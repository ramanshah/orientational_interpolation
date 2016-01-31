#!/usr/bin/python
#
# oravg_s2.py - An engine, adapted from the original oravg.py, to
# calculate orientation averaged cross sections and
# orientation-dependent spectral information starting from several
# fixed orientation spectra of a nanoparticle with axial
# symmetry. Uses interpolation on the sphere S^2 with a linear
# combination of translates of a single strictly positive definite
# basis function. The interpolation scheme is based on the "generating
# function kernel" in:
#
# Jetter K, Stoeckler J, Ward JD, Math. Computation 1999, 68, 733-747.
#
# Works with FDTD runs using x-propagating, z-polarized plane
# waves. It is assumed that the particle's axis of rotational symmetry
# is the z-axis when the quaternion is (0,0,0,1).
#
# To use: need a system with Python 2.5.2 and NumPy 1.1.0.
# Command line on Unix-like system:
#
# chmod +x ./oravg_s2.py
# ./oravg_s2.py <control_file >output_file
#
#
# CONTROL FILE: text file with a header, choice of interpolation
# parameters, and body.
#
# HEADER: define the job type and related parameters. Choices
# currently are orientation-averaged spectrum ("o"), fixed orientation
# spectrum ("f"), or polar view of cross section as particle is
# rotated ("p"). Put one of the following at the top of the control
# file:
#
# Orientation-averaged spectrum:
# o                          <- Job type
#
# Fixed orientation spectrum:
# f                          <- Job type
# 0.0 0.38268 0.0 0.92388    <- Unit quaternion (x,y,z,w) defining desired
#                               orientation
# Polar view:
# p                          <- Job type
# 0.0 0.0 0.0 1.0            <- Unit quaternion defining desired orientation
# 1.0 0.0 0.0                <- Unit vector defining axis about which particle
#                               is to be rotated
# 101                        <- Number of angular samples from 0 to 2pi
#                               (samples include 0 and 2pi)
# 667, 1.0                   <- Wavelength/frequency to sample and tolerance
#                               between given value and value in data file
#                               (in same units as in data file)
#
# CHOICE OF INTERPOLATION PARAMETERS: Put the following next in the
# control file:
#
# 0.3                        <- Peakedness z of the kernel (0 < z < 1)
# dinfh                      <- Point group of nanoparticle
#
# BODY: Put the following at the end of the control file:
#
# 0.0 0.0 0.0 1.0 30.dat F F
# 0.707106781187 0.0 0.0 0.707106781187 30a.dat F F
# 0.0 0.707106781187 0.0 0.707106781187 30b.dat F F
# 0.325057583672 0.325057583672 0.0 0.888073833977 30c.dat T T
#
# where first 4 fields represent the unit quaternion (x,y,z,w) by
# which the particle is rotated; the fifth is the data file name (by
# default in same subdirectory as this program) of the corresponding
# spectrum, and the last two are whether reflection in the xy and xz
# planes, respectively, give a different but valid orientation with an
# identical spectrum.
#
#
# DATA FILES:
# The data files must be lists of wavelength or frequency followed by the
# cross section, one entry per line, with no blank lines:
#
# 810.2498 1.6144834E-12
#
# The list of wavelengths sampled must be identical between files.
#
#
# Raman Shah, April-August 2010


# IMPORTED MODULES

from math import *
import numpy as np
from itertools import imap
from operator import mul


# MISCELLANEOUS UTILITIES

# Convert Ts and Fs in file to boolean values in Python
# Note: Text must be trimmed of any leading whitespace!
def parsebool(string):
    return string[0].upper() == 'T'

# Read standard input until a nonblank like is found, and return the nonblank
# line as a list of words. Recognize comma as a word separator.
def readlist():
    while True:
        line = raw_input().replace(',', ' ').split()
        if line != []:
            return line


# QUATERNION CLASS AND RELATED METHODS

class quaternion:
    """A self-normalized rotation quaternion (x,y,z,w)."""
    def __init__(self, xpart, ypart, zpart, wpart):
        norm = sqrt(xpart * xpart + ypart * ypart + zpart * zpart \
            + wpart * wpart)
        if norm == 0:
            raise ValueError
        self.x = xpart / norm
        self.y = ypart / norm
        self.z = zpart / norm
        self.w = wpart / norm

    def invert(self):
        # Quaternion for the inverse rotation has same axis
        # but opposite angle, so the vector components flip sign
        return quaternion(-1 * self.x, -1 * self.y, -1 * self.z, self.w)

    def xyreflect(self):
        return quaternion(-1 * self.x, -1 * self.y, self.z, self.w)

    def xzreflect(self):
        return quaternion(-1 * self.x, self.y, -1 * self.z, self.w)

    def show(self):
        print self.x, self.y, self.z, self.w

def multiply(quat1, quat2):
    prodw = quat1.w*quat2.w - quat1.x*quat2.x \
        - quat1.y*quat2.y - quat1.z*quat2.z
    prodx = quat1.w*quat2.x + quat1.x*quat2.w \
        + quat1.y*quat2.z - quat1.z*quat2.y
    prody = quat1.w*quat2.y + quat1.y*quat2.w \
        + quat1.z*quat2.x - quat1.x*quat2.z
    prodz = quat1.w*quat2.z + quat1.z*quat2.w \
        + quat1.x*quat2.y - quat1.y*quat2.x
    return quaternion(prodx,prody,prodz,prodw)


# FUNCTIONS FOR WORKING IN S^2

def point(quat):
    # Gives the unit vector (x, y, z) in R^3 resulting when (0, 0, 1)
    # is rotated by the given quaternion. The third column of the
    # quaternion rotation matrix.
    return (2 * (quat.w * quat.y + quat.x * quat.z),
            2 * (quat.y * quat.z - quat.x * quat.w),
            (quat.w * quat.w - quat.x * quat.x - \
                     quat.y * quat.y + quat.z * quat.z))

def dot(x, y):
    # Dot product of two vectors in R^3, written as tuples
    return sum(imap(mul, x, y))


# FUNCTIONS FOR SYMMETRY ADAPTATION OF DATA SET

def cinfv_symadapt(input_tuple):
    return [input_tuple]

def dinfh_symadapt(input_tuple):
    # Takes a 2-tuple (quaternion + filename) and returns a
    # list of the 2 tuples equivalent to the given tuple via the discrete
    # subgroup of D_infh.
    # Right-handed rotation by pi about the x axis
    perpc2 = quaternion(1, 0, 0, 0)
    quat = input_tuple[0]
    text = input_tuple[1]
    return [(quat, text), (multiply(quat, perpc2), text)]

def choose_point_group(point_group):
    if point_group == 'cinfv':
        return cinfv_symadapt
    if point_group == 'dinfh':
        return dinfh_symadapt


# DEFINITION OF INTERPOLATOR AND ITS INTEGRAL

def choose_interpolator(z):
    # This is the generating function kernel, normalized to 1 for
    # cos_theta = 1, and its integral over S^2.
    interpolator = (lambda cos_theta: \
        ((1 - z) ** 3.0 / (1 + z ** 2 - 2 * z * cos_theta) ** 1.5))
    integral = (1 - z) ** 2 / (1 + z)

    return (interpolator, integral)


# EXECUTION BEGINS HERE: READ CONTROL FILE TO INITIALIZE

# Initialize job using header of control file piped to standard input
job_type = readlist()[0]
if job_type == 'o':
    pass
elif job_type == 'f':
    components = map(float, readlist())
    orientation = quaternion(components[0], components[1], \
                             components[2], components[3])
elif job_type == 'p':
    components = map(float, readlist())
    orientation = quaternion(components[0], components[1], \
                             components[2], components[3])
    axis = map(float, readlist())
    stepcount = int(readlist()[0])
    wavelength, wavelength_tol = map(float, readlist())
else:
    raise IOError('Job type read as ' + job_type + ' is not valid.')


# Define interpolator, integral, rot_adapt from control file
z = float(readlist()[0])
if not ((z > 0) and (z < 1)):
    raise ValueError('Need 0 < z < 1 to give a valid kernel.')

point_group = readlist()[0]
interpolator, integral = choose_interpolator(z)
rot_adapt = choose_point_group(point_group)

# Parse input from body of control file into a control array with
# entries [quaternion, string, boolean, boolean]. Also create a
# dictionary of open data files to read later in program.
control_list = []
datafile_dict = {}
while True:
    try:
        words = readlist()
        curr_quat = quaternion(float(words[0]), float(words[1]), \
                               float(words[2]), float(words[3]))
        xybool = parsebool(words[5])
        xzbool = parsebool(words[6])
        control_list.append([curr_quat, words[4], xybool, xzbool])
        datafile_dict[words[4]] = open('./' + words[4], 'r')
    except EOFError:
        break


# ADAPT DATA SET TO SYMMETRY OF PROBLEM

# Add reflected quaternions as directed in the control list
# giving an array of tuples (quaternion, string)
reflected_list = []
for item in control_list:
    reflected_list.append((item[0], item[1]))
    if (item[2] and item[3]):
        #xy and xz both give new equivalent configurations
        reflected_list.append((item[0].xyreflect(), item[1]))
        reflected_list.append((item[0].xzreflect(), item[1]))
        reflected_list.append((item[0].xyreflect().xzreflect(), item[1]))
        continue
    if item[2]:
        #only xy gives a new equivalent configuration
        reflected_list.append((item[0].xyreflect(), item[1]))
        continue
    if item[3]:
        #only xz gives a new equivalent configuration
        reflected_list.append((item[0].xzreflect(), item[1]))
        continue

# For each quaternion, add the configurations equivalent
# under the rotational point group, giving a master array
# of tuples (quaternion, string) ready for processing
master_list = []
for item in reflected_list:
    master_list += rot_adapt(item)

# Give tuples of the quaternions and associated filenames
main_quat_list, main_filename_list = zip(*master_list)


# NUMERICAL SECTION: FUNCTION EVALUATION AND MATRIX INVERSION

# Instantiate a matrix of the dot products between corresponding points in
# R^3 as well as a matrix of the interpolator evaluations. The former is useful
# for measuring the separation distance of the data set.
# a_ij = x_i dot x_j
# A_ij = Phi_z(a_ij)
dimension = len(main_quat_list)
a = np.matrix(np.zeros((dimension, dimension)))
A = np.matrix(np.zeros((dimension, dimension)))
for i in range(dimension):
    for j in range(dimension):
        a[i, j] = dot(point(main_quat_list[i]), point(main_quat_list[j]))
        A[i, j] = interpolator(a[i, j])

# # Some verbose info
#
# # Find separation distance of the set of quaternions
# biggest_costheta = 0
# for i in range(dimension):
#     for j in range(dimension):
#         if (i != j) and a[i,j] > biggest_costheta:
#             biggest_costheta = a[i,j]
# print 'Smallest separation angle (in rad) is ' + \
#    str(acos(biggest_costheta)) + '.'
#
# print 'Dimension is ' + str(dimension) + '.'
#
#print 'Points sampled are:'
#for quat in main_quat_list:
#    print point(quat)
#
# for i in range(dimension):
#     for j in range(dimension):
#         if ((i != j) and A[i,j] == 1):
#             print 'Off-diagonal element ' + str(i+1) + ', ' + str(j+1) + \
#                 ' is 1.'
#         if A[i,j] == 0:
#             print 'Element ' + str(i+1) + ', ' + str(j+1) + ' is 0.'

# Invert the matrix A to get the matrix that converts data
# to coefficients
Ainv = np.linalg.inv(A)


# CALCULATION AND OUTPUT OF RESULTS

# Spectrum-like jobs: main loop increments line number for data files
if job_type in ['o', 'f']:
    if job_type == 'f':
        eval_row = np.matrix(np.zeros((1, dimension)))
        for j in range(dimension):
            eval_row[0,j] = interpolator(dot(point(orientation), \
                                             point(main_quat_list[j])))

    while True: # Loop runs until a blank line is seen in data file
        # Initialize the iteration
        cs_vector = np.matrix(np.zeros((dimension, 1)))
        coef_vector = np.matrix(np.zeros((dimension, 1)))
        curr_line_dict = {}

        # Read a new line, parsed into floats, to give a dictionary
        for filename in datafile_dict:
            curr_line_dict[filename] = \
                map(float, datafile_dict[filename].readline().split())

        if curr_line_dict.values()[0] == []:
            break

        # Look up the wavelength and cross sections for this iteration
        curr_wavelength = curr_line_dict.values()[0][0]
        for i in range(dimension):
            cs_vector[i] = curr_line_dict[main_filename_list[i]][1]

        # Compute the vector of coefficients and the sum of
        # coefficients
        coef_vector = Ainv * cs_vector

        if job_type == 'o':
            print curr_wavelength, np.sum(coef_vector) * integral

        if job_type == 'f':
            print curr_wavelength, (eval_row * coef_vector)[0,0]

# Monochromatic jobs: main loop increments parameter; one line per file used
if job_type == 'p':
    # Look up data for the desired wavelength
    while True:
        curr_line_dict = {}

        # Read a new line, parsed into floats, to give a dictionary
        for filename in datafile_dict:
            curr_line_dict[filename] = \
                map(float, datafile_dict[filename].readline().split())

        if curr_line_dict.values()[0] == []:
            raise IOError('Wavelength/frequency ' + \
                              str(wavelength) + ' not found!')

        curr_wavelength = float(curr_line_dict.values()[0][0])

        if abs(curr_wavelength - wavelength) < wavelength_tol:
            # Comment out following line to get pure stream of numbers
            print '# Wavelength/frequency sampled is: ' + str(curr_wavelength)
            break

    # Find vector of coefficients just once for the chosen wavelength
    cs_vector = np.matrix(np.zeros((dimension, 1)))
    for i in range(dimension):
        cs_vector[i] = curr_line_dict[main_filename_list[i]][1]

    coef_vector = Ainv * cs_vector

    angle_list = np.linspace(0, 2*pi, stepcount, True)

    # Main loop
    for angle in angle_list:
        curr_quat = multiply(quaternion(axis[0] * sin(angle / 2), \
                                        axis[1] * sin(angle / 2), \
                                        axis[2] * sin(angle / 2), \
                                        cos (angle / 2)), \
                             orientation)

        eval_row = np.matrix(np.zeros((1, dimension)))
        for j in range(dimension):
            eval_row[0,j] = interpolator(dot(point(curr_quat), \
                                             point(main_quat_list[j])))
        print angle, (eval_row * coef_vector)[0,0]


# FINAL I/O CLEANUP

for filename in datafile_dict:
    datafile_dict[filename].close()
