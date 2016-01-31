#!/usr/bin/python
#
# oravg.py - An engine to calculate orientation averaged cross sections
# and orientation-dependent spectral information starting from several
# fixed orientation spectra of a nanoparticle without azimuthal symmetry
# but ideally with some discrete rotational symmetries. Uses interpolation
# on SO(3) with a linear combination of translates of a single positive
# definite basis function. The basis functions used have closed form values
# for their Haar integrals over SO(3). The interpolation scheme is based on
# Filbir F, Schmid D, J. Approx. Theory 2008, 153, 170-183.
#
#
# Works with FDTD runs using x-propagating, z-polarized plane waves.
#
# To use: need a system with Python 2.5.2 and NumPy 1.1.0.
# Command line on Unix-like system:
#
# chmod +x ./oravg.py
# ./oravg.py <control_file >output_file
#
#
# CONTROL FILE: text file with a header, choice of interpolator and point
# group, and body.
#
# HEADER: define the job type and parameters. Choices currently are
# orientation-averaged spectrum ("o"), fixed orientation spectrum ("f"),
# or polar view of cross section as particle is rotated ("p"). Put one of
# the following at the top of the control file:
#
# Orientation-averaged spectrum:
# o                          <- Job type
#
# Fixed orientation spectrum:
# f                          <- Job type
# 0.0 0.0 0.0 1.0            <- Unit quaternion (x,y,z,w) defining desired
#                               orientation
#
# Polar view:
# p                          <- Job type
# 0.0 0.0 0.0 1.0            <- Unit quaternion defining desired orientation
# 1.0 0.0 0.0                <- Unit vector defining axis about which particle
#                               is to be rotated
# 101                        <- Number of angular samples from 0 to 2pi
#                               (samples include 0 and 2pi)
# 672E-9, 1E-9               <- Wavelength/frequency to sample and tolerance
#                               between given value and value in data file
#                               (in same units as in data file)
#
# CHOICE OF INTERPOLATOR AND POINT GROUP: Put the following next in the
# control file:
#
# symm_interp                <- Family of interpolator to use
# 4                          <- Degree of interpolator in family
# d5                         <- Rotational point group of nanoparticle
#
# BODY: Put the following at the end of the control file:
#
# 0.0 0.0 0.0 1.0 76d.dat F F
# 0.7071 0.0 0.0 0.7071 75.dat F F
# 0.0 0.7071 0.0 0.7071 76f.dat T F
# -0.2706 0.2706 0.0 0.9239 76h.dat T T
#
# where first 4 fields represent the unit quaternion (x,y,z,w) by which
# the particle is rotated; the fifth is the file name (by default in same
# subdirectory as this program) of the corresponding spectrum, and the
# last two are whether reflection in the xy and xz planes, respectively,
# give a different but valid orientation with an identical spectrum.
#
#
# DATA FILES:
# The data files must be lists of wavelength or frequency followed by the
# cross section, one entry per line, with no blank lines:
#
# 7.49481103E-07  551.72174
#
# The list of wavelengths sampled must be identical between files.
#
#
# Raman Shah, April-May 2010
#
# Modified for minor bug fixes 8 Jun 2010, 20 Aug 2010


# IMPORTED MODULES

from math import *
import numpy as np


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

def chain_multiply(quat_list):
    # Recursively multiplies an arbitrary list of quaternions,
    # returning a quaternion.
    # q1 * q2 * q3 * q4 evaluated as ((q1 * q2) * q3) * q4, etc.
    if len(quat_list) == 0:
        raise IndexError
    if len(quat_list) == 1:
        return quat_list[0]
    if len(quat_list) == 2:
        return multiply(quat_list[0], quat_list[1])
    else:
        right_term = quat_list.pop()
        return multiply(chain_multiply(quat_list), right_term)


# FUNCTIONS FOR SYMMETRY ADAPTATION OF DATA SET

def d5symadapt(input_tuple):
    # Takes a 2-tuple (quaternion + filename) and returns a
    # list of the 10 tuples equivalent to the given tuple
    # under D_5 symmetry (first output is the parent).
    # Right-handed rotation by 2pi/5 about the z axis
    c5 = quaternion(0, 0, sin(pi/5), cos(pi/5))
    # Right-handed rotation by pi about the x axis
    perpc2 = quaternion(1, 0, 0, 0)
    quat = input_tuple[0]
    text = input_tuple[1]
    adapted_list = [ \
        (quat, text), \
        (chain_multiply([quat, c5]), text), \
        (chain_multiply([quat, c5, c5]), text), \
        (chain_multiply([quat, c5, c5, c5]), text), \
        (chain_multiply([quat, c5, c5, c5, c5]), text), \
        (chain_multiply([quat, perpc2]), text), \
        (chain_multiply([quat, c5, perpc2]), text), \
        (chain_multiply([quat, c5, c5, perpc2]), text), \
        (chain_multiply([quat, c5, c5, c5, perpc2]), text), \
        (chain_multiply([quat, c5, c5, c5, c5, perpc2]), text)]
    return adapted_list

def choose_point_group(point_group):
    if point_group == 'd5':
        return d5symadapt


# DEFINITION OF INTERPOLATORS AND THEIR INTEGRALS

def choose_interpolator(int_choice, n):
    # The Filbir family of interpolating basis functions on SO(3)
    # Filbir F, Schmid D, J. Approx. Theory 153 (2008) 170.
    if int_choice == 'filbir':
        interpolator = (lambda quat: \
            ((quat.w * quat.w - quat.z * quat.z + 1) / 2) ** n)
        # integral[n] gives the Haar integral of (...)^n
        # for n up to 15. Note that the first entry is integral[0]!
        integral = [1.0, \
                    1.0/2, \
                    7.0/24, \
                    3.0/16, \
                    83.0/640, \
                    73.0/768, \
                    523.0/7168, \
                    119.0/2048, \
                    14051.0/294912, \
                    13103.0/327680, \
                    98601.0/2883584, \
                    15565.0/524288, \
                    1423159.0/54525952, \
                    1361617.0/58720256, \
                    10461043.0/503316480, \
                    1259743.0/67108864]
        return (interpolator, integral[n])

    # Maybe a better choice? Symmetrically uses all three
    # functions phi_1,0,0, phi_0,1,0, phi_0,0,1 in the paper,
    # adds 2, and analogously normalizes to 1
    if int_choice == 'symm_interp':
        interpolator = (lambda quat: \
            ((4 * quat.w * quat.w + 1) / 5) ** n)
        integral = [1.0, \
                    2.0/5, \
                    1.0/5, \
                    3.0/25, \
                    51.0/625, \
                    188.0/3125, \
                    731.0/15625, \
                    118.0/3125, \
                    2447.0/78125, \
                    51822.0/1953125, \
                    223191.0/9765625, \
                    974427.0/48828125, \
                    860529.0/48828125, \
                    767244.0/48828125, \
                    17242377.0/1220703125, \
                    78049611.0/6103515625]
        return (interpolator, integral[n])


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
int_choice = readlist()[0]
n = int(readlist()[0])
point_group = readlist()[0]
interpolator, integral = choose_interpolator(int_choice, n)
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

# Instantiate array of quaternion arguments to be passed to the
# interpolator. This can also be used to look at the separation
# quality of the data set.
# a_ij = q_j^-1 * q_i
dimension = len(main_quat_list)
arg_array = []
for i in range(dimension):
    row = []
    for j in range(dimension):
        row.append(multiply(main_quat_list[j].invert(), \
                                main_quat_list[i]))
    arg_array.append(row)

# Instantiate the main matrix of function evaluations:
# A_ij = Phi_n(a_ij)
A = np.matrix(np.zeros((dimension, dimension)))
for i in range(dimension):
    for j in range(dimension):
        A[i,j] = interpolator(arg_array[i][j])

## Some verbose info
#
## Find separation distance of the set of quaternions
#biggest_w = 0.0
#for i in range(dimension):
#    for j in range(dimension):
#        if (i != j) and arg_array[i][j].w > biggest_w:
#            biggest_w = arg_array[i][j].w
#print 'Biggest w is ' + str(biggest_w)
#print 'Smallest separation angle (in rad) is ' + str(2*acos(biggest_w))
#
#print 'Dimension is ' + str(dimension) + '.'
#
#for i in range(dimension):
#    for j in range(dimension):
#        if ((i != j) and A[i,j] == 1):
#            print 'Off-diagonal element ' + str(i+1) + ', ' + str(j+1) + \
#                ' is 1.'
#        if A[i,j] == 0:
#            print 'Element ' + str(i+1) + ', ' + str(j+1) + ' is 0.'

# Invert the matrix A to get the matrix that converts data
# to coefficients
Ainv = np.linalg.inv(A)


# CALCULATION AND OUTPUT OF ORIENTATION-AVERAGED SPECTRUM

# Spectrum-like jobs: main loop increments line number for data files
if job_type in ['o', 'f']:
    if job_type == 'f':
        eval_row = np.matrix(np.zeros((1, dimension)))
        for j in range(dimension):
            eval_row[0,j] = interpolator(multiply(main_quat_list[j].invert(), \
                                                      orientation))
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
            eval_row[0,j] = interpolator(multiply(main_quat_list[j].invert(), \
                                                      curr_quat))

        print angle, (eval_row * coef_vector)[0,0]


# FINAL I/O CLEANUP

for filename in datafile_dict:
    datafile_dict[filename].close()
