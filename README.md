# Orientational Interpolation

This code implements (in Python 2, using NumPy) the orientational
interpolation (OI) method described in the following paper:

R. A. Shah, P. Guyot-Sionnest, and
S. K. Gray. [Orientational
Interpolation of the Optical Spectra of Nonspherical
Nanoparticles](http://dx.doi.org//10.1021/jp300621q). *J. Phys. Chem. C*, 2012, 116(23), 12712-12724.

OI was developed for connecting fixed-orientation simulations of the
optical spectra of a symmetric metal nanoparticle (e.g., with
finite-difference time-domain methods) to the orientation-averaged
extinction spectra that one measures in the lab when synthesizing
them. However, the method is much more general than that. I'm
optimistic that it could be a useful tool in any situation where one
is interested in an orientation-dependent measurement of a symmetric
object that one can sample only sparsely.

If you do use this code, I'd be delighted to hear about it in a GitHub
issue or at the gmail address in the git log. My co-authors and I
would be thankful if you cite the above paper in any related work you
publish.

## Not under active development

This is a pretty spare proof of concept code covering the D_5 and
D_infinity rotational point groups in the examples we showed in the
above paper. It also reflects my programming practices in 2010 when I
was considerably less sophisticated about software. Rather than
attempt to refactor and expand this code, at the time of writing
(2016) I'm hoping to build a Haskell implementation of OI that covers
all point groups and offers more niceties.

Of course, such goals are a moving target, and I'm willing to
reconsider if you want to use OI for something cool. Keep in mind that
I'm no longer a physical scientist by trade, so this is a hobby
project for me.

## Acknowledgment

We thank V. Kitaev for helpful synthetic advice and J. C.  Mitchell
for discussions regarding grids on SO(3). R.A.S. was supported by an
NDSEG Fellowship from the U.S. Department of Defense, Air Force Office
of Scientific Research, 32 CFR 168a, as well as a National Science
Foundation Graduate Research Fellowship. Facilities for synthesis and
characterization of nanoparticles were provided by NSF-MRSEC.
Synthesis was supported by NSF under Grant NSF-CHE
1111799. Simulations were performed under proposals CNM- 21175 and
CNM-24645 at Argonne National Laboratory supported by the
U.S. Department of Energy, Office of Science, Office of Basic Energy
Sciences, under Contract DEAC02-06CH11357.
