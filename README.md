#CoLoRe - Cosmological Lognormal Realizations


##1 Methods.

CoLoRe is a parallel C code for generating fast mock realizations
of a given galaxy sample using a lognormal model for the matter density.
The process is as follows:
1. Generate a Gaussian realization of the linearized density field at z=0, as well as the corresponding linear radial velocity field. This is done in a Cartesian grid.
2. Calculate the redshift of each grid point and linearly evolve the density and velocity to that redshift. Include a linear galaxy bias (redshift-dependent) in the evolution of the overdensity field.
3. Log-normalize the density field and poisson-sample it using an input N(z). Each source is randomly placed inside its cell.
4. Compute the cosmological redshift and angular coordinates for each source, and introduce redshift-space distortions based on the local value of the velocity field.
5. Write sources to file.

The source code can be found in the folder "src/"

When in doubt, bear in mind that by default CoLoRe uses the following units:
 - Lenghts: Mpc/h
 - Angles: degrees


##2 Compilation and usage.

To compile CoLoRe, open the Makefile and edit it according to your
system. The default options (except for the paths to the external
libraries) should work for most systems. CoLoRe may be very memory-
demanding. To minimize the memory overhead use single precision
floating point (USE_PRECISION = yes).

OpenMP parallelization is enabled by setting the option USE_OMP
to "yes".

MPI parallelization is enabled by setting the option USE_MPI
to "yes".

CoLoRe uses 4 external packages:
 - GSL. The GNU Scientific Library (tested for versions 3.*)
 - FFTW. The Fastest Fourier Transform of the West (versions 3.*)
 - CFITSIO. FITS format library. This package is optional.
 - HDF5. HDF5 format library. This package is optional.
The paths to the corresponding headers and libraries should be correctly
set in the Makefile.

Once the Makefile has been editted, typing 'make' should generate
the executable 'CoLoRe'. To run CoLoRe just type

> mpirun -np <number-of-nodes> ./CoLoRe <param_file>

where <param_file> is the path to the parameter file described in
section 3.


##3 Parameter file.

The behaviour of CoLoRe is mainly controlled by the input param file. The
param file is basically a set of name-value pairs. Any blank lines, and
anything beyond a #-symbol will be ignored. We provide a sample param
file (param_sample.ini) that includes all the input parameters needed by
CoLoRe. The comments included in this file explain the meaning and
functionality of these parameters.


##4 Output.

The main output of CoLoRe is a catalogue of sources written either
as ASCII or FITS files. Each source is characterized by 4 quantities:
 - Z0 -> cosmological redshift (without RSDs)
 - RA -> right ascension
 - DEC -> declination
 - RSD -> RSD contribution to the redshift
          (i.e. true redshift = Z0 + RSD).
A python script (read_colore.py) is also provided to illustrate how
to read in CoLoRe's output.


##5 License:
CoLoRe is distributed under the GPL license (see COPYING in the root
directory). We kindly ask you to report the program's website
"https://github.com/damonge/CoLoRe" when using it.


##6 Contact:
Regarding bugs, suggestions, questions or petitions, feel free to contact
the author:
    David Alonso: david.alonso@astro.ox.ac.uk
