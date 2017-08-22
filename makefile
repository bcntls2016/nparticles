#
# Makefile for program '4hetddft-anisotropic'
# 3D dynamic (real-time evolution) 4-He Density Functional Theory for AN-isotropic PES's
#

COMP = ifort
CFLAGS = -c -O3 -xAVX -align array64byte -qopenmp\
		 -parallel -qopt-matmul -unroll
LD_FLAGS = -threads -I${MKLROOT}/include/fftw -mkl=parallel -qopt-matmul

# Name of the program
PROGNAME = nparticles

#   Fortran objects
objs = nparticles.o

.SUFFIXES: .f90 .f .o
$(PROGNAME): $(objs)
	$(COMP) -o $(PROGNAME) $(objs) $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS) -o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS) -o $(@) $<;

clean:
	rm -f *.o *.bak *.lst *.mod;
distclean:
	rm -f *.o *.bak *.lst *.mod $(PROGNAME);
