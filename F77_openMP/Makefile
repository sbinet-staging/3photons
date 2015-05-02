#
#   Makefile for e+e- --> yyy calculations
#

FC         = ifort
FC         = gfortran

COMMONFLAGS=
COMMONFLAGS= -O3 -march=native
COMMONFLAGS= -O3 -march=native -fopenmp
VECTOPT    = -xHost -vec-report2
VECTOPT    =
FFLAGS     = $(COMMONFLAGS) $(VECTOPT)
FFLAGS     = -fno-range-check $(COMMONFLAGS) $(VECTOPT)
LDFLAGS    = $(COMMONFLAGS)

LINKER     = $(FC)

# pour debuguer
LIB        = -L/usr/lib/i386-linux-gnux
LIB        =
LIBS       = -lgenlib -lkernlib -lpacklib
LIBS       = -lkernlib -lpacklib
LIBS       =
MONTEC     = mc

MCOBJ      = mc.o amps.o mat.o rambo.o hasa.o mt19937_mod.o cut.o fawzi.o eric.o
SOURCE    := $(MCOBJ:%.o=%.f)

.PHONY: clean

$(MONTEC): $(MCOBJ)
	$(LINKER) $(LDFLAGS) -o $@ $^ $(LIB) $(LIBS)

check: $(SOURCE)
	ftnchek $(SOURCE)

clean:
	/bin/rm -f $(MCOBJ)

# I DID figure out the dependencies on the header files.
# (Via mkmf)
amps.o  : spinor.for
book.o  : paw.for ppp.for angle.for xtrpro.for
cut.o   : ppp.for angle.for cutpar.for
eric.o  : param.for ppp.for resfin.for
fawzi.o : cutpar.for param.for ppp.for resfin.for
hnorm.o : paw.for
mat.o   : scalar.for result.for param.for ppp.for spinor.for angle.for
mc.o    : ppp.for spinor.for scalar.for angle.for result.for cutpar.for \
          resfin.for param.for paw.for xtrpro.for
