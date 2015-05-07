#
#   Makefile for e+e- --> yyy calculations
#

F90        =g95
F90        =ifort
F90        =gfortran

LINKER     =$(F90)

COMMONFLAGS=
COMMONFLAGS=-O3 -march=native -fno-range-check -funroll-loops -finline# # 6.1519999504089354E-006 s = f77 * 1.77639176363769144091
COMMONFLAGS=-O3 -march=native -fno-range-check -funroll-loops -fopenmp
COMMONFLAGS=-O3 -march=native -fno-range-check -funroll-loops#          # 5.5640001296997072E-006 s = f77 * 1.60660664547316982055

VECTOPT    =-xHost -vec-report2
VECTOPT    =
FFLAGS     =$(COMMONFLAGS) $(VECTOPT)
LDFLAGS    =$(COMMONFLAGS)

LIB        =
LIBS       =-lgenlib -lkernlib -lpacklib
LIBS       =-lkernlib -lpacklib
LIBS       =

MONTEC     =mc

MCOBJ      = precision.o angle.o cutpar.o param.o ppp.o resfin.o result.o scalar.o spinor.o xtrpro.o \
             mc.o amps.o mat.o rambo.o hasa.o cut.o fawzi.o eric.o # book.o
SOURCE    :=$(MCOBJ:%.o=%.f90)
MCMOD     :=$(MCOBJ:%.o=%.mod)

.SUFFIXES: .f90 .o

.f90.o:
	$(F90) $(FFLAGS) $(CFLAGS) -c $< -o $@
#.mod.o:
#	$(F90) $(FFLAGS) $(CFLAGS) $*.f90

.PHONY: clean

$(MONTEC): $(MCOBJ)
	$(LINKER) $(LDFLAGS) -o $@ $^ $(LIB) $(LIBS)

check: $(SOURCE)
	ftnchek $(SOURCE)

clean:
	/bin/rm -f $(MCOBJ) $(MCMOD)
	/bin/rm -f *.o *.mod

# I DID figure out the dependencies on the header files.
# (Via mkmf)
amps.o	: spinor.f90 precision.f90
book.o	: paw.f90 ppp.f90 angle.f90 xtrpro.f90 precision.f90
cut.o	: ppp.f90 angle.f90 cutpar.f90 precision.f90
eric.o	: param.f90 ppp.f90 resfin.f90 precision.f90
fawzi.o	: cutpar.f90 param.f90 ppp.f90 resfin.f90 precision.f90
hnorm.o	: paw.f90
mat.o	: scalar.f90 result.f90 param.f90 ppp.f90 spinor.f90 angle.f90 precision.f90
mc.o	: ppp.f90 spinor.f90 scalar.f90 angle.f90 result.f90 cutpar.f90 resfin.f90 \
	param.f90 paw.f90 xtrpro.f90 precision.f90
