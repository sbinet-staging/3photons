#
#   Makefile for e+e- --> yyy calculations
#

#FC         = fort77
#FC         = f77
FC         = ifort
FC         = gfortran
FC         = g77
FC         = fort77

#FFLAGS	= -G -K +ppu +E1 -C -R8 -c +U77 -w66 +O3
#LDFLAGS	= -G -K +ppu +E1 -C +U77 +O3
# pour debuguer
#FFLAGS	   = -g -K +ppu +E1 -C -R8 -c +U77 -w66 +O3
#LDFLAGS   = -g -K +ppu +E1 -C +U77 +O3
# pour optimiser la source
# +OP3P    # Optimizer
FFLAGS	   = -c
FFLAGS	   = -c -O3 -m32
FFLAGS	   = -c -O3
LDFLAGS   = 
LDFLAGS   = -melf_i386 -O3 -m32
LDFLAGS   = -O3

LINKER	   = $(FC)

# pour debuguer
#LIB        = -L/cern/pro/lib -L/usr/lib/end.o
LIB        = -L/cern/pro/lib
LIB        = -L/usr/lib/i386-linux-gnu
LIB        =
LIBS       = -lgenlib -lkernlib -lpacklib
LIBS       = -lkernlib -lpacklib

MONTEC	= mc

MCOBJ	= mc.o amps.o mat.o rambo.o hasa.o cut.o book.o fawzi.o eric.o

$(MONTEC):	$(MCOBJ)
		$(LINKER) $(LDFLAGS) -o $(MONTEC) $(MCOBJ) $(LIB) $(LIBS)

# I DID figure out the dependencies on the header files.
# (Via mkmf)
amps.o	: spinor.i
book.o	: paw.i ppp.i angle.i xtrpro.i
cut.o	: ppp.i angle.i cutpar.i
eric.o	: param.i ppp.i resfin.i
fawzi.o	: cutpar.i param.i ppp.i resfin.i
hnorm.o	: paw.i
mat.o	: scalar.i result.i param.i ppp.i spinor.i angle.i
mc.o	: ppp.i spinor.i scalar.i angle.i result.i cutpar.i resfin.i \
	param.i paw.i xtrpro.i
