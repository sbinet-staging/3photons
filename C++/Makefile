#
#   Makefile for e+e- --> yyy calculations
#

CXX        =icc
CXX        =g++

LINKER     =$(CXX)

COMMONFLAGS=
COMMONFLAGS=-O3 -march=native -fopenmp
COMMONFLAGS=-O3 -march=native

VECTOPT    =-xHost -vec-report2
VECTOPT    =
CXXFLAGS   =$(COMMONFLAGS) $(VECTOPT)
CFLAGS     =$(COMMONFLAGS) $(VECTOPT)
LDFLAGS    =$(COMMONFLAGS)

LIB        =
LIBS       =-lstdc++
LIBS       =

MONTEC     =mc

MCOBJ      =mc.o angle.o ppp.o spinor.o result.o resfin.o hasa.o unixtime.o
SOURCE    :=$(MCOBJ:%.o=%.cpp)

.PHONY: clean

$(MONTEC): $(MCOBJ)
	$(LINKER) $(LDFLAGS) -o $@ $^ $(LIB) $(LIBS)

check: $(SOURCE)
	cppcheck $(SOURCE)

clean:
	/bin/rm -f $(MCOBJ)
