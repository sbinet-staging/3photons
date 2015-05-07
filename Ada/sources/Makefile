#
#   Makefile for e+e- --> yyy calculations
#

COMMONFLAGS=
COMMONFLAGS=-O3 -march=native
COMMONFLAGS=-O3 -march=native -funroll-loops                   # 4.76000022888184E-06 s = f77 * 1.37445144175223291913
COMMONFLAGS=-O3 -march=native -funroll-loops -gnatp            # 4.28000020980835E-06 s = f77 * 1.23585129752248468481
VECTOPT    =-xHost -vec-report2
VECTOPT    =
ADAFLAGS   =$(COMMONFLAGS) $(VECTOPT)
LDFLAGS    =$(COMMONFLAGS)

LIB        =
LIBS       =

MONTEC     =mc

SOURCE     =mc.adb angle.adb cpuclock.adb hasa.adb ppp.adb resfin.adb result.adb scalar.adb spinor.adb unit.adb
OTHERSPEC  =cutpar.ads param.ads precision.ads xtrpro.ads
MCOBJ     :=$(SOURCE:%.adb=%.o) $(OTHERSPEC:%.ads=%.o) unixtime.o
MCLNK     :=$(SOURCE:%.adb=%.ali) $(OTHERSPEC:%.ads=%.ali)

.PHONY: clean

$(MONTEC): $(SOURCE) unixtime.o
	gnatmake $(ADAFLAGS) $@ -largs unixtime.o

check: $(SOURCE)
	gnatcheck $(SOURCE)

clean:
	/bin/rm -f $(MCOBJ)
	/bin/rm -f $(MCLNK)
