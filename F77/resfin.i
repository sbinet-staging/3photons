C
C   Common pour passer les resultats finaux
C
      INTEGER NRES
      PARAMETER (NRES=8)
      REAL*8 SPM2DIF(NRES),SPM2(2,NRES),VAR(2,NRES)
      COMMON /resfin/ SPM2DIF,SPM2,VAR
