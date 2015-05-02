C
C   Common pour passer les resultats finaux
C
      INTEGER NRES
      PARAMETER (NRES=8)
      double precision SPM2DIF(NRES),SPM2(2,NRES),VAR(2,NRES)
      COMMON /resfin/ SPM2,VAR
!      COMMON /resfin/ SPM2DIF,SPM2,VAR
      save /resfin/
!OMP SHARED(/resfin/)
