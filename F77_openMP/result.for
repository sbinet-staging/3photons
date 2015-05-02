C
C   Common pour passer les resultats engendrees par mat(rice)
C
      INTEGER NRESUL
      PARAMETER (NRESUL=8)
      double precision M2(2,2,2,NRESUL)
      COMMON /result/ M2
      save /result/
!$OMP THREADPRIVATE(/result/)
