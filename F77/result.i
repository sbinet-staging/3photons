C
C   Common pour passer les resultats engendrees par mat(rice)
C
      INTEGER NRESUL
      PARAMETER (NRESUL=8)
      REAL*8 M2(2,2,2,NRESUL)
      COMMON /resultat/ M2
