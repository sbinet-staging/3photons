C
C   Common pour passer les impulsions engendrees par Rambo
C
      INTEGER INP
      PARAMETER (INP=3)
      REAL*8 POUT(4,100),P1(4),P2(4)
      COMMON /ppp/ POUT,P1,P2
