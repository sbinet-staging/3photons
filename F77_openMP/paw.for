C
C   Common pour PAW
C
      INTEGER NWPAWC
      PARAMETER (NWPAWC=40000)
      REAL*4 H(NWPAWC)
      COMMON /PAWC/ H
      save /pawc/
!$OMP THREADPRIVATE(/pawc/)
