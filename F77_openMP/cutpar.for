C
C   Common pour passer les coupures
C
      double precision ACUT,BCUT,EMIN,SINCUT
      COMMON /cutpar/ ACUT,BCUT,EMIN,SINCUT
      save /cutpar/
!$OMP THREADPRIVATE(/cutpar/)
