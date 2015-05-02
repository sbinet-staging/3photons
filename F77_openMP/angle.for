C
C   Common pour passer les angles
C
      double precision COS12,COS1P1,COS2P1,COS3P1,COS13,COS23,COSN,COSAC
      COMMON /angle/ COS12,COS1P1,COS2P1,COS3P1,COS13,COS23,COSN,COSAC
      save /angle/
!$OMP THREADPRIVATE(/angle/)
