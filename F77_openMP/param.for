C
C   Common pour passer les parametres physiques
C
      double precision MZ0,GZ0,GA,GBP,GBM,POLP,POLM,POLP2,POLM2
      LOGICAL IMPR
      COMMON /param/ MZ0,GZ0,GA,GBP,GBM,POLP,POLM,POLP2,POLM2,IMPR
      save /param/
!$OMP THREADPRIVATE(/param/)
