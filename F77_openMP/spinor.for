C
C   Common pour passer les produits spinoriels engendres par mat3
C
      double complex S(5,5),T(5,5)
      double precision RAC8
      COMMON /spinor/ S,T,RAC8

      save /spinor/
!$OMP THREADPRIVATE(/spinor/)
