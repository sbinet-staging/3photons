C
C   Common pour passer les produits scalaires engendres par mat3
C
      double precision PS(5,5)
      COMMON /scalar/ PS

      save /scalar/
!$OMP THREADPRIVATE(/scalar/)
