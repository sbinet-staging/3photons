!----------------------------------------------------------------------
SUBROUTINE BOOK (PROBA, PROBA2)
!
! HBOOK les resultats
!
!----------------------------------------------------------------------
  use precision
  use ppp
  use angle
  use xtrpro
  implicit none
  INCLUDE 'paw.f90'
  REAL (pr) :: PROBA, PROBA2
  REAL*4 PROBAPAW, PROBAPAW2

  PROBAPAW=SNGL (PROBA)
  PROBAPAW2=SNGL (PROBA2)

!    STOCKE DANS UN HISTOGRAMME LES PARAMETRES DE L'EVENEMENT

  CALL HFILL (600, SNGL (COS1P1), SNGL (POUT (4, 3)), PROBAPAW2)
  CALL HFILL (601, SNGL (COS1P1), SNGL (POUT (4, 3)), PRPLUS   )
  CALL HFILL (602, SNGL (COS1P1), SNGL (POUT (4, 3)), PRMOINS  )

  CALL HFILL (610, SNGL (POUT (4, 2)), SNGL (POUT (4, 3)), PROBAPAW2)
  CALL HFILL (611, SNGL (POUT (4, 2)), SNGL (POUT (4, 3)), PRPLUS   )
  CALL HFILL (612, SNGL (POUT (4, 2)), SNGL (POUT (4, 3)), PRMOINS  )

  CALL HFILL (620, SNGL (ABS (COSN)), SNGL (POUT (4, 3)), PROBAPAW2)
  CALL HFILL (621, SNGL (ABS (COSN)), SNGL (POUT (4, 3)), PRPLUS   )
  CALL HFILL (622, SNGL (ABS (COSN)), SNGL (POUT (4, 3)), PRMOINS  )

  CALL HFILL (630, SNGL (ABS (COSN)), SNGL (POUT (4, 1)), PROBAPAW2)
  CALL HFILL (631, SNGL (ABS (COSN)), SNGL (POUT (4, 1)), PRPLUS   )
  CALL HFILL (632, SNGL (ABS (COSN)), SNGL (POUT (4, 1)), PRMOINS  )

  CALL HFILL (200, SNGL (COSAC), 0., PROBAPAW2)
  CALL HFILL (210, SNGL (COS12), 0., PROBAPAW2)
  CALL HFILL (220, SNGL (COS13), 0., PROBAPAW2)
  CALL HFILL (230, SNGL (COS23), 0., PROBAPAW2)
  CALL HFILL (240, SNGL (COS1P1), 0., PROBAPAW2)
  CALL HFILL (250, SNGL (COS2P1), 0., PROBAPAW2)
  CALL HFILL (260, SNGL (COS3P1), 0., PROBAPAW2)
  CALL HFILL (270, SNGL (POUT (4, 1)), 0., PROBAPAW2)
  CALL HFILL (280, SNGL (POUT (4, 2)), 0., PROBAPAW2)
  CALL HFILL (290, SNGL (POUT (4, 3)), 0., PROBAPAW2)
  CALL HFILL (295, SNGL (ABS (COSN)), 0., PROBAPAW2)

  CALL HFILL (201, SNGL (COSAC), 0., PRPLUS )
  CALL HFILL (211, SNGL (COS12), 0., PRPLUS )
  CALL HFILL (221, SNGL (COS13), 0., PRPLUS )
  CALL HFILL (231, SNGL (COS23), 0., PRPLUS )
  CALL HFILL (241, SNGL (COS1P1), 0., PRPLUS )
  CALL HFILL (251, SNGL (COS2P1), 0., PRPLUS )
  CALL HFILL (261, SNGL (COS3P1), 0., PRPLUS )
  CALL HFILL (271, SNGL (POUT (4, 1)), 0., PRPLUS )
  CALL HFILL (281, SNGL (POUT (4, 2)), 0., PRPLUS )
  CALL HFILL (291, SNGL (POUT (4, 3)), 0., PRPLUS )
  CALL HFILL (296, SNGL (ABS (COSN)), 0., PRPLUS )

  CALL HFILL (202, SNGL (COSAC), 0., PRMOINS)
  CALL HFILL (212, SNGL (COS12), 0., PRMOINS)
  CALL HFILL (222, SNGL (COS13), 0., PRMOINS)
  CALL HFILL (232, SNGL (COS23), 0., PRMOINS)
  CALL HFILL (242, SNGL (COS1P1), 0., PRMOINS)
  CALL HFILL (252, SNGL (COS2P1), 0., PRMOINS)
  CALL HFILL (262, SNGL (COS3P1), 0., PRMOINS)
  CALL HFILL (272, SNGL (POUT (4, 1)), 0., PRMOINS)
  CALL HFILL (282, SNGL (POUT (4, 2)), 0., PRMOINS)
  CALL HFILL (292, SNGL (POUT (4, 3)), 0., PRMOINS)
  CALL HFILL (297, SNGL (ABS (COSN)), 0., PRMOINS)

  RETURN
END SUBROUTINE BOOK
!----------------------------------------------------------------------

!----------------------------------------------------------------------
SUBROUTINE NORMA (COEFF)
!----------------------------------------------------------------------
!
! Effectue la normalisations des histogrammes 
!
!----------------------------------------------------------------------
  use precision
  implicit none
  REAL (pr) :: COEFF
  REAL*4 WGTPAW, ZERO
  PARAMETER (ZERO=0.E0)
  INCLUDE 'paw.f90'

  WGTPAW=SNGL (COEFF)

  CALL HOPERA (200, '+', 200, 200, WGTPAW, ZERO)
  CALL HOPERA (210, '+', 210, 210, WGTPAW, ZERO)
  CALL HOPERA (220, '+', 220, 220, WGTPAW, ZERO)
  CALL HOPERA (230, '+', 230, 230, WGTPAW, ZERO)
  CALL HOPERA (240, '+', 240, 240, WGTPAW, ZERO)
  CALL HOPERA (250, '+', 250, 250, WGTPAW, ZERO)
  CALL HOPERA (260, '+', 260, 260, WGTPAW, ZERO)
  CALL HOPERA (270, '+', 270, 270, WGTPAW, ZERO)
  CALL HOPERA (280, '+', 280, 280, WGTPAW, ZERO)
  CALL HOPERA (290, '+', 290, 290, WGTPAW, ZERO)
  CALL HOPERA (295, '+', 295, 295, WGTPAW, ZERO)

  CALL HOPERA (201, '+', 201, 201, WGTPAW, ZERO)
  CALL HOPERA (211, '+', 211, 211, WGTPAW, ZERO)
  CALL HOPERA (221, '+', 221, 221, WGTPAW, ZERO)
  CALL HOPERA (231, '+', 231, 231, WGTPAW, ZERO)
  CALL HOPERA (241, '+', 241, 241, WGTPAW, ZERO)
  CALL HOPERA (251, '+', 251, 251, WGTPAW, ZERO)
  CALL HOPERA (261, '+', 261, 261, WGTPAW, ZERO)
  CALL HOPERA (271, '+', 271, 271, WGTPAW, ZERO)
  CALL HOPERA (281, '+', 281, 281, WGTPAW, ZERO)
  CALL HOPERA (291, '+', 291, 291, WGTPAW, ZERO)
  CALL HOPERA (296, '+', 296, 296, WGTPAW, ZERO)

  CALL HOPERA (202, '+', 202, 202, WGTPAW, ZERO)
  CALL HOPERA (212, '+', 212, 212, WGTPAW, ZERO)
  CALL HOPERA (222, '+', 222, 222, WGTPAW, ZERO)
  CALL HOPERA (232, '+', 232, 232, WGTPAW, ZERO)
  CALL HOPERA (242, '+', 242, 242, WGTPAW, ZERO)
  CALL HOPERA (252, '+', 252, 252, WGTPAW, ZERO)
  CALL HOPERA (262, '+', 262, 262, WGTPAW, ZERO)
  CALL HOPERA (272, '+', 272, 272, WGTPAW, ZERO)
  CALL HOPERA (282, '+', 282, 282, WGTPAW, ZERO)
  CALL HOPERA (292, '+', 292, 292, WGTPAW, ZERO)
  CALL HOPERA (297, '+', 297, 297, WGTPAW, ZERO)

  CALL HOPERA (600, '+', 600, 600, WGTPAW, ZERO)
  CALL HOPERA (610, '+', 610, 610, WGTPAW, ZERO)
  CALL HOPERA (620, '+', 620, 620, WGTPAW, ZERO)
  CALL HOPERA (630, '+', 630, 630, WGTPAW, ZERO)

  CALL HOPERA (601, '+', 601, 601, WGTPAW, ZERO)
  CALL HOPERA (611, '+', 611, 611, WGTPAW, ZERO)
  CALL HOPERA (621, '+', 621, 621, WGTPAW, ZERO)
  CALL HOPERA (631, '+', 631, 631, WGTPAW, ZERO)

  CALL HOPERA (602, '+', 602, 602, WGTPAW, ZERO)
  CALL HOPERA (612, '+', 612, 612, WGTPAW, ZERO)
  CALL HOPERA (622, '+', 622, 622, WGTPAW, ZERO)
  CALL HOPERA (632, '+', 632, 632, WGTPAW, ZERO)

  RETURN
END SUBROUTINE NORMA
!----------------------------------------------------------------------

!----------------------------------------------------------------------
SUBROUTINE INIBOOK (ETOT, NBIN)
!
! Initialise HBOOK
!
!----------------------------------------------------------------------
  use precision
  implicit none
  REAL (pr) :: ETOT
  REAL*4 UN, E2, E3, E4, U2, ZERO
  INTEGER ISTA, NBIN
  PARAMETER (UN=1.E0, U2=.5E0, ZERO=0.E0)
  INCLUDE 'paw.f90'

  E2=SNGL (ETOT/2._pr)
  E3=SNGL (ETOT/3._pr)
  E4=SNGL (ETOT/4._pr)
  CALL HLIMIT (NWPAWC)
  CALL HROPEN (1, 'DON', 'hist.dat', 'NQ', 1024, ISTA)
  CALL HCDIR ('//DON', ' ')

  CALL HBOOK1 (200, 'QED:Acolinearite photons 1 et 2', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (210, 'QED:Angle photons 1 et 2', NBIN, -UN, -U2, 0.)
  CALL HBOOK1 (220, 'QED:Angle photons 1 et 3', NBIN, -UN, ZERO, 0.)
  CALL HBOOK1 (230, 'QED:Angle photons 2 et 3', NBIN, -U2, UN, 0.)
  CALL HBOOK1 (240, 'QED:Angle photon 1 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (250, 'QED:Angle photon 2 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (260, 'QED:Angle photon 3 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (270, 'QED:Energie du photon 1', NBIN, E3, E2, 0.E0)
  CALL HBOOK1 (280, 'QED:Energie du photon 2', NBIN, E4, E2, 0.E0)
  CALL HBOOK1 (290, 'QED:Energie du photon 3', NBIN, 0.E0, E3, 0.E0)
  CALL HBOOK1 (295, 'QED: Angle normale et l axe', NBIN, 0.E0, UN, 0.)

  CALL HBOOK1 (201, 'B+ :Acolinearite photons 1 et 2', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (211, 'B+ :Angle photons 1 et 2', NBIN, -UN, -U2, 0.)
  CALL HBOOK1 (221, 'B+ :Angle photons 1 et 3', NBIN, -UN, ZERO, 0.)
  CALL HBOOK1 (231, 'B+ :Angle photons 2 et 3', NBIN, -U2, UN, 0.)
  CALL HBOOK1 (241, 'B+ :Angle photon 1 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (251, 'B+ :Angle photon 2 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (261, 'B+ :Angle photon 3 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (271, 'B+ :Energie du photon 1', NBIN, E3, E2, 0.E0)
  CALL HBOOK1 (281, 'B+ :Energie du photon 2', NBIN, E4, E2, 0.E0)
  CALL HBOOK1 (291, 'B+ :Energie du photon 3', NBIN, ZERO, E3, 0.E0)
  CALL HBOOK1 (296, 'B+ : Angle normale et l axe', NBIN, ZERO, UN, 0.)

  CALL HBOOK1 (202, 'B- :Acolinearite photons 1 et 2', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (212, 'B- :Angle photons 1 et 2', NBIN, -UN, -U2, 0.)
  CALL HBOOK1 (222, 'B- :Angle photons 1 et 3', NBIN, -UN, ZERO, 0.)
  CALL HBOOK1 (232, 'B- :Angle photons 2 et 3', NBIN, -U2, UN, 0.)
  CALL HBOOK1 (242, 'B- :Angle photon 1 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (252, 'B- :Angle photon 2 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (262, 'B- :Angle photon 3 et l axe', NBIN, -UN, UN, 0.)
  CALL HBOOK1 (272, 'B- :Energie du photon 1', NBIN, E3, E2, 0.E0)
  CALL HBOOK1 (282, 'B- :Energie du photon 2', NBIN, E4, E2, 0.E0)
  CALL HBOOK1 (292, 'B- :Energie du photon 3', NBIN, ZERO, E3, 0.E0)
  CALL HBOOK1 (297, 'B- : Angle normale et l axe', NBIN, ZERO, UN, 0.)

  CALL HBOOK2 (600, 'QED: Cos (1, z)/E-', 40, -UN, UN, 40, 0.E0, E3, 0.E0)
  CALL HBOOK2 (601, 'B+ : Cos (1, z)/E-', 40, -UN, UN, 40, 0.E0, E3, 0.E0)
  CALL HBOOK2 (602, 'B- : Cos (1, z)/E-', 40, -UN, UN, 40, 0.E0, E3, 0.E0)

  CALL HBOOK2 (610, 'QED: E0/E-', 40, E4, E2, 40, 0.E0, E3, 0.E0)
  CALL HBOOK2 (611, 'B+ : E0/E-', 40, E4, E2, 40, 0.E0, E3, 0.E0)
  CALL HBOOK2 (612, 'B- : E0/E-', 40, E4, E2, 40, 0.E0, E3, 0.E0)

  CALL HBOOK2 (620, 'QED: Cos (1, n)/E-', 40, ZERO, UN, 40, 0.E0, E3, 0.E0)
  CALL HBOOK2 (621, 'B+ : Cos (1, n)/E-', 40, ZERO, UN, 40, 0.E0, E3, 0.E0)
  CALL HBOOK2 (622, 'B- : Cos (1, n)/E-', 40, ZERO, UN, 40, 0.E0, E3, 0.E0)

  CALL HBOOK2 (630, 'QED: Cos (1, n)/E+', 40, ZERO, UN, 40, E3, E2, 0.E0)
  CALL HBOOK2 (631, 'B+ : Cos (1, n)/E+', 40, ZERO, UN, 40, E3, E2, 0.E0)
  CALL HBOOK2 (632, 'B- : Cos (1, n)/E+', 40, ZERO, UN, 40, E3, E2, 0.E0)

  RETURN
END SUBROUTINE INIBOOK
!----------------------------------------------------------------------

!----------------------------------------------------------------------
SUBROUTINE NORMSUP (ETOT, C0, C1, C2)
!----------------------------------------------------------------------
!
! Effectue la normalisation supplementaire des histogrammes 
!
!----------------------------------------------------------------------
  use precision
  implicit none
  REAL (pr) :: ETOT, C0, C1, C2
  REAL*4 W0, W1, W2, ZERO, UN, D, T
  REAL*4 E2, E3, E4, D1, D2, D3, Q
  PARAMETER (ZERO=0.E0, UN=1.E0, D=2.E0, T=3.E0)
  INCLUDE 'paw.f90'

  E2=SNGL (ETOT/2._pr)
  E3=SNGL (ETOT/3._pr)
  E4=SNGL (ETOT/4._pr)
  D1=E2-E3
  D2=E2-E4
  D3=E3
  Q=SNGL (40._pr**2/100._pr)

  W0=SNGL (C0)
  W1=SNGL (C1)
  W2=SNGL (C2)

  CALL HOPERA (200, '+', 200, 200, W0/D, ZERO)
  CALL HOPERA (210, '+', 210, 210, W0*D, ZERO)
  CALL HOPERA (220, '+', 220, 220, W0  , ZERO)
  CALL HOPERA (230, '+', 230, 230, W0*D/T, ZERO)
  CALL HOPERA (240, '+', 240, 240, W0/D, ZERO)
  CALL HOPERA (250, '+', 250, 250, W0/D, ZERO)
  CALL HOPERA (260, '+', 260, 260, W0/D, ZERO)
  CALL HOPERA (270, '+', 270, 270, W0/D1, ZERO)
  CALL HOPERA (280, '+', 280, 280, W0/D2, ZERO)
  CALL HOPERA (290, '+', 290, 290, W0/D3, ZERO)
  CALL HOPERA (295, '+', 295, 295, W0, ZERO)

  CALL HOPERA (201, '+', 201, 201, W1/D, ZERO)
  CALL HOPERA (211, '+', 211, 211, W1*D, ZERO)
  CALL HOPERA (221, '+', 221, 221, W1  , ZERO)
  CALL HOPERA (231, '+', 231, 231, W1*D/T, ZERO)
  CALL HOPERA (241, '+', 241, 241, W1/D, ZERO)
  CALL HOPERA (251, '+', 251, 251, W1/D, ZERO)
  CALL HOPERA (261, '+', 261, 261, W1/D, ZERO)
  CALL HOPERA (271, '+', 271, 271, W1/D1, ZERO)
  CALL HOPERA (281, '+', 281, 281, W1/D2, ZERO)
  CALL HOPERA (291, '+', 291, 291, W1/D3, ZERO)
  CALL HOPERA (296, '+', 296, 296, W1, ZERO)

  CALL HOPERA (202, '+', 202, 202, W2/D, ZERO)
  CALL HOPERA (212, '+', 212, 212, W2*D, ZERO)
  CALL HOPERA (222, '+', 222, 222, W2  , ZERO)
  CALL HOPERA (232, '+', 232, 232, W2*D/T, ZERO)
  CALL HOPERA (242, '+', 242, 242, W2/D, ZERO)
  CALL HOPERA (252, '+', 252, 252, W2/D, ZERO)
  CALL HOPERA (262, '+', 262, 262, W2/D, ZERO)
  CALL HOPERA (272, '+', 272, 272, W2/D1, ZERO)
  CALL HOPERA (282, '+', 282, 282, W2/D2, ZERO)
  CALL HOPERA (292, '+', 292, 292, W2/D3, ZERO)
  CALL HOPERA (297, '+', 297, 297, W2, ZERO)

  CALL HOPERA (600, '+', 600, 600, W0/D/D3*Q, ZERO)
  CALL HOPERA (610, '+', 610, 610, W0/D2/D3*Q, ZERO)
  CALL HOPERA (620, '+', 620, 620, W0/D3*Q, ZERO)
  CALL HOPERA (630, '+', 630, 630, W0/D1*Q, ZERO)

  CALL HOPERA (601, '+', 601, 601, W1/D/D3*Q, ZERO)
  CALL HOPERA (611, '+', 611, 611, W1/D2/D3*Q, ZERO)
  CALL HOPERA (621, '+', 621, 621, W0/D3*Q, ZERO)
  CALL HOPERA (631, '+', 631, 631, W0/D1*Q, ZERO)

  CALL HOPERA (602, '+', 602, 602, W2/D/D3*Q, ZERO)
  CALL HOPERA (612, '+', 612, 612, W2/D2/D3*Q, ZERO)
  CALL HOPERA (622, '+', 622, 622, W0/D3*Q, ZERO)
  CALL HOPERA (632, '+', 632, 632, W0/D1*Q, ZERO)

  RETURN
END SUBROUTINE NORMSUP
!----------------------------------------------------------------------
