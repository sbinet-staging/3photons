!
!   Common pour passer les angles
!
!> \brief set of angles between outgoing particles, beam and plane defined by the outgoing particles
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module angle
  use precision
  implicit none
  REAL (pr) :: &
       COS1P1, & !< cosine of angle between highest energy photon and electron beam
       COS2P1, & !< cosine of angle between middle energy photon and electron beam
       COS3P1, & !< cosine of angle between lowest energy photon and electron beam
       COS12,  & !< cosine of angle between highest energy photon and middle energy photon
       COS13,  & !< cosine of angle between highest energy photon and lowest energy photon
       COS23,  & !< cosine of angle between middle energy photon and lowest energy photon
       COSN,   & !< cosine of angle between outgoing photons plane and electron beam
       COSAC     !< cosine of angle between outgoing photons perpendicular momenta
  !> \brief type replacing the previous common
  type angle_typ
     real (pr) :: &
          cos12,  & !< cosine of angle between highest energy photon and electron beam
          cos1p1, & !< cosine of angle between middle energy photon and electron beam
          cos2p1, & !< cosine of angle between lowest energy photon and electron beam
          cos3p1, & !< cosine of angle between highest energy photon and middle energy photon
          cos13,  & !< cosine of angle between highest energy photon and lowest energy photon
          cos23,  & !< cosine of angle between middle energy photon and lowest energy photon
          cosn,   & !< cosine of angle between outgoing photons plane and electron beam
          cosac     !< cosine of angle between outgoing photons perpendicular momenta
  end type angle_typ
end module angle
