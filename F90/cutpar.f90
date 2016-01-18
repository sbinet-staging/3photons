!
!   Common pour passer les coupures
!
!> \brief set of experimental cuts
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module cutpar
  use precision
  implicit none
  REAL (pr) :: &
       ACUT,   & !< cut on maximum cosine of (beam, photons) angle
       BCUT,   & !< cut on maximum cosine of (photon, photon) angle
       EMIN,   & !< cut on minimum photon energy
       SINCUT    !< cut on minimum cosine of (beam, normal to the photon plane) angle
  !> \brief type replacing the previous common
  type cutpar_typ
     real (pr) :: &
          acut,   & !< cut on maximum cosine of (beam, photons) angle                   
          bcut,   & !< cut on maximum cosine of (photon, photon) angle                  
          emin,   & !< cut on minimum photon energy                                     
          sincut    !< cut on minimum cosine of (beam, normal to the photon plane) angle
  end type cutpar_typ
end module cutpar
