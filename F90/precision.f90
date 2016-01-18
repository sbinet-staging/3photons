!> \brief parameter that set the overall precision for the program
!> \author Vincent C. LAFAGE
!> \date 2007-06-15 ISO
module precision
  implicit none
!  integer, parameter :: pr = selected_real_kind (15,300)
!  integer, parameter :: pr = selected_real_kind (30,3000)
!  integer, parameter :: pr = 4
  integer, parameter :: pr = 8 !< kind of real used in computation throughout the program
!  integer, parameter :: pr = 16
end module precision
