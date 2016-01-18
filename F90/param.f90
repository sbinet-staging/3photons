!
!   Common pour passer les parametres physiques
!
!> \brief array of generated momenta
!> \author Vincent C. LAFAGE
!> \date 2007-06-13 ISO
module param
  use precision
  implicit none
  REAL (pr) :: &
       MZ0,    & !< Z⁰ boson mass (GeV)
       GZ0,    & !< Z⁰ boson width (GeV)
       GA,     & !< Standard Model contribution electromagnetic coupling √(4𝜋𝛼)³
       GBP,    & !< 𝛽₊ anomalous contribution electroweak coupling
       GBM,    & !< 𝛽₋ anomalous contribution electroweak coupling
       POLP,   & !< electroweak polarisations factors for 𝛽₊ anomalous contribution
       POLM,   & !< electroweak polarisations factors for 𝛽₋ anomalous contribution
       POLP2,  & !< electroweak polarisations factors for 𝛽₊ anomalous contribution
       POLM2     !< electroweak polarisations factors for 𝛽₋ anomalous contribution
  LOGICAL IMPR   !< boolean predicate value controlling dump of result
  !> \brief type replacing the previous common
  type param_typ
     real (pr) :: &
          mz0,    & !< Z⁰ boson mass (GeV)                                            
          gz0,    & !< Z⁰ boson width (GeV)                                           
          ga,     & !< Standard Model contribution electromagnetic coupling √(4𝜋𝛼)³   
          gbp,    & !< 𝛽₊ anomalous contribution electroweak coupling                 
          gbm,    & !< 𝛽₋ anomalous contribution electroweak coupling                 
          polp,   & !< electroweak polarisations factors for 𝛽₊ anomalous contribution
          polm,   & !< electroweak polarisations factors for 𝛽₋ anomalous contribution
          polp2,  & !< electroweak polarisations factors for 𝛽₊ anomalous contribution
          polm2     !< electroweak polarisations factors for 𝛽₋ anomalous contribution
     logical impr   !< boolean predicate value controlling dump of result
  end type param_typ
end module param
