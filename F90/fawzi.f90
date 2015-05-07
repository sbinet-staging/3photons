!----------------------------------------------------------------------
SUBROUTINE FAWZI (UNE, BREPEM, CONVERS, PI, ETOT)
!----------------------------------------------------------------------
!
! Affiche les resultats analytiques de Fawzi, 
! et compare aux resultats du MC
!
!----------------------------------------------------------------------
  use precision
  use cutpar
  use param
  use ppp
  use resfin
  implicit none
  INTEGER UNE
  REAL (pr) :: BREPEM, CONVERS, PI, ETOT
  REAL (pr) :: BRA, MRE, GRE, SIG, DEL, EPS
  REAL (pr) :: F1, G1, G2, G3, FF, GG
  REAL (pr) :: SIGP, SIGM, MCP, MCM, INCRP, INCRM
  COMPLEX (pr) :: SDZ

  MRE=MZ0/ETOT
  GRE=GZ0*MZ0/ETOT**2
  SDZ=CMPLX(1._pr-MRE**2, -GRE) &
   /((1._pr-MRE**2)**2+GRE**2)
  DEL=(1._pr-BCUT)/2._pr
  EPS=2._pr*EMIN/ETOT
  BRA=MZ0/3._pr/6._pr/PI**3/16._pr/120._pr
  !SIG=12._pr*PI/MZ0**2*(BREPEM*GZ0)*BRA/ETOT**2*(ETOT/MZ0)**8 &
  !       *CDABS(SDZ)**2*CONVERS
  SIG=12._pr*PI/MZ0**2*(BREPEM*GZ0)*BRA/ETOT**2*(ETOT/MZ0)**8 * ABS (SDZ)**2*CONVERS
  
  F1=       1._pr- 15._pr*EPS**4 &
   - 9._pr/ 7._pr*(1._pr- 70._pr*EPS**4)*DEL**2 &
   + 6._pr/ 7._pr*(1._pr+ 70._pr*EPS**4)*DEL**3

  G1=       1._pr- 30._pr*EPS**4 &
   - 9._pr/ 7._pr*(1._pr- 70._pr*EPS**4)*DEL &
               - 90._pr*EPS**4 *DEL**2 &
   - 1._pr/ 7._pr*(1._pr-420._pr*EPS**4)*DEL**3

  G2=       1._pr- 25._pr*EPS**4 &
   - 6._pr/ 7._pr*(1._pr- 70._pr*EPS**4)*DEL &
   - 3._pr/ 7._pr*(1._pr+210._pr*EPS**4)*DEL**2 &
   - 8._pr/21._pr*(1._pr-105._pr/2._pr*EPS**4)*DEL**3

  G3=       1._pr-195._pr/11._pr*EPS**4 &
   -18._pr/77._pr*(1._pr- 7._pr*EPS**4)*DEL &
   - 9._pr/11._pr*(9._pr/7._pr-70._pr*EPS**4)*DEL**2 &
   - 8._pr/11._pr*(1._pr-105._pr/11._pr*EPS**4)*DEL**3

  FF=F1*(1._pr-SINCUT**3)
  GG=G1-27._pr/16._pr*G2*SINCUT+11._pr/16._pr*G3*SINCUT**3

  SIGP=SIG*(FF+2._pr*GG)
  SIGM=SIG*(FF+4._pr*GG)

  MCP=(SPM2(1, 2)+SPM2(2, 2))/4._pr
  MCM=(SPM2(1, 3)+SPM2(2, 3))/4._pr
  INCRP=SQRT((SPM2(1, 2)*VAR(1, 2))**2+(SPM2(2, 2)*VAR(2, 2))**2) &
        /ABS(SPM2(1, 2)+SPM2(2, 2))
  INCRM=SQRT((SPM2(1, 3)*VAR(1, 3))**2+(SPM2(2, 3)*VAR(2, 3))**2) &
        /ABS(SPM2(1, 3)+SPM2(2, 3))

  WRITE(UNE, *)
  WRITE(UNE, 6)'s (pb) :   Sig_cut_Th    Sig_Th      Rapport'
  WRITE(UNE, 6)'       :   Sig_Num'
  WRITE(UNE, 6)'       :   Ecart_relatif  Incertitude'
  WRITE(UNE, *)
  WRITE(UNE, 7)'s+(pb) : ', SIGP, SIG*3._pr, SIGP/SIG/3._pr
  WRITE(UNE, 8)'       : ', MCP
  WRITE(UNE, 9)'       : ', MCP/SIGP-1._pr, INCRP, (MCP/SIGP-1._pr)/INCRP
  WRITE(UNE, *)
  WRITE(UNE, 7)'s-(pb) : ', SIGM, SIG*5._pr, SIGM/SIG/5._pr
  WRITE(UNE, 8)'       : ', MCM
  WRITE(UNE, 9)'       : ', MCM/SIGM-1._pr, INCRM, (MCM/SIGM-1._pr)/INCRM
  WRITE(UNE, *)

  RETURN
6 FORMAT(A)
7 FORMAT(A, 2G14.8, G10.4)
8 FORMAT(A,  G14.8)
9 FORMAT(A, 2G14.8, G10.4)
END SUBROUTINE FAWZI
!----------------------------------------------------------------------
