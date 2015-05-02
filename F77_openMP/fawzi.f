C----------------------------------------------------------------------
      SUBROUTINE FAWZI(UNE,BREPEM,CONVERS,PI,ETOT)
C----------------------------------------------------------------------
C
C Affiche les resultats analytiques de Fawzi, 
C et compare aux resultats du MC
C
C----------------------------------------------------------------------
!     $ use OMP_LIB
      IMPLICIT NONE
      INTEGER UNE
      double precision BREPEM,CONVERS,PI,ETOT
      INCLUDE 'cutpar.for'
      INCLUDE 'param.for'
      INCLUDE 'ppp.for'
      INCLUDE 'resfin.for'
      double precision BRA,MRE,GRE,SIG,DEL,EPS
      double precision F1,G1,G2,G3,FF,GG
      double precision SIGP,SIGM,MCP,MCM,INCRP,INCRM
      double complex SDZ

      MRE=MZ0/ETOT
      GRE=GZ0*MZ0/ETOT**2
      SDZ=CMPLX(1.d0-MRE**2,-GRE)
     &/((1.d0-MRE**2)**2+GRE**2)
      DEL=(1.d0-BCUT)/2.d0
      EPS=2.d0*EMIN/ETOT
      BRA=MZ0/3.d0/6.d0/PI**3/16.d0/120.d0
      SIG=12.d0*PI/MZ0**2*(BREPEM*GZ0)*BRA/ETOT**2*(ETOT/MZ0)**8
     &      *CDABS(SDZ)**2*CONVERS

      F1=           1.d0- 15.d0*EPS**4
     &- 9.d0/ 7.d0*(1.d0- 70.d0*EPS**4)*DEL**2
     &+ 6.d0/ 7.d0*(1.d0+ 70.d0*EPS**4)*DEL**3

      G1=           1.d0- 30.d0*EPS**4
     &- 9.d0/ 7.d0*(1.d0- 70.d0*EPS**4)*DEL
     &            - 90.d0*EPS**4 *DEL**2
     &- 1.d0/ 7.d0*(1.d0-420.d0*EPS**4)*DEL**3

      G2=           1.d0- 25.d0*EPS**4
     &- 6.d0/ 7.d0*(1.d0- 70.d0*EPS**4)*DEL
     &- 3.d0/ 7.d0*(1.d0+210.d0*EPS**4)*DEL**2
     &- 8.d0/21.d0*(1.d0-105.d0/2.d0*EPS**4)*DEL**3

      G3=           1.d0-195.d0/11.d0*EPS**4
     &-18.d0/77.d0*(1.d0-        7.d0*EPS**4)*DEL
     &- 9.d0/11.d0*(9.d0/7.d0-  70.d0*EPS**4)*DEL**2
     &- 8.d0/11.d0*(1.d0-105.d0/11.d0*EPS**4)*DEL**3

      FF=F1*(1.d0-SINCUT**3)
      GG=G1-27.d0/16.d0*G2*SINCUT+11.d0/16.d0*G3*SINCUT**3

      SIGP=SIG*(FF+2.d0*GG)
      SIGM=SIG*(FF+4.d0*GG)

      MCP = (SPM2 (1, 2) + SPM2 (2, 2))/4.d0
      MCM = (SPM2 (1, 3) + SPM2 (2, 3))/4.d0
      INCRP=SQRT((SPM2(1,2)*VAR(1,2))**2+(SPM2(2,2)*VAR(2,2))**2)
     &     /ABS(SPM2(1,2)+SPM2(2,2))
      INCRM=SQRT((SPM2(1,3)*VAR(1,3))**2+(SPM2(2,3)*VAR(2,3))**2)
     &     /ABS(SPM2(1,3)+SPM2(2,3))

      WRITE(UNE,*)
      WRITE(UNE,6)'s (pb) :   Sig_cut_Th    Sig_Th      Rapport'
      WRITE(UNE,6)'       :   Sig_Num'
      WRITE(UNE,6)'       :   Ecart_relatif  Incertitude'
      WRITE(UNE,*)
      WRITE(UNE,7)'s+(pb) : ',SIGP,SIG*3.d0,SIGP/SIG/3.d0
      WRITE(UNE,8)'       : ',MCP
      WRITE(UNE,9)'       : ',MCP/SIGP-1.d0,INCRP,(MCP/SIGP-1.d0)/INCRP
      WRITE(UNE,*)
      WRITE(UNE,7)'s-(pb) : ',SIGM,SIG*5.d0,SIGM/SIG/5.d0
      WRITE(UNE,8)'       : ',MCM
      WRITE(UNE,9)'       : ',MCM/SIGM-1.d0,INCRM,(MCM/SIGM-1.d0)/INCRM
      WRITE(UNE,*)

      RETURN
 6    FORMAT(A)
 7    FORMAT(A,2G14.8,G10.4)
 8    FORMAT(A, G14.8)
 9    FORMAT(A,2G14.8,G10.4)
      END
C----------------------------------------------------------------------
c      SIGP=SIG*(3.-75.*EPS**4-18./7.*DEL
c     & -9./7.*DEL**2 -8./7.*DEL**3-27./8.*SINCUT)
c      SIGM=SIG*(5.-135.*EPS**4-36./7.*DEL
c     & -9./7.*DEL**2 -10./7.*DEL**3-27./4.*SINCUT)
c      double precision SCUT
c      SCUT=SQRT(1.-ACUT**2)
c      WRITE(UNE,*) 's+(pb) : ',SIG*3.*
c     &           (1.-SCUT*(1.+  ACUT**2/ 8.))
c      WRITE(UNE,*) 's-(pb) : ',SIG*5.*
c     &           (1.-SCUT*(1.+7*ACUT**2/20.))
