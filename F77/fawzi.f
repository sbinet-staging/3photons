C----------------------------------------------------------------------
      SUBROUTINE FAWZI(UNE,BREPEM,CONVERS,PI,ETOT)
C----------------------------------------------------------------------
C
C Affiche les resultats analytiques de Fawzi, 
C et compare aux resultats du MC
C
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER UNE
      REAL*8 BREPEM,CONVERS,PI,ETOT
      INCLUDE 'cutpar.i'
      INCLUDE 'param.i'
      INCLUDE 'ppp.i'
      INCLUDE 'resfin.i'
      REAL*8 BRA,MRE,GRE,SIG,DEL,EPS
      REAL*8 F1,G1,G2,G3,FF,GG
      REAL*8 SIGP,SIGM,MCP,MCM,INCRP,INCRM
      COMPLEX*16 SDZ

      MRE=MZ0/ETOT
      GRE=GZ0*MZ0/ETOT**2
      SDZ=CMPLX(1.-MRE**2,-GRE)
     &/((1.-MRE**2)**2+GRE**2)
      DEL=(1.-BCUT)/2.
      EPS=2.*EMIN/ETOT
      BRA=MZ0/3./6./PI**3/16./120.
      SIG=12.*PI/MZ0**2*(BREPEM*GZ0)*BRA/ETOT**2*(ETOT/MZ0)**8
     &      *CDABS(SDZ)**2*CONVERS

      F1=       1.- 15.*EPS**4
     &- 9./ 7.*(1.- 70.*EPS**4)*DEL**2
     &+ 6./ 7.*(1.+ 70.*EPS**4)*DEL**3

      G1=       1.- 30.*EPS**4
     &- 9./ 7.*(1.- 70.*EPS**4)*DEL
     &            - 90.*EPS**4 *DEL**2
     &- 1./ 7.*(1.-420.*EPS**4)*DEL**3

      G2=       1.- 25.*EPS**4
     &- 6./ 7.*(1.- 70.*EPS**4)*DEL
     &- 3./ 7.*(1.+210.*EPS**4)*DEL**2
     &- 8./21.*(1.-105./2.*EPS**4)*DEL**3

      G3=       1.-195./11.*EPS**4
     &-18./77.*(1.- 7.*EPS**4)*DEL
     &- 9./11.*(9./7.-70.*EPS**4)*DEL**2
     &- 8./11.*(1.-105./11.*EPS**4)*DEL**3

      FF=F1*(1.-SINCUT**3)
      GG=G1-27./16.*G2*SINCUT+11./16.*G3*SINCUT**3

      SIGP=SIG*(FF+2.*GG)
      SIGM=SIG*(FF+4.*GG)

      MCP=(SPM2(1,2)+SPM2(2,2))/4.D0
      MCM=(SPM2(1,3)+SPM2(2,3))/4.D0
      INCRP=SQRT((SPM2(1,2)*VAR(1,2))**2+(SPM2(2,2)*VAR(2,2))**2)
     &     /ABS(SPM2(1,2)+SPM2(2,2))
      INCRM=SQRT((SPM2(1,3)*VAR(1,3))**2+(SPM2(2,3)*VAR(2,3))**2)
     &     /ABS(SPM2(1,3)+SPM2(2,3))

      WRITE(UNE,*)
      WRITE(UNE,6)'s (pb) :   Sig_cut_Th    Sig_Th      Rapport'
      WRITE(UNE,6)'       :   Sig_Num'
      WRITE(UNE,6)'       :   Ecart_relatif  Incertitude'
      WRITE(UNE,*)
      WRITE(UNE,7)'s+(pb) : ',SIGP,SIG*3.,SIGP/SIG/3.
      WRITE(UNE,8)'       : ',MCP
      WRITE(UNE,9)'       : ',MCP/SIGP-1.,INCRP,(MCP/SIGP-1.)/INCRP
      WRITE(UNE,*)
      WRITE(UNE,7)'s-(pb) : ',SIGM,SIG*5.,SIGM/SIG/5.
      WRITE(UNE,8)'       : ',MCM
      WRITE(UNE,9)'       : ',MCM/SIGM-1.,INCRM,(MCM/SIGM-1.)/INCRM
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
c      REAL*8 SCUT
c      SCUT=SQRT(1.-ACUT**2)
c      WRITE(UNE,*) 's+(pb) : ',SIG*3.*
c     &           (1.-SCUT*(1.+  ACUT**2/ 8.))
c      WRITE(UNE,*) 's-(pb) : ',SIG*5.*
c     &           (1.-SCUT*(1.+7*ACUT**2/20.))
