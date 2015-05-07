!--------------------------------------------------------------------
PROGRAM RIEMANN
!--------------------------------------------------------------------
!$ use omp_lib
  use precision
  use ppp
  use spinor
  use scalar
  use angle
  use result
  use cutpar
  use resfin
  use param
  use xtrpro
  implicit none 
  REAL (pr) :: ETOT, FLUX, NORM, POIDS, WTEV, SIGMA, PREC, VARIANCE
  REAL (pr) :: POIDS2, PROBA
  REAL (pr) :: MFIN(100)
  INTEGER ITOT, NTOT, NBIN
  INTEGER I, K
  INTEGER Sp, L1, L2, L3
  REAL (pr) :: ALPHA, CONVERS, SIN2W, COS2W, FACTCOM, E2, ALPHAZ, E2Z, GZR
  REAL (pr) :: BETAPLUS, BETAMOINS, BETAMIN
  REAL (pr) :: SSP, SSM, INCSSP, INCSSM
  REAL (pr) :: PAA, PAB, PBB
  REAL (pr) :: CAA, CAB, CBB
  REAL (pr) :: DZETA, PROPAG, ECARTPIC
  REAL (pr) :: BREPEM
  LOGICAL CUT, PLOT
  INTEGER ICYCLE
  INCLUDE 'paw.f90'
  CHARACTER*80 A
  INTEGER UNL, UNE
  PARAMETER (UNL=10, UNE=20)
  REAL*4 etime, tarray(2)
  CHARACTER*8 :: TIMESTR
  CHARACTER*9 :: TODAY = ''
  REAL (pr), parameter :: PI = 4._pr * ATAN (1._pr)

!$    integer nb_taches
!$    integer rang
!$    integer iii
! $    integer OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!   LECTURE DES PARAMETRES
  OPEN(UNIT=UNL, FILE='valeurs', STATUS='OLD')
  READ(UNL, *) ITOT, A
  READ(UNL, *) ETOT, A
  READ(UNL, *) ACUT, A
  READ(UNL, *) BCUT, A
  READ(UNL, *) EMIN, A
  READ(UNL, *) SINCUT, A
  READ(UNL, *) ALPHA, A
  READ(UNL, *) ALPHAZ, A
  READ(UNL, *) CONVERS, A
  READ(UNL, *) MZ0, A
  READ(UNL, *) GZ0, A
  READ(UNL, *) SIN2W, A
  READ(UNL, *) BREPEM, A
  READ(UNL, *) BETAPLUS, A
  READ(UNL, *) BETAMOINS, A
  READ(UNL, *) NBIN, A
  READ(UNL, *) IMPR, A
  READ(UNL, *) PLOT, A
  CLOSE(UNIT=UNL)
!      WRITE(*, *) ETOT
!   INITIALISATION DE PAW
!  CALL INIBOOK(ETOT, NBIN)
!   MASSES DES PARTICULES DANS L'ETAT FINAL
  MFIN(1)=0._pr
  MFIN(2)=0._pr
  MFIN(3)=0._pr
!   CALCUL DU FACTEUR DE FLUX (=1/2s pour 2 particules initiales sans masse)
  FLUX=1._pr/(2._pr*ETOT**2)
!   INCLUE LA NORMALISATION TOTALE DE L'ESPACE DES PHASES
!  PI=4._pr*ATAN(1._pr)
  NORM=((2._pr*PI)**(4-3*INP))/DBLE(ITOT)
!   FACTEUR COMMUN NON-MOYENNE SUR LES SPINS=1
!                   /(FACTEUR DE SYMETRIE)
!                   *(FACTEUR DE CONVERSION GeV^-2->pb)
!   POUR MOYENNER SUR LES SPINS, RAJOUTER :
!                   /(NOMBRE D'HELICITES INCIDENTES)
  FACTCOM =1._pr/6._pr*CONVERS
  E2      =4._pr*PI*ALPHA
  E2Z     =4._pr*PI*ALPHAZ
  COS2W=1._pr-SIN2W
  GZR=GZ0/MZ0
!   LE FACTEUR RAC8 VIENT DES SQRT(2) DE CHAQUE VECTEUR DE POLARISATION
!   DANS LA METHODE DES AMPLITUDES D'HELICITE
  RAC8=SQRT(8._pr) ! dorénavant, c'est un parametre... nan, pas vraiment
!   COUPLAGES
  GA =-SQRT(E2)**3
  GBP=-SQRT(E2Z/(4*COS2W*SIN2W))/MZ0**4
  GBM=GBP
!   FACTEURS DU A LA SOMME SUR LES POLARISATIONS
  POLP=-2*SIN2W
  POLM=1-2*SIN2W
  POLP2=POLP**2
  POLM2=POLM**2
  PAA=2
  PAB=1-4*SIN2W
  PBB=1-4*SIN2W+8*SIN2W**2
!   COEFFICIENT POUR HOMOGENEISER
  CAA=FACTCOM*PAA
  CAB=FACTCOM*PAB/MZ0**2
  CBB=FACTCOM*PBB/MZ0**4
!   POIDS D'UN EVENEMENT D'ESPACE DES PHASES A 3 PHOTONS
  WTEV=PI**2/8._pr*ETOT**2
!   PASSAGE EN VARIABLE A-DIMENSIONEE...
  DZETA=(ETOT/MZ0)**2
  ECARTPIC=(DZETA-1.)/GZR
  PROPAG=1./(ECARTPIC**2+1.0)
!   INITIALISATION DES IMPULSIONS ENTRANTES
  P1(4)=ETOT/2._pr
  P2(4)=ETOT/2._pr
  P1(1)=-ETOT/2._pr
  P2(1)=ETOT/2._pr
  P1(2)=0._pr
  P1(3)=0._pr
  P2(2)=0._pr
  P2(3)=0._pr

!   INITIALISATION DES CUMULANTS
  NTOT=0
  DO Sp=1, 2
     DO K=1, NRESUL
        SPM2(Sp, K)=0._pr
        VAR(Sp, K)=0._pr
     ENDDO
  ENDDO
  SIGMA=0._pr
  VARIANCE=0._pr

  !$OMP PARALLEL  DEFAULT(FIRSTPRIVATE) &
  !$OMP SHARED(SIGMA, VARIANCE, NTOT, SPM2DIF, SPM2, VAR)
  ! $ OMP COPYIN(cos12, cos1p1, cos2p1, cos3p1, cos13, cos23, cosn, cosac, acut, bcut, emin, sincut, &
  ! $ OMP mz0, gz0, ga, gbp, gbm, polp, polm, polp2, polm2, s, t, rac8, pout, p1, p2)
  ! $ OMP COPYIN(/angle/, /cutpar/, /param/, /spinor/, /PAWC/, /ppp/)
  !$     nb_taches = OMP_GET_NUM_THREADS()
  !$     rang = OMP_GET_THREAD_NUM ()
  !$     write (*, *) 'Nbtaches', nb_taches, rang
  !$OMP DO REDUCTION (+:SIGMA, VARIANCE, NTOT, SPM2, VAR)
!   DEBUT DE LA BOUCLE D'INTEGRATION
  DO I=1, ITOT

     CALL RAMBO(INP, ETOT, MFIN, POUT, WTEV)
     !write (*, '(Z16)') POUT
     WTEV=WTEV*NORM
!   TRIE LES PHOTONS SORTANTS PAR ENERGIE
     CALL TRI
!   CALCUL DES PRODUITS SPINORIELS, PRODUITS SCALAIRES ET ANGLES
     CALL PRECALC
!   CALCULE LE POIDS TOTAL DE L'EVENEMENT AVEC L'ELEMENT DE MATRICE
!   (COUPURES S'IL Y A LIEU)
     IF (.NOT.CUT(0.)) THEN
        CALL MATRX(ETOT)
        DO K=1, NRESUL
           SPM2DIF(K)=0._pr
           DO L1=1, 2
              DO L2=1, 2
                 DO L3=1, 2
                    SPM2DIF(K)=SPM2DIF(K)+M2(L1, L2, L3, K)
                 ENDDO
              ENDDO
           ENDDO
           SPM2(1, K)=SPM2(1, K)+SPM2DIF(K)
           VAR(1, K) =VAR(1, K) +SPM2DIF(K)**2
        ENDDO
        PROBA= &
             CAA*SPM2DIF(1) &
             +CBB*(BETAPLUS**2*SPM2DIF(2) &
             +BETAMOINS**2*SPM2DIF(3)) &
             /GZR**2*PROPAG &
             +CAB*2*BETAPLUS*(ECARTPIC*SPM2DIF(4) &
             -SPM2DIF(5)) &
             /GZR*PROPAG
        POIDS=PROBA*WTEV/4._pr
        SIGMA=SIGMA+POIDS
        VARIANCE=VARIANCE+POIDS**2
!    STOCKE DANS UN HISTOGRAMME LES PARAMETRES DE L'EVENEMENT
!        IF(PLOT) THEN
!           POIDS2=CAA*SPM2DIF(1)*WTEV/4._pr
!           PRPLUS =SNGL(CBB*BETAPLUS**2*SPM2DIF(2) &
!                /GZR**2*PROPAG &
!                *WTEV/4._pr)
!           PRMOINS=SNGL(CBB*BETAMOINS**2*SPM2DIF(3) &
!                /GZR**2*PROPAG &
!                *WTEV/4._pr)
!           CALL BOOK(POIDS, POIDS2)
!        ENDIF
        NTOT=NTOT+1
     ELSE
        POIDS=0._pr
     ENDIF

  ENDDO
!$OMP END DO
!$OMP END PARALLEL
! $ OMP END PARALLEL DO

! FIN DE LA BOUCLE D'ECHANTILLONNAGE

! CALCUL DES INCERTITUDES RELATIVES
  DO K=1, NRESUL
     VAR(1, K)=(VAR(1, K)-SPM2(1, K)**2/DBLE(ITOT))/DBLE(ITOT-1)
     VAR(1, K)=SQRT(VAR(1, K)/DBLE(ITOT))/ABS(SPM2(1, K) &
          /DBLE(ITOT))
  ENDDO
! COPIE POUR LES SPINS OPPOSEES
  DO K=1, NRESUL
     SPM2(2, K)=SPM2(1, K)
     VAR(2, K)=VAR(1, K)
  ENDDO
! POLARISATIONS
  DO K=2, 3
     SPM2(1, K)=SPM2(1, K)*POLM2
     SPM2(2, K)=SPM2(2, K)*POLP2
  ENDDO
  DO K=4, 5
     SPM2(1, K)=SPM2(1, K)*POLM
     SPM2(2, K)=SPM2(2, K)*POLP
  ENDDO
! COEFFICIENTS PHYSIQUES et PROPAGATEUR DU Z0
  DO Sp=1, 2
     DO K=1, NRESUL
        SPM2(Sp, K)=SPM2(Sp, K)*FACTCOM*FLUX*WTEV
     ENDDO
     SPM2(Sp, 1)=SPM2(Sp, 1)
     SPM2(Sp, 2)=SPM2(Sp, 2)/GZR**2/MZ0**4*PROPAG
     SPM2(Sp, 3)=SPM2(Sp, 3)/GZR**2/MZ0**4*PROPAG
     SPM2(Sp, 4)=SPM2(Sp, 4)/GZR/MZ0**2*PROPAG*ECARTPIC
     SPM2(Sp, 5)=SPM2(Sp, 5)/GZR/MZ0**2*PROPAG
  ENDDO
  BETAMIN=SQRT((SPM2(1, 1)+SPM2(2, 1))/(SPM2(1, 2)+SPM2(2, 2)))
  SSP=(SPM2(1, 2)+SPM2(2, 2))/SQRT(SPM2(1, 1)+SPM2(2, 1))/2._pr
  SSM=(SPM2(1, 3)+SPM2(2, 3))/SQRT(SPM2(1, 1)+SPM2(2, 1))/2._pr
  INCSSP= &
       SQRT((SPM2(1, 2)*VAR(1, 2))**2 &
       +(SPM2(2, 2)*VAR(2, 2))**2)/ABS(SPM2(1, 2)+SPM2(2, 2)) &
       +SQRT((SPM2(1, 1)*VAR(1, 1))**2 &
       +(SPM2(2, 1)*VAR(2, 1))**2)/ABS(SPM2(1, 1)+SPM2(2, 1))/2._pr
  INCSSM= &
       SQRT((SPM2(1, 3)*VAR(1, 3))**2 &
       +(SPM2(2, 3)*VAR(2, 3))**2)/ABS(SPM2(1, 3)+SPM2(2, 3)) &
       +SQRT((SPM2(1, 1)*VAR(1, 1))**2 &
       +(SPM2(2, 1)*VAR(2, 1))**2)/ABS(SPM2(1, 1)+SPM2(2, 1))/2._pr

  VARIANCE=(VARIANCE-SIGMA**2/DBLE(ITOT))/DBLE(ITOT-1)
  PREC=SQRT(VARIANCE/DBLE(ITOT))/ABS(SIGMA/DBLE(ITOT))
  SIGMA=SIGMA*FLUX
!   NORMALISATION DES HISTOGRAMMES
!  CALL NORMA(FLUX*NBIN)
!  CALL NORMSUP(ETOT, 4./(SPM2(1, 1)+SPM2(2, 1)), &
!                    4./(SPM2(1, 2)+SPM2(2, 2)), &
!                    4./(SPM2(1, 3)+SPM2(2, 3)))
!   STOCKAGE DES RESULTATS DE PAW
!  CALL HROUT(0, ICYCLE, ' ')
!  CALL HREND('DON')
!   STOCKAGE DES RESULTATS NUMERIQUES
!  CALL DATE(TODAY)
!  CALL TIME(TIMESTR)
  OPEN(UNIT=UNE, FILE='res.dat', STATUS='UNKNOWN')
  WRITE(UNE, *) TODAY, '   ', TIMESTR
  WRITE(UNE, *) 
  WRITE(UNE, *) 'Nombre d''evenements            : ', ITOT
  WRITE(UNE, *) '... apres coupure              : ', NTOT
  WRITE(UNE, *) 'energie dans le CdM      (GeV) : ', ETOT
  WRITE(UNE, *) 'coupure / cos(photon, faisceau) : ', ACUT
  WRITE(UNE, *) 'coupure / cos(photon, photon)   : ', BCUT
  WRITE(UNE, *) 'coupure / sin(normale, faisceau): ', SINCUT
  WRITE(UNE, *) 'coupure sur l''energie    (GeV) : ', EMIN
  WRITE(UNE, *) '1/(constante de structure fine): ', 1._pr/ALPHA
  WRITE(UNE, *) '1/(structure fine au pic)      : ', 1._pr/ALPHAZ
  WRITE(UNE, *) 'facteur de conversion GeV-2/pb : ', CONVERS
  ! WRITE(UNE, *) 'Volume d''espace des phases     : ', VOLUME/NTOT
  WRITE(UNE, *) 'Masse du Z0              (GeV) : ', MZ0
  WRITE(UNE, *) 'Largeur du Z0            (GeV) : ', GZ0
  WRITE(UNE, *) 'Sinus^2 Theta Weinberg         : ', SIN2W
  WRITE(UNE, *) 'Taux de branchement Z--->e+e-  : ', BREPEM
  WRITE(UNE, *) 'Beta plus                      : ', BETAPLUS
  WRITE(UNE, *) 'Beta moins                     : ', BETAMOINS
  WRITE(UNE, *) '---------------------------------------------'
  WRITE(UNE, *) 'Section Efficace          (pb) : ', SIGMA
  WRITE(UNE, *) 'Ecart-Type                (pb) : ', SIGMA*PREC
  WRITE(UNE, *) 'Precision Relative             : ', PREC
  WRITE(UNE, *) '---------------------------------------------'
  WRITE(UNE, *) 'Beta minimum                   : ', BETAMIN
  WRITE(UNE, *) 'Stat. Significance  B+(pb-1/2) : ', SSP
  WRITE(UNE, *) 'Incert. Stat. Sign. B+(pb-1/2) : ', SSP*INCSSP
  WRITE(UNE, *) 'Stat. Significance  B-(pb-1/2) : ', SSM
  WRITE(UNE, *) 'Incert. Stat. Sign. B+(pb-1/2) : ', SSM*INCSSM
  WRITE(UNE, *) 'Temps ecoule                   : ', etime(tarray)
  WRITE(UNE, *) 'Temps ecoule utilisateur       : ', tarray(1)
  WRITE(UNE, *) 'Temps ecoule systeme           : ', tarray(2)
  WRITE(UNE, *) 'Temps ecoule par evenement     : ', &
       tarray(1)/DBLE(ITOT)
  WRITE(UNE, *)
  DO Sp=1, 2
     DO K=1, NRESUL
        WRITE(UNE, 930) Sp, K, SPM2(Sp, K), ABS(SPM2(Sp, K))*VAR(Sp, K), &
             VAR(Sp, K)
     ENDDO
     WRITE(UNE, *)
  ENDDO
  DO K=1, NRESUL
     WRITE(UNE, 940) K, (SPM2(1, K)+SPM2(2, K))/4._pr, &
          SQRT((SPM2(1, K)*VAR(1, K))**2 &
          +(SPM2(2, K)*VAR(2, K))**2)/4._pr, &
          SQRT((SPM2(1, K)*VAR(1, K))**2 &
          +(SPM2(2, K)*VAR(2, K))**2)/ABS(SPM2(1, K)+SPM2(2, K))
  ENDDO
  CALL ERIC(UNE, BREPEM, CONVERS, PI)
  CALL FAWZI(UNE, BREPEM, CONVERS, PI, ETOT)
  CLOSE(UNIT=UNE)

  OPEN(UNIT=UNE, ACCESS='APPEND', FILE='pil.mc', STATUS='unknown')
  WRITE(UNE, *) TODAY, '   ', TIMESTR
  WRITE(UNE, 960) ETOT, &
       (SPM2(1, 1)+SPM2(2, 1))/4._pr, &
       (SPM2(1, 2)+SPM2(2, 2))/4._pr*BETAPLUS**2, &
       (SPM2(1, 3)+SPM2(2, 3))/4._pr*BETAMOINS**2, &
       (SPM2(1, 4)+SPM2(2, 4))/4._pr*BETAPLUS, &
       ((SPM2(1, 1)+SPM2(2, 1)) &
       +(SPM2(1, 2)+SPM2(2, 2))*BETAPLUS**2 &
       +(SPM2(1, 3)+SPM2(2, 3))*BETAMOINS**2 &
       +2.*(SPM2(1, 4)+SPM2(2, 4))*BETAPLUS)/4._pr, &
       SIGMA
  CLOSE(UNIT=UNE)

!900 FORMAT(I5, 3E21.14)
!920 FORMAT(I2, 2I3, 8E8.2)
930 FORMAT(2I3, 3E15.7)
940 FORMAT(3X, I3, 3E15.7)
960 FORMAT(G12.4, 6E15.7)

  STOP
END PROGRAM RIEMANN
!--------------------------------------------------------------------

!--------------------------------------------------------------------
SUBROUTINE TRI
!--------------------------------------------------------------------
  use precision
  use ppp
  implicit none 
  INTEGER I
  REAL (pr) :: PE(4)

!   TRIE LES PHOTONS SORTANTS PAR ENERGIE
  IF (POUT(4, 2).GT.POUT(4, 1)) THEN
     DO I=1, 4
        PE(I)=POUT(I, 1)
        POUT(I, 1)=POUT(I, 2)
        POUT(I, 2)=PE(I)
     ENDDO
  ENDIF
  IF (POUT(4, 3).GT.POUT(4, 1)) THEN
     DO I=1, 4
        PE(I)=POUT(I, 1)
        POUT(I, 1)=POUT(I, 3)
        POUT(I, 3)=PE(I)
     ENDDO
  ENDIF
  IF (POUT(4, 3).GT.POUT(4, 2)) THEN
     DO I=1, 4
        PE(I)=POUT(I, 2)
        POUT(I, 2)=POUT(I, 3)
        POUT(I, 3)=PE(I)
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE TRI
!--------------------------------------------------------------------
