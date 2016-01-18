!>    @mainpage About the 3 photons simple Monte-Carlo
!>      @section intro1 Introduction (for the physicist)
!>       This small computational program computes cross-section for the particle physics process
!> electron + positron gives three photons.
!> It distinguishes a classical Standard Model contribution, of purely Quantum ElectroDynamic origin
!> and an hypothetic, beyond the Standard Model, New Physics contribution, phenomenologically
!> described by two effectives operators.
!> It was designed in the LEP era, so these new interactions occurs between the Z⁰ boson and the three photons.
!>  The effective operator can be related to specific models, among which magnetic monopoles that run in a four points loop.
!>  The two operators exhibit different
!>      @section intro2 Introduction (for the numerical guy)
!> The physicist want to compute a (multidimensional) integral, so we chose a Monte Carlo algorithm
!>      @section intro3 Introduction (for the computer guy)
!> this program started in a purely procedural style:
!> read in parameters and initialise counters
!> loop over (random) event,
!>           determining their geometrical and energy configuration,
!>           their phase space weight,
!>           their transition probability for each polarisation/helicity configuration,
!>                   depending on coupling strength
!>           sum it up
!>  then display / store the result.
!>  The use of common (for the original Fortran) or struct (in C) or record types (in Ada) or classes (in C++)
!> illustrates an object oriented design.
!>  The fact that we can plug each phase's output as the input of the next phase lend to a functionnal approach.
!>
!> @brief A Monte-Carlo generator for three-photons production at LEP with New Physics Anomalous Couplings
!> reads in parameters (physical as well as computational), and run the main Monte-Carlo loop
!> @author Vincent C. LAFAGE
PROGRAM RIEMANN
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
  real (pr) :: ETOT      !< total center-of-mass-energy
  real (pr) :: FLUX      !< flux factor
  real (pr) :: NORM      !< phase space normalisation factor
  real (pr) :: POIDS     !< event by event contribution
  real (pr) :: WTEV      !< phase space weight of the event
  real (pr) :: SIGMA     !< total cross-section
  real (pr) :: PREC      !< relative accuracy of cross-section estimator
  real (pr) :: VARIANCE  !< variance of total cross-section estimator
  real (pr) :: POIDS2    !< event by event pure QED contribution
  real (pr) :: PROBA     !< probability of given event
  real (pr), dimension (100) :: MFIN !< array of masses for final particles
  integer ITOT           !< non-cut number of samples
  integer NTOT           !< total number of samples
  integer NBIN           !< number of histogramming bin
  integer I              !< likely dummy loop variables ?
  integer K              !< likely dummy loop variables ?
  integer Sp             !< dummy loop variables on spin and helicities
  integer L1             !< dummy loop variables on spin and helicities
  integer L2             !< dummy loop variables on spin and helicities
  integer L3             !< dummy loop variables on spin and helicities

  real (pr) :: ALPHA     !< QED coupling constant squared over \f$4\pi\f$
  real (pr) :: CONVERS   !< \f$\mathrm{GeV}^{-2}\to\mathrm{pb}\f$ conversion factor
  real (pr) :: SIN2W     !< \f$\sin^2\theta_W\f$ squared sine of electroweak mixing angle
  real (pr) :: COS2W     !< \f$\cos^2\theta_W\f$ squared cosine of electroweak mixing angle
  real (pr) :: FACTCOM   !< conversion factor taking permutation into account
  real (pr) :: E2        !< QED coupling constant squared
  real (pr) :: ALPHAZ    !< QED coupling constant squared over \f$4\pi\f$, at \f$Z^0\f$ scale
  real (pr) :: E2Z       !< QED coupling constant squared, at \f$Z^0\f$ scale
  real (pr) :: GZR       !< adimensional relative \f$Z^0\f$ width
  real (pr) :: BETAPLUS  !< dual lagrangien \f$\beta_+\f$ adimensional coupling constant (input parameter)
  real (pr) :: BETAMOINS !< dual lagrangien \f$\beta_-\f$ adimensional coupling constant (input parameter)
  real (pr) :: BETAMIN   !< ...
  real (pr) :: SSP       !< ...
  real (pr) :: SSM       !< ...
  real (pr) :: INCSSP    !< ...
  real (pr) :: INCSSM    !< ...
  real (pr) :: PAA       !< electroweak prefactors of various component of amplitude squared
  real (pr) :: PAB       !< electroweak prefactors of various component of amplitude squared
  real (pr) :: PBB       !< electroweak prefactors of various component of amplitude squared
  real (pr) :: CAA       !< prefactors of various component of amplitude
  real (pr) :: CAB       !< prefactors of various component of amplitude
  real (pr) :: CBB       !< prefactors of various component of amplitude
  real (pr) :: DZETA     !< adimensional coefficient for \f$Z^0\f$ propagator modeling
  real (pr) :: PROPAG    !< adimensional coefficient for \f$Z^0\f$ propagator modeling
  real (pr) :: ECARTPIC  !< adimensional coefficient for \f$Z^0\f$ propagator modeling
  real (pr) :: BREPEM    !< \f$Z^0\to e^+e^-\f$ Branching Ratio
  logical PLOT           !< selector flag for histogramming
  integer ICYCLE    !< ...
  include 'paw.f90'
  character*80 A    !< ...
  integer, parameter :: UNL=10       !< ...
  integer, parameter :: UNE=20       !< ...
  LOGICAL CUT

  REAL*4 etime, tarray(2)
  CHARACTER*8 :: TIMESTR
  CHARACTER*9 :: TODAY = ''
  real (pr), parameter :: PI = 4 * atan (1._pr)        !< 𝜋, Archimede's constant

!$    integer nb_taches
!$    integer rang
!$    integer iii
! $    integer OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!   LECTURE DES PARAMETRES
!< Reading input parameters
  open (UNIT=UNL, FILE='valeurs', STATUS='OLD')
  READ (UNL, *) ITOT, A
  READ (UNL, *) ETOT, A
  READ (UNL, *) ACUT, A
  READ (UNL, *) BCUT, A
  READ (UNL, *) EMIN, A
  READ (UNL, *) SINCUT, A
  READ (UNL, *) ALPHA, A
  READ (UNL, *) ALPHAZ, A
  READ (UNL, *) CONVERS, A
  READ (UNL, *) MZ0, A
  READ (UNL, *) GZ0, A
  READ (UNL, *) SIN2W, A
  READ (UNL, *) BREPEM, A
  READ (UNL, *) BETAPLUS, A
  READ (UNL, *) BETAMOINS, A
  READ (UNL, *) NBIN, A
  READ (UNL, *) IMPR, A
  READ (UNL, *) PLOT, A
  CLOSE(UNIT=UNL)

!   PAW initialisation
!  CALL INIBOOK(ETOT, NBIN)
!< Sets final state particles masses
  MFIN (1)=0._pr
  MFIN (2)=0._pr
  MFIN (3)=0._pr
!<   flux factor (=1/2s for 2 initial massless particles)
  FLUX=1._pr/(2*ETOT**2)
!<   Includes total phase space normalisation
  NORM=((2*PI)**(4-3*INP))/real (ITOT, pr)
!<  Common factor, non-averaged over spins=1
!<                  /(symmetry factor)
!<                  *(conversion factor GeV^-2->pb)
!<  To average over spins, add :
!<                  /(number of incoming helicities)
  FACTCOM = 1._pr / 6._pr *CONVERS
  E2      = 4._pr * PI *ALPHA
  E2Z     = 4._pr * PI *ALPHAZ
  COS2W   = 1._pr - SIN2W
  GZR = GZ0 / MZ0
!>   RAC8 factor arise from SQRT (2) normalisation in each polarisation vector
!>   in helicity amplitudes method
  RAC8=SQRT(8._pr) ! dorénavant, c'est un parametre... nan, pas vraiment
!>   Couplings
  GA  = -SQRT (E2)**3
  GBP = -SQRT (E2Z/(4*COS2W*SIN2W))/MZ0**4
  GBM = GBP
!>   sum over polarisations factors
  POLP  =    -2 * SIN2W
  POLM  = 1 - 2 * SIN2W
  POLP2 = POLP**2
  POLM2 = POLM**2
  PAA=2
  PAB=1-4*SIN2W
  PBB=1-4*SIN2W+8*SIN2W**2
!>   Homogeneity coefficient
  CAA = FACTCOM * PAA
  CAB = FACTCOM * PAB / MZ0**2
  CBB = FACTCOM * PBB / MZ0**4
!>   Weight of a 3 photons phase space event
  WTEV = PI**2/8._pr*ETOT**2
!   PASSAGE EN VARIABLE A-DIMENSIONEE...
  DZETA    = (ETOT/MZ0)**2
  ECARTPIC = (DZETA-1._pr)/GZR
  PROPAG   = 1._pr/(ECARTPIC**2+1._pr)
!> Incoming momenta initialisation
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
!>   Start of integration loop
  DO I=1, ITOT

     CALL RAMBO(INP, ETOT, MFIN, POUT, WTEV)
     !write (*, '(Z16)') POUT
     WTEV=WTEV*NORM
!> Sort outgoing photons by energy
     CALL TRI
!>    Spinor inner product, scalar product and center-of-mass frame angles computation
     CALL PRECALC
!>  Compute event total weight including Matrix element
!>  (with cuts if needed)
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
!>    Store event parameters in an histogram
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
!> end of sampling loop

!> Computing the relative uncertainties
  DO K=1, NRESUL
     VAR(1, K)=(VAR(1, K)-SPM2(1, K)**2/REAL (ITOT, pr))/REAL (ITOT-1, pr)
     VAR(1, K)=SQRT(VAR(1, K)/REAL (ITOT, pr))/ABS(SPM2(1, K) &
          /REAL (ITOT, pr))
  ENDDO
!> Copy for opposite spins
  DO K=1, NRESUL
     SPM2(2, K)=SPM2(1, K)
     VAR(2, K)=VAR(1, K)
  ENDDO
!> Polarisations
  DO K=2, 3
     SPM2(1, K)=SPM2(1, K)*POLM2
     SPM2(2, K)=SPM2(2, K)*POLP2
  ENDDO
  DO K=4, 5
     SPM2(1, K)=SPM2(1, K)*POLM
     SPM2(2, K)=SPM2(2, K)*POLP
  ENDDO
!> Physical coefficients and Z⁰ propagator
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

  VARIANCE=(VARIANCE-SIGMA**2/REAL (ITOT, pr))/REAL (ITOT-1, pr)
  PREC=SQRT(VARIANCE/REAL (ITOT, pr))/ABS(SIGMA/REAL (ITOT, pr))
  SIGMA=SIGMA*FLUX
!>   Histograms normalisation
!  CALL NORMA(FLUX*NBIN)
!  CALL NORMSUP(ETOT, 4./(SPM2(1, 1)+SPM2(2, 1)), &
!                    4./(SPM2(1, 2)+SPM2(2, 2)), &
!                    4./(SPM2(1, 3)+SPM2(2, 3)))
!   PAW results storage
!  CALL HROUT(0, ICYCLE, ' ')
!  CALL HREND('DON')
!>   Numerical results storage
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
       tarray(1)/REAL (ITOT, pr)
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

!> @brief A small permutation routine that sorts the outgoing photons according to their energy
!--------------------------------------------------------------------
SUBROUTINE TRI
!--------------------------------------------------------------------
  use precision
  use ppp
  implicit none
  INTEGER I
  REAL (pr) :: PE(4)

!   TRIE LES PHOTONS SORTANTS PAR ENERGIE
!> sort outgoing photons by energy
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
