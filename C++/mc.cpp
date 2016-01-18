#include "intel_compatibility.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <iosfwd>
#include "precision.h"
#include "ppp.h"
#include "spinor.h"
#include "scalar.h"
#include "angle.h"
#include "result.h"
#include "cutpar.h"
#include "resfin.h"
#include "param.h"
#include "xtrpro.h"
extern "C" {
#include "unixtime.h"
}
#include "time.h"

// g++ mc.cxx angle.cxx ppp.cxx spinor.cxx result.cxx resfin.cxx fawzi.cxx

//!    \mainpage About the 3 photons simple Monte-Carlo
//!      \section intro1 Introduction (for the physicist)
//!       This small computational program computes cross-section for the particle physics process
//! electron + positron gives three photons (e⁺e⁻ → 𝛾𝛾𝛾).
//! It distinguishes a classical Standard Model contribution, of purely Quantum ElectroDynamic origin
//! and an hypothetic, beyond the Standard Model, New Physics contribution, phenomenologically
//! described by two effectives operators.
//! It was designed in the LEP era, so these new interactions occurs between the Z⁰ boson and the three photons.
//!  The effective operator can be related to specific models, among which magnetic monopoles that run in a four points loop.
//!  The two operators exhibit different
//!      \section intro2 Introduction (for the numerical guy)
//! The physicist want to compute a (multidimensional) integral, so we chose a Monte Carlo algorithm
//!      \section intro3 Introduction (for the computer guy)
//! this program started in a purely procedural style:
//! read in parameters and initialise counters
//! loop over (random) event,
//!           determining their geometrical and energy configuration,
//!           their phase space weight,
//!           their transition probability for each polarisation/helicity configuration,
//!                   depending on coupling strength
//!           sum it up
//!  then display / store the result.
//!  The use of common (for the original Fortran) or struct (in C) or record types (in Ada) or classes (in C++)
//! illustrates an object oriented design.
//!  The fact that we can plug each phase's output as the input of the next phase lend to a functionnal approach.

//! \file mc.cpp function for the 3 photons simple MonteCarlo
//! \brief A Monte-Carlo generator for three-photons production at LEP with New Physics Anomalous Couplings
//! reads in parameters (physical as well as computational), and run the main Monte-Carlo loop
//! \author Vincent C. LAFAGE
//! \date 2007-08-10 ISO
int main ()
{
  double ETOT, FLUX, PI, NORM, POIDS, WTEV, SIGMA, PREC, VARIANCE;
  double POIDS2, PROBA;
  double MFIN [100];
  int ITOT, NTOT, NBIN;
  int I;
  double ALPHA, CONVERS, SIN2W, COS2W, FACTCOM, E2, ALPHAZ, E2Z, GZR;
  double BETAPLUS, BETAMOINS, BETAMIN;
  double SSP, SSM, INCSSP, INCSSM;
  double PAA, PAB, PBB;
  double CAA, CAB, CBB;
  double DZETA, PROPAG, ECARTPIC;
  double BREPEM;
  bool CUT, PLOT;
  int ICYCLE;
  //INCLUDE 'paw.f90'
  std::string A; // char* A;
  //int UNL=10, UNE=20;
  float etime, tarray [2];
  //std::string TIMESTR = ""; // char* TIMESTR;
  //std::string TODAY   = ""; // char* TODAY = '';
  char TODAY_s[10];
  char TIMESTR_s[9];

  cutpar oCutpar;
  param  oParam;
  //spinor oSpinor;
  //ppp    oPpp;
  //result oResult;
  resfin oResfin;
  xtrpro oXtrpro;

  int Saved_Time = unixtime ();
//   LECTURE DES PARAMETRES
//!< Reading input parameters
  std::ifstream fluxFic ("valeurs", std::ios::in);

  fluxFic >> ITOT;           getline (fluxFic, A);
  if (fluxFic.fail ()) {fluxFic.clear (); std::cout << "FAIL\n";}
  std::cout << ITOT << " | " << A << std::endl;
  fluxFic >> ETOT;           getline (fluxFic, A);
  if (fluxFic.fail ()) {fluxFic.clear (); std::cout << "FAIL\n";}
  std::cout << ETOT << " | " << A << std::endl;
  fluxFic >> oCutpar.ACUT;   getline (fluxFic, A);
  if (fluxFic.fail ()) {fluxFic.clear (); std::cout << "FAIL\n";}
  std::cout << oCutpar.ACUT << " | " << A << std::endl;
  fluxFic >> oCutpar.BCUT;   getline (fluxFic, A);
  if (fluxFic.fail ()) {fluxFic.clear (); std::cout << "FAIL\n";}
  std::cout << oCutpar.BCUT << " | " << A << std::endl;
  fluxFic >> oCutpar.EMIN;   getline (fluxFic, A);
  fluxFic >> oCutpar.SINCUT; getline (fluxFic, A);
  fluxFic >> ALPHA;          getline (fluxFic, A);
  fluxFic >> ALPHAZ;         getline (fluxFic, A);
  fluxFic >> CONVERS;        getline (fluxFic, A);
  fluxFic >> oParam.MZ0;     getline (fluxFic, A);
  fluxFic >> oParam.GZ0;     getline (fluxFic, A);
  fluxFic >> SIN2W;          getline (fluxFic, A);
  fluxFic >> BREPEM;         getline (fluxFic, A);
  fluxFic >> BETAPLUS;       getline (fluxFic, A);
  fluxFic >> BETAMOINS;      getline (fluxFic, A);
  fluxFic >> NBIN;           getline (fluxFic, A);
  if (fluxFic.fail ()) {fluxFic.clear (); std::cout << "FAIL: NBIN\n";}
  fluxFic >> oParam.IMPR;    getline (fluxFic, A);
  if (fluxFic.fail ()) {fluxFic.clear (); std::cout << "FAIL: oParam.IMPR\n";}
  fluxFic >> PLOT;           getline (fluxFic, A);
  if (fluxFic.fail ()) {fluxFic.clear (); std::cout << "FAIL: PLOT\n";}

  fluxFic.close ();

  std::cout << "ITOT           : " << ITOT           << std::endl;
  std::cout << "ETOT           : " << ETOT           << std::endl;
  std::cout << "oCutpar.ACUT   : " << oCutpar.ACUT   << std::endl;
  std::cout << "oCutpar.BCUT   : " << oCutpar.BCUT   << std::endl;
  std::cout << "oCutpar.EMIN   : " << oCutpar.EMIN   << std::endl;
  std::cout << "oCutpar.SINCUT : " << oCutpar.SINCUT << std::endl;
  std::cout << "ALPHA          : " << ALPHA          << std::endl;
  std::cout << "ALPHAZ         : " << ALPHAZ         << std::endl;
  std::cout << "CONVERS        : " << CONVERS        << std::endl;
  std::cout << "oParam.MZ0     : " << oParam.MZ0     << std::endl;
  std::cout << "oParam.GZ0     : " << oParam.GZ0     << std::endl;
  std::cout << "SIN2W          : " << SIN2W          << std::endl;
  std::cout << "BREPEM         : " << BREPEM         << std::endl;
  std::cout << "BETAPLUS       : " << BETAPLUS       << std::endl;
  std::cout << "BETAMOINS      : " << BETAMOINS      << std::endl;
  std::cout << "NBIN           : " << NBIN           << std::endl;
  std::cout << "oParam.IMPR    : " << oParam.IMPR    << std::endl;
  std::cout << "PLOT           : " << PLOT           << std::endl;

  oParam.IMPR    = false;
  PLOT           = false;

  std::cout << "2oParam.IMPR   : " << oParam.IMPR    << std::endl;
  std::cout << "2PLOT          : " << PLOT           << std::endl;
  std::cout << "1/alpha  " << 1.0/ALPHA  << " alpha  " << ALPHA  << std::endl ;
  std::cout << "1/alphaz " << 1.0/ALPHAZ << " alphaz " << ALPHAZ << std::endl ;

  ppp    oPpp (ETOT);

//  PAW initialisation
//  CALL INIBOOK(ETOT, NBIN)
//!< Sets final state particles masses
  MFIN [0] = 0.0;
  MFIN [1] = 0.0;
  MFIN [2] = 0.0;
//!<   flux factor (=1/2s for 2 initial massless particles)
  FLUX = 1.0 / (2.0*ETOT*ETOT);
//!<   Includes total phase space normalisation
  PI = 4.0 * atan (1.0);
  NORM = pow (2.0*PI, 4-3*INP) / ((double) ITOT);
//!<  Common factor, non-averaged over spins=1
//!<                  /(symmetry factor)
//!<                  *(conversion factor GeV^-2->pb)
//!<  To average over spins, add :
//!<                  /(number of incoming helicities)
  FACTCOM = 1.0/6.0*CONVERS;
  E2      = 4.0*PI*ALPHA;
  E2Z     = 4.0*PI*ALPHAZ;
  COS2W   = 1.0-SIN2W;
  GZR     = oParam.GZ0 / oParam.MZ0;

//!>   RAC8 factor arise from SQRT (2) normalisation in each polarisation vector
//!>   in helicity amplitudes method
  //oSpinor.RAC8 = sqrt (8.0); //! dorénavant, c'est un paramètre... nan, pas vraiment
//!   Couplings
  oParam.GA  = -pow (sqrt (E2), 3);
  oParam.GBP = -sqrt (E2Z/(4*COS2W*SIN2W)) / pow (oParam.MZ0, 4);
  oParam.GBM = oParam.GBP;
//!   sum over polarisations factors
  oParam.POLP  =     - 2.0 * SIN2W;
  oParam.POLM  = 1.0 - 2.0 * SIN2W;
  oParam.POLP2 = pow (oParam.POLP, 2);
  oParam.POLM2 = pow (oParam.POLM, 2);
  PAA   = 2.0;
  PAB   = 1.0-4.0*SIN2W;
  PBB   = 1.0-4.0*SIN2W+8.0*pow (SIN2W, 2);
//!   Homogeneity coefficient
  CAA   = FACTCOM*PAA;
  CAB   = FACTCOM*PAB / pow (oParam.MZ0, 2);
  CBB   = FACTCOM*PBB / pow (oParam.MZ0, 4);
//!   Weight of a 3 photons phase space event
  WTEV  = pow (PI, 2)/8.0 * pow (ETOT, 2);
//!   Passage en variable a-dimensionée...
  DZETA = pow ((ETOT/oParam.MZ0), 2);
  ECARTPIC = (DZETA-1.0)/GZR;
  PROPAG = 1.0/(pow (ECARTPIC, 2)+1.0);
//!   Incoming momenta initialisation

  //Facteur 3.81971807742127067698

//!   Initialisation des cumulants
  NTOT = 0;
  for (int Sp=0; Sp < 2; Sp++) {
    for (int K=0; K < NRESUL; K++) {
      oResfin.SPM2 [Sp] [K] = 0.0;
      oResfin.VAR  [Sp] [K] = 0.0;
    }
  }
  SIGMA    = 0.0;
  VARIANCE = 0.0;

//!$OMP PARALLEL DO SHARED (SPM2DIF, SPM2, VAR, SIGMA, VARIANCE, NTOT)
//!     rang = OMP_GET_THREAD_NUM ()
//!>  Start of integration loop
  for (I=0; I < ITOT; I++) {
    //RAMBO (INP, ETOT, MFIN, oPpp.POUT, WTEV);
    oPpp.RAMBO (INP, ETOT, MFIN, WTEV);
    WTEV = WTEV * NORM;
//!  Sort outgoing photons by energy
    oPpp.TRI ();
//!   Spinor inner product, scalar product and center-of-mass frame angles computation
    spinor oSpinor (oPpp);
    scalar oScalar (oSpinor);
    angle  oAngle  (oPpp, oScalar);
//!  Compute event total weight including Matrix element
//!  (with cuts if needed)
    if (!oAngle.CUT (oPpp, oCutpar)) {
      result oResult (oParam, oSpinor, ETOT);
      oResult.display (oParam);
      for (int K=0; K < NRESUL; K++) {
        oResfin.SPM2DIF [K] = 0.0;
        for (int L1=0; L1 < 2; L1++) {
          for (int L2=0; L2 < 2; L2++) {
            for (int L3=0; L3 < 2; L3++) {
              oResfin.SPM2DIF [K] = oResfin.SPM2DIF [K] + oResult.M2 [L1] [L2] [L3] [K];
	      //std::cout << K << " " << L1 << " " << L2 << " " << L3 << " | " << oResult.M2 [L1] [L2] [L3] [K] << std::endl ;
            }
          }
        }
        oResfin.SPM2 [0] [K] = oResfin.SPM2 [0] [K] + oResfin.SPM2DIF [K];
        oResfin.VAR [0] [K]  = oResfin.VAR [0] [K]  + pow (oResfin.SPM2DIF [K], 2);
      }
      PROBA =
          CAA * oResfin.SPM2DIF [0]
        + CBB * (+pow (BETAPLUS,  2) * oResfin.SPM2DIF [1]
		 +pow (BETAMOINS, 2) * oResfin.SPM2DIF [2])
        / pow (GZR, 2) * PROPAG
        + CAB * 2 * BETAPLUS * (ECARTPIC * oResfin.SPM2DIF [3]
                                - oResfin.SPM2DIF [4])
        / GZR * PROPAG;
      POIDS = PROBA * WTEV / 4.0;
      SIGMA = SIGMA + POIDS;
      VARIANCE = VARIANCE + pow (POIDS, 2);

      //std::cout << "SPM2DIF [0] " << oResfin.SPM2DIF [0] << " SPM2DIF [1] " << oResfin.SPM2DIF [1] << " SPM2DIF [2] " << oResfin.SPM2DIF [2] << " SPM2DIF [3] " << oResfin.SPM2DIF [3] << " SPM2DIF [4] " << oResfin.SPM2DIF [4] << std::endl ;
      //std::cout << " " <<  << "  " <<  << "  " <<  << "  " <<  << std::endl ;
      //std::cout << "CAA   " << CAA   << " CBB  " << CBB  << " CAB   " << CAB   << " GZR   " << GZR << std::endl ;
      //std::cout << "PROBA " << PROBA << " WTEV " << WTEV << " POIDS " << POIDS << " SIGMA " << SIGMA << std::endl ;
//!  Store event parameters in an histogram
      if (PLOT) {
        POIDS2 = CAA * oResfin.SPM2DIF [0] * WTEV / 4.0;
        oXtrpro.PRPLUS  = ((float) CBB * pow (BETAPLUS,  2) * oResfin.SPM2DIF [1]
			   /pow (GZR, 2) * PROPAG
			   *WTEV/4.0);
        oXtrpro.PRMOINS = ((float) CBB * pow (BETAMOINS, 2) * oResfin.SPM2DIF [2]
			   /pow (GZR, 2) * PROPAG
			   *WTEV/4.0);
        //BOOK (POIDS, POIDS2);
      }
      NTOT = NTOT + 1;
    } else {
      POIDS = 0.0;
    };
  };
//!$OMP END PARALLEL DO
//! end of sampling loop

//! Computing the relative uncertainties
  for (int K=0; K < NRESUL; K++) {
    oResfin.VAR [0] [K] = (oResfin.VAR [0] [K]-pow (oResfin.SPM2 [0] [K], 2)/ ((double) ITOT)) / ((double) ITOT-1);
    oResfin.VAR [0] [K] = sqrt (oResfin.VAR [0] [K] / ((double) ITOT)) / fabs (oResfin.SPM2 [0] [K] / (double) ITOT);
  }
//! Copy for opposite spins
  for (int K=0; K < NRESUL; K++) {
    oResfin.SPM2 [1] [K] = oResfin.SPM2 [0] [K];
    oResfin.VAR  [1] [K] = oResfin.VAR  [0] [K];
  }
//! Polarisations
  for (int K=1; K <= 2;  K++) {
    oResfin.SPM2 [0] [K] = oResfin.SPM2 [0] [K] * oParam.POLM2;
    oResfin.SPM2 [1] [K] = oResfin.SPM2 [1] [K] * oParam.POLP2;
  }
  for (int K=3; K<= 4; K++) {
    oResfin.SPM2 [0] [K] = oResfin.SPM2 [0] [K] * oParam.POLM;
    oResfin.SPM2 [1] [K] = oResfin.SPM2 [1] [K] * oParam.POLP;
  }
//! Physical coefficients and Z⁰ propagator
  for (int Sp=0; Sp < 2; Sp++) {
     for (int K=0; K < NRESUL; K++) {
       oResfin.SPM2 [Sp] [K] = oResfin.SPM2 [Sp] [K] * FACTCOM * FLUX * WTEV;
     }
     oResfin.SPM2 [Sp] [0] = oResfin.SPM2 [Sp] [0];
     oResfin.SPM2 [Sp] [1] = oResfin.SPM2 [Sp] [1]/pow (GZR, 2)/pow (oParam.MZ0, 4)*PROPAG;
     oResfin.SPM2 [Sp] [2] = oResfin.SPM2 [Sp] [2]/pow (GZR, 2)/pow (oParam.MZ0, 4)*PROPAG;
     oResfin.SPM2 [Sp] [3] = oResfin.SPM2 [Sp] [3]/GZR/pow (oParam.MZ0, 2)*PROPAG*ECARTPIC;
     oResfin.SPM2 [Sp] [4] = oResfin.SPM2 [Sp] [4]/GZR/pow (oParam.MZ0, 2)*PROPAG;
  }
  BETAMIN = sqrt ((oResfin.SPM2 [0] [0]+oResfin.SPM2 [1] [0])/(oResfin.SPM2 [0] [1]+oResfin.SPM2 [1] [1]));
  SSP=(oResfin.SPM2 [0] [1]+oResfin.SPM2 [1] [1])/sqrt (oResfin.SPM2 [0] [0]+oResfin.SPM2 [1] [0])/2.0;
  SSM=(oResfin.SPM2 [0] [2]+oResfin.SPM2 [1] [2])/sqrt (oResfin.SPM2 [0] [0]+oResfin.SPM2 [1] [0])/2.0;
  INCSSP=
    +sqrt (+pow (oResfin.SPM2 [0] [1]*oResfin.VAR [0] [1], 2)
	   +pow (oResfin.SPM2 [1] [1]*oResfin.VAR [1] [1], 2))
    / fabs (oResfin.SPM2 [0] [1]+oResfin.SPM2 [1] [1])
    +sqrt (+pow (oResfin.SPM2 [0] [0]*oResfin.VAR [0] [0], 2)
           +pow (oResfin.SPM2 [1] [0]*oResfin.VAR [1] [0], 2))
    / fabs (oResfin.SPM2 [0] [0]+oResfin.SPM2 [1] [0])/2.0;
  INCSSM=
    +sqrt (+pow (oResfin.SPM2 [0] [2]*oResfin.VAR [0] [2], 2)
           +pow (oResfin.SPM2 [1] [2]*oResfin.VAR [1] [2], 2))
    / fabs (oResfin.SPM2 [0] [2]+oResfin.SPM2 [1] [2])
    +sqrt (+pow (oResfin.SPM2 [0] [0]*oResfin.VAR [0] [0], 2)
           +pow (oResfin.SPM2 [1] [0]*oResfin.VAR [1] [0], 2))
    / fabs (oResfin.SPM2 [0] [0]+oResfin.SPM2 [1] [0])/2.0;

  VARIANCE = (VARIANCE-pow (SIGMA, 2) / ((double) ITOT)) / ((double) ITOT-1);
  PREC = sqrt (VARIANCE / ((double) ITOT)) / fabs (SIGMA/((double) ITOT));
  SIGMA = SIGMA * FLUX;
//!  Histograms normalisation
  // PAW en C++, pas tout de suite//  NORMA (FLUX*NBIN); // CALL NORMA(FLUX*NBIN)
  // PAW en C++, pas tout de suite//  NORMSUP (ETOT, // CALL NORMSUP(ETOT,
  // PAW en C++, pas tout de suite//           4./(oResfin.SPM2 [1, 1)+oResfin.SPM2 [2, 1)),
  // PAW en C++, pas tout de suite//           4./(oResfin.SPM2 [1, 2)+oResfin.SPM2 [2, 2)),
  // PAW en C++, pas tout de suite//           4./(oResfin.SPM2 [1, 3)+oResfin.SPM2 [2, 3)));
//!   PAW results storage
  // PAW en C++, pas tout de suite//  HROUT (0, ICYCLE, " "); // CALL HROUT(0, ICYCLE, " ')
  // PAW en C++, pas tout de suite//  HREND("DON"); // CALL HREND('DON')
//!   Numerical results storage
  time_t t;
  time (&t);
  //size_t strftime (char *s, size_t max, const char *format, const struct tm *tm);
  size_t size_tmp1 = strftime (TODAY_s,  10, "%d-%b-%y", gmtime (&t));
  size_t size_tmp2 = strftime (TIMESTR_s, 9, "%T", gmtime (&t));
  std::string TIMESTR = std::string (TIMESTR_s); // char* TIMESTR;
  std::string TODAY   = std::string (TODAY_s);   // char* TODAY = '';

  std::ofstream UNE;
  UNE.open ("res.dat"); //, ios::in|ios::nocreate); // OPEN(UNIT=UNE, FILE="res.dat", STATUS="UNKNOWN");
  UNE << " " << TODAY << "   " << TIMESTR << std::endl; // WRITE(UNE, *) TODAY, "   ", TIMESTR
  UNE << std::endl;
  //UNE.setf (std::ios::showpoint);
  UNE << std::setw (16);
  UNE << std::setprecision (14);
  //UNE.setf (ios::fixed);
  UNE << " Nombre d'evenements            : " << ITOT           << std::endl;
  UNE << " ... apres coupure              : " << NTOT           << std::endl;
  UNE << " energie dans le CdM      (GeV) : " << ETOT           << std::endl;
  UNE << " coupure / cos(photon,faisceau) : " << oCutpar.ACUT   << std::endl;
  UNE << " coupure / cos(photon,photon)   : " << oCutpar.BCUT   << std::endl;
  UNE << " coupure / sin(normale,faisceau): " << oCutpar.SINCUT << std::endl;
  UNE << " coupure sur l'energie    (GeV) : " << oCutpar.EMIN   << std::endl;
  UNE << " 1/(constante de structure fine): " << 1.0/ALPHA       << std::endl;
  UNE << " 1/(structure fine au pic)      : " << 1.0/ALPHAZ      << std::endl;
  UNE << " facteur de conversion GeV-2/pb : " << CONVERS         << std::endl;
//UNE << " Volume d'espace des phases     : " << VOLUME/NTOT     << std::endl;
  UNE << " Masse du Z0              (GeV) : " << oParam.MZ0      << std::endl;
  UNE << " Largeur du Z0            (GeV) : " << oParam.GZ0      << std::endl;
  UNE << " Sinus^2 Theta Weinberg         : " << SIN2W           << std::endl;
  UNE << " Taux de branchement Z--->e+e-  : " << BREPEM          << std::endl;
  UNE << " Beta plus                      : " << BETAPLUS        << std::endl;
  UNE << " Beta moins                     : " << BETAMOINS       << std::endl;
  UNE << " ---------------------------------------------"        << std::endl;
  UNE << " Section Efficace          (pb) : " << SIGMA           << std::endl;
  UNE << " Ecart-Type                (pb) : " << SIGMA*PREC      << std::endl;
  UNE << " Precision Relative             : " << PREC            << std::endl;
  UNE << " ---------------------------------------------"        << std::endl;
  UNE << " Beta minimum                   : " << BETAMIN         << std::endl;
  UNE << " Stat. Significance  B+(pb-1/2) : " << SSP             << std::endl;
  UNE << " Incert. Stat. Sign. B+(pb-1/2) : " << SSP*INCSSP      << std::endl;
  UNE << " Stat. Significance  B-(pb-1/2) : " << SSM             << std::endl;
  UNE << " Incert. Stat. Sign. B+(pb-1/2) : " << SSM*INCSSM      << std::endl;
  UNE << " Temps ecoule                   : " << std::endl;//<< etime (tarray) << std::endl;
  //UNE << " Temps ecoule utilisateur       : " << (systtime () - Saved_Time)/sysconf_sc_clk_tck () << std::endl;//<< tarray (1)  << std::endl;
  UNE << " Temps ecoule utilisateur       : " << ((double) (unixtime () - Saved_Time))/ ((double) sysconf_sc_clk_tck ()) << std::endl;//<< tarray (1)  << std::endl;
  UNE << " Temps ecoule systeme           : " << std::endl;//<< tarray (2)  << std::endl;
  UNE << " Temps ecoule par evenement     : " << ((double) (unixtime () - Saved_Time)) / ((double) sysconf_sc_clk_tck ()) / ((double) ITOT) << std::endl;//<< tarray (1)/((double) ITOT) << std::endl;
  UNE.setf (std::ios::scientific);
  UNE << std::endl;
  UNE.precision (7);
  UNE.scientific;
  //  1  1  0.1412055E+01  0.3094754E-01  0.2191666E-01
  //  1  1 3.69675182e-018.10204612e-03 2.19166623e-02
  //  1  1  3.6967518e-018.1020461e-03  2.1916662e-02
  for (int Sp=0; Sp < 2; Sp++) {
     for (int K=0; K < NRESUL; K++) { // format 930
       //UNE << std::right << std::setw (3) << Sp+1 << std::setw (3) << K+1 << std::setw (15) << std::left << oResfin.SPM2 [Sp] [K] << std::setw (15) << std::left << fabs (oResfin.SPM2 [Sp] [K])*oResfin.VAR [Sp] [K] << std::setw (15) << std::left << oResfin.VAR [Sp] [K] << std::endl;
       UNE << std::right << std::setw (3) << Sp+1 << std::setw (3) << K+1 << std::setw (15) << oResfin.SPM2 [Sp] [K] << std::setw (15) << fabs (oResfin.SPM2 [Sp] [K])*oResfin.VAR [Sp] [K] << std::setw (15) << oResfin.VAR [Sp] [K] << std::endl;
     }
     UNE << std::endl;
  }
  for (int K=0; K < NRESUL; K++) { // format 940
    UNE << std::right << "   " << std::setw (3) << K+1 << std::setw (15) << (oResfin.SPM2 [0] [K]+oResfin.SPM2 [1] [K])/4.0 << std::setw (15) <<
    sqrt (+pow (oResfin.SPM2 [0] [K]*oResfin.VAR [0] [K], 2)
          +pow (oResfin.SPM2 [1] [K]*oResfin.VAR [1] [K], 2))/4.0 << std::setw (15) <<
    sqrt (+pow (oResfin.SPM2 [0] [K]*oResfin.VAR [0] [K], 2)
          +pow (oResfin.SPM2 [1] [K]*oResfin.VAR [1] [K], 2))/fabs (oResfin.SPM2 [0] [K]+oResfin.SPM2 [1] [K]) << std::endl;
  }
  oResfin.ERIC (BREPEM, CONVERS, PI, oParam);
  oResfin.FAWZI (BREPEM, CONVERS, PI, ETOT, oParam, oCutpar);
  UNE.close ();

  UNE.open ("pil.mc"); //, ios::in|ios::nocreate); // OPEN(UNIT=UNE, ACCESS="APPEND", FILE="pil.mc", STATUS="unknown")
  UNE << TODAY << "   " << TIMESTR << std::endl;

  UNE << ETOT <<
    (oResfin.SPM2 [0] [0]+oResfin.SPM2 [1] [0])/4.0 <<
    (oResfin.SPM2 [0] [1]+oResfin.SPM2 [1] [1])/4.0*pow (BETAPLUS,  2) <<
    (oResfin.SPM2 [0] [2]+oResfin.SPM2 [1] [2])/4.0*pow (BETAMOINS, 2) <<
    (oResfin.SPM2 [0] [3]+oResfin.SPM2 [1] [3])/4.0*BETAPLUS <<
   ((oResfin.SPM2 [0] [0]+oResfin.SPM2 [1] [0]) +
    (oResfin.SPM2 [0] [1]+oResfin.SPM2 [1] [1])*pow (BETAPLUS,  2) +
    (oResfin.SPM2 [0] [2]+oResfin.SPM2 [1] [2])*pow (BETAMOINS, 2) +
 2.*(oResfin.SPM2 [0] [3]+oResfin.SPM2 [1] [3])*BETAPLUS)/4.0 <<
    SIGMA;
  UNE.close ();

}
