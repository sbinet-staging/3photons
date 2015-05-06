#include "intel_compatibility.h"
#include "precision.h"
#include "ppp.h"
#include "hasa.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>

//mc.cxx:143: erreur: no matching function for call to 'ppp::RAMBO(const int&, double&, double [100], double&)'
//ppp.h:28: note: candidats sont: void ppp::RAMBO(int, double, double*, double**, double&)


//!--------------------------------------------------------------------------
//void ppp::RAMBO (const int N, double ET, double* XM, double** P, double& WT) {
void ppp::RAMBO (const int N, double ET, double* XM, double& WT) {
//!--------------------------------------------------------------------------
//!
//!                       RAMBO
//!
//!    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
//!
//! ****** VERSION LIGHT & LIGHTSPEED ******
//!
//!    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
//!    AUTHORS@D  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
//!    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
//!
//!    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
//!    ET = TOTAL CENTRE-OF-MASS ENERGY
//!    XM = PARTICLE MASSES ( DIM=100 )
//!    P  = PARTICLE MOMENTA ( DIM=(4, 100) )
//!    WT = WEIGHT OF THE EVENT
//!
//!--------------------------------------------------------------------------
  int ITMAX = 6;
  const int NP = 100;
  static int IBEGIN = 0, NM;
  static double TWOPI, PO2LOG, Z [NP], XMT;
  double C, S, F, RMAS, G, A, X, BQ;
  double ACC = 1.0e-14;
  //double XM [NP];
  //double P [4] [NP];
  double Q [4] [NP], R [4], B [3];
  int IWARN [5] = {0, 0, 0, 0, 0}; // DATA IWARN/5*0/;
  //double RN;

//! INITIALIZATION STEP@D FACTORIALS FOR THE PHASE SPACE WEIGHT
  if (IBEGIN == 0) {
    std::cout << "IBegin" << std::endl;
    IBEGIN = 1;
    TWOPI  = 8.0 * atan (1.0);
//!        NORM=(TWOPI**(4-3*INP))/DBLE(ITOT)
    PO2LOG = log (TWOPI / 4.0);
    Z [1] = PO2LOG;
    for (int K = 2; K < NP; K++) {
      //Z [K] = Z [K-1] + PO2LOG - 2.0 * log ((double) (K-2));
      Z [K] = Z [K-1] + PO2LOG - 2.0 * log ((double) (K-1));
     } // ENDDO
    for (int K = 2; K < NP; K++) {
      //Z [K] = (Z [K] - log ((double) (K-1)));
      Z [K] = (Z [K] - log ((double) K));
    } // ENDDO

//! CHECK ON THE NUMBER OF PARTICLES
    if ((N > 1) && (N < 101)) goto l104;
    //PRINT 1001, N;
    //stop;
    exit (0);

//! CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
    l104:
    XMT=0.0;
    NM=0;
    for (int I = 0; I < N; I++) {
      if (XM [I] != 0.0) NM++;
      XMT = XMT + fabs (XM [I]);
    } // ENDDO
    if (XMT <= ET) goto l201;
    //PRINT 1002, XMT, ET;
    //STOP;
    exit (0);
  } // ENDIF

//! THE PARAMETER VALUES ARE NOW ACCEPTED

//! GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
  l201:
  for (int I=0; I < N; I++) {
    //C = 2.0 * RN [1] - 1.0;
    //C = 2.0 * ((double) random () / (double) RAND_MAX) - 1.0;
    C = 2.0 * RN () - 1.0;
    S = sqrt (1.0 - C*C);
    F = TWOPI * RN ();
    //F = TWOPI * RN [2];
    //F = TWOPI * ((double) random () / (double) RAND_MAX);
    Q [3] [I] = -log (RN () * RN ());
    //Q [4] [I] = -log (RN [3] * RN [4]);
    //Q [3] [I] = -log (((double) random () / (double) RAND_MAX) * ((double) random () / (double) RAND_MAX));
    Q [2] [I] = Q [3] [I] * C;
    Q [1] [I] = Q [3] [I] * S * cos (F);
    Q [0] [I] = Q [3] [I] * S * sin (F);
    //!write (*, *) 'C, F, Q (4, I): ', C, F, Q [4] [I]
    //std::cout << "TST: " << C << " " << F << " " << Q [3] [I] << " " << std::endl;
  } // ENDDO

//! CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
  for (int I = 0; I < 4; I++) {
    R [I] = 0.0;
  } // ENDDO
  for (int I = 0; I < N; I++) {
    for (int K = 0; K < 4; K++) {
      R [K] = R [K] + Q [K] [I];
    } // ENDDO
  } // ENDDO
  RMAS = sqrt (pow (R [3], 2) - pow (R [2], 2) - pow (R [1], 2) - pow (R [0], 2));
  for (int K = 0; K < 3; K++) {
    B [K] = -R [K] / RMAS;
  } // ENDDO
  G = R [3] / RMAS;
  A = 1.0 /(1.0 +G);
  X = ET/RMAS;

//! TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
  for (int I=0; I < N; I++) {
    BQ = B [0] * Q [0] [I] + B [1] * Q [1] [I] + B [2] * Q [2] [I];
    for (int K=0; K < 3; K++) {
      POUT [K] [I] = X * (Q [K] [I] + B [K] * (Q [3] [I] + A * BQ));
    } // ENDDO
    POUT [3] [I] = X * (G * Q [3] [I] + BQ);
  } // ENDDO

//! CALCULE LE POIDS ET LES AVERTISSEMENTS EVENTUELS
  WT = PO2LOG;
  //if (N != 2) WT = (2.0 * N - 4.0) * log (ET) + Z [N];
  if (N != 2) WT = (2.0 * N - 4.0) * log (ET) + Z [N-1];
  if (WT >= -180.0 ) goto l208;
  //if (IWARN [1] <= 5) PRINT 1004, WT;
  IWARN [1]++;
    l208:
  if (WT <= 174.0) goto l209;
  //if (IWARN [2] <= 5) PRINT 1005, WT;
  IWARN [2]++;

//! RENVOIE LES IMPULSIONS SANS MASSES PONDEREES
  l209:
  WT = exp (WT);
  return;

//1001 FORMAT(" RAMBO FAILS@D # OF PARTICLES =", I5, " IS NOT ALLOWED")
//1002 FORMAT(" RAMBO FAILS@D TOTAL MASS =", D15.6, " IS NOT", &
//            " SMALLER THAN TOTAL ENERGY =", D15.6)
//1004 FORMAT(" RAMBO WARNS@D WEIGHT = EXP(", F20.9, ") MAY UNDERFLOW")
//1005 FORMAT(" RAMBO WARNS@D WEIGHT = EXP(", F20.9, ") MAY  OVERFLOW")

  } // END SUBROUTINE RAMBO
