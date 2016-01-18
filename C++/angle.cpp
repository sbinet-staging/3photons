#include "intel_compatibility.h"
#include "precision.h"
#include "angle.h"

#include "ppp.h"
#include <math.h>

angle::angle (ppp oPpp, scalar oScalar) {
  double P [5] [4];

  for (int J = 0; J < 4; J++){
    P [0] [J] = oPpp.P1 [J];
    P [1] [J] = oPpp.P2 [J];
    P [2] [J] = oPpp.POUT [J] [0];
    P [3] [J] = oPpp.POUT [J] [1];
    P [4] [J] = oPpp.POUT [J] [2];
  }

  COS1P1 = 1.0 - oScalar.PS [0] [2] / (P [0] [3] * P [2] [3]);
  COS2P1 = 1.0 - oScalar.PS [0] [3] / (P [0] [3] * P [3] [3]);
  COS3P1 = 1.0 - oScalar.PS [0] [4] / (P [0] [3] * P [4] [3]);

  COS12  = 1.0 - oScalar.PS [2] [3]/ (oPpp.POUT [3] [0] * oPpp.POUT [3] [1]);
  COS13  = 1.0 - oScalar.PS [2] [4]/ (oPpp.POUT [3] [0] * oPpp.POUT [3] [2]);
  COS23  = 1.0 - oScalar.PS [3] [4]/ (oPpp.POUT [3] [1] * oPpp.POUT [3] [2]);

  //std::cout << "COS1P1 " << COS1P1 << std::endl ;
  //std::cout << "COS2P1 " << COS2P1 << std::endl ;
  //std::cout << "COS3P1 " << COS3P1 << std::endl ;
  //std::cout << "COS12  " << COS12  << std::endl ;
  //std::cout << "COS13  " << COS13  << std::endl ;
  //std::cout << "COS23  " << COS23  << std::endl ;

  double NX, NY, NZ, NN;
  NX = oPpp.POUT [1] [0] * oPpp.POUT [2] [1] - oPpp.POUT [1] [1] * oPpp.POUT [2] [0];
  NY = oPpp.POUT [2] [0] * oPpp.POUT [0] [1] - oPpp.POUT [2] [1] * oPpp.POUT [0] [0];
  NZ = oPpp.POUT [0] [0] * oPpp.POUT [1] [1] - oPpp.POUT [0] [1] * oPpp.POUT [1] [0];
  NN = sqrt (pow (NX, 2) + pow (NY, 2) + pow (NZ, 2));

  COSN  =  (NX*P [0] [0]+NY*P [0] [1]+NZ*P [0] [2])/P [0] [3]/NN;

  COSAC =  (P [2] [1] * P [3] [1] +P [2] [2]*P [3] [2])/sqrt (
           (P [2] [1] * P [2] [1] +P [2] [2]*P [2] [2])* 
           (P [3] [1] * P [3] [1] +P [3] [2]*P [3] [2]));
  //std::cout << "COSN   " << COSN   << std::endl ;
  //std::cout << "COSAC  " << COSAC  << std::endl ;
  //std::cout << "       " <<           std::endl ;
}

//! Determine whether an event passes cut or not
bool angle::CUT (ppp oPpp, cutpar oCutpar) {
  bool xCut;

  xCut = false;
  for (int I=0; I < INP; I++) {
    xCut = xCut || (oPpp.POUT [3] [I] < oCutpar.EMIN);
  }
  xCut = xCut ||
    (fabs (COS1P1) > oCutpar.ACUT) ||
    (fabs (COS2P1) > oCutpar.ACUT) ||
    (fabs (COS3P1) > oCutpar.ACUT) ||
    (COS12 > oCutpar.BCUT) ||
    (COS13 > oCutpar.BCUT) ||
    (COS23 > oCutpar.BCUT) ||
    (fabs (COSN) < oCutpar.SINCUT);

  return (xCut);
}
