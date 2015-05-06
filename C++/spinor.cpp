#include "intel_compatibility.h"
#include "precision.h"
#include "spinor.h"
#include <complex>

std::complex<double> spinor::APPM (int k1, int k2, int k3) {
  return (-RAC8*S [0] [1] * pow (S [0] [k3], 2)/S [0] [k1]/S [0] [k2]/S [1] [k1]/S [1] [k2]);
}

std::complex<double> spinor::APMM (int k1, int k2, int k3) {
  return (-RAC8*T [0] [1] * pow (T [1] [k1], 2)/T [1] [k2]/T [1] [k3]/T [0] [k2]/T [0] [k3]);
}

std::complex<double> spinor::BPPM (int k1, int k2, int k3) {
  return (-RAC8*T [0] [1] * pow (T [k1] [k2] * S [k3] [0], 2));
}

std::complex<double> spinor::BPMM (int k1, int k2, int k3) {
  return (-RAC8*S [0] [1] * pow (T [k1] [1] * S [k2] [k3], 2));
}

std::complex<double> spinor::BPPP (int k1, int k2, int k3) {
  return (-RAC8*S [0] [1] * (
       + pow (T [k1] [k2] * T [k3] [1], 2)
       + pow (T [k1] [k3] * T [k2] [1], 2)
       + pow (T [k2] [k3] * T [k1] [1], 2)
       ));
}

std::complex<double> spinor::BMMM (int k1, int k2, int k3) {
  return (-RAC8*T [0] [1] * (
       + pow (S [k1] [0] * S [k2] [k3], 2)
       + pow (S [k2] [0] * S [k1] [k3], 2)
       + pow (S [k3] [0] * S [k1] [k2], 2)
       ));
}

spinor::spinor (ppp oPpp) {
  double P [5] [4];
  double TX, XX [5];
  std::complex<double> FX [5], CX;

  for (int J = 0; J < 4; J++){
    P [0] [J] = oPpp.P1 [J];
    P [1] [J] = oPpp.P2 [J];
    P [2] [J] = oPpp.POUT [J] [0];
    P [3] [J] = oPpp.POUT [J] [1];
    P [4] [J] = oPpp.POUT [J] [2];
  }
  for (int K = 0; K < 5; K++){
    S [K] [K] = 0.0;
    TX = sqrt (P [K] [3]+P [K] [2]);
    XX [K] = TX;
    FX [K] = std::complex<double> (P [K] [0], P [K] [1]) / TX;
    //std::cout << "K " << K << " FX " << FX [K] << std::endl;
  }
  //std::cout << std::endl;

  for (int J = 0; J < 4; J++){
    for (int K = J+1; K < 5; K++){
//! PRODUIT SPINORIEL DE M.MANGANO,  S.PARKE
      CX = FX [J] * XX [K] - FX [K] * XX [J];
      S [J] [K] =  CX;
      S [K] [J] = -CX;
      T [K] [J] =  conj (CX);
      T [J] [K] = -conj (CX);
      //std::cout << "J " << J << " K " << K << " CX " << CX << " bar " << conj (CX) << " abs " << abs (CX) << std::endl;
    }
  }
  //std::cout << std::endl;
}
