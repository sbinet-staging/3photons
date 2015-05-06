#include "intel_compatibility.h"
//! Calcule l'élément de matrice e+e- en 3 photons
#include "precision.h"
#include "result.h"

#include "param.h"
#include "spinor.h"
#include <math.h>
#include <complex>
#include <iostream>

result::result (param  oParam, spinor oSpinor, double ETOT) {
  int l1, l2, l3;
  std::complex<double> A, BP, BM;

//! LES - VIENNENT DES Photons SORTANTS
  for (int LL1 = 0; LL1 < 2; LL1++) {
    l1 = - (2*LL1-1);
    for (int LL2 = 0; LL2 < 2; LL2++) {
      l2 = - (2*LL2-1);
      for (int LL3 = 0; LL3 < 2; LL3++) {
	l3 = - (2*LL3-1);
//! CALCULS DES AMPLITUDES D'HELICITE
	if ((l1 == +1) && (l2 == +1) && (l3 == +1)) {
	  A  = 0.0;
	  BP = 0.0;
	  BM = oSpinor.BPPP (2, 3, 4);
	  //std::cout << "+++ " << A  << " " << BP << " " << BM << " " << std::endl;
	}
	if ((l1 == -1) && (l2 == -1) && (l3 == -1)) {
	  A  = 0.0;
	  BP = 0.0;
	  BM = oSpinor.BMMM (2, 3, 4);
	  //std::cout << "--- " << A  << " " << BP << " " << BM << " " << std::endl;
	}
	if ((l1 == +1) && (l2 == +1) && (l3 == -1)) {
	  A  = oSpinor.APPM (2, 3, 4);
	  BP = oSpinor.BPPM (2, 3, 4);
	  BM = 0.0;
	  //std::cout << "++- " << A  << " " << BP << " " << BM << " " << std::endl;
	}
	if ((l1 == +1) && (l2 == -1) && (l3 == +1)) {
	  A  = oSpinor.APPM (4, 2, 3);
	  BP = oSpinor.BPPM (4, 2, 3);
	  BM = 0.0;
	  //std::cout << "+-+ " << A  << " " << BP << " " << BM << " " << std::endl;
	}
	if ((l1 == -1) && (l2 == +1) && (l3 == +1)) {
	  A  = oSpinor.APPM (3, 4, 2);
	  BP = oSpinor.BPPM (3, 4, 2);
	  BM = 0.0;
	  //std::cout << "-++ " << A  << " " << BP << " " << BM << " " << std::endl;
	}
	if ((l1 == +1) && (l2 == -1) && (l3 == -1)) {
	  A  = oSpinor.APMM (2, 3, 4);
	  BP = oSpinor.BPMM (2, 3, 4);
	  BM = 0.0;
	  //std::cout << "+-- " << A  << " " << BP << " " << BM << " " << std::endl;
	}
	if ((l1 == -1) && (l2 == -1) && (l3 == +1)) {
	  A  = oSpinor.APMM (4, 2, 3);
	  BP = oSpinor.BPMM (4, 2, 3);
	  BM = 0.0;
	  //std::cout << "--+ " << A  << " " << BP << " " << BM << " " << std::endl;
	}
	if ((l1 == -1) && (l2 == +1) && (l3 == -1)) {
	  A  = oSpinor.APMM (3, 2, 4);
	  BP = oSpinor.BPMM (3, 2, 4);
	  BM = 0.0;
	  //std::cout << "-+- " << A  << " " << BP << " " << BM << " " << std::endl;
	}
//! COUPLAGES
	A  = A  * oParam.GA;
	BP = BP * oParam.GBP;
	BM = BM * oParam.GBM;
	//std::cout << "" << A  << " " << BP << " " << BM << " " << std::endl;
//! LES DIFFERENTS TERMES DE L'ELEMENT DE MATRICE AU CARRE
	M2 [LL1] [LL2] [LL3] [0] = pow (abs (A ), 2);
	M2 [LL1] [LL2] [LL3] [1] = pow (abs (BP), 2);
	M2 [LL1] [LL2] [LL3] [2] = pow (abs (BM), 2);
	M2 [LL1] [LL2] [LL3] [3] = 2.0 * real (A * conj (BP));
	M2 [LL1] [LL2] [LL3] [4] = 2.0 * imag (A * conj (BP));
	M2 [LL1] [LL2] [LL3] [5] = 0.0;
	M2 [LL1] [LL2] [LL3] [6] = 0.0;
	M2 [LL1] [LL2] [LL3] [7] = 0.0;
      }
    }
  }
}

void result::display (param oParam) {
  if (oParam.IMPR) {
    for (int K = 0; K < NRESUL; K++) {
      std::cout << K << std::endl;
      std::cout << "  ---     --+     -+-     +--     -++     +-+     ++-     +++" << std::endl;
      //WRITE (*, 900)
      //900 FORMAT (8E8.2)
      std::cout << 
	M2 [0] [0] [0] [K] << " " << M2 [0] [0] [1] [K] << " " << 
	M2 [0] [1] [0] [K] << " " << M2 [1] [0] [0] [K] << " " << 
	M2 [0] [1] [1] [K] << " " << M2 [1] [0] [1] [K] << " " << 
	M2 [1] [1] [0] [K] << " " << M2 [1] [1] [1] [K] << " " << std::endl;
    }
  }
}
