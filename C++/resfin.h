#ifndef RESFIN_H
#define RESFIN_H
#include "precision.h"
#include "param.h"
#include "cutpar.h"

const int NRES = 8 ;

//!> \brief array of final results: differential cross-section, sum and variance
//!> \author Vincent C. LAFAGE
//!> \date 2007-08-22 ISO
class resfin
{

public:
  double
    SPM2DIF [NRES],  //!< element of cross-section
    SPM2 [2] [NRES], //!< total cross-section
    VAR [2] [NRES];  //!< variance of the sum

	// Constructeurs
  resfin () {
    for (int Sp=0; Sp < 2; Sp++) {
      for (int K=0; K < NRES; K++) {
	SPM2 [Sp] [K] = 0.0;
	VAR  [Sp] [K] = 0.0;
      }
    }
    for (int K=0; K < NRES; K++) {
      SPM2DIF [K] = 0.0;
    }
  }
	// Destructeur
	// Accesseurs
  void ERIC (double BREPEM, double CONVERS, double PI, param oParam);
  //void ERIC (ofstream UNE, double BREPEM, double CONVERS, double PI);
  //void FAWZI (ofstream UNE, double BREPEM, double CONVERS, double PI, double ETOT, param oParam, cutpar oCutpar);
  void FAWZI (double BREPEM, double CONVERS, double PI, double ETOT, param oParam, cutpar oCutpar);
	// Méthodes
};
#endif // RESFIN_H
