#ifndef RESULT_H
#define RESULT_H
#include "precision.h"

#include "param.h"
#include "spinor.h"

static const int NRESUL = 8;
class result
{
//!   Common pour passer les resultats engendrees par mat(rice)
public:
  double M2 [2] [2] [2] [NRESUL];

	// Constructeurs
  result (param oParam, spinor oSpinor, double ETOT);
	// Accesseurs
  void display (param oParam);
	// Destructeur
	// Méthodes
};
#endif // RESULT_H
