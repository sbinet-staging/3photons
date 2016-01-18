#ifndef XTRPRO_H
#define XTRPRO_H
#include "precision.h"

//!> \brief set of extra probabilities (for histogramming)
//!> \author Vincent C. LAFAGE
//!> \date 2007-08-04 ISO
class xtrpro
{

public:
  float
    PRPLUS,  //!<  extra probabilities (for histogramming): not used in C++, no PAW
    PRMOINS; //!<  extra probabilities (for histogramming): not used in C++, no PAW
  double EE1, EE2;
// present dans mc.f (naturellement)
// present dans book.f

	// Accesseurs
	// Constructeurs
	// Destructeur
	// Méthodes
};
#endif // XTRPRO_H
