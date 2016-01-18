#ifndef RESULT_H
#define RESULT_H
#include "precision.h"

#include "param.h"
#include "spinor.h"

static const int NRESUL = 8; //!< number of type of results
//! \brief array of squared matrix elements contribution with detail of helicities
//! \author Vincent C. LAFAGE
//! \date 2007-08-11 ISO
class result
{
//!   Array to pass results computed by mat(rix)
public:
  double M2 [2] [2] [2] [NRESUL];  //!< array of squared matrix elements NRESUL contribution with detail of outgoing helicities configuration

	// Constructeurs
  result (param oParam, spinor oSpinor, double ETOT);
	// Accesseurs
  //! \brief ordered dump of squared modulus transition Matrix: all contribution, all helicity configuration
  void display (param oParam);
	// Destructeur
	// Méthodes
};
#endif // RESULT_H
