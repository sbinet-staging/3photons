#ifndef PARAM_H
#define PARAM_H
#include "precision.h"

//! \brief physical constants and parameter
//! \author Vincent C. LAFAGE
//! \date 2007-08-04 ISO
class param
{

public:
  double
    MZ0,     //!< Z⁰ boson mass (GeV)
    GZ0,     //!< Z⁰ boson width (GeV)
    GA,      //!< Standard Model contribution electromagnetic coupling √(4𝜋𝛼)³
    GBP,     //!< 𝛽₊ anomalous contribution electroweak coupling
    GBM,     //!< 𝛽₋ anomalous contribution electroweak coupling
    POLP,    //!< electroweak polarisations factors for 𝛽₊ anomalous contribution
    POLM,    //!< electroweak polarisations factors for 𝛽₋ anomalous contribution
    POLP2,   //!< electroweak polarisations factors for 𝛽₊ anomalous contribution
    POLM2;   //!< electroweak polarisations factors for 𝛽₋ anomalous contribution
  bool IMPR; //!< boolean predicate value controlling dump of result

	// Accesseurs
	// Constructeurs
	// Destructeur
	// Méthodes
};
#endif // PARAM_H
