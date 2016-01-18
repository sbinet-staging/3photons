#ifndef SPINOR_H
#define SPINOR_H
#include "precision.h"
#include "ppp.h"
#include <math.h>
#include <complex>

const double RAC8 =2.0 * sqrt (2.0);
//! \brief array of massless momenta spinor products
//! \author Vincent C. LAFAGE
//! \date 2007-08-23 ISO
class spinor
{
public:
  std::complex<double>
    S [5] [5],  //!< massless momenta spinor inner products Gram matrix
    T [5] [5];  //!< massless momenta conjugate spinor inner products Gram matrix

	// Accesseurs
	// Constructeurs
  spinor (ppp oPpp);
	// Destructeur
	// Méthodes
  std::complex<double> APPM (int k1, int k2, int k3); //!< computes standard amplitude for helicities ++-
  std::complex<double> APMM (int k1, int k2, int k3); //!< computes standard amplitude for helicities +--
  std::complex<double> BPPM (int k1, int k2, int k3); //!< computes anomalous amplitude for helicities ++-
  std::complex<double> BPMM (int k1, int k2, int k3); //!< computes anomalous amplitude for helicities +--
  std::complex<double> BPPP (int k1, int k2, int k3); //!< computes anomalous amplitude for helicities +++
  std::complex<double> BMMM (int k1, int k2, int k3); //!< computes anomalous amplitude for helicities ---
};
#endif // SPINOR_H
