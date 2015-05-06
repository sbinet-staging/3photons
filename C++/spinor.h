#ifndef SPINOR_H
#define SPINOR_H
#include "precision.h"
#include "ppp.h"
#include <math.h>
#include <complex>

const double RAC8 =2.0 * sqrt (2.0);
class spinor
{
//!   Common pour passer les produits spinoriels
public:
  std::complex<double> S [5] [5], T [5] [5];

	// Accesseurs
	// Constructeurs
  spinor (ppp oPpp);
	// Destructeur
	// Méthodes
  std::complex<double> APPM (int k1, int k2, int k3);
  std::complex<double> APMM (int k1, int k2, int k3);
  std::complex<double> BPPM (int k1, int k2, int k3);
  std::complex<double> BPMM (int k1, int k2, int k3);
  std::complex<double> BPPP (int k1, int k2, int k3);
  std::complex<double> BMMM (int k1, int k2, int k3);
};
#endif // SPINOR_H
