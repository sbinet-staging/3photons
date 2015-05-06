//
//   Common pour passer les produits scalaires engendres par mat3
//
#ifndef SCALAR_H
#define SCALAR_H
#include "precision.h"
#include "spinor.h"
#include <complex>

class scalar
{
public:
  double PS [5][5];

	// Accesseurs
	// Constructeurs
  scalar (spinor oSpinor) {
    std::complex<double> CX;

    for (int K = 0; K < 5; K++){
      PS [K] [K] = 0.0;
    }

    //! PRODUIT SCALAIRE DE LORENTZ
    for (int J = 0; J < 4; J++){
      for (int K = J+1; K < 5; K++){
	CX = oSpinor.S [J] [K];
	PS [J] [K] = real (CX*conj (CX))/2.0;
	PS [K] [J] = PS [J] [K];
      }
    }
  }
	// Destructeur
	// Méthodes
};
#endif // SCALAR_H
