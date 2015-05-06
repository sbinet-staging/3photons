#ifndef PPP_H
#define PPP_H
#include "precision.h"
#include <iostream>

static const int INP = 3;    // Number of impulsions generated
class ppp
{
//!   Common pour passer les impulsions engendrees par Rambo

public:
  double POUT [4][100], P1 [4], P2 [4];

	// Accesseurs
	// Constructeurs
  ppp (double ETOT) {   //!   Initialisation des impulsions entrantes
    P1 [3] =  ETOT/2.0;
    P2 [3] =  ETOT/2.0;
    P1 [0] = -ETOT/2.0;
    P2 [0] =  ETOT/2.0;
    P1 [1] = 0.0;
    P1 [2] = 0.0;
    P2 [1] = 0.0;
    P2 [2] = 0.0;
  }

	// Destructeur
	// Méthodes
  void RAMBO (const int N, double ET, double* XM, double& WT);
  //void RAMBO (const int N, double ET, double* XM, double** P, double& WT);

  void display () {
    for (int I = 0; I < 4; I++) {
      std::cout << I << " "<< POUT [I] [0] << " " << POUT [I] [1] << " " << POUT [I] [2] << std::endl;
    }
    std::cout << std::endl;
  };

  //!   TRIE LES PHOTONS SORTANTS PAR ENERGIE
  void TRI () {
    double PE [4];

    if (POUT [3] [1] > POUT [3] [0]) {
      for (int I = 0; I < 4; I++) {
	PE [I] = POUT [I] [0];
	POUT [I] [0] = POUT [I] [1];
	POUT [I] [1] = PE [I];
      }
    }
    if (POUT [3] [2] > POUT [3] [0]) {
      for (int I = 0; I < 4; I++) {
	PE [I] = POUT [I] [0];
	POUT [I] [0] = POUT [I] [2];
	POUT [I] [2] = PE [I];
      }
    }
    if (POUT [3] [2] > POUT [3] [1]) {
      for (int I = 0; I < 4; I++) {
	PE [I]=POUT [I] [1];
	POUT [I] [1] = POUT [I] [2];
	POUT [I] [2] = PE [I];
      }
    }
  }

};
#endif // PPP_H
