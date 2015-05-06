#ifndef ANGLE_H
#define ANGLE_H
#include "precision.h"
#include "ppp.h"
#include "scalar.h"
#include "cutpar.h"

//   Common pour passer les angles
class angle
{
  double COS12, COS1P1, COS2P1, COS3P1, COS13, COS23, COSN, COSAC;
public:

	// Accesseurs
	// Constructeurs
  angle (ppp oPpp, scalar oScalar);

	// Destructeur
	// Méthodes
  bool CUT (ppp oPpp, cutpar oCutpar);
};
#endif // ANGLE_H
