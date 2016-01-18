#ifndef ANGLE_H
#define ANGLE_H
#include "precision.h"
#include "ppp.h"
#include "scalar.h"
#include "cutpar.h"

//   Common pour passer les angles
//! \brief set of angles between outgoing particles, beam and plane defined by the outgoing particles
//! \author Vincent C. LAFAGE
//! \date 2007-08-11 ISO
class angle
{
  double
    COS1P1, //!< cosine of angle between highest energy photon and electron beam
    COS2P1, //!< cosine of angle between middle energy photon and electron beam
    COS3P1, //!< cosine of angle between lowest energy photon and electron beam
    COS12,  //!< cosine of angle between highest energy photon and middle energy photon
    COS13,  //!< cosine of angle between highest energy photon and lowest energy photon
    COS23,  //!< cosine of angle between middle energy photon and lowest energy photon
    COSN,   //!< cosine of angle between outgoing photons plane and electron beam
    COSAC;  //!< cosine of angle between outgoing photons perpendicular momenta

public:

	// Accesseurs
	//! Constructor
  angle (ppp oPpp, scalar oScalar);

	// Destructeur
	// Méthodes
  bool CUT (ppp oPpp, cutpar oCutpar);
};
#endif // ANGLE_H
