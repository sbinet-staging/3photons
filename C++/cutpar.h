//
//   Common pour passer les coupures
//
#ifndef CUTPAR_H
#define CUTPAR_H
#include "precision.h"
 
//! \brief set of experimental cuts
//! \author Vincent C. LAFAGE
//! \date 2007-08-04 ISO
class cutpar
{
public:
  double
    ACUT,   //!< cut on maximum cosine of (beam, photons) angle
    BCUT,   //!< cut on maximum cosine of (photon, photon) angle
    EMIN,   //!< cut on minimum photon energy
    SINCUT; //!< cut on minimum cosine of (beam, normal to the photon plane) angle
};
#endif // CUTPAR_H
