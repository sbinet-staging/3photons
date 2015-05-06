#include "intel_compatibility.h"
//g++ -Wall -pedantic hasa.cxx -o tsthasa
#include "precision.h"
#include "hasa.h"

static const long MODULO = 1000000000;
long static NCALL = 0;
long static MCALL = 55;

//--------------------------------------------------------------------------
double RN () //FUNCTION RN (IDMY)
//--------------------------------------------------------------------------
{
//    ***************************************************************
//    * RANDOM NUMBER FUNCTION TAKEN FROM KNUTH RANF                *
//    * (SEMINUMERICAL ALGORITHMS).                                 *
//    * METHOD IS X (N)=MOD (X (N-55)-X (N-24), 1/FMODUL)           *
//    * NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.          *
//    *                                                             *
//    * RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.               *
//    * IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.   *
//    * IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE. *
//    ***************************************************************
  long static IA [56];
  double const FMODUL = 1.0e-9;

  if (NCALL==0) {
    IN55 (IA, (long) 234612947);
    NCALL = 1;
  }
  if (MCALL==0) {
    IRN55 (IA);
    MCALL = 55;
  }
  MCALL = MCALL - 1;
  return ((double) IA [MCALL+1]) * FMODUL;
}

//--------------------------------------------------------------------------
void IN55 (long IA[], long IX) //SUBROUTINE IN55 (IA, IX)
//--------------------------------------------------------------------------
{
  IA [55]=IX;
  long J  = IX;
  long K  = 1;
  for (int I=1; I <= 54; I++) {
    long II = (21*I) % 55;
    IA [II] = K;
    K = J - K;
    if (K < 0) K += MODULO;
    J = IA [II];
  }
  for (int I=0; I < 10; I++) {
    IRN55 (IA);
  }
}

//--------------------------------------------------------------------------
void IRN55 (long IA[])
//--------------------------------------------------------------------------
{
  long J;
  for (int I=1; I <= 24; I++) {
    J = IA [I] - IA [I+31];
    if (J < 0) J += MODULO;
    IA [I] = J;
  }
  for (int I=25; I <= 55; I++) {
    J = IA [I] - IA [I-24];
    if (J < 0) J += MODULO;
    IA [I] = J;
  }
}
