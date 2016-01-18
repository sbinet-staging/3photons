#include "intel_compatibility.h"
#include "precision.h"
#include "resfin.h"
#include <iostream>
#include <math.h>

#include <complex>

#include "ppp.h"


//! Affiche les resultats dans le parametrage d'Eric
//void ERIC (int UNE, double BREPEM, double CONVERS, double PI) {
//void resfin::ERIC (ofstream UNE, double BREPEM, double CONVERS, double PI) {
void resfin::ERIC (double BREPEM, double CONVERS, double PI, param oParam) {
  double MUTH;

  MUTH = BREPEM/(8*9*5*pow (PI, 2) * oParam.MZ0 * oParam.GZ0)*CONVERS;
  std::cout << std::endl;
  std::cout <<     "       :        -          +" << std::endl;
  //950 FORMAT(A, 2E12.4) // pour toute la suite
  std::cout << "sigma0  : " <<   SPM2 [0] [0]/2.0 << " | " <<  SPM2 [1] [0]/2.0 << std::endl;
  std::cout << "alpha0  : " <<   SPM2 [0] [4]/2.0 << " | " <<  SPM2 [1] [4]/2.0 << std::endl;
  std::cout << "beta0   : " <<  -SPM2 [0] [3]/2.0 << " | " << -SPM2 [1] [3]/2.0 << std::endl;
  std::cout << "lambda0 : " << (-SPM2 [0] [1]+SPM2 [0] [2])/2.0 << " | " << (-SPM2 [1] [1]+SPM2 [1] [2])/2.0 << std::endl;
  std::cout << "mu0     : " <<  (SPM2 [0] [1]+SPM2 [0] [2])/2.0 << " | " <<  (SPM2 [1] [1]+SPM2 [1] [2])/2.0 << std::endl;
  std::cout << "mu/lamb : " <<  (SPM2 [0] [1]+SPM2 [0] [2])/(-SPM2 [0] [1]+SPM2 [0] [2]) << " | " << (SPM2 [1] [1]+SPM2 [1] [2])/(-SPM2 [1] [1]+SPM2 [1] [2]) << std::endl;
  std::cout << "mu (num): " <<  (SPM2 [0] [1]+SPM2 [0] [2]+SPM2 [1] [1]+SPM2 [1] [2])/4.0 << std::endl;
  std::cout << "rapport : " <<  (SPM2 [0] [1]+SPM2 [0] [2]+SPM2 [1] [1]+SPM2 [1] [2])/4.0/MUTH << std::endl;
  std::cout << "mu (th) : " << MUTH << std::endl;
}

//void resfin::FAWZI (ofstream UNE, double BREPEM, double CONVERS, double PI, double ETOT, param oParam, cutpar oCutpar) {
void resfin::FAWZI (double BREPEM, double CONVERS, double PI, double ETOT, param oParam, cutpar oCutpar) {
//! Displays Fawzi's analytical results, 
//! and compare them to MonteCarlo results

  double BRA, MRE, GRE, SIG, DEL, EPS;
  double F1, G1, G2, G3, FF, GG;
  double SIGP, SIGM, MCP, MCM, INCRP, INCRM;
  std::complex<double> SDZ;

  MRE = oParam.MZ0/ETOT;
  GRE = oParam.GZ0*oParam.MZ0/pow (ETOT, 2);
  SDZ = std::complex<double> (1.0-pow (MRE, 2), -GRE)/(pow ((1.0-pow (MRE, 2)), 2)+pow (GRE, 2));
  DEL = (1.0-oCutpar.BCUT)/2.0;
  EPS = 2.0*oCutpar.EMIN/ETOT;
  BRA = oParam.MZ0/3.0/6.0/pow (PI, 3)/16.0/120.0;
  //!SIG=12.0*PI/oParam.MZ0**2*(BREPEM*oParam.GZ0)*BRA/ETOT**2*(ETOT/oParam.MZ0)**8
  //!       *CDABS(SDZ)**2*CONVERS
  SIG=12.0*PI/pow (oParam.MZ0, 2)*(BREPEM*oParam.GZ0)*BRA/pow (ETOT, 2)*pow (ETOT/oParam.MZ0, 8) * pow (abs (SDZ), 2) * CONVERS;

  F1=           1.0 - 15.0*pow (EPS, 4)
    - 9.0/ 7.0*(1.0 - 70.0*pow (EPS, 4)) * pow (DEL, 2)
    + 6.0/ 7.0*(1.0 + 70.0*pow (EPS, 4)) * pow (DEL, 3);

  G1=           1.0 - 30.0 * pow (EPS, 4) 
    - 9.0/ 7.0*(1.0 - 70.0 * pow (EPS, 4)) * DEL
    - 90.0                 * pow (EPS, 4)  * pow (DEL, 2)
    - 1.0/ 7.0*(1.0 -420.0 * pow (EPS, 4)) * pow (DEL, 3);

  G2=           1.0 - 25.0 * pow (EPS, 4)
    - 6.0/ 7.0*(1.0 - 70.0 * pow (EPS, 4))  * DEL
    - 3.0/ 7.0*(1.0 +210.0 * pow (EPS, 4))  * pow (DEL, 2)
    - 8.0/21.0*(1.0 -105.0/2.0*pow (EPS, 4))* pow (DEL, 3);

  G3=       1.0-195.0/11.0*pow (EPS, 4)
    -18.0/77.0*(1.0- 7.0*pow (EPS, 4))       * DEL
    - 9.0/11.0*(9.0/7.0-70.0*pow (EPS, 4))   * pow (DEL, 2)
    - 8.0/11.0*(1.0-105.0/11.0*pow (EPS, 4)) * pow (DEL, 3);

  FF=F1*(1.0-pow (oCutpar.SINCUT, 3));
  GG=G1-27.0/16.0*G2*oCutpar.SINCUT+11.0/16.0*G3*pow (oCutpar.SINCUT, 3);

  SIGP=SIG*(FF+2.0*GG);
  SIGM=SIG*(FF+4.0*GG);

  MCP = (SPM2 [0] [1]+SPM2 [1] [1])/4.0;
  MCM = (SPM2 [0] [2]+SPM2 [1] [2])/4.0;
  INCRP =sqrt (pow (SPM2 [0] [1] * VAR [0] [1], 2)+
	       pow (SPM2 [1] [1] * VAR [1] [1], 2))
    /fabs (SPM2 [0] [1]+SPM2 [1] [1]);
  INCRM =sqrt (pow (SPM2 [0] [2] * VAR [0] [2], 2)+
	       pow (SPM2 [1] [2] * VAR [1] [2], 2))
             /fabs (SPM2 [0] [2] +SPM2 [1] [2]);

  std::cout << std::endl;
  std::cout << "s (pb) :   Sig_cut_Th    Sig_Th      Rapport" << std::endl;
  std::cout << "       :   Sig_Num" << std::endl;
  std::cout << "       :   Ecart_relatif  Incertitude" << std::endl;
  std::cout << std::endl;
  std::cout << "s+(pb) : " << SIGP << " | " << SIG*3.0 << " | " << SIGP/SIG/3.0 << std::endl;
  std::cout << "       : " << MCP << std::endl;
  std::cout << "       : " << MCP/SIGP-1.0 << " | " << INCRP << " | " << (MCP/SIGP-1.0)/INCRP << std::endl;
  std::cout << std::endl;
  std::cout << "s-(pb) : " << SIGM << " | " << SIG*5.0 << " | " << SIGM/SIG/5.0 << std::endl;
  std::cout << "       : " << MCM << std::endl;
  std::cout << "       : " << MCM/SIGM-1.0 << " | " << INCRM << " | " << (MCM/SIGM-1.0)/INCRM << std::endl;
  std::cout << std::endl;

}
