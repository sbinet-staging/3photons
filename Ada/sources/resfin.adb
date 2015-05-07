with Text_IO;
with Precision;
with Param;
with Cutpar;
use  Precision;
use  Param;
use  Cutpar;

package body Resfin is
   use Text_IO;
   use Specific_Complex;
   package Real_IO    is new Text_IO.Float_IO (Real);
   use  Real_IO;
   package Entier_IO is new Text_IO.Integer_IO (Integer);
   use  Entier_IO;

   procedure ERIC (UNE : in Integer; BREPEM, CONVERS, PI : in Real; OParam : in Param.Param; OResfin : in Resfin) is
      ------------------------------------------------------------------------
      -- Affiche les resultats dans le parametrage d'Eric
      ------------------------------------------------------------------------
      MUTH : Real;
   begin
      MUTH := BREPEM/(8.0*9.0*5.0*PI**2*oParam.MZ0*oParam.GZ0)*CONVERS;
      Put_Line ("");
      Put_Line ("        :        -          +");
      --write (UNE, 950)--950 format(A, 2E12.4)
      Put ("sigma0  : "); Put (  oResfin.SPM2 (1, 1)/2.0, 3, 4, 3); Put ( oResfin.SPM2 (2, 1)/2.0, 3, 4, 3); Put_Line ("");
      Put ("alpha0  : "); Put (  oResfin.SPM2 (1, 5)/2.0, 3, 4, 3); Put ( oResfin.SPM2 (2, 5)/2.0, 3, 4, 3); Put_Line ("");
      Put ("beta0   : "); Put ( -oResfin.SPM2 (1, 4)/2.0, 3, 4, 3); Put (-oResfin.SPM2 (2, 4)/2.0, 3, 4, 3); Put_Line ("");
      Put ("lambda0 : "); Put ((-oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3))/2.0, 3, 4, 3);
                          Put ((-oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/2.0, 3, 4, 3); Put_Line ("");
      Put ("mu0     : "); Put (( oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3))/2.0, 3, 4, 3);
                          Put (( oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/2.0, 3, 4, 3); Put_Line ("");
      Put ("mu/lamb : ");
      Put ((oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3))/(-oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3)), 3, 4, 3);
      Put ((oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/(-oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3)), 3, 4, 3);
      Put_Line ("");
      Put ("mu (th) : "); Put (MUTH, 3, 4, 3); Put_Line ("");
      Put ("mu (num): "); Put ((oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/4.0, 3, 4, 3); Put_Line ("");
      Put ("rapport : "); Put ((oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/4.0/MUTH, 3, 4, 3); Put_Line ("");
   end ERIC;

   procedure FAWZI (UNE : in Integer; BREPEM, CONVERS, PI, ETOT : in Real;
                                      OParam  : in Param.Param;
                                      OCutpar : in Cutpar.Cutpar;
                                      OResfin : in Resfin) is
      ------------------------------------------------------------------------
      -- Affiche les resultats analytiques de Fawzi,
      -- et compare aux resultats du MC
      ------------------------------------------------------------------------
      use Cutpar;
      use Precision;
      use Specific_Elementary;
      BRA, MRE, GRE, SIG, DEL, EPS : Real;
      F1, G1, G2, G3, FF, GG : Real;
      SIGP, SIGM, MCP, MCM, INCRP, INCRM : Real;
      SDZ : Complex;
   begin
      MRE := oParam.MZ0/ETOT;
      GRE := oParam.GZ0*oParam.MZ0/ETOT**2;
      SDZ := (1.0-MRE**2, -GRE)
           /((1.0-MRE**2)**2+GRE**2);
      DEL := (1.0-oCutpar.BCUT)/2.0;
      EPS := 2.0*oCutpar.EMIN/ETOT;
      BRA := oParam.MZ0/3.0/6.0/PI**3/16.0/120.0;
      --SIG := 12.0*PI/MZ0**2*(BREPEM*GZ0)*BRA/ETOT**2*(ETOT/MZ0)**8 &
      --       *CDABS(SDZ)**2*CONVERS
      SIG := 12.0*PI/oParam.MZ0**2*(BREPEM*oParam.GZ0)*BRA/ETOT**2*(ETOT/oParam.MZ0)**8 * (abs SDZ)**2*CONVERS;

      F1 :=            1.0- 15.0*EPS**4
        - 9.0/ 7.0*(1.0- 70.0*EPS**4)*DEL**2
        + 6.0/ 7.0*(1.0+ 70.0*EPS**4)*DEL**3;

      G1 :=         1.0- 30.0*EPS**4
        - 9.0/ 7.0*(1.0- 70.0*EPS**4)*DEL
                       - 90.0*EPS**4 *DEL**2
        - 1.0/ 7.0*(1.0-420.0*EPS**4)*DEL**3;

      G2 :=         1.0     -  25.0 * EPS**4
        - 6.0/ 7.0*(1.0     -  70.0 * EPS**4) * DEL
        - 3.0/ 7.0*(1.0     + 210.0 * EPS**4) * DEL**2
        - 8.0/21.0*(1.0 - 105.0/2.0 * EPS**4) * DEL**3;

      G3 :=         1.0 - 195.0/11.0 * EPS**4
        -18.0/77.0*(1.0 -        7.0 * EPS**4) * DEL
        - 9.0/11.0*(9.0/7.0    -70.0 * EPS**4) * DEL**2
        - 8.0/11.0*(1.0 - 105.0/11.0 * EPS**4) * DEL**3;

      FF := F1*(1.0- oCutpar.SINCUT**3);
      GG := G1-27.0/16.0*G2*oCutpar.SINCUT+11.0/16.0*G3* oCutpar.SINCUT**3;

      SIGP := SIG*(FF+2.0*GG);
      SIGM := SIG*(FF+4.0*GG);

      MCP := (oResfin.SPM2 (1, 2)+oResfin.SPM2 (2, 2))/4.0;
      MCM := (oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 3))/4.0;
      INCRP  :=  sqrt ((oResfin.SPM2 (1, 2)* oResfin.VAR (1, 2))**2+(oResfin.SPM2 (2, 2) * oResfin.VAR (2, 2))**2)
        /abs (oResfin.SPM2 (1, 2)+oResfin.SPM2 (2, 2));
      INCRM  :=  sqrt ((oResfin.SPM2 (1, 3) * oResfin.VAR (1, 3))**2+(oResfin.SPM2 (2, 3) * oResfin.VAR (2, 3))**2)
        /abs (oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 3));

      Put_Line ("");
      Put_Line ("s (pb) :   Sig_cut_Th    Sig_Th      Rapport");
      Put_Line ("       :   Sig_Num");
      Put_Line ("       :   Ecart_relatif  Incertitude");
      Put_Line ("");
      --Put ("s+(pb) : "); Put (SIGP, 2, 8, 3); Put (SIG*3.0, 2, 8, 3); Put (SIGP/SIG/3.0, 2, 4, 3); Put_Line ("");
      Put ("s+(pb) : "); Put (SIGP, 2, 7); Put (SIG*3.0, 2, 7); Put (SIGP/SIG/3.0, 2, 3); Put_Line ("");
      Put ("       : "); Put (MCP, 2, 8, 3);Put_Line ("");
      Put ("       : "); Put (MCP/SIGP-1.0, 2, 8, 3); Put (INCRP, 2, 8, 3); Put ((MCP/SIGP-1.0)/INCRP, 2, 4, 3);
      Put_Line ("");
      Put_Line ("");
      Put ("s-(pb) : "); Put (SIGM, 2, 8, 3); Put (SIG*5.0, 2, 8, 3); Put (SIGM/SIG/5.0, 2, 4, 3);Put_Line ("");
      Put ("       : "); Put (MCM, 2, 8, 3); Put_Line ("");
      Put ("       : "); Put (MCM/SIGM-1.0, 2, 8); Put (INCRM, 2, 8); Put ((MCM/SIGM-1.0)/INCRM, 2, 4);
      --Put ("       : "); Put (MCM/SIGM-1.0, 2, 8, 3); Put (INCRM, 2, 8, 3); Put ((MCM/SIGM-1.0)/INCRM, 2, 4, 3);
      Put_Line ("");
      Put_Line ("");
   end FAWZI;
end Resfin;


      --write (UNE, 950) 'sigma0  : ', oResfin.SPM2 (1, 1)/2.0, oResfin.SPM2 (2, 1)/2.0
      --write (UNE, 950) 'alpha0  : ', oResfin.SPM2 (1, 5)/2.0, oResfin.SPM2 (2, 5)/2.0
      --write (UNE, 950) 'beta0   : ', -oResfin.SPM2 (1, 4)/2.0, -oResfin.SPM2 (2, 4)/2.0
      --write (UNE, 950) 'lambda0 : ', (-oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3))/2.0, &
      --     (-oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/2.0
      --write (UNE, 950) 'mu0     : ', (oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3))/2.0, &
      --     (oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/2.0
      --write (UNE, 950) 'mu/lamb : ',  &
      --     (oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3))/(-oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3)), &
      --     (oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/(-oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))
      --write (UNE, 950) 'mu (th) : ', MUTH
      --write (UNE, 950) 'mu (num): ',  &
      --     (oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/4.0
      --write (UNE, 950) 'rapport : ',  &
      --     (oResfin.SPM2 (1, 2)+oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 2)+oResfin.SPM2 (2, 3))/4.0/MUTH


      --7   format (A, 2G14.8, G10.4)
      --8   format (A,  G14.8)
      --9   format (A, 2G14.8, G10.4)

      --Put_Line ("");
      --6   format (A)
      --write (UNE, 6)'s (pb) :   Sig_cut_Th    Sig_Th      Rapport'
      --write (UNE, 6)'       :   Sig_Num'
      --write (UNE, 6)'       :   Ecart_relatif  Incertitude'
      --Put_Line ("");
      --write (UNE, 7)'s+(pb) : ', SIGP, SIG*3.0, SIGP/SIG/3.0
      --write (UNE, 8)'       : ', MCP
      --write (UNE, 9)'       : ', MCP/SIGP-1.0, INCRP, (MCP/SIGP-1.0)/INCRP
      --Put_Line ("");
      --write (UNE, 7)'s-(pb) : ', SIGM, SIG*5.0, SIGM/SIG/5.0
      --write (UNE, 8)'       : ', MCM
      --write (UNE, 9)'       : ', MCM/SIGM-1.0, INCRM, (MCM/SIGM-1.0)/INCRM
      --Put_Line ("");
