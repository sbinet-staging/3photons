with Precision;
with Spinor;
with Ppp;

package body Angle is
   use Precision;
   use Specific_Elementary;

   function oAngle (OPpp : in Ppp.Ppp; OScalar : in Scalar.Scalar) return Angle is
      use Spinor;

      P : array (Spinor.T_INDICE_PARTICULE, T_INDICE_MINK) of Real;
      NX, NY, NZ, NN : Real;
      OAngle : Angle;
   begin

      --for J in 1..4 loop
      for J in T_INDICE_MINK loop
         P (1, J) := oPpp.P1 (J);
         P (2, J) := oPpp.P2 (J);
         P (3, J) := oPpp.POUT (J, 1);
         P (4, J) := oPpp.POUT (J, 2);
         P (5, J) := oPpp.POUT (J, 3);
      end loop;

      oAngle.Cos1p1 := 1.0-oScalar.PS (1, 3)/ (P (1, 4) * P (3, 4));
      oAngle.Cos2p1 := 1.0-oScalar.PS (1, 4)/ (P (1, 4) * P (4, 4));
      oAngle.Cos3p1 := 1.0-oScalar.PS (1, 5)/ (P (1, 4) * P (5, 4));

      oAngle.Cos12 := 1.0-oScalar.PS (3, 4)/ (oPpp.POUT (4, 1)*oPpp.POUT (4, 2));
      oAngle.Cos13 := 1.0-oScalar.PS (3, 5)/ (oPpp.POUT (4, 1)*oPpp.POUT (4, 3));
      oAngle.Cos23 := 1.0-oScalar.PS (4, 5)/ (oPpp.POUT (4, 2)*oPpp.POUT (4, 3));

      NX := oPpp.POUT (2, 1)*oPpp.POUT (3, 2)-oPpp.POUT (2, 2)*oPpp.POUT (3, 1);
      NY := oPpp.POUT (3, 1)*oPpp.POUT (1, 2)-oPpp.POUT (3, 2)*oPpp.POUT (1, 1);
      NZ := oPpp.POUT (1, 1)*oPpp.POUT (2, 2)-oPpp.POUT (1, 2)*oPpp.POUT (2, 1);
      NN := Sqrt (NX**2+NY**2+NZ**2);
      oAngle.COSN :=  (NX*P (1, 1)+NY*P (1, 2)+NZ*P (1, 3))/P (1, 4)/NN;

      oAngle.COSAC :=  (P (3, 2)*P (4, 2)+P (3, 3)*P (4, 3)) / Specific_Elementary.Sqrt (
         (P (3, 2) * P (3, 2) + P (3, 3) * P (3, 3)) *
         (P (4, 2) * P (4, 2) + P (4, 3) * P (4, 3)));
      return OAngle;
   end OAngle;

  function CUT (OPpp : Ppp.Ppp; OCutpar : Cutpar.Cutpar; OAngle : Angle) return Boolean is
     ------------------------------------------------------------------------
     -- Determine si un evenement passe ou non les coupures
     ------------------------------------------------------------------------
     CUT : Boolean := False;
  begin
     -- Coupure s'il y a lieu
     for I in 1..INP loop
        CUT := CUT or (oPpp.POUT (4, I) < oCutpar.EMIN);
     end loop;
     CUT := CUT or
       (abs oAngle.COS1P1 > oCutpar.ACUT) or
       (abs oAngle.COS2P1 > oCutpar.ACUT) or
       (abs oAngle.COS3P1 > oCutpar.ACUT) or
       (oAngle.COS12 > oCutpar.BCUT) or
       (oAngle.COS13 > oCutpar.BCUT) or
       (oAngle.COS23 > oCutpar.BCUT) or
       (abs oAngle.COSN < oCutpar.SINCUT);
     return CUT;
  end CUT;

end Angle;
