with Ppp;
with Precision;

package body Spinor is
   use Specific_Elementary;
   function oSpinor (OPpp : in Ppp.Ppp) return Spinor is
      OSpinor : Spinor;
      P : array (T_INDICE_PARTICULE, T_INDICE_MINK) of Real;
      TX : Real;
      XX : array (T_INDICE_PARTICULE) of Real;
      CX : Complex;
      FX : array (T_INDICE_PARTICULE) of Complex;
   begin
      -- Transfere dans le constructeur de l'objet Angle
      for J in T_INDICE_MINK loop --1..4
         P (1, J) := oPpp.P1 (J);
         P (2, J) := oPpp.P2 (J);
         P (3, J) := oPpp.POUT (J, 1);
         P (4, J) := oPpp.POUT (J, 2);
         P (5, J) := oPpp.POUT (J, 3);
      end loop;

      for K in T_INDICE_PARTICULE loop -- 1..5
         OSpinor.S (K, K) := (0.0, 0.0);
         TX := sqrt (P (K, 4)+P (K, 3));
         XX (K) := TX;
         FX (K) := (P (K, 1), P (K, 2))/TX;
         --FX (K) := CMPLX (P (K, 1), P (K, 2))/TX;
      end loop;

      for J in T_INDICE_PARTICULE'First..T_INDICE_PARTICULE'Pred (T_INDICE_PARTICULE'Last) loop -- 1..4
         for K in J+1..T_INDICE_PARTICULE'Last loop -- 1..5
          -- Produit spinoriel de M.Mangano, S.Parke
            CX := FX (J)*XX (K)-FX (K)*XX (J);
            OSpinor.S (J, K) :=  CX;
            OSpinor.S (K, J) := -CX;
            OSpinor.T (K, J) :=  Conjugate (CX);
            OSpinor.T (J, K) := -Conjugate (CX);
         end loop;
      end loop;
      return OSpinor;
   end OSpinor;

   function APPM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex is
   begin
      return -RAC8 * OSpinor.S (1, 2) * OSpinor.S (1, K3)**2 / OSpinor.S (1, K1)/OSpinor.S (1, K2)/OSpinor.S (2, K1) / OSpinor.S (2, K2);
   end APPM;

   function APMM (OSpinor : in Spinor; k1, k2, K3 : in Integer) return Complex is
   begin
      return -RAC8 * OSpinor.T (1, 2) * OSpinor.T (2, K1)**2 / OSpinor.T (2, K2) / OSpinor.T (2, K3) / OSpinor.T (1, K2) / OSpinor.T (1, K3);
   end APMM;

   function BPPM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex is
   begin
      return -RAC8*OSpinor.T (1, 2)* (OSpinor.T (K1, K2)*OSpinor.S (K3, 1))**2;
   end BPPM;

   function BPMM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex is
   begin
      return -RAC8*OSpinor.S (1, 2)* (OSpinor.T (K1, 2)*OSpinor.S (K2, K3))**2;
   end BPMM;

   function BPPP (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex is
   begin
      return -RAC8*OSpinor.S (1, 2)* (
         + (OSpinor.T (K1, K2)*OSpinor.T (K3, 2))**2
         + (OSpinor.T (K1, K3)*OSpinor.T (K2, 2))**2
         + (OSpinor.T (K2, K3)*OSpinor.T (K1, 2))**2
         );
   end BPPP;

   function BMMM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex is
   begin
      return -RAC8*OSpinor.T (1, 2)* (
         + (OSpinor.S (K1, 1)*OSpinor.S (K2, K3))**2
         + (OSpinor.S (K2, 1)*OSpinor.S (K1, K3))**2
         + (OSpinor.S (K3, 1)*OSpinor.S (K1, K2))**2
         );
   end BMMM;

end Spinor;
