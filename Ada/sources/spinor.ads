--   Common pour passer les produits spinoriels engendres par mat3
with Precision;
with Ppp;
package Spinor is
   use Precision;
   use Specific_Complex;
   use Specific_Elementary;
   subtype T_INDICE_PARTICULE is POSITIVE range 1..5;
   subtype T_INDICE_MINK is POSITIVE range 1..4;
   subtype T_INDICE_PARTICULE_SORTANTE is POSITIVE range 1..3;
   RAC8 : constant Real := Sqrt (Real (8.0));
   --   Le facteur RAC8 vient des Sqrt (2) de chaque vecteur de polarisation
   --   dans la methode des amplitudes d'helicite

   type Spinor_array is array (T_INDICE_PARTICULE, T_INDICE_PARTICULE) of Complex;
   type Spinor is record
      s, t :  Spinor_array;
   end record;

   -- Constructeur
   function oSpinor (OPpp : in Ppp.Ppp) return Spinor;
   function APPM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex;
   function APMM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex;
   function BPPM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex;
   function BPMM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex;
   function BPPP (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex;
   function BMMM (OSpinor : in Spinor; K1, K2, K3 : in Integer) return Complex;
end Spinor;
