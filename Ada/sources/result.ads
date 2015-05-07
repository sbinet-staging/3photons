--   Common pour passer les resultats engendrees par mat(rice)
with Spinor;
with Param;
with Precision;
package Result is
   use Precision;
   -- duplication de code entre resfin et result
   NRESUL : constant Integer := 8;
   subtype Result_Category is POSITIVE range 1..NRESUL;
   type Result_Arr is array (1..2, 1..2, 1..2, Result_Category) of Real;
   type Result is record
      M2 : Result_Arr;
   end record;
   function oResult (OSpinor : in Spinor.Spinor; OParam : in Param.Param; ETOT : Real) return Result;
end Result;
