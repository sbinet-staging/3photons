--   Common pour passer les résultats finaux
with Precision;
with Param;
with Cutpar;
use  Precision;
use  Param;
use  Cutpar;

package Resfin is
   use Specific_Complex;
   -- duplication de code entre resfin et result
   NRES : constant := 8;
   subtype T_INDICE_RES is POSITIVE range 1..NRES;
   type Tmp_Resultat_Array is array (T_INDICE_RES) of Real;
   type Resultat_Array is array (1..2, T_INDICE_RES) of Real;

   type Resfin is record
      Spm2dif : Tmp_Resultat_Array;
      Spm2    : Resultat_array;
      Var     : Resultat_array;
   end record;

   procedure ERIC  (UNE : in Integer; BREPEM, CONVERS, PI : in Real;
                                      OParam : in Param.Param;
                                      OResfin : in Resfin);

   procedure FAWZI (UNE : in Integer; BREPEM, CONVERS, PI, ETOT : in Real;
                                      OParam  : in Param.Param;
                                      OCutpar : in Cutpar.Cutpar;
                                      OResfin : in Resfin);
end Resfin;
