--   Common pour passer les produits scalaires engendres par mat3
with Spinor;
with Precision;

package Scalar is
   use Precision;
   use Spinor;
   type Scalar_Array is array (T_INDICE_PARTICULE, T_INDICE_PARTICULE) of Real;
   type Scalar is record
      Ps :  Scalar_Array;
   end record;

   -- Constructeur
   function oScalar (OSpinor : Spinor.Spinor) return Scalar;
end Scalar;
