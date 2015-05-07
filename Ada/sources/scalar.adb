with Spinor;
with Precision;

package body Scalar is

   function oScalar (OSpinor : Spinor.Spinor) return Scalar is
      use Precision;
      use Specific_Complex;
      OScalar : Scalar;
   begin
      for K in T_INDICE_PARTICULE loop
         OScalar.PS (K, K) := 0.0;
      end loop;
      for J in T_INDICE_PARTICULE'First..T_INDICE_PARTICULE'Pred (T_INDICE_PARTICULE'Last) loop
         for K in J+1..T_INDICE_PARTICULE'Last loop
            -- produit scalaire de Lorentz
            --OScalar.PS (J, K) := OSpinor.S (J, K) * CONJG (OSpinor.S (J, K)) / 2.0;
            OScalar.PS (J, K) := (abs OSpinor.S (J, K))**2 / 2.0;
            OScalar.PS (K, J) := OScalar.PS (J, K);
         end loop;
      end loop;
      return OScalar;
   end OScalar;

end Scalar;
