with Text_IO;
with Precision;

package body Result is
   use Precision;
   use Specific_Complex;

   ------------------------------------------------------------------------
   --
   -- Calcule l'élément de matrice e+e- en 3 photons
   --
   ------------------------------------------------------------------------
   function oResult (OSpinor : in Spinor.Spinor; OParam : in Param.Param; ETOT : Real) return Result is
      OResult     : Result;
      L1, L2, L3  : Integer;
      A, BP, BM   : Complex;
      M2          : Result_Arr;
      use Text_IO;
      use Spinor;
   begin
      -- LES - VIENNENT DES Photons SORTANTS
      for LL1 in 1..2 loop
         L1 := - (2*LL1-3);
         for LL2 in 1..2 loop
            L2 := - (2*LL2-3);
            for LL3 in 1..2 loop
               L3 := - (2*LL3-3);
               -- CALCULS DES AMPLITUDES D'HELICITE
               if  ((l1=1) and then (l2=1) and then (l3=1)) then
                  A  := (0.0, 0.0);
                  BP := (0.0, 0.0);
                  BM := BPPP (oSpinor, 3, 4, 5);
               end if;
               if  ((l1=-1) and then (l2=-1) and then (l3=-1)) then
                  A  := (0.0, 0.0);
                  BP := (0.0, 0.0);
                  BM := BMMM (oSpinor, 3, 4, 5);
               end if;
               if  ((l1=1) and then (l2=1) and then (l3=-1)) then
                  A  := APPM (oSpinor, 3, 4, 5);
                  BP := BPPM (oSpinor, 3, 4, 5);
                  BM := (0.0, 0.0);
               end if;
               if  ((l1=1) and then (l2=-1) and then (l3=1)) THEN
                  A  := APPM (oSpinor, 5, 3, 4);
                  BP := BPPM (oSpinor, 5, 3, 4);
                  BM := (0.0, 0.0);
               end if;
               if  ((l1=-1) and then (l2=1) and then (l3=1)) THEN
                  A  := APPM (oSpinor, 4, 5, 3);
                  BP := BPPM (oSpinor, 4, 5, 3);
                  BM := (0.0, 0.0);
               end if;
               if  ((l1=1) and then (l2=-1) and then (l3=-1)) THEN
                  A  := APMM (oSpinor, 3, 4, 5);
                  BP := BPMM (oSpinor, 3, 4, 5);
                  BM := (0.0, 0.0);
               end if;
               if  ((l1=-1) and then (l2=-1) and then (l3=1)) THEN
                  A  := APMM (oSpinor, 5, 3, 4);
                  BP := BPMM (oSpinor, 5, 3, 4);
                  BM := (0.0, 0.0);
               end if;
               if  ((l1=-1) and then (l2=1) and then (l3=-1)) THEN
                  A  := APMM (oSpinor, 4, 3, 5);
                  BP := BPMM (oSpinor, 4, 3, 5);
                  BM := (0.0, 0.0);
               end if;
-- COUPLAGES
               A  := A  * oParam.GA;
               BP := BP * oParam.GBP;
               BM := BM * oParam.GBM;
-- LES DIFFERENTS TERMES DE L'ELEMENT DE MATRICE AU CARRE
               M2 (LL1, LL2, LL3, 1) := (abs A )**2;
               M2 (LL1, LL2, LL3, 2) := (abs BP)**2;
               M2 (LL1, LL2, LL3, 3) := (abs BM)**2;
               M2 (LL1, LL2, LL3, 4) := 2.0 * Re (A * Conjugate (BP));
               M2 (LL1, LL2, LL3, 5) := 2.0 * Im (A * Conjugate (BP));
               M2 (LL1, LL2, LL3, 6) := 0.0;
               M2 (LL1, LL2, LL3, 7) := 0.0;
               M2 (LL1, LL2, LL3, 8) := 0.0;

             --oResult.M2 (LL1, LL2, LL3, 1) := (abs A )**2;
             --oResult.M2 (LL1, LL2, LL3, 2) := (abs BP)**2;
             --oResult.M2 (LL1, LL2, LL3, 3) := (abs BM)**2;
             --oResult.M2 (LL1, LL2, LL3, 4) := 2.0 * Re (A * Conjugate (BP));
             --oResult.M2 (LL1, LL2, LL3, 5) := 2.0 * Im (A * Conjugate (BP));
             --oResult.M2 (LL1, LL2, LL3, 6) := 0.0;
             --oResult.M2 (LL1, LL2, LL3, 7) := 0.0;
             --oResult.M2 (LL1, LL2, LL3, 8) := 0.0;

--           M2 (LL1, LL2, LL3, 1) = CDABS (A )**2
--           M2 (LL1, LL2, LL3, 2) = CDABS (BP)**2
--           M2 (LL1, LL2, LL3, 3) = CDABS (BM)**2
--           M2 (LL1, LL2, LL3, 4) = 2.*DREAL (A *DCONJG (BP))
--           M2 (LL1, LL2, LL3, 5) = 2.*DIMAG (A *DCONJG (BP))

            end loop;
         end loop;
      end loop;

      if oParam.IMPR then
         for K in 1..NRESUL loop
            Put_Line (Integer'Image (K));
            Put_Line ("  ---     --+     -+-     +--     -++     +-+     ++-     +++");
            --write (*, 900)
            --     oResult.M2 (1, 1, 1, K), oResult.M2 (1, 1, 2, K),
            --     oResult.M2 (1, 2, 1, K), oResult.M2 (2, 1, 1, K),
            --     oResult.M2 (1, 2, 2, K), oResult.M2 (2, 1, 2, K),
            --     oResult.M2 (2, 2, 1, K), oResult.M2 (2, 2, 2, K)
         end loop;
      end if;

--900 format (8E8.2)
      OResult.M2 := M2;
      return OResult;
   end OResult;

end Result;
