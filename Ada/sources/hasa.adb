with Precision;
use  Precision;
package body Hasa is
   -- IA : IntRanVec; --, save ::
   -- NCALL : Integer := 0; --, save ::
   -- MCALL : Integer := 55; --, save ::

   procedure KRanPool is
   begin
      for I in 1..PoolSize loop
         for J in 1..4 loop
            Ranpool (I, J) := Rn; -- Ranpool (I, J) := rn (0);
         end loop;
      end loop;
   end KRanPool;

   procedure GRanVec (Xvec : out RanVec) is
   begin
      if (Icounter = 0) then
         KRanPool;
         Icounter := 1;
      end if;
      for J in RanVec'Range loop
         Xvec (J) := Ranpool (Icounter, J);
      end loop;
      Icounter := Icounter + 1;
      if (Icounter > PoolSize) then
         Icounter := 0;
      end if;
   end GRanVec;

   -- on va construire un tableau de vecteur aleatoires
   -- tous les ... on le remplit a bloc
   -- l'accesseur nous renvoie un tableau de quadrivecteurs aleatoires

   ----------------------------------------------------------------------------
   function Rn return Real is
      ----------------------------------------------------------------------------
      --    ***************************************************************
      --    * RANDOM NUMBER FUNCTION TAKEN FROM KNUTH RANF                *
      --    * (SEMINUMERICAL ALGORITHMS).                                 *
      --    * METHOD IS X (N)=MOD (X (N-55)-X (N-24), 1/FMODUL)           *
      --    * NO PROVISION YET FOR CONTROL OVER THE SEED NUMBER.          *
      --    *                                                             *
      --    * RANF GIVES ONE RANDOM NUMBER BETWEEN 0 AND 1.               *
      --    * IRN55 GENERATES 55 RANDOM NUMBERS BETWEEN 0 AND 1/FMODUL.   *
      --    * IN55  INITIALIZES THE 55 NUMBERS AND WARMS UP THE SEQUENCE. *
      --    ***************************************************************
      --IA (1..55), NCALL, MCALL : Integer; --, save ::
      XRN : Real;
      FMODUL : constant Real := 1.0e-9;
   begin
      if (NCALL = 0) then
         IN55 (IA, 234612947);
         NCALL := 1;
      end if;
      if (MCALL = 0) then
         IRN55 (IA);
         MCALL := 55;
      end if;
      XRN := Real (IA (MCALL)) * FMODUL;
      MCALL := MCALL - 1;
      return (XRN);
   end Rn;

  ----------------------------------------------------------------------------
   procedure IN55 (IA : in out IntRanVec; IX : Integer) is
      MODULO : constant Integer := 1_000_000_000;
      J, K, II : Integer; -- I is a dummy loop variable
   begin
      IA (55) := IX;
      J := IX;
      K := 1;
      for I in 1..54 loop
         II := (21*I) mod 55;
         IA (II) := K;
         K := J - K;
         if (K < 0) then
            K := K + MODULO;
         end if;
         J := IA (II);
      end loop;
      for I in 1..10 loop
         IRN55 (IA);
      end loop;
   end IN55;

   ----------------------------------------------------------------------------
   procedure IRN55 (IA : in out IntRanVec) is
      MODULO : constant Integer := 1_000_000_000;
      J : Integer;
   begin
      for I in 1..24 loop
         J := IA (I) - IA (I+31);
         if (J < 0) then
            J := J + MODULO;
         end if;
         IA (I) := J;
      end loop;
      for I in 25..55 loop
         J := IA (I) - IA (I-24);
         if (J < 0) then
            J := J + MODULO;
         end if;
         IA (I) := J;
      end loop;
   end IRN55;
   ----------------------------------------------------------------------------
end Hasa;
----------------------------------------------------------------------------
