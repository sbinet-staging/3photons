with Precision;
use  Precision;
with Hasa;
use  Hasa;

package body Ppp is
   use Specific_Elementary;

   IBEGIN : Integer := 0; --integer, save      ::
   NM : Integer := 0; --integer, save      ::
   TWOPI : Real := 8.0 * Arctan (1.0); --real (pr),save ::
   PO2LOG : Real := Log (TWOPI / 4.0); --real (pr),save ::
   XMT : Real := 0.0; --real (pr),save ::
   Z : T_VECTEUR_PARTICULE; --real (pr),save ::
   -- NORM := (TWOPI**(4-3*INP))/REAL(ITOT)

   function oPpp (ETOT : in Real) return Ppp is
      -- Initialisation des impulsions entrantes
      OPpp : Ppp;
   begin
      OPpp.P1 (4) :=  ETOT/2.0;
      OPpp.P2 (4) :=  ETOT/2.0;
      OPpp.P1 (1) := -ETOT/2.0;
      OPpp.P2 (1) :=  ETOT/2.0;
      OPpp.P1 (2) := 0.0;
      OPpp.P1 (3) := 0.0;
      OPpp.P2 (2) := 0.0;
      OPpp.P2 (3) := 0.0;
      return (OPpp);
   end OPpp;

   procedure TRI (OPpp : in out Ppp) is
      PE : T_VECTEUR_MINKOWSKI;
   begin
      -- Trie les photons sortants par energie
      if (OPpp.POUT (4, 2) > OPpp.POUT (4, 1)) then
         for I in T_INDICE_MINKOWSKI loop
            PE (I)           := OPpp.POUT (I, 1);
            OPpp.POUT (I, 1) := OPpp.POUT (I, 2);
            OPpp.POUT (I, 2) := PE (I);
         end loop;
      end if;
      if (OPpp.POUT (4, 3) > oPpp.POUT (4, 1)) then
         for I in T_INDICE_MINKOWSKI loop
            PE (I)           := OPpp.POUT (I, 1);
            OPpp.POUT (I, 1) := OPpp.POUT (I, 3);
            OPpp.POUT (I, 3) := PE (I);
         end loop;
      end if;
      if (OPpp.POUT (4, 3) > OPpp.POUT (4, 2)) then
         for I in T_INDICE_MINKOWSKI loop
            PE (I)           := OPpp.POUT (I, 2);
            OPpp.POUT (I, 2) := OPpp.POUT (I, 3);
            OPpp.POUT (I, 3) := PE (I);
         end loop;
      end if;
   end TRI;

   procedure RAMBO (N : in Integer; ET : in Real; XM : in T_VECTEUR_PARTICULE; OPpp : in out Ppp; WT : out Real) is
      ----------------------------------------------------------------------------
      --
      --                       RAMBO
      --
      --    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
      --
      -- ****** VERSION LIGHT & LIGHTSPEED ******
      --
      --    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
      --    AUTHORS@D  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
      --    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
      ----------------------------------------------------------------------------
      -- integer, intent (in)    :: N  -- Number of particles (>1, IN THIS VERSION <101)
      -- real (pr), intent (in)  :: ET -- Total centre-of-mass energy
      -- real (pr), intent (in)  :: XM (NP)   -- Particle masses (DIM=100)
      -- type (ppp_typ), intent (inout) :: oPpp -- Particle momenta (DIM=(4, 100))
      -- real (pr), intent (out) :: WT -- Weight of the event

      C, S, F, RMAS, G, A, X, BQ : Real;
      Q : T_EVENT;
      R : T_VECTEUR_MINKOWSKI;
      B : array (1..3) of Real;
      IWARN : array (1..5) of Integer := (others => 0); --data IWARN/5*0/
      Wrong_Event : Exception;
   begin
      -- Initialization step@d factorials for the phase space weight
      if (IBEGIN = 0) then
         IBEGIN := 1;
         TWOPI  := 8.0 * Arctan (1.0);
         PO2LOG := Log (TWOPI / 4.0);
         Z (2) := PO2LOG;
         for K in 3..NP loop
            Z (K) := Z (K-1) + PO2LOG -2.0 * Log (Real (K-2));
         end loop;
         for K in 3..NP loop
            Z (K) := (Z (K) - Log (Real (K-1)));
         end loop;

-- Check on the number of particles
         if (N > 1 and N < 101) then
            goto L104;
         end if;
         Put (" RAMBO FAILS@D # OF PARTICLES ="); Put (N, 5); Put (" IS NOT ALLOWED");--print 1001, N
         Put_Line ("");
         -- 1001 format (' RAMBO FAILS@D # OF PARTICLES =', I5, ' IS NOT ALLOWED')
         raise Wrong_Event; --Stop;

-- Check whether total energy is sufficient; Count nonzero masses
         <<L104>>    XMT := 0.0;
             NM := 0;
             for I in 1..N loop
                if (XM(I) /= 0.0) then
                   NM := NM + 1;
                end if;
                XMT := XMT + abs (XM (I));
             end loop;
             if (XMT <= ET) then
                goto L201;
             end if;
             Put (" RAMBO FAILS@D TOTAL MASS =");         Put (Real'Image (XMT)); --1002,
             Put (" IS NOT SMALLER THAN TOTAL ENERGY ="); Put (Real'Image (ET)); --1002,
             Put_Line ("");
             -- 1002 FORMAT (' RAMBO FAILS@D TOTAL MASS =', D15.6, ' IS NOT', &
             --              ' SMALLER THAN TOTAL ENERGY =', D15.6)
             raise Wrong_Event; -- Stop;
      end if;

-- The parameter values are now accepted

-- Generate N massless momenta in infinite phase space
<<L201>> for I in 1..N loop
       C := 2.0* Rn -1.0; -- C := 2.0*RN (1)-1.0;
       S := Sqrt (1.0-C*C);
       F := TWOPI * Rn; -- TWOPI * RN (2);
       Q (4, I) := -Log (Rn * Rn); --  -Log (RN (3) * RN (4));
       Q (3, I) := Q (4, I) * C;
       Q (2, I) := Q (4, I) * S * Cos (F);
       Q (1, I) := Q (4, I) * S * Sin (F);
    end loop;

-- Calculate the parameters of the conformal transformation
    for I in T_INDICE_MINKOWSKI loop
       R(I) := 0.0;
    end loop;
    for I in 1..N loop
       for K in T_INDICE_MINKOWSKI loop
          R (K) := R (K) + Q (K, I);
       end loop;
    end loop;
    RMAS := Sqrt (R (4)**2 - R (3)**2 - R (2)**2 - R (1)**2);
    for K in 1..3 loop
       B (K) := -R (K) / RMAS;
    end loop;
    G := R (4) / RMAS;
    A := 1.0 / (1.0+G);
    X := ET / RMAS;

-- Transform the Q's conformally into the P's
    for I in 1..N loop
       BQ := B (1) * Q (1, I) + B (2) * Q (2, I) + B (3) * Q (3, I);
       for K in 1..3 loop
          oPpp.pout (K, I) := X * (Q (K, I) + B (K) * (Q (4, I) + A * BQ));
       end loop;
       oPpp.pout (4, I) := X * (G * Q (4, I) + BQ);
    end loop;

-- CALCULE LE POIDS ET LES AVERTISSEMENTS EVENTUELS
    WT := PO2LOG;
    if (N /= 2) then
       WT := Real (2 * N - 4) * Log (ET) + Z (N);
    end if;
    if (WT >= -180.0) then
       goto L208;
    end if;
    if (IWARN (1) <= 5) then
       Put (" RAMBO WARNS@D WEIGHT = EXP(");
       Put (Real'Image (WT)); -- 1004,
       Put (") MAY UNDERFLOW");
       Put_Line ("");
       -- 1004 format (' RAMBO WARNS@D WEIGHT = EXP(', F20.9, ') MAY UNDERFLOW')
    end if;
    IWARN (1) := IWARN (1) + 1;
    <<L208>> if (WT <= 174.0) then
      goto L209;
    end if;
    if (IWARN (2) <= 5) then
       Put (" RAMBO WARNS@D WEIGHT = EXP(");
       Put (Real'Image (WT)); -- 1005,
       Put (") MAY  OVERFLOW");
       Put_Line ("");
       -- 1005 format (' RAMBO WARNS@D WEIGHT = EXP(', F20.9, ') MAY  OVERFLOW')
    end if;
    IWARN (2) := IWARN (2) + 1;

-- Renvoie les impulsions sans masses ponderees
<<L209>> WT := exp (WT);

  end RAMBO;
end Ppp;
