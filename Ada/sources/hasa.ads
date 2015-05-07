with Precision;

package Hasa is
   use Precision;
   PoolSize : constant Integer := 55 * 3;
   Ranpool : array (1..PoolSize, 1..4) of Real;
   type RanVec is array (1..4) of Real;
   type IntRanVec is array (1..55) of Integer;

   IA : IntRanVec; --, save ::
   NCALL : Integer := 0; --, save ::
   MCALL : Integer := 55; --, save ::

   Icounter : Integer := 0; -- , save

   -- Constructeurs
   procedure KRanPool; --private :: kRanPool

   -- Accesseurs
   procedure GRanVec (Xvec : out RanVec);

   function Rn return Real;

   procedure IN55 (IA : in out IntRanVec; IX : Integer);
   procedure IRN55 (IA : in out IntRanVec);

end Hasa;
