--   Common pour passer les impulsions engendrees par Rambo
with Precision;

with Ada.Text_IO, Ada.Integer_Text_IO, Ada.Float_Text_IO;
use  Ada.Text_IO, Ada.Integer_Text_IO, Ada.Float_Text_IO;
package Ppp is
   use Precision;
   INP : constant Integer := 3; -- Number of impulsions generated
   NP  : constant Integer := 100; -- will be changed to dynamic sizing later ;
                                  -- kept here for comparison purposes

   type T_VECTEUR is array (POSITIVE range <>) of Real;
   subtype T_INDICE_MINKOWSKI  is POSITIVE range 1..4;
   subtype T_INDICE_PARTICULE  is POSITIVE range 1..NP;
   type T_VECTEUR_MINKOWSKI is new T_VECTEUR (T_INDICE_MINKOWSKI);
   type T_VECTEUR_PARTICULE is new T_VECTEUR (T_INDICE_PARTICULE);

   type T_EVENT is array (T_INDICE_MINKOWSKI, T_INDICE_PARTICULE) of Real;

   type Ppp is record
      Pout : T_EVENT;
      P1 : T_VECTEUR_MINKOWSKI;
      P2 : T_VECTEUR_MINKOWSKI;
   end record;

   -- Accesseurs
   -- Constructeurs
   function oPpp (ETOT : in Real) return Ppp;

   -- Méthodes
   procedure TRI (OPpp : in out Ppp);

   --procedure RAMBO (N, ET, XM, oPpp, WT)
   procedure RAMBO (N : in Integer; ET : in Real; XM : in T_VECTEUR_PARTICULE; OPpp : in out Ppp; WT : out Real);

end Ppp;
