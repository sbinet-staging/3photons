--   Common pour passer les probabilites supplementaires

with Precision;

package Xtrpro is
  use Precision;
   type Xtrpro is record
     Prplus, Prmoins : Real; -- real*4 ::
     Ee1, Ee2 : Real; -- real (pr) ::
  end record;
end Xtrpro;
-- present dans mc (naturellement)
-- present dans book
