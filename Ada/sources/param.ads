--   Common pour passer les parametres physiques

with Precision;

package Param is
   use Precision;
   type Param is record
     Mz0, Gz0, Ga, Gbp, Gbm, Polp, Polm, Polp2, Polm2 : Real;
     Impr : Boolean;
  end record;
end Param;
