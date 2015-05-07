--   Common pour passer les angles
with Ppp;
with Scalar;
with Cutpar;
with Precision;

package Angle is
  use Precision;
  use Ppp;
  use Scalar;
  use Cutpar;
  type Angle is record
     Cos12, Cos1p1, Cos2p1, Cos3p1, Cos13, Cos23, Cosn, Cosac : Real;
  end record;

  -- Constructeur
  function oAngle (OPpp : in Ppp.Ppp; OScalar : in Scalar.Scalar) return Angle;
  -- Methodes
  function CUT (OPpp : Ppp.Ppp; OCutpar : Cutpar.Cutpar; OAngle : Angle) return Boolean;
end Angle;
