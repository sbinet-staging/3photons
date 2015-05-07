--   Common pour passer les coupures
with Precision;
package Cutpar is
   use Precision;
   type Cutpar is record
      Acut, Bcut, Emin, Sincut : Real;
   end record;
end Cutpar;
