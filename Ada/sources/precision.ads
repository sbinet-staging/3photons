with Ada.Numerics.Generic_Complex_Types;
with Ada.Numerics.Generic_Elementary_Functions;

package Precision is
   --type Real is new Float;
   --type Real is new Long_Long_Float;
   type Real is new Long_Float;
   package Specific_Complex    is new Ada.Numerics.Generic_Complex_Types (Real => Real);
   package Specific_Elementary is new Ada.Numerics.Generic_Elementary_Functions (Float_Type => Real);
end Precision;
