with Text_IO;

package body CPUClock is

   function unixtime return Integer;
   pragma Interface (C, unixtime);

   function Usertime return Integer;
   pragma Interface (C, usertime);

   function Systtime return Integer;
   pragma Interface (C, systtime);

   function sysconf_sc_clk_tck return Integer;
   pragma Interface (C, sysconf_sc_clk_tck);

   Saved_Time: Integer;
   Saved_System_Time: Integer;

   function CPUTime return CPUSecond is
   begin
      return CPUSecond (unixtime - Saved_Time) / CPUSecond (sysconf_sc_clk_tck);
   end CPUTime;

   function CPUSysTime return CPUSecond is
   begin
      return CPUSecond (Systtime - Saved_Time) / CPUSecond (sysconf_sc_clk_tck);
      --return CPUSecond (Systtime - Saved_System_Time) / CPUSecond (sysconf_sc_clk_tck);
   end CPUSysTime;

   procedure ResetCPUTime is
   begin
      Saved_Time := Unixtime;
      Saved_System_Time := Systtime;
   end ResetCPUTime;

begin -- initialization of package
   ResetCPUTime;
end CPUClock;
