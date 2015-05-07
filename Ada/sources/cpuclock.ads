package CPUClock is

   -- specification for a package to do CPU timing of algorithms

   subtype CPUSecond is Float range 0.0..Float'Last;
   -- We make CPUSecond a Float type so the usual operations are available

   procedure ResetCPUTime;
   -- Pre:  none
   -- Post: resets a CPU timer

   function CPUTime return CPUSecond;
   -- Pre:  none
   -- Post: returns the number of CPUSeconds since the last reset

   function CPUSysTime return CPUSecond;
   -- Pre:  none
   -- Post: returns the number of System CPUSeconds since the last reset

end CPUClock;
