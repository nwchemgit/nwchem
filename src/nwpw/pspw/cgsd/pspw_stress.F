      logical function pspw_stress(rtdb)
      implicit none
      integer rtdb
      logical  cgsdv5
      external cgsdv5

      pspw_stress = cgsdv5(rtdb,3)
      return
      end

