c     $Id$
      logical function ofpw_energy(rtdb)
      implicit none
      integer rtdb
      logical  ofpw_minimizer
      external ofpw_minimizer

      write(*,*) "hello from ofpw_energy"
      ofpw_energy = ofpw_minimizer(rtdb,0)
      return
      end

