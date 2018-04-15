
*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_init              *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_init(rtdb,geom)
      implicit none
      integer rtdb,geom

#include "errquit.fh"
#include "bafdecls.fh"
#include "stdio.fh"
#include "btdb.fh"

      logical oprint

      logical cosmo_on,cosmo_last
      integer cosmo_phase
      common /nwpw_cosmo_blk/ cosmo_phase,cosmo_on,cosmo_last
c
c     == Cosmo initialization ==
      oprint = .true.
      cosmo_last = .false.
      if (.not.btdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) 
     >   cosmo_on = .false.

      if(cosmo_on) then
         call cosmo_initialize(rtdb,geom,-1,oprint)
c
c        Turn cosmo on, we want to run the calculation
c        Start with gas_phase run unless told otherwise
c
         cosmo_last = .true.
         cosmo_on = .true.
         if(.not.rtdb_get(rtdb,'cosmo_phase',mt_int,1,cosmo_phase))
     >      cosmo_phase = 1

      endif
      return
      end

      subroutine nwpw_cosmo_end(rtdb)
      implicit none
      integer rtdb

      logical cosmo_on,cosmo_last
      integer cosmo_phase
      common /nwpw_cosmo_blk/ cosmo_phase,cosmo_on,cosmo_last

c     == Clean up cosmo stuff, if necessary ==
      if (cosmo_on) then
         call cosmo_tidy(rtdb)
         cosmo_on = .false.
         cosmo_phase = 1
      endif
      return
      end
