
*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_init              *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_init(rtdb0,geom)
      implicit none
      integer rtdb0,geom

#include "errquit.fh"
#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "btdb.fh"

*     **** parallel variables ****
      logical  omaster,oprint
      integer  taskid
      integer  MASTER
      parameter(MASTER=0)

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,nefc,rtdb
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

*     ***** external functions *****
      logical  control_print
      external control_print
      integer  Parallel_threadid
      external Parallel_threadid

      rtdb = rtdb0
c
c     == Cosmo initialization ==
      call Parallel_taskid(taskid)
      omaster =  (taskid.eq.MASTER).and.(Parallel_threadid().eq.MASTER)
      oprint = omaster.and.control_print(print_medium)

      do_cosmo_smd = .false.
      nefc = 0
      gcds = 0.0d0

      if (.not.btdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) 
     >   cosmo_on = .false.

      if(cosmo_on) then
         call cosmo_initialize(rtdb,geom,-1,oprint)
c
c        Turn cosmo on, we want to run the calculation
c        Start with gas_phase run unless told otherwise
c
         cosmo_on = .true.
         if (.not.btdb_get(rtdb,'cosmo_phase',mt_int,1,cosmo_phase))
     >      cosmo_phase = 1

         if (.not.btdb_get(rtdb,'cosmo:do_cosmo_smd',
     >                     mt_log,1,do_cosmo_smd))
     >      do_cosmo_smd = .false.
     
         gcds = 0.0d0
         if (do_cosmo_smd) then
            if (.not.btdb_get(rtdb,'cosmo:gcds',mt_dbl,1,gcds)) gcds = 0.0d0
         end if

         if(.not.btdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc)) nefc = 0

      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_end               *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_end()
      implicit none
      
*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

c     == Clean up cosmo stuff, if necessary ==
      if (cosmo_on) then
         call cosmo_tidy(rtdb)
         cosmo_on = .false.
         cosmo_phase = 1
         nefc = 0
         gcds = 0.0d0
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_nq                *
*     *                                         *
*     *******************************************
      integer function nwpw_cosmo_nq()
      implicit none

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      nwpw_cosmo_nq = 0
      if (cosmo_on) then
         nwpw_cosmo_nq = nefc
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_phase             *
*     *                                         *
*     *******************************************
      integer function nwpw_cosmo_phase()
      implicit none

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      nwpw_cosmo_phase = 0
      if (cosmo_on) then
         nwpw_cosmo_phase = cosmo_phase
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_smd               *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo_smd()
      implicit none

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      nwpw_cosmo_smd = .false.
      if (cosmo_on) then
         nwpw_cosmo_smd = do_cosmo_smd
      end if
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_on                *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo_on()
      implicit none

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      nwpw_cosmo_on = .false.
      if (cosmo_on) then
         nwpw_cosmo_on = cosmo_on
      end if
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_qcoords           *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_qcoords(qcoord)
      implicit none
      real*8 qcoord(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      if (cosmo_on) then
         if(.not.btdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,qcoord))
     >     call errquit('nwpw_cosmo_qcoord: failed',0,rtdb_err)
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_qs                *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_qs(qs)
      implicit none
      real*8 qs(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      if (cosmo_on) then
         if(.not.btdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,qs))
     >     call errquit('nwpw_cosmo_qs: failed',0,rtdb_err)
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_qsurfaces         *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_qsurfaces(qs)
      implicit none
      real*8 qs(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      if (cosmo_on) then
         if(.not.btdb_get(rtdb,'cosmo:efcs',mt_dbl,nefc,qs))
     >     call errquit('nwpw_cosmo_qsurfaces: failed',0,rtdb_err)
      endif
      return
      end



*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_qa_indexing       *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_qa_indexing(indx)
      implicit none
      integer indx

#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  gcds
      common /nwpw_cosmo_blk/ gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      if (cosmo_on) then
         if(.not.btdb_get(rtdb,'cosmo:efciat',mt_int,nefc,indx))
     >     call errquit('nwpw_cosmo_qa_indexing: failed',0,rtdb_err)
      endif
      return
      end

