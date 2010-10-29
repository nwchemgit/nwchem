      subroutine grad_store ( rtdb, geom, g_force, chktime, nat, nquad, 
     $     frc, oinit, ofinish, oactive, oskel, omp2 )

c$Id$
C=======================================================================
C     Syncs and writes the partially calculated forces to rtdb if either:
C     - checkpointing is requested and the specified time has elapsed. 
C     - ofinish == .true.: gradient calculation is finished
C     oinit == .true.: just initializes timer
C=======================================================================
      implicit none
#include "errquit.fh"

#include "util.fh"
#include "msgids.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      integer rtdb, geom, g_force, chktime, nat, nquad

      double precision frc
      dimension frc (3,nat)

      logical oinit, ! initialize time and return
     $     ofinish,  ! this is the final storage of forces
     $     oactive,  ! array of active atoms (or point charges)
     $     oskel,    ! use symmetry
     $     omp2      ! this is the separable part for mp2
      dimension oactive (nat)

      integer lnquad
      double precision last, now, interval
      logical status

      character*26 rtdb_string
      data rtdb_string /'gradients:partial gradient'/

      save last

      if ( oinit ) then
        last = util_wallsec()
        return
      endif

      now = util_wallsec()

      interval = chktime
      if ( ofinish .or. ((now-last) .gt. interval) ) then
C       CREATE A CHECKPOINT

        rtdb_string = 'gradients:partial gradient'

        call ga_dgop ( msg_grad_2el, frc, 3*nat, '+' )
        if ( ga_nodeid() .eq. 0 ) then
          call ga_acc ( g_force, 1, 3, 1, nat, frc, 3, 1.D0 )
        endif
        call ga_sync()
        if ( ga_nodeid() .eq. 0 ) then
          call ga_get ( g_force, 1, 3, 1, nat, frc, 3 )
        endif

C       find out highest atom quadrupel
        lnquad = nquad
        call ga_igop ( Msg_grad_nquad, lnquad, 1, 'max' )

C       find out if some other process has finished the loops
        call grad_finished ( ofinish )
C       process zero only gets the accumulated forces
        status = rtdb_parallel (.false.)
        if (ga_nodeid().eq.0) then

          call ga_get ( g_force, 1, 3, 1, nat, frc, 3 )

          if ( util_print('information', print_high) ) then
            call util_print_centered ( luout,
     $           'writing forces to rtdb', 35, .false.)
            call util_flush(luout)
          endif

C         store forces and # of atom quadruples in the database
          status = rtdb_put ( rtdb, rtdb_string, MT_DBL, 3*nat, frc)
          if ( .not. status ) 
     $         call errquit('could not store gradients',1, RTDB_ERR)
          status = rtdb_put ( rtdb, 'gradients:ndone', MT_INT, 1, 
     $         lnquad )
        endif
C       zero the local 2-el. force array
        call dfill ( 3*nat, 0.D0, frc, 1 )
        status = rtdb_parallel (.true.)

        call ga_sync()

        last = now
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function grad_restart ( rtdb, geom, basis, nat, ndone, 
     $     frc, g_force )

C=======================================================================
C     Finds out if this is a restart of a partially finished gradient run.
C     If yes, returns the previously calculated forces (proc 0 only)
C     and the last calculated atom quadruple in grad2
C========================================================================
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"

      integer rtdb, geom, basis, nat,
     $     ndone, g_force

      double precision frc

      dimension frc(3,nat)

      logical status, yes

      logical grad_same_geom
      external grad_same_geom

      status = rtdb_get ( rtdb, 'gradients:restart', mt_log, 1, yes )

      if ( status .and. yes ) then

        status = rtdb_get ( rtdb, 'gradients:ndone', mt_int, 1, 
     $       ndone )

        if ( status ) then
C         have basis or geometry changed?
          status = grad_same_geom ( rtdb, geom, basis, nat )
        endif

        if ( status ) then
          if ( ga_nodeid() .eq. 0 ) then
C           node 0 gets the old forces
            if ( util_print('information', print_high) ) then
              call util_print_centered ( luout,
     $             'read forces from rtdb', 35, .false.)
              call util_flush(luout)
            endif

            status = rtdb_parallel (.false.)
            status = rtdb_get ( rtdb, 'gradients:partial gradient', 
     $           mt_dbl, 3*nat, frc )
            if ( .not. status ) 
     $           call errquit('could not read partial gradients',1,
     &       RTDB_ERR)
            status = rtdb_parallel (.true.)

          else  ! node <> 0
C           all other nodes simply zero their forces
CCC         call dfill ( 3*nat, 0.D0, frc, 1 )
          endif
C         successful completion
          grad_restart = .true.
          return
        else
C         just print error message and continue to default end
          if ( ga_nodeid() .eq. 0 ) then
            if ( util_print('information', print_medium) ) then
              call util_print_centered ( luout,
     $             'grad_restart: no start indices', 35, .false.)
              call util_print_centered ( luout,
     $             'starting from the beginning', 35, .true.)
            endif
          endif
        endif
      endif
      ndone = 0
      grad_restart = .false.
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine zero_forces ( force, oactive, nat )

      double precision force, zero
      logical oactive
      integer nat

      parameter ( zero = 0.D0 )

      dimension force(3,nat)
      dimension oactive(nat)

      do i=1, nat
        if ( .not.oactive(i) ) then
          force(1,i) = zero
          force(2,i) = zero
          force(3,i) = zero
        endif
      enddo

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine grad_finished ( ofinish )
C     sets ofinish true if it is true on any node
      implicit none

#include "msgids.fh"
#include "global.fh"

      logical ofinish
      integer itrue

      if ( ofinish ) then
        itrue = 1
      else
        itrue = 0
      endif

      call ga_igop ( Msg_grad_finished, itrue, 1, 'max' )

      if ( itrue .ne. 0 ) ofinish = .true.

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function grad_same_geom ( rtdb, geom, basis, nat )

#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
C     unfinished routine for automaic restart
C     check if geometry has changed since the partial forces were stored
C     put on hold for release
      grad_same_geom = .true.
      return
      end
