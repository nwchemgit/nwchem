      subroutine grad_store ( rtdb, geom, g_force, chktime, nat, nquad, 
     $     frc, oinit, ofinish, oactive, oskel )

c$Id: grad_store.F,v 1.1 1995-12-07 23:03:20 vg038 Exp $
C=======================================================================
C     Syncs and writes the partially calculated forces to rtdb if either:
C     - checkpointing is requested and the specified time has elapsed. 
C     - ofinish == .true.: gradient calculation is finished
C     oinit == .true.: just initializes timer
C=======================================================================
      implicit none

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
     $     oskel     ! use symmetry
      dimension oactive (nat)

      integer lnquad
      double precision last, now, interval
      logical status

      save last

      if ( oinit ) then
        last = util_wallsec()
        return
      endif

      now = util_wallsec()

      interval = chktime
      if ( ofinish .or. ((now-last) .gt. interval) ) then
C       CREATE A CHECKPOINT

C       accumulate local sums to global force array
        call ga_acc ( g_force, 1, 3, 1, nat, frc, 3, 1.D0 )
        call ga_sync()

C       find out highest atom quadrupel
        lnquad = nquad
        call ga_igop ( Msg_grad_nquad, lnquad, 1, 'max' )

C       find out if some other process has finished the loops
        call grad_finished ( ofinish )
C       process zero only gets the accumulated forces
        status = rtdb_parallel (.false.)
        if (ga_nodeid().eq.0) then

          call ga_get ( g_force, 1, 3, 1, nat, frc, 3 )
          if ( ofinish ) then
C           only in the end (not for checkpoint)

C           zero force contributions on inactive atoms
            call zero_forces ( frc, oactive, nat )

C           symmetrize
            if (oskel) then
              call sym_grad_symmetrize(geom, frc)
            endif
          endif

          call util_print_centered ( luout,
     $         'writing forces to rtdb', 35, .false.)
          call util_flush(luout)

C         store forces and # of atom quadruples in the database
          status = rtdb_put ( rtdb, 'scf:gradients', MT_DBL, 3*nat, frc)
          if ( .not. status ) 
     $         call errquit('could not store gradients',1)
          status = rtdb_put ( rtdb, 'gradients:ndone', MT_INT, 1, 
     $         lnquad )
        else
C         zero the forces on other nodes (node 0 has them now)
          call dfill ( 3*nat, 0.D0, frc, 1 )
        endif
        status = rtdb_parallel (.true.)

C       zero the global array (the next ga_acc from node 0 contains all)
        call ga_zero ( g_force )
        call ga_sync()

        last = now
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      logical function grad_restart ( rtdb, nat, ndone, frc, g_force )

C=======================================================================
C     Finds out if this is a restart of a partially finished gradient run.
C     If yes, returns the previously calculated forces (proc 0 only)
C     and the last calculated atom quadruple in grad2
C========================================================================
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"

      integer rtdb, nat,
     $     ndone, g_force

      double precision frc

      dimension frc(3,nat)

      logical status, yes

      status = rtdb_get ( rtdb, 'gradients:restart', mt_log, 1, yes )

      if ( status .and. yes ) then
C       delete entry (the next run e.g. in a geometry optimization 
C       is NOT a restart)
        status = rtdb_delete ( rtdb, 'gradients:restart' )

        status = rtdb_get ( rtdb, 'gradients:ndone', mt_int, 1, 
     $       ndone )

        if ( status ) then

C         delete one-electron derivatives (they come from the previous run)
          call ga_zero ( g_force )

          if ( ga_nodeid() .eq. 0 ) then
C           node 0 gets the old forces
            call util_print_centered ( luout,
     $           'read forces from rtdb', 35, .false.)
            call util_flush(luout)

            status = rtdb_parallel (.false.)
            status = rtdb_get ( rtdb, 'scf:gradients', mt_dbl, 3*nat, 
     $           frc )
            if ( .not. status ) 
     $           call errquit('could not read gradients',1)
            status = rtdb_parallel (.true.)

          else  ! node <> 0
C           all other nodes simply zero their forces
C           1el. forces are already in the rtdb
            call dfill ( 3*nat, 0.D0, frc, 1 )
          endif
C         successful completion
          grad_restart = .true.
        else
          if ( ga_nodeid() .eq. 0 ) then
            call util_print_centered ( luout,
     $           'grad_restart: no start indices', 35, .false.)
            call util_print_centered ( luout,
     $          'starting from the beginning', 35, .true.)
          endif
          ndone = 0
          grad_restart = .false.
        endif
      else
        ndone = 0
        grad_restart = .false.
      endif
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
