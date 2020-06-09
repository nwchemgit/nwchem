      integer function util_time_remaining(rtdb)
* $Id$
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
c
c     Return the minimum of the time remaining w.r.t. limits
c     .   set in the input (if any), and
c     .   set in the the batch job (if any)
c
c     If no limit information is available then return -1
c
c     NB: THIS IS A COLLECTIVE FUNCTION in order to ensure that
c     all processes get the same information and also to avoid 
c     creation of many subprocesses on machines for which it is
c     necessary to analyse the output of shell commands using perl.
c
      integer user_left, batch_left, left
      integer util_batch_job_time_remaining
      double precision util_wallsec
      external util_batch_job_time_remaining, util_wallsec
      logical rtdb_mode
c
      if (ga_nodeid() .eq. 0) then
c
         rtdb_mode = rtdb_parallel(.false.)
         if (rtdb_get(rtdb,'timelimit', mt_int, 1, user_left)) then
*            write(6,*) ' USER LIMIT ', user_left
            user_left = user_left - util_wallsec()
            if (user_left .lt. 0) user_left = 0
         else
            user_left = -1
         endif
*         write(6,*) '  USER LEFT ', user_left
         rtdb_mode = rtdb_parallel(rtdb_mode)
         batch_left = util_batch_job_time_remaining()
*         write(6,*) ' BATCH LEFT ', batch_left
c
         if (user_left .eq. -1) then ! User limit not defined
            left = batch_left
*            write(6,*) ' LEFT1 ', left
         else if (batch_left .eq. -1) then ! Batch limit not defined
            left = user_left
*            write(6,*) ' LEFT2 ', left
         else                   ! Both defined
*            left = min(batch_left,user_left)
            left = user_left    ! Seems better to allow override from input
*            write(6,*) ' LEFT3 ', left
         endif
      endif
c
      call ga_brdcst(161,left,MA_sizeof(MT_INT,1,MT_BYTE),0)
c
      util_time_remaining = left
c
      end
      logical function util_test_time_remaining(rtdb, required)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer rtdb ! [input]
      integer required ! [input] required time in seconds
c
c     Return TRUE if the required no. of seconds remain in the
c     job (minimum of user or batch job limits) or if no limits 
c     have been set.  If insufficient time remains, return FALSE.
c
c     NB: THIS IS A COLLECTIVE FUNCTION
c
      integer left, req
      integer util_time_remaining
      external util_time_remaining
c
      left = util_time_remaining(rtdb)
c
      util_test_time_remaining = .true.
c
      if (left .eq. -1) return
c
      req = required
      call ga_igop(353, req, 1, 'max')
c
      if (left.eq.0 .or. req.gt.left) then
         util_test_time_remaining = .false.
         if (ga_nodeid().eq.0) then
            write(6,1) left, required
 1          format(/' !! Terminating execution since time remaining =',
     $           i7,'s',/,' !!             is less than the ',
     $           'required time =',i7,'s'/)
            call util_flush(6)
         endif
      endif
c
      end
