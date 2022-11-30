      subroutine sifc2f( aoint2, info, ierr )
c
c  close the 2-e integral file.
c
c  input: aoint2  = unit number of the aoints2 file.
c         info(*) = info array for this file.
c
c  output: ierr = error return code. 0 for normal return.
c
c  the correct operation order in the calling program is:
c     open(unit=aoints,...)        # standard open for the 1-e file.
c     call sifo2f(..aoint2.)       # open the 2-e file.
c     call sifc2f(aoint2...)       # close the 2-e file.
c     close(unit=aoints...)        # close the 1-e file.
c
c  this routine, along with sifo2f(), properly account for cases in
c  which only one file at a time is actually used, and for FSPLIT=1
c  cases for which all integral records are on the same file.
c
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
c
      integer aoint2, info(*), ierr
c
      integer fsplit
c
      fsplit = info(1)
      ierr   = 0
c
c     # only close if the file is split.
c
      if ( fsplit .eq. 2 ) then
c
c        # 2-e records are separate.  use async i/o routines.
c        # close the file.
c
         call aiclos( aoint2 )
c        # aiclos() does not return ierr.
         ierr = 0
c
      endif
c
      return
      end
