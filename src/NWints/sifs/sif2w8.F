      subroutine sif2w8( aoint2, info, reqnum, ierr )
c
c  wait (w8) for completion of any pending i/o operations on the
c  2-e integral file of the i/o request identified by reqnum.
c
c  aoint2  = unit number.
c  info(*) = info array for this file.
c  reqnum  = i/o request number.  this value was returned by the
c            async i/o routines at the initial i/o reqest.
c
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
c
      integer aoint2, info(*), reqnum, ierr
c
      integer fsplit
c
      fsplit = info(1)
c
      ierr = 0
      if ( fsplit .eq. 2 ) then
c
c        # 2-e records are separate.  use async i/o routines.
c        # otherwise, this is just a dummy call.
c
c        # aiwait() does not use reqnum.
         call aiwait ( aoint2 )
c        # aiwait() does not return ierr.
         ierr = 0
c
      endif
c
      return
      end
