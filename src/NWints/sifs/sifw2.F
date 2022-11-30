      subroutine sifw2( aoint2, iwait, info, buffer, reqnum, ierr )
c
c  write a 2-e integral record.  buffer(*) has already been encoded.
c
c  input:
c  aoint2  = output unit number.
c  iwait   = asynchronous i/o parameter.
c          = 0  don't wait.  use asynch i/o and return without
c               waiting for i/o completion.  the calling program
c               must call sif2w8() before reusing the buffer.
c          = 1  wait for i/o completion before returning.
c  info(*) = info array for this file.
c  buffer(1:l2rec) = output buffer of length l2rec=info(4).
c
c  output:
c  reqnum = i/o request number for the i/o operatoin associated
c           with this record.
c  ierr = error return code.  0 for normal return.
c
c  all sif 2-e records should be written by this routine.  this allows
c  local conventions, such as the use of non-fortran i/o, to be
c  localized.  see also sifr2() and sif2w8().
c
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
c
      integer  aoint2, iwait,  reqnum, ierr
      integer  info(*)
      real*8   buffer(*)
c
      integer  fsplit, l2rec
c
      fsplit = info(1)
      l2rec  = info(4)
c
      ierr   = 0
      if ( fsplit .eq. 1 ) then
c
c        # must use standard fortran i/o.
c
         call seqwbf( aoint2, buffer, l2rec )
c        # seqwbf() does not return ierr.
         ierr = 0
c
      elseif ( fsplit .eq. 2 ) then
c
c        # 2-e records are separate.  use async i/o routines.
c
c        # aiwrit() and aiwait() do not use reqnum.
         reqnum = 0
         call aiwrit( aoint2, buffer, l2rec )
c        # aiwrit() does not return ierr.
         ierr = 0
         if ( iwait .eq. 1 ) call aiwait ( aoint2 )
c
      endif
c
      return
      end
