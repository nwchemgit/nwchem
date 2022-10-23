      subroutine sifskh( aoints, ierr )
c
c  skip over the header records and position the file at the
c  beginning of the 1-e integral records.
c
c  output: ierr = error return code.  0 for normal return.
c
c  20-sep-90 written by ron shepard.
c
      implicit none
c
      integer aoints, ierr
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      rewind aoints
      read( aoints, iostat=ierr )
      if ( ierr .ne. 0 ) return
c
      read( aoints, iostat=ierr )
      return
c
      end
