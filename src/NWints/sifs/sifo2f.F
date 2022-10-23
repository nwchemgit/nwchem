      subroutine sifo2f( aoints, iunit2, filnm2, info, aoint2, ierr )
c
c  determine the 2-e integral file unit number, open the file if
c  necessary, and leave the file positioned at the beginning of
c  the 2-e records.
c
c  input:
c  aoints  = unit number of the aoints file.
c            aoints is assumed to be open on entry to this routine.
c  iunit2  = unit number of the optional aoints2 file.
c            if ( aoints = iunit2 ) then the calling program should
c            assume that aoints will be closed with this call.
c  filnm2  = character filename of the aoints2 file.
c  info(*) = info array for this file.
c
c  output:
c  aoint2 = unit number to use for the 2-e integral file.
c         = aoints or iunit2 as determined by the fsplit parameter.
c  ierr   = error return code. 0 for normal return.
c
c  the correct operation order in the calling program is:
c     open(unit=aoints,...)        # standard open for the 1-e file.
c     call sifo2f(..aoint2.)       # open the 2-e file.
c     call sifc2f(aoint2...)       # close the 2-e file.
c     close(unit=aoints...)        # close the 1-e file.
c
c  this routine, along with sifc2f(), properly account for cases in
c  which only one file at a time is actually used, and for FSPLIT=1
c  cases for which all integral records are on the same file.
c
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
c
      integer  aoints, iunit2, aoint2, ierr
      integer  info(*)
      character*(*)    filnm2
c
      integer  fsplit, l2rec
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      fsplit = info(1)
      l2rec  = info(4)
c
      ierr   = 0
      if ( fsplit .eq. 1 ) then
c
c        # everything is on the same file.  assume aoints is already
c        # opened.  set aoint2, position the file, and return.
c
         aoint2 = aoints
         call sifsk1( aoint2, info, ierr )
         if ( ierr .ne. 0 ) then
c           # 1-e records were not properly written.
            return
         endif
c
      elseif ( fsplit .eq. 2 ) then
c
c        # 2-e records are separate.  use async i/o routines.
c        # set aoint2 and open the file.
c
         aoint2 = iunit2
         if ( aoints .eq. iunit2 ) then
c           # can't use the same unit for two files, so close the old
c           # one before opening the new one.
            close ( unit = aoints, iostat=ierr )
            if ( ierr .ne. 0 ) return
         endif
         call aiopen( aoint2, filnm2, l2rec )
c        # aiopen() does not return ierr.
         ierr = 0
c
      endif
c
      return
      end
