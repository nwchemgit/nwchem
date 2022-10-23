      subroutine sifsk1( aoint2, info, ierr )
c
c  skip the 1-e integrals and position the file at the beginning
c  of the 2-e integral records.
c
c  input:
c  aoint2 = input file unit number for the 2-e integral file.
c           note: this is not necessarily the same file as the 1-e
c                 integral file.
c  info(*) = info array for this file.
c
c  output:
c  ierr = error return code.  0 for normal return.
c
c  20-sep-90 sifskh() version. -rls
c  01-aug-89 dword(1) version. -rls
c  24-jul-89 written by ron shepard.
c
      implicit none
c
      integer  aoint2, ierr
      integer  info(*)
c
      integer fsplit, last,   num,    ibvtyp, ifmt,   itypeb, itypea
      integer idummy(1)
      real*8  wdummy(1),      dword(1)
      integer   nipv,   msame,   nomore,   iretbv,   l1rec,   n1max
      parameter(nipv=0, msame=0, nomore=2, iretbv=0, l1rec=1, n1max=1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr = 0
c
      fsplit = info(1)
      if ( fsplit .eq. 1 ) then
c
c        # header1, header2, 1-e records, and 2-e records are all
c        # on the same file.
c
         call sifskh( aoint2, ierr )
         if ( ierr .ne. 0 ) return
c
c        # do while ( last .ne. nomore)...
         last = msame
100      continue
         if ( last .ne. nomore ) then
c
c           # just read and decode the first word of each record.
c
            read(aoint2,iostat=ierr)dword
            if ( ierr.ne.0 ) return
c
            call sifd1(
     &       info,   nipv,   iretbv, dword,
     &       num,    last,   itypea, itypeb,
     &       ifmt,   ibvtyp, wdummy, idummy,
     &       wdummy(1), idummy, ierr )
            if ( ierr.ne.0 ) return
c
            goto 100
         endif
      elseif ( fsplit .eq. 2 ) then
c        # 2-e integrals only are on aoint2; rewind is sufficient.
c        # note: assume that a standard rewind works with async i/o.
c        #       this may need modification later. -rls
         rewind aoint2
      else
         call bummer('sifsk1: fsplit=',fsplit,wrnerr)
         ierr = -1
         return
      endif
c
      return
      end
