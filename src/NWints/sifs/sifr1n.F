      subroutine sifr1n(
     & aoints,  info,    itypea,  btypmx,
     & btypes,  buffer,  values,  labels,
     & nsym,    nbpsy,   symoff,  mapin,
     & lda,     array,   fcore,   symb,
     & kntin,   lasta,   lastb,   last,
     & nrec,    ierr )
c
c  read an arbitrary number of 1-e integral arrays.
c
c  this is a basic, no-frills, routine to read several 1-e integral
c  arrays of the same generic type with a single pass through the
c  1-e integral records.  see routines sifr1x() and sifrsh() for
c  examples of the use of this routine.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  itypea = generic integral type to be read.
c         >= 0  read the entire set of 1-e records, and leave the
c               file positioned after the last one.
c         = -1, then read the next array, of any type, from the
c               current file position and return.  The file is
c               left positioned at the beginning of the next array.
c  btypmx = maximum value of itypeb to consider.
c  btypes(0:btypmx) = vector of typeb specific arrays to be read.
c                 btypes(i) <= 0       # ignore this array type.
c                 bytpes(i) =  iarray  # accumulate into array(*,iarray)
c                 if itypea=-1, then iarray=1 always and btypes(*) is
c                 not referenced.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  symoff(*) = symmetry block offsets to use for itypea=1 or itypea=2
c              arrays .  the blocks are offset by
c              symoff( nndxf(isym) + jsym ).
c              itypea=0 arrays are addressed by nnskp(*) defined below.
c  mapin(*) = input_ao-to-ao mapping vector.
c  lda  = leading dimension of array(*,*). ( lda > 0 )
c         for fully dense arrays of itypea=1 or itypea=2, lda should be
c         at least (nbft*(nbft+1))/2.  the calling program is
c         responsible for ensuring that lda and symoff(*) are
c         consistent with the arrays that are being read.
c  array(1:lda,1:*) = initial values of the array values.
c               contributions from the 1-e records are accumulated into
c               this array.  The second dimension is determined by the
c               maximum value of the btypes(*) vector elements.
c  fcore(*) = initial values of frozen core energy contributions.
c             fcore(iarray) is associated with array(*,iarray).
c
c  output:
c  array(*,*) = updated array elements.
c  fcore(*) = updated frozen core values.
c  symb(1:nbft) = symmetry index of each basis function
c  kntin( 1 : (nsym*(nsym+1))/2 ) = number of elements read in each
c             symmetry block, referenced as ( nndxf(isym) + jsym ).
c             this array is to allow the calling program to verify that
c             the input arrays have the expected symmetry blocking of
c             elements.
c  lasta, lastb = last itypea and itypeb array values that were
c                   read and used.
c  nrec = total number of 1-e records processed.
c  ierr = error return code.
c       =  0     for normal return.
c       = -1     eof was found on aoints.
c       = -2     unsatisfied search.
c       = -3     symmetry blocking errors were detected.
ctm
c       = -4     some integral types  are not on file
ctm
c       >  0     iostat error while reading aoints.
c
c  24-apr-92 info(*) declaration fixed. -rls
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  04-oct-90 sifskh() call added. -rls
c  26-jul-89 written by ron shepard.
c
      implicit none
c
      integer   nipv,   msame,   nmsame,   nomore,   iretbv
      parameter(nipv=2, msame=0, nmsame=1, nomore=2, iretbv=0)
      integer  ibtypmx
      parameter ( ibtypmx=40)

c
      integer  fdbtype(0:ibtypmx)
      integer  aoints, itypea, btypmx, nsym,   lda,
     & lasta,  lastb,  last,   nrec,   ierr
      integer  info(*),        kntin(*),       nbpsy(nsym),
     & btypes(0:btypmx),       symoff(*),      mapin(*),       symb(*),
     & labels(nipv,*)
      real*8   buffer(*),      values(*),      array(lda,*),   fcore(*)
c
c     # local...
      integer  i,      j,      akeep,  bkeep,  ntot,   nntot,  isym,
     & ifmt,   itypbx, itypax, num,    iarray, ilab,
     & jlab,   jsym,   ijsym,  ij,     nerror, ibvtyp
      integer  nskp(8),        nnskp(8),       idummy(1)
      logical  check(msame:nomore)
      real*8   fcorex
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
      nndxf(i) = (i * (i - 1)) / 2
c
      ierr = 0
       if (btypmx.gt.ibtypmx)
     . call bummer('sifr1n: btypmx.gt.ibtypmx,btypmx=',btypmx,2)

c
c     # set up nskp(*), nnskp(*), symb(*), and kntin(*)..
c
      ntot  = 0
      nntot = 0
      do 20 isym = 1, nsym
         nskp(isym)  = ntot
         nnskp(isym) = nntot
         ntot        = ntot  + nbpsy(isym)
         nntot       = nntot + nndxf( nbpsy(isym) + 1 )
         do 10 j = (nskp(isym)+1), ntot
            symb(j) = isym
10       continue
20    continue
c
c     # initialize the output symmetry block counter.
      call izero( nndxf(nsym+1), kntin, 1 )
c
c     # position the integral file.
c
      if ( itypea .ge. 0 ) then
         call sifskh( aoints, ierr )
         if ( ierr .ne. 0 ) return
      endif
c
c     # if ( itypea .ge. 0 ) then
c     #    do while ( (last .eq. msame) .or. (last .eq. nmsame) )
c     # else
c     #    do while ( last .eq. msame )
c     # endif
c
      if ( itypea .ge. 0 ) then
         check(msame)  = .true.
         check(nmsame) = .true.
         check(nomore) = .false.
      else
         check(msame)  = .true.
         check(nmsame) = .false.
         check(nomore) = .false.
      endif
c
c     # read the 1-e integral records and accumulate the appropriate
c     # contributions.
c
      do i=0,btypmx
      fdbtype(i)=-1
      enddo

      nrec  = 0
      akeep = -1
      bkeep = -1
c
      last = msame
100   continue
      if ( check(last) ) then
         call sifrd1(
     &    aoints, info,   nipv,   iretbv,
     &    buffer, num,    last,   itypax,
     &    itypbx, ifmt,   ibvtyp, values,
     &    labels, fcorex, idummy, ierr )
         if ( ierr .ne. 0 ) return
c
         nrec = nrec + 1
c
c        # set iarray to point into array(*,iarray).
c
         if ( itypea .eq. -1 ) then
c           # single array mode. accumulate into array(*,1) always.
            iarray = 1
         elseif ( itypea .eq. itypax ) then
c           # correct generic array type.
c           # mask btypmx with btypes(*) to see if this array should
c           # be accumulated somewhere.
            if ( itypbx .le. btypmx) then
               iarray = btypes(itypbx)
               fdbtype(itypbx)=1
            else
               iarray = -1
            endif
         else
c           # ignore this array.
            iarray = -1
         endif
c
         if ( iarray .gt. 0 ) then
c
c           # a valid array has been found.
            akeep = itypax
            bkeep = itypbx
c
c           # accumulate the values from this record.
c
c           # map the input orbital labels.
            do 120 i = 1, num
               labels(1,i) = mapin( labels(1,i) )
               labels(2,i) = mapin( labels(2,i) )
120         continue
c
            if ( itypax .eq. 0 ) then
c
c              # symmetric matrix, totally symmetric operator.
c
c              # nnskp(1:nsym) offsets are referenced.
c              # diagonal symmetry blocks are returned lower-triangular
c              # packed by rows, asub( nndxf(i) + j ), where asub(*)
c              # begins at nnskp( isym ).
c
               do 140 i = 1, num
                  ilab = max( labels(1,i), labels(2,i) )
                  jlab = min( labels(1,i), labels(2,i) )
                  isym = symb(ilab)
                  jsym = symb(jlab)
                  ijsym = nndxf(isym) + jsym
                  kntin( ijsym ) = kntin( ijsym ) + 1
                  if ( isym .eq. jsym ) then
c                    # only valid array(*) elements are referenced.
                     ij   = nnskp(isym) +
     &                nndxf( ilab - nskp(isym) ) + jlab - nskp(isym)
                     array(ij,iarray) = array(ij,iarray) + values(i)
                  endif
140            continue
            elseif ( itypax .eq. 1 ) then
c
c              # symmetric matrix, nontotally symmetric operator.
c
c              # symoff(1: nndxf(nsym+1) ) is referenced.
c              # diagonal symmetry blocks of array(*,iarray) are
c              # returned lower-triangular packed.
c              # off-diagonal symmetry blocks are returned as
c              # rectangular arrays, asub(1:nbpsy(isym), 1:nbpsy(jsym)),
c              # where asub(*) begins at symoff( nndxf(isym) + jsym ).
c
               do 160 i = 1, num
                  ilab  = max( labels(1,i), labels(2,i) )
                  jlab  = min( labels(1,i), labels(2,i) )
                  isym  = symb(ilab)
                  jsym  = symb(jlab)
                  ijsym = nndxf(isym) + jsym
                  kntin( ijsym ) = kntin( ijsym ) + 1
                  if ( isym .eq. jsym ) then
c                    # diagonal blocks are stored lower-triangle packed.
                     ij = symoff(ijsym) +
     &                nndxf( ilab - nskp(isym) ) + jlab - nskp(isym)
                  else
                     ij = symoff(ijsym) +
     &                ( jlab - nskp(jsym) - 1 ) * nbpsy(isym) +
     &                ilab - nskp(isym)
                  endif
                  array(ij,iarray) = array(ij,iarray) + values(i)
160            continue
            elseif ( itypax .eq. 2 ) then
c
c              # antisymmetric matrix, nontotally symmetric operator.
c
c              # symoff(1: nndxf(nsym+1) ) is referenced.
c              # diagonal symmetry blocks of array(*,iarray) are
c              # returned lower-triangular packed.
c              # off-diagonal symmetry blocks are returned as
c              # rectangular arrays, asub(1:nbpsy(isym), 1:nbpsy(jsym)),
c              # where asub(*) begins at symoff( nndxf(isym) + jsym ).
c
c              # note that antisymmetric arrays are stored the same way
c              # as symmetric arrays.  In particular, the diagonal
c              # elements, which are zero, are explicitly stored.
c
               do 180 i = 1, num
                  ilab  = labels(1,i)
                  jlab  = labels(2,i)
                  if ( ilab .lt. jlab ) then
                     ilab      = labels(2,i)
                     jlab      = labels(1,i)
c                    # note the sign change upon transposition.
                     values(i) = -values(i)
                  endif
                  isym  = symb(ilab)
                  jsym  = symb(jlab)
                  ijsym = nndxf(isym) + jsym
                  kntin( ijsym ) = kntin( ijsym ) + 1
                  if ( isym .eq. jsym ) then
c                    # diagonal blocks are stored lower-triangle packed.
                     ij = symoff(ijsym) +
     &                nndxf( ilab - nskp(isym) ) + jlab - nskp(isym)
                  else
                     ij = symoff(ijsym) +
     &                ( jlab - nskp(jsym) - 1 ) * nbpsy(isym) +
     &                ilab - nskp(isym)
                  endif
                  array(ij,iarray) = array(ij,iarray) + values(i)
180            continue
            endif
            if ( last .ne. msame ) then
c              # last record of this array.
c              # accumulate the core contribution.
               fcore(iarray) = fcore(iarray) + fcorex
            endif
         endif
         go to 100
      endif
c
      lasta = akeep
      lastb = bkeep
      if ( akeep .eq. -1 ) then
c        # unsatisfied search.
         ierr = -2
      elseif ( akeep .eq. 0 ) then
c        # if the processed array was a totally symmetric operator,
c        # then check for nonzero off-diagonal symmetry block
c        # elements.  if found, then either the integral file is in
c        # error, or there is an inconsistency in mapin(*).
         nerror = 0
         do 240 isym = 2, nsym
            do 220 jsym = 1, (isym-1)
               nerror = nerror + kntin( nndxf(isym) + jsym )
220         continue
240      continue
c        # reset ierr if symmetry errors were detected.
         if ( nerror .ne. 0 )then
            ierr = -3
            return
         endif
      endif

c        # check whether all of the requested integral blocks
c        # have been read successfully


      if ( itypea .gt.-1) then
        do i=1,btypmx
         if( fdbtype(i) .ne.1 .and. btypes(i).gt.0 ) then
         ierr=-4
         return
        endif
        enddo
      endif

      return
      end
