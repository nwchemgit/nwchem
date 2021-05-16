*deck siffw2
      subroutine siffw2(
     & info,    nipv,    num,     last,
     & itypea,  itypeb,  ifmt,    ibvtyp,
     & values,  labels,  ibitv,   nlist,
     & ierr )
c
c  formatted write of a 2-e buffer.
c
c  input:
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 four orbital labels are packed in each labels(*) entry.
c       = 2 two orbital labels are packed in each labels(*) entry.
c       = 4 one orbital label is stored in each labels(*) entry.
c  num = number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
c                 elements of labels(*) are referenced.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n2max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = error return code. 0 for normal return.
c
c  24-apr-92 nipv added to the output record. -rls
c  16-aug-89 num=0 short return added.
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  nipv,   num,    itypea, itypeb, last,
     & ifmt,   ibvtyp, nlist,  ierr
      real*8   values(*)
      integer  info(*),        labels(nipv,*), ibitv(*)
c
      integer  l2rec,  n2max,  lenpl,  lenbv,  l2recx,
     & lab1,   ifmtv,  i,      j
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l2rec = info(4)
      n2max = info(5)
c
c     # check parameters for consistency...
c
      if(num.gt.n2max)then
         call sifbummer('siffw2: num=',num,wrnerr)
         ierr = -1
         return
      elseif(itypea.lt.0 .or. itypea.gt.7)then
         call sifbummer('siffw2: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif(itypeb.lt.0 .or. itypeb.gt.1023)then
         call sifbummer('siffw2: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif(last.lt.0 .or. last.gt.3)then
         call sifbummer('siffw2: last=',last,wrnerr)
         ierr = -1
         return
      elseif(ibvtyp.lt.0 .or. ibvtyp.gt.7)then
         call sifbummer('siffw2: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if(ifmt.eq.0)then
         lenpl=(num+1)/2
      elseif(ifmt.eq.1)then
         lenpl=num
      else
         call sifbummer('siffw2: ifmt=',ifmt,wrnerr)
         ierr = -1
         return
      endif
      if(ibvtyp.ne.0)then
         lenbv=(n2max+63)/64
      else
         lenbv=0
      endif
      l2recx=(1+num+lenpl+lenbv)
      if( l2recx .gt. l2rec )then
         call sifbummer('siffw2: (l2recx-l2rec)=',(l2recx-l2rec),wrnerr)
         ierr = -1
         return
      endif
c
      lab1 = num + 2
c
c     # write out the dword information.
c
      write(nlist,6010)
     & num,    lab1,   ibvtyp, itypea,
     & itypeb, ifmt,   last,   nipv
6010  format(1x,8i7)
c
      if ( nipv .eq. 1 ) then
         assign 6021 to ifmtv
      elseif( nipv .eq. 2 ) then
         assign 6022 to ifmtv
      elseif( nipv .eq. 4 ) then
         assign 6024 to ifmtv
      else
         call sifbummer('siffw2: nipv=',nipv,wrnerr)
         ierr = -1
         return
      endif
c
      do 10 i = 1, num
         write(nlist,ifmtv) values(i), (labels(j,i), j=1,nipv)
10    continue
6021  format(1x,1pe20.12, i11 )
6022  format(1x,1pe20.12, 2i7 )
6024  format(1x,1pe20.12, 4i4 )
c
      if(ibvtyp.ne.0)then
         write(nlist,6030) ( ibitv(i), i = 1, num )
      endif
6030  format(1x,20i2)
c
      return
      end
*deck sifmge
      subroutine sifmge( nenrgy, energy, ietype, nnew, enew, typnew )
c
c  merge two energy(*) vectors such that exact duplicates are not
c  repeated.
c
c  input: nenrgy    = number of current energy(*) entries.
c         energy(*) = current energy entries. (distinct)
c         ietype(*) = energy types.
c         nnew      = number of new candidates.
c         enew(*)   = new energy values. (distinct)
c         typnew(*) = new types.
c
c  output: nenrgy    = updated value. may be as large as
c                      (original_nenrgy + nnew).
c          energy(1:nenrgy) = updated values.
c          ietype(1:nenrgy) = updated types..
c
c  this routine should be called by programs that write integral
c  files and which read energy(*) arrays from several sources.  If
c  input energy(*) values from two different sources have some common
c  ancestor, then duplication of common energy(*) elements should be
c  eliminated.  this is not only desirable for total energ values, but
c  it is necessary for frozen-core energy contributions in order to
c  avoid overcounting.
c
c  note that tests for exact equality are used in the comparisons.
c  this means that nothing should be done to input energy(*) values
c  which are passed on to other programs other than copy operations.
c  in particular, the programmer should avoid adding and then
c  subtracting the nuclear repulsion energy from total energy values.
c  it is the original, unmodified, values which should be transferred.
c
c  18-oct-90 written by ron shepard.
c
      implicit none
c
      integer  nenrgy, ietype(*),      nnew,   typnew(*)
      real*8   energy(*),      enew(*)
c
      integer  nold,   inew,   iold
c
c     # current energy(*) values are assumed to be distinct (i.e. no
c     # repetitions of the same energy contribution) among themselves.
c
c     # new contributions must be distinct among themselves.
c
      nold = nenrgy
c
      do 20 inew = 1, nnew
         do 10 iold = 1, nold
            if ( (energy(iold) .eq. enew(inew)) .and.
     &       ( ietype(iold) .eq. typnew(inew)) ) goto 20
10       continue
c        # loop exit means a distinct new energy has been found.
         nenrgy = nenrgy + 1
         energy(nenrgy) = enew(inew)
         ietype(nenrgy) = typnew(inew)
20    continue
c
      return
      end
*deck sifo2f
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
         call sifaiopen( aoint2, filnm2, l2rec )
c        # aiopen() does not return ierr.
         ierr = 0
c
      endif
c
      return
      end
*deck sifpre
      subroutine sifpre( nlist, nenrgy, energy, ietype )
c
c  print out the energy(*) array with character labels.
c
c  24-apr-92 cnvginf values added. -rls
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit none
c
      integer  nlist,  nenrgy, ietype(*)
      real*8   energy(*)
c
      integer  i,      itypea, itypeb
      character*8      chrtyp
c
      do 10 i = 1, nenrgy
c
         itypea = ietype(i) / 1024
         itypeb = mod( ietype(i), 1024 )
         call siftyp( itypea, itypeb, chrtyp )
c
         if ( ietype(i) .ge. 0 ) then
c
c           # fcore from a 1-e hamiltonian array.
c
            write(nlist,6010) i, energy(i), ietype(i), 'fcore', chrtyp
c
         elseif ( itypea .eq. 0 ) then
c
c           # core energy value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'core', chrtyp
c
         elseif ( itypea .eq. -1 ) then
c
c           # total energy value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'total', chrtyp
c
         elseif ( itypea .eq. -2 ) then
c
c           # energy or wave function convergence value.
c
            write(nlist,6010) i, energy(i), ietype(i), 'cnvginf', chrtyp
c
         else
c
c           # undefined energy type.
c
            write(nlist,6010) i, energy(i), ietype(i), 'unknown', chrtyp
c
         endif
10    continue
c
      return
6010  format(' energy(', i2, ')=', 1pe20.12, ', ietype=', i5,
     & ', ', a7, ' energy of type: ', a )
      end
*deck sifr1n
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
c
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
      call ifill( nndxf(nsym+1), 0, kntin, 1 )
**      call izero( nndxf(nsym+1), kntin, 1 )
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
c
      return
      end
*deck sifr1x
      subroutine sifr1x(
     & aoints,  info,    itypea,  itypeb,
     & nsym,    nbpsy,   symoff,  buffer,
     & values,  labels,  mapin,   symb,
     & lena,    array,   fcore,   kntin,
     & lasta,   lastb,   last,    nrec,
     & ierr )
c
c  read the designated 1-e integral array and accumulate the matrix
c  elements into array(*).
c
c  this is a no-frills routine to read generic integral arrays.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  itypea, itypeb = generic and specific integral type to be read.
c                   if itypea>=0, then search the entire 1-e integral
c                                 file for candidates.  the file is
c                                 left positioned at the end of the
c                                 1-e integral records.
c                   if itypea<0,  then read the next single array only.
c                                 the file is not repositioned before
c                                 the read.
c                                 the file is left positioned at the end
c                                 of this array.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  symoff(*) = symmetry block offsets used for nontotally symmetric
c              operators.  This includes itypea=1 and itypea=2 arrays.
c              the symmetry blocks within array(*) are offset by
c              symoff( nndxf(isym) + jsym). this allows compact
c              storage of sparse nonsymmetric operator matrices.
c              itypea=0 arrays are referenced in the normal manner,
c              and symoff(*) is not used.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  mapin(*) = input_ao-to-ao mapping vector.
c  lena = length of array(*).
c  array(1:lena) = initial values for the integral array.
c  fcore = initial value for the frozen core contribution.
c
c  output:
c  symb(1:nbft) = symmetry index of each basis function.
c  array(*) = designated 1-e integral array (only if ierr = 0).
c             symmetry block order is determined by the symoff(*) array.
c             diagonal symmetry blocks are stored lower-triangle packed
c             by rows.  off-diagonal symmetry blocks are stored in
c             standard column order.
c  fcore = updated with core contributions from the designated arrays.
c  kntin(*) = number of elements read in each symmetry block.  this
c              array is always referenced as ( nndxf(isym) + jsym ).
c  lasta, lastb = last integral type read and used.
c  last = last parameter from the last record read.
c  ierr = error return code.
c       =  0     for normal return.
c       = -1     eof was found on aoints.
c       = -2     unsatisfied search.
c       = -n*100 if n symmetry blocking errors were detected.
c       >  0     iostat error while reading aoints..
c
c  08-oct-90 (columbus day) 1-e fcore added. sifr1n() interface
c            used.  ierr added. -rls
c  01-jul-90 multiple matrix capability added. -rls
c  15-aug-89 written by ron shepard.
c
      implicit none
c
      integer   nipv,   nmsame,   nomore,   iretbv
      parameter(nipv=2, nmsame=1, nomore=2, iretbv=0)
c
      integer  aoints, itypea, itypeb, nsym,   lena,
     & lasta,  lastb,  last,   nrec,   ierr
      integer  info(*),        nbpsy(nsym),    labels(nipv,*),
     & mapin(*),       symb(*),        symoff(*),      kntin(*)
      real*8   fcore,  buffer(*),      values(*),      array(*)
c
c     # local:
      integer    btypmx
      parameter( btypmx=20 )
      integer btypes(0:btypmx)
      real*8 fcorex(1)
c
      if ( itypea .ge. 0 ) then
         if ( (itypeb .lt. 0) .or. (itypeb .gt. btypmx) ) then
c           # unsupported itypeb value.
            ierr = -3
            return
         else
c           # setup the btypes(*) array for this integral type.
            call ifill( (itypeb+1), 0, btypes, 1 )
**            call izero( (itypeb+1), btypes, 1 )
            btypes(itypeb) = +1
         endif
      endif
c
      fcorex(1) = (0)
c
      call sifr1n(
     & aoints, info,   itypea, itypeb,
     & btypes, buffer, values, labels,
     & nsym,   nbpsy,  symoff, mapin,
     & lena,   array,  fcorex, symb,
     & kntin,  lasta,  lastb,  last,
     & nrec,   ierr )
c
c     # update the frozen core value.
      fcore = fcore + fcorex(1)
c
      return
      end
*deck sifr2
      subroutine sifr2( aoint2, iwait, info, buffer, reqnum, ierr )
c
c  read a 2-e integral record without decoding the contents.
c
c  input:
c  aoint2  = input unit number.
c  iwait   = asynchronous i/o parameter.
c          = 0  don't wait.  use asynch i/o and return without
c               waiting for i/o completion.  the calling program
c               must call sif2w8() before processing the buffer.
c          = 1  wait for i/o completion before returning.
c  info(*) = info array for this file.
c
c  output:
c  buffer(1:l2rec) = input buffer of length l2rec=info(4) to be filled.
c  reqnum = i/o request number for the i/o operation associated
c           with this record.
c  ierr = error return code. 0 for normal return.
c
c  all sif 2-e records should be read by this routine.  this allows
c  local conventions, such as the use of non-fortran i/o, to be
c  localized.  see also sifw2() and sif2w8().
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
      ierr = 0
      if ( fsplit .eq. 1 ) then
c
c        # must use standard fortran i/o.
c
         call sifseqrbf( aoint2, buffer, l2rec )
c        # sifseqrbf() does not return ierr.
         ierr = 0
c
      elseif ( fsplit .eq. 2 ) then
c
c        # 2-e records are separate.  use async i/o routines.
c
c        # airead() and aiwait() do not allow use of reqnum.
         reqnum = 0
         call sifairead( aoint2, buffer, l2rec )
c        # airead() does not return ierr.
         ierr = 0
         if ( iwait .eq. 1 ) call sifaiwait ( aoint2 )
c
      endif
c
      return
      end
*deck sifrd1
      subroutine sifrd1(
     & aoints,  info,    nipv,    iretbv,
     & buffer,  num,     last,    itypea,
     & itypeb,  ifmt,    ibvtyp,  values,
     & labels,  fcore,   ibitv,   ierr )
c
c  read and decode a 1-e integral record.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  nipv = number of integers per value to be returned.
c       = 0 only unpack dword.  values(*), labels(*), and ibitv(*)
c           are not referenced.
c       = 1 return two orbital labels packed in each labels(*) entry.
c       = 2 return one orbital label in each labels(*) entry.
c  iretbv = bit vector request type.
c     if ( iretbv=0 ) then
c         null request, don't return ibitv(*).
c     elseif ( iretbv=ibvtyp ) then
c         request return of the bit-vector of type iretbv.
c     elseif ( iretbv=-1 .and. ibvtyp<>0 ) then
c         return any type of bit-vector that is on the record.
c     else
c        error. requested bit-vector is not available in buffer(*).
c     endif
c  buffer(1:l1rec) = packed  buffer.
c
c  output:
c  num = actual number of values in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format of the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = values (referenced only if nipv.ne.0).
c  labels(1:nipv,1:num) = integral labels
c           (referenced only if nipv.ne.0).
c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if iretbv.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = return code. 0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, nipv,   iretbv, num,   last,
     & itypea, itypeb, ifmt,   ibvtyp, ierr
      integer  info(*),        labels(*),     ibitv(*)
      real*8   buffer(*),      values(*),     fcore
c
      integer  l1rec,  n1max
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l1rec = info(2)
      n1max = info(3)
c
c     # read the input file...
c
      call sifseqrbf ( aoints, buffer, l1rec )
c     # sifseqrbf() does not return error codes.
      ierr = 0
c
c     # unpack the buffer...
c
      call sifd1(
     & info,  nipv,   iretbv, buffer,
     & num,   last,   itypea, itypeb,
     & ifmt,  ibvtyp, values, labels,
     & fcore, ibitv,  ierr )
c
      return
      end
*deck sifrd2
      subroutine sifrd2(
     & aoint2,  info,    nipv,    iretbv,
     & buffer,  num,     last,    itypea,
     & itypeb,  ifmt,    ibvtyp,  values,
     & labels,  ibitv,   ierr )
c
c  read and decode a 2-e integral record.
c
c  input:
c  aoint2 = input file unit number.
c  info(*) = info array for this file.
c  nipv = number of integers per value to be returned.
c       = 0 only unpack dword.  values(*), labels(*), and ibitv(*)
c           are not referenced.
c       = 1 return four orbital labels packed in each labels(*) entry.
c       = 2 return two orbital labels packed in each labels(*) entry.
c       = 4 return one orbital label in each labels(*) entry.
c  iretbv = bit vector request type.
c     if ( iretbv=0 ) then
c         null request, don't return ibitv(*).
c     elseif ( iretbv=ibvtyp ) then
c         request return of the bit-vector of type iretbv.
c     elseif ( iretbv=-1 .and. ibvtyp<>0 ) then
c         return any type of bit-vector that is on the record.
c     else
c        error. requested bit-vector is not available in buffer(*).
c     endif
c  buffer(1:l2rec) = packed  buffer with l2rec=info(4).
c
c  output:
c  num = actual number of values in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format of the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = values (referenced only if nipv.ne.0).
c  labels(1:nipv,1:num) = integral labels
c           (referenced only if nipv.ne.0).
c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
c                 elements of labels(*) are referenced.
c  ibitv(*) = unpacked bit vector (referenced only if iretbv.ne.0).
c             note: as many as ((n2max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoint2, nipv,   iretbv, num,    last,
     & itypea, itypeb, ifmt,   ibvtyp, ierr
      integer  info(*),        labels(*),      ibitv(*)
      real*8   buffer(*),      values(*)
c
      integer  reqnum
c
      integer   iwait
      parameter(iwait=1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # read the input file.
c
      call sifr2( aoint2, iwait, info, buffer, reqnum, ierr )
c
c     # unpack the buffer...
c
      call sifd2(
     & info,  nipv,   iretbv, buffer,
     & num,   last,   itypea, itypeb,
     & ifmt,  ibvtyp, values, labels,
     & ibitv, ierr )
c
      return
      end
*deck sifrh1
      subroutine sifrh1(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    ierr )
c
c  read header_1 from the standard integral file.
c
c  input:
c  aoints = input file unit number.
c
c  output:
c  ntitle = number of titles.
c  nsym = number of symmetry blocks.
c  nbft = total number of basis functions.
c  ninfo = number of record-definition parameters.
c  nenrgy = number of core energies.  the first element must be the
c           nuclear repulsion energy.
c  nmap = number of optional map vectors.
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   ierr
c
c  vrsion = routine library version number.
c  ntitmx = maximum number of titles allowed.
c  ninchk = minimum number of info(*) elements.
c  lrecmx = maximum record length allowed.  this should be consistent
c           with dword bit-packing in the record-writing routines.
c
      integer   vrsion,   ntitmx,    ninchk
      parameter(vrsion=1, ntitmx=20, ninchk=5)

      integer   lrecmx
      parameter(lrecmx=2**16-1)
c
      integer  verin
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      read(aoints,iostat=ierr)
     & verin, ntitle, nsym, nbft, ninfo, nenrgy, nmap
c
      if ( ierr.ne.0 ) then
         return
      elseif ( verin.ne.vrsion ) then
         call sifbummer('sifrh1: (verin-vrsion)=',(verin-vrsion),wrnerr)
         ierr = -2
         return
      elseif ( ntitle.le.0 .or. ntitle.gt.ntitmx ) then
         call sifbummer('sifrh1: ntitle=',ntitle,wrnerr)
         ierr = -3
         return
      elseif ( nsym.ne.1 .and. nsym.ne.2 .and. nsym.ne.4
     &    .and. nsym.ne.8 ) then
         call sifbummer('sifrh1: nsym=',nsym,wrnerr)
         ierr = -4
         return
      elseif ( nbft.le.0 ) then
         call sifbummer('sifrh1: nbft=',nbft,wrnerr)
         ierr = -5
         return
      elseif ( ninfo.lt.ninchk ) then
         call sifbummer('sifrh1: ninfo=',ninfo,wrnerr)
         ierr = -6
         return
      elseif ( nenrgy.le.0 ) then
         call sifbummer('sifrh1: nenrgy=',nenrgy,wrnerr)
         ierr = -7
         return
      elseif ( nmap.lt.0 ) then
         call sifbummer('sifrh1: nmap=',nmap,wrnerr)
         ierr = -8
         return
      endif
c
      return
      end
*deck sifrh2
      subroutine sifrh2(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    title,
     & nbpsy,   slabel,  info,    bfnlab,
     & ietype,  energy,  imtype,  map,
     & ierr )
c
c  read header_2 from the standard integral file.
c
c  input:
c  aoints = input file unit number.
c  ntitle = number of titles.
c  nsym = number of symmetry blocks.
c  nbft = total number of basis functions.
c  ninfo = number of record-definition parameters.
c  nenrgy = number of core energies.  the first element must be the
c           nuclear repulsion energy.
c  nmap = number of optional map vectors.
c
c  output:
c  title*80(1:ntitle) = identifying titles.
c  nbpsy(1:nsym) = number of basis functions per symmetry block.
c  slabel*4(1:nsym) = symmetry labels.
c  info(1:ninfo) = record-definition parameters.
c  bfnlab*8(1:nbft) = basis function labels.
c  ietype(1:nenrgy) = core energy types.
c  energy(1:nenrgy) = core energies.
c  imtype(1:nmap) = map vector types.
c  map(1:nbft,1:nmap) = basis function map vectors.
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, ierr
      character*80 title(*)
      integer nbpsy(*)
      character*4 slabel(*)
      integer info(*)
      character*8 bfnlab(*)
      integer ietype(*)
      real*8 energy(*)
      integer imtype(*)
      integer map(*)
c
c
      integer mapdim, nbftx, isym
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      mapdim = max ( nmap, 1 )
      call sifzh2(
     & aoints, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   title,
     & nbpsy,  slabel, info,   bfnlab,
     & ietype, energy, imtype, map,
     & mapdim, ierr )
c
      if ( ierr.ne.0 ) return
c
      nbftx=0
      do 10 isym=1,nsym
         if ( nbpsy(isym) .lt. 0 ) then
            call sifbummer('sifrh2: nbpsy(isym)=',nbpsy(isym),wrnerr)
            ierr = -2
            return
         endif
         nbftx=nbftx+nbpsy(isym)
10    continue
      if ( nbftx.ne.nbft ) then
         call sifbummer('sifrh2: (nbftx-nbft)=',(nbftx-nbft),wrnerr)
         ierr = -3
         return
      endif
c
      return
      end
*deck sifrsh
      subroutine sifrsh(
     & aoints,  info,    buffer,  values,
     & labels,  nsym,    nbpsy,   mapin,
     & nnbft,   s1h1,    score,   hcore,
     & symb,    ierr )
c
c  read the overlap and 1-e hamiltonian matrices.
c
c  this is a basic, no-frills, routine to read the 1-e integral arrays
c  necessary for energy calculations.
c
c  on exit, the integral file is positioned after the
c  last 1-e integral record.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  mapin(*) = input_ao-to-ao mapping vector.
c  nnbft = leading dimension of s1h1(*,1:2).
c
c  output:
c  s1h1(*) = the overlap s1(*) matrix is returned in s1h1(*,1) and
c            the total h1(*) matrix is returned in s1h1(*,2).
c            both are symmetry-blocked lower-triangle packed by rows.
c            all 1-e contributions on the integral file are summedc
c            into this array.  consequently, the entries on the file
c            must be only the distinct array elements.
c  score = frozen core contribution. tr( s1 * dfc ) = nfrzct
c  hcore = frozen core contribution. tr( h1 * dfc ) = total_core_energy
c  symb(1:nbft) = symmetry index of each basis function
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 if no arrays were found on the integral file.
c       = -n if n symmetry blocking errors were detected.
c       >  0 for iostat error.
c
c  08-oct-90 (columbus day) 1-e fcore change.  sifr1n() interface
c            used. ierr added. -rls
c  04-oct-90 sifskh() call added. -rls
c  26-jul-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, nsym,   nnbft,  ierr
      integer  info(*),        nbpsy(nsym),    labels(2,*),
     & mapin(*),       symb(*)
      real*8   score,  hcore
      real*8   buffer(*),      values(*),      s1h1(nnbft,2)
c
c     # local...
      integer    itypea,   btypmx
      parameter( itypea=0, btypmx=6 )
c
      integer i,      nntot,  isym,   nrec,   last,   lastb,  lasta
      integer symtot(36), idummy(1), btypes(0:btypmx)
      real*8  fcore(2)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
c
c     # set the btypes(*) array:
c     #            0:s1, 1:t1, 2:v1, 3:vec, 4:vfc, 5:vref, 6:generic_h1
      data btypes/ 1,    2,    2,    2,     2,     2,      2          /
c
      nndxf(i) = (i * (i - 1)) / 2
c
c     # nntot = the actual number of elements in the arrays.
      nntot=0
      do 20 isym=1,nsym
         nntot = nntot + nndxf( nbpsy(isym) + 1 )
20    continue
c
      if ( nntot .gt. nnbft ) then
c        # inconsistent nnbft value.
         call sifbummer('sifrsh: (nntot-nnbft)=',(nntot-nnbft),wrnerr)
         ierr = -2
         return
      endif
c
c     # initialize the output arrays...
c
      call dfill( nntot, 0.0d0, s1h1(1,1), 1 )
**      call wzero( nntot, s1h1(1,1), 1 )
      call dfill( nntot, 0.0d0, s1h1(1,2), 1 )
**      call wzero( nntot, s1h1(1,2), 1 )
      fcore(1) = (0)
      fcore(2) = (0)
c
      call sifr1n(
     & aoints, info,   itypea, btypmx,
     & btypes, buffer, values, labels,
     & nsym,   nbpsy,  idummy, mapin,
     & nnbft,  s1h1,   fcore,  symb,
     & symtot, lasta,  lastb,  last,
     & nrec,   ierr )
c
c     # save the appropriate core values.
      score = fcore(1)
      hcore = fcore(2)
c
      return
      end
*deck sifsce
      function sifsce( nenrgy, energy, ietype )
c
c  select and sum the core energies in the energy(*) array.
c
c  usage: real*8   sifsce
c         external sifsce
c         total_core = sifsce(...)
c
c  input: nenrgy = number of energy(*) values.
c         ietype(1:nenrgy) = energy types.
c         energy(1:nenrgy) = energy array.
c
c  output: sifsce = total core energy such that
c                       total_potential = total_electronic + sifsce(...)
c                   is the total clamped-nucleus, born-oppenheimer
c                   potential.
c
      implicit none
c
      integer  nenrgy, ietype(nenrgy)
      real*8   sifsce, energy(nenrgy)
c
      integer  i,      itypea
      real*8   ecore
c
      ecore = (0)
      do 10 i = 1, nenrgy
c
c        # 0 <= ietype elements are 1-e hamiltonian frozen core terms.
c        # ietype < 0  elements with itypea=0 are other core terms to
c                      be added.  e.g. ietype=-1 = nuclear repulsion.
c
         itypea = ietype(i) / 1024
         if ( (ietype(i) .ge. 0 ) .or.
     &    ( (ietype(i) .lt. 0) .and. (itypea .eq. 0) ) ) then
            ecore = ecore + energy(i)
         endif
10    continue
c
      sifsce = ecore
c
      return
      end
*deck sifsk1
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
     &       wdummy, idummy, ierr )
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
         call sifbummer('sifsk1: fsplit=',fsplit,wrnerr)
         ierr = -1
         return
      endif
c
      return
      end
*deck sifskh
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
*deck sifstv
      subroutine sifstv(
     & aoints,  info,    buffer,  values,
     & labels,  nsym,    nbpsy,   mapin,
     & nnbft,   stv,     score,   tcore,
     & vcore,   symb,    ierr )
c
c  read the overlap, kinetic energy, and potential energy arrays.
c
c  this is a basic, no-frills, routine to read the 1-e integral arrays
c  necessary for energy and virial ratio calculations.
c
c  on exit, the integral file is positioned after the
c  last 1-e integral record.
c
c  input:
c  aoints = input file unit number.
c  info(*) = info array for this file.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  mapin(*) = input_ao-to-ao mapping vector.
c  nnbft = leading dimension of stv(*,1:3).
c
c  output:
c  stv(*) = the overlap s1(*) matrix is returned in stv(*,1), the total
c           kinetic matrix T1_total(*) is returned in stv(*,2), and the
c           total 1-e potential matrix V1_total(*) is returned in
c           stv(*,3).
c           all arrays are returned symmetry-blocked lower-triangle
c           packed by rows.
c           all 1-e contributions on the integral file are summed
c           into this array.  consequently, the entries on the file
c           must be only the distinct array elements.
c  score = frozen core contribution. tr( s1 * dfc ) = nfrzct
c  tcore = frozen core contribution. tr( t1 * dfc )
c  vcore = frozen core contribution. tr( v1 * dfc )
c  symb(1:nbft) = symmetry index of each basis function
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 if no arrays were found on the integral file.
c       = -n if n symmetry blocking errors were detected.
c       >  0 for iostat error.
c
c  31-oct-90 sifstv() created from sifrsh(). -rls
c  08-oct-90 (columbus day) 1-e fcore change.  sifr1n() interface
c            used. ierr added. -rls
c  04-oct-90 sifskh() call added. -rls
c  26-jul-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, nsym,   nnbft,  ierr
      integer  info(*),        nbpsy(nsym),    labels(2,*),
     & mapin(*),       symb(*)
      real*8   score,  tcore,  vcore
      real*8   buffer(*),      values(*),      stv(nnbft,3)
c
c     # local...
      integer    itypea,   btypmx
      parameter( itypea=0, btypmx=6 )
c
      integer i,      nntot,  isym,   nrec,   last,   lastb,  lasta
      integer symtot(36), idummy(1), btypes(0:btypmx)
      real*8  fcore(3)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
c
c     # set the btypes(*) array:
c     #            0:s1, 1:t1, 2:v1, 3:vec, 4:vfc, 5:vref, 6:generic_h1
      data btypes/ 1,    2,    3,    3,     3,     3,      3          /
c
c     # note that generic_h1(*) contributions are added to the total 1-e
c     # potential, and any corresponding fcore contribution is included
c     # into vcore.  This may not be best in all cases, but these
c     # terms must be included somewhere.  -rls
c
      nndxf(i) = (i * (i - 1)) / 2
c
c     # nntot = the actual number of elements in the arrays.
      nntot=0
      do 20 isym=1,nsym
         nntot = nntot + nndxf( nbpsy(isym) + 1 )
20    continue
c
      if ( nntot .gt. nnbft ) then
c        # inconsistent nnbft value.
         call sifbummer('sifrsh: (nntot-nnbft)=',(nntot-nnbft),wrnerr)
         ierr = -2
         return
      endif
c
c     # initialize the output arrays...
c
      call dfill( nntot, 0.0d0, stv(1,1), 1 )
**      call wzero( nntot, stv(1,1), 1 )
      call dfill( nntot, 0.0d0, stv(1,2), 1 )
**      call wzero( nntot, stv(1,2), 1 )
      call dfill( nntot, 0.0d0, stv(1,3), 1 )
**      call wzero( nntot, stv(1,3), 1 )
      fcore(1) = (0)
      fcore(2) = (0)
      fcore(3) = (0)
c
      call sifr1n(
     & aoints, info,   itypea, btypmx,
     & btypes, buffer, values, labels,
     & nsym,   nbpsy,  idummy, mapin,
     & nnbft,  stv,    fcore,  symb,
     & symtot, lasta,  lastb,  last,
     & nrec,   ierr )
c
c     # save the appropriate core values.
      score = fcore(1)
      tcore = fcore(2)
      vcore = fcore(3)
c
      return
      end
*deck siftdy
      subroutine siftdy( chrtdy )
c
c  return the machine_name-time-date-year character string.
c
c  output: chrtdy = character*40 representation of m-t-d-y.
c
c  it is recommended that mtdy-stamps be included on all SIFS files,
c  along with the program name that created the file.  the chrtdy
c  string returned is limited to 40 characters so that both sets of
c  information will fit within one title(*) record:
c
c     title(1:40)=program_name_and_function ; title(41:80)=chrtdy
c
c  since the generation of this data is very machine dependent, local
c  copies of this routine, or the local interface routines on which
c  it relys, may have to be maintained on each machine at a
c  particular site.
c
c  recommended format: machine_id   [wkd] hh:mm:ss [TMZ] dd-mmm-yyyy
c
c  [wkd] is optional and is the 3-char day of the week,
c  hh is the 24 hour time, [TMZ] is optional and is the local time
c  zone, and mmm is the 3-char month abbrieviation.
c  Please do not use the 2-digit month value so as to avoid confusion
c  between dd-mm and mm-dd.  The order of these date parameters is not
c  important, as long as they are unambigious.  yyyy is either the
c  4-digit gregorian year, or its 2-digit abbreviation.
c
c  03-sep-91 cray uname call added. -rls
c  13-mar-91 posix code added. -rls
c  08-oct-90 (columbus day) written by ron shepard.
c
      implicit logical(a-z)
c
      character*(*) chrtdy
c
c     # warning: case dependent code
c     # please do not change the case of the parameter constants and
c     # character data arrays in this program.
c
      character*(*) site
*mdc*if argonne
*      parameter( site = 'ANL' )
*mdc*elseif fsu
*      parameter( site = 'FSU' )
*mdc*elseif osu
*      parameter( site = 'OSU' )
*mdc*elseif wien
*      parameter( site = 'Wien' )
*mdc*elseif pnnl
       parameter (site = 'PNNL')
*mdc*else
*      parameter( site = 'Site=?' )
*mdc*endif
c
*mdc*if posix
*c
*      integer jutsname, ierr, lenn
*      integer iatime(9)
*      character*3 month(12)
*c
*      data month/
*     & 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
*     & 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /
*c
*c     # create the utsname structure.
*      call f77structcreate( 'utsname', jutsname, ierr )
*      if ( ierr .ne. 0 ) then
*c        # encode the error number and return.
*         write( chrtdy, '(a,i9)' ) 'f77structcreate() error=', ierr
*         return
*      endif
*c
*c     # set the values of the utsname structure.
*      call f77uname( jutsname, ierr )
*      if ( ierr .ne. 0 ) then
*c        # encode the error number and return.
*         write( chrtdy, '(a,i9)' ) 'f77uname() error=', ierr
*         goto 9999
*      endif
*c
*c     # extract the nodename element, truncating at 20 characters.
*      call f77strget( jutsname, 'nodename', chrtdy(1:20), lenn, ierr )
*      if ( ierr .ne. 0 ) then
*c        # encode the error number and return.
*         write( chrtdy, '(a,i9)' ) 'f77strget() error=', ierr
*         goto 9999
*      endif
*c
*c     # make sure undefined characters are set to spaces.
*      chrtdy( min(lenn,20)+1 :) = ' '
*c
*c     # get the current time and date.
*      call f77localtime( iatime, ierr )
*      if ( ierr .ne. 0 ) then
*c        # encode the error number and return.
*         write( chrtdy(21:40), '(a11,i9)' ) 'time_error=', ierr
*         goto 9999
*      endif
*c
*c        # return time in the format "hh:mm:ss dd-mmm-yyyy".
*      write( chrtdy(21:40),
*     & '(i2.2,a1,i2.2,a1,i2.2,1x,i2.2,a1,a3,a1,i4.4)' )
*     & iatime(3), ':', iatime(2), ':', 'iatime(1)',
*     & iatime(4), '-', month(iatime(5)), '-', iatime(6)
*c
*9999  continue
*c     # structure cleanup.
*      call f77structfree( jutsname, ierr )
*mdc*elseif vax
*c     # can't get the machine_id easily, so punt.
*      chrtdy = site // ' vax'
*      call time( chrtdy( 22:29) )
*      call date( chrtdy(32:) )
*mdc*elseif cray
*c     # 03-sep-91 unicos 6.0 uname() call added. -rls
*c     # 23-oct-90 cray version written by ron shepard.
*c     # can't access unicos date() from fortran.
*c     # must use clock() and date() instead.
*c     # braindamaged clock() and date() return hollerith strings. -rls
*      integer  clock, date, uname, strlen
*      external clock, date, uname, strlen
*      integer     isys, inode, imach
*      character*9 sys,  node,  mach
*      integer mm
*c     # need to convert ambigous  mm/dd/yy into dd-mmm-yy.
*      character*8 cdate8
*      character*3 month(12)
*      data month/
*     & 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
*     & 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /
*c
*c     # get the system and node name; use cdate8 for scratch.
*      call uname( sys, node, cdate8, cdate8, mach )
*c     # pack into the first 21 spaces of chrtdy(:).
*      isys  = strlen( sys )
*      inode = strlen( node )
*      imach = strlen( mach )
*      chrtdy = node(:inode) // ' ' // sys(:isys) // ' ' // mach(:imach)
*c     # pack time and date info into chrtdy(:).
*      write( chrtdy(22:30), '(a1,a8)'   ) ' ', clock()
*      write( cdate8,        '(a8)'      ) date()
*      read(  cdate8(1:2),   '(i2)'      ) mm
*      chrtdy(32:40) = cdate8(4:5) // '-'
*     & //             month(mm)   // '-' // cdate8(7:8)
*mdc*elseif delta paragon
*      integer i
*      do 2788 i = 1,40
*         chrtdy(i:i) = '?'
* 2788 continue
*mdc*elseif unix
*c
*c     # generic bsd unix version.
*c     # see also $COLUMBUS/special/unix/ for library
*c     # routines to support the following machines:
*c     #    IBM RS/6000
*c     #    Cray Y-MP, Cray 2 (optional, the cray may also use the mdc
*c     #                       block in this routine.)
*c     #    Fujitsu VP2000
*c
*      integer ierr, i
*c
*      character*24     fdate
*      integer  hostnm
*      external hostnm, fdate
*c
*      ierr = hostnm( chrtdy(1:16) )
*      if ( ierr .ne. 0 ) then
*c        # pack the error code into chrtdy(*)
*         write( chrtdy(1:16),'(a9,i7.7)' )'hostnm()=', ierr
*      endif
*c
*c     # change nulls to spaces.  this code assumes that ichar()
*c     # returns the ascii character value.
*c
*      do 10 i = 1, 16
*         if ( ichar( chrtdy(i:i) ) .eq. 0 ) then
*            chrtdy(i:i) = ' '
*         endif
*10    continue
*      chrtdy(17:) = fdate()
*mdc*if pnnl
       character*256 rjhjunk
       call util_hostname(rjhjunk)
       chrtdy(1:16) = rjhjunk
       call util_date(rjhjunk)
       chrtdy(17:) = rjhjunk
*mdc*else
*c
*c     # default code: just return a dummy string.
*c
*      chrtdy = site // 'Machine=?  ??:??:?? ??-???-?? '
*mdc*endif
c
      return
      end
*deck siftyp
      subroutine siftyp( itypea, itypeb, chrtyp )
c
c  return a character description of the integral type.
c
c  input:
c  itypea, itypeb = generic and specific integral or energy(*) types.
c
c  output:
c  chrtyp = character description. (this should be at least character*8
c           in the calling program.)
c
c  08-oct-90 total energy and convergence types added. -rls
c  01-jun-90 table version. -rls
c  11-aug-89 written by ron shepard.
c
      implicit none
c
c     # dummy...
      integer          itypea,  itypeb
      character*(*)    chrtyp
c
c     # local...
      integer  i
c
      integer  typ1e,  st1e,   end1e
      integer  typ2e,  st2e,   end2e
      integer  typc,   stc,    endc
      integer  typte,  stte,   endte
      integer  typcv,  stcv,   endcv
c
c     # typ1e = number of defined 1-e array types.
c     # typ2e = number of defined 2-e array types.
c     # typc  = number of core energy types.
c     # typte = number of total energy types.
c     # typcv = number of convergence types.
c
      parameter( typ1e = 16, st1e = 1,        end1e = typ1e        )
      parameter( typ2e =  2, st2e = end1e +1, end2e = end1e +typ2e )
      parameter( typc  =  1, stc  = end2e +1, endc  = end2e +typc  )
      parameter( typte = 14, stte = endc  +1, endte = endc  +typte )
      parameter( typcv = 10, stcv = endte +1, endcv = endte +typcv )
c
      integer    ntype
      parameter( ntype = endcv )
c
      integer ltypea(ntype), ltypeb(ntype)
      character*8 lctype(ntype)
c
c     # these array types can be in any order, so new ones can be added
c     # to the end or inserted into the middle as appropriate. -rls
c
c     # warning: case dependent code:  do not change the case of the
c     # following character strings.
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st1e, end1e ) /
     & 0,    0, 'S1(*)',     0,    1, 'T1(*)',
     & 0,    2, 'V1(*)',     0,    3, 'Veff(*)',
     & 0,    4, 'VFC(*)',    0,    5, 'Vref(*)',
     & 0,    6, 'H1(*)',     0,    7, 'D1(*)',
     & 0,    8, 'F(*)',      0,    9, 'Q(*)',
     & 1,    0, 'X(*)',      1,    1, 'Y(*)',
     & 1,    2, 'Z(*)',
     & 2,    0, 'Im(SO:x)',  2,    1, 'Im(SO:y)',
     & 2,    2, 'Im(SO:z)'  /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = st2e, end2e ) /
     & 3,    0, '1/r12',     3,    1, 'd2(*)'     /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stc, endc ) /
     & 0,   -1, 'Nuc.Rep.' /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stte, endte ) /
     & -1,   0, 'SCF',      -1,   -1, 'MCSCF',
     & -1,  -2, 'MRSDCI',   -1,   -3, 'CPF',
     & -1,  -4, 'ACPF',     -1,   -5, 'LCC-SD',
     & -1,  -6, 'MRPT',     -1,   -7, 'Bk',
     & -1,  -8, 'DV1',      -1,   -9, 'DV2',
     & -1, -10, 'EPOPLE',   -1,  -11, 'S.O. CI',
     & -1, -12, 'SR-SDCI',  -1,  -13, 'UCEPA'  /
c
      data ( ltypea(i), ltypeb(i), lctype(i), i = stcv, endcv ) /
     & -2,   0, 'SCF-D.E.', -2,   -1, 'SCF-D.D1',
     & -2,  -2, 'MC-D.E.',  -2,   -3, 'MC-Wnorm',
     & -2,  -4, 'MC-Knorm', -2,   -5, 'MC-ApxDE',
     & -2,  -6, 'Bk-Resid', -2,   -7, 'CI-Resid',
     & -2,  -8, 'CI-D.E.',  -2,   -9, 'CI-ApxDE' /
c
      do 10 i = 1, ntype
         if ( itypea .eq. ltypea(i) ) then
            if ( itypeb .eq. ltypeb(i) ) then
               chrtyp = lctype(i)
               return
            endif
         endif
10    continue
c     # loop exit means unrecognized type.
      chrtyp = 'Unknown'
c
      return
      end
*deck sifw1x
      subroutine sifw1x(
     & aoints,  info,    lstflg,  itypea,
     & itypeb,  ifmt,    mapout,  array,
     & nsym,    nbpsy,   symoff,  kntin,
     & buffer,  values,  labels,  fcore,
     & small,   kntout,  numtot,  nrec,
     & ierr )
c
c  write the designated 1-e integral array.
c
c  the output file must be correctly positioned on entry.
c
c  this is a no-frills routine to write generic 1-e integral arrays.
c
c  input:
c  aoints = output file unit number.
c  info(*) = info array for this file.
c  lstflg = flag value to write for the last record.
c  itypea, itypeb = generic and specific integral types.
c  ifmt = output format.
c  mapout(*) = bfn-to-output_bfn mapping vector.
c  array(*) = array to be output.
c             "standard" order of the symmetry blocks and elements
c             within a block is assumed.  see below for details.
c  nsym = number of symmetry blocks.
c  nbpsy(*) = number of basis functions per symmetry block.
c  symoff(*) = symmetry block offsets used for nontotally symmetric
c              arrays. the elements are referenced as
c              ( nndxf(isym) + jsym ).  this allows compact
c              storage of sparse nonsymmetric operator matrices.
c              symoff(*) is not referenced for itypea=0 totally
c              symmetric operator arrays.
c  kntin(*) = nominally the number of elements in each symmetry block.
c             this routine only uses the fact that an element is
c             nonzero to determine whether to write out the block.
c             kntin(*) is referenced as (nndxf(isym)+jsym).   when
c             combined with symoff(*), this allows compact storage
c             of sparse nonsymmetric operator matrices in which the
c             blocks are stored arbitrarily in memory.
c  buffer(1:l1rec) = record buffer.
c  values(1:n1max) = value buffer.
c  labels(1:2,1:n1max) = orbital label buffer.
c  fcore = frozen core contribution.
c  small = outout integral cutoff threshold.
c  nrec = initial aoints record count.
c
c  output:
c  array(*) = unmodified, but written to the output file.
c  kntout(1:nndx(nsym+1)) = number of output elements in each symmetry
c                           block.  referenced as (nndxf(isym)+jsym).
c  numtot = actual total number of values written to aoints.
c  nrec = updated record count.
c  ierr = error return code.  0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  15-aug-89  written by ron shepard.
c
      implicit none
c
c     # ibvtyp = 0; bit vectors are assumed to be associated with
c     #             individual integral records, and not with
c     #             the array elements. -rls
c
      integer   nipv,   msame,   ibvtyp
      parameter(nipv=2, msame=0, ibvtyp=0)
c
      integer  aoints, lstflg, itypea, itypeb, ifmt,
     & nsym,   numtot, nrec,   ierr
      integer  info(*),        mapout(*),      nbpsy(nsym),
     & symoff(*),      kntin(*),       labels(nipv,*), kntout(*)
      real*8   fcore,  small
      real*8   array(*),       buffer(*),      values(*)
c
      integer  i,      j,      ij,     ij0,    n1max,
     & num,    i2,     isym,   i1,     ijsym,  last,
     & jsym,   j1,     numtx,  skipd,  j2
      integer  idummy(1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer nndxf
      nndxf(i) = ( i * (i-1) ) / 2
c
      ierr   = 0
      n1max  = info(3)
      num    = 0
      numtot = 0
      call ifill( nndxf(nsym+1), 0, kntout, 1 )
**      call izero( nndxf(nsym+1), kntout, 1 )
c
      if ( itypea.eq.0 ) then
c
c        # symmetric diagonal-symmetry-blocked input array corresponding
c        # to a totally symmetric operator.
c
         ij0 = 0
         i2  = 0
         do 130 isym = 1, nsym
            i1 = i2 + 1
            i2 = i2 + nbpsy(isym)
c
            ij  = ij0
            ij0 = ij0 + nndxf( nbpsy(isym) + 1 )
c
            ijsym = nndxf(isym+1)
c           # check to make sure this block is to be written.
            if ( kntin(ijsym) .gt. 0 ) then
c
               do 120 i = i1, i2
                  do 110 j = i1, i
                     ij = ij + 1
                     if ( abs(array(ij)).gt.small) then
                        if ( num.eq.n1max ) then
                           numtot = numtot + num
                           last   = msame
                           call sifew1(
     &                      aoints, info,   nipv,   num,
     &                      last,   itypea, itypeb, ifmt,
     &                      ibvtyp, values, labels, fcore,
     &                      idummy, buffer, nrec,   ierr )
                           if ( ierr .ne. 0 ) return
                           numtot = numtot - num
                        endif
                        kntout(ijsym) = kntout(ijsym) + 1
                        num           = num + 1
                        values(num)   = array(ij)
                        labels(1,num) = mapout(i)
                        labels(2,num) = mapout(j)
                     endif
110               continue
120            continue
            endif
130      continue
      elseif ( itypea.eq.1 .or. itypea.eq.2 ) then
c
c        # symmetric and antisymmetric arrays of nontotally symmetric
c        # operators are both stored the same.
c
         if ( itypea .eq. 1 ) then
c           # diagonal elements will be written.
            skipd = 0
         elseif ( itypea .eq. 2 ) then
c           # diagonal elements will be ignored.
            skipd = 1
         endif
c
         i2 = 0
         do 260 isym = 1, nsym
            i1 = i2 + 1
            i2 = i2 + nbpsy(isym)
c
            j2 = 0
            do 230 jsym = 1, (isym-1)
               j1 = j2 + 1
               j2 = j2 + nbpsy(jsym)
c
               ijsym = nndxf(isym) + jsym
               if ( kntin(ijsym) .gt. 0 ) then
c
                  ij = symoff(ijsym)
                  do 220 j = j1, j2
                     do 210 i = i1, i2
                        ij = ij + 1
                        if ( abs(array(ij)) .gt. small ) then
                           if ( num.eq.n1max ) then
                              numtot = numtot + num
                              last   = msame
                              call sifew1(
     &                         aoints, info,   nipv,   num,
     &                         last,   itypea, itypeb, ifmt,
     &                         ibvtyp, values, labels, fcore,
     &                         idummy, buffer, nrec,   ierr )
                              if ( ierr .ne. 0 ) return
                              numtot = numtot - num
                           endif
                           kntout(ijsym) = kntout(ijsym) + 1
                           num           = num + 1
                           values(num)   = array(ij)
                           labels(1,num) = mapout(i)
                           labels(2,num) = mapout(j)
                        endif
210                  continue
220               continue
               endif
230         continue
c
c           # diagonal symmetry block of a nontotally symmetric array.
c
c           # note that diagonal itypea=2 array elements are ignored.
c
            ijsym = nndxf(isym+1)
            if ( kntin(ijsym) .gt. 0 ) then
c
               do 250 i = (i1+skipd), i2
                  ij = symoff(ijsym) + nndxf(i - i1 + 1)
                  do 240 j = i1, (i-skipd)
                     ij = ij + 1
                     if ( abs(array(ij)) .gt. small ) then
                        if ( num.eq.n1max ) then
                           numtot = numtot + num
                           last   = msame
                           call sifew1(
     &                      aoints, info,   nipv,   num,
     &                      last,   itypea, itypeb, ifmt,
     &                      ibvtyp, values, labels, fcore,
     &                      idummy, buffer, nrec,   ierr )
                           if ( ierr .ne. 0 ) return
                           numtot = numtot - num
                        endif
                        kntout(ijsym) = kntout(ijsym) + 1
                        num           = num + 1
                        values(num)   = array(ij)
                        labels(1,num) = mapout(i)
                        labels(2,num) = mapout(j)
                     endif
240               continue
250            continue
            endif
260      continue
      else
         call sifbummer('sifw1x: unsupported itypea=',itypea,wrnerr)
         ierr = -1
         return
      endif
c
c     # dump the last buffer.
      numtot = numtot + num
      last=lstflg
      call sifew1(
     & aoints, info,   nipv,   num,
     & last,   itypea, itypeb, ifmt,
     & ibvtyp, values, labels, fcore,
     & idummy, buffer, nrec,   ierr )
      if ( ierr .ne. 0 ) return
c
c     check for consistency between kntout(*) and numtot.
      numtx = 0
      do 320 isym = 1, nsym
         do 310 jsym = 1, isym
            numtx = numtx + kntout( nndxf(isym) + jsym )
310      continue
320   continue
      if ( numtx .ne. numtot ) then
         call sifbummer('sifw1x: (numtx-numtot)=',(numtx-numtot),wrnerr)
         ierr = -2
         return
      endif
c
      return
      end
*deck sifw2
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
         call sifseqwbf( aoint2, buffer, l2rec )
c        # sifseqwbf() does not return ierr.
         ierr = 0
c
      elseif ( fsplit .eq. 2 ) then
c
c        # 2-e records are separate.  use async i/o routines.
c
c        # aiwrit() and aiwait() do not use reqnum.
         reqnum = 0
         call sifaiwrit( aoint2, buffer, l2rec )
c        # aiwrit() does not return ierr.
         ierr = 0
         if ( iwait .eq. 1 ) call sifaiwait ( aoint2 )
c
      endif
c
      return
      end
*deck sifwh
      subroutine sifwh(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    title,
     & nbpsy,   slabel,  info,    bfnlab,
     & ietype,  energy,  imtype,  map,
     & ierr )
c
c  write the header records of the standard integral file structure.
c
c  input:
c  aoints = output file unit number.
c  ntitle = number of titles.
c  nsym = number of symmetry blocks.
c  nbft = total number of basis functions.
c  ninfo = number of record-definition parameters.
c  nenrgy = number of core energies.  the first element must be the
c           nuclear repulsion energy.
c  nmap = number of optional map vectors.
c  title*80(1:ntitle) = identifying titles.
c  nbpsy(1:nsym) = number of basis functions per symmetry block.
c  slabel*4(1:nsym) = symmetry labels.
c  info(1:ninfo) = record-definition parameters.
c  bfnlab*8(1:nbft) = basis function labels.
c  ietype(1:nenrgy) = core energy types.
c  energy(1:nenrgy) = core energies.
c  imtype(1:nmap) = map vector types.
c  map(1:nbft,1:nmap) = basis function map vectors.
c
c  output:
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89  written by ron shepard.
c
      implicit none
c
c  vrsion = routine library version number.
c  ntitmx = maximum number of titles allowed.
c  ninchk = minimum number of info(*) elements.
c  lrecmx = maximum record length allowed.  this should be consistent
c           with dword bit-packing in the record-writing routines.
c
      integer   vrsion,   ntitmx,    ninchk
      parameter(vrsion=1, ntitmx=20, ninchk=5)

      integer   lrecmx
      parameter(lrecmx=2**16-1)
c
c     # using (*) dimensions until parameters are checked...
      integer  aoints, ntitle, nsym,   nbft,   ninfo,
     & nenrgy, nmap,   ierr
      character*80     title(*)
      integer  nbpsy(*)
      character*4      slabel(*)
      integer  info(*)
      character*8      bfnlab(*)
      integer  ietype(*)
      real*8   energy(*)
      integer  imtype(*)
      integer  map(*)
c
      integer  isym,   mapdim, nbftx
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # check arguments for validity...
c
      ierr = 0
      if ( ntitle.le.0 .or. ntitle.gt.ntitmx ) then
         call sifbummer('sifwh: ntitle=',ntitle,wrnerr)
         ierr = -1
         return
      elseif ( nsym.ne.1 .and. nsym.ne.2 .and. nsym.ne.4
     &    .and. nsym.ne.8 ) then
         call sifbummer('sifwh: nsym=',nsym,wrnerr)
         ierr = -1
         return
      elseif ( nbft.le.0 ) then
         call sifbummer('sifwh: nbft=',nbft,wrnerr)
         ierr = -1
         return
      elseif ( ninfo.lt.ninchk ) then
         call sifbummer('sifwh: ninfo=',ninfo,wrnerr)
         ierr = -1
         return
      elseif ( info(1).ne.1 .and. info(1).ne.2  ) then
c        # fsplit parameter.
         call sifbummer('sifwh: fsplit=',info(1),wrnerr)
         ierr = -1
         return
      elseif ( info(2).le.0 .or. info(2).gt.lrecmx ) then
c        # l1rec.
         call sifbummer('sifwh: l1rec=',info(2),wrnerr)
         ierr = -1
         return
      elseif ( info(3).le.0 ) then
c        # n1max.
         call sifbummer('sifwh: n1max=',info(3),wrnerr)
         ierr = -1
         return
      elseif ( info(4).le.0 .or. info(4).gt.lrecmx ) then
c        # l2rec.
         call sifbummer('sifwh: l2rec=',info(4),wrnerr)
         ierr = -1
         return
      elseif ( info(5).le.0 ) then
c        # n2max.
         call sifbummer('sifwh: n2max=',info(5),wrnerr)
         ierr = -1
         return
      elseif ( nenrgy.le.0 ) then
         call sifbummer('sifwh: nenrgy=',nenrgy,wrnerr)
         ierr = -1
         return
      elseif ( ietype(1).ne.-1 ) then
         call sifbummer('sifwh: ietype(1)=',ietype(1),wrnerr)
         ierr = -1
         return
      elseif ( nmap.lt.0 ) then
         call sifbummer('sifwh: nmap=',nmap,wrnerr)
         ierr = -1
         return
      endif
c
      nbftx=0
      do 10 isym=1,nsym
         if ( nbpsy(isym) .lt. 0 ) then
            call sifbummer('sifwh: nbpsy(isym)=',nbpsy(isym),wrnerr)
            ierr = -1
            return
         endif
         nbftx=nbftx+nbpsy(isym)
10    continue
      if ( nbftx.ne.nbft ) then
         call sifbummer('sifwh: nbftx=',nbftx,wrnerr)
         ierr = -1
         return
      endif
c
      mapdim = max( nmap, 1 )
      call sifzwh(
     & aoints, vrsion, ntitle, nsym,
     & nbft,   ninfo,  nenrgy, nmap,
     & title,  nbpsy,  slabel, info,
     & bfnlab, ietype, energy, imtype,
     & map,    mapdim, ierr )
c
      return
      end
*deck sifzh2
      subroutine sifzh2(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    title,
     & nbpsy,   slabel,  info,    bfnlab,
     & ietype,  energy,  imtype,  map,
     & mapdim,  ierr )
c
c  low-level header2 reading routine.
c
c  *** this routine should not be called directly by user programs. ***
c
c  see sifrh2() for argument description.
c
c  if ( nmap.eq.0 ) then
c     mapdim=1 ;imtype(*) and map(*) are not referenced.
c  else
c     mapdim=nmap ;imtype(*) and map(*) are read.
c  endif
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, ntitle, nsym,   nbft,   ninfo,
     & nenrgy, nmap,   mapdim, ierr
      character*80     title(ntitle)
      integer  nbpsy(nsym)
      character*4      slabel(nsym)
      integer  info(ninfo)
      character*8      bfnlab(nbft)
      integer  ietype(nenrgy)
      real*8   energy(nenrgy)
      integer  imtype(mapdim)
      integer  map(nbft,mapdim)
c
      if ( nmap.eq.0 ) then
         read(aoints,iostat=ierr)
     &    title,nbpsy,slabel,info,bfnlab,ietype,energy
      else
         read(aoints,iostat=ierr)
     &    title,nbpsy,slabel,info,bfnlab,ietype,energy,imtype,map
      endif
c
      return
      end
*deck sifzwh
      subroutine sifzwh(
     & aoints,  vrsion,  ntitle,  nsym,
     & nbft,    ninfo,   nenrgy,  nmap,
     & title,   nbpsy,   slabel,  info,
     & bfnlab,  ietype,  energy,  imtype,
     & map,     mapdim,  ierr )
c
c  low-level header-writing routine.
c
c  *** this routine should not be called directly by user programs. ***
c
c  this version writes imtype(*) and map(*).
c  see sifwh() for argument description.
c
c  if ( nmap.ne.0 ) then
c     mapdim=nmap
c  else
c     mapdim=1
c  endif
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, vrsion, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   mapdim, ierr
      character*80     title(ntitle)
      integer  nbpsy(nsym)
      character*4      slabel(nsym)
      integer  info(ninfo)
      character*8      bfnlab(nbft)
      integer  ietype(nenrgy)
      real*8   energy(nenrgy)
      integer  imtype(mapdim)
      integer  map(nbft,mapdim)
c
      write(aoints,iostat=ierr)
     & vrsion, ntitle, nsym, nbft, ninfo, nenrgy, nmap
c
      if ( ierr .ne. 0 ) return
c
      if ( nmap .eq. 0 ) then
         write(aoints,iostat=ierr)
     &    title, nbpsy, slabel, info, bfnlab, ietype, energy
      else
         write(aoints,iostat=ierr)
     &    title, nbpsy, slabel, info, bfnlab, ietype, energy,
     &    imtype, map
      endif
c
      return
      end
      subroutine sifairead(unit, buffer, len)
      implicit none
#include "errquit.fh"
      integer unit, len
      double precision buffer(len)
      read(unit,err=10) buffer
      return
 10   call errquit('sifairead: error with unit ', unit, INPUT_ERR)
      end
      subroutine sifaiwrit(unit, buffer, len)
      implicit none
#include "errquit.fh"
      integer unit, len
      double precision buffer(len)
      write(unit,err=10) buffer
      return
 10   call errquit('sifaiwrit: error with unit ', unit, INPUT_ERR)
      end
      subroutine sifaiopen(unit, filename, len)
      implicit none
#include "errquit.fh"
      integer unit, len
      character*(*) filename
      open(unit,file=filename,form='unformatted',status='unknown',
     $     err=10)
      return
 10   call errquit('sifaiopen: error with unit ', unit, INPUT_ERR)
      end
      subroutine sifaiclos(unit)
      implicit none
#include "errquit.fh"
      integer unit
      close(unit,status='keep',err=10)
      return
 10   call errquit('sifaiclos: error with unit ', unit, INPUT_ERR)
      end
      subroutine sifaiwait(unit)
      integer unit
      end
*deck sifplab8
      subroutine sifplab8( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+7)/8)*8 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw8
c
      nuw8=((nuw+7)/8)*8
      if ( nuw8 .ne. 0 ) call pack(p,8,u,nuw8)
#elif defined(SUN)
c  17-apr-91 mask added to fix garbage overwrite bug. -rls/tom kovar
      integer p(2,*)
      integer u(8,*)
c
      integer nuw8
c
      integer   or, lshift, and
      intrinsic or, lshift, and
c
      integer    m8
      parameter( m8=2**8-1 )
c
      integer i, j
      integer mor
      mor(i,j) = or( i, and(j,m8) )
c
      nuw8=(nuw+7)/8
      do 10 i=1,nuw8
         p(1,i)=mor(lshift(mor(lshift(mor(lshift(
     &    u(1,i),8),u(2,i)),8),u(3,i)),8),u(4,i))
         p(2,i)=mor(lshift(mor(lshift(mor(lshift(
     &    u(5,i),8),u(6,i)),8),u(7,i)),8),u(8,i))
10    continue
#else
c  general byte-addressable 32-bit integer machines.
      logical*1 p(8,*)
      logical*1 u(4,8,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
c      
      call util_pack_info(bytes_per_int, obigendian)
      if (.not.obigendian) then
c        # ...little-endian.
         do 10 i=1,((nuw+7)/8)
            p(8,i)=u(1,1,i)
            p(7,i)=u(1,2,i)
            p(6,i)=u(1,3,i)
            p(5,i)=u(1,4,i)
            p(4,i)=u(1,5,i)
            p(3,i)=u(1,6,i)
            p(2,i)=u(1,7,i)
            p(1,i)=u(1,8,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+7)/8)
            p(1,i)=u(4,1,i)
            p(2,i)=u(4,2,i)
            p(3,i)=u(4,3,i)
            p(4,i)=u(4,4,i)
            p(5,i)=u(4,5,i)
            p(6,i)=u(4,6,i)
            p(7,i)=u(4,7,i)
            p(8,i)=u(4,8,i)
20       continue
      endif
#endif
      return
      end
*deck sifulab8
      subroutine sifulab8( p, u, nuw )
c
c  unpack integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+7)/8)*8 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  20-jul-91 cray nuw check added. -galen gawboy/rls
c  written by ron shepard.
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw8
c
      nuw8=((nuw+7)/8)*8
      if ( nuw8 .ne. 0 ) call unpack(p,8,u,nuw8)
#elif defined(SUN)
      integer p(2,*)
      integer u(8,*)
c
      integer i, nuw8
c
      integer   m8
      parameter(m8=2**8-1)
c
      integer   and, rshift
      intrinsic and, rshift
c
      nuw8=(nuw+7)/8
      do 10 i=1,nuw8
         u(1,i)=and(rshift(p(1,i),24),m8)
         u(2,i)=and(rshift(p(1,i),16),m8)
         u(3,i)=and(rshift(p(1,i), 8),m8)
         u(4,i)=and(       p(1,i),    m8)
         u(5,i)=and(rshift(p(2,i),24),m8)
         u(6,i)=and(rshift(p(2,i),16),m8)
         u(7,i)=and(rshift(p(2,i), 8),m8)
         u(8,i)=and(       p(2,i),    m8)
10    continue
#else
c  general byte-addressable 32-bit integer machines.
      logical*1 p(8,*)
      logical*1 u(4,8,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
c      
      call lfill(nuw,.false.,u,1)
**      call izero(nuw,u,1)
      call util_pack_info(bytes_per_int, obigendian)
      if (.not.obigendian) then
c        # ...little-endian.
         do 10 i=1,((nuw+7)/8)
            u(1,1,i)=p(8,i)
            u(1,2,i)=p(7,i)
            u(1,3,i)=p(6,i)
            u(1,4,i)=p(5,i)
            u(1,5,i)=p(4,i)
            u(1,6,i)=p(3,i)
            u(1,7,i)=p(2,i)
            u(1,8,i)=p(1,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+7)/8)
            u(4,1,i)=p(1,i)
            u(4,2,i)=p(2,i)
            u(4,3,i)=p(3,i)
            u(4,4,i)=p(4,i)
            u(4,5,i)=p(5,i)
            u(4,6,i)=p(6,i)
            u(4,7,i)=p(7,i)
            u(4,8,i)=p(8,i)
20       continue
      endif
#endif
      return
      end
*deck sifplab16
      subroutine sifplab16( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+3)/4)*4 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw16
c
      nuw16=((nuw+3)/4)*4
      if ( nuw16 .ne. 0 ) call pack(p,16,u,nuw16)
#else
c     # general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer*2 u(2,4,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
c      
      call util_pack_info(bytes_per_int, obigendian)
      if (.not.obigendian) then
c        # ...little-endian.
         do 10 i=1,((nuw+3)/4)
            p(4,i)=u(1,1,i)
            p(3,i)=u(1,2,i)
            p(2,i)=u(1,3,i)
            p(1,i)=u(1,4,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+3)/4)
            p(1,i)=u(2,1,i)
            p(2,i)=u(2,2,i)
            p(3,i)=u(2,3,i)
            p(4,i)=u(2,4,i)
20       continue
      endif
#endif
      return
      end
*deck sifulab16
      subroutine sifulab16( p, u, nuw )
c
c  unpack one-electron integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u(1: ((nuw+3)/4)*4 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  20-jul-91 cray nuw check added. -galen gawboy/rls
c  written by ron shepard.
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw16
c
      nuw16=((nuw+3)/4)*4
      if ( nuw16 .ne. 0 ) call unpack(p,16,u,nuw16)
#else
c     # general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer*2 u(2,4,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
c      
      call ifill2 ( nuw, 0, u, 1 )
**      call izero ( nuw, u, 1 )
      call util_pack_info(bytes_per_int, obigendian)
      if (.not.obigendian) then
c        # ...little-endian.
         do 10 i=1,((nuw+3)/4)
            u(1,1,i)=p(4,i)
            u(1,2,i)=p(3,i)
            u(1,3,i)=p(2,i)
            u(1,4,i)=p(1,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+3)/4)
            u(2,1,i)=p(1,i)
            u(2,2,i)=p(2,i)
            u(2,3,i)=p(3,i)
            u(2,4,i)=p(4,i)
20       continue
      endif
#endif
      return
      end
*deck sifplab32
      subroutine sifplab32( p, u, nuw )
c
c  pack integral labels from u(*) into p(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+1)/2)*2 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  written by ron shepard.
c  version: 14-mar-89
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw32
#ifdef SOLARIS
c     needed for -stackvar option
      save shortw, longw
#endif
c
      nuw32=((nuw+1)/2)*2
      if ( nuw32 .ne. 0 ) call pack(p,32,u,nuw32)
#else
      integer p(2,*)
      integer u(2,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
c      
      call util_pack_info(bytes_per_int, obigendian)
      if (.not.obigendian) then
c        # ...little-endian.
         do 10 i=1,((nuw+1)/2)
            p(2,i)=u(1,i)
            p(1,i)=u(2,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+1)/2)
            p(1,i)=u(1,i)
            p(2,i)=u(2,i)
20       continue
      endif
#endif
      return
      end
*deck sifulab32
      subroutine sifulab32( p, u, nuw )
c
c  unpack integral labels from p(*) into u(*,*).
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+1)/2)*2 ) are referenced.
c  nuw  = number of unpacked integral labels.
c
c  20-jul-91 cray nuw check added. -galen gawboy/rls
c  written by ron shepard.
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw32
#ifdef SOLARIS
c     needed for -stackvar option
      save shortw, longw
#endif
c
      nuw32=((nuw+1)/2)*2
      if ( nuw32 .ne. 0 ) call unpack(p,32,u,nuw32)
#else
      integer p(2,*)
      integer u(2,*)
c
      integer i
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
c      
      call util_pack_info(bytes_per_int, obigendian)
      if (.not.obigendian) then
c        # ...little-endian.
         do 10 i=1,((nuw+1)/2)
            u(2,i)=p(1,i)
            u(1,i)=p(2,i)
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+1)/2)
            u(1,i)=p(1,i)
            u(2,i)=p(2,i)
20       continue
      endif
#endif
      return
      end
*deck sifplab1
      subroutine sifplab1( p, u, nuw )
c
c  pack a bit vector.
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+63)/64)*64 ) are referenced.
c  nuw  = number of unpacked array elements.
c
c  written by ron shepard.
c  version: 5-jul-89
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw64
c
      nuw64=((nuw+63)/64)*64
      if ( nuw64 .ne. 0 ) call pack(p,1,u,nuw64)
#elif defined(VAX) || defined(FX2800)
c
c  the following code may be used for 32-bit integer, little-endian
c  machines that use vax-type bit operators.
c  30-jul-90  written by ron shepard.
c
      integer p(2,*)
      integer u(*)
c
      integer local1, local2, u0
      integer sor, i, j
      sor(i,j) = ior(ishft(i,1), iand(j,1) )
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1 = 0
         do 11 j = 1, 32
            local1 = sor(local1,u(u0+j))
11       continue
         p(2,i) = local1
c
         u0 = (i-1)*64 + 32
         local2 = 0
         do 12 j = 1, 32
            local2 = sor(local2,u(u0+j))
12       continue
         p(1,i) = local2
c
20    continue
#elif defined(TITAN) || defined(FX8) || defined(FUJITSUVP)
c
c  the following code may be used for 32-bit integer, big-endian
c  machines that use vax-type bit operators.
c  30-nov-89  written by ron shepard.
c
c  the general code breaks the titan compiler. use this version
c  instead. -rls
      integer p(2,*)
      integer u(*)
c
      integer local1, local2, u0
      integer sor, i, j
      sor(i,j) = ior(ishft(i,1), iand(j,1) )
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1 = 0
         do 11 j = 1, 32
            local1 = sor(local1,u(u0+j))
11       continue
         p(1,i) = local1
c
         u0 = (i-1)*64 + 32
         local2 = 0
         do 12 j = 1, 32
            local2 = sor(local2,u(u0+j))
12       continue
         p(2,i) = local2
c
20    continue
#else
c  general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer u(*)
c
      integer i, u0
      integer local1, local2, local3, local4
      integer*2 llocl1(2), llocl2(2), llocl3(2), llocl4(2)
#ifdef F77      
      equivalence (llocl1,local1),(llocl2,local2)
      equivalence (llocl3,local3),(llocl4,local4)
#else
      integer*2 mold
#endif
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
c      
#ifdef SOLARIS
c     needed for -stackvar option
      save local1, llocl1, local2, llocl2,local3, llocl3,
     $     local4, llocl4, shortw, longw
#endif
c
c  define bit-packing functions recursively for maximum pipelining...
      intrinsic mod
      integer lp2,lp4,lp8,lp16
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16
      lp2(i1,i2) = mod(i1,2) * 2 + mod(i2,2)
      lp4(i1,i2,i3,i4) = (lp2(i1,i2) * 4 + lp2(i3,i4))
      lp8(i1,i2,i3,i4,i5,i6,i7,i8)=
     & (lp4(i1,i2,i3,i4) * 16 + lp4(i5,i6,i7,i8))
      lp16(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16)=
     & ( lp8(i1,i2,i3,i4,i5,i6,i7,i8) * 256 +
     & lp8(i9,i10,i11,i12,i13,i14,i15,i16) )
c
      call util_pack_info(bytes_per_int, obigendian)
      if (.not.obigendian) then
c        # ...little-endian.
         do 10 i=1,((nuw+63)/64)
c
            u0=(i-1)*64
            local1=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl1=transfer(local1,mold,bytes_per_int/2)
#endif            
            p(4,i)=llocl1(1)
c
            u0=(i-1)*64+16
            local2=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl2=transfer(local2,mold,bytes_per_int/2)
#endif            
            p(3,i)=llocl2(1)
c
            u0=(i-1)*64+32
            local3=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl3=transfer(local3,mold,bytes_per_int/2)
#endif            
            p(2,i)=llocl3(1)
c
            u0=(i-1)*64+48
            local4=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl4=transfer(local4,mold,bytes_per_int/2)
#endif           
            p(1,i)=llocl4(1)
c
10       continue
      else
c        # ...big-endian.
         do 20 i=1,((nuw+63)/64)
c
            u0=(i-1)*64
            local1=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl1=transfer(local1,mold,bytes_per_int/2)
#endif            
            p(1,i)=llocl1(2)
c
            u0=(i-1)*64+16
            local2=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl2=transfer(local2,mold,bytes_per_int/2)
#endif            
            p(2,i)=llocl2(2)
c
            u0=(i-1)*64+32
            local3=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl3=transfer(local3,mold,bytes_per_int/2)
#endif            
            p(3,i)=llocl3(2)
c
            u0=(i-1)*64+48
            local4=lp16(
     &       u(u0+1 ),u(u0+2 ),u(u0+3 ),u(u0+4 ),
     &       u(u0+5 ),u(u0+6 ),u(u0+7 ),u(u0+8 ),
     &       u(u0+9 ),u(u0+10),u(u0+11),u(u0+12),
     &       u(u0+13),u(u0+14),u(u0+15),u(u0+16))
#ifndef F77      
      llocl4=transfer(local4,mold,bytes_per_int/2)
#endif            
            p(4,i)=llocl4(2)
c
20       continue
      endif
#endif
      return
      end
*deck sifulab1
      subroutine sifulab1( p, u, nuw )
c
c  unpack a bit vector.
c
c  p(*) = packed array (working precision in the calling program).
c  u(*) = unpacked array.  u( 1 : ((nuw+63)/64)*64 ) are referenced.
c  nuw  = number of unpacked array elements.
c
c  written by ron shepard.
c  version: 5-jul-89
c
      implicit logical(a-z)
      integer nuw
#if defined(CRAY)&& !defined(__crayx1)
      real*8 p(*)
      integer u(nuw)
c
      integer nuw64
c
      nuw64=((nuw+63)/64)*64
      if ( nuw64 .ne. 0 ) call unpack(p,1,u,nuw64)
#elif defined(VAX) || defined(FX2800)
c
c  the following code may be used for 32-bit integer, little-endian
c  machines that use vax-type bit operators.
c  30-jul-90  written by ron shepard.
c
      integer p(2,*)
      integer u(*)
c
      integer local1, local2, u0
      integer sand, i, j
      sand(i,j)=iand(ishft(i,-j),1)
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1=p(2,i)
         do 11 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local1,j)
11       continue
c
         u0 = (i-1)*64 + 32
         local2=p(1,i)
         do 12 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local2,j)
12       continue
c
20    continue
#elif defined(FX8) || defined(TITAN) | defined(FJUITSUVP)
c
c  the following code may be used for 32-bit integer, big-endian
c   machines that use vax-type bit operators.
c  11-jun-89  written by ron shepard.
c
c  the general code breaks the alliant compiler. use this instead. -rls
      integer p(2,*)
      integer u(*)
c
      integer local1, local2, u0
      integer sand, i, j
      sand(i,j)=iand(ishft(i,-j),1)
c
      do 20 i=1,((nuw+63)/64)
c
         u0 = (i-1)*64
         local1=p(1,i)
         do 11 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local1,j)
11       continue
c
         u0 = (i-1)*64 + 32
         local2=p(2,i)
         do 12 j=31,0,-1
            u0=u0+1
            u(u0)=sand(local2,j)
12       continue
c
20    continue
#else
c  general byte-addressable 32-bit integer machines.
      integer*2 p(4,*)
      integer u(*)
c
      integer i, ipower, j, u0
      integer local1, local2, local3, local4
      integer*2 llocl1(2), llocl2(2), llocl3(2), llocl4(2)
#ifdef F77      
      equivalence (llocl1,local1),(llocl2,local2)
      equivalence (llocl3,local3),(llocl4,local4)
#else
      integer*2 mold
#endif
c
c     # to determine the little-endian or big-endian
c     # addressing convention.
      integer bytes_per_int
      logical obigendian
#ifdef SOLARIS
c     needed for -stackvar option
      save local1, llocl1, local2, llocl2,local3, llocl3,
     $     local4, llocl4, shortw, longw
#endif
c
      intrinsic mod
c
      local1=0
      local2=0
      local3=0
      local4=0
      call util_pack_info(bytes_per_int, obigendian)
#ifndef F77      
      llocl1=transfer(local1,mold,bytes_per_int/2)
      llocl2=transfer(local2,mold,bytes_per_int/2)
      llocl3=transfer(local3,mold,bytes_per_int/2)
      llocl4=transfer(local4,mold,bytes_per_int/2)
#endif            
c
      if (.not.obigendian) then
c        # ...little-endian.
         u0=0
         do 10 i=1,((nuw+63)/64)
c
            llocl1(1)=p(4,i)
            ipower=2**15
            do 1 j=1,16
               u0=u0+1
               u(u0)=mod(local1/ipower,2)
               ipower=ipower/2
1           continue
c
            llocl2(1)=p(3,i)
            ipower=2**15
            do 2 j=1,16
               u0=u0+1
               u(u0)=mod(local2/ipower,2)
               ipower=ipower/2
2           continue
c
            llocl3(1)=p(2,i)
            ipower=2**15
            do 3 j=1,16
               u0=u0+1
               u(u0)=mod(local3/ipower,2)
               ipower=ipower/2
3           continue
c
            llocl4(1)=p(1,i)
            ipower=2**15
            do 4 j=1,16
               u0=u0+1
               u(u0)=mod(local4/ipower,2)
               ipower=ipower/2
4           continue
c
10       continue
      else
c        # ...big-endian.
         u0=0
         do 20 i=1,((nuw+63)/64)
c
            llocl1(2)=p(1,i)
            ipower=2**15
            do 11 j=1,16
               u0=u0+1
               u(u0)=mod(local1/ipower,2)
               ipower=ipower/2
11          continue
c
            llocl2(2)=p(2,i)
            ipower=2**15
            do 12 j=1,16
               u0=u0+1
               u(u0)=mod(local2/ipower,2)
               ipower=ipower/2
12           continue
c
            llocl3(2)=p(3,i)
            ipower=2**15
            do 13 j=1,16
               u0=u0+1
               u(u0)=mod(local3/ipower,2)
               ipower=ipower/2
13           continue
c
            llocl4(2)=p(4,i)
            ipower=2**15
            do 14 j=1,16
               u0=u0+1
               u(u0)=mod(local4/ipower,2)
               ipower=ipower/2
14           continue
c
20       continue
      endif
#endif
      return
      end
      subroutine sifbummer(msg,code,level)
      implicit none
#include "errquit.fh"
      character*(*) msg
      integer code, level
c
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      if (level .eq. faterr) then
         call errquit(msg, code, UNKNOWN_ERR)
      else
         write(6,*) msg, code
         call util_flush(6)
      endif
c
      end
      subroutine sifseqwbf( unit, buffer, buflen )
c
c----------------------------------------------------------
c this is a primitive i/o routine to be used by all higher
c level routines requiring i/o of this type. machine i/o
c peculiarites are meant to reside in these routines only.
c----------------------------------------------------------
c
c routine name: sifseqwbf
c version: 1.0                        date: 8/24/88
c author: eric stahlberg - ohio state university
c purpose: this routine will write a given vector to a
c          sequential file. the length of the vector is given
c          in minimum necessary real*8 words to access the entire
c          vector
c parameters:
c     unit:   unit number to be written
c     buffer: buffer to be written
c     buflen: number of real*8 words to write to file
c
      implicit integer(a-z)
c
      integer unit, buflen
      real*8 buffer ( buflen )
c
*mdc*if debug
*      print *,'sifseqwbf: unit=, buflen=',unit,buflen
*mdc*endif
*mdc*if (debug .and. osu)
*      print *, (buffer(i),i=1, min(100,buflen) )
*mdc*endif
c
      write (unit) buffer
      return
      end
      subroutine sifseqrbf( unit, buffer, buflen )
c
c----------------------------------------------------------
c this is a primitive i/o routine to be used by all higher
c level routines requiring i/o of this type. machine i/o
c peculiarites are meant to reside in these routines only.
c----------------------------------------------------------
c
c routine name: sifseqrbf
c version: 1.0                        date: 8/24/88
c author: eric stahlberg - ohio state university
c purpose: this routine will read a given vector from a
c          sequential file. the length of the vector is given
c          in minimum necessary real*8 words to access the entire
c          vector
c parameters:
c     unit:   unit number to read from
c     buffer: buffer to transfer information to
c     buflen: number of real*8 words to read from file
c
      implicit integer(a-z)
c
      integer unit, buflen
      real*8 buffer ( buflen )
c
*mdc*if debug
*      print *,'sifseqrbf: unit=, buflen=',unit,buflen
*mdc*endif
c
      read (unit) buffer
c
*mdc*if ( debug .and. osu)
*      print *, (buffer(i), i=1, min(100,buflen) )
*mdc*endif
      return
      end
c $Id$
