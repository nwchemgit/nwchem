*
* $Id$
*
*colib8.f
*colib part=8 of 9.  SIFS library routines.
*version=4.1 last modified: 24-apr-92
c
c  RCS $Revision: 1.1 $  $Date: 2000-11-02 23:53:40 $
c
c  see colib1.f for version history info.
c
*  *** this entire file is incremental ***
*deck sif2w8
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
         call sifaiwait ( aoint2 )
c        # aiwait() does not return ierr.
         ierr = 0
c
      endif
c
      return
      end
*deck sifc2f
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
         call sifaiclos( aoint2 )
c        # aiclos() does not return ierr.
         ierr = 0
c
      endif
c
      return
      end
*deck sifcfg
      subroutine sifcfg(
     & itype,   lrecal,  nbft,    ibvtyp,
     & ifmt,    lrec,    nmax,    ierr )
c
c  return a set of consistent configuration parameters for a standard
c  integral file structure.
c
c  input:
c  itype = integral type.
c          1 for 1-e integrals,
c          2 for 2-e integrals.
c  lrecal = maximum buffer length to be allocated.
c         =-1 is a special case for default output values.
c  nbft = total number of basis functions.
c  ibvtyp = 0 if no bit vector is to be stored.
c         .ne.0 if a bit vector is to be stored in the record buffer.
c
c  notes: (1) ifmt, lrec, and nmax may eventually have meaning on input
c         in future versions of this routine.  for now the input values
c         are ignored.
c         (2) for extensibility, the input variables may be passed
c         into this routine as array elements, instead of scalars, in
c         future versions of this routine.
c
c  output:
c  ifmt = ifmt parameter.
c  lrec = actual buffer length to be written.
c  nmax = maximum number of values in each record.
c  ierr = error return code.
c       =  0 for normal return.
c       = -1 for itype error.
c       = -2 for nbft error.
c       = -3 for lrec error.
c
c  08-oct-90 (columbus day) 1-e fcore change. chunk added. -rls
c  05-jul-89 ibvtyp added. -rls
c  01-aug-89  written by ron shepard.
c
      implicit none
c
      integer itype, lrecal, nbft, ibvtyp, ifmt, lrec, nmax, ierr
c
      integer nbig, chunk, nword, ibvfac, space
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # nbftmx = largest nbft consistent with label packing methods.
c     # lrecmn = minimum allowed record length.
c     # lrecmx = maximum allowed record length.  this should be
c     #          consistent with the record "dword" encoding.
c     # lrcinc = record length block size increments.
c     # ndeflt = default number of block increments. such that
c     #          the default record length = (lrecmn + ndeflt * lrcinc)
c
      integer    nbftmx
      parameter( nbftmx=2**16-1 )
      integer    lrecmn,     lrecmx,         lrcinc,      ndeflt
      parameter( lrecmn=512, lrecmx=2**16,   lrcinc=512,  ndeflt=3 )
c
c     # check the input parameters...
c
      ierr = 0
c
      if ( itype .ne. 1    .and.    itype .ne. 2 ) then
         ierr = -1
         return
      elseif ( nbft .le. 0    .or.    nbft .gt. nbftmx ) then
         ierr = -2
         return
      endif
c
c     # adjust lrec for efficient I/O.
c     # machine-dependent quirks, magic record lengths, etc. should
c     # be localized here and in the above parameter definitions..
c
      if ( lrecal .eq. -1 ) then
c        # return the default values for everything.
         lrec = lrecmn + (ndeflt * lrcinc)
      else
c        # want:     lrec = lrecmn + (n * lrcinc) <= lrecal
c        #           for the largest possible n.
         lrec = min( lrecmx, lrecal )
         lrec = lrecmn + ( ((lrec - lrecmn) / lrcinc ) * lrcinc )
      endif
      lrec = min( lrecmx, lrec )
      if ( lrec .lt. lrecmn ) then
         ierr = -3
         return
      endif
c
c     # for this record length, compute the number of values.
c
c     # 1-e records: dword // values // fcore // packed_labels // ibitv
c     # 2-e records: dword // values //          packed_labels // ibitv
c
      if ( nbft .lt. 2**8 ) then
c        # use standard 8-bit orbital label packing.
         ifmt=0
      else
c        # use standard 16-bit orbital label packing.
         ifmt=1
      endif
c
c     # compute a reasonable upper bound to nmax.
c     # chunk is used to constrain nmax so that sifulab8() or sifulab16() do
c     # not overwrite the labels(*) array.
c     #
c     # another choice would be to require the programmer to account
c     # for these overruns when allocating labels(*).  this would
c     # result in more efficient buffer(*) use, but the chunk choice
c     # is simpler.
c
c     # note that the programmer must explicitly account for ibitv(*)
c     # overruns; constraining chunk to 64 is too wasteful. -rls
c
      if ( itype .eq. 1 ) then
         if ( ifmt .eq. 0 ) then
c           # n + n/4 <= lrec
            nbig  = ( 4 * lrec ) / 5
            chunk = 4
         elseif ( ifmt .eq. 1 ) then
c           # n + n/2 <= lrec
            nbig  = ( 2 * lrec ) / 3
            chunk = 2
         endif
      elseif ( itype .eq. 2 ) then
         if ( ifmt .eq. 0 ) then
c           # n + n/2 <= lrec
            nbig  = ( 2 * lrec ) / 3
            chunk = 2
         elseif ( ifmt .eq. 1 ) then
c           # n + n <= lrec
            nbig  = lrec / 2
            chunk = 1
         endif
      endif
c
c     # round up to a higher multiple of chunk.
      nbig = ( ( nbig + 2*chunk ) / chunk ) * chunk
c
c     # account for dword and fcore.
      if ( itype .eq. 1 ) then
         nword = 2
      else
         nword = 1
      endif
c
c     # packed_bit_vector_space = ibvfac * ((n+63)/64)
      if ( ibvtyp .ne. 0 ) then
         ibvfac = 1
      else
         ibvfac = 0
      endif
c
c     # loop backwards from the upper bound until a valid nmax is found.
      do 10 nmax = nbig, 0, -chunk
         space = nword + nmax + ( ( nmax + chunk - 1) / chunk ) +
     &    ibvfac * ( ( nmax + 63 ) / 64 )
         if ( space .le. lrec ) goto 11
10    continue
11    continue
      if ( nmax .le. 0 ) then
         ierr = -4
         return
      endif
c
c     # ifmt, lrec, and nmax are all ok.
c
      return
      end
*deck sifd1
      subroutine sifd1(
     & info,   nipv,    iretbv,  buffer,
     & num,    last,    itypea,  itypeb,
     & ifmt,   ibvtyp,  values,  labels,
     & fcore,  ibitv,   ierr )
c
c  decode a 1-e buffer.
c  buffer has the form:
c    dword // packed_values(*) // fcore // packed_labels(*) //
c          //packed_bit_vector(*)
c
c  input:
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
c           note: if ifmt=0, then as many as ((nipv*n1max+7)/8)*8
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if iretbv.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c  ierr = error return code. 0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  nipv,   iretbv, num,    itypea, itypeb,
     & last,   ifmt,   ibvtyp, ierr
      real*8 buffer(*), values(*), fcore
      integer info(*), labels(*), ibitv(*)
c
      integer  l1rec,  n1max,  lab1,   lenpl,  lenbv,
     & l1recx, nuw,    bv1
      integer  unpack(4)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
c
      l1rec = info(2)
      n1max = info(3)
c
c     # unpack dword...
c
      call sifulab16( buffer(1), unpack, 4 )
      num    = unpack(1)
      lab1   = unpack(2)
      last   = mod(unpack(4),         4)
      ifmt   = mod(unpack(4)/2**2,    8)
      itypeb = mod(unpack(3),      1024)
      itypea = mod(unpack(3)/2**10,   8)
      ibvtyp = mod(unpack(3)/2**13,   8)
c
c     # if nipv=0 then only dword is unpacked...
c
      if ( nipv .eq. 0 ) return
c
      if ( ifmt .eq. 0 ) then
         lenpl=(num+3)/4
      elseif ( ifmt .eq. 1 ) then
         lenpl=(num+1)/2
      else
c        # illegal ifmt.
         ierr = -1
         return
      endif
      if ( ibvtyp .ne. 0 ) then
         lenbv=(n1max+63)/64
      else
         lenbv=0
      endif
      l1recx=(2+num+lenpl+lenbv)
      if ( l1recx .gt. l1rec ) then
c        # inconsistent l1rec.
         ierr = -2
         return
      endif
c
c     # unpack/copy the values(*)...
c
      call dcopy( num, buffer(2), 1,  values, 1 )
c
      fcore = buffer( num + 2 )
c
c     # unpack the labels(*)...
c
      if ( ifmt .eq. 0 ) then
c        # 8-bit packing of orbital labels.
         if ( nipv .eq. 1 ) then
c           # 1 integer/value output.
            nuw=num
            call sifulab16( buffer(lab1), labels, nuw )
         elseif ( nipv .eq. 2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call sifulab8( buffer(lab1), labels, nuw )
         else
c           # illegal nipv.
            ierr = -3
            return
         endif
      elseif ( ifmt .eq. 1 ) then
c        # 16-bit packing of orbital labels.
         if ( nipv .eq. 1 ) then
c           # 1 integer/value output.
            nuw=num
            call sifulab32( buffer(lab1), labels, nuw )
         elseif ( nipv .eq. 2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call sifulab16( buffer(lab1), labels, nuw )
         else
c           # illegal nipv.
            ierr = -3
            return
         endif
      endif
c
      if ( iretbv .eq. 0 ) then
c        # ignore bit-vector processing.
         continue
      elseif ( (iretbv .eq. ibvtyp) .or.
     &    ( (iretbv .eq. -1) .and. (ibvtyp .ne. 0 ) ) ) then
c        # unpack the bit vector from the end of buffer(*).
         bv1=l1rec+1-lenbv
         nuw=num
         call sifulab1( buffer(bv1), ibitv, nuw )
      elseif ( iretbv .eq. -1 ) then
c        # general bitvector request with ibvtyp=0.  not an error.
         continue
      else
c        # inconsistent ibvtyp.
         ierr = -4
         return
      endif
c
      return
      end
*deck sifd2
      subroutine sifd2(
     & info,    nipv,    iretbv,  buffer,
     & num,     last,    itypea,  itypeb,
     & ifmt,    ibvtyp,  values,  labels,
     & ibitv,   ierr )
c
c  decode a 2-e buffer.
c  buffer has the form:
c    dword // packed_values(*) // packed_labels(*) //
c          //packed_bit_vector(*)
c
c  input:
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
c  buffer(1:l2rec) = packed  buffer.
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
c  ierr = error return code. 0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  nipv,   iretbv, num,    itypea, itypeb,
     & last,   ifmt,   ibvtyp, ierr
      real*8   buffer(*),      values(*)
      integer  info(*),        labels(*),      ibitv(*)
c
      integer  l2rec,  n2max,  lab1,   lenpl,  lenbv,
     & l2recx, nuw,    bv1
      integer  unpack(4)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l2rec = info(4)
      n2max = info(5)
c
c     # unpack dword...
c
      call sifulab16( buffer(1), unpack, 4 )
      num    = unpack(1)
      lab1   = unpack(2)
      last   = mod(unpack(4),         4)
      ifmt   = mod(unpack(4)/2**2,    8)
      itypeb = mod(unpack(3),      1024)
      itypea = mod(unpack(3)/2**10,   8)
      ibvtyp = mod(unpack(3)/2**13,   8)
c
c     # if nipv=0 then only dword is unpacked...
c
      if ( nipv .eq. 0 ) return
c
      if ( ifmt.eq.0 ) then
         lenpl=(num+1)/2
      elseif ( ifmt.eq.1 ) then
         lenpl=num
      else
c        # illegal ifmt.
         ierr = -1
         return
      endif
      if ( ibvtyp.ne.0 ) then
         lenbv=(n2max+63)/64
      else
         lenbv=0
      endif
      l2recx=(1+num+lenpl+lenbv)
      if ( l2recx .gt. l2rec ) then
c        # inconsistent l2recx.
         ierr = -1
         return
      endif
c
c     # unpack/copy the values(*)...
c
      call dcopy( num, buffer(2), 1,   values, 1 )
c
c     # unpack the labels(*)...
c
      if ( ifmt.eq.0 ) then
c        # 8-bit packing of orbital labels.
         if ( nipv.eq.1 ) then
c           # 1 integer/value output.
            nuw=num
            call sifulab32( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call sifulab16( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.4 ) then
c           # 4 integers/value output.
            nuw=4*num
            call sifulab8( buffer(lab1), labels, nuw )
         else
c           # illegal nipv.
            ierr = -3
            return
         endif
      elseif ( ifmt.eq.1 ) then
c        # 16-bit packing of orbital labels.
         if ( nipv.eq.1 ) then
c           # 1 integer/value output.
c           # not allowed on 32-bit integer machines.
            ierr = -3
            return
         elseif ( nipv.eq.2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call sifulab32( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.4 ) then
c           # 4 integers/value output.
            nuw=4*num
            call sifulab16( buffer(lab1), labels, nuw )
         else
c           # illeval nipv.
            ierr = -3
            return
         endif
      endif
c
      if ( iretbv.eq.0 ) then
c        # ignore bit-vector processing.
         continue
      elseif ( (iretbv .eq. ibvtyp) .or.
     &    ( (iretbv .eq. -1) .and. (ibvtyp .ne. 0 ) ) ) then
c        # unpack the bit vector from the end of buffer(*).
         bv1=l2rec+1-lenbv
         nuw=num
         call sifulab1( buffer(bv1), ibitv, nuw )
      elseif ( iretbv .eq. -1 ) then
c        # general bitvector request with ibvtyp=0.  not an error.
         continue
      else
c        # inconsistent ibvtyp.
         ierr = -4
         return
      endif
c
      return
      end
*deck sife1
      subroutine sife1(
     & info,    nipv,    num,     last,
     & itypea,  itypeb,  ifmt,    ibvtyp,
     & values,  labels,  fcore,   ibitv,
     & buffer,  ierr )
c
c  encode a 1-e buffer.
c  buffer has the form:
c    dword // packed_values(*) // fcore // packed_labels(*) //
c          // packed_bit_vector(*)
c
c  input:
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 two orbital labels are packed in each labels(*) entry.
c       = 2 one orbital label is stored in each labels(*) entry.
c  num = number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n1max+7)/8)*8
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c
c  output:
c  num = number of values(*) and labels(*) remaining.  the calling
c        program should not assume that this is zero on return. use:
c                   numtot = numtot + num
c                   call (...num...)
c                   numtot = numtot - num
c        in the calling program to compute correctly the total
c        number of output values.  -rls
c  values(1:num) = elements that were not written on this call.
c  labels(1:nipv,1:num) = corresponding unwritten labels.
c  ibitv(1:num) = corresponding unwritten bit-vector elements.
c  buffer(1:l1rec) = packed  buffer.
c  ierr = error return code. 0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  16-aug-89 num=0 short return added. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  num,    nipv,   itypea, itypeb,
     & last,   ifmt,   ibvtyp, ierr
      real*8   buffer(*),      values(*),      fcore
      integer  info(*),        labels(nipv,*), ibitv(*)
c
      integer  l1rec,  n1max,  lenpl,  lenbv,
     & l1recx, lab1,   nuw,    bv1
      integer  unpack(4)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l1rec = info(2)
      n1max = info(3)
c
c     # check parameters for consistency...
c
      if ( num.gt.n1max ) then
         call sifbummer('sife1: num=',num,wrnerr)
         ierr = -1
         return
      elseif ( itypea.lt.0 .or. itypea.gt.7 ) then
         call sifbummer('sife1: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif ( itypeb.lt.0 .or. itypeb.gt.1023 ) then
         call sifbummer('sife1: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif ( last.lt.0 .or. last.gt.3 ) then
         call sifbummer('sife1: last=',last,wrnerr)
         ierr = -1
         return
      elseif ( ibvtyp.lt.0 .or. ibvtyp.gt.7 ) then
         call sifbummer('sife1: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if ( ifmt.eq.0 ) then
         lenpl=(num+3)/4
      elseif ( ifmt.eq.1 ) then
         lenpl=(num+1)/2
      else
         call sifbummer('sife1: ifmt=',ifmt,wrnerr)
         ierr = -1
         return
      endif
      if ( ibvtyp.ne.0 ) then
         lenbv=(n1max+63)/64
      else
         lenbv=0
      endif
      l1recx=(2+num+lenpl+lenbv)
      if ( l1recx .gt. l1rec ) then
         call sifbummer('sife1: (l1recx-l1rec)=',(l1recx-l1rec),wrnerr)
         ierr = -1
         return
      endif
c
      lab1=num+3
c
c     # pack dword...
c
      unpack(1) = num
      unpack(2) = lab1
      unpack(3) = (ibvtyp*8+itypea)*1024+itypeb
      unpack(4) = ifmt*4+last
      call sifplab16( buffer(1), unpack, 4 )
c
c     # if num=0, then don't bother with the packing.
c
      if ( num.eq.0)return
c
c     # pack/copy the values(*)...
c
      call dcopy( num, values, 1,   buffer(2), 1 )
c
      buffer( num + 2 ) = fcore
c
c     # pack the labels(*)...
c
      if ( ifmt.eq.0 ) then
c        # 8-bit packing of orbital labels.
         if ( nipv.eq.1 ) then
c           # 1 integer/value input.
            nuw=num
            call sifplab16( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.2 ) then
c           # 2 integers/value input.
            nuw=2*num
            call sifplab8( buffer(lab1), labels, nuw )
         else
            call sifbummer('sife1: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      elseif ( ifmt.eq.1 ) then
c        # 16-bit packing of orbital labels.
         if ( nipv.eq.1 ) then
c           # 1 integer/value input.
            nuw=num
            call sifplab32( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.2 ) then
c           # 2 integers/value input.
            nuw=2*num
            call sifplab16( buffer(lab1), labels, nuw )
         else
            call sifbummer('sife1: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      endif
      if ( ibvtyp.ne.0 ) then
c        # pack the bit vector at the end of buffer(*).
         bv1=l1rec+1-lenbv
         nuw=num
         call sifplab1( buffer(bv1), ibitv, nuw )
      endif
c
      return
      end
*deck sife2
      subroutine sife2(
     & info,    nipv,    num,     last,
     & itypea,  itypeb,  ifmt,    ibvtyp,
     & values,  labels,  ibitv,   buffer,
     & ierr )
c
c  encode a 2-e buffer.
c  buffer has the form:
c    dword // packed_values(*) // packed_labels(*) //
c          // packed_bit_vector(*)
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
c
c  output:
c  num = number of values(*) and labels(*) remaining.  the calling
c        program should not assume that this is zero on return. use:
c                   numtot = numtot + num
c                   call (...num...)
c                   numtot = numtot - num
c        in the calling program to compute correctly the total
c        number of output values.  -rls
c  values(1:num) = elements that were not written on this call.
c  labels(1:nipv,1:num) = corresponding unwritten labels.
c  ibitv(1:num) = corresponding unwritten bit-vector elements.
c  buffer(1:l2rec) = packed  buffer.
c  ierr = error return code. 0 for normal return.
c
c  16-aug-89 num=0 short return added.
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  nipv,   num,    itypea, itypeb, last,
     & ifmt,   ibvtyp, ierr
      real*8   buffer(*),      values(*)
      integer  info(*),        labels(nipv,*), ibitv(*)
c
      integer  n2max,  l2rec,  lab1,   lenpl,  lenbv,
     & l2recx, nuw,    bv1
      integer unpack(4)
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
         call sifbummer('sife2: num=',num,wrnerr)
         ierr = -1
         return
      elseif(itypea.lt.0 .or. itypea.gt.7)then
         call sifbummer('sife2: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif(itypeb.lt.0 .or. itypeb.gt.1023)then
         call sifbummer('sife2: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif(last.lt.0 .or. last.gt.3)then
         call sifbummer('sife2: last=',last,wrnerr)
         ierr = -1
         return
      elseif(ibvtyp.lt.0 .or. ibvtyp.gt.7)then
         call sifbummer('sife2: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if(ifmt.eq.0)then
         lenpl=(num+1)/2
      elseif(ifmt.eq.1)then
         lenpl=num
      else
         call sifbummer('sife2: ifmt=',ifmt,wrnerr)
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
         call sifbummer('sife2: (l2recx-l2rec)=',(l2recx-l2rec),wrnerr)
         ierr = -1
         return
      endif
c
      lab1=num+2
c
c     # pack dword...
c
      unpack(1) = num
      unpack(2) = lab1
      unpack(3) = (ibvtyp*8+itypea)*1024+itypeb
      unpack(4) = ifmt*4+last
      call sifplab16( buffer(1), unpack, 4 )
c
c     # if num=0, then don't bother with the packing.
c
      if(num.eq.0)return
c
c     # pack/copy the values(*)...
c
      call dcopy( num, values, 1,   buffer(2), 1 )
c
c     # pack the labels(*)...
c
      if(ifmt.eq.0)then
c        # 8-bit packing of orbital labels.
         if(nipv.eq.1)then
c           # 1 integer/value input.
            nuw=num
            call sifplab32( buffer(lab1), labels, nuw )
         elseif(nipv.eq.2)then
c           # 2 integers/value input.
            nuw=2*num
            call sifplab16( buffer(lab1), labels, nuw )
         elseif(nipv.eq.4)then
c           # 4 integers/value input.
            nuw=4*num
            call sifplab8( buffer(lab1), labels, nuw )
         else
            call sifbummer('sife2: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      elseif(ifmt.eq.1)then
c        # 16-bit packing of orbital labels.
         if(nipv.eq.1)then
c           # 1 integer/value input.
c           # *** not allowed on 32-bit integer machines. ***
            call sifbummer('sif2e: ifmt=1, nipv=',nipv,wrnerr)
            ierr = -1
            return
         elseif(nipv.eq.2)then
c           # 2 integers/value input.
            nuw=2*num
            call sifplab32( buffer(lab1), labels, nuw )
         elseif(nipv.eq.4)then
c           # 4 integers/value input.
            nuw=4*num
            call sifplab16( buffer(lab1), labels, nuw )
         else
            call sifbummer('sife2: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      endif
      if(ibvtyp.ne.0)then
c        # pack the bit vector at the end of buffer(*).
         bv1=l2rec+1-lenbv
         nuw=num
         call sifplab1( buffer(bv1), ibitv, nuw )
      endif
c
      return
      end
*deck sifew1
      subroutine sifew1(
     & aoints,  info,    nipv,    num,
     & last,    itypea,  itypeb,  ifmt,
     & ibvtyp,  values,  labels,  fcore,
     & ibitv,   buffer,  nrec,    ierr )
c
c  encode and write a 1-e integral record.
c
c  input:
c  aoints = output file unit number.
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 two orbital labels are packed in each labels(*) entry.
c       = 2 one orbital label is stored in each labels(*) entry.
c  num = actual number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n1max+7)/8)*8
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c
c  output:
c  num = reset to the number of unwritten elements in values(*).
c        if (last.ne.0) then num is always zero on return. otherwise,
c        the calling program should not assume that this is set to zero.
c        note: this provision is to allow future data-dependent
c              value(*) and labels(*) packing methods.  use:
c                   numtot = numtot + num
c                   call (...num...)
c                   numtot = numtot - num
c              in the calling program to compute correctly the total
c              number of output values.  -rls
c  values(1:num) = elements that were not written on this call.
c  labels(1:nipv,1:num) = corresponding unwritten labels.
c  ibitv(1:num) = corresponding unwritten bit-vector elements.
c  buffer(1:l1rec) = packed  buffer.
c  nrec = updated record count.
c         note: the calling program should not assume that this is
c               always incremented by 1.
c  ierr = error return code.  0 for normal return.
c
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, nipv,   num,    last,
     & itypea, itypeb, ifmt,   ibvtyp, nrec,   ierr
      integer  info(*),        labels(*),      ibitv(*)
      real*8   values(*),      buffer(*),      fcore
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
c     # pack the buffer...
c
      call sife1(
     & info,   nipv,   num,   last,
     & itypea, itypeb, ifmt,  ibvtyp,
     & values, labels, fcore, ibitv,
     & buffer, ierr )
      if ( ierr .ne. 0 ) return
c
c     # write to the output file...
c
      call sifseqwbf( aoints, buffer, l1rec )
c     # sifseqwbf() does not return ierr.
      ierr = 0
c
c     # update nrec, reset num, and move any unwritten values..
c
      nrec = nrec + 1
      num  = 0
c
      return
      end
*deck sifew2
      subroutine sifew2(
     & aoint2,  info,    nipv,    num,
     & last,    itypea,  itypeb,  ifmt,
     & ibvtyp,  values,  labels,  ibitv,
     & buffer,  iwait,   nrec,    reqnum,
     & ierr )
c
c  encode and write a 2-e integral record.
c
c  input:
c  aoint2 = output file unit number.
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 four orbital labels are packed in each labels(*) entry.
c       = 2 two orbital labels are packed in each labels(*) entry.
c       = 4 one orbital label is stored in each labels(*) entry.
c  num = actual number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n2max+7)/8)*8
c                 elements of labels(*) are referenced.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c             note: as many as ((n1max+63)/64)*64 elements of this
c                   array are referenced.
c  iwait   = asynchronous i/o parameter.
c          = 0  don't wait.  use asynch i/o and return without
c               waiting for i/o completion.  the calling program
c               must call sif2w8() before reusing the buffer or
c               before calling this routine again with the same
c               output buffer.
c          = 1  wait for i/o completion before returning.
c               no sif2w8() call is required in the calling program.
c               buffer(*) can be reused immediately upon return.
c
c  output:
c  num = reset to the number of unwritten elements in values(*).
c        if (last.ne.0) then num is always zero on return. otherwise,
c        the calling program should not assume that this is set to zero.
c        note: this provision is to allow future data-dependent
c              value(*) and labels(*) packing methods. use:
c                   numtot = numtot + num
c                   call (...num...)
c                   numtot = numtot - num
c        in the calling program to compute correctly the total
c        number of output values.  -rls
c  values(1:num) = elements that were not written on this call.
c  labels(1:nipv,1:num) = corresponding unwritten labels.
c  ibitv(1:num) = corresponding unwritten bit-vector elements.
c  buffer(1:l2rec) = packed  buffer.
c  nrec = updated record count.
c         note: the calling program should not assume that this is
c               always incremented by 1.
c  reqnum = i/o request number for the i/o operation associated with
c           this record.
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoint2, nipv,   num,    last,   itypea, itypeb,
     & ifmt,   ibvtyp, iwait,  nrec,   reqnum, ierr
      integer  info(*),        labels(*),      ibitv(*)
      real*8   values(*),      buffer(*)
c
      integer  l2rec,  n2max
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l2rec = info(4)
      n2max = info(5)
c
c     # pack the buffer...
c
      call sife2(
     & info,   nipv,   num,   last,
     & itypea, itypeb, ifmt,  ibvtyp,
     & values, labels, ibitv, buffer,
     & ierr )
      if ( ierr .ne. 0 ) return
c
c     # write to the output file.
c
      call sifw2( aoint2, iwait, info, buffer, reqnum, ierr )
      if ( ierr .ne. 0 ) return
c
c     # update nrec, reset num, and move any uncopied values..
c
      nrec = nrec + 1
      num  = 0
c
      return
      end
*deck siffr1
      subroutine siffr1(
     & ninput,  info,    nipv,    num,
     & last,    itypea,  itypeb,  ifmt,
     & ibvtyp,  values,  labels,  fcore,
     & ibitv,   ierr )
c
c  read a formatted 1-e integral record.
c  (see also routine siffw1().)
c
c  input:
c  ninput = input file unit number.
c  info(*) = info array for this file.
c
c  output:
c  nipv = number of integers per value returned.
c       = 1 return two orbital labels packed in each labels(*) entry.
c       = 2 return one orbital label in each labels(*) entry.
c  num = actual number of values in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format of the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = values.
c  labels(1:nipv,1:num) = integral labels
c  fcore = frozen core contribution.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c  ierr = return code. 0 for normal return.
c
c  24-apr-92 nipv added to the input record. -rls
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  ninput, nipv,   num,   last,
     & itypea, itypeb, ifmt,   ibvtyp, ierr
      integer  info(*),        labels(*),     ibitv(*)
      real*8   values(*),      fcore
c
      integer  i, ioff, lab1
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
c
c     # read in the "dword" information.
c
      read(ninput,*,iostat=ierr)
     & num,    lab1,   ibvtyp, itypea,
     & itypeb, ifmt,   last,   nipv
      if ( ierr .ne. 0 ) then
         call sifbummer(' siffr1() dword, ierr=', ierr, wrnerr )
         return
      endif
c
      read(ninput,*,iostat=ierr) fcore
      if ( ierr .ne. 0 ) then
         call sifbummer(' siffr1() fcore, ierr=', ierr, wrnerr )
         return
      endif
c
c     # read the values and labels.
c
      if ( nipv .eq. 1 ) then
         do 10 i = 1, num
            read(ninput,*,iostat=ierr)
     &       values(i), labels(i)
            if ( ierr .ne. 0 ) then
               call sifbummer(' siffr1() nipv=1, ierr=', ierr, wrnerr )
               return
            endif
10       continue
      elseif ( nipv .eq. 2 ) then
         ioff = 0
         do 20 i = 1, num
            read(ninput,*,iostat=ierr)
     &       values(i), labels(ioff+1), labels(ioff+2)
            if ( ierr .ne. 0 ) then
               call sifbummer(' siffr1() nipv=2, ierr=', ierr, wrnerr )
               return
            endif
            ioff = ioff + 2
20       continue
      else
         ierr = -10
         return
      endif
c
c     # read the bit vector if present.
c
      if ( ibvtyp .ne. 0 ) then
         read(ninput,*,iostat=ierr) ( ibitv(i), i = 1, num )
         if ( ierr .ne. 0 ) then
            call sifbummer(' siffr1() ibitv, ierr=', ierr, wrnerr )
            return
         endif
      endif
c
      return
      end
*deck siffr2
      subroutine siffr2(
     & ninput,  info,    nipv,    num,
     & last,    itypea,  itypeb,  ifmt,
     & ibvtyp,  values,  labels,  ibitv,
     & ierr )
c
c  read a formatted 2-e integral record.
c  (see also routine siffw2().)
c
c  input:
c  ninput = input file unit number.
c  info(*) = info array for this file.
c
c  output:
c  nipv = number of integers per value to be returned.
c       = 1 return four orbital labels packed in each labels(*) entry.
c       = 2 return two orbital labels packed in each labels(*) entry.
c       = 4 return one orbital label in each labels(*) entry.
c  num = actual number of values in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format of the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = values.
c  labels(1:nipv,1:num) = integral labels.
c  ibitv(*) = unpacked bit vector (referenced only if ibvtyp.ne.0).
c  ierr = error return code.  0 for normal return.
c
c  24-apr-92 nipv added to the input record. -rls
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  ninput, nipv,   num,    last,
     & itypea, itypeb, ifmt,   ibvtyp, ierr
      integer  info(*),        labels(*),      ibitv(*)
      real*8   values(*)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer i, ioff, lab1
c
      ierr = 0
c
c     # read the "dword" information.
c
      read(ninput,*,iostat=ierr)
     & num,    lab1,   ibvtyp, itypea,
     & itypeb, ifmt,   last,   nipv
      if ( ierr .ne. 0 ) then
         call sifbummer('siffr2(), dword ierr=', ierr, wrnerr )
         return
      endif
c
      if ( nipv .eq. 1 ) then
         do 10 i = 1, num
            read(ninput,*,iostat=ierr) values(i), labels(i)
            if ( ierr .ne. 0 ) then
               call sifbummer('siffr2(), nipv=1, ierr=', ierr, wrnerr )
               return
            endif
10       continue
      elseif ( nipv .eq. 2 ) then
         ioff = 0
         do 20 i = 1, num
            read(ninput,*,iostat=ierr)
     &       values(i), labels(ioff+1), labels(ioff+2)
            if ( ierr .ne. 0 ) then
               call sifbummer('siffr2(), nipv=2, ierr=', ierr, wrnerr )
               return
            endif
            ioff = ioff + 2
20       continue
      elseif ( nipv .eq. 4 ) then
         ioff = 0
         do 30 i = 1, num
            read(ninput,*,iostat=ierr)
     &       values(i),
     &       labels(ioff+1), labels(ioff+2),
     &       labels(ioff+3), labels(ioff+4)
            if ( ierr .ne. 0 ) then
               call sifbummer('siffr2(), nipv=4, ierr=', ierr, wrnerr )
               return
            endif
            ioff = ioff + 4
30       continue
      else
         ierr = -10
         return
      endif
c
c     # read the bit vector if present.
c
      if ( ibvtyp .ne. 0 ) then
         read(ninput,*,iostat=ierr) ( ibitv(i), i = 1, num )
         if ( ierr .ne. 0 ) then
            call sifbummer('siffr2(), ibitv ierr=', ierr, wrnerr )
            return
         endif
      endif
c
      return
      end
*deck siffw1
      subroutine siffw1(
     & info,    nipv,    num,     last,
     & itypea,  itypeb,  ifmt,    ibvtyp,
     & values,  labels,  fcore,   ibitv,
     & nlist,   ierr )
c
c  formatted write of a 1-e buffer.
c
c  input:
c  info(*) = info array for this file.
c  nipv = number of integers per value.
c       = 1 two orbital labels are packed in each labels(*) entry.
c       = 2 one orbital label is packed in each labels(*) entry.
c  num = number of values to place in the packed buffer.
c  last = integral continuation parameter.
c  itypea,itypeb = generic and specific integral types.
c  ifmt = format to use for the packed buffer.
c  ibvtyp = type of packed bit-vector.
c  values(1:num) = integral values.
c  labels(1:nipv,1:num) = integral labels
c           note: if ifmt=0, then as many as ((nipv*n2max+3)/4)*4
c                 elements of labels(*) are referenced.
c  fcore = frozen core contribution.
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
      integer  nlist,  nipv,   num,    itypea, itypeb,
     & last,   ifmt,   ibvtyp, ierr
      real*8   values(*),      fcore
      integer  info(*),        labels(nipv,*), ibitv(*)
c
      integer  l1rec,  n1max,  lenpl,  lenbv,  l1recx,
     & lab1,   ifmtv,  i,      j
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      ierr  = 0
      l1rec = info(2)
      n1max = info(3)
c
c     # check parameters for consistency...
c
      if(num.gt.n1max)then
         call sifbummer('siffw1: num=',num,wrnerr)
         ierr = -1
         return
      elseif(itypea.lt.0 .or. itypea.gt.7)then
         call sifbummer('siffw1: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif(itypeb.lt.0 .or. itypeb.gt.1023)then
         call sifbummer('siffw1: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif(last.lt.0 .or. last.gt.3)then
         call sifbummer('siffw1: last=',last,wrnerr)
         ierr = -1
         return
      elseif(ibvtyp.lt.0 .or. ibvtyp.gt.7)then
         call sifbummer('siffw1: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if(ifmt.eq.0)then
         lenpl=(num+3)/4
      elseif(ifmt.eq.1)then
         lenpl=(num+1)/2
      else
         call sifbummer('siffw1: ifmt=',ifmt,wrnerr)
         ierr = -1
         return
      endif
      if(ibvtyp.ne.0)then
         lenbv=(n1max+63)/64
      else
         lenbv=0
      endif
      l1recx=(2+num+lenpl+lenbv)
      if( l1recx .gt. l1rec )then
         call sifbummer('siffw1: (l1recx-l1rec)=',(l1recx-l1rec),wrnerr)
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
      else
         call sifbummer('siffw1: nipv=',nipv,wrnerr)
         ierr = -1
         return
      endif
c
      write(nlist,6021) fcore
c
      do 10 i = 1, num
         write(nlist,ifmtv) values(i), (labels(j,i), j=1,nipv)
10    continue
6021  format(1x,1pe20.12, i7  )
6022  format(1x,1pe20.12, 2i4 )
c
      if ( ibvtyp .ne. 0 ) then
         write(nlist,6030) ( ibitv(i), i = 1, num )
      endif
6030  format(1x,20i2)
c
      return
      end
