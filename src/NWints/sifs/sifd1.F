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
      call ulab16( buffer(1), unpack, 4 )
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
      call cdcopy( num, buffer(2), 1,  values, 1 )
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
            call ulab16( buffer(lab1), labels, nuw )
         elseif ( nipv .eq. 2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call ulab8( buffer(lab1), labels, nuw )
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
            call ulab32( buffer(lab1), labels, nuw )
         elseif ( nipv .eq. 2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call ulab16( buffer(lab1), labels, nuw )
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
         call ulab1( buffer(bv1), ibitv, nuw )
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
