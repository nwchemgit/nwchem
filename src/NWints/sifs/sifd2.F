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
cgk debug
!      write(*,*)'gk: in sifd2, l2rec = ', l2rec
      write(*,*)'n2max = ', n2max
      write(*,*)'nipv = ', nipv
cgk end
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
cgk debug
!      write(*,*)'num = ', num
!      write(*,*)'lab1 = ', lab1
!      write(*,*)'last = ', last
!      write(*,*)'ifmt = ', ifmt
!      write(*,*)'itypeb = ', itypeb
!      write(*,*)'itypea = ', itypea
!      write(*,*)'ibvtyp = ', ibvtyp
cgk end
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
      call cdcopy( num, buffer(2), 1,   values, 1 )
c
c     # unpack the labels(*)...
c
      if ( ifmt.eq.0 ) then
c        # 8-bit packing of orbital labels.
         if ( nipv.eq.1 ) then
c           # 1 integer/value output.
            nuw=num
            call ulab32( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.2 ) then
c           # 2 integers/value output.
            nuw=2*num
            call ulab16( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.4 ) then
c           # 4 integers/value output.
            nuw=4*num
            call ulab8( buffer(lab1), labels, nuw )
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
            call ulab32( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.4 ) then
c           # 4 integers/value output.
            nuw=4*num
            call ulab16( buffer(lab1), labels, nuw )
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
