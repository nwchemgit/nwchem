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
         call bummer('sife2: num=',num,wrnerr)
         ierr = -1
         return
      elseif(itypea.lt.0 .or. itypea.gt.7)then
         call bummer('sife2: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif(itypeb.lt.0 .or. itypeb.gt.1023)then
         call bummer('sife2: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif(last.lt.0 .or. last.gt.3)then
         call bummer('sife2: last=',last,wrnerr)
         ierr = -1
         return
      elseif(ibvtyp.lt.0 .or. ibvtyp.gt.7)then
         call bummer('sife2: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if(ifmt.eq.0)then
         lenpl=(num+1)/2
      elseif(ifmt.eq.1)then
         lenpl=num
      else
         call bummer('sife2: ifmt=',ifmt,wrnerr)
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
         call bummer('sife2: (l2recx-l2rec)=',(l2recx-l2rec),wrnerr)
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
      call plab16( buffer(1), unpack, 4 )
c
c     # if num=0, then don't bother with the packing.
c
      if(num.eq.0)return
c
c     # pack/copy the values(*)...
c
      call cdcopy( num, values, 1,   buffer(2), 1 )
c
c     # pack the labels(*)...
c
      if(ifmt.eq.0)then
c        # 8-bit packing of orbital labels.
         if(nipv.eq.1)then
c           # 1 integer/value input.
            nuw=num
            call plab32( buffer(lab1), labels, nuw )
         elseif(nipv.eq.2)then
c           # 2 integers/value input.
            nuw=2*num
            call plab16( buffer(lab1), labels, nuw )
         elseif(nipv.eq.4)then
c           # 4 integers/value input.
            nuw=4*num
            call plab8( buffer(lab1), labels, nuw )
         else
            call bummer('sife2: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      elseif(ifmt.eq.1)then
c        # 16-bit packing of orbital labels.
         if(nipv.eq.1)then
c           # 1 integer/value input.
c           # *** not allowed on 32-bit integer machines. ***
            call bummer('sif2e: ifmt=1, nipv=',nipv,wrnerr)
            ierr = -1
            return
         elseif(nipv.eq.2)then
c           # 2 integers/value input.
            nuw=2*num
            call plab32( buffer(lab1), labels, nuw )
         elseif(nipv.eq.4)then
c           # 4 integers/value input.
            nuw=4*num
            call plab16( buffer(lab1), labels, nuw )
         else
            call bummer('sife2: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      endif
      if(ibvtyp.ne.0)then
c        # pack the bit vector at the end of buffer(*).
         bv1=l2rec+1-lenbv
         nuw=num
         call plab1( buffer(bv1), ibitv, nuw )
      endif
c
      return
      end
