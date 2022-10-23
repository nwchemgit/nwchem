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
cgk debug
      write(*,*)'gk: in sife1'
cgk end
      ierr  = 0
      l1rec = info(2)
      n1max = info(3)
c
c     # check parameters for consistency...
c
      if ( num.gt.n1max ) then
         call bummer('sife1: num=',num,wrnerr)
         ierr = -1
         return
      elseif ( itypea.lt.0 .or. itypea.gt.7 ) then
         call bummer('sife1: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif ( itypeb.lt.0 .or. itypeb.gt.1023 ) then
         call bummer('sife1: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif ( last.lt.0 .or. last.gt.3 ) then
         call bummer('sife1: last=',last,wrnerr)
         ierr = -1
         return
      elseif ( ibvtyp.lt.0 .or. ibvtyp.gt.7 ) then
         call bummer('sife1: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if ( ifmt.eq.0 ) then
         lenpl=(num+3)/4
      elseif ( ifmt.eq.1 ) then
         lenpl=(num+1)/2
      else
         call bummer('sife1: ifmt=',ifmt,wrnerr)
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
         call bummer('sife1: (l1recx-l1rec)=',(l1recx-l1rec),wrnerr)
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
cgk debug
      write(*,*)'gk: calling plab16'
cgk end
      call plab16( buffer(1), unpack, 4 )
cgk debug
      write(*,*)'gk: back from plab16'
cgk end
c
c     # if num=0, then don't bother with the packing.
c
      if ( num.eq.0)return
c
c     # pack/copy the values(*)...
c
cgk debug
      write(*,*)'gk: calling cdcopy'
cgk end
      call cdcopy( num, values, 1,   buffer(2), 1 )
cgk debug
      write(*,*)'gk: back from cdcopy'
cgk end
c
      buffer( num + 2 ) = fcore
c
c     # pack the labels(*)...
c
cgk debug
      write(*,*)'gk: in sife1, ifmt = ', ifmt
cgk end
      if ( ifmt.eq.0 ) then
c        # 8-bit packing of orbital labels.
         if ( nipv.eq.1 ) then
c           # 1 integer/value input.
            nuw=num
            call plab16( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.2 ) then
c           # 2 integers/value input.
            nuw=2*num
            call plab8( buffer(lab1), labels, nuw )
         else
            call bummer('sife1: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      elseif ( ifmt.eq.1 ) then
c        # 16-bit packing of orbital labels.
         if ( nipv.eq.1 ) then
c           # 1 integer/value input.
            nuw=num
            call plab32( buffer(lab1), labels, nuw )
         elseif ( nipv.eq.2 ) then
c           # 2 integers/value input.
            nuw=2*num
            call plab16( buffer(lab1), labels, nuw )
         else
            call bummer('sife1: nipv=',nipv,wrnerr)
            ierr = -1
            return
         endif
      endif
cgk debug
      write(*,*)'gk: in sife1, ibvtyp = ', ibvtyp
cgk end
      if ( ibvtyp.ne.0 ) then
c        # pack the bit vector at the end of buffer(*).
         bv1=l1rec+1-lenbv
         nuw=num
         call plab1( buffer(bv1), ibitv, nuw )
      endif
cgk debug
      write(*,*)'gk:  leaving sife1'
cgk end
c
      return
      end
