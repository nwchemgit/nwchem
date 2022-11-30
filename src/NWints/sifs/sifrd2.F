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
