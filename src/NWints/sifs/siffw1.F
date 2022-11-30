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
         call bummer('siffw1: num=',num,wrnerr)
         ierr = -1
         return
      elseif(itypea.lt.0 .or. itypea.gt.7)then
         call bummer('siffw1: itypea=',itypea,wrnerr)
         ierr = -1
         return
      elseif(itypeb.lt.0 .or. itypeb.gt.1023)then
         call bummer('siffw1: itypeb=',itypeb,wrnerr)
         ierr = -1
         return
      elseif(last.lt.0 .or. last.gt.3)then
         call bummer('siffw1: last=',last,wrnerr)
         ierr = -1
         return
      elseif(ibvtyp.lt.0 .or. ibvtyp.gt.7)then
         call bummer('siffw1: ibvtyp=',ibvtyp,wrnerr)
         ierr = -1
         return
      endif
c
      if(ifmt.eq.0)then
         lenpl=(num+3)/4
      elseif(ifmt.eq.1)then
         lenpl=(num+1)/2
      else
         call bummer('siffw1: ifmt=',ifmt,wrnerr)
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
         call bummer('siffw1: (l1recx-l1rec)=',(l1recx-l1rec),wrnerr)
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
         call bummer('siffw1: nipv=',nipv,wrnerr)
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
