      subroutine esp_matinv(a,n,is,iad1,iad2,d,mdm)
      implicit real*8 (a-h,o-z)
c
c     inversion of square matrix a by means of the gauss-jordan
c     algorithm
c
c     april 72/rs9b
c
      real*8 a(mdm,mdm),d(mdm)
      integer is(2,mdm),iad1(mdm),iad2(mdm)
c
      data zero/0.0d0/, one/1.0d0/, small/1.0d-20/
c
 2000 format(' warning from inv: matrix is singular')
c
      do 1 l=1,n
      is(1,l)=0
  1   is(2,l)=0
      do 9 ima=1,n
      b=zero
      do 2 l=1,n
      do 2 m=1,n
      if(is(1,l).eq.1.or.is(2,m).eq.1) goto 2
      e=dabs(a(l,m))
      if(e.lt.b) goto 8
      i=l
      k=m
    8 b=dmax1(b,e)
  2   continue
      is(1,i)=1
      is(2,k)=1
      iad1(k)=i
      iad2(i)=k
      b=a(i,k)
c
c     pivot
c
      if(dabs(b).lt. small) goto 20
      a(i,k)=one/b
      do 6 l=1,n
      if(l.eq.k) goto 6
c
c     kellerzeile
c
      a(i,l)=-a(i,l)/b
  6   continue
      do 5 l=1,n
      do 5 m=1,n
      if(l.eq.i.or.m.eq.k) goto 5
c
c     rechteck-regel
c
      a(l,m)=a(l,m)+a(l,k)*a(i,m)
  5   continue
      do 11 l=1,n
      if(l.eq.i) goto 11
c
c     pivot-spalte
c
      a(l,k)=a(l,k)/b
  11  continue
  9   continue
c
c     permutation der zeilen, um die natuerliche ordnung wieder herzuste
c
      do 15 l=1,n
      do 13 j=1,n
      k=iad1(j)
  13  d(j)=a(k,l)
      do 14 j=1,n
  14  a(j,l)=d(j)
  15  continue
c
c     permutation der spalten
c
      do 16 l=1,n
      do 17 j=1,n
      k=iad2(j)
  17  d(j)=a(l,k)
      do 18 j=1,n
  18  a(l,j)=d(j)
  16  continue
      return
c
c     error exit: matrix is singular
c
  20  write(*,2000)
      stop 'inv in polar'
      end
