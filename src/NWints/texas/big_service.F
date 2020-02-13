c $Id$
c
c     this module contains the home-made service routines
c     like matrix manipulations etc.
c     some of these duplicate blas routines and should be replaced
c     by them later
c.... now come the homebrewn 
c..  replaced by blas routine (same calling sequence)
c     subroutine add1 (a,con,b,nn)
c     implicit real*8 (a-h,o-z)
c     dimension a(1), b(1)
c     n=nn
c     do 10 i=1,n
c  10 b(i)=b(i)+con*a(i)
c     return
c
c     end
c     subroutine tfer (a,b,n)
c     implicit real*8 (a-h,o-z)
c     dimension a(1), b(1)
c
c           transfer a to b
c
c     nn=n
c     do 10 i=1,nn
c  10 b(i)=a(i)
c     return
c
c     end
      subroutine txs_add(a,b,n)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
c...  adds a to b, result in b
c%include '/sys/ins/base.ins.ftn'
c%include '/sys/ins/vec.ins.ftn'
c      call vec_$dadd_vector(a,b,n,b)
      do 100 i=1,n
        b(i)=b(i)+a(i)
 100  continue
      end
      subroutine txs_add1(a,con,b,n)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      call daxpy(n,con,a,1,b,1)
      end
c
      subroutine tfer (a,b,n)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
c
*      if (n.lt.1000) then
         do i = 1, n
            b(i) = a(i)
         enddo
*      else
*         call dcopy(n,a,1,b,1)
*      endif
c
      end
      subroutine mxmtr(a,b,c,m,n,p)
c...  this is not a true blas routine- it is our own creation
c...  this subroutine performs the matrix multiply c=a(trspse)*b
c     both the storage dimension and the real dimension of the
c     matrices is given by a(p,m),b(p,n),c(m,n)
c     it uses an efficient blocking algorithm
c     note that ingeneral it is more efficient to perform a(tr)b
c     than a*b. transpose the matrix before if you have to
      implicit real*8 (a-h,o-z)
      integer p
      dimension  a(p,m), b(p,n),c(m,n)
      parameter (ibl=16)
c     performs a(tr)*b=c
      do 500 i=1,m,ibl
        iend=min0(i+ibl-1,m)
        do 500 j=1,n,ibl
          jend=min0(j+ibl-1,n)
          do 500 i1=i,iend 
            do 500 j1=j,jend 
          sum=0.0d0
c         sum=vec_$ddot(a(1,i1),b(1,j1),p)
          sum=ddot(p,a(1,i1),1,b(1,j1),1)
c         call scalar(a(1,i1),b(1,j1),p,sum)
c         do 400 k=1,p 
c            sum=sum+a(k,i1)*b(k,j1)
c400  continue
       c(i1,j1)=sum
 500  continue 
      end
c
c     subroutine sdiag2v(m,n,a,d,x)
      subroutine sdiag2 (m,n,a,d,x)
      implicit real*8 (a-h,o-z)
      integer*4 f4,h4,l4
ckwol parameter (mxdim=401)
ckwol parameter (mxdim=601)
      parameter (mxdim=2000)
c
c
c      computation of all eigenvalues and eigenvectors of a real
c      symmetric matrix by the method of qr transformations.
c      if the euclidean norm of the rows varies   s t r o n g l y
c      most accurate results may be obtained by permuting rows and
c      columns to give an arrangement with increasing norms of rows.
c
c      two machine constants must be adjusted appropriately,
c      eps = minimum of all x such that 1+x is greater than 1 on the
c      e     computer,
c      tol = inf / eps  with inf = minimum of all positive x represen-
c            table within the computer.
c      a dimension statement e(256) may also be changed appropriately.
c
c      input
c
c      (m)   not larger than mxdim,  corresponding value of the actual
c            dimension statement a(m,m), d(m), x(m,m),
c      (n)   not larger than (m), order of the matrix,
c      (a)   the matrix to be diagonalized, its lower triangle has to
c            be given as  ((a(i,j), j=1,i), i=1,n),
c
c      output
c
c      (d)   components d(1), ..., d(n) hold the computed eigenvalues
c            in ascending sequence. the remaining components of (d) are
c            unchanged,
c      (x)   the computed eigenvector corresponding to the j-th eigen-
c            value is stored as column (x(i,j), i=1,n). the eigenvectors
c            are normalized and orthogonal to working accuracy. the
c            remaining entries of (x) are unchanged.
c
c      array (a) is unaltered. however, the actual parameters
c      corresponding to (a) and (x)  may be identical, ''overwriting''
c      the eigenvectors on (a).
c
c      leibniz-rechenzentrum, munich 1965
c
c
      dimension a(m,m), d(m), x(m,m)
      dimension e(mxdim)
c
c     correct adjustment for ieee 64  bit floating  point numbers`
c       is about (not too sensitive)
      eps=1.0d-14
      tol=1.0d-200
c
      if (n.eq.1) go to 350
      if (n.gt.mxdim) then
         write(*,*) ' dimension too big, increase mxdim in sdiag2'
         write (*,*) ' max= ',mxdim, ' actual= ',n
          call txs_error
c         can be replaced by go to 360 if you do not have "error"
      end if
      do 10 i=1,n
      do 10 j=1,i
   10 x(i,j)=a(i,j)
c
c     householder reduction
c     simulation of loop do 150 i=n,2,(-1)
c
      do 150 ni=2,n
         ii=n+2-ni
c
c     fake loop for recursive address calculation
c
      do 150 i=ii,ii
         l=i-2
         h=0.0
         g=x(i,i-1)
#if 1
         l4=l
         if (l4) 140,140,20
#else
         if (l) 140,140,20
#endif
   20    do 30 k=1,l
   30    h=h+x(i,k)**2
         s=h+g*g
         if (s.ge.tol) go to 40
         h=0.0
         go to 140
#if 1
   40    h4=h
         if (h4) 140,140,50
#else
   40    if (h) 140,140,50
#endif
   50    l=l+1
         f=g
         g=dsqrt(s)
#if 1
         f4=f
         if (f4) 70,70,60
#else
         if (f) 70,70,60
#endif
   60    g=-g
   70    h=s-f*g
         x(i,i-1)=f-g
         f=0.0
c
         do 110 j=1,l
            x(j,i)=x(i,j)/h
            s=0.0
c           do 80 k=1,j
c  80       s=s+x(j,k)*x(i,k)
c           s=vec_$ddot_i(x(j,1),m,x(i,1),m,j)
            s=ddot(j,x(j,1),m,x(i,1),m)
            j1=j+1
            if (j1.gt.l) go to 100
c           do 90 k=j1,l
c  90       s=s+x(k,j)*x(i,k)
            j1l=l-j1+1
c           write(*,*) s,j,j1,i,m,j1l
c           s=s+vec_$ddot_i(x(j1,j),1,x(i,j1),m,j1l)
            s=s+ddot(j1l,x(j1,j),1,x(i,j1),m)
  100       e(j)=s/h
  110    f=f+s*x(j,i)
c
         f=f/(h+h)
c
c        do 120 j=1,l
c 120    e(j)=e(j)-f*x(i,j)
c        call vec_$dmult_add_i(e,1,x(i,1),m,l,-f,e,1)
         call daxpy(l,-f,x(i,1),m,e,1)
c
         do 130 j=1,l
            f=x(i,j)
            s=e(j)
c        do 130 k=1,j
c 130    x(j,k)=x(j,k)-f*e(k)-x(i,k)*s
c        call vec_$dmult_add_i(x(j,1),m,e,1,j,-f,x(j,1),m)
c        call vec_$dmult_add_i(x(j,1),m,x(i,1),m,j,-s,x(j,1),m)
         call daxpy(j,-f,e,1,x(j,1),m)
         call daxpy(j,-s,x(i,1),m,x(j,1),m)
  130    continue
c
  140    d(i)=h
  150 e(i-1)=g
c
c     accumulation of transformation matrices
c
      d(1)=x(1,1)
      x(1,1)=1.0
      do 200 i=2,n
         l=i-1
         if (d(i)) 190,190,160
  160    do 180 j=1,l
            s=0.0
c           do 170 k=1,l
c 170       s=s+x(i,k)*x(k,j)
c           s=vec_$ddot_i(x(i,1),m,x(1,j),1,l)
            s=ddot(l,x(i,1),m,x(1,j),1)
c        do 180 k=1,l
c 180    x(k,j)=x(k,j)-s*x(k,i)
c      call vec_$dmult_add(x(1,j),x(1,i),l,-s,x(1,j))
       call daxpy(l,-s,x(1,i),1,x(1,j),1)
  180  continue
  190    d(i)=x(i,i)
         x(i,i)=1.0
      do 200 j=1,l
         x(i,j)=0.0
  200 x(j,i)=0.0
c     call vec_$dzero_i(x(i,1),m,l)
c     call vec_$dzero(x(1,i),l)
c 200 continue
c
c     diagonalization of the tridiagonal matrix
c
      b=0.0
      f=0.0
      e(n)=0.0
c
      do 310 l=1,n
         h=eps*(dabs(d(l))+dabs(e(l)))
         if (h.gt.b) b=h
c
c     test for splitting
c
         do 210 j=l,n
            if (dabs(e(j)).le.b) go to 220
  210    continue
c
c     test for convergence
c
  220    if (j.eq.l) go to 310
c
c     shift from upper 2*2 minor
c
  230    p=(d(l+1)-d(l))*0.5/e(l)
         r=dsqrt(p*p+1.0)
         if (p) 240,250,250
  240    p=p-r
         go to 260
  250    p=p+r
  260    h=d(l)-e(l)/p
         do 270 i=l,n
  270    d(i)=d(i)-h
         f=f+h
c
c     qr transformation
c
         p=d(j)
         c=1.0
         s=0.0
c
c     simulation of loop do 330 i=j-1,l,(-1)
c
         j1=j-1
         do 300 ni=l,j1
            ii=l+j1-ni
c
c     fake loop for recursive address calculation
c
         do 300 i=ii,ii
            g=c*e(i)
            h=c*p
c
c     protection against underflow of exponents
c
            if (dabs(p).lt.dabs(e(i))) go to 280
            c=e(i)/p
            r=dsqrt(c*c+1.0)
            e(i+1)=s*p*r
            s=c/r
            c=1.0/r
            go to 290
  280       c=p/e(i)
            r=dsqrt(c*c+1.0)
            e(i+1)=s*e(i)*r
            s=1.0/r
            c=c/r
  290       p=c*d(i)-s*g
            d(i+1)=h+s*(c*g+s*d(i))
c        do 300 k=1,n
c           h=x(k,i+1)
c           x(k,i+1)=x(k,i)*s+h*c
c 300    x(k,i)=x(k,i)*c-h*s
c        call vec_$dcopy(x(1,i),tmp1,n)
c        call vec_$dcopy(x(1,i+1),tmp2,n)
c        call vec_$dmult_constant(tmp2,n,c,x(1,i+1))
c        call vec_$dmult_add(x(1,i+1),x(1,i),n,s,x(1,i+1))
c        call vec_$dmult_constant(tmp2,n,-s,x(1,i))
c        call vec_$dmult_add(x(1,i),tmp1,n,c,x(1,i))
         call drot(n,x(1,i+1),1,x(1,i),1,c,s)
  300    continue
c
         e(l)=s*p
         d(l)=c*p
         if (abs(e(l)).gt.b) go to 230
c
c     convergence
c
  310 d(l)=d(l)+f
c
c     ordering of eigenvalues
c
      ni=n-1
      do 340 i=1,ni
         k=i
         p=d(i)
         j1=i+1
         do 320 j=j1,n
            if (d(j).ge.p) go to 320
            k=j
            p=d(j)
  320    continue
         if (k.eq.i) go to 340
         d(k)=d(i)
         d(i)=p
         do 330 j=1,n
            p=x(j,i)
            x(j,i)=x(j,k)
  330    x(j,k)=p
  340 continue
      go to 360
c
c     special treatment of case n = 1
c
  350 d(1)=a(1,1)
      x(1,1)=1.0
  360 return
c
      end
c
c     subroutine sdiag2 (m,n,a,d,x)
      subroutine sdiag2sc (m,n,a,d,x)
      implicit real*8 (a-h,o-z)
c
c
c      computation of all eigenvalues and eigenvectors of a real
c      symmetric matrix by the method of qr transformations.
c      if the euclidean norm of the rows varies   s t r o n g l y
c      most accurate results may be obtained by permuting rows and
c      columns to give an arrangement with increasing norms of rows.
c
c      two machine constants must be adjusted appropriately,
c      eps = minimum of all x such that 1+x is greater than 1 on the
c      e     computer,
c      tol = inf / eps  with inf = minimum of all positive x represen-
c            table within the computer.
c      a dimension statement e(160) may also be changed appropriately.
c
c      input
c
c      (m)   not larger than 160,  corresponding value of the actual
c            dimension statement a(m,m), d(m), x(m,m),
c      (n)   not larger than (m), order of the matrix,
c      (a)   the matrix to be diagonalized, its lower triangle has to
c            be given as  ((a(i,j), j=1,i), i=1,n),
c
c      output
c
c      (d)   components d(1), ..., d(n) hold the computed eigenvalues
c            in ascending sequence. the remaining components of (d) are
c            unchanged,
c      (x)   the computed eigenvector corresponding to the j-th eigen-
c            value is stored as column (x(i,j), i=1,n). the eigenvectors
c            are normalized and orthogonal to working accuracy. the
c            remaining entries of (x) are unchanged.
c
c      array (a) is unaltered. however, the actual parameters
c      corresponding to (a) and (x)  may be identical, ''overwriting''
c      the eigenvectors on (a).
c
c      leibniz-rechenzentrum, munich 1965
c
c
      dimension a(m,m), d(m), x(m,m)
      dimension e(160)
      integer*4 f4,h4,l4
c
c     correct adjustment for ieee floating point numbers (64 bits)
c
      eps=2.0d-15
      tol=1.0d-292
c
      if (n.eq.1) go to 350
      do 10 i=1,n
      do 10 j=1,i
   10 x(i,j)=a(i,j)
c
c     householder reduction
c     simulation of loop do 150 i=n,2,(-1)
c
      do 150 ni=2,n
         ii=n+2-ni
c
c     fake loop for recursive address calculation
c
      do 150 i=ii,ii
         l=i-2
         h=0.0d0
         g=x(i,i-1)
#if 1
         l4=l
         if (l4) 140,140,20
#else
         if (l) 140,140,20
#endif
   20    do 30 k=1,l
   30    h=h+x(i,k)**2
         s=h+g*g
         if (s.ge.tol) go to 40
         h=0.0d0
         go to 140
#if 1
   40    h4=h
         if (h4) 140,140,50
#else
   40    if (h) 140,140,50
#endif
   50    l=l+1
         f=g
         g=dsqrt(s)
#if 1
         f4=f
         if (f4) 70,70,60
#else
         if (f) 70,70,60
#endif
   60    g=-g
   70    h=s-f*g
         x(i,i-1)=f-g
         f=0.0d0
c
         do 110 j=1,l
            x(j,i)=x(i,j)/h
            s=0.0d0
            do 80 k=1,j
   80       s=s+x(j,k)*x(i,k)
            j1=j+1
            if (j1.gt.l) go to 100
            do 90 k=j1,l
   90       s=s+x(k,j)*x(i,k)
  100       e(j)=s/h
  110    f=f+s*x(j,i)
c
         f=f/(h+h)
c
         do 120 j=1,l
  120    e(j)=e(j)-f*x(i,j)
c
         do 130 j=1,l
            f=x(i,j)
            s=e(j)
         do 130 k=1,j
  130    x(j,k)=x(j,k)-f*e(k)-x(i,k)*s
c
  140    d(i)=h
  150 e(i-1)=g
c
c     accumulation of transformation matrices
c
      d(1)=x(1,1)
      x(1,1)=1.0d0
      do 200 i=2,n
         l=i-1
         if (d(i)) 190,190,160
  160    do 180 j=1,l
            s=0.0d0
            do 170 k=1,l
  170       s=s+x(i,k)*x(k,j)
         do 180 k=1,l
  180    x(k,j)=x(k,j)-s*x(k,i)
  190    d(i)=x(i,i)
         x(i,i)=1.0d0
      do 200 j=1,l
         x(i,j)=0.0d0
  200 x(j,i)=0.0d0
c
c     diagonalization of the tridiagonal matrix
c
      b=0.0d0
      f=0.0d0
      e(n)=0.0d0
c
      do 310 l=1,n
         h=eps*(abs(d(l))+abs(e(l)))
         if (h.gt.b) b=h
c
c     test for splitting
c
         do 210 j=l,n
            if (abs(e(j)).le.b) go to 220
  210    continue
c
c     test for convergence
c
  220    if (j.eq.l) go to 310
c
c     shift from upper 2*2 minor
c
  230    p=(d(l+1)-d(l))*0.5d0/e(l)
         r=dsqrt(p*p+1.0d0)
         if (p) 240,250,250
  240    p=p-r
         go to 260
  250    p=p+r
  260    h=d(l)-e(l)/p
         do 270 i=l,n
  270    d(i)=d(i)-h
         f=f+h
c
c     qr transformation
c
         p=d(j)
         c=1.0d0
         s=0.0d0
c
c     simulation of loop do 330 i=j-1,l,(-1)
c
         j1=j-1
         do 300 ni=l,j1
            ii=l+j1-ni
c
c     fake loop for recursive address calculation
c
         do 300 i=ii,ii
            g=c*e(i)
            h=c*p
c
c     protection against underflow of exponents
c
            if (abs(p).lt.abs(e(i))) go to 280
            c=e(i)/p
            r=dsqrt(c*c+1.0d0)
            e(i+1)=s*p*r
            s=c/r
            c=1.0d0/r
            go to 290
  280       c=p/e(i)
            r=dsqrt(c*c+1.0d0)
            e(i+1)=s*e(i)*r
            s=1.0d0/r
            c=c/r
  290       p=c*d(i)-s*g
            d(i+1)=h+s*(c*g+s*d(i))
         do 300 k=1,n
            h=x(k,i+1)
            x(k,i+1)=x(k,i)*s+h*c
  300    x(k,i)=x(k,i)*c-h*s
c
         e(l)=s*p
         d(l)=c*p
         if (abs(e(l)).gt.b) go to 230
c
c     convergence
c
  310 d(l)=d(l)+f
c
c     ordering of eigenvalues
c
      ni=n-1
      do 340 i=1,ni
         k=i
         p=d(i)
         j1=i+1
         do 320 j=j1,n
            if (d(j).ge.p) go to 320
            k=j
            p=d(j)
  320    continue
         if (k.eq.i) go to 340
         d(k)=d(i)
         d(i)=p
         do 330 j=1,n
            p=x(j,i)
            x(j,i)=x(j,k)
  330    x(j,k)=p
  340 continue
      go to 360
c
c     special treatment of case n = 1
c
  350 d(1)=a(1,1)
      x(1,1)=1.0d0
  360 return
c
      end
c
      subroutine zeroit (a,n)
      implicit real*8 (a-h,o-z)
      dimension a(*)
c     call vec_$dzero (a,n)
      do  100 i=1,n
        a(i)=0.0d0
 100  continue
      end
      subroutine mult(tmat,con,num)
      implicit real*8 (a-h,o-z)
      dimension tmat(1)
c
c     this subroutine multiplies an array by a constant
c     this is jmb invention; it is the same as dscal
      call  dscal(num,con,tmat,1)
c
      end
      subroutine mave(amat,vold,vnew,num)
      implicit real*8 (a-h,o-z)
      dimension amat(1),vold(1),vnew(1)
c
c     this subroutine multiplies a triangular matrix by a column vector
c     and returns a new column vector.  vnew cannot be same place as vold
c
      n=num
      do 13 i=1,n
      ij=i*(i-1)/2
      sum=0.0
         do 20 j=1,n
         ij=ij+1
         if(j.gt.i) ij=ij+j-2
         sum=sum+amat(ij)*vold(j)
   20 continue
      vnew(i)=sum
   13 continue
      return
c
      end
      subroutine inner(x,y,sum,ne)
c   WHO ON EARTH WROTE THIS - THIS DUPLICATES SCALAR, AND USES
C   SINGLE PRECISION ZERO
      implicit real*8 (a-h,o-z)
      dimension x(1),y(1)
c     multiplies the elements in two arrays with the same index and adds
c     them similar to inner product
c
c     sum=0.0
      sum=0.0d0
      m=ne
      do 30 l=1,m
      sum=sum+x(l)*y(l)
   30 continue
      return
c
      end
      subroutine tri(a,b,m)
      implicit real*8 (a-h,o-z)
      dimension a(1), b(1)
c
c     a symmetric quadratic matrix a is changed to triangular form b
c
      mm=m
      ij=0
      do 200 i=1,mm
         iad=i*mm-mm
         do 100 j=1,i
            ij=ij+1
            b(ij)=a(iad+j)
  100 continue
  200 continue
      return
c
      end
      subroutine quadbfta(a,b,c,m)
      implicit real*8 (a-h,o-z)
      dimension a(1),b(1)
c
c     make a quadratic symmetric or antisymmetric matrix b from
c     triangular matrix a
c
      c1=c
      con=abs(c1)
      ij=0
      do 10 i=1,m
         iad=i*m-m
         do 20 j=1,i
            jad=j*m-m
            ij=ij+1
            b(iad+j)=con*a(ij)
            b(jad+i)=c1*a(ij)
   20 continue
      b(iad+i)=b(iad+i)*(c1+con)/2
   10 continue
      return
c
      end
      subroutine vecmat(amat,istcol,icol,vecin,vecout,idime)
      implicit real*8 (a-h,o-z)
      dimension amat(1), vecin(1), vecout(1)
c
c     multiplies amat(t) times a vector starting from stcol to stcol+ind
c     looks like josep bofill creation
      ist=istcol
      ic=icol
      id=idime
      indcol=ist+ic-1
      do 10 i=ist,indcol
         nad=i*id-id
         sum=0.0
         do 20 j=1,id
            sum=sum+amat(j+nad)*vecin(j)
   20 continue
         vecout(i)=sum
   10 continue
      return
c
      end
      subroutine matmat(a,ias,ia,b,ibs,ib,d,id)
      implicit real*8 (a-h,o-z)
      dimension a(1),b(1),d(1)
c
c     multiplies a(t) x b with common dimension id, result in d
c
      ja=ia
      jas=ias
      jb=ib
      jbs=ibs
      jae=ja+jas-1
      jbe=jb+jbs-1
      idim=id
      do 30 i=jas,jae
         iad=i*idim-idim
         do 20 j=jbs,jbe
            jad=j*idim-idim
            sum=0.0
            do 10 k=1,idim
               sum=sum+a(k+iad)*b(k+jad)
   10 continue
            d(jad+i)=sum
   20 continue
   30 continue
      return
c
      end
      subroutine eig (b,a,d,n)
      implicit real*8 (a-h,o-z)
      dimension a(n,n), b(1), d(1)
      ij=0
      do 10 i=1,n
      do 10 j=1,i
         ij=ij+1
         a(i,j)=b(ij)
   10 continue
      call sdiag2 (n,n,a,d,a)
      return
c
c
      end
      subroutine spur (a,b,ncf,s)
      implicit real*8 (a-h,o-z)
      dimension a(1), b(1)
c
c     ....    spur(ab) - both are synmmetrical triangular matrices
c
cccc  n=ncf
      s=zero
      ij=0
      do 10 i=1,ncf
      do 10 j=1,i
         ij=ij+1
         s=s+a(ij)*b(ij)
         if (i.eq.j) s=s-a(ij)*b(ij)*0.5d0
   10 continue
      s=2*s
      return
c
      end
      subroutine txs_trans (a,ncf)
      implicit real*8 (a-h,o-z)
      dimension a(ncf,ncf)
      do 10 i=1,ncf
      do 10 j=1,i
         s=a(i,j)
         a(i,j)=a(j,i)
         a(j,i)=s
   10 continue
      return
c
      end
c======================================================================
c linking routine between mamu and mxma which translates existing calls
c of mamu directly into calls of mxma.
c richard g.a. bone, march, 1991.
c main interpretative problem is that mamu performs b*a matrix multiplication
c whereas mxma performs a*b.
c mamu works with either square or triangular matrices according to
c arguments k,l,m, whereas mxma requires squares only.
c n.b. in pulay group mxma may not take the same array address twice.
c
c     subroutine mamu (b,a,c,k,l,m,n,w)
c     implicit double precision (a-h,o-z)
c     dimension a(1),b(1),c(1),w(1)
c     common /big/bl(1)
c     common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,
c    1 npl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),
c    2 icode(200),inpf,ioutf
c     nsq = n**2
c     ntri = n*(n+1)/2
c
c     if((k.ne.0).and.(k.ne.1)) goto 9998
c     if((l.ne.0).and.(l.ne.1)) goto 9998
c     if((m.ne.0).and.(m.ne.1)) goto 9998
c
c     call getmem(nsq,ia)
c     call getmem(nsq,ic)
c
c     if (l.eq.0) call squr(a,bl(ia),n)
c     if (l.eq.1) call tfer(a,bl(ia),nsq)
c
c     if (k.eq.0) then
c      call getmem(nsq,ib)
c      call squr(b,bl(ib),n)
c      call mxma(bl(ia),1,n,bl(ib),1,n,bl(ic),1,n,n,n,n)
c      call retmem(1)
c     endif
c     if (k.eq.1) then
c      ir=ir-1
c      call mxma(bl(ia),1,n,b,1,n,bl(ic),1,n,n,n,n)
c     endif
c
c     if (m.eq.0) call tri(bl(ic),c,n)
c     if (m.eq.1) call tfer(bl(ic),c,nsq)
c
c     call retmem(ir)
c
c     call retmem(2)
c     return
c9998 write(icond,*)
c    1  'incorrect use of mamu: integer arguments must be 1 or 0'
c     return
c     end
c
c======================================================================
      subroutine squr(t,sq,n)
      implicit real*8(a-h,o-z)
c     triangle to square
      dimension sq(1),t(1)
      ij=0
      ii=0
      do 20 i=1,n
      jj=0
c$dir no_recurrence
cdir$ ivdep
      do 10 j=1,i
      ij=ij+1
      sq(ii+j)=t(ij)
      sq(jj+i)=t(ij)
10    jj=jj+n
20    ii=ii+n
      return
      end
      subroutine mxma(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
c     this code should run at about 34 mflops for matrices
c     of dimensions between 50 and 500 on a 25 mhz rs6000
      implicit real*8 (a-h,o-z)
      dimension  a(*),b(*),r(*)
      parameter (nb=60)
c     written based on tips by  ron bell, ibm united kingdom, ltd.
c     seeibm document no. gg24-3611
c     note that this code, in spite of all efforts, does not reach
c     the numerical performance (43 mflops/s) claimed in the above
c     publication. it is assumed now that the performance figures
c     of bell refer to the inner loop only and do not reflect
c     the overhead.  the maximum performance of this code is about
c     35 mflops/s
      dimension rr(nb,nb),aa(nb,nb),bb(nb,nb)
c     performs r=a*b; the actual dimensions are naxma,maxmb,naxmb
c     the addressing function is: r(i,j)=r((i-1)*mcolr+(j-1)*mrowr+1)
c     a(i,k)=a((i-1)*mcola+(k-1)*mrowa+1)
c     b(k,j)=b((k-1)*mcolb+(j-1)*mrowb+1)
c     good example of obscure addressing in fortran
c     please do not write code like this. unfortunately, this piece of
c     code is widely used on the crays
c     the buffer size has been roughly optimized for the apollo dn10k
      ir=1
      do 3 j=1,nrow
        irr=ir
        do 2 i=1,ncol
c        ncol is the number of rows. the dutch mind is perplexing
          r(irr)=0.0d0
          irr=irr+mcolr
 2      continue
        ir=ir+mrowr
 3    continue
      do 1400 ii=1,ncol,nb
          ie=min0(ii+nb-1,ncol)
          ie1=ie-ii+1
          ie2=3*(ie1/3)
c         write(*,*) ie,ie1
          do 1300 jj=1,nrow,nb
            je=min0(jj+nb-1,nrow)
            je1=je-jj+1
            je2=3*(je1/3)
c      make sure that ie2 and je2 are divisible by 3 for the
c      loop unrolling
            ir=(ii-1)*mcolr+(jj-1)*mrowr+1
            do 200 i=1,ie1
              irr=ir
              do 100 j=1,je1
c               rr(i,j)=r(irr)
                rr(i,j)=0.0d0
                irr=irr+mrowr
 100          continue
              ir=ir+mcolr
 200        continue
            do 1000 kk=1,nlink,nb
              ke=min0(kk+nb-1,nlink)
              ke1=ke-kk+1
              ia=(ii-1)*mcola+(kk-1)*mrowa+1
              do 400 i=1,ie1
                iaa=ia
                do 300 k=1,ke1
                  aa(k,i)=a(iaa)
                  iaa=iaa+mrowa
 300            continue
                ia=ia+mcola
 400          continue
              ib=(kk-1)*mcolb+(jj-1)*mrowb+1
              do 600 j=1,je1
                ibb=ib
                do 500 k=1,ke1
                  bb(k,j)=b(ibb)
                  ibb=ibb+mcolb
 500            continue
                ib=ib+mrowb
 600          continue
              do 900 i=1,ie2,3
                do 800 j=1,je2,3
                  s00=rr(i,j)
                  s10=rr(i+1,j)
                  s20=rr(i+2,j)
                  s01=rr(i,j+1)
                  s11=rr(i+1,j+1)
                  s21=rr(i+2,j+1)
                  s02=rr(i,j+2)
                  s12=rr(i+1,j+2)
                  s22=rr(i+2,j+2)
                  do 700 k=1,ke1
                    s00=s00+aa(k,i)*bb(k,j)
                    s01=s01+aa(k,i)*bb(k,j+1)
                    s02=s02+aa(k,i)*bb(k,j+2)
                    s10=s10+aa(k,i+1)*bb(k,j)
                    s11=s11+aa(k,i+1)*bb(k,j+1)
                    s12=s12+aa(k,i+1)*bb(k,j+2)
                    s20=s20+aa(k,i+2)*bb(k,j)
                    s21=s21+aa(k,i+2)*bb(k,j+1)
                    s22=s22+aa(k,i+2)*bb(k,j+2)
 700              continue
                  rr(i,j)=s00
                  rr(i+1,j)=s10
                  rr(i+2,j)=s20
                  rr(i,j+1)=s01
                  rr(i+1,j+1)=s11
                  rr(i+2,j+1)=s21
                  rr(i,j+2)=s02
                  rr(i+1,j+2)=s12
                  rr(i+2,j+2)=s22
 800            continue
                do 820 j=je2+1,je1
                  s00=rr(i,j)
                  s10=rr(i+1,j)
                  s20=rr(i+2,j)
                  do 810 k=1,ke1
                    s00=s00+aa(k,i)*bb(k,j)
                    s10=s10+aa(k,i+1)*bb(k,j)
                    s20=s20+aa(k,i+2)*bb(k,j)
 810              continue
                  rr(i,j)=s00
                  rr(i+1,j)=s10
                  rr(i+2,j)=s20
 820            continue
 900          continue
              do 950 i=ie2+1,ie1
               do 940 j=1,je1
                do 930 k=1,ke1
                 rr(i,j)=rr(i,j)+aa(k,i)*bb(k,j)
 930            continue
 940           continue
 950          continue
 1000       continue
            ir=(ii-1)*mcolr+(jj-1)*mrowr+1
            do 1200 i=1,ie1
              irr=ir
              do 1100 j=1,je1
                r(irr)=r(irr)+rr(i,j)
                irr=irr+mrowr
 1100         continue
              ir=ir+mcolr
 1200       continue
 1300     continue
 1400 continue
      end
      subroutine nerror(noer,routine,message,n1,n2)
#include "errquit.fh"
      character*(*) routine
      character*(*) message
      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      write(iout,*) 'Error no. ',noer,' in ',routine,' ',message,n1,n2
      write(icond,*) 'Error no. ',noer,' in ',routine,' ',message,n1,n2
c100  format(' Error No.',i4,' in ',a20,2x,/,a80,/,' variables ',2i6) 
      call errquit('texas: nerror called', 0, INT_ERR)
*      stop 20
      end
      subroutine txs_message(routine,messag1,n1,n2)
      character*(*) routine
      character*(*) messag1
      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      write(iout,*) 'Message from ',routine,' ',messag1,n1,n2
      write(icond,*) 'Message from ',routine,' ',messag1,n1,n2
      end
      subroutine txs_error
#include "errquit.fh"
      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      write(iout,*) ' error found'
      write(icond,*) ' error found'
*      stop 10
      call errquit('texas: txs_error called', 0, INT_ERR)
      end
c==================================================================
      subroutine dxypz(n,dx,incx, dy,incy, dz,incz)
c
c vector_x*vector_y + vector_z ===> vector_z 
c
c modification of the daxpy blas routine. (KW)
c
c
      double precision dx(*),dy(*),dz(*) 
      integer i,incx,incy,incz, ix,iy,iz, m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1.and.incz.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      iz = 1
      if(incx.lt.0) ix= (-n+1)*incx + 1
      if(incy.lt.0) iy= (-n+1)*incy + 1
      if(incz.lt.0) iz =(-n+1)*incz + 1
c
      do 10 i = 1,n
        dz(iz) = dz(iz) + dy(iy)*dx(ix)
        ix = ix + incx
        iy = iy + incy
        iz = iz + incz
   10 continue
c
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dz(i) = dz(i) + dy(i)*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i0= mp1,n,4
        i1=i0+1
        i2=i1+1
        i3=i2+1
        dz(i0) = dz(i0) + dy(i0)*dx(i0)
        dz(i1) = dz(i1) + dy(i1)*dx(i1)
        dz(i2) = dz(i2) + dy(i2)*dx(i2)
        dz(i3) = dz(i3) + dy(i3)*dx(i3)
   50 continue
      return
      end
c==================================================================
      subroutine dxyz(n,dx,incx, dy,incy, dz,incz)
c
c vector_x*vector_y  ===> vector_z 
c
c modification of the daxpy blas routine. (KW)
c
c
      double precision dx(*),dy(*),dz(*) 
      integer i,incx,incy,incz, ix,iy,iz, m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1.and.incz.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      iz = 1
      if(incx.lt.0) ix= (-n+1)*incx + 1
      if(incy.lt.0) iy= (-n+1)*incy + 1
      if(incz.lt.0) iz =(-n+1)*incz + 1
c
      do 10 i = 1,n
        dz(iz) = dy(iy)*dx(ix)
        ix = ix + incx
        iy = iy + incy
        iz = iz + incz
   10 continue
c
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dz(i) = dy(i)*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i0= mp1,n,4
        i1=i0+1
        i2=i1+1
        i3=i2+1
        dz(i0) =  dy(i0)*dx(i0)
        dz(i1) =  dy(i1)*dx(i1)
        dz(i2) =  dy(i2)*dx(i2)
        dz(i3) =  dy(i3)*dx(i3)
   50 continue
      return
      end
c==================================================================
      subroutine dxpy(n,dx,incx, dy,incy)
c
c vector_x+vector_y  ===> vector_y 
c
c modification of the daxpy blas routine. (KW)
c
c
      double precision dx(*),dy(*)
      integer i,incx,incy, ix,iy, m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0) ix= (-n+1)*incx + 1
      if(incy.lt.0) iy= (-n+1)*incy + 1
c
      do 10 i = 1,n
        dy(iy) = dy(iy) + dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
c
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i0= mp1,n,4
        i1=i0+1
        i2=i1+1
        i3=i2+1
        dy(i0) =  dy(i0) + dx(i0)
        dy(i1) =  dy(i1) + dx(i1)
        dy(i2) =  dy(i2) + dx(i2)
        dy(i3) =  dy(i3) + dx(i3)
   50 continue
      return
      end
c==================================================================
      subroutine get_ij_half(ij, i, j)
c
c extracts indeces i & j from common index ij=i*(i-1)/2 + j
c
      implicit none
      integer ij, i, j
      intrinsic sqrt, float
c
      i = sqrt(float(ij + ij))
      j = ij - i*(i-1)/2
      if (i .lt. j) then
         i = i + 1
         j = ij - i*(i-1)/2
      endif
c
      end
c-----------------------------------------------------------------
      subroutine get_ij_full(ij,nj, i, j)
c
c     extracts indeces i and j from common index ij=(i-1)*nj + j
c
      implicit none
      integer ij,nj, i, j
c
      i=ij/nj+1
      j=ij-(i-1)*nj
c
      if(j.eq.0) then
         i=i-1
         j=ij-(i-1)*nj
      endif
c
      end
c-----------------------------------------------------------------
