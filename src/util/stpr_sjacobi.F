      subroutine stpr_sjacobi(nm,n,a,w,eivr,ierr)
c $Id: stpr_sjacobi.F,v 1.1 2000-08-28 23:27:41 edo Exp $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      dimension a(n*(n+1)/2),w(n),eivr(nm,n)
c
c     use jacobi method to emulate eispack rs
c     this routine uses a variable threshold jacobi method
c     it gives very good eigenvalues and eigenvectors
c     the routine is much faster than the old hdiag routine written
c     at m.i.t. that uses the jacobi method but not the variable
c     threshold technique that is applied here
c
c-----parameters---------------
      sq2inv = 1.d0/sqrt(2.d0)
      t1     = 1.d-12
      t2     = 1.d-12
c------------------------------
      ierr   = 0
      avgf   = dble(n*(n-1))*0.55d0
      if(n.lt.1) then
        goto 160
      else if(n.lt.2) then
        eivr(1,1)=1.0d0
        goto 160
      end if
      do 20 j=1,n
        do 10 i=1,n
          eivr(i,j)=0.0d0
   10   continue
        eivr(j,j)=1.0d0
   20 continue
c
c     Find the absolutely largest element of a.
c
      atop=0.0d0
      ij = 0
      do 40 i=1,n
        do 30 j=1,i
          ij = ij + 1
          if(atop.lt.dabs(a(ij)))atop = dabs(a(ij))
   30   continue
   40 continue
      if(atop.le.0.0d0) then
        ierr = 1
        goto 160
      end if
c
c     Calculate the stopping criterion -- dstop.
c
      d = 0.0d0
      ij = 0
      do 60 i=1,n
        do 50 j=1,i-1
          ij = ij + 1
          d = d + a(ij)**2
   50   continue
        ij = ij + 1
   60 continue
      dstop=t1*d
c
c     Calculate the threshold, thrsh.
c
      thrsh= sqrt(d/avgf)
c
c     Start a sweep.
c
   70 continue
      iflag=0
      do 150 jcol=2,n
        jcoltr = jcol*(jcol-1)/2
        jcol1=jcol-1
        do 140 irow=1,jcol1
          indxa = jcoltr + irow
          aij=a(indxa)
c
c         Compare the off-diagonal element with thrsh.
c
          daij = dabs(aij)
          if(daij.le.thrsh) goto 140
          irowtr = irow*(irow-1)/2
          indxa = irowtr + irow
          aii=a(indxa)
          indxa = jcoltr + jcol
          ajj=a(indxa)
          s=ajj-aii
          ds = dabs(s)
c
c         The chosen rotation is less than the rounding.
c         Do not rotate.
c
          if (daij.lt.t2*ds) goto 140
          iflag=1
          if(t2*daij.ge.ds)then
c
c           Rotation is very close to 45 degrees,
c           sin and cos = 1/(root 2).
c
            s = sq2inv
            c = s
          else
c
c           Rotation is not very close to 45 degrees.
c
            t = aij/s
            u = 0.25d0/sqrt(0.25d0+t*t)
            c = sqrt(0.5d0+u)
            s = 2.d0*t*u/c
          end if
c
c         Calculate new elements of matrix a.
c
          do 80 i=1,irow
            t         = a(irowtr + i)
            u         = a(jcoltr + i)
            a(irowtr + i) = c*t-s*u
            a(jcoltr + i) = s*t+c*u
   80     continue
          i2 = irow+2
          if (i2.le.jcol) then
            do 90 i=i2,jcol
              itr = (i-1)*(i-2)/2
              t           = a(jcoltr + i-1)
              u           = a(itr + irow)
              a(jcoltr + i-1) = s*u+c*t
              a(itr + irow) = c*u-s*t
   90         continue
          end if
          a(jcoltr + jcol) = s*aij+c*ajj
          a(irowtr + irow) = c*a(irowtr + irow)-s*(c*aij-s*ajj)
          do 100 j=jcol,n
            jtr = j*(j-1)/2
            t         = a(jtr + irow)
            u         = a(jtr + jcol)
            a(jtr + irow) = c*t-s*u
            a(jtr + jcol) = s*t+c*u
  100     continue
c
c         Rotation completed. see if eigenvectors are wanted by
c         user.
c
          do 110 i=1,n
            t=eivr(i,irow)
            eivr(i,irow)=c*t-eivr(i,jcol)*s
            eivr(i,jcol)=s*t+eivr(i,jcol)*c
  110     continue
c
c         Calculate the new norm d and compare with dstop.
c
          s=aij
          d=d-s*s
          if(d.lt.dstop) then
c
c           Recalculate dstop and thrsh to discard rounding errors.
c
            d=0.0d0
            ij = 0
            do 130 i=1,n
              do 120 j=1,i-1
                ij = ij + 1
                d = d + a(ij)**2
  120         continue
              ij = ij + 1
  130       continue
            dstop=t1*d
          end if
          thrsh=sqrt(d/avgf)
  140   continue
  150 continue
      if(iflag.ne.0) goto 70
  160 continue
c
c     Fill eigenvalue vector.
c
      idiag = 0
      do 170 i=1,n
        idiag = idiag + i
        w(i) = a(idiag)
  170 continue
c
c     Arrange eigenvalues & vectors in ascending order.
c
      isort=1
  180 continue
      if(isort.eq.1)then
        isort=0
        do 200 i=1,n-1
          if(w(i).gt.w(i+1))then
            temp=w(i)
            w(i)=w(i+1)
            w(i+1)=temp
            do 190 k = 1,n
              temp=eivr(k,i)
              eivr(k,i)=eivr(k,i+1)
              eivr(k,i+1)=temp
  190       continue
            isort=1
          endif
  200   continue
        goto 180
      endif
      return
      end
