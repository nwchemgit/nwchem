
!***** Double precision linear algebra *******************************

      subroutine loupd(x,n,nm,iperm,nperm)
! Lower-upper decomposition of matrix x, logical size n x n,
! stored as size nm x nm.  Returns matrix x(i,j),i<=j as upper
! triangular matrix and x(i,j),i>j as lower diagonal matrix with 1 on
! the diagonal.  The product of these two gives the original matrix x.
! Nperm is the order of the permutations of rows, either +1 or -1,
! unless the matrix is singular, in which case nperm is returned as 0
! and evaluation of the subroutine is stopped.
      implicit none
      integer n,nm,nperm
      integer iperm(nm)
      integer i,ii,j,jj,k,kk,im
      double precision x(nm,nm),temp
      double precision cmax,scl(nm),check
! scl(i) is scaling factor for row i, to normalize largest element to 1
      nperm=1
      do i=1,n
        cmax=0.d0
        do j=1,n
          if(abs(x(i,j)).gt.cmax) cmax=abs(x(i,j))
        enddo
        if (cmax.eq.0.d0) then
          nperm=0
          return
        endif
        scl(i)=1.d0/cmax
      enddo
      do j=1,n
        do i=1,j-1
          temp=x(i,j)
          do k=1,i-1
            temp=temp-x(i,k)*x(k,j)
          enddo
          x(i,j)=temp
        enddo
        cmax=0.d0
        do i=j,n
          temp=x(i,j)
          do k=1,j-1
            temp=temp-x(i,k)*x(k,j)
          enddo
          x(i,j)=temp
          check=scl(i)*abs(temp)
          if (check.ge.cmax) then
            im=i
            cmax=check
          endif
        enddo
        if (j.ne.im) then
          do k=1,n
            temp=x(im,k)
            x(im,k)=x(j,k)
            x(j,k)=temp
          enddo
          nperm=-nperm
          check=scl(j)
          scl(j)=scl(im)
          scl(im)=check
        endif
        iperm(j)=im
        if(j.ne.n) then
          if (x(j,j).ne.0) then
            temp=(1.d0,0.d0)/x(j,j)
            do i=j+1,n
              x(i,j)=x(i,j)*temp
            enddo
          else
            nperm=0
            return
          endif
        endif
      enddo
      return
      end

      subroutine baksublu(x,v,n,nm,iperm)
! Solves a linear equation x'*w=v for w.  X is the lower-upper
! decomposition of x' from subroutine loupd.  Vector v is replaced by
! the solutions w as output.
      implicit none
      integer n,nm,iperm(nm),i,j,ii,jj
      double precision x(nm,nm),v(nm),temp
      ii=0
      do i=1,n
        temp=v(iperm(i))
        v(iperm(i))=v(i)
        if (ii.ne.0) then
          do j=ii,i-1
            temp=temp-x(i,j)*v(j)
          enddo
        elseif (temp.ne.0.d0) then
          ii=i
        endif
        v(i)=temp
      enddo
      do i=n,1,-1
        temp=v(i)
        do j=i+1,n
          temp=temp-x(i,j)*v(j)
        enddo
        v(i)=temp/x(i,i)
      enddo
      return
      end


      subroutine inverse(x,y,n,nm)
! Returns matrix y as the inverse of matrix x, both
! of dimension n x n.
      implicit none
      integer n,nm
      double precision x(nm,nm),y(nm,nm),v(nm)
      integer i,j,k,ii,jj,kk
      integer iperm(nm),nperm
      do i=1,n
        do j=1,n
          if (i.eq.j) then
            y(i,j)=1.d0
          else
            y(i,j)=0.d0
          endif
        enddo
      enddo
      call loupd(x,n,nm,iperm,nperm)
      do j=1,n
        call baksublu(x,y(1,j),n,nm,iperm)
      enddo
      return
      end

!***** Complex linear algebra *********************************************

      subroutine loupdc(x,n,nm,iperm,nperm)
! Lower-upper decomposition of matrix x, logical size n x n,
! stored as size nm x nm.  Returns matrix x(i,j),i<=j as upper
! triangular matrix and x(i,j),i>j as lower diagonal matrix with 1 on
! the diagonal.  The product of these two gives the original matrix x.
! Nperm is the order of the permutations of rows, either +1 or -1,
! unless the matrix is singular, in which case nperm is returned as 0
! and evaluation of the subroutine is stopped.
      implicit none
      integer n,nm,nperm
      integer iperm(nm)
      integer i,ii,j,jj,k,kk,im
      double complex x(nm,nm),temp
      double precision cmax,scl(nm),check
! scl(i) is scaling factor for row i, to normalize largest element to 1
      nperm=1
      do i=1,n
        cmax=0.d0
        do j=1,n
          if(abs(x(i,j)).gt.cmax) cmax=abs(x(i,j))
        enddo
        if (cmax.eq.0.d0) then
          nperm=0
          return
        endif
        scl(i)=1.d0/cmax
      enddo
      do j=1,n
        do i=1,j-1
          temp=x(i,j)
          do k=1,i-1
            temp=temp-x(i,k)*x(k,j)
          enddo
          x(i,j)=temp
        enddo
        cmax=0.d0
        do i=j,n
          temp=x(i,j)
          do k=1,j-1
            temp=temp-x(i,k)*x(k,j)
          enddo
          x(i,j)=temp
          check=scl(i)*abs(temp)
          if (check.ge.cmax) then
            im=i
            cmax=check
          endif
        enddo
        if (j.ne.im) then
          do k=1,n
            temp=x(im,k)
            x(im,k)=x(j,k)
            x(j,k)=temp
          enddo
          nperm=-nperm
          check=scl(j)
          scl(j)=scl(im)
          scl(im)=check
        endif
        iperm(j)=im
        if(j.ne.n) then
          if (x(j,j).ne.0) then
            temp=(1.d0,0.d0)/x(j,j)
            do i=j+1,n
              x(i,j)=x(i,j)*temp
            enddo
          else
            nperm=0
            return
          endif
        endif
      enddo
      return
      end

      subroutine baksubluc(x,v,n,nm,iperm)
! Solves a linear equation x'*w=v for w.  X is the lower-upper
! decomposition of x' from subroutine loupdc.  Vector v is replaced by
! the solutions w as output.
      implicit none
      integer n,nm,iperm(nm),i,j,ii,jj
      double complex x(nm,nm),v(nm),temp
      ii=0
      do i=1,n
        temp=v(iperm(i))
        v(iperm(i))=v(i)
        if (ii.ne.0) then
          do j=ii,i-1
            temp=temp-x(i,j)*v(j)
          enddo
        elseif (temp.ne.0.d0) then
          ii=i
        endif
        v(i)=temp
      enddo
      do i=n,1,-1
        temp=v(i)
        do j=i+1,n
          temp=temp-x(i,j)*v(j)
        enddo
        v(i)=temp/x(i,i)
      enddo
      return
      end


      subroutine inversec(x,y,n,nm)
! Returns matrix y as the inverse of matrix x, both
! of dimension n x n.
      implicit none
      integer n,nm
      double complex x(nm,nm),y(nm,nm),v(nm)
      integer i,j,k,ii,jj,kk
      integer iperm(nm),nperm
!
      double complex ll(nm,nm),uu(nm,nm),AA(nm,nm)
      do i=1,n
        do j=1,n
          if (i.eq.j) then
            y(i,j)=1.d0
          else
            y(i,j)=0.d0
          endif
        enddo
      enddo
      call loupdc(x,n,nm,iperm,nperm)
      do j=1,n
        call baksubluc(x,y(1,j),n,nm,iperm)
      enddo
      return
      end


