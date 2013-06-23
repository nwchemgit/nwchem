      subroutine util_nwsort(n,value,vector,order)
c
c $Id$
c
c Sort the eigenvalues and eigenvectors 
c in an ascending/descending order
c Written by So Hirata, Feb 2002.
c (c) Battelle, PNNL, 2002.
c
      implicit none
      integer n
      double precision value(n)
      double precision vector(n,n)
      double precision minval,maxval,swap
      character*1 order
      integer i,j,k

      if ((order.eq.'A').or.(order.eq.'a')) then
        do i=1,n-1
          minval=value(i)
          k=0
          do j=i+1,n
            if (value(j).lt.minval) then
              k=j
              minval=value(j)
            endif
          enddo
          if (k.ne.0) then
            swap=value(i)
            value(i)=value(k)
            value(k)=swap
            do j=1,n
              swap=vector(j,i)
              vector(j,i)=vector(j,k)
              vector(j,k)=swap
            enddo
          endif
        enddo
      else if ((order.eq.'D').or.(order.eq.'d')) then
        do i=1,n-1
          maxval=value(i)
          k=0
          do j=i+1,n
            if (value(j).gt.maxval) then
              k=j
              maxval=value(j)
            endif
          enddo
          if (k.ne.0) then
            swap=value(i)
            value(i)=value(k)
            value(k)=swap
            do j=1,n
              swap=vector(j,i)
              vector(j,i)=vector(j,k)
              vector(j,k)=swap
            enddo
          endif
        enddo
      endif
      return
      end
