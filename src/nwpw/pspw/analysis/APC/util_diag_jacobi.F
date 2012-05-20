*     **** subroutine util_diag_jacobi  --  jacobi matrix diagonalization ****

*     Matrix diagonalization of a real symmetric matrix using Jacobi rotations
*
*     unit:  output unit
*     n:     logical dimension of the matrix to be diagonalized
*     np:    physical dimension of the matrix storage area
*     a:     input with the matrix to be diagonalized (real symmetric)
*     d:     returned with the eigenvalues in ascending order
*     v:     returned with the eigenvectors of the matrix
*     b:     temporary work vector
*     z:     temporary work vector

      subroutine util_diag_jacobi(unit,n,np,a,d,v,b,z)

      implicit none

      integer i,j,k,ip,iq,n,np,nrot,maxrot,unit
      real*8  sm,tresh,s,c,t,theta,tau,h,g,p
      real*8  a(np,np),d(np),v(np,np),b(np),z(np)

*     **** preliminaries ****
      maxrot = 100
      nrot = 0
      do ip = 1, n
         do iq = 1, n
            v(ip,iq) = 0.0d0
         end do
         v(ip,ip) = 1.0d0
      end do
      do ip = 1, n
         b(ip) = a(ip,ip)
         d(ip) = b(ip)
         z(ip) = 0.0d0
      end do

*     **** perform jacobi rotations ****
      do i = 1, maxrot
         sm = 0.0d0
         do ip = 1, n-1
            do iq = ip+1, n
               sm = sm + abs(a(ip,iq))
            end do
         end do
         if (sm .eq. 0.0d0)  goto 10
         if (i .lt. 4) then
            tresh = 0.2d0*sm / n**2
         else
            tresh = 0.0d0
         end if
         do ip = 1, n-1
            do iq = ip+1, n
               g = 100.0d0 * abs(a(ip,iq))
               if (i.gt.4 .and. abs(d(ip))+g.eq.abs(d(ip))
     &                    .and. abs(d(iq))+g.eq.abs(d(iq))) then
                  a(ip,iq) = 0.0d0
               else if (abs(a(ip,iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h)+g .eq. abs(h)) then
                     t = a(ip,iq) / h
                  else
                     theta = 0.5d0*h / a(ip,iq)
                     t = 1.0d0 / (abs(theta)+sqrt(1.0d0+theta**2))
                     if (theta .lt. 0.0d0)  t = -t
                  end if
                  c = 1.0d0 / sqrt(1.0d0+t**2)
                  s = t * c
                  tau = s / (1.0d0+c)
                  h = t * a(ip,iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip,iq) = 0.0d0
                  do j = 1, ip-1
                     g = a(j,ip)
                     h = a(j,iq)
                     a(j,ip) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = ip+1, iq-1
                     g = a(ip,j)
                     h = a(j,iq)
                     a(ip,j) = g - s*(h+g*tau)
                     a(j,iq) = h + s*(g-h*tau)
                  end do
                  do j = iq+1, n
                     g = a(ip,j)
                     h = a(iq,j)
                     a(ip,j) = g - s*(h+g*tau)
                     a(iq,j) = h + s*(g-h*tau)
                  end do
                  do j = 1, n
                     g = v(j,ip)
                     h = v(j,iq)
                     v(j,ip) = g - s*(h+g*tau)
                     v(j,iq) = h + s*(g-h*tau)
                  end do
                  nrot = nrot + 1
               end if
            end do
         end do
         do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
         end do
      end do

*     **** print warning if not converged **** 
   10 continue
      if (nrot .eq. maxrot) then
         write(unit,20)
   20    format (/,' util_diag_jacobi -- Diagonalization not converged')
      end if

      return

*     **** sort eigenvalues and eigenvectors in ascending order ****
c      if (lsort) then
c       do i = 1, n-1
c         k = i
c         p = d(i)
c         do j = i+1, n
c            if (d(j) .lt. p) then
c               k = j
c               p = d(j)
c            end if
c         end do
c         if (k .ne. i) then
c            d(k) = d(i)
c            d(i) = p
c            do j = 1, n
c               p = v(j,i)
c               v(j,i) = v(j,k)
c               v(j,k) = p
c            end do
c         end if
c       end do
c      end if ! lsort
c      return

      end
c $Id$
