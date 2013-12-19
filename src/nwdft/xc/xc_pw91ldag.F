      Subroutine xc_pw91ldag(rs, i, e ,d1e, d2e)
      implicit none
      double precision rs ! [in]
      integer i ! [in]
c     i=1 eps_c(rs,0)
c     i=2 eps_c(rs,1)
c     i=3 alpha_c(rs)
      double precision e ! [out]
      double precision d1e ! [out]
      double precision d2e ! [out]
c
c
      double precision h1, h2, d1h1,d2h1,d2h2,d1x,d2x,
     D     d1h2 ,x
c
c     Functional Parameters
c
      double precision A(3), alp(3), b(4,3)
      save A, alp, b
      integer j,n, initial
      save initial
      data A   / 0.0310907d0, 0.01554535d0, 0.0168869d0 /
      data alp / 0.21370d0, 0.20548d0, 0.11125d0 /
      data b   / 7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0,
     &          14.1189d0, 6.1977d0, 3.3662d0, 0.62517d0,
     &          10.357d0, 3.6231d0, 0.88026d0, 0.49671d0 /
      data initial /1/
c
c     Define miscellaneous parameters.
c
      if (initial.eq.1)then
         initial = 0
c        For convenience store -2A as A and multiply betas by 2A
         do j = 1, 3
            A(j) = -2d0*A(j)
            do n = 1, 4
               b(n,j) = -A(j)*b(n,j)
            enddo
         enddo
c        Finally, change the sign on A for spin stiffness since
c        the negative of that is fitted in the PW'91 paper.  We can't
c        just take the negative of A at the start since A also contributes
c        to the argument of the ln function.
         A(3) = -A(3)
      endif
      x = sqrt(rs)
      d1x = 0.5d0/x
      d2x = -0.5d0*d1x/rs
c
      h2 = x*(b(1,i) + x*(b(2,i) + x*(b(3,i) + x*b(4,i))))
      d1h2 = b(1,i)
     &     + x*(2d0*b(2,i) + x*(3d0*b(3,i) + 4d0*x*b(4,i)))
      
      d2h2 = 2d0*b(2,i) + x*(6d0*b(3,i) + 12d0*x*b(4,i))
      
c     
      h1 = DLOG(1d0+1d0/h2)
      d1h1 = -d1h2/(h2*(h2+1d0))
      
      d2h1 = d1h1*d1h1*(2d0*h2+1d0) - d2h2/(h2*(h2+1d0))
      
c     
      e = A(i)*(1d0+alp(i)*rs)*h1
      
      d1e = A(i)*(2d0*alp(i)*x*h1+(1d0+alp(i)*rs)*d1h1)
      
      d2e = A(i)*(2d0*alp(i)*h1+4d0*alp(i)*x*d1h1
     &     +(1d0+alp(i)*rs)*d2h1)
      
c
c           Transform derivatives wrt x to derivatives wrt rs
c     

c           Do 2nd derivative first so the x first derivative in d1e
c           is not lost
      d2e = d2e*d1x*d1x + d1e*d2x
      
      d1e = d1e*d1x

      return
      end
c $Id$
