C> \ingroup wfn1
C> @{
C>
C> \brief Test program to check that we can correctly compute parabola
C>
C> I have a couple of routines that try and fit parabola through a
C> couple of data points. If these routines work correctly I they should
C> be able to reproduce the parabola I used to generate the data points.
C>
      program wfn1_test_para
      implicit none
      double precision ai,bi,ci !< Input a, b, and c
      double precision ao1,bo1,co1 !< Output a, b, and c
      double precision ao2,bo2,co2 !< Output a, b, and c
      double precision x0,x1,x2,x3 !< Input x-coordinates
      double precision f0,f1,f2,f3 !< Input function values at
                                   !< x-coordinates
      double precision df0         !< Input the gradient at x0
      double precision t           !< Input the trust region
      double precision xm1         !< Output location of minimum
      double precision fxm1        !< Output function value at minimum
      double precision xm2         !< Output location of minimum
      double precision fxm2        !< Output function value at minimum
c
      double precision x
      double precision fx,dfx
      fx(x) = ai*x*x+bi*x+ci
      dfx(x) = 2.0d0*ai*x+bi
c
10    read(*,*,ERR=10)ai,bi,ci
      if (ai.eq.0.0d0.and.bi.eq.0.0d0.and.ci.eq.0.0d0) goto 900
20    read(*,*,ERR=20)x1,x2,x3
      x0  = 0.0d0
      f0  = fx(0.0d0)
      df0 = dfx(0.0d0)
      f1  = fx(x1)
      f2  = fx(x2)
      f3  = fx(x3)
      t   = 5.0d0
      call wfn1_f0df0f1(x1,f0,df0,f1,t,ao1,bo1,co1,xm1,fxm1)
      call wfn1_f1f2f3(x1,x2,x3,f1,f2,f3,t,ao2,bo2,co2,xm2,fxm2)
c
      write(*,*)ao1,bo1,co1,xm1,fxm1
      write(*,*)ao2,bo2,co2,xm2,fxm2
      write(*,*)
      goto 10
c
900   continue
      end
