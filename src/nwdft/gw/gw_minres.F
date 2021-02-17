      subroutine minres(A,x,b,n)

      implicit none
#include "errquit.fh"
      integer,intent(in) :: n
      double precision,dimension(n,n),intent(in) :: A
      double precision,dimension(n),intent(in)   :: b
      double precision,dimension(n),intent(inout)  :: x

      double precision :: eps
      double precision :: Anorm, Acond, Arnorm, rnorm, ynorm
      double precision :: alpha, beta, beta1, cs, dbar, delta, denom
      double precision :: diag, epsa, epsln, epsr, epsx, gam
      double precision :: gbar, gmax, gmin, oldb, oldeps, grnorm, phi
      double precision :: phibar, rhs1, rhs2, rnorml, rootl, Arnorml
      double precision :: relArnorml, s, sn, t, tnorm2, ynorm2, z
      double precision,dimension(n) :: r1, r2, sol, v, w, w1, w2, y

      integer iter, istop, i
      integer,parameter :: maxiter = 30
      double precision,parameter :: rtol = 1d-5
      double precision,parameter :: shift = 0.0d0
#ifdef USE_OPENMP
      integer iMaxThreads
      integer, external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif


      eps  = epsilon(eps)
      iter = 0
      istop = 0
      Anorm = 0d0
      Acond = 0d0
      rnorm = 0d0
      ynorm = 0d0

      ! Setup y and v for the first Lanczos vector v1

!$omp parallel do      
      do i=1,n
        y(i) = b(i)
      enddo
!$omp end parallel do      

      call dsymv('l',n,-1d0,A,n,x,1,1d0,y,1)

!$omp parallel do      
      do i=1,n
        r1(i) = y(i)
      enddo
!$omp end parallel do

      beta1 = dot_product(r1,r1)

      ! Solution already converged
      if (beta1.eq.0d0) then
        goto 1000
        istop = -1
      endif

      beta1 = dsqrt(beta1)

      call dsymv('l',n,1d0,A,n,y,1,0d0,w,1)
      Arnorml = dsqrt(dot_product(w,w))

      oldb = 0d0
      beta = beta1
      dbar = 0d0
      epsln = 0d0
      grnorm = beta1
      phibar = beta1
      rhs1 = beta1
      rhs2 = 0d0
      tnorm2 = 0d0
      ynorm2 = 0d0
      cs = -1d0
      sn = 0d0
      w(:) = 0d0
      w1(:) = 0d0
      w2(:) = 0d0
      r2(:) = r1
      sol = 0d0

      do
        iter = iter + 1
        s = 1d0/beta

!$omp   parallel do        
        do i=1,n
          v(i) = s*y(i)
        enddo
!$omp   end parallel do        

        call dsymv('l',n,1d0,A,n,v,1,0d0,y,1)

!$omp   parallel do
        do i=1,n
          y(i) = y(i) - shift*v(i)
        enddo
!$omp   end parallel do        

        if (iter.ge.2) then
!$omp     parallel do
          do i=1,n
            y(i) = y(i) - (beta/oldb)*r1(i)
          enddo
!$omp     end parallel do          
        endif

        alpha = dot_product(v,y)

!$omp   parallel do 
        do i=1,n
          r1(i) = r2(i)
          y(i) = y(i) - (alpha/beta)*r2(i)
          r2(i) = y(i)
        enddo
!$omp   end parallel do

        oldb = beta
        beta = dsqrt(dot_product(r2,y))
        tnorm2 = tnorm2 + alpha**2 + oldb**2 + beta**2

        if (iter.eq.1) then
          if (beta/beta1 .le. 10d0*eps) then
            istop = -1
          endif
          gmax = abs(alpha)
          gmin = gmax
        endif

        oldeps = epsln
        delta = cs*dbar + sn*alpha
        gbar  = sn*dbar - cs*alpha
        epsln = sn*beta
        dbar  = -cs*beta

        gam    = dsqrt( gbar**2 + beta**2 )
        gam    = max(gam,eps)
        cs     = gbar/gam
        sn     = beta/gam
        phi    = cs*phibar
        phibar = sn*phibar

!$omp   parallel do        
        do i=1,n
          w1(i) = w2(i)
          w2(i) = w(i)
          w(i) = (v(i) - oldeps*w1(i) - delta*w2(i))/gam
          sol(i) = sol(i) + phi*w(i)
        enddo
!$omp   end parallel do        

        gmax = max(gmax,gam)
        gmin = min(gmin,gam)
        z    = rhs1/gam
        ynorm2 = z**2 + ynorm2
        rhs1 = rhs2 - delta*z
        rhs2 = -epsln*z

        Anorm = dsqrt(tnorm2)
        ynorm = dsqrt(ynorm2)
        epsa = Anorm*eps
        epsx = Anorm*ynorm*eps
        epsr = Anorm*ynorm*rtol
        diag = gbar
        if (diag.eq.0) diag = epsa

        grnorm = phibar
        rnorml = rnorm
        rnorm  = grnorm
        rootl  = dsqrt(gbar**2 + dbar**2)
        Arnorml = rnorml*rootl
        relArnorml = rootl/Anorm

        Acond = gmax/gmin

        if (istop .eq. 0) then
          if (iter .ge. maxiter) istop = 5
          if (Acond .ge. 0.1d0/eps) istop=4
          if (epsx  .ge. beta1    ) istop=3
          if (grnorm .le. epsx .or. relArnorml .le. epsx) istop = 2
          if (grnorm .le. epsr .or. relArnorml .le. epsr) istop = 1
c          if (grnorm/dsqrt(dble(n)) .le. rtol ) istop = 1
        endif

        if (istop.ne.0) exit
      enddo

 1000 Arnorm = Arnorml

!$omp parallel do      
      do i=1,n
        x(i) = x(i) + sol(i)
      enddo
!$omp end parallel do      

      if (istop.eq.4) 
     $  call errquit('gw_minres: matrix is very ill-conditioned',0,0)

      if (istop.eq.5)
     $  call errquit('gw_minres: minres did not converge',0,0)


#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif      

      end subroutine





