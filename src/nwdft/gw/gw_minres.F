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
      integer,parameter :: maxiter = 10
      double precision,parameter :: rtol = 1d-5
      double precision,parameter :: shift = 0.0d0


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

c      if (istop.eq.5)
c     $  call errquit('gw_minres: minres did not converge',0,0)



      end subroutine




      subroutine gw_minres(OV1,OV2,factor,x,b,npoles,maxpoles,n,ipol,
     $                     maxiter,lprint)
      implicit none
#include "errquit.fh"
#include "global.fh"      
#include "mafdecls.fh"
#include "tcgmsg.fh"      
      logical,intent(in) :: lprint
      integer,intent(in) :: ipol
      integer,intent(in) :: npoles(ipol),maxpoles,n,maxiter
      double precision,intent(in) :: factor(maxpoles,ipol)
      double precision,intent(in) :: OV1(n,maxpoles),OV2(n,maxpoles)
      double precision,dimension(n),intent(in)   :: b
      double precision,dimension(n),intent(inout)  :: x

      integer l_pi, k_pi, l_ipiv, k_ipiv, lscr, l_scr, k_scr, info
      double precision :: eps
      double precision :: Anorm, Acond, Arnorm, rnorm, ynorm
      double precision :: alpha, beta, beta1, cs, dbar, delta, denom
      double precision :: diag, epsa, epsln, epsr, epsx, gam
      double precision :: gbar, gmax, gmin, oldb, oldeps, grnorm, phi
      double precision :: phibar, rhs1, rhs2, rnorml, rootl, Arnorml
      double precision :: relArnorml, s, sn, t, tnorm2, ynorm2, z
      double precision :: t1, t2, test1, test2, root, bnorm, dxnorm
      double precision :: test3, xnorm
      double precision,dimension(n) :: r1, r2, sol, v, w, w1, w2, y
      double precision,dimension(n) :: precon

      
      logical main,converged,oactive
      integer iter, istop, i, gy, me, lo ,hi, ld
      double precision,parameter :: shift = 0.000d0, rtol=1d-4

      me = ga_nodeid()
      main = me.eq.0
      converged = .false.

      eps  = epsilon(eps)
      iter = 0
      istop = 0
      Anorm = 0d0
      Acond = 0d0
      ynorm = 0d0
      gmin  = 1d0
      test3 = 0d0


      if (ipol.eq.1) then
        call gw_aprod(OV1,factor,x,y,-4d0,0d0,npoles,n)
      else
        call gw_aprod(OV1,factor,x,y,-2d0,0d0,npoles,n)
        call gw_aprod(OV2,factor(1,2),x,y,-2d0,1d0,npoles(2),n)
      endif
      call dgop(mt_dbl,y,n,'+')

      beta1 = 0d0
      xnorm = 0d0
!$omp parallel do reduction(+:beta1,xnorm)
      do i=1,n
        y(i) = y(i) + b(i) - x(i)
        beta1 = beta1 + y(i)**2
        xnorm = xnorm + x(i)**2
      enddo
!$omp end parallel do
      beta1 = dsqrt(beta1)
      xnorm = dsqrt(xnorm)

      ! Solution already converged
      if (beta1.lt.rtol) then
        rnorm = beta1
        istop = -1
        goto 2000
      endif

      bnorm = 0d0
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
      gmax = 0d0
      gmin = 1d0

!$omp parallel
!$omp do
      do i=1,n
        r1(i) = y(i)
      enddo
!$omp do
      do i=1,n
        r2(i) = y(i)
      enddo
!$omp do      
      do i=1,n
        w(i) = 0d0
      enddo
!$omp do      
      do i=1,n
        w1(i) = 0d0
      enddo
!$omp do      
      do i=1,n
        w2(i) = 0d0
      enddo
!$omp do      
      do i=1,n
        sol(i) = 0d0
      enddo
!$omp do reduction(+: bnorm)
      do i=1,n
        bnorm = bnorm + b(i)**2
      enddo
!$omp end parallel      
      bnorm = dsqrt(bnorm)



      do
        iter = iter + 1
        s = 1d0/beta


        v(:) = s*y(:)

        if (ipol.eq.1) then
          call gw_aprod(OV1,factor,v,y,4d0,0d0,npoles,n)
        else
          call gw_aprod(OV1,factor,v,y,2d0,0d0,npoles,n)
          call gw_aprod(OV2,factor(1,2),v,y,2d0,1d0,npoles(2),n)
        endif
        call dgop(mt_dbl,y,n,'+')

!$omp   parallel
!$omp   do        
        do i=1,n
          y(i) = y(i) + v(i)
        enddo
!$omp   end do        

        if (iter.ge.2) then
!$omp     do          
          do i=1,n
            y(i) = y(i) - beta/oldb*r1(i)
          enddo
!$omp     enddo        
        endif

!$omp   single
        alpha = 0d0
!$omp   end single

!$omp   do reduction(+:alpha)
        do i=1,n
          alpha = alpha + v(i)*y(i)
        enddo
!$omp   end do

!$omp   do
        do i=1,n
          y(i) = y(i) - alpha/beta*r2(i)
        enddo
!$omp   end do

!$omp   do
        do i=1,n
          r1(i) = r2(i)
        enddo
!$omp   end do

!$omp   do
        do i=1,n
          r2(i) = y(i)
        enddo
!$omp   end do          


!$omp   single        
        oldb = beta
        beta = 0d0
!$omp   end single


!$omp   do reduction(+:beta)
        do i=1,n
          beta = beta + r2(i)**2
        enddo
!$omp   end do

!$omp   single
        beta = dsqrt(beta)

        tnorm2 = tnorm2 + alpha**2 + oldb**2 + beta**2

        if (iter.eq.1) then
          if (beta/beta1 .le. 10d0*eps) then
            istop = -1
          endif
        endif

        oldeps = epsln
        delta = cs*dbar + sn*alpha
        gbar  = sn*dbar - cs*alpha
        epsln = sn*beta
        dbar  = -cs*beta
        root = dsqrt(gbar**2 + dbar**2)
        Arnorm = phibar*root

        if (beta.eq.0d0) then
          cs = sign(1d0,gbar)
          sn = 0d0
          gam = abs(gbar)
        elseif (gbar.eq.0d0) then
          cs = 0d0
          sn = 1d0
          gam = beta
        elseif (beta.gt.abs(gbar)) then
          t = gbar/beta
          sn = 1d0/dsqrt(1d0 + t**2)
          cs = sn*t
          gam = beta/sn
        else
          t = beta/gbar
          cs = sign(1d0,gbar)/dsqrt(1d0+t**2)
          sn = cs*t
          gam = gbar/cs
        endif
        gam = max(gam,eps)

        phi    = cs*phibar
        phibar = sn*phibar
!$omp   end single

!$omp   do
        do i=1,n
          w1(i) = w2(i)
        enddo
!$omp   end do
!$omp   do
        do i=1,n
          w2(i) = w(i)
        enddo
!$omp   end do
!$omp   do
        do i=1,n
          w(i) = (v(i) - oldeps*w1(i) - delta*w(i))/gam
        enddo
!$omp   end do
!$omp   do
        do i=1,n
          sol(i) = sol(i) + phi*w(i)
        enddo
!$omp   end do
!$omp   single
        dxnorm = 0d0
        xnorm = 0d0
!$omp   end single
!$omp   do reduction(+:dxnorm,xnorm)
        do i=1,n
          dxnorm = dxnorm + sol(i)**2
          xnorm = xnorm + (sol(i)+x(i))**2
        enddo
!$omp   end do        
!$omp   end parallel

        dxnorm = dsqrt(dxnorm)
        xnorm  = dsqrt(xnorm)

        gmax = max(gmax,gam)
        gmin = min(gmin,gam)
        z    = rhs1/gam
        rhs1 = rhs2 - delta*z
        rhs2 = -epsln*z

        Anorm = dsqrt(tnorm2)
        rnorm = phibar

        Acond = gmax/gmin
        if (istop.ne.0) exit

        epsx = (Anorm*dxnorm + beta1)*eps
        epsr = (Anorm*dxnorm + beta1)*rtol
        test1 = rnorm/(Anorm*dxnorm + bnorm)
        test2 = Arnorm/(Anorm*(rnorm + eps))
        test3 = xnorm*rnorm + rnorm**2/gmin
        test3 = max(test3, Anorm*(test3/gmin+rnorm*xnorm/gmin))
        t1 = 1d0 + test1
        t2 = 1d0 + test2

        if (iter  .ge. maxiter  ) istop=6
        if (Acond .ge. 0.1d0/eps) istop=5
        if (t2    .le. 1d0      ) istop=2
        if (t1    .le. 1d0      ) istop=1
        if (epsx  .ge. beta1    ) istop=3
        if (iter.gt.1 .and. test3 .le. rtol/10d0) istop=10
        if (test2 .le. rtol     ) istop=2
        if (test1 .le. rtol     ) istop=1
        if (rnorm .le. rtol     ) istop=8


        if (istop.ne.0) exit
      enddo

      x(:) = x(:) + sol(:)

 1000 continue
 2000 continue

      if (istop.eq.6) then
        if (.not.ma_alloc_get(mt_dbl,n**2,'pi',l_pi,k_pi))
     $    call errquit('gw_en: can''t allocate pi matrix',0,MA_ERR) 
        if (.not.ma_push_get(mt_int,n,'ipiv',l_ipiv,k_ipiv))
     $    call errquit('gw_en: can''t allocate ipiv',0,MA_ERR) 
        call dsytrf('l',n,dbl_mb(k_pi),n,int_mb(k_ipiv),
     $              x,-1,info)
        lscr = max(n,int(x(1)))
        x(:) = b(:)
        if (.not.ma_push_get(mt_dbl,lscr,'scr',l_scr,k_scr))
     $    call errquit('gw_en: can''t allocate scratch',0,MA_ERR) 

        call gw_buildpi('r',ov1,ov2,dbl_mb(k_pi),factor,n,n,npoles,
     $                   main,ipol,maxpoles)
        call ga_dgop(mt_dbl,dbl_mb(k_pi),n**2,'+')
        call dsysv('l',n,1,dbl_mb(k_pi),n,int_mb(k_ipiv),x,n,
     $              dbl_mb(k_scr),lscr,info)
#ifdef SCALAPACK
c        call gw_symmat(dbl_mb(k_pi),n,1,n)
c        call gw_lu_solve(dbl_mb(k_pi),x,n)
#endif
        if (.not.ma_chop_stack(l_ipiv))
     $    call errquit('gw_en: can''t chop stack',0,MA_ERR)
        if (.not.ma_free_heap(l_pi))
     $    call errquit('gw_en: can''t deallocate pi matrix',0,MA_ERR)   
        rnorm = 0d0
      endif


c      if (main.and.lprint) write(*,*) istop,test3

      end subroutine


      subroutine gw_aprod(ov,factor,x,y,scal,mult,npoles,n)
      implicit none
      integer,intent(in) :: npoles, n
      double precision, intent(in) :: scal, mult
      double precision, intent(in) :: ov(n,npoles),factor(npoles)
      double precision, intent(in) :: x(n)
      double precision, intent(out) :: y(n)

      integer :: i
      double precision :: tmp(npoles)

      call dgemv('t',n,npoles,scal,ov,n,x,1,0d0,tmp,1)

!$omp parallel do
      do i=1,npoles
        tmp(i) = tmp(i)*factor(i)
      enddo
!$omp end parallel do

      call dgemv('n',n,npoles,1d0,ov,n,tmp,1,mult,y,1)

      end subroutine

      subroutine gw_precon(ov,factor,precon,scal,npoles,n)
      implicit none
      integer,intent(in) :: npoles, n
      double precision,intent(in) :: scal
      double precision,intent(in) :: ov(n,npoles), factor(npoles)
      double precision,intent(out) :: precon(n)

      integer i
      double precision :: temp(n)

      temp(:) = factor(1)*ov(:,1)**2

      do i=2,npoles
        temp(:) = temp(:) + factor(i)*ov(:,i)**2
      enddo
      call daxpy(n,scal,temp,1,precon,1)

      end subroutine


