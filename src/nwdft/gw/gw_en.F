      subroutine gw_en(ov1,ov2,factor,x,b,npoles,maxpoles,n,ipol)
      implicit none
#include "mafdecls.fh"      
#include "global.fh"
#include "errquit.fh"
      integer,intent(in) :: n,npoles(*),maxpoles,ipol
      double precision,intent(in) :: ov1(n,*),ov2(n,*)
      double precision,intent(in) :: factor(maxpoles,*),b(n)
      double precision,intent(inout) :: x(n)

      integer,parameter :: maxres = 25
      integer,parameter :: maxiter = 10
      double precision,parameter :: rtol = 1d-4

      logical :: converged, main
      double precision,dimension(n) :: ahr,hr,r,w,eta,xi,precon
      double precision :: u(n,maxres),c(n,maxres),beta
      double precision :: alpha
      double precision norm, rnorm, xnorm

      integer iter, totiter, jiter, ga, gb
      integer l_pi, k_pi, l_ipiv, k_ipiv, l_scr, k_scr, lscr, info

      iter = 0
      totiter = 0
      converged = .false.
      main = ga_nodeid().eq.0


      xnorm = dot_product(x,x)

      if (xnorm.eq.0d0) then
        r(:) = b(:)
      else
        if (ipol.eq.1) then
          call gw_aprod(OV1,factor,x,r,-4d0,0d0,npoles,n)
        else
          call gw_aprod(OV1,factor,x,r,-2d0,0d0,npoles,n)
          call gw_aprod(OV2,factor(1,2),x,r,-2d0,1d0,npoles(2),n)
        endif
        call dgop(mt_dbl,r,n,'+')
        r = r + b - x
      endif

      rnorm = dsqrt(dot_product(r,r))

      if (rnorm.lt.rtol) goto 200

      if ((rnorm.gt.0.01d0)) then
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

        goto 200
      endif

      hr = 0d0

  100 continue

      iter = iter + 1
      totiter = totiter + 1

      hr(:) = r(:)

      if (ipol.eq.1) then
        call gw_aprod(OV1,factor,hr,ahr,4d0,0d0,npoles,n)
      else
        call gw_aprod(OV1,factor,hr,ahr,2d0,0d0,npoles,n)
        call gw_aprod(OV2,factor(1,2),hr,ahr,2d0,1d0,npoles(2),n)
      endif
      call dgop(mt_dbl,ahr,n,'+')
      ahr = ahr + hr

      norm = max(dot_product(ahr,ahr),1d-14)
      norm = dot_product(ahr,r)/norm

      xi(:) = r(:) - norm*ahr(:)
      eta(:) = norm*hr(:)

      do jiter=1,iter-1
        alpha = dot_product(c(:,jiter),xi)
        call daxpy(n,-alpha,c(1,jiter),1,xi,1)
        call daxpy(n,alpha,u(1,jiter),1,eta,1)
      enddo

      u(:,iter) = xi(:)


      if (ipol.eq.1) then
        call gw_aprod(OV1,factor,u(1,iter),c(1,iter),4d0,0d0,npoles,n)
      else
        call gw_aprod(OV1,factor,u(1,iter),c(1,iter),2d0,0d0,npoles,n)
        call gw_aprod(OV2,factor(1,2),u(1,iter),c(1,iter),2d0,1d0,
     $                npoles(2),n)
      endif
      call dgop(mt_dbl,c(1,iter),n,'+')
      c(:,iter) = c(:,iter) + u(:,iter)

      do jiter=1,iter-1
        beta = -dot_product(c(:,iter),c(:,jiter))
        call daxpy(n,beta,c(1,jiter),1,c(1,iter),1)
        call daxpy(n,beta,u(1,jiter),1,u(1,iter),1)
      enddo

      norm = max(dsqrt(dot_product(c(:,iter),c(:,iter))),1d-14)

      c(:,iter) = c(:,iter)/norm
      u(:,iter) = u(:,iter)/norm

      norm = dot_product(c(:,iter),xi)

      x(:) = x(:) + eta(:) + norm*u(:,iter)
      r(:) = xi(:) - norm*c(:,iter)

      rnorm = dsqrt(dot_product(r,r))


      if (rnorm.lt.rtol) goto 200
      if (totiter.ge.maxiter) goto 200
      if (iter.ge.maxres) iter = 0


      goto 100

  200 continue    


      end subroutine
