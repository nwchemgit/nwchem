      recursive subroutine gw_en(ov1,ov2,factor,x,b,npoles,maxpoles,n,
     $                           ipol,maxiter,rtol,precond)
      implicit none
#include "mafdecls.fh"      
      logical,intent(in) :: precond
      integer,intent(in) :: n,npoles(*),maxpoles,ipol,maxiter
      double precision,intent(in) :: ov1(n,*),ov2(n,*),rtol
      double precision,intent(in) :: factor(maxpoles,*),b(n)
      double precision,intent(out) :: x(n)

      integer,parameter :: maxres = 25
c      double precision,parameter :: rtol = 1d-4

      logical :: converged
      double precision,dimension(n) :: ahr,hr,r,w,eta,xi,precon
      double precision :: u(n,maxres),c(n,maxres),beta(maxres)
      double precision :: alpha(maxres)
      double precision norm, rnorm

      integer iter, totiter, jiter
      integer, external :: ga_nodeid

      iter = 0
      totiter = 0
c      hr = 0d0
c      u(:,:) = 0d0
      converged = .false.

      if (ipol.eq.1) then
        call gw_aprod(OV1,factor,x,r,-4d0,0d0,npoles,n)
      else
        call gw_aprod(OV1,factor,x,r,-2d0,0d0,npoles,n)
        call gw_aprod(OV2,factor(1,2),x,r,-2d0,1d0,npoles(2),n)
      endif
      call ga_dgop(mt_dbl,r,n,'+')

      r(:) = r(:) + b(:) - x(:)

      rnorm = dsqrt(dot_product(r,r))
      if (precond.and.(rnorm.lt.rtol)) goto 200


  100 continue

      iter = iter + 1
      totiter = totiter + 1

      if (precond) then
        hr(:) = 0d0
        call gw_en(ov1,ov2,factor,hr,r,npoles,maxpoles,n,ipol,25,
     $             1d-4,.false.)
c        call gw_minres(ov1,ov2,factor,hr,r,npoles,maxpoles,n,ipol,2,
c     $                 .false.)
      else
c        hr(:) = precon(:)*r(:)
        hr(:) = r(:)
      endif

      if (ipol.eq.1) then
        call gw_aprod(OV1,factor,hr,ahr,4d0,0d0,npoles,n)
      else
        call gw_aprod(OV1,factor,hr,ahr,2d0,0d0,npoles,n)
        call gw_aprod(OV2,factor(1,2),hr,ahr,2d0,1d0,npoles(2),n)
      endif
      call ga_dgop(mt_dbl,ahr,n,'+')
      ahr(:) = (ahr(:) + hr(:))      

      norm = max(dot_product(ahr,ahr),1d-14)
      norm = dot_product(ahr,r)/norm

      xi(:) = r(:) - norm*ahr(:)
      eta(:) = norm*hr(:)

      do jiter=1,iter-1
        alpha(jiter) = dot_product(c(:,jiter),xi)
        xi(:) = xi(:) - alpha(jiter)*c(:,jiter)
        eta(:) = eta(:) + alpha(jiter)*u(:,jiter)
      enddo

      if (precond) then
c        call gw_minres(ov1,ov2,factor,u(1,iter),xi,npoles,maxpoles,
c     $                 n,ipol,2,.false.)
        u(:,iter) = 0d0
        call gw_en(ov1,ov2,factor,u(1,iter),xi,npoles,maxpoles,
     $             n,ipol,25,1d-4,.false.)
      else
c        u(:,iter) = precon(:)*xi(:)
        u(:,iter) = xi(:)
      endif

      if (ipol.eq.1) then
        call gw_aprod(OV1,factor,u(1,iter),c(1,iter),4d0,0d0,npoles,n)
      else
        call gw_aprod(OV1,factor,u(1,iter),c(1,iter),2d0,0d0,npoles,n)
        call gw_aprod(OV2,factor(1,2),u(1,iter),c(1,iter),2d0,1d0,
     $                npoles(2),n)
      endif
      call ga_dgop(mt_dbl,c(1,iter),n,'+')
      c(:,iter) = c(:,iter) + u(:,iter)

      do jiter=1,iter-1
        beta(jiter) = -dot_product(c(:,iter),c(:,jiter))
        c(:,iter) = c(:,iter) + beta(jiter)*c(:,jiter)
        u(:,iter) = u(:,iter) + beta(jiter)*u(:,jiter)
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

c      if (precond.and.(ga_nodeid().eq.0)) write(*,*) totiter,rnorm
c      if ((ga_nodeid().eq.0)) write(*,*) totiter,rnorm


      end subroutine
