      subroutine gw_minres(OV1,OV2,factor,x,b,npoles,maxpoles,n,ipol)
C
C     Obtain the action of the inverse dielectric matrix on a given ERI
C     vector by using the MINRES solver. (See Equations (43)-(46))
C     without the explicit computation of the polarizability or
C     dielectric matrices.
C
C     If the convergence seems to be slow, the algorithm switches to
C     actually building the dielectric matrix and decompose it with the
C     LDL^T Bunch-Kaufmann factorization.
C

      use iso_c_binding
      implicit none
#include "errquit.fh"
#include "global.fh"      
#include "mafdecls.fh"
#include "util.fh"
      integer,parameter :: maxiter = 25
      integer,parameter :: maxthreads = 1024

      integer,intent(in) :: ipol
      integer,intent(in) :: npoles(ipol),maxpoles,n
      double precision,intent(in) :: factor(maxpoles,ipol)
      double precision,intent(in) :: OV1(n,maxpoles),OV2(n,maxpoles)
      double precision,dimension(n),intent(in)   :: b
      double precision,dimension(n),intent(inout)  :: x


      integer l_pi, k_pi, l_ipiv, k_ipiv, lscr, l_scr, k_scr, info
      integer mynri, nprocs, llpi, pilo
      double precision :: eps, dummy
      double precision :: Anorm, Acond, Arnorm, rnorm, ynorm
      double precision :: alpha, beta, beta1, cs, dbar, delta, denom
      double precision :: diag, epsa, epsln, epsr, epsx, gam
      double precision :: gbar, gmax, gmin, oldb, oldeps, grnorm, phi
      double precision :: phibar, rhs1, rhs2, rnorml, rootl, Arnorml
      double precision :: relArnorml, s, sn, t, tnorm2, ynorm2, z
      double precision :: t1, t2, test1, test2, root, bnorm, dxnorm
      double precision :: test3, xnorm, x0norm
      double precision :: aprod,rest,t0,tdgop

      integer(kind=c_int) :: stat_alloc
      real(kind=C_DOUBLE),allocatable,dimension(:) :: r1, r2, sol, v, w,
     $                                                w1, w2, y

      
      logical main,converged,oactive,parallel_diag
      integer iter, istop, i, gy, me, ld, nmax, nlocal, ithread
      integer lo, hi, jthread
      double precision,parameter :: shift = 0.000d0, rtol=1d-5

      call ga_sync()

      allocate(y(n), w(n), w1(n), w2(n), sol(n), r1(n), r2(n), v(n),
     $         stat=stat_alloc)
      if (stat_alloc.ne.0) then
        call errquit('gw_minres: allocation failed',0,MA_ERR)
      endif

      me = ga_nodeid()
      main = me.eq.0
      converged = .false.

      !
      ! Initialize MINRES variables
      !
      eps  = epsilon(eps)
      iter = 0
      istop = 0
      Anorm = 0d0
      Acond = 0d0
      ynorm = 0d0
      gmin  = 1d0
      test3 = 0d0
      x0norm = dot_product(x,x)

      !
      ! Check the starting residual.
      ! If the solution vector is 0, just skip the computation
      !
      if (x0norm.eq.0d0) then
!$omp   parallel do simd schedule(static)
        do i=1,n
          y(i) = b(i)
        end do
!$omp   end parallel do simd
        goto 1000
      endif

      !
      ! Obtain the action of the polarizability matrix on a vector
      ! (bottom line of Equation (47)).
      !
      if (ipol.eq.1) then
        call gw_aprod(OV1,factor,x,y,-4d0,0d0,npoles,n)
      else
        call gw_aprod(OV1,factor,x,y,-2d0,0d0,npoles,n)
        call gw_aprod(OV2,factor(1,2),x,y,-2d0,1d0,npoles(2),n)
      endif
      call ga_mask_sync(.true.,.false.)
      call ga_dgop(1038,y,n,'+')

      !
      ! Obtain the residual by completing the action of the dielectric
      ! matrix (top line of Equation (47)).
      !
!$omp parallel do simd schedule(static)
      do i=1,n
        y(i) = y(i) + b(i) - x(i)
      enddo
!$omp end parallel do simd

1000  beta1 = util_norm2(n,y)

      !
      ! Solution already converged
      !
      if (beta1.lt.rtol) then
        rnorm = beta1
        istop = -1
        goto 2000
      endif

      !
      ! Continue with MINRES
      !
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

!$omp parallel do simd schedule(static)
      do i=1,n
        r1(i) = y(i)
        r2(i) = y(i)
        w(i) = 0d0
        w1(i) = 0d0
        w2(i) = 0d0
        sol(i) = 0d0
      enddo
!$omp end parallel do simd

      bnorm = util_norm2(n,b)

      do
        iter = iter + 1
        s = 1d0/beta

!$omp   parallel do simd schedule(static)
        do i=1,n
          v(i) = s*y(i)
        enddo
!$omp   end parallel do simd

        !
        ! Obtain the action of the polarizability matrix
        ! (bottom line of Equation (47)).
        !
        if (ipol.eq.1) then
          call gw_aprod(OV1,factor,v,y,4d0,0d0,npoles,n)
        else
          call gw_aprod(OV1,factor,v,y,2d0,0d0,npoles,n)
          call gw_aprod(OV2,factor(1,2),v,y,2d0,1d0,npoles(2),n)
        endif
        call ga_mask_sync(.true.,.false.)
        call ga_dgop(1038,y,n,'+')

        !
        ! Complete the action of the dielectric matrix
        ! (top line of Equation (47)).
        !
        call daxpy(n,1d0,v,1,y,1)

        if (iter.ge.2) call daxpy(n,-beta/oldb,r1,1,y,1)
        alpha = dot_product(v,y)
        oldb = beta
        call daxpy(n,-alpha/beta,r2,1,y,1)

!$omp   parallel do simd schedule(static)
        do i=1,n
          r1(i) = r2(i)
          r2(i) = y(i)
        enddo
!$omp   end parallel do simd

        beta = util_norm2(n,r2)

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

!$omp   parallel do schedule(static)
        do i=1,n
          w1(i) = w2(i)
          w2(i) = w(i)
          w(i) = (v(i) - oldeps*w1(i) - delta*w(i))/gam
        enddo
!$omp   end parallel do

        call daxpy(n,phi,w,1,sol,1)

        dxnorm = util_norm2(n,sol)
        xnorm = dsqrt(x0norm + dxnorm**2 + 2d0*dot_product(x,sol))
        
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
        t1 = 1d0 + test1
        t2 = 1d0 + test2

        !
        ! Check how slow is MINRES convergence. If it seems to take many
        ! iterations, then explicitly build the dielectric matrix.
        !
        if (iter  .ge. maxiter  )                   istop = 101
        if (iter  .ge.  5 .and. rnorm .ge. 5.0d0)   istop = 102
        if (iter  .ge. 10 .and. rnorm .ge. 1.0d0)   istop = 103
        if (iter  .ge. 15 .and. rnorm .ge. 0.5d0)   istop = 104
        if (iter  .ge. 20 .and. rnorm .ge. 0.01d0)  istop = 105

        if (Acond .ge. 0.1d0/eps) istop=5
        if (t2    .le. 1d0      ) istop=2
        if (t1    .le. 1d0      ) istop=1
        if (epsx  .ge. beta1    ) istop=3
        if (test3 .le. rtol     ) istop=4
        if (test2 .le. rtol     ) istop=2
        if (test1 .le. rtol     ) istop=1
        if (rnorm .le. rtol     ) istop=8

        call ga_mask_sync(.true.,.false.)
        call ga_brdcst(1038,istop,ma_sizeof(mt_int,1,mt_byte),0)

        if (istop.ne.0) exit
      enddo

      !
      ! Update the solution vector
      !
      call daxpy(n,1d0,sol,1,x,1)

 2000 continue


c     if we reached this point, convergence with
c     MINRES seems rather slow. Let's better try
c     actually building the matrix and solving the
c     equation system using LAPACK.
      if ((istop.gt.100 .and. rnorm.gt.10d0*rtol)) then

c       we need an n**2 matrix
        if (.not.ma_push_get(mt_dbl,n**2,'pi',l_pi,k_pi))
     $    call errquit('gw_en: can''t allocate pi matrix',0,MA_ERR) 

        !
        ! Build the polarizability matrix and transform it
        !
        call gw_cdgw_buildpi('r',ov1,ov2,dbl_mb(k_pi),factor,n,n,npoles,
     $                   main,ipol,maxpoles)
        call ga_mask_sync(.true.,.false.)
        call ga_dgop(1038,dbl_mb(k_pi),(n*(n+1))/2,'+')

c       Pi is built in rectangular full packed format.
c       Transform it to standard format.
        call rfp2std(dbl_mb(k_pi),n)

c       Use DSYSV locally
        lscr = 64*n
        x(:) = b(:)
        if (.not.ma_push_get(mt_int,n,'ipiv',l_ipiv,k_ipiv))
     $    call errquit('gw_en: can''t allocate ipiv',0,MA_ERR) 
        if (.not.ma_push_get(mt_dbl,lscr,'scr',l_scr,k_scr))
     $    call errquit('gw_en: can''t allocate scratch',0,MA_ERR) 

        call dsysv('l',n,1,dbl_mb(k_pi),n,int_mb(k_ipiv),x,n,
     $              dbl_mb(k_scr),lscr,info)

        if (.not.ma_chop_stack(l_pi))
     $    call errquit('gw_en: can''t deallocate pi matrix',0,MA_ERR)   
        rnorm = 0d0

c     MINRES converged, broadcast result from node 0
      else
        call ga_mask_sync(.true.,.false.)
        call ga_brdcst(1038,x,ma_sizeof(mt_dbl,n,mt_byte),0)
      endif

c     Deallocate local fields
      deallocate(r1, r2, y, w, w1, w2, v, sol, stat=stat_alloc)
      if (stat_alloc.ne.0) then
        call errquit('gw_minres: error deallocating',0,MA_ERR)
      endif

      end subroutine


      subroutine gw_aprod(ov,factor,x,y,scal,mult,npoles,n)
C
C     Computes the action of the polarizability matrix on a given vector
C
      implicit none
      integer,intent(in) :: npoles, n
      double precision, intent(in) :: scal, mult
      double precision, intent(in) :: ov(n,npoles),factor(npoles)
      double precision, intent(in) :: x(n)
      double precision, intent(out) :: y(n)

      integer :: i,tid,lo,hi
      double precision :: tmp(npoles)

      !
      ! Rightmost sum on the bottom line of Equation (47)
      !
      call dgemv('t',n,npoles,scal,ov,n,x,1,0d0,tmp,1)

      !
      ! Factor in square brackets of Equation (47)
      !
!$omp parallel do simd schedule(static)
      do i=1,npoles
        tmp(i) = tmp(i)*factor(i)
      enddo
!$omp end parallel do simd

      !
      ! First sum on the bottom line of Equation (47)
      !
      call dgemv('n',n,npoles,1d0,ov,n,tmp,1,mult,y,1)

      end subroutine

      
      subroutine minres(A,x,b,n)
C
C     This version of MINRES explicitly uses the dielectric matrix A
C
      implicit none
#include "errquit.fh"
#include "global.fh"      
#include "mafdecls.fh"
#include "util.fh"
      integer,intent(in) :: n
      double precision,dimension(n),intent(in)   :: b
      double precision,dimension(n,n),intent(in) :: A
      double precision,dimension(n),intent(inout)  :: x

      integer,parameter :: maxiter = 50

      integer l_pi, k_pi, l_ipiv, k_ipiv, lscr, l_scr, k_scr, info
      integer mynri, nprocs, llpi, pilo
      double precision :: eps, dummy
      double precision :: Anorm, Acond, Arnorm, rnorm, ynorm
      double precision :: alpha, beta, beta1, cs, dbar, delta, denom
      double precision :: diag, epsa, epsln, epsr, epsx, gam
      double precision :: gbar, gmax, gmin, oldb, oldeps, grnorm, phi
      double precision :: phibar, rhs1, rhs2, rnorml, rootl, Arnorml
      double precision :: relArnorml, s, sn, t, tnorm2, ynorm2, z
      double precision :: t1, t2, test1, test2, root, bnorm, dxnorm
      double precision :: test3, xnorm, x0norm
      double precision,dimension(n) :: r1, r2, sol, v, w, w1, w2, y

      
      logical main,converged,oactive,parallel_diag,doparallel
      integer iter, istop, i, gy, me, lo ,hi, ld, nmax, nlocal, ithread
      double precision,parameter :: shift = 0.000d0, rtol=1d-5

      converged = .false.
      doparallel = n > 32768

      eps  = epsilon(eps)
      iter = 0
      istop = 0
      Anorm = 0d0
      Acond = 0d0
      ynorm = 0d0
      gmin  = 1d0
      test3 = 0d0

      x0norm = dot_product(x,x)
      bnorm = util_norm2(n,b)

      call dcopy(n,b,1,y,1)
      if (x0norm.gt.0d0) then
        call dsymv('l',n,-1d0,A,n,x,1,1d0,y,1)
        beta1 = util_norm2(n,y)
      else
        beta1 = bnorm
      end if

      ! Solution already converged
      if (beta1.lt.rtol) then
        rnorm = beta1
        istop = -1
        goto 2000
      endif

      ! Continue with MINRES
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

!$omp parallel do simd schedule(static) if (doparallel)
      do i=1,n
        r1(i) = y(i)
        r2(i) = y(i)
        w(i) = 0d0
        w1(i) = 0d0
        w2(i) = 0d0
        sol(i) = 0d0
      end do
!$omp end parallel do simd

      do
        iter = iter + 1
        s = 1d0/beta


        v(:) = s*y(:)

        call dsymv('l',n,1d0,A,n,v,1,0d0,y,1)
        if (iter.ge.2) call daxpy(n,-beta/oldb,r1,1,y,1)

        alpha = dot_product(v,y)
        call daxpy(n,-alpha/beta,r2,1,y,1)
        r1(:) = r2(:)
        r2(:) = y(:)

        oldb = beta
        beta = util_norm2(n,r2)

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

        w1(:) = w2(:)
        w2(:) = w(:)
        w(:) = (v(:) - oldeps*w1(:) - delta*w(:))/gam
        call daxpy(n,phi,w,1,sol,1)
      
        dxnorm = util_norm2(n,sol)
        xnorm = dsqrt(x0norm + dxnorm**2 + 2d0*dot_product(sol,x))

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
        t1 = 1d0 + test1
        t2 = 1d0 + test2

        if (iter  .ge. maxiter  ) istop=6
        if (Acond .ge. 0.1d0/eps) istop=5
        if (t2    .le. 1d0      ) istop=2
        if (t1    .le. 1d0      ) istop=1
        if (epsx  .ge. beta1    ) istop=3
        if (test2 .le. rtol     ) istop=2
        if (test1 .le. rtol     ) istop=1
        if (rnorm .le. rtol     ) istop=8

        if (istop.ne.0) exit
      enddo

      call daxpy(n,1d0,sol,1,x,1)

 1000 continue
 2000 continue

      end subroutine

      subroutine rfp2std(matrix,n)
C
      ! Converts a matrix from RFP packed format, to lower-triangular full
      ! format
      !
      ! Note: when the size of the matrix is even, in reality we are
      ! using a shifted RFP. For example, for N=4 the full format and
      ! RFP formats are given by:
      !
      ! [  a11  ---  ---  ---  ]        [  a33  a34  ]
      ! [  a12  a22  ---  ---  ]        [  a11  a44  ]
      ! [  a13  a23  a33  ---  ]        [  a12  a22  ]
      ! [  a14  a24  a34  a44  ]        [  a13  a23  ]
      !                                 [  a14  a24  ]
      !
      ! The format we are using here is
      !
      ! [  a11  a33  a34 ]
      ! [  a12  a22  a44 ]
      ! [  a13  a23  --- ]
      ! [  a14  a24  --- ]
C
      implicit none
      integer n
      double precision matrix(n,n)
      integer i

      if (mod(n,2).eq.0) then
        do i=n/2+1,n
          matrix(i:,i) = matrix(i-n/2,i-n/2+1:n/2+1)
        enddo
      else
        do i=n/2+2,n
          matrix(i:,i) = matrix(i-n/2-1,i-n/2:n/2+1)
        enddo
      endif
      end subroutine
