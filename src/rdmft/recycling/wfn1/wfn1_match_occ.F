C> \ingroup wfn1
C> @{
C>
C> \brief Fit occupation functions to match occupation numbers
C>
C> Assume we have a set of orthonormal occupation functions
C> \f$\{\phi_p(x)\}\f$, we now want to find a set of occupation
C> functions \f$\{\phi'_p(x)\}\f$ that minimize
C> \f{eqnarray*}{
C>   R&=&
C>   \int\left(\sum_{p'}\phi'_{p'}^2(x)-d(x)\right)^2 dx
C> \f}
C> where \f$p'\f$ means that the sum runs over just the occupied 
C> occupation functions. In addition the new occupation functions can
C> be expressed in terms of the old ones as
C> \f{eqnarray*}{
C>   \phi'_p &=& \sum_q \phi_q U_{qp}
C> \f}
C> where \f$U\f$ is a unitary matrix. We will now establish \f$U\f$
C> by minimizing \f$R\f$ as a function of \f$U\f$.
C>
C> First the gradient of \f$R\f$ wrt. \f$U\f$ is
C> \f{eqnarray*}{
C>   \frac{dR}{dU} &=&
C>   \int 2\left(\sum_{p'}\phi'^2_{p'}(x)-d(x)\right)
C>   \frac{d\sum_{q'}\phi'^2_{q'}(x)}{dU} dx\\\\
C>   &=& \int 4\left(\sum_{p'}\phi'^2_{p'}(x)-d(x)\right)
C>   \sum_{q'}\phi'_{q'}(x)\frac{d\phi'_{q'}(x)}{dU} dx\\\\
C> \f}
C> Substituting \f$\phi'_{q'}\f$ and differentiating explicitly wrt.
C> \f$U_{uw'}\f$ we get
C> \f{eqnarray*}{
C>   \frac{dR}{dU_{uw'}} 
C>   &=& \int 4\left(\sum_{p'}\phi'^2_{p'}(x)-d(x)\right)
C>   \sum_{q'rs}\phi_r(x)U_{rq'}\frac{d\phi_s(x)U_{sq'}}{dU_{uw'}} dx\\\\
C>   &=& \int 4\left(\sum_{p'}\phi'^2_{p'}(x)-d(x)\right)
C>   \sum_{r}\phi_r(x) U_{rw'}\phi_u(x) dx\\\\
C> \f}
C> Substituting \f$\phi'_{p'}\f$ we obtain
C> \f{eqnarray*}{
C>   \frac{dR}{dU_{uw'}} 
C>   &=& \int 4\left(\sum_{p'st}\phi_{s}(x)U_{sp'}\phi_{t}(x)U_{tp'}
C>    -  d(x)\right)
C>   \sum_{r}\phi_r(x) U_{rw'}\phi_u(x) dx\\\\
C> \f}
C> Rearranging
C> \f{eqnarray*}{
C>   \frac{dR}{dU_{uw'}} 
C>   &=& \int 4\left(\sum_{p'st}\phi_{s}(x)U_{sp'}U^T_{p't}\phi_{t}(x)
C>    -  d(x)\right)
C>   \sum_{r}\phi_r(x) U_{rw'}\phi_u(x) dx\\\\
C>   &=& \int 4\left(\sum_{p'st}\phi_{s}(x)U_{sp'}U^T_{p't}\phi_{t}(x)
C>    -  d(x)\right)
C>   \phi_u(x)\sum_{r}\phi_r(x) U_{rw'} dx\\\\
C>   &=& \int 4\sum_{p'st}\phi_{s}(x)U_{sp'}U^T_{p't}\phi_{t}(x)
C>       \phi_u(x)\sum_{r}\phi_r(x) U_{rw'} dx\\\\
C>   &-& \int 4d(x)\phi_u(x)\sum_{r}\phi_r(x) U_{rw'} dx\\\\
C>   &=& 0
C> \f}
C> If we write
C> \f$w(x) = 4\sum_{p'st}\phi_{s}(x)U_{sp'}U^T_{p't}\phi_{t}(x)\f$
C> then we get
C> \f{eqnarray*}{
C>   \sum_r \int w(x)\phi_u(x)\phi_r(x)dx U_{rw'}
C>   &=& \sum_r \int d(x)\phi_u(x)\phi_r(x)dx U_{rw'}
C> \f}
C> Next we can write \f$A_{ur}=\int w(x)\phi_u(x)\phi_r(x)\mathrm{d}x\f$
C> and \f$B_{ur}=\int d(x)\phi_u(x)\phi_r(x)\mathrm{d}x\f$ to summarize
C> the equation above as
C> \f{eqnarray*}{
C>   \sum_r A_{ur}U_{rw'} &=& \sum_r B_{ur}U_{rw'}
C> \f}
C> Modifying this equation a little bit we get
C> \f{eqnarray*}{
C>   \sum_r (A_{ur}-B_{ur})U_{rw'} &=& \lambda\sum_r \delta_{ur}U_{rw'}
C> \f}
C> and we can solve this problem as an eigenvalue equation. The 
C> approach has to be iterative as \f$A\f$ depends on the solution.
C> Also at convergence all eigenvalues should be zero.
C>
      subroutine wfn1_match_occ(newof,nbf,ne,docc,r_of,r_of2,aa,bb,uu,
     +                          bt,m1,m2,ww,v1,oprint)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      logical newof  !< [Input] If .TRUE. start with fresh occupation
                     !< functions, otherwise start with what is there
                     !< already
      integer nbf    !< [Input] The number of basis functions
      integer ne     !< [Input] The number of electrons
      logical oprint !< [Input] If .TRUE. print some progress data
c
      double precision docc(nbf)      !< [Input] The occupation numbers
      !< to match
      double precision r_of(nbf,nbf)  !< [In/Output] The occupation
      !< functions stored in columns, i.e.
      !< \f$\mathrm{r\_of}(i,p) = C^O_{ip}\f$
      double precision r_of2(nbf,nbf) !< [Scratch] Space to store
      !< \f$\sum_P \mathrm{r\_of}(i,p)U_{pq}\f$
c
      double precision aa(nbf,nbf) !< [Scratch] matrix \f$A\f$
      double precision bb(nbf,nbf) !< [Scratch] matrix \f$B\f$
      double precision uu(nbf,nbf) !< [Scratch] totalmatrix \f$U\f$
      double precision bt(nbf,nbf) !< [Scratch] matrix \f$U^T B U\f$
      double precision ww(nbf)     !< [Scratch] vector \f$w\f$
      double precision m1(nbf,nbf) !< [Scratch] matrix
      double precision m2(nbf,nbf) !< [Scratch] matrix
      double precision v1(nbf)     !< [Scratch] vector
c
c     Local
c
      integer ii    !< Counter over natural orbitals
      integer jj    !< Counter over natural orbitals
      integer ip    !< Counter over occupation functions
      integer iq    !< Counter over occupation functions
      integer iter  !< Iteration counter
      integer lwork !< The work space size for LAPACK
      integer l_w   !< The memory handle for LAPACK work space
      integer k_w   !< The memory offset for LAPACK work space
      integer info  !< The LAPACK error number
      integer imax  !< Index of maximum value
c
      double precision s    !< Overlap
      double precision dnrm !< The infinite order norm of the gradient
      double precision derr !< The error in the occupation numbers
      double precision dmax !< Maximum value
c
      double precision sqrt12 !< The square root of 1/2
      double precision sqrt13 !< The square root of 1/3
c
c     Functions
c
      double precision ddot, wfn1_ddot3
      external         ddot, wfn1_ddot3
c
c     Code
c
      info  = 0
      lwork = -1
      sqrt12 = sqrt(1.0d0/2.0d0)
      sqrt13 = sqrt(1.0d0/3.0d0)
c
c     Find out how much memory LAPACK needs and allocate corresponding
c     chunks
c
c     call dsyev('v','u',nbf,aa,nbf,v1,s,lwork,info)
c     lwork = s + 0.25d0
c     if (.not.ma_push_get(MT_DBL,lwork,'work dsygv',l_w,k_w))
c    +   call errquit("wfn1_match_occ: alloc work dsygv failed",
c    +                ma_sizeof(MT_DBL,lwork,MT_BYTE),MA_ERR)
c
cDEBUG
c     write(*,*)"wfn1_match_occ: newof = ",newof
cDEBUG
      if (newof) then
c
c       Generate a new set of orthonormal vectors
c       All vectors need to interlink to guarantee that the solution can
c       be reached.
c
        call dfill(nbf*nbf,0.0d0,r_of,1)
c       do ii = 1, nbf
c         r_of(ii,ii) = 1.0d0
c         r_of(ii,1)  = 1.0d0
c       enddo
        if (nbf.eq.1) then
          r_of(1,1) = 1.0d0
        else if (nbf.gt.1) then
          r_of(1    ,1  ) = +sqrt12
          r_of(2    ,1  ) = +sqrt12
          r_of(nbf-1,nbf) = -sqrt12
          r_of(nbf  ,nbf) = +sqrt12
          do ii = 2, nbf-1
            r_of(ii-1,ii) = -sqrt13
            r_of(ii  ,ii) = +sqrt13
            r_of(ii+1,ii) = +sqrt13
          enddo
        else
          call errquit("wfn1_match_occ: invalid number of basis "//
     +                 "functions",nbf,UERR)
        endif
c
      endif
c
c     Orthonormalize what ever is currently in r_of (just in case)
c
      do ip = 1, nbf
        s = ddot(nbf,r_of(1,ip),1,r_of(1,ip),1)
        if (s.lt.1.0d-15) then
          call errquit("wfn1_match_occ: normalization catastrophy",
     +                  ip,UERR)
        endif
        s = 1.0d0/sqrt(s)
        call dscal(nbf,s,r_of(1,ip),1)
      enddo
      do ip = 2, nbf
        iter = 1
        s = 2.0d0
        do while (s.gt.sqrt(2.0d0))
          do iq = 1, ip-1
            s =  ddot(nbf,r_of(1,iq),1,r_of(1,ip),1)
            call daxpy(nbf,-s,r_of(1,iq),1,r_of(1,ip),1)
          enddo
          s = ddot(nbf,r_of(1,ip),1,r_of(1,ip),1)
          if (s.lt.0.5d0.and.iter.ge.2) then
            call errquit("wfn1_match_occ: orthonormalization "//
     +                   "catastrophy",ip,UERR)
          endif
          s = 1.0d0/sqrt(s)
          call dscal(nbf,s,r_of(1,ip),1)
          iter = iter + 1
        enddo
      enddo
c
c     Initially U is going to be the unit matrix
c
      call dfill(nbf*nbf,0.0d0,uu,1)
      do ii = 1, nbf
        uu(ii,ii) = 1.0d0
      enddo
c
c     Construct matrix B in the basis of input occupation functions
c
c     do iq = 1, nbf
c       do ip = 1, nbf
c         bb(ip,iq) = 0.0d0
c         do ii = 1, nbf
c           do jj = 1, nbf
c             bb(ip,iq) = bb(ip,iq) + sqrt(docc(ii)*docc(jj))
c    +                              * r_of(ii,ip)*r_of(jj,iq)
c           enddo
c         enddo
c       enddo
c     enddo
c
      dnrm = 10.0d0
      derr = 10.0d0
      iter = 0
      do while ((dnrm.gt.1.0d-12).and.
     +          (iter.lt.200000))
        iter = iter + 1
c
c       Use current U to calculate W
c
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_of,nbf,uu,nbf,0.0d0,
     +             r_of2,nbf)
        do ii = 1, nbf
          ww(ii) = 0.0d0
          do ip = 1, ne
            ww(ii) = ww(ii) + r_of2(ii,ip)**2
          enddo
        enddo
c
c       Calculate current error in matching the occupation numbers
c
        derr = 0.0d0
        do ii = 1, nbf
          derr = derr + (ww(ii)-docc(ii))**2
          v1(ii) = ww(ii)-docc(ii)
        enddo
        derr = sqrt(derr)
cDEBUG
c       write(*,*)iter,' derr = ',derr,dnrm
        if (mod(iter,1000).eq.0.or..false.) then
          write(*,*)
          do ii = 1, nbf
            write(*,*)iter,' occ: ',ii,docc(ii),ww(ii)
          enddo
          write(*,*)iter,' derr = ',derr
        endif
cDEBUG
c
c       Construct skew symmetric matrix A in the basis of current best
c       occupation functions
c
        call dfill(nbf*nbf,0.0d0,aa,1)
        do iq = 1, nbf
          do ip = 1, iq-1
            s = wfn1_ddot3(nbf,v1,1,r_of2(1,ip),1,r_of2(1,iq),1)
            aa(ip,iq) =  s
            aa(iq,ip) = -s
          enddo
        enddo
c
c       Compute the rotation exp(A)
c
        dnrm = 0.0d0
        do ip = 1, nbf
          do iq = ip+1, nbf
            dnrm = max(dnrm,abs(aa(iq,ip)))
          enddo
        enddo
cDEBUG
        if (dnrm.gt.1.0d-20) then
          s = derr/dnrm
c         s = 1.0d0/sqrt(sqrt(dnrm))
          call dscal(nbf*nbf,s,aa,1)
        endif
c       write(*,*)'*** matrix A:'
c       call hess_hssout(aa,nbf,nbf,nbf)
cDEBUG
        call wfn1_exp(nbf,aa,m1,m2,bb)
c       call dsyev('v','u',nbf,aa,nbf,v1,dbl_mb(k_w),lwork,info)
c       if (info.lt.0) then
c         call errquit("wfn1_match_occ: dsyev failed: invalid data in"
c    +               //" argument",-info,UERR)
c       else if (info.gt.0) then
c         call errquit("wfn1_match_occ: dsyev failed to converge: no."
c    +               //" non-zero off-diag elements remain",info,UERR)
c       endif
cDEBUG
c       do ii=1, nbf
c         write(*,*)iter," :eigval: ",ii,v1(ii)
c       enddo
c       write(*,*)'*** vectors:'
c       call hess_hssout(aa,nbf,nbf,nbf)
cDEBUG
c
c       Apply this rotation to the U we have so far
c
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,uu,nbf,bb,nbf,0.0d0,m1,nbf)
        call dcopy(nbf*nbf,m1,1,uu,1)
c
c       Use current U to calculate overlap of W_i with D_occ
c
c       call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_of,nbf,uu,nbf,0.0d0,
c    +             r_of2,nbf)
c       do ii = 1, nbf
c         ww(ii) = 0.0d0
c       enddo
c       do ip = 1, nbf
c         do ii = 1, nbf
c           ww(ip) = ww(ip) + (r_of2(ii,ip)**2)*docc(ii)
c         enddo
c       enddo
c
c       Transform BT with these vectors and find the ne largest
c       eigenvalues and put the corresponding vectors at the front
c
c       call dgemm('t','n',nbf,nbf,nbf,1.0d0,aa,nbf,bt,nbf,0.0d0,m1,nbf)
c       call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,aa,nbf,0.0d0,m2,nbf)
cDEBUG
c       do ii=1, nbf
c         write(*,*)iter," :eigval bt: ",ii,m2(ii,ii)
c       enddo
cDEBUG
c       do ip = 1, ne
c         dmax = m2(ip,ip)
c         imax = ip
c         do iq = 2, nbf
c           if (m2(iq,iq).gt.dmax) then
c             dmax = m2(iq,iq)
c             imax = iq
c           endif
c         enddo
c         if (dmax.gt.m2(ip,ip)) then
c
c           Swap columns in U
c
c           do ii = 1, nbf
c             s           = aa(ii,imax)
c             aa(ii,imax) = aa(ii,ip)
c             aa(ii,ip)   = s
c           enddo
c
c           Swap entries on the diagonal of m2
c
c           s             = m2(imax,imax)
c           m2(imax,imax) = m2(ip,ip)
c           m2(ip,ip)     = s
c         endif
c       enddo
c       do ip = 1, ne
c         dmax = ww(ip)
c         imax = ip
c         do iq = 2, nbf
c           if (ww(iq).gt.dmax) then
c             dmax = ww(iq)
c             imax = iq
c           endif
c         enddo
c         if (dmax.gt.ww(ip)) then
c
c           Swap columns in U
c
c           do ii = 1, nbf
c             s           = uu(ii,imax)
c             uu(ii,imax) = uu(ii,ip)
c             uu(ii,ip)   = s
c           enddo
c
c           Swap entries on the diagonal of ww
c
c           s        = ww(imax)
c           ww(imax) = ww(ip)
c           ww(ip)   = s
c         endif
c       enddo
c
        if (oprint.and.(mod(iter,1000).eq.0)) then
          write(luout,*)'wfn1_match_occ: ',iter,dnrm,derr
        endif
c
      enddo
c
c     Construct the final vectors by transforming the original ones
c     with the optimized transformation
c
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_of,nbf,uu,nbf,0.0d0,
     +           r_of2,nbf)
      call dcopy(nbf*nbf,r_of2,1,r_of,1)
c
      end
C>
C> \brief Dot-product with non-unit weights
C>
C> For the optimization of the occupation functions a dot-product is
C> needed that includes the current error in the occupation numbers
C> as the weight function. This function implements this 3-vector
C> dot-product.
C>
C> The weighted dot-product between 2 functions \f$\phi_a\f$ and
C> \f$\phi_b\f$ is defined as
C> \f{eqnarray*}{
C>   S_{ab} &=& \int w(x)\phi_a(x)\phi_b(x) \mathrm{d}x \\\\
C>          &=& \sum_i w_i \phi_{ia}\phi_{ib}
C> \f}
C>
      double precision function wfn1_ddot3(n,w,iw,x,ix,y,iy)
      implicit none
c
      integer n  !< [Input] The length of the vectors
      integer iw !< [Input] The increment of the weight vector
      integer ix !< [Input] The increment of the `x` vector
      integer iy !< [Input] The increment of the `y` vector
c
      double precision w(1+(n-1)*iw) !< [Input] The weight vector
      double precision x(1+(n-1)*ix) !< [Input] The `x` vector
      double precision y(1+(n-1)*iy) !< [Input] The `y` vector
c
      integer i
      double precision result
c
      wfn1_ddot3 = 0.0d0
      if (n.le.0) return
c
      result = w(1)*x(1)*y(1)
      do i = 1, n-1
        result = result + w(1+i*iw)*x(1+i*ix)*y(1+i*iy)
      enddo
c
      wfn1_ddot3 = result
c
      return
      end
C> 
C> @}
