      Subroutine grid_quadv0b(
     ,     ntot_cube,ictr_buf,qxyz, qwght, nqpts,rad,grad,
     .     rtdb, g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess,  xyz,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq, rdens_atom, 
     &                       cntobfr, cetobfr, 
     ,     ibf_ao,bp,rscr,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc)
c
c$Id: grid_quadv0b.F,v 1.7 2000-04-07 16:34:26 edo Exp $
c
      implicit none
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "grid_quad.fh"
      integer  nqpts    ! [in] no. grid points  to be processed
      integer  ictr_buf ! [in] ctr of grid
      double precision acc_AO_gauss, rad, dummy(1),rho0
      integer nbf_ao_mxnbf_ce2, mbf_ao, npol, nscr
      integer lchi_ao, ichi_ao, ldelchi_ao, idelchi_ao, iscr, lscr
      integer ldmat,idmat,i, grid_nbfm,ntot_cube
      integer max_elem_ang_scr, max_component
      double precision dabsmax,wmax
      external dabsmax
      external grid_nbfm
C
      npol = (ipol*(ipol+1))/2
      acc_AO_gauss = dble(iAOacc)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(nqpts.eq.0) return
c
c     compute variables for scratch space in function evaluation
c     assumes cartesian; spherical set would only require 2*MXang + 1
c     
      max_component = (nbf_ao_mxang+1)*(nbf_ao_mxang+2)/2 
c     
c     Allocate scratch space.
c     
      if (nbf_ao_mxang .eq. 0)then
         max_elem_ang_scr = max_component * 3
      elseif (nbf_ao_mxang .le. 3)then
         max_elem_ang_scr = max_component * 9
      else                      ! general case
         max_elem_ang_scr = max_component * 28
      endif 
C     
c     Evaluate the AO basis set at each of the quad. points.
c     allocate arrays for exponents and contraction coefficients
c     within int_eval_gbsets
c     Evaluate the distances (squared) between each center and the points
c     
      mbf_ao=nbf_ao
      if(natoms.gt.1) then
         call ifill(mbf_ao, 0, ibf_ao, 1)
         mbf_ao=grid_nbfm(AO_bas_han,  natoms, 
     &        ictr_buf,rad,xyz,qxyz,nqpts,
     &        ibf_ao,  iniz, ifin, expo, 
     &        minexp,ldiff,acc_ao_gauss)
         if (mbf_ao.eq.0) return
      else
         iniz(1)=1
         ifin(1)=nbf_ao
         do i=1,nbf_ao
            ibf_ao(i)=i
         enddo
      endif
c     
      if(GRAD) then
         if (.not.MA_Push_Get(mt_dbl, 3*nqpts*mbf_ao,
     .        'delchi_ao', ldelchi_ao, idelchi_ao))
     &        call errquit('grid_quadv0a: cant alloc delchi_ao',0)
      endif
      if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_ao, 'chi_ao', 
     &     lchi_ao, ichi_ao))
     &     call errquit('grid_quadv0a: cannot allocate chi_ao',0)
c     
c     
      nscr = 3*nqpts*nbf_ao_mxcont + 
     &     max(3*nqpts*nbf_ao_mxprim,
     &     nqpts*max_elem_ang_scr) + 1
c     
      if (.not.MA_Push_Get(mt_dbl,nscr,'scrCD',lscr,iscr))
     &     call errquit('grid_quadv0a: cannot allocate scrCD',12)
c     
      call qdist(rchi_atom, rq, qxyz, xyz, nqpts, natoms)
      call new_eval_gbsets(ao_bas_han, grad, .false., 
     &     dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &     hess, ibf_ao, 0, rchi_atom, 
     &     rq, qxyz, xyz, nqpts, nqpts, 
     &     mbf_ao, natoms, dbl_mb(iscr), nscr, 
     &     expo, nbf_ao_mxprim, ccoef, ncoef_max, 
     &     iniz,ifin)
c     
      if (.not.ma_pop_stack(lscr))
     &     call errquit('grid_quadv0a: cannot pop stack',0)
c     
c     get reduced Xi(r) and dXi(r) over atoms
c     
      call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,     dbl_mb(ichi_ao),rchi_atom)
      if(grad) 
     .     call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,     dbl_mb(idelchi_ao),rdelchi_atom)
      
      call grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
     ,     rchi_atom,rdelchi_atom,rdens_atom,iniz,
     ,     dbl_mb(ichi_ao),dbl_mb(idelchi_ao),mbf_ao,grad)
      if(nqpts.eq.0) goto 2010
!         call grid_rchi3(ntot_cube+1,
!     ,        nqpts, natoms,iniz,ifin,
!     ,        dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
!     ,        rchi_atom,rdelchi_atom,
!     ,        grad)
c     
      wmax= dabsmax(nqpts,qwght)
      call xc_rhogen(tol_rho, AO_bas_han, g_dens,
     &     nbf_ao_mxnbf_ce, natoms, ipol, nqpts, nbf_ao,
     &     mbf_ao, ldelley, GRAD, Fmat, Pmat, Bmat, 
     &     dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &     ibf_ao, iniz, ifin, rho, delrho, 
     &     bp, rscr, 
     &     rchi_atom, rdelchi_atom,
     &     rdens_atom, cetobfr,wmax)
c     
c     Evaluate the integrated charge density.
c     
      nqtot = nqtot + nqpts
c     
      rho0=rho_int
      rho_int = rho_int + 
     &     ddot(nqpts, qwght, 1, rho, 1)
      if(abs(rho0-rho_int).lt.tol_rho) goto 2010
c     
      call xc_eval_fnl(rho, delrho, Amat,
c     &           Amat2, Cmat, Cmat2, nqpts, Ex, Ec, qwght,
     &     dummy, Cmat, dummy, nqpts, Ex, Ec, qwght,
     &     GRAD, .false., dummy, .false.)
      
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &     'Dmat', lDmat,iDmat))
     &     call errquit('grdi_quadv0_a: cannot allocate Dmat',0)
c     
      call xc_tabcd(tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, 
     &     Dbl_MB(iDmat), Dbl_MB(ichi_ao), 
     &     Dbl_MB(idelchi_ao), ipol, nqpts, mbf_ao, 
     &     nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD, 
     &     AO_bas_han, 
     &     natoms, iniz, ifin, g_vxc, ibf_ao,
     &     rchi_atom, rdelchi_atom,
     &     rdens_atom, cetobfr)
      
      
      if (.not.ma_pop_stack(lDmat))
     &     call errquit('grid_quadv0_a: cannot pop stack', 1)
 2010 continue
      if (.not.ma_chop_stack(lchi_ao))
     &     call errquit('grid_quadv0a: cannot pop stack', 3)
      return
      end
      subroutine grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
     ,     rchi_atom,rdelchi_atom,rdens_atom,iniz,chi,delchi,
     ,     mbf_ao,grad)
      implicit none
#include "mafdecls.fh"
      integer natoms
      integer nqpts ![in/out]
      double precision qwght(*),qxyz(3,*)! [in/out]
      double precision rchi_atom(*),rdelchi_atom(*)
      double precision chi(*),delchi(3,*)
      double precision rdens_atom(natoms,natoms,*)
      integer iniz(*),mbf_ao
      double precision tol_rho
      logical grad
c
      double precision dabsmax
      external dabsmax
c
      integer k_iscr,l_iscr,k_dscr,l_dscr
      integer nzero,i
      double precision pmax,funcmax,ff,delmax
c
      nzero=0
c
c     get pmax ang funcmax
c
      pmax=0d0
      funcmax=dabsmax(natoms,rchi_atom)
      delmax=funcmax
      if(grad) delmax=max(delmax,dabsmax(natoms,rdelchi_atom))
      do i=1,natoms
         if(iniz(i).ne.0) then
             pmax=max(pmax,dabsmax(natoms,rdens_atom(1,i,1)))
          endif
      enddo
      ff=abs(funcmax*delmax*pmax)
      do i=1,nqpts
         if(ff*qwght(i).lt.tol_rho) then
            nzero=nzero+1
            qwght(i)=0d0
         endif
      enddo
      if(nzero.ne.0) then
            if (.not.MA_Push_Get(mt_int,nqpts,
     &         'iscr', l_iscr,k_iscr))
     .        call errquit('grid_pcom: pushget failed',1)
            k_dscr=0
            if(grad) then
            if (.not.MA_Push_Get(mt_dbl,nqpts*3*mbf_ao,
     &         'dscr', l_dscr,k_dscr))
     .        call errquit('grid_pcom: pushget failed',2)
            endif
c        write(0,'(" % ",f12.6,I4,2X,1PE12.6)') 
c     .           (nzero*100d0)/nqpts,nzero,tol_rho
         call grid_pcompress(nqpts,qwght,chi,qxyz,delchi,
     ,           int_mb(k_iscr),dbl_mb(k_dscr),
     ,        mbf_ao,1d-19,grad)
            if (.not.ma_chop_stack(l_iscr))
     &           call errquit('grid_quadv0_a: cannot pop stack', 11)
      endif
      return
      end
      subroutine grid_pcompress(nq,qwght,chi,qxyz,delchi,
     ,     iscr,dscr,
     .     mbf_ao,toll,grad)
      implicit none
      integer nq,mbf_ao
      double precision qwght(*),chi(*),qxyz(3,*),delchi(*)
      double precision dscr(*)
      integer iscr(*)
      logical grad
      double precision toll ![input]
c
      integer nsign,i,q,j
c        
c       compress grid
c
      if(nq.eq.0 )return
      nsign=0
      do i=1,nq
      if(abs(qwght(i)).gt.toll) then
       nsign=nsign+1
       iscr(nsign)=i
       endif
      enddo
      if(nq.eq.nsign) return
      if(nsign.eq.0)  then
         nq=0
         return
      else
      do q=1,nsign
        qxyz(1,q)=qxyz(1,iscr(q))
        qxyz(2,q)=qxyz(2,iscr(q))
        qxyz(3,q)=qxyz(3,iscr(q))
        qwght(q)=qwght(iscr(q))
      enddo
        do j=1,mbf_ao
           do q=1,nsign
              chi(q+(j-1)*nsign)=chi(iscr(q)+(j-1)*nq)
           enddo
        enddo
        if(grad) then
           call dcopy(3*nq*mbf_ao,delchi,1,dscr,1)
           call grid_pcompres1(nsign,nq,mbf_ao,delchi,dscr,iscr)
        endif
        nq=nsign
      endif
      return
      end
      subroutine grid_pcompres1(nsign,nq,mbf_ao,delchi,dscr,iscr)
      implicit none
      integer nsign,nq,mbf_ao
      integer iscr(*)
      double precision delchi(nsign,3,*)
      double precision dscr(nq,3,*)
c
      integer j,q
c
      do j=1,mbf_ao
         do q=1,nsign
            delchi(q,1,j)=dscr(iscr(q),1,j)
            delchi(q,2,j)=dscr(iscr(q),2,j)
            delchi(q,3,j)=dscr(iscr(q),3,j)
         enddo
      enddo
      return
      end
