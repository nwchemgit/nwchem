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
c$Id: grid_quadv0b.F,v 1.2 2000-03-29 01:59:11 edo Exp $
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
      integer ldmat,idmat,i, grid_nbfm,ntot_cube,iii
      external grid_nbfm
!       write(0,*) ' v0b nrc ',ntot_cube,' nqpt ',nqpts
C
      npol = (ipol*(ipol+1))/2
      acc_AO_gauss = dble(iAOacc)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
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
         else                   ! general case
            max_elem_ang_scr = max_component * 28
         endif 
C
cedo         do iii=1,nqpts
cedo            write(0,'(A,I4,4F16.9)') 
cedo     .           ' RR ',ictr_buf,
cedo     ,           qxyz(1,iii),qxyz(2,iii),qxyz(3,iii),
cedo     .           qwght(iii)
cedo         enddo
c     
c
c           Evaluate the AO basis set at each of the quad. points.
c           allocate arrays for exponents and contraction coefficients
c           within int_eval_gbsets
c           Evaluate the distances (squared) between each center and the points
c
         mbf_ao=nbf_ao
         if(natoms.gt.1) then
         call ifill(mbf_ao, 0, ibf_ao, 1)
         mbf_ao=grid_nbfm(AO_bas_han,  natoms, 
     &        ictr_buf,rad,xyz,qxyz,nqpts,
     &        ibf_ao,  iniz, ifin, expo, 
     &        minexp,ldiff,acc_ao_gauss)
cedo         write(6,*) ' mbf_ao ',mbf_ao
         if (mbf_ao.eq.0) return
         else
            iniz(1)=1
            ifin(1)=nbf_ao
            do i=1,nbf_ao
               ibf_ao(i)=i
            enddo
         endif
c
         if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_ao, 'chi_ao', 
     &        lchi_ao, ichi_ao))
     &        call errquit('grid_quadv0a: cannot allocate chi_ao',0)
c     
         if(GRAD) then
            if (.not.MA_Push_Get(mt_dbl, 3*nqpts*mbf_ao,
     .           'delchi_ao', ldelchi_ao, idelchi_ao))
     &           call errquit('grid_quadv0a: cant alloc delchi_ao',0)
         endif
c     
         nscr = 3*nqpts*nbf_ao_mxcont + 
     &        max(3*nqpts*nbf_ao_mxprim,
     &        nqpts*max_elem_ang_scr) + 1
c     
         if (.not.MA_Push_Get(mt_dbl,nscr,'scrCD',lscr,iscr))
     &        call errquit('grid_quadv0a: cannot allocate scrCD',12)
c     
         call qdist(rchi_atom, rq, qxyz, xyz, nqpts, natoms)
         call new_eval_gbsets(ao_bas_han, grad, .false., 
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &        hess, ibf_ao, 0, rchi_atom, 
     &        rq, qxyz, xyz, nqpts, nqpts, 
     &        mbf_ao, natoms, dbl_mb(iscr), nscr, 
     &        expo, nbf_ao_mxprim, ccoef, ncoef_max, 
     &        iniz,ifin)
c     
         if (.not.ma_pop_stack(lscr))
     &        call errquit('grid_quadv0a: cannot pop stack',0)
         
c     
c        get reduced Xi(r) and dXi(r) over atoms
c     
         
         call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,           dbl_mb(ichi_ao),rchi_atom)
         if(grad) 
     .        call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,           dbl_mb(idelchi_ao),rdelchi_atom)

!         call grid_rchi3(ntot_cube+1,
!     ,        nqpts, natoms,iniz,ifin,
!     ,        dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
!     ,        rchi_atom,rdelchi_atom,
!     ,        grad)
c     
         call xc_rhogen(tol_rho, AO_bas_han, g_dens,
     &        nbf_ao_mxnbf_ce, natoms, ipol, nqpts, nbf_ao,
     &        mbf_ao, ldelley, GRAD, Fmat, Pmat, Bmat, 
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &        ibf_ao, iniz, ifin, rho, delrho, 
     &        bp, rscr, 
     &        rchi_atom, rdelchi_atom,
     &        rdens_atom, cetobfr)
c
c
c           Evaluate the integrated charge density.
c     
         nqtot = nqtot + nqpts
c
         rho0=rho_int
         rho_int = rho_int + 
     &        ddot(nqpts, qwght, 1, rho, 1)
cedo         write(6,'(A,I3,3(2X,1PE12.6))') 
cedo     .        ' rho ',nqpts,rho_int,rho(1)
c
            call xc_eval_fnl(rho, delrho, Amat,
!     &           Amat2, Cmat, Cmat2, nqpts, Ex, Ec, qwght,
     &           dummy, Cmat, dummy, nqpts, Ex, Ec, qwght,
     &           GRAD, .false., dummy, .false.)

            if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &         'Dmat', lDmat,iDmat))
     &         call errquit('grdi_quadv0_a: cannot allocate Dmat',0)
c     
            call xc_tabcd(tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, 
     &           Dbl_MB(iDmat), Dbl_MB(ichi_ao), 
     &           Dbl_MB(idelchi_ao), ipol, nqpts, mbf_ao, 
     &           nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD, 
     &           AO_bas_han, 
     &           natoms, iniz, ifin, g_vxc, ibf_ao,
     &           rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr)


            if (.not.ma_pop_stack(lDmat))
     &           call errquit('grid_quadv0_a: cannot pop stack', 1)
            if(grad) then
               if (.not.ma_pop_stack(ldelchi_ao)) call
     &              errquit('grid_quadv0a: cannot pop stack', 2)
            endif
         if (.not.ma_pop_stack(lchi_ao))
     &        call errquit('grid_quadv0a: cannot pop stack', 3)
      return
      end
