#if 0
#define QUADV0_REVERT
#else
#undef QUADV0_REVERT
#endif
      Subroutine grid_quadv0b(
     ,     ictr_buf,qxyz, qwght, nqpts,rad,grad,
     .     g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess,  xyz,iscr,nscr,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq, rdens_atom,  cetobfr, 
     ,     ibf_ao,bas_cent_info,bas_cset_info,iscratch,docset,
     ,     rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,amat2,cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat,
     ,     donewbf,l3d,triplet)
c
c$Id: grid_quadv0b.F,v 1.31 2002-03-22 02:37:37 sohirata Exp $
c
      implicit none
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "grid_quad.fh"
#include "grid_sic.fh"
c !!! temp !!! ***
#include "rtdb.fh"
#include "bgj.fh"
c !!! temp !!! ***
      integer  nqpts    ! [in] no. grid points  to be processed
      integer  ictr_buf ! [in] ctr of grid
      integer bas_cent_info(3,natoms),bas_cset_info(6,*)
      logical docset(*)
      integer iscratch(*)
      double precision acc_AO_gauss, rad, dummy(1),rho0
      integer nbf_ao_mxnbf_ce2, mbf_ao, npol, ncontrset, maxdbas
      integer lchi_ao, ichi_ao, ldelchi_ao, idelchi_ao
      integer ldmat,idmat,i, grid_nbfm
      double precision dabsmax,wmax
      external dabsmax
      external grid_nbfm
      integer iscf_rho,iscf_delrho
      logical grid_sicinit,grid_sicend
      external grid_sicinit,grid_sicend
      double precision acc_xc_gauss
      logical donewbf           !!! !!!! *** temp
      logical l3d
      logical triplet
C
      npol = (ipol*(ipol+1))/2
      acc_AO_gauss = dble(iAOacc)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(nqpts.eq.0) return
      if (GRAD) then
         maxdbas = 1
      else
         maxdbas = 0
      endif
c
c     Evaluate the AO basis set at each of the quad. points.
c     allocate arrays for exponents and contraction coefficients
c     within int_eval_gbsets
c     Evaluate the distances (squared) between each center and the points
c     
      if(xcfit) acc_xc_gauss=dble(ixcacc)
      mbf_ao=nbf_ao
      if(natoms.gt.1) then
         call ifill(mbf_ao, 0, ibf_ao, 1)
         mbf_ao=grid_nbfm(AO_bas_han,  natoms, 
     &        ictr_buf,rad,xyz,qxyz,nqpts,
     &        ibf_ao, docset, iniz, ifin, expo, 
     &        minexp,ldiff,acc_ao_gauss,iatype_pt_chg)
         if (mbf_ao.eq.0) return
      else
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in grid_quadv0b',0)
         iniz(1)=1
         ifin(1)=nbf_ao
         do i=1,ncontrset
            docset(i)=.true.
         enddo
         do i=1,nbf_ao
            ibf_ao(i)=i
         enddo
      endif
c     
      if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_ao, 'chi_ao', 
     &     lchi_ao, ichi_ao))
     &     call errquit('grid_quadv0b: cannot allocate chi_ao',0)
      if(GRAD) then
         if (.not.MA_Push_Get(mt_dbl, 3*nqpts*mbf_ao,
     .        'delchi_ao', ldelchi_ao, idelchi_ao))
     &        call errquit('grid_quadv0b: cant alloc delchi_ao',0)
      endif
      if (sic_orb_index.eq.1) then
         if(.not.grid_sicinit(nqpts,
     ,        l_vect1,i_vect1, ldelrho_sig,idelrho_sig))
     ,        call errquit(' grid_quadv0b: grid_sicinit failed',0)
      endif
c     
c     
      call qdist(rchi_atom, rq, qxyz, xyz, nqpts, natoms)
      if (donewbf) then
      call xc_eval_basis(ao_bas_han, maxdbas, dbl_mb(ichi_ao),
     &     dbl_mb(idelchi_ao), 0d0, 0d0, rq, qxyz, xyz, nqpts, natoms,
     &     iniz, docset, bas_cent_info, bas_cset_info)
      else
c !!! Can delete this and donewbf logic once confident in new basis fn code
      call new_eval_gbsets(ao_bas_han, grad, .false., 
     &     dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &     hess, ibf_ao, 0, rchi_atom, 
     &     rq, qxyz, xyz, nqpts, nqpts, 
     &     mbf_ao, natoms, iscr, nscr, 
     &     expo, nbf_ao_mxprim, ccoef, ncoef_max, 
     &     iniz,ifin)
      endif
c     
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
     ,     dbl_mb(ichi_ao),dbl_mb(idelchi_ao),mbf_ao,grad,
     ,     ipol,iscratch)
      if(nqpts.eq.0) goto 2010

      wmax= dabsmax(nqpts,qwght)
      if(do_scf) then
         if (sic_orb_index.eq.1) then
            call grid_sicrho(nqpts,mbf_ao,cetobfr,wmax,
     ,           rho,delrho,grad,ibf_ao, iniz, ifin,
     ,           dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     ,           rchi_atom, rdelchi_atom,
     ,           rdens_atom,
     ,           Fmat, Pmat, Bmat,
     ,           dbl_mb(i_vect1),dbl_mb(idelrho_sig))
         else
            call xc_rhogen(tol_rho, AO_bas_han, g_dens,
     &           nbf_ao_mxnbf_ce, natoms, ipol, nqpts, nbf_ao,
     &           mbf_ao, ldelley, GRAD, Fmat, Pmat, Bmat, 
     &           dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &           ibf_ao, iniz, ifin, rho, delrho, 
     &           rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr,wmax)
         endif
      else  if (do_cpks_l) then
#ifdef QUADV0_REVERT
            call xc_rhogen(tol_rho, AO_bas_han, g_dens(2),
     &           nbf_ao_mxnbf_ce, natoms, ipol, nqpts, nbf_ao,
     &           mbf_ao, ldelley, GRAD, Fmat, Pmat, Bmat, 
     &           dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &           ibf_ao, iniz, ifin, rho, delrho, 
     &           rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr,1d0)
#else
            call xc_rhogen_pert(tol_rho, AO_bas_han, g_dens,
     &           nbf_ao_mxnbf_ce, natoms, nmat+1,ipol, nqpts, 
     .           nbf_ao, mbf_ao, GRAD, Fmat, Pmat, Bmat, 
     &           dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     &           ibf_ao, iniz, ifin, rho, delrho, 
     &           rchi_atom, 
     .           rdelchi_atom, rdens_atom, cetobfr,wmax)
#endif
            else if (do_cpks_r) then
               call errquit('gridquadv0b: what do do for cpks rhs?',0)
            endif
c     
c     Evaluate the integrated charge density.
c     
      nqtot = nqtot + nqpts
c     
      rho0=rho_int
      rho_int = rho_int + 
     &     ddot(nqpts, qwght, 1, rho, 1)
      if(do_scf.and.abs(rho0-rho_int).lt.tol_rho) goto 2010
c     
c
c           Evaluate the XC potential and energy sampling matrices.
c           Set up pointer to the SCF density for the CPS LHS case.
c
            if (do_cpks_l) then
#ifdef QUADV0_REVERT
               iscf_rho = 1
               iscf_delrho = 1
#else
               iscf_rho = 1 + nqpts*ipol*nmat
               iscf_delrho = 1 + nqpts*3*ipol*nmat
#endif
            else
               iscf_rho = 1
               iscf_delrho = 1
            endif
      if(do_cpks_l.and.ipol.eq.2)then
c     copy beta to 3
         call dcopy(nqpts, rho(iscf_rho+nqpts), 1, 
     &        rho(iscf_rho+2*nqpts), 1)
c     copy alpha to 2
         call dcopy(nqpts, rho(iscf_rho), 1, 
     &        rho(iscf_rho+nqpts), 1)
c     add beta to 1 to get alpha+beta
         call daxpy(nqpts, 1.d0, rho(iscf_rho+2*nqpts), 1,
     &        rho(iscf_rho), 1)
      endif     
      call xc_eval_fnl(rho(iscf_rho), delrho(iscf_delrho), Amat,
     &     Amat2, Cmat, Cmat2, nqpts, Ex, Ec, qwght,
     &     GRAD, .false., dummy, do_cpks_l .or. do_cpks_r)

      if (sic_orb_index.eq.1) then
         call grid_sicpot(nqpts,qxyz,qwght,rho,delrho,amat,cmat,
     ,     dbl_mb(idelrho_sig),dbl_mb(i_vect1))
      endif

      if(xcfit) then
         call grid_xcfit(xc_bas_han,nbf_xc,tmat,grad,
     ,     iniz,ifin,ibf_ao,ccoef,expo,nscr,iscr,hess,
     ,     amat,cmat,rq,rchi_atom,nbf_xc_mxprim, ncoef_max,
     ,     xyz,qxyz,nqpts,natoms,ipol,ictr_buf,rad,
     ,     acc_xc_gauss)
      else
         if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &        'Dmat', lDmat,iDmat))
     &        call errquit('grid_quadv0b: cannot allocate Dmat',0)
c     
         if(do_scf) then
            call xc_tabcd(l3d,tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, 
     &           Dbl_MB(iDmat), Dbl_MB(ichi_ao), 
     &           Dbl_MB(idelchi_ao), ipol, nqpts, mbf_ao, 
     &           nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD, 
     &           AO_bas_han, 
     &           natoms, iniz, ifin, g_vxc, ibf_ao,
     &           rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr)
         else if (do_cpks_l) then
#ifdef QUADV0_REVERT
            call xc_tabcd(l3d,tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, 
     &           Dbl_MB(iDmat), Dbl_MB(ichi_ao), 
     &           Dbl_MB(idelchi_ao), ipol, nqpts, mbf_ao, 
     &           nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD, 
     &           AO_bas_han, 
     &           natoms, iniz, ifin, g_vxc, ibf_ao,
     &           rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr)
#else
c
c     Combine derivs of functional with deriv densities to get coeffs
c     of bfn pairs in CPKS LHS - these are written back into rho and
c     delrho
c
            call xc_cpks_coeff(rho, delrho, Amat2, Cmat2, Cmat,
     &           delrho(iscf_delrho), nmat, ipol, nqpts, GRAD,
     &           triplet)
            call xc_tabcd_cpks_lhs(tol_rho, Fmat, Pmat, rho, Bmat,
     &          delrho, Dbl_MB(iDmat), Dbl_MB(ichi_ao), 
     &          Dbl_MB(idelchi_ao), nmat, ipol, nqpts, mbf_ao,
     &          nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD, AO_bas_han,
     &          natoms, iniz, ifin, g_vxc, ibf_ao,
     &          rchi_atom, rdelchi_atom,
     &          rdens_atom, cetobfr)
#endif
         else if (do_cpks_r) then
            call errquit('gridquadv0b: what do do for cpks rhs?',0)
         endif

         if (.not.ma_pop_stack(lDmat))
     &      call errquit('grid_quadv0b: cannot pop stack', 1)
      endif
 2010 continue
      if (sic_orb_index.eq.1) then
         if(.not.grid_sicend(l_vect1,ldelrho_sig))
     ,        call errquit(' grid_quadv0b: grid_sicend failed',0)
      endif
      if(grad)then
         if (.not.ma_pop_stack(ldelchi_ao))
     &        call errquit('grid_quadv0b: cannot pop stack', 2)
      endif
      if (.not.ma_pop_stack(lchi_ao))
     &     call errquit('grid_quadv0b: cannot pop stack', 3)
      return
      end
