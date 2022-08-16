      subroutine grid_quadv0b(
     &     ictr_buf,qxyz, qwght, nqpts,rad,grad,
     &     g_dens, ncoef_max, natoms, 
     &     rho_int, nqtot, 
     &     iniz, ifin, rho, delrho, 
     &     hess,  xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq, rdens_atom,  cetobfr, 
     &     ibf_ao,bas_cent_info,bas_cset_info,iscratch,docset,
     &     rchi_atom,rdelchi_atom,
     &     Ex, Ec, Amat, Cmat, g_vxc, Amat2, Cmat2,
     &     nmat, do_scf, do_cpks_l, do_cpks_r, 
     &     do_nmrks_r, do_xdm, tmat,
     &     l3d,triplet,
     &     bas_cent_xinfo,bas_cset_xinfo,
     &     docsetx, 
     &     ttau, kske, Mmat, Mmat2,
     &     kslap, Lmat,
     &     StericEnergy,
     &     do_gxc, Amat3, Cmat3,
     Q     qbecke,
cc AJL/Begin/FDE
     &     fde_option, g_dens_fde, ncenters_fde, rdens_atom_fde,
     &     Fmat_fde, Pmat_fde, Bmat_fde, ttau_fde, rho_fde, delrho_fde,
     &     cetobfr_fde, Amat_fde, Amat2_fde, Amat3_fde,
     &     Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     &     ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &     iniz_fde, ifin_fde, ibf_fde, xyz_fde,
     &     rchi_atom_fde, rdelchi_atom_fde, rq_fde,
     &     bas_cent_info_fde, bas_cset_info_fde, docset_fde,rho_int_fde)
cc AJL/End
c
      implicit none
#include "errquit.fh"
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "grid_quad.fh"
#include "grid_sic.fh"
#include "dftps.fh"
#include "geom.fh"
c
      integer  nqpts    ! [in] no. grid points  to be processed
      integer  ictr_buf ! [in] ctr of grid
      integer bas_cent_info(3,natoms),bas_cset_info(6,*)
      logical docset(*)
      integer bas_cent_xinfo(3,natoms),bas_cset_xinfo(6,*)
      logical docsetx(*)
      integer iscratch(*)
      double precision acc_AO_gauss, rad, rho0
      integer nbf_ao_mxnbf_ce2, mbf_ao, npol, ncontrset, maxdbas
      integer lchi_ao, ichi_ao, ldelchi_ao, idelchi_ao
      integer lheschi_ao, iheschi_ao
      integer lDmat,iDmat,i, grid_nbfm,llap,ilap
      integer lGmat,iGmat
      double precision dabsmax,wmax
      external dabsmax
      external grid_nbfm
      integer iscf_rho,iscf_delrho
      integer iscf_tau,iscf_ttau
      logical grid_sicinit,grid_sicend,xc_chkdelrq
      external grid_sicinit,grid_sicend,xc_chkdelrq
      double precision acc_xc_gauss
      logical l3d
      logical triplet, kske, kslap
      integer lemat,iemat,lfmat,ifmat,k_scr,l_scr
      logical do_2nd
      double precision StericEnergy
      logical do_3rd
      double precision qbecke(*)
      double precision lapdum
c     AOR begin
      integer nxdm
      integer ixdm_v, ixdm_ml, lxdm_v, lxdm_ml
      common /xdmd/ nxdm, ixdm_v, lxdm_v, ixdm_ml, lxdm_ml
c     AOR end
      logical oprint_qbecke
cc AJL/Begin/FDE
c
      logical docset_fde(*)
      double precision bas_cent_info_fde(3,ncenters_fde)
      double precision bas_cset_info_fde(6,*)

      integer lchi_fde, ichi_fde
      integer ldelchi_fde, idelchi_fde
      integer lheschi_fde, iheschi_fde

      integer llap_fde, ilap_fde
      integer lDmat_fde, iDmat_fde

      character*32 pname
      integer mbf_fde
      integer ncontrset_fde

      double precision Ex_qm, Ex_fde, Ex_tot
      double precision Ec_qm, Ec_fde, Ec_tot
      double precision StericEnergy_qm, StericEnergy_fde
      double precision StericEnergy_tot

      double precision ddot
      external ddot

      integer basis,mxnbf_ce,mbf_xc

      if (odftps) call pstat_on(ps_qv0b)

      pname = 'grid_quadv0b: '
      mbf_fde = 0
      Ex_qm = 0.0d0
      Ex_fde = 0.0d0
      Ex_tot = 0.0d0
      Ec_qm = 0.0d0
      Ec_fde = 0.0d0
      Ec_tot = 0.0d0
      StericEnergy_qm = 0.0d0
      StericEnergy_fde = 0.0d0
      StericEnergy_tot = 0.0d0

      if (adft) then
        basis = xc_bas_han
        mbf_ao = nbf_xc
        mxnbf_ce = nbf_xc_mxnbf_ce
      else
        basis = ao_bas_han
        mbf_ao = nbf_ao
        mxnbf_ce = nbf_ao_mxnbf_ce
      endif

      npol = (ipol*(ipol+1))/2
      acc_AO_gauss = dble(iAOacc)
      nbf_ao_mxnbf_ce2 = mxnbf_ce*mxnbf_ce
      if(nqpts.eq.0) then
        if (odftps) call pstat_off(ps_qv0b)
        return
      endif
cc AJL/Begin/FDE
      if (kslap.or.do_xdm) then
         maxdbas = 2
      elseif (grad) then
         maxdbas = 1
         ilap = 0
         iheschi_ao = 0
         iheschi_fde = 0
      else
         maxdbas = 0
         ilap = 0
         idelchi_ao = 0
         iheschi_ao = 0
         idelchi_fde = 0
         iheschi_fde = 0
      endif
cc AJL/End
c
c     becke charges
c
      oprint_qbecke = util_print('becke charges', print_high)
c
c     Evaluate the AO basis set at each of the quad. points.
c     allocate arrays for exponents and contraction coefficients
c     within int_eval_gbsets
c     Evaluate the distances (squared) between each center and the points
c     
      if((.not.adft).and.xcfit) acc_xc_gauss=dble(ixcacc)

      if(natoms.gt.1) then
         call icopy(mbf_ao, 0,0, ibf_ao, 1)  
cc AJL/Begin/FDE
c         if (fde_option.lt.1) then
c No need to pre-screen basis functions if we are calculating ESP
         if (fde_option.ne.1) then
           mbf_ao=grid_nbfm(basis,  natoms, 
     &            ictr_buf,rad,xyz,qxyz,nqpts,
     &            ibf_ao, docset, iniz, ifin, expo, 
     &            minexp,ldiff,acc_ao_gauss,iatype_pt_chg)
c         else if (fde_option.eq.1) then
         else
c Calculating the ESP over the grid:
           if (.not.bas_numcont(AO_bas_han, ncontrset))
     &       call errquit('Exiting in grid_quadv0b',0, BASIS_ERR)
           do i=1,natoms
             iniz(i)=1
           enddo
c           ifin(1)=nbf_ao
           do i=1,ncontrset
             docset(i)=.true.
           enddo
c           do i=1,nbf_ao
c              ibf_ao(i)=i
c           enddo
         endif
cc AJL/End
         if (mbf_ao.eq.0) then
           if (odftps) call pstat_off(ps_qv0b)
           return
         endif
      else
c We've only got one atom
      if (.not.bas_numcont(basis, ncontrset))
     &     call errquit('Exiting in grid_quadv0b',0, BASIS_ERR)
         iniz(1)=1
         ifin(1)=nbf_ao
         if (adft) ifin(1) = nbf_xc
         do i=1,ncontrset
            docset(i)=.true.
         enddo
         do i=1,ifin(1)
            ibf_ao(i)=i
         enddo
      endif

cc AJL/Begin/FDE
c Perform the same process of screening relevant basis functions 
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then
        mbf_fde=nbf_fde

        if (ncenters_fde.gt.1) then
          call icopy(mbf_fde, 0,0, ibf_fde, 1)
          mbf_fde=grid_nbfm(FD_bas_han,ncenters_fde,
     &            ictr_buf,rad,xyz_fde,qxyz,nqpts,
     &            ibf_fde, docset_fde, 
     &            iniz_fde, ifin_fde, expo,
     &        minexp_fde,ldiff_fde,
     &        acc_ao_gauss,iatype_pt_chg_fde)
        else
          if (.not.bas_numcont(FD_bas_han, ncontrset_fde))
     &      call errquit('Exiting in grid_quadv0b',0, BASIS_ERR)
          iniz_fde(1)=1
          ifin_fde(1)=nbf_fde
          do i=1,ncontrset_fde
             docset_fde(i)=.true.
          enddo
          do i=1,nbf_fde
             ibf_fde(i)=i
          enddo
        endif

        if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_fde, 'chi_fde',
     &           lchi_fde, ichi_fde))
     &           call errquit(pname//'chi_fde',0, MA_ERR)
        if( grad ) then
           if (.not.MA_Push_Get(mt_dbl, 3*nqpts*mbf_fde,
     &        'delchi_fde', ldelchi_fde, idelchi_fde)) then
           call errquit(pname//'cant alloc delchi_fde',0,0)
           endif
        endif
c XDM and SIC are disabled for FDE.... needs further testing
c        if (do_xdm) then
c           if (.not.MA_Push_Get(mt_dbl, 6*nqpts*mbf_fde,
c     &           'heschi_fde', lheschi_fde, iheschi_fde)) then
c               call errquit(pname//'cant alloc heschi_fde',0,0)
c           endif
c        endif
      endif ! fde_option.eq.(3 or 4 or 6 or 7)
cc AJL/End

      if (.not.MA_Push_Get(mt_dbl, nqpts*mbf_ao, 'chi_ao', 
     &     lchi_ao, ichi_ao))
     &     call errquit(pname//'cannot allocate chi_ao',0,
     &       MA_ERR)
      if( grad.or.do_xdm ) then
         if (.not.MA_Push_Get(mt_dbl, 3*nqpts*mbf_ao,
     &      'delchi_ao', ldelchi_ao, idelchi_ao)) then
           call errquit(pname//'cant alloc delchi_ao',0,0)
        endif
      endif
      if (kslap.or.do_xdm) then
         if (.not.MA_Push_Get(mt_dbl, 6*nqpts*mbf_ao,
     &        'heschi_ao', lheschi_ao, iheschi_ao)) then
            call errquit(pname//'cant alloc heschi_ao',0,0)
         endif
      endif
      if (sic_orb_index.eq.1) then
         if(.not.grid_sicinit(nqpts,
     &        l_vect1,i_vect1, ldelrho_sig,idelrho_sig))
     &        call errquit(pname//'grid_sicinit failed',0,
     &       CALC_ERR)
      endif

      if (odftps) call pstat_on(ps_ebf)
      call qdist(rchi_atom, rq, qxyz, xyz, nqpts, natoms)
      call xc_eval_basis(basis, maxdbas, dbl_mb(ichi_ao),
     &     dbl_mb(idelchi_ao), dbl_mb(iheschi_ao), 0d0, rq, 
     Q     qxyz, xyz, nqpts, natoms,
     &     iniz, docset, bas_cent_info, bas_cset_info)

cc AJL/BEGIN/FDE
c Evaluate the basis over the grid for FDE
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then
c
        call qdist(rchi_atom_fde, rq_fde,
     &           qxyz, xyz_fde, nqpts, ncenters_fde)
        call xc_eval_basis(FD_bas_han, maxdbas, dbl_mb(ichi_fde),
     &     dbl_mb(idelchi_fde), dbl_mb(iheschi_fde), 0d0, 
     &     rq_fde, qxyz, xyz_fde, nqpts, 
     &     ncenters_fde, iniz_fde, docset_fde, 
     &     bas_cent_info_fde, bas_cset_info_fde)
c
      endif !fde_option
cc AJL/End
c
      if (odftps) call pstat_off(ps_ebf)
c
C AJL/BEGIN/FDE
      if (fde_option.eq.1) then
c Calculate potential at each point, embedding potential.
c Outsourced this to a subroutine: 
        call dft_frozemb_fde_esp(nqpts,qxyz,qwght,
     &       mbf_ao,dbl_mb(ichi_ao),g_dens_fde,g_vxc)
c We're finished so lets jump to the final deallocates
        goto 2020
      endif ! fde_option = 1
c AJL/END
c
c     get reduced Xi(r) and dXi(r) over atoms
c
      call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,     dbl_mb(ichi_ao),rchi_atom)
      if( grad ) 
     .     call util_rmsatbf(3*nqpts, natoms,iniz,ifin,
     ,     dbl_mb(idelchi_ao),rdelchi_atom)
cc AJL/Begin/FDE
c We pre-screen grid points to see if they meet a threshold only if the
c calculation involves addressing the QM density only. For the FDE and
c combined systems we cannot screen in this same way as it will results in
c losing useful grid points
      if (fde_option.eq.0.or.fde_option.eq.2.or.fde_option.eq.5) then
c      call grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
           if (odftps) call pstat_on(ps_fitv)
           call grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
     ,     rchi_atom,rdelchi_atom,rdens_atom,iniz,
     ,     dbl_mb(ichi_ao),dbl_mb(idelchi_ao),dbl_mb(iheschi_ao),
     M     mbf_ao,grad,kslap,
     ,     ipol,iscratch,adft)
           if (odftps) call pstat_off(ps_fitv)
      endif
      if (.not.MA_Push_Get(mt_dbl,mxnbf_ce*nqpts*3,
     &     'Dmat', lDmat,iDmat))
     &     call errquit(pname//'cannot allocate Dmat',0,
     &    MA_ERR)
      if (.not.MA_Push_get(mt_dbl,mxnbf_ce*nqpts,
     &    'Gmat',lGmat,iGmat))
     &     call errquit(pname//'cannot allocate Gmat',0,
     &    MA_ERR)
c This could probably be turned in to an "else" without harm
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then
c
c     get reduced Xi(r) and dXi(r) over atoms
c
        call util_rmsatbf(nqpts, ncenters_fde,iniz_fde,
     &       ifin_fde, dbl_mb(ichi_fde),rchi_atom_fde)
        if( grad )
     .       call util_rmsatbf(3*nqpts,ncenters_fde,iniz_fde,
     .       ifin_fde, dbl_mb(idelchi_fde), rdelchi_atom_fde)
c Not sure if want to pre-screen with the FDE basis - otherwise might lose
c valuable points from the QM environment. Only applied for FDE calcs.
c For combined systems we cannot use this method of removing grid points
        if (fde_option.eq.3.or.fde_option.eq.6)
     ,       call grid_pscreen(nqpts,qwght,qxyz,tol_rho,ncenters_fde,
     ,       rchi_atom_fde,rdelchi_atom_fde,
     ,       rdens_atom_fde,iniz_fde,
     ,       dbl_mb(ichi_fde),dbl_mb(idelchi_fde),dbl_mb(iheschi_fde),
     ,       mbf_fde,grad,do_xdm,
     ,       ipol,iscratch,adft)

        if (.not.MA_Push_Get(mt_dbl,nbf_fde_mxnbf_ce*nqpts*3,
     &       'Dmat_fde', lDmat_fde,iDmat_fde))
     &       call errquit(pname//'cannot allocate Dmat',0,
     &       MA_ERR)
       endif !fde_option
c AJL/End

      if(nqpts.eq.0) goto 2010

      if (kslap.or.do_xdm) then
         if (.not.MA_Push_Get(mt_dbl,nqpts*ipol*(ipol+1)/2,
     &        'lap', llap,ilap))
     &        call errquit('grid_quadv0b: cannot allocate lap',0,
     &        MA_ERR)
      endif

      wmax= dabsmax(nqpts,qwght)
      if(do_scf.or.do_xdm.or.do_nmrks_r) then
         if (odftps) call pstat_on(ps_xcrho)
            if (.not.adft) then
              call xc_rhogen(0,
     T          tol_rho,AO_bas_han,g_dens,nbf_ao_mxnbf_ce, 
     N          natoms, 1,1,1,
     I          ipol, nqpts, nbf_ao, mbf_ao, grad, (ipol*(ipol+1))/2,
     F          Fmat, Pmat, Bmat, dbl_mb(iDmat),
     &          dbl_mb(ichi_ao), dbl_mb(idelchi_ao), dbl_mb(iheschi_ao),
     &          ibf_ao, iniz, ifin, rho, delrho, dbl_mb(ilap),
     &          rchi_atom, rdelchi_atom,
     &          rdens_atom, cetobfr,wmax, 
     &          ttau, kske, kslap.or.do_xdm)
              if(xc_chkdelrq()) call
     X             xc_delrhodotr(nqpts,ipol,qxyz,rho,delrho,ttau)
            else
              call xc_rhogen_adft(tol_rho, natoms, ipol, nqpts, 
     .          nbf_xc, mbf_ao, mxnbf_ce, grad, npol, dbl_mb(ichi_ao),
     .          dbl_mb(idelchi_ao), dbl_mb(iheschi_ao), ibf_ao, iniz,
     .          ifin, rho, delrho, dbl_mb(ilap), rchi_atom,
     .          rdelchi_atom, cetobfr, wmax, kslap, tmat)
            endif
         if (odftps) call pstat_off(ps_xcrho)
         if (sic_orb_index.eq.1) then
            call grid_sicrho(nqpts,mbf_ao,cetobfr,wmax,
     ,        rho,delrho,grad,ibf_ao, iniz, ifin,
     ,        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     ,        rchi_atom, rdelchi_atom,
     ,        rdens_atom,
     ,        Fmat, Pmat, Bmat,
     ,        dbl_mb(i_vect1),dbl_mb(idelrho_sig),
     &        ttau, kske)
         endif
      else  if (do_cpks_l.or.do_gxc) then
         if (odftps) call pstat_on(ps_xcrho)
            call xc_rhogen(1,
     T        tol_rho, AO_bas_han, g_dens, nbf_ao_mxnbf_ce, 
     N        natoms, 1,1, nmat+1,
     I        ipol, nqpts, nbf_ao, mbf_ao, grad, ipol,
     F        Fmat, Pmat, Bmat, dbl_mb(iDmat),
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 0d0,
     &        ibf_ao, iniz, ifin, rho, delrho, lapdum,
     &        rchi_atom, rdelchi_atom, 
     &        rdens_atom, cetobfr,wmax, 
     &        ttau, kske, .false.)
         if (odftps) call pstat_off(ps_xcrho)
      else if (do_cpks_r) then
       call errquit(pname//'what do do for cpks rhs?',0, CALC_ERR)
      endif

c AJL/Begin
c Generate density for FDE on QM grid
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then

c Currently disabled. To be tested
c        if (do_xdm) then
c          if (.not.MA_Push_Get(mt_dbl,nqpts*ipol*(ipol+1)/2,
c     &        'lap', llap_fde, ilap_fde))
c     &        call errquit(pname//'cannot allocate lap',0,
c     &        MA_ERR)
c        endif

c wmax is already defined further above
c        wmax= dabsmax(nqpts,qwght)
        if(do_scf.or.do_xdm.or.do_nmrks_r) then
           if (odftps) call pstat_on(ps_xcrho)
              call xc_rhogen(0,
     T          tol_rho,FD_bas_han,g_dens_fde,nbf_fde_mxnbf_ce,
     N          ncenters_fde, 1,1,1,
     I          ipol, nqpts, nbf_fde, mbf_fde, grad, (ipol*(ipol+1))/2,
     F          Fmat_fde, Pmat_fde, Bmat_fde, dbl_mb(iDmat_fde),
     &        dbl_mb(ichi_fde), dbl_mb(idelchi_fde),dbl_mb(iheschi_fde),
c     &          dbl_mb(ichi_fde), dbl_mb(idelchi_fde), 0d0,
     &          ibf_fde, iniz_fde, ifin_fde,
     &          rho_fde, delrho_fde, dbl_mb(ilap_fde),
     &          rchi_atom_fde, rdelchi_atom_fde,
     &          rdens_atom_fde, cetobfr_fde, wmax,
     &          ttau_fde, kske, do_xdm)
           if (odftps) call pstat_off(ps_xcrho)
c Also currently disabled. To be tested 
c           if (sic_orb_index.eq.1) then
c              call grid_sicrho(nqpts,mbf_ao,cetobfr,wmax,
c     ,          rho,delrho,grad,ibf_ao, iniz, ifin,
c     ,          dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
c     ,          rchi_atom, rdelchi_atom,
c     ,          rdens_atom,
c     ,          Fmat, Pmat, Bmat,
c     ,          dbl_mb(i_vect1),dbl_mb(idelrho_sig),
c     &          ttau, kske)
c           endif
c      else  if (do_cpks_l.or.do_gxc) then
c             call xc_rhogen(1,
c    T          tol_rho, FD_bas_han, g_dens_fde, nbf_fde_mxnbf_ce,
c    N          ncenters_fde, 1,1, nmat+1,
c    I          ipol, nqpts, nbf_fde, mbf_fde, grad, ipol,
c    F          Fmat_fde, Pmat_fde, Bmat_fde, dbl_mb(iDmat_fde),
c    &          dbl_mb(ichi_fde), dbl_mb(idelchi_fde), 0d0,
c    &          ibf_fde, iniz_fde, ifin_fde, 
c    &          rho_fde, delrho_fde, lapdum,
c    &          rchi_atom_fde, rdelchi_atom_fde,
c    &          rdens_atom_fde, cetobfr_fde, wmax,
c   &          ttau_fde, kske, .false.)
        else  if (do_cpks_l.or.do_cpks_r.or.do_gxc) then
c        else if (do_cpks_r) then
         call 
     &   errquit(pname//'what do do for cpks in fde? Not
     &            implemented',0, CALC_ERR)
        endif !do_scf
c
      endif !fde_option
c
cc AJL/End
c     
c     Evaluate the integrated charge density.
c   
  
      nqtot = nqtot + nqpts
c     
      rho0=rho_int
      rho_int = rho_int + ddot(nqpts, qwght, 1, rho, 1)

      if(oprint_qbecke) then
          qbecke(ictr_buf)=qbecke(ictr_buf)+(rho_int-rho0)
      endif

      if(do_scf.and.abs(rho0-rho_int).lt.tol_rho) goto 2010
c
c     Evaluate the XC potential and energy sampling matrices.
c     Set up pointer to the SCF density for the CPS LHS case.
c
      iscf_rho = 1
      iscf_delrho = 1
      if (do_cpks_l.or.do_gxc) then
          iscf_rho = 1 + nqpts*ipol*nmat
          iscf_delrho = 1 + nqpts*3*ipol*nmat
      endif
c
      if((do_cpks_l.or.do_gxc).and.ipol.eq.2)then
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
c
      if (odftps) call pstat_on(ps_excf)
      do_2nd = (do_cpks_l .or. do_cpks_r)

cc AJL/Comment
c XDM is disabled for FDE
      if (do_xdm) then
         call xc_eval_xdm(rho,delrho,dbl_mb(ilap),nqpts,
     &        qxyz,qwght,ttau, natoms,xyz,
     &        dbl_mb(ixdm_v),dbl_mb(ixdm_ml))
      else

        do_3rd = do_gxc  ! need 3rd derivatives ?
c Daniel (2-7-13): The last false in this call is used to initiallize
c ldew2, which is used for TDDFT gradients to prevent weighting of the
c XC-derivative matrices.
c
cc AJL/Begin/FDE
c Evaluate the energy from the specific functional/density
        if (fde_option.eq.0) then
c Standard assessment of QM density using QM settings
          call xc_eval_fnl(rho(iscf_rho), delrho(iscf_delrho), 
     &         Amat, Amat2, Cmat, Cmat2, nqpts, Ex, Ec, qwght,
     &         grad, .false., 0d0, do_2nd, 
     &         ttau, kske, Mmat, Mmat2,
     &         dbl_mb(ilap), kslap, Lmat,
     &         StericEnergy,
     &         do_3rd, Amat3, Cmat3, .false.)
        elseif (fde_option.eq.2.or.fde_option.eq.4) then
c Assessment of QM contribution for non-Additive exchange-correlation
          call xc_eval_fnl_fde(rho(iscf_rho), delrho(iscf_delrho),
     &       Amat, Amat2, Cmat, Cmat2, nqpts, Ex_qm, Ec_qm, qwght,
     &       grad, .false., 0d0, do_2nd,
     &       ttau, kske, Mmat, Mmat2,
     &       StericEnergy_qm,
     &       do_3rd, Amat3, Cmat3, .false.)
        elseif (fde_option.eq.5.or.fde_option.eq.7) then
c Kinetic Energy for QM contribution (FDE)
          call ts_eval_fnl(rho(iscf_rho), delrho(iscf_delrho),
     &       Amat, Amat2, Cmat, Cmat2, nqpts, Ex_qm, qwght,
     &       grad, .false., 0d0, do_2nd, ttau, kske, Mmat, Mmat2,
     &       do_3rd, Amat3, Cmat3, .false.)
        endif !fde_option

        if (fde_option.eq.3.or.fde_option.eq.4) then
cc Assessment of MM contribution for non-Additive exchange-correlation
          call xc_eval_fnl_fde(rho_fde(iscf_rho), 
     &       delrho_fde(iscf_delrho),
     &       Amat_fde, Amat2_fde, Cmat_fde, Cmat2_fde, nqpts, Ex_fde,
     &       Ec_fde, qwght,
     &       grad, .false., 0d0, do_2nd,
     &       ttau_fde, kske, Mmat_fde, Mmat2_fde,
     &       StericEnergy_fde,
     &       do_3rd, Amat3_fde, Cmat3_fde, .false.)
        elseif (fde_option.eq.6.or.fde_option.eq.7) then
c Assessment of MM contribution for non-Additive Kinetic Energy 
          call ts_eval_fnl(rho_fde(iscf_rho), delrho_fde(iscf_delrho),
     &       Amat_fde, Amat2_fde, Cmat_fde, Cmat2_fde, nqpts, Ex_fde,
     &       qwght, grad, .false., 0d0, do_2nd,
     &       ttau_fde, kske, Mmat_fde, Mmat2_fde,
     &       do_3rd, Amat3_fde, Cmat3_fde, .false.)
        endif !fde_option

        if (fde_option.eq.4.or.fde_option.eq.7) then
c Copy QM values in to the QM/FDE arrays
          call dcopy(npol*nqpts, rho(iscf_rho), 1,
     &               rho_qmfde(iscf_rho), 1)
          call daxpy(npol*nqpts, 1.d0, rho_fde(iscf_rho), 1,
     &                rho_qmfde(iscf_rho), 1)
c
          if (grad) then
            call dcopy(3*ipol*nqpts, delrho(iscf_delrho), 1, 
     &                 delrho_qmfde(iscf_delrho), 1)
            call daxpy(3*ipol*nqpts, 1.d0, delrho_fde(iscf_delrho), 1,
     &                delrho_qmfde(iscf_delrho), 1)
          endif
c
          if (kske) then
            call dcopy(len_ttau, ttau, 1, ttau_qmfde, 1)
            call daxpy(len_ttau, 1.d0, ttau_fde, 1, ttau_qmfde, 1)
          endif
c
cc Assessment of QM/FDE contribution for non-Additive exchange-correlation
          if (fde_option.eq.4) then
          call xc_eval_fnl_fde(rho_qmfde(iscf_rho), 
     &       delrho_qmfde(iscf_delrho),
     &       Amat_fde, Amat2_fde, Cmat_fde, Cmat2_fde, nqpts, 
     &       Ex_tot, Ec_tot, qwght,
     &       grad, .false., 0d0, do_2nd,
     &       ttau_qmfde, kske, Mmat_fde, Mmat2_fde,
     &       StericEnergy_tot,
     &       do_3rd, Amat3_fde, Cmat3_fde, .false.)
          elseif (fde_option.eq.7) then
c Kinetic Energy of combined system
          call ts_eval_fnl(rho_qmfde(iscf_rho), 
     &       delrho_qmfde(iscf_delrho), 
     &       Amat_fde, Amat2_fde, Cmat_fde, Cmat2_fde, nqpts, 
     &       Ex_tot, qwght, grad, .false., 0d0, do_2nd,
     &       ttau_qmfde, kske, Mmat_fde, Mmat2_fde,
     &       do_3rd, Amat3_fde, Cmat3_fde, .false.)
          endif !fnl_eval
        endif !fde_option

c Collect all the energies correctly.
        if (fde_option.le.1) then
c For normal calculations continue, and ESP calculation as well
          continue
        else if (fde_option.eq.2.or.fde_option.eq.5) then
c QM contributions using FDE variables
          Ex = Ex + Ex_qm
          Ec = Ec + Ec_qm
          StericEnergy = StericEnergy + StericEnergy_qm
        else if (fde_option.eq.3.or.fde_option.eq.6) then
c FDE contributions using FDE variables
          Ex = Ex + Ex_fde
          Ec = Ec + Ec_fde
          StericEnergy = StericEnergy + StericEnergy_fde
          rho_int_fde = rho_int_fde + ddot(nqpts, qwght, 1, rho_fde, 1)
        else if (fde_option.eq.4.or.fde_option.eq.7) then
c Non-additive contributions using FDE variables
!          write(*,*) 'BEFORE: Ex, Ex_tot, Ex_qm, Ex_fde', 
!     &                Ex, Ex_tot, Ex_qm, Ex_fde
          Ex = Ex + (Ex_tot - Ex_qm - Ex_fde)
          Ec = Ec + (Ec_tot - Ec_qm - Ec_fde)
          StericEnergy = StericEnergy + (StericEnergy_tot -
     &                   StericEnergy_qm - StericEnergy_fde)
c          write(*,*) 'AFTER: Ex, Ex_tot, Ex_qm, Ex_fde',
c     &                Ex, Ex_tot, Ex_qm, Ex_fde
c          call util_flush(luout)
         rho_int_fde = rho_int_fde + ddot(nqpts, qwght, 1, rho_qmfde, 1)
        endif

! Fortunately, for the potentials the derivative of fde with respect to
! qm is 0, and so we only need the V_xc(p_tot) and V_xc(p_qm) terms

      endif !do_xdm
cc AJL/End

      if (odftps) call pstat_off(ps_excf)
c
      if (sic_orb_index.eq.1) then
         call grid_sicpot(nqpts,qxyz,qwght,rho,delrho,Amat,cmat,
     ,     dbl_mb(idelrho_sig),dbl_mb(i_vect1))
      endif
c
      if(adft) then
         if (odftps) call pstat_on(ps_tabcd)
         call xc_tabcd_adft(tol_rho,amat,cmat,lmat,dbl_mb(ichi_ao),
     $   dbl_mb(idelchi_ao),dbl_mb(iheschi_ao),
     $   nmat,ipol,nqpts,nbf_xc,mxnbf_ce,grad,xc_bas_han,natoms,iniz,
     $   ifin,ibf_ao,rchi_atom,rdelchi_atom,cetobfr,kslap,
     $   tmat(nbf_xc*ipol+1))
         if (odftps) call pstat_off(ps_tabcd)
      elseif(xcfit) then
         call grid_xcfit(xc_bas_han,nbf_xc,tmat,grad,
     ,     iniz,ifin,ibf_ao,ccoef,expo,nscr,iscr,hess,
     ,     amat,cmat,rq,rchi_atom,nbf_xc_mxprim, ncoef_max,
     ,     xyz,qxyz,nqpts,natoms,ipol,ictr_buf,rad,
     .     bas_cent_xinfo,bas_cset_xinfo,docsetx,
     ,     acc_xc_gauss)
      else
c     
         if(do_scf) then
c
            if (odftps) call pstat_on(ps_tabcd)

cc AJL/Begin/FDE
c Evaluate potentials for QM contributions on QM grid
            if (fde_option.eq.0
     &          .or.fde_option.eq.2.or.fde_option.eq.4.
     &          .or.fde_option.eq.5.or.fde_option.eq.7) then
              call xc_tabcd(0,l3d,
     &        tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, Dbl_MB(iDmat), 
     &        0d0,0d0,0d0,0d0,
     &        Dbl_MB(ichi_ao), dbl_MB(idelchi_ao), dbl_MB(iheschi_ao),
     &        0d0,0,1,
     &        ipol, nqpts, mbf_ao, nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, 
     &        grad, AO_bas_han, 
     &        natoms, iniz, ifin, g_vxc, ibf_ao,
     &        rchi_atom, rdelchi_atom,
     &        rdens_atom, cetobfr,
     &        kske,Mmat,
     &        kslap,dbl_MB(iGmat),Lmat,0d0)
            endif

c The derivative of FDE Energy with repsect to QM density is zero, so we
c only need to do anything here for fde if we have the qm/fde system
            if (fde_option.eq.3.or.fde_option.eq.6) then
c By definition the below value should be zero - this is just a test to see if things work
              continue
            else if (fde_option.eq.4.or.fde_option.eq.7) then
c            if (fde_option.eq.4.or.fde_option.eq.7) then

c We will use the _fde arrays where possible, as well as the already
c declared _qmfde array structures

              call xc_tabcd(0,l3d,
     &          tol_rho, Fmat, Pmat, Amat_fde, Bmat,
     &          Cmat_fde, Dbl_MB(iDmat), 0d0,0d0,0d0,0d0,
     &          Dbl_MB(ichi_ao), dbl_MB(idelchi_ao), 0d0,
     &          0d0,0,1,
     &          ipol, nqpts, mbf_ao, nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2,
     &          grad, AO_bas_han,
     &          natoms, iniz, ifin, g_fde, ibf_ao,
     &          rchi_atom, rdelchi_atom,
     &          rdens_atom, cetobfr,
     &          kske,Mmat_fde,
     &          .false.,0d0,0d0,0d0)

c To get the non-additive potential one needs to subtract G_fde from G_vxc,
c this is currently included in dft_scf as we can't manipulate GAs here.
            endif
c AJL/End

            if (odftps) call pstat_off(ps_tabcd)
c
         else if (do_nmrks_r) then
            if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts,
     &           'iEmat', lEmat,iEmat))
     &           call errquit(pname//'cannot allocate iEmat',0,
     &          MA_ERR)
            if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &           'iFmat', lFmat,iFmat))
     &           call errquit(pname//'cannot allocate iFmat',0,
     &          MA_ERR)
            if (.not.MA_Push_Get(mt_dbl,nqpts,
     &           'iFmat', l_scr,k_scr))
     &           call errquit(pname//'cannot allocate iFmat',0,
     &          MA_ERR)        

            call xc_tabcd(3,.true.,
     T           tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, Dbl_MB(iDmat), 
     E           Dbl_MB(iEmat), Dbl_MB(iFmat),  qxyz,xyz,
     &           Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     N           0d0,0,1,
     I           ipol,nqpts,mbf_ao,nbf_ao_mxnbf_ce,nbf_ao_mxnbf_ce2, 
     G           grad, AO_bas_han,natoms, iniz, ifin, 
     G           g_vxc, ibf_ao, rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr,
     &           kske, Mmat, 
     &           .false.,0d0,0d0, dbl_mb(k_scr))

c
            if (.not.ma_chop_stack(lEmat))
     &         call errquit(pname//'cannot chop Emat', 1, MA_ERR)
         else if (do_cpks_l) then
c
c     Combine derivs of functional with deriv densities to get coeffs
c     of bfn pairs in CPKS LHS --> these are written back into rho, delrho, ttau
c
            if (odftps) call pstat_on(ps_cpkscoef)
            call xc_cpks_coeff(rho, delrho, ttau, 
     &           Amat2, Cmat2, Cmat, Mmat2,
     &           delrho(iscf_delrho), nmat, ipol, nqpts, grad, kske, 
     &           triplet)
            if (odftps) call pstat_off(ps_cpkscoef)
c
            if (odftps) call pstat_on(ps_tabcd)
            call xc_tabcd(1,l3d,
     T           tol_rho,Fmat,Pmat, rho, Bmat,delrho, Dbl_MB(iDmat), 
     N           0d0,0d0,0d0,0d0,
     C           Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     N           0d0,0,nmat, 
     I           ipol,nqpts,mbf_ao,nbf_ao_mxnbf_ce,nbf_ao_mxnbf_ce2, 
     G           grad, AO_bas_han, natoms, iniz, ifin, 
     G           g_vxc, ibf_ao, rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr,
     &           kske, ttau,
     &           .false.,0d0,0d0,0d0)
            if (odftps) call pstat_off(ps_tabcd)
c
c Daniel (1-11-13): Need specific control over contraction for the 
c XC-third derivatives
         else if (do_gxc) then
c Daniel (1-14-13): This routine builds the XC-coefficient matrices,
c which are products of the 3rd derivatives and the perturbed density
c from (X+Y) for the gradients.
            call xc_3rd_deriv(rho, delrho,
     1        Amat3, Cmat3, Cmat2, delrho(iscf_delrho), nmat, ipol,
     2        nqpts, grad, triplet)
c Daniel (1-14-13): Build the XC-third derivative matrix (contracted
c with perturbed densities).
            call xc_tabcd(1, l3d,
     1        tol_rho, Fmat, Pmat, rho, Bmat, delrho, Dbl_MB(iDmat),
     2        0d0, 0d0, 0d0, 0d0,
     3        Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     4        0d0, 0, nmat,
     5        ipol, nqpts, mbf_ao, nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2,
     6        grad, AO_bas_han, natoms, iniz, ifin,
     7        g_vxc, ibf_ao, rchi_atom, rdelchi_atom,
     8        rdens_atom, cetobfr,
     9        kske, ttau,
     &        .false.,0d0,0d0,0d0)
         else if (do_cpks_r) then
            call errquit(pname//'what do do for cpks rhs?',0,
     &       CALC_ERR)
         endif
c
      endif

 2010 continue

      if (.not.ma_chop_stack(lGmat))
     &      call errquit(pname//'cannot chop Gmat', 1, MA_ERR)
      if (.not.ma_chop_stack(lDmat))
     &      call errquit(pname//'cannot chop Dmat', 1, MA_ERR)

 2020 continue

      if (sic_orb_index.eq.1) then
         if(.not.grid_sicend(l_vect1,ldelrho_sig))
     ,        call errquit(pname//'grid_sicend failed',0,
     &       CALC_ERR)
      endif

      if(kslap.or.do_xdm)then
         if(.not.ma_chop_stack(lheschi_ao))
     &       call errquit(pname//'cannot chop heschi_ao', 2, MA_ERR)
      endif

      if(grad)then
         if (.not.ma_chop_stack(ldelchi_ao))
     &        call errquit(pname//'cannot chop delchi_ao', 2, MA_ERR)
      endif

      if (.not.ma_pop_stack(lchi_ao))
     &     call errquit(pname//'cannot pop chi_ao', 3, MA_ERR)

CC AJL/BEGIN/FDE
c
c     clear memory of FDE allocations
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then
c        if (grad) then
c           if(.not.MA_pop_stack(ldelchi_fde))
c     &        call errquit(pname//"pop ldelchi_fde", 1, MA_ERR)
c        endif
        if(.not.MA_chop_stack(lchi_fde))
     &     call errquit(pname//"chop lchi_fde", 1, MA_ERR)
      endif

C
cC AJL/End
      if (odftps) call pstat_off(ps_qv0b)

      return
      end
c
c     function grid_nbfm
c
      integer function grid_nbfm(basis_hand,   mcenters, 
     &     ctr_quad,  rad_quad, xyz, qxyz,nq,
     &                  ibf, docset, iniz, ifin, zprim,
     .     minexp,ldiff,acc_gauss,iatype_pt_chg)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
      integer basis_hand
      integer mcenters ! [input]
      double precision acc_gauss ! [input]
      double precision xyz(3,*) ! [input]
      double precision rad_quad ! [input]
      integer nq ! [in] # grid pts
      double precision qxyz(3,nq) ! [input] coord grid pts
      integer ctr_quad ! grid center  [input]
      integer ibf(*) ! [output]
      logical docset(*) ! [output]
      integer ldiff(*) ! [in]
      double precision minexp(*) ! [in]
      logical iatype_pt_chg(*) ! [in]
c
#include "bas.fh"
c
c     Distance Squared between Sampling Points and Centers
c
      double precision r_q0 ! min dist
      integer  iniz(mcenters),ifin(mcenters)
      double precision zprim(*)
      integer n1, icset, ictr,  nprim, ncontr, 
     &        isphere,  l, iprimo,npt,nshbf
      double precision zmin,acc_loc,x,y,z
      integer n,ic1,ic2,m,ibfhi,ibflo
      double precision alpha,logeps
      double precision gaussian_range,r2,r_arg
      integer icontr,l0(2),l_in
      logical veryfar
      gaussian_range(n,logeps,alpha) = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
      if(acc_gauss.gt.25d0.or.acc_gauss.lt.3d0) call errquit(
     R ' grid_nbfm: silly accgauss ',nint(acc_gauss), UNKNOWN_ERR)
      acc_loc=-acc_gauss
c
      n1 = 0
      npt=0
      do 400 ictr=1,mcenters
        if(iatype_pt_chg(ictr)) goto 400
        if (.not.bas_ce2cnr(basis_hand,ictr,ic1,ic2))
     &       call errquit('Exiting in xc_signf.',11, BASIS_ERR)

        r2=0d0
        if (ictr.ne.ctr_quad) then
          x = xyz(1,ctr_quad) - xyz(1,ictr)
          y = xyz(2,ctr_quad) - xyz(2,ictr)
          z = xyz(3,ctr_quad) - xyz(3,ictr)
          r2 = sqrt(x*x + y*y + z*z)
        endif
        r_arg=0d0
        if (rad_quad.lt.r2) r_arg = (r2-rad_quad)
c
c       check on most diffuse fn
c
        veryfar=r_arg.gt.gaussian_range(ldiff(ictr),acc_loc,
     .       minexp(ictr)) 
        if(veryfar) then
           if(.not.bas_ce2bfr(basis_hand, ictr, ibflo, ibfhi))
     &          call errquit('Exiting in grid_nbfm',4, BASIS_ERR)
           nshbf=ibfhi-ibflo+1
           npt=npt+nshbf
        else
    
           r_q0=1d10
           do n=1,nq
             x = qxyz(1,n) - xyz(1,ictr)
             y = qxyz(2,n) - xyz(2,ictr)
             z = qxyz(3,n) - xyz(3,ictr)
             r_q0 = min(r_q0,sqrt(x*x + y*y + z*z))
           enddo
           do icset = ic1,ic2
              docset(icset) = .false.
c         
c             get info about current contraction set
c       
              if (.not. bas_continfo(basis_hand, icset,  l_in ,nprim,
     &            ncontr, isphere))
     &            call errquit('Exiting in grid_nbfm.',5, BASIS_ERR)
c             fix for sps
              if(ncontr.eq.2) then
                 l0(1)=0
                 l0(2)=1
              else
                 l0(1)=l_in
              endif
c         
c       
c             get exponents and contraction coefficients for this
c             contraction set
c       
              if (.not.bas_get_exponent(basis_hand, icset, zprim))
     &             call errquit('Exiting in grid_nbfm.',7, BASIS_ERR)
c       
c             Determine the minimum Gaussian exponent.
c       

              zmin = 1.D+06
              do iprimo = 1,nprim
                zmin = min(zprim(iprimo),zmin)
              enddo
              do icontr=1,ncontr
                 l=l0(icontr)
c       
c                Only include those basis functions that are "non-zero"
c                for at least one point in the sampling set.
c       
                 nshbf=(((l+1)*(l+2))/2)
                 if(isphere.eq.1) then
                    nshbf=(2*l+1)
                 endif
c
c                pre-screening for worst case (max radius)
c
                 if (r_q0.lt.gaussian_range(l0(ncontr),acc_loc,zmin))
     &           then
                       do m=1,nshbf
                          ibf(n1+m ) = npt+m
                       enddo
                       docset(icset) = .true.
                       if (iniz(ictr).eq.0) iniz(ictr)=n1+1
                       n1=n1+nshbf
                 endif
                 npt=npt+nshbf
              enddo
           enddo
c           
           ifin(ictr)= n1
        endif
c
  400 continue
c
      grid_nbfm = n1
c
      return
      end
c
c     Precompute relevant basis info for XC
c
c     BGJ - 9/00
c
C> \brief Precompute basis info for XC
C>
C> Precompute information on the particular basis set in hand.
C> The tables generated are used in the basis function evaluation.
C>
      Subroutine xc_make_basis_info(basis_hand, bas_cent_info,
     &     bas_cset_info, mcenters)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
      integer basis_hand        !< [Input] The basis set handle
      integer mcenters          !< [Input] The number of centers with basis
                                !< functions
      integer bas_cent_info(3,mcenters) !< [Output] Basis set info for each
                                        !< center
                                        !< - 1. The maximum angular momentum
                                        !< - 2. The lowest rank shell on the
                                        !<   atom
                                        !< - 3. The highest rank shell on the
                                        !<   atom
      integer bas_cset_info(6,*) !< [Output] Basis set info for each shell
                                 !< - 1. The lowest rank basis function in the
                                 !<   shell
                                 !< - 2. The highest rank basis function in the
                                 !<   shell
                                 !< - 3. The angular momentum of the shell
                                 !< - 4. The number of primitives in the shell
                                 !< - 5. The number of contractions in the shell
                                 !< - 6. Whether the shell is Cartesian or
                                 !<   spherical harmonic
c
#include "bas.fh"
c
      integer ictr, icset, ifirst, ilast, nprim, ncontr, l, isphere
      integer ic1, ic2
c
      do ictr = 1,mcenters
         bas_cent_info(1,ictr) = 0 ! max angular momentum
         if (.not.bas_ce2cnr(basis_hand,ictr,ic1,ic2))
     &        call errquit('Exiting in xc_make_basis_info',1, BASIS_ERR)
         bas_cent_info(2,ictr) = ic1
         bas_cent_info(3,ictr) = ic2
c
         do icset = ic1,ic2
c         
c     Get info about current contraction set - first and last basis function,
c     angular momentum, number of primitives, number of contractions and
c     whether spherical harmomic
c       
            if (.not. bas_cn2bfr(basis_hand, icset, ifirst, ilast))
     &           call errquit('Exiting in xc_make_basis_info',1,
     &       BASIS_ERR)
            if (.not. bas_continfo(basis_hand, icset, l, nprim,
     &           ncontr, isphere))
     &           call errquit('Exiting in xc_make_basis_info',2,
     &       BASIS_ERR)
c
            if (l.lt.0) then
#if 0
               call errquit('L code < 0 in xc_make_basis_info',1,
     &       BASIS_ERR)
#else
               l=-1
#endif
            endif
c
            bas_cent_info(1,ictr) = max(bas_cent_info(1,ictr),l)
c
            bas_cset_info(1,icset) = ifirst
            bas_cset_info(2,icset) = ilast
            bas_cset_info(3,icset) = l
            bas_cset_info(4,icset) = nprim
            bas_cset_info(5,icset) = ncontr
            bas_cset_info(6,icset) = isphere
         enddo
      enddo
c
      return
      end
