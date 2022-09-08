c
c     Generalized to do XC hessian as well as gradient; made original
c     routine a wrapper for the more general routine
c
c     BGJ - 8/98
c
      Subroutine dftg_gridv0(rtdb, geom,basis,ipol,nbf_ao,
     ,     iga_dens, force, natoms, oactive,
     &     nactive, rdens_atom, cetobfr, cd_coef, adft)
c
c     $Id$
c
      implicit none      
c
      integer rtdb,geom,ipol,nbf_ao
      integer natoms
      integer nactive
      logical oactive(natoms)
      integer iga_dens(2)
      double precision force(3*natoms)
      integer cetobfr(*)
      double precision rdens_atom(*)
      double precision cd_coef(*)
      logical adft
      integer basis
c
c Daniel (2-27-13): The last 0 in the call is the number of perturbed
c densities.  Since this routine is a wrapper for doing ground state
c gradients, there are no perturbed densities in this call.  The false
c is for consistency because these routines are also used for TDDFT
c gradients. 
      call dftgh_gridv0(rtdb, geom, basis,ipol,nbf_ao,
     .     iga_dens, force, 0d0, 0, 1, natoms,
     &                  oactive, nactive, rdens_atom, cetobfr,
     1                  0, .false.,cd_coef, adft)
c
      return
      end

c
c     Generalized to do XC gradient and hessian
c
      Subroutine dftgh_gridv0(rtdb, geom, basis,ipol,nbf_ao,
     ,     iga_dens, force, hess, g_rhs,
     &                        calc_type, natoms, oactive, nactive,
     &                        rdens_atom, cetobfr,
     1     ngxc, triplet,cd_coef,adft)
c      
c $Id$
c     
      implicit none      
#include "errquit.fh"
c
      logical adft
      integer basis
      integer rtdb,geom,ipol,nbf_ao ![in]
      integer natoms
      integer nactive
      logical oactive(natoms)
      integer iga_dens(2)
      double precision cd_coef(*)
      double precision Force(3*natoms), hess(3*natoms,3*natoms)
      integer g_rhs(*)          ! [input] handles for CPKS RHS mats
      integer calc_type         ! [input] 1 = XC gradient
                                !         2 = XC hessian
                                !         3 = XC CPKS RHS
                                !         4 = XC kernel gradient
                                !             (XC-third derivatives)
                                !         5 = XC potential gradient
                                !             (XC-second derivatives)
c
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "msgids.fh"
#include "cgridfile.fh"
c
#include "dft2drv.fh"
c Daniel (1-14-13): Header file for XC-third derivatives
#include "dft3drv.fh"
c
c
      integer cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      logical GRAD, GRADCHI, HESSCHI, DER3CHI, kske, kslap
      logical xc_chkgrad,xc_rep_close, xc_chktau, xc_chklap
      external xc_chkgrad,xc_rep_close, xc_chktau, xc_chklap
           
      logical ldew, oprint
      double precision eps
      parameter(eps=1.d-8)
      double precision rho_int
      integer me, nproc, NXprim, nbf_max, nqstep,
     &        npol, nxyz, max_at_bf,  nqtot,
     &        nd_qwght, len_Amat, len_Amat2, len_Cmat, 
     &        len_Cmat2, ncoef_max, max_at_bf2, nelem, 
     &        n1, icount, ilo, ihi, iatom, nd_func
      integer lxyz, ixyz, l_xyzm,k_xyzm,lcharge, icharge,  lAmat, iAmat, 
     &        lCmat, iCmat, lqxyz, iqxyz, lqwght, iqwght, ld_qwght, 
     &        id_qwght, lfunc, ifunc, lrq0, irq0, lrq, irq, lniz, iniz, 
     &        lfin, ifin, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, lBmat, iBmat,
     &        lhh, ihh,  landex, iandex, idelrho, ldelrho,
     &        ltags, itags,
     &        iMmat, lMmat, ittau, lttau,
     &        iLmat, lLmat
      logical do_force, do_hess, do_cpks_r
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2, len_Bmat
      integer l_Mmat2, k_Mmat2
      integer lzzz,  izzz, lff, iff, lppp, ippp
      integer libf_ao, iibf_ao,k_curatoms,l_curatoms
      integer lrchi_atom, irchi_atom,
     &        lrdelchi_atom, irdelchi_atom
      integer  l_drhonuc, k_drhonuc,l_ddelrhonuc, k_ddelrhonuc,
     ,     l_d2scr, k_d2scr
      integer nscr,icrad,iprad,iang
      integer llmax, MXcont, max_elem_ang_scr,nxang
      integer nbf_ao_mxnbf_ce
      integer lgthpmat
      integer len_Mmat,len_ttau,len_Mmat2
      integer len_Lmat
      integer lscr,iscr
      integer nxtask
      external nxtask
c Daniel (1-14-13): This is needed for evaluating the gradient of the
c XC-kernel, which can be directly evaluated from the XC-third
c derivative multiplied by the gradient of the ground state density.
      logical do_gxc
      integer len_Amat3, len_Cmat3
      integer l_Amat3, k_Amat3, l_Cmat3, k_Cmat3
      integer ngxc
      integer len_rho, len_delrho
      integer l_rho_save, k_rho_save, l_delrho_save, k_delrho_save
      integer mbf_dchiffp
c Daniel (2-6-13): Variable for doing restricted triplet gradients for
c the XC-kernel term in TDDFT gradients
      logical triplet
c Daniel (2-7-13): Variable for the alternative approach to the dVxc*P
c term in the TDDFT gradients.
      logical do_tddftvxc
c
      do_force = calc_type .eq. 1
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      do_gxc = calc_type .eq. 4
      do_tddftvxc = calc_type .eq. 5
c      write(6,*)'*** Entered dftgh_quadv0, logicals =',do_force,
c     &     do_hess, do_cpks_r
      if (.not. (do_force .or. do_hess .or. do_cpks_r 
     1          .or. do_gxc .or. do_tddftvxc)) then
c      if (.not. (do_force .or. do_hess .or. do_cpks_r 
c     1          .or. do_gxc)) then
         call errquit('dftg_gridb0: illegal calculation type',0,
     &       INPUT_ERR)
      endif
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('quadrature', print_high)

c
c     Define the number of basis functions in the AO Gaussian basis set.
c
      if (.not. bas_nbf_ce_max(basis,nbf_ao_mxnbf_ce))
     .     call errquit('dftggridv0: cemax failed ',0, BASIS_ERR)
       if (.not. bas_nprim_cn_max(basis, NXprim))
     &   call errquit('exiting in dftg_gridv0',6, BASIS_ERR)
      if (.not.bas_nbf_cn_max(basis, nbf_max))
     &   call errquit('exiting in dftg_gridv0.',7, BASIS_ERR)
c
c     Define the max number of basis functions, contractions,
c     and primitives comprising any given shell of this basis set.
c
      if (.not. rtdb_get(rtdb, 'dft:grad_nq', mt_int, 1, nqstep))
     &   nqstep=1
c
c     see if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad()
c
c     see if any Kinetic Energy Density is not null
c
      kske=xc_chktau()
c
c     see if any Laplacian dependent functional is not null
c
      kslap=xc_chklap()
c
c     weight derivatives switched on by default
c
      if (.not. rtdb_get(rtdb, 'dft:weight derivatives', mt_log, 1, 
     &   ldew))ldew=.true.
c     weight deriv bit missing for metaggas
cedo      if(kske) ldew=.false.
      if(ldew.and.(me.eq.0).and.oprint) 
     .     write(LuOut,*) '  weight derivatives on'
c
c     !!! Are these logicals ever used? Are they in a common block? !!!
      if (kslap) then
         gradchi = .true.
         hesschi = .true.
         der3chi = .true.
      elseif (grad.or.kske) then
         gradchi = .true.
         hesschi = .true.
         der3chi = do_hess
      else
         gradchi = .true.
         hesschi = do_hess
         der3chi = .false.
      endif
c      write(6,*)'*** dftgh_quadv0, hesschi der3chi =',hesschi,der3chi
c
      npol = (ipol*(ipol+1))/2
c
      nxyz = 3*natoms
c
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit('dftg_gridv0: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',l_xyzm,k_xyzm))
     &   call errquit('dftg_gridv0: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'charge',lcharge,icharge))
     &   call errquit('dftg_gridv0: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,natoms*16,'tags',ltags,itags))
     &   call errquit('dftg_gridv0: cannot allocate tags',0, MA_ERR)
c
      if(.not.geom_cart_get(geom, natoms, byte_mb(itags), 
     &          dbl_mb(ixyz), dbl_mb(icharge)))
     &     call errquit('dftg_gridv0: geomcartget failed',0, GEOM_ERR)
c
      if (.not.ma_pop_stack(ltags))
     &   call errquit('dftg_gridv0: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('dftg_gridv0: cannot pop stack',0, MA_ERR)
c     
c     get max dimens of atom block
c     
      max_at_bf = 0
      do iatom = 1, natoms
         if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     &      call errquit('dftg_gridv0: bas_ce2bfr failed', iatom,
     &       BASIS_ERR)
         max_at_bf = max(max_at_bf, ihi-ilo+1)
      enddo
      max_at_bf2 = max_at_bf*max_at_bf
c
c     
c     Divide the multi-center integrand into a set of single center integrands.
c     
      nqtot = 0
      rho_int = 0.D0
c
      len_Cmat  = 3*n_per_rec*ipol
      len_Amat  = n_per_rec*2
      len_Amat2 = n_per_rec*NCOL_AMAT2
      len_Cmat2 = n_per_rec*NCOL_CMAT2
      len_Mmat2 = n_per_rec*NCOL_MMAT2
c Daniel (1-14-13): XC-third derivative matrix lengths.
      len_Amat3 = n_per_rec*NCOL_AMAT3
      len_Cmat3 = n_per_rec*NCOL_CMAT3
c Daniel (1-14-13): Presumably this is here to initialize these
c variables to avoid compiler warnings if they aren't used...
      k_Amat2 = 0
      k_Cmat2 = 0
      k_Mmat2 = 0
      k_drhonuc = 0
      k_ddelrhonuc = 0
      k_d2scr = 0
      k_Amat3 = 0
      k_Cmat3 = 0
c
      if (.not.MA_Push_Get(mt_dbl,len_Amat,'Amat',lAmat,iAmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',8, MA_ERR)
c      if (do_hess .or. do_cpks_r) then
      if (do_hess .or. do_cpks_r .or. do_tddftvxc) then
       if (.not.MA_Push_Get(mt_dbl,len_Amat2,'Amat2',l_Amat2,k_Amat2))
     &  call errquit(' exiting in dftg_gridv0: insuff stack',108,MA_ERR)
        if (grad) then
         if (.not.MA_Push_Get(mt_dbl,len_Cmat2,'Cmat2',l_Cmat2,k_Cmat2))
     &     call errquit(' exiting in dftg_gridv0: insuff stack',109,
     &       MA_ERR)
        endif
        if (kske) then
         if (.not.MA_Push_Get(mt_dbl,len_Mmat2,'Mmat2',l_Mmat2,k_Mmat2))
     &     call errquit(' exiting in dftg_gridv0: insuff stack',109,
     &       MA_ERR)
        endif
c Daniel (1-14-13): We need 3rd derivative matrices here.  Presently we
c don't initiallize the kinetic-energy density terms, although these 
c can be added in the future.
      elseif (do_gxc) then
       if (.not.MA_Push_Get(mt_dbl,len_Amat2,'Amat2',l_Amat2,k_Amat2))
     &  call errquit(' exiting in dftg_gridv0: insuff stack',110,MA_ERR)
       if (.not.MA_Push_Get(mt_dbl,len_Amat3,'Amat3',l_Amat3,k_Amat3))
     &  call errquit(' exiting in dftg_gridv0: insuff stack',111,MA_ERR)
       if (grad) then
        if (.not.MA_Push_Get(mt_dbl,len_Cmat2,'Cmat2',l_Cmat2,k_Cmat2))
     &    call errquit(' exiting in dftg_gridv0: insuff stack',112,
     &      MA_ERR)
        if (.not.MA_Push_Get(mt_dbl,len_Cmat3,'Cmat3',l_Cmat3,k_Cmat3))
     &    call errquit(' exiting in dftg_gridv0: insuff stack',113,
     &      MA_ERR)
       endif        
      endif
c
      if (.not.MA_Push_Get(mt_dbl,len_Cmat,'Cmat',lCmat,iCmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',9, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',10,
     &       MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',11,
     &       MA_ERR)
      nd_qwght = 1
      if (ldew.or.do_gxc)nd_qwght = 3*n_per_rec*natoms
      if (.not.MA_Push_get(mt_dbl,nd_qwght,'qwght',ld_qwght,id_qwght))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',11,
     &       MA_ERR)
      nd_func=1
      if (ldew.or.do_gxc) nd_func=n_per_rec
      if (.not.MA_Push_get(mt_dbl,nd_func,'func',lfunc,ifunc))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',11,
     &       MA_ERR)
      if (.not.MA_Push_get(mt_dbl,natoms,'rq0',lrq0,irq0))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',12,
     &       MA_ERR)
      if (.not.MA_Push_get(mt_dbl,natoms*n_per_rec,'rq',lrq,irq))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
      if (.not.MA_Push_get(mt_int,natoms,'atom list',lniz,iniz))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
      if (.not.MA_Push_get(mt_int,natoms,'atom nbf',lfin,ifin))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,NXprim,'expo',lexpo,iexpo))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
      ncoef_max = nbf_max * NXprim
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
      lgthpmat=max_at_bf2
      if (do_cpks_r)  lgthpmat= lgthpmat*3
      if (.not.MA_Push_Get(mt_dbl,lgthpmat,'DENSF',lFmat,iFmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
c
c this for a single nfock pass in xc_tabcd
c      if(do_cpks_r) lgthpmat=lgthpmat*3*natoms
c
      if (.not.MA_Push_Get(mt_dbl,lgthpmat,'DENS',lPmat,iPmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',17,
     &       MA_ERR)
c Daniel (1-16-13): Need to account for having the ground state
c density matrices and perturbed density matrices for the 
c XC-kernel gradient and XC-potential gradient for TDDFT.
c      if (do_gxc) then
      if (do_gxc .or. do_tddftvxc) then
c        len_rho = npol*n_per_rec*(ngxc+1) + n_per_rec
        len_rho = ipol*n_per_rec*(ngxc+1) + n_per_rec
      else
        len_rho = npol*n_per_rec
      endif
      if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho',lrho,irho))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
c      if (.not.MA_Push_Get(MT_Dbl,npol*n_per_rec,'rho',lrho,irho))
c     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
c     &       MA_ERR)
c
c Daniel (1-16-13): Ditto to the above comment for the density gradient
c      if (do_gxc) then
      if (do_gxc .or. do_tddftvxc) then
c        len_delrho = 3*npol*n_per_rec*(ngxc+1) + 3*n_per_rec
        len_delrho = 3*ipol*n_per_rec*(ngxc+1) + 3*n_per_rec
      else
        len_delrho = 3*npol*n_per_rec
      endif
      if (.not.MA_Push_Get(MT_Dbl,len_delrho,'delrho',ldelrho,idelrho))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
c      nelem = 3*npol*n_per_rec
c      if (.not.MA_Push_Get(MT_Dbl,nelem,'delrho',ldelrho,idelrho))
c     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
c     &       MA_ERR)
c
c Daniel (1-17-13): We need extra matrices for storing the perturbed 
c electron density and its gradient, because we may overwrite old 
c matrices with new ones.  For all other calculation types, we don't
c need these matrices. 
      k_rho_save = 0
      k_delrho_save = 0
c      if (do_gxc) then
      if (do_gxc .or. do_tddftvxc) then
c        len_rho = npol*n_per_rec*(ngxc+1) + n_per_rec
        len_rho = ipol*n_per_rec*(ngxc+1) + n_per_rec
        if (.not.ma_push_get(mt_dbl, len_rho, 'rho_save', l_rho_save,
     1    k_rho_save))
     2    call errquit(
     3      'dftgh_gridv0: Could not allocate rho_save', 0, MA_ERR)
      endif
c
c      if (do_gxc) then
      if (do_gxc .or. do_tddftvxc) then
c        len_delrho = 3*npol*n_per_rec*(ngxc+1) + 3*n_per_rec
        len_delrho = 3*ipol*n_per_rec*(ngxc+1) + 3*n_per_rec
        if (.not.ma_push_get(mt_dbl, len_delrho, 'delrho_save',
     1    l_delrho_save, k_delrho_save))
     2    call errquit(
     3      'dftgh_gridv0: Could not allocate delrho_save', 0, MA_ERR)
      endif
c
      len_ttau = 1
      if (kske) len_ttau = ipol*n_per_rec
      if (.not.MA_Push_Get(MT_Dbl,len_ttau,'ttau',lttau,ittau))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13,
     &       MA_ERR)
c
      len_Mmat = 1
      if (kske) len_Mmat = ipol*n_per_rec
      if (.not.MA_Push_Get(MT_Dbl,len_Mmat,'Mmat',lMmat,iMmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',
     &   13, MA_ERR)
c
      len_Lmat = 1
      if (kslap) len_Lmat = ipol*n_per_rec
      if (.not.MA_Push_Get(MT_Dbl,len_Lmat,'Lmat',lLmat,iLmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',
     &   13, MA_ERR)
c
c     Note: hessian and cpks rhs need a larger version of Bmat
c     See xc_quadv0.F
c
      len_Bmat = n_per_rec*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(mt_dbl,len_Bmat,'Bmat',lBmat,iBmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',31,
     &       MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,3*n_per_rec*max_at_bf,
     H     'hess scratch',lhh,ihh))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',32,
     &       MA_ERR)
c
c     stuff for neighbours list
c
      if (.not.MA_Push_Get(MT_int,natoms,'iandex',landex,iandex))
     &   call errquit(' exiting in quadvxc0: insuff stack',77,
     &       MA_ERR)
c
c
c     allocate scratch array for grid compression
c
      if (.not. MA_Push_Get(MT_Int, nbf_ao, 'ibf_ao',
     &     libf_ao, iibf_ao))call errquit
     &     (' exiting in dftg_gridv0: insuff stack',9, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit('dft_scf: cannot allocate rchi_atom',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'rdelchi_atom',
     &     lrdelchi_atom,irdelchi_atom))call errquit
     &     ('dftg_gridv0: cannot allocate rdelchi_atom',0, MA_ERR)
      if (.not. MA_Push_Get(mt_dbl, nbf_ao*max_at_bf, 'PPPP', lppp, 
     &     ippp))call errquit
     &     (' dftg_gridv0: not enough stack',13, MA_ERR)
c         
c !!! BGJ - probably don't need to allocate some of this if .not.do_force
      if (.not. MA_Push_Get(mt_dbl, nbf_ao*max_at_bf, 'F', lff, iff))
     &     call errquit(' dftg_gridv0: not enough stack',13, MA_ERR)
      if (.not. MA_Push_Get(mt_dbl, 3*nbf_ao*max_at_bf, 
     Z     'zzz', lzzz, izzz))call errquit
     &     (' dftg_gridv0: not enough stack',13, MA_ERR)
      if (.not.MA_push_get(MT_INT, natoms, 'curatoms',
     &     l_curatoms, k_curatoms)) call
     &     errquit('dftg_gridv0: cant allocate curatoms',1, MA_ERR)
      if (do_hess .or. do_cpks_r) then
c     
c     Set up space for derivatives of density wrt significant nuclei
c
         if (.not.MA_push_get(MT_DBL, ipol*n_per_rec*3*natoms,
     &        'drhonuc', l_drhonuc, k_drhonuc)) call
     &        errquit('dftg_gridv0: cannot allocate drhonuc',1, MA_ERR)
         if (grad) then
            if (.not.MA_push_get(MT_DBL,
     &           ipol*3*n_per_rec*3*natoms, 'ddelrhonuc',
     &           l_ddelrhonuc, k_ddelrhonuc)) call
     &           errquit('dftg_gridv0: cant alloc ddelrhonuc',
     &           1, MA_ERR)
         endif
         if(do_hess) then
            if (.not.MA_push_get(mt_dbl, 
     1           15*n_per_rec*nbf_ao_mxnbf_ce, 'd2scr',
     &           l_d2scr, k_d2scr)) call errquit
     &           ('dftg_gridv0: not enough stack for d2scr',0, MA_ERR)
         endif
      endif
c    
c     compute scratch space for int_eval_gbsets
c    
c
c     Define the number of basis functions in the AO Gaussian basis set.
c
      if (.not. bas_ncontr_cn_max(basis, MXcont))
     &   call errquit('Exiting in dftg_quadv0.',1, BASIS_ERR)
      if (.not. bas_high_angular(basis, NXang))
     &   call errquit('exiting in dftg_quadv0.',7, BASIS_ERR)
      llMAX=((NXang+2)*(NXang+1))/2
c scratch for eval_func
      if (NXang .ge. 3)then
         max_elem_ang_scr = 28
      elseif (NXang .ge. 1)then
         max_elem_ang_scr = 9
      else 
         max_elem_ang_scr = 3
      endif 
      max_elem_ang_scr=max_elem_ang_scr*llMAX

      if (der3chi) then
c     
c     Allocate scratch space.
c     mxcont is used here instead of NXprim as in new_eval_gbsets
c     this is correct, but sometimes fails, hence the choice
c     in new_eval_gsets.  here, the length is also doubled...
c     
         icrad = 1              ! contracted radial parts
         iprad = icrad + 4*n_per_rec*MXcont ! primitive radial parts
         iang  = iprad          ! also angular parts
c     
c     note 4*NXprim for longer primitive and contracted radial
c     scratch storage
c     
         nscr=iang+max(4*n_per_rec*NXprim,max_elem_ang_scr*n_per_rec)
c     
      else if (hesschi)then
c     
c     Allocate scratch space.
c     mxcont is used here instead of NXprim as in new_eval_gbsets
c     this is correct, but sometimes fails, hence the choice
c     in new_eval_gsets.  here, the length is also doubled...
c     
         icrad = 1              ! contracted radial parts
         iprad = icrad + 3*n_per_rec*MXcont ! primitive radial parts
         iang  = iprad          ! also angular parts
c     
c     note 3*NXprim for longer primitive and contracted radial
c     scratch storage
c     
         nscr=iang+max(3*n_per_rec*NXprim,max_elem_ang_scr*n_per_rec)
c     
      elseif (gradchi)then
c     
c     Allocate scratch space.
c     
         icrad = 1
         iprad = icrad + 2*n_per_rec*NXprim
         iang  = iprad
         nscr  = iang+max(2*n_per_rec*NXprim,9*n_per_rec*llMAX) - 1
         nscr = 2*nscr
      else 
c
c                 Allocate scratch space.
         
         icrad = 1
         iprad = icrad + n_per_rec*NXprim
         iang  = iprad
         nscr  = iang  + max(n_per_rec*NXprim,6*n_per_rec*llMAX) 
     &        - 1
         nscr = 2*nscr
c     
      endif 
c     
      if (.not.adft) then
        call xc_rep_init(rtdb, geom, basis,iga_dens,iga_dens,
     &     nbf_ao,ipol,do_force,.true.)
      endif
      
c
      if(grid_written) then
c Daniel (1-14-13): We need to pass 3rd derivatives information 
c through here so we can evaluate the gradient of the XC-kernel.
      call dftg_gridv0a(rtdb,iga_dens,ldew,
     &     max_at_bf,
     &     rho_int, nqtot, int_mb(iniz), int_mb(ifin),
     &     dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ifunc),
     &     dbl_mb(id_qwght), dbl_mb(iqwght),
     &     dbl_mb(iqxyz), dbl_mb(ixyz),dbl_mb(k_xyzm), 
     &     int_mb(iandex), dbl_mb(iexpo),
     &     dbl_mb(iAmat), dbl_mb(k_Amat2),
     &     dbl_mb(iCmat), dbl_mb(k_Cmat2),
     &     dbl_mb(iBmat), 
     &     dbl_mb(iFmat), dbl_mb(iPmat),
     &     dbl_mb(irq), dbl_mb(irq0), dbl_mb(ihh),
     ,     int_mb(iibf_ao),int_mb(k_curatoms),
     ,     dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),dbl_mb(ippp),
     ,     dbl_mb(iff),dbl_mb(izzz),
     .     dbl_mb(k_drhonuc),dbl_mb(k_ddelrhonuc),dbl_mb(k_d2scr),
     .     nscr,
     &     force, hess, g_rhs, do_force, do_hess,
     &     do_cpks_r, natoms, oactive, nactive,
     &     rdens_atom,  cetobfr,
     &     dbl_mb(ittau), 
     &     dbl_mb(iMmat), dbl_mb(k_Mmat2),
     &     dbl_mb(iLmat),
     1     do_gxc, dbl_mb(k_Amat3), dbl_mb(k_Cmat3), ngxc,
     2     dbl_mb(k_rho_save), dbl_mb(k_delrho_save), triplet, 
     3     do_tddftvxc,cd_coef)
      else
       if (.not.MA_Push_Get(mt_dbl, nscr, 'scrCD', lscr, iscr))
     &     call errquit
     &     (' dftg_gridv0: insuff stack',12, MA_ERR)
c Daniel (1-14-13): We need to pass a logical variable relevant to
c 3rd derivatives through here to prevent them from happening, since
c some compilers (such as the newer iforts) default logical variables
c to true.
      call dftg_griddo(rtdb,
     &     iga_dens,ldew,
     &     ncoef_max, max_at_bf,
     &     rho_int, nqtot, int_mb(iniz), int_mb(ifin),
     &     dbl_mb(irho), dbl_mb(idelrho), 
     &     dbl_mb(ifunc), dbl_mb(id_qwght), dbl_mb(iqwght),
     &     dbl_mb(iqxyz), dbl_mb(ixyz),
     &     int_mb(iandex), dbl_mb(iexpo), dbl_mb(icoef),
     &     dbl_mb(iAmat), dbl_mb(k_Amat2),dbl_mb(iCmat), 
     &     dbl_mb(k_Cmat2), dbl_mb(iBmat), 
     &     dbl_mb(iFmat), dbl_mb(iPmat),dbl_mb(irq), 
     &     dbl_mb(irq0), dbl_mb(ihh),
     ,     int_mb(iibf_ao),int_mb(k_curatoms),dbl_mb(irchi_atom),
     ,     dbl_mb(irdelchi_atom),dbl_mb(ippp),
     ,     dbl_mb(iff),dbl_mb(izzz),
     .     dbl_mb(k_drhonuc),dbl_mb(k_ddelrhonuc),dbl_mb(k_d2scr),
     .     dbl_mb(iscr),nscr,
     &     force, hess, g_rhs, do_force, do_hess,
     &     do_cpks_r, natoms, oactive, nactive,
     &     rdens_atom,  cetobfr,
     &     dbl_mb(ittau), 
     &     dbl_mb(iMmat), dbl_mb(k_Mmat2),
     &     dbl_mb(iLmat),
     1     do_gxc, dbl_mb(k_Amat3), dbl_mb(k_Cmat3), 0,
     2     dbl_mb(k_rho_save), dbl_mb(k_delrho_save),
     3     do_tddftvxc,cd_coef)
      if (.not.ma_pop_stack(lscr))
     &   call errquit('dftg_gridv0: cannot pop stack',21, MA_ERR)
      endif
      if (.not.adft) then
        if(.not.xc_rep_close(rtdb, nbf_ao,ipol,ipol,
     D     iga_dens, .true.)) call
     .     errquit(' dftggridv0: xcrepclose failed ',0, 0)
      endif
      if (.not.ma_chop_stack(lhh))
     &   call errquit('dftg_gridv0: cannot chop stack',1, MA_ERR)
      if (.not.ma_chop_stack(lCmat))
     &   call errquit('dftg_gridv0: cannot pop stack',0, MA_ERR)
c Daniel (1-14-13): Third derivative stuff
      if (do_gxc) then
         if (grad) then
            if (.not.ma_pop_stack(l_Cmat3))
     &         call errquit('dftg_gridv0: cannot pop stack Cmat3',0, 
     &            MA_ERR)
            if (.not.ma_pop_stack(l_Cmat2))
     &         call errquit('dftg_gridv0: cannot pop stack Cmat2',0, 
     &            MA_ERR)
         endif
         if (.not.ma_pop_stack(l_Amat3))
     &      call errquit('dftg_gridv0: cannot pop stack Amat3',0, 
     &         MA_ERR)
         if (.not.ma_pop_stack(l_Amat2))
     &      call errquit('dftg_gridv0: cannot pop stack Amat2',0, 
     &         MA_ERR)
c      elseif (do_hess .or. do_cpks_r) then
      elseif (do_hess .or. do_cpks_r .or. do_tddftvxc) then
         if (grad) then
            if (.not.ma_pop_stack(l_Cmat2))
     &           call errquit('dftg_gridv0: cannot pop stack',0, MA_ERR)
         endif
         if (.not.ma_pop_stack(l_Amat2))
     &      call errquit('dftg_gridv0: cannot pop stack',0, MA_ERR)
      endif
      if (.not.ma_pop_stack(lAmat))
     &   call errquit('dftg_gridv0: cannot pop stack',0, MA_ERR)
c
      icount = 1
      n1 = nxtask(-nproc,icount)
c
      call ga_igop(Msg_nqtot,nqtot,1,'+')
      call ga_dgop(Msg_Excrho,rho_int,1,'+')
c     
      if (me.eq.0 .and. oprint)then
         write(LuOut,'(a)')'         Num. Quad. of the Electron Density'
         write(LuOut,'(a,f20.10)')' integrated electron density:  ',
     &                            rho_int
         write(LuOut,'(a,i20)')' no. of quad. points/center :  ',
     &                         nqtot/natoms
         write(LuOut,*)
      end if
c
      if (.not.ma_free_heap(lxyz))
     &   call errquit('dftg_gridv0: cannot free heap',0, MA_ERR)
      if (.not.ma_free_heap(l_xyzm))
     &   call errquit('dftg_gridv0: cannot free heap',1, MA_ERR)
c
      return 
c     
      end
      subroutine oactive_to_curatoms(oactive, natoms, curatoms,
     &                               ncuratoms)
c
c     Converts list of active logicals to indexing array, mapping
c     an atom to its position in a compressed list of atoms, or 0 if
c     atom is not active.
c
c     !!! Make more sophisticated to use cutoffs !!!
c
      implicit none
c
      integer natoms           ! Total number of atoms [input]
      logical oactive(natoms)  ! Whether atom is active [input]
      integer curatoms(natoms) ! Mapping array [output]
      integer ncuratoms        ! Number of significant atoms [output]
c
      integer i
c
      ncuratoms = 0
      do i = 1, natoms
         curatoms(i) = 0
         if (oactive(i)) then
            ncuratoms = ncuratoms + 1
            curatoms(i) = ncuratoms
         endif
c         write(*,*)'oactive_to_curatoms:',i,oactive(i),curatoms(i)
      enddo
c      write(*,*)'oactive_to_curatoms: ncuratoms',ncuratoms
c
      return
      end
