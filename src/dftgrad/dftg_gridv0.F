c
c     Generalized to do XC hessian as well as gradient; made original
c     routine a wrapper for the more general routine
c
c     BGJ - 8/98
c
      Subroutine dftg_gridv0(rtdb, geom,ao_bas_han,ipol,nbf_ao,
     ,     iga_dens, force, natoms, oactive,
     &                       nactive, rdens_atom, cetobfr)
c
c     $Id: dftg_gridv0.F,v 1.7 2000-08-25 20:05:11 edo Exp $
c
      implicit none      
c
      integer rtdb,geom,ao_bas_han,ipol,nbf_ao
      integer natoms
      integer nactive
      logical oactive(natoms)
      integer iga_dens(2)
      double precision force(3*natoms)
      integer cetobfr(*)
      double precision rdens_atom(*)
c
      call dftgh_gridv0(rtdb, geom, ao_bas_han,ipol,nbf_ao,
     .     iga_dens, force, 0d0, 0, 1, natoms,
     &                  oactive, nactive, rdens_atom, cetobfr)
c
      return
      end

c
c     Generalized to do XC gradient and hessian
c
      Subroutine dftgh_gridv0(rtdb, geom, ao_bas_han,ipol,nbf_ao,
     ,     iga_dens, force, hess, g_rhs,
     &                        calc_type, natoms, oactive, nactive,
     &                        rdens_atom, cetobfr)
c      
c $Id: dftg_gridv0.F,v 1.7 2000-08-25 20:05:11 edo Exp $
c     
      implicit none      
c
      integer rtdb,geom,ao_bas_han,ipol,nbf_ao ![in]
      integer natoms
      integer nactive
      logical oactive(natoms)
      integer iga_dens(2)
      double precision Force(3*natoms), hess(3*natoms,3*natoms)
      integer g_rhs(*)          ! [input] handles for CPKS RHS mats
      integer calc_type         ! [input] 1 = XC gradient
                                !         2 = XC hessian
                                !         3 = XC CPKS RHS
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
c
c
      integer cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      logical LResult, GRAD, GRADCHI, HESSCHI, DER3CHI
      logical xc_chkgrad
      external xc_chkgrad
      
      logical ldew, oprint
      double precision eps
      parameter(eps=1.d-8)
      double precision rho_int
      integer me, nproc, NXprim, nbf_max, nqstep,
     &        npol, nxyz, max_at_bf,  nqtot,
     &        nd_qwght, mA, mA2, mC, mC2, ncoef_max, max_at_bf2, nelem, 
     &        n1, icount, ilo, ihi, iatom, nd_func
      integer lxyz, ixyz, l_xyzm,k_xyzm,lcharge, icharge,  lAmat, iAmat, 
     &        lCmat, iCmat, lqxyz, iqxyz, lqwght, iqwght, ld_qwght, 
     &        id_qwght, lfunc, ifunc, lrq0, irq0, lrq, irq, lniz, iniz, 
     &        lfin, ifin, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, lBmat, iBmat, lDmat, iDmat,
     &        lhh, ihh,  landex, iandex, idelrho, ldelrho,
     &        ltags, itags
      logical do_force, do_hess, do_cpks_r
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2, len_Bmat, len_Dmat
      integer lzzz,  izzz, lff, iff, lppp, ippp,k_bp,l_bp
      integer libf_ao, iibf_ao,k_curatoms,l_curatoms
      integer lrchi_atom, irchi_atom,
     &        lrdelchi_atom, irdelchi_atom
      integer  l_drhonuc, k_drhonuc,l_ddelrhonuc, k_ddelrhonuc,
     ,     l_d2scr, k_d2scr,lscr,iscr,lrscr,irscr
      integer nscr,icrad,iprad,iang
      integer llmax, MXcont, max_elem_ang_scr,nxang
      integer nbf_ao_mxnbf_ce
      integer nxtask
      external nxtask
c     
      do_force = calc_type .eq. 1
      do_hess = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
c      write(6,*)'*** Entered dftgh_quadv0, logicals =',do_force,
c     &     do_hess, do_cpks_r
      if (.not. (do_force .or. do_hess .or. do_cpks_r)) then
         call errquit('dftg_gridb0: illegal calculation type',0)
      endif
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('quadrature', print_high)
c
c     weight derivatives switched on by default
c
      if (.not. rtdb_get(rtdb, 'dft:weight derivatives', mt_log, 1, 
     &   ldew))ldew=.true.
      if(ldew.and.(me.eq.0).and.oprint) 
     .     write(LuOut,*) '  weight derivatives on'
c
c     Define the number of basis functions in the AO Gaussian basis set.
c
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     .     call errquit('dftggridv0: cemax failed ',0)
       if (.not. bas_nprim_cn_max(AO_bas_han, NXprim))
     &   call errquit('exiting in dftg_gridv0',6)
      if (.not.bas_nbf_cn_max(AO_bas_han, nbf_max))
     &   call errquit('exiting in dftg_gridv0.',7)
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
c     !!! Are these logicals ever used? Are they in a common block? !!!
      if (grad) then
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
     &   call errquit('dftg_gridv0: cannot allocate xyz',0)
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',l_xyzm,k_xyzm))
     &   call errquit('dftg_gridv0: cannot allocate xyz',0)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'charge',lcharge,icharge))
     &   call errquit('dftg_gridv0: cannot allocate charge',0)
      if (.not.MA_Push_Get(MT_Byte,natoms*16,'tags',ltags,itags))
     &   call errquit('dftg_gridv0: cannot allocate tags',0)
c
      LResult = geom_cart_get(geom, natoms, byte_mb(itags), 
     &          dbl_mb(ixyz), dbl_mb(icharge))
c
      if (.not.ma_pop_stack(ltags))
     &   call errquit('dftg_gridv0: cannot pop stack',0)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('dftg_gridv0: cannot pop stack',0)
c     
c     get max dimens of atom block
c     
      max_at_bf = 0
      do iatom = 1, natoms
         if (.not. bas_ce2bfr(AO_bas_han, iatom, ilo, ihi))
     &      call errquit('dftg_gridv0: bas_ce2bfr failed', iatom)
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
      mC=3*n_per_rec*ipol
      mA = ipol*n_per_rec
      mA2 = n_per_rec*NCOL_AMAT2
      mC2 = n_per_rec*NCOL_CMAT2
      if (.not.MA_Push_Get(mt_dbl,mA,'Amat',lAmat,iAmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',8)
      if (do_hess .or. do_cpks_r) then
         if (.not.MA_Push_Get(mt_dbl,mA2,'Amat2',l_Amat2,k_Amat2))
     &      call errquit(' exiting in dftg_gridv0: insuff stack',108)
         if (grad) then
            if (.not.MA_Push_Get(mt_dbl,mC2,'Cmat2',l_Cmat2,k_Cmat2))
     &         call errquit(' exiting in dftg_gridv0: insuff stack',109)
         endif
      endif
      if (.not.MA_Push_Get(mt_dbl,mC,'Cmat',lCmat,iCmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',9)
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',10)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',11)
      nd_qwght = 1
      if (ldew)nd_qwght = 3*n_per_rec*natoms
      if (.not.MA_Push_get(mt_dbl,nd_qwght,'qwght',ld_qwght,id_qwght))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',11)
      nd_func=1
      if (ldew) nd_func=n_per_rec
      if (.not.MA_Push_get(mt_dbl,nd_func,'func',lfunc,ifunc))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',11)
      if (.not.MA_Push_get(mt_dbl,natoms,'rq0',lrq0,irq0))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',12)
      if (.not.MA_Push_get(mt_dbl,natoms*n_per_rec,'rq',lrq,irq))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
      if (.not.MA_Push_get(mt_int,natoms,'atom list',lniz,iniz))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
      if (.not.MA_Push_get(mt_int,natoms,'atom nbf',lfin,ifin))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
      if (.not.MA_Push_Get(MT_Dbl,NXprim,'expo',lexpo,iexpo))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
      ncoef_max = nbf_max * NXprim
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
      if (.not.MA_Push_Get(mt_dbl,max_at_bf2,'DENSF',lFmat,iFmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
      if (.not.MA_Push_Get(mt_dbl,max_at_bf2,'DENS',lPmat,iPmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
      if (.not.MA_Push_Get(MT_Dbl,npol*n_per_rec,'rho',lrho,irho))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
c
      nelem = 3*npol*n_per_rec
c 
      if (.not.MA_Push_Get(MT_Dbl,nelem,'delrho',ldelrho,idelrho))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',13)
c
c     Note: hessian and cpks rhs need a larger version of Bmat
c     See xc_quadv0.F
c
      if (do_force) then
         len_Bmat = n_per_rec
      else if (do_hess .or. do_cpks_r) then
         len_Bmat = n_per_rec*nbf_ao_mxnbf_ce
      endif
      if (.not.MA_Push_Get(mt_dbl,len_Bmat,'Bmat',lBmat,iBmat))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',31)
c
c     hessian and cpks rhs also need a Dmat
c
      if (do_hess .or. do_cpks_r) then
         len_Dmat = n_per_rec*3*nbf_ao_mxnbf_ce
         if (.not.MA_Push_Get(mt_dbl,len_Dmat,'Dmat',lDmat,iDmat))
     &        call errquit(' exiting in dftg_gridv0: insuff stack',131)
      else
         len_Dmat = 0
      endif
c
c     !!! Appears that this array is never used !!! ask Jeff !!!
      if (.not.MA_Push_Get(mt_dbl,3*n_per_rec,'hess scratch',lhh,ihh))
     &   call errquit(' exiting in dftg_gridv0: insuff stack',32)
c
c     stuff for neighbours list
c
      if (.not.MA_Push_Get(MT_int,natoms,'iandex',landex,iandex))
     &   call errquit(' exiting in quadvxc0: insuff stack',77)
c
c
c     allocate scratch array for grid compression
c
      if (.not. MA_Push_Get(MT_Int, nbf_ao, 'ibf_ao',
     &     libf_ao, iibf_ao))call errquit
     &     (' exiting in dftg_gridv0: insuff stack',9)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit('dft_scf: cannot allocate rchi_atom',0)
      if (.not.MA_Push_Get(MT_Dbl,natoms,'rdelchi_atom',
     &     lrdelchi_atom,irdelchi_atom))call errquit
     &     ('dftg_gridv0: cannot allocate rdelchi_atom',0)
      if (.not. MA_Push_Get(mt_dbl, nbf_ao, 'PPPP', lppp, 
     &     ippp))call errquit
     &     (' dftg_gridv0: not enough stack',13)
c         
c !!! BGJ - probably don't need to allocate some of this if .not.do_force
      if (.not. MA_Push_Get(mt_dbl, nbf_ao, 'F', lff, iff))
     &     call errquit(' dftg_gridv0: not enough stack',13)
      if (.not. MA_Push_Get(mt_dbl, 3*nbf_ao, 'zzz', lzzz, 
     &     izzz))call errquit
     &     (' dftg_gridv0: not enough stack',13)
      if (.not.MA_push_get(MT_INT, natoms, 'curatoms',
     &     l_curatoms, k_curatoms)) call
     &     errquit('dftg_gridv0: cant allocate curatoms',1)
      if (.not.MA_push_get(MT_INT, natoms*n_per_rec, 'bp',
     &     l_bp, k_bp)) call
     &     errquit('dftg_gridv0: cant allocate bp',1)
      if (do_hess .or. do_cpks_r) then
c     
c     Set up space for derivatives of density wrt significant nuclei
c
         if (.not.MA_push_get(MT_DBL, ipol*n_per_rec*3*nactive,
     &        'drhonuc', l_drhonuc, k_drhonuc)) call
     &        errquit('dftg_gridv0a: cannot allocate drhonuc',1)
         if (grad) then
            if (.not.MA_push_get(MT_DBL,
     &           ipol*3*n_per_rec*3*nactive, 'ddelrhonuc',
     &           l_ddelrhonuc, k_ddelrhonuc)) call
     &           errquit('dftg_gridv0: cant alloc ddelrhonuc',
     &           1)
         endif
         if(do_hess) then
            if (.not.MA_push_get(mt_dbl, 15*n_per_rec, 'd2scr',
     &           l_d2scr, k_d2scr)) call errquit
     &           ('dftg_gridv0: not enough stack for d2scr',0)
         endif
      endif
c    
c     compute scratch space for int_eval_gbsets
c    
c
c     Define the number of basis functions in the AO Gaussian basis set.
c
      if (.not. bas_ncontr_cn_max(ao_bas_han, MXcont))
     &   call errquit('Exiting in dftg_quadv0.',1)
      if (.not. bas_high_angular(AO_bas_han, NXang))
     &   call errquit('exiting in dftg_quadv0.',7)
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
      if (.not.MA_Push_Get(mt_dbl, nscr, 'scrCD', lscr, iscr))
     &     call errquit
     &     (' exiting in dftg_gridv0: insuff stack',12)
      if (.not.MA_Push_Get(mt_dbl,n_per_rec,'rho scr',
     &     lrscr,irscr))call errquit
     &     ('dftg_gridv0: cannot allocate rho scr',13)
      
c
      call dftg_gridv0a(iga_dens,ldew,
     &                   ncoef_max, NXprim, max_at_bf,
     &                   rho_int, nqtot, int_mb(iniz), int_mb(ifin),
     &                   dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ifunc),
     &                   dbl_mb(id_qwght), dbl_mb(iqwght),
     &                   dbl_mb(iqxyz), dbl_mb(ixyz),dbl_mb(k_xyzm), 
     &                   int_mb(iandex), dbl_mb(iexpo), dbl_mb(icoef),
     &                   dbl_mb(iAmat), dbl_mb(k_Amat2),
     &                   dbl_mb(iCmat), dbl_mb(k_Cmat2),
     &                   dbl_mb(iBmat), dbl_mb(iDmat),
     &                   dbl_mb(iFmat), dbl_mb(iPmat), 
     &                   dbl_mb(irq), dbl_mb(irq0), dbl_mb(ihh),
     ,     int_mb(iibf_ao),int_mb(k_curatoms),
     ,     dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),dbl_mb(ippp),
     ,     dbl_mb(iff),dbl_mb(izzz),dbl_mb(k_bp),
     .     dbl_mb(k_drhonuc),dbl_mb(k_ddelrhonuc),dbl_mb(k_d2scr),
     .     dbl_mb(iscr),nscr,dbl_mb(irscr),
     &                   force, hess, g_rhs, do_force, do_hess,
     &                   do_cpks_r, natoms, oactive, nactive,
     &                   rdens_atom,  cetobfr)
      if (.not.ma_chop_stack(lhh))
     &   call errquit('dftg_gridv0: cannot chop stack',1)
      if (do_hess .or. do_cpks_r) then
         if (.not.ma_pop_stack(lDmat))
     &      call errquit('dftg_gridv0: cannot pop stack',0)
      endif
      if (.not.ma_chop_stack(lCmat))
     &   call errquit('dftg_gridv0: cannot pop stack',0)
      if (do_hess .or. do_cpks_r) then
         if (grad) then
            if (.not.ma_pop_stack(l_Cmat2))
     &           call errquit('dftg_gridv0: cannot pop stack',0)
         endif
         if (.not.ma_pop_stack(l_Amat2))
     &      call errquit('dftg_gridv0: cannot pop stack',0)
      endif
      if (.not.ma_pop_stack(lAmat))
     &   call errquit('dftg_gridv0: cannot pop stack',0)
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
     &   call errquit('dftg_gridv0: cannot free heap',0)
      if (.not.ma_free_heap(l_xyzm))
     &   call errquit('dftg_gridv0: cannot free heap',1)
c
      return 
c     
      end
