      Subroutine dftg_gridv0b(nqpts,rad,ictr_buf,iga_dens,
     &                         ncoef_max, NXprim, max_at_bf,
     &                         rho_int, nqtot, iniz, ifin, rho, delrho,
     &                         func, d_qwght, qwght, qxyz, xyz, 
     &                         iandex, expo, ccoef, Amat, Amat2, Cmat,
     &                         Cmat2, Bmat, Dmat, Fmat, Pmat,  rq,
     &                         rq0, hh, 
     .     ibf_ao,curatoms,rchi_atom,rdelchi_atom,ppp,ff,zzz,
     ,     force, hess, g_rhs, do_force,
     ,     llmax,max_elem_ang_scr,nbf_ao_mxnbf_ce2,mxcont,
     ,     acc_ao_gauss, ldew,grad,gradchi,hesschi,der3chi,npol,
     &                         do_hess, do_cpks_r, natoms, oactive,
     &                         nactive, rdens_atom, 
     &                         cetobfr)
c      
c     $Id: dftg_gridv0b.F,v 1.3 2000-04-14 01:47:45 edo Exp $
c     
      implicit none
#include "bas.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "cgridfile.fh"
c
#include "dft2drv.fh"
c !!! BGJ test !!!
#include "bgj.fh"
c !!! BGJ test !!!
      integer nqpts, ictr_buf ! [in]
      double precision rad ! [in]
      integer natoms
      integer ncoef_max, NXprim, max_at_bf
      integer nqtot
      integer nactive
      logical oactive(natoms)
      logical grad, gradchi, hesschi, der3chi
      logical xc_chkgrad
      external xc_chkgrad
      logical ldew
      logical do_force, do_hess, do_cpks_r
c
      integer cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      double precision Force(3*natoms), hess(3*natoms,3*natoms)
      integer g_rhs(*)
      double precision hh(3,n_per_rec)
      double precision rho(n_per_rec,ipol*(ipol+1)/2)
      double precision delrho(n_per_rec,3,ipol)
      double precision Amat(n_per_rec,ipol), Cmat(n_per_rec,3,ipol),
     &                 Bmat(n_per_rec,*), Dmat(n_per_rec,3,max_at_bf)
      double precision Amat2(n_per_rec,NCOL_AMAT2)
      double precision Cmat2(n_per_rec,NCOL_CMAT2)
      double precision Pmat(max_at_bf*max_at_bf),
     &                 Fmat(max_at_bf*max_at_bf)
      double precision func(n_per_rec)
      double precision d_qwght(3,n_per_rec,*)
      double precision qwght(nqpts), qxyz(3,nqpts)
      double precision xyz(3,natoms)
      integer iandex(natoms)
      integer iniz(natoms), ifin(natoms)
      integer ibf_ao(nbf_ao),curatoms(natoms)
      double precision rchi_atom(natoms),rdelchi_atom(3,natoms)
      double precision ppp(*),ff(*),zzz(*)
      double precision ccoef(ncoef_max), expo(NXprim)
      double precision rq(n_per_rec,natoms)
      double precision rq0(natoms)
c
c     external
c
      double precision ddot
      integer nxtask
      External ddot, nxtask
      double precision eps
      parameter(eps=1.d-8)
c
c---> Matrix Elements of the XC Potential & Energy
c
      integer iga_dens(2)
c
c     local crap
c
      integer iang, npol
      integer mbf_ao
      integer llmax, MXcont, max_elem_ang_scr
      integer ibp, lbp
      integer  m4
      integer nscr, iscr, lscr, icrad, iprad, lheschi_ao, iheschi_ao,
     &        ld3chi_ao, id3chi_ao
      integer lrscr, irscr
      integer lchi_ao, ichi_ao, 
     &        lgradchi_ao, idelchi_ao
      integer nbf_ao_mxnbf_ce2
      double precision acc_AO_gauss 
      double precision Ex, Ec
      double precision rho_int ! [in/out]
      double precision dabsmax,wmax
      external dabsmax
      integer grid_nbfm
      external grid_nbfm
c
      integer ncuratoms, l_drhonuc, k_drhonuc,
     &        l_ddelrhonuc, k_ddelrhonuc, l_d2scr, k_d2scr
c     
         
c         Evaluate the AO basis set at each of the quad. points.
c           allocate arrays for exponents and contraction coefficients
c           within int_eval_gbsets mbf is # of orbs with nonzero 
c           points in sampling set.  iniz is a mapping array.
c           currently set up only for cartesian (6d) sets
c         
            mbf_ao = nbf_ao
            call ifill(mbf_ao, 0, ibf_ao, 1)
            mbf_ao=grid_nbfm( AO_bas_han,  natoms, 
     &           ictr_buf,rad,xyz,qxyz,nqpts,
     &           ibf_ao, iniz, ifin, expo, 
     .           minexp,ldiff, acc_AO_gauss)
            if (mbf_ao .eq. 0) return
            if(ldew) then
!               if(lssw) then
                  call dftg_gridssw(d_qwght,qwght, qxyz, xyz, 
     .                 rq, ictr_buf, natoms, iandex,
     .                 nqpts,tol_rho,whichssw)
!               else
!                     call dftg_beckew(d_qwght, qwght, qxyz, xyzm,
!     &                                Rij, rq, dbl_mb(ibp),
!     &                                dbl_mb(izeta), dbl_mb(idzeta),
!     &                                dbl_mb(idp), ictr, qcenters,
!     &                                nq_orig,nq, iandex)
!                  endif
               endif
            
            if (do_hess .or. do_cpks_r) then
c     
c     Set up space for derivatives of density wrt significant nuclei
c
               if (.not.MA_push_get(MT_DBL, ipol*nqpts*3*nactive,
     &              'drhonuc', l_drhonuc, k_drhonuc)) call
     &              errquit('dftg_quadv0a: cannot allocate drhonuc',1)
               if (grad) then
                  if (.not.MA_push_get(MT_DBL,
     &                 ipol*3*nqpts*3*nactive, 'ddelrhonuc',
     &                 l_ddelrhonuc, k_ddelrhonuc)) call
     &                 errquit('dftg_quadv0a: cant alloc ddelrhonuc',
     &                 1)
               endif
            endif
            
c
               if (.not. MA_Push_Get(MT_Dbl, nqpts*mbf_ao, 'chi_ao',
     &            lchi_ao, ichi_ao))call errquit
     &            (' exiting in dftg_quadv0: insuff stack',9)
               if (.not. MA_Push_Get(MT_Dbl, 3*nqpts*mbf_ao, 
     &            'gradchi_ao', lgradchi_ao,idelchi_ao))call errquit
     &            (' exiting in dftg_quadv0: insuff stack',9)
c         
c              hessian is used only for derivatives; set array to length 1
c              here
c         
c
               if (hesschi)then
                  m4 = 6*nqpts*mbf_ao
               else
                  m4 = 1
               endif
c
               if (.not. MA_Push_Get(MT_Dbl, m4, 'heschi_ao', 
     &            lheschi_ao, iheschi_ao))call errquit
     &            (' exiting in dftg_quadv0: insuff stack',9)
c    
               if (der3chi)then
                  m4 = 10*nqpts*mbf_ao
               else
                  m4 = 1
               endif
c
               if (.not. MA_Push_Get(MT_Dbl, m4, 'd3chi_ao', 
     &            ld3chi_ao, id3chi_ao))call errquit
     &            (' exiting in dftg_quadv0: insuff stack',19)
c    
c              compute scratch space for int_eval_gbsets
c    
               if (der3chi) then
c
c                 Allocate scratch space.
c                 mxcont is used here instead of NXprim as in new_eval_gbsets
c                 this is correct, but sometimes fails, hence the choice
c                 in new_eval_gsets.  here, the length is also doubled...
c           
                  icrad = 1 ! contracted radial parts
                  iprad = icrad + 4*nqpts*MXcont ! primitive radial parts
                  iang  = iprad ! also angular parts
c           
c                 note 4*NXprim for longer primitive and contracted radial
c                 scratch storage
c           
                  nscr=iang+max(4*nqpts*NXprim,max_elem_ang_scr*nqpts)
c
               else if (hesschi)then
c
c                 Allocate scratch space.
c                 mxcont is used here instead of NXprim as in new_eval_gbsets
c                 this is correct, but sometimes fails, hence the choice
c                 in new_eval_gsets.  here, the length is also doubled...
c           
                  icrad = 1 ! contracted radial parts
                  iprad = icrad + 3*nqpts*MXcont ! primitive radial parts
                  iang  = iprad ! also angular parts
c           
c                 note 3*NXprim for longer primitive and contracted radial
c                 scratch storage
c           
                  nscr=iang+max(3*nqpts*NXprim,max_elem_ang_scr*nqpts)
c
               elseif (gradchi)then
c
c                 Allocate scratch space.
c
                  icrad = 1
                  iprad = icrad + 2*nqpts*NXprim
                  iang  = iprad
                  nscr  = iang+max(2*nqpts*NXprim,9*nqpts*llMAX) - 1
                  nscr = 2*nscr
               else 
c
c                 Allocate scratch space.

                  icrad = 1
                  iprad = icrad + nqpts*NXprim
                  iang  = iprad
                  nscr  = iang  + max(nqpts*NXprim,6*nqpts*llMAX) 
     &                    - 1
                  nscr = 2*nscr
c
               endif 
c
               if (.not.MA_Push_Get(mt_dbl, nscr, 'scrCD', lscr, iscr))
     &            call errquit
     &            (' exiting in dftg_quadv0: insuff stack',12)
c
               call qdist(rq0, rq, qxyz, xyz, nqpts, natoms)
               if (der3chi) then
                  call new_eval_gbsets_d3(ao_bas_han,
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), dbl_mb(id3chi_ao),
     &                 ibf_ao, 0, rq0, rq,
     &                 qxyz, xyz, nqpts, nqpts, mbf_ao,
     &                 natoms, dbl_mb(iscr), nscr, expo, NXprim,
     &                 ccoef, ncoef_max, 
     &                 iniz,ifin)
               else
                  call new_eval_gbsets(ao_bas_han, gradchi, 
     &                                 hesschi, dbl_mb(ichi_ao), 
     &                                 dbl_mb(idelchi_ao), 
     &                                 dbl_mb(iheschi_ao), 
     &                                 ibf_ao, 0, rq0,
     &                                 rq, qxyz, xyz, nqpts, 
     &                                 nqpts, mbf_ao, natoms, 
     .                 dbl_mb(iscr), 
     &                                 nscr, expo, NXprim, ccoef, 
     &                                 ncoef_max, 
     &                                 iniz,ifin)
               endif
               if (.not.ma_pop_stack(lscr))
     &            call errquit('dftg_quadv0a: cannot pop stack',0)
c         
c
c              Pre-compute reduced Xi(r) and dXi(r) over atoms
c                 reduce chi(nq,mbf_ao) to chi(natoms)
c                 reduce delchi(nq,3,mbf_ao) to delchi(natoms)
         call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,           dbl_mb(ichi_ao),rchi_atom)
         if(grad) 
     .        call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,           dbl_mb(idelchi_ao),rdelchi_atom)
c         
               if (.not.MA_Push_Get(mt_dbl, nqpts*natoms, 'bp', 
     &            lbp, ibp))call errquit
     &            (' exiting in quadvxc0: insuff stack',7)
         
               if (.not.MA_Push_Get(mt_dbl,nqpts,'rho scr',
     &            lrscr,irscr))call errquit
     &            ('dftg_quadv0a: cannot allocate rho scr',12)
c
               wmax= dabsmax(nqpts,qwght)
               call xc_rhogen(tol_rho, AO_bas_han, iga_dens, max_at_bf,
     &              natoms, ipol, nqpts, nbf_ao, mbf_ao,
     &              ldelley, grad, Fmat, Pmat, Bmat,
     &              dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &              ibf_ao, iniz, ifin, rho, delrho,
     &              dbl_mb(ibp), dbl_mb(irscr),
     &              rchi_atom, rdelchi_atom,
     &              rdens_atom, cetobfr,wmax)
c
               if (do_hess .or. do_cpks_r) then
                  call oactive_to_curatoms(oactive, natoms,
     &                 curatoms, ncuratoms)
                  call dfill(ipol*nqpts*3*nactive,0.D0,
     &                       dbl_mb(k_drhonuc),1)
                  if (grad)
     &               call dfill(ipol*3*nqpts*3*nactive,0.D0,
     &                          dbl_mb(k_ddelrhonuc),1)
c                  write(*,*)'calling xc_rhogen_nucder'
#if 1
                  call xc_rhogen_nucder(tol_rho, AO_bas_han, iga_dens,
     &                 max_at_bf, natoms, curatoms, ncuratoms, 
     &                 ipol, nqpts, nbf_ao, mbf_ao, grad, 
     &                 Fmat, Pmat, Bmat, Dmat,
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), ibf_ao, iniz, ifin,
     &                 dbl_mb(k_drhonuc), dbl_mb(k_ddelrhonuc),
     &                 dbl_mb(irscr), rchi_atom,
     &                 rdelchi_atom, rdens_atom, cetobfr)
#endif
               endif
c
               if (.not.ma_pop_stack(lrscr))
     &            call errquit('xc_quadv0_a: cannot pop stack',0)
c
               if (.not.ma_pop_stack(lbp))
     &            call errquit('dftg_quadv0a: cannot pop stack',0)
c         
c              Evaluate the integrated charge density.
c         
               nqtot = nqtot + nqpts

               rho_int = rho_int + 
     &                   ddot(nqpts, qwght, 1, rho, 1)

               call xc_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &              nqpts, Ex, Ec, qwght, GRAD, ldew, func,
     &              do_hess .or. do_cpks_r)
c     
c              exchange-correlation contribution, with or without gradient
c              terms
c     
               if (do_force) then
c 
!                 call dftg_quadv0_b(ictr_buf, .false., AO_bas_han, ipol,
                 call dftg_quadv0_b(ictr_buf, ldew, AO_bas_han, ipol,
     &                            grad, nbf_ao, mbf_ao, natoms, nqpts,
     &                            dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                            dbl_mb(iheschi_ao), Amat, Cmat,
     &                            ibf_ao, iniz, ifin,
     &                            ppp, ff, 
     &                            zzz, hh, iga_dens, force,
     &                            oactive, nactive)
c
               else if (do_hess) then
c
                  if (.not.MA_push_get(mt_dbl, 15*nqpts, 'd2scr',
     &                l_d2scr, k_d2scr)) call errquit
     &                ('dftg_quadv0: not enough stack for d2scr',0)
c
c     Note that k_d2scr is passed in three times, since it is also
c     used for the scratch arrays ff and s.  This is ok since by the
c     time ff and s are used, the routine is done with the array scr.
c -- what about using Bmat? !!!
                  call xc_d2expl(tol_rho, dbl_mb(k_d2scr), Amat, Amat2,
     &                 Bmat, Cmat, Cmat2, Dmat, Fmat, Pmat,
     &                 dbl_mb(k_d2scr), dbl_mb(k_d2scr+3*nqpts),
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), dbl_mb(id3chi_ao),
     &                 curatoms, ncuratoms, ipol, nqpts,
     &                 nbf_ao, mbf_ao, grad, AO_bas_han, natoms, iniz,
     &                 ifin, dbl_mb(k_drhonuc), dbl_mb(k_ddelrhonuc),
     &                 delrho, iga_dens, hess, ibf_ao,
     &                 rchi_atom, rdelchi_atom,
     &                 rdens_atom, cetobfr)
                  if (.not.ma_pop_stack(l_d2scr))
     &               call errquit('dftg_quadv0a: cannot pop stack',0)
c
               else if (do_cpks_r) then
c
c     The CPKS RHS is conveniently done in two parts:  the part
c     involving the derivs of the density looks just like the CPKS
c     LHS, so that code is called for that part, while the part
c     involving only the derivs of the basis functions is calculated
c     with a separate routine
c
c     First part: LHS-like contribution
c
c     Combine derivs of functional with deriv densities to get coeffs
c     of bfn pairs in CPKS LHS - these are written back into rho and
c     delrho
c
                  call xc_cpks_coeff(dbl_mb(k_drhonuc),
     &                 dbl_mb(k_ddelrhonuc), Amat2, Cmat2, Cmat, delrho,
     &                 3*ncuratoms, ipol, nqpts, GRAD)
#if 1
c                  write(*,*)'cpks first part'
                  call xc_tabcd_cpks_lhs(tol_rho, Fmat, Pmat,
     &                 dbl_mb(k_drhonuc), Bmat, dbl_mb(k_ddelrhonuc),
     &                 Dmat, Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao),
     &                 3*ncuratoms, ipol, nqpts, mbf_ao,
     &                 nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD,
     &                 AO_bas_han, natoms, iniz, ifin, g_rhs,
     &                 ibf_ao, rchi_atom,
     &                 rdelchi_atom, rdens_atom, cetobfr)
c                  write(*,*)'Cmat after lhs'
c                  call output(cmat,1,nq,1,6,nq,6,1)
#endif
#if 1
c
c     Second part: Basis function derivative contribution
c
c                  write(*,*)'cpks second part'
c
c     Transform Cmat to delrho form for use in xc_cpks_tabcd_rhs
c
                  if (grad) then
                     call transform_Cmat(delrho, Cmat, ipol, nqpts)
                  endif
                  call xc_tabcd_cpks_rhs(tol_rho, Fmat, Pmat, Amat,
     &                 Bmat, Cmat, Dmat, Dbl_MB(ichi_ao),
     &                 Dbl_MB(idelchi_ao), dbl_mb(iheschi_ao),
     &                 curatoms,
     &                 ncuratoms, ipol, nqpts, mbf_ao,
     &                 nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD,
     &                 AO_bas_han, natoms, iniz, ifin, g_rhs,
     &                 ibf_ao, rchi_atom,
     &                 rdelchi_atom, rdens_atom, cetobfr)
#endif
c
               endif
c
              if (.not.ma_pop_stack(ld3chi_ao))
     &            call errquit('dftg_quadv0a: cannot pop stack',0)
               if (.not.ma_pop_stack(lheschi_ao))
     &            call errquit('dftg_quadv0a: cannot pop stack',0)
               if (.not.ma_pop_stack(lgradchi_ao))
     &            call errquit('dftg_quadv0a: cannot pop stack',0)
               if (.not.ma_pop_stack(lchi_ao))
     &            call errquit('dftg_quadv0a: cannot pop stack',0)
c
               if (do_hess .or. do_cpks_r) then
                  if (grad) then
                     if (.not.ma_pop_stack(l_ddelrhonuc))
     &               call errquit('xc_quadv0_a: cannot pop ddelrhonuc',
     &                            0)
                  endif
                  if (.not.ma_pop_stack(l_drhonuc))
     &              call errquit('xc_quadv0_a: cannot pop drhonuc',0)
               endif
c
c
               if (do_force.and.ldew)then
                  call dftg_wderiv(ictr_buf, natoms, nqpts, nqpts, 
     .                 iandex, 
     &                             natoms, func, d_qwght,
     &                             force, oactive, nactive)
               endif
c         
c
      return 
c     
      end
