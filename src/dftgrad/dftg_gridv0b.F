      Subroutine dftg_gridv0b(nqpts,rad,ictr_buf,iga_dens,
     &                         ncoef_max, NXprim, max_at_bf,
     &                         rho_int, nqtot, iniz, ifin, rho, delrho,
     &                         func, d_qwght, qwght, qxyz, xyz, xyzm,
     &                         iandex, expo, ccoef, Amat, Amat2, Cmat,
     &                         Cmat2, Bmat, Dmat, Fmat, Pmat,  rq,
     &                         rq0, hh, 
     .     ibf_ao,bas_cent_info,bas_cset_info,docset,
     .     curatoms,rchi_atom,rdelchi_atom,ppp,ff,zzz,
     .     drhonuc,ddelrhonuc,d2scr,iscr,nscr,
     ,     force, hess, g_rhs, do_force,
     ,     nbf_ao_mxnbf_ce2,
     ,     acc_ao_gauss, ldew,grad,gradchi,hesschi,der3chi,
     &                         do_hess, do_cpks_r, natoms, oactive,
     &                         nactive, rdens_atom, 
     &                         cetobfr,donewbf)
c      
c     $Id: dftg_gridv0b.F,v 1.18 2002-03-22 02:20:55 sohirata Exp $
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
c !!! temporary !!!
#include "rtdb.fh"
      integer nqpts, ictr_buf ! [in]
      double precision rad ! [in]
      integer natoms
      integer ncoef_max, NXprim, max_at_bf
      integer nqtot
      integer nactive
      logical oactive(natoms)
      logical grad, gradchi, hesschi, der3chi
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
      double precision xyz(3,natoms),xyzm(3,natoms)
      integer iandex(natoms)
      integer iniz(natoms), ifin(natoms)
      integer ibf_ao(nbf_ao),curatoms(natoms)
      integer bas_cent_info(3,natoms),bas_cset_info(6,*)
      logical docset(*)
      double precision rchi_atom(natoms),rdelchi_atom(3,natoms)
      double precision ppp(*),ff(*),zzz(*)
      double precision ccoef(ncoef_max), expo(NXprim)
      double precision rq(n_per_rec,natoms)
      double precision rq0(natoms)
      double precision drhonuc(*),ddelrhonuc(*),d2scr(*),iscr(*)
      integer nscr ! [in]
c     
c     external
c
      double precision ddot
      integer nxtask
      External ddot, nxtask
c
c---> Matrix Elements of the XC Potential & Energy
c
      integer iga_dens(2)
c
c     local crap
c
      integer mbf_ao
      integer  m4
      integer   lheschi_ao, iheschi_ao,
     &        ld3chi_ao, id3chi_ao
      integer lchi_ao, ichi_ao, 
     &        lgradchi_ao, idelchi_ao
      integer nbf_ao_mxnbf_ce2, maxdbas
      double precision acc_AO_gauss 
      double precision Ex, Ec
      double precision rho_int ! [in/out]
      double precision dabsmax,wmax
      external dabsmax
      integer grid_nbfm
      external grid_nbfm
      integer natoms_l ! atoms left after xc_sig
c
      integer ncuratoms
      logical donewbf           !!! !!!! *** temp
      integer l_Acof2, k_Acof2, l_Ccof2, k_Ccof2
c     
c           Evaluate the AO basis set at each of the quad. points.
c           allocate arrays for exponents and contraction coefficients
c           within int_eval_gbsets mbf is # of orbs with nonzero 
c           points in sampling set.  iniz is a mapping array.
c           currently set up only for cartesian (6d) sets
c         
      Ex=0
      Ec=0
      if (der3chi) then
         maxdbas = 3
      elseif (hesschi) then
         maxdbas = 2
      elseif (gradchi) then
         maxdbas = 1
      else
         maxdbas = 0
      endif
c
            mbf_ao = nbf_ao
            call ifill(mbf_ao, 0, ibf_ao, 1)
            mbf_ao=grid_nbfm( AO_bas_han,  natoms, 
     &           ictr_buf,rad,xyz,qxyz,nqpts,
     &           ibf_ao, docset, iniz, ifin, expo, 
     .           minexp,ldiff, acc_AO_gauss,iatype_pt_chg)
            if (mbf_ao .eq. 0) return
            
            
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
               call qdist(rq0, rq, qxyz, xyz, nqpts, natoms)
      if (donewbf) then
               call xc_eval_basis(ao_bas_han, maxdbas, dbl_mb(ichi_ao),
     &              dbl_mb(idelchi_ao), dbl_mb(iheschi_ao),
     &              dbl_mb(id3chi_ao), rq, qxyz, xyz, nqpts, natoms,
     &              iniz, docset, bas_cent_info, bas_cset_info)
      else
c !!! Can delete this once confident in new basis fn code
               if (der3chi) then
                  call new_eval_gbsets_d3(ao_bas_han,
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), dbl_mb(id3chi_ao),
     &                 ibf_ao, 0, rq0, rq,
     &                 qxyz, xyz, nqpts, nqpts, mbf_ao,
     &                 natoms, iscr, nscr, expo, NXprim,
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
     .                 iscr, 
     &                                 nscr, expo, NXprim, ccoef, 
     &                                 ncoef_max, 
     &                                 iniz,ifin)
               endif
      endif
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
         
c
               wmax= dabsmax(nqpts,qwght)
               call xc_rhogen(tol_rho, AO_bas_han, iga_dens, max_at_bf,
     &              natoms, ipol, nqpts, nbf_ao, mbf_ao,
     &              ldelley, grad, Fmat, Pmat, Bmat,
     &              dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &              ibf_ao, iniz, ifin, rho, delrho,
     &              rchi_atom, rdelchi_atom,
     &              rdens_atom, cetobfr,wmax)
c
               if (do_hess .or. do_cpks_r) then
                  call oactive_to_curatoms(oactive, natoms,
     &                 curatoms, ncuratoms)
                  call dfill(ipol*nqpts*3*nactive,0.D0,
     &                       drhonuc,1)
                  if (grad)
     &               call dfill(ipol*3*nqpts*3*nactive,0.D0,
     &                          ddelrhonuc,1)
c                  write(*,*)'calling xc_rhogen_nucder'
#if 1
                  call xc_rhogen_nucder(tol_rho, AO_bas_han, iga_dens,
     &                 max_at_bf, natoms, curatoms, ncuratoms, 
     &                 ipol, nqpts, nbf_ao, mbf_ao, grad, 
     &                 Fmat, Pmat, Bmat, Dmat,
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), ibf_ao, iniz, ifin,
     &                 drhonuc, ddelrhonuc,
     &                 rchi_atom,
     &                 rdelchi_atom, rdens_atom, cetobfr,1d0)
#endif
               endif
c
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
                 call dftg_quadv0_b(ictr_buf, ldew, AO_bas_han, ipol,
!                 call dftg_quadv0_b(ictr_buf, .false., AO_bas_han, ipol,
     &                            grad, nbf_ao, mbf_ao, natoms, nqpts,
     &                            dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                            dbl_mb(iheschi_ao), Amat, Cmat,
     &                            ibf_ao, iniz, ifin,
     &                            ppp, ff, 
     &                            zzz, hh, iga_dens, force,
     &                            oactive)
c
               else if (do_hess) then
c
c     Note that k_d2scr is passed in three times, since it is also
c     used for the scratch arrays ff and s.  This is ok since by the
c     time ff and s are used, the routine is done with the array scr.
c
c     Allocate scratch space to use with xc_cpks_coeff in xc_d2expl
c
                  if (.not.MA_push_get(MT_DBL, ipol*nqpts*3,
     &                 'Acof2', l_Acof2, k_Acof2)) call
     &                 errquit('dftg_gridv0b: cannot allocate Acof2',1)
                  if (grad) then
                     if (.not.MA_push_get(MT_DBL, ipol*9*nqpts, 'CCof2',
     &                    l_Ccof2, k_Ccof2)) call
     &                    errquit('dftg_gridv0b: cant alloc Ccof2',1)
                  endif
                  call xc_d2expl(tol_rho, d2scr(3*nqpts+1), Amat, Amat2,
     &                 dbl_mb(k_Acof2), Cmat, Cmat2, dbl_mb(k_Ccof2),
     &                 Fmat, Pmat, d2scr, d2scr(3*nqpts+1),
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), dbl_mb(id3chi_ao),
     &                 curatoms, ncuratoms, ipol, nqpts,
     &                 nbf_ao, mbf_ao, grad, AO_bas_han, natoms, iniz,
     &                 ifin, drhonuc, ddelrhonuc,
     &                 delrho, iga_dens, hess, ibf_ao,
     &                 rchi_atom, rdelchi_atom,
     &                 rdens_atom, cetobfr)
                  if (grad) then
                     if (.not.ma_pop_stack(l_Ccof2))
     &                    call errquit('dftg_gridv0b: cannot pop CCof2',
     &                    0)
                  endif
                  if (.not.ma_pop_stack(l_Acof2))
     &                 call errquit('dftg_gridv0b: cannot pop Acof2',0)
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
                  call xc_cpks_coeff(drhonuc,
     &                 ddelrhonuc, Amat2, Cmat2, Cmat, delrho,
     &                 3*ncuratoms, ipol, nqpts, GRAD, .false.)
#if 1
                  call xc_tabcd_cpks_lhs(tol_rho, Fmat, Pmat,
     &                 drhonuc, Bmat, ddelrhonuc,
     &                 Dmat, Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao),
     &                 3*ncuratoms, ipol, nqpts, mbf_ao,
     &                 nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, GRAD,
     &                 AO_bas_han, natoms, iniz, ifin, g_rhs,
     &                 ibf_ao, rchi_atom,
     &                 rdelchi_atom, rdens_atom, cetobfr)
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
 2010          continue
c
               if (.not.ma_chop_stack(lchi_ao))
     &            call errquit('dftg_gridv0b: cannot chop stack',0)
c
c
               if (do_force.and.ldew)then
                  natoms_l=natoms
                  call dftg_gridpart(d_qwght,qwght, qxyz, xyz,xyzm, 
     .                 rq, ictr_buf, natoms_l, iandex,
     .                 nqpts,whichssw,
     ,                 rad,expo)
                  call dftg_wderiv(ictr_buf, natoms, nqpts, nqpts, 
     .                 iandex, 
     &                 natoms_l, func, d_qwght,
     &                 force, oactive)
               endif
c         
c
      return 
c     
      end
