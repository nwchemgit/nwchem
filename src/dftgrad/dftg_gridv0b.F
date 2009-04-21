      subroutine dftg_gridv0b(nqpts,rad,ictr_buf,iga_dens,
     &     max_at_bf,
     &     rho_int, nqtot, iniz, ifin, rho, delrho,
     &     func, d_qwght, qwght, qxyz, xyz, xyzm,
     &     iandex, expo, Amat, Amat2, Cmat,
     &     Cmat2, Bmat, Fmat, Pmat, rq,
     &     rq0, hh, 
     &     ibf_ao,bas_cent_info,bas_cset_info,docset,
     &     curatoms,rchi_atom,rdelchi_atom,ppp,ff,zzz,
     &     drhonuc,ddelrhonuc,d2scr,
     &     force, hess, g_rhs, do_force,
     &     nbf_ao_mxnbf_ce2,
     &     acc_ao_gauss, ldew,grad,gradchi,hesschi,der3chi,
     &     do_hess, do_cpks_r, natoms, oactive,
     &     nactive, rdens_atom, 
     &     cetobfr, tau, dttau, kske, Mmat, Mmat2,
     &     ldbw, fact, dtype, nconstr)  ! Qin
c      
c     $Id: dftg_gridv0b.F,v 1.39 2008-08-14 06:31:58 niri Exp $
c     
      implicit none
#include "errquit.fh"
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
#include "dftps.fh"

      logical ldbw
      integer ispin, i, ii, nconstr
      double precision fact(nconstr)
      integer dtype(nconstr)          ! Qin
c
      integer nqpts, ictr_buf ! [in]
      double precision rad ! [in]
      integer natoms
      integer max_at_bf
      integer nqtot
      integer nactive
      logical oactive(natoms)
      logical grad, gradchi, hesschi, der3chi, kske
      logical ldew
      logical do_force, do_hess, do_cpks_r
c
      integer cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      double precision Force(3*natoms), hess(3*natoms,3*natoms)
      integer g_rhs(*)
      double precision hh(*)
      double precision rho(n_per_rec,ipol*(ipol+1)/2)
      double precision delrho(n_per_rec,3,ipol)
c
      double precision tau(n_per_rec,3,ipol)
      double precision dttau(*)
      double precision Mmat(*)
c
      double precision Amat(n_per_rec,ipol), Cmat(n_per_rec,3,ipol),
     &                 Bmat(n_per_rec,*)
c
      double precision Amat2(n_per_rec,NCOL_AMAT2)
      double precision Cmat2(n_per_rec,NCOL_CMAT2)
      double precision Mmat2(n_per_rec,NCOL_MMAT2)
c
      double precision Pmat(max_at_bf*max_at_bf),
     &                 Fmat(max_at_bf*max_at_bf)
c
      double precision func(n_per_rec)
      double precision d_qwght(3,n_per_rec,*)
      double precision qwght(*), qxyz(3,*)
      double precision xyz(3,natoms),xyzm(3,natoms)
      integer iandex(natoms)
      integer iniz(natoms), ifin(natoms)
      integer ibf_ao(nbf_ao),curatoms(natoms)
      integer bas_cent_info(3,natoms),bas_cset_info(6,*)
      logical docset(*)
      double precision rchi_atom(natoms),rdelchi_atom(natoms)
      double precision ppp(*),ff(*),zzz(*)
      double precision expo(*)
      double precision rq(n_per_rec,natoms)
      double precision rq0(natoms)
      double precision drhonuc(*),ddelrhonuc(*),d2scr(*)
c     
c     external
c
      double precision ddot
      External ddot
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
      double precision StericEnergy
      double precision rho_int ! [in/out]
      double precision dabsmax,wmax
      external dabsmax
      integer grid_nbfm
      external grid_nbfm
      integer natoms_l ! atoms left after xc_sig
c
      integer ncuratoms
      integer l_Acof2, k_Acof2
      integer l_Ccof2, k_Ccof2
      integer l_Mcof2, k_Mcof2
      integer i_ccdel,l_ccdel,idmat,ldmat
      logical do_2nd
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
     &            (' exiting in dftg_gridv0b: insuff stack',9, MA_ERR)
               if (.not. MA_Push_Get(MT_Dbl, 3*nqpts*mbf_ao, 
     &            'gradchi_ao', lgradchi_ao,idelchi_ao))call errquit
     &            (' exiting in dftg_gridv0b: insuff stack',9, MA_ERR)
c         
c              hessian is used only for derivatives; set array to length 1 here
c
               if (hesschi)then
                  m4 = 6*nqpts*mbf_ao
               else
                  m4 = 1
               endif
c
               if (.not. MA_Push_Get(MT_Dbl, m4, 'heschi_ao', 
     &            lheschi_ao, iheschi_ao))call errquit
     &            (' exiting in dftg_gridv0b: insuff stack',9, MA_ERR)
c    
               if (der3chi)then
                  m4 = 10*nqpts*mbf_ao
               else
                  m4 = 1
               endif
c
               if (.not. MA_Push_Get(MT_Dbl, m4, 'd3chi_ao', 
     &            ld3chi_ao, id3chi_ao))call errquit
     &            (' exiting in dftg_gridv0b: insuff stack',19, MA_ERR)
c
               if (odftps) call pstat_on(ps_excf)
               call qdist(rq0, rq, qxyz, xyz, nqpts, natoms)
               call xc_eval_basis(ao_bas_han, maxdbas, dbl_mb(ichi_ao),
     &              dbl_mb(idelchi_ao), dbl_mb(iheschi_ao),
     &              dbl_mb(id3chi_ao), rq, qxyz, xyz, nqpts, natoms,
     &              iniz, docset, bas_cent_info, bas_cset_info)
               if (odftps) call pstat_off(ps_excf)
c         
c
c              Pre-compute reduced Xi(r) and dXi(r) over atoms
c                 reduce chi(nq,mbf_ao) to chi(natoms)
c                 reduce delchi(nq,3,mbf_ao) to delchi(natoms)
         call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,           dbl_mb(ichi_ao),rchi_atom)
         if(grad.or.kske) 
     .        call util_rmsatbf(nqpts*3, natoms,iniz,ifin,
     ,           dbl_mb(idelchi_ao),rdelchi_atom)
c         
      if (do_hess .or. do_cpks_r.or.kske) then
         if (.not.MA_Push_Get(mt_dbl,nqpts*3*nbf_ao_mxnbf_ce,
     A        'Dmat',lDmat,iDmat))
     &        call errquit(' exiting in dftg_gridv0: insuff stack',131,
     &       MA_ERR)
      endif
c
               wmax= dabsmax(nqpts,qwght)
               if (odftps) call pstat_on(ps_xcrho)
               call xc_rhogen(0,
     T              tol_rho, AO_bas_han, iga_dens, max_at_bf,
     &              natoms, 1,1,1,
     I              ipol, nqpts, nbf_ao, mbf_ao, grad,
     P              (ipol*(ipol+1))/2,
     F              Fmat, Pmat, Bmat, dbl_mb(iDmat),
     &              dbl_mb(ichi_ao), dbl_mb(idelchi_ao),0d0,
     &              ibf_ao, iniz, ifin, rho, delrho,
     &              rchi_atom, rdelchi_atom,
     &              rdens_atom, cetobfr,wmax, 
     &              tau, dttau, kske)
               if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'
c
               if (do_hess .or. do_cpks_r) then
                  call oactive_to_curatoms(oactive, natoms,
     &                 curatoms, ncuratoms)
                  call dfill(ipol*nqpts*3*nactive,0.D0,
     &                       drhonuc,1)
                  if (grad)
     &               call dfill(ipol*3*nqpts*3*nactive,0.D0,
     &                          ddelrhonuc,1)
                  call xc_rhogen(2,
     &                 tol_rho, AO_bas_han, iga_dens,max_at_bf, 
     N                 natoms, curatoms, ncuratoms, 1,
     &                 ipol, nqpts, nbf_ao, mbf_ao, grad, ipol,
     &                 Fmat, Pmat, Bmat, dbl_mb(iDmat),
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), ibf_ao, iniz, ifin,
     &                 drhonuc, ddelrhonuc,
     &                 rchi_atom, rdelchi_atom, 
     &                 rdens_atom, cetobfr,1d0,
     D                 tau, dttau, kske)
               endif
               if (odftps) call pstat_off(ps_xcrho)
c
c              Evaluate the integrated charge density.
c         
               nqtot = nqtot + nqpts
               rho_int = rho_int + 
     &                   ddot(nqpts, qwght, 1, rho, 1)

c Qin for lcdft and pop1
c the idea is to add qwght to Amat, which is dfunc and will time drho,
c and rho to func, which will time dweight
c to count for the contribution of Becke population of contrained atoms

               if(ldbw) then
                do i = 1, nconstr
                 do ispin = 1, ipol
                   ii = (ispin-1)*nqpts + 1
                   if(dtype(i).eq.2.and.ispin.eq.2) then
                    call daxpy(nqpts, -fact(i), qwght, 1, Amat(ii,1),1)
                   else
                    call daxpy(nqpts, fact(i), qwght, 1, Amat(ii,1), 1)
                   endif
                 enddo
                 if(dtype(i).eq.1) then
                  call daxpy(nqpts, fact(i), rho(1,1), 1, func, 1)
                 else
                  call daxpy(nqpts, fact(i), rho(nqpts+1,1), 1, func, 1)
                  call daxpy(nqpts,-fact(i), rho(2*nqpts+1,1),1,func,1)
                 endif
                enddo
               endif
c Qin end
               if (odftps) call pstat_on(ps_ebf)
               do_2nd = (do_hess .or. do_cpks_r)
               call xc_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &            nqpts, Ex, Ec, qwght, GRAD, ldew, func, do_2nd,
     &            tau, dttau, kske, Mmat, Mmat2,
     &            StericEnergy)
               if (odftps) call pstat_off(ps_ebf)
c     
c              exchange-correlation contribution, with or without gradient terms
c     
               if (do_force) then
c 
                  if(grad) then
                     if (.not. MA_Push_Get(MT_Dbl, nqpts*mbf_ao*ipol, 
     &                    'gradchi_ao2', l_ccdel,i_ccdel))call errquit
     &                    (' exit dftggridv0b: insuff stack',9, 0)
                  endif

               if (odftps) call pstat_on(ps_qv0b)
               call dftg_quadv0_b(ictr_buf, ldew, AO_bas_han, ipol,
     &                            grad, kske, nbf_ao, mbf_ao, natoms, 
     &                            nqpts,
     &                            dbl_mb(ichi_ao), dbl_mb(idelchi_ao ),
     &                            dbl_mb(iheschi_ao), Amat, Cmat, Mmat,
     &                            dbl_mb(i_ccdel),
     &                            ibf_ao, iniz, ifin, cetobfr,
     &                            ppp, ff, 
     &                            zzz, hh, iga_dens, force,
     &                            oactive)
               if (odftps) call pstat_off(ps_qv0b)
                 if(grad) then
                    if (.not.ma_pop_stack(l_ccdel))
     &                   call errquit('dftggv0b: pop stack ',1,0)
                    endif
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
     &                 errquit('dftg_gridv0b: cannot allocate Acof2',1,
     &       MA_ERR)
                  if (grad) then
                     if (.not.MA_push_get(MT_DBL, ipol*9*nqpts, 'CCof2',
     &                    l_Ccof2, k_Ccof2)) call
     &                    errquit('dftg_gridv0b: cant alloc Ccof2',1,
     &       MA_ERR)
                  endif
                  if (kske) then
                     if (.not.MA_push_get(MT_DBL, ipol*3*nqpts, 'MCof2',
     &                    l_Mcof2, k_Mcof2)) call
     &                    errquit('dftg_gridv0b: cant alloc Mcof2',1,
     &       MA_ERR)
                  end if
c
                  call xc_d2expl(tol_rho, 
     D                 d2scr(3*nqpts*nbf_ao_mxnbf_ce+1),
     A                 Amat, Amat2, dbl_mb(k_Acof2), 
     &                 Cmat, Cmat2, dbl_mb(k_Ccof2),
     &                 Mmat, Mmat2, dbl_mb(k_Mcof2),
     &                 Fmat, Pmat, d2scr, 
     D                 d2scr(3*nqpts*nbf_ao_mxnbf_ce+1),
     &                 dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &                 dbl_mb(iheschi_ao), dbl_mb(id3chi_ao),
     &                 curatoms, ncuratoms, ipol, nqpts,
     &                 nbf_ao, mbf_ao, grad, AO_bas_han, natoms, iniz,
     &                 ifin, drhonuc, ddelrhonuc, dttau,
     &                 delrho, iga_dens, hess, ibf_ao,
     &                 rchi_atom, rdelchi_atom,
     &                 rdens_atom, cetobfr, kske)
c
                  if (grad) then
                     if (.not.ma_pop_stack(l_Ccof2))
     &                    call errquit('dftg_gridv0b: cannot pop CCof2',
     &                    0, MA_ERR)
                  endif
                  if (kske) then
                     if (.not.ma_pop_stack(l_Mcof2))
     &                    call errquit('dftg_gridv0b: cannot pop CCof2',
     &                    0, MA_ERR)
                  endif
                  if (.not.ma_pop_stack(l_Acof2))
     &                 call errquit('dftg_gridv0b: cannot pop Acof2',0,
     &       MA_ERR)
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
                  call xc_cpks_coeff(drhonuc, ddelrhonuc, dttau,
     &                 Amat2, Cmat2, Cmat, Mmat2, delrho,
     &                 3*ncuratoms, ipol, nqpts, grad, .false., ! kske = false 
     &                 .false.)
c
                  if (odftps) call pstat_on(ps_tabcd)
                  call xc_tabcd(1,.true.,
     T                 tol_rho, Fmat, Pmat, drhonuc, Bmat, ddelrhonuc,
     &                 dbl_mb(iDmat), 
     N                 0d0,0d0,0d0,0d0,
     C                 Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     &                 0d0,0,3*ncuratoms, 
     I                 ipol,nqpts,mbf_ao,
     N                 nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, 
     G                 GRAD,AO_bas_han, natoms, iniz, ifin, 
     G                 g_rhs, ibf_ao, rchi_atom, rdelchi_atom, 
     R                 rdens_atom, cetobfr,
     &                 .false., dttau, 0d0)  ! kske = false
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
                  call xc_tabcd(2,.true.,
     &                 tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, 
     D                 dbl_mb(iDmat), 
     N                 0d0,0d0,0d0,0d0,
     C                 Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 
     H                 dbl_mb(iheschi_ao),
     &                 curatoms, ncuratoms, 1,
     I                 ipol, nqpts, mbf_ao,
     &                 nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, 
     G                 GRAD, AO_bas_han, natoms, iniz, ifin, 
     G                 g_rhs, ibf_ao, rchi_atom, rdelchi_atom, 
     R                 rdens_atom, cetobfr,
     &                 .false.,Mmat,0d0)  ! kske = false
c
                  if (odftps) call pstat_off(ps_tabcd)
               endif
c
               if (.not.ma_chop_stack(lchi_ao))
     &            call errquit('dftg_gridv0b: cannot chop stack',0,
     &       MA_ERR)
c
c
               if (do_force.and.ldew)then
                  natoms_l=natoms
                  if (odftps) call pstat_on(ps_ssw)                  
                  call dftg_gridpart(d_qwght,qwght, qxyz, xyz,xyzm, 
     .                 rq, ictr_buf, natoms_l, iandex,
     .                 nqpts,whichssw,
     ,                 rad,expo)
                  call dftg_wderiv(ictr_buf, natoms, nqpts, nqpts, 
     .                 iandex, 
     &                 natoms_l, func, d_qwght,
     &                 force, oactive)
                  if (odftps) call pstat_off(ps_ssw)
               endif
c         
      return 
      end
