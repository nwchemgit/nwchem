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
     &     cetobfr, dttau, kske, Mmat, Mmat2,
     &     kslap, Lmat,
     &     ldbw, fact, dtype, nconstr,  ! Qin
     1     do_gxc, Amat3, Cmat3, ngxc, 
     2     crho, cdelrho, triplet,
     3     do_tddftvxc,cd_coef)
c      
c     $Id$
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
c Daniel (1-14-13): Added 3rd derivative header file
#include "dft3drv.fh"
c !!! BGJ test !!!
#include "bgj.fh"
c !!! BGJ test !!!
c !!! temporary !!!
#include "rtdb.fh"
#include "dftps.fh"
#include "util.fh"

      logical ldbw
      integer ispin, i, ii, nconstr
      double precision fact(nconstr)
      integer dtype(nconstr)          ! Qin
      double precision cd_coef(nbf_cd,*)
c
      integer nqpts, ictr_buf ! [in]
      double precision rad ! [in]
      integer natoms
      integer max_at_bf
      integer nqtot
      integer nactive
      logical oactive(natoms)
      logical grad, gradchi, hesschi, der3chi, kske, kslap
      logical ldew
      logical do_force, do_hess, do_cpks_r
c
      integer cetobfr(2,natoms)
      double precision rdens_atom(ipol*natoms*natoms)
      double precision Force(3*natoms), hess(3*natoms,3*natoms)
      integer g_rhs(*)
      double precision hh(*)
      double precision lapdum
c Daniel (2-12-13): Nothing appears to be broken because of this
c change. 
c      double precision rho(n_per_rec,ipol*(ipol+1)/2)
c      double precision delrho(n_per_rec,3,ipol)
      double precision rho(*)
      double precision delrho(*)
c
      double precision dttau(*)
      double precision Mmat(nqpts,ipol)
c
      double precision Lmat(nqpts,ipol)
c
      double precision Amat(nqpts,ipol), Cmat(nqpts,3,ipol),
     &                 Bmat(nqpts,*)
c
      double precision Amat2(nqpts,NCOL_AMAT2)
      double precision Cmat2(nqpts,NCOL_CMAT2)
      double precision Mmat2(nqpts,NCOL_MMAT2)
c
      double precision Pmat(max_at_bf*max_at_bf),
     &                 Fmat(max_at_bf*max_at_bf)
c
      double precision func(nqpts)
      double precision d_qwght(3,nqpts,*)
      double precision qwght(*), qxyz(3,*)
      double precision xyz(3,natoms),xyzm(3,natoms)
      integer iandex(natoms)
      integer iniz(natoms), ifin(natoms)
      integer ibf_ao(*),curatoms(natoms)
      integer bas_cent_info(3,natoms),bas_cset_info(6,*)
      logical docset(*)
      double precision rchi_atom(natoms),rdelchi_atom(natoms)
      double precision ppp(*),ff(*),zzz(*)
      double precision expo(*)
      double precision rq(nqpts,natoms)
      double precision rq0(natoms)
      double precision drhonuc(*),ddelrhonuc(*),d2scr(*)
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
      integer llap, ilap
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
c Daniel (1-14-13): Define XC-third derivatives for evaluating the 
c XC-kernel gradients.
      logical do_gxc, do_3rd
      double precision Amat3(nqpts,NCOL_AMAT3)
      double precision Cmat3(nqpts,NCOL_CMAT3)
c Daniel (1-17-13): Other required stuff for the current implementation
c of the TDDFT gradients.  The "coefficient rho" and "coefficient
c delrho" are created because we can't overwrite the perturbed density
c matrix as we do for linear response calculations.
      integer ngxc
      integer iscf_rho, iscf_delrho
c      double precision crho(n_per_rec,ipol,ngxc)
c      double precision cdelrho(n_per_rec,3,ipol,ngxc)
      double precision crho(*)
      double precision cdelrho(*)
c Daniel (1-24-13): This is needed for evaluating the grid contribution
c to the dVxc term in the TDDFT gradients
      logical do_pertdens
c Daniel (2-6-13): Logical variable for doing restricted triplet 
c TDDFT gradients from the XC-kernel.
      logical triplet
c Daniel (2-7-13): Logical variable storing if we're calculating the
c dVxc*P contribution.
      logical do_tddftvxc
      logical ldew2 ! Prevent weighting of derivative matrices

      double precision ddot
      external ddot

      integer basis
c     
c           Evaluate the AO basis set at each of the quad. points.
c           allocate arrays for exponents and contraction coefficients
c           within int_eval_gbsets mbf is # of orbs with nonzero 
c           points in sampling set.  iniz is a mapping array.
c           currently set up only for cartesian (6d) sets
c         
c Daniel (2-14-13): Fixed the indentation in this routine.
      iDmat = 0
      ichi_ao = 0
      idelchi_ao = 0
      iheschi_ao = 0
      id3chi_ao = 0
      i_ccdel = 0
      k_Acof2 = 0
      k_Ccof2 = 0
      k_Mcof2 = 0

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

      if (adft) then
        basis = cd_bas_han
        mbf_ao = nbf_cd
      else
        basis = ao_bas_han
        mbf_ao = nbf_ao
      endif
c
c Daniel (1-14-13): Fixed the indenting in this routine because it
c made the code difficult to read.
      call ifill(mbf_ao, 0, ibf_ao, 1)
      mbf_ao=grid_nbfm( basis,  natoms, 
     &     ictr_buf,rad,xyz,qxyz,nqpts,
     &     ibf_ao, docset, iniz, ifin, expo, 
     .     minexp,ldiff, acc_AO_gauss,iatype_pt_chg)
      if (mbf_ao .eq. 0) return
c
      if (.not. MA_Push_Get(MT_Dbl, nqpts*mbf_ao, 'chi_ao',
     &   lchi_ao, ichi_ao))call errquit
     &   (' exiting in dftg_gridv0b: insuff stack',9, MA_ERR)
      if (.not. MA_Push_Get(MT_Dbl, 3*nqpts*mbf_ao, 
     &   'gradchi_ao', lgradchi_ao,idelchi_ao))call errquit
     &   (' exiting in dftg_gridv0b: insuff stack',9, MA_ERR)
c     
c     hessian is used only for derivatives; set array to length 1 here
c
      if (hesschi)then
         m4 = 6*nqpts*mbf_ao
      else
         m4 = 1
      endif
c
      if (.not. MA_Push_Get(MT_Dbl, m4, 'heschi_ao', 
     &   lheschi_ao, iheschi_ao))call errquit
     &   (' exiting in dftg_gridv0b: insuff stack',9, MA_ERR)
c    
      if (der3chi)then
         m4 = 10*nqpts*mbf_ao
      else
         m4 = 1
      endif
c
      if (.not. MA_Push_Get(MT_Dbl, m4, 'd3chi_ao', 
     &   ld3chi_ao, id3chi_ao))call errquit
     &   (' exiting in dftg_gridv0b: insuff stack',19, MA_ERR)
c
      if (odftps) call pstat_on(ps_excf)
      call qdist(rq0, rq, qxyz, xyz, nqpts, natoms)
      call xc_eval_basis(basis, maxdbas, dbl_mb(ichi_ao),
     &     dbl_mb(idelchi_ao), dbl_mb(iheschi_ao),
     &     dbl_mb(id3chi_ao), rq, qxyz, xyz, nqpts, natoms,
     &     iniz, docset, bas_cent_info, bas_cset_info)
      if (odftps) call pstat_off(ps_excf)
c      
c
c           Pre-compute reduced Xi(r) and dXi(r) over atoms
c              reduce chi(nq,mbf_ao) to chi(natoms)
c              reduce delchi(nq,3,mbf_ao) to delchi(natoms)
      call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,        dbl_mb(ichi_ao),rchi_atom)
      if(grad.or.kske) 
     .     call util_rmsatbf(nqpts*3, natoms,iniz,ifin,
     ,        dbl_mb(idelchi_ao),rdelchi_atom)
c      
      if (do_hess .or. do_cpks_r.or.kske) then
         if (.not.MA_Push_Get(mt_dbl,nqpts*3*nbf_ao_mxnbf_ce,
     A        'Dmat',lDmat,iDmat))
     &        call errquit(' exiting in dftg_gridv0: insuff stack',
     &       131, MA_ERR)
      endif
c
      if (kslap) then
         m4 = nqpts*ipol*(ipol+1)/2
      else
         m4 = 1
      endif

      if (.not.MA_Push_Get(mt_dbl,m4,'lap',llap,ilap))
     &     call errquit(' exiting in dftg_gridv0b: insuff stack',
     &     0, MA_ERR)
c
      wmax= dabsmax(nqpts,qwght)
      if (odftps) call pstat_on(ps_xcrho)
c Daniel (1-16-13): Avoid this part when we do the third derivatives.
c      do_pertdens = do_gxc
      do_pertdens = (do_gxc .or. do_tddftvxc)
      if (.not. do_pertdens) then
        if (adft) then
          call xc_rhogen_adft(tol_rho,natoms,ipol,nqpts,
     $         nbf_cd,mbf_ao,max_at_bf,grad,(ipol*(ipol+1))/2,
     $         dbl_mb(ichi_ao),dbl_mb(idelchi_ao),dbl_mb(iheschi_ao),
     $         ibf_ao, iniz, ifin, rho, delrho, dbl_mb(ilap), rchi_atom,
     $         rdelchi_atom, cetobfr, wmax, kslap, cd_coef)
        else
          call xc_rhogen(0,
     T       tol_rho, AO_bas_han, iga_dens, max_at_bf,
     &       natoms, 1,1,1,
     I       ipol, nqpts, nbf_ao, mbf_ao, grad,
     P       (ipol*(ipol+1))/2,
     F       Fmat, Pmat, Bmat, dbl_mb(iDmat),
     &       dbl_mb(ichi_ao), dbl_mb(idelchi_ao), dbl_mb(iheschi_ao),
     &       ibf_ao, iniz, ifin, rho, delrho, dbl_mb(ilap),
     &       rchi_atom, rdelchi_atom,
     &       rdens_atom, cetobfr,wmax, 
     &       dttau, kske, kslap)
        endif
      else
c Daniel (1-16-13): Here, the call is identical to the call to 
c xc_rhogen in grid_quadv0b.
        call xc_rhogen(1,
     T    tol_rho, AO_bas_han, iga_dens, max_at_bf,
     N    natoms, 1, 1, ngxc+1,
     I    ipol, nqpts, nbf_ao, mbf_ao, grad, ipol,
     F    Fmat, Pmat, Bmat, dbl_mb(iDmat),
     &    dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 0d0,
     &    ibf_ao, iniz, ifin, rho, delrho, lapdum,
     &    rchi_atom, rdelchi_atom,
     &    rdens_atom, cetobfr,wmax,
     &    dttau, kske, .false.)
      endif
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2e'
c
c Daniel (2-7-13): Presumably, this is constructing the gradient of
c the ground state density along all atoms.
      if (do_hess .or. do_cpks_r) then
         call oactive_to_curatoms(oactive, natoms,
     &        curatoms, ncuratoms)
         call dfill(ipol*nqpts*3*nactive,0.D0,
     &              drhonuc,1)
         if (grad)
     &      call dfill(ipol*3*nqpts*3*nactive,0.D0,
     &                 ddelrhonuc,1)
         call xc_rhogen(2,
     &        tol_rho, AO_bas_han, iga_dens,max_at_bf, 
     N        natoms, curatoms, ncuratoms, 1,
     &        ipol, nqpts, nbf_ao, mbf_ao, grad, ipol,
     &        Fmat, Pmat, Bmat, dbl_mb(iDmat),
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &        dbl_mb(iheschi_ao), ibf_ao, iniz, ifin,
     &        drhonuc, ddelrhonuc, lapdum,
     &        rchi_atom, rdelchi_atom, 
     &        rdens_atom, cetobfr,1d0,
     D        dttau, kske, .false.)
      endif
      if (odftps) call pstat_off(ps_xcrho)
c
c              Evaluate the integrated charge density.
c         
      nqtot = nqtot + nqpts
      rho_int = rho_int + 
     &          ddot(nqpts, qwght, 1, rho, 1)

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
c         call daxpy(nqpts, fact(i), rho(1,1), 1, func, 1)
         call daxpy(nqpts, fact(i), rho(1), 1, func, 1)
        else
c         call daxpy(nqpts, fact(i), rho(nqpts+1,1), 1, func, 1)
         call daxpy(nqpts, fact(i), rho(nqpts+1), 1, func, 1)
c         call daxpy(nqpts,-fact(i), rho(2*nqpts+1,1),1,func,1)
         call daxpy(nqpts,-fact(i), rho(2*nqpts+1),1,func,1)
        endif
       enddo
      endif
c Qin end
c Daniel (1-14-13): Here we mimic what was done in grid_quadv0b to
c make sure the correct density is used to construct the XC-functional
c derivatives
      iscf_rho = 1
      iscf_delrho = 1
      if (do_pertdens) then
        iscf_rho = 1 + nqpts*ipol*ngxc
        iscf_delrho = 1 + nqpts*3*ipol*ngxc
      endif
c
      if (do_pertdens.and.ipol.eq.2) then
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
      if (odftps) call pstat_on(ps_ebf)
c      do_2nd = (do_hess .or. do_cpks_r)
      do_2nd = (do_hess .or. do_cpks_r .or. do_tddftvxc)
c Daniel (1-14-13): The call to xc_eval_fnl is adapted here for
c XC-kernel gradients in TDDFT.
      do_3rd = do_gxc
c Daniel (2-7-13): This is used to prevent the code from weighting the
c functional derivative matrices with qwght in xc_eval_fnl.  We do this
c to reduce the work we need to do later (i.e. weighting the derivatives
c and then un-weighting them to get the quadrature weight gradients). 
      ldew2 = do_tddftvxc
      call xc_eval_fnl(rho(iscf_rho), delrho(iscf_delrho), 
     &   Amat, Amat2, Cmat, Cmat2,
     &   nqpts, Ex, Ec, qwght, GRAD, ldew, func, do_2nd,
     &   dttau, kske, Mmat, Mmat2,
     &   dbl_mb(ilap),kslap,Lmat,
     &   StericEnergy,
     &   do_3rd, Amat3, Cmat3, ldew2)
c      call xc_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
c     &   nqpts, Ex, Ec, qwght, GRAD, ldew, func, do_2nd,
c     &   dttau, kske, Mmat, Mmat2,
c     &   StericEnergy,
c     &   do_3rd, Amat3, Cmat3)
      if (odftps) call pstat_off(ps_ebf)
c     
c              exchange-correlation contribution, with or without gradient terms
c     
      if (do_force) then
c 
         if (grad) then
            if (.not. MA_Push_Get(MT_Dbl,
     $           max(10*nqpts,nqpts*mbf_ao*ipol), 
     &           'gradchi_ao2', l_ccdel,i_ccdel))call errquit
     &           (' exit dftggridv0b: insuff stack',9, 0)
         endif

         if (odftps) call pstat_on(ps_qv0b)
         if (adft) then
           call dftg_quadv0_b_adft(ictr_buf,ldew,basis,ipol,grad,
     $          kslap,nbf_cd,mbf_ao,natoms,nqpts,dbl_mb(ichi_ao),
     $          dbl_mb(idelchi_ao),dbl_mb(iheschi_ao),
     $          dbl_mb(id3chi_ao),Amat,Cmat,Lmat,dbl_mb(i_ccdel),
     $          ibf_ao,iniz,ifin,cetobfr,force,oactive,cd_coef)
         else
           call dftg_quadv0_b(ictr_buf, ldew, AO_bas_han, ipol,
     &                      grad, kske, kslap, nbf_ao, mbf_ao, natoms,
     &                      nqpts,
     &                      dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
     &                      dbl_mb(iheschi_ao),dbl_mb(id3chi_ao),
     &                      Amat,Cmat,Mmat,Lmat,
     &                      dbl_mb(i_ccdel),
     &                      ibf_ao, iniz, ifin, cetobfr,
     &                      ppp, ff, 
     &                      zzz, hh, iga_dens, force,
     &                      oactive)
         endif
        if (odftps) call pstat_off(ps_qv0b)
        if (grad) then
           if (.not.ma_pop_stack(l_ccdel))
     &        call errquit('dftggv0b: pop stack ',1,0)
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
     &        'Acof2', l_Acof2, k_Acof2)) call
     &        errquit('dftg_gridv0b: cannot allocate Acof2',1,
     &                MA_ERR)
         if (grad) then
            if (.not.MA_push_get(MT_DBL, ipol*9*nqpts, 'CCof2',
     &           l_Ccof2, k_Ccof2)) call
     &           errquit('dftg_gridv0b: cant alloc Ccof2',1,
     &                   MA_ERR)
         endif
         if (kske) then
            if (.not.MA_push_get(MT_DBL, ipol*3*nqpts, 'MCof2',
     &           l_Mcof2, k_Mcof2)) call
     &           errquit('dftg_gridv0b: cant alloc Mcof2',1,
     &                   MA_ERR)
         end if
c
         call xc_d2expl(tol_rho, 
     D        d2scr(3*nqpts*nbf_ao_mxnbf_ce+1),
     A        Amat, Amat2, dbl_mb(k_Acof2), 
     &        Cmat, Cmat2, dbl_mb(k_Ccof2),
     &        Mmat, Mmat2, dbl_mb(k_Mcof2),
     &        Fmat, Pmat, d2scr, 
     D        d2scr(3*nqpts*nbf_ao_mxnbf_ce+1),
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao),
     &        dbl_mb(iheschi_ao), dbl_mb(id3chi_ao),
     &        curatoms, ncuratoms, ipol, nqpts,
     &        nbf_ao, mbf_ao, grad, AO_bas_han, natoms, iniz,
     &        ifin, drhonuc, ddelrhonuc, dttau,
     &        delrho, iga_dens, hess, ibf_ao,
     &        rchi_atom, rdelchi_atom,
     &        rdens_atom, cetobfr, kske)
c
         if (grad) then
            if (.not.ma_pop_stack(l_Ccof2))
     &           call errquit('dftg_gridv0b: cannot pop CCof2',
     &           0, MA_ERR)
         endif
         if (kske) then
            if (.not.ma_pop_stack(l_Mcof2))
     &           call errquit('dftg_gridv0b: cannot pop CCof2',
     &           0, MA_ERR)
         endif
         if (.not.ma_pop_stack(l_Acof2))
     &        call errquit('dftg_gridv0b: cannot pop Acof2',0,
     &                     MA_ERR)
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
     &        Amat2, Cmat2, Cmat, Mmat2, delrho,
     &        3*ncuratoms, ipol, nqpts, grad, .false., ! kske = false 
     &        .false.)
c
         if (odftps) call pstat_on(ps_tabcd)
c Daniel (1-24-13): The number of perturbed densities is 3*ncuratoms, 
c which correspond to the x-, y-, and z-gradients for the active atoms
c in the system
c Daniel (2-7-13): This gives the matrix elements of the XC-kernel 
c multiplied by the ground state density gradient.  The required are 
c matrices are fed into xc_tabcd using drhonuc and ddelrhonuc.  In
c xc_tabcd, drhonuc -> Amat and ddelrhonuc -> Cmat (not to be confused
c with the XC-derivative matrices from xc_eval_fnl).
         call xc_tabcd(1,.true.,
     T        tol_rho, Fmat, Pmat, drhonuc, Bmat, ddelrhonuc,
     &        dbl_mb(iDmat), 
     N        0d0,0d0,0d0,0d0,
     C        Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     &        0d0,0,3*ncuratoms, 
     I        ipol,nqpts,mbf_ao,
     N        nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, 
     G        GRAD,AO_bas_han, natoms, iniz, ifin, 
     G        g_rhs, ibf_ao, rchi_atom, rdelchi_atom, 
     R        rdens_atom, cetobfr,
     &        .false., dttau, .false., 0d0, 0d0, 0d0)  ! kske = false
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
c Daniel (2-6-13): Here, we have the Hessian of the basis functions,
c and 1 XC-matrix (1 perturbation) rather than 3*Natoms as above.
         call xc_tabcd(2,.true.,
     &        tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, 
     D        dbl_mb(iDmat), 
     N        0d0,0d0,0d0,0d0,
     C        Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 
     H        dbl_mb(iheschi_ao),
     &        curatoms, ncuratoms, 1,
     I        ipol, nqpts, mbf_ao,
     &        nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, 
     G        GRAD, AO_bas_han, natoms, iniz, ifin, 
     G        g_rhs, ibf_ao, rchi_atom, rdelchi_atom, 
     R        rdens_atom, cetobfr,
     &        .false.,Mmat,.false.,0d0,0d0,0d0)  ! kske = false
c
         if (odftps) call pstat_off(ps_tabcd)
c Daniel (1-14-13): TDDFT XC-kernel gradient terms.
      elseif (do_gxc) then
c The formulae for doing XC-kernel gradients are shown in 
c J. Chem. Phys., 124, 094107.  
c 
c ----------------------------------------------
c Quadrature weight contribution to the gradient
c ----------------------------------------------
c This is done before the gradient of the XC-kernel part, because we 
c overwrite the matrix that stores the perturbed electron density at 
c that step.
        call tddft_grad_dxcdmat(func, rho, delrho, nqpts, ipol, 
     1    ngxc, grad, rho(iscf_rho), delrho(iscf_delrho),
     3    Amat, Cmat, Amat2, Cmat2, qwght, triplet, 2)
        natoms_l = natoms
c We calculate derivatives of quadrature weights by reusing the code
c from the ground state gradients.
        call dftg_gridpart(d_qwght, qwght, qxyz, xyz, xyzm, rq,
     1    ictr_buf, natoms_l, iandex, nqpts,
     2    whichssw, rad, expo, basis)
c Construct the quadrature weight derivatives multiplied by the 
c XC-kernel.  The call to the following routine is identical to the 
c call for ground state gradients, but func contains the XC-kernel
c contracted with the perturbed density.
        call dftg_wderiv(ictr_buf, natoms, nqpts, nqpts,
     1    iandex, natoms_l, func, d_qwght, force, oactive)

c Multiply the derivative arrays with the quadrature weights, since
c this wasn't done earlier to prevent double counting them.
c Daniel (3-11-13): The .false. passed to this subroutine is used to
c initialize kske.
        call setACmat_d3(delrho(iscf_delrho), Amat, Amat2, Amat3, 
     1         Cmat, Cmat2, Cmat3, qwght, ipol, nqpts, grad, 
     2         (.not. do_3rd), .false., .false.)
c 
c Based on the formulae in the article above, it is possible to 
c separate the problem into a contribution from the perturbed
c density and the ground state density.
c
c ----------------------------------------------
c Perturbed density contribution to the gradient
c ----------------------------------------------
        call tddft_grad_dfxc(rho, delrho, Amat3, Cmat3, Amat2, 
     1    Cmat2, Amat, Cmat, delrho(iscf_delrho), ngxc, 
     2    ipol, nqpts, grad, crho, cdelrho, 2, triplet)
        if (grad) then
          if (.not.ma_push_get(mt_dbl, nqpts*mbf_ao*ipol, 'ccdel',
     1      l_ccdel, i_ccdel))
     2      call errquit(
     3      'dftg_gridv0b: Failed to allocate ccdel', 0, MA_ERR)
        end if
        call tddft_grad_quadv0b(ictr_buf, ldew, ipol,
     1    grad, nbf_ao, mbf_ao, natoms,
     2    nqpts,
     3    dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
     4    dbl_mb(iheschi_ao), crho, cdelrho,
     5    dbl_mb(i_ccdel),
     6    ibf_ao, iniz, ifin, cetobfr,
     7    ppp, ff,
     8    zzz, hh, iga_dens, force,
     9    oactive, ngxc, 2)
c -------------------------------------------------
c Ground state density contribution to the gradient
c -------------------------------------------------
c Build the matrix for the gradient of the XC-kernel, for the
c products involving the gradient of the ground state density.
        call tddft_grad_dfxc(rho, delrho,
     1    Amat3, Cmat3, Amat2, Cmat2, Amat, Cmat,
     2    delrho(iscf_delrho), ngxc, ipol, nqpts,
     3    grad, crho, cdelrho, 1, triplet)
c Calculate the gradient due to the XC-functional term, involving
c the ground state density gradient.
        call tddft_grad_quadv0b(ictr_buf, ldew, ipol,
     1    grad, nbf_ao, mbf_ao, natoms,
     2    nqpts,
     3    dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
     4    dbl_mb(iheschi_ao), crho, cdelrho,
     5    dbl_mb(i_ccdel),
     6    ibf_ao, iniz, ifin, cetobfr,
     7    ppp, ff,
     8    zzz, hh, iga_dens, force,
     9    oactive, ngxc, 1)
        if (grad) then
          if (.not.ma_pop_stack(l_ccdel))
     1      call errquit(
     2      'dftg_gridv0b: Failed to deallocate ccdel', 0, MA_ERR)
        end if
c Daniel (2-7-13): TDDFT XC-potential gradient terms.
      else if (do_tddftvxc) then
cc The formulae for doing XC-potential gradients are shown in 
cc J. Chem. Phys., 124, 094107.  
cc 
cc ----------------------------------------------
cc Quadrature weight contribution to the gradient
cc ----------------------------------------------
cc This is done before the gradient of the XC-kernel part, because we 
cc overwrite the matrix that stores the perturbed electron density at 
cc that step.
c        call tddft_grad_dxcdmat(func, rho, delrho, nqpts, ipol, 
c     1    ngxc, grad, rho(iscf_rho), delrho(iscf_delrho),
c     3    Amat, Cmat, Amat2, Cmat2, qwght, triplet, 1)
c        natoms_l = natoms
cc We calculate derivatives of quadrature weights by reusing the code
cc from the ground state gradients.
c        call dftg_gridpart(d_qwght, qwght, qxyz, xyz, xyzm, rq,
c     1    ictr_buf, natoms_l, iandex, nqpts,
c     2    whichssw, rad, expo)
cc Construct the quadrature weight derivatives multiplied by the 
cc XC-kernel.  The call to the following routine is identical to the 
cc call for ground state gradients, but func contains the XC-potential
cc contracted with the perturbed density.
c        call dftg_wderiv(ictr_buf, natoms, nqpts, nqpts,
c     1    iandex, natoms_l, func, d_qwght, force, oactive)
c
cc Multiply the derivative arrays with the quadrature weights, since
cc this wasn't done earlier to prevent double counting them.
c        call setACmat_d2(delrho(iscf_delrho), Amat, Amat2, 
c     1         Cmat, Cmat2, qwght, ipol, nqpts, grad, 
c     2         (.not. do_2nd), kske, Mmat2)
cc 
cc Based on the formulae in the article above, it is possible to 
cc separate the problem into a contribution from the perturbed
cc density and the ground state density.
cc
cc ----------------------------------------------
cc Perturbed density contribution to the gradient
cc ----------------------------------------------
c        call tddft_grad_dvxc(rho, delrho, Amat2, Cmat2, Amat, Cmat,
c     2    delrho(iscf_delrho), ngxc, ipol, nqpts, grad,
c     3    crho, cdelrho, 2)
c        if (grad) then
c          if (.not.ma_push_get(mt_dbl, nqpts*mbf_ao*ipol, 'ccdel',
c     1      l_ccdel, i_ccdel))
c     2      call errquit(
c     3      'dftg_gridv0b: Failed to allocate ccdel', 0, MA_ERR)
c        end if
c        call tddft_grad_quadv0b(ictr_buf, ldew, ipol,
c     1    grad, nbf_ao, mbf_ao, natoms,
c     2    nqpts,
c     3    dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
c     4    dbl_mb(iheschi_ao), crho, cdelrho,
c     5    dbl_mb(i_ccdel),
c     6    ibf_ao, iniz, ifin, cetobfr,
c     7    ppp, ff,
c     8    zzz, hh, iga_dens, force,
c     9    oactive, ngxc, 2)
cc -------------------------------------------------
cc Ground state density contribution to the gradient
cc -------------------------------------------------
cc Build the matrix for the gradient of the XC-potential, for the
cc products involving the gradient of the ground state density.
c        call tddft_grad_dvxc(rho, delrho, Amat2, Cmat2, Amat, Cmat,
c     2    delrho(iscf_delrho), ngxc, ipol, nqpts, grad,
c     3    crho, cdelrho, 1)
cc Calculate the gradient due to the XC-functional term, involving
cc the ground state density gradient.
c        call tddft_grad_quadv0b(ictr_buf, ldew, ipol,
c     1    grad, nbf_ao, mbf_ao, natoms,
c     2    nqpts,
c     3    dbl_mb(ichi_ao),dbl_mb(idelchi_ao),
c     4    dbl_mb(iheschi_ao), crho, cdelrho,
c     5    dbl_mb(i_ccdel),
c     6    ibf_ao, iniz, ifin, cetobfr,
c     7    ppp, ff,
c     8    zzz, hh, iga_dens, force,
c     9    oactive, ngxc, 1)
c        if (grad) then
c          if (.not.ma_pop_stack(l_ccdel))
c     1      call errquit(
c     2      'dftg_gridv0b: Failed to deallocate ccdel', 0, MA_ERR)
c        end if
      endif
c
      if (.not.ma_chop_stack(lchi_ao))
     &   call errquit('dftg_gridv0b: cannot chop stack',0,
     &     MA_ERR)
c
c
      if (do_force.and.ldew)then
         natoms_l=natoms
         if (odftps) call pstat_on(ps_ssw)                  
         call dftg_gridpart(d_qwght,qwght, qxyz, xyz,xyzm, 
     .        rq, ictr_buf, natoms_l, iandex,
     .        nqpts,whichssw,
     ,        rad,expo,basis)
         call dftg_wderiv(ictr_buf, natoms, nqpts, nqpts, 
     .        iandex, 
     &        natoms_l, func, d_qwght,
     &        force, oactive)
         if (odftps) call pstat_off(ps_ssw)
      endif
c         
      return 
      end
