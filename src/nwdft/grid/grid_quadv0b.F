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
     &     nmat, do_scf, do_cpks_l, do_cpks_r,do_nmrks_r,tmat,
     &     l3d,triplet,
     &     bas_cent_xinfo,bas_cset_xinfo,
     &     docsetx, 
     &     tau, ttau, kske, Mmat, Mmat2,
     &     StericEnergy)
c
c$Id: grid_quadv0b.F,v 1.64 2008-02-01 22:27:45 niri Exp $
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
#include "zora.fh"
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
      integer ldmat,idmat,i, grid_nbfm
      double precision dabsmax,wmax
      external dabsmax
      external grid_nbfm
      integer iscf_rho,iscf_delrho
      integer iscf_tau,iscf_ttau
      logical grid_sicinit,grid_sicend
      external grid_sicinit,grid_sicend
      double precision acc_xc_gauss
      logical l3d
      logical triplet, kske
      integer lemat,iemat,lfmat,ifmat,k_scr,l_scr
      logical do_2nd
      double precision StericEnergy
c
      npol = (ipol*(ipol+1))/2
      acc_AO_gauss = dble(iAOacc)
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if(nqpts.eq.0) return
      if (grad) then
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
         call icopy(mbf_ao, 0,0, ibf_ao, 1)
         mbf_ao=grid_nbfm(AO_bas_han,  natoms, 
     &        ictr_buf,rad,xyz,qxyz,nqpts,
     &        ibf_ao, docset, iniz, ifin, expo, 
     &        minexp,ldiff,acc_ao_gauss,iatype_pt_chg)
         if (mbf_ao.eq.0) return
      else
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in grid_quadv0b',0, BASIS_ERR)
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
     &     call errquit('grid_quadv0b: cannot allocate chi_ao',0,
     &       MA_ERR)
      if( grad ) then
         if (.not.MA_Push_Get(mt_dbl, 3*nqpts*mbf_ao,
     &      'delchi_ao', ldelchi_ao, idelchi_ao)) then
           call errquit('gridquadv0b: cant alloc delchi_ao',0,0)
        endif
      endif
      if (sic_orb_index.eq.1) then
         if(.not.grid_sicinit(nqpts,
     &        l_vect1,i_vect1, ldelrho_sig,idelrho_sig))
     &        call errquit(' grid_quadv0b: grid_sicinit failed',0,
     &       CALC_ERR)
      endif
c     
      if (odftps) call pstat_on(ps_ebf)
      call qdist(rchi_atom, rq, qxyz, xyz, nqpts, natoms)
      call xc_eval_basis(ao_bas_han, maxdbas, dbl_mb(ichi_ao),
     &     dbl_mb(idelchi_ao), 0d0, 0d0, rq, qxyz, xyz, nqpts, natoms,
     &     iniz, docset, bas_cent_info, bas_cset_info)
      if (odftps) call pstat_off(ps_ebf)
c
c     get reduced Xi(r) and dXi(r) over atoms
c
      call util_rmsatbf(nqpts, natoms,iniz,ifin,
     ,     dbl_mb(ichi_ao),rchi_atom)
      if( grad ) 
     .     call util_rmsatbf(3*nqpts, natoms,iniz,ifin,
     ,     dbl_mb(idelchi_ao),rdelchi_atom)
      call grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
     ,     rchi_atom,rdelchi_atom,rdens_atom,iniz,
     ,     dbl_mb(ichi_ao),dbl_mb(idelchi_ao),mbf_ao,grad,
     ,     ipol,iscratch)
         if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &        'Dmat', lDmat,iDmat))
     &        call errquit('grid_quadv0b: cannot allocate Dmat',0,
     &       MA_ERR)
      if(nqpts.eq.0) goto 2010

      wmax= dabsmax(nqpts,qwght)
      if(do_scf.or.do_nmrks_r) then
         if (odftps) call pstat_on(ps_xcrho)
            call xc_rhogen(0,
     T        tol_rho,AO_bas_han,g_dens,nbf_ao_mxnbf_ce, 
     N        natoms, 1,1,1,
     I        ipol, nqpts, nbf_ao, mbf_ao, grad, (ipol*(ipol+1))/2,
     F        Fmat, Pmat, Bmat, dbl_mb(idmat),
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 0d0,
     &        ibf_ao, iniz, ifin, rho, delrho, 
     &        rchi_atom, rdelchi_atom,
     &        rdens_atom, cetobfr,wmax, 
     &        tau, ttau, kske)
         if (odftps) call pstat_off(ps_xcrho)
         if (sic_orb_index.eq.1) then
            call grid_sicrho(nqpts,mbf_ao,cetobfr,wmax,
     ,        rho,delrho,grad,ibf_ao, iniz, ifin,
     ,        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 
     ,        rchi_atom, rdelchi_atom,
     ,        rdens_atom,
     ,        Fmat, Pmat, Bmat,
     ,        dbl_mb(i_vect1),dbl_mb(idelrho_sig),
     &        tau, ttau, kske)
         endif
      else  if (do_cpks_l) then
            call xc_rhogen(1,
     T        tol_rho, AO_bas_han, g_dens, nbf_ao_mxnbf_ce, 
     N        natoms, 1,1, nmat+1,
     I        ipol, nqpts, nbf_ao, mbf_ao, grad, ipol,
     F        Fmat, Pmat, Bmat, dbl_mb(idmat),
     &        dbl_mb(ichi_ao), dbl_mb(idelchi_ao), 0d0,
     &        ibf_ao, iniz, ifin, rho, delrho, 
     &        rchi_atom, rdelchi_atom, 
     &        rdens_atom, cetobfr,wmax, 
     &        tau, ttau, kske)
      else if (do_cpks_r) then
       call errquit('gridquadv0b: what do do for cpks rhs?',0, CALC_ERR)
      endif
c     
c     Evaluate the integrated charge density.
c     

      nqtot = nqtot + nqpts
c     
      rho0=rho_int
      rho_int = rho_int + ddot(nqpts, qwght, 1, rho, 1)

      if(do_scf.and.abs(rho0-rho_int).lt.tol_rho) goto 2010
c     
c
c     Evaluate the XC potential and energy sampling matrices.
c     Set up pointer to the SCF density for the CPS LHS case.
c
      iscf_rho = 1
      iscf_delrho = 1
      if (do_cpks_l) then
          iscf_rho = 1 + nqpts*ipol*nmat
          iscf_delrho = 1 + nqpts*3*ipol*nmat
      endif
c
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
c
czora...start
czora...calculate the zora correction (hartree and nuclear potentials included in correction)
      if (do_scf .and. do_zora_get_correction) then 

czora...calculate the zora correction (consolidated density => hartree piece)
        call grid_zoracorr(nqpts,qxyz,qwght,natoms,g_dens,amat)

czora...evaluate zora contributions
           call zora_tabcd(0,l3d,
     T        tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, Dbl_MB(iDmat), 
     &        0d0,0d0,0d0,0d0,
     &        Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     N        0d0,0,1,
     I        ipol,nqpts,mbf_ao,nbf_ao_mxnbf_ce,nbf_ao_mxnbf_ce2, 
     G        grad, AO_bas_han,natoms, iniz, ifin, 
     G        g_vxc, ibf_ao, rchi_atom, rdelchi_atom,
     &        rdens_atom, cetobfr,
     &        .false.,0d0,0d0)

      else  ! do_zora_get_correction
c
      if (odftps) call pstat_on(ps_excf)
      do_2nd = (do_cpks_l .or. do_cpks_r)
      call xc_eval_fnl(rho(iscf_rho), delrho(iscf_delrho), 
     &     Amat, Amat2, Cmat, Cmat2, nqpts, Ex, Ec, qwght,
     &     grad, .false., 0d0, do_2nd, 
     &     tau, ttau, kske, Mmat, Mmat2,
     &     StericEnergy)
      if (odftps) call pstat_off(ps_excf)
c
      if (sic_orb_index.eq.1) then
         call grid_sicpot(nqpts,qxyz,qwght,rho,delrho,Amat,cmat,
     ,     dbl_mb(idelrho_sig),dbl_mb(i_vect1))
      endif
c
      if(xcfit) then
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
            call xc_tabcd(0,l3d,
     &        tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, Dbl_MB(iDmat), 
     &        0d0,0d0,0d0,0d0,
     &        Dbl_MB(ichi_ao), dbl_MB(idelchi_ao), 0d0,
     &        0d0,0,1,
     &        ipol, nqpts, mbf_ao, nbf_ao_mxnbf_ce, nbf_ao_mxnbf_ce2, 
     &        grad, AO_bas_han, 
     &        natoms, iniz, ifin, g_vxc, ibf_ao,
     &        rchi_atom, rdelchi_atom,
     &        rdens_atom, cetobfr,
     &        kske,Mmat,0d0)
         if (odftps) call pstat_off(ps_tabcd)
c
         else if (do_nmrks_r) then
            if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts,
     &           'iEmat', lEmat,iEmat))
     &           call errquit('grid_quadv0b: cannot allocate iEmat',0,
     &          MA_ERR)
            if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*nqpts*3,
     &           'iFmat', lFmat,iFmat))
     &           call errquit('grid_quadv0b: cannot allocate iFmat',0,
     &          MA_ERR)
            if (.not.MA_Push_Get(mt_dbl,nqpts,
     &           'iFmat', l_scr,k_scr))
     &           call errquit('grid_quadv0b: cannot allocate iFmat',0,
     &          MA_ERR)
c
            call xc_tabcd(3,.true.,
     T           tol_rho, Fmat, Pmat, Amat, Bmat, Cmat, Dbl_MB(iDmat), 
     E           Dbl_MB(iEmat), Dbl_MB(iFmat),  qxyz,xyz,
     &           Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     N           0d0,0,1,
     I           ipol,nqpts,mbf_ao,nbf_ao_mxnbf_ce,nbf_ao_mxnbf_ce2, 
     G           grad, AO_bas_han,natoms, iniz, ifin, 
     G           g_vxc, ibf_ao, rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr,
     &           kske, Mmat, dbl_mb(k_scr))
c
            if (.not.ma_chop_stack(lEmat))
     &         call errquit('grid_quadv0b: cannot pop stack', 1, MA_ERR)
         else if (do_cpks_l) then
c
c     Combine derivs of functional with deriv densities to get coeffs
c     of bfn pairs in CPKS LHS --> these are written back into rho, delrho, ttau
c
            call xc_cpks_coeff(rho, delrho, ttau, 
     &           Amat2, Cmat2, Cmat, Mmat2,
     &           delrho(iscf_delrho), nmat, ipol, nqpts, grad, kske, 
     &           triplet)
c
            call xc_tabcd(1,l3d,
     T           tol_rho,Fmat,Pmat, rho, Bmat,delrho, Dbl_MB(iDmat), 
     N           0d0,0d0,0d0,0d0,
     C           Dbl_MB(ichi_ao), Dbl_MB(idelchi_ao), 0d0,
     N           0d0,0,nmat, 
     I           ipol,nqpts,mbf_ao,nbf_ao_mxnbf_ce,nbf_ao_mxnbf_ce2, 
     G           grad, AO_bas_han, natoms, iniz, ifin, 
     G           g_vxc, ibf_ao, rchi_atom, rdelchi_atom,
     &           rdens_atom, cetobfr,
     &           kske, ttau, 0d0)
c
         else if (do_cpks_r) then
            call errquit('gridquadv0b: what do do for cpks rhs?',0,
     &       CALC_ERR)
         endif
c
      end if ! do_zora_get_correction
      endif

 2010 continue
c
      if (.not.ma_pop_stack(lDmat))
     &      call errquit('grid_quadv0b: cannot pop stack', 1, MA_ERR)
      if (sic_orb_index.eq.1) then
         if(.not.grid_sicend(l_vect1,ldelrho_sig))
     ,        call errquit(' grid_quadv0b: grid_sicend failed',0,
     &       CALC_ERR)
      endif
      if(grad)then
         if (.not.ma_pop_stack(ldelchi_ao))
     &        call errquit('grid_quadv0b: cannot pop stack', 2, MA_ERR)
      endif
      if (.not.ma_pop_stack(lchi_ao))
     &     call errquit('grid_quadv0b: cannot pop stack', 3, MA_ERR)
c
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
C$Id: grid_quadv0b.F,v 1.64 2008-02-01 22:27:45 niri Exp $
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
      logical veryfar
      gaussian_range(n,logeps,alpha) = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      call ifill(mcenters,0,iniz,1)
      call ifill(mcenters,0,ifin,1)
      if(acc_gauss.gt.25d0.or.acc_gauss.lt.3d0) call errquit(
     ' ' grid_nbfm: silly accgauss ',nint(acc_gauss), UNKNOWN_ERR)
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
c     check on most diffuse fn
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
c       get info about current contraction set
c       
          if (.not. bas_continfo(basis_hand, icset,  l ,nprim,
     &         ncontr, isphere))
     &         call errquit('Exiting in grid_nbfm.',5, BASIS_ERR)
c         
c       
c       get exponents and contraction coefficients for this contraction set
c       
        if (.not.bas_get_exponent(basis_hand, icset, zprim))
     &       call errquit('Exiting in grid_nbfm.',7, BASIS_ERR)
c       
c       Determine the minimum Gaussian exponent.
c       

        zmin = 1.D+06
        do iprimo = 1,nprim
          zmin = min(zprim(iprimo),zmin)
        enddo
c       
c       Only include those basis functions that are "non-zero" for at least 
c       one  point in the sampling set.
c       
        nshbf=ncontr*(((l+1)*(l+2))/2)
        if(isphere.eq.1) then
           nshbf=ncontr*(2*l+1)
        endif
c
c     pre-screening for worst case (max radius)
c
              if (r_q0.lt.gaussian_range(l,acc_loc,zmin)) then
                 do m=1,nshbf
                    ibf(n1+m ) = npt+m
                 enddo
                 docset(icset) = .true.
                 if (iniz(ictr).eq.0) iniz(ictr)=n1+1
                 n1=n1+nshbf
              endif

        npt=npt+nshbf
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
      Subroutine xc_make_basis_info(basis_hand, bas_cent_info,
     &     bas_cset_info, mcenters)
c
C$Id: grid_quadv0b.F,v 1.64 2008-02-01 22:27:45 niri Exp $
c
      implicit none
#include "errquit.fh"
c      
      integer basis_hand        ! [input]
      integer mcenters          ! [input]
      integer bas_cent_info(3,mcenters) ! [output]
      integer bas_cset_info(6,*) ! [output]
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
               call errquit('L code < 0 in xc_make_basis_info',1,
     &       BASIS_ERR)
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
