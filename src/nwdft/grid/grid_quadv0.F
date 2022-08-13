      subroutine grid_quadv0(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,tmat)
c
C$Id$
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(*)
      double precision rho_n
      double precision Exc(2)
      double precision tmat(*)  ! XC fit
cc AJL/Begin/FDE
      integer g_dens_fde_dummy(2) ! Acts as a dummy variable

      call grid_quadv0_gen(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,
     ,     1,1,tmat,.false.,.false.,g_dens_fde_dummy)
cc AJL/End

      return
      end
c
      subroutine grid_quadv0_gen(rtdb, g_dens, g_vxc,
     &      nexc, rho_n, Exc, nmat, calc_type ,tmat,l3d_dum,triplet,
cc AJL/Begin/FDE
     &      g_dens_fde)
cc AJL/End

      implicit none

#include "errquit.fh"
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(4)

      Logical oprint,grad,xc_chkgrad,oprintw
      Logical xc_chktau, kske
      logical xc_chklap, kslap
      logical xc_chkdelrq, ksdelrq
      external xc_chkgrad,xc_rep_close, xc_chktau, xc_chklap,xc_chkdelrq
      logical xc_rep_close
      double precision Exc(2)
      double precision rho_n
      integer nmat              ! [input] number of XC matrices to make;
c                               !         only used for lhs cpks
      double precision tmat(*)  ! XC fit
      integer calc_type         ! what to calculate [input]
c                               !   1 = XC matrices for SCF
c                               !   2 = XC matrices for CPKS LHS
c                               !   3 = XC matrices for CPKS RHS
c                               !   4 = XC matrices for GIAO NMR RHS
c                               !   5 = XC matrices for third derivatives
c                               !   6 = post-SCF XDM energy

      logical l3d_dum           ! true if 3d GA is used
      logical triplet           ! true if triplet excitation energy calc
      double precision StericEnergyTotal
      logical l3d               ! true if 3d GA is used
c
#include "cgridfile.fh"
#include "geom.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "dft2drv.fh"
#include "steric.fh"
#include "dft3drv.fh"
c
      logical rtdb_status       
      logical grid_reopen
      external grid_reopen
      logical lsa ! Qin
c
      integer me, nproc, nxyz, nctrs, 
     &        nbf_ao_mxnbf_ce2, nqtot, 
     &        ncoef_max,
     &        npol
c
      integer lxyz, ixyz, lcharge, icharge, 
     &        lqxyz, iqxyz, lqwght, iqwght, 
     &        lniz, iniz, lexpo, iexpo, lcoef, icoef, lFmat, iFmat, 
     &        lPmat, iPmat, lrho, irho, ldelrho, idelrho, lBmat, iBmat,
     &        lhess, ihess, lrq, irq, 
     &        lfin, ifin, ltags, itags,
     &        lmmat, immat, lttau, ittau, lLmat, iLmat, iz, lz
c
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr,lamat,iamat,icmat,lcmat

      double precision nelec
c
      integer libf_ao,iibf_ao,irqbuf,lrqbuf,lbp,ibp,
     .     irchi_atom,irdelchi_atom,lrchi_atom,lrdelchi_atom
      integer iqs,lqs,nqlen,mxnrad,
     , laqbuf, iaqbuf,lRij,iRij,ndbl,l_qat,k_qat,irqlen
      integer l_Amat2, k_Amat2, l_Cmat2, k_Cmat2, l_Mmat2, k_Mmat2
      integer len_rho, len_delrho, len_Amat2, len_Cmat2, len_Mmat2
      integer len_Amat
      integer len_mmat, len_ttau
      integer len_Lmat
      logical do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm
      integer max_elem_ang_scr, max_component
      integer nexpo,nibf
      double precision rho_int, ex, ec,dummy(3),job_grid_acc,rho_err
      integer ii,isp
      equivalence (dummy(1),Ex), (dummy(2),Ec), (dummy(3),rho_int)
c
      logical use_nwxc
      logical fon ! fractional occupancy
      double precision pstrace  ! a measure of the total charge
c
      logical do_gxc
      integer l_Amat3, k_Amat3, l_Cmat3, k_Cmat3
      integer len_Amat3, len_Cmat3
      logical lcpksl_gxc

C AJL/BEGIN/FDE
c It seems a bit excessive to have to declare all these extra variables
c and pass them through the subroutines, but without some serious
c re-writing of quadv0b there is no real alternative, as we need to
c generate the density of both the QM and FDE systems where applicable
      integer ncenters_fde
      integer icetobfr_fde, lcetobfr_fde
      integer irdens_atom_fde, lrdens_atom_fde
      integer lFmat_fde, iFmat_fde
      integer lPmat_fde, iPmat_fde
      integer lBmat_fde, iBmat_fde
      integer lAmat_fde, iAmat_fde
      integer lAmat2_fde, iAmat2_fde
      integer lAmat3_fde, iAmat3_fde
      integer lCmat_fde, iCmat_fde
      integer lCmat2_fde, iCmat2_fde
      integer lCmat3_fde, iCmat3_fde
      integer lMmat_fde, iMmat_fde
      integer lMmat2_fde, iMmat2_fde
      integer lttau_fde, ittau_fde

      integer lrho_fde, irho_fde
      integer ldelrho_fde, idelrho_fde

      integer lrho_qmfde, irho_qmfde
      integer ldelrho_qmfde, idelrho_qmfde
      integer lttau_qmfde, ittau_qmfde

c Copied out of grid_quadv0b
      integer lniz_fde, iniz_fde
      integer ixyz_fde, lxyz_fde
      integer icharge_fde, lcharge_fde
      integer itags_fde, ltags_fde
      integer l_rchi_atom_fde, i_rchi_atom_fde
      integer lrdelchi_atom_fde, irdelchi_atom_fde
      integer l_rq_fde, i_rq_fde
      integer l_ibf_fde, i_ibf_fde
      integer lfin_fde, ifin_fde
c
      integer g_dens_fde(2)
      integer nbf_fde_mxnbf_ce2
      character*32 pname
      integer g_fde(4)
      integer i
      double precision rho_int_fde
      integer fde_option
      integer g_s

      integer nshells_bas,nbf_mxnbf_ce2
      integer basis

#ifdef SOLARIS
c     needed for -stackvar option
      save ex, dummy, ec, rho_int, rho_int_fde
#endif
c    
      pname = 'grid_quadv0: '
cc AJL/End
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      oprint = util_print('quadrature', print_high)
      oprintw = util_print('quadrature warning', print_medium)
      grad = xc_chkgrad()
      kske = xc_chktau()
      kslap = xc_chklap()
      ksdelrq = xc_chkdelrq()
      l3d = (ga_ndim(g_vxc).eq.3)

      if (adft) then
        basis = xc_bas_han
      else
        basis = ao_bas_han
      endif
c
      irdelchi_atom = 0
c
      npol = (ipol*(ipol+1))/2
      if (.not.rtdb_get(rtdb,'dft:largenode', mt_log, 1, largenode))
     &     largenode=.false.
c
c     Open grid pts file
c     
cng
      if (grid_written .and. n_rec_in_file .le. 0 ) then
c
c       Somebody (fock_xc?) closed the grid file, better reopen it
c       before continuing. Note, fock_xc tends to open the grid file
c       at the beginning and close it at the end.
c
        if (.not.grid_reopen(geom)) call errquit(
     &    pname//'could not open grid file',0,DISK_ERR)
      endif
cng
      if (.not. grid_written .or. n_per_rec .le. 0 ) then
         call grid_packing_info
      end if
c
      if (oprint) then      
        write(*,*) " *** "//pname//"calc_type:",calc_type 
        write(*,*) " *** "//pname//"grid_written:",grid_written 
        write(*,*) " *** "//pname//"grad:",grad 
        write(*,*) " *** "//pname//"kske:",kske 
        write(*,*) " *** "//pname//"kslap:",kslap
      end if
c
c     Set up logicals denoting calculation type
c
      do_scf = calc_type .eq. 1
      do_cpks_l = calc_type .eq. 2
      do_cpks_r = calc_type .eq. 3
      do_nmrks_r = calc_type .eq. 4
      do_gxc = calc_type .eq. 5 ! for 3rd derivatives
      do_xdm = calc_type .eq. 6
C AJL/Begin/FDE
c Defined options for FDE:
c fde_option: 1 = ESP
c             2 = QM with Exc FDE
c             3 = FDE with Exc FDE
c             4 = QM/MM with Exc FDE
c             5 = QM with Ts FDE
c             6 = FDE with Ts FDE
c             7 = QM/MM with Ts TDE
      if (calc_type .gt. 6) then
        fde_option = calc_type - 6
        if (fde_option .gt. 1) do_scf = .true.
      else 
        fde_option = 0
      endif
C AJL/End

      if (.not. (do_scf .or. do_cpks_l .or. do_cpks_r .or. 
     &    do_nmrks_r .or. do_gxc. or. do_xdm .or.
     &    (fde_option .eq. 1))) then
         call errquit(pname//'illegal calculation type',1,
     &       INPUT_ERR)
      endif

      if (kslap.and.((do_cpks_l .or. do_cpks_r .or.do_nmrks_r .or. 
     &                do_gxc. or. do_xdm .or. xcfit .or.
     &                (fde_option.ne.0)))) then
         call errquit(pname//'illegal calculation type',1,
     &       INPUT_ERR)
      endif
c
      if (oprint) then
        write(*,*) " *** "//pname//"do_scf:",do_scf 
        write(*,*) " *** "//pname//"fde_option:",fde_option
        write(*,*) " *** "//pname//"n_per_rec:",n_per_rec
      end if
c
cc AJL/Begin/FDE
c These need defining else it causes confusion when calculating the ESP
c for frozen density embedding (fde_option.eq.1)
      len_rho = 1
      len_delrho = 1
cc AJL/End
      len_mmat = 1
      len_ttau = 1
      len_Lmat = 1
      len_Amat2 = 1
      len_Cmat2 = 1
      len_Mmat2 = 1
      len_Amat3 = 1
      len_Cmat3 = 1 
      use_nwxc = util_module_avail("nwxc")
      if (use_nwxc) then 
         call nwxc_getvals("nwxc_is_on",use_nwxc)
      endif
      if (use_nwxc) then
c       the Automatic Differentiation code uses the unrestricted open-
c       shell approach for 2nd and higher order derivatives.
        len_Amat  = 2*n_per_rec
      else
        len_Amat  = ipol*n_per_rec
      endif
c
      if (do_scf .or. do_nmrks_r. or. do_xdm) then
         len_rho = npol*n_per_rec
         len_delrho = 3*ipol*n_per_rec
         if (kske.or.ksdelrq) len_ttau = npol*n_per_rec
         if (kske) len_mmat =ipol*n_per_rec  
         if (kslap) len_Lmat = ipol*n_per_rec
      elseif (do_cpks_l) then
         len_rho = ipol*n_per_rec*(nmat+1)+n_per_rec
         len_delrho = 3*ipol*n_per_rec*(nmat+1)+3*n_per_rec
         if (kske) len_ttau = ipol*n_per_rec*(nmat+1)+n_per_rec
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
         if (kske) len_Mmat2 = n_per_rec*NCOL_MMAT2
      elseif (do_cpks_r) then
         call errquit(pname//'what for do_cpks_r?',0,
     &       UNKNOWN_ERR)
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
c
      elseif (do_gxc) then
c Here, len_rho and len_delrho reflect the fact that we have ipol ground
c state density matrices and ipol*nmat perturbed density matrices, where
c nmat depends on the number of roots requested.  These should be the
c correct dimensions.  The extra n_per_rec or 3*n_per_rec accounts for
c the total ground state density (or gradient) for unrestricted 
c calculations.
         len_rho = ipol*n_per_rec*(nmat+1) + n_per_rec
         len_delrho = 3*ipol*n_per_rec*(nmat+1) + 3*n_per_rec
         len_Amat2 = n_per_rec*NCOL_AMAT2
         len_Cmat2 = n_per_rec*NCOL_CMAT2
         len_Amat3 = n_per_rec*NCOL_AMAT3
         len_Cmat3 = n_per_rec*NCOL_CMAT3
      endif
      if (.not. geom_ncent(geom, nctrs)) 
     &   call errquit(pname//'geom_ncent failed',geom, GEOM_ERR)
      nxyz = 3*nctrs
      if (.not.MA_Alloc_Get(MT_Dbl,nxyz,'xyz',lxyz,ixyz))
     &   call errquit(pname//'cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_int,nctrs,'atns',lz,iz))
     &     call errquit(pname//'cannot allocate atns',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'charge',lcharge,icharge))
     &   call errquit(pname//'cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,nctrs*16,'tags',ltags,itags))
     &   call errquit(pname//'cannot allocate tags',0, MA_ERR)

      if (.not. geom_cart_get2(geom, nctrs, 
     &     Byte_MB(itags), Dbl_MB(ixyz), Dbl_MB(icharge),int_mb(iz)))
     &   call errquit(pname//'geom_cart_get2 failed',74, GEOM_ERR)


      if (.not. geom_cart_get(geom, nctrs, Byte_MB(itags),
     &                        Dbl_MB(ixyz), Dbl_MB(icharge)))
     &   call errquit(pname//'geom_cart_get failed',74, GEOM_ERR)
c
c      write(*,*)' Coordinates from grid_den.F '
c      call output(dbl_mb(ixyz), 1, 3, 1, nctrs, 3, nctrs, 1)
c     &                    
c
      if (.not.ma_pop_stack(ltags))
     &   call errquit(pname//'cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit(pname//'cannot pop stack',0, MA_ERR)
c
      nbf_ao_mxnbf_ce2 = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (do_nmrks_r) nbf_ao_mxnbf_ce2=nbf_ao_mxnbf_ce2*3
c AJL/Begin/FDE
      nbf_fde_mxnbf_ce2 = nbf_fde_mxnbf_ce*nbf_fde_mxnbf_ce
      if (do_nmrks_r) nbf_fde_mxnbf_ce2=nbf_fde_mxnbf_ce2*3
c AJL/End
c
c     Divide the multi-center integrand into a set of single center integrands.
c
      nqtot = 0
      rho_int = 0.D0
      rho_int_fde = 0.D0
c
      if (.not.MA_Push_get(mt_dbl,3*n_per_rec,'qxyz',lqxyz,iqxyz))
     &   call errquit(pname//'cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,n_per_rec,'qwght',lqwght,iqwght))
     &   call errquit(pname//'cannot allocate qwght',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom list',lniz,iniz))
     &   call errquit(pname//'cannot allocate atom list',0, MA_ERR)
      if (.not.MA_Push_get(MT_int,nctrs,'atom nbf',lfin,ifin))
     &   call errquit(pname//'cannot allocate atom nbf',0, MA_ERR)
cc AJL/Begin/FDE
c Check for largest value of nbf
      nexpo=nbf_ao_mxprim
      if(xcfit)nexpo=max(nexpo,nbf_xc_mxprim)
      if(fde_option.gt.0)nexpo=max(nexpo,nbf_fde_mxprim)
c      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
      if (.not.MA_Push_Get(MT_Dbl,nexpo,'expo',lexpo,iexpo))
     &   call errquit(pname//'cannot allocate expo',0, MA_ERR)
cc AJL/End
c
      ncoef_max = nbf_ao_mxnbf_cn * nbf_ao_mxprim
      if(xcfit) ncoef_max=max(ncoef_max,nbf_xc_mxnbf_cn * nbf_xc_mxprim)
c
      if (.not.MA_Push_Get(MT_Dbl,ncoef_max,'coef',lcoef,icoef))
     &   call errquit(pname//'cannot allocate coef',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENSF',lFmat,iFmat))
     &   call errquit(pname//'cannot allocate DENSF',0, MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,nbf_ao_mxnbf_ce2,'DENS',lPmat,iPmat))
     &   call errquit(pname//'cannot allocate DENS',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho',lrho,irho))
     &   call errquit(pname//'cannot allocate rho',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,len_delrho,'delrho',ldelrho, idelrho))
     &   call errquit(pname//'cannot allocate delrho',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl,len_ttau,'ttau',lttau, ittau))
     &   call errquit(pname//'cannot allocate ttau',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,nbf_ao_mxnbf_ce*n_per_rec,'Bmat',
     &   lBmat,iBmat))call errquit(pname//'cannot allocate Bmat',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,1,'hesdum',lhess,ihess))
     &   call errquit(pname//'cannot allocate hesdum',0, MA_ERR)
c
c     Assume this numerical density routine can be called from anywhere.

c     -compute mapping vectors
c
      nshells_bas = nshells_ao
      if (adft) nshells_bas = nshells_xc

      if (.not.ma_push_get
     &   (mt_int,nctrs*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit(pname//'push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_bas,'cntoce map',lcntoce,icntoce))
     &   call errquit(pname//'push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_bas*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit(pname//'push_get failed', 13, MA_ERR)
c     
      call build_maps(basis, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), nctrs, nshells_bas)
      if (.not.ma_chop_stack(lcntoce))
     &   call errquit(pname//'cannot pop stack',1, MA_ERR)
c
c     -compute reduced total density matrices over atoms
c
      if (.not.MA_Push_Get(MT_Dbl,ipol*nctrs*nctrs,'rdens_atom',
     &   lrdens_atom,irdens_atom))
     &   call errquit(pname//'cannot allocate rdens_atom',0, MA_ERR)
c      
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce*(nmat+1)
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit(pname//'cannot allocate scr',0, MA_ERR)
c Daniel (1-11-13): I think we can stay similar to calc_type = 2 for
c XC-third derivatives in this portion of the code, since this part
c doesn't depend on having 3-dimensional GAs.
      if (calc_type.eq.2.or.calc_type.eq.5) then
         do isp=1,ipol
            call util_ga_mat_reduce(nbf_ao, nctrs, int_mb(icetobfr), 
     A           g_dens(1+(isp-1)*(nmat+1)),  1, 
     A           dbl_mb(irdens_atom+(isp-1)*nctrs*nctrs), 'absmax', 
     &           dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
            do ii=2,(nmat+1)
               call util_ga_mat_reduce(nbf_ao,nctrs,int_mb(icetobfr), 
     A              g_dens(ii+(isp-1)*(nmat+1)),  1, 
     A              dbl_mb(irdens_atom+(isp-1)*nctrs*nctrs), 'absmax', 
     &              dbl_mb(iscr), nbf_ao_mxnbf_ce,.false.)
            enddo
         enddo
      elseif (adft) then
        call adft_reduce(nbf_xc, nctrs, int_mb(icetobfr), tmat, ipol,
     &                   dbl_mb(irdens_atom), nbf_xc_mxnbf_ce)
      else   
      call util_ga_mat_reduce(nbf_ao, nctrs, int_mb(icetobfr), g_dens, 
     &                        ipol, dbl_mb(irdens_atom), 'absmax', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
      endif
c
      if (.not.ma_pop_stack(lscr))
     &     call errquit(pname//'cannot pop stack',0, MA_ERR)

c
c
cc AJL/BEGIN/FDE
c Duplicates of the above allocations but for the FDE environment
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then
c
        if (.not.MA_Push_Get(MT_dbl,nbf_fde_mxnbf_ce2,'DENSF_fde',
     &     lFmat_fde,iFmat_fde))
     &     call errquit(pname//'cannot allocate DENSF_fde',0, MA_ERR)
        if (.not.MA_Push_Get(MT_dbl,nbf_fde_mxnbf_ce2,'DENS_fde',
     &     lPmat_fde,iPmat_fde))
     &     call errquit(pname//'cannot allocate DENS_fde',0, MA_ERR)
c
        if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho_fde',lrho_fde,
     &     irho_fde))
     &     call errquit(pname//'cannot allocate rho_fde',0, MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,len_delrho,'delrho_fde',
     &     ldelrho_fde, idelrho_fde))
     &     call errquit(pname//'cannot allocate delrho_fde',0, MA_ERR)
c
        if (.not.MA_Push_Get(MT_Dbl,len_ttau,'ttau_fde',
     &     lttau_fde, ittau_fde))
     &     call errquit(pname//'cannot allocate ttau_fde',0, MA_ERR)
c
        if (.not.MA_Push_Get(mt_dbl,nbf_fde_mxnbf_ce*n_per_rec,
     &     'Bmat_fde',lBmat_fde,iBmat_fde))
     &     call errquit(pname//'cannot allocate Bmat_fde',0, MA_ERR)
c
c Reduce density matrix for FDE
c
c     -compute reduced total density matrices over atoms
c
        if (.not. geom_ncent(geom_fde, ncenters_fde))
     &     call errquit(pname//'geom_ncent failed',73, GEOM_ERR)
c
c     Assume this numerical density routine can be called from anywhere.
c     -compute mapping vectors
c 
        if (.not.ma_push_get
     &     (mt_int,ncenters_fde*2,'cetobfr_fde map',
     &      lcetobfr_fde,icetobfr_fde))
     &     call errquit(pname//'push_get failed', 13, MA_ERR)
c Re-using deallocated variables
        if (.not.ma_push_get
     &     (mt_int,nshells_fde,'cntoce_fde map',lcntoce,icntoce))
     &     call errquit(pname//'push_get failed', 13, MA_ERR)
        if (.not.ma_push_get
     &     (mt_int,nshells_fde*2,'cntobfr_fde map',lcntobfr,icntobfr))
     &     call errquit(pname//'push_get failed', 13, MA_ERR)
c     
        call build_maps(FD_bas_han, int_mb(icntoce), int_mb(icntobfr),
     &                  int_mb(icetobfr_fde), ncenters_fde, nshells_fde)
        if (.not.ma_chop_stack(lcntoce))
     &     call errquit(pname//'cannot pop stack',1, MA_ERR)
c
        if (.not.MA_Push_Get(MT_Dbl,ipol*ncenters_fde*ncenters_fde,
     &      'rdens_atom_fde', lrdens_atom_fde,irdens_atom_fde))
     &  call errquit(pname//'cannot allocate rdens_atom_fde',0, MA_ERR)
        nscr = nbf_fde_mxnbf_ce*nbf_fde_mxnbf_ce*(nmat+1)
        if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &     call errquit(pname//'cannot allocate scr',0, MA_ERR)
c
c Daniel (1-11-13): I think we can stay similar to calc_type = 2 for
c XC-third derivatives in this portion of the code, since this part
c doesn't depend on having 3-dimensional GAs.
c
c AJL: Not possible inside this if statement
c        if (calc_type.eq.2.or.calc_type.eq.5) then
c           do isp=1,ipol
c              call util_ga_mat_reduce(nbf_fde, ncenters_fde,
c     &             int_mb(icetobfr_fde),
c     &             g_dens_fde(1+(isp-1)*(nmat+1)),  1,
c     &             dbl_mb(irdens_atom_fde+(isp-1)*
c     &                              ncenters_fde*ncenters_fde), 'rms',
c     &             dbl_mb(iscr), nbf_fde_mxnbf_ce,.true.)
c              do ii=2,(nmat+1)
c                 call util_ga_mat_reduce(nbf_fde,ncenters_fde,
c     &                int_mb(icetobfr_fde),
c     &                g_dens_fde(ii+(isp-1)*(nmat+1)),  1,
c     &                dbl_mb(irdens_atom_fde+(isp-1)*
c     &                              ncenters_fde*ncenters_fde), 'rms',
c     &                dbl_mb(iscr), nbf_fde_mxnbf_ce,.false.)
c              enddo
c           enddo
c        else
          call util_ga_mat_reduce(nbf_fde, ncenters_fde, 
     &                          int_mb(icetobfr_fde), g_dens_fde,
     &                          ipol, dbl_mb(irdens_atom_fde), 'rms',
     &                          dbl_mb(iscr), nbf_fde_mxnbf_ce,.true.)
c        endif ! calc_type
c
        if (.not.ma_pop_stack(lscr))
     &       call errquit(pname//'cannot pop stack',0, MA_ERR)

        if (.not.MA_Push_Get(mt_dbl,ipol*n_per_rec,'Amat_fde',
     &                       lAmat_fde,iAmat_fde))
     &     call errquit(pname//'cannot allocate Amat_fde',0, MA_ERR)
c
        if (.not.MA_Push_Get(mt_dbl,len_Amat2,'Amat2_fde',
     &                       lAmat2_fde,iAmat2_fde))
     &     call errquit(pname//'cannot allocate Amat2_fde',0, MA_ERR)
        call dfill(len_Amat2, 0.0d0, dbl_mb(iAmat2_fde), 1)
        if (.not.MA_Push_Get(mt_dbl,len_Amat3,'Amat3_fde',
     &                       lAmat3_fde,iAmat3_fde))
     &    call errquit(pname//'cannot allocate Amat3_fde', 0, MA_ERR)
        call dfill(len_Amat3, 0.0d0, dbl_mb(iAmat3_fde), 1)
c
        if (.not.MA_Push_Get(
     &       mt_dbl,3*n_per_rec*ipol,'Cmat_fde',lCmat_fde,iCmat_fde))
     &     call errquit(pname//'cannot allocate Cmat_fde',0, MA_ERR)
c
        if (.not.MA_Push_Get(mt_dbl,len_Cmat2,'Cmat2_fde',
     &                       lCmat2_fde,iCmat2_fde))
     &     call errquit(pname//'cannot allocate Cmat2_fde',0, MA_ERR)
        call dfill(len_Cmat2, 0.0d0, dbl_mb(iCmat2_fde), 1)
        if (.not.MA_Push_Get(mt_dbl,len_Cmat3,'Cmat3_fde',
     &                       lCmat3_fde,iCmat3_fde))
     &    call errquit(pname//'cannot allocate Cmat3_fde', 0, MA_ERR)
        call dfill(len_Cmat3, 0.0d0, dbl_mb(iCmat3_fde), 1)
c
        if (.not.MA_Push_Get(
     &       mt_dbl,len_mmat,'Mmat',lMmat_fde,iMmat_fde))
     &     call errquit(pname//'cannot allocate Mmat',0, MA_ERR)
c
        if (.not.MA_Push_Get(mt_dbl,len_Mmat2,'Mmat2_fde',
     &                       lMmat2_fde,iMmat2_fde))
     &     call errquit(pname//'cannot allocate Mmat2',0, MA_ERR)
c QM/FDE
        if (.not.MA_Push_Get(MT_Dbl,len_rho,'rho_qmfde',
     &     lrho_qmfde,irho_qmfde))
     &     call errquit(pname//'cannot allocate rho_qmfde',0, MA_ERR)
c
        if (.not.MA_Push_Get(MT_Dbl,len_delrho,'delrho_qmfde',
     &     ldelrho_qmfde, idelrho_qmfde))
     &     call errquit(pname//'cannot allocate delrho_qmfde',0, MA_ERR)
c
        if (.not.MA_Push_Get(MT_Dbl,len_ttau,'ttau_qmfde',
     &     lttau_qmfde, ittau_qmfde))
     &     call errquit(pname//'cannot allocate ttau_qmfde',0, MA_ERR)
c
        if (fde_option.eq.4.or.fde_option.eq.7) then
          do i=1,ipol
            if (.not. ga_duplicate(g_vxc(i), g_fde(i), 'g_fde matrix'))
     &      call errquit(pname//'failed duplicate',g_vxc, GA_ERR)
            call ga_zero(g_fde(i))
          enddo
        endif ! fde_option.eq.4.or.7
c
        if (.not.MA_Push_Get(mt_int, nbf_fde, 'ibf_fde', l_ibf_fde,
     &     i_ibf_fde))
     &     call errquit(pname//'cannot allocate fde nbf',0, MA_ERR)
        if (.not.MA_Push_get(MT_int,ncenters_fde,'fde nbf',
     &      lfin_fde,ifin_fde))
     &      call errquit(pname//'cannot allocate fde nbf',0, MA_ERR)
c
        if(.not.MA_push_get(MT_int, ncenters_fde, 'iniz_fde',
     &     lniz_fde, iniz_fde))
     &     call errquit(pname//"cannot allocate iniz_fde",0, MA_ERR)
        if (.not.MA_push_Get(MT_Dbl,3*ncenters_fde,
     &                      'xyz_fde',lxyz_fde,ixyz_fde))
     &     call errquit(pname//'cannot allocate xyz_fde',0, MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,ncenters_fde,
     &                      'charge_fde',lcharge_fde,icharge_fde))
     &     call errquit(pname//'cannot allocate charge_fde',0, MA_ERR)
        if (.not.MA_Push_Get(MT_Byte,ncenters_fde*16,
     &                      'tags_fde',ltags_fde,itags_fde))
     &     call errquit(pname//'cannot allocate tags_fde',0, MA_ERR)
        if (.not. geom_cart_get(geom_fde, ncenters_fde,
     &     Byte_MB(itags_fde), Dbl_MB(ixyz_fde), Dbl_MB(icharge_fde)))
     &     call errquit(pname//'geom_cart_get failed',74, GEOM_ERR)
c
c      compute the basis functions over the grid
         if(.not.MA_Push_get(MT_dbl, ncenters_fde, 'rchi_atom_fde',
     &         l_rchi_atom_fde,i_rchi_atom_fde))
     &         call errquit(pname//"rchi_atom_fde",0, MA_ERR)
         if(grad) then
           if (.not.MA_Push_Get(MT_Dbl,ncenters_fde,'rdelchi_atom_fde',
     &        lrdelchi_atom_fde,irdelchi_atom_fde))call errquit
     &        (pname//'cannot allocate rdelchi_atom_fde',0, MA_ERR)
         endif
         if(.not.MA_Push_get(MT_dbl, n_per_rec*ncenters_fde, 'rq_fde',
     &         l_rq_fde,i_rq_fde))
     &         call errquit(pname//"rq_fde",0, MA_ERR)
c
      endif ! fde_option.eq.3.4.6.or.7
cc AJL/END
c
      if (.not.MA_Push_Get(mt_dbl,len_Amat,'Amat',lamat,iamat))
     &   call errquit(pname//'cannot allocate Amat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Amat2,'Amat2',l_Amat2,k_Amat2))
     &   call errquit(pname//'cannot allocate Amat2',0, MA_ERR)
      call dfill(len_Amat2, 0.0d0, dbl_mb(k_Amat2), 1)
      if (.not.MA_Push_Get(mt_dbl,len_Amat3,'Amat3',l_Amat3,k_Amat3))
     1  call errquit(pname//'cannot allocate Amat3', 0, MA_ERR)
      call dfill(len_Amat3, 0.0d0, dbl_mb(k_Amat3), 1)
c
      if (.not.MA_Push_Get(
     .     mt_dbl,3*n_per_rec*ipol,'Cmat',lcmat,icmat))
     &   call errquit(pname//'cannot allocate Cmat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Cmat2,'Cmat2',l_Cmat2,k_Cmat2))
     &   call errquit(pname//'cannot allocate Cmat2',0, MA_ERR)
      call dfill(len_Cmat2, 0.0d0, dbl_mb(k_Cmat2), 1)
      if (.not.MA_Push_Get(mt_dbl,len_Cmat3,'Cmat3',l_Cmat3,k_Cmat3))
     1  call errquit(pname//'cannot allocate Cmat3', 0, MA_ERR)
      call dfill(len_Cmat3, 0.0d0, dbl_mb(k_Cmat3), 1)
c
      if (.not.MA_Push_Get(
     .     mt_dbl,len_mmat,'Mmat',lmmat,immat))
     &   call errquit(pname//'cannot allocate Mmat',0, MA_ERR)
c
      if (.not.MA_Push_Get(mt_dbl,len_Mmat2,'Mmat2',l_Mmat2,k_Mmat2))
     &   call errquit(pname//'cannot allocate Mmat2',0, MA_ERR)
      call dfill(len_Mmat2, 0.0d0, dbl_mb(k_Mmat2), 1)
c
      if (.not.MA_Push_Get(mt_dbl,len_Lmat,'Lmat',lLmat,iLmat))
     &   call errquit(pname//'cannot allocate Lmat',0, MA_ERR)
c Daniel (1-11-13): Eventually we'll need to put stuff here for third order
c functional derivative arrays of the kinetic energy density, but I'm
c avoiding these because we don't yet have second derivatives for
c meta-GGAs (as far as I'm aware).
c
c     allocate scratch array for grid compression
c
      nibf=nbf_ao
      if(xcfit) nibf=max(nibf,nbf_xc)
      if (.not.MA_Push_Get(mt_int, nibf, 'ibf_ao', libf_ao,
     &     iibf_ao))
     &     call errquit(pname//'cannot allocate ibf_ao',2,
     &       MA_ERR)
      if (.not.MA_Push_get(MT_dbl, 4*buffer_size+4, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit(pname//'cannot allocate quad pt buffer', 3,
     &       MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rchi_atom',
     &     lrchi_atom,irchi_atom))
     &     call errquit(pname//'cannot allocate rchi_atom',0, MA_ERR)
      if(grad) then
      if (.not.MA_Push_Get(MT_Dbl,nctrs,'rdelchi_atom',
     &     lrdelchi_atom,irdelchi_atom))call errquit
     &        (pname//'cannot allocate rdelchi_atom',0, MA_ERR)
      endif
      irqlen=nctrs*n_per_rec
      if(.not.grid_written) irqlen=max(irqlen,nctrs*nqmax)
      if (.not.MA_Push_get(mt_dbl,irqlen,'rq',lrq,irq))
     &   call errquit(pname//'cannot allocate rq',0, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,irqlen,'bp',lbp,ibp))
     &     call errquit(pname//'cannot allocate bp',4, MA_ERR)
c
c     compute variables for scratch space in function evaluation
c     assumes cartesian; spherical set would only require 2*MXang + 1
c     
      max_component = (nbf_ao_mxang+1)*(nbf_ao_mxang+2)/2 
      if (xcfit) max_component = max(max_component,
     $  (nbf_cd_mxang+1)*(nbf_cd_mxang+2)/2) 
c     
c     Allocate scratch space.
c     
      if (nbf_ao_mxang .eq. 0)then
         max_elem_ang_scr = max_component * 3
      elseif (nbf_ao_mxang .le. 3)then
         max_elem_ang_scr = max_component * 9
      else                      ! general case
         max_elem_ang_scr = max_component * 28
      endif 

      if (xcfit) then
        if (nbf_cd_mxang.le.3) then
          max_elem_ang_scr = max(max_elem_ang_scr, max_component*9)
        else
          max_elem_ang_scr = max(max_elem_ang_scr, max_component*28)
        endif
      endif

      nscr = 3*n_per_rec*nbf_ao_mxcont + 
     &     max(3*n_per_rec*nbf_ao_mxprim,
     &     n_per_rec*max_elem_ang_scr) + 1
      if (xcfit) nscr=max(nscr, 3*n_per_rec*nbf_xc_mxcont + 
     &     max(3*n_per_rec*nbf_xc_mxprim,
     &     n_per_rec*max_elem_ang_scr) + 1)
cc AJL/Begin/FDE
      if (frozemb_fde) 
     &     nscr=max(nscr, 3*n_per_rec*nbf_fde_mxcont +
     &     max(3*n_per_rec*nbf_fde_mxprim,
     &     n_per_rec*max_elem_ang_scr) + 1)
cc AJL/End
c     
      if (.not.MA_Push_Get(mt_dbl,nscr,'func scr',lscr,iscr))
     &     call errquit(pname//'cannot allocate bp',4, MA_ERR)
c
      Ex=0d0
      Ec=0d0
      call ga_sync()
c
      if(.not.do_cpks_l.or..not.do_gxc)
     &  call xc_rep_init(rtdb, geom, ao_bas_han,g_dens,g_vxc,
     &  nbf_ao,ipol,(do_scf.or.do_xdm).and.(.not.l3d),.false.)
c      
c     Grid quadrature
c
      if(grid_written) then
c
         if (iAOacc .le. 0) then
          if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc))then
             iAOacc=-nint(log(e_conv))
          else
             iAOacc=max(iAOacc,-nint(log(e_conv)))
          endif
         end if
c
         call grid_quadv0a(g_dens,ncoef_max,nctrs,rho_int,nqtot, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     .        dbl_mb(iscr),nscr,
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        dbl_mb(irdens_atom), int_mb(icetobfr),
     ,        int_mb(iibf_ao), dbl_mb(irqbuf),
     ,        dbl_mb(irchi_atom), dbl_mb(irdelchi_atom),
     ,        ex,ec,dbl_mb(iamat), dbl_mb(icmat),g_vxc,
     ,        dbl_mb(k_amat2), dbl_mb(k_cmat2),
     ,        nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm,
     T        tmat, l3d,triplet,
     ,        dbl_mb(ittau), 
     &        dbl_mb(immat), dbl_mb(k_mmat2),
     &        dbl_mb(iLmat),
     &        StericEnergyTotal,
     1        do_gxc, dbl_mb(k_Amat3), dbl_mb(k_Cmat3),
cc AJL/Begin/FDE
     &        fde_option,
     &        g_dens_fde, ncenters_fde, dbl_mb(irdens_atom_fde),
     &        dbl_mb(iFmat_fde),dbl_mb(iPmat_fde),dbl_mb(iBmat_fde),
     &        dbl_mb(ittau_fde),dbl_mb(irho_fde), dbl_mb(idelrho_fde),
     &        int_mb(icetobfr_fde),
     &        dbl_mb(iAmat_fde), dbl_mb(iAmat2_fde), dbl_mb(iAmat3_fde),
     &        dbl_mb(iCmat_fde), dbl_mb(iCmat2_fde), dbl_mb(iCmat3_fde),
     &        dbl_mb(iMmat_fde), dbl_mb(iMmat2_fde),
     &        dbl_mb(ittau_qmfde),dbl_mb(irho_qmfde), 
     &        dbl_mb(idelrho_qmfde), len_ttau, g_fde,
     &        int_mb(iniz_fde), int_mb(ifin_fde), int_mb(i_ibf_fde), 
     &        dbl_mb(ixyz_fde), dbl_mb(i_rchi_atom_fde),
     &        dbl_mb(irdelchi_atom_fde), dbl_mb(i_rq_fde), rho_int_fde)
cc AJL/End
      else
c
         if (.not. rtdb_get(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &        call errquit(pname//'rtdb_get failed', 119, RTDB_ERR)
         nqlen = ncenters*mxnrad
         if (.not.MA_Push_get(MT_int,4*nqlen,'shell list',
     &        lqs,iqs))
     &        call errquit(pname//'cannot allocate shell list',0,
     &       MA_ERR)
         if (.not.MA_Push_get(MT_dbl, 4*nqmax,
     &        'ang quad pts buffer', laqbuf, iaqbuf))
     &        call errquit(pname//'cant allocate quad pt', 0, MA_ERR)
         ndbl = (ncenters*(ncenters+1))/2
         if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lRij,iRij))
     &        call errquit(pname//'cannot allocate Rij',0, MA_ERR)
         if (.not.MA_Push_get(mt_int,nqshells,'qshat',l_qat,k_qat))
     &        call errquit(pname//'cannot allocate qat',5, MA_ERR)
c
c Daniel (1-11-13): Added XC-third derivative stuff here also, although
c it should never be used when we need to construct the grid.  This is
c needed for preventing odd behavior in calls to the grid_quadv0b and
c xc_eval_fnl routine later on, when someone implements new functionality
c down the line for the DFT code.
         call grid_quadvw(rtdb,
     ,        int_mb(iqs),dbl_mb(iaqbuf),
     ,        dbl_mb(irij),int_mb(k_qat),
     ,        g_dens,ncoef_max,nctrs,rho_int,nqtot, 
     &        int_mb(iniz), int_mb(ifin), 
     &        dbl_mb(irho), dbl_mb(idelrho), dbl_mb(ihess), 
     &        dbl_mb(iqwght), dbl_mb(iqxyz), dbl_mb(ixyz), 
     .        dbl_mb(iscr),nscr,
     &        dbl_mb(iexpo), dbl_mb(icoef), 
     &        dbl_mb(iBmat), dbl_mb(iFmat), 
     &        dbl_mb(iPmat), dbl_mb(irq), 
     &        dbl_mb(irdens_atom), int_mb(icetobfr),
     ,        int_mb(iibf_ao),dbl_mb(ibp),
     ,        dbl_mb(irchi_atom),dbl_mb(irdelchi_atom),
     ,        ex,ec,dbl_mb(iamat),dbl_mb(icmat),g_vxc,
     ,        dbl_mb(k_amat2),dbl_mb(k_cmat2),
     ,        nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm,
     T        tmat, l3d,triplet,
     ,        dbl_mb(ittau), 
     &        dbl_mb(immat), dbl_mb(k_mmat2),
     &        dbl_mb(iLmat),
     1        do_gxc, dbl_mb(k_Amat3), dbl_mb(k_Cmat3),
     2        0, 0.0d0, 0.0d0, 
cc AJL/Begin/FDE
     &        fde_option,
     &        g_dens_fde, ncenters_fde, dbl_mb(irdens_atom_fde),
     &        dbl_mb(iFmat_fde),dbl_mb(iPmat_fde),dbl_mb(iBmat_fde),
     &        dbl_mb(ittau_fde),dbl_mb(irho_fde), dbl_mb(idelrho_fde),
     &        int_mb(icetobfr_fde),
     &        dbl_mb(iAmat_fde), dbl_mb(iAmat2_fde), dbl_mb(iAmat3_fde),
     &        dbl_mb(iCmat_fde), dbl_mb(iCmat2_fde), dbl_mb(iCmat3_fde),
     &        dbl_mb(iMmat_fde), dbl_mb(iMmat2_fde),
     &        dbl_mb(ittau_qmfde),dbl_mb(irho_qmfde),
     &        dbl_mb(idelrho_qmfde), len_ttau, g_fde,
     &        int_mb(iniz_fde), int_mb(ifin_fde), int_mb(i_ibf_fde),
     &        dbl_mb(ixyz_fde), dbl_mb(i_rchi_atom_fde),
     &        dbl_mb(irdelchi_atom_fde), dbl_mb(i_rq_fde), rho_int_fde)
cc AJL/End
c
         if (.not. MA_chop_Stack(lqs))
     &        call errquit(pname//'pop stack failed.',1, MA_ERR)
c Qin
         if (.not. rtdb_get(rtdb, 'dft:grid_make_sa', mt_log, 1, lsa))
     &     lsa = .false.
         if(lsa) then
           call grid_make_sa(rtdb)
           lsa = .false.
           if (.not. rtdb_put(rtdb, 'dft:grid_make_sa', mt_log, 1, lsa))
     &       call errquit(pname//'rtdb_put lsa failed', 0, RTDB_ERR)
         endif
c end Qin

      endif  ! grid_written

      if(.not.xc_rep_close(rtdb, nbf_ao,ipol,g_vxc,g_dens,.false.)) 
     &     call errquit(pname//'xcrepclose failed ',0, DISK_ERR)

c     
C AJL/Begin/FDE
c Calculate the non-additive potential by removing the potential for the
c FDE from the potential for the combined system. This cannot be
c achieved at a deeper level because the global arrays are only available
c here upwards
      if (fde_option.eq.4.or.fde_option.eq.7) then
        do i=1,ipol
          call ga_dadd(1d0, g_vxc(i), -1d0, g_fde(i), g_vxc(i))
          if (.not. ga_destroy(g_fde(i)))
     &    call errquit(pname//'failed destroy',g_fde, GA_ERR)
        enddo
      endif
c AJL/End
c
      if (.not.ma_chop_stack(lqxyz))
     &   call errquit(pname//'cannot pop stack',0, MA_ERR)
      call ga_sync()
      if((do_scf.or.do_xdm).and.(.not.l3d).and.(.not.adft)) then
c
c         write(0,*) ' symmetr scf l r',do_scf,do_cpks_l,do_cpks_r,
c     A        calc_type
c         call ga_print(g_vxc(1))
c
          call ga_symmetrize(g_vxc(1))
          if(ipol.gt.1) call ga_symmetrize(g_vxc(2))
      endif
c
      call ga_igop(Msg_nqtot,nqtot,1,'+')
      call ga_dgop(Msg_Excrho,dummy,3,'+')
c
cc AJL/Begin/FDE
      if (fde_option.eq.1) then
c Skip through as we are calculating ESP on the grid
        continue
      else
        if(xcfit.and.(.not.adft)) then
          call ga_dgop(Msg_Tmat,Tmat,nbf_xc*ipol,'+')
        elseif(adft) then
          call ga_dgop(Msg_Tmat,Tmat(nbf_xc*ipol+1),nbf_xc*ipol,'+')
        endif
        if (nExc.eq.1)then
           Exc(1) = Exc(1) + Ec + Ex
        else
           Exc(1) = Exc(1) + Ex
           Exc(2) = Exc(2) + Ec
        endif
c
c     == Print total steric energy ==
c
        if (lsteric) then
           if (me.eq.0) write(LuOut,2221) StericEnergyTotal
        end if
 2221   format('         Total Steric Energy =', f22.12/)
c
c     get grid accuracy for job
c
        if (.not. rtdb_get(rtdb, 'dft:job_grid_acc', mt_dbl, 1, 
     &     job_grid_acc))
     &     call errquit(pname//'rtdb_put failed', 911, RTDB_ERR)
c
        fon = .false.
c     fractional occupations?
c
        if (rtdb_get(rtdb,'dft:fon',mt_log,1,fon)) then 
          if (.not.fon) call errquit(
     &       pname//'fon stored in RTDB but not .true.', 1,
     &       RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'dft:pstrace', mt_dbl, 1, pstrace))
     &       call errquit(pname//'rtdb_get pstrace failed',
     &       1, RTDB_ERR)
        else
          fon = .false.
        end if
c
c     get number of electrons
c
        if(ipol.eq.1) then
           nelec=dble(2*noc(1))
        else
           nelec=dble(noc(1)+noc(2))
        endif
c
c     for fractional occupations the correct electron number comes from tr[PS]
        if (fon) nelec = pstrace
c         
        if(oprintw) then
           if(sic_orb_occ.eq.0) then
              rho_err=abs(rho_int-nelec)/nelec
           else
              rho_err=abs(rho_int-1d0)
           endif
c Daniel (1-11-13): I changed the logic here to stop the warning from
c being printed for either the TDDFT solver or TDDFT gradient solver.
           lcpksl_gxc = (do_cpks_l.or.do_gxc)
           if(.not.lcpksl_gxc.and.rho_err.gt.job_grid_acc.and.me.eq.0)
     &        write(LuOut,123) rho_int,job_grid_acc
c         if(.not.do_cpks_l.and.rho_err.gt.job_grid_acc.and.me.eq.0)
c     &      write(LuOut,123) rho_int,job_grid_acc
        endif ! oprintw

c     some debug output to check the density integration error:
c$$$      if (me.eq.0) then
c$$$        write(luout,'(1x,a,t15,i10)') "noc(1):",noc(1)
c$$$        write(luout,'(1x,a,t15,i10)') "noc(2):",noc(2)
c$$$        write(luout,'(1x,a,t15,i10)') "nelec:",nelec
c$$$        write(luout,'(1x,a,t15,f15.8)') "rho_int:",rho_int
c$$$        write(luout,'(1x,a,t15,f15.8)')
c$$$     &     "rhodiff: ",abs(rho_int-dble(nelec))
c$$$      end if
c
        rho_n = rho_int  ! integrated density
c
        if (me.eq.0.and.oprint)then
           write(LuOut,*)
           if (fde_option.eq.0) then
              write(LuOut,'(a,f14.6,a,i8,a,i8)')
     &        ' XC quadr. integrated density : ',
     &        rho_int,' npts: ',nqtot,' nctrs: ',nctrs
           else
              write(LuOut,'(a,f14.6,a,f14.6,a,i8,a,i8)')
     &        ' XC quadr. integrated density (QM): ', rho_int,
     &                                     ' (FDE): ',rho_int_fde,
     &        ' npts: ',nqtot,' nctrs: ',nctrs
           endif
           write(LuOut,*)
        endif
      endif !fde_option
cc AJL/End
c
      if (.not.ma_pop_stack(lz))
     &   call errquit(pname//'cannot pop stack',0, ma_err)
      if (.not.MA_Free_Heap(lxyz))
     &   call errquit(pname//'cannot free heap',0, MA_ERR)
c
 123  format(
     &   1x,'Grid integrated density: ',f20.12,/
     &   1x,'Requested integration accuracy: ',e10.2)
c
      return
      end

      subroutine adft_reduce(nbf, nat, cetobfr, coef, ipol, rdens_atom)
      implicit none
      integer nbf, nat, ipol, mxnbf
      integer cetobfr(2,nat)
      double precision rdens_atom(nat,nat,ipol),coef(nbf,ipol)
      double precision dabsmax
      external dabsmax

      integer iat,istart,iend

      rdens_atom(:,:,:) = 0d0

      do iat=1,nat
        istart = cetobfr(1,iat)
        iend   = cetobfr(2,iat)
        rdens_atom(1,iat,1) = dabsmax(iend-istart+1,coef(istart,1))
        if (ipol.gt.1) then
          rdens_atom(1,iat,2) = dabsmax(iend-istart+1,coef(istart,2))
        endif
      enddo

      end
