#define NO_FDE_ARGS 1
      Subroutine dftg_loop(
     ,     nq_out,iqsh,irsh,rad_sh,dosympts,
     ,     fnxyz,nxyz,nops,grad,
     ,     nq_cube,nqtotc,symfact,
     ,     iang_done,
     ,     gridscr,
     ,     iscr_loc,symops,gops,
     ,     qsh,aqbuf,rij, 
     ,     g_dens, ncoef_max, natoms, 
     &     rho_int,
     &     iniz, ifin, rho, delrho, 
     &     hess, qwght, qxyz, xyz,iscr,
     &     expo, Bmat, Fmat, 
     &     Pmat,  rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,do_xdm,tmat,
     ,     ttau, kske, Mmat, Mmat2,
     ,     kslap, Lmat,
     D     func,d_qwght,rq0,hh,ppp,ff,zzz,
     D     drhonuc,ddelrhonuc,d2scr,force,
     D     acc_ao_gauss,
     L     do_force,do_hess,ldew,gradchi,hesschi,der3chi,oactive,
     I     iandex,max_at_bf,curatoms,g_rhs,
     I     nbf_ao_mxnbf_ce2,nactive,
     1     do_gxc, Amat3, Cmat3, ngxc,
     2     rhosave, delrhosave, triplet,
#ifdef NO_FDE_ARGS
     3     do_tddftvxc)
#else
cc AJL/Begin/FDE
     &     do_tddftvxc, fde_option,
     &     g_dens_fde, ncenters_fde, rdens_atom_fde,
     &     Fmat_fde, Pmat_fde, Bmat_fde,
     &     ttau_fde, rho_fde, delrho_fde,
     &     cetobfr_fde,
     &     Amat_fde, Amat2_fde, Amat3_fde,
     &     Cmat_fde, Cmat2_fde, Cmat3_fde,
     &     Mmat_fde, Mmat2_fde,
     &     ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &     iniz_fde, ifin_fde, ibf_fde,
     &     xyz_fde, rchi_atom_fde, rdelchi_atom_fde, rq_fde, 
     &     rho_int_fde)
c
cc AJL/End
c
#endif
c$Id$
c
      implicit none
c
#define USE_BP
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "sym.fh"
#include "errquit.fh"
#include "bas.fh"
#define NO_FDE_INCL
#include "grid_quad.fh"
#include "grid_cube.fh"
c
      double precision rad_sh  ! [in]
      integer irsh,iqsh ! [in]
      logical dosympts ! [in]
      double precision  gridscr(*)
      integer iscr_loc(*)
      integer nxyz , fnxyz ![in]
      double precision symops(*),gops(*)
      integer qsh(3,*)! [in]
      double precision aqbuf(4,*) ! [scr]
      integer nops              ! [in]
      integer iang_done ![in/out]
      double precision symfact ![in]
      double precision rij ! [scr] for ssw
      integer nqtotc,nq_cube ! [in/out]
c
      integer nq_left ! [local]
      integer npts_per_ang
      integer iang,ictr,nq
      integer nq_out
      integer grid_nxyz
      external grid_nxyz
c
      logical kske, kslap
c
      double precision func(*),d_qwght(*),
     D     rq0(*),hh(*),ppp(*),ff(*),zzz(*),force(*)
      double precision drhonuc(*),ddelrhonuc(*),d2scr(*),
     ,     acc_ao_gauss
      integer iandex(*),max_at_bf,curatoms(*),g_rhs(*)
      logical do_force,do_hess,ldew,gradchi,hesschi,der3chi,
     ,     oactive(*)
      integer nbf_ao_mxnbf_ce2,nactive
c
      integer lbas_cent_info, ibas_cent_info,
     &        lbas_cset_info, ibas_cset_info,
     &        ldocset, idocset
      integer ncontrset,n
c Qin
      integer rtdb
      logical lcdft, ldbw
      integer ipop, nconstr, i
      integer iatms1(max_constr), iatms2(max_constr)
      integer fatms1(max_constr), fatms2(max_constr)
      integer dtype(max_constr)
      double precision multipliers(max_constr), fact(max_constr)  ! pass to dftg_gridv0b
c
      save nq_left
c
c Daniel (1-16-13): For consistency, we pass the number of perturbed
c density matrices to this routine, although they should never be
c used for anything in dftg_gridv0b
      integer ngxc
      double precision rhosave(*)
      double precision delrhosave(*)
      logical triplet
      logical do_tddftvxc ! dVxc*P contribution to the TDDFT grads
      integer basis

      if (adft) then
        basis = xc_bas_han
      else
        basis = ao_bas_han
      endif
c
c Qin
      ldbw = .false.
      lcdft=.false.
#if 0
      if (.not. rtdb_get(rtdb, 'dft:lcdft', mt_log, 1, lcdft))
     &     lcdft = .false.
      if(lcdft) then
        if (.not. rtdb_get(rtdb, 'dft:ipop', mt_int, 1, ipop))
     &     ipop = 3
        if(ipop.eq.1) then
         ldbw = .true.
         if (.not. rtdb_get(rtdb, 'dft:nconstr', mt_int, 1, nconstr))
     &     call errquit('dft_gridv0a: get nconstr failed', 0, RTDB_ERR)
        if (.not. (
     $            rtdb_get(rtdb,'dft:iatms1', mt_int, nconstr,iatms1(1))
     $      .and. rtdb_get(rtdb,'dft:iatms2', mt_int, nconstr,iatms2(1))
     $      .and. rtdb_get(rtdb,'dft:fatms1', mt_int, nconstr,fatms1(1))
     $      .and. rtdb_get(rtdb,'dft:fatms2', mt_int, nconstr,fatms2(1))
     $      ))
     &   call errquit('dft_gridv0a: atms failed', 2700, RTDB_ERR)
        if (.not.rtdb_get(rtdb, 'dft:multipliers', mt_dbl, nconstr,
     $       multipliers(1)))
     &    call errquit('dft_gridv0a: get multipliers failed', 0)
        if (.not. rtdb_get(rtdb,'dft:dtype', mt_int, nconstr, dtype(1)))
     &   call errquit('dft_gridv0a: get dtype failed', 3000, RTDB_ERR)
        endif
      endif
#endif
c Qin / end
c
      npts_per_ang=nqmax
      ictr = qsh(2,iqsh)
      iang = qsh(3,iqsh)
      nq = npts_per_ang
c
      call grid_finpts(dosympts,
     .        rad_sh,xyz, ictr, 
     .        irsh, iang, iang_done,nq_left,
     ,        aqbuf, nq, fnxyz,nxyz,
     ,        gridscr,nops,symops,gops,
     .        iscr_loc,
     .     ncenters, rij,
     .     rq,bp,expo,
     .     oskel,symfact)
      do n=1,nq
         qxyz(1,n) = aqbuf(1,n)
         qxyz(2,n) = aqbuf(2,n)
         qxyz(3,n) = aqbuf(3,n)
         qwght(n) = aqbuf(4,n)
      enddo
c
      if(nq.ne.0) then
      if (.not.bas_numcont(basis, ncontrset))
     &     call errquit('Exiting in dftg_gridv0a',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('dftg_gridv0a: cannot allocate docset',
     ,     ncontrset,0)
      if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('dftg_gridv0a: cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('dftg_gridv0a: cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(basis, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natoms)
c     
c     copy to buffer for XC stuff and compute if needed
c     
         nqtotc=nqtotc+nq
         nq_cube=nq_cube+nq
         nq_out=nq
         call dftg_gridv0b(nq_out,rad_sh,ictr,g_dens,
     &     max_at_bf,
     &     rho_int, nq, iniz, ifin, rho, delrho,
     &     func, d_qwght, qwght, qxyz, xyz, bp,
     &     iandex, expo, Amat, Amat2, Cmat,
     &     Cmat2, Bmat, Fmat, Pmat,  rq,
     &     rq0, hh, 
     .     ibf_ao,int_mb(ibas_cent_info),int_mb(ibas_cset_info),
     .     log_mb(idocset),curatoms,rchi_atom,rdelchi_atom,ppp,
     .     ff,zzz,drhonuc,ddelrhonuc,d2scr,
     ,     force, hess, g_rhs, do_force,
     ,     nbf_ao_mxnbf_ce2,
     ,     acc_ao_gauss, ldew,grad,gradchi,hesschi,der3chi,
     &     do_hess, do_cpks_r, natoms, oactive,
     &     nactive, rdens_atom, 
     &     cetobfr,
     &     ttau, kske, Mmat, Mmat2,
     &     kslap, Lmat,
     &     ldbw, fact, dtype, nconstr,  ! Qin
     &     do_gxc, Amat3, Cmat3, ngxc,
     &     rhosave, delrhosave, triplet,
     &     do_tddftvxc,tmat)

c         
         if (.not.MA_Pop_Stack(lbas_cset_info))
     &        call errquit('dftg_gridv0a: pop stack failed.',2, MA_ERR)
         if (.not.MA_Pop_Stack(lbas_cent_info))
     &        call errquit('dftg_gridv0a: pop stack failed.',3, MA_ERR)
         if (.not.MA_Pop_Stack(ldocset))
     &        call errquit('dftg_gridv0a: pop stack failed.',4, MA_ERR)
      endif
c
      return
      end
