      subroutine grid_quadvc(cube_done,nwrbuf,
     .     aqbuf,nq,nq_out,
     ,     ntot_cube,ictr_buf,qxyz, qwght, 
     ,     rad_write,grad,
     .     g_dens, ncoef_max, natoms, 
     &     rho_int, nqtot, 
     &     iniz, ifin, rho, delrho, 
     &     hess,  xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq, rdens_atom, cetobfr, 
     ,     ibf_ao,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, Amat,Cmat,g_vxc,Amat2,Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm, tmat,
     ,     l3d, triplet, 
     &     ttau, kske, Mmat, Mmat2,
     &     kslap, Lmat,
     1     do_gxc, Amat3, Cmat3, 
cc AJL/Begin/FDE
     &     fde_option, g_dens_fde, ncenters_fde, rdens_atom_fde,
     &     Fmat_fde, Pmat_fde, Bmat_fde,ttau_fde, rho_fde, delrho_fde,
     &     cetobfr_fde, Amat_fde, Amat2_fde, Amat3_fde,
     &     Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     &     ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &     iniz_fde, ifin_fde, ibf_fde, xyz_fde, 
     &     rchi_atom_fde, rdelchi_atom_fde, rq_fde, rho_int_fde)
cc AJL/End
C$Id$
      implicit none
#include "errquit.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "cgridfile.fh"
#include "grid_quad.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
c
c     copy to buffer for XC and compute XC if npts=buff_size
c
      logical cube_done ! [in]
      integer nq        ! [in]
      integer nwrbuf    ! [in/out] buffer simulation
      double precision rad_write ! [in]
      integer ntot_cube ! [in] cube identific no. for rchi
      double precision aqbuf(4,*) ! [scr]
      integer nq_out ! [in/out] incremented
      integer ictr_buf
      integer ncontrset, n
      integer lbas_cent_info, ibas_cent_info,
     &        lbas_cset_info, ibas_cset_info,
     &        ldocset, idocset,i_iscratch,l_iscratch
      integer ncontrsetx,lbas_cent_xinfo, ibas_cent_xinfo,
     &        lbas_cset_xinfo, ibas_cset_xinfo,
     .     ldocsetx, idocsetx
      logical buffwrite
      logical l3d
      logical triplet
      logical kske, kslap
      integer k_becke
      integer l_becke
      integer basis
      logical oprint_qbecke
      double precision StericEnergy
c
cc AJL/Begin/FDE
      integer ncontrset_fde
      integer lbas_cset_info_fde, ibas_cset_info_fde
      integer lbas_cent_info_fde, ibas_cent_info_fde
      integer ldocset_fde, idocset_fde
      character*32 pname

      pname = 'grid_quadvc: '
      ibas_cset_info_fde = 0
      ibas_cent_info_fde = 0
      idocset_fde = 0
cc AJL/End
c
      basis = ao_bas_han
      if(adft) basis = xc_bas_han
c     becke charges
c
      oprint_qbecke = util_print('becke charges', print_high)
      if(oprint_qbecke) 
     R    call dft_initqbecke(natoms,k_becke,l_becke)
      ibas_cent_xinfo = 0
      ibas_cset_xinfo = 0
      idocsetx = 0
c
c     Allocate and create info for new basis function evaluator
c
      if (.not.bas_numcont(basis, ncontrset))
     &     call errquit('Exiting in grid_quadv0c',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit(pname//'cannot allocate bas_cent_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit(pname//'cannot allocate bas_cset_info',0,
     &       MA_ERR)
      call xc_make_basis_info(basis, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natoms)
      if (xcfit.and.(.not.adft)) then
         if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_xinfo',
     &        lbas_cent_xinfo, ibas_cent_xinfo))
     &        call errquit(pname//'cannot alloc bas_cent_info',0,
     &       MA_ERR)
      if (.not.bas_numcont(xc_bas_han, ncontrsetx))
     &     call errquit('Exiting in grid_quadv0a',0, BASIS_ERR)
      if (.not.MA_Push_Get(mt_int, 6*ncontrsetx, 'bas_cset_xinfo',
     &     lbas_cset_xinfo, ibas_cset_xinfo))
     &     call errquit(pname//'cannot allocate bas_cset_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrsetx, 'docsetx',
     &     ldocsetx, idocsetx))
     &     call errquit(pname//'cannot allocate docset',
     .     ncontrset, MA_ERR)
      call xc_make_basis_info(xc_bas_han, int_mb(ibas_cent_xinfo),
     &     int_mb(ibas_cset_xinfo), natoms)
      endif ! xcfit
c
      if (.not.MA_Push_Get(mt_int, 6*n_per_rec, 'iscratch',
     &     l_iscratch, i_iscratch))
     &     call errquit(pname//'cannot allocate bas_cset_info',0,
     &       MA_ERR)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit(pname//'cannot allocate docset',
     .     ncontrset, MA_ERR)
cc AJL/Begin/FDE
      if (fde_option.eq.3.or.fde_option.eq.4.or.
     &    fde_option.eq.6.or.fde_option.eq.7) then
       if (.not.bas_numcont(FD_bas_han, ncontrset_fde))
     &     call errquit(pname//'bas_numcont',0, BASIS_ERR)
       if (.not.MA_Push_Get(mt_int, 3*ncenters_fde, 'bas_cent_info_fde',
     &   lbas_cent_info_fde, ibas_cent_info_fde))
     &   call errquit(pname//'cannot allocate bas_cent_info_fde',0,
     &       MA_ERR)
       if (.not.MA_Push_Get(mt_int, 6*ncontrset_fde,'bas_cset_info_fde',
     &   lbas_cset_info_fde, ibas_cset_info_fde))
     &   call errquit(pname//'cannot allocate bas_cset_info_fde',0,
     &       MA_ERR)
       call xc_make_basis_info(FD_bas_han, int_mb(ibas_cent_info_fde),
     &     int_mb(ibas_cset_info_fde), ncenters_fde)
       if (.not.MA_Push_Get(mt_log, ncontrset_fde, 'docset_fde',
     &     ldocset_fde, idocset_fde))
     &     call errquit(pname//'cannot allocate docset_fde',
     .     ncontrset_fde, MA_ERR)
      endif !fde_option.eq.3,4,6 or 7
cc AJL/End

c
c     pack aqbuf to qxyz and compute when bufsize exceeded
c
      n = 0
 30   continue
      buffwrite=((n.eq.0).and.(nwrbuf+1.ge.n_per_rec-1)
     .     .and.(nq_out.ne.0)).or.
     .   (n.eq.0).and.((n_in_rec+nq.gt.n_per_rec)
     .     .and.(n_in_rec.gt.n_per_rec*3/5)
     .     ).or.
     .     ((n.gt.0).and.(nwrbuf+1.eq.n_per_rec))
      if(buffwrite.or.cube_done) then
         if(nq_out.ne.0) 
     .    call grid_quadv0b(
     ,        ictr_buf,qxyz, qwght, nq_out,
     ,        rad_write,grad,
     .        g_dens, ncoef_max, natoms, 
     &        rho_int, nqtot, 
     &        iniz, ifin, rho, delrho, 
     &        hess,  xyz,iscr,nscr,
     &        expo, ccoef, Bmat, Fmat, 
     &        Pmat, rq,  rdens_atom, cetobfr, 
     ,        ibf_ao,int_mb(ibas_cent_info),int_mb(ibas_cset_info),
     .        int_mb(i_iscratch),
     L        log_mb(idocset),rchi_atom,rdelchi_atom,
     ,        Ex, Ec, Amat,Cmat,g_vxc,Amat2,Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm,
     T        tmat, l3d,triplet,
     .     int_mb(ibas_cent_xinfo),int_mb(ibas_cset_xinfo),
     .        log_mb(idocsetx),
     &        ttau, kske, Mmat, Mmat2, 
     &        kslap, Lmat,
     &        StericEnergy,
     1        do_gxc, Amat3, Cmat3,
     Q        dbl_mb(k_becke),
cc AJL/Begin/FDE
     &        fde_option, g_dens_fde, ncenters_fde, rdens_atom_fde,
     &        Fmat_fde, Pmat_fde, Bmat_fde,
     &        ttau_fde, rho_fde, delrho_fde, cetobfr_fde,
     &        Amat_fde, Amat2_fde, Amat3_fde,
     &        Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     &        ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &        iniz_fde, ifin_fde, ibf_fde, xyz_fde, 
     &        rchi_atom_fde, rdelchi_atom_fde, rq_fde,
     &        int_mb(ibas_cent_info_fde), int_mb(ibas_cset_info_fde),
     &        log_mb(idocset_fde), rho_int_fde)
cc AJL/End
         ntot_cube=ntot_cube+1
         if(cube_done) then
            nwrbuf=nwrbuf+1
            goto 100
         endif
         if(buffwrite) then
            nwrbuf=1
            nq_out=0
         endif
      endif
      n=n+1
      nwrbuf=nwrbuf+1
      nq_out=nq_out+1
      qxyz(1,nq_out) = aqbuf(1,n)
      qxyz(2,nq_out) = aqbuf(2,n)
      qxyz(3,nq_out) = aqbuf(3,n)
      qwght(nq_out) = aqbuf(4,n)
!            write(0,'(A,i4,4F16.9)') 
!     .           ' RE ',n+nq_out,
!     ,       qxyz(1,n),qxyz(2,n),qxyz(3,n),qwght(n)
      if(n.lt.nq) goto 30            
!      write(0,*) ' endc ',nwrbuf+1,nq,nq_out
 100  continue
      if (.not.MA_chop_Stack(lbas_cent_info))
     &     call errquit(pname//'pop stack failed.',1, MA_ERR)
      if(oprint_qbecke) then
ctoomanyprints         call dft_printqbecke(geom,ao_bas_han,ncenters,dbl_mb(k_becke))
         if (.not.ma_pop_stack(l_becke)) 
     Z        call errquit(pname//'pop stack failed.',0, MA_ERR)
      endif
c
      return
      end
