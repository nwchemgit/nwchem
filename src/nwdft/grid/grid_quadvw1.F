      Subroutine grid_quadvw1(ncu,nq_out,iqsh,irsh,rad_sh,dosympts,
     ,     fnxyz,nxyz,nops,grad,
     ,     nq_cube,nqtotc,rad_write,symfact,
     ,     iang_done,nwrbuf,
     ,     gridscr,iscr_loc,symops,gops,
     ,     qsh, aqbuf,rij, 
     ,     g_dens, ncoef_max, natoms, 
     &     rho_int, nqtot, 
     &     iniz, ifin, rho, delrho, 
     &     hess, qwght, qxyz, xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm, tmat,
     ,     l3d, triplet, ttau, kske, Mmat, Mmat2,
     ,     kslap, Lmat,
     1     do_gxc, Amat3, Cmat3,
cc AJL/Begin/FDE
     &     fde_option, g_dens_fde, ncenters_fde, rdens_atom_fde,
     &     Fmat_fde, Pmat_fde, Bmat_fde,
     &     ttau_fde, rho_fde, delrho_fde, cetobfr_fde,
     &     Amat_fde, Amat2_fde, Amat3_fde,
     &     Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     &     ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &     iniz_fde, ifin_fde, ibf_fde, xyz_fde, 
     &     rchi_atom_fde, rdelchi_atom_fde, rq_fde, rho_int_fde)
cc AJL/End
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
      double precision rad_write ! [in]
c
      integer nq_left ! [local]
      integer ntot_cube,npts_per_ang
      integer iang,ictr,nq,ncu
      integer nq_out,nwrbuf
      integer grid_nxyz,nq_in
      external grid_nxyz
      logical l3d
      logical triplet
      logical  grid_file_write
      external grid_file_write
c
      logical kske, kslap
c
      save nq_left
c
      npts_per_ang=nqmax
      ictr = qsh(2,iqsh)
      iang = qsh(3,iqsh)
      nq = npts_per_ang
c
      call grid_finpts(dosympts,
     .     rad_sh,xyz, ictr, 
     .     irsh, iang, iang_done,nq_left,
     ,     aqbuf, nq, fnxyz,nxyz,
     ,     gridscr,nops,symops,gops,
     .     iscr_loc,
     .     ncenters, rij,
     .     rq,bp,expo,
     .     oskel,symfact)
c
      if(nq.ne.0) then
c     
c     copy to buffer for XC stuff and compute if needed
c     
         nqtotc=nqtotc+nq
         nq_cube=nq_cube+nq
         nq_in=nq
         call grid_quadvc(.false.,nwrbuf,
     ,        aqbuf,nq_in,nq_out,
     ,        ntot_cube,ictr,qxyz, qwght, 
     ,        rad_write,grad,
     .        g_dens, ncoef_max, natoms, 
     &        rho_int, nqtot, 
     &        iniz, ifin, rho, delrho, 
     &        hess,  xyz,iscr,nscr,
     &        expo, ccoef, Bmat, Fmat, 
     &        Pmat, rq, rdens_atom, cetobfr, 
     ,        ibf_ao,rchi_atom,rdelchi_atom,
     ,        Ex, Ec, Amat,Cmat,g_vxc,Amat2,Cmat2,
     ,        nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm,
     T        tmat, l3d, triplet,
     &        ttau, kske, Mmat, Mmat2,
     &        kslap, Lmat,
     1        do_gxc, Amat3, Cmat3,
cc AJL/Begin
     &        fde_option, g_dens_fde, ncenters_fde, rdens_atom_fde,
     &        Fmat_fde, Pmat_fde, Bmat_fde,
     &        ttau_fde, rho_fde, delrho_fde, cetobfr_fde,
     &        Amat_fde, Amat2_fde, Amat3_fde,
     &        Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     &        ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &        iniz_fde, ifin_fde, ibf_fde, xyz_fde, rchi_atom_fde, 
     &        rdelchi_atom_fde, rq_fde, rho_int_fde)
cc AJL/FDE
c     
c     now pack this aqbuf into the write buffer 
c     (write and update task counter as needed)
c     
         if(.not.largenode)
     .        grid_written = grid_file_write(nq, ictr, 
     ,        rad_write, aqbuf,  ncu)
         
      endif
c
      return
      end
c $Id$
