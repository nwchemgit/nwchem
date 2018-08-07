      Subroutine grid_quadvw(rtdb, 
     ,     qsh,aqbuf,
     ,     rij, iqat,
     ,     g_dens, ncoef_max, natoms, 
     &     rho_int, nqtot, 
     &     iniz, ifin, rho, delrho, 
     &     hess, qwght, qxyz, xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, do_xdm,
     T     tmat, l3d,triplet, ttau, Mmat, Mmat2,
     ,     Lmat,
     1     do_gxc, Amat3, Cmat3,
     2     td_idum, td_ddum, td_ddum1,
cc AJL/Begin/FDE
     &     fde_option,
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
cc AJL/End
c
c$Id$
c
      implicit none
c
#include "errquit.fh"
#define USE_BP
#include "cgridfile.fh"
#include "stdio.fh"
#include "geom.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "grid_quad.fh"
#include "grid_cube.fh"
#include "msgids.fh"
c
      integer rtdb ![in]
      integer qsh(3,*)! [in]
      double precision aqbuf(4,*) ! [scr]
      double precision rij ! [scr] for ssw
      integer iqat(*) ! [scr]
      integer npts_per_ang
      logical xc_chkgrad, kske, xc_chktau, kslap, xc_chklap
      external xc_chkgrad, xc_chktau, xc_chklap
      integer me,nproc
      integer nops ! [out]
      logical dosympts
      integer grid_nxyz,grid_nsub4
      external grid_nxyz,grid_nsub4
      double precision grid_lcube
      external grid_lcube
      logical oprint
      logical l3d
      logical triplet
c
      double precision ddum1
      integer idum1
      logical ldum1
      integer td_idum
      double precision td_ddum, td_ddum1
c
      nproc=ga_nnodes()
      me=ga_nodeid()
      oprint = util_print('quadrature', print_high)
     
c
c     see if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad()
      kske=xc_chktau()
      kslap=xc_chklap()
      npts_per_ang=nqmax
c
c     initialize the grid
c
      call grid_init(rtdb,geom,
     ,     ncenters,
     ,     oskel,leb,dosympts,qsh)
c
      call grid_loop(.false.,nops,
     ,     oprint,grad,g_dens,
     ,     dosympts,sub4,largenode,grid_written,
     ,     ncoef_max,natoms,
     ,     nqtot, iniz, ifin, rho, delrho, 
     &     hess, qwght, qxyz, xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, Amat,Cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,
     ,     do_nmrks_r, do_xdm, tmat,
     ,     l3d,triplet,
     &     ttau, kske, Mmat, Mmat2,
     &     kslap, Lmat,
     ,     iqat,
     ,     qsh, aqbuf,rij,rho_int,
     D     ddum1,ddum1,ddum1,ddum1,ddum1,ddum1,ddum1,
     D     ddum1,ddum1,ddum1,ddum1,
     D     ddum1,
     L     ldum1,ldum1,ldum1,ldum1,ldum1,ldum1,ldum1,
     I     idum1,idum1,idum1,idum1,
c     I     idum1,idum1,idum1,
     I     idum1,idum1,
     1     do_gxc, Amat3, Cmat3,
     2     td_idum, td_ddum, td_ddum1,
     3     .false.,   !  false is for the dVxc*P contribution to the TDDFT grads
cc AJL/Begin/FDE
     &     fde_option, g_dens_fde, ncenters_fde, rdens_atom_fde,
     &     Fmat_fde, Pmat_fde, Bmat_fde,
     &     ttau_fde, rho_fde, delrho_fde, cetobfr_fde,
     &     Amat_fde, Amat2_fde, Amat3_fde,
     &     Cmat_fde, Cmat2_fde, Cmat3_fde, Mmat_fde, Mmat2_fde,
     &     ttau_qmfde, rho_qmfde, delrho_qmfde, len_ttau, g_fde,
     &     iniz_fde, ifin_fde, ibf_fde,
     &     xyz_fde, rchi_atom_fde, rdelchi_atom_fde, rq_fde,rho_int_fde)
cc AJL/End
c
      return
      end
