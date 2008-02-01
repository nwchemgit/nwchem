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
     ,     nmat, do_scf, do_cpks_l, do_cpks_r, do_nmrks_r, tmat,
     .     l3d,triplet, tau, ttau, Mmat, Mmat2)
c
c$Id: grid_quadvw.F,v 1.46 2008-02-01 22:27:45 niri Exp $
c
      implicit none
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
      integer rtdb ![in]
      integer qsh(3,*)! [in]
      double precision aqbuf(4,*) ! [scr]
      double precision rij ! [scr] for ssw
      integer iqat(*) ! [scr]
      integer npts_per_ang
      logical xc_chkgrad, kske, xc_chktau
      external xc_chkgrad, xc_chktau
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
c
      nproc=ga_nnodes()
      me=ga_nodeid()
      oprint = util_print('quadrature', print_high)
     
c
c     see if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad()
      kske=xc_chktau()
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
     ,     do_nmrks_r, tmat,
     ,     l3d,triplet,
     &     tau, ttau, kske, Mmat, Mmat2,
     ,     iqat,
     ,     qsh, aqbuf,rij,rho_int,
     D     ddum1,ddum1,ddum1,ddum1,ddum1,ddum1,ddum1,
     D     ddum1,ddum1,ddum1,ddum1,
     D     ddum1,
     L     ldum1,ldum1,ldum1,ldum1,ldum1,ldum1,ldum1,
     I     idum1,idum1,idum1,idum1,
     I     idum1,idum1,idum1)
c
      return
      end
