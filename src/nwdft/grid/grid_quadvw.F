      Subroutine grid_quadvw(rtdb, 
     ,     qsh,aqbuf,
     ,     rij, iqat,
     ,     g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess, qwght, qxyz, xyz,iscr,nscr,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat,
     .     l3d,triplet, tau,Mmat)
c
c$Id: grid_quadvw.F,v 1.38 2004-03-19 03:25:15 edo Exp $
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
      integer ntot_cube,npts_per_ang
      logical xc_chkgrad, kske, xc_chktau
      external xc_chkgrad, xc_chktau
      logical grid_file_open
      external grid_file_open
      logical grid_test_mem
      external grid_test_mem
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

      double precision tau(*), Mmat(*)
C      
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
      call grid_loop(nops,
     ,     oprint,grad,g_dens,
     ,     dosympts,sub4,largenode,grid_written,
     , ncoef_max,natoms,
     ,     nqtot,
     &                       iniz, ifin, rho, delrho, 
     &                       hess, qwght, qxyz, xyz,iscr,nscr,
     &                    expo, ccoef, Bmat, Fmat, 
     &                    Pmat, rq,  rdens_atom, cetobfr, 
     ,                    ibf_ao,bp,rchi_atom,rdelchi_atom,
     ,                    Ex, Ec, amat,cmat,g_vxc,
     ,                    Amat2, Cmat2,
     ,                    nmat, do_scf, do_cpks_l, do_cpks_r,tmat,
     ,                    l3d,triplet,
     &                    tau, kske,Mmat,
     ,     iqat,
     ,     qsh, aqbuf,rij,rho_int)

c
      return
      end
