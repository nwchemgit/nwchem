      Subroutine grid_quadvw1(ncu,nq_out,iqsh,irsh,rad_sh,dosympts,
     ,     fnxyz,nxyz,nops,grad,
     ,     nq_cube,nqtotc,rad_write,symfact,
     ,     iang_done,nwrbuf,
     ,     gridscr,iscr_loc,symops,
     ,     qsh, aqbuf,rij, 
     ,     g_dens, ncoef_max, natoms, 
     &     rho_int, nqtot, 
     &     iniz, ifin, rho, delrho, 
     &     hess, qwght, qxyz, xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq,  rdens_atom, cetobfr, 
     ,     ibf_ao,bp,rscr,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat,
     ,     donewbf)
c
c$Id: grid_quadvw1.F,v 1.6 2000-10-17 21:21:05 edo Exp $
c
      implicit none
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "sym.fh"
#include "grid_quad.fh"
#include "grid_cube.fh"
      double precision rad_sh  ! [in]
      integer irsh,iqsh ! [in]
      logical dosympts ! [in]
      double precision  gridscr(*)
      integer iscr_loc(*)
      integer nxyz , fnxyz ![in]
      double precision symops(*)
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
      logical donewbf           !!! !!!! *** temp
c
      save nq_left
c
      npts_per_ang=nqmax
      ictr = qsh(2,iqsh)
      iang = qsh(3,iqsh)
      nq = npts_per_ang
c     
c     Generate all angular quadrature points for this radial shell
c     and keep only those inside cube
c     
      if(dosympts) then
         call grid_sympts(dosympts,
     .        rad_sh,xyz(1,ictr), ictr, 
     .        irsh, iang, iang_done,nq_left,
     ,        aqbuf, nq, fnxyz,nxyz,
     ,        gridscr,nops,symops,
     .        iscr_loc)
      else
         call grid_qpts(rad_sh,xyz(1,ictr), ictr, 
     ,        irsh, iang, 
     .        aqbuf, nq, 
     *        fnxyz,nxyz)
      endif
c     
c     modify the weight to include the spatial weight (and symmetry factor)
c     
      if(ncenters.gt.1) then
         call grid_setspac(xyz, aqbuf,
     ,        rij, 
     &        rq, bp, ictr, nq,
     ,        expo, rad_sh)
      endif
      if(nq.ne.0) then
         if(oskel) call dscal(
     /        nq, symfact, aqbuf(4,1), 4)
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
     ,        ibf_ao,bp,rscr,rchi_atom,rdelchi_atom,
     ,        Ex, Ec, amat,cmat,g_vxc,amat2,cmat2,
     ,        nmat, do_scf, do_cpks_l, do_cpks_r,tmat,
     ,        donewbf)
         
c     
c     now pack this aqbuf into the write buffer 
c     (write and update task counter as needed)
c     
         call grid_file_write(nq, ictr, 
     ,        rad_write, aqbuf,  ncu)
         
      endif
      return
      end
