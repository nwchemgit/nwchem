      Subroutine grid_quadv0a(g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess, qwght, qxyz, xyz,iscr,nscr,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq, rdens_atom, cetobfr, 
     ,     ibf_ao,rqbuf,bp,rscr,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,
     ,     Amat2, Cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat)
c
c$Id: grid_quadv0a.F,v 1.22 2000-08-01 20:11:54 edo Exp $
c
      implicit none
c
#include "cgridfile.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "grid_quad.fh"
c
      double precision rqbuf(*)
      double precision rad
      
      integer iqsh
      integer nqpts, ictr_buf
      integer ncube,istep,ntot_cube
      logical grid_file_rewind
      external grid_file_rewind
      logical xc_chkgrad
      external xc_chkgrad
c
c     see if any Gradient Corrected Functional is not null
c
      grad=xc_chkgrad(xfac,cfac)
      if(n_rec_in_file.eq.0) goto 100
c
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit('grid_quadv0a: rewinding gridpts?', 0)
C
      ntot_cube=0
      do 200 iqsh = 1, n_rec_in_file
c
c     Define the current range of radial shells and integration center.
c     
         call grid_file_read(n_per_rec, nqpts, ictr_buf, 
     &        rad,rqbuf,nsubb)
         if(nqpts.gt.buffer_size) call
     '        errquit(' buffersize exceed by qpts ',nqpts)
        if(nqpts.eq.0) goto 200
        istep=0
        do  ncube=1,nsubb
c     
c        temp put buf into currently used arrays qxyz and qwght
c     
         call grid_repack(rqbuf, qxyz, qwght, nqpts,
     ,          rad,istep)
         if(nqpts.ne.0) then
            call grid_quadv0b(
     ,           ictr_buf,qxyz, qwght, nqpts,rad,grad,
     .           g_dens, ncoef_max, natoms, 
     &                       rho_int, nqtot, 
     &                       iniz, ifin, rho, delrho, 
     &                       hess,  xyz,iscr,nscr,
     &                       expo, ccoef, Bmat, Fmat, 
     &                       Pmat, rq, rdens_atom,  cetobfr, 
     ,           ibf_ao,bp,rscr,rchi_atom,rdelchi_atom,
     ,           Ex, Ec, amat,cmat,g_vxc,amat2,cmat2,
     ,           nmat, do_scf, do_cpks_l, do_cpks_r,tmat)
            ntot_cube=ntot_cube+1
         endif
      enddo
 200  continue
 100  continue
c
      call ga_sync
      return
      end
