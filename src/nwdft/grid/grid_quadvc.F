      subroutine grid_quadvc(cube_done,nwrbuf,
     .     aqbuf,nq,nq_out,
     ,     ntot_cube,ictr_buf,qxyz, qwght, 
     ,     rad_write,grad,
     .     rtdb, g_dens, ncoef_max, natoms, 
     &     rho_int, nqtot, 
     &     iniz, ifin, rho, delrho, 
     &     hess,  xyz,iscr,nscr,
     &     expo, ccoef, Bmat, Fmat, 
     &     Pmat, rq, rdens_atom, 
     &     cntobfr, cetobfr, 
     ,     ibf_ao,bp,rscr,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,amat2,cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat)
C$Id: grid_quadvc.F,v 1.6 2000-04-21 00:28:46 edo Exp $
      implicit none
#include "cdft.fh"
#include "cgridfile.fh"
#include "grid_quad.fh"
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
      integer n
      logical buffwrite
c
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
!      write(6,11) ictr_buf,ntot_cube,nq_out,nwrbuf,cube_done,
!     ,        buffwrite
! 11   format(' ictr ',i2,' nrchi ',i4,' nq_out ',i4,' nwrb ',i4,
!     ,     ' cd ',i1,' bw ',i1) 
         if(nq_out.ne.0) 
     .    call grid_quadv0b(
     ,        ntot_cube,ictr_buf,qxyz, qwght, nq_out,
     ,        rad_write,grad,
     .        rtdb, g_dens, ncoef_max, natoms, 
     &        rho_int, nqtot, 
     &        iniz, ifin, rho, delrho, 
     &        hess,  xyz,iscr,nscr,
     &        expo, ccoef, Bmat, Fmat, 
     &        Pmat, rq,  rdens_atom, 
     &        cntobfr, cetobfr, 
     ,        ibf_ao,bp,rscr,rchi_atom,rdelchi_atom,
     ,        Ex, Ec, amat,cmat,g_vxc,amat2,cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r)
         ntot_cube=ntot_cube+1
         if(cube_done) then
            nwrbuf=nwrbuf+1
cedo            nwrbuf=1
            return
         endif
         if(buffwrite) then
            nwrbuf=1
         endif
         nq_out=0
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
      return
      end
