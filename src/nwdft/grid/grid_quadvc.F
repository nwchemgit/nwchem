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
     ,     ibf_ao,bp,rscr,rchi_atom,rdelchi_atom,
     ,     Ex, Ec, amat,cmat,g_vxc,amat2,cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat)
C$Id: grid_quadvc.F,v 1.9 2000-09-21 15:23:53 bjohnson Exp $
      implicit none
#include "cdft.fh"
#include "mafdecls.fh"
#include "bas.fh"
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
      integer ncontrset, n
      integer lbas_cent_info, ibas_cent_info,
     &        lbas_cset_info, ibas_cset_info,
     &        ldocset, idocset
      logical buffwrite
c
c     Allocate and create info for new basis function evaluator
c
      if (.not.MA_Push_Get(mt_int, 3*natoms, 'bas_cent_info',
     &     lbas_cent_info, ibas_cent_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cent_info',0)
      if (.not.bas_numcont(AO_bas_han, ncontrset))
     &     call errquit('Exiting in grid_quadv0a',0)
      if (.not.MA_Push_Get(mt_int, 6*ncontrset, 'bas_cset_info',
     &     lbas_cset_info, ibas_cset_info))
     &     call errquit('grid_quadv0a: cannot allocate bas_cset_info',0)
      if (.not.MA_Push_Get(mt_log, ncontrset, 'docset',
     &     ldocset, idocset))
     &     call errquit('grid_quadv0a: cannot allocate docset',0)
      call xc_make_basis_info(AO_bas_han, int_mb(ibas_cent_info),
     &     int_mb(ibas_cset_info), natoms)
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
     ,        ictr_buf,qxyz, qwght, nq_out,
     ,        rad_write,grad,
     .        g_dens, ncoef_max, natoms, 
     &        rho_int, nqtot, 
     &        iniz, ifin, rho, delrho, 
     &        hess,  xyz,iscr,nscr,
     &        expo, ccoef, Bmat, Fmat, 
     &        Pmat, rq,  rdens_atom, cetobfr, 
     ,        ibf_ao,int_mb(ibas_cent_info),int_mb(ibas_cset_info),
     '        log_mb(idocset),bp,rscr,rchi_atom,rdelchi_atom,
     ,        Ex, Ec, amat,cmat,g_vxc,amat2,cmat2,
     ,     nmat, do_scf, do_cpks_l, do_cpks_r,tmat)
         ntot_cube=ntot_cube+1
         if(cube_done) then
            nwrbuf=nwrbuf+1
cedo            nwrbuf=1
            goto 100
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
 100  continue
      if (.not.MA_Pop_Stack(ldocset))
     &     call errquit('grid_quadvc: pop stack failed.',3)
      if (.not.MA_Pop_Stack(lbas_cset_info))
     &     call errquit('grid_quadvc: pop stack failed.',2)
      if (.not.MA_Pop_Stack(lbas_cent_info))
     &     call errquit('grid_quadvc: pop stack failed.',1)
c
      return
      end
