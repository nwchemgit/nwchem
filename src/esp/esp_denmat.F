      subroutine esp_denmat(geom,basis,g_vec,g_dns,occ)
c
c $Id: esp_denmat.F,v 1.3 1998-08-13 21:35:18 d3j191 Exp $
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "esp_params.fh"
#include "esp_common.fh"
c
      integer geom,basis
      integer g_vec,g_dns
      real*8 occ(nbf)
c
      integer ilo,ihi,jlo,jhi,lthi,lthj,i
      integer g_scr,k_scr,l_scr,icol
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      real*8 scale
      logical lresult
c
      g_scr=ga_create_atom_blocked(geom, basis, 'scr')
c
      call ga_copy(g_vec,g_scr)
      call ga_sync()
c
      call ga_distribution(g_scr,me,ilo,ihi,jlo,jhi)
      lthi=ihi-ilo+1
      lthj=jhi-jlo+1
c
      if(lthi.ne.0.and.lthj.ne.0) then
      if(lthi.lt.0.or.lthj.lt.0) 
     + call errquit('esp_denmat: distribution error',me)
c
      if(.not.ma_push_get(mt_dbl,lthi,'vector',l_scr,k_scr))
     + call errquit('esp_denmat: failed to allocate vector',me)
c
      do 1 icol=jlo,jhi
      call ga_get(g_scr,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx      write(*,6666) icol,(dbl_mb(k_scr+i),i=0,nbf-1)
cx 6666 format(//,'Vector',i5,/,(10f10.6))
      scale=occ(icol)
cx      write(*,'(a,i5,f12.6)') 'Scale ',icol,scale
      call dscal(lthi,scale,dbl_mb(k_scr),1)
      call ga_put(g_scr,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
      call ga_get(g_scr,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx      write(*,6667) (dbl_mb(k_scr+i),i=0,nbf-1)
cx 6667 format(10f10.6)
    1 continue
c
      if(.not.ma_pop_stack(l_scr))
     + call errquit('esp_denmat: failed to deallocate vector',me)
c
      endif
c
      call ga_sync()
      call ga_dgemm('n','t',nbf,nbf,nbf,1.d0,g_vec,g_scr,0.d0,g_dns)
c
      if(.not.ga_destroy(g_scr))
     +  call errquit('esp_denmat: failed to destroy g_scr',me)
c
      return
      end
