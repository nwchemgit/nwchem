      subroutine esp_denmat(geom,basis,g_vec,g_dns,nocc,nbf)
      implicit none
*
************************************************************************
*     include standard input and memory checking common structures.
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
************************************************************************
*     
*---- arguments of a dmat call
      integer geom,basis
      integer g_vec,g_dns,nbf
      real*8 nocc(nbf)
*---- local variables
      integer iproc,ilo,ihi,jlo,jhi,lthi,lthj
      integer g_scr,k_scr,l_scr,icol
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      real*8 scale
      logical lresult,debug
*
************************************************************************
*     
      debug = .false.
      iproc = ga_nodeid()
*
*---- allocate the global array
c      call gacrea(geom,basis,nbf,nbf,'scr',g_scr,'atom')
      g_scr = ga_create_atom_blocked(geom, basis, 'scr')
      call ga_copy(g_vec,g_scr)
      call ga_sync
*
*---- find the distribution of the g_scr
      call ga_distribution(g_scr,iproc,ilo,ihi,jlo,jhi)
      lthi = ihi - ilo + 1
      lthj = jhi - jlo + 1
      if (debug) then
         write(*,*)' dmat: iproc,ilo,ihi,lthi ',iproc,ilo,ihi,lthi
         write(*,*)' dmat: iproc,jlo,jhi,lthj ',iproc,jlo,jhi,lthj
      end if
      if (lthi.eq.0 . or . lthj.eq.0) goto 9999
      if (lthi.lt.0 . or . lthj.lt.0)
     &   call errquit('dmat: wrong distribution !!!',0)
*
*---- allocate the local vector
      lresult = ma_push_get(mt_dbl,lthi,'vector',l_scr,k_scr)
      call errmem(lresult,'dmat: failed to allocate vector',0)
*
*---- scale each column of g_scr by corresponding occupation number
      do icol = jlo, jhi
         call ga_get(g_scr,ilo,ihi,icol,icol,dbl_mb(k_scr),lthi)
         scale = nocc(icol)
         call dscal(lthi,scale,dbl_mb(k_scr),1)
         call ga_put(g_scr,ilo,ihi,icol,icol,dbl_mb(k_scr),lthi)
      end do
*
*---- deallocate vector
      lresult = ma_pop_stack(l_scr)
      call errmem(lresult,'dmat: failed to deallocate vector',0)
*
 9999 continue
      call ga_sync
      call ga_dgemm('n','t',nbf,nbf,nbf,1.d0,g_vec,g_scr,0.d0,g_dns)
*
c      call gadest(g_scr)
      if (.not. ga_destroy(g_scr)) call errquit
     &   ('dmat: could not destroy g_scr', 0)
*
*---- normal termination
      return
      end
