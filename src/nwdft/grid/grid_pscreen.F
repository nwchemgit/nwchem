      subroutine grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
     ,     rchi_atom,rdelchi_atom,rdens_atom,iniz,chi,delchi,
     ,     mbf_ao,grad)
c
c$Id: grid_pscreen.F,v 1.1 2000-04-11 20:49:40 edo Exp $
c
      implicit none
#include "mafdecls.fh"
      integer natoms
      integer nqpts ![in/out]
      double precision qwght(*),qxyz(3,*)! [in/out]
      double precision rchi_atom(*),rdelchi_atom(*)
      double precision chi(*),delchi(3,*)
      double precision rdens_atom(natoms,natoms,*)
      integer iniz(*),mbf_ao
      double precision tol_rho
      logical grad
c
      double precision dabsmax
      external dabsmax
c
      integer k_iscr,l_iscr,k_dscr,l_dscr
      integer nzero,i
      double precision pmax,funcmax,ff,delmax
c
      nzero=0
c
c     get pmax ang funcmax
c
      pmax=0d0
      funcmax=dabsmax(natoms,rchi_atom)
      delmax=funcmax
      if(grad) delmax=max(delmax,dabsmax(natoms,rdelchi_atom))
      do i=1,natoms
         if(iniz(i).ne.0) then
             pmax=max(pmax,dabsmax(natoms,rdens_atom(1,i,1)))
          endif
      enddo
      ff=abs(funcmax*delmax*pmax)
      do i=1,nqpts
         if(ff*qwght(i).lt.tol_rho) then
            nzero=nzero+1
            qwght(i)=0d0
         endif
      enddo
      if(nzero.ne.0) then
            if (.not.MA_Push_Get(mt_int,nqpts,
     &         'iscr', l_iscr,k_iscr))
     .        call errquit('grid_pcom: pushget failed',1)
            k_dscr=0
            if(grad) then
            if (.not.MA_Push_Get(mt_dbl,nqpts*3*mbf_ao,
     &         'dscr', l_dscr,k_dscr))
     .        call errquit('grid_pcom: pushget failed',2)
            endif
c        write(0,'(" % ",f12.6,I4,2X,1PE12.6)') 
c     .           (nzero*100d0)/nqpts,nzero,tol_rho
         call grid_pcompress(nqpts,qwght,chi,qxyz,delchi,
     ,           int_mb(k_iscr),dbl_mb(k_dscr),
     ,        mbf_ao,1d-19,grad)
            if (.not.ma_chop_stack(l_iscr))
     &           call errquit('grid_quadv0_a: cannot pop stack', 11)
      endif
      return
      end
      subroutine grid_pcompress(nq,qwght,chi,qxyz,delchi,
     ,     iscr,dscr,
     .     mbf_ao,toll,grad)
      implicit none
      integer nq,mbf_ao
      double precision qwght(*),chi(*),qxyz(3,*),delchi(*)
      double precision dscr(*)
      integer iscr(*)
      logical grad
      double precision toll ![input]
c
      integer nsign,i,q,j
c        
c       compress grid
c
      if(nq.eq.0 )return
      nsign=0
      do i=1,nq
      if(abs(qwght(i)).gt.toll) then
       nsign=nsign+1
       iscr(nsign)=i
       endif
      enddo
      if(nq.eq.nsign) return
      if(nsign.eq.0)  then
         nq=0
         return
      else
      do q=1,nsign
        qxyz(1,q)=qxyz(1,iscr(q))
        qxyz(2,q)=qxyz(2,iscr(q))
        qxyz(3,q)=qxyz(3,iscr(q))
        qwght(q)=qwght(iscr(q))
      enddo
        do j=1,mbf_ao
           do q=1,nsign
              chi(q+(j-1)*nsign)=chi(iscr(q)+(j-1)*nq)
           enddo
        enddo
        if(grad) then
           call dcopy(3*nq*mbf_ao,delchi,1,dscr,1)
           call grid_pcompres1(nsign,nq,mbf_ao,delchi,dscr,iscr)
        endif
        nq=nsign
      endif
      return
      end
      subroutine grid_pcompres1(nsign,nq,mbf_ao,delchi,dscr,iscr)
      implicit none
      integer nsign,nq,mbf_ao
      integer iscr(*)
      double precision delchi(nsign,3,*)
      double precision dscr(nq,3,*)
c
      integer j,q
c
      do j=1,mbf_ao
         do q=1,nsign
            delchi(q,1,j)=dscr(iscr(q),1,j)
            delchi(q,2,j)=dscr(iscr(q),2,j)
            delchi(q,3,j)=dscr(iscr(q),3,j)
         enddo
      enddo
      return
      end
