      subroutine grid_pscreen(nqpts,qwght,qxyz,tol_rho,natoms,
     ,     rchi_atom,rdelchi_atom,rdens_atom,iniz,chi,delchi,
     ,     mbf_ao,grad,ipol)
c
c$Id: grid_pscreen.F,v 1.4 2001-06-16 00:00:47 edo Exp $
c
      implicit none
#include "mafdecls.fh"
      integer natoms
      integer nqpts ![in/out]
      double precision qwght(*),qxyz(3,*)! [in/out]
      double precision rchi_atom(*),rdelchi_atom(*)
      double precision chi(*),delchi(*)
      double precision rdens_atom(natoms,natoms,*)
      integer iniz(*),mbf_ao
      double precision tol_rho
      logical grad
      integer ipol
c
      double precision dabsmax
      external dabsmax
c
      integer nzero,i,j,ispin
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
      do ispin=1,ipol
         do i=1,natoms
            if(iniz(i).ne.0) 
     .           pmax=max(pmax,dabsmax(natoms,rdens_atom(1,i,ispin)))
         enddo
      enddo
      ff=abs(funcmax*delmax*pmax)
      do i=1,nqpts
         if(abs(ff*qwght(i)).lt.tol_rho) then
            nzero=nzero+1
            qwght(i)=0d0
c
c     invalidate all the functions
c
            do j=1,mbf_ao
               chi(i+(j-1)*nqpts)=0d0
               if(grad)then
                  delchi(i+(j-1)*3*nqpts)=0d0
                  delchi(i+(j-1)*3*nqpts+nqpts)=0d0
                  delchi(i+(j-1)*3*nqpts+nqpts*2)=0d0
               endif
            enddo
         endif
      enddo
      if(nzero.ne.0) then
         call grid_pcompress(nqpts,qwght,chi,qxyz,delchi,
     ,        mbf_ao,1d-19,grad)
      endif
      return
      end
      subroutine grid_pcompress(nq,qwght,chi,qxyz,delchi,
     .     mbf_ao,toll,grad)
      implicit none
      integer nq,mbf_ao
      double precision qwght(*),chi(*),qxyz(3,*),delchi(*)
      logical grad
      double precision toll ![input]
c
      integer i,q,nqxmbf,nzero,nq_orig,j,qj
c        
c       compress grid
c
      if(nq.eq.0 )return
      nq_orig=nq
      nqxmbf=nq*mbf_ao
c
c     BFs
c
      nzero=0
      qj=0
      do j=1,mbf_ao
         do q=1,nq
            qj=qj+1
            if(abs(qwght(q)).lt.toll) then
               nzero=nzero+1
            else
               if(nzero.ne.0) chi(qj-nzero)=chi(qj)
            endif
         enddo
      enddo
      if(grad) then
         nzero=0
         qj=0
         do j=1,mbf_ao*3
            do q=1,nq
               qj=qj+1
               if(abs(qwght(q)).lt.toll) then
                  nzero=nzero+1
               else
                  if(nzero.ne.0) delchi(qj-nzero)=delchi(qj)
               endif
            enddo
         enddo
      endif
c
c     grid pts
c
      i=1
      nzero=0
      do q=1,nq
         if(abs(qwght(q)).lt.toll) then
            nzero=nzero+1
         else
            if(nzero.ne.0) then
               qwght(q-nzero)=qwght(q)
               qxyz(1,q-nzero)=qxyz(1,q)
               qxyz(2,q-nzero)=qxyz(2,q)
               qxyz(3,q-nzero)=qxyz(3,q)
            endif
         endif
      enddo
      nq=nq-nzero
#if defined(DEBUG)
      if(nqxmbf/mbf_ao.ne.nq) then
      write(0,*) ' NONO: bf ',nqxmbf/mbf_ao,' ww ',nq
      endif
#endif
      return
      end
