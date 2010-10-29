      subroutine diis_refvec_so(nbf,g_err,rowre,colre,rowim, colim)
C     $Id$
      implicit none
      
      integer nbf
      integer g_err(2)
      double precision rowre(nbf),colre(nbf),rowim(nbf),colim(nbf)
      double precision diag

#include "global.fh"
#include "tcgmsg.fh"
c     
      integer nproc,igrad,next,i,j
      integer nxtask
      external nxtask
c     
c      write(*,*)"g_err first"
      nproc=ga_nnodes() 
      if(ga_nodeid().eq.0)then 
         call ga_put(g_err(1),1,1,1,1,0.d0,1)
         call ga_get(g_err(2),1,1,1,1,diag,1)
         diag=-2.0d0*diag
         call ga_put(g_err(2),1,1,1,1,diag,1)
      endif
      igrad=1
      next=nxtask(nproc,igrad)+2
      do i = 2,nbf
         if(i.eq.next) then
            call ga_get(g_err(1),1,i-1,i,i,colre,1)
            call ga_get(g_err(1),i,i,1,i-1,rowre,1)
            call ga_get(g_err(2),1,i,i,i,colim,1)
            call ga_get(g_err(2),i,i,1,i-1,rowim,1)
            rowre(i)=0.d0
            colre(i)=0.d0
            colim(i)=-2.d0*colim(i)
            do j=1,i-1
               rowre(j)=-rowre(j)+colre(j)
               colre(j)=-rowre(j)
               rowim(j)=-rowim(j)-colim(j)
               colim(j)=rowim(j)
            enddo
            call ga_put(g_err(1),1,i,i,i,colre,1)
            call ga_put(g_err(1),i,i,1,i-1,rowre,1)
            call ga_put(g_err(2),1,i,i,i,colim,1)
            call ga_put(g_err(2),i,i,1,i-1,rowim,1)
            next=nxtask(nproc,igrad)+2
         endif
      enddo
      next=nxtask(-nproc,igrad)
c      write(*,*)"g_err anti"
c     
      return
      end
