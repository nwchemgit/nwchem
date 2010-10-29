      subroutine diis_core_so(derr,pulay,
     &     mxhist,nhist,icall,nbf,g_diis,g_err)
C$Id$
      implicit none
c
      integer mxhist
      integer nbf
      integer g_diis(2)
      integer g_err(2)
      double precision pulay((mxhist+1)*(mxhist+2)/2)
      Integer nhist, icall
c***********************************************************************
c
c     subroutine diis
c
c     Direct inversion of the iterative subspace.
c
c     Purpose:  Construct Pulays B matrix from the trace of the inner 
c               product of error vectors.
c                              t                      
c               B_ij = trace(Ei ,Ej)
c           
c
c     References:
c                 Peter Pulay, Convergence Acceleration of Iterative
c                 Sequences.  The Case of SCF Iterations.
c                 Chem. Phys. Lett.  73, 393 (1980)
c               
c                 and 
c
c                 Peter Pulay, Improved SCF Convergence Acceleration,
c                 J. Comp. Chem. 3, 566, (1982).
c
c                 
c     Variables:
c                 errvec:  error vector.
c                 mxhist:  maximum history.
c                 nhist:   current number of objects in history.
c                 icall:   number of calls made to this routine.
c                 pulay:   Pulays B matrix.
c                 scr:     scratch array
c                 maxhist: hardwired upper-bound to the maximum allowed
c                          history.
c
c***********************************************************************
c
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
c
      integer irow,ijcanon,i,j
      integer me,nproc,i1,ii,istep,ijm1,ijlast,jj,ij
      double precision derr
      logical oprint
c
c
      irow(j) = mod(j-1,mxhist)+1
      ijcanon(i,j) = i*(i-1)/2 +j
      oprint = util_print('convergence', print_default)
      me=ga_nodeid()
      nproc=ga_nnodes()
c
      if(icall.eq.1) then
c
c     Initialize the Pulay matrix.  Store -1 on the 1st row and column.
c     (1,1) element = 0.
c     
c     if(me.eq.0)write(LUout,*)' Initializing DIIS; ispin = ',ispin

         pulay(1) = 0.0d0
c     
         do i = 2,mxhist+1
            i1=ijcanon(i,1)
            pulay(i1) = -1.0d0
         enddo 
c     
         nhist=0
      endif 
      ii=mxhist + irow(icall)
      istep=(ii-1)*nbf
      call ga_sync
c     
C     Store new errvec
c     
      call ga_copy_patch('N',
     *     g_err(1),1,nbf,1,nbf, 
     *     g_diis(1),1,nbf,istep+1,istep+nbf)
      call ga_copy_patch('N',
     *     g_err(2),1,nbf,1,nbf, 
     *     g_diis(2),1,nbf,istep+1,istep+nbf)
      
c     
c     If more than mxhist histories are stored, shift the
c     rows of the Pulays B matrix.
c     
      if(nhist.eq.mxhist)then
c     
         do i = 3,nhist+1
            do j = 3,i
c     
               ij = ijcanon(i,j)
               ijm1 = ijcanon(i-1,j-1)
               pulay(ijm1)=pulay(ij)
            enddo 
         enddo 
c     
      endif 
c     
c     Determine the size of the iterative subspace history.
c     
      nhist = min(nhist+1,mxhist)
c     
c     Construct the last column of the B matrix.
c     
      derr=ga_ddot(g_err(1),g_err(1))+ga_ddot(g_err(2),g_err(2))
      do jj = 1,nhist
c     
         ii=mxhist+irow(icall-jj+1)
         istep=(ii-1)*nbf
         call ga_sync
         ijlast = ijcanon(nhist+1,nhist+2-jj)
         pulay(ijlast)=
     &        ga_ddot_patch(g_diis(1),'N',1,nbf,istep+1,istep+nbf,
     *        g_err(1),'N',1,nbf,1,nbf)+
     &        ga_ddot_patch(g_diis(2),'N',1,nbf,istep+1,istep+nbf,
     *        g_err(2),'N',1,nbf,1,nbf)
c     
      enddo
      return
C     1111 format(/,15x,'Diis error vector norm:   ',f20.10,/)
      end
      
      
