      subroutine diis_core(derr,pulay,
     &     mxhist,nhist,icall,ispin,ipol, g_diis,g_err)
C$Id$
      implicit none
c
      integer mxhist
      integer ispin
      integer ipol
      integer g_diis(mxhist)
      integer g_err
      double precision pulay((mxhist+1)*(mxhist+2)/2)
      Integer nhist(2), icall
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
      integer me,nproc,i1,ii,ijm1,ijlast,jj,ij
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
         if(ispin.eq.1) then
c
c       Initialize the Pulay matrix.  Store -1 on the 1st row and column.
c       (1,1) element = 0.
c
c        if(me.eq.0)write(LUout,*)' Initializing DIIS; ispin = ',ispin

         call dfill((mxhist+1)*(mxhist+2)/2,0d0,pulay,1)
!        pulay(1) = 0.0d0
c
        do i = 2,mxhist+1
          i1=ijcanon(i,1)
          pulay(i1) = -1.0d0
        enddo 
      endif
c
        nhist(ispin)=0
      endif 
      ii=(ipol+ispin-1)*mxhist + irow(icall)
c
C     Store new errvec
c
      call ga_copy(g_err,g_diis(ii))
c
c     If more than mxhist histories are stored, shift the
c     rows of the Pulays B matrix.
c     
      if(nhist(1).eq.mxhist.and.ispin.eq.1)then
c
        do i = 3,nhist(1)+1
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
      nhist(ispin) = min(nhist(ispin)+1,mxhist)
c
c     Construct the last row of the B matrix.
c

      derr=ga_ddot(g_err,g_err)
      if(ispin.eq.1) then
         do jj = 1,nhist(1)
            ijlast = ijcanon(nhist(1)+1,nhist(1)+2-jj)
            pulay(ijlast)=0d0
         enddo
      endif
      do jj = 1,nhist(1)
        ii=(ipol+ispin-1)*mxhist+irow(icall-jj+1)
        ijlast = ijcanon(nhist(1)+1,nhist(1)+2-jj)
        pulay(ijlast)=pulay(ijlast)+ga_ddot(g_diis(ii),g_err)
      enddo
      return
C 1111 format(/,15x,'Diis error vector norm:   ',f20.10,/)
      end

