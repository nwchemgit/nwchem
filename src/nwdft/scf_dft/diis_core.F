      subroutine diis_core(pulay,mxhist,nhist,icall,nbf,ispin,ipol,
     &     iga_diis,iga_err)
C$Id: diis_core.F,v 1.1 1995-04-27 03:16:22 og845 Exp $
      implicit real*8 (a-h,o-z) 
c
      dimension pulay((mxhist+1)*(mxhist+2)/2,ipol)
      Integer nhist(2), icall(2)
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
c
#include "mafdecls.h"
#include "stdio.fh"
#include <global.fh>
#include <tcgmsg.fh>
c
      irow(j) = mod(j-1,mxhist)+1
      ijcanon(i,j) = i*(i-1)/2 +j
      me=ga_nodeid()
      nproc=ga_nnodes()
c
      if(icall(ispin).eq.1) then
c
c       Initialize the Pulay matrix.  Store -1 on the 1st row and column.
c       (1,1) element = 0.
c
        pulay(1,ispin) = 0.0d0
c
        do i = 2,mxhist+1
          i1=ijcanon(i,1)
          pulay(i1,ispin) = -1.0d0
        enddo 
c
        nhist(ispin)=0
      endif 
      ii=(ipol+ispin-1)*mxhist + irow(icall(ispin))
      istep=(ii-1)*nbf
      call ga_sync
c
C     Store new errvec
c

      call ga_copy_patch('N',
     *     iga_err,1,nbf,1,nbf, 
     *     iga_diis,1,nbf,istep+1,istep+nbf)
c
c     If more than mxhist histories are stored, shift the
c     rows of the Pulays B matrix.
c     
      if(nhist(ispin).eq.mxhist)then
c
        do i = 3,nhist(ispin)+1
          do j = 3,i
c
            ij = ijcanon(i,j)
            ijm1 = ijcanon(i-1,j-1)
            pulay(ijm1,ispin)=pulay(ij,ispin)
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

      temp=ga_ddot(iga_err,iga_err)
      if(me.eq.0)write(LuOut,1111)temp
      do jj = 1,nhist(ispin)
c
        ii=(ipol+ispin-1)*mxhist+irow(icall(ispin)-jj+1)
        istep=(ii-1)*nbf
        call ga_sync

        ijlast = ijcanon(nhist(ispin)+1,nhist(ispin)+2-jj)
        pulay(ijlast,ispin)=ga_ddot_patch(
     *       iga_diis,'N',1,nbf,istep+1,istep+nbf,
     *       iga_err,'N',1,nbf,1,nbf)
c        
      enddo
      return
 1111 format(/,15x,'Diis error vector norm:   ',f20.10,/)
      end

