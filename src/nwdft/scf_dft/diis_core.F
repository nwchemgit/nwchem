      subroutine diis_core(derr,pulay,
     &     mxhist,nhist,icall,ispin,ipol, g_diis,g_err,
     &     os_constr)
C$Id: diis_core.F,v 1.12 2003-03-12 02:54:36 edo Exp $
      implicit none
c
      integer mxhist
      integer ispin
      integer ipol
      integer g_diis(mxhist)
      integer g_err
      double precision pulay((mxhist+1)*(mxhist+2)/2,*)
      Integer nhist, icall(2)
      logical os_constr
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
      integer iptr
c
c
      irow(j) = mod(j-1,mxhist)+1
      ijcanon(i,j) = i*(i-1)/2 +j
      oprint = util_print('convergence', print_default)
      me=ga_nodeid()
      nproc=ga_nnodes()
c
      iptr=1
      if(os_constr) iptr=ispin
c
      if(icall(iptr).eq.1) then
         if(ispin.eq.1.or.os_constr) then
c
c       Initialize the Pulay matrix.  Store -1 on the 1st row and column.
c       (1,1) element = 0.
c

         call dfill((mxhist+1)*(mxhist+2)/2,0d0,pulay(1,iptr),1)
c
        do i = 2,mxhist+1
          i1=ijcanon(i,1)
          pulay(i1,iptr) = -1.0d0
        enddo 
      endif
c
        nhist=0
      endif 
      ii=(ipol+ispin-1)*mxhist + irow(icall(iptr))
c
C     Store new errvec
c
      call ga_copy(g_err,g_diis(ii))
c
c     If more than mxhist histories are stored, shift the
c     rows of the Pulays B matrix.
c     
      if(nhist.eq.mxhist.and.(ispin.eq.1.or.os_constr))then
c
        do i = 3,nhist+1
          do j = 3,i
c
            ij = ijcanon(i,j)
            ijm1 = ijcanon(i-1,j-1)
            pulay(ijm1,iptr)=pulay(ij,iptr)
          enddo 
        enddo 
c       
      endif 
c
c     Determine the size of the iterative subspace history.
c
      nhist = min(nhist+1,mxhist)
c
c     Construct the last row of the B matrix.
c

      derr=ga_ddot(g_err,g_err)
      if(ispin.eq.1.or.os_constr) then
         do jj = 1,nhist
            ijlast = ijcanon(nhist+1,nhist+2-jj)
            pulay(ijlast,iptr)=0d0
         enddo
      endif
      do jj = 1,nhist
        ii=(ipol+ispin-1)*mxhist+irow(icall(iptr)-jj+1)
        ijlast = ijcanon(nhist+1,nhist+2-jj)
        pulay(ijlast,iptr)=pulay(ijlast,iptr)+
     +       ga_ddot(g_diis(ii),g_err)
      enddo
      return
C 1111 format(/,15x,'Diis error vector norm:   ',f20.10,/)
      end

