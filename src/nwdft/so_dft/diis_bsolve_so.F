      subroutine diis_bsolve_so(pulay,wght,mxhist,nhist)

C$Id$

c     
c***********************************************************************
c     
c     subroutine bsolve
c
c     Purpose:  Solve the following equation:
c           
c               B*W = C,
c
c     where B is Pulays B matrix constructed from the inner product
c     of error vectors; W is a weighting vector used to mix previous
c     guesses.  
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
c     Variables:
c                 wght(out):   weight vector, (W)
c                 c(in):       c(i) = 0 for i > 1, and c(1) = -1.
c                 mxhist(in):  maximum history.
c                 nhist(in):   current number of objects in history.
c                 pulay(in):   Pulays B matrix.
c                 maxhist(in): hardwired upper-bound to the maximum 
c                              allowed history.
c     
c***********************************************************************
c     
c     
      implicit none
      integer mxhist
      integer nhist
      double precision wght(mxhist+1),pulay((mxhist+1)*(mxhist+2)/2)
c     
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
c     
      logical LResult
      integer ijcanon,i,j
      integer ntemp,info,ij,ijm1
      integer lipiv,iipiv,liwork,iiwork,lwork,iwork,lcc,icc,
     ,     lafp,iafp
      double precision  rcond, ferr, berr
c

      ijcanon(i,j) = i*(i-1)/2 +j
      
 1000 if (nhist.lt.2) then 
        return
      endif 
c
      ntemp = (nhist+1)*(nhist+2)/2
      LResult = MA_Push_Get(MT_int,nhist+1,'Lapack 5',lipiv,iipiv)
      LResult = MA_Push_Get(MT_int,nhist+1,'Lapack 4',liwork,iiwork)
      LResult = MA_Push_Get(MT_Dbl,3*(nhist+1),'Lapack 3',lwork,iwork)
      LResult = MA_Push_Get(MT_Dbl,nhist+1,'Lapack 2',lcc,icc)
      LResult = MA_Push_Get(MT_Dbl, ntemp,'Lapack 1', lafp,iafp)
c
c     Define the c array.
c
      call dcopy(nhist+1,0d0,0,Dbl_MB(icc), 1)
      Dbl_MB(icc) = -1.0d0
c      if(ga_nodeid().eq.0)then
c        write(6,*)' pulay matrix '
c        call stpr_tprntd(pulay,nhist+1)
c        write(6,*)' rhs '
c        call output(Dbl_MB(icc),1,nhist+1,1,1,nhist+1,1,1)
c      endif
      call dspsvx( 'N', 'U', nhist+1, 1, pulay, Dbl_MB(iafp), 
     &             Int_MB(iipiv), Dbl_MB(icc), nhist+1, wght, nhist+1, 
     &             rcond, ferr, berr, Dbl_MB(iwork), Int_MB(iiwork), 
     &             info)
c      if(ga_nodeid().eq.0)then
c        write(6,*)' wght '
c        call output(wght,1,nhist+1,1,1,nhist+1,1,1)
c      endif
      LResult = MA_Pop_Stack(lafp)
      LResult = MA_Pop_Stack(lcc)
      LResult = MA_Pop_Stack(lwork)
      LResult = MA_Pop_Stack(liwork)
      LResult = MA_Pop_Stack(lipiv)
c
c     Check for singularity.  If linear dependency occurs, decrease nhist
c     and modify the Pulay matrix accordingly.
c     
      if(info.ne.0)then
        if(ga_nodeid().eq.0 .and. util_print('warning',print_default))
     &        write(LuOut,*)' Singularity in Pulay matrix. ',
     &       'Error and Fock matrices removed. '
        do i = 3,nhist+1
          do j = 3,i
c           
            ij = ijcanon(i,j)
            ijm1 = ijcanon(i-1,j-1)
            pulay(ijm1)=pulay(ij)
c            
          enddo 
        enddo 
        nhist=nhist-1
c
c Removed by work of Sherwood/Apra (8/1/96)
c        icall=icall-1
c
        goto 1000
      endif 
c
      return
 9000 format(10x,'DIIS skipped this iteration:',/,
     &       10x,'Number of error vectors in history ',
     &           'less than two.')
      end
