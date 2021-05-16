      subroutine diis_bsolve(pulay,wght,mxhist,nhist)

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
      integer lipiv,iipiv,liwork,iiwork,lwork,iwork,l_cc,i_cc,
     ,     lafp,iafp
      double precision  rcond, ferr, berr
      logical oprint
c
      ijcanon(i,j) = i*(i-1)/2 +j
      oprint= util_print('diis_debug',print_high).and.
     A     ga_nodeid().eq.0
      
 1000 continue
      if(oprint)then
        write(6,*)'nhist ',nhist,' pulay matrix '
        do i = 1,nhist+1
           write(6,123) (pulay(ijcanon(i,j)),j=1,i)
 123       format(10(1x,1pe12.4))
        enddo 
      endif
      if (nhist.lt.2) return
c
      ntemp = (nhist+1)*(nhist+2)/2
      LResult = MA_Push_Get(MT_int,nhist+1,'Lapack 5',lipiv,iipiv)
      LResult = MA_Push_Get(MT_int,nhist+1,'Lapack 4',liwork,iiwork)
      LResult = MA_Push_Get(MT_Dbl,3*(nhist+1),'Lapack 3',lwork,iwork)
      LResult = MA_Push_Get(MT_Dbl,nhist+1,'Lapack 2',l_cc,i_cc)
      LResult = MA_Push_Get(MT_Dbl, ntemp,'Lapack 1', lafp,iafp)
c
c     Define the c array.
c
      call dcopy(nhist+1,0d0,0,Dbl_MB(i_cc), 1)
      Dbl_MB(i_cc) = -1.0d0
      if(oprint)then
        write(6,*)' rhs '
        write(6,123) (Dbl_MB(i_cc+j-1),j=1,nhist+1)
      endif
      info=0
      call dspsvx( 'N', 'U', nhist+1, 1, pulay, Dbl_MB(iafp), 
     &             Int_MB(iipiv), Dbl_MB(i_cc), nhist+1, wght, nhist+1, 
     &             rcond, ferr, berr, Dbl_MB(iwork), Int_MB(iiwork), 
     &             info)
      if(oprint)then
        write(6,*)' wght '
        call output(wght,1,nhist+1,1,1,nhist+1,1,1)
      endif
      LResult = MA_chop_Stack(lipiv)
c
c     Check for singularity.  If linear dependency occurs, decrease nhist
c     and modify the Pulay matrix accordingly.
c     
      if(info.ne.0)then
c        write(6,*)' info ',info
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
