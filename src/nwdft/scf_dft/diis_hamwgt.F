c$Id: diis_hamwgt.F,v 1.9 2000-11-14 20:26:09 edo Exp $
c***********************************************************************
c
c     subroutine hamwght
c
c     Construct new Hamiltonian from weighted average of Hamiltonians
c     constructed in previous iterations.  Used in conjunction with the
c     scferrv andd diis routines.
c
c     H  = w(n)*H(n) + w(n-1)*H(n-1) + w(n-2)*H(n-2) + . . .
c
c
c***********************************************************************
c
      subroutine diis_hamwgt(wght,mxhist,nhist,icall,ipol,
     &     g_focks,g_diis)
      implicit none
c
      integer mxhist
      integer ipol
      integer ispin
      integer g_focks(2)
      integer g_diis(mxhist)
      double precision wght(mxhist+1)
c
      Integer nhist(2), icall(2)
c
#include "mafdecls.fh"
c
      integer ii,jj
      double precision xwght,anum
      integer irow,j
c
      irow(j) = mod(j-1,mxhist)+1
      do ispin=1,ipol
         ii=(ispin-1)*mxhist + irow(icall(1))
         call ga_copy(g_focks(ispin),g_diis(ii))

c     
c
c     Construct new Hamiltonian from weighted sum of previous Hamiltonians.
c
         if(nhist(1).gt.1) then
            do jj = 1,nhist(1)
               xwght = wght(nhist(1)-jj+2)
               ii=(ispin-1)*mxhist + irow(icall(1)-jj+1)
               anum=1.d0
               if(jj.eq.1) anum=0.d0
               call ga_dadd(xwght,g_diis(ii),anum,
     .              g_focks(ispin),g_focks(ispin))
            enddo 
         endif
      enddo
      return
      end


