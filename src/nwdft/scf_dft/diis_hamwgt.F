c$Id: diis_hamwgt.F,v 1.5 2000-03-16 00:06:27 edo Exp $
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
      subroutine diis_hamwgt(wght,mxhist,nhist,nbf,icall,ispin,ipol,
     &     g_fock,g_diis,g_tmp)
      implicit none
c
      integer mxhist
      integer nbf
      integer ispin
      integer ipol
      integer g_fock
      integer g_diis(mxhist)
      integer g_tmp
      double precision wght(mxhist+1)
c
      Integer nhist(2), icall(2)
c
#include "mafdecls.fh"
c
      integer ii,istep,jj
      double precision xwght,anum
      integer irow,j
c
      irow(j) = mod(j-1,mxhist)+1
      ii=(ispin-1)*mxhist + irow(icall(ispin))
      istep=(ii-1)*nbf
cedo      call ga_sync
c
cedo      call ga_copy_patch('N',
cedo     *     g_fock,1,nbf,1,nbf, 
cedo     *     g_diis,1,nbf,istep+1,istep+nbf)
      call ga_copy(g_fock,g_diis(ii))

c     
      if (nhist(ispin).eq.1) then
        return
      endif 
c
c     Construct new Hamiltonian from weighted sum of previous Hamiltonians.
c
      do jj = 1,nhist(ispin)
        xwght = wght(nhist(ispin)-jj+2)
        ii=(ispin-1)*mxhist + irow(icall(ispin)-jj+1)
        istep=(ii-1)*nbf
        anum=1.d0
        if(jj.eq.1) anum=0.d0
cedo        call ga_dadd_patch(
cedo     *       xwght,g_diis,1,nbf,istep+1,istep+nbf,
cedo     *        anum,g_fock,1,nbf,1,nbf,
cedo     *            g_tmp,1,nbf,1,nbf)
        call ga_dadd(xwght,g_diis(ii),anum,g_fock,g_tmp)
        call ga_copy(g_tmp,g_fock)
      enddo 
      return
      end


