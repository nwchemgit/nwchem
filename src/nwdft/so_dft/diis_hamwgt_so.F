c$Id$
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
      subroutine diis_hamwgt_so(wght,mxhist,nhist,nbf,icall,
     &     g_fock,g_diis,g_tmp)
      implicit none
c
      integer mxhist
      integer nbf
      integer g_fock(2)
      integer g_diis(2)
      integer g_tmp(2)
      double precision wght(mxhist+1)
c
      Integer nhist, icall
c
#include "mafdecls.fh"
c
      integer ii,istep,jj
      double precision anum
      double precision xwght
      integer irow,j
c
      irow(j) = mod(j-1,mxhist)+1
      ii=irow(icall)
      istep=(ii-1)*nbf
      call ga_sync
c
      call ga_copy_patch('N',
     *     g_fock(1),1,nbf,1,nbf, 
     *     g_diis(1),1,nbf,istep+1,istep+nbf)
      call ga_copy_patch('N',
     *     g_fock(2),1,nbf,1,nbf, 
     *     g_diis(2),1,nbf,istep+1,istep+nbf)
c     
      if (nhist.eq.1) then
        return
      endif 
c     
c     Construct new Hamiltonian from weighted sum of previous Hamiltonians.
c     
      do jj = 1,nhist
         xwght = wght(nhist-jj+2)
         ii=irow(icall-jj+1)
         istep=(ii-1)*nbf
         anum=1.d0
         if(jj.eq.1) anum=0.d0
c         write(*,*)"jj=", jj, xwght
c         call ga_print_patch(g_diis(1),1,nbf,istep+1,istep+nbf,1)
         call ga_dadd_patch(
     *        xwght,g_diis(1),1,nbf,istep+1,istep+nbf,
     *        anum,g_fock(1),1,nbf,1,nbf,
     *        g_tmp(1),1,nbf,1,nbf)
         call ga_copy(g_tmp(1),g_fock(1))
         call ga_dadd_patch(
     *        xwght,g_diis(2),1,nbf,istep+1,istep+nbf,
     *        anum,g_fock(2),1,nbf,1,nbf,
     *        g_tmp(2),1,nbf,1,nbf)
         call ga_copy(g_tmp(2),g_fock(2))
      enddo 
      return
      end




