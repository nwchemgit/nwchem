c$Id: diis_hamwgt.F,v 1.2 1995-10-17 08:55:00 d3g681 Exp $
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
     &     iga_fock,iga_diis,iga_tmp)
      implicit real*8 (a-h,o-z) 
      dimension wght(mxhist+1)
c
      Integer nhist(2), icall(2)
c
#include "mafdecls.fh"
c
      irow(j) = mod(j-1,mxhist)+1
      ii=(ispin-1)*mxhist + irow(icall(ispin))
      istep=(ii-1)*nbf
      call ga_sync
c
#ifdef DATA_SERVER
        call ga_release(iga_fock,1,nbf,1,nbf)
        call ga_release(iga_diis,1,nbf,istep+1,istep+nbf)
#endif  
      call ga_copy_patch('N',
     *     iga_fock,1,nbf,1,nbf, 
     *     iga_diis,1,nbf,istep+1,istep+nbf)

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
#ifdef DATA_SERVER
        call ga_release(iga_tmp,1,nbf,1,nbf)
        call ga_release(iga_fock,1,nbf,1,nbf)
        call ga_release(iga_diis,1,nbf,istep+1,istep+nbf)
#endif  
        call ga_dadd_patch(
     *       xwght,iga_diis,1,nbf,istep+1,istep+nbf,
     *        anum,iga_fock,1,nbf,1,nbf,
     *            iga_tmp,1,nbf,1,nbf)
        call ga_copy(iga_tmp,iga_fock)
      enddo 
      return
      end


