      subroutine diis_builds12(
     *     isvals,iga_svecs,iga_sout,iga_stmp,nbf,iwhich)
c***********************************************************************
c     sbuild:  Build S, S**(-1/2) or S**(+1/2) from S eigenvalues and
c     S eigenvectors.
c     iwhich = 1; build S
c     iwhich = 2; build S**(-1/2)
c     iwhich = 3; build S**(+1/2)
c***********************************************************************
C$Id: diis_bld12.F,v 1.1 1995-12-13 01:35:52 d3g681 Exp $
      implicit real*8 (a-h,o-z) 
      logical LResult
      parameter (toll=1.d-6)
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c     
      me=ga_nodeid()
      nproc=ga_nnodes()
      call ga_sync
      LResult = MA_Push_Get(MT_Dbl,nbf,'tmpm',
     &     ltmpm, itmpm)

      call ga_copy(iga_svecs,iga_stmp)
      
      if(iwhich.eq.2)then
c     
c       Build S**(-1/2).
c     
        do i=me+1,nbf,nproc
          aaz = Dbl_MB(isvals+i-1)
          if (aaz.lt.toll)  then
            write(LuOut,*) ' singular eigenvalue in S matrix ',aaz
            aaz=0.d0
          else
            aaz=1.d0/dsqrt(aaz)
          endif
c          call ga_get(iga_svecs,
c     &         1,nbf,i,i,DBL_MB(itmpm),1)
          call get_col(iga_svecs,nbf,i,DBL_MB(itmpm))
          call dscal(nbf,aaz,DBL_MB(itmpm),1)
          call put_col(iga_stmp,nbf,i,DBL_MB(itmpm))
        enddo

      elseif(iwhich.eq.3)then
c     
c     Build S**(1/2).
c     
        do i=me+1,nbf,nproc
          aaz = Dbl_MB(isvals+i-1)
          if (aaz.lt.toll)  then
C            write(LuOut,*) ' singular eigenvalue in S matrix ',aaz
            aaz=0.d0
          else
            aaz=dsqrt(aaz)
          endif
          call get_col(iga_svecs,nbf,i,DBL_MB(itmpm))
          call dscal(nbf,aaz,DBL_MB(itmpm),1)
          call put_col(iga_stmp,nbf,i,DBL_MB(itmpm))
        enddo

      endif
      LResult = MA_Pop_Stack(ltmpm)

      call ga_sync

      call ga_dgemm('N','T',nbf,nbf,nbf,1.d0,
     &     iga_stmp,iga_svecs,0.d0,iga_sout)

      return
      end
