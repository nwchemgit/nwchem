      subroutine diis_builds12(
     *     isvals,iga_svecs,iga_sout,iga_stmp,nbf,iwhich)
c***********************************************************************
c     sbuild:  Build S, S**(-1/2) or S**(+1/2) from S eigenvalues and
c     S eigenvectors.
c     iwhich = 1; build S
c     iwhich = 2; build S**(-1/2)
c     iwhich = 3; build S**(+1/2)
c***********************************************************************
C$Id: diis_bld12.F,v 1.3 1996-04-16 22:09:44 og845 Exp $
      implicit none
      integer isvals    ! MA dbl handle for S evals [input]
      integer iga_svecs ! GA handle for S evecs [input]
      integer iga_sout  ! GA handle for S^(n) [ouput]
      integer iga_stmp  ! GA handle for scratch [input]
      integer nbf       ! no. basis fns [input]
      integer iwhich    ! 2->S(-1/2) 3->S(1/2) 
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c     
      logical LResult
      integer me,nproc,i
      integer ltmpm,itmpm
      double precision toll,aaz
      parameter (toll=1.d-5)
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
            aaz=1.d0/sqrt(aaz)
          endif
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
