      subroutine diis_bld12_so(toll_s, svals, iga_svecs, iga_sout, 
     &                      iga_stmp, nbf, iwhich)
c***********************************************************************
c     sbuild:  Build S, S**(-1/2) or S**(+1/2) from S eigenvalues and
c     S eigenvectors.
c     iwhich = 1; build S
c     iwhich = 2; build S**(-1/2)
c     iwhich = 3; build S**(+1/2)
c     iwhich = 4; build S**(-1)
c***********************************************************************
C$Id$
c
      implicit none
#include "errquit.fh"
c
      double precision svals(*) !  S evals [input]
      integer iga_svecs ! GA handle for S evecs [input]
      integer iga_sout  ! GA handle for S^(n) [ouput]
      integer iga_stmp  ! GA handle for scratch [input]
      integer nbf       ! no. basis fns [input]
      integer iwhich    ! 1:S, 2:S(-1/2), 3:S(1/2), 4:S(-1) 
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "rtdb.fh"
c     
      logical LResult
      integer me, nproc, i
      integer ltmpm, itmpm
      double precision toll_s, aaz
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      call ga_sync
c
      if (.not. MA_Push_Get(MT_Dbl, nbf, 'tmpm', ltmpm, itmpm))
     &   call errquit('diis_bld12_so: failed to alloc tmpm',0, MA_ERR)
c
      call ga_copy(iga_svecs, iga_stmp)
c      
c     Build S,S^(-1/2),S^(1/2),S^(-1)
c     
      call ga_zero(iga_stmp)
      do i = me+1, nbf, nproc
         aaz = svals(i)
         if (aaz.ge.toll_s)  then
c
            if (iwhich.eq.1) aaz = 1.d0*aaz      ! S
            if (iwhich.eq.2) aaz = 1.d0/sqrt(aaz)! S^(-1/2)
            if (iwhich.eq.3) aaz = dsqrt(aaz)    ! S^(1/2)
            if (iwhich.eq.4) aaz = 1.d0/aaz      ! S^(-1)
c
            call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
            call dscal(nbf, aaz, DBL_MB(itmpm), 1)
            call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
         else
            call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
            call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
         endif
      enddo
c
      LResult = MA_Pop_Stack(ltmpm)
c
      call ga_sync
c
      call ga_dgemm('N', 'T', nbf, nbf, nbf, 1.d0, 
     &              iga_stmp, iga_svecs, 0.d0, iga_sout)
c
      return
      end
