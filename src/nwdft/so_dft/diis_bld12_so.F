      subroutine diis_bld12_so(toll_s, svals, iga_svecs, iga_sout, 
     &                      iga_stmp, nbf, iwhich)
c***********************************************************************
c     sbuild:  Build S, S**(-1/2) or S**(+1/2) from S eigenvalues and
c     S eigenvectors.
c     iwhich = 1; build S
c     iwhich = 2; build S**(-1/2)
c     iwhich = 3; build S**(+1/2)
c***********************************************************************
C$Id: diis_bld12_so.F,v 1.1 2000-04-17 21:45:57 zzhang Exp $
c
      implicit none
c
      double precision svals(*) !  S evals [input]
      integer iga_svecs ! GA handle for S evecs [input]
      integer iga_sout  ! GA handle for S^(n) [ouput]
      integer iga_stmp  ! GA handle for scratch [input]
      integer nbf       ! no. basis fns [input]
      integer iwhich    ! 2->S(-1/2) 3->S(1/2) 
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
     &   call errquit('diis_bld12: failed to alloc tmpm',0)
c
      call ga_copy(iga_svecs, iga_stmp)
c      
      if (iwhich.eq.1)then
c     
c        Build S.
c     
         call ga_zero(iga_stmp)
         do i = me+1, nbf, nproc
            aaz = svals(i)
            if (aaz.ge.toll_s)  then
               call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
               call dscal(nbf, aaz, DBL_MB(itmpm), 1)
               call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
            else
c               write(LuOut,*)' Small eig of S in build of S',
c     &           i, aaz, toll_s
               call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
               call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
            endif
         enddo
c
      elseif (iwhich.eq.2)then
c     
c        Build S**(-1/2).
c     
         call ga_zero(iga_stmp)
         do i = me+1, nbf, nproc
            aaz = svals(i)
            if (aaz.ge.toll_s)  then
               aaz = 1.d0/sqrt(aaz)
               call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
               call dscal(nbf, aaz, DBL_MB(itmpm), 1)
               call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
            else
c               write(LuOut,*)' Small eig of S in build of S^-1/2',
c     &           i, aaz, toll_s
               call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
               call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
            endif
         enddo
c
      elseif (iwhich.eq.3)then
c     
c        Build S**(1/2).
c     
         call ga_zero(iga_stmp)
         do i = me+1, nbf, nproc
            aaz = svals(i)
            if (aaz.ge.toll_s)  then
               aaz = dsqrt(aaz)
               call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
               call dscal(nbf, aaz, DBL_MB(itmpm), 1)
               call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
            else
c               write(LuOut,*)' Small eig of S in build of S^1/2',
c     &           i, aaz, toll_s
               call get_col(iga_svecs, nbf, i, DBL_MB(itmpm))
               call put_col(iga_stmp, nbf, i, DBL_MB(itmpm))
            endif
         enddo
      endif
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
