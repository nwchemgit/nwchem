      subroutine dft_canorg( nmo, svals, g_svecs, g_xij)
c
c     Canonical Orthogonalization
c     Szabo & Ostlund page 144 (Dover Ed.)
c     
c     $Id: dft_canorg.F,v 1.3 2001-10-05 22:27:12 edo Exp $
c     
      implicit none
c
      double precision svals(*) !  S evals [input]
      integer g_svecs ! GA handle for S evecs [input]
      integer g_xij  ! GA [ouput]
      integer nmo ! [out]
c
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
c     
      integer me, nproc, i
      integer ltmpm, itmpm,ndep
      double precision aaz
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      ndep=nbf-nmo
      
c
      if (.not. MA_Push_Get(MT_Dbl, nbf_ao, 'tmpm', ltmpm, itmpm))
     &   call errquit('diis_bld12: failed to alloc tmpm',0)
c     
c        Build Xij
c     
      call ga_zero(g_xij)
      do i = me+1, nmo, nproc
         aaz = svals(i+ndep)
         call get_col(g_svecs, nbf_ao, i+ndep, DBL_MB(itmpm))
         aaz = 1.d0/sqrt(aaz)
         call dscal(nbf_ao, aaz, DBL_MB(itmpm), 1)
         call put_col(g_xij, nbf_ao, i, DBL_MB(itmpm))
      enddo
      if(.not.MA_Pop_Stack(ltmpm)) call errquit(
     .     'failed popstack in dftcanorg ',0)
c     
      return
      end
      integer function dft_lindep(nbf,toll_s,svals,oprint,luout)
      implicit none
      integer nbf ! [in]
      double precision toll_s ! [in]
      double precision svals(*) ! [in]
      logical oprint
      integer luout
c
      double precision smax,smin
      integer i,nmo
c
c     get linear dep
c
      smax=-1d30
      smin=1d30
      nmo = nbf
      do i=1,nbf
         if(svals(i).lt.toll_s) then
            nmo=nmo-1
            smax=max(smax,svals(i))
            smin=min(smin,svals(i))
         endif
      enddo
      if(oprint.and.nbf.ne.nmo) then
         write(luout,*)
         write(luout,'(A,i5,A)') ' WARNING : Found ',nbf-nmo,
     /        ' linear dependencies '
         write(luout,11) ' S eigenvalue threshold:',toll_s
 11      format(A,1pe16.5)
         write(luout,11) ' Smallest S eigenvalue :',smin
         write(luout,11) ' Largest  S eigenvalue :',smax
         write(luout,*)
      endif
      dft_lindep=nmo
      return
      end
