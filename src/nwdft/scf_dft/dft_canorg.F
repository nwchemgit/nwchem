      subroutine dft_canorg( nmo, svals, g_svecs, g_xij)
c
c     Canonical Orthogonalization
c     Szabo & Ostlund page 144 (Dover Ed.)
c     
c     $Id: dft_canorg.F,v 1.4 2001-11-01 01:26:25 edo Exp $
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
      subroutine dft_guessin(movecs_in,ldmix,ncanorg,fon,
     ,     vecs_or_dens,
     ,     ipol,nbf_ao,g_movecs,g_gmovecs,
     ,     toll_s,svals)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
      character*(*) movecs_in
      logical ldmix ! [out]
      integer ncanorg ! [out]
      logical fon ! [in]
      character*(*) vecs_or_dens ! [out]
      integer ipol,g_movecs(2),g_gmovecs(2),nmo(2)
      integer nbf_ao
      double precision svals(*),toll_s
c
      external dft_lindep
      integer dft_lindep,nmo_out,me
c
      me=ga_nodeid()
      vecs_or_dens='density'
      ldmix=.false.
      ncanorg=0
      if((movecs_in.ne.'atomic'.and.
     *     movecs_in.ne.'fragment'.and.movecs_in.ne.'hcore')
     .     .or.fon) then
         vecs_or_dens='vectors'
         ldmix=.true.
         ncanorg=1
      else
      endif
      nmo_out=dft_lindep(nbf_ao,toll_s,
     .        svals,me.eq.0,luout)
      if(nmo_out.ne.nbf_ao) then
         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'gvecs', 
     &        nbf_ao, 0, g_gmovecs(1))) 
     .        call errquit(' dftscf: ga_create failed ',0)
         if(ipol.eq.2) then
            if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'gvecs', 
     &           nbf_ao, 0, g_gmovecs(2))) 
     .           call errquit(' dftscf: ga_create failed ',0)
         endif
         call ga_copy(g_movecs(1),g_gmovecs(1))
         if(ipol.eq.2)
     /        call ga_copy(g_movecs(2),g_gmovecs(2))
      else
         g_gmovecs(1)=g_movecs(1)
         g_gmovecs(2)=g_movecs(2)
      endif
      if(movecs_in.eq.'atomic') then
         nmo(1)=nmo_out
         nmo(2)=nmo(1)
      endif
      return
      end
      subroutine dft_guessout(nmo,nbf_ao,g_gmovecs,g_movecs,ipol)
      implicit none
#include "global.fh"
      integer nmo(2),nbf_ao,g_gmovecs(2),g_movecs(2),ipol
c
      if(nmo(1).ne.nbf_ao) then
         call ga_zero(g_movecs(1))
            write(0,*) ' nbf nm1 ',nbf_ao,nmo(1)
         call ga_copy_patch('n', g_gmovecs(1), 1, nbf_ao, 
     .        1, nmo(1),
     ,        g_movecs(1), 1,nbf_ao,1,nmo(1))
      if (.not. ga_destroy(g_gmovecs(1))) call errquit
     &     ('dft_scf: could not destroy g_denso', 0)
         if(ipol.eq.2) then
            write(0,*) ' nbf nm2 ',nbf_ao,nmo(2)
            nmo(2)=nmo(1)
            call ga_zero(g_movecs(2))
            call ga_copy_patch('n', g_gmovecs(2), 1, nbf_ao, 
     .           1, nmo(2),
     ,           g_movecs(2), 1,nbf_ao,1,nmo(2))
         if (.not. ga_destroy(g_gmovecs(2))) call errquit
     &        ('dft_scf: could not destroy g_denso', 0)
         endif
      endif
      return
      end
