c-----------------------------------------------------------------------
c
      subroutine wfn1_ortho(wfn1_mat,wfn1_wave)
      implicit none
#include "wfn1_matcnst.fh"
#include "wfn1_wfn.fh"
C>
C>    \brief Orthonormalize the vectors 
C>
C>    This subroutine ortho-normalizes all vectors making up the 
C>    wavefunction, i.e. the natural orbitals as well as the correlation
C>    functions. 
C>
      type(wfn1_matcnst), intent(in) :: wfn1_mat ! matrix constants
      type(wfn1_wfn),  intent(inout) :: wfn1_wave ! the wave function
c
      integer :: idum  ! dummy variable
c
      idum = 0
      call ga_orthog(wfn1_wave%g_ca,idum,.false.)
      call ga_orthog(wfn1_wave%g_cb,idum,.false.)
      call ga_orthog(wfn1_wave%g_na,wfn1_mat%g_s,.true.)
      call ga_orthog(wfn1_wave%g_nb,wfn1_mat%g_s,.true.)
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_tst_ortho(wfn1_mat,wfn1_wave,g_bfmo,g_momo)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_wfn.fh"
C>
C>    \brief Test routine for wfn1_ortho
C>
C>    This routine simply calls wfn1_ortho but then checks the results
C>    to make sure that the orbitals returned are actually orthonormal.
C>
      type(wfn1_matcnst), intent(in) :: wfn1_mat ! matrix constants
      type(wfn1_wfn),  intent(inout) :: wfn1_wave ! the wave function
c
      integer, intent(inout) :: g_bfmo ! temporary matrix nbf x nmo
      integer, intent(inout) :: g_momo ! temporary matrix nmo x nmo
c
      integer :: nbf   ! the number of basis functions
      integer :: nmo   ! the number of MOs
      integer :: itype ! the GA type
      integer :: iproc ! the current processor rank
      integer :: ii    ! counter
c
      integer :: ilo, ihi ! local matrix block limits
      integer :: jlo, jhi ! local matrix block limits
c
      double precision :: dnormT  ! norm of temporary matrix
c
      call wfn1_ortho(wfn1_mat,wfn1_wave)
c
      iproc = ga_nodeid()
      call ga_inquire(wfn1_wave%g_na,itype,nbf,nmo)
c
c     Check alpha natural orbitals Na^T S Na = I
c
      call ga_dgemm('n','n',nbf,nmo,nbf,1.0d0,wfn1_mat%g_s,
     &              wfn1_wave%g_na,0.0d0,g_bfmo)
      call ga_dgemm('t','n',nmo,nmo,nbf,1.0d0,wfn1_wave%g_na,
     &              g_bfmo,0.0d0,g_momo)
      call ga_distribution(g_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_momo,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(g_momo,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Na: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(wfn1_wave%g_na)
c       call ga_print(g_momo)
      endif
c
c     Check beta natural orbitals Nb^T S Nb = I
c
      call ga_dgemm('n','n',nbf,nmo,nbf,1.0d0,wfn1_mat%g_s,
     &              wfn1_wave%g_nb,0.0d0,g_bfmo)
      call ga_dgemm('t','n',nmo,nmo,nbf,1.0d0,wfn1_wave%g_nb,
     &              g_bfmo,0.0d0,g_momo)
      call ga_distribution(g_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_momo,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(g_momo,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Nb: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(wfn1_wave%g_nb)
c       call ga_print(g_momo)
      endif
c
c     Check alpha correlation functions Ca^T Ca = I
c
      call ga_dgemm('t','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,
     &              wfn1_wave%g_ca,0.0d0,g_momo)
      call ga_distribution(g_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_momo,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(g_momo,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Ca: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(wfn1_wave%g_ca)
c       call ga_print(g_momo)
      endif
c
c     Check beta correlation functions Cb^T Cb = I
c
      call ga_dgemm('t','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,
     &              wfn1_wave%g_cb,0.0d0,g_momo)
      call ga_distribution(g_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(g_momo,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(g_momo,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Cb: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(wfn1_wave%g_cb)
c       call ga_print(g_momo)
      endif
      end

c
c-----------------------------------------------------------------------
