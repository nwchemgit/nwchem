c-----------------------------------------------------------------------
c
      subroutine wfn1_2el_wfn1(wfn1_param,wfn1_inst,wfn1_wave,g_2aa)
      implicit none
C>
C> \brief Calculate the 2-electron density matrix according to
C> WFN1 theory (should be the as the one of Hartree-Fock theory
C> if the correlation functions are unit vectors) for debugging.
C>
C> The two electron energy is given by 
C> sum_ij (ii|jj) - (ij|ji)
C> where the integrals are given in natural orbital basis.
C>
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
c
      type(wfn1_prm),     intent(in) :: wfn1_param
      type(wfn1_prminst), intent(in) :: wfn1_inst
      type(wfn1_wfn),     intent(in) :: wfn1_wave
c
      integer, intent(inout) :: g_2aa
c
      double precision, allocatable :: d2el(:,:,:,:)
      double precision, allocatable :: bufw(:,:) ! wave function
      double precision, allocatable :: bufo(:,:) ! occupation numbers
c
      integer :: nbf
      integer :: nmo
      integer :: nea
      integer :: ii, jj
      integer :: ia, ib, ic, id
c
      double precision occ
c
      nbf = wfn1_param%nbf
      nmo = wfn1_param%nmo
      nea = wfn1_param%nea
      call wfn1_calc_occ2(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                    g_2aa)
      allocate(d2el(nbf,nbf,nbf,nbf))
      allocate(bufw(nbf,nmo),bufo(nmo,nmo))
      call ga_get(wfn1_wave%g_na,1,nbf,1,nmo,bufw,nbf)
      call ga_get(g_2aa,1,nmo,1,nmo,bufo,nmo)
      d2el = 0.0d0
c
c     Coulomb term
c
      do ii = 1, nmo
        do jj = 1, nmo
          occ = 0.5d0*bufo(ii,jj)
          do ia = 1, nbf
            do ib = 1, nbf
              do ic = 1, nbf
                do id = 1, nbf
                  d2el(ia,ib,ic,id) = d2el(ia,ib,ic,id)
     &                              + bufw(ia,ii)*bufw(ib,ii)*
     &                                bufw(ic,jj)*bufw(id,jj)*occ
     &                              + bufw(ia,jj)*bufw(ib,jj)*
     &                                bufw(ic,ii)*bufw(id,ii)*occ
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
c     exchange term
c
      do ii = 1, nmo
        do jj = 1, nmo
          occ = 0.5d0*bufo(ii,jj)
          do ia = 1, nbf
            do ib = 1, nbf
              do ic = 1, nbf
                do id = 1, nbf
                  d2el(ia,ib,ic,id) = d2el(ia,ib,ic,id)
     &                              - bufw(ia,ii)*bufw(ib,jj)*
     &                                bufw(ic,jj)*bufw(id,ii)*occ
     &                              - bufw(ia,jj)*bufw(ib,ii)*
     &                                bufw(ic,ii)*bufw(id,jj)*occ
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      write(*,*)'HvD: 2-electron density WFN1'
      write(*,'(6x,$)')
      do id = 1, nbf
        do ic = 1, nbf
          write(*,'(i12,$)')id
        enddo
      enddo
      write(*,*)
      write(*,'(6x,$)')
      do id = 1, nbf
        do ic = 1, nbf
          write(*,'(i12,$)')ic
        enddo
      enddo
      write(*,*)
      do ia = 1, nbf
        do ib = 1, nbf
          write(*,'(2i3,$)')ia,ib
          do id = 1, nbf
            do ic = 1, nbf
              write(*,'(f12.5,$)')d2el(ia,ib,ic,id)
            enddo
          enddo
          write(*,*)
        enddo
      enddo
c
      deallocate(bufo,bufw,d2el)
      end
