c-----------------------------------------------------------------------
c
      subroutine wfn1_2el_hf(wfn1_param,wfn1_inst,wfn1_wave)
      implicit none
C>
C> \brief Calculate the 2-electron density matrix according to
C> Hartree-Fock theory for debugging
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
      double precision, allocatable :: d2el(:,:,:,:)
      double precision, allocatable :: buf(:,:)
c
      integer :: nbf
      integer :: ii, jj
      integer :: ia, ib, ic, id
c
      nbf = wfn1_param%nbf
      allocate(d2el(nbf,nbf,nbf,nbf))
      allocate(buf(nbf,nbf))
      call ga_get(wfn1_wave%g_na,1,nbf,1,nbf,buf,nbf)
      d2el = 0.0d0
c
c     Coulomb term
c
      do ii = 1, nbf
        do jj = 1, nbf
          do ia = 1, nbf
            do ib = 1, nbf
              do ic = 1, nbf
                do id = 1, nbf
                  d2el(ia,ib,ic,id) = d2el(ia,ib,ic,id)
     &                              + buf(ia,ii)*buf(ib,ii)*
     &                                buf(ic,jj)*buf(id,jj)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      do ii = 1, nbf
        do jj = 1, nbf
          do ia = 1, nbf
            do ib = 1, nbf
              do ic = 1, nbf
                do id = 1, nbf
                  d2el(ia,ib,ic,id) = d2el(ia,ib,ic,id)
     &                              - buf(ia,ii)*buf(ib,jj)*
     &                                buf(ic,jj)*buf(id,ii)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c
      write(*,*)'HvD: 2-electron density Hartree-Fock'
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
      deallocate(buf,d2el)
      end
