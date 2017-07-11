c-----------------------------------------------------------------------
c
      subroutine wfn1_calc_occ(g_c,nmo,ne,icnf,dcnt,occ)
      implicit none
C>
C> \brief Calculate occupation numbers from a set of correlation
C> functions
C>
#include "global.fh"
#include "mafdecls.fh"
c
      integer, intent(in) :: g_c      ! the correlation functions
      integer, intent(in) :: nmo      ! the number of molecular orbitals
      integer, intent(in) :: ne       ! the number of electrons
      integer, intent(in) :: icnf(nmo)! the electron configuration
c
      double precision, intent(in)  :: dcnt(nmo)! the electron fraction
c
      double precision, intent(out) :: occ(nmo) ! the occupation numbers
c
      integer iproc ! the processor rank
      integer nproc ! the number of processor ranks
c
      integer ilo, ihi, jlo, jhi ! patch limits
      integer ii, ir             ! counters
c
      double precision, allocatable :: buf(:,:)
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      occ = 0.0d0
      call ga_distribution(g_c,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      call ga_get(g_c,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      do ir = 1, ne
        if(icnf(ir).ge.jlo.and.icnf(ir).le.jhi)
     &  then
          do ii = ilo, ihi
            occ(ii) = occ(ii) + buf(ii,icnf(ir)) * buf(ii,icnf(ir)) *
     &                          dcnt(ir)
          enddo
        endif
      enddo
      deallocate(buf)
      call ga_dgop(MT_DBL,occ,nmo,"+")
c
      end subroutine wfn1_calc_occ
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_calc_occ2(g_c,nmo,ne,icnf,dcnt,g_occ)
      implicit none
C>
C> \brief Calculate the 2-electron occupation numbers from a set of
C> correlation functions
C>
#include "global.fh"
#include "mafdecls.fh"
c
      integer, intent(in) :: g_c       ! the correlation functions
      integer, intent(in) :: nmo       ! the number of molecular orbitals
      integer, intent(in) :: ne        ! the number of electrons
      integer, intent(in) :: icnf(1:ne)! the electron configuration
      integer, intent(inout) :: g_occ  ! the 2-electron occupation 
                                       ! numbers
c
      double precision, intent(in) :: dcnt(1:ne)! the electron fraction
c
      integer :: iolo, iohi, jolo, johi ! block limits of occupation GA
c
      double precision, allocatable :: bocc(:,:)  ! buffer for 
                                                  ! occupation numbers
      double precision, allocatable :: bcori(:,:) ! buffer for 
                                                  ! correlation function
      double precision, allocatable :: bcorj(:,:) ! buffer for 
                                                  ! correlation function
c
      integer iproc ! the current processor rank
      integer nproc ! the number of processor ranks
c
      integer is, it ! counters for the correlation functions
      integer ii, ij ! counters for the natural orbitals
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      call ga_distribution(g_occ,iproc,iolo,iohi,jolo,johi)
      allocate(bocc(iolo:iohi,jolo:johi))
      allocate(bcori(iolo:iohi,1:nmo))
      allocate(bcorj(jolo:johi,1:nmo))
      bocc = 0.0d0
      call ga_get(g_c,iolo,iohi,1,nmo,bcori,iohi-iolo+1)
      call ga_get(g_c,jolo,johi,1,nmo,bcorj,johi-jolo+1)
      do ij = jolo, johi
        do ii = iolo, iohi
          do it = 1, ne
            do is = 1, ne
              bocc(ii,ij) = bocc(ii,ij) 
     +                    + bcori(ii,icnf(is))*bcori(ii,icnf(is))
     +                    * bcorj(ij,icnf(it))*bcorj(ij,icnf(it))
     +                    * dcnt(is)*dcnt(it)
     +                    - bcori(ii,icnf(is))*bcori(ii,icnf(it))
     +                    * bcorj(ij,icnf(it))*bcorj(ij,icnf(is))
     +                    * dcnt(is)*dcnt(it)
            enddo
          enddo
        enddo
      enddo
      call ga_put(g_occ,iolo,iohi,jolo,johi,bocc,iohi-iolo+1)
c
      deallocate(bcorj)
      deallocate(bcori)
      deallocate(bocc)
c
      end subroutine wfn1_calc_occ2
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_calc_docx2(g_c,nmo,ne,icnf,dcnt,g_ocx)
      implicit none
C>
C> \brief Calculate the derivative of the 2-electron exchange 
C> occupation numbers
C>
C> These quantities are the sum_s Cjs Cks terms in the derivative of 
C> the 2-electron density matrix.
C>
#include "global.fh"
#include "mafdecls.fh"
c
      integer, intent(in) :: g_c       ! the correlation functions
      integer, intent(in) :: nmo       ! the number of molecular orbitals
      integer, intent(in) :: ne        ! the number of electrons
      integer, intent(in) :: icnf(1:ne)! the electron configuration
c
      double precision, intent(in) :: dcnt(1:ne)! the electron fraction
c
      integer, intent(inout) :: g_ocx  ! the 2-electron occupation 
                                       ! numbers
c
      integer :: iolo, iohi, jolo, johi ! block limits of occupation GA
c
      double precision, allocatable :: bocc(:,:)  ! buffer for 
                                                  ! occupation numbers
      double precision, allocatable :: bcori(:,:) ! buffer for 
                                                  ! correlation function
      double precision, allocatable :: bcorj(:,:) ! buffer for 
                                                  ! correlation function
c
      integer iproc ! the current processor rank
      integer nproc ! the number of processor ranks
c
      integer is     ! counters for the correlation functions
      integer ii, ij ! counters for the natural orbitals
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      call ga_distribution(g_ocx,iproc,iolo,iohi,jolo,johi)
      allocate(bocc(iolo:iohi,jolo:johi))
      allocate(bcori(iolo:iohi,1:nmo))
      allocate(bcorj(jolo:johi,1:nmo))
      bocc = 0.0d0
      call ga_get(g_c,iolo,iohi,1,nmo,bcori,iohi-iolo+1)
      call ga_get(g_c,jolo,johi,1,nmo,bcorj,johi-jolo+1)
      do ij = jolo, johi
        do ii = iolo, iohi
          do is = 1, ne
            bocc(ii,ij) = bocc(ii,ij) 
     +                    + bcori(ii,icnf(is))*bcorj(ij,icnf(is))
     +                    * dcnt(is)
          enddo
        enddo
      enddo
      call ga_put(g_ocx,iolo,iohi,jolo,johi,bocc,iohi-iolo+1)
c
      deallocate(bcorj)
      deallocate(bcori)
      deallocate(bocc)
c
      end subroutine wfn1_calc_docx2
c
c-----------------------------------------------------------------------
