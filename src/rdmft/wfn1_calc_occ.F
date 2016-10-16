c-----------------------------------------------------------------------
c
      subroutine wfn1_calc_occ(g_c,nmo,ne,icnf,occ)
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
      double precision, intent(out) :: occ(nmo) ! the occupation numbers
c
      integer iproc ! the processor rank
      integer nproc ! the number of processor ranks
c
      integer ilo, ihi, jlo, jhi ! patch limits
      integer ii, jj             ! counters
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
      do jj = 1, ne
        if(icnf(jj).ge.jlo.and.icnf(jj).le.jhi)
     &  then
          do ii = ilo, ihi
            occ(ii) = occ(ii) + buf(ii,icnf(jj)) * buf(ii,icnf(jj))
          enddo
        endif
      enddo
      deallocate(buf)
      call ga_dgop(MT_DBL,occ,nmo,"+")
c
      end
c
c-----------------------------------------------------------------------
