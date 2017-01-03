c-----------------------------------------------------------------------
c
      subroutine wfn1_print_occ(nmo,occa,occb)
      implicit none
C>
C> \brief Print the occupation numbers
C>
#include "stdio.fh"
#include "global.fh"
c
      integer, intent(in) :: nmo ! the number of orbitals
c
      double precision, intent(in) :: occa(nmo) ! a-occupation numbers
      double precision, intent(in) :: occb(nmo) ! b-occupation numbers
c
      integer :: ilo ! lower index
      integer :: ihi ! higher index
c
      integer :: iproc ! this processor
c
      integer :: ii  ! counter
c
      double precision, parameter :: tol = 1.0d-4
c
      iproc = ga_nodeid()
c
      if (iproc.eq.0) then
        ilo = 1
        ihi = nmo
        do ii = 1, nmo
          if (occa(ii)+occb(ii).gt.2.0d0-tol) ilo = ii
          if (occa(nmo-ii+1)+occb(nmo-ii+1).lt.tol) ihi = nmo-ii+1
        enddo
c
 600    format(//,6x,"Natural orbital occupation numbers")
        write(LuOut,600)
        write(LuOut,'(/,6x,4x,"Alpha",5x,"Beta")')
        write(LuOut,'(  6x,3x,"------",3x,"------")')
        do ii = ilo, ihi
          write(LuOut,'(i6,2f9.4)')ii,occa(ii),occb(ii)
        enddo
      endif
c
      end
c
c-----------------------------------------------------------------------
