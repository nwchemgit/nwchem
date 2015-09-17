c
c     Print calculation information
c     HVD, NG: 2012
c
      subroutine tddft_grad_print_info(luout,nroots,iroots,
     +           nsroots,isroots,ntroots,itroots,title)
c
      implicit none
c
#include "global.fh"
c
c     Inputs:
      character*(*) title      ! the job title
c
      integer luout            ! the output unit
      integer nroots           ! the number of unrestricted roots
      integer iroots(nroots)   ! the unrestricted roots
      integer nsroots          ! the number of singlet roots
      integer isroots(nsroots) ! the singlet roots
      integer ntroots          ! the number of triplet roots
      integer itroots(ntroots) ! the triplet roots
c
c     Local:
c
      integer i
      character*32 pname
c
      pname="tddft_grad_print_info: "
c
      if (ga_nodeid().ne.0) return
      if (title .ne. ' ') then
        call util_print_centered(luout, title, 40, .false.)
        write(luout,*)
        write(luout,*)
      endif
      write(luout,*)'Calculated gradients of:'
      if (nroots.gt.0) then
        write(luout,"(' Number of Unrestricted roots ',i3)")nroots
        write(luout,"(' Unrestricted roots ',$)")
        do i = 1, nroots
          write(luout,"(i4,$)")iroots(i)
        enddo
        write(luout,*)
      else
        if (nsroots+ntroots.gt.0) then
          if (nsroots.gt.0) then
            write(luout,"(' Number of Singlet roots ',i3)")nsroots
            write(luout,"(' Singlet roots ',$)")
            do i = 1, nsroots
              write(luout,"(i4,$)")isroots(i)
            enddo
            write(luout,*)
          endif
          if (ntroots.gt.0) then
            write(luout,"(' Number of Triplet roots ',i3)")ntroots
            write(luout,"(' Triplet roots ',$)")
            do i = 1, ntroots
              write(luout,"(i4,$)")itroots(i)
            enddo
            write(luout,*)
          endif
        else
          write(luout,*)'No roots whatsoever ...'
        endif
c
      endif  ! ga_nodeid()
c
      end
c $Id$
