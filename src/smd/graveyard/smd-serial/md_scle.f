c
c $Id$
c

      SUBROUTINE md_scle(ntshel)

      implicit none

      integer i,ntshel

      if(ntshel.gt.0)then

       do i=1,4

        call com_motion()
        call scle_velo()
        call shll_qnch(ntshel)

       enddo

      else

       call com_motion()
       call scle_velo()

      endif

      return

      END
