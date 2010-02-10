c
c $Id: md_scle.f,v 1.1 2008-10-01 22:32:28 marat Exp $
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
