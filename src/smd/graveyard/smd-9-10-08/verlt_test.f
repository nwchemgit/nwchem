c
c $Id$
c

      SUBROUTINE verlt_test(tstep,ivv,lupdate)

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'cm_atom.inc'
      include 'cm_cuto.inc'

      integer i,exceed

      real*8 tstep,tstepsq,ivv
      real*8 dispmax,dispxsq,dispysq,dispzsq,disprsq 

      logical lnew,lupdate

      dimension ivv(mxatms,3) 

      data lnew/.true./

      save lnew

      tstepsq=tstep**2

      if(lnew)then

       lupdate=.true.
       lnew=.false.

       do i=1,natms

        ivv(i,1)=0.0
        ivv(i,2)=0.0
        ivv(i,3)=0.0

       enddo

      else

       lupdate=.false.

       dispmax=((vcut-rcut)/2.0)**2

       do i=1,natms

        ivv(i,1)=ivv(i,1)+vvv(i,1)
        ivv(i,2)=ivv(i,2)+vvv(i,2)
        ivv(i,3)=ivv(i,3)+vvv(i,3)

       enddo

       exceed=0

       do i=1,natms

        dispxsq=ivv(i,1)**2
        dispysq=ivv(i,2)**2
        dispzsq=ivv(i,3)**2
        disprsq=tstepsq*(dispxsq+dispysq+dispzsq)
        if(disprsq.gt.dispmax)exceed=exceed+1
        if(exceed.ge.2)lupdate=.true.

       enddo

       if(lupdate)then

        do i=1,natms

         ivv(i,1)=0
         ivv(i,2)=0
         ivv(i,3)=0

        enddo

       endif

      endif

      return

      END
