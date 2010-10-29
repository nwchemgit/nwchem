c
c $Id$
c

      SUBROUTINE com_motion()

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'cm_atom.inc'

      integer i,iatm

      real*8 commass,comxvv,comyvv,comzvv

      comxvv=0.d0
      comyvv=0.d0
      comzvv=0.d0
      commass=0.d0

      do i=1,natms

       iatm=atmtype(i)
       comxvv=comxvv+vvv(i,1)*typmass(iatm)
       comyvv=comyvv+vvv(i,2)*typmass(iatm)
       comzvv=comzvv+vvv(i,3)*typmass(iatm)
       commass=commass+typmass(iatm)

      enddo

      comxvv=comxvv/commass
      comyvv=comyvv/commass
      comzvv=comzvv/commass

      do i=1,natms

       vvv(i,1)=vvv(i,1)-comxvv
       vvv(i,2)=vvv(i,2)-comyvv
       vvv(i,3)=vvv(i,3)-comzvv

      enddo

      return

      END
