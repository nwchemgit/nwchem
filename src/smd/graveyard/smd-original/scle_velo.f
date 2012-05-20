      SUBROUTINE scle_velo()

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_temp.inc'

      integer i,iatm

      real*8 instanke,xscale

      instanke=0.d0

      do i=1,natms

       iatm=atmtype(i)
       instanke=instanke
     $          +typmass(iatm)*((vvv(i,1)**2+vvv(i,2)**2+vvv(i,3)**2))

      enddo

      instanke=0.5*instanke
      xscale=sqrt(targetke/instanke)

      do i=1,natms

       vvv(i,1)=xscale*vvv(i,1)
       vvv(i,2)=xscale*vvv(i,2)
       vvv(i,3)=xscale*vvv(i,3)

      enddo

      return

      END
c $Id$
