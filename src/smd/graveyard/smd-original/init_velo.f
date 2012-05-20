      SUBROUTINE init_velo(iseed)

      implicit none

      include 'p_array.inc'
      include 'cm_atom.inc'
      include 'cm_temp.inc'

      integer i,iatm,iseed

      real*8 x
      real*8 commass,comxvv,comyvv,comzvv
      real*8 instanke,xscale

      iseed=620419483

      comxvv=0.d0
      comyvv=0.d0
      comzvv=0.d0
      commass = 0.0d0

      do i=1,natms

       iatm=atmtype(i)
       call tool_randm(iseed,x)
       vvv(i,1)=(x-0.5)/sqrt(typmass(iatm))
       call tool_randm(iseed,x)
       vvv(i,2)=(x-0.5)/sqrt(typmass(iatm))
       call tool_randm(iseed,x)
       vvv(i,3)=(x-0.5)/sqrt(typmass(iatm))

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

      instanke=0.d0

      do i=1,natms

       iatm=atmtype(i)
       instanke=instanke
     $         +typmass(iatm)*((vvv(i,1)**2+vvv(i,2)**2+vvv(i,3)**2))

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

      SUBROUTINE write_velo(un)

      implicit none

      integer un
      include 'p_array.inc'
      include 'cm_atom.inc'
      include 'cm_temp.inc'

      integer i

      do i=1,natms

       write(un,*) vvv(i,1),vvv(i,2),vvv(i,3)

      enddo
      return

      END
c $Id$
