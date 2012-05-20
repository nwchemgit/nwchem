      SUBROUTINE inte_leapf(tstep,ekin)

      implicit none

      include 'p_array.inc'
      include 'cm_atom.inc'
      include 'cm_latt.inc'

      integer i,iatm

      real*8  ekin,tmpvx,tmpvy,tmpvz,tstep

      ekin=0.0

      do i=1,natms

       iatm=atmtype(i)

       tmpvx=vvv(i,1)
       tmpvy=vvv(i,2)
       tmpvz=vvv(i,3)

       vvv(i,1)=vvv(i,1)+fff(i,1)*tstep/typmass(iatm)
       vvv(i,2)=vvv(i,2)+fff(i,2)*tstep/typmass(iatm)
       vvv(i,3)=vvv(i,3)+fff(i,3)*tstep/typmass(iatm)

       tmpvx=0.5*(tmpvx+vvv(i,1))
       tmpvy=0.5*(tmpvy+vvv(i,2))
       tmpvz=0.5*(tmpvz+vvv(i,3))

       ekin=ekin+typmass(iatm)*(tmpvx**2+tmpvy**2+tmpvz**2)

       ccc(i,1)=ccc(i,1)+tstep*vvv(i,1)
       ccc(i,2)=ccc(i,2)+tstep*vvv(i,2)
       ccc(i,3)=ccc(i,3)+tstep*vvv(i,3)

      enddo 

      call tool_rebox(natms,mxatms,latt,rlatt,ccc) 

      ekin=0.5*ekin

      return

      END
c $Id$
