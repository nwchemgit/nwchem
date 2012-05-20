
      subroutine smd_leapf(natms,tstep,ekin,mass,fff,vvv,ccc)

      implicit none

      integer natms
      double precision tstep
      double precision ekin
      double precision mass(natms)
      double precision fff(natms,3)
      double precision vvv(natms,3)
      double precision ccc(natms,3)
c
      integer i,iatm

      double precision  tmpvx,tmpvy,tmpvz

      ekin=0.0

      do i=1,natms
       tmpvx=vvv(i,1)
       tmpvy=vvv(i,2)
       tmpvz=vvv(i,3)

       vvv(i,1)=vvv(i,1)+fff(i,1)*tstep/mass(i)
       vvv(i,2)=vvv(i,2)+fff(i,2)*tstep/mass(i)
       vvv(i,3)=vvv(i,3)+fff(i,3)*tstep/mass(i)

       tmpvx=0.5*(tmpvx+vvv(i,1))
       tmpvy=0.5*(tmpvy+vvv(i,2))
       tmpvz=0.5*(tmpvz+vvv(i,3))

       ekin=ekin+mass(i)*(tmpvx**2+tmpvy**2+tmpvz**2)

       ccc(i,1)=ccc(i,1)+tstep*vvv(i,1)
       ccc(i,2)=ccc(i,2)+tstep*vvv(i,2)
       ccc(i,3)=ccc(i,3)+tstep*vvv(i,3)
       
      enddo 

      ekin=0.5*ekin

      return

      end

c $Id$
