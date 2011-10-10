
      subroutine smd_leapf(natms,tstep,ekin,mass,fff,vvv,ccc)

      implicit none

      integer natms
      double precision tstep
      double precision ekin
      double precision mass(natms)
      double precision fff(3,natms)
      double precision vvv(3,natms)
      double precision ccc(3,natms)
c
      integer i,iatm

      double precision  tmpvx,tmpvy,tmpvz

      ekin=0.0

      do i=1,natms
       tmpvx=vvv(1,i)
       tmpvy=vvv(2,i)
       tmpvz=vvv(3,i)

       vvv(1,i)=vvv(1,i)+fff(1,i)*tstep/mass(i)
       vvv(2,i)=vvv(2,i)+fff(2,i)*tstep/mass(i)
       vvv(3,i)=vvv(3,i)+fff(3,i)*tstep/mass(i)

       tmpvx=0.5*(tmpvx+vvv(1,i))
       tmpvy=0.5*(tmpvy+vvv(2,i))
       tmpvz=0.5*(tmpvz+vvv(3,i))

       ekin=ekin+mass(i)*(tmpvx**2+tmpvy**2+tmpvz**2)

       ccc(1,i)=ccc(1,i)+tstep*vvv(1,i)
       ccc(2,i)=ccc(2,i)+tstep*vvv(2,i)
       ccc(3,i)=ccc(3,i)+tstep*vvv(3,i)
       
      enddo 

      ekin=0.5*ekin

      return

      end

c $Id$
