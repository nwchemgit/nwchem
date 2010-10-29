c
c $Id$
c

      SUBROUTINE ewald_recp(ewald2)

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_ewld.inc'
      include 'cm_latt.inc'

      integer i,iatm,k,kx,ky,kz,kminx,kminy,kminz

      real*8 rksq,rx,ry,rz,rkx,rky,rkz
      real*8 kcoeff,factor,force,ewald2

      complex rhosum
      complex eikr(mxatms)
      complex eikx(1:mxatms,       0:mxkvect)
      complex eiky(1:mxatms,-mxkvect:mxkvect)
      complex eikz(1:mxatms,-mxkvect:mxkvect)

      do i=1,natms

       eikx(i,0)=(1.0,0.0)
       eiky(i,0)=(1.0,0.0)
       eikz(i,0)=(1.0,0.0)
       rx=rlatt(1,1)*ccc(i,1)+rlatt(1,2)*ccc(i,2)+rlatt(1,3)*ccc(i,3)
       ry=rlatt(2,1)*ccc(i,1)+rlatt(2,2)*ccc(i,2)+rlatt(2,3)*ccc(i,3)
       rz=rlatt(3,1)*ccc(i,1)+rlatt(3,2)*ccc(i,2)+rlatt(3,3)*ccc(i,3)
       eikx(i,1)=dcmplx(dcos(twopi*rx),dsin(twopi*rx))
       eiky(i,1)=dcmplx(dcos(twopi*ry),dsin(twopi*ry))
       eikz(i,1)=dcmplx(dcos(twopi*rz),dsin(twopi*rz))
       eiky(i,-1)=conjg(eiky(i,1))
       eikz(i,-1)=conjg(eikz(i,1))

      enddo

      do i=1,natms

       do k=2,kmaxx
        eikx(i,k)=eikx(i,k-1)*eikx(i,1)
       enddo
       do k=2,kmaxy
        eiky(i,k)=eiky(i,k-1)*eiky(i,1)
        eiky(i,-k)=conjg(eiky(i,k))
       enddo
       do k=2,kmaxz
        eikz(i,k)=eikz(i,k-1)*eikz(i,1)
        eikz(i,-k)=conjg(eikz(i,k))
       enddo

      enddo

      kminx=0
      kminy=-kmaxy
      kminz=-kmaxz

      do kx=kminx,kmaxx

       if(kx.eq.0)then
        factor=1.0
       else
        factor=2.0
       endif

       do ky=kminy,kmaxy

        do kz=kminz,kmaxz

         rkx=real(kx)*rlatt(1,1)+real(ky)*rlatt(1,2)+real(kz)*rlatt(1,3)
         rky=real(kx)*rlatt(2,1)+real(ky)*rlatt(2,2)+real(kz)*rlatt(2,3)
         rkz=real(kx)*rlatt(3,1)+real(ky)*rlatt(3,2)+real(kz)*rlatt(3,3)
         rkx=twopi*rkx
         rky=twopi*rky
         rkz=twopi*rkz
         rksq=rkx*rkx+rky*rky+rkz*rkz

          if(rksq.lt.rksqmax.and.rksq.ne.0.0)then

           rhosum=(0.0,0.0)

           do i=1,natms

            iatm=atmtype(i)
            eikr(i)=eikx(i,kx)*eiky(i,ky)*eikz(i,kz)
            rhosum=rhosum+typchge(iatm)*eikr(i)

           enddo

           kcoeff=exp(rksq*ralphsq)/rksq
           ewald2=ewald2+factor*kcoeff*conjg(rhosum)*rhosum

           do i=1,natms

            iatm=atmtype(i)
            force=-factor*2.0*twopi*convfct1/vol*kcoeff*
     $            aimag(rhosum*conjg(eikr(i)))*typchge(iatm)
            fff(i,1)=fff(i,1)+convfct2*rkx*force
            fff(i,2)=fff(i,2)+convfct2*rky*force
            fff(i,3)=fff(i,3)+convfct2*rkz*force

           enddo

          endif

         enddo

        enddo 

       enddo 

       ewald2=twopi*ewald2*convfct1/vol

       return

       END
