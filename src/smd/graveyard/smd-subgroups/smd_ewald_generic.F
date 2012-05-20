
      subroutine smd_ewald_recip_generic(
     >                          na,
     >                          nk,
     >                          eikr,
     >                          eikx,
     >                          eiky,
     >                          eikz,
     >                          ralphsq,
     >                          rksqmax,
     >                          vol,
     >                          rlatt,
     >                          kmax,
     >                          ccc,
     >                          q,
     >                          fff,
     >                          ewald2)

      implicit none
#include "smd_const_data.fh"

      integer na
      integer nk
      double precision ralphsq
      double precision rksqmax
      double precision vol
      double complex eikr(na)
      double complex eikx(1:na,0:nk)
      double complex eiky(1:na,-nk:nk)
      double complex eikz(1:na,-nk:nk)
      double precision rlatt(3,3)
      double precision ccc(na,3)
      double precision fff(na,3)
      double precision q(na)
      integer kmax(3)



      integer i,k,kx,ky,kz,kminx,kminy,kminz

      real*8 rksq,rx,ry,rz,rkx,rky,rkz
      real*8 kcoeff,factor,force,ewald2

      double complex rhosum

      ewald2 = 0.0d0
      do i=1,na

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

      do i=1,na

       do k=2,kmax(1)
        eikx(i,k)=eikx(i,k-1)*eikx(i,1)
       enddo
       do k=2,kmax(2)
        eiky(i,k)=eiky(i,k-1)*eiky(i,1)
        eiky(i,-k)=conjg(eiky(i,k))
       enddo
       do k=2,kmax(3)
        eikz(i,k)=eikz(i,k-1)*eikz(i,1)
        eikz(i,-k)=conjg(eikz(i,k))
       enddo

      enddo

      kminx=0
      kminy=-kmax(2)
      kminz=-kmax(3)

      do kx=kminx,kmax(1)

       if(kx.eq.0)then
        factor=1.0
       else
        factor=2.0
       endif

       do ky=kminy,kmax(2)

        do kz=kminz,kmax(3)

         rkx=real(kx)*rlatt(1,1)+real(ky)*rlatt(1,2)+real(kz)*rlatt(1,3)
         rky=real(kx)*rlatt(2,1)+real(ky)*rlatt(2,2)+real(kz)*rlatt(2,3)
         rkz=real(kx)*rlatt(3,1)+real(ky)*rlatt(3,2)+real(kz)*rlatt(3,3)
         rkx=twopi*rkx
         rky=twopi*rky
         rkz=twopi*rkz
         rksq=rkx*rkx+rky*rky+rkz*rkz

          if(rksq.lt.rksqmax.and.rksq.ne.0.0)then

           rhosum=(0.0,0.0)

           do i=1,na

            eikr(i)=eikx(i,kx)*eiky(i,ky)*eikz(i,kz)
            rhosum=rhosum+q(i)*eikr(i)

           enddo

           kcoeff=exp(rksq*ralphsq)/rksq
           ewald2=ewald2+factor*kcoeff*conjg(rhosum)*rhosum
           do i=1,na

            force=-factor*2.0*twopi*convfct1/vol*kcoeff*
     $            dimag(rhosum*dconjg(eikr(i)))*q(i)
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


      subroutine smd_ewald_excl_generic(na,
     >                                  nl,
     >                                  alpha,
     >                                  rcutsq,
     >                                  latt,
     >                                  rlatt,
     >                                  q,
     >                                  ccc,
     >                                  fff,
     >                                  epoint,
     >                                  elist,
     >                                  e)

      implicit none

#include "smd_const_data.fh"


      integer na
      integer nl
      double precision alpha
      double precision rcutsq
      double precision rlatt(3,3),latt(3,3)
      integer epoint(na)
      double precision q(na)
      integer elist(nl)
      double precision ccc(na,3),fff(na,3)
      double precision e
c
      integer i,j,k,jnab
      integer jbeg,jend

      double precision  dr,ar,rsq
      double precision  erfxc,force

      double precision x,y,z

      double precision  ssx,ssy,ssz,xss,yss,zss

      e=0

      do i=1,na-1

       jbeg=epoint(i)
       jend=epoint(i+1)-1

       write(*,*) "i,jbeg,jend",i,jbeg,jend

      do jnab=jbeg,jend

       j=elist(jnab)

       x=ccc(i,1)-ccc(j,1)
       y=ccc(i,2)-ccc(j,2)
       z=ccc(i,3)-ccc(j,3)
c
c      reboxing here
c      ------------
       ssx=(rlatt(1,1)*x+rlatt(1,2)*y+rlatt(1,3)*z)
       ssy=(rlatt(2,1)*x+rlatt(2,2)*y+rlatt(2,3)*z)
       ssz=(rlatt(3,1)*x+rlatt(3,2)*y+rlatt(3,3)*z)

       xss=ssx-nint(ssx)
       yss=ssy-nint(ssy)
       zss=ssz-nint(ssz)

       x=(latt(1,1)*xss+latt(1,2)*yss+latt(1,3)*zss)
       y=(latt(2,1)*xss+latt(2,2)*yss+latt(2,3)*zss)
       z=(latt(3,1)*xss+latt(3,2)*yss+latt(3,3)*zss)
c      done reboxing

       rsq=x*x+y*y+z*z
       if(rsq.lt.rcutsq)then

        dr=sqrt(rsq)
        ar=alpha*dr

        e=e-convfct1*q(i)*q(j)
     $       *(1-erfxc(ar))/dr

        force=-convfct1*q(i)*q(j)*
     $       ((1-erfxc(ar))-2*ar/sqrpi*exp(-ar*ar))
     $       /(dr*rsq)

        fff(i,1)=fff(i,1)+convfct2*force*x
        fff(i,2)=fff(i,2)+convfct2*force*y
        fff(i,3)=fff(i,3)+convfct2*force*z

        fff(j,1)=fff(j,1)-convfct2*force*x
        fff(j,2)=fff(j,2)-convfct2*force*y
        fff(j,3)=fff(j,3)-convfct2*force*z

       endif

      end do
      end do

      return

      END

      subroutine smd_ewald_real_generic(na,
     >                                  nl,
     >                                  alpha,
     >                                  rcutsq,
     >                                  q,
     >                                  ccc,
     >                                  fff,
     >                                  point,
     >                                  list,
     >                                  e)

      implicit none

#include "smd_const_data.fh"


      integer na
      integer nl
      double precision alpha
      double precision rcutsq
      integer point(na)
      double precision q(na)
      integer list(nl)
      double precision ccc(nl,3)
      double precision fff(na,3)
      double precision e
c
      integer i,j,k,jnab
      integer jbeg,jend
      integer nlist

      double precision  dr,ar,rsq
      double precision  erfxc,force

      double precision x,y,z

      e=0

      nlist = 0
      do i=1,na-1

       jbeg=point(i)
       jend=point(i+1)-1


      do jnab=jbeg,jend

       j=list(jnab)

       nlist = nlist + 1
       x=ccc(nlist,1)
       y=ccc(nlist,2)
       z=ccc(nlist,3)

       rsq=x*x+y*y+z*z

       if(rsq.lt.rcutsq)then

        dr=sqrt(rsq)
        ar=alpha*dr

        e=e+convfct1*q(i)*q(j)
     $       *erfxc(ar)/dr
        write(*,*) q(i),q(j),ar,dr

        force=convfct1*q(i)*q(j)
     $       *(erfxc(ar)+2*ar/sqrpi*exp(-ar*ar))/(dr*rsq)

        fff(i,1)=fff(i,1)+convfct2*force*x
        fff(i,2)=fff(i,2)+convfct2*force*y
        fff(i,3)=fff(i,3)+convfct2*force*z

        fff(j,1)=fff(j,1)-convfct2*force*x
        fff(j,2)=fff(j,2)-convfct2*force*y
        fff(j,3)=fff(j,3)-convfct2*force*z

       endif

      end do
      end do

      return

      END

c $Id$
