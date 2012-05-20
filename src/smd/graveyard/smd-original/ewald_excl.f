      SUBROUTINE ewald_excl(iii,rij,rijsq,jbeg,jend,ewald3)

      implicit none

      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_ewld.inc'
      include 'cm_cuto.inc'
      include 'cm_elst.inc'

      integer iii,i,j,k
      integer jbeg,jend,iatm,jatm

      real*8 rij,drij,arij,rijsq
      real*8 ewald3,erfxc,force

      dimension rij(mxnlist,3),rijsq(mxnlist)

      ewald3=0.0
      k=0
      iatm=atmtype(iii)

      do i=jbeg,jend

       j=elist(i)
       jatm=atmtype(j)
       k=k+1

       if(rijsq(k).lt.rcutsq)then

        drij=sqrt(rijsq(k))
        arij=alpha*drij

        ewald3=ewald3-convfct1*typchge(iatm)*typchge(jatm)
     $       *(1-erfxc(arij))/drij

        force=-convfct1*typchge(iatm)*typchge(jatm)*
     $       ((1-erfxc(arij))-2*arij/sqrpi*exp(-arij*arij))
     $       /(drij*rijsq(k))

        fff(iii,1)=fff(iii,1)+convfct2*force*rij(k,1)
        fff(iii,2)=fff(iii,2)+convfct2*force*rij(k,2)
        fff(iii,3)=fff(iii,3)+convfct2*force*rij(k,3)

        fff(j,1)=fff(j,1)-convfct2*force*rij(k,1)
        fff(j,2)=fff(j,2)-convfct2*force*rij(k,2)
        fff(j,3)=fff(j,3)-convfct2*force*rij(k,3)

       endif

      enddo

      return

      END
c $Id$
