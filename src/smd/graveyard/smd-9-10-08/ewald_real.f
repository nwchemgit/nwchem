c
c $Id$
c

      SUBROUTINE ewald_real(iii,rij,rijsq,jbeg,jend,ewald4)

      implicit none

      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_ewld.inc'
      include 'cm_cuto.inc'
      include 'cm_vlst.inc'

      integer iii,i,j,k,iatm,jatm
      integer jbeg,jend

      real*8 rij,drij,arij,rijsq
      real*8 ewald4,erfxc,force

      dimension rij(mxnlist,3),rijsq(mxnlist)

      ewald4=0.0
      k=0
      iatm=atmtype(iii)

      do i=jbeg,jend

       j=list(i)
       jatm=atmtype(j)
       k=k+1

       if(rijsq(k).lt.rcutsq)then

        drij=sqrt(rijsq(k))
        arij=alpha*drij

        ewald4=ewald4+convfct1*typchge(iatm)*typchge(jatm)
     $        *erfxc(arij)/drij

        force=convfct1*typchge(iatm)*typchge(jatm)*
     $       (erfxc(arij)+2*arij/sqrpi*exp(-arij*arij))/(drij*rijsq(k))

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
