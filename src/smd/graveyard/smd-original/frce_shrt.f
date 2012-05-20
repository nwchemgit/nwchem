      SUBROUTINE frce_shrt(iii,ntype,rij,rijsq,jbeg,jend,evdw)

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_vlst.inc'
      include 'cm_pote.inc'
      include 'cm_cuto.inc'

      integer iii,i,j,k,l
      integer ntype,potindex
      integer jbeg,jend,imin,imax,idif

      real*8  evdw,rij,drij,rijsq,force,etmp

      dimension rij(mxnlist,3),rijsq(mxnlist)

      evdw=0.0
      k=0

      do i=jbeg,jend

       etmp = 0.0d0
       j=list(i)
       k=k+1


       if(rijsq(k).lt.rcutsq)then

        drij=sqrt(rijsq(k))

        imin=min(atmtype(iii),atmtype(j))
        imax=max(atmtype(iii),atmtype(j))
        idif=imax-imin

        potindex=0
        do l=1,imin-1
         potindex=potindex+(ntype-l+1)
        enddo
        potindex=potindex+idif+1


        if(potkey(potindex).eq.1)then

         evdw=evdw+(potpar(potindex,1)/drij**12
     $             -potpar(potindex,2)/drij**6)
         etmp = (potpar(potindex,1)/drij**12
     $             -potpar(potindex,2)/drij**6)

         force=(12*potpar(potindex,1)/drij**12
     $          -6*potpar(potindex,2)/drij**6)/rijsq(k)


         fff(iii,1)=fff(iii,1)+convfct2*force*rij(k,1)
         fff(iii,2)=fff(iii,2)+convfct2*force*rij(k,2)
         fff(iii,3)=fff(iii,3)+convfct2*force*rij(k,3)

         fff(j,1)=fff(j,1)-convfct2*force*rij(k,1)
         fff(j,2)=fff(j,2)-convfct2*force*rij(k,2)
         fff(j,3)=fff(j,3)-convfct2*force*rij(k,3)

        elseif(potkey(potindex).eq.2)then

         evdw=evdw+(potpar(potindex,1)
     $              *exp(-drij/potpar(potindex,2))
     $              -potpar(potindex,3)/drij**6)

         force=(drij*potpar(potindex,1)
     $         *exp(-drij/potpar(potindex,2))
     $         /potpar(potindex,2)-6*potpar(potindex,3)/drij**6)
     $         /rijsq(k)

         fff(iii,1)=fff(iii,1)+convfct2*force*rij(k,1)
         fff(iii,2)=fff(iii,2)+convfct2*force*rij(k,2)
         fff(iii,3)=fff(iii,3)+convfct2*force*rij(k,3)

         fff(j,1)=fff(j,1)-convfct2*force*rij(k,1)
         fff(j,2)=fff(j,2)-convfct2*force*rij(k,2)
         fff(j,3)=fff(j,3)-convfct2*force*rij(k,3)

        endif

       endif

      enddo

      return

      END
c $Id$
