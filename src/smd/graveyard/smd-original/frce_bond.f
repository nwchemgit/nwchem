      SUBROUTINE frce_bond(ntbond,ebond)

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_bond.inc'
      include 'cm_latt.inc'

      integer i,iatm1,iatm2,ibond
      integer ntbond

      real*8 ebond,rij,dij0,rdij,dijt,force,kcst

      dimension rij(mxbond2,3)

      ebond=0.0

      do i=1,ntbond

       iatm1=bonlist(i,1)
       iatm2=bonlist(i,2)

       rij(i,1)=ccc(iatm1,1)-ccc(iatm2,1)
       rij(i,2)=ccc(iatm1,2)-ccc(iatm2,2)
       rij(i,3)=ccc(iatm1,3)-ccc(iatm2,3)

      enddo

      call tool_rebox(ntbond,mxbond2,latt,rlatt,rij)

      do i=1,ntbond

       iatm1=bonlist(i,1)
       iatm2=bonlist(i,2)
       ibond=bonlist(i,3)

       dij0=bonddist(ibond)
       kcst=bondfrce(ibond)

       dijt=sqrt(rij(i,1)**2+rij(i,2)**2+rij(i,3)**2)
       rdij=1.0/dijt

       ebond=ebond+0.5*kcst*(dijt-dij0)**2

       force=-kcst*(dijt-dij0)*rdij

       fff(iatm1,1)=fff(iatm1,1)+convfct2*force*rij(i,1)
       fff(iatm1,2)=fff(iatm1,2)+convfct2*force*rij(i,2)
       fff(iatm1,3)=fff(iatm1,3)+convfct2*force*rij(i,3)

       fff(iatm2,1)=fff(iatm2,1)-convfct2*force*rij(i,1)
       fff(iatm2,2)=fff(iatm2,2)-convfct2*force*rij(i,2)
       fff(iatm2,3)=fff(iatm2,3)-convfct2*force*rij(i,3)

      enddo

      return

      END
c $Id$
