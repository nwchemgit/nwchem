c
c $Id$
c

      SUBROUTINE frce_shel(ntshel,eshel)

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_latt.inc'
      include 'cm_shel.inc'

      integer i,iatm1,iatm2,ishel,ntshel

      real*8  eshel,rij,dijt,force,kcst

      dimension rij(mxshel2,3)

      eshel=0.0

      do i=1,ntshel

       iatm1=shllist(i,1)
       iatm2=shllist(i,2)

       rij(i,1)=ccc(iatm1,1)-ccc(iatm2,1)
       rij(i,2)=ccc(iatm1,2)-ccc(iatm2,2)
       rij(i,3)=ccc(iatm1,3)-ccc(iatm2,3)

      enddo

      call tool_rebox(ntshel,mxshel2,latt,rlatt,rij)

      do i=1,ntshel

       iatm1=shllist(i,1)
       iatm2=shllist(i,2)
       ishel=shllist(i,3)

       kcst=shelfrce(ishel)

       dijt=sqrt(rij(i,1)**2+rij(i,2)**2+rij(i,3)**2)

       eshel=eshel+0.5*kcst*dijt**2

       force=-kcst

       fff(iatm1,1)=fff(iatm1,1)+convfct2*force*rij(i,1)
       fff(iatm1,2)=fff(iatm1,2)+convfct2*force*rij(i,2)
       fff(iatm1,3)=fff(iatm1,3)+convfct2*force*rij(i,3)

       fff(iatm2,1)=fff(iatm2,1)-convfct2*force*rij(i,1)
       fff(iatm2,2)=fff(iatm2,2)-convfct2*force*rij(i,2)
       fff(iatm2,3)=fff(iatm2,3)-convfct2*force*rij(i,3)

      enddo

      return

      END
