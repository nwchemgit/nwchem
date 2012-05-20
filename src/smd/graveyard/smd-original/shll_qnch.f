      SUBROUTINE shll_qnch(ntshel)

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'p_const.inc'
      include 'cm_atom.inc'
      include 'cm_temp.inc'
      include 'cm_shel.inc'

      integer i,iatm1,iatm2,iatma,iatmb
      integer ntshel

      real*8 cske,rmu,vij,xscale
      real*8 tmx,tmy,tmz

      dimension vij(mxshel2,3)

      cske=boltzmann*temp*1.d-4

      do i=1,ntshel

       iatm1=shllist(i,1)
       iatm2=shllist(i,2)
       iatma=atmtype(iatm1)
       iatmb=atmtype(iatm2)
       rmu=(typmass(iatma)*typmass(iatmb))/
     $     (typmass(iatma)+typmass(iatmb))

       vij(i,1)=vvv(iatm2,1)-vvv(iatm1,1)
       vij(i,2)=vvv(iatm2,2)-vvv(iatm1,2)
       vij(i,3)=vvv(iatm2,3)-vvv(iatm1,3)

       xscale=sqrt(cske/(rmu*(vij(i,1)**2+vij(i,2)**2+vij(i,3)**2)))

       tmx=typmass(iatma)*vvv(iatm1,1)+typmass(iatmb)*vvv(iatm2,1)
       tmy=typmass(iatma)*vvv(iatm1,2)+typmass(iatmb)*vvv(iatm2,2)
       tmz=typmass(iatma)*vvv(iatm1,3)+typmass(iatmb)*vvv(iatm2,3)

       vvv(iatm1,1)=tmx/(typmass(iatma)+typmass(iatmb))
     $             -xscale*rmu*vij(i,1)/typmass(iatma)
       vvv(iatm2,1)=tmx/(typmass(iatma)+typmass(iatmb))
     $             +xscale*rmu*vij(i,1)/typmass(iatmb)
       vvv(iatm1,2)=tmy/(typmass(iatma)+typmass(iatmb))
     $             -xscale*rmu*vij(i,2)/typmass(iatma)
       vvv(iatm2,2)=tmy/(typmass(iatma)+typmass(iatmb))
     $             +xscale*rmu*vij(i,2)/typmass(iatmb)
       vvv(iatm1,3)=tmz/(typmass(iatma)+typmass(iatmb))
     $             -xscale*rmu*vij(i,3)/typmass(iatma)
       vvv(iatm2,3)=tmz/(typmass(iatma)+typmass(iatmb))
     $             +xscale*rmu*vij(i,3)/typmass(iatmb)

      enddo
 
      return

      END
c $Id$
