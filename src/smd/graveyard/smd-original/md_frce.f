      SUBROUTINE md_frce(ntype,ecoul,eshrt,ebond,ntbond,eshel,ntshel,
     $     ewald1)

      implicit none

      include 'p_array.inc'
      include 'p_input.inc'
      include 'cm_atom.inc'
      include 'cm_latt.inc'
      include 'cm_vlst.inc'
      include 'cm_elst.inc'

      integer i,j,k
      integer jbeg,jend,jnab
      integer ntype,ntbond,ntshel

      real*8 rij,rijsq
      real*8 ewald1,ewald2,ewald3,ewald4
      real*8 evdw,ecoul,eshrt,ebond,eshel
      double precision ewald_excluded

      dimension rij(mxnlist,3),rijsq(mxnlist)

      ewald_excluded = 0.0d0
      do i=1,natms
       fff(i,1)=0.0
       fff(i,2)=0.0
       fff(i,3)=0.0
      enddo

      ecoul=ewald1
      eshrt=0.0
      ebond=0.0
      eshel=0.0
      evdw=0.0
      ewald2=0.0
      ewald3=0.0
      ewald4=0.0

c Reciprocal space sum
      call ewald_recp(ewald2)
      write(*,*) "ewald reciprocal",ewald2
      ecoul=ecoul+ewald2

      do i=1,natms-1

c Excluded electrostatic interactions
       jbeg=epoint(i)
       jend=epoint(i+1)-1
       k=0

       if(jbeg.le.jend)then

        do jnab=jbeg,jend

         j=elist(jnab)
         k=k+1

         if(k.gt.mxnlist)then
          write(output,*)'k greater than mxnlist'
          stop
         endif

         rij(k,1)=ccc(i,1)-ccc(j,1)
         rij(k,2)=ccc(i,2)-ccc(j,2)
         rij(k,3)=ccc(i,3)-ccc(j,3)

        enddo

        call tool_rebox(k,mxnlist,latt,rlatt,rij)

        k=0
        do jnab=jbeg,jend

         k=k+1
         rijsq(k)=rij(k,1)*rij(k,1)+rij(k,2)*rij(k,2)+rij(k,3)*rij(k,3)

        enddo

        call ewald_excl(i,rij,rijsq,jbeg,jend,ewald3)
        ewald_excluded = ewald_excluded + ewald3
        ecoul=ecoul+ewald3

       endif

       jbeg=point(i)
       jend=point(i+1)-1
       k=0

       if(jbeg.le.jend)then

        do jnab=jbeg,jend

         j=list(jnab)
         k=k+1

         if(k.gt.mxnlist)then
          write(output,*)'k greater than mxnlist'
          stop
         endif

         rij(k,1)=ccc(i,1)-ccc(j,1)
         rij(k,2)=ccc(i,2)-ccc(j,2)
         rij(k,3)=ccc(i,3)-ccc(j,3)

        enddo

        call tool_rebox(k,mxnlist,latt,rlatt,rij)

        k=0
        do jnab=jbeg,jend

         k=k+1
         rijsq(k)=rij(k,1)*rij(k,1)+rij(k,2)*rij(k,2)+rij(k,3)*rij(k,3)

        enddo

c Real space sum
        call ewald_real(i,rij,rijsq,jbeg,jend,ewald4)
        ecoul=ecoul+ewald4

c Non-bonded interactions
        call frce_shrt(i,ntype,rij,rijsq,jbeg,jend,evdw)
        eshrt=eshrt+evdw
        write(*,*) "evdw",evdw

       endif

      enddo

      write(*,*) "ewald exluded",ewald_excluded
c Bonded interactions
c      if(ntbond.gt.0)call frce_bond(ntbond,ebond)

c Core-shell interactions
c      if(ntshel.gt.0)call frce_shel(ntshel,eshel)

      return

      END
c $Id$
