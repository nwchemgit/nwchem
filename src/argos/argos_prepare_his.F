      logical function argos_prepare_his(num,mat,
     + csa,isat,isgm,imol,ifra,xs,msa,nsa)
c
c $Id$
c
      implicit none
c
#include "util.fh"
#include "mafdecls.fh"
c
      character*2 argos_prepare_atnam
      real*8 argos_prepare_atsiz
      integer argos_prepare_merge
      external argos_prepare_atnam,argos_prepare_atsiz,
     + argos_prepare_merge
c
      integer mat,msa,nsa
      integer num(mat),isat(msa),isgm(msa),imol(msa),ifra(msa)
      character*16 csa(msa)
      real*8 xs(3,msa)
c
      real*8 d,dmin
      integer i,j,jmin,k,l
c
      do 1 i=1,nsa
      if((csa(i)(1:3).eq.'HID'.or.csa(i)(1:3).eq.'HIE'.or.
     + csa(i)(1:3).eq.'HIP'.or.csa(i)(1:3).eq.'HIN').and.
     + (csa(i)(11:14).eq.' ND1'.or.csa(i)(11:14).eq.' HD1'.or.
     + csa(i)(11:14).eq.' NE2'.or.csa(i)(11:14).eq.' HE2')) then
      jmin=0
      dmin=9.9d9
      do 2 j=1,nsa
      if(isgm(i).ne.isgm(j).and.(num(isat(j)).eq.7.or.
     + num(isat(j)).eq.7.or.num(isat(j)).eq.8.or.
     + num(isat(j)).eq.16)) then
      d=sqrt((xs(1,i)-xs(1,j))**2+(xs(2,i)-xs(2,j))**2+
     + (xs(3,i)-xs(3,j))**2)
      if(d.lt.dmin) then
      dmin=d
      jmin=j
      endif
      endif
    2 continue
      if(jmin.gt.0) then
c      if(csa(i)(1:3).eq.'HID'
      write(*,10) isgm(i),csa(i),isgm(jmin),csa(jmin),d
   10 format(i5,a,i5,a,f12.6)
      endif
      endif
    1 continue
c
      argos_prepare_his=.true.
c
      return
      end
