      subroutine obftrans(ovbf,ovbftrn,tr,nsy,nchan,nscat,nsch,nlm,
     1 ngauss,ngch,nbfmax,lmtop,nchnl,iprint,ovbftrnp,ovbfp,istat,ioft)
       implicit real*8 (a-h,o-z)
       character*8 istat,ioft
c
c  transform the bound index of the bound-free overlaps
c
      real*8 tr(nbfmax,nbfmax)
      complex*16 ovbf(lmtop,nbfmax,nchnl),ovbfp(lmtop,nbfmax,nchnl)
      complex*16 ovbftrn(lmtop,nbfmax,nchnl)
      complex*16 ovbftrnp(lmtop,nbfmax,nchnl)
      integer nscat(nchnl),ngauss(nchnl),nlm(nchnl),nsch(nbfmax,nchnl)
      integer ngch(nbfmax,nchnl)
      do 100 ic=1,nchan
      nlmic=nlm(ic)
      nsic=nscat(ic)
      do 110 i=1,nlmic
      do 110 j=1,nsic
      ovbftrnp(i,j,ic)=0.0
110   ovbftrn(i,j,ic)=0.0
      ngic=ngauss(ic)
      if(istat.eq.ioft)then
      do 102 kkc=1,ngic
      do 102 ilm=1,nlmic
      do 102 isc=1,nsic
 102     ovbftrnp(ilm,isc,ic)=ovbftrnp(ilm,isc,ic) +
     1   ovbfp(ilm,kkc,ic)*tr(ngch(kkc,ic),nsch(isc,ic))
      endif
      do 100 kkc=1,ngic
      do 100 ilm=1,nlmic
      do 100 isc=1,nsic
100   ovbftrn(ilm,isc,ic)=ovbftrn(ilm,isc,ic) +
     1   ovbf(ilm,kkc,ic)*tr(ngch(kkc,ic),nsch(isc,ic))
      if(iprint.ne.0) then
      do 200 ic=1,nchan
      write(6,107) ic
107   format(//' transformed bound free overlaps for channel:',i4)
      nlmic=nlm(ic)
      nsic=nscat(ic)
      do 200 ilm=1,nlmic
200   write(6,101) ilm,(ovbftrn(ilm,j,ic),j=1,nsic)
101   format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &     (4x,3("(",f8.5,3x,f8.5,")",3x)))
      if(istat.eq.ioft)then
      do 201 ic=1,nchan
      write(6,108) ic
 108  format(//' transformed bound free (p) overlaps for channel:',i4)
      nlmic=nlm(ic)
      nsic=nscat(ic)
      do 201 ilm=1,nlmic
 201     write(6,101) ilm,(ovbftrnp(ilm,j,ic),j=1,nsic)     
      endif
      endif
      return
      end
c $Id$
