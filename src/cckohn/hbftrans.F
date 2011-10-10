      subroutine hbftrans(hbf,hbftrn,tr,nsy,nchan,nscat,nsch,nlm,
     1 ngauss,ngch,nbfmax,lmtop,nchnl,iprint,hbfp,hbftrnp,istat,ioft)
       implicit real*8 (a-h,o-z)
c
c  transform the bound index of the bound-free hamiltonian
c
      real*8 tr(nbfmax,nbfmax)
      complex*16 hbf(lmtop,nbfmax,nchnl**2),hbfp(lmtop,nbfmax,nchnl**2)
      complex*16 hbftrn(lmtop,nbfmax,nchnl**2),hbftrnp(lmtop,nbfmax
     1 ,nchnl**2)
      integer nscat(nchnl),ngauss(nchnl),nlm(nchnl),nsch(nbfmax,nchnl)
      integer ngch(nbfmax,nchnl)
      character*8 istat,ioft
      do 102 ic=1,nchan
      nlmic=nlm(ic)
      do 102 jc=1,nchan
      nsjc=nscat(jc)
      ngjc=ngauss(jc)
      icc=nchan*(ic-1) + jc
      do 80 i=1,nlmic
      do 80 j=1,nsjc
      hbftrnp(i,j,icc) = 0.
80    hbftrn(i,j,icc) = 0.
      do 100 kkc=1,ngjc
      do 100 ilm=1,nlmic
      do 100 jsc=1,nsjc
100    hbftrn(ilm,jsc,icc)=hbftrn(ilm,jsc,icc) +
     1   hbf(ilm,kkc,icc)*tr(ngch(kkc,jc),nsch(jsc,jc))
      if(istat.eq.ioft)then
      do 103 kkc=1,ngjc
      do 103 ilm=1,nlmic
      do 103 jsc=1,nsjc
103    hbftrnp(ilm,jsc,icc)=hbftrnp(ilm,jsc,icc) +
     1   hbfp(ilm,kkc,icc)*tr(ngch(kkc,jc),nsch(jsc,jc))
      endif
 102  continue
      if(iprint.ne.0) then
      do 200 ic=1,nchan
      nlmic=nlm(ic)
      do 200 jc=1,nchan
      write(6,107)ic,jc
107   format(//,' transformed bound free ham. matrix for channels:',2i4)
      nsjc=nscat(jc)
      icc=nchan*(ic-1)+jc
      do 200 ilm=1,nlmic
200   write(6,101) ilm,(hbftrn(ilm,j,icc),j=1,nsjc)
101   format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &     (4x,3("(",f8.5,3x,f8.5,")",3x)))
      if(istat.eq.ioft)then
      do 201 ic=1,nchan
      nlmic=nlm(ic)
      do 201 jc=1,nchan
      write(6,108)ic,jc
 108  format(//,' transformed bound free (p) ham. matrix for channels:'
     1 ,2i4)
      nsjc=nscat(jc)
      icc=nchan*(ic-1)+jc
      do 201 ilm=1,nlmic
 201     write(6,101) ilm,(hbftrn(ilm,j,icc),j=1,nsjc)
      endif
      endif
      return
      end
c $Id$
