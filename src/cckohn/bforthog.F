c
      subroutine bforthog(ovbftrn,lmtop,nbfmax,nchnl,nstate,hpvbtrn,
     1 hbbtrn,nchan,nlm,nscat,nbas,bforth,iprint,bforthp,
     2 obbtrn,kchan,eground,ovbftrnp,hbbtrne,hpvbtrnp,istat,ioft,
     $  nsmall)
c
       implicit real*8 (a-h,o-z)
c
c construct free-bound matrix of (h-e) for free functions
c orthogonalized to bound basis
c
      complex*16 bforth(lmtop,nbfmax,nchnl**2)
      complex*16 bforthp(lmtop,nbfmax,nchnl**2)
      complex*16 hpvbtrn(lmtop,nbfmax,nchnl**2)
      complex*16 hpvbtrnp(lmtop,nbfmax,nchnl**2)
      complex*16 ovbftrn(lmtop,nbfmax,nchnl)
      complex*16 ovbftrnp(lmtop,nbfmax,nchnl)
      real*8 hbbtrn(nbfmax,nbfmax,nstate)
      real*8 hbbtrne(nbfmax,nbfmax,nstate)
      real*8 obbtrn(nbfmax,nbfmax,nstate), kchan(nchnl)
      integer nlm(nchnl),nscat(nchnl),nbas(nchnl),nsmall(nchnl)
      character*8 istat,ioft
c
c zero the bound-free output matrix
c
      do 100 ic=1,nchan
      nlmic = nlm(ic)
      do 100 jc=1,nchan
      nsjc = nscat(jc)
      icc = nchan*(ic-1) + jc
      do 100 i=1,nlmic
      do 100 j=1,nsjc
100   bforth(i,j,icc) = hpvbtrn(i,j,icc)
      if(istat.eq.ioft)then
      do 102 ic=1,nchan
      nlmic = nlm(ic)
      do 102 jc=1,nchan
      nsjc = nscat(jc)
      icc = nchan*(ic-1) + jc
      do 102 i=1,nlmic
      do 102 j=1,nsjc
102   bforthp(i,j,icc) = hpvbtrnp(i,j,icc)
      endif
c
c*******note********
c this routine assumes target energies are not added to hbb
c**************
c     e=eground + kchan(1)**2/2
c
c do the orthogonalization
c
      do 400 ic=1,nchan
      nlmic=nlm(ic)
      nbkc=nbas(ic)
      do 400 jc=1,nchan
      nsjc=nscat(jc)
      icc=nchan*(ic-1) + jc
      ii=max0(ic,jc)
      jj=ic+jc-ii
c*************
      e=kchan(ii)**2/2.
      ist=ii*(ii-1)/2+jj
      if(ic.gt.jc)then
      do 402 kbc=1,nbkc
      do 402 ilm=1,nlmic
      do 402 jsc=1,nsjc
402   bforth(ilm,jsc,icc) = bforth(ilm,jsc,icc)
     1 -ovbftrn(ilm,kbc,ic)*(hbbtrn(kbc,jsc,ist)-e*obbtrn(kbc,jsc,ist))
      else
      do 401 kbc=1,nbkc
      do 401 ilm=1,nlmic
      do 401 jsc=1,nsjc
401   bforth(ilm,jsc,icc) = bforth(ilm,jsc,icc)
     1 -ovbftrn(ilm,kbc,ic)*(hbbtrn(jsc,kbc,ist)-e*obbtrn(jsc,kbc,ist))
      endif
400   continue
      if(istat.eq.ioft)then
      do 403 ic=1,nchan
      nlmic=nlm(ic)
c      nbkc=nbas(ic)
      nbkc=nsmall(ic)
      do 403 jc=1,nchan
      nsjc=nscat(jc)
      icc=nchan*(ic-1) + jc
      ii=max0(ic,jc)
      jj=ic+jc-ii
c*************
      ist=ii*(ii-1)/2+jj
      if(ic.gt.jc)then
      do 404 kbc=1,nbkc
      do 404 ilm=1,nlmic
      do 404 jsc=1,nsjc
 404     bforthp(ilm,jsc,icc) = bforthp(ilm,jsc,icc)
     1 +ovbftrnp(ilm,kbc,ic)*hbbtrne(kbc,jsc,ist)
      nbasic=nbas(ic)
      if(nbkc.lt.nbasic)then
      istart=nbkc+1
      do 414 kbc=istart,nbasic
      do 414 ilm=1,nlmic
      do 414 jsc=1,nsjc
 414     bforthp(ilm,jsc,icc) = bforthp(ilm,jsc,icc)
     1 -ovbftrnp(ilm,kbc,ic)*(hbbtrn(kbc,jsc,ist)-e*obbtrn(kbc,jsc,ist))
      endif
      else
      do 405 kbc=1,nbkc
      do 405 ilm=1,nlmic
      do 405 jsc=1,nsjc
 405     bforthp(ilm,jsc,icc) = bforthp(ilm,jsc,icc)
     1 +ovbftrnp(ilm,kbc,ic)*hbbtrne(jsc,kbc,ist)
      nbasic=nbas(ic)
      if(nbkc.lt.nbasic)then
      istart=nbkc+1
      do 415 kbc=istart,nbasic
      do 415 ilm=1,nlmic
      do 415 jsc=1,nsjc
 415     bforthp(ilm,jsc,icc) = bforthp(ilm,jsc,icc)
     1 -ovbftrnp(ilm,kbc,ic)*(hbbtrn(jsc,kbc,ist)-e*obbtrn(jsc,kbc,ist))
      endif
      endif
 403  continue
      endif
c
      if(iprint.ne.0) then
      do 200 ic=1,nchan
      nlmic=nlm(ic)
      do 200 jc=1,nchan
      write(6,107)ic,jc
107   format(//,' bound-(orth)free ham. matrix for channels:',2i4)
      nsjc=nscat(jc)
      icc=nchan*(ic-1)+jc
      do 200 ilm=1,nlmic
200   write(6,101) ilm,(bforth(ilm,j,icc),j=1,nsjc)
c101   format(1x,i3,6e12.5,/,(4x,6e12.5))
101   format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &     (4x,3("(",f8.5,3x,f8.5,")",3x)))
      if(istat.eq.ioft)then
      do 201 ic=1,nchan
      nlmic=nlm(ic)
      do 201 jc=1,nchan
      write(6,108)ic,jc
108   format(//,' bound-freep ham. matrix for channels:',2i4)
      nsjc=nscat(jc)
      icc=nchan*(ic-1)+jc
      do 201 ilm=1,nlmic
201   write(6,101) ilm,(bforthp(ilm,j,icc),j=1,nsjc)
      endif
      endif
      return
      end
c $Id$
