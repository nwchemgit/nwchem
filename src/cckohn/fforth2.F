      subroutine fforth2(ovbftrn,lmtop,nbfmax,nchnl,nstate,
     1 hpvbtrn,bforth,nchan,nlm,nbas,fpfmorth,iprint,
     2 ibcondx,hpvbtrnp,fpfmporth,istat,ioft,nscat,ovbftrnp)
       implicit real*8 (a-h,o-z)
c
c construct free-free matrix element:
c h+(orthogonalized) (h-e) h-(orthogonalized)
c
c  ibcondx = 0 gives S-matrix boundary conditions
c  ibcondx = 1 gives T-matrix boundary conditions
c
c  Note: the code FFBF puts out <h+|V|j> in the T-matrix case
c  but it is always called hpvhm -- whether S- or T-matrix boundary
c  conditions are in effect
c
      complex*16 fpfmorth(lmtop,lmtop,nchnl**2),fpfmporth(lmtop,lmtop
     1 ,nchnl**2)
      complex*16 ovbftrn(lmtop,nbfmax,nchnl)
      complex*16 ovbftrnp(lmtop,nbfmax,nchnl)
      complex*16 hpvbtrn(lmtop,nbfmax,nchnl**2),hpvbtrnp(lmtop,nbfmax
     1 ,nchnl**2)
      complex*16 bforth(lmtop,nbfmax,nchnl**2)
      integer nbas(nchnl),nlm(nchnl),nscat(nchnl)
c
c orthogonalize wrt fcns in channel ic
c
      do 200 ic=1,nchan
      nlmic=nlm(ic)
      nbic=nbas(ic)
      do 200 jc=1,nchan
      nlmjc=nlm(jc)
      icc=nchan*(ic-1) + jc
c note indexing subtlety: we use hpvbtrn (transformed
c  <h+|V|bound> ) matrix to represent <bound|V|h-> by first
c  complex conjugating and then permuting the channel indices in
c  the S-matrix case.  In the T-matrix case we do the same, but by
c  permuting the indices and taking the imaginary part of
c  <h+|V|bound> to get <bound|V|j>.
c That is why we are using both icc and jcic in the following.
      jcic = nchan*(jc-1) + ic
      if(ibcondx.eq.0) then
      do 198 kbc=1,nbic
      do 198 ilm=1,nlmic
      do 198 jlm=1,nlmjc
198   fpfmorth(ilm,jlm,icc) = fpfmorth(ilm,jlm,icc)
     1 -ovbftrn(ilm,kbc,ic)*conjg(hpvbtrn(jlm,kbc,jcic))
      if(ioft.eq.istat)then
      do 201 kbc=1,nbic
      do 201 ilm=1,nlmic
      do 201 jlm=1,nlmjc
201   fpfmporth(ilm,jlm,icc) = fpfmporth(ilm,jlm,icc)
     1 -ovbftrn(ilm,kbc,ic)*conjg(hpvbtrnp(jlm,kbc,jcic))
      endif
      else
      do 199 kbc=1,nbic
      do 199 ilm=1,nlmic
      do 199 jlm=1,nlmjc
199   fpfmorth(ilm,jlm,icc) = fpfmorth(ilm,jlm,icc)
     1 -ovbftrn(ilm,kbc,ic)*imag(hpvbtrn(jlm,kbc,jcic))
      if(ioft.eq.istat)then
      do 202 kbc=1,nbic
      do 202 ilm=1,nlmic
      do 202 jlm=1,nlmjc
202   fpfmporth(ilm,jlm,icc) = fpfmporth(ilm,jlm,icc)
     1 -ovbftrn(ilm,kbc,ic)*imag(hpvbtrnp(jlm,kbc,jcic))
      endif
      endif
200   continue
c
c orthogonalize wrt functions in channel jc
c
c  There is no indexing subtlety in orthogonalizing wrt to
c  functions in channel jc analogous to the one in channel ic.
c  This is because we are constructing
c  <h+ already orthogonalized|V|h- orthogonalized>
c  and the channel indices are always ic,jc in every component
c  in this section.
      do 300 ic=1,nchan
      nlmic=nlm(ic)
      do 300 jc=1,nchan
      nlmjc=nlm(jc)
      nbjc=nbas(jc)
      nbsc=nscat(jc)
      istart=nbsc+1
      icc=nchan*(ic-1) + jc
      if(ibcondx.eq.0) then
      do 298 kbc=1,nbjc
      do 298 ilm=1,nlmic
      do 298 jlm=1,nlmjc
298   fpfmorth(ilm,jlm,icc) = fpfmorth(ilm,jlm,icc)
     1 -bforth(ilm,kbc,icc)*conjg(ovbftrn(jlm,kbc,jc))
      if(ioft.eq.istat.and.istart.le.nbjc)then
      do 398 kbc=istart,nbjc
      do 398 ilm=1,nlmic
      do 398 jlm=1,nlmjc
398   fpfmporth(ilm,jlm,icc) = fpfmporth(ilm,jlm,icc)
     1 -bforth(ilm,kbc,icc)*conjg(ovbftrnp(jlm,kbc,jc))
      endif
      else
      do 299 kbc=1,nbjc
      do 299 ilm=1,nlmic
      do 299 jlm=1,nlmjc
299   fpfmorth(ilm,jlm,icc) = fpfmorth(ilm,jlm,icc)
     1 -bforth(ilm,kbc,icc)*imag(ovbftrn(jlm,kbc,jc))
      if(ioft.eq.istat.and.istart.le.nbjc)then
      do 399 kbc=istart,nbjc
      do 399 ilm=1,nlmic
      do 399 jlm=1,nlmjc
399   fpfmporth(ilm,jlm,icc) = fpfmporth(ilm,jlm,icc)
     1 -bforth(ilm,kbc,icc)*imag(ovbftrnp(jlm,kbc,jc))
      endif
      endif
300   continue
c
c write it out
c
      if(iprint.ne.0) then
      write(6,107)
107   format(//,' h+(h-e)h- matrix with orthogonalized free fcns')
      do 500 ic=1,nchan
      nlmic=nlm(ic)
      do 500 jc=1,nchan
      nlmjc=nlm(jc)
      icc=nchan*(ic-1)+jc
      do 500 jlm=1,nlmjc
500   write(6,108) jlm,(fpfmorth(ilm,jlm,icc),ilm=1,nlmic)
108   format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &     (4x,3("(",f8.5,3x,f8.5,")",3x)))
      if(istat.eq.ioft)then
      write(6,109)
 109  format(//,' h+(h-e)h-p matrix with orthogonalized free fcns')
      do 501 ic=1,nchan
      nlmic=nlm(ic)
      do 501 jc=1,nchan
      nlmjc=nlm(jc)
      icc=nchan*(ic-1)+jc
      do 501 jlm=1,nlmjc
501   write(6,108) jlm,(fpfmporth(ilm,jlm,icc),ilm=1,nlmic)
      endif
      endif
      return
      end
c $Id$
