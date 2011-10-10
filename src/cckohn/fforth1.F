      subroutine fforth1(ovbftrn,lmtop,nbfmax,nchnl,nstate,
     1 hpvbtrn,bforth,nchan,nlm,nbas,fpfporth,iprint)
       implicit real*8 (a-h,o-z)
c
c construct free-free matrix element:
c h+(orthogonalized) (h-e) h+(orthogonalized)
c
 
      complex*16 fpfporth(lmtop,lmtop,nstate)
      complex*16 ovbftrn(lmtop,nbfmax,nchnl)
      complex*16 hpvbtrn(lmtop,nbfmax,nchnl**2)
      complex*16 bforth(lmtop,nbfmax,nchnl**2)
      integer nbas(nchnl),nlm(nchnl)
c
c orthogonalize wrt fcns in channel ic
c
      do 200 ic=1,nchan
      nlmic=nlm(ic)
      nbic=nbas(ic)
      do 200 jc=1,ic
      nlmjc=nlm(jc)
      ist=ic*(ic-1)/2 + jc
c Note indexing subtlety:  channel indices for hpvbtrn must
c be permuted here because we are using <h+|V|bound> to represent <bound|V|h+>
      jcic = nchan*(jc-1) + ic
      do 200 kbc=1,nbic
      do 200 ilm=1,nlmic
      do 200 jlm=1,nlmjc
200   fpfporth(ilm,jlm,ist) = fpfporth(ilm,jlm,ist)
     1 -ovbftrn(ilm,kbc,ic)*hpvbtrn(jlm,kbc,jcic)
c
c orthogonalize wrt functions in channel jc
c
      do 300 ic=1,nchan
      nlmic=nlm(ic)
      do 300 jc=1,ic
      nlmjc=nlm(jc)
      nbjc=nbas(jc)
      ist=ic*(ic-1)/2 + jc
      icc=nchan*(ic-1) + jc
      do 300 kbc=1,nbjc
      do 300 ilm=1,nlmic
      do 300 jlm=1,nlmjc
300   fpfporth(ilm,jlm,ist) = fpfporth(ilm,jlm,ist)
     1 -bforth(ilm,kbc,icc)*ovbftrn(jlm,kbc,jc)
c
c write it out
c
      if(iprint.ne.0) then
      write(6,107)
107   format(//,' h+(h-e)h+ matrix with orthogonalized free fcns')
      do 500 ic=1,nchan
      nlmic=nlm(ic)
      do 500 jc=1,ic
      nlmjc=nlm(jc)
      ist=ic*(ic-1)/2 + jc
      do 500 jlm=1,nlmjc
500   write(6,108) jlm,(fpfporth(ilm,jlm,ist),ilm=1,nlmic)
108   format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &     (4x,3("(",f8.5,3x,f8.5,")",3x)))
      endif
      return
      end
c $Id$
