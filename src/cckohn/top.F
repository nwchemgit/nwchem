      subroutine top(htop,nbig,nsmall,hpvhm,lmtop,nstate,hpvbtrn,
     1     nbfmax,nchnl,nchan,nscat,nlm,iprint,ntot,nfree,ibcondx
     2     ,htopp,hpvhmp,hpvbtrnp,itell,ihow,ovbftrnp,istat,ioft,cov,
     $     nbas,nopen)
      implicit real*8 (a-h,o-z)
c
c     construct numerator (h-e) from various pieces
c     
      complex*16 htop(nbig,nsmall),hpvhm(lmtop,lmtop,nchnl**2)
      complex*16 cov(nbig,nsmall)
      complex*16 htopp(nbig,nsmall),hpvhmp(lmtop,lmtop,nchnl**2)
      complex*16 hpvbtrn(lmtop,nbfmax,nchnl**2)
      complex*16 hpvbtrnp(lmtop,nbfmax,nchnl**2)
      complex*16 ovbftrnp(lmtop,nbfmax,nchnl)
      integer nscat(nchnl),nlm(nchnl),nbas(nchnl)
      character*8 itell,ihow,istat,ioft
      istart=0
      do 21 ic=1,nchan
21    istart=istart+nscat(ic)
c
c     free free part
c
      iprev=0
c     do 31 ic=1,nchan
c.....tnr...9-94.......
      do 31 ic=1,nopen
         nlmic=nlm(ic)
         jprev=0
         do 30 jc=1,nopen
            nlmjc=nlm(jc)
            ist=nchan*(ic-1) + jc
            do ilm=1,nlmic
               do jlm=1,nlmjc
                  isub=ilm+iprev+istart
                  jsub=jlm+jprev
                  htop(isub,jsub)=hpvhm(ilm,jlm,ist)
                  htopp(isub,jsub)=hpvhmp(ilm,jlm,ist)
               enddo
            enddo
            jprev=jprev+nlmjc
 30      continue
         iprev=iprev+nlmic
 31   continue
c
c     new nall determination  ..tnr 6/20/94
c     
      nall=isub
c     
c     bound-free part
c     
      iprev=0
      do 41 ic=1,nopen
         nlmic=nlm(ic)
         jprev=0
         do 40 jc=1,nchan
            nsjc=nscat(jc)
            icc=nchan*(ic-1) + jc
            do 39 ilm=1,nlmic
               if(ibcondx.eq.0) then
                  do 37 jsc=1,nsjc
                     isub=ilm+iprev
                     jsub=jsc+jprev
                     htopp(jsub,isub) = conjg( hpvbtrnp(ilm,jsc,icc))
                     htop(jsub,isub) = conjg( hpvbtrn(ilm,jsc,icc))
 37               continue
               else
                  do 38 jsc=1,nsjc
                     isub=ilm+iprev
                     jsub=jsc+jprev
                     htopp(jsub,isub) = imag( hpvbtrnp(ilm,jsc,icc))
                     htop(jsub,isub) = imag( hpvbtrn(ilm,jsc,icc))
 38               continue
               endif
 39         continue
            jprev=jprev+nsjc
 40      continue
         iprev=iprev+nlmic
 41   continue
c     ... tnr 6/20/94
      nfree=iprev
c     nball=jprev
c     
      if(ihow.eq.itell)then
c     ... tnr 6/20/94
c     nall=nfree+jprev do not write (7)
c
c     write(7)((htopp(i,j),i=1,nall),j=1,nfree)
c     
c     load and write out bound-free overlap matrix
c     
         do i=1,nbig
            do j=1,nfree
               cov(i,j)=0.
            enddo
         enddo
         if(istat.eq.ioft)then
            nlmprev=0
            nbprev=0
c.....tnr...9-94.......
c     do 11 i=1,nchan
            do 11 i=1,nopen
               nlmic=nlm(i)
               nsic=nbas(i)      
               do 12 ilm=1,nlmic
                  ill=ilm+nlmprev
                  if(ibcondx.eq.0) then
                     do 13 ib=1,nsic
                        ibb=ib+nbprev
                        cov(ibb,ill)=conjg(ovbftrnp(ilm,ib,i))
 13                  continue
                  else
                     do 14 ib=1,nsic
                        ibb=ib+nbprev
                        cov(ibb,ill)=imag(ovbftrnp(ilm,ib,i))
 14                  continue
                  endif
 12            continue
               nlmprev=nlmprev+nlmic
               nbprev=nbprev+nsic
 11         continue
         endif
c         write(7)nbprev
c         write(7)((cov(i,j),i=1,nbprev),j=1,nfree)
c         write(7)((htop(i,j),i=1,nall),j=1,nfree)
      endif
c
      if(iprint.ne.0) then
         write(6,107)
 107     format(//' numerator matrix of (h-e)')
         do 60 i=1,nfree
            write(6,108) i,(htop(j,i),j=1,ntot)
 60      continue
 108     format(1x,i3,6f12.5,/,(4x,6f12.5))
      endif
      return
      end
c $Id$
