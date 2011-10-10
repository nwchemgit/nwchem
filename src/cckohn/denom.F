      subroutine denom(hdenom,nbig,hpvhp,lmtop,nstate,hpvbtrn,nbfmax
     1     ,nchnl,hbbtrn,nchan,nscat,nlm,kchan,eground,obbtrn,iprint,
     2     ntot,vopt,iopt,nsch,hpptemp,nmotot,hbbtrnd,istat,itell,ihow,
     3     isym,jsym, hbbtrne,nopen)
      implicit real*8 (a-h,o-z)
c     
c     construct denominator (h-e) from various pieces
c     
      complex*16 hdenom(nbig,nbig),hpvhp(lmtop,lmtop,nstate)
      complex*16 hpvbtrn(lmtop,nbfmax,nchnl**2)
      real*8 hbbtrn(nbfmax,nbfmax,nstate)
      real*8 hbbtrne(nbfmax,nbfmax,nstate)
      real*8 hpptemp(nbig,nbig),  vopt(nbig,nbig)
      integer nsch(nbfmax,nchan), nscat(nchnl),nlm(nchnl)
      real*8 kchan(nchnl), obbtrn(nbfmax,nbfmax,nstate)
      dimension hbbtrnd(nbfmax,nbfmax,nstate)
      character*8 istatic,istat,itell,ihow,isym,jsym
      data istatic/8hstatic  /
c     
c     note hpptemp is a real array which is EQUIVALENCED  to hdenom
c     
c     read size of Hpp and Hopt together with number of MO's in
c     target space from mesa file
c     
      if(istat.eq.istatic) go to 603
c     if(ihow.eq.itell)go to 608
c     if(istat.ne.istatic) go to 603
      write(6,194) 
 194  format(' this is not a static calculation')
c     
c     read(9) npvec,nsmall
c     
      nsmall = 7                ! number of small space 
      npvec = 9                 ! number of bound p space vectors
      write(6,877)npvec,nsmall
 877  format(/' npvec,nsmall from mesa:(hardwired for now)',2i5)
c     
c     read hpp-E in with mesa indexing
c     call rdbinsqr(hpptemp,nbig,npvec,9)
c     read hopt in with mesa indexing
c     call rdbinsqr(vopt,nbig,npvec,9)
c     add hpp-E and hopt together, still in mesa indexing
c     
c     hpptemp calculated in selci_hpp, which is called form cknew.F
c     
      write(*,*)"there is no vopt in this case, what is in hpptemp?"
      if(iopt.eq.1) then
         do ip = 1,npvec
            do jp=1,npvec
               vopt(ip,jp)=vopt(ip,jp)+hpptemp(ip,jp)
            enddo
         enddo
      else
c     
c     if iopt says to turn hopt off, load just hpp-E
c     
         do ip=1,npvec
            do jp=1,npvec
               vopt(ip,jp) = hpptemp(ip,jp)
            enddo
         enddo
      endif
      go to 605
 603  continue
      nsmall=0
      ist=0
      iprev=0
      do 606 ic=1,nchan
         nsic=nscat(ic)
         jprev=0
         do 607 jc=1,ic
            en=0.0
            if(ic.eq.jc) en=0.5*kchan(ic)**2
            ist=ist+1
            nsjc=nscat(jc)
            do isc=1,nsic
               do jsc=1,nsjc
                  iref=isc+iprev
                  jref=jsc+jprev
                  hdenom(iref,jref)=hbbtrnd(isc,jsc,ist)
     &                 -en*obbtrn(isc,jsc,ist)
                  hdenom(jref,iref)=hdenom(iref,jref)
               enddo
            enddo
            jprev=jprev+nsjc
 607     continue
         iprev=iprev+nsic
 606  continue
      if(ihow.eq.itell)then
         do i=1,iprev
            do j=1,iprev
               hpptemp(i,j)=hdenom(i,j)
            enddo
         enddo
c     call wrbinsqr(hpptemp,nbig,iprev,99)
      endif
      go to 608
 605  continue
      do i=1,nbig
         do j=1,nbig
            hdenom(i,j) = 0.0
         enddo
      enddo
c     
c     symmetry case
c     
      if(isym.eq.jsym)then
         do i=1,npvec
            do j=1,npvec
               hdenom(i,j)=vopt(i,j)
            enddo
         enddo
         iprev=npvec
         jprev=npvec
      else
c     
c     move heff=hpp+hopt - E into hdenom, reindexing to conform
c     with assignment of MO's to each channel.
c     mesa indexing will have all MO's associated with every channel
c     .
         iprev=0
         do 11 ic=1,nchan
            jprev=0
            nsic=nscat(ic)
            do 10 jc=1,ic
               nsjc=nscat(jc)
               ist=ic*(ic-1)/2+jc
c
c     This is the crucial computation of indices in mesa convention.
c     nsmall = number of MO's contributing to target configurations
c     nmotot = total number of MO's in transformation matrix
c     so (nmotot-nsmall) = number of MO's associated with every channel
c     by mesa, which does not build symmetry adapted configurations.
c     nsch(isc,ic) = index in transformation matrix of scattering orbital
c     (MO) isc in channel ic, but nsmall target MO's precede the first scattering
c     orbital in that matrix.  So we substract nsmall to get the mesa index
c     relative to the beginning of the entries for that channel.
c     
               do isc=1,nsic
                  isub=isc+iprev
                  do jsc=1,nsjc
                     jsub=jsc+jprev
                     iref=nsch(isc,ic)-nsmall+(ic-1)*(nmotot-nsmall)
                     jref=nsch(jsc,jc)-nsmall+(jc-1)*(nmotot-nsmall)
                     hdenom(isub,jsub) = vopt(iref,jref)
                     hdenom(jsub,isub) = hdenom(isub,jsub)
                 enddo
              enddo
              jprev= jprev + nsjc
 10        continue
           iprev=iprev + nsic
 11     continue
      endif
 608  continue
      istart=0
      do ic=1,nchan
         istart=istart+nscat(ic)
      enddo
      if(ihow.eq.itell)go to 609
      if(istart.ne.iprev.or.istart.ne.jprev) then
         write(6,101)
 101     format(//' indexing error in kohnopt routine hdenom')
         stop
      endif
 609  continue
c     
c     free free part
c     
      iprev=0
c.....tnr...9-94.......
c     do 31 ic=1,nchan
      do 31 ic=1,nopen
         nlmic=nlm(ic)
         jprev=0
         do 30 jc=1,ic
            nlmjc=nlm(jc)
            ist=ic*(ic-1)/2 + jc
            do ilm=1,nlmic
               isub=ilm+iprev+istart
               do jlm=1,nlmjc
                  jsub=jlm+jprev+istart
                  hdenom(isub,jsub)=hpvhp(ilm,jlm,ist)
                  hdenom(jsub,isub)=hdenom(isub,jsub)
               enddo
            enddo
            jprev=jprev+nlmjc
 30      continue
         iprev=iprev+nlmic
 31   continue
c     
c     bound-free part
c     
      iprev=istart
c.....tnr...9-94.......
c     do 41 ic=1,nchan
      do 41 ic=1,nopen
         nlmic=nlm(ic)
         jprev=0
         do 40 jc=1,nchan
            nsjc=nscat(jc)
            icc=nchan*(ic-1) + jc
            do ilm=1,nlmic
               isub=ilm+iprev
               do jsc=1,nsjc
                  jsub=jsc+jprev
                  hdenom(isub,jsub) = hpvbtrn(ilm,jsc,icc)
                  hdenom(jsub,isub) =hdenom(isub,jsub)
               enddo
            enddo
            jprev=jprev+nsjc
 40      continue
         iprev=iprev+nlmic
 41   continue
      ntot=iprev
      if(ihow.eq.itell)then
         iw=istart+1
         if(isym.ne.jsym)then
            iw=1
         endif
c     
c     do not write hdenom
c     
c     write(7)((hdenom(i,j),i=1,ntot),j=iw,ntot)
      endif
c
      if(iprint.ne.0) then
         write(6,107)
 107     format(//' denominator matrix of (h-e)')
         do 60 i=1,ntot
            write(6,108) i,(hdenom(j,i),j=1,ntot)
 60      continue
 108     format(1x,i3,3("(",f8.5,3x,f8.5,")",3x),/,
     &        (4x,3("(",f8.5,3x,f8.5,")",3x)))
      endif
      return
      end
      

c $Id$
