      subroutine cknew(ovbf,hpvb,hpvhp,hpvhm)
      implicit none
c     
c     parameters, for explanation see cckohn.F
c     
      integer ltop, mtop, lmtop, nchnl, nstate, nsymst, nbfmax, nbfcmx,
     &     maxene
      parameter (ltop=6, mtop=4, lmtop=43, nchnl=20, nstate=nchnl**2, 
     &     nsymst=nchnl*(nchnl+1)/2, nbfmax=130, nbfcmx=130, maxene=20)
      double precision echan(nchnl),energy(maxene),kchan(nchnl)
c     
      double complex 
     1     ovbf(lmtop,nbfcmx,nchnl),
     2     hpvb(lmtop,nbfcmx,nstate),
     3     hpvhp(lmtop,lmtop,nsymst),
     4     hpvhm(lmtop,lmtop,nstate)
c
c*****************modified april1995********************************************
c     includes pseudopotential option
c     
c*****************modified september1994*****************************************
c     modified to handle closed channels explicitly
c     for closed channels, only L-2 functions appear and free functions are dropped
c
c*****************modified august 1992*****************************************
c incorporates either symmetry in denom or standard no symmetry runs in mesa
c also has an off-shell option and a static option
c
c   set istat = exchange   for an exact exchange run
c             = static     for a static run
c             = offshell   for an off-shell run
c             = other      for a standard run
c
c   set ihow  = mesa       for a mesa setup, with no solution of linear equations
c             = other      for a full run
c
c   set isym  = symmetry   if mesa was run in symmetry
c             = other      if mesa was not run in symmetry
c   
c*************modified to handle explicit inclusion of exchange*****************
c complex kohn scattering code started  8-8-88
c
c  reads bound-bound, bound-free, and free-free
c  matrices from the sequence of codes:
c newgrid, ass1, bespline, basis, statnu, ffbf, bndnew2, optpot
c
c changed 5/1/90 to read nmotot BHL
c
c changed to read in channel-dependent bound-bound matrices 12/21/88
c
c changed 12/22/88 to handle multi-channel bound-bound matrices
c correctly
c
c  addition of option to orthogonalize free functions to bound fcns
c  8-10-88
c
c  addition of eigenphase computation 10-3-88
c  note this uses routine cgeev for complex diagonalization from
c  SLATEC library
c
c optical potential added 1/19/89
c
c debugging for 3 channel case:  fixed fforth1 indexing 3/18/89
c                                fixed fforth2 indexing 3/18/89
c
c added total cross section  for multichannel case 3/18/89
c
c  ### MESA interface implemented 7/12/89 ####
c  reads file bndmat created my MESA code m950
c
      parameter (nblok=1000)
      parameter (nbig=700)
      parameter (nsmall=nchnl*lmtop)
c
      integer fmt(10),iclosed(nchnl)
      integer nlm(nchnl)
      integer ngauss(nchnl),ngch(nbfmax,nchnl)
      integer nscat(nchnl),nsch(nbfmax,nchnl)
      integer lch(lmtop,nchnl),mch(lmtop,nchnl)
      integer ntarg(nchnl),ntch(nbfmax,nchnl)
      integer nbas(nchnl),nbch(nbfmax,nchnl),kpvt(nbig)
      character*1 jobvl, jobvr
      double precision kchan(nchnl),kchanp(nchnl)
      double precision scratch(nbfmax,nbfmax,nstate)
      double complex unity(nsmall,nsmall)
      double precision tr(nbfmax,nbfmax)
      double precision hbb(nbfmax**2),hbbe(nbfmax**2)
      double precision hbbtrnd(nbfmax,nbfmax,nstate),hbbtrne(nbfmax,nbfmax,nstate)
      double precision obbtrn(nbfmax,nbfmax,nstate),hbbtrn(nbfmax,nbfmax,nstate)
      double precision vopt(nbig,nbig)
      double precision xsecmat(nchnl,nchnl)
      double precision hpptemp(nbig,nbig)
      double complex ovbf(lmtop,nbfmax,nchnl),ovbfp(lmtop,nbfmax,nchnl)
      double complex hpvb(lmtop,nbfmax,nchnl**2)
      double complex hpvbp(lmtop,nbfmax,nchnl**2)
      double complex hpvbtrn(lmtop,nbfmax,nchnl**2)
      double complex hpvbtrnp(lmtop,nbfmax,nchnl**2)
      double complex ovbftrn(lmtop,nbfmax,nchnl)
      double complex ovbftrnp(lmtop,nbfmax,nchnl)
      double complex hpvhp(lmtop,lmtop,nstate)
      double complex hpvhm(lmtop,lmtop,nchnl**2)
      double complex hpvhmp(lmtop,lmtop,nchnl**2)
      double complex hdenom(nbig,nbig),cov(nbig,nsmall)
      double complex htop(nbig,nsmall)
      double complex bforth(lmtop,nbfmax,nchnl**2)
      double complex bforthp(lmtop,nbfmax,nchnl**2)
      double complex hpvhmps(lmtop,lmtop,nchnl)
      double complex zdotc,zdotu,ai, det, htop2(nbig,nsmall)
      double complex smat(nsmall,nsmall),work(nbig)
      double complex seig(nsmall),svec(nsmall,nsmall)
      double complex swork(2*nsmall)
      double precision srwork(2*nsmall)
      equivalence (vopt,hbb)
      equivalence (unity, hdenom),(hbbtrne,hbbe)
      equivalence (hpptemp,hdenom)
      equivalence (ovbf, hpvb, hpvhm)
      equivalence (ovbfp, hpvbp, hpvhmp)
      equivalence (hpvhp,scratch),(hpvbtrnp,cov)
c      equivalence (fpfmorth,hdenom),(fpfporth,smat),(bforth,htop)
      equivalence (bforth,htop),(smat,ovbftrn)
     1 ,(bforthp,htop2)
c
c     variables for options
c      
      integer iprint, iorth, iopt,ipseudo
      character*8 ioft,iex,itell,jsym,istatic,ihow,istat,isym
      data ioft /"offshell"/
      data iex/"exchange"/, itell/"mesa    "/, jsym/"symmetry"/,
     &     istatic/"static  "/
c
c..unicos
c      call link("unit8=intsbf,unit5=inkohn,unit6=(outkohn,hc,create),
c     1 unit10=intsff,unit9=bndmat
c     2 ,unit7=(punkohn,create)
c     3 //")
c..unicos
c
      ai=cmplx(0.,1.)
czz
czz      open(5,file='inkohn',status='unknown')
czz      open(6,file='outkohn',status='unknown')
czz      read(5,*)iprint,iorth,iopt,ipseudo
czz      read(5,734) istat
czz
      write(*,*)"options are hardwired, 
     &     iprint=1, iorth=1, iopt=0, ipseudo=0, istat=nex, ihow =bmesa"
      iprint = 1
      iorth=1
      iopt=0
      ipseudo=0
      istat="nex"
      ihow = "bmesa"
      isym = "symmetry" 
      ibcondx=1                 ! tmatrix boundary condition
      write(*,"(2a8)")"istat=  ", istat
c
      if(istat.eq.ioft)then     ! ioft = "offshell"
         open(77,file='intsffof',form='unformatted',status='old') 
         open(66,file='intsbfof',form='unformatted',status='old') 
      endif
      if(istat.eq.iex.or.istat.eq.ioft)then ! iex= exchange
         open(15,file='exchmat',form='unformatted',status='old')
      endif
      if(iopt.eq.1)then
         write(6,595)
      endif
      if(ipseudo.ne.0)then
         write(6,596)
      else
         write(6,597)
      endif
c     
c     read(5,734)ihow
c     
      if(ihow.eq.itell)then     ! itell = "mesa    "
         open(7,file='bmesa',form='unformatted',status='unknown')
         if(istat.eq.istatic)then ! istatic = static
            open(99,file='bstatic',form='unformatted',status='unknown')
         endif
      else
         write(*,*)"input originally from files tmats=88 and punkohn=7"
c     open(88,file='tmats',status='unknown')
c     open(7,file='punkohn',form='unformatted',status='unknown')
      endif
c     
c     read(5,734)isym
c     write(6,734)isym
c     
      
      write(*,*)"files intsbf, bndmat,intsff eliminated"
c     
c      open(8,file='intsbf',form='unformatted',status='old') 
c      open(9,file='bndmat',form='unformatted',status='old')
c      open(10,file='intsff',form='unformatted',status='old') 
c
      if(ipseudo.ne.0)then
         open(11,file='intsbfps',form='unformatted',status='old') 
         open(12,file='intsffps',form='unformatted',status='old') 
      endif
c     
c     read(10) ibcondx
c     
      if(ibcondx.eq.0) then     ! smatrix
         write(6,731)
      else if(ibcondx.eq.1) then
         write(6,732)
      else
         write(6,733)
         stop
      endif
c     
c      read(8) nener,nchan,(nlm(ic),ic=1,nchan)
c
      if(nchan.gt.nchnl) then
         write(6,*)' nchan gt nchnl ',nchan,nchnl
         stop 11
      end if
c
c      read(8) ((lch(j,ic),mch(j,ic),j=1,nlm(ic)),ic=1,nchan)
c      read(8) eground
c      read(8) (ngauss(ic),ic=1,nchan)
c      read(8) ((ngch(j,ic),j=1,ngauss(ic)),ic=1,nchan)
c
      nbmax=0
      nlmax=0
      do 851 i=1,nchan
         if(nlm(i).gt.nlmax)nlmax=nlm(i)
         if(ngauss(i).gt.nbmax)nbmax=ngauss(i)
 851  continue
      if(nlmax.gt.lmtop)then
         write(6,*)' nlmax gt lmtop', nlmax,lmtop
         stop
      endif
      if(nbmax.gt.nbfmax)then
         write(6,*)' nbmax gt nbfmax', nbmax,nbfmax
         stop
      endif
c
c      read(10)
c      read(10)
c      read(10)
c
      if(ipseudo.ne.0)then
         read(11)
         read(11)
         read(11)
         read(11)
         read(11)
         read(12)
         read(12)
         read(12)
         read(12)
      endif
c     
c     read assignment of scattering basis to channels
c     
      read(5,*) (nscat(ic),ic=1,nchan)
      write(6,138) (nscat(ic),ic=1,nchan)
      read(5,*) ((nsch(j,ic),j=1,nscat(ic)),ic=1,nchan)
      do 875 ic=1,nchan
         write(6,139) ic, (nsch(j,ic),j=1,nscat(ic))
 875  continue
c     
c     read assignment of target orbitals to channels
c     
      if(iorth.ne.0) then
         read(5,*) (ntarg(ic),ic=1,nchan)
         itarg=0
         write(6,137) (ntarg(ic),ic=1,nchan)
         do  i=1,nchan
            itarg = itarg + ntarg(i)
         enddo
         if(itarg.ne.0) then
            read(5,*) ((ntch(j,ic),j=1,ntarg(ic)),ic=1,nchan)
            do 876 ic=1,nchan
               write(6,139) ic, (ntch(j,ic),j=1,ntarg(ic))
 876        continue
         endif
      else
         do ic=1,nchan
            ntarg(ic) = 0
         enddo
      endif
c     
c     construct total orbital list for orthogonalization of free fcns.
c     if iorth was zero this list is the same as nscat list
c     
      nntot=0
      nbtot=0
      mfree=0
      do 4 ic=1,nchan
         nntot=nntot+nlm(ic)+nscat(ic)    
         nbtot=nbtot+nscat(ic)
         mfree=mfree+nlm(ic)
         ntic=ntarg(ic)
         nsic=nscat(ic)
         nbas(ic)=ntic+nsic
         do  i=1,nsic
            nbch(i,ic)=nsch(i,ic)
         enddo
         if (ntic.eq.0.or.itarg.eq.0) go to 4
         do  i=1,ntic
            nbch(nsic+i,ic) = ntch(i,ic)
         enddo
 4    continue
      if(nntot.gt.nbig)then
         write(6,*)' stopping because ntot gt nbig', nntot,nbig
         stop
      else
         write(6,*)' total dimension of hamiltonian is ',nntot
      endif
c     
      if(ihow.eq.itell)then ! itell = mesa
         write(7)nchan,nbtot,mfree,(nlm(i),i=1,nchan)
      endif
      if(iorth.ne.0) then
         write(6,787)
         do ic=1,nchan
            write(6,789) ic,(nbch(i,ic),i=1,nbas(ic))
         enddo
         write(6,791)
      endif
c     
c     skip some preliminary info on file bndmat made by m950
c     
c      read(9)
c      read(9)
c      read(9)
c      read(9)
c
c
c     istat = nex
c
      if(istat.eq.iex.or.istat.eq.ioft)then ! iex = exchange ioft = offshell
         read(15)
      endif
c
c     read in transformation matrix
c     
c..   bhl 5/1/90
c     read(9) nbftot,nmotot
c     
      if(nbftot.gt.nbfmax) then
         write(6,*)' nbftot gt nbfmax ',nbftot,nbfmax
         stop 11
      end if
c..   bhl      nmotot=nbftot
      nsy = nbftot
c
c      call rdbinsqr(tr,nbfmax,nbftot,9)
c
      write(6,99)
      do 1 i=1,nmotot
         write(6,100)i,(tr(j,i),j=1,nbftot)
 1    continue
c
c     read in bound-bound overlap matrix
c     
      call rdbintri(hbb,nsy,9)
c     
c     transform overlap matrix as a test
c     
      call ovtrans(hbb,obbtrn,nsy,nchan,nbas,nbch,tr,scratch,nbfmax
     1     ,nstate,iprint)
c     
c     read in energy-independent bound-bound direct hamiltonian
c     and transform hbb direct
c     
      call bbtrans(hbb,hbbtrnd,nsy,nchan,nbas,nbch,tr,scratch,nbfmax
     1     ,nstate,iprint,hbbe,hbbtrne,istat)
c     
c     open a loop on energies
c     
      if(ihow.eq.itell)then     ! itell = mesa
         write(7)nener
      endif
      do 1000 iene=1,nener
c
c         read(8) (kchan(ic),ic=1,nchan)
c         read(8) (iclosed(ic),ic=1,nchan)
c         read(10)
c
         if(ipseudo.ne.0)then
            read(11)
            read(12)
         endif
         write(6,111) (kchan(i),i=1,nchan)
c     
c     compute the number open open channels at this energy
c     
         nopen=nchan
         do  iopen=1,nchan
            nopen=nopen-iclosed(iopen)
         enddo
         if(ihow.eq.itell)then  ! itell = mesa
            write(7)(kchan(ic),ic=1,nchan)
         endif
         if(istat.eq.ioft)then  ! ioft = offshell
            read(66)(kchanp(i),i=1,nchan)
            read(77)(kchanp(i),i=1,nchan)
            write(6,112)(kchanp(i),i=1,nchan)
         endif
         do 507 ic=1,nchan
            nlmic=nlm(ic)
            ngic=ngauss(ic)
c
c            read(8)((ovbf(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
c
            if(istat.eq.ioft)then ! ioft = offshell
               read(66)((ovbfp(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
            endif
 507     continue
c     
c     transform the bound-free overlap matrices
c     
         call obftrans(ovbf,ovbftrn,tr,nsy,nchan,nbas,nbch,nlm,ngauss,
     1        ngch,nbfmax,lmtop,nchnl,iprint,ovbftrnp,ovbfp,istat,ioft)
         do ic=1,nchan
            do jc=1,nchan
               ist=nchan*(ic-1)+jc
               nlmic=nlm(ic)
               ngjc=ngauss(jc)
c               
c               read(8)((hpvb(ilm,ig,ist),ilm=1,nlmic),ig=1,ngjc)
c
               if(istat.eq.ioft)then ! ioft = offshell
                  read(66)((hpvbp(ilm,ig,ist),ilm=1,nlmic),ig=1,ngjc)
               endif
            enddo
         enddo
c     
c     add pseudopotential terms to bound-free integrals
c     
         if(ipseudo.ne.0)then
            do ic=1,nchan
               nlmic=nlm(ic)
               ngic=ngauss(ic)
c     
c     read(11)((hpvbtrn(ilm,ig,ic),ilm=1,nlmic),ig=1,ngic)
c     
               ist=nchan*(ic-1)+ic
               do ilm=1,nlmic
                  do ig=1,ngic
                     hpvb(ilm,ig,ist)=
     &                    hpvb(ilm,ig,ist)+hpvbtrn(ilm,ig,ic)
                  enddo
               enddo
            enddo
         endif
c
c     transform the bound-free (h-e) matrices
c     
         call hbftrans(hpvb,hpvbtrn,tr,nsy,nchan,nbas,nbch,nlm,ngauss,
     &        ngch,nbfmax,lmtop,nchnl,iprint,hpvbp,hpvbtrnp,istat,ioft)
c     
c     read the free free integrals file
c     
         do ic=1,nchan
            do jc=1,ic
               ist = ic*(ic-1)/2 +jc
               nlmic=nlm(ic)
               nlmjc=nlm(jc)
c     
c     read(10) ((hpvhp(ilm,jlm,ist),ilm=1,nlmic),jlm=1,nlmjc)
c     
            enddo
         enddo
c     
c     add pseudopotential terms to free-free integrals
c     
         if(ipseudo.ne.0)then
            do ic=1,nchan
               nlmic=nlm(ic)
c     
c     read(12) ((hpvhm(ilm,jlm,ic),ilm=1,nlmic),jlm=1,nlmic)
c     
               ist=ic*(ic+1)/2
               do ilm=1,nlmic
                  do jlm=1,nlmic
                     hpvhp(ilm,jlm,ist)=hpvhp(ilm,jlm,ist)+
     $                    hpvhm(ilm,jlm,ic)
                  enddo
               enddo
            enddo
         endif
         do ic=1,nchan
            do jc=1,nchan
               ist=nchan*(ic-1)+jc
               nlmic=nlm(ic)
               nlmjc=nlm(jc)
               read(10) ((hpvhm(ilm,jlm,ist),ilm=1,nlmic),jlm=1,nlmjc)
               if(istat.eq.ioft)then ! ioft = offshell
                  read(77) ((hpvhmp(ilm,jlm,ist),ilm=1,nlmic),
     &                 jlm=1,nlmjc)
               endif
            enddo
         enddo
c     
c     add pseudopotential terms to free-free integrals
c     
         if(ipseudo.ne.0)then
            do ic=1,nchan
               nlmic=nlm(ic)
               read(12) ((hpvhmps(ilm,jlm,ic),ilm=1,nlmic),jlm=1,nlmic)
               ist=nchan*(ic-1)+ic
               do ilm=1,nlmic
                  do jlm=1,nlmic
                     hpvhm(ilm,jlm,ist)=hpvhm(ilm,jlm,ist)+
     $                    hpvhmps(ilm,jlm,ic)
                  enddo
               enddo
            enddo
         endif
c
c     ***orthogonalization section***
c     ***must use for static exchange***
c     
c     construct bound-free (h-e) matrices for free fcns orthogonalized
c     bound functions
c     
         if(iorth.ne.0) then    ! orthogonalize
            call bforthog(ovbftrn,lmtop,nbfmax,nchnl,nstate,hpvbtrn,
     1           hbbtrnd,nchan,nlm,nbas,nbas,bforth,iprint,bforthp,
     2           obbtrn,kchan,eground,ovbftrnp,hbbtrne,hpvbtrnp,istat,
     3           ioft,nscat)
c     
c     construct free-free (h-e) matrices for free functions orthogonalized
c     to bound functions
c     
            call fforth1(ovbftrn,lmtop,nbfmax,nchnl,nstate,
     1           hpvbtrn,bforth,nchan,nlm,nbas,hpvhp,iprint)
            call fforth2(ovbftrn,lmtop,nbfmax,nchnl,nstate,
     1           hpvbtrn,bforth,nchan,nlm,nbas,hpvhm,iprint,
     2           ibcondx,hpvbtrnp,hpvhmp,istat,ioft,nscat,ovbftrnp)
c     
c     load matrices with orthogonalized free functions back into
c     the arrays which held the same matrices with unorthogonalized free functions
c     so that remainder of the code is independent of orthogonalization procedure
c     
            if(istat.ne.ioft)then ! ioft = offshell
               do ic=1,nchan
                  nlmic=nlm(ic)
                  do jc=1,nchan
                     nlmjc=nlm(jc)
                     icc=nchan*(ic-1)+jc
                     do ilm=1,nlmic
                        do jlm=1,nlmjc
                           hpvhmp(ilm,jlm,icc)=hpvhm(ilm,jlm,icc)
                        enddo
                     enddo
                  enddo
               enddo
            endif
            do ic=1,nchan
               nlmic=nlm(ic)
               do jc=1,nchan
                  nsjc=nscat(jc)
                  icc=nchan*(ic-1)+jc
                  do ilm=1,nlmic
                     do jsc=1,nsjc
                        if(istat.ne.ioft)then ! ioft = offshell
                           hpvbtrnp(ilm,jsc,icc)=bforth(ilm,jsc,icc)
                        else
                           hpvbtrnp(ilm,jsc,icc)=bforthp(ilm,jsc,icc)
                        endif
                        hpvbtrn(ilm,jsc,icc)=bforth(ilm,jsc,icc)
                     enddo
                  enddo
               enddo
            enddo     
         else                   ! do not orthogonalize
            if(istat.ne.ioft)then ! ioft = offshell
               do ic=1,nchan
                  nlmic=nlm(ic)
                  do jc=1,nchan
                     nlmjc=nlm(jc)
                     icc=nchan*(ic-1)+jc
                     do ilm=1,nlmic
                        do jlm=1,nlmjc
                           hpvhmp(ilm,jlm,icc)=hpvhm(ilm,jlm,icc)
                        enddo
                     enddo
                  enddo
               enddo
               do ic=1,nchan
                  nlmic=nlm(ic)
                  do jc=1,nchan
                     nsjc=nscat(jc)
                     icc=nchan*(ic-1)+jc
                     do ilm=1,nlmic
                        do jsc=1,nsjc
                           hpvbtrnp(ilm,jsc,icc)=hpvbtrn(ilm,jsc,icc)
                        enddo
                     enddo
                  enddo
               enddo
            endif
         endif
c     
c     *** end of orthogonalization section ***
c     
c     construct denominator
c     
         call denom(hdenom,nbig,hpvhp,lmtop,nstate,hpvbtrn,nbfmax,
     1        nchnl,hbbtrn,nchan,nscat,nlm,kchan,eground,obbtrn,iprint,
     2        ntot,vopt,iopt,nsch,hpptemp,nmotot,hbbtrnd,istat,itell,
     3        ihow,isym,jsym,hbbtrne,nopen)
c
c     cnstruct numerator
c     
         call top(htop,nbig,nsmall,hpvhm,lmtop,nstate,hpvbtrn,nbfmax,
     1        nchnl,nchan,nscat,nlm,iprint,ntot,nfree,ibcondx,htop2,
     2        hpvhmp,hpvbtrnp,itell,ihow,ovbftrnp,istat,ioft,cov,
     3        nbas,nopen)
         if(mfree.ne.nfree)then
            write(6,*)' note that nfree and mfree are different', 
     &           mfree,nfree
         endif
c     
c     if this is a mesa run, write out born terms and quit
c     
         if(itell.eq.ihow)then ! itell = mesa
            iprev = 0
            do 651 ic=1,nopen
               nlmic=nlm(ic)
               jprev = 0
               do 650 jc=1,nopen
                  nlmjc = nlm(jc)
                  ii=max0(ic,jc)
                  jj=ic+jc-ii
                  ist=ii*(ii-1)/2 + jj
                  icc=nchan*(ic-1) + jc
                  do 640 ilm=1,nlmic
                     isub=ilm+iprev
                     if(ibcondx.eq.0) then ! smatrix boundary condition
c
c     ******note**************
c     this is incorrect for off-shell s-matrix 
c     because hpvhpp was never calculated
c     
                        do jlm=1,nlmjc
                           jsub = jlm+jprev
                           smat(isub,jsub) =conjg(hpvhp(ilm,jlm,ist))
                        enddo
                     else       ! tmatrix boundary condition
                        do jlm=1,nlmjc
                           jsub = jlm+jprev
                           smat(isub,jsub) =imag(hpvhmp(ilm,jlm,icc))
                        enddo
                     endif
 640              continue
                  jprev=jprev+nlmjc
 650           continue
               iprev=iprev+nlmic
 651        continue
            write(7)((smat(i,j),i=1,nfree),j=1,nfree)
         endif
         if(ihow.eq.itell)go to 1000 ! itell = mesa
c     
c     save numerator matrix
c     
         nterm=nbig*nsmall
c
c     call ccopy(nterm,htop,1,htop2,1)
c     
c     solve linear equations
c     
c     call matinv(hdenom,ntot,htop2,nfree,det,nbig)
         call cgefs(hdenom,nbig,ntot,htop2(1,1),1,ind,work,kpvt)
         write(6,*)' solving linear equations'
         do i=2,nfree
            call cgefs(hdenom,nbig,ntot,htop2(1,i),2,ind,work,kpvt)
         enddo
c     
         if(iprint.ne.0) then
            write(6,107)
            do  i=1,ntot
               write(6,108) i,(htop2(i,j),j=1,nfree)
            enddo
         endif
         do ilm=1,nfree
            do jlm=1,nfree
               smat(ilm,jlm)=zdotu(ntot,htop2(1,jlm),1,htop(1,ilm),1)
            enddo
         enddo
         if(iprint.ne.0) then
            write(6,109)
            do i=1,nfree
               write(6,108) i,(smat(i,j),j=1,nfree)
            enddo
         endif
c     
c     add born terms
c     
         iprev = 0
         do 551 ic=1,nopen
            nlmic=nlm(ic)
            jprev = 0
            do 550 jc=1,nopen
               nlmjc = nlm(jc)
               ii=max0(ic,jc)
               jj=ic+jc-ii
               ist=ii*(ii-1)/2 + jj
               icc=nchan*(ic-1) + jc
               do 540 ilm=1,nlmic
                  isub=ilm+iprev
                  if(ibcondx.eq.0) then ! smatrix
c     
c     ******note**************
c     this is incorrect for off-shell s-matrix because hpvhpp 
c     was never calculated
c     
                     do jlm=1,nlmjc
                        jsub = jlm+jprev
                        smat(isub,jsub) =ai*(-smat(isub,jsub)+
     &                       conjg(hpvhp(ilm,jlm,ist)))
                     enddo
                  else
                     do jlm=1,nlmjc
                        jsub = jlm+jprev
                        smat(isub,jsub) =-2.0*(-smat(isub,jsub)+
     &                       imag(hpvhmp(ilm,jlm,icc)))
                     enddo
                  endif
 540           continue
               jprev=jprev+nlmjc
 550        continue
            iprev=iprev+nlmic
 551     continue
c     
c     write out s matrix
c     
         if(ibcondx.eq.0) then  ! smatrix
            write(6,776) (kchan(i),i=1,nopen)
         else                   ! tmatrix
c     
c     Here we have the T-matrix.  Compute the cross sections from it.
c     
            istart = 0
            do 879 ic=1,nopen
               jstart = 0
               ni = nlm(ic)
               do 878 jc=1,nopen
                  nj=nlm(jc)
                  ist=istart+1
                  ifin=istart+ni
                  jst=jstart+1
                  jfin=jstart+nj
                  write(7)ic,jc,ni,nj,kchan(ic),kchan(jc)
                  write(88,177)ic,jc,ni,nj,kchan(ic),kchan(jc)
                  write(7)((smat(ii,jj),ii=ist,ifin),jj=jst,jfin)
                  write(88,277)((smat(ii,jj),ii=ist,ifin),jj=jst,jfin)
                  summod = 0.
                  do ilm=1,ni
                     do jlm=1,nj
                        isub = istart + ilm
                        jsub = jstart + jlm
                        summod = summod + abs(smat(isub,jsub))**2
                     enddo
                  enddo
                  xsecmat(ic,jc) = 4.0*3.141592654*summod/kchan(ic)**2
                  jstart = jstart + nlm(jc)
 878           continue
               istart = istart + nlm(ic)
 879        continue
            write(6,874)
            write(6,871) (i,i=1,nopen)
            do i=1,nopen
               write(6,872) i, (xsecmat(i,j),j=1,nopen)
            enddo
c     
            write(6,778) (kchan(i),i=1,nopen)
         endif
         do i=1,nfree
            write(6,777) i,(smat(j,i),j=1,nfree)
         enddo
c
c     test unitarity
c     
c     watch out; unity(i,j) is equivalenced to denom(i,j)
c     
         if(ibcondx.ne.0) then
            do i=1,nfree
               do j=1,nfree
                  smat(i,j) =  (0.0,2.0)*smat(i,j)
               enddo
            enddo
            do i=1,nfree
               smat(i,i) = smat(i,i) + 1.0
            enddo
         endif
c     
c     we now own the S-matrix -- no matter what type of boundary
c     we used in the calculation
c     first test unitarity and then compute eigenphases
c     
         do i=1,nfree
            do j=1,nfree
               unity(i,j)=zdotc(nfree,smat(1,i),1,smat(1,j),1)
            enddo
         enddo
         write(6,793)
         do i=1,nfree
            write(6,777) i,(unity(j,i),j=1,nfree)
         enddo
c
c     job = 0
c     call cgeev(smat,nsmall,nfree,seig,svec,nsmall,swork,job,inform)
c     
         write(6,*)' S-matrix'
         do i=1,nfree
            write(6,777) i,(smat(j,i),j=1,nfree)
         enddo
         jobvl='n'
         jobvr='v'
         lwork=2*nbig
         call zgeev(jobvl,jobvr,nfree,smat,nsmall,seig,svec,nsmall,
     &        svec,nsmall,swork,lwork,srwork,inform)
         if(inform.ne.0) write(6,1299) inform
         phasesum = 0.0
         do 573 i = 1,nfree
            arg1=imag(seig(i))
            arg2=dble(seig(i))
            phase = atan2(arg1,arg2)
            phase=phase*0.5
            rtest = abs(seig(i))
            phasesum=phasesum + phase
            write(6,574) i, phase,rtest
 573     continue
         write(6,576)phasesum,(kchan(i),i=1,nopen)
         if (ibcondx.ne.0) then
            eincidnt = kchan(1)**2/2.
         endif
 1000 continue
c     
 595  format(/" this is an optical potential calculation"/)
 596  format(/" this is an effective core potential calculation"/)
 597  format(/" this is not an effective core potential calculation"/)
 731  format(//,' S-matrix boundary conditions')
 732  format(//,' T-matrix boundary conditions')
 733  format(//,'  boundary cond. on file intsff is nonsense')
 138  format(' number of scattering orbitals in each channel:', 6i3)
 139  format(' channel',i3,' orbs:',20i3,/,(17x,20i3))
 137  format(' number of target orbitals in each channel:', 6i3)
 787  format(//,' *****orthogonalization is on*****',//,
     &     ' free functions will be orthogonalized to following orbitals')
 789  format(' in channel',i3,' :  ',10i3/(18x,10i3))
 791  format(//,1x)
 111  format(///,' ********************************',
     #     ' new energy ********************************',
     &     //,' channel momenta:',6e12.5)
 99   format(" transformation matrix")
 100  format(i5,6f12.7/(5x,6f12.7))
 112  format(' off-shell momenta:',6e12.4)
 107  format(//' solution matrix ')
 108  format(1x,i3,6e12.5,/,(4x,6e12.5))
 109  format(//' scattered wave part  ')
 776  format(//,' S-matrix for channel momenta:',3e12.6)
 177  format(4i5,2f20.10)
 277  format(4e20.10)
 874  format(//,' total cross sections: row index = initial chnl,',
     &     ' column index = final chnl')
 871  format(4x,6(6x,i2,4x))
 872  format(1x,i2,1x,6e12.5/(4x,6e12.5))
 778  format(//,' T-matrix for channel momenta:',3e12.6)
 777  format(1x,i3,6e12.5,/,(4x,6e12.5))
 793  format(//,' unitarity check: Hermitian conjg(S) * S')
 1299 format(//,' *** S-matrix diagonalization failed, inform =',i5)
 574  format(' eigenphase #',i3,' =',e12.5,' modulus =',e12.5)
 576  format(' eigenphase sum =',e12.5,
     &     ' for channel momenta:', (6e12.5))
 884  format(' incident energy =',f15.10)
c
      return
      end
