
      block data feffbd

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header

      character*10 shole(0:9)
      character*8 sout(0:6)
      common /labels/ shole, sout


      character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
      common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch

c     character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
c     common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch

      data shole /'no hole',    'K shell',
     1            'LI shell',   'LII shell',
     2            'LIII shell', 'MI shell',
     3            'MII shell',  'MIII shell',
     4            'MIV shell',  'MV shell'/
      data sout /'H-L exch', 'D-H exch', 'Gd state', 
     1           'DH - HL ', 'DH + HL ', 'HLnoimag', 'Gs HL   '/

c                   123456789012
      data vfeff  /'  Feff 6.01l'/
      data vpotph /'  potph 4.12'/
      data vpaths /'  paths 3.05'/
      data vgenfm /' genfmt 1.44'/
      data vff2ch /' ff2chi 2.01'/

c     6.01l EXAFS only lite version 10/02 jjr
c     5.05a is current working version
c     5.05j is jjr's version 6/93
c     6.00 Alexey's polarization and XANES 
c     6.01 Release version of FEFF6 including bug fixes ala and jjr
c     4.04 Major code reorganization.  Muffin tin finder modified -- now
c     uses average of all possible muffin tin radii instead of minimum.
c     26 March, 1991   Steven Zabinsky
c     4.05 Yet another improvement to muffin tin finder, now averages
c     based on volume of lense-shaped overlapping region, April, 1991
c     4.06 Bug fix in sumax, april 1991
c     4.07 Several minor changes involving non-standard F77 6/6/91, siz
c     4.08 ION card added 7/24/91, siz
c     4.08a, bug in header for ION card fixed 9/10/91, siz
c     4.09, quinn correction added to imhl, interstitial calculation
c           corrected, rmt modified to handle too few neighbors and
c           error msg in phase about hard test in fovrg modified,
c           folp card added
c     POTPH 4.1  Same as feff4.09, but version hacked to work with
c     module potph of feff5, Mar 1992, siz
c
c     new version common added, siz, Mar 1992
c     feff 5.03, first 'real' release, lots of little changes.
c                4 criteria added is the big change.  siz, April 1992
c     feffx 5.04, intermediate intermittent version of code with
c                 background, xsect, xmu, timereversal, lots
c                 of input cards, xanes, etc.  July 1992, siz
c     e REQUIRE card removed, Oct 92, siz
c     f, and paths 3.04, new crits, 9 points. Oct 92
c     g: major bug in xsect -  ixc not passed to xcpot, beginning with
c        5.04g, it's fixed.
c     h use gs for xsect (hard coded)
c     i fixed init and final state mixup in xsect
c     Feff 5.05, release version with all of the above in it.  XANES
c        is turned off in RDINP for the release -- turn it back on
c        there for development.
c     Feff 6 includes polarization (Alexey) and XANES (Steve Z.)
c     Feff 6.01 is the first release version of FEFF6.
c     Feff 6.01l EXAFS only lite version 10/02 jjr

      end
c     code: relativistic atom code (relativistic hartree fock slater)
c     modified desclaux code -- partially translated from the french
c
c     modified by: r. c. albers (from previously modified code from
c                  j. e. muller who in turn got it from the danes)
c                  j. j. rehr and s. i. zabinsky for inclusion in feff 
c
c     special features:  renormalizes charge density at wigner-seitz
c                        radius
c
c     version 2 (30 september 87): renormalized coulomb potential and
c     renormalized charge density are produced to be used in XAFS
c     calculations by cphase program. j.j. rehr, j. mustre  university
c     of washington., a.djaoui university of essex.
c     please acknowledge use. r. c. albers  (los alamos national lab)
c     j.j. rehr (university of washington),
c
c     Subroutine calling hierarchy siz 1/8/90
c     ATOM
c        INDATA
c           GETORB
c           FPOT
c        DIRAC
c           INOUH
c           INTH
c        POTSL
c        SOMM
c        TOTALE
c           SOMM
c        CDSLD
c           SOMM
c           YKDIR
c        RENORM
c           POTSLW
c
c     Version 1/11/90:  Input and output re-organized to work
c                       easily with overlapped potential code
c                       in FEFF.
c
c     Version Aug 1990: Minor modification to work more easily with
c                       FEFF4, cluster version.  SRHO no longer has
c                       factor of r**2.  INDATA uses rr function to
c                       set r grid.
c     Version Dec 1990: Writes to atom.dat restored
c     Version Feb 1991: Unit 16 opened in atom if necessary
c     June 1992  dirac upper and lower components and total energy
c                passed out for use with matrix element calculations
c
c     Input:   title    title, max 40 characters
c              ifr      index of free atom, used for output labels
c              iz       atomic number of atom
c              ihole    location of electron hole
c              rws      Wigner-Seitz radius
c              ionin    ionicity
c              iprint   print flag, passed through commom /print/
c              ispinr   0, do not save dirac spinors, else save for
c                       orbital ispinr
c
c     Output:  vcoul(251)  coulomb potential (no factor r**2)
c              srho(251)   electron density in form
c                          4*pi*density (formerly 4*pi*density*r**2)
c              dgc0(251)   large component (set if ispinr.ne.0)
c              dpc0(251)   small component (set if ispinr.ne.0)
c              eatom       total energy in rydbergs
c
c     All data is on a grid r(i) = exp (-8.8 + (i-1)*0.05)

      subroutine feff_atom(title,ifr,iz,ihole,rws,ionin,vcoul,srho,
     1                 ispinr, dgc0, dpc0, eatom)

      implicit double precision (a-h,o-z)
      save

c     Save central atom dirac components, see comments below.
      dimension dgc0(251), dpc0(251)

      character*(*)  title
      dimension vcoul(251)
      dimension srho(251)
      common /print/ iprint

      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco
      integer*4 nstop,nuc
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc

      common /ps2/ dexv, dexe, dcop, test, teste,
     1             testy, testv, niter, ion, icut, iprat, irnorm

      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1              dpc(251,30)

      character*40 ttl
      character*2  titre
      common /char2/ titre(30), ttl

      dimension tden(30)
      character*30 fname

      data harryd /2./
      character*72 header
      common /header_common/ header


      if (iprint .ge. 3)  then
c        prepare file for atom output
         write(fname,14)  ifr
   14    format('atom', i2.2, '.dat')
         open (unit=16, file=trim(header)//fname, 
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//fname, 'atom')
c        call head (16)
         write(16,*)  ' Free atom ', ifr
      endif

      ttl = title

      nstop=1
      mark=0

      call indata (iz, ihole, rws, ionin)
      iter=1
      do 30 i=1,np
      do 30 j=1,norb
         dgc(i,j)=0.0
         dpc(i,j)=0.0
   30 continue

      if (iprint .ge. 3)  write(16,40) ttl
   40 format (1h1,40x,a40)
      n=-(ion+1)

   60 continue
      do 70 i=1,np
         d(i)=0.0
   70 continue
      tets=test
      ymax=0.0
      vmax=0.0
      emax=0.0

c resolution of the dirac equation for each orbital
      do 150 j=1,norb
         de=den(j)
   80    call feff_dirac (nqn(j),nql(j),nk(j),imax,den(j),
     D        dfl(j),dq1(j),j)
            if (nstop.eq.0) go to 110
            if (nstop.ne.362.or.iter.ge.10.or.tets.gt.test) go to 90
            tets=testv
         go to 80
   90    if (iprint .ge. 3)  write(16,100) nstop,nqn(j),titre(j)
  100    format ('  nstop=',i4,'  for the orbital',i3,a2)
         write(77,*) ' Fatal error.'
         write(77,*) ' Wigner-Seitz or muffin tin radius may be',
     1              ' too small.'
         go to 999

  110    val=abs((den(j)-de)/de)
         if (val.gt.emax) emax=val
         nmax(j)=imax
         do 140 i=1,np
            val=dgc(i,j)-dp(i)
            if (abs(dp(i)).gt.1.0d0) val=val/dp(i)
            if (abs(val).lt.abs(ymax)) go to 120
               ymax=val
               y=dp(i)
               yn=dgc(i,j)
  120       val=dpc(i,j)-dq(i)
            if (abs(dq(i)).gt.1.0d0) val=val/dq(i)
            if (abs(val).lt.abs(ymax)) go to 130
               ymax=val
               y=dq(i)
               yn=dpc(i,j)
  130       dgc(i,j)=dp(i)
            dpc(i,j)=dq(i)
  140    d(i)=d(i)+nel(j)*(dp(i)*dp(i)+dq(i)*dq(i))
  150 continue

c     dgc and dpc are set in loop above, only referenced in remainder
c     of code, so save them into dgc0 and dpc0 here.  Note: np=251,
c     set in indata.  dgc0 is large component
c                     dpc0 is small
      if (ispinr .ne. 0)  then
         do 152  i = 1, np
            dgc0(i) = dgc(i,ispinr)
            dpc0(i) = dpc(i,ispinr)
  152    continue
      endif

      if (mark.eq.0) go to 280

c  This is case mark .ne. 0
c  d is the core electron density resulting from the renormalized pot.
      dval=0.0
      do 160 j=1,norb
  160    dval=dval+nel(j)*den(j)

      dval=dval*2.0
c jm-- core charge density commented away in unit 6 appears in unit 3--
      if (iprint .ge. 3)  write(16,170) dval
  170 format (1h ,' core energy = ',e15.8)

c jm- renormalized potential

c     note conversion to rydbergs using constant harryd
c     passvt is part of old system to pass data directly from
c     ATOM to PHASE
c      do 200 ixx=1,251
c  200    passvt(ixx)=harryd*dr(ixx)*dr(ixx)*dv(ixx)


c  d is the core electron density resulting from the renormalized pot.

c  next write renormalized electron density for each shell
      do 270 j=1,norb
         do 240 i=1,np
            d(i)=dgc(i,j)*sqrt(12.56637062d0)
  240    continue
  270 continue
      go to 750

c     mark .eq. 0 case
  280 continue

      call potsl (dc,d,dp,dr,dpas,dexv,z,np,ion,icut,dvn)
      if (nuc.le.0) go to 300
         do 290 i=1,nuc
            dc(i)=dc(i)+z/dr(i)+z*((dr(i)/dr(nuc))**2-3.0d0) /
     1            (dr(nuc)+dr(nuc))
  290    continue
  300 continue
      do 310 i=1,np
         dval=abs(dc(i)-dv(i))
         if ((dr(i)*dc(i)).le.n) dval=-dval/dc(i)
         if (dval.le.vmax) go to 310
            vmax=dval
            j=i
  310 continue

c     print 320, iter,vmax,dr(j),dv(j),dc(j),emax,ymax,yn,y
c 320 format (i5,1pe11.2,3(1pe16.6),2(1pe11.2),2(1pe16.6))

      if (tets.le.test.and.emax.le.teste.and.vmax.le.testv.and.ymax.le
     1 .testy) go to 430
      if (mark.eq.1) go to 430
      iter=iter+1
      if (iter.le.niter) go to 340
      if (iprint .ge. 3)  write(16,330) niter
  330 format (' number of iterations greater than',i4)
      nstop=2
c      print*, ' ATOM-Fatal error, too many iterations.'
c      print*, '   iter, niter ', iter, niter
      write(77,*) ' ATOM-Fatal error, too many iterations.'
      write(77,*) '   iter, niter ', iter, niter
      go to 999
c potential for the following iteration

  340 continue
      if (iter.eq.2) go to 350
      if (iprat.eq.0) then
      do 400 i=1,np
      dval=dalp(dvn(i),dvf(i),dv(i),dc(i))
      dvn(i)=dv(i)
      dvf(i)=dc(i)
  400 dv(i)=dval*dv(i)+(1.0d0-dval)*dc(i)
      else
  350 dval=1.0-dcop
      do  i=1,np
      dvn(i)=dv(i)
      dvf(i)=dc(i)
      dv(i)=dval*dv(i)+dcop*dc(i)
      enddo
      endif
      go to 60

  430 if (iprint .ge. 3)  write(16,40) ttl
      if (iprint .ge. 3)  write(16,460)
  460 format (12x,'energie',12x,'(r4)',14x,'(r2)',14x,'(r)',15x,'(r-1)',
     1 13x,'(r-3)'/)

c valeurs moyennes de r
      do 470 i=1,np
      dvf(i)=dc(i)
  470 dq(i)=0.0
      dval=0.0
      do 560 i=1,norb
      im=nmax(i)
      dval=dval+nel(i)*den(i)
      do 480 j=1,im
  480 dc(j)=dgc(j,i)*dgc(j,i)+dpc(j,i)*dpc(j,i)
      l=5
      if (iabs(nk(i)).eq.1) l=l-1
      do 550 j=1,l
      dp(j)=dfl(i)+dfl(i)
c      if (j-2) 490,500,510
      if(j.lt.2) then
  490 n=4
      elseif(j.eq.2) then
  500 n=2
      else
c  510 if (j-4) 520,530,540
      if(j.lt.4) then
  520 n=1
      elseif(j.eq.4) then
  530 n=-1
      else
  540 n=-3
      endif
      endif
  550 call somm (dr,dc,dq,dpas,dp(j),n,im)
  560 if (iprint .ge. 3)  write(16,570) nqn(i),titre(i),
     1                                   den(i),(dp(j),j=1,l)
  570 format (i3,a2,6(1pe18.7))

      if (dexv.eq.0.0) go to 650

c energie totale en moyenne spherique
      do 580 i=1,norb
  580 tden(i)=-2.0d0*den(i)

      dc(1)=1
      do 600 i=1,np
  600 dp(i)=d(i)/dr(i)
      if (nuc.le.0) go to 620
      do 610 i=1,nuc
  610 dp(i)=d(i)*(3.0d0-dr(i)*dr(i)/(dr(nuc)*dr(nuc)))/(dr(nuc)+dr(nuc))
      dc(1)=4
  620 call somm (dr,dp,dq,dpas,dc(1),0,np)
      do 630 i=1,np
      dp(i)=d(i)*dvf(i)
  630 d(i)=d(i)*((d(i)*dr(i))**(1.0d0/3.0d0))
      dc(2)=3
      dc(3)=1
      if (nuc.ne.0) dc(3)=4
      call somm (dr,dp,dq,dpas,dc(3),0,np)
      call somm (dr,d,dq,dpas,dc(2),-1,np)
      dc(2)=-3.0d0*dc(2)/(105.27578d0**(1.0d0/3.0d0))
      dc(1)=-z*dc(1)
      dc(4)=dval-dc(3)
      dval=dval+(dc(1)-dc(3)+(dexe-dexv)*dc(2))/2.0d0
      dc(3)=(dc(3)-dc(1)-dexv*dc(2))/2.0d0
      dc(2)=dc(2)*dexe/2.0d0
      if (iprint .ge. 3)  write(16,640) dval,dc(4),dc(3),dc(2),dc(1)
  640 format (1h0,5x,'et=',1pe14.7,5x,'ec=',1pe14.7,5x,'ee=',1pe14.7,5x,
     1 'ex=',1pe14.7,5x,'en=',1pe14.7)
      go to 660
  650 call totale (dval)
  660 continue

c     pass out eatom (total energy) (factor of 2 is to put energy in
c     rydberg units)
      eatom = 2 * dval

      if (norb.eq.1) go to 710
      if (iprint .ge. 3)  write(16,40) ttl
      if (iprint .ge. 3)  write(16,670)
  670 format (1h0,47x,'overlap integrals         '/)

c overlap integrals
      do 700 i=2,norb
      k=i-1
      do 700 j=1,k
      if (nql(i).ne.nql(j).or.nk(i).ne.nk(j)) go to 700
      im=nmax(j)
      if (nmax(i).lt.im) im=nmax(i)
      do 680 l=1,im
      dq(l)=dpc(l,i)*dpc(l,j)
  680 dc(l)=dgc(l,i)*dgc(l,j)
      dval=dfl(i)+dfl(j)
      call somm (dr,dc,dq,dpas,dval,0,im)
      if (iprint .ge. 3)  write(16,690) nqn(i),titre(i),
     1                                   nqn(j),titre(j),dval
  690 format (34x,i1,a2,i3,a2,f19.7)
  700 continue
  710 call cdsld


      if (irnorm.eq.1) then
         call renorm (dexv, vcoul, srho)
      endif
      do 720 i=1,np
  720 dc(i)=harryd*dv(i)*dr(i)**2
      if (irnorm.ne.1) stop 0000
      norb=norbco
      if (norbco.eq.0) go to 750
      if (mark.eq.1) go to 750
      mark=1
      go to 60

  750 continue

c     return srho as 4*pi*density instead of 4*pi*density*r**2
      do 760  i = 1, 251
         srho(i) = srho(i) / (dr(i)**2)
  760 continue

      if (iprint .ge. 3)  close(unit=16)

      return


  999 continue
      stop 'ATOM-1'
      end
      subroutine besjn (x, jl, nl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 30 (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this array nl = abramowitz yl.
c       jl and nl must be dimensioned 
c            complex*16 jl(ltot+2), nl(ltot+2), with ltot defined in 
c            dim.h.
c
c     notes:  jl and nl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c
c-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)

      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      complex*16 x
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 cjl(ltot+2), sjl(ltot+2), cnl(ltot+2), snl(ltot+2)

      complex*16 xjl,xnl,asx,acx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      parameter (xcut = 1.0d0, xcut1 = 7.51d0, xcut2 = 5.01d0)

      if (dble(x) .le. 0)  stop 'Re(x) is .le. zero in besjn'

      lmaxp1 = ltot+2

      if (dble(x) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 il = 1,lmaxp1
            l = il-1
            ifl = 0
            call bjnser (x,l,xjl,xnl,ifl)
            jl(il) = xjl
            nl(il) = xnl
   10    continue

      elseif (dble(x) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmaxp1-1,xjl,xnl,1)
         jl(lmaxp1) = xjl

         call bjnser (x,lmaxp1-2,xjl,xnl,1)
         jl(lmaxp1-1) = xjl

         if (dble(x) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1.0d0 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lx = 3,lmaxp1
            lp1 = lmaxp1+1-lx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(lp1) = tlxp3 * jl(lp1+1) / x  -  jl(lp1+2)
   60    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1.0d0 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3.0d0*xi3 - xi
         sjl(4) = 15.0d0*xi4 - 6.0d0*xi2
         sjl(5) = 105.0d0*xi5 - 45.0d0*xi3 + xi
         sjl(6) = 945.0d0*xi6 - 420.0d0*xi4 + 15.0d0*xi2
         sjl(7) = 10395.0d0*xi7 - 4725.0d0*xi5 + 210.0d0*xi3 - xi
         sjl(8) = 135135.0d0*xi8 - 62370.0d0*xi6 + 3150.0d0*xi4 
     >             - 28.0d0*xi2
         sjl(9) = 2027025.0d0*xi9 - 945945.0d0*xi7 + 51975.0d0*xi5 
     1            - 630.0d0*xi3 + xi
         sjl(10) = 34459425.0d0*xi10 - 16216200.0d0*xi8 +945945.0d0*xi6
     1            - 13860.0d0*xi4 + 45.0d0*xi2
         sjl(11) = 654729075.0d0*xi11 - 310134825.0d0*xi9 
     >            + 18918900.0d0*xi7 
     1            - 315315.0d0*xi5 + 1485.0d0*xi3 - xi
         cjl(1) = 0.0d0
         cjl(2) = -xi
         cjl(3) = -3.0d0*xi2
         cjl(4) = -15.0d0*xi3 + xi
         cjl(5) = -105.0d0*xi4 + 10.0d0*xi2
         cjl(6) = -945.0d0*xi5 + 105.0d0*xi3 - xi
         cjl(7) = -10395.0d0*xi6 + 1260.0d0*xi4 - 21.0d0*xi2
         cjl(8) = -135135.0d0*xi7 + 17325.0d0*xi5 - 378.0d0*xi3 + xi
         cjl(9) = -2027025.0d0*xi8 + 270270.0d0*xi6 - 6930.0d0*xi4 
     >            + 36.0d0*xi2
         cjl(10) = -34459425.0d0*xi9 + 4729725.0d0*xi7 - 135135.0d0*xi5 
     1             + 990.0d0*xi3 - xi
         cjl(11) = -654729075.0d0*xi10 + 91891800.0d0*xi8 
     >             - 2837835.0d0*xi6 
     1             + 25740.0d0*xi4 - 55.0d0*xi2
         do 80 ie = 1,11
            snl(ie) = cjl(ie)
            cnl(ie) = -sjl(ie)
   80    continue
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = dble(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
            snl(lp1) = tlxp1*xi*snl(lp1-1)-snl(lp1-2)
            cnl(lp1) = tlxp1*xi*cnl(lp1-1)-cnl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         do 110 lp1 = 1,lmaxp1
            jl(lp1) = asx*sjl(lp1)+acx*cjl(lp1)
            nl(lp1) = asx*snl(lp1)+acx*cnl(lp1)
  110    continue
      endif

      return
      end
      subroutine bjnser (x, l, jl, nl, ifl)

c-----------------------------------------------------------------------
c
c     subroutine: bjnser (x,l,jl,nl,ifl)
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c
c     arguments:
c       x = argument of jl and nl
c       l = l value calculated (no offset)
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c       ifl = 0 return both jl and nl
c             1 return jl only
c             2 return nl only
c
c     notes:  jl and nl are calculated by a series
c             expansion according to 10.1.2 and 10.1.3
c             in abramowitz and stegun (ninth printing),
c             page 437
c
c             error msgs written with PRINT statements.
c
c     first coded by r. c. albers on 26 jan 83
c
c     version 2
c
c     last modified: 27 jan 83 by r. c. albers
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      complex*16 x,u,ux,del,pj,pn
      complex*16 jl,nl

      parameter (niter = 20, tol = 1.d-15)

      if (l .lt. 0) then
         write(77,*) 'l .lt. 0 in bjnser'
         stop 'bjnser 1'
      endif
   20 if (dble(x).lt. 0.0d0) then
         write(77,30) x
   30    format (/, ' x = ', 1p, 2e14.6, ' is .le. 0 in bjnser')
         stop 'bjnser 2'
      endif

      lp1 = l+1
      u = x**2 / 2.0d0

c     make djl = 1 * 3 * 5 * ... * (2*l+1),
c          dnl = 1 * 3 * 5 * ... * (2*l-1)
      djl = 1
      fac = -1
      do 50 il = 1, lp1
         fac = fac + 2
         djl = fac * djl
   50 continue
      dnl = djl / (2*l+1)


      if (ifl .eq. 2)   goto 90
c     make jl
c     pj is term in { } in 10.1.2, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pj = 1
      nf = 1
      nfac = 2*l + 3
      den = nfac
      sgn = -1
      ux = u
      do 60 il = 1, niter
         del = sgn*ux / den
         pj = pj + del
         trel = abs (del / pj)
         if (trel .le. tol)  goto 80
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
   60 continue
      stop  'jl does not converge in bjnser'
   80 jl = pj * (x**l) / djl

   90 if (ifl.eq.1) return
c     make nl
c     pn is term in { } in 10.1.3, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pn = 1
      nf = 1
      nfac = 1 - 2*l
      den = nfac
      sgn = -1
      ux = u
      do 100  il = 1, niter
         del = sgn * ux / den
         pn = pn + del
         trel = abs (del / pn)
         if (trel .le. tol) goto 120
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
  100 continue
      stop  'nl does not converge in bjnser'
  120 nl = -pn * dnl / (x**lp1)

      return
      end
      subroutine ccrit(npat, ipat, ckspc,
     1    fbetac, rmax, pcrith, pcritk, nncrit, ipotnn, ipot,
     2    rpath, lheap, lkeep, xcalcx)
      implicit double precision (a-h, o-z)

c     lheap to add to heap, lkeep if keep path at output.
c     NB, if lheap is false, lkeep is not used (since path
c     won't be in the heap).


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      logical lheap, lkeep
      dimension ipat(npatx)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)

c     local variables
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)

c     mrb is efficient way to get only ri and beta
c     note that beta is cos(beta)
      call mrb (npat, ipat, ri, beta)

      rpath = 0.0d0
      do 300  i = 1, npat+1
         rpath = rpath + ri(i)
  300 continue

c     If we can decide only on rpath, do it here...
      if (rpath .gt. rmax)  then
         lheap = .false.
         lkeep = .false.
         return
      endif

c     If last atom central atom, do put in heap, don't use it
c     as an actual path at output
      if (ipat(npat).eq.0)  then
         lheap = .true.
         lkeep = .false.
         return
      endif

c     Make index into fbetac array (this is nearest cos(beta) grid 
c     point, code is a bit cute [sorry!], see prcrit for grid).
      do 290  i = 1, npat+1
         tmp = abs(beta(i))
         n = tmp / 0.025d0
         del = tmp - n*0.025d0
         if (del .gt. 0.0125d0)  n = n+1
         if (beta(i) .lt. 0.0d0)  n = -n
         indbet(i) = n
  290 continue

c     Decide if we want the path added to the heap if necessary.
c     (Not necessary if no pcrith in use.)
      if (pcrith .gt. 0)  then

         call mcrith(npat, ipat, ri, indbet,
     1                ipot, nncrit, fbetac, ckspc, xheap)

c        xheap = -1 if not defined for this path (too few legs, etc.)
         if (xheap .ge. 0  .and.  xheap .lt. pcrith)  then
c           Do not want path in heap
            lheap = .false.
            lkeep = .false.
            return
         endif
      endif
c     Keep this path in the heap
      lheap = .true.

c     We may want path in heap so that other paths built from this
c     path will be considered, but do not want this path to be
c     written out for itself.  Decide that now and save the flag
c     in the heap, so we won't have to re-calculate the mpprm
c     path parameters later.

c     Skip calc if pcritk < 0
      if (pcritk .le. 0)  then
         lkeep = .true.
         return
      endif

c     Make xout, output inportance factor.
      call mcritk (npat, ipat, ri, beta, indbet,
     1             ipot, nncrit, fbetac, ckspc, xout, xcalcx)

c     See if path wanted for output
c     Do not want it if last atom is central atom (xout = -1) or
c     if xout is too small
      lkeep = .false.
      if (xout .ge. pcritk)  lkeep = .true.

      return
      end
      subroutine cdsld

      implicit double precision (a-h,o-z)
      save
      common /print/ iprint
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco

      integer*4 nstop,nuc
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)

c titre = identification of the wave functions  s,p*,p,........
      character*40 ttl
      character*2  titre
      integer*4 n1,n1sum
      common /char2/ titre(30), ttl

c  -- This read commented out to make input easier, not used for
c     PHASE calculations
      irm  = 0
      ins  = 0
      npun = 0
      nmfg = 0
      nmrk = 0
c     read (5,10) irm,ins,npun,nmfg,nmrk
   10 format (8i3)

c valeurs moyennes de r**j  if irm non-zero
c tabulation of the wave functions if ins non-zero
c the potential multiplied by r is perfore if npun non-zero
      if (irm.eq.0) go to 200
      if (iprint .ge. 5)  write(16,20) ttl
   20 format (1h1,40x,a40,/)
   30 read (5,10) j,l,n1,l1,j1,n2,l2,j2
      if (l.eq.0) go to 200

c valeur moyenne of (p1*p2+q1*q2)*r**j  if l positive
c valeur moyenne of (p1*q2+p2*q1)*r**j  if l negative
      if (n1.gt.0) go to 40
      if (((n1+1)*(n1+2)).ne.0) go to 60
      i1=1
      i2=1
      go to 80
   40 i1=0
      i2=0
      do 50 i=1,norb
      if (nqn(i).eq.n1.and.nql(i).eq.l1.and.(j1-1).eq.(-nk(i)/iabs(nk(i)
     1 ))) i1=i
      if (nqn(i).eq.n2.and.nql(i).eq.l2.and.(j2-1).eq.(-nk(i)/iabs(nk(i)
     1 ))) i2=i
   50 continue
      if (i1.ne.0.and.i2.ne.0) go to 80
   60 if (iprint .ge. 5)  write(16,70) j,l,n1,l1,j1,n2,l2,j2
   70 format (1h0,'    error for the card     ',8i3)
      go to 30
   80 dval=dfl(i1)+dfl(i2)
      if ((dval+j).gt.-1.0) go to 90
      if (n1) 170,170,60
   90 im=nmax(i1)
      if (nmax(i2).lt.im) im=nmax(i2)
      if (l.lt.0) go to 110
      do 100 i=1,im
      dv(i)=dgc(i,i1)*dgc(i,i2)
  100 dq(i)=dpc(i,i1)*dpc(i,i2)
      go to 130
  110 do 120 i=1,im
      dv(i)=dgc(i,i1)*dpc(i,i2)
  120 dq(i)=dgc(i,i2)*dpc(i,i1)
  130 call somm (dr,dv,dq,dpas,dval,j,im)
      if (l.lt.0) go to 150
      if (iprint .ge. 5)  write(16,140) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dval
  140 format (24x,'(p1p2+q1q2)r**',i2,' for  ',i1,a2,i3,a2,5x,'=',1pe14.
     1 7,/)
      go to 170
  150 if (iprint .ge. 5)  write(16,160) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dval
  160 format (24x,'(p1q2+q1p2)r**',i2,' for  ',i1,a2,i3,a2,5x,'=',1pe14.
     1 7,/)
  170 n1sum=n1+1
      if (n1sum) 190,180,30
  180 i1=i1+1
      i2=i1
      n1sum=i1-norb
      if (n1sum) 80,80,30
  190 i2=i2+1
      n1sum=i2-norb
      if (n1sum) 80,80,180
  200 if (ins.eq.0) go to 260
      do 250 i=1,norb,3
      j=i+2
      if (j.gt.norb) j=norb
      im=0
      do 210 l=i,j
      if (nmax(l).gt.im) im=nmax(l)
  210 continue
      do 230 k=1,im
      if (((k-1)*(k-48*(k/48))).ne.0) go to 230
      if (iprint .ge. 5)  write(16,20) ttl
      if (iprint .ge. 5)  write(16,220) (nqn(l),titre(l),nqn(l),
     1                                     titre(l),l=i,j)
  220 format (9x,'r',14x,3(i1,a2,'g.c.',i11,a2,'p.c.',10x))
  230 if (iprint .ge. 5)  write(16,240) dr(k),
     1                                   (dgc(k,l),dpc(k,l),l=i,j)
  240 format (7(1pe17.7))
  250 continue
  260 if (npun.eq.0) go to 300
      do 270 i=1,np
  270 dp(i)=dvf(i)*dr(i)
c     write(8,280) ttl
  280 format (a40)
c     write(8,290) (dp(i),i=1,np)
  290 format (8f9.4)
  300 do 310 i=1,np
  310 d(i)=0.0
      nag=1
      if (nmfg.eq.0) go to 470
      if (iprint .ge. 5)  write(16,20)
      if (iprint .ge. 5)  write(16,320)
  320 format (/,30x,'integrales magnetiques directes et d echange'//)
  330 read (5,10) i1,i2,n1
      if (i1.le.0) go to 470
      if (i2.gt.0) go to 350
      if (((i2+1)*(i2+2)).ne.0) go to 340
      if (n1.le.0) n1=1
      i1=n1
      n1=i2
      i2=i1
      go to 360
  340 if (iprint .ge. 5)  write(16,70) i1,i2,n1
      go to 330
  350 if (i1.gt.norb.or.i2.gt.norb) go to 340
      n1=1
  360 j1=2*iabs(nk(i1))-1
      j2=2*iabs(nk(i2))-1
      kma=min0(j1,j2)
      nm=nmax(i2)
      do 380 j=1,kma,2
      call ykdir (i1,i1,j,nag)
      do 370 i=1,nm
  370 dp(i)=dq(i)*dgc(i,i2)*dpc(i,i2)
      dval=j+1
      call somm (dr,d,dp,dpas,dval,-1,nm)
  380 if (iprint .ge. 5)  write(16,390) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dval
  390 format (20x,'fm',i2,' (',i1,a2,',',i1,a2,') =',1pe14.7)
      if (i1.eq.i2) go to 440
      j1=(iabs(1-2*nk(i1))-1)/2
      j2=(iabs(1-2*nk(i2))-1)/2
      kma=max0(nql(i1)+j2,nql(i2)+j1)
      j1=iabs(nql(i2)-j1)
      j2=iabs(nql(i1)-j2)
      kmi=min0(j1,j2)
      j1=kmi+nql(i1)+nql(i2)
      j1=j1-2*(j1/2)
      if (j1.eq.0) kmi=kmi+1
      nm=min0(nmax(i1),nmax(i2))
      do 420 j=kmi,kma,2
      call ykdir (i1,i2,j,nag)
      do 400 i=1,nm
      dp(i)=dq(i)*dgc(i,i1)*dpc(i,i2)
  400 dc(i)=dq(i)*dgc(i,i2)*dpc(i,i1)
      dval=j+1
      dvalp=dval
      dvalm=dval
      call somm (dr,d,dp,dpas,dvalp,-1,nm)
      call somm (dr,d,dc,dpas,dval,-1,nm)
      call ykdir (i2,i1,j,nag)
      do 410 i=1,nm
  410 dp(i)=dq(i)*dgc(i,i2)*dpc(i,i1)
      call somm (dr,d,dp,dpas,dvalm,-1,nm)
  420 if (iprint .ge. 5)  write(16,430) j,nqn(i1),titre(i1),nqn(i2),
     1                                   titre(i2),dvalm,dval,dvalp
  430 format (' gm',i2,' (',i1,a2,',',i1,a2,')',5x,'(-1)=',1pe14.7,5x,'(
     10)=',1pe14.7,5x,'(+1)=',1pe14.7)
  440 n1sum=n1+1
      if (n1sum) 460,450,330
  450 i1=i1+1
      i2=i1
      n1sum=i1-norb
      if (n1sum) 360,360,330
  460 i2=i2+1
      n1sum=i2-norb
      if (n1sum) 360,360,450
  470 if (nmrk.eq.0) go to 530
      if (iprint .ge. 5)  write(16,20)
      if (iprint .ge. 5)  write(16,480)
  480 format (/,20x,'integrales magnetiques rk=integrale de p1(1)*q2(1)*
     1uk(1,2)*p3(2)*q4(2)'//)
  490 read (5,10) i1,i2,i3,i4,k
      if (i1.le.0) go to 530
      if (i1.le.norb.and.i2.gt.0.and.i2.le.norb.and.i3.gt.0.and.i3.le
     1 .norb.and.i4.gt.0.and.i4.le.norb.and.k.ge.0) go to 500
      if (iprint .ge. 5)  write(16,70) i1,i2,i3,i4,k
      go to 490
  500 call ykdir (i1,i2,k,nag)
      do 510 i=1,np
  510 dp(i)=dq(i)*dgc(i,i3)*dpc(i,i4)
      dval=k+1
      call somm (dr,d,dp,dpas,dval,-1,np)
      if (iprint .ge. 5)  write(16,520) k,nqn(i1),titre(i1),nqn(i2),
     1              titre(i2),nqn(i3),titre(i3),nqn(i4),titre(i4),dval
  520 format (20x,'rm',i2,' (',i1,a2,',',i1,a2,',',i1,a2,',',i1,a2,') ='
     1 ,1pe14.7)
      go to 490
  530 return
      end
      subroutine chopen (ios, fname, mod)
      implicit double precision (a-h, o-z)
c     Writes error msg and stops if error in ios flag from open
c     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod

c     open successful
      if (ios .le. 0)  return

c     error opening file, tell user and die.
      write(77,100) fname, mod

  100 format (' ERROR opening file, ', /,
     1        ' filename:  ', a, /,
     2        ' in module: ', a)

      write(77,*) 'Fatal error'
      stop 'CHOPEN'
      end
      subroutine cpl0 (x, pl0, lmaxp1)
      implicit double precision (a-h, o-z)

c-----------------------------------------------------------------------
c
c     cpl0:  Calculate associated legendre polynomials p_l0(x)
c            by recursion.
c            Adapted from aslgndr.
c
c     first written: (25 june 86) by j. j. rehr
c
c     version 1 (25 june 86) (aslgndr)
c     version 2 (March, 1992) siz
c
c-----------------------------------------------------------------------

      dimension pl0 (lmaxp1)

      lmax = lmaxp1-1

c     calculate legendre polynomials p_l0(x) up to l=lmax
      pl0(1) = 1.0d0
      pl0(2) = x
      do 10  il = 2, lmax
         l = il-1
         pl0(il+1) = ( (2*l+1)*x*pl0(il) - l*pl0(l) ) / il
   10 continue

      return
      end
c Copyright Notice: FEFF6 is copyright protected software and users must
c obtain a license from the University of Washington Office of
c Technology Transfer for its use; see section V of FEFF document.

c Main Authors of FEFF5: please contact us concerning problems.
c A. L. Ankudinov, alex@phys.washington.edu      (206) 543 0435
c S. I. Zabinsky, zabinsky@phys.washington.edu   (206) 543 0435
c J. J. Rehr,     jjr@phys.washington.edu        (206) 543 8593
c R. C. Albers,   rca@nidhug.lanl.gov            (505) 665 0417

c Citations: Please cite at least one of the following articles if 
c FEFF is used in published work: 
c    1) Multiple scattering
c       J.J. Rehr and R.C. Albers, Phys. Rev. B41, 8139 (1990).
c       J.J. Rehr, S.I. Zabinsky and R.C. Albers, 
c          Phys. Rev. Let. 69, 3397 (1992).
c    2) General reference
c       J.J. Rehr, J. Mustre de Leon, S.I. Zabinsky, and R.C. Albers,
c          J. Am. Chem. Soc. 113, 5135 (1991).
c    3) Technical reference
c       J. Mustre de Leon, J.J. Rehr, S.I.  Zabinsky, and R.C. Albers,
c          Phys. Rev. B44, 4146 (1991).


      subroutine csomm (dr,dp,dq,dpas,da,m,np)
c Modified to use complex p and q.  SIZ 4/91
c integration by the method of simpson of (dp+dq)*dr**m from 
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),dq(*),da,dc
      mm=m+1
      d1=da+mm
      da=0.0d0
      db=0.0d0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      da=da+dc
      dc=dq(i)*dl
      da=da+dc
   70 continue
      da=dpas*da/3.0d0
      dd=exp(dpas)-1.0d0
      db=d1*(d1+1.0d0)*dd*exp((d1-1.0d0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0d0+1.0d0/(dd*(d1+1.0d0)))/d1
      da=da+dd*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine cubic (xk0, wp, alph, rad, qplus, qminus)

c     input:  xk0, wp, alph
c     output: rad, qplus, qminus

      implicit double precision (a-h, o-z)
      complex*16 s1,s13
      parameter (three = 3.0d0)
      parameter (third = 1.0d0/three)

c     this subroutine finds the roots of the equation
c     4xk0 * q^3  +  (alph-4xk0^2) * q^2  +  wp^2 = 0
c     see abramowitz and stegun pg 17 for formulae.

      a2 = (alph / (4.0d0*xk0**2)  -  1.0d0) * xk0
      a0 = wp**2 / (4.0d0*xk0)
      a1 = 0.0d0
      q = a1/3.0d0 - a2**2/9.0d0
      r = (a1*a2 - 3.0d0*a0)/6.0d0  -  a2**3/27.0d0
      rad = q**3 + r**2
      if (rad .gt. 0.0d0) then
         qplus = 0.0d0
         qminus = 0.0d0
         return
      endif

      s13 = dcmplx (r, sqrt(-rad))
      s1 = s13 ** third
      qz1 = 2.0d0*s1 - a2/3.0d0
      qz2 = -(s1 + sqrt(three)*dimag(s1) + a2/3.0d0)
      qz3 = -(s1 - sqrt(three)*dimag(s1) + a2/3.0d0)
      qplus = qz1
      qminus = qz3

      return
      end
      double precision function dalp (d1,d2,d3,d4)
      implicit double precision (a-h,o-z)
      save
c
c procedure of pratt to accelerate the convergence
c d1=initial (n-1);   d2=final (n-1);   d3=initial (n);   d4=final (n);
c **********************************************************************
      if ((d1+d4).eq.(d2+d3)) go to 10
      d=(d4-d2)/((d1+d4)-(d2+d3))
      if (d.lt.0.0d0) go to 20
      if (d.lt.0.5d0) go to 30
   10 d=0.5d0
      go to 30
   20 d=0.0d0
   30 dalp=d
      return
      end
      subroutine feff_diff(v, dx, n, vm)
      implicit double precision (a-h,o-z)
      complex*16 v(n), vm(n)
      vm(1)=((6.0d0*v(2)+6.66666666667d0*v(4)+1.2d0*v(6))-(2.45d0*v(1)
     > +7.0d0
     1 5*v(3)+3.75d0*v(5)+.166666666667d0*v(7)))/dx
      vm(2)=((6.0d0*v(3)+6.66666666667d0*v(5)+1.2d0*v(7))-(2.45d0*v(2)
     > +7.0d0
     1 5*v(4)+3.75d0*v(6)+.166666666667d0*v(8)))/dx
      nm2=n-2
      do 10 i=3,nm2
   10 vm(i)=((v(i-2)+8.0d0*v(i+1))-(8.0d0*v(i-1)+v(i+2)))/12.0d0/dx
      vm(n-1)=(v(n)-v(n-2))/(2.0d0*dx)
      vm(n)=(v(n-2)*0.5d0-2.0d0*v(n-1)+1.5d0*v(n))/dx
      return
      end
      subroutine feff_dirac (nqn,nql,nk,imax,de,dfl,dq1,jc)
c
c solution of the dirac equation
c nqn=principal quantum number; nql=orbital quantum number
c nk=kappa quantum number;  imax=the last tabulated point of the
c wave function; de=energy;   dfl=power of the first term of the
c developpement limite; dq1=slope at the origin of dp or dq
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      integer*4 nstop,nuc
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, test,
     1              z, nstop, nes, np, nuc
c
c dv=potential in a.u. and negative;  dr=radial mesh
c dp=large component;    dq=small component;    dpas=exponential step;
c nes=number of attempts to adjust the energy
c z=atomic number; nstop controls the numeric integration
c test=precision obtained in the energies; np=maximum number of points
c finite nuclear size if nuc is non-zero
c **********************************************************************
      common /ps1/ dep(5), deq(5), db, dvc, dsal, dk, dm
c
c dep,deq=derivatives of op and dq;  db=energie/dvc;
c dvc=speed of light in a.u.; dsal=2.*dvc;  dk=kappa quantum number
c dm=exponential step/720., dkoef=1./720.
c **********************************************************************
      common /trois/ dpno(4,30), dqno(4,30)
      data dkoef /0.1388888888888888d-2/
      nstop=0
      dvc=137.0373d0
      dsal=dvc+dvc
      imm=0
      ies=0
      dk=nk
      lll=(nql*(nql+1))/2
      nd=0
      noeud=nqn-nql
      if (lll.ne.0) go to 10
      elim=-z*z/(1.5d0*nqn*nqn)
      go to 40
   10 elim=dv(1)+lll/(dr(1)*dr(1))
      do 20 i=2,np
      val=dv(i)+lll/(dr(i)*dr(i))
      if (val.le.elim) elim=val
   20 continue
      if (elim) 40,30,30
   30 nstop=17
c 2*v+l*(l+1)/r**2 is everywhere positive
c **********************************************************************
      return
   40 if (de.le.elim) de=elim*0.5d0
   50 if (imm.eq.1) go to 80
      do 60 i=7,np,2
      imat=np+1-i
      if ((dv(imat)+lll/(dr(imat)*dr(imat))-de).le.0.0d0) go to 70
   60 continue
   70 if (imat.gt.5) go to 80
      de=de*0.5d0
      if (de.lt.-test.and.nd.le.noeud) go to 50
      nstop=28
c 2*v+l*(l+1)/r**2-2*e is everywhere positive
c **********************************************************************
      return
c initial value for the outward integration
c **********************************************************************
   80 db=de/dvc
      call inouh (dp,dq,dr,dq1,dfl,dv(1),z,test,nuc,nstop,jc)
      if (nstop) 310,90,310
c     nstop=45
c the expansion at the origin does not converge
c **********************************************************************
   90 nd=1
      do 110 i=1,5
      dval=dr(i)**dfl
      if (i.eq.1) go to 100
      if (dp(i-1).eq.0.0d0) go to 100
      if ((dp(i)/dp(i-1)).gt.0.0d0) go to 100
      nd=nd+1
  100 dp(i)=dp(i)*dval
      dq(i)=dq(i)*dval
      dep(i)=dep(i)*dval
  110 deq(i)=deq(i)*dval
      k=-1+2*(noeud-2*(noeud/2))
      if ((dp(1)*k).gt.0.0d0) go to 130
  120 nstop=53
c error in the expansion at the origin
c **********************************************************************
      return
  130 if ((k*nk*dq(1)).lt.0.0d0) go to 120
      dm=dpas*dkoef
c outward integration
c **********************************************************************
      do 140 i=6,imat
      dp(i)=dp(i-1)
      dq(i)=dq(i-1)
      call inth (dp(i),dq(i),dv(i),dr(i))
      if (dp(i-1).eq.0.0d0) go to 140
      if ((dp(i)/dp(i-1)).gt.0.0d0) go to 140
      nd=nd+1
      if (nd.gt.noeud) go to 150
  140 continue
      if (nd.eq.noeud) go to 160
      de=0.8d0*de
      if (de.lt.-test) go to 50
      nstop=206
c the number of nodes is too small
c **********************************************************************
      return
  150 de=1.2d0*de
      if (de.gt.elim) go to 50
      nstop=210
c the number of nodes is too big
c **********************************************************************
      return
c initial values for the inward integration
c **********************************************************************
  160 dqm=dq(imat)
      dpm=dp(imat)
      if (imm.eq.1) go to 180
      do 170 i=1,np,2
      imax=np+1-i
      if(((dv(imax)-de)*dr(imax)*dr(imax)).le.300.0d0) go to 180
  170 continue
  180 dd=sqrt(-de*(2.0d0+db/dvc))
      dpq=-dd/(dsal+db)
      dm=-dm
      do 190 i=1,5
      j=imax+1-i
      dp(j)=exp(-dd*dr(j))
      dep(i)=-dd*dp(j)*dr(j)
      dq(j)=dpq*dp(j)
  190 deq(i)=dpq*dep(i)
      m=imax-5
c inward integration
c***********************************************************************
      do 200 i=imat,m
      j=m+imat-i
      dp(j)=dp(j+1)
      dq(j)=dq(j+1)
  200 call inth (dp(j),dq(j),dv(j),dr(j))
c joining of the large components
c **********************************************************************
      dval=dpm/dp(imat)
      if (dval.gt.0.0d0) go to 210
      nstop=312
c error in the sign of the large component
c **********************************************************************
      return
  210 do 220 i=imat,imax
      dp(i)=dp(i)*dval
  220 dq(i)=dq(i)*dval
c calculation of the norm
c **********************************************************************
      dsum=3.0d0*dr(1)*(dp(1)**2+dq(1)**2)/(dpas*(dfl+dfl+1.0d0))
      do 230 i=3,imax,2
  230 dsum=dsum+dr(i)*(dp(i)**2+dq(i)**2)
     > +4.0d0*dr(i-1)*(dp(i-1)**2+dq(i-
     1 1)**2)+dr(i-2)*(dp(i-2)**2+dq(i-2)**2)
      dsum=dpas*(dsum+dr(imat)*(dqm*dqm-dq(imat)*dq(imat)))*0.3333333333
     1 333333d0
c modification of the energy
c **********************************************************************
      dbe=dp(imat)*(dqm-dq(imat))*dvc/dsum
      imm=0
      val=abs(dbe/de)
      if (val.le.test) go to 260
  240 dval=de+dbe
      if (dval.lt.0.0d0) go to 250
      dbe=dbe*0.5d0
      val=val*0.5d0
      if (val.gt.test) go to 240
      nstop=345
c energie nulle
c **********************************************************************
      return
  250 de=dval
      if (val.le.0.1d0) imm=1
      ies=ies+1
      if (ies.le.nes) go to 50
      nstop=362
c number of iterations too big
c **********************************************************************
      return
  260 dsum=sqrt(dsum)
      dq1=dq1/dsum
      do 270 i=1,imax
      dp(i)=dp(i)/dsum
  270 dq(i)=dq(i)/dsum
      do 280 i=1,4
      dpno(i,jc)=dpno(i,jc)/dsum
  280 dqno(i,jc)=dqno(i,jc)/dsum
      if (imax.eq.np) go to 300
      j=imax+1
      do 290 i=j,np
      dp(i)=0.0d0
  290 dq(i)=0.0d0
  300 nstop=0
  310 return
      end
      double precision function feff_dist(r0, r1)
c     find distance between cartesian points r0 and r1
      implicit double precision (a-h, o-z)
      dimension r0(3), r1(3)
      feff_dist = 0.0d0
      do 10  i = 1, 3
         feff_dist = feff_dist + (r0(i) - r1(i))**2
   10 continue
      feff_dist = dsqrt(feff_dist)
      return
      end
c***********************************************************************
c
c     this subroutine calculates the ' energy dependent
c     exchange-correlation potential' (or 'dirac- hara potential')
c     ref.: paper by s.h.chou, j.j.rehr, e.a.stern, e.r.davidson (1986)
c
c     inputs:    rs in a.u.
c                xk momentum in a.u.
c                vi0 constant imaginary part in rydbergs
c     outputs:   vr --- dirac potential (Hartrees)
c                vi --- constant imag part of the potential (Hartrees)
c     written by j. mustre 8/31/87
c**********************************************************************

      subroutine edp (rs, xk, vi0, vr, vi)
      implicit double precision (a-h, o-z)

      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      xf = fa / rs

c     p = sqrt (k^2 + kf^2) is the local momentum, and x = p / kf
c     Reference formula 23 in Role of Inelastic effects in EXAFS
c     by Rehr and Chou. EXAFS1 conference editted by Bianconi.
c     x is local momentum in units of fermi momentum

      x = xk / xf
      x = x + 1.0d-5
c     set to fermi level if below fermi level
      if (x .lt. 1.00001d0) x = 1.00001d0
      c = abs( (1+x) / (1-x) )
      c = log(c)
      vr = - (xf/pi) * (1.0d0 + c * (1.0d0-x**2) / (2*x))

c     Note vi=vi0/2 to have both real and imaginary part in hartrees
c     to be consistent with  other subroutines.
      vi = vi0 / 2.0d0
      return
      end
      double precision function exchan (d,dr,dexv)
      implicit double precision (a-h,o-z)
      save
c  dexv=0.0, hedin-barth corr. and exch. potential
c  dexv.ne. 0.0, dexv*slater exchange potential
c  d=4pi*rho*r^2 , radial density for r=dr
c  this function calculates exch=-r*Vexch
c  105.27578=32*(pi^2)/3
c  comments added by j. mustre 8/27/87
      if (dexv.eq.0.0d0) go to 10
      exchan=3.0d0*dexv*((dr*d/105.27578d0)**(1.0d0/3.0d0))
      return
   10 continue
      rrs=(d/(3.0d0*dr**2))**0.33333333333d0
      exchan=+0.5d0*(1.22177412d0*rrs
     >        +0.0504d0*log(30.0d0*rrs+1.0d0))*dr
      return
      end
      double precision function exchee (d,dr)
      implicit double precision (a-h,o-z)
      save
c jm if density= 0,make exchange energy equal to zero
      if (d .eq. 0.0d0) then
      exchee=0.0d0
      else
      x=(3.0d0*dr**2/d)**0.333333333333d0/30.0d0
      rx=1.0d0/x
      exchee=0.02520d0*(x**3*log(1.0d0+rx)+x*0.50d0
     > -x**2-1.0d0/3.0d0-0.2020129d0
     1 2*rx)
      endif
      return
      end
      subroutine exjlnl (z, l, jl, nl)

c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0, 1, 2 or 3 using exact analytic expression
c
c     arguments:
c       z = argument of jl and nl
c       l = integer order of spherical bessel function
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this nl = abramowitz yl.
c
c       analytic expressions from abramowitz 10.1.11 and 10.1.12
c       recurrence relation to get analytic j4,n4  eqns 10.1.19-22

      implicit double precision (a-h, o-z)

      complex*16 z, jl, nl

      complex*16 cosz, sinz

c     Exact formulae unstable for very small z, so use series
c     expansion there.  Limit of .3 chosen for 9 digit agreement.
      if (abs(z) .lt. 0.3d0)  then
         call bjnser (z, l, jl, nl, 0)
      else
c        use analytic formulae
         cosz = cos(z)
         sinz = sin(z)

         if (l .eq. 0)  then
            jl =  sinz / z
            nl = -cosz / z

         elseif (l .eq. 1)  then
            jl =  sinz/z**2 - cosz/z
            nl = -cosz/z**2 - sinz/z

         elseif (l .eq. 2)  then
            jl = ( 3.0d0/z**3 - 1.0d0/z)*sinz - 3.0d0*cosz/z**2
            nl = (-3.0d0/z**3 + 1.0d0/z)*cosz - 3.0d0*sinz/z**2

         elseif (l .eq. 3)  then
            jl = ( 15.0d0/z**4 - 6.0d0/z**2)*sinz 
     >          + (-15.0d0/z**3 + 1.0d0/z)*cosz
            nl = (-15.0d0/z**4 + 6.0d0/z**2)*cosz 
     >          + (-15.0d0/z**3 + 1.0d0/z)*sinz

         else
            stop 'exjlnl, l out of range'

         endif
      endif

      return
      end
      subroutine fermi (rhoint, vint, xmu, rs, xf)

      implicit double precision (a-h, o-z)


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


c     calculate fermi level of the system (mu) according to formula
c     mu=vcoulomb(interstitial)+vxc(interstitial)+kf(interstitial)^2
c     formula  2.13 in lee and beni, phys. rev. b15,2862(1977)

c     note that vint includes both coulomb and ground state
c     exchange-correlation potentials

c     den is the interstitial density
c     rs is the density parameter
c     xf is the interstital fermi momentum
c     xmu is the fermi level in rydbergs

      den = rhoint / (4.0d0*pi)
      rs = (3.0d0 / (4.0d0*pi*den)) ** third
      xf = fa / rs
      xmu = vint + xf**2

      return
      end
      subroutine ff2chi (ipr4, critcw, s02, tk, thetad, icsig,
     1                   vrcorr, vicorr)
c     modified for feff6l by jjr
      implicit double precision (a-h, o-z)


      character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
      common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (delk = 0.05d0)
      parameter (eps = 1.0d-10)
      parameter (eps4 = 1.0d-4)
c     e (eV) = bohr**2 * ryd * k**2 (k in invA), b2r ~=3.81
      parameter (b2r = bohr**2 * ryd)

c     This is set in dim.h for other parts of the code
      parameter (nex = 100)

c     Max number of points on fine k grid for chi output
      parameter (nfinex = 601)

      dimension achi(nex), achix(nex)
      dimension xk(nex), cdelta(nex), afeff(nex), phfeff(nex),
     1          redfac(nex), xlam(nex), rep(nex)

      dimension emxs(nex), omega(nex), xkxs(nex), xsec(nex)

      complex*16 p2, pp2
      complex*16 ck(nex), dw
      complex*16 cchi(nfinex), ccc, ccpath(nfinex)

c     head is headers from files.dat, hdxs is headers from xsect.bin
      parameter (nheadx = 30)
      character*80 head(nheadx), hdxs(nheadx)
      dimension lhead(nheadx), lhdxs(nheadx)

      parameter (nlegx = 10)
      dimension rat(3,0:nlegx), iz(0:nlegx)

      character*80  line
      parameter (nwordx = 4)
      character*50  words(nwordx), fname

c     do (or don't) correlated debye model dw factor
      logical dwcorr
c     write xmu file only if xsect.bin exists
      logical wxmu
      character*72 header
      common /header_common/ header


c     icsig 0, use real    momentum for debye waller factor
c           1, use complex momentum for debye waller factor

c     NB: code units for this module are Ang, Ang**-1, eV, etc.
      vrcorr = vrcorr * ryd
      vicorr = vicorr * ryd

      do 22  i = 1, nfinex
         cchi(i) = 0
   22 continue

c     Keep stats on total paths and paths used to make chi
      ntotal = 0
      nused = 0

c     open files.dat
      open (unit=2, file=trim(header)//'files.dat', 
     >      status='old', iostat=ios)
      call chopen (ios, trim(header)//'files.dat', 'ff2chi')
      nhead = nheadx
      call rdhead (2, nhead, head, lhead)
c     header from rdhead includes carriage control
c     skip a label line
      read(2,*)

      dwcorr = .false.
      if (tk .gt. 1.0d-1)  dwcorr = .true.

c     Open chi.dat and xmu.dat (output) and start header
      open (unit=3, file=trim(header)//'chi.dat',
     >      status='unknown', iostat=ios)
      call chopen (ios, trim(header)//'chi.dat', 'ff2chi')
c      open (unit=8, file='xsect.bin', status='old', iostat=ios)
      wxmu = .false.
      if (ios .le. 0)  wxmu = .true.
      if(wxmu) then
c        read xsect.bin
         nhdxs = nheadx
c        skip label
         edge0 = (emxs(1)/ryd + xkxs(1)**2*bohr**2)*ryd

      endif

      do 14  ihead = 1, nhead
         if (lhead(ihead) .gt. 0)  then
            write(3,12) head(ihead)(1:lhead(ihead))
         endif
   12    format ('#',a)
   14 continue
      if (dwcorr)  then
         write(3,800)  s02, tk, thetad, vfeff, vff2ch
  800    format ('# S02', f7.3, '   Temp', f8.2, '  Debye temp', f8.2,
     1           t57, 2a12)
      else
         write(3,801)  s02, vfeff, vff2ch
  801    format ('# S02', f7.3, t57, 2a12)
      endif

      if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
         write(3,802) vrcorr, vicorr
         write(77,802) vrcorr, vicorr
      endif
  802 format ('# Energy zero shift, vr, vi ', 1p, 2e14.5)


      if (critcw .gt. 0)  write(3,15) critcw
   15 format ('# Curved wave amplitude ratio filter ', f7.3, '%')
      write(3,16)
   16 format ('#    file           sig2  cw amp ratio   deg',
     1        '  nlegs  r effective')

c     Open sig2.dat if necessary (output) and start header
      if (ipr4 .ge. 1)  then
         open (unit=4, file=trim(header)//'sig2.dat',
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//'sig2.dat', 'ff2chi')
         do 514  ihead = 1, nhead
            if (lhead(ihead) .gt. 0)
     1            write(4,12) head(ihead)(1:lhead(ihead))
  514    continue
         if (dwcorr)  then
            write(4,800)  s02, tk, thetad, vfeff, vff2ch
         else
            write(4,801)  s02, vfeff, vff2ch
         endif
         write(4,16)
      endif
      write(77,515) critcw
  515 format ('    Use all paths with cw amplitude ratio', f7.2, '%')
      if (dwcorr)  then
         write(77,516) s02, tk, thetad
      else
         write(77,517) s02
      endif
  516 format('    Use correlated Debye model.  S02', f7.3,
     1        '  Temp', f8.2, '  Debye temp', f8.2)
  517 format('    Use Debye-Waller factors from files.dat.  S02', f7.3)

   10 continue
         read(2,11,end=399)  line
   11    format (a)
         call triml (line)
         nwords = nwordx
         call bwords (line, nwords, words)
c        if line was blank, skip it and go on to next line
         if (nwords .lt. 1)  goto 10

         ntotal = ntotal+1
c        word 1 - feff.dat file name
c             2 - sig2 for path
c             3 - amplitude ratio, full k range

         read(words(2),40,err=900)  sig2
         read(words(3),40,err=900)  crit
   40    format (bn, f15.0)
c        Skip un-important path

c        Write output if path is important enough (ie, path is
         if (crit .lt. critcw)  then
            write(77,17) words(1)(1:15), crit, '   (not used)  '
   17       format (4x, a, f10.4, a)
            goto 10
         endif

c        Read feff.dat file
         nused = nused+1
         write(77,17) words(1)(1:15), crit
         fname = words(1)
         open (unit=1, file=trim(header)//words(1),
     >         status='old', iostat=ios)
         call chopen (ios, trim(header)//words(1), 'ff2chi')
         nhead = nheadx
         call rdhead (1, nhead, head, lhead)
         read(1,*)  nleg, deg, reff, rs, edge
         if (abs(vrcorr) .gt. eps4) edge = edge-vrcorr
         if (nleg .gt. nlegx)  stop 'too many legs'
c        skip label
         read(1,*)
         do 30  ileg = 0, nleg-1
            read(1,*) (rat(j,ileg),j=1,3), ipot, iz(ileg)
   30    continue
c        skip label
         read(1,*)
         do 20  j = 1, 3
            rat(j,nleg) = rat(j,0)
   20    continue
         iz(nleg) = iz(0)

c        Get sig2 from correlated debye model if required
         if (dwcorr)  then
c           replace sig2 from files.dat
            call sigms (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
         endif

c        Put path into chi.dat header, sig2.dat as required
         write(3,110)  words(1)(1:15), sig2, crit,
     1                 deg, nleg, reff
         if (ipr4 .ge. 1)  then
            write(4,110)  words(1)(1:15), sig2, crit,
     1                    deg, nleg, reff
         endif
  110    format('#',1x, a, f8.5, 2f10.2, i6, f9.4)

c        read data
         i = 1
  120    read(1,*,end=130)  xk(i), cdelta(i), afeff(i),
     1             phfeff(i), redfac(i), xlam(i), rep(i)

c           make complex momentum
c           add correction to imag part of energy to xlam here

c           use atomic units for this
            viryd = vicorr / ryd
            preal = rep(i) * bohr
            xlamb = xlam(i) / bohr
            pimag = 1 / xlamb
c           p2 is p**2, pp2 is p' **2 (p prime squared, new p)
            p2 = (preal + coni*pimag)**2
            pp2 = p2 + coni*viryd
            ck(i) = sqrt (pp2)
            xlam(i) = 1 / dimag(ck(i))
            rep(i) = dble(ck(i))
c           put everything back into Ang and invAng
            ck(i) = ck(i) / bohr
            xlam(i) = xlam(i) * bohr
            rep(i) = rep(i) / bohr

            npts = i
            i = i+1
         goto 120
  130    continue
         close(unit=1)

c        Make chi, note that |feff| at k=0 is zero.  Must interpolate
c        or extrapolate to find it.  Can interpolate when we have 
c        data for k<0, but just extrapolate for all cases for now.
         iextr = 0
         do 300  i = 1, npts

c           extrapolate chi when k=0, otherwise calculate it
c           achi has no 2kr term
            dw = exp(-2*sig2*ck(i)**2)
            phdw = atan2 (dimag(dw), dble(dw))
            if (abs(xk(i)) .lt. 0.01d0)  then
               iextr = i
            else
               achi(i) = afeff(i) * deg * abs(dw) *
     1             exp(-2*reff/xlam(i)) * redfac(i) * s02 / 
     2             (abs(xk(i))*reff**2)
            endif
            achix(i) = cdelta(i) + phfeff(i) + phdw
  300    continue
c        fill in achi where extrapolation necessary
         if (iextr .gt. 0)  then
            achi(iextr) = 2*achi(iextr+1) - achi(iextr+2)
         endif

c        make sure no 2pi jumps in phase
         do 310  i = 2, npts
            call pijump (achix(i), achix(i-1))
  310    continue

c        Decide on fine grid -- need k' if vrcorr /= 0
         if (abs(vrcorr) .gt. eps4)  then
            xkpmin = xk2xkp (xk(1), vrcorr)
            n = xkpmin / delk
c           need 1st int ABOVE xkpmin/delk
            if (xkpmin .gt. 0.0d0)  n = n+1
c           First k grid point moved by vrcorr
            xkmin = n * delk
         else
c           Use unmodified grid
            xkmin = xk(1)
         endif

c        sum chi on fine k grid
         nkx = nfinex
         do 330  i = 1, nfinex
c           xkout is k value for output, xk0 is k value used for 
c           interpolation and reconstruction of chi with original grid.
c           If vrcorr=0, xk0 and xkout will be the same.
            xkout = delk * (i-1) + xkmin
            xk0 = xkp2xk (xkout, vrcorr)

c           find end of data, eps4 is to handle round-off (we've been 
c           reading files with limited precision)
            if (xk0 .gt. xk(npts)+eps4)  then
               nkx = i-1
               goto 331
            endif
            call terp (xk, achi,  npts, xk0, achi0)
            call terp (xk, achix, npts, xk0, achix0)
            cchi(i) = cchi(i) + achi0 *
     1                exp (coni * (2*xk0*reff + achix0))
            ccpath(i) = achi0 * exp (coni * (2*xk0*reff + achix0))
  330    continue
  331    continue

c        write out a chinnnn.dat for this path, if necessary.  Headers
c        later...
         if (ipr4 .ge. 2)  then
c           Assume file is form  feffnnnn.whatever, change it to
c                                chipnnnn.whatever.  Other filenames
c           will turn out wierdly
            fname(1:4) = 'chip'
            open (unit=9, file=trim(header)//fname, status='unknown')
            do 370  ihead = 1, nhead
               if (lhead(ihead) .gt. 0)  then
                  write(9,12) head(ihead)(1:lhead(ihead))
               endif
  370       continue
            if (dwcorr)  then
               write(9,800)  s02, tk, thetad, vfeff, vff2ch
            else
               write(9,801)  s02, vfeff, vff2ch
            endif

            if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
               write(9,802)  vrcorr, vicorr
            endif
            write(9,*) 'Debye-waller factor ', sig2

            write(9,407)
            write(9,338)
  338       format ('       k         chi           mag          ',
     1              'phase        phase-2kr  @#')
            do 340  i = 1, nkx
               xk0 = delk * (i-1) + xkmin
               ccc = ccpath(i)
               phase=0
               if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
               if (i .gt. 1)  call pijump (phase, phase0)
               phase0 = phase
               write(9,410)  xk0, dimag(ccc), abs(ccc), phase,
     1                       phase-2*xk0*reff
  340       continue
         endif

      goto 10
  399 continue
      close (unit=2)

c     Write it out
      write(3,405)  nused, ntotal
  405 format ('#',1x, i4, '/', i4, ' paths used')
      write(3,407)
  407 format ('#',1x, 79('-'))
      write(3,406)
  406 format ( '#     k          chi          mag           phase @#')
      do 420  i = 1, nkx
         xk0 = delk * (i-1) + xkmin
         ccc = cchi(i)
         phase=0
         if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
         if (i .gt. 1)  call pijump (phase, phase0)
         phase0 = phase
         write(3,410)  xk0, dimag(ccc), abs(ccc), phase
  410    format (1x, f10.4, 3x, 4(1pe13.6,1x))
  420 continue
      close (unit=3)


      write(77,500) nused, ntotal
  500 format (' ff2chi done, ', i4, '/', i4, ' paths used.')
      return

  900 stop 'Error reading files.dat importance factors'
      end

c     following functions use invA and eV as input and output,
c     internal workings in atomic units

      double precision function xk2xkp (xk, vrcorr)
      implicit double precision (a-h, o-z)

      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)

      xk0 = xk*bohr
      vr = vrcorr / ryd
      xksign = sign (one, xk0)
      e = xksign*xk0**2 + vr
      xk2xkp = getxk(e) / bohr
      return
      end

      double precision function xkp2xk (xkp, vrcorr)
      implicit double precision (a-h, o-z)

      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)

      xkp0 = xkp*bohr
      vr = vrcorr / ryd
      xkpsgn = sign (one, xkp0)
      e = xkpsgn*xkp0**2 - vr
      xkp2xk = getxk(e) / bohr
      return
      end
      double precision function ffq(q, ef, xk, wp, alph)
      implicit double precision (a-h,o-z)

c     input:  q, wp, alph, ef, xk
c             q is dimensionless, normalized to fermi momentum
c             xk is momentum in invBohrs
c     output: ffq only

      wq = sqrt (wp**2 + alph*q**2 + q**4)
      ffq = (wp+wq)/(q**2) + alph/(2.0d0*wp)
      ffq = ((ef*wp) / (4.0d0*xk))  * log(ffq)

      return
      end
      subroutine fixvar (rmt, edens, vtot,
     1                   vint, rhoint, nr, dx, x0, ri,
     2                   vtotph, rhoph)

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header



      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      dimension edens(nrptx), vtot (nrptx)
      dimension vtotph(nr), rhoph(nr)
      dimension ri(nr)

c     PHASE needs
c     vtot = total potential including gs xcorr, no r**2
c     edens = rho, charge density, no factor of 4*pi, no r**2
c     From overlapping, vtot = potential only, ok as is
c                       edens = density*4*pi, so fix this here.

c     If new grid is different from old one, be sure to interpolate
c     somehow...

c     Only values inside the muffin tin are used, except that XCPOT
c     (in PHASE) uses values at imt+1 and requires these to be the
c     interstitial values.  So set the last part of the arrays to
c     interstitial values...

      imt = ii(rmt)

      do 190  i = 1, imt
         vtotph(i) = vtot(i)
         rhoph(i) = edens(i)/(4.0d0*pi)
  190 continue
      do 200  i = imt+1, nrptx
         vtotph(i) = vint
         rhoph(i) = rhoint/(4.0d0*pi)
  200 continue

      return
      end
      subroutine fmtrxi (lam1x, lam2x, ie, iterm, ileg, ilegp)
      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


c     all commons except for /fmat/ are inputs

c     inputs:
c       lam1x, lam2x:  limits on lambda and lambda'
c       ie:  energy grid points
c       iterm = 1 if we're doing the termination matrix M,
c              -1 otherwise
c       ileg, ilegp: leg and leg'
c
c     Inputs from common:
c        phases, use ph(ie,...,ilegp), and lmax(ie,ilegp)
c        lambda arrays
c        rotation matrix for ilegp
c        clmz for ileg and ilegp
c        path data, eta(ilegp) and ipot(ilegp)
c        xnlm array
c
c     Output:  fmati(...,ilegp) in common /fmatrx/ is set for
c              current energy point.

c     calculate scattering amplitude matrices
c     f(lam,lam') = sum_l tl gam(l,m,n)dri(l,m,m',ileg)gamt(l,m',n')
c                 *cexp(-i*m*eta),  eta = gamma+alpha'
c     lam lt lam1x, lam' lt lam2x such that m(lam) lt l0, n(lam) lt l0
c     gam = (-)**m c_l,n+m*xnlm, gamt = (2l+1)*c_ln/xnlm,
c     gamtl = gamt*tl


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      save /nlm/
      common /nlm/ xnlm(ltot+1,mtot+1)


      common /lambda/  
     4   mlam(lamtot), 	!mu for each lambda
     5   nlam(lamtot),	!nu for each lambda
     1   lamx, 		!max lambda in problem
     2   laml0x, 	!max lambda for vectors involving absorbing atom
     3   mmaxp1, nmax 	!max mu in problem + 1, max nu in problem


      save /clmz/
      complex*16 clmi
      common /clmz/ clmi(ltot+1,mtot+ntot+1,legtot)


      complex*16 fmati
      common /fmatrx/ fmati(lamtot,lamtot,legtot)


      save /rotmat/
      common /rotmat/ dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines


      complex*16 cam, camt, cterm, tltl
      complex*16 gam(ltot+1,mtot+1,ntot+1),
     1           gamtl(ltot+1,mtot+1,ntot+1), tl

c     calculate factors gam and gamtl
      iln = 1
      ilx = lmax(ie,ipot(ilegp)) + 1
      if (iterm .gt. 0)  then
         iln = il0
         ilx = il0
      endif
      do 30  il = iln, ilx
         tltl = 2.0d0*il - 1.0d0
         if (iterm .lt. 0)  then
            tl = (exp(2.0d0*coni*ph(ie,il,ipot(ilegp))) - 1.0d0)
     >           / (2.0d0*coni)
            tltl = tltl * tl
         endif
         lam12x = max (lam1x, lam2x)
         do 20  lam = 1, lam12x
            m = mlam(lam)
            if (m .lt. 0)  goto 20
            im = m+1
            if (im .gt. il)  goto 20
            in = nlam(lam) + 1
            imn = in + m
            if (lam .gt. lam1x)  goto 10
            cam = xnlm(il,im) * (-1)**m
            if (imn .le. il)  gam(il,im,in) = cam * clmi(il,imn,ileg)
            if (imn .gt. il)  gam(il,im,in) = 0
   10       if (lam .gt. lam2x) goto 20
            camt = tltl / xnlm(il,im)
            gamtl(il,im,in) = camt * clmi(il,in,ilegp)
   20    continue
   30 continue

      do 60 lam1 = 1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do 60  lam2 = 1, lam2x
            m2 = mlam(lam2)
            in2 = nlam(lam2) + 1
            iam2 = iabs(m2) + 1
            imn1 = iam1 + in1 - 1
            cterm = 0.0d0
            ilmin = max (iam1, iam2, imn1, in2, iln)
            do 40  il = ilmin, ilx
c              skip terms with mu > l (NB il=l+1, so mu=il is mu>l)
               if (abs(m1).ge.il .or. abs(m2).ge.il)  goto 40
               m1d = m1 + mtot+1
               m2d = m2 + mtot+1

               cterm = cterm + gam(il,iam1,in1)*gamtl(il,iam2,in2)
     1                         *dri(il,m1d,m2d,ilegp)

   40       continue
            if (eta(ileg) .ne. 0.0d0) then
               m1 = mlam(lam1)
               cterm = cterm * exp(-coni*eta(ileg)*m1)
            endif
c           Above was org coding, change to use eta(ilegp) as test
c           based on algebra check.  July 20, 1992, siz&jjr
c           Changed back with redifinition of eta(see rdpath.f)
c           which is more convinient in polarization case.
c           August 8,1993, ala.
c           if (eta(ilegp) .ne. 0.0) then
c              m1 = mlam(lam1)
c              cterm = cterm * exp(-coni*eta(ilegp)*m1)
c           endif
            fmati(lam1,lam2,ilegp) = cterm
   60 continue

c     test of fmati(lam,lam',ileg)
c     plot fmat(lam,lam') = csqrt((z/2)**(m1-m2))*fmat

      return
      end
      subroutine fovrg (il, ihard, rmt, xmt, jri, e, nr, dx, ri, v, dny,
     1                  pu, qu, p, q, ps, qs, vm)
      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


c     Input:
c        il      ang mom number + 1
c        ihard   number of times convergence test fails
c        rmt     muffin tin radius
c        xmt     x such that rmt = exp ((x-1)*dx - 8.8)
c        jri     first interstitial grid point (imt + 1)
c        e       current complex energy
c        nr      number of points in r grid
c        dx      dx in Loucks' grid (usually .05)
c        ri(nr)  Loucks' position grid, r = exp ((i-1)*dx - 8.8)
c        v(nr)   total complex potential including energy dep xc
c                v is in the form  pot*r**2
c
c     Work space:
c        complex*16 p(nr), q(nr), ps(nr), qs(nr), vm(nr)
c        Must be dimensioned in calling program.  Coded like this
c        to make using different r-grids with different nrmax easy.
c
c     Output:
c        ihard   incremented each time convergence test fails
c        dny     r*g'/g, see loucks (4-85), q/p = cf/g (eq 4-86)
c        pu, qu  upper and lower components at muffin tin
c        q and q arrays  upper and lower components (see comments)

      complex*16 v(nr), e
      dimension ri(nr)
      complex*16 dny, pu, qu
      complex*16 p(nr), q(nr), ps(nr), qs(nr), vm(nr)


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)

      parameter (c = clight)
      parameter (csq = c**2)

      double precision lp1, ldcsq
      complex*16 c1,c2,c3,pc,qc,dp1,dq1,dp2,dq2,dp3,dq3,dp4,dq4
      complex*16 vh,vmh,vmnp1,psn,qsn,psnm1,qsnm1,psnm2,qsnm2
      complex*16 psnm3,qsnm3,psnm4,qsnm4,pp,qp,psnp1,qsnp1,prel,qrel
      complex*16 psu,vu,dummy
      complex*16 vn,vmn
      integer*4 n1sum

c     test=1.e+04 value in loucks
      test=1.d+05
      nrk=6

      expdxh=exp(dx/2.0d0)
      dxd4=dx/4.0d0
      dxd8=dx/8.0d0
      a1=dx*3.30d0
      a2=-dx*4.20d0
      a3=dx*7.80d0
      a4=dx*14.0d0/45.0d0
      a5=dx*64.0d0/45.0d0
      a6=dx*24.0d0/45.0d0
      call feff_diff (v,dx,jri,vm)
      twoz=-dble (v(1))/ri(1)
      l=il-1
      lp1=l+1.0d0
      ldcsq=l/csq
      ie=1
      r=ri(1)
      vn=v(1)
      vmn=vm(1)
cv    p(1)=1.0
      p(1)=1.d-20
      q(1)=-e/(2.0d0*l+3.0d0)*r*p(1)
      beta=lp1
      if (twoz.eq.0.0d0) go to 10
      beta=sqrt(lp1*l+1.0d0-(twoz/c)**2)
      sb0=(beta-lp1)*csq/twoz
      sa1=(3.0d0*beta-(twoz/c)**2)/(2.0d0*beta+1.0d0)
      sb1=csq/twoz*((beta-l)*sa1-1.0)-sb0
      sa2=((beta+3.0*lp1)*sa1-3.0d0*l+twoz/csq*(beta+lp1+3.0d0)*sb1)/
     1 (beta+1.0d0)/4.0d0
      sb2=(csq/twoz*(2.0d0*l*(beta+2.0d0-lp1)-l-(twoz/c)**2)*sa1-3.0d0*l
     1 *csq/twoz*(beta+2.0d0-lp1)
     > +(beta+3.0d0-2.0d0*lp1-(twoz/c)**2)*sb1)/
     2 (beta+1.0)/4.0d0
      delta=r*csq/twoz
      q(1)=(sb0+delta*(sb1+delta*sb2))/
     >     (1.0d0+delta*(sa1+delta*sa2))*p(1)
   10 continue
c     runge kutta method  (see loucks)
      c1=vn/r**2-e
      c2=1.0d0-c1/csq
      c3=(vmn-2.0d0*vn)/c2/c2*ldcsq
      ps(1)=r*c2*q(1)+lp1*p(1)
      qs(1)=-lp1*q(1)+(r*c1-c3/r**3)*p(1)
      n=1
   20 continue
      pc=p(n)
      qc=q(n)
      dp1=dx*(r*c2*qc+lp1*pc)
      dq1=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      pc=pc+0.50d0*dp1
      qc=qc+0.50d0*dq1
      r=r*expdxh
      vnp1=v(n+1)
      vmnp1=vm(n+1)
      vh=(vn+vnp1)*0.50d0+(vmn-vmnp1)*dxd8
      vmh=(1.50d0*(vnp1-vn)-(vmn+vmnp1)*dxd4)/dx
      c1=vh/r/r-e
      c2=1.0d0-c1/csq
      c3=(vmh-2.0d0*vh)/c2/c2*ldcsq
      dp2=dx*(r*c2*qc+lp1*pc)
      dq2=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      pc=pc+0.50d0*(dp2-dp1)
      qc=qc+0.50d0*(dq2-dq1)
      dp3=dx*(r*c2*qc+lp1*pc)
      dq3=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      pc=pc+dp3-0.50d0*dp2
      qc=qc+dq3-0.50d0*dq2
      n=n+1
      r=ri(n)
      c1=vnp1/r/r-e
      c2=1.0d0-c1/csq
      c3=(vmnp1-2.0d0*vnp1)/c2/c2*ldcsq
      dp4=dx*(r*c2*qc+lp1*pc)
      dq4=dx*(-lp1*qc+(r*c1-c3/r**3)*pc)
      p(n)=p(n-1)+(dp1+2.0d0*(dp2+dp3)+dp4)/6.0d0
      q(n)=q(n-1)+(dq1+2.0d0*(dq2+dq3)+dq4)/6.0d0
      ps(n)=r*c2*q(n)+lp1*p(n)
      qs(n)=-lp1*q(n)+(r*c1-c3/r**3)*p(n)
      vn=vnp1
      vmn=vmnp1
      n1sum=n-nrk
      if (n1sum) 20,30,30
   30 if (n.ge.jri) go to 120
      psn=ps(nrk)
      qsn=qs(nrk)
      psnm1=ps(nrk-1)
      qsnm1=qs(nrk-1)
      psnm2=ps(nrk-2)
      qsnm2=qs(nrk-2)
      psnm3=ps(nrk-3)
      qsnm3=qs(nrk-3)
      psnm4=ps(nrk-4)
      qsnm4=qs(nrk-4)
c     milne method
   40 r=ri(n+1)
      c1=v(n+1)/r/r-e
      c2=1.0d0-c1/csq
      c3=(vm(n+1)-2.0d0*v(n+1))/c2/c2*ldcsq
      pp=p(n-5)+a1*(psn+psnm4)+a2*(psnm1+psnm3)+a3*psnm2
      qp=q(n-5)+a1*(qsn+qsnm4)+a2*(qsnm1+qsnm3)+a3*qsnm2
      nit=0
   50 psnp1=r*c2*qp+lp1*pp
      qsnp1=-lp1*qp+(r*c1-c3/r**3)*pp
      pc=p(n-3)+a4*(psnp1+psnm3)+a5*(psn+psnm2)+a6*psnm1
      qc=q(n-3)+a4*(qsnp1+qsnm3)+a5*(qsn+qsnm2)+a6*qsnm1
      if (abs(test*(pc-pp))-abs(pc)) 60,60,70
   60 if (abs(test*(qc-qp))-abs(qc)) 110,110,70
   70 n1sum=nit-40
      if (n1sum) 100,80,100
c  70 if (nit-5) 100,80,100 value in loucks
   80 prel=(pc-pp)/pc
      qrel=(qc-qp)/qc
c     count times hard test fails
      ihard = ihard + 1
c     print90, il,ie,n,prel,qrel
   90 format (' hard test in fovrg il=',i2,' ie=',i1,' n=',i3,' prel='
     1 ,e16.8,' qrel=',e16.8,' **********')
      go to 110
  100 nit=nit+1
      pp=pc
      qp=qc
      go to 50
  110 n=n+1
      p(n)=pc
      q(n)=qc
      ps(n)=psnp1
      qs(n)=qsnp1
      psnm4=psnm3
      psnm3=psnm2
      psnm2=psnm1
      psnm1=psn
      psn=psnp1
      qsnm4=qsnm3
      qsnm3=qsnm2
      qsnm2=qsnm1
      qsnm1=qsn
      qsn=qsnp1
c     introduce scale factor to prevent overflow on vax jjr
      if(abs(pc).lt.1.d+20) go to 119
      scale=1.d-20
      do 112 mm=1,6
      nm=n-mm+1
      p(nm)=scale*p(nm)
      q(nm)=scale*q(nm)
      ps(nm)=scale*ps(nm)
      qs(nm)=scale*qs(nm)
  112 continue
      psnm4=scale*psnm4
      psnm3=scale*psnm3
      psnm2=scale*psnm2
      psnm1=scale*psnm1
      psn=scale*psn
      qsnm4=scale*qsnm4
      qsnm3=scale*qsnm3
      qsnm2=scale*qsnm2
      qsnm1=scale*qsnm1
      qsn=scale*qsn
  119 n1sum=n-jri
      if (n1sum) 40,120,120
  120 jm=jri-1
      x=dx*(xmt-jm)
      call intpol (zero,dx,p(jm),p(jri),ps(jm),ps(jri),x,pu,psu)
      call intpol (zero,dx,q(jm),q(jri),qs(jm),qs(jri),x,qu,dummy)
      call intpol (zero,dx,v(jm),v(jri),vm(jm),vm(jri),x,vu,dummy)
      dny=rmt*(1.0-(vu/rmt**2-e)/csq)*qu/pu+l
c dny is r*g'/g, see loucks (4-85), q/p = cf/g (eq 4-86)
c (watch for factors of rmt)
      return
      end
      double precision function fpot(r,z,wa)
      implicit double precision (a-h,o-z)
      save
c
c thomas fermi potential at the point r; z=atomic number
c wa=number of electrons-z-1
c **********************************************************************
      wc=sqrt((r*(z+wa)**(1.0d0/3.0d0))/0.88530d0)
      wd=wc*(0.601120d0*wc+1.810610d0)+1.0d0
      we=wc*(wc*(wc*(wc*(0.04793d0*wc+0.21465d0)+0.77112d0)+1.39515d0)
     > +1.81061d0)+1.0d0
      wc=(z+wa)*(wd/we)**2-wa
      fpot=-wc/r
      return
      end
      subroutine frnrm (rho, iz, rnrm)
      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header

      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      dimension rho(nrptx)

      real*8 sum,fr,fl

c     finds norman radius

c     Need overlapped densities.  We'll get them in the form
c     4*pi*density = rho.  Also need z of atom

c     Then integrate out to the point where the integral of
c     4*pi*density*r**2 is equal to iz
      sum = 0.0d0
      do 10  i = 1, nrptx-1
         fr = rho(i+1) * rr(i+1)**3
         fl = rho(i)   * rr(i)**3
         sumsav = sum
         sum = sum + 0.025d0*(fr+fl)
         if (sum .ge. iz)  then
            inrm = i+1
            goto 20
         endif
   10 continue
      write(77,*) ' FRNRM Could not integrate enough charge to reach'
      write(77,*) '       required z.'
      write(77,*) "error sum,iz=",sum,iz
      stop 'FRNRM-1'
   20 continue
c     inrm is too big, subtract one from irnm and interpolate
c     to get correct value
      inrm = inrm - 1
      deltaq = iz - sumsav
      fr = rho(inrm+1) * rr(inrm+1)**3
      fl = rho(inrm)   * rr(inrm)**3
c     dipas is delta i * 0.05
      dipas = 2*deltaq / (fl + fr)
      rnrm = rr(inrm)*(1.0d0 + dipas)

      return
      end
      subroutine genfmt (ipr3, critcw, sig2g, iorder)
      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      save /clmz/
      complex*16 clmi
      common /clmz/ clmi(ltot+1,mtot+ntot+1,legtot)


      complex*16 fmati
      common /fmatrx/ fmati(lamtot,lamtot,legtot)


      common /lambda/  
     4   mlam(lamtot), 	!mu for each lambda
     5   nlam(lamtot),	!nu for each lambda
     1   lamx, 		!max lambda in problem
     2   laml0x, 	!max lambda for vectors involving absorbing atom
     3   mmaxp1, nmax 	!max mu in problem + 1, max nu in problem


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines


      save /nlm/
      common /nlm/ xnlm(ltot+1,mtot+1)


      save /rotmat/
      common /rotmat/ dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)



      character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
      common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch


c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola


      complex*16  rho(legtot), pmati(lamtot,lamtot,2)
      complex*16  pllp, ptrac, srho, prho, cdel1, cfac
      complex*16  cchi(nex), cfms, mmati
      dimension   mmati(-mtot:mtot,-mtot:mtot)
      dimension   t3j(-mtot-1:mtot+1,-1:1)
      dimension   xk(nex), ckmag(nex)
      complex*16  ck(nex)
      dimension   ffmag(nex)

      character*12 fname

      logical done

c     Input flags:
c     iorder, order of approx in f-matrix expansion (see setlam)
c             (normal use, 2.  Do ss exactly regardless of iorder)

c     used for divide-by-zero and trig tests
      parameter (eps = 1.0d-16)

c     Read phase calculation input, data returned via commons
      open (unit=1, file=trim(header)//'phase.bin', status='old',
     1      access='sequential', form='unformatted', iostat=ios)
      call chopen (ios, trim(header)//'phase.bin', 'genfmt')
      call rphbin (1)
      close (unit=1)

c     Open path input file (unit in) and read title.  Use unit 1.
      ntitle = 5
      open (unit=1, file=trim(header)//'paths.dat',
     >      status='old', iostat=ios)
      call chopen (ios, trim(header)//'paths.dat', 'genfmt')
      call rdhead (1, ntitle, title, ltitle)
      if (ntitle .le. 0)  then
         title(1) = ' '
         ltitle(1) = 1
      endif

c     cgam = gamma in mean free path calc (eV).  Set to zero in this
c     version.  Set it to whatever you want if you need it.
c     cgam = 0
c     cgam = cgam / ryd
c     add cnst imag part to eref
c     do 20  ie = 1, ne
c        eref(ie) = eref(ie) - coni*cgam/2
c  20 continue

   50 format (a)
   60 format (1x, a)
   70 format (1x, 79('-'))

c     Save filenames of feff.dat files for use by ff2chi
      open (unit=2, file=trim(header)//'files.dat',
     >      status='unknown', iostat=ios)
      call chopen (ios, trim(header)//'files.dat', 'genfmt')
c     Put phase header on top of files.dat
      do 100  itext = 1, ntext
         write(2,60)  text(itext)(1:ltext(itext))
  100 continue
      write(2,70)
      write(2,120)
  120 format ('    file        sig2   amp ratio    ',
     1        'deg    nlegs  r effective')

c     Set crit0 for keeping feff.dat's
      if (ipr3 .le. 0)  crit0 = 2*critcw/3
c     Make a header for the running messages.
      write(77,130) critcw
  130 format ('    Curved wave chi amplitude ratio', f7.2, '%')
      if (ipr3 .le. 0)  write(77,131) crit0
  131 format ('    Discard feff.dat for paths with cw ratio <',
     1         f7.2, '%')
      write(77,132)
  132 format ('    path  cw ratio     deg    nleg  reff')

c     Set nlm factors in common /nlm/ for use later
      call snlm (ltot+1, mtot+1)

      if (pola) then
c        Make 3j factors in t3j  (multiplied by sqrt(3*(2l0+1)) for
c        further convinience - the same expression for chi)
c        l0 - final momentum, initial momentum = l0-1.
         do 140  m0 = -l0+1,l0-1
            t3j(m0, 1) = (-1)**(l0+1+m0)*sqrt(3.0d0*(l0+m0)*(l0+m0+1)
     1                /(2*l0)/(2*l0-1))
            t3j(m0, 0) = (-1)**(l0+m0)*sqrt(3.0d0*(l0*l0-m0*m0)/
     1                l0/(2*l0-1))
  140    continue
         do 145  m0 = -l0+1,l0-1
            t3j(m0,-1) = t3j(-m0,1)
  145    continue
      endif

c     While not done, read path, find feff.
      open (unit=4,file=trim(header)//'nstar.dat',
     >      status='unknown', iostat=ios)
      write(4,198, iostat=ios) evec
  198 format('polarization  ',3f8.4)
      write(4,199, iostat=ios)
  199 format('npath  nstar')
      npath = 0
      ntotal = 0
      nused = 0
      xportx = -1
  200 continue

c        Read current path
         call rdpath (1, pola, done,xstar)
         icalc = iorder
         if (done)  goto  1000
         npath = npath + 1
         ntotal = ntotal + 1

         write (4,201,iostat=ios) npath, xstar
  201    format (i5, f8.4)

c        Need reff
         reff = 0
         do 220  i = 1, nleg
            reff = reff + ri(i)
  220    continue
         reff = reff/2

c        Set lambda for low k
         call setlam (icalc, 1)

c        Calculate and store rotation matrix elements
c        Only need to go to (il0, il0, ...) for isc=nleg and
c        nleg+1 (these are the paths that involve the 'z' atom
         call rot3i (il0, il0, nleg)
         do 400  isc = 1, nsc
            call rot3i (lmaxp1, mmaxp1, isc)
  400    continue
         if (pola) then
c           one more rotation in polarization case
            call rot3i (il0, il0, nleg+1)
            call mmtr(t3j,mmati)
         endif 


c        Big energy loop
         do 800  ie = 1, ne

c           real momentum (k)
            xk(ie) = getxk (em(ie) - edge)

c           complex momentum (p)
            ck(ie) = sqrt (em(ie) - eref(ie))
            ckmag(ie) = abs(ck(ie))
c           complex rho
            do 420  ileg = 1, nleg
               rho(ileg) = ck(ie) * ri(ileg)
  420       continue

c           if ck is zero, xafs is undefined.  Make it zero and jump
c           to end of calc part of loop.
            if (abs(ck(ie)) .le. eps)  then
               cchi(ie) = 0
               goto 620
            endif

c           Calculate and store spherical wave factors c_l^(m)z^m/m!
c           in a matrix clmi(il,im,ileg), ileg=1...nleg.
c           Result is that common /clmz/ is updated for use by fmtrxi.

c           zero clmi arrays
            do 440  ileg = 1, legtot
            do 440  il = 1, ltot+1
            do 440  im = 1, mtot+ntot+1
               clmi(il,im,ileg) = 0
  440       continue

            mnmxp1 = mmaxp1 + nmax

            lxp1 = max (lmax(ie,ipot(1))+1, l0+1)
            mnp1 = min (lxp1, mnmxp1)
            call sclmz (rho, lxp1, mnp1, 1)

            lxp1 = max (lmax(ie,ipot(nsc))+1, l0+1)
            mnp1 = min (lxp1, mnmxp1)
            call sclmz (rho, lxp1, mnp1, nleg)

            do 460  ileg = 2, nleg-1
               isc0 = ileg-1
               isc1 = ileg
               lxp1 = max (lmax(ie,ipot(isc0))+1, lmax(ie,ipot(isc1))+1)
               mnp1 = min (lxp1, mnmxp1)
               call sclmz (rho, lxp1, mnp1, ileg)
  460       continue

c           Calculate and store scattering matrices fmati.

            if (pola) then
c              Polarization version, make new m matrix
c              this will fill fmati(...,nleg) in common /fmtrxi/
               call mmtrxi (laml0x, mmati, ie, 1, nleg)
            else 
c              Termination matrix, fmati(...,nleg)
               iterm = 1
               call fmtrxi (laml0x, laml0x, ie, iterm, 1, nleg)
            endif

            iterm = -1
c           First matrix
            call fmtrxi (lamx, laml0x, ie, iterm, 2, 1)
c           Last matrix if needed
            if (nleg .gt. 2)  then
               call fmtrxi (laml0x, lamx, ie, iterm, nleg, nleg-1)
            endif
c           Intermediate scattering matrices
            do 480  ilegp = 2, nsc-1
               ileg = ilegp + 1
               call fmtrxi (lamx, lamx, ie, iterm, ileg, ilegp)
  480       continue

c           Big matrix multiplication loops.
c           Calculates trace of matrix product
c           M(1,N) * f(N,N-1) * ... * f(3,2) * f(2,1), as in reference.
c           We will (equivalently) calculate the trace over lambda_N of
c           f(N,N-1) * ... * f(3,2) * f(2,1) * M(1,N), working from
c           right to left.
c           Use only 2 pmati arrays, alternating indp (index p)
c           1 and 2.

c           f(2,1) * M(1,N) -> pmat(1)
            indp = 1
            do 520  lmp = 1, laml0x
            do 520  lm = 1, lamx
               pllp = 0
               do 500  lmi = 1, laml0x
                  pllp = pllp + fmati(lm,lmi,1) * fmati(lmi,lmp,nleg)
  500          continue
               pmati(lm,lmp,indp)=pllp
  520       continue

c           f(N,N-1) * ... * f(3,2) * [f(2,1) * M(1,N)]
c           Term in [] is pmat(1)
            do 560 isc = 2, nleg-1
c              indp is current p matrix, indp0 is previous p matrix
               indp = 2 - mod(isc,2)
               indp0 = 1 + mod(indp,2)
               do 550  lmp = 1, laml0x
               do 550  lm = 1, lamx
                  pllp=0
                  do 540 lmi = 1, lamx
                     pllp = pllp +
     1                      fmati(lm,lmi,isc)*pmati(lmi,lmp,indp0)
  540             continue
  550          pmati(lm,lmp,indp) = pllp
  560       continue

c           Final trace over matrix
            ptrac=0
            do 580  lm = 1, laml0x
               ptrac = ptrac + pmati(lm,lm,indp)
  580       continue

c           Calculate xafs
c           srho=sum pr(i), prho = prod pr(i)
            srho=0
            prho=1
            do 600  ileg = 1, nleg
               srho = srho + rho(ileg)
               prho = prho * rho(ileg)
  600       continue
c           Complex chi (without 2kr term)
c           ipot(nleg) is central atom
            cdel1 = exp(2*coni*ph(ie,il0,ipot(nleg)))
            cfac = cdel1 * exp(coni*(srho-2*xk(ie)*reff)) / prho

            cchi(ie) = ptrac * cfac/(2*l0+1)

c           When ck(ie)=0, xafs is set to zero.  Calc above undefined.
c           Jump to here from ck(ie)=0 test above.
  620       continue

c        end of energy loop
  800    continue

c        Make importance factor, deg*(integral (|chi|*d|p|))
c        make ffmag (|chi|)
c        xport   importance factor
         do 810  ie = 1, ne
               ffmag(ie) = abs(cchi(ie))
  810    continue

c        integrate from edge (ik0) to ne
         nemax = ne - ik0 + 1
         call feff_trap (ckmag(ik0), ffmag(ik0), nemax, xport)
         xport = abs(deg*xport)
         if (xport .gt. xportx)  xportx = xport
         crit = 100 * xport / xportx

c        Write output if path is important enough (ie, path is

c        Write feff.dat if we need it.
         if (ipr3 .ge. 1  .or.  crit .ge. crit0)  then
c           Prepare output file feffnnnn.dat (unit 3)
            write(fname,241)  ipath
  241       format ('feff', i4.4, '.dat')
            open (unit=3, file=trim(header)//fname,
     >            status='unknown', iostat=ios)
            call chopen (ios, trim(header)//fname, 'genfmt')
c           put header on feff.dat
            do 245  itext = 1, ntext
               write(3,60)  text(itext)(1:ltext(itext))
  245       continue
            write(3,250) ipath, icalc, vfeff, vgenfm
  250       format (' Path', i5, '      icalc ', i7, t57, 2a12)
            write(3,70)
            write(3,290)  nleg, deg, reff*bohr, rnrmav, edge*ryd
  290       format (1x, i3, f8.3, f9.4, f10.4, f11.5, 
     1              ' nleg, deg, reff, rnrmav(bohr), edge')
            write(3,300)
  300       format ('        x         y         z   pot at#')
            write(3,310)  (rat(j,nleg)*bohr,j=1,3), ipot(nleg),
     1                    iz(ipot(nleg)), potlbl(ipot(nleg))
  310       format (1x, 3f10.4, i3, i4, 1x, a6, '   absorbing atom')
            do 330  ileg = 1, nleg-1
               write(3,320)  (rat(j,ileg)*bohr,j=1,3), ipot(ileg),
     1                       iz(ipot(ileg)), potlbl(ipot(ileg))
  320          format (1x, 3f10.4, i3, i4, 1x, a6)
  330       continue

            write(3,340)
  340       format    ('    k   real[2*phc]   mag[feff]  phase[feff]',
     1                 ' red factor   lambda      real[p]@#')

c           Make the feff.dat stuff and write it to feff.dat
            do 900  ie = 1, ne
c              Consider chi in the standard XAFS form.  Use R = rtot/2.
               xlam = 1.0d10
               if (dabs(dimag(ck(ie))) .gt. eps) 
     >            xlam = 1.0d0/dimag(ck(ie))
               redfac = exp(-2 * dimag (ph(ie,il0,ipot(nleg))))
               cdelt = 2*dble(ph(ie,il0,ipot(nleg)))
               cfms = cchi(ie) * xk(ie) * reff**2 *
     1              exp(2*reff/xlam) / redfac
               if (abs(cchi(ie)) .lt. eps)  then
                  phff = 0
               else
                  phff = atan2(dimag(cchi(ie)), dble(cchi(ie)))
               endif
c              remove 2 pi jumps in phases
               if (ie .gt. 1)  then
                  call pijump (phff, phffo)
                  call pijump (cdelt, cdelto)
               endif
               phffo = phff
               cdelto = cdelt

c              write 1 k, momentum wrt fermi level k=sqrt(p**2-kf**2)
c                    2 central atom phase shift (real part),
c                    3 magnitude of feff,
c                    4 phase of feff,
c                    5 absorbing atom reduction factor,
c                    6 mean free path = 1/(Im (p))
c                    7 real part of local momentum p

               write(3,640)
     1            xk(ie)/bohr,
     2            cdelt + l0*pi,
     3            abs(cfms) * bohr,
     4            phff - cdelt - l0*pi,
     5            redfac,
     6            xlam * bohr,
     7            dble(ck(ie))/bohr
  640          format (1x, f6.3, 1x, 3(1pe11.4,1x),0pe11.4,1x,
     1                               2(1pe11.4,1x))
  900       continue

c           Done with feff.dat
            close (unit=3)

c           Put feff.dat and stuff into files.dat
            write(2,820) fname, sig2g, crit, deg,
     1                   nleg, reff*bohr
  820       format(1x, a, f8.5, 2f10.3, i6, f9.4)

c           Tell user about the path we just did
            write(77,210) ipath, crit, deg, nleg, reff*bohr
  210       format (3x, i4, 2f10.3, i6, f9.4)
            nused = nused+1

         else
c           path unimportant, tell user
            write(77,211) ipath, crit, deg, nleg, reff*bohr
  211       format (3x, i4, 2f10.3, i6, f9.4, ' neglected')
         endif

c        Do next path
         goto 200

c     Done with loop over paths
 1000 continue
c     close paths.dat, files.dat
      close (unit=1)
      close (unit=2)
      close (unit=4)
      write(77,1010) nused, ntotal
 1010 format (1x, i4, ' paths kept, ', i4, ' examined.')

      return
      end
      subroutine getorb (iz, ihole, ion, norb, norbco,
     1                  den, nqn, nk, nel)

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


c     Save internal variables in case this gets re-entered
      save

c     Gets orbital data for chosen element.  Input is iz, atomic number
c     of desired element, other arguments are output.

c     Written by Steven Zabinsky, July 1989
c
c     last modified (20 aug 1989)  table increased to at no 95

c     Table for each element has occupation of the various levels.
c     The order of the levels in each array is:

c     element  level     principal qn (nqn), kappa qn (nk)
c           1  1s        1  -1
c           2  2s        2  -1
c           3  2p1/2     2   1
c           4  2p3/2     2  -2
c           5  3s        3  -1
c           6  3p1/2     3   1
c           7  3p3/2     3  -2
c           8  3d3/2     3   2
c           9  3d5/2     3  -3
c          10  4s        4  -1
c          11  4p1/2     4   1
c          12  4p3/2     4  -2
c          13  4d3/2     4   2
c          14  4d5/2     4  -3
c          15  4f5/2     4   3
c          16  4f7/2     4  -4
c          17  5s        5  -1
c          18  5p1/2     5   1
c          19  5p3/2     5  -2
c          20  5d3/2     5   2
c          21  5d5/2     5  -3
c          22  5f5/2     5   3
c          23  5f7/2     5  -4
c          24  6s        6  -1
c          25  6p1/2     6   1
c          26  6p3/2     6  -2
c          27  6d3/2     6   2
c          28  6d5/2     6  -3
c          29  7s        7  -1

      dimension den(30), nqn(30), nk(30), nel(30)
      dimension kappa (29)
      dimension iocc (95, 29)
      dimension nnum (29)
c     dimension ncore(95)

c     kappa quantum number for each orbital
c     k = - (j + 1/2)  if l = j - 1/2
c     k = + (j + 1/2)  if l = j + 1/2
      data kappa /-1,-1, 1,-2,-1,   1,-2, 2,-3,-1,   1,-2, 2,-3, 3,
     1            -4,-1, 1,-2, 2,  -3, 3,-4,-1, 1,  -2, 2,-3,-1/

c     principal quantum number (energy eigenvalue)
      data nnum  /1,2,2,2,3,  3,3,3,3,4,  4,4,4,4,4,
     1            4,5,5,5,5,  5,5,5,6,6,  6,6,6,7/

c     number of core orbitals for z = 1 to 95
c     data ncore
c    1  /0, 0, 1, 1, 1,  1, 1, 1, 1, 1,  4, 4, 4, 4, 4,  4, 4, 4, 4, 4,
c    2   4, 4, 4, 4, 4,  4, 4, 4, 9, 9,  9, 9, 9, 9, 9,  9, 9, 9, 9, 9,
c    3   9, 9, 9, 9, 9,  9, 9, 9, 9, 9,  9, 9, 9, 9, 9,  9, 9, 9, 9, 9,
c    4   9, 9, 9, 9, 9,  9, 9, 9, 9, 9, 16,16,16,16,16, 16,16,16,16,16,
c    5  16,16,16,16,16, 16,16,16,16,16, 16,16,16,16,16/

c     occupation of each level for z = 1, 95
      data (iocc( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 3,i),i=1,29)  /2,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 4,i),i=1,29)  /2,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 5,i),i=1,29)  /2,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 6,i),i=1,29)  /2,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 7,i),i=1,29)  /2,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 8,i),i=1,29)  /2,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 9,i),i=1,29)  /2,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(10,i),i=1,29)  /2,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(11,i),i=1,29)  /2,2,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(12,i),i=1,29)  /2,2,2,4,2,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(13,i),i=1,29)  /2,2,2,4,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(14,i),i=1,29)  /2,2,2,4,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(15,i),i=1,29)  /2,2,2,4,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(16,i),i=1,29)  /2,2,2,4,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(17,i),i=1,29)  /2,2,2,4,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(18,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(19,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(20,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(21,i),i=1,29)  /2,2,2,4,2,  2,4,1,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(22,i),i=1,29)  /2,2,2,4,2,  2,4,2,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(23,i),i=1,29)  /2,2,2,4,2,  2,4,3,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(24,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(25,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(26,i),i=1,29)  /2,2,2,4,2,  2,4,4,2,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(27,i),i=1,29)  /2,2,2,4,2,  2,4,4,3,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(28,i),i=1,29)  /2,2,2,4,2,  2,4,4,4,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(29,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(30,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(31,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(32,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(33,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(34,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(35,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(36,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(37,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(38,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(39,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(40,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(41,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(42,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(43,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(44,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(45,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(46,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(47,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(48,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(49,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(50,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(51,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(52,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(53,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(54,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(55,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (iocc(56,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(57,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(58,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,2,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(59,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,3,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(60,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,4,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(61,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,5,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(62,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(63,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(64,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(65,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           3,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(66,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           4,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(67,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           5,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(68,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           6,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(69,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           7,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(70,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(71,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(72,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (iocc(73,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  0,0,0,2,0,  0,0,0,0/
      data (iocc(74,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  0,0,0,2,0,  0,0,0,0/
      data (iocc(75,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  1,0,0,2,0,  0,0,0,0/
      data (iocc(76,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  2,0,0,2,0,  0,0,0,0/
      data (iocc(77,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  3,0,0,2,0,  0,0,0,0/
      data (iocc(78,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  5,0,0,1,0,  0,0,0,0/
      data (iocc(79,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,1,0,  0,0,0,0/
      data (iocc(80,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,0,  0,0,0,0/
      data (iocc(81,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,1,  0,0,0,0/
      data (iocc(82,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  0,0,0,0/
      data (iocc(83,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  1,0,0,0/
      data (iocc(84,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  2,0,0,0/
      data (iocc(85,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  3,0,0,0/
      data (iocc(86,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,0/
      data (iocc(87,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,1/
      data (iocc(88,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,2/
      data (iocc(89,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,1,0,2/
      data (iocc(90,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,2,0,2/
      data (iocc(91,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,2,0,2,2,  4,1,0,2/
      data (iocc(92,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,3,0,2,2,  4,1,0,2/
      data (iocc(93,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,4,0,2,2,  4,1,0,2/
      data (iocc(94,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,0,2,2,  4,0,0,2/
      data (iocc(95,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,1,2,2,  4,0,0,2/

      if (iz .lt. 1  .or.  iz .gt. 95)  then
         write(77,*)  ' Atomic number ', iz, ' not available.'
         stop
      endif

      index = iz - ion
      if (ihole .gt. 0)  then
         index = index + 1
c        remove an electron from the level specified by ihole
         if (iocc(index,ihole) .lt. 1)  then
            write(77,*) ' Cannot remove an electron from this level'
            stop 'GETORB-1'
         endif
         iocc(index,ihole) = iocc(index,ihole) - 1
      endif

      norb = 0
      do 10  i = 1, 29
         if (iocc(index,i) .ne. 0)  then
            norb = norb + 1
            nqn(norb) = nnum(i)
            nk(norb)  = kappa(i)
            nel(norb) = iocc(index,i)
            den(norb) = 0.0d0
         endif
   10 continue

c     restore iocc array for neatness
      if (ihole .gt. 0)  then
         iocc(index,ihole) = iocc(index,ihole) + 1
      endif

      norbco = norb

      return
      end
      double precision function getxk(e)
      implicit double precision (a-h, o-z)

c     Make xk from energy as
c          k =  sqrt( e)  for e > 0  (above the edge)
c          k = -sqrt(-e)  for e < 0  (below the edge)

      getxk = sqrt(abs(e))
      if (e .lt. 0.0d0)  getxk = - getxk
      return
      end
      subroutine sthead (ntitle, title, ltitle, nph, iz, rmt, rnrm,
     1                  ion, ifrph, ihole, ixc,
     2                  vr0, vi0, rs0, gamach, xmu, xf, vint, rs,
     3                  nhead, lhead, head)

c     SeT HEAD
c     This routine makes the file header, returned in head array.
c     header lines do not include a leading blank.
c     Last header line is not --------- end-of-header line

c     title lines coming into sthead include carriage control, since
c     they were read from potph.dat

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      dimension ifrph(0:nphx)
      dimension ion(0:nfrx)
      dimension iz(0:nfrx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

      character*80 title(ntitle)
      parameter (nheadx = 30)
      character*80 head(nheadx)
      dimension lhead(nheadx), ltitle(ntitle)

      character*80 heada(nheadx)
      dimension lheada(nheadx)
      save nheada, lheada, heada
c     heada, etc., are saved for use by entry wthead

      character*10 shole(0:9)
      character*8  sout(0:6)
      common /labels/ shole, sout

      character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
      common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch

c     character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
c     common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch

c     FiLl head array with HEADer
c     Fills head arrray, n = number of lines used.
c     Does not include line of dashes at the end.

      nhead = 1
      if (ntitle .ge. 1  .and.  ltitle(1).gt.1)  then
         write(head(nhead),100)  title(1)(2:), vfeff, vpotph
      else
         write(head(nhead),102)  vfeff, vpotph
      endif
  100 format(a55, t56, 2a12)
  102 format(t56, 2a12)
      do 120  ititle = 2, ntitle
         if (ltitle(ititle).le.1)  goto 120
         nhead = nhead+1
         write(head(nhead),110) title(ititle)(2:)
  110    format(a79)
  120 continue
      if (ion(0) .ne. 0)  then
         nhead = nhead+1
         write(head(nhead),130)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, ion(0), shole(ihole)
      else
         nhead = nhead+1
         write(head(nhead),140)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, shole(ihole)
      endif
  130 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',i2,1x,a10)
  140 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3, 1x,a10)

      do 150  iph = 1, nph
         ifr = ifrph(iph)
         if (ion(ifr) .ne. 0)  then
            nhead = nhead+1
            write(head(nhead),160)  iph, iz(ifr),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr, ion(ifr)
         else
            nhead = nhead+1
            write(head(nhead),170)  iph, iz(ifr),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr
         endif
  150 continue
  160 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',i2)
  170 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3)
      if (abs(vi0) .gt. 1.0d-8 .or. abs(vr0) .gt. 1.0d-8)  then
         nhead = nhead+1
         write(head(nhead),180)  gamach*ryd, sout(ixc), vi0*ryd,
     1                           vr0*ryd
      else
         nhead = nhead+1
         write(head(nhead),190)  gamach*ryd, sout(ixc)
      endif
      nhead = nhead+1
  180 format('Gam_ch=',1pe9.3, 1x,a8, ' Vi=',1pe10.3, ' Vr=',1pe10.3)
  190 format('Gam_ch=',1pe9.3, 1x,a8)
  200 format('Mu=',1pe10.3, ' kf=',1pe9.3, ' Vint=',1pe10.3,
     x        ' Rs_int=',0pf6.3)
      write(head(nhead),200)  xmu*ryd, xf/bohr, vint*ryd, rs
      if (ixc .eq. 4)  then 
          nhead = nhead+1
          write(head(nhead),210)  rs0
  210     format ('Experimental DH-HL exch, rs0 = ', 1pe14.6)
      endif
      do 220  i = 1, nhead
         lhead(i) = istrln(head(i))
         heada(i) = head(i)
         lheada(i) = lhead(i)
  220 continue
      nheada = nhead

      return

      entry wthead (io)
c     Dump header to unit io, which must be open.  Add carraige control
c     to head array, which doesn't have it.

      do 310 i = 1, nheada
         ll = lheada(i)
         write(io,300)  heada(i)(1:ll)
  300    format (1x, a)
  310 continue
      end
c     These heap routines maintain a heap (array h) and an index
c     array (array ih) used to keep other data associated with the heap
c     elements.

      subroutine hup (h, ih, n)
      implicit double precision (a-h, o-z)
c     heap is in order except for last element, which is new and must
c     be bubbled through to its proper location
c     new element is at i, j = index of parent
      integer  n,i,j
      integer  ih(n)
      dimension h(n)


      i = n

   10 j = i/2
c     if no parent, we're at the top of the heap, and done
      if (j .eq. 0)  return
      if (h(i) .lt. h(j))  then
         call swapfeff (h(i), h(j))
         call iswapfeff (ih(i), ih(j))
         i = j
         goto 10
      endif
      return
      end

      subroutine hdown (h, ih, n)
      implicit double precision (a-h, o-z)
c     h is in order, except that 1st element has been replaced.
c     Bubble it down to its proper location.  New element is i,
c     children are j and k.

      integer  n,i,j,k
      integer  ih(n)
      dimension h(n)

      i = 1

   10 continue
      j = 2*i
      k = j + 1

c     if j > n, new element is at bottom, we're done
      if (j .gt. n)  return
c     handle case where new element has only one child
      if (k .gt. n)  k = j

      if (h(j) .gt. h(k))  j = k
c     j is now index of smallest of children

      if (h(i) .gt. h(j))  then
         call swapfeff (h(i), h(j))
         call iswapfeff (ih(i), ih(j))
         i = j
         goto 10
      endif

      return
      end

      subroutine swapfeff (a, b)
      implicit double precision (a-h, o-z)
      t = a
      a = b
      b = t
      return
      end

      subroutine iswapfeff (i, j)
      implicit double precision (a-h, o-z)
      integer  i,j,k
      k = i
      i = j
      j = k
      return
      end
      subroutine imhl (rs, xk, eim, icusp)
      implicit double precision (a-h,o-z)

c     what is xk?  k**2 - mu + kf**2?

c written by j. mustre (march 1988)
c code is based on analytical expression derived by john rehr.
c it leaves the real part, calculated in rhl unchanged.
c
c modified by j. rehr  (oct 1991) - adds quinn approximation for
c losses due to electron-hole pairs below the plasmon turn on
c see new subroutine quinn.f, which incorporates r. albers coding of
c j.j. quinn's approximations for details.


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)

c     alph is Hedin-Lundquist parameter
      parameter (alph = 4.0d0 / 3.0d0)
      external ffq

      integer icount
      save icount
      data icount /0/

      icusp=0
      xf = fa / rs
      ef = xf**2 / 2

c     xk0 is xk normalized by k fermi.
      xk0 = xk/xf
c     set to fermi level if below fermi level
      if (xk0 .lt. 1.00001d0) then
         xk0 = 1.00001d0
      endif

c     wp is given in units of the fermi energy in the formula below.
      wp = sqrt (3 / rs**3) / ef
      xs = wp**2 - (xk0**2 - 1)**2

      eim = 0
      if (xs .lt. 0.0d0)  then
         q2 = sqrt ( (sqrt(alph**2-4*xs) - alph) / 2 )
         qu = min (q2, (1+xk0))
         d1 = qu - (xk0 - 1)
         if (d1 .gt. 0)  then
            eim = ffq (qu,ef,xk,wp,alph) - ffq (xk0-1,ef,xk,wp,alph)
         endif
      endif
      call cubic (xk0, wp, alph, rad, qplus, qminus)

      if (rad .le. 0) then
         d2 = qplus - (xk0 + 1)
         if (d2 .gt. 0)  then
            eim = eim + ffq (qplus,ef,xk,wp,alph) - 
     1                  ffq (xk0+1,ef,xk,wp,alph)
         endif
         d3 = (xk0-1) - qminus
         if (d3 .gt. 0)  then
            eim = eim + ffq (xk0-1,ef,xk,wp,alph) - 
     1                  ffq (qminus,ef,xk,wp,alph)
c           beginning of the imaginary part and position of the cusp x0
            icusp = 1
         endif
      endif

      call quinn (xk0, rs, wp, ef, ei)
      if (eim .ge. ei)  eim = ei

      icount = icount+1
      return
      end
c     major revision, input now comes from main program feff
c     input data is passed here to indata for processing

      subroutine indata (iz, ihole, wsin, ionin)

      implicit double precision (a-h, o-z)
      save

c     logical unit from which to read input
      parameter (linp = 1)

      common /print/ iprint
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco

      integer*4 nstop,nuc
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /ps2/ dexv, dexe, dcop, test, teste,
     1             testy, testv, niter, ion, icut, iprat, irnorm

      character*40 ttl
      character*2  titre
      common /char2/ titre(30), ttl

      character*2  ttire(9)
      data ttire /'s ', 'p*', 'p ', 'd*', 'd ', 'f*', 'f ','g*', 'g '/

c following variables fixed as data by jm 4/20/87
      data i /0/
      data j /0/
      data k /0/
      data l /0/

      idep   = 0
      icut   = 0
c     Normal use, iprat = 1
      iprat  = 1
      irnorm = 1
      iex    = 1
      nuc    = 0

c idep=0 starting potential = thomas-fermi potential
c idep=1 starting potential read in from cards
c if icut is zero one corrects the potential by -(ion+1)/r
c if iprat is zero the pratt procedure is used
c if iex is zero one uses the unmodified slater exchange
c l=0 standard option for the bloc ofs points and their precision
c finite nuclear size option if nuc is positive
c if irnorm=1 renormalize potential to wigner-seitz radius

      dvc=137.0373d0
      dsal=dvc+dvc
      iz1=0
      ion1=0
      nuc1=-1
      dpas=0.05d0
      dr1=0.01d0
      nes=15

      niter=50

c     orig values:  teste 5.e-6, testy 1.e-5, testv 1.e-5, test 1.e-7
c     JM used teste 5.0e-5 to treat negative ion,
c     SZ changed teste to 1.0e-4 for selenium only to avoid convergence
c     problems with this particular atom.
c     teste set to 1.0e-4 to reduce run time (sz and jjr)
      teste = 1.0d-4
      testy=1.d-04
      testv=1.d-04
      test=1.d-07

      np=251
      nstop=30

c     Set dexv to zero for use with exafs model
      dexv = 0.0d0

      dexe=1.5d0
      dcop=0.3d0

c     i, j, k set to zero when old read statements removed
      i=0
      j=0
      k=0

c iz     = atomic number
c ion    = iz-number of electrons
c norb   = number of orbitals
c idep   = should be either 0 or 1
c i      = number of points for the integration = 251 by default
c j      = number of attempts to adjust the energy = 15 by default
c k      = number of iterations = 50 by default
c norbco = number of core orbitals

c put input data passed from feff into the necessary variables
      ws  = wsin
      ion = ionin
c     given iz, find norb, norbco, then den, nqn, nk and nel for
c     each orbital.
      call getorb (iz, ihole, ion, norb, norbco,
     1            den, nqn, nk, nel)

      if (norb .gt. nstop)  then
         if (iprint .ge. 5)  write(16,44) norb
         write(77,44) norb
   44    format (' norb=',i3,'too big')
         goto 999
      endif

c dexv = exchange coefficient for the potential: dexv=1. for slater
c dexe = exchange energy coefficient
c dexv should be equal to 2.*dexe/3. in order to satisfy the virial theo
c dexv=0.0 and iex=1, hedin-barth exchange and correlation is used

c dpas  = exponential step;  dr1 defines the first point = dr1/iz
c test  = energy precision criteria in dirac
c teste = self-consistency criteria for the energies of all the electron
c testy = self-consistency criteria for the wavefunctions
c testv = self-consistency criteria for the potential
      z=iz

      if (nuc .gt. 0)  then
         write(77,118)
  118    format(' enter atomic mass ')
         read (linp,*,end=900) dval
c        dval = atomic mass if nuc positive

         dval=z*(dval**(1.0d0/3.0d0))*2.267700d-05/exp(4.0d0*dpas)
         if (dval .le. dr1)  then
            dr1=dval
            nuc=5
         else
            dval=dval*exp(4.0d0*dpas)
            do 170 i=6,np
               d1=dr1*exp((i-1)*dpas)
               if (d1.ge.dval) goto 190
  170       continue
            write(77,180)
            if (iprint .ge. 5)  write(16,180)
  180       format (' error for the atomic mass')
            goto 999

  190       nuc=i
            dr1=dr1*dval/d1
         endif
      endif

      if (iprint .ge. 5)  write(16,210) ttl,niter,teste,testy,testv
  210 format (1h1,40x,A40,//,5x,'number of iterations',i4,//,
     1        5x,'precision of the energies',1pe9.2,//,
     2        23x,'wave functions  ',1pe9.2,//,
     3        23x,'potential',1pe9.2,/)

      xtmp = 8.8d0
      dr1=z*exp(-xtmp)

      if (iprint .ge. 5)  write(16,220) np,dr1,iz,dpas
  220 format (' the integration is made on ', i3,
     1        ' points-the first is equal to ' ,f7.4, '/', i2,/,
     2        ' and the step-size pas = ',f7.4,/)
      if (iprint .ge. 5)  write(16,230) test,nes,idep,icut,iprat
  230 format (' dans le sous programme resld la precision relative a',
     1        ' obtenir sur l energie est ', 1pe9.2,
     2        ' et le nombre d essais ',i3, //,
     3        'idep=', i3, 5x, 'icut=', i3, 5x, 'iprat=', i3, /)
      if (iprint .ge. 5)  write(16,240) dexv,dexe
  240 format ('  dexv=', 1pe14.7, '     dexe=' ,1pe14.7,
     1        ' if dexv=0.0 hedin-barth corr. and exchan. is used'/)
      k=0
      dval=z*z/(dvc*dvc)


      if (nuc.gt.0) then
         if (iprint .ge. 5)  write(16,250)
  250    format (1h0,30x,'finite nucleus case used'/)
      endif

      do 350 i=1,norb
c        den = orbital energy in atomic units and negative
c        nqn = principal quantum number; nk = kappa quantum number
c        nel = occupation of the orbital

         k=k+nel(i)
         if (den(i) .ge. 0.0)  den(i) = -z*z / (4.0*nqn(i)*nqn(i))

         nql(i)=iabs(nk(i))

         if (nk(i).lt.0) nql(i)=nql(i)-1
         if (nuc .le. 0)  then
            dfl(i)=nk(i)*nk(i)
            dfl(i)=sqrt(dfl(i)-dval)
         else
            dfl(i)=iabs(nk(i))
         endif
         l=2*iabs(nk(i))


         if (nql(i).lt.nqn(i)  .and.  nel(i).le.l  .and.
     1       nqn(i).gt.0       .and.  nql(i).le.4)   goto 340
            write(77,330) den(i),nqn(i),nql(i),j,nel(i)
            if (iprint .ge. 5)  write(16,330) den(i),nqn(i),nql(i),
     1                                         j,nel(i)
  330       format (' error in the card    ',e15.8,i2,3i2)
            goto 999
  340    continue
         j=nql(i)+iabs(nk(i))
         titre(i)=ttire(j)
         if (iprint .ge. 5)  write(16,345) nqn(i),titre(i),nel(i),
     1                                      den(i)
  345    format (7x,i1,a2,i16,1pe23.7)
  350 continue

      if (iprint .ge. 5)  write(16,370) norbco
  370 format (' no. of core orbitals=',i3)
      if (k.eq.(iz-ion)) goto 390
         write(77,380)
         if (iprint .ge. 5)  write(16,380)
  380    format (' error for the number of electrons')
         goto 999
  390 continue

      if (iprat .eq. 0)  then
         if (iprint .ge. 5)  write(16,410)
  410    format (1h0,'  the pratt procedure is used'/)
      else
         if (iprint .ge. 5)  write(16,430) ws
  430    format (1h0,'  wigner-seitz radius = ',0pf10.6,/)
      endif

      if (nuc .eq. nuc1)  then
         if (iz.eq.iz1.and.ion.eq.ion1) goto 600
         if (iz.eq.iz1) goto 470
      endif

c     dr(1)=dr1/z
c     do 460 i=2,np
c        dr(i)=dr(1)*exp((i-1)*dpas)
c 460 continue
c     Let's make this consistant with grid in other routines
c     dr array commeted out above
c     SIZ  December 1990
      do 461  i = 1, 251
         dr(i) = rr(i)
  461 continue

c starting potential

  470 val=-ion-1

c     Following code is a block, block ends at line 600
      if (idep .eq. 1)  then

c        read in starting potential (in a.u. and negative) if idep=1
         read (linp,480,end=900) (dv(i),i=1,np)
  480    format (8f9.4)

         if (iprint .ge. 5)  write(16,490) TTL,(dv(i),i=1,np)
  490    format (1h1, 40x, A40, //,
     1           5x, 'starting potential multiplied by r ' /,
     2           10(2x, f9.4))
         dval = -z/dv(1)
         if (nuc.gt.0)  dval = 1.0d0
         do 500 i=1,np
            dv(i)=dv(i)*dval/dr(i)
  500    continue

      else

         if (idep .ne. 0)  then
            write(77,510)
            if (iprint .ge. 5)  write(16,510)
  510       format (' error for idep')
            goto 999
         endif

         if (iz.ne.iz1  .or .  ion.le.ion1  .or.   nuc.ne.nuc1)  then
            do 520 i=1,np
               r=dr(i)
               dv(i)=fpot(r,z,val)
  520       continue
            if (nuc .gt. 0)  then
               do 530 i=1,nuc
                  dv(i) = dv(i) + z/dr(i) +
     1                    z*((dr(i)/dr(nuc))**2-3.0)/(dr(nuc)+dr(nuc))
  530          continue
            endif
            goto 600
         endif
      endif
      if (icut .eq. 0)  then
         do 540 i=1,np
            if ((dr(i)*dv(i)).gt.val)  dv(i)=val/dr(i)
  540    continue
      endif
      val=z+dv(1)*dr(1)
      if (nuc.gt.0)  val=z+dv(nuc)*dr(nuc)
      if (abs(val) .ge. 0.1d0)  then
         write(77,550)
         if (iprint .ge. 5)  write(16,550)
  550    format (' error for the potential ')
         goto 999
      endif

  600 continue
c     End of block above


      if (norb .ne. 1)  then
         do 620 i=2,norb
            k=i-1
            do 620 j=1,k
            if (nqn(i).eq.nqn(j)  .and. nk(i).eq.nk(j))   then
               write(77,610)
               if (iprint .ge. 5)  write(16,610)
  610          format (' standard configuration')
               goto 999
            endif
  620    continue
      endif

  630 iz1=iz
      ion1=ion
      nuc1=nuc
      do 660 i=1,norb
         nmax(i)=np
         l=1
         j=nqn(i)-nql(i)
         if ((j-2*(j/2)).eq.0) l=-l
         dq1(i)=l*nk(i)/iabs(nk(i))
         if (nuc .ne. 0  .and.  nk(i) .lt. 0)  then
            dq1(i)=dq1(i)*(nk(i)-dfl(i))*dvc/z
         endif
  660 continue


c  -- Normal return
      return


c  -- Error condition, stop program

c     Unexpected end of file during read -- stop program
  900 continue
      write(77,910)
  910 format (' Unexpected end of file')

c     Fatal error, stop gracefully (sic)
  999 continue
      stop 'INDATA-1'
      end
      subroutine inouh (dp,dq,dr,dq1,dfl,dv,z,test,nuc,nstop,jc)
c
c initial values for the outward integration
c dp=large component;     dq=small component;     dr=radial mesh
c dq1=slope at the origin of dp or dq;  dfl=power of the first term
c du=developpement limite;  dv=potential at the first point
c z=atomic number      test=test of the precision
c finite nuclear size if nuc is non-zero
c nstop controls the convergence  du developpement limite
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      integer*4 nstop,nuc
      common /ps1/ dep(5), deq(5), dd, dvc, dsal, dk, dm
c
c dep,deq=derivatives of dp and dq; dd=energy/dvc;
c dvc=speed of light in a.u.;
c dsal=2.*dvc   dk=kappa quantum number
c dm=exponential step/720.
c **********************************************************************
      common /trois/ dpno(4,30), dqno(4,30)
      dimension dp(251), dq(251), dr(251)
      do 10 i=1,10
      dp(i)=0.0
   10 dq(i)=0.0
      if (nuc) 20,20,60
   20 dval=z/dvc
      deva1=-dval
      deva2=dv/dvc+dval/dr(1)-dd
      deva3=0.0
      if (dk) 30,30,40
   30 dbe=(dk-dfl)/dval
      go to 50
   40 dbe=dval/(dk+dfl)
   50 dq(10)=dq1
      dp(10)=dbe*dq1
      go to 90
   60 dval=dv+z*(3.0d0-dr(1)*dr(1)/(dr(nuc)*dr(nuc)))/(dr(nuc)+dr(nuc))
      deva1=0.0d0
      deva2=(dval-3.0d0*z/(dr(nuc)+dr(nuc)))/dvc-dd
      deva3=z/(dr(nuc)*dr(nuc)*dr(nuc)*dsal)
      if (dk) 70,70,80
   70 dp(10)=dq1
      go to 90
   80 dq(10)=dq1
   90 do 100 i=1,5
      dp(i)=dp(10)
      dq(i)=dq(10)
      dep(i)=dp(i)*dfl
  100 deq(i)=dq(i)*dfl
      m=1
  110 dm=m+dfl
      dsum=dm*dm-dk*dk+deva1*deva1
      dqr=(dsal-deva2)*dq(m+9)-deva3*dq(m+7)
      dpr=deva2*dp(m+9)+deva3*dp(m+7)
      dval=((dm-dk)*dqr-deva1*dpr)/dsum
      dsum=((dm+dk)*dpr+deva1*dqr)/dsum
      j=-1
      do 130 i=1,5
      dpr=dr(i)**m
      dqr=dsum*dpr
      dpr=dval*dpr
      if (m.eq.1) go to 120
  120 dp(i)=dp(i)+dpr
      dq(i)=dq(i)+dqr
      if (abs(dpr/dp(i)).le.test.and.abs(dqr/dq(i)).le.test) j=1
      dep(i)=dep(i)+dpr*dm
  130 deq(i)=deq(i)+dqr*dm
      if (j.eq.1) go to 140
      dp(m+10)=dval
      dq(m+10)=dsum
      m=m+1
      if (m.le.20) go to 110
      nstop=45
  140 do 150 i=1,4
      dpno(i,jc)=dp(i+9)
  150 dqno(i,jc)=dq(i+9)
      return
      end
      subroutine inth (dp,dq,dv,dr)
c
c integration by the 5-point method of adams for the large
c component dp and the small component dq at the point dr;
c dv being the potential at this point
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      common /ps1/ dep(5), deq(5), db, dvc, dsal, dk, dm
c
c dep,deq the derivatives of dp and dq; db=energy/dvc;
c dvc=speed of light in atomic units; dsal=2.*dvc; dk=kappa quantum numb
c dm=exponential step/720.
c dkoef1=405./502., dkoef2=27./502.
c **********************************************************************
      data dkoef1 /0.9462151394422310d0/, dkoef2 /0.5378486055776890d-1/
      dpr=dp+dm*((251.0d0*dep(1)+2616.0d0*dep(3)
     > +1901.0d0*dep(5))-(1274.0d0
     1 *dep(2)+2774.0d0*dep(4)))
      dqr=dq+dm*((251.0d0*deq(1)+2616.0d0*deq(3)
     >   +1901.0d0*deq(5))-(1274.0d0
     1 *deq(2)+2774.0d0*deq(4)))
      do 10 i=2,5
      dep(i-1)=dep(i)
   10 deq(i-1)=deq(i)
      dsum=(db-dv/dvc)*dr
      dep(5)=-dk*dpr+(dsal*dr+dsum)*dqr
      deq(5)=dk*dqr-dsum*dpr
      dp=dp+dm*((106.0d0*dep(2)+646.0d0*dep(4)
     >   +251.0d0*dep(5))-(19.0d0*dep(1
     1 )+264.0d0*dep(3)))
      dq=dq+dm*((106.0d0*deq(2)+646.0d0*deq(4)
     >  +251.0d0*deq(5))-(19.0d0*deq(1
     1 )+264.0d0*deq(3)))
      dp=dkoef1*dp+dkoef2*dpr
      dq=dkoef1*dq+dkoef2*dqr
      dep(5)=-dk*dp+(dsal*dr+dsum)*dq
      deq(5)=dk*dq-dsum*dp
      return
      end
      subroutine intpol (a,b,fa,fb,fma,fmb,x,fx,fmx)
      implicit double precision (a-h,o-z)
c     Only output is fx, fmx
      complex*16 fa,fb,fma,fmb,fx,fmx
      dx=b-a
      d=(x-a)/dx
c     if (d*(1.0-d).lt.0.0) stop 'Died in intpol'
      if (d*(1.0d0-d).lt.0.0d0) then
         write(77,*) 'a, b, dx'
         write(77,*) a, b, dx
         write(77,*) 'x, x-a'
         write(77,*) x, x-a
         write(77,*) 'd, d*(1-d)'
         write(77,*) d, d*(1-d)
         stop 'Died in intpol'
      endif
      c2=3.0d0*(fb-fa)-(fmb+2.0d0*fma)*dx
      c3=2.0d0*(fa-fb)+(fma+fmb)*dx
      fx=fa+d*(dx*fma+d*(c2+d*c3))
      fmx=fma+d*(2.0d0*c2+3.0d0*c3*d)/dx
      return
      end
      subroutine ipack (iout, n, ipat)
      implicit double precision (a-h, o-z)

c     Input:  n          number of things to pack, nmax=8
c             ipat(1:n)  integers to pack
c     Output: iout(3)    packed version of n and ipat(1:n)
c
c     Packs n and ipat(1:n) into 3 integers, iout(1:3).  Algorithm
c     packs three integers (each between 0 and 1289 inclusive) into a
c     single integer.  Single integer must be INT*4 or larger, we assume
c     that one bit is wasted as a sign bit so largest positive int
c     is 2,147,483,647 = (2**31 - 1).
c     This version is specifically for the path finder and
c     degeneracy checker.

      dimension iout(3), ipat(n)
      dimension itmp(8)
      parameter (ifac1 = 1290, ifac2 = 1290**2)

      if (n .gt. 8)  stop 'ipack n too big'

      do 10  i = 1, n
         itmp(i) = ipat(i)
   10 continue
      do 20  i = n+1, 8
         itmp(i) = 0
   20 continue

      iout(1) = n       + itmp(1)*ifac1 + itmp(2)*ifac2
      iout(2) = itmp(3) + itmp(4)*ifac1 + itmp(5)*ifac2
      iout(3) = itmp(6) + itmp(7)*ifac1 + itmp(8)*ifac2

      return
      end
      subroutine upack (iout, n, ipat)
      implicit double precision (a-h, o-z)

c     retrieve n and ipat from iout
c     Input:  iout(3)  packed integers
c             n        max number to get, must be .le. 8
c     Output: n        number unpacked
c             ipat(1:n) unpacked integers

      dimension iout(3), ipat(n)
      dimension itmp(8)
      parameter (ifac1 = 1290, ifac2 = 1290**2)

      nmax = n
      if (nmax .gt. 8)  stop 'nmax .gt. 8 in upack'

      n = mod (iout(1), ifac1)
      if (n .gt. nmax)  stop 'nmax in upack too small'

      itmp(1) = mod (iout(1), ifac2) / ifac1
      itmp(2) = iout(1) / ifac2
      itmp(3) = mod (iout(2), ifac1)
      itmp(4) = mod (iout(2), ifac2) / ifac1
      itmp(5) = iout(2) / ifac2
      itmp(6) = mod (iout(3), ifac1)
      itmp(7) = mod (iout(3), ifac2) / ifac1
      itmp(8) = iout(3) / ifac2

      do 10  i = 1, n
         ipat(i) = itmp(i)
   10 continue

      return
      end
      subroutine istprm(nph, nat, iphat, rat, iatph, xnatph,
     1                   novr, iphovr, nnovr, rovr, folp, edens,
     2                   vclap, vtot, imt, inrm, rmt, rnrm, 
     2                   rhoint,
     3                   vint, rs, xf, xmu, rnrmav, intclc)

c     Finds interstitial parameters, rmt, vint, etc.
      implicit double precision (a-h, o-z)


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension folp(0:nphx)
      dimension edens(nrptx,0:nphx)
      dimension vclap(nrptx,0:nphx)
      dimension vtot (nrptx,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

c     intclc = 0, average evenly over all atoms
c              1, weight be lorentzian, 1 / (1 + 3*x**2), x = r/rnn,
c                 r   = distance to central atom,
c                 rnn = distance of near neighbor to central atom

c Find muffin tin radii.  We'll find rmt based on norman prescription,
c ie, rmt(i) = R * folp * rnrm(i) / (rnrm(i) + rnrm(j)),
c a simple average
c based on atoms i and j.  We average the rmt's from each pair of
c atoms, weighting by the volume of the lense shape formed by the
c overlap of the norman spheres.
c NB, if folp=1, muffin tins touch without overlap, folp>1 gives
c overlapping muffin tins.
c
c rnn is distance between sphere centers
c rnrm is the radius of the norman sphere
c xl_i is the distance to the plane containing the circle of the
c    intersection
c h_i  = rnrm_i - xl_i is the height of the ith atom's part of
c    the lense
c vol_i = (pi/3)*(h_i**2 * (3*rnrm_i - h_i))
c
c xl_i = (rnrm_i**2 - rnrm_j**2 + rnn**2) / (2*rnn)

      do 140  iph = 0, nph
         voltot = 0
         rmtavg = 0
         if (novr(iph) .gt. 0)  then
c           Overlap explicitly defined by overlap card
            do 124  iovr = 1, novr(iph)
               rnn  = rovr(iovr,iph)
               inph = iphovr(iovr,iph)
c              Don't avg if norman spheres don't overlap
               if (rnrm(iph)+rnrm(inph) .le. rnn)  goto 124
               voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
               voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
               rmttmp = rnn * folp(iph) * rnrm(iph) /
     1                  (rnrm(iph) + rnrm(inph))
               ntmp = nnovr(iovr,iph)
               rmtavg = rmtavg + rmttmp*voltmp*ntmp
               voltot = voltot + voltmp*ntmp
  124       continue
         else
            iat = iatph(iph)
            do 130  inat = 1, nat
               if (inat .eq. iat)  goto 130
               rnn = feff_dist(rat(1,inat), rat(1,iat))
               inph = iphat(inat)
c              Don't avg if norman spheres don't overlap
               if (rnrm(iph)+rnrm(inph) .lt. rnn)  goto 130
               voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
               voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
               rmttmp = rnn * folp(iph) * rnrm(iph) /
     1                  (rnrm(iph) + rnrm(inph))
               rmtavg = rmtavg + rmttmp*voltmp
               voltot = voltot + voltmp
  130       continue
         endif
         if (rmtavg .le. 0.0d0)  then
            write(77,132) iat, iph
  132       format (' WARNING: NO ATOMS CLOSE ENOUGH TO OVERLAP ATOM',
     1              i5, ',  UNIQUE POT', i5, '!!', /,
     2              '          Rmt set to Rnorman.  May be error in ',
     3              'input file.')
            rmt(iph) = rnrm(iph)
         else
            rmt(iph) = rmtavg / voltot
         endif
  140 continue

c     Need potential with ground state xc, put it into vtot
      do 160  iph = 0, nph
         call sidx (edens(1,iph), 250, rmt(iph), rnrm(iph),
     1              imax, imt(iph), inrm(iph))
         do 150  i = 1, imax
            rs = (edens(i,iph)/3)**(-third)
c           vhedbr from Von Barth Hedin paper, 1971
            vhedbr = -1.22177412d0/rs - 0.0504d0*log(30.0d0/rs + 1)
            vtot(i,iph) = vclap(i,iph) + vhedbr
  150    continue
  160 continue

c     What to do about interstitial values?
c     Calculate'em for all atoms, print'em out for all unique pots along
c     with derivative quantities, like fermi energy, etc.
c     Interstitial values will be average over all atoms in problem.

c     rnrmav is averge norman radius,
c     (4pi/3)rnrmav**3 = (sum((4pi/3)rnrm(i)**3)/n, sum over all atoms
c     in problem
      rnrmav = 0.0d0
      xn = 0.0d0
      rs = 0.0d0
      vint   = 0.0d0
      rhoint = 0.0d0
c     volint is total interstitial volume
      volint = 0

      do 170  iph = 0, nph
c        Use all atoms
         call istval(vtot(1,iph), edens(1,iph), rmt(iph), imt(iph),
     2                rnrm(iph), inrm(iph), vintx, rhintx, ierr)
c        if no contribution to interstitial region, skip this unique pot
         if (ierr .ne. 0)  goto 170
         call fermi (rhintx, vintx, xmu, rs, xf)
c        (factor 4pi/3 cancel in numerator and denom, so leave out)
         volx = (rnrm(iph)**3 - rmt(iph)**3)
         if (volx .le. 0)  goto 170
         volint = volint + volx * xnatph(iph)
         vint   = vint   + vintx * volx * xnatph(iph)
         rhoint = rhoint + rhintx* volx * xnatph(iph)
  170 continue
c     If no contribution to interstitial from any atom, die.
      if (volint .le. 0)  then
         write(77,*) ' No interstitial density.  Check input file.'
         stop 'ISTPRM'
      endif
      vint   = vint   / volint
      rhoint = rhoint / volint
      call fermi (rhoint, vint, xmu, rs, xf)
      do 180  iph = 0, nph
         rnrmav = rnrmav + xnatph(iph) * rnrm(iph)**3
         xn = xn + xnatph(iph)
  180 continue
      rnrmav = (rnrmav/xn) ** third


      return
      end

      double precision function calcvl(r1, r2, r)
      implicit double precision (a-h, o-z)

      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)

      xl = (r1**2 - r2**2 + r**2) / (2*r)
      h = r1 - xl
      calcvl = (pi/3) * h**2 * (3*r1 - h)
      return
      end
      subroutine istval (vtot, rholap, rmt, imt, rws, iws, vint, rhoint,
     1                   ierr)

c     This subroutine calculates interstitial values of v and rho
c     for an overlapped atom.  Inputs are everything except vint and
c     rhoint, which are returned.  vtot includes ground state xc.
c     rhoint is form density*4*pi, same as rholap
c
c     ierr = 0, normal exit
c          =-1, rmt=rws, no calculation possible

      implicit double precision (a-h, o-z)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      parameter (delta = 0.050000000000000d0)

      dimension vtot (nrptx)
      dimension rholap (nrptx)

c     Integrations are done in x (r = exp(x), see Louck's grid)
c     Trapezoidal rule, end caps use linear interpolation.
c     imt is grid point immediately below rmt, etc.
c     We will integrate over spherical shell and divide by volume of
c     shell, so leave out factor 4pi, vol = r**3/3, not 4pi*r**3/3,
c     similarly leave out 4pi in integration.

c     If rmt and rws are the same, cannot contribute to interstitial
c     stuff, set error flag
      vol = (rws**3 - rmt**3) / 3.0d0
      if (vol .le. 0.0d0)  then
         ierr = -1
         return
      endif
      ierr = 0

c     Calculation of vint including exchange correlation
c     Trapezoidal rule from imt+1 to iws
      vint = 0.0d0
      do 100  i = imt, iws-1
         fr = rr(i+1)**3 * vtot(i+1)
         fl = rr(i)**3   * vtot(i)
         vint = vint + (fr+fl)*delta/2.0d0
  100 continue
c     End cap at rws (rr(iws) to rws)
      xws = log (rws)
      xiws = xx(iws)
      g = xws - xiws
      fr = rr(iws+1)**3 * vtot(iws+1)
      fl = rr(iws)**3   * vtot(iws)
      vint = vint + (g/2.0d0) * ( (2.0d0-(g/delta))*fl + (g/delta)*fr)
c     End cap at rmt (rmt to rr(imt+1))
      xmt = log (rmt)
      ximt = xx(imt)
      g = xmt - ximt
      fr = rr(imt+1)**3 * vtot(imt+1)
      fl = rr(imt)**3   * vtot(imt)
      vint = vint - (g/2.0d0) * ( (2.0d0-(g/delta))*fl + (g/delta)*fr)
      vint = vint / vol

c     Calculation of rhoint
c     Trapezoidal rule from imt+1 to iws
      rhoint = 0
      do 200  i = imt, iws-1
         fr = rr(i+1)**3 * rholap(i+1)
         fl = rr(i)**3   * rholap(i)
         rhoint = rhoint + (fr+fl)*delta/2.0d0
  200 continue
c     End cap at rws (rr(iws) to rws)
      xws = log (rws)
      xiws = xx(iws)
      g = xws - xiws
      fr = rr(iws+1)**3 * rholap(iws+1)
      fl = rr(iws)**3   * rholap(iws)
      rhoint = rhoint + (g/2.0d0) 
     >       * ( (2.0d0-(g/delta))*fl + (g/delta)*fr)
c     End cap at rmt (rmt to rr(imt+1))
      xmt = log (rmt)
      ximt = xx(imt)
      g = xmt - ximt
      fr = rr(imt+1)**3 * rholap(imt+1)
      fl = rr(imt)**3   * rholap(imt)
      rhoint = rhoint - (g/2.0d0) 
     >   * ( (2.0d0-(g/delta))*fl + (g/delta)*fr)
      rhoint = rhoint / vol

      return
      end
      subroutine mcrith (npat, ipat, ri, indbet,
     1                   ipot, nncrit, fbetac, ckspc, xheap)
      implicit double precision (a-h, o-z)


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      dimension ipat(npatx)
      dimension ri(npatx+1), indbet(npatx+1)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)

c     Decide if we want the path added to the heap.

      if (ipat(npat) .eq. 0 .or. npat.le.2)  then
c        Partial path is used for xheap, not defined for ss and
c        triangles.  Special case: central atom added to end of path 
c        necessary for complete tree, but not a real path, again,
c        xheap not defined.  Return -1 as not-defined flag.
         xheap = -1
      else
c        Calculate xheap and see if we want to add path to heap.
c        Factor for comparison is sum over nncrit of
c        f(beta1)*f(beta2)*..*f(beta npat-2)/(rho1*rho2*..*rho npat-1).
c        Compare this to sum(1/p), multiply by 100 so we can think 
c        in percent.  Allow for degeneracy when setting crit.
         xheap = 0
         spinv = 0
         do 340  icrit = 1, nncrit
            x = ckspc(icrit) ** (-(npat-1)) * ri(npat-1)
            do 320  i = 1, npat-2
               ipot0 = ipot(ipat(i))
               x = x * fbetac(indbet(i),ipot0,icrit) / ri(i)
  320       continue
            spinv = spinv + 1/ckspc(icrit)
            xheap = xheap + x
  340    continue
         xheap = 100 * xheap / spinv

c        Factor for comparison is sum over nncrit of
c        New xheap:
c        Full chi is
c f(beta1)*f(beta2)*..*f(beta npat)cos(beta0)/(rho1*rho2*..*rho nleg).
c Some of this stuff may change when the path is modified --
c we can't use rho nleg or nleg-1, beta0, beta(npat) or beta(npat-1).
c We DO want to normalize wrt first ss path, f(pi)/(rho nn)**2.
c
c So save f(pi)/(rho nn)**2, 
c calculate 
c f(beta1)*f(beta2)*..*f(beta npat-2)/(rho1*rho2*..*rho npat-1).
c divide nn ss term by stuff we left out -- beta(npat), beta(npat-1),
c cos(beta0), rho nleg, rho nleg-1.
c
c Sum this over nncrit and try it out.
*
c        Sum over nncrit of
c        1/(rho1+rho2+..+rho npat-1).
*        reff = 0
*        do 350  i = 1, npat-1
*           reff = reff + ri(i)
* 350    continue
*        xss = 0
*        do 360  icrit = 1, nncrit
*           rho = ckspc(icrit) * reff
*           xss = xss + 1/rho
* 360    continue
*        xheap = 100 * xheap / xss
      endif

      return
      end
      subroutine mcritk (npat, ipat, ri, beta, indbet,
     1                   ipot, nncrit, fbetac, ckspc, xout, xcalcx)
      implicit double precision (a-h, o-z)


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      dimension ipat(npatx)
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)

c     xcalcx is max xcalc encountered so far.  Set to -1 to reset it --
c     otherwise it gets passed in and out as mcritk gets called.

c     We may want path in heap so that other paths built from this
c     path will be considered, but do not want this path to be
c     written out for itself.  Decide that now and save the flag
c     in the heap, so we won't have to re-calculate the mpprm
c     path parameters later.

c     Do not want it for output if last atom is central atom,
c     use xout = -1 as flag for undefined, don't keep it.
      if (ipat(npat) .eq. 0)  then
         xout = -1
         return
      endif

c     Make xout, output inportance factor.  This is sum over p of
c     (product of f(beta)/rho for the scatterers) * 
c                                 (cos(beta0)/rho(npat+1).
c     Compare this to xoutx, max xout encountered so far.
c     Multiply by 100 so we can think in percent.
      xcalc = 0
      do 460  icrit = 1, nncrit
         rho = ri(npat+1) * ckspc(icrit)
c        when beta(0)=90 degrees, get zero, so fudge with cos=.2
         x = max (abs(beta(npat+1)), 0.2d0) / rho
         do 420  iat = 1, npat
            rho = ri(iat) * ckspc(icrit)
            ipot0 = ipot(ipat(iat))
            x = x * fbetac(indbet(iat),ipot0,icrit) / rho
  420    continue
         xcalc = xcalc + x
  460 continue
      if (xcalc .gt. xcalcx)  xcalcx = xcalc
      xout = 100 * xcalc / xcalcx
      return
      end
      subroutine mkptz
c     makes polarization temsor ptz if necessary
      implicit double precision (a-h, o-z)

c     all input and output through common area /pol/

c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


c     addittonal local stuff to create polarization tensor ptz(i,j)
      real*8 e2(3)
      complex*16  e(3),eps,epc
      dimension eps(-1:1),epc(-1:1)


c     Begin to make polarization tensor
c     Normalize polarization vector
      x = sqrt(evec(1)**2 + evec(2)**2 + evec(3)**2)
      if (x .eq. 0.0d0) then
         write(77,*) 'STOP  Polarization vector of zero length'
         stop
      endif
      do 290  i = 1, 3
         evec(i) = evec(i) / x
  290 continue
      if (elpty .eq. 0.0d0) then
c        run linear polarization code
         do 291 i = 1, 3
            ivec(i) = 0.0d0
  291    continue
      endif
      x = sqrt (ivec(1)**2 + ivec(2)**2 + ivec(3)**2)
      if (x .gt. 0) then
c        run elliptical polarization code
         do 293  i = 1, 3
            ivec(i) = ivec(i) / x
  293    continue
         x = evec(1)*ivec(1)+evec(2)*ivec(2)+evec(3)*ivec(3)
         if (abs(x) .gt. 0.9d0) then
            write(77,*) 
     1         'STOP polarization almost parallel to the incidence'
            write(77,*) ' polarization',(evec(i), i=1,3)
            write(77,*) ' incidence   ',(ivec(i), i=1,3)
            write(77,*) ' dot product ', x
            stop
         endif
         if (x .ne. 0.0d0) then
c          if ivec not normal to evec then make in normal, keeping the
c          plane based on two vectors
           do 294 i = 1,3
              ivec(i) = ivec(i) - x*evec(i)
  294      continue
           x = sqrt (ivec(1)**2 + ivec(2)**2 + ivec(3)**2)
           do 295  i = 1, 3
              ivec(i) = ivec(i) / x
  295      continue
         endif
      else
         elpty = 0.0
      endif 
     
      e2(1) = ivec(2)*evec(3)-ivec(3)*evec(2)
      e2(2) = ivec(3)*evec(1)-ivec(1)*evec(3)
      e2(3) = ivec(1)*evec(2)-ivec(2)*evec(1)
      do 296  i = 1,3
        e(i) = (evec(i)+elpty*e2(i)*coni)
  296 continue 
      eps(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
      eps(0)  =   e(3)
      eps(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
      do 297  i = 1,3
        e(i) = (evec(i)-elpty*e2(i)*coni)
  297 continue 
      epc(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
      epc(0)  =   e(3)
      epc(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
      do 298 i = -1,1
      do 298 j = -1,1
c        ptz(i,j) = ((-1.0)**i)*epc(-i)*eps(j)/(1+elpty**2)
c       above - true polarization tensor for given ellipticity, 
c       below - average over left and right in order to have
c       path reversal simmetry
        ptz(i,j) = ((-1.0d0)**i)*(epc(-i)*eps(j)+eps(-i)*epc(j))
     1               /(1+elpty**2)/2.0d0
  298 continue
c     end of making polarization tensor

      return
      end
      subroutine mmtr(t3j,mmati)
c     calculates the part of matrix M which does not depend on energy
c     point.( see Rehr and Albers paper)

      implicit double precision (a-h, o-z)

c     all commons are inputs
c     inputs:
c        t3j: appropriate table of the 3j symbols
c     Inputs from common:
c        rotation matrix for ilegp
c        path data, eta(ilegp) and ipot(ilegp)
c        mtot,l0
c     Output:  mmati(...) 


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola


      save /rotmat/
      common /rotmat/ dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines


      complex*16 mmati
      dimension mmati(-mtot:mtot,-mtot:mtot),t3j(-mtot-1:mtot+1,-1:1)

      do 10 i = -mtot,mtot
      do 10 j = -mtot,mtot
         mmati(i,j)=0
  10  continue
      li = l0-1
c     l0 is final orb. momentum. Thus here we need to change code
c     in case when initial momemtum larger than final one.
      lx = min(mtot,l0)

      do 60 mu1 = -lx,lx
         mu1d = mu1+mtot+1
         do 50 mu2 = -lx,lx
            mu2d = mu2+mtot+1
            do 35  m0 = -li,li 
               do 34 i = -1,1
               do 34 j = -1,1
                  m1 = m0-j
                  m2 = m0-i
                  m1d = m1 + mtot+1
                  m2d = m2 + mtot+1
                  if (abs(m1).gt.lx .or. abs(m2).gt.lx)  goto 34
                  mmati(mu1,mu2) = mmati(mu1,mu2) + 
     1              dri(il0,mu1d,m1d,nsc+2)*dri(il0,m2d,mu2d,nleg)
     2              *exp(-coni*(eta(nsc+2)*m2+eta(0)*m1))
     3              *t3j(-m0,i)*t3j(-m0,j)*ptz(i,j)

c           dri(nsc+2)  is angle between z and leg1
c           dri(nsc+1)  is angle between last leg and z
c           eta(nsc+3)  is gamma between eps and rho1,
c           eta(nsc+2)  is alpha between last leg and eps
c           t3j(m0,i)    are 3j symbols multiplied by sqrt(3) 
   34          continue
   35       continue
            mmati(mu1,mu2) = mmati(mu1,mu2)*exp(-coni*eta(1)*mu1)
   50    continue
   60  continue

      return
      end
      subroutine mmtrxi (lam1x, mmati, ie, ileg, ilegp)
c     calculates matrix M in Rehr,Albers paper.
c     in polarization case
      implicit double precision (a-h, o-z)

c     all commons except for /fmat/ are inputs

c     inputs:
c       lam1x:  limits on lambda and lambda'
c       ie:  energy grid points
c       ileg, ilegp: leg and leg'
c
c     Inputs from common:
c        phases, use ph(ie,...,ilegp), and lmax(ie,ilegp)
c        lambda arrays
c        rotation matrix for ilegp
c        clmz for ileg and ilegp
c        path data, eta(ilegp) and ipot(ilegp)
c        xnlm array
c
c     Output:  fmati(...,ilegp) in common /fmatrx/ is set for
c              current energy point.

c     calculate scattering amplitude matrices
c     f(lam,lam') = sum_l tl gam(l,m,n)dri(l,m,m',ileg)gamt(l,m',n')
c                 *cexp(-i*m*eta),  eta = gamma+alpha'
c     lam lt lam1x, lam' lt lam2x such that m(lam) lt l0, n(lam) lt l0
c     gam = (-)**m c_l,n+m*xnlm, gamt = (2l+1)*c_ln/xnlm,
c     gamtl = gamt*tl


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      save /nlm/
      common /nlm/ xnlm(ltot+1,mtot+1)


      common /lambda/  
     4   mlam(lamtot), 	!mu for each lambda
     5   nlam(lamtot),	!nu for each lambda
     1   lamx, 		!max lambda in problem
     2   laml0x, 	!max lambda for vectors involving absorbing atom
     3   mmaxp1, nmax 	!max mu in problem + 1, max nu in problem


      save /clmz/
      complex*16 clmi
      common /clmz/ clmi(ltot+1,mtot+ntot+1,legtot)


c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola


      complex*16 fmati
      common /fmatrx/ fmati(lamtot,lamtot,legtot)


      save /rotmat/
      common /rotmat/ dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines


      complex*16 cam, camt, tltl,mmati
      dimension mmati(-mtot:mtot,-mtot:mtot)
      complex*16 gam(ltot+1,mtot+1,ntot+1),
     1           gamtl(ltot+1,mtot+1,ntot+1)

c     calculate factors gam and gamtl
      iln = il0
      ilx = il0
      do 30  il = iln, ilx
         tltl = 2*il - 1
         do 20  lam = 1, lam1x
            m = mlam(lam)
            if (m .lt. 0)  goto 20
            im = m+1
            if (im .gt. il)  goto 20
            in = nlam(lam) + 1
            imn = in + m
            if (lam .gt. lam1x)  goto 10
            cam = xnlm(il,im) * (-1)**m
            if (imn .le. il)  gam(il,im,in) = cam * clmi(il,imn,ileg)
            if (imn .gt. il)  gam(il,im,in) = 0
   10       if (lam .gt. lam1x) goto 20
            camt = tltl / xnlm(il,im)
            gamtl(il,im,in) = camt * clmi(il,in,ilegp)
   20    continue
   30 continue

      do 60 lam1 = 1,lam1x
         m1 = mlam(lam1)
         in1 = nlam(lam1) + 1
         iam1 = abs(m1) + 1
         do 50  lam2 = 1, lam1x
            m2 = mlam(lam2)
            in2 = nlam(lam2) + 1
            iam2 = iabs(m2) + 1
            imn1 = iam1 + in1 - 1
            fmati(lam1,lam2,ilegp) = mmati(m1,m2)*
     1                       gam(il0,iam1,in1)*gamtl(il0,iam2,in2)
   50    continue
   60 continue

      return
      end
      subroutine mpprmd (npat, ipat, ri, beta, eta)
      implicit double precision (a-h, o-z)
c     double precision version so angles come out right
c     for output...

c     Used with pathsd, a single precision code, so BE CAREFUL!!
c     No implicit, all variables declared explicitly.

c     make path parameters, ie, ri, beta, eta for each leg for a given
c     path.

c     Input is list of atoms (npat, ipat(npat)), output is
c     ri(npat+1), beta, eta.

      dimension ipat(npat)

      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     /atoms/ is single precision from pathsd
      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      complex*16  coni
      parameter (coni = (0,1))

      complex*16  alph(npatx+1), gamm(npatx+2), eieta
      double precision beta(npatx+1)
      double precision ri(npatx+1), eta(npatx+1)

      double precision x, y, z
      double precision ct, st, cp, sp, ctp, stp, cpp, spp
      double precision cppp, sppp

      n = npat + 1
      do 100  j = 1, n

c        get the atoms in this path
c        we actually have them already via the ipat array
c        remember that we'll want rat(,npat+1)=rat(,0) and
c                                 rat(,npat+2)=rat(,1) later on
c        make alpha, beta, and gamma for point i from 1 to N
c        NB: N is npat+1, since npat is number of bounces and N is
c            number of legs, or think of N=npat+1 as the central atom
c            that is the end of the path.
c
c        for euler angles at point i, need th and ph (theta and phi)
c        from rat(i+1)-rat(i)  and  thp and php
c        (theta prime and phi prime) from rat(i)-rat(i-1)
c
c        Actually, we need cos(th), sin(th), cos(phi), sin(phi) and
c        also for angles prime.  Call these  ct,  st,  cp,  sp   and
c                                            ctp, stp, cpp, spp.
c
c        We'll need angles from n-1 to n to 1,
c        so use rat(n+1) = rat(1), so we don't have to write code
c        later to handle these cases.

c        i = ipat(j)
c        ip1 = ipat(j+1)
c        im1 = ipat(j-1)
c        except for special cases...
         if (j .eq. n)  then
c           j central atom, j+1 first atom, j-1 last path atom
            i = 0
            ip1 = ipat(1)
            im1 = ipat(npat)
         elseif (j .eq. npat)  then
c           j last path atom, j+1 central, j-1 next-to last atom
c              unless only one atom, then j-1 central
            i = ipat(j)
            ip1 = 0
            if (npat .eq. 1)  then
               im1 = 0
            else
               im1 = ipat(npat-1)
            endif
         elseif (j .eq. 1)  then
c           j first atom, j+1 second unless only one,
c           then j+1 central, j-1 central
            i = ipat(j)
            if (npat .eq. 1)  then
               ip1 = 0
            else
               ip1 = ipat (j+1)
            endif
            im1 = 0
         else
            i = ipat(j)
            ip1 = ipat(j+1)
            im1 = ipat(j-1)
         endif

         x = rat(1,ip1) - rat(1,i)
         y = rat(2,ip1) - rat(2,i)
         z = rat(3,ip1) - rat(3,i)
         call strigd (x, y, z, ct, st, cp, sp)
         x = rat(1,i) - rat(1,im1)
         y = rat(2,i) - rat(2,im1)
         z = rat(3,i) - rat(3,im1)
         call strigd (x, y, z, ctp, stp, cpp, spp)

c        cppp = cos (phi prime - phi)
c        sppp = sin (phi prime - phi)
         cppp = cp*cpp + sp*spp
         sppp = spp*cp - cpp*sp

c        alph = exp**(i alpha)  in ref eqs 18
c        beta = cos(beta)
c        gamm = exp**(i gamma)
         alph(j) = st*ctp - ct*stp*cppp - coni*stp*sppp
         beta(j) = ct*ctp + st*stp*cppp
c        Watch out for roundoff errors
         if (beta(j) .lt. -1)  beta(j) = -1
         if (beta(j) .gt.  1)  beta(j) =  1
         gamm(j) = st*ctp*cppp - ct*stp + coni*st*sppp
         ri(j) = sdist (rat(1,i), rat(1,im1))
  100 continue

c     Make eta(i) = alpha(i) + gamma(i+1).  We only really need
c     exp(i*eta)=eieta, so that's what we'll calculate.
c     We'll need gamm(N+1)=gamm(npat+2)=gamm(1)
      gamm(npat+2) = gamm(1)
      do 150  j = 1, npat+1
         eieta = alph(j) * gamm(j+1)
         call sargd(eieta, eta(j))
  150 continue

c     Return beta as an angle, ie, acos(beta).  Check for beta >1 or
c     beta <1 (roundoff nasties)
      do 160  j = 1, npat+1
         if (beta(j) .gt.  1)  beta(j) =  1
         if (beta(j) .lt. -1)  beta(j) = -1
         beta(j) = dacos(beta(j))
  160 continue

      return
      end
      subroutine strigd (x, y, z, ct, st, cp, sp)
      implicit double precision (a-h, o-z)
      double precision x, y, z, ct, st, cp, sp, r, rxy
c     returns cos(theta), sin(theta), cos(phi), sin(ph) for (x,y,z)
c     convention - if x=y=0, phi=0, cp=1, sp=0
c                - if x=y=z=0, theta=0, ct=1, st=0
      parameter (eps = 1.0d-6)
      r = sqrt (x**2 + y**2 + z**2)
      rxy = sqrt (x**2 + y**2)
      if (r .lt. eps)  then
         ct = 1
         st = 0
      else
         ct = z/r
         st = rxy/r
      endif
      if (rxy .lt. eps)  then
         cp = 1
         sp = 0
      else
         cp = x / rxy
         sp = y / rxy
      endif

      return
      end
      subroutine sargd (c, th)
      implicit double precision (a-h, o-z)

      double precision x, y, th
      complex*16  c
      parameter (eps = 1.0d-6)
      x = dble(c)
      y = dimag(c)
      if (abs(x) .lt. eps)  x = 0
      if (abs(y) .lt. eps)  y = 0
      if (abs(x) .lt. eps  .and.  abs(y) .lt. eps)  then
         th = 0
      else
         th = atan2 (y, x)
      endif
      return
      end
      subroutine mpprmp (npat, ipat, xp, yp, zp)
      implicit double precision (a-h, o-z)

c     make path parameters,  xp, yp,zp for each atom for a given
c     path.

c     Input is list of atoms (npat, ipat(npat)), output are
c     x,y,z coord. of path in standard frame of reference
c     (see comments in timrep.f or here below)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola

      double precision  ro2, norm, zvec, xvec, yvec, ri, xp1, yp1, zp1
      dimension ipat(npatx+1), zvec(3), xvec(3), yvec(3)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      dimension xp(npatx), yp(npatx), zp(npatx)
      dimension xp1(npatx), yp1(npatx), zp1(npatx)
      dimension ri(3,npatx)

      parameter (eps4 = 1.0E-4)

c        get the atoms in this path
c        we actually have them already via the ipat array

c     initialize staff
      do 10 j = 1, npatx
         xp(j) = 0
         yp(j) = 0
         zp(j) = 0
         xp1(j) = 0
         yp1(j) = 0
         zp1(j) = 0
   10 continue
      nleg = npat + 1
      do 20  j = 1, npat
      do 20  i = 1, 3
         ri(i,j) = rat(i,ipat(j)) - rat(i,0)
   20 continue
      do 30  j = nleg, npatx
      do 30  i = 1, 3
         ri(i,j) = 0
   30 continue
      do 40 i =1, 3
         xvec(i) = 0.0
         yvec(i) = 0.0
         zvec(i) = 0.0
   40 continue

      if (.not. pola) then
c        z-axis along first leg
         norm = ri(1,1)*ri(1,1)+ri(2,1)*ri(2,1)+ri(3,1)*ri(3,1)
         norm = sqrt(norm)
         do 140 i = 1, 3
           zvec(i) = ri(i,1)/norm
  140    continue
      else
c        z-axis in direction of polarization
         do 120 i = 1, 3
           zvec(i) = evec(i)
  120    continue
      endif

      do 160 j = 1,npat
      do 160 i = 1, 3
        zp1(j) = zp1(j) + zvec(i)*ri(i,j)
  160 continue

      num = 1
      if (.not. pola) then
c        first nonzero z-coord. is already positive
         goto 240
      endif
  200 continue
      if (abs(zp1(num)) .gt. eps4) then
         if (zp1(num) .lt. 0.0) then
c           inverse all z-coordinates and zvec, if 
c           first nonzero z-coordinate is negative 
            do 210 j = 1, 3
               zvec(j) = - zvec(j)
  210       continue
            do 220 j = 1, npat
               zp1(j) = - zp1(j)
  220       continue
         endif
         goto 240
      endif
      num = num +1
      if (num .lt. nleg) then
         goto 200
      endif
c     here first nonzero z-coordinate is positive
  240 continue

      num = 1
  300 continue
      ro2 = 0.0
      do 310 i =1, 3
         ro2 = ro2 + ri(i,num)*ri(i,num)
  310 continue
c     looking for first atom which is not on z-axis
      ro2 = ro2 - zp1(num)*zp1(num)
      ro2 = sqrt(abs(ro2))
      if (ro2 .ge. eps4) then
c     if atom not on the z-axis then
         if (elpty .eq. 0.0) then
c           if not elliptical polarization then
c           choose x-axis so that x-coord. positive and y=0.
            do 320 i = 1, 3
               xvec(i) = ri(i,num) - zvec(i)*zp1(num)
  320       continue
            do 330 i = 1, 3
               xvec(i) = xvec(i)/ro2
  330       continue
         else
c           if elliptical polarization then
c           choose x-axis along incident beam
            do 350 i =1, 3
               xvec(i) = ivec(i)
  350       continue
         endif
         yvec(1) = zvec(2)*xvec(3) - zvec(3)*xvec(2)
         yvec(2) = zvec(3)*xvec(1) - zvec(1)*xvec(3)
         yvec(3) = zvec(1)*xvec(2) - zvec(2)*xvec(1)
         goto 390
      endif
      num = num + 1
      if (num .lt. nleg) then
         goto 300
      endif
  390 continue

c     calculate x,y coord for each atom in chosen frame of reference
      do 400 j = 1, npat
      do 400 i =1,3
         xp1(j) = xp1(j) + xvec(i)*ri(i,j)
         yp1(j) = yp1(j) + yvec(i)*ri(i,j)
  400 continue

      if ( elpty .ne. 0.0) then
c        if no polarization or linear polarization then first nonzero
c        x-coordinate is already positive, no need to check it.
         num = 1
  500    continue
         if (abs(xp1(num)) .ge. eps4) then
            if (xp1(num) .lt. 0.0) then
               do 510 j = 1, npat
                  xp1(j) = - xp1(j)
  510          continue
            endif
            goto 520
         endif
         num = num + 1
         if (num .lt. nleg) then
            goto 500
         endif
  520    continue
      endif

      num = 1
  570 continue
c     inverse all y-coordinates if first nonzero y-coord is negative
      if (abs(yp1(num)) .ge. eps4) then
         if (yp1(num) .lt. 0.0) then
            do 580 j = 1, npat
               yp1(j) = - yp1(j)
  580       continue
         endif
         goto 590
      endif
      num = num + 1
      if (num .lt. nleg) then
         goto 570
      endif
  590 continue

      do 595 j = 1, npat
        xp(j) = xp1(j)
        yp(j) = yp1(j)
        zp(j) = zp1(j)
  595 continue
c     now xp,yp,zp represent the path in standard order
      return
      end
      subroutine mrb (npat, ipat, ri, beta)
      implicit double precision (a-h, o-z)

c     Make ri, beta and rpath path parameters for crit calculations.

c     Input is list of atoms (npat, ipat(npat)), output is
c     ri(npat+1), beta, eta.


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      dimension ipat(npatx)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      dimension beta(npatx+1), ri(npatx+1), ipat0(npatx+1)

      nleg = npat+1
c     central atom is atom 0 in rat array
c     need local ipat0 array since we use ipat0(npat+1), final atom
c     in path (final atom is, of course, the central atom)
      do 10  i = 1, npat
         ipat0(i) = ipat(i)
   10 continue
      ipat0(nleg) = 0

      do 30  ileg = 1, nleg
c        make beta and ri for point i from 1 to N
c        NB: N is npat+1, since npat is number of bounces and N is
c            number of legs, or think of N=npat+1 as the central atom
c            that is the end of the path.
c
c        We'll need angles from n-1 to n to 1,
c        so use rat(n+1) = rat(1), so we don't have to write code
c        later to handle these cases.

c        Work with atom j
c        jp1 = (j+1)
c        jm1 = (j-1)
         j = ileg
         jm1 = j-1
         jp1 = j+1
c        Fix special cases (wrap around when j is near central atom,
c        also handle ss and triangular cases).
         if (jm1 .le.    0)  jm1 = nleg
         if (jp1 .gt. nleg)  jp1 = 1

         jat = ipat0(j)
         jm1at = ipat0(jm1)
         jp1at = ipat0(jp1)

         ri(ileg) = sdist (rat(1,jat), rat(1,jm1at))

c        Make cos(beta) from dot product
         call dotcos(rat(1,jm1at), rat(1,jat), rat(1,jp1at),
     1               beta(ileg))
   30 continue

      rpath = 0
      do 60  ileg = 1, nleg
         rpath = rpath + ri(ileg)
   60 continue

      return
      end
      subroutine dotcos (rm1, r, rp1, cosb)
      implicit double precision (a-h, o-z)
      dimension rm1(3), r(3), rp1(3)

      parameter (eps = 1.0d-8)

      cosb = 0
      do 100  i = 1, 3
         cosb = cosb + (r(i)-rm1(i)) * (rp1(i)-r(i))
  100 continue

c     if denom is zero (and it may be if 2 atoms are in the same place,
c     which will happen when last path atom is central atom), set
c     cosb = 0, so it won't be undefined.

      denom = (sdist(r,rm1) * sdist(rp1,r))
      if (denom .gt. eps)  then
         cosb = cosb / denom
      else
         cosb = 0
      endif
      return
      end
      subroutine outcrt (npat, ipat, ckspc,
     1    nncrit, fbetac, ne, ik0, cksp, fbeta, ipotnn, ipot,
     1    xport, xheap, xheapr,
     1    xout, xcalcx)
      implicit double precision (a-h, o-z)

c     This make pw importance factor for pathsd, also recalculates
c     pathfinder criteria for output.  Pathfinder recalculation
c     is hacked from ccrit, so be sure to update this if ccrit
c     is changed.


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      dimension ipat(npatx)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)
      dimension fbeta(-nbeta:nbeta,0:npotx,nex), cksp(nex)

c     local variables
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)
      dimension xporti(nex)
      parameter (eps = 1.0d-6)

c     Space for variables for time reversed path (used in xheapr
c     calculation below)
      dimension ipat0(npatx)
      dimension ri0(npatx+1), indbe0(npatx+1)

c     mrb is 'efficient' way to get only ri and beta
c     note that beta is cos(beta)
      call mrb (npat, ipat, ri, beta)

c     Make index into fbeta array (this is nearest cos(beta) grid point,
c     code is a bit cute [sorry!], see prcrit for grid).
      do 290  i = 1, npat+1
         tmp = abs(beta(i))
         n = tmp / 0.025d0
         del = tmp - n*0.025d0
         if (del .gt. 0.0125d0)  n = n+1
         if (beta(i) .lt. 0)  n = -n
         indbet(i) = n
  290 continue

c     Make pw importance factor by integrating over all points
c     above the edge
c     Path importance factor is integral d|p| of
c        (product of f(beta)/rho for the scatterers) * cos(beta0)/rho0
      do 560  ie = ik0, ne
         rho = ri(npat+1) * cksp(ie)
         crit = max (abs(beta(npat+1)), 0.2d0) / rho
         do 520  iat = 1, npat
            rho = ri(iat) * cksp(ie)
            ipot0 = ipot(ipat(iat))
            crit = crit * fbeta(indbet(iat),ipot0,ie) / rho
  520    continue
         xporti(ie) =  abs(crit)
  560 continue
c     integrate from ik0 to ne
      nmax = ne - ik0 + 1
      call strap (cksp(ik0), xporti(ik0), nmax, xport)

c     Stuff for  output.
c     Heap crit thing (see ccrit and mcrith for comments)
c     If a path got time reversed, its xheap may be smaller than
c     it was before it got time-reversed.  So calculate it both
c     ways.
c     xheap for path, xheapr for time-reversed path

      xheap  = -1
      xheapr = -1
      call mcrith (npat, ipat, ri, indbet,
     1             ipot, nncrit, fbetac, ckspc, xheap)

c     Prepare arrays for time reversed path and make xheapr
c     See timrev.f for details on indexing here.

      nleg = npat+1
c     ri
      do 200  i = 1, nleg
         ri0(i) = ri(nleg+1-i)
  200 continue
c     indbet  and ipat
      indbe0(nleg) = indbet(nleg)
      do 210  i = 1, nleg-1
         indbe0(i) = indbet(nleg-i)
         ipat0(i) = ipat(nleg-i)
  210 continue

      call mcrith(npat, ipat0, ri0, indbe0,
     1             ipot, nncrit, fbetac, ckspc, xheapr)

c     Keep crit thing (see mcritk for comments)
      call mcritk (npat, ipat, ri, beta, indbet,
     1             ipot, nncrit, fbetac, ckspc, xout, xcalcx)
c     print*, npat, xout, xcalcx

      return
      end
      subroutine ovrlp (iph, iphat, rat, iatph, ifrph, novr,
     1                  iphovr, nnovr, rovr, iz, nat, rho, vcoul,
     2                  edens, vclap, rnrm)

c     Overlaps coulomb potentials and electron densities for current
c     unique potential
      implicit double precision (a-h, o-z)


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension ifrph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension iz(0:nfrx)
      dimension rho(251,0:nfrx)
      dimension vcoul(251,0:nfrx)
      dimension edens(nrptx,0:nphx)
      dimension vclap(nrptx,0:nphx)
      dimension rnrm(0:nphx)

c     find out which free atom we're dealing with
      ifr = ifrph(iph)

c     start with free atom values for current atom
      do 100  i = 1, 250
         vclap(i,iph) = vcoul(i,ifr)
         edens(i,iph) = rho  (i,ifr)
  100 continue

      if (novr(iph) .gt. 0)  then
         do 104  iovr = 1, novr(iph)
            rnn  = rovr(iovr,iph)
            ann  = nnovr(iovr,iph)
            infr = ifrph(iphovr(iovr,iph))
            call sumax (250, rnn, ann, vcoul(1,infr), vclap(1,iph))
            call sumax (250, rnn, ann, rho  (1,infr), edens(1,iph))
  104    continue
      else
c        Do overlapping from geometry with model atom iat
         iat = iatph(iph)

c        overlap with all atoms within r overlap max (rlapx)
c        12 au = 6.35 ang  This number pulled out of a hat...
         rlapx = 12
c        inat is Index of Neighboring ATom
         do 110  inat = 1, nat
c           don't overlap atom with itself
            if (inat .eq. iat)  goto 110

c           if neighbor is too far away, don't overlap it
            rnn = feff_dist(rat(1,inat), rat(1,iat))
            if (rnn .gt. rlapx)  goto 110

            infr = ifrph(iphat(inat))
            call sumax (250, rnn, one, vcoul(1,infr), vclap(1,iph))
            call sumax (250, rnn, one, rho  (1,infr), edens(1,iph))
  110       continue
      endif

c     set norman radius
      call frnrm (edens(1,iph), iz(ifr), rnrm(iph))

      return
      end
      subroutine paths(ckspc, fbetac, pcritk, pcrith, nncrit,
     1                  rmax, nlegxx, ipotnn)
 
      implicit double precision (a-h, o-z)

c     finds multiple scattering paths
c     This is single precision, units are Angstroms.  BE CAREFUL!

c     pcrith is cut-off fraction used when building paths
c            (path criterion for heap)
c     pcritk is cut-off fraction used on output
c            (path criterion for keeping)

c     ipotnn is output, used by pathsd to duplicate paths criteria,
c     which are used only for diagnostic output.


      character*72 header
      common /header_common/ header


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)


      character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
      common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch


c     This common in pathsd, mpprm
      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      dimension m(-1:natx,0:natx)
      dimension mindex(natx+1)
c     Used for packed integers
      dimension iout(3)

c     ok true if all paths to rmax found.  If heap full, npx exceeded,
c     etc., last general shell may be incomplete, set ok=.false.
      logical ok

c     Heap data structure:
c     index is the pointer to the element of the data structure.
c     Each element contains
c        r        total path length
c                 Note that r is sorted along with index -- this keeps
c                 the heap maintenance routines fast.
c        mi, mj   m matrix elements used to place last atom in this path
c        npat     number of atoms in this path
c        ipat(npatx) indices of atoms in this path
c     next is the index of the next data structure element available.
c     If an element is freed, npat is the index of the free element
c     to use after using current next element.
c     nx is max number in heap
      integer    nx
      parameter (nx = 10000)
c     parameter (nx = 60 000)
c     r also used in making m matrix, must have nx >= natx+1
      integer   index(nx), npx, np, n, ip, i
c     parameter (npx = 100 000)
      parameter (npx = 4000000)
      dimension r(nx), mi(nx), mj(nx)
      dimension npat(nx)
      dimension ipat (npatx,nx)
c     Keep this path on output
      logical keep1(nx), kp1tmp

c     Used with ipack, so need ipat(8)
      dimension ipat0(8)

c     paths are typically about 10 or 20 Ang
      parameter (big = 1.0d3)

      parameter (nheadx = 30)
      character*80  head(nheadx)
      character*80  title
      dimension lhead(nheadx)

c     Returned from criterion checker, false if path fails criterion
      logical keep

c     read input
c     header...
c     i, x, y, z, ipot, i1b   of nat+1 atoms (i=0 is central atom)
      open (1, file=trim(header)//'geom.dat', status='old', iostat=ios)
      call chopen (ios, trim(header)//'geom.dat', 'paths')
      nhead = nheadx
      call rdhead (1, nhead, head, lhead)
c     header from geom.dat includes carriage control...
c     nlegxx is max number of legs user wants to consider.
c     nlegs = npat+1, so set npatxx = min (npatx, nlegxx-1)
      npatxx = min (npatx, nlegxx-1)
c     Input rmax is one-way distances
      rmax = rmax*2
      nat = -1
c     ratx is distance to most distant atom, used to check rmax
      ratx = 0
   10 continue
         nat = nat+1
         if (nat .gt. natx)  then
            write(77,*) ' nat, natx ', nat, natx
            stop 'Bad input'
         endif
         read(1,*,end=20)  idum, (rat(j,nat),j=1,3), ipot(nat), i1b(nat)
         rtmp = sdist(rat(1,nat),rat(1,0))
         if (rtmp .gt. ratx)  ratx = rtmp
      goto 10
   20 continue
      nat = nat-1
      close (unit=1)

c     Warn user if rmax > dist to most distant atom
      if (rmax/2.0d0 .gt. ratx+0.02d0)  then
        write(77,*) '   WARNING:  rmax > distance to most distant atom.'
        write(77,*) '             Some paths may be missing.'
        write(77,*) '             rmax, ratx ', rmax/2, ratx
      endif

c     Count number of 1st bounce atoms (at least 1 required).
      n1b = 0
      do 30  i = 1, nat
         if (i1b(i) .gt. 0)  n1b = n1b + 1
   30 continue
      if (n1b .lt. 1) stop 'At least one 1st bounce atoms required.'

      if (rmax .ge. big)  stop 'Hey, get real with rmax!'

c     Make title for this run, include carriage control because head
c     (read above) includes carriage control.
      write(title,32)  rmax/2, pcritk, pcrith, vfeff, vpaths
   32 format(' Rmax', f8.4, ',  keep limit', f7.3,
     1       ', heap limit', f7.3, t57, 2a12)

      write(77,34) rmax/2, pcritk, pcrith
   34 format ('    Rmax', f8.4,
     1        '  keep and heap limits', 2f12.7)

      write(77,36) '   Preparing neighbor table'
   36 format (1x, a)
c     prepare table telling distance from atom i to atom j and then
c     back to central atom
c     First bounce is m(-1,...), m(0,...) is bounces from central
c     atom that are not first bounces.
      do 60  i = -1, nat
         ir = i
         if (i .eq. -1)  ir = 0
         do 40  j = 0, nat
c           r begins with element 1 so sort routine later will work
            r(j+1) = sdist (rat(1,ir), rat(1,j))
            r(j+1) = r(j+1) + sdist (rat(1,j), rat(1,0))
c           we don't need m(i,i), since this will be = shortest
c           of the r(j), so just set it to something very big,
c           it will sort to the end of this row and it won't
c           bother us
            if (j .eq. ir)  r(j+1) = big
c           If we're doing first bounce, use only the allowed first
c           bounce paths.
            if (i .eq. -1)  then
               if (i1b(j) .le. 0)  r(j+1) = big
            endif
   40    continue

c        prepare row i of m table
c        m is a distance table ordered such that distance from
c               i to m(i,0) to 0 <
c               i to m(i,1) to 0 <
c               i    m(i,2)    0 <
c               :    :    :
c               i    m(i,nat)  0
c
c        That is, m(i,0) is index of atom that gives shortest path,
c                 m(i,1)                        next shortest path, etc.
c        Note that m(0,0) is shortest single bounce path.

c        Again, r and mindex go from 1 to nat+1, m goes from 0 to nat
         call sortir (nat+1, mindex, r)
         do 50  j = 0, nat
            m(i,j) = mindex(j+1)-1
   50    continue
   60 continue

      write(77,61)
   61 format ('    nfound  nheap  nheapx  nsc    r')

c     initialize heap data space next pointers
      do 70  i = 1, nx-1
         npat(i) = i+1
   70 continue
      npat(nx) = -1
c     initial condition:  make the first path
c     n    number in heap
c     nna  number skipped counter
c     nhx  number used in heap max, a counter
      n = 1
      nna = 0
      nhx = n
      nwrote = 0
      index(n) = 1
      ip = index(n)
      next = 2
      mi(ip) = -1
      mj(ip) = 0
      npat(ip) = 1
      ipat(npat(ip),1) = m(mi(ip),mj(ip))

c     near neighbor is atom ipat(npat(ip),1) for first path into heap
      ipotnn = ipot(ipat(npat(ip),1))

c     Someday change keep and keep1 to lkeep and lheap to match
c     ccrit variable names.
c     Initialize keep criterion
      xcalcx = -1
      call ccrit(npat(ip), ipat(1,ip), ckspc,
     1    fbetac, rmax, pcrith, pcritk, nncrit, ipotnn, ipot,
     2    r(n), keep, keep1(ip), xcalcx)

      open (file=trim(header)//'paths.bin', unit=3, access='sequential',
     1      form='unformatted', status='unknown', iostat=ios)
      call chopen (ios, trim(header)//'paths.bin', 'paths')
c     These strings are all char*80 and include carriage control
      write(3) nhead+1
      do 88  ihead = 1, nhead
         write(3) head(ihead)
         write(3) lhead(ihead)
   88 continue
      write(3) title
      write(3) istrln(title)
      write(3)  nat
      do 90  i = 0, nat
         write(3) (rat(j,i),j=1,3), ipot(i), i1b(i)
   90 continue

c     r is the heap, index is the pointer to the rest of the data
c     np is the number of paths found and saved
      np = 0
c     nbx  mpat max (Number of Bounces maX)
      nbx = 0

c     done if path at top of heap is longer than longest path we're
c        interested in
c     done if max number of paths we want have been found
c     begin 'while not done' loop
      ok = .false.
  800 continue
         if (r(1) .gt. rmax  .or.  np .ge. npx .or. n.le.0)  then
c           n=0 means heap is empty
            if (n.le.0)  ok=.true.
c           if (n.le.0)  print*, '   Heap empty'
            goto 2000
         endif

c        save element at top of heap in arrays labeled 0
c        dump to unit 3 (unformatted)
         ip = index(1)
         npat0 = npat(ip)
         do 100  i = 1, npat0
            ipat0(i) = ipat(i,ip)
  100    continue
         r0 = r(1)

c        Don't write out path if last atom is central atom, or
c        if it doesn't meet pcritk
         if (ipat0(npat0).eq.0 .and. keep1(ip)) then
            write(77,*) ipat0(npat0), keep1(ip), ' odd case...'
         endif
         if (ipat0(npat0).ne.0 .and. keep1(ip))  then
            np = np+1
c           pack integers
            call ipack (iout, npat0, ipat0)
            write(3)  r0, iout
            nwrote = nwrote+1
c           write status report to screen
            if (mod(np,1000) .eq. 0)  then
               write(77,132) np, n, nhx, nbx, r0/2
  132          format (4x, i6, i7, i8, i4, f10.4)
            endif
         endif

         if (np .ge. npx)  then
            write(77,*) np, ' paths found.  (np .ge. npx)'
            goto 2000
         endif

c        Make new path by replacing last atom in path from top of heap,
c        put this path on top of heap and buble it down.  If row is
c        finished, or new path is too long, don't add it, instead
c        move last path in heap to the top.
c        If working on row mi=-1 (first bounce atoms), don't
c        use them if not allowed 1st bounce atoms.
         mj(ip) = mj(ip) + 1
         if (mi(ip).eq.-1  .and.  i1b(m(mi(ip),mj(ip))).le.0)  then
c           not allowed first bounce atom
            r(1) = big
            keep = .false.
c           print*, '1st bounce limit!'
         elseif (mj(ip) .ge. nat)  then
c           we've finished a row of m matrix
            r(1) = big
            keep = .false.
         else
c           new path has same indices, etc.  Only need to replace
c           last atom.
            ipat(npat(ip),ip) = m(mi(ip),mj(ip))
            call ccrit (npat(ip), ipat(1,ip), ckspc,
     1                  fbetac, rmax, pcrith, pcritk, nncrit,
     1                  ipotnn, ipot,
     2                  r(1), keep, keep1(ip), xcalcx)
         endif

c        If r is bigger than rmax or keep=false, remove element from
c        heap by taking the last element in the heap and moving it to
c        the top.  Then bubble it down.  When removing an element
c        from the heap, be sure to save the newly freed up index.
c        r(1) and index(1) are new path, set above
         if (r(1).gt.rmax .and. keep)  then
            write(77,*) 'odd case rmax...'
         endif
         if (r(1).gt.rmax .or. .not.keep)  then
            index(1) = index(n)
            r(1) = r(n)
c           use npat as pointer to next free location
            npat(ip) = next
            next = ip
            n = n-1
c           nna is Number Not Added to heap
            nna = nna + 1
c           Maybe heap may be empty here, but that's alright
         endif
         if (npat(index(1)).gt.nbx .and. n.gt.0)  nbx = npat(index(1))

c        If heap is empty, don't call hdown.
         if (n.gt.0)  call hdown (r, index, n)

c        and make a new path by adding an atom onto the end of the path
c        we saved, put this at the end of the heap and bubble it up.
c        Do this only if it won't be too many bounces.
         if (npat0+1 .le. npatxx)  then
            ip = next
            if (ip .lt. 0)  then
c              print*, '   Heap full'
               goto 2000
            endif
            next0 = npat(ip)
            do 200  i = 1, npat0
               ipat(i,ip) = ipat0(i)
  200       continue
            mi(ip) = ipat0(npat0)
            mj(ip) = 0
            npat(ip) = npat0+1
            ipat(npat(ip),ip) = m(mi(ip),mj(ip))
            call ccrit (npat(ip), ipat(1,ip), ckspc,
     1                  fbetac, rmax, pcrith, pcritk, nncrit,
     1                  ipotnn, ipot,
     2                  rtmp, keep, kp1tmp, xcalcx)
            if (rtmp .gt. rmax  .and.  keep)  then
               write(77,*) 'odd case rmax and tmp...'
            endif
            if (rtmp .gt. rmax  .or.  .not.keep)  then
               npat(ip) = next0
               nna = nna+1
            else
c              add it to the heap
               next = next0
               n = n+1
               if (n .gt. nhx)  nhx = n
               index(n) = ip
               r(n) = rtmp
               keep1(ip) = kp1tmp
               if (npat(index(n)) .gt. nbx)  nbx = npat(index(n))
               call hup (r, index, n)
            endif
         endif

      goto 800
 2000 continue
c     end of 'while not done' loop
      if (.not. ok)  then
         write(77,*) '   Internal path finder limit exceeded -- ',
     1           'path list may be incomplete.'
      endif
      close (unit=3)
      write(77,2010) np, nhx, nbx
 2010 format ('    Paths found', i9, 3x,
     1        '(nheapx, nbx', i8, i4, ')')

      end
      subroutine pathsd(ckspc, fbetac, ne, ik0, cksp, fbeta,
     1                   critpw, ipotnn, ipr2, 
     1                   pcritk, pcrith, nncrit, potlbl)

      implicit double precision (a-h, o-z)
c     New degeneracy checker, cute and hopefully fast for large
c     problems

c     pcritk and pcrith used only for analysis after outcrt

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

c     np1x  number of paths to consider at 1 time
      parameter (np1x = 12 000)
c     parameter (np1x = 60 000)
      dimension iout(3,np1x), iout0(3)

      dimension index(np1x)
      double precision dhash(np1x), dcurr, ddum
      dimension rx(npatx), ry(npatx), rz(npatx), ipat(npatx+1)
      dimension rx0(npatx), ry0(npatx), rz0(npatx), ipat0(npatx+1)
      double precision rid(npatx+1), betad(npatx+1), etad(npatx+1)

      parameter (nheadx = 40)
      character*80 head(nheadx)
      dimension lhead(nheadx)

      character*6  potlbl(0:npotx)

c     eps5 for rtotal range, eps3 for individual leg parameters.
c     eps3 large since code single precision and don't want round-off
c     error to reduce degeneracy.
      parameter (eps5 = 2.0d-5)
      parameter (eps3 = 2.0d-3)

      logical ldiff, last
      parameter (necrit=9, nbeta=40)
      real*8 fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)
      real*8 fbeta(-nbeta:nbeta,0:npotx,nex), cksp(nex)

      write(77,30) critpw
   30 format ('    Plane wave chi amplitude filter', f7.2, '%')

c     Read atoms info
      open (file=trim(header)//'paths.bin', unit=3, access='sequential',
     1      form='unformatted', status='old', iostat=ios)
      call chopen (ios, trim(header)//'paths.bin', 'pathsd')
      read(3) nhead
      do 40  ihead = 1, nhead
         read(3)  head(ihead)
         read(3)  lhead(ihead)
   40 continue
c     Header lines above include carriage control
      read(3)  nat
      do 50  i = 0, nat
         read(3) (rat(j,i),j=1,3), ipot(i), i1b(i)
   50 continue

c     Initialize stuff...
c     nptot  number of total paths, incl all degeneracies
c     nuptot number of unique paths for which must calc xafs
c     ngs    number of generalized shells (unique distances)
      nptot = 0
      nuptot = 0
      ngs = 0
      xportx = eps5
      ndegx = 1
      c0lim = 1.0d10
      c1lim = 1.0d10
c     Initialize keep criterion
      xcalcx = -1

c     write output to paths.dat
      if (ipr2 .ne. 5)  then
         open (unit=1, file=trim(header)//'paths.dat',
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//'paths.dat', 'pathsd')
         do 60  ihead = 1, nhead
            write(1,58)  head(ihead)(1:lhead(ihead))
   58       format(a)
   60    continue
         write(1,61)  critpw
   61    format (' Plane wave chi amplitude filter', f7.2, '%')
         write(1,62)
   62    format (1x, 79('-'))
      endif

c     Write crit.dat (criteria information)
      if (ipr2 .ge. 1)  then
         open (unit=4, file=trim(header)//'crit.dat',
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//'crit.dat', 'pathsd')
         do 65  ihead = 1, nhead
            write(4,58)  head(ihead)(1:lhead(ihead))
   65    continue
         write(4,61)  critpw
         write(4,62)
         write(4,80)
   80    format (' ipath nleg ndeg     r       pwcrit    ',
     1           'xkeep   accuracy   xheap    accuracy')
      endif

c     Read path data for each total path length range

c     Prepare for first path.
      read(3,end=999)  r0, iout0

c     Begin next total path length range
      last = .false.
  100 continue
      ngs = ngs+1
      rcurr = r0
      np = 1
      do 110  i = 1,3
         iout(i,np) = iout0(i)
  110 continue
  120 read(3,end=140)  r0, iout0
         if (abs(r0-rcurr) .lt. eps3)  then
            np = np+1
            if (np .gt. np1x) then
               write(77,*) ' np, np1x ', np, np1x
               stop 'np > np1x'
            endif
            do 130  i = 1, 3
               iout(i,np) = iout0(i)
  130       continue
         else
c           r0 is the rtot for the next set
c           iout0 is the packed atom list for the first path of the
c           next set
            goto 200
         endif
      goto 120
  140 continue
c     Get here only if end-of-file during read
      last = .true.

  200 continue

      nupr = 0
c     variable nuprtt was nuprtot, changed to be six chars, SIZ 12/93
      nuprtt = 0

c     Hash each path into an integer
      iscale = 1000
      do 230  ip = 1, np

         npat = npatx
         call upack (iout(1,ip), npat, ipat)

c        Get hash key for this path.
c        If two paths are the same, except time-reversed, the xafs
c        will be the same, so check for this type of degeneracy.
c        We do this by choosing a 'standard order' for a path --
c        if it's the other-way-around, we time-reverse here.
         call timrep (npat, ipat, rx, ry, rz, dhash(ip))

  230 continue

c     Do a heap sort on these things
      call sortid (np, index, dhash)

c     Find beginning and end of range with same hash key
c     i0 is beginning of hash range, i1 is end of the range

      i0 = 1
  300 continue
         i1 = np + 1
         dcurr = dhash(index(i0))
         do 310  ip = i0+1, np
            if (dhash(index(ip)) .ne. dcurr)  then
c              end of a hash range
               i1 = ip
               goto 311
            endif
  310    continue
  311    continue
         i1 = i1-1

c        At this point, i0 is the first path and i1 the last
c        of a hash range.  Do whatever you want with them!

c        Sum degeneracy, including degeneracy from 1st bounce atom.
c        Check this range to see if all of the paths are actually 
c        degenerate.  Make sure time-ordering is standard.
         npat0 = npatx
         call upack (iout(1,index(i0)), npat0, ipat0)
         call timrep (npat0, ipat0, rx0, ry0, rz0, ddum)

         ndeg = 0
         do 430  ii = i0, i1
            npat = npatx
            call upack (iout(1,index(ii)), npat, ipat)
c           Note that if path gets time-reversed, we lose 1st bounce 
c           flag (since first atom is now last...), so save path deg
            ndpath = i1b(ipat(1))
            call timrep (npat, ipat, rx, ry, rz, ddum)
c           Sum degeneracy here.
            ndeg = ndeg + ndpath
c           Check for hash collisons begins here.
            ldiff = .false.
            if (npat .ne. npat0)  then
               ldiff = .true.
               goto 430
            endif
            do 320  iat = 1, npat
               if (ipot(ipat(iat)) .ne. ipot(ipat0(iat)))  then
                  ldiff = .true.
                  goto 400
               endif
  320       continue
            do 330  ileg = 1, npat
               if (abs(rx(ileg)-rx0(ileg)) .gt. eps3  .or.
     1             abs(ry(ileg)-ry0(ileg)) .gt. eps3  .or.
     2             abs(rz(ileg)-rz0(ileg)) .gt. eps3)  then
                  ldiff = .true.
                  goto 400
               endif
  330       continue
  400       continue
            if (ldiff)  then
              write(77,*) 'WARNING!!  Two non-degenerate paths hashed ',
     1                 'to the same hash key!!'
               write(77,*) dhash(index(i0)), dhash(index(ii))
               write(77,*) npat0, npat, '  npat0, npat'
               write(77,*) ' iat, ipot0, ipot, ipat0, ipat'
               do 410  iat = 1, npat
                  write(77,*) iat, ipot(ipat0(iat)), ipot(ipat(iat)),
     1                         ipat0(iat), ipat(iat)
  410          continue
               write(77,*) 'ileg, rx0,ry0,rz0,  rx1,ry1,rz1'
               do 420  ileg = 1, npat
                  write(77,*) ileg, rx0(ileg), rx(ileg)
                  write(77,*) ileg, ry0(ileg), ry(ileg)
                  write(77,*) ileg, rz0(ileg), rz(ileg)
  420          continue
               stop 'hash error'
            endif
  430    continue

c        Find path pw importance factors, and recalculate 
c        pathfinder crits for output
         call outcrt (npat0, ipat0, ckspc,
     1                nncrit, fbetac, ne, ik0, cksp, fbeta, 
     1                ipotnn, ipot,
     1                xport, xheap, xheapr, xkeep, xcalcx)

         if (xport*ndeg .gt. xportx*ndegx)  then
            xportx = xport
c           ndegx is degeneracy of path that makes xportx, used for
c           testing new path keep crit
            ndegx = ndeg
         endif
c        frac is fraction of max importance to use for test
         frac = 100*ndeg*xport/(ndegx*xportx)

c        Write output if path is important enough (ie, path is
c        at least critpw % important as most important path found
c        so far.)
         if (frac .ge. critpw)  then
            nupr = nupr+1
            nuprtt = nuprtt+ndeg
            nptot = nptot + ndeg
            nuptot = nuptot + 1

c           Write path info to paths.dat
c           mpprmd is double precision, used to get angles
c           180.000 instead of 179.983, etc.
            call mpprmd (npat0, ipat0, rid, betad, etad)
c           skip paths.dat if not necessary
            if (ipr2 .eq. 5)  goto 576
            write(1,500) nuptot, npat0+1, real(ndeg),
     1              rcurr/2
  500       format (1x, 2i5, f8.3,
     1             '  index, nleg, degeneracy, r=', f8.4)
            write(1,502)
  502       format ('      x           y           z     ipot  ',
     1              'label      rleg      beta        eta')
            do 510  i = 1, npat0
               iat = ipat0(i)
               write(1,506)  rat(1,iat), rat(2,iat),
     1                  rat(3,iat), ipot(iat), potlbl(ipot(iat)),
     1                  rid(i), betad(i)*raddeg, etad(i)*raddeg
  506          format (3f12.6, i4, 1x, '''', a6, '''', 1x, 3f10.4)
  510       continue
            write(1,506)  rat(1,0), rat(2,0), rat(3,0), ipot(0), 
     1         potlbl(ipot(0)),
     1         rid(npat0+1), betad(npat0+1)*raddeg, etad(npat0+1)*raddeg
c           End of paths.dat writing for this path

c           Write to crit.dat here (unit 4, opened above)
  576       continue

c           cmpk is degeneracy corrected xkeep, should equal frac
            cmpk = xkeep*ndeg/ndegx
c           cmpk is accuracy of xkeep, 100 is perfect
            cmpk = 100.0d0 - 100.0d0*(abs(frac-cmpk)/frac)

c           cmph is same thing for xheap
            if (xheap .lt. 0.0d0)  then
               cmph = 100.0d0
            else
               cmph = xheap*ndeg/ndegx
               cmph = 100.0d0 - 100.0d0*(abs(frac-cmph)/frac)
            endif

            if (ipr2 .ge. 1)  then
               write(4,560)  nuptot, npat0+1, ndeg, rcurr/2, frac,
     1             xkeep, cmpk, xheap, cmph
  560          format (i6, i4, i6, 3f10.4, f8.2, f10.4, 1pe14.3)
            endif

c           write out fraction error between xkeep and critpw
         endif

c        And do next ihash range
         i0 = i1+1
      if (i0 .le. np)  goto 300

c     print 600,  ngs, rcurr, nupr
  600 format (1x, i5, f12.6, i7, ' igs, rcurr, nupr')
c     write(80,601)  ngs, rcurr/2, nupr, nuprtt
  601 format (1x, i8, f12.6, 2i9)

      if (.not. last) goto 100

      if (ipr2 .ne. 5)  close (unit=1)
c     delete paths.bin when done...
      close (unit=3, status='delete')
      close (unit=4)

      write(77,620) nuptot, nptot
  620 format ('    Unique paths', i7, ',  total paths', i8)

c     Do not let user accidently fill up their disk
      if (nuptot .gt. 3200)  then
      write(77,*) 'You have found more than 1200 paths.  Genfmt'
      write(77,*) 'could require a lot of time and more than 6 meg of'
      write(77,*) 'storage.  Suggest a larger critpw to reduce number'
      write(77,*) 'of paths.  To continue this calculation, restart'
      write(77,*) 'with current paths.dat and module genfmt (3rd module'
      write(77,*) 'on CONTROL card).'
      stop 'User must verify very large run.'
      endif
      return
  999 stop 'no input'
      end
c     Periodic table of the elements
c     Written by Steven Zabinsky, Feb 1992.  Deo Soli Gloria

c     atwts(iz)  single precision fn, returns atomic weight
c     atwtd(iz)  double precision fn, returns atomic weight
c     atsym(iz)  character*2 fn, returns atomic symbol

      double precision function atwtd(iz)
      implicit double precision (a-h, o-z)
      double precision weight
      save /atwtco/
      common /atwtco/ weight(103)
      atwtd = weight(iz)
      return
      end

      real*8 function atwts(iz)
      implicit double precision (a-h, o-z)
      double precision weight
      save /atwtco/
      common /atwtco/ weight(103)
      atwts = weight(iz)
      return
      end

      character*2 function atsym (iz)
      implicit double precision (a-h, o-z)
      character*2 sym
      save /atsyco/
      common /atsyco/ sym(103)
      atsym = sym(iz)
      return
      end

      block data prtbbd
c     PeRiodic TaBle Block Data

c     Atomic weights from inside front cover of Ashcroft and Mermin.

      double precision weight
      save /atwtco/
      common /atwtco/ weight(103)

      character*2 sym
      save /atsyco/
      common /atsyco/ sym(103)

      data weight /
     1   1.0079d0, 4.0026d0, 6.941d0,  9.0122d0, 10.81d0,  12.01d0,
     2   14.007d0, 15.999d0, 18.998d0, 20.18d0,  22.9898d0, 24.305d0,
     3   26.982d0, 28.086d0, 30.974d0, 32.064d0, 35.453d0, 39.948d0,
     4   39.09d0,  40.08d0,  44.956d0, 47.90d0,  50.942d0, 52.00d0,
     5   54.938d0, 55.85d0,  58.93d0,  58.71d0,  63.55d0,  65.38d0,
     6   69.72d0,  72.59d0,  74.922d0, 78.96d0,  79.91d0,  83.80d0,
     7   85.47d0,  87.62d0,  88.91d0,  91.22d0,  92.91d0,  95.94d0,
     8   98.91d0,  101.07d0, 102.90d0, 106.40d0, 107.87d0, 112.40d0,
     9   114.82d0, 118.69d0, 121.75d0, 127.60d0, 126.90d0, 131.30d0,
     x   132.91d0, 137.34d0, 138.91d0, 140.12d0, 140.91d0, 144.24d0,
     1   145.0d0,  150.35d0, 151.96d0, 157.25d0, 158.92d0, 162.50d0,
     2   164.93d0, 167.26d0, 168.93d0, 173.04d0, 174.97d0, 178.49d0,
     3   180.95d0, 183.85d0, 186.2d0,  190.20d0, 192.22d0, 195.09d0,
     4   196.97d0, 200.59d0, 204.37d0, 207.19d0, 208.98d0, 210.0d0,
     5   210.0d0,  222.0d0,  223.0d0,  226.0d0, 227.0d0,   232.04d0,
     6   231.0d0,  238.03d0, 237.05d0, 244.0d0, 243.0d0,   247.0d0,
     7   247.0d0,  251.0d0,  254.0d0, 257.0d0, 256.0d0,    254.0d0,
     8   257.0d0/

      data sym /  'H', 'He','Li','Be','B', 'C', 'N', 'O', 'F', 'Ne',
     1            'Na','Mg','Al','Si','P', 'S', 'Cl','Ar','K', 'Ca',
     2            'Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y', 'Zr',
     4            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5            'Sb','Te','I', 'Xe','Cs','Ba','La','Ce','Pr','Nd',
     6            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7            'Lu','Hf','Ta','W', 'Te','Os','Ir','Pt','Au','Hg',
     8            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9            'Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     x            'Md','No','Lw'/

      end
      subroutine phase (iph, nr, dx, x0, ri, ne, em, edge,
     1                  index, rmt, xmu, vi0, rs0, gamach,
     2                  vtot, edens,
     3                  eref, ph, lmax)

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


c     INPUT
c     iph          unique pot index (used for messages only)
c     nr, dx, x0, ri(nr)
c                  Loucks r-grid, ri=exp((i-1)*dx-x0)
c     ne, em(ne)   number of energy points, real energy grid
c     edge         energy for k=0 (note, edge=xmu-vr0)
c     index        0  Hedin-Lunqist + const real & imag part
c                  1  Dirac-Hara + const real & imag part
c                  2  ground state + const real & imag part
c                  3  Dirac-Hara + HL imag part + const real & imag part
c                  4, 5, 6, see rdinp or xcpot
c     rmt          r muffin tin
c     xmu          fermi level
c     vi0          const imag part to add to complex potential
c     rs0          user input density cutoff, used only with ixc=4
c     gamach       core hole lifetime
c     vtot(nr)     total potential, including gsxc
c     edens(nr)    density
c
c     OUTPUT
c     eref(ne)     complex energy reference including energy dep xc
c     ph(nex,ltot+1) complex scattering phase shifts
c     lmax         max l (lmax = kmax*rmt)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)

      dimension   ri(nr), em(nex), vtot(nr), edens(nr)
      complex*16  eref(nex)
      complex*16  ph(nex,ltot+1)

c     work space for xcpot
      dimension   vxcrmu(nrptx), vxcimu(nrptx)
c     work space for fovrg
      complex*16 p(nrptx), q(nrptx), ps(nrptx), qs(nrptx), vm(nrptx)

      complex*16  p2, xkmt, temp, dny, pu, qu
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 v(nrptx)
      external besjn

      ihard = 0
c     zero phase shifts (some may not be set below)
      do 100  ie = 1, ne
         do 90  il = 1, ltot+1
            ph(ie,il) = dcmplx(0.0d0,0.0d0)
   90    continue
  100 continue

c     limit l, lmax = kmax * rmt
c     lmax = rmt * sqrt(em(ne)-edge)
c     Use kmax = 20 so we get enough l-points even if kmax is small
      lmax = rmt * (20 * bohr)
      lmax = min (lmax, ltot)

c     set imt and jri (use general Loucks grid)
c     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = (log(rmt) + x0) / dx  +  1
      jri = imt+1
      if (jri .gt. nr)  stop 'jri .gt. nr in phase'
c     xmt is floating point version of imt, so that
c     rmt = (exp (x-1)*dx - x0).  xmt used in fovrg
      xmt = (log(rmt) + x0) / dx  +  1

      ifirst = 0
c     calculate phase shifts
      do 220 ie = 1, ne

         call xcpot (iph, ie, nr, index, ifirst, jri,
     1               em(ie), xmu, vi0, rs0, gamach,
     2               vtot, edens,
     3               eref(ie), v,
     4               vxcrmu, vxcimu)

c        fovrg needs v in form pot*r**2
         do 120  i = 1, jri
            v(i) = v(i) * ri(i)**2
  120    continue

c        p2 is (complex momentum)**2 referenced to energy dep xc
         p2 = em(ie) - eref(ie)
         xkmt = rmt * sqrt (p2)
         call besjn (xkmt, jl, nl)

         do 210  il = 1, lmax+1
            l = il - 1

            call fovrg(il, ihard, rmt, xmt, jri, p2, 
     1                  nr, dx, ri, v, dny,
     1                  pu, qu, p, q, ps, qs, vm)


            temp = (jl(il)*(dny-l) + xkmt*jl(il+1))  /
     1             (nl(il)*(dny-l) + xkmt*nl(il+1))
            xx = dble (temp)
            yy = dimag(temp)
            if (xx .ne. 0)  then
               alph = (1 - xx**2 - yy**2)
               alph = sqrt(alph**2 + 4*xx**2) - alph
               alph = alph / (2 * xx)
               alph = atan (alph)
            else
               alph = 0
            endif
            beta = (xx**2 + (yy+1)**2) /
     1             (xx**2 + (yy-1)**2)
            beta = log(beta) / 4

            ph(ie,il) = dcmplx (alph, beta)

c           cut phaseshift calculation if they become too small
            if (abs(ph(ie,il)) .lt. 1.0d-6)  goto 220

  210    continue

  220 continue


c     Warn user if fovrg failed ihard test.
      if (ihard .gt. 0)  then
         write(77,*) ' Hard test failed in fovrg ', ihard, ' times.'
         write(77,*) ' Muffin-tin radius may be too large;',
     1               ' coordination number too small.'
      endif

      return
      end
      subroutine phash (npat, ipat, rx, ry, rz, dhash)
c     hashes a path into double precision real dhash
      implicit double precision (a-h, o-z)


      character*72 header
      common /header_common/ header


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      double precision dhash
      dimension rx(npatx), ry(npatx), rz(npatx), ipat(npatx+1)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      double precision xx

      parameter (iscale = 1000)
      parameter (factor = 16.12345678d0)

c     Hashing scheme: Assume about 15 significant digits in a double 
c     precision number.  This is 53 bit mantissa and 11 bits for sign 
c     and exponent, vax g_floating and probably most other machines.
c     With max of 9 legs, 47**9 = 1.12e15, so with a number less than 
c     47, we can use all these digits, scaling each leg's data by 
c     47**(j-1).  Actually, since our numbers can go up to about 10,000,
c     we should keep total number < 1.0e11, 17**9 = 1.18e11, which means
c     a factor a bit less than 17.  Choose 16.12345678, a non-integer,
c     to help avoid hash collisions.

c     iscale and 'int' below are to strip off trailing digits, which
c     may contain roundoff errors

      dhash = 0
      do 210  j = 1, npat
         xx = factor**(j-1)
         dhash = dhash + xx * (nint(rx(j)*iscale) +
     1               nint(ry(j)*iscale)*0.894375 +
     2               nint(rz(j)*iscale)*0.573498)
  210 continue
      do 220  j = 1, npat
         xx = factor**(j-1)
         dhash = dhash + xx * ipot(ipat(j))
  220 continue
      dhash = dhash + npat * 40 000 000

      return
      end
c     make e and r mesh for phase
c     input:  nr, dx, x0, nemax, iprint,
c             ixanes, edge, xmu, vint, vr0, imt, edens, nph
c             edge, xmu... used only with ixanes = 1
c     output: ri(nr), ne, em(ne), ik0 [grid point with k=0]
c
c     set nemax = nex (from dim.h) for max number of points

      subroutine phmesh (nr, dx, x0, nemax, iprint,
     1                   ixanes, edge, xmu, vint, vr0,
     1                   imt, edens, nph,
     2                   ri, ne, em, ik0)
      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      dimension ri(nr), em(nex)

c     edens       overlapped density*4*pi
c     imt         r mesh index just inside rmt
c     see arrays.h
      dimension edens(nrptx,0:nphx)
      dimension imt(0:nphx)

c     r mesh
      do 100  i = 1, nr
         ri(i) = rr(i)
  100 continue

c     xkmin needed only with ixanes
      if (ixanes .gt. 0)  then
c        Need xf**2 min for all unique potentials, take rho(imt) as
c        min rho
         xf2int = xmu-vint
         xf2min = xf2int
         do 400  i = 0, nph
            rs = (3 / edens(imt(i),i)) ** third
            xf2 = (fa / rs) ** 2
            if (xf2 .le. xf2min) xf2min = xf2
  400    continue

         xkmin2 = xf2min - vr0
         if (xkmin2 .lt. 0)  then
            write(77,*) ' xf2min, vr0, xkmin2'
            write(77,*) xf2min, vr0, xkmin2
            write(77,*) 'bad vr0 in phmesh'
            stop 'bad vr0 in phmesh'
         endif

         delk = bohr/5.0d0
         xkmin = sqrt (xkmin2)
         n = int(xkmin/delk) - 1
      else
         xkmin = 0.0d0
         n = 0
      endif

c     energy mesh
c      n pts (-2 le k lt 0,  delk=0.2 ang(-1) ) (only if xanes)
c     30 pts (0 le k le 5.8, delk=0.2 ang(-1) )
c      9 pts (6 le k le 10., delk=0.5 ang(-1) )
c     10 pts (11 le k le 20.0, delk=1.0 ang(-1) )
      ne = 0
      delk = bohr/5.0d0
      if (ixanes .gt. 0)  then
         xkmin = n*delk
         do 110 i=1,n
            tempk=-xkmin+(i-1)*delk
            ne = ne+1
            em(ne)=-tempk**2+edge
  110    continue
      endif
      delk = bohr/5
      do 112 i=1,30
         tempk=(i-1)*delk
         ne = ne+1
         em(ne)=tempk**2+edge
         if (i.eq.1)  ik0 = ne
  112 continue
      delk = bohr/2
      do 113 i=1,9
         tempk=6.*bohr + (i-1)*delk
         ne = ne+1
         em(ne)=tempk**2+edge
  113 continue
      delk=bohr
      do 114 i=1,10
         tempk=11.0d0*bohr + (i-1)*delk
         ne = ne+1
         em(ne)=tempk**2+edge
  114 continue

c     print*, 'phmesh: ne, nex, nemax before setting ne ',
c    1                 ne, nex, nemax
      ne = min (ne, nemax)
c     print*, 'phmesh: ne, nex, nemax after  setting ne ',
c    1                 ne, nex, nemax


      if (iprint .ge. 3)  then
         open (unit=44, file=trim(header)//'emesh.dat')
         write(44,*) 'edge, bohr, edge*ryd ', edge, bohr, edge*ryd
         write(44,*) 'ixanes, ik0 ', ixanes, ik0
         write(44,*) vint, xkmin, n, ' vint, xkmin, n'
         write(44,*) 'ie, em(ie), xk(ie)'
         do 230  ie = 1, ne
            write(44,220)  ie, em(ie), getxk(em(ie)-edge)/bohr
  220       format (i5, 2f20.5)
  230    continue
         close (unit=44)
      endif

      return
      end
      subroutine pijump (ph, old)
      implicit double precision (a-h, o-z)


      character*72 header
      common /header_common/ header

c     removes jumps of 2*pi in phases

c     ph = current value of phase (may be modified on output, but
c          only by multiples of 2*pi)
c     old = previous value of phase


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)

      parameter (twopi = 2 * pi)
      dimension xph(3)

      xph(1) = ph - old
      jump =  (abs(xph(1))+ pi) / twopi
      xph(2) = xph(1) - jump*twopi
      xph(3) = xph(1) + jump*twopi


      xphmin = min (abs(xph(1)), abs(xph(2)), abs(xph(3)))
      isave = 0
      do 10  i = 1, 3
         if (abs (xphmin - abs(xph(i))) .le. 0.01)  isave = i
   10 continue
      if (isave .eq. 0)  then
         write(77,*) 'isave ', isave
         write(77,*) xph(1)
         write(77,*) xph(2)
         write(77,*) xph(3)
         stop 'pijump'
      endif

      ph = old + xph(isave)

      return
      end
      subroutine potph (isporb)

c     Cluster code -- multiple shell single scattering version of FEFF
c     This program (or subroutine) calculates potentials and phase
c     shifts for unique potentials specifed by atoms and overlap cards.
c
c     Input files:  potph.inp    input data, atoms, overlaps, etc.
c     Output:       phases.bin   phase shifts for use by the rest of the
c                                program
c                   xxx.dat      various diagnostics

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt



c     Notes:
c        nat	number of atoms in problem
c        nph	number of unique potentials
c        nfr	number of unique free atoms
c        ihole	hole code of absorbing atom
c        iph=0 for central atom
c        ifr=0 for central atom

c     Specific atom input data
      dimension iphat(natx)	!given specific atom, which unique pot?
      dimension rat(3,natx)	!cartesian coords of specific atom

c     Unique potential input data
      dimension iatph(0:nphx)	!given unique pot, which atom is model?
				!(0 if none specified for this unique pot)
      dimension ifrph(0:nphx)	!given unique pot, which free atom?
      dimension xnatph(0:nphx)	!given unique pot, how many atoms are there 
				!of this type? (used for interstitial calc)
      character*6 potlbl(0:nphx)	!label for user convienence

      dimension folp(0:nphx)	!overlap factor for rmt calculation
      dimension novr(0:nphx)	!number of overlap shells for unique pot
      dimension iphovr(novrx,0:nphx)	!unique pot for this overlap shell
      dimension nnovr(novrx,0:nphx)	!number of atoms in overlap shell
      dimension rovr(novrx,0:nphx)	!r for overlap shell

c     Free atom data
      dimension ion(0:nfrx)	!ionicity, input
      dimension iz(0:nfrx)	!atomic number, input

c     ATOM output
c     Note that ATOM output is dimensioned 251, all other r grid
c     data is set to nrptx, currently 250
      dimension rho(251,0:nfrx)		!density*4*pi
      dimension vcoul(251,0:nfrx)	!coulomb potential

c     Overlap calculation results
      dimension edens(nrptx,0:nphx)	!overlapped density*4*pi
      dimension vclap(nrptx,0:nphx) 	!overlapped coul pot
      dimension vtot (nrptx,0:nphx)	!overlapped total potential

c     Muffin tin calculation results
      dimension imt(0:nphx)	!r mesh index just inside rmt
      dimension inrm(0:nphx)	!r mesh index just inside rnorman
      dimension rmt(0:nphx)	!muffin tin radius
      dimension rnrm(0:nphx)	!norman radius

c     PHASE output
      complex*16 eref(nex)		!interstitial energy ref
      complex*16 ph(nex,ltot+1,0:nphx)	!phase shifts
      dimension lmax(0:nphx)		!number of ang mom levels

      common /print/ iprint

      parameter (nheadx = 30)
      character*80 head(nheadx)
      dimension lhead(nheadx)

c     head0 is header from potph.dat, include carriage control
      character*80 head0(nheadx)
      dimension lhead0(nheadx)

      dimension em(nex)
      dimension dgc0(251), dpc0(251)
      dimension xsec(nex), xsatan(nex)

c     nrx = max number of r points for phase r grid
      parameter (nrx = 250)
      dimension ri(nrptx), vtotph(nrx), rhoph(nrx)

   10 format (4x, a, i5)

c     Read input from file potph.inp
      open (unit=1, file=trim(header)//'potph.dat',
     >       status='old', iostat=ios)
      call chopen (ios, trim(header)//'potph.dat', 'potph')
      nhead0 = nheadx
      call rpotph (1, nhead0, head0, lhead0, nat, nph,
     1             nfr, ihole, gamach, iafolp, intclc,
     1             ixc, vr0, vi0, rs0, iphat, rat, iatph, ifrph, 
     1             xnatph, novr,
     2             iphovr, nnovr, rovr, folp, ion, iz, iprint, 
     2             ixanes, nemax, xkmin, xkmax, potlbl)
      close (unit=1)

c     Free atom potentials and densities
c     NB wsatom is needed in SUMAX, if changed here, change it there
      wsatom = 15
c     do not save spinors
      ispinr = 0
      do 20  ifr = 0, nfr
         itmp = 0
         if (ifr .eq. 0)  itmp = ihole
         write(77,10) 'free atom potential and density for atom type', ifr
         call feff_atom(head0(1)(1:40), ifr, iz(ifr), itmp, wsatom,
     1              ion(ifr), vcoul(1,ifr), rho(1,ifr),
     2              ispinr, dgc0, dpc0, et)
c        etfin is absorbing atom final state total energy
c        etinit is absorbing atom initial state (no hole)
         if (ifr .eq. 0)  etfin = et
   20 continue
      if (ixanes .gt. 0)  then
         write(77,10) 'initial state energy'
c        save spinor for core hole orbital
         ispinr = ihole
c        if no hole, use orbital from isporb
         if (ihole .eq. 0)  ispinr = isporb
         itmp = 0
         call feff_atom (head0(1)(1:40), 0, iz(0), itmp, wsatom,
     1              ion(0), vcoul(1,nfr+1), rho(1,nfr+1),
     2              ispinr, dgc0, dpc0, etinit)
      endif
c     Need etfin if xanes and no hole, use K shell for this
      if (ixanes .gt. 0 .and. ihole .eq. 0)  then
c        K hole
         itmp = 1
         ispinr = 0
         call feff_atom (head0(1)(1:40), 0, iz(0), itmp, wsatom,
     1              ion(0), vcoul(1,nfr+1), rho(1,nfr+1),
     2              ispinr, dgc0, dpc0, etfin)
      endif

c     Overlap potentials and densitites
      do 40  iph = 0, nph
         write(77,10)
     1    'overlapped potential and density for unique potential', iph
         call ovrlp (iph, iphat, rat, iatph, ifrph, novr,
     1               iphovr, nnovr, rovr, iz, nat, rho, vcoul,
     2               edens, vclap, rnrm)
   40 continue

c     Find muffin tin radii, add gsxc to potentials, and find
c     interstitial parameters
      write(77,10) 'muffin tin radii and interstitial parameters'
      call istprm (nph, nat, iphat, rat, iatph, xnatph,
     1             novr, iphovr, nnovr, rovr, folp, edens,
     2             vclap, vtot, imt, inrm, rmt, rnrm, rhoint,
     3             vint, rs, xf, xmu, rnrmav, intclc)

c     Automatic max reasonable overlap
      if (iafolp .eq. 1)  then
         write(77,10) 'automatic overlapping'
         write(77,*) 'iph, rnrm(iph)*bohr, rmt(iph)*bohr, folp(iph)'
         do 400  iph = 0, nph
            folp(iph) = 1 + 0.7*(rnrm(iph)/rmt(iph) - 1)
            write(77,*) iph, rnrm(iph)*bohr, rmt(iph)*bohr, folp(iph)
  400    continue
         call istprm (nph, nat, iphat, rat, iatph, xnatph,
     1                novr, iphovr, nnovr, rovr, folp, edens,
     2                vclap, vtot, imt, inrm, rmt, rnrm, rhoint,
     3                vint, rs, xf, xmu, rnrmav, intclc)
      endif

c     Initialize header routine and write misc.dat
      call sthead (nhead0, head0, lhead0, nph, iz, rmt, rnrm,
     1             ion, ifrph, ihole, ixc,
     2             vr0, vi0, rs0, gamach, xmu, xf, vint, rs,
     3             nhead, lhead, head)
      if (iprint .ge. 1)  then
         open (unit=1, file=trim(header)//'misc.dat',
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//'misc.dat', 'potph')
         call wthead(1)
         close (unit=1)
      endif

      if (iprint .ge. 2)  then
         call wpot (nph, edens, ifrph, imt, inrm,
     1              rho, vclap, vcoul, vtot)
      endif

c     Phase shift calculation
c     Make energy mesh and position grid
      nr = 250
      dx = .05
      x0 = 8.8
      edge = xmu - vr0
      call phmesh (nr, dx, x0, nemax, iprint,
     1             ixanes, edge, xmu, vint, vr0,
     1             imt, edens, nph,
     2             ri, ne, em, ik0)

c     Cross section calculation, use phase mesh for now
c     remove xanes calculation in feff6l

      do 60  iph = 0, nph
         write(77,10) 'phase shifts for unique potential', iph
c        fix up variable for phase
         call fixvar (rmt(iph), edens(1,iph), vtot(1,iph),
     1                vint, rhoint, nr, dx, x0, ri,
     2                vtotph, rhoph)

         call phase (iph, nr, dx, x0, ri, ne, em, edge,
     1               ixc, rmt(iph), xmu, vi0, rs0, gamach,
     2               vtotph, rhoph,
     3               eref, ph(1,1,iph), lmax(iph))
   60 continue

      if (iprint .ge. 2)  then
         call wphase (nph, em, eref, lmax, ne, ph)
      endif

c     Write out phases for genfmt
c     May need stuff for use with headers only
      open (unit=1, file=trim(header)//'phase.bin', access='sequential',
     1      form='unformatted', status='unknown', iostat=ios)
      call chopen (ios, trim(header)//'phase.bin', 'potph')
      write(1) nhead
      do 62  i = 1, nhead
         write(1) head(i)
         write(1) lhead(i)
   62 continue
      write(1) ne, nph, ihole, rnrmav, xmu, edge, ik0
      write(1) (em(ie),ie=1,ne)
      write(1) (eref(ie),ie=1,ne)
      do 80  iph = 0, nph
         write(1) lmax(iph), iz(ifrph(iph))
         write(1) potlbl(iph)
         do 70  ie = 1, ne
            write(1)  (ph(ie,ll,iph), ll=1,lmax(iph)+1)
   70    continue
   80 continue
      close (unit=1)

      return
      end
      subroutine potsl (dv,d,dp,dr,dpas,dexv,z,np,ion,icut,dvn)
c
c coulomb potential uses a 4-point integration method
c dv=potential;  d=density;  dp=bloc de travail; dr=radial mesh;
c dpas=exponential step; dexv=multiplicative coefficient for the exchang
c z=atomic number;  np=number of points; ion=z-number of electrons
c if icut is zero one corrects the potential by -(ion+1)/r
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      dimension dv(251), d(251), dp(251), dr(251), dvn(251)
      das=dpas/24.0d0
      do 10 i=1,np
   10 dv(i)=d(i)*dr(i)
      dlo=exp(dpas)
      dlo2=dlo*dlo
      dp(2)=dr(1)*(d(2)-d(1)*dlo2)/(12.0*(dlo-1.0))
      dp(1)=dv(1)/3.0-dp(2)/dlo2
      dp(2)=dv(2)/3.0-dp(2)*dlo2
      j=np-1
      do 20 i=3,j
   20 dp(i)=dp(i-1)+das*(13.0*(dv(i)+dv(i-1))-(dv(i-2)+dv(i+1)))
      dp(np)=dp(j)
      dv(j)=dp(j)
      dv(np)=dp(j)
      do 30 i=3,j
      k=np+1-i
   30 dv(k)=dv(k+1)/dlo+das*(13.0*(dp(k+1)/dlo+dp(k))-(dp(k+2)/dlo2+dp
     1 (k-1)*dlo))
      dv(1)=dv(3)/dlo2+dpas*(dp(1)+4.0*dp(2)/dlo+dp(3)/dlo2)/3.0
      dlo=-(ion+1)
      do 40 i=1,np
      dvn(i)=dv(i)/dr(i)
      dv(i)=dv(i)-(z+exchan(d(i),dr(i),dexv))
      if (icut.ne.0) go to 40
      if (dv(i).gt.dlo) dv(i)=dlo
   40 dv(i)=dv(i)/dr(i)
      return
      end
      subroutine potslw (dv,d,dp,dr,dpas,np)
c
c coulomb potential uses a 4-point integration method
c dv=potential;  d=density;  dp=bloc de travail; dr=radial mesh
c dpas=exponential step;
c np=number of points
c **********************************************************************

      implicit double precision (a-h,o-z)
      save
      dimension dv(251), d(251), dp(251), dr(251)
      das=dpas/24.0d0
      do 10 i=1,np
   10 dv(i)=d(i)*dr(i)
      dlo=exp(dpas)
      dlo2=dlo*dlo
      dp(2)=dr(1)*(d(2)-d(1)*dlo2)/(12.0*(dlo-1.0))
      dp(1)=dv(1)/3.0d0-dp(2)/dlo2
      dp(2)=dv(2)/3.0d0-dp(2)*dlo2
      j=np-1
      do 20 i=3,j
   20 dp(i)=dp(i-1)+das*(13.0d0*(dv(i)+dv(i-1))-(dv(i-2)+dv(i+1)))
      dp(np)=dp(j)
      dv(j)=dp(j)
      dv(np)=dp(j)
      do 30 i=3,j
      k=np+1-i
   30 dv(k)=dv(k+1)/dlo+das*(13.0d0*(dp(k+1)/dlo+dp(k))-(dp(k+2)/dlo2+dp
     1 (k-1)*dlo))
      dv(1)=dv(3)/dlo2+dpas*(dp(1)+4.0d0*dp(2)/dlo+dp(3)/dlo2)/3.0d0
      do 40 i=1,np
   40 dv(i)=dv(i)/dr(i)
      return
      end
      subroutine prcrit (neout, nncrit, ik0out, cksp, fbeta, ckspc, 
     1                   fbetac, potlb0)
      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


c     Prepare fbeta arrays, etc., for pathfinder criteria
c
c     Note that path finder is single precision, so be sure that
c     things are correct precision in calls and declarations!
c     See declarations below for details.
c     
c     Inputs:  Reads phase.bin
c     Output:  neout   'ne', number of energy grid points
c              ik0out  index of energy grid with k=0
c              cksp    |p| at each energy grid point in single precision
c              fbeta   |f(beta)| for each angle, npot, energy point, sp
c              ckspc   |p| at each necrit point in single precision
c              fbetac  |f(beta)| for each angle, npot, nncrit point, sp
c              potlb0  unique potential labels


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines


c     Output variables SINGLE PRECISION for use with path finder.
c     BE CAREFUL!!
      parameter (necrit=9, nbeta=40)
      real*8 fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)
      real*8 fbeta(-nbeta:nbeta,0:npotx,nex), cksp(nex)
      character*6  potlb0(0:npotx)

c     Local variables
      complex*16 cfbeta, tl
      dimension dcosb(-nbeta:nbeta)
      dimension pl(ltot+1)
      dimension iecrit(necrit)


c     Need stuff from phase.bin
c     Read phase calculation input, data returned via commons
      open (unit=1, file=trim(header)//'phase.bin', status='old',
     1      access='sequential', form='unformatted', iostat=ios)
      call chopen (ios, trim(header)//'phase.bin', 'prcrit')
      call rphbin (1)
      close (unit=1)
c     Pass out ne, ik0, potlbl (from rphbin via /pdata/)
      neout = ne
      ik0out = ik0
      do 40  i = 0, npotx
         potlb0(i) = potlbl(i)
   40 continue

c     |p| at each energy point (path finder uses invA, convert here)
      do 100  ie = 1, ne
         cksp(ie) = abs (sqrt (em(ie) - eref(ie))) / bohr
  100 continue

c     Make the cos(beta)'s
c     Grid is from -40 to 40, 81 points from -1 to 1, spaced .025
      do 200  ibeta = -nbeta, nbeta
         dcosb(ibeta) = 0.025d0 * ibeta
  200 continue
c     watch out for round-off error
      dcosb(-nbeta) = -1
      dcosb(nbeta)  =  1

c     make fbeta (f(beta) for all energy points
      do 280  ibeta = -nbeta, nbeta
         call cpl0 (dcosb(ibeta), pl, lmaxp1)
         do 260  iii = 0, npot
            do 250  ie = 1, ne
               cfbeta = 0
               do 245  il = 1, lmax(ie,iii)+1
                  tl = (exp(2.0d0*coni*ph(ie,il,iii)) - 1.0d0)/(2*coni)
                  cfbeta = cfbeta + tl*pl(il)*(2*il-1)
  245          continue
               fbeta(ibeta,iii,ie) = abs(cfbeta)
  250       continue
  260    continue
  280 continue

c     Make similar arrays for only the icrit points

c     Use 9 points at k=0,1,2,3,4,6,8,10,12 invA
c     See phmesh for energy gid definition.  These seem to work fine, 
c     and results aren't too sensitive to choices of k.  As few as 4
c     points work well (used 0,3,6,9), but time penalty for 9 points
c     is small and increased safety seems to be worth it.
      iecrit(1) = ik0
      iecrit(2) = ik0 + 5
      iecrit(3) = ik0 + 10
      iecrit(4) = ik0 + 15
      iecrit(5) = ik0 + 20
      iecrit(6) = ik0 + 30
      iecrit(7) = ik0 + 34
      iecrit(8) = ik0 + 38
      iecrit(9) = ik0 + 40

c     make sure that we have enough energy grid points to use all
c     9 iecrits
      nncrit = 0
      do 290  ie = 1, necrit
         if (iecrit(ie) .gt. ne)  goto 295
         nncrit = ie
  290 continue
  295 continue
      if (nncrit .eq. 0) stop 'bad nncrit in prcrit'
      write(77,*) ' nncrit in prcrit ', nncrit
            

      do 320  icrit = 1, nncrit
         ie = iecrit(icrit)
         ckspc(icrit) = cksp(ie)
         do 310  ibeta = -nbeta, nbeta
            do 300  iii = 0, npot
               fbetac(ibeta,iii,icrit) = fbeta(ibeta,iii,ie)
  300       continue
  310    continue
  320 continue

      return
      end
      subroutine quinn (x, rs, wp, ef, ei)
      implicit double precision (a-h, o-z)

c     input  x, rs, wp, ef
c     output ei

c***********************************************************************
c
c     quinn: calculates low energy gamma (approx. proportional to e**2)
c             formula taken from john j. quinn, phys. rev. 126,
c             1453 (1962); equation (7).
c             a cut-off is set up at quinn's cutoff + ef = ekc; it is a
c             rounded inverted step function (a fermi function)
c             theta = 1/( 1 + exp((e-ekc)/gam)) )
c             where the rounding factor gam is set to be about 0.3 ekc.
c     modified by j. rehr (oct 1991) based on coding of r. albers
c     subroutines quinn.f and quinnc.f
c
c     variables:
c        x  = p/pf
c        rs = ws density parameter
c        ei = imaginary self energy
c        pfqryd = quinn's prefactor in atomic-rydberg units
c        wkc = quinn's plasmon threshold
c
c***********************************************************************

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (alphaq = 1.0d0/ fa)

c     calculate quinn prefactor in atomin Hartree units
      pisqrt = sqrt(pi)
      pfq = pisqrt / (32.0d0 * (alphaq*rs)**1.5d0)
      temp1 = atan (sqrt (pi / (alphaq*rs)))
      temp2 = sqrt(alphaq*rs/pi) / (1 + alphaq*rs/pi)
      pfq = pfq * (temp1 + temp2)

c     calculate quinn cutoff
c     wkc = quinn's plasmon threshold
c     wkc is cut-off of quinn, pr126, 1453, 1962, eq. (11)
c     in formulae below wp=omegap/ef
      wkc = (sqrt(1+wp) - 1)**2
      wkc = (1 + (6.0d0/5.0d0) * wkc / wp**2) * wp * ef

c     we add fermi energy to get correct energy for
c     plasma excitations to turn on
      ekc = wkc + ef

c     calculate gamma
c     gamryd = 2 * (pfqryd/x) * (x**2-1)**2
      gam = (pfq/x) * (x**2-1)**2

c     put in fermi function cutoff
      eabs = ef * x**2
      arg = (eabs-ekc) / (0.3d0*ekc)
      f = 0
      if (arg .lt. 80)  f = 1.0d0 / (1.0d0 + exp(arg))

      ei = -gam * f / 2.0d0

      return
      end
      subroutine rdhead (io, nhead, head, lhead)
      implicit double precision (a-h, o-z)


      character*72 header
      common /header_common/ header

c     Reads title line(s) from unit io.  Returns number of lines
c     read.  If more than nheadx lines, skips over them.  End-of-header
c     marker is a line of 1 blank, 79 '-'s.
c     lhead is length of each line w/o trailing blanks.
c     header lines returned will have 1st space on line blank for
c     carriage control

      character*(*) head(nhead)
      dimension lhead(nhead)
      character*80  line

      n = 0
      nheadx = nhead
      nhead = 0
   10 read(io,20)  line
   20    format(a)
         if (line(4:11) .eq. '--------')  goto 100
         n = n+1
         if (n .le. nheadx)  then
            head(n) = line
            lhead(n) = istrln(head(n))
            nhead = n
         endif
      goto 10
  100 continue
      return
      end
      subroutine rdinp(mphase, mpath, mfeff, mchi, ms,
     1                  ntitle, title, ltit,
     2                  critcw,
     1                  ipr2, ipr3, ipr4,
     1                  s02, tk, thetad, sig2g,
     1                  nlegxx,
     1                  rmax, critpw, pcritk, pcrith, nncrit,
     2                  icsig, iorder, vrcorr, vicorr, isporb)

c     Read input for multiple scattering feff
      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola


c     Following passed to pathfinder, which is single precision.
c     Be careful to always declare these!
      real*8 rmax, critpw, pcritk, pcrith

c     Data for potph (see arrays.h for comments)
      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension ifrph(0:nphx)
      dimension xnatph(0:nphx)
      dimension folp(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension ion(0:nfrx)
      dimension iz(0:nfrx)

      character*6  potlbl(0:nphx)

c     Local stuff
      character*150  line
      parameter (nwordx = 12)
      character*15 words(nwordx)

      parameter (ntitx = 10)
      character*79  title(ntitx)
      dimension ltit(ntitx)
      dimension ionph(0:nphx), izph(0:nphx)
      logical iscomm
      parameter (nssx = 16)
      dimension indss(nssx), iphss(nssx)
      dimension degss(nssx), rss(nssx)
      logical nogeom

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

c     initialize things

      ihole = 1
      ntitle = 0
      ixc = 0
      vr0 = 0
      vi0 = 0
      rs0 = 0
      rmax = -1
      tk = 0
      thetad = 0
      sig2g = 0
      rmult = 1
      s02 = 1
      mphase = 1
      mpath = 1
      mfeff = 1
      mchi = 1
      ms = 0
      ipr1 = 0
      ipr2 = 0
      ipr3 = 0
      ipr4 = 0
      nlegxx = 10
      xkmin = 0
      xkmax = 20
      critcw = 4.0d0
      critpw = 2.5d0
      pcritk = 0
      pcrith = 0
      nogeom = .false.
      icsig = 1
      iorder = 2
      ixanes = 0
      vrcorr = 0
      vicorr = 0
      iafolp = 0
      intclc = 0
      nemax = nex
      isporb = -1

c     average over polarization by default
      pola = .false.
      elpty = 0
      do 50 i = 1, 3 
         evec(i) = 0
         ivec(i) = 0
  50  continue 

c     nncrit is number of necrit points to use.  necrit is
c     currently 9, this was at once an input used for testing.
      nncrit = 9

      nat = 0
      do 100  iat = 1, natx
         iphat(iat) = -1
  100 continue

      nss = 0
      do 102  iss = 1, nssx
         indss(iss) = 0
         iphss(iss) = 0
         degss(iss) = 0
         rss(iss) = 0
  102 continue

      nph = 0
      do 110  iph = 0, nphx
         iatph(iph) = 0
         ifrph(iph) = -1
         xnatph(iph) = 0
         folp(iph) = 1
         novr(iph) = 0
         ionph(iph) = 0
         izph(iph) = 0
         potlbl(iph) = ' '
  110 continue

      nfr = 0
      do 120  ifr = 0, nfrx
         ion(ifr) = 0
         iz(ifr) = 0
  120 continue

c     Open feff.inp, the input file we're going to read
      open (unit=1, file=trim(header)//'feff.inp',
     >      status='old', iostat=ios)
      call chopen (ios, trim(header)//'feff.inp', 'rdinp')

c     tokens  0 if not a token
c             1 if ATOM (ATOMS)
c             2 if HOLE
c             3 if OVER (OVERLAP)
c             4 if CONT (CONTROL)
c             5 if EXCH (EXCHANGE)
c             6 if ION
c             7 if TITL (TITLE)
c             8 if FOLP
c             9 if RMAX
c            10 if DEBY (DEBYE)
c            11 if RMUL (RMULTIPLIER)
c            12 if SS
c            13 if PRIN (PRINT)
c            14 if POTE (POTENTIALS)
c            15 if NLEG
c            16 if REQU (REQUIRE), now dead
c            17 if KLIM (KLIMIT)
c            18 if CRIT (CRITERIA)
c            19 if NOGEOM
c            20 if CSIG
c            21 if IORDER
c            22 if PCRI (PCRITERIA)
c            23 if SIG2
c            24 if XANE (XANES), disabled for current release
c            25 if CORR (CORRECTIONS)
c            26 if AFOL (AFOLP)
c            27 if NEMA (NEMAX)
c            28 if INTCALC
c            29 if POLA (POLARIZATION)
c            30 if ELLI (ELLIPTICITY) 
c            31 if ISPO (ISPORB)
c            -1 if END  (end)
c     mode flag  0 ready to read a keyword card
c                1 reading atom positions
c                2 reading overlap instructions for unique pot
c                3 reading unique potential definitions

      mode = 0
  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  line='END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1))

c        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
c              ATOM
c              Following lines are atom postions, one per line
               mode = 1
            elseif (itok .eq. 2)  then
c              HOLE     1  1.0
c                   holecode s02
               read(words(2),20,err=900)  ihole
               read(words(3),30,err=900)  s02
               mode = 0
            elseif (itok .eq. 3)  then
c              OVERLAP iph
c                  iph  n  r
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               mode = 2
            elseif (itok .eq. 4)  then
c              CONTROL  mphase, mpath, mfeff, mchi
c               0 - do not run modules, 1 - run module
               read(words(2),20,err=900)  mphase
               read(words(3),20,err=900)  mpath
               read(words(4),20,err=900)  mfeff
               read(words(5),20,err=900)  mchi
               mode = 0
            elseif (itok .eq. 5)  then
c              EXCHANGE  ixc  vr0  vi0
c              ixc=0  Hedin-Lunqvist + const real & imag part
c              ixc=1  Dirac-Hara + const real & imag part
c              ixc=2  ground state + const real & imag part
c              ixc=3  Dirac-Hara + HL imag part + const real & imag part
c              ixc=4  DH below rs0 + HL above rs0 + const real
c                     & imag part, form is
c                     EXCHANGE  4  vr0  vi0  rs0
c              vr0 is const imag part of potential
c              vi0 is const imag part of potential
c              Default is HL.
               read(words(2),20,err=900)  ixc
               read(words(3),30,err=900)  vr0
               read(words(4),30,err=900)  vi0
               if (ixc .eq. 4) read(words(5),30,err=900)  rs0
               if (ixc .ge. 3)  call warnex(1)
               mode = 0
            elseif (itok .eq. 6)  then
c              ION  iph ionph(iph)
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),20,err=900)  ionph(iph)
               mode = 0
            elseif (itok .eq. 7)  then
c              TITLE title...
               ntitle = ntitle + 1
               if (ntitle .le. ntitx)  then
                  title(ntitle) = line(6:)
                  call triml (title(ntitle))
               else
                  write(77,*) 'Too many title lines, title ignored'
                  write(77,*) line(1:79)
               endif
               mode = 0
            elseif (itok .eq. 8)  then
c              FOLP iph folp (overlap factor, default 1)
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),30,err=900)  folp(iph)
               mode = 0
            elseif (itok .eq. 9)  then
c              RMAX  rmax (max r for ss and pathfinder)
               read(words(2),30,err=900)  rmax
               mode = 0
            elseif (itok .eq. 10)  then
c              DEBYE  temp debye-temp
c                   temps in kelvin
c                   if tk and thetad > 0, use these instead of sig2g
               read(words(2),30,err=900)  tk
               read(words(3),30,err=900)  thetad
               mode = 0
            elseif (itok .eq. 11)  then
c              RMULTIPLIER  rmult
c              Multiples atom coord, rss, overlap and rmax distances by
c              rmult (default 1).  DOES NOT modify sig2g
               read(words(2),30,err=900)  rmult
               mode = 0
            elseif (itok .eq. 12)  then
c              SS index ipot deg rss
               nss = nss + 1
               if (nss .gt. nssx)  then
                  write(77,*)
     >             'Too many ss paths requested, max is ', nssx
                  stop 'RDINP'
               endif
               read(words(2),20,err=900)  indss(nss)
               read(words(3),20,err=900)  iphss(nss)
               read(words(4),30,err=900)  degss(nss)
               read(words(5),30,err=900)  rss(nss)
               mode = 0
            elseif (itok .eq. 13)  then
c              PRINT  ipr1  ipr2  ipr3  ipr4
c              print flags for various modules
c              ipr1 potph  0 phase.bin only
c                          1 add misc.dat
c                          2 add pot.dat, phase.dat
c                          5 add atom.dat
c                          6 add central atom dirac stuff
c                          7 stop after doing central atom dirac stuff
c              ipr2 pathfinder  0 paths.dat only
c                               1 add crit.dat
c                               2 keep geom.dat
c                               3 add fbeta files
c                               5 special magic code, crit&geom only
c                                 not paths.dat.  Use for path studies
c              ipr3 genfmt 0 files.dat, feff.dats that pass 2/3 of
c                            curved wave importance ratio
c                          1 keep all feff.dats
c              ipr4 ff2chi 0 chi.dat
c                          1 add sig2.dat with debye waller factors
c                          2 add chipnnnn.dat for each path
               read(words(2),20,err=900)  ipr1
               read(words(3),20,err=900)  ipr2
               read(words(4),20,err=900)  ipr3
               read(words(5),20,err=900)  ipr4
               mode = 0
            elseif (itok .eq. 14)  then
c              POTENTIALS
c              Following lines are unique potential defs, 1 per line
               mode = 3
            elseif (itok .eq. 15)  then
c              NLEG nlegmax (for pathfinder)
               read(words(2),20,err=900)  nlegxx
               mode = 0
            elseif (itok .eq. 16)  then
c              REQUIRE rreq, ipot (for pathfinder, require than ms paths
c                            length >rreq contain atom ipot)
               write(77,*) 'REQUIRE no longer available'
               stop
            elseif (itok .eq. 17)  then
c              KLIMIT xkmin, xkmax
               write(77,*) 'KLIMIT no longer available, run continues.'
               mode = 0
            elseif (itok .eq. 18)  then
c              CRIT critcw critpw
               read(words(2),30,err=900)  critcw
               read(words(3),30,err=900)  critpw
               mode = 0
            elseif (itok .eq. 19)  then
c              NOGEOM (do not write geom.dat)
               nogeom = .true.
               mode = 0
            elseif (itok .eq. 20)  then
c              CSIG (use complex momentum with debye waller factor)
c              note: this is always on anyway, so this card unnecessary
               icsig = 1
               mode = 0
            elseif (itok .eq. 21)  then
c              IORDER  iorder (used in genfmt, see setlam for meaning)
               read(words(2),20,err=900)  iorder
               call warnex(2)
               mode = 0
            elseif (itok .eq. 22)  then
c              PCRIT  pcritk pcrith
c                     (keep and heap criteria for pathfinder)
               read(words(2),30,err=900)  pcritk
               read(words(3),30,err=900)  pcrith
               mode = 0
            elseif (itok .eq. 23)  then
c              SIG2  sig2g   global sig2 written to files.dat
               read(words(2),30,err=900)  sig2g
               mode = 0
            elseif (itok .eq. 24)  then
c              XANES
c              Use extended k range for xanes
               ixanes = 1
c              to avoid problems with debye waller factors below the
c              edge, always use complex p for debye waller
               icsig = 1
               call warnex(3)
               write(77,212)
  212          format ( ' CORRECTIONS and other cards may be needed.',
     1            '  See FEFF6 document for', /,
     2            ' details and a discussion of approximations.')
               mode = 0
            elseif (itok .eq. 25)  then
c              CORRECTIONS  e0-shift, lambda correction
c              e0 shift is in eV, edge will be edge-e0
c              lambda corr is a const imag energy in eV
c              e0 and lambda corr same as vr0 and vi0 in EXCH card
               read(words(2),30,err=900)  vrcorr
               read(words(3),30,err=900)  vicorr
               mode = 0
            elseif (itok .eq. 26)  then
c              AFOLP use generalized automatic folp
               iafolp = 1
               mode =0
            elseif (itok .eq. 27)  then
c              NEMAX  nemax for energy grid
               read(words(2),20,err=900)  nemax
               call warnex(4)
               if (nemax .gt. nex)  then
                  write(77,*) 'nemax too big, nemax, nex, ', nemax, nex
                  nemax = nex
                  write(77,*) 'nemax reset to ', nemax
               endif
               mode = 0
            elseif (itok .eq. 28)  then
c              INTCALC  intclc
c              0  use average over all atoms
c              1  use current experimental method 1
c              2  use current experimental method 2
c              read(words(2),20,err=900)  intclc
               write(77,*) 'INTCALC not implemented -- card ignored.'
               mode = 0
            elseif (itok .eq. 29)  then
c              POLARIZATION  X Y Z
               pola = .true.
c              run polarization code if 'pola' is true
c              run usual feff otherwise
               read(words(2),30,err=900)  evec(1)
               read(words(3),30,err=900)  evec(2)
               read(words(4),30,err=900)  evec(3)
               mode = 0
            elseif (itok .eq. 30)  then
c              ELLIPTICITY  E incident direction
               read(words(2),30,err=900)  elpty
               read(words(3),30,err=900)  ivec(1)
               read(words(4),30,err=900)  ivec(2)
               read(words(5),30,err=900)  ivec(3)
               mode = 0
            elseif (itok .eq. 31)  then
c              ISPORB  isporb
               read(words(2),20,err=900)  isporb
               write(77,*) ' isporb set ', isporb
               mode = 0
            elseif (itok .eq. -1)  then
c              END
               goto 220
            else
               write(77,*) line(1:70)
               write(77,*) words(1)
               write(77,*) 'Token ', itok
               write(77,*) 'Keyword unrecognized.'
               write(77,*) 'See FEFF document -- some old features'
               write(77,*) 'are no longer available.'
               stop 'RDINP-2'
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
c              We're done reading atoms.
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            nat = nat+1
            if (nat .gt. natx)  then
               write(77,*) 'Too many atoms, max is ', natx
               stop 'RDINP-3'
            endif
            read(words(1),30,err=900)  rat(1,nat)
            read(words(2),30,err=900)  rat(2,nat)
            read(words(3),30,err=900)  rat(3,nat)
            read(words(4),20,err=900)  iphat(nat)
         elseif (mode .eq. 2)  then
            if (itok .ne. 0)  then
c              We're done reading these overlap instructions.
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            novr(iph) = novr(iph)+1
            iovr = novr(iph)
            if (iovr .gt. novrx)  then
               write(77,*) 'Too many overlap shells, max is ', novrx
               stop 'RDINP-5'
            endif
            read(words(1),20,err=900) iphovr(iovr,iph)
            read(words(2),20,err=900) nnovr(iovr,iph)
            read(words(3),30,err=900) rovr(iovr,iph)
         elseif (mode .eq. 3)  then
            if (itok .ne. 0)  then
c              We're done reading unique potential definitions
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900)  iph
            if (iph .lt. 0  .or.  iph .gt. nphx)  then
               write(77,*) 'Unique potentials must be between 0 and ',
     1                 nphx
               write(77,*) iph, ' not allowed'
               write(77,*) line(1:79)
               stop 'RDINP'
            endif
            read(words(2),20,err=900)  izph(iph)
c           No potential label if user didn't give us one
c           Default set above is potlbl=' '
            if (nwords .ge. 3)  potlbl(iph) = words(3)
         else
            write(77,*) 'Mode unrecognized, mode ', mode
            stop 'RDINP-6'
         endif
      goto 200
  220 continue

c     We're done reading the input file, close it.
      close (unit=1)

c     Fix up defaults, error check limits, figure out free atoms, etc.

      if (pola) then
c        make polarization tensor
         call mkptz
      endif

c     Find out how many unique potentials we have
      nph = 0
      do 300  iph = nphx, 0, -1
         if (izph(iph) .gt. 0)  then
            nph = iph
            goto 301
         endif
  300 continue
  301 continue
c     Must have central atom
      if (izph(0) .le. 0)  then
       write(77,*) 'Absorbing atom, unique potential 0, is not defined.'
       stop 'RDINP'
      endif

c     Then find model atoms for unique pots that have them
      do 330  iph = 0, nphx
c        Use first atom in atom list that is of unique pot iph
         do 320  iat = 1, nat
            if (iph .eq. iphat(iat))  then
               iatph(iph) = iat
               goto 321
            endif
  320    continue
  321    continue
  330 continue
c     if iatph > 0, a model atom has been found.

c     No gaps allowed in unique pots.  Make sure we have enough
c     to overlap all unique pots 0 to nph.
      do 340  iph = 0, nph
         if (iatph(iph) .le. 0  .and.  novr(iph) .le. 0)  then
c           No model atom, no overlap cards, can't do this unique pot
            write(77,*) ' No atoms or overlap cards for unique pot ', iph
            write(77,*) ' Cannot calculate potentials, etc.'
            stop 'RDINP-'
         endif
  340 continue

c     Need number of atoms of each unique pot, count them.  If none,
c     set to one.
      do 350  iph = 0, nph
         xnatph(iph) = 0
         do 346  iat = 1, nat
            if (iphat(iat) .eq. iph)  xnatph(iph) = xnatph(iph)+1
  346    continue
         if (xnatph(iph) .le. 0)  xnatph(iph) = 1
  350 continue

c     Do the free atom shuffling, do central atom as special case
      iz(0) = izph(0)
      ion(0) = ionph(0)
      ifrph(0) = 0
      nfr = 0
      do 390  iph = 1, nph
         ifrph(iph) = -1
         do 380  ifr = 1, nfr
            if (iz(ifr).eq.izph(iph) .and. ion(ifr).eq.ionph(iph)) then
               ifrph(iph) = ifr
               goto 381
            endif
  380    continue
  381    continue
c        add free atom type if necessary
         if (ifrph(iph) .lt. 0)  then
            nfr = nfr+1
            if (nfr .gt. nfrx)  then
               write(77,*) ' Too many free atoms, max is ', nfrx
               stop 'RDINP10'
            endif
            ion(nfr) = ionph(iph)
            iz(nfr) = izph(iph)
            ifrph(iph) = nfr
         endif
  390 continue

c     Find central atom (only 1 permitted)
      iatabs = -1
      do 400  iat = 1, nat
         if (iphat(iat) .eq. 0)  then
            if (iatabs .lt. 0)  then
               iatabs = iat
            else
               write(77,*) 'More than one absorbing atom (potential 0)'
               write(77,*) 'Only one absorbing atom allowed'
               stop 'RDINP'
            endif
         endif
  400 continue

c     Find distance to nearest and most distant atom (use overlap card
c     if no atoms specified.)
      if (iatabs .lt. 0  .or.  nat .lt. 2)  then
         ratmin = rovr(1,0)
         ratmax = rovr(novr(0),0)
      else
         ratmax = 0
         ratmin = 1.0d10
         do 412  iat = 1, nat
c           skip absorbing atom
            if (iat .eq. iatabs)  goto 412
            tmp = feff_dist(rat(1,iat), rat(1,iatabs))
            if (tmp .gt. ratmax)  ratmax = tmp
            if (tmp .lt. ratmin)  ratmin = tmp
  412    continue
      endif

c     Set rmax if necessary
      if (rmax.le.0 .and. nss.le.0)  then
c        set to min (2+ times ratmin, ratmax)
         rmax = min (2.001 * ratmin, ratmax)
      endif

c     Set core hole lifetime (central atom quantity)
      ifr = ifrph(0)
      call setgam (iz(ifr), ihole, gamach)

c     Set s02 if necessary
      if (s02 .le. 1.0d-10)  s02 = 1

c     Convert everything to code units, and use rmult factor
c     rmax is for pathfinder, so leave it in Ang.
      rmax = rmax * rmult
      vr0 = vr0 / ryd
      vi0 = vi0 / ryd
      vrcorr = vrcorr / ryd
      vicorr = vicorr / ryd
      xkmin = xkmin * bohr
      xkmax = xkmax * bohr
      do 430  iat = 1, nat
         do 420  i = 1, 3
            rat(i,iat) = rat(i,iat) * rmult / bohr
  420    continue
  430 continue
      do 460  iph = 0, nph
         do 450  iovr = 1, novr(iph)
            rovr(iovr,iph) = rovr(iovr,iph) * rmult / bohr
  450    continue
  460 continue
      do 462  iss = 1, nss
c        rss used only to make paths.dat, so leave it in Angstroms.
         rss(iss) = rss(iss) * rmult
  462 continue

c     Check if 2 atoms are closer together than 1.75 ryd (~.93 Ang)
      ratmin = 1.0d20
      do 480  iat = 1, nat
         do 470  jat = iat+1, nat
            rtmp = feff_dist(rat(1,iat),rat(1,jat))
            if (rtmp .lt. ratmin)  ratmin = rtmp
            if (rtmp .lt. 1.75d0)  then
c           if (dist(rat(1,iat),rat(1,jat)) .lt. 1.5)  then
               write(77,*) 'WARNING:  TWO ATOMS VERY CLOSE TOGETHER.',
     1                 '  CHECK INPUT.'
               write(77,*) ' atoms ', iat, jat
               write(77,*) iat, (rat(i,iat)*bohr,i=1,3)
               write(77,*) jat, (rat(i,jat)*bohr,i=1,3)
               write(77,*) 'Run continues in case you really meant it.'
            endif
  470    continue
  480 continue

c     default to k shell
      if (isporb .lt. 0)  isporb = 1

c     Clean up control flags
      if (mphase .ne. 0)  mphase = 1
      if (mpath  .ne. 0)  mpath = 1
      if (mfeff  .ne. 0)  mfeff = 1
      if (mchi   .ne. 0)  mchi = 1
      if (nss    .le. 0)  ms = 1

      if (ntitle .le. 0)  then
         ntitle = 1
         title(i) = 'No title input'
      endif
      do 490  i = 1, ntitle
         ltit(i) = istrln (title(i))
  490 continue

c     Write output files

c     For potph...
      if (mphase .eq. 1)  then
         open (unit=1, file=trim(header)//'potph.dat',
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//'potph.dat', 'rdinp')
         do 705  i = 1, ntitle
            write(1,700)  title(i)(1:ltit(i))
  700       format (1x, a)
  705    continue
         write(1,706)
  706    format (1x, 79('-'))
         write(1,709) ihole, gamach, ipr1, iafolp, intclc
  709    format(i5, 1p, e14.6, 3i4, 
     1         ' ihole, gamach, iprint, iafolp, intclc')
         write(1,702)  ixc, vr0, vi0, rs0
  702    format (i5, 1p, 3e14.6, ' ixc, vr0, vi0, rs0')
         write(1,701)  ixanes, nemax, xkmin, xkmax
  701    format (2i5, 1p, 2e14.6, 
     1           ' ixanes, nemax, xkmin, xkmax (inv bohr)')
         write(1,707) nfr, '  nfr'
  707    format (i5, a)
         do 710  ifr = 0, nfr
            write(1,708)  ifr, iz(ifr), ion(ifr)
  708       format (3i5, ' ifr, iz, ion')
  710    continue
         write(1,707) nat, '  nat.   iat, iph, x, y, z'
         do 720  iat = 1, nat
            write(1,715) iat, iphat(iat), (rat(j,iat),j=1,3)
  715       format (2i5, 3f12.6)
  720    continue
         write(1,707) nph, '  nph'
         do 740  iph = 0, nph
            write(1,722) iph, iatph(iph), ifrph(iph), xnatph(iph),
     1                   folp(iph), novr(iph),
     2                   ' iph, iat, ifr, xnat, folp, novr'
  722       format (3i5, 2f12.6, i5, a)
            write(1,723) potlbl(iph)
  723       format (' ''', a6, '''  potlbl')
            do 730  iovr = 1, novr(iph)
               write(1,724) iphovr(iovr,iph), nnovr(iovr,iph),
     1                      rovr(iovr,iph),
     2                      ' ovr...  iph, n, r'
  724       format (2i5, f12.6, a)
  730       continue
  740    continue
         close (unit=1)
      endif

c     Single scattering paths for genfmt
      if (nss .gt. 0  .and.  mpath .eq. 1)  then
         open (unit=1, file=trim(header)//'paths.dat',
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//'paths.dat', 'rdinp')
         do 750  i = 1, ntitle
            write(1,748)  title(i)(1:ltit(i))
  748       format (1x, a)
  750    continue
         write(1,751)
  751    format (' Single scattering paths from ss lines cards',
     1           ' in feff input')
         write(1,706)
         do 760  iss = 1, nss
            if (rmax.le.0  .or.  rss(iss).le.rmax)  then
c              NB, rmax and rss are in angstroms
               write(1,752) indss(iss), 2, degss(iss),
     2              rss(iss)
  752          format ( 2i4, f8.3,
     1             '  index,nleg,degeneracy,r=', f8.4)
               write(1,766)
  766          format (' single scattering')
               write(1,754) rss(iss)*bohr, zero, zero, iphss(iss),
     1                      potlbl(iphss(iss))
               write(1,753) zero, zero, zero, 0, potlbl(0)
  753          format (3f12.6, i4,  1x, '''', a6, '''', '  x,y,z,ipot')
  754          format (3f12.6, i4,  1x, '''', a6, '''')
            endif
  760    continue
         close (unit=1)
      endif

c     Atoms for the pathfinder
      if (nss.le.0  .and.  mpath.eq.1  .and.  nat.gt.0)  then
         if (iatabs .le. 0)  then
            write(77,*) 'Absorbing atom coords not specified.'
            write(77,*) 'Cannot find multiple scattering paths.'
            stop 'RDINP'
         endif
c        if user doesn't want geom.dat, don't do it
         if (nogeom)  goto 792
         open (unit=1, file=trim(header)//'geom.dat',
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//'geom.dat', 'rdinp')
c        Echo title cards to geom.dat
         do 770  i = 1, ntitle
            write(1,700)  title(i)(1:ltit(i))
  770    continue
         write(1,706)
c        Central atom first
         ii = 0
         write(1,780)  ii, (rat(j,iatabs)*bohr,j=1,3), 0, 1
c        Rest of the atoms (skip central atom)
         do 790   iat = 1, nat
            if (iat .eq. iatabs)  goto 790
            ii = ii+1
            write(1,780)  ii, (rat(j,iat)*bohr,j=1,3), iphat(iat), 1
  780       format (i4, 3f12.6, 2i4)
  790    continue
         close (unit=1)
      endif
  792 continue

      return

  900 continue
      write(77,*) 'Error reading input, bad line follows:'
      write(77,*) line(1:79)
      stop 'RDINP fatal error.'

      end

      function itoken (word)
      implicit double precision (a-h, o-z)
c     chars in word assumed upper case, left justified
c     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w

      w = word(1:4)
      if     (w .eq. 'ATOM')  then
         itoken = 1
      elseif (w .eq. 'HOLE')  then
         itoken = 2
      elseif (w .eq. 'OVER')  then
         itoken = 3
      elseif (w .eq. 'CONT')  then
         itoken = 4
      elseif (w .eq. 'EXCH')  then
         itoken = 5
      elseif (w .eq. 'ION ')  then
         itoken = 6
      elseif (w .eq. 'TITL')  then
         itoken = 7
      elseif (w .eq. 'FOLP')  then
         itoken = 8
      elseif (w .eq. 'RMAX')  then
         itoken = 9
      elseif (w .eq. 'DEBY')  then
         itoken = 10
      elseif (w .eq. 'RMUL')  then
         itoken = 11
      elseif (w .eq. 'SS  ')  then
         itoken = 12
      elseif (w .eq. 'PRIN')  then
         itoken = 13
      elseif (w .eq. 'POTE')  then
         itoken = 14
      elseif (w .eq. 'NLEG')  then
         itoken = 15
      elseif (w .eq. 'REQU')  then
         itoken = 16
      elseif (w .eq. 'KLIM')  then
         itoken = 17
      elseif (w .eq. 'CRIT')  then
         itoken = 18
      elseif (w .eq. 'NOGE')  then
         itoken = 19
      elseif (w .eq. 'CSIG')  then
         itoken = 20
      elseif (w .eq. 'IORD')  then
         itoken = 21
      elseif (w .eq. 'PCRI')  then
         itoken = 22
      elseif (w .eq. 'SIG2')  then
         itoken = 23
      elseif (w .eq. 'XANE')  then
         itoken = 24
      elseif (w .eq. 'CORR')  then
         itoken = 25
      elseif (w .eq. 'AFOL')  then
         itoken = 26
      elseif (w .eq. 'NEMA')  then
         itoken = 27
      elseif (w .eq. 'INTC')  then
         itoken = 28
      elseif (w .eq. 'POLA')  then
         itoken = 29
      elseif (w .eq. 'ELLI')  then
         itoken = 30
      elseif (w .eq. 'ISPO')  then
         itoken = 31
      elseif (w .eq. 'END ')  then
         itoken = -1
      else
         itoken = 0
      endif
      return
      end
      logical function iscomm (line)
      implicit double precision (a-h, o-z)
c     returns true if line is a comment or blank line, false otherwise
      character*(*) line
      iscomm = .false.
      if (istrln(line).le.0  .or.  line(1:1).eq.'*')  iscomm = .true.
      return
      end
      subroutine phstop (iph,line)
      implicit double precision (a-h, o-z)
      character*(*) line

      character*72 header
      common /header_common/ header


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      if (iph .lt. 0  .or.  iph .gt. nphx)  then
         write(77,10) iph, nphx, line
   10    format (' Unique potential index', i5, ' out of range.', /,
     1           ' Must be between 0 and', i5, '.  Input line:', /,
     2           1x, a)
         stop 'RDINP - PHSTOP'
      endif
      return
      end
      subroutine warnex (i)
      implicit double precision (a-h, o-z)
c     This prints a warning message if the user is using an
c     expert option.
c     i    expert option card
c     1    EXCHANGE with code >= 3
c     2    IORDER
c     3    XANES
c     4    NEMAX
c     5    INTCALC

c     message max of 22 characters to keep warning on 80 char line.
  100 format (1x, a, 
     1   ': Expert user option, please read documentation', /,
     2   ' carefully and check your results.')

      if (i .eq. 1)  then
         write(77,100) 'EXCHANGE code >= 3'
      elseif (i .eq. 2)  then
         write(77,100) 'IORDER'
      elseif (i .eq. 3)  then
         write(77,100) 'XANES'
      elseif (i .eq. 4)  then
         write(77,100) 'NEMAX'
      elseif (i .eq. 5)  then
         write(77,100) 'INTCALC'
      endif
      return
      end
      subroutine rdpath (in, pol, done,xstar)
      implicit double precision (a-h, o-z)
      logical done, pol


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines


c     global polarization data
      logical  pola
      double precision evec,ivec,elpty
      complex*16 ptz
      common /pol/ evec(3), ivec(3), elpty, ptz(-1:1,-1:1), pola


      complex*16  alph, gamm
      dimension  alpha(0:legtot), gamma(legtot)

      read(in,*,end=200)  ipath, nleg, deg
      if (nleg .gt. legtot)  then
         write(77,*) 'nleg .gt. legtot, nleg, legtot ', nleg, legtot
         write(77,*) 'ERROR'
         goto 200
      endif
c     skip label (x y z ipot rleg beta eta)
      read(in,*)
      do 20  ileg = 1, nleg
         read(in,*,end=999)  (rat(j,ileg),j=1,3), ipot(ileg), 
     1                       potlbl(ipot(ileg))
c        convert to code units
         do 10  j = 1, 3
            rat(j,ileg) = rat(j,ileg)/bohr
   10    continue
         if (ipot(ileg) .gt. npot)  then
            write(77,*) 'ipot(ileg) too big, ipot, ileg, npot ',
     1               ipot(ileg), ileg, npot
            write(77,*) 'ERROR'
            goto 200
         endif
   20 continue
      nsc = nleg-1

c     We need the 'z' atom so we can use it below.  Put
c     it in rat(nleg+1).  No physical significance, just a handy
c     place to put it.
      if (pol) then
         rat(1,nleg+1) = rat(1,nleg)
         rat(2,nleg+1) = rat(2,nleg)
         rat(3,nleg+1) = rat(3,nleg) + 1.0d0
      endif

c     add rat(0) and ipot(0) (makes writing output easier)
      do 22 j = 1, 3
         rat(j,0) = rat(j,nleg)
   22 continue
      ipot(0) = ipot(nleg)

c     beginnnig of calculating nstar=deg*cos(eps r1)*cos(eps rN)
      x1 = 0.0
      do 23 j = 1,3
         x1 = x1 + evec(j) * ( rat(j,1) - rat(j,0) )
   23 continue
      xnorm = 0.0
      do 24 j = 1,3
         xnorm = xnorm + (rat(j,1) - rat(j,0))**2
   24 continue
      x1 = x1/sqrt(xnorm)
      x2 = 0.0
      do 25 j = 1,3
         x2 = x2 + evec(j) * ( rat(j,nleg-1) - rat(j,0) )
   25 continue
      xnorm = 0.0
      do 26 j = 1,3
         xnorm = xnorm + (rat(j,nleg-1) - rat(j,0))**2
   26 continue
      x2 = x2/sqrt(xnorm)
      xstar = deg* abs(x1*x2)
c     end of calculating nstar

      nangle = nleg
      if (pol) then 
c        in polarization case we need one more rotation
         nangle = nleg + 1
      endif
      do 100  j = 1, nangle

c        for euler angles at point i, need th and ph (theta and phi)
c        from rat(i+1)-rat(i)  and  thp and php
c        (theta prime and phi prime) from rat(i)-rat(i-1)
c
c        Actually, we need cos(th), sin(th), cos(phi), sin(phi) and
c        also for angles prime.  Call these  ct,  st,  cp,  sp

c        i = (j)
c        ip1 = (j+1)
c        im1 = (j-1)
c        except for special cases...
         ifix = 0
         if (j .eq. nsc+1)  then
c           j+1 'z' atom, j central atom, j-1 last path atom
            i = 0
            ip1 = 1
            if (pol) then
               ip1 = nleg+1
            endif
            im1 = nsc

         elseif (j .eq. nsc+2)  then
c           j central atom, j+1 first path atom, j-1 'z' atom
            i = 0
            ip1 = 1
            im1 = nleg+1
            ifix = 1
         else
            i = j
            ip1 = j+1
            im1 = j-1
         endif

         x = rat(1,ip1) - rat(1,i)
         y = rat(2,ip1) - rat(2,i)
         z = rat(3,ip1) - rat(3,i)
         call trig (x, y, z, ctp, stp, cpp, spp)
         x = rat(1,i) - rat(1,im1)
         y = rat(2,i) - rat(2,im1)
         z = rat(3,i) - rat(3,im1)
         call trig (x, y, z, ct, st, cp, sp)

c        Handle special case, j=central atom, j+1 first
c        path atom, j-1 is 'z' atom.  Need minus sign
c        for location of 'z' atom to get signs right.
         if (ifix .eq. 1)  then
            x = 0
            y = 0
            z = 1.0
            call trig (x, y, z, ct, st, cp, sp)
            ifix = 0
         endif

c        cppp = cos (phi prime - phi)
c        sppp = sin (phi prime - phi)
         cppp = cp*cpp + sp*spp
         sppp = spp*cp - cpp*sp
         phi  = atan2(sp,cp)
         phip = atan2(spp,cpp)

c        alph = exp(i alpha)  in ref eqs 18
c        beta = cos (beta)         
c        gamm = exp(i gamma)
         alph = -(st*ctp - ct*stp*cppp - coni*stp*sppp)
         beta(j) = ct*ctp + st*stp*cppp
c        watch out for roundoff errors
         if (beta(j) .lt. -1) beta(j) = -1
         if (beta(j) .gt.  1) beta(j) =  1
         gamm = -(st*ctp*cppp - ct*stp + coni*st*sppp)
         call feff_arg(alph,phip-phi,alpha(j))
         beta(j) = acos(beta(j))
         call feff_arg(gamm,phi-phi,gamma(j))
c       Convert from the rotation of FRAME used before to the rotation 
c       of VECTORS used in ref.
         dumm = alpha(j)
         alpha(j) =  pi- gamma(j)
         gamma(j) =  pi- dumm

         if (j .le. nleg)  then
            ri(j) = feff_dist(rat(1,i), rat(1,im1))
         endif
  100 continue

c     Make eta(i) = alpha(i-1) + gamma(i). 
c     We'll need alph(nangle)=alph(0)
      alpha(0) = alpha(nangle)
      do 150  j = 1, nleg
         eta(j) = alpha(j-1) + gamma(j)
  150 continue
      if (pol) then
         eta(0) = gamma(nleg+1)
         eta(nleg+1) = alpha(nleg)
      endif

c     eta and beta in radians at this point.
      done = .false.
      return

c     If no more data, tell genfmt we're done
  200 continue
      done = .true.
      return

c     If unexpected end of file, die
  999 continue
      write(77,*) 'Unexpected end of file'
      stop 'ERROR'
      end
      subroutine trig (x, y, z, ct, st, cp, sp)
      implicit double precision (a-h, o-z)
c     returns cos(theta), sin(theta), cos(phi), sin(ph) for (x,y,z)
c     convention - if x=y=0 and z>0, phi=0, cp=1, sp=0
c                  if x=y=0 and z<0, phi=180, cp=-1,sp=0
c                - if x=y=z=0, theta=0, ct=1, st=0
      parameter (eps = 1.0d-6)
      r = sqrt (x**2 + y**2 + z**2)
      rxy = sqrt (x**2 + y**2)
      if (r .lt. eps)  then
         ct = 1
         st = 0
      else
         ct = z/r
         st = rxy/r
      endif
      if (rxy .lt. eps)  then
         cp = 1
         if (ct .lt. 0) cp = -1
         sp = 0
      else
         cp = x / rxy
         sp = y / rxy
      endif
      return
      end
      subroutine feff_arg(c,fi,th)
      implicit double precision (a-h, o-z)
      complex*16  c
      parameter (eps = 1.0d-6)
      x = dble(c)
      y = dimag(c)
      if (abs(x) .lt. eps) x = 0
      if (abs(y) .lt. eps) y = 0
      if (abs(x) .lt. eps  .and.  abs(y) .lt. eps) then
        th = fi
      else
        th = atan2(y,x)
      endif
      return
      end
      subroutine renorm (dexv, vcoul, srho)

      implicit double precision (a-h,o-z)
      save

      common /print/ iprint
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco
      integer*4 nstop
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)

c     vcoul is the coulomb potential (no factor of r**2) (output)
      dimension vcoul(251)
c     srho is charge density in form 4*pi*density*r**2 output)
      dimension srho(251)
c jm  9/23/87 added srho renormalized charge density to be used
c     in cphase

      do 10 i=1,np
         dv(i)=0.0
         d(i)=0.0
   10 continue
      ddjri=log(ws/dr(1))/dpas
      jri=1.0+ddjri
      jr1=jri
      ddjr1=ddjri-jr1+1.0

      if (jri-2*(jri/2).ne.0) go to 20
         jri=jri+1
   20 continue

      ddjri=ddjri-jri+1.0
c  ddjri = (log(ws)-dri)/dpas
c  dri  =  log(dr(jri))

      da=0.0
      do 30 j=1,norb
      do 30 i=1,np
   30    d(i)=d(i)+nel(j)*(dgc(i,j)**2+dpc(i,j)**2)

      do 50 i=jri,np
         dl=dr(i)
         if (i.eq.jri.or.i.eq.np) go to 40
            dl=dl+dl
            if ((i-2*(i/2)).eq.0) dl=dl+dl
   40    dd=d(i)*dl
         da=da+dd
   50 continue

      da=dpas*da/3.0
      dfo=dr(jri-1)*d(jri-1)
      df1=dr(jri)*d(jri)
      df2=dr(jri+1)*d(jri+1)
      dcor=-dpas*(df1*ddjri+(df2+dfo-2.0*df1)*ddjri**3/6.0+(df2-dfo)
     1 *ddjri**2*.25)
      da=da+dcor
      if (iprint .ge. 5)  write(16,60) da
   60 format (1h ,' no. of electrons outside the ws-radius',e16.8)
      db=0.0

      do 80 i=jri,np
         dl=1.0
         if (i.eq.jri.or.i.eq.np) go to 70
            dl=dl+dl
            if ((i-2*(i/2)).eq.0) dl=dl+dl
   70    dd=d(i)*dl
         db=db+dd
   80 continue

      db=dpas*db/3.0
      df0=d(jri-1)
      df1=d(jri)
      df2=d(jri+1)
      dcor=-dpas*(df1*ddjri+(df2+df0-2.0*df1)*ddjri**3/6.0+(df2-df0)
     1 *ddjri**2*.25)
      db=db+dcor
      if (iprint .ge. 5)  write(16,90) db
   90 format (1h ,' db= ',e16.8)

      call potslw (dvn,d,dp,dr,dpas,np)

      du=da*3.0/(ws**3)

      do 120 i=1,np
         if (i.gt.jr1+1) then
            srho(i)=0.0
            go to 100
         endif
            d(i)=d(i)+du*dr(i)**2
            srho(i)=d(i)
  100    continue
         dumm=-exchan(d(i),dr(i),dexv)/dr(i)
         dvf(i)=dumm
         if (i.gt.jr1) go to 110
            dvn(i)=dvn(i)-z/dr(i)+da*(1.50/ws-.50*dr(i)**2/ws**3)-db
            go to 120
  110    continue
            dvn(i)=0.0
  120 dv(i)=dvn(i)+dumm

c ad1 write the mt index and radius
      if (iprint .ge. 5)  write(16,55)jr1,dr(jr1)
  55  format(' jr1 = ',i10,10x,'wigner-seitz radius = ',e16.8)

c ad1 output 2.*dvn*r**2 for use in phase (dvn = normalised coulomb)
c     write(17,200)((2.0*dvn(i)*dr(i)*dr(i)),i=1,np)
c 200 format(1p5e16.8)
c      passvc formerly used to pass data directly to PHASE
c      do 151  i = 1, np
c         passvc (i) = 2.0 * dvn(i) * dr(i) * dr(i)
c  151 continue
c
c     passvc above is vcoul*r**2
      do 151  i = 1, np
         vcoul(i) = 2 * dvn(i)
  151 continue


c jm  output renormalized charge density for use in cphase
c                                          (d=4pi*rho*r^2)
c     write(18,200) srho

cjm write out rs as function of r
c     do 8934 i=1,jr1
c     xxrs=(3*dr(i)*dr(i)/srho(i))**.33333333
c8934 write(29,140) dr(i), xxrs
      return
      end
      subroutine rhl (rs, xk, erl, eim)
      implicit double precision (a-h, o-z)

c     input:  rs, xk
c     output: erl, eim

c     This is a new hl subroutine, using interpolation for the
c     real part while the imaginary part is calculated analytically.
c     It uses hl to calculate values at the mesh points for the inter-
c     polation of the real part. The imaginary part is calculated
c     using subroutine imhl.
c
c     written by jose mustre
c     polynomial in rs has a 3/2 power term. j.m.


c     for the right branch the interpolation has the form:
c     hl(rs,x) = e/x + f/x**2 + g/x**3
c     where e is known and
c        f = sum (i=1,3) ff(i) rs**(i+1)/2
c        g = sum (i=1,3) gg(i) rs**(i+1)/2
c
c
c     lrs=number of rs panels, in this case one has 4 panels
c     nrs=number of standard rs values, also order of rs expansion
c     if you change nrs you need to change the expansion of hl
c     in powers of rs that only has 3 terms!
c     nleft=number of coefficients for x<x0
c     nright=number of coefficients for x>x0

      parameter (lrs=4, nrs=3, nleft=4, nright=2)

      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      dimension cleft(nleft), cright(nright)

      save rcfl, rcfr
      dimension rcfl(lrs,nrs,nleft), rcfr(lrs,nrs,nright)
      data rcfr/-0.173963d+00,-0.173678d+00,-0.142040d+00,-0.101030d+00,
     1     -0.838843d-01,-0.807046d-01,-0.135577d+00,-0.177556d+00,
     2     -0.645803d-01,-0.731172d-01,-0.498823d-01,-0.393108d-01,
     3     -0.116431d+00,-0.909300d-01,-0.886979d-01,-0.702319d-01,
     4      0.791051d-01,-0.359401d-01,-0.379584d-01,-0.419807d-01,
     5     -0.628162d-01, 0.669257d-01, 0.667119d-01, 0.648175d-01/
      data rcfl/ 0.590195d+02, 0.478860d+01, 0.812813d+00, 0.191145d+00,
     1     -0.291180d+03,-0.926539d+01,-0.858348d+00,-0.246947d+00,
     2      0.363830d+03, 0.460433d+01, 0.173067d+00, 0.239738d-01,
     3     -0.181726d+03,-0.169709d+02,-0.409425d+01,-0.173077d+01,
     4      0.886023d+03, 0.301808d+02, 0.305836d+01, 0.743167d+00,
     5     -0.110486d+04,-0.149086d+02,-0.662794d+00,-0.100106d+00,
     6      0.184417d+03, 0.180204d+02, 0.450425d+01, 0.184349d+01,
     7     -0.895807d+03,-0.318696d+02,-0.345827d+01,-0.855367d+00,
     8      0.111549d+04, 0.156448d+02, 0.749582d+00, 0.117680d+00,
     9     -0.620411d+02,-0.616427d+01,-0.153874d+01,-0.609114d+00,
     1      0.300946d+03, 0.109158d+02, 0.120028d+01, 0.290985d+00,
     2      -0.374494d+03,-0.535127d+01,-0.261260d+00,-0.405337d-01/

c
c     calculate hl using interpolation coefficients
      rkf = fa/rs
      ef  = rkf**2/2
      wp  = sqrt (3/rs**3)
      call imhl (rs, xk, eim, icusp)

c     eim already has a factor of ef in it j.m.
c     eim also gives the position of the cusp

      xx = xk / rkf
c     set to fermi level if below fermi level
      if (xx .lt. 1.00001) then
          xx = 1.00001
      endif
c     calculate right hand side coefficients
      if (rs .lt. 0.2) then
         mrs=1
      elseif (rs .lt. 1.0) then
         mrs=2
      elseif (rs .lt. 5.0) then
         mrs=3
      else
         mrs=4
      endif

      do 210 j=1,nright
         cright(j) = rcfr(mrs,1,j)*rs + rcfr(mrs,2,j)*rs*sqrt(rs)
     1               + rcfr(mrs,3,j)*rs**2
  210 continue
      eee=-pi*wp/(4*rkf*ef)

      if (icusp .ne. 1) then
         do 230 j=1,nleft
            cleft(j) = rcfl(mrs,1,j)*rs + rcfl(mrs,2,j)*rs**1.5
     1                 + rcfl(mrs,3,j)*rs**2
  230    continue
         erl=cleft(1)
         do 250 j=2,nleft
            erl=erl+cleft(j)*xx**(j-1)
  250    continue
      else
c        right branch
         erl=eee/xx
         do 280 j=1,nright
            erl=erl+cright(j)/xx**(j+1)
  280    continue
      endif

      erl = erl * ef

      return
      end
      subroutine rot3i (lxp1, mxp1, ileg)
      implicit double precision (a-h,o-z)

c     input:  lxp1, mxp1, ileg (lmax+1, mmax+1)
c             also beta(ileg) used from common /pdata/
c     output: dri(...ileg) in common /rotmat/

c     subroutine rot3 calculates rotation matrices for l = 0,lxp1-1

c     subroutine rot3 calculates the beta dependence of rotation
c     matrix elements using recursion of an iterated version of
c     formula (4.4.1) in edmonds.
c
c     first written:(september 17,1986) by j. mustre
c     version 2  (17 sep 86)
c     version 3  (22 feb 87) modified by j. rehr
c     version for genfmt, modified by s. zabinsky, Sept 1991
c     Initialized dri0.  Some elements may be used before being
c        initialized elsewhere -- rot3i needs to be carefully
c        checked.  S. Zabinsky, April 1993
c
c******************** warning******************************************
c     ltot must be at least lxp1 or overwriting will occur
c     nmax must be at least nm or overwriting will occur
c----------------------------------------------------------------------
c     notation dri0(l,m,n) =  drot_i(l'm'n')
c     l = l'+1, n' = n-l, m' = m-l, primes denoting subscripts
c     thus dri0(1,1,1) corresponds to the rotation matrix with
c     l' = 0, and n' and m' = 0; dri0(3,5,5) : l' = 2,n' = 2,m' = 2.
c--------------------------------------------------------------------


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      save /rotmat/
      common /rotmat/ dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines

c     dri0 is larger than needed for genfmt, but necessary for
c     this calculation algorithm.  Copy result into smaller
c     dri arrays (in common) at end of this routine.
      dimension  dri0 (ltot+1, 2*ltot+1, 2*ltot+1)

c     initialize dri0
      do 200 il = 1, ltot+1
         do 200 im = 1, 2*ltot+1
            do 200 in = 1, 2*ltot+1
               dri0(il,im,in) = 0
  200 continue

      nm = mxp1
      ndm = lxp1+nm-1
      xc = cos(beta(ileg)/2)
      xs = sin(beta(ileg)/2)
      s = sin(beta(ileg))
      dri0(1,1,1) = 1
      dri0(2,1,1) = xc**2
      dri0(2,1,2) = s/sqrt(2.0d0)
      dri0(2,1,3) = xs**2
      dri0(2,2,1) = -dri0(2,1,2)
      dri0(2,2,2) = cos(beta(ileg))
      dri0(2,2,3) = dri0(2,1,2)
      dri0(2,3,1) = dri0(2,1,3)
      dri0(2,3,2) = -dri0(2,2,3)
      dri0(2,3,3) = dri0(2,1,1)
      do 30  l = 3, lxp1
         ln = 2*l - 1
         lm = 2*l - 3
         if (ln .gt. ndm)  ln = ndm
         if (lm .gt. ndm)  lm = ndm
         do 20  n = 1, ln
            do 10  m = 1, lm
               t1 = (2*l-1-n) * (2*l-2-n)
               t = (2*l-1-m) * (2*l-2-m)
               f1 = sqrt (t1/t)
               f2 = sqrt ((2*l-1-n) * (n-1) / t)
               t3 = (n-2) * (n-1)
               f3 = sqrt(t3/t)
               dlnm = f1 * xc**2 * dri0(l-1,n,m)
               if (n-1 .gt. 0) dlnm = dlnm - f2*s*dri0(l-1,n-1,m)
               if (n-2 .gt. 0) dlnm = dlnm + f3*xs**2*dri0(l-1,n-2,m)
               dri0(l,n,m) = dlnm
               if (n .gt. (2*l-3))
     1            dri0(l,m,n) = (-1)**(n-m) * dri0(l,n,m)
   10       continue
            if (n .gt. (2*l-3)) then
               dri0(l,2*l-2,2*l-2) = dri0(l,2,2)
               dri0(l,2*l-1,2*l-2) = -dri0(l,1,2)
               dri0(l,2*l-2,2*l-1) = -dri0(l,2,1)
               dri0(l,2*l-1,2*l-1) = dri0(l,1,1)
            endif
   20    continue
   30 continue
   40 continue

c-----test sum rule on d
c     open (29,file='rotmat.dat',status='new',carriagecontrol='list')
c     write(29,*)  ' l, m, sum'
c     write(29,*) ' (dri0(il,im,in),in = 1,ln)'
c     do 70 il = 1,lxp1
c        l = il-1
c        ln = 2*l+1
c        if(ln.gt.ndm) ln = ndm
c        do 37 im = 1,ln
c           sum = 0
c           do 50 in = 1,ln
c              m = im-il
c              term = dri0(il,im,in)
c  50       sum = sum+term**2
c           write(29,60) l,m,sum
c           write(29,62) (dri0(il,im,in),in = 1,ln)
c  60       format(2i3,e30.20)
c  62       format(5e14.6)
c  70 continue
c     close(29)
c-----end test------------------------

c     Copy result into dri(...ileg) in /rotmat/ (zero it first...)
      do 90  il = 1, ltot+1
         do 90  m1 = 1, 2*mtot+1
            do 90  m2 = 1, 2*mtot+1
               dri(il,m1,m2,ileg) = 0
   90 continue

      do 120  il = 1, lxp1
         mx = min (il-1, mxp1-1)
         do 110  m1 = -mx, mx
            do 100  m2 = -mx, mx
               dri(il,m1+mtot+1,m2+mtot+1,ileg)=dri0(il,m1+il,m2+il)
  100       continue
  110    continue
  120 continue

      return
      end
      subroutine rphbin (in)
      implicit double precision (a-h, o-z)

c     Reads input from unit in.  Returns (via /pdata/)
c       energy mesh (ne, em and eref),
c       ph (npot, lmax, lmaxp1, ph),
c       final state (l0, il0)
c
c     phmin is min value to use for |phase shift|


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines


      parameter (phmin = 1.0d-8)

c     These header lines do not include carriage control
      read(in) ntext
      do 62  i = 1, ntext
         read(in) text(i)
         read(in) ltext(i)
   62 continue
      read(in) ne, npot, ihole, rnrmav, xmu, edge, ik0
      read(in) (em(ie),ie=1,ne)
      read(in) (eref(ie),ie=1,ne)
      lmaxp1 = 0
      do 80  iph = 0, npot
         read(in) lmax0, iz(iph)
         read(in) potlbl(iph)
         do 70  ie = 1, ne
            read(in)  (ph(ie,ll,iph), ll=1,lmax0+1)
            lmax(ie,iph) = 0
c           Set lmax to include only non-zero phases
            do 60  il = 1, lmax0+1
               if (abs(ph(ie,il,iph)) .lt. phmin)  goto 61
               lmax(ie,iph) = il-1
   60       continue
   61       continue
            if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
   70    continue
   80 continue

c-----l0 is angular momentum of final state
c     Selection rule says that final state has angmom = l_init+1
c     ihole  initial state from ihole         final state
c     1      K    1s      L=0 -> linit=0   L0=1 -> lfinal=1
c     2      LI   2s      L=0 -> linit=0   L0=1 -> lfinal=1
c     3      LII  2p 1/2  L=1 -> linit=1   L0=2 -> lfinal=2
c     4      LIII 2p 3/2  L=1 -> linit=1   L0=2 -> lfinal=2
c     5+     M -- think about this later...
      if (ihole .le. 2)  then
c        hole in s state (1s or 2s)
         linit = 0
         lfinal = 1
      elseif (ihole .le. 4)  then
c        hole in p state (2p 1/2  or  2p 3/2)
         linit = 1
         lfinal = 2
      else
c        some m hole, n=3, could go to d state
         stop 'Can not handle M shell.'
      endif
      l0 = lfinal
      il0 = l0 + 1

      return
      end
      subroutine rpotph (io, nhead0, head0, lhead0,
     1             nat, nph, nfr, ihole, gamach, iafolp, intclc,
     1             ixc, vr0, vi0, rs0, iphat, rat, iatph, ifrph, 
     1             xnatph, novr,
     2             iphovr, nnovr, rovr, folp, ion, iz, iprint, 
     2             ixanes, nemax, xkmin, xkmax, potlbl)
      implicit double precision (a-h, o-z)

c     Notes:
c        nat   number of atoms in problem
c        nph   number of unique potentials
c        nfr   number of unique free atoms
c        ihole hole code of absorbing atom
c        iph=0 for central atom
c        ifr=0 for central atom
c        xkmin, xkmax  min and max energy mesh points to consider


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      character*(*) head0(nhead0)
      dimension lhead0(nhead0)

c     End of line comments removed -- see include file arrays.h for
c     comments.
c     Specific atom input data
      dimension iphat(natx)
      dimension rat(3,natx)

c     Unique potential input data
      dimension iatph(0:nphx)
      dimension ifrph(0:nphx)
      dimension xnatph(0:nphx)
      character*6  potlbl(0:nphx)

      dimension folp(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)

c     Free atom data
      dimension ion(0:nfrx)
      dimension iz(0:nfrx)

c     read and save header from old file, has carriage control char
      head0(1) = ' '
      call rdhead (io, nhead0, head0, lhead0)
      read(io,*) ihole, gamach, iprint, iafolp, intclc
      read(io,*) ixc, vr0, vi0, rs0
      read(io,*) ixanes, nemax, xkmin, xkmax
      read(io,*) nfr
      do 710  ifr = 0, nfr
         read(io,*)  index, iz(ifr), ion(ifr)
  710 continue
      read(io,*) nat
      do 720  iat = 1, nat
         read(io,*) index, iphat(iat), (rat(j,iat),j=1,3)
  720 continue
      read(io,*) nph
      do 740  iph = 0, nph
         read(io,*) index, iatph(iph), ifrph(iph), xnatph(iph),
     1                folp(iph), novr(iph)
         read(io,*) potlbl(iph)
         do 730  iovr = 1, novr(iph)
            read(io,*) iphovr(iovr,iph), nnovr(iovr,iph),
     1                   rovr(iovr,iph)
  730    continue
  740 continue

      return
      end
      subroutine sclmz (rho, lmaxp1, mmaxp1, ileg)
      implicit double precision (a-h, o-z)

c     Set CLM(Z) for current leg.
c     Makes clm(z) (eq B11).  Fills array clmi in /clmz/ for ileg,
c     elements clm(0,0) -> clm(lmax+1,mmax+1).
c     If mmaxp1 > lmaxp1, fills m only to lmaxp1.

c     calculates energy dependent factors
c     c(il,im) = c_l^(m)z**m/m! = c_lm    by recursion
c     c_l+1,m = c_l-1,m-(2l+1)z(c_l,m-c_l,m-1, l ne m
c     c_m,m = (-z)**m (2m)!/(2**m m!) with z = 1/i rho
c
c     To test pw approx, set z = 0


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      save /clmz/
      complex*16 clmi
      common /clmz/ clmi(ltot+1,mtot+ntot+1,legtot)


      complex*16 rho(legtot)
      complex*16 z, cmm

      cmm = 1
      z = -coni / rho(ileg)

      clmi(1,1,ileg) = (1,0)
      clmi(2,1,ileg) = clmi(1,1,ileg) - z

      lmax = lmaxp1-1

      do 10  il = 2, lmax
         clmi(il+1,1,ileg) =
     1           clmi(il-1,1,ileg) - z*(2*il-1)*clmi(il,1,ileg)
   10 continue
      mmxp1 = min (mmaxp1, lmaxp1)
      do 20  im = 2, mmxp1
         m = im-1
         imp1 = im+1
         cmm = -cmm * (2*m-1) * z
         clmi(im,im,ileg) = cmm
         clmi(imp1,im,ileg) = cmm * (2*m+1) * (1-im*z)
         do 20  il = imp1, lmax
            l = il-1
            clmi(il+1,im,ileg) = clmi(l,im,ileg) -
     1          (2*l+1) * z * (clmi(il,im,ileg) + clmi(il,m,ileg))
   20 continue

      return
      end
      double precision function sdist(r0, r1)
      implicit double precision (a-h, o-z)
c     find distance squared between cartesian points r0 and r1
c     single precision
      dimension r0(3), r1(3)
      sdist = 0
      do 10  i = 1, 3
         sdist = sdist + (r0(i) - r1(i))**2
   10 continue
      sdist = sqrt(sdist)
      return
      end
      subroutine setgam (iz, ihole, gamach)

c     Sets gamach, core hole lifetime.  Data comes from graphs in
c     K. Rahkonen and K. Krause,
c     Atomic Data and Nuclear Data Tables, Vol 14, Number 2, 1974.

      implicit double precision (a-h, o-z)

      dimension gamk(6), zk(6),  famk(6)
      dimension gaml1(6), zl1(6),faml1(6)
      dimension gaml2(6), zl2(6),faml2(6)
      parameter (ryd  = 13.6058)

      save ienter

c     Note that 0.99 replaces 1.0, 95.1 replaces 95.0 to avoid roundoff
c     trouble.
c     Gam arrays contain the gamma values.
c     We will take log10 of the gamma values so we can do linear
c     interpolation from a log plot.

      data  zk   / 0.99d0,  10.0d0, 20.0d0,  40.0d0,  60.0d0,   95.1d0/
c      data  gamk / 0.07,   0.3,  0.75,  5.0,  20.0,  100.0/
      data  famk / 0.07d0,   0.3d0,  0.75d0,  5.0d0,  20.0d0,  100.0d0/

      data  zl1   / 0.99d0,  20.0d0, 35.0d0, 50.0d0,  75.0d0,  95.1d0/
c      data  gaml1 / 0.07,   4.0,  7.0,  4.0,   8.0,  19.0/
      data  faml1 / 0.07d0,   4.0d0,  7.0d0,  4.0d0,   8.0d0,  19.0d0/

      data  zl2   / 0.99d0,  26.0d0, 31.0d0, 60.0d0,  80.0d0,  95.1d0/
c      data  gaml2 / 0.001,  1.7,  0.8,  3.5,   5.0,  10.0/
      data  faml2 / 0.001d0,  1.7d0,  0.8d0,  3.5d0,   5.0d0,  10.0d0/

      data ienter /0/

c     Call this only once, if it gets called a second time the gamma
c     values will be messed up by repeated taking of log10

c      if (ienter .gt. 0)  then
c         write(77,*) ' Re-entered SETGAM'
c         stop 'SETGAM-1'
c      endif
c      ienter = 1

      if (ihole .le. 0)  then
         gamach = 0
         write(77,*) 'No hole in SETGAM, gamach = ', gamach
         return
      endif
      if (ihole .gt. 4)  then
         write(77,*) ' This version of FEFF only handles through L III',
     1              ' shell absorption.'
         stop 'SETGAM-2'
      endif

      zz = iz
      if (ihole .le. 1)  then
         do 10  i = 1, 6
c            gamk(i) = log10 (gamk(i))
            gamk(i) = log10 (famk(i))
   10    continue
         call terp (zk, gamk, 6, zz, gamach)
      else if (ihole .le. 2)  then
         do 20  i = 1, 6
c            gaml1(i) = log10 (gaml1(i))
            gaml1(i) = log10 (faml1(i))
   20    continue
         call terp (zl1, gaml1, 6, zz, gamach)
      else if (ihole .le. 4)  then
c        note that LII and LIII have almost exactly the same
c        core hole lifetimes
         do 30  i = 1, 6
c            gaml2(i) = log10 (gaml2(i))
            gaml2(i) = log10 (faml2(i))
   30    continue
         call terp (zl2, gaml2, 6, zz, gamach)
      endif

c     Change from log10 (gamma) to gamma
      gamach = 10.0 ** gamach

c     Table values are in eV, code requires atomic units
      gamach = gamach / ryd

      return
      end
      subroutine setlam (icalc, ie)
      implicit double precision (a-h, o-z)

c     Set lambda array based on icalc and ie
c     icalc  what to do
c      0     i0, ss exact
c      1     i1, ss exact
c      2     i2, ss exact
c     10     cute algorithm
c     <0     do exactly as told, decode as:
c               icalc = -(nmax + 100*mmax + 10 000*(iord+1))
c               Note that iord=0 <=> nmax=mmax=0, so use
c                  icalc = -10 000 for this case.
c               iord = 2*nmax + mmax, so if you want iord to control,
c               set nmax and mmax large enough-- if you want nmax and
c               mmax to control, set iord = 2*nmax + mmax...

c     inputs: ie used for cute algorithm
c             nsc used from /pdata/ to recognize ss paths
c     output: variables in /lambda/ set


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      common /lambda/  
     4   mlam(lamtot), 	!mu for each lambda
     5   nlam(lamtot),	!nu for each lambda
     1   lamx, 		!max lambda in problem
     2   laml0x, 	!max lambda for vectors involving absorbing atom
     3   mmaxp1, nmax 	!max mu in problem + 1, max nu in problem


c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines

      dimension mlam0(lamtot), nlam0(lamtot)

c     one degree in radians
      parameter (onedeg = .01745329252d0)

c     Set iord, nmax and mmax based on icalc
      if (icalc .lt. 0)  then
c        decode it and do what user wants
         icode = -icalc
         nmax = mod(icode,100)
         mmax = mod(icode,10000)/100
         iord = icode/10000 -1
      elseif (nsc .eq. 1)  then
         mmax = il0-1
         nmax = il0-1
         iord = 2*nmax + mmax
      elseif (icalc .lt. 10)  then
         iord = icalc
         mmax = iord
         nmax = iord/2
      elseif (icalc .eq. 10)  then
c        do cute algorithm
c        set mmax = L0 if straight line path, otherwise set mmax = 3
         mmax = il0-1
         do 10  ileg = 1, nleg
            mag1 = abs(beta(ileg))
            mag2 = abs(mag1 - pi)
c           if beta is not 0 or pi, path is non-linear
            if (mag1.gt.onedeg .and. mag2.gt.onedeg) mmax = 3
   10    continue
c        Set nmax based on ie and l0.
c        k <= 12 invA (ie=41)  nmax = L0
c        k >= 13 invA (ie=42)  nmax =  9
         nmax = il0-1
         if (ie .ge. 42)  nmax = 9
         iord = 2*nmax + mmax
      else
         write(77,*) 'undefined icalc ', icalc
         stop 'setlam'
      endif

c-----construct index lambda (lam), (mu, nu) = mlam(lam), nlam(lam)
c     lamtot, ntot, mtot are maximum lambda, mu and nu to consider
c     Use ...0 for making indices, then sort into arrays with no
c     trailing 0 so laml0x is minimimized. (note: this is a crude
c     n**2 sort -- can 'improve' to nlog_2(n) if necessary)
      lam = 0
      do 20 in = 1, nmax+1
         n = in - 1
         do 20  im = 1, mmax+1
            m = im-1
            jord = 2*n+m
            if (jord .gt. iord)  goto 20
            if (lam .ge. lamtot)  then
               write(77,*) 'Lambda array filled, some order lost'
               goto 21
            endif
            lam = lam+1
            mlam0(lam) = -m
            nlam0(lam) = n
            if (m .eq. 0)  goto 20
            if (lam .ge. lamtot)  then
               write(77,*) 'Lambda array filled, some order lost'
               goto 21
            endif
            lam = lam+1
            mlam0(lam) = m
            nlam0(lam) = n
   20 continue
   21 continue
      lamx=lam
c     lamx must be less than lamtot
      if (lamx .gt. lamtot) stop 'SETLAM lamx > lamtot'

c     laml0x is biggest lam for non-zero fmatrix, also set mmax and nmax
c     Sort mlam0 and nlam0 to use min possible laml0x
      lam = 0
      do 30  lam0 = 1, lamx
         if ((nlam0(lam0).le.l0) .and. (iabs(mlam0(lam0)).le.l0)) then
            lam = lam+1
            nlam(lam) = nlam0(lam0)
            mlam(lam) = mlam0(lam0)
            nlam0(lam0) = -1
         endif
   30 continue
      laml0x = lam
      do 40  lam0 = 1, lamx
         if (nlam0(lam0) .ge. 0)  then
            lam = lam+1
            nlam(lam) = nlam0(lam0)
            mlam(lam) = mlam0(lam0)
         endif
   40 continue

      mmaxp1 = 0
      nmax = 0
      do 50  lam = 1, lamx
         if (mlam(lam)+1 .gt. mmaxp1)  mmaxp1 = mlam(lam)+1
         if (nlam(lam) .gt. nmax)  nmax = nlam(lam)
   50 continue

      if (nmax.gt.ntot .or. mmaxp1.gt.mtot+1)  then
         write(77,*) 'mmaxp1, nmax, mtot, ntot ',
     1            mmaxp1, nmax, mtot, ntot
         write(77,*) 'icalc ', icalc
         stop 'setlam'
      endif

      return
      end
      subroutine sidx (rholap, npts, rmt, rnrm, imax, imt, inrm)
      implicit double precision (a-h, o-z)
      dimension rholap (npts)

      imt = ii (rmt)
      inrm = ii (rnrm)

c     Set imax (last non-zero rholap data)
      do 220  i = 1, npts
         if (rholap(i) .le. 1.0d-5)  goto 230
         imax = i
  220 continue
  230 continue

c     We need data up to the norman radius, so move norman
c     radius if density is zero inside rnrm.
      if (inrm .gt. imax)  then
         inrm = imax
         rnrm = rr (inrm)
         write(77,*) ' Moved rnrm.  New rnrm (au) ', rnrm
      endif
      if (imt .gt. imax)  then
         imt = imax
         rmt = rr (imt)
         write(77,*) ' Moved rmt.  New rmt (au) ', rmt
      endif
      return
      end
c---------------------------------------------------------------------
c     program sigms.f
c
c     calculates debye-waller factors for each multiple
c     scattering path using Debye-Model correlations
c
c     files:  input  pathd_all.dat  multiple scattering path data
c             output fort.3  sig**2 vs path
c                    fort.2  long output
c
c     version 1  (29 july 91)
c
c     coded by j. rehr
c     path data from s. zabinsky
c
c     modified to use pdata.inp, Dec 1991, siz
c     Subroutine version, Dec 1991, siz
c
c---------------------------------------------------------------------

      subroutine sigms (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
c               tk temperature in degrees K
c               thetad debye temp in degrees K
c               rs=wigner seitz or norman radius in bohr, averaged
c                  over entire problem
c                  (4pi/3)*rs**3 = sum( (4pi/3)rnrm**3 ) / N
c                  (sum is over all atoms in the problem)
c               nlegx used in dimensions of rat and iz
c               nleg nlegs in path
c               rat positions of each atom in path (in bohr)
c               iz atomic number of each atom in path
c               NB Units of distance in this routine
c                  are angstroms, including sig**2
c               sig2 is output, debye waller factor in bohr**-2

      implicit double precision (a-h,o-z)


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


c     nlegx is max number of atoms in any one path
      dimension rat(3,0:nlegx)
      dimension iz(0:nlegx)

c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)

c     Note:  There are nleg atoms including the central atom
c            index 0 and index nleg both refer to central atom,
c            which makes special code unnecessary later.
      sum = 0.0d0
      ntot = 0

      sigtot=0
      do 800 il=1,nleg
      do 800 jl=il,nleg

c        calculate r_i-r_i-1 and r_j-r_j-1

         rij = feff_dist (rat(1,il), rat(1,jl))
         call corrfn (rij, cij, thetad, tk, iz(il), iz(jl), rs)
         sig2ij=cij

         rimjm = feff_dist (rat(1,il-1), rat(1,jl-1))
         call corrfn (rimjm, cimjm, thetad, tk, iz(il-1), iz(jl-1), rs)
         sig2ij=sig2ij+cimjm

         rijm = feff_dist (rat(1,il), rat(1,jl-1))
         call corrfn (rijm, cijm, thetad, tk, iz(il), iz(jl-1), rs)
         sig2ij=sig2ij-cijm

         rimj = feff_dist (rat(1,il-1), rat(1,jl))
         call corrfn (rimj, cimj, thetad, tk, iz(il-1), iz(jl), rs)
         sig2ij=sig2ij-cimj

         riim = feff_dist (rat(1,il), rat(1,il-1))
         rjjm = feff_dist (rat(1,jl), rat(1,jl-1))

         ridotj=(rat(1,il)-rat(1,il-1))*(rat(1,jl)-rat(1,jl-1))+
     1          (rat(2,il)-rat(2,il-1))*(rat(2,jl)-rat(2,jl-1))+
     2          (rat(3,il)-rat(3,il-1))*(rat(3,jl)-rat(3,jl-1))
         ridotj=ridotj/(riim*rjjm)

c        double count i .ne. j  terms
         if(jl.ne.il) sig2ij=2*sig2ij
         sig2ij=sig2ij*ridotj
         sigtot=sigtot+sig2ij

  800 continue
      sig2=sigtot/4.0d0

c     sig2 is in bohr**2, just as we wanted for ff2chi
      return
      end



      subroutine corrfn(rij,cij,thetad,tk,iz1,iz2,rsavg)
c     subroutine calculates correlation function
c     c(ri,rj)=<xi xj> in the Debye approximation
c
c             =(1/N)sum_k exp(ik.(Ri-Rj))(1/sqrt(mi*mj))*
c              (hbar/2w_k)*coth(beta hbar w_k/2)
c             = (3kT/mu w_d**2)*sqrt(mu**2/mi*mj)*I
c
c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)
c
c      solution by numerical integration
c
      implicit double precision (a-h, o-z)
      common /xy/ x, yinv


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508459393094d0)

c     external fn
c     rij=2.55
c     tk=295
c     thetad=315
c     ami=amj=63.55 at wt for Cu
c     rs=2.7

      ami=atwtd(iz1)
      amj=atwtd(iz2)
      rs=rsavg
c     thetad in degrees K, t temperature in degrees K
c     y=thetad/tk
      yinv=tk/thetad
      xkd=(9.0d0*pi/2.0d0)**(third)/(rs*bohr)
      fac=(3.0d0/2.0d0)*con/(thetad*sqrt(ami*amj))
      rj=rij
      x=xkd*rj
c     call numerical integration
      call bingrt (grater, eps, nx)
      cij=fac*grater
      return
      end
      double precision function fn(w)
      implicit double precision (a-h,o-z)
      common/xy/x,yinv
c     fn=(sin(wx)/x)*coth(wy/2)
c     change code to allow t=0 without bombing
c     fn=2/y
      fn=2.0d0*yinv
      if(w.lt.1.d-20) return
      fac=w
      if(x.gt.0.0d0) fac=sin(w*x)/x
      emwy=0.0d0
      if(yinv.gt.0.0125d0) emwy=exp(-w/yinv)
      emwy=exp(-w/yinv)
      fn=fac*(1.0d0+emwy)/(1.0d0-emwy)
      return
      end
c-----------------------------------------------
      subroutine bingrt (b, eps, n)
c     subroutine calculates integrals between [0,1]
c      b = int_0^1 f(z) dz
c     by trapezoidal rule and binary refinement
c     (romberg integration)
c     coded by j rehr (10 Feb 92)
c     see, e.g., numerical recipes for discussion
c     and a much fancier version
c-----------------------------------------------
c     del=dz  itn=2**n tol=1.e-5
c     starting values
      implicit double precision (a-h,o-z)
      common /xy/x,yinv
c     external fn
c     error is approximately 2**(-2n) ~ 10**(-.6n)
c     so nmax=10 implies an error of 1.e-6
      parameter(nmax = 10, tol = 1.d-5)
      parameter(zero=0, one=1)
      n=0
      itn=1
      del=1.0d0
      bn=(fn(zero)+fn(one))/2.0d0
      bo=bn
 10   continue
c     nth iteration
c     b_n+1=(b_n)/2+deln*sum_0^2**n f([2n-1]deln)
      n=n+1
      if(n.gt.nmax) go to 40
      del=del/2.0d0
      sum=0.0d0
      do 20 i=1, itn
      zi=(2*i-1)*del
 20   sum=sum+fn(zi)
c     bnp1=b_n+1 is current value of integral
      bnp1=bn/2.0d0+del*sum
c     cancel leading error terms b=[4b-bn]/3
c     note: this is the first term in the
c     neville table - remaining errors were
c     found too small to justify the added code
      b=(4*bnp1-bn)/3.0d0
      eps=abs((b-bo)/b)
      if(eps.lt.tol) goto 60
      bn=bnp1
      bo=b
      itn=itn*2
      goto 10
 40   write(77,50) n,itn, b,eps
 50   format(' not converged, n,itn,b,eps=',
     1  2i4,2e14.6)
      return
 60   continue
c     print70, n, itn, b, eps
c70   format(' n,itn,b,eps=' 2i4,2e16.8)
      return
      end
      subroutine snlm (lmaxp1, mmaxp1)
      implicit double precision(a-h,o-z)

c     Set nlm, legendre normalization factors, xnlm in common /nlm/
c     Calculates legendre norm factors
c     xnlm= sqrt ((2l+1)(l-m)!/(l+m)!)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      save /nlm/
      common /nlm/ xnlm(ltot+1,mtot+1)


c     flg(i) = i! * afac**i, set in factst
      dimension flg(0:210)

      call factst (afac, flg)

c     initialize xnlm explicitly
      do 5  il = 1, ltot+1
      do 5  im = 1, mtot+1
         xnlm(il,im) = 0
    5 continue

      do 10  il = 1, lmaxp1
         mmxp1 = min (mmaxp1, il)
         do 10  im = 1, mmxp1
            l = il-1
            m = im-1
            cnlm = (2*l+1) * flg(l-m) / flg(l+m)
            cnlm = sqrt(cnlm) * afac**m
            xnlm(il,im) = cnlm
   10 continue

      return
      end
      subroutine factst (afac, flg)
      implicit double precision (a-h,o-z)

c     FACTorial SeT, flg(i) = i! * afac**i
      dimension flg(0:210)

c     afac = 1/64 works with double precision on a VAX
      afac = 1.0d0/64.0d0

      flzero = 1
      flg(0) = 1
      flg(1) = afac

      do 10  i = 2, 210
   10 flg(i) = flg(i-1) * i * afac

      return
      end
      subroutine somm (dr,dp,dq,dpas,da,m,np)
c
c integration by the method of simpson of (dp+dq)*dr**m from
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      dimension dr(251), dp(251), dq(251)
      mm=m+1
      d1=da+mm
      da=0.0d0
      db=0.0d0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      if (dc) 20,40,30
   20 db=db+dc
      go to 40
   30 da=da+dc
   40 dc=dq(i)*dl
      if (dc) 50,70,60
   50 db=db+dc
      go to 70
   60 da=da+dc
   70 continue
      da=dpas*(da+db)/3.0d0
      dc=exp(dpas)-1.0d0
      db=d1*(d1+1.0d0)*dc*exp((d1-1.0d0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dc=(dr(1)**mm)*(1.0d0+1.0d0/(dc*(d1+1.0d0)))/d1
      da=da+dc*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine sortir (n, index, r)
      implicit double precision (a-h, o-z)

c     SORT by rearranges Indices, keys are Real numbers
c     Heap sort, following algorithm in Knuth using r as key
c     Knuth, The Art of Computer Programming,
c     Vol 3 / Sorting and Searching, pp 146-7
c     Array r is not modified, instead array index is returned
c     ordered so that r(index(1)) is smallest, etc.
c     rr is temporary r storage (Knuth's R), irr is index of stored r

      dimension r(n), index(n)

c     Initialize index array
      do 10  i = 1, n
         index(i) = i
   10 continue
c     only 1 element is already sorted
      if (n .eq. 1)  return

c     H1: initialize
      l = n/2 + 1
      ir = n

c     H2: Decrease l or ir
   20 continue
      if (l .gt. 1)  then
         l = l-1
         irr = index(l)
         rr = r(irr)
      else
         irr = index(ir)
         rr = r(irr)
         index(ir) = index(1)
         ir = ir-1
         if (ir .eq. 1) then
            index(1) = irr
            return
         endif
      endif

c     H3: Prepare for sift-up
      j = l

c     H4: Advance downward
   40 continue
      i = j
      j = 2 * j
      if (j .eq. ir)  goto 60
      if (j .gt. ir)  goto 80

c     H5: Find larger son of i
      if (r(index(j)) .lt. r(index(j+1)))  j = j+1

c     H6: Son larger than rr?
   60 continue
      if (rr .ge. r(index(j)))  goto 80

c     H7: Move son up
      index(i) = index(j)
      goto 40

c     H8: Store rr in it's proper place
   80 continue
      index(i) = irr
      goto 20

      end
      subroutine sortii (n, index, k)
      implicit double precision (a-h, o-z)

c     SORT by rearranges Indices, keys are Integers
c     Heap sort, following algorithm in Knuth using r as key
c     Knuth, The Art of Computer Programming,
c     Vol 3 / Sorting and Searching, pp 146-7
c     Array r is not modified, instead array index is returned
c     ordered so that r(index(1)) is smallest, etc.
c     rr is temporary r storage (Knuth's R), irr is index of stored r

      dimension k(n)
      dimension index(n)

c     Initialize index array
      do 10  i = 1, n
         index(i) = i
   10 continue
c     only 1 element is already sorted
      if (n .eq. 1)  return

c     H1: initialize
      l = n/2 + 1
      ir = n

c     H2: Decrease l or ir
   20 continue
      if (l .gt. 1)  then
         l = l-1
         irr = index(l)
         kk = k(irr)
      else
         irr = index(ir)
         kk = k(irr)
         index(ir) = index(1)
         ir = ir-1
         if (ir .eq. 1) then
            index(1) = irr
            return
         endif
      endif

c     H3: Prepare for sift-up
      j = l

c     H4: Advance downward
   40 continue
      i = j
      j = 2 * j
      if (j .eq. ir)  goto 60
      if (j .gt. ir)  goto 80

c     H5: Find larger son of i
      if (k(index(j)) .lt. k(index(j+1)))  j = j+1

c     H6: Son larger than kk?
   60 continue
      if (kk .ge. k(index(j)))  goto 80

c     H7: Move son up
      index(i) = index(j)
      goto 40

c     H8: Store kk in it's proper place
   80 continue
      index(i) = irr
      goto 20

      end
      subroutine sortid (n, index, r)

c     SORT by rearranges Indices, keys are Double precision numbers
c     Heap sort, following algorithm in Knuth using r as key
c     Knuth, The Art of Computer Programming,
c     Vol 3 / Sorting and Searching, pp 146-7
c     Array r is not modified, instead array index is returned
c     ordered so that r(index(1)) is smallest, etc.
c     rr is temporary r storage (Knuth's R), irr is index of stored r

      implicit double precision (a-h, o-z)
      dimension r(n), index(n)

c     Initialize index array
      do 10  i = 1, n
         index(i) = i
   10 continue
c     only 1 element is already sorted
      if (n .eq. 1)  return

c     H1: initialize
      l = n/2 + 1
      ir = n

c     H2: Decrease l or ir
   20 continue
      if (l .gt. 1)  then
         l = l-1
         irr = index(l)
         rr = r(irr)
      else
         irr = index(ir)
         rr = r(irr)
         index(ir) = index(1)
         ir = ir-1
         if (ir .eq. 1) then
            index(1) = irr
            return
         endif
      endif

c     H3: Prepare for sift-up
      j = l

c     H4: Advance downward
   40 continue
      i = j
      j = 2 * j
      if (j .eq. ir)  goto 60
      if (j .gt. ir)  goto 80

c     H5: Find larger son of i
      if (r(index(j)) .lt. r(index(j+1)))  j = j+1

c     H6: Son larger than rr?
   60 continue
      if (rr .ge. r(index(j)))  goto 80

c     H7: Move son up
      index(i) = index(j)
      goto 40

c     H8: Store rr in it's proper place
   80 continue
      index(i) = irr
      goto 20

      end
C FUNCTION ISTRLN (STRING)  Returns index of last non-blank
C                           character.  Returns zero if string is
C                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING

C  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

C  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING (I:I) .NE. ' ')  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
C SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*200  TMP

      JLEN = ISTRLN (STRING)

C  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

C  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING (I:I) .NE. ' ')  GOTO 20
   10 CONTINUE
   20 CONTINUE

C  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

C  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE BWORDS (S, NWORDS, WORDS)
C
C     Breaks string into words.  Words are seperated by one or more
C     blanks, or a comma and zero or more blanks.
C
C     ARGS        I/O      DESCRIPTION
C     ----        ---      -----------
C     S            I       CHAR*(*)  String to be broken up
C     NWORDS      I/O      Input:  Maximum number of words to get
C                          Output: Number of words found
C     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
C                          Contains words found.  WORDS(J), where J is
C                          greater then NWORDS found, are undefined on
C                          output.
C
C      Written by:  Steven Zabinsky, September 1984
C
C**************************  Deo Soli Gloria  **************************

C  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA
      PARAMETER (BLANK = ' ', COMMA = ',')

C  -- BETW    .TRUE. if between words
C     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

C  -- Maximum number of words allowed
      WORDSX = NWORDS

C  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

C  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

C  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSEIF (S(I:I) .EQ. COMMA)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S(BEGC : I-1)
               BETW = .TRUE.
            ELSEIF (COMFND)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = BLANK
            ENDIF
            COMFND = .TRUE.
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END
      subroutine strap (x, y, n, sum)
      implicit double precision (a-h, o-z)

c     Trapeziodal integration of y(x), result in sum
c     SINGLE PRECISION

      dimension x(n), y(n)

      sum = y(1) * (x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * (x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * (x(n) - x(n-1))
      sum = sum/2.0d0

      return
      end
c SUBROUTINE SUMAX (NPTS, RN, ANN, AA2, AASUM)
c This is a version of the subroutine sumax found on page 110 of
c Louck's book.  It performs eq 3.22, using simpson's rule and
c taking advantage of the logarithmic grid so that sum f(r)*dr becomes
c sum over f(r)*r*(0.05).  Linear interpolation is used at the end
c caps.  This version does not sum over 14 shells of identical
c atoms, instead it averages the contribution of one or more atoms
c of type 2 at the location of atom 1.  Louck's description (except
c for his integration algorithm) is very clear.
c
c input:  npts      number of points to consider
c         rn        distance from atom 1 to atom 2 in au
c         ann       number of type 2 atoms to add to atom 1, can
c                   be fractional
c         aa2(i)    potential or density at atom 2
c output: aasum(i)  spherically summed contribution added into this
c                   array so that sumax can be called repeatedly
c                   and the overlapped values summed into aasum
c
c Note that this routine requires that all position data be on a
c grid  rr(j) = exp (-8.8d0 + (j-1)*0.05d0), which is the grid
c used by Louck, and also used by ATOM if nuclear options not used.
c
c Coded by Steven Zabinsky, December 1989
c Modified for FEFF cluster code, August 1990, siz
c Bug fixed, May 1991, SIZ
c Another bug fixed, Mar 1992, SIZ
c
c T.L.Louck, Augmented Plane Wave Method, W.A.Benjamin, Inc., 1967

      subroutine sumax (npts, rn, ann, aa2, aasum)
      implicit double precision (a-h, o-z)
      parameter (nptx=250)
      dimension aa2(nptx), aasum(nptx)
      dimension stor(nptx)

c     jjchi     index beyond which aa2 is zero
c     jtop      index just below distance to neighbor
c               aasum is calculated only up to index jtop

c     Wigner-Seitz radius is set to 15 in ATOM.
      rws = 15.0d0
      jjchi = ii(rws)
      jtop  = ii(rn)

      topx = xx(jjchi)

      do 120  i = 1, jtop
         x = xx(i)
         xint = 0.0d0
         et = exp(x)
         blx = log(rn-et)
         if (blx .ge. topx)  goto 119
         jbl = 2.0d0+20.0d0*(blx+8.8d0)
         if (jbl .lt. 1)  jbl=1
         if (jbl .ge. 2)  then
c           use linear interp to make end cap near center of neighbor
            xjbl = jbl
            xbl = 0.05d0 * (xjbl-1.0d0) - 8.8d0
            g = xbl-blx
            xint =xint+0.5d0*g*(aa2(jbl)*(2.0d0-20.0d0*g)*exp(2.0d0*xbl)
     1             +20.0d0*g*aa2(jbl-1)*exp(2.0d0*(xbl-0.05d0)))
         endif
         tlx = log(rn+et)
         if (tlx .ge. topx)  then
            jtl = jjchi
            go to 90
         endif
         jtl = 1.0d0 + 20.0d0*(tlx+8.8d0)
         if (jtl .lt. jbl)  then
c           handle peculiar special case at center of atom 1
            fzn = aa2(jtl)*exp(2.0d0*(xbl-0.05d0))
            fz3 = aa2(jbl)*exp(2.0d0*xbl)
            fz2 = fzn+20.0d0*(fz3-fzn)*(tlx-xbl+0.05d0)
            fz1 = fzn+20.0d0*(fz3-fzn)*(blx-xbl+0.05d0)
            xint = 0.5d0*(fz1+fz2)*(tlx-blx)
            go to 119
         endif
         xjtl = jtl
         xtl = 0.05d0*(xjtl-1.0d0)-8.8d0
         c = tlx-xtl
         xint = xint+0.5d0*c*(aa2(jtl)*(2.0d0-20.0d0*c)
     1         *exp(2.0d0*xtl)+aa2(jtl+1)*20.0d0*c
     2         *exp(2.0d0*(xtl+0.05d0)))

   90    if (jtl .gt. jbl)  then
  100       xint = xint+0.5d0*(aa2(jbl)*exp(2.0d0*xbl)+aa2(jbl+1)
     1             *exp(2.0d0*(xbl+0.05d0)))*0.05d0
            jbl = jbl+1
            if (jbl .lt. jtl) then
               xbl = xbl+0.05d0
               go to 100
            endif
         endif
  119    stor(i) = 0.5d0*xint*ann/(rn*et)
  120 continue

      do 190  i = 1, jtop
         aasum(i) = aasum(i) + stor(i)
  190 continue

      return
      end
c     Linear interpolation and extrapolation.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terp (x, y, n, x0, y0)
      implicit double precision (a-h, o-z)

      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
c     if i < 1, set i=1, if i > n-1, set i=n-1
      i = max (i, 1)
      i = min (i, n-1)

      if (x(i+1) - x(i) .eq. 0)  stop 'TERP-1'

      y0 = y(i) +  (x0 - x(i)) * (y(i+1) - y(i)) / (x(i+1) - x(i))

      return
      end

      integer function locat(x, n, xx)
      implicit double precision (a-h, o-z)
      double precision x, xx(n)
      integer  u, m, n

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat = 0
      u = n+1

   10 if (u-locat .gt. 1)  then
         m = (u + locat) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat = m
         endif
         goto 10
      endif

      return
      end
      subroutine timrep (npat, ipat, rx, ry, rz, dhash)
      implicit double precision (a-h, o-z)

c     subroutine timrev(...) is modified for polarization case 
c     Time-orders path and returns path in standard order,
c     standard order defined below.
c     Input:  npat, ipat
c     Output: ipat in standard order (time reversed if necessary)
c             rx, ry, rz   contain x,y,z coordinates of the path atoms,
c             where z-axis is along polarization vector or first leg, if
c               running usual feff,
c             x-axis is chosen so that first atom, which does not lie on
c               z-axis, lies in xz-plane,
c               for elliptically polarized light, x-axis is along the
c               incidence direction
c             y-axis is cross product of two previos unit vectors
c             Standarrd order is defined so that first nonzero x,y and z
c             coords are positive.(Otherwise we use the inversion of
c             the corresponding unit vector)
c             dhash double precision hash key for path in standard
c                order


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt

      common /atoms/ rat(3,0:natx), ipot(0:natx), ilb(0:natx)
      dimension ipat(npatx+1), rx(npatx), ry(npatx), rz(npatx)
      dimension ipat0(npatx+1), rx0(npatx), ry0(npatx), rz0(npatx)

      double precision dhash, dhash0

c     Time reverses path if time reversing it will put it
c     in standard order.  Standard order is defined by min hash
c     number, using path hash algorithm developed for the path
c     degeneracy checker.  See subroutine phash for details.
c     Symmetrical paths are, of course, always standard ordered.
c     Also returns hash number for standard ordered path.

c     Use suffix 0 for (') in variable names

c     If no time-reversal standard ordering needed, make hash number
c     and return.  No timrev needed if 2 leg path (symmetrical).
      nleg = npat + 1
      ipat(nleg) = 0
      do 10 i = 1, npatx
         rx(i)   = 0.0d0
         ry(i)   = 0.0d0
         rz(i)   = 0.0d0
         rx0(i)   = 0.0d0
         ry0(i)   = 0.0d0
         rz0(i)   = 0.0d0
   10 continue
      call mpprmp(npat, ipat, rx, ry, rz)
      call phash (npat, ipat, rx, ry, rz, dhash)

      if (npat .le. 1)  then
         return
      endif

c     Make time reversed path

      ipat0(nleg) = ipat(nleg)
      do 210  i = 1, npat
         ipat0(i) = ipat(nleg-i)
  210 continue
      call mpprmp(npat, ipat0, rx0, ry0, rz0)
      call phash (npat, ipat0, rx0, ry0, rz0, dhash0)

c     Do the comparison using hash numbers
c     Want representation with smallest hash number
      if (dhash0 .lt. dhash)  then
c        time reversed representation is smaller, so return
c        that version of the path
         dhash = dhash0
         do 300  i = 1, npat
            ipat(i) = ipat0(i)
            rx(i)   = rx0(i)
            ry(i)   = ry0(i)
            rz(i)   = rz0(i)
  300    continue
      endif

      return
      end
      subroutine totale (dval)
      implicit double precision (a-h,o-z)
      save
      common /print/ iprint
      integer*4 nstop
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)
      dc(1)=1.0d0
      do 10 i=1,np
   10 dp(i)=d(i)/dr(i)
      if (nuc.le.0) go to 30
      do 20 i=1,nuc
   20 dp(i)=d(i)*(3.0d0-dr(i)*dr(i)/(dr(nuc)*dr(nuc)))/(dr(nuc)+dr(nuc))
      dc(1)=4.0d0
   30 call somm (dr,dp,dq,dpas,dc(1),0,np)
      dc(1)=-z*dc(1)
      do 40 i=1,np
      dp(i)=d(i)*dvf(i)
      dvn(i)=d(i)*dvn(i)
   40 d(i)=d(i)*exchee(d(i),dr(i))
      dc(2)=2.0d0
      dc(3)=1.0d0
      dc(5)=2.0d0
      if (nuc.ne.0) dc(3)=4.0d0
      call somm (dr,dp,dq,dpas,dc(3),0,np)
      call somm (dr,dvn,dq,dpas,dc(5),0,np)
      call somm (dr,d,dq,dpas,dc(2),0,np)
      dc(4)=dval-dc(3)
      dval=dval-0.50d0*dc(5)-dc(2)
      dc(2)=dc(3)-dc(1)-dc(5)-dc(2)
      dc(3)=0.50d0*dc(5)
      if (iprint .ge. 5)  write(16,50) dval,dc(4),dc(3),dc(2),dc(1)
   50 format (1h0,5x,'et=',1pe14.7,5x,'ec=',1pe14.7,5x,'ee=',1pe14.7,5x,
     1 'ex=',1pe14.7,5x,'en=',1pe14.7)
      return
      end
      subroutine feff_trap (x, y, n, sum)
      implicit double precision (a-h, o-z)

c     Trapeziodal integration of y(x), result in sum

      dimension x(n), y(n)

      sum = y(1) * (x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * (x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * (x(n) - x(n-1))
      sum = sum/2

      return
      end
      subroutine wphase (nph, em, eref, lmax, ne, ph)

c     Writes phase data to file PHASExx.DAT for each shell

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      complex*16 eref(nex)
      complex*16 ph(nex,ltot+1,0:nphx)
      dimension em(nex)
      dimension lmax(0:nphx)
      character*30  fname

c     Dump phase data, eref and complex phase for each shell
      do 260  iph = 0, nph
c        prepare file for shell's phase data
         write(fname,242)  iph
  242    format('phase', i2.2, '.dat')
         open (unit=1, file=trim(header)//fname,
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//fname, 'wphase')
         call wthead  (1)
c        write out unique pot and lmax
         write(1,244)   iph, lmax(iph), ne
  244    format (1x, 3i4, '   unique pot,  lmax, ne ')
c        for each energy
c        ie, em, eref, p=sqrt(em-eref)
c        ph array to ltot+1, 5 values per line
         do 250  ie = 1, ne
            xp = sqrt(em(ie) - eref(ie))
            write(1,246)  ie, em(ie), eref(ie), sqrt(em(ie)-eref(ie))
  246       format ('   ie        energy      re(eref)',
     1              '      im(eref)',
     2              '         re(p)         im(p)', /,
     3              1x, i4, 1p, 5e14.6)
            write(1,248)  (ph(ie,ll,iph), ll=1,lmax(iph)+1)
  248       format (1x, 1p, 4e14.6)
  250    continue
         close(unit=1)
  260 continue

      return
      end
      subroutine wpot (nph, edens, ifrph, imt, inrm,
     1                 rho, vclap, vcoul, vtot)

c     Writes potentials to file name POTxx.DAT for each unique pot.

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      dimension ifrph(0:nphx)
      dimension rho(251,0:nfrx)
      dimension vcoul(251,0:nfrx)
      dimension edens(nrptx,0:nphx)
      dimension vclap(nrptx,0:nphx)
      dimension vtot (nrptx,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)

      character*30 fname

c     note units --
c     potentials in rydbergs, so that v * 13.6 -> eV
c     density in #/(bohr)**3, so rho * e / (.529)**3 -> e/(Ang)**3

      do 180  iph = 0, nph
         ifr = ifrph(iph)
c        prepare file for unique potential data
         write(fname,172)  iph
  172    format('pot', i2.2, '.dat')
         open (unit=1, file=trim(header)//fname,
     >         status='unknown', iostat=ios)
         call chopen (ios, trim(header)//fname, 'wpot')
         call wthead(1)
         write(1,173)  iph, imt(iph), inrm(iph)
  173    format (1x, 3i4, '  Unique potential, I_mt, I_norman.',
     1          '    Following data in atomic units.')
         write(1,*) ' ifr ', ifr
         write(1,174)
  174    format ('   i      r         vcoul        rho',
     1           '     ovrlp vcoul  ovrlp vtot  ovrlp rho')
         do 178  i = 1, nrptx
            write(1,176) i, rr(i), vcoul(i,ifr), rho(i,ifr)/(4*pi),
     1                vclap(i,iph), vtot(i,iph), edens(i,iph)/(4*pi)
  176       format (1x, i3, 1p, 6e12.4)
  178    continue
         close(unit=1)
  180 continue

      return
      end
      subroutine xcpot (iph, ie, nr, index, ifirst, jri,
     1                  em, xmu, vi0, rs0, gamach,
     2                  vr, densty,
     3                  eref, v,
     4                  vxcrmu, vxcimu)

      implicit double precision (a-h, o-z)

      character*72 header
      common /header_common/ header


c     INPUT
c     iph, ie used only for debug and labels.
c     nr          number of points in current Loucks r-grid
c     index       0  Hedin-Lunqvist + const real & imag part
c                 1  Dirac-Hara + const real & imag part
c                 2  ground state + const real & imag part
c                 3  Dirac-Hara + HL imag part + const real & imag part
c                 4  See rdinp for comment
c     ifirst      first entry flag, set to zero before first call for
c                 each unique potential, see vxcrmu and vxcimu below
c     jri         index of first interstitial point in current
c                 Loucks r grid
c     em          current energy grid point
c     xmu         fermi level
c     vi0         const imag part to subtract from potential
c     rs0         user input density cutoff, index=4 only
c     gamach      core hole lifetime
c     vr(nr)      total potential (coulomb and gs exchange corr)
c     densty(nr)  electron density
c
c     OUTPUT
c     eref        complex energy reference for current energy
c     v(nr)       complex potential including energy dep xc
c
c     WORKSPACE
c     vxcrmu and vxcimu are calculated only on first entry for a
c     particular unique potential, re-used on subsequent entries.
c     vxcrmu(nr)  real part of xc at fermi level
c     vxcimu(nr)  imag part of xc at fermi level
c
c     This subroutine uses atomic (hartree) units for energy,
c     phase uses rydbergs.  All inputs to and outputs from xcpot are
c     in rydbergs.  (Factor of 2 to convert from one to the other.)



      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      dimension   vr(nr), densty(nr)
      complex*16  eref, v(nr)
      dimension   vxcrmu(nr), vxcimu(nr)

      complex*16  delta

c     First calculate vxc to correct the local momentum dispersion
c     relation, delta = vxc(e,k) - vxc(mu,k), and
c               p^2 = k^2 -mu + kf^2 - delta.
c     In jr theory, v(e,r) = vcoul(r) + vxc(e,r) =
c                          = vcoul(r) + vxcgs(r) + delta(e,r).

      if (index .eq. 2)  then
c        Ground state exchange, no self energy calculation
         do 10  i = 1, jri
            v(i) = vr(i)
   10    continue
      else
c        Add the self energy correction
         do 20  i = 1, jri
            rs = (3 / (4*pi*densty(i))) ** third
c           xf = 1.9191.../rs
            xf = fa / rs

c           xk2 is the local momentum squared, p^2 = k^2 - mu + kf^2,
c           k^2 represents energy measured from vacuum.
c           See formula 2.15 in Lee and Beni's paper with the last 2
c           terms neglected.  (complete reference?)
            xk2 = em + xf**2 - xmu

            if (xk2 .lt. 0)  then
               write(77,*) 'i, jri'
               write(77,*) i, jri
               write(77,*) 'rs, densty(i)'
               write(77,*) rs, densty(i)
               write(77,*) 'xf, fa'
               write(77,*) xf, fa
               write(77,*) 'em, xmu, xk2'
               write(77,*) em, xmu, xk2
               stop 'XCPOT-1'
            endif
            xk = sqrt(xk2)
            if (index .eq. 0)  call rhl(rs,xk,vxcr,vxci)
            if (index .eq. 1)  call edp(rs,xk,vi0,vxcr,vxci)
            if (index .eq. 3)  then
               call edp(rs,xk,vi0,vxcr,vxci)
               call imhl(rs,xk,vxci,icusp)
            elseif (index .eq. 4)  then
               rstmp = (1.0d0/rs**3 - 1.0d0/rs0**3) ** (-third)
               call edp(rstmp,xk,vi0,vxcr1,vxci1)
               call rhl(rs0,xk,vxcr2,vxci2)
               vxcr = vxcr1 + vxcr2
               vxci = vxci1 + vxci2
            endif

            if (ifirst .eq. 0)  then
c              vxc_mu indep of energy, calc only once
c              Calculate vxc at fermi level e = mu, j.m. 1/12/89
               xk = xf * 1.00001d0
               if (index .eq. 0) call rhl(rs,xk,vxcrmu(i),vxcimu(i))
               if (index .eq. 1) call edp(rs,xk,vi0,vxcrmu(i),vxcimu(i))
               if (index .eq. 3) then
                  call edp(rs,xk,vi0,vxcrmu(i),vxcimu(i))
                  call imhl (rs,xk,vxcimu(i),icusp)
               elseif (index .eq. 4)  then
                  rstmp = (1.0d0/rs**3 - 1.0d0/rs0**3) ** (-third)
                  call edp(rstmp,xk,vi0,vxcr1,vxci1)
                  call rhl(rs0,xk,vxcr2,vxci2)
                  vxcrmu(i) = vxcr1 + vxcr2
                  vxcimu(i) = vxci1 + vxci2
               endif
            endif

            delta = dcmplx (vxcr-vxcrmu(i), vxci-vxcimu(i))

c           Correct local momentum according to the formula
c           p^2 = k^2 - mu + kf^2 - delta.  Note that imag part
c           of delta is ignored, since xk2 is a real quantity.
            xk2 = em + xf**2 - xmu - delta
            if (xk2 .lt. 0)  then
               write(77,*) xk2, i, ie, iph, ' xk2, i, ie, iph'
               write(77,*) 'em, xf**2, xmu, delta'
               write(77,*) em, xf**2, xmu, delta
               stop 'XCPOT-2'
            endif
            xk = sqrt (xk2)

c           recalculate vxc(e,k) and vxc(mu,k) with the corrected
c           local momentum
            if (index .eq. 0)  call rhl(rs,xk,vxcr, vxci)
            if (index .eq. 1)  call edp(rs,xk,vi0,vxcr,vxci)
            if (index .eq. 3)  then
               call edp(rs,xk,vi0,vxcr,vxci)
               call imhl (rs,xk,vxci,icusp)
            elseif (index .eq. 4)  then
               rstmp = (1.0d0/rs**3 - 1.0d0/rs0**3) ** (-third)
               call edp(rstmp,xk,vi0,vxcr1,vxci1)
               call rhl(rs0,xk,vxcr2,vxci2)
               vxcr = vxcr1 + vxcr2
               vxci = vxci1 + vxci2
            endif

c           delta corrected calculated with new local momentum
            delta = dcmplx (vxcr-vxcrmu(i), vxci-vxcimu(i))

c           Note multiplication by 2 in the exchange correlation part to
c           to convert it to rydberg units.
   19       continue
            v(i) = vr(i) + 2*delta

   20    continue
      endif

c     Reference the potential with respect to mt potential, ie,
c     first interstitial point.  v(jri) = 0

c     Note that the reference does not contain the core hole lifetime
c     since the total atomic potential should have it. However in the
c     perturbation  deltav = v - vmt it cancels out.
c     ( deltav = vat - igamma - (vatmt-igamma) ).

      eref = v(jri)
      do 11  i = 1, jri
         v(i) = v(i) - eref
   11 continue

c     igamma added to the reference so that k^2 = E - Eref, where
c     Eref = Vat(mt) - igamma / 2
      eref = eref - coni * gamach / 2

c     Add const imag part
      eref = eref - coni * vi0

      ifirst = 1
      return
      end
      double precision function xx (j)
      implicit double precision (a-h, o-z)
c     x grid point at index j, x = log(r), r=exp(x)
      parameter (delta = 0.050000000000000d0)
      parameter (c88   = 8.800000000000000d0)
c     xx = -8.8 + (j-1)*0.05
      xx = -c88 + (j-1)*delta
      return
      end

      double precision function rr(j)
      implicit double precision (a-h, o-z)
c     r grid point at index j
      rr = exp (xx(j))
      return
      end

      integer function ii(r)
      implicit double precision (a-h, o-z)
c     index of grid point immediately below postion r
      parameter (delta = 0.050000000000000d0)
      parameter (c88   = 8.800000000000000d0)
c     ii = (log(r) + 8.8) / 0.05 + 1
      ii = (log(r) + c88) / delta + 1
      return
      end
      subroutine ykdir (ia,ib,nk1,nag)

      implicit double precision (a-h,o-z)
      save
      common /atomco/ den(30), dq1(30), dfl(30), ws, nqn(30), nql(30),
     1                nk(30), nmax(30), nel(30), norb, norbco
      integer*4 nstop,n1sum
      common /dira/ dv(251), dr(251), dp(251), dq(251), dpas, tets,
     1              z, nstop, nes, np, nuc
      common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30),
     1 dpc(251,30)
      common /trois/ dpno(4,30), dqno(4,30)
      dimension dpn1(4)
      dpah=exp(dpas)
      dpyk=dpas/24.0d0
      id=min0(nmax(ia)+2,nmax(ib)+2,np)
      idm1=id-1
      if (nag.ne.0) go to 30
      do 10 i=1,id
   10 dq(i)=dr(i)*(dgc(i,ia)*dgc(i,ib)+dpc(i,ia)*dpc(i,ib))
      do 20 i=1,4
      dpn1(i)=0.0d0
      do 20 j=1,i
   20 dpn1(i)=dpn1(i)+dpno(j,ia)*dpno(i+1-j,ib)+dqno(j,ia)*dqno(i+1-j,ib
     1 )
      go to 60
   30 do 40 i=1,id
   40 dq(i)=dr(i)*dgc(i,ia)*dpc(i,ib)
      do 50 i=1,4
      dpn1(i)=0.0d0
      do 50 j=1,i
   50 dpn1(i)=dpn1(i)+dpno(j,ia)*dqno(i+1-j,ib)
   60 di=dfl(ia)+dfl(ib)+nk1
      dp(1)=0.0d0
      dp(2)=0.0d0
      do 70 i=1,4
      di=di+1.0d0
      dp(1)=dp(1)+(dr(1)**di)*dpn1(i)/di
   70 dp(2)=dp(2)+(dr(2)**di)*dpn1(i)/di
      dm=dpah**(-nk1)
      dim2=-dpyk*dm*dm
      dim1=13.0d0*dpyk*dm
      di=13.0d0*dpyk
      dip1=-dpyk/dm
      do 80 i=3,idm1
   80 dp(i)=dp(i-1)*dm+dim2*dq(i-2)+dip1*dq(i+1)+dim1*dq(i-1)+di*dq(i)
      dq(id-2)=dp(id-2)
      do 90 i=idm1,np
   90 dq(i)=dq(i-1)*dm
      i=nk1+nk1+1
      dm=dm/dpah
      dim2=i*dim2/(dpah*dpah)
      dim1=i*dim1/dpah
      di=i*di
      dip1=i*dip1*dpah
      i=id-3
  100 dq(i)=dq(i+1)*dm+dim2*dp(i+2)+dip1*dp(i-1)+dim1*dp(i+1)+di*dp(i)
      i=i-1
      n1sum=i-1
      if (n1sum) 110,110,100
  110 dq(1)=dq(3)*dm*dm+8.0d0*((di*dp(1)
     >  +4.0d0*dim1*dp(2))/13.0d0-dim2*dp(3))
      return
      end


      subroutine feff6(header_in)
c     EXAFS only lite version of FEFF6 
c     see LICENSE for copying details
      implicit double precision (a-h, o-z)
      character*(*) header_in


      character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
      common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch


      parameter (nphx = 7)	!max number of unique potentials (potph)
      parameter (npotx = nphx)	!max number of unique potentials (genfmt, paths)
      parameter (nfrx = nphx)	!max number of free atom types
      parameter (novrx = 8)	!max number of overlap shells
      parameter (natx = 250)	!max number of atoms in problem
      parameter (ltot = 24)	!max number of ang mom (arrays 1:ltot+1)
      parameter (nrptx = 250)	!Loucks r grid used through overlap
      parameter (nex = 100)	!Number of energy points genfmt, etc.

      parameter (lamtot=15)	!Max number of distinct lambda's for genfmt
 				!15 handles iord 2 and exact ss
      parameter (mtot=4, ntot=2) !vary mmax and nmax independently
      parameter (legtot=9)	!matches path finder, used in GENFMT
      parameter (npatx = 8)	!max number of path atoms, used in path
				!finder, NOT in genfmt


      parameter (pi = 3.1415926535897932384626433d0)
      parameter (one = 1, zero = 0)
      parameter (third = 1.0d0/3.0d0)
      parameter (raddeg = 180.0d0 / pi)
      complex*16 coni
      parameter (coni = (0.0d0,1.0d0))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919158292677512811d0)

      parameter (bohr = 0.529177249d0, ryd  = 13.605698d0)
      parameter (alpinv = 137.03598956d0)
c     fine structure alpha
      parameter (alphfs = 1.0d0 / alpinv)
c     speed of light in louck's units (rydbergs?)
      parameter (clight = 2 * alpinv)


      parameter (ntitx = 10)
      character*79  title(ntitx)
      dimension ltit(ntitx)
      character*12 tmpstr
      character*30 fname

c     Following passed to pathfinder, which is single precision.
c     Be careful to always declare these!
      parameter (necrit=9, nbeta=40)
      real*8 fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)
      real*8 fbeta(-nbeta:nbeta,0:npotx,nex), cksp(nex)
      real*8 rmax, critpw, pcritk, pcrith
      character*6  potlbl(0:npotx)

      character*72 header
      common /header_common/ header

   10 format (1x, a)

      header = header_in

      open (unit=77,file=trim(header)//'feff.stdout',status='unknown')

      tmpstr = vfeff
      call triml (tmpstr)
      write(77,10) tmpstr
      call rdinp(mphase, mpath, mfeff, mchi, ms,
     1            ntitle, title, ltit,
     2            critcw, 
     1            ipr2, ipr3, ipr4,
     1            s02, tk, thetad, sig2g,
     1            nlegxx,
     1            rmax, critpw, pcritk, pcrith, nncrit,
     2            icsig, iorder, vrcorr, vicorr, isporb)

      do 20  i = 1, ntitle
         write(77,10) title(i)(1:ltit(i))
   20 continue

      if (mphase .eq. 1)  then
         write(77,10) 'Calculating potentials and phases...'
         call potph (isporb)
         open (unit=1, file=trim(header)//'potph.dat',
     >         status='old', iostat=ios)
         call chopen (ios, trim(header)//'potph.dat', 'feff')
         close (unit=1, status='delete')
      endif

      if (ms.eq.1  .and.  mpath.eq.1)  then

         write(77,10) 'Preparing plane wave scattering amplitudes...'
         call prcrit(ne, nncrit, ik0, cksp, fbeta, ckspc, 
     1                fbetac, potlbl)

c        Dump out fbetac for central atom and first pot
         if (ipr2 .ge. 3 .and. ipr2.ne.5)  then
            do 260  ipot = 0, 1
               do 250  ie = 1, nncrit
                  write(fname,200)  ie, ipot
  200             format ('fbeta', i1, 'p', i1, '.dat')
                  open (unit=1, file=trim(header)//fname)
                  write(1,210)  ipot, ie, ckspc(ie)
  210             format ('# ipot, ie, ckspc(ie) ', 2i5, 1pe20.6, /
     1                    '#  angle(degrees), fbeta/|p|,  fbeta')
                  do 240  ibeta = -nbeta, nbeta
                     cosb = .025 * ibeta
                     if (cosb .gt.  1)  cosb =  1
                     if (cosb .lt. -1)  cosb = -1
                     angle = acos (cosb)
                     write(1,230)  angle*raddeg, 
     1                  fbetac(ibeta,ipot,ie)/ckspc(ie),
     2                  fbetac(ibeta,ipot,ie)
  230                format (f10.4, 1p, 2e15.6)
  240             continue
                  close (unit=1)
  250          continue
  260       continue
         endif

         write(77,10) 'Searching for paths...'
         call paths(ckspc, fbetac, pcritk, pcrith, nncrit,
     1               rmax, nlegxx, ipotnn)

         write(77,10) 'Eliminating path degeneracies...'
         call pathsd(ckspc, fbetac, ne, ik0, cksp, fbeta,
     1                critpw, ipotnn, ipr2, 
     1                pcritk, pcrith, nncrit, potlbl)

         if (ipr2 .lt. 2)  then
            open (unit=1, file=trim(header)//'geom.dat', status='old')
            call chopen (ios, trim(header)//'geom.dat', 'feff')
            close (unit=1, status='delete')
         endif
      endif

      if (mfeff .eq. 1)  then
         write(77,10) 'Calculating EXAFS parameters...'
         call genfmt (ipr3, critcw, sig2g, iorder)
      endif

      if (mchi .eq. 1)  then
         write(77,10) 'Calculating chi...'
         call ff2chi (ipr4, critcw, s02, tk, thetad, icsig,
     1                vrcorr, vicorr)
      endif

      write(77,500)
  500 format (1x, 'Feff done.  Have a nice day.')

      close(77)

      return
      end


      character*2 function upperlower(string)
      implicit none
      character*2 string
      character*2 item
      character t1,t2
      integer uca,ucz,lca,lcz,shift
      uca = ichar('A')
      ucz = ichar('Z')
      lca = ichar('a')
      lcz = ichar('z')
      shift = lca - uca
      item = '  '
      t1 = string(1:1)
      t2 = string(2:2)
      if ((ichar(t1).ge.lca).and.(ichar(t1).le.lcz))
     >   t1 = char(ichar(t1)-shift)
      if ((ichar(t2).ge.uca).and.(ichar(t2).le.ucz))
     >   t2 = char(ichar(t2)+shift)
      item(1:1) = t1
      item(2:2) = t2
      upperlower = item
      return
      end


      logical function str_compare(a,b)
      implicit none
      character*(*) a
      character*(*) b
      character t1,t2
      integer la,lb,uca,ucz,lca,i,shift
      str_compare = .false.
      la = len_trim(a)
      lb = len_trim(b)
      if (la.eq.lb) then
         uca = ichar('A')
         ucz = ichar('Z')
         lca = ichar('a')
         shift = lca-uca
         do i=1,la
            t1 = a(i:i)
            t2 = b(i:i)
            if ((ichar(t1).ge.uca).and.(ichar(t1).le.ucz))
     >         t1 = char(ichar(t1)+shift)
            if ((ichar(t2).ge.uca).and.(ichar(t2).le.ucz))
     >         t2 = char(ichar(t2)+shift)
            if (t1.ne.t2) return
         end do
         str_compare = .true.
      end if

      return
      end



      subroutine feff_codeversion(version)
      implicit none
      character*(*) version

      character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
      common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch

      version = vfeff
      return
      end



      subroutine feff_serial(dict_in,outtype,dict_out,nkf,kf,chi)
      implicit none
      character*(*) dict_in
      integer outtype
      character*(*) dict_out
      integer nkf
      real*8 kf(*),chi(*)


*     **** local variables ****
      character*2 symbols(112)
      data symbols/
     $     'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     $     'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     $     'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     $     'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     $     'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     $     'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     $     'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     $     'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     $     'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     $     'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/

      logical done,found
      character*100 buf
      character*80 header
      character*80 spectroscopy

      integer      nedge
      character*80 edge

      integer nabsorber
      character*80 absorber
      character*2 item
      integer ncenter
      integer center(20)
      character*12 cnum

      real*8 rmax 

      integer nkatm
      integer katm(1000)
      integer ipot(50),zkatm(50)
      integer nion
      character*30 cnum1
      real*8  x,y,z,dist,rion(3,1000)
      integer zi,zion(1000)

      integer ind,ind0,ind1,ind2,ind3,ind4,ind4a,ind4b
      integer i,j,ii,ia,ip

*     **** external functions ****
      character*2 upperlower
      external    upperlower
      logical  str_compare
      external str_compare


c     **** parse "scratch_dir": json item****
      header = ' '
      ind = index(dict_in,"""scratch_dir"":")
      if (ind.gt.0) then
         ind2 = ind+16
         ind3 = ind+14+index(dict_in(ind+16:),"""")
         header = dict_in(ind2:ind3)
      else
         header = ' '
      end if


c     **** parse "spectroscopy": json item ****
      spectroscopy = ' '
      ind = index(dict_in,"""spectroscopy"":")
      if (ind.gt.0) then
         ind2 = ind+17
         ind3 = ind+16+index(dict_in(ind+18:),"""")
         spectroscopy = dict_in(ind2:ind3)
      else
         spectroscopy = "exafs"
      end if


c     **** parse "edge": json item ****
      edge = ' '
      ind = index(dict_in,"""edge"":")
      if (ind.gt.0) then
         ind2 = ind+9
         ind3 = ind+7+index(dict_in(ind+9:),"""")
         edge = dict_in(ind2:ind3)
      else
         edge = "k"
      end if

c     **** parse "rmax": json item ****
      ind = index(dict_in,"""rmax"":")
      if (ind.gt.0) then
         ind1 = ind + 7
         ind4a = index(dict_in(ind1:),",")
         ind4b = index(dict_in(ind1:),"}")
         if ((ind4a.lt.ind4b).and.(ind4a.gt.0)) then
            ind2 = ind1+ind4a
         else
            ind2 = ind1+ind4b
         end if
         cnum1 = trim(adjustl(dict_in(ind1:ind2-2)))
         read(cnum1,*) rmax
      else
         rmax = 10.0d0
      end if


c     **** parse "absorber": json item ****
      absorber = repeat(' ',80)
      nabsorber = 0
      ind0 = index(dict_in,"""absorber"":")
      if (ind0.gt.0) then
         ind1 = ind0 + index(dict_in(ind0:),"[")
         done = .false.
         do while (.not.done)
            ind2 = ind1 + index(dict_in(ind1:),"""")
            ind3 = ind2 + index(dict_in(ind2:),"""")-2
            item = '  '
            item = dict_in(ind2:ind3)
            absorber(1+2*nabsorber:2+2*nabsorber) = item
            nabsorber = nabsorber + 1

            ind4a = index(dict_in(ind3:),",")
            ind4b = index(dict_in(ind3:),"]")
            if ((ind4a.lt.ind4b).and.(ind4a.gt.0)) then
               ind4a = ind4a + ind3
               ind1 = ind4a
            else
               ind4b = ind4b + ind3
               ind1 = ind4b
               done = .true.
            end if
         end do
      end if


c     **** parse "center": json item ****
      ncenter = 0
      ind0 = index(dict_in,"""center"":")
      if (ind0.gt.0) then
         ind1 = ind0 + index(dict_in(ind0:),"[")
         done = .false.
         do while (.not.done)
            ind4a = index(dict_in(ind1:),",")
            ind4b = index(dict_in(ind1:),"]")
            if ((ind4a.lt.ind4b).and.(ind4a.gt.0)) then
               ind4a = ind4a + ind1
               cnum = trim(adjustl(dict_in(ind1:ind4a-2)))
               ind1 = ind4a
            else
               ind4b = ind4b + ind1
               cnum = trim(adjustl(dict_in(ind1:ind4b-2)))
               ind1 = ind4b
               done = .true.
            end if
            ncenter = ncenter + 1
            read(cnum,'(I12)') center(ncenter)
         end do
      end if
      if (ncenter.lt.1) then
         ncenter = 1
         center(1) = 1
      end if

c     **** parse "geometry": json item ****
      nion = 0
      ind0 = index(dict_in,"""geometry"":")
      if (ind0.gt.0) then
         ind0 = ind0 + index(dict_in(ind0:),"[")
         done = .false.
         do while (.not.done)
            ind0 = ind0 + index(dict_in(ind0:),"[")

            ind1 = ind0 + index(dict_in(ind0:),",")
            ind2 = ind1 + index(dict_in(ind1:),",")
            ind3 = ind2 + index(dict_in(ind2:),",")
            ind4 = ind3 + index(dict_in(ind3:),"]")

            item = '  '
            item = dict_in(ind0+1:ind1-3)
            item = upperlower(item)
            zi = -1
            do ii=1,112
               if (item.eq.symbols(ii)) zi = ii
            end do
      
            cnum1 = trim(adjustl(dict_in(ind0:ind1-2)))

            cnum1 = trim(adjustl(dict_in(ind1:ind2-2)))
            read(cnum1,*) x

            cnum1 = trim(adjustl(dict_in(ind2:ind3-2)))
            read(cnum1,*) y

            cnum1 = trim(adjustl(dict_in(ind3:ind4-2)))
            read(cnum1,*) z
           

            ind4a = index(dict_in(ind4:),",")
            ind4b = index(dict_in(ind4:),"]")
            if ((ind4a.lt.ind4b).and.(ind4a.gt.0)) then
               ind0 = ind4+ind4a
            else
               ind0 = ind4+ind4b
               done = .true.
            end if 
            nion = nion + 1
            rion(1,nion) = x
            rion(2,nion) = y
            rion(3,nion) = z
            zion(nion)   = zi
         end do
      end if

      nkatm = 0
      do ii=1,nion
         found = .false.
         do j=1,nkatm
            if (zion(ii).eq.zkatm(j)) then
               found = .true.
               ia = j
            end if
         end do
         if (found) then
            katm(ii) = ia
         else
            nkatm = nkatm + 1
            zkatm(nkatm) = zion(ii)
            katm(ii) = nkatm
         end if
      end do

      do ia=1,nkatm
         ipot(ia) = -1
      end do
      ip = 1
      do ii=1,nion
         if (ii.ne.center(1)) then
            ia = katm(ii)
            if (ipot(ia).eq.-1) then
               ipot(ia) = ip
               ip = ip + 1
            end if
         end if
      end do
     
      nedge = 1
      if (str_compare(edge,"k"))   nedge = 1
      if (str_compare(edge,"l1"))  nedge = 2
      if (str_compare(edge,"l2"))  nedge = 3
      if (str_compare(edge,"l3"))  nedge = 4
      if (str_compare(edge,"m1"))  nedge = 5
      if (str_compare(edge,"m2"))  nedge = 6
      if (str_compare(edge,"m3"))  nedge = 7
      if (str_compare(edge,"m4"))  nedge = 8
      if (str_compare(edge,"m5"))  nedge = 9


      open (unit=76,file=trim(header)//'feff.inp',status='unknown')
      write(76,'("TITLE ...")')
      write(76,*) 
      write(76,'("HOLE ",I2," 1.0")') nedge
      write(76,*) 
      write(76,*) "*   mphase,mpath,mfeff,mchi"
      write(76,'("CONTROL  1   1   1   1")')
      write(76,'("PRINT    1   0   0   0")')
      write(76,*) 
      write(76,'("RMAX  ",F10.3)') rmax
      write(76,*) 
      if (nkatm.gt.0) then
         write(76,'("POTENTIALS")')
         write(76,*) "*  ipot       Z    element"
         write(76,100) 0,zion(center(1)),symbols(zion(center(1)))
         do ia=1,nkatm
            if (ipot(ia).ne.-1)
     >         write(76,100) ipot(ia), zkatm(ia),symbols(zkatm(ia))
         end do
  100 format(I8,I8,9x,A2)
      end if
      if (nion.gt.0) then
         write(76,*) 
         write(76,'("ATOMS")')
         write(76,200) 'x','y','z','ipot','tag','distance'
         do ii=1,nion
            ia = katm(ii)
            ip = ipot(ia)
            if (ii.eq.center(1)) ip = 0
            x = rion(1,ii)-rion(1,center(1))
            y = rion(2,ii)-rion(2,center(1))
            z = rion(3,ii)-rion(3,center(1))
            dist = sqrt(x**2 + y**2 + z**2)
            write(76,201) rion(1,ii),rion(2,ii),rion(3,ii),
     >                   ip,symbols(zion(ii)),dist
         end do
 200  format('*',A19,2A20,A5,2x,A4,A20)
 201  format(3E20.9,I5,2x,A4,E20.9)
      end if
      close(76)

      !*** call feff6 ****
      call feff6(header)

*     **** just cat "chi.dat" to dict_out ****
      if (outtype.eq.1) then
         dict_out  = " "
         open(15,file=trim(header)//'chi.dat')
         do
            read(15,'(A)',end=311,err=311) buf
            !dict_out = trim(dict_out)//trim(buf)//NEW_LINE('A')
         end do
 311     close(15)
      else
         nkf = 0
         open(15,file=trim(header)//'chi.dat')
         do
            read(15,'(f10.4,3e13.6)',end=411,err=410) x,y,z,dist
            nkf = nkf+1
            kf(nkf)  = x
            chi(nkf) = y
 410        continue
         end do
 411     close(15)
      end if

      return
      end




      subroutine feff_fortran(header,spectroscopy,absorption,edge,
     >                        center,rmax,e0,s0,
     >                        nkatm,katm,zkatm,nion,zion,rion,
     >                        nohydrogen,
     >                        nkf,kf,chi)
      implicit none
      character*(*) header
      character*(*) spectroscopy
      character*(*) absorption
      character*(*) edge
      integer center
      real*8  rmax,e0,s0
      integer nkatm
      integer katm(*)
      integer zkatm(*)
      integer nion
      integer zion(*)
      real*8  rion(3,*)
      logical nohydrogen

      integer nkf
      real*8 kf(*),chi(*)

*     **** local variables ****
      real*8 autoang
      parameter (autoang = 0.529177d0)

      character*2 symbols(112)
      data symbols/
     $     'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',
     $     'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',
     $     'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     $     'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     $     'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     $     'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     $     'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     $     'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     $     'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     $     'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     $     'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
     $     'Rg', 'Cn'/

      logical done,found
      integer nedge
      integer ipot(50)
      real*8  x,y,z,dist
      integer i,j,ii,ia,ip

*     **** external functions ****
      character*2 upperlower
      external    upperlower
      logical  str_compare
      external str_compare


      if (rmax.le.0.0d0) rmax   = 10.0d0
      if (center.lt.1)   center = 1
     
      if (str_compare(edge,"k"))   nedge = 1
      if (str_compare(edge,"l1"))  nedge = 2
      if (str_compare(edge,"l2"))  nedge = 3
      if (str_compare(edge,"l3"))  nedge = 4
      if (str_compare(edge,"m1"))  nedge = 5
      if (str_compare(edge,"m2"))  nedge = 6
      if (str_compare(edge,"m3"))  nedge = 7
      if (str_compare(edge,"m4"))  nedge = 8
      if (str_compare(edge,"m5"))  nedge = 9


      do ia=1,nkatm
         ipot(ia) = -1
      end do

      if (nohydrogen) then
         do ia=1,nkatm
            if (zkatm(ia).eq.1) ipot(ia) = -2
         end do 
      end if

      ip = 1
      do ii=1,nion
         if (ii.ne.center) then
            ia = katm(ii)
            if (ipot(ia).eq.-1) then
               ipot(ia) = ip
               ip = ip + 1
            end if
         end if
      end do

      open (unit=76,file=trim(header)//'feff.inp',status='unknown')
      write(76,'("TITLE ...")')
      write(76,*) 
      write(76,'("HOLE ",I2,F10.3)') nedge,s0
      write(76,*) 
      write(76,*) "*   mphase,mpath,mfeff,mchi"
      write(76,'("CONTROL  1   1   1   1")')
      write(76,'("PRINT    1   0   0   0")')
      write(76,*) 
      write(76,'("RMAX  ",F10.3)') rmax
      if (dabs(e0).gt.1.0e-3) write(76,'("CORRECTIONS  ",F10.3)') e0
      write(76,*) 
      if (nkatm.gt.0) then
         write(76,'("POTENTIALS")')
         write(76,*) "*  ipot       Z    element"
         write(76,100) 0,zion(center),symbols(zion(center))
         do ia=1,nkatm
            if (ipot(ia).gt.-1)
     >         write(76,100) ipot(ia), zkatm(ia),symbols(zkatm(ia))
         end do
  100 format(I8,I8,9x,A2)
      end if
      if (nion.gt.0) then
         write(76,*) 
         write(76,'("ATOMS")')
         write(76,200) 'x','y','z','ipot','tag','distance'
         do ii=1,nion
            ia = katm(ii)
            ip = ipot(ia)
            if (ii.eq.center) ip = 0
            x = rion(1,ii)-rion(1,center)
            y = rion(2,ii)-rion(2,center)
            z = rion(3,ii)-rion(3,center)
            dist = sqrt(x**2 + y**2 + z**2)
            if (ip.gt.-1) then
               write(76,201) rion(1,ii)*autoang,
     >                       rion(2,ii)*autoang,
     >                       rion(3,ii)*autoang,
     >                       ip,symbols(zion(ii)),dist*autoang
            end if
         end do
 200  format('*',A19,2A20,A5,2x,A4,A20)
 201  format(3F20.9,I5,2x,A4,E20.9)
      end if
      close(76)

      !*** call feff6 ****
      call feff6(header)

*     **** read chi.dat to generate nkf, kf, and chi ****
      nkf = 0
      open(unit=15,file=trim(header)//'chi.dat')
      do
         read(15,*,end=411,err=410) x,y,z,dist
         nkf = nkf+1
         kf(nkf)  = x
         chi(nkf) = y
 410     continue
      end do
 411  close(15)

      return
      end

