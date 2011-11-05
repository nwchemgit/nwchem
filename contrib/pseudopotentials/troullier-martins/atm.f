      program atm
c               ______________
c              /              \
c             /   Main - ATM   \
c             \     START      /
c              \______________/
c                      |
c                ______V_____
c               |            |          __________          _______      ______
c               |   Zesec    |         |          |        |       |    /      \
c               |____________|    |--->|  Charge  |-error->|  Ext  |-->|  STOP  |
c                      |          |    |__________|        |_______|    \______/
c                ______V______    |
c               /             /<--|   __________
c              /             /       |          |
c             /    Input    /<------>|  Zedate  |
c            /             /         |__________|
c           /_____________/------|         _______       _____
c                    |           |        |       |     /      \
c   |---->---------->|           |-error->|  Ext  |--->|  STOP  |
c   |                |                    |_______|     \______/
c   |               _V_                                        ____________
c   |             /     \        _________      _______       /            \
c   |            / More  \      |         |    |       |     /  Main - ATM  \
c   |           <         >-no->|  Zesec  |--->|  Ext  |--->|                |
c   A            \ Data? /      |_________|    |_______|     \     STOP     /
c   |             \ ___ /                                     \____________/
c   |                | yes
c   |               _V_
c   |             /     \
c   A            /Config-\
c   |           < uration >yes->|
c   |            \ Test? /      |
c   |             \ ___ /       |
c   |                | no       |
c   |           _____V_____     |
c   A          |           |    |
c   |          |  Vionic   |    |        __________
c   |          |___________|    |       |          |
c   |                |          V  |--->|  Splift  |
c   |                |<---------|  |    |__________|
c   |           _____V_____        |
c   |          |           |<------|   _______       ______
c   A          |           |          |       |     /      \
c   |          |  Velect   |--error-->|  Ext  |--->|  STOP  |
c   |          |           |          |_______|     \______/
c   |          |___________|<------|     _________
c   |                |             |    |         |
c   |                |             |--->|  Spliq  |
c   A   |--->------->|                  |_________|
c   |   A            |                                   __________
c   |   |           _V_                                 |          |
c   |   |         / 1'st\          __________     |---->|  Tridib  |
c   |   |        /or 2'nd\        |          |<---|     |__________|
c   |   |       <  Itera- >--yes->|  Dsolv1  |           __________
c   |   |        \ tion? /        |__________|<---|     |          |
c   A   |         \ ___ /               |         |---->|  Tinvit  |
c   |   A            |no                |               |__________|
c   |   |            |                  V
c   |   |            |                  |--->------>----->----->---->----->-----
c   |   |            V               ___
c   |   |            |             /     \          __________         _______
c   |   |            |            / Rela- \        |          |       |       |
c   |   |       _____V____   |--><  tivis- ><-yes->|  Difrel  |-error>|  Ext  |
c   |   A      |          |<-|    \  tic? /        |__________|       |_______|
c   |   |      |  Dsolv2  |        \ ___ /                                |
c   A   |      |__________|<-|        A                                 __V___
c   |   |            |       |        | no                             /      \
c   |   |            |       |        |                               |  STOP  |
c   |   |            |       |        |                                \______/
c   |   |            |       |   _____V____          _______      ______
c   |   A            V       |  |          |        |       |    /      \
c   |   |            |       |  |  Difnrl  |-error->|  Ext  |-->|  STOP  |
c   |   |            |       |  |__________|        |_______|    \______/
c   A   |            |       |      ___
c   |   |            |       |    /     \          _________
c   |   |            V       |   / Con-  \        |         |
c   |   |            |       |-><   verg  ><-yes->|  Orban  |
c   |   A            |           \  ed?  /        |_________|
c   |   |            |            \ ___ /
c   |   |            |
c   |   |            |<------<-----<-------<-------<-------<------<------<-----<
c   |   |            |                   __________
c   |   |            |                  |          |
c   |   A            |             |--->|  Splift  |
c   |   |       _____V_____        |    |__________|
c   |   |      |           |<------|   _______       ______
c   A   |      |           |          |       |     /      \
c   |   |      |  Velect   |--error-->|  Ext  |--->|  STOP  |
c   |   |      |           |          |_______|     \______/
c   |   A      |___________|<------|     _________
c   |   |            |             |    |         |
c   |   |            |             |--->|  Spliq  |
c   |   |            |                  |_________|
c   A   |           _V_                            ___
c   |   A         /     \       __________       /Val- \         __________
c   |   |        / Con-  \     |          |     / ence  \       |          |
c   |   |       <   verg  >--->|  Etotal  |--->< Modify? >-yes->|  Vionic  |
c   |   |        \  ed?  /     |__________|     \       /       |__________|
c   |   |         \ ___ /                        \ ___ /              |
c   A   A            | no                           |no               |
c   |   |        ____V____                          V                 V
c   |   |       |         |            |<------------<--------------<-|
c   |   |       |  Dmixp  |            |           __________________________
c   |   |       |_________|           _V_         |              0)Pseudo    |
c   |   |            |              /     \       |              1)Pseudk    |
c   |   A           _V_            /Pseudo \      |  Pseudo-     2)Pseudt    |
c   |   |         /Pass \         < Generate>-yes>|  Potential   3)Pseudv    |
c   A   |        / Max.  \         \   ?   /      |  Generation  4)Datout    |
c   |   |<--no--< Itera-  >         \ ___ /       |  Block       5)Pseudb    |
c   |            \ tion? /             |no        |              6)Pseud2    |
c   |             \ ___ /              V          |__________________________|
c   |                | yes             |---<----------<-------|
c   |             ___V___             _V_
c   A            |       |          /     \        __________
c   |            |  Ext  |         /Config-\      |          |
c   |            |_______|        < uration >-yes>|  Prdiff  |
c   |                |             \ Test? /      |__________|
c   |             ___V__            \ ___ /            |
c   |            /      \              |no             |
c   |           |  STOP  |             |               V
c   A            \______/              |<---------<----|
c   |                                  V
c   |---------<------------<-----------|
c
c
c  *************************************************************
c  *     Program for atomic calculations                       *
c  *     Copyright Norman J. Troullier Jr &                    *
c  *     Jose Luis Martins                                     *
c  *     Written by Norman J. Troullier Jr., Sept 89           *
c  *     while at U of Minn, from a Sverre Froyen              *
c  *     UC Berkeley code.  Program input/output is            *
c  *     compatible with earlier Berkeley versions.            *
c  *                                                           *
c  *     Send comments/suggestions/bug reports to:             *
c  *     troullie@csfsa.cs.umn.edu                             *
c  *     612 625-0392                                          *
c  *                                                           *
c  *     Version 5.06, Dated Oct. 19, 1990                     *
c  *                                                           *
c  *************************************************************
c
c    Some parameters are set inside the program,
c  the most important ones are:
c  1)the tolerance for selfconsistency in the screening
c    potential (set in the main-atm program-named tol),
c  2)the accuracy in the eigenvalue for a given potential
c    (set in difnrl-named tol or difrel-named tol),
c  3)the dimensions of the work space used: nrmax,
c    norbmx, lmax(needs to be set only in main-atm),
c  4)the machine precision - MACHEP, for use in the
c    eispack group of fuctions: tinvit, and tridib.
c    (The current value is ok for this application.)
c  5)the machine precision exp(-2*expzer), set in difnrl
c    and difrel for the predictor-corrector methods
c    (both named expzer),
c
c    For traditional applications the default values
c  should be enough for 1-4.  For 5, expzer should be
c  found for the system in use.
c    NOTE: that for most VAX systems the value of expzer
c  can be very very small !!
c
c    The subroutine orban is called once for each orbital
c  after selfconsistency has been achieved.
c  Any additional analysis of the orbitals should therefore
c  be placed in orban.  Note that arrays ar and br have
c  different meaning for non-relativistic (ar-wave,
c  br-d(wave)/dj) and relativistic (ar-major, br-minor)
c  calculations.
c
c    There are six ways to generate the pseudopotential :
c  ikerk = 6 Improved Troullier and Martins
c  ikerk = 5 Bachelet, Hamann, and Schluter
c  ikerk = 4 generates data file another pseudopotential
c   generation program.
c  ikerk = 3 Vanderbilt
c  ikerk = 2 Troullier and Martins
c  ikerk = 1 Kerker
c  ikerk = 0 Hamann Schluter and Chiang
c
c      This main - atm routine has had extremly major
c    modifications with respect to the Berkeley version.
c    However, all input and output files are still compatible
c    with earlier Berkeley versions of this program.
c
c    1)Machine dependent timing calls were placed
c      in the program so it could be used as a machine
c      perfomance indicatior.  The user will either have
c      to change these calls for his machine or
c      comment them out.  Corresponding subroutines
c      are included for the Apollo, Cray, Sun, and Vax
c      computers.
c    2)The plot.dat file is now opened as a formatted file,
c      this is user/plot method dependent.  The atom.job
c      file is no longer used.  Note that the Apollo
c      system does not use standard Fortran methods to
c      open unformatted files.
c    3)The charge density startup is scaled with
c      an empirical formula that reduces the
c      number of iterations needed for the screening
c      potential convergence.
c    4)The screening potential mixing parameter is
c      an empirical function of the nuclear charge.
c      Larger atoms require a slower convergence
c      then smaller atoms.
c    5)The screening potential is intially mixed with a
c      percentage of old and new for the first itsm
c      iterations. This brings the potential to a stable
c      region after which an Anderson's extrapolation scheme
c      is used.
c    6)The files pseudo.dat and plot.dat files are closed
c      and new ones are opened before the generation of a
c      pseudopotential.  This allows the generation of
c      multiple pseudopotentials(up to 99).
c    7)The pseudopotentail generation scheme of Troullier
c      and Martins is added - pseudt.  The pseudopotential
c      scheme of Vanderbilt has been added - pseudv.
c      The improved pseudopotential scheme of Troullier
c      and Martins has been added - pseud2.
c      The datout routine generates a data file for use
c      in external pseudopotential generation programs.
c      The user may wish to alter for his own use or ignore.
c    8)Only the major modifications(not programming style)
c      to the algorithms of each subroutine are commented
c      in that routine.
c    9)Cray(and other machines) conversions are indicated
c      at the begining of each routine.
c   10)The difrel and difnrl allow for the calculation of
c      a nonbound state(zero eigenvalue).  These states
c      can only be used in the pseudt, pseudk and
c      pseud2 generation
c      routines.  The pseudo, pseudb and pseudv will fail with
c      a zero eigenvalue due to the generation method.
c      The user should be very careful in using a nonbound
c      state, and should always  compare the resulting pseudopotential
c      to a bound state pseudopotential calculation.
c   11)What comes free comes with no guarantee!!
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Make sure the 2 open(unit=1) statements are
c  ###      a non-recl (non-Apollo type) format.
c  ###    3)Switch the 2 double precision parameters
c  ###      to single precision parameter statements.
c  ###  Cray conversions
c  njtj
c
c  tolerance for self-consistency
c
      implicit double precision (a-h,o-z)
      parameter (tol=1.E-8)
c
      parameter (zero=0.0,one=1.0)
c
      parameter(lmax=5,nrmax=1000,norbmx=40)
c
c
      dimension r(nrmax),rab(nrmax),no(norbmx),lo(norbmx),
     1 so(norbmx),zo(norbmx),cdd(nrmax),cdu(nrmax),cdc(nrmax),
     2 viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     3 vod(nrmax),vou(nrmax),vn1d(nrmax),vn1u(nrmax),
     4 vn11d(nrmax),vn11u(nrmax),vn2d(nrmax),vn2u(nrmax),
     5 vn22d(nrmax),vn22u(nrmax),ev(norbmx),evi(norbmx),ek(norbmx),
     6 ep(norbmx),wk1(nrmax),wk2(nrmax),wk3(nrmax),wk4(nrmax),
     7 wk5(nrmax),wk6(nrmax),wk7(nrmax),wk8(nrmax),wk9(nrmax),
     8 wkb(6*nrmax)
c
      dimension econf(100),etot(10)
c
      character*1 ispp
      character*2 naold,icold,icorr,nameat
      character*10 plotfile
      character*12 pseudofile
c
c      CALL DROPFILE(0)
c  njtj  ***  machine call  ***
c    Call to machine dependent cpu time routine.
c    User may want to comment out timing calls -
c    here and at exit of main - atm
c
      CALL ZESEC(t1)
c
c  njtj  ***  machine call  ***
c
c  Startup values for doing multiple input data sets.
c
      naold = '  '
      icold = '  '
      zsold = zero
      nconf = 0
c
c      open files
c
      isize = 8*nrmax
c
c  Note that the open(unit=1,...,recl..) is needed
c  by the Apollo systems.  For other systems the
c  Cray statements(standard Fortran 77) should be
c  the ones used.
c
          open(unit=1,file='pseudo.dat',form='unformatted',
     1 status='unknown')
c
c  njtj  ***  modification  start ***
c    The plot.dat file is now opened as a formatted file.
c    This is user/plot method dependent.  The atom.job
c    file is no longer used.
c
      open(unit=3,file='plot.dat',status='new',form='formatted')
      open(unit=5,file='atom.dat',status='old',form='formatted')
      open(unit=6,file='atom.out',status='new',form='formatted')
c
c  njtj  ***  modification end  ***
c
c   Start of loop over configuration.
c   Read the input data.
c

 20   nr = nrmax
      norb = norbmx
      call input(itype,ikerk,icorr,ispp,zsh,rsh,
     1 nr,a,b,r,rab,nameat,norb,ncore,no,lo,so,zo,
     2 znuc,zel,evi)
c
c  njtj  ***  machine call  ***
c  Stop - no more data in input file,
c  Find time taken for total calculation.
c  second - machine dependent routine
c
      if (itype .lt. 0) then
        CALL ZESEC(t2)
        write(6,2000)t2-t1
 2000 format(//,' The total time for the calculation is ',
     1 f12.5,' seconds')
        call ext(0)
      endif
c
c  njtj  *** machine call ***
c
c  Jump charge density 'set up' and ionic data input if
c  configuration test.
c
      itsm=znuc/9+3
      if (zsold .eq. zsh .and. naold .eq. nameat
     1 .and. itype .lt. 1 ) then
      else
        if (itype .lt. 4) then
c
c  Set up the initial charge density.
c  cdd and cdu  =  (4 pi r**2 rho(r))/2
c
c  njtj  ***  modification  ***
c    The charge density setup (aa) is scaled with
c    an empirical formula that reduces the
c    number of iterations needed for the screening
c    potential convergence.
c
          aa = sqrt(sqrt(znuc))/2+one
          do 30 i=1,nr
            cdd(i) = zel*aa**3*exp(-aa*r(i))*r(i)**2/4
            cdu(i) = cdd(i)
 30       continue
        endif
c
c  njtj  ***  modification end  ***
c
c  set up ionic potentials
c
        call vionic(ispp,itype,icorr,ifcore,zsh,rsh,
     1   lmax,nr,a,b,r,rab,nameat,ncore,znuc,
     2   cdd,cdu,cdc,viod,viou)
      endif
c
c   Set up the electronic potential.
c
      call velect(0,0,icorr,ispp,ifcore,
     1 nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,wk1,wk2,
     2 wk3,wk4,wk5,wkb)
c
      do 50 i=1,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 50   continue
c
c   Start the iteration loop for electronic convergence.
c
      iconv = 0
      icon2 = 0
      maxit = 100
c
c  njtj  ***  modification start  ***
c    The screening potential mixing parameter is
c    an empirical function of the nuclear charge.
c    Larger atoms require a slower convergence
c    then smaller atoms.
c
      xmixo = one/log(znuc+7*one)
c
c  njtj  ***  modifications end  ***
c
      do 100 iter=1,maxit
        if (iter .eq. maxit) iconv=1
c
c  compute orbitals
c
        if (icon2 .lt. 2) then
          call dsolv1(lmax,nr,a,b,r,rab,norb,ncore,
     1     no,lo,so,zo,cdd,cdu,viod,viou,vid,viu,ev,
     2     wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wkb)
        else
          call dsolv2(iter,iconv,ispp,ifcore,lmax,nr,
     1     a,b,r,rab,norb,ncore,no,lo,so,zo,znuc,cdd,
     2     cdu,cdc,viod,viou,vid,viu,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,evi)
        endif
c
c  set up output electronic potential from charge density
c
        call velect(iter,iconv,icorr,ispp,ifcore,
     1   nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,wk1,wk2,
     2   wk3,wk4,wk5,wkb)
c
c  check for convergence
c
        if (iconv .gt. 0) goto 120
        dvmax = zero
        do 60 i=1,nr
          dv = (vod(i)-vid(i))/(1.D0+vod(i)+vou(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
          dv = (vou(i)-viu(i))/(1.D0+vou(i)+vod(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
 60     continue
        icon2 = icon2+1
        if (dvmax .le. tol) iconv=1
c
c  Mix the input and output electronic potentials.
c
c  njtj  ***  modification  start  ***
c    The screening potential is initially mixed with a
c    percentage of old and new for itsm iterations.
c    This brings the potential to a stable region
c    after which an Anderson's extrapolation scheme
c    is used.
c
        if (iter .lt. itsm) then
          iiter=2
        else
          iiter=iter-itsm+3
        endif
        call dmixp(vod,vid,xmixo,iiter,3,nr,wk1,wk2,
     1   vn1d,vn11d,vn2d,vn22d)
        call dmixp(vou,viu,xmixo,iiter,3,nr,wk1,wk2,
     1   vn1u,vn11u,vn2u,vn22u)
 100  continue
c
c   End of iteration of electronic convergence loop.
c
      write(6,110) dvmax,xmixo
 110  format(/,' potential not converged - dvmax =',e10.4,
     1 ' xmixo =',f5.3)
      call ext(1)
c
c  njtj  ***  modification end  ***
c
c  Find the total energy.
c
 120  write(6,121)icon2
 121  format(/,'Total number of iterations needed for',
     1 ' electron screening potential is ',i2,/)
      call etotal(itype,zsh,nameat,norb,
     1 no,lo,so,zo,etot,ev,ek,ep)
c
c   Replace the valence charge density.
c
      if (itype .eq. 5) call vionic(ispp,6,icorr,
     1 ifcore,zsh,rsh,lmax,nr,a,b,r,rab,nameat,
     2 ncore,znuc,cdd,cdu,cdc,viod,viou)
c
c  Pseudopotential generation.
c
c  njtj  ***  modification  start  ***
c    Current pseudo.dat and plot.dat files are closed
c    and new ones are opened.  This allows the
c    generation of multiple pseudopotentials(up to 99).
c
      if (itype .ge.1 .and. itype .le. 3) then
        if (ikerk .ne. 4 ) then
          close(unit=1)
          close(unit=3)
          if (nconf .le. 8) then
            write(pseudofile,8000)nconf+1
            write(plotfile,8002)nconf+1
          else
            write(pseudofile,8001)nconf+1
            write(plotfile,8003)nconf+1
          endif
          write(6,8004)nconf+1
 8000 format('pseudo.dat0',i1)
 8001 format('pseudo.dat',i2)
 8002 format('plot.dat0',i1)
 8003 format('plot.dat',i2)
 8004 format(//,' Pseudopotential generation file number ',i2)

c
c  Note that the open(unit17,...,recl..) is needed
c  by the Apollo systems.  For other systems the
c  Cray statements(standard Fortran 77) should be the
c  ones used.
c
           open(unit=1,file=pseudofile,form='unformatted')
          open(unit=3,file=plotfile,status='new',
     1     form='formatted')
        endif
c
c  njtj  ***  modification  end  ***
c
c  njtj  ***  modification start  ***
c    The pseudopotentail generation scheme of Troullier
c    and Martins is added - pseudt.  The pseudopotential
c    scheme of Vanderbilt is added - pseudv.  The
c    pseudopotential scheme of BHS is added - pseudb.
c    The improved pseudopotential scheme of Troullier
c    and Martins is added - pseud2.  The
c    dataout routine generates a data file for use in
c    external pseudopotential generation programs.
c    The user may wish to alter for their own use or ignore.
c
        if(ikerk.eq.0) then
          call pseudo(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,vn2d,vn2u,
     4     vn11d,wkb,evi)
        elseif (ikerk .eq. 1) then
          call pseudk(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,wkb,evi)
        elseif (ikerk .eq. 2) then
          call pseudt(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,wkb,evi)
        elseif (ikerk .eq. 3) then
          call pseudv(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,vn2d,vn2u,
     4     vn11d,wkb,evi)
        elseif (ikerk .eq. 4) then
          call datout(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdc,
     2     viod,viou,vid,viu,ev)
        elseif (ikerk .eq. 5) then
          call pseudb(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,vn2d,vn2u,
     4     vn11d,wkb,evi)
        elseif (ikerk .eq. 6) then
          call pseud2(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1     nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2     viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3     wk3,wk4,wk5,wk6,wk7,wk8,wk9,vn1d,vn1u,wkb,evi)
        endif
      endif
c
c  njtj   ***  modification end  ***
c
c  printout difference from first configuration
c
      nconf = nconf + 1
      econf(nconf) = etot(10)
      if(naold .eq. nameat .and. icold .eq. icorr .and. nconf .ne. 1
     1 .and. (itype .lt. 1 .or. itype .gt. 3)) then
        call prdiff(nconf,econf)
        write(6,130) etot(10)-econf(1)
      endif
      write(6,135)
 130  format(//,' excitation energy         =',f18.8,/)
 135  format(//,60('%'),//)
      naold = nameat
      icold = icorr
      zsold = zsh
c
c   End loop of configuration.
c
      goto 20
      end
C
C
C
       double precision function charge(name)
Cray        function charge(name)
c
c      function determines the nuclear charge of an element
c
c  njtj  ***  modifications  ***
c    All elements from H to Lr are included
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Switch double precision function to function.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
        parameter (one=1.0)
c
       character*2 name
c
       if (name .eq. 'H ' .or. name .eq. ' H') then
         charge = 1*one
       elseif (name .eq. 'He') then
         charge = 2*one
       elseif (name .eq. 'Li') then
         charge = 3*one
       elseif (name .eq. 'Be') then
         charge = 4*one
       elseif (name .eq. 'B ' .or. name .eq. ' B') then
         charge = 5*one
       elseif (name .eq. 'C ' .or. name .eq. ' C') then
         charge = 6*one
       elseif (name .eq. 'N ' .or. name .eq. ' N') then
         charge = 7*one
       elseif (name .eq. 'O ' .or. name .eq. ' O') then
         charge = 8*one
       elseif (name .eq. 'F ' .or. name .eq. ' F') then
         charge = 9*one
       elseif (name .eq. 'Ne') then
         charge = 10*one
       elseif (name .eq. 'Na') then
         charge = 11*one
       elseif (name .eq. 'Mg') then
         charge = 12*one
       elseif (name .eq. 'Al') then
         charge = 13*one
       elseif (name .eq. 'Si') then
         charge = 14*one
       elseif (name .eq. 'P ' .or. name .eq. ' P') then
         charge = 15*one
       elseif (name .eq. 'S ' .or. name .eq. ' S') then
         charge = 16*one
       elseif (name .eq. 'Cl') then
         charge = 17*one
       elseif (name .eq. 'Ar') then
         charge = 18*one
       elseif (name .eq. 'K ' .or. name .eq. ' K') then
         charge = 19*one
       elseif (name .eq. 'Ca') then
         charge = 20*one
       elseif (name .eq. 'Sc') then
         charge = 21*one
       elseif (name .eq. 'Ti') then
         charge = 22*one
       elseif (name .eq. 'V ' .or. name .eq. ' V') then
         charge = 23*one
       elseif (name .eq. 'Cr') then
         charge = 24*one
       elseif (name .eq. 'Mn') then
         charge = 25*one
       elseif (name .eq. 'Fe') then
         charge = 26*one
       elseif (name .eq. 'Co') then
         charge = 27*one
       elseif (name .eq. 'Ni') then
         charge = 28*one
       elseif (name .eq. 'Cu') then
         charge = 29*one
       elseif (name .eq. 'Zn') then
         charge = 30*one
       elseif (name .eq. 'Ga') then
         charge = 31*one
       elseif (name .eq. 'Ge') then
         charge = 32*one
       elseif (name .eq. 'As') then
         charge = 33*one
       elseif (name .eq. 'Se') then
         charge = 34*one
       elseif (name .eq. 'Br') then
         charge = 35*one
       elseif (name .eq. 'Kr') then
         charge = 36*one
       elseif (name .eq. 'Rb') then
         charge = 37*one
       elseif (name .eq. 'Sr') then
         charge = 38*one
       elseif (name .eq. 'Y ' .or. name .eq. ' Y') then
         charge = 39*one
       elseif (name .eq. 'Zr') then
         charge = 40*one
       elseif (name .eq. 'Nb') then
         charge = 41*one
       elseif (name .eq. 'Mo') then
         charge = 42*one
       elseif (name .eq. 'Tc') then
         charge = 43*one
       elseif (name .eq. 'Ru') then
         charge = 44*one
       elseif (name .eq. 'Rh') then
         charge = 45*one
       elseif (name .eq. 'Pd') then
         charge = 46*one
       elseif (name .eq. 'Ag') then
         charge = 47*one
       elseif (name .eq. 'Cd') then
         charge = 48*one
       elseif (name .eq. 'In') then
         charge = 49*one
       elseif (name .eq. 'Sn') then
         charge = 50*one
       elseif (name .eq. 'Sb') then
         charge = 51*one
       elseif (name .eq. 'Te') then
         charge = 52*one
       elseif (name .eq. 'I ' .or. name .eq. ' I') then
         charge = 53*one
       elseif (name .eq. 'Xe') then
         charge = 54*one
       elseif (name .eq. 'Cs') then
         charge = 55*one
       elseif (name .eq. 'Ba') then
         charge = 56*one
       elseif (name .eq. 'La') then
         charge = 57*one
       elseif (name .eq. 'Ce') then
         charge = 58*one
       elseif (name .eq. 'Pr') then
         charge = 59*one
       elseif (name .eq. 'Nd') then
         charge = 60*one
       elseif (name .eq. 'Pm') then
         charge = 61*one
       elseif (name .eq. 'Sm') then
         charge = 62*one
       elseif (name .eq. 'Eu') then
         charge = 63*one
       elseif (name .eq. 'Gd') then
         charge = 64*one
       elseif (name .eq. 'Tb') then
         charge = 65*one
       elseif (name .eq. 'Dy') then
         charge = 66*one
       elseif (name .eq. 'Ho') then
         charge = 67*one
       elseif (name .eq. 'Er') then
         charge = 68*one
       elseif (name .eq. 'Tm') then
         charge = 69*one
       elseif (name .eq. 'Yb') then
         charge = 70*one
       elseif (name .eq. 'Lu') then
         charge = 71*one
       elseif (name .eq. 'Hf') then
         charge = 72*one
       elseif (name .eq. 'Ta') then
         charge = 73*one
       elseif (name .eq. 'W ' .or. name .eq. ' W') then
         charge = 74*one
       elseif (name .eq. 'Re') then
         charge = 75*one
       elseif (name .eq. 'Os') then
         charge = 76*one
       elseif (name .eq. 'Ir') then
         charge = 77*one
       elseif (name .eq. 'Pt') then
         charge = 78*one
       elseif (name .eq. 'Au') then
         charge = 79*one
       elseif (name .eq. 'Hg') then
         charge = 80*one
       elseif (name .eq. 'Tl') then
         charge = 81*one
       elseif (name .eq. 'Pb') then
         charge = 82*one
       elseif (name .eq. 'Bi') then
         charge = 83*one
       elseif (name .eq. 'Po') then
         charge = 84*one
       elseif (name .eq. 'At') then
         charge = 85*one
       elseif (name .eq. 'Rn') then
         charge = 86*one
       elseif (name .eq. 'Fr') then
         charge = 87*one
       elseif (name .eq. 'Ra') then
         charge = 88*one
       elseif (name .eq. 'Ac') then
         charge = 89*one
       elseif (name .eq. 'Th') then
         charge = 90*one
       elseif (name .eq. 'Pa') then
         charge = 91*one
       elseif (name .eq. ' U' .or. name .eq. 'U ') then
         charge = 92*one
       elseif (name .eq. 'Np') then
         charge = 93*one
       elseif (name .eq. 'Pu') then
         charge = 94*one
       elseif (name .eq. 'Am') then
         charge = 95*one
       elseif (name .eq. 'Cm') then
         charge = 96*one
       elseif (name .eq. 'Bk') then
         charge = 97*one
       elseif (name .eq. 'Cf') then
         charge = 98*one
       elseif (name .eq. 'Es') then
         charge = 99*one
       elseif (name .eq. 'Fm') then
         charge = 100*one
       elseif (name .eq. 'Md') then
         charge = 101*one
       elseif (name .eq. 'No') then
         charge = 102*one
       elseif (name .eq. 'Lr') then
         charge = 103*one
       else
         write(6,100) name
 100   format(//,'element ',a2,' unknown')
         call ext(200)
       endif
       return
       end
C
C
C
      subroutine datout(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdc,viod,viou,
     2 vid,viu,ev)
c
c  *********************************************************
c  *
c  *  njtj
c  *  The routine writes needed data to file 'datafile.dat'
c  *  for latter use in a minimization program.
c  *  Users may want to remove or modify this routine
c  *  depending on their needs.
c  *  njtj
c  *
c  ***********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Make sure the open(unit=7) is a non-recl
c  ###      (non-Apollo type) format.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdc(nr),viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),
     2 ev(norb)
c
      character*1 ispp
      character*2 icorr,nameat
c
c  Open and write out data to current file datafile.dat.
c  Note that the open(unit=7,...,recl..) is needed
c  by the Apollo systems.  For other systems the
c  Cray statements(standard Fortran 77) should be used.
c
       open (unit=7,file='datafile.dat',status='new',
     1 form='unformatted')
      write(7)itype,icorr,ispp,nr,a,b
      write(7)(r(i),i=1,nr)
      write(7)(rab(i),i=1,nr)
      write(7)lmax,nameat,norb,ncore
      write(7)(no(i),i=1,norb)
      write(7)(lo(i),i=1,norb)
      write(7)(so(i),i=1,norb)
      write(7)(zo(i),i=1,norb)
      write(7)znuc,zel
      write(7)(cdc(i),i=1,nr)
      do 1,j=1,lmax
        write(7)(viod(j,i),i=1,nr)
 1    continue
      do 2,j=1,lmax
        write(7)(viou(j,i),i=1,nr)
 2    continue
      write(7)(vid(i),i=1,nr)
      write(7)(viu(i),i=1,nr)
      write(7)(ev(i),i=1,norb)
      close (unit=7)
c
      return
      end
c
c
c
      subroutine difnrl(iter,iorb,v,ar,br,lmax,
     1 nr,a,b,r,rab,norb,no,lo,so,znuc,viod,viou,
     2 vid,viu,ev,iflag,rab2,fa,fb,evi)
c
c    difnrl integrates the Schroedinger equation
c    if finds the eigenvalue ev, the wavefunction ar
c    and the derivative br = d(ar)/dr
c
c  njtj  ***  modifications  ***
c    This routine has had major modifications.  Some
c    of the data used inside the main loop has been
c    calculated outside the main loop to reduce the number
c    of operations(uses extra array space to gain speed)
c    and are passed as work arrays form the main.
c    The predictor-corrector functions have been put
c    into a array.
c    The iflag variable was added to indicate nonconvergence
c    for other programs.  It has no use in the atom program
c    and can be removed by the user.
c    All output from the routine is compatible to
c    the Berkeley/Sverre Froyen version.
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch the double precision parameter statements
c  ###      to single precision parameter statements.
c  ###  Cray conversions
c  njtj
c
c  njtj
c  &&&  Machine dependent Parameter
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for
c  &&&    the machine in use from the list, or find
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c
      implicit double precision (a-h,o-z)
       parameter(zero=0.0,pnine=0.9,two=2.0,etol=-1.E-7)
c
c  Tolerence
c
       parameter(tol=1.E-10,five=5.0)
c
c  Integration coefficients
c
       parameter(abc1=190.1/72,abc2=-138.7/36,abc3=10.9/3,
     1 abc4=-63.7/36,abc5=25.1/72,amc0=25.1/72,amc1=32.3/36,
     2 amc2=-1.1/3,amc3=5.3/36,amc4=-1.9/72)
c
c
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr),no(norb),
     1 lo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),evi(norb)
c
c  njtj  *** start modification  ***
c    Arrays added to gain speed.
c
      dimension rabrlo(5),rlp(5),rab2(nr),fa(nr),fb(nr)
c
c  njtj  ***  end modification  ***
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
      expzer = 3.7D2
cVax      expzer = 44.D0
Cray       expzer =  2.8E3
c
c  njtj  *** major modification start  ***
c    Loop data calculated outside loop to gain speed.
c
      itmax = 100
      iflag = 0
      lp = lo(iorb)+1
      ar(1) = zero
      if (lo(iorb) .eq. 0) then
        br(1) = b*a
      else
        br(1) = zero
      endif
      do 1 j=2,nr
        ar(j) = zero
 1    continue
      do 2 j=2,nr
        br(j) =zero
 2    continue
      do 4 j=2,5
        rlp(j)=r(j)**lp
 4    continue
      do 5 j=2,5
        rabrlo(j)=rab(j)*r(j)**lo(iorb)
 5    continue
      do 6 j=1,nr
        rab2(j)=rab(j)*rab(j)
 6    continue
c
c   set underflow trap, error from Berkeley version,
c   fixed by Troy Barbee sqrt(expzer) should be expzer/2
c   4/17/90
c
      juflow=1
      do 42 j=2,nr
        if (lp*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue
c
c  njtj  *** end major modification  ***
c
c   determine effective charge and vzero for startup of
c   outward integration
c   ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
c   aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
c
      zeff = zero
      if (so(iorb) .lt. 0.1 .and. viod(lp,2) .lt. -0.1) zeff=znuc
      if (so(iorb) .gt. 0.1 .and. viou(lp,2) .lt. -0.1) zeff=znuc
      aa = -zeff/lp
      vzero = -2*zeff*aa
      if (zeff .eq. zero) then
        if (so(iorb) .lt. 0.1 ) then
          vzero=vzero+viod(lp,2)/r(2)
        else
          vzero=vzero+viou(lp,2)/r(2)
        endif
      endif
      if (so(iorb) .lt. 0.1) then
        vzero=vzero+vid(2)
      else
        vzero=vzero+viu(2)
      endif
      var0 = zero
      if (lo(iorb) .eq. 0) var0=-2*zeff
      if (lo(iorb) .eq. 1) var0=two
c
      emax = zero
      emin = -two*100000
      if (ev(iorb) .gt. emax) ev(iorb) = emax
 10   if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15   format(' iorb =',i3,' iter =',i3,' ev =',e18.10,' nodes =',i2)
      if (itmax .eq. 0) then
        iflag =1
        return
      endif
      if (ev(iorb) .gt. zero) then
        write(6,1000)iorb
        call ext(620+iorb)
      endif
 1000 format(//,' error in difnrl - ev(',i2,
     1 ') greater then v(infinty)')
c
c   find practical infinity ninf and classical turning
c   point nctp for orbital
c
      icount=0
 20   icount=icount+1
      do 22 j=nr,2,-1
        temp = v(j) -ev(iorb)
        if (temp .lt. zero) temp = zero
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
 22   continue
 23   ninf=j
      nctp = ninf - 5
      do 25 j=2,ninf-5
        if (v(j) .lt. ev(iorb)) nctp = j
 25   continue
      if (ev(iorb) .ge. etol*10) nctp=ninf-5
      if (ev(iorb) .ge. etol) ev(iorb)=zero
      if (evi(iorb) .ne. zero) then
        ev(iorb) = evi(iorb)
        do 26 j=1,nr
          if (r(j) .lt. five) nctp=j
 26     continue
      endif
      if (nctp .le. 6) then
        ev(iorb) = pnine*ev(iorb)
        if (icount .gt. 100) then
          write(6,1010)iorb
          call ext(650+iorb)
        endif
        goto 20
      endif
 1010 format(//,'error in difnrl - cannot find the classical '
     1 ,/' turning point for orbital ',i2)
c
c   outward integration from 1 to nctp
c   startup
c
      bb = (vzero-ev(iorb))/(4*lp+2)
      do 35 j=2,5
        ar(j) = rlp(j) * (1+(aa+bb*r(j))*r(j))
        br(j) = rabrlo(j) * (lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
 35   continue
c
c  njtj  ***  start major modification  ***
c    Predictor-corrector array added.
c
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-ev(iorb))*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-ev(iorb))*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-ev(iorb))*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-ev(iorb))*ar(5)
c
c   intergration loop
c
      nodes = 0
      do 40 j=6,nctp
c
c   predictor (Adams-Bashforth)
c
        j1=j-1
        j2=j-2
        j3=j-3
        j4=j-4
        j5=j-5
        vev=v(j)-ev(iorb)
        arp = ar(j1) + abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5)
        brp = br(j1) + abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5)
        fb1 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) + amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4)
        brc = br(j1) + amc0*fb1+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4)
        fb0 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc + amc0*(brc-brp)
        br(j) = brc + amc0*(fb0-fb1)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
c
c   count nodes - if no underflow
c
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.zero)nodes=nodes+1
 40   continue
c
c  njtj  ***  end major modification  ***
c
      arctp = ar(nctp)
      brctp = br(nctp)
c
c   end outward integration
c
c   if number of nodes correct, start inward integration
c   else modify energy stepwise and try again
c
      if (evi(iorb) .ne. zero) goto 111
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
        if (nodes .lt. no(iorb)-lo(iorb)-1) then
c
c  too few nodes; increase ev
c
          if (ev(iorb) .gt. emin) emin = ev(iorb)
          ev(iorb) = ev(iorb) - ev(iorb)/10
        else
c
c  too many nodes; decrease ev
c
          if (ev(iorb) .lt. emax) emax = ev(iorb)
          ev(iorb) = ev(iorb) + ev(iorb)/10
        endif
        itmax = itmax-1
        goto 10
      endif
c
c   inward integration from ninf to nctp
c   startup
c
      do 71 j=ninf,ninf-4,-1
        alf = v(j) - ev(iorb)
        if (alf .lt. zero) alf = zero
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = -rab(j)*alf*ar(j)
 71   continue
c
c  njtj  ***  start major modification  ***
c    Array for predictor-corrector added.
c
      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*
     1 (v(ninf)-ev(iorb))*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*
     1       (v(ninf1)-ev(iorb))*ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*
     1       (v(ninf2)-ev(iorb))*ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*
     1       (v(ninf3)-ev(iorb))*ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*
     1       (v(ninf4)-ev(iorb))*ar(ninf4)
c
c   integration loop
c
      istop = ninf - nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
c
c   predictor (Adams-Bashforth)
c
        j1 = j + 1
        j2 = j + 2
        j3 = j + 3
        j4 = j + 4
        j5 = j + 5
        vev = v(j)-ev(iorb)
        arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5))
        brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5))
        fb0 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4))
        brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4))
c
        fb1 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc - amc0*(brc-brp)
        br(j) = brc - amc0*(fb1-fb0)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
 80   continue
c
c   end inward integration
c
c  njtj  *** end major modification  ***
c
c   rescale ar and br outside nctp to match ar(nctp) from
c   outward integration
c
  222 factor = arctp/ar(nctp)
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
c
c   find normalizing factor
c
      factor = zero
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*ar(j)*ar(j)*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
c
c   modify eigenvalue ev
c
      dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev(iorb)) dev=dsign(ev(iorb),dev)/5
      itmax = itmax-1
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) / 2
      if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) / 2
      if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c   normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
c
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j) / rab(j)
 110  continue
 111  continue
      if (evi(iorb) .ne. zero) then
        factor = zero
        ll = 4
        do 112 j=2,nctp
          factor = factor + ll*ar(j)*ar(j)*rab(j)
          ll = 6 - ll
 112    continue
        factor = factor / 3
        factor = 1 / sqrt(factor)
        do 113 j=1,nctp
          ar(j) = factor*ar(j)
          br(j) = factor*br(j) / rab(j)
 113    continue
      endif
      return
      end
C
C
C
      subroutine difrel(iter,iorb,v,ar,br,lmax,nr,a,b,r,rab,
     1 norb,no,lo,so,znuc,viod,viou,vid,viu,ev,rabkar,
     2 rabai,fa,fb,evi)
c
c  difrel integrates the relativistic Dirac equation
c  it finds the eigenvalue ev, the major and minor component
c  of the wavefunction, ar and br.  It uses an intial guess
c  for the eigenvalues from dsolv1
c
c  njtj  ***  modifications  ***
c    This routine has major modifications.
c    1)The data needed inside the loops has been calculated
c    outside the main loop(increases speed for non-opt
c    compiliers, i.e. dumb compiliers).
c    2)The predict/correct values are placed in an array.
c    Output is unchanged
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch the 3 double precision parameter
c  ###      to single precision parameter statements.
c  ###  Cray conversions
c  njtj
c
c  njtj
c  &&&  Machine dependent Parameter
c  &&&    The value of expzer is machine dependent.
c  &&&    The user must switch in the correct value for
c  &&&    the machine in use from the list, or find
c  &&&    it for their machine.
c  &&&  Machine dependent Parameter
c  njtj
c
      implicit double precision (a-h,o-z)
      parameter (zero=0.0,pnine=0.9,one=1.0,ai=2*137.0360411)
      parameter (etol=-1.E-7)
c
c  Tolernce
c
      parameter (tol = 1.D-10,five=5.0D0)
Cray      parameter (tol = 1.E-10,five=5.0)
c
c  Integration coefficients
c
      parameter(abc1=190.1/72,abc2=-138.7/36,abc3=10.9/3,
     1 abc4=-63.7/36,abc5=25.1/72,amc0=25.1/72,amc1=32.3/36,
     2 amc2=-1.1/3,amc3=5.3/36,amc4=-1.9/72)
c
c
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr),
     1 no(norb),lo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),rabkar(nr),rabai(nr),
     3 fa(nr),fb(nr),evi(norb)
c
      dimension rs(5)
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
cApollo      expzer = 3.7D2
cSun      expzer = 3.7D2
      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
      itmax = 100
      ai2 = ai * ai
      az = znuc/(2*ai)
      ka = lo(iorb)+1
      if (so(iorb) .lt. 0.1 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c  determine effective charge and vzero for startup of
c  outward integration
c  ar = r**s * (1  + a1 r + a2 r**2 + ... )
c  br = r**s * (b0 + b1 r + b2 r**2 + ... )
c  s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
c  an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
c        / (n ai (2 s + n))
c  bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))
c
      s = sqrt(ka*ka-az*az)
      if (ka .gt. 0) then
        b0 = -az/(s+ka)
      else
        b0 = (s-ka)/az
      endif
      if (so(iorb) .lt. 0.1) then
        vzero=vid(2)
      else
        vzero=viu(2)
      endif
c
c  njtj  ***  start major modification  ***
c    Loop data calculated only once.
c    Set ar() and br() to zero.
c
      do 1 j=1,nr
        ar(j) = zero
        br(j) = zero
 1    continue
      do 3 j=2,nr
        rabkar(j)=rab(j)*ka/r(j)
 3    continue
      do 4 j=2,nr
        rabai(j)=rab(j)/ai
 4    continue
      do 5 j=2,5
        rs(j)=r(j)**s
 5    continue
c
c  set the underflow trap, error from Berkeley version,
c  fixed by Troy Barbee, sqrt(expzer) should be expzer/2,
c  4/17/90.
c
      juflow=1
      do 42 j=2,nr
        if (s*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue
c  njtj *** end major modification  ***
c
      emax = zero
      emin = -one*100000
      if (ev(iorb) .gt. emax) ev(iorb) = emax
 10   if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15   format(' iorb =',i3,' iter =',i3,' ev =',e18.10,' nodes =',i2)
      if (itmax .eq. 0) return
      if (ev(iorb) .gt. zero) then
        write(6,1000)iorb
        call ext(620+iorb)
      endif
 1000 format(//,' error in difrel - ev(',i2,
     1 ') greater then v(infinty)')
c
c  Find practical infinity ninf and classical turning
c  point nctp for orbital.
c
      icount=0
 20   icount=icount+1
      do 22 j=nr,2,-1
        temp = v(j) - ev(iorb)
        if (temp .lt. zero) temp = zero
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
 22   continue
 23   ninf=j
      nctp = ninf - 5
      do 25 j=2,ninf-5
        if (v(j) .lt. ev(iorb)) nctp = j
 25   continue
      if (ev(iorb) .ge. etol*100) nctp=ninf-5
      if (ev(iorb) .ge. etol) ev(iorb)=zero
      if (evi(iorb) .ne. zero) then
        ev(iorb)=evi(iorb)
        do 26 j=2,nr
          if (r(j) .lt. five) nctp=j
 26     continue
      endif
      if (nctp .le. 6) then
        ev(iorb) = pnine*ev(iorb)
        if (icount .gt. 100) then
          write(6,1010)iorb
          call ext(650+iorb)
        endif
        goto 20
      endif
 1010 format(//,'error in difrel - cannot find classical',
     1 /,'turning point in orbital ',i2)
c
c  Outward integration from 1 to nctp, startup.
c
      a1 = (az*(vzero-ev(iorb))-(s+1+ka)*(vzero-ev(iorb)-ai2)*b0)
     1   / (ai*(2*s+1))
      b1 = ((vzero-ev(iorb))-2*znuc*a1) / (ai*(s+1+ka))
      a2 = (az*(vzero-ev(iorb))*a1-(s+2+ka)*(vzero-ev(iorb)-ai2)*b1)
     1   / (2*ai*(2*s+2))
      b2 = ((vzero-ev(iorb))*a1-2*znuc*a2) / (ai*(s+2+ka))
      do 35 j=2,5
        ar(j) = rs(j) * (1 +(a1+a2*r(j))*r(j))
        br(j) = rs(j) * (b0+(b1+b2*r(j))*r(j))
 35   continue
      fa(1) = zero
      fb(1) = zero
      fa(2) = rabkar(2)*ar(2)+(ev(iorb)-v(2)+ai2)*br(2)*rabai(2)
      fb(2) = -rabkar(2)*br(2)-(ev(iorb)-v(2))*ar(2)*rabai(2)
      fa(3) = rabkar(3)*ar(3)+(ev(iorb)-v(3)+ai2)*br(3)*rabai(3)
      fb(3) = -rabkar(3)*br(3)-(ev(iorb)-v(3))*ar(3)*rabai(3)
      fa(4) = rabkar(4)*ar(4)+(ev(iorb)-v(4)+ai2)*br(4)*rabai(4)
      fb(4) = -rabkar(4)*br(4)-(ev(iorb)-v(4))*ar(4)*rabai(4)
      fa(5) = rabkar(5)*ar(5)+(ev(iorb)-v(5)+ai2)*br(5)*rabai(5)
      fb(5) = -rabkar(5)*br(5)-(ev(iorb)-v(5))*ar(5)*rabai(5)
c
c  Intergration loop.
c
      nodes = 0
      do 40 j=6,nctp
c
c  Predictor (Adams-Bashforth).
c
        evvai2=ev(iorb)-v(j)+ai2
        evv=ev(iorb)-v(j)
        arp = ar(j-1) + abc1*fa(j-1)+abc2*fa(j-2)+abc3*fa(j-3)
     1   +abc4*fa(j-4)+abc5*fa(j-5)
        brp = br(j-1) + abc1*fb(j-1)+abc2*fb(j-2)+abc3*fb(j-3)
     1   +abc4*fb(j-4)+abc5*fb(j-5)
        fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
        fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
        arc = ar(j-1) + amc0*fa(j)+amc1*fa(j-1)+amc2*fa(j-2)
     1   +amc3*fa(j-3)+amc4*fa(j-4)
        brc = br(j-1) + amc0*fb(j)+amc1*fb(j-1)+amc2*fb(j-2)
     1   +amc3*fb(j-3)+amc4*fb(j-4)
        faj = rabkar(j)*arc+evvai2*brc*rabai(j)
        fbj = -rabkar(j)*brc-evv*arc*rabai(j)
c
c  Error reduction step.
c
        ar(j) = arc + amc0*(faj-fa(j))
        br(j) = brc + amc0*(fbj-fb(j))
        fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
        fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)
c
c  Count nodes - if no underflow.
c
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.zero)nodes=nodes+1
 40   continue
       arout = ar(nctp)
       arpout = fa(nctp)
c
c  End outward integration.
c  If number of nodes correct, start inward integration
c  else modify energy stepwise and try again.
c
      if (evi(iorb) .ne. zero) goto 111
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
c
c  too many nodes decrease ev
c
        if (nodes .gt. no(iorb)-lo(iorb)-1) then
          if (ev(iorb) .lt. emax) emax = ev(iorb)
          ev(iorb) = ev(iorb) + ev(iorb)/10
c
c  too few nodes increase ev
c
        else
          if (ev(iorb) .gt. emin) emin = ev(iorb)
          ev(iorb) = ev(iorb) - ev(iorb)/10
        endif
        itmax = itmax-1
        goto 10
      endif
c
c  Inward integration from ninf to nctp startup.
c
      do 70 j=ninf,ninf-4,-1
        alf = v(j) - ev(iorb)
        if (alf .lt. zero) alf = zero
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = ai*(alf+ka/r(j))*ar(j)/(v(j)-ev(iorb)-ai2)
 70   continue
      fa(ninf) = rabkar(ninf)*ar(ninf)+
     1    (ev(iorb)-v(ninf)+ai2)*br(ninf)*rabai(ninf)
      fb(ninf) = -rabkar(ninf)*br(ninf)
     1    -(ev(iorb)-v(ninf))*ar(ninf)*rabai(ninf)
      fa(ninf-1) = rabkar(ninf-1)*ar(ninf-1)+
     1    (ev(iorb)-v(ninf-1)+ai2)*br(ninf-1)*rabai(ninf-1)
      fb(ninf-1) = -rabkar(ninf-1)*br(ninf-1)
     1    -(ev(iorb)-v(ninf-1))*ar(ninf-1)*rabai(ninf-1)
      fa(ninf-2) = rabkar(ninf-2)*ar(ninf-2)
     1    +(ev(iorb)-v(ninf-2)+ai2)*br(ninf-2)*rabai(ninf-2)
      fb(ninf-2) = -rabkar(ninf-2)*br(ninf-2)
     1    -(ev(iorb)-v(ninf-2))*ar(ninf-2)*rabai(ninf-2)
      fa(ninf-3) = rabkar(ninf-3)*ar(ninf-3)
     1    +(ev(iorb)-v(ninf-3)+ai2)*br(ninf-3)*rabai(ninf-3)
      fb(ninf-3) = -rabkar(ninf-3)*br(ninf-3)
     1    -(ev(iorb)-v(ninf-3))*ar(ninf-3)*rabai(ninf-3)
      fa(ninf-4) = rabkar(ninf-4)*ar(ninf-4)
     1    +(ev(iorb)-v(ninf-4)+ai2)*br(ninf-4)*rabai(ninf-4)
      fb(ninf-4) = -rabkar(ninf-4)*br(ninf-4)
     1    -(ev(iorb)-v(ninf-4))*ar(ninf-4)*rabai(ninf-4)
c
c  Integration loop.
c
      istop = ninf-nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
c
c  Predictor (Adams-Bashforth).
c
        evvai2=ev(iorb)-v(j)+ai2
        evv=ev(iorb)-v(j)
        arp = ar(j+1)-(abc1*fa(j+1)+abc2*fa(j+2)+abc3*fa(j+3)
     1   +abc4*fa(j+4)+abc5*fa(j+5))
        brp = br(j+1)-(abc1*fb(j+1)+abc2*fb(j+2)+abc3*fb(j+3)
     1   +abc4*fb(j+4)+abc5*fb(j+5))
        fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
        fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
        arc = ar(j+1)-(amc0*fa(j)+amc1*fa(j+1)+amc2*fa(j+2)
     1   +amc3*fa(j+3)+amc4*fa(j+4))
        brc = br(j+1)-(amc0*fb(j)+amc1*fb(j+1)+amc2*fb(j+2)
     1   +amc3*fb(j+3)+amc4*fb(j+4))
        faj = rabkar(j)*arc+evvai2*brc*rabai(j)
        fbj = -rabkar(j)*brc-evv*arc*rabai(j)
c
c  Error reduction step.
c
        ar(j) = arc + amc0*(faj-fa(j))
        br(j) = brc + amc0*(fbj-fb(j))
        fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
        fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)
 80   continue
 222  arin = ar(nctp)
      arpin = fa(nctp)
c
c  End inward integration
c  Rescale ar and br outside nctp to match ar(nctp) from
c  outward integration.
c
      factor = arout/arin
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
      arpin = factor * arpin
c
c  Find the normalizing factor.
c
      factor = zero
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
c
c  Modify the eigenvalue ev.
c
      dev = arout * (arpout-arpin) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev(iorb)) dev=dsign(ev(iorb),dev)/5
      itmax = itmax-1
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) .gt. emax) then
        ev(iorb) = (evold + emax) / 2
      elseif (ev(iorb) .lt. emin) then
        ev(iorb) = (evold + emin) / 2
      endif
      if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c  Normalize the wavefunction.
c
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j)
 110  continue
 111  continue
      if (evi(iorb) .ne. zero) then
        factor = zero
        ll = 4
        do 112 j=2,nctp
          factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
          ll = 6 - ll
 112    continue
        factor = factor / 3
        factor = 1 / sqrt(factor)
        do 113 j=1,nctp
          ar(j) = factor*ar(j)
          br(j) = factor*br(j)
 113    continue
      endif
      return
      end
C
C
C
      SUBROUTINE DMIXP(A,B,BETA,ICY,ID,NMSH,
     1                 C,D,VN1,VN12,VN2,VN22)
C*    ADAPTED FROM K.C.PANDEY
C*    USING ANDERSON'S EXTRAPOLATION SCHEME
C*    EQS 4.1-4.9,4.15-4.18 OF
C*    D.G.ANDERSON J.ASSOC.COMPUTING MACHINERY,12,547(1965)
C*    COMPUTES A NEW VECTOR IN A ITERATIVE SCHEME
C*    INPUT A=NEWPOT B=OLDPOT
C*    OUTPUT A=A-B B=NEWPOT
C*    BETA=MIXING,IN=ITER. NUMBER
C*    ID=1,2 OR 3 DIFF CONV METH.
C*    ICY CYCLE NUMBER ,ICY=1 ON FIRST/ZEROTH CALL
C*    C,D WORK ARRAYS OF SIZE NMSH
C*    VN1,VN12,VN2,VN22 STORAGE ARRAYS OF SIZE NMSH
C
C  njtj
C  ###  Cray conversions
C  ###    1)Comment out implicit double precision.
C  ###    2)Switch double precision parameter
C  ###      to single precision parameter statement.
C  ###  Cray conversions
C  njtj
C
C
      implicit double precision (a-h,o-z)
      PARAMETER (UZE=0.0D0,UM=1.0D0,DETOL=1.D-9)
Cray      PARAMETER (UZE=0.0,UM=1.0,DETOL=1.E-9)
C
      DIMENSION A(NMSH),B(NMSH),C(NMSH),D(NMSH)
      DIMENSION VN1(NMSH),VN12(NMSH),VN2(NMSH),VN22(NMSH)
      IN=ICY-1
      IF(IN.EQ.0) THEN
        CALL TRNSVV(B,A,UM,NMSH)
        RETURN
      ENDIF
      CALL TRNSVV(A,B,-UM,NMSH)
      CALL DOTTVV(A,A,R2,NMSH)
      IF(ID.EQ.1) THEN
        CALL TRNSVV(B,A,BETA,NMSH)
        RETURN
      ENDIF
      IF(IN.EQ.1) THEN
        DO 100 I=1,NMSH
          VN1(I)=A(I)
 100    CONTINUE
        DO 105 I=1,NMSH
          VN2(I)=B(I)
 105    CONTINUE
        CALL TRNSVV(B,A,BETA,NMSH)
        RETURN
      ENDIF
      DO 110 I=1,NMSH
        C(I)=VN1(I)
 110  CONTINUE
      IF(ID.EQ.3.AND.IN.GT.2) THEN
        DO 115 I=1,NMSH
          D(I)=VN12(I)
 115    CONTINUE
      ENDIF
      DO 120 I=1,NMSH
        VN1(I)=A(I)
 120  CONTINUE
      IF(ID.GT.2.AND.IN.GT.1) THEN
        DO 125 I=1,NMSH
          VN12(I)=C(I)
 125    CONTINUE
      ENDIF
      CALL TRNSVV(C,A,-UM,NMSH)
      CALL DOTTVV(C,C,D11,NMSH)
      CALL DOTTVV(A,C,RD1M,NMSH)
      IF(IN.LE.2.OR.ID.LE.2) THEN
        T1=-RD1M/D11
        X=UM-T1
        BT1=BETA*T1
        DO 5 I=1,NMSH
          A(I)=BETA*A(I)
 5      CONTINUE
        CALL TRNSVV(A,C,BT1,NMSH)
        DO 130 I=1,NMSH
          D(I)=VN2(I)
 130    CONTINUE
        CALL TRNSVV(A,D,T1,NMSH)
        DO 135 I=1,NMSH
          VN2(I)=B(I)
 135    CONTINUE
        IF(ID.GT.2.AND.IN.EQ.2) THEN
          DO 140 I=1,NMSH
            VN22(I)=D(I)
 140      CONTINUE
        ENDIF
        DO 10 I=1,NMSH
          B(I)=X*B(I)+A(I)
 10     CONTINUE
        RETURN
      ENDIF
      CALL TRNSVV(D,A,-UM,NMSH)
      CALL DOTTVV(D,D,D22,NMSH)
      CALL DOTTVV(C,D,D12,NMSH)
      CALL DOTTVV(A,D,RD2M,NMSH)
      A2=D11*D22
      DET=A2-D12*D12
      DETT=DET/A2
      IF(ABS(DETT).GE.DETOL) THEN
        T1=(-RD1M*D22+RD2M*D12)/DET
        T2=( RD1M*D12-RD2M*D11)/DET
      ELSE
        T1=-RD1M/D11
        T2=UZE
      ENDIF
      X=UM-T1-T2
      BT1=BETA*T1
      BT2=BETA*T2
      DO 15 I=1,NMSH
        A(I)=BETA*A(I)
 15   CONTINUE
      CALL TRNSVV(A,C,BT1,NMSH)
      CALL TRNSVV(A,D,BT2,NMSH)
      CALL TRNSVV(A,VN2,T1,NMSH)
      CALL TRNSVV(A,VN22,T2,NMSH)
      DO 145 I=1,NMSH
        VN22(I)=VN2(I)
 145  CONTINUE
      DO 155 I=1,NMSH
        VN2(I)=B(I)
 155  CONTINUE
      DO 20 I=1,NMSH
        B(I)=X*B(I)+A(I)
 20   CONTINUE
      RETURN
      END
C
C
C
      subroutine dottvv(a,b,c,n)
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch 1 function line from double
c  ###    precision to a single precision statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      dimension a(n),b(n)
c
       c=0.0
c
      do 10 i=1,n
        c=c+a(i)*b(i)
 10   continue
      return
      end
C
C
C
       subroutine dsolv1(lmax,nr,a,b,r,rab,norb,ncore,
     1  no,lo,so,zo,cdd,cdu,viod,viou,vid,viu,
     3  ev,dk,d,sd,sd2,rv1,rv2,rv3,rv4,rv5,z)
c
c   dsolv1 finds the (non)-relativistic wave function
c   using finite differences and matrix diagonalization.
c   An initial guess for the eigenvalues need not be supplied.
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,one=1.D0,pone=0.1D0,opf=1.5D0)
Cray      parameter (zero=0.0,one=1.0,pone=0.1,opf=1.5)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),
     1 zo(norb),cdd(nr),cdu(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),dk(nr),d(nr),sd(nr),sd2(nr),
     2 z(6*nr),rv1(nr),rv2(nr),rv3(nr),rv4(nr),rv5(nr)
c
      dimension nmax(2,5),e(10),ind(10)
c
c   Initialize the charge density arrays.
c
       do 10 i=1,nr
         cdd(i) = zero
         cdu(i) = zero
 10    continue
c
c   Find the max n given l and s.
c   Zero spin is treated as down.
c
      do 20 i=1,2
        do 20 j=1,lmax
          nmax(i,j) = 0
          do 20 k=1,norb
            if (no(k) .le. 0) goto 20
            if (lo(k) .ne. j-1) goto 20
            if ((so(k)-pone)*(i-opf) .lt. zero) goto 20
            nmax(i,j)=no(k)
 20   continue
c
c   Set up hamiltonian matrix for kinetic energy.
c   Only the diagonal depends on the potential.
c
      c2 = -one/b**2
      c1 = -2*one*c2 + one/4
      dk(1)  = c1 / (r(2)+a)**2
      sd(1)  = zero
      sd2(1) = zero
      do 30 i=3,nr
        dk(i-1)  = c1 / (r(i)+a)**2
        sd(i-1)  = c2 / ((r(i)+a)*(r(i-1)+a))
        sd2(i-1) = sd(i-1)**2
 30   continue
c
c   Start loop over spin down=1 and spin up=2.
c
      nrm = nr - 1
      do 80 i=1,2
c
c   Start loop over s p d... states.
c
        do 80 j=1,lmax
          if (nmax(i,j) .eq. 0) goto 80
          llp = j*(j-1)
          do 40 k=2,nr
            if (i .eq. 1) then
              d(k-1)=dk(k-1)+(viod(j,k)+llp/r(k))/r(k)+vid(k)
            else
              d(k-1)=dk(k-1)+(viou(j,k)+llp/r(k))/r(k)+viu(k)
            endif
 40       continue
c
c   Diagonalize the matrix.
c
          eps = -one
          call tridib(nrm,eps,d,sd,sd2,bl,bu,1,
     1     nmax(i,j),e,ind,ierr,rv4,rv5)
          if (ierr .ne. 0) write(6,50) ierr
 50   format(/,' error in tridib ****** ierr =',i3,/)
          call tinvit(nrm,nrm,d,sd,sd2,nmax(i,j),e,ind,z,ierr,
     1     rv1,rv2,rv3,rv4,rv5)
          if (ierr .ne. 0) write(6,55) ierr
 55   format(/,' error in tinvit ****** ierr =',i3,/)
c
c   Save the energy levels and add to charge density.
c
          ki = 1
          kn = 0
          do 70 k=1,norb
            if (no(k) .le. 0) goto 70
            if (lo(k) .ne. j-1) goto 70
            if ((so(k)-pone)*(i-opf) .lt. zero) goto 70
            ev(k) = e(ki)
            do 60 l=2,nr
            denr = zo(k) * z(kn+l-1)**2 / rab(l)
            if (i .eq. 1) then
              cdd(l) = cdd(l) + denr
            else
              cdu(l) = cdu(l) + denr
            endif
 60       continue
          ki = ki + 1
          kn = kn + nrm
 70     continue
 80   continue
c
c   End loop over s p and d states.
c
      return
      end
C
C
C
      subroutine dsolv2(iter,iconv,ispp,ifcore,lmax,nr,a,b,r,
     1 rab,norb,ncore,no,lo,so,zo,znuc,cdd,cdu,cdc,viod,
     2 viou,vid,viu,ev,ek,ep,wk1,wk2,wk3,wk4,v,ar,br,evi)
c
c  dsolv2 finds the (non) relativistic wave function using
c  difnrl to intgrate the Scroedinger equation or
c  difrel to intgrate the Dirac equation.
c  The energy level from the previous iteration is used
c  as initial guess, and it must therefore be reasonable
c  accurate.
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      character*1 ispp
c
      parameter (zero=0.D0,smev=1.D-4)
Cray      parameter (zero=0.0,smev=1.E-4)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),ek(norb),ep(norb),evi(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),v(nr),ar(nr),br(nr)
c
c  Initialize arrays for charge density.
c
      do 5 i=1,nr
        cdd(i) = zero
 5    continue
      do 10 i=1,nr
        cdu(i) = zero
 10   continue
      if (ifcore .eq. 0) then
        do 15 i=1,nr
          cdc(i)= zero
 15     continue
      endif
c
c  Start the loop over orbitals.
c  Note that spin zero is treated as down.
c
      do 50 i=1,norb
        if (no(i) .le. 0) goto 50
        if (zo(i) .eq. 0.0 .and. iconv .eq. 0) goto 50
        if (ev(i) .ge. 0.0) ev(i)=-smev
c
c  Set up the potential, set the wave functionc array to zero-ar.
c
        lp  = lo(i)+1
        llp = lo(i)*lp
        do 17 j=1,nr
          ar(j)=zero
 17     continue
        if (so(i) .lt. 0.1) then
          do 18 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 18       continue
        else
          do 19 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 19       continue
        endif
        if (ispp .ne. 'r') then
          do 20 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 20       continue
        endif
c
c  Call the integration routine.
c
        if (ispp .ne. 'r') then
          call difnrl(iter,i,v,ar,br,lmax,nr,a,b,r,
     1     rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2     ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(iter,i,v,ar,br,lmax,nr,a,b,r,
     1     rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2     ev,wk1,wk2,wk3,wk4,evi)
        endif
c
c  Add to the charge density.
c
       if (ispp .eq. 'r') then
         if (so(i) .lt. 0.1) then
           do 30 j=1,nr
             denr = zo(i) *(br(j) * br(j) + ar(j) * ar(j))
             cdd(j) = cdd(j) + denr
 30        continue
         else
           do 31 j=1,nr
             denr = zo(i) *(br(j) * br(j) + ar(j) * ar(j))
             cdu(j) = cdu(j) + denr
 31        continue
         endif
       else
         if (so(i) .lt. 0.1) then
           do 32 j=1,nr
             denr = zo(i) * ar(j) * ar(j)
             cdd(j) = cdd(j) + denr
 32        continue
         else
           do 33 j=1,nr
             denr = zo(i) * ar(j) * ar(j)
             cdu(j) = cdu(j) + denr
 33        continue
         endif
       endif
       if (ifcore .eq. 0 .and. i .le. ncore) then
         do 34 j=1,nr
           denr = zo(i) * ar(j) * ar(j)
           cdc(j)=cdc(j)+denr
 34      continue
       endif
c
c  Compute various quantitities if last iteration.
c
        if (iconv .eq. 1) call orban(ispp,i,ar,br,
     1   lmax,nr,a,b,r,rab,norb,no,lo,zo,so,viod,viou,
     2   vid,viu,ev,ek,ep)
 50   continue
c
c  End loop over orbitals.
c
      return
      end
C
C
C
      subroutine etotal(itype,zsh,nameat,norb,
     1 no,lo,so,zo,etot,ev,ek,ep)
c
c  etotal computes the total energy from the
c  electron charge density.
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
      parameter (zero=0.D0)
Cray      parameter (zero=0.0)
c
c
      character*1 il(5)
      character*2 nameat
c
      dimension no(norb),lo(norb),so(norb),zo(norb),
     1 etot(10),ev(norb),ek(norb),ep(norb)
c
c      etot(i)    i=1,10 contains various contributions to the total
c                 energy.
c                 (1)   sum of eigenvalues ev
c                 (2)   sum of orbital kinetic energies ek
c                 (3)   el-ion interaction from sum of orbital
c                       potential energies ep
c                 (4)   electrostatic el-el interaction  (from velect)
c                 (5)   vxc (exchange-correlation) correction to sum
c                       of eigenvalues                   (from velect)
c                 (6)   3 * vc - 4 * ec
c                       correction term for virial theorem
c                       when correlation is included     (from velect)
c                 (7)   exchange and correlation energy  (from velect)
c                 (8)   kinetic energy from eigenvalues  (1,3,4,5)
c                 (9)   potential energy
c                 (10)  total energy
c
c
c      sum up eigenvalues ev, kinetic energies ek, and
c      el-ion interaction ep
c
      etot(1) = zero
      etot(2) = zero
      etot(3) = zero
      do 10 i=1,norb
        etot(1) = etot(1) + zo(i)*ev(i)
        etot(2) = etot(2) + zo(i)*ek(i)
        etot(3) = etot(3) + zo(i)*ep(i)
 10   continue
c
c   kinetic energy
c
      etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)
c
c   potential energy
c
      etot(9) = etot(3) + etot(4) + etot(7)
c
c      total energy
c
      etot(10) = etot(1) - etot(4) - etot(5) + etot(7)
c
c   printout
c
      il(1) = 's'
      il(2) = 'p'
      il(3) = 'd'
      il(4) = 'f'
      il(5) = 'g'
      write(6,20) nameat
 20   format(//,1x,a2,' output data for orbitals',/,1x,28('-'),//,
     1 ' nl    s      occ',9x,'eigenvalue',4x,'kinetic energy',
     2 6x,'pot energy',/)
      do 40 i=1,norb
        write(6,30) no(i),il(lo(i)+1),so(i),zo(i),ev(i),ek(i),ep(i)
 30   format(1x,i1,a1,f6.1,f10.4,3f17.8)
 40   continue
      write(6,50) (etot(i),i=1,10)
 50   format(//,' total energies',/,1x,14('-'),/,
     1 /,' sum of eigenvalues        =',f18.8,
     2 /,' kinetic energy from ek    =',f18.8,
     3 /,' el-ion interaction energy =',f18.8,
     4 /,' el-el  interaction energy =',f18.8,
     5 /,' vxc    correction         =',f18.8,
     6 /,' virial correction         =',f18.8,
     7 /,' exchange + corr energy    =',f18.8,
     8 /,' kinetic energy from ev    =',f18.8,
     9 /,' potential energy          =',f18.8,/,1x,45('-'),
     X /,' total energy              =',f18.8)
       if (itype .ge. 4 .or. abs(zsh) .gt. 0.00001) return
c
c   virial theorem
c
       vsum = 2*etot(8) + etot(9) + etot(6)
       write(6,60) 2*etot(8),etot(9),etot(6),vsum
 60    format(//,' virial theorem(nonrelativistic)',/,1x,14('-'),/,
     1 /,' kinetic energy  *  2      =',f18.8,
     2 /,' potential energy          =',f18.8,
     3 /,' virial correction         =',f18.8,/,1x,45('-'),
     4 /,' virial sum                =',f18.8)
       return
       end
C
C
C
      subroutine ext(i)
c
c  Stops program in case of errors or completion.
c
c  i is a stop parameter
c   000-099 main (0 is normal exit)
c   100-199 input
c   200-299 charge
c   300-399 vionic
c   400-499 velect
c   500-599 dsolv1
c   600-699 dsolv2 (including difnrl and difrel)
c   700-799 etotal
c   800-899 pseudo, pseudk, pseudt and pseudv
c
      if (i .ne. 0) write(6,10) i
 10   format('stop parameter =',i3)
      close (unit=1)
      close (unit=3)
      close (unit=5)
      close (unit=6)
      call exit
      end
C
C
C
      subroutine gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,gamma,
     2 alpha,alpha1,alpha2,alpha3,alpha4,v0pp,ar)
c
c *********************************************************
c *                                                       *
c *  njtj                                                 *
c *   Retuns the values of delta, alpha, alpha1, alpha2,  *
c *   alpha3, and alpha4 given a fixed value of gamma.    *
c *   Returns V"(0) for the braketing and bisection       *
c *   routines.  Subroutine used in pseudtk routine.      *
c *  njtj                                                 *
c *                                                       *
c *********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      dimension r(jrc),rab(jrc),aj(5,5),bj(5),ar(jrc)
c
      parameter (zero=0.D0,pfive=0.5D0,one=1.D0,errmin=1.D-12)
Cray      parameter (zero=0.0,pfive=0.5,one=1.0,errmin=1.E-12)
c
      rc9  = rc8*rc1
      rc10 = rc8*rc2
      rc11 = rc8*rc3
      rc12 = rc8*rc4
      delta=zero
      bj(1)=log(arc/rc1**lp)-gamma*rc2
      bj1=bj(1)
      bj(2)=brc-lp/rc1-2*gamma*rc1
      bj2a=bj(2)+2*gamma*rc1
      bj2=bj(2)
      bj(3)=vrc-ev-2*lp/rc1*bj2a-bj2a**2-2*gamma
      bj3=bj(3)
      bj3a=bj(3)+2*gamma
      bj(4)=vap+2*lp/rc2*bj2a-2*lp/rc1*bj3a-2*bj2a*bj3a
      bj4=bj(4)
      bj(5)=vapp-4*lp/rc3*bj2a+4*lp/rc2*bj3a-2*lp/rc1*bj4-2*bj3a**2
     1 -2*bj2a*bj4
      bj5=bj(5)
      aj(1,1)=rc4
      aj(1,2)=rc6
      aj(1,3)=rc8
      aj(1,4)=rc10
      aj(1,5)=rc12
      aj(2,1)=4*rc3
      aj(2,2)=6*rc5
      aj(2,3)=8*rc7
      aj(2,4)=10*rc9
      aj(2,5)=12*rc11
      aj(3,1)=12*rc2
      aj(3,2)=30*rc4
      aj(3,3)=56*rc6
      aj(3,4)=90*rc8
      aj(3,5)=132*rc10
      aj(4,1)=24*rc1
      aj(4,2)=120*rc3
      aj(4,3)=336*rc5
      aj(4,4)=720*rc7
      aj(4,5)=1320*rc9
      aj(5,1)=24*one
      aj(5,2)=360*rc2
      aj(5,3)=1680*rc4
      aj(5,4)=5040*rc6
      aj(5,5)=11880*rc8
      call gaussj(aj,5,5,bj,1,1)
      alpha=bj(1)
      alpha1=bj(2)
      alpha2=bj(3)
      alpha3=bj(4)
      alpha4=bj(5)
c
c   start iteration loop to find delta(with gamma fixed)
c
      do 550 j=1,200
c
c   generate pseudo wavefunction-note missing factor exp(delta)
c
        do 560 k=1,jrc
          rp=r(k)
          r2=rp*rp
          polyr = r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+
     1     alpha1)*r2+ alpha)*r2+gamma)
          ar(k) = rp**lp * exp(polyr)
 560    continue
c
c   integrate pseudo charge density from r = 0 to rc
c
        ll = 2
        cdps = - ar(jrc) * ar(jrc) * rab(jrc)
        if (jrc .ne. 2*(jrc/2)) then
          do 120 k=jrc,1,-1
            cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
            ll = 6 - ll
 120      continue
        else
          do 121 k=jrc,4,-1
            cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
            ll = 6 - ll
 121      continue
          cdps = cdps - ar(4) * ar(4) * rab(4)
          cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1     3 * ar(2) *ar(2) * rab(2) +
     2     3 * ar(3) *ar(3) * rab(3) +
     3     ar(4) * ar(4) * rab(4))/8
        endif
        cdps = cdps/3
c
c   Calculate new delta(with gamma fixed), uses false position
c
        fdnew = log(cdrc/cdps) - 2*delta
        if (abs(fdnew) .lt. errmin) then
          v0pp=8*((2*one*(lp-one)+5*one)*alpha+gamma**2)
          return
        endif
        if (j .eq. 1) then
          ddelta=-pfive
        else
          ddelta = - fdnew * ddelta / (fdnew-fdold)
        endif
        delta = delta + ddelta
        bj(1)=bj1-delta
        bj(2)=bj2
        bj(3)=bj3
        bj(4)=bj4
        bj(5)=bj5
        aj(1,1)=rc4
        aj(1,2)=rc6
        aj(1,3)=rc8
        aj(1,4)=rc10
        aj(1,5)=rc12
        aj(2,1)=4*rc3
        aj(2,2)=6*rc5
        aj(2,3)=8*rc7
        aj(2,4)=10*rc9
        aj(2,5)=12*rc11
        aj(3,1)=12*rc2
        aj(3,2)=30*rc4
        aj(3,3)=56*rc6
        aj(3,4)=90*rc8
        aj(3,5)=132*rc10
        aj(4,1)=24*rc1
        aj(4,2)=120*rc3
        aj(4,3)=336*rc5
        aj(4,4)=720*rc7
        aj(4,5)=1320*rc9
        aj(5,1)=24*one
        aj(5,2)=360*rc2
        aj(5,3)=1680*rc4
        aj(5,4)=5040*rc6
        aj(5,5)=11880*rc8
        call gaussj(aj,5,5,bj,1,1)
        alpha=bj(1)
        alpha1=bj(2)
        alpha2=bj(3)
        alpha3=bj(4)
        alpha4=bj(5)
        fdold = fdnew
 550  continue
      write(6,1000)
 1000 format(//, 'error in gamfind - delta not found')
      call ext(860+lp)
      end
C
C
C
      subroutine gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,gamma,
     2 alpha,alpha1,alpha2,alpha3,alpha4,v0pp,ar)
c
c *********************************************************
c *                                                       *
c *  njtj                                                 *
c *   Retuns the values of delta, alpha, alpha1, alpha2,  *
c *   alpha3, and alpha4 given a fixed value of gamma.    *
c *   Returns V"(0) for the braketing and bisection       *
c *   routines.  Subroutine used in pseudtk routine.      *
c *  njtj                                                 *
c *                                                       *
c *********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      dimension r(jrc),rab(jrc),aj(5,5),bj(5),ar(jrc)
c
      parameter (zero=0.D0,pfive=0.5D0,one=1.D0,errmin=1.D-12)
Cray      parameter (zero=0.0,pfive=0.5,one=1.0,errmin=1.E-12)
c
      delta=zero
      bj(1)=log(arc/rc1**lp)-gamma*rc2
      bj(2)=brc-lp/rc1-2*gamma*rc1
      bj(3)=vrc-ev+(lp/rc1)**2-brc**2-2*gamma
      vt=vrc-ev+lp*(lp-1)/rc2
      bj(4)=vap-2*(vt*brc+lp**2/rc3-brc**3)
      bj(5)=vapp-2*((vap-2*lp*(lp-1)/rc3)*brc+(vt-3*brc**2)*
     1 (vt-brc**2)-3*lp**2/rc4)
      aj(1,1)=rc4
      aj(1,2)=rc5
      aj(1,3)=rc6
      aj(1,4)=rc7
      aj(1,5)=rc8
      aj(2,1)=4*rc3
      aj(2,2)=5*rc4
      aj(2,3)=6*rc5
      aj(2,4)=7*rc6
      aj(2,5)=8*rc7
      aj(3,1)=12*rc2
      aj(3,2)=20*rc3
      aj(3,3)=30*rc4
      aj(3,4)=42*rc5
      aj(3,5)=56*rc6
      aj(4,1)=24*rc1
      aj(4,2)=60*rc2
      aj(4,3)=120*rc3
      aj(4,4)=210*rc4
      aj(4,5)=336*rc5
      aj(5,1)=24*one
      aj(5,2)=120*rc1
      aj(5,3)=360*rc2
      aj(5,4)=840*rc3
      aj(5,5)=1680*rc4
      call gaussj(aj,5,5,bj,1,1)
      alpha=bj(1)
      alpha1=bj(2)
      alpha2=bj(3)
      alpha3=bj(4)
      alpha4=bj(5)
c
c   start iteration loop to find delta(with gamma fixed)
c
      do 550 j=1,200
c
c   generate pseudo wavefunction-note missing factor exp(delta)
c
        do 560 k=1,jrc
          rp=r(k)
          r2=rp*rp
          polyr = r2*(((((alpha4*rp+alpha3)*rp+alpha2)*rp+
     1     alpha1)*rp+ alpha)*r2+gamma)
          ar(k) = rp**lp * exp(polyr)
 560    continue
c
c   integrate pseudo charge density from r = 0 to rc
c
        ll = 2
        cdps = - ar(jrc) * ar(jrc) * rab(jrc)
        if (jrc .ne. 2*(jrc/2)) then
          do 120 k=jrc,1,-1
            cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
            ll = 6 - ll
 120      continue
        else
          do 121 k=jrc,4,-1
            cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
            ll = 6 - ll
 121      continue
          cdps = cdps - ar(4) * ar(4) * rab(4)
          cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1     3 * ar(2) *ar(2) * rab(2) +
     2     3 * ar(3) *ar(3) * rab(3) +
     3     ar(4) * ar(4) * rab(4))/8
        endif
        cdps = cdps/3
c
c   Calculate new delta(with gamma fixed), uses false position
c
        fdnew = log(cdrc/cdps) - 2*delta
        if (abs(fdnew) .lt. errmin) then
          v0pp=8*((2*one*(lp-one)+5*one)*alpha+gamma**2)
          return
        endif
        if (j .eq. 1) then
          ddelta=-pfive
        else
          ddelta = - fdnew * ddelta / (fdnew-fdold)
        endif
        delta = delta + ddelta
        bj(1)=log(arc/rc1**lp)-delta-gamma*rc2
        bj(2)=brc-lp/rc1-2*gamma*rc1
        bj(3)=vrc-ev+(lp/rc1)**2-brc**2-2*gamma
        vt=vrc-ev+lp*(lp-1)/rc2
        bj(4)=vap-2*(vt*brc+lp**2/rc3-brc**3)
        bj(5)=vapp-2*((vap-2*lp*(lp-1)/rc3)*brc+(vt-3*brc**2)*
     1   (vt-brc**2)-3*lp**2/rc4)
        aj(1,1)=rc4
        aj(1,2)=rc5
        aj(1,3)=rc6
        aj(1,4)=rc7
        aj(1,5)=rc8
        aj(2,1)=4*rc3
        aj(2,2)=5*rc4
        aj(2,3)=6*rc5
        aj(2,4)=7*rc6
        aj(2,5)=8*rc7
        aj(3,1)=12*rc2
        aj(3,2)=20*rc3
        aj(3,3)=30*rc4
        aj(3,4)=42*rc5
        aj(3,5)=56*rc6
        aj(4,1)=24*rc1
        aj(4,2)=60*rc2
        aj(4,3)=120*rc3
        aj(4,4)=210*rc4
        aj(4,5)=336*rc5
        aj(5,1)=24*one
        aj(5,2)=120*rc1
        aj(5,3)=360*rc2
        aj(5,4)=840*rc3
        aj(5,5)=1680*rc4
        call gaussj(aj,5,5,bj,1,1)
        alpha=bj(1)
        alpha1=bj(2)
        alpha2=bj(3)
        alpha3=bj(4)
        alpha4=bj(5)
        fdold = fdnew
 550  continue
      write(6,1000)
 1000 format(//, 'error in gamfind - delta not found')
      call ext(860+lp)
      end
C
C
C
      subroutine gaussj(a,n,np,b,m,mp)
c
c ****************************************************************
c *                                                              *
c *  njtj                                                        *
c *    Gauss-Jordan routine used by pseudt to find polynominal   *
c *  constants.  Taken from Numerical Recipes, page 28.          *
c *  njtj                                                        *
c *                                                              *
c ****************************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (nmax=50,zero=0.D0,one=1.D0)
Cray      parameter (nmax=50,zero=0.0,one=1.0)
c
      dimension a(np,np),b(np,mp),ipiv(nmax),indxr(nmax),indxc(nmax)
c
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=zero
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(6,100)
                call ext(890)
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.zero) then
          write(6,100)
          call ext(891)
        endif
        pivinv=one/a(icol,icol)
        a(icol,icol)=one
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=zero
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
 100  format(//,'Singular matrix, stopped in gaussj')
      end
C
C
C
      subroutine input(itype,ikerk,icorr,ispp,zsh,rsh,
     1 nr,a,b,r,rab,nameat,norb,ncore,no,lo,so,zo,
     2 znuc,zel,evi)
c
c  subroutine to read input parameters
c
c  njtj ***  modifications  ***
c    The input and output variables passed have been changed.
c    There are five new pseudopotential generation options
c    The input variables znuc,zsh,rsh,rmax,aa,bb are
c    compared to a small positive value - eliminates
c    floating point comparisions errors(zero is
c    not always zero).
c  njtj ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (one=1.D0,zero=0.D0,pfive=0.5D0)
Cray      parameter (one=1.0,zero=0.0,pfive=0.5)

      character*1 ispp
      character*2 type,icorr,nameat
      character*3 name,kerker
      character*10 iray(5),ititle(5)
c
c  dimension of transfered data
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),
     1 zo(norb),evi(norb)
c
c  dimensions of data used in routine
c
      dimension nc(15),lc(15),nomin(5)
c
c  data for orbitals, 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p
c
      data nc /1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      data lc /0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
c
      do 5 i=1,5
        nomin(i)=10
 5    continue
      do 6 i=1,norb
        no(i)=0
        lo(i)=0
        so(i)=zero
        zo(i)=zero
        evi(i)=zero
 6    continue
c
c  read the type of calculation and title card
c   itype =
c   ae = 0 all electron calculation
c   pg = 1 pseudopotential generation w/o core correction
c   pe = 2 pseudopotential generation w/  core correction exchange
c   ph = 3 pseudopotential generation w/  core correction hartree/exc
c   pt = 4 pseudopotential test
c   pm = 5 pseudopotential test + valence charge modify
c
      read(5,10) type,ititle
 10   format(3x,a2,5a10)
c
c  if type = ' ' , no more data, program ends
c
      if (type .eq. 'ae') then
        itype=0
      elseif (type .eq. 'pg') then
        itype=1
      elseif (type .eq. 'pe') then
        itype=2
      elseif (type .eq. 'ph') then
        itype=3
      elseif (type .eq. 'pt') then
        itype=4
      elseif (type .eq. 'pm') then
        itype=5
      else
        itype=-1
        return
      endif
c
c  njtj  ***  major modification  start  ***
c  There are seven ways to generate the pseudopotential :
c    kerker = van Vanderbilt
c    kerker = tam Troullier and Martins
c    kerker = ker (yes) Kerker
c    kerker = hsc (no)  Hamann Schluter and Chiang
c    kerker = min (oth) datafile made for minimization
c    kerker = bhs Bachelet, Hamann and Schluter
c    kerker = tm2 Improved Troullier and Martins
c
      if (itype.gt.0) then
        read(5,11)kerker
   11 format(8x,a3)
        if(kerker .eq. 'tm2' .or. kerker .eq. 'TM2') then
          ikerk = 6
        elseif(kerker .eq. 'bhs' .or. kerker .eq. 'BHS') then
          ikerk = 5
        elseif(kerker .eq. 'oth' .or. kerker .eq. 'OTH' .or.
     1   kerker .eq. 'min' .or. kerker .eq. 'MIN') then
          ikerk = 4
        elseif (kerker .eq. 'van' .or. kerker .eq.'VAN') then
          ikerk = 3
        elseif (kerker .eq. 'tbk' .or. kerker .eq. 'TBK'
     1   .or. kerker .eq. 'tam' .or. kerker .eq. 'TAM') then
          ikerk = 2
        elseif (kerker .eq. 'yes' .or. kerker .eq. 'YES' .or.
     1   kerker .eq. 'ker' .or. kerker .eq. 'KER') then
          ikerk = 1
        elseif (kerker .eq. 'no ' .or. kerker .eq. ' no' .or.
     1   kerker .eq. 'NO ' .or. kerker .eq. ' NO' .or. kerker
     2   .eq. 'hsc' .or. kerker .eq. 'HSC') then
          ikerk = 0
        else
          write(6,1000)kerker
          call ext(150)
        endif
      endif
 1000 format(//,'error in input - kerker =',a3,' unknown')
c  njtj  ***  major modification end  ***
c
c   read element name and correlation type
c   ispp = ' ' - nonspin calculation
c   ispp = s  - spin polarized calculation
c   ispp = r  - relativistic calculation
c
      read(5,15) nameat,icorr,ispp
 15   format(3x,a2,3x,a2,a1)
      if (ispp .ne. 's' .and. ispp .ne. 'r') ispp=' '
      if (ispp .eq. 's' .and. icorr .eq. 'xa') ispp=' '
      if (ispp .eq. 's' .and. icorr .eq. 'wi') ispp=' '
      if (ispp .eq. 's' .and. icorr .eq. 'hl') ispp=' '
c
c  njtj   ***  major modification start  ***
c   Floating point comparison error modification.
c   Read the atomic number (nuclear charge),
c   shell charge and radius (added to the nuclear potential),
c   and radial grid parameters.
c
      read(5,20) znuc,zsh,rsh,rmax,aa,bb
 20   format(6f10.3)
      if (abs(znuc) .le. 0.00001) znuc=charge(nameat)
      if (itype .lt. 4) then
c
c   set up grid
c
        if (abs(rmax) .lt. 0.00001) rmax=80*one
        if (abs(aa) .lt. 0.00001) aa=6*one
        if (abs(bb) .lt. 0.00001) bb=40*one
        a = exp(-aa)/znuc
        b = 1/bb
        do 30 i=1,nr
          if (i .eq. nr) then
            write(6,50)
            call ext(100)
          endif
          r(i) = a*(exp(b*(i-1))-1)
          rab(i) = (r(i)+a)*b
          if (r(i) .gt. rmax) goto 60
 30     continue
 60     nr = i-1
      endif
 50   format(/,' error in input - arraylimits',
     1 ' for radial array exceeded',/)
c  njtj  ***  major modification end  ***
c
c   read the number of core and valence orbitals
c

      read(5,70) ncore,nval
 70   format(2i5,4f10.3)
      if (ncore .gt. 15) then
        write(6,1010)
        call ext(101)
      endif
 1010 format(//,'error in input - max number of core orbitals',
     1 'is 15')
c
c   compute occupation numbers and orbital energies for the core
c
      zcore = zero
      if (ncore .eq. 0) goto 85
      sc = zero
      if (ispp .ne. ' ') sc=-pfive
      norb = 0
      do 80 i=1,ncore
        do 80 j=1,2
          if (ispp .eq. ' ' .and. j .eq. 2) goto 80
          norb = norb + 1
          no(norb) = nc(i)
          lo(norb) = lc(i)
          so(norb) = sc
          zo(norb) = 2*lo(norb)+1
          if (ispp .eq. ' ') zo(norb) = 2*zo(norb)
          if (ispp .eq. 'r') zo(norb) = 2*(lo(norb)+sc)+1
          zcore = zcore + zo(norb)
          if (abs(zo(norb)) .lt. 0.1) norb=norb-1
          if (ispp .ne. ' ') sc=-sc
 80   continue
      ncore = norb
c
c   for the valence orbitals
c
 85   if (itype .ge. 4) ncore =0
      norb = ncore
      zval = zero
      if (nval .eq. 0) goto 105
      do 90 i=1,nval
        read(5,70) ni,li,zd,zu,evd
        si = zero
        if (ispp .ne. ' ') si=pfive
        do 90 j=1,2
          if (ispp .eq. ' ' .and. j .eq. 2) goto 90
          norb = norb + 1
          if (ispp .ne. ' ') si=-si
          no(norb) = ni
          lo(norb) = li
          so(norb) = si
          zo(norb) = zd+zu
          if (zo(norb) .eq. zero) evi(norb)=evd
          if (ispp .eq. 's') then
            if (si .lt. 0.1) then
              zo(norb) = zd
            else
              zo(norb) = zu
            endif
          elseif (ispp .eq. 'r') then
            zo(norb)=zo(norb)*(2*(li+si)+1)/(4*li+2)
          endif
          zval = zval + zo(norb)
          if (ispp .eq. 'r' .and. li+si .lt. zero) norb=norb-1
          if (norb .eq. 0) goto 90
          if (nomin(lo(norb)+1) .gt. no(norb))
     1     nomin(lo(norb)+1)=no(norb)
 90   continue
c
c   abort if two orbitals are equal
c
      nval = norb - ncore
      do 100 i=1,norb
        do 100 j=1,norb
          if (i .le. j) goto 100
          if (no(i) .ne. no(j)) goto 100
          if (lo(i) .ne. lo(j)) goto 100
          if (abs(so(i)-so(j)) .gt. 0.001) goto 100
          write(6,1020)i
          call ext(110+i)
 100  continue
 1020 format(//,'error in input - orbital ',i2,
     1 'is already occupied')
c
c   reduce n quantum number if pseudoatom
c
      if (itype .ge. 4) then
        do 103 i=1,nval
          no(i) = no(i)-nomin(lo(i)+1)+lo(i)+1
 103    continue
      endif
 105  zion = znuc - zcore - zval
      zel = zval
      if (itype .lt. 4) then
        zel=zel+zcore
      else
        znuc=znuc-zcore
      endif
c
c   find jobname and date and printout, zedate is a machine dependent
c   routine
c
      iray(1)='atom-lda  '
      call zedate(iray(2))
c
c   printout
c
      write(6,110) iray(1),iray(2),ititle
 110  format(1x,a10,a10,5x,5a10,/,21('*'),/)
      if (itype .eq. 0) then
        write(6,120) nameat
      elseif (itype .lt. 4) then
        write(6,121) nameat
      elseif (itype .eq. 4) then
        write(6,124) nameat
      elseif (itype .eq. 5) then
        write(6,125) nameat
      endif
 120  format(1x,a2,' all electron calculation ',/,1x,27('-'),/)
 121  format(1x,a2,' pseudopotential generation',/,1x,29('-'),/)
 124  format(1x,a2,' pseudopotential test',/,1x,23('-'),/)
 125  format(1x,a2,' pseudo test + charge mod ',/,1x,27('-'),/)
      if (ispp .eq. 'r') then
        write(6,150)
 150  format(' r e l a t i v i s t i c ! !',/)
        name = '   '
      elseif (ispp .eq. ' ') then
        name = 'non'
      else
        name = '   '
      endif
      write(6,160) icorr,name
 160  format(' correlation = ',a2,3x,a3,'spin-polarized',/)
      write(6,170) znuc,ncore,nval,zel,zion
 170  format(' nuclear charge             =',f10.6,/,
     1       ' number of core orbitals    =',i3,/,
     2       ' number of valence orbitals =',i3,/,
     3       ' electronic charge          =',f10.6,/,
     4       ' ionic charge               =',f10.6,//)
      if (zsh .gt. 0.00001) write(6,175) zsh,rsh
 175  format(' shell charge =',f6.2,' at radius =',f6.2,//)
      write(6,180)
 180  format(' input data for orbitals',//,
     1 '  i    n    l    s     j     occ',/)
      xji = zero
      do 200 i=1,norb
        if (ispp .eq. 'r') xji = lo(i) + so(i)
        write(6,190) i,no(i),lo(i),so(i),xji,zo(i)
 190  format(1x,i2,2i5,2f6.1,f10.4)
 200  continue
      if (itype .lt. 4) write(6,210) r(2),nr,r(nr),aa,bb
 210  format(//,' radial grid parameters',//,
     1 ' r(1) = .0 , r(2) =',e8.2,' , ... , r(',i3,') =',f6.2,
     2 /,' a =',f5.2,'  b =',f6.2,/)
      return
      end
C
C
C
      subroutine orban(ispp,iorb,ar,br,lmax,nr,a,b,r,rab,
     1 norb,no,lo,zo,so,viod,viou,vid,viu,ev,ek,ep)
c
c  orban is used to analyze and printout data
c  about the orbital.
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (ai=2*137.0360411D0,zero=0.D0)
Cray      parameter (ai=2*137.0360411,zero=0.0)
c
      character*1 ispp
      character*10 name
c
      dimension ar(nr),br(nr),r(nr),rab(nr),no(norb),
     1 lo(norb),zo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),ek(norb),ep(norb)
c
      dimension rzero(10),rextr(10),aextr(10),bextr(10)
c
c      dimension wk1(1000),wk2(1000),wk3(1000),v(1000)
c
      ka = lo(iorb)+1
      lp = ka
      if (so(iorb) .lt. 0.1 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c      compute zeroes and extrema
c
      nzero = 0
      nextr = 0
      rzero(1) = zero
      arp = br(2)
      if (ispp .eq. 'r') then
        if (so(iorb) .lt. 0.1) then
          arp = ka*ar(2)/r(2) + (ev(iorb) - viod(lp,2)/r(2)
     1     - vid(2) + ai*ai) * br(2) / ai
        else
          arp = ka*ar(2)/r(2) + (ev(iorb) - viou(lp,2)/r(2)
     1     - viu(2) + ai*ai) * br(2) / ai
        endif
      endif
      do 20 i=3,nr
        if (nextr .ge. no(iorb)-lo(iorb)) goto 30
        if (ar(i)*ar(i-1) .gt. zero) goto 10
c
c   zero
c
        nzero = nzero + 1
        rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i)) / (ar(i)-ar(i-1))
 10     arpm = arp
        arp = br(i)
        if (ispp .eq. 'r') then
          if ( so(iorb) .lt. 0.1) then
            arp = ka*ar(i)/r(i) + (ev(iorb) - viod(lp,i)/r(i)
     1       - vid(i) + ai*ai) * br(i) / ai
          else
            arp = ka*ar(i)/r(i) + (ev(iorb) - viou(lp,i)/r(i)
     1       - viu(i) + ai*ai) * br(i) / ai
          endif
        endif
        if (arp*arpm .gt. zero) goto 20
c
c   extremum
c
        nextr = nextr + 1
        rextr(nextr) = (arp*r(i-1)-arpm*r(i)) / (arp-arpm)
        aextr(nextr) = (ar(i)+ar(i-1))/2
     1   - (arp**2+arpm**2) * (r(i)-r(i-1)) / (4*(arp-arpm))
        bextr(nextr) = br(i)
 20   continue
c
c   find orbital kinetic and potential energy
c   the potential part includes only the interaction with
c   the nuclear part
c
 30   ek(iorb) = br(1)*br(1)*rab(1)
      ep(iorb) = zero
      sa2 = zero
      lp = lo(iorb)+1
      llp = lo(iorb)*lp
      ll = 2
      if (2*(nr/2) .eq. nr) ll=4
      i90=nr
      i99=nr
      do 40 i=nr,2,-1
        ar2 = ar(i)*ar(i)
        br2 = br(i)*br(i)
        deni = ar2
        if (ispp .eq. 'r') deni=deni+br2
        ek(iorb) = ek(iorb) + ll * (br2 + ar2*llp/r(i)**2)*rab(i)
        if (so(iorb) .lt. 0.1) then
          ep(iorb) = ep(iorb) + ll * deni*viod(lp,i)*rab(i)/r(i)
        else
          ep(iorb) = ep(iorb) + ll * deni*viou(lp,i)*rab(i)/r(i)
        endif
        ll = 6 - ll
        if (sa2 .gt. 0.1) goto 40
        sa2 = sa2 + deni*rab(i)
        if (sa2 .le. 0.01) i99 = i
        i90 = i
 40   continue
      ek(iorb) = ek(iorb) / 3
      ep(iorb) = ep(iorb) / 3
      if (ispp .eq. 'r') ek(iorb) = zero
c
c   printout
c
      write(6,80) no(iorb),lo(iorb),so(iorb)
 80   format(/,' n =',i2,'  l =',i2,'  s =',f4.1)
      name = 'a extr    '
      write(6,100) name,(aextr(i),i=1,nextr)
      name = 'b extr    '
      if (ispp .eq. 'r') write(6,100) name,(bextr(i),i=1,nextr)
      name = 'r extr    '
      write(6,100) name,(rextr(i),i=1,nextr)
      name = 'r zero    '
      write(6,100) name,(rzero(i),i=1,nzero)
      name = 'r 90/99 % '
      write(6,100) name,r(i90),r(i99)
      if (ev(iorb) .eq. zero) then
        if (zo(iorb) .ne. zero) then
          write(6,110)zo(iorb)
        else
          write(6,120)
        endif
      endif
 100  format(8x,a10,2x,8f8.3)
 110  format(8x,'WARNING: This orbital is not bound',
     1 ' and contains ',f6.4,' electrons!!')
 120  format(8x,'WARNING:  This orbital is not bound!')
c
c  njtj  ***  plotting routines  ***
c    Save plotting information to current plot.dat file
c  (unit = 3),  User must specify what orbital
c   is to be saved(or all).
c
c      iorbplot=3
c       ist=1
c       if (ar(nr-80) .lt. 0.0) ist=-1
c       call potrw(ar,r,nr-85,lo(iorb),1,ist)
c       call wtrans(ar,r,nr,rab,lo(iorb),ist,wk1)
c       do 125 i=2,nr
c         v(i)=viod(lo(iorb)+1,i)/r(i)
c 125   continue
c       zion=4
c       call potran(lo(iorb)+1,v,r,nr,zion,wk1,wk2,wk3)
c       call potrv(v,r,nr-120,lo(iorb))
c
c  njtj  ***  user should adjust for their needs  ***
c
       return
       end
C
C
C
      subroutine polcoe(x,y,n,cof)
c
c ************************************************
c *  njtj                                        *
c *  Returns the coefficients of a polynominal.  *
c *  Taken from numerical recipes, page 93.      *
c *  njtj                                        *
c ************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (nmax=10,zero=0.D0,one=1.D0)
Cray      parameter (nmax=10,zero=0.0,one=1.0)
c
      dimension x(n),y(n),cof(n),s(nmax)
      do 11 i=1,n
        s(i)=zero
        cof(i)=zero
11    continue
      s(n)=-x(1)
      do 13 i=2,n
        do 12 j=n+1-i,n-1
          s(j)=s(j)-x(i)*s(j+1)
12      continue
        s(n)=s(n)-x(i)
13    continue
      do 16 j=1,n
        phi=n
        do 14 k=n-1,1,-1
          phi=k*s(k+1)+x(j)*phi
14      continue
        ff=y(j)/phi
        b=one
        do 15 k=n,1,-1
          cof(k)=cof(k)+b*ff
          b=s(k)+x(j)*b
15      continue
16    continue
      return
      end
C
C
C
      subroutine potran(i,vd,r,nr,zion,a,b,c)
c
c ***********************************************************
c *                                                         *
c *    This is a plotting routine; the user should adjust   *
c *  for their own needs.  The potential is fitted with a   *
c *  second degree polynomial, which is muliplied with the  *
c *  appropriate functions and then integrated by parts     *
c *  to find the fourier transform.  The result is then     *
c *  printed to the current plot.dat file (unit=3) for      *
c *  later plotting.  A marker(marker fn#) is placed at     *
c *  the end of each set of data.                           *
c *                                                         *
c ***********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,one=1.D0)
Cray      parameter (zero=0.0,one=1.0)
c
      dimension vd(nr),r(nr),a(nr),b(nr),c(nr),vql(100)
c
c  The potential times r is fitted to the polynominal
c  a + bx + cx^2 at every other point.
c
      rm=zero
      vm=2*zion
      do 130 k=2,nr,2
        r0=r(k)
        v0=r0*vd(k)+2*zion
        rp=r(k+1)
        vp=rp*vd(k+1)+2*zion
        d1=1/((rp-rm)*(r0-rm))
        d2=1/((rp-r0)*(rm-r0))
        d3=1/((r0-rp)*(rm-rp))
        a(k)=vm*d1+v0*d2+vp*d3
        b(k)=-vm*(r0+rp)*d1-v0*(rm+rp)*d2-vp*(rm+r0)*d3
        c(k)=vm*r0*rp*d1+v0*rm*rp*d2+vp*rm*r0*d3
        rm=rp
        vm=vp
 130  continue
c
c  Find the fourier transform q^2/4pi/zion*vql. Everything is
c  rescaled  by zion.
c
      do 150 j=1,94
        q=one/4*j
        q2=q*q
        vql(j)=zero
        rm=zero
        do 140 k=2,nr-1,2
          rp=r(k+1)
          vql(j)=vql(j)+(2*a(k)*rp+b(k))/q*sin(q*rp)
     1     -((a(k)*rp+b(k))*rp+c(k)-2*a(k)/q2)*cos(q*rp)
     2     -(2*a(k)*rm+b(k))/q*sin(q*rm)
     3     +((a(k)*rm+b(k))*rm+c(k)-2*a(k)/q2)*cos(q*rm)
          rm=rp
 140    continue
        vql(j)=vql(j)/2/zion-one
 150  continue
c
c  Print out the transforms( really q^2/(4pi*zion)*v(q) ) to
c  the current plot.dat file (unit=3) for latter plotting.
c
      do 170 j=1,48
        write(3,6000)one/4*j,vql(j)
 170  continue
      write(3,6008)i
      return
c
c  format statements
c
 6000 format(1x,f7.4,3x,f10.6)
 6008 format(1x,'marker fn',i1)
c
      end
C
C
C
      subroutine potrv(vd,r,nr,k)
c
c ***********************************************************
c *                                                         *
c *    This is a plotting routine; the user should          *
c *  adjust for their own needs.  Prints                    *
c *  out the potential to the current plot.dat              *
c *  file (unit=3) for later ploting.  A marker (marker)    *
c *  is placed at the end of each group of data.            *
c *                                                         *
c ***********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      character*3 marker
c
      dimension vd(nr),r(nr)
c
c  Step size of 0.05 is adjustable as seen fit to give
c  a reasonalble plot.
c
      step=0.0
      do 150,j=5,nr
        if (r(j) .ge. step) then
          write(3,6000)r(j),vd(j)
          step=step+0.05
        endif
 150  continue
      if (k .eq. 0) then
        marker='vns'
      elseif (k .eq. 1) then
        marker='vnp'
      elseif (k .eq. 2) then
        marker='vnd'
      elseif (k .eq. 3) then
        marker='vnf'
      elseif (k .eq. 4) then
        marker='vng'
      endif
      write(3,6001)marker
      return
c
c  Format statements
c
 6000 format(1x,f7.4,3x,f10.5)
 6001 format(1x,'marker ',a3)
      end
C
C
C
      subroutine potrw(vd,r,nr,k,kj,ist)
c
c ***********************************************************
c *                                                         *
c *    This is a plotting routine; the user should          *
c *  adjust/eliminatebfor their own needs.  Prints          *
c *  out the wave functions to the current plot.dat         *
c *  file (unit=3) for later ploting.  A marker (marker)    *
c *  is placed at the end of each group of data.            *
c *                                                         *
c ***********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter statement
c  ###    to single precision statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
      parameter (zero=0.D0,pzf=0.05D0)
Cray      parameter (zero=0.0,pzf=0.05)
c
c
      character*3 marker
c
      dimension vd(nr),r(nr)
c
c  Step size of 0.05 is adjustable as seen fit to give
c  a reasonalble plot.
c
      step=zero
      do 150,j=2,nr
        if (r(j) .ge. step) then
          write(3,6000)r(j),vd(j)*ist
          step=step+pzf
        endif
 150  continue
      if (kj .eq. 0) then
        if (k .eq. 0) then
          marker='wsp'
        elseif (k .eq. 1) then
          marker='wpp'
        elseif (k .eq. 2) then
          marker='wdp'
        elseif (k .eq. 3) then
          marker='wfp'
        elseif (k .eq. 4) then
          marker='wgp'
        endif
      else
        if (k .eq. 0) then
          marker='wst'
        elseif (k .eq. 1) then
          marker='wpt'
        elseif (k .eq. 2) then
          marker='wdt'
        elseif (k .eq. 3) then
          marker='wft'
        elseif (k .eq. 4) then
          marker='wgt'
        endif
      endif
      write(3,6001)marker
      return
c
c  Format statements
c
 6000 format(1x,f7.4,3x,f18.14)
 6001 format(1x,'marker ',a3)
      end
C
C
C
      subroutine prdiff(nconf,econf)
c
c   Prints out the energy differences between
c   different atomic configurations.
c
c   njtj  ***  modifications  ***
c     econf is able to handle larger numbers
c     of configurations.
c   njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      dimension econf(100)
c
      write(6,10) (i,i=1,nconf)
      do 30 i=1,nconf
        write(6,20) i,(econf(i)-econf(j),j=1,i)
 30   continue
 10   format(/,' total energy difference',//,2x,9i9)
 20   format(1x,i2,1x,9f9.4)
      return
      end
C
C
C
      subroutine pseud2(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2 viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3 wk3,wk4,wk5,wk6,wk7,nops,v,ar,br,wkb,evi)
c
c *************************************************************
c *                                                           *
c *     This routine was written by Norman J. Troullier Jr.   *
c *   April 1990, while at the U. of Minnesota, all           *
c *   comments concerning this routine should be directed     *
c *   to him.                                                 *
c *                                                           *
c *     troullie@128.101.224.101                              *
c *     troullie@csfsa.cs.umn.edu                             *
c *     612 625-0392                                          *
c *                                                           *
c *     pseud2 generates a pseudopotential using the          *
c *   improved scheme of N. Troullier and J. L. Martins.      *
c *   The general format of this routine is the same as the   *
c *   pseudo and pseudk routines.  Output/input is            *
c *   compatible.                                             *
c *                                                           *
c *************************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,one=1.D0,tpfive=2.5D0,ecuts=1.0D-3)
      parameter (small=1.D-12,pnine=0.9D0,ai=2*137.0360411D0,sml=0.1D0)
Cray      parameter (zero=0.0,one=1.0,tpfive=2.5,ecuts=1.0E-3)
Cray      parameter (small=1.E-12,pnine=0.9,ai=2*137.0360411,sml=0.1)
c
      character*1 ispp,blank,il(5)
      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 iray(6),ititle(7)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),vod(nr),vou(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),wk5(nr),wk6(nr),wk7(nr),
     4 wkb(3*nr),nops(norb),v(nr),ar(nr),br(nr),evi(norb)
c
      dimension indd(5),indu(5),rc(5),rcut(10),vstore(1000),
     1 etot(10),aa(7),rr(7),coe(7),aj(5,5),bj(5)
c
      data il/'s','p','d','f','g'/
      if (ncore .eq. norb) return
      ifcore = itype-1
      pi = 4*atan(one)
      do 3 i=1,5
        indd(i)=0
        indu(i)=0
 3    continue
      do 4 i=1,40
        nops(i) = 0
 4    continue
c
c  read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be
c  set to greater then zero.
c
      read(5,10) (rc(i),i=1,5),cfac,rcfac
 10   format(7f10.5)
      if (cfac .eq. 0.D0) cfac=one
c
c  Reset vod and vou to zero,
c  they are here used to store the pseudo valence charge density.
c
      do 15 i=1,nr
        vod(i) = zero
 15   continue
      do 16 i=1,nr
        vou(i) = zero
 16   continue
c
c  print heading
c
      write(6,20) nameat
 20   format(//,1x,a2,' pseudopotential generation using the ',
     1 'Improved Troullier and Martins method',/,1x,60('-'),//,
     2 ' nl    s    eigenvalue',6x,'rc',10x,'cdrc',7x,'delta',/)
c
c  Start loop over valence orbitals, only one orbital for each
c  angular momentum and spin can exist.
c
      ncp = ncore+1
      do 190 i=ncp,norb
        lp = lo(i) + 1
        llp = lo(i)*lp
        if (so(i) .lt. 0.1) then
          if (indd(lp) .ne. 0) then
            write(6,1000)lp-1
            call ext(800+lp)
          else
            indd(lp) = i
          endif
        else
          if (indu(lp) .ne. 0) then
            write(6,1010)lp-1
            call ext(810+lp)
          else
            indu(lp) = i
          endif
        endif
 1000 format(//,'error in pseud2 - two down spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
 1010 format(//,'error in pseud2 - two up spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
c
c  Find the all electron wave function.
c
        do 29 j=1,nr
          ar(j) = zero
 29     continue
        if (so(i) .lt. 0.1) then
          do 30 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 30       continue
        else
          do 31 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 31       continue
        endif
        if (ispp .ne. 'r') then
          do 32 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 32       continue
        endif
c
c  The parameter iflag has been added as a nonconvegence
c  indicator for auxillary routines.  Its value does
c  not change its operation.  iflag is a returned value,
c  set to 1 for none convergence.
c
        if (ispp .ne. 'r') then
          iflag=0
          call difnrl(0,i,v,ar,br,lmax,nr,a,b,
     1     r,rab,norb,no,lo,so,znuc,viod,viou,
     2     vid,viu,ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(0,i,v,ar,br,lmax,nr,a,b,r,
     1     rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2     ev,wk1,wk2,wk3,wk4,evi)
         endif
c
c  Find last zero and extremum
c
        ka = lo(i)+1
        if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
        nextr = no(i)-lo(i)
        rzero = zero
        arp = br(2)
c
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp = ka*ar(2)/r(2) + (ev(i) - viod(lp,2)/r(2)
     1       - vid(2) + ai*ai) * br(2) / ai
          else
            arp = ka*ar(2)/r(2) + (ev(i) - viou(lp,2)/r(2)
     1       - viu(2) + ai*ai) * br(2) / ai
          endif
        endif
c
        do 40 j=3,nr-7
          if (nextr .eq. 0) goto 50
          if (ar(j-1)*ar(j) .le. zero .and. evi(i) .eq. zero)
     1     rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
          arpm = arp
          arp = br(j)
c
          if (ispp .eq. 'r') then
            if(so(i) .lt. 0.1) then
              arp = ka*ar(j)/r(j) + (ev(i) - viod(lp,j)/r(j)
     1         - vid(j) + ai*ai) * br(j) / ai
            else
              arp = ka*ar(j)/r(j) + (ev(i) - viou(lp,j)/r(j)
     1         - viu(j) + ai*ai) * br(j) / ai
            endif
          endif
c
          if (arp*arpm .gt. zero) goto 40
          rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
          nextr = nextr - 1
 40     continue
 50     if (rzero .lt. r(2)) rzero = r(2)
c
c  Check rc if inside rzero,
c  reset to .9 between rmax and rzero if inside
c  if rc(lp) is negative, rc(lp) is percent of way
c  betweeen rzero and rmax.
c
        if (rc(lp) .gt. rzero) then
        elseif(rc(lp) .ge. zero) then
          rc(lp) = rzero + pnine*(rextr-rzero)
        else
          rc(lp) = rzero - rc(lp)*(rextr-rzero)
        endif
c
c  Find the index for odd grid point closest to rc.
c
        do 70 j=1,nr
          if (r(j) .gt. rc(lp)) goto 80
 70     continue
 80     jrc=j-1
        rc(lp)=r(jrc)
c
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c
        ist=1
        if (ar(jrc) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),1,ist)
        do 41 j=1,nr
          ar(j)=ar(j)*ist
          br(j)=br(j)*ist
 41     continue
c
c  njtj  ***  user should adjust for their needs  ***
c
c
c  Reset n quantum numbers.
c
        nops(i) = lp
c
c  Find the integrated charge inside rc(1-charge outside).
c
        ll = 2
        if (ispp .eq. 'r') then
          cdrc = -(ar(jrc)*ar(jrc)+br(jrc)*br(jrc))*rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 102 k=jrc,1,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 102        continue
          else
            do 103 k=jrc,4,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 103        continue
            cdrc = cdrc-(ar(4)*ar(4)+br(4)*br(4))*rab(4)
            cdrc = cdrc+9*((ar(1)*ar(1)+br(1)*br(1))*rab(1)+
     1       3*(ar(2)*ar(2)+br(2)*br(2))*rab(2)+
     2       3*(ar(3)*ar(3)+br(3)*br(3))*rab(3)+
     3       (ar(4)*ar(4)+br(4)*br(4))*rab(4))/8
          endif
          cdrc = cdrc/3
        else
          cdrc = - ar(jrc) * ar(jrc) * rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 100 k=jrc,1,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 100        continue
          else
            do 101 k=jrc,4,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 101        continue
            cdrc = cdrc - ar(4) * ar(4) * rab(4)
            cdrc = cdrc + 9 * ( ar(1) * ar(1) * rab(1) +
     1       3 * ar(2) *ar(2) * rab(2) +
     2       3 * ar(3) *ar(3) * rab(3) +
     3       ar(4) * ar(4) * rab(4))/8
          endif
          cdrc = cdrc/3
        endif
c
c  Find the values for wave(arc), d(wave)/dr(arp), potential(vrc),
c  d(potential)/dr(vrp), and d2(potential)/dr2(vrpp)
c
        rc1 = r(jrc)
        rc2 = rc1 * rc1
        rc3 = rc2 * rc1
        rc4 = rc2 * rc2
        rc5 = rc4 * rc1
        rc6 = rc4 * rc2
        rc7 = rc4 * rc3
        rc8 = rc4 * rc4
        rc9 = rc4 * rc5
        rc10= rc4 * rc6
        arc = ar(jrc)
        arp = br(jrc)
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp=ka*ar(jrc)/r(jrc) + (ev(i) - viod(lp,jrc)/r(jrc)
     1       - vid(jrc) + ai*ai) * br(jrc)/ai
          else
            arp=ka*ar(jrc)/r(jrc) + (ev(i) - viou(lp,jrc)/r(jrc)
     1       - viu(jrc) + ai*ai) * br(jrc)/ai
          endif
        endif
        arp =arp
        brc = arp / arc
c
        if (so(i) .lt. 0.1) then
          vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
          aa(1)=viod(lp,jrc-3)/r(jrc-3) + vid(jrc-3)
          aa(2)=viod(lp,jrc-2)/r(jrc-2) + vid(jrc-2)
          aa(3)=viod(lp,jrc-1)/r(jrc-1) + vid(jrc-1)
          aa(4)=vrc
          aa(5)=viod(lp,jrc+1)/r(jrc+1) + vid(jrc+1)
          aa(6)=viod(lp,jrc+2)/r(jrc+2) + vid(jrc+2)
          aa(7)=viod(lp,jrc+3)/r(jrc+3) + vid(jrc+3)
       else
          vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
          aa(1)=viou(lp,jrc-3)/r(jrc-3) + viu(jrc-3)
          aa(2)=viou(lp,jrc-2)/r(jrc-2) + viu(jrc-2)
          aa(3)=viou(lp,jrc-1)/r(jrc-1) + viu(jrc-1)
          aa(4)=vrc
          aa(5)=viou(lp,jrc+1)/r(jrc+1) + viu(jrc+1)
          aa(6)=viou(lp,jrc+2)/r(jrc+2) + viu(jrc+2)
          aa(7)=viou(lp,jrc+3)/r(jrc+3) + viu(jrc+3)
        endif
        rr(1)=r(jrc-3)-r(jrc)
        rr(2)=r(jrc-2)-r(jrc)
        rr(3)=r(jrc-1)-r(jrc)
        rr(4)=zero
        rr(5)=r(jrc+1)-r(jrc)
        rr(6)=r(jrc+2)-r(jrc)
        rr(7)=r(jrc+3)-r(jrc)
        call polcoe(rr,aa,7,coe)
        vap   = coe(2)
        vapp  = 2*coe(3)
c
c   Set up matrix without the d2(potential(0)/dr2=0 condition
c   to find an intial guess for gamma.
c
        delta=zero
        bj(1)=log(arc/rc1**lp)
        bj1=bj(1)
        bj(2)=brc-lp/rc1
        bj2=bj(2)
        bj(3)=vrc-ev(i)-2*lp/rc1*bj2-bj2**2
        bj3=bj(3)
        bj(4)=vap+2*lp/rc2*bj2-2*lp/rc1*bj3-2*bj2*bj3
        bj4=bj(4)
        bj(5)=vapp-4*lp/rc3*bj2+4*lp/rc2*bj3-2*lp/rc1*bj4-2*bj3**2
     1   -2*bj2*bj4
        bj5=bj(5)
        aj(1,1)=rc2
        aj(1,2)=rc4
        aj(1,3)=rc6
        aj(1,4)=rc8
        aj(1,5)=rc10
        aj(2,1)=2*rc1
        aj(2,2)=4*rc3
        aj(2,3)=6*rc5
        aj(2,4)=8*rc7
        aj(2,5)=10*rc9
        aj(3,1)=2*one
        aj(3,2)=12*rc2
        aj(3,3)=30*rc4
        aj(3,4)=56*rc6
        aj(3,5)=90*rc8
        aj(4,1)=zero
        aj(4,2)=24*rc1
        aj(4,3)=120*rc3
        aj(4,4)=336*rc5
        aj(4,5)=720*rc7
        aj(5,1)=zero
        aj(5,2)=24*one
        aj(5,3)=360*rc2
        aj(5,4)=1680*rc4
        aj(5,5)=5040*rc6
        call gaussj(aj,5,5,bj,1,1)
        gamma=bj(1)
        alpha=bj(2)
        alpha1=bj(3)
        alpha2=bj(4)
        alpha3=bj(5)
c
c  Start iteration loop to find delta, uses false postion.
c
        do 150 j=1,50
c
c  Generate pseudo wavefunction-note missing factor exp(delta).
c
          do 110 k=1,jrc
            rp=r(k)
            r2=rp*rp
            polyr = r2*((((alpha3*r2+alpha2)*r2+
     1       alpha1)*r2+ alpha)*r2+gamma)
            ar(k) = rp**lp * exp(polyr)
 110      continue
c
c  Integrate pseudo charge density from r = 0 to rc.
c
          ll = 2
          cdps = - ar(jrc) * ar(jrc) * rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 120 k=jrc,1,-1
              cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 120        continue
          else
            do 121 k=jrc,4,-1
              cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 121        continue
            cdps = cdps - ar(4) * ar(4) * rab(4)
            cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1       3 * ar(2) *ar(2) * rab(2) +
     2       3 * ar(3) *ar(3) * rab(3) +
     3       ar(4) * ar(4) * rab(4))/8
          endif
          cdps = cdps/3
c
c   Calculate new delta
c
          fdnew = log(cdrc/cdps) - 2*delta
          if (abs(fdnew) .lt. small) goto 160
          if (j .eq. 1) then
            ddelta=-one/2
          else
            ddelta = - fdnew * ddelta / (fdnew-fdold)
          endif
          delta = delta + ddelta
          bj(1)=bj1-delta
          bj(2)=bj2
          bj(3)=bj3
          bj(4)=bj4
          bj(5)=bj5
          aj(1,1)=rc2
          aj(1,2)=rc4
          aj(1,3)=rc6
          aj(1,4)=rc8
          aj(1,5)=rc10
          aj(2,1)=2*rc1
          aj(2,2)=4*rc3
          aj(2,3)=6*rc5
          aj(2,4)=8*rc7
          aj(2,5)=10*rc9
          aj(3,1)=2*one
          aj(3,2)=12*rc2
          aj(3,3)=30*rc4
          aj(3,4)=56*rc6
          aj(3,5)=90*rc8
          aj(4,1)=zero
          aj(4,2)=24*rc1
          aj(4,3)=120*rc3
          aj(4,4)=336*rc5
          aj(4,5)=720*rc7
          aj(5,1)=zero
          aj(5,2)=24*one
          aj(5,3)=360*rc2
          aj(5,4)=1680*rc4
          aj(5,5)=5040*rc6
          call gaussj(aj,5,5,bj,1,1)
          gamma=bj(1)
          alpha=bj(2)
          alpha1=bj(3)
          alpha2=bj(4)
          alpha3=bj(5)
          fdold = fdnew
 150    continue
c
c  End iteration loop for delta.
c
        write(6,1020)lp-1
        call ext(820+lp)
 1020 format(//,'error in pseud2 - nonconvergence in finding',
     1 /,' starting delta for angular momentum ',i1)
c
c  Bracket the correct gamma, use gamma and -gamma
c  from above as intial brackets, expands brackets
c  until a root is found..
c
 160    alpha4=zero
        x1=gamma
        x2=-gamma
c
        call zrbac2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1   rc8,lp,arc,brc,vrc,vap,vapp,ev(i),cdrc,r,rab,
     2   jrc,delta,gamma,alpha,alpha1,alpha2,alpha3,
     3   alpha4,ar)
c
c  Iteration loop to find correct gamma, uses
c  bisection to find gamma.
c
        call rtbis2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1   rc8,lp,arc,brc,vrc,vap,vapp,ev(i),cdrc,r,rab,jrc,delta,
     2   gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c  Augment charge density and invert schroedinger equation
c  to find new potential.
c
 645    expd = exp(delta)
        if (so(i) .lt. 0.1) then
          do 169 j=1,jrc
            r2=r(j)*r(j)
            poly=r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*
     1       r2+alpha)*r2+gamma)
            ar(j) = r(j)**lp * expd * exp(poly)
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
            xlamda=((((12*alpha4*r2+10*alpha3)*r2+8*alpha2)*r2+
     1       6*alpha1)*r2+4*alpha)*r2+2*gamma
            vj = ev(i) + xlamda * (2 * lp + xlamda * r2)
     1       +((((132*alpha4*r2+90*alpha3)*r2+56*alpha2)*r2+30*alpha1)*
     2       r2+12*alpha)*r2+2*gamma
            viod(lp,j) = (vj-vid(j)) * r(j)
 169      continue
          do 168 j=jrc+1,nr
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
 168      continue
        else
          do 170 j=1,jrc
            r2=r(j)*r(j)
            poly=r2*(((((alpha4*r2+alpha3)*r2+alpha2)*r2+alpha1)*
     1       r2+alpha)*r2+gamma)
            ar(j) = r(j)**lp * expd * exp(poly)
c
c bug fix Alberto Garcia 5/11/90
c
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
            xlamda=((((12*alpha4*r2+10*alpha3)*r2+8*alpha2)*r2+
     1       6*alpha1)*r2+4*alpha)*r2+2*gamma
            vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1       +((((132*alpha4*r2+90*alpha3)*r2+56*alpha2)*r2+30*alpha1)*
     2       r2+12*alpha)*r2+2*gamma
            viou(lp,j) = (vj-viu(j)) * r(j)
 170      continue
          do 171 j=jrc+1,nr
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
 171      continue
        endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
        ist=1
        call potrw(ar,r,nr-85,lo(i),0,ist)
        if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
        call wtrans(ar,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c
        write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cdrc,delta
 180  format(1x,i1,a1,f6.1,5f12.6)
 190  continue
c
c  End loop over valence orbitals.
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 200 i=ncp,norb
        nops(i) = lo(i) + 1
        zval = zval + zo(i)
 200  continue
      zion = zval+znuc-zel
      if (zval .ne. zero) zratio=zion/zval
      do 210 i=1,nr
        cdd(i) = vod(i)
 210  continue
      do 211 i=1,nr
        cdu(i) = vou(i)
 211  continue
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used
c
      if (ifcore .ne. 0) then
        ac = zero
        bc = zero
        icore = 1
        if (cfac .le. zero .or. zratio .eq. zero) then
          write(6,280) r(icore),ac,bc
        else
          if (rcfac .le. zero) then
            do 220 i=nr,2,-1
              if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220        continue
          else
            do 221 i=nr,2,-1
              if (r(i) .le. rcfac ) goto 230
 221        continue
          endif
 230      icore = i
          cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
          tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
          rbold = tpfive
          do 240 i=1,50
            rbnew = pi+atan(tanb*rbold)
            if (abs(rbnew-rbold) .lt. .00001) then
              bc = rbnew / r(icore)
              ac = cdc(icore) / (r(icore)*sin(rbnew))
              do 260 j=1,icore
                cdc(j) = ac*r(j)*sin(bc*r(j))
 260          continue
              write(6,280) r(icore),ac,bc
              goto 290
            else
              rbold=rbnew
            endif
 240      continue
          write(6,1030)
          call ext(830)
        endif
      endif
 280  format(//,' core correction used',/,
     1 ' pseudo core inside r =',f6.3,/,' ac =',f6.3,' bc =',f6.3,/)
 1030 format(//,' error in pseud2 - noncovergence in finding ',
     1 /,'pseudo-core values')
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
 290  if (ispp .eq. 's') then
        blank='s'
      else
        blank=' '
      endif
      zval2=zval
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
      if (ifcore .eq. 2) zion=zion+zval-zval2
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      ecut=ecuts
      do 315 i=ncp,norb
        lp = lo(i)+1
        if (so(i) .lt. 0.1) then
          do 300 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 300      continue
          rcut(i-ncore) = r(jcut)
          do 310 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 310      continue
          do 311 j=2,nr
            v(j) = viod(lp,j)/r(j)
 311      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        else
          do 312 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 312      continue
          rcut(i-ncore) = r(jcut)
          do 313 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 313      continue
          do 314 j=2,nr
            v(j) = viou(lp,j)/r(j)
 314      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-110,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        endif
 315  continue
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a
c   marker 'zio' for latter plotting
c
      write(3,4559)
      write(3,4560) zion
 4559 format(1x,'marker zio')
 4560 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (ispp .eq. 's' ) then
        do 500 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 505 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 505        continue
          else
            do 506 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 506        continue
          endif
 500    continue
      endif
c
      do 320 i=2,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 320  continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,blank,ifcore,lmax,nr,a,b,r,rab,
     1 norb-ncore,0,nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 znuc,cdd,cdu,cdc,viod,viou,vid,viu,ev(ncp),ek(ncp),
     3 ep(ncp),wk1,wk2,wk3,wk4,wk5,wk6,wk7,evi(ncp))
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,326) (ev(i),i=ncp,norb)
 325  format(//,' test of eigenvalues',//,' rcut =',8(2x,a1,f7.2))
 326  format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,330)
 330  format(///,' l    vps(0)    vpsmin      at r',/)
      do 370 i=1,lmax
        if (indd(i)+indu(i) .eq. 0) goto 370
        if (indd(i) .ne. 0) then
          vpsdm = zero
          do 350 j=2,nr
            if (r(j) .lt. .00001) goto 350
            vps = viod(i,j)/r(j)
            if (vps .lt. vpsdm) then
              vpsdm = vps
              rmind = r(j)
            endif
 350      continue
          write(6,360) il(i),viod(i,2)/r(2),vpsdm,rmind
        endif
        if (indu(i) .ne. 0) then
          vpsum = zero
          do 351 j=2,nr
            if (r(j) .lt. .00001) goto 351
            vps = viou(i,j)/r(j)
            if (vps .lt. vpsum) then
              vpsum = vps
              rminu = r(j)
            endif
 351      continue
          write(6,360) il(i),viou(i,2)/r(2),vpsum,rminu
        endif
 360  format(1x,a1,3f10.3)
 370  continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,
     1 nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1) = 'atom-lda  '
      call zedate(iray(2))
      iray(3) = '  Improved'
      iray(4) = ' Troullier'
      iray(5) = ' - Martins'
      iray(6) = ' potential'
c
c  Encode the title array.
c
      do 390 i=1,7
        ititle(i) = '          '
 390  continue
      do 420 i=1,lmax
        if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
        zelu = zero
        zeld = zero
        if (indd(i) .ne. 0) then
          noi = no(indd(i))
          zeld = zo(indd(i))
        endif
        if (indu(i) .ne. 0) then
          noi = no(indu(i))
          zelu = zo(indu(i))
        endif
        zelt = zeld + zelu
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt
         write(ititle(2*i),401)ispp,rc(i)
 400     format(i1,a1,'(',f6.2,')')
 401     format(a1,' rc=',f5.2)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld
         write(ititle(2*i),411)zelu,ispp,rc(i)
 410     format(i1,a1,'  (',f4.2,',')
 411     format(f4.2,')',a1,f4.2)
        endif
 420  continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
        if (indu(1) .eq. 0) goto 429
        indd(1)=indu(1)
        indu(1)=0
        do 428 j=2,nr
          viod(1,j) = viou(1,j)
          viou(1,j) = zero
 428    continue
 429    do 431 i=2,lmax
          if (indd(i) .eq. 0 .or. indu(i) .eq. 0) goto 431
          do 430 j=2,nr
            viodj = viod(i,j)
            viouj = viou(i,j)
            viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
            viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430      continue
 431    continue
      endif
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 450 i=1,lmax
        if (indd(i) .ne. 0) npotd=npotd+1
        if (indu(i) .ne. 0) npotu=npotu+1
 450  continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
        if (ifull .eq. 0) then
          nicore = 'pcec'
        else
          nicore = 'fcec'
        endif
      elseif (ifcore .eq. 2) then
        if (ifull .eq. 0) then
          nicore = 'pche'
        else
          nicore = 'fche'
        endif
      else
        nicore = 'nc  '
      endif
      if (ispp .eq. 's') then
        irel='isp'
      elseif (ispp .eq. 'r') then
        irel='rel'
      else
        irel = 'nrl'
      endif
      rewind 1
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat
c  file (unit=1).
c
      do 460 i=1,lmax
        if (indd(i) .eq. 0) goto 460
        write(1) i-1,(viod(i,j),j=2,nr)
 460  continue
      do 465 i=1,lmax
        if (indu(i) .eq. 0) goto 465
        write(1) i-1,(viou(i,j),j=2,nr)
 465  continue
c
c  Write the charge densities to the current pseudo.dat
c  file (unit=1).
c
      if (ifcore .eq. 0) then
        write(1) (zero,i=2,nr)
      else
        write(1) (cdc(i),i=2,nr)
      endif
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
      end
C
C
C
      subroutine pseudb(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2 viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,wk3,
     3 wk4,wk5,wk6,wk7,f,g,nops,v,ar,br,arps,wkb,evi)
c
c *************************************************************
c *                                                           *
c *    pseudo generates the pseudo potential using            *
c *  the scheme of Bachelet, Hamann, and Schluter -           *
c *  Phys. Rev. B. 26, 4199.                                  *
c *                                                           *
c *************************************************************
c
c  njtj  *** modifications  ***
c    The only major modifications are in the spin-polarized
c    treatment of the el-el unscreening of the pseudopotential
c    A spin-polarized pseudopotential is unscreened
c    with a spin-polarized valence charge.  This was not done
c    in pseudo or pseudok in earlier versions of this
c    program.
c  njtj  *** modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter(zero=0.D0,ecuts=1.0D-3,tpfive=2.5D0,one=1.D0)
      parameter(small=1.D-13,small2=1.D-10,small3=1.D-18,pzfive=.05D0)
      parameter(pfive=0.5D0,small4=1.D-6,ai=2*137.0360411D0)
Cray       parameter(zero=0.0,ecuts=1.0E-3,tpfive=2.5,one=1.0)
Cray       parameter(small=1.E-13,small2=1.E-10,small3=1.E-18,pzfive=.05)
Cray       parameter(pfive=0.5,small4=1.E-6,ai=2*137.0360411)
c
      character*1 ispp,blank,il(5)
      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 ititle(7),iray(6)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),vod(nr),vou(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),wk5(nr),wk6(nr),wk7(nr),
     4 wkb(6*nr),f(nr),g(nr),nops(norb),v(nr),
     5 ar(nr),br(nr),arps(nr),evi(norb)
c
      dimension etot(10),indd(5),indu(5),rc(5),rcut(10)
c
      data il/'s','p','d','f','g'/
      do 3 i=1,5
        indd(i)=0
        indu(i)=0
 3    continue
      if (ncore .eq. norb) return
      if (itype .ne. 1 .and. itype .ne. 2 .and. itype .ne. 3) return
      ifcore = itype - 1
      pi = 4*atan(one)
c
c  Spin-polarized potentails should be unscreened with
c  a spin-polarized valence charge.  This was not
c  done in pseudo and pseudk in earlier versions
c  of this program.
c
      if (ispp .eq. 's' ) then
        blank = 's'
      else
        blank = ' '
      endif
c
c  read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be
c  set to greater then zero.
c
      read(5,10) (rc(i),i=1,5),cfac,rcfac
 10   format(7f10.5)
      if (cfac .eq. zero) cfac=one
c
c   Reset vod and vou to zero.  They are here used to store
c   the pseudo valence charge density.
c
      do 15 i=1,nr
        vod(i) = zero
        vou(i) = zero
 15   continue
c
c  Print the heading.
c
      write(6,20) nameat
 20   format(//,a2,' Pseudopotential BHS generation',/,1x,35('-'),//,
     1 ' nl    s    eigenvalue',6x,'rc',4x,6x,'cl',9x,'gamma',
     2 7x,'delta',/)
c
c      start loop over valence orbitals
c
      ncp = ncore+1
      do 190 i=ncp,norb
        lp = lo(i) + 1
        llp = lo(i)*lp
        if (so(i) .lt. 0.1) then
          if (indd(lp) .ne. 0) then
            write(6,1000)lp-1
            call ext(800+lp)
          else
            indd(lp) = i
          endif
        else
          if (indu(lp) .ne. 0) then
            write(6,1010)lp-1
            call ext(810+lp)
          else
            indu(lp) = i
          endif
        endif
 1000 format(//,'error in pseudb - two down spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
 1010 format(//,'error in pseudb - two up spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
c
c      find all electron wave function
c
        do 25 j=1,nr
          ar(j)=zero
 25     continue
        if (so(i) .lt. 0.1) then
          do 27 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 27       continue
        else
          do 30 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 30       continue
        endif
        if (ispp .ne. 'r') then
          do 32 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 32       continue
        endif
        if (ispp .ne. 'r') then
          call difnrl(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,wk1,wk2,wk3,wk4,evi)
        endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c
        ist=1
        if (ar(nr-85) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),1,ist)
c
c  njtj  ***  user should adjust for their needs  ***
c
c  Find the last zero and extremum.
c
        ka = lo(i)+1
        if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
        nextr = no(i)-lo(i)
        rzero = zero
        arp = br(2)
c
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp = ka*ar(2)/r(2) + (ev(i) - viod(lp,2)/r(2)
     1       - vid(2) + ai*ai) * br(2) / ai
          else
            arp = ka*ar(2)/r(2) + (ev(i) - viou(lp,2)/r(2)
     1       - viu(2) + ai*ai) * br(2) / ai
          endif
        endif
c
        do 40 j=3,nr-7
          if (nextr .eq. 0) goto 50
          if (ar(j-1)*ar(j) .le. zero)
     1     rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
          arpm = arp
          arp = br(j)
c
          if (ispp .eq. 'r') then
            if (so(i) .lt. 0.1) then
              arp = ka*ar(j)/r(j) + (ev(i) - viod(lp,j)/r(j)
     1         - vid(j) + ai*ai) * br(j) / ai
            else
              arp = ka*ar(j)/r(j) + (ev(i) - viou(lp,j)/r(j)
     1         - viu(j) + ai*ai) * br(j) / ai
            endif
          endif
c
          if (arp*arpm .gt. zero) goto 40
          rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
          nextr = nextr - 1
 40     continue
c
c  Check rc, if outside bounds reset.
c
 50     if (rzero .lt. r(2)) rzero = r(2)
        if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) goto 60
        if (rc(lp) .ge. rzero) then
          write(6,2001)rc(lp),rextr
          goto 60
        endif
 2001   format(/,'Warning, the Core radius =',f5.2,
     1   /,' is larger then wave function',
     1   ' extrema position =',f5.2,/)
        if (rc(lp) .lt. zero) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c  Reset the n quantum numbers.
c
 60     do 70 j=1,norb
          nops(j) = 0
 70     continue
        nops(i) = lp
c
c  njtj  ***  modification start  ***
c    Sset up the functions f(r/rc) and g(r/rc) and
c  modify the ionic potential.
c
        aa = (7*one)/2
        dcl = -6*one*lp
        cl = dcl
c
        do 80 j=1,nr
          rrc = r(j)/rc(lp)
          rra = rrc**aa
          f(j) = zero
          if (rra .lt. 88*one) f(j)=exp(-rra)
          g(j) = rrc**lp * f(j)
          fjm1 = one-f(j)
          if (fjm1 .lt. small4) fjm1=(one-pfive*rra)*rra
          if (so(i) .lt. 0.1) then
            viod(lp,j)=fjm1*viod(lp,j)-f(j)*r(j)*vid(j)+dcl*r(j)*f(j)
          else
c
c bug fix Alberto Garcia 5/11/90
c
            viou(lp,j)=fjm1*viou(lp,j)-f(j)*r(j)*viu(j)+dcl*r(j)*f(j)
          endif
          if (rrc .lt. 3*one) j3rc = j
 80     continue
        dcl=dcl/2
c
c   Start the iteration loop to find cl.
c
        eviae = ev(i)
        devold = zero
        do 130 j=1,100
          call dsolv2(j,2,blank,ifcore,lmax,
     1     nr,a,b,r,rab,norb,ncore,nops,lo,so,zo,znuc,cdd,cdu,cdc,
     2     viod,viou,vid,viu,ev,ek,ep,wk1,wk2,wk3,wk4,wk5,wk6,
     3     wk7,evi)
          dev = eviae-ev(i)
c
c    The abs(dev-devold) condition was added to eliminate
c   division by zero errors in the calculation of
c   dcl = -dev*dcl / (dev-devold).
c
          if ((abs(dev) .lt. small2 .or. abs(dev-devold)
     1     .lt. small3) .and. j .ne. 1) then
            goto 140
          else
            if (j  .gt. 20 .or. abs(dev) .lt. 0.001) then
c
c   Use newton raphson iteration to change cl.
c
              dcl = -dev*dcl / (dev-devold)
            else
              if (dev*dcl .lt. zero) then
                dcl=-dcl/3
              endif
            endif
          endif
c
c  njtj  ***  modification end  ***
c
c  Find the new potential.
c
 100      if (so(i) .lt. 0.1) then
            do 110 k=2,nr
              viod(lp,k) = viod(lp,k) + dcl*r(k)*f(k)
 110        continue
          else
            do 111 k=2,nr
              viou(lp,k) = viou(lp,k) + dcl*r(k)*f(k)
 111        continue
          endif
          cl = cl + dcl
          devold = dev
 130    continue
c
c  End the iteration loop for cl.
c
        call ext(820+lp)
c
c   Find the pseudo-wavefunction.
c
 140    if (so(i) .lt. 0.1) then
          do 150 j=2,nr
            v(j) = (viod(lp,j)+llp/r(j))/r(j) + vid(j)
 150      continue
        else
          do 151 j=2,nr
            v(j) = (viou(lp,j)+llp/r(j))/r(j) + viu(j)
 151      continue
        endif
        call difnrl(0,i,v,arps,br,lmax,nr,a,b,r,rab,norb,
     1   nops,lo,so,znuc,viod,viou,vid,viu,ev,iflag,wk1,
     2   wk2,wk3,evi)
c
c  Compute delta and gamma.
c
        gamma = abs(ar(j3rc)/arps(j3rc)+ar(j3rc+1)/arps(j3rc+1))/2
        ag = zero
        gg = zero
        ll = 4
        do 160 j=2,nr
          ag = ag + ll*arps(j)*g(j)*rab(j)
          gg = gg + ll*g(j)*g(j)*rab(j)
          ll = 6 - ll
 160    continue
        ag = ag/3
        gg = gg/3
        delta = sqrt((ag/gg)**2+(1/gamma**2-1)/gg) - ag/gg
c
c     Modify the pseudo-wavefunction and pseudo-potential and
c   add to charge density.
c
        if (so(i) .lt. 0.1) then
          do 170 j=2,nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vod(j)=vod(j)+zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j)=small
            rrp = r(j)/rc(lp)
            gpp=(llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)
     1       *g(j)/r(j)**2
            viod(lp,j) = viod(lp,j)+gamma*delta*((ev(i)-
     1       v(j))*g(j)+gpp)*r(j)/arps(j)
 170      continue
        else
          do 171 j=2,nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vou(j)=vou(j)+zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j)=small
            rrp = r(j)/rc(lp)
            gpp=(llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)
     1       *g(j)/r(j)**2
            viou(lp,j) = viou(lp,j)+gamma*delta*((ev(i)-
     1       v(j))*g(j)+gpp)*r(j)/arps(j)
 171      continue
        endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
        ist=1
        if (arps(nr-85) .lt. zero) ist=-1
        call potrw(arps,r,nr-85,lo(i),0,ist)
        if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
        call wtrans(arps,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c
        write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cl,gamma,delta
 180    format(1x,i1,a1,f6.1,5f12.6)
 190  continue
c
c  End loop over valence orbitals.
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 200 i=ncp,norb
        nops(i) = lo(i) + 1
        zval = zval + zo(i)
 200  continue
      zion = zval+znuc-zel
      if (zval .ne. zero) zratio=zion/zval
      do 210 i=1,nr
        cdd(i) = vod(i)
 210  continue
      do 211 i=1,nr
        cdu(i) = vou(i)
 211  continue
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used
c
      if (ifcore .ne. 0) then
        ac = zero
        bc = zero
        icore = 1
        if (cfac .le. zero .or. zratio .eq. zero) then
          write(6,280) r(icore),ac,bc
        else
          if (rcfac .le. zero) then
            do 220 i=nr,2,-1
              if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220        continue
          else
            do 221 i=nr,2,-1
              if (r(i) .le. rcfac ) goto 230
 221        continue
          endif
 230      icore = i
          cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
          tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
          rbold = tpfive
          do 240 i=1,50
            rbnew = pi+atan(tanb*rbold)
            if (abs(rbnew-rbold) .lt. .00001) then
              bc = rbnew / r(icore)
              ac = cdc(icore) / (r(icore)*sin(rbnew))
              do 260 j=1,icore
                cdc(j) = ac*r(j)*sin(bc*r(j))
 260          continue
              write(6,280) r(icore),ac,bc
              goto 290
            else
              rbold=rbnew
            endif
 240      continue
          write(6,1030)
          call ext(830)
        endif
      endif
 280  format(//,' core correction used',/,
     1 ' pseudo core inside r =',f6.3,/,' ac =',f6.3,' bc =',f6.3,/)
 1030 format(//,' error in pseudb - noncovergence in finding ',
     1 /,'pseudo-core values')
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
 290  if (ispp .eq. 's') then
        blank='s'
      else
        blank=' '
      endif
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      ecut=ecuts
      do 315 i=ncp,norb
        lp = lo(i)+1
        if (so(i) .lt. 0.1) then
          do 300 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 300      continue
          rcut(i-ncore) = r(jcut)
          do 310 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 310      continue
          do 311 j=2,nr
            v(j) = viod(lp,j)/r(j)
 311      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        else
          do 312 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 312      continue
          rcut(i-ncore) = r(jcut)
          do 313 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 313      continue
          do 314 j=2,nr
            v(j) = viou(lp,j)/r(j)
 314      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        endif
 315  continue
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a
c   marker 'zio' for latter plotting
c
      write(3,4559)
      write(3,4560) zion
 4559 format(1x,'marker zio')
 4560 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (ispp .eq. 's' ) then
        do 500 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 505 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 505        continue
          else
            do 506 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 506        continue
          endif
 500    continue
      endif
c
      do 320 i=2,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 320  continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,blank,ifcore,lmax,nr,a,b,r,rab,
     1 norb-ncore,0,nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 znuc,cdd,cdu,cdc,viod,viou,vid,viu,ev(ncp),ek(ncp),
     3 ep(ncp),wk1,wk2,wk3,wk4,wk5,wk6,wk7,evi(ncp))
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,326) (ev(i),i=ncp,norb)
 325  format(//,' test of eigenvalues',//,' rcut =',8(2x,a1,f7.2))
 326  format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,330)
 330  format(///,' l    vps(0)    vpsmin      at r',/)
      do 370 i=1,lmax
        if (indd(i)+indu(i) .eq. 0) goto 370
        if (indd(i) .ne. 0) then
          vpsdm = zero
          do 350 j=2,nr
            if (r(j) .lt. .00001) goto 350
            vps = viod(i,j)/r(j)
            if (vps .lt. vpsdm) then
              vpsdm = vps
              rmind = r(j)
            endif
 350      continue
          write(6,360) il(i),viod(i,2)/r(2),vpsdm,rmind
        endif
        if (indu(i) .ne. 0) then
          vpsum = zero
          do 351 j=2,nr
            if (r(j) .lt. .00001) goto 351
            vps = viou(i,j)/r(j)
            if (vps .lt. vpsum) then
              vpsum = vps
              rminu = r(j)
            endif
 351      continue
          write(6,360) il(i),viou(i,2)/r(2),vpsum,rminu
        endif
 360  format(1x,a1,3f10.3)
 370  continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,
     1 nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1) = 'atom-lda  '
      call zedate(iray(2))
      iray(3) = 'Bachelet, '
      iray(4) = 'Hamann,and'
      iray(5) = ' Schluter '
      iray(6) = ' potential'
c
c  Encode the title array.
c
      do 390 i=1,7
        ititle(i) = '          '
 390  continue
      do 420 i=1,lmax
        if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
        zelu = zero
        zeld = zero
        if (indd(i) .ne. 0) then
          noi = no(indd(i))
          zeld = zo(indd(i))
        endif
        if (indu(i) .ne. 0) then
          noi = no(indu(i))
          zelu = zo(indu(i))
        endif
        zelt = zeld + zelu
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt
         write(ititle(2*i),401)ispp,rc(i)
 400     format(i1,a1,'(',f6.2,')')
 401     format(a1,' rc=',f5.2)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld
         write(ititle(2*i),411)zelu,ispp,rc(i)
 410     format(i1,a1,'  (',f4.2,',')
 411     format(f4.2,')',a1,f4.2)
        endif
 420  continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
        if (indu(1) .eq. 0) goto 429
        indd(1)=indu(1)
        indu(1)=0
        do 428 j=2,nr
          viod(1,j) = viou(1,j)
          viou(1,j) = zero
 428    continue
 429    do 431 i=2,lmax
          if (indd(i) .eq. 0 .or. indu(i) .eq. 0) goto 431
          do 430 j=2,nr
            viodj = viod(i,j)
            viouj = viou(i,j)
            viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
            viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430      continue
 431    continue
      endif
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 450 i=1,lmax
        if (indd(i) .ne. 0) npotd=npotd+1
        if (indu(i) .ne. 0) npotu=npotu+1
 450  continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
        if (ifull .eq. 0) then
          nicore = 'pcec'
        else
          nicore = 'fcec'
        endif
      elseif (ifcore .eq. 2) then
        if (ifull .eq. 0) then
          nicore = 'pche'
        else
          nicore = 'fche'
        endif
      else
        nicore = 'nc  '
      endif
      if (ispp .eq. 's') then
        irel='isp'
      elseif (ispp .eq. 'r') then
        irel='rel'
      else
        irel = 'nrl'
      endif
      rewind 1
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat
c  file (unit=1).
c
      do 460 i=1,lmax
        if (indd(i) .eq. 0) goto 460
        write(1) i-1,(viod(i,j),j=2,nr)
 460  continue
      do 465 i=1,lmax
        if (indu(i) .eq. 0) goto 465
        write(1) i-1,(viou(i,j),j=2,nr)
 465  continue
c
c  Write the charge densities to the current pseudo.dat
c  file (unit=1).
c
      if (ifcore .eq. 0) then
        write(1) (zero,i=2,nr)
      else
        write(1) (cdc(i),i=2,nr)
      endif
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
      end
C
C
C
      subroutine pseudk(itype,icorr,ispp,lmax,nr,a,b,
     1 r,rab,nameat,norb,ncore,no,lo,so,zo,znuc,zel,
     2 cdd,cdu,cdc,viod,viou,vid,viu,vod,vou,etot,ev,
     3 ek,ep,wk1,wk2,wk3,wk4,wk5,wk6,wk7,nops,v,ar,br,
     4 wkb,evi)
c
c *************************************************************
c *                                                           *
c *     pseudk generates the pseudo potential using the       *
c *   scheme of G. P. Kerker, J. Phys. C13, L189 (1980).      *
c *                                                           *
c *************************************************************
c
c  njtj  ***  modifications  ***
c    The only major modification is in the spin-polarization
c    treatment of the unscreening of the pseudopotential.
c    Spin-polarized potentails should be unscreend with
c    spin-polarized valence charge.  This was not
c    done in pseudo and pseudk in earlier Berkeley/Froyen
c    versions of this program.
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,tpfive=2.5D0,pfive=0.5D0,smtol=1.D-12)
      parameter (one=1.D0,ai=2*137.0360411D0,ecuts=1.0D-3)
Cray       parameter (zero=0.0,tpfive=2.5,pfive=0.5,smtol=1.E-12)
Cray       parameter (one=1.0,ai=2*137.0360411,ecuts=1.0D-3)
c
      character*1 ispp,blank,il(5)
      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 iray(6),ititle(7)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),vod(nr),vou(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),wk5(nr),wk6(nr),wk7(nr),
     4 wkb(3*nr),nops(norb),v(nr),ar(nr),br(nr),evi(norb)
c
      dimension indd(5),indu(5),rc(5),rcut(10),etot(10)
c
      data il/'s','p','d','f','g'/
      do 3 i=1,5
        indd(i)=0
        indu(i)=0
 3    continue
      if (ncore .eq. norb) return
      if (itype .lt. 1 .or. itype .gt. 3) return
      ifcore = itype-1
      pi = 4*atan(one)
c
c  read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be
c  set to greater then zero.
c
      read(5,10) (rc(i),i=1,5),cfac,rcfac
 10   format(7f10.5)
      if (cfac .eq. zero) cfac=one
c
c    Reset vod and vou to zero, they are used to store the pseudo
c  valence charge density.
c
      do 15 i=1,nr
        vod(i) = zero
        vou(i) = zero
 15   continue
c
c  Print the heading.
c
      write(6,20) nameat
 20   format(//,a2,' Pseudopotential generation using the method',
     1 ' of Kerker',/,1x,60('-'),//,
     2 ' nl    s    eigenvalue',6x,'rc',4x,6x,'cdrc',7x,'delta',
     3 7x,/)
c
c      start loop over valence orbitals
c
      ncp = ncore+1
      do 190 i=ncp,norb
        lp = lo(i) + 1
        llp = lo(i)*lp
        if (so(i) .lt. 0.1 .and. indd(lp) .ne. 0) call ext(800+lp)
        if (so(i) .gt. 0.1 .and. indu(lp) .ne. 0) call ext(810+lp)
        if (so(i) .lt. 0.1) then
          indd(lp) = i
        else
          indu(lp) = i
        endif
c
c  Find the all-electron wave function.
c
        if (so(i) .lt. 0.1) then
          do 30 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 30       continue
        else
          do 31 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 31       continue
        endif
        if (ispp .ne. 'r') then
          do 32 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 32       continue
        endif
        if (ispp .ne. 'r') then
          call difnrl(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,wk1,wk2,wk3,wk4,evi)
        endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c
        ist=1
        if (ar(nr-85) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),1,ist)
c
c  njtj  ***  user should adjust for their needs  ***
c
c  Find the last zero and extremum point.
c
        ka = lo(i)+1
        if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
        nextr = no(i)-lo(i)
        rzero = zero
        arp = br(2)
c
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp=ka*ar(2)/r(2)+(ev(i)-viod(lp,2)/r(2)-vid(2)+ai*ai)*br(2)/ai
          else
            arp=ka*ar(2)/r(2)+(ev(i)-viou(lp,2)/r(2)-viu(2)+ai*ai)*br(2)/ai
          endif
        endif
c
        do 40 j=3,nr-7
          if (nextr .eq. 0) goto 50
          if (ar(j-1)*ar(j) .le. zero)
     1     rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
          arpm = arp
          arp = br(j)
c
          if (ispp .eq. 'r') then
            if(so(i) .lt. 0.1) then
              arp=ka*ar(j)/r(j)+(ev(i)-viod(lp,j)/r(j)-
     1         vid(j)+ai*ai)*br(j)/ai
            else
              arp=ka*ar(j)/r(j)+(ev(i)-viou(lp,j)/r(j)-
     1         viu(j)+ai*ai)*br(j)/ai
            endif
          endif
c
          if (arp*arpm .gt. zero) goto 40
          rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
          nextr = nextr - 1
 40     continue
c
c   Check rc, if outside bounds reset.
c
 50     if (rzero .lt. r(2)) rzero = r(2)
        if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) goto 60
        if (rc(lp) .lt. zero) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c   Find index for grid point closest to rc.
c
 60     do 70 j=1,nr
          if (r(j) .gt. rc(lp)) goto 80
          jrc = j
 70     continue
c
c   Reset the n quantum numbers.
c
 80     rc(lp)=r(jrc)
        do 90 j=1,norb
          nops(j) = 0
 90     continue
        nops(i) = lp
c
c  Find the integrated charge inside rc.
c
        ll = 2
        if (ispp .eq. 'r') then
          cdrc = -(ar(jrc)*ar(jrc)+br(jrc)*br(jrc))*rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 102 k=jrc,1,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 102        continue
          else
            do 103 k=jrc,4,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 103        continue
            cdrc = cdrc-(ar(4)*ar(4)+br(4)*br(4))*rab(4)
            cdrc = cdrc+9*((ar(1)*ar(1)+br(1)*br(1))*rab(1)+
     1       3*(ar(2)*ar(2)+br(2)*br(2))*rab(2)+
     2       3*(ar(3)*ar(3)+br(3)*br(3))*rab(3)+
     3       (ar(4)*ar(4)+br(4)*br(4))*rab(4))/8
          endif
          cdrc = cdrc/3
        else
          cdrc = - ar(jrc) * ar(jrc) * rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 100 k=jrc,1,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 100        continue
          else
            do 101 k=jrc,4,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 101        continue
            cdrc = cdrc - ar(4) * ar(4) * rab(4)
            cdrc = cdrc + 9 * ( ar(1) * ar(1) * rab(1) +
     1       3 * ar(2) *ar(2) * rab(2) +
     2       3 * ar(3) *ar(3) * rab(3) +
     3       ar(4) * ar(4) * rab(4))/8
          endif
          cdrc = cdrc/3
        endif
c
c   The initial values for alpha, beta, gamma and delta.
c
         rc2 = r(jrc) * r(jrc)
         rc3 = r(jrc) * rc2
         rc4 = r(jrc) * rc3
         iswtch = 1
         if (ar(jrc) .lt. zero) iswtch = -1
         arc = iswtch * ar(jrc)
         arp = br(jrc)
c
         if (ispp .eq. 'r') then
           if(so(i) .lt. 0.1) then
             arp=ka*ar(jrc)/r(jrc)+(ev(i)-viod(lp,jrc)/r(jrc)-
     1        vid(jrc) + ai*ai) * br(jrc)/ai
           else
             arp=ka*ar(jrc)/r(jrc)+(ev(i)-viou(lp,jrc)/r(jrc)-
     1        viu(jrc) + ai*ai) * br(jrc)/ai
           endif
         endif
c
         brc = arp / ar(jrc)
         if (so(i) .lt. 0.1) then
           vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
         else
           vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
         endif
         alpha = ( 3*log(arc/r(jrc)**lp) - 2*(r(jrc)*brc-lp)
     1    + (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))/2 ) / rc4
         beta  = (-8*log(arc/r(jrc)**lp) + 5*(r(jrc)*brc-lp)
     1    - (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))   ) / rc3
         gamma = ( 6*log(arc/r(jrc)**lp) - 3*(r(jrc)*brc-lp)
     1    + (rc2*vrc+lp*lp-rc2*(ev(i)+brc*brc))/2 ) / rc2
         delta = zero
c
c  Start the iteration loop to find delta.
c
         do 150 j=1,50
c
c  Generate the pseudo-wavefunction (note missing factor exp(delta)).
c
           do 110 k=1,jrc
             polyr=r(k)*r(k)*((alpha*r(k)+beta)*r(k)+gamma)
             ar(k) = iswtch * r(k)**lp * exp(polyr)
 110       continue
c
c  Integrate  the pseudo charge density from r = 0 to rc.
c
           ll = 2
           cdps = - ar(jrc) * ar(jrc) * rab(jrc)
           if (jrc .ne. 2*(jrc/2)) then
             do 120 k=jrc,1,-1
               cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
               ll = 6 - ll
 120         continue
           else
             do 121 k=jrc,4,-1
               cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
               ll = 6 - ll
 121         continue
             cdps = cdps - ar(4) * ar(4) * rab(4)
             cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1        3 * ar(2) *ar(2) * rab(2) +
     2        3 * ar(3) *ar(3) * rab(3) +
     3        ar(4) * ar(4) * rab(4))/8
           endif
           cdps = cdps/3
c
c  Find the new delta.
c
           fdnew = log(cdrc/cdps) - 2*delta
           if (abs(fdnew) .lt. smtol) goto 160
           if (j .eq. 1) then
             ddelta = pfive
           else
             ddelta = - fdnew * ddelta / (fdnew-fdold)
           endif
           alpha = alpha - 3*ddelta/rc4
           beta  = beta  + 8*ddelta/rc3
           gamma = gamma - 6*ddelta/rc2
           delta = delta + ddelta
           fdold = fdnew
 150     continue
c
c  End the iteration loop for delta.
c
         call ext(820+lp)
c
c    Augment the charge density and invert schroedinger equation
c  to find new potential.
c
 160     expd = exp(delta)
         if (so(i) .lt. 0.1) then
           do 170 j=1,jrc
             ar(j) = expd * ar(j)
             xlamda=(4*alpha*r(j)+3*beta)*r(j)+2*gamma
             vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1        + (12 * alpha * r(j) + 6 * beta) * r(j) + 2 * gamma
             viod(lp,j) = (vj - vid(j)) * r(j)
             vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
 170       continue
           do 171 j=jrc+1,nr
             vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
 171       continue
         else
           do 175 j=1,jrc
             ar(j) = expd * ar(j)
             xlamda=(4*alpha*r(j)+3*beta)*r(j)+2*gamma
             vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1        + (12 * alpha * r(j) + 6 * beta) * r(j) + 2 * gamma
             viou(lp,j) = (vj - viu(j)) * r(j)
             vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
 175       continue
           do 176 j=jrc+1,nr
             vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
 176       continue
         endif
         write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cdrc,delta
 180     format(1x,i1,a1,f6.1,5f12.6)
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
         ist=1
         if (ar(nr-85) .lt. zero) ist=-1
         call potrw(ar,r,nr-85,lo(i),0,ist)
         if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
         call wtrans(ar,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c

 190   continue
c
c  End loop over valence orbitals.
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 200 i=ncp,norb
        nops(i) = lo(i) + 1
        zval = zval + zo(i)
 200  continue
      zion = zval+znuc-zel
      if (zval .ne. zero) zratio=zion/zval
      do 210 i=1,nr
        cdd(i) = vod(i)
 210  continue
      do 211 i=1,nr
        cdu(i) = vou(i)
 211  continue
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used
c
      if (ifcore .ne. 0) then
        ac = zero
        bc = zero
        icore = 1
        if (cfac .le. zero .or. zratio .eq. zero) then
          write(6,280) r(icore),ac,bc
        else
          if (rcfac .le. zero) then
            do 220 i=nr,2,-1
              if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220        continue
          else
            do 221 i=nr,2,-1
              if (r(i) .le. rcfac ) goto 230
 221        continue
          endif
 230      icore = i
          cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
          tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
          rbold = tpfive
          do 240 i=1,50
            rbnew = pi+atan(tanb*rbold)
            if (abs(rbnew-rbold) .lt. .00001) then
              bc = rbnew / r(icore)
              ac = cdc(icore) / (r(icore)*sin(rbnew))
              do 260 j=1,icore
                cdc(j) = ac*r(j)*sin(bc*r(j))
 260          continue
              write(6,280) r(icore),ac,bc
              goto 290
            else
              rbold=rbnew
            endif
 240      continue
          write(6,1030)
          call ext(830)
        endif
      endif
 280  format(//,' core correction used',/,
     1 ' pseudo core inside r =',f6.3,/,' ac =',f6.3,' bc =',f6.3,/)
 1030 format(//,' error in pseudk - noncovergence in finding ',
     1 /,'pseudo-core values')
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
 290  if (ispp .eq. 's') then
        blank='s'
      else
        blank=' '
      endif
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      ecut=ecuts
      do 315 i=ncp,norb
        lp = lo(i)+1
        if (so(i) .lt. 0.1) then
          do 500 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 500      continue
          rcut(i-ncore) = r(jcut)
          do 510 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 510      continue
          do 511 j=2,nr
            v(j) = viod(lp,j)/r(j)
 511      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        else
          do 512 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 512      continue
          rcut(i-ncore) = r(jcut)
          do 513 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 513      continue
          do 514 j=2,nr
            v(j) = viou(lp,j)/r(j)
 514      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        endif
 315  continue
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a
c   marker 'zio' for latter plotting
c
      write(3,4559)
      write(3,4560) zion
 4559 format(1x,'marker zio')
 4560 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (ispp .eq. 's' ) then
        do 700 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 705 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 705        continue
          else
            do 706 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 706        continue
          endif
 700    continue
      endif
c
      do 320 i=2,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 320  continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,blank,ifcore,lmax,nr,a,b,r,rab,
     1 norb-ncore,0,nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 znuc,cdd,cdu,cdc,viod,viou,vid,viu,ev(ncp),ek(ncp),
     3 ep(ncp),wk1,wk2,wk3,wk4,wk5,wk6,wk7,evi(ncp))
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,326) (ev(i),i=ncp,norb)
 325  format(//,' test of eigenvalues',//,' rcut =',8(2x,a1,f7.2))
 326  format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,330)
 330  format(///,' l    vps(0)    vpsmin      at r',/)
      do 370 i=1,lmax
        if (indd(i)+indu(i) .eq. 0) goto 370
        if (indd(i) .ne. 0) then
          vpsdm = zero
          do 350 j=2,nr
            if (r(j) .lt. .00001) goto 350
            vps = viod(i,j)/r(j)
            if (vps .lt. vpsdm) then
              vpsdm = vps
              rmind = r(j)
            endif
 350      continue
          write(6,360) il(i),viod(i,2)/r(2),vpsdm,rmind
        endif
        if (indu(i) .ne. 0) then
          vpsum = zero
          do 351 j=2,nr
            if (r(j) .lt. .00001) goto 351
            vps = viou(i,j)/r(j)
            if (vps .lt. vpsum) then
              vpsum = vps
              rminu = r(j)
            endif
 351      continue
          write(6,360) il(i),viou(i,2)/r(2),vpsum,rminu
        endif
 360  format(1x,a1,3f10.3)
 370  continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,
     1 nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1) = 'atom-lda  '
      call zedate(iray(2))
      iray(3) = '   Kerker-'
      iray(4) = 'potential '
      do 380 i=5,6
        iray(i) = '          '
 380  continue
c
c  Encode the title array.
c
      do 390 i=1,7
        ititle(i) = '          '
 390  continue
      do 420 i=1,lmax
        if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
        zelu = zero
        zeld = zero
        if (indd(i) .ne. 0) then
          noi = no(indd(i))
          zeld = zo(indd(i))
        endif
        if (indu(i) .ne. 0) then
          noi = no(indu(i))
          zelu = zo(indu(i))
        endif
        zelt = zeld + zelu
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt
         write(ititle(2*i),401)ispp,rc(i)
 400     format(i1,a1,'(',f6.2,')')
 401     format(a1,' rc=',f5.2)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld
         write(ititle(2*i),411)zelu,ispp,rc(i)
 410     format(i1,a1,'  (',f4.2,',')
 411     format(f4.2,')',a1,f4.2)
        endif
 420  continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
        if (indu(1) .eq. 0) goto 429
        indd(1)=indu(1)
        indu(1)=0
        do 428 j=2,nr
          viod(1,j) = viou(1,j)
          viou(1,j) = zero
 428    continue
 429    do 431 i=2,lmax
          if (indd(i) .eq. 0 .or. indu(i) .eq. 0) goto 431
          do 430 j=2,nr
            viodj = viod(i,j)
            viouj = viou(i,j)
            viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
            viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430      continue
 431    continue
      endif
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 450 i=1,lmax
        if (indd(i) .ne. 0) npotd=npotd+1
        if (indu(i) .ne. 0) npotu=npotu+1
 450  continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
        if (ifull .eq. 0) then
          nicore = 'pcec'
        else
          nicore = 'fcec'
        endif
      elseif (ifcore .eq. 2) then
        if (ifull .eq. 0) then
          nicore = 'pche'
        else
          nicore = 'fche'
        endif
      else
        nicore = 'nc  '
      endif
      if (ispp .eq. 's') then
        irel='isp'
      elseif (ispp .eq. 'r') then
        irel='rel'
      else
        irel = 'nrl'
      endif
      rewind 1
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat
c  file (unit=1).
c
      do 460 i=1,lmax
        if (indd(i) .eq. 0) goto 460
        write(1) i-1,(viod(i,j),j=2,nr)
 460  continue
      do 465 i=1,lmax
        if (indu(i) .eq. 0) goto 465
        write(1) i-1,(viou(i,j),j=2,nr)
 465  continue
c
c  Write the charge densities to the current pseudo.dat
c  file (unit=1).
c
      if (ifcore .eq. 0) then
        write(1) (zero,i=2,nr)
      else
        write(1) (cdc(i),i=2,nr)
      endif
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
      end
C
C
C
      subroutine pseudo(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2 viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,wk3,
     3 wk4,wk5,wk6,wk7,f,g,nops,v,ar,br,arps,wkb,evi)
c
c *************************************************************
c *                                                           *
c *    pseudo generates the pseudo potential using            *
c *  the scheme of Hamann, Schluter and Chiang -              *
c *  Phys. Rev. Lett. 43, 1494 (1979).                        *
c *                                                           *
c *************************************************************
c
c  njtj  *** modifications  ***
c    The only major modifications are in the spin-polarized
c    treatment of the el-el unscreening of the pseudopotential
c    A spin-polarized pseudopotential is unscreened
c    with a spin-polarized valence charge.  This was not done
c    in pseudo or pseudok in earlier versions of this
c    program.
c  njtj  *** modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter(zero=0.D0,ecuts=1.0D-3,tpfive=2.5D0,one=1.D0)
      parameter(small=1.D-13,small2=1.D-10,small3=1.D-18,pzfive=.05D0)
      parameter(pfive=0.5D0,small4=1.D-6,ai=2*137.0360411D0)
Cray       parameter(zero=0.0,ecuts=1.0E-3,tpfive=2.5,one=1.0)
Cray       parameter(small=1.E-13,small2=1.E-10,small3=1.E-18,pzfive=.05)
Cray       parameter(pfive=0.5,small4=1.E-6,ai=2*137.0360411)
c
      character*1 ispp,blank,il(5)
      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 ititle(7),iray(6)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),vod(nr),vou(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),wk5(nr),wk6(nr),wk7(nr),
     4 wkb(6*nr),f(nr),g(nr),nops(norb),v(nr),
     5 ar(nr),br(nr),arps(nr),evi(norb)
c
      dimension etot(10),indd(5),indu(5),rc(5),rcut(10)
c
      data il/'s','p','d','f','g'/
      do 3 i=1,5
        indd(i)=0
        indu(i)=0
 3    continue
      if (ncore .eq. norb) return
      if (itype .ne. 1 .and. itype .ne. 2 .and. itype .ne. 3) return
      ifcore = itype - 1
      pi = 4*atan(one)
c
c  Spin-polarized potentails should be unscreened with
c  a spin-polarized valence charge.  This was not
c  done in pseudo and pseudk in earlier versions
c  of this program.
c
      if (ispp .eq. 's' ) then
        blank = 's'
      else
        blank = ' '
      endif
c
c  read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be
c  set to greater then zero.
c
      read(5,10) (rc(i),i=1,5),cfac,rcfac
 10   format(7f10.5)
      if (cfac .eq. zero) cfac=one
c
c   Reset vod and vou to zero.  They are here used to store
c   the pseudo valence charge density.
c
      do 15 i=1,nr
        vod(i) = zero
        vou(i) = zero
 15   continue
c
c  Print the heading.
c
      write(6,20) nameat
 20   format(//,a2,' Pseudopotential HSC generation',/,1x,35('-'),//,
     1 ' nl    s    eigenvalue',6x,'rc',4x,6x,'cl',9x,'gamma',
     2 7x,'delta',/)
c
c      start loop over valence orbitals
c
      ncp = ncore+1
      do 190 i=ncp,norb
        lp = lo(i) + 1
        llp = lo(i)*lp
        if (so(i) .lt. 0.1) then
          if (indd(lp) .ne. 0) then
            write(6,1000)lp-1
            call ext(800+lp)
          else
            indd(lp) = i
          endif
        else
          if (indu(lp) .ne. 0) then
            write(6,1010)lp-1
            call ext(810+lp)
          else
            indu(lp) = i
          endif
        endif
 1000 format(//,'error in pseudo - two down spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
 1010 format(//,'error in pseudo - two up spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
c
c      find all electron wave function
c
        do 25 j=1,nr
          ar(j)=zero
 25     continue
        if (so(i) .lt. 0.1) then
          do 27 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 27       continue
        else
          do 30 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 30       continue
        endif
        if (ispp .ne. 'r') then
          do 32 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 32       continue
        endif
        if (ispp .ne. 'r') then
          call difnrl(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1     znuc,viod,viou,vid,viu,ev,wk1,wk2,wk3,wk4,evi)
        endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c
        ist=1
        if (ar(nr-85) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),1,ist)
c
c  njtj  ***  user should adjust for their needs  ***
c
c  Find the last zero and extremum.
c
        ka = lo(i)+1
        if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
        nextr = no(i)-lo(i)
        rzero = zero
        arp = br(2)
c
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp = ka*ar(2)/r(2) + (ev(i) - viod(lp,2)/r(2)
     1       - vid(2) + ai*ai) * br(2) / ai
          else
            arp = ka*ar(2)/r(2) + (ev(i) - viou(lp,2)/r(2)
     1       - viu(2) + ai*ai) * br(2) / ai
          endif
        endif
c
        do 40 j=3,nr-7
          if (nextr .eq. 0) goto 50
          if (ar(j-1)*ar(j) .le. zero)
     1     rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
          arpm = arp
          arp = br(j)
c
          if (ispp .eq. 'r') then
            if (so(i) .lt. 0.1) then
              arp = ka*ar(j)/r(j) + (ev(i) - viod(lp,j)/r(j)
     1         - vid(j) + ai*ai) * br(j) / ai
            else
              arp = ka*ar(j)/r(j) + (ev(i) - viou(lp,j)/r(j)
     1         - viu(j) + ai*ai) * br(j) / ai
            endif
          endif
c
          if (arp*arpm .gt. zero) goto 40
          rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
          nextr = nextr - 1
 40     continue
c
c  Check rc, if outside bounds reset.
c
 50     if (rzero .lt. r(2)) rzero = r(2)
        if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) goto 60
        if (rc(lp) .ge. rzero) then
          write(6,2001)rc(lp),rextr
        endif
 2001   format(/,'Warning, the Core radius =',f5.2,
     1   /,' is larger then wave function',
     1   ' extrema position =',f5.2,/)
        if (rc(lp) .lt. zero) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c  Reset the n quantum numbers.
c
 60     do 70 j=1,norb
          nops(j) = 0
 70     continue
        nops(i) = lp
c
c  njtj  ***  modification start  ***
c    Sset up the functions f(r/rc) and g(r/rc) and
c  modify the ionic potential.
c
        aa = 4*one
        dcl = -6*one*lp
        cl = dcl
c
        do 80 j=1,nr
          rrc = r(j)/rc(lp)
          rra = rrc**aa
          f(j) = zero
          if (rra .lt. 88*one) f(j)=exp(-rra)
          g(j) = rrc**lp * f(j)
          fjm1 = one-f(j)
          if (fjm1 .lt. small4) fjm1=(one-pfive*rra)*rra
          if (so(i) .lt. 0.1) then
            viod(lp,j)=fjm1*viod(lp,j)-f(j)*r(j)*vid(j)+dcl*r(j)*f(j)
          else
c
c bug fix Alberto Garcia 5/11/90
c
            viou(lp,j)=fjm1*viou(lp,j)-f(j)*r(j)*viu(j)+dcl*r(j)*f(j)
          endif
          if (rrc .lt. 3*one) j3rc = j
 80     continue
        dcl=dcl/2
c
c   Start the iteration loop to find cl.
c
        eviae = ev(i)
        devold = zero
        do 130 j=1,100
          call dsolv2(j,2,blank,ifcore,lmax,
     1     nr,a,b,r,rab,norb,ncore,nops,lo,so,zo,znuc,cdd,cdu,cdc,
     2     viod,viou,vid,viu,ev,ek,ep,wk1,wk2,wk3,wk4,wk5,wk6,
     3     wk7,evi)
          dev = eviae-ev(i)
c
c    The abs(dev-devold) condition was added to eliminate
c   division by zero errors in the calculation of
c   dcl = -dev*dcl / (dev-devold).
c
          if ((abs(dev) .lt. small2 .or. abs(dev-devold)
     1     .lt. small3) .and. j .ne. 1) then
            goto 140
          else
            if (j  .gt. 20 .or. abs(dev) .lt. 0.001) then
c
c   Use newton raphson iteration to change cl.
c
              dcl = -dev*dcl / (dev-devold)
            else
              if (dev*dcl .lt. zero) then
                dcl=-dcl/3
              endif
            endif
          endif
c
c  njtj  ***  modification end  ***
c
c  Find the new potential.
c
 100      if (so(i) .lt. 0.1) then
            do 110 k=2,nr
              viod(lp,k) = viod(lp,k) + dcl*r(k)*f(k)
 110        continue
          else
            do 111 k=2,nr
              viou(lp,k) = viou(lp,k) + dcl*r(k)*f(k)
 111        continue
          endif
          cl = cl + dcl
          devold = dev
 130    continue
c
c  End the iteration loop for cl.
c
        call ext(820+lp)
c
c   Find the pseudo-wavefunction.
c
 140    if (so(i) .lt. 0.1) then
          do 150 j=2,nr
            v(j) = (viod(lp,j)+llp/r(j))/r(j) + vid(j)
 150      continue
        else
          do 151 j=2,nr
            v(j) = (viou(lp,j)+llp/r(j))/r(j) + viu(j)
 151      continue
        endif
        call difnrl(0,i,v,arps,br,lmax,nr,a,b,r,rab,norb,
     1   nops,lo,so,znuc,viod,viou,vid,viu,ev,iflag,wk1,
     2   wk2,wk3,evi)
c
c  Compute delta and gamma.
c
        gamma = abs(ar(j3rc)/arps(j3rc)+ar(j3rc+1)/arps(j3rc+1))/2
        ag = zero
        gg = zero
        ll = 4
        do 160 j=2,nr
          ag = ag + ll*arps(j)*g(j)*rab(j)
          gg = gg + ll*g(j)*g(j)*rab(j)
          ll = 6 - ll
 160    continue
        ag = ag/3
        gg = gg/3
        delta = sqrt((ag/gg)**2+(1/gamma**2-1)/gg) - ag/gg
c
c     Modify the pseudo-wavefunction and pseudo-potential and
c   add to charge density.
c
        if (so(i) .lt. 0.1) then
          do 170 j=2,nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vod(j)=vod(j)+zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j)=small
            rrp = r(j)/rc(lp)
            gpp=(llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)
     1       *g(j)/r(j)**2
            viod(lp,j) = viod(lp,j)+gamma*delta*((ev(i)-
     1       v(j))*g(j)+gpp)*r(j)/arps(j)
 170      continue
        else
          do 171 j=2,nr
            arps(j) = gamma*(arps(j)+delta*g(j))
            vou(j)=vou(j)+zo(i)*arps(j)*arps(j)
            if (arps(j) .lt. small .and. r(j) .gt. one) arps(j)=small
            rrp = r(j)/rc(lp)
            gpp=(llp-aa*(2*lp+aa-1)*rrp**aa+(aa*rrp**aa)**2)
     1       *g(j)/r(j)**2
            viou(lp,j) = viou(lp,j)+gamma*delta*((ev(i)-
     1       v(j))*g(j)+gpp)*r(j)/arps(j)
 171      continue
        endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
        ist=1
        if (arps(nr-85) .lt. zero) ist=-1
        call potrw(arps,r,nr-85,lo(i),0,ist)
        if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
        call wtrans(arps,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c
        write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cl,gamma,delta
 180    format(1x,i1,a1,f6.1,5f12.6)
 190  continue
c
c  End loop over valence orbitals.
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 200 i=ncp,norb
        nops(i) = lo(i) + 1
        zval = zval + zo(i)
 200  continue
      zion = zval+znuc-zel
      if (zval .ne. zero) zratio=zion/zval
      do 210 i=1,nr
        cdd(i) = vod(i)
 210  continue
      do 211 i=1,nr
        cdu(i) = vou(i)
 211  continue
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used
c
      if (ifcore .ne. 0) then
        ac = zero
        bc = zero
        icore = 1
        if (cfac .le. zero .or. zratio .eq. zero) then
          write(6,280) r(icore),ac,bc
        else
          if (rcfac .le. zero) then
            do 220 i=nr,2,-1
              if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220        continue
          else
            do 221 i=nr,2,-1
              if (r(i) .le. rcfac ) goto 230
 221        continue
          endif
 230      icore = i
          cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
          tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
          rbold = tpfive
          do 240 i=1,50
            rbnew = pi+atan(tanb*rbold)
            if (abs(rbnew-rbold) .lt. .00001) then
              bc = rbnew / r(icore)
              ac = cdc(icore) / (r(icore)*sin(rbnew))
              do 260 j=1,icore
                cdc(j) = ac*r(j)*sin(bc*r(j))
 260          continue
              write(6,280) r(icore),ac,bc
              goto 290
            else
              rbold=rbnew
            endif
 240      continue
          write(6,1030)
          call ext(830)
        endif
      endif
 280  format(//,' core correction used',/,
     1 ' pseudo core inside r =',f6.3,/,' ac =',f6.3,' bc =',f6.3,/)
 1030 format(//,' error in pseudo - noncovergence in finding ',
     1 /,'pseudo-core values')
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
 290  if (ispp .eq. 's') then
        blank='s'
      else
        blank=' '
      endif
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      ecut=ecuts
      do 315 i=ncp,norb
        lp = lo(i)+1
        if (so(i) .lt. 0.1) then
          do 300 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 300      continue
          rcut(i-ncore) = r(jcut)
          do 310 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 310      continue
          do 311 j=2,nr
            v(j) = viod(lp,j)/r(j)
 311      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        else
          do 312 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 312      continue
          rcut(i-ncore) = r(jcut)
          do 313 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 313      continue
          do 314 j=2,nr
            v(j) = viou(lp,j)/r(j)
 314      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        endif
 315  continue
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a
c   marker 'zio' for latter plotting
c
      write(3,4559)
      write(3,4560) zion
 4559 format(1x,'marker zio')
 4560 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c

c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (ispp .eq. 's' ) then
        do 500 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 505 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 505        continue
          else
            do 506 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 506        continue
          endif
 500    continue
      endif
c
      do 320 i=2,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 320  continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,blank,ifcore,lmax,nr,a,b,r,rab,
     1 norb-ncore,0,nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 znuc,cdd,cdu,cdc,viod,viou,vid,viu,ev(ncp),ek(ncp),
     3 ep(ncp),wk1,wk2,wk3,wk4,wk5,wk6,wk7,evi(ncp))
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,326) (ev(i),i=ncp,norb)
 325  format(//,' test of eigenvalues',//,' rcut =',8(2x,a1,f7.2))
 326  format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,330)
 330  format(///,' l    vps(0)    vpsmin      at r',/)
      do 370 i=1,lmax
        if (indd(i)+indu(i) .eq. 0) goto 370
        if (indd(i) .ne. 0) then
          vpsdm = zero
          do 350 j=2,nr
            if (r(j) .lt. .00001) goto 350
            vps = viod(i,j)/r(j)
            if (vps .lt. vpsdm) then
              vpsdm = vps
              rmind = r(j)
            endif
 350      continue
          write(6,360) il(i),viod(i,2)/r(2),vpsdm,rmind
        endif
        if (indu(i) .ne. 0) then
          vpsum = zero
          do 351 j=2,nr
            if (r(j) .lt. .00001) goto 351
            vps = viou(i,j)/r(j)
            if (vps .lt. vpsum) then
              vpsum = vps
              rminu = r(j)
            endif
 351      continue
          write(6,360) il(i),viou(i,2)/r(2),vpsum,rminu
        endif
 360  format(1x,a1,3f10.3)
 370  continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,
     1 nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1) = 'atom-lda  '
      call zedate(iray(2))
      iray(3) = '   Hamann,'
      iray(4) = ' Schluter '
      iray(5) = 'and Chiang'
      iray(6) = ' potential'
c
c  Encode the title array.
c
      do 390 i=1,7
        ititle(i) = '          '
 390  continue
      do 420 i=1,lmax
        if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
        zelu = zero
        zeld = zero
        if (indd(i) .ne. 0) then
          noi = no(indd(i))
          zeld = zo(indd(i))
        endif
        if (indu(i) .ne. 0) then
          noi = no(indu(i))
          zelu = zo(indu(i))
        endif
        zelt = zeld + zelu
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt
         write(ititle(2*i),401)ispp,rc(i)
 400     format(i1,a1,'(',f6.2,')')
 401     format(a1,' rc=',f5.2)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld
         write(ititle(2*i),411)zelu,ispp,rc(i)
 410     format(i1,a1,'  (',f4.2,',')
 411     format(f4.2,')',a1,f4.2)
        endif
 420  continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
        if (indu(1) .eq. 0) goto 429
        indd(1)=indu(1)
        indu(1)=0
        do 428 j=2,nr
          viod(1,j) = viou(1,j)
          viou(1,j) = zero
 428    continue
 429    do 431 i=2,lmax
          if (indd(i) .eq. 0 .or. indu(i) .eq. 0) goto 431
          do 430 j=2,nr
            viodj = viod(i,j)
            viouj = viou(i,j)
            viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
            viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430      continue
 431    continue
      endif
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 450 i=1,lmax
        if (indd(i) .ne. 0) npotd=npotd+1
        if (indu(i) .ne. 0) npotu=npotu+1
 450  continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
        if (ifull .eq. 0) then
          nicore = 'pcec'
        else
          nicore = 'fcec'
        endif
      elseif (ifcore .eq. 2) then
        if (ifull .eq. 0) then
          nicore = 'pche'
        else
          nicore = 'fche'
        endif
      else
        nicore = 'nc  '
      endif
      if (ispp .eq. 's') then
        irel='isp'
      elseif (ispp .eq. 'r') then
        irel='rel'
      else
        irel = 'nrl'
      endif
      rewind 1
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat
c  file (unit=1).
c
      do 460 i=1,lmax
        if (indd(i) .eq. 0) goto 460
        write(1) i-1,(viod(i,j),j=2,nr)
 460  continue
      do 465 i=1,lmax
        if (indu(i) .eq. 0) goto 465
        write(1) i-1,(viou(i,j),j=2,nr)
 465  continue
c
c  Write the charge densities to the current pseudo.dat
c  file (unit=1).
c
      if (ifcore .eq. 0) then
        write(1) (zero,i=2,nr)
      else
        write(1) (cdc(i),i=2,nr)
      endif
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
      end
C
C
C
      subroutine pseudt(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2 viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,
     3 wk3,wk4,wk5,wk6,wk7,nops,v,ar,br,wkb,evi)
c
c *************************************************************
c *                                                           *
c *     This routine was written by Norman J. Troullier Jr.   *
c *   Sept. 1989, while at the U. of Minnesota, all           *
c *   comments concerning this routine should be directed     *
c *   to him.                                                 *
c *                                                           *
c *     troullie@128.101.224.101                              *
c *     troullie@csfsa.cs.umn.edu                             *
c *     612 625-0392                                          *
c *                                                           *
c *     pseudt generates a pseudopotential using the          *
c *   scheme of N. Troullier and J. L. Martins.               *
c *   The general format of this routine is the same as the   *
c *   pseudo and pseudk routines.  Output/input is            *
c *   compatible.                                             *
c *                                                           *
c *************************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,one=1.D0,tpfive=2.5D0,ecuts=1.0D-3)
      parameter (small=1.D-12,pnine=0.9D0,ai=2*137.0360411D0)
Cray      parameter (zero=0.0,one=1.0,tpfive=2.5,ecuts=1.0E-3)
Cray      parameter (small=1.E-12,pnine=0.9,ai=2*137.0360411)
c
      character*1 ispp,blank,il(5)
      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 iray(6),ititle(7)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),vod(nr),vou(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),wk5(nr),wk6(nr),wk7(nr),
     4 wkb(3*nr),nops(norb),v(nr),ar(nr),br(nr),evi(norb)
c
      dimension indd(5),indu(5),rc(5),rcut(10),
     1 etot(10),aa(7),rr(7),coe(7),aj(5,5),bj(5)
c
      data il/'s','p','d','f','g'/
      if (ncore .eq. norb) return
      ifcore = itype-1
      pi = 4*atan(one)
      do 3 i=1,5
        indd(i)=0
        indu(i)=0
 3    continue
      do 4 i=1,40
        nops(i) = 0
 4    continue
c
c  read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be
c  set to greater then zero.
c
      read(5,10) (rc(i),i=1,5),cfac,rcfac
 10   format(7f10.5)
      if (cfac .eq. 0.D0) cfac=one
c
c  Reset vod and vou to zero,
c  they are here used to store the pseudo valence charge density.
c
      do 15 i=1,nr
        vod(i) = zero
 15   continue
      do 16 i=1,nr
        vou(i) = zero
 16   continue
c
c  print heading
c
      write(6,20) nameat
 20   format(//,1x,a2,' pseudopotential generation using the ',
     1 'Troullier and Martins method',/,1x,60('-'),//,
     2 ' nl    s    eigenvalue',6x,'rc',10x,'cdrc',7x,'delta',/)
c
c  Start loop over valence orbitals, only one orbital for each
c  angular momentum and spin can exist.
c
      ncp = ncore+1
      do 190 i=ncp,norb
        lp = lo(i) + 1
        llp = lo(i)*lp
        if (so(i) .lt. 0.1) then
          if (indd(lp) .ne. 0) then
            write(6,1000)lp-1
            call ext(800+lp)
          else
            indd(lp) = i
          endif
        else
          if (indu(lp) .ne. 0) then
            write(6,1010)lp-1
            call ext(810+lp)
          else
            indu(lp) = i
          endif
        endif
 1000 format(//,'error in pseudt - two down spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
 1010 format(//,'error in pseudt - two up spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
c
c  Find the all electron wave function.
c
        do 29 j=1,nr
          ar(j) = zero
 29     continue
        if (so(i) .lt. 0.1) then
          do 30 j=2,nr
            v(j) = viod(lp,j)/r(j) + vid(j)
 30       continue
        else
          do 31 j=2,nr
            v(j) = viou(lp,j)/r(j) + viu(j)
 31       continue
        endif
        if (ispp .ne. 'r') then
          do 32 j=2,nr
            v(j) = v(j) + llp/r(j)**2
 32       continue
        endif
c
c  The parameter iflag has been added as a nonconvegence
c  indicator for auxillary routines.  Its value does
c  not change its operation.  iflag is a returned value,
c  set to 1 for none convergence.
c
        if (ispp .ne. 'r') then
          iflag=0
          call difnrl(0,i,v,ar,br,lmax,nr,a,b,
     1     r,rab,norb,no,lo,so,znuc,viod,viou,
     2     vid,viu,ev,iflag,wk1,wk2,wk3,evi)
        else
          call difrel(0,i,v,ar,br,lmax,nr,a,b,r,
     1     rab,norb,no,lo,so,znuc,viod,viou,vid,viu,
     2     ev,wk1,wk2,wk3,wk4,evi)
         endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c
        ist=1
        if (ar(nr-85) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),1,ist)
c
c  njtj  ***  user should adjust for their needs  ***
c
c
c  Find last zero and extremum
c
        ka = lo(i)+1
        if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
        nextr = no(i)-lo(i)
        rzero = zero
        arp = br(2)
c
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp = ka*ar(2)/r(2) + (ev(i) - viod(lp,2)/r(2)
     1       - vid(2) + ai*ai) * br(2) / ai
          else
            arp = ka*ar(2)/r(2) + (ev(i) - viou(lp,2)/r(2)
     1       - viu(2) + ai*ai) * br(2) / ai
          endif
        endif
c
        do 40 j=3,nr-7
          if (nextr .eq. 0) goto 50
          if (ar(j-1)*ar(j) .le. zero)
     1     rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
          arpm = arp
          arp = br(j)
c
          if (ispp .eq. 'r') then
            if(so(i) .lt. 0.1) then
              arp = ka*ar(j)/r(j) + (ev(i) - viod(lp,j)/r(j)
     1         - vid(j) + ai*ai) * br(j) / ai
            else
              arp = ka*ar(j)/r(j) + (ev(i) - viou(lp,j)/r(j)
     1         - viu(j) + ai*ai) * br(j) / ai
            endif
          endif
c
          if (arp*arpm .gt. zero) goto 40
          rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
          nextr = nextr - 1
 40     continue
 50     if (rzero .lt. r(2)) rzero = r(2)
c
c  Check rc if inside rzero,
c  reset to .9 between rmax and rzero if inside
c  if rc(lp) is negative, rc(lp) is percent of way
c  betweeen rzero and rmax.
c
        if (rc(lp) .gt. rzero) then
        elseif(rc(lp) .ge. zero) then
          rc(lp) = rzero + pnine*(rextr-rzero)
        else
          rc(lp) = rzero - rc(lp)*(rextr-rzero)
        endif
c
c  Find the index for odd grid point closest to rc.
c
        do 70 j=1,nr
          if (r(j) .gt. rc(lp)) goto 80
 70     continue
 80     jrc=j-1
        rc(lp)=r(jrc)
c
c  Reset n quantum numbers.
c
        nops(i) = lp
c
c  Find the integrated charge inside rc(1-charge outside).
c
        ll = 2
        if (ispp .eq. 'r') then
          cdrc = -(ar(jrc)*ar(jrc)+br(jrc)*br(jrc))*rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 102 k=jrc,1,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 102        continue
          else
            do 103 k=jrc,4,-1
              cdrc = cdrc+ll*(ar(k)*ar(k)+br(k)*br(k))*rab(k)
              ll = 6 - ll
 103        continue
            cdrc = cdrc-(ar(4)*ar(4)+br(4)*br(4))*rab(4)
            cdrc = cdrc+9*((ar(1)*ar(1)+br(1)*br(1))*rab(1)+
     1       3*(ar(2)*ar(2)+br(2)*br(2))*rab(2)+
     2       3*(ar(3)*ar(3)+br(3)*br(3))*rab(3)+
     3       (ar(4)*ar(4)+br(4)*br(4))*rab(4))/8
          endif
          cdrc = cdrc/3
        else
          cdrc = - ar(jrc) * ar(jrc) * rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 100 k=jrc,1,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 100        continue
          else
            do 101 k=jrc,4,-1
              cdrc = cdrc +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 101        continue
            cdrc = cdrc - ar(4) * ar(4) * rab(4)
            cdrc = cdrc + 9 * ( ar(1) * ar(1) * rab(1) +
     1       3 * ar(2) *ar(2) * rab(2) +
     2       3 * ar(3) *ar(3) * rab(3) +
     3       ar(4) * ar(4) * rab(4))/8
          endif
          cdrc = cdrc/3
        endif
c
c  Find the values for wave(arc), d(wave)/dr(arp), potential(vrc),
c  d(potential)/dr(vrp), and d2(potential)/dr2(vrpp)
c
        rc1 = r(jrc)
        rc2 = rc1 * rc1
        rc3 = rc2 * rc1
        rc4 = rc2 * rc2
        rc5 = rc4 * rc1
        rc6 = rc4 * rc2
        rc7 = rc4 * rc3
        rc8 = rc4 * rc4
        iswtch = 1
        if (ar(jrc) .lt. zero) iswtch = -1
        arc = iswtch * ar(jrc)
        arp = br(jrc)
        if (ispp .eq. 'r') then
          if (so(i) .lt. 0.1) then
            arp=ka*ar(jrc)/r(jrc) + (ev(i) - viod(lp,jrc)/r(jrc)
     1       - vid(jrc) + ai*ai) * br(jrc)/ai
          else
            arp=ka*ar(jrc)/r(jrc) + (ev(i) - viou(lp,jrc)/r(jrc)
     1       - viu(jrc) + ai*ai) * br(jrc)/ai
          endif
        endif
        arp =arp *iswtch
        brc = arp / arc
c
        if (so(i) .lt. 0.1) then
          vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
          aa(1)=viod(lp,jrc-3)/r(jrc-3) + vid(jrc-3)
          aa(2)=viod(lp,jrc-2)/r(jrc-2) + vid(jrc-2)
          aa(3)=viod(lp,jrc-1)/r(jrc-1) + vid(jrc-1)
          aa(4)=vrc
          aa(5)=viod(lp,jrc+1)/r(jrc+1) + vid(jrc+1)
          aa(6)=viod(lp,jrc+2)/r(jrc+2) + vid(jrc+2)
          aa(7)=viod(lp,jrc+3)/r(jrc+3) + vid(jrc+3)
       else
          vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
          aa(1)=viou(lp,jrc-3)/r(jrc-3) + viu(jrc-3)
          aa(2)=viou(lp,jrc-2)/r(jrc-2) + viu(jrc-2)
          aa(3)=viou(lp,jrc-1)/r(jrc-1) + viu(jrc-1)
          aa(4)=vrc
          aa(5)=viou(lp,jrc+1)/r(jrc+1) + viu(jrc+1)
          aa(6)=viou(lp,jrc+2)/r(jrc+2) + viu(jrc+2)
          aa(7)=viou(lp,jrc+3)/r(jrc+3) + viu(jrc+3)
        endif
        rr(1)=r(jrc-3)-r(jrc)
        rr(2)=r(jrc-2)-r(jrc)
        rr(3)=r(jrc-1)-r(jrc)
        rr(4)=zero
        rr(5)=r(jrc+1)-r(jrc)
        rr(6)=r(jrc+2)-r(jrc)
        rr(7)=r(jrc+3)-r(jrc)
        call polcoe(rr,aa,7,coe)
        vap = coe(2)
        vapp= 2*coe(3)
c
c   Set up matrix without the d2(potential(0)/dr2=0 condition
c   to find an intial guess for gamma.
c
        delta=zero
        bj(1)=log(arc/rc1**lp)
        bj(2)=brc-lp/rc1
        bj(3)=vrc-ev(i)+(lp/rc1)**2-brc**2
        vt=vrc-ev(i)+lp*(lp-1)/rc2
        bj(4)=vap-2*(vt*brc+lp**2/rc3-brc**3)
        bj(5)=vapp-2*((vap-2*lp*(lp-1)/rc3)*brc+(vt-3*brc**2)*
     1  (vt-brc**2)-3*lp**2/rc4)
        aj(1,1)=rc2
        aj(1,2)=rc4
        aj(1,3)=rc5
        aj(1,4)=rc6
        aj(1,5)=rc7
        aj(2,1)=2*rc1
        aj(2,2)=4*rc3
        aj(2,3)=5*rc4
        aj(2,4)=6*rc5
        aj(2,5)=7*rc6
        aj(3,1)=2*one
        aj(3,2)=12*rc2
        aj(3,3)=20*rc3
        aj(3,4)=30*rc4
        aj(3,5)=42*rc5
        aj(4,1)=zero
        aj(4,2)=24*rc1
        aj(4,3)=60*rc2
        aj(4,4)=120*rc3
        aj(4,5)=210*rc4
        aj(5,1)=zero
        aj(5,2)=24*one
        aj(5,3)=120*rc1
        aj(5,4)=360*rc2
        aj(5,5)=840*rc3
        call gaussj(aj,5,5,bj,1,1)
        gamma=bj(1)
        alpha=bj(2)
        alpha1=bj(3)
        alpha2=bj(4)
        alpha3=bj(5)
c
c  Start iteration loop to find delta, uses false postion.
c
        do 150 j=1,50
c
c  Generate pseudo wavefunction-note missing factor exp(delta).
c
          do 110 k=1,jrc
            rp=r(k)
            r2=rp*rp
            polyr = r2*((((alpha3*rp+alpha2)*rp+
     1       alpha1)*rp+ alpha)*r2+gamma)
            ar(k) = iswtch * rp**lp * exp(polyr)
 110      continue
c
c  Integrate pseudo charge density from r = 0 to rc.
c
          ll = 2
          cdps = - ar(jrc) * ar(jrc) * rab(jrc)
          if (jrc .ne. 2*(jrc/2)) then
            do 120 k=jrc,1,-1
              cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 120        continue
          else
            do 121 k=jrc,4,-1
              cdps = cdps +  ll * ar(k) * ar(k) * rab(k)
              ll = 6 - ll
 121        continue
            cdps = cdps - ar(4) * ar(4) * rab(4)
            cdps = cdps + 9 * ( ar(1) * ar(1) * rab(1) +
     1       3 * ar(2) *ar(2) * rab(2) +
     2       3 * ar(3) *ar(3) * rab(3) +
     3       ar(4) * ar(4) * rab(4))/8
          endif
          cdps = cdps/3
c
c   Calculate new delta
c
          fdnew = log(cdrc/cdps) - 2*delta
          if (abs(fdnew) .lt. small) goto 160
          if (j .eq. 1) then
            ddelta=-one/2
          else
            ddelta = - fdnew * ddelta / (fdnew-fdold)
          endif
          delta = delta + ddelta
          bj(1)=log(arc/rc1**lp)-delta
          bj(2)=brc-lp/rc1
          bj(3)=vrc-ev(i)+(lp/rc1)**2-brc**2
          vt=vrc-ev(i)+lp*(lp-1)/rc2
          bj(4)=vap-2*(vt*brc+lp**2/rc3-brc**3)
          bj(5)=vapp-2*((vap-2*lp*(lp-1)/rc3)*brc+(vt-3*brc**2)*
     1     (vt-brc**2)-3*lp**2/rc4)
          aj(1,1)=rc2
          aj(1,2)=rc4
          aj(1,3)=rc5
          aj(1,4)=rc6
          aj(1,5)=rc7
          aj(2,1)=2*rc1
          aj(2,2)=4*rc3
          aj(2,3)=5*rc4
          aj(2,4)=6*rc5
          aj(2,5)=7*rc6
          aj(3,1)=2*one
          aj(3,2)=12*rc2
          aj(3,3)=20*rc3
          aj(3,4)=30*rc4
          aj(3,5)=42*rc5
          aj(4,1)=zero
          aj(4,2)=24*rc1
          aj(4,3)=60*rc2
          aj(4,4)=120*rc3
          aj(4,5)=210*rc4
          aj(5,1)=zero
          aj(5,2)=24*one
          aj(5,3)=120*rc1
          aj(5,4)=360*rc2
          aj(5,5)=840*rc3
          call gaussj(aj,5,5,bj,1,1)
          gamma=bj(1)
          alpha=bj(2)
          alpha1=bj(3)
          alpha2=bj(4)
          alpha3=bj(5)
          fdold = fdnew
 150    continue
c
c  End iteration loop for delta.
c
        write(6,1020)lp-1
        call ext(820+lp)
 1020 format(//,'error in pseudt - nonconvergence in finding',
     1 /,' starting delta for angular momentum ',i1)
c
c  Bracket the correct gamma, use gamma and -gamma
c  from above as intial brackets, expands brackets
c  until a root is found..
c
 160    x1=gamma
        x2=-gamma
        alpha4=zero
c
        call zrbact(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1   rc8,lp,arc,brc,vrc,vap,vapp,ev(i),cdrc,r,rab,
     2   jrc,delta,gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c  Iteration loop to find correct gamma, uses
c  bisection to find gamma.
c
        call rtbist(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1   rc8,lp,arc,brc,vrc,vap,vapp,ev(i),cdrc,r,rab,jrc,delta,
     2   gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c  Augment charge density and invert schroedinger equation
c  to find new potential.
c
        expd = exp(delta)
        if (so(i) .lt. 0.1) then
          do 169 j=1,jrc
            poly = r(j)*r(j)*(((((alpha4*r(j)+alpha3)
     1       *r(j)+alpha2)*r(j)+alpha1)*r(j)+alpha)*r(j)**2+gamma)
            ar(j) = iswtch * r(j)**lp * expd * exp(poly)
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
            xlamda=((((8*alpha4*r(j)+7*alpha3)*r(j)
     1       +6*alpha2)*r(j)+5*alpha1)*r(j)+4*alpha)*r(j)**2+
     2       2*gamma
            vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1       +((((56*alpha4*r(j)+42*alpha3)*r(j)
     2       +30*alpha2)*r(j)+20*alpha1)*r(j)+12*alpha)*r(j)**2
     3       +2*gamma
            viod(lp,j) = (vj-vid(j)) * r(j)
 169      continue
          do 168 j=jrc+1,nr
            vod(j) = vod(j) + zo(i)*ar(j)*ar(j)
 168      continue
        else
          do 170 j=1,jrc
            poly = r(j)*r(j)*(((((alpha4*r(j)+alpha3)
     1       *r(j)+alpha2)*r(j)+alpha1)*r(j)+alpha)*r(j)**2+gamma)
            ar(j) = iswtch * r(j)**lp * expd * exp(poly)
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
            xlamda=((((8*alpha4*r(j)+7*alpha3)*r(j)
     1       +6*alpha2)*r(j)+5*alpha1)*r(j)+4*alpha)*r(j)**2+
     2       2*gamma
            vj = ev(i) + xlamda * (2 * lp + xlamda * r(j)**2)
     1       +((((56*alpha4*r(j)+42*alpha3)*r(j)
     2       +30*alpha2)*r(j)+20*alpha1)*r(j)+12*alpha)*r(j)**2
     3       +2*gamma
            viou(lp,j) = (vj-viu(j)) * r(j)
 170      continue
          do 171 j=jrc+1,nr
            vou(j) = vou(j) + zo(i)*ar(j)*ar(j)
 171      continue
        endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
        ist=1
        if (ar(nr-85) .lt. zero) ist=-1
        call potrw(ar,r,nr-85,lo(i),0,ist)
        if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
        call wtrans(ar,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c
        write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cdrc,delta
 180  format(1x,i1,a1,f6.1,5f12.6)
 190  continue
c
c  End loop over valence orbitals.
c
c  Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 200 i=ncp,norb
        nops(i) = lo(i) + 1
        zval = zval + zo(i)
 200  continue
      zion = zval+znuc-zel
      if (zval .ne. zero) zratio=zion/zval
      do 210 i=1,nr
        cdd(i) = vod(i)
 210  continue
      do 211 i=1,nr
        cdu(i) = vou(i)
 211  continue
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used
c
      if (ifcore .ne. 0) then
        ac = zero
        bc = zero
        icore = 1
        if (cfac .le. zero .or. zratio .eq. zero) then
          write(6,280) r(icore),ac,bc
        else
          if (rcfac .le. zero) then
            do 220 i=nr,2,-1
              if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220        continue
          else
            do 221 i=nr,2,-1
              if (r(i) .le. rcfac ) goto 230
 221        continue
          endif
 230      icore = i
          cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
          tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
          rbold = tpfive
          do 240 i=1,50
            rbnew = pi+atan(tanb*rbold)
            if (abs(rbnew-rbold) .lt. .00001) then
              bc = rbnew / r(icore)
              ac = cdc(icore) / (r(icore)*sin(rbnew))
              do 260 j=1,icore
                cdc(j) = ac*r(j)*sin(bc*r(j))
 260          continue
              write(6,280) r(icore),ac,bc
              goto 290
            else
              rbold=rbnew
            endif
 240      continue
          write(6,1030)
          call ext(830)
        endif
      endif
 280  format(//,' core correction used',/,
     1 ' pseudo core inside r =',f6.3,/,' ac =',f6.3,' bc =',f6.3,/)
 1030 format(//,' error in pseudt - noncovergence in finding ',
     1 /,'pseudo-core values')
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
 290  if (ispp .eq. 's') then
        blank='s'
      else
        blank=' '
      endif
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      ecut=ecuts
      do 315 i=ncp,norb
        lp = lo(i)+1
        if (so(i) .lt. 0.1) then
          do 300 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 300      continue
          rcut(i-ncore) = r(jcut)
          do 310 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 310      continue
          do 311 j=2,nr
            v(j) = viod(lp,j)/r(j)
 311      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        else
          do 312 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 312      continue
          rcut(i-ncore) = r(jcut)
          do 313 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 313      continue
          do 314 j=2,nr
            v(j) = viou(lp,j)/r(j)
 314      continue
c
c  njtj  ***  plotting routines ***
c
          call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
          call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
        endif
 315  continue
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a
c   marker 'zio' for latter plotting
c
      write(3,4559)
      write(3,4560) zion
 4559 format(1x,'marker zio')
 4560 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c

c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (ispp .eq. 's' ) then
        do 500 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 505 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 505        continue
          else
            do 506 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 506        continue
          endif
 500    continue
      endif
c
      do 320 i=2,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 320  continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
      call dsolv2(0,1,blank,ifcore,lmax,nr,a,b,r,rab,
     1 norb-ncore,0,nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 znuc,cdd,cdu,cdc,viod,viou,vid,viu,ev(ncp),ek(ncp),
     3 ep(ncp),wk1,wk2,wk3,wk4,wk5,wk6,wk7,evi(ncp))
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,326) (ev(i),i=ncp,norb)
 325  format(//,' test of eigenvalues',//,' rcut =',8(2x,a1,f7.2))
 326  format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,330)
 330  format(///,' l    vps(0)    vpsmin      at r',/)
      do 370 i=1,lmax
        if (indd(i)+indu(i) .eq. 0) goto 370
        if (indd(i) .ne. 0) then
          vpsdm = zero
          do 350 j=2,nr
            if (r(j) .lt. .00001) goto 350
            vps = viod(i,j)/r(j)
            if (vps .lt. vpsdm) then
              vpsdm = vps
              rmind = r(j)
            endif
 350      continue
          write(6,360) il(i),viod(i,2)/r(2),vpsdm,rmind
        endif
        if (indu(i) .ne. 0) then
          vpsum = zero
          do 351 j=2,nr
            if (r(j) .lt. .00001) goto 351
            vps = viou(i,j)/r(j)
            if (vps .lt. vpsum) then
              vpsum = vps
              rminu = r(j)
            endif
 351      continue
          write(6,360) il(i),viou(i,2)/r(2),vpsum,rminu
        endif
 360  format(1x,a1,3f10.3)
 370  continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,
     1 nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1) = 'atom-lda  '
      call zedate(iray(2))
      iray(3) = ' Troullier'
      iray(4) = ' - Martins'
      iray(5) = ' potential'
      iray(6) = '          '
c
c  Encode the title array.
c
      do 390 i=1,7
        ititle(i) = '          '
 390  continue
      do 420 i=1,lmax
        if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
        zelu = zero
        zeld = zero
        if (indd(i) .ne. 0) then
          noi = no(indd(i))
          zeld = zo(indd(i))
        endif
        if (indu(i) .ne. 0) then
          noi = no(indu(i))
          zelu = zo(indu(i))
        endif
        zelt = zeld + zelu
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt
         write(ititle(2*i),401)ispp,rc(i)
 400     format(i1,a1,'(',f6.2,')')
 401     format(a1,' rc=',f5.2)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld
         write(ititle(2*i),411)zelu,ispp,rc(i)
 410     format(i1,a1,'  (',f4.2,',')
 411     format(f4.2,')',a1,f4.2)
        endif
 420  continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
        if (indu(1) .eq. 0) goto 429
        indd(1)=indu(1)
        indu(1)=0
        do 428 j=2,nr
          viod(1,j) = viou(1,j)
          viou(1,j) = zero
 428    continue
 429    do 431 i=2,lmax
          if (indd(i) .eq. 0 .or. indu(i) .eq. 0) goto 431
          do 430 j=2,nr
            viodj = viod(i,j)
            viouj = viou(i,j)
            viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
            viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430      continue
 431    continue
      endif
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 450 i=1,lmax
        if (indd(i) .ne. 0) npotd=npotd+1
        if (indu(i) .ne. 0) npotu=npotu+1
 450  continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
        if (ifull .eq. 0) then
          nicore = 'pcec'
        else
          nicore = 'fcec'
        endif
      elseif (ifcore .eq. 2) then
        if (ifull .eq. 0) then
          nicore = 'pche'
        else
          nicore = 'fche'
        endif
      else
        nicore = 'nc  '
      endif
      if (ispp .eq. 's') then
        irel='isp'
      elseif (ispp .eq. 'r') then
        irel='rel'
      else
        irel = 'nrl'
      endif
      rewind 1
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat
c  file (unit=1).
c
      do 460 i=1,lmax
        if (indd(i) .eq. 0) goto 460
        write(1) i-1,(viod(i,j),j=2,nr)
 460  continue
      do 465 i=1,lmax
        if (indu(i) .eq. 0) goto 465
        write(1) i-1,(viou(i,j),j=2,nr)
 465  continue
c
c  Write the charge densities to the current pseudo.dat
c  file (unit=1).
c
      if (ifcore .eq. 0) then
        write(1) (zero,i=2,nr)
      else
        write(1) (cdc(i),i=2,nr)
      endif
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
      end
C
C
C
      subroutine pseudv(itype,icorr,ispp,lmax,nr,a,b,r,rab,
     1 nameat,norb,ncore,no,lo,so,zo,znuc,zel,cdd,cdu,cdc,
     2 viod,viou,vid,viu,vod,vou,etot,ev,ek,ep,wk1,wk2,wk3,
     3 wk4,wk5,wk6,wk7,f,g,nops,v,ar,br,arps,wkb,evi)
c
c *************************************************************
c *                                                           *
c *     This routine was written by Norman J. Troullier Jr.   *
c *   Nov. 1989, while at the U. of Minnesota, all            *
c *   comments concerning this routine should be directed     *
c *   to him.                                                 *
c *                                                           *
c *     troullie@128.101.224.101                              *
c *     troullie@csfsa.cs.umn.edu                             *
c *     612 625-0392                                          *
c *                                                           *
c *     pseudv generates a pseudopotential using the          *
c *   scheme of D. Vanderbilt, ref. Physical Review B,        *
c *   vol. 32, num 12, page 8412.                             *
c *   The general format of this routine is the same as the   *
c *   pseudo, pseudk and pseudt routines.  Output/input is    *
c *   compatible.                                             *
c *                                                           *
c *************************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
       parameter(zero=0.D0,deltas=1.D-3,tpfive=2.5D0,one=1.D0,two=2.D0)
       parameter(small=1.D-32,small2=1.D-8,small3=1.D-16,pzfive=0.05D0)
       parameter(pfive=0.5D0,small4=1.D-6,ai=2*137.0360411D0)
       parameter(onepf=1.5D0,oneh=100.D0)
Cray       parameter(zero=0.0,deltas=1.E-3,tpfive=2.5,one=1.0,two=2.D0)
Cray       parameter(small=1.E-32,small2=1.E-8,small3=1.E-16,pzfive=0.5)
Cray       parameter(pfive=0.5,small4=1.E-6,ai=2*137.0360411)
Cray       parameter(onepf=1.5,oneh=100.0)
c
       character*1 ispp,blank,il(5)
       character*2 icorr,nameat
       character*3 irel
       character*4 nicore
       character*10 ititle(7),iray(6)
c
      dimension r(nr),rab(nr),no(norb),lo(norb),so(norb),zo(norb),
     1 cdd(nr),cdu(nr),cdc(nr),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),vod(nr),vou(nr),ev(norb),ek(norb),ep(norb),
     3 wk1(nr),wk2(nr),wk3(nr),wk4(nr),wk5(nr),wk6(nr),wk7(nr),
     4 wkb(6*nr),f(nr),g(nr),nops(norb),v(nr),ar(nr),br(nr),
     5 arps(nr),evi(norb)
c
      dimension etot(10),indd(5),indu(5),rc(5),rcut(10),ab(5),
     1 rr(5),coe(5),bj(3),aj(3,3)
c
       data il/'s','p','d','f','g'/
       do 3 i=1,5
         indd(i)=0
         indu(i)=0
 3     continue
       if (ncore .eq. norb) return
       ifcore = itype-1
       pi = 4*atan(one)
c
c  Spin-polarized potentails should be unscreened with
c  a spin-polarized valence charge.  This was not
c  done in pseudo and pseudk in earlier versions
c  of this program.
c
      if (ispp .eq. 's' ) then
        blank = 's'
      else
        blank = ' '
      endif
c
c  read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac
c
c    cfac is used for the pseudocore - the pseudocore stops where
c  the core charge density equals cfac times the renormalized
c  valence charge density (renormalized to make the atom neutral).
c  If cfac is input as negative, the full core charge is used,
c  if cfac is input as zero, it is set equal to one.
c    rcfac is used for the pseudocore cut off radius.  If set
c  to less then or equal to zero cfac is used.  cfac must be
c  set to greater then zero.
c
      read(5,10) (rc(i),i=1,5),cfac,rcfac
 10   format(7f10.5)
      if (cfac .eq. zero) cfac=one
c
c   Reset vod and vou to zero.  They are here used to store
c   the pseudo valence charge density.
c
       do 15 i=1,nr
         vod(i) = zero
         vou(i) = zero
 15    continue
c
c  Print the heading.
c
       write(6,20) nameat
 20    format(//,a2,' Pseudopotential Vanderbilt generation',/,1x,
     1  50('-'),//,' nl    s    eigenvalue',6x,'rc',4x,6x,'cl',
     2  9x,'gamma',7x,'delta',/)
c
c      start loop over valence orbitals
c
       ncp = ncore+1
       do 190 i=ncp,norb
         lp = lo(i) + 1
         llp = lo(i)*lp
         if (so(i) .lt. 0.1) then
           if (indd(lp) .ne. 0) then
             write(6,1000)lp-1
             call ext(800+lp)
           else
             indd(lp) = i
           endif
         else
           if (indu(lp) .ne. 0) then
             write(6,1010)lp-1
             call ext(810+lp)
           else
             indu(lp) = i
           endif
         endif
 1000 format(//,'error in pseudv - two down spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
 1010 format(//,'error in pseudv - two up spin orbitals of the same ',
     1 /,'angular momentum (',i1,') exist')
c
c      find all electron wave function
c
         do 25 j=1,nr
           ar(j)=zero
 25      continue
         if (so(i) .lt. 0.1) then
           do 27 j=2,nr
             v(j) = viod(lp,j)/r(j) + vid(j)
 27        continue
         else
           do 30 j=2,nr
             v(j) = viou(lp,j)/r(j) + viu(j)
 30        continue
         endif
         if (ispp .ne. 'r') then
           do 32 j=2,nr
             v(j) = v(j) + llp/r(j)**2
 32        continue
         endif
         if (ispp .ne. 'r') then
           call difnrl(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1      znuc,viod,viou,vid,viu,ev,iflag,wk1,wk2,wk3,evi)
         else
           call difrel(0,i,v,ar,br,lmax,nr,a,b,r,rab,norb,no,lo,so,
     1      znuc,viod,viou,vid,viu,ev,wk1,wk2,wk3,wk4,evi)
         endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save an usefull number of points
c  of the wave function to make a plot.  The info is
c  written to the current plot.dat file.
c
         ist=1
         if (ar(nr-85) .lt. zero) ist=-1
         call potrw(ar,r,nr-85,lo(i),1,ist)
c
c  njtj  ***  user should adjust for their needs  ***
c
c  Find the last zero and extremum.
c
         ka = lo(i)+1
         if (so(i) .lt. 0.1 .and. lo(i) .ne. 0) ka=-lo(i)
         nextr = no(i)-lo(i)
         rzero = zero
         arp = br(2)
c
         if (ispp .eq. 'r') then
           if (so(i) .lt. 0.1) then
             arp = ka*ar(2)/r(2) + (ev(i) - viod(lp,2)/r(2)
     1        - vid(2) + ai*ai) * br(2) / ai
           else
             arp = ka*ar(2)/r(2) + (ev(i) - viou(lp,2)/r(2)
     1        - viu(2) + ai*ai) * br(2) / ai
           endif
         endif
c
         do 40 j=3,nr-7
           if (nextr .eq. 0) goto 50
           if (ar(j-1)*ar(j) .le. zero)
     1      rzero = (ar(j)*r(j-1)-ar(j-1)*r(j)) / (ar(j)-ar(j-1))
           arpm = arp
           arp = br(j)
c
           if (ispp .eq. 'r') then
             if (so(i) .lt. 0.1) then
               arp = ka*ar(j)/r(j) + (ev(i) - viod(lp,j)/r(j)
     1          - vid(j) + ai*ai) * br(j) / ai
             else
               arp = ka*ar(j)/r(j) + (ev(i) - viou(lp,j)/r(j)
     1          - viu(j) + ai*ai) * br(j) / ai
             endif
           endif
c
           if (arp*arpm .gt. zero) goto 40
           rextr = (arp*r(j-1)-arpm*r(j)) / (arp-arpm)
           nextr = nextr - 1
 40      continue
c
c  Check rc, if outside bounds reset.
c
 50      if (rzero .lt. r(2)) rzero = r(2)
         if (rc(lp) .gt. rzero .and. rc(lp) .lt. rextr) goto 60
         if (rc(lp) .ge. rzero) write(6,2001)rc(lp),rextr
 2001  format(/,'Warning, the Core radius =',f5.2,
     1  /,' is larger then wave function',
     1  ' extrema position =',f5.2,/)
         if (rc(lp) .lt. zero) rc(lp) = rzero - rc(lp)*(rextr-rzero)
c
c  Find the index for grid point closest to 1.5*rc.
c  Find the index for 3*rc which is used for matching norms.
c
 60      rcopf= onepf*rc(lp)
         do 71 j=1,nr
           if (r(j) .le. rcopf) then
             jrc=j
           endif
           if (r(j) .lt. 3*rc(lp)) then
             j3rc = j
           endif
 71      continue
c
c  Reset the n quantum numbers.
c
         do 70 j=1,norb
           nops(j) = 0
 70      continue
         nops(i) = lp
c
c   Set up potential vl1, first find true potential,
c   its first and second derivative at rc.  Store new
c   potential(unscreen it first, screening added back
c   in dsolv2).
c
         if (so(i) .lt. 0.1) then
           vrc = viod(lp,jrc)/r(jrc) + vid(jrc)
           ab(1)=viod(lp,jrc-2)/r(jrc-2) + vid(jrc-2)
           ab(2)=viod(lp,jrc-1)/r(jrc-1) + vid(jrc-1)
           ab(3)=vrc
           ab(4)=viod(lp,jrc+1)/r(jrc+1) + vid(jrc+1)
           ab(5)=viod(lp,jrc+2)/r(jrc+2) + vid(jrc+2)
         else
           vrc = viou(lp,jrc)/r(jrc) + viu(jrc)
           ab(1)=viou(lp,jrc-2)/r(jrc-2) + viu(jrc-2)
           ab(2)=viou(lp,jrc-1)/r(jrc-1) + viu(jrc-1)
           ab(3)=vrc
           ab(4)=viou(lp,jrc+1)/r(jrc+1) + viu(jrc+1)
           ab(5)=viou(lp,jrc+2)/r(jrc+2) + viu(jrc+2)
         endif
         rr(1)=r(jrc-2)-r(jrc)
         rr(2)=r(jrc-1)-r(jrc)
         rr(3)=zero
         rr(4)=r(jrc+1)-r(jrc)
         rr(5)=r(jrc+2)-r(jrc)
         call polcoe(rr,ab,5,coe)
         vap = coe(2)
         vapp= 2*coe(3)
         bj(1)=vrc
         bj(2)=vap
         bj(3)=vapp
         aj(1,1)=one
         aj(2,1)=zero
         aj(3,1)=zero
         aj(1,2)=r(jrc)**2
         aj(2,2)=2*r(jrc)
         aj(3,2)=2*one
         aj(1,3)=r(jrc)**4
         aj(2,3)=4*r(jrc)**3
         aj(3,3)=12*r(jrc)**2
         call gaussj(aj,3,3,bj,1,1)
         b0=bj(1)
         b2=bj(2)
         b4=bj(3)
         if (so(i) .lt. 0.1) then
           do 82 j=1,jrc
             viod(lp,j)=((b0+b2*r(j)**2+b4*r(j)**4)-vid(j))*r(j)
 82        continue
         else
           do 83 j=1,jrc
             viou(lp,j)=((b0+b2*r(j)**2+b4*r(j)**4)-viu(j))*r(j)
 83        continue
         endif
c
c  Set up the functions f(r/rc) and g(r/rc) and  modify the ionic potential.
c
         if (lp .eq. 1) then
           dcl = sqrt(znuc)
         else
           dcl=-2*one*lp*llp
         endif
         cl=dcl
         sinhb2=(sinh(one))**2
c
         do 80 j=1,nr
           rrc = r(j)/rc(lp)/onepf
           f(j)=oneh**(-((sinh(rrc))**2)/sinhb2)
           if (f(j) .lt. small2) f(j)=zero
           g(j) =  f(j)
 80      continue
         if (so(i) .lt. 0.1) then
           do 81 j=2,nr
             viod(lp,j)=viod(lp,j)+dcl*f(j)*r(j)
 81        continue
         else
           do 84 j=2,nr
             viou(lp,j)=viou(lp,j)+dcl*f(j)*r(j)
 84        continue
         endif
         dcl=dcl/2
c
c   Start the iteration loop to find cl.
c
         eviae = ev(i)
         devold = zero
         do 130 j=1,100
           call dsolv2(j,2,blank,ifcore,lmax,
     1      nr,a,b,r,rab,norb,ncore,nops,lo,so,zo,znuc,cdd,cdu,cdc,
     2      viod,viou,vid,viu,ev,ek,ep,wk1,wk2,wk3,wk4,wk5,wk6,
     3      wk7,evi)
           dev = eviae-ev(i)
c
c    The abs(dev-devold) condition was added to eliminate
c    division by zero errors in the calculation of
c    dcl = -dev*dcl / (dev-devold).
c
       if ((abs(dev) .lt. small2 .or. abs(dev-devold) .lt. small3)
     1  .and. j .ne. 1) then
         goto 140
       else
         if (j  .gt. 15 .or. abs(dev) .lt. 0.001) then
c
c  Use newton raphson iteration to change cl.
c
           dcl = -dev*dcl / (dev-devold)
         else
           if (dev*dcl .le. zero) then
             dcl=-dcl/4
           endif
         endif
       endif
c
c  Find the new potential.
c
       if (so(i) .lt. 0.1) then
         do 110 k=2,nr
           viod(lp,k) = viod(lp,k) + dcl*r(k)*f(k)
 110     continue
       else
         do 111 k=2,nr
           viou(lp,k) = viou(lp,k) + dcl*r(k)*f(k)
 111     continue
       endif
       cl = cl + dcl
       devold = dev
 130   continue
c
c  End the iteration loop for cl.
c
       call ext(820+lp)
c
c  Find the new pseudo-wavefunction.
c
 140   if (so(i) .lt. 0.1) then
         do 150 j=2,nr
           v(j) = (viod(lp,j)+llp/r(j))/r(j) + vid(j)
 150     continue
       else
         do 151 j=2,nr
           v(j) = (viou(lp,j)+llp/r(j))/r(j) + viu(j)
 151     continue
       endif
       do 152 j=1,nr
         arps(j)=zero
 152   continue
       call difnrl(0,i,v,arps,br,lmax,nr,a,b,r,rab,norb,
     1  nops,lo,so,znuc,viod,viou,vid,viu,ev,iflag,wk1,
     2  wk2,wk3,evi)
c
c  Compute yl store in g, store ln(arps) in br.
c
       do 155 j=2,nr
         g(j)=arps(j)*f(j)
 155   continue
       do 157 j=2,nr
         br(j)=log(arps(j)+small)
 157   continue
c
c  Compute delta and gamma.
c
       gamma = abs(ar(j3rc)/arps(j3rc)+ar(j3rc+1)/arps(j3rc+1))/2
       ag = zero
       gg = zero
       ll = 4
       do 160 j=2,nr
         ag = ag + ll*arps(j)*g(j)*rab(j)
         gg = gg + ll*g(j)*g(j)*rab(j)
         ll = 6 - ll
 160   continue
       ag = ag/3
       gg = gg/3
       delta = sqrt((ag/gg)**2+(one/gamma**2-one)/gg) - ag/gg
c
c  Modify the pseudo-wavefunction.
c
       do 171 j=2,nr
         arps(j) = gamma*arps(j)*(one+delta*f(j))
 171   continue
c
c     Find d(ln(wl)/dr and store in g().  Note the use of additional
c   given information of the Vanderbilt method, i.e. the use of
c   d(ln(wl)/dr to improve stability.
c
       do 172 j=4,nr-2
         ab(1) = br(j-2)
         ab(2) = br(j-1)
         ab(3) = br(j)
         ab(4) = br(j+1)
         ab(5) = br(j+2)
         rr(1)=r(j-2)-r(j)
         rr(2)=r(j-1)-r(j)
         rr(3)=zero
         rr(4)=r(j+1)-r(j)
         rr(5)=r(j+2)-r(j)
         call polcoe(rr,ab,5,coe)
         g(j)=coe(2)
 172   continue
       g(nr-1)=g(nr-2)
       g(nr)=g(nr-2)
       ab(1) = g(4)
       ab(2) = g(5)
       ab(3) = g(6)
       ab(4) = g(7)
       ab(5) = g(8)
       rr(1)=r(4)-r(3)
       rr(2)=r(5)-r(3)
       rr(3)=r(6)-r(3)
       rr(4)=r(7)-r(3)
       rr(5)=r(8)-r(3)
       call polcoe(rr,ab,5,coe)
       g(3)=coe(1)
       ab(1) = g(3)
       ab(2) = g(4)
       ab(3) = g(5)
       ab(4) = g(6)
       ab(5) = g(7)
       rr(1)=r(3)-r(2)
       rr(2)=r(4)-r(2)
       rr(3)=r(5)-r(2)
       rr(4)=r(6)-r(2)
       rr(5)=r(7)-r(2)
       call polcoe(rr,ab,5,coe)
       g(2)=coe(1)
c
c   Find constants for inversion.
c
       c3=log(oneh)/onepf/rc(lp)/sinhb2
       c2=2/onepf/rc(lp)*c3
       c1=c3**2
c
c    Modify potential and find total charge density.
c
       if (so(i) .lt. 0.1) then
         do 173 j=2,nr
           vod(j)=vod(j)+zo(i)*arps(j)*arps(j)
 173     continue
       else
         do 174 j=2,nr
           vou(j)=vou(j)+zo(i)*arps(j)*arps(j)
 174     continue
       endif
       if (so(i) .lt. 0.1) then
         do 175 j=2,nr
           xr=two*r(j)/rc(lp)/onepf
           sinhxr=sinh(xr)
           coshxr=cosh(xr)
           viod(lp,j)=viod(lp,j)+delta*f(j)/(one+delta*f(j))*
     1      (c1*(sinhxr)**2-c2*coshxr-2*c3*sinhxr*g(j))*r(j)
 175     continue
       else
         do 176 j=2,nr
           xr=two*r(j)/rc(lp)/onepf
           sinhxr=sinh(xr)
           coshxr=cosh(xr)
           viou(lp,j)=viou(lp,j)+delta*f(j)/(one+delta*f(j))*
     1      (c1*(sinhxr)**2-c2*coshxr-2*c3*sinhxr*g(j))*r(j)
 176     continue
       endif
c
c  njtj  ***  plotting routines ***
c  potrw is called to save a usefull number of points
c  of the pseudowave function to make a plot.  The
c  info is written to the current plot.dat file.
c  wtrans is called to fourier transform the the pseudo
c  wave function and save it to the current plot.dat file.
c
        ist=1
        if (arps(nr-85) .lt. zero) ist=-1
        call potrw(arps,r,nr-85,lo(i),0,ist)
        if (ev(i) .eq. zero .or. evi(i) .ne. zero) ist=2
        call wtrans(arps,r,nr,rab,lo(i),ist,wk1)
c
c  njtj  ***  user should adjust for their needs  ***
c
       write(6,180) nops(i),il(lp),so(i),ev(i),rc(lp),cl,gamma,delta
 180   format(1x,i1,a1,f6.1,2f12.6,f12.3,2f12.4)
 190   continue
c
c   End the loop over the valence orbitals.
c
c    Reset the n quantum numbers to include all valence orbitals.
c  Compute the ratio between the valence charge present and the
c  valence charge of a neutral atom.
c  Transfer pseudo valence charge to charge array
c
      zval = zero
      zratio = zero
      do 200 i=ncp,norb
        nops(i) = lo(i) + 1
        zval = zval + zo(i)
 200  continue
      zion = zval+znuc-zel
      if (zval .ne. zero) zratio=zion/zval
      vod(1)=zero
      vou(1)=zero
      do 210 i=1,nr
        cdd(i) = vod(i)
 210  continue
      do 211 i=1,nr
        cdu(i) = vou(i)
 211  continue
c
c  If a core correction is indicated construct pseudo core charge
c  cdc(r) = ac*r * sin(bc*r) inside r(icore)
c  if cfac < 0 or the valence charge is zero the full core is used
c
      if (ifcore .ne. 0) then
        ac = zero
        bc = zero
        icore = 1
        if (cfac .le. zero .or. zratio .eq. zero) then
          write(6,280) r(icore),ac,bc
        else
          if (rcfac .le. zero) then
            do 220 i=nr,2,-1
              if (cdc(i) .gt. cfac*zratio*(cdd(i)+cdu(i))) goto 230
 220        continue
          else
            do 221 i=nr,2,-1
              if (r(i) .le. rcfac ) goto 230
 221        continue
          endif
 230      icore = i
          cdcp = (cdc(icore+1)-cdc(icore)) / (r(icore+1)-r(icore))
          tanb = cdc(icore) / (r(icore)*cdcp-cdc(icore))
          rbold = tpfive
          do 240 i=1,50
            rbnew = pi+atan(tanb*rbold)
            if (abs(rbnew-rbold) .lt. .00001) then
              bc = rbnew / r(icore)
              ac = cdc(icore) / (r(icore)*sin(rbnew))
              do 260 j=1,icore
                cdc(j) = ac*r(j)*sin(bc*r(j))
 260          continue
              write(6,280) r(icore),ac,bc
              goto 290
            else
              rbold=rbnew
            endif
 240      continue
          write(6,1030)
          call ext(830)
        endif
      endif
 280  format(//,' core correction used',/,
     1 ' pseudo core inside r =',f6.3,/,' ac =',f6.3,' bc =',f6.3,/)
 1030 format(//,' error in pseudv - noncovergence in finding ',
     1 /,'pseudo-core values')
c
c  End the pseudo core charge.
c  Compute the potential due to pseudo valence charge.
c
c  njtj  ***  NOTE  ***
c  Spin-polarized potentails should be unscreend with
c  spin-polarized valence charge.  This was not
c  done in pseudo and pseudok in earlier versions
c  of this program.
c  njtj  ***  NOTE  ***
c
 290  if (ispp .eq. 's') then
        blank='s'
      else
        blank=' '
      endif
      call velect(0,1,icorr,blank,ifcore,nr,r,rab,zval,
     1 cdd,cdu,cdc,vod,vou,etot,wk1,wk2,wk3,wk4,wk5,wkb)
c
c  Construct the ionic pseudopotential and find the cutoff,
c  ecut should be adjusted to give a reassonable ionic cutoff
c  radius, but should not alter the pseudopotential, ie.,
c  the ionic cutoff radius should not be inside the pseudopotential
c  cutoff radius
c
      ecut=deltas
      do 315 i=ncp,norb
        lp = lo(i)+1
        if (so(i) .lt. 0.1) then
          do 300 j=2,nr
            viod(lp,j)=viod(lp,j) + (vid(j)-vod(j))*r(j)
            vp2z = viod(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 300      continue
          rcut(i-ncore) = r(jcut)
          do 310 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viod(lp,j) = - 2*zion + fcut * (viod(lp,j)+2*zion)
 310      continue
          do 311 j=2,nr
            v(j) = viod(lp,j)/r(j)
 311      continue
        else
          do 312 j=2,nr
            viou(lp,j)=viou(lp,j)+ (viu(j)-vou(j))*r(j)
            vp2z = viou(lp,j) + 2*zion
            if (abs(vp2z) .gt. ecut) jcut = j
 312      continue
          rcut(i-ncore) = r(jcut)
          do 313 j=jcut,nr
            fcut = exp(-5*(r(j)-r(jcut)))
            viou(lp,j) = - 2*zion + fcut * (viou(lp,j)+2*zion)
 313      continue
          do 314 j=2,nr
            v(j) = viou(lp,j)/r(j)
 314      continue
        endif
c
c  njtj  ***  plotting routines ***
c
        call potran(lo(i)+1,v,r,nr,zion,wk1,wk2,wk3)
        call potrv(v,r,nr-120,lo(i))
c
c  njtj  ***  user should adjust for their needs  ***
c
 315  continue
c
c  njtj  ***  plotting routines ***
c   The calls to 1)potran take the fourier transform of
c   the potential and saves it in the current plot.dat file,
c   2)potrv saves the potential in the current plot.dat file
c   3)zion is saved to the current plot.dat file wtih a
c   marker 'zio' for latter plotting
c
      write(3,4559)
      write(3,4560) zion
 4559 format(1x,'marker zio')
 4560 format(2x,f5.2)
c
c  njtj  ***  user should adjust for their needs  ***
c
c   Convert spin-polarized potentials back to nonspin-polarized
c   by occupation weight(zo).  Assumes core polarization is
c   zero, ie. polarization is only a valence effect.
c
      if (ispp .eq. 's' ) then
        do 500 i=ncp,norb,2
          lp = lo(i)+1
          zot=zo(i)+zo(i+1)
          if (zot .ne. zero) then
            do 505 j=2,nr
              viod(lp,j)=(viod(lp,j)*zo(i)+viou(lp,j)
     1         *zo(i+1))/zot
              viou(lp,j)=viod(lp,j)
 505        continue
          else
            do 506 j=2,nr
              viod(lp,j)=viod(lp,j)/2+viou(lp,j)/2
              viou(lp,j)=viod(lp,j)
 506        continue
          endif
 500    continue
      endif
c
      do 320 i=1,nr
        vid(i) = vod(i)
        viu(i) = vou(i)
 320  continue
c
c   Test the pseudopotential self consistency.  Spin-polarized
c   is tested as spin-polarized(since up/down potentials are
c   now the same)
c
       call dsolv2(0,1,blank,ifcore,lmax,
     1  nr,a,b,r,rab,norb,ncore,nops,lo,so,zo,znuc,cdd,cdu,cdc,
     2  viod,viou,vid,viu,ev,ek,ep,wk1,wk2,wk3,wk4,wk5,wk6,
     3  wk7,evi)
c
c  Printout the pseudo eigenvalues after cutoff.
c
      write(6,325) (il(lo(i)+1),rcut(i-ncore),i=ncp,norb)
      write(6,326) (ev(i),i=ncp,norb)
 325  format(//,' test of eigenvalues',//,' rcut =',8(2x,a1,f7.2))
 326  format(' eval =',8(2x,f8.5))
c
c  Printout the data for potentials.
c
      write(6,330)
 330  format(///,' l    vps(0)    vpsmin      at r',/)
      do 370 i=1,lmax
        if (indd(i)+indu(i) .eq. 0) goto 370
        if (indd(i) .ne. 0) then
          vpsdm = zero
          do 350 j=2,nr
            if (r(j) .lt. .00001) goto 350
            vps = viod(i,j)/r(j)
            if (vps .lt. vpsdm) then
              vpsdm = vps
              rmind = r(j)
            endif
 350      continue
          write(6,360) il(i),viod(i,2)/r(2),vpsdm,rmind
        endif
        if (indu(i) .ne. 0) then
          vpsum = zero
          do 351 j=2,nr
            if (r(j) .lt. .00001) goto 351
            vps = viou(i,j)/r(j)
            if (vps .lt. vpsum) then
              vpsum = vps
              rminu = r(j)
            endif
 351      continue
          write(6,360) il(i),viou(i,2)/r(2),vpsum,rminu
        endif
 360  format(1x,a1,3f10.3)
 370  continue
c
c   Print out the energies from etotal.
c
      call etotal(itype,one,nameat,norb-ncore,
     1 nops(ncp),lo(ncp),so(ncp),zo(ncp),
     2 etot,ev(ncp),ek(ncp),ep(ncp))
c
c  Find the jobname and date, date is a machine
c  dependent routine and must be chosen/written/
c  comment in/out in the zedate section.
c
      iray(1)='atom-lda  '
      call zedate(iray(2))
      iray(3) = 'Vanderbilt'
      iray(4) = ' Pseudo - '
      iray(5) = 'potential '
      iray(6) = 'generation'
c
c  Encode the title array.
c
      do 390 i=1,7
        ititle(i) = '          '
 390  continue
      do 420 i=1,lmax
        if (indd(i) .eq. 0 .and. indu(i) .eq. 0) goto 420
        zelu = zero
        zeld = zero
        if (indd(i) .ne. 0) then
          noi = no(indd(i))
          zeld = zo(indd(i))
        endif
        if (indu(i) .ne. 0) then
          noi = no(indu(i))
          zelu = zo(indu(i))
        endif
        zelt = zeld + zelu
       if (ispp .ne. 's') then
         write(ititle(2*i-1),400) noi,il(i),zelt
         write(ititle(2*i),401)ispp,rc(i)
 400     format(' ',i1,a1,'(',f5.2,')')
 401     format(a1,' rc=',f5.2)
       else
         write(ititle(2*i-1),410) noi,il(i),zeld
         write(ititle(2*i),411)zelu,ispp,rc(i)
 410     format(i1,a1,'  (',f4.2,',')
 411     format(f4.2,')',a1,f4.2)
        endif
 420  continue
c
c  Construct relativistic sum and difference potentials.
c
      if (ispp .eq. 'r') then
        if (indu(1) .eq. 0) goto 429
        indd(1)=indu(1)
        indu(1)=0
        do 428 j=2,nr
          viod(1,j) = viou(1,j)
          viou(1,j) = zero
 428    continue
 429    do 431 i=2,lmax
          if (indd(i) .eq. 0 .or. indu(i) .eq. 0) goto 431
          do 430 j=2,nr
            viodj = viod(i,j)
            viouj = viou(i,j)
            viod(i,j) = ((i-1)*viodj + i*viouj) / (2*i-1)
            viou(i,j) = 2 * (viouj - viodj) / (2*i-1)
 430      continue
 431    continue
      endif
c
c  Determine the number of  potentials.  Coded them as
c  two digits, where the first digit is the number
c  of down or sum potentials and the second the number of
c  up or difference potentials.
c
      npotd = 0
      npotu = 0
      do 450 i=1,lmax
        if (indd(i) .ne. 0) npotd=npotd+1
        if (indu(i) .ne. 0) npotu=npotu+1
 450  continue
c
c  Write the heading to the current pseudo.dat
c  file (unit=1).
c
      ifull = 0
      if (cfac .le. zero .or. zratio .eq. zero) ifull = 1
      if (ifcore .eq. 1) then
        if (ifull .eq. 0) then
          nicore = 'pcec'
        else
          nicore = 'fcec'
        endif
      elseif (ifcore .eq. 2) then
        if (ifull .eq. 0) then
          nicore = 'pche'
        else
          nicore = 'fche'
        endif
      else
        nicore = 'nc  '
      endif
      if (ispp .eq. 's') then
        irel='isp'
      elseif (ispp .eq. 'r') then
        irel='rel'
      else
        irel = 'nrl'
      endif
      rewind 1
      write(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     1 (ititle(i),i=1,7),npotd,npotu,nr-1,a,b,zion
      write(1) (r(i),i=2,nr)
c
c  Write the potentials to the current pseudo.dat
c  file (unit=1).
c
      do 460 i=1,lmax
        if (indd(i) .eq. 0) goto 460
        write(1) i-1,(viod(i,j),j=2,nr)
 460  continue
      do 465 i=1,lmax
        if (indu(i) .eq. 0) goto 465
        write(1) i-1,(viou(i,j),j=2,nr)
 465  continue
c
c  Write the charge densities to the current pseudo.dat
c  file (unit=1).
c
      if (ifcore .eq. 0) then
        write(1) (zero,i=2,nr)
      else
        write(1) (cdc(i),i=2,nr)
      endif
      write(1) (zratio*(cdd(i)+cdu(i)),i=2,nr)
c
      return
      end
C
C
C
      subroutine rtbis2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1 rc8,lp,arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2 gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c *************************************************************
c *  njtj
c *  Finds the value of gamma for the v"(0)=0 criteria.
c *  The method used is bisection.  This routine
c *  was taken from Numerical Recipes, page 247.
c *  njtj
c *************************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out the implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (jmax=80,pfive=0.5D0,zero=0.D0,xacc=1.D-10)
Cray      parameter (jmax=80,pfive=0.5,zero=0.0,xacc=1.E-10)
c
      dimension r(jrc),rab(jrc),ar(jrc)
c
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f,ar)
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2 alpha,alpha1,alpha2,alpha3,alpha4,fmid,ar)
      if(f*fmid.ge.zero) then
        write(6,4000)
        call ext(840+lp)
      endif
      if(f.lt.zero)then
        gamma=x1
        dx=x2-x1
      else
        gamma=x2
        dx=x1-x2
      endif
      do 11 j=1,jmax
        dx=dx*pfive
        xmid=gamma+dx
        call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1   arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2   xmid,alpha,alpha1,alpha2,alpha3,alpha4,fmid,ar)
        if(fmid.lt.zero)gamma=xmid
        if(abs(dx).lt.xacc .or. fmid.eq. zero) return
11    continue
      write(6,4001)
      call ext(850+lp)
 4000 format(' error in bisection method(rtbistk)',
     1 ' - root must be bracketed.',
     2 /,'a b o r t i n g   p r o g r a m')
 4001 format(' error in bisection method(rtbistk)',
     1 ' - too many bisections used',
     2 /,'a b o r t i n g   p r o g r a m')
      end
C
C
C
      subroutine rtbist(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1 rc8,lp,arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2 gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c *************************************************************
c *  njtj
c *  Finds the value of gamma for the v"(0)=0 criteria.
c *  The method used is bisection.  This routine
c *  was taken from Numerical Recipes, page 247.
c *  njtj
c *************************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out the implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (jmax=80,pfive=0.5D0,zero=0.D0,xacc=1.D-10)
Cray      parameter (jmax=80,pfive=0.5,zero=0.0,xacc=1.E-10)
c
      dimension r(jrc),rab(jrc),ar(jrc)
c
      call gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f,ar)
      call gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2 alpha,alpha1,alpha2,alpha3,alpha4,fmid,ar)
      if(f*fmid.ge.zero) then
        write(6,4000)
        call ext(840+lp)
      endif
      if(f.lt.zero)then
        gamma=x1
        dx=x2-x1
      else
        gamma=x2
        dx=x1-x2
      endif
      do 11 j=1,jmax
        dx=dx*pfive
        xmid=gamma+dx
        call gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1   arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2   xmid,alpha,alpha1,alpha2,alpha3,alpha4,fmid,ar)
        if(fmid.lt.zero)gamma=xmid
        if(abs(dx).lt.xacc .or. fmid.eq. zero) return
11    continue
      write(6,4001)
      call ext(850+lp)
 4000 format(' error in bisection method(rtbistk)',
     1 ' - root must be bracketed.',
     2 /,'a b o r t i n g   p r o g r a m')
 4001 format(' error in bisection method(rtbistk)',
     1 ' - too many bisections used',
     2 /,'a b o r t i n g   p r o g r a m')
      end
C
C
C
       DOUBLE PRECISION FUNCTION SBESSJ(N,X)
       implicit double precision(a-h, o-z)
       PARAMETER(ONE=1.D0,TWO=2.D0,THREE=3.D0,ZERO=0.D0)
       PARAMETER( FIVE = 5.0D0 , TEN = 10.0D0 , FOURTN = 14.0D0 )
C      SPHERICAL BESSEL FUNCTION OF THE FIRST KIND
C

       IF(ABS(X) .GT. 0.001) THEN
         SB0 = SIN(X)/X
       ELSE
         X2 = X*X/TWO
         SB0 = ONE - (X2/THREE)*(ONE - X2/TEN)
       ENDIF
       IF(N .EQ. 0) THEN
         SBESSJ = SB0
       ELSE
         IF(ABS(X) .GT. 0.001) THEN
           SB1 = (SIN(X)/X - COS(X)) / X
         ELSE
           X2 = X*X/TWO
           SB1 = (X/THREE)*(ONE - (X2/FIVE)*(1.0 - X2/FOURTN))
         ENDIF
         IF(N .EQ. 1) THEN
           SBESSJ = SB1
         ELSEIF(X .EQ. ZERO) THEN
           SBESSJ = ZERO
         ELSE
           BY = SB1
           BYM = SB0
           UX = ONE / X
           DO 10 J=1,N-1
             BYP = REAL(2*J+1)*UX*BY - BYM
             BYM = BY
             BY = BYP
 10        CONTINUE
           SBESSJ = BY
         ENDIF
       ENDIF
       RETURN
       END

       SUBROUTINE SPLIFT (X,Y,YP,YPP,N,W,IERR,ISX,A1,B1,AN,BN)
C
      implicit double precision(a-h,o-z)

      PARAMETER (FOUR=4.D0)
CRAY      PARAMETER (FOUR=4.0)
C
C  NJTJ
C  ###  CRAY CONVERSIONS
C  ###    1)Comment out the implicit double precision.
C  ###    2)Switch double precision parameter
C  ###      to single precision parameter
C  ###  CRAY CONVERSIONS
C  NJTJ
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
C         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
C         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
C         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
C         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
C         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
C
C       --INPUT--
C
C         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
C         Y    - ARRAY OF ORDINATES OF DATA
C         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
C                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
C         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
C                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
C                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
C                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
C                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
C                SET TO ONE FOR FASTER EXECUTION.
C         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
C                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
C                OF THE SPLINE AT THE END POINTS (SEE YPP).
C                THE END CONDITION CONSTRAINTS ARE
C                        YPP(1) = A1*YPP(2) + B1
C                AND
C                        YPP(N) = AN*YPP(N-1) + BN
C                WHERE
C                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
C
C                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
C                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
C                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
C                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
C                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
C                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
C                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
C                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
C                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
C                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
C                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
C
C                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
C                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
C                FOLLOWING CONDITIONS.
C
C                A1=-0.5
C
C                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
C
C                AN=-0.5
C
C                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
C
C       --OUTPUT--
C
C         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
C         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
C         IERR - A STATUS CODE
C              --NORMAL CODE
C                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
C              --ABNORMAL CODES
C                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
C                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.
C
C       --WORK--
C
C         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
       DIMENSION X(N),Y(N),YP(N),YPP(N),W(N,3)
C
       IF (N.LT.4) THEN
         IERR = 2
         RETURN
       ENDIF
       NM1  = N-1
       NM2  = N-2
       IF (ISX.GT.0) GO TO 40
       DO 5 I=2,N
         IF (X(I)-X(I-1) .LE. 0) THEN
           IERR = 3
           RETURN
         ENDIF
 5     CONTINUE
C
C     DEFINE THE TRIDIAGONAL MATRIX
C
       W(1,3) = X(2)-X(1)
       DO 10 I=2,NM1
         W(I,2) = W(I-1,3)
         W(I,3) = X(I+1)-X(I)
 10      W(I,1) = 2*(W(I,2)+W(I,3))
       W(1,1) = FOUR
       W(1,3) =-4*A1
       W(N,1) = FOUR
       W(N,2) =-4*AN
C
C     L U DECOMPOSITION
C
       DO 30 I=2,N
         W(I-1,3) = W(I-1,3)/W(I-1,1)
 30    W(I,1) = W(I,1) - W(I,2)*W(I-1,3)
C
C     DEFINE *CONSTANT* VECTOR
C
 40   YPP(1) = 4*B1
      DOLD = (Y(2)-Y(1))/W(2,2)
      DO 50 I=2,NM2
        DNEW   = (Y(I+1) - Y(I))/W(I+1,2)
        YPP(I) = 6*(DNEW - DOLD)
        YP(I)  = DOLD
 50   DOLD = DNEW
      DNEW = (Y(N)-Y(N-1))/(X(N)-X(N-1))
      YPP(NM1) = 6*(DNEW - DOLD)
      YPP(N) = 4*BN
      YP(NM1)= DOLD
      YP(N) = DNEW
C
C     FORWARD SUBSTITUTION
C
      YPP(1) = YPP(1)/W(1,1)
      DO 60 I=2,N
 60   YPP(I) = (YPP(I) - W(I,2)*YPP(I-1))/W(I,1)
C
C     BACKWARD SUBSTITUTION
C
       DO 70 J=1,NM1
         I = N-J
   70 YPP(I) = YPP(I) - W(I,3)*YPP(I+1)
C
C     COMPUTE FIRST DERIVATIVES
C
      YP(1) = (Y(2)-Y(1))/(X(2)-X(1)) - (X(2)-X(1))*(2*YPP(1)
     1  + YPP(2))/6
      DO 80 I=2,NM1
 80   YP(I) = YP(I) + W(I,2)*(YPP(I-1) + 2*YPP(I))/6
      YP(N) = YP(N) + (X(N)-X(NM1))*(YPP(NM1) + 2*YPP(N))/6
C
      IERR = 1
      RETURN
      END
C
C
C
       SUBROUTINE SPLINT (X,Y,YPP,N,XI,YI,YPI,YPPI,NI,KERR)
       implicit double precision (a-h,o-z)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     WRITTEN BY RONDALL E. JONES
C
C     ABSTRACT
C
C         SPLINT EVALUATES A CUBIC SPLINE AND ITS FIRST AND SECOND
C         DERIVATIVES AT THE ABSCISSAS IN XI.  THE SPLINE (WHICH
C         IS DEFINED BY X, Y, AND YPP) MAY HAVE BEEN DETERMINED BY
C         SPLIFT OR SMOO OR ANY OTHER SPLINE FITTING ROUTINE THAT
C         PROVIDES SECOND DERIVATIVES.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YPP(N), XI(NI), YI(NI), YPI(NI), YPPI(NI)
C
C       --INPUT--
C
C         X   - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C         Y   - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C         YPP - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C         N   - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C               THE ARRAYS X, Y, AND YPP MUST BE DIMENSIONED AT LEAST N
C               N MUST BE GREATER THAN OR EQUAL TO 2.
C         XI  - THE ABSCISSA OR ARRAY OF ABSCISSAS (IN ARBITRARY ORDER)
C               AT WHICH THE SPLINE IS TO BE EVALUATED.
C               EACH XI(K) THAT LIES BETWEEN X(1) AND X(N) IS A CASE OF
C               INTERPOLATION.  EACH XI(K) THAT DOES NOT LIE BETWEEN
C               X(1) AND X(N) IS A CASE OF EXTRAPOLATION.  BOTH CASES
C               ARE ALLOWED.  SEE DESCRIPTION OF KERR.
C         NI  - THE NUMBER OF ABSCISSAS AT WHICH THE SPLINE IS TO BE
C               EVALUATED.  IF NI IS GREATER THAN 1, THEN XI, YI, YPI,
C               AND YPPI MUST BE ARRAYS DIMENSIONED AT LEAST NI.
C               NI MUST BE GREATER THAN OR EQUAL TO 1.
C
C       --OUTPUT--
C
C         YI  - ARRAY OF VALUES OF THE SPLINE (ORDINATES) AT XI.
C         YPI - ARRAY OF VALUES OF THE FIRST DERIVATIVE OF SPLINE AT XI
C         YPPI- ARRAY OF VALUES OF SECOND DERIVATIVES OF SPLINE AT XI.
C         KERR- A STATUS CODE
C             --NORMAL CODES
C                1 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI USING ONLY INTERPOLATION.
C                2 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
C                  IN XI, BUT AT LEAST ONE EXTRAPOLATION WAS PERFORMED.
C             -- ABNORMAL CODE
C                3 MEANS THAT THE REQUESTED NUMBER OF EVALUATIONS, NI,
C                  WAS NOT POSITIVE.
C
       DIMENSION X(N),Y(N),YPP(N),XI(NI),YI(NI),YPI(NI),YPPI(NI)
C
C     CHECK INPUT
C
      IF (NI) 1,1,2
 1    CONTINUE
C    1 CALL ERRCHK(67,67HIN SPLINT,  THE REQUESTED NUMBER OF INTERPOLATI
C     1NS WAS NOT POSITIVE)
      KERR = 3
      RETURN
    2 KERR = 1
      NM1= N-1
C
C     K IS INDEX ON VALUE OF XI BEING WORKED ON.  XX IS THAT VALUE.
C     I IS CURRENT INDEX INTO X ARRAY.
C
       K  = 1
       XX = XI(1)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IL = 1
       IR = N
C
C     BISECTION SEARCH
C
   10 I  = (IL+IR)/2
       IF (I.EQ.IL) GO TO 100
       IF (XX-X(I)) 20,100,30
   20 IR = I
       GO TO 10
   30 IL = I
       GO TO 10
C
C     LINEAR FORWARD SEARCH
C
   50 IF (XX-X(I+1)) 100,100,60
   60 IF (I.GE.NM1) GO TO 80
       I  = I+1
       GO TO 50
C
C     EXTRAPOLATION
C
   80 KERR = 2
      I  = NM1
      GO TO 100
   90 KERR = 2
      I  = 1
C
C     INTERPOLATION
C
  100 H  = X(I+1) - X(I)
       H2 = H*H
       XR = (X(I+1)-XX)/H
       XR2= XR*XR
       XR3= XR*XR2
       XL = (XX-X(I))/H
       XL2= XL*XL
       XL3= XL*XL2
       YI(K) = Y(I)*XR + Y(I+1)*XL
     1       -H2*(YPP(I)*(XR-XR3) + YPP(I+1)*(XL-XL3))/6.0D0
       YPI(K) = (Y(I+1)-Y(I))/H
     1 +H*(YPP(I)*(1.0D0-3.0D0*XR2)-YPP(I+1)*(1.0D0-3.0D0*XL2))/6.0D0
       YPPI(K) = YPP(I)*XR + YPP(I+1)*XL
C
C     NEXT POINT
C
       IF (K.GE.NI) RETURN
       K = K+1
       XX = XI(K)
       IF (XX.LT.X(1)) GO TO 90
       IF (XX.GT.X(N)) GO TO 80
       IF (XX-XI(K-1)) 110,100,50
  110 IL = 1
       IR = I+1
       GO TO 10
C
       END
       SUBROUTINE SPLIQ(X,Y,YP,YPP,N,XLO,XUP,NUP,ANS,IERR)
C
C
C  NJTJ
C  ###  CRAY CONVERSIONS
C  ###    1)Comment out implicit double precision.
C  ###  CRAY CONVERSIONS
C  NJTJ
C
       implicit double precision (a-h,o-z)
       DIMENSION X(N),Y(N),YP(N),YPP(N),XUP(NUP),ANS(NUP)
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2613
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                    ISSUED BY SANDIA LABORATORIES
C  *                   A PRIME CONTRACTOR TO THE
C  *                UNITED STATES DEPARTMENT OF ENERGY
C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
C  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
C  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
C  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
C  * OWNED RIGHTS.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
C  * PART IS SAND77-1441.
C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     THIS ROUTINE WAS WRITTEN BY M. K. GORDON
C
C     ABSTRACT
C
C     SUBROUTINE SPLIQ INTEGRATES A CUBIC SPLINE (GENERATED BY
C     SPLIFT, SMOO, ETC.) ON THE INTERVALS (XLO,XUP(I)), WHERE XUP
C     IS A SEQUENCE OF UPPER LIMITS ON THE INTERVALS OF INTEGRATION.
C     THE ONLY RESTRICTIONS ON XLO AND XUP(*) ARE
C                XLO .LT. XUP(1),
C                XUP(I) .LE. XUP(I+1)   FOR EACH I .
C     ENDPOINTS BEYOND THE SPAN OF ABSCISSAS ARE ALLOWED.
C     THE SPLINE OVER THE INTERVAL (X(I),X(I+1)) IS REGARDED
C     AS A CUBIC POLYNOMIAL EXPANDED ABOUT X(I) AND IS INTEGRATED
C     ANALYTICALLY.
C
C     DESCRIPTION OF ARGUMENTS
C         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C         E.G.  X(N), Y(N), YP(N), YPP(N), XUP(NUP), ANS(NUP)
C
C      --INPUT--
C
C        X    - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
C               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
C        Y    - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
C               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
C        YP   - ARRAY OF FIRST DERIVATIVES OF THE SPLINE AT ABSCISSAS.
C               USUALLY YP IS THE SAME AS YP IN SPLIFT OR R1 IN SMOO.
C        YPP  - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
C               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
C        N    - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
C        XLO  - LEFT ENDPOINT OF INTEGRATION INTERVALS.
C        XUP  - RIGHT ENDPOINT OR ARRAY OF RIGHT ENDPOINTS OF
C               INTEGRATION INTERVALS IN ASCENDING ORDER.
C        NUP  - THE NUMBER OF RIGHT ENDPOINTS.  IF NUP IS GREATER THAN
C               1, THEN XUP AND ANS MUST BE DIMENSIONED AT LEAST NUP.
C
C      --OUTPUT--
C
C        ANS -- ARRAY OF INTEGRAL VALUES, THAT IS,
C               ANS(I) = INTEGRAL FROM XLO TO XUP(I)
C        IERR -- ERROR STATUS
C                = 1 INTEGRATION SUCCESSFUL
C                = 2 IMPROPER INPUT - N.LT.4 OR NUP.LT.1
C                = 3 IMPROPER INPUT - ABSCISSAS NOT IN
C                        STRICTLY ASCENDING ORDER
C                = 4 IMPROPER INPUT - RIGHT ENDPOINTS XUP NOT
C                        IN ASCENDING ORDER
C                = 5 IMPROPER INPUT - XLO.GT.XUP(1)
C                = 6 INTEGRATION SUCCESSFUL BUT AT LEAST ONE ENDPOINT
C                        NOT WITHIN SPAN OF ABSCISSAS
C              ** NOTE.  ERRCHK PROCESSES DIAGNOSTICS FOR CODES 2,3,4,5
C
C   CHECK FOR IMPROPER INPUT
C
       IERR = 2
       IF(N .LT. 4  .OR.  NUP .LT. 1) THEN
         RETURN
       ENDIF
       NM1 = N-1
       NM2 = N-2
       IERR = 3
       DO 2 I = 1,NM1
         IF(X(I) .GE. X(I+1)) THEN
           RETURN
         ENDIF
 2     CONTINUE
       IF(NUP .NE. 1) THEN
         IERR = 4
         DO 3 I = 2,NUP
           IF(XUP(I-1) .GT. XUP(I)) THEN
             RETURN
           ENDIF
 3       CONTINUE
       ENDIF
       IERR = 5
       IF(XLO .GT. XUP(1)) THEN
         RETURN
       ENDIF
       IERR = 1
       IF(XLO .LT. X(1)  .OR.  XUP(NUP) .GT. X(N)) IERR = 6
C
C   LOCATE XLO IN INTERVAL (X(I),X(I+1))
C
       DO 10 I = 1,NM2
         IF(XLO .LT. X(I+1)) GO TO 20
 10      CONTINUE
       I = NM1
 20    HLO = XLO-X(I)
       HLO2 = HLO*HLO
       HI = X(I+1)-X(I)
       HI2 = HI*HI
       DO 30 J = 1,NUP
         IF(XUP(J) .GT. X(I+1)  .AND.  XLO .LT. X(NM1)) GO TO 40
C
C   COMPUTE SPECIAL CASES OF XUP IN INTERVAL WITH XLO
C
         HUP = XUP(J)-X(I)
         HSUM = HUP+HLO
         HDIFF = HUP-HLO
         HUP2 = HUP*HUP
         SUM = (YPP(I+1)-YPP(I))*HSUM*HDIFF*(HUP2+HLO2)/(24*HI)
         SUM = SUM + YPP(I)*HDIFF*(HUP2+HLO*HUP+HLO2)/6
         SUM = SUM + YP(I)*HDIFF*HSUM/2
         SUM = SUM + Y(I)*HDIFF
 30    ANS(J) = SUM
       RETURN
C
C   COMPUTE INTEGRAL BETWEEN XLO AND X(I+1) AS FOUR TERMS IN TAYLOR
C   POLYNOMIAL AND ADVANCE I TO I+1
C
 40    HDIFF = HI-HLO
       HSUM = HI+HLO
       SUM0 = Y(I)*HDIFF
       SUM1 = YP(I)*HDIFF*HSUM
       SUM2 = YPP(I)*HDIFF*(HI2+HI*HLO+HLO2)
       SUM3 = (YPP(I+1)-YPP(I))*HDIFF*HSUM*(HI2+HLO2)/HI
       I = I+1
C
C   LOCATE EACH XUP(M) IN INTERVAL (X(I),X(I+1))
C
       DO 80 M = J,NUP
 50      IF(XUP(M) .LT. X(I+1)  .OR.  I .EQ. NM1) GO TO 60
C
C   AUGMENT INTEGRAL BETWEEN ABSCISSAS TO INCLUDE INTERVAL
C   (X(I),X(I+1)) AND ADVANCE I TO I+1
C
         HI = X(I+1)-X(I)
         HI2 = HI*HI
         HI3 = HI2*HI
         SUM0 = SUM0 + Y(I)*HI
         SUM1 = SUM1 + YP(I)*HI2
         SUM2 = SUM2 + YPP(I)*HI3
         SUM3 = SUM3 + (YPP(I+1)-YPP(I))*HI3
         I = I+1
         GO TO 50
C
C   INTEGRAL BETWEEN X(I) AND XUP(M) IS ZERO
C
 60      IF(XUP(M) .NE. X(I)) THEN
C
C   COMPUTE INTEGRAL BETWEEN X(I) AND XUP(M) AND EVALUATE
C   TAYLOR POLYNOMIAL IN REVERSE ORDER
C
           HUP = XUP(M)-X(I)
           HUP2 = HUP*HUP
           HUP3 = HUP2*HUP
           HUP4 = HUP3*HUP
           HI = X(I+1)-X(I)
           PSUM0 = Y(I)*HUP
           PSUM1 = YP(I)*HUP2
           PSUM2 = YPP(I)*HUP3
           PSUM3 = (YPP(I+1)-YPP(I))*HUP4/HI
           SUM = (SUM3+PSUM3)/24 + (SUM2+PSUM2)/6
           SUM = SUM + (SUM1+PSUM1)/2
           SUM = SUM + (SUM0+PSUM0)
         ELSE
           SUM = ((SUM3/24 + SUM2/6) + SUM1/2) + SUM0
         ENDIF
 80    ANS(M) = SUM
       RETURN
       END
C
C
C
c
c  ********************************************************
c  *                                                      *
c  *   njtj                                               *
c  *     These are machine dependent routines.            *
c  *   Included are routine for Apollo, Sun,               *
c  *   Vax, and Cray systems.  The user must              *
c  *   1)compile with their systems lines uncommented     *
c  *   or 2)supply their own                              *
c  *   or 3)remove-comment out all references to          *
c  *   these calls in the program.                        *
c  *                                                      *
c  ********************************************************
c
c  ****************Apollo start***********************
c
C
C **************Cray start***********************
C
Cray       SUBROUTINE ZESEC(T)
C
C   GETS CPU TIME IN SECONDS
C   CRAY-2 VERSION
C
Cray       T = SECOND()
Cray       RETURN
Cray       END
C
Cray       SUBROUTINE ZEDATE(BDATE)
C
C    GETS THE DATE (DAY-MONTH-YEAR)
C    CRAY-2 VERSION
C
Cray       CHARACTER*10 BDATE
Cray       CHARACTER*8 ADATE
Cray       CHARACTER*3 MONTH(12)
Cray       CHARACTER*1 DASH,DUM1,DUM2
Cray       DATA DASH/'-'/
Cray       DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
Cray     2  'OCT','NOV','DEC'/
Cray
Cray       WRITE(ADATE,100) DATE()
Cray       READ(ADATE,101) LMONTH,DUM1,LDAY,DUM2,LYEAR
Cray       WRITE(BDATE,102) LDAY,DASH,MONTH(LMONTH),DASH,LYEAR
Cray  100  FORMAT(A8)
Cray  101  FORMAT(I2,A1,I2,A1,I2)
Cray  102  FORMAT(I2,A1,A3,A1,I2,' ')
Cray       RETURN
Cray       END
C
C  *****************Cray end***********************
C
C  *****************Vax start**********************
C
cVax      SUBROUTINE ZESEC(T)
C
C   CALCULATES THE ELAPSED CPU TIME SINCE
C   THE FIRST CALL IN A VAX/VMS SYSTEM
C
cVax      REAL*8 T
cVax      COMMON/ZESEC/IFLAG
cVax      DATA IFLAG /0/
cVax      IF(IFLAG.EQ.0) THEN
cVax        CALL LIB$INIT_TIMER
cVax        IFLAG=1
cVax        T=0.0
cVax      ELSE
cVax        CALL LIB$STAT_TIMER(2,ITS)
cVax        T=0.01*FLOAT(ITS)
cVax      ENDIF
cVax      RETURN
cVax      END
C
cVax      SUBROUTINE ZEDATE(BDATE)
C
C   Gets the data (DAY-MONTH-YEAR)
C   VAX version
C
cVax       CHARACTER*10 BDATE
cVax       CHARACTER*9 ADATE
cVax       CALL DATE(ADATE)
cVax       WRITE(BDATE,100) ADATE
cVax 100   FORMAT(A9,' ')
cVax       RETURN
cVax       END
C
C  ********************Vax end***********************
C
C  ********************Sun start ********************
C
       SUBROUTINE ZESEC(TBACK)
C
C   GETS CPU TIME IN SECONDS
C   Sun version
C
       REAL TARRAY(2)
       DOUBLE PRECISION TBACK
       T=ETIME(TARRAY)
       T=TARRAY(1)
       TBACK=T
       RETURN
       END
C
       SUBROUTINE ZEDATE(BDATE)
C
C   GETS THE DATE (DAY-MONTH-YEAR)
C   Sun version
C
       CHARACTER*1 BDATE(10)
       CHARACTER*1 LOCTIM(24)
       CALL FDATE(LOCTIM)
       DO 101 I = 11, 20
       II = I - 10
       BDATE(II) = LOCTIM(I)
101    CONTINUE
       RETURN
       END
C
C *****************Sun end *********************
C
C
C
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
C
c  njtj
c  ###  Cray conversions
c  ###    1)Switch double precision to real.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
C
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      INTEGER IND(M)
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
Cray      REAL D(N),E(N),E2(N),W(M),Z(NM,M),
Cray     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,
     X       NORM,ORDER,MACHEP
Cray      REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,ONEM3=1.D-3,TWO=2.D0)
Cray      PARAMETER(ZERO=0.0,ONE=1.0,ONEM3=1.E-3,TWO=2.0)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE,
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES,
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER,
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT-
C
C        ALL INPUT ARRAYS ARE UNALTERED,
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS,
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      MACHEP = TWO**(-40)
C
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = ONE - E2(1)
      Q = 0
C     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. ZERO) GO TO 140
  120 CONTINUE
C     ********** FIND VECTORS BY INVERSE ITERATION **********
  140 TAG = TAG + 1
      S = 0
C
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
C     ********** CHECK FOR ISOLATED ROOT **********
         XU = ONE
         IF (P .NE. Q) GO TO 490
         RV6(P) = ONE
         GO TO 870
  490    NORM = ABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = NORM + ABS(D(I)) + ABS(E(I))
C     ********** EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW **********
         EPS2 = ONEM3 * NORM
         EPS3 = MACHEP * NORM
         UK = REAL(Q-P+1)
         EPS4 = UK * EPS3
         UK = EPS4 / SQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
C     ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********
  510    IF (ABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. ZERO) X1 = X0 + ORDER * EPS3
C     ********** ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR **********
  520    V = ZERO
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
C     ********** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY **********
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = ZERO
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = ZERO
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. ZERO) U = EPS3
         RV1(Q) = U
         RV2(Q) = ZERO
         RV3(Q) = ZERO
C     ********** BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- **********
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP **********
         IF (GROUP .EQ. 0) GO TO 700
         J = R
C
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = ZERO
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = ZERO
C
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. ONE) GO TO 840
C     ********** FORWARD SUBSTITUTION **********
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. ZERO) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     ********** ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE **********
  780    DO 820 I = IP, Q
            U = RV6(I)
C     ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS **********
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     ********** SET ERROR -- NON-CONVERGED EIGENVECTOR **********
  830    IERR = -R
         XU = ZERO
         GO TO 870
C     ********** NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER **********
  840    U = ZERO
C
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
C
         XU = ONE / SQRT(U)
C
  870    DO 880 I = 1, N
  880    Z(I,R) = ZERO
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
      IF (Q .LT. N) GO TO 100
 1001 RETURN
C     ********** LAST CARD OF TINVIT **********
      END
C
C
C
      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
c
c  njtj
c  ###  Cray conversions
c  ###    1)Switch double precision to real.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      INTEGER IND(M)
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
Cray      REAL D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
Cray      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,PFIVE=0.5D0)
Cray      PARAMETER(ZERO=0.0,ONE=1.0,TWO=2.0,PFIVE=0.5)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
C     USING BISECTION.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY,
C
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
C          EIGENVALUES,
C
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
C
C     ON OUTPUT-
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE,
C
C        D AND E ARE UNALTERED,
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO,
C
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
C          EIGENVALUES,
C
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER,
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      MACHEP = TWO**(-40)
C
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = ZERO
C     ********** LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES **********
      DO 40 I = 1, N
         X1 = U
         U = ZERO
         IF (I .NE. N) U = ABS(E(I+1))
         XU = MIN(D(I)-(X1+U),XU)
         X0 = MAX(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         IF (ABS(E(I)) .GT. MACHEP * (ABS(D(I)) + ABS(D(I-1))))
     X      GO TO 40
   20    E2(I) = ZERO
   40 CONTINUE
C
      X1 = MAX(ABS(XU),ABS(X0)) * MACHEP * REAL(N)
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     ********** DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES **********
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
C     ********** ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS **********
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = ZERO
C
      DO 120 Q = P, N
         X1 = U
         U = ZERO
         V = ZERO
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = MIN(D(Q)-(X1+U),XU)
         X0 = MAX(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0) GO TO 140
  120 CONTINUE
C
  140 X1 = MAX(ABS(XU),ABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     ********** CHECK FOR ISOLATED ROOT WITHIN INTERVAL **********
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * REAL(Q-P+1)
      LB = MAX(T1,XU-X1)
      UB = MIN(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     ********** FIND ROOTS BY BISECTION **********
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     ********** LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) **********
      K = M2
  250    XU = LB
C     ********** FOR I=K STEP -1 UNTIL M1 DO -- **********
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     ********** NEXT BISECTION STEP **********
  300    X1 = (XU + X0) * PFIVE
         IF ((X0 - XU) .LE. (TWO * MACHEP *
     X      (ABS(XU) + ABS(X0)) + ABS(EPS1))) GO TO 420
C     ********** IN-LINE PROCEDURE FOR STURM SEQUENCE **********
  320    S = P - 1
         U = ONE
C
         DO 340 I = P, Q
            IF (U .NE. ZERO) GO TO 325
            V = ABS(E(I)) / MACHEP
            IF (E2(I) .EQ. ZERO) V = ZERO
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. ZERO) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     ********** REFINE INTERVALS **********
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     ********** K-TH EIGENVALUE FOUND **********
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     ********** ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS **********
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     ********** SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES **********
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
C     ********** LAST CARD OF TRIDIB **********
      END
C
C
C
      subroutine trnsvv(a,b,c,n)
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      dimension a(n),b(n)
c
      do 10 i=1,n
        a(i)=a(i)+c*b(i)
 10   continue
      return
      end

      subroutine velect(iter,iconv,icorr,ispp,ifcore,
     1 nr,r,rab,zel,cdd,cdu,cdc,vod,vou,etot,y,yp,
     2 ypp,s1,s2,w)
c
c    velect generates the electronic output potential from
c    the electron charge density.  The ionic part is
c    added in dsolv1/dsolv2.
c
c  njtj  ***  modifications  ***
c    The only major modiication is that the constants for the
c    ceperly-alder 'ca' method are placed in parameter
c    statements, this was done so non-opt compiliers
c    would minimize the number of calculations.
c  njtj  ***  modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter statements
c  ###      to single precision parameter statements.
c  ###  Cray conversions
c  njtj
c
       implicit double precision (a-h,o-z)
c
       character*1 ispp
       character*2 icorr
c
c  njtj  *** modification start  ***
c
       parameter (zero=0.D0,one=1.D0,pfive=.5D0,opf=1.5D0,pnn=.99D0)
       parameter (pthree=0.3D0,psevf=0.75D0,c0504=0.0504D0)
       parameter (c0254=0.0254D0,c014=0.014D0,c0406=0.0406D0)
       parameter (c15p9=15.9D0,c0666=0.0666D0,c11p4=11.4D0)
       parameter (c045=0.045D0,c7p8=7.8D0,c88=0.88D0,c20p592=20.592D0)
       parameter (c3p52=3.52D0,c0311=0.0311D0,c0014=0.0014D0)
       parameter (c0538=0.0538D0,c0096=0.0096D0,c096=0.096D0)
       parameter (c0622=0.0622D0,c004=0.004D0,c0232=0.0232D0)
       parameter (c1686=0.1686D0,c1p3981=1.3981D0,c2611=0.2611D0)
       parameter (c2846=0.2846D0,c1p0529=1.0529D0,c3334=0.3334D0)
Cray       parameter (zero=0.0,one=1.0,pfive=0.5,opf=1.5,pnn=0.99)
Cray       parameter (pthree=0.3,psevf=0.75,c0504=0.0504)
Cray       parameter (c0254=0.0254,c014=0.014,c0406=0.0406)
Cray       parameter (c15p9=15.9,c0666=0.0666,c11p4=11.4)
Cray       parameter (c045=0.045,c7p8=7.8,c88=0.88,c20p592=20.592)
Cray       parameter (c3p52=3.52,c0311=0.0311,c0014=0.0014)
Cray       parameter (c0538=0.0538,c0096=0.0096,c096=0.096)
Cray       parameter (c0622=0.0622,c004=0.004,c0232=0.0232)
Cray       parameter (c1686=0.1686,c1p3981=1.3981,c2611=0.2611)
Cray       parameter (c2846=0.2846,c1p0529=1.0529,c3334=0.3334)
c
c    Ceperly-Alder 'ca' constants
c
       parameter (con1=1.D0/6, con2=0.008D0/3, con3=0.3502D0/3)
       parameter (con4=0.0504D0/3, con5=0.0028D0/3, con6=0.1925D0/3)
       parameter (con7=0.0206D0/3, con8=9.7867D0/6, con9=1.0444D0/3)
       parameter (con10=7.3703D0/6, con11=1.3336D0/3)
Cray       parameter (con1=1.0/6, con2=0.008/3, con3=0.3502/3)
Cray       parameter (con4=0.0504/3, con5=0.0028/3, con6=0.1925/3)
Cray       parameter (con7=0.0206/3, con8=9.7867/6, con9=1.0444/3)
Cray       parameter (con10=7.3703/6, con11=1.3336/3)
c
c  njtj  ***  modification end  ***
c
      dimension r(nr),rab(nr),cdd(nr),cdu(nr),cdc(nr),
     1 vod(nr),vou(nr),etot(10),y(nr),yp(nr),ypp(nr),
     2 s1(nr),s2(nr),w(3*nr)
c
       pi=4*atan(one)
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c
Csun      expzer = 3.7D2
cApollo      expzer = 3.7D2
       expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
c      fit cd/r by splines
c
       y(1) = zero
       do 10 i=2,nr
         y(i) = (cdd(i)+cdu(i))/r(i)
 10    continue
       if (ifcore .eq. 2) then
         do 11 i=2,nr
           y(i) = y(i) + cdc(i)/r(i)
 11      continue
       endif
       isx = 0
       a1 = zero
       an = zero
       b1 = zero
       bn = zero
       nrm=nr
       call splift(r,y,yp,ypp,nrm,w,ierr,isx,a1,b1,an,bn)
       if(ierr.ne.1) then
         write(6,20000)ierr
         call ext(420+ierr)
       endif
20000  format(1x,'****** Error in splift ierr =',i2)
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
       xlo = zero
       call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s2,ierr)
       if(ierr.ne.1) then
         write(6,20001)ierr
         call ext(440+ierr)
       endif
20001  format(1x,'****** Error in spliq ierr =',i2)
       do 20 i=1,nr
         ypp(i) = r(i)*ypp(i) + 2*yp(i)
         yp(i)  = r(i)*yp(i)  + y(i)
         y(i)   = r(i)*y(i)
 20    continue
       call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s1,ierr)
       if(ierr.ne.1) then
         write(6,20002)ierr
         call ext(460+ierr)
       endif
20002  format(1x,'****** Error in spliq ierr =',i2)
c
c      check normalization
c
       xnorm = zero
       if (ifcore .eq. 2 .and. iter .eq. 0 ) zel=s1(nr)
       if (zel .ne. zero) xnorm = zel/s1(nr)
       if (iter .gt. 3 .and. abs(zel-s1(nr)) .gt. 0.01) then
         if (zel .lt. s1(nr)+1.0 ) then
           write(6,24) iter,xnorm
 24    format(/,' warning *** charge density rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         else
           xnorm=pnn*xnorm
           write(6,25) iter,xnorm
 25    format(/,' warning *** charge density partially rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         endif
       endif
c
c      compute new hartree potential
c      renormalize the charge density
c
       do 30 i=2,nr
         vod(i) = 2 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
         vou(i) = vod(i)
         cdd(i) = xnorm*cdd(i)
         cdu(i) = xnorm*cdu(i)
 30    continue
c
c      compute hartree contribution to total energy
c
       if (iconv .eq. 1) then
         ehart = zero
         ll = 4
         do 40 i=2,nr
           ehart = ehart+ll*(cdd(i)+cdu(i))*vod(i)*rab(i)
           ll = 6 - ll
 40      continue
         ehart = ehart / 6
       endif
c
c      add exchange and correlation
c
       trd = one/3
       ftrd = 4*trd
       tftm = 2**ftrd-2
       a0 = (4/(9*pi))**trd
c
c      set x-alpha
c
       alp = one
       if (icorr .ne. 'xa') alp = 2 * trd
       vxc = zero
       vc  = zero
       exc = zero
       ec  = zero
c
c      start loop
c
       ll = 4
       do 210 i=2,nr
         cdsum = cdd(i) + cdu(i)
         if (ifcore .ge. 1) cdsum=cdsum+cdc(i)
         if (cdsum .le. zero) goto 210
c
c  Vax bug fix.  Troy Barbee - 4/17/90
c
         if (log(3*r(i)**2/cdsum) .gt. 2*expzer) goto 210
         rs = (3*r(i)**2/cdsum)**trd
         z = zero
         fz = zero
         fzp = zero
         if (ispp .eq. 's') then
           z = (cdd(i)-cdu(i)) / cdsum
           fz = ((1+z)**ftrd+(1-z)**ftrd-2)/tftm
           fzp = ftrd*((1+z)**trd-(1-z)**trd)/tftm
         endif
c
c      exchange (only use (xa))
c
         vxp = -3*alp/(pi*a0*rs)
         exp = 3*vxp/4
         if (ispp .eq. 'r') then
           beta = c014/rs
           sb = sqrt(1+beta*beta)
           alb = log(beta+sb)
           vxp = vxp * (-pfive + opf * alb / (beta*sb))
           exp = exp *(one-opf*((beta*sb-alb)/beta**2)**2)
         endif
 65      vxf = 2**trd*vxp
         exf = 2**trd*exp
         vcp = zero
         ecp = zero
         vcf = zero
         ecf = zero
         if (icorr .eq. 'ca') then
c          ceperly-alder (ca)
c          The Perdew-Zunger parameterization is used.
c          See Phys. Rev. B 23 5075 (1981).
           if (rs .gt. one) then
             sqrs=sqrt(rs)
             te = one+con10*sqrs+con11*rs
             be = one+c1p0529*sqrs+c3334*rs
             ecp = -c2846/be
             vcp = ecp*te/be
             te = one+con8*sqrs+con9*rs
             be = one+c1p3981*sqrs+c2611*rs
             ecf = -c1686/be
             vcf = ecf*te/be
           else
             rslog=log(rs)
             ecp=(c0622+c004*rs)*rslog-c096-c0232*rs
             vcp=(c0622+con2*rs)*rslog-con3-con4*rs
             ecf=(c0311+c0014*rs)*rslog-c0538-c0096*rs
             vcf=(c0311+con5*rs)*rslog-con6-con7*rs
           endif
         elseif (icorr .eq. 'xa') then
c          correlation
         elseif (icorr .eq. 'wi') then
c          wigner (wi)
           vcp = -(c3p52*rs+c20p592)/(3*(rs+c7p8)**2)
           ecp = -c88/(rs+c7p8)
         elseif (icorr .eq. 'hl') then
c          hedin-lundqvist (hl)
           x = rs/21
           aln = log(1+1/x)
           vcp = -c045*aln
           ecp = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecp=((con1/x-pthree)/x+psevf)/x
           ecp = -c045*ecp
         elseif (icorr .eq. 'gl') then
c          gunnarson-lundqvist-wilkins (gl)
           x = rs/c11p4
           aln = log(1+1/x)
           vcp = -c0666*aln
           ecp = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecp=((con1/x-pthree)/x+psevf)/x
           ecp = -c0666*ecp
           x = rs/c15p9
           aln = log(1+1/x)
           vcf = -c0406*aln
           ecf = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecf=((con1/x-pthree)/x+psevf)/x
           ecf = -c0406*ecf
         elseif (icorr .eq. 'bh') then
c          von barth - hedin (bh)
           x = rs/30
           aln = log(1+1/x)
           vcp = -c0504*aln
           ecp = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecp=((con1/x-pthree)/x+psevf)/x
           ecp = -c0504*ecp
           x = rs/75
           aln = log(1+1/x)
           vcf = -c0254*aln
           ecf = aln+(x**3*aln-x*x)+x/2-trd
           if (x .gt. 500*one) ecf=((con1/x-pthree)/x+psevf)/x
           ecf = -c0254*ecf
         else
           write(6,70) icorr
           call ext(400)
         endif
 70   format('error in velect - icorr =',a2,' not implemented')
         vxcp = vxp + vcp
         vxcf = vxf + vcf
         vxcd = vxcp
         vxcu = vxcp
         excp = exp + ecp
         excf = exf + ecf
         vcd = vcp
         vcu = vcp
         exct = excp
         ect = ecp
         if (z .ne. zero) then
           vxcd = vxcd + fz*(vxcf-vxcp) + (1-z)*fzp*(excf-excp)
           vxcu = vxcu + fz*(vxcf-vxcp) - (1+z)*fzp*(excf-excp)
           vcd = vcd + fz*(vcf-vcp) + (1-z)*fzp*(ecf-ecp)
           vcu = vcu + fz*(vcf-vcp) - (1+z)*fzp*(ecf-ecp)
           exct = exct + fz*(excf-excp)
           ect = ect + fz*(ecf-ecp)
         endif
         vod(i) = vod(i) + vxcd
         vou(i) = vou(i) + vxcu
         vxc = vxc + ll * (cdd(i)*vxcd + cdu(i)*vxcu) * rab(i)
         vc  = vc  + ll * (cdd(i)*vcd  + cdu(i)*vcu ) * rab(i)
         exc = exc + ll * cdsum * exct * rab(i)
         ec  = ec  + ll * cdsum * ect  * rab(i)
         ll = 6 - ll
 210   continue
       etot(4) = ehart
       etot(5) = vxc / 3
       etot(6) = (3*vc - 4*ec) / 3
       etot(7) = exc / 3
       vod(1) = vod(2) - (vod(3)-vod(2))*r(2)/(r(3)-r(2))
       vou(1) = vou(2) - (vou(3)-vou(2))*r(2)/(r(3)-r(2))
       return
       end
C
C
C
      subroutine vionic(ispp,itype,icorr,ifcore,zsh,rsh,
     1 lmax,nr,a,b,r,rab,nameat,ncore,znuc,
     2 cdd,cdu,cdc,viod,viou)
c
c  Vionic sets up the ionic potential.
c  Note that viod/viou is the ionic potential times r.
c
c  njtj ***  major modifications  ***
c    If a potential does not exist, it is approximated
c    by an existing potential.
c    A nonspin or spin-polarized pseudo test, uses the
c    down(nonspin generation), weighted average(spin-
c    polarized), or averaged(relativistic) potentials.
c    A relativistic pseudo test, must use relativistic
c    generated potentials.  The Schroedinger equation is
c    used to integrate a relativistic pseudo test,
c    not the Dirac equation.
c  njtj  ***  major modifications  ***
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0)
Cray      parameter (zero=0.0)
c
      character*1 ispp
      character*2 icorr,icorrt,nameat,namet
      character*3 irel
      character*4 nicore
      character*10 iray(6),ititle(7)

      dimension r(nr),rab(nr),cdd(nr),cdu(nr),cdc(nr),
     1 viod(lmax,nr),viou(lmax,nr),npd(5),npu(5)
c
c  2*znuc part
c
      ifcore = 0
      if (itype .lt. 4) then
        do 10 i=1,lmax
          do 12 j=1,nr
            viod(i,j) = -2*znuc
            viou(i,j) = -2*znuc
 12       continue
 10     continue
      else
c
c  read pseudopotentials from tape1
c
        rewind 1
        read(1) namet,icorrt,irel,nicore,(iray(i),i=1,6),
     1   (ititle(i),i=1,7),npotd,npotu,nrm,a,b,zion
        if(nicore.eq.'fcec'.or.nicore.eq.'pcec') ifcore = 1
        if(nicore.eq.'fche'.or.nicore.eq.'pche') ifcore = 2
        nr = nrm+1
        read(1) (r(i),i=2,nr)
        r(1) = zero
c
c   down potentials (or average relativistic potentials)
c
c njtj  ***  major start  ***
c   if a potential does not exist, it is replaced by the
c   next existing lower angular momentum potential or
c   the next existing higher if no lower exist.
c
        do 15 i=1,lmax
          npd(i)=0
 15     continue
        do 20 i=1,npotd
          read(1) loi,(viod(loi+1,j),j=2,nr)
          viod(loi+1,1) = zero
          npd(loi+1)=1
 20     continue
        if (npd(1) .eq. 0) then
          do 25 i=2,lmax
            if (npd(i) .gt. 0) then
              do 24 j=1,nr
                viod(1,j)=viod(i,j)
 24           continue
              goto 30
            endif
 25       continue
        endif
 30     do 33 i=2,lmax
          if (npd(i) .eq. 0) then
            do 32 j=1,nr
              viod(i,j)=viod(i-1,j)
 32         continue
          endif
 33     continue
c
c   up potentials (or spin orbit potentials)
c
        if (npotu .le. 0) goto 49
        do 35 i=1,lmax
          npu(i)=0
 35     continue
        do 37 i=1,npotu
          read(1) loi,(viou(loi+1,j),j=2,nr)
          viou(loi+1,1) = zero
          npu(loi+1)=1
 37     continue
        if (npu(1) .eq. 0) then
          do 38 i=2,lmax
            if (npu(i) .gt. 0) then
              do 39 j=1,nr
                viou(1,j)=viou(i,j)
 39           continue
              goto 40
            endif
 38       continue
        endif
 40     do 45 i=2,lmax
          if (npu(i) .eq. 0) then
            do 43 j=1,nr
              viou(i,j)=viou(i-1,j)
 43         continue
          endif
 45     continue
c
c  njtj  ***  major end  ***
c
c
c  core and valence charges
c
 49     read(1) (cdc(i),i=2,nr)
        cdc(1) = zero
c
c  replace valence charge on tape(valence charge modify)
c
        if (itype .eq. 6) then
          write(1) (cdd(i)+cdu(i),i=2,nr)
          return
        endif
        read(1) (cdd(i),i=2,nr)
        cdd(1) = zero
c
c  njtj  ***   major start  ***
c   distribute charge as up and down charge
c   generate radial intergration grid
c   set up potentials equal to down potentials for
c   spin-polarized pseudo test of nonspin and relativistic
c   generated potentails.  Construct spin-orbit potentials
c   from relativistic sum and difference potentials and
c   change ispp='r' to ispp=' '.
c
        do 50 i=1,nr
          rab(i) = (r(i)+a)*b
          cdd(i) = cdd(i)/2
          cdu(i) = cdd(i)
 50     continue
        if (ispp .eq. 's' .and. irel .ne. 'isp') then
          do 51 i=1,lmax
            do 52 j=1,nr
              viou(i,j) = viod(i,j)
 52         continue
 51       continue
        endif
        if (ispp .eq. 'r') then
          ispp=' '
          if (irel .ne. 'rel') then
            write(6,130)irel
 130  format(//,'Pseudopotentail is not relativistic!!!!',/
     1 ' setting up potentials equal to down!!!',//)
            do 53 i=1,lmax
              do 54 j=1,nr
                viou(i,j) = viod(i,j)
 54           continue
 53         continue
          else
            do 57 j=1,nr
              viou(1,j)=viod(1,j)
 57         continue
            do 58 i=2,lmax
              do 56 j=1,nr
                vsum=viod(i,j)
                vdiff=viou(i,j)
                viod(i,j)=vsum-i*vdiff/2
                viou(i,j)=vsum+(i-1)*vdiff/2
 56           continue
 58         continue
          endif
        endif
c
c   njtj  ***  major end   ***
c
c
c   printout
c
        write(6,60) namet,icorrt,irel,nicore,(iray(i),i=1,6),
     1   (ititle(i),i=1,7)
 60   format(//,1x,a2,2x,a2,2x,a3,2x,a4,
     1 '  pseudopotential read from tape',
     2 /,1x,2a10,5x,4a10,/,1x,7a10,//)
        if (nameat .ne. namet) write(6,70) nameat,namet
 70   format(' input element ',a2,
     1 ' not equal to element on tape ',a2,//)
        if (icorr .ne. icorrt) write(6,80) icorr,icorrt
 80   format(' input correlation ',a2,
     1 ' not equal to correlation from tape ',a2,//)
        write(6,90) r(2),nr,r(nr)
 90   format(' radial grid parameters',//,
     1 ' r(1) = .0 , r(2) =',e8.2,' , ... , r(',i3,') =',
     2 f6.2,//)
      endif
c
c   add potential from shell charge
c
      if (abs(zsh) .gt. 0.e-5) then
        do 110 i=1,lmax
          do 120 j=1,nr
            if (r(j) .ge. rsh) then
              viod(i,j) = viod(i,j) - 2*zsh
              viou(i,j) = viou(i,j) - 2*zsh
            else
              viod(i,j) = viod(i,j) - 2*zsh*r(i)/rsh
              viou(i,j) = viou(i,j) - 2*zsh*r(i)/rsh
            endif
 120      continue
 110    continue
       endif
       return
       end
C
C
C
      subroutine wtrans(vd,r,nr,rab,l,ist,b)
c
c **********************************************************
c *
c *    This is a plotting routine; the user should adjust
c *  for their own needs.  The result
c *  is then printed to the current plot.dat file (unit=3)
c *  for later plotting of the data.  A marker (marker fw#)
c *  is placed at the end of each set of data.
c *
c **********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.D0,one=1.D0,big=17280.0D0,p5=.05D0)
Cray      parameter (zero=0.0,one=1.0,big=17280.0,p5=.05)
c
      dimension vd(nr),r(nr),rab(nr),b(nr),vql(48),vql2(48),
     1 a(2000),vdpp(2000),r2(2000),v(2000),w(4000)
c
      do 1 i=1,48
        vql(i)=zero
 1    continue
c
c  The wavefuncion(rR) times r times rab.
c
      if (abs(ist) .eq. 2) goto 400
      pi4=16*atan(one)
      do 10 k=2,nr
        if (r(k)-r(k-1).gt. p5) then
          nr2=k
          goto 20
        endif
 10   continue
 20   nr2=7*(nr2/7)+1
      nr3=nr2-7
      do 130 k=2,nr2
        b(k)=vd(k)*r(k)*rab(k)
 130  continue
      do 150 k=nr2,nr
        a(k-nr2+1)=vd(k)*r(k)
 150  continue
      isx = 0
      a1 = -p5*10
      an = -p5*10
      b1 = zero
      bn = zero
      nrm=nr-nr2+1
      call splift(r(nr2),a,r2,vdpp,nrm,w,ierr,isx,a1,b1,an,bn)
      if(ierr.ne.1) then
        call exit
      endif
      nr4=0
      do 155 ak=r(nr2),100.0D0,0.05D0
        nr4=nr4+1
        r2(nr4)=ak
 155  continue
      call splint(r(nr2),a,vdpp,nrm,r2,v,w,w(2000),nr4,kerr)
c
c  Find the fourier transform-vql.
c
      do 140 j=1,48
        q=one/4*j
        vql(j)=zero
        a(1)=zero
        do 135 k=2,nr2
          a(k)=b(k)*sbessj(l,q*r(k))
 135    continue
c
c  Due to the high number of occilations in the intagrand,
c  an eight point Newton-Cotes intagration method is used.
c  See  Abramowitz and Stegun Eq. 25.4.17
c
        do 145 k=1,nr3,7
          vql(j)=vql(j)+751*(a(k)+a(k+7))+3577*(a(k+1)+a(k+6))+
     1     1323*(a(k+2)+a(k+5))+2989*(a(k+3)+a(k+4))
 145    continue
        vql(j)=pi4*7*vql(j)/big
        do 160 k=1,nr4
          a(k)=v(k)*sbessj(l,q*r2(k))
 160    continue
        vql2(j)=zero
        do 165 kk=8,nr4,7
          k=kk-7
          vql2(j)=vql2(j)+751*(a(k)+a(k+7))+3577*(a(k+1)+
     1     a(k+6))+1323*(a(k+2)+a(k+5))+2989*(a(k+3)+a(k+4))
 165    continue
        vql2(j)=0.35D0*pi4*vql2(j)/big
        vql(j)=vql(j)+vql2(j)
 140  continue
c
c  Print out the transform vql(q) to the current plot.dat
c  file (unit=3) for latter plotting.
c
 400  do 170 j=1,48
        write(3,6000)one/4*j,ist*vql(j)
 170  continue
      write(3,6001)l
      return
c
c  format statements
c
 6000 format(1x,f7.4,3x,f10.6)
 6001 format(1x,'marker fw',i1)
      end
C
C
C
      subroutine zrbac2(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1 rc8,lp,arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2 gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c **********************************************************
c *  njtj
c *    Routine brackets the root of the given function.
c *    Taken from Numerical Recipes page 245.
c *  njtj
c **********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (factor=1.6D0,ntry=50)
Cray      parameter (factor=1.6,ntry=50)
c
      dimension r(jrc),rab(jrc),ar(jrc)
c
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f1,ar)
      call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f2,ar)
c
      do 11 j=1,ntry
        if(f1*f2.lt.0.0)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+factor*(x1-x2)
          call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1     arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2     alpha,alpha1,alpha2,alpha3,alpha4,f1,ar)
        else
          x2=x2+factor*(x2-x1)
          call gamfn2(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1     arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2     alpha,alpha1,alpha2,alpha3,alpha4,f2,ar)
        endif
11    continue
c
c  failure, abort program
c
      write(6,1000)lp
      call ext(830+lp)
 1000 format(//,'error in zbractk - can not bracket orbital ',i2)
      return
      end
C
C
C
      subroutine zrbact(x1,x2,rc1,rc2,rc3,rc4,rc5,rc6,rc7,
     1 rc8,lp,arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,
     2 gamma,alpha,alpha1,alpha2,alpha3,alpha4,ar)
c
c **********************************************************
c *  njtj
c *    Routine brackets the root of the given function.
c *    Taken from Numerical Recipes page 245.
c *  njtj
c **********************************************************
c
c  njtj
c  ###  Cray conversions
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
      implicit double precision (a-h,o-z)
c
      parameter (factor=1.6D0,ntry=50)
Cray      parameter (factor=1.6,ntry=50)
c
      dimension r(jrc),rab(jrc),ar(jrc)
c
      call gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f1,ar)
      call gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1 arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2 alpha,alpha1,alpha2,alpha3,alpha4,f2,ar)
c
      do 11 j=1,ntry
        if(f1*f2.lt.0.0)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+factor*(x1-x2)
          call gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1     arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x1,
     2     alpha,alpha1,alpha2,alpha3,alpha4,f1,ar)
        else
          x2=x2+factor*(x2-x1)
          call gamfnd(rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,lp,
     1     arc,brc,vrc,vap,vapp,ev,cdrc,r,rab,jrc,delta,x2,
     2     alpha,alpha1,alpha2,alpha3,alpha4,f2,ar)
        endif
11    continue
c
c  failure, abort program
c
      write(6,1000)lp
      call ext(830+lp)
 1000 format(//,'error in zbractk - can not bracket orbital ',i2)
      return
      end
C
C
C
c $Id$
