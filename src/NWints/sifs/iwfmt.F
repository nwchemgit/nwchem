*iwfmt.f
*iwfmt part=1 of 1.  formatted integral file utility program
*version=1.0 last modified: 04-may-92
*version 5.0
c
      program iwfmt
c
c  read an integral file and write a formatted version.
c
c  version 1.0b5 04-may-92
c
c  version history:
c  04-may-92 minor ftnchek-related cleanup. -rls
c  21-oct-90 written by ron shepard, based on program istat.
c  istat history:
c  08-oct-90 (columbus day) 1-e fcore change. -rls
c  06-oct-90 fsplit=2 capability added. -rls
c  11-aug-89 siftyp() call added. -rls
c  08-aug-89 SIFS version. -rls
c  written by ron shepard.
c
      implicit none
c
      integer aoints, buf, i, ibitv, idummy, ierr, ilab, itotal,
     & lenbuf, nbft, nenrgy, ninfo, nlist, nmap, nmax, nsym, ntitle,
     & val, nin
      character*80 title(20)
      integer nbpsy(8)
c
      integer   nbfmxp,     nengmx,    nmapmx
      parameter(nbfmxp=255, nengmx=20, nmapmx=10)
c
      integer ietype(nengmx), info(10)
      integer imtype(nmapmx), map(nbfmxp*nmapmx)
      real*8 energy(nengmx)
      character*4 slabel(8)
      character*8 bfnlab(nbfmxp)
      character*60 fname
c
      integer   lencor
      parameter(lencor=200000)
      real*8 core(lencor)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
#if !defined(NONAMELIST)
      namelist /input/ fname
#endif 
c
      integer  forbyt, atebyt
      external forbyt, atebyt
c
      nlist = 6
      nin = 5
      open( unit=nin,   file='iwfmtin',  status='unknown'     )
      aoints= 10
c
#if defined(FUTURE)
c  open nlist here if necessary.
#else 
c  assume preconnected unit.
#endif 
c
      call ibummr(nlist)
c
      fname = 'aoints'
      read( nin, input, end=3 )
3     continue
      call trnfln( 1, fname )
      open(unit=aoints,file=fname,form='unformatted',status='old')
c
c     # read the header info.
c
      call sifrh1( aoints, ntitle, nsym, nbft,
     & ninfo, nenrgy, nmap, ierr )
      if ( ierr .ne. 0 ) then
         call bummer('iwfmt: ierr=',ierr,faterr)
      elseif ( ntitle .gt. 20 ) then
         call bummer('iwfmt: ntitle=',ntitle,faterr)
      elseif ( nbft .gt. nbfmxp ) then
         call bummer('iwfmt: nbft=',nbft,faterr)
      elseif ( ninfo .gt. 10 ) then
         call bummer('iwfmt: ninfo=',ninfo,faterr)
      elseif ( nenrgy .gt. nengmx ) then
         call bummer('iwfmt: nenrgy=',nenrgy,faterr)
      elseif( nmap .gt. nmapmx ) then
         call bummer('iwfmt: nmap=',nmap,faterr)
      endif
c
      write(nlist,6100) ntitle, nsym, nbft, ninfo, nenrgy, nmap
6100  format(1x,10i8)
c
      call sifrh2(
     & aoints, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   title,
     & nbpsy,  slabel, info,   bfnlab,
     & ietype, energy, imtype, map,
     & ierr )
      if ( ierr.ne.0 ) call bummer('iwfmt: ierr=',ierr,faterr)
c
      write(nlist,6020) (title(i),i=1,ntitle)
6020  format(1x,a)
c
      write(nlist,6100) (nbpsy(i),i=1,nsym)
c
      write(nlist,6030) (slabel(i),i=1,nsym)
6030  format(1x,8a5)
c
      write(nlist,6100) (info(i),i=1,ninfo)
c
      write(nlist,6040) (bfnlab(i),i=1,nbft)
6040  format(1x,8a10)
c
      write(nlist,6100) (ietype(i),i=1,nenrgy)
c
      write(nlist,6050) (energy(i),i=1,nenrgy)
6050  format(1x,1p4e20.12)
c
      if ( nmap .gt. 0 ) then
         write(nlist,6100) ( imtype(i), i = 1, nmap )
         write(nlist,6060) ( map(i), i = 1, nmap*nbft )
      endif
6060  format(1x,20i4)
c
      if ( (info(1) .ne. 1) .and. (info(1) .ne. 2) ) then
         call bummer('iwfmt: fsplit=',info(1),faterr)
      endif
c
c     # print the 1-e integrals.
c     # core=>1:buffer,2:ilab1,3:val1,4:ibitv
c
      lenbuf = info(2)
      nmax   = info(3)
c
      buf    = 1
      ilab   = buf   + atebyt( lenbuf )
      val    = ilab  + forbyt( 2*nmax )
      ibitv  = val   + atebyt( nmax )
      itotal = ibitv + forbyt( ((nmax+63)/64)*64 ) -1
      if ( itotal .gt. lencor ) then
         call bummer( 'iwfmt: prt1e itotal=',itotal,faterr)
      endif
c
      call prt1e(
     & nlist,      aoints,    info,        core(buf),
     & core(ilab), core(val), core(ibitv) )
c
c     # open the 2-e file.
      fname = 'aoints2'
      call trnfln( 1, fname )
      call sifo2f( aoints, aoints, fname, info, idummy, ierr )
      if ( ierr .ne. 0 ) then
         call bummer('iwfmt: from sifo2f, ierr=',ierr,faterr)
      endif
c
c     # print the 2-e integrals.
c     # core=>1:buffer,2:ilab1,3:val1,4:ibitv
c
      lenbuf = info(4)
      nmax   = info(5)
c
      buf    = 1
      ilab   = buf   + atebyt( lenbuf )
      val    = ilab  + forbyt( 4*nmax )
      ibitv  = val   + atebyt( nmax )
      itotal = ibitv + forbyt( ((nmax+63)/64)*64 ) -1
      if ( itotal .gt. lencor ) then
         call bummer( 'iwfmt: prt2e itotal=',itotal,faterr)
      endif
c
      call prt2e(
     & nlist,      aoints,    info,        core(buf),
     & core(ilab), core(val), core(ibitv) )
c
      call sifc2f( aoints, info, ierr )
      if ( ierr .ne. 0 ) then
         call bummer('iwfmt: from sifc2f, ierr=',ierr,faterr)
      endif
      close (unit = aoints)
c
      call bummer('normal termination',0,3)
      stop
      end
*deck prt1e
      subroutine prt1e(
     & nlist,  ntape,  info,   buf,
     & ilab,   val,    ibitv )
c
c  read and dump the 1-e integral records.
c
c  21-oct-90 written by ron shepard.
c
      implicit none
c
      integer   nipv,   msame,   nmsame,   nomore
      parameter(nipv=2, msame=0, nmsame=1, nomore=2)
c
      integer nlist, ntape
      real*8 buf(*), val(*)
      integer info(*), ilab(nipv,*), ibitv(*)
c
      integer ibvtyp, ierr, ifmt, itypea, itypeb, last, num
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      integer   iretbv
      parameter(iretbv=-1)
c
      real*8 fcore
c
      last = msame
200   if ( last .ne. nomore ) then
c
         call sifrd1( ntape, info, nipv, iretbv,
     &    buf, num, last, itypea,
     &    itypeb, ifmt, ibvtyp, val,
     &    ilab, fcore, ibitv, ierr )
c
         if ( ierr.ne.0 ) call bummer(
     &    'iwfmt: sifrd1 ierr=',ierr,faterr)
c
         call siffw1(
     &    info,   nipv,   num,    last,
     &    itypea, itypeb, ifmt,   ibvtyp,
     &    val,    ilab,   fcore,  ibitv,
     &    nlist,  ierr  )
c
         if ( ierr .ne. 0 ) call bummer(
     &    'iwfmt: siffw1 ierr=',ierr,faterr)
c
         goto 200
      endif
c
      return
      end
*deck prt2e
      subroutine prt2e(
     & nlist,  ntape,  info,   buf,
     & ilab,   val,    ibitv )
c
c  read and dump the 2-e integral records.
c
c  21-oct-90 written by ron shepard.
c
      implicit none
c
      integer   nipv,   msame,   nmsame,   nomore
      parameter(nipv=4, msame=0, nmsame=1, nomore=2)
c
      integer nlist, ntape
      real*8 buf(*), val(*)
      integer info(*), ilab(nipv,*), ibitv(*)
c
      integer ibvtyp, ierr, ifmt, itypea, itypeb, last, num
      integer   iretbv
      parameter(iretbv=-1)
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      last = msame
200   if ( last .ne. nomore ) then
c
         call sifrd2(
     &    ntape,  info,   nipv,   iretbv,
     &    buf,    num,    last,   itypea,
     &    itypeb, ifmt,   ibvtyp, val,
     &    ilab,   ibitv,  ierr )
c
         if ( ierr.ne.0 ) call bummer(
     &    'iwfmt: sifrd2 ierr=',ierr,faterr)
c
         call siffw2(
     &    info, nipv, num, last,
     &    itypea, itypeb, ifmt, ibvtyp,
     &    val, ilab, ibitv, nlist,
     &    ierr )
c
         if ( ierr.ne.0 ) call bummer(
     &    'iwfmt: siffw2 ierr=',ierr,faterr)
c
         goto 200
      endif
c
      return
      end
