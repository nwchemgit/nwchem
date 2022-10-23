      subroutine sifwh(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    title,
     & nbpsy,   slabel,  info,    bfnlab,
     & ietype,  energy,  imtype,  map,
     & ierr )
c
c  write the header records of the standard integral file structure.
c
c  input:
c  aoints = output file unit number.
c  ntitle = number of titles.
c  nsym = number of symmetry blocks.
c  nbft = total number of basis functions.
c  ninfo = number of record-definition parameters.
c  nenrgy = number of core energies.  the first element must be the
c           nuclear repulsion energy.
c  nmap = number of optional map vectors.
c  title*80(1:ntitle) = identifying titles.
c  nbpsy(1:nsym) = number of basis functions per symmetry block.
c  slabel*4(1:nsym) = symmetry labels.
c  info(1:ninfo) = record-definition parameters.
c  bfnlab*8(1:nbft) = basis function labels.
c  ietype(1:nenrgy) = core energy types.
c  energy(1:nenrgy) = core energies.
c  imtype(1:nmap) = map vector types.
c  map(1:nbft,1:nmap) = basis function map vectors.
c
c  output:
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89  written by ron shepard.
c
      implicit none
c
c  vrsion = routine library version number.
c  ntitmx = maximum number of titles allowed.
c  ninchk = minimum number of info(*) elements.
c  lrecmx = maximum record length allowed.  this should be consistent
c           with dword bit-packing in the record-writing routines.
c
      integer   vrsion,   ntitmx,    ninchk
      parameter(vrsion=1, ntitmx=20, ninchk=5)

      integer   lrecmx
      parameter(lrecmx=2**16-1)
c
c     # using (*) dimensions until parameters are checked...
      integer  aoints, ntitle, nsym,   nbft,   ninfo,
     & nenrgy, nmap,   ierr
      character*80     title(*)
      integer  nbpsy(*)
      character*4      slabel(*)
      integer  info(*)
      character*8      bfnlab(*)
      integer  ietype(*)
      real*8   energy(*)
      integer  imtype(*)
      integer  map(*)
c
      integer  isym,   mapdim, nbftx
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
c     # check arguments for validity...
c
cgk debug
      write(6,*)'gk: in sifwh'
      write(*,*)'gk: in sifwh'
cgk end
      ierr = 0
      if ( ntitle.le.0 .or. ntitle.gt.ntitmx ) then
         call bummer('sifwh: ntitle=',ntitle,wrnerr)
         ierr = -1
         return
      elseif ( nsym.ne.1 .and. nsym.ne.2 .and. nsym.ne.4
     &    .and. nsym.ne.8 ) then
         call bummer('sifwh: nsym=',nsym,wrnerr)
         ierr = -1
         return
      elseif ( nbft.le.0 ) then
         call bummer('sifwh: nbft=',nbft,wrnerr)
         ierr = -1
         return
      elseif ( ninfo.lt.ninchk ) then
         call bummer('sifwh: ninfo=',ninfo,wrnerr)
         ierr = -1
         return
      elseif ( info(1).ne.1 .and. info(1).ne.2  ) then
c        # fsplit parameter.
         call bummer('sifwh: fsplit=',info(1),wrnerr)
         ierr = -1
         return
      elseif ( info(2).le.0 .or. info(2).gt.lrecmx ) then
c        # l1rec.
         call bummer('sifwh: l1rec=',info(2),wrnerr)
         ierr = -1
         return
      elseif ( info(3).le.0 ) then
c        # n1max.
         call bummer('sifwh: n1max=',info(3),wrnerr)
         ierr = -1
         return
      elseif ( info(4).le.0 .or. info(4).gt.lrecmx ) then
c        # l2rec.
         call bummer('sifwh: l2rec=',info(4),wrnerr)
         ierr = -1
         return
      elseif ( info(5).le.0 ) then
c        # n2max.
         call bummer('sifwh: n2max=',info(5),wrnerr)
         ierr = -1
         return
      elseif ( nenrgy.le.0 ) then
         call bummer('sifwh: nenrgy=',nenrgy,wrnerr)
         ierr = -1
         return
      elseif ( ietype(1).ne.-1 ) then
         call bummer('sifwh: ietype(1)=',ietype(1),wrnerr)
         ierr = -1
         return
      elseif ( nmap.lt.0 ) then
         call bummer('sifwh: nmap=',nmap,wrnerr)
         ierr = -1
         return
      endif
c
      nbftx=0
      do 10 isym=1,nsym
         if ( nbpsy(isym) .lt. 0 ) then
            call bummer('sifwh: nbpsy(isym)=',nbpsy(isym),wrnerr)
            ierr = -1
            return
         endif
         nbftx=nbftx+nbpsy(isym)
10    continue
      if ( nbftx.ne.nbft ) then
         call bummer('sifwh: nbftx=',nbftx,wrnerr)
         ierr = -1
         return
      endif
c
      mapdim = max( nmap, 1 )
cgk debug
      write(*,*)'gk: in sifwh calling sifzwh'
      write(*,*)'gk:    ntitle= ', ntitle
cgk end
      call sifzwh(
     & aoints, vrsion, ntitle, nsym,
     & nbft,   ninfo,  nenrgy, nmap,
     & title,  nbpsy,  slabel, info,
     & bfnlab, ietype, energy, imtype,
     & map,    mapdim, ierr )
cgk debug
      write(*,*)'gk: back from sifzwh'
cgk end
      return
      end
