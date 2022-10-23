      subroutine sifrh2(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    title,
     & nbpsy,   slabel,  info,    bfnlab,
     & ietype,  energy,  imtype,  map,
     & ierr )
c
c  read header_2 from the standard integral file.
c
c  input:
c  aoints = input file unit number.
c  ntitle = number of titles.
c  nsym = number of symmetry blocks.
c  nbft = total number of basis functions.
c  ninfo = number of record-definition parameters.
c  nenrgy = number of core energies.  the first element must be the
c           nuclear repulsion energy.
c  nmap = number of optional map vectors.
c
c  output:
c  title*80(1:ntitle) = identifying titles.
c  nbpsy(1:nsym) = number of basis functions per symmetry block.
c  slabel*4(1:nsym) = symmetry labels.
c  info(1:ninfo) = record-definition parameters.
c  bfnlab*8(1:nbft) = basis function labels.
c  ietype(1:nenrgy) = core energy types.
c  energy(1:nenrgy) = core energies.
c  imtype(1:nmap) = map vector types.
c  map(1:nbft,1:nmap) = basis function map vectors.
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer aoints, ntitle, nsym, nbft, ninfo, nenrgy, nmap, ierr
      character*80 title(*)
      integer nbpsy(*)
      character*4 slabel(*)
      integer info(*)
      character*8 bfnlab(*)
      integer ietype(*)
      real*8 energy(*)
      integer imtype(*)
      integer map(*)
c
c
      integer mapdim, nbftx, isym
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      mapdim = max ( nmap, 1 )
c     write(*,*) 'sifrh2:1: nenrgy=',nenrgy
      call sifzh2(
     & aoints, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   title,
     & nbpsy,  slabel, info,   bfnlab,
     & ietype, energy, imtype, map,
     & mapdim, ierr )
c     write(*,*) 'sifrh2:2: nenrgy=',nenrgy
c
      if ( ierr.ne.0 ) return
c
      nbftx=0
      do 10 isym=1,nsym
         if ( nbpsy(isym) .lt. 0 ) then
            call bummer('sifrh2: nbpsy(isym)=',nbpsy(isym),wrnerr)
            ierr = -2
            return
         endif
         nbftx=nbftx+nbpsy(isym)
10    continue
      if ( nbftx.ne.nbft ) then
         call bummer('sifrh2: (nbftx-nbft)=',(nbftx-nbft),wrnerr)
         ierr = -3
         return
      endif
c     write(*,*) 'sifrh2:3: nenrgy=',nenrgy
c
      return
      end
