      subroutine sifrh1(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    ierr )
c
c  read header_1 from the standard integral file.
c
c  input:
c  aoints = input file unit number.
c
c  output:
c  ntitle = number of titles.
c  nsym = number of symmetry blocks.
c  nbft = total number of basis functions.
c  ninfo = number of record-definition parameters.
c  nenrgy = number of core energies.  the first element must be the
c           nuclear repulsion energy.
c  nmap = number of optional map vectors.
c  ierr = error return code.  0 for normal return.
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   ierr
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
      integer  verin
c
c     # bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      read(aoints,iostat=ierr)
     & verin, ntitle, nsym, nbft, ninfo, nenrgy, nmap
c
      if ( ierr.ne.0 ) then
         return
      elseif ( verin.ne.vrsion ) then
         call bummer('sifrh1: (verin-vrsion)=',(verin-vrsion),wrnerr)
         ierr = -2
         return
      elseif ( ntitle.le.0 .or. ntitle.gt.ntitmx ) then
         call bummer('sifrh1: ntitle=',ntitle,wrnerr)
         ierr = -3
         return
      elseif ( nsym.ne.1 .and. nsym.ne.2 .and. nsym.ne.4
     &    .and. nsym.ne.8 ) then
         call bummer('sifrh1: nsym=',nsym,wrnerr)
         ierr = -4
         return
      elseif ( nbft.le.0 ) then
         call bummer('sifrh1: nbft=',nbft,wrnerr)
         ierr = -5
         return
      elseif ( ninfo.lt.ninchk ) then
         call bummer('sifrh1: ninfo=',ninfo,wrnerr)
         ierr = -6
         return
      elseif ( nenrgy.le.0 ) then
         call bummer('sifrh1: nenrgy=',nenrgy,wrnerr)
         ierr = -7
         return
      elseif ( nmap.lt.0 ) then
         call bummer('sifrh1: nmap=',nmap,wrnerr)
         ierr = -8
         return
      endif
c
      return
      end
