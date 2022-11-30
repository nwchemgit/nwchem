      subroutine sifzh2(
     & aoints,  ntitle,  nsym,    nbft,
     & ninfo,   nenrgy,  nmap,    title,
     & nbpsy,   slabel,  info,    bfnlab,
     & ietype,  energy,  imtype,  map,
     & mapdim,  ierr )
c
c  low-level header2 reading routine.
c
c  *** this routine should not be called directly by user programs. ***
c
c  see sifrh2() for argument description.
c
c  if ( nmap.eq.0 ) then
c     mapdim=1 ;imtype(*) and map(*) are not referenced.
c  else
c     mapdim=nmap ;imtype(*) and map(*) are read.
c  endif
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, ntitle, nsym,   nbft,   ninfo,
     & nenrgy, nmap,   mapdim, ierr
      character*80     title(ntitle)
      integer  nbpsy(nsym)
      character*4      slabel(nsym)
      integer  info(ninfo)
      character*8      bfnlab(nbft)
      integer  ietype(nenrgy)
      real*8   energy(nenrgy)
      integer  imtype(mapdim)
      integer  map(nbft,mapdim)
c
c     write(*,*) 'sifzh2:1: nenrgy=',nenrgy
c    & , ntitle,nsym,nbft,mapdim
      if ( nmap.eq.0 ) then
c     write(*,*) 'sifzh2:2: nenrgy=',nenrgy
         read(aoints,iostat=ierr)
     &    title,nbpsy,slabel,info,bfnlab,ietype,energy
c     write(*,*) 'sifzh2:3: nenrgy=',nenrgy
      else
c     write(*,*) 'sifzh2:4: nenrgy=',nenrgy
c    & , ntitle,nsym,nbft,mapdim
         read(aoints,iostat=ierr)
     &    title,nbpsy,slabel,info,bfnlab,ietype,energy,imtype,map
c     write(*,*) 'sifzh2:5: nenrgy=',nenrgy
c    & , ntitle,nsym,nbft,mapdim
      endif
c
      return
      end
