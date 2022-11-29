      subroutine sifzwh(
     & aoints,  vrsion,  ntitle,  nsym,
     & nbft,    ninfo,   nenrgy,  nmap,
     & title,   nbpsy,   slabel,  info,
     & bfnlab,  ietype,  energy,  imtype,
     & map,     mapdim,  ierr )
c
c  low-level header-writing routine.
c
c  *** this routine should not be called directly by user programs. ***
c
c  this version writes imtype(*) and map(*).
c  see sifwh() for argument description.
c
c  if ( nmap.ne.0 ) then
c     mapdim=nmap
c  else
c     mapdim=1
c  endif
c
c  26-jun-89 written by ron shepard.
c
      implicit none
c
      integer  aoints, vrsion, ntitle, nsym,   nbft,
     & ninfo,  nenrgy, nmap,   mapdim, ierr
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
cgk debug
*      write(6,*)'gk: in sifzwh, aoints=',aoints
*      write(6,*)'gk:      vrsion=',vrsion
*      write(6,*)'gk:      ntitle=',ntitle
*      write(6,*)'gk:      nsym=',nsym
*      write(6,*)'gk:      nbft=',nbft
*      write(6,*)'gk:      ninfo=',ninfo
*      write(6,*)'gk:      nenrgy=',nenrgy
*      write(6,*)'gk:      nmap=',nmap
cgk end
      write(aoints,iostat=ierr)
     & vrsion, ntitle, nsym, nbft, ninfo, nenrgy, nmap
cgk debug
*      write(6,*)
*     & vrsion, ntitle, nsym, nbft, ninfo, nenrgy, nmap
cgk end
c
      if ( ierr .ne. 0 ) return
c
      if ( nmap .eq. 0 ) then
         write(aoints,iostat=ierr)
     &    title, nbpsy, slabel, info, bfnlab, ietype, energy
      else
         write(aoints,iostat=ierr)
     &    title, nbpsy, slabel, info, bfnlab, ietype, energy,
     &    imtype, map
      endif
c
      return
      end
