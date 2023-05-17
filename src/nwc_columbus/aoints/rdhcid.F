c
c**********************************************************************
c
*deck rdhcid
      subroutine rdhcid(denfl,infocd,ntitle,title,nsym,nmot,
     &                  nbpsy)
c
c  read the header info from the ci density matrix file.
c
c  written by ron shepard.
c
      implicit none
c
      character*80 title(*)
      integer infocd(*)
      integer ntitle,nsym,nmot,ninfo,nenrgy,ierr,nmap
      integer nmotx
      parameter (nmotx=255)
      integer i, ntitmx, denfl,idummy, nsymx, nengmx
      parameter( ntitmx = 30, nsymx=1, nengmx=20)
      real*8 cienrg(nengmx)
      character*4 slabel(nsymx)
      character*8 cidlab(nmotx)
      integer nbpsy(*),cietyp(nengmx)
c
      real*8    one,     zero
      parameter(one=1d0, zero=0d0)
c
      integer nmotn, nsymn
c
c  bummer error types.
      integer   wrnerr,  nfterr,  faterr
      parameter(wrnerr=0,nfterr=1,faterr=2)
c
      rewind denfl
      call sifrh1(
     1 denfl, ntitle, nsym, nmot,
     2 ninfo, nenrgy, nmap, ierr)
c
cgk debug
!       write(*,*)'gk: in rdhcid, ntitle=',ntitle
!       write(*,*)'gk:  nsym=',nsym
!       write(*,*)'gk:  nmot = ', nmot
!       write(*,*)'gk:  ninfo = ', ninfo
!       write(*,*)'gk:  nenrgy = ', nenrgy
!       write(*,*)'gk:  nmap = ', nmap
cgk end
      if (ierr .ne. 0) then
         call bummer('rdhcid: from sifrh1, ierr=',ierr,faterr)
      elseif (nenrgy .gt. nengmx) then
         call bummer('rdhcid: from sifrh1, nenrgy=',ntitle,faterr)
      elseif (nsym .gt. nsymx) then
         call bummer('rdhcid: from sifrh1, nsymx=',nsymx,faterr)
      elseif (ntitle .gt. ntitmx) then
         call bummer('rdhcid: from sifrh1, ntitle=',ntitle,faterr)
      elseif (ninfo .gt. 10) then
         call bummer('rdhcid: from sifrh1, ninfo=',ninfo,wrnerr)
      endif
c     #ignore map
      nmap = 0
      call sifrh2(
     1 denfl,   ntitle,  nsym,    nmot,
     2 ninfo,   nenrgy,  nmap,    title,
     3 nbpsy,   slabel,  infocd,  cidlab,
     4 cietyp,  cienrg,  idummy,  idummy,
     5 ierr)
cgk debug
!       write(*,*)'gk: in rdhcid, ntitle=',ntitle
!       write(*,*)'gk:  nsym=',nsym
!       write(*,*)'gk:  nmot = ', nmot
!       write(*,*)'gk:  ninfo = ', ninfo
!       write(*,*)'gk:  nenrgy = ', nenrgy
!       write(*,*)'gk:  nmap = ', nmap
!       write(*,*)'gk:  title = ', title(1:6)
!       write(*,*)'gk:  nbpsy = ', nbpsy(1)
!       write(*,*)'gk:  slabel = ', slabel
!       write(*,*)'gk:  infocd = ', infocd(1:5)
!       write(*,*)'gk:  cidlab = ', cidlab
!       write(*,*)'gk:  cietyp = ', cietyp
!       write(*,*)'gk:  cienrg = ', cienrg
cgk end
c      write(*,*)'ntitle=',ntitle
c
20    continue
c
      write(*,6020)(title(i),i=1,ntitle)
6020  format(/' ci density file header information:'/(1x,a))
      write(*,6040)(i,nbpsy(i),i=1,nsym)
6040  format(' (isym:nbpsy)',8(i3,':',i3))
c
      return
      end
