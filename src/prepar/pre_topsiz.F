      logical function pre_topsiz(lfntop,filtop,lfnout,
     + nat,nsa,nwa,nsb,nqu,nseq)
c
c $Id$
c
      implicit none
c
      integer lfntop,lfnout,nsa,nwa,nat,nqu
      character*255 filtop
      character*1 cdummy
      character*80 card
c
      integer nwb,nwh,nwd,nwo,nwt,nwn,nsb,nsh,nsd,nso,nst,nsx
      integer i,j,k,len,npar,nseq
c
      len=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:len),form='formatted',
     + status='old',err=9999)
c
      do 1 i=1,4
      read(lfntop,1000) cdummy
 1000 format(a1)
    1 continue
      read(lfntop,1001) npar
      read(lfntop,1001) nat
      read(lfntop,1001) nqu
      read(lfntop,1001) nseq
      read(lfntop,1000) cdummy
 1001 format(i5)
      do 2 i=1,nat*npar
      read(lfntop,1000) cdummy
    2 continue
      do 3 i=1,nat
      do 4 j=i,nat
      do 5 k=1,npar
      read(lfntop,1000) cdummy
    5 continue
    4 continue
    3 continue
      do 6 i=1,nqu*npar
      read(lfntop,1000) cdummy
    6 continue
      do 7 i=1,nseq
      read(lfntop,1000) cdummy
    7 continue
      read(lfntop,1002) nwa,nwb,nwh,nwd,nwo,nwt,nwn
      read(lfntop,1002) nsa,nsb,nsh,nsd,nso,nst,nsx
 1002 format(7i7)
c
      close(unit=lfntop)
c
      pre_topsiz=.true.
      return
c
 9998 continue
      close(unit=lfntop)
 9999 continue
      pre_topsiz=.false.
      return
      end
