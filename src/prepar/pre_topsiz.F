      logical function pre_topsiz(lfntop,filtop,lfnout,
     + nat,nsa,nwa,nsb,nqu)
c
c $Id: pre_topsiz.F,v 1.5 1998-11-06 21:18:09 d3j191 Exp $
c
      implicit none
c
      integer lfntop,lfnout,nsa,nwa,nat,nqu
      character*255 filtop
      character*1 cdummy
c
      integer nwb,nwh,nwd,nwo,nwt,nwn,nsb,nsh,nsd,nso,nst,nsx
      integer i,j,len
c
      len=index(filtop,' ')-1
      open(unit=lfntop,file=filtop(1:len),form='formatted',
     + status='old',err=9999)
c
      do 1 i=1,4
      read(lfntop,1000) cdummy
 1000 format(a1)
    1 continue
      read(lfntop,1001) nat,nqu
 1001 format(i5,41x,i5)
      do 2 i=1,nat
      read(lfntop,1000) cdummy
    2 continue
      do 3 i=1,nat
      do 4 j=i,nat
      read(lfntop,1000) cdummy
      read(lfntop,1000) cdummy
    4 continue
    3 continue
      do 5 i=1,nqu
      read(lfntop,1000) cdummy
    5 continue
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
