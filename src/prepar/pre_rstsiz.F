      logical function pre_rstsiz(lfnrst,filrst,nwm,nwa,nsa)
c
c $Id: pre_rstsiz.F,v 1.1 1998-07-15 21:18:11 d3j191 Exp $
c
      implicit none
c
      integer lfnrst,nwm,nwa
      character*255 filrst
c
      character*1 cdummy
      integer i,nsa,npbtyp
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',err=9999)
      rewind(lfnrst)
c
      do 1 i=1,4
      read(lfnrst,1000) cdummy
 1000 format(a1)
    1 continue
      read(lfnrst,1001) npbtyp
 1001 format(i5)
      do 2 i=1,2
      read(lfnrst,1000) cdummy
    2 continue
c
      read(lfnrst,1002) nwm,nwa,nsa
 1002 format(2i10,10x,i10)
c
      close(unit=lfnrst)
c
      pre_rstsiz=.true.
      return
c
 9999 continue
      pre_rstsiz=.false.
      return
      end
