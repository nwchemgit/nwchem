      logical function pre_natyps(lfnout,catm,matm,natm,natt,nats)
c
c $Id: pre_natyps.F,v 1.5 2001-10-23 20:46:28 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,matm,natm,natt,nats
      character*6 catm(4,matm)
c
      integer i,j,ii,jj,iii
c
      natt=0
      nats=0
c
      if(util_print('topology',print_high)) then
      write(lfnout,1000)
 1000 format(' Atom types',/)
      endif
c
      do 1 i=1,natm
      do 2 ii=1,3
      do 3 j=1,i-1
      do 4 jj=1,3
      if(catm(ii+1,i).eq.catm(jj+1,j)) goto 2
    4 continue
    3 continue
      do 5 iii=1,ii-1
      if(catm(ii+1,i).eq.catm(iii+1,i)) goto 2
    5 continue
      natt=natt+1
      if(util_print('topology',print_high)) then
      write(lfnout,1001) natt,catm(ii+1,i)
 1001 format(1x,i5,1x,a6)
      endif
    2 continue
    1 continue
c
      if(util_print('topology',print_high)) then
      write(lfnout,1002)
 1002 format(' ')
      endif
c
      do 6 i=1,natm
      do 7 j=1,i-1
      if(catm(2,i).eq.catm(2,j).and.catm(3,i).eq.catm(3,j).and.
     + catm(4,i).eq.catm(4,j)) goto 6
    7 continue
      nats=nats+1
    6 continue
c
      pre_natyps=.true.
      return
      end
