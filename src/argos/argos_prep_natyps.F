      logical function argos_prep_natyps(lfnout,nparms,mparms,
     + catm,matm,natm,natt,nats)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,matm,natm,natt,nats,nparms,mparms
      character*6 catm(mparms,matm)
c
      integer i,j,k,ii,jj,iii
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
      do 2 ii=1,nparms
      do 3 j=1,i-1
      do 4 jj=1,nparms
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
      ii=0
      do 8 k=1,nparms
      if(catm(k+1,i).eq.catm(k+1,j)) ii=ii+1
    8 continue
      if(ii.eq.nparms) goto 6
    7 continue
      nats=nats+1
    6 continue
c
      argos_prep_natyps=.true.
      return
      end
