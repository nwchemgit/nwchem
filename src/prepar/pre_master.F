      logical function pre_master(lfnfrg,filfrg,lfnout)
c
      implicit none
c
#include "util.fh"
#include "pre_common.fh"
c
      integer lfnfrg,lfnout
      character*255 filfrg,card
c
      integer length
c
      length=index(filfrg,' ')-1
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filfrg(1:length)
 2000 format(' Attempt reading fragment',t40,a,/)
      endif
c
      open(unit=lfnfrg,file=filfrg(1:length),form='formatted',
     + status='old',err=9999)
c
    1 continue
c
      read(lfnfrg,1000,end=999,err=9999) card
 1000 format(a)
c
      if(card(1:8).eq.'#MASTER ') then
      filmst=card(9:255)
      else
      goto 1
      endif
c
  999 continue
c
      close(lfnfrg)
c
      pre_master=.true.
      return
c
 9999 continue
      pre_master=.false.
      return
      end
