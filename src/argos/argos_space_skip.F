      logical function argos_space_skip(lfntri,nskip)
c
      implicit none
c
#include "argos_space_common.fh"
c
      integer lfntri,nskip
c
      integer i
      character*80 card
c
      if(me.eq.0) then
      do 1 i=1,nskip
    2 continue
      read(lfntri,1000,end=9999) card
 1000 format(a)
      if(card(1:5).ne.'frame') goto 2
    1 continue
      endif
c
      argos_space_skip=.true.
c
      return
c
 9999 continue
c
      argos_space_skip=.false.
c
      return
      end
c $Id$
