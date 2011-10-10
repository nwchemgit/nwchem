      subroutine argos_diana_initag(tag,val)
c
c     implicit none
c
#include "argos_diana_common.fh"
c
      character*24 tag(nsa,2)
      real*8 val(nsa,2)
c
      integer i
c
      do 1 i=1,nsa
      tag(i,1)='                        '
      tag(i,2)='                        '
      val(i,1)=1.0d0
      val(i,2)=1.0d0
    1 continue
c
      do 2 i=1,nwa
      wtag(i,1)='                        '
      wtag(i,2)='                        '
      wval(i,1)=1.0d0
      wval(i,2)=1.0d0
    2 continue
c
      return
      end
c $Id$
