      subroutine argos_prop_record()
c
      implicit none
c
#include "argos_prop.fh"
c
      integer i
      character*10 pdate,ptime
c
      if(me.ne.0) return
c
      if(.not.lfhdr) then
      call swatch(pdate,ptime)
      write(lfnprp,3000) nprop,pdate,ptime,np,npfft
 3000 format(i7,1x,2a10,2i5)
      write(lfnprp,3001) (pronam(ixp(i)),i=1,nprop)
 3001 format(a50)
      lfhdr=.true.
      endif
      write(lfnprp,3002)
 3002 format('frame')
      if(iprop.eq.0) then
      write(lfnprp,3003) (p(ixp(i)),i=1,nprop)
      else
      write(lfnprp,3003) (psump(ixp(i))/dble(nsump),i=1,nprop)
 3003 format(4(1pe12.5))
      endif
      nsump=0
      do 1 i=1,nprop
      psump(ixp(i))=zero
    1 continue
c
      return
      end
c $Id$
