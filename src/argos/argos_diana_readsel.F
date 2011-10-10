      subroutine argos_diana_readsel(lfn,fil,iwsel,isel)
c
      implicit none
c
#include "argos_diana_common.fh"
c
      integer lfn,iwsel(mwm),isel(msa)
      character*255 fil
c
      integer i
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='old')
      rewind(lfn)
      read(lfn,1000) (iwsel(i),i=1,nwm)
      read(lfn,1000) (isel(i),i=1,nsa)
 1000 format(40i2)
      close(unit=lfn)
      write(*,'(/,a,a)') ' Selection read from ',fil(1:index(fil,' ')-1)
      endif
c
      return
      end
c $Id$
