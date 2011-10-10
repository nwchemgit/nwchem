      subroutine argos_cafe_print_der(lfnout,number,snam,dera)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout,number
      real*8 dera(6,nsatot)
      character*16 snam(nsatot)
c
      integer i,j
c
      call ga_dgop(mcf_68,dera,5*nsatot,'+')
c
      if(me.eq.0) then
c
      write(lfnout,1000) number
 1000 format(//,' Solute MCTI derivatives from last',i8,' steps',//,
     + '     Atom',13x,'Solvent     Solvent      Solute      Solute',
     + 6x,'Solute      Solute',/,
     + 21x,'v.d.Waals   Electrost   v.d.Waals   Electrost',5x,
     + 'Bonded        Mass',/)
c
      do 1 i=1,nsatot
      write(lfnout,1001) snam(i),(dera(j,i)/dble(number),j=1,6)
 1001 format(1x,a16,6f12.3)
    1 continue
c
      endif
c
      return
      end
c $Id$
