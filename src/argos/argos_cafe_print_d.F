      subroutine argos_cafe_print_d(lfnout,snam,dec)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout
      real*8 dec(6,nsatot),decr(6),dtot
      character*16 snam(nsatot)
c
      integer i,j
c
      write(lfnout,1001)
 1001 format(//,
     + ' SOLUTE ATOM DECOMPOSITION THERMODYNAMIC INTEGRATION',//,
     + '     Atom',13x,'Solvent     Solvent      Solute      Solute',
     + 6x,'Solute      Solute       Total',/,
     + 21x,'v.d.Waals   Electrost   v.d.Waals   Electrost',5x,
     + 'Bonded        Mass',/)
      do 1 i=1,nsatot
      dtot=zero
      do 2 j=1,6
      dtot=dtot+dec(j,i)
    2 continue
      write(lfnout,1009) snam(i),(dec(j,i),j=1,6),dtot
 1009 format(1x,a16,7f12.3)
    1 continue
c
      write(lfnout,1002)
 1002 format(//,
     + ' SOLUTE SEGMENT DECOMPOSITION THERMODYNAMIC INTEGRATION',//,
     + '     Atom',13x,'Solvent     Solvent      Solute      Solute',
     + 6x,'Solute      Solute       Total',/,
     + 21x,'v.d.Waals   Electrost   v.d.Waals   Electrost',5x,
     + 'Bonded        Mass',/)
      do 3 j=1,6
      decr(j)=dec(j,1)
    3 continue
      do 4 i=2,nsatot
      if(snam(i)(11:16).eq.snam(i-1)(11:16)) then
      do 5 j=1,6
      decr(j)=decr(j)+dec(j,i)
    5 continue
      else
      dtot=zero
      do 6 j=1,6
      dtot=dtot+decr(j)
    6 continue
      write(lfnout,1010) snam(i-1)(1:5),snam(i-1)(11:16),
     + (decr(j),j=1,6),dtot
 1010 format(1x,a5,5x,a6,7f12.3)
      do 7 j=1,6
      decr(j)=dec(j,i)
    7 continue
      endif
    4 continue
      dtot=zero
      do 8 j=1,6
      dtot=dtot+decr(j)
    8 continue
      write(lfnout,1010) snam(nsatot)(1:5),snam(nsatot)(11:16),
     + (decr(j),j=1,6),dtot
c
      return
      end
c $Id$
