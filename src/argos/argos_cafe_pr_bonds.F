      subroutine argos_cafe_pr_bonds(lfnout,nsb,msb,msp,ibnd,
     + bnd,rbnd,snam)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nsb,msb,msp
      integer ibnd(msb,3)
      real*8 bnd(msb,msp,mset),rbnd(msb,2)
      character*16 snam(nsatot)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_51,rbnd,2*msb,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute bonds',/,
     + '      i      j     atom_i           atom_j        ',
     + ' c        b_0        C          b          U',/)
c
      do 1 i=1,nsb
      if(ibnd(i,3).eq.0) then
      write(lfnout,1001) (ibnd(i,j),j=1,2),
     + (snam(ibnd(i,j))(12:16)//':'//snam(ibnd(i,j))(1:11),j=1,2),
     + (bnd(i,j,iset),j=1,2),(rbnd(i,j),j=1,2)
 1001 format(2i7,2a17,7x,f10.5,e12.5,f10.5,f12.6)
      else
      write(lfnout,1002) (ibnd(i,j),j=1,2),
     + (snam(ibnd(i,j))(12:16)//':'//snam(ibnd(i,j))(1:11),j=1,2),
     + (bnd(i,j,iset),j=1,2),rbnd(i,1)
 1002 format(2i7,2a17,' cnstr ',f10.5,e12.5,f10.5,f12.6)
      endif
    1 continue
c
      return
      end
c $Id$
