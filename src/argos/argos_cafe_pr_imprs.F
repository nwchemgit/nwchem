      subroutine argos_cafe_pr_imprs(lfnout,nso,mso,msp,
     + ioop,oop,roop,snam)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nso,mso,msp
      integer ioop(mso,5)
      real*8 oop(mso,msp,mset),roop(mso,2)
      character*16 snam(nsatot)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_54,roop,2*mso,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute out of plane torsions',/,
     + '      i      j      k      l     atom_i           ',
     + 'atom_j           atom_k           atom_l           ',
     + ' c  mul    theta_0     C         theta         U',/)
c
      do 1 i=1,nso
      write(lfnout,1001) (ioop(i,j),j=1,4),
     + (snam(ioop(i,j))(12:16)//':'//snam(ioop(i,j))(1:11),j=1,4),
     + ioop(i,5),(oop(i,j,iset),j=1,3),(roop(i,j),j=1,2)
 1001 format(4i7,4a17,i7,f5.1,f10.5,e12.5,f10.5,f12.6)
    1 continue
c
      return
      end
c $Id$
