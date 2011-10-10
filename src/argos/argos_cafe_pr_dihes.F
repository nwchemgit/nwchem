      subroutine argos_cafe_pr_dihes(lfnout,nsd,msd,msp,
     + idih,dih,rdih,snam)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nsd,msd,msp
      integer idih(msd,5)
      real*8 dih(msd,msp,mset),rdih(msd,2)
      character*16 snam(nsatot)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_53,rdih,2*msd,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute torsions',/,
     + '      i      j      k      l     atom_i           ',
     + 'atom_j           atom_k           atom_l           ',
     + ' c  mul     phi_0      C          phi          U',/)
c
      do 1 i=1,nsd
      write(lfnout,1001) (idih(i,j),j=1,4),
     + (snam(idih(i,j))(12:16)//':'//snam(idih(i,j))(1:11),j=1,4),
     + idih(i,5),(dih(i,j,iset),j=1,3),(rdih(i,j),j=1,2)
 1001 format(4i7,4a17,i7,f5.1,f10.5,e12.5,f10.5,f12.6)
    1 continue
c
      return
      end
c $Id$
