      subroutine argos_cafe_pr_angles(lfnout,nsh,msh,msp,iang,
     + ang,rang,rub,snam)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer lfnout,nsh,msh,msp
      integer iang(msh,4)
      real*8 ang(msh,msp,mset),rang(msh,2),rub(msh,2)
      character*16 snam(nsatot)
c
      integer i,j
c
      if(np.gt.1) call ga_dgop(mcf_52,rang,2*msh,'+')
      if(np.gt.1.and.iffld.eq.2) call ga_dgop(mcf_81,rub,2*msh,'+')
      if(me.ne.0) return
c
      write(lfnout,1000)
 1000 format(/,' Solute angles',/,
     + '      i      j      k     atom_i           atom_j           ',
     + 'atom_k            c     th_0       C          th         U',/)
c
      do 1 i=1,nsh
      write(lfnout,1001) (iang(i,j),j=1,3),
     + (snam(iang(i,j))(12:16)//':'//snam(iang(i,j))(1:11),j=1,3),
     + iang(i,4),(ang(i,j,iset),j=1,2),(rang(i,j),j=1,2)
 1001 format(3i7,3a17,i7,f10.5,e12.5,f10.5,f12.6)
    1 continue
c
      if(iffld.eq.2) then
c
      write(lfnout,1002)
 1002 format(/,' Solute Urey-Bradley terms',/,
     + '      i      k     atom_i           atom_k           ',
     + ' c     s_0        C          s          U',/)
c
      do 2 i=1,nsh
      if(ang(i,3,iset).gt.0.0d0) then
      write(lfnout,1003) iang(i,1),iang(i,3),
     + snam(iang(i,1))(12:16)//':'//snam(iang(i,1))(1:11),
     + snam(iang(i,3))(12:16)//':'//snam(iang(i,3))(1:11),
     + iang(i,4),(ang(i,j,iset),j=3,4),(rub(i,j),j=1,2)
 1003 format(2i7,2a17,i7,f10.5,e12.5,f10.5,f12.6)
      endif
    2 continue
c
      endif
c
      return
      end
c $Id$
