      subroutine esp_wrtinp(q)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      integer i,j,jcon,lcon
      real*8 q(natoms)
c
      if(me.ne.0) return
c
      if(util_print('output',print_medium)) then
      write(lfnout,601) natoms
  601 format(' Atom parameters',//,
     + ' Number of atoms is ',t45,i12)
c
      write(lfnout,602) nbf
  602 format(' Number of basis functions is',t45,i12)
c
      write(lfnout,603) mxgrid
  603 format(//,' Grid parameters',//,
     + ' Maximum number of grid points is',t45,i12)
c
      write(lfnout,604) ngrid
  604 format(' Number of grid points is',t45,i12)
c
      if(newgrd.ne.0) then
c
      write(lfnout,605) cau2nm*rcut
  605 format(' Grid range',t45,f12.6,' nm')
c
      write(lfnout,606) cau2nm*spac
  606 format(' Grid spacing',t45,f12.6,' nm')
c
      write(lfnout,607) probe
  607 format(' Probe radius',t45,f12.6,' nm')
c
      write(lfnout,608) factor
  608 format(' Atom radius factor',t45,f12.6)
c
      if(lscrn.ne.0) then
      write(lfnout,609) scrtol
  609 format(/,' Schwarz screening with tolerance',1pe8.1)
      endif
c
      else
c
      call util_file_name('grid',.false.,.false.,grdfil)
      write(lfnout,618) grdfil(1:index(grdfil,' ')-1)
  618 format(/,' Grid read from ',a)
c
      endif
c
      if(nconst.gt.1) then
      write(lfnout,610)
  610 format(//,' Constraint parameters',/)
      jcon=0
      do 1 i=1,nconst-1
      jcon=jcon+1
      lcon=icon(jcon)
      if(lcon.gt.0) then
      write(lfnout,611) rcon(i),(icon(jcon+j),j=1,lcon)
  611 format(' Constraint charge',f12.6,' on atoms',10i5,/,(39x,10i5))
      else
      lcon=-lcon
      write(lfnout,1611) rcon(i),(icon(jcon+j),j=1,lcon)
 1611 format(' Constraint charge',f12.6,' in 2nd stage on atoms',10i5,
     + /,(39x,10i5))
      endif
      jcon=jcon+lcon
    1 continue
      endif
c
      if(irest.gt.0) then
      write(lfnout,612)
  612 format(//,' Restraint parameters',/)
      endif
c
      if(irest.eq.1) then
      write(lfnout,613) resa
  613 format(' Harmonic restraint to zero charge',/,
     + '    with scaling',t30,f12.6,/)
      endif
c
      if(irest.eq.2) then
      write(lfnout,614) resa,resb,maxit,toler
  614 format(' Hyperbolic restraint to zero charge',/,
     + '    with scaling',t30,f12.6,/,
     + '         tightness',t30,f12.6,/,
     + '         maximum iterations',t30,i5,/,
     + '         tolerance',t30,f12.6,/)
      endif
c
      if(irest.gt.0) then
      if(ihfree.gt.0) then
      write(lfnout,615)
  615 format(/,' Hydrogen atoms are excluded from the restraints')
      endif
      endif
c
      write(lfnout,616)
  616 format(/,' Atomic radii',/)
      do 2 i=1,107
      do 3 j=1,natoms
      if(int(q(j)).eq.i) then
      write(lfnout,617) i,radius(i)
  617 format(i5,f12.6)
      goto 2
      endif
    3 continue
    2 continue
c
      endif
c
      return
      end
