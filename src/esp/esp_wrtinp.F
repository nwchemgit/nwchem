      subroutine esp_wrtinp
c
c $Id: esp_wrtinp.F,v 1.4 1998-06-08 21:31:58 d3j191 Exp $
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
      write(lfnout,605) cau2nm*rcut
  605 format(' Grid range',t45,f12.6,' nm')
c
      write(lfnout,606) cau2nm*spac
  606 format(' Grid spacing',t45,f12.6,' nm')
c
      if(nconst.gt.1) then
      write(lfnout,607)
  607 format(//,' Constraint parameters',/)
      jcon=0
      do 1 i=1,nconst-1
      jcon=jcon+1
      lcon=icon(jcon)
      write(lfnout,608) rcon(i),(icon(jcon+j),j=1,lcon)
  608 format(' Constraint charge',f12.6,' on atoms',10i5,/,(39x,10i5))
      jcon=jcon+lcon
    1 continue
      endif
c
      if(irest.gt.0) then
      write(lfnout,609)
  609 format(//,' Restraint parameters',/)
      endif
c
      if(irest.eq.1) then
      write(lfnout,610) resa
  610 format(' Harmonic restraint to zero charge',/,
     + '    with scaling',t30,f12.6,/)
      endif
c
      if(irest.eq.2) then
      write(lfnout,611) resa,resb,maxit,toler
  611 format(' Hyperbolic restraint to zero charge',/,
     + '    with scaling',t30,f12.6,/,
     + '         tightness',t30,f12.6,/,
     + '         maximum iterations',t30,i5,/,
     + '         tolerance',t30,f12.6,/)
      endif
c
      if(irest.gt.0) then
      if(ihfree.gt.0) then
      write(lfnout,612)
  612 format(/,' Hydrogen atoms are excluded from the restraints')
      endif
      endif
c
      endif
c
      return
      end
