      subroutine esp_wrtinp
c
c $Id: esp_wrtinp.F,v 1.1 1997-10-06 23:00:37 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      integer i,j,jcon,lcon
c
      if(me.ne.0) return
c
      write(6,601) natoms
  601 format(' Atom parameters',//,
     + ' Number of atoms is ',t45,i12)
c
      write(6,602) nbf
  602 format(' Number of basis functions is',t45,i12)
c
      write(6,603) mxgrid
  603 format(//,' Grid parameters',//,
     + ' Maximum number of grid points is',t45,i12)
c
      write(6,604) ngrid
  604 format(' Number of grid points is',t45,i12)
c
      write(6,605) cau2nm*rcut
  605 format(' Grid range beyond van der Waals spheres',t45,f12.6,' nm')
c
      write(6,606) cau2nm*spac
  606 format(' Grid spacing',t45,f12.6,' nm')
c
      if(nconst.gt.1) then
      write(6,607)
  607 format(//,' Constraint parameters',/)
      jcon=0
      do 1 i=1,nconst-1
      jcon=jcon+1
      lcon=icon(jcon)
      write(6,608) rcon(i),(icon(jcon+j),j=1,lcon)
  608 format(' Constraint charge',f12.6,' on atoms',10i5)
      jcon=jcon+lcon
    1 continue
      endif
c
      if(irest.gt.0) then
      write(6,609)
  609 format(//,' Restraint parameters',/)
      endif
c
      if(irest.eq.1) then
      write(6,610) resa
  610 format(' Harmonic restraints fit to zero charge with weight',
     + f12.6)
      endif
c
      if(irest.eq.2) then
      write(6,611) resa,resb,maxit,toler
  611 format(' Hyperbolic restraints fit to zero charge with weights',
     + f12.6,' and ',f12.6,/,
     + t25,'maximum number of iterations',i5,/,
     + t25,'tolerance',f12.6)
      endif
c
      if(irest.gt.0) then
      if(ihfree.gt.0) then
      write(6,612)
  612 format(' Hydrogen atoms excluded from restraints')
      endif
      endif
c
      return
      end
