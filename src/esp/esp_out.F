      subroutine esp_out(tag,x,q)
c
c $Id: esp_out.F,v 1.1 1997-10-06 23:00:36 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      character*16 tag(natoms)
      real*8 x(3,natoms),q(natoms),qsum
c
      integer i,j
c
      if(me.ne.0) return
      if(irest.eq.2) then
      write(6,605) niter
  605 format(//,' Number of iterations',t59,i5,/)
      endif
c
      write(6,600)
  600 format('    Atom',t27,'Coordinates',t63,'Charge',/)
      qsum=0.0d0
      do 1 i=1,natoms
      write(6,602) i,tag(i),(cau2nm*x(j,i),j=1,3),q(i)
  602 format(i5,1x,a16,3f12.6,f12.6)
      qsum=qsum+q(i)
    1 continue
      write(6,603) qsum
  603 format(58x,'------------',/,58x,f12.6)
c
      write(6,604) dipmom
  604 format(/,' Dipole moment',t59,f12.6)
      return
      end
