      subroutine esp_out(tag,x,q,g,v)
c
c $Id: esp_out.F,v 1.2 1997-10-07 18:09:54 d3j191 Exp $
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
      real*8 g(3,mxgrid),v(mxgrid)
c
      integer i,j
c
      if(me.ne.0) return
      if(irest.eq.2) then
      write(6,601) niter
  601 format(//,' Number of iterations',t59,i5,/)
      endif
c
      write(6,602)
  602 format('    Atom',t27,'Coordinates',t63,'Charge',/)
      qsum=zero
      do 1 i=1,natoms
      write(6,603) i,tag(i),(cau2nm*x(j,i),j=1,3),q(i)
  603 format(i5,1x,a16,3f12.6,f12.6)
      qsum=qsum+q(i)
    1 continue
      write(6,604) qsum
  604 format(58x,'------------',/,58x,f12.6)
c
      write(6,605) dipmom
  605 format(/,' Dipole moment',t59,f12.6)
c
      write(6,606) rms(1)
  606 format(/,' RMS deviation',t59,f12.6,' kJ/mol')
      write(6,607) rms(2)
  607 format(/,' RMS deviation',t59,f12.6,' %')
c
      if(idump.eq.1) then
      call util_file_name('grid',.false.,.false.,grdfil)
      open(unit=lfngrd,file=grdfil,form='formatted',status='unknown')
      do 2 i=1,ngrid
      write(lfngrd,'(4f20.10)') (g(j,i),j=1,3),v(i)
    2 continue
      close(unit=lfngrd)
      endif
c
      return
      end
