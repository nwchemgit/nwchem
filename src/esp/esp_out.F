      subroutine esp_out(tag,x,qf,g,v)
c
c $Id: esp_out.F,v 1.3 1997-10-07 23:41:45 d3j191 Exp $
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
      real*8 x(3,natoms),qf(ndim,4),qsum(4)
      real*8 g(3,mxgrid),v(mxgrid)
c
      integer i,j
c
      if(me.ne.0) return
c
      write(6,600)
  600 format('    Atom',t17,'Coordinates',t71,'Charge',/)
      write(6,601) (labelf(i),i=1,nf)
      write(6,601) (labelc(i),i=1,nf)
  601 format(t50,4a12)
      write(6,602)
  602 format(' ')
      qsum(1)=zero
      qsum(2)=zero
      qsum(3)=zero
      qsum(4)=zero
      do 1 i=1,natoms
      write(6,603) i,tag(i)(1:6),(cau2nm*x(j,i),j=1,3),(qf(i,j),j=1,nf)
  603 format(i5,1x,a6,3f12.6,4f12.6)
      do 2 j=1,nf
      qsum(j)=qsum(j)+qf(i,j)
    2 continue
    1 continue
      write(6,604) ('------------',i=1,nf)
  604 format(48x,4a12)
      write(6,605) (qsum(i),i=1,nf)
  605 format(48x,4f12.6)
c
      write(6,606) (dipmom(i),i=1,nf)
  606 format(/,' Dipole moment',t49,4f12.6)
c
      write(6,607) (rms(1,i),i=1,nf)
  607 format(/,' RMS deviation kJ/mol',t49,4f12.6)
      write(6,608) (rms(2,i),i=1,nf)
  608 format(/,' RMS deviation %',t49,4f12.6)
c
      if(idump.eq.1) then
      call util_file_name('grid',.false.,.false.,grdfil)
      open(unit=lfngrd,file=grdfil,form='formatted',status='unknown')
      do 3 i=1,ngrid
      write(lfngrd,'(4f20.10)') (g(j,i),j=1,3),v(i)
    3 continue
      close(unit=lfngrd)
      endif
c
      return
      end
