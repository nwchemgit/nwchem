      subroutine esp_out(tag,x,qf,g,v)
c
c $Id: esp_out.F,v 1.9 1998-08-10 17:05:04 d3j191 Exp $
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
      character*16 tag(natoms)
      real*8 x(3,natoms),qf(ndim,4),qsum(4)
      real*8 g(3,mxgrid),v(mxgrid)
c
      integer i,j,k
      character*16 Atag
      character*80 name
      real*8 coord(3)
      logical done
c
      if(me.eq.0) then
c
      if(util_print('output',print_medium)) then
      write(lfnout,600)
  600 format('    Atom',t17,'Coordinates',t71,'Charge',/)
      write(lfnout,601) (labelf(i),i=1,nf)
      write(lfnout,601) (labelc(i),i=1,nf)
  601 format(t50,4a12)
      write(lfnout,602)
  602 format(' ')
      endif
      qsum(1)=zero
      qsum(2)=zero
      qsum(3)=zero
      qsum(4)=zero
      do 1 i=1,natoms
      if(util_print('output',print_medium)) then
      write(lfnout,603)
     + i,tag(i)(1:6),(cau2nm*x(j,i),j=1,3),(qf(i,j),j=1,nf)
  603 format(i5,1x,a6,3f12.6,4f12.6)
      endif
      do 2 j=1,nf
      qsum(j)=qsum(j)+qf(i,j)
    2 continue
    1 continue
      if(util_print('output',print_medium)) then
      write(lfnout,604) ('------------',i=1,nf)
  604 format(48x,4a12)
      write(lfnout,605) (qsum(i),i=1,nf)
  605 format(48x,4f12.6)
c
      write(lfnout,606) (dipmom(i),i=1,nf)
  606 format(/,' Dipole moment',t49,4f12.6)
c
      write(lfnout,607) (quamom(1,1,i),i=1,nf)
  607 format(/,' Quadrupole moment Qxx',t49,4f12.6)
c      write(lfnout,608) (quamom(1,2,i),i=1,nf)
c  608 format('                   Qxy',t49,4f12.6)
c      write(lfnout,609) (quamom(1,3,i),i=1,nf)
c  609 format('                   Qxz',t49,4f12.6)
      write(lfnout,610) (quamom(2,2,i),i=1,nf)
  610 format('                   Qyy',t49,4f12.6)
c      write(lfnout,611) (quamom(2,3,i),i=1,nf)
c  611 format('                   Qyz',t49,4f12.6)
      write(lfnout,612) (quamom(3,3,i),i=1,nf)
  612 format('                   Qzz',t49,4f12.6)
c
      write(lfnout,613) (rms(1,i),i=1,nf)
  613 format(/,' RMS deviation kJ/mol',t49,4f12.6)
      write(lfnout,614) (rms(2,i),i=1,nf)
  614 format(/,' RMS deviation %',t49,4f12.6)
      endif
c
      endif
c
c     put charges on rtdb
c     -------------------
c
      done=geom_rtdb_delete(irtdb,'esp:geometry')
      done=geom_rtdb_delete(irtdb,'cesp:geometry')
      done=geom_rtdb_delete(irtdb,'resp:geometry')
      done=geom_rtdb_delete(irtdb,'cresp:geometry')
      do 3 i=1,nf
      k=1
      if(labelf(i).ne.'      ESP   ') k=k+2
      if(labelc(i).ne.'            ') k=k+1
      if(k.eq.1) name='esp:geometry'
      if(k.eq.2) name='cesp:geometry'
      if(k.eq.3) name='resp:geometry'
      if(k.eq.4) name='cresp:geometry'
      do 4 j=1,natoms
      if(.not.geom_cent_get(igeom,j,Atag,coord,charge))
     + call errquit('esp_out: geom_cent_get failed',9999)
      charge=qf(j,i)
      if(.not.geom_cent_set(igeom,j,Atag,coord,charge))
     + call errquit('esp_out: geom_cent_set failed',9999)
    4 continue
      if(k.gt.0) then
      if(.not.geom_rtdb_store(irtdb,igeom,name))
     + call errquit('esp_out: geom_rtdb_store failed',9999)
      else
      call errquit('error in esp_out',9999)
      endif
    3 continue
c
      if(me.eq.0.and.idump.eq.1) then
      call util_file_name('grid',.false.,.false.,grdfil)
      open(unit=lfngrd,file=grdfil,form='formatted',status='unknown')
      do 10 i=1,ngrid
      write(lfngrd,'(4f20.10)') (g(j,i),j=1,3),v(i)
   10 continue
      close(unit=lfngrd)
      endif
c
      return
      end
