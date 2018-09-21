      subroutine esp_out(tag,x,qf,g,v)
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
      character*16 tag(natoms)
      real*8 x(3,natoms),qf(ndim,6),qsum(6)
      real*8 g(3,mxgrid),v(mxgrid)

      logical util_io_unit
      external util_io_unit
c
      integer i,j,k
      integer fn
      character*16 Atag
      character*80 name,espfilename
      real*8 coord(3)
c
      if(me.eq.0) then

      call util_file_name('esp',.false.,.false.,chgfil)
      open(unit=lfnchg,file=chgfil,form='formatted',status='unknown')
      rewind(lfnchg)
c
      call util_file_name('q',.false.,.false.,chgfil)
      open(unit=lfnchg,file=chgfil,form='formatted',status='unknown')
      rewind(lfnchg)
c
      write(lfnchg,'(i7,i3)') natoms,nf
c      
      if(.not.util_io_unit(80,90,fn))
     +  call errquit("cannot get io unit in esp_out",0,0)
      call util_file_name('esp',.false.,.false.,espfilename)
      open(unit=fn,file=espfilename,status='unknown',form="formatted")
cc
      if(util_print('output',print_medium)) then
      write(lfnout,600)
  600 format(//,'    Atom',t23,'Coordinates',t61,'Charge',/)
      write(lfnout,601) (labelf(i),i=1,nf)
      write(lfnout,601) (labelc(i),i=1,nf)
  601 format(t45,6a12)
      write(lfnout,602)
  602 format(' ')
      endif
      qsum(1)=zero
      qsum(2)=zero
      qsum(3)=zero
      qsum(4)=zero
      qsum(5)=zero
      qsum(6)=zero
      write(fn,902) natoms
      do 1 i=1,natoms
      if(util_print('output',print_medium)) then
      write(lfnout,603)
     + i,tag(i)(1:2),(cau2nm*x(j,i),j=1,3),(qf(i,j),j=1,nf)
  603 format(i5,1x,a2,3f12.6,6f12.6)
      endif
      write(fn,903)
     + tag(i)(1:2),(cau2ang*x(j,i),j=1,3),(qf(i,nf))
  902 format(i5,/)
  903 format(a2,3f12.6,6f12.6)
      write(lfnchg,703)
     + tag(i)(1:2),(cau2nm*x(j,i),j=1,3),(qf(i,j),j=1,nf) 
  703 format(a2,3f12.6,6f12.6)
      do 2 j=1,nf
      qsum(j)=qsum(j)+qf(i,j)
    2 continue
    1 continue
      close(fn)
      close(unit=lfnchg)
      if(util_print('output',print_medium)) then
      write(lfnout,604) ('------------',i=1,nf)
  604 format(44x,6a12)
      write(lfnout,605) (qsum(i),i=1,nf)
  605 format(44x,6f12.6)
c
      write(lfnout,606) (dipmom(i),i=1,nf)
  606 format(/,' Dipole moment',t45,6f12.6)
c
      write(lfnout,607) (quamom(1,1,i),i=1,nf)
  607 format(/,' Quadrupole moment Qxx',t45,6f12.6)
c      write(lfnout,608) (quamom(1,2,i),i=1,nf)
c  608 format('                   Qxy',t49,4f12.6)
c      write(lfnout,609) (quamom(1,3,i),i=1,nf)
c  609 format('                   Qxz',t49,4f12.6)
      write(lfnout,610) (quamom(2,2,i),i=1,nf)
  610 format('                   Qyy',t45,6f12.6)
c      write(lfnout,611) (quamom(2,3,i),i=1,nf)
c  611 format('                   Qyz',t49,4f12.6)
      write(lfnout,612) (quamom(3,3,i),i=1,nf)
  612 format('                   Qzz',t45,6f12.6)
c
      write(lfnout,613) (rms(1,i),i=1,nf)
  613 format(/,' RMS deviation kJ/mol',t45,6f12.6)
      write(lfnout,614) (rms(2,i),i=1,nf)
  614 format(/,' RRMS %',t45,6f12.6)
      endif
c
      endif
c
c     put charges on rtdb
c     -------------------
c
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
     + call md_abort('esp_out: geom_cent_get failed',9999)
      charge=qf(j,i)
      if(.not.geom_cent_set(igeom,j,Atag,coord,charge))
     + call md_abort('esp_out: geom_cent_set failed',9999)
    4 continue
      if(k.gt.0) then
      if(.not.geom_rtdb_store(irtdb,igeom,name))
     + call md_abort('esp_out: geom_rtdb_store failed',9999)
      else
      call md_abort('error in esp_out',9999)
      endif
      if(i.eq.nf) then
      if(.not.geom_rtdb_store(irtdb,igeom,'esp_fit:geometry'))
     + call md_abort('esp_out: geom_rtdb_store failed',9999)
      endif
    3 continue
c
      return
      end
