      subroutine smd_charge_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_charge,sp_param,sp_atom,sp_type
      character*32 tag,pname
      character*255 filename
      logical result

      pname = "smd_charge_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "charge"
      call smd_system_get_component(sp_charge,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "type"
      call smd_system_get_component(sp_type,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_charge_init(sp_charge)
      call smd_chargefile_input(filename,result)
      if(result) then
        call smd_charge_read(filename) 
      else
        call smd_charge_set()
      end if
c
      return
      end

      subroutine smd_charge_init(sp_charge)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_charge_data.fh"
c     
      character*(*) sp_charge
c
      character*32 pname
      character*255 filename
      integer na
c
      integer ndims, type, ga_charges
      integer dims(1), chunk(1)
      integer rtdb
c
      pname = "smd_charge_init"
c
c      write(*,*) "in "//pname
c
      type = MT_DBL
      ndims = 1
c
c     get total number of atoms 
c     -------------------------
      call smd_atom_ntot(na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c     create charge data structures
c     ---------------------------
      call smd_namespace_create(sp_charge)
      call smd_data_create(sp_charge,"charge",na,MT_DBL)
      return
      end

      subroutine smd_charge_set()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 pname,tag
      integer na
      logical result
      integer i_q,i_q0,i_id
      integer i,j
      double precision q0
c
      pname = "smd_charge_set"
c
c      write(*,*) "in "//pname
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      tag = "param:charge"
      call smd_get_ind(tag,i_q0,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      tag = "type:id"
      call smd_get_ind(tag,i_id,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      do i=1,na
       j=int_mb(i_id+i-1)
       q0=dbl_mb(i_q0+j-1)
       dbl_mb(i_q+i-1)=q0
      end do

      return
      end

      subroutine smd_charge_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer un
c
      character*72 sp_charge
      character*32 pname,tag
      integer na
      logical result
      integer i_q
      integer i
c
      pname = "smd_charge_print"
c
c      write(*,*) "in "//pname
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      do i=1,na
       write(un,*) i,dbl_mb(i_q+i-1)
      end do

      return
      end

      subroutine smd_charge_write(filename)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) filename
c
      integer un
      character*72 sp_charge
      character*32 pname,tag
      integer na
      logical result
      integer i_q
      integer i
c
      pname = "smd_charge_write"
c
c      write(*,*) "in "//pname
c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
      open(unit=un,status="unknown",form="formatted",file=filename)
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      do i=1,na
       write(un,*) dbl_mb(i_q+i-1)
      end do

      close(un)

      return
      end

      subroutine smd_charge_write_default()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(255) filename
      logical result
      call smd_chargefile_output(filename,result)
c
      if(result) then
        call smd_charge_write(filename) 
      end if
      return
      end

      subroutine smd_charge_read(filename)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) filename
c
      integer un
      character*72 sp_charge
      character*32 pname,tag
      integer na
      logical result
      integer i_q
      integer i
c
      pname = "smd_charge_read"
c
c      write(*,*) "in "//pname
c
      if(.not.util_get_io_unit(un))
     >   call errquit("cannot get file number",0,0)
      open(unit=un,status="old",form="formatted",file=filename)
c
      tag = "charge"
      call smd_get_ind_dim(tag,i_q,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, rtdb_err)

      do i=1,na
       read(un,*) dbl_mb(i_q+i-1)
      end do

      close(un)

      return
      end

      subroutine smd_chargefile_output(filename,result)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "global.fh"

c
      character*(*) filename
      logical result
c
      character*30 pname

      pname = "smd_chargefile_input"

      result = .true.
      call smd_rtdb_get_string("smd:charge:output",1,
     >                           filename,result)

      end

      subroutine smd_chargefile_input(filename,result)
      implicit none
#include "rtdb.fh"
#include "smd_rtdb_data.fh"
#include "mafdecls.fh"
#include "global.fh"

c
      character*(*) filename
      logical result
c
      character*30 pname

      pname = "smd_chargefile_input"

      result = .true.
      call smd_rtdb_get_string("smd:charge:input",1,
     >                           filename,result)

      end
c $Id$
