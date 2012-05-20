      block data initial_smd_energy_data
      implicit none
#include "smd_energy_data.fh"
c     
      data smd_energy_name /"total",
     c                      "kinetic",
     c                      "vdw",
     c                      "electrostatic",
     c                      "ewald self",
     c                      "ewald reciprocal",
     c                      "ewald excluded",
     c                      "ewald real"/


      end

      subroutine smd_energy_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_energy
      character*32 tag,pname
      logical result

      pname = "smd_energy_init_system"
c
      tag = "energy"
      call smd_system_get_component(sp_energy,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_energy_init(sp_energy)

      return
      end

      subroutine smd_energy_init(sp_energy)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_energy_data.fh"
c     
      character*(*) sp_energy
c
      character*32 pname
      integer na
c
      pname = "smd_energy_init"
c
c      write(*,*) "in "//pname
c
c     create energy data structures
c     ---------------------------
      call smd_namespace_create(sp_energy)
      call smd_data_create(sp_energy,"energy",smd_energy_nc,MT_DBL)
      call smd_data_create(sp_energy,"oenergy",smd_energy_nc,MT_LOG)
      return
      end

      subroutine smd_energy_set_component(aname,avalue)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_energy_data.fh"
c     
      character*(*) aname
      double precision avalue
c
      character*32 pname
      integer nrec
      integer i
      logical result
      logical ocase
      integer i_e
      integer i_oe
      character*32 sp_energy
      character*32 tag
c
      pname = "smd_energy_set_component"
c
      write(*,*) "in "//pname
c
      call smd_system_get_component(sp_energy,"energy",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no energy ',0,0)

      tag = "energy"
      call smd_data_get_index(sp_energy,tag,i_e,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      tag = "oenergy"
      call smd_data_get_index(sp_energy,tag,i_oe,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

c
c     case sensitive
      ocase = .true.
c
      nrec = smd_energy_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_energy_name,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      dbl_mb(i_e+i-1) = avalue
      log_mb(i_oe+i-1) = .true.
      
      return
      end

      subroutine smd_energy_unset_component(aname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_energy_data.fh"
c     
      character*(*) aname
c
      character*32 pname
      integer nrec
      integer i
      logical result
      logical ocase
      integer i_e
      integer i_oe
      character*32 sp_energy
      character*32 tag
cc
      pname = "smd_energy_set_component"
c
      write(*,*) "in "//pname
c
      call smd_system_get_component(sp_energy,"energy",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no energy ',0,0)

      tag = "energy"
      call smd_data_get_index(sp_energy,tag,i_e,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      tag = "oenergy"
      call smd_data_get_index(sp_energy,tag,i_oe,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
c     case sensitive
      ocase = .true.
c
      nrec = smd_energy_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_energy_name,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      log_mb(i_oe+i-1) = .false. 
      dbl_mb(i_e+i-1) = 0.0d0
      
      return
      end

      subroutine smd_energy_get_component(avalue,aname,oexist)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_energy_data.fh"
c     
      character*(*) aname
      double precision avalue
c
      character*32 pname
      integer nrec
      integer i
      logical oexist,result
      logical ocase
      integer i_e
      integer i_oe
      character*32 sp_energy
      character*32 tag
c
      pname = "smd_energy_get_component"
c
      write(*,*) "in "//pname
c
      call smd_system_get_component(sp_energy,"energy",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no energy ',0,0)

      tag = "energy"
      call smd_data_get_index(sp_energy,tag,i_e,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      tag = "oenergy"
      call smd_data_get_index(sp_energy,tag,i_oe,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
c     case sensitive
      ocase = .true.
c
      nrec = smd_energy_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_energy_name,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      avalue = dbl_mb(i_e+i-1)
      oexist = log_mb(i_oe+i-1)
      
      return
      end

      subroutine smd_energy_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_energy_data.fh"
c     
      integer un
c
      character*32 pname
      integer nrec
      integer i
      logical oexist,result
      logical ocase
      integer i_e
      integer i_oe
      character*32 sp_energy
      character*32 tag
c
      pname = "smd_energy_print"
c
      write(*,*) "in "//pname
c
      call smd_system_get_component(sp_energy,"energy",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no energy ',0,0)

      tag = "energy"
      call smd_data_get_index(sp_energy,tag,i_e,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      tag = "oenergy"
      call smd_data_get_index(sp_energy,tag,i_oe,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      do i=1,smd_energy_nc
        oexist = log_mb(i_oe+i-1)
        if(oexist)
     >     write(un,'(A16," : ", F12.6)') 
     >     smd_energy_name(i),dbl_mb(i_e+i-1)
      end do
      write(91,*) dbl_mb(i_e),dbl_mb(i_e+3),dbl_mb(i_e+2)
     
      return
      end

      subroutine smd_energy_compute()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_energy_data.fh"
c     
      integer un
c
      character*32 pname
      integer nrec
      integer i
      logical oexist,result
      logical ocase
      integer i_e
      integer i_oe
      character*32 sp_energy
      character*32 tag
c
      pname = "smd_energy_compute"
c
      write(*,*) "in "//pname
c
      call smd_system_get_component(sp_energy,"energy",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no energy ',0,0)

      tag = "energy"
      call smd_data_get_index(sp_energy,tag,i_e,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      tag = "oenergy"
      call smd_data_get_index(sp_energy,tag,i_oe,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      dbl_mb(i_e+3) = 0.0d0
      do i=5,8
          dbl_mb(i_e+3) = dbl_mb(i_e+3)+ dbl_mb(i_e+i-1)
      end do
      log_mb(i_oe+3) = .true.

      dbl_mb(i_e) = 0.0d0
      do i=2,4
          dbl_mb(i_e) = dbl_mb(i_e)+ dbl_mb(i_e+i-1)
      end do
      log_mb(i_oe) = .true.
      
      return
      end

c $Id$
