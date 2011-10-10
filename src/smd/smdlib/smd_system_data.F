      subroutine smd_comp_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      character*(smd_string_size) namespace

      namespace = "component"
      call smd_namespace_create(namespace)
      return
      end

      subroutine smd_system_set_component(aname,sp_name)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      character*(*) aname
      character*(*) sp_name
c
      character*(smd_string_size) namespace
      character*(smd_string_size) pname,tag
      integer nrec
      integer i
      logical result
      logical ocase
c
      pname = "smd_system_set_component"
c
      namespace = "component"
c
c      write(*,*) "in "//pname
c
      call smd_data_put_string(namespace,aname,1,sp_name)
c
      return
      end

      subroutine smd_system_get_component(aname,sp_name,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      character*(*) aname
      character*(*) sp_name
      logical result
c
      character*(smd_string_size) namespace
      character*(smd_string_size) pname,tag
      integer nrec
      integer i
      logical ocase
c
      pname = "smd_system_set_component"
c
      namespace = "component"
c
c      write(*,*) "in "//pname
c
      call smd_data_get_string(namespace,sp_name,1,aname,result)
c
      return
      end

      subroutine smd_system_unset_component(aname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_const_data.fh"
c     
      character*(*) aname
      logical result
c
      character*(smd_string_size) namespace
      character*(smd_string_size) pname,tag
      integer nrec
      integer i
      logical ocase
c
      pname = "smd_system_set_component"
c
      namespace = "component"
c
c      write(*,*) "in "//pname
c
      call smd_data_destroy(namespace,aname,result)
c
      return
      end

      subroutine smd_system_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
c     
      integer un
c
      integer i
c
      character*30 pname
      pname = "smd_system_print"

      call smd_data_namespace_print("component",un)
c
      
      return
      end
c $Id$
