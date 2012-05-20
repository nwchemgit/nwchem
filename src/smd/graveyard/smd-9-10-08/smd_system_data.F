      block data initial_smd_system_data
      implicit none
#include "smd_system_data.fh"
c     
      data smd_system_cname /smd_system_nc*" "/

      data smd_system_rname /"atom",
     c                       "parameters",
     c                       "lattice",
     c                       "type",
     c                       "potential",
     c                       "coordinates",
     c                       "charge",
     c                       "fragment",
     c                       "velocity",
     c                       "mass",
     c                       "energy",
     c                       "force",
     c                       "excl_list",
     c                       "verlet_list",
     c                       "shakelist",
     c                       "bondlist",
     c                       "cutoff",
     c                       "temperature",
     c                       "rtdb",
     c                       "ewald",
     c                       "bond"/


      end

      subroutine smd_system_set_component(aname,sp_name)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_system_data.fh"
c     
      character*(*) aname
      character*(*) sp_name
c
      character*32 pname,tag
      integer nrec
      integer i
      logical result
      logical ocase
c
      pname = "smd_system_set_component"
c
c      write(*,*) "in "//pname
c
c     case sensitive
      ocase = .true.
c
      nrec = smd_system_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_system_rname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      smd_system_cname(i) = sp_name
      
      return
      end

      subroutine smd_system_unset_component(aname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_system_data.fh"
c     
      character*(*) aname
c
      character*32 pname,tag
      integer nrec
      integer i
      logical result
      logical ocase
c
      pname = "smd_system_set_component"
c
c      write(*,*) "in "//pname
c
c     case sensitive
      ocase = .true.
c
      nrec = smd_system_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_system_rname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      smd_system_cname(i) = " "
      
      return
      end

      subroutine smd_system_get_component(sp_name,aname,oexist)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_system_data.fh"
c     
      character*(*) aname
      character*(*) sp_name
c
      character*32 pname,tag
      integer nrec
      integer i
      logical oexist,result
      logical ocase
c
      pname = "smd_system_get_component"
c
c      write(*,*) "in "//pname
c
c     case sensitive
      ocase = .true.
c
      nrec = smd_system_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_system_rname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      sp_name = smd_system_cname(i)
      oexist = sp_name.ne." "
      
      return
      end

      subroutine smd_system_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_system_data.fh"
c     
      integer un
c
      integer i
c
      character*30 pname
      pname = "smd_system_print"
c
      write(*,*) "in "//pname
      do i=1,smd_system_nc
        if(smd_system_cname(i).ne."") then
        write(un,'(A16,"= ",A16)') smd_system_rname(i),
     >                 smd_system_cname(i)
        end if
      end do
      
      return
      end
c $Id$
