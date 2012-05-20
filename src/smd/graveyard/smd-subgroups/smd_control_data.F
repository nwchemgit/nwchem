      block data initial_smd_control_data
      implicit none
#include "smd_control_data.fh"
c     
      data smd_control_fstatus /smd_control_nf*.false./

      data smd_control_fname /"rcut",
     c                       "rcut_verlet",
     c                       "step",
     c                       "temp"/

      data smd_control_istatus /smd_control_ni*.false./

      data smd_control_iname /"ndata",
     c                       "nequil",
     c                       "nprint",
     c                       "temp"/


      end

      subroutine smd_control_set_component(aname,sp_name)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_control_data.fh"
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
      pname = "smd_control_set_component"
c
      write(*,*) "in "//pname
c
c     case sensitive
      ocase = .true.
c
      nrec = smd_control_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_control_rname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      smd_control_cname(i) = sp_name
      
      return
      end

      subroutine smd_control_get_component(sp_name,aname,oexist)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_control_data.fh"
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
      pname = "smd_control_get_component"
c
      write(*,*) "in "//pname
c
c     case sensitive
      ocase = .true.
c
      nrec = smd_control_nc
c
      tag = aname
      result = inp_match(nrec,ocase,aname,smd_control_rname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//tag,0, RTDB_ERR)

      sp_name = smd_control_cname(i)
      oexist = sp_name.ne." "
      
      return
      end

      subroutine smd_control_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "smd_control_data.fh"
c     
      integer un
c
      integer i
c
      character*30 pname
      pname = "smd_control_print"
c
      write(*,*) "in "//pname
      do i=1,smd_control_nc
        write(un,'(A16,"=",A16)') smd_control_rname(i),
     >                 smd_control_cname(i)
      end do
      
      return
      end
c $Id$
