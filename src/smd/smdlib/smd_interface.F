      block data initial_smd_interface_data
      implicit none
#include "smd_interface_data.fh"
c     
      data smd_interface_cname /"charge",
     c                       "param:atomname",
     c                       "param:typename",
     c                       "param:charge",
     c                       "param:mass",
     c                       "param:resname",
     c                       "coords",
     c                       "type:id",
     c                       "atom:iconst",
     c                       "atom:name",
     c                       "atom:resname",
     c                       "atom:formal_charge",
     c                       "atom:resid",    
     c                       "atom:res_type",    
     c                       "atom:res_size",
     c                       "atom:res_pointer",
     c                       "nres",
     c                       "res_dipoles",
     c                       "exlist:pointer",
     c                       "exlist:list",
     c                       "vlist:pointer",
     c                       "vlist:list",
     c                       "vlist:distances",
     c                       "vlist:displacement",
     c                       "bond:i1",
     c                       "bond:i2",
     c                       "bond:type",
     c                       "bond:distance",
     c                       "bond:strength",
     c                       "shake:i1",
     c                       "shake:i2",
     c                       "shake:distance",
     c                       "pot:p1",
     c                       "pot:p2",
     c                       "pot:p3",
     c                       "pot:type",
     c                       "force",
     c                       "lat:cell",
     c                       "lat:rcell",
     c                       "lat:fconst",
     c                       "mass",
     c                       "rcut",
     c                       "rcut_verlet",
     c                       "target_temperature",
     c                       "vel"
     c                       /

      data smd_interface_rname /"charge",
     c                       "parameters",  
     c                       "parameters",
     c                       "parameters",
     c                       "parameters",
     c                       "parameters",
     c                       "coordinates",
     c                       "type",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "atom",
     c                       "excl_list",
     c                       "excl_list",
     c                       "verlet_list",
     c                       "verlet_list",
     c                       "verlet_list",
     c                       "verlet_list",
     c                       "bond",
     c                       "bond",
     c                       "bond",
     c                       "bond",
     c                       "bond",
     c                       "shakelist",
     c                       "shakelist",
     c                       "shakelist",
     c                       "potential",
     c                       "potential",
     c                       "potential",
     c                       "potential",
     c                       "force",
     c                       "lattice",
     c                       "lattice",
     c                       "lattice",
     c                       "mass",
     c                       "cutoff",
     c                       "cutoff",
     c                       "temperature",
     c                       "velocity"
     c                       /


      end
      subroutine smd_get_ind_dim(atag,ind,n,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_interface_data.fh"
      character*(*) atag
      integer ind,n
c     
      character*32 stag,sp
      character*32 pname,buffer
      logical result,ocase
      integer i,nrec

      pname = "smd_get_ind_dim"
c
      result = .false.
c      if(atag.eq."charge") then
c        stag = "charge"
c      else
c      else
c        buffer = atag
c        call smd_message_post(pname//'tag unknown '//buffer,0,0)
c        goto 200
c      end if

c     case sensitive
      ocase = .true.

      nrec = smd_interface_nc
      buffer = atag
      result = inp_match(nrec,ocase,atag,smd_interface_cname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//buffer,0, RTDB_ERR)

      stag = smd_interface_rname(i) 
 
      call smd_system_get_component(sp,stag,result)
c
      call smd_data_get_ind_dim(sp,atag,ind,n,result)
      if(.not.result) then
        call smd_message_post(pname//'error getting index for'//buffer) 
        goto 200
      end if

200   continue
      return
      end

      subroutine smd_get_size(atag,n,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_interface_data.fh"
      character*(*) atag
      integer n
c     
      character*32 stag,sp
      character*32 pname,buffer
      logical result,ocase
      integer i,nrec

      pname = "smd_get_size"
c
      result = .false.
c     case sensitive
      ocase = .true.

      nrec = smd_interface_nc
      buffer = atag
      result = inp_match(nrec,ocase,atag,smd_interface_cname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//buffer,0, RTDB_ERR)

      stag = smd_interface_rname(i) 
 
      call smd_system_get_component(sp,stag,result)
c
      call smd_data_get_size(sp,atag,n,result)
      if(.not.result) then
        call smd_message_post(pname//'error getting index for'//buffer) 
        goto 200
      end if

200   continue
      return
      end

      subroutine smd_get_ind(atag,ind,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "smd_interface_data.fh"
      character*(*) atag
      integer ind
c     
      character*32 stag,sp
      character*32 pname,buffer
      logical result,ocase
      integer i,nrec

      pname = "smd_get_ind"
c
      result = .false.
c     case sensitive
      ocase = .true.

      nrec = smd_interface_nc
      buffer = atag
      result = inp_match(nrec,ocase,atag,smd_interface_cname,i)
      if(.not. result) 
     >  call errquit(
     >       pname//'no component '//buffer,0, RTDB_ERR)

      stag = smd_interface_rname(i) 
 
      call smd_system_get_component(sp,stag,result)
c
      call smd_data_get_index(sp,atag,ind,result)
      if(.not.result) then
        write(*,*) "WARNING"
        write(*,*) "in "//pname
        write(*,*) "looking for "//atag
        write(*,*) "in "//sp
        call smd_message_post(pname//'error getting index for'//buffer) 
        goto 200
      end if

200   continue
      return
      end
c $Id$
