      subroutine smd_global_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_global
      character*32 tag,pname
      logical result

      pname = "smd_global_init_system"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call ga_sync()
      call smd_global_init(sp_global)
      write(*,*) "out of "//pname
      call ga_sync()
c
      return
      end

      subroutine smd_global_init(sp_global)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_global
c
      character*32 pname,tag,sp_coords
      integer na,np
      integer i_c
      logical result
c
      integer ndims
      integer ga_coords,ga_grads
c
      pname = "smd_global_init"
c
c      write(*,*) "in "//pname
c
c     get total number of atoms 
c     -------------------------
      call smd_atom_ntot(na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
      tag = "coordinates"
      call smd_system_get_component(sp_coords,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c
c     create data structures
c     ---------------------------
      call smd_namespace_create(sp_global)
c     coordinates
      if(.not.nga_create(mt_dbl, 1, 3*na, 
     >    'coordinates',1, ga_coords)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
      call smd_data_put_int(sp_global,"coordinates",1,ga_coords)
c     gradients
      if(.not.nga_create(mt_dbl, 1, 3*na, 
     >    'gradients',1, ga_grads)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
      call smd_data_put_int(sp_global,"gradients",1,ga_grads)
c
      call ga_sync()
c
c     get coordinate information
c     --------------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      if(ga_nodeid().eq.0)
     >  call nga_put(ga_coords,1, np, dbl_mb(i_c), 1)
      call ga_sync()
      call ga_print_distribution(ga_coords)
      call ga_print(ga_coords)


     
      return
      end
