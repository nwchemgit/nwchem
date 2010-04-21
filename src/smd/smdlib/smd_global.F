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
      call smd_global_update("coords","put")
c      call smd_global_update_res("coords","put",2)

c      call ga_print_distribution(ga_coords)
c      call ga_print(ga_coords)
c
      call smd_global_pull("gradients")
c      call ga_print_distribution(ga_grads)
c      call ga_print(ga_grads)


     
      return
      end

      subroutine smd_global_pull_res_coords(ires)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer ires
c
      character*32 sp_global
      character*32 pname,tag,sp_coords
      integer na,np
      integer i_c,i_ir
      integer nm,i0,nt
      integer i,k
      logical result
      integer ilo,ihi
c
      integer ndims
      integer ga_coords,ga_grads
c
      pname = "smd_global_init"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "coordinates"
      call smd_system_get_component(sp_coords,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c
c     get coordinate information
c     --------------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_get_ind_dim(tag,i_ir,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

       call smd_geom_monomer_index_size(ires,nt,int_mb(i_ir),nm,i0)

c       write(*,*) "coords of monomer",ires,nm
c       do i=1,3*nm
c         write(*,*) (dbl_mb(i_c+3*(i0-1)+i-1))
c       end do
c       write(*,*) "coords of monomer done"
      call smd_data_get_int(sp_global,"coordinates",1,ga_coords,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)

c      if(ga_nodeid().eq.0) then
c         call ga_summarize(.false.)
c         write(*,*) "ga coord handle",ga_coords
c         call util_flush(6)
c      end if
c      call ga_zero(ga_coords)
      if(ga_nodeid().eq.0) then
        ilo = 3*(i0-1)+1
        ihi = 3*(i0-1)+3*nm
        call nga_put(ga_coords,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
      end if
      call ga_sync()
c      write(*,*) "coords after pushing monomer",ires
c      call ga_print_distribution(ga_coords)
c      call ga_print(ga_coords)
     
      return
      end

      subroutine smd_global_pull_res(atag,ires)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) atag
      integer ires
c
      character*32 sp_global
      character*32 pname,tag,ltag,gtag
      integer na,np
      integer i_c,i_ir
      integer nm,i0,nt
      integer i,k
      logical result
      integer ilo,ihi
c
      integer ndims
      integer ga_handle
c
      pname = "smd_global_pull_res"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "atom:resid"
      call smd_get_ind_dim(tag,i_ir,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      call smd_geom_monomer_index_size(ires,nt,int_mb(i_ir),nm,i0)
c
      if(atag.eq."coords") then
        gtag = "coordinates"
        ltag = "coords"
      else
        call errquit(
     >       pname//'unknown tag',0, 0)

      end if
c
c     get local array
c     ---------------
      call smd_get_ind_dim(ltag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c     get global array
c     ---------------
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)
c
c     push data onto global array
c     ---------------------------
      if(ga_nodeid().eq.0) then
        ilo = 3*(i0-1)+1
        ihi = 3*(i0-1)+3*nm
        call nga_put(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
      end if
      call ga_sync()
     
      return
      end

      subroutine smd_global_update_res(atag,oper,ires)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) atag
      character*(*) oper
      integer ires
c
      character*32 sp_global
      character*32 pname,tag,ltag,gtag
      integer na,np
      integer i_c,i_ir
      integer nm,i0,nt
      integer i,k
      logical result
      integer ilo,ihi
c
      integer ndims
      integer ga_handle
c
      pname = "smd_global_pull_res"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "atom:resid"
      call smd_get_ind_dim(tag,i_ir,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      call smd_geom_monomer_index_size(ires,nt,int_mb(i_ir),nm,i0)
c
      if(atag.eq."coords") then
        gtag = "coordinates"
        ltag = "coords"
      else if(atag.eq."gradients") then
        gtag = "gradients"
        ltag = "force"
      else
        call errquit(
     >       pname//'unknown tag',0, 0)

      end if
c
c     get local array
c     ---------------
      call smd_get_ind_dim(ltag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c     get global array
c     ---------------
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)
c
c     put/get data on/from global array
c     ---------------------------
      if(ga_nodeid().eq.0) then
        ilo = 3*(i0-1)+1
        ihi = 3*(i0-1)+3*nm
        if(oper.eq."put") then
          call nga_put(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
        else if(oper.eq."get") then
          call nga_get(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
        else
          call errquit(
     >       pname//'unknown operation',0, 0)
        end if
      end if
      call ga_sync()
     
      return
      end

      subroutine smd_global_update(atag,oper)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) atag
      character*(*) oper
c
      character*32 sp_global
      character*32 pname,tag,ltag,gtag
      integer na,np
      integer i_c,i_ir
      integer nm,i0,nt
      integer i,k
      logical result
      integer ilo,ihi
c
      integer ndims
      integer ga_handle
c
      pname = "smd_global_pull_res"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      if(atag.eq."coords") then
        gtag = "coordinates"
        ltag = "coords"
      else if(atag.eq."gradients") then
        gtag = "gradients"
        ltag = "force"
      else 
        call errquit(
     >       pname//'unknown tag',0, 0)

      end if
c
c     get local array
c     ---------------
      call smd_get_ind_dim(ltag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c     get global array
c     ---------------
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)
c
c     put/get data on/from global array
c     ---------------------------
      if(ga_nodeid().eq.0) then
        ilo = 1
        ihi = np
        if(oper.eq."put") then
          call nga_put(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
        else if(oper.eq."get") then
          call nga_get(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
        else
          call errquit(
     >       pname//'unknown operation',0, 0)
        end if
      end if
      call ga_sync()
     
 
      return
      end

      subroutine smd_global_print(atag)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) atag
c
      character*32 sp_global
      character*32 pname,tag,ltag,gtag
      integer na,np
      integer i_c,i_ir
      integer nm,i0,nt
      integer i,k
      logical result
      integer ilo,ihi
c
      integer ndims
      integer ga_handle
c
      pname = "smd_global_pull_res"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      if(atag.eq."coords") then
        gtag = "coordinates"
        ltag = "coords"
      else if(atag.eq."gradients") then
        gtag = "gradients"
        ltag = "force"
      else 
        call errquit(
     >       pname//'unknown tag',0, 0)

      end if
c
c     get global array
c     ---------------
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)

      call ga_print(ga_handle)

      call ga_sync()
     
 
      return
      end

      subroutine smd_global_pull(atag)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) atag
c
      character*32 sp_global
      character*32 pname,tag,ltag,gtag
      integer na,np
      integer i_c,i_ir
      integer nm,i0,nt
      integer i,k
      logical result
      integer ilo,ihi
c
      integer ndims
      integer ga_handle
c
      pname = "smd_global_pull_res"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      if(atag.eq."coords") then
        gtag = "coordinates"
        ltag = "coords"
      else if(atag.eq."gradients") then
        gtag = "gradients"
        ltag = "force"
      else 
        call errquit(
     >       pname//'unknown tag',0, 0)

      end if
c
c     get local array
c     ---------------
      call smd_get_ind_dim(ltag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c     get global array
c     ---------------
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)
c
c     push data onto global array
c     ---------------------------
      if(ga_nodeid().eq.0) then
        ilo = 1
        ihi = np
        call nga_put(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
      end if
      call ga_sync()
     
      return
      end

      subroutine smd_global_pull_coords()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_global
      character*32 pname,tag,sp_coords
      integer na,np
      integer i_c,i_ir
      integer nm,i0,nt
      integer i,k
      logical result
      integer ilo,ihi
c
      integer ndims
      integer ga_coords,ga_grads
c
      pname = "smd_global_init"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "coordinates"
      call smd_system_get_component(sp_coords,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c
c     get coordinate information
c     --------------------------
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,np,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      if(ga_nodeid().eq.0) then
        ilo = 1
        ihi = np
        call nga_put(ga_coords,ilo,ihi, dbl_mb(i_c), 1)
      end if
      call ga_sync()
     
      return
      end
