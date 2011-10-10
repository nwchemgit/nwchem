C>      smd global data structure stores global arrays
C>      and handles data exchange with their local counterparts
      subroutine smd_global_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "msgids.fh"
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
      character*32 pname,tag
      integer na,np
      integer nres
      integer i_c
      logical result
c
      integer ndims
      integer ga_coords,ga_grads,ga_dip,ga_charges
c
      pname = "smd_global_init"
c
c
c     get total number of atoms and residues
c     --------------------------------------
      call smd_atom_ntot(na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
      call smd_system_nres(nres)
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
c     dipoles
      if(.not.nga_create(mt_dbl, 1, 3*nres, 
     >    'residue dipoles',1, ga_dip)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
      call smd_data_put_int(sp_global,"res-dipoles",1,ga_dip)
c      call ga_print(ga_dip)
c     charges
      if(.not.nga_create(mt_dbl, 1, na, 
     >    'charges',1, ga_charges)) 
     >  call errquit(
     >       pname//'nga_create error ',0, RTDB_ERR)
      call smd_data_put_int(sp_global,"charges",1,ga_charges)
ccc
      call ga_sync()
c
c     get coordinate information
c     --------------------------
      call smd_global_update("coords","put")
      call smd_global_update("gradients","put")
      call smd_global_update("charges","put")
c
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
#include "msgids.fh"
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
      integer i_pr,i_mr
c
      integer ndims
      integer ga_handle
c
      pname = "smd_global_update_res"
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
c      
      tag = "atom:res_pointer"
      call smd_get_ind(tag,i_pr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:res_size"
      call smd_get_ind(tag,i_mr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c     ----------------------------------------------------------------
c     ndims is the size of 2nd dimension
c     which is 3 for gradients, and 1 for scalar charges
c     i0 is the beginning index of the residue ires in the atomic array
c     nm is the size (number of atoms) in the residue ires
c     ----------------------------------------------------------------
      if(atag.eq."coords") then
        ndims = 3
        i0=int_mb(i_pr+ires-1)
        nm=int_mb(i_mr+ires-1)
        gtag = "coordinates"
        ltag = "coords"
      else if(atag.eq."dipoles") then
        ndims = 3
        i0=ires
        nm=1
        gtag = "res-dipoles"
        ltag = "res_dipoles"
      else if(atag.eq."charges") then
        ndims = 1
        i0=int_mb(i_pr+ires-1)
        nm=int_mb(i_mr+ires-1)
        gtag = "charges"
        ltag = "charge"
      else if(atag.eq."gradients") then
        ndims = 3
        i0=int_mb(i_pr+ires-1)
        nm=int_mb(i_mr+ires-1)
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
        ilo = ndims*(i0-1)+1
        ihi = ndims*(i0-1)+ndims*nm
        if(oper.eq."put") then
          if(ga_nodeid().eq.0) then
            call nga_put(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
          end if
        else if(oper.eq."get") then
          if(ga_nodeid().eq.0) then
            call nga_get(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
          end if
          call ga_brdcst(msg_smd,dbl_mb(i_c+ilo-1),
     >     ndims*nm*ma_sizeof(mt_dbl,1,mt_byte),0)
        else
          call errquit(
     >       pname//'unknown operation',0, 0)
        end if
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
#include "msgids.fh"
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
      else if(atag.eq."charges") then
        gtag = "charges"
        ltag = "charge"
      else if(atag.eq."dipoles") then
        gtag = "res-dipoles"
        ltag = "res_dipoles"
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
      ilo = 1
      ihi = np
      if(oper.eq."put") then
        if(ga_nodeid().eq.0) then
          call nga_put(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
        end if
      else if(oper.eq."get") then
        if(ga_nodeid().eq.0) then
          call nga_get(ga_handle,ilo,ihi, dbl_mb(i_c+ilo-1), 1)
        end if
        call ga_brdcst(msg_smd,dbl_mb(i_c+ilo-1),
     >     np*ma_sizeof(mt_dbl,1,mt_byte),0)
      else
        call errquit(
     >     pname//'unknown operation',0, 0)
      end if
     
c      call ga_print(ga_handle) 
      return
      end

      subroutine smd_global_create_copy(atag,gh)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) atag
      integer gh
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
      pname = "smd_global_get_copy"
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
      else if(atag.eq."charges") then
        gtag = "charges"
        ltag = "charge"
      else if(atag.eq."dipoles") then
        gtag = "res-dipoles"
        ltag = "res_dipoles"
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
c
      if(.not.ga_duplicate(ga_handle,gh,"copy of "//gtag))
     >  call errquit(
     >       pname//'ga_duplicate error',0, RTDB_ERR)

      call ga_copy(ga_handle,gh)
 
      return
      end

      subroutine smd_global_diff(atag,gh,eps)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) atag
      integer gh
      double precision eps
c
      character*32 sp_global
      character*32 pname,tag,ltag,gtag
      integer na,np
      integer i_c,i_ir
      logical result
c
      integer gh0
      integer type,ndim,dims
      integer type0,ndim0,dims0
c
      pname = "smd_global_get_copy"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      if(atag.eq."coords") then
        gtag = "coordinates"
      else if(atag.eq."charges") then
        gtag = "charges"
      else if(atag.eq."dipoles") then
        gtag = "res-dipoles"
      else if(atag.eq."gradients") then
        gtag = "gradients"
      else 
        call errquit(
     >       pname//'unknown tag',0, 0)

      end if
c
c     get global array
c     ---------------
      call smd_data_get_int(sp_global,gtag,1,gh0,result)
      if(.not.result) call errquit(
     >       pname//'no coordindates ',0,0)
c
c      write(*,*) "check gh"
c      call ga_print(gh)
c      write(*,*) "check gh0"
c      call ga_print(gh0)
c
      call ga_add(1.0d0,gh,-1.0d0,gh0,gh)
      call ga_norm_infinity(gh,eps)
      call ga_copy(gh0,gh)

      call nga_inquire(gh,type,ndim,dims)
      call nga_inquire(gh0,type0,ndim0,dims0)
     
      if(type.ne.type0) then
        call errquit(pname//'incomp types ',0,0) 
      end if
      if(ndim.ne.ndim0) then
        call errquit(pname//'incomp ndim ',0,0) 
      end if
      if(dims.ne.dims0) then
        call errquit(pname//'incomp dims',0,0) 
      end if

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
      else if(atag.eq."charges") then
        gtag = "charges"
        ltag = "charge"
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

     
 
      return
      end

      subroutine smd_global_end()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
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
      pname = "smd_global_end"
c
      tag = "global"
      call smd_system_get_component(sp_global,tag,result)
      if(.not.result) return
c
c     get global array
c     ---------------
      gtag="coordinates"
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(result) then
       if(.not.ga_destroy(ga_handle)) call errquit(
     >       pname//'cannot destroy '//gtag,0,0)
      end if

      gtag="gradients"
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(result) then
       if(.not.ga_destroy(ga_handle)) call errquit(
     >       pname//'cannot destroy '//gtag,0,0)
      end if

      gtag="charges"
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(result) then
       if(.not.ga_destroy(ga_handle)) call errquit(
     >       pname//'cannot destroy '//gtag,0,0)
      end if

      gtag="residue dipoles"
      call smd_data_get_int(sp_global,gtag,1,ga_handle,result)
      if(result) then
       if(.not.ga_destroy(ga_handle)) call errquit(
     >       pname//'cannot destroy '//gtag,0,0)
      end if
      
      call smd_namespace_destroy(sp_global)
 
      return
      end

c $Id$
