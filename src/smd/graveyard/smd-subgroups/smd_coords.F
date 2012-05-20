      subroutine smd_coords_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_coords,sp_atom
      character*32 tag,pname
      logical result

      pname = "smd_coords_init_system"
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "coordinates"
      call smd_system_get_component(sp_coords,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_coords_init(sp_coords,sp_atom)
      
      call smd_coords_read(sp_coords)
c
      return
      end

      subroutine smd_coords_init(sp_coords,sp_atom)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_coords
      character*(*) sp_atom
c
      character*32 pname
      integer na
c
      pname = "smd_coords_init"
c
c      write(*,*) "in "//pname
c
c     get total number of atoms 
c     -------------------------
      call smd_atom_ntot(sp_atom,na)
      if(na.le.0)
     >  call errquit(
     >       pname//'no atoms ',0, RTDB_ERR)
c
c     create coords data structures
c     ---------------------------
      call smd_namespace_create(sp_coords)
      call smd_data_create(sp_coords,"coords",3*na,MT_DBL)

      return
      end

      subroutine smd_coords_read(sp_coords)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_coords
      integer rtdb
c
      character*32 pname
      character*72 tag
      integer i_c
      integer na
      logical result
c
      pname = "smd_coords_read"
c
      write(*,*) "in "//pname
c
c     fill in coordinates from pdb file if any
c     ----------------------------------------
      tag = "coords"
      call smd_data_get_index(sp_coords,tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      call smd_data_get_size(sp_coords,tag,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      if(na.le.0)
     >call errquit(
     >     pname//'no atoms ',0, RTDB_ERR)
      na = na/3


      call smd_coordfile_read_coords(na,
     +                          dbl_mb(i_c))


      return
      end

      subroutine smd_coords_rebox()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_atom_rebox"

c
c     get atomic coordinates
c     ----------------------
      call smd_system_get_component(sp_coords,"coordinates",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coordinates ',0,0)

      tag = "coords"
      call smd_data_get_index(sp_coords,tag,i_c,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      call smd_data_get_size(sp_coords,tag,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      na = na/3

      call smd_lat_rebox(na,
     >                    dbl_mb(i_c))

      return
      end

      subroutine smd_coords_print(un)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer un
      character*32 sp_coords
c
      character*72 tag
      character*30 pname
      integer na
      integer i_c,i_lrc,i_lc
      logical result

      pname = "smd_atom_rebox"

c
c     get atomic coordinates
c     ----------------------
      call smd_system_get_component(sp_coords,"coordinates",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no coordinates ',0,0)

      tag = "coords"
      call smd_data_get_index(sp_coords,tag,i_c,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      call smd_data_get_size(sp_coords,tag,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      na = na/3

      
      call smd_util_print_force_array(un,na,
     >                           dbl_mb(i_c))  

      return
      end

      subroutine smd_coords_update()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "smd_system.fh"
c     
      character*72 sp_vel
      character*72 sp_coords
      character*72 sp_mass
      character*72 sp_force
      character*72 sp_shakelist
      character*32 pname
      character*80 tag
      integer nt,na
      integer i_f
      integer i_m,i_c,i_v
      logical result
      integer rtdb
      double precision ekin, tstep
      logical oshake
      integer i_is1,i_is2,i_ds,ns
      integer h_ncc,i_ncc
      integer h_nvv,i_nvv
      integer h_dcc,i_dcc
      integer h_nrij,i_nrij
      integer h_orij,i_orij
c
      pname = "smd_coords_update"
c
c      write(*,*) "in "//pname
c
      call  smd_rtdb_get_handle(rtdb)
c
c     get components
c     --------------
      tag = "velocity"
      call smd_system_get_component(sp_vel,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "mass"
      call smd_system_get_component(sp_mass,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "force"
      call smd_system_get_component(sp_force,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "coordinates"
      call smd_system_get_component(sp_coords,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c
c     get velocity array
c     ------------------
      tag = "vel"
      call smd_data_get_index(sp_vel,tag,i_v,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      call  smd_data_get_size(sp_vel,tag,na,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
      na = na/3

c
c     get mass array
c     ------------------
      tag = "mass"
      call smd_data_get_index(sp_mass,tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)
c
c     get force array
c     ---------------
      tag = "force"
      call smd_data_get_index(sp_force,tag,i_f,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

c
c     get coord array
c     ---------------
      tag = "coords"
      call smd_data_get_index(sp_coords,tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, 0)

c
c     get time step
c     -------------
      tag="smd:step"
        if (.not.rtdb_get(rtdb,tag,mt_dbl,1,tstep)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
 

      oshake = smd_system_shake()

      if(oshake) then

        tag = "shakelist"
        call smd_system_get_component(sp_shakelist,tag,result)
        if(.not.result)
     >    call errquit(
     >         pname//'no component '//tag,0,0)

        tag = "shake:i1"
        call smd_data_get_index(sp_shakelist,tag,i_is1,result)
        if(.not. result) 
     >    call errquit(
     >         pname//'error getting index for '//tag,0, 0)

        tag = "shake:i2"
        call smd_data_get_index(sp_shakelist,tag,i_is2,result)
        if(.not. result) 
     >    call errquit(
     >         pname//'error getting index for '//tag,0, 0)

        tag = "shake:distance"
        call smd_data_get_index(sp_shakelist,tag,i_ds,result)
        if(.not. result) 
     >    call errquit(
     >         pname//'error getting index for '//tag,0, 0)

        call smd_data_get_size(sp_shakelist,tag,ns,result)
        if(.not. result) 
     >    call errquit(
     >         pname//'error getting size for '//tag,0, 0)

        if(.not.ma_push_get(mt_dbl,na*3,'i_ncc',h_ncc,i_ncc))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,na*3,'i_nvv',h_nvv,i_nvv))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,na*3,'i_dcc',h_dcc,i_dcc))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,ns*3,'i_nrij',h_nrij,i_nrij))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

        if(.not.ma_push_get(mt_dbl,ns*3,'i_orij',h_orij,i_orij))
     +   call errquit(pname//'Failed to allocate memory',
     +   0, MA_ERR)

      end if
      if(oshake) then
            call smd_leapf_shake(na,
     >                           ns,
     >                           tstep,
     >                           ekin,
     >                           dbl_mb(i_m),
     >                           int_mb(i_is1),
     >                           int_mb(i_is2),
     >                           dbl_mb(i_ds),
     >                           dbl_mb(i_ncc),
     >                           dbl_mb(i_nvv),
     >                           dbl_mb(i_dcc),
     >                           dbl_mb(i_nrij),
     >                           dbl_mb(i_orij),
     >                           dbl_mb(i_f),
     >                           dbl_mb(i_v),
     >                           dbl_mb(i_c))


      else
      call smd_leapf(na,
     >                 tstep,
     >                 ekin,
     >                 dbl_mb(i_m),
     >                 dbl_mb(i_f),
     >                 dbl_mb(i_v),
     >                 dbl_mb(i_c))

      end if

      call smd_coords_rebox()

      if(oshake) then

        if(.not.ma_pop_stack(h_orij))
     &   call errquit(pname//'Failed to deallocate stack h_orij',0,
     &         MA_ERR)

        if(.not.ma_pop_stack(h_nrij))
     &   call errquit(pname//'Failed to deallocate stack h_nrij',0,
     &         MA_ERR)

        if(.not.ma_pop_stack(h_dcc))
     &   call errquit(pname//'Failed to deallocate stack h_dcc',0,
     &         MA_ERR)

        if(.not.ma_pop_stack(h_nvv))
     &   call errquit(pname//'Failed to deallocate stack h_nvv',0,
     &         MA_ERR)


        if(.not.ma_pop_stack(h_ncc))
     &   call errquit(pname//'Failed to deallocate stack h_ncc',0,
     &         MA_ERR)

      end if

      return
      end
c $Id$
