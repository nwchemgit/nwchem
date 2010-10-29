*
* $Id$
*
      subroutine smd_pot_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_pot,sp_param
      character*32 tag,pname
      logical result

      pname = "smd_pot_init_system"
c
      tag = "potential"
      call smd_system_get_component(sp_pot,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_pot_init(sp_pot,sp_param,result)

      if(.not.result) then
       tag = "potential"
       call smd_system_unset_component(tag)
      end if
c
      return
      end

      subroutine smd_pot_init(sp_pot,sp_param,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_pot
      character*(*) sp_param
      logical result
c
      character*32 pname
      character*80 tag
      character*255 filename
      integer nt,np
      integer i_p1,i_p2,i_p3
      integer i_pt
      integer i_it
      integer i_tt,i_tr
c
      pname = "smd_pot_init"
c
c      write(*,*) "in "//pname
c
c     get total number of parameters
c     ------------------------------
      call smd_param_ntot(sp_param,nt)
c
c     total number of possible pairwise interactions
c     ----------------------------------------------
      np=nt*(nt+1)
      np=np/2
c
c     create pot data structures
c     ---------------------------
      call smd_namespace_create(sp_pot)
      call smd_data_create(sp_pot,"pot:p1",np,MT_DBL)
      call smd_data_create(sp_pot,"pot:p2",np,MT_DBL)
      call smd_data_create(sp_pot,"pot:p3",np,MT_DBL)
      call smd_data_create(sp_pot,"pot:type",np,MT_INT)
c
c     fill in parameters from parameter file
c     -------------------------------------------
      tag = "pot:p1"
      call smd_data_get_index(sp_pot,tag,i_p1,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "pot:p2"
      call smd_data_get_index(sp_pot,tag,i_p2,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "pot:p3"
      call smd_data_get_index(sp_pot,tag,i_p3,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "pot:type"
      call smd_data_get_index(sp_pot,tag,i_pt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "param:typename"
      call smd_data_get_index(sp_param,tag,i_tt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      tag = "param:resname"
      call smd_data_get_index(sp_param,tag,i_tr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)


      call smd_paramfile_read_pot(np,
     +                         int_mb(i_pt),
     +                         dbl_mb(i_p1),
     +                         dbl_mb(i_p2),
     +                         dbl_mb(i_p3),
     +                         nt,
     +                         byte_mb(i_tt),
     +                         byte_mb(i_tr))


      if(np.eq.0) then
       result = .false.
       call smd_namespace_destroy(sp_pot)
      end if
      return
      end

      subroutine smd_pot_compute()

      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_const_data.fh"
      double precision  e
c
      character*30 pname
      character*32 sp_pot
      character*32 sp_lattice
      character*32 sp_force
      character*32 sp_vlist
      character*32 sp_type
      character*32 sp
      character*32 tag
      integer i_p1,i_p2,i_p3,i_pt
      integer nt,np
      integer i_f,i_c,i_rlat
      integer i_lat
      integer i_l,i_p
      integer i_it
      integer nl,na
      double precision vol
      double precision alpha
      double precision rcutsq
      integer kmax(3)
      integer ntmp,nk
      logical result

      pname = "smd_pot_compute"


      e = 0.0d0
c
c     get components
c     --------------
      call smd_system_get_component(sp_pot,"potential",result)
      if(.not.result) return

      call smd_system_get_component(sp_force,"force",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no force ',0,0)

      call smd_system_get_component(sp_vlist,"verlet_list",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no verlet_list ',0,0)

      call smd_system_get_component(sp_type,"type",result)
      if(.not.result)
     >  call errquit(
     >       pname//'no types ',0,0)


c
c     get verlet list data
c     --------------------
      sp = sp_vlist
      tag = "vlist:pointer"
      call smd_data_get_index(sp,tag,i_p,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "vlist:list"
      call smd_data_get_index(sp,tag,i_l,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      call smd_data_get_size(sp,tag,nl,result)

      tag = "vlist:distances"
      call smd_data_get_index(sp_vlist,tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c
c     get type data
c     ------------
      sp = sp_type
      tag = "type:id"
      call smd_data_get_index(sp,tag,i_it,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      call smd_data_get_size(sp,tag,nt,result)
c
c     get pot data
c     ------------
      sp = sp_pot
      tag = "pot:p1"
      call smd_data_get_index(sp,tag,i_p1,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      call smd_data_get_size(sp,tag,np,result)
      tag = "pot:p2"
      call smd_data_get_index(sp,tag,i_p2,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      tag = "pot:p3"
      call smd_data_get_index(sp,tag,i_p3,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      tag = "pot:type"
      call smd_data_get_index(sp,tag,i_pt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
c     get force data
c     --------------
      tag = "force"
      sp = sp_force
      call smd_data_get_index(sp,tag,i_f,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      call smd_data_get_size(sp,tag,na,result)
      na = na/3
c
c     get cutoff data
c     ---------------
      call smd_cutoff_get_rcut(rcutsq)
      rcutsq = rcutsq*rcutsq

      call       smd_pot_compute_generic(np,
     >                                  nl,
     >                                  na,
     >                                  rcutsq,
     >                                  int_mb(i_it),
     >                                  dbl_mb(i_p1),
     >                                  dbl_mb(i_p2),
     >                                  dbl_mb(i_p3),
     >                                  int_mb(i_pt),
     >                                  dbl_mb(i_c),
     >                                  dbl_mb(i_f),
     >                                  int_mb(i_p),
     >                                  int_mb(i_l),
     >                                  e)

      call smd_energy_set_component("vdw",e)
      return

      END



