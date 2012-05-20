      subroutine smd_bond_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_bond,sp_param
      character*32 tag,pname
      logical result

      pname = "smd_bond_init_system"
c
      tag = "bond"
      call smd_system_get_component(sp_bond,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "parameters"
      call smd_system_get_component(sp_param,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)


      call smd_bond_init(sp_bond,sp_param,result)

      if(.not.result) then
       tag = "bond"
       call smd_system_unset_component(tag)
      end if
c
      return
      end

      subroutine smd_bond_init(sp_bond,sp_param,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_bond
      character*(*) sp_param
      logical result
c
      character*32 pname
      character*80 tag
      character*255 filename
      integer i,nb,nt
      integer i_i1,i_i2,i_d,i_k,i_it
      integer i_ta,i_tr
c
      pname = "smd_bond_init"
c
c      write(*,*) "in "//pname
c
c     get total number of shake parameters
c     ------------------------------------
      call smd_paramfile_nbond(nb)
c      write(*,*) "number of shake parameters",nb
      if(nb.eq.0) then
        result = .false.
        goto 100
      end if
c
c     create parameter data structures
c     ---------------------------
      call smd_namespace_create(sp_bond)
      tag = "bond:i1"
      call smd_data_create_get(sp_bond,tag,nb,MT_INT,i_i1)
      tag = "bond:i2"
      call smd_data_create_get(sp_bond,tag,nb,MT_INT,i_i2)
      tag = "bond:type"
      call smd_data_create_get(sp_bond,tag,nb,MT_INT,i_it)
      tag="bond:distance"
      call smd_data_create_get(sp_bond,tag,nb,MT_DBL,i_d)
      tag = "bond:strength"
      call smd_data_create_get(sp_bond,"bond:strength",nb,MT_DBL,i_k)
c
c     get array of type names
c     -----------------------
      tag = "param:atomname"
      call smd_get_ind_dim(tag,i_ta,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

      tag = "param:resname"
      call smd_get_ind(tag,i_tr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
 

      call smd_paramfile_read_bond(nb,nt,
     +                         int_mb(i_i1),
     +                         int_mb(i_i2),
     +                         int_mb(i_it),
     +                         dbl_mb(i_d),
     +                         dbl_mb(i_k),
     +                         byte_mb(i_ta),
     +                         byte_mb(i_tr))


      do i=1,nb
       write(*,*) "BOND",
     >            int_mb(i_i1+i-1),
     >            int_mb(i_i2+i-1),
     >            int_mb(i_it+i-1),
     >            dbl_mb(i_d+i-1),
     >            dbl_mb(i_k+i-1)
      end do
100   continue
      return
      end

      subroutine smd_bond_nshake(sp_bond,ns)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
c     
      integer ns
      character*(*) sp_bond
c
      character*30 tag
      character*32 pname
      integer i,nb,nt
      integer i_it
      logical result
c
      pname = "smd_bond_nshake"
c
      write(*,*) "in "//pname
c
c     get bond arrays
c     ----------------
      tag = "bond:type"
      call smd_get_ind_dim(tag,i_it,nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)


      ns = 0 
      do i=1,nb
       if(int_mb(i_it+i-1).eq.1)
     >    ns = ns + 1
 
      end do
100   continue
      return
      end
c $Id$
