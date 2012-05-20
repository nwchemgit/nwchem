      subroutine smd_geom_monomer_create(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "bq.fh"
c
      integer im
c     
      character*32 sp_coords,sp_atom,sp_charge,sp_mass
      character*32 tag,pname
      logical result
      integer i_ir,i_ta,i_q,i_c,i_m,i_qf
      integer i_ct,h_ct
      integer i_cd,h_cd
      integer i,j
      integer nt,nm,i0
      integer rtdb
      integer igeom
      integer bq_handle
      integer nbq,i_qbq,i_cbq
      double precision scale
      double precision cm(3)
c
      logical geom_cart_set1
      external geom_cart_set1
      pname = "smd_geom_monomer_create"
c
      call smd_rtdb_get_handle(rtdb)
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
c
      tag = "charge"
      call smd_system_get_component(sp_charge,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "mass"
      call smd_system_get_component(sp_mass,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c     
      tag = "coords"
      call smd_data_get_index(sp_coords,tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)

      tag = "charge"
      call smd_data_get_index(sp_charge,tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:name"
      call smd_data_get_index(sp_atom,tag,i_ta,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:formal_charge"
      call smd_data_get_index(sp_atom,tag,i_qf,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_data_get_index(sp_atom,tag,i_ir,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      call smd_data_get_size(sp_atom,tag,nt,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c
      tag = "mass"
      call smd_data_get_index(sp_mass,tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c    
      call smd_geom_monomer_index_size(im,nt,int_mb(i_ir),nm,i0)
c
c      write(*,*) "monomer size",i0,nm
c
      if(.not.ma_push_get(mt_dbl,3*nm,'tmp',h_cd,i_cd))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nm, MA_ERR)
 
c
      if(.not.ma_push_get(mt_dbl,3*nt,'tmp',h_ct,i_ct))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nm, MA_ERR)
 
c
c     shift box so that monomer is in the center
c     -----------------------------------------
      do i=1,3*nt
        dbl_mb(i_ct+i-1) = dbl_mb(i_c+i-1)
      end do

      call smd_geom_recenter(nt,nm,i0,dbl_mb(i_cd),
     >                       dbl_mb(i_ct))

      if(.not.geom_create(igeom,"geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

     
      if(.not.geom_cart_set1(igeom,i0,nm,nt,byte_mb(i_ta),
     + dbl_mb(i_ct),dbl_mb(i_qf)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)

      if(.not.geom_masses_set(igeom,nm,dbl_mb(i_m+i0-1)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(igeom)
c
c      if (.not. geom_print(igeom))
c     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
c
      if(.not.geom_rtdb_store(rtdb,igeom,"geometry"))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c      call smd_coords_print(55)
c
c     now do the charges
c
      nbq = nt - nm
      if(.not.bq_create("monomer charges",bq_handle))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)

      if(.not.bq_alloc(bq_handle,nbq))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)


      if(.not.bq_index_coord(bq_handle,i_cbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      if(.not.bq_index_charge(bq_handle,i_qbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      scale = 1.88972598858d0

      do i=1,i0-1
         do j=1,3
            dbl_mb(i_cbq+3*(i-1)+j-1)=dbl_mb(i_ct+nt*(j-1)+i-1)*scale
         end do
         dbl_mb(i_qbq+i-1)=dbl_mb(i_q+i-1)
      end do

      do i=i0,nbq
         do j=1,3
            dbl_mb(i_cbq+3*(i-1)+j-1)=dbl_mb(i_ct+nt*(j-1)+i-1+nm)
     >                                *scale
         end do
         dbl_mb(i_qbq+i-1)=dbl_mb(i_q+nm+i-1)
      end do

c      call bq_print_info(bq_handle)

      if(.not.bq_activate(bq_handle))
     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)
c


      if(.not.ma_pop_stack(h_ct))
     & call errquit(pname//'Failed to deallocate stack h_l',nm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_cd))
     & call errquit(pname//'Failed to deallocate stack h_l',nm,
     &       MA_ERR)
c


      return
      end
c
      subroutine smd_geom_monomer_end()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "bq.fh"
c
      integer bq_handle

      character*32 pname

      pname = "smd_geom_end"
      if(.not.bq_get_active(bq_handle))
     + call errquit(pname//'Failed bq_get_active',0,CALC_ERR)

      if(.not.bq_destroy(bq_handle))
     + call errquit(pname//'Failed bq_destroy',0,CALC_ERR)


      return
      end
c
      subroutine smd_geom_recenter(nt,nm,i0,cm,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer nt
      integer nm
      integer i0
      double precision cm(nm,3)
      double precision c(nt,3)
c
      integer i,j
      double precision cshift(3)
      double precision com(3)
c
c     assemble a complete residue
c     ---------------------------
      do j=1,3
         cshift(j)=c(i0,j)
      end do

      do i=1,nm
      do j=1,3
         cm(i,j)=c(i0+i-1,j)-cshift(j)
      end do
      end do
 
      call smd_lat_rebox(nm,cm) 

      do i=1,nm
      do j=1,3
         cm(i,j)=cm(i,j)+cshift(j)
      end do
      end do
c
c     shift the entire box so that
c     monomer geom. center is at the 
c     center
c     ------------------------------ 
      com(1) = 0.0d0
      com(2) = 0.0d0
      com(3) = 0.0d0
      
      do i=1,nm
      do j=1,3
         com(j)=com(j)+cm(i,j)
      end do
      end do

      com(1) = com(1)/3.0d0 
      com(2) = com(2)/3.0d0 
      com(3) = com(3)/3.0d0 

      do i=1,nt
      do j=1,3
         c(i,j)=c(i,j)-com(j)
      end do
      end do
 
      call smd_lat_rebox(nt,c) 

c
      return
c
  
      end
c
      subroutine smd_geom_push_charges(im)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
c
      integer im
c     
      character*32 sp_coords,sp_atom,sp_charge,sp_mass
      character*32 tag,pname
      logical result
      integer i_ir,i_ta,i_q,i_c,i_m,i_qf
      integer i_ct,h_ct
      integer i,j
      integer nt,nm,i0,nqm
      integer rtdb
      integer igeom
c
      logical geom_cart_set1
      external geom_cart_set1
      logical geom_rtdb_ncent
      external geom_rtdb_ncent
c
      pname = "smd_geom_push_charges"
c
      call smd_rtdb_get_handle(rtdb)
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
c
      tag = "charge"
      call smd_system_get_component(sp_charge,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "mass"
      call smd_system_get_component(sp_mass,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

c     
      tag = "coords"
      call smd_data_get_index(sp_coords,tag,i_c,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "charge"
      call smd_data_get_index(sp_charge,tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:name"
      call smd_data_get_index(sp_atom,tag,i_ta,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:formal_charge"
      call smd_data_get_index(sp_atom,tag,i_qf,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
      tag = "atom:resid"
      call smd_data_get_index(sp_atom,tag,i_ir,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      call smd_data_get_size(sp_atom,tag,nt,result)
      if(.not. result)
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c
      tag = "mass"
      call smd_data_get_index(sp_mass,tag,i_m,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)

c    
      call smd_geom_monomer_index_size(im,nt,int_mb(i_ir),nm,i0)
c
      if(.not.ma_push_get(mt_dbl,3*nm,'tmp',h_ct,i_ct))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nm, MA_ERR)

      if(.not.geom_create(igeom,"geometry"))
     + call errquit('qmmm: Failed to create geometry',0, GEOM_ERR)

c
c     pick the best available charge geometry
c     ---------------------------------------
      tag = "cresp:geometry"
      result = geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      if(.not.result) then
        tag = "resp:geometry"
        result =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if
      if(.not.result) then
        tag = "cesp:geometry"
        result =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if
      if(.not.result) then
        tag = "esp:geometry"
        result =  geom_rtdb_ncent(rtdb,tag(1:inp_strlen(tag)),nqm)
      end if

      if(.not.result) then
         call errquit(pname//'cannot find charges')
      end if

      if(nqm.ne.nm)
     + call errquit(pname//'mismatch in sizes'//tag,0, GEOM_ERR)

      if(.not.geom_rtdb_load(rtdb,igeom,tag(1:inp_strlen(tag))))
     + call errquit(pname//'Failed to load geometry'//tag,0, GEOM_ERR)

      if(.not.geom_cart_get_charges(igeom,nm,dbl_mb(i_q+i0-1)))
     + call errquit(pname//'Failed to get charges',0, GEOM_ERR)

       if(.not.geom_destroy(igeom))
     + call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)

c
      if(.not.ma_pop_stack(h_ct))
     & call errquit(pname//'Failed to deallocate stack h_l',nm,
     &       MA_ERR)
c
      return
      end
c
      subroutine smd_geom_monomer_index_size(im,nt,ir,nm,i0)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer im
      integer rtdb
      integer nt
      integer ir(nt)
      integer nm
      integer i0
c
      integer i,j
c
c     count number of atoms in residue im
c     contigious storage is presumed
c     -----------------------------------
      nm = 0
      do i=1,nt
       j = i
       if(im.eq.ir(i)) then
100     continue
        nm = nm + 1
        j = j + 1
        if(im.eq.ir(j)) then
         goto 100
        else
         goto 101 
        end if
       end if
      end do
101   continue
      i0=i
      return
c
  
      end

c $Id$
