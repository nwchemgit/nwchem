      subroutine smd_shakelist_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_bond,sp_type,sp_shakelist
      character*32 tag,pname
      logical result

      pname = "smd_shakelist_init_system"
c
      tag = "bond"
      call smd_system_get_component(sp_bond,tag,result)
      if(.not.result) goto 200

      tag = "type"
      call smd_system_get_component(sp_type,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "shakelist"
      call smd_system_get_component(sp_shakelist,tag,result)
      if(.not.result) 
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_shakelist_init(sp_shakelist,sp_bond,sp_type,result)

200   continue
      if(.not.result) then
       tag = "shakelist"
       call smd_system_unset_component(tag)
      end if
c
      return
      end

      subroutine smd_shakelist_init(sp_shakelist,sp_bond,sp_type,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_type
      character*(*) sp_bond
      character*(*) sp_shakelist
      integer rtdb
      logical result
c
      character*32 pname
      character*80 tag
      character*255 filename
      integer na,nb,ns
      integer i_it
      integer i_ib1,i_ib2,i_db,i_itb
      integer h_is1t,i_is1t
      integer h_is2t,i_is2t
      integer h_dst,i_dst
      integer i_is1
      integer i_is2
      integer i_ds
      integer i
c
      pname = "smd_shakelist_init"
c
c      write(*,*) "in "//pname
c
c     get array of types
c     ------------------
      tag = "type:id"
      call smd_data_get_index(sp_type,tag,i_it,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      call smd_data_get_size(sp_type,tag,na,result)
c
c     get bond arrays
c     ---------------
      tag = "bond:i1"
      call smd_data_get_index(sp_bond,tag,i_ib1,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "bond:i2"
      call smd_data_get_index(sp_bond,tag,i_ib2,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "bond:distance"
      call smd_data_get_index(sp_bond,tag,i_db,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "bond:type"
      call smd_data_get_index(sp_bond,tag,i_itb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)
      call smd_data_get_size(sp_bond,tag,nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)


c
c     allocate initial storage for shake list
c     ---------------------------------------
      ns = na
      if(.not.ma_push_get(mt_int,ns,'tmp i1',h_is1t,i_is1t))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_int,ns,'tmp i2',h_is2t,i_is2t))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'tmp d',h_dst,i_dst))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)


      call smd_shakelist_set(ns,nb,na,
     >                       int_mb(i_is1t),
     >                       int_mb(i_is2t),
     >                       dbl_mb(i_dst),
     >                       int_mb(i_ib1),
     >                       int_mb(i_ib2),
     >                       int_mb(i_itb),
     >                       dbl_mb(i_db),
     >                       int_mb(i_it))
c
c     create shake list structure
c     ---------------------------
      if(ns.eq.0) then
        result = .false.
        goto 200
      end if
      call smd_namespace_create(sp_shakelist)
      call smd_data_create(sp_shakelist,"shake:i1",ns,MT_INT)
      call smd_data_create(sp_shakelist,"shake:i2",ns,MT_INT)
      call smd_data_create(sp_shakelist,"shake:distance",ns,MT_DBL)
c
c     get shake arrays
c     ----------------
      tag = "shake:i1"
      call smd_data_get_index(sp_shakelist,tag,i_is1,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "shake:i2"
      call smd_data_get_index(sp_shakelist,tag,i_is2,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "shake:distance"
      call smd_data_get_index(sp_shakelist,tag,i_ds,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      do i=1,ns
       int_mb(i_is1+i-1) = int_mb(i_is1t+i-1)
       int_mb(i_is2+i-1) = int_mb(i_is2t+i-1)
       dbl_mb(i_ds+i-1)  = dbl_mb(i_dst+i-1)
      end do

200   continue
      if(.not.ma_pop_stack(h_dst))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_is2t))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_is1t))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      return
      end

      subroutine smd_shakelist_set(ns,nb,na,
     >                       is1,
     >                       is2,
     >                       ds,
     >                       ib1,
     >                       ib2,
     >                       itb,
     >                       db,
     >                       it)
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      integer ns,nb,na
      integer is1(ns)
      integer is2(ns)
      double precision  ds(ns)
      integer ib1(nb)
      integer ib2(nb)
      integer itb(nb)
      double precision  db(ns)
      integer it(na)
c
      integer i,i1,i2,j,nlist
c
      nlist = 0
      do i=1,nb
       if(itb(i).eq.1) then
        i1=0
        i2=0
        do j=1,na
         if(it(j).eq.ib1(i)) i1=j
         if(it(j).eq.ib2(i)) i2=j
         if(i1*i2.ne.0) then
          nlist = nlist + 1
          is1(nlist) = min(i1,i2)
          is2(nlist) = max(i1,i2)
          ds(nlist)  = db(i)
          i1=0
          i2=0
         end if
        end do
       end if
      end do
      ns = nlist

      return
      end

c $Id$
