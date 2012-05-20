      subroutine smd_bondlist_init_system()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*32 sp_bond,sp_type,sp_bondlist
      character*32 tag,pname
      logical result

      pname = "smd_bondlist_init_system"
c
      tag = "bond"
      call smd_system_get_component(sp_bond,tag,result)
      if(.not.result) goto 200

      tag = "type"
      call smd_system_get_component(sp_type,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      tag = "bondlist"
      call smd_system_get_component(sp_bondlist,tag,result)
      if(.not.result) 
     >  call errquit(
     >       pname//'no component '//tag,0,0)

      call smd_bondlist_init(sp_bondlist,result)

200   continue
      if(.not.result) then
       tag = "bondlist"
       call smd_system_unset_component(tag)
      end if
c
      return
      end

      subroutine smd_bondlist_init(sp_bondlist,result)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_bondlist
      integer rtdb
c
      character*32 pname
      character*80 tag
      character*255 filename
      integer na,nb,ns
      integer i_it
      integer i_kb
      integer i_ib1,i_ib2,i_db,i_itb
      integer h_il1t,i_il1t
      integer h_il2t,i_il2t
      integer h_dlt,i_dlt
      integer h_klt,i_klt
      integer h_tlt,i_tlt
      integer i_il1
      integer i_il2
      integer i_dl
      integer i_kl
      integer i_tl
      logical result
      integer i
c
      pname = "smd_bondlist_init"
c
c      write(*,*) "in "//pname
c
c     get array of types
c     ------------------
      tag = "type:id"
      call smd_get_ind_size(tag,i_it,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c     get bond arrays
c     ---------------
      tag = "bond:i1"
      call smd_get_ind(tag,i_ib1,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "bond:i2"
      call smd_get_ind(tag,i_ib2,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "bond:distance"
      call smd_get_ind(tag,i_db,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "bond:strength"
      call smd_get_ind(tag,i_kb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)


      tag = "bond:type"
      call smd_get_ind_size(tag,i_itb,nb,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

c
c     allocate initial storage for bond list
c     ---------------------------------------
      ns = na
      if(.not.ma_push_get(mt_int,ns,'tmp i1',h_il1t,i_il1t))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_int,ns,'tmp i2',h_il2t,i_il2t))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'tmp d',h_dlt,i_dlt))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'tmp k',h_klt,i_klt))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_int,ns,'tmp t',h_tlt,i_tlt))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)


      call smd_bondlist_set(ns,nb,na,
     >                       int_mb(i_il1t),
     >                       int_mb(i_il2t),
     >                       int_mb(i_tlt),
     >                       dbl_mb(i_dlt),
     >                       dbl_mb(i_klt),
     >                       int_mb(i_ib1),
     >                       int_mb(i_ib2),
     >                       int_mb(i_itb),
     >                       dbl_mb(i_db),
     >                       dbl_mb(i_kb),
     >                       int_mb(i_it))
c
c     create bond list structure
c     ---------------------------
      if(ns.eq.0) then
        result = .false.
        goto 200
      end if
      call smd_namespace_create(sp_bondlist)
      tag = "bond:i1"
      call smd_data_create_get(sp_bondlist,tag,ns,MT_INT,i_il1)
      tag = "bond:i2"
      call smd_data_create_get(sp_bondlist,tag,ns,MT_INT,i_il2)
      tag = "bond:distance"
      call smd_data_create_get(sp_bondlist,tag,ns,MT_DBL,i_dl)
      tag = "bond:strength"
      call smd_data_create_get(sp_bondlist,tag,ns,MT_DBL,i_kl)
      tag = "bond:type"
      call smd_data_create_get(sp_bondlist,tag,ns,MT_INT,i_tl)

c
      do i=1,ns
       int_mb(i_il1+i-1) = int_mb(i_il1t+i-1)
       int_mb(i_il2+i-1) = int_mb(i_il2t+i-1)
       int_mb(i_tl+i-1)  = int_mb(i_tlt+i-1)
       dbl_mb(i_dl+i-1)  = dbl_mb(i_dlt+i-1)
       dbl_mb(i_kl+i-1)  = dbl_mb(i_klt+i-1)
      end do

200   continue

      if(.not.ma_pop_stack(h_tlt))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_klt))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_dlt))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_il2t))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_il1t))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      return
      end

      subroutine smd_bondlist_set(ns,nb,na,
     >                       il1,
     >                       il2,
     >                       itl,
     >                       dl,
     >                       kl,
     >                       ib1,
     >                       ib2,
     >                       itb,
     >                       db,
     >                       kb,
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
      integer il1(ns)
      integer il2(ns)
      integer itl(ns)
      double precision  dl(ns)
      double precision  kl(ns)
      integer ib1(nb)
      integer ib2(nb)
      integer itb(nb)
      double precision  db(ns)
      double precision  kb(ns)
      integer it(na)
c
      integer i,i1,i2,j,nlist
c
      nlist = 0
      do i=1,nb
        i1=0
        i2=0
        do j=1,na
         if(it(j).eq.ib1(i)) i1=j
         if(it(j).eq.ib2(i)) i2=j
         if(i1*i2.ne.0) then
          nlist = nlist + 1
          il1(nlist) = min(i1,i2)
          il2(nlist) = max(i1,i2)
          itl(nlist) = itb(i)
          kl(nlist)  = kb(i)
          dl(nlist)  = db(i)
          i1=0
          i2=0
         end if
        end do
      end do
      ns = nlist

      return
      end

c $Id$
