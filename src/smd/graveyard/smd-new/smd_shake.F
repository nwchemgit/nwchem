      subroutine smd_shake_init(sp_shake,sp_param,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_shake
      character*(*) sp_param
      integer rtdb
c
      character*32 pname
      character*80 tag
      character*255 filename
      integer na,nt
      integer i_i1,i_i2,i_d,i_tt
      integer i_i1t,h_i1t
      integer i_i2t,h_i2t
      integer i_dt,h_dt
      logical result
      integer i
c
      pname = "smd_shake_init"
c
      write(*,*) "in "//pname
c
c     get array of type names
c     -----------------------
      tag = "param:typename"
      call smd_data_get_index(sp_param,tag,i_tt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      call smd_data_get_size(sp_param,tag,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting size for'//tag,0, RTDB_ERR)
      nt = nt/16
c
       na = nt*(nt-1)
       na = na/2
c
c     allocate initial storage for shake params
c     ----------------------------------------
      if(.not.ma_push_get(mt_int,na,'tmp i1',h_i1t,i_i1t))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_int,na,'tmp i2',h_i2t,i_i2t))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      if(.not.ma_push_get(mt_dbl,na,'tmp d',h_dt,i_dt))
     + call errquit(pname//'Failed to allocate memory',
     + 0, MA_ERR)

      call smd_paramfile_read_shake(na,nt,
     +                         int_mb(i_i1t),
     +                         int_mb(i_i2t),
     +                         dbl_mb(i_dt),
     +                         byte_mb(i_tt))

c
c     create parameter data structures
c     ---------------------------
      call smd_data_create(sp_shake,"shake:i1",na,MT_INT)
      call smd_data_create(sp_shake,"shake:i2",na,MT_INT)
      call smd_data_create(sp_shake,"shake:distance",na,MT_DBL)

c
c     get shake arrays
c     ----------------
      tag = "shake:i1"
      call smd_data_get_index(sp_shake,tag,i_i1,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "shake:i2"
      call smd_data_get_index(sp_shake,tag,i_i2,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      tag = "shake:distance"
      call smd_data_get_index(sp_shake,tag,i_d,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      do i=1,na
       int_mb(i_i1+i-1) = int_mb(i_i1t+i-1)
       int_mb(i_i2+i-1) = int_mb(i_i2t+i-1)
       dbl_mb(i_d+i-1)  = dbl_mb(i_dt+i-1)
      end do

      if(.not.ma_pop_stack(h_dt))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_i2t))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_i1t))
     & call errquit(pname//'Failed to deallocate stack',0,
     &       MA_ERR)

      return
      end

      subroutine smd_shake_ntot(sp_shake,na)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) sp_shake
      integer na
c
      logical result
      character*30 tag,pname
c
      pname = "smd_shake_ntot"
c
      tag = "shake:i1"
      call smd_data_get_index(sp_shake,tag,na,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, 0)

      return
      end

c $Id$
