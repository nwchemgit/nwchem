      subroutine smd_bq_init()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "bq.fh"
#include "smd_const_data.fh"
c
      integer im
c     
      character*32 tag,pname,sp_geom,s_bq
      logical result
      integer i_q,i_c
      integer i,j,ip,k
      integer nt,nm,i0,nr
      integer rtdb
      integer bq_handle
      integer nbq,i_qbq,i_cbq
      double precision scale
      double precision cm(3)
      integer i_pr
      integer i_nb,h_nb,nb
      integer ir
      integer i_lbq,h_lbq,i_rbq,h_rbq
      logical nobq
      integer ng
      integer i_im,h_im
      integer i_nm,h_nm
c
      pname = "smd_geom_bq_init"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
c
      tag = "bq"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'no component '//tag,0,0)
      call smd_namespace_create(sp_bq)
c
      call smd_rtdb_get_handle(rtdb)
c     
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      nt=nt/3
c
      tag = "charge"
      call smd_get_ind(tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
      tag = "atom:res_pointer"
      call smd_get_ind_dim(tag,i_pr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c      
c     geometry data
c     -------------
      tag = "position"
      call smd_data_get_index_dim(sp_geom,tag,i_im,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "size"
      call smd_data_get_index_dim(sp_geom,tag,i_nm,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "neighbor_index"
      call smd_data_get_index_dim(sp_geom,tag,i_nb,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
c      get some memory for bq construction
c      -----------------------------------
       if(.not.ma_push_get(mt_int,nb+1,'tmp',h_rbq,i_rbq))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nb, MA_ERR)
       if(.not.ma_push_get(mt_int,nb+1,'tmp',h_lbq,i_lbq))
     + call errquit(pname//'Failed to allocate memory for tmp l',
     + nb, MA_ERR)
      
      int_mb(i_lbq) = 1
      int_mb(i_rbq+nb) = nt
      j=1
      do i=1,nb
        ir=int_mb(i_nb+i-1)
        write(44,*) "residue ir",ir,int_mb(i_pr+ir-1),
     +               int_mb(i_pr+ir)   
        int_mb(i_rbq+i-1) = int_mb(i_pr+ir-1)-1
        int_mb(i_lbq+i) = int_mb(i_pr+ir)
      end do
      call smd_util_sort(nb+1,int_mb(i_rbq))
      call smd_util_sort(nb+1,int_mb(i_lbq))
      write(44,*) "boundaries for bqs"
      do i=1,nb+1
        write(44,*) i,int_mb(i_lbq+i-1),
     +              int_mb(i_rbq+i-1)
      end do
c
c     now do the charges
c
      nbq=0
      do i=1,nb+1
        nbq=nbq+int_mb(i_rbq+i-1)-
     +              int_mb(i_lbq+i-1)+1
      end do
      write(*,*) "nbq=",nbq
c      nbq = nt - ng
c      if(.not.bq_create("monomer charges",bq_handle))
c     + call errquit(pname//'Failed bq_create',0,CALC_ERR)
c
c      if(.not.bq_alloc(bq_handle,nbq))
c     + call errquit(pname//'Failed bq_create',0,CALC_ERR)
c
c      if(.not.bq_index_coord(bq_handle,i_cbq))
c     + call errquit(pname//'Failed bq_index',0,CALC_ERR)
c
c      if(.not.bq_index_charge(bq_handle,i_qbq))
c     + call errquit(pname//'Failed bq_index',0,CALC_ERR)
c
c      scale = 1.88972598858d0
c
c      j=0
c      do ip=1,nb+1
c      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
c         j = j+1
c         do k=1,3
c            dbl_mb(i_cbq+3*(j-1)+k-1)=dbl_mb(i_ct+3*(i-1)+k-1)
c     >                                *scale
c            write(*,*) dbl_mb(i_cbq+3*(j-1)+k-1),
c     >                 dbl_mb(i_ct+3*(i-1)+k-1)
c         end do
c         dbl_mb(i_qbq+j-1)=dbl_mb(i_q+i-1)
c      end do
c      end do
c
c      call bq_print_xyz(bq_handle,98)
c
c      if(.not.bq_activate(bq_handle))
c     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)
cc
c
cc        call bq_print_info(bq_handle)
      if(.not.ma_pop_stack(h_lbq))
     & call errquit(pname//'Failed to deallocate stack h_lbq',nm,
     &       MA_ERR)
      if(.not.ma_pop_stack(h_rbq))
     & call errquit(pname//'Failed to deallocate stack h_rbq',nm,
     &       MA_ERR)

c


      return
      end


