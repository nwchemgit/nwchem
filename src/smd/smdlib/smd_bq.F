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
      character*(smd_string_size) tag,pname
      character*(smd_string_size) sp_geom,sp_bq
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
      integer i_sht
c
      pname = "smd_geom_bq_init"
c
      tag = "geom"
      call smd_system_get_component(sp_geom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'need geom component '//tag,0,0)
c
      tag = "bq"
      call smd_system_get_component(sp_bq,tag,result)
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
      call smd_data_get_ind_dim(sp_geom,tag,i_im,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "size"
      call smd_data_get_ind_dim(sp_geom,tag,i_nm,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "neighbor_index"
      call smd_data_get_ind_dim(sp_geom,tag,i_nb,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )

      tag = "shift vector"
      call smd_data_get_index(sp_geom,tag,i_sht,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
      tag = "right index"
      call smd_data_create_get(sp_bq,tag,nb+1,MT_INT,i_rbq)
      tag = "left index"
      call smd_data_create_get(sp_bq,tag,nb+1,MT_INT,i_lbq)
c
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
      if(.not.bq_create("monomer charges test",bq_handle))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)

      if(.not.bq_alloc(bq_handle,nbq))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)

      if(.not.bq_index_coord(bq_handle,i_cbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      if(.not.bq_index_charge(bq_handle,i_qbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      scale = 1.88972598858d0
      j=0
      do ip=1,nb+1
      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
         j = j+1
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=dbl_mb(i_c+3*(i-1)+k-1)
     >                                -dbl_mb(i_sht+k-1)       
         write(102,*) dbl_mb(i_c+3*(i-1)+k-1),dbl_mb(i_sht+k-1),
     >                dbl_mb(i_cbq+3*(j-1)+k-1)
         end do
         dbl_mb(i_qbq+j-1)=dbl_mb(i_q+i-1)
      end do
      end do
      call smd_lat_rebox(nbq,dbl_mb(i_cbq))

      do i=1,3*nbq
        dbl_mb(i_cbq+i-1) = dbl_mb(i_cbq+i-1)*scale
      end do
      call bq_print_xyz(bq_handle,101)

      if(.not.bq_activate(bq_handle))
     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)
c

c        call bq_print_info(bq_handle)




      return
      end


