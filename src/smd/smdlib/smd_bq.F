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
#include "stdio.fh"
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
      integer nb
      integer i_lbq,h_lbq,i_rbq,h_rbq
      logical nobq
      integer ng
      integer i_im,h_im
      integer i_nm,h_nm
      integer i_sht
c
      pname = "smd_geom_bq_init"
c
      call smd_rtdb_get_handle(rtdb)
c
      if(.not.rtdb_get(rtdb,'smd:nobq',mt_log,1,nobq))
     +     nobq = .false.
      if(nobq) return
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
c     set the boundary regions for bq 
c     -------------------------------
      int_mb(i_lbq) = 1
      int_mb(i_rbq+nb) = nt
      do i=1,nb
        int_mb(i_rbq+i-1) = int_mb(i_im+i-1)-1
        int_mb(i_lbq+i) = int_mb(i_im+i-1)+
     >                    int_mb(i_nm+i-1)   
      end do
      call smd_util_sort(nb+1,int_mb(i_rbq))
      call smd_util_sort(nb+1,int_mb(i_lbq))
c      write(44,*) "boundaries for bqs"
c      do i=1,nb+1
c        write(44,*) i,int_mb(i_lbq+i-1),
c     +              int_mb(i_rbq+i-1)
c      end do
c
c     now do the charges
c
      nbq=0
      do i=1,nb+1
        nbq=nbq+int_mb(i_rbq+i-1)-
     +              int_mb(i_lbq+i-1)+1
      end do
      if(nbq.eq.0) return
      if(ga_nodeid().eq.0) 
     +  write(luout,*) "Total number of Bq's",nbq
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
         end do
         dbl_mb(i_qbq+j-1)=dbl_mb(i_q+i-1)
      end do
      end do
      call smd_lat_rebox(nbq,dbl_mb(i_cbq))

      do i=1,3*nbq
        dbl_mb(i_cbq+i-1) = dbl_mb(i_cbq+i-1)*scale
      end do
c      call bq_print_xyz(bq_handle,101)

      if(.not.bq_activate(bq_handle))
     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)
c

c        call bq_print_info(bq_handle)




      return
      end

      subroutine smd_bq_dipole_init()
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
#include "stdio.fh"
c
      integer im
c     
      character*(smd_string_size) tag,pname
      character*(smd_string_size) sp_geom,sp_bq,sp_atom
      logical result
      integer i_q,i_c
      integer i,j,ip,k
      integer nt,nm,i0
      integer i_pr,i_mr,nr
      integer rtdb
      integer bq_handle
      integer nbq,i_qbq,i_cbq
      double precision scale
      double precision cm(3)
      integer i_nb,nb
      integer i_lbq,h_lbq,i_rbq,h_rbq
      logical nobq
      integer ng
      integer i_im,h_im
      integer i_nm,h_nm
      integer i_sht
      integer h_fd,i_fd
      double precision coc(3)
      double precision shift(3)
      double precision rd
      double precision qtmp
c
      pname = "smd_bq_dipole_init"
c
      rd=1.0d-2
c
      call smd_rtdb_get_handle(rtdb)
c
      if(.not.rtdb_get(rtdb,'smd:nobq',mt_log,1,nobq))
     +     nobq = .false.
      if(nobq) return
c      
c      write(luout,*) "reading dipoles"
c      call smd_atom_dipole_read(.true.)
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'need component '//tag,0,0)
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
c     
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      nt=nt/3
c
      tag = "atom:formal_charge"
      call smd_get_ind(tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
c     residue data
c     ------------
      tag = "atom:res_pointer"
      call smd_get_ind_dim(tag,i_pr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      write(*,*) "i_pr",(int_mb(i_pr+i-1),i=1,nr)
      tag = "atom:res_size"
      call smd_get_ind_dim(tag,i_mr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c
c     dipole   data
c     -------------
      if(.not.ma_push_get(mt_dbl,3*nr,'tmp i_fd',h_fd,i_fd))
     + call errquit(pname//'Failed to allocate memory for h_fd',0,0)
      tag = "res_dipoles"
      call smd_data_get_dbl1(sp_atom,tag,1,3*nr,dbl_mb(i_fd))
c
c     geometry data
c     -------------
      tag = "neighbor_index"
      call smd_data_get_ind_dim(sp_geom,tag,i_nb,nb,result)

      tag = "position"
      call smd_data_get_ind_dim(sp_geom,tag,i_im,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
      
      tag = "size"
      call smd_data_get_ind_dim(sp_geom,tag,i_nm,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "shift vector"
      call smd_data_get_index(sp_geom,tag,i_sht,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
      call smd_data_get_dbl(sp_geom,tag,3,shift,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
      tag = "right index"
      call smd_data_create_get(sp_bq,tag,nb+1,MT_INT,i_rbq)
      tag = "left index"
      call smd_data_create_get(sp_bq,tag,nb+1,MT_INT,i_lbq)
c
c     set the boundary regions for bq 
c     -------------------------------
      int_mb(i_lbq) = 1
      int_mb(i_rbq+nb) = nr
      do i=1,nb
        int_mb(i_rbq+i-1) = int_mb(i_nb+i-1)-1
        int_mb(i_lbq+i) = int_mb(i_nb+i-1)+1
      end do
      call smd_util_sort(nb+1,int_mb(i_rbq))
      call smd_util_sort(nb+1,int_mb(i_lbq))
      write(6,*) "boundaries for bqs"
      do i=1,nb+1
        write(6,*) i,int_mb(i_lbq+i-1),
     +              int_mb(i_rbq+i-1)
      end do
      if(ga_nodeid().eq.0) then
      write(*,*) "printing dipoles"
      do i=1,nr
         write(*,*) (dbl_mb(i_fd+3*(i-1)+k-1),k=1,3) 
      end do
      end if
c
c
      nbq=0
      do i=1,nb+1
        nbq=nbq+int_mb(i_rbq+i-1)-
     +              int_mb(i_lbq+i-1)+1
      end do
      nbq = nbq*2
      if(nbq.eq.0) return
      if(ga_nodeid().eq.0) 
     +  write(luout,*) "Total number of Bq's",nbq
      if(.not.bq_create("monomer charges test",bq_handle))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)

      if(.not.bq_alloc(bq_handle,nbq))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)

      if(.not.bq_index_coord(bq_handle,i_cbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      if(.not.bq_index_charge(bq_handle,i_qbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      if(ga_nodeid().eq.0) 
     +   call smd_data_namespace_print(sp_geom,203)
      scale = 1.88972598858d0
      j=0
      do ip=1,nb+1
      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
         call smd_geom_coc(nt,int_mb(i_mr+i-1),int_mb(i_pr+i-1),
     +                  dbl_mb(i_q),dbl_mb(i_c),coc)
         qtmp=0.0d0 
         do k=1,3
            write(*,*) "i_fd",i,k,dbl_mb(i_fd+3*(i-1)+k-1)
            qtmp=qtmp+dbl_mb(i_fd+3*(i-1)+k-1)**2
         end do
         qtmp = SQRT(qtmp)
         if(qtmp.ne.0) then
         do k=1,3
            dbl_mb(i_fd+3*(i-1)+k-1)=dbl_mb(i_fd+3*(i-1)+k-1)/qtmp
         end do
         end if
         j = j+1
         dbl_mb(i_qbq+j-1)=qtmp/rd
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=coc(k)
     >                                -dbl_mb(i_sht+k-1)       
         end do
         j=j+1
         dbl_mb(i_qbq+j-1)=-qtmp/rd
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=coc(k)
     >                                -dbl_mb(i_sht+k-1)       
         end do
      end do
      end do
      call smd_lat_rebox(nbq,dbl_mb(i_cbq))
      j=0
      do ip=1,nb+1
      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
         j = j+1
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=dbl_mb(i_cbq+3*(j-1)+k-1)+
     +                             0.5d0*rd*dbl_mb(i_fd+3*(i-1)+k-1)
         end do
         j=j+1
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=dbl_mb(i_cbq+3*(j-1)+k-1)-
     +                             0.5d0*rd*dbl_mb(i_fd+3*(i-1)+k-1)
         end do
      end do
      end do

      j=0
      do ip=1,nb+1
      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
         do k=1,3
           coc(k) = 0.0d0
         end do
         j = j+1
         do k=1,3
            coc(k) = coc(k) + dbl_mb(i_cbq+3*(j-1)+k-1)*
     +                        dbl_mb(i_qbq+j-1)
         end do
         j=j+1
         do k=1,3
            coc(k) = coc(k) + dbl_mb(i_cbq+3*(j-1)+k-1)*
     +                        dbl_mb(i_qbq+j-1)
         end do
         if(ga_nodeid().eq.0) then
         write(*,*) "dipole comparision"
         do k=1,3
           write(*,*) coc(k),dbl_mb(i_fd+3*(i-1)+k-1)
         end do
         end if
      end do
      end do

      do i=1,3*nbq
        dbl_mb(i_cbq+i-1) = dbl_mb(i_cbq+i-1)*scale
      end do
      call bq_print_xyz(bq_handle,201)

      if(.not.bq_activate(bq_handle))
     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)
c

c        call bq_print_info(bq_handle)


      if(.not.ma_chop_stack(h_fd))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)


      return
      end

      subroutine smd_bq_dipole_init1()
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
#include "stdio.fh"
c
      integer im
c     
      character*(smd_string_size) tag,pname
      character*(smd_string_size) sp_geom,sp_bq,sp_atom
      logical result
      integer i_q,i_c
      integer i,j,ip,k
      integer nt,nm,i0
      integer i_pr,i_mr,nr
      integer rtdb
      integer bq_handle
      integer nbq,i_qbq,i_cbq
      double precision scale
      double precision cm(3)
      integer i_nb,nb
      integer i_lbq,h_lbq,i_rbq,h_rbq
      logical nobq
      integer ng
      integer i_im,h_im
      integer i_nm,h_nm
      integer i_sht
      integer h_fd,i_fd
      double precision coc(3)
      double precision shift(3)
      double precision rd
      double precision qtmp
c
      pname = "smd_bq_dipole_init"
c
      rd=1.0d-2
c
      call smd_rtdb_get_handle(rtdb)
c
      if(.not.rtdb_get(rtdb,'smd:nobq',mt_log,1,nobq))
     +     nobq = .false.
      if(nobq) return
c      
c      write(luout,*) "reading dipoles"
c      call smd_atom_dipole_read(.true.)
c
      tag = "atom"
      call smd_system_get_component(sp_atom,tag,result)
      if(.not.result)
     >  call errquit(
     >       pname//'need component '//tag,0,0)
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
c     
      tag = "coords"
      call smd_get_ind_dim(tag,i_c,nt,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
      nt=nt/3
c
      tag = "atom:formal_charge"
      call smd_get_ind(tag,i_q,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for '//tag,0, RTDB_ERR)
c
c     residue data
c     ------------
      tag = "atom:res_pointer"
      call smd_get_ind_dim(tag,i_pr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
      write(*,*) "i_pr",(int_mb(i_pr+i-1),i=1,nr)
      tag = "atom:res_size"
      call smd_get_ind_dim(tag,i_mr,nr,result)
      if(.not. result) 
     >  call errquit(
     >       pname//'error getting index for'//tag,0, RTDB_ERR)
c
c
c     dipole   data
c     -------------
      if(.not.ma_push_get(mt_dbl,3*nr,'tmp i_fd',h_fd,i_fd))
     + call errquit(pname//'Failed to allocate memory for h_fd',0,0)
      tag = "res_dipoles"
      call smd_data_get_dbl1(sp_atom,tag,1,3*nr,i_fd)
c
c     geometry data
c     -------------
      tag = "neighbor_index"
      call smd_data_get_ind_dim(sp_geom,tag,i_nb,nb,result)

      tag = "position"
      call smd_data_get_ind_dim(sp_geom,tag,i_im,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
      
      tag = "size"
      call smd_data_get_ind_dim(sp_geom,tag,i_nm,nb,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c      
      tag = "shift vector"
      call smd_data_get_index(sp_geom,tag,i_sht,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
      call smd_data_get_dbl(sp_geom,tag,3,shift,result)
      if(.not. result)  call errquit(
     >       pname//'error getting index for'//tag,0,0 )
c
      tag = "right index"
      call smd_data_create_get(sp_bq,tag,nb+1,MT_INT,i_rbq)
      tag = "left index"
      call smd_data_create_get(sp_bq,tag,nb+1,MT_INT,i_lbq)
c
c     set the boundary regions for bq 
c     -------------------------------
      int_mb(i_lbq) = 1
      int_mb(i_rbq+nb) = nr
      do i=1,nb
        int_mb(i_rbq+i-1) = int_mb(i_nb+i-1)-1
        int_mb(i_lbq+i) = int_mb(i_nb+i-1)+1
      end do
      call smd_util_sort(nb+1,int_mb(i_rbq))
      call smd_util_sort(nb+1,int_mb(i_lbq))
      write(6,*) "boundaries for bqs"
      do i=1,nb+1
        write(6,*) i,int_mb(i_lbq+i-1),
     +              int_mb(i_rbq+i-1)
      end do
      if(ga_nodeid().eq.0) then
      write(*,*) "printing dipoles"
      do i=1,nr
         write(*,*) (dbl_mb(i_fd+3*(i-1)+k-1),k=1,3) 
      end do
      end if
c
c
      nbq=0
      do i=1,nb+1
        nbq=nbq+int_mb(i_rbq+i-1)-
     +              int_mb(i_lbq+i-1)+1
      end do
      nbq = nbq*2
      if(nbq.eq.0) return
      if(ga_nodeid().eq.0) 
     +  write(luout,*) "Total number of Bq's",nbq
      if(.not.bq_create("monomer charges test",bq_handle))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)

      if(.not.bq_alloc(bq_handle,nbq))
     + call errquit(pname//'Failed bq_create',0,CALC_ERR)

      if(.not.bq_index_coord(bq_handle,i_cbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      if(.not.bq_index_charge(bq_handle,i_qbq))
     + call errquit(pname//'Failed bq_index',0,CALC_ERR)

      if(ga_nodeid().eq.0) 
     +   call smd_data_namespace_print(sp_geom,203)
      scale = 1.88972598858d0
      j=0
      do ip=1,nb+1
      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
         call smd_geom_coc(nt,int_mb(i_mr+i-1),int_mb(i_pr+i-1),
     +                  dbl_mb(i_q),dbl_mb(i_c),coc)
         qtmp=0.0d0 
         do k=1,3
            write(*,*) "i_fd",i,k,dbl_mb(i_fd+3*(i-1)+k-1)
            qtmp=qtmp+dbl_mb(i_fd+3*(i-1)+k-1)**2
         end do
         qtmp = SQRT(qtmp)/rd
         write(*,*) "qtmp=",qtmp
         j = j+1
         dbl_mb(i_qbq+j-1)=qtmp
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=coc(k)
     >                                -dbl_mb(i_sht+k-1)       
         end do
         j=j+1
         dbl_mb(i_qbq+j-1)=-qtmp
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=coc(k)
     >                                -dbl_mb(i_sht+k-1)       
         end do
      end do
      end do
      call smd_lat_rebox(nbq,dbl_mb(i_cbq))
      j=0
      do ip=1,nb+1
      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
         j = j+1
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=dbl_mb(i_cbq+3*(j-1)+k-1)+
     +                             0.5d0*dbl_mb(i_fd+3*(i-1)+k-1)/
     +                                dbl_mb(i_qbq+j-1)
         end do
         j=j+1
         do k=1,3
            dbl_mb(i_cbq+3*(j-1)+k-1)=dbl_mb(i_cbq+3*(j-1)+k-1)+
     +                             0.5d0*dbl_mb(i_fd+3*(i-1)+k-1)/
     +                                dbl_mb(i_qbq+j-1)
         end do
      end do
      end do

      j=0
      do ip=1,nb+1
      do i=int_mb(i_lbq+ip-1),int_mb(i_rbq+ip-1)
         do k=1,3
           coc(k) = 0.0d0
         end do
         j = j+1
         do k=1,3
            coc(k) = coc(k) + dbl_mb(i_cbq+3*(j-1)+k-1)*
     +                        dbl_mb(i_qbq+j-1)
         end do
         j=j+1
         do k=1,3
            coc(k) = coc(k) + dbl_mb(i_cbq+3*(j-1)+k-1)*
     +                        dbl_mb(i_qbq+j-1)
         end do
         if(ga_nodeid().eq.0) then
         write(*,*) "dipole comparision"
         do k=1,3
           write(*,*) coc(k),dbl_mb(i_fd+3*(i-1)+k-1)
         end do
         end if
      end do
      end do

      do i=1,3*nbq
        dbl_mb(i_cbq+i-1) = dbl_mb(i_cbq+i-1)*scale
      end do
      call bq_print_xyz(bq_handle,201)

      if(.not.bq_activate(bq_handle))
     + call errquit(pname//'Failed bq_activate',0,CALC_ERR)
c

c        call bq_print_info(bq_handle)


      if(.not.ma_chop_stack(h_fd))
     & call errquit(pname//'Failed to deallocate stack h_l',0,
     &       MA_ERR)


      return
      end

