c $Id: cons_utils.F,v 1.1 2004-03-29 20:34:16 marat Exp $
      subroutine cons_add_spring(i,j,k,r,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer i,j
      double precision k
      double precision r
      integer rtdb
c
      logical status
      integer nhb0
      integer nhb
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      
      if(.not.rtdb_get(rtdb,"cons:hbond_n",mt_int,1,nhb0)) then
         nhb0 = 0
      end if
      nhb = nhb0 + 1

      if ( .not. ma_push_get( MT_INT, 2*nhb, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_k',
     &      h_tmp_k, i_tmp_k) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_r',
     &      h_tmp_r, i_tmp_r) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)

       if(nhb0 .ge. 1) then
        status = rtdb_get(rtdb,"cons:hbond_ij",
     >     mt_int,2*nhb0,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,"cons:hbond_kij",
     >     mt_dbl,nhb0,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,"cons:hbond_rij",
     >     mt_dbl,nhb0,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)
       status =              rtdb_delete(rtdb,"cons:hbond_n")
       status = status .and. rtdb_delete(rtdb,"cons:hbond_ij")
       status = status .and. rtdb_delete(rtdb,"cons:hbond_kij")
       status = status .and. rtdb_delete(rtdb,"cons:hbond_rij")
       if(.not.status) call errquit(
     >     'cons_add_spring: unable to delete prior harm bonds',
     >      0, MA_ERR)
       end if
       
       int_mb(i_tmp_id + 2*nhb0)     = i
       int_mb(i_tmp_id + 2*nhb0 + 1) = j
       dbl_mb(i_tmp_k + nhb0)=k
       dbl_mb(i_tmp_r + nhb0)=r

       status = rtdb_put(rtdb,"cons:hbond_n",mt_int,1,nhb)
       status = status . and. rtdb_put(rtdb,"cons:hbond_ij",
     >     mt_int,2*nhb,int_mb(i_tmp_id)) 
        status = status . and. rtdb_put(rtdb,"cons:hbond_kij",
     >     mt_dbl,nhb,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_put(rtdb,"cons:hbond_rij",
     >     mt_dbl,nhb,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_r) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_r',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_k',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_spring: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_add_fixed_atoms(nat,id,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer nat
      integer id(nat)
      integer rtdb
c
      logical status
      integer i
      integer n0
      integer n
      integer h_tmp_id,i_tmp_id 
      
      character tag_id*(*) 
      parameter (tag_id = 'cons:fxd_at_id') 
      character tag_n*(*) 
      parameter (tag_n = 'cons:fxd_at_n') 

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n0)) then
         n0 = 0
      end if
      n = n0 + nat

      if ( .not. ma_push_get( MT_INT, n, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      n, MA_ERR)

       if(n0 .gt. 1) then
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,n0,int_mb(i_tmp_id)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
        status =  rtdb_delete(rtdb,tag_id)
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to delete prior fixed atoms',
     >      0, MA_ERR)
       end if
       
       do i=1,nat
         int_mb(i_tmp_id + n0+i-1)=id(i)
       end do

       status = rtdb_put(rtdb,tag_n,mt_int,1,n)
       status = status . and. rtdb_put(rtdb,tag_id,
     >     mt_int,n,int_mb(i_tmp_id)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior harm bonds',
     >      n, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_info_fixed_atoms(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer rtdb
c
      logical status
      integer i
      integer n
      integer h_tmp_id,i_tmp_id 
      
      character tag_id*(*) 
      parameter (tag_id = 'cons:fxd_at_id') 
      character tag_n*(*) 
      parameter (tag_n = 'cons:fxd_at_n') 

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n)) then
         n = 0
      end if
      if (n .eq. 0 ) then
         write(*,*) "no fixed atoms found"
         return
      end if

      if ( .not. ma_push_get( MT_INT, n, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_info_fixed_atoms: unable to allocate scratch space',
     &      n, MA_ERR)

       status = rtdb_get(rtdb,tag_id,
     >     mt_int,n,int_mb(i_tmp_id)) 
       if(.not.status) call errquit(
     >     'cons_info_fixed_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
       
       write(*,*) "fixed atoms",(int_mb(i_tmp_id +i-1),i=1,n)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_info_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_info_springs(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c
      integer rtdb
c     
      logical status
      integer i
      integer nhb
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      
      if(.not.rtdb_get(rtdb,"cons:hbond_n",mt_int,1,nhb)) then
         nhb = 0
      end if
       
      if (nhb .eq. 0 ) then
         write(*,*) "no springs found"
         return
      end if

      if ( .not. ma_push_get( MT_INT, 2*nhb, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      2*nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_k',
     &      h_tmp_k, i_tmp_k) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)
      if ( .not. ma_push_get( MT_DBL, nhb, 'h_tmp_r',
     &      h_tmp_r, i_tmp_r) ) call errquit(
     &      'cons_add_spring: unable to allocate scratch space',
     &      nhb, MA_ERR)

       if(nhb .ge. 1) then
        status = rtdb_get(rtdb,"cons:hbond_ij",
     >     mt_int,2*nhb,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,"cons:hbond_kij",
     >     mt_dbl,nhb,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,"cons:hbond_rij",
     >     mt_dbl,nhb,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)
       end if

      write(*,*) "sppprings"
      do i=1,nhb
      write(*,*) int_mb(i_tmp_id+2*(i-1)),
     >           int_mb(i_tmp_id+2*(i-1)+1),
     >           dbl_mb(i_tmp_k+i-1),
     >           dbl_mb(i_tmp_r+i-1)
      end do
      
      if (.not. ma_pop_stack(h_tmp_r) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_r',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_k) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_k',
     &    0, MA_ERR)
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_info_springs: unable to free h_tmp_id',
     &    0, MA_ERR)
      end

