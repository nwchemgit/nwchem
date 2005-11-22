c $Id: cons_rtdb.F,v 1.4 2005-11-22 20:28:02 marat Exp $
      subroutine cons_create(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      
      call cons_set_namespace(namespace)
      call cons_set_irtdb(rtdb)

c     return if space already exists
      if(cons_inquire(namespace,rtdb)) return

      call cons_tag(namespace,tag)

c     contraints will be active by default
      if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit('cons_create', 0,
     >        RTDB_ERR)
      
c     no fixed atoms by default
      call cons_fxd_atoms_n_tag(namespace,tag)
      if(.not.rtdb_put(rtdb,tag,mt_int,1,0)) 
     >      call errquit('cons_create', 0,
     >        rtdb_err)

c     no springs by default
      call cons_hbond_n_tag(namespace,tag)
      if(.not.rtdb_put(rtdb,tag,mt_int,1,0)) 
     >      call errquit('cons_create', 0,
     >        rtdb_err)


      end

      subroutine cons_activate(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      
      if(.not.cons_inquire(namespace,rtdb)) then
           call errquit('nonexistent constraints', 0,
     >        0)

      end if

      if(.not. rtdb_cput(rtdb,"constraints" , 1, namespace)) then
         call errquit('qmmm: Failed constraints',0,0)
      end if
      

      end

      subroutine cons_get_namespace(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons.fh"
#include "cons_params.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      
      if(.not. rtdb_cget(rtdb,"constraints" , 1, namespace)) then
         namespace=cons_default
      end if
      

      end

      function cons_inquire(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_inquire
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag
      integer n

      call cons_fxd_atoms_n_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_int,1,n)
      if(cons_inquire) return

      call cons_hbond_n_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_int,1,n)
      if(cons_inquire) return
      
      return
      end

      subroutine cons_add_fixed_atoms(namespace,nat,id,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*(*) namespace
      integer nat
      integer id(nat)
      integer rtdb
c
      logical status
      integer i,j
      integer n0
      integer n
      integer h_tmp_id,i_tmp_id 
      
      character*255 tag_id
      character*255 tag_n


      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)

      call cons_sort_i(nat,id)
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n0)) then
         n0 = 0
      end if
      n = n0 + nat

      if ( .not. ma_push_get( MT_INT, n, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      n, MA_ERR)

       if(n0 .gt. 0) then
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
       do i=1,n0
        do j=1,nat
           if(int_mb(i_tmp_id + i-1).eq.id(j)) then
             id(j)=-1
           end if
        end do
       end do
       
       j=0
       do i=1,nat
         if(id(i) .ne. -1) then
          int_mb(i_tmp_id + n0+j)=id(i)
          j=j+1
         end if
       end do
       n=n0+j

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

      subroutine cons_add_spring(namespace,id1,id2,k,r,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*(*) namespace
      integer id1,id2
      double precision k
      double precision r
      integer rtdb
c
      logical status
      integer i
      integer nhb0
      integer nhb
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      integer id01,id02
      double precision k0
      double precision r0
      character*255 tag_id
      character*255 tag_n
      character*255 tag_r
      character*255 tag_k

      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb0)) then
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
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,2*nhb0,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,tag_k,
     >     mt_dbl,nhb0,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,tag_r,
     >     mt_dbl,nhb0,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)
       status =              rtdb_delete(rtdb,tag_n)
       status = status .and. rtdb_delete(rtdb,tag_id)
       status = status .and. rtdb_delete(rtdb,tag_k)
       status = status .and. rtdb_delete(rtdb,tag_r)
       if(.not.status) call errquit(
     >     'cons_add_spring: unable to delete prior harm bonds',
     >      0, MA_ERR)
       end if
       
      status = .true.
      do i=1,nhb0
       id01 = int_mb(i_tmp_id+2*(i-1))
       id02 = int_mb(i_tmp_id+2*(i-1)+1)
       k0 = dbl_mb(i_tmp_k+i-1)
       r0 = dbl_mb(i_tmp_r+i-1)
       if(id01.eq.id1 .and. id02.eq.id2) then
         write(*,*) "cons_add_spring: replacing existing spring"
         dbl_mb(i_tmp_k+i-1)=k
         dbl_mb(i_tmp_r+i-1)=r
         status = .false.
       end if
      end do
        
      if(status) then
       int_mb(i_tmp_id + 2*nhb0)     = id1
       int_mb(i_tmp_id + 2*nhb0 + 1) = id2
       dbl_mb(i_tmp_k + nhb0)=k
       dbl_mb(i_tmp_r + nhb0)=r
      else
       nhb=nhb0
      end if



       status = rtdb_put(rtdb,tag_n,mt_int,1,nhb)
       status = status . and. rtdb_put(rtdb,tag_id,
     >     mt_int,2*nhb,int_mb(i_tmp_id)) 
        status = status . and. rtdb_put(rtdb,tag_k,
     >     mt_dbl,nhb,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_put(rtdb,tag_r,
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

      subroutine cons_unfix_atoms(namespace,nat,id,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*(*) namespace
      integer nat
      integer id(nat)
      integer rtdb
c
      logical status
      integer i,j
      integer n0
      integer n
      integer h_tmp_id,i_tmp_id 

      character*255 tag_id
      character*255 tag_n

      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n0)) then
         n=0
      end if
      if(n.eq.0) then
         return
      end if
      n = n0

      if ( .not. ma_push_get( MT_INT, n, 'h_tmp_id',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_add_fixed_atoms: unable to allocate scratch space',
     &      n, MA_ERR)

        status = rtdb_get(rtdb,tag_id,
     >     mt_int,n0,int_mb(i_tmp_id)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
        status =  rtdb_delete(rtdb,tag_id)
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to delete prior fixed atoms',
     >      0, MA_ERR)

       do i=1,n0
        do j=1,nat
           if(int_mb(i_tmp_id + i-1).eq.id(j)) then
             int_mb(i_tmp_id + i-1)=-1
           end if
        end do
       end do

       call cons_sort_i(n0,int_mb(i_tmp_id))
       
       j=0
       do i=1,n0
         if(int_mb(i_tmp_id + i-1).eq.-1) then
           j=j+1
         end if
       end do


       status = rtdb_put(rtdb,tag_n,mt_int,1,n-j)
       status = status . and. rtdb_put(rtdb,tag_id,
     >     mt_int,n-j,int_mb(i_tmp_id+j)) 
        if(.not.status) call errquit(
     >     'cons_add_fixed_atoms: unable to get prior harm bonds',
     >      n, MA_ERR)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_info_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical status
      integer i
      integer n
      integer h_tmp_id,i_tmp_id 

      character*255 tag_id
      character*255 tag_n

      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n)) then
         n = 0
      end if
      if (n .eq. 0 ) then
         if(ga_nodeid().eq.0)
     >     write(*,*) "no fixed atoms found"
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
       
       if(ga_nodeid().eq.0) then
          write(*,'(" fixed atoms",T31,":",(T32,4(I4,2X),:))')
     >         (int_mb(i_tmp_id +i-1),i=1,n)
          write(*,'(/)')
       end if

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_info_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

      end

      subroutine cons_info_springs(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c
      character*(*) namespace
      integer rtdb
c     
      logical status
      integer i
      integer nhb
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      character*255 tag_id
      character*255 tag_n
      character*255 tag_k
      character*255 tag_r

      call cons_hbond_id_tag(namespace,tag_id)
      call cons_hbond_n_tag(namespace,tag_n)
      call cons_hbond_k_tag(namespace,tag_k)
      call cons_hbond_r_tag(namespace,tag_r)

      
      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) then
         nhb = 0
      end if
       
      if (nhb .eq. 0 ) then
         if(ga_nodeid().eq.0)    
     >     write(*,*) "no springs found"
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
        status = rtdb_get(rtdb,tag_id,
     >     mt_int,2*nhb,int_mb(i_tmp_id)) 
        status = status . and. rtdb_get(rtdb,tag_k,
     >     mt_dbl,nhb,dbl_mb(i_tmp_k)) 
        status = status . and. rtdb_get(rtdb,tag_r,
     >     mt_dbl,nhb,dbl_mb(i_tmp_r)) 
        if(.not.status) call errquit(
     >     'cons_add_spring: unable to get prior harm bonds',
     >      nhb, MA_ERR)
       end if

     
      if(ga_nodeid().eq.0) then   
      do i=1,nhb
      write(*,FMT=200) int_mb(i_tmp_id+2*(i-1)),
     >           int_mb(i_tmp_id+2*(i-1)+1),
     >           dbl_mb(i_tmp_k+i-1),
     >           dbl_mb(i_tmp_r+i-1)
      end do
      end if
200   FORMAT(" spring bond (i,j,k,r0)",T29,I6,2X,I6,2X,F6.3,2X,F6.3)
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

      subroutine cons_print_gen_info(namespace,rtdb)
      implicit none
#include "util.fh"
#include "global.fh"
      character*(*) namespace
      integer rtdb

      if(ga_nodeid().eq.0) then
        call cons_print_header()
        call util_print_centered(6,
     >     "General Information",
     >     1,.true.)
      end if

      call cons_info_fixed_atoms(namespace,rtdb)
      call cons_info_springs(namespace,rtdb)

      end

      function cons_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_fixed_atoms
      character*(*) namespace
      integer rtdb
c
      integer nfa
      character*255 tag_n

      call cons_fxd_atoms_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nfa)) nfa=0
     
      if(nfa.eq.0) then
        cons_fixed_atoms=.false.
      else
        cons_fixed_atoms=.true.
      end if

      return
      end

      function cons_hbonds(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_hbonds
      character*(*) namespace
      integer rtdb
c
      integer nhb
      character*255 tag_n

      call cons_hbond_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) nhb=0
     
      if(nhb.eq.0) then
        cons_hbonds=.false.
      else
        cons_hbonds=.true.
      end if

      return
      end

      function cons_status(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_status
      character*(*) namespace
      integer rtdb
c
      character*255 tag

      call cons_status_tag(namespace,tag)

      if(.not.rtdb_get(rtdb,tag,mt_log,1,cons_status)) then
        cons_status=.true.
      end if
     
      return
      end

      subroutine cons_set_status(namespace,cons_status,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_status
      character*(*) namespace
      integer rtdb
c
      character*255 tag

      call cons_status_tag(namespace,tag)

      if(.not.rtdb_put(rtdb,tag,mt_log,1,cons_status)) then
        call errquit("cons_set_status",0,RTDB_ERR)
      end if
     
      return
      end

      subroutine cons_delete_fixed_atoms(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
c     
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character*255 tag_id
      character*255 tag_n
c
      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)
c
      ignore = rtdb_delete(rtdb,tag_n)
      ignore = rtdb_delete(rtdb,tag_id)
c
      end

