c $Id: cons_utils.F,v 1.12 2004-04-20 18:07:21 marat Exp $
      function cons_print()
      implicit none
#include "global.fh"
c     
      logical cons_print
      cons_print=(ga_nodeid().eq.0)
      return
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
      character(255) tag

      call cons_tag(namespace,tag)
      cons_inquire=rtdb_get(rtdb,tag,mt_log,1,ignore)
      return
      end

      function cons_activate(namespace,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c     
      logical cons_activate
      character*(*) namespace
      integer rtdb
c
      logical ignore
      character(255) tag
      
      call cons_tag(namespace,tag)
      cons_activate=rtdb_get(rtdb,tag,mt_log,1,ignore)

      if(cons_activate) call cons_set_namespace(namespace)

      return
      end

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
      character(255) tag
      
      call cons_tag(namespace,tag)

      if(cons_inquire(namespace,rtdb)) return

      if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit('cons_create', 0,
     >        RTDB_ERR)
      
      end

      subroutine cons_add_spring(namespace,i,j,k,r,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      character*(*) namespace
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
      integer i0,j0
      double precision k0
      double precision r0
      character(255) tag_id
      character(255) tag_n
      character(255) tag_r
      character(255) tag_k

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
       i0 = int_mb(i_tmp_id+2*(i-1))
       j0 = int_mb(i_tmp_id+2*(i-1)+1)
       k0 = dbl_mb(i_tmp_k+i-1)
       r0 = dbl_mb(i_tmp_r+i-1)
       if(i0.eq.i .and. j0.eq.j) then
         write(*,*) "cons_add_spring: replacing existing spring"
         dbl_mb(i_tmp_k+i-1)=k
         dbl_mb(i_tmp_r+i-1)=r
         status = .false.
       end if
      end do
        
      if(status) then
       int_mb(i_tmp_id + 2*nhb0)     = i
       int_mb(i_tmp_id + 2*nhb0 + 1) = j
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
      
      character(255) tag_id
      character(255) tag_n


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
       write(*,*) "cons_add_fixed_atoms",namespace
       write(*,*) "cons_add_fixed_atoms",(int_mb(i_tmp_id+i-1),i=1,n)

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_add_fixed_atoms: unable to free h_tmp_id',
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

      character(255) tag_id
      character(255) tag_n

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

      subroutine cons_get_active_atoms(nat,nactive,oactive,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
      integer nat
      logical oactive(nat)
      integer nactive
      integer rtdb
c
      logical status
      integer i
      integer n
      integer id
      integer h_tmp_id,i_tmp_id 

      character(255) tag_id
      character(255) tag_n
      character(255) namespace

      call cons_get_namespace(namespace)
      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)

      
      do i=1,nat
        oactive(i)=.true. 
      end do
      nactive = nat

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,n)) then
         n = 0
      end if
      if (n .eq. 0 ) then
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
       
       do i=1,n
         id=int_mb(i_tmp_id +i-1)
         if(id .gt. nat) goto 1
         oactive(id)=.false.
         nactive = nactive - 1
       end do

1     continue
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_info_fixed_atoms: unable to free h_tmp_id',
     &    0, MA_ERR)

       write(*,*) "cons_get_active_atoms",nactive,(oactive(i),i=1,nat)

      end

      subroutine cons_info_fixed_atoms(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
c
      logical status
      integer i
      integer n
      integer h_tmp_id,i_tmp_id 

      character(255) tag_id
      character(255) tag_n
      character(255) namespace

      call cons_get_namespace(namespace)
      call cons_fxd_atoms_id_tag(namespace,tag_id)
      call cons_fxd_atoms_n_tag(namespace,tag_n)
      
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
       
       if(ga_nodeid().eq.0) then
       write(*,*) "fixed atoms          ",(int_mb(i_tmp_id +i-1),i=1,n)
       end if

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
#include "global.fh"
c
      integer rtdb
c     
      logical status
      integer i
      integer nhb
      integer h_tmp_id,i_tmp_id 
      integer h_tmp_k,i_tmp_k 
      integer h_tmp_r,i_tmp_r 
      character(255) tag_id
      character(255) tag_n
      character(255) tag_k
      character(255) tag_r
      character(255) namespace

      call cons_get_namespace(namespace)
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

      subroutine cons_sort_i(n,a)
c     
      integer n
      integer a(n)
c
c     local variables:
      integer i
      integer pass  
      integer sorted 
      integer temp

      pass = 1
      sorted = 0
      do while(sorted .eq. 0) 
        sorted = 1
        do 2 i = 1,n-pass
          if(a(i) .gt. a(i+1)) then
            temp = a(i)
            a(i) = a(i+1)
            a(i+1) = temp
            sorted = 0
          endif
 2      continue
        pass = pass +1
      end do
      do i=1,n-1
       if(a(i).eq.a(i+1)) a(i)=-1 
      end do

      return
      end

      subroutine cons_print_header()
      implicit none
#include "util.fh"

      call util_print_centered(6,
     >     "NWChem Constraints Module",
     >     36,.true.)

      end

      subroutine cons_print_gen_info(rtdb)
      implicit none
#include "util.fh"
#include "global.fh"
      integer rtdb

      if(ga_nodeid().eq.0) then
        call cons_print_header()
        call util_print_centered(6,
     >     "General Information",
     >     1,.true.)
      end if

      call cons_info_fixed_atoms(rtdb)
      call cons_info_springs(rtdb)

      end
c
      subroutine cons_fxd_atoms_id_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character(255) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":fxd_atom_id"

      end

      subroutine cons_fxd_atoms_n_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character(255) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":fxd_atom_n"

      end
c
      subroutine cons_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character(255) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)

      end

      subroutine cons_hbond_id_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character(255) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hbond_id"

      end

      subroutine cons_hbond_n_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character(255) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hbond_n"

      end

      subroutine cons_hbond_r_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character(255) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hbond_r"

      end

      subroutine cons_hbond_k_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character(255) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hbond_k"

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
      character(255) tag_n

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
      character(255) tag_n

      call cons_hbond_n_tag(namespace,tag_n)

      if(.not.rtdb_get(rtdb,tag_n,mt_int,1,nhb)) nhb=0
     
      if(nhb.eq.0) then
        cons_hbonds=.false.
      else
        cons_hbonds=.true.
      end if

      return
      end

