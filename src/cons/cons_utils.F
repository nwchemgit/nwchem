c $Id: cons_utils.F,v 1.18 2004-11-29 22:16:58 marat Exp $
      function cons_print()
      implicit none
#include "global.fh"
c     
      logical cons_print
      cons_print=(ga_nodeid().eq.0)
      return
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

      character*255 tag_id
      character*255 tag_n
      character*255 namespace

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
     &      'cons_get_active_atoms: unable to allocate scratch space',
     &      n, MA_ERR)

       status = rtdb_get(rtdb,tag_id,
     >     mt_int,n,int_mb(i_tmp_id)) 
       if(.not.status) call errquit(
     >     'cons_get_active_atoms: unable to get prior fixed atoms',
     >      0, MA_ERR)
       
       do i=1,n
         id=int_mb(i_tmp_id +i-1)
         if(id .gt. nat) goto 1
         oactive(id)=.false.
         nactive = nactive - 1
       end do

1     continue
      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_get_active_atoms: unable to free h_tmp_id',
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
 
      subroutine cons_fxd_atoms_id_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
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
      character*(*) cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":fxd_atom_n"

      end
 
      subroutine cons_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)

      end

      subroutine cons_status_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":status"

      end

      subroutine cons_hbond_id_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
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
      character*255 cons_name
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
      character*255 cons_name
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
      character*255 cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hbond_k"

      end
c


