c $Id$
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
      subroutine cons_hdihed_id_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hdihed_id"

      end

      subroutine cons_hdihed_n_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hdihed_n"

      end

      subroutine cons_hdihed_r_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hdihed_r"

      end

      subroutine cons_hdihed_k_tag(cons_name,tag)
      implicit none
#include "util.fh"
#include "inp.fh"

      character*(*) tag
c
      integer n
      character*255 cons_name
c
      n=inp_strlen(cons_name)
      tag="cons:"//cons_name(1:n)//":hdihed_k"

      end
c
      subroutine cons_dihed(r1,r2,r3,r4,phi)
      implicit none
#include "util.fh"
#include "inp.fh"

      double precision r1(3),r2(3),r3(3),r4(3)
      double precision phi
c
      double precision r12(3),r32(3),r34(3)
      double precision r52(3),r63(3)
      double precision s52,s63,s1263,s5263
      integer k
      double precision phi0
c    
      do k=1,3
       r12(k)=r1(k)-r2(k)
       r32(k)=r3(k)-r2(k)
       r34(k)=r3(k)-r4(k)
      end do
      
      call cons_vdot(r12,r32,r52)
      call cons_vdot(r32,r34,r63)
      call cons_ddot(r52,r52,s52)
      call cons_ddot(r63,r63,s63)
      call cons_ddot(r52,r32,s5263)

      phi0 = acos(s5263/sqrt(s52*s63))

      call cons_ddot(r12,r63,s1263)

      phi = sign(phi0,s1263)

      end
c      
      subroutine cons_dihed_deriv(r1,r2,r3,r4,d1,d2,d3,d4)
      implicit none
#include "util.fh"
#include "inp.fh"

      double precision r1(3),r2(3),r3(3),r4(3)
      double precision d1(3),d2(3),d3(3),d4(3)
c
      double precision r12(3),r32(3),r34(3)
      double precision r52(3),r63(3)
      double precision s32,s52,s63,s1232,s3432
      integer k
      double precision a1,a2,a3,a4
      double precision a31,a42,a31_1,a42_2
c    
      do k=1,3
       r12(k)=r1(k)-r2(k)
       r32(k)=r3(k)-r2(k)
       r34(k)=r3(k)-r4(k)
      end do
      
      call cons_vdot(r12,r32,r52)
      call cons_vdot(r32,r34,r63)

      call cons_ddot(r32,r32,s32)
      call cons_ddot(r52,r52,s52)
      call cons_ddot(r63,r63,s63)
      call cons_ddot(r12,r32,s1232)
      call cons_ddot(r34,r32,s3432)
     
      a1=sqrt(s32)/s52
      a2=sqrt(s32)/s63
      a3=s1232/s32
      a4=s3432/s32
      a31 = a3*a1
      a31_1 = a31-a1
      a42 = a4*a2
      a42_2 = a42-a2

      do k=1,3
        d1(k)= a1*r52(k)
        d2(k)= a31_1*r52(k)+a42*r63(k)
        d3(k)= a42_2*r63(k)+a31*r52(k)
        d4(k)=-a4*r63(k)
      end do

      end
c      
      subroutine cons_vdot(r1,r2,r12)
      implicit none
#include "util.fh"
#include "inp.fh"

      double precision r1(3),r2(3),r12(3)
c
      r12(1) = r1(2)*r2(3)-r1(3)*r2(2)
      r12(2) = r1(3)*r2(1)-r1(1)*r2(3)
      r12(3) = r1(1)*r2(2)-r1(1)*r2(2)
      end
c
      subroutine cons_ddot(r1,r2,s)
      implicit none
#include "util.fh"
#include "inp.fh"

      double precision r1(3),r2(3),s
c
      integer k
      s = 0.d0
      do k=1,3
        s = s+r1(k)*r2(k)
      end do
      s = sqrt(s)
      end
