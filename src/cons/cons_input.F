c $Id: cons_input.F,v 1.14 2004-09-21 19:12:09 marat Exp $
      subroutine cons_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_params.fh"
#include "cons.fh"
c     
      integer rtdb
c
      character*255 test
      character*255 namespace
      integer max_fp
      integer i_tmp_id,h_tmp_id
      integer nfp
      logical status
      integer id(2)
      double precision k,r
c
      logical cons_harm_bond_input
      external cons_harm_bond_input

      max_fp = min( max_fixed_atoms*4, ma_inquire_avail(MT_INT))
      if ( .not. ma_push_get( MT_INT, max_fp, 
     &      'cons scratch fixed atoms',
     &      h_tmp_id, i_tmp_id) ) call errquit(
     &      'cons_input: failed cons scratch space for fixed atoms',
     &      max_fp, MA_ERR)

c     set namespace for constraints if any
c     ------------------------------------
      if(inp_n_field().eq.1) then
       namespace=cons_default
      else if(inp_n_field().eq.2) then
       call inp_set_field(1)
       if (.not. inp_a(test)) then
         call errquit('cons_input: inp_read failed',0,0)
       else
         namespace=test
       end if
      end if 
      
 100  if (.not. inp_read()) 
     >      call errquit('cons_input: inp_read failed',0,0)
c
      if (.not. inp_a(test))
     $     call errquit('cons_input: failed to read keyword', 0,
     &       0)
c     
c     print level 
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, "cons")
c
      goto 100
c
c     reset/clear
c     ----------------------------
      else if (inp_compare(.false.,'clear', test)) then
         call cons_delete_hbonds(namespace,rtdb)
         call cons_delete_fixed_atoms(namespace,rtdb)

      goto 100

c      enable
c     ----------------------------
      else if (inp_compare(.false.,'enable', test)) then
         call cons_set_status(namespace,.true.,rtdb)

      goto 100
c
c     disable 
c     ----------------------------
      else if (inp_compare(.false.,'disable', test)) then
         call cons_set_status(namespace,.false.,rtdb)

      goto 100
c
c     fixed degress of freedom 
c     ----------------------------
      else if (inp_compare(.false.,'fix', test)) then
        if (.not. inp_a(test))
     >      call errquit('cons_input: failed to read keyword', 0,
     >        0)
         if (inp_compare(.false.,'atom', test)) then
           status = inp_ilist(max_fp,int_mb(i_tmp_id),nfp)
           if(status) then
             call cons_add_fixed_atoms(namespace,nfp,
     >                     int_mb(i_tmp_id),rtdb)
           else
             call errquit('cons_input: failed fixed atom input', 0,
     >        0)
           end if
         end if
c
      goto 100
c
c     unfixed degress of freedom 
c     ----------------------------
      else if (inp_compare(.false.,'unfix', test)) then
        if (.not. inp_a(test))
     >      call errquit('cons_input: failed to read keyword', 0,
     >        0)
         if (inp_compare(.false.,'atom', test)) then
           status = inp_ilist(max_fp,int_mb(i_tmp_id),nfp)
           if(status) then
             call cons_unfix_atoms(namespace,nfp,int_mb(i_tmp_id),rtdb)
           else
             call errquit('cons_input: failed unfixed atom input', 0,
     >        0)
           end if
         end if
c
      goto 100
c
c     spring contraints
c     ----------------------------
      else if (inp_compare(.false.,'spring', test)) then
        if (.not. inp_a(test))
     >      call errquit('cons_input: failed to read keyword', 0,
     >        0)
         if (inp_compare(.false.,'bond', test)) then
           if(cons_harm_bond_input(id(1),id(2),k,r)) then 
             call cons_add_spring(namespace,id(1),id(2),k,r,rtdb)
           else
             call errquit('cons_input: failed harm bond input', 0,
     >        0)
           end if
         end if
c 
      goto 100
c
c     end of input block
c     -------------------
      else if (inp_compare(.false.,'end', test)) then
       goto 101
      endif
c     
c
101   continue

      if (.not. ma_pop_stack(h_tmp_id) ) call errquit(
     &    'cons_input_springs: unable to free h_tmp_id',
     &    0, MA_ERR)
        
      call cons_info_fixed_atoms(namespace,rtdb)
      call cons_info_springs(namespace,rtdb)
      return
      end
c
      function cons_harm_bond_input(i,j,sk,sr)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      logical cons_harm_bond_input
      integer i,j
      double precision sk
      double precision sr

      integer itest
      double precision ftest

      if(.not.inp_i(itest))
     & call errquit('cons_harm_input: not an integer field',0,0)
      i=itest
      if(.not.inp_i(itest))
     & call errquit('cons_harm_input: not an integer field',0,0)
      j=itest
      if(.not.inp_f(ftest))
     & call errquit('cons_harm_input: not a real field',0,0)
      sk=ftest
      if(.not.inp_f(ftest))
     & call errquit('cons_harm_input: not a real field',0,0)
      sr=ftest
      
c     do some checks
      if(sr .le. 0 ) then
        cons_harm_bond_input=.false.
        return
      else if(sk .le. 0) then
        cons_harm_bond_input=.false.
        return
      else if(i.eq.j) then
        cons_harm_bond_input=.false.
        return
      end if

      if (j .lt. i) then
        itest  = j
        j = i
        i = itest
      end if 

      cons_harm_bond_input = .true.
      return
      end

