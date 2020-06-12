c $Id$
      subroutine cons_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
#include "cons.fh"
#include "geom.fh"
c     
      integer rtdb
c
      character*255 test
      character*255 namespace
      integer max_fp
      integer i_tmp_id,h_tmp_id
      integer nfp
      logical status
      integer id(4)
      double precision k,r,K0,gcut0,gamma0,c0,coef(20)
      integer n,indx(20)
      integer geom, ncenters
c
      logical  cons_harm_bond_input
      external cons_harm_bond_input
      logical  cons_harm_dihed_input
      external cons_harm_dihed_input
      logical  cons_harm_bondings_input
      external cons_harm_bondings_input
      logical  cons_penalty_bondings_input
      external cons_penalty_bondings_input
      logical  cons_active
      external cons_active
c
c     get no. atoms
c
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_input: geom_create?',70, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_input: geom_rtdb?',71, GEOM_ERR)
      if (.not. geom_ncent(geom, ncenters))
     &     call errquit('cons_input: geom_ncent failed',73, GEOM_ERR)
      if (.not. geom_destroy(geom))
     &     call errquit('cons_input: geom_destroy fails',7, GEOM_ERR)
      max_fp=ncenters*4
      if ( .not. ma_push_get( MT_INT, ncenters*4, 
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
         call cons_delete_hbondings(namespace,rtdb)
         call cons_delete_pbondings(namespace,rtdb)
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
             call errquit('cons_input: failed harm bond input',0,0)
           end if
         else if (inp_compare(.false.,'dihedral', test)) then
           write(*,*) "found dihedral"
           if(cons_harm_dihed_input(id,k,r)) then 
             call cons_add_dihed(namespace,id,k,r,rtdb)
           else
             call errquit('cons_input: failed harm dihedral input',0,0)
           end if
         else if (inp_compare(.false.,'bondings', test)) then
           if(cons_harm_bondings_input(n,indx,coef,K0,gamma0)) then
            call cons_add_bondings(namespace,n,indx,coef,K0,gamma0,rtdb)
           else
            call errquit('cons_input: failed harm bondings input',0,0)
           end if
         end if
c 
      goto 100

c
c     penalty functions
c     ----------------------------
      else if (inp_compare(.false.,'penalty', test)) then
        if (.not. inp_a(test))
     >      call errquit('cons_input: failed to read keyword', 0,
     >        0)
         if (inp_compare(.false.,'bondings', test)) then
           if(cons_penalty_bondings_input(n,indx,coef,
     >                                    K0,gcut0,gamma0,c0))
     >     then
            call cons_add_pbondings(namespace,n,indx,coef,K0,
     >                              gcut0,gamma0,c0,rtdb)
           else
            call errquit('cons_input:failed penalty bondings input',0,0)
           end if
         end if

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
      call cons_info_bondings(namespace,rtdb)
      call cons_info_pbondings(namespace,rtdb)
      call cons_info_dihed(namespace,rtdb)
      return
      end
c
      function cons_harm_dihed_input(id,sk,sr)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      logical cons_harm_dihed_input
      integer id(4)
      double precision sk
      double precision sr

      integer itest
      double precision ftest
      character*30 pname
      character*72 message
      integer k,m

      pname = "cons_harm_dihed_input"

      do k=1,4
        if(.not.inp_i(id(k)))
     &   call errquit(pname//': not an integer field',0,0)
      end do
      if(.not.inp_f(ftest))
     & call errquit(pname//': not a real field',0,0)
      sk=ftest
      if(.not.inp_f(ftest))
     & call errquit(pname//': not a real field',0,0)
      sr=ftest*deg_to_rad
      
c     do some checks
      do k=1,4
       do m=k+1,4
         if(id(m).eq.id(k)) then
           message = " index input error "
           goto 10
         end if
       end do
      end do
      if(sk .le. 0 ) then
        cons_harm_dihed_input=.false.
        message = " use positive K "
        goto 10 
      end if
      cons_harm_dihed_input = .true.
      return
c     
c     error processing
c     ---------------
10    continue
      call errquit(pname//message,0,0)
      end
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


*     ****************************************
*     *                                      *
*     *       cons_harm_bondings_input       *
*     *                                      *
*     ****************************************
*
* bondings K0 gamma0 -1.0 1 2 1.0 2 3
*
      logical function cons_harm_bondings_input(n,indx,coef,K0,gamma0)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      logical ok
      integer n,indx(*)
      double precision coef(*),K0,gamma0

      integer itest
      double precision ftest

      ok = .true.

      !*** enter K0 and gamma0 ***
      if (.not.inp_f(ftest)) ok = .false.
      K0 = ftest
      if (ok.and.(.not.inp_f(ftest))) ok = .false.
      gamma0 = ftest

      !*** coef and indx ***
      n = 0
      do while (ok.and.inp_f(ftest))
         coef(n+1) = ftest
         if(ok.and.(.not.inp_i(itest))) ok = .false.
         indx(2*n+1) = itest
         if(ok.and.(.not.inp_i(itest))) ok = .false.
         indx(2*n+2) = itest
         n = n + 1
      end do
      if (n.eq.0) ok = .false.

      cons_harm_bondings_input = ok
      return
      end




*     ****************************************
*     *                                      *
*     *       cons_penalty_bondings_input    *
*     *                                      *
*     ****************************************
*
* penalty bondings K0 gcut0 gamma0 c0 -1.0 1 2 1.0 2 3
*
      logical function cons_penalty_bondings_input(n,indx,coef,K0,
     >                                             gcut0,gamma0,c0)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      logical ok
      integer n,indx(*)
      double precision coef(*),K0,gcut0,gamma0,c0

      integer itest
      double precision ftest

      ok = .true.

      !*** enter K0, gcut0, gamma0, and c0 ***
      if (.not.inp_f(ftest)) ok = .false.
      K0 = ftest
      if (ok.and.(.not.inp_f(ftest))) ok = .false.
      gcut0 = ftest
      if (ok.and.(.not.inp_f(ftest))) ok = .false.
      gamma0 = ftest
      if (ok.and.(.not.inp_f(ftest))) ok = .false.
      c0 = ftest

      !*** coef and indx ***
      n = 0
      do while (ok.and.inp_f(ftest))
         coef(n+1) = ftest
         if(ok.and.(.not.inp_i(itest))) ok = .false.
         indx(2*n+1) = itest
         if(ok.and.(.not.inp_i(itest))) ok = .false.
         indx(2*n+2) = itest
         n = n + 1
      end do
      if (n.eq.0) ok = .false.

      cons_penalty_bondings_input = ok
      return
      end


