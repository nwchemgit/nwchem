c $Id: cons_input.F,v 1.1 2004-01-28 01:30:59 marat Exp $
      subroutine cons_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb
      character*255 test

      write(*,*) "found cons"
c
 100  if (.not. inp_read()) call errquit('cons_input: inp_read failed',0,
     &       0)
c
      if (.not. inp_a(test))
     $     call errquit('cons_input: failed to read keyword', 0,
     &       0)
c     
c     print level is read
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, "cons")
c
      goto 100
c
c     fixed degress of freedom 
c     ----------------------------
c     else if (inp_compare(.false.,'fixed', test)) then
c        call cons_fixed_input(rtdb)
  
      goto 100
c
c     harmonic contraints
c     ----------------------------
      else if (inp_compare(.false.,'harm', test)) then
         call cons_harm_input(rtdb)
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
      
      return
      end
c
c
c
c
      subroutine cons_fixed_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb
      character*255 test
      integer itest
      integer ncent
      integer h_fixed,i_fixed,n_fixed
      logical status
c
c     check if we have integers here
c     ------------------------------
      if(.not.inp_i(itest))
     &  call errquit('cons_fixed_input: not an integer field',0,0) 
      call inp_prev_field()
c
c     allocate memory for fixed atom list 
      status = inp_ilist_size(ncent)
      if (.not. status) call errquit(
     &    'cons_fixed_input: failed reading list ', 0, INPUT_ERR)
      if ( .not. ma_push_get( MT_INT, ncent, 'input scratch',
     &      h_fixed, i_fixed) ) call errquit(
     &      'cons_fixed_input: unable to allocate input scratch space',
     &      ncent, MA_ERR)
c
c     read a list of fixed atoms
c     --------------------------
      status = inp_ilist(ncent, int_mb(i_fixed), n_fixed)
      if (.not. status) call errquit(
     &    'cons_fixed_input: failed reading list ', 0, INPUT_ERR)

      status = rtdb_put(rtdb, 'cons:fixed_atoms', MT_INT, n_fixed,
     &      int_mb(i_fixed))
      if (.not. status) call errquit(
     &    'cons_fixed_input: failed storing list ', 0, RTDB_ERR)

      if (.not. ma_pop_stack(h_fixed) ) call errquit(
     &    'cons_fixed_input: unable to free scratch space',
     &    h_fixed, MA_ERR)

      return
      end
c
c
      subroutine cons_harm_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cons_data.fh"
c     
      integer rtdb
      character*255 test
      integer itest
      double precision ftest
      integer space_avail
      integer i
      integer ma_type
      logical status

      write(*,*) "found harm bonds"
c
c     allocate scratch memory buffer for bonds
c     ------------------------------
      space_avail = min( max_harm_bonds*4, ma_inquire_avail(MT_DBL) )
      if ( .not. ma_push_get( MT_DBL, space_avail, 'cons scratch',
     &      h_bond, i_bond) ) call errquit(
     &      'input_set: unable to allocate cons scratch space',
     &      space_avail, MA_ERR)

c
c     allocate scratch memory buffer for angles
c     ------------------------------
      space_avail = min( max_harm_angles*4, ma_inquire_avail(MT_DBL) )
      if ( .not. ma_push_get( MT_DBL, space_avail, 'cons scratch',
     &      h_angle, i_angle) ) call errquit(
     &      'input_set: unable to allocate cons scratch space',
     &      space_avail, MA_ERR)

c
      n_bond=0
      n_angle=0
100   if (.not. inp_read()) 
     & call errquit('cons_harm_input: inp_read failed',0,
     &       0)
c
      if(.not.inp_a(test))
     &  call errquit('cons_harm_input: not a string field',0,0) 

      if (inp_compare(.false.,'bond', test)) then
           i=i_bond+n_bond*4
           if(.not.inp_i(itest))
     &     call errquit('cons_harm_input: not an integer field',0,0)
           dbl_mb(i)=itest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           if(.not.inp_i(itest))
     &     call errquit('cons_harm_input: not an integer field',0,0)
           dbl_mb(i)=itest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           if(.not.inp_f(ftest))
     &     call errquit('cons_harm_input: not a real field',0,0)
           dbl_mb(i)=ftest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           if(.not.inp_f(ftest))
     &     call errquit('cons_harm_input: not a real field',0,0)
           dbl_mb(i)=ftest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           n_bond=n_bond+1
           goto 100
c
      else if (inp_compare(.false.,'angle', test)) then
           i=i_angle+n_angle*4
           if(.not.inp_i(itest))
     &     call errquit('cons_harm_input: not an integer field',0,0)
           dbl_mb(i)=itest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           if(.not.inp_i(itest))
     &     call errquit('cons_harm_input: not an integer field',0,0)
           dbl_mb(i)=itest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           if(.not.inp_i(itest))
     &     call errquit('cons_harm_input: not an integer field',0,0)
           dbl_mb(i)=itest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           if(.not.inp_f(ftest))
     &     call errquit('cons_harm_input: not a real field',0,0)
           dbl_mb(i)=ftest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           if(.not.inp_f(ftest))
     &     call errquit('cons_harm_input: not a real field',0,0)
           dbl_mb(i)=ftest
           write(*,*) "dbl_mb(i)",dbl_mb(i)
           i=i+1
           n_angle=n_angle+1
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

        status = rtdb_put(rtdb, 'harm_cons_test', mt_log, 1, .true.)
      if(n_bond.ne.0 .or. n_angle.ne.0) then
        write(*,*) "constraints are inthe database"
        status = rtdb_put(rtdb, 'harm_cons', mt_log, 1, .true.)
        if (.not. status) call errquit(
     &    'cons_harm_input: failed storing  ', 0, RTDB_ERR)
      end if

      if(rtdb_get(rtdb, 'harm_cons', mt_log, 1, status)) then
        write(6,*) "found constraints"
      else
        write(6,*) "not found constraints"
      end if


      if(n_bond.ne.0) then
        status = rtdb_put(rtdb, 'cons:harm_bonds', MT_DBL, n_bond*4,
     &      dbl_mb(i_bond))
        if (.not. status) call errquit(
     &    'cons_harm_input: failed storing list ', 0, RTDB_ERR)
      end if

      if(n_angle.ne.0) then
        status = rtdb_put(rtdb, 'cons:harm_angles', MT_DBL, n_angle*5,
     &      dbl_mb(i_angle))
        if (.not. status) call errquit(
     &    'cons_harm_input: failed storing list ', 0, RTDB_ERR)
      end if


      if (.not. ma_pop_stack(h_angle) ) call errquit(
     &    'cons_bond_input: unable to free scratch space',
     &    h_bond, MA_ERR)

      if (.not. ma_pop_stack(h_bond) ) call errquit(
     &    'cons_bond_input: unable to free scratch space',
     &    h_bond, MA_ERR)

      return
      end
