c
c $Id$
c

*     *************************************************
*     *                                               *
*     *                frac_occ_input                 *
*     *                                               *
*     *************************************************
      subroutine frac_occ_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 8)
      character*22 dirs(num_dirs)
      data dirs / 'occ', 'filling', 'occupation',
     >            'occ:','filling:','occupation:',
     >            'extra_orbitals',
     >            'end'/


      character*255     test
      character*30     rtdb_name
      integer          num_states,k,l
      double precision socc
      integer          state,norbs
      integer          s(2),  s_new(2)
      integer          occ(2),occ_new(2)


c     
c     **** intialize stuff ****
      num_states = 0
      norbs = 0
c
c
 10   if (.not. inp_read()) 
     >     call errquit(
     >           'frac_occ_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'frac_occ_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'frac_occ_input: unknown directive', 0,
     &       INPUT_ERR)


      goto ( 200,200,200,200,200,200,300,
     >      9999) ind
      call errquit(
     >      'frac_occ_input: unimplemented directive', ind,
     &       INPUT_ERR)


c
c occupation
c
 200  if (.not. inp_f(socc))  socc  = 1.0d0
      if (.not. inp_i(state)) state = 1
     
      num_states = num_states + 1
      if  ((.not. BA_alloc_get(mt_dbl,(num_states), 
     >        'occ_new',occ_new(2),occ_new(1)))
     > .or.(.not. BA_alloc_get(mt_int,(num_states), 
     >        's_new',s_new(2),s_new(1))))
     >     call errquit(
     >          'frac_occ_input: heap failed 1', 0,
     &       INPUT_ERR)

      if (num_states.gt.1) then
        do k=1,(num_states-1)
         dbl_mb(occ_new(1) + (k-1)) = dbl_mb(occ(1) + (k-1))
         int_mb(s_new(1)   + (k-1)) = int_mb(s(1)   + (k-1))

        end do
        if ((.not.BA_free_heap(s(2))).or.(.not.BA_free_heap(occ(2))))
     >     call errquit(
     >          'frac_occ_input: heap failed 2', 0, MA_ERR)

      end if
      dbl_mb(occ_new(1) + (num_states-1)) = socc
      int_mb(s_new(1)   + (num_states-1)) = state
      occ(1) = occ_new(1)
      occ(2) = occ_new(2)
      s(1)   = s_new(1)
      s(2)   = s_new(2)


      goto 10

c
c extra_orbitals
c
 300  if (.not. inp_i(norbs))  norbs = 1

      rtdb_name = 'nwpw:frac_occ:extra_orbitals'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,norbs))
     >  call errquit(
     >  'frac_occ_input: rtdb_put(extra_orbitals) failed', 0,
     >       RTDB_ERR)

      goto 10

*     ***** add frac_occ to rtdb ****
 9999 continue


      rtdb_name = 'nwpw:frac_occ:number_states'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_states))
     >  call errquit(
     >  'frac_occ_input: rtdb_put(number_states) failed', 0,
     >       RTDB_ERR)

      rtdb_name = 'nwpw:frac_occ:occupations'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                   (num_states),
     >                    dbl_mb(occ(1))))
     >  call errquit(
     >  'frac_occ_input: rtdb_put(occupations) failed', 0,
     >       RTDB_ERR)

      rtdb_name = 'nwpw:frac_occ:states'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,
     >                   (num_states),
     >                    int_mb(s(1))))
     >  call errquit(
     >  'frac_occ_input: rtdb_put(states) failed', 0,
     >       RTDB_ERR)

      if ((.not.BA_free_heap(s(2))).or.(.not.BA_free_heap(occ(2))))
     >  call errquit(
     >       'frac_occ_input: heap failed 3', 0, MA_ERR)

      return
      end

