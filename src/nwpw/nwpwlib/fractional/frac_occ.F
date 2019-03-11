c
c $Id$
c

*     *************************************************
*     *                                               *
*     *                frac_occ_init                  *
*     *                                               *
*     *************************************************
      subroutine frac_occ_init(rtdb,ispin,ne)
      implicit none
      integer rtdb
      integer ispin
      integer ne(2)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "frac_occ.fh"


*     **** local variables ****
      character*50 rtdb_name
      integer i,j
      integer num_states,occ(2),s(2)


      fractional = .false.

      rtdb_name = 'nwpw:frac_occ:number_states'
      if (.not. btdb_get(rtdb,rtdb_name,mt_int,1,num_states))
     >  num_states = 0

      if (num_states.gt.0) then
        fractional = .true.

*       **** allocate from heap memory ****
        if (.not.BA_alloc_get(mt_dbl,(ne(1)+ne(2)),
     >                        'fweight',fweight(2),fweight(1)))
     >     call errquit(
     >          'frac_occ_init: heap failed 1', 0,MA_ERR)

*       **** allocate from stack memory ****
        if((.not.BA_push_get(mt_dbl,(num_states),'occ',occ(2),occ(1)))
     > .or.(.not. BA_push_get(mt_int,(num_states),'s',s(2),s(1))))
     >     call errquit(
     >          'frac_occ_init: stack failed 1', 0,MA_ERR)


*       **** get info from rtdb ****
        rtdb_name = 'nwpw:frac_occ:occupations'
        if (.not. btdb_get(rtdb,rtdb_name,mt_dbl,
     >                   (num_states),
     >                    dbl_mb(occ(1))))
     >    call errquit(
     >    'frac_occ_init: btdb_get(occupation) failed', 0,RTDB_ERR)
        rtdb_name = 'nwpw:frac_occ:states'
        if (.not. btdb_get(rtdb,rtdb_name,mt_int,
     >                   (num_states),
     >                    int_mb(s(1))))
     >    call errquit(
     >    'frac_occ_init: btdb_get(states) failed', 0,RTDB_ERR)


*       **** set the fractional weights ****
        call dcopy((ne(1)+ne(2)),1.0d0,0,dbl_mb(fweight(1)),1)
        do j=1,num_states
          i = int_mb(s(1)+j-1) 
          dbl_mb(fweight(1)+i-1) = dbl_mb(occ(1)+j-1)
        end do


*       **** pop stack memory ****
        if ((.not.BA_pop_stack(s(2))).or.(.not.BA_pop_stack(occ(2))))
     >  call errquit(
     >       'frac_occ_init: stack failed', 0, MA_ERR)

      end if

      return
      end

*     *************************************************
*     *                                               *
*     *                frac_occ_end                   *
*     *                                               *
*     *************************************************
      subroutine frac_occ_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "frac_occ.fh"

      if (fractional) then
      if (.not.BA_free_heap(fweight(2)))
     >  call errquit('frac_occ_end: heap failed',0,MA_ERR)
      end if
      return
      end

*     *************************************************
*     *                                               *
*     *                frac_occ_set                   *
*     *                                               *
*     *************************************************

*     *************************************************
      subroutine frac_occ_set(rtdb,ispin,ne,occ)
      implicit none
      integer rtdb
      integer ispin
      integer ne(2)
      real*8 occ(*)

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"


*     **** local variables ****
      logical      assending
      character*50 rtdb_name
      integer i,j
      integer num_states,focc(2),s(2)

      rtdb_name = 'nwpw:frac_occ:number_states'
      if (.not. btdb_get(rtdb,rtdb_name,mt_int,1,num_states))
     >  num_states = 0

      if (num_states.gt.0) then
         assending =  occ(1).gt.occ(ne(1))

*        **** allocate from stack memory ****
        if((.not.BA_push_get(mt_dbl,(num_states),'occ',focc(2),focc(1)))
     >  .or.(.not. BA_push_get(mt_int,(num_states),'s',s(2),s(1))))
     >     call errquit(
     >          'frac_occ_set: stack failed 1', 0,MA_ERR)

*        **** get info from rtdb ****
         rtdb_name = 'nwpw:frac_occ:occupations'
         if (.not. btdb_get(rtdb,rtdb_name,mt_dbl,
     >                     (num_states),
     >                      dbl_mb(focc(1))))
     >    call errquit('frac_occ_set: btdb_get(occupation) failed',
     >                 0,RTDB_ERR)
         rtdb_name = 'nwpw:frac_occ:states'
         if (.not. btdb_get(rtdb,rtdb_name,mt_int,
     >                     (num_states),
     >                      int_mb(s(1))))
     >    call errquit('frac_occ_set: btdb_get(states) failed',
     >                 0,RTDB_ERR)


*        **** reset the fractional weights ****
         if (assending) then
            do j=1,num_states
               i = int_mb(s(1)+j-1)
               occ(ne(1)+ne(2) - i + 1) = dbl_mb(focc(1)+j-1)
            end do
         else
            do j=1,num_states
               i = int_mb(s(1)+j-1)
               occ(i) = dbl_mb(focc(1)+j-1)
            end do
         end if

*        **** pop stack memory ****
         if ((.not.BA_pop_stack(s(2))).or.(.not.BA_pop_stack(focc(2))))
     >   call errquit(
     >        'frac_occ_set: stack failed', 0, MA_ERR)

      end if

      return
      end

