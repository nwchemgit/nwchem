c
c $Id: frac_occ.F,v 1.2 2004-03-08 22:51:25 bylaska Exp $
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

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "frac_occ.fh"


*     **** local variables ****
      character*30 rtdb_name
      integer i,j
      integer num_states,occ(2),s(2)


      fractional = .false.

      rtdb_name = 'nwpw:frac_occ:number_states'
      if (.not. rtdb_get(rtdb,rtdb_name,mt_int,1,num_states))
     >  num_states = 0

      if (num_states.gt.0) then
        fractional = .true.

*       **** allocate from heap memory ****
        if (.not.MA_alloc_get(mt_dbl,(ne(1)+ne(2)),
     >                        'fweight',fweight(2),fweight(1)))
     >     call errquit(
     >          'frac_occ_init: heap failed 1', 0,MA_ERR)

*       **** allocate from stack memory ****
        if((.not.MA_push_get(mt_dbl,(num_states),'occ',occ(2),occ(1)))
     > .or.(.not. MA_push_get(mt_int,(num_states),'s',s(2),s(1))))
     >     call errquit(
     >          'frac_occ_init: stack failed 1', 0,MA_ERR)


*       **** get info from rtdb ****
        rtdb_name = 'nwpw:frac_occ:occupations'
        if (.not. rtdb_get(rtdb,rtdb_name,mt_dbl,
     >                   (num_states),
     >                    dbl_mb(occ(1))))
     >    call errquit(
     >    'frac_occ_init: rtdb_get(occupation) failed', 0,RTDB_ERR)
        rtdb_name = 'nwpw:frac_occ:states'
        if (.not. rtdb_get(rtdb,rtdb_name,mt_int,
     >                   (num_states),
     >                    int_mb(s(1))))
     >    call errquit(
     >    'frac_occ_init: rtdb_get(states) failed', 0,RTDB_ERR)


*       **** set the fractional weights ****
        call dcopy((ne(1)+ne(2)),1.0d0,0,dbl_mb(fweight(1)),1)
        do j=1,num_states
          i = int_mb(s(1)+j-1) 
          dbl_mb(fweight(1)+i-1) = dbl_mb(occ(1)+j-1)
        end do


*       **** pop stack memory ****
        if ((.not.MA_pop_stack(s(2))).or.(.not.MA_pop_stack(occ(2))))
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "frac_occ.fh"

      if (fractional) then
      if (.not.MA_free_heap(fweight(2)))
     >  call errquit('frac_occ_end: heap failed',0,MA_ERR)
      end if
      return
      end
