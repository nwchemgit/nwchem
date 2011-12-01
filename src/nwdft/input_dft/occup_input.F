c     fractional occupation handler
c
      subroutine occup_input(rtdb)
c
C     $Id$
C     Adapted from geom_input
c     Occupations format:
c     occupations [nfocc_in]
c     [focc-1]
c     [focc-2]
c      ...
c     [focc-nfocc_in]
c     end
c     For example:
c     occupations 5
c     0.25
c     0.1
c     1.7
c     1.1
c     0.175
c     end
c     RESTRICTIONS:
c     2. [focc-1]   is integer or float and positive
c     3. [nfocc_in] is integer and positive and < nmo*2
c        Note.- In so-dft the MOs are arranged as 
c        {nmo-alpha-1 nmo-beta-1}{nmo-alpha-2 nmo-beta-2} ...
c        As an example: In CH2 which is triplet we have for 6-311G*
c        nmo=25
c        Default occupations are: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
c        This means 8 pairs of 1's.
c        This would be defined in our occupation module scheme as
c        This scheme is designed for so-dft calculation.
c        Note.- For other type of calculation is not done yet.
c        occupations 5 3
c        1.0 1.0
c        1.0 1.0
c        1.0 1.0
c        1.0
c        1.0
c        end
c     FA-02-10-11
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "global.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer rtdb              ! [input]     
      character*255 field       ! for character input
      integer ac,noccA,noccB    ! counts no. of centers as input
      logical status            ! scratch for return codes
      integer i,ind_min,ind_max
      integer nfocc_in(2),nfocc ! input of nfocc
      integer l_occ,k_occ
      integer switch_focc
c     
c     read occupations from the input deck 
c     and output it to the rtdb.
c     
c     current input line should begin 'occupations ...'
c     
      if (ga_nodeid() .ne. 0) return
c     
c     Check that this is indeed a occupations line
c     
      call inp_set_field(0)     ! goto start of line
      if (.not. inp_a(field))
     $     call errquit('occup_input: no input present', 0, INPUT_ERR)
      if (.not. inp_compare(.false., 'occup', field))
     $     call errquit('occup_input: not occup input',
     &                   0,INPUT_ERR)   
      status = .true.
      status = status .and. inp_i(nfocc_in(1))
      status = status .and. inp_i(nfocc_in(2))
      if (status) then ! --------nfocc-in---START      
c       write(*,2) nfocc_in(1),nfocc_in(2)
c2      format('nfocc_in=(',i5,',',i5,')')  
        nfocc=nfocc_in(1)+nfocc_in(2)
        ind_min=2
        ind_max=1
        if (nfocc_in(1).lt.nfocc_in(2)) then
         ind_min=1
         ind_max=2
        endif
        if (.not. ma_push_get(mt_dbl,nfocc,'occup', l_occ,k_occ))
     &  call errquit('occ_input ma_push_get fail k_occ',
     &               0,MA_ERR)  
c----- Start reading <occ-A> <occ-B> sets
       noccA=1
       noccB=1
       ac = 0
 20    if (inp_read()) then   ! --- if-read-focc-START
         status = inp_a(field)
         if (inp_compare(.false., 'end', field)) then
            goto 30
         else
            if (noccA .gt. nfocc_in(ind_max)) call errquit
     $           ('occup_input: too many occupations?', 
     &            ac, INPUT_ERR)
c        
c     style input ... <occ-A>  <occ-B> 
c     
            call inp_set_field(0)
            status = .true.
c ------- Case 1: Reading occ-A occ-B
          if (noccA .le. nfocc_in(ind_min)) then
            status = status .and. inp_f(dbl_mb(k_occ+ac))
            status = status .and. inp_f(dbl_mb(k_occ+ac+1))
            if (.not. status) call errquit
     $          ('occup_input: <occ>',ac+1,INPUT_ERR)
            ac = ac + 2
            noccA=noccA+1
            noccB=noccB+1
          else if (ac/2+1 .gt. nfocc_in(ind_min)) then
c ------- Case 2: Reading occ-A (unpaired electrons)
            status = status .and. inp_f(dbl_mb(k_occ+ac))
            if (.not. status) call errquit
     $          ('occup_input: <occ>',ac+1,INPUT_ERR)
            if (ind_max.eq.1) noccA=noccA+1
            if (ind_max.eq.2) noccB=noccB+1
            ac = ac + 1
          endif
         end if
         goto 20
       else
         call errquit('occup_input: premature end of file',
     &                0,INPUT_ERR)
       end if                 ! --- if-read-focc-END
      else
       call errquit(
     &        'occup_input: missing nfocc-A -B in occupations',
     &         0,INPUT_ERR) 
      endif ! --------nfocc-in---END    
 30   continue    
c     write(*,5) noccA-1,noccB-1
c5    format('(noccA,noccB)=(',i5,',',i5,')')
      if (noccA-1.ne.nfocc_in(1) .or. noccB-1.ne.nfocc_in(2)) then
         call errquit('occup_input: mismatch noccA-ith or noccB-ith',
     &                0,INPUT_ERR)      
      endif
c
c ---- check that int_mb(k_nocc) > 0 and dbl_mb(k_occ) > 0
c      do i=1,nfocc
c       if (dbl_mb(k_occ+i-1).le.0.0d0) then
c        call errquit('occup_input: Problem focc input ith val<0',
c     &                0,INPUT_ERR)    
c       endif 
c      enddo
c
c     store in rtdb 
      switch_focc=1 ! 1 means using occupations keyword
      if (.not. rtdb_put(rtdb,'focc:occup_switch',
     &                   mt_int,1,switch_focc))  
     &   call errquit('occup_input: rtdb_put failed', 1201, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'focc:occup',mt_int,2,nfocc_in))
     &   call errquit('occup_input: rtdb_put failed', 1201, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'focc:occup_list',
     &                   mt_dbl,nfocc,dbl_mb(k_occ)))
     $      call errquit('occup_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not.ma_pop_stack(l_occ)) ! Free memory
     &      call errquit('occ input: ma_pop_stack failed k_occ',
     &                    0,MA_ERR)
c
50    continue
      end
