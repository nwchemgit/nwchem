c     fractional occupation handler
c
      subroutine occup_input(rtdb)
c
C     $Id$

c     occup 5 3
c      1.0 1.0
c      1.0 1.0
c      1.0 1.0
c      1.0
c      1.0
c     end
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
#include "global.fh"
#include "rtdb.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      integer rtdb              ! [input]     
      character*255 field       ! for character input
      integer ac,noccA,noccB    ! counts no. of centers as input
      logical status            ! scratch for return codes
      integer i,ind_min,ind_max
      integer nfocc_in(2),nfocc ! input of nfocc
      integer l_occ,k_occ
      integer switch_focc
c     
c     read occup from the input deck and persist it in the rtdb
c     
      if (ga_nodeid() .ne. 0) return
c     
c     Check that this is indeed a occupations line
c     
      call inp_set_field(0)     ! goto start of line
      if (.not. inp_a(field))
     $ call errquit('occup_input: occup input error', 0, INPUT_ERR)
      if (.not. inp_compare(.false., 'occup', field))
     $ call errquit('occup_input: occup input error', 0,INPUT_ERR)   
c
      nfocc_in(1) = 0  ! spin-1
      nfocc_in(2) = 0  ! spin-2
      status = .true.
      status = status .and. inp_i(nfocc_in(1))
      status = status .and. inp_i(nfocc_in(2))
c
      if (status) then
        nfocc=nfocc_in(1)+nfocc_in(2) ! total fractional occupation
        ind_min=2
        ind_max=1   ! spin-1 is dominant
        if (nfocc_in(1).lt.nfocc_in(2)) then
          ind_min=1
          ind_max=2  ! spin-2 is dominant
        endif
c
c       set the memory based on the number of fractional occupied states
        if (.not. ma_push_get(mt_dbl,nfocc,'occup', l_occ,k_occ))
     &    call errquit('occ_input ma_push_get fail k_occ',0,MA_ERR)  

c----- Start reading <occ-A> <occ-B> sets
       noccA=1
       noccB=1
       ac = 0
 20    if (inp_read()) then
        status = inp_a(field)
        if (inp_compare(.false., 'end', field)) then
          goto 30
        else
          if (noccA .gt. nfocc_in(ind_max)) 
     &        call errquit('occup_input: too many occupations?', 
     &            ac, INPUT_ERR)
c        
c           input format  <occ-A>  <occ-B> 
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
        call errquit('occup_input: premature end of file',0,INPUT_ERR)
       end if
      else
       call errquit('occup_input: missing alpha and beta occupations',
     &         0,INPUT_ERR) 
      endif  ! status
c
 30   continue    
c
c     consistency check
      if ((noccA-1 .ne. nfocc_in(1)) .or. (noccB-1.ne.nfocc_in(2))) then
         write(luout,*) "noccA-1: ",noccA-1
         write(luout,*) "nfocc_in(1): ",nfocc_in(1)
         write(luout,*) "noccB-1: ",noccB-1
         write(luout,*) "nfocc_in(2): ",nfocc_in(2)
         call errquit('occup_input: occupancy mismatch',0,INPUT_ERR)
      end if ! consistency check
c
c     store in rtdb 
      switch_focc=1 ! 1 means using occup keyword
c
      if (.not. rtdb_put(rtdb,'focc:occup_switch',mt_int,1,switch_focc))
     &  call errquit('occup_input: rtdb_put failed', 1201, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'focc:occup',mt_int,2,nfocc_in))
     &  call errquit('occup_input: rtdb_put failed', 1201, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'focc:occup_list',
     &    mt_dbl,nfocc,dbl_mb(k_occ)))
     &  call errquit('occup_input: rtdb_put failed', 0, RTDB_ERR)
      if (.not.ma_pop_stack(l_occ)) ! Free memory
     &      call errquit('occ input: ma_pop_stack failed k_occ',
     &                    0,MA_ERR)
c
50    continue
      end
