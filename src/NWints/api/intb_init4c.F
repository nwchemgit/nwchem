      logical function intb_init4c(brain, icl, jcl, ketin, kcl, lcl,
     &    num_q, lscr, scr, l_erilab, block_eff)
c $Id: intb_init4c.F,v 1.5 1996-09-03 06:02:50 d3e129 Exp $  
      implicit none
*
c
c blocking initialization routien for 2 electron 4 center integrals
c
*
      integer brain      ! [input] basis set handle for bra basis
      integer ketin      ! [input] basis set handle for ket basis
      integer num_q      ! [input] number of quartets
      integer icl(num_q) ! [input] i-contraction labels for quartets
      integer jcl(num_q) ! [input] j-contraction labels for quartets
      integer kcl(num_q) ! [input] k-contraction labels for quartets
      integer lcl(num_q) ! [input] l-contraction labels for quartets
      integer l_erilab   ! [input] size of eri and label arrays that will 
*................................. be used in intb_2e4c.F
      integer lscr       ! [input] length of scratch array
      double precision scr(lscr)  ! [input] scratch array
      double precision block_eff  ! [output] blocking efficiency parameter
#include "numb_qP.fh"
c
      logical status
      integer i_q
c
      logical cando_txs
      external cando_txs
c
      num_q_start = 1
      num_q_last  = -1
c
      more_ints_txs = .false. 
      intb_init4c = .true.
      block_eff = -1.0d00
*
*.... check shells for txs only after intb_init is called and established        
      txs_b2e4c = .true.
      status = .true.
*.... check first 20 quartets but every other one
      do i_q = 1,min(20,num_q),2
        status = status.and. (cando_txs(brain,icl(i_q),jcl(i_q)))
        status = status.and. (cando_txs(ketin,kcl(i_q),lcl(i_q)))
        if (.not.status) then
          txs_b2e4c = .false.
          return
        endif
      enddo
      end
