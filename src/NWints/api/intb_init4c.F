      logical function intb_init4c(brain, icl, jcl, ketin, kcl, lcl,
     &    num_q, lscr, scr, block_eff)
c $Id: intb_init4c.F,v 1.1 1995-05-09 18:01:01 d3e129 Exp $  
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
      integer lscr       ! [input] length of scratch array
      double precision scr(lscr)  ! [input] scratch array
      double precision block_eff  ! [output] blocking efficiency parameter
#include "numb_qP.fh"
      integer i_q, ii, jj, kk, ll
c
      num_q_start = 1
      num_q_last  = -1
c
      do 00100 i_q = 1,num_q
        ii = icl(i_q)
        jj = jcl(i_q)
        kk = kcl(i_q)
        ll = lcl(i_q)
        call int_canon(ii,jj,kk,ll,
     &      icl(i_q),jcl(i_q),kcl(i_q),lcl(i_q))
00100 continue
      intb_init4c = .true.
      block_eff = 0.0d00
      end
