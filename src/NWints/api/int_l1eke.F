      subroutine int_l1eke(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, lke, Ke, lscr, scr, numgen) 
c $Id: int_l1eke.F,v 1.3 1996-01-22 21:46:28 d3e129 Exp $
      implicit none
#include "basP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c
c::external subroutines used
c errquit
c::function
      logical int_chk_sh
      external int_chk_sh
c::passed
      integer i_basis, j_basis, ish, jsh, lke, numgen, lscr
      integer ilab(lke),jlab(lke)
      double precision zerotol, scr(lscr)
      double precision Ke(lke)
c::local
      integer ibas, jbas, icount, i, j
      logical ijbas
      logical shells_ok
      logical nonzero
      integer numint, newlscr
c
#include "bas_ibs_sfn.fh"
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_l1eke: invalid contraction/shell',0)
c      
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
      ijbas = ibas .eq. jbas
c
      numint =         (sf_ibs_cn2bfr(2,ish,ibas))-
     &      (sf_ibs_cn2bfr(1,ish,ibas))+1
      numint = numint*((sf_ibs_cn2bfr(2,jsh,jbas))-
     &      (sf_ibs_cn2bfr(1,jsh,jbas))+1)
      newlscr = lscr - numint
      call int_1eke(i_basis, ish, j_basis, jsh,
     &    newlscr,scr(numint+1),numint,scr)
c
      numgen = 0
      icount = 0
      do 00100 i = (sf_ibs_cn2bfr(1,ish,ibas)),
     &      (sf_ibs_cn2bfr(2,ish,ibas))
        do 00200 j = (sf_ibs_cn2bfr(1,jsh,jbas)),
     &        (sf_ibs_cn2bfr(2,jsh,jbas))
          icount = icount + 1
          nonzero = abs(scr(icount)).ge.zerotol
          if (nonzero) then
            if ((.not.ijbas).or.i.ge.j) then
              numgen = numgen + 1
              Ke(numgen) = scr(icount)
              ilab(numgen) = i
              jlab(numgen) = j
            endif
          endif
c
00200   continue
00100 continue
c
      end
