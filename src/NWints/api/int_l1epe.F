      subroutine int_l1epe(i_basis, ish, j_basis, jsh, zerotol,
     &       ilab, jlab, lpe, Pe, lscr, scr, numgen) 
c $Id: int_l1epe.F,v 1.4 1996-07-18 23:58:08 d3e129 Exp $
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
      integer i_basis, j_basis, ish, jsh, lpe, numgen, lscr
      integer ilab(lpe),jlab(lpe)
      double precision zerotol, scr(lscr)
      double precision Pe(lpe)
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
     &       call errquit('int_l1epe: invalid contraction/shell',0)
c
c  check if spherical
c
      call int_nospherical_check(i_basis,'int_l1epe')
      call int_nospherical_check(j_basis,'int_l1epe')
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
* compute pe integrals
      call int_1epe(i_basis, ish, j_basis, jsh,
     &    newlscr, scr(numint+1), numint,scr)
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
              Pe(numgen) = scr(icount)
              ilab(numgen) = i
              jlab(numgen) = j
            endif
          endif
c
00200   continue
00100 continue
c
      end
