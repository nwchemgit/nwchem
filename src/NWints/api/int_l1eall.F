      subroutine int_l1eall(i_basis, j_basis, ish, jsh, zerotol,
     &       ilab, jlab, lstv, S, T, V, lscr, scr, numstv) 
c $Id: int_l1eall.F,v 1.3 1994-04-09 15:46:27 d3e129 Exp $
#include "basP.fh"
#include "geobasmapP.fh"
      implicit none
c::external subroutines used
c dcopy
c errquit
c::function
      logical int_chk_init
      external int_chk_init
      logical int_chk_sh
      external int_chk_sh
c::passed
      integer i_basis, j_basis, ish, jsh, lstv, numstv, lscr
      integer ilab(lstv),jlab(lstv)
      double precision zerotol, scr(lscr)
      double precision S(lstv), T(lstv), V(lstv)
c::local
      integer ibas, jbas, iscrS, iscrT, iscrV, icount, i, j
      logical ijbas, shells_ok
      logical nonzero
c
c check initialization and shells
c
      if (.not.int_chk_init('int_l1eall'))
     &       call errquit('int_l1eall: int_init was not called' ,0)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_l1eall: invalid contraction/shell',0)
c
c check canoncialization of shells
c
      ijbas = i_basis .eq. j_basis
      if (.not.(ijbas.and.ish.ge.jsh)) then
        write(6,*)' shells not in canonical order '
        write(6,*)' ish = ',ish
        write(6,*)' jsh = ',jsh
        call errquit('int_l1eall: canonical error',0)
      endif
c
      if ((3*lstv).gt.lscr) then
        write(6,*)'int_l1eall: need more scratch space '
        write(6,*)'int_l1eall: have :',lscr
        write(6,*)'int_l1eall: need :',(3*lstv)
        call errquit('int_l1eall: scratch space error ',0)
      endif
c
      ibas  = i_basis + BASIS_HANDLE_OFFSET
      jbas  = j_basis + BASIS_HANDLE_OFFSET
c
      iscrS = 1
      call dcopy(lstv,S,1,scr(iscrS),1)
      iscrT = iscrS + lstv
      call dcopy(lstv,T,1,scr(iscrT),1)
      iscrV = iscrT + lstv
      call dcopy(lstv,V,1,scr(iscrV),1)
c
      numstv = 0
      icount = 0
      do 00100 i = ibs_cn2bfr(1,ish,ibas), ibs_cn2bfr(2,ish,ibas)
        do 00200 j = ibs_cn2bfr(1,jsh,jbas), ibs_cn2bfr(2,jsh,jbas)
          nonzero = abs(scr((icount+iscrS))).ge.zerotol
          nonzero = nonzero .and.(abs(scr((icount+iscrT))).ge.zerotol)
          nonzero = nonzero .and.(abs(scr((icount+iscrV))).ge.zerotol)
          if (nonzero) then
            if ((.not.ijbas).or.i.ge.j) then
              numstv = numstv + 1
              S(numstv) = scr((icount+iscrS))
              T(numstv) = scr((icount+iscrT))
              V(numstv) = scr((icount+iscrV))
              ilab(numstv) = i
              ilab(numstv) = j
            endif
          endif
c
          icount = icount + 1
00200   continue
00100 continue
c
      end
