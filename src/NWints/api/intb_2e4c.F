      logical function intb_2e4c(brain, icl, jcl, ketin, kcl, lcl,
     &    num_q, zerotol, canonicalize, ilab, jlab, klab, llab, eri,
     &    l_erilab, nint, lscr, scr)
* $Id: intb_2e4c.F,v 1.3 1995-10-31 09:25:59 d3e129 Exp $
      implicit none
#include "numb_qP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
c::functions
      logical cando_txs
      external cando_txs
c:: passed
      integer brain      ! [input]  basis set handle for bra functions
      integer ketin      ! [input]  basis set handle for ket functions
      integer num_q      ! [input]  number of quartets input
      integer icl(num_q) ! [input]  i-contraction quartet labels
      integer jcl(num_q) ! [input]  j-contraction quartet labels
      integer kcl(num_q) ! [input]  k-contraction quartet labels
      integer lcl(num_q) ! [input]  l-contraction quartet labels
      integer l_erilab   ! [input]  length of eri and integral label arrays
      integer ilab(l_erilab) ! [output] integral labels for "i" gaussian
      integer jlab(l_erilab) ! [output] integral labels for "j" gaussian
      integer klab(l_erilab) ! [output] integral labels for "k" gaussian
      integer llab(l_erilab) ! [output] integral labels for "l" gaussian
      integer nint           ! [output] number of integrals stored in output arrays
      double precision eri(l_erilab) ! [output] integral values
      integer lscr                ! [input] length of scratch array
      double precision scr(lscr)  ! [input] scratch array
      double precision zerotol    ! [input] zero integral threshold
      logical canonicalize        ! [input] Canonicalize integral labels on exit??
c::local
      double precision roff(3)
      integer bra, ket
      integer ip_int, leri_left
      integer numint, numinti, numintj, numintk, numintl, nint_for_q
      integer i_q, i_q_done
      logical status
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      intb_2e4c = .true.
*.... using texas
      status = .true.
      status = (bra .eq. ket) .and. status
      do i_q = num_q_start,(num_q_start+10)
        status = status.and. (cando_txs(brain,icl(i_q),jcl(i_q)))
        status = status.and. (cando_txs(ketin,kcl(i_q),lcl(i_q)))
        if (.not.status) goto 00001
      enddo
*... checking all is costly check first 10
*rak:      do i_q = num_q_start,num_q
*rak:        status = status.and. (cando_txs(brain,icl(i_q),jcl(i_q)))
*rak:        status = status.and. (cando_txs(ketin,kcl(i_q),lcl(i_q)))
*rak:        if (.not.status) goto 00001
*rak:      enddo
      call dfill(3,0.0d00,roff,1)
      call texas_hf2_m(brain,icl,jcl,ketin,kcl,lcl,num_q,
     $     roff,roff,roff,roff,.false.,
     $     eri,l_erilab, ilab, jlab, klab, llab, nint, .true.,
     &    more_ints_txs, scr, lscr)
      if (nint.gt.l_erilab) then
        write(luout,*)' memory trashed.'
        write(luout,*)' number of integrals from texas   :',nint
        write(luout,*)' size of integral and label buffer:',l_erilab
        write(luout,*)' you must increase eri buffer/label size '
        write(luout,*)' '
*        write(luout,*)' Tell Kendall to fix the memory routines!!'
        write(luout,*)' '
        call errquit(' oops ',911)
      endif
      intb_2e4c = more_ints_txs
      i_q_done = num_q  ! all done by texas now with split via more_ints_txs
      return
*.... using nwchem
00001 continue
      intb_2e4c = .true.
      ip_int = 1
      leri_left = l_erilab
      i_q_done = num_q_start - 1
      if (num_q_start.eq.num_q_last) then
        write(luout,*)' l_erilab: too small ',l_erilab
        call errquit(' intb_2e4c: infinite blocking loop ',911)
      endif
      num_q_last = num_q_start
      nint = 0
      do 00100 i_q = num_q_start, num_q
        numinti=ibs_cn2bfr(2,icl(i_q),bra)-ibs_cn2bfr(1,icl(i_q),bra)+1
        numintj=ibs_cn2bfr(2,jcl(i_q),bra)-ibs_cn2bfr(1,jcl(i_q),bra)+1
        numintk=ibs_cn2bfr(2,kcl(i_q),ket)-ibs_cn2bfr(1,kcl(i_q),ket)+1
        numintl=ibs_cn2bfr(2,lcl(i_q),ket)-ibs_cn2bfr(1,lcl(i_q),ket)+1
        numint=numinti*numintj*numintk*numintl
        if (numint.gt.leri_left) then
          goto 00002
        else
          call int_l2e4c(brain, icl(i_q), jcl(i_q),
     &        ketin, kcl(i_q), lcl(i_q),
     &        zerotol, canonicalize, leri_left, eri(ip_int),
     &        nint_for_q,
     &        ilab(ip_int), jlab(ip_int), klab(ip_int), llab(ip_int),
     &        lscr, scr)
          i_q_done = i_q_done + 1
          ip_int = ip_int + nint_for_q
          leri_left = leri_left - nint_for_q
          nint = nint + nint_for_q
        endif
00100 continue
      intb_2e4c = .false.
00002 continue
      num_q_start = i_q_done + 1
      end
