      logical function intb_2e4c(brain, icl, jcl, ketin, kcl, lcl,
     $     num_q, q4, use_q4, zerotol, canonicalize,
     $     ilab, jlab, klab, llab, eri,
     $     l_erilab, nint, lscr, scr)
*     $Id: intb_2e4c.F,v 1.15 1997-02-20 09:08:52 d3e129 Exp $
      implicit none
#include "numb_qP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c::   functions
      logical intb_nw_2e4c
      external intb_nw_2e4c
c::   passed
      integer brain             ! [input]  basis set handle for bra functions
      integer ketin             ! [input]  basis set handle for ket functions
      integer num_q             ! [input]  number of quartets input
      integer icl(num_q)        ! [input]  i-contraction quartet labels
      integer jcl(num_q)        ! [input]  j-contraction quartet labels
      integer kcl(num_q)        ! [input]  k-contraction quartet labels
      integer lcl(num_q)        ! [input]  l-contraction quartet labels
      logical use_q4            ! [input]  logical stating whether to use q4
      double precision q4(num_q) ! [input] symmetry prefactors for integrals
*     
      integer l_erilab          ! [input]  length of eri and integral label arrays
      integer ilab(l_erilab)    ! [output] integral labels for ``i'' gaussian
      integer jlab(l_erilab)    ! [output] integral labels for ``j'' gaussian
      integer klab(l_erilab)    ! [output] integral labels for ``k'' gaussian
      integer llab(l_erilab)    ! [output] integral labels for ``l'' gaussian
      integer nint              ! [output] number of integrals stored in output arrays
      double precision eri(l_erilab) ! [output] integral values
      integer lscr              ! [input] length of scratch array
      double precision scr(lscr) ! [input] scratch array
      double precision zerotol  ! [input] zero integral threshold
      logical canonicalize      ! [input] Canonicalize integral labels on exit??
c::   local
      double precision roff(3)
c     
c     Cannot do texas.
c     
      if (brain.ne.ketin .or. (.not.txs_b2e4c)) then
         intb_2e4c = intb_nw_2e4c(
     $        brain, icl, jcl, ketin, kcl, lcl,
     $        num_q, q4, use_q4, zerotol, canonicalize,
     $        ilab, jlab, klab, llab, eri,
     $        l_erilab, nint, lscr, scr)
         return
      endif
c     
c     If faster to use SP than texas
c     
      if (nsp .gt. 0) then
         intb_2e4c = intb_nw_2e4c(
     $        brain, icl, jcl, ketin, kcl, lcl,
     $        nsp, q4, use_q4, zerotol, canonicalize,
     $        ilab, jlab, klab, llab, eri,
     $        l_erilab, nint, lscr, scr)
         if ((nint.eq.0) .or. (.not. intb_2e4c)) then
c     
c     Must skip past those already done with SP
c     
            num_q_start = nsp + 1
            nsp = 0
            if (num_q_start .gt. num_q) then
               intb_2e4c = .false.
               num_q_start = -99999999
            else
               intb_2e4c = .true.
            endif
            return
         else
            return
         endif
      endif
c     
      intb_2e4c = .true.
c     
      call dfill(3,0.0d00,roff,1)
      call texas_hf2_m(
     $     brain, icl(num_q_start), jcl(num_q_start),
     $     ketin, kcl(num_q_start), lcl(num_q_start), 
     $     num_q-num_q_start+1,
     $     q4(num_q_start), use_q4,
     $     roff,roff,roff,roff,.false.,
     $     eri,l_erilab, ilab, jlab, klab, llab, nint, .true.,
     $     more_ints_txs, scr, lscr, zerotol)
      if (nint.gt.l_erilab) then
         write(luout,*)' number of integrals from texas   :',nint
         write(luout,*)' size of integral and label buffer:',l_erilab
         write(luout,*)' you must increase eri buffer/label size '
         call errquit(' integral buffer overflow ',911)
      endif
      intb_2e4c = more_ints_txs
c     
      end
      logical function intb_nw_2e4c(brain, icl, jcl, ketin, kcl, lcl,
     $     num_q, q4, use_q4, zerotol, canonicalize,
     $     ilab, jlab, klab, llab, eri,
     $     l_erilab, nint, lscr, scr)
*     $Id: intb_2e4c.F,v 1.15 1997-02-20 09:08:52 d3e129 Exp $
      implicit none
#include "numb_qP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::   passed
      integer brain             ! [input]  basis set handle for bra functions
      integer ketin             ! [input]  basis set handle for ket functions
      integer num_q             ! [input]  number of quartets input
      integer icl(num_q)        ! [input]  i-contraction quartet labels
      integer jcl(num_q)        ! [input]  j-contraction quartet labels
      integer kcl(num_q)        ! [input]  k-contraction quartet labels
      integer lcl(num_q)        ! [input]  l-contraction quartet labels
      logical use_q4            ! [input]  logical stating whether to use q4
      double precision q4(num_q) ! [input] symmetry prefactors for integrals
*     
      integer l_erilab          ! [input]  length of eri and integral label arrays
      integer ilab(l_erilab)    ! [output] integral labels for ``i'' gaussian
      integer jlab(l_erilab)    ! [output] integral labels for ``j'' gaussian
      integer klab(l_erilab)    ! [output] integral labels for ``k'' gaussian
      integer llab(l_erilab)    ! [output] integral labels for ``l'' gaussian
      integer nint              ! [output] number of integrals stored in output arrays
      double precision eri(l_erilab) ! [output] integral values
      integer lscr              ! [input] length of scratch array
      double precision scr(lscr) ! [input] scratch array
      double precision zerotol  ! [input] zero integral threshold
      logical canonicalize      ! [input] Canonicalize integral labels on exit??
c::   local
      integer bra, ket
      integer ip_int, leri_left
      integer numint, numinti, numintj, numintk, numintl, nint_for_q
      integer i_q, i_q_done
c     
#include "bas_ibs_sfn.fh"
c     
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c     
      intb_nw_2e4c = .true.
c     
      ip_int = 1
      leri_left = l_erilab
      i_q_done = num_q_start - 1
#if defined(DEBUG_TXS)
      write(6,*)' num_q_start ',num_q_start
      write(6,*)' num_q_last  ',num_q_last
#endif
      if (num_q_start.eq.num_q_last) then
         write(luout,*)' l_erilab: too small ',l_erilab
         call errquit(' intb_nw_2e4c: infinite blocking loop ',911)
      endif
      num_q_last = num_q_start
      nint = 0
*     this will have to change for spherical code
      do 00100 i_q = num_q_start, num_q
         numinti=(sf_ibs_cn2bfr(2,icl(i_q),bra))-
     $        (sf_ibs_cn2bfr(1,icl(i_q),bra))+1
         numintj=(sf_ibs_cn2bfr(2,jcl(i_q),bra))-
     $        (sf_ibs_cn2bfr(1,jcl(i_q),bra))+1
         numintk=(sf_ibs_cn2bfr(2,kcl(i_q),ket))-
     $        (sf_ibs_cn2bfr(1,kcl(i_q),ket))+1
         numintl=(sf_ibs_cn2bfr(2,lcl(i_q),ket))-
     $        (sf_ibs_cn2bfr(1,lcl(i_q),ket))+1
         numint=numinti*numintj*numintk*numintl
         if (numint.gt.leri_left) then
#if defined(DEBUG_TXS)
           write(6,*)' shells ',icl(i_q),jcl(i_q),kcl(i_q),lcl(i_q)
           write(6,*)' numinti  ',numinti
           write(6,*)' numintj  ',numintj
           write(6,*)' numintk  ',numintk
           write(6,*)' numintl  ',numintl
           write(6,*)' numint   ',numint
           write(6,*)' leri_left',leri_left
#endif
            goto 00002
         else
            call int_l2e4c(brain, icl(i_q), jcl(i_q),
     $           ketin, kcl(i_q), lcl(i_q),
     $           zerotol, canonicalize, leri_left, eri(ip_int),
     $           nint_for_q,
     $           ilab(ip_int), jlab(ip_int), klab(ip_int), llab(ip_int),
     $           lscr, scr)
            if (use_q4 .and. nint_for_q.gt.0)
     $           call dscal(nint_for_q,q4(i_q),eri(ip_int),1)
            i_q_done = i_q_done + 1
#if defined(DEBUG_TXS)
            write(6,*)' numinti  ',numinti
            write(6,*)' numintj  ',numintj
            write(6,*)' numintk  ',numintk
            write(6,*)' numintl  ',numintl
            write(6,*)' i_q_done = ',i_q_done,
#endif
     &          '   shells ',icl(i_q),jcl(i_q),kcl(i_q),lcl(i_q),
     &          '   numint ',numint,' nint_for_q ',nint_for_q


            ip_int = ip_int + nint_for_q
            leri_left = leri_left - nint_for_q
            nint = nint + nint_for_q
         endif
00100 continue
      intb_nw_2e4c = .false.
00002 continue
      num_q_start = i_q_done + 1
#if defined(DEBUG_TXS)
      write(6,*)' num_q_start on next pass is ',num_q_start
#endif
c     
      end
