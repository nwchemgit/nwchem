      subroutine int_2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       lscr, scr, leri, eri)
c $Id: int_2e4c.F,v 1.22 1996-07-18 23:57:52 d3e129 Exp $
      implicit none
c
c basic api routine to generate a block of two electron integrals
c eri = <bra_g(ish).bra_g(jsh) | ket_g(ksh).ket_g(lsh)>
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical cando_sp
      logical cando_txs
      external cando_nw
      external cando_sp
      external cando_txs
#define PERF_NORMAL
#if defined(PERF_NORMAL)
      logical int_chk_sh
      logical int_chk_init
      external int_chk_sh
      external int_chk_init
#endif
c:: passed
      integer brain, ketin, lscr, leri
      integer ish,jsh,ksh,lsh
      double precision scr(lscr), eri(leri)
c:: local
      integer bra, ket 
      integer ab_geom, cd_geom, ucont
      integer La, a_prim, a_gen, a_iexp, a_icfp, a_cent
      integer Lb, b_prim, b_gen, b_iexp, b_icfp, b_cent
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent
c.txs
      integer num_quart, dummy_lab(2)
      integer txs_i,txs_j,txs_k,txs_l
      double precision roff(3)
      double precision q4
      integer nint
      logical dum_log
      logical status_sp, status_nw, status_txs, status_gen
      integer texas_ang_limit
c
      logical any_spherical
      integer a_nbf,   b_nbf,   c_nbf,   d_nbf
      integer a_nbf_s, b_nbf_s, c_nbf_s, d_nbf_s
c
#if defined(PERF_NORMAL)
      logical shells_ok
      logical inline_chk_sh
c
*debug:      integer ilo, jlo, klo, llo, ihi, jhi, khi, lhi
*debug:      integer aa,bb,cc,dd,abcd
      integer WarnP
      save WarnP
      data WarnP /0/
c... statement function for int_chk_sh
      inline_chk_sh(bra,ish) =
     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(bra)))
#endif
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
#if defined(PERF_NORMAL)
c
c check initialization
c
      if (.not.int_chk_init('int_2e4c'))
     &       call errquit('int_2e4c: int_init was not called' ,0)
c
c     define center information required 
c
      if (.not. bas_check_handle(brain,'int_2e4c:bra')) call errquit
     $     ('int_2e4c: invalid bra basis handle', brain)
      if (.not. bas_check_handle(ketin,'int_2e4c:ket')) call errquit
     $     ('int_2e4c: invalid ket basis handle', ketin)
#endif
*      call int_nospherical_check(brain,'int_2e4c')
*      call int_nospherical_check(ketin,'int_2e4c')
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
#if defined(PERF_NORMAL)
c
c     check input shell ranges
c
      shells_ok = inline_chk_sh(bra,ish)
      shells_ok = shells_ok .and. inline_chk_sh(bra,jsh)
      shells_ok = shells_ok .and. inline_chk_sh(ket,ksh)
      shells_ok = shells_ok .and. inline_chk_sh(ket,lsh)
      if (.not.shells_ok)
     &       call errquit('int_2e4c: invalid contraction/shell',0)
c
      if (ab_geom.ne.cd_geom.and.WarnP.eq.0) then
        write(6,*)'int_2e4c: WARNING: possible geometry inconsistency'
        write(6,*)'bra geometry handle:',ab_geom
        write(6,*)'ket geometry handle:',cd_geom
        WarnP = 1
      endif
c
#endif
      a_cent  = (sf_ibs_cn2ce(ish,bra))
      b_cent  = (sf_ibs_cn2ce(jsh,bra))
      c_cent  = (sf_ibs_cn2ce(ksh,ket))
      d_cent  = (sf_ibs_cn2ce(lsh,ket))
c
c
      any_spherical = bas_spherical(bra).or.bas_spherical(ket)
c
c... stat sp
      status_sp = cando_sp(brain,ish,jsh).and.cando_sp(ketin,ksh,lsh)

c... stat texas
      texas_ang_limit = 2
      ucont   = (sf_ibs_cn2ucn(ish,bra))
      La      = infbs_cont(CONT_TYPE ,ucont,bra)
      a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
      ucont   = (sf_ibs_cn2ucn(jsh,bra))
      Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
      b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
      ucont   = (sf_ibs_cn2ucn(ksh,ket))
      Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
      c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
      ucont   = (sf_ibs_cn2ucn(lsh,ket))
      Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
      d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
      status_txs = cando_txs(brain,ish,jsh).and.cando_txs(ketin,ksh,lsh)
      status_txs = status_txs .and. brain.eq.ketin
      status_gen = max(a_gen,b_gen,c_gen,d_gen).gt.1
      if (.not.status_gen) then
        status_txs = status_txs .and.
     &        ( (La.gt.texas_ang_limit).or.
     &        (Lb.gt.texas_ang_limit).or.
     &        (Lc.gt.texas_ang_limit).or.
     &        (Ld.gt.texas_ang_limit)  )
      endif
      status_txs = status_txs .and. .not.any_spherical
c... stat nw
      status_nw  = cando_nw(brain,ish,jsh).and.cando_nw(ketin,ksh,lsh)
      if (status_sp) then
c
        call genr70(
     &         brain,ish,coords(1,a_cent,ab_geom),
     &               jsh,coords(1,b_cent,ab_geom),
     &         ketin,ksh,coords(1,c_cent,cd_geom),
     &               lsh,coords(1,d_cent,cd_geom),
     &         eri)
c        
      elseif(status_txs)  then
        num_quart = 1
        dummy_lab(1) = 0
        dummy_lab(2) = 0
        roff(1) = 0.0d00
        roff(2) = 0.0d00
        roff(3) = 0.0d00
        dum_log = .false.
*...check symmetry here
        q4 = 1.0d00
        txs_i = ish
        txs_j = jsh
        txs_k = ksh
        txs_l = lsh
        call texas_hf2_m(
     &        brain,txs_i,txs_j,
     &        ketin,txs_k,txs_l,num_quart,
     &        q4,.false.,
c............................... use roff set false
     &        roff,roff,roff,roff,.false.,
     &        eri, leri, dummy_lab, dummy_lab, dummy_lab, dummy_lab,
c............... gen labs .. more_integrals
     &        nint, .false., dum_log, scr, lscr)

      elseif(status_nw) then
c
        ucont   = (sf_ibs_cn2ucn(ish,bra))
        La      = infbs_cont(CONT_TYPE ,ucont,bra)
        a_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        a_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        a_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = (sf_ibs_cn2ucn(jsh,bra))
        Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
        b_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        b_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        b_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c      
        ucont   = (sf_ibs_cn2ucn(ksh,ket))
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c
        ucont   = (sf_ibs_cn2ucn(lsh,ket))
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c      
        call hf2(
     &         coords(1,a_cent,ab_geom), dbl_mb(mb_exndcf(a_iexp,bra)),
     &         dbl_mb(mb_exndcf(a_icfp,bra)), a_prim, a_gen, La,
     &         coords(1,b_cent,ab_geom), dbl_mb(mb_exndcf(b_iexp,bra)),
     &         dbl_mb(mb_exndcf(b_icfp,bra)), b_prim, b_gen, Lb,
     &         coords(1,c_cent,cd_geom), dbl_mb(mb_exndcf(c_iexp,ket)),
     &         dbl_mb(mb_exndcf(c_icfp,ket)), c_prim, c_gen, Lc,
     &         coords(1,d_cent,cd_geom), dbl_mb(mb_exndcf(d_iexp,ket)),
     &         dbl_mb(mb_exndcf(d_icfp,ket)), d_prim,d_gen,Ld,
c....................... canAB    canCD    canPQ
     &         eri, leri, .false., .false., .false.,
c............ dryrun
     &         .false., scr, lscr)
c
c
c eri has cartesian block of integrals (llo:lhi,klo:khi,jlo:jhi,ilo:ihi)
c
        if (any_spherical)then
*debug:          write(6,*)' for shells:',ish,jsh,ksh,lsh
*debug:          if (.not.bas_cn2bfr(brain,ish,ilo,ihi))
*debug:     &        stop 'cn2bfr error i'
*debug:          if (.not.bas_cn2bfr(brain,jsh,jlo,jhi))
*debug:     &        stop 'cn2bfr error j'
*debug:          if (.not.bas_cn2bfr(ketin,ksh,klo,khi))
*debug:     &        stop 'cn2bfr error k'
*debug:          if (.not.bas_cn2bfr(ketin,lsh,llo,lhi))
*debug:     &        stop 'cn2bfr error l'
*debug:          write(6,*)' nbf low   :',ilo,jlo,klo,llo
*debug:          write(6,*)' nbf hi    :',ihi,jhi,khi,lhi
          a_nbf = int_nbf_x(La)
          b_nbf = int_nbf_x(Lb)
          c_nbf = int_nbf_x(Lc)
          d_nbf = int_nbf_x(Ld)
*debug:          write(6,*)'bra/ket-cartesian integrals'
*debug:          abcd = 0
*debug:          do aa = 1,a_nbf
*debug:            do bb = 1,b_nbf
*debug:              do cc = 1,c_nbf
*debug:                do dd = 1,d_nbf
*debug:                  abcd = abcd + 1
*debug:                  write(6,*)aa,bb,cc,dd,eri(abcd)
*debug:                  write(79,*)aa,bb,cc,dd,eri(abcd)
*debug:                enddo
*debug:              enddo
*debug:            enddo
*debug:          enddo
          if(bas_spherical(bra)) then
            a_nbf_s = int_nbf_s(La)
            b_nbf_s = int_nbf_s(Lb)
            call spcart_bra2etran(eri,scr,
     &            b_nbf,a_nbf,b_nbf_s,a_nbf_s,
     &            Lb, La, b_gen, a_gen,
     &            (c_nbf*d_nbf),.false.)
            a_nbf = a_nbf_s
            b_nbf = b_nbf_s
          endif
*debug:          write(6,*)'bra-shperical/ket-cartesian integrals'
*debug:          abcd = 0
*debug:          do aa = 1,a_nbf
*debug:            do bb = 1,b_nbf
*debug:              do cc = 1,c_nbf
*debug:                do dd = 1,d_nbf
*debug:                  abcd = abcd + 1
*debug:                  write(6,*)aa,bb,cc,dd,eri(abcd)
*debug:                  write(80,*)aa,bb,cc,dd,eri(abcd)
*debug:                enddo
*debug:              enddo
*debug:            enddo
*debug:          enddo
          if(bas_spherical(ket)) then
            c_nbf_s = int_nbf_s(Lc)
            d_nbf_s = int_nbf_s(Ld)
            call spcart_ket2etran(eri,scr,
     &            d_nbf,c_nbf,d_nbf_s,c_nbf_s,
     &            Ld, Lc, d_gen, c_gen,
     &            (a_nbf*b_nbf),.false.)
            c_nbf = c_nbf_s
            d_nbf = d_nbf_s
          endif
        endif
*debug:        write(6,*)'bra-shperical/ket-spherical integrals'
*debug:        abcd = 0
*debug:        do aa = 1,a_nbf
*debug:          do bb = 1,b_nbf
*debug:            do cc = 1,c_nbf
*debug:              do dd = 1,d_nbf
*debug:                abcd = abcd + 1
*debug:                write(6,*)aa,bb,cc,dd,eri(abcd)
*debug:                write(80,*)aa,bb,cc,dd,eri(abcd)
*debug:              enddo
*debug:            enddo
*debug:          enddo
*debug:        enddo
      else
        write(6,*)'int_2e4c: could not do sp or nw integrals'
        write(6,*)' brain :',brain
        write(6,*)' ketin :',ketin
        write(6,*)' ish   :',ish
        write(6,*)' jsh   :',jsh
        write(6,*)' ksh   :',ksh
        write(6,*)' lsh   :',lsh
        call errquit('int_2e4c: fatal error ',0)
      endif
      end
