      subroutine int_2e3c(brain, ish, ketin, jsh, ksh,
     &       lscr, scr, leri, eri)
c $Id: int_2e3c.F,v 1.4 1996-01-05 01:59:52 d3e129 Exp $
      implicit none
c
c basic api routine to generate a block of 3 center two electron integrals
c eri = <bra_g(ish)|ket_g(jsh).ket_g(ksh)>
c
#include "apiP.fh"
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical int_chk_sh
      logical int_chk_init
      external cando_nw
      external int_chk_sh
      external int_chk_init
#if defined(USE_TEXAS)
      logical cando_txs
      external cando_txs
#endif
c:: passed
      integer brain, ketin, lscr, leri
      integer ish,jsh,ksh
      double precision scr(lscr), eri(leri)
c:: local
      logical shells_ok
      integer bra, ket 
      integer p_geom, cd_geom, ucont
      integer Lp, p_prim, p_gen, p_iexp, p_icfp, p_cent
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent
      logical status_nw
#if defined(USE_TEXAS)
      logical status_txs
      integer num_quart, dummy_lab(2)
      double precision roff(3), q4 
      integer nint
      logical dum_log
#endif
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
c check initialization
c
      if (.not.int_chk_init('int_2e3c'))
     &       call errquit('int_2e3c: int_init was not called' ,0)
c
c     check input shell ranges
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(ketin,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      if (.not.shells_ok)
     &       call errquit('int_2e3c: invalid contraction/shell',0)
c
c set integral code status
c
      status_nw = cando_nw(brain,ish,0).and.cando_nw(ketin,jsh,ksh)
#if defined(USE_TEXAS)
      status_txs = cando_txs(brain,ish,0).and.cando_txs(ketin,jsh,ksh)
#endif
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      p_geom  = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
c
      if (p_geom.ne.cd_geom.and.WarnP.eq.0) then
        write(6,*)'int_2e3c: WARNING: possible geometry inconsistency'
        write(6,*)'bra geometry handle:',p_geom
        write(6,*)'ket geometry handle:',cd_geom
        WarnP = 1
      endif
c
      p_cent  = ibs_cn2ce(ish,bra)
      c_cent  = ibs_cn2ce(jsh,ket)
      d_cent  = ibs_cn2ce(ksh,ket)
c
#if defined(USE_TEXAS)
      if (status_txs) then
        num_quart = 1
        dummy_lab(1) = 0
        dummy_lab(2) = 0
        roff(1) = 0.0d00
        roff(2) = 0.0d00
        roff(3) = 0.0d00
        dum_log = .false.
        q4 = 1.0d00
        nint = 0
        call texas_hf2_m(brain,ish,0,ketin,jsh,ksh,num_quart,
     &      q4,.false.,roff,roff,roff,roff,.false.,
     &      eri, leri, dummy_lab, dummy_lab, dummy_lab, dummy_lab,
     &      nint, .false., dum_log, scr, lscr)

      else if(status_nw) then
#else
      if(status_nw) then
#endif
c
        ucont   = ibs_cn2ucn(ish,bra)
        Lp      = infbs_cont(CONT_TYPE ,ucont,bra)
        p_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        p_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        p_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        p_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = ibs_cn2ucn(jsh,ket)
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c
        ucont   = ibs_cn2ucn(ksh,ket)
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c      
        call hf2(
     &         coords(1,p_cent,p_geom), exndcf(p_iexp,bra),
     &         exndcf(p_icfp,bra), p_prim, p_gen, Lp,
     &         coords(1,p_cent,p_geom), DCexp,
     &         DCcoeff           ,      1,     1, 0,
     &         coords(1,c_cent,cd_geom), exndcf(c_iexp,ket),
     &         exndcf(c_icfp,ket), c_prim, c_gen, Lc,
     &         coords(1,d_cent,cd_geom), exndcf(d_iexp,ket),
     &         exndcf(d_icfp,ket), d_prim,d_gen,Ld,
c......................... canAB    canCD    canPQ
     &         eri, leri, .false., .false., .false.,
c............ dryrun
     &         .false., scr, lscr)
c
      else
        write(6,*)'int_2e3c: could not do nw integrals'
        write(6,*)' brain :',brain
        write(6,*)' ketin :',ketin
        write(6,*)' ish   :',ish
        write(6,*)' jsh   :',jsh
        write(6,*)' ksh   :',ksh
        call errquit('int_2e3c: fatal error ',0)
      endif
      end
