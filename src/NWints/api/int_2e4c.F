      subroutine int_2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &       lscr, scr, leri, eri)
c $Id: int_2e4c.F,v 1.6 1994-04-09 15:46:25 d3e129 Exp $
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
c
c::external subroutines used
c errquit
c::functions
      logical cando_nw
      logical cando_sp
      logical int_chk_sh
      logical int_chk_init
      external cando_nw
      external cando_sp
      external int_chk_sh
      external int_chk_init
c:: passed
      integer brain, ketin, lscr, leri
      integer ish,jsh,ksh,lsh
      double precision scr(lscr), eri(leri)
c:: local
      logical shells_ok
      integer bra, ket 
      integer ab_geom, cd_geom, ucont
      integer La, a_prim, a_gen, a_iexp, a_icfp, a_cent
      integer Lb, b_prim, b_gen, b_iexp, b_icfp, b_cent
      integer Lc, c_prim, c_gen, c_iexp, c_icfp, c_cent
      integer Ld, d_prim, d_gen, d_iexp, d_icfp, d_cent
c
      integer WarnP
      data WarnP /0/
      save WarnP
c
c check initialization
c
      if (.not.int_chk_init('int_2e4c'))
     &       call errquit('int_2e4c: int_init was not called' ,0)
c
c     check input shell ranges
c
      shells_ok = int_chk_sh(brain,ish)
      shells_ok = shells_ok .and. int_chk_sh(brain,jsh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,ksh)
      shells_ok = shells_ok .and. int_chk_sh(ketin,lsh)
      if (.not.shells_ok)
     &       call errquit('int_l2e4c: invalid contraction/shell',0)
c
c     define center information required 
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
      ab_geom = ibs_geom(bra)
      cd_geom = ibs_geom(ket)
c
      if (ab_geom.ne.cd_geom.and.WarnP.eq.0) then
        write(6,*)'int_2e4c: WARNING: possible geometry inconsistency'
        write(6,*)'bra geometry handle:',ab_geom
        write(6,*)'ket geometry handle:',cd_geom
        WarnP = 1
      endif
c
      a_cent  = ibs_cn2ce(ish,bra)
      b_cent  = ibs_cn2ce(jsh,bra)
      c_cent  = ibs_cn2ce(ksh,ket)
      d_cent  = ibs_cn2ce(lsh,ket)
c
      if (cando_sp(brain,ish,jsh).and.cando_sp(ketin,ksh,lsh)) then
c
        call genr70(
     &         brain,ish,coords(1,a_cent,bra),
     &               jsh,coords(1,b_cent,bra),
     &         ketin,ksh,coords(1,c_cent,ket),
     &               lsh,coords(1,d_cent,ket),
     &         eri)
c        
      elseif(cando_nw(brain,ish,jsh).and.cando_nw(ketin,ksh,lsh)) then
c
        ucont   = ibs_cn2ucn(ish,bra)
        La      = infbs_cont(CONT_TYPE ,ucont,bra)
        a_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        a_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        a_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        a_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c
        ucont   = ibs_cn2ucn(jsh,bra)
        Lb      = infbs_cont(CONT_TYPE ,ucont,bra)
        b_prim  = infbs_cont(CONT_NPRIM,ucont,bra)
        b_gen   = infbs_cont(CONT_NGEN ,ucont,bra)
        b_iexp  = infbs_cont(CONT_IEXP ,ucont,bra)
        b_icfp  = infbs_cont(CONT_ICFP ,ucont,bra)
c      
        ucont   = ibs_cn2ucn(ksh,ket)
        Lc      = infbs_cont(CONT_TYPE ,ucont,ket)
        c_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        c_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        c_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        c_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c
        ucont   = ibs_cn2ucn(lsh,ket)
        Ld      = infbs_cont(CONT_TYPE ,ucont,ket)
        d_prim  = infbs_cont(CONT_NPRIM,ucont,ket)
        d_gen   = infbs_cont(CONT_NGEN ,ucont,ket)
        d_iexp  = infbs_cont(CONT_IEXP ,ucont,ket)
        d_icfp  = infbs_cont(CONT_ICFP ,ucont,ket)
c      
        call hf2(
     &         coords(1,a_cent,ab_geom), exndcf(a_iexp,bra),
     &         exndcf(a_icfp,bra), a_prim, a_gen, La,
     &         coords(1,b_cent,ab_geom), exndcf(b_iexp,bra),
     &         exndcf(b_icfp,bra), b_prim, b_gen, Lb,
     &         coords(1,c_cent,cd_geom), exndcf(c_iexp,ket),
     &         exndcf(c_icfp,ket), c_prim, c_gen, Lc,
     &         coords(1,d_cent,cd_geom), exndcf(d_iexp,ket),
     &         exndcf(d_icfp,ket), d_prim,d_gen,Ld,
c....................... canAB    canCD    canPQ
     &         eri, leri, .false., .false., .false.,
c............ dryrun
     &         .false., scr, lscr)
c
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
