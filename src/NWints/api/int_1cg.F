      subroutine int_1cg(i_basis,ish,lscr,scr,l1g,G1e)
c $Id: int_1cg.F,v 1.3 1996-01-22 21:46:14 d3e129 Exp $
      implicit none
#include "apiP.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_sp
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_sp
      external cando_nw
c::passed
      integer i_basis, ish, lscr, l1g
      double precision scr(lscr), G1e(l1g)
c::local
      logical shells_ok
      integer i_geom, ibas, ucont
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_1cg'))
     &       call errquit('int_1cg: int_init was not called' ,0)
c
      shells_ok = int_chk_sh(i_basis,ish)
      if (.not.shells_ok)
     &       call errquit('int_1cg: invalid contraction/shell',0)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
c
      if (cando_sp(i_basis,ish,0)) then
        call int_hf1sp(
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li,
     &        coords(1,i_cent,i_geom), DCexp, DCcoeff, 1, 1,
     &        0,coords(1,1,i_geom),charge(1,i_geom),
     &        ncenter(i_geom),
c............................. doS    doT     doV     canonical
     &        G1e,scr,scr,l1g,.true.,.false.,.false.,.false.,
c.............. dryrun
     &        .false.,scr,lscr,'int_1cg')
      elseif (cando_nw(i_basis,ish,0)) then
        call hf1(
     &        coords(1,i_cent,i_geom),
     &        dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li,
     &        coords(1,i_cent,i_geom),DCexp,
     &        DCcoeff, 1, 1, 0,
     &        coords(1,1,i_geom),charge(1,i_geom),ncenter(i_geom),
c............................. doS    doT     doV     canonical
     &         G1e,scr,scr,l1g,.true.,.false.,.false.,.false.,
c.............. dryrun
     &         .false.,scr,lscr)
      else
        call errquit('int_1cg: could not do sp or nw integrals',0)
      endif
c
      end
