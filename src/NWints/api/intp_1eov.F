      subroutine intp_1eov(i_basis,ish,j_basis,jsh,R,lscr,scr,lov,Ov)
c $Id: intp_1eov.F,v 1.2 1994-12-10 15:51:19 d3e129 Exp $
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
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
      integer i_basis, j_basis, ish, jsh, lscr, lov
      double precision scr(lscr), Ov(lov)
      double precision R(3)  ! translational vector for periodic
c::local
      logical shells_ok
      integer i_geom, j_geom, ibas, jbas, ucont
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
      double precision xyz_frac(3) ! fractional coordinates for jsh function center
      double precision xyz_new(3)  ! new coordinates for jsh function center
c::debug
      integer jjj
c
      logical inline_chk_sh
      integer WarnP
      data WarnP /0/
      save WarnP
c... statement function for int_chk_sh
      inline_chk_sh(ibas,ish) =
     $     ((ish.gt.0) .and. (ish.le.ncont_tot_gb(ibas)))
c
c check initialization and shells
c
      if (.not.int_chk_init('intp_1eov'))
     &       call errquit('intp_1eov: int_init was not called' ,0)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      shells_ok = inline_chk_sh(ibas,ish)
      shells_ok = shells_ok .and. inline_chk_sh(jbas,jsh)
      if (.not.shells_ok)
     &       call errquit('intp_1eov: invalid contraction/shell',0)
c
      ucont   = ibs_cn2ucn(ish,ibas)
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = ibs_cn2ce(ish,ibas)
      i_geom  = ibs_geom(ibas)
c
      ucont   = ibs_cn2ucn(jsh,jbas)
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = ibs_cn2ce(jsh,jbas)
      j_geom  = ibs_geom(jbas)
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(6,*)'int_1eall: WARNING: possible geometry inconsistency'
        write(6,*)'i_basis geometry handle:',i_geom
        write(6,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
c.. zero local arrays
      call dfill(3,0.0d00,xyz_new,1)
      call dfill(3,0.0d00,xyz_frac,1)
c
c... compute j center fractional coordinates
c....    xyz_frac = amatrix_inv*coords(jsh)
c
      call dgemv('n',3,3,1.0d00,amatrix_inv(1,1,j_geom),3,
     &    coords(1,j_cent,j_geom),1,
     &    0.0d00,xyz_frac,1)
c
c... compute j translated fractional coordinates
c.........  xyz_frac = xyz_frac + R
c
      do 00100 jjj=1,3
        xyz_frac(jjj) = xyz_frac(jjj) + R(jjj)
00100 continue
c
c... computer j translated cart. coordinates.
c........ xyz_new = amatrix*xyz_frac
c
      call dgemv('n',3,3,1.0d00,amatrix(1,1,j_geom),3,
     &    xyz_frac,1,
     &    0.0d00,xyz_new,1)
c
      if (cando_sp(i_basis,ish,0).and.cando_sp(j_basis,jsh,0)) then
        call int_hf1sp(
     &         coords(1,i_cent,i_geom),exndcf(i_iexp,ibas),
     &         exndcf(i_icfp,ibas), i_prim, i_gen, Li,
     &         xyz_new,exndcf(j_iexp,jbas),
     &         exndcf(j_icfp,jbas), j_prim, j_gen, Lj,
     &         coords(1,1,i_geom),charge(1,i_geom),ncenter(i_geom),
c............................. doS    doT     doV     canonical
     &         Ov,scr,scr,lov,.true.,.false.,.false.,.false.,
c.............. dryrun
     &         .false.,scr,lscr,'intp_1eov')
      elseif (cando_nw(i_basis,ish,0).and.cando_nw(j_basis,jsh,0)) then
        call hf1(
     &         coords(1,i_cent,i_geom),exndcf(i_iexp,ibas),
     &         exndcf(i_icfp,ibas), i_prim, i_gen, Li,
     &         xyz_new,exndcf(j_iexp,jbas),
     &         exndcf(j_icfp,jbas), j_prim, j_gen, Lj,
     &         coords(1,1,i_geom),charge(1,i_geom),ncenter(i_geom),
c............................. doS    doT     doV     canonical
     &         Ov,scr,scr,lov,.true.,.false.,.false.,.false.,
c.............. dryrun
     &         .false.,scr,lscr)
      else
        call errquit('intp_1eov: could not do sp or nw integrals',0)
      endif
c
      end
