      subroutine int_1e3ov (i_basis, ish, j_basis, jsh, k_basis, ksh,
     &       lscr, scr, lov3, OV3)
c
c routine to calculate the 3 center overlap integrals
c
      implicit none
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
c:: external subroutines used
c..  errquit
c::functions
      logical int_chk_sh, int_chk_init
      external int_chk_sh, int_chk_init
c::passed
      integer i_basis            ! [input] basis set handle for ish
      integer ish                ! [input] i shell/contraction
      integer j_basis            ! [input] basis set handle for jsh
      integer jsh                ! [input] j shell/contraction
      integer k_basis            ! [input] basis set handle for ksh
      integer ksh                ! [input] k shell/contraction
      integer lscr               ! [input] length of scratch array
      double precision scr(lscr) ! [input] scratch array
      integer lov3               ! [input/output] length of 3c overlap integrals array
      double precision OV3(lov3) ! [output] 3c overlap integrals
c::local
      logical shells_ok
      integer ucont
      integer ibas, jbas, kbas
      integer i_geom, j_geom, k_geom
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
      integer Lk, k_prim, k_gen, k_iexp, k_icfp, k_cent
c
      integer WarnP
      save WarnP
      data WarnP /0/
c
c check initialization and shells
c
      if (.not.int_chk_init('int_1e3ov'))
     &       call errquit('int_1e3ov: int_init was not called' ,0)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      shells_ok = shells_ok .and. int_chk_sh(k_basis,ksh)
      if (.not.shells_ok)
     &       call errquit('int_1e3ov: invalid contraction/shell',0)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      kbas = k_basis + BASIS_HANDLE_OFFSET
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
      ucont   = ibs_cn2ucn(ksh,kbas)
      Lk      = infbs_cont(CONT_TYPE ,ucont,kbas)
      k_prim  = infbs_cont(CONT_NPRIM,ucont,kbas)
      k_gen   = infbs_cont(CONT_NGEN ,ucont,kbas)
      k_iexp  = infbs_cont(CONT_IEXP ,ucont,kbas)
      k_icfp  = infbs_cont(CONT_ICFP ,ucont,kbas)
      k_cent  = ibs_cn2ce(ksh,kbas)
      k_geom  = ibs_geom(kbas)
c
      if ((i_geom.ne.j_geom.or.j_geom.ne.k_geom).and.WarnP.eq.0) then
        write(6,*)'int_1e3ov: WARNING: possible geometry inconsistency'
        write(6,*)'i_basis geometry handle:',i_geom
        write(6,*)'j_basis geometry handle:',j_geom
        write(6,*)'k_basis geometry handle:',k_geom
        WarnP = 1
      endif
c
      call hf3OIs(
     &       coords(1,i_cent,i_geom),exndcf(i_iexp,ibas),
     &       exndcf(i_icfp,ibas),i_prim, Li,
     &       coords(1,j_cent,j_geom),exndcf(j_iexp,jbas),
     &       exndcf(j_icfp,jbas),j_prim, Lj,
     &       coords(1,k_cent,k_geom),exndcf(k_iexp,kbas),
     &       exndcf(k_icfp,kbas),k_prim, Lk,
c.....................TriDiag..DryRun
     &       OV3,lov3,.false.,.false.,scr,lscr)
      end
