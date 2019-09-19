C> \ingroup nwint
C> @{
C>
C> \brief Compute 1-electron Gaussian periodic overlap
C> integral derivatives
C>
C> See [1] for details.
C>
C> [1] JE Jaffe, AC Hess,
C>     <i>"Gaussian basis density functional theory for systems
C>     periodic in two or three dimensions: Energy and forces"</i>,
C>    J.Chem.Phys. <b>105</b>, 10983-10998 (1996), DOI:
C>    <a href="https://doi.org/10.1063/1.472866">
C>    10.1063/1.472866</a>
C>
      subroutine intpd_1eov(i_basis,ish,j_basis,jsh,R,lscr,scr,
     &    lOva,Ova,idatom)
*
* $Id$
*
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_nw_1e
      external cando_nw
      integer int_nint_cart
      external int_nint_cart
c::passed
      integer i_basis           !< [Input] basis set handle for ish functions
      integer j_basis           !< [Input] basis set handle for jsh functions
      integer ish               !< [Input] lexical contraction/shell index
      integer jsh               !< [Input] lexical contraction/shell index
      integer lscr              !< [Input] length of the scratch array
      integer lOva              !< [Input] length of the overlap integral array
      double precision Ova(lOva) !< [Output] overlap integral array
      double precision scr(lscr) !< [Scratch] scratch array
      double precision R(3)     !< [Input] translational vector fractional coordinates
      integer idatom(*)         !< [Output] array identifying centers for derivatives
c                       ! e.g., the first nint*3  derivatives go to center idatom(1)
c                       !       the second nint*3 derivatives go to center idatom(2)
c::local
      logical shells_ok
      integer i_geom, j_geom, ibas, jbas, ucont, mynint
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
      double precision xyz_new(3) ! new coordinates for jsh function center
*rak:      integer jjj
c
      logical inline_chk_sh
      integer WarnP
      save WarnP
      data WarnP /0/
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c... statement function for int_chk_sh
      inline_chk_sh(ibas,ish) =
     $    ((ish.gt.0) .and. (ish.le.ncont_tot_gb(ibas)))
c
c check initialization and shells
c
      if (.not.int_chk_init('intpd_1eov'))
     &    call errquit('intpd_1eov: int_init was not called' ,0,
     &            INT_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
c
      shells_ok = inline_chk_sh(ibas,ish)
      shells_ok = shells_ok .and. inline_chk_sh(jbas,jsh)
      if (.not.shells_ok)
     &    call errquit('intpd_1eov: invalid contraction/shell',0,
     &            BASIS_ERR)
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
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
c
      mynint = int_nint_cart(i_basis,ish,j_basis,jsh,0,0,0,0)
      if (i_cent.eq.j_cent) then
*        write(luout,*)' automatic zero '
        call ifill(2,0,idatom,1)
        call dcopy((mynint*3*2),0.0d00,0,Ova,1)
        return
      endif
      idatom(1) = i_cent
      idatom(2) = j_cent
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'intpd_1eov: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
c.. translate coordinates based on R
      call intp_txyz(j_cent,j_geom,R,xyz_new)
c
      if (cando_nw(i_basis,ish,0).and.cando_nw(j_basis,jsh,0)) then
        call hf1d(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li, i_cent,
     &      xyz_new,dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)), j_prim, j_gen, Lj, j_cent,
     &      coords(1,1,i_geom),charge(1,i_geom),
     &      geom_invnucexp(1,i_geom),ncenter(i_geom),
c.............................. doS    doT     doV     canonical
     &      Ova,scr,scr,mynint,.true.,.false.,.false.,.false.,
c........... dryrun
     &      .false.,scr,lscr)
        if (bas_spherical(ibas).or.bas_spherical(jbas)) then
          if (Li.eq.-1) i_gen = 1
          if (Lj.eq.-1) j_gen = 1
          call spcart_2cBtran(Ova,scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        (3*2),.false.)
        endif
      else
        call errquit('intpd_1eov: could not do sp or nw integrals',0,
     &                INT_ERR)
      endif
c
      end
C> @}
