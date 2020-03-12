C> \ingroup nwint
C> @{
C>
C> \brief Compute 1-electron Gaussian periodic nuclear attraction
C> integrals
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
      subroutine intp_1epe(i_basis,ish,Ri,j_basis,jsh,Rj,
     &    lscr,scr,lpe,Pe)
c $Id$
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
c errquit
c::functions
      logical int_chk_init
      logical int_chk_sh
      logical cando_nw_1e
      logical cando_nw
      external int_chk_init
      external int_chk_sh
      external cando_nw_1e
      external cando_nw
c::passed
      integer i_basis !< [Input] basis set handle for ish functions
      integer j_basis !< [Input] basis set handle for jsh functions
      integer ish     !< [Input] lexical contraction/shell index
      integer jsh     !< [Input] lexical contraction/shell index
      integer lscr    !< [Input] length of the scratch array
      integer lpe     !< [Input] length of potential energy integral array
      double precision Pe(lpe)   !< [Output] nuclear attraction integral array
      double precision scr(lscr) !< [Scratch] scratch array
      double precision Ri(3)  !< [Input] translation vector for ish center (fractional coordinates)
      double precision Rj(3)  !< [Input] translation vector for jsh center (fractional coordinates)
c::local
      logical shells_ok
      integer i_geom, j_geom, ibas, jbas, ucont
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
      double precision xyz_new_i(3)  ! new coordinates for ish function center
      double precision xyz_new_j(3)  ! new coordinates for jsh function center
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
      if (.not.int_chk_init('intp_1epe'))
     &       call errquit('intp_1epe: int_init was not called' ,0,
     &        INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('intp_1epe: invalid contraction/shell',0,
     &         BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
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
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'intp_1epe: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c.. translate ish center coordinates based on Ri
      call intp_txyz(i_cent,i_geom,Ri,xyz_new_i)
c.. translate jsh center coordinates based on Rj
      call intp_txyz(j_cent,j_geom,Rj,xyz_new_j)
c
      if (cando_nw_1e(i_basis,ish,0).and.cando_nw_1e(j_basis,jsh,0))
     &    then
        call int_hf1sp(
     &         xyz_new_i,dbl_mb(mb_exndcf(i_iexp,ibas)),
     &         dbl_mb(mb_exndcf(i_icfp,ibas)),
     &         i_prim, i_gen, Li, i_cent,
     &         xyz_new_j,dbl_mb(mb_exndcf(j_iexp,jbas)),
     &         dbl_mb(mb_exndcf(j_icfp,jbas)),
     &         j_prim, j_gen, Lj, j_cent,
     &         coords(1,1,i_geom),charge(1,i_geom),
     &         geom_invnucexp(1,i_geom),ncenter(i_geom),
c............................. doS     doT     doV    canonical
     &         scr,scr,Pe,lpe,.false.,.false.,.true.,.false.,
c.............. dryrun
     &         .false.,scr,lscr,'intp_1epe')
        if (bas_spherical(ibas).or.bas_spherical(jbas)) then
          if (Li.eq.-1) i_gen = 1
          if (Lj.eq.-1) j_gen = 1
          call spcart_2ctran(Pe,scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &        .false.)          
        endif
      else
        call errquit('intp_1epe: could not do sp or nw integrals',0,
     &           INT_ERR)
      endif
c
      end
C> @}
