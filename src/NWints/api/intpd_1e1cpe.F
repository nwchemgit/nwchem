C> \ingroup nwint
C> @{
C>
C> \brief Compute 1-electron 1-center Gaussian periodic nuclear attraction
C> integrals derivatives
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
      subroutine intpd_1e1cpe(i_basis,ish,Ri,lscr,scr,lpea,Pea)
*
* $Id$
*
      implicit none
#include "apiP.fh"
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
      integer int_nint_cart
      external int_chk_init
      external int_chk_sh
      external cando_nw_1e
      external cando_nw
      external int_nint_cart
c::passed
      integer i_basis           !< [Input] basis set handle for ish functions
      integer ish               !< [Input] lexical contraction/shell index
      integer lscr              !< [Input] length of the scratch array
      integer lpea              !< [Input] length of potential energy derivative integral array
      double precision Pea(lpea) !< [Output] potential energy derivative integral array
      double precision scr(lscr) !< [Scratch] scratch array
      double precision Ri(3)    !< [Input] translation vector for ish center (fractional coordinates)
c::local
      logical shells_ok
      integer mynint
      integer i_geom, ibas, ucont, mynat
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      double precision xyz_new_i(3) ! new coordinates for ish function center
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('intpd_1e1cpe'))
     &    call errquit('intpd_1e1cpe: int_init was not called' ,0, 
     &            INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      if (.not.shells_ok)
     &    call errquit('intpd_1e1cpe: invalid contraction/shell',0,
     &             BASIS_ERR)
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
      mynint  = int_nint_cart(i_basis,ish, 0,0, 0,0, 0,0)
      mynat = ncenter(i_geom)
      if (lpea*3*mynat.lt.mynint) stop ' internal error '
c
c.. translate ish center coordinates based on Ri
      call intp_txyz(i_cent,i_geom,Ri,xyz_new_i)
c
      if (cando_nw(i_basis,ish,0)) then        
        call hf1d(
     &      xyz_new_i,dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li, i_cent,
     &      xyz_new_i,DCexp,DCcoeff, 1, 1, 0, i_cent,
     &      coords(1,1,i_geom),charge(1,i_geom),
     &      geom_invnucexp(1,i_geom),ncenter(i_geom),
c.............................. doS     doT     doV    canonical
     &      scr,scr,Pea,mynint,.false.,.false.,.true.,.false.,
c........... dryrun
     &      .false.,scr,lscr)
        if (bas_spherical(ibas)) then
          if (Li.eq.-1) i_gen = 1
          call spcart_2cBtran(Pea,scr,lscr,
     &        int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &        int_nbf_x( 0),int_nbf_x( 0), 0,    1,.false.,
     &        (3*ncenter(i_geom)),.false.)
        endif
      else
        call errquit('intpd_1e1cpe: could not do sp or nw integrals',0,
     &            INT_ERR)
      endif
c
      end
C> @}
