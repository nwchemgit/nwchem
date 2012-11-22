C> \ingroup nwint
C> @{
C>
C> \brief Computes the electrostatic potential at a point due to a charge
C> distribution
C>
C> Computes the electrostatic potential at a set of points due to the charge
C> distribution given by pairs of basis functions. In other words, compute the
C> integrals
C> \f{eqnarray*}{
C>    P_e(i,\mu,\nu) &=& \int \frac{\chi_\mu(r')\chi_\nu(r')q_i}{|r_i-r'|}dr'
C> \f}
C> where \f$ \chi(r) \f$ is a Gaussian basis function, the \f$ r_i \f$ are
C> the point positions of charges \f$ q_i \f$.
C>
      subroutine int_1er(r,npt,charg,
     &     i_basis,ish,j_basis,jsh,lscr,scr,lpe,Pe)
c $Id$
c
c     computes - for a given set of points r - the integral
c
c     \int dr' chi(mu;r') chi(nu;r')/|r-r'|
c
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "util.fh"
c::external subroutines used
c errquit
c::functions
      logical int_chk_init
      logical int_chk_sh
      logical cando_sp
      logical cando_nw
      external int_chk_init
      external int_chk_sh
      external cando_sp
      external cando_nw
c::passed
      integer i_basis !< [Input] the basis set handle for ish
      integer ish     !< [Input] the i-th shell/contraction
      integer j_basis !< [Input] the basis set handle for jsh
      integer jsh     !< [Input] the j-th shell/contraction
      integer lscr    !< [Input] the length of the scratch space
      integer lpe     !< [Input] the number of charge distributions
      integer npt     !< [Input] the number of points
      double precision scr(lscr)   !< [Scratch] the scratch space
      double precision Pe(npt,lpe) !< [Output] the integrals
      double precision r(3,npt)    !< [Input] the center positions
      double precision charg(npt)  !< [Input] the center charges
c::local
      logical shells_ok
      integer i_geom, j_geom, ibas, jbas, ucont
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
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
      if (.not.int_chk_init('int_1epe'))
     &       call errquit('int_1epe: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1epe: invalid contraction/shell',0,
     &       BASIS_ERR)
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
        write(6,*)'int_1epe: WARNING: possible geometry inconsistency'
        write(6,*)'i_basis geometry handle:',i_geom
        write(6,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
        call hf1_er(
     &         coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &         dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li,
     &         coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &         dbl_mb(mb_exndcf(j_icfp,jbas)), j_prim, j_gen, Lj,
c               coord  charge      no. of centers
     &         r     , charg,             npt ,
c............................. doS     doT     doV    canonical
     &         scr,scr,Pe,lpe,.false.,.false.,.true.,.false.,
c.............. dryrun
     &         .false.,scr,lscr)
c
      end
C> @}
