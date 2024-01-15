c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute 1-electron magnetic hyperfine integrals used with x2c
C>
C> Compute 1-electron magnetic hyperfine integrals as defined
C> by (see e.g. [1], Equations (2.60) and (2.61) and [2], Equation
C> (15) and the text follwing Eq. (15)):
C> \f{eqnarray*}{
C>   I_{x,x'} &=&  \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})
C>   \frac{(X_B-r_1)_x}{|X_B-r_1|^3}
C>   \frac{\partial}{\partial r_1^{x'}}
C>   g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> Above, \f$ B \f$ is the atom for which the magnetic hyperfine
C> integrals are evaluated, \f$ x \f$ and \f$ x' \f$ are the corresponding
C> Cartesian coordinates.
C>
C> [1] M. Dupuis,
C>     <i>"New integral transforms for molecular properties and applications
C>      to a massively parallel GIAO-SCF implementation"</i>,
C>     Comp. Phys. Comm. <b>134</b>, 150-166 (2001), DOI:
C>     <a href="https://doi.org/10.1016/S0010-4655(00)00195-8">
C>     10.1016/S0010-4655(00)00195-8</a>.
C> 
C> [2] J. Autschbach,
C>     <i>"Relativistic Effects on Electron−Nucleus Hyperﬁne Coupling 
C>      Studied with an Exact 2-Component (X2C) Hamiltonian"</i>, DOI:
C>     <a href="https://doi.org/10.1021/acs.jctc.6b01014">
C>     10.1021/acs.jctc.6b01014</a>.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_ops}
c:tex-This routine computes PSO-like integrals, but all i-j combinations
c:tex-of (r_N_i/r_N**3) (d/d r_j)
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_ops(i_basis,ish,j_basis,jsh,lscr,scr,lh11,h11,
     &                       xyzpt,nat)
c:tex-\end{verbatim}
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
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e_prp
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e_prp
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lh11               !< [Input] length of h11 buffer
      double precision h11(lh11) !< [Output] h11 integrals
      integer nat    !< [Input] number of atoms under consideration
      double precision xyzpt(3,nat) !< [Input] coords of atoms under consideration
c:tex-\end{verbatim}
c::local
      integer igeom, jgeom, ibas, jbas, ucont
      integer itype, inp, igen, iexp, icent, icf, iatom
      integer jtype, jnp, jgen, jexp, jcent, jcf, jatom
c
      logical any_spherical, trani, tranj, shells_ok
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
      integer ipts, ncartint ,i,j
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_ops'))
     &       call errquit('int_ops: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_ops: invalid contraction/shell',0,
     &       INT_ERR)
c
c  check if gencont
c
      call int_nogencont_check(i_basis,'int_ops:i_basis')
      call int_nogencont_check(j_basis,'int_ops:j_basis')
c
      ibas = i_basis + basis_handle_offset
      jbas = j_basis + basis_handle_offset
c
      ucont = (sf_ibs_cn2ucn(ish,ibas))
      itype = infbs_cont(CONT_TYPE ,ucont,ibas)
      inp   = infbs_cont(CONT_NPRIM,ucont,ibas)
      igen  = infbs_cont(CONT_NGEN ,ucont,ibas)
      iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      icf   = infbs_cont(CONT_ICFP ,ucont,ibas)
      iatom = (sf_ibs_cn2ce(ish,ibas))
      igeom = ibs_geom(ibas)
c
      ucont = (sf_ibs_cn2ucn(jsh,jbas))
      jtype = infbs_cont(CONT_TYPE ,ucont,jbas)
      jnp   = infbs_cont(CONT_NPRIM,ucont,jbas)
      jgen  = infbs_cont(CONT_NGEN ,ucont,jbas)
      jexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      jcf   = infbs_cont(CONT_ICFP ,ucont,jbas)
      jatom = (sf_ibs_cn2ce(jsh,jbas))
      jgeom = ibs_geom(jbas)
c
      if (igeom.ne.jgeom) then
        write(luout,*)'int_ops: two different geometries for',
     &         ' properties?'
        call errquit('int_ops: geom error ',911, GEOM_ERR)
      endif
c
c     Determine # of cartesian integrals in block
c
      ncartint = int_nbf_x(itype)*int_nbf_x(jtype)
c
      call hnd_ops(
     &       coords(1,iatom,igeom),
     &       dbl_mb(mb_exndcf(iexp,ibas)),
     &       dbl_mb(mb_exndcf(icf,ibas)),
     &       inp,igen,itype,
c
     &       coords(1,jatom,jgeom),
     &       dbl_mb(mb_exndcf(jexp,jbas)),
     &       dbl_mb(mb_exndcf(jcf,jbas)),
     &       jnp,jgen,jtype,
c
     &       ncartint,h11,scr,lscr,
c
     &       xyzpt,nat)        
c
c
c     h11 now has nine blocks, for each point/atom
c      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
      i_nbf_x = int_nbf_x(itype)
      j_nbf_x = int_nbf_x(jtype)
c
c... assume we need to transform both i and j integrals
c
      trani = .true.
      tranj = .true. 
*.. do not tranform i component
      if (.not.bas_spherical(ibas)) trani = .false.
*.. do not tranform j component
      if (.not.bas_spherical(jbas)) tranj = .false.
c        
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4.
c
      if (itype.eq.-1) igen = 1
      if (jtype.eq.-1) jgen = 1
      call spcart_2cBtran(h11,scr,lscr,
     &    j_nbf_x,int_nbf_s(jtype),jtype,jgen,tranj,
     &    i_nbf_x,int_nbf_s(itype),itype,igen,trani,
     &    9,.false.)
c
c     We now have the integrals in array (nsph_ints,9)
c
      return
      end
C> @}
