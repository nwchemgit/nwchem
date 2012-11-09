c $Id$
*
C> \ingroup nwint
C> @{
C> \brief Calculates 3-center overlap integrals
C>
C> This routine computes the 3 center overlap integral:
C> \f{eqnarray*}{
C> ({\mu}{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})g_{\lambda}(X_{\lambda},r_{1})dr_{1}
C> \f}
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1e3ov}
c:tex-This routine computes the 3 center overlap integral:
c:tex-\begin{eqnarray*}
c:tex-({\mu}{\nu}{\lambda}) = \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})g_{\lambda}(X_{\lambda},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-{\it Syntax:} 
c:tex-\begin{verbatim}
      subroutine int_1e3ov (i_basis, ish, j_basis, jsh, k_basis, ksh,
     &       lscr, scr, lov3, OV3)
c
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "util.fh"
c:: external subroutines used
c..  errquit
c::functions
      logical int_chk_sh, int_chk_init
      external int_chk_sh, int_chk_init
c::passed
c:tex-\begin{verbatim} 
      integer i_basis            !< [Input] basis set handle for ish
      integer ish                !< [Input] i shell/contraction
      integer j_basis            !< [Input] basis set handle for jsh
      integer jsh                !< [Input] j shell/contraction
      integer k_basis            !< [Input] basis set handle for ksh
      integer ksh                !< [Input] k shell/contraction
      integer lscr               !< [Input] length of scratch v
      double precision scr(lscr) !< [Scratch] scratch array
      integer lov3               !< [Input] length of 3c overlap buffer
      double precision OV3(lov3) !< [Output] 3c overlap integrals
c:tex-\end{verbatim} 
c::local
      logical shells_ok
      integer ucont
      integer ibas, jbas, kbas
      integer i_geom, j_geom, k_geom
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
      integer Lk, k_prim, k_gen, k_iexp, k_icfp, k_cent
      integer nint_ck
      integer int_nint_cart
      external int_nint_cart
*     integer nintx, ninta, nintb
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
      if (.not.int_chk_init('int_1e3ov'))
     &       call errquit('int_1e3ov: int_init was not called' ,0,
     &       INT_ERR)
      nint_ck = int_nint_cart(i_basis,ish,j_basis,jsh,k_basis,ksh,
     ,     k_basis,0)
      if (nint_ck.gt.lov3) then
        write(luout,*) 'nint = ',nint_ck
        write(luout,*) 'leri     = ',lov3
        call errquit('int_1e3ov: nint>leri error',911, INT_ERR)
      endiF
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      shells_ok = shells_ok .and. int_chk_sh(k_basis,ksh)
      if (.not.shells_ok)
     &    call errquit('int_1e3ov: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      call int_nogencont_check(i_basis,'int_1e3ov:i_basis')
      call int_nogencont_check(j_basis,'int_1e3ov:j_basis')
      call int_nogencont_check(k_basis,'int_1e3ov:k_basis')
      call int_nospshell_check(i_basis,'int_1e3ov:i_basis')
      call int_nospshell_check(j_basis,'int_1e3ov:j_basis')
      call int_nospshell_check(k_basis,'int_1e3ov:k_basis')
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      kbas = k_basis + BASIS_HANDLE_OFFSET
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
      ucont   = (sf_ibs_cn2ucn(ksh,kbas))
      Lk      = infbs_cont(CONT_TYPE ,ucont,kbas)
      k_prim  = infbs_cont(CONT_NPRIM,ucont,kbas)
      k_gen   = infbs_cont(CONT_NGEN ,ucont,kbas)
      k_iexp  = infbs_cont(CONT_IEXP ,ucont,kbas)
      k_icfp  = infbs_cont(CONT_ICFP ,ucont,kbas)
      k_cent  = (sf_ibs_cn2ce(ksh,kbas))
      k_geom  = ibs_geom(kbas)
c
      if ((i_geom.ne.j_geom.or.j_geom.ne.k_geom).and.WarnP.eq.0) then
        write(luout,*)
     &      'int_1e3ov: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        write(luout,*)'k_basis geometry handle:',k_geom
        WarnP = 1
      endif
c
      call hf3OIs(
     &      coords(1,i_cent,i_geom),
     &      dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      i_prim, Li,
     &      coords(1,j_cent,j_geom),
     &      dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      j_prim, Lj,
     &      coords(1,k_cent,k_geom),
     &      dbl_mb(mb_exndcf(k_iexp,kbas)),
     &      dbl_mb(mb_exndcf(k_icfp,kbas)),
     &      k_prim, Lk,
c.....................TriDiag..DryRun
     &       OV3,lov3,.false.,.false.,scr,lscr)
c
*  OV3 has the cartesian block (klo:khi,jlo:jhi,ilo:ihi)
*
      if (.not.(bas_spherical(ibas).or.
     &          bas_spherical(jbas).or.
     &          bas_spherical(kbas))) return  ! no spherical

*... when general contraction reset sp function ngen to 1
*... since int_nbf_* arrays are set to the appropriate size.
*
      if (Li.eq.-1)  i_gen = 1
      if (Lj.eq.-1)  j_gen = 1
      if (Lk.eq.-1)  k_gen = 1
      call spcart_3ctran(OV3,scr,lscr,
     &    int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &    int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &    int_nbf_x(Lk),int_nbf_s(Lk),Lk,k_gen,bas_spherical(kbas),
     &    .false.)
c
      end
C> @}
