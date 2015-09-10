c $Id$
*
C> \ingroup nwint
C> @{
C> \brief Compute any 1-electron integrals
C>
C> This is an internal routine that most of the external 1-electron
C> routines call.  This is the actual workhorse routine.
C> This routine computes the 1 electron integrals \f$S\f$, \f$T\f$, and \f$V\f$:
C> \f{eqnarray*}{
C> S & = & ({\mu}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> T & = & ({\mu}|-\frac{1}{2}\nabla^{2}|{\nu}) \\\\
C>   & = & -\frac{1}{2}\int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\nabla^{2}(r_{1})\\
C>         g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
C> {-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1} 
C> \f}
C> 
C> If an ECP is defined then the ECP integral contributions are summed
C> directly into the V integrals.
C> 
C> If a relativistic basis is defined then the 1-electron integrals
C> for the case where both shells are relativistic are modified to
C> \f{eqnarray*}{
C> S & = & ({\mu^L}|{\nu^L})
C>       - ({\mu^S}|\frac{\alpha^2}{4}{\nabla^{2}}|{\nu^S}) \\\\
C> T & = & -\frac{1}{2} ({\mu^L}|{\nabla^{2}}|{\nu^S})
C>       -  \frac{1}{2} ({\mu^S}|{\nabla^{2}}|{\nu^L})
C>       +  \frac{1}{2} ({\mu^S}|{\nabla^{2}}|{\nu^S}) \\\\
C> V & = & ({\mu^L}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu^L})
C>     - \frac{\alpha^2}{4} ({\mu^S}|\nabla\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}\cdot\nabla|{\nu^S}) \\\\
C> \f}
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1estv}
c:tex-This is an internal routine that most of the external 1 electron
c:tex-routines call.  This is the actual workhorse routine.
c:tex-This routine computes the 1 electron integrals S, T, and V:
c:tex-\begin{eqnarray*}
c:tex-S & = & ({\mu}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-T & = & ({\mu}|-\frac{1}{2}\nabla^{2}|{\nu}) \\
c:tex-  & = & -\frac{1}{2}\int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\nabla^{2}(r_{1})
c:tex-        g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
c:tex-{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-\end{eqnarray*}
c:tex-
c:tex-If an ECP is defined then the ECP integral contributions are summed
c:tex-directly into the V integrals.
c:tex-
c:tex-If a relativistic basis is defined then the one-electron integrals
c:tex-for the case where both shells are relativistic are modified to
c:tex-\begin{eqnarray*}
c:tex-S & = & ({\mu^L}|{\nu^L})
c:tex-      - ({\mu^S}|\frac{\alpha^2}{4}{\nabla^{2}}|{\nu^S}) \\
c:tex-T & = & -\frac{1}{2} ({\mu^L}|{\nabla^{2}}|{\nu^S})
c:tex-      -  \frac{1}{2} ({\mu^S}|{\nabla^{2}}|{\nu^L})
c:tex-      +  \frac{1}{2} ({\mu^S}|{\nabla^{2}}|{\nu^S}) \\
c:tex-V & = & ({\mu^L}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu^L})
c:tex-    - \frac{\alpha^2}{4} ({\mu^S}|\nabla\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}\cdot\nabla|{\nu^S}) \\
c:tex-\end{eqnarray*}
c:tex- 
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1estv(i_basis,ish,j_basis,jsh,lscr,scr,lstv,S,T,V,
cc AJL/Begin
     &    doS,doT,doV,msg)
cc AJL/End
c:tex-\end{verbatim}
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geom.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
#include "util.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e
      external cando_nw_1e
      external cando_nw
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      integer lstv               !< [Input] length of integral buffer
      double precision scr(lscr) !< [Scratch] scratch array
      double precision S(lstv)   !< [Output] overlap integrals
      double precision T(lstv)   !< [Output] kinetic energy integrals
      double precision V(lstv)   !< [Output] potential energy integrals
      logical doS                !< [Input] flag for overlap integrals
      logical doT                !< [Input] flag for kinetic energy integrals
      logical doV                !< [Input] flag for potential energy integrals
c:tex-\end{verbatim}
c::local
      logical ohnd_ok, onw_ok, onw1e_ok
      logical shells_ok, orel, oirel, ojrel, oNR, canAB
      integer i_geom, j_geom, ibas, jbas, ucont, uconts
      integer lbas, sbas, abas, isbas, jsbas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_icfpS
c
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
c
      integer WarnP
      save WarnP
      data WarnP /0/
cc AJL/Begin
cc Soft-coded the value of msg
      character*(*) msg
cc AJL/End
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_1eov'))
     &       call errquit('int_1eov: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1eov: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
      canAB = .false.
c
      if (dyall_mod_dir) then
c
c     get basis set handles; relativistic integral option only valid
c     if both ibas and jbas are the ao basis.
c
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        orel = ibas .eq. abas .and. jbas .eq. abas
      end if
c
c   i shell
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1estv: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          ibas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        isbas = ibas
      end if
c
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS = infbs_cont(CONT_ICFP ,uconts,isbas)
c
c   j shell
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,jbas) .ne. 0)) then
        ojrel = .true.
        jsbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1estv: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          jbas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        jsbas = jbas
      end if
c
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_1eov: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif

      ohnd_ok = cando_hnd_1e(i_basis,ish,0)
     &    .and. cando_hnd_1e(j_basis,jsh,0)
     &    .and. (.not.geom_any_finuc (i_geom))
     &    .and. (.not.geom_any_finuc (j_geom))
      onw_ok = cando_nw(i_basis,ish,0) .and. cando_nw(j_basis,jsh,0)
      onw1e_ok = cando_nw_1e(i_basis,ish,0)
     &    .and. cando_nw_1e(j_basis,jsh,0)
      if (orel) then
        call rel_onel (
     &      coords(1,i_cent,i_geom),
     &      dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfpS,isbas)),i_prim,i_gen,Li,
     &      coords(1,j_cent,j_geom),
     &      dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfpS,jsbas)),j_prim,j_gen,Lj,
     &      coords(1,1,i_geom),charge(1,i_geom),
     &      geom_invnucexp(1,i_geom),ncenter(i_geom),
     &      S,T,V,lstv,doS,doT,doV,canAB,onw_ok,ohnd_ok,oNR,.false.,
     &      scr,lscr,rel_dbg,rel_typ)
      else if (ohnd_ok) then
        call hnd_stvint(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      j_prim, j_gen, Lj,
     &      coords(1,1,i_geom),charge(1,i_geom),ncenter(i_geom),
     &      S,T,V,lstv,doS,doT,doV,scr,lscr)
c
      elseif (onw1e_ok) then
        call int_hf1sp(
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li, i_cent,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),
     &        j_prim, j_gen, Lj, j_cent,
     &        coords(1,1,i_geom),charge(1,i_geom),
     &        geom_invnucexp(1,i_geom),ncenter(i_geom),
     &        S,T,V,lstv,doS,doT,doV,canAB,.false.,
cc AJL/Begin
c     &        scr,lscr,'int_1eov')
     &        scr,lscr,msg)
cc AJL/End
      elseif (onw_ok) then
        call hf1(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)), j_prim, j_gen, Lj,
     &      coords(1,1,i_geom),charge(1,i_geom),
     &      geom_invnucexp(1,i_geom),ncenter(i_geom),
     &      S,T,V,lstv,doS,doT,doV,canAB,.false.,
     &      scr,lscr)
      else
        call errquit('int_1eov: could not do hnd, sp or nw integrals',
     &                0, INT_ERR)
      endif
c
*     We now have the cartesian integral block(s)  (jlo:jhi,ilo:ihi)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
c
      if (li.eq.-1) i_gen = 1
      if (lj.eq.-1) j_gen = 1
c
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*... transform both i and j integrals
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
c
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = j_nbf_x
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
      else
        call errquit(
     &      'int_1eov: should never reach transform blocked else',911,
     &       INT_ERR)
      endif
      return
c
      end
C> @}
c----------------------------------------------------------------------
c
c     Wrapper routines to compute different classes of integrals. 
c     All of these call int_1estv.
c
c----------------------------------------------------------------------
*
C> \ingroup nwint
C> @{
C> \brief Compute the overlap integrals
C>
C> Compute the overlap integrals between two shells. These integrals
C> are defined as:
C> \f{eqnarray*}{
C> S & = & ({\mu}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1eov}
c:tex-This routine computes the 1 electron overlap integrals ($S$):
c:tex-\begin{eqnarray*}
c:tex-S & = & ({\mu}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1eov(i_basis,ish,j_basis,jsh,lscr,scr,lov,Ov)
#include "util.fh"
c:tex-\end{verbatim}
c
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lov                !< [Input] length of Ov buffer
      double precision Ov(lov)   !< [Output] overlap integrals
c:tex-\end{verbatim}
c
      call int_1estv(i_basis,ish,j_basis,jsh,lscr,scr,lov,Ov,scr,scr,
cc AJL/Begin
c     &    .true.,.false.,.false.)
     &    .true.,.false.,.false.,'int_1eov')
c          doS    doT     doV     msg
cc AJL/End
      return
      end
c----------------------------------------------------------------------
*
C> \brief Calculates the 1-electron kinetic energy integrals
C>
C> This routine calculates the 1-electron kinetic energy integrals 
C> between a pair of specified shells. The integrals are defined as:
C> \f{eqnarray*}{
C> T & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}\nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1eke}
c:tex-This routine computes the 1 electron kinetic energy integrals, ($T$).:
c:tex-\begin{eqnarray*}
c:tex-T & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}\nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1eke(i_basis,ish,j_basis,jsh,lscr,scr,lke,Ke)
#include "util.fh"
c:tex-\end{verbatim}
c
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lke                !< [Input] length of Ke buffer
      double precision Ke(lke)   !< [Output] kinetic energy integrals
c:tex-\end{verbatim}
c
      call int_1estv(i_basis,ish,j_basis,jsh,lscr,scr,lke,scr,Ke,scr,
cc AJL/Begin
c     &    .false.,.true.,.false.)
     &    .false.,.true.,.false.,'int_1eke')
c          doS     doT    doV     msg
cc AJL/End
      return
      end
c----------------------------------------------------------------------
*
C> \brief Computes the nuclear attraction integrals
C>
C> This routine calculates the 1-electron nuclear attraction integrals
C> for a pair of shells. If any ECPs are present then the ECP terms 
C> are directly summed into these integrals.
C>
C> The integrals are defined as:
C> \f{eqnarray*}{
C> V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
C> {-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1epe}
c:tex-This routine computes the 1 electron potential integrals, ($V$):
c:tex-\begin{eqnarray*}
c:tex-V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
c:tex-{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-If an ECP is defined then the ECP integral contributions are summed
c:tex-directly into the V integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1epe(i_basis,ish,j_basis,jsh,lscr,scr,lpe,Pe)
#include "util.fh"
c:tex-\end{verbatim}
c
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lpe                !< [Input] length of Pe buffer
      double precision Pe(lpe)   !< [Output] kinetic energy integrals
c:tex-\end{verbatim}
c
      call int_1estv(i_basis,ish,j_basis,jsh,lscr,scr,lpe,scr,scr,Pe,
cc AJL/Begin
c     &    .false.,.false.,.true.)
     &    .false.,.false.,.true.,'int_1epe')
c          doS     doT     doV    msg
cc AJL/End
      return
      end
C> @}
cc AJL/Begin/SPIN ECPs
c----------------------------------------------------------------------
*
C> \brief Computes the nuclear attraction integrals
C>
C> This routine calculates the 1-electron nuclear attraction integrals
C> for a pair of shells. If any ECPs are present then the ECP terms 
C> are directly summed into these integrals.
C>
C> The integrals are defined as:
C> \f{eqnarray*}{
C> V & = &
C({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}
Cg_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
C> {-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1epe}
c:tex-This routine computes the 1 electron potential integrals, ($V$):
c:tex-\begin{eqnarray*}
c:tex-V & = &({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
c:tex-{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-If a *BETA* ECP is defined then the ECP integral contributions are summed
c:tex-directly into the V integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1epe_beta(i_basis,ish,j_basis,jsh,lscr,scr,lpe,Pe)
#include "util.fh"
c:tex-\end{verbatim}
c
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lpe                !< [Input] length of Pe buffer
      double precision Pe(lpe)   !< [Output] kinetic energy integrals
c:tex-\end{verbatim}
c
      call int_1estv(i_basis,ish,j_basis,jsh,lscr,scr,lpe,scr,scr,Pe,
     &    .false.,.false.,.true.,'int_1epe_beta')
c          doS     doT     doV    msg
      return
      end
cc AJL/End
C> @}
c----------------------------------------------------------------------
*
C> \ingroup nwint
C> @{
C> \brief Compute the 1-electron Hamiltonian
C>
C> This routine computes the 1 electron hamiltonian \f$H1\f$.
C> \f{eqnarray*}{
C> H1 & = & T + V      \\\\
C> T  & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\\\
C>    & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}
C> \nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> V  & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\\\
C>    & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
C> {-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> 
C> If an ECP is defined then the ECP integral contributions are summed
C> directly into the \f$H1\f$ integrals.
C> 
C> If a relativistic basis is defined then the one-electron integrals for
C> the case where both shells are relativistic are the modified integrals.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1eh1}
c:tex-This routine computes the 1 electron hamiltonian, ($H1$).
c:tex-\begin{eqnarray*}
c:tex-H1 & = & T + V      \\
c:tex-T  & = & ({\mu}|\frac{-1}{2}\nabla^{2}|{\nu}) \\
c:tex-   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}
c:tex-\nabla^{2}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-V  & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-   & = & \int_{-\infty}^{\infty} g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac
c:tex-{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-If an ECP is defined then the ECP integral contributions are summed
c:tex-directly into the $H1$ integrals.
c:tex-
c:tex-If a relativistic basis is defined then the one-electron integrals for
c:tex-the case where both shells are relativistic are the modified integrals.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1eh1(i_basis,ish,j_basis,jsh,lscr,scr,lh1,H1)
#include "util.fh"
c:tex-\end{verbatim}
#include "apiP.fh"
c
c::functions
      logical cando_hnd_1e
      external cando_hnd_1e
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lh1                !< [Input] length of H1 buffer.
      double precision H1(lh1)   !< [Output] one electron
c:tex-\end{verbatim}
c
      call int_1estv(i_basis,ish,j_basis,jsh,lscr-lh1,scr(lh1+1),lh1,
cc AJL/Begin
c     &    scr,scr,H1,.false.,.true.,.true.)
     &    scr,scr,H1,.false.,.true.,.true.,'int_1eh1')
c                     doS     doT    doV    msg
cc AJL/End
      do i = 1,lh1
        H1(i) = H1(i)+scr(i)
      end do
      return
      end
c----------------------------------------------------------------------
*
C> \brief Compute all 1-electron integrals
C>
C> This routine computes the 1 electron integrals \f$S\f$, \f$T\f$, and \f$V\f$:
C> \f{eqnarray*}{
C> S & = & ({\mu}|{\nu}) \\\\
C>   & = & \int_{{-}\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> T & = & ({\mu}|-\frac{1}{2}{\nabla^{2}}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}{\nabla^{2}}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\\\
C> V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-Z_{\alpha}}
C> {|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C>
c:tex-
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1eall}
c:tex-This routine computes the 1 electron integrals S, T, and V:
c:tex-\begin{eqnarray*}
c:tex-S & = & ({\mu}|{\nu}) \\
c:tex-  & = & \int_{{-}\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-T & = & ({\mu}|-\frac{1}{2}{\nabla^{2}}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\frac{-1}{2}{\nabla^{2}}(r_{1})g_{\nu}(X_{\nu},r_{1})dr_{1} \\
c:tex-V & = & ({\mu}|\sum_{\alpha}\frac{-Z_{\alpha}}{|r_{1}-R_{\alpha}|}|{\nu}) \\
c:tex-  & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-Z_{\alpha}}
c:tex-{|r_{1}-R_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
c:tex-\end{eqnarray*}
c:tex-
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim} 
      subroutine int_1eall(i_basis,ish,j_basis,jsh,lscr,scr,lstv,S,T,V)
#include "util.fh"
c:tex-\end{verbatim}
c
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lstv               !< [Input] length of one electron buffers
      double precision T(lstv)   !< [Output] kinetic integral buffer
      double precision V(lstv)   !< [Output] potential integral buffer
      double precision S(lstv)   !< [Output] overlap integral buffer
c:tex-\end{verbatim}
c
      call int_1estv(i_basis,ish,j_basis,jsh,lscr,scr,lstv,S,T,V,
cc AJL/Begin
c     &    .true.,.true.,.true.)
     &    .true.,.true.,.true.,'int_1eall')
c          doS    doT    doV    msg
cc AJL/End
      return
      end
c----------------------------------------------------------------------
C> \brief Compute the 1-electron COSMO embedding potential 
C>
C> This routine computes the 1 electron integrals \f$V\f$:
C> \f{eqnarray*}{
C> V & = & ({\mu}|\sum_{\alpha}\frac{-q_{\alpha}}{|r_{1}-r_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-q_{\alpha}}
C> {|r_{1}-r_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> where \f$ q_\alpha \f$ is a point charge situated at position \f$ r_\alpha \f$.
C> The COSMO embedding point charges are accessed in int_1eefc.
C>
      subroutine int_1epot(i_basis,ish,j_basis,jsh,lscr,scr,lve,Ve)
#include "util.fh"
c
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lve                !< [Input] length of Ve buffer
      double precision Ve(lve)   !< [Output] potential energy integrals
c
      call int_1eefc(i_basis,ish,j_basis,jsh,lscr,scr,lve,scr,scr,Ve,
     &    .false.,.false.,.true.)
c          doS     doT     doV
      return
      end
c----------------------------------------------------------------------
C> \brief Compute the 1-electron integrals but the potential is due to
C> COSMO embedding charges
C>
C> This routine computes the 1 electron integrals \f$V\f$:
C> \f{eqnarray*}{
C> V & = & ({\mu}|\sum_{\alpha}\frac{-q_{\alpha}}{|r_{1}-r_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-q_{\alpha}}
C> {|r_{1}-r_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> where \f$ q_\alpha \f$ is a point charge situated at position \f$ r_\alpha \f$.
C> The location of the COSMO embedding point charges is specified in the common
C> block include from prop.fh.
C>
      subroutine int_1eefc(i_basis,ish,j_basis,jsh,lscr,scr,lstv,S,T,V,
     &    doS,doT,doV)
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bq.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
#include "util.fh"
C
C The following include block contains the geometry handle for the
C geometry and charges of the solvent for COSMO. This is a clear violation
C of the integral API, but for now the cleanest solution
C
#include "prop.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e
      external cando_nw_1e
      external cando_nw
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      integer lstv               !< [Input] length of integral buffer
      double precision scr(lscr) !< [Scratch] scratch array
      double precision S(lstv)   !< [Output] overlap integrals
      double precision T(lstv)   !< [Output] kinetic energy integrals
      double precision V(lstv)   !< [Output] potential energy integrals
      logical doS                !< [Input] flag for overlap integrals
      logical doT                !< [Input] flag for kinetic energy integrals
      logical doV                !< [Input] flag for potential energy integrals
c:tex-\end{verbatim}
c::local
      logical ohnd_ok, onw_ok, onw1e_ok
      logical shells_ok, orel, oirel, ojrel, oNR, canAB
      integer i_geom, j_geom, ibas, jbas, ucont, uconts
      integer lbas, sbas, abas, isbas, jsbas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_icfpS
      integer iefc_c, iefc_q, nefc
      integer k_efc_inv
c
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
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
      if (.not.int_chk_init('int_1eefc'))
     &       call errquit('int_1eefc: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1eefc: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
      canAB = .false.
c
      if (dyall_mod_dir) then
c
c     get basis set handles; relativistic integral option only valid
c     if both ibas and jbas are the ao basis.
c
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        orel = ibas .eq. abas .and. jbas .eq. abas
      end if
c
c   i shell
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1eefc: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          ibas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        isbas = ibas
      end if
c
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS = infbs_cont(CONT_ICFP ,uconts,isbas)
c
c   j shell
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,jbas) .ne. 0)) then
        ojrel = .true.
        jsbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1eefc: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          jbas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        jsbas = jbas
      end if
c
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_1eefc: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
c     ----- current cosmo restriction ... -----
c
      if(i_geom.ne.j_geom) then
         write(luout,*) 'i_geom andr j_geom is noteq'
         call errquit('int_1eefc: i_geom and j_geom must be the same',0,
     &         INT_ERR)
      endif
c
      if (.not.bq_ncenter(cosmo_bq_efc,nefc))
     +  call errquit("int_1eefc: bq_ncenter failed",nefc,uerr)
      if (.not.bq_index_coord(cosmo_bq_efc,iefc_c))
     +  call errquit("int_1eefc: bq_index_coord failed",iefc_c,uerr)
      if (.not.bq_index_charge(cosmo_bq_efc,iefc_q))
     +  call errquit("int_1eefc: bq_index_charge failed",iefc_q,uerr)
c
c     Pick up the array we created and zeroed in cosmo_initialize
c
      if (.not.bq_index_charge(cosmo_bq_invnuc,k_efc_inv))
     +  call errquit("int_1eefc: bq_index_charge failed",iefc_q,uerr)
c
      ohnd_ok = cando_hnd_1e(i_basis,ish,0)
     &    .and. cando_hnd_1e(j_basis,jsh,0)
     &    .and. (.not.geom_any_finuc (i_geom))
     &    .and. (.not.geom_any_finuc (j_geom))
      onw_ok = cando_nw(i_basis,ish,0) .and. cando_nw(j_basis,jsh,0)
      onw1e_ok = cando_nw_1e(i_basis,ish,0)
     &    .and. cando_nw_1e(j_basis,jsh,0)
c
      if (orel) then
        call rel_onel (
     &      coords(1,i_cent,i_geom),
     &      dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfpS,isbas)),i_prim,i_gen,Li,
     &      coords(1,j_cent,j_geom),
     &      dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfpS,jsbas)),j_prim,j_gen,Lj,
     &      dbl_mb(iefc_c),dbl_mb(iefc_q),
     &      dbl_mb(k_efc_inv),nefc,
     &      S,T,V,lstv,doS,doT,doV,canAB,onw_ok,ohnd_ok,oNR,.false.,
     &      scr,lscr,rel_dbg,rel_typ)
      else if (ohnd_ok) then
        call hnd_stvint(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      j_prim, j_gen, Lj,
     &      dbl_mb(iefc_c),dbl_mb(iefc_q),nefc,
     &      S,T,V,lstv,doS,doT,doV,scr,lscr)
c
      elseif (onw1e_ok) then
          call int_hf1sp(
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li, i_cent,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),
     &        j_prim, j_gen, Lj, j_cent,
     &        dbl_mb(iefc_c),dbl_mb(iefc_q),
     &        dbl_mb(k_efc_inv),nefc,
     &        S,T,V,lstv,doS,doT,doV,canAB,.false.,
     &        scr,lscr,'int_1eefc')
      elseif (onw_ok) then
        call hf1(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)), j_prim, j_gen, Lj,
     &      dbl_mb(iefc_c),dbl_mb(iefc_q),
     &      dbl_mb(k_efc_inv),nefc,
     &      S,T,V,lstv,doS,doT,doV,canAB,.false.,
     &      scr,lscr)
      else
        call errquit('int_1eefc: could not do hnd, sp or nw integrals',
     &                0, INT_ERR)
      endif
c
*     We now have the cartesian integral block(s)  (jlo:jhi,ilo:ihi)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
c
      if (li.eq.-1) i_gen = 1
      if (lj.eq.-1) j_gen = 1
c
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*... transform both i and j integrals
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
c
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = j_nbf_x
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
      else
        call errquit(
     &      'int_1eefc: should never reach transform blocked else',911,
     &          INT_ERR)
      endif
      return
c
      end

c----------------------------------------------------------------------
C> \brief Compute the 1-electron Bq embedding potential 
C>
C> This routine computes the 1 electron integrals \f$V\f$:
C> \f{eqnarray*}{
C> V & = & ({\mu}|\sum_{\alpha}\frac{-q_{\alpha}}{|r_{1}-r_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-q_{\alpha}}
C> {|r_{1}-r_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> where \f$ q_\alpha \f$ is a point charge situated at position \f$ r_\alpha \f$.
C> The Bq embedding point charges are accessed in int_1eefc1.
C>
      subroutine int_1epot1(i_basis,ish,j_basis,jsh,lscr,scr,lpe,Pe)
#include "util.fh"
c
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lpe                !< [Input] length of Pe buffer
      double precision Pe(lpe)   !< [Output] potential energy integrals
c
      call int_1eefc1(i_basis,ish,j_basis,jsh,lscr,scr,lpe,scr,scr,Pe,
     &    .false.,.false.,.true.)
c          doS     doT     doV
      return
      end
c----------------------------------------------------------------------
C> \brief Compute the 1-electron integrals but the potential is due to
C> Bq embedding charges
C>
C> This routine computes the 1 electron integrals \f$V\f$:
C> \f{eqnarray*}{
C> V & = & ({\mu}|\sum_{\alpha}\frac{-q_{\alpha}}{|r_{1}-r_{\alpha}|}|{\nu}) \\\\
C>   & = & \int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-q_{\alpha}}
C> {|r_{1}-r_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> where \f$ q_\alpha \f$ is a point charge situated at position \f$ r_\alpha \f$.
C> The location of the Bq embedding point charges is specified in a 
C> Bq instance that has been activated elsewhere.
C>
      subroutine int_1eefc1(i_basis,ish,j_basis,jsh,lscr,scr,lstv,S,T,V,
     &    doS,doT,doV)
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geom.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
#include "util.fh"
C
C The following include block contains the geometry handle for the
C geometry and charges of the solvent for COSMO. This is a clear violation
C of the integral API, but for now the cleanest solution
C
#include "prop.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e
      external cando_nw_1e
      external cando_nw
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      integer lstv               !< [Input] length of integral buffer
      double precision scr(lscr) !< [Scratch] scratch array
      double precision S(lstv)   !< [Output] overlap integrals
      double precision T(lstv)   !< [Output] kinetic energy integrals
      double precision V(lstv)   !< [Output] potential energy integrals
      logical doS                !< [Input] flag for overlap integrals
      logical doT                !< [Input] flag for kinetic energy integrals
      logical doV                !< [Input] flag for potential energy integrals
c:tex-\end{verbatim}
c::local
      logical ohnd_ok, onw_ok, onw1e_ok
      logical shells_ok, orel, oirel, ojrel, oNR, canAB
      integer i_geom, j_geom, ibas, jbas, ucont, uconts
      integer lbas, sbas, abas, isbas, jsbas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_icfpS
c
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
c
      integer WarnP
      save WarnP
      data WarnP /0/
c 
      integer bq_ncent
      integer i_qbq,i_cbq,h_xbq0,i_xbq0
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c check initialization and shells
c
      if (.not.int_chk_init('int_1eefc1'))
     &       call errquit('int_1eefc1: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1eefc1: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
      canAB = .false.
c
      if (dyall_mod_dir) then
c
c     get basis set handles; relativistic integral option only valid
c     if both ibas and jbas are the ao basis.
c
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        orel = ibas .eq. abas .and. jbas .eq. abas
      end if
c
c   i shell
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1eefc1: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          ibas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        isbas = ibas
      end if
c
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS = infbs_cont(CONT_ICFP ,uconts,isbas)
c
c   j shell
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,jbas) .ne. 0)) then
        ojrel = .true.
        jsbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1eefc1: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          jbas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        jsbas = jbas
      end if
c
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_1eefc1: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
c     ----- current cosmo restriction ... -----
c
      if(i_geom.ne.j_geom) then
         write(luout,*) 'i_geom andr j_geom is noteq'
         call errquit('int_1eefc1: i_geom and j_geom must be the same',
     &                0,INT_ERR)
      endif

      ohnd_ok = cando_hnd_1e(i_basis,ish,0)
     &    .and. cando_hnd_1e(j_basis,jsh,0)
     &    .and. (.not.geom_any_finuc (i_geom))
     &    .and. (.not.geom_any_finuc (j_geom))
      onw_ok = cando_nw(i_basis,ish,0) .and. cando_nw(j_basis,jsh,0)
      onw1e_ok = cando_nw_1e(i_basis,ish,0)
     &    .and. cando_nw_1e(j_basis,jsh,0)
c
c     get external charges here (MV)
c     ------------------------------
      if(.not.geom_extbq_on())
     >   call errquit('int_1eefc1:no active bqs',0,0) 
      bq_ncent = geom_extbq_ncenter()
      i_cbq = geom_extbq_coord()
      i_qbq = geom_extbq_charge()
      if(.not.ma_push_get(mt_dbl,bq_ncent,'xbq',h_xbq0,i_xbq0))
     +     call errquit( 'int_1eefc1',0,0)
      call dfill(bq_ncent,0.d0,dbl_mb(i_xbq0),1)
c
      if (orel) then
        call rel_onel (
     &      coords(1,i_cent,i_geom),
     &      dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfpS,isbas)),i_prim,i_gen,Li,
     &      coords(1,j_cent,j_geom),
     &      dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfpS,jsbas)),j_prim,j_gen,Lj,
c     &      coords(1,1,cosmo_geom_efc),charge(1,cosmo_geom_efc),
c     &      geom_invnucexp(1,i_geom),ncenter(cosmo_geom_efc),
     &      dbl_mb(i_cbq),dbl_mb(i_qbq),
     &      dbl_mb(i_xbq0),bq_ncent,
     &      S,T,V,lstv,doS,doT,doV,canAB,onw_ok,ohnd_ok,oNR,.false.,
     &      scr,lscr,rel_dbg,rel_typ)
      else if (ohnd_ok) then
        call hnd_stvint(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      j_prim, j_gen, Lj,
c     &      coords(1,1,cosmo_geom_efc),charge(1,cosmo_geom_efc),
c     &      ncenter(cosmo_geom_efc),
     &      dbl_mb(i_cbq),dbl_mb(i_qbq),
     &      bq_ncent,
     &      S,T,V,lstv,doS,doT,doV,scr,lscr)
c
      elseif (onw1e_ok) then
          call int_hf1sp(
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li, i_cent,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),
     &        j_prim, j_gen, Lj, j_cent,
c     &        coords(1,1,cosmo_geom_efc),charge(1,cosmo_geom_efc),
c     &        geom_invnucexp(1,i_geom),ncenter(cosmo_geom_efc),
     &        dbl_mb(i_cbq),dbl_mb(i_qbq),
     &        dbl_mb(i_xbq0),bq_ncent,
     &        S,T,V,lstv,doS,doT,doV,canAB,.false.,
     &        scr,lscr,'int_1eefc')
      elseif (onw_ok) then
        call hf1(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)), j_prim, j_gen, Lj,
c     &      coords(1,1,cosmo_geom_efc),charge(1,cosmo_geom_efc),
c     &      geom_invnucexp(1,i_geom),ncenter(cosmo_geom_efc),
     &      dbl_mb(i_cbq),dbl_mb(i_qbq),
     &      dbl_mb(i_xbq0),bq_ncent,
     &      S,T,V,lstv,doS,doT,doV,canAB,.false.,
     &      scr,lscr)
      else
        call errquit('int_1eefc1: could not do hnd, sp or nw integrals',
     &                0, INT_ERR)
      endif
c
      if(.not.ma_pop_stack(h_xbq0))
     +     call errquit( 'int_1eefc1',0,0)
c
*     We now have the cartesian integral block(s)  (jlo:jhi,ilo:ihi)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
c
      if (li.eq.-1) i_gen = 1
      if (lj.eq.-1) j_gen = 1
c
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*... transform both i and j integrals
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
c
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = j_nbf_x
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
      else
        call errquit(
     &      'int_1eefc1: should never reach transform blocked else',911,
     &          INT_ERR)
      endif
      return
c
      end
C>
C> @}
cc AJL/BEGIN/FDE
c----------------------------------------------------------------------
C> \brief Compute the 1-electron Bq embedding potential 
C>
C> This routine computes the 1 electron integrals \f$V\f$:
C> \f{eqnarray*}{
C> V & = &
C({\mu}|\sum_{\alpha}\frac{-q_{\alpha}}{|r_{1}-r_{\alpha}|}|{\nu}) \\\\
C>   & = &
C\int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-q_{\alpha}}
C> {|r_{1}-r_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> where \f$ q_\alpha \f$ is a point charge situated at position \f$
Cr_\alpha \f$.
C> The Bq embedding point charges are accessed in int_1eefc1.
C>
      subroutine int_1efde1(i_basis,ish,j_basis,jsh,fde_basis,lscr,scr,
     &                      lpe,Pe)
#include "util.fh"
c
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer fde_basis !< [Input] FDE basis set handle
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lpe                !< [Input] length of Pe buffer
      double precision Pe(lpe)   !< [Output] potential energy integrals
c
c      call int_1epe(i_basis,ish,j_basis,jsh,lscr,scr,lpe,scr,
c     &     scr,Pe,.false.,.false.,.true.,'pot')
      call int_1efde(i_basis,ish,j_basis,jsh,fde_basis,lscr,scr,lpe,scr,
     &     scr,Pe,.false.,.false.,.true.)
c                   doS     doT     doV
      return
      end
c----------------------------------------------------------------------
C> \brief Compute the 1-electron integrals but the potential is due to
C> Bq embedding charges
C>
C> This routine computes the 1 electron integrals \f$V\f$:
C> \f{eqnarray*}{
C> V & = &
C({\mu}|\sum_{\alpha}\frac{-q_{\alpha}}{|r_{1}-r_{\alpha}|}|{\nu}) \\\\
C>   & = &
C\int_{-\infty}^{\infty}g_{\mu}(X_{\mu},r_{1})\sum_{\alpha}\frac{-q_{\alpha}}
C> {|r_{1}-r_{\alpha}|}g_{\nu}(X_{\nu},r_{1})dr_{1}
C> \f}
C> where \f$ q_\alpha \f$ is a point charge situated at position \f$
Cr_\alpha \f$.
C> The location of the Bq embedding point charges is specified in a 
C> Bq instance that has been activated elsewhere.
C>
      subroutine int_1efde(i_basis,ish,j_basis,jsh,fde_basis,lscr,scr,
     &    lstv,S,T,V,doS,doT,doV)
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geom.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "rel_nwc.fh"
#include "util.fh"
C
C The following include block contains the geometry handle for the
C geometry and charges of the solvent for COSMO. This is a clear
C violation
C of the integral API, but for now the cleanest solution
C
#include "prop.fh"
c::external subroutines used
c... errquit
c::functions
      logical cando_hnd_1e
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_hnd_1e
      external cando_nw_1e
      external cando_nw
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer fde_basis !< [Input] FDE basis set handle
      integer lscr    !< [Input] length of scratch array
      integer lstv               !< [Input] length of integral buffer
      double precision scr(lscr) !< [Scratch] scratch array
      double precision S(lstv)   !< [Output] overlap integrals
      double precision T(lstv)   !< [Output] kinetic energy integrals
      double precision V(lstv)   !< [Output] potential energy integrals
      logical doS                !< [Input] flag for overlap integrals
      logical doT                !< [Input] flag for kinetic energy integrals
      logical doV                !< [Input] flag for potential energy integrals
c:tex-\end{verbatim}
c::local
      logical ohnd_ok, onw_ok, onw1e_ok
      logical shells_ok, orel, oirel, ojrel, oNR, canAB
      integer i_geom, j_geom, ibas, jbas, ucont, uconts
      integer lbas, sbas, abas, isbas, jsbas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_icfpS
c
      integer fde_bas, fde_geom
c
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s
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
      if (.not.int_chk_init('int_1efde'))
     &       call errquit('int_1efde: int_init was not called' ,0,
     &       INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1efde: invalid contraction/shell',0,
     &       BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      fde_bas = fde_basis + BASIS_HANDLE_OFFSET
      oNR = .true.
      oirel = .false.
      ojrel = .false.
      orel = .false.
      canAB = .false.
c
      if (dyall_mod_dir) then
c
c     get basis set handles; relativistic integral option only valid
c     if both ibas and jbas are the ao basis.
c
        lbas = lc_bsh + BASIS_HANDLE_OFFSET
        sbas = sc_bsh + BASIS_HANDLE_OFFSET
        abas = ao_bsh + BASIS_HANDLE_OFFSET
        orel = ibas .eq. abas .and. jbas .eq. abas
      end if
c
c   i shell
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,ibas) .ne. 0)) then
        oirel = .true.
        isbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1efde: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          ibas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        isbas = ibas
      end if
c
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS = infbs_cont(CONT_ICFP ,uconts,isbas)
c
c   j shell
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
c
c     check for relativistic shell
c
      if (orel .and. (infbs_cont(CONT_RELLS ,ucont,jbas) .ne. 0)) then
        ojrel = .true.
        jsbas = sbas
        uconts = ao_to_ls(ucont)
        if (uconts .eq. 0) call errquit (
     &      'int_1efde: no relativistic pointer',911, INT_ERR)
        if (nesc_1e_approx) then
          jbas = lbas
          ucont = uconts
        end if
      else
        uconts = ucont
        jsbas = jbas
      end if
c
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS = infbs_cont(CONT_ICFP ,uconts,jsbas)
c
      oNR = .not.(oirel.and.ojrel)
      orel = oirel.or.ojrel
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_1efde: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
c     ----- current cosmo restriction ... -----
c
      if(i_geom.ne.j_geom) then
         write(luout,*) 'i_geom and j_geom is not eq'
         call errquit('int_1efde: i_geom and j_geom must be the same',
     &                0,INT_ERR)
      endif

      ohnd_ok = cando_hnd_1e(i_basis,ish,0)
     &    .and. cando_hnd_1e(j_basis,jsh,0)
     &    .and. (.not.geom_any_finuc (i_geom))
     &    .and. (.not.geom_any_finuc (j_geom))
      onw_ok = cando_nw(i_basis,ish,0) .and. cando_nw(j_basis,jsh,0)
      onw1e_ok = cando_nw_1e(i_basis,ish,0)
     &    .and. cando_nw_1e(j_basis,jsh,0)
c
c     get external charges here (MV)
c     ------------------------------
c      if(.not.geom_extbq_on())
c     >   call errquit('int_1eefc1:no active bqs',0,0) 
c      bq_ncent = geom_extbq_ncenter()
c      i_cbq = geom_extbq_coord()
c      i_qbq = geom_extbq_charge()
c      if(.not.ma_push_get(mt_dbl,bq_ncent,'xbq',h_xbq0,i_xbq0))
c     +     call errquit( 'int_1eefc1',0,0)
c      call dfill(bq_ncent,0.d0,dbl_mb(i_xbq0),1)
c
      fde_geom  = ibs_geom(fde_bas)

c      write(*,*) 'FDE_GEOM ', fde_geom
c      write(*,*) 'COORDS ', coords(1,1,fde_geom)
c      write(*,*) 'CHARGE ', charge(1,fde_geom)
c      write(*,*) 'INVNUCLEARCHARGE', geom_invnucexp(1,fde_geom)
c      write(*,*) 'NCENTER ', ncenter(fde_geom)

c      call util_flush(luout)

      if (orel) then
        call rel_onel (
     &      coords(1,i_cent,i_geom),
     &      dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfpS,isbas)),i_prim,i_gen,Li,
     &      coords(1,j_cent,j_geom),
     &      dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfpS,jsbas)),j_prim,j_gen,Lj,
     &      coords(1,1,fde_geom),charge(1,fde_geom),
     &      geom_invnucexp(1,fde_geom),ncenter(fde_geom),
     &      S,T,V,lstv,doS,doT,doV,canAB,onw_ok,ohnd_ok,oNR,.false.,
     &      scr,lscr,rel_dbg,rel_typ)
      else if (ohnd_ok) then
        call hnd_stvint(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)),
     &      i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)),
     &      j_prim, j_gen, Lj,
     &      coords(1,1,fde_geom),charge(1,fde_geom),
     &      ncenter(fde_geom),
     &      S,T,V,lstv,doS,doT,doV,scr,lscr)
c
      elseif (onw1e_ok) then
          call int_hf1sp(
     &        coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &        dbl_mb(mb_exndcf(i_icfp,ibas)),
     &        i_prim, i_gen, Li, i_cent,
     &        coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &        dbl_mb(mb_exndcf(j_icfp,jbas)),
     &        j_prim, j_gen, Lj, j_cent,
     &        coords(1,1,fde_geom),charge(1,fde_geom),
     &        geom_invnucexp(1,fde_geom),ncenter(fde_geom),
     &        S,T,V,lstv,doS,doT,doV,canAB,.false.,
     &        scr,lscr,'int_1eefc')
      elseif (onw_ok) then
        call hf1(
     &      coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &      dbl_mb(mb_exndcf(i_icfp,ibas)), i_prim, i_gen, Li,
     &      coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &      dbl_mb(mb_exndcf(j_icfp,jbas)), j_prim, j_gen, Lj,
     &      coords(1,1,fde_geom),charge(1,fde_geom),
     &      geom_invnucexp(1,fde_geom),ncenter(fde_geom),
     &      S,T,V,lstv,doS,doT,doV,canAB,.false.,
     &      scr,lscr)
      else
        call errquit('int_1efde: could not do hnd, sp or nw integrals',
     &                0, INT_ERR)
      endif
c
c      if(.not.ma_pop_stack(h_xbq0))
c     +     call errquit( 'int_1eefc1',0,0)
c
*     We now have the cartesian integral block(s)  (jlo:jhi,ilo:ihi)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4. Since int_nbf_* arrays are set to the appropriate size.
c
      if (li.eq.-1) i_gen = 1
      if (lj.eq.-1) j_gen = 1
c
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*... transform both i and j integrals
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
c
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = j_nbf_x
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,0,j_gen,
     &      j_nbf_s,i_nbf_s,Li,i_gen,
     &      .false.)
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
        if (doS) call spcart_tran1e(S,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doT) call spcart_tran1e(T,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
        if (doV) call spcart_tran1e(V,scr,
     &      j_nbf_x,i_nbf_x,Lj,j_gen,
     &      j_nbf_s,i_nbf_s,0,i_gen,
     &      .false.)
      else
        call errquit(
     &      'int_1efde: should never reach transform blocked else',911,
     &          INT_ERR)
      endif
      return
c
      end
cc AJL/End
C>
C> @}
