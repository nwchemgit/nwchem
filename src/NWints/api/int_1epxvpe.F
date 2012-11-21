c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Computes the 1-electron double derivative potential energy integrals
C>
C> Computes the 1-electron double derivative potential energy integrals of
C> the form \f$pxVp\f$.
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{int\_1epxvpe}
c:tex-This routine computes the 1 electron double derivative potential energy integrals pxVp.
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_1epxvpe(i_basis,ish,j_basis,jsh,
     $                       lscr,scr,lpvpe,pvpe)
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
      double precision scr(lscr)     !< [Scratch] scratch array
      integer lpvpe                  !< [Input] length of pvpe buffer
      double precision pvpe(lpvpe,3) !< [Output] pxVp integrals
c:tex-\end{verbatim}
c::local
      logical shells_ok
      integer i_geom, j_geom, ibas, jbas, ucont
      integer isbas, jsbas
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent, i_icfpS
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent, j_icfpS
c
      logical any_spherical
      integer i_nbf_x, j_nbf_x
      integer i_nbf_s, j_nbf_s, icomp
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
      if (.not.int_chk_init('int_1epxvpe'))
     &       call errquit('int_1epxvpe: int_init was not called' ,0,
     &            INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1epxvpe: invalid contraction/shell',0,
     &        BASIS_ERR)
c
      ibas = i_basis + BASIS_HANDLE_OFFSET
      jbas = j_basis + BASIS_HANDLE_OFFSET
      isbas = ibas
      jsbas = jbas
c
      ucont   = (sf_ibs_cn2ucn(ish,ibas))
      Li      = infbs_cont(CONT_TYPE ,ucont,ibas)
      i_prim  = infbs_cont(CONT_NPRIM,ucont,ibas)
      i_gen   = infbs_cont(CONT_NGEN ,ucont,ibas)
      i_iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
      i_icfp  = infbs_cont(CONT_ICFP ,ucont,ibas)
      i_cent  = (sf_ibs_cn2ce(ish,ibas))
      i_geom  = ibs_geom(ibas)
      i_icfpS  = infbs_cont(CONT_ICFP ,ucont,isbas)
c
      ucont   = (sf_ibs_cn2ucn(jsh,jbas))
      Lj      = infbs_cont(CONT_TYPE ,ucont,jbas)
      j_prim  = infbs_cont(CONT_NPRIM,ucont,jbas)
      j_gen   = infbs_cont(CONT_NGEN ,ucont,jbas)
      j_iexp  = infbs_cont(CONT_IEXP ,ucont,jbas)
      j_icfp  = infbs_cont(CONT_ICFP ,ucont,jbas)
      j_cent  = (sf_ibs_cn2ce(jsh,jbas))
      j_geom  = ibs_geom(jbas)
      j_icfpS  = infbs_cont(CONT_ICFP ,ucont,jsbas)
c
      if (i_geom.ne.j_geom.and.WarnP.eq.0) then
        write(luout,*)
     &      'int_1epxvpe: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
      call rel_pvp (
     &     coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &     dbl_mb(mb_exndcf(i_icfp,ibas)),
     &     i_prim, i_gen, Li, i_cent,
     &     coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &     dbl_mb(mb_exndcf(j_icfp,jbas)),
     &     j_prim, j_gen, Lj, j_cent,
     &     coords(1,1,i_geom),charge(1,i_geom),
     &     geom_invnucexp(1,i_geom),ncenter(i_geom),
     &     pvpe,lpvpe,.false.,.false.,scr,lscr,'int_1epxvpe',0,3)
c
*     Ke now has the cartesian integral block  (jlo:jhi,ilo:ihi)
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
c    Transformation of each of the three components to sphericals if needed
c
      if (bas_spherical(ibas).and.bas_spherical(jbas)) then
*... transform both i and j integrals
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
c        
        do icomp = 1, 3
           call spcart_tran1e(pvpe(1,icomp),scr,
     &           j_nbf_x,i_nbf_x,Lj,j_gen,
     &           j_nbf_s,i_nbf_s,Li,i_gen,
     &           .false.)
        enddo
      else if (bas_spherical(ibas)) then
*.. transform on i component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = int_nbf_s(Li)
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = j_nbf_x
        do icomp = 1, 3
           call spcart_tran1e(pvpe(1,icomp),scr,
     &           j_nbf_x,i_nbf_x,0,j_gen,
     &           j_nbf_s,i_nbf_s,Li,i_gen,
     &           .false.)
        enddo
      else if (bas_spherical(jbas)) then
*.. transform on j component
        i_nbf_x = int_nbf_x(Li)
        i_nbf_s = i_nbf_x
        j_nbf_x = int_nbf_x(Lj)
        j_nbf_s = int_nbf_s(Lj)
        do icomp = 1, 3
           call spcart_tran1e(pvpe(1,icomp),scr,
     &           j_nbf_x,i_nbf_x,Lj,j_gen,
     &           j_nbf_s,i_nbf_s,0,i_gen,
     &           .false.)
        enddo
      else
        call errquit(
     &        'int_1epxvpe: should never reach transform blocked else',
     &        911, INT_ERR)
      endif
c
      end
C> @}
