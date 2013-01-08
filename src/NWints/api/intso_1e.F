c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Compute the 1-electron spin-orbit potential integrals
C>
c:tex-% this is part of the API Standard Integral routines.
c:tex-\subsection{intso\_1e}
c:tex-This routine computes the 1 electron spin orbit potential integrals,
c:tex-
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine intso_1e(i_basis,ish,j_basis,jsh,lscr,scr,lso,SO)
c:tex-\end{verbatim}
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
      logical cando_nw_1e
      logical cando_nw
      logical int_chk_init
      logical int_chk_sh
      external int_chk_init
      external int_chk_sh
      external cando_nw_1e
      external cando_nw
c::passed
c:tex-\begin{verbatim}
      integer i_basis !< [Input] basis set handle for ish
      integer ish     !< [Input] i shell/contraction
      integer j_basis !< [Input] basis set handle for jsh
      integer jsh     !< [Input] j shell/contraction
      integer lscr    !< [Input] length of scratch array
      double precision scr(lscr) !< [Scratch] scratch array
      integer lso                !< [Input] length of SO buffer
      double precision SO(lso)   !< [Output] spin orbit potential 
                                 !< integral buffer
c:tex-\end{verbatim}
c::local
      logical shells_ok
      integer i_geom, j_geom, ibas, jbas, ucont
      integer Li, i_prim, i_gen, i_iexp, i_icfp, i_cent
      integer Lj, j_prim, j_gen, j_iexp, j_icfp, j_cent
c
      logical any_spherical
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
      if (.not.int_chk_init('int_1eso'))
     &       call errquit('int_1eso: int_init was not called' ,0,
     &           INT_ERR)
c
      shells_ok = int_chk_sh(i_basis,ish)
      shells_ok = shells_ok .and. int_chk_sh(j_basis,jsh)
      if (.not.shells_ok)
     &       call errquit('int_1eso: invalid contraction/shell',0,
     &            BASIS_ERR)
c
      call int_nospshell_check(i_basis,'intso_1e:i_basis')
      call int_nospshell_check(j_basis,'intso_1e:j_basis')
*
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
     &      'int_1eso: WARNING: possible geometry inconsistency'
        write(luout,*)'i_basis geometry handle:',i_geom
        write(luout,*)'j_basis geometry handle:',j_geom
        WarnP = 1
      endif
c
      if (cando_nw_1e(i_basis,ish,0).and.cando_nw_1e(j_basis,jsh,0))
     &    then
        call intso_hf1(
     &         coords(1,i_cent,i_geom),dbl_mb(mb_exndcf(i_iexp,ibas)),
     &         dbl_mb(mb_exndcf(i_icfp,ibas)),
     &         i_prim, i_gen, Li, i_cent,
     &         coords(1,j_cent,j_geom),dbl_mb(mb_exndcf(j_iexp,jbas)),
     &         dbl_mb(mb_exndcf(j_icfp,jbas)),
     &         j_prim, j_gen, Lj, j_cent,
     &         SO, lso, scr,lscr,.false.) ! false>dryrun
      else
        call errquit('int_1eso: could not do sp or nw integrals',0,
     &          INT_ERR)
      endif
c
*     SO now has the cartesian integral block  (jlo:jhi,ilo:ihi,3)
*      
      any_spherical = bas_spherical(ibas).or.bas_spherical(jbas)
      if (.not.any_spherical) return
c
c ... reset general contractions for sp shells to 1 since they are handled
c     as a block of 4.
c
      if (li.eq.-1) i_gen = 1
      if (lj.eq.-1) j_gen = 1
      call spcart_2cBtran(SO,scr,lscr,
     &    int_nbf_x(Li),int_nbf_s(Li),Li,i_gen,bas_spherical(ibas),
     &    int_nbf_x(Lj),int_nbf_s(Lj),Lj,j_gen,bas_spherical(jbas),
     &    3,.false.)
c
      end
C> @}
