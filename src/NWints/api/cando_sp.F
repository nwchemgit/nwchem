c $Id$
*
C> \ingroup nwint
C> @{
C>
C> \brief Check whether integrals can be done with the sp-rotated axis
C> code
C>
C> This routine returns `.true.` to the calling integral API 
C> routine if the two electron integrals over `ISH` and `JSH` can be 
C> computed using the sp-rotated axis code.
C> 
C> \return Return .true. if the combination of shells can be handled
C> by the sp-rotated axis code, and .false. otherwise
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{cando\_sp}
c:tex-This routine returns ``.true.'' to the calling integral API 
c:tex-routine if the two electron integrals over ISH and JSH can be 
c:tex-computed using the sp-rotated axis code.
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function cando_sp(basisin,ish,jsh)
c:tex-\end{verbatim}
      implicit none
c
c this routine does not check the shell and int_init
c because it is only called by routines that do this.
c 
c
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "candoP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
#include "util.fh"
c::passed
c:tex-\begin{verbatim}
      integer basisin !< [Input] basis set handle for identified shells 
      integer ish     !< [Input] lexical shell index for I
      integer jsh     !< [Input] lexical shell index for J
c:tex-\end{verbatim}
c:tex-A zero [i|j]sh means ignore this shell. this allows one function 
c:tex-for all 2e permutations and multiple basis sets.
c::local
#if defined(USE_SP_CODE)
      integer basis
      integer i_ucont
      integer j_ucont
#endif
*_old_code:      integer mytype, mycont
*_old_code:      integer ucont
c
      logical cando_sp_u
      external cando_sp_u
c
#include "bas_ibs_sfn.fh"
c
#if defined(USE_SP_CODE)
      basis = basisin + BASIS_HANDLE_OFFSET
c.. get ish unique shell index
      i_ucont = (sf_ibs_cn2ucn(ish,basis))
c.. get jsh unique shell index
      j_ucont = (sf_ibs_cn2ucn(jsh,basis))
c
      cando_sp = cando_sp_u(basisin, i_ucont, j_ucont)
#else
      cando_sp = .false.
#endif
c
      return
      end
*_old_code:c
*_old_code:c.... check ish
*_old_code:      if (ish.gt.0) then
*_old_code:        ucont = 
*_old_code:        mytype = infbs_cont(CONT_TYPE, ucont,basis)
*_old_code:        mycont = infbs_cont(CONT_NPRIM,ucont,basis)
*_old_code:        cando_sp = cando_sp .and.
*_old_code:     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1)) .and.
*_old_code:     $         (mycont.le.25)
*_old_code:      endif
*_old_code:c
*_old_code:c.... check jsh
*_old_code:      if (jsh.gt.0) then
*_old_code:        ucont = (sf_ibs_cn2ucn(jsh,basis))
*_old_code:        mytype = infbs_cont(CONT_TYPE, ucont,basis)
*_old_code:        mycont = infbs_cont(CONT_NPRIM,ucont,basis)
*_old_code:        cando_sp = cando_sp .and.
*_old_code:     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1)) .and.
*_old_code:     $         (mycont.le.25)
*_old_code:      endif
*_old_code:c
*_old_code:      end
***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsection{cando\_nw}
c:tex-This routine returns ``.true.'' to the calling integral API 
c:tex-routine if the two electron integrals over ISH and JSH can be 
c:tex-computed using the McMurchie-Davidson code.
***********************************************************************
c:tex-{\it Syntax:}
C> \brief Checks whether the 1-electron shell pair integrals can be
C> evaluated using the McMurchie-Davidson code
C>
C> Considering the specified basis set and the shells `ISH` and `JSH`
C> checks whether those integrals can be evaluated using the
C> McMurchie-Davidson code.
C>
C> \return Return .true. if the integrals can be evaluated using the 
C> McMurchie-Davidson code, and .false. otherwise
c
      logical function cando_nw(basisin,ish,jsh)
      implicit none
c 
c:tex-For two electron cases then there are two calls to cando_nw
c:tex-that are required with each shell pair and cooresponding 
c:tex-basis set handle (bra or ket). 
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
c
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
      integer basisin !< [Input] the basis set handle for identified
                      !< shells 
      integer ish     !< [Input] the lexical shell index for I
      integer jsh     !< [Input] the lexical shell index for J
c:tex-A zero [i|j]sh means ignore this shell. this allows one function 
c:tex-for all 2e permutations and multiple basis sets.
c::local
      integer basis, mytype, ucont
c
#include "bas_ibs_sfn.fh"
c
      if (user_cando_nw) then
        cando_nw = def_cando_nw
        return
      endif
c
      cando_nw = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c check shell range
      if (ish.lt.0 .or. ish.gt.ncont_tot_gb(basis))
     $     call errquit('cando_nw: bad ish =', ish, UNKNOWN_ERR)
      if (jsh.lt.0 .or. jsh.gt.ncont_tot_gb(basis))
     $     call errquit('cando_nw: bad jsh =', jsh, UNKNOWN_ERR)
c.... check ish
      if (ish.gt.0) then
        ucont  = (sf_ibs_cn2ucn(ish,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        cando_nw = cando_nw .and.(mytype.gt.-1) 
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        ucont  = (sf_ibs_cn2ucn(jsh,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        cando_nw = cando_nw .and.(mytype.gt.-1)
      endif
c
      end
***********************************************************************
C> \brief Checks whether the unique shell combination can be handled by
C> the sp-rotated axis code
C>
C> This routine returns `.true.` to the calling integral API 
C> routine if the two electron integrals over unique shells `ISH`
C> and `JSH` can be computed using the sp-rotated axis code.
C>
C> \return Return `.true.` if the combination of unique shells can be
C> handled by the sp-rotated axis code, and `.false.` otherwise.
c:tex-% part of API Internal Routines
c:tex-\subsection{cando\_sp\_u}
c:tex-This routine returns ``.true.'' to the calling integral API 
c:tex-routine if the two electron integrals over unique shells ISH 
c:tex-and JSH can be computed using the sp-rotated axis code.
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function cando_sp_u(basisin,ish_u,jsh_u)
c:tex-\end{verbatim}
      implicit none
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
c::passed
c:tex-\begin{verbatim}
      integer basisin !< [Input] basis set handle for identified shells 
      integer ish_u   !< [Input] unique shell index for I
      integer jsh_u   !< [Input] unique shell index for J
c:tex-\end{verbatim}
c:tex-A zero [i|j]sh means ignore this shell. this allows one function 
c:tex-for all 2e permutations and multiple basis sets.  This routine 
c:tex-{\tt cando\_sp\_u} works for unique shells only.
c::local
#if defined(USE_SP_CODE)
      integer basis, mytype, mycont, mygen
c
      if (user_cando_sp) then
        cando_sp_u = def_cando_sp
        return
      endif
c
      cando_sp_u = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c check shell range
      if (ish_u.lt.0 .or. ish_u.gt.infbs_head(Head_NCONT,basis))
     $     call errquit('cando_sp_u: bad ish_u =', ish_u, UNKNOWN_ERR)
      if (jsh_u.lt.0 .or. jsh_u.gt.infbs_head(Head_NCONT,basis))
     $     call errquit('cando_sp_u: bad jsh_u =', jsh_u, UNKNOWN_ERR)
c.... check ish_u
      if (ish_u.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ish_u,basis)
        mycont = infbs_cont(CONT_NPRIM,ish_u,basis)
        mygen  = infbs_cont(CONT_NGEN,ish_u,basis)
c
c     SP shell or singly contracted S or P, nprim <=25
c
        cando_sp_u = cando_sp_u .and.
     &       ((mytype.eq.-1).or.
     $       (((mytype.eq.0).or.(mytype.eq.1)).and.mygen.eq.1)) .and.
     $         (mycont.le.25)
      endif
c
c.... check jsh_u
      if (jsh_u.gt.0) then
        mytype = infbs_cont(CONT_TYPE,jsh_u,basis)
        mycont = infbs_cont(CONT_NPRIM,jsh_u,basis)
        mygen  = infbs_cont(CONT_NGEN,jsh_u,basis)
        cando_sp_u = cando_sp_u .and.
     &       ((mytype.eq.-1).or.
     $       (((mytype.eq.0).or.(mytype.eq.1)).and.mygen.eq.1)) .and.
     $         (mycont.le.25)
      endif
#else
      cando_sp_u = .false.
#endif
c
      end
***********************************************************************
C> \brief Checks whether the unique shell combination can be handled by
C> the texas axis code
C>
C> This routine returns `.true.` to the calling integral API 
C> routine if the two electron integrals over unique shells `ISH`
C> and `JSH` can be computed using the texas axis code.
C>
C> \return Return `.true.` if the combination of unique shells can be
C> handled by the texas axis code, and `.false.` otherwise.
c:tex-% part of API Internal Routines
c:tex-\subsection{cando\_txs}
c:tex-This routine returns ``.true.'' to the calling integral API 
c:tex-routine if the two electron integrals over ISH and JSH can be 
c:tex-computed using the texas integral code.
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function cando_txs(basisin,ish,jsh)
c:tex-\end{verbatim}
      implicit none
c:tex-\begin{verbatim}
      integer basisin !< [Input] basis set handle for identified shells 
      integer ish     !< [Input] lexical shell index for I
      integer jsh     !< [Input] lexical shell index for J
c:tex-\end{verbatim}
c:tex-A zero [i|j]sh means ignore this shell. This allows one function 
c:tex-for all 2e permutations.
#define USE_TEXAS
#if defined(USE_TEXAS)
      integer mytype, basis
#endif
c
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c
#include "bas_ibs_sfn.fh"
#if defined(USE_TEXAS)
      if (user_cando_txs) then
        cando_txs = def_cando_txs
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c
      cando_txs = .true.
c
c ish/jsh != 0 for texas right now.
c   this will change once texas can do 2 and 3 center integrals
c
c once available change .le. to .lt.
c
      if (ish.lt.0 .or. ish.gt.ncont_tot_gb(basis))
     $     call errquit('cando_txs: ish bad', ish, UNKNOWN_ERR)
      if (jsh.lt.0 .or. jsh.gt.ncont_tot_gb(basis))
     $     call errquit('cando_txs: jsh bad', jsh, UNKNOWN_ERR)
c
c check for spherical false until spherical debuged in texas
*      cando_txs = cando_txs .and. (.not.bas_spherical(basis))
      if (ish.gt.0) then
        mytype = infbs_cont(CONT_TYPE,sf_ibs_cn2ucn(ish,basis),basis)
        cando_txs = cando_txs.and. (mytype.ge.-1)
     $             .and. (mytype.le.5)
      endif
      if (jsh.gt.0) then
        mytype = infbs_cont(CONT_TYPE,sf_ibs_cn2ucn(jsh,basis),basis)
        cando_txs = cando_txs.and. (mytype.ge.-1) 
     $             .and. (mytype.le.5)
      endif
#else
      cando_txs = .false.
#endif
      end
***********************************************************************
C> \brief Checks whether the shell combination can be handled by
C> the McMurchie-Davidson 1-electron integral code
C>
C> This routine returns `.true.` to the calling integral API 
C> routine if the one electron integrals over `ISH` and `JSH` can be 
C> computed using the McMurchie-Davidson code.
C>
C> \return Return `.true.` if the combination of shells can be
C> handled by the 1-electron McMurchi-Davidson code, and `.false.`
C> otherwise.
c:tex-% part of API Internal Routines
c:tex-\subsection{cando\_nw_\1e}
c:tex-This routine returns ``.true.'' to the calling integral API 
c:tex-routine if the one electron integrals over ISH and JSH can be 
c:tex-computed using the McMurchie-Davidson code.
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function cando_nw_1e(basisin,ish,jsh)
c:tex-\end{verbatim}
      implicit none
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
c:tex-\begin{verbatim}
      integer basisin !< [Input] basis set handle for identified shells 
      integer ish     !< [Input] lexical shell index for I
      integer jsh     !< [Input] lexical shell index for J
c:tex-\end{verbatim}
c:tex-A zero [i|j]sh means ignore this shell. this allows one function 
c:tex-for all 1e permutations and multiple basis sets.
c::local
      integer basis, mytype, ucont
c
#include "bas_ibs_sfn.fh"
c
      cando_nw_1e = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c check shell range
      if (ish.lt.0 .or. ish.gt.ncont_tot_gb(basis))
     $     call errquit('cando_nw_1e: bad ish =', ish, UNKNOWN_ERR)
      if (jsh.lt.0 .or. jsh.gt.ncont_tot_gb(basis))
     $     call errquit('cando_nw_1e: bad jsh =', jsh, UNKNOWN_ERR)
c.... check ish
      if (ish.gt.0) then
        ucont  = (sf_ibs_cn2ucn(ish,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        cando_nw_1e = cando_nw_1e .and.(mytype.ge.-1)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        ucont  = (sf_ibs_cn2ucn(jsh,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        cando_nw_1e = cando_nw_1e .and.(mytype.ge.-1)
      endif
c
      end
***********************************************************************
C> \brief Checks whether the 1-electron shell pair integrals can be
C> evaluated using the HONDO code
C>
C> Considering the specified basis set and the shells `ISH` and `JSH`
C> checks whether those integrals can be evaluated using the HONDO
C> code.
C>
C> \return Return `.true.` if the integrals can be evaluated using the 
C> HONDO code, and `.false.` otherwise
c
      logical function cando_hnd_1e(basisin,ish,jsh)
      implicit none
c
c this routine does not check the shell and int_init
c because it is only called by routines that do this.
c
#include "apiP.fh"
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
      integer basisin  !< [Input] the basis set handle for identified
                       !< shells
      integer ish      !< [Input] the lexical shell index for I
      integer jsh      !< [Input] the lexical shell index for J
c:tex-A zero [i|j]sh means ignore this shell. this allows one function
c:tex-for all 1e permutations and multiple basis sets.
c::local
      integer basis, mytype, ucont, ngen
c
c
#include "bas_ibs_sfn.fh"
c
      if (user_cando_hnd) then
        cando_hnd_1e = def_cando_hnd
        return
      endif
c
      cando_hnd_1e = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c check shell range
c also check for general contractions which the code does not
c currently handle
c
      if (ish.lt.0 .or. ish.gt.ncont_tot_gb(basis))
     $     call errquit('cando_hnd_1e: bad ish =', ish, UNKNOWN_ERR)
      if (jsh.lt.0 .or. jsh.gt.ncont_tot_gb(basis))
     $     call errquit('cando_hnd_1e: bad jsh =', jsh, UNKNOWN_ERR)
c.... check ish
      if (ish.gt.0) then
        ucont  = (sf_ibs_cn2ucn(ish,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        ngen = infbs_cont(CONT_NGEN,ucont,basis)
        cando_hnd_1e = cando_hnd_1e .and. (mytype.ge.0) .and.
     $     (mytype.le.4) .and. (ngen.eq.1)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        ucont  = (sf_ibs_cn2ucn(jsh,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        ngen = infbs_cont(CONT_NGEN,ucont,basis)
        cando_hnd_1e = cando_hnd_1e .and. (mytype.ge.0) .and.
     $     (mytype.le.4) .and. (ngen.eq.1)
      endif
c
c for now the code does not work with ECPs
c this may be added later.
c
      if (any_ecp) cando_hnd_1e = .false.
c
      end
***********************************************************************
C> \brief Turn the Texas integral code off
C>
C> This routine stores any current user setting for the texas 
C> integral code and then sets the cando_txs flag on the RTDB to
C> false.  This is for applications that have problems using the
C> texas integrals.  Right now this is the dft module.  This
C> routine should be called right before the call to `int_init`.
c:tex-% part of API application interface
c:tex-\subsection{int\_app\_set\_no\_texas}
c:tex-This routine stores any current user setting for the texas 
c:tex-integral code and then sets the cando_txs flag on the RTDB to
c:tex-false.  This is for applications that have problems using the
c:tex-texas integrals.  Right now this is the dft module.  This
c:tex-routine should be called right before the call to ``int_init''
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_app_set_no_texas(rtdb)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "candoP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer rtdb !< [Input] Run Time Data Base (RTDB) handle
c:tex-\end{verbatim}
c
      logical user_status
c
      if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,user_status)) then
        app_stored_txs = user_status
      else
        app_stored_txs = .true.
      endif
      if (.not.rtdb_put(rtdb,'int:cando_txs',mt_log,1,.false.))
     &    call errquit('int_app_set_no_texas:rtdb_put failed',911,
     &       RTDB_ERR)
      rtdbIused = rtdb
      end
***********************************************************************
C> \brief Turn the Texas integral code on
C>
C> This routine resets the current user setting for the texas 
C> integral code and that was stored by `int_app_set_no_texas` by 
C> setting the cando_txs flag on the RTDB to the state prior to 
C> turning texas off.  This is for applications that have problems 
C> using the texas integrals.  Right now this is the dft module.  
C> This routine should be called right after the call to 
C> `int_terminate`.
c:tex-% part of API application interface
c:tex-\subsection{int\_app\_unset\_no\_texas}
c:tex-This routine resets the current user setting for the texas 
c:tex-integral code and that was stored by int\_app\_set\_no\_texas by 
c:tex-setting the cando_txs flag on the RTDB to the state prior to 
c:tex-turning texas off.  This is for applications that have problems 
c:tex-using the texas integrals.  Right now this is the dft module.  
c:tex-Note: This routine should be called right after the call to 
c:tex-``int_terminate''
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_app_unset_no_texas(rtdb)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "candoP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer rtdb !< [Input] Run Time Data Base (RTDB) handle
c:tex-\end{verbatim}
c
      if (app_stored_txs) then
        if(.not.rtdb_delete(rtdb,'int:cando_txs')) 
     &      call errquit
     &      ('int_app_unset_no_texas:rtdb_delete failed',911, RTDB_ERR)
      else
        if (.not.rtdb_put(rtdb,'int:cando_txs',mt_log,1,app_stored_txs))
     &      call errquit('int_app_unset_no_texas:rtdb_put failed',911, 
     &      RTDB_ERR)
      endif
      end
*_old_code:      subroutine int_set_cando_txs(user_value)
*_old_code:      implicit none
*_old_code:c
*_old_code:#include "candoP.fh"
*_old_code:c
*_old_code:      logical oldvalue_user, oldvalue_def
*_old_code:      common /cndo_txs/ oldvalue_user, oldvalue_def
*_old_code:      logical user_value
*_old_code:      
*_old_code:      oldvalue_user = user_cando_txs
*_old_code:      oldvalue_def  = def_cando_txs
*_old_code:      user_cando_txs = .true.
*_old_code:      def_cando_txs  = user_value
*_old_code:      if (def_cando_txs) call errquit
*_old_code:     &('int_set_cando_txs: logic error with user cando_texas settings',
*_old_code:     &    911)
*_old_code:      end
*_old_code:      subroutine int_unset_cando_txs()
*_old_code:      implicit none
*_old_code:c
*_old_code:#include "candoP.fh"
*_old_code:c
*_old_code:      logical oldvalue_user, oldvalue_def
*_old_code:      common /cndo_txs/ oldvalue_user, oldvalue_def
*_old_code:c
*_old_code:      user_cando_txs = oldvalue_user
*_old_code:      def_cando_txs  = oldvalue_def
*_old_code:c
*_old_code:      end
***********************************************************************
C> \brief Turns the SP-rotated axis integral code off
C>
C> This routine stores any current user setting for the spint 
C> integral code and then sets the cando_sp flag on the RTDB to
C> false.  This is for applications that have problems using the
C> sp integrals.  This routine should be called right before the 
C> call to `int_init`.
c:tex-% part of API application interface
c:tex-\subsection{int\_app\_set\_no\_spint}
c:tex-This routine stores any current user setting for the spint 
c:tex-integral code and then sets the cando_sp flag on the RTDB to
c:tex-false.  This is for applications that have problems using the
c:tex-sp integrals.  This routine should be called right before the 
c:tex-call to ``int_init''
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_app_set_no_spint(rtdb)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "candoP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer rtdb !< [Input] Run Time Data Base (RTDB) handle
c:tex-\end{verbatim}
c
      logical user_status
c
      if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,user_status)) then
        app_stored_spint = user_status
      else
        app_stored_spint = .true.
      endif
      if (.not.rtdb_put(rtdb,'int:cando_sp',mt_log,1,.false.))
     &    call errquit('int_app_set_no_spint:rtdb_put failed',911,
     &       RTDB_ERR)
      rtdbIused = rtdb
      end
***********************************************************************
C> \brief Turns the SP-rotated axis integral code on
C>
C> This routine resets the current user setting for the sp
C> integral code and that was stored by int_app_set_no_spint by 
C> setting the cando_sp flag on the RTDB to the state prior to 
C> turning sp integrals off.  This is for applications that have 
C> problems using the sp integrals.  
C> This routine should be called right after the call to 
C> `int_terminate`.
c:tex-% part of API application interface
c:tex-\subsection{int\_app\_unset\_no\_spint}
c:tex-This routine resets the current user setting for the sp
c:tex-integral code and that was stored by int\_app\_set\_no\_spint by 
c:tex-setting the cando_sp flag on the RTDB to the state prior to 
c:tex-turning sp integrals off.  This is for applications that have 
c:tex-problems using the sp integrals.  
c:tex-Note: This routine should be called right after the call to 
c:tex-``int_terminate''
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_app_unset_no_spint(rtdb)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "candoP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer rtdb !< [Input] Run Time Data Base (RTDB) handle
c:tex-\end{verbatim}
c
      if (app_stored_spint) then
        if(.not.rtdb_delete(rtdb,'int:cando_sp')) 
     &      call errquit
     &      ('int_app_unset_no_spint:rtdb_delete failed',911, RTDB_ERR)
      else
        if (.not.rtdb_put(rtdb,'int:cando_sp',
     &      mt_log,1,app_stored_spint))
     &      call errquit('int_app_unset_no_spint:rtdb_put failed',911, 
     &      RTDB_ERR)
      endif
      end
***********************************************************************
C> \brief Turn the Hondo integral code off
C>
C> This routine stores any current user setting for the hondo
C> integral code and then sets the cando_hnd flag on the RTDB to
C> false.  This is for applications that have problems using the
C> hondo integrals.  This routine should be called right before the
C> call to `int_init`.
c:tex-% part of API application interface
c:tex-\subsection{int\_app\_set\_no\_hnd}
c:tex-This routine stores any current user setting for the hondo
c:tex-integral code and then sets the cando_hnd flag on the RTDB to
c:tex-false.  This is for applications that have problems using the
c:tex-hondo integrals.  This routine should be called right before the
c:tex-call to ``int_init''
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_app_set_no_hnd(rtdb)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "candoP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer rtdb !< [Input] Run Time Data Base (RTDB) handle
c:tex-\end{verbatim}
c
      logical user_status
c
      if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,user_status)) then
        app_stored_hnd = user_status
      else
        app_stored_hnd = .true.
      endif
      if (.not.rtdb_put(rtdb,'int:cando_hnd',mt_log,1,.false.))
     &    call errquit('int_app_set_no_hnd:rtdb_put failed',911,
     &       RTDB_ERR)
      rtdbIused = rtdb
      end
***********************************************************************
C> \brief Turn the Hondo integral code on
C>
C> This routine resets the current user setting for the hondo
C> integral code that was stored by `int_app_set_no_hnd` by
C> setting the cando_hnd flag on the RTDB to the state prior to
C> turning hondo integrals off.  This is for applications that have
C> problems using the hondo integrals.
C> This routine should be called right after the call to
C> `int_terminate`.
c:tex-% part of API application interface
c:tex-\subsection{int\_app\_unset\_no\_hnd}
c:tex-This routine resets the current user setting for the hondo
c:tex-integral code that was stored by int\_app\_set\_no\_hnd by
c:tex-setting the cando_hnd flag on the RTDB to the state prior to
c:tex-turning hondo integrals off.  This is for applications that have
c:tex-problems using the hondo integrals.
c:tex-Note: This routine should be called right after the call to
c:tex-``int_terminate''
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      subroutine int_app_unset_no_hnd(rtdb)
c:tex-\end{verbatim}
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "candoP.fh"
#include "util.fh"
c:tex-\begin{verbatim}
      integer rtdb !< [Input] Run Time Data Base (RTDB) handle
c:tex-\end{verbatim}
c
      if (app_stored_hnd) then
        if(.not.rtdb_delete(rtdb,'int:cando_hnd'))
     &      call errquit
     &      ('int_app_unset_no_hnd:rtdb_delete failed',911, RTDB_ERR)
      else
        if (.not.rtdb_put(rtdb,'int:cando_hnd',
     &      mt_log,1,app_stored_hnd))
     &      call errquit('int_app_unset_no_hnd:rtdb_put failed',911, 
     &      RTDB_ERR)
      endif
      end
***********************************************************************
C> \brief Check whether the Hondo code can calculate 1-electron 
C> second derivative integrals for the given shells
C>
C> This routine returns ``.true.'' to the calling integral API routine
C> if the second derivatives of the one electron integrals over `ISH`
C> and `JSH` can be computed using the HONDO code.
C>
C> \return Return `.true.` if the Hondo code can calculate the
C> integrals, and `.false.` otherwise
c:tex-% part of API Internal Routines
c:tex-\subsection{cando\_hnd_\1edd}
c:tex-This routine returns ``.true.'' to the calling integral API
c:tex-routine if the second derivatives of the one electron integrals over ISH and JSH can be
c:tex-computed using the HONDO code.
***********************************************************************
c:tex-{\it Syntax:}
c:tex-\begin{verbatim}
      logical function cando_hnd_1edd(basisin,ish,jsh)
c:tex-\end{verbatim}
      implicit none
c
c this routine does not check the shell and int_init
c because it is only called by routines that do this.
c
#include "errquit.fh"
#include "apiP.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
c:tex-\begin{verbatim}
      integer basisin !< [Input] basis set handle for identified shells
      integer ish     !< [Input] lexical shell index for I
      integer jsh     !< [Input] lexical shell index for J
c:tex-\end{verbatim}
c:tex-A zero [i|j]sh means ignore this shell. this allows one function
c:tex-for all 1e permutations and multiple basis sets.
c::local
      integer basis, mytype, ucont, ngen
c
c
#include "bas_ibs_sfn.fh"
c
      if (user_cando_hnd) then
        cando_hnd_1edd = def_cando_hnd
        return
      endif
c
      cando_hnd_1edd = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c check shell range
c also check for general contractions which the code does not
c currently handle
c
      if (ish.lt.0 .or. ish.gt.ncont_tot_gb(basis))
     $     call errquit('cando_hnd_1edd: bad ish =', ish, BASIS_ERR)
      if (jsh.lt.0 .or. jsh.gt.ncont_tot_gb(basis))
     $     call errquit('cando_hnd_1edd: bad jsh =', jsh, BASIS_ERR)
c.... check ish
      if (ish.gt.0) then
        ucont  = (sf_ibs_cn2ucn(ish,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        ngen = infbs_cont(CONT_NGEN,ucont,basis)
        cando_hnd_1edd = cando_hnd_1edd .and. (mytype.ge.0) .and.
     $     (mytype.le.6) .and. (ngen.eq.1)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        ucont  = (sf_ibs_cn2ucn(jsh,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        ngen = infbs_cont(CONT_NGEN,ucont,basis)
        cando_hnd_1edd = cando_hnd_1edd .and. (mytype.ge.0) .and.
     $     (mytype.le.6) .and. (ngen.eq.1)
      endif
c
      end
C> \return Return .true. if the integrals can be evaluated using the 
C> Simint code, and .false. otherwise
c
      logical function cando_sim(basisin,ish,jsh)
      implicit none
c 
c:tex-For two electron cases then there are two calls to cando_nw
c:tex-that are required with each shell pair and cooresponding 
c:tex-basis set handle (bra or ket). 
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
c
#include "bas.fh"
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
      integer basisin !< [Input] the basis set handle for identified
                      !< shells 
      integer ish     !< [Input] the lexical shell index for I
      integer jsh     !< [Input] the lexical shell index for J
c:tex-A zero [i|j]sh means ignore this shell. this allows one function 
c:tex-for all 2e permutations and multiple basis sets.
c::local
      integer basis, mytype, ucont
      logical nwcsim_initdone
      external nwcsim_initdone
c
#include "bas_ibs_sfn.fh"
#ifdef USE_SIMINT
c
      cando_sim = .false.
      if(.not.nwcsim_initdone()) return
      if (user_cando_sim) then
        cando_sim = def_cando_sim
        return
      endif
c
      cando_sim = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c check shell range
      if (ish.lt.0 .or. ish.gt.ncont_tot_gb(basis))
     $     call errquit('cando_sim: bad ish =', ish, UNKNOWN_ERR)
      if (jsh.lt.0 .or. jsh.gt.ncont_tot_gb(basis))
     $     call errquit('cando_sim: bad jsh =', jsh, UNKNOWN_ERR)
c.... check ish
      if (ish.gt.0) then
        ucont  = (sf_ibs_cn2ucn(ish,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        cando_sim = cando_sim .and.(mytype.gt.-1) 
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        ucont  = (sf_ibs_cn2ucn(jsh,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        cando_sim = cando_sim .and.(mytype.gt.-1)
      endif
c
c     add check if ish and jsh are stored in c_simint_shell
c
c
#else
      cando_sim = .false.
#endif
      end
C> @}
