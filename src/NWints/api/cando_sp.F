      logical function cando_sp(basisin,ish,jsh)
      implicit none
c $Id: cando_sp.F,v 1.18 1996-06-27 22:33:26 d3g681 Exp $
c
c this routine does not check the shell and int_init
c because it is only called by routines that do this.
c 
***********************************************************************
c-tex-% part of API Internal Routines
c-tex-\subsubsection{cando\_sp}
c-tex-This routine returns ``.true.'' to the calling integral API 
c-tex-routine the integrals over ISH and JSH can be computed using 
c-tex-the sp-rotated axis code.
***********************************************************************
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "candoP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
c-tex-\begin{verbatim}
c-tex-      logical function cando_sp(basisin,ish,jsh)
      integer basisin   ! [input] basis set handle for identified shells 
      integer ish, jsh  ! [input] lexical shell index for I and J
c-tex-\end{verbatim}
c-tex-A zero [i|j]sh means ignore this shell. this allows one function 
c-tex-for all 2e permutations and multiple basis sets
c::local
      integer basis, mytype, mycont
      integer ucont
#include "bas_ibs_sfn.fh"
c
      if (user_cando_sp) then
        cando_sp = def_cando_sp
        return
      endif
c
      cando_sp = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c.... check ish
      if (ish.gt.0) then
        ucont = (sf_ibs_cn2ucn(ish,basis))
        mytype = infbs_cont(CONT_TYPE, ucont,basis)
        mycont = infbs_cont(CONT_NPRIM,ucont,basis)
        cando_sp = cando_sp .and.
     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1)) .and.
     $         (mycont.le.25)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        ucont = (sf_ibs_cn2ucn(jsh,basis))
        mytype = infbs_cont(CONT_TYPE, ucont,basis)
        mycont = infbs_cont(CONT_NPRIM,ucont,basis)
        cando_sp = cando_sp .and.
     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1)) .and.
     $         (mycont.le.25)
      endif
c
      end
      logical function cando_nw(basisin,ish,jsh)
      implicit none
c
c a zero [i|j]sh means ignore this shell
c this allows one function for all 2e permutations
c 
c for two electron cases then there are two calls to cando_nw
c that are required
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
      integer basisin,ish,jsh
c::local
      integer basis, mytype, mygen, ucont
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
c.... check ish
      if (ish.gt.0) then
        ucont  = (sf_ibs_cn2ucn(ish,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        mygen  = infbs_cont(CONT_NGEN,ucont,basis)
        cando_nw = cando_nw .and.(mytype.gt.-1) 
        cando_nw = cando_nw .and.(mygen.eq.1)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        ucont  = (sf_ibs_cn2ucn(jsh,basis))
        mytype = infbs_cont(CONT_TYPE,ucont,basis)
        mygen  = infbs_cont(CONT_NGEN,ucont,basis)
        cando_nw = cando_nw .and.(mytype.gt.-1)
        cando_nw = cando_nw .and.(mygen.eq.1)
      endif
c
      end
      logical function cando_sp_u(basisin,ish_u,jsh_u)
      implicit none
c
c a zero [i|j]sh means ignore this shell/contraction
c this allows one function for all 2e permutations
c 
c for two electron cases then there are two calls to cando_sp_u
c that are required
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
c this routine cando_sp_u works for unique shells only
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
c::passed
      integer basisin,ish_u,jsh_u
c::local
      integer basis, mytype, mycont
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
c.... check ish_u
      if (ish_u.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ish_u,basis)
        mycont = infbs_cont(CONT_NPRIM,ish_u,basis)
        cando_sp_u = cando_sp_u .and.
     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1)) .and.
     $         (mycont.le.25)
      endif
c
c.... check jsh_u
      if (jsh_u.gt.0) then
        mytype = infbs_cont(CONT_TYPE,jsh_u,basis)
        mycont = infbs_cont(CONT_NPRIM,jsh_u,basis)
        cando_sp_u = cando_sp_u .and.
     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1)) .and.
     $         (mycont.le.25)
      endif
c
      end
      logical function cando_txs(basisin,ish,jsh)
      implicit none
      integer basisin
      integer ish,jsh
#if defined(USE_TEXAS)
      integer mytype, basis
#endif
c
#include "bas.fh"
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
*old      if(ish.eq.0.or.jsh.eq.0) cando_txs=.false.
crjh debug
      if (ish.le.0 .or. ish.gt.ncont_tot_gb(basis))
     $     call errquit('cando_txs: ish bad', ish)
      if (jsh.le.0 .or. jsh.gt.ncont_tot_gb(basis))
     $     call errquit('cando_txs: jsh bad', jsh)
crjh enddebug
      if (ish.gt.0) then
        mytype = infbs_cont(CONT_TYPE,sf_ibs_cn2ucn(ish,basis),basis)
        cando_txs = cando_txs.and. (mytype.ge.-1)
      endif
      if (jsh.gt.0) then
        mytype = infbs_cont(CONT_TYPE,sf_ibs_cn2ucn(jsh,basis),basis)
        cando_txs = cando_txs.and. (mytype.ge.-1)
      endif
#else
      cando_txs = .false.
#endif
      end
      subroutine int_set_cando_txs(user_value)
      implicit none
c
#include "candoP.fh"
c
      logical oldvalue_user, oldvalue_def
      common /cndo_txs/ oldvalue_user, oldvalue_def
      logical user_value
      
      oldvalue_user = user_cando_txs
      oldvalue_def  = def_cando_txs
      user_cando_txs = .true.
      def_cando_txs  = user_value
      if (def_cando_txs) call errquit
     &('int_set_cando_txs: logic error with user cando_texas settings',
     &    911)
      end
      subroutine int_unset_cando_txs()
      implicit none
c
#include "candoP.fh"
c
      logical oldvalue_user, oldvalue_def
      common /cndo_txs/ oldvalue_user, oldvalue_def
c
      user_cando_txs = oldvalue_user
      def_cando_txs  = oldvalue_def
c
      end
      logical function cando_nw_1e(basisin,ish,jsh)
      implicit none
c
c a zero [i|j]sh means ignore this shell
c this allows one function for all 1e permutations
c 
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
      integer basisin,ish,jsh
c::local
      integer basis, mytype, ucont
c
#include "bas_ibs_sfn.fh"
c
      cando_nw_1e = .true.
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
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
