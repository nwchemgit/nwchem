      logical function cando_sp(basisin,ish,jsh)
      implicit none
c $Id: cando_sp.F,v 1.12 1996-01-05 01:59:49 d3e129 Exp $
c
c a zero [i|j]sh means ignore this shell
c this allows one function for all 1e and 2e permutations
c 
c for two electron cases then there are two calls to cando_sp
c that are required
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "candoP.fh"
#include "geobasmapP.fh"
c::passed
      integer basisin,ish,jsh
c::local
      integer basis, mytype, mycont
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
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(ish,basis),basis)
        mycont = infbs_cont(CONT_NPRIM,ibs_cn2ucn(ish,basis),basis)
        cando_sp = cando_sp .and.
     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1)) .and.
     $         (mycont.le.25)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(jsh,basis),basis)
        mycont = infbs_cont(CONT_NPRIM,ibs_cn2ucn(jsh,basis),basis)
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
c this allows one function for all 1e and 2e permutations
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
c::passed
      integer basisin,ish,jsh
c::local
      integer basis, mytype
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
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(ish,basis),basis)
        cando_nw = cando_nw .and.(mytype.gt.-1)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(jsh,basis),basis)
        cando_nw = cando_nw .and.(mytype.gt.-1)
      endif
c
      end
      logical function cando_sp_u(basisin,ish_u,jsh_u)
      implicit none
c
c a zero [i|j]sh means ignore this shell/contraction
c this allows one function for all 1e and 2e permutations
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
      integer mytype, basis
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "candoP.fh"
c
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
      if (ish.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(ish,basis),basis)
        cando_txs = cando_txs.and. (mytype.ge.-1)
      endif
      if (jsh.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(jsh,basis),basis)
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
