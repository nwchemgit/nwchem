***********************************************************************
c:tex-% part of API Internal Routines
c:tex-\subsubsection{cando\_sp}
c:tex-This routine returns ``.true.'' to the calling integral API 
c:tex-routine the integrals over ISH and JSH can be computed using 
c:tex-the sp-rotated axis code.
***********************************************************************
c:tex-\begin{verbatim}
      logical function cando_sp(basisin,ish,jsh)
c:tex-\end{verbatim}
      implicit none
c $Id$
c
c this routine does not check the shell and int_init
c because it is only called by routines that do this.
c 
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "candoP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_ibs_dec.fh"
c::passed
c:tex-\begin{verbatim}
      integer basisin   ! [input] basis set handle for identified shells 
      integer ish, jsh  ! [input] lexical shell index for I and J
c:tex-\end{verbatim}
c:tex-A zero [i|j]sh means ignore this shell. this allows one function 
c:tex-for all shell permutations and multiple basis sets
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
      basis = basisin + BASIS_HANDLE_OFFSET  !:tex-A basis set handle computed by BASIS_HANDLE_OFFSET
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

