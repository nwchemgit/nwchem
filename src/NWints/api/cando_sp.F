      logical function cando_sp(basisin,ish,jsh)
      implicit none
c
c a zero [i|j]sh means ignore this shell
c this allows one function for all 1e and 2e permutations
c 
c for two electron cases then there are two calles to cando_sp
c that are required
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
c::passed
      integer basisin,ish,jsh
c::local
      integer basis, mytype
c
      cando_sp = bas_check_handle(basisin,'cando_sp')
      if (.not.cando_sp) then
        write(6,*)' bad basis handle '
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c.... check ish
      if (ish.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(ish,basis),basis)
        cando_sp = cando_sp .and.
     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1))
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(jsh,basis),basis)
        cando_sp = cando_sp .and.
     &         ((mytype.eq.-1).or.(mytype.eq.0).or.(mytype.eq.1))
      endif
c
#if defined(DEBUGNW)
      if(cando_sp) cando_sp = .false.
#endif
      end
      logical function cando_nw(basisin,ish,jsh)
      implicit none
c
c a zero [i|j]sh means ignore this shell
c this allows one function for all 1e and 2e permutations
c 
c for two electron cases then there are two calles to cando_nw
c that are required
c
c this routine does not check the shell and int_init
c because it only called by routines that do this.
c
c this routine will disappear after the api has been modified 
c to handle sp functions using the NWints code
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
c::passed
      integer basisin,ish,jsh
c::local
      integer basis, mytype
c
      cando_nw = bas_check_handle(basisin,'cando_nw')
      if (.not.cando_nw) then
        write(6,*)' bad basis handle '
        return
      endif
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c.... check ish
      if (ish.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(ish,basis),basis)
        cando_nw = cando_nw .and.(mytype.ne.-1)
      endif
c
c.... check jsh
      if (jsh.gt.0) then
        mytype = infbs_cont(CONT_TYPE,ibs_cn2ucn(jsh,basis),basis)
        cando_nw = cando_nw .and.(mytype.ne.-1)
      endif
c
      end
