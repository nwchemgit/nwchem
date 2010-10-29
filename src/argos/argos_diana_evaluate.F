      subroutine argos_diana_evaluate()
c
c $Id$
c
c     evaluate the defined properties
c     -------------------------------
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer itval
      integer ilt,iht,jlt,jht,ndxt,ldt
c
      filord=card(8:80)
      fildat=filord(1:index(filord,' ')-1)//'.dat'
      filord=filord(1:index(filord,' ')-1)//'.ord'
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
c     allocate global array for the properties
c     ----------------------------------------
c
      if(ntval.gt.0) then
      if(.not.ga_destroy(ga_tval))
     + call md_abort('Failed to destroy ga_tval',0)
      endif
      ntval=ndist+nangl+ntors+nimps
c
      if(ntval.gt.0) then
c
      if(.not.ga_create(mt_dbl,ntval,nfrtot,'tval',ntval,nfrdim,
     + ga_tval)) call md_abort('Failed to create ga_tval',0)
c
      itval=1
c
      call ga_distribution(ga_tval,me,ilt,iht,jlt,jht)
      call ga_access(ga_tval,ilt,iht,jlt,jht,ndxt,ldt)
c
c     evaluate distances
c     ------------------
c
      if(ndist.gt.0)
     + call argos_diana_evaldist(dbl_mb(ndx),dbl_mb(ndxt),itval)
      itval=itval+ndist
c
c     evaluate angles
c     ---------------
c
      if(nangl.gt.0)
     + call argos_diana_evalangl(dbl_mb(ndx),dbl_mb(ndxt),itval)
      itval=itval+nangl
c
c     evaluate torsions
c     -----------------
c
      if(ntors.gt.0)
     + call argos_diana_evaltors(dbl_mb(ndx),dbl_mb(ndxt),itval)
      itval=itval+ntors
c
      endif
c
      if(nord.gt.0)
     + call argos_diana_evalord(dbl_mb(ndx))
c
      call argos_diana_export()
c
c     reset
c
      ndist=0
      nangl=0
      ntors=0
      nord=0
c
      return
      end
