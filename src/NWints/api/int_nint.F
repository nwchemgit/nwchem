      integer function int_nint(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
c $Id: int_nint.F,v 1.1 1994-06-08 22:39:39 d3e129 Exp $
      implicit none
c
c calculates the number of integrals for a given shell/contraction 
c grouping if an input shell is zero then it ignores this shell.  This 
c routine will work for all integral routines.  
c
#include "bas.fh"
c::passed
      integer ibasin   ! [input] basis set handle for icnt
      integer icnt     ! [input] contraction index (e.g., ish)
      integer jbasin   ! [input] basis set handle for jcnt
      integer jcnt     ! [input] contraction index (e.g., jsh)
      integer kbasin   ! [input] basis set handle for kcnt
      integer kcnt     ! [input] contraction index (e.g., ksh)
      integer lbasin   ! [input] basis set handle for lcnt
      integer lcnt     ! [input] contraction index (e.g., lsh)
c::local
      integer hi,lo,nbf
c
      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_nint = 0
        return
      endif
c
      int_nint = 1
c
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_cn2bfr(ibasin,icnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_cn2bfr(jbasin,jcnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_cn2bfr(kbasin,kcnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_cn2bfr(lbasin,lcnt,lo,hi))
     &         call errquit('int_nint: bas_cn2bfr failure',911)
        nbf = hi - lo + 1
        int_nint = int_nint*nbf
      endif
      end

