      integer function int_nint(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
c $Id: int_nint.F,v 1.3 1996-07-18 23:58:27 d3e129 Exp $
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
      integer function int_unint(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
      implicit none
c
c calculates the number of integrals for a given shell/contraction 
c grouping if an input shell is zero then it ignores this shell.  This 
c routine will work for all integral routines.  
c
c::functions
      integer nbf_from_ucont
      external nbf_from_ucont
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
      integer nbf
c
      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_unint = 0
        return
      endif
c
      int_unint = 1
c
c..icnt/ish
      if (icnt.gt.0) then
        nbf = nbf_from_ucont(icnt,ibasin)
        int_unint = int_unint*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        nbf = nbf_from_ucont(jcnt,jbasin)
        int_unint = int_unint*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        nbf = nbf_from_ucont(kcnt,kbasin)
        int_unint = int_unint*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        nbf = nbf_from_ucont(lcnt,lbasin)
        int_unint = int_unint*nbf
      endif
      end
      integer function int_nint_raw(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
      implicit none
#include "bas.fh"
#include "int_nbf.fh"
c
c::passed
      integer ibasin   ! [input] basis set handle for icnt
      integer icnt     ! [input] contraction index (e.g., ish)
      integer jbasin   ! [input] basis set handle for jcnt
      integer jcnt     ! [input] contraction index (e.g., jsh)
      integer kbasin   ! [input] basis set handle for kcnt
      integer kcnt     ! [input] contraction index (e.g., ksh)
      integer lbasin   ! [input] basis set handle for lcnt
      integer lcnt     ! [input] contraction index (e.g., lsh)
c:local
      integer type, nprim, ngen, spcart
      integer nbf

      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_nint_raw = 0
        return
      endif

      int_nint_raw = 1

      
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_continfo(ibasin,icnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_raw: bas_continfo failure',911)
        nbf = int_nbf_x(type)
        int_nint_raw = int_nint_raw*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_continfo(jbasin,jcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_raw: bas_continfo failure',911)
        nbf = int_nbf_x(type)
        int_nint_raw = int_nint_raw*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_continfo(kbasin,kcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_raw: bas_continfo failure',911)
        nbf = int_nbf_x(type)
        int_nint_raw = int_nint_raw*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_continfo(lbasin,lcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_raw: bas_continfo failure',911)
        nbf = int_nbf_x(type)
        int_nint_raw = int_nint_raw*nbf
      endif
      end
      integer function int_unint_raw(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
      implicit none
#include "int_nbf.fh"
c
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
      integer type, nprim, ngen, spcart
      integer nbf
c::functions
      logical  bas_ucontinfo
      external bas_ucontinfo

      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_unint_raw = 0
        return
      endif

      int_unint_raw = 1

      
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_ucontinfo(ibasin,icnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)
        int_unint_raw = int_unint_raw*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_ucontinfo(jbasin,jcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)
        int_unint_raw = int_unint_raw*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_ucontinfo(kbasin,kcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)
        int_unint_raw = int_unint_raw*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_ucontinfo(lbasin,lcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_raw: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)
        int_unint_raw = int_unint_raw*nbf
      endif
      end
      integer function int_nint_cart(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
      implicit none
#include "bas.fh"
#include "int_nbf.fh"
c
c::passed
      integer ibasin   ! [input] basis set handle for icnt
      integer icnt     ! [input] contraction index (e.g., ish)
      integer jbasin   ! [input] basis set handle for jcnt
      integer jcnt     ! [input] contraction index (e.g., jsh)
      integer kbasin   ! [input] basis set handle for kcnt
      integer kcnt     ! [input] contraction index (e.g., ksh)
      integer lbasin   ! [input] basis set handle for lcnt
      integer lcnt     ! [input] contraction index (e.g., lsh)
c:local
      integer type, nprim, ngen, spcart
      integer nbf

      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_nint_cart = 0
        return
      endif

      int_nint_cart = 1

      
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_continfo(ibasin,icnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_continfo(jbasin,jcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_continfo(kbasin,kcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_continfo(lbasin,lcnt,type,nprim,ngen,spcart))
     &         call errquit('int_nint_cart: bas_continfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_nint_cart = int_nint_cart*nbf
      endif
      end
      integer function int_unint_cart(ibasin,icnt,jbasin,jcnt,
     &       kbasin,kcnt,lbasin,lcnt)
      implicit none
#include "int_nbf.fh"
c
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
      integer type, nprim, ngen, spcart
      integer nbf
c::functions
      logical  bas_ucontinfo
      external bas_ucontinfo

      if (
     &    icnt.eq.jcnt.and.
     &    jcnt.eq.kcnt.and.
     &    kcnt.eq.lcnt.and.
     &    icnt.eq.0
     &   ) then
        int_unint_cart = 0
        return
      endif

      int_unint_cart = 1

      
c..icnt/ish
      if (icnt.gt.0) then
        if (.not.bas_ucontinfo(ibasin,icnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
c
c..jcnt/jsh
      if (jcnt.gt.0) then
        if (.not.bas_ucontinfo(jbasin,jcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
c
c..kcnt/ksh
      if (kcnt.gt.0) then
        if (.not.bas_ucontinfo(kbasin,kcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
c
c..lcnt/lsh
      if (lcnt.gt.0) then
        if (.not.bas_ucontinfo(lbasin,lcnt,type,nprim,ngen,spcart))
     &         call errquit('int_unint_cart: bas_ucontinfo failure',911)
        nbf = int_nbf_x(type)*ngen
        int_unint_cart = int_unint_cart*nbf
      endif
      end
