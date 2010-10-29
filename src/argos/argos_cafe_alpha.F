      subroutine argos_cafe_alpha
c $Id$
      implicit none
c
#include "argos_cafe_common.fh"
#include "argos_cafe_funcs_dec.fh"
c
      real*8 a,alow,ahigh,b
      integer i,n
c
#include "argos_cafe_funcs_sfn.fh"
c
      a=half
      i=0
    1 continue
      a=two*a
      i=i+1
      b=(erfc(a*rlong))/rlong
      if(b.ge.pmetol) goto 1
c
      n=i+60
      alow=zero
      ahigh=a
      do 2 i=1,n
      a=half*(alow+ahigh)
      b=(erfc(a*rlong))/rlong
      if(b.ge.pmetol) then
      alow=a
      else
      ahigh=a
      endif
    2 continue
      ealpha=a
c
      return
      end
