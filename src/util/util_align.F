c $Id$
      subroutine util_align(in,alignval)
      implicit none
      integer in,alignval
      if(mod(in,alignval).ne.0) in=(1+in/alignval)*alignval
      in=in+alignval
      return
      end
