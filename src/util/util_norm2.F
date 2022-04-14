      function util_norm2(n,v)
      implicit none
      double precision :: util_norm2
      integer ,intent(in) :: n
      double precision,intent(in) :: v(n)
#if defined(XLFLINUX) || (defined(__GNUC__) && !defined(GCC46))
      util_norm2=dsqrt(dot_product(v,v))
#else      
      util_norm2=norm2(v)
#endif      
      return
      end function util_norm2
