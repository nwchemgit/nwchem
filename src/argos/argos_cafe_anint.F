#if defined(CRAY_T3D) || defined(CRAY_T3E)
      real*8 function argos_cafe_anint(x)
      real*8 x
      argos_cafe_anint=anint(x)
      return
      end
#endif
