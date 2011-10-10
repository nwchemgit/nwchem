      subroutine argos_cafe_free()
c
      implicit none
c
      call argos_cafe_lfree()
      call argos_cafe_sfree()
c
      return
      end
c $Id$
