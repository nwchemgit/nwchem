      subroutine riscf_tidy
C$Id: riscf_tidy.F,v 1.2 1995-12-29 02:04:15 vg038 Exp $
      implicit none
c
c     this is just a wrapper for del_3c_ints, to free ALL arrays of 3c 
c     integrals
c
      logical status, riscf_del_3c_ints

      status = riscf_del_3c_ints(-1)
      if (.not.status)
     $     call errquit('riscf_tidy: failed to free 3c ints', 0)

      return
      end
