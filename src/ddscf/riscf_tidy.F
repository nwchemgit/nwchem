      subroutine riscf_tidy
C$Id: riscf_tidy.F,v 1.1 1995-05-22 22:28:29 vg038 Exp $
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
