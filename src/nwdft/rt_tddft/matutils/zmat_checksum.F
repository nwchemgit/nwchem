C
C     zmat_checksum.F
C
C     Computes the "checksum" of a complex GA (used to confirm old and
C     new matrices are the same).
C
C
      double precision function zmat_checksum (g_z)
      
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in) :: g_z        ! matrix to output


C     == Parameters ==
      character(*), parameter :: pname = "zmat_checksum: "

      
C     == Variables ==
      integer g_re, g_im
      integer dtype, n1, n2
      logical ok
      double precision sum1, sum2


C
C     Check the GA.
C
      call ga_check_handle (g_z,
     $     "first argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_z, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matrices", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")



C
C     Allocate real and imag parts.
C
      ok = .true.
      ok = ok.and.ga_create(mt_dbl, n1, n1, "Re part", 0, 0, g_re)
      ok = ok.and.ga_create(mt_dbl, n1, n1, "Im part", 0, 0, g_im)

      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      call convert_z2d (g_z, g_re, g_im)
      

C
C     Just some random checksum I came up with...
C
      call ga_norm1 (g_re, sum1)
      call ga_norm_infinity (g_re, sum2)

      zmat_checksum = 1.67*sum1 + 0.12*sqrt (abs(sum2))

      
C
C     Clean up.
C
      ok = .true.
      ok = ok .and. ga_destroy(g_re)
      ok = ok .and. ga_destroy(g_im)
      if (.not.ok) call errquit(pname//"failed to destroy GAs", 0,0)
      
      end function
      
c $Id$
