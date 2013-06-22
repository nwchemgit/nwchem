C====================================================================
C
C     Wrapper routine for inverting generic complex matrix g_za.
C
      subroutine zmat_inv (g_za, g_zainv)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"


C     == Inputs ==
      integer, intent(in) :: g_za              !complex matrix to invert


C     == Outputs ==
      integer, intent(in) :: g_zainv           !inverse of matrix


C     == Parameter ==
      character(*), parameter :: pname = "zmat_inv: "


C     == Variables ==
      integer dim1, dim2, dtype
      integer n


C      
C     Get dims of GAs and check that they are correct types
C
C     Check the matrix (input 1).
C
      call ga_check_handle (g_za,
     $     "first argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_za, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as first argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)

C      
C     The size of all matricies must be n x n.
C
      n = dim1

      
C
C     Check the inverse (output) matrix.
C
      call ga_check_handle (g_zainv,
     $     "second argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zainv, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1.ne.n)
     $     call errquit (pname//"size of ainv must match size of a")


CXXX  [KAL]: only serial inverse routine works for now
      call zmat_inv_serial (n, g_za, g_zainv)
      
      end subroutine


C====================================================================
      subroutine zmat_inv_serial (n, g_za, g_zainv)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"


C     == Inputs ==
      integer, intent(in) :: n                 !size of mats
      integer, intent(in) :: g_za              !input matrix


C     == Outputs ==      
      integer, intent(in) :: g_zainv           !inverse of matrix


C     == Parameter ==
      character(*), parameter :: pname = "zmat_inv_serial: "


C     == Variables ==
      logical ok
      integer i, j
      integer me
      integer info
      integer g_zident
      integer izident, lzident
      integer lza, iza
      integer lzainv, izainv
      integer lpivot, ipivot

      
      call ga_sync ()
      
C     
C     Do all work on all processors (wasteful but works).
C     
      me = ga_nodeid ()


C     
C     Allocation.
C      
      ok=.true.
      ok=ok.and.ma_push_get(mt_dcpl, n*n, "complex mat A", lza,iza)
      ok=ok.and.ma_push_get(mt_dcpl, n*n, "complex mat Ainv",
     $     lzainv, izainv)
      ok=ok.and.ma_push_get(mt_dcpl, n*n, "complex ident mat",
     $     lzident, izident)
      ok=ok.and.ma_push_get(mt_int, n, "pivot indices", lpivot, ipivot)
      
      ok = ok .and. ga_duplicate (g_za, g_zident, "ident")
      if (.not.ok) call errquit (pname//"alloc failed", 0, 0)
      
      
C
C     Unpack buffer to GA
C
      call pack_ga2buffer_dcpl (g_za, dcpl_mb(iza))

      
C
C     LU factorization
C
      call zgetrf (n, n, dcpl_mb(iza), n, int_mb(ipivot), info)

      if (info .ne. 0)
     $     call errquit (pname//"LU factorization failed", 0, 0)


C
C     Make complex identity matrix buffer.
C
      call ga_zero (g_zident)
      call mat_set_ident (g_zident)
      call pack_ga2buffer_dcpl (g_zident, dcpl_mb(izident))

      
C
C     Compute inverse by solving AX = I (I is the nxn ident mat)
C
C      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      call zgetrs ("N", n, n, dcpl_mb(iza), n, int_mb(ipivot),
     $     dcpl_mb(izident), n, info)
      
      if (info .ne. 0)
     $     call errquit (pname//"Inversion failed", 0, 0)


C
C     Pack solution (stored in zident buffer) into output GA.
C
      call pack_buffer2ga_dcpl (dcpl_mb(izident), g_zainv)


C
C     Clean up.
C      
      ok = .true.
      ok = ok .and. ma_chop_stack (lza)
      ok = ok .and. ga_destroy (g_zident)
      if (.not.ok) call errquit (pname//"free failed", 0, 0)
      
      end subroutine


C====================================================================
      subroutine zmat_inv_check (g_za, g_zainv)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"

C     == Inputs ==
      integer, intent(in) :: g_za              !complex matrix to invert
      integer, intent(in) :: g_zainv           !inverse of matrix


C     == Parameter ==
      character(*), parameter     :: pname = "zmat_inv_check: "
      double complex, parameter   :: z1 = (1d0, 0d0)
      double complex, parameter   :: z0 = (0d0, 0d0)
      double precision, parameter :: thresh = 1d-8


C     == Variables ==
      integer dim1, dim2, dtype
      integer n
      integer g_zprod


C      
C     Get dims of GAs and check that they are correct types
C
C     Check the matrix (input 1).
C
      call ga_check_handle (g_za,
     $     "first argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_za, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as first argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)

C      
C     The size of all matricies must be n x n.
C
      n = dim1

      
C
C     Check the inverse (output) matrix.
C
      call ga_check_handle (g_zainv,
     $     "second argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zainv, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1.ne.n)
     $     call errquit (pname//"size of ainv must match size of a")


C
C     Allocation
C
      if (.not. ga_duplicate (g_za, g_zprod, "prod"))
     $     call errquit (pname//"failed to alloc prod", 0, GA_ERR)



C
C     Check that A^-1 A = I
C      
      call ga_zgemm ("N", "N", n, n, n, z1, g_zainv, g_za, z0, g_zprod)
      if (.not. mat_is_ident (g_zprod, thresh))
     $     call errquit (pname//"inversion check failed", 0, 0)
      

C
C     Clean up
C     
      if (.not. ga_destroy (g_zprod))
     $     call errquit (pname//"destroy failed", 0, 0)

      end subroutine


C====================================================================
      subroutine zmat_inv_example_driver ()
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"


C     == Parameters ==
      character(*), parameter :: pname = "zmat_inv_example_driver: "


C     == Variables ==
      integer me
      integer g_zmat, g_zmatinv
      double complex val


      me = ga_nodeid ()


C
C     EXAMPLE ONE
C
C     octave:1> A=[[1,2,0];[0,3,0];[2,-4,2]]
C     A =
C
C     1   2   0
C     0   3   0
C     2  -4   2
C
C     octave:2> inverse(A)
C     ans =
C
C     1.00000  -0.66667   0.00000
C     0.00000   0.33333   0.00000
C     -1.00000   1.33333   0.50000
C

      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout, "Inversion example one",
     $        20, .true.)
      endif
      

C
C     Allocation.
C
      if (.not. ga_create(mt_dcpl,3,3,"matrix" , 0, 0, g_zmat))
     $     call errquit ("failed to create mat", 0, 0)

      if (.not. ga_create(mt_dcpl,3,3,"matrix inv" , 0, 0, g_zmatinv))
     $     call errquit ("failed to create mat inv", 0, 0)
      
      call ga_zero (g_zmat)

      val = dcmplx (1d0, 0d0)
      call ga_put (g_zmat, 1, 1, 1, 1, val, 1)

      val = dcmplx (2d0, 0d0)
      call ga_put (g_zmat, 1, 1, 2, 2, val, 1)

      val = dcmplx (3d0, 0d0)
      call ga_put (g_zmat, 2, 2, 2, 2, val, 1)

      val = dcmplx (2d0, 0d0)
      call ga_put (g_zmat, 3, 3, 1, 1, val, 1)

      val = dcmplx (-4d0, 0d0)
      call ga_put (g_zmat, 3, 3, 2, 2, val, 1)

      val = dcmplx (2d0, 0d0)
      call ga_put (g_zmat, 3, 3, 3, 3, val, 1)


C
C     Inversion
C
      call zmat_inv (g_zmat, g_zmatinv)
      call zmat_inv_check (g_zmat, g_zmatinv)
      call ga_print (g_zmat)
      call ga_print (g_zmatinv)


C
C     EXAMPLE TWO
C
C     octave:1> A=[[1.2+i,1,i];[0,3,-0.98*i];[2+3.2*i,-4,2.5+8*i]]
C     A =
C
C     1.20000 + 1.00000i   1.00000 + 0.00000i   0.00000 + 1.00000i
C     0.00000 + 0.00000i   3.00000 + 0.00000i  -0.00000 - 0.98000i
C     2.00000 + 3.20000i  -4.00000 + 0.00000i   2.50000 + 8.00000i
C     
C     octave:2> inverse(A)
C     ans =
C     
C     0.867525 - 0.256532i  -0.512594 + 0.069857i  -0.167565 - 0.011740i
C     -0.073269 - 0.137812i   0.419917 + 0.105576i   0.046621 + 0.044729i
C     -0.421875 + 0.224292i   0.323190 - 0.265053i   0.136924 - 0.142717i
C

      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout, "Inversion example two",
     $        20, .true.)
      endif

      call ga_zero (g_zmat)

      val = dcmplx (1.2d0, 1d0)
      call ga_put (g_zmat, 1, 1, 1, 1, val, 1)

      val = dcmplx (1.0d0, 0d0)
      call ga_put (g_zmat, 1, 1, 2, 2, val, 1)

      val = dcmplx (0d0, 1d0)
      call ga_put (g_zmat, 1, 1, 3, 3, val, 1)

      val = dcmplx (0d0, 0d0)
      call ga_put (g_zmat, 2, 2, 1, 1, val, 1)

      val = dcmplx (3d0, 0d0)
      call ga_put (g_zmat, 2, 2, 2, 2, val, 1)

      val = dcmplx (0d0, -0.98d0)
      call ga_put (g_zmat, 2, 2, 3, 3, val, 1)

      val = dcmplx (2d0, 3.2d0)
      call ga_put (g_zmat, 3, 3, 1, 1, val, 1)

      val = dcmplx (-4d0, 0d0)
      call ga_put (g_zmat, 3, 3, 2, 2, val, 1)

      val = dcmplx (2.5d0, 8d0)
      call ga_put (g_zmat, 3, 3, 3, 3, val, 1)


C
C     Inversion
C
      call zmat_inv (g_zmat, g_zmatinv)
      call zmat_inv_check (g_zmat, g_zmatinv)
      call ga_print (g_zmat)
      call ga_print (g_zmatinv)


C
C     Clean up.
C     
      if (.not. ga_destroy (g_zmat))
     $     call errquit (pname//"failed to destroy g_zmat",0,0)
      if (.not. ga_destroy (g_zmatinv))
     $     call errquit (pname//"failed to destroy g_zmatinv",0,0)

      end subroutine
      
      
c $Id$
