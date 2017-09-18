C
C     matutils.F
C
C     Various matrix utility rouintes.
C      
C
C     Manifest:
C-----------------------------
C
C     logical mat_is_symmetric (g_a, symant, tol)
C     logical mat_is_idempotent (g_z/a, tol)
C     logical mat_is_hermitian (g_z, tol)
C     logical mat_is_tracex (g_z/a, x, tol)
C     logical mat_is_ident (g_z/a, tol)
C     logical mat_is_same (g_z/a, g_z/a, tol)
C     logical mat_is_purereal (g_z, tol)
C     logical mat_is_unitary (g_z, tol)
C     logical mat_is_posneg (g_z/a, posneg, tol)
C
C     double precision mat_norm (g_z/a)
C     subroutine mat_dtrace (g_a, dtr)
C     subroutine mat_ztrace (g_z, ztr)
C
C     subroutine mat_set_ident (g_z/a)
C     subroutine mat_hermconj (g_z, g_zhc)
C      
C====================================================================
      logical function mat_is_symmetric (g_a, symant, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in)          :: g_a       !real matrix to analyze
      character*1, intent(in)      :: symant    !"S" check symmetric; "A" check antisymm
      double precision, intent(in) :: tol

C     == External ==
      double precision, external :: mat_norm
      logical, external          :: mat_is_same
      

C     == Variables ==      
      integer g_tmp
      integer dtype, n1, n2
      double precision norm
      character(len=20), parameter :: pname = "mat_is_symmetric:"

      
      mat_is_symmetric = .false.

C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)
      
      if (dtype .ne. mt_dbl) call errquit (
     $     pname//" only valid for real matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")


      if (.not. ga_duplicate (g_a, g_tmp, "work matrix"))
     $     call errquit (pname//"cannot create g_tmp",0,0)

      if (symant.eq."S") then   !check symmetric
         call ga_transpose (g_a, g_tmp)
         mat_is_symmetric = mat_is_same (g_a, g_tmp, tol)
         
      elseif (symant.eq."A") then !check antisymmetric
         call ga_transpose (g_a, g_tmp)
         call ga_scale (g_tmp, -1d0)
         mat_is_symmetric = mat_is_same (g_a, g_tmp, tol)
         
      else
         call errquit (pname//"bad flag '"//symant//"'", 0, 0)
      endif

      if (.not. ga_destroy (g_tmp))
     $     call errquit (pname//"cannot destroy g_tmp",0,0)
      
      end function
      
C====================================================================
C
C     Checks that A^2 = A; works for both real and complex.
C
C      
      logical function mat_is_idempotent (g_a, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"


C     == Inputs ==
      integer, intent(in)          :: g_a    !real or complex ga
      double precision, intent(in) :: tol

C     == External ==
      double precision, external :: mat_norm

C     == Variables ==
      integer g_tmp
      integer dtype, n1, n2
      double precision norm
      character(len=20), parameter :: pname = "mat_is_idempotent:"
      double complex, parameter :: zone = (1d0, 0d0)
      double complex, parameter :: zzero = (0d0, 0d0)

C     == External ==
      logical, external :: mat_is_same

      mat_is_idempotent = .false.

C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)

      if (.not. ga_duplicate (g_a, g_tmp, "work matrix"))
     $     call errquit (pname//"cannot create g_tmp",0,0)

      call ga_zero (g_tmp)

      call ga_inquire (g_a, dtype, n1, n2)

C     If matrix is square, certainly not idempotent.
      if (n1.ne.n2) then
         mat_is_idempotent = .false.
         return
      endif

C     Take matmul of the two, subtract original matrix, take norm.
C     Note we have dbl and dcpl versions

      if (dtype .eq. mt_dbl) then !real 
         call ga_dgemm ("n","n",n1,n1,n1, 1d0, g_a, g_a, 0d0, g_tmp)
         mat_is_idempotent = mat_is_same (g_a, g_tmp, tol)

      elseif (dtype .eq. mt_dcpl) then !complex
         call ga_zgemm ("n","n",n1,n1,n1, zone, g_a, g_a, zzero, g_tmp)
         mat_is_idempotent = mat_is_same (g_a, g_tmp, tol)

      else
         call errquit (pname//"unsupported data type", 0, 0)
         
      endif

      if (.not. ga_destroy (g_tmp))
     $     call errquit (pname//"cannot destroy g_tmp",0,0)
      
      end function
C====================================================================

      logical function mat_is_hermitian (g_z, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in)          :: g_z     ! complex mat to check
      double precision, intent(in) :: tol

C     == Parameters ==
      character*18, parameter :: pname = "mat_is_hermitian"

C     == Variables ==
      integer g_re, g_im
      integer n1, n2, dtype
      logical herm

C     == External ==
      logical mat_is_symmetric
      external mat_is_symmetric

C     == Preliminary checks ==
      call ga_check_handle (g_z,
     $     "first argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_z, dtype, n1, n2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as first argument", 0, 0)
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")
      

C     == Allocate and extract real and im parts ==
      if (.not.ga_create(mt_dbl, n1, n1,"Real part", 0, 0, g_re))
     $     call errquit (pname//"couldnt create g_re", 0, 0)
      if (.not.ga_create(mt_dbl, n1, n1,"Imag part", 0, 0, g_im))
     $     call errquit (pname//"couldnt create g_im", 0, 0)

      call convert_z2d (g_z, g_re, g_im)

C
C     If g_z is Hermitian, the real part must be symmetric
C     and the imaginary part must be antisymmetric.
C     
      herm = .true.
      herm = herm .and. mat_is_symmetric (g_re, "S", tol)
      herm = herm .and. mat_is_symmetric (g_im, "A", tol)

C     == Clean up ==
      if (.not.ga_destroy(g_re))
     $     call errquit (pname//"couldnt destroy g_re",0,0)
      if (.not.ga_destroy(g_im))
     $     call errquit (pname//"couldnt destroy g_im",0,0)
      
      mat_is_hermitian = herm

      end function

C====================================================================
C
C     Computes a kind of "norm" of a matrix, equivalent to f90:
C     sqrt(sum(abs(a**2))); works for both real and complex mats.
CXXX  [KAL]: uses build in norm infty instead!
C
      double precision function mat_norm (g_a)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      
C     == Inputs ==
      integer, intent(in) :: g_a               !real or complex ga

C     == Parameters ==
      character*16, parameter  :: pname = "mat_abs_sqrd: "

C     == Variables ==
      integer dtype, n1, n2
      integer g_tmp
      double precision ds2
      double complex zs2
      double precision n1sqinv
      
C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"() is not a valid GA")

      call ga_inquire (g_a, dtype, n1, n2)


#if 0
C      if (n1 .ne. n2)
C     $     call errquit (pname//"n1 must equal n2")

      if (.not. ga_duplicate (g_a, g_tmp, "work matrix"))
     $     call errquit (pname//"cannot create g_tmp",0,0)

C
C     Compute the norm, note we have real and complex versions,
C     but the norm is always real.
C
      call ga_copy (g_a, g_tmp)
      call ga_abs_value (g_tmp)

C     "Normalize" by dividing by N^2:
      n1sqinv = 1d0 / dble ((n1*n1))
      
      if (dtype.eq.mt_dbl) then !real matrix
         ds2 = ga_ddot (g_tmp, g_tmp)
         mat_norm = dsqrt (ds2) * n1sqinv
      elseif (dtype.eq.mt_dcpl) then !complex matrix
         zs2 = ga_zdot (g_tmp, g_tmp)
         mat_norm = dsqrt (abs (zs2)) * n1sqinv
      else
         call errquit (pname//"unsupported data type", 0, 0)         
      endif
#endif


CXXX  [KAL]: CHECKING USING BUILT-IN

      call ga_norm_infinity (g_a, mat_norm)


C      if (.not. ga_destroy (g_tmp))
C     $     call errquit (pname//"cannot destroy g_tmp",0,0)
      
      end function

C====================================================================
C
C     Calculates the trace of a real matrix.
C
      subroutine mat_dtrace (g_a, dtr)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs == 
      integer, intent(in) :: g_a               !real matrix


C     == Outputs ==
      double precision, intent(out) :: dtr

      
C     == Parameters ==
      character(len=15), parameter :: pname = "mat_dtrace:"

      
C     == Variables ==
      integer dtype, n1, n2
      integer g_ident

      
C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)
      
      if (dtype .ne. mt_dbl) call errquit (
     $     pname//" only valid for real matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")


C
C     Stupid way of calculating trace:  take the ddot between
C     the matrix and the identity matrix.
C
      if (.not. ga_create (mt_dbl, n1, n1, "ident mat", 0, 0, g_ident))
     $     call errquit (pname//"cannot create g_ident",0,0)

      call mat_set_ident (g_ident)
      dtr = ga_ddot (g_a, g_ident)

      if (.not. ga_destroy (g_ident))
     $     call errquit (pname//"cannot destroy g_ident",0,0)

      end subroutine


C====================================================================
C
C     Calculates the trace of a complex matrix.
C
      subroutine mat_ztrace (g_z, ztr)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"


C     == Inputs == 
      integer, intent(in) :: g_z               !complex matrix


C     == Outputs ==
      double complex, intent(out) :: ztr        !trace
      

C     == Parameters ==
      character(len=15), parameter :: pname = "mat_ztrace:"

      
C     == Variables ==
      integer dtype, n1, n2
      integer g_ident

      
C
C     Check the GA.
C
      call ga_check_handle (g_z,
     $     "first argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_z, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")


C     Stupid way of calculating trace:  take the ddot between
C     the matrix and the identity matrix
      if (.not. ga_create (mt_dcpl, n1, n1, "ident mat", 0, 0, g_ident))
     $     call errquit (pname//"cannot create g_ident",0,0)
      
      call mat_set_ident (g_ident)
      ztr = ga_zdot (g_z, g_ident)
      
      if (.not. ga_destroy (g_ident))
     $     call errquit (pname//"cannot destroy g_ident",0,0)
      
      end subroutine


C====================================================================
C
C     Sets the input GA to be the identity matrix.
C     Works for both complex and real GAs.
C
      subroutine mat_set_ident (g_a)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"


C     == Inputs == 
      integer, intent(in) :: g_a               !real or complex matrix

      
C     == Parameters ==
      character(len=*), parameter :: pname = "mat_set_trace: "
      double complex, parameter   :: z1 = (1d0, 0d0)

      
C     == Variables ==
      integer dtype, n1, n2
      integer me
      integer i
      double precision one
      double complex zone


      me = ga_nodeid ()

C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "first argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")


      if (dtype .eq. mt_dbl) then !real ident mat

#define NEWIDENT 1

CXXX  [KAL]: silly way to make ident mat
CXXX  [KAL]: diagonal GA routines scramble???

#if NEWIDENT
         one = 1d0
         call ga_zero (g_a)
         if (me.eq.0) then
            do i = 1, n1
               call ga_put (g_a, i,i,i,i, one, 1)
            enddo
         endif
         call ga_sync ()
#else
         call ga_zero (g_a)
         call ga_shift_diagonal (g_a, 1d0)
#endif
         
      elseif (dtype .eq. mt_dcpl) then !complex ident mat

#if NEWIDENT
         zone = dcmplx (1d0, 0d0)
         call ga_zero (g_a)
         if (me.eq.0) then
            do i = 1, n1
               call ga_put (g_a, i,i,i,i, zone, 1)
            enddo
         endif
         call ga_sync ()
#else
         call ga_zero (g_a)
         call ga_shift_diagonal (g_a, z1)
#endif
         
      else
         call errquit (pname//"unsupported data type")         
      endif   
      
      end subroutine

C====================================================================
C
C     Checks if matrix has trace = x, where x is a real number.
C
      logical function mat_is_tracex (g_a, x, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"

C     == Inputs ==
      integer, intent(in)          :: g_a      !real or complex matrix
      double precision, intent(in) :: x        !desired trace
      double precision, intent(in) :: tol
      

C     == Parameters ==
      character(len=20), parameter :: pname = "mat_is_tracex: "


C     == Variables ==
      integer dtype, n1, n2
      double precision dtrace_n1 !Tr[A] - 1d0
      double complex ztrace_n1


      mat_is_tracex = .false.

C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)

      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")

      
C
C     Find trace of the real or complex matrix.
C     Check that is is ~x.
C    
      if (dtype.eq.mt_dbl) then !real matrix
         call mat_dtrace (g_a, dtrace_n1)
         dtrace_n1 = dtrace_n1 - x
         mat_is_tracex = abs(dtrace_n1).lt.tol
         
      elseif (dtype.eq.mt_dcpl) then !complex matrix
         call mat_ztrace (g_a, ztrace_n1)
         ztrace_n1 = ztrace_n1 - dcmplx (x, 0d0)
         mat_is_tracex = abs(ztrace_n1).lt.tol

      else
         call errquit (pname//"unsupported data type", 0, 0)
      endif


      end function 


C====================================================================
C
C     Compares two matricies to see if they differ.
C     Works for mixed real and/or complex matricies.
C
      logical function mat_is_same (g_a, g_b, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

C     == Inputs ==
      integer, intent(in)          :: g_a, g_b
      double precision, intent(in) :: tol


C     == Variables ==
      integer dtype_a, dtype_b
      integer n1_a, n2_a, n1_b, n2_b
      double precision norm
      logical same
      integer g_tmp
      integer me

C     == Parameters ==
      character(*), parameter      :: pname = "mat_is_same: "
      double complex, parameter    :: zone = (1d0, 0d0)
      double complex, parameter    :: nzone = (-1d0, 0d0)
      
C     == External ==
      double precision, external :: mat_norm

      me = ga_nodeid ()
      
C     Check GAs and get their data types
      call ga_check_handle (g_a, "1st arg of "//pname//"not valid GA")
      call ga_check_handle (g_b, "2nd arg of "//pname//"not valid GA")

      call ga_inquire (g_a, dtype_a, n1_a, n2_a)
      call ga_inquire (g_b, dtype_b, n1_b, n2_b)

      if ((n1_a .ne. n1_b).or.(n2_a .ne. n2_b))
     $     call errquit (pname//"mats must be same size", 0, 0)

      same = .false.
      
      if (dtype_a .eq. dtype_b)  then !they are the same data type

         if (dtype_a .eq. mt_dbl) then ! real data
            if (.not. ga_create (mt_dbl,n1_a,n1_a, "tmp", 0, 0, g_tmp))
     $           call errquit (pname//"cannot create g_tmp",0,0)

            call ga_zero (g_tmp)
            call ga_dadd (1d0, g_a, -1d0, g_b, g_tmp)
            norm = mat_norm (g_tmp)

            if (norm < tol) then
               same = .true.
            else
               same = .false.
               if (me.eq.0) write (luout, *)
     $              pname//"|| A - B || = ", norm
               call util_flush (luout)
            endif

            if (.not. ga_destroy (g_tmp))
     $           call errquit (pname//"cannot destroy g_tmp",0,0)

         elseif (dtype_a .eq. mt_dcpl) then ! complex data
            if (.not. ga_create (mt_dcpl,n1_a,n1_a, "tmp", 0, 0, g_tmp))
     $           call errquit (pname//"cannot create g_tmp",0,0)

            call ga_zero (g_tmp)
C            call ga_zadd (zone, g_a, nzone, g_b, g_tmp)
            call ga_add (zone, g_a, nzone, g_b, g_tmp)
            norm = mat_norm (g_tmp)

            if (norm < tol) then
               same = .true.
            else
               same = .false.
               if (me.eq.0) write (luout, *)
     $              pname//"|| A - B || = ", norm
               call util_flush (luout)
            endif

            if (.not. ga_destroy (g_tmp))
     $           call errquit (pname//"cannot destroy g_tmp",0,0)

      else
         call errquit (pname//"unsupported data type", 0, 0)
      endif                     ! ends real/complex data check
            
      else                      ! they are different data types
         call errquit (pname//"both must be same data type", 0, 0)
         
      endif                     ! ends same data check
      
      mat_is_same = same

      end function


C====================================================================
C
C     Checks if matrix is the identity.
C     
C
      logical function mat_is_ident (g_a, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

C     == Inputs ==
      integer, intent(in)          :: g_a        !matrix to check
      double precision, intent(in) :: tol


C     == Parameters ==
      character(len=*), parameter :: pname = "mat_is_ident:"

      
C     == External ==
      logical, external :: mat_is_same


C     == Variables ==
      integer g_ident
      integer dtype, n1, n2


      mat_is_ident = .false.
      
C
C     Check the GA.
C
      call ga_check_handle (g_a,
     $     "argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)
      
      if (n1 .ne. n2) then
         mat_is_ident = .false.
         return
      endif

      if (dtype .eq. mt_dbl) then !real matrix
         if (.not. ga_create (mt_dbl, n1, n1, "I", 0, 0, g_ident))
     $        call errquit (pname//"cannot create g_ident",0,0)

         call mat_set_ident (g_ident)
         mat_is_ident = mat_is_same (g_a, g_ident, tol)

         if (.not. ga_destroy (g_ident))
     $        call errquit (pname//"cannot destroy g_ident",0,0)

         
      elseif (dtype .eq. mt_dcpl) then !complex matrix
         if (.not. ga_create (mt_dcpl, n1, n1, "I", 0, 0, g_ident))
     $        call errquit (pname//"cannot create g_ident",0,0)

         call mat_set_ident (g_ident)
         mat_is_ident = mat_is_same (g_a, g_ident, tol)

         if (.not. ga_destroy (g_ident))
     $        call errquit (pname//"cannot destroy g_ident",0,0)

      else
         call errquit (pname//"unsupported data type",0,0)
      endif
      
      end function



C====================================================================
C
C     Checks if complex matrix has no imaginary part.
C     
C
      logical function mat_is_purereal (g_za, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

C     == Inputs ==
      integer, intent(in)          :: g_za          !matrix to check
      double precision, intent(in) :: tol


C     == Parameters ==
      character(len=*), parameter :: pname = "mat_is_purereal: "
      
C     == External ==
      double precision, external :: mat_norm


C     == Variables ==
      integer dtype, n1, n2
      double precision norm
      integer g_are, g_aim


C
C     Check the GA.
C
      call ga_check_handle (g_za,
     $     "argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_za, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matricies", 0, 0)

C     == Allocation == 
      if (.not.ga_create(mt_dbl,n1,n2,"Re[A]", 0, 0, g_are))
     $     call errquit ("failed to create g_are", 0, GA_ERR)
      if (.not.ga_create(mt_dbl,n1,n2,"Im[A]", 0, 0, g_aim))
     $     call errquit ("failed to create g_aim", 0, GA_ERR)


C     == Extract imaginary part and see that it is small ==
      call convert_z2d (g_za, g_are, g_aim)
      norm = mat_norm (g_aim)
      
      if (norm < tol) then
         mat_is_purereal = .true.
      else
         mat_is_purereal = .false.
      endif            

C     == Clean up ==
      if (.not.ga_destroy(g_are))
     $     call errquit ("failed to destroy g_are", 0, GA_ERR)
      if (.not.ga_destroy(g_aim))
     $     call errquit ("failed to destroy g_aim", 0, GA_ERR)

      end function

C====================================================================
C
C     Returns Hermitian conjugate A^+ of matrix.
C     
      subroutine mat_hermconj (g_z, g_zhc)
      implicit none
      
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

      
C     == Inputs ==
      integer, intent(in) :: g_z
      

C     == Outputs ==
      integer, intent(in) :: g_zhc             !hermitian conjugate of g_z


C     == Parameters ==
      character(*), parameter :: pname = "mat_hermconj:"
      double precision, parameter :: znone = (-1d0, 0d0)

      
C     == Variables ==      
      integer dtype, n1, n2
      integer g_dscr1, g_dscr2
      logical ok
      
 
C
C     Check the GA.
C
      call ga_check_handle (g_z,
     $     "first argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_z, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matricies", 0, 0)

      ok = .true.
      ok = ok .and. ga_create(mt_dbl, n2, n1, "dscr1", 0, 0, g_dscr1)
      ok = ok .and. ga_create(mt_dbl, n2, n1, "dscr2", 0, 0, g_dscr2)
      if (.not.ok) call errquit (pname//"alloc failed", 0, GA_ERR)

      call ga_zero (g_dscr1)
      call ga_zero (g_dscr2)
      call ga_zero (g_zhc)

      
      call ga_transpose (g_z, g_zhc)
      call convert_z2d (g_zhc, g_dscr1, g_dscr2)
      call convert_d2z (1d0, g_dscr1, -1d0, g_dscr2, g_zhc)

      
      ok = .true.
      ok = ok .and. ga_destroy(g_dscr1)
      ok = ok .and. ga_destroy(g_dscr2)
      if (.not.ok) call errquit (pname//"destroy failed", 0, GA_ERR)

      end subroutine

C====================================================================
      logical function mat_is_unitary (g_z, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

      
C     == Inputs ==
      integer, intent(in)          :: g_z
      double precision, intent(in) :: tol

      
C     == Parameters ==
      character(*), parameter :: pname = "mat_is_unitary:"
      double complex, parameter :: zone = (1d0, 0d0)
      double complex, parameter :: zzero = (0d0, 0d0)
      
C     == External ==
      logical, external :: mat_is_ident

      
C     == Variables ==      
      integer dtype, n1, n2
      integer g_zscr1, g_zscr2
      

C
C     Check the GA.
C
      call ga_check_handle (g_z,
     $     "first argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_z, dtype, n1, n2)
      
      if (dtype .ne. mt_dcpl) call errquit (
     $     pname//" only valid for complex matricies", 0, 0)
      
      if (n1 .ne. n2)
     $     call errquit (pname//"n1 must equal n2")

      if (.not.ga_duplicate (g_z, g_zscr1, "zscr1"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      if (.not.ga_duplicate (g_z, g_zscr2, "zscr2"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      
C
C     Check if A^+ A = I.
C     
      call ga_zero (g_zscr1)
      call mat_hermconj (g_z, g_zscr1)

      call ga_zero (g_zscr2)
      call ga_zgemm ("N", "N", n1, n1, n1,
     $     zone, g_zscr1, g_z, zzero, g_zscr2)

      mat_is_unitary = mat_is_ident (g_zscr2, tol)

      if (.not.ga_destroy (g_zscr1))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not.ga_destroy (g_zscr2))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      end function



C====================================================================
C
C     Checks if matrix is purely positive or purely negative.
C     
C
      logical function mat_is_posneg (g_a, posneg, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

C     == Inputs ==
      integer, intent(in)          :: g_a        !matrix to check
      character, intent(in)        :: posneg     !"P" for positive, "N" for negative
      double precision, intent(in) :: tol


C     == Parameters ==
      character(len=*), parameter :: pname = "mat_is_posneg:"

      
C     == Variables ==
      integer dtype, n1, n2
      integer me
      integer i, j
      double precision val
      logical stat
C      double complex zval

      
      me = ga_nodeid ()

C
C     Check the GA and params
C
      call ga_check_handle (g_a,
     $     "argument of "//pname//"() is not a valid GA")
      
      call ga_inquire (g_a, dtype, n1, n2)

      if (tol .lt. 0)
     $     call errquit (pname//"invalid negative tol",0,0)

      
C     (not square is ok)
c$$$      if (n1 .ne. n2) then
c$$$         mat_is_posneg = .false.
c$$$         return
c$$$      endif

      stat = .true.
      if (me.eq.0) then
         do i = 1, n1
            do j = 1, n2

               if (dtype .eq. mt_dbl) then
                  call ga_get (g_a, i, i, j, j, val, 1)
               elseif (dtype .eq. mt_dcpl) then
                  call errquit (pname//"dcpl not implemented yet",0,0)
               else
                  call errquit (pname//"unsupported data type",0,0)
               endif

               if (posneg .eq. "P") then !check that mat is purely positive
                  if (val .lt. -1d0*tol) then
                     write (luout, "(a,i0,a,i0,a,e12.3)")
     $                    pname//" Element ", i, ", ", j,
     $                    " not positive: ", val
                     stat = .false.
                  endif
                  
               elseif (posneg .eq. "N") then !check that mat is purely negative
                  if (val .gt. 1d0*tol) then
                     write (luout, "(a,i0,a,i0,a,e12.3)")
     $                    pname//" Element ", i, ", ", j,
     $                    " not negative: ", val
                     stat = .false.
                  endif
                  
               else
                  call errquit (pname//
     $                 "invalid posneg switch: "//posneg,0,0)
               endif

            enddo
         enddo
      endif
      
      call ga_sync ()
      call ga_brdcst (mt_log, stat, 1, 0)
      
      mat_is_posneg = stat

      end function


c $Id$
