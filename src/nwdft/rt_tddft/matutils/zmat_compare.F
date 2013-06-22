C====================================================================
      
      subroutine zmat_compare_dump (g_re, g_im, basename)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "matutils.fh"
      
C     == Inputs ==
      integer, intent(in)           :: g_re, g_im
      character(len=*), intent(in)  :: basename

      
C     == Parameters ==
      character(len=*), parameter :: pname = "zmat_compare_dump: "


C     == Variables ==
      character(len=255)  :: fname

      call util_file_name ("zmat_compare."//trim(basename)//".re",
     $     .false., .false., fname)
      
      if (.not. dmat_io_dump (g_re, fname))
     $     call errquit(pname//"dump of g_re failed: "//trim(fname),0,0)

      call util_file_name ("zmat_compare."//trim(basename)//".im",
     $     .false., .false., fname)
      
      if (.not. dmat_io_dump (g_im, fname))
     $     call errquit(pname//"dump of g_im failed: "//trim(fname),0,0)

      if (ga_nodeid().eq.0) then
         write (luout, *) pname//'saved complex matrix to "'
     $        //trim(fname)//'"'
      endif
      
      end subroutine 



C====================================================================
      
      subroutine zmat_compare_dump_purereal (g_re, basename)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "matutils.fh"
      
C     == Inputs ==
      integer, intent(in)           :: g_re
      character(len=*), intent(in)  :: basename

      
C     == Parameters ==
      character(len=*), parameter :: pname =
     $     "zmat_compare_dump_purereal: "


C     == Variables ==
      character(len=255)  :: fname
      integer g_zero

      if (.not. ga_duplicate (g_re, g_zero, "g_zero"))
     $     call errquit (pname//"duplicate failed", 0, 0)

      call ga_zero (g_zero)

      call util_file_name ("zmat_compare."//trim(basename)//".re",
     $     .false., .false., fname)
      if (.not. dmat_io_dump (g_re, fname))
     $     call errquit (pname//"dump of g_re failed", 0, 0)

      call util_file_name ("zmat_compare."//trim(basename)//".im",
     $     .false., .false., fname)
      if (.not. dmat_io_dump (g_zero, fname))
     $     call errquit (pname//"dump of g_im failed", 0, 0)  ! pure 0 imag part

      if (.not. ga_destroy (g_zero))
     $     call errquit (pname//"destroy failed", 0, 0)

      if (ga_nodeid().eq.0) then
         write (luout, *) pname//'saved real matrix to "'
     $        //trim(fname)//'"'
      endif

      end subroutine 

      


C====================================================================
      
      logical function zmat_compare (n, basename, g_zmat, tol)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "matutils.fh"
      
C     == Inputs ==
      integer, intent(in)           :: n         ! size of all matrices
      character(len=*), intent(in)  :: basename  ! base filename for mat to compare against
      integer, intent(in)           :: g_zmat    ! matrix to compare against
      double precision, intent(in)  :: tol


C     == Parameters ==
      character(len=*), parameter :: pname = "zmat_compare: "

      
C     == Variables
      integer g_tmpre, g_tmpim, g_ztmp
      character(len=255)  :: fname
      integer me
      integer type, dim1, dim2

      
      me = ga_nodeid ()

      if (.not. ga_create(mt_dbl, n, n, "tmpre", 0, 0, g_tmpre))
     $     call errquit (pname//"couldnt create g_tmpre", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, n, n, "tmpim", 0, 0, g_tmpim))
     $     call errquit (pname//"couldnt create g_tmpim", 0, GA_ERR)


      
      if (.not. ga_create(mt_dcpl, n, n, "ztmp", 0, 0, g_ztmp))
     $     call errquit (pname//"couldnt create g_ztmp", 0, GA_ERR)


C     (read in real part)
      call ga_zero (g_tmpre)
      call util_file_name ("zmat_compare."//trim(basename)//".re",
     $     .false., .false., fname)
      if (.not. dmat_io_read (g_tmpre, fname))
     $     call errquit (pname//"read re failed: " //trim(fname), 0, 0)

C     (read in imag part)
      call ga_zero (g_tmpim)
      call util_file_name ("zmat_compare."//trim(basename)//".im",
     $     .false., .false., fname)
      if (.not. dmat_io_read (g_tmpim, fname))
     $     call errquit (pname//"read im failed: "//trim(fname), 0, 0)


C     (pack into complex and compare with input matrix)
      call ga_zero (g_ztmp)
      call convert_d2z (1d0, g_tmpre, 1d0, g_tmpim, g_ztmp)

      if (mat_is_same (g_ztmp, g_zmat, tol)) then
         zmat_compare = .true.
      else
         if (me.eq.0) then
            write (luout, "(1x,a,es8.2)")
     $           pname//trim(basename)//": "//
     $           "Matrices differ by more than ", tol
            call util_flush (luout)
         endif

C         call elem_compare (g_ztmp, g_zmat, tol*10)  !dont do this--output too verbose

C         call ga_print (g_ztmp)
C         call ga_print (g_zmat)
         zmat_compare = .false.
      endif

      
      if (.not. ga_destroy (g_tmpre))
     $     call errquit (pname//"destroy g_tmpre failed",0,0)
      if (.not. ga_destroy (g_tmpim))
     $     call errquit (pname//"destroy g_tmpim failed",0,0)
      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"destroy g_ztmp failed",0,0)

      end function
      


C====================================================================

      subroutine zmat_compare_and_print (n, g_a, matname, tol, desc)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "matutils.fh"

      
C     == Inputs ==
      integer, intent(in)               :: n       ! size of matrices
      integer, intent(in)               :: g_a     ! GA to compare against
      character(len=*), intent(in)      :: matname ! name of other matrix to compare against (uses zmat_compare)
      double precision, intent(in)      :: tol     ! tolerance
      character(len=*), intent(in)      :: desc    ! description (for output)


C     == External ==
      logical, external  :: zmat_compare
      

C     == Parameters ==
      character(len=*), parameter :: pname = "zmat_compare_and_print: "
      character(len=*), parameter :: fmt = "(a,es8.2,a)"


C     == Variables ==
      integer :: me
      

C
C     Checks
C
      if (n < 1) call errquit (pname//"invalid n < 1", 0, 0)

      if (tol < 0d0) call errquit (pname//"invalid negative tol", 0, 0)


      me = ga_nodeid ()

      
      if (.not. zmat_compare (n, matname, g_a, tol)) then
         if (me.eq.0)
     $        write (luout, fmt) "*** COMPARISON FAILED: " //
     $        trim(desc) // " vs tolerance ", tol, " ***" 
      else
         if (me.eq.0)
     $        write (luout, fmt) "Comparison passed: " //
     $        trim(desc) // " vs tolerance ", tol, ""
      endif


      end subroutine
c $Id$
