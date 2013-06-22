C====================================================================
C
C     Wrapper routine for diagonalizing generic complex matrix g_za.
C     Note that checks will probably fail if you reorder...
C
      subroutine zmat_diag (g_za, g_zevecs, zevals, lsort, lcheck)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"


C     == Inputs ==
      integer, intent(in) :: g_za               !complex matrix to diagonalize
      logical, intent(in) :: lsort              !sort by eigenvalue or not
      logical, intent(in) :: lcheck             !check result


C     == Outputs ==
      integer, intent(out)          :: g_zevecs    !complex eigenvectors
      double complex, intent(inout) :: zevals(*)   !complex eigenvalues


C     == Parameter ==
      character(*), parameter :: pname = "zmat_diag: "


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
C     Check the eigenvector matrix.
C
      call ga_check_handle (g_zevecs,
     $     "second argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zevecs, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1.ne.n)
     $     call errquit (pname//"size of evecs must match size of mat")


CXXX  [KAL]: only serial diagonalizer works for now
      call zmat_diag_serial (n, g_za, g_zevecs, zevals)


      if (lsort) call zmat_diag_sort (g_zevecs, zevals)
      if (lcheck) call zmat_diag_check (g_za, g_zevecs, zevals, 1d-6)  !CHECK WILL FAIL WITH REORDERD NON-HERMITIAN MATS
      
      end subroutine
      

C====================================================================
C
C     Serial complex diagonalizer using stock zgeev blas routine.
C
C     We use zgeev (for non-symmetric mat) to compute the *right*
C     eigenvectors for a non-Hermitian matricies, and zheev for
C     Hermitian matricies.
C
C     A X = a X .
C

      subroutine zmat_diag_serial (n, g_za, g_zevecs, zevals)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"


C     == Inputs ==
      integer, intent(in) :: n                 !size of mats
      integer, intent(in) :: g_za


C     == Outputs ==      
      integer, intent(in)         :: g_zevecs
      double complex, intent(out) :: zevals(n)


C     == Parameter ==
      character(*), parameter :: pname = "zmat_diag_serial: "


C     == Variables ==
      logical ok
      integer me
      integer lza, iza
      integer levecs_right, ievecs_right, ievecs_left, levecs_left
      integer levals_re, ievals_re
      integer llwork, ilwork
      integer lrwork, irwork
      integer nlwork, nrwork
      integer info
      integer i
      double complex zval
      character*100 aname

      logical lhermitian


      call ga_sync ()


      
      me = ga_nodeid ()


      call ga_inquire_name (g_za, aname)

      

      if (mat_is_hermitian (g_za, 1d-5)) then
         lhermitian = .true.
         if (me.eq.0) write (luout, *) pname//
     $        'Matrix "'//trim(aname)//
     $        '" is Hermitian; using zheev'
      else
         lhermitian = .false.
         if (me.eq.0) write (luout, *) pname//
     $        'Matrix "'//trim(aname)//
     $        '" is not Hermitian; using zgeev'

         if (me.eq.0) write (luout, *) pname//
     $        "CAREFUL: non-Hermitian may be buggy"
      endif   
      

      if (lhermitian) then

C
C     Allocate buffers
C
         nlwork = max(1, 2*n-1)
         nrwork = max(1, 3*n-2)
         
         ok=.true.
         ok=ok.and.ma_push_get(mt_dcpl, n*n, "complex mat A", lza,iza)
         ok=ok.and.ma_push_get(mt_dbl, n, "real evals",
     $        levals_re, ievals_re)  !note we need to make this evals array passed to routine input is complex
         ok=ok.and.ma_push_get(mt_dcpl, nlwork, "Lwork", llwork, ilwork)
         ok=ok.and.ma_push_get(mt_dbl, nrwork, "Rwork", lrwork, irwork)
         
         
C     Load input matrix GA into local buffer.
         call pack_ga2buffer_dcpl (g_za, dcpl_mb(iza))

         call ga_sync ()
         call zheev ("V", "U", n, dcpl_mb(iza), n,
     $        dbl_mb(ievals_re), dcpl_mb(ilwork), nlwork,
     $        dbl_mb(irwork), info)
         call ga_sync ()

         if (info.lt.0) call errquit (pname//"bad args to zheev",0,0)
         if (info.gt.0) call errquit (pname//"zheev failed",0,0)

C     load eigenvectors and values for shipping out
         call ga_zero (g_zevecs)
         call pack_buffer2ga_dcpl (dcpl_mb(iza), g_zevecs)
         
C     hermitian matrix therefore all evals pure real:
         do i = 1, n
            zevals(i) = dcmplx ( dbl_mb(ievals_re + i - 1), 0d0)
         enddo
         

      else ! general (non-Hermitian) case
         
C     
C     Allocate buffers.
C     
         nlwork = max(1, 2*n)
         nrwork = 2*n

         ok=.true.
         ok=ok.and.ma_push_get(mt_dcpl, n*n, "complex mat A", lza,iza)
         ok=ok.and.ma_push_get(mt_dcpl, n*n, "evecs_right",
     $        levecs_right, ievecs_right)
         ok=ok.and.ma_push_get(mt_dcpl, n*n, "evecs_left",
     $        levecs_left, ievecs_left)
         ok=ok.and.ma_push_get(mt_dcpl, nlwork, "Lwork", llwork, ilwork)
         ok=ok.and.ma_push_get(mt_dbl, nrwork, "Rwork", lrwork, irwork)
         if (.not.ok) call errquit (pname//"alloc failed", 0, MA_ERR)

         
C     Load input matrix GA into local buffer.
         call pack_ga2buffer_dcpl (g_za, dcpl_mb(iza))

C      
C     Call lapack generic complex diagonalizer.
C
         call ga_sync ()
         call zgeev ("V", "V", n, dcpl_mb(iza), n, zevals,
     $        dcpl_mb(ievecs_left), n, dcpl_mb(ievecs_right), n, 
     $        dcpl_mb(ilwork), nlwork, dbl_mb(irwork), info)
         call ga_sync ()
c         
         if (info.lt.0) call errquit (pname//"bad args to zgeev",0,0)
         if (info.gt.0) call errquit (pname//"zgeev failed",0,0)
C     
C     Load eigenvectors into GA (output); note we want the "right"
C     eigenvectors: A v_j = a_j v_j
C     
         call ga_zero (g_zevecs)
         call pack_buffer2ga_dcpl (dcpl_mb(ievecs_right), g_zevecs)
      
      endif

      if (.not.ma_chop_stack(lza))
     $     call errquit (pname//"failed to chop stack", 0, MA_ERR)
      
      end subroutine


C====================================================================
C
C     Sort eigenvalues in acending order (and corresponding
C     eigenvectors).  Uses a simple bubble sort on proc 0 (very slow but
C     it works).
C
C
      subroutine zmat_diag_sort (g_zevecs, zevals)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"

C     == Inputs ==
      integer, intent(inout)          :: g_zevecs
      double complex, intent(inout)   :: zevals(*)


C     == Parameters ==
      character(*), parameter   :: pname = "zmat_diag_sort: "
      double complex, parameter :: zzero = (0d0, 0d0)


C     == Variables ==
      integer dim1, dim2, dtype
      integer n
      integer me
      integer i
      logical swapped

      double complex v1, v2
      integer ivec1, lvec1, ivec2, lvec2
      

      
      me = ga_nodeid()
      

C
C     Check the evec matrix (input 1).
C
      call ga_check_handle (g_zevecs,
     $     "first argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zevecs, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as first argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)

      
C      
C     The size of all matricies must be n x n.
C
CXXX  [KAL]: should prbably check that zevecs array is size n as well...
C
      n = dim1


      
C
C     Allocation
C      
      if (.not. ma_push_get (mt_dcpl,n,"vec1",
     $     lvec1, ivec1))
     $     call errquit (pname//"failed to push vec1", 0, MA_ERR)

      if (.not. ma_push_get (mt_dcpl,n,"vec2",
     $     lvec2, ivec2))
     $     call errquit (pname//"failed to push vec2", 0, MA_ERR)



C
C     Sort by eigenvalues in acending order by REAL part
C
      if (me.eq.0) then
         do while (.true.)
            
            swapped = .false.
            do i = 1, n - 1
               v1 = zevals(i)
               v2 = zevals(i+1)
               
               if (dble(v1) .gt. dble(v2)) then !swap eigenvalues and eigenvectors

                  zevals(i)   = v2
                  zevals(i+1) = v1

                  call ga_get(g_zevecs, 1,n, i,   i  , dcpl_mb(ivec1),n)
                  call ga_get(g_zevecs, 1,n, i+1, i+1, dcpl_mb(ivec2),n)
                  
                  call ga_put(g_zevecs, 1,n, i,   i  , dcpl_mb(ivec2),n)
                  call ga_put(g_zevecs, 1,n, i+1, i+1, dcpl_mb(ivec1),n)

                  swapped = .true.
               endif
            enddo
            
            if (.not. swapped) exit
         enddo
      endif
      call ga_sync ()


C
C     Clean up
C
      if (.not. ma_chop_stack (lvec1))
     $     call errquit (pname//"failed to chop stack", 0, MA_ERR)
      
      end subroutine 

      

C====================================================================
C
C     Complex parallel diagonalizer.
C
CXXX  [KAL]: STUB ONLY
      subroutine zmat_diag_parallel ()
      implicit none
      
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"

      call errquit ("zmat_diag_parallel not implemented yet",0,0)

      end subroutine


C====================================================================
C
C     Check diagonalization.
C      
      subroutine zmat_diag_check (g_za, g_zevecs, zevals, tol)
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "matutils.fh"

C     == Inputs ==
      integer, intent(in)          :: g_za
      integer, intent(in)          :: g_zevecs
      double complex, intent(in)   :: zevals(*)
      double precision, intent(in) :: tol


C     == Parameter ==
      character(*), parameter   :: pname = "zmat_diag_check: "
      double complex, parameter :: zone = (1d0, 0d0)
      double complex, parameter :: zzero = (0d0, 0d0)
      double complex, parameter :: znone = (-1d0, 0d0)
      double complex, parameter :: zjunk = (-23d0, 5.6d0)


C     == Variables ==
      integer dim1, dim2, dtype
      integer n
      integer solve_result
      integer g_zx, g_zvec, g_zavec
      integer ivec, lvec, iavec, lavec
      integer i, j, me
      double precision diff

      
      me = ga_nodeid()


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
C     Check the eigenvector matrix.
C
      call ga_check_handle (g_zevecs,
     $     "second argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zevecs, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1.ne.n)
     $     call errquit (pname//"size of evecs must match size of mat")



C
C     Allocation
C      
      if (.not. ga_duplicate (g_za, g_zx, "zX"))
     $     call errquit (pname//"failed to alloc zX", 0, GA_ERR)

      if (.not. ga_create (mt_dcpl, n, 1, "g_zvec", 0, 0, g_zvec))
     $     call errquit (pname//"failed to alloc g_zvec", 0, GA_ERR)

      if (.not. ga_create (mt_dcpl, n, 1, "g_zavec", 0, 0, g_zavec))
     $     call errquit (pname//"failed to alloc g_zavec", 0, GA_ERR)
      
      if (.not. ma_push_get (mt_dcpl, n, "vec", lvec, ivec))
     $     call errquit (pname//"failed to push vec", 0, MA_ERR)

      if (.not. ma_push_get (mt_dcpl, n, "avec", lavec, iavec))
     $     call errquit (pname//"failed to push avec", 0, MA_ERR)

      
      do i = 1, n               !loop over eigenvectors
C
C     Get v_i.
C
         call ga_zero (g_zvec)

         call ga_copy_patch ("N", g_zevecs, 1, n, i, i,
     $        g_zvec, 1, n, 1, 1)

         
C
C     Check that A v_i = a_i v_i .
C
C
C     A v_i
C
         call ga_zero (g_zavec)
         call ga_zgemm ("N", "N", n, 1, n,
     $        zone, g_za, g_zvec, zzero, g_zavec)


C
C     A v_i - a_i v_i
C
         call ga_add (zone, g_zavec, -1d0*zevals(i), g_zvec, g_zvec)

         
C
C     Norm of the difference vector
C         
C         diff = mat_norm (g_zvec)
         diff = sqrt (ga_zdot (g_zvec, g_zvec))
         
         if (diff.gt.tol)
     $        call errquit (pname//"diag check failed", 0, 0)
         
         call ga_sync ()
      enddo                     !end loop over eigenvectors

      
      if (.not.ga_destroy(g_zx))
     $     call errquit (pname//"failed to destroy g_zx", 0, GA_ERR)

      if (.not.ga_destroy(g_zvec))
     $     call errquit (pname//"failed to destroy g_zvec", 0, GA_ERR)

      if (.not.ga_destroy(g_zavec))
     $     call errquit (pname//"failed to destroy g_zavec", 0, GA_ERR)

      if (.not. ma_chop_stack (lvec))
     $     call errquit (pname//"failed to chop stack", 0, MA_ERR)
     
      end subroutine
      


C====================================================================
      subroutine zmat_diag_print_eigeninfo (g_za, g_zevecs, zevals)
      implicit none
      
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"


C     == Inputs ==
      integer, intent(in)        :: g_za       !complex matrix to diagonalize
      integer, intent(in)        :: g_zevecs   !complex eigenvectors
      double complex, intent(in) :: zevals(*)  !complex eigenvalues


C     == Parameter ==
      character(*), parameter :: pname = "zmat_diag_print_eigeninfo: "


C     == Variables ==
      integer i
      integer me
      integer dim1, dim2, dtype
      integer n


      
      me = ga_nodeid ()


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
C     Check the eigenvector matrix.
C
      call ga_check_handle (g_zevecs,
     $     "second argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zevecs, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1.ne.n)
     $     call errquit (pname//"size of evecs must match size of mat")



C
C     Print info to stdout
C
      if (me.eq.0) then 
         write (luout, *) ""
         call util_print_centered (luout, "Eigensystem information",
     $        40, .true.)
      endif


      if (me.eq.0) then 
         write (luout, *) ""
         call util_print_centered (luout, "Original matrix", 0, .true.)
      endif
      call ga_print (g_za)

C     (original matrix)
      if (me.eq.0) then 
         write (luout, *) ""
         call util_print_centered (luout, "Eigenvalues", 0, .true.)
      endif

C     (eigenvalues and abs eigenvalues)
      do i = 1, n
         if (me.eq.0)
     $        write (luout,"(i8,a,1es16.8e3,a,1es16.8e3,a,1es16.8e3,a)")
     $        i, ": " , dble (zevals(i)), " + ",
     $        aimag (zevals(i)), " i   (", abs(zevals(i)), " )"
         call util_flush (luout)
      enddo

C     (eigenvectors)
      if (me.eq.0) then 
         write (luout, *) ""
         call util_print_centered (luout, "Eigenvector matrix",0,.true.)
      endif
      call ga_print (g_zevecs)

      if (me.eq.0) write (luout, *) ""
      if (me.eq.0) write (luout, *) ""

      end subroutine



C====================================================================
      subroutine zmat_diag_print_eigenvals (g_za, g_zevecs, zevals)
      implicit none
      
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"


C     == Inputs ==
      integer, intent(in)        :: g_za       !complex matrix to diagonalize
      integer, intent(in)        :: g_zevecs   !complex eigenvectors
      double complex, intent(in) :: zevals(*)  !complex eigenvalues


C     == Parameter ==
      character(*), parameter :: pname = "zmat_diag_print_eigeninfo: "


C     == Variables ==
      integer i
      integer me
      integer dim1, dim2, dtype
      integer n


      
      me = ga_nodeid ()


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
C     Check the eigenvector matrix.
C
      call ga_check_handle (g_zevecs,
     $     "second argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zevecs, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)
      if (dim1.ne.n)
     $     call errquit (pname//"size of evecs must match size of mat")



C
C     Print info to stdout
C
c$$$      if (me.eq.0) then 
c$$$         write (luout, *) ""
c$$$         call util_print_centered (luout, "Eigensystem information",
c$$$     $        40, .true.)
c$$$      endif
c$$$
c$$$
c$$$      if (me.eq.0) then 
c$$$         write (luout, *) ""
c$$$         call util_print_centered (luout, "Original matrix", 0, .true.)
c$$$      endif
c$$$      call ga_print (g_za)

      if (me.eq.0) then 
         write (luout, *) ""
         call util_print_centered (luout, "Eigenvalues", 0, .true.)
      endif

C     (eigenvalues and abs eigenvalues)
      do i = 1, n
         if (me.eq.0)
     $        write (luout,"(i8,a,1es16.8e3,a,1es16.8e3,a,1es16.8e3,a)")
     $        i, ": " , dble (zevals(i)), " + ",
     $        aimag (zevals(i)), " i   (", abs(zevals(i)), " )"
         call util_flush (luout)
      enddo

C     (eigenvectors)
c$$$      if (me.eq.0) then 
c$$$         write (luout, *) ""
c$$$         call util_print_centered (luout, "Eigenvector matrix",0,.true.)
c$$$      endif
c$$$      call ga_print (g_zevecs)
c$$$
c$$$      if (me.eq.0) write (luout, *) ""
      if (me.eq.0) write (luout, *) ""

      end subroutine




C====================================================================
      subroutine zmat_diag_example_driver ()
      implicit none

#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"


C     == Parameters ==
      character(*), parameter :: pname = "zmat_diag_example_driver: "


C     == Variables ==
      integer me
      integer g_zevecs, g_zmat
      integer lvals, ivals
      integer i
      double complex val


      me = ga_nodeid ()


C      
C     EXAMPLE ONE
C
C     octave:3> A=[[1,2,0];[0,3,0];[2,-4,2]]
C     A =
C     
C     1   2   0
C     0   3   0
C     2  -4   2
C     
C     octave:4> [vec,val]=eig(A)
C     vec =
C     
C     0.00000   0.44721   0.40825
C     0.00000   0.00000   0.40825
C     1.00000  -0.89443  -0.81650
C     
C     val =
C     
C     Diagonal Matrix
C     
C     2   0   0
C     0   1   0
C     0   0   3
C     
C     
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout, "Diagonalization example one",
     $        20, .true.)
      endif
      

C
C     Allocation.
C
      if (.not. ga_create(mt_dcpl,3,3,"matrix" , 0, 0, g_zmat))
     $     call errquit ("failed to create mat", 0, 0)
      
      if (.not. ga_create(mt_dcpl,3,3,"eigenvectors" , 0, 0, g_zevecs))
     $     call errquit ("failed to create vecs", 0, 0)
      
      if (.not.ma_push_get(mt_dcpl,3,'vals',lvals,ivals))
     &     call errquit(pname//'cannot allocate vals',0, MA_ERR)

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
C     Diagonalization
C
      call zmat_diag(g_zmat, g_zevecs, dcpl_mb(ivals), .false., .false.)
      call zmat_diag_check (g_zmat, g_zevecs, dcpl_mb(ivals), 1d-6)
      call zmat_diag_print_eigeninfo (g_zmat, g_zevecs, dcpl_mb(ivals))



      
C
C     EXAMPLE TWO
C
C     octave:5> A=[[1.2+i,1,i];[0,3,-0.98*i];[2+3.2*i,-4,2.5+8*i]]
C     A =
C     
C     1.20000 + 1.00000i   1.00000 + 0.00000i   0.00000 + 1.00000i
C     0.00000 + 0.00000i   3.00000 + 0.00000i  -0.00000 - 0.98000i
C     2.00000 + 3.20000i  -4.00000 + 0.00000i   2.50000 + 8.00000i
C     
C     octave:6> [vec,val]=eig(A)
C     vec =
C     
C     0.12030 + 0.04509i   0.82804 + 0.00000i   0.60683 + 0.06457i
C     -0.11377 - 0.00080i  -0.05384 - 0.23655i   0.70499 + 0.00000i
C     0.98517 + 0.00000i  -0.48480 + 0.14304i  -0.27740 - 0.23160i
C     
C     val =
C     
C     Diagonal Matrix
C
C     3.05971 + 8.48554i                    0                    0
C                      0   0.96224 + 0.12885i                    0
C                      0                    0   2.67805 + 0.38561i
C
      if (me.eq.0) then
         write (luout, *) ""
         write (luout, *) ""
         call util_print_centered (luout, "Diagonalization example two",
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
C     Diagonalization
C
      call zmat_diag(g_zmat, g_zevecs, dcpl_mb(ivals), .false., .false.)
      call zmat_diag_check (g_zmat, g_zevecs, dcpl_mb(ivals), 1d-6)
      call zmat_diag_print_eigeninfo (g_zmat, g_zevecs, dcpl_mb(ivals))
      


C
C     Clean up.
C     
      if (.not. ga_destroy (g_zmat))
     $     call errquit (pname//"failed to destroy g_zmat",0,0)
      if (.not. ga_destroy (g_zevecs))
     $     call errquit (pname//"failed to destroy g_zmat",0,0)
      if (.not. ma_pop_stack (lvals))
     $     call errquit (pname//"failed to pop lvals",0,0)

      end subroutine


      
c $Id$
