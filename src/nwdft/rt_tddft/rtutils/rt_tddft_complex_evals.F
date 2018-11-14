C
C     Compute and print complex eigenvalues of a potentially
C     non-hermitian matrix.
C
      subroutine rt_tddft_complex_evals (params,tt, g_za, name)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "matutils.fh"      
#include "rtdb.fh"


      
C     == Inputs ==
      type (rt_params_t), intent(in) :: params
      integer, intent(in)            :: g_za !matrix to diagonalize and print eigenvalues for
      character(len=*), intent(in)   :: name
      double precision , intent(in)            :: tt
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_complex_evals: "


C     == Variables ==
      integer ivals, lvals
      integer g_zevecs
      integer dtype, dim1, dim2
      integer i
      double complex zval, eval_print
c      double complex zocc
      integer j, n
!      double precision tt       ! time
!      integer it
C
C     Check the matrix and create arrays.
C
      call ga_check_handle (g_za,
     $     "first argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_za, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)

      if (.not.ma_push_get(mt_dcpl, dim1,'vals',lvals,ivals))
     &     call errquit(pname//'cannot allocate vals',0, MA_ERR)   !xxxx size???

      if (.not. ga_duplicate (g_za, g_zevecs, "zevecs"))
     $     call errquit (pname//"failed to create zevecs", 0, GA_ERR)


C
C     Diagonalize and print to stdout
C
      call zmat_diag (g_za, g_zevecs, dcpl_mb(ivals), .true., .false.) !true = reorder evals by energy, false = non hermitian

C      call ga_print (dcpl_mb(ivals))
!      tt = params%tmin + (it-1)*params%dt
!      tt = params%nt  

!      do it = 1, params%nt
!         tt = params%tmin + (it-1)*params%dt

!     enddo

      
      if (ga_nodeid().eq.0) then
         write(luout, *) ""
         call util_print_centered (luout,
     $        "Complex eigenvalues of "//trim(name), 40, .true.)
         
         do i = 1, dim1
            zval = dcpl_mb(ivals + i - 1)
C     call ga_get (g_za, i, i, i, i, zocc, 1)  !on diagonal part (occup for MO dens matrix)
C     write(luout, "(1x,a,i8,3x, 2e16.6)")
C     $           "Vector", i, tt,  dble(zval), aimag(zval)
            write(luout, *)
     $           "Vector", i, tt,  dble(zval), aimag(zval)
         enddo
         
      endif
      
      call ga_sync()
      
C     
C     Clean up
C

      if (.not. ma_chop_stack (lvals))
     &   call errquit(pname//'cannot chop stack',0, MA_ERR)

      if (.not. ga_destroy (g_zevecs))
     $     call errquit (pname//"failed to destrory zevecs", 0, GA_ERR)

      end subroutine



C
C     Compare the abs diff of complex eigenvalues between two matricies
C     which are potentially non-Hermitian.
C
      subroutine rt_tddft_complex_evals_compare (params,
     $     g_za, g_zb, diff, imax)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "matutils.fh"      
 
C     == Inputs ==
      type (rt_params_t), intent(in) :: params
      integer, intent(in)            :: g_za
      integer, intent(in)            :: g_zb
      double precision               :: diff  ! the max abs diff between two elements
      integer, intent(out)           :: imax  ! which eigenvalue has max diff

      
C     == Parameters ==
      character(*), parameter :: pname =
     $     "rt_tddft_complex_evals_compare: "


C     == Variables ==
      integer ivals_a, lvals_a
      integer ivals_b, lvals_b
      integer g_zevecs_a, g_zevecs_b
      integer dtype, dim1, dim2
      integer i
      double complex zval
      integer me
      integer msglen
      double precision :: diff_one, diff_max
      integer imax_local
c      double complex zocc

      me = ga_nodeid()
      

C
C     Check the matrix and create arrays.
C
      call ga_check_handle (g_za,
     $     "first argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_za, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)

      if (.not.ma_push_get(mt_dcpl, dim1,'vals_a',lvals_a,ivals_a))
     &     call errquit(pname//'cannot allocate vals_a',0, MA_ERR)   !xxxx size???

      if (.not. ga_duplicate (g_za, g_zevecs_a, "zevecs_a"))
     $     call errquit (pname//"failed to create zevecs_a", 0, GA_ERR)

      call ga_check_handle (g_zb,
     $     "first argument of "//pname//"() is not a valid GA")
      call ga_inquire (g_zb, dtype, dim1, dim2)
      if (dtype .ne. mt_dcpl) call errquit (pname//
     $     "expecting complex-valued GA as second argument", 0, 0)
      if (dim1 .ne. dim2)
     $     call errquit (pname//"dim1 must equal dim2", 0, 0)

      if (.not.ma_push_get(mt_dcpl, dim1,'vals_b',lvals_b,ivals_b))
     &     call errquit(pname//'cannot allocate vals_b',0, MA_ERR)   !xxxx size???

      if (.not. ga_duplicate (g_zb, g_zevecs_b, "zevecs_b"))
     $     call errquit (pname//"failed to create zevecs_b", 0, GA_ERR)

C     XXX should check that A and B are consistent

C
C     Diagonalize and print to stdout
C
C      call zmat_diag (g_za, g_zevecs_a, dcpl_mb(ivals_a),.true.,.false.) !true = reorder evals_a by energy, false = non hermitian
C      call zmat_diag (g_zb, g_zevecs_b, dcpl_mb(ivals_b),.true.,.false.) !true = reorder evals_a by energy, false = non hermitian

C      call ga_print (dcpl_mb(ivals_a))
C      call ga_print (dcpl_mb(ivals_b))
C      call ga_print ( dcpl_mb(ivals_b+i-1))
C
C     Print eigenvalues to screen (messy)
C
      if (me.eq.0) then
         do i = 1, dim1
            write(luout,*) "EVAL ",
     $           i, dcpl_mb(ivals_a+i-1), dcpl_mb(ivals_b+i-1) 
         enddo
      endif
      
c$$$
c$$$
c$$$      write(luout,*) "INSIDE COMPLEX EVALS a: ",
c$$$     $     me, dcpl_mb(ivals_a + 0),
c$$$     $     dcpl_mb(ivals_a + 9), dcpl_mb(ivals_a + 21)
c$$$
c$$$      write(luout,*) "INSIDE COMPLEX EVALS b: ",
c$$$     $     me, dcpl_mb(ivals_b + 0),
c$$$     $     dcpl_mb(ivals_b + 9), dcpl_mb(ivals_b + 21)

      
C
C     Compute sum abs diff of eigenvalues
C
      call ga_sync()
      if (ga_nodeid().eq.0) then
         diff_max = 0d0
         diff = 0d0
         imax_local = -99
         do i = 1, dim1
            diff_one = abs(dcpl_mb(ivals_a + i - 1) -
     $           dcpl_mb(ivals_b + i - 1))
            write(luout,*) "DIFFS: ", i, diff_one

            if (diff_one .gt. diff_max) then 
               diff_max = diff_one
               imax_local = i
C               write(luout,*) "imax_local is now ", imax_local
            endif
            
C            diff = diff + diff_one  !old way: sum of diffs
         enddo
         diff = diff_max
         imax = imax_local
      endif
      
      call ga_sync()
      msglen = MA_sizeof(MT_DBL,1,MT_BYTE)
      call ga_brdcst(100001, diff, msglen, 0)
      
      call ga_sync()
      msglen = MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_brdcst(100002, imax, msglen, 0)
      
C
C     Clean up
C
      if (.not. ma_chop_stack (lvals_a))
     &   call errquit(pname//'cannot chop stack',0, MA_ERR)

      if (.not. ga_destroy (g_zevecs_a))
     $     call errquit (pname//"failed to destrory zevecs_a",0,GA_ERR)
      if (.not. ga_destroy (g_zevecs_b))
     $     call errquit (pname//"failed to destrory zevecs_b",0,GA_ERR)

      end subroutine
