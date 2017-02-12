! In any Gaussian based QM code the performance of the 2-electron
! integrals is likely to be critical (unless the wave function is
! represented on a grid). Recent reviews on integral evaluations is
! [1] and [2].
! Wim Klopper suggested that the evaluation of the integrals could be
! performed faster by approximating the natural exponent functions (EXP)
! as a Taylor series [3]. To see to what extent using a Taylor series
! may improve the performance this kernel was written. Given the order
! of the Taylor series all other quantities are automatically
! initialized. We also use the fact that we are interested only in
! exp(x) for ln(eps) <= x <= 0 (where eps is the machine precision).
! For x < ln(eps) we assume exp(x)=0.
! If nord = 4 then the Taylor series approximation is 1.5 times faster
! than the exp function on a Intel(R) Xeon(R) CPU E5-2670 0 @ 2.60GHz.
! The interpolation table requires 1.4 MB of memory. The memory
! requirements mean that the algorithm could be badly impacted by cache
! misses. Knights Landing does not have such big caches. Also the Volta
! GPU is unlikely going to have such big caches. In addition Knights
! Landing has added AVX-512ER instructions to perform exponentiations
! of the form 2^x. These might help boost the performance of the exp
! function well beyond what a Taylor series brings. These instructions
! are not available on all Intel processors though.
!
! [1] Simen Reine, Trygve Helgaker, Roland Lindh, "Multi-electron
!     integrals", WIREs Comput. Mol. Sci. 2012, 2:290-303,
!     doi:10.1002/wcms.78
! [2] Justin T. Fermann, Edward F. Valeev, "Fundamentals of Molecular
!     Integrals Evaluation", 
!     https://www.valeevgroup.chem.vt.edu/docs/ints.07032015.pdf
! [3] Wim Klopper, "SCF methods, basis sets, and integrals:
!     Lecture IV: Integrals", ESQC 2011, Torre Normanna, 9/18 –
!     10/01/2011,
!     http://esqc.ups-tlse.fr/11/lectures/KlopperLecture4-1x2.pdf
!     Slide: "Calculating the Boys function".
!
      module taylor
      double precision, allocatable :: ts(:) ! table of exp(x_i)
      double precision :: h  ! step size
      double precision :: xm ! x_min (exp(x) is needed only in the range
                             ! x_min <= x <= 0)
      double precision :: xmh ! = xm+h/2
      integer, parameter :: nord = 4 ! the maximum order
      double precision, parameter :: eps = 1.0d-16 ! machine precision
      contains
      subroutine t_init
        ! Initialize the interpolation table
        implicit none
        integer :: npt ! number of points
        integer :: ipt ! current point
        double precision :: a, x
        xm  = dlog(eps)
        h   = 2.0d0*(eps**(1.0d0/dble(nord)))
        xmh = xm+0.5d0*h
        npt = (-xm/h)+1
        allocate(ts(npt))
        do ipt = 1, npt
          x=xmh+h*(ipt-1)
          ts(ipt) = dexp(x)
        enddo
        write(*,*)'xm ',xm
        write(*,*)'h  ',h
        write(*,*)'npt',npt,' requires ',dble(8*npt)/1024,' KB'
      end subroutine t_init
      double precision function t_exp(x)
        ! Use Taylor series to evaluate exp(x)
        implicit none
        double precision, intent(in) :: x
        double precision :: xi, xd
        double precision :: r
        integer ipt
        integer iord
        ipt = (x-xmh)/h
        xi  = xmh+h*(ipt-1)
        xd  = x-xi
        r   = xd/dble(nord)+1.0d0
        do iord = nord-1, 1, -1
          r = r*xd/dble(iord)+1.0d0
        enddo
        t_exp = r*ts(ipt)
      end function t_exp
      end module taylor

      program performance
      use taylor
      implicit none
#include "mpif.h"
      integer :: ierr, nmax, i
      double precision :: a, b
      double precision :: tt, te, t0
      nmax = 10000000
      a = -30.0d0
      call mpi_init(ierr)
      call t_init
      t0 = -mpi_wtime()
      b =   0.0d0
      do i = 1, nmax
        b = b + a/i
      enddo
      t0 = t0 + mpi_wtime()
      write(*,*)'Null time:  ',t0,b
      tt = -mpi_wtime()
      b =   0.0d0
      do i = 1, nmax
        b = b + t_exp(a/i)
      enddo
      tt = tt + mpi_wtime()
      write(*,*)'t-exp time: ',tt-t0,b
      te = -mpi_wtime()
      b =   0.0d0
      do i = 1, nmax
        b = b + exp(a/i)
      enddo
      te = te + mpi_wtime()
      write(*,*)'Exp time:   ',te-t0,b
      write(*,*)'Exp/t-exp:  ',(te-t0)/(tt-t0)
!     do i = 1, 1000
!       write(*,*)a/i,exp(a/i),t_exp(a/i)
!     enddo
      call mpi_finalize(ierr)
      end

      


