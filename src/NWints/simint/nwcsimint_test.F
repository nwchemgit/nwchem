      subroutine nwcsim_f90test()
      use SimintFortran
      implicit none
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_shell), target :: sh(2)
      type(c_simint_multi_shellpair), target :: msh

      double precision, pointer :: p1(:), p2(:)
      double precision :: alpha(3), coef(3)
      double precision :: integrals(1000)
      double precision, allocatable :: work(:)
      integer :: i, ncomputed
      integer :: worksize
      integer :: num_der, max_ang
      integer l_wrk,k_wrk

      call simint_init()

      num_der=0
      max_ang=1

      worksize = simint_eri_worksize(num_der, max_ang)

c      allocate(work(worksize))
      if (.not. MA_Push_Get(MT_Dbl,worksize,'simintw',l_wrk,k_wrk))
     C     call errquit(' simint wrk ',worksize,0)


      alpha(1) = 130.7093200d0
      coef(1) = 0.15432897d0
      alpha(2) = 23.8088610d0
      coef(2) = 0.53532814d0
      alpha(3) = 6.4436083d0
      coef(3) = 0.44463454d0
      
      call simint_initialize_shell(sh(1))
      call simint_initialize_shell(sh(2))
      call simint_initialize_multi_shellpair(msh)
      
      call simint_create_shell(3, 0, 0.0d0, 0.0d0, 0.0d0, 
     &                          alpha, coef, sh(1)) 
      alpha(1)=0.20d0
      coef(1)=1d0
      call simint_create_shell(1, 0, 0.0d0, 0.0d0, 0.0d0, 
     &                          alpha, coef, sh(2)) 
      call simint_normalize_shells(2, sh)
      
      call C_F_POINTER(sh(1)%alpha, p1, shape=[sh(1)%nprim])
      call C_F_POINTER(sh(1)%coef, p2, shape=[sh(1)%nprim])
      
      call simint_create_multi_shellpair(2, sh, 2, sh, msh, 0)
      
      call C_F_POINTER(msh%alpha, p1, shape=[msh%nprim])
      call C_F_POINTER(msh%prefac, p2, shape=[msh%nprim])
      
      write(*,*) "Shell Pair info"
      do i = 1, msh%nprim
        write(*,*) p1(i), p2(i)
      end do
      
      ncomputed = simint_compute_eri(msh, msh, 0.0d0, 
     C     dbl_mb(k_wrk), integrals)
      write(*,*) integrals(1),integrals(ncomputed)
      
c      deallocate(work)
      if(.not.ma_pop_stack(l_wrk)) call errquit(
     E     ' pop stack failed ',l_wrk,0)
      call simint_free_shell(sh(1))
      call simint_free_shell(sh(2))
      call simint_free_multi_shellpair(msh)
      
      call simint_finalize()
      
      
      return
      end
      
