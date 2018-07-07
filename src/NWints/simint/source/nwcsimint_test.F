      subroutine nwcsim_f90test()
      use SimintFortran
      implicit none
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_shell), target :: sh(2), zero_sh(1,3)
      type(c_simint_multi_shellpair), target :: bra_msh
      type(c_simint_multi_shellpair), target :: ket_msh

      double precision, pointer :: p1(:), p2(:)
      double precision :: alpha(3), coef(3)
      double precision :: integrals(1000)
      double precision, allocatable :: work(:)
      integer :: i, ncomputed
      integer :: worksize
      integer :: num_der, max_ang
      integer l_wrk,k_wrk
      integer*8 smnt_screen_method
      parameter(smnt_screen_method=2)
      integer*4 smnt_deriv
      parameter(smnt_deriv=1)

      call simint_init()

      num_der=0
#ifdef SIMINT_GRADIENT
      num_der=1
#endif
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
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)
      
      call simint_create_shell(3, 0, 0.0d0, 0.0d0, 0.0d0, 
     &                          alpha, coef, sh(1)) 
      alpha(1)=0.20d0
      coef(1)=1d0
      call simint_create_shell(1, 0, 0.0d0, 0.0d0, 0.0d0, 
     &                          alpha, coef, sh(2)) 
      call simint_normalize_shells(2, sh)
      
      call C_F_POINTER(sh(1)%alpha, p1, shape=[sh(1)%nprim])
      call C_F_POINTER(sh(1)%coef, p2, shape=[sh(1)%nprim])
      
      call simint_create_multi_shellpair(2, sh, 2, sh, bra_msh, 0)
      call simint_create_multi_shellpair(2, sh, 2, sh, ket_msh, 0) 
      
      call C_F_POINTER(bra_msh%alpha, p1, shape=[bra_msh%nprim])
      call C_F_POINTER(bra_msh%prefac, p2, shape=[bra_msh%nprim])
      
      write(*,*) "Shell Pair info"
      do i = 1, bra_msh%nprim
        write(*,*) p1(i), p2(i)
      end do
      
#ifdef SIMINT_GRADIENT
      write(6,*) ' using derivatives 1 code '
      ncomputed = simint_compute_eri_deriv(smnt_deriv,bra_msh, ket_msh,
     s     1d-15,
     C     dbl_mb(k_wrk), integrals)
#else
      ncomputed = simint_compute_eri(bra_msh, ket_msh, 1d-15,
     C     dbl_mb(k_wrk), integrals)
#endif
      write(*,*) ' 2e4c', integrals(1),integrals(ncomputed)
c     try 2e3c
      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
#if 0
      call simint_create_zero_shell(zero_sh(1,1)) 
#else
      alpha(1)=0d0
      coef(1)=1d0
      call simint_create_shell(1, 0, 0.0d0, 0.0d0, 0.0d0, 
     &                          alpha, coef, zero_sh(1,1)) 
#endif
      call C_F_POINTER(zero_sh(1,1)%alpha, p1, 
     S     shape=[zero_sh(1,1)%nprim])
      call C_F_POINTER(zero_sh(1,1)%coef, p2, 
     S     shape=[zero_sh(1,1)%nprim])
      
      write(*,*) "Shell info"
      do i = 1, zero_sh(1,1)%nprim
        write(*,*) p1(i), p2(i)
      end do

      call simint_create_multi_shellpair(1, sh, 1, zero_sh(1,1), 
     B     bra_msh, 0)
      call simint_create_multi_shellpair(1, sh, 1, sh, ket_msh, 0)

      ncomputed = simint_compute_eri(bra_msh, ket_msh, 0.0d0, 
     C     dbl_mb(k_wrk), integrals)

      write(*,*) ' 2e3c',integrals(1),integrals(ncomputed)
      call C_F_POINTER(bra_msh%alpha, p1, shape=[bra_msh%nprim])
      call C_F_POINTER(bra_msh%prefac, p2, shape=[bra_msh%nprim])
      
      write(*,*) "Shell Pair info"
      do i = 1, bra_msh%nprim
        write(*,*) p1(i), p2(i)
      end do
c     try 2e2c
      call simint_free_multi_shellpair(ket_msh)
      call simint_create_multi_shellpair(1, sh, 1, zero_sh(1,1), 
     B     ket_msh, 0)
      ncomputed = simint_compute_eri(bra_msh, ket_msh, 0.0d0, 
     C     dbl_mb(k_wrk), integrals)
      write(*,*) ' 2e2c',integrals(1),integrals(ncomputed)
      
      
c      deallocate(work)
      if(.not.ma_pop_stack(l_wrk)) call errquit(
     E     ' pop stack failed ',l_wrk,0)
      call simint_free_shell(sh(1))
      call simint_free_shell(sh(2))
      call simint_free_shell(zero_sh(1,1))
      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
      
      call simint_finalize()
      
      
      return
      end
      
