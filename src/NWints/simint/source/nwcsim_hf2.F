#ifdef USE_SIMINT
      subroutine nwcsim_hf2(ij_basis,i_sh,j_sh,
     J     kl_basis,k_sh,l_sh,
     O     ncomputed, eri,leri,scr,lscr)
c      use SimintFortran
      use nwcsim_shell
      implicit none
      integer ij_basis        ! basis for bra [in]
      integer kl_basis       ! basis for ket [in]
      integer i_sh,j_sh,k_sh,l_sh
      integer ncomputed
      double precision eri(*)
      integer leri
      double precision scr(*)
      integer lscr
#include "errquit.fh"
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_multi_shellpair), target :: bra_msh, ket_msh
#define DEBUG_ 1
#ifdef DEBUG
      double precision, pointer :: p1(:), p2(:),p3(:),p4(:)
#endif
      integer i
      character*10 pname
      integer*4 zerod
      parameter(zerod=0)
      pname="nwcsim_hf2"
c
c      match  ijkl shell to simint shell

c
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)
      call simint_create_multi_shellpair(1, smnt_sh(i_sh,ij_basis), 1, 
     S     smnt_sh(j_sh,ij_basis), bra_msh, smnt_screen_method)
      call simint_create_multi_shellpair(1,smnt_sh(k_sh,kl_basis),
     S     1,smnt_sh(l_sh,kl_basis), ket_msh, smnt_screen_method)
#ifdef DEBUG
      write(6,'(a,2i10)') ' basis tags ',ij_basis,kl_basis
      write(6,'(4i10)') i_sh,j_sh,k_sh,l_sh
      call C_F_POINTER(bra_msh%alpha, p1, shape=[bra_msh%nprim])
      call C_F_POINTER(bra_msh%prefac, p2, shape=[bra_msh%nprim])
      call C_F_POINTER(bra_msh%ab_x, p3, shape=[bra_msh%nprim])
      call C_F_POINTER(bra_msh%pa_x, p4, shape=[bra_msh%nprim])
      
      write(*,*) "Shell Pair info"
      do i = 1, bra_msh%nprim
        write(*,*) p1(i), p2(i),p3(i),p4(i)
      end do
#endif

#ifdef SIMINT_GRADIENT
      ncomputed = simint_compute_eri_deriv(zerod,bra_msh, ket_msh,
     c     smnt_screen_tol,
     C     scr, eri)
#else
      ncomputed = simint_compute_eri(bra_msh, ket_msh,
     C     smnt_screen_tol,
     C     scr, eri)
#endif

      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
#ifdef DEBUG
      write(6,*) ncomputed,' eri ', eri(1),eri(ncomputed)
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
#endif
      if(ncomputed.le.0) call errquit(
     P     pname//'return value from simint ',ncomputed, BASIS_ERR)
cwrong      if(ncomputed.gt.leri) call errquit(
c     P     pname//'ncomputed gt leri ',ncomputed*10000+leri, BASIS_ERR)
 
      return
      end
      subroutine nwcsim_hf2_3c(i_basis,i_sh,jk_basis,j_sh,k_sh,
     O     ncomputed, eri,leri,scr,lscr)
c      use SimintFortran
      use nwcsim_shell
      implicit none
      integer i_basis        ! basis for bra [in]
      integer jk_basis       ! basis for ket [in]
      integer i_sh,j_sh,k_sh ! [in]
      integer ncomputed
      double precision eri(*)
      integer leri
      double precision scr(*)
      integer lscr
#include "errquit.fh"
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_multi_shellpair), target :: bra_msh, ket_msh
      integer i
       character*10 pname
      pname="nwcsim_hf2_3c"
c
c      match  ijkl shell to simint shell

c
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)

      call simint_create_multi_shellpair(1,smnt_sh(i_sh,i_basis),
     S     1,zero_sh(j_sh,i_basis), bra_msh, smnt_screen_method)
      call simint_create_multi_shellpair(1, smnt_sh(j_sh,jk_basis), 1, 
     S     smnt_sh(k_sh,jk_basis), ket_msh, smnt_screen_method)
      ncomputed = simint_compute_eri(bra_msh, ket_msh, smnt_screen_tol, 
     C     scr, eri)
      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
#ifdef DEBUG
      write(6,*) ' eri ', eri(1),eri(ncomputed)
      if (.not.ma_verify_allocator_stuff()) stop ' ma ck 2b'
#endif
      if(ncomputed.gt.leri) call errquit(
     P     pname//'ncomputed gt leri ',ncomputed*10000+leri, BASIS_ERR)
 
      return
      end
      subroutine nwcsim_hf2_2c(i_basis,i_sh,j_basis,j_sh,
     O     ncomputed, eri,leri,scr,lscr)
c      use SimintFortran
      use nwcsim_shell
      implicit none
      integer i_basis        ! basis for bra [in]
      integer j_basis       ! basis for ket [in]
      integer i_sh,j_sh ! [in]
      integer ncomputed
      double precision eri(*)
      integer leri
      double precision scr(*)
      integer lscr
#include "errquit.fh"
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_multi_shellpair), target :: bra_msh, ket_msh
#define DEBUG_ 1
      integer i
       character*10 pname
      pname="nwcsim_hf2_2c"
c
c      match  ijkl shell to simint shell

c
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)

      call simint_create_multi_shellpair(1,smnt_sh(i_sh,i_basis),
     S     1,zero_sh(i_sh,i_basis), bra_msh, smnt_screen_method)
      call simint_create_multi_shellpair(1, smnt_sh(j_sh,j_basis), 1, 
     S     zero_sh(j_sh,j_basis), ket_msh, smnt_screen_method)
      ncomputed = simint_compute_eri(bra_msh, ket_msh,
     S     smnt_screen_tol,
     C     scr, eri)
      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
      if(ncomputed.le.0) call errquit(
     P     pname//'return value from simint ',ncomputed, BASIS_ERR)
c      if(ncomputed.gt.leri) call errquit(
c     P     pname//'ncomputed gt leri ',ncomputed*10000+leri, BASIS_ERR)
 
      return
      end
#ifdef SIMINT_GRADIENT
      subroutine nwcsim_hf2d(ij_basis,i_sh,j_sh,
     J     kl_basis,k_sh,l_sh,
     O     ncomputed, eri,leri,scr,lscr)
c      use SimintFortran
      use nwcsim_shell
      implicit none
      integer ij_basis        ! basis for bra [in]
      integer kl_basis       ! basis for ket [in]
      integer i_sh,j_sh,k_sh,l_sh
      integer ncomputed
      double precision eri(*)
      integer leri
      double precision scr(*)
      integer lscr
#include "errquit.fh"
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_multi_shellpair), target :: bra_msh, ket_msh
      integer i
       character*10 pname
      pname="nwcsim_hf2d"
c
c      match  ijkl shell to simint shell

c
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)

      call simint_create_multi_shellpair(1, smnt_sh(i_sh,ij_basis), 1, 
     S     smnt_sh(j_sh,ij_basis), bra_msh, smnt_screen_method)
      call simint_create_multi_shellpair(1,smnt_sh(k_sh,kl_basis),
     S     1,smnt_sh(l_sh,kl_basis), ket_msh, smnt_screen_method)

      ncomputed = simint_compute_eri_deriv(smnt_deriv,bra_msh, ket_msh,
     S     smnt_screen_tol,
     C     scr, eri)
      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
      if(ncomputed.gt.leri) call errquit(
     P     pname//'ncomputed gt leri ',ncomputed*10000+leri, BASIS_ERR)
 
      return
      end
      subroutine nwcsim_hf2d_2c(i_basis,i_sh,
     J     j_basis,j_sh,
     O     ncomputed, eri,leri,scr,lscr)
c      use SimintFortran
      use nwcsim_shell
      implicit none
      integer i_basis        ! basis for bra [in]
      integer j_basis       ! basis for ket [in]
      integer i_sh,j_sh
      integer ncomputed
      double precision eri(*)
      integer leri
      double precision scr(*)
      integer lscr
#include "errquit.fh"
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_multi_shellpair), target :: bra_msh, ket_msh
      integer i
       character*10 pname
      pname="nwcsim_hf2d_2c"
c
c      match  ijkl shell to simint shell

c
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)

      call simint_create_multi_shellpair(1, smnt_sh(i_sh,i_basis), 1, 
     S     zero_sh(i_sh,i_basis), bra_msh, smnt_screen_method)
      call simint_create_multi_shellpair(1, smnt_sh(j_sh,j_basis), 1,
     S     zero_sh(j_sh,j_basis), ket_msh, smnt_screen_method)

      ncomputed = simint_compute_eri_deriv(smnt_deriv,bra_msh, ket_msh,
     S     smnt_screen_tol,
     C     scr, eri)
      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
      if(ncomputed.gt.leri) call errquit(
     P     pname//'ncomputed gt leri ',ncomputed*10000+leri, BASIS_ERR)
 
      return
      end
      subroutine nwcsim_hf2d_3c(i_basis,i_sh,
     J     jk_basis,j_sh,k_sh,
     O     ncomputed, eri,leri,scr,lscr)
c      use SimintFortran
      use nwcsim_shell
      implicit none
      integer i_basis        ! basis for bra [in]
      integer jk_basis       ! basis for ket [in]
      integer i_sh,j_sh,k_sh
      integer ncomputed
      double precision eri(*)
      integer leri
      double precision scr(*)
      integer lscr
#include "errquit.fh"
#include "mafdecls.fh"
c      use iso_c_binding
      type(c_simint_multi_shellpair), target :: bra_msh, ket_msh
      integer i
       character*10 pname
      pname="nwcsim_hf2d_3c"
c
c      match  ijkl shell to simint shell

c
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)

      call simint_create_multi_shellpair(1, smnt_sh(i_sh,i_basis), 1, 
     S     zero_sh(i_sh,i_basis), bra_msh, smnt_screen_method)
      call simint_create_multi_shellpair(1,smnt_sh(j_sh,jk_basis), 1,
     S     smnt_sh(k_sh,jk_basis), ket_msh, smnt_screen_method)

      ncomputed = simint_compute_eri_deriv(smnt_deriv,bra_msh, ket_msh,
     s     smnt_screen_tol,
     C     scr, eri)
      call simint_free_multi_shellpair(bra_msh)
      call simint_free_multi_shellpair(ket_msh)
      if(ncomputed.gt.leri) call errquit(
     P     pname//'ncomputed gt leri ',ncomputed*10000+leri, BASIS_ERR)
 
      return
      end
#endif
#else
      subroutine nwcsim_hf2()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
      subroutine nwcsim_hf2_3c()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
      subroutine nwcsim_hf2_2c()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
      subroutine nwcsim_hf2d()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
      subroutine nwcsim_hf2d_3c()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
      subroutine nwcsim_hf2d_2c()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
#endif
