#ifdef USE_SIMINT
      subroutine nwcsim_hf2(i_sh,j_sh,k_sh,l_sh,
     O     ncomputed, eri,leri,scr,lscr)
c      use SimintFortran
      use nwcsim_shell
      implicit none
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
      pname="nwcsim_hf2"
c
c      match  ijkl shell to simint shell

c
      call simint_initialize_multi_shellpair(bra_msh)
      call simint_initialize_multi_shellpair(ket_msh)
      call simint_create_multi_shellpair(1,smnt_sh(i_sh),
     S     1,smnt_sh(j_sh), bra_msh, 0)
#ifdef DEBUG
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
      call simint_create_multi_shellpair(1, smnt_sh(k_sh), 1, 
     S     smnt_sh(l_sh), ket_msh, 0)
      ncomputed = simint_compute_eri(bra_msh, ket_msh, 0.0d0, 
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
#else
      subroutine nwcsim_hf2()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
#endif
