c
c     Compute unrelaxed difference density matrix T
c     HvD 9/2006, NG 11/2012   
c
      subroutine tddft_grad_compute_t(tda,ipol,naoc,nav,nroot,
     +      g_xpy,g_xmy,g_tp)
c
c     Eq. 19: Furche & Ahlrich JCP (2002)
c     Tijs = -1/2 Sum_a [(X+Y)ias(X+Y)jas+(X-Y)ias(X-Y)jas]
c     Tabs =  1/2 Sum_i [(X+Y)ias(X+Y)ibs+(X-Y)ias(X-Y)ibs]
c     Tias = Tais = 0 (but we won't actually initialise these)
c
      implicit none
c
#include "global.fh"
#include "tddft_grad_util.fh"
#include "stdio.fh"
c
c     Input:
c
      logical tda       ! .true. if Tamm-Dancoff approximation is used
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer naoc(2)   ! the number of active occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nroot     ! the number of states to consider
      integer g_xpy(2)  ! global array handle for (X+Y)
      integer g_xmy(2)  ! global array handle for (X-Y)
c
c     Output:
c
      integer g_tp(2)   ! global array handle for T
c
c     Local:
c
      integer alo(3)    ! lower chunk limits on A
      integer ahi(3)    ! upper chunk limits on A
      integer blo(3)    ! lower chunk limits on B
      integer bhi(3)    ! upper chunk limits on B
      integer clo(3)    ! lower chunk limits on C
      integer chi(3)    ! upper chunk limits on C
      integer ip        ! counter on polarization
      integer ir        ! counter on roots
      logical oroot
c
      character*32 pname
      pname="tddft_grad_compute_t: "
c
c     Loop over spins...
c
      do ip = 1, ipol
c
        call ga_zero(g_tp(ip))
c
c       Do the occupied-occupied block
c
        do ir = 1, nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = naoc(ip)
          alo(3) = 1
          ahi(3) = nav(ip)
          blo(1) = ir
          bhi(1) = ir
          blo(2) = 1
          bhi(2) = nav(ip)
          blo(3) = 1
          bhi(3) = naoc(ip)
          clo(1) = ir
          chi(1) = ir
          clo(2) = 1
          chi(2) = naoc(ip)
          clo(3) = 1
          chi(3) = naoc(ip)
c Daniel (12-3-12): Rather than calling nga_matmul_patch followed by
c nga_scale_patch, it is probably better to check whether the
c calculation is CIS first.
          if (tda) then
c For CIS, we have:
c (X+Y)ias = Xias , (X+Y)jas = Xjas
c (X-Y)ias = Xias , (X-Y)jas = Xjas
c Tij = - sum_a Xias*Xjas 
            call nga_matmul_patch('n','t',-0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
            call nga_matmul_patch('n','t',-0.5d0,1.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          else
c (X+Y)ias(X+Y)jas 
            call nga_matmul_patch('n','t',-0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
c (X-Y)ias(X-Y)jas
            call nga_matmul_patch('n','t',-0.5d0,1.0d0,
     +                            g_xmy(ip),alo,ahi,
     +                            g_xmy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          endif ! tda
        enddo  ! ir = 1, nroot
c
c       Do the virtual-virtual block
c
        do ir = 1, nroot
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = nav(ip)
          alo(3) = 1
          ahi(3) = naoc(ip)
          blo(1) = ir
          bhi(1) = ir
          blo(2) = 1
          bhi(2) = naoc(ip)
          blo(3) = 1
          bhi(3) = nav(ip)
          clo(1) = ir
          chi(1) = ir
          clo(2) = naoc(ip)+1
          chi(2) = naoc(ip)+nav(ip)
          clo(3) = naoc(ip)+1
          chi(3) = naoc(ip)+nav(ip)
c Daniel (12-3-12): Rather than calling nga_matmul_patch followed by
c nga_scale_patch, it is probably better to check whether the
c calculation is CIS first.
          if (tda) then
c For CIS, we have:
c (X+Y)ias = Xias , (X+Y)ibs = Xibs
c (X-Y)ias = Xias , (X-Y)ibs = Xibs
c Tab = sum_i Xias*Xibs 
            call nga_matmul_patch('t','n',0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
            call nga_matmul_patch('t','n',0.5d0,1.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          else
c (X+Y)ias(X+Y)ibs 
            call nga_matmul_patch('t','n',0.5d0,0.0d0,
     +                            g_xpy(ip),alo,ahi,
     +                            g_xpy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
c (X-Y)ias(X-Y)ibs 
            call nga_matmul_patch('t','n',0.5d0,1.0d0,
     +                            g_xmy(ip),alo,ahi,
     +                            g_xmy(ip),blo,bhi,
     +                            g_tp(ip),clo,chi)
          endif  ! tda
        enddo  ! ir = 1,nroot
c
      enddo ! ip = 1, ipol
c
      if (tddft_grad_util_print('tddft grad t',print_debug)) then
        oroot = ga_nodeid().eq.0
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'T'
        call tddft_grad_print_array(ipol,nroot,g_tp,dble(ipol))
      endif
c
      end
c $Id$
