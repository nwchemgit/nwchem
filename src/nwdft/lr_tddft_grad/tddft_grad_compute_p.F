c
c     Add Zia to Tij+Tab to form P  (Eq. 23, Furche & Ahlrichs)
c     P = one-particle difference density matrix
c     HvD 10/2007, NG 11/2012
c
      subroutine tddft_grad_compute_p(ipol,naoc,nav,nroots,g_tp,g_z)
c
      implicit none
c
c     Note that the subroutine grad1 that evaluates the 1-electron
c     derivative integrals assumes that the density matrices are
c     symmetric. If they are not the routine returns rubbish. So 
c     instead of constructing
c
c              | Tij  Zia |
c        Ppq = |          |
c              |  0   Tab |
c
c     we have to construct
c
c              | Tij       Zia/2 |
c        Ppq = |                 |  where t(A) is the transpose of A.
c              | t(Zia/2)  Tab   |
c
c     This way we get a symmetric matrix that has the same effect
c     as the unsymmetric one apart from the fact that it will get the
c     right results from grad1.
c
#include "global.fh"
#include "stdio.fh"
#include "tddft_grad_util.fh"
c
c     Input:
c
      integer ipol    ! =1 (restricted), =2 (unrestricted)
      integer naoc(2) ! the number of active occupied orbitals
      integer nav(2)  ! the number of active virtual orbitals
      integer nroots  ! the number of roots
c
      integer g_z(2)  ! the global arrays holding Z
c
c     Input/Output:
c
      integer g_tp(2) ! on input: T, on output P
c
c     Local:
c
      integer ip      ! counter over spin components
      integer alo(3)  ! lower limits on A
      integer ahi(3)  ! upper limits on A
      integer blo(3)  ! lower limits on B
      integer bhi(3)  ! upper limits on B
      logical oroot
c
      character*32 pname
c
      pname="tddft_grad_compute_p: "
c
c     Copy Zia to Tia
c
      do ip = 1, ipol
        alo(1) = 1
        ahi(1) = nroots
        alo(2) = 1
        ahi(2) = naoc(ip)
        alo(3) = 1
        ahi(3) = nav(ip)
        blo(1) = 1
        bhi(1) = nroots
        blo(2) = 1
        bhi(2) = naoc(ip)
        blo(3) = naoc(ip)+1
        bhi(3) = naoc(ip)+nav(ip)
        call nga_copy_patch('n',g_z(ip),alo,ahi,g_tp(ip),blo,bhi)
c
c       scale the patch by half because of required symmetrization for
c       grad1
c
        call nga_scale_patch(g_tp(ip),blo,bhi,0.5d0)
      enddo
c
c     Copy Zia to Tai
c
      do ip = 1, ipol
        alo(1) = 1
        ahi(1) = nroots
        alo(2) = 1
        ahi(3) = nav(ip)
        alo(3) = 1
        ahi(2) = naoc(ip)
        blo(1) = 1
        bhi(1) = nroots
        blo(2) = naoc(ip)+1
        bhi(2) = naoc(ip)+nav(ip)
        blo(3) = 1
        bhi(3) = naoc(ip)
        call nga_copy_patch('t',g_z(ip),alo,ahi,g_tp(ip),blo,bhi)
c
c       scale the patch by half because of required symmetrization for
c       grad1
c
        call nga_scale_patch(g_tp(ip),blo,bhi,0.5d0)
      enddo
c
      if (util_print(pname,print_debug)) then
        oroot = ga_nodeid().eq.0
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'P'
        call tddft_grad_print_array(ipol,nroots,g_tp,dble(ipol))
      endif
c
      return
      end
c $Id$
