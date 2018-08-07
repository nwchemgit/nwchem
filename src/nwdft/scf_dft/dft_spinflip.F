      subroutine dft_spinflip(g_dens,g_s,basis,c_a,c_b,
     S     no_sflip,at_flip)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
      integer g_dens(2)
      integer g_s
      integer basis,geom
      double precision c_a,c_b
      integer no_sflip
      integer at_flip(no_sflip)
c
      integer at,ii,nbf
      integer ilo,ihi,ld
      integer l_dm,k_dm
      double precision pstrace_in,pstrace_out
c
c g_dens(2) has -B
c     compute spin DM and store it in g_dens(2)
c
!      call ga_add(1d0, g_dens(1), -1d0, g_dens(2), g_dens(2))
clast      call ga_add(1d0, g_dens(1), -2d0, g_dens(2), g_dens(2))
      call ga_add(c_a, g_dens(1), c_b, g_dens(2), g_dens(2))
c
c     check on trace(spin_dm *S)
c
      pstrace_in=ga_ddot(g_dens(2),g_s)
      if(ga_nodeid().eq.0) then
         write(6,*) ' spinflip: input pstrace ',pstrace_in
         if (.not. bas_numbf(basis, nbf)) call errquit
     $        ('dft_spinflip: bas_numbf ?', 0, BASIS_ERR)
         
         do ii=1,no_sflip
           at=at_flip(ii)
            write(6,*) '  spinflip: flipping at ',at
c     
c     grab atomic block and invert sign
c
            if (.not. bas_ce2bfr(basis, at, ilo, ihi)) call errquit
     $           ('dft_spinflip: bas_ce2bfr ?', 0, BASIS_ERR)
            ld=ihi-ilo+1
            if (ld .gt. 0) then
               if (.not.MA_Push_Get(MT_Dbl,ld*ld,'dm',
     F              l_dm,k_dm))
     &              call errquit('dft_spinf: cannot allocate',0, MA_ERR)
               call ga_get(g_dens(2),ilo,ihi,ilo,ihi,dbl_mb(k_dm),ld)
               call dscal(ld*ld,-1d0,dbl_mb(k_dm),1)
               call ga_put(g_dens(2),ilo,ihi,ilo,ihi,dbl_mb(k_dm),ld)
               if (.not.ma_pop_stack(l_dm))
     &              call errquit('dft_spinf: cannot pop stack',3,MA_ERR)
            else
               call errquit('dftspinflip: no basis for atom ',
     B              at,BASIS_ERR)
            endif
         enddo
      endif
      call ga_sync()
      pstrace_out=ga_ddot(g_dens(2),g_s)
      if(ga_nodeid().eq.0)  then
         write(6,*) ' spinflip: output pstrace ',pstrace_out
         if(abs(pstrace_in-pstrace_out).gt.1d-3) write(6,*)
     W        'WARNING: large change in pstrace!'
      endif
      
c     restore - beta in g_dens(2)
c       A - (A -B) = B
!      call ga_add(1d0, g_dens(1), -1d0, g_dens(2), g_dens(2))
clast      call ga_add(0.5d0, g_dens(1), -0.5d0, g_dens(2), g_dens(2))
      call ga_add(-c_a/c_b, g_dens(1), 1d0/c_b, g_dens(2), g_dens(2))
      return
      end
