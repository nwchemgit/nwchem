      subroutine tddft_grad_mat_to_vec(g_m, ilo, ihi, g_v, klo, col,
     $           scale, op)
c
c     $Id$
c
c     This subroutine is a copy of ga_mat_to_vec modified to handle
c     3D global arrays as they are used in the TDDFT Gradients. 
c     I realize that this introduces some redundancy but personally
c     I would like to see the need for this sort of routine to 
c     vanish altogether. The latter may be harder to achieve than I
c     think...
c
c     Huub van Dam, June 2006.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
c     v(klo:*,col:*)=scale*m(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3))
c
c     l = col
c     do m = ilo(1),ihi(1)
c     .  k = klo
c     .  do j = jlo,jhi
c     .     do i = ilo,ihi
c     .        if (op .eq. '=') then
c     .          v(k,l) = scale*m(i,j)
c     .        else if (op .eq. '+') then
c     .          v(k,l) = v(k,l) + scale*m(i,j)
c     .        else
c     .          error
c     .        k = k + 1
c     .     enddo
c     .  enddo
c     .  l = l + 1
c     enddo
c
      integer g_m, g_v
      integer ilo(3), ihi(3), klo, col
      character*1 op
      double precision scale
c
c     Local:
c
      integer vlo(2), vhi(2) 
c     integer iplo, iphi, jplo, jphi, k_m   ! local patch of m
c     integer ioff, j, numi, numj, k, ld
c
cifdef  BAD_GACCESS
      vlo(1) = klo
      vhi(1) = klo+(ihi(2)-ilo(2)+1)*(ihi(3)-ilo(3)+1)-1
      vlo(2) = col
      vhi(2) = col+(ihi(1)-ilo(1)+1)-1
      if (op .eq. '=') then
         call nga_copy_patch('n', g_m, ilo, ihi, g_v, vlo, vhi) 
         call nga_scale_patch(g_v, vlo, vhi, scale)
      else if (op .eq. '+') then
         call nga_add_patch(scale, g_m, ilo, ihi, 1.0d0, g_v, vlo, vhi,
     +                      g_v, vlo, vhi)
      else
         call errquit('tddft_grad_mat_to_vec: bad op',0, 0)
      endif
celse
c     call ga_distribution(g_m, ga_nodeid(), iplo, iphi, jplo, jphi)
c
c     call ga_sync()
c     if (iplo.gt.0 .and. iplo.le.iphi) then
c        numi =  min(ihi,iphi) - max(ilo,iplo) + 1
c        numj =  min(jhi,jphi) - max(jlo,jplo) + 1
c        if (numi.gt.0 .and. numj.gt.0) then
c          call ga_access(g_m,iplo,iphi,jplo,jphi,k_m,ld)
c          do j = max(jlo,jplo), min(jhi,jphi)
c             k = klo + (j-jlo)*(ihi-ilo+1) + max(ilo,iplo) - ilo
c             ioff = k_m + (j-jplo)*ld + max(ilo,iplo) - iplo
c             if (op .eq. '=') then
c               call ga_put(g_v, k, k+numi-1, col, col, dbl_mb(ioff), 1)
c             else if (op .eq. '+') then
c               call ga_acc(g_v, k, k+numi-1, col, col, dbl_mb(ioff), 1,
c    $                      scale)
c             else
c               call errquit('ga_mat_to_vec: bad op', 0, GA_ERR)
c             endif
c          enddo
c          call ga_release(g_m,iplo,iphi,jplo,jphi)
c        endif
c     endif
cendif
c     call ga_sync()
c
      end
      subroutine tddft_grad_vec_to_mat(g_m, ilo, ihi, g_v, klo, col)
c
c     $Id$
c
c     This subroutine is a copy of ga_vec_to_mat modified to handle
c     3D global arrays as they are used in the TDDFT Gradients. 
c     I realize that this introduces some redundancy but personally
c     I would like to see the need for this sort of routine to 
c     vanish altogether. The latter may be harder to achieve than I
c     think...
c
c     Huub van Dam, June 2006.
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
c     m(ilo(1):ihi(1),ilo(2):ihi(2),ilo(3):ihi(3)) = v(klo:*,col:*)
c
c     l = col
c     do m = ilo(1),ihi(1)
c     .  k = klo
c     .  do j = jlo,jhi
c     .     do i = ilo,ihi
c     .        m(i,j) = v(k,l)
c     .        k = k + 1
c     .     enddo
c     .  enddo
c     .  l = l + 1
c     enddo
c
      integer g_m, g_v
      integer ilo(3), ihi(3), klo, col
c
      integer vlo(2), vhi(2) 
c     integer iplo, iphi, jplo, jphi, k_m   ! local patch of m
c     integer ioff, j, numi, numj, k, ld
c
cifdef  BAD_GACCESS
      vlo(1) = klo
      vhi(1) = klo+(ihi(2)-ilo(2)+1)*(ihi(3)-ilo(3)+1)-1
      vlo(2) = col
      vhi(2) = col+(ihi(1)-ilo(1)+1)-1
      call nga_copy_patch('n', 
     +        g_v, vlo, vhi, 
     +        g_m, ilo, ihi) 
celse 
c     call ga_distribution(g_m, ga_nodeid(), iplo, iphi, jplo, jphi)
c
c     call ga_sync()
c     if (iplo.gt.0 .and. iplo.le.iphi) then
c        numi =  min(ihi,iphi) - max(ilo,iplo) + 1
c        numj =  min(jhi,jphi) - max(jlo,jplo) + 1
c        if (numi.gt.0 .and. numj.gt.0) then
c          call ga_access(g_m,iplo,iphi,jplo,jphi,k_m,ld)
c          do j = max(jlo,jplo), min(jhi,jphi)
c             k = klo + (j-jlo)*(ihi-ilo+1) + max(ilo,iplo) - ilo
c             ioff = k_m + (j-jplo)*ld + max(ilo,iplo) - iplo
c             call ga_get(g_v, k, k+numi-1, col, col, dbl_mb(ioff), 1)
c          enddo
c          call ga_release(g_m,iplo,iphi,jplo,jphi)
c        endif
c     endif
cendif
c     call ga_sync()
c
      end
