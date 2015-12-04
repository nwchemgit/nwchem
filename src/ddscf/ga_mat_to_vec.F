      subroutine ga_mat_to_vec(g_m, ilo, ihi, jlo, jhi, g_v, klo, col,
     $           scale, op)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
c     v(klo:*,col) = m(ilo:ihi,jlo:jhi)
c
c     k = klo
c     do j = jlo,jhi
c     .  do i = ilo,ihi
c     .     if (op .eq. '=') then
c     .       v(k,col) = m(i,j)
c     .     else if (op .eq. '+') then
c     .       v(k,col) = v(k,col) + scale*m(i,j)
c     .     else
c     .       error
c     .     k = k + 1
c     .  enddo
c     enddo
c
      integer g_m, g_v
      integer ilo, ihi, jlo, jhi, klo, col
      character*1 op
      double precision scale
c
      integer iplo, iphi, jplo, jphi, k_m   ! local patch of m
      integer ioff, j, numi, numj, k, ld
c
#ifdef  BAD_GACCESS
      if (op .eq. '=') then
         call ga_copy_patch('n', 
     $        g_m, ilo, ihi, jlo, jhi, 
     $        g_v, klo, klo+(ihi-ilo+1)*(jhi-jlo+1)-1, col, col) 
         
      else if (op .eq. '+') then
         call ga_dadd_patch(scale, 
     $        g_m, ilo, ihi, jlo, jhi, 
     $        1d0,g_v, klo, klo+(ihi-ilo+1)*(jhi-jlo+1)-1, col, col,
     $        g_v, klo, klo+(ihi-ilo+1)*(jhi-jlo+1)-1, col, col) 
      else
         call errquit('ga_mat_to_vec: bad op',0, 0)
      endif
#else
      call ga_distribution(g_m, ga_nodeid(), iplo, iphi, jplo, jphi)
c
      call ga_sync()
      if (iplo.gt.0 .and. iplo.le.iphi) then
         numi =  min(ihi,iphi) - max(ilo,iplo) + 1
         numj =  min(jhi,jphi) - max(jlo,jplo) + 1
         if (numi.gt.0 .and. numj.gt.0) then
           call ga_access(g_m,iplo,iphi,jplo,jphi,k_m,ld)
           do j = max(jlo,jplo), min(jhi,jphi)
              k = klo + (j-jlo)*(ihi-ilo+1) + max(ilo,iplo) - ilo
              ioff = k_m + (j-jplo)*ld + max(ilo,iplo) - iplo
              if (op .eq. '=') then
                call ga_put(g_v, k, k+numi-1, col, col, dbl_mb(ioff), 
     N                numi)
              else if (op .eq. '+') then
                call ga_acc(g_v, k, k+numi-1, col, col, dbl_mb(ioff), 
     N                numi, scale)
              else
                call errquit('ga_mat_to_vec: bad op', 0, GA_ERR)
              endif
           enddo
           call ga_release(g_m,iplo,iphi,jplo,jphi)
         endif
      endif
#endif
      call ga_sync()
c
      end
      subroutine ga_vec_to_mat(g_m, ilo, ihi, jlo, jhi, g_v, klo, col)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c
c     m(ilo:ihi,jlo:jhi) = v(klo:*,col) 
c
c     k = klo
c     do j = jlo,jhi
c     .  do i = ilo,ihi
c     .     m(i,j) = v(k,col)
c     .     k = k + 1
c     .  enddo
c     enddo
c
      integer g_m, g_v
      integer ilo, ihi, jlo, jhi, klo, col
c
      integer iplo, iphi, jplo, jphi, k_m   ! local patch of m
      integer ioff, j, numi, numj, k, ld
c
#ifdef  BAD_GACCESS
      call ga_copy_patch('n', 
     $     g_v, klo, klo+(ihi-ilo+1)*(jhi-jlo+1)-1, col, col,
     $     g_m, ilo, ihi, jlo, jhi)
#else 
      call ga_distribution(g_m, ga_nodeid(), iplo, iphi, jplo, jphi)
c
      call ga_sync()
      if (iplo.gt.0 .and. iplo.le.iphi) then
         numi =  min(ihi,iphi) - max(ilo,iplo) + 1
         numj =  min(jhi,jphi) - max(jlo,jplo) + 1
         if (numi.gt.0 .and. numj.gt.0) then
           call ga_access(g_m,iplo,iphi,jplo,jphi,k_m,ld)
           do j = max(jlo,jplo), min(jhi,jphi)
              k = klo + (j-jlo)*(ihi-ilo+1) + max(ilo,iplo) - ilo
              ioff = k_m + (j-jplo)*ld + max(ilo,iplo) - iplo
              call ga_get(g_v, k, k+numi-1, col, col, dbl_mb(ioff),numi)
           enddo
           call ga_release(g_m,iplo,iphi,jplo,jphi)
         endif
      endif
#endif
      call ga_sync()
c
      end
