      double precision function potential(basis, g_dens, x, y, z)
      implicit none
#include "bas.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      integer basis, g_dens
      double precision x, y, z
c
c     Compute the electrostatic potential due to the density
c     at the given point.
c
      integer nshell
      integer ishell, jshell, mem1, max1e
      integer ilo, ihi, jlo, jhi, idim, jdim
      integer l_buf, l_scr
      integer k_buf, k_scr
      double precision pot, r(3)
c
      r(1) = x
      r(2) = y
      r(3) = z
c
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('potential: bas_numcont failed for basis', basis)
c     
c     l_buf ... buffer to hold shell block of matrix
c     l_scr ... workspace for integral routines
c     
      call int_mem_1e(max1e, mem1)
      if (.not. MA_push_get(MT_DBL,max1e,'potential:buf',l_buf,k_buf))
     $     call errquit('potential: ma failed', max1e)
      if (.not. MA_push_get(MT_DBL, mem1,'potential:scr',l_scr,k_scr))
     $     call errquit('potential: ma failed', mem1)
c     
**      call ga_print(g_dens)
      pot = 0.0d0
      do jshell = 1, nshell
         if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     $        call errquit('potential: bas_cn2bfr ?', basis)
         jdim = jhi - jlo + 1
         do ishell = 1, nshell
            if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) 
     $           call errquit('potential: bas_cn2bfr ?', basis)
c
            idim = ihi - ilo + 1
            call int_1er (r,1,1.0d0, basis, jshell, basis, ishell, 
     $           mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
*            call int_1eov (basis, jshell, basis, ishell, 
*     $           mem1, dbl_mb(k_scr), max1e, dbl_mb(k_buf))
            call ga_get(g_dens,ilo,ihi,jlo,jhi,dbl_mb(k_scr),idim)
c
            pot = pot + ddot(idim*jdim,dbl_mb(k_scr),1,dbl_mb(k_buf),1)
         enddo
      enddo
c
c$$$      call ga_zero(g_dens)
c$$$      call int_1e_ga(basis, basis, g_dens, 'overlap', .false.)
c$$$      call ga_print(g_dens)
c$$$      stop
      if (.not. MA_pop_stack(l_scr)) call errquit
     $     ('potential: pop failed', 0)
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('potential: pop failed', 0)
c
      potential = pot
c
      end
