      Subroutine dftg_getxc(rtdb, natom, iga_dens, force, oactive, 
     &                      nactive, s34)
*
* $Id: dftg_getxc.F,v 1.10 1997-11-04 10:12:14 d3e129 Exp $
*
      implicit none
c      
c     compute XC component of energy gradient.  
c
      integer rtdb
      integer natom
      integer iga_dens(2)
      double precision force(3*natom), s34(natom)
      integer nactive
      logical oactive(natom)
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "sym.fh"
#include "cdft.fh"
c     
      double precision zero, one, onem
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      logical status
c
      integer nproc, max_at_bf, lbuf, max1e, max2e,
     &        lscratch, mscratch_1e, mscratch_2e, i, lsqa
      integer l_buf, k_buf, l_scr, k_scr, l_ik, i_ik, l_il, i_il,
     &        l_jk, i_jk, l_jl, i_jl, ilo, ihi
      double precision tol2e, tot
c
      nproc = ga_nnodes()
c
      if (abs(xfac(1)).gt.1e-8)then
c       
c        Compute the exact exchange potential (as in Hartree-Fock 
c        calculations).
c       
         tol2e=10.d0**(-itol2e)
c
         max_at_bf = 0
         do i = 1, natom
            if (.not. bas_ce2bfr(ao_bas_han, i, ilo, ihi))
     &         call errquit('dftg_getxc: bas_ce2bfr failed', i)
            max_at_bf = max(max_at_bf, ihi-ilo+1)
         enddo
c
         call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
         lbuf = max(max1e, max2e)
         lscratch = max(mscratch_1e, mscratch_2e)
         status = MA_Push_get(MT_DBL, lbuf, 'deriv buffer', l_buf,
     &            k_buf)
         if (.not.status)call errquit('could not allocate buffer',1)
         status = MA_Push_get(MT_DBL, lscratch, 'deriv scratch', 
     &            l_scr, k_scr)
         if (.not.status)call errquit('could not allocate deriv buf',1)
c
c        allocate local density matrix block
c
         lsqa = max_at_bf * max_at_bf
c
         status = status.and.MA_Push_get(MT_DBL,lsqa,
     &                       'local_w_density', l_ik, i_ik)
         status = status.and.MA_Push_get(MT_DBL,lsqa,
     &                       'local_w_density', l_il, i_il)
         status = status.and.MA_Push_get(MT_DBL,lsqa,
     &                       'local_w_density', l_jk, i_jk)
         status = status.and.MA_Push_get(MT_DBL,lsqa,
     &                       'local_w_density', l_jl, i_jl)
         if (.not.status) 
     &      call errquit('force_dft: could not allocate scratch',1)
c
c        store total DM in ga_dens(1)
c
         if (ipol .eq. 2)then
            call ga_dadd (one,iga_dens(1),one,iga_dens(2),iga_dens(1))
         endif 
c
         call dftg_hf_x(xfac(1), dbl_mb(i_ik), dbl_mb(i_jl),
     &                  dbl_mb(i_il), dbl_mb(i_jk), force, iga_dens,
     &                  max_at_bf, AO_bas_han, nproc, natom, 
     &                  dbl_mb(k_scr), lscratch, dbl_mb(k_buf), lbuf, 
     &                  tol2e, oactive, nactive, geom, s34, oskel)
c
         status = MA_Pop_Stack (l_jl)
         status = MA_Pop_Stack (l_jk)
         status = MA_Pop_Stack (l_il)
         status = MA_Pop_Stack (l_ik)
         status = MA_Pop_Stack (l_scr)
         status = MA_Pop_Stack (l_buf)
         if (.not.status) 
     &      call errquit('force_dft: could not pop stack',1)
c
c        restore alpha DM in g_dens(1)
c
         if (ipol .eq. 2)then
            call ga_dadd (one,iga_dens(1),onem,iga_dens(2),iga_dens(1))
         endif 
c 
         if (.not.status) 
     &      call errquit(' dftg_getcx: failed to deallocate',11)
c
      endif
      call ga_sync
      tot = -xfac(1)
      do i = 1, 15
         tot = tot + xfac(i) + cfac(i)
      enddo
c
      if (abs(tot).gt.1e-8)then
c     
c        Compute via direct numerical quadrature.
c     
         call dftg_quadv0(rtdb, iga_dens, force, natom, oactive, 
     &                    nactive) 
      endif
c
      return
      end
