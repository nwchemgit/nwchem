      Subroutine dftg_getxc(rtdb, geom, AO_bas_han, ipol, natom, 
     &                      iga_dens, force, oactive, nactive,
     &                      s34, oskel)
c
      implicit none
c      
c     compute XC component of energy gradient.  
c
      integer rtdb, geom, AO_bas_han
      integer ipol
      integer natom
      integer iga_dens(2)
      double precision force(3*natom), s34(natom)
      integer nactive
      logical oactive(natom), oskel
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "sym.fh"
c     
      double precision zero, one, onem
      parameter(zero = 0.0d0, one = 1.0d0, onem = -1.0d0)
      logical status
      double precision xfac(15), cfac(15)
c
      integer nproc, nbf_ao, itol2e, max_at_bf, lbuf, max1e, max2e,
     &        lscratch, mscratch_1e, mscratch_2e, i, lsqa
      integer l_buf, k_buf, l_scr, k_scr, l_ik, i_ik, l_il, i_il,
     &        l_jk, i_jk, l_jl, i_jl
      double precision tol2e, tot
c
      nproc = ga_nnodes()
c
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, 15, cfac))
     &   call errquit('xc_getv: rtdb_put failed', 0)
      if (.not. rtdb_get(rtdb, 'dft:xfac', mt_dbl, 15, xfac))
     &   call errquit('xc_getv:rtdb_put failed', 0)
      if (.not. bas_numbf(AO_bas_han,nbf_ao))then
         call errquit('Exiting in getdexc.',1)
      endif
c
      if (abs(xfac(1)).gt.1e-8)then
c       
c        Compute the exact exchange potential (as in Hartree-Fock 
c        calculations).
c       
         if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1, itol2e))
     &      call errquit('fitcd: rtdb_get failed', 0)
         tol2e=10.d0**(-itol2e)
c
         status = bas_nbf_ce_max(AO_bas_han, max_at_bf)
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
         call dftg_quadv0(rtdb, geom, AO_bas_han, ipol, iga_dens, 
     &                    force, natom, oactive, nactive, oskel) 
      endif
c
      return
      end
