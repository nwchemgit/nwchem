c     $Id: et_2e.F,v 1.9 2002-03-27 22:16:51 d3k958 Exp $
      subroutine et_2e(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &                 g_pa,g_pb,omega2)
      implicit none
c
c     calculates the two-electron contribution to Vba, omega2
c
c     nbf, basis, geom                                          [input]
c     max2e, mem2 (for int_2e4c)                                [input]
c     g_pa, g_pb (handles for alpha and beta densities)         [input]
c     omega2 (two electron contrib to Vab)                      [output]
c
#include "rtdb.fh"
#include "schwarz.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "msgids.fh"   
c
      integer rtdb
      integer g_pa,g_pb
      integer nbf, geom, basis,max2e, mem2
      integer l_gJ, k_gJ, l_scr, k_scr
      integer ish, jsh, ksh, lsh, nsh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer i,j,nproc, num
      integer maxlsh, shellnum, ng
      double precision omega2, omega2_val
      integer nxtask, next
      external nxtask
      double precision start_cpu, end_cpu, elapsed
      double precision smax, sij, skl, Sba
      double precision pa(nbf,nbf), pb(nbf,nbf)
      double precision density(nbf*nbf)
      double precision norm, G(10000), dij
      double precision tol2e, maxPa, maxPb
      double precision density_shell, dabsmax, denmax
      external  density_shell
c     
c     get densities and keep in memory
c     ================================
c
c
      do i = 1,nbf
      do j = 1,nbf
        call ga_get(g_pa, i, i, j, j, pa(i,j),1)
        call ga_get(g_pb, i, i, j, j, pb(i,j),1)
      enddo
      enddo
c
      call ga_sync()
c
c     get the maximum value of the alpha and beta densities
c     ====================================================
c
      num = 1 
      do i = 1,nbf
      do j = 1,nbf
        density(num) = pa(i,j)
        num=num+1
      enddo
      enddo
c
      maxPa = dabsmax(nbf*nbf,density)
c
      num = 1 
      do i = 1,nbf
      do j = 1,nbf
        density(num) = pb(i,j)
        num=num+1
      enddo
      enddo
c
      maxPb = dabsmax(nbf*nbf,density)
c
      denmax = max(maxPa,maxPb)
      denmax = denmax*denmax
c
c     initialize omega2 and establish 2e and scr arrays
c     ================================================
c     
      omega2 = 0.0d0
      elapsed = 0.0d0
c     
      if (.not.ma_push_get(MT_DBL, max2e, '2e buffer', l_gJ, k_gJ))
     $     call errquit('et_2e:ma_push_get failed for 2e J buffer',k_gJ)
c     
      if (.not.ma_push_get(MT_DBL, mem2, '2e scr', l_scr, k_scr))
     $     call errquit('et_2e:ma_push_get failed for 2e scratch',k_scr)
c     
c     get nsh, the number of shells
c     =============================
c     
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('et_2e: bas_numcont failed', basis)
c     
c     read or set screening parameters
c     ================================
c
      call schwarz_init(geom, basis)
      smax = schwarz_max() 
c
c       if (.not. rtdb_get(rtdb,'et:Srp', MT_DBL, 1, Sba))
c     $        call errquit('et_calc: rtdb_put of Srp failed',0)
c      Sba = Sba*1.0d-7
cc
c      if (.not. rtdb_get(rtdb, 'et:tol2e', MT_DBL, 1, tol2e)) then
c       tol2e    =  min(1.0d-7,dabs(Sba))
c       tol2e    =  max(1.0d-12,tol2e)
c      endif
c
 1212 format (2x,a30,1x,a1,1pE14.4)
      if (ga_nodeid().eq.0) then
      write(6,*)
      write(6,*)'Start two-electron contribution,  H2(rp)...'
      write(6,1212)'Density screening tolerance   ',':',tol2e
      write(6,1212)'Largest value of alpha density',':',maxPa
      write(6,1212)'Largest value of beta density ',':',maxPb
      write(6,1212)'Largest two-electron integral ',':',smax
      write(6,*)
      endif
      call util_flush(6) 
c     
c     begin loop over shells
c     ======================
c     
      shellnum = 0
      nproc = ga_nnodes()
      next = nxtask(nproc,1)
c
      do ish = nsh, 1, -1
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('et_2e: bas_cn2bfr', ish)
c     
         do jsh = ish, 1, -1
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('et_2e: bas_cn2bfr', jsh)
c
             dij =  density_shell(ilo,ihi,jlo,jhi,ish,jsh,pa,pb,nbf) 
c
              sij = schwarz_shell(ish,jsh)
                if (sij*smax*dij .gt. tol2e) then
c     
               do ksh = ish, 1, -1
                  if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $                 call errquit('et_2e: bas_cn2bfr', ksh)
c     
                  maxlsh = ksh
                  if (ksh.EQ.ish) maxlsh = jsh
                  do lsh = maxlsh, 1, -1
                     if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $                 call errquit('et_2e: bas_cn2bfr', lsh)
c
                   if (shellnum.eq.next) then
c
                        skl = schwarz_shell(ksh,lsh)      
c
                        call dens_prod(nbf,
     $                    ilo,ihi,jlo,jhi,klo,khi,llo,lhi, 
     $                    ish,jsh,ksh,lsh, pa,pb, norm,G,ng) 
c
                        if (sij*skl*smax*norm.gt.tol2e) then
c
                           call int_2e4c(basis, ish, jsh, basis, 
     $                          ksh, lsh, mem2, dbl_mb(k_scr), max2e, 
     $                          dbl_mb(k_gJ))
c     
                           call omega2_mult(ilo,ihi,jlo,jhi,
     &                          klo,khi,llo,lhi,
     &                          ish,jsh,ksh,lsh,
     &                          G,ng,dbl_mb(k_gJ),omega2_val)
c
                           omega2 = omega2 + omega2_val
c     
                       endif  !norm screening
                   next = nxtask(nproc,1)
                   endif   !shellnum
                   shellnum = shellnum + 1
                  enddo
               enddo
            endif  !sij x smax
         enddo
      enddo
c     
      next = nxtask(-nproc,1)
      call ga_sync()
c     
      call ga_dgop(msg_et_2e,omega2,1,'+')
c     
c     destroy arrays
c     ==============
      call schwarz_tidy()
c     
      if (.not. ma_pop_stack(l_scr)) call errquit('et_2e: pop',0)
      if (.not. ma_pop_stack(l_gJ))  call errquit('et_2e: pop',1)
c
c
      return
      end



