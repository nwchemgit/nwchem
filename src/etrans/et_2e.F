c     $Id: et_2e.F,v 1.7 2002-03-06 15:59:59 d3k958 Exp $
      subroutine et_2e(nbf,basis,geom,max2e,mem2,g_pa,g_pb,omega2)
      implicit none
c
c     calculates the two-electron contribution to Vba, omega2
c
c     nbf, basis, geom                                          [input]
c     max2e, mem2 (for int_2e4c)                                [input]
c     g_pa, g_pb (handles for alpha and beta densities)         [input]
c     omega2 (two electron contrib to Vab)                      [output]
c
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
      integer g_pa,g_pb
      integer g_paa,g_pbb
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
      double precision smax, sij, skl
      double precision pa(nbf,nbf), pb(nbf,nbf)
      double precision paa(nbf,nbf), pbb(nbf,nbf)
      double precision pt2(nbf,nbf)
      double precision norm, G(10000), dij, dkl
      double precision tol2e, tol_dij, tol_dkl
      double precision ij_density_shell, kl_density_shell
      external  ij_density_shell, kl_density_shell
c     
c     calculate alpha-alpha and beta-beta density products
c     ====================================================
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'paa',0,0,g_paa))
     $        call errquit('et_2e: ga_create failed', 0)    
c
      if(.not.ga_create(mt_dbl,nbf,nbf,'pbb',0,0,g_pbb))
     $        call errquit('et_2e: ga_create failed', 0)
c
      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_pa,g_pa,0.0d0,g_paa)
      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_pb,g_pb,0.0d0,g_pbb)
c
c     get densities and keep in memory
c     ================================
c
c      call ga_get(g_pa, 1, nbf, 1, nbf, pa, 1)
c      call ga_get(g_pb, 1, nbf, 1, nbf, pb, 1)
c
      do i = 1,nbf
      do j = 1,nbf
        call ga_get(g_pa, i, i, j, j, pa(i,j),1)
        call ga_get(g_pb, i, i, j, j, pb(i,j),1)
        call ga_get(g_paa, i, i, j, j, paa(i,j),1)
        call ga_get(g_pbb, i, i, j, j, pbb(i,j),1)
      enddo
      enddo
c
      call ga_sync()
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
c     begin loop over shells
c     ======================
c     
      call schwarz_init(geom, basis)
      smax = schwarz_max() 
      shellnum = 0
      nproc = ga_nnodes()
      next = nxtask(nproc,1)
c     
      tol2e    = 1.0d-8
      tol_dij  = 1.0d-9
      tol_dkl  = 1.0d-9
c
      if (ga_nodeid().eq.0) then
      write(6,*)
      write(6,*)'Start two-electron contribution,  H2(rp)...'
      write(6,*)'Density screening tolerances:'
      write(6,*)'tol2e   = ', tol2e
      write(6,*)'tol_dij = ', tol_dij
      write(6,*)'tol_dkl = ', tol_dkl
      endif
      call util_flush(6) 
c     
      do ish = nsh, 1, -1
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('et_2e: bas_cn2bfr', ish)
c     
         do jsh = ish, 1, -1
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('et_2e: bas_cn2bfr', jsh)
c
             dij =  ij_density_shell(ilo,ihi,jlo,jhi,pa,pb,paa,pbb,nbf) 
             if (dij.gt.tol_dij) then
c
              sij = schwarz_shell(ish,jsh)
                if (sij*smax .gt. tol2e) then
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
             dkl =  kl_density_shell(klo,khi,llo,lhi,pa,pb,nbf) 
             if (dkl.gt.tol_dkl) then
c
                   if (shellnum.eq.next) then
c
                       call dens_prod(nbf,
     $                    ilo,ihi,jlo,jhi,klo,khi,llo,lhi, 
     $                    ish,jsh,ksh,lsh, pa,pb, norm,G,ng) 
c
                        skl = schwarz_shell(ksh,lsh)      
c
                        if (sij*skl.gt.tol2e.AND.norm.gt.tol2e) then
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
c                write(6,1111) dij,dkl,norm,dabs(omega2_val)
c     
                       endif  !norm x schwarz screening
                   next = nxtask(nproc,1)
                   endif   !shellnum
                   shellnum = shellnum + 1
                   endif   !dkl
                  enddo
               enddo
            endif  !sij x smax
            endif  !dij
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
c     
      if (.not.ga_destroy(g_paa))
     &  call errquit('could not destroy paa array',g_paa)
      if (.not.ga_destroy(g_pbb))
     &  call errquit('could not destroy pbb array',g_pbb)
c
      if (.not. ma_pop_stack(l_scr)) call errquit('et_2e: pop',0)
      if (.not. ma_pop_stack(l_gJ))  call errquit('et_2e: pop',1)
c
 1111 format(3F20.16,3x,F20.16)     
      return
      end



