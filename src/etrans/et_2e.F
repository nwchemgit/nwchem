c     $Id: et_2e.F,v 1.4 2002-02-19 21:38:23 d3k958 Exp $
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
      double precision smax, sij, skl, tol
      double precision pa(nbf,nbf), pb(nbf,nbf)
      double precision pt2(nbf,nbf)
      double precision norm, G(10000), dij, dkl
      double precision density_shell
      external  density_shell
c     

c     get densities and keep in memory
c     ================================
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
c      write(6,*)'nsh',nsh
c     
c     begin loop over shells
c     ======================
c     
      smax = schwarz_max() 
      shellnum = 0
      nproc = ga_nnodes()
      next = nxtask(nproc,100)
c     
      tol = 1.0d-9
c     
      do ish = 1, nsh
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('et_2e: bas_cn2bfr', ish)
c     
         do jsh = 1, ish
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('et_2e: bas_cn2bfr', jsh)
c
              sij = schwarz_shell(ish,jsh)
              dij = density_shell(ilo,ihi,jlo,jhi,pa,pb,nbf)
                if (sij*smax .gt. tol) then
c     
               do ksh = 1, ish
                  if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $                 call errquit('et_2e: bas_cn2bfr', ksh)
c     
                  maxlsh = ksh
                  if (ksh.EQ.ish) maxlsh = jsh
                  do lsh = 1, maxlsh
                     if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $                 call errquit('et_2e: bas_cn2bfr', lsh)
c
                     call dens_prod(nbf,
     $                  ilo,ihi,jlo,jhi,klo,khi,llo,lhi, 
     $                  ish,jsh,ksh,lsh, pa,pb, norm,G,ng) 
c
                        skl = schwarz_shell(ksh,lsh)      
                        dkl = density_shell(klo,khi,llo,lhi,pa,pb,nbf)
c
                        if (norm*sij*skl.gt. tol) then
c
                        if (shellnum.eq.next) then
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
c                           write(6,1113) ish,jsh,ksh,lsh,dij,omega2_val
                           omega2 = omega2 + omega2_val
c     
                           next = nxtask(nproc,100)
                        endif
                        shellnum = shellnum + 1
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
c     
      next = nxtask(-nproc,1)
      call ga_sync()
c      write(6,*)'et_2e time: ',ga_nodeid(),elapsed
c     
      call ga_dgop(msg_et_2e,omega2,1,'+')
c     
c     destroy arrays
c     ==============
c     
      if (.not. ma_pop_stack(l_scr)) call errquit('et_2e: pop',0)
      if (.not. ma_pop_stack(l_gJ))  call errquit('et_2e: pop',1)
c     
 1111 format(a18,1x,i4,1x,i4,1x,i4,1x,i4,3(1x,F24.12))
 1113 format(4i3,1x,2(F20.16),a10)
 1112 format (a31,1x,4F16.12)
c     
      return
      end



