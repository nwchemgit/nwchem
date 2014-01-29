c     
c
c  $Id$
c
c
c  This routine returns the Coulomb and exchange integral
c  operator matrices for the range of MO-indices as mo_indx_hi, mo_indx_lo
c  The g_coul, g_exch global arrays are ordered as
c
c
c
c               ij
c  (ij|ab) = ( J  )  = g_coul[ ij : (a-1)*N2 + b ] = g_coul [ ij : (b-1)*N2 + a ]
c                  ab
c
c               ij
c  (ia|jb) = ( K  )  = g_exch[ ij : (a-1)*N2 + b ]
c                  ab
c
c -----------------
c Note: This routine differs from the standard in that
c it computes the integrals 6 times as opposed to 2 times
c in the standard. However, it has lower communication scaling
c for K, O(N^4), compared to O(N^5) in the standard. Therefore
c this routine is more appropriate for very high parallelism
c -----------------
c
c
c
       subroutine moints_build_6x(basis, ousesym,
     $                            occ_start, mo1_lo, mo1_hi,
     $                            mo2_lo, mo2_hi,
     $                            g_movecs,
     $                            g_coul, ocoul,
     $                            g_exch, oexch,
     $                            blen, oblk )
       implicit none
#include "errquit.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "util.fh"
#include "sym.fh"
#include "schwarz.fh"
#include "msgids.fh"
c
c Arguments
c
       integer basis                          ! Basis handle
       logical ousesym                        ! Symmetry toggle
       integer occ_start                      ! Offset from frozen core
       integer mo1_lo, mo1_hi                 ! 1st Pair Index range
       integer mo2_lo, mo2_hi                 ! 2nd Pair Index range
       integer g_movecs                       ! MO coefficients
       integer g_coul                         ! Coulomb operator
       integer g_exch                         ! Exchange operator
       logical ocoul,oexch                    ! Type selection
       integer blen                           ! Blocking length
       logical oblk                           ! Toggle AO blocking
c
c Local variables
c
       integer geom, nmo, nmo1, nmo2, nbf, nsh, maxbfsh
       integer qlo, qhi
       integer nblk, bsize, ngrp
       integer ii, jj, iz, kz, jz
       integer l_gmap, k_gmap
       integer ish, jsh, ilen, jlen
       integer ibflo,ibfhi,jbflo,jbfhi,kbflo,kbfhi,lbflo,lbfhi
       integer kshlo, kshhi, lshlo, lshhi
       integer kblen, lblen
       integer kgr, lgr
       integer i, nmixed
       integer l_ssbb, k_ssbb, l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_mo, k_mo, l_xmo, k_xmo
       integer l_shmap, k_shmap, l_bfmap, k_bfmap, l_rbfmap, k_rbfmap
       integer l_glo, k_glo, l_ghi, k_ghi
       integer l_sym, k_sym
       integer n_ssbb, n_ssbb1, n_ssni, n_hlp, n_ijni
       integer mem2, max2e
       integer num_nodes, ploop, next
       double precision tol2e, scale, scale0, schw_ij, integ_acc
       double precision tz, thalf, tint, t1jidx, t2jidx, t1kidx, t2kidx
       double precision t34kidx, t34jidx, thalf0, ttotal, tsynch
       double precision ttask, ttaskmax, ttaskmin, ttaskagg
       double precision flop1, flop2k, flop2j, q2
       integer tottask
       logical status, st, osym, odoit
#include "moints_stats.fh"
#include "moints.fh"
c
c
c
c      integer moints_numgr, gr_len, moints_nxttask
c      external moints_numgr, gr_len, moints_nxttask
       integer  gr_len
       external gr_len
c
c
c
       data tol2e/1.d-12/
c
c
c
       integ_acc = min(1d-9, max(0.01d0*tol2e,1d-17))
       call int_acc_set(integ_acc)
c
c  General basis info
c
       ttotal = util_cpusec()
       num_nodes = ga_nnodes()
       if (.not. bas_geom(basis, geom)) call errquit
     $      ('moints: cannot get geometry', 0, GEOM_ERR)
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not.status) call errquit('moints: cannot get basis info',0,
     &       BASIS_ERR)
       qlo  = min(occ_start,mo2_lo,mo1_lo)
       qhi  = max(mo1_hi,mo2_hi)
       nmo  = qhi - qlo + 1
       nmo1 = mo1_hi - mo1_lo + 1
       nmo2 = mo2_hi - mo2_lo + 1
c     
c     Symmetry adapt the MOs and renumber irreps to start at zero
c     
       osym = (ousesym.and.(sym_number_ops(geom).gt.0))
       if (osym) then
          if (.not. ma_push_get(MT_INT, nbf, 'movec syms',l_sym, k_sym))
     $         call errquit('moints_6x: no memory for syms?',0,
     &       MA_ERR)
          call sym_movecs_adapt(basis, 1d-8, g_movecs, int_mb(k_sym),
     $                          nmixed)
          if (nmixed .ne. 0)
     $         call errquit('moints_6x: symmetry contamination', nmixed,
     &       INPUT_ERR)
          do i =0, nbf-1
             int_mb(k_sym+i) = int_mb(k_sym+i) - 1
          enddo
          if (util_print('orbital symmetry',print_debug)) then
            write(6,887) 
 887        format('Symmetry of MOs')
            write(6,888) (int_mb(k_sym+i),i=0,nbf-1)
 888        format(16i3)
          endif
       endif
c
c  Initialize AO disk caching
c       
       status = status .and. moints_aodisk_prep()
c
c  Local MO coefficients
c
       status = status.and.ma_push_get(MT_DBL,(nbf*nmo),
     $                                 'movecs cols',l_mo,k_mo)
       call ga_get(g_movecs, 1, nbf, qlo, qhi, dbl_mb(k_mo), nbf)
c
c  Integrals allocation
c
       if (oblk) then
         call intb_mem_2e4c( max2e, mem2 )
         max2e = max( max2e, maxbfsh*maxbfsh*blen*blen )
       else
         call int_mem_2e4c( max2e, mem2)
       endif
       status = ma_push_get(MT_DBL, max2e,'moints: buf', l_eri, k_eri)
       status = ma_push_get(MT_DBL, mem2, 'moints: scr', l_iscr, k_iscr)
c
c  Shell group mapping
c
       nblk = moints_numgr( basis, blen )
       status = ma_push_get(MT_INT,(nblk*4),'shell group map',
     $                      l_gmap, k_gmap)
       call moints_grmap( basis, blen, nblk, int_mb(k_gmap))
c
c  Reorder and group shells by type
c
       status = ma_push_get(MT_INT,nsh,'shell order map',
     $                      l_shmap, k_shmap)
       status = ma_push_get(MT_INT,nsh,'group lo', l_glo, k_glo )
       status = ma_push_get(MT_INT,nsh,'group hi', l_ghi, k_ghi)
       status = ma_push_get(MT_INT,nbf,'basis map',
     $                      l_bfmap, k_bfmap)
       status = ma_push_get(MT_INT,nbf,'rev basis map',
     $                      l_rbfmap, k_rbfmap)
       call moints_shorder( basis, nsh, nbf, blen, ngrp,
     $                      int_mb(k_glo), int_mb(k_ghi),
     $                      int_mb(k_shmap),
     $                      int_mb(k_bfmap), int_mb(k_rbfmap) )
c
c  Copy of MO coefficients with reordered rows
c
       status = ma_push_get(MT_DBL,(nbf*nmo),'reorder mos',
     $                      l_xmo, k_xmo)
       call row_exch( nbf, nmo, int_mb(k_rbfmap), dbl_mb(k_mo),
     $                       dbl_mb(k_xmo) )
c
c  Temporary partially-transformed arrays
c     
       bsize = max(blen,maxbfsh)
       n_ssbb = maxbfsh*maxbfsh*bsize*bsize
       n_ssbb1 = max((nmo1*nmo1),n_ssbb)
       n_hlp = max((bsize*maxbfsh*maxbfsh*nmo1),(maxbfsh*nbf))
       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
       status = ma_push_get(MT_DBL,n_ssbb1,'ssbb block',l_ssbb,k_ssbb)
       status = ma_push_get(MT_DBL,n_ssbb,'ssbbt block',l_ssbbt,k_ssbbt)
       status = ma_push_get(MT_DBL,n_hlp,'hlp block',l_hlp,k_hlp)
       status = ma_push_get(MT_DBL,n_ssni,'ssni block',l_ssni,k_ssni)
       if (.not.(status)) call errquit
     $      ('moints_6x: cannot allocate local memory',0, MA_ERR)
c
c Initialize
c
       flop1 = 0.d0
       flop2k = 0.d0
       flop2j = 0.d0
       tint = 0.d0
       t1jidx = 0.d0
       t2jidx = 0.d0
       t1kidx = 0.d0
       t2kidx = 0.d0
       t34jidx = 0.d0
       t34kidx = 0.d0
       tottask = 0
       ttaskmax = 0.d0
       ttaskmin = 1.d24
       ttaskagg = 0.d0
       thalf0 = util_cpusec()
       ploop = 0
       next = -1
       if (ocoul) call ga_zero(g_coul)
       if (oexch) call ga_zero(g_exch)
c
c  Exchange part
c  Four fold shell loop
c
c
       if (oexch) then
         if (next.lt.0) next = moints_nxttask(num_nodes)

         do ii=1,nsh
           ish = int_mb(k_shmap+ii-1)
           status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
           ilen = ibfhi - ibflo + 1
           do jj=1,ii
             jsh = int_mb(k_shmap+jj-1)
             status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
             jlen = jbfhi - jbflo + 1
             scale0 = 1.d0
             schw_ij = schwarz_shell(ish,jsh)
             if (ish.eq.jsh) scale0 = 0.5d0
             if (next.eq.ploop) then
               ttask = util_cpusec()
               tottask = tottask + 1
               call dfill(n_ssni,0.d0,dbl_mb(k_ssni),1)
               do kgr=1,ngrp
                 kshlo = int_mb(k_glo+kgr-1)
                 kshhi = int_mb(k_ghi+kgr-1)
                 st = bas_cn2bfr(basis,int_mb(k_shmap+kshlo-1),iz,kz)
                 st = bas_cn2bfr(basis,int_mb(k_shmap+kshhi-1),kz,jz)
                 kbflo = int_mb(k_rbfmap+iz-1)
                 kbfhi = int_mb(k_rbfmap+jz-1)
                 kblen = kbfhi - kbflo + 1
                 do lgr=1,kgr
                   lshlo = int_mb(k_glo+lgr-1)
                   lshhi = int_mb(k_ghi+lgr-1)
                   st = bas_cn2bfr(basis,int_mb(k_shmap+lshlo-1),iz,kz)
                   st = bas_cn2bfr(basis,int_mb(k_shmap+lshhi-1),kz,jz)
                   lbflo = int_mb(k_rbfmap+iz-1)
                   lbfhi = int_mb(k_rbfmap+jz-1)
                   lblen = lbfhi - lbflo + 1
                   tz = util_cpusec()
                   scale = 1.d0
                   if (kgr.eq.lgr) scale = 0.5d0
                   call moints6x_gblkk( basis, ish, jsh, kshlo, kshhi,
     $                                  lshlo, lshhi, int_mb(k_shmap),
     $                                  int_mb(k_rbfmap),
     $                                  schw_ij, tol2e, osym, oblk,
     $                                  max2e, dbl_mb(k_eri),
     $                                  mem2, dbl_mb(k_iscr),
     $                                  ibflo, ibfhi, jbflo, jbfhi,
     $                                  kbflo, kbfhi, lbflo, lbfhi,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt) )
                   tint = tint + util_cpusec() - tz
                   flop1 = flop1 + 4*lblen*ilen*jlen*kblen*nmo1
                   tz = util_cpusec()
                   call moints6x_trf1K(nbf, qlo, qhi, mo1_lo, mo1_hi, 
     $                                 ilen, jlen, kbflo, kbfhi, lbflo,
     $                                 lbfhi, scale, dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt), dbl_mb(k_xmo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp))
                   t1kidx = t1kidx + util_cpusec() - tz
                 enddo
               enddo
               
               tz = util_cpusec()
               call moints6x_trf2K( nbf, qlo, qhi, 
     $                              occ_start, mo1_lo, mo1_hi,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              dbl_mb(k_ssni),
     $                              dbl_mb(k_hlp), dbl_mb(k_ssbb),
     $                              dbl_mb(k_hlp), dbl_mb(k_ssbb),
     $                              dbl_mb(k_xmo), g_exch )

                       
               t2kidx = t2kidx + util_cpusec() - tz
               next = moints_nxttask(num_nodes)
               ttask = util_cpusec() - ttask
               ttaskmax = max(ttaskmax,ttask)
               ttaskmin = min(ttaskmin,ttask)
               ttaskagg = ttaskagg + ttask
             endif
             ploop = ploop + 1
           enddo
         enddo
       endif
c
       if (oexch) then
         if (util_print('exchange half integral',print_debug)) then
           call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_exch)
         endif
       endif
c
c  Coulomb part
c  4-fold shell loop
c
       if (ocoul) then
         if (next.lt.0) next = moints_nxttask(num_nodes)

         do ii=1,nsh
           do jj=1,ii
           ish = max(int_mb(k_shmap+ii-1),int_mb(k_shmap+jj-1))
           jsh = min(int_mb(k_shmap+ii-1),int_mb(k_shmap+jj-1))
           status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
           status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
           ilen = ibfhi - ibflo + 1
           jlen = jbfhi - jbflo + 1
           schw_ij = schwarz_shell(ish,jsh)
           scale = 1.d0
           if (ish.eq.jsh) scale = scale*0.5d0
           odoit = schw_ij*schwarz_max().ge.tol2e
           if (odoit .and. osym) then
              odoit = sym_shell_pair(basis, ish, jsh, q2)
           endif

           if (odoit) then
             if (next.eq.ploop) then
               ttask = util_cpusec()
               tottask = tottask + 1
                 n_ijni = ilen*jlen*nbf*nmo1
                 call dfill(n_ijni,0.d0,dbl_mb(k_ssni),1)

               do kgr=1,ngrp
                 kshlo = int_mb(k_glo+kgr-1)
                 kshhi = int_mb(k_ghi+kgr-1)
                 st = bas_cn2bfr(basis,int_mb(k_shmap+kshlo-1),iz,kz)
                 st = bas_cn2bfr(basis,int_mb(k_shmap+kshhi-1),kz,jz)
                 kbflo = int_mb(k_rbfmap+iz-1)
                 kbfhi = int_mb(k_rbfmap+jz-1)
                 kblen = kbfhi - kbflo + 1
                 do lgr=1,kgr
                   lshlo = int_mb(k_glo+lgr-1)
                   lshhi = int_mb(k_ghi+lgr-1)
                   st = bas_cn2bfr(basis,int_mb(k_shmap+lshlo-1),iz,kz)
                   st = bas_cn2bfr(basis,int_mb(k_shmap+lshhi-1),kz,jz)
                   lbflo = int_mb(k_rbfmap+iz-1)
                   lbfhi = int_mb(k_rbfmap+jz-1)
                   lblen = lbfhi - lbflo + 1
                   tz = util_cpusec()
                   call moints_gblk( basis, ish, jsh,
     $                               kshlo, kshhi, lshlo, lshhi,
     $                               int_mb(k_shmap),int_mb(k_rbfmap),
     $                               schw_ij, tol2e, osym, oblk,
     $                               max2e, dbl_mb(k_eri),
     $                               mem2, dbl_mb(k_iscr),
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               dbl_mb(k_ssbb) )
                   tint = tint + util_cpusec() - tz
                   flop1 = flop1 + 4*lblen*ilen*jlen*kblen*nmo1
                   if (lgr.ne.kgr) then
                     tz = util_cpusec()
                     call moints_blktr(ilen, jlen, kblen, lblen,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt))
                     call moints_trf1(nbf, qlo, qhi, mo1_lo, mo1_hi, 
     $                                ilen, jlen, kbflo, kbfhi,
     $                                lbflo, lbfhi, 1.d0,
     $                                dbl_mb(k_ssbb),
     $                                dbl_mb(k_ssbbt),
     $                                dbl_mb(k_xmo),
     $                                dbl_mb(k_ssni),
     $                                dbl_mb(k_hlp))
                   else
                     tz = util_cpusec()
                     call moints_trf1(nbf, qlo, qhi, mo1_lo, mo1_hi,
     $                                ilen, jlen, kbflo, kbfhi, 
     $                                lbflo, lbfhi, 0.5d0, 
     $                                dbl_mb(k_ssbb),
     $                                dbl_mb(k_ssbb),
     $                                dbl_mb(k_xmo),
     $                                dbl_mb(k_ssni),
     $                                dbl_mb(k_hlp))
                   endif
                   t1jidx = t1jidx + util_cpusec() - tz
                 enddo
               enddo
                 flop2j = flop2j+2*nbf*ilen*jlen*((nmo1*(nmo1+1))/2)
                 tz = util_cpusec()
                 call moints_trf2J(nbf, qlo, qhi, 
     $                             occ_start, mo1_lo, mo1_hi,
     $                             ibflo, ibfhi, jbflo, jbfhi,
     $                             dbl_mb(k_ssni), dbl_mb(k_hlp),
     $                             dbl_mb(k_ssbb), dbl_mb(k_xmo),
     $                             g_coul)
                 t2jidx = t2jidx + util_cpusec() - tz
                 next = moints_nxttask(num_nodes)
                 ttask = util_cpusec() - ttask
                 ttaskmax = max(ttaskmax,ttask)
                 ttaskmin = min(ttaskmin,ttask)
                 ttaskagg = ttaskagg + ttask
               endif
               ploop = ploop + 1
             endif
           enddo
         enddo
       endif
c
       if (ocoul) then
         if (util_print('coulomb half integral',print_debug)) then
           call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_coul)
         endif
       endif
c
c Synchronize
c
       tz = util_cpusec()
       next = moints_nxttask(-num_nodes)
       tsynch = util_cpusec() - tz
       thalf = util_cpusec() - thalf0
       call moints_aodisk_tidy()
c
c Clean-up
c
       if (.not. ma_pop_stack(l_ssni))
     $     call errquit('moints: failed to pop', l_ssni, MA_ERR)
       if (.not. ma_pop_stack(l_hlp))
     $     call errquit('moints: failed to pop', l_hlp, MA_ERR)
       if (.not. ma_pop_stack(l_ssbbt))
     $     call errquit('moints: failed to pop', l_ssbbt, MA_ERR)
       if (.not. ma_pop_stack(l_ssbb))
     $     call errquit('moints: failed to pop', l_ssbb, MA_ERR)
       if (.not. ma_pop_stack(l_xmo))
     $     call errquit('moints: failed to pop', l_xmo, MA_ERR)
       if (.not. ma_pop_stack(l_rbfmap))
     $     call errquit('moints: failed to pop', l_rbfmap, MA_ERR)
       if (.not. ma_pop_stack(l_bfmap))
     $     call errquit('moints: failed to pop', l_bfmap, MA_ERR)
       if (.not. ma_pop_stack(l_ghi))
     $     call errquit('moints: failed to pop', l_ghi, MA_ERR)
       if (.not. ma_pop_stack(l_glo))
     $     call errquit('moints: failed to pop', l_glo, MA_ERR)
       if (.not. ma_pop_stack(l_shmap))
     $     call errquit('moints: failed to pop', l_shmap, MA_ERR)
       if (.not. ma_pop_stack(l_gmap))
     $     call errquit('moints: failed to pop', l_gmap, MA_ERR)
       if (.not. ma_pop_stack(l_iscr))
     $     call errquit('moints: failed to pop', l_iscr, MA_ERR)
       if (.not. ma_pop_stack(l_eri))
     $     call errquit('moints: failed to pop', l_eri, MA_ERR)
c
c 
c
       status = ma_push_get(MT_DBL,(nbf*nbf),'hlp',l_hlp,k_hlp)
c
c  Exchange 3 & 4 qtr transformation
c
       if (oexch) then
         tz = util_cpusec()
         call moints_Ktrf34(nbf, qlo, qhi, occ_start, mo1_lo, mo1_hi, 
     $                      mo2_lo, mo2_hi, .false., dbl_mb(k_mo),
     $                      dbl_mb(k_hlp), osym, int_mb(k_sym),
     $                      g_exch)
         t34kidx = util_cpusec() - tz
         if (util_print('exchange integral',print_debug)) then
           call moints_op_print(occ_start,mo1_lo,mo1_hi,nmo2,g_exch)
         endif
C         call MOINTS_OP_PRINT(OCC_START,MO1_LO,MO1_HI,NMO2,G_EXCH)
       endif
c
c  Coulomb 3 & 4 qtr transformation
c

       if (ocoul) then
         tz = util_cpusec()
         call moints_Jtrf34(nbf, qlo, qhi, occ_start, mo1_lo, mo1_hi,
     $                      mo2_lo, mo2_hi, dbl_mb(k_mo),
     $                      dbl_mb(k_hlp), osym, int_mb(k_sym),
     $                      g_coul)
         if (util_print('coulomb integral',print_debug)) then
           call moints_op_print(occ_start,mo1_lo,mo1_hi,nmo2,g_coul)
         endif
         t34jidx = util_cpusec() - tz
       endif
       call ga_sync()
c
c  Clean-up
c
       if (.not. ma_pop_stack(l_hlp))
     $     call errquit('moints: failed to pop', l_hlp, MA_ERR)
       if (.not. ma_pop_stack(l_mo))
     $     call errquit('moints: failed to pop', l_mo, MA_ERR)
       if (osym) then
          if (.not. ma_pop_stack(l_sym))
     $         call errquit('moints_2x: memory corrupt',0, MA_ERR)
       endif
c
c
c
#ifdef BLOCK_TRANSF
       if (ga_nodeid().eq.0) write(6,965)
 965   format(/,10x,'**** BLOCK TRANSFER ENABLED ****')
#endif
#ifdef NOCOMMS
       if (ga_nodeid().eq.0) write(6,334)
 334   format(/,10x,'**** COMMUNICATION DISABLED ****')
#endif


c
c Timings and Statistics bookkeeping
c
       flop1 = flop1*1.d-6
       flop2k = flop2k*1.d-6
       flop2j = flop2j*1.d-6
       ttotal = util_cpusec() - ttotal
       mi_npass = mi_npass + 1.d0
       mi_ttotal = mi_ttotal + ttotal
       mi_thalf = mi_thalf + thalf
       mi_tint = mi_tint + tint
       mi_t1j = mi_t1j + t1jidx
       mi_t2j = mi_t2j + t2jidx
       mi_t1k = mi_t1k + t1kidx
       mi_t2k = mi_t2k + t2kidx
       mi_t34k = mi_t34k + t34kidx
       mi_t34j = mi_t34j + t34jidx
       mi_flop1 = mi_flop1 + flop1
       mi_synch = mi_synch + tsynch
       mi_maxsynch = max(mi_maxsynch,tsynch)
       mi_minsynch = min(mi_minsynch,tsynch)
       mi_aggsynch = mi_aggsynch + tsynch
       mi_nsynchs = mi_nsynchs + 1
       mi_maxtask = max(mi_maxtask,ttaskmax)
       mi_mintask = min(mi_mintask,ttaskmin)
       mi_aggtask = mi_aggtask + ttaskagg
       mi_ntasks = mi_ntasks + tottask

       return
       end










c
c               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                 OLD VERSIONS WITH NO REORDERING
c
c               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c  Generate AO integrals in suitable form
c  exchange transformation:
c  Branch according to block or non-block
c  scheme
c
c

#ifdef NOCOMPILE
      subroutine moints6x_gblkk_old( basis, ish, jsh,
     $                           kshlo, kshhi, lshlo, lshhi,
     $                           schw_ij, tol2e, erilen, eri,
     $                           scrlen, iscr,
     $                           ibflo, ibfhi, jbflo, jbfhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           sbsb, sbsbt, oblk )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision sbsb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      double precision sbsbt(kblo:kbhi,lblo:lbhi,
     $                       jbflo:jbfhi,ibflo:ibfhi)
      logical oblk
c
      integer mxshq, imem
      integer k_tmp, l_tmp
      integer k_iv, k_jv, k_kv, k_lv
      integer k_il, k_jl, k_kl, k_ll
c


      if (.not.(oblk)) then
        call moints6x_gblkk_sq_old( basis, ish, jsh,
     $                          kshlo, kshhi, lshlo, lshhi,
     $                          schw_ij, tol2e, erilen, eri,
     $                          scrlen, iscr,
     $                          ibflo, ibfhi, jbflo, jbfhi,
     $                          kblo, kbhi, lblo, lbhi,
     $                          sbsb, sbsbt )
      else
        mxshq = (kshhi - kshlo + 1)*(lshhi - lshlo + 1)
        imem = 4*mxshq + 4*erilen
        if (.not. ma_push_get(MT_INT, imem, 'gblk tmp',l_tmp, k_tmp))
     $     call errquit('moints6x_gblkk: no memory for labels',imem,
     &       MA_ERR)
        k_iv = k_tmp
        k_jv = k_tmp + mxshq
        k_kv = k_tmp + 2*mxshq
        k_lv = k_tmp + 3*mxshq
        k_il = k_tmp + 4*mxshq
        k_jl = k_tmp + 4*mxshq + erilen
        k_kl = k_tmp + 4*mxshq + 2*erilen
        k_ll = k_tmp + 4*mxshq + 3*erilen
        call moints6x_gblkk_mq_old( basis, ish, jsh,
     $                          kshlo, kshhi, lshlo, lshhi,
     $                          schw_ij, tol2e, erilen, eri,
     $                          scrlen, iscr,
     $                          ibflo, ibfhi, jbflo, jbfhi,
     $                          kblo, kbhi, lblo, lbhi,
     $                          sbsb, sbsbt, mxshq,
     $                          int_mb(k_iv), int_mb(k_jv),
     $                          int_mb(k_kv), int_mb(k_lv),
     $                          int_mb(k_il), int_mb(k_jl),
     $                          int_mb(k_kl), int_mb(k_ll) )
        if (.not.ma_pop_stack(l_tmp))
     $    call errquit('moints6x_gblkk: failed to pop',0, MA_ERR)
      endif
      return
      end





c
c This version uses single shell 
c quartet AO integrals.
c
      subroutine moints6x_gblkk_sq_old( basis, ish, jsh,
     $                              kshlo, kshhi, lshlo, lshhi,
     $                              schw_ij, tol2e, erilen, eri,
     $                              scrlen, iscr,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              kblo, kbhi, lblo, lbhi,
     $                              sbsb, sbsbt )
      implicit none
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision sbsb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      double precision sbsbt(kblo:kbhi,lblo:lbhi,
     $                       jbflo:jbfhi,ibflo:ibfhi)
c
      integer ilen, jlen
      integer ksh, lsh, kbflo, kbfhi, lbflo, lbfhi, ltop
      integer klen, llen, kblen, lblen, bsize, k1, k2, j, i
      logical status
c
c
c
      ilen = ibfhi - ibflo + 1
      jlen = jbfhi - jbflo + 1
      kblen = kbhi - kblo + 1
      lblen = lbhi - lblo + 1
      bsize = kblen*lblen*ilen*jlen
      call dfill(bsize,0.d0,sbsb,1)
      call dfill(bsize,0.d0,sbsbt,1)
      do ksh=kshlo,kshhi
        status = bas_cn2bfr(basis,ksh,kbflo,kbfhi)
        klen = kbfhi - kbflo + 1
        ltop = lshhi
        if (kshlo.eq.lshlo) ltop = ksh
        do lsh=lshlo,ltop
C        do lsh=lshlo,lshhi
          status = bas_cn2bfr(basis,lsh,lbflo,lbfhi)
          llen = lbfhi - lbflo + 1
          if (schwarz_shell(ish,ksh)*
     $        schwarz_shell(jsh,lsh).ge.tol2e) then
            call int_2e4c(basis, ish, ksh, basis, jsh, lsh,
     $                    scrlen, iscr, erilen, eri )
            call moints6x_eri2blkk(ilen, jlen, klen, llen,
     $                             eri, sbsb(lbflo,kbflo,jbflo,ibflo),
     $                             lblen, kblen )
          endif
          if (schwarz_shell(ish,lsh)*
     $        schwarz_shell(jsh,ksh).ge.tol2e) then
            call int_2e4c(basis, ish, lsh, basis, jsh, ksh,
     $                    scrlen, iscr, erilen, eri )
            call moints6x_eri2blkk(ilen, jlen, llen, klen,
     $                             eri, sbsbt(kbflo,lbflo,jbflo,ibflo),
     $                             kblen, lblen )
          endif
        enddo
      enddo
      if (kshlo.eq.lshlo) then
        do i=ibflo,ibfhi
          do j=jbflo,jbfhi
            do k1=kblo,kbhi
              do k2=kblo,k1-1
                sbsb(k1,k2,j,i) = sbsbt(k1,k2,j,i)
                sbsbt(k2,k1,j,i) = sbsb(k2,k1,j,i)
              enddo
            enddo
          enddo
        enddo
      endif
c
c
c
      return
      end






c
c  This uses multiple quartets at a time
c  AO integrals with labels
c
c
c
      subroutine moints6x_gblkk_mq_old( basis, ish, jsh,
     $                              kshlo, kshhi, lshlo, lshhi,
     $                              schw_ij, tol2e, erilen, eri,
     $                              scrlen, iscr,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              kblo, kbhi, lblo, lbhi,
     $                              sbsb, sbsbt, mxshq,
     $                              ishv, jshv, kshv, lshv,
     $                              ilab, jlab, klab, llab )
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision sbsb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      double precision sbsbt(kblo:kbhi,lblo:lbhi,
     $                       jbflo:jbfhi,ibflo:ibfhi)
      integer mxshq
      integer ishv(mxshq), jshv(mxshq), kshv(mxshq), lshv(mxshq)
      integer ilab(erilen), jlab(erilen), klab(erilen), llab(erilen)
c
      integer bsize, ieri, neri, neritotal, nshq
      integer ii, jj, kk, ll
      double precision beffect, q4
      logical use_q4, more
c
      logical intb_2e4c, intb_init4c
      external intb_2e4c, intb_init4c
c
      data use_q4/.false./
      data q4/1.d0/
c
c
      neritotal = 0
      bsize = (kbhi - kblo + 1)*(lbhi - lblo + 1)*
     $        (ibfhi - ibflo + 1)*(jbfhi - jbflo + 1)
      call dfill(bsize,0.d0,sbsb,1)
      call dfill(bsize,0.d0,sbsbt,1)
c
c
c  First time around   ( I K | J L ) type
c
c  Fill out shell arrays with interacting quartet labels
c        
      nshq = 0
      do kk=kshlo,kshhi
        do ll=lshlo,lshhi
          if ((schwarz_shell(ish,kk)*schwarz_shell(jsh,ll))
     $                .ge.tol2e) then
            if (nshq.eq.mxshq)
     $         call errquit('moints6x_gblkk: internal error',0,
     &       UNKNOWN_ERR)
            nshq = nshq + 1
            kshv(nshq) = kk
            lshv(nshq) = ll
          endif
        enddo
      enddo
      call ifill( nshq, ish, ishv, 1 )
      call ifill( nshq, jsh, jshv, 1 )
c
c  Prepare texas for shell block
c
      if (.not. intb_init4c( basis, ishv, kshv, basis, jshv, lshv,
     $     nshq, q4, use_q4, 
     $     scrlen, iscr, erilen, beffect ))
     $   call errquit('moints6x_gblkk: intb_init4c failed',0,
     &       INT_ERR)
c
c  Get some batch of integrals
c
 100  more = intb_2e4c( basis, ishv, kshv, basis, jshv, lshv,
     $                  nshq, q4, use_q4, tol2e, .false.,
     $                  ilab, klab, jlab, llab,
     $                  eri, erilen, neri, scrlen, iscr )
c
c  Unpack labels and assign to integrals
c
      if (neri.gt.0) then
        neritotal = neritotal + neri
        do ieri=1,neri
          ii = ilab(ieri)
          kk = klab(ieri)
          jj = jlab(ieri)
          ll = llab(ieri)
          sbsb(ll,kk,jj,ii) = eri(ieri)
        enddo
      endif
      if (more) goto 100
c
c
c  Second time around   ( I L | J K ) type
c
c  Fill out shell arrays with interacting quartet labels
c        
      nshq = 0
      do kk=kshlo,kshhi
        do ll=lshlo,lshhi
          if ((schwarz_shell(ish,ll)*schwarz_shell(jsh,kk))
     $                .ge.tol2e) then
            if (nshq.eq.mxshq)
     $         call errquit('moints6x_gblkk: internal error',0, INT_ERR)
            nshq = nshq + 1
            kshv(nshq) = kk
            lshv(nshq) = ll
          endif
        enddo
      enddo
      call ifill( nshq, ish, ishv, 1 )
      call ifill( nshq, jsh, jshv, 1 )
c
c  Prepare texas for shell block
c
      if (.not. intb_init4c( basis, ishv, lshv, basis, jshv, kshv,
     $     nshq, q4, use_q4,
     $     scrlen, iscr, erilen, beffect ))
     $   call errquit('moints6x_gblkk: intb_init4c failed',0, INT_ERR)
c
c  Get some batch of integrals
c
 200  more = intb_2e4c( basis, ishv, lshv, basis, jshv, kshv,
     $                  nshq, q4, use_q4, tol2e, .false.,
     $                  ilab, llab, jlab, klab,
     $                  eri, erilen, neri, scrlen, iscr )
c
c  Unpack labels and assign to integrals
c
      if (neri.gt.0) then
        neritotal = neritotal + neri
        do ieri=1,neri
          ii = ilab(ieri)
          ll = llab(ieri)
          jj = jlab(ieri)
          kk = klab(ieri)
          sbsbt(kk,ll,jj,ii) = eri(ieri)
        enddo
      endif
      if (more) goto 200
c
c
c
      return
      end






      subroutine moints6x_eri2blkk( ilen, jlen, klen, llen,
     $                              eri, blk, lblen, kblen )
      implicit none
      integer ilen, jlen, klen, llen, lblen, kblen
      double precision blk(lblen, kblen, jlen, ilen )
      double precision eri(llen, jlen, klen, ilen )
      integer k,l,i,j

      do i=1,ilen
        do j=1,jlen
          do k=1,klen
            do l=1,llen
              blk(l,k,j,i) = eri(l,j,k,i)
            enddo
          enddo
        enddo
      enddo
      
      return
      end

#endif















c
c               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                 NEW VERSIONS WITH SHELL REORDERING
c
c               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c  Generate AO integrals in suitable form
c  exchange transformation:
c  Branch according to block or non-block
c  scheme
c
c
      subroutine moints6x_gblkk( basis, ish, jsh, kshlo, kshhi,
     $                           lshlo, lshhi, shmap, rbfmap,
     $                           schw_ij, tol2e, osym, oblk,
     $                           erilen, eri, scrlen, iscr,
     $                           ibflo, ibfhi, jbflo, jbfhi,
     $                           kblo, kbhi, lblo, lbhi,
     $                           sbsb, sbst )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer shmap(*), rbfmap(*)
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      logical osym, oblk
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision sbsb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      double precision sbst(kblo:kbhi,lblo:lbhi,jbflo:jbfhi,ibflo:ibfhi)
c
c
      integer mxshq, imem
      integer k_tmp, l_tmp
      integer k_iv, k_jv, k_kv, k_lv
      integer k_il, k_jl, k_kl, k_ll
      integer blklen
      logical status
      integer blkid, blkidt
c
      logical moints_gblk_fromdisk, moints_gblk_todisk
      external moints_gblk_fromdisk, moints_gblk_todisk
c
      blklen = (ibfhi-ibflo+1)*(jbfhi-jbflo+1)*
     $         (kbhi-kblo+1)*(lbhi-lblo+1)
c
c  Check to see if we can retrieve from disk
c  For diagonal blocks, we just need to copy since (ik|jl) = (jl|ik)
c
      blkid = kshlo*1000 + lshlo
      blkidt = (kshlo*1000 + lshlo)*1000
      if (moints_gblk_fromdisk( blkid, ish, jsh, kshlo, lshlo,
     $                          ibflo, ibfhi, jbflo, jbfhi, 
     $                          kblo, kbhi, lblo, lbhi,
     $                          sbsb )) then
        if (kshlo.ne.lshlo) then
          status = moints_gblk_fromdisk( blkidt, ish, jsh, 
     $                                   kshlo, lshlo, 
     $                                   ibflo, ibfhi, jbflo, jbfhi, 
     $                                   lblo, lbhi, kblo, kbhi,
     $                                   sbst )
          if (.not.(status)) call errquit(
     $        'moints_6x_gblkk: inconsistent AO disk records',0,
     &       INT_ERR)
        else
          call dcopy( blklen, sbsb, 1, sbst, 1 )
        endif
        return
      endif
c
c  Otherwise compute directly
c
      if (oblk) then
        mxshq = (kshhi - kshlo + 1)*(lshhi - lshlo + 1)
        imem = 4*mxshq + 4*erilen
        if (.not. ma_push_get(MT_INT, imem, 'gblk tmp',l_tmp, k_tmp))
     $     call errquit('moints6x_gblkk: no memory for labels',imem,
     &       INT_ERR)
        k_iv = k_tmp
        k_jv = k_tmp + mxshq
        k_kv = k_tmp + 2*mxshq
        k_lv = k_tmp + 3*mxshq
        k_il = k_tmp + 4*mxshq
        k_jl = k_tmp + 4*mxshq + erilen
        k_kl = k_tmp + 4*mxshq + 2*erilen
        k_ll = k_tmp + 4*mxshq + 3*erilen
        call moints6x_gblkk_mq( basis, ish, jsh,
     $                          kshlo, kshhi, lshlo, lshhi,
     $                          shmap, rbfmap,
     $                          schw_ij, tol2e, erilen, eri,
     $                          scrlen, iscr,
     $                          ibflo, ibfhi, jbflo, jbfhi,
     $                          kblo, kbhi, lblo, lbhi,
     $                          sbsb, sbst, mxshq,
     $                          int_mb(k_iv), int_mb(k_jv),
     $                          int_mb(k_kv), int_mb(k_lv),
     $                          int_mb(k_il), int_mb(k_jl),
     $                          int_mb(k_kl), int_mb(k_ll) )
        if (.not.ma_pop_stack(l_tmp))
     $    call errquit('moints6x_gblkk: failed to pop',0, MA_ERR)
      else
        call moints6x_gblkk_sq( basis, ish, jsh,
     $                          kshlo, kshhi, lshlo, lshhi,
     $                          shmap, rbfmap,
     $                          schw_ij, tol2e, erilen, eri,
     $                          scrlen, iscr,
     $                          ibflo, ibfhi, jbflo, jbfhi,
     $                          kblo, kbhi, lblo, lbhi,
     $                          sbsb, sbst )
      endif
c
c  Cache this block to disk if required
c  Only one copy for diagonal blocks
c     
      status = moints_gblk_todisk( blkid, ish, jsh, kshlo, lshlo,
     $                             ibflo, ibfhi, jbflo, jbfhi,
     $                             kblo, kbhi, lblo, lbhi,
     $                             sbsb )
      if (kshlo.ne.lshlo) then
        status = moints_gblk_todisk( blkidt, ish, jsh, kshlo, lshlo,
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               lblo, lbhi, kblo, kbhi,
     $                               sbst )
      endif


      return
      end





c
c This version uses single shell 
c quartet AO integrals.
c
      subroutine moints6x_gblkk_sq( basis, ish, jsh,
     $                              kshlo, kshhi, lshlo, lshhi,
     $                              shmap, rbfmap,
     $                              schw_ij, tol2e, erilen, eri,
     $                              scrlen, iscr,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              kblo, kbhi, lblo, lbhi,
     $                              sbsb, sbst )
      implicit none
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer shmap(*), rbfmap(*)
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision sbsb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      double precision sbst(kblo:kbhi,lblo:lbhi,jbflo:jbfhi,ibflo:ibfhi)
c
      integer ilen, jlen
      integer ksh, lsh, kbflo, kbfhi, lbflo, lbfhi, ltop, kk, ll
      integer klen, llen, kblen, lblen, bsize, k1, k2, j, i
      logical status
c
c
c
      ilen = ibfhi - ibflo + 1
      jlen = jbfhi - jbflo + 1
      kblen = kbhi - kblo + 1
      lblen = lbhi - lblo + 1
      bsize = kblen*lblen*ilen*jlen
      call dfill(bsize,0.d0,sbsb,1)
      call dfill(bsize,0.d0,sbst,1)
c
      do kk=kshlo,kshhi
        ksh = shmap(kk)
        status = bas_cn2bfr(basis,ksh,kbflo,kbfhi)
        klen = kbfhi - kbflo + 1
        ltop = lshhi
        if (kshlo.eq.lshlo) ltop = kk
        do ll=lshlo,ltop
          lsh = shmap(ll)
          status = bas_cn2bfr(basis,lsh,lbflo,lbfhi)
          llen = lbfhi - lbflo + 1
          if (schwarz_shell(ish,ksh)*
     $        schwarz_shell(jsh,lsh).ge.tol2e) then
            call int_2e4c(basis, ish, ksh, basis, jsh, lsh,
     $                    scrlen, iscr, erilen, eri )
            call moints6x_eri2blkkx( rbfmap, ilen, jlen,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               kblo, kbhi, lblo, lbhi,
     $                               eri, sbsb )
          endif
          if (schwarz_shell(ish,lsh)*
     $        schwarz_shell(jsh,ksh).ge.tol2e) then
            call int_2e4c(basis, ish, lsh, basis, jsh, ksh,
     $                    scrlen, iscr, erilen, eri )
            call moints6x_eri2blkkx( rbfmap, ilen, jlen,
     $                               lbflo, lbfhi, kbflo, kbfhi,
     $                               lblo, lbhi, kblo, kbhi,
     $                               eri, sbst )
          endif
        enddo
      enddo
c
      if (kshlo.eq.lshlo) then
        do i=ibflo,ibfhi
          do j=jbflo,jbfhi
            do k1=kblo,kbhi
              do k2=kblo,k1-1
                sbsb(k1,k2,j,i) = sbst(k1,k2,j,i)
                sbst(k2,k1,j,i) = sbsb(k2,k1,j,i)
              enddo
            enddo
          enddo
        enddo
      endif
c
      return
      end






c
c  This uses multiple quartets at a time
c  AO integrals with labels
c
      subroutine moints6x_gblkk_mq( basis, ish, jsh,
     $                              kshlo, kshhi, lshlo, lshhi,
     $                              shmap, rbfmap,
     $                              schw_ij, tol2e, erilen, eri,
     $                              scrlen, iscr,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              kblo, kbhi, lblo, lbhi,
     $                              sbsb, sbst, mxshq,
     $                              ishv, jshv, kshv, lshv,
     $                              ilab, jlab, klab, llab )
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer shmap(*), rbfmap(*)
      integer erilen, scrlen
      double precision schw_ij, tol2e, eri(erilen), iscr(scrlen)
      integer ibflo, ibfhi, jbflo, jbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision sbsb(lblo:lbhi,kblo:kbhi,jbflo:jbfhi,ibflo:ibfhi)
      double precision sbst(kblo:kbhi,lblo:lbhi,jbflo:jbfhi,ibflo:ibfhi)
      integer mxshq
      integer ishv(mxshq), jshv(mxshq), kshv(mxshq), lshv(mxshq)
      integer ilab(erilen), jlab(erilen), klab(erilen), llab(erilen)
c
      integer bsize, ieri, neri, neritotal, nshq
      integer ii, jj, kk, ll, ksh, lsh
      double precision beffect, q4
      logical use_q4, more
c
      logical intb_2e4c, intb_init4c
      external intb_2e4c, intb_init4c
c
      data use_q4/.false./
      data q4/1.d0/
c
c
      neritotal = 0
      bsize = (kbhi - kblo + 1)*(lbhi - lblo + 1)*
     $        (ibfhi - ibflo + 1)*(jbfhi - jbflo + 1)
      call dfill(bsize,0.d0,sbsb,1)
      call dfill(bsize,0.d0,sbst,1)
c
c
c  First time around   ( I K | J L ) type
c
c  Fill out shell arrays with interacting quartet labels
c        
      nshq = 0
      do kk=kshlo,kshhi
        ksh = shmap(kk)
        do ll=lshlo,lshhi
          lsh = shmap(ll)
          if ((schwarz_shell(ish,ksh)*schwarz_shell(jsh,lsh))
     $                .ge.tol2e) then
            if (nshq.eq.mxshq)
     $         call errquit('moints6x_gblkk: internal error',0, INT_ERR)
            nshq = nshq + 1
            kshv(nshq) = ksh
            lshv(nshq) = lsh
          endif
        enddo
      enddo
      call ifill( nshq, ish, ishv, 1 )
      call ifill( nshq, jsh, jshv, 1 )
c
c  Prepare texas for shell block
c
      if (.not. intb_init4c( basis, ishv, kshv, basis, jshv, lshv,
     $     nshq, q4, use_q4,
     $     scrlen, iscr, erilen, beffect ))
     $   call errquit('moints6x_gblkk: intb_init4c failed',0, INT_ERR)
c
c  Get some batch of integrals
c
 100  more = intb_2e4c( basis, ishv, kshv, basis, jshv, lshv,
     $                  nshq, q4, use_q4, tol2e, .false.,
     $                  ilab, klab, jlab, llab,
     $                  eri, erilen, neri, scrlen, iscr )
c
c  Unpack labels and assign to integrals
c
      if (neri.gt.0) then
        neritotal = neritotal + neri
        do ieri=1,neri
          ii = ilab(ieri)
          jj = jlab(ieri)
          kk = rbfmap(klab(ieri))
          ll = rbfmap(llab(ieri))
          sbsb(ll,kk,jj,ii) = eri(ieri)
        enddo
      endif
      if (more) goto 100
c
c
c  Second time around   ( I L | J K ) type
c
c  Fill out shell arrays with interacting quartet labels
c        
      nshq = 0
      do kk=kshlo,kshhi
        ksh = shmap(kk)
        do ll=lshlo,lshhi
          lsh = shmap(ll)
          if ((schwarz_shell(ish,lsh)*schwarz_shell(jsh,ksh))
     $                .ge.tol2e) then
            if (nshq.eq.mxshq)
     $         call errquit('moints6x_gblkk: internal error',0, INT_ERR)
            nshq = nshq + 1
            kshv(nshq) = ksh
            lshv(nshq) = lsh
          endif
        enddo
      enddo
      call ifill( nshq, ish, ishv, 1 )
      call ifill( nshq, jsh, jshv, 1 )
c
c  Prepare texas for shell block
c
      if (.not. intb_init4c( basis, ishv, lshv, basis, jshv, kshv,
     $     nshq, q4, use_q4, 
     $     scrlen, iscr, erilen, beffect ))
     $   call errquit('moints6x_gblkk: intb_init4c failed',0, INT_ERR)
c
c  Get some batch of integrals
c
 200  more = intb_2e4c( basis, ishv, lshv, basis, jshv, kshv,
     $                  nshq, q4, use_q4, tol2e, .false.,
     $                  ilab, llab, jlab, klab,
     $                  eri, erilen, neri, scrlen, iscr )
c
c  Unpack labels and assign to integrals
c
      if (neri.gt.0) then
        neritotal = neritotal + neri
        do ieri=1,neri
          ii = ilab(ieri)
          jj = jlab(ieri)
          ll = rbfmap(llab(ieri))
          kk = rbfmap(klab(ieri))
          sbst(kk,ll,jj,ii) = eri(ieri)
        enddo
      endif
      if (more) goto 200
c
c
c
      return
      end






      subroutine moints6x_eri2blkkx( rbfmap, ilen, jlen,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               kblo, kbhi, lblo, lbhi,
     $                               eri, blk )
      implicit none
      integer rbfmap(*)
      integer ilen, jlen
      integer kbflo, kbfhi, lbflo, lbfhi
      integer kblo, kbhi, lblo, lbhi
      double precision blk(lblo:lbhi,kblo:kbhi,jlen,ilen)
      double precision eri(lbflo:lbfhi,jlen,kbflo:kbfhi,ilen)
      integer kk,ll,k,l,i,j

      do i=1,ilen
        do j=1,jlen
          do k=kbflo,kbfhi
            kk = rbfmap(k)
            do l=lbflo,lbfhi
              ll = rbfmap(l)
              blk(ll,kk,j,i) = eri(l,j,k,i)
            enddo
          enddo
        enddo
      enddo
      
      return
      end














      subroutine moints6x_pushK( nbf, ilen, jlen, klo, khi,
     $                           nmo, x, ssni )
      implicit none
      integer nbf, ilen, jlen, klo, khi, nmo
      double precision x(klo:khi,jlen,ilen,nmo)
      double precision ssni(nbf,jlen,ilen,nmo)
      integer a,i,j,k
      
      do a=1,nmo
        do i=1,ilen
          do j=1,jlen
            do k=klo,khi
              ssni(k,j,i,a) = ssni(k,j,i,a) + x(k,j,i,a)
            enddo
          enddo
        enddo
      enddo
      return
      end







      subroutine moints6x_trf1K(nbf, qlo, qhi, molo, mohi, 
     $                          ilen, jlen, klo, khi, llo, lhi, 
     $                          scale, ssbb, ssbbt, c, ssni, hlp )
      implicit none
      integer nbf, qlo, qhi
      integer molo, mohi, ilen, jlen, klo, khi, llo, lhi
      double precision scale
      double precision ssbb(llo:lhi,klo:khi,jlen,ilen)
      double precision ssbbt(klo:khi,llo:lhi,jlen,ilen)
      double precision c(nbf,qlo:qhi)
      double precision ssni(nbf,jlen,ilen,molo:mohi)
      double precision hlp(*)
c
c
      integer nmo, llen, klen, kjilen, ljilen
c
c      
      nmo = mohi - molo + 1
      llen = lhi - llo + 1
      klen = khi - klo + 1
      kjilen = klen*jlen*ilen
      ljilen = llen*jlen*ilen
      call dgemm( 't', 'n', kjilen, nmo, llen, scale, ssbb, llen,
     $            c(llo,molo), nbf, 0.d0, hlp, kjilen )
      call moints6x_pushK( nbf, ilen, jlen, klo, khi,
     $                     nmo, hlp, ssni )
      call dgemm( 't', 'n', ljilen, nmo, klen, scale, ssbbt, klen,
     $            c(klo,molo), nbf, 0.d0, hlp, ljilen )
      call moints6x_pushK( nbf, ilen, jlen, llo, lhi,
     $                     nmo, hlp, ssni )

      return
      end









      subroutine moints6x_trf2K(nbf, qlo, qhi, ostart, olo, ohi, 
     $                          ilo, ihi, jlo, jhi, ssni, h1, h2, 
     $                          g1, g2, c, g_exch )
      implicit none
      integer nbf, qlo, qhi
      integer ostart, olo, ohi, ilo, ihi, jlo, jhi
      double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
#ifdef SCATTER_TRANSF
      double precision h1(*)
      integer iv(100),jv(100)
#endif
      double precision h1(nbf,ilo:ihi)
      double precision h2(jlo:jhi,ilo:ihi)
      double precision g1(nbf,jlo:jhi)
      double precision g2(ilo:ihi,jlo:jhi)
      double precision c(nbf,qlo:qhi)
      integer g_exch
c
      integer ilen, jlen, ijlen
      integer nni, nnj, ijlo, ijhi, jilo, jihi
      integer aoff, ofroz, aa, bb, a, b, ab, i, j
c
      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      ijlen = ilen*jlen
      nni = ilen*nbf
      nnj = jlen*nbf

#ifdef BLOCK_TRANSF
      ijlo = (ilo-1)*nbf + 1
      ijhi = ihi*nbf
#endif
      do a=olo,ohi
        do b=ostart,a
          call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,b),
     $               nbf,c(1,a),nbf,0.d0,h2,ijlen)
#ifndef NOCOMMS
          aa = a - ofroz
          bb = b - ofroz
          ab = (aa*(aa-1))/2 + bb - aoff
#ifdef BLOCK_TRANSF
          call dfill(nni,0.d0,h1,1)
          do i=ilo,ihi
            do j=jlo,jhi
              h1(j,i) = h2(j,i)
            enddo
          enddo
          call ga_acc(g_exch,ijlo,ijhi,ab,ab,h1,nni,1.d0)
#else
#if SCATTER_TRANSF
          ij = 0
          do i=ilo,ihi
            ii = (i-1)*nbf
            do j=jlo,jhi
              ij = ij + 1
              h1(ij) = h2(j,i)
              iv(ij) = ii + j
              jv(ij) = ab
            enddo
          enddo
          call ga_dscatter(g_exch,h1,iv,jv,ijlen)
#else
          do i=ilo,ihi
            ijlo = (i-1)*nbf + jlo
            ijhi = (i-1)*nbf + jhi
            call ga_acc(g_exch,ijlo,ijhi,ab,ab,h2(jlo,i),ilen,1.d0)
          enddo
#endif
#endif
#endif
        enddo
      enddo



      if (ilo.eq.jlo) goto 120
#ifdef BLOCK_TRANSF
      jilo = (jlo-1)*nbf + 1
      jihi = jhi*nbf
#endif
      do a=olo,ohi
        do b=ostart,a
          call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,a),
     $               nbf,c(1,b),nbf,0.d0,h2,ijlen)
#ifndef NOCOMMS
          aa = a - ofroz
          bb = b - ofroz
          ab = (aa*(aa-1))/2 + bb - aoff
#ifdef BLOCK_TRANSF
          call dfill(nnj,0.d0,g1,1)
          do j=jlo,jhi
            do i=ilo,ihi
              g1(i,j) = h2(j,i)
            enddo
          enddo
          call ga_acc(g_exch,jilo,jihi,ab,ab,g1,nnj,1.d0)
#else
#if SCATTER_TRANSF
          ij = 0
          do j=jlo,jhi
            jj = (j-1)*nbf
            do i=ilo,ihi
              ij = ij + 1
              h1(ij) = h2(j,i)
              iv(ij) = jj + i
              jv(ij) = ab
            enddo
          enddo
          call ga_dscatter(g_exch,h1,iv,jv,ijlen)
#else
          do j=jlo,jhi
            do i=ilo,ihi
              g1(i,j) = h2(j,i)
            enddo
          enddo
          do j=jlo,jhi
            jilo = (j-1)*nbf + ilo
            jihi = (j-1)*nbf + ihi
            call ga_acc(g_exch,jilo,jihi,ab,ab,g1(ilo,j),ilen,1.d0)
          enddo
#endif
#endif
#endif
        enddo
      enddo

 120  continue
      return
      end
