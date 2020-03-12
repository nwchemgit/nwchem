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
c
c
c
       subroutine moints_build_2x(basis, ohalf, ousesym,
     $                            occ_start, mo1_lo, mo1_hi,
     $                            mo2_lo, mo2_hi,
     $                            g_movecs,
     $                            g_coul, ocoul,
     $                            g_exch, oexch,
     $                            blen, oblk )
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "util.fh"
#include "schwarz.fh"
#include "msgids.fh"
#include "sym.fh"
c
c  Arguments
c
       integer basis                          ! Basis handle
       logical ohalf                          ! Toggle alternative algorithm for half-ints intermediates
       logical ousesym                        ! Use symmetry if present
       integer occ_start                      ! Offset from frozen core
       integer mo1_lo, mo1_hi                 ! 1st Pair Index range
       integer mo2_lo, mo2_hi                 ! 2nd Pair Index range
       integer g_movecs                       ! MO coefficients
       integer g_coul                         ! Coulomb operator
       integer g_exch                         ! Exchange operator
       logical ocoul,oexch                    ! Type selection
       integer blen                           ! Blocking length
       logical oblk                           ! Toggle blocked AO ints
c
c  Local variables
c

       integer nmo, nmo1, nmo2, nbf, nsh, maxbfsh
       integer qlo, qhi
       integer bsize, ngrp
       integer ii, jj, ish, jsh, ilen, jlen
       integer ibflo, ibfhi, jbflo, jbfhi
       integer kbflo, kbfhi, lbflo, lbfhi
       integer kshlo, kshhi, lshlo, lshhi
       integer kblen, lblen
       integer kgr, lgr
       integer i
       integer iz, jz, kz
       integer l_ssbb, k_ssbb, l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni, l_hlp2, k_hlp2
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_mo, k_mo, l_xmo, k_xmo
       integer l_shmap, k_shmap, l_bfmap, k_bfmap, l_rbfmap, k_rbfmap
       integer l_glo, k_glo, l_ghi, k_ghi
       integer k_sym, l_sym
       integer n_ssbb, n_ssbb1, n_ssni, n_ssai, n_hlp, n_hlp2, n_ijni
       integer mem2, max2e
       integer num_nodes, ploop, next
       integer nmixed
       double precision tol2e, scale, schw_ij
*       double precision integ_acc
       double precision tz, thalf, tint, t1idx, t2kidx, t3kidx, t4kidx
       double precision t34jidx, ttotal, t2jidx, tsynch
       double precision ttask, ttaskmax, ttaskmin, ttaskagg
       double precision flop1, q2
       integer tottask
       integer geom
       logical status, osym, st
       logical odoit,doloop
c
#include "moints_stats.fh"
#include "moints.fh"
c
       integer clo, chi, rlo, rhi,nodesloop
c      integer moints_numgr, gr_len, moints_nxttask
c      external moints_numgr, gr_len, moints_nxttask
       integer  gr_len
       external gr_len
       data k_sym/0/
c
c     Set accuracy of integrals from required precision
c     ... this should come via the argument list but that can wait
c
       tol2e = 1d-16             ! No database !!!!
*       integ_acc = min(1d-9, max(0.01d0*tol2e,1d-17))
*       call int_acc_set(integ_acc)
c
c  General basis info
c
       ttotal = util_cpusec()
       num_nodes = ga_nnodes()
       if (.not. bas_geom(basis, geom)) call errquit
     $      ('moints: cannot get geometry', 0, BASIS_ERR)
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not.status) call errquit('moints: cannot get basis info',0,
     &       BASIS_ERR)
       qlo  = min(occ_start,mo2_lo,mo1_lo)    ! min MO index
       qhi  = max(mo1_hi,mo2_hi)              ! max MO index
       nmo  = qhi - qlo + 1          
       nmo1 = mo1_hi - mo1_lo + 1
       nmo2 = mo2_hi - mo2_lo + 1
c     
c  Symmetry adapt the MOs and renumber irreps to start at zero
c     
       osym = (ousesym.and.(sym_number_ops(geom).gt.0)) .and.
     $      sym_abelian_group(geom)
       if (osym) then
          if (.not. ma_push_get(MT_INT, nbf, 'movec syms',l_sym, k_sym))
     $         call errquit('moints_2x: no memory for syms?',0,
     &       MA_ERR)
          call sym_movecs_adapt(basis, tol2e, g_movecs, int_mb(k_sym),
     $                          nmixed)
          if (nmixed.ne.0 .and. ga_nodeid().eq.0) then
             write(6,*) ' !! warning: moints 2x: ',nmixed,
     $            ' MO vectors were symmetry contaminated'
             call util_flush(6)
          endif
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
c  Perpare AO disk caching
c       
       status = status .and. moints_aodisk_prep()
c
c  Local MO coefficients
c  Subset of columns                    ( qlo, ...., qhi )
c  excludes frozen core & virtual
c
       status = status.and.ma_push_get(MT_DBL,(nbf*nmo),
     $                                 'movecs cols',l_mo,k_mo)
       call ga_sync()
#ifdef CRAYXT
       if(ga_nodeid().eq.0)
     G      call ga_get(g_movecs, 1, nbf, qlo, qhi, dbl_mb(k_mo), nbf)
       call ga_brdcst(65537, dbl_mb(k_mo), 
     G      MA_sizeof(MT_DBL,1,MT_BYTE)*nbf*(qhi-qlo+1), 0)
#else
       call ga_get(g_movecs, 1, nbf, qlo, qhi, dbl_mb(k_mo), nbf)
#endif
c
c  AO ints allocation
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
       n_hlp2 = max((maxbfsh*maxbfsh*nmo2),(maxbfsh*nbf))
       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
       n_ssai = maxbfsh*maxbfsh*nmo2*nmo1
       status = ma_push_get(MT_DBL,n_ssbb1,'ssbb block',l_ssbb,k_ssbb)
       status = ma_push_get(MT_DBL,n_ssbb,'ssbbt block',l_ssbbt,k_ssbbt)
       status = ma_push_get(MT_DBL,n_hlp,'hlp block',l_hlp,k_hlp)
       status = ma_push_get(MT_DBL,n_hlp2,'hlp block',l_hlp2,k_hlp2)
       status = ma_push_get(MT_DBL,n_ssni,'ssni block',l_ssni,k_ssni)
       if (.not.(status)) call errquit
     $      ('moints_2x: cannot allocate local memory',0, MA_ERR)
c
c  Initialize
c
       flop1 = 0.d0
       tint = 0.d0
       t1idx = 0.d0
       t2jidx = 0.d0
       t34jidx = 0.d0
       t2kidx = 0.d0
       t3kidx = 0.d0
       t4kidx = 0.d0
       tottask = 0
       ttaskmax = 0.d0
       ttaskmin = 1.d24
       ttaskagg = 0.d0
       if (oexch) call ga_zero(g_exch)
       if (ocoul) call ga_zero(g_coul)
       ploop = 0
       thalf = util_cpusec()
#ifdef CRAYXT
cedo
c     skip if we do not own any patch of g_exch to avoid slowdowns
c     
       call ga_distribution(g_exch, ga_nodeid(), rlo, rhi, clo, chi )
       doloop=.not.(rlo.eq.0.and.rhi.eq.-1.and.clo.eq.0.and.chi.eq.-1)
       nodesloop=0
       if(doloop) nodesloop=1
       call ga_igop(65536,nodesloop,1,'+')
cedo better use at last 2k procs since they 
cedo we should not have any congestion at 2k
       if(nodesloop.lt.2048) then
          nodesloop=min(2048,ga_nnodes())
          doloop=(ga_nodeid().lt.nodesloop) 
       endif
       if(ga_nodeid().eq.0) then
          write(6,*) ' using ',nodesloop,
     C      ' nodes out of ',ga_nnodes()
          call util_flush(6)
       endif
       if (doloop) then
#else
       next = moints_nxttask(num_nodes)
#endif
c
c  4-fold shell loop
c
       do ii=1,nsh
#ifdef ACC_STRIPS
          if(ga_nodeid().eq.0) then
             write(6,'(i5,a,i6,a,i6,a,f16.2,a)') ga_nodeid(),
     S            ' shell ',ii,
     c            ' out of ',nsh, ' in ',util_wallsec(),' secs'
             call util_flush(6)
          endif
#endif

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
#ifdef CRAYXT
             if (mod(ploop,nodesloop).eq.ga_nodeid()) then
#else
             if (next.eq.ploop) then
#endif
               tottask = tottask + 1
               ttask = util_cpusec()
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
                   tz = util_cpusec()
                   if (lgr.ne.kgr) then
                     call moints_blktr( ilen, jlen, kblen, lblen,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt))
                     call moints_trf1( nbf, qlo, qhi, mo1_lo, mo1_hi, 
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 1.d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt), dbl_mb(k_xmo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )

                   else
                     call moints_trf1( nbf, qlo, qhi, mo1_lo, mo1_hi, 
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 0.5d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbb), dbl_mb(k_xmo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   endif
                   t1idx = t1idx + util_cpusec() - tz
                 enddo
               enddo
               if (ocoul) then
                 tz = util_cpusec()
                 call moints_trf2J( nbf, qlo, qhi, 
     $                              occ_start, mo1_lo, mo1_hi,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              dbl_mb(k_ssni), dbl_mb(k_hlp),
     $                              dbl_mb(k_ssbb), dbl_mb(k_xmo),
     $                              g_coul )
                 t2jidx = t2jidx + util_cpusec() - tz
               endif
               if ((oexch).and.(ohalf)) then
                 tz = util_cpusec()
                 call moints2x_trf2K( nbf, qlo, qhi, 
     $                                occ_start, mo1_lo, mo1_hi,
     $                                ibflo, ibfhi, jbflo, jbfhi,
     $                                int_mb(k_bfmap),
     $                                scale, dbl_mb(k_ssni),
     $                                dbl_mb(k_hlp), dbl_mb(k_hlp),
     $                                dbl_mb(k_hlp2),
     $                                dbl_mb(k_mo), g_exch )
                 t2kidx = t2kidx + util_cpusec() - tz
               endif
               if ((oexch).and.(.not.(ohalf))) then
                 tz = util_cpusec()
                 call moints2xv_trf23K( nbf, qlo, qhi, occ_start,
     $                                  mo1_lo, mo1_hi, mo2_lo, mo2_hi,
     $                                  ibflo, ibfhi, jbflo, jbfhi, 
     $                                  scale, dbl_mb(k_ssni),
     $                                  dbl_mb(k_hlp2),
     $                                  dbl_mb(k_hlp), dbl_mb(k_hlp),
     $                                  dbl_mb(k_mo), dbl_mb(k_xmo),
     $                                  g_exch )
                 t3kidx = t3kidx + util_cpusec() - tz
               endif
#ifndef CRAYXT
               next = moints_nxttask(num_nodes)
#endif
               ttask = util_cpusec() - ttask
               ttaskmax = max(ttaskmax,ttask)
               ttaskmin = min(ttaskmin,ttask)
               ttaskagg = ttaskagg + ttask
             endif
             ploop = ploop + 1
           endif
         enddo
       enddo
#ifdef CRAYXT
       endif!doloop
#endif
       tsynch = util_cpusec()
       call ga_sync()
#ifndef CRAYXT
       next = moints_nxttask(-num_nodes)
#endif
       tsynch = util_cpusec() - tsynch
       thalf = util_cpusec() - thalf
       call moints_aodisk_tidy()
c
*       write(6,*) ' halfway g_exch ', ga_ddot(g_exch, g_exch)
c
c  Clean-up
c
       if (.not. ma_pop_stack(l_ssni))
     $     call errquit('moints: failed to pop', l_ssni, MA_ERR)
       if (.not. ma_pop_stack(l_hlp2))
     $     call errquit('moints: failed to pop', l_hlp2, MA_ERR)
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
       if (.not. ma_pop_stack(l_iscr))
     $     call errquit('moints: failed to pop', l_iscr, MA_ERR)
       if (.not. ma_pop_stack(l_eri))
     $     call errquit('moints: failed to pop', l_eri, MA_ERR)
c
c  2nd half transformation
c
       status = ma_push_get(MT_DBL,(nbf*nbf),'hlp',l_hlp,k_hlp)
c
c  Exchange term
c
       if ((oexch).and.(ohalf)) then
         tz = util_cpusec()
         if (util_print('exchange half integral',print_debug))
     $     call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_exch)
         call moints_Ktrf34( nbf, qlo, qhi, occ_start, mo1_lo, mo1_hi,
     $                       mo2_lo, mo2_hi, .true.,
     $                       dbl_mb(k_mo), dbl_mb(k_hlp),
     $                       osym, int_mb(k_sym), g_exch)
         t4kidx = t4kidx + util_cpusec() - tz
C	 CALL MOINTS_OP_PRINT(OCC_START,MO1_LO,MO1_HI,NMO2,G_EXCH)
         if (util_print('exchange integral',print_debug))
     $      call moints_op_print(occ_start,mo1_lo,mo1_hi,nmo2,g_exch)
       endif
       if ((oexch).and.(.not.(ohalf))) then
         tz = util_cpusec()
         call moints2xv_Ktrf4(nbf, qlo, qhi, occ_start, mo1_lo, mo1_hi,
     $                        mo2_lo, mo2_hi, dbl_mb(k_mo),
     $                        dbl_mb(k_hlp), osym, int_mb(k_sym),
     $                        g_exch )
         t4kidx = t4kidx + util_cpusec() - tz
C	 CALL MOINTS_OP_PRINT(OCC_START,MO1_LO,MO1_HI,NMO2,G_EXCH)
         if (util_print('exchange integral',print_debug))
     $      call moints_op_print(occ_start,mo1_lo,mo1_hi,nmo2,g_exch)
       endif
c
c  Coulomb term
c
       if (ocoul) then
         tz = util_cpusec()
         if (util_print('coulomb half integral',print_debug))
     $      call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_coul)
         call moints_Jtrf34(nbf, qlo, qhi, occ_start, mo1_lo, mo1_hi,
     $                      mo2_lo, mo2_hi, dbl_mb(k_mo),
     $                      dbl_mb(k_hlp), osym, int_mb(k_sym),
     $                      g_coul )
         t34jidx = util_cpusec() - tz
C	 CALL MOINTS_OP_PRINT(OCC_START,MO1_LO,MO1_HI,NMO2,G_COUL)
         if (util_print('coulomb integral',print_debug))
     $      call moints_op_print(occ_start,mo1_lo,mo1_hi,nmo2,g_coul)
       endif
c
       call ga_sync
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
#ifdef NOCOMMS
       if (ga_nodeid().eq.0) write(6,334)
 334   format(/,10x,'**** COMMUNICATION DISABLED ****')
#endif

c
c  Timings and Statistics bookkeeping
c
       flop1 = flop1*1.d-6
       ttotal = util_cpusec() - ttotal
       mi_npass = mi_npass + 1.d0
       mi_ttotal = mi_ttotal + ttotal
       mi_thalf = mi_thalf + thalf
       mi_tint = mi_tint + tint
       mi_t1 = mi_t1 + t1idx
       mi_t2k = mi_t2k + t2kidx
       mi_t3k = mi_t3k + t3kidx
       mi_t4k = mi_t4k + t4kidx
       mi_t2j = mi_t2j + t2jidx
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
c
c
       call ga_sync()
*       write(6,*) ' fullway g_exch ', ga_ddot(g_exch, g_exch)
c
       return
       end









c
c  ------ old version -------
c  Exchange 2nd Index transformation to occupied
c  Used when half-transformed integrals
c  are required.
c
       subroutine moints2x_trf2K_old( nbf, ostart, olo, ohi, ilo, ihi,
     $                           jlo, jhi, scale, ssni, h1, h2,
     $                           c, g_exch )
       implicit none
       integer nbf, ostart, olo, ohi, ilo, ihi, jlo, jhi
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision h1(nbf,ilo:ihi), h2(nbf,jlo:jhi)
       double precision c(nbf,nbf)
       integer g_exch
c
       integer nni,nnj,ab,ijlo,ijhi,ofroz
       integer aoff,a,b,i,j,iz
c
       ofroz = ostart - 1
       aoff = ((olo-ofroz)*(olo-ofroz-1))/2
       nni = (ihi - ilo + 1)*nbf
       nnj = (jhi - jlo + 1)*nbf
       do a=olo,ohi
         do b=ostart,a
           call dfill(nni,0.d0,h1,1)
           do i=ilo,ihi
             do j=jlo,jhi
               do iz=1,nbf
                 h1(iz,i) = h1(iz,i) + scale*c(j,b)*ssni(iz,j,i,a)
               enddo
             enddo
           enddo
#ifndef NOCOMMS
           ab = ((a-ofroz)*(a-ofroz-1))/2 + b - ofroz - aoff
           ijlo = (ilo-1)*nbf + 1
           ijhi = ihi*nbf
           call ga_acc(g_exch,ijlo,ijhi,ab,ab,h1,nni,1.d0)
#endif
         enddo
       enddo
c
       do a=olo,ohi
         do b=ostart,a
           call dfill(nnj,0.d0,h2,1)
           do j=jlo,jhi
             do i=ilo,ihi
               do iz=1,nbf
                 h2(iz,j) = h2(iz,j) + scale*c(i,b)*ssni(iz,j,i,a)
               enddo
             enddo
           enddo
#ifndef NOCOMMS
           ab = ((a-ofroz)*(a-ofroz-1))/2 + b - ofroz - aoff
           ijlo = (jlo-1)*nbf + 1
           ijhi = jhi*nbf
           call ga_acc(g_exch,ijlo,ijhi,ab,ab,h2,nnj,1.d0)
#endif
         enddo
       enddo
       return
       end







c
c  ----- old version -----
c  Exchange 2nd index & 3rd index transformations
c  Note that the 2nd index is transformed to virtuals
c  initially.
c
       subroutine moints2xv_trf23K_old( nbf, ostart, olo, ohi, nvir,
     $                              vlo, vhi, ilo, ihi, jlo, jhi,
     $                              scale, ssni, ssv, h1, h2, c,
     $                              g_exch )
       implicit none
       integer nbf, ostart, ohi, olo, nvir, vlo, vhi
       integer ilo, ihi, jlo, jhi
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision ssv(nvir,jlo:jhi,ilo:ihi)
       double precision h1(nvir,ilo:ihi), h2(nvir,jlo:jhi)
       double precision c(nbf,nbf)
       integer g_exch

c
       integer aoff,nni,nnij,nnj,ab,ivlo,ivhi,jvlo,jvhi,ofroz
       integer ssij
       integer aa,bb,a,b,i,j,v
       double precision cc
c     
       ofroz = ostart - 1
       ssij = (ihi - ilo + 1)*(jhi - jlo + 1)
       nni = (ihi - ilo + 1)*nvir
       nnij = (jhi - jlo + 1)*nni
       ivlo = (ilo-1)*nvir + 1
       ivhi = ihi*nvir
       aoff = ((olo-ofroz)*(olo-ofroz-1))/2
       nnj = (jhi - jlo + 1)*nvir
       jvlo = (jlo-1)*nvir + 1
       jvhi = jhi*nvir

       do a=olo,ohi
         call dgemm('t', 'n', nvir, ssij, nbf, 1.d0,
     $               c(1,vlo), nbf, ssni(1,jlo,ilo,a),
     $               nbf, 0.d0, ssv, nvir )
         do b=ostart,a
           call dfill(nni,0.d0,h1,1)
           do j=jlo,jhi
             cc = c(j,b)
             do i=ilo,ihi
               do v=1,nvir
                 h1(v,i) = h1(v,i) + scale*cc*ssv(v,j,i)
               enddo
             enddo
           enddo
#ifndef NOCOMMS
           aa = a - ofroz
           bb = b - ofroz
           ab = (aa*(aa-1))/2 + bb - aoff
           call ga_acc(g_exch,ivlo,ivhi,ab,ab,h1,nni,1.d0)
#endif
         enddo
         do b=ostart,a
           call dfill(nnj,0.d0,h2,1)
           do i=ilo,ihi
            cc = c(i,b)
            do j=jlo,jhi
              do v=1,nvir
                h2(v,j) = h2(v,j) + scale*cc*ssv(v,j,i)
              enddo
            enddo
          enddo
#ifndef NOCOMMS
          aa = a - ofroz
          bb = b - ofroz
          ab = (aa*(aa-1))/2 + bb - aoff
          call ga_acc(g_exch,jvlo,jvhi,ab,ab,h2,nnj,1.d0)
#endif
         enddo
       enddo
c
c
       return
       end







c
c  Transform last index of K to virtual
c
      subroutine moints2xv_Ktrf4( nbf, qlo, qhi, ostart, olo, ohi, 
     $                            vlo, vhi, c, tmp, osym, syms, 
     $                            g_exch )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nbf, qlo, qhi
      integer ostart, olo, ohi, vlo, vhi
      double precision c(nbf,qlo:qhi)
      double precision tmp(vlo:vhi,vlo:vhi)
      logical osym
      integer syms(nbf)
      integer g_exch
c
      integer nvir, clo, chi, rlo, rhi, my_id
      integer ofroz, a1, a2, aa1, aa2, aa, aoff, k_local, ld
#ifdef BAD_GACCESS
      integer l_local
#endif
c
      ofroz = ostart - 1
      nvir = vhi - vlo + 1
      my_id = ga_nodeid()
      call ga_distribution(g_exch, my_id, rlo, rhi, clo, chi )
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      do a1=olo,ohi
        aa1 = a1 - ofroz
        do a2=ostart,a1
          aa2 = a2 - ofroz
          aa = (aa1*(aa1-1))/2 + aa2 - aoff
          if ((aa.ge.clo).and.(aa.le.chi)) then
#ifdef BAD_GACCESS
      ld=rhi-rlo+1
      if(.not.ma_push_get(MT_DBL,ld,
     $  'scratch buff2x', l_local, k_local)) call
     $  errquit('moints2x: pushget failed',0,0)
            call ga_get(g_exch,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
#else
            call ga_access(g_exch,rlo,rhi,aa,aa,k_local,ld)
#endif
            call dgemm('t', 't', nvir, nvir, nbf, 1.d0, c(1,vlo), nbf,
     $                  dbl_mb(k_local), nvir, 0.d0, tmp, nvir )
            if (osym) call moints_sym_zero( nbf, vlo, vhi, syms,
     $                                      a1, a2, tmp)
            call dcopy((nvir*nvir),tmp,1,dbl_mb(k_local),1)
#ifdef BAD_GACCESS
        call ga_put(g_exch,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('moints2x: popstack failed',0,0)
#else
            call ga_release(g_exch,rlo,rhi,aa,aa)
#endif
          endif
        enddo
      enddo
      return
      end

















c
c  New version using reordered shells
c  ----------------------------------
c
c  2nd Index transformation to occupied for Exchange integrals
c  Used when half-transformed integrals only
c  are required.
c
       subroutine moints2x_trf2K( nbf, qlo, qhi, ostart, olo, ohi, 
     $                            ilo, ihi, jlo, jhi, bfmap, 
     $                            scale, ssni, h1, h2, h3, c, 
     $                            g_exch )
       implicit none
#include "global.fh"
       integer nbf, qlo, qhi
       integer ostart, olo, ohi, ilo, ihi, jlo, jhi
       integer bfmap(nbf)
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision h1(nbf,ilo:ihi)
       double precision h2(nbf,jlo:jhi)
       double precision h3(nbf,*)
       double precision c(nbf,qlo:qhi)
       integer g_exch
c
       integer nni,nnj,ab,ijlo,ijhi,ofroz,ilen,jlen
       integer aoff,a,b,i,j,iz
       integer ain,olb
       double precision xx
c
       ofroz = ostart - 1
       aoff = ((olo-ofroz)*(olo-ofroz-1))/2
       ilen = ihi - ilo + 1
       jlen = jhi - jlo + 1
       nni = ilen*nbf
       nnj = jlen*nbf
cedo       do a=olo,ohi
cstaggering to avoid congestion in ga_acc
       olb=ohi-olo+1
       do ain=ga_nodeid()+1,ga_nodeid()+olb
          a=mod(ain,olb)+olo
         do b=ostart,a
           call dfill(nni,0.d0,h1,1)
           do i=ilo,ihi
             do j=jlo,jhi
               xx = scale*c(j,b)
               do iz=1,nbf
                 h1(iz,i) = h1(iz,i) + xx*ssni(iz,j,i,a)
               enddo
             enddo
           enddo
           call row_exch(nbf, ilen, bfmap, h1, h3 )
c h3(nbf,ilo:ihi)
#ifndef NOCOMMS
           ab = ((a-ofroz)*(a-ofroz-1))/2 + b - ofroz - aoff
           ijlo = (ilo-1)*nbf + 1
           ijhi = ihi*nbf
#ifdef ACC_STRIPS
           call util_gacc_strips(g_exch,ijlo,ijhi,ab,ab,h3,nni,1.d0,
     2          nbf)
#else
           call ga_acc(g_exch,ijlo,ijhi,ab,ab,h3,nni,1.d0)
#endif
#endif
         enddo
       enddo
c
cedo      do a=olo,ohi
cstaggering to avoid congestion in ga_acc
       olb=ohi-olo+1
       do ain=ga_nodeid()+1,ga_nodeid()+olb
          a=mod(ain,olb)+olo
         do b=ostart,a
           call dfill(nnj,0.d0,h2,1)
           do j=jlo,jhi
             do i=ilo,ihi
               xx = scale*c(i,b)
               do iz=1,nbf
                 h2(iz,j) = h2(iz,j) + xx*ssni(iz,j,i,a)
               enddo
             enddo
           enddo
           call row_exch(nbf, jlen, bfmap, h2, h3 )
c h3(1:nbf,jlo:jhi)
#ifndef NOCOMMS
           ab = ((a-ofroz)*(a-ofroz-1))/2 + b - ofroz - aoff
           ijlo = (jlo-1)*nbf + 1
           ijhi = jhi*nbf
#ifdef ACC_STRIPS
           call util_gacc_strips(g_exch,ijlo,ijhi,ab,ab,h3,nnj,1.d0,
     2          nbf          )
#else
           call ga_acc(g_exch,ijlo,ijhi,ab,ab,h3,nnj,1.d0)
#endif
#endif
         enddo
       enddo
       return
       end




c
c  Exchange 2nd index & 3rd index transformations
c  Note that the 2nd index is transformed to virtuals
c  initially.
c
       subroutine moints2xv_trf23K( nbf, qlo, qhi, ostart, olo, ohi, 
     $                              vlo, vhi, ilo, ihi, jlo, jhi, scale, 
     $                              ssni, ssv, h1, h2, c, xc, g_exch )
       implicit none
#include "global.fh"
       integer nbf, qlo, qhi
       integer ostart, ohi, olo, nvir, vlo, vhi
       integer ilo, ihi, jlo, jhi
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision ssv(vlo:vhi,jlo:jhi,ilo:ihi)
       double precision h1(vlo:vhi,ilo:ihi), h2(vlo:vhi,jlo:jhi)
       double precision c(nbf,qlo:qhi), xc(nbf,qlo:qhi)
       integer g_exch

c
       integer aoff,nni,nnij,nnj,ab,ivlo,ivhi,jvlo,jvhi,ofroz
       integer ssij
       integer aa,bb,a,b,i,j,v
       integer ain,olb
       double precision cc
c     
       ofroz = ostart - 1
       nvir = vhi - vlo + 1
       ssij = (ihi - ilo + 1)*(jhi - jlo + 1)
       nni = (ihi - ilo + 1)*nvir
       nnij = (jhi - jlo + 1)*nni
       ivlo = (ilo-1)*nvir + 1
       ivhi = ihi*nvir
       aoff = ((olo-ofroz)*(olo-ofroz-1))/2
       nnj = (jhi - jlo + 1)*nvir
       jvlo = (jlo-1)*nvir + 1
       jvhi = jhi*nvir

cedo       do a=olo,ohi
cstaggering to avoid congestion in ga_acc
       olb=ohi-olo+1
       do ain=ga_nodeid()+1,ga_nodeid()+olb
          a=mod(ain,olb)+olo
         call dgemm('t', 'n', nvir, ssij, nbf, 1.d0,
     $               xc(1,vlo), nbf, ssni(1,jlo,ilo,a),
     $               nbf, 0.d0, ssv, nvir )
         do b=ostart,a
           call dfill(nni,0.d0,h1,1)
           do j=jlo,jhi
             cc = c(j,b)*scale
             if (abs(cc).gt.1.d-16) then
               do i=ilo,ihi
                 do v=vlo,vhi
                   h1(v,i) = h1(v,i) + cc*ssv(v,j,i)
                 enddo
               enddo
             endif
           enddo
#ifndef NOCOMMS
           aa = a - ofroz
           bb = b - ofroz
           ab = (aa*(aa-1))/2 + bb - aoff
           call ga_acc(g_exch,ivlo,ivhi,ab,ab,h1,nni,1.d0)
#endif
         enddo
         do b=ostart,a
           call dfill(nnj,0.d0,h2,1)
           do i=ilo,ihi
             cc = c(i,b)*scale
             if (abs(cc).gt.1.d-16) then
               do j=jlo,jhi
                 do v=vlo,vhi
                   h2(v,j) = h2(v,j) + cc*ssv(v,j,i)
                 enddo
               enddo
             endif
           enddo
#ifndef NOCOMMS
           aa = a - ofroz
           bb = b - ofroz
           ab = (aa*(aa-1))/2 + bb - aoff
           call ga_acc(g_exch,jvlo,jvhi,ab,ab,h2,nnj,1.d0)
#endif
         enddo
       enddo
c
c
       return
       end
