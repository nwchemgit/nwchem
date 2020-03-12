
      integer function moints_trp_mem(basis, nmo1, blen)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bas.fh"
#include "mafdecls.fh"
      integer basis, nmo1, blen
c
c     return memory required for moints_trp with given arguments
c
      integer nbf, nsh, maxbfsh, max2e, mem2
      integer n_ssni, n_integ, n_ssbb, n_hlp, bsize
      logical status
c
      status = bas_numbf(basis,nbf)
      status = status.and.bas_numcont(basis,nsh)
      status = status.and.bas_nbf_cn_max(basis,maxbfsh)
      if (.not.status) call errquit('moints_trp: cannot get basis info',0,
     &       BASIS_ERR)
      bsize = max(blen,maxbfsh)
      call intb_mem_2e4c(max2e, mem2) ! Determine mem2 = scratch space
      max2e = max(max2e,min(50*maxbfsh**4,21**4)) ! Enuf room for 1 cartesian H shell
      n_integ = max2e+mem2+ma_sizeof(mt_int, 4*max2e, mt_dbl)
      n_ssbb = maxbfsh*maxbfsh*bsize*bsize
      n_ssni = maxbfsh*maxbfsh*nbf*nmo1
      n_hlp = max((bsize*maxbfsh*maxbfsh*nmo1),(maxbfsh*nbf))
c
      moints_trp_mem = 
     $     ma_sizeof(mt_int, 3*(nbf+nsh), mt_dbl) +
     $     nbf**2 + n_ssni + n_ssbb + max(n_integ, n_hlp+n_ssbb)
c
      end
      

c  $Id$
c
c
c     No-one knows what this stuff is doing ... are the next comments correct?
c
c  This routine returns the Coulomb and exchange integral
c  operator matrices for the range of MO-indices as mo_indx_hi, mo_indx_lo
c  The g_coul, g_exch global arrays are ordered as
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

       subroutine moints_trp( rtdb, basis, ohalf,
     $                        occ_start, mo1_lo, mo1_hi,
     $                        mo2_lo, mo2_hi,
     $                        g_coul, ocoul,
     $                        g_exch, oexch,
     $                        blen, cmo, ncor, nocc, nvir, ndel, Tol2E)
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "sym.fh"
#include "util.fh"
#include "schwarz.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "ccsdps.fh"
#include "itri.fh"
c
c Arguments
c
       integer rtdb
       integer basis                          ! Basis handle
       logical ohalf                          ! Toggle alternative algorithm for half-ints intermediates
       integer occ_start                      ! Offset from frozen core
       integer mo1_lo, mo1_hi                 ! 1st Pair Index range
       integer mo2_lo, mo2_hi                 ! 2nd Pair Index range
       integer g_coul                         ! Coulomb operator
       integer g_exch                         ! Exchange operator
       logical ocoul,oexch                    ! Type selection
       integer blen                           ! Blocking length
       integer ncor, nocc, nvir, ndel
       double precision cmo(*)
       Double Precision Tol2E                 ! Integral tolerance
c
c Local variables
c
       integer geom, nmo1, nmo2, nbf, nsh, maxbfsh
       integer bsize, ngrp, kgr, lgr
       integer ish0, jsh0, ish, jsh, ilen, jlen
       integer ibflo,ibfhi,jbflo,jbfhi,kbflo,kbfhi,lbflo,lbfhi
       integer kshlo, kshhi, lshlo, lshhi
       integer kblen, lblen
       integer l_ssbb, k_ssbb, l_ssbbt, k_ssbbt, l_xmo, k_xmo
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_hlp2, k_hlp2, l_hlp3, k_hlp3
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_shmap, k_shmap, l_bfmap, k_bfmap, l_rbfmap, k_rbfmap
       integer l_glo, k_glo, l_ghi, k_ghi
       integer n_ssbb, n_ssni, n_hlp, n_hlp2, n_hlp3, n_ijni
       integer mem2, max2e
       integer num_nodes, ploop, next, nmixed, i, iz, jz, kz
       integer l_sym, k_sym
       integer g_tmp                       
       integer nmo
       double precision scale, schw_ij, q2
*       double precision d1, d2
       logical status, osym, odapt,odoit, oblk, st
c
       integer moints_numgr, gr_len, nxtask
       external moints_numgr, gr_len, nxtask
       double precision wall
       integer n_progr,pct_progr,npair
       parameter(n_progr=20)
       logical i_progr(n_progr+1)
c
       double precision half
       data half/0.5d00/
       data oblk/.true./    ! toggle whether to use Texas blocking or shell quartet AO ints
c
      if (occsdps) then
         call pstat_on(ps_trpmos)
      else
         call qenter('trpmos',0)
      endif

c
c  General basis info
c
       num_nodes = ga_nnodes()
       if (.not. bas_geom(basis, geom)) call errquit
     $      ('moints: cannot get geometry', 0, GEOM_ERR)
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not.status) 
     $  call errquit('moints_trp: cannot get basis info',0,BASIS_ERR)
       nmo1 = mo1_hi - mo1_lo + 1
       nmo2 = mo2_hi - mo2_lo + 1
c
       nmo = ncor + nocc + nvir + ndel
c     
c  Symmetry adapt the MOs and renumber irreps to start at zero
c     
      if (.not.rtdb_get(rtdb, 'scf:adapt',MT_LOG, 1, osym)) then
         osym = sym_number_ops(geom) .gt. 0
      endif
c
c Check if point group is abelian, if not turn off symmetry
c
       osym = osym.and.sym_abelian_group(geom)
       if(ga_nodeid().eq.0) write(6,*) ' symmetry use ',osym
       if (osym) then
          if (.not. ma_push_get(MT_INT, nbf, 'movec syms',l_sym, k_sym))
     $         call errquit('moints_trp: no memory for syms?',0, MA_ERR)
*ga:1:0
          if (.not.ga_create(MT_DBL,nbf,nmo,'tmp',nbf,0,g_tmp))
     $      call errquit('moints_trp: cannot ga_create tmp',0, GA_ERR)
          if (ga_nodeid() .eq. 0) call ga_put(g_tmp,1,nbf,1,nmo,cmo,nbf)
	  call ga_sync
          call sym_movecs_adapt( basis, 1d-8, g_tmp,
     $                            int_mb(k_sym), nmixed )
          call ga_get(g_tmp,1,nbf,1,nmo,cmo,nbf)
          if (.not.ga_destroy(g_tmp))
     $      call errquit('cannot destroy tmp',0, GA_ERR)
          if (nmixed .ne. 0)
     $       call errquit('moints_trp: symmetry contamination', nmixed,
     &       UNKNOWN_ERR)
          do i =0, nmo-1
             int_mb(k_sym+i) = int_mb(k_sym+i) - 1
          enddo
          if (util_print('orbital symmetry',print_debug)) then
            write(6,887) 
 887        format('Symmetry of MOs')
            write(6,888) (int_mb(k_sym+i),i=0,nmo-1)
 888        format(16i3)
          endif
       endif
c
c  Integrals allocation
c
*      call int_mem_2e4c(max2e, mem2)
       call intb_mem_2e4c(max2e, mem2) ! Determine mem2 = scratch space
       max2e = max(max2e,min(50*maxbfsh**4,21**4)) ! Enuf room for 1 cartesian H shell
c
c  Reorder shells by descending shell-length
c
       status = ma_push_get(MT_INT,nsh,'shell order map',
     $                      l_shmap, k_shmap)
       status = ma_push_get(MT_INT,nsh,'group lo', l_glo, k_glo )
     $       .and. status
       status = ma_push_get(MT_INT,nsh,'group hi', l_ghi, k_ghi)
     $       .and. status
       status = ma_push_get(MT_INT,nbf,'basis map',
     $                      l_bfmap, k_bfmap)
     $       .and. status
       status = ma_push_get(MT_INT,nbf,'rev basis map',
     $                      l_rbfmap, k_rbfmap)
     $       .and. status

       if (.not. status) call errquit('moints_trp: map memory?',0,
     &       MA_ERR)

       call moints_shorder( basis, nsh, nbf, blen, ngrp,
     $                      int_mb(k_glo), int_mb(k_ghi),
     $                      int_mb(k_shmap),
     $                      int_mb(k_bfmap), int_mb(k_rbfmap) )
c
c  Copy of MO coefficients with reordered rows
c
       status = ma_push_get(MT_DBL,(nbf*nbf),'reorder mos',
     $                      l_xmo, k_xmo)

       call row_exch( nbf, nbf, int_mb(k_rbfmap), cmo,
     $                       dbl_mb(k_xmo) )
c
c  Temporary partially-transformed arrays
c     
       bsize = max(blen,maxbfsh)
       n_ssbb = maxbfsh*maxbfsh*bsize*bsize
*       n_ssbb1 = max((nmo1*nmo1),n_ssbb)
       n_hlp = max((bsize*maxbfsh*maxbfsh*nmo1),(maxbfsh*nbf))
       n_hlp2 = maxbfsh*maxbfsh*nmo2
       n_hlp3 = maxbfsh*maxbfsh
       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
       status = ma_push_get(MT_DBL,n_ssni,'ssni block',l_ssni,k_ssni)
       if (.not.(status)) call errquit('cannot allocate local memory',0,
     &       MA_ERR)
c
c Initialize
c
       if (oexch) call ga_zero(g_exch)
       if (ocoul) call ga_zero(g_coul)

       call ga_sync
       do ish0=1,n_progr+1
          i_progr(ish0)=.true.
       enddo
      wall=-util_wallsec()

       ploop = 0
       next = nxtask(num_nodes, 1)
c
c  4-fold shell loop
c
       do ish0=1,nsh
         do jsh0=1,ish0
           ish = max(int_mb(k_shmap+ish0-1),int_mb(k_shmap+jsh0-1))
           jsh = min(int_mb(k_shmap+ish0-1),int_mb(k_shmap+jsh0-1))
           status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
           status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
           ilen = ibfhi - ibflo + 1
           jlen = jbfhi - jbflo + 1
           schw_ij = schwarz_shell(ish,jsh)
           scale = 1.d0
           if (ish.eq.jsh) scale = scale*half
           odoit = schw_ij*schwarz_max().ge.tol2e
           if (odoit .and. osym) then
              odoit = sym_shell_pair(basis, ish, jsh, q2)
           endif
           if (odoit) then
              if (next.eq.ploop) then
               n_ijni = ilen*jlen*nbf*nmo1
               call dfill(n_ijni,0.d0,dbl_mb(k_ssni),1)

               if (.not. ma_push_get(MT_DBL,n_ssbb,'ssbb block',
     $                               l_ssbb,k_ssbb)) 
     $           call errquit('moints_trp: push ssbb',n_ssbb,MA_ERR)

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
                   
                   if (.not. ma_push_get(MT_DBL, max2e,'moints: buf',
     $                                   l_eri, k_eri))
     $               call errquit('moints_trp: push buf',max2e,MA_ERR)
                   if (.not. ma_push_get(MT_DBL, mem2, 'moints: scr', 
     $                                   l_iscr, k_iscr)) 
     $               call errquit('moints_trp: push scr',mem2,MA_ERR)

                   call moints_gblk( basis, ish, jsh,
     $                               kshlo, kshhi, lshlo, lshhi,
     $                               int_mb(k_shmap),int_mb(k_rbfmap),
     $                               schw_ij, tol2e, osym, oblk,
     $                               max2e, dbl_mb(k_eri),
     $                               mem2, dbl_mb(k_iscr),
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               dbl_mb(k_ssbb) )

                   if (.not. ma_pop_stack(l_iscr))
     $               call errquit('moints_trp: pop iscr',l_iscr,MA_ERR)
                   if (.not. ma_pop_stack(l_eri))
     $               call errquit('moints_trp: pop eri',l_eri,MA_ERR)

                   if (.not. ma_push_get(MT_DBL,n_hlp,'hlp block',
     $                                   l_hlp,k_hlp)) 
     $               call errquit('moints_trp: push hlp',n_hlp,MA_ERR)

                   if (lgr.ne.kgr) then

                      if (.not. ma_push_get(MT_DBL,n_ssbb,'ssbbt block',
     $                                      l_ssbbt,k_ssbbt)) 
     $                  call errquit('moints_trp: push ssbbt',n_ssbb,
     $                               MA_ERR)

                     call moints_blktr( ilen, jlen, kblen, lblen,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt))

                     call moints_trf1( nbf, 1, nbf, mo1_lo, mo1_hi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 1.d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt), dbl_mb(k_xmo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )

                     if (.not. ma_pop_stack(l_ssbbt))
     $                    call errquit('moints_trp: pop ssbbt',l_ssbbt,
     $                                 MA_ERR)

                   else
                     call moints_trf1( nbf, 1, nbf, mo1_lo, mo1_hi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, half,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbb), dbl_mb(k_xmo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   endif

                   if (.not. ma_pop_stack(l_hlp))
     $               call errquit('moints_trp: pop hlp',l_hlp,MA_ERR)

                 enddo !lgr
               enddo ! kgr

               if (.not. ma_pop_stack(l_ssbb))
     $           call errquit('moints_trp: pop ssbb',l_ssbb,MA_ERR)

               if (.not. ma_push_get(MT_DBL,n_hlp,'hlp block',
     $                               l_hlp,k_hlp)) 
     $           call errquit('moints_trp: push hlp',n_hlp,MA_ERR)

               if (oexch) then

                  if (.not. ma_push_get(MT_DBL,n_hlp2,'hlp2 block',
     $                                  l_hlp2,k_hlp2))
     $              call errquit('moints_trp: push hlp2',n_hlp2,MA_ERR)

                 call mo_trp_trf23K( nbf, occ_start,
     $                               mo1_lo, mo1_hi, nmo2,
     $                               mo2_lo, mo2_hi, ibflo, ibfhi,
     $                               jbflo, jbfhi, scale,
     $                               dbl_mb(k_ssni),dbl_mb(k_hlp2),
     $                               dbl_mb(k_hlp), dbl_mb(k_hlp),
     $                               dbl_mb(k_xmo), cmo, 
     $                               g_exch, ncor, nocc)

                 if (.not. ma_pop_stack(l_hlp2))
     $             call errquit('moints_trp: pop hlp2',l_hlp2,MA_ERR)

               endif
               if (ocoul) then

                  if (.not. ma_push_get(MT_DBL,n_hlp3,'hlp3 block',
     $                                  l_hlp3,k_hlp3))
     $              call errquit('moints_trp: push hlp3',n_hlp3,MA_ERR)

                 call mo_trp_trf2J( nbf, occ_start, mo1_lo, mo1_hi,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              dbl_mb(k_ssni),
     $                              dbl_mb(k_hlp3), dbl_mb(k_xmo),
     $                              g_coul, ncor, nocc)

                 if (.not. ma_pop_stack(l_hlp3))
     $             call errquit('moints_trp: pop hlp3',l_hlp3,MA_ERR)

               endif
               if(ga_nodeid().eq.0) then
                  npair=itri(ish0,jsh0)+1
               pct_progr=(npair*n_progr)/(nsh*(nsh+1)/2)+1
               if(i_progr(pct_progr)) then
                  i_progr(pct_progr)=.false.
                  write(6,'(a,i9,a,i9,a,i3,a,f10.1,a)')
     &                  'task ',npair,' out of ',nsh*(nsh+1)/2,
     D                 ' done ',
     &                  int((npair*100d0)/(nsh*(nsh+1)/2)),'%  at',
     &                  wall+util_wallsec(),' secs'
                  call util_flush(6)
               endif
            endif
               if (.not. ma_pop_stack(l_hlp))
     $           call errquit('moints_trp: pop hlp',l_hlp,MA_ERR)
               next = nxtask(num_nodes, 1)
             endif
             ploop = ploop + 1
           endif
         enddo
       enddo
       next = nxtask(-num_nodes, 1)
c
c Clean-up
c
       if (.not. ma_chop_stack(l_shmap))
     $     call errquit('moints_trp: failed to pop', l_shmap, MA_ERR)
c
       status = ma_push_get(MT_DBL,(nbf*nbf),'hlp',l_hlp,k_hlp)
       if (oexch) then
           call mo_trp_Ktrf4( g_exch, occ_start, mo1_lo, mo1_hi,
     $                        mo2_lo, mo2_hi, nbf, cmo,
     $                        dbl_mb(k_hlp), ncor, nocc,
     $                        osym, int_mb(k_sym) )
       endif
       if (ocoul) then
         call mo_trp_Jtrf34( g_coul, occ_start, mo1_lo, mo1_hi,
     $                       mo2_lo, mo2_hi, nbf, cmo,
     $                       dbl_mb(k_hlp), ncor, nocc, nvir,
     $                       osym, int_mb(k_sym) )
       endif

       call ga_sync

*       d1 = ga_ddot(g_coul,g_coul)
*       d2 = ga_ddot(g_exch,g_exch)
*       if (ga_nodeid() .eq. 0) then
*         write(6,*) ' GACOUL1 ', d1
*         write(6,*) ' GAEXCH1 ', d2
*       endif
c
c  Clean up
c
       if (.not. ma_pop_stack(l_hlp))
     $   call errquit('moints_trp: pop hlp', l_hlp, MA_ERR)
       if (osym) then
         if (.not. ma_pop_stack(l_sym))
     $     call errquit('moints_trp: pop sym',0, MA_ERR)
       endif
c
c
c
#ifdef NOCOMMS
       if (ga_nodeid().eq.0) write(6,334)
 334   format(/,10x,'**** COMMUNICATION DISABLED ****')
#endif
c
       if (occsdps) then
          call pstat_off(ps_trpmos)
       else
          call qexit('trpmos',0)
       endif

       return
       end





       subroutine mo_trp_trf23K( nbf, ostart, olo, ohi, nvir,
     $                           vlo, vhi, ilo, ihi, jlo, jhi,
     $                           scale, ssni, ssv, h1, h2, xc, c,
     $                           g_exch, ncor, nocc )
       implicit none
#include "global.fh"
#include "yflop.fh"
       integer nbf, ostart, ohi, olo, nvir, vlo, vhi
       integer ilo, ihi, jlo, jhi
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision ssv(nvir,jlo:jhi,ilo:ihi)
       double precision h1(nvir,ilo:ihi), h2(nvir,jlo:jhi)
       double precision xc(nbf,nbf), c(nbf,nbf)
       integer g_exch,ncor,nocc

c
       integer aoff,nni,nnj,ab,ivlo,ivhi,jvlo,jvhi,ofroz
       integer ssij,ain,olb
       integer a,b,i,j,v
       double precision cc
c
       ofroz = ostart - 1
       ssij = (ihi - ilo + 1)*(jhi - jlo + 1)
       nni = (ihi - ilo + 1)*nvir
       ivlo = (ilo-1)*nvir + 1
       ivhi = ihi*nvir
       aoff = ((olo-ofroz)*(olo-ofroz-1))/2
       nnj = (jhi - jlo + 1)*nvir
       jvlo = (jlo-1)*nvir + 1
       jvhi = jhi*nvir

cedo       do a=olo,ohi
cstaggering to avoid congestion in ga-acc
       olb=ohi-olo+1
       do ain=ga_nodeid()+1,ga_nodeid()+olb
          a=mod(ain,olb)+olo
         call dgemm('t', 'n', nvir, ssij, nbf, 1.d0,
     $               xc(1,vlo), nbf, ssni(1,jlo,ilo,a),
     $               nbf, 0.d0, ssv, nvir )
! scale=1.0 is the common case, scale=0.5 is exception on i=j
         do b=ostart,ostart+nocc-1
           call dfill(nni,0.d0,h1,1)
! this could be DGEMV if ssv(v,j,i) is transposed to ssv2(v,i,j)
! in the current form, we can use a loop over range(i) DAXPYs
#if USE_BLAS_NOT_LOOPS
             do j=jlo,jhi
               do i=ilo,ihi
                  call daxpy(nvir,scale*c(j,b),ssv(1,j,i),1,h1(1,i),1)
               enddo
             enddo
#elif USE_NEW_LOOPS
           this code is not correct
           if (abs(scale-1.0d0).lt.(1.e-13)) then
             do j=jlo,jhi
               do i=ilo,ihi
                 do v=1,nvir
                   h1(v,i) = h1(v,i) + c(j,b)*ssv(v,j,i)
                 enddo
               enddo
             enddo
           else ! scale .ne. 1.0
             do j=jlo,jhi
               do i=ilo,ihi
                 do v=1,nvir
                   h1(v,i) = h1(v,i) + scale*c(j,b)*ssv(v,j,i)
                 enddo
               enddo
             enddo
           endif
#else
#ifdef USE_YFLOP
      flops_ycount = flops_ycount + nnj*(ihi-ilo+1)*3
#endif
           do j=jlo,jhi
             cc = c(j,b)
             do i=ilo,ihi
               do v=1,nvir
                 h1(v,i) = h1(v,i) + scale*cc*ssv(v,j,i)
               enddo
             enddo
           enddo
#endif
#ifndef NOCOMMS
           ab = (a-olo)*nocc+b-ostart+1
#ifdef ACC_STRIPS
           call util_gacc_strips(g_exch,ivlo,ivhi,ab,ab,h1,nni,1.d0,
     X          nvir)
#else           
           call ga_acc(g_exch,ivlo,ivhi,ab,ab,h1,nni,1.d0)
#endif
#endif
         enddo
         do b=ostart,ostart+nocc-1
           call dfill(nnj,0.d0,h2,1)
! this is a DGEMV call of the form
! h2(vj) += cb(i) * ssv(vj,i)
! where we fuse vj=j*vdim+v
! and set cb(:) = c(:,b)
#if USE_BLAS_NOT_LOOPS
! TODO: actually use DGEMV here
          do i=ilo,ihi
            do j=jlo,jhi
              do v=1,nvir
                h2(v,j) = h2(v,j) + scale*c(i,b)*ssv(v,j,i)
              enddo
            enddo
          enddo
#elif USE_NEW_LOOPS
          this code is not correct
          if (abs(scale-1.0d0).lt.(1.e-13)) then
            do i=ilo,ihi
              do j=jlo,jhi
                do v=1,nvir
                  h2(v,j) = h2(v,j) + c(i,b)*ssv(v,j,i)
                enddo
              enddo
            enddo
          else ! scale .ne. 1.0
            do i=ilo,ihi
              do j=jlo,jhi
                do v=1,nvir
                  h2(v,j) = h2(v,j) + scale*c(i,b)*ssv(v,j,i)
                enddo
              enddo
            enddo
          endif
#else
#ifdef USE_YFLOP
      flops_ycount = flops_ycount + nnj*(ihi-ilo+1)*3
#endif
          do i=ilo,ihi
            cc = c(i,b)
            do j=jlo,jhi
              do v=1,nvir
                h2(v,j) = h2(v,j) + scale*cc*ssv(v,j,i)
              enddo
            enddo
          enddo
#endif
#ifndef NOCOMMS
          ab = (a-olo)*nocc+b-ostart+1
#ifdef ACC_STRIPS
          call util_gacc_strips(g_exch,jvlo,jvhi,ab,ab,h2,nnj,1.d0,
     X         nvir)
#else           
       call ga_acc(g_exch,jvlo,jvhi,ab,ab,h2,nnj,1.d0)
#endif
#endif
         enddo
       enddo
c
c
       return
       end
      subroutine mo_trp_Ktrf4( g_exch, ostart, olo, ohi,
     $                         vlo, vhi, nbf, c, tmp ,ncor, nocc,
     $                         osym, syms )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_exch
      integer nbf, ostart, olo, ohi, vlo, vhi,ncor,nocc
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)
      logical osym
      integer syms(*)
      integer nvir, clo, chi, rlo, rhi, my_id
      integer ofroz, a1, a2, aa1, aa, aoff, k_local, ld
#ifdef BAD_GACCESS
      integer l_local
#endif

      ofroz = ostart - 1
      nvir = vhi - vlo + 1
      call ga_sync
      my_id = ga_nodeid()
      call ga_distribution(g_exch, my_id, rlo, rhi, clo, chi )
#ifdef BAD_GACCESS
      ld=rhi-rlo+1
      if(.not.ma_push_get(MT_DBL,ld,
     $     'scratch buff2x', l_local, k_local)) call
     $     errquit('moints2x: pushget failed',0,0)
#endif
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      do a1=olo,ohi
        aa1 = a1 - ofroz
        do a2=ostart,ostart+nocc-1
          aa = (a1-olo)*nocc+a2-ostart+1
          if ((aa.ge.clo).and.(aa.le.chi)) then
#ifdef BAD_GACCESS
            call ga_get(g_exch,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
#else 
            call ga_access(g_exch,rlo,rhi,aa,aa,k_local,ld)
#endif
            call moints2xv_Ktrf4a( vlo, vhi, nbf, c,
     $                            dbl_mb(k_local), tmp )
            if (osym) call moints_sym_zero( nbf, vlo, vhi, syms,
     $                                      a1, a2, tmp)
            call dcopy((nvir*nvir),tmp,1,dbl_mb(k_local),1)
#ifdef BAD_GACCESS
            call ga_put(g_exch,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
            if(.not.ma_pop_stack(l_local)) call 
     $           errquit('moints2x: popstack failed',0,0)
#else 
            call ga_release(g_exch,rlo,rhi,aa,aa)
#endif
          endif
        enddo
      enddo
#ifdef BAD_GACCESS
      if(.not.ma_pop_stack(l_local)) call 
     $     errquit('moints2x: popstack failed',0,0)
#endif
      call ga_sync
      return
      end
      subroutine mo_trp_trf2J( nbf, ostart, olo, ohi, ilo, ihi,
     $     jlo, jhi, ssni, h2, c, g_coul,
     &     ncor, nocc )
      implicit none
#include "global.fh"
#include "yflop.fh"
      integer nbf, ostart, olo, ohi, ilo, ihi, jlo, jhi
      double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
      double precision h1(21)   ! cartesian H shell 
      double precision h2(jlo:jhi,ilo:ihi)
      double precision c(nbf,nbf)
      integer g_coul,ncor,nocc
      integer ain,olb
c     
      double precision tmp
      integer nni, ijlo, ijhi, ilen, jlen, ijlen
      integer ab, a, b, i, j, aoff, ofroz, u
c     
c     
      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      nni = ilen*nbf
      ijlen = ilen*jlen
cedo      do a=olo,ohi
cstaggering
       olb=ohi-olo+1
#ifdef USE_YFLOP
      flops_ycount = flops_ycount + olb*nocc*ijlen*nbf*2
#endif
       do ain=ga_nodeid()+1,ga_nodeid()+olb
          a=mod(ain,olb)+olo
         do b=ostart,ostart+nocc-1
            do i = ilo,ihi
               do j = jlo,jhi
                  tmp = 0.0d0
                  do u = 1,nbf
                     tmp = tmp + ssni(u,j,i,a)*c(u,b)
                  enddo
                  h2(j,i) = tmp
               enddo
            enddo
*            call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,a),
*     $           nbf,c(1,b),nbf,0.d0,h2,ijlen)
#ifndef NOCOMMS
            ab = (a-olo)*nocc+b-ostart+1
            if(jlo.le.ilo)then
               do i=ilo,ihi
*                  do j=jlo,min(i,jhi)
*                     ijlo = i*(i-1)/2 + j
*                     ijhi = i*(i-1)/2 + j
*                     call ga_acc(g_coul,ijlo,ijhi,ab,ab,h2(j,i),1,1.d0)
*                  enddo
                  ijlo = i*(i-1)/2 + jlo
                  ijhi = i*(i-1)/2 + min(i,jhi)
                  if (ijhi .ge. ijlo) call ga_acc(g_coul, ijlo, ijhi,
     $                 ab, ab, h2(jlo,i), 1, 1.d0)
               enddo
            else
               do j=jlo,jhi
*     do i=ilo,min(j,ihi)
*     ijlo = j*(j-1)/2 + i
*     ijhi = j*(j-1)/2 + i
*     call ga_acc(g_coul,ijlo,ijhi,ab,ab,h2(j,i),1,1.d0)
*     enddo
                  ijlo = j*(j-1)/2 + ilo
                  ijhi = j*(j-1)/2 + min(j,ihi)
                  if (ijhi .ge. ijlo)  then
                     do i = ilo,min(j,ihi)
                        h1(i-ilo+1) = h2(j,i)
                     enddo
#ifdef ACC_STRIPS
                     call util_gacc_strips(g_coul, ijlo, ijhi, ab, ab,
     J                    h1,1,1.d0,0)
#else
                     call ga_acc(g_coul, ijlo, ijhi, ab, ab, 
     $                    h1,1,1.d0)
#endif
                  endif
               enddo
            endif
#endif
         enddo
      enddo
      return
      end
      subroutine mo_trp_Jtrf34( g_coul, ostart, olo, ohi, vlo, vhi,
     $                          nbf, c, tmp, ncor, nocc, nvir,
     $                          osym, syms )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_coul
      integer ostart, olo, ohi, vlo, vhi, nbf, ncor, nocc, nvir
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)
      logical osym
      integer syms(*)

      integer clo, chi, rlo, rhi, my_id
      integer ofroz, aoff, a1, a2, aa, k_local, ld
#ifdef BAD_GACCESS
      integer l_local
#endif

      call ga_sync

      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      my_id = ga_nodeid()
      call ga_distribution(g_coul, my_id, rlo, rhi, clo, chi )
#ifdef BAD_GACCESS
      ld=rhi-rlo+1
      if(.not.ma_push_get(MT_DBL,ld,'scratch buff2x', 
     &                     l_local, k_local)) 
     &  call errquit('moints2x: pushget failed',0,MA_ERR)
#endif
      do a1=olo,ohi
        do a2=ostart,ostart+nocc-1
          aa = (a1-olo)*nocc+a2-ostart+1
          if ((aa.ge.clo).and.(aa.le.chi)) then
#ifdef BAD_GACCESS
            call ga_get(g_coul,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
#else
            call ga_access(g_coul,rlo,rhi,aa,aa,k_local,ld)
#endif
            call mo_trp_Jtrf34a( nbf, c, dbl_mb(k_local), tmp, 
     &                           ncor,nocc,nvir)
            if (osym) call moints_sym_zero(nbf, vlo, vhi, syms, a1, a2,
     $                                     dbl_mb(k_local))
#ifdef BAD_GACCESS
            call ga_put(g_coul,rlo,rhi,aa,aa,dbl_mb(k_local),ld)
#else 
            call ga_release(g_coul,rlo,rhi,aa,aa)
#endif
          endif
        enddo
      enddo
#ifdef BAD_GACCESS
      if(.not.ma_pop_stack(l_local)) 
     &       call errquit('moints2x: popstack failed',0,MA_ERR)
#endif
      call ga_sync()
      return
      end
      subroutine mo_trp_Jtrf34a( nbf, c, x, tmp, ncor, nocc, nvir)
      implicit none
      integer nbf,ncor,nocc,nvir
      double precision c(nbf,nbf)
      double precision x(nbf*nbf)
      double precision tmp(nbf,nbf)

      integer i, j, icnt

      icnt=0
      do i=1,nbf
        do j=1,i
          ! icnt = i*(i-1)/2 + j
          icnt=icnt+1
          tmp(j,i) = x(icnt)
          tmp(i,j) = x(icnt)
        enddo
      enddo

      call dgemm('t','n',nvir,nbf,nbf,1.d0,c(1,ncor+nocc+1),nbf,tmp,nbf,
     $           0.d0,x,nvir)
      call dgemm('n','n',nvir,nvir,nbf,1.d0,x,nvir,c(1,ncor+nocc+1),nbf,
     $           0.d0,tmp,nvir)

      call dcopy(nvir*nvir,tmp,1,x,1)

      return
      end
