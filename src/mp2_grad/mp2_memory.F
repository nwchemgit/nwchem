      subroutine mp2_memory(basis,geom,rtdb,tol2e,oskel,
     $     nbf,noa,nva,nob,nvb,
     $     batch_size,max_ma,max_ga,max_file,task,scftype)
*     
*     $Id: mp2_memory.F,v 1.22 2003-10-17 22:55:29 carlfahl Exp $
*     
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "util.fh"
c     
c     Return the maximum local memory, global memory per node and 
c     file space per node required by any algorithm step for the 
c     specified task.
c     
c     If (ga_uses_ma) then
c     .   for each step add GA/P to MA and return 0 for GA
c     else
c     .   GA and MA are in separate spaces so return both
c     endif
c     
c     A 2 GB limit is imposed upon on a single file.
c     
      integer nbf,basis,rtdb,geom
      integer noa,nva
      double precision tol2e
      logical oskel
      integer nob,nvb
      integer max_ma
      integer max_ga
      double precision max_file
      integer nnbf, nnbf_sparse
      integer nva_local, nsh, ilo, ihi, jlo, jhi, ish,jsh
      character*(*)task,scftype
      double precision k_file_size,t_file_size,p_file_size
      integer shdim(1),shlo(1),shpairs(1),shpairslocal(1)
      integer map(1024)         ! No. of processors
      integer lengd,lenscrd
      integer nshpairlocal,nshpair,batch_size,shmax
      integer nbfpairlocal,nbfpair
c
      integer ga_1,ga_2,ga_3,ga_4,ga_5,ga_6,ga_7
      integer ma_1,ma_2,ma_3,ma_4,ma_5,ma_6,ma_7,ma_1a,ma_6a
      integer max2e, mem2
      logical oga_uses_ma
c     
      nnbf=nbf*(nbf+1)/2
      nva_local=(nva-1)/ga_nnodes()+1
c     
c     Count no. of non-zero basis function pairs in same way as 
c     moints-semi to estimate the file size
c     
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('mp2_memory: bas_numcont failed', 0, BASIS_ERR)
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      nnbf_sparse = 0
      do ish = 1+ga_nodeid(), nsh, ga_nnodes()
         if (.not. bas_cn2bfr(basis,ish,ilo,ihi)) call errquit
     $        ('mp2_memory: bas_cntobfr failed', ish, BASIS_ERR)
         do jsh = 1, ish
            if (.not. bas_cn2bfr(basis,jsh,jlo,jhi)) call errquit
     $           ('mp2_memory: bas_cntobfr failed', jsh, BASIS_ERR)
            if (schwarz_shell(ish,jsh)*schwarz_max().ge.tol2e) then
               if (ish.ne. jsh) then
                  nnbf_sparse = nnbf_sparse+(ihi-ilo+1)*(jhi-jlo+1)
               else
                  nnbf_sparse = nnbf_sparse+(ihi-ilo+1)*(ihi-ilo+2)/2
               endif
            endif
         enddo
      enddo
      call ga_igop(543, nnbf_sparse, 1, '+')
      if (nnbf_sparse .gt. nnbf) call errquit('mp2mem: confused',0,
     &       MEM_ERR)
c     
      k_file_size=dble(nnbf_sparse)*batch_size*nva_local + nnbf + 2
      t_file_size = 0
      p_file_size = 0
#ifndef EXT_INT
      if (k_file_size .gt. 2.4d8) k_file_size = 1d99
#endif
*      if (task .eq. 'gradient') then
         call mp2_backt_info(basis,tol2e,oskel,
     $        .false.,nshpair,nshpairlocal,
     $        nbfpair,nbfpairlocal,shmax,
     $        shdim,shlo,shpairs,shpairslocal,map)
      if (task .eq. 'gradient') then
         t_file_size=dble((nob*nvb+noa*nva))*batch_size*nva_local
         p_file_size=dble(nbfpairlocal)*nbf*batch_size
#ifndef EXT_INT
         if (t_file_size .gt. 2.4d8) t_file_size = 1d99
         if (p_file_size .gt. 2.4d8) p_file_size = 1d99
#endif
      endif
*      endif
c
      call intb_mem_2e4c(max2e, mem2) ! Determine mem2 = scratch space
      max2e = max(max2e,min(50*shmax**4,21**4)) ! Enuf room for 1 cartesian H shell
c
      call schwarz_tidy()
      call int_terminate()
c
      if (task .eq. 'gradient') then
         call intd_init(rtdb,1,basis)
         call intb_mem_2e4c(lengd, lenscrd) ! blocking algorithm
         lengd = lengd/12
         lengd = max(lengd,15**4)          ! MUST MATCH MP_NONSEP
         call int_terminate
      endif
c
      oga_uses_ma = ga_uses_ma() .and. (.not. ga_memory_limited())
c
      ma_1 = 0
      ma_2 = 0
      ma_3 = 0
      ma_4 = 0
      ma_5 = 0
      ma_6 = 0
      ma_7 = 0
      ga_1 = 0
      ga_2 = 0
      ga_3 = 0
      ga_4 = 0
      ga_5 = 0
      ga_6 = 0
      ga_7 = 0
c     
c     Setup and persistent data
c
      ma_1 = 0
      ga_1 = 2*nbf**2
      if (task .eq. 'gradient') ga_1 = ga_1 + 4*nbf**2
      ga_1 = ga_1 / ga_nnodes()
      if (oga_uses_ma) then
         ma_1 = ma_1 + ga_1
         ga_1 = 0
      endif
c     
c     Moints-semi
c
      ma_2 = ma_1 + (nbf**2 + nsh**2)/ma_sizeof(mt_dbl,1,mt_int) +
     $     nbf**2 + 
     $     shmax**2*batch_size*nbf + 
     $     max(shmax**2*32**2+3*max2e+mem2,2*nva*shmax**2+nbf**2)
      ga_2 = ga_1 + 
     $     min(nbf**2,shmax**2*ga_nnodes())*batch_size*nva / ga_nnodes()
      if (oga_uses_ma) then
         ma_2 = ma_2 + ga_2
         ga_2 = 0
      endif
      if (task .eq. 'energy') goto 1111
c
c     Make-T
c
      ma_1a = ma_1 + nbf**2 + noa*nva + nob*nvb + 
     $     nnbf_sparse/ma_sizeof(mt_dbl,1,mt_int)
      if (scftype .eq. 'UHF') ma_1 = ma_1 + nbf**2
c
      ma_3 = ma_1a + 2*nbf**2 + max(noa*nva,nob*nvb)
      ga_3 = ga_1 + (noa*nva*nva + nob*nvb*nva) / ga_nnodes()
      if (oga_uses_ma) then
         ma_3 = ma_3 + ga_3
         ga_3 = 0
      endif
c
c     Make P and W
c
      ma_4 = ma_1a + 2*nbf**2
      ga_4 = ga_1
      if (oga_uses_ma) then
         ma_4 = ma_4 + ga_4
         ga_4 = 0
      endif
c
c     Make L
c
      ma_5 = ma_1a + nvb**2 + noa**2 + 2*noa*nva + 2*nob*nvb + 
     $     max(noa*nva,nob*nvb) + 2*nbf**2
      ga_5 = ga_1
      if (oga_uses_ma) then
         ma_5 = ma_5 + ga_5
         ga_5 = 0
      endif
c
c     Back transform
c
      ma_1a = ma_1a - 
     $     (noa*nva+nob*nvb + nnbf_sparse/ma_sizeof(mt_dbl,1,mt_int))
      ma_6a= (nshpair*2 + nshpairlocal*2)/ma_sizeof(mt_dbl,1,mt_int)
      ma_6 = ma_1a + ma_6a + 2*max(nbf*nbf, nbf*shmax*shmax) + nbf**2
      ga_6 = ga_1 + nva*nbfpair / ga_nnodes()
      if (oga_uses_ma) then
         ma_6 = ma_6 + ga_6
         ga_6 = 0
      endif
c
c     Non-separable 2-e gradient
c
      ma_7 = ma_1 + ma_6a + shmax*shmax*nbf*batch_size + 
     $     shmax*shmax*nbf + batch_size*nbf + 12*lengd + shmax**4 +
     $     4*lengd/ma_sizeof(mt_dbl,1,mt_int) + lenscrd
*      if (ga_nodeid().eq.0)
*     $write(6,*) ma_1 , ma_6a , shmax*shmax*nbf*batch_size , 
*     $     shmax*shmax*nbf , batch_size*nbf , 12*lengd , shmax**4 ,
*     $     4*lengd/ma_sizeof(mt_dbl,1,mt_int) , lenscrd
      ga_7 = ga_1
      if (oga_uses_ma) then
         ma_7 = ma_7 + ga_7
         ga_7 = 0
      endif
c
 1111 max_ga = max(ga_1,ga_2,ga_3,ga_4,ga_5,ga_6,ga_7) * 1.02d0
      max_ma = max(ma_1,ma_2,ma_3,ma_4,ma_5,ma_6,ma_7) * 1.02d0
      max_file = max(k_file_size+t_file_size,t_file_size+p_file_size)
     $     * 1.02d0
      call ga_igop(542, max_file, 1, 'min')
c
      if (ga_nodeid().eq.0 .and. 
     $     util_print('mpass debug',print_debug)) then
         write(6,1) batch_size,
     $        ma_1, ga_1, ma_2, ga_2, ma_3, ga_3, ma_4, ga_4, 
     $        ma_5, ga_5, ma_6, ga_6, ma_7, ga_7, 
     $        max_ma, max_ga, k_file_size, 
     $        t_file_size, p_file_size, max_file
 1       format(/'  MP2 multi-pass memory management'/
     $        '  --------------------------------'/
     $        '  Batch size ', i3/
     $        '  Step                       MA        GA '/
     $        '  ----------------------  --------  ---------'/
     $        '  Setup                  ',i9,2x,i9/
     $        '  Forward transformation ',i9,2x,i9/
     $        '  Making amplitudes      ',i9,2x,i9/
     $        '  Making densities       ',i9,2x,i9/
     $        '  Making Lagrangian      ',i9,2x,i9/
     $        '  Back transformation    ',i9,2x,i9/
     $        '  Non-separable gradient ',i9,2x,i9/
     $        '  ----------------------  --------  ---------'/
     $        '  Maximum memory + 2%    ',i9,2x,i9/
     $        '  ----------------------  --------  ---------'/
     $        '  K-file                 ',1p,d9.2,0p/
     $        '  T-file                 ',1p,d9.2,0p/
     $        '  P-file                 ',1p,d9.2,0p/
     $        '  ----------------------  --------'/
     $        '  Maximum filesize + 2%  ',1p,d9.2,0p/)
         call util_flush(6)
      endif
c     
      end
