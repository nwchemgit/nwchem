c
c $Id: shell_fock_build.F,v 1.1 2000-08-08 20:46:06 windus Exp $
c
      subroutine shell_fock_build(geom, basis, nder, nfock, jfac, 
     $     kfac, tol2e, g_dens, g_fock, oactive)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
      integer geom, basis, nder, nfock 
      integer g_dens(nfock), g_fock(nder, nfock)
      double precision jfac(*), kfac(*), tol2e
      logical oactive(*)        ! [input] T if atom is active
c
      integer nsh, nbf, max_at_bf, max_sh_bf, nshblocks
      integer blen
      integer lend, leri, lscr, max2e
      integer l_shmap, l_shglo, l_shghi, l_bfmap, l_rbfmap,
     $     l_bfglo, l_bfghi, l_labels, l_list, l_q4, l_bftoat
      integer k_shmap, k_shglo, k_shghi, k_bfmap, k_rbfmap,
     $     k_bfglo, k_bfghi, k_labels, k_list, k_q4, k_bftoat
      integer l_shbflo, l_shbfhi, l_eri, l_scr
      integer k_shbflo, k_shbfhi, k_eri, k_scr
      integer i, iat
c
      integer 
     $     k_d_ij,  k_d_kl,  k_d_ik,  k_d_jl,  k_d_il,  k_d_jk, 
     $     k_f_ij,  k_f_kl,  k_f_ik,  k_f_jl,  k_f_il,  k_f_jk,
     $     l_d_ij,  l_d_kl,  l_d_ik,  l_d_jl,  l_d_il,  l_d_jk, 
     $     l_f_ij,  l_f_kl,  l_f_ik,  l_f_jl,  l_f_il,  l_f_jk
c
      integer maxblen           ! maximum value for blen = max no. bf in group
      integer maxsh             ! maximum no. of shells in a group
      integer maxq              ! max quartets in a request
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)
c
      if (nfock .ne. 1) call errquit('NFOCK', nfock)
      call dscal(nfock, 0.5d0, kfac, 1)
c      
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('new_fock: could not get nbf',0)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('new_fock: could not get nsh',0)
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('new_fock: could not get max_at_bf',0)
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('new_fock: could not get max_sh_bf',0)
c
      blen = min(nbf,maxblen,6*max_sh_bf) ! d(6)*6=36, 6**4=1296 quartets
c
      lend = blen*blen*nfock
      if (.not. ma_push_get(mt_dbl, lend, 'd_ij', l_d_ij, k_d_ij)) 
     $     call errquit('new_fock: d_ij', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'd_kl', l_d_kl, k_d_kl)) 
     $     call errquit('new_fock: d_kl', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'd_ik', l_d_ik, k_d_ik)) 
     $     call errquit('new_fock: d_ik', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'd_il', l_d_il, k_d_il)) 
     $     call errquit('new_fock: d_il', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jk', l_d_jk, k_d_jk)) 
     $     call errquit('new_fock: d_jk', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jl', l_d_jl, k_d_jl)) 
     $     call errquit('new_fock: d_jl', lend)
c
      lend = lend*nder
      if (.not. ma_push_get(mt_dbl, lend, 'f_ij', l_f_ij, k_f_ij)) 
     $     call errquit('new_fock: f_ij', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'f_kl', l_f_kl, k_f_kl)) 
     $     call errquit('new_fock: f_kl', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'f_ik', l_f_ik, k_f_ik)) 
     $     call errquit('new_fock: f_ik', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'f_il', l_f_il, k_f_il)) 
     $     call errquit('new_fock: f_il', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jk', l_f_jk, k_f_jk)) 
     $     call errquit('new_fock: f_jk', lend)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jl', l_f_jl, k_f_jl)) 
     $     call errquit('new_fock: f_jl', lend)
c
      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     $  call errquit('shell_fock_build:could not allocate shmap',nsh)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $  call errquit('shell_fock_build:could not allocate blo',nsh)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $  call errquit('shell_fock_build:could not allocate bhi',nsh)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $  call errquit('shell_fock_build:could not allocate bflo',nsh)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $  call errquit('shell_fock_build:could not allocate bfhi',nsh)
      if (.not. ma_push_get(mt_int, nbf, 'bfglo', l_bfglo, k_bfglo))
     $  call errquit('shell_fock_build:could not allocate blo',nbf)
      if (.not. ma_push_get(mt_int, nbf, 'bfghi', l_bfghi, k_bfghi))
     $  call errquit('shell_fock_build:could not allocate bhi',nbf)
      if (.not. ma_push_get(mt_int, nbf, 'bfmap', l_bfmap, k_bfmap))
     $  call errquit('shell_fock_build:could not allocate bfmap',nbf)
      if (.not. ma_push_get(mt_int, nbf, 'rbfmap', l_rbfmap, k_rbfmap))
     $  call errquit('shell_fock_build:could not allocate rbfmap',nbf)
      if (.not. ma_push_get(mt_int, nbf, 'bftoat', l_bftoat, k_bftoat))
     $  call errquit('shell_fock_build:could not allocate bftoat',nbf)
c 
      call grad_shorder(basis, nsh, nbf, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap),
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))
c
      call intb_mem_2e4c(max2e, lscr) ! blocking algorithm
      leri = max(max2e,1296*100) ! 100 D quartets
c 
      if (.not. ma_push_get(mt_dbl,leri,'eri',l_eri,k_eri)) 
     $     call errquit('new_fock:could not allocate buffer',leri)
      if (.not. ma_push_get(mt_dbl,lscr,'scr',l_scr,k_scr))
     $     call errquit('new_fock: scratch alloc failed', lscr)
      if (.not. ma_push_get(mt_int,4*leri,'labels',l_labels,k_labels))
     $     call errquit('new_fock: could not allocate labels',4*leri)
c     
      do i = 1, nbf
         if (.not. bas_bf2ce(basis, i, iat)) call errquit('bf2ce',i)
         int_mb(k_bftoat+i-1) = iat
      end do
c
      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit('new_fock:could not allocate list',4*maxq)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('new_fock:could not allocate q4',maxq)
c
      do i = 1, nfock
         call ga_reorder(g_dens(i), .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
      end do
c
      call newfock(
     $     basis, geom, nder, nfock,
     $     g_dens, g_fock,
     $     dbl_mb(k_d_ij), dbl_mb(k_d_kl), dbl_mb(k_d_ik),  
     $     dbl_mb(k_d_jl), dbl_mb(k_d_il), dbl_mb(k_d_jk), 
     $     dbl_mb(k_f_ij), dbl_mb(k_f_kl), dbl_mb(k_f_ik),  
     $     dbl_mb(k_f_jl), dbl_mb(k_f_il), dbl_mb(k_f_jk), 
     $     lscr, dbl_mb(k_scr), leri, dbl_mb(k_eri), int_mb(k_labels), 
     $     maxq, int_mb(k_list), dbl_mb(k_q4), tol2e,
     $     nsh,  nbf, blen,
     $     nshblocks, int_mb(k_shmap), int_mb(k_shglo), int_mb(k_shghi),
     $     int_mb(k_bfglo), int_mb(k_bfghi), int_mb(k_bfmap), 
     $     int_mb(k_rbfmap),
     $     int_mb(k_bftoat), int_mb(k_shbflo), int_mb(k_shbfhi), 
     $     oactive,
     $     jfac, kfac)
c
      do i = 1, nfock
         call ga_reorder(g_dens(i), .true., int_mb(k_bfmap),
     $        .true., int_mb(k_bfmap))
c        call ga_reorder(g_fock(i), .true., int_mb(k_bfmap),
c    $        .true., int_mb(k_bfmap))
      end do
c
      call dscal(nfock, 2.0d0, kfac, 1)
c
c     do i = 1, nfock
c        call ga_dscal(g_fock(i), 4.0d0)
c        call ga_symmetrize(g_fock(i))
c     end do
      call ga_dscal(g_fock, 4.0d0)
c        call ga_symmetrize(g_fock(i))
c
c  Clean up memory allocated in this routine
c
      if (.not. ma_chop_stack(l_d_ij))
     $  call errquit('shell_fock_build:failed chopping MA stack',555)
c
      end
c
      subroutine newfock(
     $     basis, geom, nder, nfock,
     $     g_dens, g_fock,
     $     d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $     f_ij,  f_kl,  f_ik,  f_jl,  f_il,  f_jk, 
     $     lscr, scr, leri, eri, labels, 
     $     maxq, list, q4, tol2e,
     $     nsh,  nbf, blen,
     $     nshblocks, shmap, shglo, shghi, bfglo, bfghi, bfmap, rbfmap,
     $     bftoat, shbflo, shbfhi, oactive, jfac, kfac)
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"
c
      integer basis, geom       ! [input] familiar handles
      integer nder              ! [input] No. of derivatives (1 = energy)
      integer nfock             ! [input] No. of fock matrices
      integer g_dens, g_fock
      integer blen              ! [input] max dimension of density block
      double precision          ! [scratch] 6 blocks per possible density
     $     d_ij(nfock,blen,blen),  d_kl(nfock,blen,blen),  
     $     d_ik(nfock,blen,blen),  d_jl(nfock,blen,blen),  
     $     d_il(nfock,blen,blen),  d_jk(nfock,blen,blen)
      double precision          ! [scratch] 6 blocks per possible fock
     $     f_ij(nder,nfock,blen,blen),  f_kl(nder,nfock,blen,blen),  
     $     f_ik(nder,nfock,blen,blen),  f_jl(nder,nfock,blen,blen),  
     $     f_il(nder,nfock,blen,blen),  f_jk(nder,nfock,blen,blen)
c     
      integer lscr              ! [input] Size of scratch for integrals
      double precision scr(lscr) ![scratch] Scratch space for integrals
      integer leri              ! [input] Size of eri buffer
      double precision eri(*) ! [scratch] Derivative integral buffer
      integer labels(leri,4)    ! [scratch] Labels for integrals
      integer maxq              ! [input] Max. no. of shell quartets
      integer list(maxq,4)      ! [scratch] Buffer for shell labels
      double precision q4(maxq) ! [scratch] Buffer for quartet prefactors
      double precision tol2e    ! [input] Integral screening threshold
      integer nsh               ! [input] Number of shells
      integer nbf               ! [input] No. of basis functions (MCSCF)
      integer nshblocks         ! [input] No. of shell blocks
      integer shmap(nsh)        ! [input] Map from new to old shell order
      integer shglo(nshblocks)  ! [input] First new shell in block
      integer shghi(nshblocks)  ! [input] Last new shell in block
      integer bfglo(nshblocks)  ! [input] First new basis func in block
      integer bfghi(nshblocks)  ! [input] Last new basis func in block
      integer bfmap(nbf)        ! [input] Map from new to old basis func order
      integer rbfmap(nbf)       ! [input] Map from old to new basis func order
      integer bftoat(nbf)       ! [input] Map from old basis func to atom
      integer shbflo(nsh)       ! [input] First new bf in new shell
      integer shbfhi(nsh)       ! [input] Last new bf in new shell
      logical oactive(*)        ! [input] T if atom is active
      double precision jfac(*), kfac(*) ! [input] DFT scaling for Coulomb/Exchange
c
      integer  next, nint, ijklblock
      double precision scale
c
      integer nxtask, task_size
      external nxtask
c
      integer iilo, jjlo, kklo, lllo
      integer iihi, jjhi, kkhi, llhi
      integer ish, jsh, ksh, lsh, idim, jdim, kdim, ldim
      integer ibflo, ibfhi, jbflo, jbfhi, kbflo, kbfhi, lbflo, lbfhi
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer lblockhi, jshtop, kshtop, lshtop
      integer ishblock, jshblock, kshblock, lshblock
      integer iish, jjsh, kksh, llsh
      integer nq, atoms(4)
      double precision smax, sij, sijkl, block_eff
c      
      logical oij, okl, oikjl, omore
      logical intb_init4c, intb_2e4c
      logical intbd_init4c, intbd_2e4c
      external intbd_init4c, intbd_2e4c
      external intb_init4c, intb_2e4c
c
      integer nproc
c     
      smax  = schwarz_max()
      call int_acc_std()
c
      nproc = ga_nnodes()
      task_size = nshblocks*(nshblocks+1)/2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))
      write(6,*) 'nshblocks and task_size ',nshblocks,task_size
c     
c     Parallel loop over quartets of reordered shell blocks
c     
      ijklblock = 0
      next = nxtask(nproc,task_size)
      do ishblock = nshblocks, 1, -1
         ishlo = shglo(ishblock)
         ishhi = shghi(ishblock)
         ibflo = bfglo(ishblock)
         ibfhi = bfghi(ishblock)
         idim  = ibfhi - ibflo + 1
         do jshblock = 1, ishblock
            jshlo = shglo(jshblock)
            jshhi = shghi(jshblock)
            jbflo = bfglo(jshblock)
            jbfhi = bfghi(jshblock)
            jdim  = jbfhi - jbflo + 1
            do kshblock = ishblock, 1, -1
               kshlo = shglo(kshblock)
               kshhi = shghi(kshblock)
               kbflo = bfglo(kshblock)
               kbfhi = bfghi(kshblock)
               kdim  = kbfhi - kbflo + 1
               lblockhi = kshblock
               if (ishblock .eq. kshblock) lblockhi = jshblock
               do lshblock = 1, lblockhi
                  lshlo = shglo(lshblock)
                  lshhi = shghi(lshblock)
                  lbflo = bfglo(lshblock)
                  lbfhi = bfghi(lshblock)
                  ldim  = lbfhi - lbflo + 1
                  if (next .eq. ijklblock) then
c     
c     Get blocks of the one-particle densities
c     
                     call new_get_dens_block(g_dens, nfock, blen, 
     $                    ibflo, ibfhi, jbflo, jbfhi, d_ij, jfac)
                     call new_get_dens_block(g_dens, nfock, blen, 
     $                    kbflo, kbfhi, lbflo, lbfhi, d_kl, jfac)
                     call new_get_dens_block(g_dens, nfock, blen, 
     $                    ibflo, ibfhi, kbflo, kbfhi, d_ik, kfac)
                     call new_get_dens_block(g_dens, nfock, blen, 
     $                    jbflo, jbfhi, lbflo, lbfhi, d_jl, kfac)
                     call new_get_dens_block(g_dens, nfock, blen, 
     $                    ibflo, ibfhi, lbflo, lbfhi, d_il, kfac)
                     call new_get_dens_block(g_dens, nfock, blen, 
     $                    jbflo, jbfhi, kbflo, kbfhi, d_jk, kfac)
                     if (.not. ma_verify_allocator_stuff())
     $                    call errquit(' ',0)
c
                     call dfill(nder*nfock*blen*blen, 0d0, f_ij, 1)
                     call dfill(nder*nfock*blen*blen, 0d0, f_kl, 1)
                     call dfill(nder*nfock*blen*blen, 0d0, f_ik, 1)
                     call dfill(nder*nfock*blen*blen, 0d0, f_jl, 1)
                     call dfill(nder*nfock*blen*blen, 0d0, f_il, 1)
                     call dfill(nder*nfock*blen*blen, 0d0, f_jk, 1)
                     if (.not. ma_verify_allocator_stuff())
     $                    call errquit(' ',0)
c     
c     Build the list of integral shell quartets in NWChem shell labelling
c     
                     oij = ishblock.eq.jshblock
                     okl = kshblock.eq.lshblock
                     oikjl = ishblock.eq.kshblock .and. 
     $                    jshblock.eq.lshblock
c     
                     nq = 0
c     
                     do iish = ishlo, ishhi
                        atoms(1)=bftoat(bfmap(shbflo(iish)))
                        iilo = shbflo(iish)-ibflo+1
                        iihi = shbfhi(iish)-ibflo+1
                        ish = shmap(iish)
                        jshtop = jshhi
                        if (oij) jshtop = iish
                        do jjsh = jshlo, jshtop
                           atoms(2)=bftoat(bfmap(shbflo(jjsh)))
                           jjlo = shbflo(jjsh)-jbflo+1
                           jjhi = shbfhi(jjsh)-jbflo+1
                           jsh = shmap(jjsh)
                           sij = schwarz_shell(ish,jsh)
                           if (sij*smax*8d0 .lt. tol2e) goto 100 ! Next jjsh
                           kshtop = kshhi
                           if (oikjl) kshtop = iish
                           do kksh = kshlo, kshtop
                              atoms(3)=bftoat(bfmap(shbflo(kksh)))
                              kklo = shbflo(kksh)-kbflo+1
                              kkhi = shbfhi(kksh)-kbflo+1
                              ksh = shmap(kksh)
                              lshtop = lshhi
                              if (okl) lshtop = kksh
                              if (oikjl .and. iish.eq.kksh) 
     $                             lshtop = jjsh
                              do llsh = lshlo, lshtop
                                 atoms(4)=bftoat(bfmap(shbflo(llsh)))
                                 lllo = shbflo(llsh)-lbflo+1
                                 llhi = shbfhi(llsh)-lbflo+1
                                 lsh = shmap(llsh)
c     
c     Check on sparsity, symmetry, active list, single center
c     
                                 sijkl = sij*schwarz_shell(ksh,lsh)
c     
                                 if (sijkl*8d0 .lt. tol2e) goto 200 ! next llsh
                                 if (nder.gt.1) then
                                   if (.not. (
     $                                oactive(atoms(1)) .or.
     $                                oactive(atoms(2)) .or.
     $                                oactive(atoms(3)) .or.
     $                                oactive(atoms(4)))) goto 200
                                   if (
     $                                atoms(1).eq.atoms(2) .and.
     $                                atoms(2).eq.atoms(3) .and.
     $                                atoms(3).eq.atoms(4)) goto 200
                                 endif
c     
c     Scale according to permutation symmetry of the shell labels
c     
                                 scale = 1d0
                                 if (ish.eq.jsh) scale = scale*0.5d0
                                 if (ksh.eq.lsh) scale = scale*0.5d0
                                 if (ish.eq.ksh .and. jsh.eq.lsh)
     $                                scale = scale*0.5d0
c     
                                 nq = nq + 1
                                 q4(nq) = scale
                                 list(nq,1) = ish
                                 list(nq,2) = jsh
                                 list(nq,3) = ksh
                                 list(nq,4) = lsh
c                    if (.not. ma_verify_allocator_stuff())
c    $                    call errquit(' ',0)
c     
 200                          end do ! next llsh
                           end do
 100                    end do  ! next jjsh
                     end do
c
                     if (nder.eq.1) then
                       if (.not. intb_init4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri, 
     $                    block_eff)) call errquit
     $                    ('twodd_coul_ex:txs init?',nq)
                     else
                       if (.not. intbd_init4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri,
     $                    block_eff)) call errquit
     $                    ('twodd_coul_ex:txs init?',nq)                        
                     endif
c                    if (.not. ma_verify_allocator_stuff())
c    $                    call errquit(' ',0)
c     
 10                  continue
                     if (nder.eq.1) then
                       omore = intb_2e4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., tol2e, .false.,
     $                    labels(1,1),labels(1,2), 
     $                    labels(1,3), labels(1,4), 
     $                    eri, leri, nint, lscr, scr)
                     else
                       omore = intbd_2e4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., tol2e, .false.,
     $                    labels(1,1),labels(1,2), 
     $                    labels(1,3), labels(1,4), 
     $                    eri, leri, nint, lscr, scr)
                     endif
c                    if (.not. ma_verify_allocator_stuff())
c    $                    call errquit(' ',0)
c     
                     if (nint .gt. 0) then
                       if (nder.eq.1) then
                          call new_fock_doit(nint,labels,leri,eri,tol2e,
     $                    d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $                    nder, nfock, blen, rbfmap,
     $                    ibflo, jbflo, kbflo, lbflo)
                        else
                          call new_fock_doit_der(nint,labels,leri,
     $                    eri,tol2e,
     $                    d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $                    nder, nfock, blen, rbfmap, bftoat,
     $                    ibflo, jbflo, kbflo, lbflo)
                        endif
                     endif

                     if (.not. ma_verify_allocator_stuff())
     $                    call errquit(' ',0)

                     if (omore) goto 10
c                                      
                     call new_acc_fock_block(g_fock, nder, nfock, 
     $                 blen, ibflo, ibfhi, jbflo, jbfhi, f_ij)
                     call new_acc_fock_block(g_fock, nder, nfock, 
     $                 blen, kbflo, kbfhi, lbflo, lbfhi, f_kl)
                     call new_acc_fock_block(g_fock, nder, nfock, 
     $                 blen, ibflo, ibfhi, kbflo, kbfhi, f_ik)
                     call new_acc_fock_block(g_fock, nder, nfock, 
     $                 blen, jbflo, jbfhi, lbflo, lbfhi, f_jl)
                     call new_acc_fock_block(g_fock, nder, nfock, 
     $                 blen, ibflo, ibfhi, lbflo, lbfhi, f_il)
                     call new_acc_fock_block(g_fock, nder, nfock, 
     $                 blen, jbflo, jbfhi, kbflo, kbfhi, f_jk)
                     if (.not. ma_verify_allocator_stuff())
     $                    call errquit(' ',0)
c                                         
                     next = nxtask(nproc,task_size)
                  end if
c     
                  ijklblock = ijklblock + 1
c     
               end do
            end do
         end do
      end do
c
      call int_acc_std()
c
      next = nxtask(-nproc,task_size)
      call ga_sync()
c     
*      write(6,*) ' numq ', numq
c     
      return
      end

      subroutine new_get_dens_block(
     $     g_dens, nfock, blen, ibflo, ibfhi, jbflo, jbfhi, buf, fac)
      implicit none
      integer  nfock, g_dens(nfock), blen, ibflo, ibfhi, jbflo, jbfhi
      integer blen2
      double precision buf(nfock), fac(nfock)
c
c     integer lo(3), hi(3), ld(3)
      integer ifock
c
c     lo(1) = 1
c     hi(1) = nfock
c     ld(1) = nfock
c     lo(2) = ibflo
c     hi(2) = ibfhi
c     ld(2) = blen
c     lo(3) = jbflo
c     hi(3) = jbfhi
c     ld(3) = blen
      blen2 = blen*blen
      call dfill(nfock*blen*blen,0d0,buf,1)
***      call nga_get(g_dens, lo, hi, buf, ld) 
      do ifock = 1, nfock
         call ga_get(g_dens(ifock),ibflo,ibfhi,jbflo,jbfhi,
     *     buf((ifock-1)*blen2+1),blen)
c        call dscal(blen2, fac(ifock), buf((ifock-1)*blen2+1), nfock)
         call dscal(blen2, fac(ifock), buf((ifock-1)*blen2+1), 1)
      end do
c
      end
      subroutine new_acc_fock_block(
     $   g_fock, nder, nfock, blen, ibflo, ibfhi, jbflo, jbfhi, buf)
      implicit none
      integer  g_fock, nder, nfock, blen, ibflo, ibfhi, jbflo, jbfhi
      double precision buf(*)
c
      integer lo(4), hi(4), ld(3)
c
      lo(1) = 1
      hi(1) = nder
      ld(1) = nder
      lo(2) = 1
      hi(2) = nfock
      ld(2) = nfock
      lo(3) = ibflo
      hi(3) = ibfhi
      ld(3) = blen
      lo(4) = jbflo
      hi(4) = jbfhi
      call nga_acc(g_fock, lo, hi, buf, ld, 1d0) 
c     call ga_acc(g_fock,ibflo,ibfhi,jbflo,jbfhi,buf,blen,1d0)
c
      end
      subroutine new_fock_doit(nint,labels,leri,eri,tol2e,
     $     d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $     nder, nfock, blen, rbfmap, ibflo, jbflo, kbflo, lbflo)
      implicit none
c
      integer nder, nfock, blen
      integer nint
      integer leri
      integer labels(leri,4)
      double precision eri(*)
      double precision tol2e
      double precision
     $     d_ij(nfock,blen,blen), d_kl(nfock,blen,blen), 
     $     d_ik(nfock,blen,blen), d_jl(nfock,blen,blen), 
     $     d_il(nfock,blen,blen), d_jk(nfock,blen,blen)
      double precision
     $     f_ij(1,nfock,blen,blen), f_kl(1,nfock,blen,blen), 
     $     f_ik(1,nfock,blen,blen), f_jl(1,nfock,blen,blen), 
     $     f_il(1,nfock,blen,blen), f_jk(1,nfock,blen,blen)
      integer rbfmap(*)
      integer ibflo, jbflo, kbflo, lbflo
c
      integer integ, i, j, k, l, v
      double precision g
c
      do integ = 1, nint
         g = eri(integ)
         if (abs(g) .gt. tol2e) then
            i = labels(integ,1)
            j = labels(integ,2)
            k = labels(integ,3)
            l = labels(integ,4)
            i = rbfmap(i) - ibflo + 1
            j = rbfmap(j) - jbflo + 1
            k = rbfmap(k) - kbflo + 1
            l = rbfmap(l) - lbflo + 1
            do v = 1, nfock
               f_ij(1,v,i,j) = f_ij(1,v,i,j) + g*d_kl(v,k,l)
               f_kl(1,v,k,l) = f_kl(1,v,k,l) + g*d_ij(v,i,j)
               f_ik(1,v,i,k) = f_ik(1,v,i,k) + g*d_jl(v,j,l)
               f_il(1,v,i,l) = f_il(1,v,i,l) + g*d_jk(v,j,k)
               f_jl(1,v,j,l) = f_jl(1,v,j,l) + g*d_ik(v,i,k)
               f_jk(1,v,j,k) = f_jk(1,v,j,k) + g*d_il(v,i,l)
            end do
         end if
      end do
c
      end
c
      subroutine new_fock_doit_der(nint,labels,leri,eri,tol2e,
     $     d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $     nder, nfock, blen, rbfmap, bftoat, 
     $     ibflo, jbflo, kbflo, lbflo)
      implicit none
c
      integer nder, nfock, blen
      integer nint
      integer leri
      integer labels(leri,4)
      double precision eri(*)
      double precision tol2e
      double precision
     $     d_ij(nfock,blen,blen), d_kl(nfock,blen,blen), 
     $     d_ik(nfock,blen,blen), d_jl(nfock,blen,blen), 
     $     d_il(nfock,blen,blen), d_jk(nfock,blen,blen)
      double precision
     $     f_ij(nder,nfock,blen,blen), f_kl(nder,nfock,blen,blen), 
     $     f_ik(nder,nfock,blen,blen), f_jl(nder,nfock,blen,blen), 
     $     f_il(nder,nfock,blen,blen), f_jk(nder,nfock,blen,blen)
      integer rbfmap(*)
      integer bftoat(*)       ! [input] Map from old basis func to atom
      integer ibflo, jbflo, kbflo, lbflo
c
      integer integ, i, j, k, l, v
      integer nbuf, iiat, iat
      integer atoms(4),nfockder,xyz
      double precision g
c
      nbuf = 1
      do integ = 1, nint
c        if (abs(g) .gt. tol2e) then
            i = labels(integ,1)
            j = labels(integ,2)
            k = labels(integ,3)
            l = labels(integ,4)
            atoms(1) = bftoat(i)
            atoms(2) = bftoat(j)
            atoms(3) = bftoat(k)
            atoms(4) = bftoat(l)
            i = rbfmap(i) - ibflo + 1
            j = rbfmap(j) - jbflo + 1
            k = rbfmap(k) - kbflo + 1
            l = rbfmap(l) - lbflo + 1
            do iiat = 1, 4
              iat = atoms(iiat)
              nfockder = (iat-1)*3 + 1
              do xyz = nfockder, nfockder+2
              g = eri(nbuf)
              do v = 1, nfock
                f_ij(xyz,v,i,j) = f_ij(xyz,v,i,j) + g*d_kl(v,k,l)
                f_kl(xyz,v,k,l) = f_kl(xyz,v,k,l) + g*d_ij(v,i,j)
                f_ik(xyz,v,i,k) = f_ik(xyz,v,i,k) + g*d_jl(v,j,l)
                f_il(xyz,v,i,l) = f_il(xyz,v,i,l) + g*d_jk(v,j,k)
                f_jl(xyz,v,j,l) = f_jl(xyz,v,j,l) + g*d_ik(v,i,k)
                f_jk(xyz,v,j,k) = f_jk(xyz,v,j,k) + g*d_il(v,i,l)
              end do
              nbuf = nbuf + 1
              end do
            end do
c        end if
      end do
c
      end
