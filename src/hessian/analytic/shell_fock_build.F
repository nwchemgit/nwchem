c
c $Id: shell_fock_build.F,v 1.22 2002-01-16 17:21:38 windus Exp $
c
      subroutine shell_fock_build(geom,basis, nder, ndens,
     $     jfac, kfac, tol2e, g_dens, g_fock, osym)
c
C Took oactive out for the time being, but will probably need to
c bring it back - tlw
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
      integer basis, nder, ndens, nfock
      integer geom
      integer g_dens, g_fock
      double precision jfac(*), kfac(*), tol2e
c      logical oactive(*)        ! [input] T if atom is active
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
      integer i, j, iat
      integer ityp, ndim, dims(3)
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
      logical osym
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)
      integer ga_create_atom_blocked
      logical xc_gotxc
      external ga_create_atom_blocked,xc_gotxc
c
c Quietly modify the values of jfac and kfac so that the original
c values coming in can be consistent with the original fock_2e
c
      do i= 1, ndens
         jfac(i) = jfac(i)*4.0d00
         kfac(i) = kfac(i)*2.0d00
      enddo
c
c Get the dimensions of g_fock and eventually g_dens.
c
      call nga_inquire(g_fock, ityp, ndim, dims)
      nfock = dims(1)
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
      lend = blen*blen*ndens
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
      lend = blen*blen*nfock
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
      call nga_reorder(g_dens, .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
      call nga_reorder(g_fock, .true., int_mb(k_rbfmap),
     $        .true., int_mb(k_rbfmap))
c
      call newfock(
     $     basis, nder, ndens, nfock,
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
c     $     oactive,
     $     jfac, kfac)
c
      call nga_reorder(g_dens, .true., int_mb(k_bfmap),
     $     .true., int_mb(k_bfmap))
      call nga_reorder(g_fock, .true., int_mb(k_bfmap),
     $     .true., int_mb(k_bfmap))
c
c  Clean up memory allocated in this routine
c
      if (.not. ma_chop_stack(l_d_ij))
     $  call errquit('shell_fock_build:failed chopping MA stack',555)
c
c Quietly modify the values of jfac and kfac back to the original
c
      do i= 1, ndens
         jfac(i) = jfac(i)/4.0d00
         kfac(i) = kfac(i)/2.0d00
      enddo
c
      if(xc_gotxc().and.nder.eq.0)then
         call xc_newfock(geom,basis,jfac,kfac,nbf,
     %        nder, ndens, nfock, g_fock,g_dens)
      endif
c
      if (osym) call nga_symmetrize(g_fock)
c
      return
      end
c
      subroutine newfock(
     $     basis, nder, ndens, nfock,
     $     g_dens, g_fock,
     $     d_ij,  d_kl,  d_ik,  d_jl,  d_il,  d_jk, 
     $     f_ij,  f_kl,  f_ik,  f_jl,  f_il,  f_jk, 
     $     lscr, scr, leri, eri, labels, 
     $     maxq, list, q4, tol2e,
     $     nsh,  nbf, blen,
     $     nshblocks, shmap, shglo, shghi, bfglo, bfghi, bfmap, rbfmap,
     $     bftoat, shbflo, shbfhi, 
c     $     oactive, 
     $     jfac, kfac)
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"
c
      integer basis             ! [input] familiar handles
      integer nder              ! [input] No. of derivatives (0 = energy)
      integer ndens             ! [input] No. of density matrices (for UHF)
      integer nfock             ! [input] No. of fock matrices
      integer g_dens, g_fock
      integer blen              ! [input] max dimension of density block
      double precision          ! [scratch] 6 blocks per possible density
     $     d_ij,  d_kl,  
     $     d_ik,  d_jl,  
     $     d_il,  d_jk
      double precision          ! [scratch] 6 blocks per possible fock
     $     f_ij,  f_kl,  
     $     f_ik,  f_jl,  
     $     f_il,  f_jk
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
c      logical oactive(*)        ! [input] T if atom is active
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
      integer id, jd
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
c     write(6,*) 'nshblocks and task_size ',nshblocks,task_size
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
                     call new_get_dens_block(g_dens, ndens, blen, 
     $                    ibflo, ibfhi, jbflo, jbfhi, d_ij, jfac)
c                    call print_dens_bit(d_ij, ndens, ibflo, ibfhi,
c    $                    jbflo, jbfhi, blen)
                     call new_get_dens_block(g_dens, ndens, blen, 
     $                    kbflo, kbfhi, lbflo, lbfhi, d_kl, jfac)
                     call new_get_dens_block(g_dens, ndens, blen, 
     $                    ibflo, ibfhi, kbflo, kbfhi, d_ik, kfac)
                     call new_get_dens_block(g_dens, ndens, blen, 
     $                    jbflo, jbfhi, lbflo, lbfhi, d_jl, kfac)
                     call new_get_dens_block(g_dens, ndens, blen, 
     $                    ibflo, ibfhi, lbflo, lbfhi, d_il, kfac)
                     call new_get_dens_block(g_dens, ndens, blen, 
     $                    jbflo, jbfhi, kbflo, kbfhi, d_jk, kfac)
                     if (.not. ma_verify_allocator_stuff())
     $                    call errquit(' ',0)
c
                     call dfill(nfock*blen*blen, 0d0, f_ij, 1)
                     call dfill(nfock*blen*blen, 0d0, f_kl, 1)
                     call dfill(nfock*blen*blen, 0d0, f_ik, 1)
                     call dfill(nfock*blen*blen, 0d0, f_jl, 1)
                     call dfill(nfock*blen*blen, 0d0, f_il, 1)
                     call dfill(nfock*blen*blen, 0d0, f_jk, 1)
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
                                 if (nder.gt.0) then
c                                   if (.not. (
c     $                                oactive(atoms(1)) .or.
c     $                                oactive(atoms(2)) .or.
c     $                                oactive(atoms(3)) .or.
c     $                                oactive(atoms(4)))) goto 200
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
                     if (nder.eq.0) then
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
                     if (nder.eq.0) then
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
c
c     write(6,*) ' nint ', nint, omore
c
                     endif
c                    if (.not. ma_verify_allocator_stuff())
c    $                    call errquit(' ',0)
c     
                     if (nint .gt. 0) then
                       if (nder.eq.0) then
                          call new_fock_doit(nint,labels,leri,eri,tol2e,
     $                    d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $                    ndens, nfock, blen, rbfmap,
     $                    ibflo, jbflo, kbflo, lbflo)
                        else
                          call new_fock_doit_der(nint,labels,leri,
     $                    eri,tol2e,
     $                    d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $                    ndens, nfock/ndens, blen, rbfmap, bftoat,
     $                    ibflo, jbflo, kbflo, lbflo)
c                                      
                        endif
                     endif

                     if (.not. ma_verify_allocator_stuff())
     $                    call errquit(' ',0)

                     if (omore) goto 10
c                                      
                     call new_acc_fock_block(g_fock, nfock, 
     $                 blen, ibflo, ibfhi, jbflo, jbfhi, f_ij)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 blen, kbflo, kbfhi, lbflo, lbfhi, f_kl)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 blen, ibflo, ibfhi, kbflo, kbfhi, f_ik)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 blen, jbflo, jbfhi, lbflo, lbfhi, f_jl)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 blen, ibflo, ibfhi, lbflo, lbfhi, f_il)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 blen, jbflo, jbfhi, kbflo, kbfhi, f_jk)
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
c
      subroutine new_get_dens_block(
     $     g_dens, ndens, blen, ibflo, ibfhi, jbflo, jbfhi, buf, fac)
      implicit none
      integer  ndens, g_dens, blen, ibflo, ibfhi, jbflo, jbfhi
      integer blen2
      double precision buf(*), fac(ndens)
c
      integer lo(3), hi(3), ld(3)
      integer idens
c
      blen2 = blen*blen
      call dfill(ndens*blen*blen,0d0,buf,1)
      lo(1) = 1
      hi(1) = ndens
      lo(2) = ibflo
      hi(2) = ibfhi
      lo(3) = jbflo
      hi(3) = jbfhi
      ld(1) = ndens
      ld(2) = blen
      call nga_get(g_dens, lo, hi, buf, ld) 
c
      call scale_dens_block(buf,ndens,blen,ibflo,ibfhi,
     *   jbflo,jbfhi,fac)
c
      return
      end
c
      subroutine scale_dens_block(buf,ndens,blen,ibflo,ibfhi,
     *   jbflo,jbfhi,fac)
      implicit none
c
      integer ndens, blen, ibflo, ibfhi, jbflo, jbfhi
      double precision buf(ndens,blen,blen), fac(ndens)
c
      integer idens, i, j
      double precision factor
c
      do idens = 1, ndens
        factor = fac(idens)
        do i = 1, ibfhi-ibflo+1
          do j = 1, jbfhi-jbflo+1
            buf(idens,i,j) = buf(idens,i,j) * factor
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine new_acc_fock_block(
     $   g_fock, nfock, blen, ibflo, ibfhi, jbflo, jbfhi, buf)
      implicit none
      integer nfock, g_fock, blen 
      integer ibflo, ibfhi, jbflo, jbfhi
      double precision buf(*)
c
      integer i, index
      integer lo(3), hi(3), ld(2)
c
      lo(1) = 1
      hi(1) = nfock
      lo(2) = ibflo
      hi(2) = ibfhi
      lo(3) = jbflo
      hi(3) = jbfhi
      ld(1) = nfock
      ld(2) = blen
      call nga_acc(g_fock, lo, hi, buf, ld, 1d0)
c
      end
C
      subroutine new_fock_doit(nint,labels,leri,eri,tol2e,
     $     d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $     ndens, nfock, blen, rbfmap,
     $     ibflo, jbflo, kbflo, lbflo)
      implicit none
c
      integer ndens, nfock, blen
      integer nint
      integer leri
      integer labels(leri,4)
      double precision eri(*)
      double precision tol2e
      double precision
     $     d_ij(ndens,blen,blen), d_kl(ndens,blen,blen), 
     $     d_ik(ndens,blen,blen), d_jl(ndens,blen,blen), 
     $     d_il(ndens,blen,blen), d_jk(ndens,blen,blen)
      double precision
     $     f_ij(nfock,blen,blen), f_kl(nfock,blen,blen), 
     $     f_ik(nfock,blen,blen), f_jl(nfock,blen,blen), 
     $     f_il(nfock,blen,blen), f_jk(nfock,blen,blen)
      integer rbfmap(*)
      integer ibflo, jbflo, kbflo, lbflo
c
      integer integ, i, j, k, l, nf
      double precision g
c
c  There is an assumption that the number of fock matrices is the same
c  as the number of density matrices.  If this is not true, we bomb for
c  now.  If there is a situation where should be different, ...
c
      if (nfock.ne.ndens) call errquit(
     $  'new_fock_doit: nfock not equal to ndens',555)
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
            do nf = 1, nfock
               f_ij(nf,i,j)=f_ij(nf,i,j)+g*d_kl(nf,k,l)
               f_kl(nf,k,l)=f_kl(nf,k,l)+g*d_ij(nf,i,j)
               f_ik(nf,i,k)=f_ik(nf,i,k)+g*d_jl(nf,j,l)
               f_il(nf,i,l)=f_il(nf,i,l)+g*d_jk(nf,j,k)
               f_jl(nf,j,l)=f_jl(nf,j,l)+g*d_ik(nf,i,k)
               f_jk(nf,j,k)=f_jk(nf,j,k)+g*d_il(nf,i,l)
            end do
         end if
      end do
c
      end
      subroutine new_fock_doit_der(nint,labels,leri,eri,tol2e,
     $     d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $     ndens, nfock, blen, rbfmap, bftoat, 
     $     ibflo, jbflo, kbflo, lbflo)
      implicit none
c
c Note that nfock is the number of alpha or beta fock matrices for UHF
c and not the full number of fock matrices (as in the rest of the code)
c
      integer ndens, nfock, blen
      integer nint
      integer leri
      integer labels(leri,4)
      double precision eri(*)
      double precision tol2e
      double precision
     $     d_ij(ndens,blen,blen), d_kl(ndens,blen,blen), 
     $     d_ik(ndens,blen,blen), d_jl(ndens,blen,blen), 
     $     d_il(ndens,blen,blen), d_jk(ndens,blen,blen)
      double precision
     $     f_ij(nfock*ndens,blen,blen), f_kl(nfock*ndens,blen,blen), 
     $     f_ik(nfock*ndens,blen,blen), f_jl(nfock*ndens,blen,blen), 
     $     f_il(nfock*ndens,blen,blen), f_jk(nfock*ndens,blen,blen)
      integer rbfmap(*)
      integer bftoat(*)       ! [input] Map from old basis func to atom
      integer ibflo, jbflo, kbflo, lbflo
c
      integer integ, i, j, k, l, v
      integer nbuf, iiat, iat, nstart
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
c
c At this point we have the atom information and will combine
c integrals with density and put into local fock storage.  The loops
c below reflect that the same integrals go into both alpha and beta
c matrices for UHF.  So we keep track of where we are in the integral 
c buffer, add contributions to the alpha matrices, go back to the
c original place in the integral buffer and add the contribution to the
c beta matrices.  An alternative would be to put the density loop outside
c of the integral loop, but I didn't see any real advantage to this
c at the moment.
c
              nstart = nbuf
              do v = 1, ndens
              nbuf = nstart
              do xyz=nfockder+(nfock*(v-1)),nfockder+2+(nfock*(v-1))
               g = eri(nbuf)
               f_ij(xyz,i,j)=f_ij(xyz,i,j)+g*d_kl(v,k,l)
               f_kl(xyz,k,l)=f_kl(xyz,k,l)+g*d_ij(v,i,j)
               f_ik(xyz,i,k)=f_ik(xyz,i,k)+g*d_jl(v,j,l)
               f_il(xyz,i,l)=f_il(xyz,i,l)+g*d_jk(v,j,k)
               f_jl(xyz,j,l)=f_jl(xyz,j,l)+g*d_ik(v,i,k)
               f_jk(xyz,j,k)=f_jk(xyz,j,k)+g*d_il(v,i,l)
               nbuf = nbuf + 1
              end do
              end do
            end do
c        end if
      end do
c
      end
      subroutine xc_newfock(geom,basis,jfac,kfac,nbf,
     %     nder, ndens, nfock,g_fock,g_dens)
      implicit none
#include "global.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
      integer geom
      integer basis
      double precision jfac(*),kfac(*)
      integer nbf
      integer nder,nfock,ndens
      integer g_fock,g_dens
c
      integer dims(3),chunk(3)
      integer g_xc
      integer ifock
      integer i
      logical oprint
      character*4 scftype
c
      scftype='RHF'
      if(ndens.eq.2) scftype='UHF'
      oprint= util_print('fock_xc',print_debug) 
      if(oprint) then
         write(luout,*) ' newfockxc nder ',nder, ' nfock ',nfock,
     ,        ' ndens ',ndens
         call util_flush(luout)
      endif
      dims(1) = nfock
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'nfockxc:3dxc',chunk, 
     &     g_xc))
     &        call errquit('fock_2e: could not allocate g_xc',555)
      call ga_zero(g_xc)

      call fock_xc(geom, nbf,basis,
     ,      ndens,  g_dens, g_xc,scftype,.true.)
      call ga_add(1d0, g_xc, 1.0d0, g_fock,
     &     g_fock)
      if (oprint) then
         write(luout,*)'*** fock_2e: fock matrix',ifock
         call util_flush(luout)
         call ga_print(g_fock)
      endif
      if (.not.ga_destroy(g_xc)) then
         call errquit('fock_2e: problem destroying test array',1)
      endif

      return
      end
c
      subroutine uhf_fock_setup(g_dens, g_fock, jfac, kfac, ndens,
     &   nfock, nbf, oprint)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
c  This routine takes the multidimensional density and the fock matrices
c  that are set up for alpha and beta and expand them to be alpha and
c  beta/coulomb and exchange (2 -> 4) so that the correct summation can
c  be performed after the shell_fock_build call.
c
      integer g_dens       ! [input/output] density
      integer g_fock       ! [input/output] Fock matrices
      integer ndens        ! [input/output] number of density matrices
      integer nfock        ! [input/output] number of Fock matrices
      integer nbf          ! [input] number of basis functions
      double precision jfac(ndens)  ! [input/output] Coulomb scale factors
      double precision kfac(ndens)  ! [input/output] Exchange scale factors
      logical oprint       ! [input] Print output
c
      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3)
      integer g_dens2, g_fock2
      integer i, noffset, noffset2, noffset3
      double precision zero
      parameter (zero = 0.0D0)
c
c First set up jfac and kfac correctly.  This assumes that there is
c already enough space in jfac and kfac to expand from 2 to 4.
c
      noffset = ndens/2
      noffset2 = noffset*2
      noffset3 = noffset*3
      do i = 1, noffset
        jfac(noffset3+i) = zero
        jfac(noffset2+i) = jfac(noffset+i)
        jfac(noffset+i)  = zero
        kfac(noffset3+i) = kfac(noffset+i)
        kfac(noffset2+i) = zero
        kfac(noffset+i)  = kfac(i)
        kfac(i)          = zero
      enddo
c
c Next expand the density.  We create a new ga and substitute it
c for the original.
c
      dims(1) = ndens * 2
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'Density',chunk,
     &     g_dens2)) call errquit
     &     ('uhf_fock_setup: could not allocate g_dens2',555)
c
c Copy "backwards" so that we don't loose any data
c
      alo(1) = noffset + 1 ! Beta
      ahi(1) = ndens
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = noffset*3 + 1
      bhi(1) = noffset*4
      blo(2) = 1
      bhi(2) = nbf
      blo(3) = 1
      bhi(3) = nbf
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      blo(1) = noffset*2 + 1
      bhi(1) = noffset*3
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      alo(1) = 1       ! Alpha
      ahi(1) = noffset
      blo(1) = noffset + 1
      bhi(1) = ndens
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      blo(1) = 1
      bhi(1) = noffset
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
c
c  Now get rid of original density
c
      if (.not. ga_destroy(g_dens)) call errquit
     *   ('uhf_fock_setup: failed to free g_dens', ndens)
c
c  Assign g_dens to be the new matrix and update ndens
c
      g_dens = g_dens2
      ndens = ndens * 2
c
c Now expand the fock matrices doing what we did above for
c the densities.
c
      dims(1) = nfock * 2
      chunk(1) = dims(1)
      if (.not. nga_create (MT_DBL, 3, dims, 'Fock matrices',
     &    chunk, g_fock2)) call errquit
     &     ('uhf_fock_setup: could not allocate g_fock2',555)
c
      noffset = nfock/2
      alo(1) = noffset + 1  ! Beta
      ahi(1) = nfock
      blo(1) = noffset*3 + 1
      bhi(1) = noffset*4
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      blo(1) = noffset*2 + 1
      bhi(1) = noffset*3
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      alo(1) = 1       ! Alpha
      ahi(1) = noffset
      blo(1) = noffset + 1
      bhi(1) = nfock
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
      blo(1) = 1
      bhi(1) = noffset
      call nga_copy_patch('N',g_fock,alo,ahi,g_fock2,blo,bhi)
c
      if (.not. ga_destroy(g_fock)) call errquit
     *   ('uhf_fock_setup: failed to free g_fock', nfock)
c
      g_fock = g_fock2
      nfock = nfock * 2
c
      if (oprint) then
        if (ga_nodeid().eq.0) write(LuOut,*)
     &     'Density and Fock matrices after uhf_fock_setup'
        call ga_print(g_dens)
        call ga_print(g_fock)
      endif
c
      return
      end
c
      subroutine uhf_fock_finish(g_dens,g_fock,jfac, kfac, ndens,
     &       nfock, nbf, oprint)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
c  This routine takes the results of shell_fock_build and forms the
c  complete alpha and beta Fock matrices (4 -> 2).
c
      integer g_dens       ! [input/output] density
      integer g_fock       ! [input/output] Fock matrices
      integer ndens        ! [input/output] number of density matrices
      integer nfock        ! [input/output] number of Fock matrices
      integer nbf          ! [input] number of basis functions
      double precision jfac(ndens)  ! [input/output] Coulomb scale factors
      double precision kfac(ndens)  ! [input/output] Exchange scale factors
      logical oprint       ! [input] Print output
c
      integer dims(3), chunk(3)
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer g_dens2, g_fock2
      integer i, noffset, noffset2, noffset3
      double precision zero, one, mone
      parameter (zero =  0.0D0)
      parameter (one  =  1.0D0)
      parameter (mone = -1.0D0)
c
c  Put jfac and kfac back to their original form.
c
      ndens = ndens/2
      noffset = ndens/2
      noffset2 = noffset*2
      noffset3 = noffset*3
      do i = 1, noffset
        jfac(noffset+i)  = jfac(noffset2+i)
        kfac(i)          = kfac(noffset+i)
        kfac(noffset+i)  = kfac(noffset3+i)
      enddo
c
c  Put the density back to its original.  We create a new ga and 
c  substitute it for the original.
c
      dims(1) = ndens
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'Density',chunk,
     &     g_dens2)) call errquit
     &     ('uhf_fock_finish: could not allocate g_dens2',555)
c
c Copy data
c
      alo(1) = noffset*3 + 1  ! Beta
      ahi(1) = noffset*4
      alo(2) = 1
      ahi(2) = nbf
      alo(3) = 1
      ahi(3) = nbf
      blo(1) = noffset + 1
      bhi(1) = ndens
      blo(2) = 1
      bhi(2) = nbf
      blo(3) = 1
      bhi(3) = nbf
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
      alo(1) = 1   ! Alpha
      ahi(1) = noffset
      blo(1) = 1
      bhi(1) = noffset
      call nga_copy_patch('N',g_dens,alo,ahi,g_dens2,blo,bhi)
c
c  Now get rid of original density
c
      if (.not. ga_destroy(g_dens)) call errquit
     *   ('uhf_fock_finish: failed to free g_dens', ndens)
c
c  Assign g_dens to be the new matrix and update ndens
c
      g_dens = g_dens2
c
c  Now get the Fock matrices set up.
c
      nfock = nfock/2
      dims(1) = nfock
      chunk(1) = dims(1)
      if (.not. nga_create (MT_DBL, 3, dims, 'Fock matrices',
     &    chunk, g_fock2)) call errquit
     &     ('uhf_fock_setup: could not allocate g_fock2',555)
c
      noffset = nfock/2
      noffset2 = noffset*2
      noffset3 = noffset*3
      alo(1) = 1  ! Coulomb part (add alpha and beta)
      ahi(1) = noffset
      blo(1) = noffset*2 + 1
      bhi(1) = noffset*3
      clo(1) = 1
      chi(1) = noffset
      clo(2) = 1
      chi(2) = nbf
      clo(3) = 1
      chi(3) = nbf
      call nga_add_patch(one,g_fock,alo,ahi,
     &                    one,g_fock,blo,bhi,
     &                    g_fock2,clo,chi)
      alo(1) = noffset + 1
      ahi(1) = nfock
      call nga_copy_patch('N',g_fock2,clo,chi,
     &    g_fock2,alo,ahi)
c
c  Add exchange of same spin
c
      alo(1) = noffset + 1   ! Alpha
      ahi(1) = nfock
      blo(1) = 1
      bhi(1) = noffset
      call nga_add_patch(mone,g_fock,alo,ahi,
     &                   one,g_fock2,blo,bhi,
     &                    g_fock2,blo,bhi)
c
      alo(1) = noffset*3 + 1 ! Beta
      ahi(1) = noffset*4
      blo(1) = noffset + 1
      bhi(1) = nfock
      call nga_add_patch(mone,g_fock,alo,ahi,
     &                   one,g_fock2,blo,bhi,
     &                    g_fock2,blo,bhi)
c
c  Now get rid of g_fock and assign the new one (g_fock2) to g_fock
c
      if (.not. ga_destroy(g_fock)) call errquit
     &   ('uhf_fock_setup: failed to free g_fock', nfock)
c
      g_fock = g_fock2
c
      if (oprint) then
        if (ga_nodeid().eq.0) write(LuOut,*)
     &     'Density and Fock matrices after uhf_fock_finish'
        call ga_print(g_dens)
        call ga_print(g_fock)
      endif
c
      return
      end
c
      subroutine print_dens_bit(d_ij, ndens, ibflo, ibfhi,
     $                    jbflo, jbfhi, blen)
      implicit none
c
      integer ndens, ibflo, ibfhi, jbflo, jbfhi, blen
      double precision d_ij(ndens,blen,blen)
      integer ndns, id, jd
c
      write(6,*) 'density bit'
      do ndns = 1, ndens
        write(6,*) 'ndns is ',ndns
        do id = 1, 1+ ibfhi -ibflo
          do jd = 1, 1+ jbfhi -jbflo
            write(6,*) id,jd,d_ij(ndns,id,jd)
          enddo
        enddo
      enddo
c
      return
      end
