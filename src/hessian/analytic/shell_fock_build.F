c#define NBACC 1
c
c $Id$
c
C> \ingroup hess
C> @{
C>
C> \file shell_fock_build.F
C> Fock matrix derivatives
C>
C> \brief Construct Fock matrix derivatives
C>
C> This routine constructs Fock matrix derivatives of one or more
C> Fock matrices. The Fock matrices may be symmetric or anti-symmetric.
C> The latter occurs, for example, when the Fock matrix is the
C> imaginary part of a matrix.
C>
      subroutine shell_fock_build(geom,basis, nder, ndens,
     $     jfac, kfac, tol2e, g_dens, g_fock, osym)
c
c Took oactive out for the time being, but will probably need to
c bring it back - tlw
c
c osym=.false. now means that the density matrix is anti-symmetric
c and so, at the end of this subroutine, the output fock matrix
c will be anti-symmetrized - sh, 2/15/02
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "prop.fh"
#include "stdio.fh"
c
      integer basis  !< [Input] The basis set handle
      integer nder   !< [Input] The number of derivatives
                     !< - nder = 0: Energy
                     !< - nder = 3: E.g. electric field
                     !< - nder = 3*natoms: Hessian CPHF RHS
      integer ndens  !< [Input] The number of density matrices
      integer nfock  !< [Input] The number of Fock matrices
      integer geom   !< [Input] The geometry handle
      integer g_dens !< [Input] The GA handle for the density matrices
      integer g_fock !< [Input] The GA handle for the Fock matrices
c
      double precision jfac(*) !< [Input] The Coulomb scale factor
      double precision kfac(*) !< [Input] The exchange scale factor
      double precision tol2e   !< [Input] The 2-electron integral
                               !< tolerance
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
      logical osym !< [Input] .True. if the Fock matrices are symmetric,
                   !< otherwise the Fock matrices are anti-symmetric
      parameter (maxblen=36, maxsh=10)
      parameter (maxq=maxsh**4)
      integer ga_create_atom_blocked
      logical xc_gotxc
      external ga_create_atom_blocked,xc_gotxc
c
c     Quietly modify the values of jfac and kfac so that the original
c     values coming in can be consistent with the original fock_2e
c
      do i= 1, ndens
         jfac(i) = jfac(i)*4.0d00
         kfac(i) = kfac(i)*2.0d00
      enddo
c
c     Get the dimensions of g_fock and eventually g_dens.
c
      call nga_inquire(g_fock, ityp, ndim, dims)
      nfock = dims(1)
c
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('new_fock: could not get nbf',0, BASIS_ERR)
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('new_fock: could not get nsh',0, BASIS_ERR)
      if (.not. bas_nbf_ce_max(basis,max_at_bf))
     $     call errquit('new_fock: could not get max_at_bf',0,
     &       BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('new_fock: could not get max_sh_bf',0,
     &       BASIS_ERR)
c
      blen = min(nbf,maxblen,6*max_sh_bf) ! d(6)*6=36, 6**4=1296 quartets
c
      lend = blen*blen*ndens
      if (.not. ma_push_get(mt_dbl, lend, 'd_ij', l_d_ij, k_d_ij)) 
     $     call errquit('new_fock: d_ij', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_kl', l_d_kl, k_d_kl)) 
     $     call errquit('new_fock: d_kl', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_ik', l_d_ik, k_d_ik)) 
     $     call errquit('new_fock: d_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_il', l_d_il, k_d_il)) 
     $     call errquit('new_fock: d_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jk', l_d_jk, k_d_jk)) 
     $     call errquit('new_fock: d_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'd_jl', l_d_jl, k_d_jl)) 
     $     call errquit('new_fock: d_jl', lend, MA_ERR)
c
      lend = blen*blen*nfock
      if (.not. ma_push_get(mt_dbl, lend, 'f_ij', l_f_ij, k_f_ij)) 
     $     call errquit('new_fock: f_ij', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_kl', l_f_kl, k_f_kl)) 
     $     call errquit('new_fock: f_kl', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_ik', l_f_ik, k_f_ik)) 
     $     call errquit('new_fock: f_ik', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_il', l_f_il, k_f_il)) 
     $     call errquit('new_fock: f_il', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jk', l_f_jk, k_f_jk)) 
     $     call errquit('new_fock: f_jk', lend, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lend, 'f_jl', l_f_jl, k_f_jl)) 
     $     call errquit('new_fock: f_jl', lend, MA_ERR)
c
      if (.not. ma_push_get(mt_int, nsh, 'shmap', l_shmap, k_shmap))
     $  call errquit('shell_fock_build:could not allocate shmap',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shglo', l_shglo, k_shglo))
     $  call errquit('shell_fock_build:could not allocate blo',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shghi', l_shghi, k_shghi))
     $  call errquit('shell_fock_build:could not allocate bhi',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbflo', l_shbflo, k_shbflo))
     $  call errquit('shell_fock_build:could not allocate bflo',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'shbfhi', l_shbfhi, k_shbfhi))
     $  call errquit('shell_fock_build:could not allocate bfhi',nsh,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfglo', l_bfglo, k_bfglo))
     $  call errquit('shell_fock_build:could not allocate blo',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfghi', l_bfghi, k_bfghi))
     $  call errquit('shell_fock_build:could not allocate bhi',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bfmap', l_bfmap, k_bfmap))
     $  call errquit('shell_fock_build:could not allocate bfmap',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'rbfmap', l_rbfmap, k_rbfmap))
     $  call errquit('shell_fock_build:could not allocate rbfmap',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nbf, 'bftoat', l_bftoat, k_bftoat))
     $  call errquit('shell_fock_build:could not allocate bftoat',nbf,
     &       MA_ERR)
c 
      call grad_shorder(basis, nsh, nbf, maxsh, blen,
     $     nshblocks, int_mb(k_shglo), int_mb(k_shghi),
     &     int_mb(k_shmap),
     $     int_mb(k_bfmap), int_mb(k_rbfmap), int_mb(k_bfglo),
     $     int_mb(k_bfghi), int_mb(k_shbflo), int_mb(k_shbfhi))
      call intb_mem_2e4c(max2e, lscr) ! blocking algorithm
      if (nder.eq.0) then
        leri = max(max2e,1296*100) ! 100 D quartets
      else
        leri = max(max2e*12,1296*100*12) ! 100 D quartets
      endif
c 
      if (.not. ma_push_get(mt_dbl,leri,'eri',l_eri,k_eri)) 
     $     call errquit('new_fock:could not allocate buffer',leri,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lscr,'scr',l_scr,k_scr))
     $     call errquit('new_fock: scratch alloc failed', lscr,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int,4*leri,'labels',l_labels,k_labels))
     $     call errquit('new_fock: could not allocate labels',4*leri,
     &       MA_ERR)
c     
      do i = 1, nbf
         if (.not. bas_bf2ce(basis, i, iat)) call errquit('bf2ce',i,
     &       MA_ERR)
         int_mb(k_bftoat+i-1) = iat
      end do
c
      if (.not. ma_push_get(mt_int, 4*maxq, 'list', l_list, k_list))
     $     call errquit('new_fock:could not allocate list',4*maxq,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('new_fock:could not allocate q4',maxq,
     &       MA_ERR)
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
     $     jfac, kfac)
c
      call nga_reorder(g_dens, .true., int_mb(k_bfmap),
     $     .true., int_mb(k_bfmap))
      call nga_reorder(g_fock, .true., int_mb(k_bfmap),
     $     .true., int_mb(k_bfmap))
c
c     Clean up memory allocated in this routine
c
      if (.not. ma_chop_stack(l_d_ij))
     $  call errquit('shell_fock_build:failed chopping MA stack',555,
     &       MA_ERR)
c
c     Quietly modify the values of jfac and kfac back to the original
c
      do i= 1, ndens
         jfac(i) = jfac(i)/4.0d00
         kfac(i) = kfac(i)/2.0d00
      enddo
c
c     Check the theory type to make sure dft is active 
c     in addition to checking xc_active()
c
      if(use_theory.eq.'dft'.and.xc_gotxc().and.nder.eq.0)then
         call xc_newfock(geom,basis,jfac,kfac,nbf,
     &        nder, ndens, nfock, g_fock,g_dens)
      endif
c
      if (osym) then
         call ga_symmetrize(g_fock)
      else
         call ga_antisymmetrize(g_fock)
      endif
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
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "schwarz.fh"
#include "util.fh"
#include "bas.fh"
#include "rtdb.fh"
c
      integer basis             ! [input] familiar handles
      integer nder              ! [input] No. of derivatives (0 = energy)
      integer ndens             ! [input] No. of density matrices (for UHF)
      integer nfock             ! [input] No. of fock matrices
      integer g_dens, g_fock
      integer blen              ! [input] max dimension of density block
      double precision          ! [scratch] 6 blocks per possible density
     $     d_ij(*),  d_kl(*),  
     $     d_ik(*),  d_jl(*),  
     $     d_il(*),  d_jk(*)
      double precision          ! [scratch] 6 blocks per possible fock
     $     f_ij(*),  f_kl(*),  
     $     f_ik(*),  f_jl(*),  
     $     f_il(*),  f_jk(*)
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
      integer  next, nint, ijklblock,ijkblock,ijblock
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
      integer i0, j0, k0, l0
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
      integer nproc,bledo,bledoij,bledoijk,bledo_dm
      logical firstij,firstkl,firstik,firstjl,firstil,firstjk
      integer handij,handkl,handik,handjl,handil,handjk
      logical dorepl
      integer g_densrep,g_dens0,g_fockrep,g_fock0,dorep_glob
      logical docopy,dozero,dorepd,dorepf,dorepon
      logical util_mirrmat
      external util_mirrmat
#define IJ 1
#define IK 2
#define JK 3
#define IL 4
#define JL 5
#define KL 6

      integer cache_ij(6)
      integer ptr_ij(6)
      integer
     $     nbh_ij,  nbh_kl,  
     $     nbh_ik,  nbh_jl,  
     $     nbh_il,  nbh_jk
      integer cache_hitmiss(2,6)
      integer cache_hit,cache_miss
      integer bgj_get_rtdb_handle
      external bgj_get_rtdb_handle
      integer rtdb_out
      integer jj
      logical skipk
c     
      firstij=.true.
      firstkl=.true.
      firstik=.true.
      firstjl=.true.
      firstil=.true.
      firstjk=.true.
      smax  = schwarz_max()
      skipk=.true.
c     this allows to skip Exchange K Terms
c     ik jk il jl
      do jj=1,ndens
         if(abs(kfac(jj)).gt.1d-8) skipk=.false.
      enddo
      call int_acc_std()
      g_dens0 = 0
      g_fock0 = 0
c
      nproc = ga_nnodes()
      task_size = nshblocks*(nshblocks+1)/2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))
c
c     replicated DM and Fock Matrix
c
c     turned off for nbf gt 3000 since becomes a memory hog
c
      rtdb_out=bgj_get_rtdb_handle()
      dorepon=.true.
      if (.not.rtdb_get(rtdb_out,'fock:mirrmat',mt_log,1,dorepon)) then
         dorepon=nbf.lt.3000
      endif

      dorepl=.false.
      if(ga_cluster_nnodes().gt.1.and.dorepon) then
         docopy=.true.
         dozero=.false.
         dorepd=util_mirrmat(1,g_dens,g_densrep,docopy,dozero)
         docopy=.false.
         dozero=.true.
         dorepf=util_mirrmat(1,g_fock,g_fockrep,docopy,dozero)
         dorepl=dorepd.and.dorepf
         dorep_glob=0
         if(dorepl) dorep_glob=1
         call ga_igop(375,dorep_glob,1, '+')
         dorepl=dorep_glob.eq.ga_nnodes()
         if(dorepl) then
            g_dens0=g_dens
            g_dens=g_densrep
            g_fock0=g_fock
            g_fock=g_fockrep
         else
            if(dorepd) call util_mirrstop(g_densrep)
            if(ga_nodeid().eq.0) then
               write(6,*) ' no mirroring in shellfockbld'
               call util_flush(6)
            endif
         endif
      endif
c      if(ga_nodeid().eq.0) 
c     W     write(6,*) ' dorepon dorepd dorepl ',dorepon,dorepd,dorepl
c     write(6,*) 'nshblocks and task_size ',nshblocks,task_size
c
c     caching to reduce get dm blocks data movement
c
      nbh_ij=0
      nbh_ik=0
      nbh_jk=0
      nbh_il=0
      nbh_jl=0
      nbh_kl=0
      do jj=1,6
         cache_ij(jj)=-1
c first row = hit
c 2nd   row = miss
         cache_hitmiss(1,jj)=0
         cache_hitmiss(2,jj)=0
c  fill ptr with something
         ptr_ij(jj)=0
      enddo
c   get max blk lentgh
      bledo_dm=0
      do ishblock = nshblocks, 1, -1
         ibflo = bfglo(ishblock)
         ibfhi = bfghi(ishblock)
         idim  = ibfhi - ibflo + 1
         bledo_dm=max(bledo_dm,idim)
      enddo
c     
c     Parallel loop over quartets of reordered shell blocks
c     
      ijklblock = 0
      ijkblock=0
      ijblock=0
#define NXT_ 1
#define DISTR_IJK 1
#define DISTR_IJ_ 1
#ifdef NXT
      task_size=1
      ijkblock=-1
#if 1
c     bcast dij i=nshblock,j=1
         ishblock=nshblocks
         ibflo = bfglo(ishblock)
         ibfhi = bfghi(ishblock)
         idim  = ibfhi - ibflo + 1
         jshblock=1
         jbflo = bfglo(jshblock)
         jbfhi = bfghi(jshblock)
         jdim  = jbfhi - jbflo + 1
      if(ga_nodeid().eq.0) then
         call new_get_dens_block(g_dens, ndens, bledo_dm, 
     $        ibflo, ibfhi, jbflo, jbfhi, d_ij, jfac)
      endif
      CALL ga_brdcst(987651,d_ij,
     L     ndens*bledo_dm*jdim*ma_sizeof(mt_dbl,1,mt_byte),0)
      next = nxtask(nproc,task_size)
      ptr_ij(IJ)=(jshblock-1)*nshblocks+ishblock
      cache_ij(IJ)=ptr_ij(IJ)
      cache_hitmiss(2,IJ)=cache_hitmiss(2,IJ)+1
#endif      
      do ishblock = nshblocks, 1, -1
#else
c      do ishblock = 1, nshblocks
      do i0 = ga_nodeid(), nshblocks+ga_nodeid()-1
         ishblock = mod(i0,nshblocks)+1
#endif
         ishlo = shglo(ishblock)
         ishhi = shghi(ishblock)
         ibflo = bfglo(ishblock)
         ibfhi = bfghi(ishblock)
         idim  = ibfhi - ibflo + 1
#ifdef NXT
         do jshblock = 1, ishblock
#else
         do j0 = ga_nodeid(), ishblock+ga_nodeid()-1
            jshblock=mod(j0,ishblock)+1
#endif
          ptr_ij(IJ)=(jshblock-1)*nshblocks+ishblock
            jshlo = shglo(jshblock)
            jshhi = shghi(jshblock)
            jbflo = bfglo(jshblock)
            jbfhi = bfghi(jshblock)
            jdim  = jbfhi - jbflo + 1
            bledoij=max(idim,jdim)
#ifdef NXT
            do kshblock = ishblock, 1, -1
#else
               ijblock= (ishblock*(ishblock-1))/2+jshblock-1
c     distribute work on ij
#ifdef DISTR_IJ
               if (mod(ijblock,ga_nnodes()).eq.ga_nodeid()) then
#endif
            do k0 = ga_nodeid(), ishblock+ga_nodeid()-1
               kshblock=mod(k0,ishblock)+1
#endif
               ptr_ij(IK)=(kshblock-1)*nshblocks+ishblock
               ptr_ij(JK)=(kshblock-1)*nshblocks+jshblock
               kshlo = shglo(kshblock)
               kshhi = shghi(kshblock)
               kbflo = bfglo(kshblock)
               kbfhi = bfghi(kshblock)
               kdim  = kbfhi - kbflo + 1
               bledoijk=max(bledoij,kdim)
               lblockhi = kshblock
#ifndef NXT
               ijkblock=(jshblock-1)*ishblock+kshblock-1
               do jj=1,ishblock-1
                  ijkblock=ijkblock+jj*jj
               enddo
#else
               ijkblock=ijkblock+1
               if (next .eq. ijkblock) then
cstatic load balancing
c     distribute work on ijk
#endif
#ifdef DISTR_IJK
               if (mod(ijkblock,ga_nnodes()).eq.ga_nodeid()) then
#endif

               if (ishblock .eq. kshblock) lblockhi = jshblock
c     
                     if(ptr_ij(IJ).eq.cache_ij(IJ)) then
                        cache_hitmiss(1,IJ)=cache_hitmiss(1,IJ)+1
                     else
                        call new_get_dens_block(g_dens, ndens, bledo_dm, 
     $                       ibflo, ibfhi, jbflo, jbfhi, d_ij, jfac)
                        cache_ij(IJ)=ptr_ij(IJ)
                        cache_hitmiss(2,IJ)=cache_hitmiss(2,IJ)+1
                     endif
                     if(.not.skipk) then
                     if(ptr_ij(IK).eq.cache_ij(IK)) then
                        cache_hitmiss(1,IK)=cache_hitmiss(1,IK)+1
                     else
                        call new_get_dens_block(g_dens, ndens, bledo_dm, 
     $                       ibflo, ibfhi, kbflo, kbfhi, d_ik, kfac)
                        cache_ij(IK)=ptr_ij(IK)
                        cache_hitmiss(2,IK)=cache_hitmiss(2,IK)+1
                     endif
                     if(ptr_ij(JK).eq.cache_ij(JK)) then
                        cache_hitmiss(1,JK)=cache_hitmiss(1,JK)+1
                     else
                        call new_get_dens_block(g_dens, ndens, bledo_dm, 
     $                       jbflo, jbfhi, kbflo, kbfhi, d_jk, kfac)
                        cache_ij(JK)=ptr_ij(JK)
                        cache_hitmiss(2,JK)=cache_hitmiss(2,JK)+1
                     endif
                     endif !skipk
c               do lshblock = 1, lblockhi
               do l0 = ga_nodeid(), lblockhi+ga_nodeid()-1
                  lshblock=mod(l0,lblockhi)+1
                  ptr_ij(KL)=(lshblock-1)*nshblocks+kshblock
                  ptr_ij(JL)=(lshblock-1)*nshblocks+jshblock
                  ptr_ij(IL)=(lshblock-1)*nshblocks+ishblock
                  lshlo = shglo(lshblock)
                  lshhi = shghi(lshblock)
                  lbflo = bfglo(lshblock)
                  lbfhi = bfghi(lshblock)
                  ldim  = lbfhi - lbflo + 1
c                  if (next .eq. ijklblock) then
                     oij = ishblock.eq.jshblock
                     okl = kshblock.eq.lshblock
                     oikjl = ishblock.eq.kshblock .and. 
     $                    jshblock.eq.lshblock
                     bledo=max(bledoijk,ldim)
c     
c     Get blocks of the one-particle densities
                     if(.not.skipk) then
                     if(ptr_ij(IL).eq.cache_ij(IL)) then
                        cache_hitmiss(1,IL)=cache_hitmiss(1,IL)+1
                     else
                        call new_get_dens_block(g_dens, ndens, bledo_dm, 
     $                       ibflo, ibfhi, lbflo, lbfhi, d_il, kfac)
                        cache_ij(IL)=ptr_ij(IL)
                        cache_hitmiss(2,IL)=cache_hitmiss(2,IL)+1
                     endif
                     if(ptr_ij(JL).eq.cache_ij(JL)) then
                        cache_hitmiss(1,JL)=cache_hitmiss(1,JL)+1
                     else
                        call new_get_dens_block(g_dens, ndens, bledo_dm, 
     $                       jbflo, jbfhi, lbflo, lbfhi, d_jl, kfac)

                        cache_ij(JL)=ptr_ij(JL)
                        cache_hitmiss(2,JL)=cache_hitmiss(2,JL)+1
                     endif
                     endif ! skipk
                     if(ptr_ij(KL).eq.cache_ij(KL)) then
                        cache_hitmiss(1,KL)=cache_hitmiss(1,KL)+1
!                     elseif(.false.) then
                     elseif(oikjl) then
c   ij.eq.kl   copy dm_ij to dm_kl
                        call dcopy(ndens*bledo_dm*ldim,
     D                       d_ij,1,d_kl,1)
                        cache_hitmiss(1,KL)=cache_hitmiss(1,KL)+1
                        cache_ij(KL)=ptr_ij(KL)
                     else
#define NBG_KL 1
#ifdef NBG_KL
                     call new_nbget_dens_block(g_dens, ndens, bledo_dm, 
     $                       kbflo, kbfhi, lbflo, lbfhi, 
     D                       d_kl, jfac,nbh_kl)
#else
                        call new_get_dens_block(g_dens, ndens, bledo_dm, 
     $                       kbflo, kbfhi, lbflo, lbfhi, d_kl, jfac)
#endif
                        cache_ij(KL)=ptr_ij(KL)
                        cache_hitmiss(2,KL)=cache_hitmiss(2,KL)+1
                     endif

                     call hess_zerofij(f_ij,bledo,nfock,idim,jdim)
                     call hess_zerofij(f_kl,bledo,nfock,kdim,ldim)
                     if(.not.skipk) then
                        call hess_zerofij(f_jl,bledo,nfock,jdim,ldim)
                        call hess_zerofij(f_ik,bledo,nfock,idim,kdim)
                        call hess_zerofij(f_il,bledo,nfock,idim,ldim)
                        call hess_zerofij(f_jk,bledo,nfock,jdim,kdim)
                     endif

c
c     
c     Build the list of integral shell quartets in NWChem shell labelling
c     
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
c    $                    call errquit(' ',0, MA_ERR)
c     
 200                             continue
                              end do ! next llsh
                           end do
 100                       continue
                        end do  ! next jjsh
                     end do
c
                     if (nder.eq.0) then
                       if (.not. intb_init4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri, 
     $                    block_eff)) call errquit
     $                    ('twodd_coul_ex:txs init?',nq, UNKNOWN_ERR)
                     else
                       if (.not. intbd_init4c(
     $                    basis, list(1,1), list(1,2),
     $                    basis, list(1,3), list(1,4),
     $                    nq, q4, .true., lscr, scr, leri,
     $                    block_eff)) call errquit
     $                    ('twodd_coul_ex:txs init?',nq, UNKNOWN_ERR)
                     endif
c                    if (.not. ma_verify_allocator_stuff())
c    $                    call errquit(' ',0, MA_ERR)
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
c    $                    call errquit(' ',0, MA_ERR)
c     
#ifdef NBG_KL
                     if(nbh_kl.ne.0) then
                        call ga_nbwait(nbh_kl)
                        call scale_dens_block(d_kl,ndens,bledo_dm,
     $                       kbflo, kbfhi, lbflo, lbfhi, 
     *                       jfac)
                     endif
                     nbh_kl=0
#endif
                     if (nint .gt. 0) then
                       if (nder.eq.0) then
                          call new_fock_doit(nint,labels,leri,eri,tol2e,
     $                    d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $                    ndens, nfock, bledo, bledo_dm,
     R                         rbfmap,
     $                    ibflo, jbflo, kbflo, lbflo, skipk)
                        else
                          call new_fock_doit_der(nint,labels,leri,
     $                    eri,tol2e,
     $                    d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $                    f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $                    ndens, nfock/ndens, bledo, bledo_dm,
     R                          rbfmap, bftoat,
     $                    ibflo, jbflo, kbflo, lbflo, skipk)
c                                      
                        endif
                     endif

c                     if (.not. ma_verify_allocator_stuff())
c     $                    call errquit(' ',0, MA_ERR)

                     if (omore) goto 10
c                                      
                     call new_acc_fock_block(g_fock, nfock, 
     $                 bledo, ibflo, ibfhi, jbflo, jbfhi, f_ij,
     F                    handij, firstij)
                     if(.not.skipk) then
                     call new_acc_fock_block(g_fock, nfock, 
     $                 bledo, jbflo, jbfhi, lbflo, lbfhi, f_jl,
     F                    handjl, firstjl)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 bledo, jbflo, jbfhi, kbflo, kbfhi, f_jk,
     F                    handjk, firstjk)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 bledo, ibflo, ibfhi, lbflo, lbfhi, f_il,
     F                    handil, firstil)
                     call new_acc_fock_block(g_fock, nfock, 
     $                 bledo, ibflo, ibfhi, kbflo, kbfhi, f_ik,
     F                    handik, firstik)
                     endif
                     call new_acc_fock_block(g_fock, nfock, 
     $                 bledo, kbflo, kbfhi, lbflo, lbfhi, f_kl,
     F                    handkl, firstkl)
c
c     
                  ijklblock = ijklblock + 1
c     
               end do ! lshblock
#ifdef DISTR_IJK
                  end if ! work distr on ijk
#endif
#ifdef NXT
                  next = nxtask(nproc,task_size)
                  end if ! work distr on nxtask
#endif
            end do ! kshblock
#ifdef DISTR_IJ
            end if ! work distr on ij
#endif
         end do
      end do
c
      call int_acc_std()
c
#ifdef NXT
      next = nxtask(-nproc,task_size)
#endif
      call ga_igop(65331, cache_hitmiss, 12, '+')
      cache_hit=0
      cache_miss=0
      do jj=1,6
         cache_hit=cache_hit+cache_hitmiss(1,jj)
         cache_miss=cache_miss+cache_hitmiss(2,jj)
      enddo
#ifdef PRINT_CACHE
      if(ga_nodeid().eq.0) then
         if((cache_hit+cache_miss).ne.0) then
            write(6,14) ' cache hit: ',100d0*cache_hit/
     C           (cache_hit+cache_miss), '% out of ',
     C           cache_hit+cache_miss, ' ga_gets'
         else
            write(6,*) 'cm ',cache_hitmiss
         endif
         if((cache_hitmiss(1,IJ)+cache_hitmiss(2,IJ)).gt.0)
     W        write(6,14) ' cacij hit: ',100d0*cache_hitmiss(1,IJ)/
     C     (cache_hitmiss(1,IJ)+cache_hitmiss(2,IJ)), '% out of ',
     C     cache_hitmiss(1,IJ)+cache_hitmiss(2,IJ), ' ga_gets'
         if((cache_hitmiss(1,IK)+cache_hitmiss(2,IK)).gt.0)
     W        write(6,14) ' cacik hit: ',100d0*cache_hitmiss(1,IK)/
     C     (cache_hitmiss(1,IK)+cache_hitmiss(2,IK)), '% out of ',
     C     cache_hitmiss(1,IK)+cache_hitmiss(2,IK), ' ga_gets'
         if((cache_hitmiss(1,JK)+cache_hitmiss(2,JK)).gt.0)
     W        write(6,14) ' cacjk hit: ',100d0*cache_hitmiss(1,3)/
     C     (cache_hitmiss(1,3)+cache_hitmiss(2,3)), '% out of ',
     C     cache_hitmiss(1,3)+cache_hitmiss(2,3), ' ga_gets'
         if((cache_hitmiss(1,IL)+cache_hitmiss(2,IL)).gt.0)
     W        write(6,14)  ' cacil hit: ',100d0*cache_hitmiss(1,4)/
     C     (cache_hitmiss(1,4)+cache_hitmiss(2,4)), '% out of ',
     C     cache_hitmiss(1,4)+cache_hitmiss(2,4), ' ga_gets'
         if((cache_hitmiss(1,JL)+cache_hitmiss(2,JL)).gt.0)
     W        write(6,14)  ' cacjl hit: ',100d0*cache_hitmiss(1,5)/
     C     (cache_hitmiss(1,5)+cache_hitmiss(2,5)), '% out of ',
     C     cache_hitmiss(1,5)+cache_hitmiss(2,5), ' ga_gets'
         if((cache_hitmiss(1,KL)+cache_hitmiss(2,KL)).gt.0)
     W        write(6,14)  ' cackl hit: ',100d0*cache_hitmiss(1,6)/
     C     (cache_hitmiss(1,6)+cache_hitmiss(2,6)), '% out of ',
     C     cache_hitmiss(1,6)+cache_hitmiss(2,6), ' ga_gets'
      endif
 14   format(a,f8.2,a,i12,a)
#endif
      call ga_sync()
      if(dorepl) then
         call util_mirrstop(g_densrep)
         g_dens=g_dens0
         call util_mirrmrgadd(g_fockrep,g_fock0)
         g_fock=g_fock0
      endif
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
      double precision buf(*), fac(ndens)
c
      integer lo(3), hi(3), ld(3)
      integer idens
c
cedo      call dcopy(ndens*blen*blen,0d0,0,buf,1)
      lo(1) = 1
      hi(1) = ndens
      lo(2) = ibflo
      hi(2) = ibfhi
      lo(3) = jbflo
      hi(3) = jbfhi
      ld(1) = ndens
      ld(2) = blen
c
      call nga_get(g_dens, lo, hi, buf, ld) 
      call scale_dens_block(buf,ndens,blen,ibflo,ibfhi,
     *   jbflo,jbfhi,fac)
c
      return
      end
      subroutine new_nbget_dens_block(
     $     g_dens, ndens, blen, ibflo, ibfhi, jbflo, jbfhi, buf, fac,
     N     nbhandle)
      implicit none
      integer  ndens, g_dens, blen, ibflo, ibfhi, jbflo, jbfhi
      double precision buf(*)
      double precision fac(ndens) ! [ignored]
      integer nbhandle ! [out]
c
      integer lo(3), hi(3), ld(3)
      integer idens
c
cedo      call dcopy(ndens*blen*blen,0d0,0,buf,1)
      lo(1) = 1
      hi(1) = ndens
      lo(2) = ibflo
      hi(2) = ibfhi
      lo(3) = jbflo
      hi(3) = jbfhi
      ld(1) = ndens
      ld(2) = blen
c
!      call nga_get(g_dens, lo, hi, buf, ld)
      call nga_nbget(g_dens, lo, hi, buf, ld,nbhandle) 
!                  call ga_nbwait(nbhandle)
c      call scale_dens_block(buf,ndens,blen,ibflo,ibfhi,
c     *   jbflo,jbfhi,fac)
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
      if(ndens.eq.1) then
        factor = fac(1)
        do j = 1, jbfhi-jbflo+1
           do i = 1, ibfhi-ibflo+1
              buf(1,i,j) = buf(1,i,j) * factor
           enddo
        enddo
      elseif(ndens.eq.2) then
         do j = 1, jbfhi-jbflo+1
            do i = 1, ibfhi-ibflo+1
               buf(1,i,j) = buf(1,i,j) * fac(1)
               buf(2,i,j) = buf(2,i,j) * fac(2)
            enddo
         enddo
      else
         do j = 1, jbfhi-jbflo+1
            do i = 1, ibfhi-ibflo+1
               do idens = 1, ndens
                  factor = fac(idens)
                  buf(idens,i,j) = buf(idens,i,j) * fac(idens)
               enddo
            enddo
         enddo
      endif
c
      return
      end
c
      subroutine new_acc_fock_block(
     $   g_fock, nfock, blen, ibflo, ibfhi, jbflo, jbfhi, buf,
     H     hand, first)
      implicit none
      integer nfock, g_fock, blen 
      integer ibflo, ibfhi, jbflo, jbfhi
      double precision buf(nfock,blen,*)
c
      integer i, index
      integer lo(3), hi(3), ld(2)
      integer hand
      logical first
      double precision eps
#define FOCK0 1
      integer j,nf
      double precision maxf

      parameter(eps=1d-14)
c     
      integer fock_0,fock_1
      double precision dabsmax
      external dabsmax
c
#ifdef FOCK0
      fock_1=0
      fock_0=nfock+1
      do nf=1,nfock
         maxf=-1d0
c         if(dabsmax(blen*(jbfhi-jbflo+1),buf(i,1,1)).gt.eps)  then
         do j=1,jbfhi-jbflo+1
            do i=1,ibfhi-ibflo+1
               maxf=max(maxf,abs(buf(nf,i,j)))
            enddo
         enddo
         if(maxf.gt.eps)  then
            fock_0=min(nf,fock_0)
            fock_1=max(nf,fock_1)
         endif
      enddo
      if(fock_1.eq.0.and.fock_0.eq.nfock+1) return
      lo(1) = fock_0
      hi(1) = fock_1
#else
      maxf=-1d0
      do j=1,jbfhi-jbflo+1
         do i=1,ibfhi-ibflo+1
            maxf=max(maxf,dabsmax(nfock,buf(1,i,j)))
         enddo
      enddo
      if(maxf.lt.eps)
     R     return
      lo(1) = 1
      hi(1) = nfock
#endif
      lo(2) = ibflo
      hi(2) = ibfhi
      lo(3) = jbflo
      hi(3) = jbfhi
      ld(1) = nfock
      ld(2) = blen
#ifdef NBACC
      if(.not.first) then
         call ga_nbwait(hand)
      endif
      first=.false.
#ifdef FOCK0
      call nga_nbacc(g_fock,lo,hi,buf(fock_0,1,1),ld,1d0,hand)
#else
      call nga_nbacc(g_fock,lo,hi,buf,ld,1d0,hand)
#endif
#else
#ifdef FOCK0
      call nga_acc(g_fock, lo, hi, buf(fock_0,1,1), ld, 1d0)
#else
      call nga_acc(g_fock, lo, hi, buf, ld, 1d0)
#endif
#endif
c
      end
C
      subroutine new_fock_doit(nint,labels,leri,eri,tol2e,
     $     d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $     ndens, nfock, blen, blen_dm, rbfmap,
     $     ibflo, jbflo, kbflo, lbflo, skipk)
      implicit none
#include "errquit.fh"
c
      integer ndens, nfock, blen, blen_dm
      integer nint
      integer leri
      integer labels(leri,4)
      double precision eri(*)
      double precision tol2e
      double precision
     $     d_ij(ndens,blen_dm,*), d_kl(ndens,blen_dm,*), 
     $     d_ik(ndens,blen_dm,*), d_jl(ndens,blen_dm,*), 
     $     d_il(ndens,blen_dm,*), d_jk(ndens,blen_dm,*)
      double precision
     $     f_ij(nfock,blen,*), f_kl(nfock,blen,*), 
     $     f_ik(nfock,blen,*), f_jl(nfock,blen,*), 
     $     f_il(nfock,blen,*), f_jk(nfock,blen,*)
      integer rbfmap(*)
      integer ibflo, jbflo, kbflo, lbflo
      logical skipk ![in]
c
      integer integ, i, j, k, l, nf
      double precision g
c
c  There is an assumption that the number of fock matrices is the same
c  as the number of density matrices.  If this is not true, we bomb for
c  now.  If there is a situation where should be different, ...
c
      if (nfock.ne.ndens) call errquit(
     $  'new_fock_doit: nfock not equal to ndens',555, INPUT_ERR)
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
            if(skipk) then
            do nf = 1, nfock
               f_ij(nf,i,j)=f_ij(nf,i,j)+g*d_kl(nf,k,l)
               f_kl(nf,k,l)=f_kl(nf,k,l)+g*d_ij(nf,i,j)
            end do
            else
            do nf = 1, nfock
               f_ij(nf,i,j)=f_ij(nf,i,j)+g*d_kl(nf,k,l)
               f_kl(nf,k,l)=f_kl(nf,k,l)+g*d_ij(nf,i,j)
               f_ik(nf,i,k)=f_ik(nf,i,k)+g*d_jl(nf,j,l)
               f_il(nf,i,l)=f_il(nf,i,l)+g*d_jk(nf,j,k)
               f_jl(nf,j,l)=f_jl(nf,j,l)+g*d_ik(nf,i,k)
               f_jk(nf,j,k)=f_jk(nf,j,k)+g*d_il(nf,i,l)
            end do
            endif
         end if
      end do
c
      end
      subroutine new_fock_doit_der(nint,labels,leri,eri,tol2e,
     $     d_ij, d_kl, d_ik, d_jl, d_il, d_jk, 
     $     f_ij, f_kl, f_ik, f_jl, f_il, f_jk, 
     $     ndens, nfock, blen, blen_dm, rbfmap, bftoat, 
     $     ibflo, jbflo, kbflo, lbflo, skipk)
      implicit none
c
c Note that nfock is the number of alpha or beta fock matrices for UHF
c and not the full number of fock matrices (as in the rest of the code)
c
      integer ndens, nfock, blen, blen_dm
      integer nint
      integer leri
      integer labels(leri,4)
      double precision eri(*)
      double precision tol2e
      double precision
     $     d_ij(ndens,blen_dm,*), d_kl(ndens,blen_dm,*), 
     $     d_ik(ndens,blen_dm,*), d_jl(ndens,blen_dm,*), 
     $     d_il(ndens,blen_dm,*), d_jk(ndens,blen_dm,*)
      double precision
     $     f_ij(nfock*ndens,blen,*), f_kl(nfock*ndens,blen,*), 
     $     f_ik(nfock*ndens,blen,*), f_jl(nfock*ndens,blen,*), 
     $     f_il(nfock*ndens,blen,*), f_jk(nfock*ndens,blen,*)
      integer rbfmap(*)
      integer bftoat(*)       ! [input] Map from old basis func to atom
      integer ibflo, jbflo, kbflo, lbflo
      logical skipk ! [in]
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
              if(skipk) then
              do xyz=nfockder+(nfock*(v-1)),nfockder+2+(nfock*(v-1))
               g = eri(nbuf)
               nbuf = nbuf + 1
               f_ij(xyz,i,j)=f_ij(xyz,i,j)+g*d_kl(v,k,l)
               f_kl(xyz,k,l)=f_kl(xyz,k,l)+g*d_ij(v,i,j)
              end do
              else
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
              endif
              end do
            end do
c        end if
      end do
c
      end
      subroutine xc_newfock(geom,basis,jfac,kfac,nbf,
     %     nder, ndens, nfock,g_fock,g_dens)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "bgj.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
      integer geom
      integer basis
      double precision jfac(*),kfac(*),tdum,Exc(2)
      integer nExc
      integer nbf
      integer nder,nfock,ndens
      integer g_fock,g_dens
c
      integer dims(3),chunk(3),alo(3),ahi(3)
      integer g_xc
      integer ifock
      integer i
      logical oprint
      integer ilo(3),ihi(3)
      integer ga_create3d_atom_blocked
      external ga_create3d_atom_blocked
c
      oprint= util_print('fock_xc',print_debug) 
      if(oprint) then
         write(luout,*) ' xcnewfock nder ',nder, ' nfock ',nfock,
     ,        ' ndens ',ndens
         call util_flush(luout)
      endif
#if 1
      dims(1) = nfock
      dims(2) = nbf
      dims(3) = nbf
      chunk(1) = dims(1)
      chunk(2) = -1
      chunk(3) = -1
      if (.not. nga_create (MT_DBL, 3, dims, 'nfockxc:3dxc',chunk, 
     &     g_xc))
     &        call errquit('xcnewfock: could not allocate g_xc',555,
     &       GA_ERR)
#else
      g_xc=ga_create3d_atom_blocked(geom, basis,'nfockxc:3dxc', nfock)
      call nga_distribution(g_xc,
     .     ga_nodeid(), ilo,ihi)
      write(6,123) ga_nodeid(),' gxc distr ilo',ilo
      write(6,123) ga_nodeid(),' gxc distr ihi',ihi
 123  format(i2,a,3i5)
#endif

      call ga_zero(g_xc)

      Exc(1) = 0.0d0
      Exc(2) = 0.0d0
      nExc = 1
      call fock_xc(geom, nbf,basis,
     ,             ndens,  g_dens, g_xc,Exc,nExc,.true.)
      call ga_add(1.0d0, g_xc, 1.0d0, g_fock,
     &     g_fock)
      if (oprint) then
         write(luout,*)'*** xcnewfock: fock matrix after add'
         call util_flush(luout)
         alo(2)=1
         ahi(2)=nbf
         alo(3)=1
         ahi(3)=nbf
         do i=1,nfock
            alo(1)=i
            ahi(1)=i
            tdum=nga_ddot_patch(g_fock, 'N', alo, ahi, 
     .           g_fock, 'N', alo, ahi) 
            if(abs(tdum).gt.1d-9) call nga_print_patch(g_fock,alo,ahi,0)
         enddo
      endif
      if (.not.ga_destroy(g_xc)) then
         call errquit('xcnewfock: problem destroying test array',1,
     &       GA_ERR)
      endif

      return
      end
      subroutine hess_zerofij(f_ij,bledo,nfock,idim,jdim)
      implicit none
      integer bledo,nfock,idim,jdim
      double precision f_ij(nfock,bledo,*)
c
      integer j
c
      do j=1,jdim
         call dcopy(nfock*idim, 0d0, 0,f_ij(1,1,j), 1)
      enddo
      return
      end
C> @}
