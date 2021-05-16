      subroutine tce_ao2e_cholesky_2d(rtdb,              ! input
     &                               d_ao2e)            ! output (first full AO matrix, then replaced by full Cholesky factor)
c
c $Id$
c
c     Returns a SF handle for the AO 2-e integral matrix.
c
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sf.fh"
#include "errquit.fh"
c#include "tce.fh"
#include "tce_main.fh"
      integer rtdb                    ! Run-time database
      logical nodezero                ! True if node zero
      integer d_ao2e                  ! AO integral GA handle
      integer work1,work2             ! Work array sizes
      integer l_work1,k_work1         ! Work array 1
      integer l_work2,k_work2         ! Work array 2
      integer l_ijkl,k_ijkl           ! AO integral patch
      integer ish,jsh,ksh,lsh         ! Integral shell indices
      integer i,j,k,l,kl              ! Integral indices
      integer ikjl_offset,ijkl_offset ! 
      integer ilo,ihi,jlo,jhi         ! Integral shell index range
      integer klo,khi,llo,lhi         ! Integral shell index range
      integer irange,jrange           ! Integral shell index range size
      integer krange,lrange           ! Integral shell index range size
      integer kloff,joff,loff         ! Offset temporaries
      integer nxtask
      external nxtask
      integer next,count,nprocs
c      character*255 filename
c
c     New variables: formerly from tce_main.fh, now defined locally.
c     This is not the multigrid code so I can use the global variables - if I want.
c
c      integer geom                          ! determined in this subroutine
c      character*255 basis_name              ! determined in this subroutine
c      integer ao_bas_han                  ! determined in this subroutine
c      integer nbf                           ! determined in this subroutine
c      integer nshells                       ! determined in this subroutine
c      double precision tol2e
c      logical int_normalize
c      external int_normalize
c
c     New variables for the Cholesky code
c
      integer size_aoints
      integer rank_aoints
      integer ijkl_patch_size
      integer ikjl_patch_size
      integer ij_sh_offset_lo,ij_sh_offset_hi
      integer kl_sh_offset_lo,kl_sh_offset_hi
c
c     =========
c     Basis set
c     =========
c
c#ifdef DEBUG_PRINT
c      print*,'initialize basis set'
c#endif
c      basis_name = 'ao basis'
c      if (.not.bas_create(ao_bas_han,basis_name))
c     1  call errquit('tce_ao2e_cholesky_2d: bas_create failed',
c     2               0,BASIS_ERR)
c      if (.not.bas_rtdb_load(rtdb,geom,ao_bas_han,basis_name))
c     1  call errquit('tce_ao2e_cholesky_2d: basis not found',
c     2               0,BASIS_ERR)
c      if (.not.int_normalize(rtdb,ao_bas_han))
c     1  call errquit('tce_ao2e_cholesky_2d: int_normalize failed',
c     2               0,INT_ERR)
c      if (.not.bas_numbf(ao_bas_han,nbf))
c     1  call errquit('tce_ao2e_cholesky_2d: failed to read nbf',
c     2               0,BASIS_ERR)
c      if (.not.bas_numcont(ao_bas_han,nshells))
c     1  call errquit('tce_ao2e_cholesky_2d: failed to read nshells',
c     2               0,BASIS_ERR)
c
c     ===========
c     Allocate GA
c     ===========
c
      nodezero=(ga_nodeid().eq.0)

#ifdef DEBUG_PRINT
      print*,'starting GA allocation'
#endif
c      if (.not.(ioalg.eq.2)) then ! ioalg not accessible variable right now
c        call errquit('tce_ao2e_cholesky_2d: you must use ioalg=2 (GA);
c     1                your chosen ioalg is',ioalg,GA_ERR)
c      endif
      rank_aoints = nbf*nbf
      size_aoints = rank_aoints*rank_aoints
      if (nodezero.and.util_print('ao2e',print_default)) then
        write(LuOut,*)
        write(LuOut,9000) rank_aoints
        write(LuOut,9001) size_aoints
      endif
 9000 format(1x,'rank of AO integral matrix = ',i8)
 9001 format(1x,'size of AO integral matrix = ',i8)
      if (.not.ga_create(mt_dbl,rank_aoints,rank_aoints,'ao2e',-1,-1,
     1  d_ao2e)) then
        write(LuOut,*) ' available GA memory ',
     1    ga_memory_avail(),' bytes' 
        call errquit('tce_ao2e_cholesky_2d: ga_create failed; 
     1                filesize in MB is',size_aoints,GA_ERR)
      endif
c
c     =======================
c     Initialize integral API
c     =======================
c
c#ifdef DEBUG_PRINT
c      print*,'initialize integral API'
c#endif
      call int_mem_2e4c(work1,work2)
      if (.not.ma_push_get(mt_dbl,work1,'work1',l_work1,k_work1))
     1  call errquit('tce_ao2e_cholesky_2d: MA problem work1',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,work2,'work2',l_work2,k_work2))
     1  call errquit('tce_ao2e_cholesky_2d: MA problem work2',1,MA_ERR)
c
c     ================
c     Loop over shells
c     ================
c
#ifdef DEBUG_PRINT
      print*,'begin loop over shells'
#endif
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
c
c     BEGIN KL LOOPS
c
#ifdef DEBUG_PRINT
      print*,'nshells = ',nshells
      do ish=1,nshells
        if (.not.bas_cn2bfr(ao_bas_han,ish,ilo,ihi))
     1    call errquit('tce_ao2e_cholesky_2d: basis fn range problem i',
     2                 ish,BASIS_ERR)
          irange = ihi - ilo + 1
          write(LuOut,9900) 'i',ish,ilo,ihi,irange
      enddo
c      call errquit('manual stop',0,CALC_ERR)
#endif
c
      do lsh=1,nshells
        do ksh=1,nshells
          if (next.eq.count) then
            if (util_print('ao2e',print_debug))
     1        write(LuOut,9190) ksh,lsh,ga_nodeid()
c
            if (.not.bas_cn2bfr(ao_bas_han,ksh,klo,khi))
     1        call errquit('tce_ao2e_cholesky_2d: 
     2        basis fn range problem K',ksh,BASIS_ERR)
            krange = khi - klo + 1
c
            if (.not.bas_cn2bfr(ao_bas_han,lsh,llo,lhi))
     1        call errquit('tce_ao2e_cholesky_2d: 
     2        basis fn range problem L',lsh,BASIS_ERR)
            lrange = lhi - llo + 1
c
            ijkl_patch_size = rank_aoints*krange*lrange
            if (.not.ma_push_get(mt_dbl,ijkl_patch_size,'ao2e',
     1        l_ijkl,k_ijkl)) call errquit('tce_ao2e_cholesky_2d: 
     2        MA problem k_ijkl',ijkl_patch_size,MA_ERR)
c
            do ijkl_offset = 1,ijkl_patch_size
              dbl_mb(k_ijkl+ijkl_offset-1) = 0.0d0
            enddo
c
c           BEGIN IJ LOOPS
c
            do jsh=1,nshells
              do ish=1,nshells
                tol2e = 0.0d0 ! no screening right now
                if (schwarz_shell(ish,ksh)*schwarz_shell(jsh,lsh) 
     1            .ge. tol2e) then
c
                  if (util_print('ao2e',print_debug))
     1              write(LuOut,9191) ish,jsh,ga_nodeid()
c
                  if (.not.bas_cn2bfr(ao_bas_han,ish,ilo,ihi))
     1              call errquit('tce_ao2e_cholesky_2d: 
     2              basis fn range problem I',ish,BASIS_ERR)
                  irange = ihi - ilo + 1
c
                  if (.not.bas_cn2bfr(ao_bas_han,jsh,jlo,jhi))
     1              call errquit('tce_ao2e_cholesky_2d: 
     2              basis fn range problem J',jsh,BASIS_ERR)
                  jrange = jhi - jlo + 1
c
c                 dbl_mb(k_work1) = eri
c                 eri = <bra_g(ish).bra_g(jsh) | ket_g(ksh).ket_g(lsh)>
c
c#ifdef DEBUG_PRINT
c                  print*,'calling int_2e4c'
c#endif
                  call int_2e4c(ao_bas_han,ish,ksh,ao_bas_han,jsh,
     1              lsh,work2,dbl_mb(k_work2),work1,dbl_mb(k_work1))
c
c =============== BEGIN INNER LOOPS ==============
c
c This should be replaced by a tce_sort-like routine.
c
c#ifdef DEBUG_PRINT
c                  print*,'start transpose'
c#endif
                  kl = 0
                  do l = llo,lhi
                    do k = klo,khi
                      kl = kl + 1
                      kloff = (kl-1)*rank_aoints
                      do j = jlo,jhi
                        joff = (j-1)*nbf
                        do i = ilo,ihi
                          ikjl_offset = (
     1                                   (
     2                                    (i-ilo)*krange+(k-klo)
     3                                   )*jrange+(j-jlo)
     4                                  )*lrange+(l-llo)
                          dbl_mb(k_ijkl+kloff+joff+i-1) = 
     1                      dbl_mb(k_work1+ikjl_offset)
                        enddo              ! i
                      enddo                ! j
                    enddo                  ! k
                  enddo                    ! l
c#ifdef DEBUG_PRINT
c                  print*,'end transpose'
c#endif
c
c =============== END INNER LOOPS ==============
c
                endif                      ! schwarz_shell
c
c             Put local stuff into GA
c
              ij_sh_offset_lo = 0
              ij_sh_offset_hi = 0
              kl_sh_offset_lo = 0
              kl_sh_offset_hi = 0
              ijkl_patch_size = 0
#ifdef DEBUG_PRINT
                  write(LuOut,9900) 'i',ish,ilo,ihi,irange
                  write(LuOut,9900) 'j',jsh,jlo,jhi,jrange
                  write(LuOut,9900) 'k',ksh,klo,khi,krange
                  write(LuOut,9900) 'l',lsh,llo,lhi,lrange
c                  print*,'call ga_put'
#endif
c              call ga_put(d_ao2e, ij_sh_offset_lo, ij_sh_offset_hi, 
c     1                    kl_sh_offset_lo, kl_sh_offset_hi, 
c     2                    dbl_mb(k_ijkl), ijkl_patch_size)
              enddo                        ! ish
            enddo                          ! jsh
c
c           END IJ LOOPS
c
            if (.not.ma_pop_stack(l_ijkl)) 
     1        call errquit('tce_ao2e_cholesky_2d: 
     2                     ma_pop_stack of l_ijkl failed',16,MA_ERR)
            next = nxtask(nprocs,1)
          endif
          count = count + 1
        enddo                              ! ksh
      enddo                                ! lsh
c
c     END KL LOOPS
c
      next = nxtask(-nprocs,1)
      call ga_sync() ! only necessary in parallel but that variable not available right now
c
c     ==========
c     Deallocate
c     ==========
c
      if (.not.ma_pop_stack(l_work2)) 
     1  call errquit('tce_ao2e: MA problem',14,MA_ERR)
      if (.not.ma_pop_stack(l_work1)) 
     1  call errquit('tce_ao2e: MA problem',15,MA_ERR)
c
c     ======
c     Return
c     ======
c
 8000 format(1x,'cpu & wall time in seconds for ',A,2f20.1)
 9190 format(/,'k- and l-shells',2i3,' have been assigned to node ',i3)
 9191 format(/,'i- and j-shells',2i3,' have been assigned to node ',i3)
 9200 format(/,'k and l',2i3)
 9900 format(a1,'-shell number ',i3,': lo = ',i3,' hi = ',i3,
     &                                             ' range = ',i3)
      return
      end
