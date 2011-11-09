c===============================================================c
c                                                               c
c       NCC - a New Coupled-Cluster code for NWChem             c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c               Karol Kowalski                                  c
c               Environmental Molecular Sciences Laboratory     c
c               Pacific Northwest National Laboratory           c
c               karol.kowalski@pnl.gov                          c
c                                                               c
c               Marta Włoch                                     c
c               Department of Chemistry                         c
c               Michigan Technological University               c
c               wloch@mtu.edu                                   c
c                                                               c
c===============================================================c
c
#ifdef VACUUM
#define DEBUG_PRINT
      program nwchem
      implicit none
      integer x2info(2)
      integer nocc, nvir, tszocc, tszvir
      character*(20) label
c
      nocc = 13
      nvir = 57
      tszocc = 7
      tszvir = 17
c
      call ncc_doubles_create_aa(x2info,   ! output array
     1                           label,    ! character name for this array
     2                           nocc,     ! number of occupied orbitals
     3                           nvir,     ! number of virtual orbitals
     4                           tszocc,   ! tilesize of occupied orbitals
     5                           tszvir)   ! tilesize of virtual orbitals
c
      end
c
c this removes the dependency on stdio.fh
c
#define LuOut 6
c
c this removes the dependency on errquit.fh
c
#define errquit errbomb
#define GA_ERR 100
#define MA_ERR 200
c
      subroutine errbomb(string, icode, errcode)
      implicit none
      character*(*) string
      integer icode
      integer errcode
      call ga_error(string, icode)
      end
#endif

#ifdef DEBUG_PRINT
#define DP(x) print*,' x = ', x
#else
#define DP(x)
#endif
c
      subroutine ncc_doubles_create_aa(x2info,   ! output array
     1                                 label,    ! character name for this array
     2                                 nocc,     ! number of occupied orbitals
     3                                 nvir,     ! number of virtual orbitals
     4                                 tszocc,   ! tilesize of occupied orbitals
     5                                 tszvir)   ! tilesize of virtual orbitals
c
c $Id$
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#ifndef VACUUM
#include "errquit.fh"
#include "stdio.fh"
#endif
c
c     interface variables
c
      integer x2info(*)  ! array containing doubles descriptor:
c
c                       x2info(1) = GA handle ~ g_x2
c                       x2info(2) = GA size   ~ s_x2
c                       x2info(3) = 
c
c     other handles, such as for check-pointing, should be 
c     added to this array, so be careful to not hard-core
c     the length too many places
c
      character*(*) label  ! label
      integer nocc         ! number of occupied orbitals
      integer nvir         ! number of virtual orbitals
      integer tszocc       ! tilesize of occupied orbitals
      integer tszvir       ! tilesize of virtual orbitals
c
c     internal variables
c
      integer g_x2         ! GA handle
      integer s_x2         ! GA size
c
      integer nftocc       ! number of full occupied tiles (1D)
      integer nftvir       ! number of full virtual tiles (1D)
c
      integer ntocc       ! total number of occupied tiles (1D)
      integer ntvir       ! total number of virtual tiles (1D)
c
      integer etszocc      ! end tilesize of occupied orbitals
      integer etszvir      ! end tilesize of virtual orbitals
c
      integer ifetocc      ! 1 if end tile exists, 0 otherwise
      integer ifetvir      ! 1 if end tile exists, 0 otherwise
c
c     virtual 2D tile configuration info
c
      integer p1vv_num, p1vv_dim, p1vv_tot
      integer p2vv_num, p2vv_dim, p2vv_tot
      integer p3vv_num, p3vv_dim, p3vv_tot
      integer p4vv_num, p4vv_dim, p4vv_tot
c
c     occupied 2D tile configuration info
c
      integer p1oo_num, p1oo_dim, p1oo_tot
      integer p2oo_num, p2oo_dim, p2oo_tot
      integer p3oo_num, p3oo_dim, p3oo_tot
      integer p4oo_num, p4oo_dim, p4oo_tot
c
c     GA info
c
      integer dims(2)      ! GA dimensions
      integer chunk(2)     ! GA chunking
      integer ndim         ! using 2D GA storage
      parameter (ndim = 2)
      integer gatype       ! numerical type for GA (always double)
      parameter (gatype = MT_DBL)
      integer pgroup       ! GA processor group handle
c
c     irregular distribution variables
c
      integer dist_block(2)             ! blocks per dimension
      integer k_blksz_list,l_blksz_list ! MA handles for map helper
      integer k_dist_map,l_dist_map     ! MA handles for map
      integer b, bsum
c
c     function declarations
c
      integer ncc_anti, ncc_symm
      external ncc_anti, ncc_symm
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_doubles_create_aa'
      DP(tszocc)
      DP(tszvir)
#endif
c
c     determine number of full tiles (dimension tszocc^2 * tszvir^2)
c
      nftocc = nocc / tszocc
      nftvir = nvir / tszvir
c
c     determine size of end tilesize
c
      etszocc = mod(nocc,tszocc)
      etszvir = mod(nvir,tszvir)
c
c     determine total tile number from full tiles + end tile (if exist)
c
      ntocc = nftocc
      ntvir = nftvir
      if (etszocc.gt.0) ntocc = ntocc+1
      if (etszvir.gt.0) ntvir = ntvir+1
      DP(nocc)
      DP(nvir)
      DP(tszocc)
      DP(tszvir)
      DP(etszocc)
      DP(etszvir)
      DP(nftocc)
      DP(nftvir)
      DP(ntocc)
      DP(ntvir)
c
c     configure tiling in four parts - virtuals
c
c       B_last is the end tile index, if it exists, otherwise it is ntvir+1
c
c       part 1 - A < B < B_last (all cases)
c       part 2 - A < B = B_last (if end tile exists)
c       part 3 - A = B < B_last (all cases)
c       part 4 - A = B = B_last (if end tile exists)
c
c     now we tabulate the number of 2D VV tiles,
c     their individual 2D dimension, and their total 2D dimension
c
c     part 1 - A < B < B_last (all cases)
c
      p1vv_num = ncc_anti(nftvir)
      p1vv_dim = tszvir*tszvir
      p1vv_tot = p1vv_num * p1vv_dim
      DP(p1vv_num)
      DP(p1vv_dim)
      DP(p1vv_tot)
c
c     part 2 - A < B = B_last (if end tile exists)
c
      if (etszocc.gt.0) then
          p2vv_num = nftvir
          p2vv_dim = tszvir*etszvir
      else
          p2vv_num = 0
          p2vv_dim = 0
      endif
      p2vv_tot = p2vv_num * p2vv_dim
      DP(p2vv_num)
      DP(p2vv_dim)
      DP(p2vv_tot)
c
c     part 3 - A = B < B_last (all cases)
c
      p3vv_num = nftvir
      p3vv_dim = ncc_anti(tszvir)
      p3vv_tot = p3vv_num * p3vv_dim
      DP(p3vv_num)
      DP(p3vv_dim)
      DP(p3vv_tot)
c
c     part 4 - A = B = B_last (if end tile exists)
c
      if (etszocc.gt.0) then
          p4vv_num = 1
          p4vv_dim = ncc_anti(etszvir)
      else
          p4vv_num = 0
          p4vv_dim = 0
      endif
      p4vv_tot = p4vv_num * p4vv_dim
      DP(p4vv_num)
      DP(p4vv_dim)
      DP(p4vv_tot)
c
c     sanity check, since tiled total 2D dimension should
c     equal the untiled 2D dimension
c
      if ( (p1vv_tot + p2vv_tot + p3vv_tot + p4vv_tot)
     1     .ne. ( ncc_anti(nvir) ) ) then
          if (ga_nodeid().eq.0) write(LuOut,100) ncc_anti(nvir),
     1              (p1vv_tot + p2vv_tot + p3vv_tot + p4vv_tot)
 100  format(1x,'virtual tiling wrong!!! ',/,
     1       1x,'untiled VV dimension = ',i16,/,
     2       1x,'  tiled VV dimension = ',i16)
      endif
      DP(ncc_anti(nvir))
c
c     configure tiling in four parts - occupied
c
c       J_last is the end tile index, if it exists, otherwise it is ntocc+1
c
c       part 1 - I < J < J_last (all cases)
c       part 2 - I < J = J_last (if end tile exists)
c       part 3 - I = J < J_last (all cases)
c       part 4 - I = J = J_last (if end tile exists)
c
c     now we tabulate the number of 2D VV tiles,
c     their individual 2D dimension, and their total 2D dimension
c
c     part 1 - I < J < J_last (all cases)
c
      p1oo_num = ncc_anti(nftocc)
      p1oo_dim = tszocc*tszocc
      p1oo_tot = p1oo_num * p1oo_dim
      DP(p1oo_num)
      DP(p1oo_dim)
      DP(p1oo_tot)
c
c     part 2 - I < J = J_last (if end tile exists)
c
      if (etszocc.gt.0) then
          p2oo_num = nftocc
          p2oo_dim = tszocc*etszocc
      else
          p2oo_num = 0
          p2oo_dim = 0
      endif
      p2oo_tot = p2oo_num * p2oo_dim
      DP(p2oo_num)
      DP(p2oo_dim)
      DP(p2oo_tot)
c
c     part 3 - I = J < J_last (all cases)
c
      p3oo_num = nftocc
      p3oo_dim = ncc_anti(tszocc)
      p3oo_tot = p3oo_num * p3oo_dim
      DP(p3oo_num)
      DP(p3oo_dim)
      DP(p3oo_tot)
c
c     part 4 - I = J = J_last (if end tile exists)
c
      if (etszocc.gt.0) then
          p4oo_num = 1
          p4oo_dim = ncc_anti(etszocc)
      else
          p4oo_num = 0
          p4oo_dim = 0
      endif
      p4oo_tot = p4oo_num * p4oo_dim
      DP(p4oo_num)
      DP(p4oo_dim)
      DP(p4oo_tot)
c
c     sanity check, since tiled total 2D dimension should
c     equal the untiled 2D dimension
c
      if ( (p1oo_tot + p2oo_tot + p3oo_tot + p4oo_tot)
     1     .ne. ( ncc_anti(nocc) ) ) then
          if (ga_nodeid().eq.0) write(LuOut,200) ncc_anti(nocc),
     1              (p1oo_tot + p2oo_tot + p3oo_tot + p4oo_tot)
 200  format(1x,'occupied tiling wrong!!! ',/,
     1       1x,'untiled OO dimension = ',i16,/,
     2       1x,'  tiled OO dimension = ',i16)
      endif
      DP(ncc_anti(nocc))
c
c     end of tiling configuration
c
c     the GA will be 2D nocc*(nocc-1)/2 * nvir*(nvir-1)/2
c
      s_x2 = ncc_anti(nvir) * ncc_anti(nocc)
      x2info(2) = s_x2
      if (ga_nodeid().eq.0) write(LuOut,300) s_x2
 300  format(1x,'creating GA of ',i16,' doubles')
c
c     GA handle creation and labeling
c
      g_x2 = ga_create_handle()
      call ga_set_array_name(g_x2, label)
c
c     GA processor group configuration
c        to begin, we'll use the world group
c
      pgroup = ga_pgroup_get_world()
      call ga_set_pgroup(g_x2, pgroup)
c
c     GA dimensions
c
      dims(1) = ncc_anti(nvir) ! leading dimension is "fast" in Fortran GA API
      dims(2) = ncc_anti(nocc)
      call ga_set_data(g_x2, ndim, dims, gatype)
c
c     GA distribution
c        all 4D tiles should be contiguous on a single host
c        and distributed evenly across all nodes
c
c     dist_block(2) - number of blocks each dimension is divided into
c
      dist_block(1) = p1vv_num + p2vv_num + p3vv_num + p4vv_num ! VV first
      dist_block(2) = p1oo_num + p2oo_num + p3oo_num + p4oo_num ! OO second
      DP(dist_block(1))
      DP(dist_block(2))
c
c     dist_map - starting index for each block; the size s is a sum of all elements of nblock array
c
      if (.not.ma_push_get(mt_int,dist_block(1)+dist_block(2),
     1                     'blksz_list',l_blksz_list,k_blksz_list)) then
          call errquit ('ncc_doubles_create_aa: ma_push_get blksz_list',
     1                  dist_block(1)+dist_block(2),MA_ERR)
      endif
c
      if (.not.ma_push_get(mt_int,dist_block(1)+dist_block(2),
     1                     'dist_map',l_dist_map,k_dist_map)) then
          call errquit ('ncc_doubles_create_aa: ma_push_get dist_map',
     1                  dist_block(1)+dist_block(2),MA_ERR)
      endif
c
      do b = 1, p1vv_num
          int_mb(k_blksz_list+b-1) = p1vv_dim
      enddo
      do b = 1,p2vv_num
          int_mb(k_blksz_list+p1vv_num+b-1) = p2vv_dim
      enddo
      do b = 1,p3vv_num
          int_mb(k_blksz_list+p1vv_num+p2vv_num+b-1) = p3vv_dim
      enddo
      do b = 1,p4vv_num
          int_mb(k_blksz_list+p1vv_num+p2vv_num+p3vv_num+b-1) = p4vv_dim
      enddo
c
      bsum = 1
      int_mb(k_dist_map) = bsum
      do b = 2,p1vv_num+p2vv_num+p3vv_num+p4vv_num
          bsum = bsum + int_mb(k_blksz_list+b-2)
          int_mb(k_dist_map+b-1) = bsum
      enddo
c
#ifdef DEBUG_PRINT
      print*,'========================================'
      do b = 1,dist_block(1)
          print*,'blksz_list',b,int_mb(k_blksz_list+b-1)
      enddo
c
      do b = 1,dist_block(1)
          print*,'dist_map',b,int_mb(k_dist_map+b-1)
      enddo
      print*,'========================================'
#endif
c
      do b = 1, p1oo_num
          int_mb(k_blksz_list+dist_block(1)+b-1) = p1oo_dim
      enddo
      do b = 1,p2oo_num
          int_mb(k_blksz_list+dist_block(1)+p1oo_num+b-1) = p2oo_dim
      enddo
      do b = 1,p3oo_num
          int_mb(k_blksz_list+dist_block(1)+
     1           p1oo_num+p2oo_num+b-1) = p3oo_dim
      enddo
      do b = 1,p4oo_num
          int_mb(k_blksz_list+dist_block(1)+
     1           p1oo_num+p2oo_num+p3oo_num+b-1) = p4oo_dim
      enddo
c
      bsum = 1
      int_mb(k_dist_map+dist_block(1)) = bsum
      do b = 2,p1oo_num+p2oo_num+p3oo_num+p4oo_num
          bsum = bsum + int_mb(k_blksz_list+dist_block(1)+b-2)
          int_mb(k_dist_map+dist_block(1)+b-1) = bsum
      enddo
c
#ifdef DEBUG_PRINT
      print*,'========================================'
      do b = dist_block(1)+1,dist_block(1)+dist_block(2)
          print*,'blksz_list',b,int_mb(k_blksz_list+b-1)
      enddo
c
      do b = dist_block(1)+1,dist_block(1)+dist_block(2)
          print*,'dist_map',b,int_mb(k_dist_map+b-1)
      enddo
      print*,'========================================'
#endif
c
c     send the irregular distribution to GA
c
!      call ga_set_irreg_distr(g_x2, int_mb(k_dist_map), dist_block)
c
c     don't need the dist_map anymore since it is inside of GA now
c
      if (.not.ma_pop_stack(l_dist_map)) then
        call errquit('ncc_doubles_create_aa: ma_pop_stack dist_map ',
     1               0,MA_ERR)
      endif
c
      if (.not.ma_pop_stack(l_blksz_list)) then
        call errquit('ncc_doubles_create_aa: ma_pop_stack blksz_list ',
     1               0,MA_ERR)
      endif
c
c     Regular distribution with chunking for now
c
      chunk(1) = -1
      chunk(2) = -1
      call ga_set_chunk(g_x2, chunk)
c
c     GA allocation
c
      if (.not. ga_allocate(g_x2) ) then
          call errquit ('ncc_doubles_create_aa: ga_allocate',g_x2,
     1                  GA_ERR)
      endif
c
      x2info(1) = g_x2
c
      call ga_zero(g_x2)
c
#ifdef DEBUG_PRINT
      call ga_print_distribution(g_x2)
#endif
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_doubles_create_aa'
#endif
c
      return
      end




      subroutine ncc_doubles_destroy(x2info)
c
      implicit none
c#include "mafdecls.fh"
#include "global.fh"
#ifndef VACUUM
#include "errquit.fh"
#include "stdio.fh"
#endif
c
c     interface variables
c
      integer x2info(*)  ! array containing doubles descriptor:
c
c                       x2info(1) = GA handle ~ g_x2
c                       x2info(2) = GA size   ~ s_x2
c
c     other handles, such as for check-pointing, should be 
c     added to this array, so be careful to not hard-core
c     the length too many places
c
c     internal variables
c
      integer g_x2 ! GA handle
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_doubles_destroy'
#endif
c
      g_x2 = x2info(1)
c
      if (.not. ga_destroy(g_x2) ) then
          call errquit ('ncc_doubles_destroy: ga_destroy',g_x2,GA_ERR)
      endif
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_doubles_destroy'
#endif
c
      return
      end