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
      subroutine ncc_doubles_create(x2info,   ! output array
     1                              label,    ! character name for this array
     2                              nocc,     ! number of occupied orbitals
     3                              tszocc,   ! tilesize of occupied orbitals
     4                              nvir,     ! number of virtual orbitals
     5                              tszvir)   ! tilesize of virtual orbitals
c
c $Id: ncc_doubles.F,v 1.0 2009/20/10 23:48:58 jhammond Exp $
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
c
c     interface variables
c
      integer x2info(*)  ! array containing doubles descriptor:
c
c                       x2info(1) = 
c                       x2info(2) = 
c                       x2info(3) = 
c
c     other handles, such as for check-pointing, should be 
c     added to this array, so be careful to not hard-core
c     the length too many places
c
      character*(*) label  ! label
      integer nocc         ! number of occupied orbitals
      integer tszocc       ! tilesize of occupied orbitals
      integer nvir         ! number of virtual orbitals
      integer tszvir       ! tilesize of virtual orbitals
c
c     internal variables
c
      integer g_x2 ! GA handle
      integer s_x2 ! GA size
c
      integer ntocc    ! number of occupied tiles (1D)
      integer ntvir    ! number of virtual tiles (1D)
      integer ntall    ! total number of tiles (4D)
c
      integer tilesize ! size of 4D tile
      integer dims(2)  ! GA dimensions
      integer chunk(2) ! GA chunking
      integer ndim     ! using 2D GA storage
      parameter (ndim = 2)
      integer gatype   ! numerical type for GA (always double)
      parameter (gatype = MT_DBL)
      integer pgroup   ! GA processor group handle
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_doubles_create'
#endif
c
c     determine number of tiles necessary to hold occ/vir orbitals
c
      ntocc = ceiling( real(nocc) / real(tszocc) )
      ntvir = ceiling( real(nvir) / real(tszvir) )
c
c     doubles amplitudes X^{a,b}_{i,j} satisfy i<j, a<b
c     but need to store tiles for i<=j, a<=b
c
      ntall = ntocc*(ntocc+1)/2 * ntvir*(ntvir+1)/2
c
c     the uniform 4D tile is non-symmetric
c
      tilesize = tszocc * tszocc * tszvir * tszvir
c
c     the GA will be 2D: ntall x tilesize
c
      s_x2 = ntall * tilesize
      if (ga_nodeid().eq.0) write(LuOut,100) s_x2
 100  format(1x,'creating GA of ',i16,' doubles')
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
      dims(1) = tilesize ! leading dimension is "fast" in Fortran
      dims(2) = ntall
      call ga_set_data(g_x2, ndim, dims, gatype)
c
c     GA chunking
c        all 4D tiles should be contiguous on a single host
c        and distributed evenly across all nodes
c
c     this choice probably induces round-robin allocation
c      chunk(1) = 1
c     this choice probably clusters adjacent tiles
      chunk(1) = nint( real(ntall) / real(ga_nnodes()) )
      chunk(2) = tilesize
      call ga_set_chunk(g_x2, chunk)
c
c     GA allocation
c
      if (.not. ga_allocate(g_x2) ) then
          call errquit ('ncc_doubles_create: ga_allocate',g_x2,GA_ERR)
      endif
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_doubles_create'
#endif
c
      return
      end




      subroutine ncc_doubles_destroy(x2info)
c
      implicit none
#include "mafdecls.fh"
c#include "global.fh"
#include "errquit.fh"
c#include "stdio.fh"
c
c     interface variables
c
      integer x2info(*)  ! array containing doubles descriptor:
c
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
      g_x2 = x2info(3)
c

c
#ifdef DEBUG_PRINT
      print*,'end of ncc_doubles_destroy'
#endif
c
      return
      end
