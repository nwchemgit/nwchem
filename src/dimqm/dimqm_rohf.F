      subroutine dimqm_rohf_wrap(rtdb, geom, basis, nbf, g_dens,
     $                           g_dim)
c
c     Author: Justin Moore
c
c     Called from: rohf_fock.F
c
c     Calls: dimqm_seeding.F dimqm_indDipoles.F dimqm_fock_xc.F
c
c     Wrapper routine for calculating the DIM dipoles and DIM potential
c     from the ground state SCF in order to build the unperturbed part
c     of the Fock matrix for response.
c
      implicit none
#include "dimqm.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
c     Input Variables
      integer rtdb
      integer geom
      integer basis
      integer nbf
      integer g_dens
c
c     Ouput Variables
      integer g_dim
c
c     Local Variables
      double precision Edim
c      
c     Global Variables from dimqm.fh
c
c     logical lfirst

c     Calculate dipoles.  Mark this is as the first DIM cycle
c     since the arrays were cleaned up after ending the ground state
c     SCF and need reinitialized.  Also mark this as the last DIM cycle
c     since these arrays won't be needed in response.
      lfirst = .true.

      if(.not. rtdb_put(rtdb, 'dimqm:lfirst', mt_log, 1, lfirst))
     $   call errquit('dimqm_rohf_wrap: lfirst put failed', 1, RTDB_ERR)

      call dimqm_main(rtdb, basis, geom, g_dens, .true.)
      if (.not.ga_create(MT_DBL, nbf, nbf, 'DIM FOCK', -1, -1, g_dim))
     $  call errquit('dimqm_rohf_wrap: create dim fock failed', 1,
     $               GA_ERR)
      call ga_zero(g_dim)
      Edim = 0.0d0
      call dim_grid_quadv0(rtdb, g_dim, 1, Edim)

      end subroutine dimqm_rohf_wrap

      subroutine dimqm_rohf_calcDIM(rtdb, geom, basis, nbf, nfock,
     $                              g_dens, g_dim, luhf_fd)
c
c     Author: Justin Moore
c
c     Called from: rohf_hessv2.F
c
c     Calls: dimqm_seeding.F dimqm_indDipoles.F dimqm_fock_xc.F
c
c     Wrapper routine for calculating the DIM dipoles and DIM potential
c     within rohf_hessv_2e2
c

      implicit none
#include "dimqm.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

c     Input Variables
      integer rtdb  
      integer geom
      integer basis
      integer nbf
      integer nfock
      integer g_dens
c
c     Output Variables
      integer g_dim
c
c     Local Variables
      integer rmax
      integer dims(3), chunk(3)
      logical luhf_fd
c
c     Global Variables from dimqm.fh
c
c     logical ldimqm, lrsp, lfirst
c     

c
      ! Check to see DIM is on for response
      if (ldimqm .and. lrsp) then

        ! Pull in residual and check seeding status
        if (.not.lfirst) then
          if (.not.rtdb_get(rtdb, 'lkain:rmax', mt_dbl, 1, rmax))
     $      call errquit('dimqm_rohf: rmax get failed', 1, RTDB_ERR)
          call dimqm_check_dipoles(1.0d0, rmax)
        end if
        
        ! Calculate DIM dipoles
!  jbecca START
        if (dimqm_on .and. luhf_fd ) 
     $    call dimqm_indDipoles(rtdb, nbf, geom, basis, g_dens, .true.)
        if (dimqm_on .and. (.not. luhf_fd))
     $    call dimqm_indDipoles(rtdb, nbf, geom, basis, g_dens, .false.)
        
!        if (dimqm_on) call dimqm_indDipoles(rtdb, nbf, geom, basis,
!     $                                      g_dens)
!jbecca END

        ! Build potential from dipoles
        dims(1)  = nfock
        dims(2)  = nbf
        dims(3)  = nbf
        chunk(1) = dims(1)
        chunk(2) = -1
        chunk(3) = -1
        if (.not. nga_create (MT_DBL, 3, dims, 'dimqm:3ddim',chunk,
     $                        g_dim))
     $    call errquit('dimqm_rohf: could not allocate g_dim',555,
     $                 GA_ERR)
        call ga_zero(g_dim)  
        call fock_dim(geom, nbf, basis, nfock, g_dim, 0, 1)
        call ga_symmetrize(g_dim)
        call ga_sync()
        
        if (.not. rtdb_put(rtdb, 'dimqm:g_dim', mt_int, 1, g_dim))
     $    call errquit('dimqm_rohf: put g_dim failed', 1, RTDB_ERR)
      end if

      end subroutine dimqm_rohf_calcDIM


      subroutine dimqm_rohf_addDIM(rtdb, g_ax, ivec, xoff, voff,
     $                             nclosed, nvir, nmo, nbf,
     $                             g_dcv, g_movec, g_tmp1, g_tmp2,
     $                             alo, ahi, blo, bhi)
c
c     Author: Justin Moore
c
c     Called from: rohf_hessv2.F
c
c     Calls: None
c
c     Wrapper routine for adding the DIM potential to the Fock matrix
c     generated in rohf_hessv_2e2
c

      implicit none
#include "dimqm.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "dimqm_constants.fh"

c     Input Variables
      integer rtdb
      integer g_ax
      integer ivec
      integer xoff
      integer voff
      integer nclosed
      integer nvir
      integer nmo
      integer nbf
      integer g_dcv
      integer g_movec
      integer g_tmp1
      integer g_tmp2
      integer alo
      integer ahi
      integer blo
      integer bhi

c     Local Variables
      integer g_dim

c     Global Variables from dimqm.fh
c
c     logical ldimqm, lrsp

      if (ldimqm .and. lrsp) then
        call nga_copy_patch('N',g_dim,alo,ahi,g_dcv,blo,bhi)
        call ga_scale(g_dcv, TWO)
        call ga_matmul_patch('n', 'n', TWO, zero,
     $                       g_dcv,   1, nbf, 1, nbf,
     $                       g_movec, 1, nbf, 1, nclosed,
     $                       g_tmp1,  1, nbf, 1, nclosed)
        call ga_matmul_patch('t', 'n', one, zero,
     $                       g_movec, voff, nmo, 1, nbf,
     $                       g_tmp1, 1, nbf,  1, nclosed,
     $                       g_tmp2, 1, nvir, 1, nclosed)
        call ga_sync()
        call ga_mat_to_vec(g_tmp2, 1, nvir, 1, nclosed, g_ax,
     $                     xoff, ivec, FOUR, '+')
      end if

      end subroutine dimqm_rohf_addDIM

      subroutine dimqm_rohf_cleanup(rtdb)
c
c     Author: Justin Moore
c
c     Called from: rohf_hessv2.F
c
c     Calls: None
c
c     Destroys g_dim when done in rohf_hessv_2e2
c
      implicit none
#include "dimqm.fh"
#include "rtdb.fh"
#include "mafdecls.fh"

c     Input Variables
      integer rtdb



      end subroutine dimqm_rohf_cleanup
