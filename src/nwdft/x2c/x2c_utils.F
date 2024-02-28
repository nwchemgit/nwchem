      subroutine x2c_guess(rtdb,nexc,ipol,
     t     toll_s,
     d     do_x2c,
     g     g_vxc,g_rel,g_dens_at,
     g     g_x2c_u,g_x2c_u2c,
     N     nbf_ao)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
      integer rtdb,nexc,ipol
      double precision toll_s
      logical do_x2c
      integer g_vxc(*) ! [input]
      integer g_rel(*)          ! [output]
      integer g_dens_at(*)
      integer g_x2c_u(*),g_x2c_u2c,g_x2c_c2u
      integer nbf_ao
c
      integer me
c      
c     -------------------------------------------------
c     set up an exact decoupling Hamiltonian (optional)
c     -------------------------------------------------
      me=ga_nodeid()
      if (do_x2c) then
        if (me.eq.0) write (luout,*) 'allocating g_x2c',nbf_ao

        if(.not.ga_duplicate(g_vxc(1),g_rel(1),'x2c 1'))
     &     call errquit('dft_scf: ga_duplicate failed',666, GA_ERR)
        call ga_zero(g_rel(1))

        if (me.eq.0) write (luout,*) 'calling calc_x2c_1e_scalar'

c       note: g_vxc is used as temp array in the numerical
c       integration routines.
c       g_x2c_u2c is created in the x2c routine and survives
c       if keepU=.T. Same for the transformation matrices g_xc_u(:).
c       optionally, g_x2c_c2u is also allocated and kept.

        call calc_x2c_1e_scalar_ga(rtdb, 
     &    nexc, g_dens_at, g_vxc, g_rel(1), g_x2c_u,
     &    g_x2c_u2c, g_x2c_c2u, toll_s)


        if (ipol.gt.1) then
          if(.not.ga_duplicate(g_vxc(2),g_rel(2),'x2c 2'))
     &     call errquit('dft_scf: ga_duplicate failed',1, GA_ERR)
          call ga_copy(g_rel(1),g_rel(2))
        end if
        

      end if                    ! do_x2c

      return
      end

      subroutine x2czora_atguess(rtdb,geom,ao_bas_han,
     i     ipol, nbf_ao, oskel, ntotel, nclosed,
     d     do_zora, do_x2c,
     g     g_dens_at)
      implicit none
      integer rtdb,geom,ao_bas_han
      integer ipol, nbf_ao, oskel, ntotel, nclosed
      logical do_zora, do_x2c
      integer g_dens_at(*)
c
      logical x2c_noatm
      external x2c_noatm

c
c     -------------------------------------------
c     calculate superposition of atomic densities
c     if ZORA or X2C functionality is requested
c     -------------------------------------------

      if ((do_zora .or. do_x2c) .and. .not.x2c_noatm()) then
        call ga_zero(g_dens_at(1))
        if (ipol.gt.1) call ga_zero(g_dens_at(2))
        call guess_dens(rtdb, geom, ao_bas_han, g_dens_at)
        call dft_quickguess(ipol,nbf_ao,g_dens_at,oskel,
     N     ntotel,nclosed,
     D     .false.)
      else if (do_x2c .and. x2c_noatm()) then
        call ga_zero(g_dens_at(1))
        if (ipol.gt.1) call ga_zero(g_dens_at(2))
      end if                    ! do_zora .or. do_x2c
      return
      end
