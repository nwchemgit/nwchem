czora...
czora...Scale zora eigenvalues and energy
czora...
          subroutine dft_zora_scale_so(geom, 
     &                   ao_bas_han,
     &                   nbf,
     &                   nbf_ao,
     &                   nbf_mo,
     &                   g_dens,
     &                   g_s,
     &                   g_moso,
     &                   g_zora_scal_sf,
     &                   g_zora_scal_so,
     &                   evals,
     &                   focc,
     &                   noc,
     &                   ipol,
     &                   ener_scal)

       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer g_dens(2)
      integer g_moso(2)
      integer g_orb(2)
      integer g_dens_so(2)
      integer g_scr(2)
      integer g_s
      integer g_zora_scal_sf(2)
      integer g_zora_scal_so(3)

      double precision numelecs, ener_kin

      integer l_vecsre, k_vecsre
      integer l_vecsim, k_vecsim

      integer iorb
      double precision eval_scal
      double precision ener_scal
      double precision ener_sf
      double precision ener_so_x
      double precision ener_so_y
      double precision ener_so_z
      double precision ener_tot
      integer noc
      integer ispin
      integer ipol

      integer geom
      integer ao_bas_han
      integer nbf
      integer nbf_ao
      integer nbf_mo

      double precision focc(*)  ! occupation no.
      double precision evals(*) ! eigenvalues
c
c     allocate memory
c
      if (.not.MA_Push_Get(MT_Dbl, nbf_mo, 'real vec aux',
     &          l_vecsre, k_vecsre))
     & call errquit('dft_zora_scale_so: cannot allocate vec',0, MA_ERR)
c
      if (.not.MA_Push_Get(MT_Dbl, nbf_mo, 'imag vec aux',
     &          l_vecsim, k_vecsim))
     & call errquit('dft_zora_scale_so: cannot allocate vec',0, MA_ERR)
c 
c     spin-orbit vector - real
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'orbs re',0,0,g_orb(1)))
     & call errquit('dft_zora_scale_so: orb real error',0, GA_ERR)
      call ga_zero(g_orb(1))
c
c     spin-orbit vector - imag
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'orbs im',0,0,g_orb(2)))
     & call errquit('dft_zora_scale_so: orb imag 2 error',0, GA_ERR)
      call ga_zero(g_orb(2))
c
c     spin-orbit density matrix - real
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'denmxre',0,0,
     & g_dens_so(1)))
     & call errquit('dft_zora_scale_so: dens real error',0, GA_ERR)
      call ga_zero(g_dens_so(1))
c
c     spin-orbit density matrix - imag
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'denmxim',0,0,
     & g_dens_so(2)))
     & call errquit('dft_zora_scale_so: dens imag error',0, GA_ERR)
      call ga_zero(g_dens_so(2))
c
c     scratch array
       if(.not.ga_duplicate(g_dens(1),g_scr(1),'scratch 1'))
     &  call errquit('dft_zora_scale_so: ga_duplicate failed',1, GA_ERR)
        call ga_zero(g_scr(1))
       if(.not.ga_duplicate(g_dens(2),g_scr(2),'scratch 2'))
     &  call errquit('dft_zora_scale_so: ga_duplicate failed',1, GA_ERR)
        call ga_zero(g_scr(2))
c
czora...debug
c      call dft_densm_so(g_dens_so, g_moso, nbf_ao, noc)
c      call ga_zero(g_scratch(1))
c      call ga_zero(g_scratch(2))
c      call ga_dens_sf(g_scratch, g_dens_so(1), nbf_ao)
c      ener_kin = ga_ddot(g_scratch(1),g_zora_scal_sf(1))
c      write(luout,*) "Kinetic energy spin up again:",ener_kin
c      if (ipol.gt.1) then
c        ener_kin = ga_ddot(g_scratch(2),g_zora_scal_sf(2))
c        write(luout,*) "Kinetic energy spin down again:",ener_kin
c      end if
czora...debug

      ener_scal = 0.d0
      ener_tot = 0.d0
      do iorb = 1, nbf_mo
c      write(luout,*) iorb, evals(iorb)
       call ga_get(g_moso(1),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsre),1)
       call ga_zero(g_orb(1))
       call ga_put(g_orb(1),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsre),1)

       call ga_get(g_moso(2),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsim),1)
       call ga_zero(g_orb(2))
       call ga_put(g_orb(2),1,nbf_mo,iorb,iorb,dbl_mb(k_vecsim),1)

       call dft_densm_so(g_dens_so, g_orb, nbf_ao, noc)
       call ga_zero(g_scr(1))
       call ga_zero(g_scr(2))
       call ga_dens_sf(g_scr, g_dens_so(1), nbf_ao)

       ener_sf = ga_ddot(g_scr(1),g_zora_scal_sf(1)) +
     &        ga_ddot(g_scr(2),g_zora_scal_sf(2))

       call ga_zero(g_scr(1))
       call ga_dens_so(g_scr(1), g_dens_so, nbf_ao, 'z')
       ener_so_z = ga_ddot(g_scr(1), g_zora_scal_so(1)) 

       call ga_zero(g_scr(1))
       call ga_dens_so(g_scr(1), g_dens_so, nbf_ao, 'y')
       ener_so_y = ga_ddot(g_scr(1), g_zora_scal_so(2)) 

       call ga_zero(g_scr(1))
       call ga_dens_so(g_scr(1), g_dens_so, nbf_ao, 'x')
       ener_so_x = ga_ddot(g_scr(1), g_zora_scal_so(3)) 

       ener_tot = ener_sf + ener_so_z + ener_so_y + ener_so_x
       eval_scal = evals(iorb)
       eval_scal = eval_scal/(1.d0 + ener_tot) 
       ener_scal = ener_scal - eval_scal*ener_tot 

c       ener_scal = ener_scal 
c     &    + ener_sf + ener_so_z + ener_so_y + ener_so_x

       evals(iorb) = eval_scal
c      write(luout,*) iorb, evals(iorb), eval_scal

      end do 
c      write(luout,*) "ener_scal:",ener_scal

c     deallocate memory
      if (.not. ma_chop_stack(l_vecsim)) 
     & call errquit('dft_zora_scale_so:l_vecsim', 0, MA_ERR)
      if (.not. ma_chop_stack(l_vecsre))
     & call errquit('dft_zora_scale_so:l_vecsre', 0, MA_ERR)

      if (.not. ga_destroy(g_orb(1))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_orb(2))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)

      if (.not. ga_destroy(g_dens_so(1))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_dens_so(2))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)

      if (.not. ga_destroy(g_scr(1))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)
      if (.not. ga_destroy(g_scr(2))) 
     & call errquit('dft_zora_scale_so: ga_destroy failed',0, GA_ERR)

      return
      end
