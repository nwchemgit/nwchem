czora...
czora...Scale zora eigenvalues and energy
czora...
          subroutine dft_zora_scale(geom, 
     &                   ao_bas_han,
     &                   nbf,
     &                   nbf_ao,
     &                   g_dens,
     &                   g_s,
     &                   g_movecs,
     &                   g_zora_scal_sf,
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
      integer g_movecs(2)
      integer g_s
      integer g_zora_scal_sf(4)

      double precision numelecs, ener_kin
      integer g_orb
      integer g_orb_dens
      integer orb_mem_hndl
      integer orb_mem_indx
      integer iorb
      double precision eval_scal
      double precision ener_scal
      integer noc(2)
      integer ispin
      integer ipol

      integer geom
      integer ao_bas_han
      integer nbf,nbf_ao

      double precision focc(nbf*ipol)  ! occupation no.
      double precision evals(ipol*nbf) ! eigenvalues
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     debugging
c     separate the density
c      call ga_dadd(1.d0,g_dens(1),-1.d0,g_dens(2),g_dens(1))
c
c     electron number test
c      numelecs = ga_ddot(g_dens(1),g_s)
c      write(luout,*) "Number of spin up electrons:",numelecs
c      if (ipol.gt.1) then
c       numelecs = ga_ddot(g_dens(2),g_s)
c       write(luout,*) "Number of spin down electrons:",numelecs
c      end if
c
c     kinetic energy test
c      ener_kin = ga_ddot(g_dens(1),g_zora_scal_sf(1))   
c      write(luout,*) "Kinetic energy spin up:",ener_kin
c      if (ipol.gt.1) then
c        ener_kin = ga_ddot(g_dens(2),g_zora_scal_sf(2))
c        write(luout,*) "Kinetic energy spin down:",ener_kin
c      end if
c
c     consolidate the density before leaving
c     call ga_dadd(1.d0,g_dens(1),1.d0,g_dens(2),g_dens(1))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     allocate the necessary arrays
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'vec aux',
     &            orb_mem_hndl, orb_mem_indx))
     &     call errquit('xc_sic: cannot allocate vec aux',0, MA_ERR)
      g_orb = ga_create_atom_blocked(geom, ao_bas_han, 'orbs')
      g_orb_dens = ga_create_atom_blocked(geom, ao_bas_han,'orbs dens')
c
c     scale the eigenvalues and energy
      ener_scal = 0.d0
      do ispin=1,ipol
       do iorb = 1,nbf

        call ga_get(g_movecs(ispin), 1, nbf, iorb, iorb,
     &                      dbl_mb(orb_mem_indx), nbf)
        call ga_zero(g_orb)
        call ga_put(g_orb, 1, nbf, iorb, iorb,
     &                    dbl_mb(orb_mem_indx), nbf)

        call ga_zero(g_orb_dens)
        call ga_dgemm('n', 't', nbf, nbf, noc(ispin), 1.0d00,
     &                      g_orb, g_orb, 0.d00, g_orb_dens)

        eval_scal = evals(iorb)
        if (ispin.gt.1) eval_scal = evals(iorb + nbf)

        eval_scal 
     &    = eval_scal/(1.d0 + ga_ddot(g_orb_dens,g_zora_scal_sf(ispin)))

        ener_scal = ener_scal
     &     - eval_scal*ga_ddot(g_orb_dens,g_zora_scal_sf(ispin))

c       ener_scal = ener_scal + ga_ddot(g_orb_dens,g_zora_scal_sf(ispin))

        if (ispin.le.1) then
            evals(iorb) = eval_scal
c            write(luout,*) iorb,evals(iorb),eval_scal
        else
            evals(iorb+nbf) = eval_scal
c            write(luout,*) iorb,evals(iorb+nbf),eval_scal
        end if

       end do ! orbital loop
      end do ! polarization loop
c
c     double the energy for closed-shell calculations
      if (ipol.eq.1) ener_scal = 2.d0*ener_scal
      
c      write(luout,*) "scaled_energy:",ener_scal

c     deallocate the memory
      if (.not.ma_pop_stack(orb_mem_hndl))
     &   call errquit('dft_scf: cannot pop stack',0, MA_ERR)
      if (.not. ga_destroy(g_orb)) call errquit(
     &  'zora_scale_evals: ga_destroy failed ',0, GA_ERR)
      if (.not. ga_destroy(g_orb_dens)) call errquit(
     &  'zora_scale_evals: ga_destroy failed ',0, GA_ERR)

      return
      end
