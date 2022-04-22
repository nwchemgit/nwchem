      subroutine dft_fockbld(rtdb,natoms,ntotel,
     ,     g_fock,g_focks,g_dens,g_vxc,g_movecs,g_xcinv,
     ,     g_tmp,g_nq,g_wght,g_xyz,
     ,     ivcoul_Opt,nexc,ivxc_opt,wght_GA,
     ,     n_batch,n3c_int,n3c_dbl,IOLGC,
     ,     evals, rdens_atom,dd3ceri,ii3cwhat,
     ,     cetobfr, iter, ecore, ecoul,exc, rho_n,ldmix,
     ,     test_sic,exact_pot,ecoul_aux,dif_lamda,totsic,
     ,     i_degen,g_rho_exact,g_tmp_exact, g_addit_exact,
     ,     n_levels, act_levels, oprint_time,
cc AJL/Begin/FDE
c     ,     g_rel_sf,g_frozemb)
     ,     g_rel_sf,g_frozemb,g_vxc_fde,g_xcinv_fde,g_dens_fde,eelec )
cc AJL/End
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
#include "cdft.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "dftps.fh"
#include "cosmo.fh"
#include "zora.fh"
#include "x2c.fh"
#include "case.fh"
c#include "dimqm.fh"
c
      integer g_fock,g_dens(2),g_vxc(4),g_movecs(2) ! [handle in/values out]
      integer g_focks(2)
      integer g_xcinv
      integer g_tmp,g_nq,g_wght,g_xyz
      integer g_rho_exact,g_tmp_exact, g_addit_exact
      integer rtdb ! [in]
      integer natoms ! [in]
      logical wght_GA
c     double precision evals (*)
      integer evals(2)
      double precision rdens_atom(*)
      double precision dd3ceri(*)
      integer ii3cwhat(*)
      integer cetobfr(*)
      logical oprint_time ![in]
      integer iter ! [in]
      integer n_levels(2), act_levels
      double precision ecore, ecoul,exc(2) ! [out]
      integer ntotel,n_batch,n3c_int,n3c_dbl
      logical IOLGC, ldmix
      integer test_sic ! [in]
      integer exact_pot, iter_ini
      integer i_degen(2)  ! [in]
      double precision ecoul_aux(2),dif_lamda,totsic
      double precision rho_n ![out]
      double precision fant_a,fant_d,dft_fant
      external dft_fant
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      integer ivcoul_Opt,nexc,ivxc_opt
c
      integer me,nscr,lscr,iscr,lcd_coef,icd_coef
c
      double precision Ebq
c
      integer g_rel_sf(*)       ! spin-free zora/x2c ga
      double precision Erel_sf  ! spin-free zora/x2c energy
c
      integer g_frozemb       ! frozen embedding ga
      double precision Efrozemb  ! energy

      integer g_cosmo

      logical no2e
c
c     timings
c
      double precision time1_ke,time2_ke
      double precision time1_pot,time2_pot
c  jbecca START-- making NWChem great again with efield
      double precision efield(3)
      double precision eelec
      logical lescale
c  jbecca END
c
c MN solvation models -->
c
      double precision gcds
      integer do_cosmo_vem
      logical do_cosmo_smd
c
c DIM/QM variables
c
      double precision g_vdim
      double precision Edim
      logical ldimqm, dimqm_go,use_dimqm
c
cc AJL/Begin/SPIN ECPs
      integer ecp_channels
cc AJL/Begin/FDE 
      integer i
      double precision Exc_fde(2)
      integer g_xcinv_fde
      integer g_vxc_fde(4)
c      double precision rho_n_fde ![out]
      integer g_dens_fde(2)
cc AJL/End
      double precision t0
      integer g_dmr
c
      do_cosmo_smd = .false.
      do_cosmo_vem = 0
      gcds = 0.D0
c
      if (cosmo_on.and.cosmo_phase.eq.2) then
       if (.not. rtdb_get(
     $ rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $ call errquit(
     $ 'dft_fockbld: cannot get do_cosmo_smd from rtdb',
     $ 0,rtdb_err)
       if (.not. rtdb_get(
     $ rtdb,'cosmo:do_cosmo_vem',mt_int,1,do_cosmo_vem))
     $ call errquit(
     $ 'dft_fockbld: cannot get do_cosmo_vem from rtdb',
     $ 0,rtdb_err)
       if (do_cosmo_smd) then
        if (.not. rtdb_get
     $   (rtdb,'cosmo:gcds',mt_dbl,1,gcds))
     $   call errquit(
     $   'dft_fockbld: cannot get gcds from rtdb',
     $   0,rtdb_err)
        endif
      endif
c
c <-- MN solvation models
c
      use_dimqm = util_module_avail("dimqm")
      if (use_dimqm) then
         call dimqm_init(rtdb)
         call dimqm_used(ldimqm)
         call dimqm_getgo(dimqm_go)
      else
         ldimqm=.false.
         dimqm_go=.false.
      endif
c     =========================================================================

      me=ga_nodeid()
      if (odftps) call pstat_on(ps_fockb)

c ... jochen: to run tests on one-electron systems, we
c     simply disable all two-electron contributions if
c     requested by an input switch
      if(.not.rtdb_get(rtdb,'dft:no2e',MT_LOG,1,no2e))
     &   no2e = .false.
      if (no2e.and.me.eq.0) write (luout,*)
     &   'WARNING: dft:no2e set. Disabling Coulomb and XC Fock terms'


c
      Ecoul  = 0d0
      Exc(1) = 0d0
      Exc(2) = 0d0

      icd_coef = 1

c Jorge Garza/Begin
      if (ldmix) then
        iter_ini = 1
      else
        iter_ini = 3
      end if
      if (test_sic.eq.2.or.test_sic.eq.4) then
        if (iter.ge.iter_ini) then
           call xc_sicdeg(i_degen, evals, n_levels,act_levels)
        end if
      end if
c Jorge Garza/End
c     
c     Accumulate core hamiltonian into Fock matrix; compute core energy
c     
      if (oprint_time)
     &     call dft_tstamp(' calling 1-e int ')
      if (odftps) call pstat_on(ps_onel)
      call ga_zero(g_fock)
c
      if(util_print('dft timings', print_high))
     &      time1_ke=util_cpusec()   ! start kinetic time
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'kinetic', oskel)
      if(util_print('dft timings', print_high))
     &      time2_ke=util_cpusec()   ! end kinetic time
c
c     print kinetic energy matrix time
      if(util_print('dft timings', print_high)) then
        if (me.eq.0) then
           write(*,"(4x,'Kinetic Matrix Time:',F13.1,'s')") 
     &        time2_ke-time1_ke
        endif
      end if
c
      if(util_print('kinetic_energy', print_high)) then
         Ecore = ga_ddot(g_dens(1), g_fock)
         if(ga_nodeid().eq.0) then
            write(luout,222) Ecore
 222        format(/
     &     '           Kinetic energy =', f20.12/)
         endif
      endif
c
cc AJL/Begin/SPIN ECPs
cc    Read this value from rtdb vvvv
      if (.not.rtdb_get(rtdb, 'dft:spin_polarised_ecps',
     &  mt_int, 1, ecp_channels)) then
        ecp_channels = 1
      endif

      if (ecp_channels.eq.2.and.ipol.eq.1)
     & call errquit('dft_fockbld: spin_polarised_ecps with mult = 1', 
     &               0, BASIS_ERR)

c     == Copy Fock matrix to individual components for spin ECPs ==
      call ga_copy(g_fock, g_focks(1))
      if (ecp_channels.gt.1) call ga_copy(g_fock, g_focks(2))
c
c     if requested separate calculation of bq-electron energy (MV)
c     ----------------------------------------------------------- 
      if (rtdb_get(rtdb, 'dft:bq_energy', mt_dbl, 1,Ebq)) then
c        Ecore = ga_ddot(g_dens(1), g_fock)
c        call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 
c     &                 'bq_pot', oskel)

        Ecore = ga_ddot(g_dens(1), g_focks(1))
        call int_1e_ga(ao_bas_han, ao_bas_han, g_focks(1),
     &                 'bq_pot', oskel)
c Copy to g_focks(2), as this is identical for both channels
        if (ecp_channels.gt.1) call ga_copy(g_focks(1), g_focks(2))
c
c        Ebq = ga_ddot(g_dens(1), g_fock)-Ecore
c        call int_1e_ga(ao_bas_han, ao_bas_han, g_fock,
c     &                 'potential0', oskel)

        Ebq = ga_ddot(g_dens(1), g_focks(1))-Ecore
        call int_1e_ga(ao_bas_han, ao_bas_han, g_focks(1),
     &                 'potential0', oskel)
c Calculate Beta channel 1e potential
        if (ecp_channels.gt.1)
     &    call int_1e_ga(ao_bas_han, ao_bas_han, g_focks(2),
     &                   'potential_beta0', oskel)
cc AJL/End
c
        if (.not. rtdb_put(rtdb,'dft:bq_energy',mt_dbl,1,Ebq)) 
     $     call errquit('dft: writing bq_energy failed', 
     $     0, RTDB_ERR)
      else
        if(util_print('dft timings', print_high))
     &        time1_pot=util_cpusec()   ! start time
cc AJL/Begin/SPIN ECPs
c        call int_1e_ga(ao_bas_han, ao_bas_han, g_fock,
        call int_1e_ga(ao_bas_han, ao_bas_han, g_focks(1), 
     &                 'potential', oskel)
        if (ecp_channels.gt.1)
     &    call int_1e_ga(ao_bas_han, ao_bas_han, g_focks(2),
     &                   'potential_beta', oskel)
cc AJL/End
        if(util_print('dft timings', print_high))
     &        time2_pot=util_cpusec()   ! end time
c
c     print potential matrix time
        if(util_print('dft timings', print_high)) then
         if (me.eq.0) then
           write(*,"(4x,'Potential Matrix Time:',F13.1,'s')")
     &        time2_pot-time1_pot
         endif
        end if
      end if
c
c     ----- Add cosmo potential when needed -----
c
cc AJL/Begin/SPIN ECPs
      if(cosmo_on.and.cosmo_phase.eq.2) then
c         call int_1e_ga(ao_bas_han,ao_bas_han,g_fock,
c         t0 = util_wallsec()
         if (.not.ga_duplicate(g_focks(1),g_cosmo,'cosmo'))
     &     call errquit('dft_fockbld: create g_cosmo',0,GA_ERR)
         call ga_zero(g_cosmo)
         call cosmo_fock(ao_bas_han, g_cosmo)
         call ga_add(1d0, g_cosmo, 1d0, g_focks(1), g_focks(1))
         if (ecp_channels.gt.1)
     &     call ga_add(1d0, g_cosmo, 1d0, g_focks(2), g_focks(2))
         if (.not.ga_destroy(g_cosmo))
     &     call errquit('dft_fockbld: destroy g_cosmo',0,GA_ERR)    
c         if (me.eq.0) then
c           write(*,'("cos_chg_pot: ",F8.3)') util_wallsec()-t0
c         endif
      endif

c     == Copy Fock matrix to individual components for normal calcs ==
      if (ecp_channels.eq.1.and.ipol.gt.1) 
     &  call ga_copy(g_focks(1), g_focks(2))

      call ga_sync()
c
cc      Ecore = ga_ddot(g_dens(1), g_fock)
cc AJL: Ecore = ga_ddot(g_dens(1), g_focks(1))
cc      write(6,*)'Ecore(1:T):',ga_ddot(g_dens(1), g_focks(1))

cc    With spin-polarised ECPs Ecore will be spin dependent
cc    See Szabo and Ostlund pg. 215
cc        
      if (ecp_channels.gt.1) then
c     
c       Restore alpha and beta densities to calculate energy
        call ga_dadd(1d0, g_dens(1), -1d0, g_dens(2), g_dens(1))

        Ecore = ga_ddot(g_dens(1), g_focks(1))
c        write(6,*)'Ecore(1:A):',Ecore
        Ecore = Ecore + ga_ddot(g_dens(2), g_focks(2))
c        write(6,*)'Ecore(2:B):',ga_ddot(g_dens(2), g_focks(2))
c
c       And then recombine the densities for now. 
        call ga_dadd(1d0, g_dens(1), 1d0, g_dens(2), g_dens(1))
c
c       call ga_sync()   
      else

        Ecore = ga_ddot(g_dens(1), g_focks(1))
c        write(6,*)'Ecore(1:T):',ga_ddot(g_dens(1), g_focks(1)) 
 
      endif
cc AJL/End
c
c     add Bonacic-Fantucci repulsive term
c
      if (.not.rtdb_get(rtdb, 'dft:fant_d', mt_dbl, 1,
     &   fant_d)) fant_d=-1d0
      if (.not.rtdb_get(rtdb, 'dft:fant_a', mt_dbl, 1,
     &   fant_a)) fant_a=-1d0
      if(fant_a.ne.-1d0.and.fant_d.ne.-1d0)
     A     Ecore=Ecore+dft_fant(geom,natoms,fant_a,fant_d)
      if (odftps) call pstat_off(ps_onel)
c
c     Pre-compute reduced total density matrices over atoms
c 
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit('dft_scf: cannot allocate scr',0, MA_ERR)
      call util_ga_mat_reduce(nbf_ao, natoms, cetobfr, g_dens, 
     &                     ipol, rdens_atom, 'absmax', 
     &                     dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
      if (.not.ma_pop_stack(lscr))
     &   call errquit('dft_scf: cannot pop stack',0, MA_ERR)
c
      if (CDFIT)then
c
c        Fit the electron charge density.
c     
         if (.not.MA_Push_Get(MT_Dbl,nbf_cd,'cd_coef',lcd_coef,
     &        icd_coef))
     &        call errquit('dft_scf: cannot allocate cd_coef',0, MA_ERR)
         if (oprint_time)
     &        call dft_tstamp(' Before call to FITCD.   ')
         if (odftps) call pstat_on(ps_fitcd)
         call case_setflags(.false.)  ! no attenuation for Coulomb
         call dft_fitcd(1,Dbl_MB(icd_coef),dd3cERI,Ecoul,
     &                  g_dens, nTotEl, n_batch, n3c_int,
     &                  ii3cwhat, n3c_dbl, IOLGC,
     .        natoms,
     &                  .false., 0d0, .false.)
         if (odftps) call pstat_off(ps_fitcd)
      endif
c     
      if (oprint_time)
     &     call dft_tstamp(' Before call to GETVCOUL.')
      if (odftps) call pstat_on(ps_vcoul)
      call dft_getvc(Dbl_MB(icd_coef), dd3cERI, Ecoul,
     &               g_tmp, iVcoul_opt, n_batch, 
     &               n3c_int, ii3cwhat, n3c_dbl,
     &               IOLGC,
     &               .false., 1)
      if (odftps) call pstat_off(ps_vcoul)

      if (no2e) call ga_zero(g_tmp)
c     
c     Add V coul to Fock Matrix
c     
cc AJL/Begin/SPIN ECPs
c      call ga_dadd(1d0, g_tmp, 1d0, g_fock, g_fock)
      do i=1,ipol
        call ga_dadd(1d0, g_tmp, 1d0, g_focks(i), g_focks(i))
      enddo
cc AJL/End
c
      if (CDFIT)then
         if (.not.ma_pop_stack(lcd_coef))
     &        call errquit('dft_scf: cannot pop stack',0, MA_ERR)
      endif
c
c     
      if (ipol.gt.1)then
c     
c        Restore alpha and beta densities.
c     
         call ga_dadd(1d0, g_dens(1), -1d0, g_dens(2), g_dens(1))
c     
c        Note that g_dens(1) now contains the alpha density
c        matrix and g_dens(2) contains the beta
c
c        Pre-compute reduced alpha and beta density matrices over atoms
c 
         call dfill(ipol*natoms*natoms, 0.0d0, rdens_atom, 1)
         nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
         if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &      call errquit('dft_scf: cannot allocate scr',0, MA_ERR)
         call util_ga_mat_reduce(nbf_ao, natoms, cetobfr, 
     &                       g_dens, ipol, rdens_atom, 
     &                       'absmax', dbl_mb(iscr), nbf_ao_mxnbf_ce,
     A        .true.)
         if (.not.ma_pop_stack(lscr))
     &      call errquit('dft_scf: cannot pop stack',0, MA_ERR)
      endif
c
c     DIM/QM JEM
c     Calculate DIM potential and energy
c     We calculate the DIM potential at each point from the total density
c     This potential is used to perturb both alpha and beta spins, so only
c     one call to dim_grid_quadv0 is needed
      if (ldimqm .and. dimqm_go) then
        if(.not.ga_duplicate(g_vxc(1), g_vdim, 'dim:1'))
     $    call errquit("dft_fockbld: can't create dim:1", 99, GA_ERR)
        call ga_zero(g_vdim)
        call dim_grid_quadv0(rtdb, g_vdim)
c        call ga_dadd(1.d0, g_vxc(1), 1.d0, g_fock, g_fock)
        Edim = ga_ddot(g_dens(1), g_vdim)
        if(ipol > 1) then
          Edim = Edim + ga_ddot(g_dens(2), g_vdim)
        end if
        call ga_sync()
      end if
c     
c     == Compute the XC potential and energy ==
      g_vxc(1) = g_tmp
      call ga_zero(g_vxc(1))
      rho_n = 0.0d0
cc AJL/Begin
cc Can't remember why I implemented this? To be tested.
cc      if (ipol.gt.1) call ga_zero(g_vxc(2))
      if (ipol.gt.1) call ga_copy(g_vxc(1), g_vxc(2))
cc AJL/End
      if (oprint_time)call dft_tstamp(' Before call to GETVXC.  ')
c
      if (odftps) call pstat_on(ps_xc)
c
c     == Calculate the HF exchange (exact Coulomb), DFT exchange-correlation ==
      call xc_getv(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens, 
     &             g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &             rdens_atom, cetobfr, natoms)
c
      if (odftps) call pstat_off(ps_xc)
c
      if (no2e) then
        call ga_zero(g_vxc(1))
        call ga_zero(g_vxc(2))
      end if
c
c     add relativistic energy contributions
      if (do_zora .or. do_x2c) then
        Erel_sf =0.d0
        if (ipol.gt.1) then
          Erel_sf = ga_ddot(g_dens(1),g_rel_sf(1))
     &             + ga_ddot(g_dens(2),g_rel_sf(2))
        else
          Erel_sf = ga_ddot(g_dens(1),g_rel_sf(1))
        end if
        Ecore = Ecore + Erel_sf 
c
        call ga_dadd(1.d0,g_vxc(1),1.d0,g_rel_sf(1),g_vxc(1))
        if (ipol.gt.1) 
     &   call ga_dadd(1.d0,g_vxc(2),1.d0,g_rel_sf(2),g_vxc(2))
c
      end if  ! do_zora .or. do_x2c
c
c     == if frozen embedding ==
      Efrozemb = 0.d0
      if (frozemb) then

        do i=1,ipol
          Efrozemb = Efrozemb + ga_ddot(g_dens(i),g_frozemb)
          call ga_dadd(1.d0,g_vxc(i),1.d0,g_frozemb,g_vxc(i))
        enddo

        Ecore = Ecore + Efrozemb   ! add in linearized frozen embedding energy to core
c
cc AJL/Begin/FDE
        if (frozemb_fde) then
c Add in contributions for Non-additive potential here.
c This will be for XC:  Vxc(P_tot) - Vxc(P_qm)
c And the same for Ts:  Ts(P_tot)  - Ts(P_qm)
c
c Order:
c -> Call grid_quadv0 equivalent, but with XC set to FDE option
c -> Then recalculate values for P_tot
c -> And then finally subtract P_qm from P_tot values.. et voila.
c
c This is all done in grid_quadv0b, except for the final manipulation of
c the Hamiltonian contributions as we cannot manipulate GAs at this
c lower level (I think. I haven't tried. So subtraction is below)

cc For testing of outputs
c
cc Exc(QM)
c        call dft_frozemb_fde_getvxc(rtdb, Exc_fde, nExc, iVxc_opt,
c     &       g_xcinv_fde, g_dens, g_vxc_fde, IOLGC, 2, g_dens_fde)
cc Exc(FDE)
c        call dft_frozemb_fde_getvxc(rtdb, Exc_fde, nExc, iVxc_opt,
c     &       g_xcinv_fde, g_dens, g_vxc_fde, IOLGC, 3, g_dens_fde)
cc Exc(Non-add)
          call dft_frozemb_fde_getvxc(rtdb, Exc_fde, nExc, iVxc_opt,
     &         g_xcinv_fde, g_dens, g_vxc_fde, IOLGC, 4, g_dens_fde)

c     == Add FDE non-additive Vxc to the original matrix ==      
          do i=1,ipol 
            call ga_dadd(1.d0, g_vxc(i), 1.d0, g_vxc_fde(i), g_vxc(i))
          enddo
c     == Add FDE non-additive Exc in to the Exc ==
          do i=1,nExc 
            Exc(i) = Exc(i) + Exc_fde(i)
          enddo
c
c Repeat for Non-additive Ts; set nExc explicitly to 1 as not
c interested in separating our exchange and correlation!
c
cc Ts(QM)
c        call dft_frozemb_fde_getvxc(rtdb, Exc_fde, 1, iVxc_opt,
c     &       g_xcinv_fde, g_dens, g_vxc_fde, IOLGC, 5, g_dens_fde)
cc Ts(FDE)
c        call dft_frozemb_fde_getvxc(rtdb, Exc_fde, 1, iVxc_opt,
c     &       g_xcinv_fde, g_dens, g_vxc_fde, IOLGC, 6, g_dens_fde)
cc Ts(Non-add)
          call dft_frozemb_fde_getvxc(rtdb, Exc_fde, 1, iVxc_opt, 
     &         g_xcinv_fde, g_dens, g_vxc_fde, IOLGC, 7, g_dens_fde)

c     == Add FDE non-additive VTs to the original matrix ==      
          do i=1,ipol 
            call ga_dadd(1.d0, g_vxc(i), 1.d0, g_vxc_fde(i), g_vxc(i))
          enddo
c     == Add FDE non-additive Ts energy in to the 1e energy ==
          Ecore = Ecore + Exc_fde(1)
c
        end if ! frozemb_fde 
cc AJL/End
c
      end if
c
c     == Add exact pot contributions ==
      if (exact_pot.eq.1) then
        call xc_exact_pot(Ecoul_aux, g_dens, g_vxc, g_rho_exact,
     &                    dif_lamda, g_tmp_exact, g_addit_exact)
      end if
c
c  Jorge Garza/Begin
      if (ipol.eq.2.and.noc(2).eq.0) call ga_zero(g_vxc(2))
      totsic=0.0d00
      if (test_sic.eq.2.or.test_sic.eq.4) then
        if (iter.ge.iter_ini) then
          if (test_sic.eq.4) then
            call localize_sic(g_movecs, evals)
          end if
          call xc_sic(rtdb, nExc,iVxc_opt, g_dens, 
     &                g_vxc,  g_wght, g_xyz, g_nq, wght_GA,
     &                rdens_atom, cetobfr, natoms,
     &                g_movecs, totsic, i_degen, n_levels)
          Exc(1) = Exc(1) + totsic
        end if
      end if
c  Jorge Garza/End
c
cc AJL/Begin/SPIN ECPs
cc     == Copy Fock matrix to individual components ==
c      call ga_copy(g_fock, g_focks(1))
c      if (ipol.gt.1) 
c     &     call ga_copy(g_fock, g_focks(2))
cc AJL/End
c
c     == Add exchange-correlation to the Fock matrix ==       
      call ga_dadd(1.d0, g_vxc(1), 1.d0, g_focks(1), g_focks(1))
      if (ipol.gt.1) 
     &     call ga_dadd(1.d0, g_vxc(2), 1.d0, g_focks(2), g_focks(2))

c jbecca START - damnit, quit taking this out
C
C     Compute -D.E contribution due to external applied electric field
C     (if dft:efield is stored in rtdb).  Note that we can potentially
C     apply different fields to the alpha and beta spins this way.
C
C
      if (rtdb_get (rtdb, "dft:efield", mt_dbl, 3, efield)) then

         if (.not. rtdb_get (rtdb, "dft:efield_scaled",
     $        mt_log, 1, lescale))
     $        lescale = .false.

         call ga_zero (g_tmp)  !careful, this previously stored g_vxc
         call dft_fock_efield (rtdb, efield, g_tmp, lescale)  !g_tmp now has -D.E
c Daniel (10-9-12): This line causes problems.
c         call ga_norm_infinity (g_tmp, enorm)

         if (ipol.eq.1) then
            call ga_dadd (1d0, g_focks(1), 1d0, g_tmp, g_focks(1))
            eelec = ga_ddot (g_dens(1), g_tmp)
         elseif (ipol.eq.2) then
            call ga_dadd (1d0, g_focks(1), 1d0, g_tmp, g_focks(1))
            call ga_dadd (1d0, g_focks(2), 1d0, g_tmp, g_focks(2))
            eelec = ga_ddot (g_dens(1), g_tmp) +
     $           ga_ddot (g_dens(2), g_tmp)
         else
            call errquit ("dft_fockbld: bad ipol", ipol, 0)
         endif

c$$$         if (me.eq.0) write (luout,"(1x,a,3es12.4,a,1es12.4,a,1es12.4)")
c$$$     $        "[dft_fockbld] External E-field:", efield,
c$$$     $        "; interaction energy =", eelec, "; |-D.E| =", enorm

         if (me.eq.0) then
c            write (luout,"(1x,a,3es12.4,a,a,l,a,a,1es12.4,a,a,1es12.4)")
c     $           "Applied E-field = (", efield, " ), ",
c     $           "scaled = ", lescale, ", ",
c     $           "||D.E|| =", enorm, ", ",
c     $           "interaction energy =", eelec
            write (luout,"(1x,a,3es12.4,a,a,l,a,a,1es12.4)")
     $           "Applied E-field = (", efield, " ), ",
     $           "scaled = ", lescale, ", ",
     $           "interaction energy =", eelec
c$$$     $           "External E-field: ",
c$$$     $           "scaled = ", lescale,
c$$$     $           "; E = (", efield, " )"
         endif
      endif
c  jbecca END

c
c     DIM/QM JEM
c     Add DIM potential to the Fock matrix
      if (ldimqm .and. dimqm_go) then
c  jbecca START: scaling potential here to test some stuff.
c              originally 1.d0
        call ga_dadd(1.d0, g_vdim, 1.d0, g_focks(1), g_focks(1))
c       DIM potential is from total density, so same for both alpha and beta
        if (ipol.gt.1)
     &    call ga_dadd(1.d0, g_vdim, 1.d0, g_focks(2), g_focks(2))
        if (.not. ga_destroy(g_vdim))
     &    call errquit("destroy vdim failed", 9, GA_ERR)
      end if
      if (oprint_time)
     &     call dft_tstamp(' End of parallel region. ')
c     
      if (odftps) call pstat_off(ps_fockb)
c
cc AJL/Begin/SPIN ECPs (Debug)
c      if(util_print('Intermediate Fock Matrix', print_debug)) then
c        if (ga_nodeid() .eq. 0) then
c          write(6,*)
c          call util_print_centered(6,
c     $         'Intermediate Fock Matrix',40,.true.)
c          write(6,*)
c          call util_flush(6)
c        end if
c        call ga_print(g_focks(1))
c        if (ipol.eq.2) call ga_print(g_focks(2))
cc AJL/End
c
c MN solvation models -->
c
c     ----- Add SMD-CDS energy contribution when needed -----
c
c      if(cosmo_on.and.cosmo_phase.eq.2.and.do_cosmo_smd) then
c        Ecore=Ecore+gcds
c      endif
c
c <-- MN solvation models
c
      return
      end
