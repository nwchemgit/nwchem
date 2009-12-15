      subroutine dft_fockbld(rtdb,natoms,ntotel,
     ,     g_fock,g_focks,g_dens,g_vxc,g_movecs,g_xcinv,
     ,     g_tmp,g_nq,g_wght,g_xyz,
     ,     iwhat_max,ivcoul_Opt,nexc,ivxc_opt,wght_GA,
     ,     n_batch,n3c_int,n3c_dbl,n_semi_bufs,fd,IOLGC,
     ,     evals, rdens_atom,dd3ceri,ii3cwhat,
     ,     cetobfr, iter, ecore, ecoul,exc,noc, rho_n,ldmix,
     ,     test_sic,exact_pot,ecoul_aux,dif_lamda,totsic,
     ,     i_degen,g_rho_exact,g_tmp_exact, g_addit_exact,
     ,     n_levels, act_levels, oprint_time,
     ,     g_zora_sf )
c     
c     $Id: dft_fockbld.F,v 1.36 2008-11-11 00:40:33 niri Exp $
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
#include "case.fh"
c
      integer g_fock,g_dens(2),g_vxc(4),g_movecs(2) ! [handle in/values out]
      integer g_focks(2)
      integer g_xcinv
      integer g_tmp,g_nq,g_wght,g_xyz
      integer g_rho_exact,g_tmp_exact, g_addit_exact
      integer rtdb ! [in]
      integer natoms ! [in]
      logical wght_GA
      integer iwhat_max
c     double precision evals (*)
      integer evals(2)
      double precision rdens_atom(*)
      double precision dd3ceri(*)
      integer ii3cwhat(*)
      integer cetobfr(*)
      logical oprint_time ![in]
      integer iter ! [in]
      integer noc(2), n_levels(2), act_levels
      double precision ecore, ecoul,exc(2) ! [out]
      integer ntotel,n_batch,n3c_int,n3c_dbl,n_semi_bufs,fd
      logical IOLGC, ldmix
      integer test_sic ! [in]
      integer exact_pot, iter_ini
      integer i_degen(2)  ! [in]
      double precision ecoul_aux(2),dif_lamda,totsic
      double precision rho_n ![out]
      double precision fant_a,fant_d,dft_fant
      external dft_fant
      integer ivcoul_Opt,nexc,ivxc_opt
c
      integer me,nscr,lscr,iscr,lcd_coef,icd_coef
c
      double precision Ebq
c
czora...start
      integer g_zora_sf(*)
      double precision Ezora_sf
czora...end
c
c     timings
c
      double precision time1_ke,time2_ke
      double precision time1_pot,time2_pot
c
      me=ga_nodeid()
      if (odftps) call pstat_on(ps_fockb)
c
      Ecoul  = 0d0
      Exc(1)= 0d0
      Exc(2) = 0d0

c Jorge Garza/Begin
      if (ldmix) then
        iter_ini = 1
      else
        iter_ini = 3
      end if
      if (test_sic.eq.2.or.test_sic.eq.4) then
        if (iter.ge.iter_ini) then
           call xc_sicdeg(i_degen, evals, noc, n_levels,act_levels)
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
c     if requested separate calculation of bq-electron energy (MV)
c     ----------------------------------------------------------- 
      if (rtdb_get(rtdb, 'dft:bq_energy', mt_dbl, 1,Ebq)) then
        Ecore = ga_ddot(g_dens(1), g_fock)
        call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 
     &                 'bq_pot', oskel)
        Ebq = ga_ddot(g_dens(1), g_fock)-Ecore
        call int_1e_ga(ao_bas_han, ao_bas_han, g_fock,
     &                 'potential0', oskel)
        if (.not. rtdb_put(rtdb,'dft:bq_energy',mt_dbl,1,Ebq)) 
     $     call errquit('dft: writing bq_energy failed', 
     $     0, RTDB_ERR)
      else
        if(util_print('dft timings', print_high))
     &        time1_pot=util_cpusec()   ! start time
        call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 
     &                 'potential', oskel)
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
      if(cosmo_on.and.cosmo_phase.eq.2) then
         call int_1e_ga(ao_bas_han,ao_bas_han,g_fock,
     &                  'cos_chg_pot',.false.)
      endif
      call ga_sync()
c
      Ecore = ga_ddot(g_dens(1), g_fock)
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
     &                        ipol, rdens_atom, 'rms', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
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
     &                  ii3cwhat, n3c_dbl, iwhat_max, 
     &                  n_semi_bufs, fd, IOLGC, 
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
     &               n3c_int, ii3cwhat, n3c_dbl, iwhat_max,
     &               n_semi_bufs, fd, IOLGC,
     &               .false., 1)
      if (odftps) call pstat_off(ps_vcoul)
c     
c     Add V coul to Fock Matrix
c     
      call ga_dadd(1d0, g_tmp, 1d0, g_fock, g_fock)
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
     &                           g_dens, ipol, rdens_atom, 
     &                           'rms', dbl_mb(iscr), nbf_ao_mxnbf_ce,
     A        .true.)
         if (.not.ma_pop_stack(lscr))
     &      call errquit('dft_scf: cannot pop stack',0, MA_ERR)
      endif
c     
c     == Compute the XC potential and energy ==
      g_vxc(1) = g_tmp
      call ga_zero(g_vxc(1))
      rho_n = 0.0d0
      if (ipol.eq.2)call ga_zero(g_vxc(2))
      if (oprint_time)call dft_tstamp(' Before call to GETVXC.  ')
c
      if (odftps) call pstat_on(ps_xc)
c
czora...We have the correction, turn off the correction
czora...so_term is set via a common block (zora.fh)
czora...so_term = 0 (spin-free and the matrices are symmetric)
      do_zora_get_correction = .false.
      so_term = 0
c
c     == Calculate the HF exchange (exact Coulomb), DFT exchange-correlation ==
      call xc_getv(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens, 
     &             g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &             rdens_atom, cetobfr, natoms) 
c
      if (odftps) call pstat_off(ps_xc)
c
c     == Add ZORA contributions ==
      if (do_zora) then
        Ezora_sf =0.d0
        if (ipol.gt.1) then
          Ezora_sf = ga_ddot(g_dens(1),g_zora_sf(1))
     &      + ga_ddot(g_dens(2),g_zora_sf(2))
        else
          Ezora_sf = ga_ddot(g_dens(1),g_zora_sf(1))
        end if
        Ecore = Ecore + Ezora_sf 
c
        call ga_dadd(1.d0,g_vxc(1),1.d0,g_zora_sf(1),g_vxc(1))
        if (ipol.gt.1) 
     &   call ga_dadd(1.d0,g_vxc(2),1.d0,g_zora_sf(2),g_vxc(2))
c
      end if  ! do_zora
c
c     == Add exact pot contributions ==
      if (exact_pot.eq.1) then
        call xc_exact_pot(Ecoul_aux, g_dens, g_vxc, noc, g_rho_exact,
     &                    dif_lamda, g_tmp_exact, g_addit_exact)
      end if
c
c  Jorge Garza/Begin
      if (ipol.eq.2.and.noc(2).eq.0) call ga_zero(g_vxc(2))
      totsic=0.0d00
      if (test_sic.eq.2.or.test_sic.eq.4) then
        if (iter.ge.iter_ini) then
          if (test_sic.eq.4) then
            call localize_sic(g_movecs, noc, evals)
          end if
          call xc_sic(rtdb, nExc,iVxc_opt, g_dens, 
     &                g_vxc,  g_wght, g_xyz, g_nq, wght_GA,
     &                rdens_atom, cetobfr, natoms,
     &                g_movecs, noc, totsic, i_degen, n_levels)
          Exc(1) = Exc(1) + totsic
        end if
      end if
c  Jorge Garza/End
c
c     == Copy Fock matrix to individual components ==
      call ga_copy(g_fock, g_focks(1))
      if (ipol.gt.1) 
     &     call ga_copy(g_fock, g_focks(2))
c
c     == Add exchange-correlation to the Fock matrix ==       
      call ga_dadd(1.d0, g_vxc(1), 1.d0, g_focks(1), g_focks(1))
      if (ipol.gt.1) 
     &     call ga_dadd(1.d0, g_vxc(2), 1.d0, g_focks(2), g_focks(2))
c     
      if (oprint_time)
     &     call dft_tstamp(' End of parallel region. ')
c     
      if (odftps) call pstat_off(ps_fockb)
c
      return
      end
