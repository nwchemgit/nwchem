      subroutine dft_fockbld(rtdb,natoms,ntotel,
     ,     g_fock,g_dens,g_vxc,g_movecs,g_xcinv,
     ,     g_tmp,g_nq,g_wght,g_xyz,
     ,     iwhat_max,ivcoul_Opt,nexc,ivxc_opt,wght_GA,
     ,     n_batch,n3c_int,n3c_dbl,n_semi_bufs,fd,IOLGC,
     ,     evals, rdens_atom,dd3ceri,ii3cwhat,
     ,     cetobfr,
     ,     iter, ecore, ecoul,exc,noc, rho_n,
     ,     test_sic,exact_pot,ecoul_aux,dif_lamda,totsic,
     ,     g_rho_exact,g_tmp_exact, g_addit_exact,
     ,     oprint_time)
c     
c     $Id: dft_fockbld.F,v 1.4 2000-04-24 19:56:16 edo Exp $
c     
      implicit none
#include "cdft.fh"
#include "global.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      integer g_fock,g_dens(2),g_vxc(4),g_movecs(2) ! [handle in/values out]
      integer g_xcinv
      integer g_tmp,g_nq,g_wght,g_xyz
      integer g_rho_exact,g_tmp_exact, g_addit_exact
      integer rtdb ! [in]
      integer natoms ! [in]
      logical wght_GA
      integer iwhat_max
      double precision evals (*)
      double precision rdens_atom(*)
      double precision dd3ceri(*)
      integer ii3cwhat(*)
      integer cetobfr(*)
      logical oprint_time ![in]
      integer iter ! [in]
      integer noc(2)
      double precision ecore, ecoul,exc(2) ! [out]
      integer ntotel,n_batch,n3c_int,n3c_dbl,n_semi_bufs,fd
      logical IOLGC
      integer test_sic ! [in]
      integer exact_pot
      double precision ecoul_aux(2),dif_lamda,totsic
      double precision rho_n ![out]
      integer ivcoul_Opt,nexc,ivxc_opt
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer l_degen, i_degen
      integer me,nscr,lscr,iscr,lcd_coef,icd_coef
c
      me=ga_nodeid()
c

      Ecoul  = 0d0
      Exc(1)= 0d0
      Exc(2) = 0d0
c Jorge Garza/Begin
      if (test_sic.ne.0) then
        if (iter.ge.3) then
           call xc_sicdeg(l_degen, i_degen, evals, noc)
        end if
      end if
c Jorge Garza/End
c     
c     Accumulate core hamiltonian into Fock matrix; compute core energy
c     
      if (oprint_time)
     &     call dft_tstamp(' calling 1-e int ')
      call ga_zero(g_fock)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'kinetic', oskel)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'potential', oskel)
c
      Ecore = ga_ddot(g_dens(1), g_fock)
c
c     Pre-compute reduced total density matrices over atoms
c 
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit('dft_scf: cannot allocate scr',0)
      call util_ga_mat_reduce(nbf_ao, natoms, cetobfr, g_dens, 
     &                        ipol, rdens_atom, 'rms', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce)
      if (.not.ma_pop_stack(lscr))
     &   call errquit('dft_scf: cannot pop stack',0)
#if defined(DFT_TIMER)
      call dft_second(tp3)
#endif
c
      if (CDFIT)then
c     
c        Fit the electron charge density.
c     
         if (.not.MA_Push_Get(MT_Dbl,nbf_cd,'cd_coef',lcd_coef,
     &        icd_coef))
     &        call errquit('dft_scf: cannot allocate cd_coef',0)
         if (oprint_time)
     &        call dft_tstamp(' Before call to FITCD.   ')
         call dft_fitcd(1,Dbl_MB(icd_coef),dd3cERI,Ecoul,
     &                  g_dens, nTotEl, n_batch, n3c_int,
     &                  ii3cwhat, n3c_dbl, iwhat_max, 
     &                  n_semi_bufs, fd, IOLGC, 
     .        natoms,
     &                  .false., 0d0, .false.)
      endif
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(6) = '   FitCD  '
      call tcllct(t_pgm(1,6),tp4,tp3)
      call dft_second(tp3)
#endif
      if (oprint_time)
     &     call dft_tstamp(' Before call to GETVCOUL.')
      call dft_getvc(Dbl_MB(icd_coef), dd3cERI, Ecoul,
     &               g_tmp, iVcoul_opt, n_batch, 
     &               n3c_int, ii3cwhat, n3c_dbl, iwhat_max,
     &               n_semi_bufs, fd, IOLGC,
     &               .false., 1)
c     
c     Add V coul to Fock Matrix
c     
      call ga_dadd(1d0, g_tmp, 1d0, g_fock, g_fock)
      if (CDFIT)then
         if (.not.ma_pop_stack(lcd_coef))
     &        call errquit('dft_scf: cannot pop stack',0)
      endif
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(7) = '   GtVcoul'
      call tcllct(t_pgm(1,7),tp4,tp3)
      call dft_second(tp3)
#endif
      if (ipol.gt.1)then
c     
c        Restore alpha and beta densities.
c     
         call ga_dadd(1d0, g_dens(1), 1d0, g_dens(2), g_dens(1))
c     
c        Note that g_dens(1) now contains the alpha density
c        matrix and g_dens(2) contains the beta
c
c        Pre-compute reduced alpha and beta density matrices over atoms
c 
         call dfill(ipol*natoms*natoms, 0.0d0, rdens_atom, 1)
         nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
         if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &      call errquit('dft_scf: cannot allocate scr',0)
         call util_ga_mat_reduce(nbf_ao, natoms, cetobfr, 
     &                           g_dens, ipol, rdens_atom, 
     &                           'rms', dbl_mb(iscr), nbf_ao_mxnbf_ce)
         if (.not.ma_pop_stack(lscr))
     &      call errquit('dft_scf: cannot pop stack',0)
      endif
c     
c     Compute the XC potential and energy.
c     
      g_vxc(1) = g_tmp
      call ga_zero(g_vxc(1))
      rho_n = 0.0d0
      if (ipol.eq.2)call ga_zero(g_vxc(2))
      if (oprint_time)call dft_tstamp(' Before call to GETVXC.  ')
c
      call xc_getv(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens, 
     &             g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &             rdens_atom, cetobfr, natoms)
c
      if (exact_pot.eq.1) then
        call xc_exact_pot(Ecoul_aux, g_dens, g_vxc, noc, g_rho_exact,
     &                    dif_lamda, g_tmp_exact, g_addit_exact)
      end if
c
c  Jorge Garza/Begin
c       
      if (ipol.eq.2.and.noc(2).eq.0) call ga_zero(g_vxc(2))
      totsic=0.0d00
      if (test_sic.eq.2.or.test_sic.eq.4) then
        if (iter.ge.3) then
          if (test_sic.eq.4) then
            call localize_sic(g_movecs, noc)
          end if
          call xc_sic(rtdb, nExc,iVxc_opt,g_xcinv, g_dens, 
     &                g_vxc,  g_wght, g_xyz, g_nq, wght_GA,
     &                rdens_atom, cetobfr, natoms,
     &                g_movecs, noc, totsic, i_degen)
          Exc(1) = Exc(1) + totsic
        end if
      end if
c
c  Jorge Garza/End
c
      call ga_dadd(1.d0, g_vxc(1), 1.d0, g_fock, g_fock)
      if (ipol.gt.1)then
c     
c     Put (-alpha+beta) in g_vxc(2)
c     
         call ga_dadd(-1.d0, g_vxc(1), 1.d0, g_vxc(2), g_vxc(2))
      endif
c     
      if (oprint_time)
     &     call dft_tstamp(' End of parallel region. ')
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(8) = '   GetVXC '
      call tcllct(t_pgm(1,8),tp4,tp3)
#endif

      return
      end
