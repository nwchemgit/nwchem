C>
C> \brief Calculate Fukui indeces
C>
C> Calculate Fukui indeces for open shell wavefunctions [1,2].
C>
C> [1] E. Chamorro, P. P&eacute;rez,
C>     <i>"Condensed-to-atoms electronic Fukui functions within the
C>     framework of spin-polarized density-functional theory"</i>,
C>     J. Chem. Phys. <b>123</b>, 114107 (2005), DOI:
C>     <a href="https://doi.org/10.1063/1.2033689">
C>     10.1063/1.2033689</a>.
C>
C> [2] R.G. Parr, W. Yang,
C>     <i>"Density functional approach to the frontier-electron theory
C>     of chemical reactivity"</i>,
C>     J. Am. Chem. Soc. <b>106</b>, 4049-4050 (1984), DOI:
C>     <a href="https://doi.org/10.1021/ja00326a036">
C>     10.1021/ja00326a036</a>.
C>
      Subroutine fukui(g_movecs, k_eval, tol2e, rtdb,
     &                 nExc, iVxc_opt, g_xcinv, IOLGC, g_wght,
     &                 g_xyz, g_nq, wght_GA, rho_n, irdens_atom,
     &                 icetobfr, natoms)
c     
c     $Id$
c     
      implicit none
#include "errquit.fh"
      integer g_dens_HOMO(2), g_dens_LUMO(2), g_orb, g_dens_ss, ik,
     &        g_movecs(2), isp, k_eval(2), g_s,
     &        me, l_temp_vec, i_temp_vec,
     &        g_dens(2), noc_aux_1, noc_aux_2, noc_test
      integer irdens_atom, icetobfr
c
      integer nExc
      integer iVxc_opt
      integer g_xcinv, g_vxc(4), g_wght, g_xyz,g_nq
      integer natoms
      logical IOLGC, wght_GA
      integer rtdb
c
      double precision eig_lumo(2), eig_homo(2), jfac(4), kfac(4),
     &                 mu_n_mas, mu_n_men, mu_n_cer, tol2e,
     &                 mu_s_mas, mu_s_men, mu_s_cer, 
     &                 int_HaHa, int_HbHb, int_LaLa, int_LbLb,
     &                 int_HaLb, int_LaHb, Exc(2), ecoul, rho_n,
     &                 int_vxc_H(2), int_vxc_L(2), Exc_zero,
     &                 Exc_pert, diff_Exc, ion_pot, ele_afi,
     &                 high_mult, low_mult, e_orbital, e_coul,
     &                 e_xc
c
c #include "stdio.fh"
c #include "mafdecls.fh"
c #include "global.fh"
c #include "util.fh"
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "util.fh"
#include "sym.fh"
c
#include "stdio.fh"
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical oprint_fukui
      me = ga_nodeid()
      oprint_fukui = util_print('Fukui information', print_high)
      g_orb = ga_create_atom_blocked(geom, AO_bas_han, 'ga_orb')
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'temp vec',
     &   l_temp_vec, i_temp_vec))
     &   call errquit('fukui: cannot allocate temp vec',0, MA_ERR)
      do isp = 1, 2
         g_dens_HOMO(isp) = ga_create_atom_blocked(geom,
     &                      AO_bas_han, 'ga_dens_orb_homo')
         g_dens_LUMO(isp) = ga_create_atom_blocked(geom,
     &                      AO_bas_han, 'ga_dens_orb_lumo')
         call ga_zero(g_dens_HOMO(isp))
         call ga_zero(g_dens_LUMO(isp))
      end do
c
      g_dens_ss = ga_create_atom_blocked(geom, AO_bas_han,
     &            'ga_dens_orb_ss')
c
      do isp = 1, ipol
         if (noc(isp).eq.0) then
           noc_aux_1 = 1
           noc_aux_2 = 1
         else
           noc_aux_1 = noc(isp)
           noc_aux_2 = noc(isp) + 1
         end if
C        do ik = noc(isp), (noc(isp)+1)
         do ik = noc_aux_1, noc_aux_2
            call ga_zero(g_orb)
            call ga_get(g_movecs(isp), 1, nbf_ao, ik, ik,
     &                  Dbl_MB(i_temp_vec), nbf_ao)
            call ga_put(g_orb, 1, nbf_ao, ik, ik, 
     &                  Dbl_MB(i_temp_vec), nbf_ao)
            if (ik.eq.noc(isp)) then
               eig_homo(isp) = dbl_mb(k_eval(isp) + 
     &                         noc(isp) - 1)
               call ga_dgemm('N', 'T', nbf_ao, nbf_ao, nbf_ao,
     &                       1.0d00, g_orb, g_orb, 0.d00, 
     &                       g_dens_HOMO(isp))
            else
               eig_lumo(isp) = dbl_mb(k_eval(isp) + noc(isp))
               call ga_dgemm('N', 'T', nbf_ao, nbf_ao, nbf_ao,
     &                       1.0d00, g_orb, g_orb, 0.d00,
     &                       g_dens_LUMO(isp))
            end if
         end do
      end do
c
      if (.not.ma_pop_stack(l_temp_vec))
     &   call errquit('fukui: cannot pop stack',0, MA_ERR)
c      
      if (ipol.eq.1) then
         call ga_copy(g_dens_HOMO(1), g_dens_HOMO(2))
         call ga_copy(g_dens_LUMO(1), g_dens_LUMO(2))
         eig_homo(2) = eig_homo(1)
         eig_lumo(2) = eig_lumo(1)
      end if
c
      if (me.eq.0.and.oprint_fukui)
     &    call dft_header(' Reactivity Parameters ')
c
      mu_n_mas = 0.5d00*(eig_lumo(1) + eig_lumo(2))
      mu_n_men = 0.5d00*(eig_homo(1) + eig_homo(2))
      mu_n_cer = 0.5*(mu_n_mas + mu_n_men)
      mu_s_mas = 0.5d00*(eig_lumo(1) - eig_homo(2))
      mu_s_men = 0.5d00*(eig_homo(1) - eig_lumo(2))
      mu_s_cer = 0.5*(mu_s_mas + mu_s_men)
c
      g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      call ga_zero(g_s)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_s,'overlap',.false.)
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_LUMO(2), 0.5d00,
     &            g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header('  Condensed Fukui function [fnn(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s)
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(2), 0.5d00,
     &             g_dens_HOMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fnn(-)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s)
c
      call ga_zero(g_dens_ss)
      call ga_add(-0.5d00, g_dens_LUMO(2), 0.5d00,
     &             g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fsn(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s)
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(2), 0.5d00,
     &             g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fss(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s)
c
      call ga_zero(g_dens_ss)
      call ga_add(-0.5d00, g_dens_HOMO(2), 0.5d00,
     &             g_dens_LUMO(1), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fns(+)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s)
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(1), 0.5d00,
     &             g_dens_LUMO(2), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fss(-)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s)
c
      call ga_zero(g_dens_ss)
      call ga_add(0.5d00, g_dens_HOMO(1), -0.5d00,
     &             g_dens_LUMO(2), g_dens_ss)
      if (me.eq.0.and.oprint_fukui)
     &   call dft_header(' Condensed Fukui function [fns(-)]')
      call mull_pop_fuk(geom, AO_bas_han, g_dens_ss, g_s)
c
      mu_n_men = 0.5d00*(eig_homo(1) + eig_homo(2))
      mu_n_cer = 0.5*(mu_n_mas + mu_n_men)
      mu_s_mas = 0.5d00*(eig_lumo(1) - eig_homo(2))
      mu_s_men = 0.5d00*(eig_homo(1) - eig_lumo(2))
      mu_s_cer = 0.5*(mu_s_mas + mu_s_men)
      if (me.eq.0.and.oprint_fukui)then
         write(LuOut,*) ' ------------------------------------'
         write(LuOut,*) '    mu_n(+)    mu_n(-)    mu_n(0)'
         write(LuOut,'(3f11.4)')  mu_n_mas, mu_n_men, mu_n_cer
         write(LuOut,*) ' ------------------------------------'
         write(LuOut,*) '    mu_s(+)    mu_s(-)    mu_s(0)'
         write(LuOut,'(3f11.4)')  mu_s_mas, mu_s_men,mu_s_cer
         write(LuOut,*) ' ------------------------------------'
      endif
c
cc    Evaluating Coulomb integrals for HOMO, LUMO and differences
c
      kfac(1) = 0.d00
      jfac(1) = 1.0d0
      jfac(2) = 1.0d0
      kfac(2) = 0.0d0
      call ga_zero(g_orb)
      call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens_HOMO(1), g_orb, .false.)
      int_HaHa = ga_ddot(g_dens_HOMO(1), g_orb)
c
      int_HaLb = ga_ddot(g_dens_LUMO(2), g_orb)
c
      call ga_zero(g_orb)
      call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens_HOMO(2), g_orb, .false.)
      int_HbHb = ga_ddot(g_dens_HOMO(2), g_orb)
c
      int_LaHb = ga_ddot(g_dens_LUMO(1), g_orb)
c
      call ga_zero(g_orb)
      call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens_LUMO(1), g_orb, .false.)
      int_LaLa = ga_ddot(g_dens_LUMO(1), g_orb)
c
      call ga_zero(g_orb)
      call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens_LUMO(2), g_orb, .false.)
      int_LbLb = ga_ddot(g_dens_LUMO(2), g_orb)
c
c
cc     Evaluating exchange-correlation integrals
c
      g_dens(1) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                       'ga_dens(1)')
      g_vxc(1) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                       'g_vxc(1)')
      call ga_zero(g_dens(1))
      call ga_zero(g_vxc(1))
      if (ipol.eq.2) then
        g_dens(2) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                       'ga_dens(2)')
        g_vxc(2) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                       'g_vxc(2)')
        call ga_zero(g_dens(2))
        call ga_zero(g_vxc(2))
      endif
      do isp=1,ipol
        call ga_dgemm('N', 'T', nbf_ao, nbf_ao, 
     &                noc(isp), 2d0/dble(ipol), g_movecs(isp),
     &                g_movecs(isp), 0.0d00, g_dens(isp))
      enddo
      Exc(1) = 0.0d00
      Exc(2) = 0.0d00
      Ecoul = 0.0d00
      call xc_pot(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens,
     &             g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &             dbl_mb(irdens_atom), 
     &             int_mb(icetobfr), natoms)
      Exc_zero = Exc(1)
      call ga_zero(g_orb)
      call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens_HOMO(1), g_orb, .false.)
      do isp = 1, ipol
        int_vxc_H(isp) = ga_ddot(g_dens_HOMO(isp), g_vxc(isp))
        int_vxc_L(isp) = ga_ddot(g_dens_LUMO(isp), g_vxc(isp))
      end do
      if(ipol .eq.1) then
        int_vxc_H(2) = int_vxc_H(1)
        int_vxc_L(2) = int_vxc_L(1)
      endif
c
cc    Approximating the ionization potential
c
      noc_aux_1 = noc(1)
      noc_aux_2 = noc(2)
C    Think about it.  How can you excite only one electron and still be RKS?
C    Code produces rubbish with RKS - compare to UKS/ROKS and see
C    This code assumes once closed shell, then always closed shell
      if(ipol .eq.1) then
        noc_aux_2 = noc_aux_1
        call errquit('fukui: UKS or ROKS required. No RKS', 0, 0)
      end if

      noc_test = noc_aux_1 - noc_aux_2
      
      if (noc_test.eq.0) then
        noc(2) = noc_aux_2 - 1
        e_orbital = -eig_homo(2)
        e_coul = 0.5d00*int_HbHb
        e_xc = int_vxc_H(2)
      else
        noc(1) = noc_aux_1 - 1
        e_orbital = -eig_homo(1)
        e_coul = 0.5d00*int_HaHa
        e_xc = int_vxc_H(1)
      end if
      ion_pot = e_orbital + e_coul + e_xc
      if (.not.rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc)) 
     &   call errquit('fukui: rtdb_put of noc failed', 0, RTDB_ERR)
      do isp=1,ipol
        call ga_zero(g_dens(isp))
        call ga_dgemm('N', 'T', nbf_ao, nbf_ao,
     &                noc(isp), 2d0/dble(ipol), g_movecs(isp),
     &                g_movecs(isp), 0.0d00, g_dens(isp))
      enddo
     
      Exc(1) = 0.0d00
      Exc(2) = 0.0d00
      Ecoul = 0.0d00
      do isp = 1, ipol
        call ga_zero(g_vxc(isp))
      end do
      call xc_pot(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens,
     &             g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &             dbl_mb(irdens_atom),
     &             int_mb(icetobfr), natoms)
      Exc_pert = Exc(1)
      diff_Exc = Exc_pert - Exc_zero
      ion_pot = ion_pot + diff_Exc
      if (me.eq.0.and.oprint_fukui)then
         write(LuOut,'("  Alpha = ",i3," ; Beta = ",i3)')
     &                   noc(1), noc(2)
         write(LuOut,*) ' Contributions in atomic units:'
         write(LuOut,'("      Orbital Energy = ",f10.4)') e_orbital
         write(LuOut,'("      Coulomb Integral = ",f10.4)') e_coul
         write(LuOut,'("      XC Integral = ",f10.4)') e_xc
         write(LuOut,'("      XC Diff. Energy = ",f10.4)') diff_Exc
         write(LuOut,'("      Ionization potential = ",f7.4," a.u.")') 
     &                                           ion_pot
         write(LuOut,'("                           = ",f7.2," eV")') 
     &                                           ion_pot*27.211
      end if
c
cc    Approximating the electron affinity
c
      if (noc_test.eq.0) then
        noc(1) = noc_aux_1 + 1
        noc(2) = noc_aux_2
        e_orbital = -eig_lumo(1)
        e_coul = - 0.5d00*int_LaLa
        e_xc = int_vxc_L(1)
      else
        noc(1) = noc_aux_1
        noc(2) = noc_aux_2 + 1
        e_orbital = -eig_lumo(2)
        e_coul = - 0.5d00*int_LbLb
        e_xc = int_vxc_L(2)
      end if
      ele_afi = e_orbital + e_coul + e_xc
      if (.not.rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('fukui: rtdb_put of noc failed', 0, RTDB_ERR)
      do isp=1,ipol
        call ga_zero(g_dens(isp))
        call ga_dgemm('N', 'T', nbf_ao, nbf_ao,
     &                noc(isp), 2d0/dble(ipol), g_movecs(isp),
     &                g_movecs(isp), 0.0d00, g_dens(isp))
      enddo

      Exc(1) = 0.0d00
      Exc(2) = 0.0d00
      Ecoul = 0.0d00
      do isp = 1, ipol
        call ga_zero(g_vxc(isp))
      end do
      call xc_pot(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens,
     &             g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &             dbl_mb(irdens_atom),
     &             int_mb(icetobfr), natoms)
      Exc_pert = Exc(1)
      diff_Exc = Exc_pert - Exc_zero
      ele_afi = ele_afi - diff_Exc
      if (me.eq.0.and.oprint_fukui)then
         write(LuOut,*) ' ------------------------------------'
         write(LuOut,'("  Alpha = ",i3," ; Beta = ",i3)')
     &                   noc(1), noc(2)
         write(LuOut,*) ' Contributions in atomic units:'
         write(LuOut,'("      Orbital Energy = ",f10.4)') e_orbital
         write(LuOut,'("      Coulomb Integral = ",f10.4)') e_coul
         write(LuOut,'("      XC Integral = ",f10.4)') e_xc
         write(LuOut,'("      XC Diff. Energy = ",f10.4)') diff_Exc
         write(LuOut,'("      Electron Affinity = ",f7.4," a.u.")')
     &                                        ele_afi
         write(LuOut,'("                        = ",f7.2," eV")')
     &                                        ele_afi*27.211
      end if
      if (me.eq.0.and.oprint_fukui)then
         write(LuOut,*) ' ------------------------------------'
         write(LuOut,'("  Electronegativity (I+A)/2 = ",f7.2," eV")') 
     &                        0.5d00*(ion_pot + ele_afi)*27.211
         write(LuOut,'("  Hardness (I-A) = ",f7.2," eV")') 
     &                        (ion_pot - ele_afi)*27.211
      end if
c
cc    Energy difference for high multiplicity
c
      if (noc_aux_2.gt.0) then
        noc(1) = noc_aux_1 + 1
        noc(2) = noc_aux_2 - 1
        e_orbital = eig_lumo(1) - eig_homo(2)
        e_coul = 0.5d00*(int_LaLa + int_HbHb) - int_LaHb
        e_xc = int_vxc_H(2) - int_vxc_L(1)
        high_mult = e_orbital + e_coul + e_xc
        if (.not.rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     &     call errquit('fukui: rtdb_put of noc failed', 0, RTDB_ERR)
        do isp=1,ipol
          call ga_zero(g_dens(isp))
          call ga_dgemm('N', 'T', nbf_ao, nbf_ao,
     &                  noc(isp), 2d0/dble(ipol), g_movecs(isp),
     &                  g_movecs(isp), 0.0d00, g_dens(isp))
        enddo
        Exc(1) = 0.0d00
        Exc(2) = 0.0d00
        Ecoul = 0.0d00
        do isp = 1, ipol
          call ga_zero(g_vxc(isp))
        end do
        call xc_pot(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens,
     &               g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &               dbl_mb(irdens_atom), 
     &               int_mb(icetobfr), natoms)
        Exc_pert = Exc(1)
        diff_Exc = Exc_pert - Exc_zero
        high_mult = high_mult + diff_Exc
        if (me.eq.0.and.oprint_fukui)then
           write(LuOut,*) ' ------------------------------------'
           write(LuOut,'("  Alpha = ",i3," ; Beta = ",i3)')
     &                     noc(1), noc(2)
           write(LuOut,*) ' Contributions in atomic units:'
           write(LuOut,'("      Orbital Energy = ",f10.4)') e_orbital
           write(LuOut,'("      Coulomb Integrals = ",f10.4)') e_coul
           write(LuOut,'("      XC Integrals = ",f10.4)') e_xc
           write(LuOut,'("      XC Diff. Energy = ",f10.4)') diff_Exc
           write(LuOut,'("      High Multiplicity = ",f7.4," a.u.")')
     &                                          high_mult
           write(LuOut,'("                        = ",f7.2," eV")')
     &                                          high_mult*27.211
        end if
      end if
c
cc    Energy difference for low multiplicity
c
      if (noc_test.ge.2) then
        noc(1) = noc_aux_1 - 1
        noc(2) = noc_aux_2 + 1
        low_mult = eig_lumo(2) - eig_homo(1) +
     &              0.5d00*(int_LbLb + int_HaHa) - int_HaLb +
     &              int_vxc_H(1) - int_vxc_L(2)
        e_orbital = eig_lumo(2) - eig_homo(1)
        e_coul = 0.5d00*(int_LbLb + int_HaHa) - int_HaLb
        e_xc = int_vxc_H(1) - int_vxc_L(2)
        low_mult = e_orbital + e_coul + e_xc
        if (.not.rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     &     call errquit('fukui: rtdb_put of noc failed', 0, RTDB_ERR)
        do isp=1,ipol
          call ga_zero(g_dens(isp))
          call ga_dgemm('N', 'T', nbf_ao, nbf_ao,
     &                  noc(isp), 2d0/dble(ipol), g_movecs(isp),
     &                  g_movecs(isp), 0.0d00, g_dens(isp))
        enddo

        Exc(1) = 0.0d00
        Exc(2) = 0.0d00
        Ecoul = 0.0d00
        do isp = 1, ipol
          call ga_zero(g_vxc(isp))
        end do
        call xc_pot(rtdb, Exc, Ecoul,nExc, iVxc_opt, g_xcinv, g_dens,
     &               g_vxc, IOLGC, g_wght, g_xyz, g_nq,wght_GA, rho_n,
     &               dbl_mb(irdens_atom), 
     &               int_mb(icetobfr), natoms)
        Exc_pert = Exc(1)
        diff_Exc = Exc_pert - Exc_zero
        low_mult = low_mult + diff_Exc
        if (me.eq.0.and.oprint_fukui)then
           write(LuOut,*) ' ------------------------------------'
           write(LuOut,'("  Alpha = ",i3," ; Beta = ",i3)')
     &                     noc(1), noc(2)
           write(LuOut,*) ' Contributions:'
           write(LuOut,'("      Orbital Energy = ",f10.4)') e_orbital
           write(LuOut,'("      Coulomb Integrals = ",f10.4)') e_coul
           write(LuOut,'("      XC Integrals = ",f10.4)') e_xc
           write(LuOut,'("      XC Diff. Energy = ",f10.4)') diff_Exc
           write(LuOut,'("      Low Multiplicity = ",f7.4," a.u.")')
     &                                         low_mult
           write(LuOut,'("                       = ",f7.2," eV")')
     &                                         low_mult*27.211
           write(LuOut,*) ' ------------------------------------'
        end if
      end if
cc
c
      noc(1) = noc_aux_1
      noc(2) = noc_aux_2
      if (.not.rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('fukui: rtdb_put of noc failed', 0, RTDB_ERR)
      if (.not. ga_destroy(g_dens(1))) call errquit
     &    ('fukui: could not destroy g_dens(1)',0, GA_ERR)
      if (.not. ga_destroy(g_vxc(1))) call errquit
     &    ('fukui: could not destroy g_vxc(1)',0, GA_ERR)
      if (ipol.eq.2) then
        if (.not. ga_destroy(g_dens(2))) call errquit
     &    ('fukui: could not destroy g_dens(2)',0, GA_ERR)
        if (.not. ga_destroy(g_vxc(2))) call errquit
     &    ('fukui: could not destroy g_vxc(2)',0, GA_ERR)
      end if
c
      if (.not. ga_destroy(g_orb)) call errquit
     &   ('fukui: could not destroy g_orb', 0, GA_ERR)
      if (.not. ga_destroy(g_dens_HOMO(1))) call errquit
     &   ('fukui: could not destroy g_dens_HOMO', 0, GA_ERR)
      if (.not. ga_destroy(g_dens_LUMO(1))) call errquit
     &   ('fukui: could not destroy g_dens_LUMO', 0, GA_ERR)
      if (.not. ga_destroy(g_dens_HOMO(2))) call errquit
     &   ('fukui: could not destroy g_dens_HOMO_bet', 0, GA_ERR)
      if (.not. ga_destroy(g_dens_LUMO(2))) call errquit
     &   ('fukui: could not destroy g_dens_LUMO_bet', 0, GA_ERR)
      if (.not. ga_destroy(g_dens_ss)) call errquit
     &   ('fukui: could not destroy g_dens_ss', 0, GA_ERR)
      if (.not. ga_destroy(g_s)) call errquit
     &   ('fukui: could not destroy g_s', 0, GA_ERR)
      return
      end

      Subroutine mull_pop_fuk( geom, basis, iga_dens,iga_s)

C$Id$
      Implicit none
#include "errquit.fh"
      integer geom,basis
      integer iga_s              ! overlap   GA handle
      integer iga_dens           ! dens. mat GA handle
      integer iga_ps             ! product   GA handle

      integer natoms,nshells
      integer lPSmat,iPSmat,lqatom,iqatom,lqshell,iqshell
      integer iatom,ilo,ihi,nbf,max_at_bf2
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      logical status

#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "cscfps.fh"
#include "mafdecls.fh"

      if (oscfps) call pstat_on(ps_mull)
      

c*****************************************************************************

c     
      if(.not.geom_ncent(geom, natoms))
     &     call errquit(' exiting in mull_pop',0, GEOM_ERR)
      if( .not. bas_numcont(basis,nshells) )
     &     call errquit(' exiting in mull_pop',1, BASIS_ERR)
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in mull_op',1, BASIS_ERR)
      max_at_bf2 = 0
      do iatom = 1, natoms
        if (.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $       call errquit('mul_pop: bas_ce2bfr failed', iatom,
     &       BASIS_ERR)
        max_at_bf2 = max(max_at_bf2, ihi-ilo+1)
      enddo
      max_at_bf2 = max_at_bf2*max_at_bf2

      if(.not.MA_Push_Get(mt_dbl,max_at_bf2,'PS',lPSmat,iPSmat))
     &     call errquit(' exiting in mull_pop: insuff stack',21, MA_ERR)
      if(.not.MA_Push_Get(mt_dbl,natoms,'q atom',lqatom,iqatom))
     &     call errquit(' exiting in mull_pop: insuff stack',22, MA_ERR)
      if(.not.MA_Push_Get(mt_dbl,nshells,'q shell',lqshell,iqshell))
     &     call errquit(' exiting in mull_pop: insuff stack',3, MA_ERR)

      iga_PS=ga_create_atom_blocked(geom, basis, 'PS product')

      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,
     &     iga_dens,iga_s,0.d0,iga_PS)
      call mull_calc_fuk(basis,natoms, nshells,max_at_bf2,iga_PS,
     &     dbl_mb(iqatom),dbl_mb(iqshell),dbl_mb(iPSmat))
      
      call ga_sync

      if(.not.MA_Pop_Stack(lqshell))
     &     call errquit(' exiting in mull_pop',33, MA_ERR)
      if(.not.MA_Pop_Stack(lqatom))
     &     call errquit(' exiting in mull_pop',34, MA_ERR)
      if(.not.MA_Pop_Stack(lPSmat))
     &     call errquit(' exiting in mull_pop',35, MA_ERR)
      status= ga_destroy(iga_PS)
c
      if (oscfps) call pstat_off(ps_mull)
c
      return
      end
c
c
c
      Subroutine mull_calc_fuk(basis, natoms, nshells,max_at_bf2,iga_PS,
     &     qatom,qshell,PSmat)

      Implicit none
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "inp.fh"
      integer basis
      integer natoms,nshells
      integer iga_PS             ! product   GA handle
      integer ifirst,ilast,nbf_at,max_at_bf2
      integer ish1,ish2,ish,nn,iat,mu
      integer me,nproc, geom
      double precision psmu, coord(3), qnuc
      double precision qatom(natoms),qshell(nshells),PSmat(max_at_bf2)
      character*2 symbol
      character*16 tag
      character*32 element
      integer atn
c     
      me=ga_nodeid()
      nproc=ga_nnodes()
      
      call dfill(natoms,0.D0,qatom,1)
      call dfill(nshells,0.D0,qshell,1)

      if (.not. bas_geom(basis, geom)) call errquit
     $     ('mull_pop: bas_geom failed',basis, BASIS_ERR)

      do iat=me+1,natoms,nproc
        if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &       call errquit(' exiting in mull_pop',4, BASIS_ERR)
        call get_atom_block(iga_PS, basis,
     $       iat, iat, PSmat, nbf_at, nbf_at)
        mu=0
        do ish=ish1,ish2
          if (.not. bas_cn2bfr(basis,ish,ifirst,ilast))
     &         call errquit(' exiting in mull_pop.',5, BASIS_ERR)
          do nn=ifirst,ilast
            mu=mu+1
            psmu=PSmat((mu-1)*nbf_at+mu)
            qshell(ish)=qshell(ish)+psmu
          enddo
          qatom(iat)=qatom(iat)+qshell(ish)
        enddo
      enddo
      call ga_sync
      call ga_dgop(Msg_Mull1,qatom,natoms,'+')
      call ga_dgop(Msg_Mull2,qshell,nshells,'+')
      if(me.eq.0) then
        write(LuOut,1)
    1   format(/'    Atom       Condensed Fukui ')
        write(luout,11)
 11     format( ' -----------   ----------------')
        do iat=1,natoms
          if (.not.bas_ce2cnr(basis,iat,ish1,ish2))
     &         call errquit(' exiting in mull_pop',4, BASIS_ERR)
c
          if (.not. geom_cent_get(geom, iat, tag, coord, qnuc)) 
     $         call errquit('mull_pop: geom_cent_tag failed',0,
     &       GEOM_ERR)
c
          if (.not. geom_tag_to_element(tag, symbol, element, atn))
     $         symbol = 'X'
          if (ish2.ge.ish1) then
             write(LuOut,2) iat,symbol,nint(qnuc),qatom(iat) 
 2           format(1x,i4,1x,a2,i4,1x,f10.4)
          endif
        enddo
      endif
c
      call ga_sync

      return
      end

      Subroutine xc_pot(rtdb, Exc, ecoul,nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, g_wght, g_xyz,g_nq, 
     &                   wght_GA, rho_n, rdens_atom, 
     &                   cetobfr, natoms)
c
C$Id$
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
c      
      integer nExc
      integer iVxc_opt
      integer g_xcinv, g_dens(2), g_vxc(4), g_wght, g_xyz,g_nq
      integer natoms
      logical IOLGC, wght_GA
      integer rtdb
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "util.fh"
#include "sym.fh"
c
      integer cetobfr(2,natoms)
      double precision rho_n, rdens_atom(ipol*natoms*natoms)
      double precision jfac(4),kfac(4)
      integer g_jk(4), g_d(4)
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
c--> XC Energy
c
      double precision Exc(2)
      double precision ecoul ! [output]
c
c This driver routine solves for the XC energy and potential (Vxc) via 
c numerical quadrature methods. The results are obtained either by direct 
c numerical integration or by means of a LSQ fit of the Vxc to a set of 
c Gaussian functions. This fitted function can be used to evaluate Vxc
c via a summation of a series of 3-center overlap integrals (3OIs). The
c algorithms are formulated in terms of matrix products. See subsequent
c subroutines for further explanation.
c 
c              XC Energy and Potential Index Key, Vxc(pq,i)
c
c              Value of     |     Definition of index "i"
c            ipol     nExc  |    1        2        3       4
c           --------------------------------------------------
c              1        1   |   Vxc      
c              2        1   |   Vxc^up   Vxc^dw   
c              1        2   |   Vxc      
c              2        2   |   Vxc^up   Vxc^dw   
c
c           nTcols = ipol 
c
      integer me,nproc,i,nTrows,nTcols
      integer lTmat,iTmat,g_oep
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time, grid_on_file
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e,tot
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
      Exc(1)=0.d0
      Exc(2)=0.d0
      me=ga_nodeid()
      nproc=ga_nnodes()
c
      if (oprint_intermediate_xc)then
         write(*,*)' rtdb, Exc, nExc, iVxc_opt, g_xcinv: ',
     &               rtdb, Exc, nExc, iVxc_opt, g_xcinv
         write(*,*)' g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA: ',
     &               g_dens(1),g_vxc(1),IOLGC,g_wght,g_xyz,wght_GA
         write(*,*)' Fock XC matrix entering xc_getv: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
      endif
c
      if (abs(xfac(1)).gt.1e-8 .or. (.not. CDFIT))then
c       
c        Compute the exact exchange potential (as in Hartree-Fock calculations).
c       
         tol2e=10.d0**(-itol2e)
         call ga_sync
         if (oprint_time)call dft_tstamp(' Before call to fock_2e. ')
         if (ipol.eq.1) then
            kfac(1) = -0.5d0*xfac(1)
            jfac(1)=0.0d0
            if (.not. CDFIT) then
             jfac(2) = 1.0d0
             kfac(2) = 0d0
              g_vxc(2) = ga_create_atom_blocked(geom,ao_bas_han,'jk')
              g_dens(2)=g_dens(1)
              call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &             tol2e, oskel, g_dens(1), g_vxc(1), .false.)
              Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
              call ga_zero(g_vxc(2))
              ecoul = 0.5d0*ga_ddot(g_dens(1),g_vxc(2)) 
              call ga_dadd(1d0,g_vxc(1),1d0,g_vxc(2),g_vxc(1))
              if (.not. ga_destroy(g_vxc(2))) call errquit
     $             ('xc_getv: ga corrupt?',0, GA_ERR)
            else
              call fock_2e(geom, AO_bas_han, 1, jfac, kfac,
     &             tol2e, oskel, g_dens(1), g_vxc(1), .false.)
            Exc(1) = Exc(1)+0.5d0*ga_ddot(g_dens(1),g_vxc(1))
            endif
         else
            if (CDFIT) then
               jfac(1)=0.d0
               jfac(2)=0.d0
               kfac(1)=-1.0d0*xfac(1)
               kfac(2)=-1.0d0*xfac(1)
               call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &              tol2e, oskel, g_dens, g_vxc, .false.)
               Exc(1) = Exc(1)+0.5d0*(ga_ddot(g_dens(1),g_vxc(1)) +
     &              ga_ddot(g_dens(2),g_vxc(2)))
            else
               jfac(1) = 1.0d0
               jfac(2) = 0.0d0
               jfac(3) = 1.0d0
               jfac(4) = 0.0d0
               kfac(1) = 0.0d0
               kfac(2) = 1.0d0
               kfac(3) = 0.0d0
               kfac(4) = 1.0d0
               g_jk(1) = g_vxc(1) ! This assignment is assumed
               g_jk(2) = g_vxc(2)
               g_jk(3) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
               g_jk(4) = ga_create_atom_blocked(geom, ao_bas_han, 'jk')
               call ga_zero(g_jk(3))
               call ga_zero(g_jk(4))
               g_d(1)  = g_dens(1)
               g_d(2)  = g_dens(1)
               g_d(3)  = g_dens(2)
               g_d(4)  = g_dens(2)
               call fock_2e(geom, AO_bas_han, 4, jfac, kfac,
     &              tol2e, oskel, g_d(1), g_jk(1),.false. )
               call ga_zero(g_jk(1))
               call ga_zero(g_jk(3))
               ecoul = 0.5d0*( ! Alpha coulomb energy
     $              ga_ddot(g_dens(1),g_jk(1)) + 
     $              ga_ddot(g_dens(1),g_jk(3)))
               ecoul = ecoul + 0.5d0*( ! Beta coulomb energy
     $              ga_ddot(g_dens(2),g_jk(1)) + 
     $              ga_ddot(g_dens(2),g_jk(3)))
               exc(1) = exc(1) - xfac(1)*0.5d0*( ! All exchange energy
     $              ga_ddot(g_dens(1),g_jk(2)) +
     $              ga_ddot(g_dens(2),g_jk(4)))
               call ga_dadd(1.0d0, g_jk(1), 1.0d0, g_jk(3), g_jk(1))
               call ga_copy(g_jk(1), g_jk(3))
               call ga_dadd(1.0d0, g_jk(1), -xfac(1), g_jk(2),
     $              g_jk(1))
               call ga_dadd(1.0d0, g_jk(3), -xfac(1), g_jk(4),
     $              g_jk(2))
               if (.not. ga_destroy(g_jk(3))) call errquit
     $              ('xc_getv: ga corrupt?',0, GA_ERR)
               if (.not. ga_destroy(g_jk(4))) call errquit
     $              ('xc_getv: ga corrupt?',1, GA_ERR)
            endif
         endif
         if (oprint_time)call dft_tstamp('  After call to fock_2e. ')
         call ga_sync
c
c        Symmetrize Vxc?
c
c         if (oskel)then
c            call sym_symmetrize(geom, AO_bas_han, .false., g_vxc(1))
c            if (ipol.gt.1)then
c               call sym_symmetrize(geom, AO_bas_han, .false., 
c     &              g_vxc(2))
c            endif
c         endif
c       
c        Compute the exact exchange energy.
c       
      endif
c
      tot=-xfac(1)
      do i=1,numfunc
        tot=tot+xfac(i)+cfac(i)
      enddo
c
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &     grid_on_file))then
         grid_on_file = .false.
      endif
      if (abs(tot).gt.1e-8) then
         if(xcfit) then
            if (.not. bas_numbf(XC_bas_han,nbf_xc) )then
               call errquit('Exiting in getvxc.',1, BASIS_ERR)
            endif
            nTrows = nbf_xc
            nTcols = ipol
c         
c           Allocate scratch space for the "T" matrix.
c         
            if (.not.ma_push_get(MT_Dbl,nTrows*nTcols,'Tmat',lTmat,
     &         iTmat))call errquit('xc_getv: cannot allocate Tmat',0,
     &       MA_ERR)
            call dfill(nTrows*nTcols,0.D0,dbl_mb(iTmat),1)
         endif
            call grid_quadv0(rtdb, g_dens, g_vxc, nexc,rho_n,  Exc,
     .           dbl_mb(itmat))
      if(xcfit) then
c     
c           symmetrize the "T" vector
c     
         if (oskel)then
            call sym_vec_symmetrize(geom, xc_bas_han, Dbl_MB(iTmat))
            if (ipol.gt.1)then
               call sym_vec_symmetrize(geom, xc_bas_han, 
     &              Dbl_MB(iTmat+nbf_xc))
            endif
         endif
c      
            call xc_fitv(rtdb,Dbl_MB(iTmat), nTrows, nTcols,
     &                   g_vxc, g_xcinv, g_oep, IOLGC)
            if (oprint_time)call dft_tstamp(' After call to xc_fitv.  ')
            if (.not.ma_pop_stack(lTmat))
     &         call errquit('xc_getv: cannot pop stack',0, MA_ERR)
c
         endif
      endif
c
      if (oprint_intermediate_xc)then
         write(*,*)' Fock XC matrix leaving xc_pot: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
      endif
c
c
      return
      end

