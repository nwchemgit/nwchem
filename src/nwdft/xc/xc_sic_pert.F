      Subroutine xc_sic_pert(rtdb, nExc, iVxc_opt, g_xcinv,
     &                       g_dens, g_vxc, IOLGC, g_wght, g_xyz, g_nq, 
     &                       wght_GA, rho_n, rdens_atom, cntoce,
     &                       cntobfr, cetobfr, natoms,
     &                       g_movecs, noc, k_eval, totsic)
c
C$Id: xc_sic_pert.F,v 1.1 1999-08-10 21:21:14 jolguin Exp $
c
      implicit none
c      
      integer nExc, iVxc_opt, natoms, rtdb, g_xcinv, g_dens(2),
     &        g_vxc(4), g_wght, g_xyz, g_nq, ik, noc(2), isp,
     &        isp, n_states, naux_ipol, i, g_dens_orb(2), g_orb,
     &        g_vxc_orb(4), g_movecs(2), k_eval(2), g_coul_orb,
     &        i_temp_vec, l_temp_vec, fac_oc
      double precision Exc_orb(2), tot_sic_xc, tot_sic_coul, ecoul_orb,
     &                 totsic
      logical IOLGC, wght_GA, oprint_sic
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"
c
      integer cntoce(nshells_ao), cntobfr(2,nshells_ao), 
     &        cetobfr(2,natoms)
      double precision rho_n, rdens_atom(ipol*natoms*natoms)
      double precision jfac(4),kfac(4)
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer me,nproc,i
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e,tot
c
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
      oprint_sic = util_print('SIC information', print_high)
      me=ga_nodeid()
      nproc=ga_nnodes()
c
      tot=-xfac(1)
      do i=1,10
        tot=tot+xfac(i)+cfac(i)
      enddo
c
      if (abs(tot).gt.1e-8) then
         if (iVxc_opt.eq.0) then
           if (me.eq.0.and.oprint_sic) write(LuOut,*) 
     &                 ' Starting SIC perturbative approximation...'
           g_vxc_orb(1) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                         'Vxc_orb')
           g_vxc_orb(2) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                         'Vxc_orb')
           g_orb = ga_create_atom_blocked(geom, ao_bas_han, 'Orbsic')
           g_dens_orb(1) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                         'Densic')
           g_dens_orb(2) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                         'Densic')
           g_coul_orb = ga_create_atom_blocked(geom, ao_bas_han,
     &                                         'V_coul_orb')
           tot_sic_xc=0.d0
           tot_sic_coul=0.d0
           if (ipol.eq.2) then
             naux_ipol=2
             fac_oc=1
           else
             naux_ipol=1
             fac_oc=2
           end if
c
           if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'temp vec',
     &        l_temp_vec, i_temp_vec))
     &        call errquit('xc_sic_pert: cannot allocate temp vec',0)
c
           do isp = 1,ipol
             n_states = noc(isp)
ccccc   Start loop of total states by spin        
             do ik = 1,n_states
               if (me.eq.0.and.oprint_sic) write(LuOut,*) ' State ',ik
               call ga_zero(g_dens_orb(1))
               call ga_zero(g_dens_orb(2))
               call ga_zero(g_orb)
               call ga_zero(g_vxc_orb(1))
               call ga_zero(g_vxc_orb(2))
               call ga_zero(g_coul_orb)
               call ga_get(g_movecs(isp), 1, nbf_ao, ik, ik,
     &                     DBl_MB(i_temp_vec), nbf_ao)
               call ga_put(g_orb, 1, nbf_ao, ik, ik, 
     &                     DBl_MB(i_temp_vec), nbf_ao)
               call ga_dgemm('n', 't', nbf_ao, nbf_ao, nbf_ao, 1.0d00,
     &                       g_orb, g_orb, 0.d00, g_dens_orb(1))
               if (abs(xfac(1)).gt.1e-8 )then
                 if (me.eq.0.and.oprint_sic) write(LuOut,*) 
     &                                       ' no sic in coulomb part'
               else
cccc g_dens_orb(1) is the orbital density of orb(ik) states
                 kfac(1) = 0.d00
                 jfac(1) = 1.0d0
                 jfac(2) = 1.0d0
                 kfac(2) = 0d0
                 call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &                        tol2e, oskel, g_dens_orb(1), g_vxc_orb(1))
                 ecoul_orb = 0.0d0
                 ecoul_orb = 0.5d0*ga_ddot(g_dens_orb(1), g_vxc_orb(1))
                 tot_sic_coul = tot_sic_coul - dble(fac_oc)*ecoul_orb
cccc ecoul_orb is the orbital coulomb energy 
                 call ga_copy(g_vxc_orb(1), g_coul_orb)
               end if
               ipol = 2
               call ga_zero(g_vxc_orb(1))
               Exc_orb(1) = 0.0d0
               Exc_orb(2) = 0.0d0
               if (oprint_time)call dft_tstamp(
     &                                   ' Before call to xc_quadv0.')
               sic_orb_spin = isp
               sic_orb_occ = 0
               sic_orb_index = 0
               call xc_quadv0(rtdb, nExc, g_vxc_orb, g_dens_orb,
     &                            g_wght, g_xyz, g_nq, wght_GA,
     &                            Exc_orb, rho_n, rdens_atom, cntoce,
     &                            cntobfr, cetobfr, natoms)
               tot_sic_xc = tot_sic_xc - dble(fac_oc)*Exc_orb(1)
cccc Exc_orb(1) is the orbital exchange-correlation energy
               if (naux_ipol.eq.1) ipol = 1
               if (oprint_time)call dft_tstamp(
     &                                  ' After call to xc_quadv0.')
             end do
           end do
c
           if (.not.ma_pop_stack(l_temp_vec))
     &        call errquit('xc_sic_pert: cannot pop stack',0)
c
           totsic = tot_sic_xc + tot_sic_coul
c
           if (me.eq.0.and.oprint_sic) write(LuOut,*) 
     &                                            ' tot_sic=',totsic
c
           if (.not. ga_destroy(g_vxc_orb(1))) call errquit
     &      ('xc_sic_pert: could not destroy g_vxc_orb(1)', 0)
           if (.not. ga_destroy(g_vxc_orb(2))) call errquit
     &      ('xc_sic_pert: could not destroy g_vxc_orb(2)', 0)
           if (.not. ga_destroy(g_orb)) call errquit
     &      ('xc_sic_pert: could not destroy g_orb', 0)
           if (.not. ga_destroy(g_dens_orb(1))) call errquit
     &      ('xc_sic_pert: could not destroy g_dens_orb(1)', 0)
           if (.not. ga_destroy(g_dens_orb(2))) call errquit
     &      ('xc_sic_pert: could not destroy g_dens_orb(2)', 0)
           if (.not. ga_destroy(g_coul_orb)) call errquit
     &      ('xc_sic_pert: could not destroy g_coul_orb', 0)
c
         elseif (iVxc_opt.eq.1 )then
c
           call errquit
     &      ('xc_sic_pert: SIC + XC fitting not implemented', 0)
c
         endif
c
      endif
c
      return
      end

