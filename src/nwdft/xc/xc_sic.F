      Subroutine xc_sic(rtdb, nExc, iVxc_opt, g_xcinv,
     &                  g_dens, g_vxc,  g_wght, g_xyz,g_nq, 
     &                  wght_GA, rdens_atom, 
     &                  cetobfr, natoms,
     &                  g_movecs, noc,  totsic, i_degen)
c
C$Id: xc_sic.F,v 1.12 2000-04-27 20:29:11 zzhang Exp $
c
      implicit none
c      
      integer iVxc_opt, nExc, g_xcinv, g_dens(*), g_vxc(4), g_wght,
     &        g_xyz,g_nq
      integer ik, ij, noc(2), isp, n_states, naux_ipol, size_mat, 
c    &        g_sic_dens(2),g_orb,g_vxc_orb(4),g_movecs(2),
     &        g_orb,g_vxc_orb(4),g_movecs(2),
     &        g_matm, g_vectu,  g_coul_orb, g_vxc_sig(2),
     &        l_vec_aux, i_vec_aux, j, 
     &        l_vec_u, i_vec_u, l_consti, i_consti, l_degen,
     &        aux_ij,  i_degen(2), i_gauxsig, l_gauxsig
      double precision Exc_orb(2), tot_sic_xc, tot_sic_coul, ecoul_orb,
     &                 u_array, res_int,const_ci, totsic
      integer natoms
      logical  wght_GA
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
#include "stdio.fh"
c
      integer cetobfr(2,natoms)
      double precision rho_n, rdens_atom(ipol*natoms*natoms)
      double precision jfac(4),kfac(4)
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer me,nproc
c     integer me,nproc,i,nTrows,nTcols
c     integer lTmat,iTmat
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time, oprint_sic
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e
c******************************************************************************
c
c Compute the matrix elements for the XC potential of the OEP and the
c energy of the SIC approximation
c
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
      oprint_sic = util_print('SIC information', print_high)
c
      me=ga_nodeid()
      nproc=ga_nnodes()
c
      if (iVxc_opt.eq.0) then
c         
c           Option A:  Compute via direct numerical quadrature.
c         
        if (ipol.eq.2) then
          if (noc(2).eq.0) call ga_zero(g_vxc(2))
        end if
c
        if (me.eq.0.and.oprint_sic) write(LuOut,*)
     &               ' Starting SIC by orbital..'
        g_vxc_orb(1) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                      'Vxc_orb_alpha')
        g_vxc_orb(2) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                        'Vxc_orb_beta')
        g_orb = ga_create_atom_blocked(geom, ao_bas_han, 'Orbsic')
        g_sic_dens(1) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                          'Densic_alpha')
        g_sic_dens(2) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                           'Densic_beta')
        g_coul_orb= ga_create_atom_blocked(geom, ao_bas_han,
     &                                       'V_coul_orb')
        g_vxc_sig(1) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                          'Vxc_sig_alpha')
        if (ipol.eq.2) then
          g_vxc_sig(2) = ga_create_atom_blocked(geom, ao_bas_han,
     &                            'Vxc_sig_beta')
        end if
c
        tot_sic_xc = 0.d0
        tot_sic_coul = 0.d0
c
        if (ipol.eq.2) then
          naux_ipol = 2
          sic_orb_occ = 1
        else
          naux_ipol = 1
          sic_orb_occ = 2
        end if
c
        do 300 isp = 1,ipol
c
          n_states = noc(isp)
          size_mat = 0
          if (noc(isp).ne.0) then
            do j = 1,noc(isp)
              l_degen = Int_MB(i_degen(isp) + j - 1)
              if (l_degen.eq.0) then
                size_mat = size_mat + 1
              else
                if (me.eq.0.and.oprint_sic) write(LuOut,*)
     &            ' HOMO = orbital ',j
              end if
            end do
          end if
          call ga_zero(g_vxc_sig(isp))
          if (size_mat.gt.0) then
            if (.not. ga_create(mt_dbl, size_mat,size_mat,'Mat_m',
     &          size_mat, size_mat, g_matm))
     &        call errquit('dft_main0d: error creating g_matm',0)
            if (.not. ga_create(mt_dbl, size_mat, 1, 'Vet_u',
     &           size_mat, 1, g_vectu))
     &        call errquit('dft_main0d: error creating Vet_u',0)
            if (.not.MA_Push_Get(MT_Dbl, size_mat, 'const_i',
     &               l_consti, i_consti))
     &       call errquit('xc_sic: cannot allocate vec aux',0)
          end if
c
          if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'vec aux',
     &                 l_vec_aux, i_vec_aux))
     &         call errquit('xc_sic: cannot allocate vec aux',0)
          if (.not.MA_Push_Get(MT_Dbl, n_states, 'vec u',
     &                 l_vec_u, i_vec_u))
     &         call errquit('xc_sic: cannot allocate vec u',0)
          if (.not.ma_push_get(mt_int, n_states, 'g_aux_sig',
     &                         l_gauxsig, i_gauxsig))
     &         call errquit('xc_sic:push_get failed', 0)
c
c   Start loop of total states by spin        
c
cc
          do 100 ik = 1,n_states
            if (me.eq.0.and.oprint_sic) write(LuOut,*) ' Orbital ',ik
            aux_ik = Int_MB(i_degen(isp) + ik - 1)
            call ga_zero(g_sic_dens(1))
            call ga_zero(g_sic_dens(2))
            call ga_zero(g_orb)
            call ga_zero(g_vxc_orb(1))
            call ga_zero(g_vxc_orb(2))
            call ga_zero(g_coul_orb)
            call ga_get(g_movecs(isp), 1, nbf_ao, ik, ik,
     &                    Dbl_MB(i_vec_aux), nbf_ao)
            call ga_put(g_orb, 1, nbf_ao, ik, ik, 
     &                    Dbl_MB(i_vec_aux), nbf_ao)
            call ga_dgemm('n', 't', nbf_ao, nbf_ao, nbf_ao, 1.0d00,
     &                      g_orb, g_orb, 0.d00, g_sic_dens(1))
c
c    g_sic_dens(1) is the orbital density of the state ik
c
            if (abs(xfac(1)).gt.1e-8 )then
              if (me.eq.0.and.oprint_sic) write(LuOut,*) 
     &                                    ' NO SIC in Coulomb term'
            else
              if (me.eq.0.and.oprint_sic) write(LuOut,*)
     &                            ' SIC approximation in Coulomb term'
              kfac(1) = 0.d00
              jfac(1) = 1.0d0
              jfac(2) = 1.0d0
              kfac(2) = 0d0
              call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &             tol2e, oskel, g_sic_dens(1), g_vxc_orb(1),.false.)
              ecoul_orb = 0.0d0
              ecoul_orb = 0.5d0*ga_ddot(g_sic_dens(1), g_vxc_orb(1))
c
c   ecoul_orb is the orbital coulomb energy 
c
              tot_sic_coul = tot_sic_coul - dble(sic_orb_occ)*ecoul_orb
              call ga_copy(g_vxc_orb(1), g_coul_orb)
            end if
c
            ipol = 2
            call ga_zero(g_vxc_orb(1))
            Exc_orb(1) = 0.0d0
            Exc_orb(2) = 0.0d0
            if (oprint_time)call dft_tstamp
     &                           ('Before call to xc_quadv0.')
            sic_orb_spin = isp
            sic_orb_index = 0
            call xc_quadv0(rtdb, nExc, g_vxc_orb, g_sic_dens,
     &                     g_wght, g_xyz, g_nq, wght_GA, Exc_orb,
     &                     rho_n, rdens_atom,
     &                     cetobfr, natoms)
c
c   Exc_orb(1) is the exchange-correlation energy by orbital
c
            tot_sic_xc = tot_sic_xc - dble(sic_orb_occ)*Exc_orb(1)
c
            if (aux_ik.ne.1) then
              call ga_zero(g_orb)
              call ga_add(1.0d0, g_coul_orb, 1.0d0, g_vxc_orb(1),
     &                    g_orb)      
c
c  g_orb contains the coulomb + exchange-correlation potential evaluated with
c        the ik-orbital
c
              u_array = -ga_ddot(g_sic_dens(1), g_orb)
              Dbl_MB(i_vec_u + ik - 1) = u_array
            end if

            call ga_zero(g_vxc_orb(1))
            call ga_zero(g_vxc_orb(2))
            sic_orb_spin = isp
            sic_orb_index = 1
            call xc_quadv0(rtdb, nExc, g_vxc_orb, g_dens,
     &                     g_wght, g_xyz,g_nq, wght_GA, Exc_orb,
     &                     rho_n, rdens_atom, 
     &                     cetobfr, natoms)
            sic_orb_index = 0
            do ij=1,noc(isp)
              aux_ij=Int_MB(i_degen(isp) + ij - 1)
              if (aux_ij.ne.1.and.aux_ik.ne.1) then
                call ga_zero(g_orb)
                call ga_zero(g_sic_dens(2))
                call ga_get(g_movecs(isp), 1, nbf_ao, ij, ij, 
     &                      Dbl_MB(i_vec_aux) ,nbf_ao)
                call ga_put(g_orb, 1, nbf_ao, ij, ij, 
     &                      Dbl_MB(i_vec_aux) ,nbf_ao)
                call ga_dgemm('n', 't', nbf_ao, nbf_ao, nbf_ao,
     &                        1.0d00, g_orb, g_orb, 0.d00, 
     &                        g_sic_dens(2))
                res_int = -ga_ddot(g_sic_dens(2), g_vxc_orb(2))
                if (ij.eq.ik) res_int = 1.0d00 + res_int
                call ga_fill_patch(g_matm, ik, ik, ij, ij, res_int)
              end if
            end do
c
            if (size_mat.gt.0.and.aux_ik.ne.1) then
              if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao,'i_gauxsig',
     &                            nbf_ao, nbf_ao,
     &                            int_mb(i_gauxsig + ik - 1)))
     &             call errquit(
     '              ' xc_sic: error creating i_gauxsig',0)
              call ga_copy(g_vxc_orb(2), int_mb(i_gauxsig + ik - 1))
            end if

            call ga_add(-1.0d0, g_vxc_orb(1), 1.0d00, 
     &                  g_vxc_sig(isp), g_vxc_sig(isp))
            if (naux_ipol.eq.1) then
               ipol=1
            end if
            if (oprint_time)call dft_tstamp
     &                             (' After call to xc_quadv0.')
  100     continue
c
          do ik=1,noc(isp)
            aux_ik = Int_MB(i_degen(isp) + ik - 1)
            if (aux_ik.ne.1) then
              call ga_zero(g_orb)
              call ga_zero(g_sic_dens(2))
              call ga_get(g_movecs(isp), 1, nbf_ao, ik, ik, 
     &                    Dbl_MB(i_vec_aux) ,nbf_ao)
              call ga_put(g_orb, 1, nbf_ao, ik, ik, 
     &                    Dbl_MB(i_vec_aux), nbf_ao)
              call ga_dgemm('n', 't', nbf_ao, nbf_ao, nbf_ao, 
     &                      1.0d00, g_orb, g_orb, 0.d00,
     &                      g_sic_dens(2))
              res_int = ga_ddot(g_sic_dens(2), g_vxc_sig(isp))
              Dbl_MB(i_consti + ik - 1) = res_int-Dbl_MB(i_vec_u +
     &                                                   ik - 1)
            end if
          end do
c
          if (size_mat.gt.0) then
            call ga_zero(g_vectu)
            call ga_put(g_vectu, 1, size_mat, 1, 1, Dbl_MB(i_consti),
     &                  size_mat)
            call ga_lu_solve_seq('n', g_matm, g_vectu)
            call ga_get(g_vectu, 1, size_mat, 1, 1, Dbl_MB(i_consti),
     &                  size_mat)
            call ga_zero(g_orb)
            do ik=1,noc(isp)
              aux_ik = Int_MB(i_degen(isp) + ik - 1)
              if (aux_ik.ne.1) then
                const_ci= Dbl_MB(i_consti + ik - 1)
                call ga_add(1.0d0, g_orb, const_ci,
     &                      int_mb(i_gauxsig + ik - 1), g_orb)
              end if
            end do
            call ga_add(1.0d0, g_orb, 1.0d00, g_vxc_sig(isp),
     &                  g_vxc_sig(isp))
            do ik = 1, size_mat
              if (.not. ga_destroy(int_mb(i_gauxsig + ik - 1))) 
     &             call errquit ('xc_sic: could not destroy i_gauxsig',
     &                            0)
            end do
            if (.not.ma_pop_stack(l_gauxsig))
     &            call errquit('xc_sic: cannot pop stack',0)
            if (.not.ma_pop_stack(l_vec_u))
     &            call errquit('xc_sic: cannot pop stack',0)
            if (.not.ma_pop_stack(l_vec_aux))
     &            call errquit('xc_sic: cannot pop stack',0)
c
c  g_vxc_sig is the exchange-correlation potential associate to the OEP
c
          end if
          call ga_add(1.0d0, g_vxc_sig(isp), 1.0d00, g_vxc(isp),
     &                g_vxc(isp))
c
          if (size_mat.gt.0) then
            if (.not.ma_pop_stack(l_consti))
     &          call errquit('xc_sic: cannot pop stack',0)
            if (.not. ga_destroy(g_matm)) call errquit
     &         ('xc_getv: could not destroy g_matm', 0)
            if (.not. ga_destroy(g_vectu)) call errquit
     &         ('xc_getv: could not destroy g_vectu', 0)
          end if     
  300   continue
        totsic = tot_sic_xc + tot_sic_coul
        sic_orb_occ = 0
        if (me.eq.0.and.oprint_sic) write(LuOut,*) 
     &    ' tot_sic_coul, tot_sic_xc, tot_sic: ',tot_sic_coul,
     &                               tot_sic_xc, totsic
c
        if (.not. ga_destroy(g_vxc_orb(1))) call errquit
     &       ('xc_sic: could not destroy g_vxc_orb(1)', 0)
        if (.not. ga_destroy(g_vxc_orb(2))) call errquit
     &       ('xc_sic: could not destroy g_vxc_orb(2)', 0)
        if (.not. ga_destroy(g_orb)) call errquit
     &       ('xc_sic: could not destroy g_orb', 0)
        if (.not. ga_destroy(g_sic_dens(1))) call errquit
     &       ('xc_sic: could not destroy g_sic_dens(1)', 0)
        if (.not. ga_destroy(g_sic_dens(2))) call errquit
     &       ('xc_sic: could not destroy g_sic_dens(2)', 0)
        if (.not. ga_destroy(g_coul_orb)) call errquit
     &       ('xc_sic: could not destroy g_coul_orb', 0)
        if (.not. ga_destroy(g_vxc_sig(1))) call errquit
     &       ('xc_sic: could not destroy g_vxc_sig', 0)
        if (ipol.eq.2) then
          if (.not. ga_destroy(g_vxc_sig(2))) call errquit
     &     ('xc_sic: could not destroy g_vxc_sig', 0)
        end if
c
       elseif (iVxc_opt.eq.1 )then
           call errquit
     &      ('xc_sic_pert: SIC + XC fitting not implemented', 0)
         endif
c
      if (oprint_intermediate_xc)then
         write(*,*)' Fock XC matrix leaving xc_sic: '
         call ga_print(g_vxc(1))
         if(ipol.eq.2)call ga_print(g_vxc(2))
      endif
c
      return
      end

