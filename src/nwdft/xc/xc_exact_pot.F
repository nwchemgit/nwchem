      Subroutine open_xc_exact_pot(g_rho_exact, g_tmp_exact,
     &                             g_addit_exact, g_movecs_aux,g_dens,
     &                             dif_lamda, lamda, p_lamda,
     &                             lamda_old, Ecoul_aux)
c
      implicit none
c
#include "global.fh"
#include "cdft.fh"
      integer g_rho_exact(2), g_tmp_exact(2), g_addit_exact(2),
     &        g_movecs_aux(2), g_dens(2)
      integer i
      double precision dif_lamda, lamda, p_lamda, lamda_old,
     &                 Ecoul_aux(2)
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      g_rho_exact(1) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                       'ga_rho_exact_1')
      g_tmp_exact(1) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                       'ga_tmp_exact_1')
      g_addit_exact(1) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                     'ga_addit_exact_1')
      g_movecs_aux(1) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                        'ga_movecs_aux_1')
      if (ipol.eq.2) then
        g_rho_exact(2) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                         'ga_rho_exact_2')
        g_tmp_exact(2) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                         'ga_tmp_exact_2')
        g_addit_exact(2) = ga_create_atom_blocked(geom,
     &                          AO_bas_han, 'ga_addit_exact_2')
        g_movecs_aux(2) = ga_create_atom_blocked(geom, AO_bas_han,
     &                                        'ga_movecs_aux_2')
      end if
      do i = 1, ipol
        call ga_zero(g_rho_exact(i))
        call ga_zero(g_tmp_exact(i))
        call ga_zero(g_addit_exact(i))
        call ga_zero(g_movecs_aux(i))
        call ga_copy(g_dens(i),g_rho_exact(i))
        call ga_scale(g_rho_exact(i),-1.0d00)
      end do
      do i = 1, numfunc
        xfac(i) = 0.0d00
        cfac(i) = 0.0d00
      end do
      dif_lamda = 0.d00
      lamda = 0.d00
      p_lamda = 0.d00
      lamda_old = 0.d00
      do i = 1 , ipol
        Ecoul_aux(i) = 100.0d0
      end do
      return
      end
c
cccc
      Subroutine xc_exact_pot(ecoul_aux, g_dens, g_vxc, 
     &                        g_rho_exact, dif_lamda, g_tmp_exact,
     &                        g_addit_exact)
c
      implicit none
#include "errquit.fh"
c      
c     integer noc(2)
      integer g_dens(2), g_vxc(4),  
     &        g_vxc_aux(4), g_rho_exact(2), g_rho_diff(2),
     &        g_addit_exact(2), g_tmp_exact(2), g_rho_zero(2)
      integer i_dif
c
#include "global.fh"
#include "cdft.fh"
c
      double precision jfac(4),kfac(4)
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      double precision fer_amal, dif_lamda
      double precision ecoul, ecoul_aux(2), ocup
c
      integer me,nproc
c     double precision zero,one,onem
c     parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e
c
      me=ga_nodeid()
      nproc=ga_nnodes()
c       
      g_vxc_aux(1) = ga_create_atom_blocked(geom,ao_bas_han,
     &                                      'g_vxc_aux_1')
      g_rho_diff(1) = ga_create_atom_blocked(geom,ao_bas_han,
     &                                      'g_rho_diff_1')
      tol2e=10.d0**(-itol2e)
      call ga_sync
      kfac(1) = 0.0d00
      jfac(1) = 0.0d0
      jfac(2) = 1.0d0
      kfac(2) = 0d0
      ocup = 2.0d00
      if (ipol.eq.2) ocup = 1.0d00
      do i_dif = 1 , ipol
        fer_amal = -1.0d00/(ocup*noc(i_dif))
        g_vxc_aux(2) = ga_create_atom_blocked(geom,ao_bas_han,'jk')
        g_rho_zero(1) = ga_create_atom_blocked(geom, ao_bas_han,
     &                                        'g_rho_zero')
        call ga_copy(g_dens(i_dif),g_rho_zero(1))
        g_rho_zero(2) = g_rho_zero(1)
        call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &               tol2e, oskel, g_rho_zero(1), g_vxc_aux(1),.false.)
        ecoul = 0.5d00*ga_ddot(g_rho_zero(1),g_vxc_aux(2))
        call ga_scale(g_vxc_aux(2),fer_amal)
        call ga_dadd(1.0d00, g_vxc(i_dif), 1.0d00, g_vxc_aux(2),
     &                 g_vxc(i_dif))
        if (dif_lamda.gt.0.0d00) then
           call ga_add(1.0d00, g_vxc(i_dif), 1.0d00,
     &                 g_tmp_exact(i_dif), g_vxc(i_dif))
           call ga_zero(g_vxc_aux(1))
           call ga_zero(g_vxc_aux(2))
           call ga_dadd(1.0d00, g_dens(i_dif), 1.0d00,
     &                  g_rho_exact(i_dif), g_rho_diff(1))
           g_rho_diff(2) = g_rho_diff(1)
           call fock_2e(geom, AO_bas_han, 2, jfac, kfac,
     &          tol2e, oskel, g_rho_diff(1), g_vxc_aux(1),.false.)
           ecoul_aux(i_dif) = 0.5d0*ga_ddot(g_rho_diff(1),
     &                                      g_vxc_aux(2))
           call ga_scale(g_vxc_aux(2),dif_lamda)
           call ga_copy(g_vxc_aux(2),g_addit_exact(i_dif))
           call ga_dadd(1.0d00, g_vxc(i_dif),1.0d00,
     &                  g_addit_exact(i_dif), g_vxc(i_dif))
        end if
        if (.not. ga_destroy(g_rho_zero(1))) call errquit
     $    ('xc_exact_pot: could not detroy g_rho_zero',0, GA_ERR)
        if (.not. ga_destroy(g_vxc_aux(2))) call errquit
     $        ('xc_exact_pot: could not detroy jk',0, GA_ERR)
      end do
      if (.not. ga_destroy(g_vxc_aux(1))) call errquit
     $             ('xc_exact_pot: could not detroy g_vxc_aux_1',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_rho_diff(1))) call errquit
     $             ('xc_exact_pot: could not detroy g_rho_diff_1',0,
     &       GA_ERR)
      return
      end
c
ccccc
      Subroutine parlam_xc_exact_pot(lamda, lamda_old, dif_lamda,
     &                               p_lamda, incre, g_tmp_exact,
     &                               g_addit_exact, g_movecs_aux,
     &                               g_movecs, iter, k_eval,
     &                               Ecoul_aux, split)
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
#include "util.fh"
#include "stdio.fh"
      integer g_addit_exact(2), g_movecs_aux(2), g_movecs(2),
     &        g_tmp_exact(2), k_eval(2)
c     integer noc(2)
      integer split, i, incre, iter, me, nproc, isp
      double precision reference, lamda, p_lamda, dif_lamda,
     &                 lamda_old, Ecoul_aux(2)
      me=ga_nodeid()
      nproc=ga_nnodes()
      if (iter.ge.iterations) then
        if (me.eq.0) write(LuOut,*) ' Lamda, NO SCF CONVERGENCE',
     &                                  lamda
      else
        if (me.eq.0)then
          write(LuOut,*) ' Lamda = ',lamda
          do isp = 1, ipol
            write(LuOut,*) ' Restriction = ',Ecoul_aux(isp)
            write(LuOut,*) (dbl_mb(k_eval(isp)+i-1),i=1,noc(1)+3)
          end do
        end if
      end if
      reference = 64.0d00
      if (iter.lt.iterations) then
        split = 0
        call ga_copy(g_movecs(1),g_movecs_aux(1))
        if (ipol.eq.2) call ga_copy(g_movecs(2),g_movecs_aux(2))
        if (lamda.lt.reference) then
          lamda_old = lamda
          p_lamda = p_lamda + 1.0d00
          lamda = 2.0d00**p_lamda
          incre = 0
        else
          lamda_old = lamda
          incre = incre + 1
          lamda = reference + dble(incre)*32.0d00
          if (lamda.eq.896.d00) lamda = 900.0d00
        end if 
      else
        split = split + 1
        do i = 1, ipol
          call ga_zero(g_addit_exact(i))
          call ga_copy(g_movecs_aux(i),g_movecs(i))
        end do
        if (lamda.le.reference) then
          if (split.eq.1) p_lamda = p_lamda - 1.0d00
        else
          if (split.eq.1) incre = incre - 1
        end if
        lamda = lamda_old + dif_lamda/dble(2*split)
      end if
      dif_lamda = lamda - lamda_old
      do i = 1, ipol
        call ga_dadd(1.0d0, g_tmp_exact(i), 1.0d0, g_addit_exact(i),
     &               g_tmp_exact(i))
      end do
      return
      end
c
cccccc
      Subroutine close_xc_exact_pot(g_rho_exact, g_tmp_exact,
     &                             g_addit_exact, g_movecs_aux)
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "cdft.fh"
      integer g_rho_exact(2), g_tmp_exact(2), g_addit_exact(2),
     &        g_movecs_aux(2)

      if (.not. ga_destroy(g_rho_exact(1))) call errquit
     &       ('xc_exact_pot: could not destroy g_rho_exact(1)', 0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp_exact(1))) call errquit
     &       ('xc_exact_pot: could not destroy g_tmp_exact_1', 0,
     &       GA_ERR)
      if (.not. ga_destroy(g_addit_exact(1))) call errquit
     &       ('xc_exact_pot: could not destroy g_addit_exact_1', 0,
     &       GA_ERR)
      if (.not. ga_destroy(g_movecs_aux(1))) call errquit
     &       ('xc_exact_pot: could not destroy g_movecs_aux_1', 0,
     &       GA_ERR)
      if (ipol.eq.2) then
        if (.not. ga_destroy(g_rho_exact(2))) call errquit
     &     ('xc_exact_pot: could not destroy g_rho_exact(2)', 0,
     &       GA_ERR)
        if (.not. ga_destroy(g_tmp_exact(2))) call errquit
     &     ('xc_exact_pot: could not destroy g_tmp_exact_2', 0,
     &       GA_ERR)
        if (.not. ga_destroy(g_addit_exact(2))) call errquit
     &     ('xc_exact_pot: could not destroy g_addit_exact_2', 0,
     &       GA_ERR)
        if (.not. ga_destroy(g_movecs_aux(2))) call errquit
     &       ('xc_exact_pot: could not destroy g_movecs_aux_2', 0,
     &       GA_ERR)
      end if
      return
      end
      



c $Id$
