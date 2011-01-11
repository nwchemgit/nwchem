      Subroutine localize_sic(g_movecs, evals)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
c     integer noc(2)
      integer isp, x, me, i_degen, l_degen, n_levels,
     &        evals(2), ij, ik, n_orbitals, i_level, i_orb, j,
     &        aux_levels
      double precision start
     
      integer g_movecs(2)
      integer g_uc(4),g_sc,g_over,l_sc,k_sc,l_c,k_c,nocc
c     integer ncore,nloc,liloc,iiloc,i
      integer nloc,liloc,iiloc
      logical geom_num_core
      external geom_num_core
c
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical oprint_sic
c
c******************************************************************************
c
c
      me=ga_nodeid()
      oprint_sic = util_print('SIC information', print_high)
c
      if (.not.ma_push_get
     &   (MT_Int,nbf,'map of orb to loc',liloc,iiloc))
     &   call errquit('localize_sic:push_get failed', 911, MA_ERR)
c
      do isp=1,ipol
c
         nocc=noc(isp)
         if (.not.MA_Push_Get(MT_Int, nocc, 'levels',
     &        l_degen, i_degen))
     &        call errquit('localize_sic: cannot allocate levels',0,
     &       MA_ERR)
         if (nocc.ne.0) then
           aux_levels = 1
           n_orbitals = 1
           Int_MB(i_degen + aux_levels - 1) = n_orbitals
           do j = noc(isp),2,-1
             start = (dbl_mb(evals(isp) + j - 1) -
     &                dbl_mb(evals(isp) + j - 2))
             if (start.le.1.0d-04) then
               n_orbitals = n_orbitals + 1
             else
               n_orbitals = 1
               aux_levels = aux_levels + 1
             end if
             Int_MB(i_degen + aux_levels - 1) = n_orbitals
           end do
           n_levels = aux_levels
         else
           n_levels = 0
         end if
c
         if (nocc.gt.0) then
           do x= 1, 4
             if (.not. ga_create(MT_DBL, nbf, nbf, 'uc',
     $             nbf, 1, g_uc(x))) call errquit('g_uc_loc_sic',x,
     &       GA_ERR)
           end do
c
           call int_dip_ga(AO_bas_han, AO_bas_han, g_uc(1), g_uc(2),
     &                     g_uc(3))
c
           if (.not. ga_create(MT_DBL, nbf, nbf, 'sc',
     $          nbf, 1, g_sc)) call errquit('g_sc_loc_sic',0, GA_ERR)
           do x = 1, 3
             call ga_dgemm('n', 'n', nbf, nbf, nbf, 
     $                     1.0d0, g_uc(x), g_movecs(isp), 0.0d0, g_sc)
             call ga_copy_patch('n',g_sc,1,nbf,1,nbf,g_uc(x),1,
     &                           nbf,1,nbf)
           end do
           g_over  = ga_create_atom_blocked(geom, AO_bas_han,
     &                                    'g_over_loc_sic')
           call ga_zero(g_over)
           call int_1e_ga(AO_bas_han, AO_bas_han, g_over, 'overlap',
     &                    .false.)
           call ga_dgemm('n', 'n', nbf, nbf, nbf, 
     $                   1.0d0, g_over, g_movecs(isp), 0.0d0, g_uc(4))
           if (.not. ma_push_get(mt_dbl, 8*nbf, 'sc', l_sc, k_sc))
     $          call errquit('ma for sc', 0, MA_ERR)
           if (.not. ma_push_get(mt_dbl, 8*nbf, 'c', l_c, k_c))
     $          call errquit('ma for c', 0, MA_ERR)
c
c     Localize the occupied orbitals
c
c          if (.not. geom_num_core(geom,ncore)) ncore = 0
c          if (ncore .gt. 0) then
c             do i = 1, ncore
c                int_mb(iiloc+i-1) = i
c             end do
c             nloc = ncore
c             call localizeFB(AO_bas_han, dbl_mb(k_c), dbl_mb(k_sc), 
c    $                        nloc, int_mb(iiloc), nbf, nbf, 
c    $                        g_movecs(isp), g_uc)
c          end if
c
c     Localized the occupied
c
c          do i = ncore+1, nocc
c             int_mb(iiloc + i - ncore - 1) = i
c          end do
c          nloc = nocc - ncore
c
c          ik = noc(isp) -  Int_MB(i_degen) + 1   ! If HOMO is not involved
           ik = noc(isp) + 1
           ij = 0
c          do 100 i_level = 2,n_levels            ! If HOMO is not involved
           do 100 i_level = 1,n_levels
             n_orbitals = Int_MB(i_degen + i_level - 1)
             do 200 i_orb = 1,n_orbitals
               ik = ik - 1
               ij = ij + 1
               int_mb(iiloc + ij - 1) = ik
  200        continue
  100      continue
           nloc = ij
           if (nloc.gt.0) then
             if (me.eq.0) 
     $        call util_print_centered(6,
     $    'Foster-Boys orbital localization for the SIC approximation',
     &           40, .true.)
              if (me.eq.0.and.oprint_sic) write(LuOut,*) 
     $                           ' Number of Orbitals to localize',nloc
              call localizeFB(AO_bas_han, dbl_mb(k_c), dbl_mb(k_sc),
     $                        nloc, int_mb(iiloc), nbf, nbf, 
     $                        g_movecs(isp), g_uc)
           end if
c
           do x=1,4
             if (.not. ga_destroy(g_uc(x))) call errquit
     &        ('localize_sic: could not destroy g_uc', 0, GA_ERR)
           end do
           if (.not. ga_destroy(g_sc)) call errquit
     &        ('localize_sic: could not destroy g_sc', 0, GA_ERR)
           if (.not. ga_destroy(g_over)) call errquit
     &        ('localize_sic: could not destroy g_sc', 0, GA_ERR)
           if (.not.ma_pop_stack(l_c))
     &        call errquit('localize_sic: cannot pop stack',0, MA_ERR)
           if (.not.ma_pop_stack(l_sc))
     &        call errquit('localize_sic: cannot pop stack',0, MA_ERR)
          end if
          if (.not.ma_pop_stack(l_degen))
     &              call errquit('localize_sic: cannot pop stack',0,
     &       MA_ERR)
        end do
        if (.not.ma_pop_stack(liloc))
     &        call errquit('localize_sic: cannot pop stack',0,
     &       MA_ERR)
        return
      end

