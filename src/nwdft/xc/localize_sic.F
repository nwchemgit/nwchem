      Subroutine localize_sic(g_movecs, noc) 
c
C$Id: localize_sic.F,v 1.3 2000-04-28 17:54:38 jolguin Exp $
c
      implicit none
c      
      integer noc(2),isp,i,x,me
      integer g_movecs(2)
      integer g_uc(4),g_sc,g_over,l_sc,k_sc,l_c,k_c,nocc
      integer ncore,nloc,liloc,iiloc
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
     &   call errquit('localize_sic:push_get failed', 911)
c
      do isp=1,ipol
         if (me.eq.0) 
     $    call util_print_centered(6,
     $    'Foster-Boys orbital localization for the SIC approximation',
     &         40, .true.)
         nocc=noc(isp)
         if (me.eq.0.and.oprint_sic) write(LuOut,*) 
     $                         ' Number of Orbitals to localize',nocc
         if (nocc.gt.0) then
           do x= 1, 4
             if (.not. ga_create(MT_DBL, nbf, nbf, 'uc',
     $             nbf, 1, g_uc(x))) call errquit('g_uc_loc_sic',x)
           end do
c
           call int_dip_ga(AO_bas_han, AO_bas_han, g_uc(1), g_uc(2),
     &                     g_uc(3))
c
           if (.not. ga_create(MT_DBL, nbf, nbf, 'sc',
     $          nbf, 1, g_sc)) call errquit('g_sc_loc_sic',0)
           do x = 1, 3
             call ga_dgemm('n', 'n', nbf, nbf, nbf, 
     $                     1.0d0, g_uc(x), g_movecs, 0.0d0, g_sc)
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
     $          call errquit('ma for sc', 0)
           if (.not. ma_push_get(mt_dbl, 8*nbf, 'c', l_c, k_c))
     $          call errquit('ma for c', 0)
c
c     Localize the occupied orbitals
c
           if (.not. geom_num_core(geom,ncore)) ncore = 0
           if (ncore .gt. 0) then
              do i = 1, ncore
                 int_mb(iiloc+i-1) = i
              end do
              nloc = ncore
              call localizeFB(AO_bas_han, dbl_mb(k_c), dbl_mb(k_sc), 
     $                        nloc, int_mb(iiloc), nbf, nbf, 
     $                        g_movecs(isp), g_uc)
           end if
c
c     Localized the occupied
c
           do i = ncore+1, nocc
              int_mb(iiloc + i - ncore - 1) = i
           end do
           nloc = nocc - ncore
           call localizeFB(AO_bas_han, dbl_mb(k_c), dbl_mb(k_sc),
     $                     nloc, int_mb(iiloc), nbf, nbf, 
     $                     g_movecs(isp), g_uc)
c
           do x=1,4
             if (.not. ga_destroy(g_uc(x))) call errquit
     &        ('localize_sic: could not destroy g_uc', 0)
           end do
           if (.not. ga_destroy(g_sc)) call errquit
     &        ('localize_sic: could not destroy g_sc', 0)
           if (.not. ga_destroy(g_over)) call errquit
     &        ('localize_sic: could not destroy g_sc', 0)
           if (.not.ma_pop_stack(l_c))
     &        call errquit('localize_sic: cannot pop stack',0)
           if (.not.ma_pop_stack(l_sc))
     &        call errquit('localize_sic: cannot pop stack',0)
          end if
        end do
        if (.not.ma_pop_stack(liloc))
     &        call errquit('localize_sic: cannot pop stack',0)
        return
      end

