      subroutine rjhtest(rtdb)
      implicit none
      integer rtdb
      call print_integs(rtdb, 'geometry','mo basis')
*      call scf(rtdb)
      end
      subroutine print_integs(rtdb, geom_name, basis_name)
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "inp.fh"
#include "global.fh"
      integer rtdb
      character*(*) geom_name, basis_name
      integer geom
      integer bases(2)
      integer basis
      integer ri_basis, g_j, g_vinv, g_dens_dumm, g_j_real
      integer ncenter
      integer icent, jcent, kcent, lcent
      integer ish, jsh, ksh, lsh
      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
      integer i, j, k, l, ijkl
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      double precision eri(256), xi(3), xj(3), xk(3), xl(3)
      double precision qi, qj, qk, ql, sum
      character*16 tagi, tagj, tagk, tagl
      logical int_normalize
      integer ga_create_atom_blocked
      external int_normalize
      external ga_create_atom_blocked
c
*      call ieeetrap
c
      if (.not. geom_create(geom, geom_name))
     $     call errquit('print_integs: geom_create failed?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
     $     call errquit('print_integs: geom_load failed', 0)
      if (.not. geom_print(geom))
     $     call errquit('print_integs: geom_print failed', 0)
      if (.not. bas_create(basis, basis_name))
     $     call errquit('print_integs: basis create failed', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
     $     call errquit('print_integs: basis load failed', 0)
      write(6,*) ' basis before normalization '
      if (.not. bas_print(basis))
     $     call errquit('print_integs: basis print failed', 0)
      if (.not. int_normalize(basis))
     $     call errquit('print_integs: normalize failed', 0)
      write(6,*) ' basis after normalization '
      if (.not. bas_print(basis))
     $     call errquit('print_integs: basis print failed', 0)
      if (.not. gbs_map_print(basis))
     $     call errquit('print_integs: gbs_map_print failed', 0)
c
c... open ri basis set
      if (.not.bas_create(ri_basis, 'ri-scf basis'))
     &    call errquit('rjhtest: bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,ri_basis,'ri-scf basis'))
     &    call errquit('rjhtest: bas_rtdb_load failed',911)
      if (.not. bas_print(ri_basis))
     $     call errquit('print_integs: ri_basis print failed', 0)
      if (.not. int_normalize(ri_basis))
     $     call errquit('print_integs: normalize failed', 0)
      write(6,*) ' ri_basis after normalization '
      if (.not. bas_print(ri_basis))
     $     call errquit('print_integs: ri_basis print failed', 0)
      if (.not. gbs_map_print(ri_basis))
     $     call errquit('print_integs: gbs_map_print failed', 0)
      bases(1) = basis
      bases(2) = ri_basis
      call int_init(2,bases)
      g_dens_dumm = ga_create_atom_blocked(geom,basis,' fake density')
      g_j = ga_create_atom_blocked(geom,basis,'ri coulomb matrix ')
      g_j_real =
     &       ga_create_atom_blocked(geom,basis,'real coulomb matrix ')
      call ga_zero(g_j)
      call schwarz_init(geom, basis)
      call guess_dens(geom,basis,g_dens_dumm)
      write(6,*)' density matrix'
      call ga_print(g_dens_dumm)
      call ri_invert_v(rtdb,ri_basis,geom,g_vinv)
      write(6,*)' inverse v '
      call ga_print(g_vinv)
c
      call ri_drive(rtdb,basis,ri_basis,geom,1.0d00,
     &       g_dens_dumm,g_vinv,g_j)
      write(6,*)' from ri coulomb matrix'
      call ga_print(g_j)
c
c get real coulomb energy 
c
      call ga_zero(g_j_real)
      call rhf_fock_2e(geom,basis,g_dens_dumm,g_j_real,
     &       1.0d-09,.true.,.false.)
      write(6,*)' real coulomb matrix'
      call ga_print(g_j_real)
      call ga_dadd(1.0d00,g_j,-1.0d00,g_j_real, g_j_real)
      write(6,*)'ddot of diff',ga_ddot(g_j_real,g_j_real)
      write(6,*)' difference!!!! '
      call ga_print(g_j_real)
      return
      write(6,2) geom_name(1:inp_strlen(geom_name)),
     $     basis_name(1:inp_strlen(basis_name))
 2    format(/' Integrals for geometry "',a,'" with basis "',a,'"'/)
c
c     Partially redundant four-fold loop thru centers
c
      if (.not. geom_ncent(geom, ncenter))
     $     call errquit('print_integs: geom_ncent failed', 0)
      sum = 0.0d0
      do icent = 1, ncenter
         do jcent = 1, icent
            do kcent = 1, icent
               do lcent = 1, kcent
c
c     Figure out coords and contractions of each center
c
                  if (.not. geom_cent_get(geom, icent, tagi, xi, qi))
     $                 call errquit('print_integs: cent_get ?', icent)
                  if (.not. geom_cent_get(geom, jcent, tagj, xj, qj))
     $                 call errquit('print_integs: cent_get ?', jcent)
                  if (.not. geom_cent_get(geom, kcent, tagk, xk, qk))
     $                 call errquit('print_integs: cent_get ?', kcent)
                  if (.not. geom_cent_get(geom, lcent, tagl, xl, ql))
     $                 call errquit('print_integs: cent_get ?', lcent)
c                  
                  if (.not. bas_ce2cnr(basis, icent, ishlo, ishhi))
     $                 call errquit('print_integs: ce2cnr ?', icent)
                  if (.not. bas_ce2cnr(basis, jcent, jshlo, jshhi))
     $                 call errquit('print_integs: ce2cnr ?', jcent)
                  if (.not. bas_ce2cnr(basis, kcent, kshlo, kshhi))
     $                 call errquit('print_integs: ce2cnr ?', kcent)
                  if (.not. bas_ce2cnr(basis, lcent, lshlo, lshhi))
     $                 call errquit('print_integs: ce2cnr ?', lcent)
c
c     Partially redundant four-fold loop thru contractions
c
                  do ish = ishlo, ishhi
                     do jsh = jshlo, jshhi
                        do ksh = kshlo, kshhi
                           do lsh = lshlo, lshhi
c
c     Now compute the SP integrals over these shells
c
                              call genr70(
     $                             basis, ish, xi, jsh, xj,
     $                             basis, ksh, xk, lsh, xl,
     $                             eri)
c
c     Print the integrals out
c 
                              if (.not. bas_cn2bfr(basis, ish, ilo,ihi))
     $                             call errquit('print_integ: cn2bf?',0)
                              if (.not. bas_cn2bfr(basis, jsh, jlo,jhi))
     $                             call errquit('print_integ: cn2bf?',0)
                              if (.not. bas_cn2bfr(basis, ksh, klo,khi))
     $                             call errquit('print_integ: cn2bf?',0)
                              if (.not. bas_cn2bfr(basis, lsh, llo,lhi))
     $                             call errquit('print_integ: cn2bf?',0)
c
                              ijkl = 0
                              do i = ilo, ihi
                                 do j = jlo, jhi
                                    do k = klo, khi
                                       do l = llo, lhi
                                          ijkl = ijkl + 1
                                          if (abs(eri(ijkl)) .gt. 1e-6)
     $                                         write(6,1) i,j,k,l,
     $                                         eri(ijkl)
                                          sum = sum + abs(eri(ijkl))
                                       enddo
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
 1    format(4i6,f20.9)

      write(6,*) ' SUM of integrals ', sum

      if (.not. geom_destroy(geom))
     $     call errquit('print_integs: geom_destroy failed', 0)
      if (.not.bas_destroy(basis))
     $     call errquit('print_integs: basis destroy failed', 0)
c
      end
