      subroutine rjhtest(rtdb)
C$Id: rjhtest.F,v 1.12 1995-02-02 23:12:36 d3g681 Exp $
c$$$      implicit none
c$$$      integer rtdb
c$$$      call print_integs(rtdb, 'geometry','mo basis')
c$$$*      call scf(rtdb)
c$$$      end
c$$$      subroutine print_integs(rtdb, geom_name, basis_name)
c$$$      implicit none
c$$$#include "rtdb.fh"
c$$$#include "geom.fh"
c$$$#include "bas.fh"
c$$$#include "inp.fh"
c$$$#include "global.fh"
c$$$      integer rtdb
c$$$      character*(*) geom_name, basis_name
c$$$      integer geom
c$$$      integer bases(2)
c$$$      integer basis
c$$$      integer ri_basis, g_j, g_vinv, g_dens_dumm, g_j_real
c$$$      integer ncenter
c$$$      integer icent, jcent, kcent, lcent
c$$$      integer ish, jsh, ksh, lsh
c$$$      integer ishlo, ishhi, jshlo, jshhi, kshlo, kshhi, lshlo, lshhi
c$$$      integer i, j, k, l, ijkl
c$$$      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
c$$$      double precision eri(256), xi(3), xj(3), xk(3), xl(3)
c$$$      double precision qi, qj, qk, ql, sum
c$$$      character*16 tagi, tagj, tagk, tagl
c$$$      logical int_normalize
c$$$      integer ga_create_atom_blocked
c$$$      external int_normalize
c$$$      external ga_create_atom_blocked
c$$$      Double precision result
c$$$c
c$$$*      call ieeetrap
c$$$c
c$$$      if (.not. geom_create(geom, geom_name))
c$$$     $     call errquit('print_integs: geom_create failed?', 0)
c$$$      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
c$$$     $     call errquit('print_integs: geom_load failed', 0)
c$$$      If ( GA_NodeID() .eq. 0) then
c$$$         if (.not. geom_print(geom))
c$$$     $      call errquit('print_integs: geom_print failed', 0)
c$$$      EndIf
c$$$      if (.not. bas_create(basis, basis_name))
c$$$     $     call errquit('print_integs: basis create failed', 0)
c$$$      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
c$$$     $     call errquit('print_integs: basis load failed', 0)
c$$$      If ( GA_NodeID() .eq. 0) then
c$$$         write(6,*) ' basis before normalization '
c$$$         if (.not. bas_print(basis))
c$$$     $      call errquit('print_integs: basis print failed', 0)
c$$$      EndIf
c$$$      if (.not. int_normalize(basis))
c$$$     $     call errquit('print_integs: normalize failed', 0)
c$$$      If ( GA_NodeID() .eq. 0) then
c$$$         write(6,*) ' basis after normalization '
c$$$         if (.not. bas_print(basis))
c$$$     $      call errquit('print_integs: basis print failed', 0)
c$$$         if (.not. gbs_map_print(basis))
c$$$     $      call errquit('print_integs: gbs_map_print failed', 0)
c$$$      EndIf
c$$$c
c$$$c... open ri basis set
c$$$      if (.not.bas_create(ri_basis, 'ri-scf basis'))
c$$$     &    call errquit('rjhtest: bas_create failed',911)
c$$$      if (.not.bas_rtdb_load(rtdb,geom,ri_basis,'ri-scf basis'))
c$$$     &    call errquit('rjhtest: bas_rtdb_load failed',911)
c$$$      If ( GA_NodeID() .eq. 0) then
c$$$         if (.not. bas_print(ri_basis))
c$$$     $      call errquit('print_integs: ri_basis print failed', 0)
c$$$      EndIf
c$$$      if (.not. int_normalize(ri_basis))
c$$$     $     call errquit('print_integs: normalize failed', 0)
c$$$      If ( GA_NodeID() .eq. 0) then
c$$$         write(6,*) ' ri_basis after normalization '
c$$$         if (.not. bas_print(ri_basis))
c$$$     $      call errquit('print_integs: ri_basis print failed', 0)
c$$$         if (.not. gbs_map_print(ri_basis))
c$$$     $      call errquit('print_integs: gbs_map_print failed', 0)
c$$$      EndIf
c$$$      bases(1) = basis
c$$$      bases(2) = ri_basis
c$$$      call int_init(rtdb,2,bases)
c$$$      g_dens_dumm = ga_create_atom_blocked(geom,basis,' fake density')
c$$$      g_j = ga_create_atom_blocked(geom,basis,'ri coulomb matrix ')
c$$$      g_j_real =
c$$$     &       ga_create_atom_blocked(geom,basis,'real coulomb matrix ')
c$$$      call ga_zero(g_j)
c$$$      call schwarz_init(geom, basis)
c$$$      call guess_dens(geom,basis,g_dens_dumm)
c$$$c$$$      write(6,*)' density matrix'
c$$$c$$$      call ga_print(g_dens_dumm)
c$$$      call ri_invert_v(rtdb,ri_basis,geom,g_vinv)
c$$$c$$$      write(6,*)' inverse v '
c$$$c$$$      call ga_print(g_vinv)
c$$$c
c$$$      call ri_drive(rtdb,basis,ri_basis,geom,1.0d00,
c$$$     &       g_dens_dumm,g_vinv,g_j)
c$$$c$$$      write(6,*)' from ri coulomb matrix'
c$$$c$$$      call ga_print(g_j)
c$$$c
c$$$c get real coulomb energy 
c$$$c
c$$$      call ga_zero(g_j_real)
c$$$      call rhf_fock_2e(geom,basis,g_dens_dumm,g_j_real,
c$$$     &       1.0d-09,.true.,.false., .false.)
c$$$c$$$      write(6,*)' real coulomb matrix'
c$$$c$$$      call ga_print(g_j_real)
c$$$      Call GA_Sync
c$$$      call ga_dadd(1.0d00,g_j,-1.0d00,g_j_real, g_j_real)
c$$$      Result = ga_ddot(g_j_real,g_j_real)
c$$$      If ( GA_NodeID() .eq. 0 )
c$$$     $   write(6,*)'ddot of diff', result
c$$$c$$$      write(6,*)' difference!!!! '
c$$$c$$$      call ga_print(g_j_real)
c$$$      Call GA_Sync
c$$$      return
c$$$      write(6,2) geom_name(1:inp_strlen(geom_name)),
c$$$     $     basis_name(1:inp_strlen(basis_name))
c$$$ 2    format(/' Integrals for geometry "',a,'" with basis "',a,'"'/)
c$$$c
c$$$c     Partially redundant four-fold loop thru centers
c$$$c
c$$$      if (.not. geom_ncent(geom, ncenter))
c$$$     $     call errquit('print_integs: geom_ncent failed', 0)
c$$$      sum = 0.0d0
c$$$      do icent = 1, ncenter
c$$$         do jcent = 1, icent
c$$$            do kcent = 1, icent
c$$$               do lcent = 1, kcent
c$$$c
c$$$c     Figure out coords and contractions of each center
c$$$c
c$$$                  if (.not. geom_cent_get(geom, icent, tagi, xi, qi))
c$$$     $                 call errquit('print_integs: cent_get ?', icent)
c$$$                  if (.not. geom_cent_get(geom, jcent, tagj, xj, qj))
c$$$     $                 call errquit('print_integs: cent_get ?', jcent)
c$$$                  if (.not. geom_cent_get(geom, kcent, tagk, xk, qk))
c$$$     $                 call errquit('print_integs: cent_get ?', kcent)
c$$$                  if (.not. geom_cent_get(geom, lcent, tagl, xl, ql))
c$$$     $                 call errquit('print_integs: cent_get ?', lcent)
c$$$c                  
c$$$                  if (.not. bas_ce2cnr(basis, icent, ishlo, ishhi))
c$$$     $                 call errquit('print_integs: ce2cnr ?', icent)
c$$$                  if (.not. bas_ce2cnr(basis, jcent, jshlo, jshhi))
c$$$     $                 call errquit('print_integs: ce2cnr ?', jcent)
c$$$                  if (.not. bas_ce2cnr(basis, kcent, kshlo, kshhi))
c$$$     $                 call errquit('print_integs: ce2cnr ?', kcent)
c$$$                  if (.not. bas_ce2cnr(basis, lcent, lshlo, lshhi))
c$$$     $                 call errquit('print_integs: ce2cnr ?', lcent)
c$$$c
c$$$c     Partially redundant four-fold loop thru contractions
c$$$c
c$$$                  do ish = ishlo, ishhi
c$$$                     do jsh = jshlo, jshhi
c$$$                        do ksh = kshlo, kshhi
c$$$                           do lsh = lshlo, lshhi
c$$$c
c$$$c     Now compute the SP integrals over these shells
c$$$c
c$$$                              call genr70(
c$$$     $                             basis, ish, xi, jsh, xj,
c$$$     $                             basis, ksh, xk, lsh, xl,
c$$$     $                             eri)
c$$$c
c$$$c     Print the integrals out
c$$$c 
c$$$                              if (.not. bas_cn2bfr(basis, ish, ilo,ihi))
c$$$     $                             call errquit('print_integ: cn2bf?',0)
c$$$                              if (.not. bas_cn2bfr(basis, jsh, jlo,jhi))
c$$$     $                             call errquit('print_integ: cn2bf?',0)
c$$$                              if (.not. bas_cn2bfr(basis, ksh, klo,khi))
c$$$     $                             call errquit('print_integ: cn2bf?',0)
c$$$                              if (.not. bas_cn2bfr(basis, lsh, llo,lhi))
c$$$     $                             call errquit('print_integ: cn2bf?',0)
c$$$c
c$$$                              ijkl = 0
c$$$                              do i = ilo, ihi
c$$$                                 do j = jlo, jhi
c$$$                                    do k = klo, khi
c$$$                                       do l = llo, lhi
c$$$                                          ijkl = ijkl + 1
c$$$                                          if (abs(eri(ijkl)) .gt. 1e-6)
c$$$     $                                         write(6,1) i,j,k,l,
c$$$     $                                         eri(ijkl)
c$$$                                          sum = sum + abs(eri(ijkl))
c$$$                                       enddo
c$$$                                    enddo
c$$$                                 enddo
c$$$                              enddo
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$c
c$$$ 1    format(4i6,f20.9)
c$$$
c$$$      write(6,*) ' SUM of integrals ', sum
c$$$
c$$$      if (.not. geom_destroy(geom))
c$$$     $     call errquit('print_integs: geom_destroy failed', 0)
c$$$      if (.not.bas_destroy(basis))
c$$$     $     call errquit('print_integs: basis destroy failed', 0)
c$$$c
      end
