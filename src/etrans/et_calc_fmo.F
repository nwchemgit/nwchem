      logical function et_calc_fmo(rtdb,Vba)
      implicit none
c
c     This function calculates the electron transfer matrix elements.
c     So far only frontier couplings are calculated (total 4).
c     Only closed shells are possible at the moment.
*
#include "errquit.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "cscf.fh"
#include "stdio.fh"
#include "numerical_constants.fh"
c
      integer rtdb, me
      integer i, j, k, ilo, ihi, jlo, jhi
      integer g_s_f1, g_mos_f1, ma_occ_f1(2), ma_evls_f1(2), nhomo_f1
      integer g_s_f2, g_mos_f2, ma_occ_f2(2), ma_evls_f2(2), nhomo_f2
      integer g_s_d,  g_mos_d,  ma_occ_d(2),  ma_evls_d(2),  nhomo_d
      integer g_s_p, g_mos_p
      integer g_evls_d, g_mos_d_o, g_mos_p_o
      integer et_nmo(2),nbf_f1,nbf_f2,nbf_d,nset
      integer g_c, g_fock, g_fock_tmp
c
      double precision Vba, E1, E2, toev, ths
c
      character*255 basis_name
c      character*2 pname(2)

      logical movecs_read_header
      external movecs_read_header

      parameter (toev=27.2114d+00)
      parameter (ths=1000.0d+00)

      et_calc_fmo=.true.

      me=ga_nodeid()

      if (me.eq.0) then
      Call Util_Print_Centered(LuOut,'FMO ET Calculation'
     &     ,40,.true.)
      write(LuOut,*)
      endif
c
c for printing:
c
      call util_print_push
      call util_print_rtdb_load(rtdb,"et")
c
c     FRAGMENT 1
c
      if(.not.movecs_read_header('frag1.movecs',title,basis_name,
     &   scftype,nbf_f1,nset,et_nmo,1))
     &   call errquit('et_calc: movecs_read_header failed for frag1',
     &       me, DISK_ERR)

      if (nset.ne.1)
     &   call errquit('et_calc: only HF/DFT orbitals allowed frag1',
     &   nset,INPUT_ERR)
c
      call et_fmo_getdata(rtdb,'frag1',g_s_f1,g_mos_f1,
     &   ma_evls_f1, ma_occ_f1)
      Call GA_Sync
c
      nhomo_f1=0
      do i = 1,nbf_f1
        nhomo_f1 = nhomo_f1 + dbl_mb(ma_occ_f1(2)+i-1)
      enddo
      nhomo_f1 = nhomo_f1/2
c
c order of pop matters
      if (.not.ma_pop_stack(ma_occ_f1(1))) 
     &   call errquit('could not destroy MA',ma_occ_f1(1),MA_ERR)
c
      if (.not.ma_pop_stack(ma_evls_f1(1))) 
     &   call errquit('could not destroy MA',ma_evls_f1(1),MA_ERR)
c
c     FRAGMENT 2
c
      if(.not.movecs_read_header('frag2.movecs',title,basis_name,
     &   scftype,nbf_f2,nset,et_nmo,1))
     &   call errquit('et_calc: movecs_read_header failed for frag2',
     &       me, DISK_ERR)
c
      if (nset.ne.1)
     &   call errquit('et_calc: only HF/DFT orbitals allowed frag2',
     &   nset,INPUT_ERR)
c
      call et_fmo_getdata(rtdb,'frag2',g_s_f2,g_mos_f2,
     &   ma_evls_f2, ma_occ_f2)
      Call GA_Sync
c
      nhomo_f2=0
      do i = 1,nbf_f2
         nhomo_f2 = nhomo_f2 + dbl_mb(ma_occ_f2(2)+i-1)
      enddo
      nhomo_f2 = nhomo_f2/2
c
      if (.not.ma_pop_stack(ma_occ_f2(1))) 
     &   call errquit('could not destroy MA',ma_occ_f2(1),MA_ERR)
c
      if (.not.ma_pop_stack(ma_evls_f2(1))) 
     &   call errquit('could not destroy MA',ma_evls_f2(1),MA_ERR)
c
c
c     DIMER
c
      if(.not.movecs_read_header('dimer.movecs',title,basis_name,
     &   scftype,nbf_d,nset,et_nmo,1))
     &   call errquit('et_calc: movecs_read_header failed for dimer',
     &       me, DISK_ERR)
c
      if (nset.ne.1)
     &   call errquit('et_calc: only HF/DFT orbitals allowed dimer',
     &   nset,INPUT_ERR)
c
      call et_fmo_getdata(rtdb,'dimer',g_s_d,g_mos_d,ma_evls_d,ma_occ_d)
      Call GA_Sync
c
      nhomo_d=0
      do i = 1,nbf_d
         nhomo_d = nhomo_d + dbl_mb(ma_occ_d(2)+i-1)
      enddo
      nhomo_d = nhomo_d/2
c
      if(.not.ga_duplicate(g_mos_d, g_evls_d, 'd_evls'))
     &   call errquit('et_calc: ga_duplicate failed evls_d', 0, GA_ERR)
      call ga_zero(g_evls_d)
c
c copied from ga_trace_diag
      call ga_distribution(g_evls_d, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         ilo = max(ilo,jlo)
         ihi = min(ihi,jhi)
         do i = ilo,ihi
            call ga_put(g_evls_d, i, i, i, i,
     &         dbl_mb(ma_evls_d(2)+i-1)*toev, 1)
         enddo
      endif
c
      if (.not.ma_pop_stack(ma_occ_d(1))) 
     &   call errquit('could not destroy MA',ma_occ_d(1),MA_ERR)
c
      if (.not.ma_pop_stack(ma_evls_d(1))) 
     &   call errquit('could not destroy MA',ma_evls_d(1),MA_ERR)
c
c     PAIR
c
      if(.not.ga_duplicate(g_mos_d, g_mos_p, 'p_mos'))
     &   call errquit('et_calc: ga_duplicate failed p_mos', 0, GA_ERR)
      call ga_zero(g_mos_p)
c
      if(.not.ga_duplicate(g_mos_d, g_s_p, 'p_s'))
     &   call errquit('et_calc: ga_duplicate failed p_s', 0, GA_ERR)
      call ga_zero(g_s_p)
c
      call ga_copy_patch('n', g_mos_f1, 1, nbf_f1, 1, nbf_f1,
     &   g_mos_p, 1, nbf_f1, 1, nbf_f1)
c
      call ga_copy_patch('n', g_mos_f2, 1, nbf_f2, 1, nbf_f2,
     &   g_mos_p, nbf_f1+1, nbf_d, nbf_f1+1, nbf_d)
c
      call ga_copy_patch('n', g_s_f1, 1, nbf_f1, 1, nbf_f1,
     &   g_s_p, 1, nbf_f1, 1, nbf_f1)
c
      call ga_copy_patch('n', g_s_f2, 1, nbf_f2, 1, nbf_f2,
     &   g_s_p, nbf_f1+1, nbf_d, nbf_f1+1, nbf_d)
c
      call ga_matpow(g_s_d, FP_Half, 1d-12)
      call ga_matpow(g_s_p, FP_Half, 1d-12)
c
      if(.not.ga_duplicate(g_mos_d, g_mos_d_o, 'mos_d_o'))
     &   call errquit('et_calc: ga_duplicate failed mos_d_o', 0, GA_ERR)
c
      call ga_dgemm('n','n',nbf_d,nbf_d,nbf_d,1.0d0,
     &   g_s_d, g_mos_d, 0.0d0, g_mos_d_o)
c
      if(.not.ga_duplicate(g_mos_d, g_mos_p_o, 'mos_p_o'))
     &   call errquit('et_calc: ga_duplicate failed mos_p_o', 0, GA_ERR)
c
      call ga_dgemm('n','n',nbf_d,nbf_d,nbf_d,1.0d0,
     &   g_s_p,g_mos_p, 0.0d0, g_mos_p_o)
c
      if(.not.ga_duplicate(g_mos_d, g_c, 'c'))
     &   call errquit('et_calc: ga_duplicate failed c', 0, GA_ERR)
      call ga_zero(g_c)
c
      call ga_dgemm('t','n',nbf_d,nbf_d,nbf_d,1.0d0,
     &   g_mos_p_o,g_mos_d_o, 0.0d0, g_c)
c
      if(.not.ga_duplicate(g_mos_d, g_fock, 'fock'))
     &   call errquit('et_calc: ga_duplicate failed fock', 0, GA_ERR)
      call ga_zero(g_fock)

      if(.not.ga_duplicate(g_mos_d, g_fock_tmp, 'fock_tmp'))
     &  call errquit('et_calc: ga_duplicate failed fock_tmp', 0, GA_ERR)
      call ga_zero(g_fock_tmp)
c
      call ga_dgemm('n','n',nbf_d,nbf_d,nbf_d,1.0d0,
     &   g_c,g_evls_d, 0.0d0, g_fock_tmp)
c
      call ga_dgemm('n','t',nbf_d,nbf_d,nbf_d,1.0d0,
     &   g_fock_tmp,g_c, 0.0d0, g_fock)

      if (.not.ga_destroy(g_fock_tmp)) 
     &  call errquit('ga_destroy g_fock_tmp failed',g_fock_tmp,
     &       GA_ERR)
c
      if (me.eq.0) then
        write(LuOut,9995) 'Frag1 HOMO: ', nhomo_f1,
     &   ' Frag2 HOMO: ', nhomo_f2, ' Dimer HOMO: ', nhomo_d
        write(LuOut,*)
        write(*,"(5A10)") 'F1 #Mo','F2 #Mo','V, meV','E1, meV','E2, meV'
c
        do i=nhomo_f1, nhomo_f1+1
          do j=nhomo_f2+nbf_f1, nhomo_f2+nbf_f1+1
            Call GA_Get(g_fock, j, j, i, i, Vba, 1)
            Call GA_Get(g_fock, i, i, i, i, E1, 1)
            Call GA_Get(g_fock, j, j, j, j, E2, 1)
            write(LuOut,9996) i, j-nbf_f1, Vba*ths, E1*ths, E2*ths
          enddo
        enddo
      endif
c
 9995 format(A13,I3,A13,I3,A13,I3)
 9996 format(2I10,3F10.3)
c
      call util_print_pop
c
c     destroying
c
      if (.not.ga_destroy(g_s_f1))
     &  call errquit('could not destroy s_f1',g_s_f1, GA_ERR)
c
      if (.not.ga_destroy(g_mos_f1))
     &  call errquit('could not destroy mos_f1',g_mos_f1, GA_ERR)

      if (.not.ga_destroy(g_s_f2))
     &  call errquit('could not destroy s_f2',g_s_f2, GA_ERR)
c
      if (.not.ga_destroy(g_mos_f2))
     &  call errquit('could not destroy mos_f2',g_mos_f2, GA_ERR)

      if (.not.ga_destroy(g_s_d))
     &  call errquit('could not destroy s_d',g_s_d, GA_ERR)
c
      if (.not.ga_destroy(g_mos_d))
     &  call errquit('could not destroy mos_d',g_mos_d, GA_ERR)
c
      if (.not.ga_destroy(g_s_p))
     &  call errquit('could not destroy s_p',g_s_p, GA_ERR)
c
      if (.not.ga_destroy(g_mos_p))
     &  call errquit('could not destroy mos_p',g_mos_p, GA_ERR)
c
      if (.not.ga_destroy(g_mos_d_o))
     &  call errquit('could not destroy mos_d_o',g_mos_d_o, GA_ERR)
c
      if (.not.ga_destroy(g_mos_p_o))
     &  call errquit('could not destroy mos_p_o',g_mos_p_o, GA_ERR)
c
      if (.not.ga_destroy(g_c))
     &  call errquit('could not destroy c',g_c, GA_ERR)
c
      if (.not.ga_destroy(g_fock))
     &  call errquit('could not destroy fock',g_fock, GA_ERR)
c
      if (.not.ga_destroy(g_evls_d))
     &  call errquit('could not destroy evls_d',g_evls_d, GA_ERR)
c
      et_calc_fmo=.true.
      return
      end
c $Id$
