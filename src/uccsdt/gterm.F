      subroutine uccsdt_makeg(g_omega,d_g,spini,spinj)
c
c     Pure spin:  g(a,i) = Sum(u,v,j) t(u,v,i,j)C(u,a,spini)C(v,j,spinj)
c     Mixed spin: g(a,j) = Sum(u,v,i) t(u,v,i,j)C(u,a,spini)C(v,i,spini)
c     where u & v are SO's (no spin), e,m,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
      integer d_g          ! [output] DRA handle for g(a,i)
      integer g_omega      ! [input] GA handle for Omega(u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer g_g          ! Temporary GA for g(u,i) and g(a,i)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer k_c,l_c      ! Temporary MA for MO coeffs
      integer lenui,lenu,lenv,lenuv,lena,leni,lenai
      integer symui,symi,symu,symij,symuv,symj,symuj,syma,symv
      integer ioff(nw_max_nbf_small,0:7)
      integer ij,ijlo,ijhi
      integer vlo,alo
      integer i,j
      integer dummy
c
c     Addressing & allocate temporary GA
c
      lenui = 0
      do symui = 0,7
       do symi = 0,7
        symu = ieor(symui,symi)
        do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
         ioff(i,symui) = lenui
         lenui = lenui + bf_per_ir(symu)
        end do
       end do
      end do
      if (.not.ga_create(mt_dbl,lenui,1,'g', 
     $ -1,1,g_g)) call errquit
     $ ('uccsdt_makeg: room for g?',lenui)
      call ga_zero(g_g)
c
c     Transformed MO coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf*no(spini), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spini))
      call ga_get(g_hole(spini),1,nbf,1,no(spini),
     $ dbl_mb(k_c),nbf)
c
c     Data parallel transform 
c     Pure spin:  Omega(u,v,i,j)C(v,i) = tmp(u,j)
c     Mixed spin: Omega(u,v,i,j)C(v,j) = tmp(u,i)
c
      call ga_distribution(g_omega,ga_nodeid(),ijlo,ijhi,dummy,dummy)
      ij = 0
      do symij = 0,7
       symuv = symij
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
          ij = ij + 1
          if (ij.ge.ijlo .and. ij.le.ijhi) then
           do symv = 0,7
            symu = ieor(symuv,symv)
            symui = ieor(symu,symi)
            symuj = ieor(symu,symj)
            lenu = bf_per_ir(symu)
            lenv = bf_per_ir(symv)
            vlo = bf_per_ir_cum(symv) + 1
            lenuv = lenu * lenv
            if (lenuv.gt.0) then
             if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $        l_buf, k_buf)) 
     $       call errquit('ma? lenuv',lenuv)
             if (.not. ma_push_get(mt_dbl, lenu, 'tmp',
     $        l_tmp, k_tmp)) 
     $       call errquit('ma? lenu',lenu)
             call ga_get(g_omega,1,lenuv,ij,ij,dbl_mb(k_buf),lenu)
c
c     (uv) -> (u) 
c
             if (spini .eq. spinj) then
              write(*,*) 'correct???'
              call dgemm('n','n',lenu,1,lenv,
     $         1.0d0,dbl_mb(k_buf),lenu,
     $         dbl_mb(k_c),nbf,
     $         0.0d0,dbl_mb(k_tmp),lenu)
              call ga_acc(g_g,
     $         ioff(j,symuj)+1,ioff(j,symuj)+lenu,
     $         1,1,dbl_mb(k_tmp),lenu,1.0d0)
             else
              call dgemm('n','n',lenu,1,lenv,
     $         1.0d0,dbl_mb(k_buf),lenu,
     $         dbl_mb(k_c),nbf,
     $         0.0d0,dbl_mb(k_tmp),lenu)
              call ga_acc(g_g,
     $         ioff(i,symui)+1,ioff(i,symui)+lenu,
     $         1,1,dbl_mb(k_tmp),lenu,1.0d0)
             end if
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf)
            end if
           end do
          end if
         end do
        end do
       end do
      end do
c
c     Sequential transform (ui) -> (ai)
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
      if (.not. ma_push_get(mt_dbl, nbf*nv(spini), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*nv',nbf*nv(spini))
      call ga_get(g_part(spini),1,nbf,no(spini)+1,nmo,
     $ dbl_mb(k_c),nbf)
      do symui = 0,7
       do symi = 0,7
        symu = ieor(symui,symi)
        syma = symu
        lenu = bf_per_ir(symu)
        lena = nv_sym(syma,spini)
        leni = no_sym(symi,spini)
        lenui = lenu * leni
        lenai = lena * leni
        vlo = bf_per_ir_cum(symu) + 1
        alo = v_sym(1,symu,spini)
        if (.not. ma_push_get(mt_dbl, lenai, 'tmp',
     $   l_buf, k_buf)) 
     $   call errquit('ma? lenui',lenai)
        if (.not. ma_push_get(mt_dbl, lenui, 'tmp',
     $   l_tmp, k_tmp)) 
     $   call errquit('ma? lenui',lenui)
        call ga_get(g_g,
     $   ioff(i,symui)+1,ioff(i,symui)+lenui,
     $   1,1,dbl_mb(k_tmp),lenu)
        call dgemm('t','n',lena,leni,lenu,
     $   1.0d0,dbl_mb(k_c),nbf,
     $   dbl_mb(k_tmp),lenu,
     $   0.0d0,dbl_mb(k_buf),lena)
        call ga_put(g_g,
     $   ioff(i,symui)+1,ioff(i,symui)+lenai,
     $   1,1,dbl_mb(k_buf),lena)
       end do
      end do
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
c
c     Write one symmetry block of G to disk
c
csh   if (.not. uccsdt_ampfile_add_t1(d_g, spini, g_g))
csh  $ call errquit('write_t1 failed',0)
      if (.not. ga_destroy(g_g)) call errquit('GA 999',0)
c
      return
      end
