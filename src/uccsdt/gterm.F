      subroutine uccsdt_makeg(g_omega,g_g_big,spini,spinj)
csh   subroutine uccsdt_makeg(g_omega,d_g,spini,spinj)
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
csh   integer d_g          ! [output] DRA handle for g(a,i)
csh
      integer g_g_big      ! [output] DRA handle for g(a,i)
csh
      integer g_omega      ! [input] GA handle for Omega(u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer spina
      integer g_g          ! Temporary GA for g(u,i) and g(a,i)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer k_c,l_c      ! Temporary MA for MO coeffs
      integer lenui,lenu,lenv,lenuv,lena,leni,lenai
      integer symi,symu,symij,symuv,symj,syma,symv
      integer symioff(0:7),ioff(nbf,0:7)
      integer ij,ijlo,ijhi
      integer i,j,jlo,jhi
      integer dummy
csh
      logical odebug
      odebug = .false.
csh
      spina = spini
c
c     Addressing & allocate temporary GA
c
      lenai = 0
      do symi = 0,7
       syma = symi
       symioff(symi) = lenai
       if (odebug) write(*,*) 'symi,symioff = ',
     $  symi,symioff(symi)
       do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
        lenai = lenai + nv_sym(syma,spina)
       end do
      end do
      if (odebug) write(*,*) 'lenai = ',lenai
      if (.not.ga_create(mt_dbl,lenai,1,'g_big',-1,1,g_g_big))
     $ call errquit ('uccsdt_makeg: room for g_big?',lenai)
c
c     Addressing & allocate temporary GA
c
      lenui = 0
      do symi = 0,7
       symu = symi
       do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
        ioff(i,symi) = lenui
        if (odebug) write(*,*) 'i,symi,ioff = ',i,symi,ioff(i,symi)
        lenui = lenui + bf_per_ir(symu)
       end do
      end do
      if (odebug) write(*,*) 'lenui = ',lenui
      if (.not.ga_create(mt_dbl,lenui,1,'g',-1,1,g_g))
     $ call errquit ('uccsdt_makeg: room for g?',lenui)
      call ga_zero(g_g)
c
c     Transformed MO coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf*no(spini), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spini),0)
      call ga_get(g_part(spini),1,nbf,1,no(spini),
     $ dbl_mb(k_c),nbf)
      if (odebug) call ma_print(dbl_mb(k_c),nbf,no(spini),'MO C(P)')
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
        symu = symi
        symv = symj
        jlo = o_sym(1,symj,spinj)
        jhi = o_sym(2,symj,spinj)
        do j = jlo,jhi
         do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
          ij = ij + 1
          if (ij.ge.ijlo .and. ij.le.ijhi) then
           if (odebug) write(*,*) 'i,j,symj,symij = ',i,j,symj,symij
           lenu = bf_per_ir(symu)
           lenv = bf_per_ir(symv)
           lenuv = lenu * lenv
           if (lenuv.gt.0) then
            if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $       l_buf, k_buf)) 
     $      call errquit('ma? lenuv',lenuv,0)
            if (.not. ma_push_get(mt_dbl, lenu, 'tmp',
     $       l_tmp, k_tmp)) 
     $      call errquit('ma? lenu',lenu,0)
            call ga_get(g_omega,1,lenuv,ij,ij,dbl_mb(k_buf),lenu)
c
c     (uv) -> (u) 
c
            if (odebug) call ma_print(dbl_mb(k_buf),lenu,lenv,'UV')
            if (odebug) call ma_print(dbl_mb(k_c+(j-jlo)*nbf),nbf,1,'C')
            call dgemm('n','n',lenu,1,lenv,
     $       1.0d0,dbl_mb(k_buf),lenu,
     $       dbl_mb(k_c+(j-jlo)*nbf),nbf,
     $       0.0d0,dbl_mb(k_tmp),lenu)
            if (odebug) call ma_print(dbl_mb(k_tmp),lenu,1,'U')
            if (odebug) write(*,*) 'accumulating ... in ',
     $       ioff(i,symi)+1,'~',ioff(i,symi)+lenu,',1~1'
            call ga_acc(g_g,
     $       ioff(i,symi)+1,ioff(i,symi)+lenu,
     $       1,1,dbl_mb(k_tmp),lenu,1.0d0)
            if (.not. ma_pop_stack(l_tmp))
     $       call errquit('ma_pop?',l_tmp,0)
            if (.not. ma_pop_stack(l_buf))
     $       call errquit('ma_pop?',l_buf,0)
           end if
          end if
         end do
        end do
       end do
      end do
c
c     Sequential transform (ui) -> (ai)
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c,0)
      if (.not. ma_push_get(mt_dbl, nbf*nv(spina), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*nv',nbf*nv(spina),0)
      call ga_get(g_part(spina),1,nbf,no(spina)+1,nmo,
     $ dbl_mb(k_c),nbf)
      do symi = 0,7
       if (odebug) write(*,*) 'symi = ',symi
       symu = symi
       syma = symu
       lenu = bf_per_ir(symu)
       lena = nv_sym(syma,spini)
       leni = no_sym(symi,spini)
       lenui = lenu * leni
       lenai = lena * leni
       if ((lenai.gt.0).and.(lenui.gt.0)) then
        if (.not. ma_push_get(mt_dbl, lenai, 'buf',
     $   l_buf, k_buf)) 
     $   call errquit('ma? lenai',lenai,0)
        if (.not. ma_push_get(mt_dbl, lenui, 'tmp',
     $   l_tmp, k_tmp)) 
     $   call errquit('ma? lenui',lenui,0)
        call ga_get(g_g,
     $   ioff(i,symi)+1,ioff(i,symi)+lenui,
     $   1,1,dbl_mb(k_tmp),lenu)
        if (odebug) call ma_print(dbl_mb(k_tmp),lenu,leni,'UI')
        call dgemm('t','n',lena,leni,lenu,
     $   1.0d0,dbl_mb(k_c),nbf,
     $   dbl_mb(k_tmp),lenu,
     $   0.0d0,dbl_mb(k_buf),lena)
        if (odebug) call ma_print(dbl_mb(k_buf),lena,leni,'AI')
        if (odebug) write(*,*) 'putting ... in ',
     $   symioff(symi)+1,'~',symioff(symi)+lenai,',1~1'
        call ga_put(g_g_big,
     $   symioff(symi)+1,symioff(symi)+lenai,
     $   1,1,dbl_mb(k_buf),lena)
        if (.not. ma_pop_stack(l_tmp))
     $   call errquit('ma_pop?',l_tmp,0)
        if (.not. ma_pop_stack(l_buf))
     $   call errquit('ma_pop?',l_buf,0)
       end if
      end do
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c,0)
c
c     Write one symmetry block of G to disk
c
csh   if (.not. uccsdt_ampfile_add_t1(d_g, spini, g_g))
csh  $ call errquit('write_t1 failed',0)
      if (.not. ga_destroy(g_g)) call errquit('GA 999',0,0)
      if (odebug) call ga_print(g_g_big)
c
      return
      end
c $Id$
