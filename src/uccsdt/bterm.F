      subroutine uccsdt_bxgterms(d_t,d_b,d_x,d_g)
c $Id: bterm.F,v 1.2 2002-08-27 21:08:11 sohirata Exp $
c
c     Compute B, X, and G terms defined as follows:
c
c     Pure spin
c
c     B(a,b,i,j) = 1/2 <ab|^|cd> t(c,d,i,j)
c
c     ... t(u,v,i,j) = C(u,c)C(v,d)t(c,d,i,j)
c
c     ... Omega(u,v,i,j) = {<uv||wl> t(w,l,i,j)}/(1+delta_wl) ; w>=l
c
c     ... B(a,b,i,j) = Omega(u,v,i,j)C(u,a)C(v,b)
c
c     X(i,j,m,n) = <ij|^|ab> t(a,b,m,n)
c
c     ... X(i,j,m,n) = Omega(u,v,m,n)C(u,i)C(v,j)
c
c     G(e,m) = 1/2 <ie|^|ab> t(e,b,i,j)
c                    _    _      _   _
c                + <ei|^ ab> t(a,b,m,i)
c
c     ... G(e,m) = Omega(u,v,i,m)C(u,e)C(v,i) 
c                              _           _
c     ...        + Omega(u,v,m,i)C(u,e)C(v,i)
c
c     Mixed spin
c         _   _      _    _      _   _
c     B(a,b,i,j) = <ab|^ cd> t(c,d,i,j)
c                 _              _     _   _
c     ... t(u,v,i,j) = C(u,c)C(v,d)t(c,d,i,j)
c                      _            _            _
c     ... t[+/-](u,v,i,j) = t(u,v,i,j)+/-t(v,u,i,j)
c                     _                            _
c     ... Omega(u,v,i,j) = 1/2 {<uv||wl>[+]t(w,l,i,j)[+]
c                                                  _
c     ...                     + <uv||wl>[-]t(w,l,i,j)[-]}/(1+delta_wl) ; w>=l
c             _   _                _           _
c     ... B(a,b,i,j) = Omega(u,v,i,j)C(u,a)C(v,b)
c         _   _      _    _      _   _
c     X(i,j,m,n) = <ij|^ ab> t(a,b,m,n)
c             _   _                _           _
c     ... X(i,j,m,n) = Omega(u,v,m,n)C(u,i)C(v,j)
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "global.fh"
      integer d_t           ! [input]  DRA handle for t(a,b,i,j)
      integer d_b           ! [output] DRA handle for b(a,b,i,j)
      integer d_x           ! [output] DRA handle for x(a,b,i,j)
      integer d_g           ! [output] DRA handle for g(a,i)
      integer g_uv_p        ! GA handle for t[+](u,v,i,j)
      integer g_uv_m        ! GA handle for t[-](u,v,i,j)
      integer g_omega       ! GA handle for Omega(u,v,i,j)
      integer spini,spinj
c
c     Spin loop
c
      if (.not. uccsdt_ampfile_zero(d_g)) 
     $ call errquit('uccsdt_bxgterms: ampfile_zero?',999)
      do spini = 0,1
       do spinj = spini,1
c
c     SO half transform
c
        call uccsdt_tosopair(d_t,g_uv_p,g_uv_m,spini,spinj)
c
c     Contract
c
        call uccsdt_b_contract(g_uv_p,g_uv_m,g_omega,spini,spinj)
c
c     Make B
c
        call uccsdt_makeb(g_omega,d_b,spini,spinj)
c
c     Make X
c
        call uccsdt_makex(g_omega,d_x,spini,spinj)
c
c     Make G
c
        call uccsdt_makeg(g_omega,d_g,spini,spinj)
c
       end do
      end do
      return
      end
      subroutine uccsdt_tosopair(d_t,g_uv_p,g_uv_m,spini,spinj)
c
c     t(u,v,i,j) = Sum(a,b) C(u,a,spini)C(v,b,spinj)t(a,b,i,j)
c     where u & v are SO's (no spin), a,b,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
      integer d_t          ! [input] DRA handle for t(a,b,i,j)
      integer g_uv_p       ! [output] GA handle for t[+](u,v,i,j)
      integer g_uv_m       ! [output] GA handle for t[-](u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer spina, spinb ! Spins of a and b
      integer g_t2         ! Temporary GA t(a,b,i,j)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer k_p1,l_p1    ! Temporary MA for tranformation
      integer k_p2,l_p2    ! Temporary MA for tranformation
      integer k_m1,l_m1    ! Temporary MA for tranformation
      integer k_m2,l_m2    ! Temporary MA for tranformation
      integer lenij,lenub,lenuv
      integer lena,lenb,leni,lenu,lenv,lenia,leniu
      integer maxlenuv
      integer symij,symj,symi
      integer symub,symb,symu,symv,symuv
      integer symjb,symia,syma
      integer i,j,a,b,u,v,ia,iu,jb,ij
      integer joff(nw_max_nbf_small,0:7)
      integer boff(nw_max_nbf_small,0:7)
      integer symvoff(0:7,0:7)
      integer blo,bhi,alo,ilo,ihi,ulo,uhi,vlo
      integer jblo,jbhi,ijlo,ijhi
      double precision x,y
c
      spina = spini
      spinb = spinj
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         joff(j,symij) = lenij
         lenij = lenij + no_sym(symi,spini)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      do symub = 0,7
       lenub = 0
       do symb = 0,7
        symu = ieor(symub,symb)
        do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
         boff(b,symub) = lenub
         lenub = lenub + bf_per_ir(symu)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      maxlenuv = 0
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
c
c     Allocate g_uv_p and g_uv_m for t(u,b,i,j) and t[+/-](u,v,i,j)
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair', 
     $ maxlenuv,-1,g_uv_m)) call errquit
     $ ('uccsdt_tosopair: room for g_uv_m?',maxlenuv*lenij)
      if (spini .ne. spinj) then
       if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair', 
     $  maxlenuv,-1,g_uv_p)) call errquit
     $  ('uccsdt_tosopair: room for g_uv_p?',maxlenuv*lenij)
      end if     
c
c     Perform quarter transformation from t(a,b,i,j) to t(u,b,i,j)
c     t(i,a,syma,j,b,symb,symjb)
c
      do symjb = 0, 7
       symia = symjb
       do symb = 0, 7
        symj = ieor(symjb,symb)
c
c     Read one symmetry block of T2s from disk
c
        blo = v_sym(1,symb,spinb)
        bhi = v_sym(2,symb,spinb)
        if (.not. uccsdt_ampfile_read_t2(d_t, 
     $   spini, spina, spinj, spinb, symia, 
     $   blo, bhi, g_t2, .true., 'column'))
     $   call errquit('read_t2 failed',0)
c
c     Data-parallel transformation of a
c
        call ga_distribution(g_t2, ga_nodeid(), jblo, jbhi)
        jb = 0
        do b = blo, bhi
         do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
          jb = jb + 1
          if (jb.ge.jblo .and. jb.le.jbhi) then
           ia = 0
           iu = 0
c     Compare the next line with the comments in Amplitude.F
c     I think the comments are not accurate
           do syma = 0, 7
            symi = ieor(symia,syma)
            symij = ieor(symi,symj)
            symu = syma
            lena = nv_sym(syma,spina)
            leni = no_sym(symi,spini)
            lenu = bf_per_ir(symu)
            lenia = leni*lena
            leniu = leni*lenu
            alo = v_sym(1,syma,spina)
            ilo = o_sym(1,symi,spini)
            ihi = o_sym(2,symi,spini)
            ulo = bf_per_ir_cum(symu) + 1
            uhi = ulo + lenu - 1
            if (lenia.gt.0 .and. lenu.gt.0) then
             if (.not. ma_push_get(mt_dbl, lenia, 'buf',
     $        l_buf, k_buf)) 
     $        call errquit('ma? lenia',lenia)
             if (.not. ma_push_get(mt_dbl, leniu, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? leniu',leniu)
c
c     (ia) -> (iu) for a (virtual) & u (SO) of same irrep 
c
             call ga_get(g_t2,ia+1,ia+lenia,jb,jb,
     $        dbl_mb(k_buf),leni)
             call dgemm('n','t',lenu,leni,lena,
     $        1.0d0,c(ulo,alo,spina),nbf,
     $        dbl_mb(k_buf),leni,
     $        0.0d0,dbl_mb(k_tmp),lenu)
c     For c(), see cctrans3
             ij = joff(j,symij)
             do i = ilo,ihi
              ij = ij + 1
              call ga_put(g_uv_m,
     $         boff(b,symub)+1,boff(b,symub)+lenu,
     $         ij,ij,dbl_mb(k_tmp+lenu*(i-ilo)),lenu)
             end do
             ia = ia + lenia
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf)
            end if
           end do
          end if
         end do
        end do
        if (.not. ga_destroy(g_t2)) call errquit('GA 999',0)
       end do
      end do
c
c     Now have t(u,b,i,j) ... data parallel transform b in place 
c     (g_uv_m was allocated big enough to do this).
c
      call ga_distribution(g_uv_m, ga_nodeid(), ijlo, ijhi)
      ij = 0
      do symij = 0,7
       symuv = symij
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
          ij = ij + 1
          if (ij.ge.ijlo .and. ij.le.ijhi) then
           do symb = 0,7
            symu = ieor(symij,symb)
            symv = symb
            lenb = nv_sym(symb,spinb)
            lenu = bf_per_ir(symu)
            lenv = bf_per_ir(symv)
            lenub = lenu * lenb
            lenuv = lenu * lenv
            blo = v_sym(1,symb,spinb)
            vlo = bf_per_ir_cum(symv) + 1
            if (lenuv.gt.0 .and. lenb.gt.0) then
             if (.not. ma_push_get(mt_dbl, lenub, 'buf',
     $        l_buf, k_buf)) 
     $       call errquit('ma? lenub',lenub)
             if (.not. ma_push_get(mt_dbl, lenuv, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenuv',lenuv)
c
c     (ub) -> (uv) for b (virtual) & v (SO) of same irrep 
c
             call ga_get(g_uv_m,1,lenub,ij,ij,
     $        dbl_mb(k_buf),lenu)
             call dgemm('n','t',lenu,lenv,lenb,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        c(vlo,blo,spinb),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
c     For c(), see cctrans3
             call ga_put(g_uv_m,
     $        symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $        ij,ij,dbl_mb(k_tmp),lenu)
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
c     Symmetrize if spini .ne. spinj
c
      if (spini .ne. spinj) then
       call ga_copy(g_uv_m,g_uv_p)
       call ga_distribution(g_uv_m,ga_nodeid(),ijlo,ijhi)
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
c
c     symu >= symv to avoid referencing the already altered blocks
c
             if (symu .ge. symv) then
              lenu = bf_per_ir(symu)
              lenv = bf_per_ir(symv)
              lenuv = lenu * lenv
              if (lenuv .gt. 0) then
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'buf', l_buf, k_buf)) 
     $          call errquit('ma? lenuv',lenuv)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'tmp', l_tmp, k_tmp)) 
     $          call errquit('ma? lenuv',lenuv)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'plus', l_p1, k_p1)) 
     $          call errquit('ma? lenuv',lenuv)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'minus', l_m1, k_m1)) 
     $          call errquit('ma? lenuv',lenuv)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'plus', l_p2, k_p2)) 
     $          call errquit('ma? lenuv',lenuv)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'minus', l_m2, k_m2)) 
     $          call errquit('ma? lenuv',lenuv)
               call ga_get(g_uv_m,
     $          symvoff(symu,symuv)+1,symvoff(symu,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_tmp),lenv)
               call ga_get(g_uv_m,
     $          symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_buf),lenu)
               do v = 1,lenv
                do u = 1,lenu
                 x = dbl_mb(k_buf+(v-1)*lenu+u-1)
                 y = dbl_mb(k_buf+(u-1)*lenv+v-1)
                 dbl_mb(k_p1+(v-1)*lenu+u-1) = 0.5d0*(x+y)
                 dbl_mb(k_m1+(v-1)*lenu+u-1) = 0.5d0*(x-y)
                 dbl_mb(k_p2+(u-1)*lenv+v-1) = 0.5d0*(x+y)
                 dbl_mb(k_m2+(u-1)*lenv+v-1) = 0.5d0*(x-y)
                end do
               end do
c
c     we fill the whole g_uv_p and g_uv_m, since the u>=v restriction
c     in the contraction step is based on tiling rather than sym blocks
c
               call ga_put(g_uv_p,
     $          symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_p1),lenu)
               call ga_put(g_uv_m,
     $          symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_m1),lenu)
               call ga_put(g_uv_p,
     $          symvoff(symu,symuv)+1,symvoff(symu,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_p2),lenv)
               call ga_put(g_uv_m,
     $          symvoff(symu,symuv)+1,symvoff(symu,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_m2),lenv)
               if (.not. ma_pop_stack(l_m2))
     $          call errquit('ma_pop?',l_m2)
               if (.not. ma_pop_stack(l_p2))
     $          call errquit('ma_pop?',l_p2)
               if (.not. ma_pop_stack(l_m1))
     $          call errquit('ma_pop?',l_m1)
               if (.not. ma_pop_stack(l_p1))
     $          call errquit('ma_pop?',l_p1)
               if (.not. ma_pop_stack(l_tmp))
     $          call errquit('ma_pop?',l_tmp)
               if (.not. ma_pop_stack(l_buf))
     $          call errquit('ma_pop?',l_buf)
              end if
             end if
            end do
           end if
          end do
         end do
        end do
       end do
      end if
c
      return
      end
      subroutine uccsdt_b_contract(g_uv_p,g_uv_m,g_omega,
     $ spini,spinj)
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
c
c     Perform contraction 
c
c     for u>=v, i,j, w>=l
c     Omega(u,v,i,j) = {<uv||wl>[-] t[-](w,l,i,j)               (puer spin/mixed spin)
c                    +  <uv||wl>[+] t[+](w,l,i,j)}/(1+delta_wl) (mixed spin)
c
c     Omega(u,v,i,j) = 1/2 Omega(u,v,i,j) (mixed spin)
c
      integer g_uv_p          ! [input] GA handle for t[+](w,l,i,j)
      integer g_uv_m          ! [input] GA handle for t[-](w,l,i,j)
      integer g_omega         ! [output] GA handle for Omega(u,v,i,j)
      integer spini, spinj    ! [input] Spins of i and j
      integer l_gx,k_gx       ! Temporary MA for integrals
      integer l_gy,k_gy       ! Temporary MA for integrals
      integer l_omega,k_omega ! Temporary MA for Omega
      integer l_t,k_t         ! Temporary MA for t
      integer l_m,k_m         ! Temporary MA for Omega[-]
      integer l_p,k_p         ! Temporary MA for Omega[+]
      integer lenij,lenwl,lenlv,lenwu,lenuv
      integer lenu,lenv,lenl,lenw,lenwv,lenlu
      integer lengx,lengy
      integer symij,symj,symi,symwl,syml,symw
      integer symuv,symu,symv,symlv,symwu
      integer symjoff(0:7,0:7)
      integer symvoff(0:7,0:7)
      integer loff2(nw_max_nbf_small,0:7)
      integer maxlenuv
      integer sh_blk_lo,sh_blk_hi
      integer nsh_block
      integer sh_block(2,nw_max_shells_small)
      integer count,next,task_size
      integer ush_blk,vsh_blk,wsh_blk,lsh_blk
      integer urange(2,0:7),vrange(2,0:7),wrange(2,0:7),lrange(2,0:7)
      integer ush_blk_lo,ush_blk_hi
      integer vsh_blk_lo,vsh_blk_hi
      integer wsh_blk_lo,wsh_blk_hi
      integer lsh_blk_lo,lsh_blk_hi
      integer uoff(nw_max_nbf_small,0:7)
      integer voff(nw_max_nbf_small,0:7)
      integer woff(nw_max_nbf_small,0:7)
      integer loff(nw_max_nbf_small)
      integer symwuoff(0:7)
      integer wlo,whi,llo,lhi
      integer u,v,w,l,uvwl
      integer nxtask
      external nxtask
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        symjoff(symj,symij) = lenij
        lenij = lenij + no_sym(symi,spini)*no_sym(symj,spinj)
       end do
      end do
c
c     Addressing & sizes
c
      maxlenuv = 0
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
c
c     Addressing & sizes
c
      do symwl = 0,7
       lenwl = 0
       do syml = 0,7
        symw = ieor(symwl,syml)
        do l = bf_per_ir_cum(syml)+1,
     $         bf_per_ir_cum(syml)+bf_per_ir(syml)
         loff2(l,symwl) = lenwl
         lenwl = lenwl + bf_per_ir(symw)
        end do
       end do
      end do
c
c     Allocate g_omega_p and g_omega_m for Omega[+/-](u,v,i,j)
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair', 
     $ maxlenuv,-1,g_omega)) call errquit
     $ ('uccsdt_tosopair: room for tmp?',maxlenuv*lenij)
      call ga_zero(g_omega)
c     
c     Block the shells so that there are roughly no more than XX
c     functions in each block (better would be to look at the memory
c     available)
c     
      sh_blk_lo = 1
      sh_blk_hi = sh_n_uq
      nproc = ga_nnodes()
      task_size = sh_n_uq*(sh_n_uq+1)/2
      task_size = task_size*(task_size+1)/2
      task_size = max(1,task_size/(20*nproc))
c
      call cctrans_block(nsh_block,sh_block,sh_blk_lo,sh_blk_hi)
c
      count = 0
      next = nxtask(nproc,task_size)
      do ush_blk = 1, nsh_block
       ush_blk_lo = sh_block(1,ush_blk)
       ush_blk_hi = sh_block(2,ush_blk)
       call soints_shblock_info(basis,ush_blk_lo,ush_blk_hi,urange)
       do vsh_blk = 1, ush_blk ! (V<=U)
        vsh_blk_lo = sh_block(1,vsh_blk)
        vsh_blk_hi = sh_block(2,vsh_blk)
        call soints_shblock_info(basis,vsh_blk_lo,vsh_blk_hi,vrange)
        do wsh_blk = 1, ush_blk ! (W,L <= U,V)
         wsh_blk_lo = sh_block(1,wsh_blk)
         wsh_blk_hi = sh_block(2,wsh_blk)
         call soints_shblock_info(basis,wsh_blk_lo,wsh_blk_hi,wrange)
         do lsh_blk = 1, wsh_blk ! (L<=W)
          lsh_blk_lo = sh_block(1,lsh_blk)
          lsh_blk_hi = sh_block(2,lsh_blk)
          call soints_shblock_info(basis,lsh_blk_lo,lsh_blk_hi,lrange)
          if (next .eq. count) then
c     
c     Compute g(l,v,w,u) = (lv|wu) = <lw|vu> 
c     
           call soints_block_info(basis, nbf,
     $      lsh_blk_lo, lsh_blk_hi, loff, lrange,
     $      vsh_blk_lo, vsh_blk_hi, voff, vrange,
     $      wsh_blk_lo, wsh_blk_hi, woff, wrange,
     $      ush_blk_lo, ush_blk_hi, uoff, urange,
     $      lenlv, lenwu, symwuoff, lengx)
c     
           if (.not. ma_push_get
     $      (mt_dbl, lengx, 'b_contract', l_gx, k_gx))
     $      call errquit('b_contract: no room for g', lengx)
c     
           call soints_block(basis, nbf,
     $      lsh_blk_lo, lsh_blk_hi, loff, lrange,
     $      vsh_blk_lo, vsh_blk_hi, voff, vrange,
     $      wsh_blk_lo, wsh_blk_hi, woff, wrange,
     $      ush_blk_lo, ush_blk_hi, uoff, urange,
     $      1d-12, dbl_mb(k_gx))
c     
c     Compute g(l,v,w,u) = (wv|lu) = <wl|vu> 
c     
           call soints_block_info(basis, nbf,
     $      wsh_blk_lo, wsh_blk_hi, woff, wrange,
     $      vsh_blk_lo, vsh_blk_hi, voff, vrange,
     $      lsh_blk_lo, lsh_blk_hi, loff, lrange,
     $      ush_blk_lo, ush_blk_hi, uoff, urange,
     $      lenwv, lenlu, symwuoff, lengy)
c     
           if (.not. ma_push_get
     $      (mt_dbl, lengy, 'b_contract', l_gy, k_gy))
     $      call errquit('b_contract: no room for g', lengy)
c     
           call soints_block(basis, nbf,
     $      wsh_blk_lo, wsh_blk_hi, woff, wrange,
     $      vsh_blk_lo, vsh_blk_hi, voff, vrange,
     $      lsh_blk_lo, lsh_blk_hi, loff, lrange,
     $      ush_blk_lo, ush_blk_hi, uoff, urange,
     $      1d-12, dbl_mb(k_gy))
c
c     Contract
c
           do symuv = 0,7
            symij = symuv
            symwl = symuv
            lenuv = 0
            do symv = 0,7
             symu = ieor(symuv,symv)
             symvoff(symv,symuv) = lenuv
             lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
            end do
            do symj = 0,7
             symi = ieor(symij,symj)
             lenij = no_sym(symi,spini) * no_sym(symj,spinj)
             if (lenij .gt. 0) then
              do symv = 0,7
               symu = ieor(symuv,symv)
               lenv = vrange(2,symv) - vrange(1,symv) + 1
               lenu = urange(2,symu) - urange(1,symu) + 1
               lenuv = lenu * lenv
               if (lenuv .gt. 0) then
                if (.not. ma_push_get(mt_dbl, 
     $           lenuv*lenij, 'omega', l_omega, k_omega)) 
     $          call errquit('ma? lenuv*lenij',lenuv*lenij)
                do syml = 0,7
                 symw = ieor(symwl,syml)
                 symwu = ieor(symw,symu)
                 symlv = ieor(syml,symv)
                 lenl = lrange(2,syml) - lrange(1,syml) + 1
                 lenw = wrange(2,symw) - wrange(1,symw) + 1
                 lenwl = lenw * lenl
                 wlo = wrange(1,symw) - bf_per_ir_cum(symw)
                 whi = wrange(2,symw) - bf_per_ir_cum(symw)
                 llo = lrange(1,syml) - bf_per_ir_cum(syml)
                 lhi = lrange(2,syml) - bf_per_ir_cum(syml)
                 if (lenwl .gt. 0) then
                  if (.not. ma_push_get(mt_dbl, 
     $             lenwl*lenij, 't', l_t, k_t)) 
     $             call errquit('ma? lenwl',lenwl*lenij)
                  do l = lrange(1,syml),lrange(2,syml)
                   call ga_get(g_uv_m,
     $             wlo+loff2(l,symwl),whi+loff2(l,symwl), 
     $             symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $             dbl_mb(k_t),lenwl)
                  end do
                  if (.not. ma_push_get(mt_dbl, 
     $             lenuv*lenwl, 'minus', l_m, k_m)) 
     $             call errquit('ma? lenuv*lenwl',lenuv*lenwl)
                  uvwl = 0
                  do l = lrange(1,syml), lrange(2,syml)
                   do w = wrange(1,symw), wrange(2,symw)
                    do v = vrange(1,symv), vrange(2,symv)
                     do u = urange(1,symu), urange(2,symu)
                      uvwl = uvwl + 1
                      dbl_mb(k_m+uvwl) = 
     $                dbl_mb(k_gx+loff(l)+voff(v,symlv)
     $                 +woff(w,symwu)+uoff(u,symwu))
     $                -dbl_mb(k_gy+loff(l)+voff(v,symlv)
     $                 +woff(w,symwu)+uoff(u,symwu))
                     end do
                    end do
                   end do
                  end do
                  call dgemm('n','n',lenuv,lenij,lenwl,
     $             1.0d0,dbl_mb(k_m),lenuv,
     $             dbl_mb(k_t),lenwl,
     $             0.0d0,dbl_mb(k_omega),lenuv)
                  if (symw .eq. syml) then
                   do l = lrange(1,syml),lrange(2,syml)
                    call ga_acc(g_omega,
     $              wlo+loff2(l,symwl),whi+loff2(l,symwl), 
     $              symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $              dbl_mb(k_omega),lenuv,0.5d0)
                   end do
                  else
                   do l = lrange(1,syml),lrange(2,syml)
                    call ga_acc(g_omega,
     $              wlo+loff2(l,symwl),whi+loff2(l,symwl), 
     $              symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $              dbl_mb(k_omega),lenuv,1.0d0)
                   end do
                   do w = wrange(1,symw),wrange(2,symw)
                    call ga_acc(g_omega,
     $              llo+loff2(w,symwl),lhi+loff2(w,symwl), 
     $              symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $              dbl_mb(k_omega),lenuv,1.0d0)
                   end do
                  end if
                  if (.not. ma_pop_stack(l_m))
     $             call errquit('ma_pop?',l_m)
                  if (spini .ne. spinj) then
                   if (.not. ma_push_get(mt_dbl, 
     $              lenuv*lenwl, 'plus', l_p, k_p)) 
     $              call errquit('ma? lenuv*lenwl',lenuv*lenwl)
                   uvwl = 0
                   do l = lrange(1,syml), lrange(2,syml)
                    do w = wrange(1,symw), wrange(2,symw)
                     do v = vrange(1,symv), vrange(2,symv)
                      do u = urange(1,symu), urange(2,symu)
                       uvwl = uvwl + 1
                       dbl_mb(k_p+uvwl) = 
     $                 dbl_mb(k_gx+loff(l)+voff(v,symlv)
     $                  +woff(w,symwu)+uoff(u,symwu))
     $                 +dbl_mb(k_gy+loff(l)+voff(v,symlv)
     $                  +woff(w,symwu)+uoff(u,symwu))
                      end do
                     end do
                    end do
                   end do
                   call dgemm('n','n',lenuv,lenij,lenwl,
     $              1.0d0,dbl_mb(k_p),lenuv,
     $              dbl_mb(k_t),lenwl,
     $              0.0d0,dbl_mb(k_omega),lenuv)
                   if (symw .eq. syml) then
                    do l = lrange(1,syml),lrange(2,syml)
                     call ga_acc(g_omega,
     $               wlo+loff2(l,symwl),whi+loff2(l,symwl), 
     $               symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $               dbl_mb(k_omega),lenuv,0.5d0)
                    end do
                   else
                    do l = lrange(1,syml),lrange(2,syml)
                     call ga_acc(g_omega,
     $               wlo+loff2(l,symwl),whi+loff2(l,symwl), 
     $               symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $               dbl_mb(k_omega),lenuv,1.0d0)
                    end do
                    do w = wrange(1,symw),wrange(2,symw)
                     call ga_acc(g_omega,
     $               llo+loff2(w,symwl),lhi+loff2(w,symwl), 
     $               symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $               dbl_mb(k_omega),lenuv,1.0d0)
                    end do
                   end if
                  end if
                  if (.not. ma_pop_stack(l_t))
     $             call errquit('ma_pop?',l_t)
                 end if
                end do
                if (.not. ma_pop_stack(l_omega))
     $           call errquit('ma_pop?',l_omega)
               end if
              end do
             end if
            end do
           end do
c
           if (.not. ma_pop_stack(l_gy))
     $      call errquit('ma_pop?',l_gy)
           if (.not. ma_pop_stack(l_gx))
     $      call errquit('ma_pop?',l_gx)
           next = nxtask(nproc,task_size)
          end if
          count = count + 1
         end do
        end do
       end do
      end do
c
c     Barrier
c
      next = nxtval(-nproc)
      call ga_sync()
      if (spini .ne. spinj) call ga_scale(g_omega,0.5d0)
c
      return
      end
      subroutine uccsdt_makeb(g_omega,d_b,spini,spinj)
c
c     b(a,b,i,j) = Sum(u,v) t(u,v,i,j)C(u,a,spini)C(v,b,spinj)
c     where u & v are SO's (no spin), a,b,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
      integer d_b          ! [output] DRA handle for b(a,b,i,j)
      integer g_omega      ! [input] GA handle for Omega(u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer spina, spinb ! Spins of a and b
      integer g_tmp        ! Temporary GA t(u,b,i,j) & t(u,v,i,j)
      integer g_b          ! Temporary GA b(a,b,i,j)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer lenij,lenub,lenuv,lenu,leniu
      integer lenv,lenb,lenia,lenjb,leni,lena
      integer maxlenuv
      integer symij,symj,symi,symub,symb
      integer symuv,symu,symv,symjb,symia,syma
      integer joff(nw_max_nbf_small,0:7)
      integer boff(nw_max_nbf_small,0:7)
      integer symvoff(0:7,0:7)
      integer symboff(0:7,0:7)
      integer ij,jb,ia,iu
      integer ijlo,ijhi,jblo,jbhi
      integer i,j,b
      integer blo,bhi,vlo,alo,ilo,ihi,ulo,uhi
c
      spina = spini
      spinb = spinj
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         joff(j,symij) = lenij
         lenij = lenij + no_sym(symi,spini)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      do symub = 0,7
       lenub = 0
       do symb = 0,7
        symu = ieor(symub,symb)
        do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
         boff(b,symub) = lenub
         lenub = lenub + bf_per_ir(symu)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      do symub = 0,7
       lenub = 0
       do symb = 0,7
        symu = ieor(symub,symb)
        symboff(symb,symub) = lenub
        lenub = bf_per_ir(symu)*nv_sym(symb,spinb)
       end do
      end do
c
c     Addressing & sizes
c
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
c
c     Allocate temporary GA
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'tmp', 
     $ maxlenuv,-1,g_tmp)) call errquit
     $ ('uccsdt_makeb: room for tmp?',maxlenuv*lenij)
c
c     Data parallel transform Omega(u,v,i,j)C(v,b) = tmp(u,b,i,j)
c     (tmp was allocated big enough to do this).
c
      call ga_distribution(g_omega, ga_nodeid(), ijlo, ijhi)
      ij = 0
      do symij = 0,7
       symuv = symij
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
          ij = ij + 1
          if (ij.ge.ijlo .and. ij.le.ijhi) then
           do symb = 0,7
            symu = ieor(symij,symb)
            symv = symb
            lenb = nv_sym(symb,spinb)
            lenu = bf_per_ir(symu)
            lenv = bf_per_ir(symv)
            lenub = lenu * lenb
            lenuv = lenu * lenv
            blo = v_sym(1,symb,spinb)
            vlo = bf_per_ir_cum(symv) + 1
            if (lenuv.gt.0 .and. lenb.gt.0) then
             if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $        l_buf, k_buf)) 
     $        call errquit('ma? lenub',lenub)
             if (.not. ma_push_get(mt_dbl, lenub, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenuv',lenuv)
c
c     (uv) -> (ub) for b (virtual) & v (SO) of same irrep 
c
             call ga_get(g_omega,1,lenuv,ij,ij,
     $        dbl_mb(k_buf),lenu)
             call dgemm('n','n',lenu,lenb,lenv,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        c(vlo,blo,spinb),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
c     For c(), see cctrans3
             call ga_put(g_tmp,
     $        symboff(symb,symub)+1,symboff(symb,symub)+lenub,
     $        ij,ij,dbl_mb(k_tmp),lenu)
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
c     Perform quarter transformation from tmp(u,b,i,j) to t(a,b,i,j)
c     t(i,a,syma,j,b,symb,symjb)
c
      do symjb = 0,7
       symia = symjb
       do symb = 0,7
        symj = ieor(symjb,symb)
c
c     Allocate B block
c
        blo = v_sym(1,symb,spinb)
        bhi = v_sym(2,symb,spinb)
        lenia = 0
        do syma = 0,7
         symi = ieor(symia,syma)
         lenia = lenia + no_sym(symi,spini) * nv_sym(syma,spina)
        end do
        lenjb = (bhi-blo+1)*no_sym(symj,spinj)
        if (.not.ga_create(mt_dbl,lenia,lenjb,'b', 
     $   lenia,-1,g_b)) call errquit
     $   ('uccsdt_makeb: room for b?',lenia*lenjb)
c
c     Data-parallel transformation of a
c
        call ga_distribution(g_tmp, ga_nodeid(), jblo, jbhi)
        jb = 0
        do b = blo, bhi
         do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
          jb = jb + 1
          if (jb.ge.jblo .and. jb.le.jbhi) then
           ia = 0
           iu = 0
c     Compare the next line with the comments in Amplitude.F
c     I think the comments are not accurate
           do syma = 0, 7
            symi = ieor(symia,syma)
            symij = ieor(symi,symj)
            symu = syma
            lena = nv_sym(syma,spina)
            leni = no_sym(symi,spini)
            lenu = bf_per_ir(symu)
            lenia = leni*lena
            leniu = leni*lenu
            alo = v_sym(1,syma,spina)
            ilo = o_sym(1,symi,spini)
            ihi = o_sym(2,symi,spini)
            ulo = bf_per_ir_cum(symu) + 1
            uhi = ulo + lenu - 1
            if (lenia.gt.0 .and. lenu.gt.0) then
             if (.not. ma_push_get(mt_dbl, leniu, 'buf',
     $        l_buf, k_buf)) 
     $        call errquit('ma? leniu',leniu)
             if (.not. ma_push_get(mt_dbl, lenia, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenia',lenia)
c
c     (iu) -> (ia) for a (virtual) & u (SO) of same irrep 
c
             ij = joff(j,symij)
             do i = ilo,ihi
              ij = ij + 1
              call ga_get(g_tmp,
     $         boff(b,symub)+1,boff(b,symub)+lenu,
     $         ij,ij,dbl_mb(k_buf+lenu*(i-ilo)),lenu)
             end do
             call dgemm('t','t',lena,leni,lenu,
     $        1.0d0,c(ulo,alo,spina),nbf,
     $        dbl_mb(k_tmp),leni,
     $        0.0d0,dbl_mb(k_buf),lenu)
c     For c(), see cctrans3
             call ga_put(g_b,ia+1,ia+lenia,jb,jb,
     $        dbl_mb(k_buf),leni)
             ia = ia + lenia
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf)
            end if
           end do
          end if
         end do
        end do
c
c     Write one symmetry block of B to disk
c
        if (.not. uccsdt_ampfile_write_t2(d_b, 
     $   spini, spina, spinj, spinb, symia, 
     $   blo, bhi, g_b))
     $   call errquit('write_t2 failed',0)
        if (.not. ga_destroy(g_b)) call errquit('GA 999',0)
       end do
      end do
      if (.not. ga_destroy(g_tmp)) call errquit('GA 999',0)
c
      return
      end
      subroutine uccsdt_makex(g_omega,d_x,spini,spinj)
c
c     x(k,l,i,j) = Sum(u,v) t(u,v,i,j)C(u,k,spini)C(v,l,spinj)
c     where u & v are SO's (no spin), k,b,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
      integer d_x          ! [output] DRA handle for X(k,l,i,j)
      integer g_omega      ! [input] GA handle for Omega(u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer spink, spinl ! Spins of k and l
      integer g_tmp        ! Temporary GA t(u,l,i,j) & t(u,v,i,j)
      integer g_x          ! Temporary GA x(k,l,i,j)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer lenij,lenul,lenuv,lenu,lenv,lenl,lenik,lenjl,leni,lenk
      integer maxlenuv
      integer symij,symj,symi,symul,syml
      integer symuv,symu,symv,symjl,symik,symk
      integer joff(nw_max_nbf_small,0:7)
      integer loff(nw_max_nbf_small,0:7)
      integer symvoff(0:7,0:7)
      integer symloff(0:7,0:7)
      integer ij,jl,ik,iu
      integer ijlo,ijhi,jllo,jlhi
      integer i,j,l
      integer llo,lhi,vlo,klo,ilo,ihi,ulo,uhi
c
      spink = spini
      spinl = spinj
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         joff(j,symij) = lenij
         lenij = lenij + no_sym(symi,spini)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      do symul = 0,7
       lenul = 0
       do syml = 0,7
        symu = ieor(symul,syml)
        do l = o_sym(1,syml,spinl),o_sym(2,syml,spinl)
         loff(l,symul) = lenul
         lenul = lenul + bf_per_ir(symu)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      do symul = 0,7
       lenul = 0
       do syml = 0,7
        symu = ieor(symul,syml)
        symloff(syml,symul) = lenul
        lenul = bf_per_ir(symu)*no_sym(syml,spinl)
       end do
      end do
c
c     Addressing & sizes
c
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
c
c     Allocate temporary GA
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair', 
     $ maxlenuv,-1,g_tmp)) call errquit
     $ ('uccsdt_makex: room for tmp?',maxlenuv*lenij)
c
c     Data parallel transform Omega(u,v,i,j)C(v,l) = tmp(u,l,i,j)
c     (tmp was allocated big enough to do this).
c
      call ga_distribution(g_omega, ga_nodeid(), ijlo, ijhi)
      ij = 0
      do symij = 0,7
       symuv = symij
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
          ij = ij + 1
          if (ij.ge.ijlo .and. ij.le.ijhi) then
           do syml = 0,7
            symu = ieor(symij,syml)
            symv = syml
            lenl = no_sym(syml,spinl)
            lenu = bf_per_ir(symu)
            lenv = bf_per_ir(symv)
            lenul = lenu * lenl
            lenuv = lenu * lenv
            llo = o_sym(1,syml,spinl)
            vlo = bf_per_ir_cum(symv) + 1
            if (lenuv.gt.0 .and. lenl.gt.0) then
             if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $        l_buf, k_buf)) 
     $       call errquit('ma? lenul',lenul)
             if (.not. ma_push_get(mt_dbl, lenul, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenuv',lenuv)
c
c     (uv) -> (ul) for l (virtual) & v (SO) of same irrep 
c
             call ga_get(g_omega,1,lenuv,ij,ij,
     $        dbl_mb(k_buf),lenu)
             call dgemm('n','n',lenu,lenl,lenv,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        c(vlo,llo,spinl),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
c     For c(), see cctrans3
             call ga_put(g_tmp,
     $        symloff(syml,symul)+1,symloff(syml,symul)+lenul,
     $        ij,ij,dbl_mb(k_tmp),lenu)
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
c     Perform quarter transformation from tmp(u,l,i,j) to t(k,l,i,j)
c     t(i,k,symk,j,l,syml,symjl)
c
      do symjl = 0,7
       symik = symjl
       do syml = 0,7
        symj = ieor(symjl,syml)
c
c     Allocate X block
c
        llo = o_sym(1,syml,spinl)
        lhi = o_sym(2,syml,spinl)
        lenik = 0
        do symk = 0,7
         symi = ieor(symik,symk)
         lenik = lenik + no_sym(symi,spini) * no_sym(symk,spink)
        end do
        lenjl = (lhi-llo+1)*no_sym(symj,spinj)
        if (.not.ga_create(mt_dbl,lenik,lenjl,'l', 
     $   lenik,-1,g_x)) call errquit
     $   ('uccsdt_makex: room for x?',lenik*lenjl)
c
c     Data-parallel transformation of k
c
        call ga_distribution(g_tmp, ga_nodeid(), jllo, jlhi)
        jl = 0
        do l = llo, lhi
         do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
          jl = jl + 1
          if (jl.ge.jllo .and. jl.le.jlhi) then
           ik = 0
           iu = 0
c     Compare the next line with the comments in Amplitude.F
c     I think the comments are not accurate
           do symk = 0, 7
            symi = ieor(symik,symk)
            symij = ieor(symi,symj)
            symu = symk
            lenk = no_sym(symk,spink)
            leni = no_sym(symi,spini)
            lenu = bf_per_ir(symu)
            lenik = leni*lenk
            leniu = leni*lenu
            klo = o_sym(1,symk,spink)
            ilo = o_sym(1,symi,spini)
            ihi = o_sym(2,symi,spini)
            ulo = bf_per_ir_cum(symu) + 1
            uhi = ulo + lenu - 1
            if (lenik.gt.0 .and. lenu.gt.0) then
             if (.not. ma_push_get(mt_dbl, leniu, 'buf',
     $        l_buf, k_buf)) 
     $        call errquit('ma? lenia',leniu)
             if (.not. ma_push_get(mt_dbl, lenik, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? leniu',lenik)
c
c     (iu) -> (ik) for k (occupied) & u (SO) of same irrep 
c
             ij = joff(j,symij)
             do i = ilo,ihi
              ij = ij + 1
              call ga_get(g_tmp,
     $         loff(l,symul)+1,loff(l,symul)+lenu,
     $         ij,ij,dbl_mb(k_buf+lenu*(i-ilo)),lenu)
             end do
             call dgemm('t','t',lenk,leni,lenu,
     $        1.0d0,c(ulo,klo,spink),nbf,
     $        dbl_mb(k_tmp),leni,
     $        0.0d0,dbl_mb(k_buf),lenu)
c     For c(), see cctrans3
             call ga_put(g_t2,ik+1,ik+lenik,jl,jl,
     $        dbl_mb(k_buf),leni)
             ik = ik + lenik
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf)
            end if
           end do
          end if
         end do
        end do
c
c     Write one symmetry block of X to disk
c
        if (.not. uccsdt_ampfile_write_t2(d_b, 
     $   spini, spink, spinj, spinl, symik, 
     $   llo, lhi, g_x))
     $   call errquit('write_t2 failed',0)
        if (.not. ga_destroy(g_x)) call errquit('GA 999',0)
       end do
      end do
      if (.not. ga_destroy(g_tmp)) call errquit('GA 999',0)
c
      return
      end
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
      integer lenui,lenu,lenv,lenuv,lena,leni,lenai
      integer symui,symi,symu,symij,symuv,symj,symuj,syma
      integer ioff(nw_max_nbf_small,0:7)
      integer ij,ijlo,ijhi
      integer vlo,alo
      integer i,j
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
c     Data parallel transform 
c     Pure spin:  Omega(u,v,i,j)C(v,i) = tmp(u,j)
c     Mixed spin: Omega(u,v,i,j)C(v,j) = tmp(u,i)
c
      call ga_distribution(g_omega, ga_nodeid(), ijlo, ijhi)
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
              call dgemm('n','n',lenu,1,lenv,
     $         1.0d0,dbl_mb(k_buf),lenu,
     $         c(vlo,i,spini),nbf,
     $         0.0d0,dbl_mb(k_tmp),lenu)
              call ga_acc(g_g,
     $         ioff(j,symuj)+1,ioff(j,symuj)+lenu,
     $         1,1,dbl_mb(k_tmp),lenu)
             else
              call dgemm('n','n',lenu,1,lenv,
     $         1.0d0,dbl_mb(k_buf),lenu,
     $         c(vlo,j,spinj),nbf,
     $         0.0d0,dbl_mb(k_tmp),lenu)
              call ga_acc(g_g,
     $         ioff(i,symui)+1,ioff(i,symui)+lenu,
     $         1,1,dbl_mb(k_tmp),lenu)
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
     $   1.0d0,c(vlo,alo,spini),nbf,
     $   dbl_mb(k_tmp),lenu,
     $   0.0d0,dbl_mb(k_buf),lena)
        call ga_put(g_g,
     $   ioff(i,symui)+1,ioff(i,symui)+lenai,
     $   1,1,dbl_mb(k_buf),lena)
       end do
      end do
c
c     Write one symmetry block of G to disk
c
      if (.not. uccsdt_ampfile_add_t1(d_g, spini, g_g))
     $ call errquit('write_t1 failed',0)
      if (.not. ga_destroy(g_g)) call errquit('GA 999',0)
c
      return
      end
