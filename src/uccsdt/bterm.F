      subroutine shmocoeff(sh_nmo,sh_nbf,pa,ha,pb,hb)
      implicit none
#include "cuccsdtP.fh"
      integer sh_nmo,sh_nbf
      double precision pa(sh_nmo,sh_nbf)
      double precision ha(sh_nmo,sh_nbf)
      double precision pb(sh_nmo,sh_nbf)
      double precision hb(sh_nmo,sh_nbf)
      double precision buffer
      integer mo,bf
c
      do mo = 1,sh_nmo
       do bf = 1,sh_nbf
        buffer = pa(mo,bf)
        call ga_put(g_part(1),bf,bf,mo,mo,buffer,1) 
        buffer = ha(mo,bf)
        call ga_put(g_hole(1),bf,bf,mo,mo,buffer,1) 
        buffer = pb(mo,bf)
        call ga_put(g_part(2),bf,bf,mo,mo,buffer,1) 
        buffer = hb(mo,bf)
        call ga_put(g_hole(2),bf,bf,mo,mo,buffer,1) 
       end do
      end do
csh   call ga_print(g_part(1))
csh   call ga_print(g_hole(1))
csh   call ga_print(g_part(2))
csh   call ga_print(g_hole(2))
c
      return
      end
      subroutine cctrans_block(nvblock, vblock, vsh_blk_lo, vsh_blk_hi)
      implicit none
#include "cuccsdtP.fh"
      integer nvblock, vblock(2,*), vsh_blk_lo, vsh_blk_hi
c
      integer ninblock, lo, vsh
c
      nvblock = 0
      ninblock = 0
      lo = vsh_blk_lo
      do vsh = vsh_blk_lo, vsh_blk_hi
         ninblock = ninblock + sh_uq_bf(2,vsh) - sh_uq_bf(1,vsh) + 1
         if ( (ninblock.ge.32) .or. (vsh.eq.vsh_blk_hi) ) then
            ninblock = 0
            nvblock = nvblock + 1
            vblock(1,nvblock) = lo
            vblock(2,nvblock) = vsh
            lo = vsh+1
         end if
      end do
c
      return
      end
      subroutine uccsdt_tosopair(g_t2_big,g_uv_p,g_uv_m,spini,spinj)
csh   subroutine uccsdt_tosopair(d_t,g_uv_p,g_uv_m,spini,spinj)
c
c     t(u,v,i,j) = Sum(a,b) C(u,a,spini)C(v,b,spinj)t(a,b,i,j)
c     where u & v are SO's (no spin), a,b,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
      integer g_t2_big     ! [input] DRA handle for t(a,b,i,j)
csh   integer d_t          ! [input] DRA handle for t(a,b,i,j)
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
      integer k_c,l_c      ! Temporary MA for MO coeffs
      integer lenij,lenub,lenuv,lenj
      integer lena,lenb,leni,lenu,lenv,lenia,leniu
      integer maxlenuv
      integer symij,symj,symi
      integer symub,symb,symu,symv,symuv
      integer symjb,symia,syma
      integer i,j,b,u,v,ia,iu,jb,ij
      integer joff(nw_max_nbf_small,0:7)
      integer boff(nw_max_nbf_small,0:7)
      integer symvoff(0:7,0:7)
      integer blo,bhi,alo,ilo,ihi,ulo,uhi,vlo
      integer jblo,jbhi,ijlo,ijhi
      integer dummy
      double precision x,y
csh
c     logical int_chk_init
c     external int_chk_init
c     logical temp
c     temp = int_chk_init('b term sopair')
csh 
c
      spina = spini
      spinb = spinj
csh
c     write(*,*) 'no = ',no(1),no(2)
c     write(*,*) 'nv = ',nv(1),nv(2)
c     write(*,*) 'no_sym = ',(no_sym(i,1),i=0,7),(no_sym(i,2),i=0,7)
c     write(*,*) 'nv_sym = ',(nv_sym(i,1),i=0,7),(nv_sym(i,2),i=0,7)
c     write(*,*) 'o_sym(a) = ',(o_sym(1,i,1),i=0,7),(o_sym(2,i,1),i=0,7)
c     write(*,*) 'o_sym(b) = ',(o_sym(1,i,2),i=0,7),(o_sym(2,i,2),i=0,7)
c     write(*,*) 'v_sym(a) = ',(v_sym(1,i,1),i=0,7),(v_sym(2,i,1),i=0,7)
c     write(*,*) 'v_sym(b) = ',(v_sym(1,i,2),i=0,7),(v_sym(2,i,2),i=0,7)
c     write(*,*) 'bf_per_ir = ',(bf_per_ir(i),i=0,7)
csh
c
c     Addressing & sizes
c
c     write(*,*) "j-offsets"
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         joff(j,symij) = lenij
         lenij = lenij + no_sym(symi,spini)
c        write(*,*) j,symij,joff(j,symij)
        end do
       end do
      end do
c     write(*,*) "lenij = ",lenij
c
c     Addressing & sizes
c
c     write(*,*) "b-offsets"
      do symub = 0,7
       lenub = 0
       do symb = 0,7
        symu = ieor(symub,symb)
        do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
         boff(b,symub) = lenub
         lenub = lenub + bf_per_ir(symu)
c        write(*,*) b,symub,boff(b,symub)
        end do
       end do
      end do
c
c     Addressing & sizes
c
c     write(*,*) "symv-offsets"
      maxlenuv = 0
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
c       write(*,*) symv,symub,symvoff(symv,symuv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
c     write(*,*) "maxlenuv = ",maxlenuv
c
c     Transformed MO coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf*nv(spina), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*nv(spina))
      call ga_get(g_part(spina),1,nbf,no(spina)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
c     call ma_print(dbl_mb(k_c),nbf,nv(spina),'MO coeff')
csh
c
c     Allocate g_uv_p and g_uv_m for t(u,b,i,j) and t[+/-](u,v,i,j)
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair', 
     $ maxlenuv,-1,g_uv_m)) call errquit
     $ ('uccsdt_tosopair: room for g_uv_m?',maxlenuv*lenij)
csh
c     call ga_fill(g_uv_m,9.9999d0)
csh
      if (spini .ne. spinj) then
       if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair', 
     $  maxlenuv,-1,g_uv_p)) call errquit
     $  ('uccsdt_tosopair: room for g_uv_p?',maxlenuv*lenij)
csh
c     call ga_fill(g_uv_p,9.9999d0)
csh
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
        lenb = bhi - blo + 1
        lenj = nv_sym(symj,spinj)
        if (lenb.gt.0 .and. lenj.gt.0) then
c        write(*,*) "geting a slice for symb, symjb = ",symb,symjb
         call sh_get_t2(g_t2,symb,symjb,spini,spinj,g_t2_big,
     $   "new")
c        call ga_print(g_t2)
csh      if (.not. uccsdt_ampfile_read_t2(d_t, 
csh  $    spini, spina, spinj, spinb, symia, 
csh  $    blo, bhi, g_t2, .true., 'column'))
csh  $    call errquit('read_t2 failed',0)
c
c     Data-parallel transformation of a
c
         call ga_distribution(g_t2,ga_nodeid(),jblo,jbhi,dummy,dummy)
c        write(*,*) 'jblo,jbhi = ',jblo,jbhi
         jb = 0
         do b = blo, bhi
          do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
c          write(*,*) 'b,j = ',b,j
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
             symub = ieor(symu,symb)
csh
c            write(*,*) 'symi,syma,symj,symb,symjb = ',
c    $       symi,syma,symj,symb,symjb
csh
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
csh
c             write(*,*) 'leni,lena,lenu = ',
c    $        leni,lena,lenu
csh
              if (.not. ma_push_get(mt_dbl, lenia, 'buf',
     $         l_buf, k_buf)) 
     $         call errquit('ma? lenia',lenia)
              if (.not. ma_push_get(mt_dbl, leniu, 'tmp',
     $         l_tmp, k_tmp)) 
     $         call errquit('ma? leniu',leniu)
c
c     (ia) -> (iu) for a (virtual) & u (SO) of same irrep 
c
c             write(*,*) 'read range = ',ia+1,'~',ia+lenia,
c    $                   'and',jb,'~',jb
              call ga_get(g_t2,ia+1,ia+lenia,jb,jb,
     $         dbl_mb(k_buf),leni)
csh
c             call ma_print(dbl_mb(k_buf),leni,lena,'ia pair')
csh
              call dgemm('n','t',lenu,leni,lena,
     $         1.0d0,dbl_mb(k_c),nbf,
     $         dbl_mb(k_buf),leni,
     $         0.0d0,dbl_mb(k_tmp),lenu)
csh
c             call ma_print(dbl_mb(k_tmp),lenu,leni,'ui pair')
csh
              ij = joff(j,symij)
              do i = ilo,ihi
               ij = ij + 1
c              write(*,*) 'i,j,b = ',i,j,b
c              write(*,*) 'g_uv_m range boff,ij = ',boff(b,symub),ij
               call ga_put(g_uv_m,
     $          boff(b,symub)+1,boff(b,symub)+lenu,
     $          ij,ij,dbl_mb(k_tmp+lenu*(i-ilo)),lenu)
              end do
csh
c             call ga_print(g_uv_m)
csh
              ia = ia + lenia
              if (.not. ma_pop_stack(l_tmp))
     $         call errquit('ma_pop?',l_tmp)
              if (.not. ma_pop_stack(l_buf))
     $         call errquit('ma_pop?',l_buf)
             end if
            end do
           end if
          end do
         end do
         if (.not. ga_destroy(g_t2)) call errquit('GA 999',0)
        endif
       end do
      end do
      write(*,*) '**************************'
      write(*,*) 'transformation a -> u done'
      write(*,*) '**************************'
c
c     Now have t(u,b,i,j) ... data parallel transform b in place 
c     (g_uv_m was allocated big enough to do this).
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
      if (.not. ma_push_get(mt_dbl, nbf*nv(spinb), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*nv(spinb))
      call ga_get(g_part(spinb),1,nbf,no(spinb)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
c     call ma_print(dbl_mb(k_c),nbf,nv(spinb),'MO coeff')
csh
      call ga_distribution(g_uv_m,ga_nodeid(),ijlo,ijhi,dummy,dummy)
csh
c     write(*,*) 'ijlo,ijhi = ',ijlo,ijhi
csh
      ij = 0
      do symij = 0,7
       symuv = symij
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
          ij = ij + 1
          if (ij.ge.ijlo .and. ij.le.ijhi) then
c          write(*,*) 'symij,symj,ij = ',symij,symj,ij
           do symb = 0,7
c           write(*,*) 'symb = ',symb
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
c            write(*,*) 'lenu,lenv,lenb = ',lenu,lenv,lenb
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
csh
c            call ma_print(dbl_mb(k_buf),lenu,lenb,'ub pair')
csh
             call dgemm('n','t',lenu,lenv,lenb,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        dbl_mb(k_c),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
csh
c            call ma_print(dbl_mb(k_tmp),lenu,lenv,'uv pair')
c            write(*,*) 'symvoff,ij = ',symvoff(symv,symuv),ij
csh
             call ga_put(g_uv_m,
     $        symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $        ij,ij,dbl_mb(k_tmp),lenu)
csh
c            call ga_print(g_uv_m)
csh
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
      write(*,*) '**************************'
      write(*,*) 'transformation b -> v done'
      write(*,*) '**************************'
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
c
c     Symmetrize if spini .ne. spinj
c
      if (spini .ne. spinj) then
csh
       write(*,*) 'spini .NEQ. spinj ... antisymmetrize'
csh
       call ga_copy(g_uv_m,g_uv_p)
       call ga_distribution(g_uv_m,ga_nodeid(),ijlo,ijhi,dummy,dummy)
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
      subroutine uccsdt_b_contract(rtdb,g_uv_p,
     $g_uv_m,g_omega,
csh   subroutine uccsdt_b_contract(g_uv_p,g_uv_m,g_omega,
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
      integer lenij,lenwl,lenuv
      integer lenu,lenv,lenl,lenw
      integer lengx,lengy
      integer symij,symj,symi,symwl,syml,symw
      integer symuv,symu,symv,symlv,symwu,symlu,symwv
      integer symjoff(0:7,0:7)
      integer symvoff(0:7,0:7)
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
      integer uoff1(nbf,0:7)
      integer voff1(nbf,0:7)
      integer woff1(nbf,0:7)
      integer loff1(nbf)
      integer uoff2(nbf,0:7)
      integer voff2(nbf,0:7)
      integer loff2(nbf,0:7)
      integer woff2(nbf)
      integer lenlv(0:7),lenwu(0:7),lenwv(0:7),lenlu(0:7)
      integer symwuoff(0:7)
      integer wlo,whi,llo,lhi
      integer u,v,w,l,uvwl
      integer nproc
      integer nxtask
      external nxtask
csh
      double precision dummy
      integer rtdb
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
csh   call sym_abelian_bas_info???? SO MANY ARGUMENTS!!!!!!
csh
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        symjoff(symj,symij) = lenij
        write(*,*) 'symj,symij,symjoff = ',
     $  symj,symij,symjoff(symj,symij)
        lenij = lenij + no_sym(symi,spini)*no_sym(symj,spinj)
       end do
      end do
      write(*,*) 'lenij = ',lenij
c
c     Addressing & sizes
c
      maxlenuv = 0
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        write(*,*) 'symv,symuv,symvoff = ',
     $  symv,symuv,symvoff(symv,symuv)
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
      write(*,*) 'maxlenuv = ',maxlenuv
c
c     Allocate g_omega_p and g_omega_m for Omega[+/-](u,v,i,j)
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'omega', 
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
csh
      do ush_blk = 1,nsh_block
      write(*,*) 'shell block No.',ush_blk,'  ',
     $sh_block(1,ush_blk),':',sh_block(2,ush_blk)
      enddo
csh
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
     $      lsh_blk_lo, lsh_blk_hi, loff1, lrange,
     $      vsh_blk_lo, vsh_blk_hi, voff1, vrange,
     $      wsh_blk_lo, wsh_blk_hi, woff1, wrange,
     $      ush_blk_lo, ush_blk_hi, uoff1, urange,
     $      lenlv, lenwu, symwuoff, lengx)
csh
          write(*,*) 'shell block info retrieved for',
     $    ush_blk,vsh_blk,wsh_blk,lsh_blk
          write(*,*) 'ush block range ',ush_blk_lo,':',ush_blk_hi
          write(*,*) 'vsh block range ',vsh_blk_lo,':',vsh_blk_hi
          write(*,*) 'wsh block range ',wsh_blk_lo,':',wsh_blk_hi
          write(*,*) 'lsh block range ',lsh_blk_lo,':',lsh_blk_hi
          write(*,*) 'ush range ',urange
          write(*,*) 'vsh range ',vrange
          write(*,*) 'wsh range ',wrange
          write(*,*) 'lsh range ',lrange
          write(*,*) 'u offset ',uoff1
          write(*,*) 'v offset ',voff1
          write(*,*) 'w offset ',woff1
          write(*,*) 'l offset ',loff1
          write(*,*) 'lenlv ',lenlv
          write(*,*) 'lenwu ',lenwu
          write(*,*) 'lengx ',lengx
csh
c     
           if (.not. ma_push_get
     $      (mt_dbl, lengx, 'b_contract', l_gx, k_gx))
     $      call errquit('b_contract: no room for g', lengx)
c     
           call soints_block(basis, nbf,
     $      lsh_blk_lo, lsh_blk_hi, loff1, lrange,
     $      vsh_blk_lo, vsh_blk_hi, voff1, vrange,
     $      wsh_blk_lo, wsh_blk_hi, woff1, wrange,
     $      ush_blk_lo, ush_blk_hi, uoff1, urange,
     $      1d-12, dbl_mb(k_gx))
          write(*,*) 'shell block soints computed'
c     
c     Compute g(l,v,w,u) = (wv|lu) = <wl|vu> 
c     
           call soints_block_info(basis, nbf,
     $      wsh_blk_lo, wsh_blk_hi, woff2, wrange,
     $      vsh_blk_lo, vsh_blk_hi, voff2, vrange,
     $      lsh_blk_lo, lsh_blk_hi, loff2, lrange,
     $      ush_blk_lo, ush_blk_hi, uoff2, urange,
     $      lenwv, lenlu, symwuoff, lengy)
csh
          write(*,*) 'shell block info retrieved for',
     $    ush_blk,vsh_blk,wsh_blk,lsh_blk
          write(*,*) 'ush block range ',ush_blk_lo,':',ush_blk_hi
          write(*,*) 'vsh block range ',vsh_blk_lo,':',vsh_blk_hi
          write(*,*) 'lsh block range ',lsh_blk_lo,':',lsh_blk_hi
          write(*,*) 'wsh block range ',wsh_blk_lo,':',wsh_blk_hi
          write(*,*) 'ush range ',urange
          write(*,*) 'vsh range ',vrange
          write(*,*) 'lsh range ',lrange
          write(*,*) 'wsh range ',wrange
          write(*,*) 'u offset ',uoff2
          write(*,*) 'v offset ',voff2
          write(*,*) 'l offset ',loff2
          write(*,*) 'w offset ',woff2
          write(*,*) 'lenlu ',lenlu
          write(*,*) 'lenwv ',lenwv
          write(*,*) 'lengy ',lengy
csh
           if (.not. ma_push_get
     $      (mt_dbl, lengy, 'b_contract', l_gy, k_gy))
     $      call errquit('b_contract: no room for g', lengy)
c     
           call soints_block(basis, nbf,
     $      wsh_blk_lo, wsh_blk_hi, woff2, wrange,
     $      vsh_blk_lo, vsh_blk_hi, voff2, vrange,
     $      lsh_blk_lo, lsh_blk_hi, loff2, lrange,
     $      ush_blk_lo, ush_blk_hi, uoff2, urange,
     $      1d-12, dbl_mb(k_gy))
          write(*,*) 'shell block soints computed'
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
csh
                 write(*,*) 'symuv,symi,symj,symu,symv,symw,syml = ',
     $           symuv,symi,symj,symu,symv,symw,syml
csh
                 symwu = ieor(symw,symu)
                 symlv = ieor(syml,symv)
                 symlu = ieor(syml,symu)
                 symwv = ieor(symw,symv)
                 lenl = lrange(2,syml) - lrange(1,syml) + 1
                 lenw = wrange(2,symw) - wrange(1,symw) + 1
                 lenwl = lenw * lenl
                 wlo = wrange(1,symw) - bf_per_ir_cum(symw)
                 whi = wrange(2,symw) - bf_per_ir_cum(symw)
                 llo = lrange(1,syml) - bf_per_ir_cum(syml)
                 lhi = lrange(2,syml) - bf_per_ir_cum(syml)
csh
                 write(*,*) 'lenij,lenuv,lenwl = ',
     $           lenij,lenuv,lenwl
csh
                 if (lenwl .gt. 0) then
                  if (.not. ma_push_get(mt_dbl, 
     $             lenwl*lenij, 't', l_t, k_t)) 
     $             call errquit('ma? lenwl',lenwl*lenij)
c
c     g_uv_m(w,l,syml,i,j,symj,symij) -> t(w,l,i,j)
c
csh
                  write(*,*) 'getting t2 in the range ...',
     $            symvoff(syml,symwl)+1,'~',symvoff(syml,symwl)+lenwl,
     $            symjoff(symj,symij)+1,'~',symjoff(symj,symij)+lenij
csh
                  call ga_get(g_uv_m,
     $            symvoff(syml,symwl)+1,symvoff(syml,symwl)+lenwl,
     $            symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $            dbl_mb(k_t),lenwl)
csh
                  call ma_print(dbl_mb(k_t),lenwl,lenij,'SO t-amp')
csh
                  if (.not. ma_push_get(mt_dbl, 
     $             lenuv*lenwl, 'minus', l_m, k_m)) 
     $             call errquit('ma? lenuv*lenwl',lenuv*lenwl)
                  uvwl = 0
                  do l = lrange(1,syml), lrange(2,syml)
                   do w = wrange(1,symw), wrange(2,symw)
                    do v = vrange(1,symv), vrange(2,symv)
                     do u = urange(1,symu), urange(2,symu)
                      uvwl = uvwl + 1
                      dbl_mb(k_m+uvwl-1) =
     $                dbl_mb(k_gx+loff1(l)+voff1(v,symlv)
     $                 +woff1(w,symwu)+uoff1(u,symwu)-1)
     $                -dbl_mb(k_gy+woff2(w)+voff2(v,symwv)
     $                 +loff2(l,symlu)+uoff2(u,symlu)-1)
                     end do
                    end do
                   end do
                  end do
csh
                  write(*,*) 'lenuv,lenwl = ',lenuv, lenwl
                  call ma_print(dbl_mb(k_m),lenuv,lenwl,'SO int')
csh
                  call dgemm('n','n',lenuv,lenij,lenwl,
     $             1.0d0,dbl_mb(k_m),lenuv,
     $             dbl_mb(k_t),lenwl,
     $             0.0d0,dbl_mb(k_omega),lenuv)
csh
                  call ma_print(dbl_mb(k_omega),lenuv,lenij,'Omega')
csh
                  if (symw .eq. syml) then
csh
                   write(*,*) 'accumulating in ...',
     $             symvoff(symv,symuv)+1,'~',symvoff(symv,symuv)+lenuv,
     $             symjoff(symj,symij)+1,'~',symjoff(symj,symij)+lenij
csh
                   call ga_acc(g_omega,
     $             symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $             symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $             dbl_mb(k_omega),lenuv,0.5d0)
                  else
                   call ga_acc(g_omega,
     $             symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $             symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $             dbl_mb(k_omega),lenuv,1.0d0)
                   call ga_acc(g_omega,
     $             symvoff(symu,symuv)+1,symvoff(symu,symuv)+lenuv,
     $             symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $             dbl_mb(k_omega),lenuv,1.0d0)
                  end if
                  if (.not. ma_pop_stack(l_m))
     $             call errquit('ma_pop? l_m',l_m)
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
                       dbl_mb(k_p+uvwl-1) = 
     $                 dbl_mb(k_gx+loff1(l)+voff1(v,symlv)
     $                  +woff1(w,symwu)+uoff1(u,symwu)-1)
     $                 +dbl_mb(k_gy+woff2(w)+voff2(v,symwv)
     $                  +loff2(l,symlu)+uoff2(u,symlu)-1)
                      end do
                     end do
                    end do
                   end do
                   call dgemm('n','n',lenuv,lenij,lenwl,
     $              1.0d0,dbl_mb(k_p),lenuv,
     $              dbl_mb(k_t),lenwl,
     $              0.0d0,dbl_mb(k_omega),lenuv)
                   if (symu .eq. symv) then
                    call ga_acc(g_omega,
     $              symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $              symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $              dbl_mb(k_omega),lenuv,0.5d0)
                   else
                    call ga_acc(g_omega,
     $              symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $              symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $              dbl_mb(k_omega),lenuv,1.0d0)
                    call ga_acc(g_omega,
     $              symvoff(symu,symuv)+1,symvoff(symu,symuv)+lenuv,
     $              symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $              dbl_mb(k_omega),lenuv,1.0d0)
                   end if
                   if (.not. ma_pop_stack(l_p))
     $              call errquit('ma_pop? l_p',l_p)
                  end if
                  if (.not. ma_pop_stack(l_t))
     $             call errquit('ma_pop? l_t',l_t)
                 end if
                end do
                if (.not. ma_pop_stack(l_omega))
     $           call errquit('ma_pop? l_omega',l_omega)
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
      next = nxtask(-nproc,task_size)
      call ga_sync()
      if (spini .ne. spinj) call ga_scale(g_omega,0.5d0)
c
csh
      call int_terminate()
      call schwarz_tidy()
csh
      return
      end
      subroutine uccsdt_makeb(g_omega,g_b_big,spini,spinj)
csh   subroutine uccsdt_makeb(g_omega,d_b,spini,spinj)
c
c     b(a,b,i,j) = Sum(u,v) t(u,v,i,j)C(u,a,spini)C(v,b,spinj)
c     where u & v are SO's (no spin), a,b,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
csh   integer d_b          ! [output] DRA handle for b(a,b,i,j)
      integer g_b_big      ! [output] DRA handle for b(a,b,i,j)
      integer g_omega      ! [input] GA handle for Omega(u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer spina, spinb ! Spins of a and b
      integer g_tmp        ! Temporary GA t(u,b,i,j) & t(u,v,i,j)
      integer g_b          ! Temporary GA b(a,b,i,j)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer k_c,l_c      ! Temporary MA for MO coeffs
      integer lenij,lenub,lenuv,lenu,leniu
      integer lenv,lenb,lenia,lenjb,leni,lena
      integer maxlenuv
      integer symij,symj,symi,symub,symb
      integer symuv,symu,symv,symjb,symia,syma
      integer joff(nw_max_nbf_small,0:7)
      integer boff(nw_max_nbf_small,0:7)
      integer symvoff(0:7,0:7)
      integer symboff(0:7,0:7)
      integer ij,jb,ia
      integer ijlo,ijhi,jblo,jbhi
      integer i,j,a,b
      integer blo,bhi,vlo,alo,ilo,ihi,ulo,uhi
      integer dummy
c
      spina = spini
      spinb = spinj
csh ***** We have to create an empty g_b_big
      dummy = 0
      do symjb = 0,7
      do symb = 0,7
      symj = ieor(symjb,symj)
      do b = 1,nv_sym(symb,spinb)
      do j = 1,no_sym(symj,spinj)
      do syma = 0,7
      symi = ieor(symjb,syma)
      do a = 1,nv_sym(syma,spina)
      do i = 1,no_sym(symi,spini)
      dummy = dummy + 1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      write(*,*) 'Size of g_b_big = ',dummy
      if (.not.ga_create(mt_dbl,dummy,1,'b_big',
     $ -1,1,g_b_big)) call errquit
     $ ('uccsdt_makeb: room for b_big?',dummy)
csh
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         joff(j,symij) = lenij
         write(*,*) 'symij,j,joff = ',
     $   symij,j,joff(j,symij)
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
         write(*,*) 'symub,b,boff = ',
     $   symub,b,boff(b,symub)
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
        write(*,*) 'symub,symb,symboff = ',
     $  symub,symb,symboff(symb,symub)
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
        write(*,*) 'symuv,symv,symvoff = ',
     $  symuv,symv,symvoff(symv,symuv)
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
csh
      write(*,*) 'maxlenuv = ',maxlenuv
csh
c
c     Transformed MO coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf*nv(spinb), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*nv',nbf*nv(spinb))
      call ga_get(g_part(spinb),1,nbf,no(spinb)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
      call ma_print(dbl_mb(k_c),nbf,nv(spinb),'MO coeff')
csh
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
csh
            write(*,*) 'symij,symj,j,i,symb = ',
     $                  symij,symj,j,i,symb
csh
            if (lenuv.gt.0 .and. lenb.gt.0) then
csh
             write(*,*) 'lenub,lenuv = ',lenub,lenuv
csh
             if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $        l_buf, k_buf)) 
     $        call errquit('ma? lenuv',lenuv)
             if (.not. ma_push_get(mt_dbl, lenub, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenub',lenub)
c
c     (uv) -> (ub) for b (virtual) & v (SO) of same irrep 
c
             call ga_get(g_omega,1,lenuv,ij,ij,
     $        dbl_mb(k_buf),lenu)
csh
             call ma_print(dbl_mb(k_buf),lenu,lenv,'UV')
csh
             call dgemm('n','n',lenu,lenb,lenv,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        dbl_mb(k_c),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
csh
             call ma_print(dbl_mb(k_tmp),lenu,lenb,'UB')
             write(*,*) 'ga_put in the range ... ',
     $       symboff(symb,symub)+1,'~',symboff(symb,symub)+lenub,
     $       ij,'~',ij
csh
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
csh
      write(*,*) '**************************'
      write(*,*) 'transformation v -> b done'
      write(*,*) '**************************'
csh
c
c     Perform quarter transformation from tmp(u,b,i,j) to t(a,b,i,j)
c     t(i,a,syma,j,b,symb,symjb)
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
      if (.not. ma_push_get(mt_dbl, nbf*nv(spina), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*nv',nbf*nv(spina))
      call ga_get(g_part(spina),1,nbf,no(spina)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
      call ma_print(dbl_mb(k_c),nbf,nv(spina),'MO coeff')
csh
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
        if ((lenia.gt.0).and.(lenjb.gt.0)) then
csh
         write(*,*) 'size of g_b: lenia, lenjb = ',lenia,lenjb
csh
         if (.not.ga_create(mt_dbl,lenia,lenjb,'b', 
     $    lenia,-1,g_b)) call errquit
     $    ('uccsdt_makeb: room for b?',lenia*lenjb)
c
c     Data-parallel transformation of a
c
         call ga_distribution(g_b,ga_nodeid(),jblo,jbhi,dummy,dummy)
         jb = 0
         do b = blo, bhi
          do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
           jb = jb + 1
           if (jb.ge.jblo .and. jb.le.jbhi) then
            ia = 0
c     Compare the next line with the comments in Amplitude.F
c     I think the comments are not accurate
            do syma = 0, 7
             symi = ieor(symia,syma)
             symij = ieor(symi,symj)
             symu = syma
             symub = ieor(symu,symb)
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
csh
             write(*,*) 'symjb,symb,b,j,syma = ',
     $                   symjb,symb,b,j,syma
csh
             if (lenu.gt.0) then
              if (.not. ma_push_get(mt_dbl, leniu, 'buf',
     $         l_buf, k_buf)) 
     $         call errquit('ma? leniu',leniu)
              if (.not. ma_push_get(mt_dbl, lenia, 'tmp',
     $         l_tmp, k_tmp)) 
     $         call errquit('ma? lenia',lenia)
c
c     (iu) -> (ia) for a (virtual) & u (SO) of same irrep 
c
              ij = joff(j,symij)
              do i = ilo,ihi
               ij = ij + 1
csh
               write(*,*) 'i = ',i
               write(*,*) 'get from range ... ',
     $          boff(b,symub)+1,'~',boff(b,symub)+lenu,
     $          ij,'~',ij,
     $          ' to k_buf+',lenu*(i-ilo)
csh
               call ga_get(g_tmp,
     $          boff(b,symub)+1,boff(b,symub)+lenu,
     $          ij,ij,dbl_mb(k_buf+lenu*(i-ilo)),lenu)
              end do
csh
              call ma_print(dbl_mb(k_buf),lenu,leni,'B(UI)')
csh
              call dgemm('t','n',leni,lena,lenu,
     $         1.0d0,dbl_mb(k_buf),lenu,
     $         dbl_mb(k_c),nbf,
     $         0.0d0,dbl_mb(k_tmp),leni)
csh
              call ma_print(dbl_mb(k_tmp),leni,lena,'B(IA)')
              dummy = 0
              do a = 1,nv_sym(syma,spina)
              do i = 1,no_sym(symi,spini)
                dummy = dummy + 1
                write(*,8) a,b-no(spinb),i,j,dbl_mb(k_tmp+dummy-1)
    8           format(4i3,f15.10)
              enddo
              enddo
              write(*,*) 'writing to range ... ',
     $        ia+1,'~',ia+lenia,jb,'~',jb
csh
              call ga_put(g_b,ia+1,ia+lenia,jb,jb,
     $         dbl_mb(k_tmp),leni)
              ia = ia + lenia
              if (.not. ma_pop_stack(l_tmp))
     $         call errquit('ma_pop?',l_tmp)
              if (.not. ma_pop_stack(l_buf))
     $         call errquit('ma_pop?',l_buf)
             end if
            end do
           end if
          end do
         end do
c
c     Write one symmetry block of B to disk
c
         call ga_print(g_b)
         call sh_put_t2(g_b,symb,symjb,spini,spinj,g_b_big)
csh      if (.not. uccsdt_ampfile_write_t2(d_b, 
csh  $    spini, spina, spinj, spinb, symia, 
csh  $    blo, bhi, g_b))
csh  $    call errquit('write_t2 failed',0)
         if (.not. ga_destroy(g_b)) call errquit('GA 999',0)
        endif
       end do
      end do
csh
      write(*,*) '**************************'
      write(*,*) 'transformation u -> a done'
      write(*,*) '**************************'
csh
      if (.not. ga_destroy(g_tmp)) call errquit('GA 999',0)
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
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
      integer k_c,l_c      ! Temporary MA for MO coeffs
      integer lenij,lenul,lenuv,lenu,lenv,lenl
      integer lenik,lenjl,leni,lenk,leniu
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
      integer dummy
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
c     Transformed MO coefficients
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
      if (.not. ma_push_get(mt_dbl, nbf*no(spinl), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spinl))
      call ga_get(g_hole(spinl),1,nbf,1,no(spinl),
     $ dbl_mb(k_c),nbf)
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
     $        dbl_mb(k_c),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
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
      if (.not. ma_push_get(mt_dbl, nbf*no(spink), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spink))
      call ga_get(g_hole(spink),1,nbf,1,no(spink),
     $ dbl_mb(k_c),nbf)
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
        call ga_distribution(g_tmp,ga_nodeid(),jllo,jlhi,dummy,dummy)
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
     $        1.0d0,dbl_mb(k_c),nbf,
     $        dbl_mb(k_tmp),leni,
     $        0.0d0,dbl_mb(k_buf),lenu)
             call ga_put(g_x,ik+1,ik+lenik,jl,jl,
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
        if (.not. uccsdt_ampfile_write_t2(d_x, 
     $   spini, spink, spinj, spinl, symik, 
     $   llo, lhi, g_x))
     $   call errquit('write_t2 failed',0)
        if (.not. ga_destroy(g_x)) call errquit('GA 999',0)
       end do
      end do
      if (.not. ga_destroy(g_tmp)) call errquit('GA 999',0)
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
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
