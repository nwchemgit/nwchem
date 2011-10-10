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
#include "errquit.fh"
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
      integer joff(nbf,0:7)
      integer boff(nbf,0:7)
      integer symvoff(0:7,0:7)
      integer blo,bhi,alo,ilo,ihi,ulo,uhi,vlo
      integer jblo,jbhi,ijlo,ijhi
      integer dummy
      double precision x,y
csh
      logical odebug
c     logical int_chk_init
c     external int_chk_init
c     logical temp
c     temp = int_chk_init('b term sopair')
      odebug = .false.
csh 
c
      spina = spini
      spinb = spinj
csh
      if (odebug) then
      write(*,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      write(*,*) 'X The addressing of k_c must be wrong for sym=on X'
      write(*,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      write(*,*) 'no = ',no(1),no(2)
      write(*,*) 'nv = ',nv(1),nv(2)
      write(*,*) 'no_sym = ',(no_sym(i,1),i=0,7),(no_sym(i,2),i=0,7)
      write(*,*) 'nv_sym = ',(nv_sym(i,1),i=0,7),(nv_sym(i,2),i=0,7)
      write(*,*) 'o_sym(a) = ',(o_sym(1,i,1),i=0,7),(o_sym(2,i,1),i=0,7)
      write(*,*) 'o_sym(b) = ',(o_sym(1,i,2),i=0,7),(o_sym(2,i,2),i=0,7)
      write(*,*) 'v_sym(a) = ',(v_sym(1,i,1),i=0,7),(v_sym(2,i,1),i=0,7)
      write(*,*) 'v_sym(b) = ',(v_sym(1,i,2),i=0,7),(v_sym(2,i,2),i=0,7)
      write(*,*) 'bf_per_ir = ',(bf_per_ir(i),i=0,7)
      endif
csh
c
c     Addressing & sizes
c
      if (odebug) write(*,*) "j-offsets"
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         joff(j,symij) = lenij
         lenij = lenij + no_sym(symi,spini)
         if (odebug) write(*,*) j,symij,joff(j,symij)
        end do
       end do
      end do
      if (odebug) write(*,*) "lenij = ",lenij
c
c     Addressing & sizes
c
      if (odebug) write(*,*) "b-offsets"
      do symub = 0,7
       lenub = 0
       do symb = 0,7
        symu = ieor(symub,symb)
        do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
         boff(b,symub) = lenub
         lenub = lenub + bf_per_ir(symu)
         if (odebug) write(*,*) b,symub,boff(b,symub)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      if (odebug) write(*,*) "symv-offsets"
      maxlenuv = 0
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
        if (odebug) write(*,*) symv,symub,symvoff(symv,symuv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
      if (odebug) write(*,*) "maxlenuv = ",maxlenuv
c
c     Transformed MO coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf*nv(spina), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*nv(spina), MA_ERR)
      call ga_get(g_hole(spina),1,nbf,no(spina)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
      if (odebug)
     $call ma_print(dbl_mb(k_c),nbf,nv(spina),'MO coeff')
csh
c
c     Allocate g_uv_p and g_uv_m for t(u,b,i,j) and t[+/-](u,v,i,j)
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair -', 
     $ maxlenuv,-1,g_uv_m)) call errquit
     $ ('uccsdt_tosopair: room for g_uv_m?',maxlenuv*lenij, GA_ERR)
csh
      if (odebug)
     $call ga_fill(g_uv_m,9.9999d0)
csh
      if (spini .ne. spinj) then
       if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair +', 
     $  maxlenuv,-1,g_uv_p)) call errquit
     $  ('uccsdt_tosopair: room for g_uv_p?',maxlenuv*lenij, GA_ERR)
csh
      if (odebug)
     $call ga_fill(g_uv_p,9.9999d0)
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
         if (odebug)
     $   write(*,*) "geting a slice for symb, symjb = ",symb,symjb
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
         if (odebug) write(*,*) 'jblo,jbhi = ',jblo,jbhi
         jb = 0
         do b = blo, bhi
          do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
           if (odebug) write(*,*) 'b,j = ',b,j
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
             if (odebug)
     $       write(*,*) 'symi,syma,symj,symb,symjb = ',
     $       symi,syma,symj,symb,symjb
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
              if (odebug)
     $        write(*,*) 'leni,lena,lenu = ',
     $        leni,lena,lenu
csh
              if (.not. ma_push_get(mt_dbl, lenia, 'buf',
     $         l_buf, k_buf)) 
     $         call errquit('ma? lenia',lenia, MA_ERR)
              if (.not. ma_push_get(mt_dbl, leniu, 'tmp',
     $         l_tmp, k_tmp)) 
     $         call errquit('ma? leniu',leniu, MA_ERR)
c
c     (ia) -> (iu) for a (virtual) & u (SO) of same irrep 
c
c             write(*,*) 'read range = ',ia+1,'~',ia+lenia,
c    $                   'and',jb,'~',jb
              call ga_get(g_t2,ia+1,ia+lenia,jb,jb,
     $         dbl_mb(k_buf),leni)
csh
              if (odebug)
     $        call ma_print(dbl_mb(k_buf),leni,lena,'ia pair')
csh
              call dgemm('n','t',lenu,leni,lena,
     $         1.0d0,dbl_mb(k_c),nbf,
     $         dbl_mb(k_buf),leni,
     $         0.0d0,dbl_mb(k_tmp),lenu)
csh
              if (odebug)
     $        call ma_print(dbl_mb(k_tmp),lenu,leni,'ui pair')
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
              if (odebug)
     $        call ga_print(g_uv_m)
csh
              ia = ia + lenia
              if (.not. ma_pop_stack(l_tmp))
     $         call errquit('ma_pop?',l_tmp, MA_ERR)
              if (.not. ma_pop_stack(l_buf))
     $         call errquit('ma_pop?',l_buf, MA_ERR)
             end if
            end do
           end if
          end do
         end do
         if (.not. ga_destroy(g_t2)) call errquit('GA 999',0, GA_ERR)
        endif
       end do
      end do
      if (odebug) then
      write(*,*) '**************************'
      write(*,*) 'transformation a -> u done'
      write(*,*) '**************************'
      endif
c
c     Now have t(u,b,i,j) ... data parallel transform b in place 
c     (g_uv_m was allocated big enough to do this).
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nv(spinb), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*nv(spinb), MA_ERR)
      call ga_get(g_hole(spinb),1,nbf,no(spinb)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
      if (odebug)
     $call ma_print(dbl_mb(k_c),nbf,nv(spinb),'MO coeff')
csh
      call ga_distribution(g_uv_m,ga_nodeid(),ijlo,ijhi,dummy,dummy)
csh
      if (odebug)
     $write(*,*) 'ijlo,ijhi = ',ijlo,ijhi
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
           if (odebug)
     $     write(*,*) 'symij,symj,ij = ',symij,symj,ij
           do symb = 0,7
           if (odebug)
     $      write(*,*) 'symb = ',symb
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
             if (odebug)
     $       write(*,*) 'lenu,lenv,lenb = ',lenu,lenv,lenb
             if (.not. ma_push_get(mt_dbl, lenub, 'buf',
     $        l_buf, k_buf)) 
     $       call errquit('ma? lenub',lenub, MA_ERR)
             if (.not. ma_push_get(mt_dbl, lenuv, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenuv',lenuv, MA_ERR)
c
c     (ub) -> (uv) for b (virtual) & v (SO) of same irrep 
c
             call ga_get(g_uv_m,1,lenub,ij,ij,
     $        dbl_mb(k_buf),lenu)
csh
             if (odebug)
     $       call ma_print(dbl_mb(k_buf),lenu,lenb,'ub pair')
csh
             call dgemm('n','t',lenu,lenv,lenb,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        dbl_mb(k_c),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
csh
             if (odebug) then
             call ma_print(dbl_mb(k_tmp),lenu,lenv,'uv pair')
             write(*,*) 'symvoff,ij = ',symvoff(symv,symuv),ij
             endif
csh
             call ga_put(g_uv_m,
     $        symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $        ij,ij,dbl_mb(k_tmp),lenu)
csh
             if (odebug)
     $       call ga_print(g_uv_m)
csh
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp, MA_ERR)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf, MA_ERR)
            end if
           end do
          end if
         end do
        end do
       end do
      end do
      if (odebug) then
      write(*,*) '**************************'
      write(*,*) 'transformation b -> v done'
      write(*,*) '**************************'
      endif
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c, MA_ERR)
c
c     Symmetrize if spini .ne. spinj
c
      if (spini .ne. spinj) then
csh
       if (odebug)
     $ write(*,*) 'spini .NEQ. spinj ... antisymmetrize'
       if (odebug) call ga_print(g_uv_m)
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
              if (odebug) write(*,*) 'symij,symj,j,i,symv,symu = ',
     $        symij,symj,j,i,symv,symu
              lenu = bf_per_ir(symu)
              lenv = bf_per_ir(symv)
              lenuv = lenu * lenv
              if (lenuv .gt. 0) then
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'buf', l_buf, k_buf)) 
     $          call errquit('ma? lenuv',lenuv, MA_ERR)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'tmp', l_tmp, k_tmp)) 
     $          call errquit('ma? lenuv',lenuv, MA_ERR)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'plus', l_p1, k_p1)) 
     $          call errquit('ma? lenuv',lenuv, MA_ERR)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'minus', l_m1, k_m1)) 
     $          call errquit('ma? lenuv',lenuv, MA_ERR)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'plus', l_p2, k_p2)) 
     $          call errquit('ma? lenuv',lenuv, MA_ERR)
               if (.not. ma_push_get(mt_dbl, 
     $          lenuv, 'minus', l_m2, k_m2)) 
     $          call errquit('ma? lenuv',lenuv, MA_ERR)
               if (odebug) write(*,*) 'getting range ... ',
     $         symvoff(symu,symuv)+1,'~',symvoff(symu,symuv)+lenuv,
     $          ij,'~',ij
               call ga_get(g_uv_m,
     $          symvoff(symu,symuv)+1,symvoff(symu,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_tmp),lenv)
               if (odebug) call ma_print(dbl_mb(k_tmp),lenv,lenu,'VU')
               if (odebug) write(*,*) 'getting range ... ',
     $         symvoff(symv,symuv)+1,'~',symvoff(symv,symuv)+lenuv,
     $          ij,'~',ij
               call ga_get(g_uv_m,
     $          symvoff(symv,symuv)+1,symvoff(symv,symuv)+lenuv,
     $          ij,ij,dbl_mb(k_buf),lenu)
               if (odebug) call ma_print(dbl_mb(k_buf),lenu,lenv,'UV')
               do v = 1,lenv
                do u = 1,lenu
                 x = dbl_mb(k_tmp+(u-1)*lenv+v-1)
                 y = dbl_mb(k_buf+(v-1)*lenu+u-1)
                 if (odebug) write(*,*) 'U,V, and X,Y = ',u,v,x,y
                 dbl_mb(k_p1+(v-1)*lenu+u-1) = 0.5d0*(x+y)
                 dbl_mb(k_m1+(v-1)*lenu+u-1) = 0.5d0*(x-y)
                 dbl_mb(k_p2+(u-1)*lenv+v-1) = 0.5d0*(x+y)
                 dbl_mb(k_m2+(u-1)*lenv+v-1) = 0.5d0*(x-y)
                end do
               end do
               if (odebug) call ma_print(dbl_mb(k_p1),lenu,lenv,'UV+')
               if (odebug) call ma_print(dbl_mb(k_m1),lenu,lenv,'UV-')
               if (odebug) call ma_print(dbl_mb(k_p2),lenv,lenu,'VU+')
               if (odebug) call ma_print(dbl_mb(k_m2),lenv,lenu,'VU-')
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
     $          call errquit('ma_pop?',l_m2, MA_ERR)
               if (.not. ma_pop_stack(l_p2))
     $          call errquit('ma_pop?',l_p2, MA_ERR)
               if (.not. ma_pop_stack(l_m1))
     $          call errquit('ma_pop?',l_m1, MA_ERR)
               if (.not. ma_pop_stack(l_p1))
     $          call errquit('ma_pop?',l_p1, MA_ERR)
               if (.not. ma_pop_stack(l_tmp))
     $          call errquit('ma_pop?',l_tmp, MA_ERR)
               if (.not. ma_pop_stack(l_buf))
     $          call errquit('ma_pop?',l_buf, MA_ERR)
              end if
             end if
            end do
           end if
          end do
         end do
        end do
       end do
       if (odebug) call ga_print(g_uv_p)
       if (odebug) call ga_print(g_uv_m)
      end if
c
      return
      end
      subroutine uccsdt_b_contract(rtdb,g_uv_p,
     $g_uv_m,g_omega,
csh   subroutine uccsdt_b_contract(g_uv_p,g_uv_m,g_omega,
     $ spini,spinj)
      implicit none
#include "errquit.fh"
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
      integer rtdb
      logical odebug
      odebug = .false.
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
csh   call sym_abelian_bas_info???? SO MANY ARGUMENTS!!!!!!
      if (odebug) write(*,*) '******************************'
      if (odebug) write(*,*) '* What about wl <= uv ?????? *'
      if (odebug) write(*,*) '******************************'
csh
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        symjoff(symj,symij) = lenij
        if (odebug) write(*,*) 'symj,symij,symjoff = ',
     $  symj,symij,symjoff(symj,symij)
        lenij = lenij + no_sym(symi,spini)*no_sym(symj,spinj)
       end do
      end do
      if (odebug) write(*,*) 'lenij = ',lenij
c
c     Addressing & sizes
c
      maxlenuv = 0
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        if (odebug) write(*,*) 'symv,symuv,symvoff = ',
     $  symv,symuv,symvoff(symv,symuv)
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
      if (odebug) write(*,*) 'maxlenuv = ',maxlenuv
c
c     Allocate g_omega_p and g_omega_m for Omega[+/-](u,v,i,j)
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'omega', 
     $ maxlenuv,-1,g_omega)) call errquit
     $ ('uccsdt_tosopair: room for tmp?',maxlenuv*lenij, GA_ERR)
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
      if (odebug) write(*,*) 'shell block No.',ush_blk,'  ',
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
          if (odebug) then
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
          endif
csh
c     
           if (.not. ma_push_get
     $      (mt_dbl, lengx, 'b_contract', l_gx, k_gx))
     $      call errquit('b_contract: no room for g', lengx, MA_ERR)
c     
           call soints_block(basis, nbf,
     $      lsh_blk_lo, lsh_blk_hi, loff1, lrange,
     $      vsh_blk_lo, vsh_blk_hi, voff1, vrange,
     $      wsh_blk_lo, wsh_blk_hi, woff1, wrange,
     $      ush_blk_lo, ush_blk_hi, uoff1, urange,
     $      1d-12, dbl_mb(k_gx))
          if (odebug) write(*,*) 'shell block soints computed'
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
          if (odebug) then
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
          endif
csh
           if (.not. ma_push_get
     $      (mt_dbl, lengy, 'b_contract', l_gy, k_gy))
     $      call errquit('b_contract: no room for g', lengy, MA_ERR)
c     
           call soints_block(basis, nbf,
     $      wsh_blk_lo, wsh_blk_hi, woff2, wrange,
     $      vsh_blk_lo, vsh_blk_hi, voff2, vrange,
     $      lsh_blk_lo, lsh_blk_hi, loff2, lrange,
     $      ush_blk_lo, ush_blk_hi, uoff2, urange,
     $      1d-12, dbl_mb(k_gy))
          if (odebug) write(*,*) 'shell block soints computed'
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
     $          call errquit('ma? lenuv*lenij',lenuv*lenij, MA_ERR)
                do syml = 0,7
                 symw = ieor(symwl,syml)
csh
                 if (odebug)
     $           write(*,*) 'symuv,symi,symj,symu,symv,symw,syml = ',
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
                 if (odebug)
     $           write(*,*) 'lenij,lenuv,lenwl = ',
     $           lenij,lenuv,lenwl
csh
                 if (lenwl .gt. 0) then
                  if (.not. ma_push_get(mt_dbl, 
     $             lenwl*lenij, 't', l_t, k_t)) 
     $             call errquit('ma? lenwl',lenwl*lenij, MA_ERR)
c
c     g_uv_m(w,l,syml,i,j,symj,symij) -> t(w,l,i,j)
c
csh
                  if (odebug) write(*,*) 'getting t2 in the range ...',
     $            symvoff(syml,symwl)+1,'~',symvoff(syml,symwl)+lenwl,
     $            symjoff(symj,symij)+1,'~',symjoff(symj,symij)+lenij
csh
                  call ga_get(g_uv_m,
     $            symvoff(syml,symwl)+1,symvoff(syml,symwl)+lenwl,
     $            symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $            dbl_mb(k_t),lenwl)
csh
         if (odebug) call ma_print(dbl_mb(k_t),lenwl,lenij,'SO t-amp -')
csh
                  if (.not. ma_push_get(mt_dbl, 
     $             lenuv*lenwl, 'minus', l_m, k_m)) 
     $             call errquit('ma? lenuv*lenwl',lenuv*lenwl, MA_ERR)
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
       if (odebug) write(*,*) 'lenuv,lenwl = ',lenuv, lenwl
       if (odebug) call ma_print(dbl_mb(k_m),lenuv,lenwl,'SO int')
       if (odebug) call ma_print(dbl_mb(k_t),lenwl,lenij,'t (-)')
csh
                  call dgemm('n','n',lenuv,lenij,lenwl,
     $             1.0d0,dbl_mb(k_m),lenuv,
     $             dbl_mb(k_t),lenwl,
     $             0.0d0,dbl_mb(k_omega),lenuv)
      if (odebug) call ma_print(dbl_mb(k_omega),lenuv,lenij,'UVIJ -')
                  if (.not. ma_pop_stack(l_m))
     $             call errquit('ma_pop? l_m',l_m, MA_ERR)
c
c     Mixed spin
c
                  if (spini .ne. spinj) then
                   call ga_get(g_uv_p,
     $             symvoff(syml,symwl)+1,symvoff(syml,symwl)+lenwl,
     $             symjoff(symj,symij)+1,symjoff(symj,symij)+lenij,
     $             dbl_mb(k_t),lenwl)
csh
         if (odebug) call ma_print(dbl_mb(k_t),lenwl,lenij,'SO t-amp +')
csh
                   if (.not. ma_push_get(mt_dbl, 
     $              lenuv*lenwl, 'plus', l_p, k_p)) 
     $              call errquit('ma? lenuv*lenwl',lenuv*lenwl, MA_ERR)
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
     $              1.0d0,dbl_mb(k_omega),lenuv)
      if (odebug)  call ma_print(dbl_mb(k_omega),lenuv,lenwl,'UVIJ +-')
                   if (.not. ma_pop_stack(l_p))
     $              call errquit('ma_pop? l_p',l_p, MA_ERR)
                  end if
csh
       if (odebug) call ma_print(dbl_mb(k_omega),lenuv,lenij,'Omega')
csh
                  if (symw .eq. syml) then
csh
                   if (odebug) write(*,*) 'accumulating in ...',
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
                  if (.not. ma_pop_stack(l_t))
     $             call errquit('ma_pop? l_t',l_t, MA_ERR)
                 end if
                end do
                if (.not. ma_pop_stack(l_omega))
     $           call errquit('ma_pop? l_omega',l_omega, MA_ERR)
               end if
              end do
             end if
            end do
           end do
c
           if (.not. ma_pop_stack(l_gy))
     $      call errquit('ma_pop?',l_gy, MA_ERR)
           if (.not. ma_pop_stack(l_gx))
     $      call errquit('ma_pop?',l_gx, MA_ERR)
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
#include "errquit.fh"
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
      integer joff(nbf,0:7)
      integer boff(nbf,0:7)
      integer symvoff(0:7,0:7)
      integer symboff(0:7,0:7)
      integer ij,jb,ia
      integer ijlo,ijhi,jblo,jbhi
      integer i,j,a,b
      integer blo,bhi,vlo,alo,ilo,ihi,ulo,uhi
      integer dummy
csh
      logical odebug
      odebug = .false.
c
      if (odebug) then
      write(*,*) 'no = ',no(1),no(2)
      write(*,*) 'nv = ',nv(1),nv(2)
      write(*,*) 'no_sym = ',(no_sym(i,1),i=0,7),(no_sym(i,2),i=0,7)
      write(*,*) 'nv_sym = ',(nv_sym(i,1),i=0,7),(nv_sym(i,2),i=0,7)
      write(*,*) 'o_sym(a) = ',(o_sym(1,i,1),i=0,7),(o_sym(2,i,1),i=0,7)
      write(*,*) 'o_sym(b) = ',(o_sym(1,i,2),i=0,7),(o_sym(2,i,2),i=0,7)
      write(*,*) 'v_sym(a) = ',(v_sym(1,i,1),i=0,7),(v_sym(2,i,1),i=0,7)
      write(*,*) 'v_sym(b) = ',(v_sym(1,i,2),i=0,7),(v_sym(2,i,2),i=0,7)
      write(*,*) 'bf_per_ir = ',(bf_per_ir(i),i=0,7)
      endif
csh
      spina = spini
      spinb = spinj
csh ***** We have to create an empty g_b_big
      dummy = 0
      do symjb = 0,7
      do symb = 0,7
      symj = ieor(symjb,symb)
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
      if (odebug) write(*,*) 'Size of g_b_big = ',dummy
      if (.not.ga_create(mt_dbl,dummy,1,'b_big',
     $ -1,1,g_b_big)) call errquit
     $ ('uccsdt_makeb: room for b_big?',dummy, GA_ERR)
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
         if (odebug) write(*,*) 'symij,j,joff = ',
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
         if (odebug) write(*,*) 'symub,b,boff = ',
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
        if (odebug) write(*,*) 'symub,symb,symboff = ',
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
        if (odebug) write(*,*) 'symuv,symv,symvoff = ',
     $  symuv,symv,symvoff(symv,symuv)
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
csh
      if (odebug) write(*,*) 'maxlenuv = ',maxlenuv
csh
c
c     Transformed MO coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf*nv(spinb), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*nv',nbf*nv(spinb), MA_ERR)
      call ga_get(g_part(spinb),1,nbf,no(spinb)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
      if (odebug) call ma_print(dbl_mb(k_c),nbf,nv(spinb),'MO coeff')
csh
c
c     Allocate temporary GA
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'tmp', 
     $ maxlenuv,-1,g_tmp)) call errquit
     $ ('uccsdt_makeb: room for tmp?',maxlenuv*lenij, GA_ERR)
c
c     Data parallel transform Omega(u,v,i,j)C(v,b) = tmp(u,b,i,j)
c     (tmp was allocated big enough to do this).
c
      call ga_distribution(g_omega,ga_nodeid(),ijlo,ijhi,dummy,dummy)
      ij = 0
      do symij = 0,7
       symuv = symij
       symub = symij
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
            if (odebug) write(*,*) 'symij,symj,j,i,symb = ',
     $                  symij,symj,j,i,symb
csh
            if (lenuv.gt.0 .and. lenb.gt.0) then
csh
             if (odebug) write(*,*) 'lenub,lenuv = ',lenub,lenuv
csh
             if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $        l_buf, k_buf)) 
     $        call errquit('ma? lenuv',lenuv, MA_ERR)
             if (.not. ma_push_get(mt_dbl, lenub, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenub',lenub, MA_ERR)
c
c     (uv) -> (ub) for b (virtual) & v (SO) of same irrep 
c
             call ga_get(g_omega,1,lenuv,ij,ij,
     $        dbl_mb(k_buf),lenu)
csh
             if (odebug) call ma_print(dbl_mb(k_buf),lenu,lenv,'UV')
csh
             call dgemm('n','n',lenu,lenb,lenv,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        dbl_mb(k_c),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
csh
             if (odebug) call ma_print(dbl_mb(k_tmp),lenu,lenb,'UB')
             if (odebug) write(*,*) 'ga_put in the range ... ',
     $       symboff(symb,symub)+1,'~',symboff(symb,symub)+lenub,
     $       ij,'~',ij
csh
             call ga_put(g_tmp,
     $        symboff(symb,symub)+1,symboff(symb,symub)+lenub,
     $        ij,ij,dbl_mb(k_tmp),lenu)
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp, MA_ERR)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf, MA_ERR)
            end if
           end do
          end if
         end do
        end do
       end do
      end do
csh
      if (odebug) write(*,*) '**************************'
      if (odebug) write(*,*) 'transformation v -> b done'
      if (odebug) write(*,*) '**************************'
csh
c
c     Perform quarter transformation from tmp(u,b,i,j) to t(a,b,i,j)
c     t(i,a,syma,j,b,symb,symjb)
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nv(spina), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*nv',nbf*nv(spina), MA_ERR)
      call ga_get(g_part(spina),1,nbf,no(spina)+1,nmo,
     $ dbl_mb(k_c),nbf)
csh
      if (odebug) call ma_print(dbl_mb(k_c),nbf,nv(spina),'MO coeff')
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
         if (odebug) write(*,*) 'size of g_b: lenia, lenjb = ',
     $   lenia,lenjb
csh
         if (.not.ga_create(mt_dbl,lenia,lenjb,'b', 
     $    lenia,-1,g_b)) call errquit
     $    ('uccsdt_makeb: room for b?',lenia*lenjb, GA_ERR)
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
             if (odebug) write(*,*) 'symjb,symb,b,j,syma = ',
     $                   symjb,symb,b,j,syma
csh
             if (lenu.gt.0) then
              if (.not. ma_push_get(mt_dbl, leniu, 'buf',
     $         l_buf, k_buf)) 
     $         call errquit('ma? leniu',leniu, MA_ERR)
              if (.not. ma_push_get(mt_dbl, lenia, 'tmp',
     $         l_tmp, k_tmp)) 
     $         call errquit('ma? lenia',lenia, MA_ERR)
c
c     (iu) -> (ia) for a (virtual) & u (SO) of same irrep 
c
              ij = joff(j,symij)
              do i = ilo,ihi
               ij = ij + 1
csh
               if (odebug) write(*,*) 'i = ',i
               if (odebug) write(*,*) 'get from range ... ',
     $          boff(b,symub)+1,'~',boff(b,symub)+lenu,
     $          ij,'~',ij,
     $          ' to k_buf+',lenu*(i-ilo)
csh
               call ga_get(g_tmp,
     $          boff(b,symub)+1,boff(b,symub)+lenu,
     $          ij,ij,dbl_mb(k_buf+lenu*(i-ilo)),lenu)
              end do
csh
              if (odebug) call ma_print(dbl_mb(k_buf),lenu,leni,'B(UI)')
csh
              call dgemm('t','n',leni,lena,lenu,
     $         1.0d0,dbl_mb(k_buf),lenu,
     $         dbl_mb(k_c),nbf,
     $         0.0d0,dbl_mb(k_tmp),leni)
csh
              if (odebug) call ma_print(dbl_mb(k_tmp),leni,lena,'B(IA)')
              if (odebug) write(*,*) 'writing to range ... ',
     $        ia+1,'~',ia+lenia,jb,'~',jb
csh
              call ga_put(g_b,ia+1,ia+lenia,jb,jb,
     $         dbl_mb(k_tmp),leni)
              ia = ia + lenia
              if (.not. ma_pop_stack(l_tmp))
     $         call errquit('ma_pop?',l_tmp, MA_ERR)
              if (.not. ma_pop_stack(l_buf))
     $         call errquit('ma_pop?',l_buf, MA_ERR)
             end if
            end do
           end if
          end do
         end do
c
c     Write one symmetry block of B to disk
c
         if (odebug) call ga_print(g_b)
         call sh_put_t2(g_b,symb,symjb,spini,spinj,g_b_big)
csh      if (.not. uccsdt_ampfile_write_t2(d_b, 
csh  $    spini, spina, spinj, spinb, symia, 
csh  $    blo, bhi, g_b))
csh  $    call errquit('write_t2 failed',0)
         if (.not. ga_destroy(g_b)) call errquit('GA 999',0, GA_ERR)
        endif
       end do
      end do
csh
      if (odebug) write(*,*) '**************************'
      if (odebug) write(*,*) 'transformation u -> a done'
      if (odebug) write(*,*) '**************************'
csh
      if (.not. ga_destroy(g_tmp)) call errquit('GA 999',0, GA_ERR)
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c, MA_ERR)
c
      return
      end
c $Id$
