      subroutine uccsdt_bterm()
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "global.fh"
c
      integer uqsuper(2,nw_max_shell_small), nuqsuper
      integer usub(2,nw_max_shell_small), nusub
      integer vsub(2,nw_max_shell_small), nvsub
      integer urange(2,0:7), vrange(0:7)
      integer l_sopairmap, k_sopairmap
c
      integer supermaxnbf, submaxnbf, shmaxnbf
      integer nij(0:7)
c
c     For now just hardwire the super-blocking for out-of-core solution
c     For production code needs to be determined by memory.
c
      supermaxnbf = min(200,nbf)
      call soints_uqshblock(1, sh_n_uq, supermaxnbf, uqsuper, nuqsuper)
      if (.not. bas_nbf_cn_max(basis, shmaxnbf)) call errquit
     $     ('uccsdt_bterm: bas_nbf_cn_max?', 0)
c
c     For the subblocking determine the absolute minimum size
c     and then try to increase it subject to available parallelism
c     (aim for roughly 10 tasks per processor).
c
c     Also need to take into account memory per node ... NOT DONE YET!!!!!
c
      submaxnbf = 0
      do ishuq = 1, sh_n_uq
         submaxnbf = max(submaxnbf, sh_uq_bf(2,shuq)-sh_uq_bf(1,shuq)+1)
      end do
      submaxnbf = max(submaxnbf, submaxnbf/sqrt(20*P))
      write(6,*) ' super/sub maxnbf ', supermaxnbf, submaxnbf
c
c     Count the no. of symmetrized pairs of occupied orbitals per irrep
c
      do symij = 0,7
         nij(symij) = 0
         do symj = 0,7
            symi = ieor(symij,symj)
            if (symi .gt. symj) then
               nij(symij) = nij(symij) + no_sym(symi)*no_sym(symj)
            else if (symi .eq. symj) then
               nij(symij) = nij(symij) + 
     $              no_sym(symi)*(no_sym(symi)+1)/2
            end if
         end do
      end do
      write(6,*) ' #symmetrized pairs of occupied orbitals per irrep '
      write(6,*) (nij(symij),symij=0,7)
c
c     Generate addressing & map info for the SO pairs
c
      if (.not. ma_push_get(mt_int, nbf*(nbf+1)/2, 'sopairmap',
     $     l_sopairmap, k_sopairmap)) call errquit
     $     ('uccsdt_bterm: insufficient memory for sopairmap?',0)
c
      open the files for the SO amplitudes and result vector
c
c     Generate taa(ij,uv), tbb(uv,ij), tabp(uv,ij), tabm(uv,ij)
c     for i>=j, u>=v (zero diagonal if antisymmetric).
c
      call uccsdt_tosopair(file info)
c
      
      do usuper = 1, nuqsuper
         u_uqlo = uqsuper(1,usuper)
         u_uqhi = uqsuper(2,usuper)
         call soints_uqshblock(u_uqlo, u_uqhi, submaxnbf, usub, n_usub)
*         call soints_shblock_info(basis, u_uqlo, u_uqhi, urange)

         do vsuper = 1, usuper
            v_uqlo = uqsuper(1,vsuper)
            v_uqhi = uqsuper(2,vsuper)
            call soints_uqshblock(v_uqlo, v_uqhi, submaxnbf,vsub,n_vsub)
*            call soints_shblock_info(basis, v_uqlo, v_uqhi, vrange)
c
c     Transform the amplitudes back to the hole SO basis
c     and symmetry adapt the external pairs, storing them
c     blocked on disk for ready access.
c
         end do
      end do
c
      do usuper = 1, nuqsuper
         u_uqlo = uqsuper(1,usuper)
         u_uqhi = uqsuper(2,usuper)
         call soints_uqshblock(u_uqlo, u_uqhi, submaxnbf, usub, n_usub)
*         call soints_shblock_info(basis, u_uqlo, u_uqhi, urange)

         do vsuper = 1, usuper
            v_uqlo = uqsuper(1,vsuper)
            v_uqhi = uqsuper(2,vsuper)
            call soints_uqshblock(v_uqlo, v_uqhi, submaxnbf,vsub,n_vsub)
*            call soints_shblock_info(basis, v_uqlo, v_uqhi, vrange)
c
            read amplitudes in super-blocks
c
            do parallel blocks of uv >= wl
               SO integrals
               symmetrize integrals
               local MxM
               accumulate result
            end do
c
            accumulate results to disk
c
         end do
c
         accumulate results to disk
c
      end do
c
      de-antisymmetrize and form X, G, B in the MO basis storing on disk
c
      end do

      end
      subroutine uccsdt_sopair_info(nuqsuper, uqsuper, 
     $     submaxnbf, sopairmap, )
      implicit none
#include "cuccsdtP.fh"
      integer nuqsuper          ! [input]
      integer uqsuper(2,*)      ! [input]
      integer submaxnbf         ! [input] 
      integer sopairmap(*)      ! [output]
c
c     Analyze the SO pair basis arising from each pair of blocks of
c     symmetry unique AO shells
c     ... how many pairs of functions in each irrep
c     ... map from natural SO pair index to symmetry packed pair index
c     ... file addressing information.
c
c     Antisymmetric pairs are stored the same as symmetric ones but
c     with a zero diagonal (to unify addressing).
c
c     The symmetry packed SO pairs are stored in a complicated blocking
c     as follows.
c     
c     .   For each superblock of u>=v pairs
c     .      For each subblock of u>=v pairs
c     .          For each symuv -> symij
c     .            For each symu -> symv
c     .               For u of symu
c     .                  For v<=u of symv
c
c     On disk we have taa(ij,uv), tbb(uv,ij), tabp(uv,ij), tabm(uv,ij)
c     for i>=j, u>=v (zero diagonal if antisymmetric).  Also result vector.
c     File addressing is as follows
c     
c     For each of taa(ij,uv), tbb(uv,ij), tabp(uv,ij), tabm(uv,ij)
c     .   For each superblock of u>=v pairs
c     .      For each subblock of u>=v pairs
c     .          For each symuv -> symij
c     .            For each symu -> symv
c     .               For u of symu
c     .                  For v<=u of symv
c     .                     all ij of symij
c
      
      end
      subroutine uccsdt_tosopair(spini, spinj, file info, ijinfo
      implicit none
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
c
c     Transform the amplitudes from the MO to the SO basis
c     and symmetrize if necessary (mixed spin).  Write the
c     transformed amplitudes to disk.
c
c     This version assumes it can hold all of one spin-block
c     of the amplitudes in core
c
c     Out-of-core needs a little more work.
c
      integer lenuv, lenub, lenij, boff(nw_max_nbf_small), 
     $     voff(nw_max_nbf_small)
c
      spina = spini
      spinb = spinj
c
c     Addressing & sizes for the partially transformed amplitudes
c
      lenub = 0
      do symub = 0, 7
         do symb = 0, 7
            symu = ieor(symub,symb)
            do b = v_sym(1,symb,spinb),v_sym(2,symb,spinb)
               boff(b) = lenub
               lenub = lenub + bf_per_ir(symu)
            end do
         end do
      end do
c

      boff ?????

      joff ?????


      lenij = 0
      do symij = 0, 7
         do symj = 0, 7
            symi = ieor(symij,symj)
            do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
               joff(j) = lenij
               lenij = lenij + no_sym(symi)
            end do
         end do
      end do
c
      lenuv = 0
      do symuv = 0, 7
         do symv = 0, 7
            symu = ieor(symuv,symv)
            do v = bf_per_ir_cum(symv)+1,
     $           bf_per_ir_cum(symv)+bf_per_ir(symv)
               voff(v) = lenuv
               lenuv = lenuv + bf_per_ir(symu)
            end do
         end do
      end do
c
      write(6,*) ' LENUV, LENUB, LENIJ ', lenuv, lenub, lenij
c
      if (.not. ga_create(mt_dbl, lenuv, lenij, 'sopair', 
     $     lenub, 0, g_tmp)) call errquit
     $     ('uccsdt_tosopair: room for tmp?', 0)
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
            if (.not. uccsdt_ampfile_read_t2(d_amp, 
     $           spini, spina, spinj, spinb, symia, 
     $           blo, bhi, g_t2, .true., 'column'))
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
                     do syma = 0, 7
                        symi = ieor(symia,syma)
                        symu = syma
                        leni = no_sym(symi,spini)
                        lena = nv_sym(syma,spina)
                        lenu = bf_per_ir(syma)
                        lenia = leni*lena
                        leniu = leni*lenu
                        alo = v_sym(1,syma,spina)
                        ulo = bf_per_ir_cum(symu) + 1
                        uhi = ulo + lenu - 1
                        if (lenia.gt.0 .and. lenu.gt.0) then
                           if (.not. ma_push_get(mt_dbl, lenia, 'buf',
     $                          l_buf, k_buf)) 
     $                          call errquit('ma? lenia',lenia)
                           if (.not. ma_push_get(mt_dbl, leniu, 'tmp',
     $                          l_tmp, k_tmp)) 
     $                          call errquit('ma? leniu',leniu)
c
c     (ia) -> (iu) for a (virtual) & u (SO) of same irrep 
c
                           call ga_get(g_t2, ia+1, ia+lenia, jb, jb,
     $                          dbl_mb(k_buf), leni)
                           call dgemm('n','t',leni,lenu,lena,
     $                          1d0,dbl_mb(k_buf),leni,SO(ulo,alo),nbf,
     $                          0d0,dbl_mb(k_tmp),leni)
                           do i = o_sym(1,symi,spini),
     $                          o_sym(2,symi,spini)
                              ij = joff(j) + (i - o_sym(1,0,spini) + 1)
                              call ga_put(g_tmp,ulo+boff(b),uhi_boff(b),
     $                             ij,ij,dbl_mb(k_tmp),1)
                           end do
                        end if
                     end do
                  end if
               end do
            end do
c
            if (.not. ga_destroy(g_t2)) call errquit('GA 999',0)
         end do
      end do
c
c     Now have t(u,b,i,j) ... data parallel transform b in place 
c     (tmp was allocated big enough to do this).
c
      call ga_distribution(g_tmp, ga_nodeid(), ijlo, ijhi)
      ij = 0
      do symij = 0, 7
         do symj = 0, 7
            symi = ieor(symij,symj)
            do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
               do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
                  ij = ij + 1
                  if (ij.ge.ijlo .and. ij.le.ijhi) then
                     call ga_get(
