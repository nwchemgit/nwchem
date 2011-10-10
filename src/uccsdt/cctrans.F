      subroutine cctrans()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
c     Generate the following lists of integrals.  Each list
c     is stored just like the amplitudes, but with SO's 
c     instead of virtuals.  Note lists 13 and 14 contribute
c     to lists 9 and 12, and that list 7 is not used (nor made).
c
c     i,j = particle-transformed occupied orbitals
c     m,n = hole-transformed occupied orbitals
c     u,v = symmetry-adapted occupied orbitals
c
c     u - the SOs arising from the unique AO shells in the
c     .   range ushuqlo, ushuqhi
c     All other indices span the complete range.
c
c     #  Integral        Spins           Storage  Used for 
c     -- --------------- --------------- -------- --------------
c     1. (ui|vj)-(uj|vi) i=alpha j=alpha I(jv,iu) Z1,Z6
c     2. (ui|vj)         i=alpha j=beta           Z2,Z7
c     3. (ui|vj)         i=beta  j=alpha          Z5,Z8
c     4. (ui|vj)-(uj|vi) i=beta  j=beta           Z3,Z4
c
c     5. (um|vn)-(un|vm) m=alpha n=alpha I(nv,mu) aaA and F
c     6. (um|vn)         m=alpha n=beta           abA and F
c    *7. (um|vn)         m=beta  n=alpha          *Not used
c     8. (um|vn)-(un|vm) m=beta  n=beta           bbA and F
c
c     9. (ui|vn)-(uv|in) i=alpha n=alpha I(nv,iu) aaC1,abC4
c     10.(ui|vn)         i=alpha n=beta           abC1,bbC2
c     11.(ui|vn)         i=beta  n=alpha          abC5,aaC2
c     12.(ui|vn)-(uv|in) i=beta  n=beta           abC1,abC2
c
c     13.(uv|in)         i=alpha n=alpha I(nv,iu) abC3
c     14.(uv|in)         i=beta  n=beta           abC5
c
c     The occ-SO pairs are stored as follows
c
c     fill oso_off with -99999999
c     ind = 0
c     do symiu
c     .  do u in natural order
c     .     -> symu and symi
c     .     oso_off(u,symiu) = ind
c     .     do i of symi
c     .        pair(1+ind) = T(i,u)
c     .        ind = ind + 1
c     
c     Can also address pairs as
c     pair(1 + i-o_sym(1,symi,spini) + oso_off(u,symiu)) = T(i,u)
c
c     The lists will be used in large matrix multiplications
c     as follows:
c
c     Z(iu,nw) = <ij|uv>*T(jv,nw)
c     do sym(iu) -> sym(jv) -> sym(nw)
c     .   Read T, allocate Z
c     .   Z(iu,nw) <- I(jv,iu)*T(jv,nw)
c     end do
c
c     For best performance we need the symmetry sub-blocks distributed
c     across the whole machine as separate dense arrays.  Thus, each
c     list is stored with a separate GA for each pair symmetry.
c     
c     Need the occ*SO addressing arrays ... IN WHAT SO LABELLING?
c   
c     The initial parallel algorithm is not entirely satisfactory
c     since in order to keep memory feasible it does a lot of
c     communication.  It also may be complicated to balance
c     available parallelism and memory usage without doing
c     even more communication.
c
c     do blocks of u 
c     .  do blocks of w 
c     .     goal is to make (uw|vi) and (uw|vn) both alpha and beta
c     .     next 4 loops in parallel
c     .     do u in block
c     .        do w in block
c     .           do blocks of v
c     .              do blocks of l
c     .                 make (uw|vl)
c     .                 -> (uw|vi) and (uw|li) (ditto for n)
c     .     sync
c     .     transform the second index
c     
c     Input is
c
c     ush_uqlo, ush_uqhi - range of symmetry unique AO shells for U
c     
c     output is
c
c     14*nirrep global arrays
c
      integer urange(2,0:7), wrange(2,0:7)
      integer i, ir
      integer lists(0:7,14)
      integer nush_block, ush_block(2,nw_max_shells_small)
      integer nwsh_block, wsh_block(2,nw_max_shells_small)

      integer nxtask
      external nxtask
c
c     Someone has already made the addressing lists needed 
c     here for the results, and in the C-term and B-term
c
c     iu_len(symiu,spini) -> #pairs iu of given symmetry and spin with u
c     an SO in the second labelling scheme restricted to the range of u
c
c     iu_off(u,symiu,spini) -> offset to (i,u) with u an SO in the
c     second labelling scheme restricted to the range of u
c
c     iso_len(symiv,spini)-> #pairs iv for v any SO in the second
c     labelling scheme over all SOs
c
c     iso_off(v,symiv,spini) -> offset to (i,v) with v any SO in the
c     second labelling scheme over all SOs
c
c     Allocate the global arrays for the results.
c     How should they be best distributed?
c
      call ifill(14*8, -1, lists, 1)
      do symiu = 0,nir-1
         lena = iso_len(symiu,0)
         lenb = iso_len(symiu,1)
         ulena = iu_len(symiu,0)
         ulenb = iu_len(symiu,1)
         if (.not. ga_create(mt_dbl, lena, ulena, 'list 1', 0, 0, 
     $        lists(symiu,1)) 
     $        call errquit('cctrans: failed allocating list', 1, GA_ERR)
         if (.not. ga_create(mt_dbl, lenb, ulena, 'list 2', 0, 0, 
     $        lists(symiu,2)) 
     $        call errquit('cctrans: failed allocating list', 2, GA_ERR)
         if (.not. ga_create(mt_dbl, lena, ulenb, 'list 3', 0, 0, 
     $        lists(symiu,3)) 
     $        call errquit('cctrans: failed allocating list', 3, GA_ERR)
         if (.not. ga_create(mt_dbl, lenb, ulenb, 'list 4', 0, 0, 
     $        lists(symiu,4)) 
     $        call errquit('cctrans: failed allocating list', 4, GA_ERR)
         if (.not. ga_create(mt_dbl, lena, ulena, 'list 5', 0, 0, 
     $        lists(symiu,5)) 
     $        call errquit('cctrans: failed allocating list', 5, GA_ERR)
         if (.not. ga_create(mt_dbl, lenb, ulena, 'list 6', 0, 0, 
     $        lists(symiu,6)) 
     $        call errquit('cctrans: failed allocating list', 6, GA_ERR)
c$$$         if (.not. ga_create(mt_dbl, lena, ulenb, 'list 7', 0, 0, 
c$$$     $        lists(symiu,7)) 
c$$$     $        call errquit('cctrans: failed allocating list', 7)
         if (.not. ga_create(mt_dbl, lenb, ulenb, 'list 8', 0, 0, 
     $        lists(symiu,8)) 
     $        call errquit('cctrans: failed allocating list', 8, GA_ERR)
         if (.not. ga_create(mt_dbl, lena, ulena, 'list 9', 0, 0, 
     $        lists(symiu,9)) 
     $        call errquit('cctrans: failed allocating list', 9, GA_ERR)
         if (.not. ga_create(mt_dbl, lenb, ulena, 'list 10', 0, 0, 
     $        lists(symiu,10)) 
     $        call errquit('cctrans: failed allocating list', 10,
     &       GA_ERR)
         if (.not. ga_create(mt_dbl, lena, ulenb, 'list 11', 0, 0, 
     $        lists(symiu,11)) 
     $        call errquit('cctrans: failed allocating list', 11,
     &       GA_ERR)
         if (.not. ga_create(mt_dbl, lenb, ulenb, 'list 12', 0, 0, 
     $        lists(symiu,12)) 
     $        call errquit('cctrans: failed allocating list', 12,
     &       GA_ERR)
         if (.not. ga_create(mt_dbl, lena, ulena, 'list 13', 0, 0, 
     $        lists(symiu,13)) 
     $        call errquit('cctrans: failed allocating list', 13,
     &       GA_ERR)
         if (.not. ga_create(mt_dbl, lenb, ulenb, 'list 14', 0, 0, 
     $        lists(symiu,14)) 
     $        call errquit('cctrans: failed allocating list', 14,
     &       GA_ERR)
      end do
      do i = 1, 14
         do symiu = 0, nir-1
            if (lists(symiu,i) .ne. -1) call ga_zero(lists(symiu,i))
         end do
      end do
c
c     Block the u and w shells to control global memory.  For maximum
c     parallelism we want just one block of each.  For now force to one
c     block
c
      nush_block = 1
      ush_block(1,1) = ush_uqlo
      ush_block(2,1) = ush_uqhi
      nwsh_block = 1
      wsh_block(1,1) = 1
      wsh_block(2,1) = sh_n_uq
c
      do ush_blk = 1, nush_block
         ush_blk_lo = ush_block(1,ush_blk)
         ush_blk_hi = ush_block(2,ush_blk)
         call soints_shblock_info(basis, ush_blk_lo, ush_blk_hi, urange)
         write(6,*) ' urange ', ((urange(i,ir),i=1,2),ir=0,nir-1)
c
         do wsh_blk = 1, nwsh_block
            wsh_blk_lo = wsh_block(1,wsh_blk)
            wsh_blk_hi = wsh_block(2,wsh_blk)
            call soints_shblock_info(basis, wsh_blk_lo, wsh_blk_hi, 
     $           wrange)
            write(6,*) ' wrange ', ((wrange(i,ir),i=1,2),ir=0,nir-1)
c
c     For u and w in the given ranges of symmetry unique shells (-> SOs)
c     we will make <vw|iu> for i an occupied orbital
c     ({hole,particle}*{alpha,beta} -> 4 sets) and a SOs v in the second
c     or symmetry blocked labelling.  Since we end up doing lots of
c     small accumulates into this temporary array (to limit local memory
c     consumption) we may eventually bundle all 4 sets into one array,
c     but for now we have one array per set.
c
c     The symmetry of vw pairs is predetermined by symiu = symvw
c     and they are stored as follows.
c     do symw -> symv
c     .   do w in irrep
c     .      do v in irrep
c     
c     Figure out the size of the global arrays <vw|iu>
c
            lena = 0
            lenb = 0
            do symiu = 0, nir-1
               symvw = symiu
               do symw = 0, nir-1
                  symv = ieor(symvw,symw)
                  dimw = wrange(2,symw)-wrange(1,symw)+1
                  dimv = bf_per_ir(symv)
                  lena = lena + dimv*dimw*iu_len(symiu,1)
                  lenb = lenb + dimv*dimw*iu_len(symiu,2)
               end do
            end do
            if (.not. ga_create(mt_dbl, lena, 1, '<vw|iu>1', 0, 0, 
     $           g_vwiu(1)) 
     $           call errquit('cctrans: failed allocating <vw|iu>', 1,
     &       GA_ERR)
            if (.not. ga_create(mt_dbl, lena, 2, '<vw|iu>2', 0, 0, 
     $           g_vwiu(2)) 
     $           call errquit('cctrans: failed allocating <vw|iu>', 1,
     &       GA_ERR)
            if (.not. ga_create(mt_dbl, lenb, 3, '<vw|iu>3', 0, 0, 
     $           g_vwiu(3)) 
     $           call errquit('cctrans: failed allocating <vw|iu>', 3,
     &       GA_ERR)
            if (.not. ga_create(mt_dbl, lenb, 4, '<vw|iu>4', 0, 0, 
     $           g_vwiu(4)) 
     $           call errquit('cctrans: failed allocating <vw|iu>', 4,
     &       GA_ERR)
c
c     Need to determine the blocksize for the v and l loops again to
c     control local memory and increase parallelism ... for now hardwire
c     to 2 for debug.
c
c     Next four loops run in parallel.
c
            vsh_blk_size = 2
            count = 0
            next = nxtask(nproc, 1)
            do ush = ush_blk_lo, ush_blk_hi
               do wsh = wsh_blk_lo, wsh_blk_hi
                  do vsh_blk_lo = 1, sh_n_uq, vsh_blk_size
                     vsh_blk_hi = 
     $                    min(vsh_blk_lo+vsh_blk_size-1,sh_n_uq)
                     do lsh_blk_lo = 1, vsh_blk_lo, vsh_blk_size ! v>=l
                        lsh_blk_hi = 
     $                       min(lsh_blk_lo+vsh_blk_size-1,sh_n_uq)
                        if (next .eq. count) then
                           call cctrans2(ush, wsh,
     $                          vsh_blk_lo, vsh_blk_hi,
     $                          lsh_blk_lo, lsh_blk_hi,
     $                          info about <vw|iu>,
     $                          addressing arrays, 
     $                          )
                           next = nxtask(nproc, 1)
                        end if
                        count = count + 1
                     end do
                  end do
               end do
            end do
c
c     Barrier
c
            next = nxtask(-nproc, 1)
            call ga_sync()
c
c     Data parallel transformation of v or w
c
            do iset = 1, nset
               do symiu = 0, nir-1
                  call ga_distribution(g_vwiu(iset), me, 
     $                 vwlo, vwhi, iulo, iuhi)
                  iu = 1
                  do symiu = 0, nir-1
                     do symu = 0, nir-1
                        symi = ieor(symiu,symu)
                        do u = urange(1,symu),urange(2,symu)
                           do i = ilo(symi,iset),ihi(symi,iset) ??? irange ??
                              if (iu.ge.iulo .and. iu.le.iuhi) then
                                 ----> get (vw|iu)
                                 ----> transform and accumulate results
                              end if
                              iu = iu + 1
                           end do
                        end do
                     end do
                  end do
               end do
            end do
c
            call ga_sync()
c
            free up <vw|iu> and other work space
c            
         end do                 ! blocks of wsh
      end do                    ! blocks of ush
c
c     Done
c
      end
      subroutine cctrans2(ush, wsh,
     $     vsh_blk_lo, vsh_blk_hi,
     $     lsh_blk_lo, lsh_blk_hi,
     $     ) 
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "global.fh"
#include "cuccsdtP.fh"
      integer ush, vsh, vsh_blk_lo, vsh_blk_hi, lsh_blk_lo, lsh_blk_hi
c     
c     For given ranges of unique AO shells (-> range of SOs) with
c     restriction vsh>=lsh, transform (uw|vl) -> (uw|vi) (for up to 4
c     different i's) accumulating results into the global array.
c     
c     Have to further sub-block vsh and lsh in order to keep memory
c     tractable, but we should have enough room to accumulate (uw|vi)
c     and (uw|li) before adding into the global array (if we don't
c     do this there will be too much communication).  The local
c     result is stored as I(v,i,w,u) and I(l,i,w,u).
c     
      integer nvblock, vblock(2,nw_max_shells_small)
      integer nlblock, lblock(2,nw_max_shells_small)
      integer ninblock, lo
      integer vblk, lblk
c     
c     Block the shells so that there are roughly no more than XX
c     functions in each block (better would be to look at the memory
c     available)
c     
      call cctrans_block(nvblock, vblock, vsh_blk_lo, vsh_blk_hi)
      call cctrans_block(nlblock, lblock, lsh_blk_lo, lsh_blk_hi)
c     
      allocate memory for the temporary stuff
c     
      do vblk = 1, nvblock
         lblktop = nlbock
         if (vsh_blk_lo .eq. lsh_blk_lo) lblktop = vblk
         do lblk = 1, nlblock
            call cctrans3(ush, wsh, 
     $           vblock(1,vblk), vblock(2,vblk),
     $           lblock(1,lblk), lblock(2,lblk),
            ...)
         end do
      end do
c     
c     Accumulate to the GA ... this is the bit that is doing too
c     much communication if this routine is called more than once
c     due to exterior blocking of v and l.   The four target
c     GAs could be combined ---- 4x message no. reduction
c     
      ind = 0
      do symwu = 0, nir-1
         symvi = symwu
         do symu = 0, nir-1
            symw = ieor(symwu,symu)
            do u = urange(1,symu), urange(2,symu)
               do w = wrange(1,symw), wrange(2,symw)
                  do symi = 0, nir-1
                     symv = ieor(symvi,symi)
                     vlo = vrange(1,symv)
                     vhi = vrange(2,symv)
                     vdim = vhi - vlo + 1
                     if (vdim .gt. 0) then
                        do iset = 1, nset ! Inside in prep for bundling
                           ptr = ind + k_viwu(iset)
                           do i = ilo(symi,iset),ihi(symi,iset)
c
c     GA is <vw|iu> stored linearly.
c
                              vwiu = vlo
                              iu = 1 + i-ilo(symi,iset) + 
     $                             iu_off(u,symiu,spini)
                              vw = 1 + 
                           end do
                        end do
                     end if
c
                     if (v_blk_lo != l_blk_lo) ditto for l block
c
                  end do
               end do
            end do
         end do
      end do
                           
c     
      free locals
c     
c     Done
c     
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
            vblocks(1,nvblock) = lo
            vblocks(2,nvblock) = vsh
            lo = vsh+1
         end if
      end do
c
      end
      subroutine cctrans3(ush, wsh, 
     $     vsh_blk_lo, vsh_blk_hi, lsh_blk_lo, lsh_blk_hi,
     $     fullvrange, fulllrange, irange, nset, c, ldc, 
     $     k_viwu, k_liwu) 
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "global.fh"
#include "cuccsdtP.fh"
      integer ush, vsh, vsh_blk_lo, vsh_blk_hi, lsh_blk_lo, lsh_blk_hi
      integer fullvrange(2,0:7), fulllrange(2,0:7)
      integer nset
      integer irange(2,0:7,nset)
      double precision c(nbf,ldc,nset)
      integer k_viwu(nset), k_liwu(nset)
c     
c     Actually get to do something in this routine.
c     
c     For symmetry unique shells u and w, and ranges of v and l
c     compute the integrals and transform v and l to i.
c     
c     k_viwu(iset) is the offset in dbl_mb for (vi|wu) which
c     is stored as (ditto for (li|wu))
c     do symwu
c     .  symvi = symwu
c     .  do symu
c     .     -> symw
c     .     do u = SOs of symu from u-shell
c     .        do w = SOs of symw from w-shell
c     .           do symi
c     .              -> symv
c     .              do i occupied orbitals of symi
c     .                 do v of symv from v-shell-block
c     .                    (vi|wu)
c     
c     The lambda (l) arrays will only be referenced if the 
c     nu (v) and lambda blocks are not the same.
c     
      integer
     $     loff(nw_max_nbf_small),     lrange(2,0:7), 
     $     voff(nw_max_nbf_small,0:7), vrange(2,0:7),
     $     woff(nw_max_nbf_small,0:7), wrange(2,0:7),
     $     uoff(nw_max_nbf_small,0:7), urange(2,0:7),
     $     lenlv(0:7), lenwu(0:7), leng, symwuoff(0:7)
      integer lvwu, symwu, symu, symw, u, w, symv, syml, symi,
     $     vlo, vhi, vdim, llo, lhi, ldim, 
     $     fullvlo, fullvhi, fullvdim, fullllo, fulllhi, fullldim,
     $     iset, k_g, l_g
      integer vioff(0:7,0:7,4), lioff(0:7,0:7,4)
c     
c     Compute g(l,v,w,u) = (lv|wu) 
c     
      call soints_block_info(basis, nbf,
     $     lsh_blk_lo, lsh_blk_hi, loff, lrange,
     $     vsh_blk_lo, vsh_blk_hi, voff, vrange,
     $     wsh, wsh, woff, wrange,
     $     ush, ush, uoff, urange,
     $     lenlv, lenwu, symwuoff, leng)
c     
      if (.not. ma_push_get(mt_dbl, leng, 'cctrans3', l_g, k_g))
     $     call errquit('cctrans3: no room for g', leng, MA_ERR)
c     
      call soints_block(basis, nbf,
     $     lsh_blk_lo, lsh_blk_hi, loff, lrange,
     $     vsh_blk_lo, vsh_blk_hi, voff, vrange,
     $     wsh, wsh, woff, wrange,
     $     ush, ush, uoff, urange,
     $     1d-12, dbl_mb(k_g))
c
c     Convenience adressing arrays
c
      li = 0
      vi = 0
      do symvi = 0, nir-1
         do symi = 0, nir-1
            symv = ieor(symvi,symi)
            vioff(symv,symi) = vi
            lioff(symv,symi) = li
            vi = vi + fullvdim(symv)*idim(symi,iset)
            li = li + fullldim(symv)*idim(symi,iset)
      end do
c     
c     Transform v and l to i.  Loop thru the SO integrals in the order
c     generated.
c     
      do iset = 1, nset
         viwu = k_viwu(iset)
         liwu = k_liwu(iset)
         lvwu = k_g
         do symwu = 0, nir-1
            symlv = symwu
            do symu = 0, nir-1
               symw = ieor(symwu,symu)
               do u = urange(1,symu), urange(2,symu)
                  do w = wrange(1,symw), wrange(2,symw)
                     do symv = 0, nir-1
                        syml = ieor(symlv,symv)
c     Info about the SO integrals just computed
                        vlo = vrange(1,symv)
                        vhi = vrange(2,symv)
                        vdim = vhi - vlo + 1
                        llo = lrange(1,syml)
                        lhi = lrange(2,syml)
                        ldim = lhi - llo + 1
c     Info about the results being accumulated
                        fullvlo = fullvrange(1,symv)
                        fullvhi = fullvrange(2,symv)
                        fullvdim = fullvhi - fullvlo + 1
                        fullllo = fulllrange(1,syml)
                        fulllhi = fulllrange(2,syml)
                        fullldim = fulllhi - fullllo + 1
                        if (vdim*ldim .gt. 0) then
c     (lv|wu) -> (li|wu) ... v -> i
                           ptr = liwu + lioff(syml,symi)+llo-fullllo
                           call dgemm('n','n',
     $                          ldim,idim(symi,iset),vdim,
     $                          1.0d0,dbl_mb(lvwu),ldim,
     $                          c(vlo,ilo(symi,iset),iset),nbf,
     $                          1.0d0,dbl_mb(ptr),fullldim)
c     (lv|wu) -> (vi|wu) ... transpose and l -> i
                           ptr = viwu + vioff(symv,symi)+vlo-fullvlo

                           need to do a little more to get the correct
                           results address for l since are not going
                           thru the l symmetry in order !!!!

                           symi = syml
                           call dgemm('t','n',
     $                          vdim,idim(symi,iset),ldim,
     $                          1.0d0,dbl_mb(lvwu),ldim,
     $                          c(llo,ilo(symi,iset),iset),nbf,
     $                          1.0d0,dbl_mb(ptr),fullvdim)
                        end if
                        lvwu = lvwu + vdim*ldim
                     end do
                     do symi = 0, nir-1
                        symv = ieor(symvi,symi)
                        viwu = viwu + fullvdim(symv)*idim(symi,iset)
                        liwu = liwu + fullldim(symv)*idim(symi,iset)
                     end do
                  end do
               end do
            end do
         end do
      end do
c     
      if (.not. ma_pop_stack(l_g)) call errquit
     $     ('cctrans3: freeing g?', 0, MA_ERR)
c     
      end
c $Id$
