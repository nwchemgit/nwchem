      subroutine uccsdt_do_triples_dummy_amplitudes(
     $     e_sym, eeps,
     $     ncblock, cblock, ceps,
     $     niblock, iblock, ieps,
     $     nkblock, kblock, keps,
     $     moints, nbf, k_t2, spine, spinf,
     $     pointers, g_list, phase, oasym)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer e_sym(2,0:7)
      integer ncblock, cblock(3,*)
      integer niblock, iblock(3,*)
      integer nkblock, kblock(3,*)
      integer nbf, spine, spinf
      double precision eeps(*), ceps(*), ieps(*), keps(*)
      double precision moints(nbf, nbf, nbf, nbf)
      integer pointers(2,ncblock,niblock,nkblock)
      integer g_list, k_t2
      double precision phase
      logical oasym, flip
c
      integer cblk, iblk, kblk, c, i, k, e, syme, ptr
      double precision t
c
c     Fill in the dummy T's using the following template
c
c     14.  t(e,c,i,k) spin(e)=spin(a) mixed
c
c     THIS IS CURRENTLY MOSTLY HARDWIRED TO CLOSED SHELL
c
      flip = (spine.ne.spinf).and.(spine.ne.1)
c      write(6,*)' k_t2 ', k_t2
      do kblk = 1, nkblock
         do iblk = 1, niblock
            do cblk = 1, ncblock
               ptr = pointers(1,cblk,iblk,kblk)
               syme = ieor(ieor(cblock(3,cblk),iblock(3,iblk)),
     $              kblock(3,kblk))
               do k = kblock(1,kblk),kblock(2,kblk)
                  do i = iblock(1,iblk),iblock(2,iblk)
                     do c = cblock(1,cblk),cblock(2,cblk)
                        do e = e_sym(1,syme),e_sym(2,syme)
                           if (flip)then
                              t = moints(c,e,k,i)
                           else
                              t = moints(e,c,i,k)
                           endif
c                          if (oasym) t = t - moints(e,c,k,i)
c                          t = t / (eeps(e)+ceps(c)-ieps(i)-keps(k))
c                          t = t * phase
c                          call uccsdt_check_t2(dbl_mb(k_t2),nbf,
c    &                                          e,c,i,k,t,oasym)
                           call ga_put(g_list,ptr,ptr,1,1,t,1)
                           ptr = ptr + 1
                        end do
                     end do
                  end do
               end do
               if (ptr .ne. (pointers(1,cblk,iblk,kblk)+
     $              pointers(2,cblk,iblk,kblk))) call errquit
     $              ('bad length for t record',ptr, DISK_ERR)
            end do
         end do
      end do
c
      end
      subroutine uccsdt_check_t2(
     $     buf, nbf, e, c, i, k, t, oasym)
      implicit none
      integer nbf, e, c, i, k
      double precision buf(nbf,nbf,nbf,nbf), t
      logical oasym
      if (abs(buf(e,c,i,k)-t).gt.1e-6) then
         write(6,1) e,c,i,k,t,buf(e,c,i,k),buf(e,c,k,i),
     &                        buf(c,e,i,k),buf(c,e,k,i),oasym
 1       format(4i3,5(1x,f12.7),1x,l1)
      endif
      return
      end
      subroutine uccsdt_do_check_list(
     $     buf,
     $     ilo, ihi,
     $     jlo, jhi,
     $     klo, khi,
     $     llo, lhi,
     $     correct, nbf, oasym)
      implicit none
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi, nbf
      double precision buf(ilo:ihi,jlo:jhi,klo:khi,llo:lhi)
      double precision correct(nbf,nbf,nbf,nbf)
      logical oasym
c
c     Print out integrals that differ from the correct ones
c     
c     oasym = .true. if correct integral should be antisymmetrized
c
      integer i, j, k, l
      double precision ok
c
      do l = llo, lhi
         do k = klo, khi
            do j = jlo, jhi
               do i = ilo, ihi
                  ok = correct(i,j,k,l)
                  if (oasym) ok = ok - correct(i,j,l,k)
                  if (abs(ok-buf(i,j,k,l)) .gt. 1d-8)
     $                 write(6,1) i,j,k,l,ok,buf(i,j,k,l)
 1                format(' wrong ', 4i5,2f20.10)
               end do
            end do
         end do
      end do
c
      end
      subroutine uccsdt_triples_get_list(list, 
     $        blk1, blk2, blk3, blk4, buf, lenbuf)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "mafdecls.fh"
      integer list, blk1, blk2, blk3, blk4, lenbuf
      double precision buf(lenbuf)
c     
c     blk1 is only used for list 12, 12 & 13
c
      integer offset,i,j,k,l,ptr,reclen,ind
c
c     Offset maps a 4-d array into a 1-d array.  It is used to look up
c     ptr(ind,i,j,k,l) where ptr is dimensioned ptr(2,dim1,dim2,dim3,dim4) 
c     with dim1-4 = listinfo(6-9,list)
c
      offset(ind,i,j,k,l,list) = 
     $     ind-1 + 
     $     2*(i-1 + 
     $     listinfo(6,list)*(j-1 + 
     $     listinfo(7,list)*(k-1 + 
     $     listinfo(8,list)*(l-1))))
c
      if (list.lt.1 .or. list.gt.23) call errquit
     $     ('uccsdt_triples_get_list: bad list', list, DISK_ERR)
      if (blk1.lt.1 .or. blk1.gt.listinfo(6,list)) call errquit
     $     ('uccsdt_triples_get_list: block 1 bad', list, DISK_ERR)
      if (blk2.lt.1 .or. blk2.gt.listinfo(7,list)) call errquit
     $     ('uccsdt_triples_get_list: block 2 bad', list, DISK_ERR)
      if (blk3.lt.1 .or. blk3.gt.listinfo(8,list)) call errquit
     $     ('uccsdt_triples_get_list: block 3 bad', list, DISK_ERR)
      if (blk4.lt.1 .or. blk4.gt.listinfo(9,list)) call errquit
     $     ('uccsdt_triples_get_list: block 4 bad', list, DISK_ERR)
c
      ptr = int_mb(listinfo(2,list) + 
     $     offset(1,blk1,blk2,blk3,blk4,list))
      reclen = int_mb(listinfo(2,list) +
     $     offset(2,blk1,blk2,blk3,blk4,list))
c
      if (reclen .gt. lenbuf) call errquit
     $     ('uccsdt_triples_get_list: lenbuf too small', reclen,
     &       DISK_ERR)
c
      call ga_get(listinfo(5,list),ptr,ptr+reclen-1,1,1,buf,1)
c
c$$$      if (list .eq. 1) then ! <ei||ab>
c$$$      else if (list .eq. 2) then ! <ek|ac>
c$$$      else if (list .eq. 3) then ! <ek|bc>
c$$$      else if (list .eq. 4) then ! <ei|ca>
c$$$      else if (list .eq. 5) then ! <ei|cb>
c$$$      else if (list .eq. 6) then ! <mc|ik>
c$$$      else if (list .eq. 7) then ! <ma||ij>
c$$$      else if (list .eq. 8) then ! <mb||ij>
c$$$      else if (list .eq. 9) then ! <ma|ki>
c$$$      else if (list .eq. 10) then ! <mb|ki>
c$$$      else if (list .eq. 11) then ! <ab||ij>
c$$$      else if (list .eq. 12) then ! <ac|ik>
c$$$      else if (list .eq. 13) then ! <bc|ik>
c$$$      else if (list .eq. 14) then ! t(e,c,i,k)
c$$$      else if (list .eq. 15) then ! t(e,a,i,j)
c$$$      else if (list .eq. 16) then ! t(e,b,i,j)
c$$$      else if (list .eq. 17) then ! t(e,a,k,i)
c$$$      else if (list .eq. 18) then ! t(e,b,k,i)
c$$$      else if (list .eq. 19) then ! t(m,i,a,b)
c$$$      else if (list .eq. 20) then ! t(m,k,a,c)
c$$$      else if (list .eq. 21) then ! t(m,k,b,c)
c$$$      else if (list .eq. 22) then ! t(m,i,c,a)
c$$$      else if (list .eq. 23) then ! t(m,i,c,b)
c$$$      else
c$$$         call errquit('uccsdt_triples_get_list: bad list', list)
c$$$      end if
c
      end
      subroutine uccsdt_triples_moints_init()
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "util.fh"
#include "eaf.fh"
c
      double precision availkw
      double precision util_scratch_dir_avail_for_me
      character*255 filename
c
      availkw = util_scratch_dir_avail_for_me()*.25d0
      write(6,*) ' AVAILKW ', availkw , 1.1d0*dble(nbf)**4/(4000d0*nir)
      if (availkw .gt. 1.1d0*dble(nbf)**4/(4000d0*nir)) then
         call util_file_name('.trpao', .true., .true., filename)
         if (eaf_open(filename, eaf_rw, aoints_file) .ne. 0) call
     $        errquit('uccsdt_triples: failed to open ao file',0,
     &       DISK_ERR)
         aoints_file_read = .false. 
         aoints_file_write = .true.
      else
         aoints_file = -1
         aoints_file_read = .false.
         aoints_file_write = .false.
      end if
c
      end
      subroutine uccsdt_triples_moints_tidy()
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "util.fh"
#include "eaf.fh"
c
      character*255 filename
c
      if (aoints_file .ne. -1) then
         if (eaf_close(aoints_file) .ne. 0) call errquit
     $        ('uccsdt_triples: failed to close aoints file',0,
     &       DISK_ERR)
         call util_file_name('.trpao', .true., .true., filename)
         if (eaf_delete(filename) .ne. 0) call errquit
     $        ('uccsdt_triples: failed to delete aoints file',0,
     &       DISK_ERR)
      end if
      aoints_file = -1
      aoints_file_read = .false.
      aoints_file_write = .false.
c
      end
      subroutine uccsdt_triples_moints()
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "cuccsdtps.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
c
c     Assumes that super- and sub-block info is already set
c     up in the common cuccsdtP.fh including the final global
c     arrays for block-sorted results.
c
      integer ushlo, ushhi, ush, vsh, ushblocksize, junk, list
      integer ublocklo, ublockhi, ublockdim, ublk
      parameter (ushblocksize=50) ! FOR DEBUG
      integer g_ai, g_bi, g_ck
      integer spini, spink, spine, spinm
      integer ilo, ihi, klo, khi, vlo, vhi, ulo, uhi, elo, ehi, mlo, mhi
      integer adim, bdim, cdim, idim, kdim, vdim, udim
      integer l_wivu, l_wkvu, l_work, l_buf, l_tmp, l_res
      integer k_wivu, k_wkvu, k_work, k_buf, k_tmp, k_res
      integer lenwork, lenbuf, dummy(3)
      integer me, next, nproc, count, maxdim, maxsh, shmaxbf
      integer shlist(2,nw_max_shells_small), nshblock, need
      double precision ga_need, ma_need
      logical ok
      integer nxtask
      external nxtask
      double precision q2, tol
      logical odoit, oskel
      parameter (tol = 1d-12, oskel = .true.)
      logical sym_shell_pair
      external sym_shell_pair
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c     This is written for the mixed-spin case in which 
c     spin(a)=spin(b)=spin(i)=spin(j) and spin(c)=spin(k) 
c     might be different.
c
c     The pure-spin case is handled by the same code except that
c     all of the integrals are antisymmetrized.
c
      spini = asuper_spin
      spink = csuper_spin
c
      ilo = nc(spini) + 1
      ihi = nc(spini) + no(spini)
      klo = nc(spink) + 1
      khi = nc(spink) + no(spink)
      adim = asuper(2) - asuper(1) + 1
      bdim = bsuper(2) - bsuper(1) + 1
      cdim = csuper(2) - csuper(1) + 1
      idim = no(spini)
      kdim = no(spink)
c
c     Determine the blocking for ush 
c     ... large blocks increase parallelism and minimize communication
c     ... need to leave enough memory for efficient integral evaluation
c
c     ublockdim = max. no. of basis functions in a shell block
c     maxsh  = max. no. of shells in a shell block
c     shmaxbf= max. no of functions in a shell
c
c     SYMMETRY IS USED TO REDUCE THE FLOP COUNT BUT IT SHOULD
c     ALSO BE USED TO REDUCE THE MEMORY REQUIRED.  NEED TO ADJUST
c     THE GLOBAL ARRAY SO THAT IT STORES SYMMETRY UNIQUE USH ONLY.
c
c     !!! IF BUFFERING ON FILE NEED TO USE SAME TASK DISTRIBUTION
c     !!! AFTER FIRST PASS WRITING THE FILE AND ALSO THE SAME
c     !!! BLOCKING FACTORS FOR THE AO INTEGRALS ... THIS IS NOT
c     !!! YET BEING DONE.
c
      if (aoints_file_write .and. nproc.gt.1) call errquit
     $     ('moints_triples: see comments above here!!', 0, DISK_ERR)
c
      if (.not. bas_nbf_cn_max(basis, shmaxbf)) call errquit
     $     ('uccsdt: bas_nbf_cn_max?', 0, BASIS_ERR)
      do maxdim = nbf, 1, -1
         call moints_shell_block(basis, maxdim, shlist, nshblock)         
         maxsh = 0
         do ublk = 1, nshblock
            maxsh = max(maxsh, shlist(2,ublk)-shlist(1,ublk)+1)
         end do
c
c     Leave room for a minimum AO shell blocking of 20 shells
c     ... could be more intelligent about this.
c
         call moints_aoblock_memory(1, 1, 20, 20, need)
         ma_need = need + shmaxbf*shmaxbf*
     $        (maxdim*maxdim + nbf*(idim+kdim))
         ga_need = dble(adim*idim + bdim*idim + cdim*kdim)*
     $        dble(maxdim*nbf)
         call util_ma_ga_test(ma_need, ga_need, ok)
         if (ok) goto 10
      end do
      call errquit('moints_shpair_one_index: insufficient memory',0,
     &       MEM_ERR)
 10   continue
c
      call moints_aoblock_rewind()
c
      do ublk = 1, nshblock
         ushlo = shlist(1,ublk)
         ushhi = shlist(2,ublk)
         if (.not. bas_cn2bfr(basis, ushlo, ublocklo, junk))
     $        call errquit('moints_shpair: basis?', 0, BASIS_ERR)
         if (.not. bas_cn2bfr(basis, ushhi, junk, ublockhi))
     $        call errquit('moints_shpair: basis?', 0, BASIS_ERR)
         ublockdim = ublockhi-ublocklo+1
c
c     First parallel step generates (u*|ai), (u*|bi) and (u*|ck)
c
         if (.not. ga_create(mt_dbl, ublockdim*nbf, adim*idim,
     $        '(u*|ai)', ublockdim*nbf, 0, g_ai)) call errquit
     $        ('moints_triples: failed allocaing (u*|ai)', 0, GA_ERR)
         if (.not. ga_create(mt_dbl, ublockdim*nbf, bdim*idim,
     $        '(u*|bi)', ublockdim*nbf, 0, g_bi)) call errquit
     $        ('moints_triples: failed allocaing (u*|bi)', 0, GA_ERR)
         if (.not. ga_create(mt_dbl, ublockdim*nbf, cdim*kdim,
     $        '(u*|ck)', ublockdim*nbf, 0, g_ck)) call errquit
     $        ('moints_triples: failed allocaing (u*|ck)', 0, GA_ERR)
         call ga_zero(g_ai)
         call ga_zero(g_bi)
         call ga_zero(g_ck)
c
c     Next two loops in parallel with dynamic load balancing
c
         count = 0
         next = nxtask(nproc, 1)
         do ush = ushlo, ushhi
            do vsh = 1, ush     ! NOTE RESTRICTION VSH <= USH
               odoit = schwarz_shell(ush,vsh)*schwarz_max().gt.tol
               if (odoit .and. oskel) 
     $              odoit = sym_shell_pair(basis, ush, vsh, q2)
               if (odoit) then
                  if (next .eq. count) then
                     if (.not. bas_cn2bfr(basis, ush, ulo, uhi))
     $                    call errquit('moints_shpair: basis?', 0,
     &       BASIS_ERR)
                     if (.not. bas_cn2bfr(basis, vsh, vlo, vhi))
     $                    call errquit('moints_shpair: basis?', 0,
     &       BASIS_ERR)
                     udim = uhi - ulo + 1
                     vdim = vhi - vlo + 1
c     
c     (wl|vu) -> (wi|vu), (wk|vu)
c     
                     if (.not. ma_push_get(mt_dbl, udim*vdim*nbf*idim,
     $                    '(wi|vu)', l_wivu, k_wivu)) call errquit
     $                    ('moints_triples: memory for (vu|wi)?', 
     $                    udim*vdim*nbf*idim, MA_ERR)
                     if (.not. ma_push_get(mt_dbl, udim*vdim*nbf*kdim,
     $                    '(wk|vu)', l_wkvu, k_wkvu)) call errquit
     $                    ('moints_triples: memory for (vu|wk)?', 
     $                    udim*vdim*nbf*kdim, MA_ERR)
c     
                     if (ouccsdtps) 
     $                    call pstat_on(ps_moints_shpair_one_index)
c     
                     call moints_shpair_one_index(
     $                    basis,
     $                    ush, ulo, uhi,
     $                    vsh, vlo, vhi,
     $                    nbf, nmo, 2,
     $                    oskel,
     $                    tol,
     $                    ilo, ihi, klo, khi,
     $                    dbl_mb(k_movecs(spini)),
     $                    dbl_mb(k_movecs(spink)),
     $                    dbl_mb(k_wivu), dbl_mb(k_wkvu))
c     
                     if (ouccsdtps) 
     $                    call pstat_off(ps_moints_shpair_one_index)
c     
c     Transform the second index and transpose on the fly
c     
                     lenwork = max(idim*adim,idim*bdim,kdim*cdim)
                     if (.not. ma_push_get(mt_dbl, lenwork, 'work',
     $                    l_work, k_work)) call errquit
     $                    ('moints_triples: memory for work', lenwork,
     &       MA_ERR)
c     
c     (wi|vu) -> (ai|vu) -> (vu|ai)
c     
c     
                     if (ouccsdtps) 
     $                    call pstat_on(ps_moints_second_index)
c     
                     call moints_triples_second_index(
     $                    ublocklo, ulo, uhi, vlo, vhi, 
     $                    ilo, ihi, asuper(1), asuper(2), 
     $                    nbf, nmo, dbl_mb(k_work), 
     $                    dbl_mb(k_movecs(spini)), 
     $                    dbl_mb(k_wivu), g_ai)
c
c     (wi|vu) -> (bi|vu) -> (vu|bi)
c     
                     call moints_triples_second_index(
     $                    ublocklo, ulo, uhi, vlo, vhi, 
     $                    ilo, ihi, bsuper(1), bsuper(2),
     $                    nbf, nmo, dbl_mb(k_work), 
     $                    dbl_mb(k_movecs(spini)), 
     $                    dbl_mb(k_wivu), g_bi)
c     
c     (wk|vu) -> (ck|vu) -> (vu|ck)
c     
                     call moints_triples_second_index(
     $                    ublocklo, ulo, uhi, vlo, vhi, 
     $                    klo, khi, csuper(1), csuper(2),
     $                    nbf, nmo, dbl_mb(k_work), 
     $                    dbl_mb(k_movecs(spink)), 
     $                    dbl_mb(k_wkvu), g_ck)
c
                     if (ouccsdtps) 
     $                    call pstat_off(ps_moints_second_index)
c     
                     if (.not. ma_pop_stack(l_work)) call errquit
     $                    ('moints_triples:popping work?', 0, MA_ERR)
                     if (.not. ma_pop_stack(l_wkvu)) call errquit
     $                    ('moints_triples:popping wkvu?', 0, MA_ERR)
                     if (.not. ma_pop_stack(l_wivu)) call errquit
     $                    ('moints_triples:popping wivu?', 0, MA_ERR)
c     
                     next = nxtask(nproc, 1)
                  end if
                  count = count + 1
               end if
            end do
         end do
c
c     Barrier
c
         next = nxtask(-nproc, 1)
         call ga_sync()
c
c     For all v, and u/a/b/c in their respective blocks we have
c     (vu|ai), (vu|bi) and (vu|ck) distributed by the last pair
c     of indices with all vu local.
c
c     For each of the 13 target integral cases transform the
c     vu indices and accumulate into the final locations.
c     Note that we transform the u index even though we might not
c     have all u ... this means that we do n_u_block times more
c     communication but we halve the memory requirement.
c     Hence the importance of using as large a u block as possible.
c
c     In the comments, the Dirac integrals denote the final storage
c     order, but the mulliken integrals (12|34) are written in the order
c     generated by the code.  The indices on the end of each line
c     indicate what order those mulliken indices should end up when
c     put in dirac order.
c
c     E.g., consider <ei|ab> = (ea|bi)
c     .              <14|23> = (12|34)
c
c     Thus, dirac index 1 = mulliken index 1
c     .     dirac index 2 = mulliken index 4
c     .     dirac index 3 = mulliken index 2
c     .     dirac index 4 = mulliken index 3
c
c     Hence, <ei|ab> = (ea|bi) = <14|23>
c
c     Also, in the pure spin case we need to antisymmetrize
c     all of the lists ... the integrals only needed in this
c     case are described in [...].
c
         lenbuf = nbf*
     $            max(ublockdim,adim,bdim,cdim,idim,kdim)
         if (.not. ma_push_get(mt_dbl, lenbuf, 'buf',
     $        l_buf, k_buf)) call errquit
     $        ('moints_triples: memory for buf', lenbuf, MA_ERR)
         if (.not. ma_push_get(mt_dbl, lenbuf, 'tmp',
     $        l_tmp, k_tmp)) call errquit
     $        ('moints_triples: memory for temp', lenbuf, MA_ERR)
         if (.not. ma_push_get(mt_dbl, lenbuf, 'res',
     $        l_res, k_res)) call errquit
     $        ('moints_triples: memory for temp', lenbuf, MA_ERR)
                  
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_moints_triples_final)
c     
c
c     1.  <ei||ab>  = (ea|bi) - (eb|ai) = <14|23> - <14|32>
c     
         list = 1
         spine = asuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_1a)
c     
         call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        nablock, ablock, asuper(1), asuper(2),
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|23>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_1a)
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_1b)
c     
         call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ai, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|32>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_1b)
c
c
c     2.  <ek|ac>   = (ea|ck) = <14|23>   [-<ek|ca> = -(ec|ai) = -<14|32>]
c
         list = 2
         spine = asuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_2a)
c     
         call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        nablock, ablock, asuper(1), asuper(2),
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|23>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_2a)
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_2b)
c     
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ai, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|32>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_2b)
c
c     3.  <ek|bc>   = (eb|ck) = <14|23>  [-<ek|cb> = -(ec|bk) = -<14|32>]
c
         list = 3
         spine = asuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_3a)
c     
         call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|23>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_3a)
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_3b)
c     
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|32>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_3b)
c
c     4.  <ei|ca>   = (ec|ai) = <14|23>  [-<ei|ac> = -(ea|ci) = -<14|32>]
c
         list = 4
         spine = csuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_4a)
c     
         call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ai, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(csuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|23>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_4a)
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_4b)
c     
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        nablock, ablock, asuper(1), asuper(2),
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(csuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|32>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_4b)
c
c     5.  <ei|cb>   = (ec|bi) = <14|23>  [-<ei|bc> = -(eb|ci) = -<14|32>]
c
         list = 5
         spine = csuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_5a)
c     
         call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(csuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|23>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_5a)
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_5b)
c     
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .false., v_sym(1,0,spine), 
     $        1, dummy, elo, ehi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spine)), dbl_mb(k_movecs(csuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<14|32>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_5b)
c
c     6.  <mc|ik>   = (mi|ck) = <13|24>  [-<mc|ki> = -(mk|ci) = -<13|42>]
c
         list = 6
         spinm = spini
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_6)
c     
         call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spini)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spini)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_6)
c
c     7.  <ma||ij>  = (mi|aj) - (mj|ai) = <13|24> - <13|42>
c
         list = 7
         spinm = spini
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_7)
c     
         call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ai, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spini)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ai, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spini)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_7)
c
c     8.  <mb||ij>  = (mi|bj) - (mj|bi) = <13|24> - <13|42>
c
         list = 8
         spinm = spini
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_8)
c     
         call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spini)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spini)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_8)
c
c     9.  <ma|ki>   = (mk|ai) = <13|24>  [-<ma|ik> = -(mi|ak) = -<13|42>]
c
         list = 9
         spinm = spink
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_9)
c     
         call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ai, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spink)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_ai, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spink)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_9)
c
c     10. <mb|ki>   = (mk|bi) = <13|24>  [-<mb|ik> = -(mi|bk) = -<13|42>]
c
         list = 10
         spinm = spink
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_10)
c     
         call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spink)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .false., o_sym(1,0,spinm), 
     $        1, dummy, mlo, mhi,
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(spinm)), dbl_mb(k_movecs(spink)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_10)
c
c     11. <ab||ij>  = (ai|bj) - (aj|bi) = <13|24> - <13|42>
c
         list = 11
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_11)
c     
         call moints_triples_final(
     $        .true., dummy,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(asuper_spin)), 
     $        dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         call moints_triples_final(
     $        .true., dummy,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        g_bi, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(asuper_spin)), 
     $        dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_11)
c
c     12. <ac|ik>   = (ai|ck) = <13|24>  [-<ac|ki> = -(ak|ci) = -<13|42>]
c
         list = 12
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_12)
c     
         call moints_triples_final(
     $        .true., dummy,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(asuper_spin)), 
     $        dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .true., dummy,
     $        nablock, ablock, asuper(1), asuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(asuper_spin)), 
     $        dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_12)
c
c     13. <bc|ik>   = (bi|ck) = <13|24>  [-<bc|ki> = -(bk|ci) = -<13|42>]
c
         list = 13
c     
         if (ouccsdtps) 
     $        call pstat_on(ps_mo_fin_13)
c     
         call moints_triples_final(
     $        .true., dummy,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(asuper_spin)), 
     $        dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|24>')
         if (asuper_spin .eq. csuper_spin) call moints_triples_final(
     $        .true., dummy,
     $        nbblock, bblock, bsuper(1), bsuper(2),
     $        noblock(spini), oblock(1,1,spini), ilo, ihi,
     $        ncblock, cblock, csuper(1), csuper(2),
     $        noblock(spink), oblock(1,1,spink), klo, khi,
     $        g_ck, dbl_mb(k_buf), dbl_mb(k_tmp), dbl_mb(k_res),
     $        dbl_mb(k_movecs(asuper_spin)), 
     $        dbl_mb(k_movecs(asuper_spin)),
     $        ublocklo, ublockhi, nbf, 
     $        int_mb(listinfo(2,list)), 
     $        listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $        listinfo(5,list), '<13|42>')
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_mo_fin_13)
c     
         if (ouccsdtps) 
     $        call pstat_off(ps_moints_triples_final)
c
         if (.not. ma_pop_stack(l_res))
     $        call errquit('moints_triple: popping res?', 0, MA_ERR)
         if (.not. ma_pop_stack(l_tmp))
     $        call errquit('moints_triple: popping tmp?', 0, MA_ERR)
         if (.not. ma_pop_stack(l_buf))
     $        call errquit('moints_triple: popping buf?', 0, MA_ERR)
c
         if (.not. ga_destroy(g_ck)) call errquit
     $        ('moints_triples: failed to destroy g_ck',0, GA_ERR)
         if (.not. ga_destroy(g_bi)) call errquit
     $        ('moints_triples: failed to destroy g_bi',0, GA_ERR)
         if (.not. ga_destroy(g_ai)) call errquit
     $        ('moints_triples: failed to destroy g_ai',0, GA_ERR)
      end do
c
c      if (ocheck) then
c         do list = 1, 13
c            call uccsdt_check_list(list)
c         end do
c      endif
c
*      call ma_summarize_allocated_blocks()
c
      if (aoints_file_write) then
         aoints_file_write = .false.
         aoints_file_read  = .true.
      end if
c
      end
      subroutine moints_triples_final(
     $     oblocke, e_sym, 
     $     neblock, eblock, elo, ehi,
     $     nablock, ablock, alo, ahi,
     $     ncblock, cblock, clo, chi,
     $     nkblock, kblock, klo, khi,
     $     g_ck, buf, temp, result,
     $     evecs, avecs,
     $     ublocklo, ublockhi, nbf, pointers, dim1, dim2, dim3,
     $     g_list,
     $     order)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "cuccsdtps.fh"
      logical oblocke           ! [input] If true sub-block e index.
      integer e_sym(2,0:7)      ! [input] Range of e-orbitals in each irrep
      integer neblock, eblock(3,*) ![input] Sub-blocking info for e-orbitals
      integer elo, ehi          ! [input] Total range of e-orbitals
      integer nablock, ablock(3,*) ![input] Sub-blocking info for a-orbitals
      integer alo, ahi          ! [input] Total range of a-orbitals
      integer ncblock, cblock(3,*) ![input] Sub-blocking info for c-orbitals
      integer clo, chi          ! [input] Total range of c-orbitals
      integer nkblock, kblock(3,*) ![input] Sub-blocking info for k-orbitals
      integer klo, khi          ! [input] Total range of k-orbitals
      integer g_ck              ! [input] Array holding (vu|ck)
      integer nbf               ! [input] Total no. of basis functions
      double precision buf(*)   ! [scratch] ublockdim*nbf
      double precision temp(*)  ! [scratch] ublockdim*nbf
      double precision result(elo:ehi,alo:ahi) ! [scratch] edim*adim
      double precision evecs(nbf,*) ! [input] MO vectors for e-orbitals
      double precision avecs(nbf,*) ! [input] MO vectors for a-orbitals
      integer ublocklo, ublockhi ! [input] Range of functions in ush block
      integer dim1, dim2, dim3
      integer pointers(2,dim1,dim2,dim3,*) ! [input] Offset into g_list
      integer g_list            ! [input] Array to accumulate into
      character*(*) order
c
c     In g_ck we have (vu|ck) with all v and u in block local, and 
c     distributed by ck.  For the ck that are local, transform
c     (vu) -> (ea) and store the result according to order
c     which specifies the order of indices on the output Dirac integral
c     if the Mulliken indices are (ea|ck) = (12|34) (1=e,2=a,3=c,4=k)
c
c     Note that a phase is also included
c
c     1. <14|23> -> +<ek|ac>
c     2. <14|32> -> -<ek|ca>
c     3. <13|24> -> +<ec|ak>
c     4. <13|42> -> -<ec|ka>
c
c     There are two options for how the e index is handled
c
c     if (oblocke) then
c     .  The e index is blocked as are the other indices.
c     .  This applies to lists 11, 12, 13.  
c     .  neblock and eblock() are used.
c     .  e_sym() is not used.
c     else
c     .  The e index is not blocked and all e that are non-zero
c     .  by symmetry are transformed.
c     .  e_sym() is used.
c     .  neblock should be set to 1
c     .  eblock() is not used.
c     endif
c     
c     Loop thru all c & k in sub-block order BUT only process those
c     that are local to this processor
c
      integer kblk, k, cblk, c, ablk, a, ck, eblk, ebase
      integer edim, adim, cdim, ublockdim, nc, na, nk, ne
      integer vulo, vuhi, cklo, ckhi, me, syme
      integer ptr
      integer case
      double precision phase
c
      if ((.not. oblocke) .and. (neblock.ne.1))
     $     call errquit('moints_triples_final: neblock?',neblock,
     &       DISK_ERR)
c
c     Convert order to an integer to avoid character comparisons
c     in the inner loop
c
      if (order .eq. '<14|23>') then
         case = 1
         phase = 1d0
      else if (order .eq. '<14|32>') then
         case = 2
         phase = -1d0
      else if (order .eq. '<13|24>') then
         case = 3
         phase = 1d0
      else if (order .eq. '<13|42>') then
         case = 4
         phase = -1d0
      else
         case = 0               ! Just to keep "smart" compiler quiet
         call errquit('moints_triples_final: bad order', 0, DISK_ERR)
      end if
c
      ublockdim = (ublockhi-ublocklo+1)
      edim = (ehi-elo+1)        ! Total no. of e-orbitals
      cdim = (chi-clo+1)        ! Total no. of c-orbitals
      adim = (ahi-alo+1)
      me = ga_nodeid()
c     
      if (ouccsdtps) 
     $     call pstat_on(ps_final_ga_dist)
c     
      call ga_distribution(g_ck, me, vulo, vuhi, cklo, ckhi)
      if (vulo.ne.1 .and. vuhi.ne.ublockdim*nbf)
     $     call errquit('moints_triples_final: distribution?', 0,
     &       DISK_ERR)
c     
      if (ouccsdtps) 
     $     call pstat_off(ps_final_ga_dist)
c     
c      
      do kblk = 1, nkblock
         nk = kblock(2,kblk) - kblock(1,kblk) + 1
         do k = kblock(1,kblk), kblock(2,kblk)
            do cblk = 1, ncblock
               nc = cblock(2,cblk) - cblock(1,cblk) + 1
               do c = cblock(1,cblk),cblock(2,cblk)
                  ck = 1 + c-clo + (k-klo)*cdim
                  if (ck.ge.cklo .and. ck.le.ckhi) then
c     
                    if (ouccsdtps) 
     $                   call pstat_on(ps_final_ga_get)
c                   
                     call ga_get(g_ck, 1, ublockdim*nbf, ck, ck,
     $                    buf, 1)
c     
                    if (ouccsdtps) 
     $                   call pstat_off(ps_final_ga_get)
c
c     C(v,e)*buf(v,u) -> temp(e,u)
c     
                    if (ouccsdtps) 
     $                   call pstat_on(ps_final_dgemm1)
c                   
                     call dgemm('t','n',edim,ublockdim,nbf,
     $                    1.0d0, evecs(1,elo), nbf, buf, nbf,
     $                    0.0d0, temp, edim)
c     
                    if (ouccsdtps) 
     $                   call pstat_off(ps_final_dgemm1)
c
c     phase*temp(e,u)*C(u,a) -> result(e,a)
c     
                    if (ouccsdtps) 
     $                   call pstat_on(ps_final_dgemm2)
c                   
                     call dgemm('n','n',edim,adim,ublockdim,
     $                    phase, temp, edim, avecs(ublocklo,alo), nbf,
     $                    0.0d0, result, edim)
c     
                    if (ouccsdtps) 
     $                   call pstat_off(ps_final_dgemm2)
c
c     We just did (vu) -> (eu) -> (ea), but since we are restricting
c     u>=v we also need to do (vu) -> (au) -> (ae).
c
                    if (ouccsdtps) 
     $                   call pstat_on(ps_final_u_trans)
c                   
                     call util_transpose(buf, nbf, temp, ublockdim,
     $                    nbf, ublockdim)
c
                    if (ouccsdtps) 
     $                   call pstat_off(ps_final_u_trans)
c                   
c     temp(u,v)*C(v,a) -> buf(u,a)
c
                    if (ouccsdtps) 
     $                   call pstat_on(ps_final_dgemm3)
c                   
                     call dgemm('n', 'n', ublockdim, adim, nbf,
     $                    1.0d0, temp, ublockdim, avecs(1,alo), nbf,
     $                    0.0d0, buf, ublockdim)
c     
                    if (ouccsdtps) 
     $                   call pstat_off(ps_final_dgemm3)
c
c     phase*C(u,e)*buf(u,a) -> result(e,a)
c
                    if (ouccsdtps) 
     $                   call pstat_on(ps_final_dgemm4)
c                   
                     call dgemm('t','n',edim,adim,ublockdim,
     $                    phase, evecs(ublocklo,elo),nbf,buf,ublockdim,
     $                    1.0d0, result, edim)
c     
                    if (ouccsdtps) 
     $                   call pstat_off(ps_final_dgemm4)
c
c     Now have (ea|ck) in result ... note that if symmetry was used
c     that only elements that are non-zero by symmetry will be correct.
c     Accumulate the elements for each block into the appropriate place
c     inside a sub-block of integrals.
c     
                     do ablk = 1, nablock
                        na = ablock(2,ablk) - ablock(1,ablk) + 1
                        syme = ieor(ablock(3,ablk),
     $                       ieor(cblock(3,cblk),kblock(3,kblk)))
                        do a = ablock(1,ablk),ablock(2,ablk)
                           do eblk = 1, neblock
                              if (oblocke) then
                                 if (syme.eq.eblock(3,eblk)) then
                                    ebase = eblock(1,eblk)
                                    ne = eblock(2,eblk)-eblock(1,eblk)+1
                                 else
                                    ne = 0
                                 end if
                              else
                                 ebase = e_sym(1,syme)
                                 ne = e_sym(2,syme) - e_sym(1,syme) + 1
                              end if
                              if (ne .gt. 0) then
                                 if (case .eq. 1) then ! +<ek|ac>
                                    ptr=pointers(1,eblk,kblk,ablk,cblk)+ 
     $                                   ne*(k-kblock(1,kblk) + 
     $                                   nk*(a-ablock(1,ablk) + 
     $                                   na*(c-cblock(1,cblk))))
                                 else if (case .eq. 2) then ! -<ek|ca>
                                    ptr=pointers(1,eblk,kblk,cblk,ablk)+ 
     $                                   ne*(k-kblock(1,kblk) + 
     $                                   nk*(c-cblock(1,cblk) + 
     $                                   nc*(a-ablock(1,ablk))))
                                 else if (case .eq. 3) then ! +<ec|ak>
                                    ptr=pointers(1,eblk,cblk,ablk,kblk)+ 
     $                                   ne*(c-cblock(1,cblk) + 
     $                                   nc*(a-ablock(1,ablk) + 
     $                                   na*(k-kblock(1,kblk))))
                                 else if (case .eq. 4) then ! +<ec|ka>
                                    ptr=pointers(1,eblk,cblk,kblk,ablk)+ 
     $                                   ne*(c-cblock(1,cblk) + 
     $                                   nc*(k-kblock(1,kblk) + 
     $                                   nk*(a-ablock(1,ablk))))
                                 endif
c     
                    if (ouccsdtps) 
     $                   call pstat_on(ps_final_ga_acc)
c                   
                                 call ga_acc(g_list, ptr, ptr+ne-1,
     $                                1, 1, result(ebase,a), 
     $                                1, 1.0d0)
c     
                    if (ouccsdtps) 
     $                   call pstat_off(ps_final_ga_acc)
c                   
                              end if
                           end do
                        end do
                     end do
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine moints_triples_second_index(
     $     ublocklo, 
     $     ulo, uhi,
     $     vlo, vhi,
     $     ilo, ihi, 
     $     alo, ahi,
     $     nbf, nmo,
     $     work,
     $     c, g, g_vuai)
      implicit none
#include "mafdecls.fh"
      integer ublocklo          ! [input]
      integer ulo, uhi          ! [input]
      integer vlo, vhi          ! [input]
      integer ilo, ihi          ! [input]
      integer alo, ahi          ! [input]
      integer nbf, nmo          ! [input]
      double precision work(alo:ahi,ilo:ihi) ! [scratch]
      double precision c(nbf,*), g(nbf,ilo:ihi,vlo:vhi,ulo:uhi) ! [input]
      integer g_vuai            ! [input]
c
c     In g(w,i,v,u) we have (wi|vu) (a one-index transformed integral).
c     Transform w->a to get (ai|vu) and then store the result
c     transposed into the global array as GA(v,u,a,i)
c
      integer u, v, uv, adim, idim
c
      adim = ahi - alo + 1
      idim = ihi - ilo + 1
      do u = ulo, uhi
         do v = vlo, vhi
            call dgemm('t','n',adim,idim,nbf,
     $           1.0d0, c(1,alo), nbf, g(1,ilo,v,u), nbf,
     $           0.0d0, work, adim)
            uv = v + (u-ublocklo)*nbf
            call ga_put(g_vuai, uv, uv, 1, idim*adim, work, 1)
         end do
      end do
c
      end
      
      subroutine moints_shpair_one_index(
     $     basis,
     $     ush, ulo, uhi,
     $     vsh, vlo, vhi,
     $     nbf, nmo, nset,
     $     oskel,
     $     tol,
     $     ilo1, ihi1, ilo2, ihi2, 
     $     c1, c2,
     $     g1, g2)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      integer basis             ! [input]
      integer ush, ulo, uhi     ! [input] shell & bf range
      integer vsh, vlo, vhi     ! [input] shell & bf range
      integer nbf               ! [input] no. of basis functions
      integer nmo               ! [input] no. of MOs (& dim. of c)
      integer nset              ! [input] no. of sets of MOs
      logical oskel             ! [input] if true use petite list & q4
      double precision tol      ! [input] Screening for integrals
      integer ilo1, ihi1        ! [input] range of MOs to transform
      integer ilo2, ihi2        ! [input] range of MOs to transform
      double precision c1(nbf,nmo) ! [input] MO coeffs
      double precision c2(nbf,nmo) ! [input] MO coeffs
      double precision g1(nbf,ilo1:ihi1,vlo:vhi,ulo:uhi)
      double precision g2(nbf,ilo2:ihi2,vlo:vhi,ulo:uhi)
c     .                         ! [output] g(w,i,v,u,iset) = (wi|vu)
c     (wi|vu) = sum(l=1,nbf) (wl|vu) Cli
c
c     For given AO shells ush & vsh return in g all integrals with one 
c     index transformed into the MO basis for multiple sets of 
c     MO coefficients ... nset=1 or 2 ... extend to handle more.
c
c     If (oskel) then restrict the generated AO integrals to be
c     the symmetry unique set (with no permutation restriction on the
c     indices) and weight with the appropriate q4 factor.  The resulting
c     integrals may be transformed into a symmetry adapted basis
c     (note that the integrals that should be zero by symmetry must
c     be ignored).
c
c     Integrals and schwarz must already be initialized
c
c     Could add the use of symmetry but in initial application the
c     symmetry is being handled by using unique ush.
c
      integer idim1, idim2, udim, vdim
      integer nsh
      integer wshlo, wshhi, lshlo, lshhi
      integer wlo, whi, llo, lhi, junk, l_work, k_work, lenwork
      integer avail_ma, need, wblk, lblk
c
      integer shlist(2,nw_max_shells_small), nshblock, maxdim, maxsh
c
      if (nset.lt.1 .or. nset.gt.2) call errquit
     $     ('moints_shpair_one_index: bad nset', nset, DISK_ERR)
      udim = uhi - ulo + 1
      vdim = vhi - vlo + 1
      idim1= ihi1 - ilo1 + 1
      call dfill(nbf*idim1*udim*vdim, 0.0d0, g1, 1)
      if (nset .eq . 2) then
         idim2 = ihi2 - ilo2 + 1
         call dfill(nbf*idim2*udim*vdim, 0.0d0, g2, 1)
      end if
      if (schwarz_shell(ush,vsh)*schwarz_max() .lt. tol) return
c
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('moints_shpair_one_index: basis?', basis,
     &       BASIS_ERR)
      if (.not. bas_numbf(basis,nbf))
     $     call errquit('moints_shpair_one_index: basis?', basis,
     &       BASIS_ERR)
c
c     Figure out how big we can make the shell blocks
c
c     maxdim is the maximum no. of basis functions in each shell block
c     maxsh  is the maximum no. of shells in each block
c
      avail_ma = ma_inquire_avail(mt_dbl)
      avail_ma = nint(dble(avail_ma) * 0.95d0)
      do maxdim = nbf, 1, -1
         call moints_shell_block(basis, maxdim, shlist, nshblock)         
         maxsh = 0
         do wblk = 1, nshblock
            maxsh = max(maxsh, shlist(2,wblk)-shlist(1,wblk)+1)
         end do
         call moints_aoblock_memory(1, 1, maxsh, maxsh, need)
         need = need + udim*vdim*maxdim*maxdim
         if (avail_ma .gt. need) goto 10
      end do
      call errquit('moints_shpair_one_index: insufficient memory',0,
     &       MEM_ERR)
 10   continue
c
      lenwork = udim*vdim*maxdim*maxdim
      if (.not. ma_push_get(mt_dbl, lenwork, 'work', 
     $     l_work, k_work)) 
     $     call errquit('moints_shpair: memory estimate?',lenwork,
     &       MA_ERR)
c
c     Loop over blocks of shells wsh & lsh in order to control memory
c     ... smaller blocks = less memory but less efficient integral
c     evaluation.  
c
c     Consider switching to texas 93 since this will give larger
c     blocks than 95.
c
      do wblk = 1, nshblock
         wshlo = shlist(1,wblk)
         wshhi = shlist(2,wblk)
         if (.not. bas_cn2bfr(basis, wshlo, wlo, junk))
     $        call errquit('moints_shpair: basis?', 0, BASIS_ERR)
         if (.not. bas_cn2bfr(basis, wshhi, junk, whi))
     $        call errquit('moints_shpair: basis?', 0, BASIS_ERR)
         do lblk = 1, wblk
            lshlo = shlist(1,lblk)
            lshhi = shlist(2,lblk)
            if (.not. bas_cn2bfr(basis, lshlo, llo, junk))
     $           call errquit('moints_shpair: basis?', 0, BASIS_ERR)
            if (.not. bas_cn2bfr(basis, lshhi, junk, lhi))
     $           call errquit('moints_shpair: basis?', 0, BASIS_ERR)
            call moints_shpair_block(
     $           basis,
     $           lshlo, lshhi, llo, lhi,
     $           wshlo, wshhi, wlo, whi,
     $           vsh, vlo, vhi,
     $           ush, ulo, uhi,
     $           nbf, nmo, nset,
     $           oskel,
     $           tol,
     $           ilo1, ihi1, ilo2, ihi2, 
     $           c1, c2,
     $           g1, g2,
     $           dbl_mb(k_work))
         end do
      end do
      if (.not. ma_pop_stack(l_work)) 
     $     call errquit('moints_shpair: ma corrupt?',0, MA_ERR)
c
      end
      subroutine moints_shpair_block(
     $           basis,
     $           lshlo, lshhi, llo, lhi,
     $           wshlo, wshhi, wlo, whi,
     $           vsh, vlo, vhi,
     $           ush, ulo, uhi,
     $           nbf, nmo, nset,
     $           oskel,
     $           tol,
     $           ilo1, ihi1, ilo2, ihi2, 
     $           c1, c2,
     $           g1, g2,
     $           ints)
      implicit none
      integer basis             ! [input]
      integer ush, ulo, uhi     ! [input] shell & bf range
      integer vsh, vlo, vhi     ! [input] shell & bf range
      integer wshlo, wshhi      ! [input] shell range for w block
      integer wlo, whi          ! [input] corresponding bf range
      integer lshlo, lshhi      ! [input] shell range for l block
      integer llo, lhi          ! [input] corresponding bf range
      integer nbf               ! [input] no. of basis functions
      integer nmo               ! [input] total no. of MOs
      integer nset              ! [input] no. of sets of MOs
      logical oskel             ! [input] use symmetry
      double precision tol    ! [input] Screening for integrals
      integer ilo1, ihi1        ! [input] range of MOs to transform
      integer ilo2, ihi2        ! [input] range of MOs to transform
      double precision c1(nbf,nmo) ! [input] MO coeffs
      double precision c2(nbf,nmo) ! [input] MO coeffs
      double precision g1(nbf,ilo1:ihi1,vlo:vhi,ulo:uhi)
      double precision g2(nbf,ilo2:ihi2,vlo:vhi,ulo:uhi)
c     .                         ! [output] g(w,i,v,u,iset) = (wi|vu)
      double precision ints(llo:lhi,wlo:whi,vlo:vhi,ulo:uhi)
c
      integer u, v, idim1, idim2, wdim, ldim
      integer aodim
      logical moints_aoblock_read
c
      idim1 = ihi1 - ilo1 + 1
      if (nset .eq. 2) idim2 = ihi2 - ilo2 + 1
      wdim = whi - wlo + 1
      ldim = lhi - llo + 1
      aodim = ldim*wdim*(vhi-vlo+1)*(uhi-ulo+1)
c      
c     Form in ints(l,w,v,u) the integrals (lw|vu) with l<=w
c
      if (.not. moints_aoblock_read(ush,vsh,aodim,ints)) then
         call moints_aoblock(
     $        basis, oskel,
     $        lshlo, lshhi, llo, lhi,
     $        wshlo, wshhi, wlo, whi,
     $        vsh, vsh, vlo, vhi,
     $        ush, ush, ulo, uhi,
     $        tol,
     $        ints)
         call moints_aoblock_write(ush,vsh,aodim,ints)
      endif
c
c     Now have all of the AO integrals ... transform as necessary.
c     Since we are restricting w>=l we have to do l->i and w->i.
c
      do u = ulo, uhi
         do v = vlo, vhi
            call dgemm('t','n',wdim,idim1,ldim, ! (lw|vu)*Cli->(wi|vu)
     $           1.0d0, ints(llo,wlo,v,u), ldim, c1(llo,ilo1), nbf,
     $           1.0d0, g1(wlo, ilo1, v, u), nbf)
            call dgemm('n','n',ldim,idim1,wdim, ! (lw|vu)*Cwi->(li|vu)
     $           1.0d0, ints(llo,wlo,v,u), ldim, c1(wlo,ilo1), nbf,
     $           1.0d0, g1(llo, ilo1, v, u), nbf)
            if (nset .eq. 2) then
               call dgemm('t','n',wdim,idim2,ldim, ! (lw|vu)*Cli->(wi|vu)
     $              1.0d0, ints(llo,wlo,v,u), ldim, c2(llo,ilo2), nbf,
     $              1.0d0, g2(wlo, ilo2, v, u), nbf)
               call dgemm('n','n',ldim,idim2,wdim, ! (lw|vu)*Cwi->(li|vu)
     $              1.0d0, ints(llo,wlo,v,u), ldim, c2(wlo,ilo2), nbf,
     $              1.0d0, g2(llo, ilo2, v, u), nbf)
            end if
         end do
      end do
c
      end
      subroutine moints_aoblock_memory(nish, njsh, nksh, nlsh, memory)
      implicit none
#include "mafdecls.fh"
      integer nish, njsh, nksh, nlsh ! [input]
      integer memory            ! [output]
c
c     Return the amount of memory in doubles that will be internally 
c     allocated by moints_aoblock for the given no. of shells.
c
      integer maxq, leneri, lenscr
c
      maxq = nish * njsh * nksh * nlsh
      call intb_mem_2e4c(leneri, lenscr)
      leneri = max(leneri,100*(6**4)) ! Room for 100 (dd|dd)
c
      memory = maxq + leneri + lenscr + 
     $     ma_sizeof(mt_int, 4*maxq + 4*leneri, mt_dbl)
c
      end
      subroutine moints_aoblock(
     $     basis, oskel,
     $     ishlo, ishhi, ilo, ihi,
     $     jshlo, jshhi, jlo, jhi,
     $     kshlo, kshhi, klo, khi,
     $     lshlo, lshhi, llo, lhi,
     $     tol,
     $     ints)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "schwarz.fh"
#include "nwc_const.fh"
#include "bas.fh"
c
c     For the shells in the given range compute the AO integrals
c     returning them in the array ints(i,j,k,l) = (ij|kl) 
c     with ish<=jsh and ksh<=lsh and appropriate weighting 
c     on the diagonals (1/2 if ish=jsh, 1/2 if ksh=lsh).
c
      integer
     $     basis,
     $     ishlo, ishhi, ilo, ihi,
     $     jshlo, jshhi, jlo, jhi,
     $     kshlo, kshhi, klo, khi,
     $     lshlo, lshhi, llo, lhi
      logical oskel             ! If true use petite list symmetry
      double precision tol
      double precision ints(ilo:ihi,jlo:jhi,klo:khi,llo:lhi)
c
      integer ish, jsh, ksh, lsh
      integer i, j, k, l
      integer l_q4, l_scr, l_eri, l_lab, l_q
      integer k_q4, k_scr, k_eri, k_lab, k_q
      integer k_iq, k_jq, k_kq, k_lq
      integer k_ilab, k_jlab, k_klab, k_llab
      integer leneri, lenscr
      integer maxq, idim, jdim, kdim, ldim
      logical more, intb_init4c, intb_2e4c
      integer nq, neri, ieri
      logical odoit
      double precision block_eff
      double precision scalekl, scale, smax, skl, q4
c
      integer atom(nw_max_shells_small), geom, nsh
      logical sym_atom_gen_quartet
      external sym_atom_gen_quartet
      logical sym_shell_pair
      external sym_shell_pair
c
      idim = ihi - ilo + 1
      jdim = jhi - jlo + 1
      kdim = khi - klo + 1
      ldim = lhi - llo + 1
      call dfill(ldim*kdim*jdim*idim, 0.0d0, ints, 1)
c
      if (oskel) then
         if (.not. bas_geom(basis, geom)) call errquit
     $     ('moints_aoblock: bas_geom?', 0, BASIS_ERR)
         if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('moints_aoblock: bas_numcont?', 0, BASIS_ERR)
         do ish = 1, nsh
            if (.not. bas_cn2ce(basis, ish, atom(ish))) 
     $           call errquit('moints_aoblock: bas_cn2ce', ish,
     &       BASIS_ERR)
         end do
      end if
c
c     Naively try to compute all of the quartests at once!
c
      maxq = (kshhi-kshlo+1)*(lshhi-lshlo+1)*
     $     (ishhi-ishlo+1)*(jshhi-jshlo+1)
c
      call intb_mem_2e4c(leneri, lenscr)
      leneri = max(leneri,100*(6**4)) ! Room for 100 (dd|dd)
c
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('moints_aoblock: ma?',maxq, MA_ERR)
      if (.not. ma_push_get(mt_int, 4*maxq, 'quartets', l_q, k_q))
     $     call errquit('moints_aoblock: ma?', 4*maxq, MA_ERR)
      if (.not. ma_push_get(mt_int, 4*leneri, 'labels', l_lab, k_lab))
     $     call errquit('moints_aoblock: ma?', 4*leneri, MA_ERR)
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri))
     $     call errquit('moints_aoblock: ma?', leneri, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr))
     $     call errquit('moints_aoblock: ma?', lenscr, MA_ERR)
      k_iq = k_q          
      k_jq = k_q + maxq   
      k_kq = k_q + maxq*2 
      k_lq = k_q + maxq*3 
      k_ilab = k_lab       
      k_jlab = k_lab + leneri  
      k_klab = k_lab + leneri*2
      k_llab = k_lab + leneri*3
c
      nq = 0
      q4 = 1.0d0                ! If not using symmetry
      smax = schwarz_max()
      do lsh = lshlo, lshhi
         do ksh = kshlo, min(lsh,kshhi)
            skl = schwarz_shell(lsh,ksh)
            if (skl*smax .gt. tol) then
               scalekl = 1.0d0
               if (ksh .eq. lsh) scalekl = 0.5d0
               do jsh = jshlo, jshhi
                  do ish = ishlo, min(jsh,ishhi)
                     odoit = schwarz_shell(jsh,ish)*skl .gt. tol
c
c     The call to sym_atom_gen_quartet must have lsh, ksh first
c     if these are the outermost loops and they are being tested 
c     with sym_shell_pair(lsh,ksh)
c     
                     if (odoit .and. oskel) 
     $                    odoit = sym_atom_gen_quartet(geom, 
     $                    atom(lsh), atom(ksh), atom(jsh), atom(ish),q4)
                     if (odoit) then
                        int_mb(k_iq + nq) = ish
                        int_mb(k_jq + nq) = jsh
                        int_mb(k_kq + nq) = ksh
                        int_mb(k_lq + nq) = lsh
                        scale = scalekl
                        if (ish .eq. jsh) scale = 0.5d0*scale
                        dbl_mb(k_q4 + nq) = scale*q4
                        nq = nq + 1
                     endif
                  end do
               end do
            end if
         end do
      end do
      if (nq .gt. 0) then
         if (.not. intb_init4c(
     $        basis, int_mb(k_iq), int_mb(k_jq), 
     $        basis, int_mb(k_kq), int_mb(k_lq), 
     $        nq, dbl_mb(k_q4), .true., 
     $        lenscr, dbl_mb(k_scr), leneri, block_eff))
     $        call errquit('intb_init?',nq, UNKNOWN_ERR)
c
 10      more = intb_2e4c(
     $        basis, int_mb(k_iq), int_mb(k_jq), 
     $        basis, int_mb(k_kq), int_mb(k_lq), 
     $        nq, dbl_mb(k_q4), .true., tol, .false.,
     $        int_mb(k_ilab), int_mb(k_jlab),
     $        int_mb(k_klab), int_mb(k_llab),
     $        dbl_mb(k_eri), leneri, neri, 
     $        lenscr, dbl_mb(k_scr))
c     
         do ieri = 0, neri-1
            i = int_mb(k_ilab + ieri)
            j = int_mb(k_jlab + ieri)
            k = int_mb(k_klab + ieri)
            l = int_mb(k_llab + ieri)
            ints(i,j,k,l) = dbl_mb(k_eri + ieri)
         end do
c
      if (more) goto 10
c
      end if
c
c
      if (.not. ma_chop_stack(l_q4)) call errquit
     $     ('moints_aoblock: ma corrupted?',0, MA_ERR)
c
      end
      subroutine uccsdt_block_range(lo, hi, maxb, nblock, blocks, 
     $     ld, maxnblock)
      implicit none
#include "errquit.fh"
      integer lo, hi, maxb, maxnblock, ld ! [input]
      integer nblock, blocks(ld,maxnblock) ! [input]
c
c     Given a range [lo:hi] divide it into blocks so that each
c     block is not bigger than maxb but, if possible, we don't
c     have small trailing blocks.
c
c     nblock .. returns the no. of blocks
c     blocks(1,i) ... first index in i'th block
c     blocks(2,i) ... last  index in i'th block
c     
      integer n, b, i
c
      if (hi .lt. lo) then 
         nblock = 0
         return
      endif
c
      n = (hi - lo + 1)
      nblock = (n-1)/maxb + 1
      if (nblock .gt. maxnblock) call errquit
     $     ('uccsdt_block_range: maxnblock too small', nblock, DISK_ERR)
c
c     Compute the number of elements in each block.
c     All blocks get b elements.  (n-nblock*b) blocks get one more.
c
      b = n/nblock
      do i = 1, nblock
         blocks(2,i) = b
      end do
      do i = 1, (n-nblock*b)
         blocks(2,i) = blocks(2,i) + 1
      end do
c
c     Fill in the ranges
c
      blocks(1,1) = lo
      blocks(2,1) = lo + blocks(2,1) - 1
      do i = 2, nblock
         blocks(1,i) = blocks(2,i-1) + 1
         blocks(2,i) = blocks(2,i) + blocks(1,i) - 1
      end do
c
c$$$      write(6,*) ' BLOCKING ', lo, hi
c$$$      write(6,1) nblock, (blocks(1,i),blocks(2,i),i=1,nblock)
c$$$ 1    format(i5,2x,100(2i3,3x))
c
      end
      subroutine moints_shell_block(basis, maxdim, list, nblock)
      implicit none
#include "errquit.fh"
#include "bas.fh"
      integer basis             ! [input]
      integer maxdim            ! [input]
      integer list(2,*)         ! [output]
      integer nblock            ! [output]
c
c     Block the shells so that there are no more than maxdim
c     functions in any block.   
c
c     list(1,block) = first shell in block
c     list(2,block) = last shell in block
c
      integer nsh, ish, count, ilo, ihi, idim
c
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('moints_shell_block: bad basis?',0, BASIS_ERR)
c
      count = 0                 ! No. of functions in current block
      nblock = 1
      list(1,1) = 1
c
      do ish = 1, nsh
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('moints_shell_block: cn2bfr?', ish,
     &       BASIS_ERR)
         idim = ihi - ilo + 1
         if (maxdim .lt. idim) call errquit
     $        ('moints_shell_block: maxdim < idim', maxdim, BASIS_ERR)
c
         if ((count+idim).gt.maxdim) then
            list(2,nblock) = ish-1
            count = 0
            nblock = nblock + 1
            list(1,nblock) = ish
         end if
         count = count + idim
      end do
      if (count .ne. 0) then
         list(2,nblock) = ish-1
         nblock = nblock + 1
      end if
c
      nblock = nblock - 1
c
c$$$      write(6,*) '  shell blocking ', maxdim, nblock
c$$$      write(6,1) (list(1,ish),list(2,ish),ish=1,nblock)
c$$$ 1    format(20(2i3,3x))
c
      end
      subroutine util_ma_ga_avail(avail_ma, avail_ga)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      double precision avail_ma, avail_ga ! [output]
c
c     Returns the available MA and GA on the calling node taking into 
c     account the various interactions of memory
c
c     Units are doubles.
c
c     Within NWChem 
c     
c     if (ga_uses_ma()) then
c     .   GA is initialized unlimited so that can freely switch
c     .   memory between MA and GA
c     else
c     .   GA is initialized limited
c     end if
c
c     but the logic below should take care of the general
c     homgenous case
c
      avail_ma = ma_inquire_avail(mt_dbl)
      if (ga_uses_ma() .and. ga_memory_limited()) 
     $     avail_ma = avail_ma - ga_memory_avail()
      if(ga_memory_limited())then
         avail_ga=ga_memory_avail() / 8.0d0 ! into doubles
      else
         avail_ga=0
      end if
c
      end
      subroutine util_ma_ga_test(need_ma, need_ga, ok)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      double precision need_ma, need_ga ! [input]
      logical ok                ! [output]
c
c     COLLECTIVE
c
c     Given the maximum simultaneous usage of MA and GA **PER NODE**
c     determine if enough memory is available
c
c     Units are doubles.
c
      double precision avail_ma, avail_ga
      integer iok
c
      call util_ma_ga_avail(avail_ma, avail_ga)
      avail_ma = avail_ma * 0.98
      avail_ma = avail_ma * 0.98
c
      if (ga_uses_ma() .and. (.not.ga_memory_limited())) then
         ok = (need_ma + need_ga) .le. avail_ma
      else if (ga_memory_limited())then
         ok = (need_ma .le. avail_ma) .and.
     $        (need_ga .le. avail_ga)
      else 
         ok = (need_ma .le. avail_ma)
      end if
c
      if (ok) then
         iok = 0
      else
         iok = 1
      end if
c
      call igop(33, iok, 1, '+')
      ok = iok .eq. 0
c
      end
      subroutine moints_aoblock_rewind()
      implicit none
#include "cuccsdtP.fh"
c
      aoints_file_offset = 0.0d0
c
      end
      logical function moints_aoblock_read(ush,vsh,aodim,ints)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "eaf.fh"
      integer ush, vsh, aodim
      double precision ints(aodim)
c
      integer err
      character*80 msg
c
c     Needs logic adding to 
c     .  check that the correct data is being read 
c     .  compress the data
c     .  use a GA instead of a file if there is room
c
      if (aoints_file_read) then
c
*         write(6,*) ' READING  ', ush, vsh, aoints_file_offset, aodim
c
         err = eaf_read(aoints_file,aoints_file_offset,ints,8*aodim) 
         if (err .ne. 0) then
            call eaf_errmsg(err, msg)
            write(6,*) ' !!!! EAF error ', msg
            call errquit('moints_aoblock_read: eaf_read?', 0, DISK_ERR)
         end if
         aoints_file_offset = aoints_file_offset + 8*aodim
         moints_aoblock_read = .true.
      else
         moints_aoblock_read = .false.
      end if
c
      end
      subroutine moints_aoblock_write(ush,vsh,aodim,ints)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "eaf.fh"
      integer ush, vsh, aodim
      double precision ints(aodim)
c
      integer err
      character*80 msg
c
c     Needs logic adding to 
c     .  enable a check that the correct data is being read 
c     .  compress the data
c     .  use a GA instead of a file if there is room
c
      if (aoints_file_write) then
c
*         write(6,*) ' WRITING ', ush, vsh, aoints_file_offset, aodim
c
         err = eaf_write(aoints_file,aoints_file_offset,ints,8*aodim) 
         if (err .ne. 0) then
            call eaf_errmsg(err, msg)
            write(6,*) ' !!!! EAF error ', msg
            call errquit('moints_aoblock_write: eaf_read?', 0, DISK_ERR)
         end if
         aoints_file_offset = aoints_file_offset + 8*aodim
      end if
c
      end
C****      subroutine uccsdt_check_list(list)
C****      implicit none
C****#include "errquit.fh"
C****#include "cuccsdtP.fh"
C****#include "mafdecls.fh"
C****      integer list
C****c
C****      integer offset,i,j,k,l,ptr,reclen,ind,l_buf, k_buf
C****      integer blk1, blk2, blk3, blk4
C****      integer spini, spinj, spink, spina, spinb, spinc, spine, spinm, 
C****     $     syme, symm
C****      logical oasym
C****c
C****c     Offset maps a 4-d array into a 1-d array.  It is used to look up
C****c     ptr(ind,i,j,k,l) where ptr is dimensioned ptr(2,dim1,dim2,dim3,dim4) 
C****c     with dim1-4 = listinfo(6-9,list)
C****c
C****      offset(ind,i,j,k,l,list) = 
C****     $     ind-1 + 
C****     $     2*(i-1 + 
C****     $     listinfo(6,list)*(j-1 + 
C****     $     listinfo(7,list)*(k-1 + 
C****     $     listinfo(8,list)*(l-1))))
C****c
C****      if (list.lt.1 .or. list.gt.23) call errquit
C****     $     ('uccsdt_check_list: bad list', list)
C****      write(6,1) list, listname(list), (listinfo(i,list),i=1,9)
C**** 1    format(' list       ', i8, '       name       ', a/
C****     $       ' ptr len    ', i8, '       ptr index  ', i10/
C****     $       ' ptr handle ', i8, '       GA len     ', i8/
C****     $       ' GA handle  ', i8, '       block dims ', 4i5/)
C****c
C****      spina = asuper_spin
C****      spinb = spina
C****      spinc = csuper_spin
C****      spini = spina
C****      spinj = spina
C****      spink = spinc
C****c
C****      oasym = spina .eq. spinc  ! Pure spin integrals are asymmetric
C****c      
C****      do blk4 = 1, listinfo(9,list)
C****         do blk3 = 1, listinfo(8,list)
C****            do blk2 = 1, listinfo(7,list)
C****               do blk1 = 1, listinfo(6,list)
C****                  ptr = int_mb(listinfo(2,list) +
C****     $                 offset(1,blk1,blk2,blk3,blk4,list))
C****                  reclen = int_mb(listinfo(2,list) +
C****     $                 offset(2,blk1,blk2,blk3,blk4,list))
C****                  if (reclen .gt. 0) then
C****                     if (.not. ma_push_get(mt_dbl, reclen, 'buf',
C****     $                    l_buf, k_buf)) call errquit
C****     $                    ('moints_triples_check_list: ma?', reclen)
C****                     call uccsdt_triples_get_list(list,
C****     $                    blk1, blk2, blk3, blk4, dbl_mb(k_buf), reclen)
C****c     
C*****                     write(6,2) blk1, blk2, blk3, blk4, ptr, reclen
C***** 2                   format(' blocks=(',4i5,') ptr=',i8,' reclen=',i8)
C****c
C****                     if (list .eq. 1) then      ! <ei||ab>
C****                        spine = spina
C****                        syme = ieor(oblock(3,blk2,spini),
C****     $                       ieor(ablock(3,blk3),bblock(3,blk4)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
C****     $                       oblock(1,blk2,spini), oblock(2,blk2,spini),
C****     $                       ablock(1,blk3), ablock(2,blk3),
C****     $                       bblock(1,blk4), bblock(2,blk4),
C****     $                       dbl_mb(k_moint), nbf, .true.)
C****                     else if (list .eq. 2) then ! <ek|ac>
C****                        spine = spina
C****                        syme = ieor(oblock(3,blk2,spink),
C****     $                       ieor(ablock(3,blk3),cblock(3,blk4)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
C****     $                       oblock(1,blk2,spink), oblock(2,blk2,spink),
C****     $                       ablock(1,blk3), ablock(2,blk3),
C****     $                       cblock(1,blk4), cblock(2,blk4),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 3) then ! <ek|bc>
C****                        spine = spina
C****                        syme = ieor(oblock(3,blk2,spink),
C****     $                       ieor(bblock(3,blk3),cblock(3,blk4)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
C****     $                       oblock(1,blk2,spink), oblock(2,blk2,spink),
C****     $                       bblock(1,blk3), bblock(2,blk3),
C****     $                       cblock(1,blk4), cblock(2,blk4),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 4) then ! <ei|ca>
C****                        spine = spinc
C****                        syme = ieor(oblock(3,blk2,spini),
C****     $                       ieor(cblock(3,blk3),ablock(3,blk4)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
C****     $                       oblock(1,blk2,spini), oblock(2,blk2,spini),
C****     $                       cblock(1,blk3), cblock(2,blk3),
C****     $                       ablock(1,blk4), ablock(2,blk4),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 5) then ! <ei|cb>
C****                        spine = spinc
C****                        syme = ieor(oblock(3,blk2,spini),
C****     $                       ieor(cblock(3,blk3),bblock(3,blk4)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
C****     $                       oblock(1,blk2,spini), oblock(2,blk2,spini),
C****     $                       cblock(1,blk3), cblock(2,blk3),
C****     $                       bblock(1,blk4), bblock(2,blk4),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 6) then ! <mc|ik>
C****                        spinm = spina
C****                        symm  = ieor(oblock(3,blk4,spink),
C****     $                       ieor(oblock(3,blk3,spini),cblock(3,blk2)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
C****     $                       cblock(1,blk2), cblock(2,blk2),
C****     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
C****     $                       oblock(1,blk4,spink), oblock(2,blk4,spink),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 7) then ! <ma||ij>
C****                        spinm = spina
C****                        symm  = ieor(oblock(3,blk4,spini),
C****     $                       ieor(oblock(3,blk3,spini),ablock(3,blk2)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
C****     $                       ablock(1,blk2), ablock(2,blk2),
C****     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
C****     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
C****     $                       dbl_mb(k_moint), nbf, .true.)
C****                     else if (list .eq. 8) then ! <mb||ij>
C****                        spinm = spina
C****                        symm  = ieor(oblock(3,blk4,spini),
C****     $                       ieor(oblock(3,blk3,spini),bblock(3,blk2)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
C****     $                       bblock(1,blk2), bblock(2,blk2),
C****     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
C****     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
C****     $                       dbl_mb(k_moint), nbf, .true.)
C****                     else if (list .eq. 9) then ! <ma|ki>
C****                        spinm = spink
C****                        symm  = ieor(oblock(3,blk4,spini),
C****     $                       ieor(oblock(3,blk3,spink),ablock(3,blk2)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
C****     $                       ablock(1,blk2), ablock(2,blk2),
C****     $                       oblock(1,blk3,spink), oblock(2,blk3,spink),
C****     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 10) then ! <mb|ki>
C****                        spinm = spink
C****                        symm  = ieor(oblock(3,blk4,spini),
C****     $                       ieor(oblock(3,blk3,spink),bblock(3,blk2)))
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
C****     $                       bblock(1,blk2), bblock(2,blk2),
C****     $                       oblock(1,blk3,spink), oblock(2,blk3,spink),
C****     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 11) then ! <ab||ij>
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       ablock(1,blk1), ablock(2,blk1),
C****     $                       bblock(1,blk2), bblock(2,blk2),
C****     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
C****     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
C****     $                       dbl_mb(k_moint), nbf, .true.)
C****                     else if (list .eq. 12) then ! <ac|ik>
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       ablock(1,blk1), ablock(2,blk1),
C****     $                       cblock(1,blk2), cblock(2,blk2),
C****     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
C****     $                       oblock(1,blk4,spink), oblock(2,blk4,spink),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else if (list .eq. 13) then ! <bc|ik>
C****                        call uccsdt_do_check_list(
C****     $                       dbl_mb(k_buf), 
C****     $                       bblock(1,blk1), bblock(2,blk1),
C****     $                       cblock(1,blk2), cblock(2,blk2),
C****     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
C****     $                       oblock(1,blk4,spink), oblock(2,blk4,spink),
C****     $                       dbl_mb(k_moint), nbf, oasym)
C****                     else
C****                        call errquit('cannot check this list',0)
C****                     endif
C****                     if (.not. ma_pop_stack(l_buf)) call errquit
C****     $                    ('uccsdt_check_list: ma corrupt?',0)
C****                  end if
C****               end do
C****            end do
C****         end do
C****      end do
C****c$$$      if (list .eq. 1) then ! <ei||ab>
C****c$$$      else if (list .eq. 2) then ! <ek|ac>
C****c$$$      else if (list .eq. 3) then ! <ek|bc>
C****c$$$      else if (list .eq. 4) then ! <ei|ca>
C****c$$$      else if (list .eq. 5) then ! <ei|cb>
C****c$$$      else if (list .eq. 6) then ! <mc|ik>
C****c$$$      else if (list .eq. 7) then ! <ma||ij>
C****c$$$      else if (list .eq. 8) then ! <mb||ij>
C****c$$$      else if (list .eq. 9) then ! <ma|ki>
C****c$$$      else if (list .eq. 10) then ! <mb|ki>
C****c$$$      else if (list .eq. 11) then ! <ab||ij>
C****c$$$      else if (list .eq. 12) then ! <ac|ik>
C****c$$$      else if (list .eq. 13) then ! <bc|ik>
C****c$$$      else if (list .eq. 14) then ! t(e,c,i,k)
C****c$$$      else if (list .eq. 15) then ! t(e,a,i,j)
C****c$$$      else if (list .eq. 16) then ! t(e,b,i,j)
C****c$$$      else if (list .eq. 17) then ! t(e,a,k,i)
C****c$$$      else if (list .eq. 18) then ! t(e,b,k,i)
C****c$$$      else if (list .eq. 19) then ! t(m,i,a,b)
C****c$$$      else if (list .eq. 20) then ! t(m,k,a,c)
C****c$$$      else if (list .eq. 21) then ! t(m,k,b,c)
C****c$$$      else if (list .eq. 22) then ! t(m,i,c,a)
C****c$$$      else if (list .eq. 23) then ! t(m,i,c,b)
C****c$$$      else
C****c$$$         call errquit('uccsdt_triples_check_list: bad list', list)
C****c$$$      end if
C****c      
C****      end
C****      subroutine uccsdt_triples_dummy_amplitudes()
C****      implicit none
C****#include "errquit.fh"
C****#include "cuccsdtP.fh"
C****#include "mafdecls.fh"
C****c
C****c     Set the dummy t2 amplitudes to 
C****c
C****c     t(a,b,i,j) = <ab||ij> / (eps(i)+eps(j)-eps(a)-eps(b))
C****c
C****      integer spina, spinc, spini, spink, spine, spinm, list
C****      integer k_ptr(2:4), k_t2
C****      logical oasym
C****c
C****      spina = asuper_spin
C****      spinc = csuper_spin
C****      spini = spina
C****      spink = spinc
C****      oasym = spina .eq. spinc
C****c
C****      k_ptr(2) = k_t2aa
C****      k_ptr(3) = k_t2ab
C****      k_ptr(4) = k_t2bb
C****c
C****c      write(6,*)' k_t2aa, k_t2ab, k_t2bb, k_ptr(2), k_ptr(3), k_ptr(4)',
C****c     &            k_t2aa, k_t2ab, k_t2bb, k_ptr(2), k_ptr(3), k_ptr(4)
C****c
C****c      call jan_debug_print('k_t2aa', dbl_mb(k_t2aa), nbf, nbf,  
C****c     $                      nbf, nbf)
C****c
C****c      call jan_debug_print('k_t2ab', dbl_mb(k_t2ab), nbf, nbf,  
C****c     $                      nbf, nbf)
C****c
C****c      call jan_debug_print('k_t2bb', dbl_mb(k_t2bb), nbf, nbf,  
C****c     $                      nbf, nbf)
C****c
C****c      write(6,*)' spina, spinc, spini, spink ',
C****c     &            spina, spinc, spini, spink
C****c
C****cccccccc      if (spina .ne. spinc) call errquit('dummy_amplitudes: no mixed',0)
C****c
C****c     14.  t(e,c,i,k) spin(e)=spin(a) mixed
C****c
C****      list = 14
C****      spine = spina
C****      k_t2 = k_ptr(spine+spinc)
C****c      write(6,*)' list, spine, spinc, k_t2 ', 
C****c     &            list, spine, spinc, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     v_sym(1,0,spine), eps(1,spine),
C****     $     ncblock, cblock,  eps(1,spinc),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
C****     $     dbl_mb(k_t2), nbf, k_t2, spine, spinc,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, oasym)
C****c
C****c     15.  t(e,a,i,j) spin(e)=spin(a) pure
C****c
C****      list = 15
C****      spine = spina
C****      k_t2 = k_ptr(spine+spina)
C****c      write(6,*)' list, spine, spina, k_t2 ', 
C****c     &            list, spine, spina, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     v_sym(1,0,spine), eps(1,spine),
C****     $     nablock, ablock,  eps(1,spina),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     dbl_mb(k_t2), nbf, k_t2, spine, spina,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, .true.)
C****c
C****c     16.  t(e,b,i,j) spin(e)=spin(a) pure
C****c
C****      list = 16
C****      spine = spina
C****      k_t2 = k_ptr(spine+spina)
C****c      write(6,*)' list, spine, spina, k_t2 ', 
C****c     &            list, spine, spina, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     v_sym(1,0,spine), eps(1,spine),
C****     $     nbblock, bblock,  eps(1,spina),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     dbl_mb(k_t2), nbf, k_t2, spine, spina,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, .true.)
C****c
C****c     17.  t(e,a,k,i) spin(e)=spin(c) mixed
C****c
C****      list = 17
C****      spine = spinc
C****      k_t2 = k_ptr(spine+spina)
C****c      write(6,*)' list, spine, spina, k_t2 ', 
C****c     &            list, spine, spina, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     v_sym(1,0,spine), eps(1,spine),
C****     $     nablock, ablock,  eps(1,spina),
C****     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     dbl_mb(k_t2), nbf, k_t2, spine, spina,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, oasym)
C****c
C****c     18.  t(e,b,k,i) spin(e)=spin(c) mixed
C****c
C****      list = 18
C****      spine = spinc
C****      k_t2 = k_ptr(spine+spina)
C****c      write(6,*)' list, spine, spina, k_t2 ', 
C****c     &            list, spine, spina, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     v_sym(1,0,spine), eps(1,spine),
C****     $     nbblock, bblock,  eps(1,spina),
C****     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     dbl_mb(k_t2), nbf, k_t2, spine, spina,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, oasym)
C****c
C****c     19.  t(m,i,a,b) spin(m)=spin(a) pure
C****c
C****      list = 19
C****      spinm = spina
C****      k_t2 = k_ptr(spinm+spini)
C****c      write(6,*)' list, spinm, spini, k_t2 ', 
C****c     &            list, spinm, spini, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     o_sym(1,0,spinm), eps(1,spinm),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     nablock, ablock,  eps(1,spina),
C****     $     nbblock, bblock,  eps(1,spina),
C****     $     dbl_mb(k_t2), nbf, k_t2, spinm, spini,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, .true.)
C****c
C****c     20.  t(m,k,a,c) spin(m)=spin(a) mixed
C****c
C****      list = 20
C****      spinm = spina
C****      k_t2 = k_ptr(spinm+spink)
C****c      write(6,*)' list, spinm, spink, k_t2 ', 
C****c     &            list, spinm, spink, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     o_sym(1,0,spinm), eps(1,spinm),
C****     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
C****     $     nablock, ablock,  eps(1,spina),
C****     $     ncblock, cblock,  eps(1,spinc),
C****     $     dbl_mb(k_t2), nbf, k_t2, spinm, spink,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
C****c
C****c     21.  t(m,k,b,c) spin(m)=spin(a) mixed
C****c
C****      list = 21
C****      spinm = spina
C****      k_t2 = k_ptr(spinm+spink)
C****c      write(6,*)' list, spinm, spink, k_t2 ', 
C****c     &            list, spinm, spink, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     o_sym(1,0,spinm), eps(1,spinm),
C****     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
C****     $     nbblock, bblock,  eps(1,spina),
C****     $     ncblock, cblock,  eps(1,spinc),
C****     $     dbl_mb(k_t2), nbf, k_t2, spinm, spink,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
C****c
C****c     22.  t(m,i,c,a) spin(m)=spin(c) mixed
C****c
C****      list = 22
C****      spinm = spinc
C****      k_t2 = k_ptr(spinm+spini)
C****c      write(6,*)' list, spinm, spini, k_t2 ', 
C****c     &            list, spinm, spini, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     o_sym(1,0,spinm), eps(1,spinm),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     ncblock, cblock,  eps(1,spinc),
C****     $     nablock, ablock,  eps(1,spina),
C****     $     dbl_mb(k_t2), nbf, k_t2, spinm, spini,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
C****c
C****c     23.  t(m,i,c,b) spin(m)=spin(c) mixed
C****c
C****      list = 23
C****      spinm = spinc
C****      k_t2 = k_ptr(spinm+spini)
C****c      write(6,*)' list, spinm, spini, k_t2 ', 
C****c     &            list, spinm, spini, k_t2
C****      call uccsdt_do_triples_dummy_amplitudes(
C****     $     o_sym(1,0,spinm), eps(1,spinm),
C****     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
C****     $     ncblock, cblock,  eps(1,spinc),
C****     $     nbblock, bblock,  eps(1,spina),
C****     $     dbl_mb(k_t2), nbf, k_t2, spinm, spini,
C****     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
C****c
C****      end

         
c $Id$
