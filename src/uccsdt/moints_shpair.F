      subroutine uccsdt_check_list(list)
      implicit none
#include "cuccsdtP.fh"
#include "mafdecls.fh"
      integer list
c
      integer offset,i,j,k,l,ptr,reclen,ind,l_buf, k_buf
      integer blk1, blk2, blk3, blk4
      integer spini, spinj, spink, spina, spinb, spinc, spine, spinm, 
     $     syme, symm
      logical oasym
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
     $     ('uccsdt_check_list: bad list', list)
      write(6,1) list, listname(list), (listinfo(i,list),i=1,9)
 1    format(' list       ', i8, '       name       ', a/
     $       ' ptr len    ', i8, '       ptr index  ', i10/
     $       ' ptr handle ', i8, '       GA len     ', i8/
     $       ' GA handle  ', i8, '       block dims ', 4i5/)
c
      spina = asuper_spin
      spinb = bsuper_spin
      spinc = csuper_spin
      spini = spina
      spinj = spina
      spink = spinc
c
      oasym = spina .eq. spinc  ! Pure spin integrals are asymmetric
c      
      do blk4 = 1, listinfo(9,list)
         do blk3 = 1, listinfo(8,list)
            do blk2 = 1, listinfo(7,list)
               do blk1 = 1, listinfo(6,list)
                  ptr = int_mb(listinfo(2,list) +
     $                 offset(1,blk1,blk2,blk3,blk4,list))
                  reclen = int_mb(listinfo(2,list) +
     $                 offset(2,blk1,blk2,blk3,blk4,list))
                  if (reclen .gt. 0) then
                     if (.not. ma_push_get(mt_dbl, reclen, 'buf',
     $                    l_buf, k_buf)) call errquit
     $                    ('moints_triples_check_list: ma?', reclen)
                     call uccsdt_triples_get_list(list,
     $                    blk1, blk2, blk3, blk4, dbl_mb(k_buf), reclen)
c     
*                     write(6,2) blk1, blk2, blk3, blk4, ptr, reclen
* 2                   format(' blocks=(',4i5,') ptr=',i8,' reclen=',i8)
c
                     if (list .eq. 1) then      ! <ei||ab>
                        spine = spina
                        syme = ieor(oblock(3,blk2,spini),
     $                       ieor(ablock(3,blk3),bblock(3,blk4)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
     $                       oblock(1,blk2,spini), oblock(2,blk2,spini),
     $                       ablock(1,blk3), ablock(2,blk3),
     $                       bblock(1,blk4), bblock(2,blk4),
     $                       dbl_mb(k_moint), nbf, .true.)
                     else if (list .eq. 2) then ! <ek|ac>
                        spine = spina
                        syme = ieor(oblock(3,blk2,spink),
     $                       ieor(ablock(3,blk3),cblock(3,blk4)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
     $                       oblock(1,blk2,spink), oblock(2,blk2,spink),
     $                       ablock(1,blk3), ablock(2,blk3),
     $                       cblock(1,blk4), cblock(2,blk4),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 3) then ! <ek|bc>
                        spine = spina
                        syme = ieor(oblock(3,blk2,spink),
     $                       ieor(bblock(3,blk3),cblock(3,blk4)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
     $                       oblock(1,blk2,spink), oblock(2,blk2,spink),
     $                       bblock(1,blk3), bblock(2,blk3),
     $                       cblock(1,blk4), cblock(2,blk4),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 4) then ! <ei|ca>
                        spine = spinc
                        syme = ieor(oblock(3,blk2,spini),
     $                       ieor(cblock(3,blk3),ablock(3,blk4)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
     $                       oblock(1,blk2,spini), oblock(2,blk2,spini),
     $                       cblock(1,blk3), cblock(2,blk3),
     $                       ablock(1,blk4), ablock(2,blk4),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 5) then ! <ei|cb>
                        spine = spinc
                        syme = ieor(oblock(3,blk2,spini),
     $                       ieor(cblock(3,blk3),bblock(3,blk4)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       v_sym(1,syme,spine), v_sym(2,syme,spine),
     $                       oblock(1,blk2,spini), oblock(2,blk2,spini),
     $                       cblock(1,blk3), cblock(2,blk3),
     $                       bblock(1,blk4), bblock(2,blk4),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 6) then ! <mc|ik>
                        spinm = spina
                        symm  = ieor(oblock(3,blk4,spink),
     $                       ieor(oblock(3,blk3,spini),cblock(3,blk2)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
     $                       cblock(1,blk2), cblock(2,blk2),
     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
     $                       oblock(1,blk4,spink), oblock(2,blk4,spink),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 7) then ! <ma||ij>
                        spinm = spina
                        symm  = ieor(oblock(3,blk4,spini),
     $                       ieor(oblock(3,blk3,spini),ablock(3,blk2)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
     $                       ablock(1,blk2), ablock(2,blk2),
     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
     $                       dbl_mb(k_moint), nbf, .true.)
                     else if (list .eq. 8) then ! <mb||ij>
                        spinm = spina
                        symm  = ieor(oblock(3,blk4,spini),
     $                       ieor(oblock(3,blk3,spini),bblock(3,blk2)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
     $                       bblock(1,blk2), bblock(2,blk2),
     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
     $                       dbl_mb(k_moint), nbf, .true.)
                     else if (list .eq. 9) then ! <ma|ki>
                        spinm = spink
                        symm  = ieor(oblock(3,blk4,spini),
     $                       ieor(oblock(3,blk3,spink),ablock(3,blk2)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
     $                       ablock(1,blk2), ablock(2,blk2),
     $                       oblock(1,blk3,spink), oblock(2,blk3,spink),
     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 10) then ! <mb|ki>
                        spinm = spink
                        symm  = ieor(oblock(3,blk4,spini),
     $                       ieor(oblock(3,blk3,spink),bblock(3,blk2)))
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       o_sym(1,symm,spinm), o_sym(2,symm,spinm),
     $                       bblock(1,blk2), bblock(2,blk2),
     $                       oblock(1,blk3,spink), oblock(2,blk3,spink),
     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 11) then ! <ab||ij>
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       ablock(1,blk1), ablock(2,blk1),
     $                       bblock(1,blk2), bblock(2,blk2),
     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
     $                       oblock(1,blk4,spini), oblock(2,blk4,spini),
     $                       dbl_mb(k_moint), nbf, .true.)
                     else if (list .eq. 12) then ! <ac|ik>
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       ablock(1,blk1), ablock(2,blk1),
     $                       cblock(1,blk2), cblock(2,blk2),
     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
     $                       oblock(1,blk4,spink), oblock(2,blk4,spink),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else if (list .eq. 13) then ! <bc|ik>
                        call uccsdt_do_check_list(
     $                       dbl_mb(k_buf), 
     $                       bblock(1,blk1), bblock(2,blk1),
     $                       cblock(1,blk2), cblock(2,blk2),
     $                       oblock(1,blk3,spini), oblock(2,blk3,spini),
     $                       oblock(1,blk4,spink), oblock(2,blk4,spink),
     $                       dbl_mb(k_moint), nbf, oasym)
                     else
                        call errquit('cannot check this list',0)
                     endif
                     if (.not. ma_pop_stack(l_buf)) call errquit
     $                    ('uccsdt_check_list: ma corrupt?',0)
                  end if
               end do
            end do
         end do
      end do
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
c$$$         call errquit('uccsdt_triples_check_list: bad list', list)
c$$$      end if
c      
      end
      subroutine uccsdt_triples_dummy_amplitudes()
      implicit none
#include "cuccsdtP.fh"
#include "mafdecls.fh"
c
c     Set the dummy t2 amplitudes to 
c
c     t(a,b,i,j) = <ab||ij> / (eps(i)+eps(j)-eps(a)-eps(b))
c
      integer spina, spinc, spini, spink, spine, spinm, list
      logical oasym
c
      spina = asuper_spin
      spinc = csuper_spin
      spini = spina
      spink = spini
      oasym = spina .eq. spinc
c
cccccccc      if (spina .ne. spinc) call errquit('dummy_amplitudes: no mixed',0)
c
c     14.  t(e,c,i,k) spin(e)=spin(a) mixed
c
      list = 14
      spine = spina
      call uccsdt_do_triples_dummy_amplitudes(
     $     v_sym(1,0,spine), eps(1,spine),
     $     ncblock, cblock,  eps(1,spinc),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, oasym)
c
c     15.  t(e,a,i,j) spin(e)=spin(a) pure
c
      list = 15
      spine = spina
      call uccsdt_do_triples_dummy_amplitudes(
     $     v_sym(1,0,spine), eps(1,spine),
     $     nablock, ablock,  eps(1,spina),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, .true.)
c
c     16.  t(e,b,i,j) spin(e)=spin(a) pure
c
      list = 16
      spine = spina
      call uccsdt_do_triples_dummy_amplitudes(
     $     v_sym(1,0,spine), eps(1,spine),
     $     nbblock, bblock,  eps(1,spina),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, .true.)
c
c     17.  t(e,a,k,i) spin(e)=spin(c) mixed
c
      list = 17
      spine = spinc
      call uccsdt_do_triples_dummy_amplitudes(
     $     v_sym(1,0,spine), eps(1,spine),
     $     nablock, ablock,  eps(1,spina),
     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, oasym)
c
c     18.  t(e,b,k,i) spin(e)=spin(c) mixed
c
      list = 18
      spine = spinc
      call uccsdt_do_triples_dummy_amplitudes(
     $     v_sym(1,0,spine), eps(1,spine),
     $     nbblock, bblock,  eps(1,spina),
     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), -1d0, oasym)
c
c     19.  t(m,i,a,b) spin(m)=spin(a) pure
c
      list = 19
      spinm = spina
      call uccsdt_do_triples_dummy_amplitudes(
     $     o_sym(1,0,spinm), eps(1,spinm),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     nablock, ablock,  eps(1,spina),
     $     nbblock, bblock,  eps(1,spina),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, .true.)
c
c     20.  t(m,k,a,c) spin(m)=spin(a) mixed
c
      list = 20
      spinm = spina
      call uccsdt_do_triples_dummy_amplitudes(
     $     o_sym(1,0,spinm), eps(1,spinm),
     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
     $     nablock, ablock,  eps(1,spina),
     $     ncblock, cblock,  eps(1,spinc),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
c
c     21.  t(m,k,b,c) spin(m)=spin(a) mixed
c
      list = 21
      spinm = spina
      call uccsdt_do_triples_dummy_amplitudes(
     $     o_sym(1,0,spinm), eps(1,spinm),
     $     noblock(spink), oblock(1,1,spink),  eps(1,spink),
     $     nbblock, bblock,  eps(1,spina),
     $     ncblock, cblock,  eps(1,spinc),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
c
c     22.  t(m,i,c,a) spin(m)=spin(c) mixed
c
      list = 22
      spinm = spinc
      call uccsdt_do_triples_dummy_amplitudes(
     $     o_sym(1,0,spinm), eps(1,spinm),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     ncblock, cblock,  eps(1,spinc),
     $     nablock, ablock,  eps(1,spina),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
c
c     23.  t(m,i,c,b) spin(m)=spin(c) mixed
c
      list = 23
      spinm = spinc
      call uccsdt_do_triples_dummy_amplitudes(
     $     o_sym(1,0,spinm), eps(1,spinm),
     $     noblock(spini), oblock(1,1,spini),  eps(1,spini),
     $     ncblock, cblock,  eps(1,spinc),
     $     nbblock, bblock,  eps(1,spina),
     $     dbl_mb(k_moint), nbf,
     $     int_mb(listinfo(2,list)), listinfo(5,list), 1d0, oasym)
c
      end
      subroutine uccsdt_do_triples_dummy_amplitudes(
     $     e_sym, eeps,
     $     ncblock, cblock, ceps,
     $     niblock, iblock, ieps,
     $     nkblock, kblock, keps,
     $     moints, nbf,
     $     pointers, g_list, phase, oasym)
      implicit none
      integer e_sym(2,0:7)
      integer ncblock, cblock(3,*)
      integer niblock, iblock(3,*)
      integer nkblock, kblock(3,*)
      integer nbf
      double precision eeps(*), ceps(*), ieps(*), keps(*)
      double precision moints(nbf, nbf, nbf, nbf)
      integer pointers(2,ncblock,niblock,nkblock)
      integer g_list
      double precision phase
      logical oasym
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
                           t = moints(e,c,i,k)
                           if (oasym) t = t - moints(e,c,k,i)
                           t = t / (eeps(e)+ceps(c)-ieps(i)-keps(k))
                           t = t * phase
                           call ga_put(g_list,ptr,ptr,1,1,t,1)
                           ptr = ptr + 1
                        end do
                     end do
                  end do
               end do
               if (ptr .ne. (pointers(1,cblk,iblk,kblk)+
     $              pointers(2,cblk,iblk,kblk))) call errquit
     $              ('bad length for t record',ptr)
            end do
         end do
      end do
c
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
     $     ('uccsdt_triples_get_list: bad list', list)
      if (blk1.lt.1 .or. blk1.gt.listinfo(6,list)) call errquit
     $     ('uccsdt_triples_get_list: block 1 bad', list)
      if (blk2.lt.1 .or. blk2.gt.listinfo(7,list)) call errquit
     $     ('uccsdt_triples_get_list: block 2 bad', list)
      if (blk3.lt.1 .or. blk3.gt.listinfo(8,list)) call errquit
     $     ('uccsdt_triples_get_list: block 3 bad', list)
      if (blk4.lt.1 .or. blk4.gt.listinfo(9,list)) call errquit
     $     ('uccsdt_triples_get_list: block 4 bad', list)
c
      ptr = int_mb(listinfo(2,list) + 
     $     offset(1,blk1,blk2,blk3,blk4,list))
      reclen = int_mb(listinfo(2,list) +
     $     offset(2,blk1,blk2,blk3,blk4,list))
c
      if (reclen .gt. lenbuf) call errquit
     $     ('uccsdt_triples_get_list: lenbuf too small', reclen)
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
         
      subroutine uccsdt_triples_set_super_block
     $     (alo, ahi, blo, bhi, clo, chi, spina, spinb, spinc,
     $     oprint)
      implicit none
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer alo, ahi, blo, bhi, clo, chi ! [input]
      integer spina, spinb, spinc
      logical oprint
c
c     The outer-3-most loops of the triples are over blocks of virtual 
c     indices.  The blocking is determined outside of this routine
c     so as to minimize recomputation of the integrals within memory 
c     constraints.
c
c     This routine, given the specified super blocks and spins of 
c     abc, sets up all of the blocking info in the common block
c
c     1. The sub-blocking info for both virtual and occupied indices
c     2. The addressing info for indexing into the 5 global
c     .  arrays by sub-block indices
c     3. Allocation of the global arrays themselves
c
      integer ir, a1, a2, b1, b2, c1, c2, spin, list, i
      integer spini, spinj, spink
c
      asuper(1) = alo
      asuper(2) = ahi
      bsuper(1) = blo
      bsuper(2) = bhi
      csuper(1) = clo
      csuper(2) = chi
      asuper_spin = spina
      bsuper_spin = spinb
      csuper_spin = spinc
c
c     Determine the optimal sub-blocking for the occupied
c     and virtual indices ... these should be determined by the
c     machine characteristics and the performance model ... 
c     for debug purposes just hardwire these.
c
c     Must also ensure that there are not so many sub-blocks as
c     to overflow the dimensions of the adressing arrays.
c
      oblock_max = 2
      vblock_max = 3
c
c     Sub-block the 'a' super-block ... sub-blocks will be restricted
c     to orbitals of the same symmetry whereas super-blocks can contain
c     multiple symmetries and split symmetry blocks.
c
      nablock = 0
      do ir = 0, nir-1
         do a1 = max(alo,v_sym(1,ir,spina)),
     $           min(ahi,v_sym(2,ir,spina)), vblock_max
            a2 = a1 + vblock_max - 1
            a2 = min(a2,min(ahi,v_sym(2,ir,spina)))
            nablock = nablock + 1
            ablock(1,nablock) = a1
            ablock(2,nablock) = a2
            ablock(3,nablock) = ir
         end do
      end do
c
c     Ditto for b
c
      nbblock = 0
      do ir = 0, nir-1
         do b1 = max(blo,v_sym(1,ir,spinb)),
     $           min(bhi,v_sym(2,ir,spinb)), vblock_max
            b2 = b1 + vblock_max - 1
            b2 = min(b2,min(bhi,v_sym(2,ir,spinb)))
            nbblock = nbblock + 1
            bblock(1,nbblock) = b1
            bblock(2,nbblock) = b2
            bblock(3,nbblock) = ir
         end do
      end do
c
c     Ditto for c
c
      ncblock = 0
      do ir = 0, nir-1
         do c1 = max(clo,v_sym(1,ir,spinc)),
     $           min(chi,v_sym(2,ir,spinc)), vblock_max
            c2 = c1 + vblock_max - 1
            c2 = min(c2,min(chi,v_sym(2,ir,spinc)))
            ncblock = ncblock + 1
            cblock(1,ncblock) = c1
            cblock(2,ncblock) = c2
            cblock(3,ncblock) = ir
         end do
      end do
c
c     The blocking of the occupied never changes but just re-do it 
c     every time for convenience
c
      do spin = 1, 2
         noblock(spin) = 0
         do ir = 0, nir-1
            do i = o_sym(1,ir,spin),o_sym(2,ir,spin),oblock_max
               noblock(spin) = noblock(spin) + 1
               oblock(1,noblock(spin),spin) = i
               oblock(2,noblock(spin),spin) = 
     $              min(o_sym(2,ir,spin),i+oblock_max-1)
               oblock(3,noblock(spin),spin) = ir
            end do
         end do
      end do
c
      if (oprint) then
         write(6,*) ' A super block ', alo, ahi, spina, nablock
         write(6,1) (i,ablock(1,i),ablock(2,i),ablock(3,i),i=1,nablock)
         write(6,*) ' B super block ', blo, bhi, spinb, nbblock
         write(6,1) (i,bblock(1,i),bblock(2,i),bblock(3,i),i=1,nbblock)
         write(6,*) ' C super block ', clo, chi, spinc, ncblock
         write(6,1) (i,cblock(1,i),cblock(2,i),cblock(3,i),i=1,ncblock)
         write(6,*) ' Alpha occupied blocking'
         write(6,1) (i,oblock(1,i,1),oblock(2,i,1),oblock(3,i,1),
     $        i=1,noblock(1))
         write(6,*) ' Beta  occupied blocking'
         write(6,1) (i,oblock(1,i,2),oblock(2,i,2),oblock(3,i,2),
     $        i=1,noblock(2))
 1       format(4x,i5,4x,3i5)
      endif
c
c     Now set up the addressing arrays for the by-sub-block
c     addressing of the 1-d global arrays with contiguous sub-blocks
c
c     Mixed spin case has
c     .   spin(a)=spin(b)=spin(i)=spin(j)
c     .   spin(c)=spin(k)
c
c     Pure spin case has all spins the same but the same lists
c     are created as for mixed spin, except that all integrals
c     are antisymmetrized.
c
c     3-x integrals
c     1.  <ei||ab>  spin(e)=spin(a)  pure
c     2.  <ek|ac>   spin(e)=spin(a)  mixed 
c     3.  <ek|bc>   spin(e)=spin(a)  mixed
c     4.  <ei|ca>   spin(e)=spin(c)  mixed
c     5.  <ei|cb>   spin(e)=spin(c)  mixed
c
c     1-x integrals
c     6.  <mc|ik>  spin(m)=spin(a)  mixed
c     7.  <ma||ij> spin(m)=spin(a)  pure
c     8.  <mb||ij> spin(m)=spin(a)  pure
c     9.  <ma|ki>  spin(m)=spin(c)  mixed
c     10. <mb|ki>  spin(m)=spin(c)  mixed
c
c     2-x integrals (for V)
c     11. <ab||ij> pure
c     12. <ac|ik>  mixed
c     13. <bc|ik>  mixed
c
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
c     14.  t(e,c,i,k) spin(e)=spin(a) mixed
c     15.  t(e,a,i,j) spin(e)=spin(a) pure
c     16.  t(e,b,i,j) spin(e)=spin(a) pure
c     17.  t(e,a,k,i) spin(e)=spin(c) mixed
c     18.  t(e,b,k,i) spin(e)=spin(c) mixed
c
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
c     19.  t(m,i,a,b) spin(m)=spin(a) pure
c     20.  t(m,k,a,c) spin(m)=spin(a) mixed
c     21.  t(m,k,b,c) spin(m)=spin(a) mixed
c     22.  t(m,i,c,a) spin(m)=spin(c) mixed
c     23.  t(m,i,c,b) spin(m)=spin(c) mixed
c
      spini = spina
      spinj = spina
      spink = spinc
c
c     Setup the names of the lists for easier printing & debugging
c
c     3-x integrals
      listname(1) = '<ei||ab>'
      listname(2) = '<ek|ac>'
      listname(3) = '<ek|bc>'
      listname(4) = '<ei|ca>'
      listname(5) = '<ei|cb>'
c     1-x integrals
      listname(6) = '<mc|ik>'
      listname(7) = '<ma||ij>'
      listname(8) = '<mb||ij>'
      listname(9) = '<ma|ki>'
      listname(10) = '<mb|ki>'
c     2-x integrals (for V)
      listname(11) = '<ab||ij>'
      listname(12) = '<ac|ik>'
      listname(13) = '<bc|ik>'
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
      listname(14) = 't(e,c,i,k)'
      listname(15) = 't(e,a,i,j)'
      listname(16) = 't(e,b,i,j)'
      listname(17) = 't(e,a,k,i)'
      listname(18) = 't(e,b,k,i)'
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
      listname(19) = 't(m,i,a,b)'
      listname(20) = 't(m,k,a,c)'
      listname(21) = 't(m,k,b,c)'
      listname(22) = 't(m,i,c,a)'
      listname(23) = 't(m,i,c,b)'
c
c     Compute the size of each dimension of the pointer arrays
c     and then computer their size so they can be MA'd
c     
c     3-x integrals
      listinfo(6,1) = 1                         ! 1.  <ei||ab> 
      listinfo(7,1) = noblock(spini)
      listinfo(8,1) = nablock
      listinfo(9,1) = nbblock                   
      listinfo(6,2) = 1                         ! 2.  <ek|ac>  
      listinfo(7,2) = noblock(spink)
      listinfo(8,2) = nablock
      listinfo(9,2) = ncblock                   
      listinfo(6,3) = 1                         ! 3.  <ek|bc>  
      listinfo(7,3) = noblock(spink)
      listinfo(8,3) = nbblock
      listinfo(9,3) = ncblock                   
      listinfo(6,4) = 1                         ! 4.  <ei|ca>  
      listinfo(7,4) = noblock(spini)
      listinfo(8,4) = ncblock
      listinfo(9,4) = nablock                   
      listinfo(6,5) = 1                         ! 5.  <ei|cb>  
      listinfo(7,5) = noblock(spini)
      listinfo(8,5) = ncblock
      listinfo(9,5) = nbblock                   
c     1-x integrals
      listinfo(6,6) = 1                         ! 6.  <mc|ik>
      listinfo(7,6) = ncblock
      listinfo(8,6) = noblock(spini)
      listinfo(9,6) = noblock(spink)            
      listinfo(6,7) = 1                         ! 7.  <ma||ij>
      listinfo(7,7) = nablock
      listinfo(8,7) = noblock(spini)
      listinfo(9,7) = noblock(spinj)            
      listinfo(6,8) = 1                         ! 8.  <mb||ij>
      listinfo(7,8) = nbblock
      listinfo(8,8) = noblock(spini)
      listinfo(9,8) = noblock(spinj)            
      listinfo(6,9) = 1                         ! 9.  <ma|ki>
      listinfo(7,9) = nablock
      listinfo(8,9) = noblock(spink)
      listinfo(9,9) = noblock(spini)            
      listinfo(6,10)= 1                         ! 10. <mb|ki> 
      listinfo(7,10)= nbblock
      listinfo(8,10)= noblock(spink)
      listinfo(9,10)= noblock(spini)            
c     2-x integrals (for V)
      listinfo(6,11)= nablock                   ! 11. <ab||ij>
      listinfo(7,11)= nbblock
      listinfo(8,11)= noblock(spini)
      listinfo(9,11)= noblock(spinj)            
      listinfo(6,12)= nablock                   ! 12. <ac|ik>
      listinfo(7,12)= ncblock
      listinfo(8,12)= noblock(spini)
      listinfo(9,12)= noblock(spink)            
      listinfo(6,13)= nbblock                   ! 13. <bc|ik>
      listinfo(7,13)= ncblock
      listinfo(8,13)= noblock(spini)
      listinfo(9,13)= noblock(spink)            
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
      listinfo(6,14)= 1                         ! 14.  t(e,c,i,k)
      listinfo(7,14)= ncblock
      listinfo(8,14)= noblock(spini)
      listinfo(9,14)= noblock(spink)            
      listinfo(6,15)= 1                         ! 15.  t(e,a,i,j)
      listinfo(7,15)= nablock
      listinfo(8,15)= noblock(spini)
      listinfo(9,15)= noblock(spinj)            
      listinfo(6,16)= 1                         ! 16.  t(e,b,i,j)
      listinfo(7,16)= nbblock
      listinfo(8,16)= noblock(spini)
      listinfo(9,16)= noblock(spinj)            
      listinfo(6,17)= 1                         ! 17.  t(e,a,k,i)
      listinfo(7,17)= nablock
      listinfo(8,17)= noblock(spink)
      listinfo(9,17)= noblock(spini)            
      listinfo(6,18)= 1                         ! 18.  t(e,b,k,i)
      listinfo(7,18)= nbblock
      listinfo(8,18)= noblock(spink)
      listinfo(9,18)= noblock(spini)            
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
      listinfo(6,19)= 1                         ! 19.  t(m,i,a,b)
      listinfo(7,19)= noblock(spini)
      listinfo(8,19)= nablock
      listinfo(9,19)= nbblock                   
      listinfo(6,20)= 1                         ! 20.  t(m,k,a,c)
      listinfo(7,20)= noblock(spink)
      listinfo(8,20)= nablock
      listinfo(9,20)= ncblock                   
      listinfo(6,21)= 1                         ! 21.  t(m,k,b,c)
      listinfo(7,21)= noblock(spink)
      listinfo(8,21)= nbblock
      listinfo(9,21)= ncblock                   
      listinfo(6,22)= 1                         ! 22.  t(m,i,c,a)
      listinfo(7,22)= noblock(spini)
      listinfo(8,22)= ncblock
      listinfo(9,22)= nablock                   
      listinfo(6,23)= 1                         ! 23.  t(m,i,c,b)
      listinfo(7,23)= noblock(spini)
      listinfo(8,23)= ncblock
      listinfo(9,23)= nbblock                   
c
c     Allocate and then create the pointer/reclen arrays
c
      do list = 1, 23
         listinfo(1,list) = listinfo(6,list)*listinfo(7,list)*
     $        listinfo(8,list)*listinfo(9,list)
         if (.not. ma_push_get(mt_int,2*listinfo(1,list),listname(list), 
     $        listinfo(3,list), listinfo(2,list))) call errquit
     $        ('uccsdt_triples_set_super_block: ma of list?', list)
      end do
c
c     1.  <ei||ab>  spin(e)=spin(a)  pure
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     noblock(spini), oblock(1,1,spini), nablock, ablock, 
     $     nbblock, bblock, int_mb(listinfo(2,1)), listinfo(4,1))
c
c     2.  <ek|ac>   spin(e)=spin(a)  mixed 
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nablock, ablock, 
     $     ncblock, cblock, int_mb(listinfo(2,2)), listinfo(4,2))
c
c     3.  <ek|bc>   spin(e)=spin(a)  mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nbblock, bblock, 
     $     ncblock, cblock, int_mb(listinfo(2,3)), listinfo(4,3))
c      
c     4.  <ei|ca>   spin(e)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock, 
     $     nablock, ablock, int_mb(listinfo(2,4)), listinfo(4,4))
c
c     5.  <ei|cb>   spin(e)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock, 
     $     nbblock, bblock, int_mb(listinfo(2,5)), listinfo(4,5))
c
c     1-x integrals
c
c     6.  <mc|ik>  spin(m)=spin(a)  mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     ncblock, cblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,6)), 
     $     listinfo(4,6))
c
c     7.  <ma||ij> spin(m)=spin(a)  pure
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     nablock, ablock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,7)), 
     $     listinfo(4,7))
c
c     8.  <mb||ij> spin(m)=spin(a)  pure
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     nbblock, bblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,8)), 
     $     listinfo(4,8))
c
c     9.  <ma|ki>  spin(m)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     nablock, ablock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,9)), 
     $     listinfo(4,9))
c
c     10. <mb|ki>  spin(m)=spin(c)  mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     nbblock, bblock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,10)), 
     $     listinfo(4,10))
c
c     2-x integrals (for V)
c
c     11. <ab||ij> pure
c
      call uccsdt_2x_pointers(nablock, ablock, nbblock, bblock,
     $     noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,11)), 
     $     listinfo(4,11))
c
c     12. <ac|ik>  mixed
c
      call uccsdt_2x_pointers(nablock, ablock, ncblock, cblock,
     $     noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,12)), 
     $     listinfo(4,12))
c
c     13. <bc|ik>  mixed
c
      call uccsdt_2x_pointers(nbblock, bblock, ncblock, cblock,
     $     noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,13)), 
     $     listinfo(4,13))
c
c     Amplitudes multiplying 3-x integrals (note similarity to 1-x integrals)
c
c     14.  t(e,c,i,k) spin(e)=spin(a) mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     ncblock, cblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spink), oblock(1,1,spink), int_mb(listinfo(2,14)), 
     $     listinfo(4,14))
c
c     15.  t(e,a,i,j) spin(e)=spin(a) pure
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     nablock, ablock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,15)), 
     $     listinfo(4,15))
c      
c     16.  t(e,b,i,j) spin(e)=spin(a) pure
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spina), 
     $     nbblock, bblock, noblock(spini), oblock(1,1,spini), 
     $     noblock(spinj), oblock(1,1,spinj), int_mb(listinfo(2,16)), 
     $     listinfo(4,16))
c
c     17.  t(e,a,k,i) spin(e)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     nablock, ablock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,17)), 
     $     listinfo(4,17))
c
c     18.  t(e,b,k,i) spin(e)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(nv_sym(0,spinc), 
     $     nbblock, bblock, noblock(spink), oblock(1,1,spink), 
     $     noblock(spini), oblock(1,1,spini), int_mb(listinfo(2,18)), 
     $     listinfo(4,18))
c
c     Amplitudes multiplying 1-x integrals (note similarity to 3-x integrals)
c
c     19.  t(m,i,a,b) spin(m)=spin(a) pure
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     noblock(spini), oblock(1,1,spini), nablock, ablock,
     $     nbblock, bblock, int_mb(listinfo(2,19)), listinfo(4,19))
c
c     20.  t(m,k,a,c) spin(m)=spin(a) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nablock, ablock,
     $     ncblock, cblock, int_mb(listinfo(2,20)), listinfo(4,20))
c
c     21.  t(m,k,b,c) spin(m)=spin(a) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spina), 
     $     noblock(spink), oblock(1,1,spink), nbblock, bblock,
     $     ncblock, cblock, int_mb(listinfo(2,21)), listinfo(4,21))
c
c     22.  t(m,i,c,a) spin(m)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock,
     $     nablock, ablock, int_mb(listinfo(2,22)), listinfo(4,22))
c
c     23.  t(m,i,c,b) spin(m)=spin(c) mixed
c
      call uccsdt_1x_3x_pointers(no_sym(0,spinc), 
     $     noblock(spini), oblock(1,1,spini), ncblock, cblock,
     $     nbblock, bblock, int_mb(listinfo(2,23)), listinfo(4,23))
c
c     Print a summary of the lists
c
      if (oprint) then
         write(6,2)
 2       format('  Summary of block-sorted list of integrals',
     $        ' and amplitudes',//,
     $        '  list    name    #blocks  #elements   block',
     $        ' dimensions',/
     $        ' ----- ---------- ------- ----------',
     $        ' --------------------')
         do list = 1, 23
            write(6,3) list, listname(list), listinfo(1,list), 
     $           listinfo(4,list), 
     $           listinfo(6,list), listinfo(7,list), listinfo(8,list), 
     $           listinfo(9,list)
 3          format(1x,i5,1x,a10,i8,i11,1x,4i5)
         end do
      endif
c     
c     Finally, we can allocate the arrays for the lists themselves
c
      do list = 1, 23
         if (.not. ga_create(mt_dbl, listinfo(4,list), 1, 
     $        listname(list), 0, 0, listinfo(5,list))) call errquit
     $        ('uccsdt_triples_set_super_block: ga_create', list)
         call ga_zero(listinfo(5,list))
      end do
c
      end
      subroutine uccsdt_1x_3x_pointers(
     $     ne_sym,
     $     niblock, iblock,
     $     nablock, ablock, 
     $     nbblock, bblock,
     $     pointer, length)
      implicit none
      integer ne_sym(0:7)       ! [input] No. of e orbitals by symmetry
      integer niblock, iblock(3,*) ! [input] Blocking for i
      integer nablock, ablock(3,*) ! [input] Blocking for a
      integer nbblock, bblock(3,*) ! [input] Blocking for b
      integer pointer(2,niblock, nablock, nbblock) ! [output]
      integer length            ! [output] Total length of the array
c
c     Generate the pointer to the start of a block of <ei|ab>
c     for all e (of appropriate symmetry) and i, a, b blocked
c     as specified
c
c     pointer(iblock,ablock,bblock) = start of block of <ei|ab>
c
c     Eventually the lists must be stored with randomization to avoid
c     hot spts.  Also, eventually the lists must store double precision
c     pointers to use the long GA interface for > 16 GB arrays.
c
      integer ptr, reclen
      integer a, b, i
      integer syma, symb, symi, syme
      integer dima, dimb, dimi, dime
c
      ptr = 1
      do b = 1, nbblock
         dimb = bblock(2,b) - bblock(1,b) + 1
         symb = bblock(3,b)
         do a = 1, nablock
            dima = ablock(2,a) - ablock(1,a) + 1
            syma = ablock(3,a)
            do i = 1, niblock
               dimi = iblock(2,i) - iblock(1,i) + 1
               symi = iblock(3,i)
               syme = ieor(symi,ieor(syma,symb))
               dime = ne_sym(syme)
c
               reclen = dime*dimi*dima*dimb
               pointer(1,i,a,b) = ptr
               pointer(2,i,a,b) = reclen
               ptr = ptr + reclen
            end do
         end do
      end do
      length = ptr - 1
c
      end
      subroutine uccsdt_2x_pointers(
     $     nablock, ablock,
     $     nbblock, bblock,
     $     niblock, iblock, 
     $     njblock, jblock,
     $     pointer, length)
      implicit none
      integer nablock, ablock(3,*) ! [input] Blocking for a
      integer nbblock, bblock(3,*) ! [input] Blocking for b
      integer niblock, iblock(3,*) ! [input] Blocking for i
      integer njblock, jblock(3,*) ! [input] Blocking for i
      integer pointer(2,nablock, nbblock, niblock, njblock) ! [output]
      integer length            ! [output] Total length of the array
c
c     Generate the pointer to the start of a block of <ab|ij>
c     for i, j, a, b blocked as specified
c
c     pointer(ablock,bblock,iblock,jblock) = start of block of <ab|ij>
c
c     Eventually the lists must be stored with randomization to avoid
c     hot spts.  Also, eventually the lists must store double precision
c     pointers to use the long GA interface for > 16 GB arrays.
c
c     THESE POINTER LISTS GET TOO BIG ... NEED TO REDUCE TO JUST
C     LAST 3 BLOCK DIMENSIONS AND USE LOGIC TO FIGURE OUT FIRST.
c
      integer ptr, reclen
      integer a, b, i, j
      integer syma, symb, symi, symj
      integer dima, dimb, dimi, dimj
c
      ptr = 1
      do j = 1, njblock
         dimj = jblock(2,j) - iblock(1,j) + 1
         symj = jblock(3,j)
         do i = 1, niblock
            dimi = iblock(2,i) - iblock(1,i) + 1
            symi = iblock(3,i)
            do b = 1, nbblock
               dimb = bblock(2,b) - bblock(1,b) + 1
               symb = bblock(3,b)
               do a = 1, nablock
                  dima = ablock(2,a) - ablock(1,a) + 1
                  syma = ablock(3,a)
c
                  if (ieor(ieor(ieor(symi,symj),syma),symb) .eq. 0) then
                     reclen = dimi*dimj*dima*dimb
                     pointer(1,a,b,i,j) = ptr
                     pointer(2,a,b,i,j) = reclen
                     ptr = ptr + reclen
                  else
                     pointer(1,a,b,i,j) = -999999999 ! Zero by symmetry
                     pointer(2,a,b,i,j) = -999999999 ! Zero by symmetry
                  end if
               end do
            end do
         end do
      end do
      length = ptr - 1
c
      end
      subroutine uccsdt_triples_unset_super_block()
      implicit none
#include "cuccsdtP.fh"
#include "mafdecls.fh"
#include "global.fh"
c     
c     Clean up any MA/GA allocated by set_super_block()
c
      integer list
c
      do list = 23, 1, -1
         if (.not. ma_pop_stack(listinfo(3,list))) call errquit
     $        ('uccsdt_triples_unset_super_block: ma overwrite?',list)
         if (.not. ga_destroy(listinfo(5,list))) call errquit
     $        ('uccsdt_triples_unset_super_block: ga destroy?', list)
      end do
c
      end
      subroutine uccsdt_triples_moints(rtdb)
      implicit none
#include "cuccsdtP.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "rtdb.fh"
      integer rtdb
c
c     Assumes that super- and sub-block info is already set
c     up in the common cuccsdtP.fh including the final global
c     arrays for block-sorted results.
c
      integer ushlo, ushhi, ush, vsh, ushblocksize, junk, list
      integer ublocklo, ublockhi, ublockdim
      parameter (ushblocksize=50) ! FOR DEBUG
      integer g_ai, g_bi, g_ck
      integer spini, spink, spine, spinm
      integer ilo, ihi, klo, khi, vlo, vhi, ulo, uhi, elo, ehi, mlo, mhi
      integer adim, bdim, cdim, idim, kdim, vdim, udim
      integer l_wivu, l_wkvu, l_work, l_buf, l_tmp, l_res
      integer k_wivu, k_wkvu, k_work, k_buf, k_tmp, k_res
      integer lenwork, lenbuf, dummy(3)
      integer me, next, nproc, count
      integer nxtval
      external nxtval
      double precision q2, tol
      logical odoit, oskel, ocheck
      parameter (tol = 1d-12, oskel = .true.)
      logical sym_shell_pair
      external sym_shell_pair
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      if (.not. rtdb_get(rtdb, 'uccsdt:ocheck', mt_log, 1, ocheck))
     $     ocheck = .false.
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
c     Need a more intelligent blocking of ush based upon
c     available memory, required parallelism, load-balance, ...
c     We want as large a block size as possible.
c
      do ushlo = 1, nsh, ushblocksize
         ushhi = min(ushlo+ushblocksize-1,nsh)
         if (.not. bas_cn2bfr(basis, ushlo, ublocklo, junk))
     $        call errquit('moints_shpair: basis?', 0)
         if (.not. bas_cn2bfr(basis, ushhi, junk, ublockhi))
     $        call errquit('moints_shpair: basis?', 0)
         ublockdim = ublockhi-ublocklo+1
c
c     First parallel step generates (u*|ai), (u*|bi) and (u*|ck)
c
         if (.not. ga_create(mt_dbl, ublockdim*nbf, adim*idim,
     $        '(u*|ai)', ublockdim*nbf, 0, g_ai)) call errquit
     $        ('moints_triples: failed allocaing (u*|ai)', 0)
         if (.not. ga_create(mt_dbl, ublockdim*nbf, bdim*idim,
     $        '(u*|bi)', ublockdim*nbf, 0, g_bi)) call errquit
     $        ('moints_triples: failed allocaing (u*|bi)', 0)
         if (.not. ga_create(mt_dbl, ublockdim*nbf, cdim*kdim,
     $        '(u*|ck)', ublockdim*nbf, 0, g_ck)) call errquit
     $        ('moints_triples: failed allocaing (u*|ck)', 0)
         call ga_zero(g_ai)
         call ga_zero(g_bi)
         call ga_zero(g_ck)
c
c     Next two loops in parallel with dynamic load balancing
c
         count = 0
         next = nxtval(nproc)
         do ush = ushlo, ushhi
            do vsh = 1, ush     ! NOTE RESTRICTION VSH <= USH
               odoit = schwarz_shell(ush,vsh)*schwarz_max().gt.tol
               if (odoit .and. oskel) 
     $              odoit = sym_shell_pair(basis, ush, vsh, q2)
               if (odoit) then
                  if (next .eq. count) then
                     if (.not. bas_cn2bfr(basis, ush, ulo, uhi))
     $                    call errquit('moints_shpair: basis?', 0)
                     if (.not. bas_cn2bfr(basis, vsh, vlo, vhi))
     $                    call errquit('moints_shpair: basis?', 0)
                     udim = uhi - ulo + 1
                     vdim = vhi - vlo + 1
c     
c     (wl|vu) -> (wi|vu), (wk|vu)
c     
                     if (.not. ma_push_get(mt_dbl, udim*vdim*nbf*idim,
     $                    '(wi|vu)', l_wivu, k_wivu)) call errquit
     $                    ('moints_triples: memory for (vu|wi)?', 
     $                    udim*vdim*nbf*idim)
                     if (.not. ma_push_get(mt_dbl, udim*vdim*nbf*kdim,
     $                    '(wk|vu)', l_wkvu, k_wkvu)) call errquit
     $                    ('moints_triples: memory for (vu|wk)?', 
     $                    udim*vdim*nbf*kdim)
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
c     Transform the second index and transpose on the fly
c     
                     lenwork = max(idim*adim,idim*bdim,kdim*cdim)
                     if (.not. ma_push_get(mt_dbl, lenwork, 'work',
     $                    l_work, k_work)) call errquit
     $                    ('moints_triples: memory for work', lenwork)
c     
c     (wi|vu) -> (ai|vu) -> (vu|ai)
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
                     if (.not. ma_pop_stack(l_work)) call errquit
     $                    ('moints_triples:popping work?', 0)
                     if (.not. ma_pop_stack(l_wkvu)) call errquit
     $                    ('moints_triples:popping wkvu?', 0)
                     if (.not. ma_pop_stack(l_wivu)) call errquit
     $                    ('moints_triples:popping wivu?', 0)
c     
                     next = nxtval(nproc)
                  end if
                  count = count + 1
               end if
            end do
         end do
c
c     Barrier
c
         call ga_sync()
         next = nxtval(-nproc)
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
     $        ('moints_triples: memory for buf', lenbuf)
         if (.not. ma_push_get(mt_dbl, lenbuf, 'tmp',
     $        l_tmp, k_tmp)) call errquit
     $        ('moints_triples: memory for temp', lenbuf)
         if (.not. ma_push_get(mt_dbl, lenbuf, 'res',
     $        l_res, k_res)) call errquit
     $        ('moints_triples: memory for temp', lenbuf)
                  
c
c     1.  <ei||ab>  = (ea|bi) - (eb|ai) = <14|23> - <14|32>
c     
         list = 1
         spine = asuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
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
c     2.  <ek|ac>   = (ea|ck) = <14|23>   [-<ek|ca> = -(ec|ai) = -<14|32>]
c
         list = 2
         spine = asuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
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
c     3.  <ek|bc>   = (eb|ck) = <14|23>  [-<ek|cb> = -(ec|bk) = -<14|32>]
c
         list = 3
         spine = asuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
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
c     4.  <ei|ca>   = (ec|ai) = <14|23>  [-<ei|ac> = -(ea|ci) = -<14|32>]
c
         list = 4
         spine = csuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
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
c     5.  <ei|cb>   = (ec|bi) = <14|23>  [-<ei|bc> = -(eb|ci) = -<14|32>]
c
         list = 5
         spine = csuper_spin
         elo = nc(spine) + no(spine) + 1
         ehi = nv(spine) + elo
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
c     6.  <mc|ik>   = (mi|ck) = <13|24>  [-<mc|ki> = -(mk|ci) = -<13|42>]
c
         list = 6
         spinm = spini
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
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
c     7.  <ma||ij>  = (mi|aj) - (mj|ai) = <13|24> - <13|42>
c
         list = 7
         spinm = spini
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
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
c     8.  <mb||ij>  = (mi|bj) - (mj|bi) = <13|24> - <13|42>
c
         list = 8
         spinm = spini
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
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
c     9.  <ma|ki>   = (mk|ai) = <13|24>  [-<ma|ik> = -(mi|ak) = -<13|42>]
c
         list = 9
         spinm = spink
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
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
c     10. <mb|ki>   = (mk|bi) = <13|24>  [-<mb|ik> = -(mi|bk) = -<13|42>]
c
         list = 10
         spinm = spink
         mlo = nc(spinm) + 1
         mhi = no(spinm) + mlo
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
c     11. <ab||ij>  = (ai|bj) - (aj|bi) = <13|24> - <13|42>
c
         list = 11
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
c     12. <ac|ik>   = (ai|ck) = <13|24>  [-<ac|ki> = -(ak|ci) = -<13|42>]
c
         list = 12
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
c     13. <bc|ik>   = (bi|ck) = <13|24>  [-<bc|ki> = -(bk|ci) = -<13|42>]
c
         list = 13
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
         if (.not. ma_pop_stack(l_res))
     $        call errquit('moints_triple: popping res?', 0)
         if (.not. ma_pop_stack(l_tmp))
     $        call errquit('moints_triple: popping tmp?', 0)
         if (.not. ma_pop_stack(l_buf))
     $        call errquit('moints_triple: popping buf?', 0)
c
         if (.not. ga_destroy(g_ck)) call errquit
     $        ('moints_triples: failed to destroy g_ck',0)
         if (.not. ga_destroy(g_bi)) call errquit
     $        ('moints_triples: failed to destroy g_bi',0)
         if (.not. ga_destroy(g_ai)) call errquit
     $        ('moints_triples: failed to destroy g_ai',0)
      end do
c
      if (ocheck) then
         do list = 1, 13
            call uccsdt_check_list(list)
         end do
      endif
c
*      call ma_summarize_allocated_blocks()
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
#include "global.fh"
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
     $     call errquit('moints_triples_final: neblock?',neblock)
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
         call errquit('moints_triples_final: bad order', 0)
      end if
c
      ublockdim = (ublockhi-ublocklo+1)
      edim = (ehi-elo+1)        ! Total no. of e-orbitals
      cdim = (chi-clo+1)        ! Total no. of c-orbitals
      adim = (ahi-alo+1)
      me = ga_nodeid()
      call ga_distribution(g_ck, me, vulo, vuhi, cklo, ckhi)
      if (vulo.ne.1 .and. vuhi.ne.ublockdim*nbf)
     $     call errquit('moints_triples_final: distribution?', 0)
c      
      do kblk = 1, nkblock
         nk = kblock(2,kblk) - kblock(1,kblk) + 1
         do k = kblock(1,kblk), kblock(2,kblk)
            do cblk = 1, ncblock
               nc = cblock(2,cblk) - cblock(1,cblk) + 1
               do c = cblock(1,cblk),cblock(2,cblk)
                  ck = 1 + c-clo + (k-klo)*cdim
                  if (ck.ge.cklo .and. ck.le.ckhi) then
                     call ga_get(g_ck, 1, ublockdim*nbf, ck, ck,
     $                    buf, 1)
c
c     C(v,e)*buf(v,u) -> temp(e,u)
c
                     call dgemm('t','n',edim,ublockdim,nbf,
     $                    1.0d0, evecs(1,elo), nbf, buf, nbf,
     $                    0.0d0, temp, edim)
c
c     phase*temp(e,u)*C(u,a) -> result(e,a)
c
                     call dgemm('n','n',edim,adim,ublockdim,
     $                    phase, temp, edim, avecs(ublocklo,alo), nbf,
     $                    0.0d0, result, edim)
c
c     We just did (vu) -> (eu) -> (ea), but since we are restricting
c     u>=v we also need to do (vu) -> (au) -> (ae).
c
                     call util_transpose(buf, nbf, temp, ublockdim,
     $                    nbf, ublockdim)
c
c     temp(u,v)*C(v,a) -> buf(u,a)
c
                     call dgemm('n', 'n', ublockdim, adim, nbf,
     $                    1.0d0, temp, ublockdim, avecs(1,alo), nbf,
     $                    0.0d0, buf, ublockdim)
c
c     phase*C(u,e)*buf(u,a) -> result(e,a)
c
                     call dgemm('t','n',edim,adim,ublockdim,
     $                    phase, evecs(ublocklo,elo),nbf,buf,ublockdim,
     $                    1.0d0, result, edim)
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
                                 call ga_acc(g_list, ptr, ptr+ne-1,
     $                                1, 1, result(ebase,a), 
     $                                1, 1.0d0)
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
#include "bas.fh"
#include "schwarz.fh"
#include "mafdecls.fh"
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
c
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
      integer nsh, nsh_block
      integer wshlo, wshhi, lshlo, lshhi
      integer wlo, whi, llo, lhi, junk, l_work, k_work, lenwork
c
      if (.not. bas_numcont(basis,nsh))
     $     call errquit('sym_abelian_bas_info: basis?', basis)
c
      if (nset.lt.1 .or. nset.gt.2) call errquit
     $     ('moints_shpair_one_index: bad nset', nset)
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
c     Loop over blocks of shells wsh & lsh in order to control memory
c     ... smaller blocks = less memory but less efficient integral
c     evaluation.  
c
c     Consider switching to texas 93 since this will give larger
c     blocks than 95.
c
      nsh_block = 30
      do wshlo = 1, nsh, nsh_block
         wshhi = min(wshlo+nsh_block-1,nsh)
         if (.not. bas_cn2bfr(basis, wshlo, wlo, junk))
     $        call errquit('moints_shpair: basis?', 0)
         if (.not. bas_cn2bfr(basis, wshhi, junk, whi))
     $        call errquit('moints_shpair: basis?', 0)
         do lshlo = 1, wshlo, nsh_block
            lshhi = min(lshlo+nsh_block-1,nsh)
            if (.not. bas_cn2bfr(basis, lshlo, llo, junk))
     $           call errquit('moints_shpair: basis?', 0)
            if (.not. bas_cn2bfr(basis, lshhi, junk, lhi))
     $           call errquit('moints_shpair: basis?', 0)
c
            lenwork = (lhi-llo+1)*(whi-wlo+1)*udim*vdim
            if (.not. ma_push_get(mt_dbl, lenwork, 'work', 
     $           l_work, k_work)) 
     $           call errquit('moints_shpair: ma?',lenwork)
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
            if (.not. ma_pop_stack(l_work)) 
     $           call errquit('moints_shpair: ma corrupt?',0)
         end do
      end do
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
c      
c     Form in ints(l,w,v,u) the integrals (lw|vu) with l<=w
c
      call moints_aoblock(
     $     basis, oskel,
     $     lshlo, lshhi, llo, lhi,
     $     wshlo, wshhi, wlo, whi,
     $     vsh, vsh, vlo, vhi,
     $     ush, ush, ulo, uhi,
     $     tol,
     $     ints)
c
c     Now have all of the AO integrals ... transform as necessary.
c     Since we are restricting w>=l we have to do l->i and w->i.
c
      idim1 = ihi1 - ilo1 + 1
      if (nset .eq. 2) idim2 = ihi2 - ilo2 + 1
      wdim = whi - wlo + 1
      ldim = lhi - llo + 1
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
      subroutine moints_aoblock(
     $     basis, oskel,
     $     ishlo, ishhi, ilo, ihi,
     $     jshlo, jshhi, jlo, jhi,
     $     kshlo, kshhi, klo, khi,
     $     lshlo, lshhi, llo, lhi,
     $     tol,
     $     ints)
      implicit none
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
      integer atom(nw_max_shells), geom, nsh
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
     $     ('moints_aoblock: bas_geom?', 0)
         if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('moints_aoblock: bas_numcont?', 0)
         do ish = 1, nsh
            if (.not. bas_cn2ce(basis, ish, atom(ish))) 
     $           call errquit('moints_aoblock: bas_cn2ce', ish)
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
     $     call errquit('moints_aoblock: ma?',maxq)
      if (.not. ma_push_get(mt_int, 4*maxq, 'quartets', l_q, k_q))
     $     call errquit('moints_aoblock: ma?', 4*maxq)
      if (.not. ma_push_get(mt_int, 4*leneri, 'labels', l_lab, k_lab))
     $     call errquit('moints_aoblock: ma?', 4*leneri)
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri))
     $     call errquit('moints_aoblock: ma?', leneri)
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr))
     $     call errquit('moints_aoblock: ma?', lenscr)
      k_iq = k_q          
      k_jq = k_q + maxq   
      k_kq = k_q + maxq*2 
      k_lq = k_q + maxq*3 
      k_ilab = k_lab       
      k_jlab = k_lab + leneri  
      k_klab = k_lab + leneri*2
      k_llab = l_lab + leneri*3
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
     $        call errquit('intb_init?',nq)
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
     $     ('moints_aoblock: ma corrupted?',0)
c
      end
