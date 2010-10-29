      logical function bas_vec_info(basis, info)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "inp.fh"
#include "geom.fh"
#include "stdio.fh"
      integer basis             ! [input] 
      character*16 info(*)      ! [output]
c
c     Return in info(i), i=1,nbf a simple character string
c     describing the i'th basis function
c     
c     info(i) = "atom number" "element symbol" "angular momentum" "component"
c     Note fixed declaration of info as character *16 to enable automatic
c     allocation in character *1 array from MA.
c
      integer nbf, nsh, ish, ilo, ihi, geom, i, iat, shdim, igen, ibf
      character*8 cbf
      logical status
      integer type, nprim, ngen, sphcart, atn, nat
      character symbol*2, tag*16, element*20
      character*1 shtypes(-1:7)
*
      character*8 cartbf(-1:7,1:36)
*
      data shtypes/'l','s', 'p', 'd', 'f', 'g', 'h', 'i', 'k'/
c
c
c     Sleazy initialization of cartesian bf info
c
      cartbf(-1,1) = 's'
      cartbf(-1,2) = 'x'
      cartbf(-1,3) = 'y'
      cartbf(-1,4) = 'z'

      cartbf(0, 1) = ' '
      
      cartbf(1, 1) = 'x'
      cartbf(1, 2) = 'y'
      cartbf(1, 3) = 'z'
      
      cartbf(2, 1) = 'xx'
      cartbf(2, 2) = 'xy'
      cartbf(2, 3) = 'xz'
      cartbf(2, 4) = 'yy'
      cartbf(2, 5) = 'yz'
      cartbf(2, 6) = 'zz'
      
      cartbf(3, 1) = 'xxx'
      cartbf(3, 2) = 'xxy'
      cartbf(3, 3) = 'xxz'
      cartbf(3, 4) = 'xyy'
      cartbf(3, 5) = 'xyz'
      cartbf(3, 6) = 'xzz'
      cartbf(3, 7) = 'yyy'
      cartbf(3, 8) = 'yyz'
      cartbf(3, 9) = 'yzz'
      cartbf(3,10) = 'zzz'
      
      cartbf(4, 1) = 'xxxx'
      cartbf(4, 2) = 'xxxy'
      cartbf(4, 3) = 'xxxz'
      cartbf(4, 4) = 'xxyy'
      cartbf(4, 5) = 'xxyz'
      cartbf(4, 6) = 'xxzz'
      cartbf(4, 7) = 'xyyy'
      cartbf(4, 8) = 'xyyz'
      cartbf(4, 9) = 'xyzz'
      cartbf(4,10) = 'xzzz'
      cartbf(4,11) = 'yyyy'
      cartbf(4,12) = 'yyyz'
      cartbf(4,13) = 'yyzz'
      cartbf(4,14) = 'yzzz'
      cartbf(4,15) = 'zzzz'
      
      cartbf(5, 1) = 'xxxxx'
      cartbf(5, 2) = 'xxxxy'
      cartbf(5, 3) = 'xxxxz'
      cartbf(5, 4) = 'xxxyy'
      cartbf(5, 5) = 'xxxyz'
      cartbf(5, 6) = 'xxxzz'
      cartbf(5, 7) = 'xxyyy'
      cartbf(5, 8) = 'xxyyz'
      cartbf(5, 9) = 'xxyzz'
      cartbf(5,10) = 'xxzzz'
      cartbf(5,11) = 'xyyyy'
      cartbf(5,12) = 'xyyyz'
      cartbf(5,13) = 'xyyzz'
      cartbf(5,14) = 'xyzzz'
      cartbf(5,15) = 'xzzzz'
      cartbf(5,16) = 'yyyyy'
      cartbf(5,17) = 'yyyyz'
      cartbf(5,18) = 'yyyzz'
      cartbf(5,19) = 'yyzzz'
      cartbf(5,20) = 'yzzzz'
      cartbf(5,21) = 'zzzzz'
      
      cartbf(6, 1) = 'xxxxxx'
      cartbf(6, 2) = 'xxxxxy'
      cartbf(6, 3) = 'xxxxxz'
      cartbf(6, 4) = 'xxxxyy'
      cartbf(6, 5) = 'xxxxyz'
      cartbf(6, 6) = 'xxxxzz'
      cartbf(6, 7) = 'xxxyyy'
      cartbf(6, 8) = 'xxxyyz'
      cartbf(6, 9) = 'xxxyzz'
      cartbf(6,10) = 'xxxzzz'
      cartbf(6,11) = 'xxyyyy'
      cartbf(6,12) = 'xxyyyz'
      cartbf(6,13) = 'xxyyzz'
      cartbf(6,14) = 'xxyzzz'
      cartbf(6,15) = 'xxzzzz'
      cartbf(6,16) = 'xyyyyy'
      cartbf(6,17) = 'xyyyyz'
      cartbf(6,18) = 'xyyyzz'
      cartbf(6,19) = 'xyyzzz'
      cartbf(6,20) = 'xyzzzz'
      cartbf(6,21) = 'xzzzzz'
      cartbf(6,22) = 'yyyyyy'
      cartbf(6,23) = 'yyyyyz'
      cartbf(6,24) = 'yyyyzz'
      cartbf(6,25) = 'yyyzzz'
      cartbf(6,26) = 'yyzzzz'
      cartbf(6,27) = 'yzzzzz'
      cartbf(6,28) = 'zzzzzz'
      
      cartbf(7, 1) = 'xxxxxxx'
      cartbf(7, 2) = 'xxxxxxy'
      cartbf(7, 3) = 'xxxxxxz'
      cartbf(7, 4) = 'xxxxxyy'
      cartbf(7, 5) = 'xxxxxyz'
      cartbf(7, 6) = 'xxxxxzz'
      cartbf(7, 7) = 'xxxxyyy'
      cartbf(7, 8) = 'xxxxyyz'
      cartbf(7, 9) = 'xxxxyzz'
      cartbf(7,10) = 'xxxxzzz'
      cartbf(7,11) = 'xxxyyyy'
      cartbf(7,12) = 'xxxyyyz'
      cartbf(7,13) = 'xxxyyzz'
      cartbf(7,14) = 'xxxyzzz'
      cartbf(7,15) = 'xxxzzzz'
      cartbf(7,16) = 'xxyyyyy'
      cartbf(7,17) = 'xxyyyyz'
      cartbf(7,18) = 'xxyyyzz'
      cartbf(7,19) = 'xxyyzzz'
      cartbf(7,20) = 'xxyzzzz'
      cartbf(7,21) = 'xxzzzzz'
      cartbf(7,22) = 'xyyyyyy'
      cartbf(7,23) = 'xyyyyyz'
      cartbf(7,24) = 'xyyyyzz'
      cartbf(7,25) = 'xyyyzzz'
      cartbf(7,26) = 'xyyzzzz'
      cartbf(7,27) = 'xyzzzzz'
      cartbf(7,28) = 'xzzzzzz'
      cartbf(7,29) = 'yyyyyyy'
      cartbf(7,30) = 'yyyyyyz'
      cartbf(7,31) = 'yyyyyzz'
      cartbf(7,32) = 'yyyyzzz'
      cartbf(7,33) = 'yyyzzzz'
      cartbf(7,34) = 'yyzzzzz'
      cartbf(7,35) = 'yzzzzzz'
      cartbf(7,36) = 'zzzzzzz'
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('bas_vec_info: nbf?',0, BASIS_ERR)
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('bas_vec_info: nsh?',0, BASIS_ERR)
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('bas_vec_info: geom?', 0, BASIS_ERR)
      if (.not. geom_ncent(geom, nat)) call errquit
     $     ('bas_vec_info: ncent?',0, BASIS_ERR)
c
      i = 0
      do ish = 1, nsh
         if (.not. bas_continfo(basis, ish, type, nprim, ngen, sphcart))
     $        call errquit('bas_vec_info: continfo',ish, BASIS_ERR)
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('bas_vec_info: cn2bfr', ish, BASIS_ERR)
         if (.not. bas_cn2ce(basis, ish, iat))
     $        call errquit('bas_vec_info: cn2ce', ish, BASIS_ERR)
         if (.not. geom_cent_tag(geom, iat, tag))
     $        call errquit('bas_vec_info: tag', iat, BASIS_ERR)
         status = geom_tag_to_element(tag, symbol, element, atn)
c
         if (type.gt.7 .or. type.lt.-1) call errquit
     $        ('bas_vec_info: unsupported basis type?',type,
     &       BASIS_ERR)
c
         shdim = (ihi-ilo+1)
         if (type.ge.0) then 
            shdim = shdim / ngen
         else     ! Fudge for SP shells
            ngen = 1
         endif
         if (type .le. 1) sphcart = 0
c         
         do igen = 1, ngen
            do ibf = 1, shdim
               if (sphcart .eq. 0) then
                  cbf = cartbf(type,ibf)
               else
                  cbf = ' '
                  write(cbf,'(i3)') ibf-shdim/2-1
               endif
c
               i = i + 1
               info(i) = ' '
c
c     easier than writing the format to write the format
c
               if (nat.lt.10) then
                  write(info(i), '(i1,1x,a2,1x,a1,a7)')
     $                 iat, symbol, shtypes(type), cbf(1:7)
               else if (nat.lt.100) then
                  write(info(i), '(i2,1x,a2,1x,a1,a7)')
     $                 iat, symbol, shtypes(type), cbf(1:7)
               else if (nat.lt.1000) then
                  write(info(i), '(i3,1x,a2,1x,a1,a7)')
     $                 iat, symbol, shtypes(type), cbf(1:7)
               else
                  write(info(i), '(i4,1x,a2,1x,a1,a7)')
     $                 iat, symbol, shtypes(type), cbf(1:7)
               endif
c
*               write(LuOut,*) i, ' "', info(i), '"'
c
            enddo
         enddo
      enddo
c
      if (i.ne.nbf) call errquit('bas_vec_info: function mismatch?',i,
     &       BASIS_ERR)
c
      bas_vec_info = .true.
c
      end
c$$$      subroutine bas_test_vec_info(basis)
c$$$      implicit none
c$$$#include "errquit.fh"
c$$$      integer basis
c$$$      character*16 info(500)
c$$$      logical bas_vec_info
c$$$      external bas_vec_info
c$$$c
c$$$      if (.not. bas_vec_info(basis, info)) call errquit('u',0)
c$$$c
c$$$      end
      subroutine bas_print_labels(basis)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "stdio.fh"
      logical  bas_vec_info
      external bas_vec_info
      integer basis
c
      integer i, j, k_bftags, l_bftags, nbf
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('bas_print_labels: invalid basis handle',basis,
     &       BASIS_ERR)
      if (.not. ma_push_get(mt_byte,nbf*16,'bftags',
     $     l_bftags,k_bftags))
     $     call errquit('bas_print_labels: ma failed', nbf*16,
     &       BASIS_ERR)
      if (.not.bas_vec_info(basis, byte_mb(k_bftags)))
     &    call errquit('bas_print_labels: bas_vec_info failed',911,
     &       BASIS_ERR)
      write(LuOut,*)
      call util_print_centered(LuOut, 'Basis function labels',
     $     13, .true.)
      write(LuOut,780)
 780  format(4x,'   No.  Atm El Type'/
     $       4x,'   ---  --- -- ----')
      do i = 1, nbf
         write(LuOut,781) i, (byte_mb(k_bftags+(i-1)*16-1+j),j=1,16)
 781     format(4x,i5,4x,16a1)
      enddo
      if (.not. ma_pop_stack(l_bftags)) call errquit
     $     ('bas_print_labels: ma corrupt after bas func labs',0,
     &       BASIS_ERR)
c
      end
