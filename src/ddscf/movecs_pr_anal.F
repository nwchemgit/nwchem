      subroutine movecs_print_anal(basis, ilo, ihi, thresh,
     $     g_vecs, title, 
     $     oevals, evals, oirs, irs, oocc, occ)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
c
      integer basis
      integer ilo, ihi          ! [input] Range of vectors to print
      double precision thresh   ! [input] Print coeffs with absval >= thresh
      integer g_vecs
      character*(*) title
      logical oevals            ! [input] If true print eigenvalues
      double precision evals(*)
      logical oirs              ! [input] If true print irreps
      integer irs(*)
      logical oocc              ! [input] If true print occupations
      double precision occ(*)
c
c     Print a summary of the MO vectors in the specified range.
c
      integer l_vecs, k_vecs, i, j, k_tags, l_tags, k_list, l_list
      integer n, i, j, k, m, klo, khi, nbf, ibuf
      character buf*80
c  
      integer maxop, maxireps, geom
      parameter (maxop = 120, maxireps=20)
      integer nop, nir,  class_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      character*20 zname
      double precision chars(maxireps*maxireps)
      logical sym_char_table
      external sym_char_table
c
      call ga_sync
c
      if (ga_nodeid() .eq. 0) then
         if (.not. bas_numbf(basis, nbf)) call errquit
     $        ('movecs_anal: basis bad?',basis)
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_anal1',l_vecs,k_vecs))
     $        call errquit('movecs_anal: ma 1 failed', nbf)
         if (.not. ma_push_get(mt_int,nbf,'movecs_anal2',l_list,k_list))
     $        call errquit('movecs_anal: ma 2 failed', nbf)
         if (.not. ma_push_get(mt_byte,nbf*16,'movecs_anal3',
     $        l_tags,k_tags))
     $        call errquit('movecs_anal: ma 3 failed', nbf*16)
c
         if (oirs) then
            if (.not. bas_geom(basis, geom)) call errquit
     $           ('movecs_anal: bas geom', basis)
            call sym_group_name(geom, zname)
            if (.not. sym_char_table(zname, nop, nir, class_dim,
     $           zir, zclass, chars)) call errquit
     $           ('movecs_anal: no char table available ',geom)
         endif
c     
         call bas_vec_info(basis, byte_mb(k_tags))
c
         write(6,*)
         call util_print_centered(6,title, 40, .true.)
         write(6,*)
c

         do i = ilo, ihi
c
            call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), 1)
c
c     Identify significant coefficients and sort by size
c
            n = 0
            do j = 0, nbf-1
               if (abs(dbl_mb(k_vecs+j)).ge.thresh) then
                  int_mb(k_list + n) = j
                  n = n + 1
               endif
            enddo
            do j = 0, n-1
               do k = 0, j
                  if (abs(dbl_mb(k_vecs+int_mb(k_list+k))).lt.
     $                abs(dbl_mb(k_vecs+int_mb(k_list+j)))) then
                     m = int_mb(k_list+j)
                     int_mb(k_list+j) = int_mb(k_list+k)
                     int_mb(k_list+k) = m
                  endif
               enddo
            enddo
c
c     Construct optional output line
c
            ibuf = 1
            buf = ' '
            if (oocc) then
               write(buf(ibuf:),'(''Occ='',1p,d12.6)'), occ(i)
               ibuf = ibuf + 18
            endif
            if (oevals) then
               write(buf(ibuf:),'(''E='',1p,d13.6)'), evals(i)
               ibuf = ibuf + 17
            endif
            if (oirs) then
               write(buf(ibuf:),'(''Symmetry='',a4)'), zir(irs(i))
               ibuf = ibuf + 18
            endif
            write(6,1) i, buf(1:max(inp_strlen(buf),1))
 1          format(' Vector',i5,2x,a) 
c
c     Output the analysis
c
            write(6,22)
 22         format(1x,2('  Bfn.  Coefficient  Atom+Function  ',5x))
            write(6,23)
 23         format(1x,2(' ----- ------------  ---------------',5x))
            do klo = 0, n-1, 2
               khi = min(klo+1,n-1)
               write(6,2) (
     $              int_mb(k_list+k)+1, 
     $              dbl_mb(k_vecs+int_mb(k_list+k)),
     $              (byte_mb(k_tags+int_mb(k_list+k)*16+m),m=0,15), 
     $              k = klo,khi)
 2             format(1x,2(i5,2x,f12.6,2x,16a1,4x))
            enddo
            write(6,*)
         enddo
         call util_flush(6)
         if (.not. ma_chop_stack(l_vecs)) call errquit
     $        ('bas_vec_info: ma pop?', 0)
      endif
c
      call ga_sync()
c
      end


