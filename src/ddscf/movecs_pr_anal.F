      subroutine movecs_print_anal(basis, ilo, ihi, thresh,
     $     g_vecs, title, 
     $     oevals, evals, oirs, irs, oocc, occ)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "stdio.fh"
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
      logical or2               ! [input] If true print orbital center and r^2
c
c     Print a summary of the MO vectors in the specified range.
c
      integer l_vecs, k_vecs, i, j, k_tags, l_tags, k_list, l_list
      integer n, k, m, klo, khi, ibuf
      character buf*80
      double precision cur_thresh
c  
      integer type, nbf, nmo
      integer len_r_r2, l_r_r2, k_r_r2
      integer maxop, maxireps, geom
      parameter (maxop = 120, maxireps=20)
      integer nop, nir,  class_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps)
      character*20 zname
      double precision chars(maxireps*maxireps)
      logical sym_char_table
      external sym_char_table
      double precision maxv
      integer howmany9
c
      or2 = .true.
      if (.not. bas_cando_mpoles(basis)) or2 = .false.
c     
      call ga_sync
      if (oscfps) call pstat_on(ps_moanal)
      call ga_inquire(g_vecs, type, nbf, nmo)
      if (or2) then
c
c        local array for x, y, z, x^2, y^2, z^2 and r^2 for each MO
c
         len_r_r2 = 7*nmo
         if (.not. ma_push_get(mt_dbl, len_r_r2, 'nmo:r_and_r2', 
     &      l_r_r2, k_r_r2))call errquit
     &      ('movecs_pr_anal: cannot allocate r_r2', len_r_r2,
     &       MA_ERR)
         call mo_r_and_r2(basis, g_vecs, dbl_mb(k_r_r2))
      endif
c
      if (ga_nodeid() .eq. 0) then
         if (.not. bas_numbf(basis, nbf)) call errquit
     $        ('movecs_anal: basis bad?',basis, BASIS_ERR)
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_anal1',l_vecs,k_vecs))
     $        call errquit('movecs_anal: ma 1 failed', nbf, MA_ERR)
         if (.not. ma_push_get(mt_int,nbf,'movecs_anal2',l_list,k_list))
     $        call errquit('movecs_anal: ma 2 failed', nbf, MA_ERR)
         if (.not. ma_push_get(mt_byte,nbf*16,'movecs_anal3',
     $        l_tags,k_tags))
     $        call errquit('movecs_anal: ma 3 failed', nbf*16, MA_ERR)
c
         if (oirs) then
            if (.not. bas_geom(basis, geom)) call errquit
     $           ('movecs_anal: bas geom', basis, BASIS_ERR)
            call sym_group_name(geom, zname)
            if (.not. sym_char_table(zname, nop, nir, class_dim,
     $           zir, zclass, chars)) call errquit
     $           ('movecs_anal: no char table available ',geom,
     &       UNKNOWN_ERR)
         endif
c     
         call bas_vec_info(basis, byte_mb(k_tags))
c 
c if there is no title sent, don't print it.
         if(title.eq.' ') goto 110
         write(LuOut,*)
         call util_print_centered(LuOut,title, 40, .true.)
         write(LuOut,*)
c
 110     continue
         howmany9=9
c     if thresh small, we want lots of eigenv
         if(thresh.lt.1d-3) howmany9=99
         do i = ilo, ihi
c
            call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), nbf)
c
c     Identify significant coefficients and sort by size
c
            n = 0
            cur_thresh = thresh
 111        do j = 0, nbf-1
               if (abs(dbl_mb(k_vecs+j)).ge.cur_thresh) then
                  int_mb(k_list + n) = j
                  n = n + 1
               endif
            enddo
            if (n.eq.0) then
c
c     found nothing. get max +10% of lowe
c
               maxv=1d-99
c     max
               do j = 0, nbf-1
                  if (abs(dbl_mb(k_vecs+j)).ge.maxv) then
                     int_mb(k_list) = j
                     maxv=abs(dbl_mb(k_vecs+j))
                  endif
               enddo
c     get within 10% range
               maxv=0.9d0*maxv
               do j = 0, nbf-1
                  if (abs(dbl_mb(k_vecs+j)).ge.maxv) then
                     int_mb(k_list + n) = j
                     n=n+1
                  endif
               enddo
            endif
            do j = 0, n-1
               do k = 0, j
                  if ((abs(dbl_mb(k_vecs+int_mb(k_list+k)))*0.995.lt.
     $                 abs(dbl_mb(k_vecs+int_mb(k_list+j)))).and.
     $                (abs(dbl_mb(k_vecs+int_mb(k_list+k)))*1.005.lt.
     $                 abs(dbl_mb(k_vecs+int_mb(k_list+j))))) then
                     m = int_mb(k_list+j)
                     int_mb(k_list+j) = int_mb(k_list+k)
                     int_mb(k_list+k) = m
                  elseif((abs(dbl_mb(k_vecs+int_mb(k_list+k)))*0.995.lt.
     $                    abs(dbl_mb(k_vecs+int_mb(k_list+j))))
     &                   .and..not.
     $                   (abs(dbl_mb(k_vecs+int_mb(k_list+k)))*1.005.lt.
     $                    abs(dbl_mb(k_vecs+int_mb(k_list+j))))) then
                     if (int_mb(k_list+j).lt.int_mb(k_list+k)) then
                       m = int_mb(k_list+j)
                       int_mb(k_list+j) = int_mb(k_list+k)
                       int_mb(k_list+k) = m
                     endif
                  else
c                    the ordering is OK so do nothing
                  endif
               enddo
            enddo
c
c     Construct optional output line
c
            ibuf = 1
            buf = ' '
            if (oocc) then
               write(buf(ibuf:),'(''Occ='',1p,d12.6)') occ(i)
               ibuf = ibuf + 18
            endif
            if (oevals) then
               write(buf(ibuf:),'(''E='',1p,d13.6)') evals(i)
               ibuf = ibuf + 17
            endif
            if (oirs) then
               write(buf(ibuf:),'(''Symmetry='',a4)') zir(irs(i))
               ibuf = ibuf + 18
            endif
            write(LuOut,1) i, buf(1:max(inp_strlen(buf),1))
 1          format(' Vector',i5,2x,a) 
            if (or2) then
c
c              Construct optional 2nd output line
c
               ibuf = 1
               buf = ' '
               write(buf(ibuf:),'(''MO Center='',1p,3(1x,d8.1,'',''))') 
     &         dbl_mb(k_r_r2+(i-1)*7),
     &         dbl_mb(k_r_r2+(i-1)*7+1),
     &         dbl_mb(k_r_r2+(i-1)*7+2)
               ibuf = ibuf + 41
               write(buf(ibuf:),'(''r^2='',1p,d8.1)')
     &         dbl_mb(k_r_r2+(i-1)*7+6)
               ibuf = ibuf + 14
               write(LuOut,3) buf(1:max(inp_strlen(buf),1))
 3             format('       ',7x,a) 
            endif
c
c     Output the analysis
c
            write(LuOut,22)
 22         format(1x,2('  Bfn.  Coefficient  Atom+Function  ',5x))
            write(LuOut,23)
 23         format(1x,2(' ----- ------------  ---------------',5x))
            do klo = 0, min(n-1,howmany9), 2
               khi = min(klo+1,n-1)
               write(LuOut,2) (
     $              int_mb(k_list+k)+1, 
     $              dbl_mb(k_vecs+int_mb(k_list+k)),
     $              (byte_mb(k_tags+int_mb(k_list+k)*16+m),m=0,15), 
     $              k = klo,khi)
 2             format(1x,2(i5,2x,f12.6,2x,16a1,4x))
            enddo
            write(LuOut,*)
         enddo
         call util_flush(LuOut)
         if (.not. ma_chop_stack(l_vecs)) call errquit
     $        ('bas_vec_info: ma pop?', 0, BASIS_ERR)
      endif
      if (or2) then
         if (.not. ma_pop_stack(l_r_r2)) call errquit('movecs_anal?',0,
     &       MA_ERR)
      endif
c
      if (oscfps) call pstat_off(ps_moanal)
      call ga_sync()
c
      end
      subroutine mo_r_and_r2(basis, g_movecs, r_and_r2)
c
      implicit none
c
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer basis                    ! [input] basis 
      double precision r_and_r2(7,*) ! [output] x,y,z,x^2,y^2,z^2,r^2 for each mo
      integer g_movecs                 ! [input] GA mo vectors
c
      double precision center(3)
      integer i, j
      integer type, nbf, nmo
      integer g_xlm
      data center/3*0.0d0/
c
c     compute r and r^2 for each MO at the origin
c
      call ga_inquire(g_movecs, type, nbf, nmo)
c
c     create a global array to store x, y, z and x^2, y^2, z^2 for each AO
c
      if (.not. ga_create(mt_dbl, 6*nbf, nbf, 'GXLM',
c     $     32,32,g_xlm)) call errquit('mo_r_and_r2: g_xlm',6*nbf*nbf,
     $     nbf,32,g_xlm)) call errquit('mo_r_and_r2: g_xlm',6*nbf*nbf,
     &       GA_ERR)
c
c     compute dipoles and quadrupole components for each AO
c
      call xlm_make_ao_poles(basis, center, g_xlm)
c
c     compute dipoles and quadrupole components for each MO
c
      call xlm_ao_to_mo(g_movecs, g_xlm)
c
c     put global data into local array on node 0 (node 0
c     will be doing the printing of the desired info)
c
      call dfill(nmo*7, 0.0d0, r_and_r2, 1)
      if (ga_nodeid() .eq. 0) then
         do i = 1, nmo
            call ga_get(g_xlm,(i-1)*6+1,(i-1)*6+6,i,i,
     $           r_and_r2(1,i), 6)
         enddo
c
c         write(LuOut,*)' x, y, z, x^2, y^2, z^2 at (0,0,0)'
c         call output(r_and_r2, 1, 7, 1, nmo, 7, nmo, 1)
c
c        shift and sum
c
         do i = 1, nmo
           do j = 1, 3
             r_and_r2(j+3,i) = r_and_r2(j+3,i) - r_and_r2(j,i)**2
           enddo
           r_and_r2(7,i) = r_and_r2(4,i) + r_and_r2(5,i) + r_and_r2(6,i)
         enddo
c         write(LuOut,*)' x, y, z, x^2, y^2, z^2, r^2 shifted'
c         call output(r_and_r2, 1, 7, 1, nmo, 7, nmo, 1)
      endif
c      
      if (.not. ga_destroy(g_xlm)) call errquit('mo_r_and_r2: ga?',0,
     &       GA_ERR)
c
      return
      end
      subroutine xlm_make_ao_poles(basis, center, g_xlm)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
#include "util_params.fh"
#include "util.fh"
c     
      integer basis             ! [input] basis 
      double precision center(3) ! [input] the expansion center
      integer g_xlm              ! [input] GA that will return the mpoles
c
      double precision one, two
      parameter (one=1.d0, two=2.d0)
      double precision autoang2
      parameter (autoang2 = cau2ang*cau2ang)
c
      integer mcart, l_xlm, k_xlm
      integer geom
      integer nshell, noperators, maxscr, me, nproc
      integer nbf_max, lmpmax, ishell, jshell, ijshell
      integer ilo, ihi, jlo, jhi, idim, jdim, ind, i, j, l, ioff
      integer l_scr, k_scr, l_mp, k_mp
      integer lmax              ! Maximum value of L = 2
      integer clo, chi, rlo, rhi
      integer type, dim1, dim2
      logical doloop, oprint
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('multipole: bad basis', 0, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $     ('xlm_pole: bas_numcont failed for basis', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit(
     &     'xlm_pole: bas_nbf_cn_max failed',20, BASIS_ERR)
c
c
c     note lmax is hardwired to 2
c
      lmax = 2
c     
c     length of int_mpole integral output for full square list
c     includes l_pole = 0,...,lmax, where l_pole = 0 is simply
c     the 2-c overlap matrix.  (cartesian or sphericalcomponents).
c     
      noperators = (lmax+1)*(lmax+2)*(lmax+3)/6
      call int_mem_dipole(lmpmax,maxscr,basis,basis,lmax)
      maxscr = max(100000,maxscr)
c     
c     allocate necessary local temporary arrays on the stack
c     
      if(.not. ma_push_get(mt_dbl, lmpmax, 'mult:mp', l_mp, k_mp))
     &     call errquit('xlm_pole: cannot allocate mp', lmpmax, MA_ERR)
      if(.not. ma_push_get(mt_dbl, lmpmax, 'mult:xlm', l_xlm, k_xlm))
     &     call errquit('xlm_pole: cannot allocate xlm', lmpmax, MA_ERR)
      if(.not. ma_push_get(mt_dbl, maxscr, 'mult:scr', l_scr, k_scr))
     &     call errquit('xlm_pole: cannot allocate scratch', maxscr,
     &       MA_ERR)
c     
      call ga_distribution(g_xlm,
     .     ga_nodeid(), rlo, rhi, clo, chi)
      doloop=.not.(rlo.eq.0.and.rhi.eq.-1.and.clo.eq.0.and.chi.eq.-1)
      call ga_inquire(g_xlm, type, dim1, dim2)
      oprint=(ga_nodeid() .eq. 0 .and. 
     $     util_print('debug_xlm_make', print_debug))
      if(oprint) then
         write(luout,*) ' dim1 dim2 ',dim1,dim2
         write(luout,*) ' rlo rhi ',rlo,rhi
         write(luout,*) ' clo chi ',clo,chi
      endif
      call ga_zero(g_xlm)
       nproc=0
       if(doloop) nproc=1
       call ga_igop(65536,nproc,1,'+')
       if(oprint) then
          write(luout,*) ' using ',nproc,
     C      ' nodes out of ',ga_nnodes()
          call util_flush(6)
       endif
c
      if (doloop) then
      ijshell = -1
      me = ga_nodeid()
      do ishell = 1, nshell
         if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &        ('xlm_pole: bas_cn2bfr failed for basis', basis,
     &       BASIS_ERR)
         idim = ihi - ilo + 1

         do jshell = 1, nshell
            ijshell = ijshell + 1
            if (mod(ijshell,nproc) .eq. me) then
               if (.not. bas_cn2bfr(basis, jshell, jlo, jhi))
     &              call errquit('xlm_pole: bas_cn2bfr', basis,
     &       BASIS_ERR)
               jdim = jhi - jlo + 1
c     
               call int_mpole(basis, ishell, basis, jshell,
     &              lmax, center, maxscr, dbl_mb(k_scr),
     &              lmpmax, dbl_mb(k_mp))
c     
c     output from int_mpole is: overlap, dipole, q-pole, ...
c     within a multipole block, the order is <i|m|j>  j fastest, 
c     then m, then i ... we must put m first 
c     
               call dfill(6*idim*jdim,0.0d0,dbl_mb(k_xlm),1)
c
c               write(LuOut,*)' ishell, jshell = ', ishell, jshell
               ind = k_mp
               do l = 0, lmax
                  do i = 1, idim
                     do mcart = 1, ((l+1)*(l+2))/2
                        do j = 1, jdim
                           ioff = k_xlm + 6*(j-1 + jdim*(i-1))
c       write(LuOut,*)' l, i, mcart, j, ind-k_mp, dbl_mb(ind) ',
c     &             l, i, mcart, j, ind-k_mp, dbl_mb(ind)
                           if (l.eq.1.and.mcart.eq.1)then
                              dbl_mb(ioff) = dbl_mb(ind)*cau2ang
c                              write(LuOut,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.1.and.mcart.eq.2)then
                              ioff = ioff + 1
                              dbl_mb(ioff) = dbl_mb(ind)*cau2ang
c                              write(LuOut,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.1.and.mcart.eq.3)then
                              ioff = ioff + 2
                              dbl_mb(ioff) = dbl_mb(ind)*cau2ang
c                              write(LuOut,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.2.and.mcart.eq.1)then
                              ioff = ioff + 3
                              dbl_mb(ioff) = dbl_mb(ind)*autoang2
c                              write(LuOut,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.2.and.mcart.eq.4)then
                              ioff = ioff + 4
                              dbl_mb(ioff) = dbl_mb(ind)*autoang2
c                              write(LuOut,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.2.and.mcart.eq.6)then
                              ioff = ioff + 5
                              dbl_mb(ioff) = dbl_mb(ind)*autoang2
c                              write(LuOut,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           ind = ind + 1
                        end do
                     end do
                  end do
               end do
c
               call ga_put(g_xlm,1+(jlo-1)*6,jhi*6,ilo,ihi,
     $              dbl_mb(k_xlm), 6*jdim)
c
            end if
         end do
      end do
c
      endif
      call ga_sync
c
c      write(LuOut,*) ' THE AO MPOLES '
c      call ga_print(g_xlm)
c     
c     clean up stack
c     
      if (.not. ma_pop_stack(l_scr)) call errquit('xlm_pole: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_xlm)) call errquit('xlm_pole: ma?',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_mp)) call errquit('xlm_pole: ma?',0,
     &       MA_ERR)
c
      end
      subroutine xlm_ao_to_mo(g_vecs, g_xlm)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer lmax, g_vecs, g_xlm
c
c     Transform the multipoles from AO to MO basis overwriting
c     the input AO set.
c
      integer type, nbf, nmo, g_tmp1, g_tmp2, k_tmp, l_tmp, 
     $     i, lm
c
c     note lmax is hardwired to 2
c
      lmax = 2
c
      call ga_inquire(g_vecs, type, nbf, nmo)
      if (.not. ga_create(mt_dbl, nbf, nbf, 'aomotmp',32,32,g_tmp1)) 
     &   call errquit('xlm_ao_mo: ga_create failed: g_tmp1',nbf*nbf,
     &       GA_ERR)
      if (.not. ga_create(mt_dbl, nbf, nmo, 'aomotmp2',32,32,g_tmp2)) 
     &   call errquit('xlm_ao_mo: ga_create failed: g_tmp2',nmo*nbf,
     &       GA_ERR)
      if (.not. ma_push_get(mt_dbl,nbf,'xlmtpm',l_tmp, k_tmp))
     $   call errquit('xlm_ao_mo: ma_push_get failed: tmp',nbf,MA_ERR)
c
c     Must transform the LHS index one mpole at a time so 
c     might as well do both at the same time since this will
c     use less memory.
c
      do lm = 1, 6
         call ga_sync()
         do i = 1+ga_nodeid(), nbf, ga_nnodes()
            call ga_get(g_xlm,lm+(i-1)*6,lm+(i-1)*6,1,nbf,
     $           dbl_mb(k_tmp), 1)
            call ga_put(g_tmp1, i, i, 1, nbf, dbl_mb(k_tmp), 1)
         end do
         call ga_sync()
         call ga_dgemm('n','n',nbf,nmo,nbf,1.0d0,g_tmp1,g_vecs,
     $        0.0d0,g_tmp2)
         call ga_dgemm('t','n',nmo,nmo,nbf,1.0d0,g_vecs,g_tmp2,
     $        0.0d0,g_tmp1)
         call ga_sync()
         do i = 1+ga_nodeid(), nmo, ga_nnodes()
            call ga_get(g_tmp1, i, i, 1, nmo, dbl_mb(k_tmp), 1)
            call ga_put(g_xlm,lm+(i-1)*6,lm+(i-1)*6,1,nmo,
     $           dbl_mb(k_tmp), 1)
         end do
         call ga_sync()
      end do
c
c      write(LuOut,*) ' THE MO MPOLES'
c      call ga_print(g_xlm)
c
      if (.not. ga_destroy(g_tmp1)) call errquit('xlm_ao_mo?',1,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('xlm_ao_mo?',2,
     &       GA_ERR)
      if (.not. ma_pop_stack(l_tmp)) call errquit('xlm_ao_mo?',3,
     &       MA_ERR)
c
      end
