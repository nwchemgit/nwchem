      subroutine movecs_anal_so(basis, ilo, ihi, thresh,
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
c
      integer basis
      integer ilo, ihi          ! [input] Range of vectors to print
      double precision thresh   ! [input] Print coeffs with absval >= thresh
      double precision norm, normn 
      integer g_vecs(2)
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
      integer l_vecsre, k_vecsre, i, j, k_tags, l_tags, k_list, l_list
      integer l_vecsim, k_vecsim
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
      logical sym_char_table_so
      external sym_char_table_so 
c
      or2 = .true.
      if (.not. bas_cando_mpoles(basis)) or2 = .false.
c     
      call ga_sync
c      call ga_summarize(.true.)
      if (oscfps) call pstat_on(ps_moanal)
      call ga_inquire(g_vecs, type, nbf, nmo)
      if (or2) then
c
c        local array for x, y, z, x^2, y^2, z^2 and r^2 for each MO
c
         len_r_r2 = 7*nmo
         if (.not. ma_push_get(mt_dbl, len_r_r2, 'nmo:r_and_r2', 
     &      l_r_r2, k_r_r2))call errquit
     &      ('movecs_anal_so: cannot allocate r_r2', len_r_r2, MA_ERR)
         call so_r_and_r2(basis, g_vecs, dbl_mb(k_r_r2))
      endif
c
      if (ga_nodeid() .eq. 0) then
         if (.not. bas_numbf(basis, nbf/2)) call errquit
     $        ('movecs_anal: basis bad?',basis, BASIS_ERR)
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_anal1',l_vecsre,
     &        k_vecsre))
     $        call errquit('movecs_anal: ma 1 failed', nbf, MA_ERR)
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_anal1',l_vecsim,
     &        k_vecsim))
     $        call errquit('movecs_anal: ma 1.1 failed', nbf, MA_ERR)
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
            if (.not. sym_char_table_so(zname, nop, nir, class_dim,
     $           zir, zclass, chars)) call errquit
     $           ('movecs_anal: no char table available ',geom,
     &       GEOM_ERR)
         endif
c     
         call bas_vec_info(basis, byte_mb(k_tags))
         call bas_vec_info(basis, byte_mb(k_tags+16*nbf/2))
c
         write(6,*)
         call util_print_centered(6,title, 40, .true.)
         write(6,*)
c

         do i = ilo, ihi
c
            call ga_get(g_vecs(1), 1, nbf, i, i, dbl_mb(k_vecsre), 1)
            call ga_get(g_vecs(2), 1, nbf, i, i, dbl_mb(k_vecsim), 1)
c
c     Identify significant coefficients and sort by size
c
            n = 0
            cur_thresh = thresh
 111        do j = 0, nbf-1
               norm=sqrt(dbl_mb(k_vecsre+j)**2+dbl_mb(k_vecsim+j)**2)
               if (norm.ge.cur_thresh) then
                  int_mb(k_list + n) = j
                  n = n + 1
               endif
            enddo
            if (n.eq.0 .and. cur_thresh.le.64*thresh) then
               cur_thresh = cur_thresh*2
               goto 111         ! Go back if found nothing to print
            endif
            do j = 0, n-1
               do k = 0, j
                  norm = sqrt(dbl_mb(k_vecsre+int_mb(k_list+k))**2
     &                 + dbl_mb(k_vecsim+int_mb(k_list+k))**2)
                  normn = sqrt(dbl_mb(k_vecsre+int_mb(k_list+j))**2
     &                 + dbl_mb(k_vecsim+int_mb(k_list+j))**2)
                  if (norm.lt.normn) then
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
            write(6,1) i, buf(1:max(inp_strlen(buf),1))
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
               write(6,3) buf(1:max(inp_strlen(buf),1))
 3             format('       ',7x,a) 
            endif
c
c     Output the analysis
c
            write(6,22)
 22         format(1x,2(' Bfn.         Coefficient        Function  ',
     $                  4x,4x))
            write(6,23)
 23         format(1x,2(' ----     -------------------  ------------',
     $                  4x,4x))
            do klo = 0, min(n-1,9), 2
               khi = min(klo+1,n-1)
               write(6,2) (
     $              int_mb(k_list+k)+1,
     $              dbl_mb(k_vecsre+int_mb(k_list+k)),
     $              dbl_mb(k_vecsim+int_mb(k_list+k)),
     $              (byte_mb(k_tags+int_mb(k_list+k)*16+m),m=0,15), 
     $              k = klo,khi)
 2             format(1x,2(i5,2x,f11.6,f11.6,2x,16a1,4x))
            enddo
            write(6,*)
         enddo
         call util_flush(6)
         if (.not. ma_chop_stack(l_vecsim)) call errquit
     $        ('bas_vec_info: ma pop?:l_vecsim', 0, MA_ERR)
         if (.not. ma_chop_stack(l_vecsre)) call errquit
     $        ('bas_vec_info: ma pop?l_vecsre', 0, MA_ERR)
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
      subroutine so_r_and_r2 (basis, g_movecs, r_and_r2)
c
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
c
      integer basis                    ! [input] basis 
      double precision r_and_r2(7,*) ! [output] x,y,z,x^2,y^2,z^2,r^2 for each mo
      integer g_movecs(2)                 ! [input] GA mo vectors
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
     $     32,32,g_xlm)) call errquit('so_r_and_r2 : g_xlm',6*nbf*nbf,
     &       GA_ERR)
c
c     compute dipoles and quadrupole components for each AO
c
      call xlm_so_ao_poles(basis, center, g_xlm)
c
c     compute dipoles and quadrupole components for each MO
c
      call xlm_ao_to_mo_so(g_movecs, g_xlm)
c
c     put global data into local array on node 0 (node 0
c     will be doing the printing of the desired info)
c
      call dfill(nmo*7, 0.0d0, r_and_r2, 1)
      if (ga_nodeid() .eq. 0) then
         do i = 1, nmo
            call ga_get(g_xlm,(i-1)*6+1,(i-1)*6+6,i,i,
     $           r_and_r2(1,i), 1)
         enddo
c
c         write(6,*)' x, y, z, x^2, y^2, z^2 at (0,0,0)'
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
c         write(6,*)' x, y, z, x^2, y^2, z^2, r^2 shifted'
c         call output(r_and_r2, 1, 7, 1, nmo, 7, nmo, 1)
      endif
c      
      if (.not. ga_destroy(g_xlm)) call errquit('so_r_and_r2 : ga?',0,
     &       GA_ERR)
c
      return
      end
      subroutine xlm_so_ao_poles(basis, center, g_xlm)
c
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
c     
      integer basis             ! [input] basis 
      double precision center(3) ! [input] the expansion center
      integer g_xlm              ! [input] GA that will return the mpoles
c
      double precision one, two
      parameter (one=1.d0, two=2.d0)
      double precision autoang, autoang2
      parameter (autoang = 0.529177249)
      parameter (autoang2 = autoang*autoang)
c
      integer mcart, l_xlm, k_xlm
      integer geom
      integer nshell, noperators, maxscr, me, nproc
      integer nbf_max, lmpmax, ishell, jshell, ijshell
      integer ilo, ihi, jlo, jhi, idim, jdim, ind, i, j, l, ioff
      integer l_scr, k_scr, l_mp, k_mp
      integer lmax              ! Maximum value of L = 2
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('multiplole: bad basis', 0, BASIS_ERR)
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
      call ga_zero(g_xlm)
c
      ijshell = -1
      me = ga_nodeid()
      nproc = ga_nnodes()
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
c               write(6,*)' ishell, jshell = ', ishell, jshell
               ind = k_mp
               do l = 0, lmax
                  do i = 1, idim
                     do mcart = 1, ((l+1)*(l+2))/2
                        do j = 1, jdim
                           ioff = k_xlm + 6*(j-1 + jdim*(i-1))
c       write(6,*)' l, i, mcart, j, ind-k_mp, dbl_mb(ind) ',
c     &             l, i, mcart, j, ind-k_mp, dbl_mb(ind)
                           if (l.eq.1.and.mcart.eq.1)then
                              dbl_mb(ioff) = dbl_mb(ind)*autoang
c                              write(6,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.1.and.mcart.eq.2)then
                              ioff = ioff + 1
                              dbl_mb(ioff) = dbl_mb(ind)*autoang
c                              write(6,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.1.and.mcart.eq.3)then
                              ioff = ioff + 2
                              dbl_mb(ioff) = dbl_mb(ind)*autoang
c                              write(6,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.2.and.mcart.eq.1)then
                              ioff = ioff + 3
                              dbl_mb(ioff) = dbl_mb(ind)*autoang2
c                              write(6,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.2.and.mcart.eq.4)then
                              ioff = ioff + 4
                              dbl_mb(ioff) = dbl_mb(ind)*autoang2
c                              write(6,*)' ioff-k_xlm ', ioff-k_xlm
                           endif
                           if (l.eq.2.and.mcart.eq.6)then
                              ioff = ioff + 5
                              dbl_mb(ioff) = dbl_mb(ind)*autoang2
c                              write(6,*)' ioff-k_xlm ', ioff-k_xlm
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
      call ga_sync
c
c      write(6,*) ' THE AO MPOLES '
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
      subroutine xlm_ao_to_mo_so(g_vecs, g_xlm)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer lmax, g_vecs(2), g_xlm
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
      if (.not. ga_create(mt_dbl, nbf, nbf, 'aomotmp',
     $     32,32,g_tmp1)) call errquit('xlm_ao_mo: tmp1',nbf*nbf,
     &       GA_ERR)
      if (.not. ga_create(mt_dbl, nbf, nmo, 'aomotmp2',
     $     32,32,g_tmp2)) call errquit('xlm_ao_mo: tmp2',nmo*nbf,
     &       GA_ERR)
      if (.not. ma_push_get(mt_dbl,nbf,'xlmtpm',l_tmp, k_tmp))
     $     call errquit('xlm_ao_mo: tmp', nbf, MA_ERR)
c
c     Must transform the LHS index one mpole at a time so 
c     might as well do both at the same time since this will
c     use less memory.
c
      do lm = 1, 6
         call ga_sync
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
         call ga_sync
         do i = 1+ga_nodeid(), nmo, ga_nnodes()
            call ga_get(g_tmp1, i, i, 1, nmo, dbl_mb(k_tmp), 1)
            call ga_put(g_xlm,lm+(i-1)*6,lm+(i-1)*6,1,nmo,
     $           dbl_mb(k_tmp), 1)
         end do
         call ga_sync()
      end do
c
c      write(6,*) ' THE MO MPOLES'
c      call ga_print(g_xlm)
c
      if (.not. ga_destroy(g_tmp1)) call errquit('xlm_ao_mo?',1, GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit('xlm_ao_mo?',2, GA_ERR)
      if (.not. ma_pop_stack(l_tmp)) call errquit('xlm_ao_mo?',3,
     &       MA_ERR)
c
      end
