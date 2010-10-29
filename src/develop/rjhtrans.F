      logical function jvltest2(rtdb)
*     
*     $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb
c
      integer basis, geom, nbf
      character*255 movecs      ! Name of movector file
      character*80 title, name_of_basis, scftype
      integer nbf_file, nsets, nmo_file(2)
      logical cct_uhf
      logical movecs_read, movecs_read_header
      external movecs_read, movecs_read_header
c
      integer mult(8,8)
c
      integer nmo, nalpha, nbeta, g_tmp, l_occ, k_occ,
     $     l_aeval, k_aeval, l_beval, k_beval, l_mos_t, k_mos_t,
     $     l_eval, k_eval, l_tmp, k_tmp, nocca, noccb,
     $     nvirta, nvirtb,
     $     l_irs, k_irs
c
      integer nfrozen, nactive
      logical status
c
      logical int_normalize
      external int_normalize
c
      call ifill(8*8, 99999999, mult, 1)
      mult(1,1) = 1             ! All we are using now!
c
      write(6,*)
      write(6,*)
      call util_print_centered(6,'Joop''s cool triples program',
     $     40, .true.)
      write(6,*)
      write(6,*)
      title = ' '
      status = rtdb_cget(rtdb,'title',1,title)
      if (title .ne. ' ') then
         call util_print_centered(6,title,40, .false.)
         write(6,*)
         write(6,*)
      endif
c     
c     load the geometry/basis set and get info
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0, RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $        call errquit('scf_init: no ao basis set', 0, RTDB_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0, BASIS_ERR)
      if (.not.int_normalize(rtdb,basis))
     $           call errquit('scf:int_normalize failed', 0, INT_ERR)
c
c     Read the MO vectors and evals from a UHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('jvltest2: failed to read movecs header',0, DISK_ERR)
c
      if (nsets.ne.2) call errquit('Joop says: UHF ONLY!',0, INPUT_ERR)
c
      write(6,*) ' Read movecs header from ', 
     $     movecs(1:inp_strlen(movecs))
      write(6,*) ' Job title :                ', 
     $     title(1:inp_strlen(title))
      write(6,*) ' Basis name:                ', 
     $     name_of_basis(1:inp_strlen(name_of_basis))
      nmo = nmo_file(1)
      if (nmo_file(1).ne.nmo_file(2)) call errquit('Uh?',0, INPUT_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nalpha', mt_int, 1, nalpha))
     $     call errquit('nalpha?',0, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nbeta', mt_int, 1, nbeta))
     $     call errquit('nbeta?',0, RTDB_ERR)
c
c     Frozen core info
c
      if (rtdb_get(rtdb,'ccsd:frozen core:freeze by atoms',mt_log, 1, 
     $     status)) then
         if (.not. geom_num_core(geom, nfrozen)) 
     $        call errquit('Joop says: geom_num_core?',0, GEOM_ERR)
      else if (rtdb_get(rtdb, 'ccsd:frozen core', MT_INT, 1,
     $        nfrozen)) then
      else
         nfrozen = 0
      endif
c
      nactive = nmo - nfrozen
      nocca   = nalpha - nfrozen
      noccb   = nbeta  - nfrozen
      nvirta  = nactive - nocca
      nvirtb  = nactive - noccb
c
      write(6,*) ' No. of alpha electrons:    ', nalpha
      write(6,*) ' No. of beta  electrons:    ', nbeta
      write(6,*) ' No. of frozen core:        ', nfrozen
      write(6,*) ' No. of molecular orbitals: ', nmo
      write(6,*) ' No. of active MOs:         ', nactive
      write(6,*) ' No. of basis functions:    ', nbf
c
c     Generate local copies of the TRANSPOSED, ACTIVE molecular orbitals.
c
c     The active MOs and eigenvalues will be combined to form single
c     arrays with the first nactive entries being alpha and the next beta.
c     
      if (.not. ma_push_get(mt_dbl, nbf*nactive*2,'mos_t', 
     $     l_mos_t, k_mos_t)) call errquit('Joop: mos', 2*nbf*nactive,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, nactive*2,'evals', 
     $     l_eval, k_eval)) call errquit('Joop: evals', 2*nactive,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nactive*2, 'irs', 
     $     l_irs, k_irs)) call errquit('Joop: irs', 2*nactive,
     &       MA_ERR)
      call ifill(2*nactive, 1, int_mb(k_irs), 1)
c
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'aeval',l_aeval, k_aeval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'beval',l_beval, k_beval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nactive,'tmp', 
     $     l_tmp, k_tmp)) call errquit('Joop: tmp', nbf*nactive, MA_ERR)
c
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 0, 0, g_tmp))
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_aeval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, nfrozen+1, nmo, dbl_mb(k_tmp), nbf)
      call util_transpose(dbl_mb(k_tmp),nbf,dbl_mb(k_mos_t),nactive,
     $     nbf,nactive)
      call dcopy(nactive, dbl_mb(k_aeval+nfrozen), 1, dbl_mb(k_eval), 1)
c
      if (.not. movecs_read(movecs, 2, dbl_mb(k_occ), dbl_mb(k_beval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_tmp), nbf)
      call util_transpose(dbl_mb(k_tmp), nbf,
     $     dbl_mb(k_mos_t+nactive*nbf), nactive, nbf, nactive)
      call dcopy(nactive, dbl_mb(k_beval+nfrozen), 1, 
     $     dbl_mb(k_eval+nactive), 1)
c
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0, GA_ERR)
      if (.not. ma_chop_stack(l_occ)) call errquit(' ma chop?', 0,
     &       MA_ERR)
c
      write(6,*) ' Alpha active eigenvalues '
      call output(dbl_mb(k_eval),1,nactive,1,1,nactive,1,1)
      write(6,*) ' Beta active eigenvalues '
      call output(dbl_mb(k_eval+nactive),1,nactive,1,1,nactive,1,1)
      write(6,*) ' Alpha MOs transposed'
      call output(dbl_mb(k_mos_t),1,nactive,1,nbf,nactive,nbf,1)
      write(6,*) ' Beta  active MOs transposed'
      call output(dbl_mb(k_mos_t+nbf*nactive),1,nactive,1,nbf,
     $     nactive,nbf,1)
c
c     Now make the darn integral files
c
      call rjh_init_file_info()
      call rjh_make_integral_files(
     $     rtdb, basis, 
     $     nactive, nbf, nocca, noccb, 
     $     dbl_mb(k_mos_t), dbl_mb(k_eval))
c
c     Wake up Joop
c
      call cct_uhf_address(1, mult, int_mb(k_irs), 
     1     nocca, noccb, nactive)
c
c     Trying reading something from a file ... all 3-c integrals
c
      if (.not. ga_create(mt_dbl, nvirta*nocca, nvirtb*nvirtb, 
     $     'tmp', 0, 0, g_tmp))
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
c
      call get_cct(g_tmp, '<ba||ek>(bkea)', 
     $     1, nvirta*nocca, 1, 1, 
     $     nactive+noccb+1, nactive+nactive, 2, 
     $     nactive+noccb+1, nactive+nactive, 2)
c
      write(6,*) ' spins b k a e = ', 1, 1, 2, 2
      call ga_print(g_tmp)
c
      call get_cct(g_tmp, '<ba||ek>(bkea)', 
     $     1, nvirta*nocca, 2, 2, 
     $     noccb+1, nactive, 1, 
     $     noccb+1, nactive, 1)
c
      write(6,*) ' spins b k a e = ', 2, 2, 1, 1
      call ga_print(g_tmp)
c
      if (.not. ga_destroy(g_tmp)) call errquit('ga?',0, GA_ERR)
c
      status = cct_uhf(rtdb)
c
c     Tidy up
c
      if (.not. bas_destroy(basis)) call errquit(' bas ?',0, BASIS_ERR)
      if (.not. geom_destroy(geom)) call errquit(' geom ?',0, GEOM_ERR)
      if (.not. ma_chop_stack(l_mos_t)) call errquit(' ma ? ',0, MA_ERR)
c
      call rjh_tidy_file_info
c
      jvltest2 = status
c
      end
c$$$      subroutine rjh_test_trans(rtdb, basis, amos, bmos,
c$$$     $     nmo, nalpha, nbeta, nbf)
c$$$      implicit none
c$$$#include "errquit.fh"
c$$$#include "mafdecls.fh"
c$$$c
c$$$      integer rtdb, basis
c$$$      integer nmo, nalpha, nbeta, nbf
c$$$      double precision amos(nbf,nmo), bmos(nbf,nmo), sum
c$$$
c$$$c
c$$$      integer lenfull, k_full, l_full, k_full2, l_full2
c$$$      integer l_aoint, k_aoint, i
c$$$c
c$$$      if (nmo .ne. nbf) call errquit(' test hack ',0)
c$$$c
c$$$      lenfull = nbf**4
c$$$      if (.not. ma_push_get(mt_dbl,lenfull,'full',l_full, k_full))
c$$$     $     call errquit('ma full', lenfull)
c$$$      if (.not. ma_push_get(mt_dbl,lenfull,'full2',l_full2, k_full2))
c$$$     $     call errquit('ma full2', lenfull)
c$$$      if (.not. ma_push_get(mt_dbl,lenfull,'aoint',l_aoint, k_aoint))
c$$$     $     call errquit('ma aoint', lenfull)
c$$$c
c$$$      call util_inplace_transpose(nbf, amos) ! TRANSPOSE
c$$$      call util_inplace_transpose(nbf, bmos)
c$$$c
c$$$      call rjh_full_transform(rtdb, basis, nbf, nbf, nbf, nbf,
c$$$     $     nmo, nmo, nmo, nmo,
c$$$     $     amos, amos, amos, amos, dbl_mb(k_full), 'chargecloud')
c$$$c
c$$$      call rjh_make_integral_files(rtdb, basis, 
c$$$     $     nmo, nalpha, nbeta, amos, bmos)
c$$$      call rjh_tidy_file_info()
c$$$c
c$$$      call util_inplace_transpose(nbf, amos) !UNDO TRANSPOSE
c$$$      call util_inplace_transpose(nbf, bmos)
c$$$c
c$$$c$$$      call rjh_debug_print('full',
c$$$c$$$     $     dbl_mb(k_full), nbf, nbf, nbf, nbf)
c$$$c
c$$$      call rjh_all_ao_integrals(rtdb,basis,nbf,'mulliken',
c$$$     $     dbl_mb(k_aoint))
c$$$c
c$$$      call rjh_uhf_trans_lo_mem(nbf, nbf, int_mb(k_full2),
c$$$     $     amos, dbl_mb(k_aoint), dbl_mb(k_full2))
c$$$c$$$      call rjh_debug_print('full2',
c$$$c$$$     $     dbl_mb(k_full2), nbf, nbf, nbf, nbf)
c$$$c
c$$$      sum = 0.0d0
c$$$      do i = 0, lenfull-1
c$$$         sum = sum + abs(dbl_mb(k_full2+i)-dbl_mb(k_full2+i))
c$$$      enddo
c$$$      write(6,*) ' Difference between old and new ', sum
c$$$c
c$$$      if (.not. ma_chop_stack(l_full)) call errquit('ma99',0)
c$$$c
c$$$      end
      subroutine rjh_debug_print(string,full, n1, n2, n3, n4)
      implicit none
      character*(*) string
      integer n1, n2, n3, n4
      double precision full(n1, n2, n3, n4)
c
      integer p, q, r, s
      write(6,*) ' DEBUG FOR ', string, n1, n2, n3, n4
      do s = 1, n4
         do r = 1, n3
            do q = 1, n2
               do p = 1, n1
                  if (abs(full(p,q,r,s)).gt.1e-6) then
                     write(6,1) p,q,r,s,full(p,q,r,s)
 1                   format(4i5,2x,f12.6)
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine rjh_init_file_info()
      implicit none
#include "rjhfileinfo.fh"
c
      integer i
c
      do i = 1, maxfile
         fds(i) = -1
      enddo
c
      end
      subroutine rjh_tidy_file_info()
      implicit none
#include "errquit.fh"
#include "rjhfileinfo.fh"
#include "eaf.fh"
c
      integer i
c
      do i = 1, maxfile
         if (fds(i).ne.-1) then
            if (eaf_close(fds(i)) .ne. 0) call errquit
     $           ('rjh_tidy_file_info: eaf_close failed',i, UNKNOWN_ERR)
         endif
      enddo
c
      end
      subroutine rjh_make_integral_files(rtdb, basis, 
     $     nmo, nbf, nalpha, nbeta, mos_t, evals)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rjhfileinfo.fh"
      integer k_all, norb
      common /c_all/norb, k_all
c
      integer rtdb, basis
      integer nmo, nalpha, nbeta, nbf
      double precision mos_t(nmo,nbf,2), evals(nmo,2)
c
c     <ba||ek> = <ba|ek> - <ba|ke>
c     Write the next 6 as file(b,k,e,a)
c
c     1.  <ba||ek>           b=a=e=k=alpha         baek = 1111
c     2.  <ba||ek>           b=a=e=k=beta          baek = 2222
c     3.  <ba|ek>            b=e=alpha, a=k=beta   baek = 1212
c     4.  <ba|ek>            b=e=beta,  a=k=alpha  baek = 2121
c     5. -<ba|ke> = (bk|ae)  b=k=alpha, a=e=beta   baek = 1221 sign=-1
c     6. -<ba|ke> = (bk|ae)  b=k=beta,  a=e=alpha  baek = 2112 sign=-1
c
c     <ij||ab>
c     Write the next set as (i,j,a,b) restricting i>=j, a>=b if possible.
c
c     7.  <ij||ab>           i=j=a=b=alpha         ijab = 1111
c     8.  <ij||ab>           i=j=a=b=beta          ijab = 2222
c     9.  <ij|ab>            i=a=alpha j=b=beta    ijab = 1212
c
c     T(ij,ab) written same as (i,j,a,b)
c     Write the next set as (i,j,a,b) restricting i>=j, a>=b if possible.
c
c     10. <ij||ab>           i=j=a=b=alpha         ijab = 1111
c     11. <ij||ab>           i=j=a=b=beta          ijab = 2222
c     12. <ij|ab>            i=a=alpha j=b=beta    ijab = 1212
c
c     No use of spatial symmetry at present.  Assume can hold
c     everything in memory.
c
      integer lenfull, l_full, k_full
      integer occa, occb        ! No. of occupied of each spin
      integer virta, virtb      ! No. of virtuals of each spin
c
      call rjh_init_file_info()
      occa  = nalpha
      occb  = nbeta
      virta = nmo-occa
      virtb = nmo-occb
c
*      lenfull = max(occa,occb)*max(virta,virtb)**3
      lenfull = nmo**4
c
c     Slowly move to have the following using Joop-like spin cases.
c
      if (.not. ma_push_get(mt_dbl,lenfull,'full',l_full, k_full))
     $     call errquit('ma full', lenfull, MA_ERR)

*      call rjh_full_transform(rtdb, basis,
*     $     nmo, nmo, nmo, nmo,
*     $     nmo, nmo, nmo, nmo,
*     $     mos_t(1,1,1), mos_t(1,1,1), 
*     $     mos_t(1,1,1), mos_t(1,1,1),
*     $     dbl_mb(k_full), 'Chargecloud')
*      k_all = k_full
*      norb = nmo

      if (.not. ma_push_get(mt_dbl,lenfull,'full',l_full, k_full))
     $     call errquit('ma full', lenfull, MA_ERR)

c
c     1. <ba||ek>           b=a=e=k=alpha         baek = 1111
c
      call rjh_full_transform(rtdb, basis,
     $     virta, virta, virta, occa,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(occa+1,1,1), mos_t(occa+1,1,1), 
     $     mos_t(occa+1,1,1), mos_t(1,1,1),
     $     dbl_mb(k_full), 'LeftAsymDirac')
      call rjh_integ_write(dbl_mb(k_full), virta, virta, virta, occa,
     $     '1432', '<ba||ek>', 1111, 1.0d0)
c
c     2. <ba||ek>           b=a=e=k=beta          baek = 2222
c
      call rjh_full_transform(rtdb, basis,
     $     virtb, virtb, virtb, occb,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(occb+1,1,2), mos_t(occb+1,1,2), 
     $     mos_t(occb+1,1,2), mos_t(1,1,2),
     $     dbl_mb(k_full), 'LeftAsymDirac')
      call rjh_integ_write(dbl_mb(k_full), virtb, virtb, virtb, occb,
     $     '1432', '<ba||ek>', 2222, 1.0d0)
c
c     3. <ba|ek>            b=e=alpha, a=k=beta   baek = 1212
c
      call rjh_full_transform(rtdb, basis,
     $     virta, virtb, virta, occb,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(occa+1,1,1), mos_t(occb+1,1,2), 
     $     mos_t(occa+1,1,1), mos_t(1,1,2),
     $     dbl_mb(k_full), 'Dirac')
      call rjh_integ_write(dbl_mb(k_full), virta, virtb, virta, occb,
     $     '1432', '<ba||ek>', 1212, 1.0d0)
c
c     4. <ba|ek>            b=e=beta,  a=k=alpha  baek = 2121
c
      call rjh_full_transform(rtdb, basis,
     $     virtb, virta, virtb, occa,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(occb+1,1,2), mos_t(occa+1,1,1), 
     $     mos_t(occb+1,1,2), mos_t(1,1,1),
     $     dbl_mb(k_full), 'Dirac')
      call rjh_integ_write(dbl_mb(k_full), virtb, virta, virtb, occa,
     $     '1432', '<ba||ek>', 2121, 1.0d0)
c
c     5. <ba|ke> = (bk|ae)  b=k=alpha, a=e=beta   baek = 1221 sign=-1
c
      call rjh_full_transform(rtdb, basis,
     $     virta, occa, virtb, virtb,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(occa+1,1,1), mos_t(1,1,1), 
     $     mos_t(occb+1,1,2), mos_t(occb+1,1,2),
     $     dbl_mb(k_full), 'Chargecloud')
      call rjh_integ_write(dbl_mb(k_full), virta, occa, virtb, virtb,
     $     '1234', '<ba||ek>', 1221, -1.0d0)
c
c     6. <ba|ke> = (bk|ae)  b=k=beta,  a=e=alpha  baek = 2112 sign=-1
c
      call rjh_full_transform(rtdb, basis,
     $     virtb, occb, virta, virta,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(occb+1,1,2), mos_t(1,1,2), 
     $     mos_t(occa+1,1,1), mos_t(occa+1,1,1),
     $     dbl_mb(k_full), 'Chargecloud')
      call rjh_integ_write(dbl_mb(k_full), virtb, occb, virta, virta,
     $     '1234', '<ba||ek>', 2112, -1.0d0)
c
c     7.  <ij||ab>           i=j=a=b=alpha         ijab = 1111
c     i>=j, a>=b
c
      call rjh_full_transform(rtdb, basis,
     $     occa, occa, virta, virta,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(1,1,1), mos_t(1,1,1), 
     $     mos_t(occa+1,1,1), mos_t(occa+1,1,1),
     $     dbl_mb(k_full), 'LeftAsymDirac')
      call rjh_integ_write(dbl_mb(k_full), occa, occa, virta, virta,
     $     '1<=234', '<ij||ab>', 1111, 1.0d0)
c
c     8.  <ij||ab>           i=j=a=b=beta          ijab = 2222
c
      call rjh_full_transform(rtdb, basis,
     $     occb, occb, virtb, virtb,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(1,1,2), mos_t(1,1,2), 
     $     mos_t(occb+1,1,2), mos_t(occb+1,1,2),
     $     dbl_mb(k_full), 'LeftAsymDirac')
      call rjh_integ_write(dbl_mb(k_full), occb, occb, virtb, virtb,
     $     '1<=234', '<ij||ab>', 2222, 1.0d0)
c
c     9.  <ij|ab>            i=a=alpha j=b=beta    ijab = 1212
c     
      call rjh_full_transform(rtdb, basis,
     $     occa, occb, virta, virtb,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(1,1,1), mos_t(1,1,2), 
     $     mos_t(occa+1,1,1), mos_t(occb+1,1,2),
     $     dbl_mb(k_full), 'Dirac')
      call rjh_integ_write(dbl_mb(k_full), occa, occb, virta, virtb,
     $     '1234', '<ij||ab>', 1212, 1.0d0)
c
c     10. T(ij,ab)           i=j=a=b=alpha         ijab = 1111
c     i>=j, a>=b
c
      call rjh_full_transform(rtdb, basis,
     $     occa, occa, virta, virta,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(1,1,1), mos_t(1,1,1), 
     $     mos_t(occa+1,1,1), mos_t(occa+1,1,1),
     $     dbl_mb(k_full), 'LeftAsymDirac')
      call rjh_denoms(dbl_mb(k_full), 
     $     occa, occa, virta, virta, 
     $     evals(1,1), evals(1,1), evals(occa+1,1), evals(occa+1,1))
      call rjh_integ_write(dbl_mb(k_full), occa, occa, virta, virta,
     $     '1<=234', 'T(ij,ec)', 1111, 1.0d0)
c
c     11. T(ij,ab)           i=j=a=b=beta          ijab = 2222
c
      call rjh_full_transform(rtdb, basis,
     $     occb, occb, virtb, virtb,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(1,1,2), mos_t(1,1,2), 
     $     mos_t(occb+1,1,2), mos_t(occb+1,1,2),
     $     dbl_mb(k_full), 'LeftAsymDirac')
      call rjh_denoms(dbl_mb(k_full), 
     $     occb, occb, virtb, virtb, 
     $     evals(1,2), evals(1,2), evals(occb+1,2), evals(occb+1,2))
      call rjh_integ_write(dbl_mb(k_full), occb, occb, virtb, virtb,
     $     '1<=234', 'T(ij,ec)', 2222, 1.0d0)
c
c     12. T(ij,ab)            i=a=alpha j=b=beta    ijab = 1212
c     
      call rjh_full_transform(rtdb, basis,
     $     occa, occb, virta, virtb,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(1,1,1), mos_t(1,1,2), 
     $     mos_t(occa+1,1,1), mos_t(occb+1,1,2),
     $     dbl_mb(k_full), 'Dirac')
      call rjh_denoms(dbl_mb(k_full), 
     $     occa, occb, virta, virtb, 
     $     evals(1,1), evals(1,2), evals(occa+1,1), evals(occb+1,2))
      call rjh_integ_write(dbl_mb(k_full), occa, occb, virta, virtb,
     $     '1234', 'T(ij,ec)', 1212, 1.0d0)
c
c     13. T(ij,ab)            i=b=alpha j=a=beta    ijab = 1221
c     
      call rjh_full_transform(rtdb, basis,
     $     occa, occb, virtb, virta,
     $     nmo, nmo, nmo, nmo,
     $     mos_t(1,1,1), mos_t(1,1,2), 
     $     mos_t(occb+1,1,2), mos_t(occa+1,1,1),
     $     dbl_mb(k_full), 'Dirac')
      call rjh_denoms(dbl_mb(k_full), 
     $     occa, occb, virta, virtb, 
     $     evals(1,1), evals(1,2), evals(occb+1,2), evals(occa+1,1))
      call rjh_integ_write(dbl_mb(k_full), occa, occb, virtb, virta,
     $     '1234', 'T(ij,ec)', 1221, 1.0d0)
c
      if (.not. ma_chop_stack(l_full)) call errquit
     $     ('rjh_write_integ: ma chop?',0, MA_ERR)
c
      end
      subroutine rjh_denoms(full, dim1, dim2, dim3, dim4,
     $     e1, e2, e3, e4)
      implicit none
      integer dim1, dim2, dim3, dim4
      double precision full(dim1, dim2, dim3, dim4)
      double precision e1(*), e2(*), e3(*), e4(*)
c
      integer i, j, a, b
      double precision energy
c
      energy = 0.0d0
      do b = 1, dim4
         do a = 1, dim3
            do j = 1, dim2
               do i = 1, dim1
                  full(i,j,a,b) = full(i,j,a,b) / 
     $                 (e1(i) + e2(j) - e3(a) - e4(b))
                  energy = energy + full(i,j,a,b)**2*
     $                 (e1(i) + e2(j) - e3(a) - e4(b))
               enddo
            enddo
         enddo
      enddo
c
      write(6,*) ' ENERGY from T2(1) before writing ', energy
c
      end
      subroutine rjh_integ_write(full, dim1, dim2, dim3, dim4,
     $     order, partstub, spin, sign)
      implicit none
#include "errquit.fh"
#include "eaf.fh"
#include "inp.fh"
#include "rjhfileinfo.fh"
      integer dim1, dim2, dim3, dim4, spin
      double precision full(dim1, dim2, dim3, dim4), sign
      character*(*) order, partstub
c
      integer fd, ierr, i2, i3, i4, filenum
      double precision offset
      character*255 filename, errmsg
      character*16 stub
c
      stub = ' '
      write(stub,'(a,''-'',i4)') partstub, spin
      call util_file_name(stub, .true., .false., filename)
c
      call dscal(dim1*dim2*dim3*dim4, sign, full, 1)
      if (eaf_open(filename, eaf_rw, fd) .ne. 0)
     $     call errquit('rjh_integ_write: opening file?',0, DISK_ERR)
      offset = 0.0d0
c
**      call rjh_debug_print(stub,full,dim1,dim2,dim3,dim4)
c    
      if (order .eq. '1234') then ! full(i1,i2,i3,i4) -> file(i1,i2,i3,i4)
         do i4 = 1, dim4
            do i3 = 1, dim3
               do i2 = 1, dim2
                  ierr = eaf_write(fd, offset, full(1,i2,i3,i4), 8*dim1)
                  if (ierr .ne. 0) goto 100
                  offset = offset + dim1*8
               enddo
            enddo
         enddo
      else if (order .eq. '1<=234') then
         if (dim1.ne.dim2) call errquit('rjh_integ_write: order?',0,
     &       INT_ERR)
         do i4 = 1, dim4
            do i3 = 1, dim3
               do i2 = 1, dim2
                  ierr = eaf_write(fd, offset, full(1,i2,i3,i4), 8*i2)
                  if (ierr .ne. 0) goto 100
                  offset = offset + 8*i2
               enddo
            enddo
         enddo
      else if (order .eq. '1432') then ! full(i1,i2,i3,i4) -> file(i1,i4,i3,i2)
         do i2 = 1, dim2
            do i3 = 1, dim3
               do i4 = 1, dim4
                  ierr = eaf_write(fd, offset, full(1,i2,i3,i4), 8*dim1)
                  if (ierr .ne. 0) goto 100
                  offset = offset + 8*dim1
               enddo
            enddo
         enddo
      else
         call errquit('rjh_integ_write: unknown order',0, INT_ERR)
      endif
c     
c     Save all of the info in common
c
      do filenum = 1, maxfile
         if (fds(filenum).eq.-1) goto 10
      enddo
      call errquit('rjh_integ_write: too many open files?', 0, INT_ERR)
 10   fds(filenum)       = fd
      stubs(filenum)     = stub
      filenames(filenum) = filename
      write(6,11) filenum, stub(1:inp_strlen(stub)), 
     $     filename(1:inp_strlen(filename)), order, dim1, dim2, 
     $     dim3, dim4, offset
 11   format(' wrote file=',i2,2x,a,2x,a,' ',a,2x,4i5,f10.1)
c
      return
c
 100  call eaf_errmsg(ierr, errmsg)
      write(6,*) ' filename  ', filename(1:inp_strlen(filename))
      write(6,*) ' IO offset ', offset
      write(6,*) ' IO error message ',errmsg(1:inp_strlen(errmsg))
      call errquit('rij_integ_write: write failed',0, DISK_ERR)
c
      end
      subroutine get_cct(g_a, type, 
     $     i1i2_lo, i1i2_hi, i1_spin, i2_spin, 
     $     i3_lo, i3_hi, i3_spin, 
     $     i4_lo, i4_hi, i4_spin)
      implicit none
#include "errquit.fh"
#include "eaf.fh"
#include "rjhfileinfo.fh"
#include "cct_table_uhf.fh"
#include "mafdecls.fh"
#include "inp.fh"

      integer k_all, norb
      common /c_all/norb, k_all

c
      integer g_a, i1i2_lo, i1i2_hi, i1_spin, i2_spin,
     $     i3_lo, i3_hi, i3_spin, i4_lo, i4_hi, i4_spin
      character*(*) type
c
c     Read into the global array which is dimensioned 
c     .   ga(i1i2_lo:i1i2_hi, (i3_lo:i3_hi)*(i4_lo:i4_hi))
c     from the list described by type which is logically stored
c     as an array (i1, i2, i3, i4) as follows for given spatial 
c     symmetries of i3 and i4.
c
c     do i4 = i4_lo, i4_hi         ! All the same symmetry
c     .  do i3 = i3_lo, i3_hi      ! All the same symmetry
c     .     do = sym of i2
c     .        -> sym of i1
c     .        do i2 in symmetry block
c     .           do i1 in symmetry block
c     .              ...
c
c     Spatial symmetry is not yet used.
c
c     Indices 1 and 2 will be summed over and must be weighted.
c     Index 3 will always be a non-weighted index (e/f/m/n)
c
c     The input ranges are over SPIN orbitals
c
      integer filenum, fd
      character*8 order
      character*80 errmsg
      integer dim1, dim2, dim3, dim12, i4_lo_spatial, i4_hi_spatial,
     $     i3_lo_spatial, i3_hi_spatial, i1_lo, i2_lo, i1_hi, i2_hi,
     $     reclen, i1, i2, i3, i4, i3i4, ierr
      integer i1_occ, i2_occ, i3_occ, i4_occ, ind, l_buff, k_buff
      integer i1_top, i3_base, i4_base
      double precision offset, xxx, test
c     
c     0) determine which physical file
c     1) determine order of indices in file on disk
c     2) determine dimensions of each index
c     3) read the data
c     4) apply denominators if reading first order T2
c     5) apply weights to summation indices
c
      if (type .eq. '<ba||ek>(bkea)') then
c
c     i1=b,   i2=k,  i3=e,  i4=a
c
         call rjh_which_file('<ba||ek>',
     $        i1_spin, i4_spin, i3_spin, i2_spin, filenum)
         order = '1234'
         i1_occ = 2             ! b
         i2_occ = 1             ! k
         i3_occ = 2             ! e
         i4_occ = 2             ! a
      else if (type .eq. 'T(ij,ce)(ijec)') then   ! ? ORDERING OF INDICES
c
c     i1=i,   i2=j,  i3=e,  i4=c
c
         call rjh_which_file('T(ij,ec)',
     $        i1_spin, i2_spin, i3_spin, i4_spin, filenum)
         if (i1_spin .eq. i2_spin) then
            order = '1<=234'
         else
            order = '1234'
         endif
         i1_occ = 1             ! i
         i2_occ = 1             ! j
         i3_occ = 2             ! c/e
         i4_occ = 2             ! c/e
      else
         call errquit('get_cct: unknown type',0, INPUT_ERR)
      endif
c
      fd = fds(filenum)
c
      dim12 = lenxy(1,i1_spin,i2_spin,i1_occ+i2_occ)
      dim1  = n_occ_virt(1,i1_spin,i1_occ)
      dim2  = n_occ_virt(1,i2_spin,i2_occ)
      dim3  = n_occ_virt(1,i3_spin,i3_occ)
c
      i4_lo_spatial = i4_lo
      if (i4_spin .eq. 2) i4_lo_spatial = i4_lo - orb
      i4_hi_spatial = i4_lo_spatial + (i4_hi - i4_lo)
      i3_lo_spatial = i3_lo
      if (i3_spin .eq. 2) i3_lo_spatial = i3_lo - orb
      i3_hi_spatial = i3_lo_spatial + (i3_hi - i3_lo)
c
      i1_lo = s_occ_virt(1,i1_spin,i1_occ)
      i1_hi = s_occ_virt(2,i1_spin,i1_occ)-1
      i2_lo = s_occ_virt(1,i2_spin,i2_occ)
      i2_hi = s_occ_virt(2,i2_spin,i2_occ)-1
      i3_base = s_occ_virt(1,i3_spin,i3_occ)
      if (i3_spin .eq. 2) i3_base = i3_base - orb
      i4_base = s_occ_virt(1,i4_spin,i4_occ)
      if (i4_spin .eq. 2) i4_base = i4_base - orb

c
      reclen = i1i2_hi - i1i2_lo + 1
      if (.not. ma_push_get(mt_dbl,reclen,'buffer',l_buff, k_buff))
     $     call errquit('ma buff', reclen, MA_ERR)
c
*      write(6,*) ' GET_CCT: ', type, i1i2_lo, i1i2_hi, i1_spin, i2_spin, 
*     $     i3_lo, i3_hi, i3_spin, 
*     $     i4_lo, i4_hi, i4_spin, ' ',order
c
*      write(6,*) ' FAC'
*      call output(fac,1,orb,1,2,orb,2,1)
c
      if (order .eq. '1234' .or. order .eq. '1<=234') then
         do i4 = i4_lo_spatial, i4_hi_spatial
            do i3 = i3_lo_spatial, i3_hi_spatial
               offset = (i1i2_lo - 1 + dim12*(i3-i3_base + 
     $              dim3*(i4-i4_base)))*8.0d0
               ierr = eaf_read(fd, offset, dbl_mb(k_buff), 8*reclen)
               if (ierr .ne. 0) goto 100
*               write(6,*) ' THIS IS WHAT I READ ', offset, reclen
*               call output(dbl_mb(k_buff), 1, reclen, 1, 1, reclen,1,1)
c
               if (i1i2_lo .ne. 1) stop 99999
               ind = 0
               do i2 = i2_lo, i2_hi
                  i1_top = i1_hi
                  if (order .eq. '1<=234') i1_top = i2
                  do i1 = i1_lo, i1_top
                     dbl_mb(k_buff+ind) = dbl_mb(k_buff+ind) *
     $                    fac(i1) * fac(i2)
                     ind = ind + 1
                  enddo
               enddo
c
               if (ind .ne. reclen) call errquit('ind?',ind,
     &       UNKNOWN_ERR)
c
               i3i4 = 1 + i3-i3_lo_spatial + 
     $              (i3_hi-i3_lo+1)*(i4-i4_lo_spatial)
c
               call ga_put(g_a, 1, reclen, i3i4, i3i4, 
     $              dbl_mb(k_buff), 1)
c
*               write(6,*) ' ga_put ', 1, reclen, i3i4, i3i4
*               call output(dbl_mb(k_buff), 1, reclen, 1, 1, reclen,1,1)
c
            enddo
         enddo
      endif
c
      if (i1_spin .eq. 2) then
         i1_lo = i1_lo - norb
         i1_hi = i1_hi - norb
      endif
      if (i2_spin .eq. 2) then
         i2_lo = i2_lo - norb
         i2_hi = i2_hi - norb
      endif
c$$$      if (type .eq. '<ba||ek>(bkea)' ) then
c$$$         write(6,*) ' SPINS ', i1_spin, i2_spin, i3_spin, i4_spin
c$$$         write(6,*) ' LOS ', i1_lo, i2_lo, i3_lo_spatial, 
c$$$     $        i4_lo_spatial
c$$$         write(6,*) ' HIS ', i1_hi, i2_hi, i3_hi_spatial, i4_hi_spatial
c$$$         do i4 = i4_lo_spatial, i4_hi_spatial
c$$$            do i3 = i3_lo_spatial, i3_hi_spatial
c$$$               i3i4 = 1 + i3-i3_lo_spatial + 
c$$$     $              (i3_hi-i3_lo+1)*(i4-i4_lo_spatial)
c$$$               ind = 1
c$$$               do i2 = i2_lo, i2_hi
c$$$                  i1_top = i1_hi
c$$$                  if (order .eq. '1<=234') i1_top = i2
c$$$                  do i1 = i1_lo, i1_top
c$$$                     call ga_get(g_a,ind,ind,i3i4,i3i4,xxx,1)
c$$$                     ind = ind + 1
c$$$c     
c$$$c     XXX should be <i1 i4||i3 i2> = (i1 i3|i4 i2) - (i1 i2|i4 i3)
c$$$c     .              b  a   e  k   =  b  e  a  k      b  k  a  e
c$$$c
c$$$                     test = 0.0d0
c$$$                     if (i1_spin .eq. i3_spin) test = test + 
c$$$     $                    dbl_mb(k_all + 
c$$$     $                    i1-1+norb*(i3-1+norb*(i2-1+norb*(i4-1))))
c$$$                     if (i1_spin .eq. i2_spin) test = test -
c$$$     $                    dbl_mb(k_all + 
c$$$     $                    i1-1+norb*(i2-1+norb*(i4-1+norb*(i3-1))))
c$$$                     if (abs(test-xxx).gt.1d-8) then
c$$$                        write(6,1) i1,i2,i3,i4, test, xxx
c$$$ 1                      format(4i5,2f12.8)
c$$$                     else if (abs(test).gt.1d-8) then
c$$$                        write(6,2) i1,i2,i3,i4
c$$$ 2                      format(4i5,'  ok')
c$$$                     endif
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      endif
c
      if (.not. ma_pop_stack(l_buff)) call errquit
     $     ('get_cct: ma corrupted', 0, MA_ERR)
c
      return
c
 100  call eaf_errmsg(ierr, errmsg)
      write(6,*) ' filename  ', 
     $     filenames(filenum)(1:inp_strlen(filenames(filenum)))
      write(6,*) ' IO offset ', offset
      write(6,*) ' IO error message ',errmsg(1:inp_strlen(errmsg))
      call errquit('get_cct: read failed',0, DISK_ERR)
c      
      end
      subroutine rjh_which_file(partstub, 
     $     i1_spin, i2_spin, i3_spin, i4_spin, 
     $     filenum)
      implicit none
#include "errquit.fh"
#include "rjhfileinfo.fh"
#include "inp.fh"
c
      integer i1_spin, i2_spin, i3_spin, i4_spin, filenum
      character*(*) partstub
c
c     Determine which physical file needs to be read by examining
c     the stubs and spin labels
c
      character*16 stub
      integer ls, spin, i
c
      spin = i1_spin*1000 + i2_spin*100 + i3_spin*10 + i4_spin
      write(stub,'(a,''-'',i4)') partstub, spin
      ls = inp_strlen(stub)
c
      do i = 1, maxfile
         if (fds(i) .ne. -1) then
            if (stub(1:ls).eq.stubs(i)(1:ls)) then
               filenum = i
               write(6,*) stub, ' -> ', i, filenames(i)
               return
            endif
         endif
      enddo
c
      write(6,*) ' STUB ', stub
      call errquit(' rjh_which_file: failed to find ', 0, DISK_ERR)
c
      end
      subroutine rjh_full_transform(
     $     rtdb, basis, 
     $     n1, n2, n3, n4,
     $     ld1, ld2, ld3, ld4,
     $     c1t, c2t, c3t, c4t, 
     $     full, order)
      implicit none
#include "errquit.fh"
#include "schwarz.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "inp.fh"
c
      integer rtdb
      integer basis             ! AO basis handle
      integer n1, n2, n3, n4    ! Dimension of each MO set
      integer ld1, ld2, ld3, ld4
      double precision c1t(ld1,*), c2t(ld2,*), ! Transposed MO coeffs
     $     c3t(ld3,*), c4t(ld4,*)
      double precision full(n1,n2,n3,n4)
      character*(*) order
c
c     Generate the specified block of MO integrals with
c     no assumptions of equivalence between the sets of coefficients.
c
c     Order can be either 
c     .    ChargeCloud -> full(p,q,r,s) = (pq|rs)
c     or 
c     .          Dirac -> full(p,q,r,s) = <pq|rs>
c     or
c     .  LeftAsymDirac -> full(p,q,r,s) = <pq|rs>-<qp|rs>
c     .  (must have c1t=c2t, n1=n2)
c     or
c     . RightAsymDirac -> full(p,q,r,s) = <pq|rs>-<pq|sr>
c     .  (must have c3t=c4t, n3=n4)
c
c     Presently the antisymmetrization is done at the top level
c     and the storage of full is not reduced to use the symmetry.
c
c     Memory requirements are 
c     .  n1*n2*n3*n4 + S*n2*n3*n4 + S*S*n3*n4 + S*S*S*n4 + 
c     .  maxs2 + maxg2*(1 + 4*integer)
c
c     Index 4 is the first transformed so there is advantage in
c     making it the smallest range.
c
      double precision tol2e
      parameter (tol2e = 1d-12)
C
      integer nsh, maxbfsh, lenhalf, lenthird, geom
      integer l_half, k_half, l_third, k_third
      integer lsh, ksh, llo, lhi, klo, khi
      integer p, q, r, s
      logical ochargecloud, oasym
      character*8 side
c
      if (inp_compare(.false.,order,'chargecloud')) then
         ochargecloud = .true.
         oasym = .false.
         side  = ' '
      else if (inp_compare(.false.,order,'dirac')) then
         ochargecloud = .false.
         oasym = .false.
         side  = ' '
      else if (inp_compare(.false.,order,'leftasymdirac')) then
         ochargecloud = .false.
         oasym = .true.
         side  = 'left'
      else if (inp_compare(.false.,order,'rightasymdirac')) then
         ochargecloud = .false.
         oasym = .true.
         side  = 'right'
      else
         call errquit('rjh_full_trans: unkown integral option',0,
     &       INPUT_ERR)
      endif
c
c     Initialize integrals and Schwarz screening
c
      if (.not. bas_geom(basis, geom))
     $     call errquit('rjh_transform: basis ', basis, BASIS_ERR)
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
c      
      if (.not. bas_numcont(basis, nsh)) call errquit(
     $     'rjh_transform: bas_numcont', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,maxbfsh)) call errquit(
     $     'rjh_transform: bas_nbf_cn_max', basis, BASIS_ERR)
c
      lenhalf = n3*n4*maxbfsh**2
      lenthird= n2*n3*n4*maxbfsh
c
      if (.not. ma_push_get(mt_dbl,lenhalf,'half',l_half, k_half))
     $     call errquit('ma half', lenhalf, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenthird,'third',l_third, k_third))
     $     call errquit('ma third', lenthird, MA_ERR)
c
      call dfill(n1*n2*n3*n4, 0.0d0, full, 1)
      do ksh = 1, nsh
         if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $        call errquit('rjh_transform: bas_cn2bfr',basis, BASIS_ERR)
         call dfill(n2*n3*n4*(khi-klo+1), 0.0d0, dbl_mb(k_third), 1)
         do lsh = 1, nsh
            if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $           call errquit('rjh_transform: bas_cn2bfr',basis,
     &       BASIS_ERR)
            if (schwarz_shell(ksh,lsh)*schwarz_max()
     $           .gt. tol2e) then
c     
c     Make (rs|kl) all rs (indices 3 and 4) given shells k and l
c     
               call rjh_half_transform(basis, ksh, lsh, n3, n4,
     $              c3t, c4t, ld3, ld4,
     $              dbl_mb(k_half), ochargecloud, tol2e)

*               write(6,*) ' ksh, lsh ', ksh, lsh
*               call rjh_debug_print('half',
*     $              dbl_mb(k_half), n3, n4, khi-klo+1, lhi-llo+1)
               
c     
               call rjh_third_transform(llo, lhi, klo, khi, 
     $              n2, n3, n4, dbl_mb(k_half), dbl_mb(k_third), 
     $              c2t, ld2, tol2e)
            end if
         end do
*         write(6,*) ' lsh ', lsh
*         call rjh_debug_print('third',
*     $        dbl_mb(k_third), n2, n3, n4, lhi-llo+1)
         call rjh_final_transform(klo, khi, n1, n2, n3, n4, 
     $        dbl_mb(k_third), full, c1t, ld1, tol2e)
      end do
c
      if (oasym) call rjh_asym_trans(full,n1,n2,n3,n4,side)
c
      do s = 1, n4
         do r = 1, n3
            do q = 1, n2
               do p = 1, n1
                  if (abs(full(p,q,r,s)).lt.1d-10)
     $                 full(p,q,r,s) = 0.0d0
               end do
            end do
         end do
      end do
c
      if (.not. ma_pop_stack(l_third)) call errquit('ma third',0,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_half)) call errquit('ma half',0, MA_ERR)
c
      call schwarz_tidy()
      call int_terminate
c
      end
      subroutine rjh_asym_trans(full,n1,n2,n3,n4,side)
      implicit none
#include "errquit.fh"
c
      integer n1, n2, n3, n4
      double precision full(n1,n2,n3,n4)
      character*(*) side
c
      integer p, q, r, s
      double precision tmp
c
      if (side .eq. 'left') then
         if (n1 .ne. n2) call errquit('rjh_asym_trans: left', n1,
     &       UNKNOWN_ERR)
         do s = 1, n4
            do r = 1, n3
               do q = 1, n2
                  do p = 1, q
                     tmp = full(p,q,r,s) - full(q,p,r,s)
                     full(p,q,r,s) = tmp
                     full(q,p,r,s) =-tmp
                  enddo
               enddo
            enddo
         enddo
      else
         if (n3 .ne. n4) call errquit('rjh_asym_trans: right', n3,
     &       UNKNOWN_ERR)
         do s = 1, n4
            do r = 1, s
               do q = 1, n2
                  do p = 1, n1
                     tmp = full(p,q,r,s) - full(q,p,s,r)
                     full(p,q,r,s) = tmp
                     full(p,q,s,r) =-tmp
                  enddo
               enddo
            enddo
         enddo
      endif
c
      end
      subroutine rjh_final_transform(klo, khi, 
     $     n1, n2, n3, n4, third, full, c1t, ld1, tol2e)
      implicit none
c
      integer klo, khi, n1, n2, n3, n4, ld1
      double precision third(n2,n3,n4,klo:khi)
      double precision full(n1,n2,n3,n4)
      double precision c1t(ld1,*)
      double precision tol2e
c
      integer k, s, r, q, p
      double precision g
c
      do k = klo, khi
         do s = 1, n4
            do r = 1, n3
               do q = 1, n2
                  g = third(q,r,s,k)
                  if (abs(g) .gt. tol2e) then
                     do p = 1, n1
                        full(p,q,r,s) = full(p,q,r,s) + g*c1t(p,k)
                     end do
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine rjh_third_transform(llo, lhi, klo, khi, 
     $     n2, n3, n4, half, third, c2t, ld2, tol2e)
      implicit none
c
      integer llo, lhi, klo, khi, n2, n3, n4, ld2
      double precision half(n3,n4,klo:khi,llo:lhi)
      double precision third(n2,n3,n4,klo:khi)
      double precision c2t(ld2,*)
      double precision tol2e
c
      integer k, l, s, r, q
      double precision g
c
      do l = llo, lhi
         do k = klo, khi
            do s = 1, n4
               do r = 1, n3
                  g = half(r,s,k,l)
                  if (abs(g) .gt. tol2e) then
                     do q  = 1, n2
                        third(q,r,s,k) = third(q,r,s,k) + g*c2t(q,l)
                     end do
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine rjh_half_transform(basis, ksh, lsh, n1, n2, 
     $     c1t, c2t, ld1, ld2, half, ochargecloud, tol2e)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "mafdecls.fh"
      integer basis, ksh, lsh, n1, n2, ld1, ld2
      double precision c1t(*), c2t(*) ! Transposed MO coeffs
      double precision half(*)  ! n1*n2*kdim*ldim (pq|kl)
      logical ochargecloud
      double precision tol2e
c
c     For a pair of shells k and l fill 
c
c     .   half(p,q,k,l) = (pq|kl)
c
c     for all p, q, and k, l within their respective shells
c     where p and q are transformed into the new bases and
c     k and l are AO indices.
c
c     Eventually exploiting sparsity and abelian symmetry
c     ... should also eventually use the texas integrals
c
c     Assumes that integrals and schwarz have been initialized.
c
      integer nbf, nsh, k_i, k_j, k_k, k_l, l_i, l_j, l_k, l_l,
     $     maxg2, maxs2, k_buf, l_buf, k_scr, l_scr, maxbfsh
      integer llo, lhi, klo, khi, lenaobuf, l_aobuf, k_aobuf
c
c     Get dimensions and required scratch space info
c
      if (.not. bas_numcont(basis, nsh)) call errquit(
     $     'rjh_transform: bas_numcont', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,maxbfsh)) call errquit(
     $     'rjh_transform: bas_nbf_cn_max', basis, BASIS_ERR)
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('rjh_transform: nbf',basis, BASIS_ERR)
      if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $     call errquit('rjh_transform: bas_cn2bfr',basis, BASIS_ERR)
      if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $     call errquit('rjh_transform: bas_cn2bfr',basis, BASIS_ERR)
      call int_mem_2e4c(maxg2,maxs2)
      lenaobuf = (khi-klo+1)*(lhi-llo+1)*maxbfsh*n2 ! (iq|kl)
c
c     Allocate scratch space for integrals and buffers for
c     transformation
c
      if (.not. ma_push_get(mt_dbl,maxs2,'scr',l_scr, k_scr))
     $   call errquit('ma scr',maxs2, MA_ERR)
      if (.not. ma_push_get(mt_dbl,maxg2,'buf',l_buf, k_buf))
     $   call errquit('ma buf',maxg2, MA_ERR)
      if (.not. ma_push_get(mt_int,maxg2,'i',l_i, k_i))
     $   call errquit('ma ibuf',maxg2, MA_ERR)
      if (.not. ma_push_get(mt_int,maxg2,'j',l_j, k_j))
     $   call errquit('ma jbuf',maxg2, MA_ERR)
      if (.not. ma_push_get(mt_int,maxg2,'k',l_k, k_k))
     $   call errquit('ma kbuf',maxg2, MA_ERR)
      if (.not. ma_push_get(mt_int,maxg2,'l',l_l, k_l))
     $   call errquit('ma lbuf',maxg2, MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenaobuf,'aobuf',l_aobuf, k_aobuf))
     $   call errquit('ma aobuf',lenaobuf, MA_ERR)
c
      call rjh_do_half_transform(
     $     basis, 
     $     dbl_mb(k_buf), dbl_mb(k_scr),
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     maxg2, maxs2, 
     $     nbf, nsh, maxbfsh, n1, n2, tol2e,
     $     ksh, lsh, klo, khi, llo, lhi,
     $     half, c1t, c2t, ld1, ld2, dbl_mb(k_aobuf), ochargecloud)
c
      if (.not. ma_chop_stack(l_scr)) call errquit
     $     ('rjh_transform: chopping stack', 0, MA_ERR)
c
      end
      subroutine rjh_do_half_transform(
     $     basis, 
     $     buf, scr, ilab, jlab, klab, llab, maxg2, maxs2, 
     $     nbf, nsh, maxbfsh, n1, n2, tol2e, 
     $     ksh, lsh, klo, khi, llo, lhi,
     $     half, c1t, c2t, ld1, ld2, aobuf, ochargecloud)
      implicit none
#include "errquit.fh"
#include "schwarz.fh"
#include "bas.fh"
      integer basis
      integer maxg2, maxs2
      integer nbf, nsh, n1, n2, maxbfsh, ld1, ld2
      double precision buf(maxg2), scr(maxs2)
      double precision c1t(ld1,nbf), c2t(ld2,nbf)
      integer ilab(maxg2), jlab(maxg2), klab(maxg2), llab(maxg2) 
      integer ksh, lsh, klo, khi, llo, lhi
      double precision aobuf(n2,maxbfsh,klo:khi,llo:lhi)
      double precision half(n1,n2,klo:khi,llo:lhi)
      double precision tol2e
      logical ochargecloud
c
c     For a pair of shells k and l, fill aobuf with integrals (ij|kl) 
c     for all i>=j ... eventually exploiting sparsity and abelian symmetry
c     ... should also use the texas integrals
c
      double precision skl, g
      integer ish, jsh, i, j, k, l, p, q, ijkl, nint, ilo, ihi, 
     $     idim, kdim, ldim
c
      kdim = khi - klo + 1
      ldim = lhi - llo + 1
      skl = schwarz_shell(ksh,lsh)
      call dfill(n1*n2*kdim*ldim,0.0d0,half,1)
c
      do ish = 1, nsh
         if (.not. bas_cn2bfr(basis,ish,ilo,ihi)) 
     $        call errquit('rjh_do_half_transform',ish, BASIS_ERR)
         idim = ihi-ilo+1
         do l = llo,lhi
            do k = klo, khi
               do i = 1, idim
                  do q = 1, n2
                     aobuf(q,i,k,l) = 0.0d0
                  end do
               end do
            end do
         end do
c
         do jsh = 1, nsh
            if (ochargecloud) then ! (ij|kl)
               if (schwarz_shell(ish,jsh)*skl .gt. tol2e) then
                  call int_l2e4c(basis, ish, jsh, basis, ksh, lsh,
     &                 tol2e, .false., maxg2, buf, nint, 
     $                 ilab, jlab, klab, llab, maxs2, scr)
                  do ijkl = 1, nint
                     i = ilab(ijkl)-ilo+1
                     j = jlab(ijkl)
                     k = klab(ijkl)
                     l = llab(ijkl)
                     g = buf(ijkl)
                     if (abs(g) .gt. tol2e) then
                        do q = 1, n2
                           aobuf(q,i,k,l) = aobuf(q,i,k,l) + g*c2t(q,j)
                        end do
                     end if
                  end do
               end if
            else                ! <ij|kl> = (ik|jl)
               if (schwarz_shell(ish,ksh)*schwarz_shell(jsh,lsh) 
     $              .gt. tol2e) then
                  call int_l2e4c(basis, ish, ksh, basis, jsh, lsh,
     &                 tol2e, .false., maxg2, buf, nint, 
     $                 ilab, klab, jlab, llab, maxs2, scr)
                  do ijkl = 1, nint
                     i = ilab(ijkl)-ilo+1
                     j = jlab(ijkl)
                     k = klab(ijkl)
                     l = llab(ijkl)
                     g = buf(ijkl)
                     if (abs(g) .gt. tol2e) then
                        do q = 1, n2
                           aobuf(q,i,k,l) = aobuf(q,i,k,l) + g*c2t(q,j)
                        end do
                     end if
                  end do
               end if
            endif
         end do
         do l = llo, lhi
            do k = klo, khi
               do i = ilo, ihi
                  do q = 1, n2
                     g = aobuf(q,i-ilo+1,k,l)
                     if (abs(g) .gt. tol2e) then
                        do p = 1, n1
                           half(p,q,k,l) = half(p,q,k,l) + g*c1t(p,i)
                        end do
                     end if
                  end do
               end do
            end do
         end do
      end do
c     
      end
