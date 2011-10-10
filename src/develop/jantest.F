      logical function jantest(rtdb)
*     
*     $$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
      integer rtdb
c
      integer basis, geom, nbf
      integer l_aoint, k_aoint        ! nbf**4 array of AO integrals
      integer l_moint, k_moint        ! nbf**4 array of AO integrals
      character*255 movecs      ! Name of movector file
      character*80 title, name_of_basis, scftype
      integer nbf_file, nsets, nmo_file(2)
      logical movecs_read, movecs_read_header
      external movecs_read, movecs_read_header
c
      integer nmo, g_tmp, l_occ, k_occ,
     $     l_eval, k_eval, l_mos, k_mos, l_most, k_most
      integer nocc, nvirt, nopen, nclosed, nso, noso
      integer l_t, k_t, l_t2, k_t2, l_f, k_f
      integer l_t1, k_t1, l_r1, k_r1
      integer l_tg, k_tg, l_tf, k_tf
      integer l_w, k_w, l_v, k_v
      logical int_normalize
      external int_normalize
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
      if (.not.int_normalize(rtdb,basis))
     $           call errquit('scf:int_normalize failed', 0, INT_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0, BASIS_ERR)
c
c     Read the MO vectors and evals from a RHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('jantest: failed to read movecs header',911, DISK_ERR)
      write(6,*) ' Read movecs header from ', movecs
      write(6,*) ' Job title :                ', 
     $     title(1:inp_strlen(title))
      write(6,*) ' Basis name:                ', 
     $     name_of_basis(1:inp_strlen(name_of_basis))
      nmo = nmo_file(1)
      if (.not. rtdb_get(rtdb, 'scf:nclosed', mt_int, 1, nocc))
     $     call errquit('nocc?',0, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nopen', mt_int, 1, nopen))
     $     call errquit('nopen?',0, RTDB_ERR)
      if (nopen .ne. 0) call errquit('asjdlfkadjsl',0, UNKNOWN_ERR)
      nvirt= nmo - nocc
      write(6,*) ' No. of closed shells       ', nocc
      write(6,*) ' No. of molecular orbitals: ', nmo
      write(6,*) ' No. of basis functions:    ', nbf
c
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 0, 0, g_tmp))
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'eval',l_eval, k_eval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'mos', l_mos, k_mos))
     $     call errquit('ma mos', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'mos', l_most, k_most))
     $     call errquit('ma mos', nbf*nbf, MA_ERR)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_eval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_mos), nbf)
      call util_transpose(dbl_mb(k_mos),nbf,dbl_mb(k_most),nmo,
     $     nbf,nmo)
c
      write(6,*) ' Orbital eigenvalues '
      call output(dbl_mb(k_eval),1,nmo,1,1,nmo,1,1)
      write(6,*) ' MOs'
      call output(dbl_mb(k_mos),1,nbf,1,nmo,nbf,nmo,1)
      write(6,*) ' MOs T'
      call output(dbl_mb(k_most),1,nmo,1,nbf,nmo,nbf,1)
c
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0, GA_ERR)
c
c     Make all AO integrals
c
      if (.not. ma_push_get(mt_dbl,nbf**4,'aoint',l_aoint,k_aoint))
     $     call errquit('allocation of AO integrals failed',nbf**4,
     &       MA_ERR)
      call jan_all_ao_integrals(rtdb,basis,nbf,'dirac',dbl_mb(k_aoint))
c     call jan_debug_print('AOINTS',dbl_mb(k_aoint), nbf,  nbf,  nbf,
c    $     nbf)
c
c     Make all MO integrals in Dirac order
c
      if (.not. ma_push_get(mt_dbl,nmo**4,'moint',l_moint,k_moint))
     $     call errquit('allocation of MO integrals failed',nbf**4,
     &       MA_ERR)
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most), 
     $     dbl_mb(k_moint), 'Dirac')
c      call jan_debug_print('MOINTS',dbl_mb(k_moint), nmo,  nmo,  nmo,
c     $     nmo)
c
c     do some incore cc
c     set nso=2*nmo, noso=2*nocc
c
      nso=2*nmo
      noso=2*nocc
      if (.not. ma_push_get(mt_dbl,nso**4,'t amps',l_t,k_t))
     $     call errquit('allocation of t amplitudes failed',nso**4,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,nso**4,'t2 amps',l_t2,k_t2))
     $     call errquit('allocation of t2 amplitudes failed',nso**4,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,nmo**2,'Fock Matrix',l_f,k_f))
     $     call errquit('allocation of t2 amplitudes failed',nmo**2,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,nso**2,'t1 amps',l_t1,k_t1))
     $     call errquit('allocation of t1 amplitudes failed',nso**2,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,nso**2,'t1 resi',l_r1,k_r1))
     $     call errquit('allocation of t1 residual failed',nso**2,
     &       MA_ERR)
c
      call ccsd_incore(rtdb, basis, dbl_mb(k_moint), 
     &                 dbl_mb(k_eval), dbl_mb(k_t), dbl_mb(k_t2),
     &                 dbl_mb(k_f), dbl_mb(k_t1), dbl_mb(k_r1),
     &                 nbf, nmo, nocc, nso, noso)
c
      if (.not. ma_push_get(mt_dbl,noso**3*nso**3,'t3 tg',l_tg,k_tg))
     $     call errquit('allocation of t3 tg failed',nso**6,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,noso**3*nso**3,'t3 tf',l_tg,k_tf))
     $     call errquit('allocation of t3 tf failed',nso**6, MA_ERR)
      if (.not. ma_push_get(mt_dbl,noso**3*nso**3,'t3 w',l_tg,k_w))
     $     call errquit('allocation of t3 w failed',nso**6, MA_ERR)
      if (.not. ma_push_get(mt_dbl,noso**3*nso**3,'t3 v',l_tg,k_v))
     $     call errquit('allocation of t3 v failed',nso**6, MA_ERR)
c
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 0, 0, g_tmp))
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_eval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_mos), nbf)
      call util_transpose(dbl_mb(k_mos),nbf,dbl_mb(k_most),nmo,
     $     nbf,nmo)
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0, GA_ERR)
c
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most), 
     $     dbl_mb(k_moint), 'Dirac')
c
      call triples_incore(rtdb, basis,
     &                 dbl_mb(k_moint), 
     &                 dbl_mb(k_eval), dbl_mb(k_t2),
     &                 dbl_mb(k_t1),
     &                 dbl_mb(k_tg), dbl_mb(k_tf),
     &                 dbl_mb(k_w), dbl_mb(k_v),
     &                 nbf, nmo, nocc, nso, noso)
c
c     Tidy up
c
      if (.not. ma_chop_stack(l_occ)) call errquit(' ma chop?', 0,
     &       MA_ERR)
      if (.not. bas_destroy(basis)) call errquit(' bas ?',0, BASIS_ERR)
      if (.not. geom_destroy(geom)) call errquit(' geom ?',0, GEOM_ERR)
      jantest = .true.
c
      end
      subroutine jan_all_ao_integrals(rtdb, basis, nbf, order, ao)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
      integer rtdb, basis, nbf
      double precision ao(nbf,nbf,nbf,nbf)
      character*(*) order
c
      integer nsh, k_i, k_j, k_k, k_l, l_i, l_j, l_k, l_l,
     $     maxg2, maxs2, k_buf, l_buf, k_scr, l_scr
c
      call int_init(rtdb, 1, basis)
      if ( .not. bas_numcont(basis, nsh) ) call errquit(
     $     'ao_fock_2e: problem with call to bas_numcont', basis,
     &       BASIS_ERR)
      call int_mem_2e4c(maxg2,maxs2)
      if (.not. ma_push_get(mt_dbl,maxs2,'scr',l_scr, k_scr))
     $   call errquit('ma scr',maxg2, MA_ERR)
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
c
      call jan_do_all_ao_integrals(basis, dbl_mb(k_buf), dbl_mb(k_scr),
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     maxg2, maxs2, nbf, nsh, order, ao)
c
      if (.not. ma_chop_stack(l_scr)) call errquit('janallao: ma?',0,
     &       MA_ERR)
c
      call int_terminate
c
      end
      subroutine jan_do_all_ao_integrals(
     $     basis, buf, scr, ilab, jlab, klab,
     $  llab, maxg2, maxs2, nbf, nsh, order, ao)
      implicit none
#include "errquit.fh"
c
      integer basis, nbf, nsh, maxg2, maxs2
      double precision buf(maxg2), scr(maxs2)
      integer ilab(maxg2), jlab(maxs2), klab(maxs2), llab(maxs2)
      integer i, j, k ,l, ish, jsh, ksh, lsh, ijkl, nint
      character*(*) order
      double precision ao(nbf,nbf,nbf,nbf)
      double precision zerotol
      logical omulliken
c
      omulliken = .false.   ! avoids compiler warning
      if (order .eq. 'mulliken') then
         omulliken = .true.
      else if (order .eq. 'dirac') then
         omulliken = .false.
      else
         call errquit(' unknown order',0, UNKNOWN_ERR)
      end if
c
      call dfill(nbf**4, 0.0d0, ao, 1)
      zerotol = 1d-12
c  
      do ish = 1, nsh
         do jsh = 1, nsh
            do ksh = 1, nsh
               do lsh = 1, nsh
                  call int_l2e4c(basis, ish, jsh, basis, ksh, lsh,
     &                 zerotol, .false., maxg2, buf, nint, 
     $                 ilab, jlab, klab, llab, maxs2, scr)
                  do ijkl = 1, nint
                     i = ilab(ijkl)
                     j = jlab(ijkl)
                     k = klab(ijkl)
                     l = llab(ijkl)
                     if (omulliken) then
                        ao(i,j,k,l) = buf(ijkl)
                     else
                        ao(i,k,j,l) = buf(ijkl)
                     end if
                  end do
               end do
            end do
         end do
      end do
c
c$$$      write(6,*)
c$$$      write(6,*) ' AO integrals '
c$$$      write(6,*)
c$$$      do i = 1, nbf
c$$$         do j = 1, nbf
c$$$            do k = 1, nbf
c$$$               do l = 1, nbf
c$$$                  if ( abs(ao(i,j,k,l)) .gt. 1e-6 ) 
c$$$     $                 write(6,7) i,j,k,l,ao(i,j,k,l)
c$$$ 7                format(1x,4i5,2x,f12.6)
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c
      end
      subroutine jan_full_transform(
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
      ochargecloud = .false.   ! avoids compiler warning
      oasym = .false.          ! avoids compiler warning
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
         call errquit('jan_full_trans: unkown integral option',0,
     &       INT_ERR)
      endif
c
c     Initialize integrals and Schwarz screening
c
      if (.not. bas_geom(basis, geom))
     $     call errquit('jan_transform: basis ', basis, BASIS_ERR)
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
c      
      if (.not. bas_numcont(basis, nsh)) call errquit(
     $     'jan_transform: bas_numcont', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,maxbfsh)) call errquit(
     $     'jan_transform: bas_nbf_cn_max', basis, BASIS_ERR)
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
     $        call errquit('jan_transform: bas_cn2bfr',basis, BASIS_ERR)
         call dfill(n2*n3*n4*(khi-klo+1), 0.0d0, dbl_mb(k_third), 1)
         do lsh = 1, nsh
            if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $           call errquit('jan_transform: bas_cn2bfr',basis,
     &       BASIS_ERR)
            if (schwarz_shell(ksh,lsh)*schwarz_max()
     $           .gt. tol2e) then
c     
c     Make (rs|kl) all rs (indices 3 and 4) given shells k and l
c     
               call jan_half_transform(basis, ksh, lsh, n3, n4,
     $              c3t, c4t, ld3, ld4,
     $              dbl_mb(k_half), ochargecloud, tol2e)

*               write(6,*) ' ksh, lsh ', ksh, lsh
*               call jan_debug_print('half',
*     $              dbl_mb(k_half), n3, n4, khi-klo+1, lhi-llo+1)
               
c     
               call jan_third_transform(llo, lhi, klo, khi, 
     $              n2, n3, n4, dbl_mb(k_half), dbl_mb(k_third), 
     $              c2t, ld2, tol2e)
            end if
         end do
*         write(6,*) ' lsh ', lsh
*         call jan_debug_print('third',
*     $        dbl_mb(k_third), n2, n3, n4, lhi-llo+1)
         call jan_final_transform(klo, khi, n1, n2, n3, n4, 
     $        dbl_mb(k_third), full, c1t, ld1, tol2e)
      end do
c
      if (oasym) call jan_asym_trans(full,n1,n2,n3,n4,side)
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
      if (.not. ma_pop_stack(l_half)) call errquit('ma half',0,
     &       MA_ERR)
c
      call schwarz_tidy()
      call int_terminate
c
      end
      subroutine jan_full_transform_noinit(
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
      ochargecloud = .false.   ! avoids compiler warning
      oasym = .false.          ! avoids compiler warning
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
         call errquit('jan_full_trans: unkown integral option',0,
     &       INT_ERR)
      endif
c
c     Initialize integrals and Schwarz screening
c
      if (.not. bas_geom(basis, geom))
     $     call errquit('jan_transform: basis ', basis, BASIS_ERR)
*      call int_init(rtdb, 1, basis)
*      call schwarz_init(geom, basis)
c      
      if (.not. bas_numcont(basis, nsh)) call errquit(
     $     'jan_transform: bas_numcont', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,maxbfsh)) call errquit(
     $     'jan_transform: bas_nbf_cn_max', basis, BASIS_ERR)
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
     $        call errquit('jan_transform: bas_cn2bfr',basis, BASIS_ERR)
         call dfill(n2*n3*n4*(khi-klo+1), 0.0d0, dbl_mb(k_third), 1)
         do lsh = 1, nsh
            if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $           call errquit('jan_transform: bas_cn2bfr',basis,
     &       BASIS_ERR)
            if (schwarz_shell(ksh,lsh)*schwarz_max()
     $           .gt. tol2e) then
c     
c     Make (rs|kl) all rs (indices 3 and 4) given shells k and l
c     
               call jan_half_transform(basis, ksh, lsh, n3, n4,
     $              c3t, c4t, ld3, ld4,
     $              dbl_mb(k_half), ochargecloud, tol2e)

*               write(6,*) ' ksh, lsh ', ksh, lsh
*               call jan_debug_print('half',
*     $              dbl_mb(k_half), n3, n4, khi-klo+1, lhi-llo+1)
               
c     
               call jan_third_transform(llo, lhi, klo, khi, 
     $              n2, n3, n4, dbl_mb(k_half), dbl_mb(k_third), 
     $              c2t, ld2, tol2e)
            end if
         end do
*         write(6,*) ' lsh ', lsh
*         call jan_debug_print('third',
*     $        dbl_mb(k_third), n2, n3, n4, lhi-llo+1)
         call jan_final_transform(klo, khi, n1, n2, n3, n4, 
     $        dbl_mb(k_third), full, c1t, ld1, tol2e)
      end do
c
      if (oasym) call jan_asym_trans(full,n1,n2,n3,n4,side)
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
      if (.not. ma_pop_stack(l_half)) call errquit('ma half',0,
     &       MA_ERR)
c
*      call schwarz_tidy()
*      call int_terminate
c
      end
      subroutine jan_asym_trans(full,n1,n2,n3,n4,side)
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
         if (n1 .ne. n2) call errquit('jan_asym_trans: left', n1,
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
         if (n3 .ne. n4) call errquit('jan_asym_trans: right', n3,
     &       UNKNOWN_ERR)
         do s = 1, n4
            do r = 1, s
               do q = 1, n2
                  do p = 1, n1
                     full(p,q,r,s) = full(p,q,r,s) - full(p,q,s,r)
                  enddo
               enddo
               if (r .ne. s) then
                  do q = 1, n2
                     do p = 1, n1
                        full(p,q,s,r) = -full(p,q,r,s)
                     enddo
                  enddo
               else
                  do q = 1, n2
                     do p = 1, n1
                        full(p,q,s,r) = 0.0d0
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif
c
      end
      subroutine jan_final_transform(klo, khi, 
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
      subroutine jan_third_transform(llo, lhi, klo, khi, 
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
      subroutine jan_half_transform(basis, ksh, lsh, n1, n2, 
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
     $     'jan_transform: bas_numcont', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,maxbfsh)) call errquit(
     $     'jan_transform: bas_nbf_cn_max', basis, BASIS_ERR)
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('jan_transform: nbf',basis, BASIS_ERR)
      if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $     call errquit('jan_transform: bas_cn2bfr',basis, BASIS_ERR)
      if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $     call errquit('jan_transform: bas_cn2bfr',basis, BASIS_ERR)
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
      call jan_do_half_transform(
     $     basis, 
     $     dbl_mb(k_buf), dbl_mb(k_scr),
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     maxg2, maxs2, 
     $     nbf, nsh, maxbfsh, n1, n2, tol2e,
     $     ksh, lsh, klo, khi, llo, lhi,
     $     half, c1t, c2t, ld1, ld2, dbl_mb(k_aobuf), ochargecloud)
c
      if (.not. ma_chop_stack(l_scr)) call errquit
     $     ('jan_transform: chopping stack', 0, MA_ERR)
c
      end
      subroutine jan_do_half_transform(
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
     $        call errquit('jan_do_half_transform',ish, BASIS_ERR)
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
      subroutine jan_debug_print(string,full, n1, n2, n3, n4)
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
 1                   format(4i5,2x,f14.8)
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine ccsd_incore(rtdb, basis, g_mo, e, r2, t2, fock, 
     &                       t1, r1, nbf, nmo, nocc, nso, noso)
      implicit none
      integer rtdb, basis, nbf, nmo, nocc, nso, noso
c
c     nbf = number of basis functions
c     nmo = number of molecular orbitals
c     nocc = number of occupied orbitals
c     nso = number of spin orbitals (for now set to 2*nmo)
c     noso = number of occupied spin orbitals (for now set to 2*nocc)
c
c     toy coupled cluster program based on the green book
c
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision e(nmo), fock(nmo,nmo)
      double precision r2(nso,nso,nso,nso)
      double precision t2(nso,nso,nso,nso)
      double precision t1(nso,nso)
      double precision r1(nso,nso)
      integer n
c
c     cc stuff (greek indices denote occupied spin-orbitals)
c     cc stuff (r, s, t, u, ... denote unoccupied spin-orbitals)
c
      write(6,*)' nbf, nmo, nocc, nso, noso: ',
     &            nbf, nmo, nocc, nso, noso
c
c     Initially fill the Fock matrix with diagonal elements only 
c
      call dfill(nmo**2, 0.0d0, fock, 1) 
      do n = 1, nmo
        fock(n,n) = e(n)
      enddo
c      write(6,*) ' Fock Matrix'
c      call output(fock,1,nmo,1,nmo,nmo,nmo,1)
c
c     First generate an initial guess for t2
c
      call dfill(nso**4, 0.0d0, r2, 1) 
      call dfill(nso**4, 0.0d0, t2, 1) 
      call dfill(nso**2, 0.0d0, r1, 1) 
      call dfill(nso**2, 0.0d0, t1, 1) 
c
c      call t2_init(g_mo, e, t, nmo, nso, noso)
c      write(6,*)' T2 initial guess. '
c      call writet2(t,nso)
c      call correlation(g_mo,t,nmo,nso,noso)
c
c
c     next put this guess into t2 expression keeping all terms linear in t2
c
c      call t2_l_gb(g_mo, e, t, t2, nmo, nso, noso)
c      do n = 1, 3
c      call dfill(nso**4, 0.0d0, t2, 1) 
c      call t2_l_rjh(g_mo, t, t2, fock, nmo, nso, noso)
c      call t2_update_rjh(g_mo, t, t2, fock, nmo, nso, noso)
c      write(6,*)' T2 with linear terms. '
c      call writet2(t,nso)
c      call correlation(g_mo,t,nmo,nso,noso)
c      enddo
c
c     next put this into t2 expression keeping all terms
c
      do n = 1, 50
         call dfill(nso**4, 0.0d0, r2, 1) 
         call dfill(nso**2, 0.0d0, r1, 1) 
c         call t2_l_q_gb(g_mo, e, t, t2, nmo, nso, noso)
c         call t2_l_rjh(g_mo, t, t2, fock, nmo, nso, noso)
c
c        Generate delta T1 
c
         call t1_rjh(g_mo, r1, t2, fock, nmo, nso, noso)
c
c        Generate delta T2
c
         call t2_l_q_rjh(g_mo, t2, r2, fock, nmo, nso, noso)
c
c        Update T1
c
         call t1_update_rjh(g_mo, t1, r1, fock, nmo, nso, noso)
c         write(6,*)' T1: '
c         call writet1(t1,nso)
c
c        Update T2
c
         call t2_update_rjh(g_mo, t2, r2, fock, nmo, nso, noso)
c         write(6,*)' T2 with quadratic terms. '
c         call writet2(t2,nso)
c
c        Compute Energy expression
c
         call correlation(g_mo,t1,t2,nmo,nso,noso)
c
c        Use T1 to transform 1e- and 2e- ints
c
         call get_new_f_g(rtdb, basis, fock, g_mo, t1, nso, nmo)
      enddo
c
      return
      end
      subroutine triples_incore(rtdb, basis, g_mo, ener, t2,
     &                          t1, tg, tf, w, v, 
     &                          nbf, nmo, nocc, nso, noso)
      implicit none
      integer rtdb, basis, nbf, nmo, nocc, nso, noso
c
c     nbf = number of basis functions
c     nmo = number of molecular orbitals
c     nocc = number of occupied orbitals
c     nso = number of spin orbitals (for now set to 2*nmo)
c     noso = number of occupied spin orbitals (for now set to 2*nocc)
c
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision ener(nmo)
      double precision t2(nso,nso,nso,nso)
      double precision t1(nso,nso)
      double precision tg(nso,nso,nso,noso,noso,noso)
      double precision tf(nso,nso,nso,noso,noso,noso)
      double precision w(nso,nso,nso,noso,noso,noso)
      double precision v(nso,nso,nso,noso,noso,noso)
      double precision g, et
c occupied orbitals
      integer i, j, k, l, m, n
      integer i_orb, j_orb, k_orb, l_orb, m_orb, n_orb
      integer i_spin, j_spin, k_spin, l_spin, m_spin, n_spin
c unoccupied orbitals
      integer a, b, c, d, e, f
      integer a_orb, b_orb, c_orb, d_orb, e_orb, f_orb
      integer a_spin, b_spin, c_spin, d_spin, e_spin, f_spin
c
      write(6,*)' nbf, nmo, nocc, nso, noso: ',
     &            nbf, nmo, nocc, nso, noso
c
      call dfill(noso**3*nso**3, 0.0d0, tg, 1)
      call dfill(noso**3*nso**3, 0.0d0, tf, 1)
      call dfill(noso**3*nso**3, 0.0d0, w, 1)
      call dfill(noso**3*nso**3, 0.0d0, v, 1)
c
      do a = noso+1, nso
         a_orb = (a-1)/2 + 1
         a_spin = 0
         if ((a-2*a_orb).ne.0)a_spin = 1
         do b = noso+1, nso
            b_orb = (b-1)/2 + 1
            b_spin = 0
            if ((b-2*b_orb).ne.0)b_spin = 1
            do c = noso+1, nso
               c_orb = (c-1)/2 + 1
               c_spin = 0
               if ((c-2*c_orb).ne.0)c_spin = 1
c
               do i = 1, noso
                  i_orb = (i-1)/2 + 1
                  i_spin = 0
                  if ((i-2*i_orb).ne.0)i_spin = 1
                  do j = 1, noso
                     j_orb = (j-1)/2 + 1
                     j_spin = 0
                     if ((j-2*j_orb).ne.0)j_spin = 1
                     do k = 1, noso
                        k_orb = (k-1)/2 + 1
                        k_spin = 0
                        if ((k-2*k_orb).ne.0)k_spin = 1
c
                        do e = noso+1, nso
                           e_orb = (e-1)/2 + 1
                           e_spin = 0
                           if ((e-2*e_orb).ne.0)e_spin = 1
c
                           g=0.0d0
                           if (a_spin.eq.e_spin.and.
     &                         b_spin.eq.k_spin)then
                              g = g_mo(a_orb,b_orb,e_orb,k_orb)
                           endif
                           if (a_spin.eq.k_spin.and.
     &                         b_spin.eq.e_spin)then
                              g = g - g_mo(a_orb,b_orb,k_orb,e_orb)
                           endif
c                           write(6,*)' g:',g
                           tg(a,b,c,i,j,k) = tg(a,b,c,i,j,k) +
     &                        g*t2(c,e,i,j)
                        enddo
c
                        do m = 1, noso
                           m_orb = (m-1)/2 + 1
                           m_spin = 0
                           if ((m-2*m_orb).ne.0)m_spin = 1
c
                           g=0.0d0
                           if (c_spin.eq.i_spin.and.
     &                         m_spin.eq.j_spin)then
                              g = g_mo(c_orb,m_orb,i_orb,j_orb)
                           endif
                           if (c_spin.eq.j_spin.and.
     &                         m_spin.eq.i_spin)then
                              g = g - g_mo(c_orb,m_orb,j_orb,i_orb)
                           endif
c                           write(6,*)' g:',g
                           tg(a,b,c,i,j,k) = tg(a,b,c,i,j,k) -
     &                        g*t2(a,b,m,k)
                        enddo
c
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.j_spin.eq.b_spin)then
                           g = g_mo(i_orb,j_orb,a_orb,b_orb)
                        endif
                        if (i_spin.eq.b_spin.and.j_spin.eq.a_spin)then
                           g = g - g_mo(i_orb,j_orb,b_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        tf(a,b,c,i,j,k) = tg(a,b,c,i,j,k) +
     &                     g*t1(c,k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      do a = noso+1, nso
         a_orb = (a-1)/2 + 1
         a_spin = 0
         if ((a-2*a_orb).ne.0)a_spin = 1
         do b = noso+1, nso
            b_orb = (b-1)/2 + 1
            b_spin = 0
            if ((b-2*b_orb).ne.0)b_spin = 1
            do c = noso+1, nso
               c_orb = (c-1)/2 + 1
               c_spin = 0
               if ((c-2*c_orb).ne.0)c_spin = 1
c
               do i = 1, noso
                  i_orb = (i-1)/2 + 1
                  i_spin = 0
                  if ((i-2*i_orb).ne.0)i_spin = 1
                  do j = 1, noso
                     j_orb = (j-1)/2 + 1
                     j_spin = 0
                     if ((j-2*j_orb).ne.0)j_spin = 1
                     do k = 1, noso
                        k_orb = (k-1)/2 + 1
                        k_spin = 0
                        if ((k-2*k_orb).ne.0)k_spin = 1
c
                        w(a,b,c,i,j,k) = w(a,b,c,i,j,k) +
     &                              tg(a,b,c,i,j,k) - 
     &                              tg(a,b,c,k,j,i) - 
     &                              tg(a,b,c,i,k,j) - 
     &                              tg(c,b,a,i,j,k) - 
     &                              tg(a,c,b,i,j,k) + 
     &                              tg(c,b,a,k,j,i) + 
     &                              tg(a,c,b,k,j,i) + 
     &                              tg(c,b,a,i,k,j) + 
     &                              tg(a,c,b,i,k,j)  
c
                        v(a,b,c,i,j,k) = v(a,b,c,i,j,k) +
     &                              tf(a,b,c,i,j,k) - 
     &                              tf(a,b,c,k,j,i) - 
     &                              tf(a,b,c,i,k,j) - 
     &                              tf(c,b,a,i,j,k) - 
     &                              tf(a,c,b,i,j,k) + 
     &                              tf(c,b,a,k,j,i) + 
     &                              tf(a,c,b,k,j,i) + 
     &                              tf(c,b,a,i,k,j) + 
     &                              tf(a,c,b,i,k,j)  
c                        
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      Et = 0.0d0
c
      do a = noso+1, nso
         a_orb = (a-1)/2 + 1
         a_spin = 0
         if ((a-2*a_orb).ne.0)a_spin = 1
         do b = noso+1, nso
            b_orb = (b-1)/2 + 1
            b_spin = 0
            if ((b-2*b_orb).ne.0)b_spin = 1
            do c = noso+1, nso
               c_orb = (c-1)/2 + 1
               c_spin = 0
               if ((c-2*c_orb).ne.0)c_spin = 1
c
               do i = 1, noso
                  i_orb = (i-1)/2 + 1
                  i_spin = 0
                  if ((i-2*i_orb).ne.0)i_spin = 1
                  do j = 1, noso
                     j_orb = (j-1)/2 + 1
                     j_spin = 0
                     if ((j-2*j_orb).ne.0)j_spin = 1
                     do k = 1, noso
                        k_orb = (k-1)/2 + 1
                        k_spin = 0
                        if ((k-2*k_orb).ne.0)k_spin = 1
c
                        Et = Et - 1.0d0/36.0d0*w(a,b,c,i,j,k)*
     &                                         v(a,b,c,i,j,k)/
     &                        ( ener(a_orb) +
     &                          ener(b_orb) +
     &                          ener(c_orb) -
     &                          ener(i_orb) -
     &                          ener(j_orb) -
     &                          ener(k_orb) )

                        
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      write(6,*)' (t) energy = ',et
c
c      do a = noso+1, nso
c         a_orb = (a-1)/2 + 1
c         a_spin = 0
c         if ((a-2*a_orb).ne.0)a_spin = 1
c         do b = noso+1, nso
c            b_orb = (b-1)/2 + 1
c            b_spin = 0
c            if ((b-2*b_orb).ne.0)b_spin = 1
c            do c = noso+1, nso
c               c_orb = (c-1)/2 + 1
c               c_spin = 0
c               if ((c-2*c_orb).ne.0)c_spin = 1
cc
c               do i = 1, noso
c                  i_orb = (i-1)/2 + 1
c                  i_spin = 0
c                  if ((i-2*i_orb).ne.0)i_spin = 1
c                  do j = 1, noso
c                     j_orb = (j-1)/2 + 1
c                     j_spin = 0
c                     if ((j-2*j_orb).ne.0)j_spin = 1
c                     do k = 1, noso
c                        k_orb = (k-1)/2 + 1
c                        k_spin = 0
c                        if ((k-2*k_orb).ne.0)k_spin = 1
cc
c                        if(dabs(w(a,b,c,i,j,k)).gt.1d-10)then
c                          write(6,*)' i, j, k, a, b, c, w, v: ',
c     &                                i, j, k, a, b, c, 
c     &                                w(a,b,c,i,j,k),
c     &                                v(a,b,c,i,j,k)
c                        endif
cc
c                     enddo
c                  enddo
c               enddo
c            enddo
c         enddo
c      enddo
c
      return
      end
      subroutine t2_init(g_mo, e, t, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision e(nmo)
      double precision t(nso,nso,nso,nso)
      double precision g, denom
c
c     cc stuff (greek indices denote occupied spin-orbitals)
c     cc stuff (r, s, t, u, ... denote unoccupied spin-orbitals)
c
c occupied orbitals
      integer alpha, beta
      integer alpha_orb, beta_orb
      integer alpha_spin, beta_spin
c unoccupied orbitals
      integer m, n
      integer m_orb, n_orb
      integer m_spin, n_spin
c
c     First generate an initial guess for t2
c
      call dfill(nso**4, 0.0d0, t, 1) 
      do m = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         m_orb = (m-1)/2 + 1
         m_spin = 0
         if ((m-2*m_orb).ne.0)m_spin = 1
c         write(6,*)' m, m_orb, m_spin: ',
c     &               m, m_orb, m_spin
         do n = noso+1, nso
            n_orb = (n-1)/2 + 1
            n_spin = 0
            if ((n-2*n_orb).ne.0)n_spin = 1
c            write(6,*)' n, n_orb, n_spin: ',
c     &                  n, n_orb, n_spin
            do alpha = 1, noso
               alpha_orb = (alpha-1)/2 + 1
               alpha_spin = 0
               if ((alpha-2*alpha_orb).ne.0)alpha_spin = 1
c               write(6,*)' alpha, alpha_orb, alpha_spin: ',
c     &                     alpha, alpha_orb, alpha_spin
               do beta = 1, noso
                  beta_orb = (beta-1)/2 + 1
                  beta_spin = 0
                  if ((beta-2*beta_orb).ne.0)beta_spin = 1
c                  write(6,*)' beta, beta_orb, beta_spin: ',
c     &                        beta, beta_orb, beta_spin
c
                  denom = -e(m_orb)-e(n_orb)+e(alpha_orb)+e(beta_orb)
c                  write(6,*)' denom: ', denom
                  g=0.0d0
                  if (m_spin.eq.alpha_spin.and.n_spin.eq.beta_spin)then
                     g = g_mo(m_orb,n_orb,alpha_orb,beta_orb)
                  endif
                  if (m_spin.eq.beta_spin.and.n_spin.eq.alpha_spin)then
                     g = g - g_mo(m_orb,n_orb,beta_orb,alpha_orb)
                  endif
c                  write(6,*)' g:',g
                  t(m,n,alpha,beta) = g/denom
c                  write(6,*)' m,n,alpha,beta,t(m,n,alpha,beta):',
c     &                        m,n,alpha,beta,t(m,n,alpha,beta)
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine t2_l_gb(g_mo, e, t, t2, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision e(nmo)
      double precision t(nso,nso,nso,nso)
      double precision t2(nso,nso,nso,nso)
      double precision g, denom
c
c     cc stuff (greek indices denote occupied spin-orbitals)
c     cc stuff (r, s, t, u, ... denote unoccupied spin-orbitals)
c
c occupied orbitals
      integer alpha, beta, delta, gamma
      integer alpha_orb, beta_orb, delta_orb, gamma_orb
      integer alpha_spin, beta_spin, delta_spin, gamma_spin
c unoccupied orbitals
      integer m, n, p, q
      integer m_orb, n_orb, p_orb, q_orb
      integer m_spin, n_spin, p_spin, q_spin
c
c     put t into t2 expression keeping all terms linear in t2
c
      do m = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         m_orb = (m-1)/2 + 1
         m_spin = 0
         if ((m-2*m_orb).ne.0)m_spin = 1
         do n = noso+1, nso
            n_orb = (n-1)/2 + 1
            n_spin = 0
            if ((n-2*n_orb).ne.0)n_spin = 1
            do alpha = 1, noso
               alpha_orb = (alpha-1)/2 + 1
               alpha_spin = 0
               if ((alpha-2*alpha_orb).ne.0)alpha_spin = 1
               do beta = 1, noso
                  beta_orb = (beta-1)/2 + 1
                  beta_spin = 0
                  if ((beta-2*beta_orb).ne.0)beta_spin = 1
                  denom = e(m_orb)+e(n_orb)-e(alpha_orb)-e(beta_orb)
c                  write(6,*)' denom: ', denom
                  g=0.0d0
                  if (m_spin.eq.alpha_spin.and.n_spin.eq.beta_spin)then
                     g = g_mo(m_orb,n_orb,alpha_orb,beta_orb)
                  endif
                  if (m_spin.eq.beta_spin.and.n_spin.eq.alpha_spin)then
                     g = g - g_mo(m_orb,n_orb,beta_orb,alpha_orb)
                  endif
c                  write(6,*)' g:',g
                  t2(m,n,alpha,beta) = g
c
                  do p = noso+1, nso
                     p_orb = (p-1)/2 + 1
                     p_spin = 0
                     if ((p-2*p_orb).ne.0)p_spin = 1
                     do q = noso+1, p-1
                        q_orb = (q-1)/2 + 1
                        q_spin = 0
                        if ((q-2*q_orb).ne.0)q_spin = 1
c
                        g=0.0d0
                        if (m_spin.eq.p_spin.and.n_spin.eq.q_spin)then
                           g = g_mo(m_orb,n_orb,p_orb,q_orb)
                        endif
                        if (m_spin.eq.q_spin.and.n_spin.eq.p_spin)then
                           g = g - g_mo(m_orb,n_orb,q_orb,p_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(p,q,alpha,beta)
c
                     enddo
                  enddo
c
                  do gamma = 1, noso
                     gamma_orb = (gamma-1)/2 + 1
                     gamma_spin = 0
                     if ((gamma-2*gamma_orb).ne.0)gamma_spin = 1
                     do delta = 1, gamma-1
                        delta_orb = (delta-1)/2 + 1
                        delta_spin = 0
                        if ((delta-2*delta_orb).ne.0)delta_spin = 1
c
                        g=0.0d0
                        if (gamma_spin.eq.alpha_spin.and.
     &                      delta_spin.eq.beta_spin)then
                           g = 
     &                     g_mo(gamma_orb,delta_orb,alpha_orb,beta_orb)
                        endif
                        if (gamma_spin.eq.beta_spin.and.
     &                      delta_spin.eq.alpha_spin)then
                           g = g - 
     &                     g_mo(gamma_orb,delta_orb,beta_orb,alpha_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(m,n,gamma,delta)
c
                     enddo
                  enddo
c
                  do p = noso+1, nso
                     p_orb = (p-1)/2 + 1
                     p_spin = 0
                     if ((p-2*p_orb).ne.0)p_spin = 1
                     do gamma = 1, noso
                        gamma_orb = (gamma-1)/2 + 1
                        gamma_spin = 0
                        if ((gamma-2*gamma_orb).ne.0)gamma_spin = 1
c
                        g=0.0d0
                        if (gamma_spin.eq.beta_spin.and.
     &                      n_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,n_orb,beta_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      n_spin.eq.beta_spin)then
                           g = g - g_mo(gamma_orb,n_orb,p_orb,beta_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) + 
     &                                       g*t(m,p,alpha,gamma)
c
                        g=0.0d0
                        if (gamma_spin.eq.beta_spin.and.
     &                      m_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,m_orb,beta_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      m_spin.eq.beta_spin)then
                           g = g - g_mo(gamma_orb,m_orb,p_orb,beta_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(n,p,alpha,gamma)
c
                        g=0.0d0
                        if (gamma_spin.eq.alpha_spin.and.
     &                      n_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,n_orb,alpha_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      n_spin.eq.alpha_spin)then
                           g = g - g_mo(gamma_orb,n_orb,p_orb,alpha_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(m,p,beta,gamma)
c
                        g=0.0d0
                        if (gamma_spin.eq.alpha_spin.and.
     &                      m_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,m_orb,alpha_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      m_spin.eq.alpha_spin)then
                           g = g - g_mo(gamma_orb,m_orb,p_orb,alpha_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) + 
     &                                       g*t(n,p,beta,gamma)
                     enddo
                  enddo
                  t2(m,n,alpha,beta) = t2(m,n,alpha,beta)/denom
c                  write(6,*)' t2(m,n,alpha,beta):',
c     &                        t2(m,n,alpha,beta)
               enddo
            enddo
         enddo
      enddo
c      write(6,*)' T2 with linear terms. '
c      call writet2(t2,nso)
c      call correlation(g_mo,t2,nmo,nso,noso)
      return
      end
      subroutine t2_l_q_gb(g_mo, e, t, t2, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision e(nmo)
      double precision t(nso,nso,nso,nso)
      double precision t2(nso,nso,nso,nso)
      double precision g, denom
c
c     cc stuff (greek indices denote occupied spin-orbitals)
c     cc stuff (r, s, t, u, ... denote unoccupied spin-orbitals)
c
c occupied orbitals
      integer alpha, beta, delta, gamma
      integer alpha_orb, beta_orb, delta_orb, gamma_orb
      integer alpha_spin, beta_spin, delta_spin, gamma_spin
c unoccupied orbitals
      integer m, n, p, q
      integer m_orb, n_orb, p_orb, q_orb
      integer m_spin, n_spin, p_spin, q_spin
c
c     put t into t2 expression keeping all terms (linear and quadratic)
c
      do m = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         m_orb = (m-1)/2 + 1
         m_spin = 0
         if ((m-2*m_orb).ne.0)m_spin = 1
         do n = noso+1, nso
            n_orb = (n-1)/2 + 1
            n_spin = 0
            if ((n-2*n_orb).ne.0)n_spin = 1
            do alpha = 1, noso
               alpha_orb = (alpha-1)/2 + 1
               alpha_spin = 0
               if ((alpha-2*alpha_orb).ne.0)alpha_spin = 1
               do beta = 1, noso
                  beta_orb = (beta-1)/2 + 1
                  beta_spin = 0
                  if ((beta-2*beta_orb).ne.0)beta_spin = 1
                  denom = e(m_orb)+e(n_orb)-e(alpha_orb)-e(beta_orb)
c                  write(6,*)' denom: ', denom
                  g=0.0d0
                  if (m_spin.eq.alpha_spin.and.n_spin.eq.beta_spin)then
                     g = g_mo(m_orb,n_orb,alpha_orb,beta_orb)
                  endif
                  if (m_spin.eq.beta_spin.and.n_spin.eq.alpha_spin)then
                     g = g - g_mo(m_orb,n_orb,beta_orb,alpha_orb)
                  endif
c                  write(6,*)' g:',g
                  t2(m,n,alpha,beta) = g
c
                  do p = noso+1, nso
                     p_orb = (p-1)/2 + 1
                     p_spin = 0
                     if ((p-2*p_orb).ne.0)p_spin = 1
                     do q = noso+1, p-1
                        q_orb = (q-1)/2 + 1
                        q_spin = 0
                        if ((q-2*q_orb).ne.0)q_spin = 1
c
                        g=0.0d0
                        if (m_spin.eq.p_spin.and.n_spin.eq.q_spin)then
                           g = g_mo(m_orb,n_orb,p_orb,q_orb)
                        endif
                        if (m_spin.eq.q_spin.and.n_spin.eq.p_spin)then
                           g = g - g_mo(m_orb,n_orb,q_orb,p_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(p,q,alpha,beta)
c
                     enddo
                  enddo
c
                  do gamma = 1, noso
                     gamma_orb = (gamma-1)/2 + 1
                     gamma_spin = 0
                     if ((gamma-2*gamma_orb).ne.0)gamma_spin = 1
                     do delta = 1, gamma-1
                        delta_orb = (delta-1)/2 + 1
                        delta_spin = 0
                        if ((delta-2*delta_orb).ne.0)delta_spin = 1
c
                        g=0.0d0
                        if (gamma_spin.eq.alpha_spin.and.
     &                      delta_spin.eq.beta_spin)then
                           g = 
     &                     g_mo(gamma_orb,delta_orb,alpha_orb,beta_orb)
                        endif
                        if (gamma_spin.eq.beta_spin.and.
     &                      delta_spin.eq.alpha_spin)then
                           g = g - 
     &                     g_mo(gamma_orb,delta_orb,beta_orb,alpha_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(m,n,gamma,delta)
c
                     enddo
                  enddo
c
                  do p = noso+1, nso
                     p_orb = (p-1)/2 + 1
                     p_spin = 0
                     if ((p-2*p_orb).ne.0)p_spin = 1
                     do gamma = 1, noso
                        gamma_orb = (gamma-1)/2 + 1
                        gamma_spin = 0
                        if ((gamma-2*gamma_orb).ne.0)gamma_spin = 1
c
                        g=0.0d0
                        if (gamma_spin.eq.beta_spin.and.
     &                      n_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,n_orb,beta_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      n_spin.eq.beta_spin)then
                           g = g - g_mo(gamma_orb,n_orb,p_orb,beta_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) + 
     &                                       g*t(m,p,alpha,gamma)
c
                        g=0.0d0
                        if (gamma_spin.eq.beta_spin.and.
     &                      m_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,m_orb,beta_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      m_spin.eq.beta_spin)then
                           g = g - g_mo(gamma_orb,m_orb,p_orb,beta_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(n,p,alpha,gamma)
c
                        g=0.0d0
                        if (gamma_spin.eq.alpha_spin.and.
     &                      n_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,n_orb,alpha_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      n_spin.eq.alpha_spin)then
                           g = g - g_mo(gamma_orb,n_orb,p_orb,alpha_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) - 
     &                                       g*t(m,p,beta,gamma)
c
                        g=0.0d0
                        if (gamma_spin.eq.alpha_spin.and.
     &                      m_spin.eq.p_spin)then
                           g = g_mo(gamma_orb,m_orb,alpha_orb,p_orb)
                        endif
                        if (gamma_spin.eq.p_spin.and.
     &                      m_spin.eq.alpha_spin)then
                           g = g - g_mo(gamma_orb,m_orb,p_orb,alpha_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) + 
     &                                       g*t(n,p,beta,gamma)
                     enddo
                  enddo
c                  
                  do gamma = 1, noso
                     gamma_orb = (gamma-1)/2 + 1
                     gamma_spin = 0
                     if ((gamma-2*gamma_orb).ne.0)gamma_spin = 1
                     do delta = 1, gamma-1
                        delta_orb = (delta-1)/2 + 1
                        delta_spin = 0
                        if ((delta-2*delta_orb).ne.0)delta_spin = 1
                        do p = noso+1, nso
                           p_orb = (p-1)/2 + 1
                           p_spin = 0
                           if ((p-2*p_orb).ne.0)p_spin = 1
                           do q = noso+1, p-1
                              q_orb = (q-1)/2 + 1
                              q_spin = 0
                              if ((q-2*q_orb).ne.0)q_spin = 1
c
                              g=0.0d0
                              if (gamma_spin.eq.p_spin.and.
     &                           delta_spin.eq.q_spin)then
                              g = g_mo(gamma_orb,delta_orb,p_orb,q_orb)
                              endif
                              if (gamma_spin.eq.q_spin.and.
     &                            delta_spin.eq.p_spin)then
                          g = g - g_mo(gamma_orb,delta_orb,q_orb,p_orb)
                              endif
c                              write(6,*)' g:',g
c                                 
                              t2(m,n,alpha,beta) = t2(m,n,alpha,beta) + 
     &                      g*(t(p,q,alpha,beta)*t(m,n,gamma,delta) - 
     &                  2.0d0*(t(m,p,alpha,beta)*t(n,q,gamma,delta) + 
     &                         t(n,q,alpha,beta)*t(m,p,gamma,delta))- 
     &                  2.0d0*(t(m,n,alpha,gamma)*t(p,q,beta,delta) + 
     &                         t(p,q,alpha,gamma)*t(m,n,beta,delta))+ 
     &                  4.0d0*(t(m,p,alpha,gamma)*t(n,q,beta,delta) + 
     &                         t(n,q,alpha,gamma)*t(m,p,beta,delta)))
                           enddo
                        enddo
                     enddo
                  enddo
                  t2(m,n,alpha,beta) = t2(m,n,alpha,beta)/denom
c                  write(6,*)' t2(m,n,alpha,beta):',
c     &                        t2(m,n,alpha,beta)
               enddo
            enddo
         enddo
      enddo
c      write(6,*)' T2 with all terms. '
c      call writet2(t2,nso)
c      call correlation(g_mo,t2,nmo,nso,noso)
      return
      end
      subroutine writet2(t,len)
      implicit none
      integer len, i, j, k, l
      double precision t(len,len,len,len)
      do i = 1, len
         do j = 1, len
            do k = 1, len
               do l = 1, len
                  if (abs(t(i,j,k,l)).gt.1d-6)
     &               write(6,*)i,j,k,l,t(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine writet1(t,len)
      implicit none
      integer len, i, j
      double precision t(len,len)
      do i = 1, len
         do j = 1, len
            if (abs(t(i,j)).gt.1d-6)
     &         write(6,*)i,j,t(i,j)
         enddo
      enddo
      return
      end
      subroutine correlation(g_mo,t1,t2,nmo,nso,noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision t2(nso,nso,nso,nso)
      double precision t1(nso,nso)
      double precision g, e2
c occupied orbitals
      integer alpha, beta
      integer alpha_orb, beta_orb
      integer alpha_spin, beta_spin
c unoccupied orbitals
      integer m, n
      integer m_orb, n_orb
      integer m_spin, n_spin
c
      e2 = 0.0d0
      do m = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         m_orb = (m-1)/2 + 1
         m_spin = 0
         if ((m-2*m_orb).ne.0)m_spin = 1
         do n = noso+1, nso
            n_orb = (n-1)/2 + 1
            n_spin = 0
            if ((n-2*n_orb).ne.0)n_spin = 1
            do alpha = 1, noso
               alpha_orb = (alpha-1)/2 + 1
               alpha_spin = 0
               if ((alpha-2*alpha_orb).ne.0)alpha_spin = 1
               do beta = 1, noso
                  beta_orb = (beta-1)/2 + 1
                  beta_spin = 0
                  if ((beta-2*beta_orb).ne.0)beta_spin = 1
c
                  g=0.0d0
                  if (m_spin.eq.beta_spin.and.n_spin.eq.alpha_spin)then
                     g = g_mo(beta_orb,alpha_orb,m_orb,n_orb)
                  endif
                  if (m_spin.eq.alpha_spin.and.n_spin.eq.beta_spin)then
                     g = g - g_mo(beta_orb,alpha_orb,n_orb,m_orb)
                  endif
c                 write(6,*)' g:',g
                  e2 = e2 - 0.25d0*g*(t2(m,n,alpha,beta)+
     &                                t1(m,alpha)*t1(n,beta) -
     &                                t1(n,alpha)*t1(m,beta))
c
               enddo
            enddo
         enddo
      enddo
      write(6,*)' e2 = ',e2
c
      return
      end
      subroutine t2_l_rjh(g_mo, t, t2, fock, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision fock(nmo,nmo)
      double precision t(nso,nso,nso,nso)
      double precision t2(nso,nso,nso,nso)
      double precision g, f_t, denom
c occupied orbitals
      integer i, j, k, l, m, n
      integer i_orb, j_orb, k_orb, l_orb, m_orb, n_orb
      integer i_spin, j_spin, k_spin, l_spin, m_spin, n_spin
c unoccupied orbitals
      integer a, b, c, d, e, f
      integer a_orb, b_orb, c_orb, d_orb, e_orb, f_orb
      integer a_spin, b_spin, c_spin, d_spin, e_spin, f_spin
c
c     put t into t2 expression keeping all terms linear in t2
c
      do e = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         e_orb = (e-1)/2 + 1
         e_spin = 0
         if ((e-2*e_orb).ne.0)e_spin = 1
         do f = noso+1, nso
            f_orb = (f-1)/2 + 1
            f_spin = 0
            if ((f-2*f_orb).ne.0)f_spin = 1
            do m = 1, noso
               m_orb = (m-1)/2 + 1
               m_spin = 0
               if ((m-2*m_orb).ne.0)m_spin = 1
               do n = 1, noso
                  n_orb = (n-1)/2 + 1
                  n_spin = 0
                  if ((n-2*n_orb).ne.0)n_spin = 1
c
                  g=0.0d0
                  if (e_spin.eq.m_spin.and.f_spin.eq.n_spin)then
                     g = g_mo(e_orb,f_orb,m_orb,n_orb)
                  endif
                  if (e_spin.eq.n_spin.and.f_spin.eq.m_spin)then
                     g = g - g_mo(e_orb,f_orb,n_orb,m_orb)
                  endif
c                  write(6,*)' g:',g
                  t2(e,f,m,n) = g
c
                  do a = noso+1, nso
                     a_orb = (a-1)/2 + 1
                     a_spin = 0
                     if ((a-2*a_orb).ne.0)a_spin = 1
                     f_t = 0.0d0
                     if (f_spin.eq.a_spin)then
                        f_t = fock(f_orb,a_orb)
                     endif
                     t2(e,f,m,n) = t2(e,f,m,n) + 
     &                          f_t*t(e,a,m,n)
                     f_t = 0.0d0
                     if (e_spin.eq.a_spin)then
                        f_t = fock(e_orb,a_orb)
                     endif
                     t2(e,f,m,n) = t2(e,f,m,n) - 
     &                          f_t*t(f,a,m,n)
                  enddo
c
                  do i = 1, noso
                     i_orb = (i-1)/2 + 1
                     i_spin = 0
                     if ((i-2*i_orb).ne.0)i_spin = 1
                     f_t = 0.0d0
                     if (i_spin.eq.n_spin)then
                        f_t = fock(i_orb,n_orb)
                     endif
                     t2(e,f,m,n) = t2(e,f,m,n) -
     &                          f_t*t(e,f,m,i)
                     f_t = 0.0d0
                     if (i_spin.eq.m_spin)then
                        f_t = fock(i_orb,m_orb)
                     endif
                     t2(e,f,m,n) = t2(e,f,m,n) +
     &                          f_t*t(e,f,n,i)
                  enddo
c
                  do a = noso+1, nso
                     a_orb = (a-1)/2 + 1
                     a_spin = 0
                     if ((a-2*a_orb).ne.0)a_spin = 1
                     do i = 1, noso
                        i_orb = (i-1)/2 + 1
                        i_spin = 0
                        if ((i-2*i_orb).ne.0)i_spin = 1
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.f_spin.eq.n_spin)then
                           g = g_mo(i_orb,f_orb,a_orb,n_orb)
                        endif
                        if (i_spin.eq.n_spin.and.a_spin.eq.f_spin)then
                           g = g - g_mo(i_orb,f_orb,n_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(e,f,m,n) = t2(e,f,m,n) +
     &                             g*t(e,a,m,i)
c ef flip
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.e_spin.eq.n_spin)then
                           g = g_mo(i_orb,e_orb,a_orb,n_orb)
                        endif
                        if (i_spin.eq.n_spin.and.a_spin.eq.e_spin)then
                           g = g - g_mo(i_orb,e_orb,n_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(e,f,m,n) = t2(e,f,m,n) -
     &                             g*t(f,a,m,i)
c mn flip
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.f_spin.eq.m_spin)then
                           g = g_mo(i_orb,f_orb,a_orb,m_orb)
                        endif
                        if (i_spin.eq.m_spin.and.a_spin.eq.f_spin)then
                           g = g - g_mo(i_orb,f_orb,m_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(e,f,m,n) = t2(e,f,m,n) -
     &                             g*t(e,a,n,i)
c ef+mn flip
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.e_spin.eq.m_spin)then
                           g = g_mo(i_orb,e_orb,a_orb,m_orb)
                        endif
                        if (i_spin.eq.m_spin.and.a_spin.eq.e_spin)then
                           g = g - g_mo(i_orb,e_orb,m_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(e,f,m,n) = t2(e,f,m,n) +
     &                             g*t(f,a,n,i)
                     enddo
                  enddo
c
                  do i = 1, noso
                     i_orb = (i-1)/2 + 1
                     i_spin = 0
                     if ((i-2*i_orb).ne.0)i_spin = 1
                     do j = 1, noso
                        j_orb = (j-1)/2 + 1
                        j_spin = 0
                        if ((j-2*j_orb).ne.0)j_spin = 1
                        g=0.0d0
                        if (i_spin.eq.m_spin.and.j_spin.eq.n_spin)then
                           g = g_mo(i_orb,j_orb,m_orb,n_orb)
                        endif
                        if (i_spin.eq.n_spin.and.j_spin.eq.m_spin)then
                           g = g - g_mo(i_orb,j_orb,n_orb,m_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(e,f,m,n) = t2(e,f,m,n) +
     &                             0.5d0*g*t(e,f,i,j)
                     enddo
                  enddo
                  do a = noso+1, nso
                     a_orb = (a-1)/2 + 1
                     a_spin = 0
                     if ((a-2*a_orb).ne.0)a_spin = 1
                     do b = noso+1, nso
                        b_orb = (b-1)/2 + 1
                        b_spin = 0
                        if ((b-2*b_orb).ne.0)b_spin = 1
                        g=0.0d0
                        if (e_spin.eq.a_spin.and.f_spin.eq.b_spin)then
                           g = g_mo(e_orb,f_orb,a_orb,b_orb)
                        endif
                        if (e_spin.eq.b_spin.and.f_spin.eq.a_spin)then
                           g = g - g_mo(e_orb,f_orb,b_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        t2(e,f,m,n) = t2(e,f,m,n) +
     &                             0.5d0*g*t(a,b,m,n)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c      write(6,*)' t2(e,f,m,n) L residual:'
c      call writet2(t2,nso)
      return
      end
      subroutine t2_l_q_rjh(g_mo, t2, r2, fock, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision fock(nmo,nmo)
      double precision r2(nso,nso,nso,nso)
      double precision t2(nso,nso,nso,nso)
      double precision g, f_t, denom
c occupied orbitals
      integer i, j, k, l, m, n
      integer i_orb, j_orb, k_orb, l_orb, m_orb, n_orb
      integer i_spin, j_spin, k_spin, l_spin, m_spin, n_spin
c unoccupied orbitals
      integer a, b, c, d, e, f
      integer a_orb, b_orb, c_orb, d_orb, e_orb, f_orb
      integer a_spin, b_spin, c_spin, d_spin, e_spin, f_spin
c
c     put t2 into r2 expression keeping all terms quadratic in t2
c
      do e = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         e_orb = (e-1)/2 + 1
         e_spin = 0
         if ((e-2*e_orb).ne.0)e_spin = 1
         do f = noso+1, nso
            f_orb = (f-1)/2 + 1
            f_spin = 0
            if ((f-2*f_orb).ne.0)f_spin = 1
            do m = 1, noso
               m_orb = (m-1)/2 + 1
               m_spin = 0
               if ((m-2*m_orb).ne.0)m_spin = 1
               do n = 1, noso
                  n_orb = (n-1)/2 + 1
                  n_spin = 0
                  if ((n-2*n_orb).ne.0)n_spin = 1
c
                  g=0.0d0
                  if (e_spin.eq.m_spin.and.f_spin.eq.n_spin)then
                     g = g_mo(e_orb,f_orb,m_orb,n_orb)
                  endif
                  if (e_spin.eq.n_spin.and.f_spin.eq.m_spin)then
                     g = g - g_mo(e_orb,f_orb,n_orb,m_orb)
                  endif
c                  write(6,*)' g:',g
                  r2(e,f,m,n) = g
c
                  do a = noso+1, nso
                     a_orb = (a-1)/2 + 1
                     a_spin = 0
                     if ((a-2*a_orb).ne.0)a_spin = 1
                     f_t = 0.0d0
                     if (f_spin.eq.a_spin)then
                        f_t = fock(f_orb,a_orb)
                     endif
                     r2(e,f,m,n) = r2(e,f,m,n) + 
     &                          f_t*t2(e,a,m,n)
                     f_t = 0.0d0
                     if (e_spin.eq.a_spin)then
                        f_t = fock(e_orb,a_orb)
                     endif
                     r2(e,f,m,n) = r2(e,f,m,n) - 
     &                          f_t*t2(f,a,m,n)
                  enddo
c
                  do i = 1, noso
                     i_orb = (i-1)/2 + 1
                     i_spin = 0
                     if ((i-2*i_orb).ne.0)i_spin = 1
                     f_t = 0.0d0
                     if (i_spin.eq.n_spin)then
                        f_t = fock(i_orb,n_orb)
                     endif
                     r2(e,f,m,n) = r2(e,f,m,n) -
     &                          f_t*t2(e,f,m,i)
                     f_t = 0.0d0
                     if (i_spin.eq.m_spin)then
                        f_t = fock(i_orb,m_orb)
                     endif
                     r2(e,f,m,n) = r2(e,f,m,n) +
     &                          f_t*t2(e,f,n,i)
                  enddo
c
                  do a = noso+1, nso
                     a_orb = (a-1)/2 + 1
                     a_spin = 0
                     if ((a-2*a_orb).ne.0)a_spin = 1
                     do i = 1, noso
                        i_orb = (i-1)/2 + 1
                        i_spin = 0
                        if ((i-2*i_orb).ne.0)i_spin = 1
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.f_spin.eq.n_spin)then
                           g = g_mo(i_orb,f_orb,a_orb,n_orb)
                        endif
                        if (i_spin.eq.n_spin.and.a_spin.eq.f_spin)then
                           g = g - g_mo(i_orb,f_orb,n_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        r2(e,f,m,n) = r2(e,f,m,n) +
     &                             g*t2(e,a,m,i)
c ef flip
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.e_spin.eq.n_spin)then
                           g = g_mo(i_orb,e_orb,a_orb,n_orb)
                        endif
                        if (i_spin.eq.n_spin.and.a_spin.eq.e_spin)then
                           g = g - g_mo(i_orb,e_orb,n_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        r2(e,f,m,n) = r2(e,f,m,n) -
     &                             g*t2(f,a,m,i)
c mn flip
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.f_spin.eq.m_spin)then
                           g = g_mo(i_orb,f_orb,a_orb,m_orb)
                        endif
                        if (i_spin.eq.m_spin.and.a_spin.eq.f_spin)then
                           g = g - g_mo(i_orb,f_orb,m_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        r2(e,f,m,n) = r2(e,f,m,n) -
     &                             g*t2(e,a,n,i)
c ef+mn flip
                        g=0.0d0
                        if (i_spin.eq.a_spin.and.e_spin.eq.m_spin)then
                           g = g_mo(i_orb,e_orb,a_orb,m_orb)
                        endif
                        if (i_spin.eq.m_spin.and.a_spin.eq.e_spin)then
                           g = g - g_mo(i_orb,e_orb,m_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        r2(e,f,m,n) = r2(e,f,m,n) +
     &                             g*t2(f,a,n,i)
                     enddo
                  enddo
c
                  do i = 1, noso
                     i_orb = (i-1)/2 + 1
                     i_spin = 0
                     if ((i-2*i_orb).ne.0)i_spin = 1
                     do j = 1, noso
                        j_orb = (j-1)/2 + 1
                        j_spin = 0
                        if ((j-2*j_orb).ne.0)j_spin = 1
                        g=0.0d0
                        if (i_spin.eq.m_spin.and.j_spin.eq.n_spin)then
                           g = g_mo(i_orb,j_orb,m_orb,n_orb)
                        endif
                        if (i_spin.eq.n_spin.and.j_spin.eq.m_spin)then
                           g = g - g_mo(i_orb,j_orb,n_orb,m_orb)
                        endif
c                        write(6,*)' g:',g
                        r2(e,f,m,n) = r2(e,f,m,n) +
     &                             0.5d0*g*t2(e,f,i,j)
                     enddo
                  enddo
                  do a = noso+1, nso
                     a_orb = (a-1)/2 + 1
                     a_spin = 0
                     if ((a-2*a_orb).ne.0)a_spin = 1
                     do b = noso+1, nso
                        b_orb = (b-1)/2 + 1
                        b_spin = 0
                        if ((b-2*b_orb).ne.0)b_spin = 1
                        g=0.0d0
                        if (e_spin.eq.a_spin.and.f_spin.eq.b_spin)then
                           g = g_mo(e_orb,f_orb,a_orb,b_orb)
                        endif
                        if (e_spin.eq.b_spin.and.f_spin.eq.a_spin)then
                           g = g - g_mo(e_orb,f_orb,b_orb,a_orb)
                        endif
c                        write(6,*)' g:',g
                        r2(e,f,m,n) = r2(e,f,m,n) +
     &                             0.5d0*g*t2(a,b,m,n)
                     enddo
                  enddo
c
                  do i = 1, noso
                     i_orb = (i-1)/2 + 1
                     i_spin = 0
                     if ((i-2*i_orb).ne.0)i_spin = 1
                     do j = 1, noso
                        j_orb = (j-1)/2 + 1
                        j_spin = 0
                        if ((j-2*j_orb).ne.0)j_spin = 1
                        do a = noso+1, nso
                           a_orb = (a-1)/2 + 1
                           a_spin = 0
                           if ((a-2*a_orb).ne.0)a_spin = 1
                           do b = noso+1, nso
                              b_orb = (b-1)/2 + 1
                              b_spin = 0
                              if ((b-2*b_orb).ne.0)b_spin = 1
                              g=0.0d0
                              if (i_spin.eq.a_spin.and.
     &                            j_spin.eq.b_spin)then
                                 g = g_mo(i_orb,j_orb,a_orb,b_orb)
                              endif
                              if (i_spin.eq.b_spin.and.
     &                            j_spin.eq.a_spin)then
                                 g = g - g_mo(i_orb,j_orb,b_orb,a_orb)
                              endif
c                              write(6,*)' g:',g
                              r2(e,f,m,n) = r2(e,f,m,n) +
     &                             0.5d0*g* ( t2(f,b,j,n)*t2(e,a,i,m) -
     &                                        t2(e,b,j,n)*t2(f,a,i,m) -
     &                                        t2(f,b,j,m)*t2(e,a,i,n) +
     &                                        t2(e,b,j,m)*t2(f,a,i,n) -
     &                                        t2(f,b,m,n)*t2(e,a,i,j) +
     &                                        t2(e,b,m,n)*t2(f,a,i,j) -
     &                                        t2(e,f,m,i)*t2(a,b,n,j) +
     &                                        t2(e,f,n,i)*t2(a,b,m,j) +
     &                                  0.5d0*t2(e,f,i,j)*t2(a,b,m,n) )
 
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c      write(6,*)' t2(e,f,m,n) L+Q residual:'
c      call writet2(r2,nso)
      return
      end
      subroutine t2_update_rjh(g_mo, t2, r2, fock, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision fock(nmo,nmo)
      double precision t2(nso,nso,nso,nso)
      double precision r2(nso,nso,nso,nso)
      double precision g, denom, r2norm
c occupied orbitals
      integer i, j, k, l, m, n
      integer i_orb, j_orb, k_orb, l_orb, m_orb, n_orb
      integer i_spin, j_spin, k_spin, l_spin, m_spin, n_spin
c unoccupied orbitals
      integer a, b, c, d, e, f
      integer a_orb, b_orb, c_orb, d_orb, e_orb, f_orb
      integer a_spin, b_spin, c_spin, d_spin, e_spin, f_spin
      double precision ddot
      external ddot
c
c     update t with delta t (in t2)
c
      do e = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         e_orb = (e-1)/2 + 1
         e_spin = 0
         if ((e-2*e_orb).ne.0)e_spin = 1
         do f = noso+1, nso
            f_orb = (f-1)/2 + 1
            f_spin = 0
            if ((f-2*f_orb).ne.0)f_spin = 1
            do m = 1, noso
               m_orb = (m-1)/2 + 1
               m_spin = 0
               if ((m-2*m_orb).ne.0)m_spin = 1
               do n = 1, noso
                  n_orb = (n-1)/2 + 1
                  n_spin = 0
                  if ((n-2*n_orb).ne.0)n_spin = 1
c
                  denom = fock(e_orb,e_orb) + fock(f_orb,f_orb) -
     &                    fock(m_orb,m_orb) - fock(n_orb,n_orb)
c                 write(6,*)' denom: ', denom
c
                  t2(e,f,m,n) = t2(e,f,m,n) - r2(e,f,m,n)/denom
               enddo
            enddo
         enddo
      enddo
      r2norm =  dsqrt(ddot(nso**4,r2,1,r2,1))
      write(6,*)' r2norm = ', r2norm
      return
      end
      subroutine t1_rjh(g_mo, r1, t2, fock, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision fock(nmo,nmo)
      double precision r1(nso,nso)
      double precision t2(nso,nso,nso,nso)
      double precision g, f_t, denom
c occupied orbitals
      integer i, j, k, l, m, n
      integer i_orb, j_orb, k_orb, l_orb, m_orb, n_orb
      integer i_spin, j_spin, k_spin, l_spin, m_spin, n_spin
c unoccupied orbitals
      integer a, b, c, d, e, f
      integer a_orb, b_orb, c_orb, d_orb, e_orb, f_orb
      integer a_spin, b_spin, c_spin, d_spin, e_spin, f_spin
c
c     put t into t2 expression keeping all terms linear in t2
c
c      write(6,*)' debug inside t1_rjh '
c      call jan_debug_print('MOINTS', g_mo, nmo, nmo,  
c     $                      nmo, nmo)
c      write(6,*)' T2 with quadratic terms. '
c      call writet2(t2,nso)
c      write(6,*) ' Fock Matrix'
c      call output(fock,1,nmo,1,nmo,nmo,nmo,1)
c
      do e = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         e_orb = (e-1)/2 + 1
         e_spin = 0
         if ((e-2*e_orb).ne.0)e_spin = 1
         do m = 1, noso
            m_orb = (m-1)/2 + 1
            m_spin = 0
            if ((m-2*m_orb).ne.0)m_spin = 1
c
            r1(e,m) = 0.0d0
            if (e_spin.eq.m_spin)then
               r1(e,m) = fock(e_orb,m_orb)
            endif
c
            do j = 1, noso
               j_orb = (j-1)/2 + 1
               j_spin = 0
               if ((j-2*j_orb).ne.0)j_spin = 1
               do b = noso+1, nso
                  b_orb = (b-1)/2 + 1
                  b_spin = 0
                  if ((b-2*b_orb).ne.0)b_spin = 1
c
                  f_t = 0.0d0
                  if (j_spin.eq.b_spin)then
                     f_t = fock(j_orb,b_orb)
                  endif
c
                  r1(e,m) = r1(e,m) + f_t*t2(e,b,m,j)
c
                  do i = 1, noso
                     i_orb = (i-1)/2 + 1
                     i_spin = 0
                     if ((i-2*i_orb).ne.0)i_spin = 1
c
                     g=0.0d0
                     if (i_spin.eq.m_spin.and.j_spin.eq.b_spin)then
                        g = g_mo(i_orb,j_orb,m_orb,b_orb)
                     endif
                     if (i_spin.eq.b_spin.and.j_spin.eq.m_spin)then
                        g = g - g_mo(i_orb,j_orb,b_orb,m_orb)
                     endif

c                     write(6,*)' e,m,j,b,i,g:',e,m,j,b,i,g
c
                     r1(e,m) = r1(e,m) - 0.5d0*g*t2(e,b,i,j)
c
                  enddo
               enddo
            enddo
            do b = noso+1, nso
               b_orb = (b-1)/2 + 1
               b_spin = 0
               if ((b-2*b_orb).ne.0)b_spin = 1
c
               do i = 1, noso
                  i_orb = (i-1)/2 + 1
                  i_spin = 0
                  if ((i-2*i_orb).ne.0)i_spin = 1
c
                  do a = noso+1, nso
                     a_orb = (a-1)/2 + 1
                     a_spin = 0
                     if ((a-2*a_orb).ne.0)a_spin = 1
                     g=0.0d0
                     if (i_spin.eq.a_spin.and.e_spin.eq.b_spin)then
                        g = g_mo(i_orb,e_orb,a_orb,b_orb)
                     endif
                     if (i_spin.eq.b_spin.and.e_spin.eq.a_spin)then
                        g = g - g_mo(i_orb,e_orb,b_orb,a_orb)
                     endif
c
c                     write(6,*)' g:',g
c
                     r1(e,m) = r1(e,m) + 0.5d0*g*t2(a,b,i,m)
c
                  enddo
               enddo
            enddo
         enddo
      enddo
c      write(6,*)' t1(e,m) residual:'
c      call writet1(r1,nso)
      return
      end
      subroutine t1_update_rjh(g_mo, t1, r1, fock, nmo, nso, noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision fock(nmo,nmo)
      double precision r1(nso,nso)
      double precision t1(nso,nso)
      double precision g, denom, r1norm
c occupied orbitals
      integer i, j, k, l, m, n
      integer i_orb, j_orb, k_orb, l_orb, m_orb, n_orb
      integer i_spin, j_spin, k_spin, l_spin, m_spin, n_spin
c unoccupied orbitals
      integer a, b, c, d, e, f
      integer a_orb, b_orb, c_orb, d_orb, e_orb, f_orb
      integer a_spin, b_spin, c_spin, d_spin, e_spin, f_spin
      double precision ddot
      external ddot
c
c     update t with delta t (in t2)
c
      do e = noso+1, nso
c get spatial orbital and whether alpha or beta spin
         e_orb = (e-1)/2 + 1
         e_spin = 0
         if ((e-2*e_orb).ne.0)e_spin = 1
         do m = 1, noso
            m_orb = (m-1)/2 + 1
            m_spin = 0
            if ((m-2*m_orb).ne.0)m_spin = 1
c
            denom = fock(e_orb,e_orb) - fock(m_orb,m_orb)
c            write(6,*)' denom: ', denom
c
            t1(e,m) = t1(e,m) - r1(e,m)/denom
         enddo
      enddo
      r1norm =  dsqrt(ddot(nso**2,r1,1,r1,1))
      write(6,*)' r1norm = ', r1norm
      return
      end
      subroutine jan_h(
     $     rtdb, basis,
     $     n1, n2,
     $     lda1, lda2,
     $     c1t, c2t,
     $     h)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb, basis
      integer lda1, lda2, n1, n2
      double precision c1t(lda1, n1), c2t(lda2, n2)
      double precision h(n1, n2)
c
c     Return hij = sum(kl) c1t(i,k) hAO(k,l) c2t(j,l)
c
c     Note that the transposed MO coeffs are passed in 
c     to be compatible with jan_full_transform
c
      integer geom, nbf, g_tmp
      integer l_tmp1, k_tmp1, l_tmp2, k_tmp2
c
      call int_init(rtdb, 1, basis)
      if (.not. bas_geom(basis, geom))
     $     call errquit('jan_transform: basis ', basis, BASIS_ERR)
      call schwarz_init(geom, basis)
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('jan_h: nbf',basis, BASIS_ERR)
c
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'tmp1', l_tmp1, k_tmp1))
     $     call errquit('tmp1', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'tmp2', l_tmp2, k_tmp2))
     $     call errquit('tmp2', nbf*nbf, MA_ERR)
c
      if (.not. ga_create(mt_dbl, nbf, nbf, 'tmp', 1, 1, g_tmp))
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
      call ga_zero(g_tmp)
      call int_1e_ga(basis, basis, g_tmp, 'kinetic', .false.)
      call int_1e_ga(basis, basis, g_tmp, 'potential', .false.)
      call ga_get(g_tmp, 1, nbf, 1, nbf, dbl_mb(k_tmp1), nbf)
      if (.not. ga_destroy(g_tmp)) call errquit('ga?',0, GA_ERR)
c
      call dgemm('n', 'n', n1, nbf, nbf, 
     $     1.0d0, c1t, lda1, dbl_mb(k_tmp1), nbf, 
     $     0.0d0, dbl_mb(k_tmp2), n1)
      call dgemm('n', 't', n1, n2, nbf,
     $     1.0d0, dbl_mb(k_tmp2), n1, c2t, lda2,
     $     0.0d0, h, n1)
c
c      write(6,*) ' Transformed H'
c      call output(h, 1, n1, 1, n2, n1, n2, 1)
c
      call schwarz_tidy()
      call int_terminate
c
      if (.not. ma_chop_stack(l_tmp1)) call errquit('ma',0, MA_ERR)
c
      return
      end
      subroutine uccsdtest_lambda(nbf, nmo, no, nv, 
     $        t1, ct, pt, ht, work)
      implicit none
c
      integer nbf, nmo, no, nv
      double precision t1(nv, no)
      double precision ct(nmo, nbf)
      double precision pt(nmo, nbf)
      double precision ht(nmo, nbf)
      double precision work(nmo, nmo)
c
      integer i, a
c
c     Particle transformation
c
      call dfill(nmo*nmo, 0.0d0, work, 1)
      call dfill(nmo, 1.0d0, work, nmo+1)
      do i = 1, no
         do a = 1, nv
            work(i,a+no) = work(i,a+no) - t1(a,i)
         end do
      end do
      call dgemm('t','n',nmo,nbf,nmo,
     $     1.0d0, work, nmo, ct, nmo, 
     $     0.0d0, pt, nmo)
c
c     Hole transformation
c
      call dfill(nmo*nmo, 0.0d0, work, 1)
      call dfill(nmo, 1.0d0, work, nmo+1)
      do i = 1, no
         do a = 1, nv
            work(a+no,i) = work(a+no,i) + t1(a,i)
         end do
      end do
      call dgemm('t','n',nmo,nbf,nmo,
     $     1.0d0, work, nmo, ct, nmo, 
     $     0.0d0, ht, nmo)
c
      end
      subroutine get_new_f_g(rtdb, basis, fock, g_mo, t1, nso, nmo)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
      integer rtdb
c
      integer nmo, g_tmp, l_occ, k_occ,
     $     l_eval, k_eval, l_mos, k_mos, l_most, k_most
      integer nocc, nvirt, nopen, nclosed, nso, noso
      integer l_t, k_t, l_t2, k_t2, l_f, k_f
      integer l_t1, k_t1, l_r1, k_r1
      integer l_pmost, k_pmost, l_hmost, k_hmost
      integer l_work, k_work, l_h, k_h
      integer l_t1_vomos, k_t1_vomos
c
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision fock(nmo,nmo)
      double precision t1(nso,nso)
      integer basis, nbf
      character*255 movecs      ! Name of movector file
      character*80 title, name_of_basis, scftype
      integer nbf_file, nsets, nmo_file(2)
      logical movecs_read, movecs_read_header
      external movecs_read, movecs_read_header
c
      logical int_normalize
      external int_normalize
c
c     Read the MO vectors and evals from a RHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('jantest: failed to read movecs header',911,
     &       DISK_ERR)
c      write(6,*) ' Read movecs header from ', movecs
c      write(6,*) ' Job title :                ', 
c     $     title(1:inp_strlen(title))
c      write(6,*) ' Basis name:                ', 
c     $     name_of_basis(1:inp_strlen(name_of_basis))
      nmo = nmo_file(1)
      if (.not. rtdb_get(rtdb, 'scf:nclosed', mt_int, 1, nocc))
     $     call errquit('nocc?',0, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nopen', mt_int, 1, nopen))
     $     call errquit('nopen?',0, RTDB_ERR)
      if (nopen .ne. 0) call errquit('asjdlfkadjsl',0, UNKNOWN_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0, BASIS_ERR)
      nvirt= nmo - nocc
c      write(6,*) ' No. of closed shells       ', nocc
c      write(6,*) ' No. of molecular orbitals: ', nmo
c      write(6,*) ' No. of basis functions:    ', nbf
c
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 0, 0, g_tmp))
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'eval',l_eval, k_eval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'mos', l_mos, k_mos))
     $     call errquit('ma mos', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'mos', l_most, k_most))
     $     call errquit('ma mos', nbf*nbf, MA_ERR)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_eval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_mos), nbf)
      call util_transpose(dbl_mb(k_mos),nbf,dbl_mb(k_most),nmo,
     $     nbf,nmo)
c
c      write(6,*) ' Orbital eigenvalues '
c      call output(dbl_mb(k_eval),1,nmo,1,1,nmo,1,1)
c      write(6,*) ' MOs'
c      call output(dbl_mb(k_mos),1,nbf,1,nmo,nbf,nmo,1)
c      write(6,*) ' MOs T'
c      call output(dbl_mb(k_most),1,nmo,1,nbf,nmo,nbf,1)
c
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0, GA_ERR)
c
c     Transform MO integrals in Dirac order to new T1 bases
c
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'pmost', l_pmost, k_pmost))
     $     call errquit('ma pmost', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'hmost', l_hmost, k_hmost))
     $     call errquit('ma hmost', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'work', l_work, k_work))
     $     call errquit('ma work', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'1e-', l_h, k_h))
     $     call errquit('ma h', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nvirt*nocc,'t1_vomos', 
     $     l_t1_vomos, k_t1_vomos))
     $     call errquit('ma t1_vomos', nvirt*nocc, MA_ERR)
      call pack_t1(t1, dbl_mb(k_t1_vomos), nso, nmo, nocc, nvirt) 
      call uccsdtest_lambda(nbf, nmo, nocc, nvirt, 
     $                      dbl_mb(k_t1_vomos), dbl_mb(k_most), 
     $                      dbl_mb(k_pmost), 
     $                      dbl_mb(k_hmost), dbl_mb(k_work))
c
c      write(6,*) ' P mos'
c      call output(dbl_mb(k_pmost), 1, nmo, 1, nmo, nmo, nmo, 1)
c
c      write(6,*) ' H mos'
c      call output(dbl_mb(k_hmost), 1, nmo, 1, nmo, nmo, nmo, 1)
c
      call jan_h(rtdb, basis, nmo, nmo, nmo, nmo,
     $           dbl_mb(k_pmost), dbl_mb(k_hmost), dbl_mb(k_h))
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_pmost), dbl_mb(k_pmost), dbl_mb(k_hmost),
     $     dbl_mb(k_hmost), g_mo, 'Dirac')
c      call jan_debug_print('MOINTS', g_mo, nmo, nmo,  
c     $                      nmo, nmo)
      call build_f(fock, g_mo, dbl_mb(k_h), nmo, nocc)
c
c     Tidy up
c
      if (.not. ma_chop_stack(l_occ)) call errquit(' ma chop?', 0,
     &       MA_ERR)
c
      return
c
      end
      subroutine pack_t1(t1, t1_vomos, nso, nmo, nocc, nvirt)
      implicit none
c
      integer nso, nmo, nocc, nvirt
      integer i, j, i_orb, j_orb, i_so, j_so
      double precision t1_vomos(nvirt, nocc)
      double precision t1(nso,nso)
c
c      write(6,*)' t1 from pack_t1 '
c      do i = 1, nso
c         do j = 1, nso
c            if (abs(t1(i,j)).gt.1d-6)
c     &         write(6,*)i,j,t1(i,j)
c         enddo
c      enddo
c
      call dfill(nvirt*nocc, 0.0d0, t1_vomos, 1)
c
      do i = 1, nvirt
         do j = 1, nocc
            i_orb = nocc + i
            j_orb = j
            i_so = 2*(i_orb-1)+1
            j_so = 2*(j_orb-1)+1
            t1_vomos(i,j) = t1(i_so,j_so)
         enddo
      enddo
c      write(6,*)' t1_vomos '
c      do i = 1, nvirt
c         do j = 1, nocc
c            if (abs(t1_vomos(i,j)).gt.1d-6)
c     &         write(6,*)i,j,t1_vomos(i,j)
c         enddo
c      enddo
      return
      end
      subroutine build_f(fock, g_mo, h_mo, nmo, nocc)
      implicit none
c
      integer nmo, nocc
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision fock(nmo,nmo), h_mo(nmo,nmo)
c
c     fij = hij + 2*<ik|jk> - <ik|kj>
c
      integer i, j, k
c
      call dfill(nmo*nmo, 0.0d0, fock, 1)
c
      do j = 1, nmo
         do i = 1, nmo
            fock(i,j) = h_mo(i,j)
            do k = 1, nocc
               fock(i,j) = fock(i,j) + 2.0d0*g_mo(i,k,j,k)
     &                               - g_mo(i,k,k,j)
            end do
         end do
      end do
c
c      write(6,*) ' FOCK'
c      call output(fock, 1, nmo, 1, nmo, nmo, nmo, 1)
c
      return
      end

c $Id$
