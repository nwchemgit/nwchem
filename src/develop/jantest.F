      logical function jantest(rtdb)
*     
*     $$
*
      implicit none
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
      integer l_t, k_t, l_t2, k_t2
      logical int_normalize
      external int_normalize
c     
c     load the geometry/basis set and get info
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $        call errquit('scf_init: no ao basis set', 0)
      if (.not.int_normalize(basis))
     $           call errquit('scf:int_normalize failed', 0)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0)
c
c     Read the MO vectors and evals from a RHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('jantest: failed to read movecs header',911)
      write(6,*) ' Read movecs header from ', movecs
      write(6,*) ' Job title :                ', 
     $     title(1:inp_strlen(title))
      write(6,*) ' Basis name:                ', 
     $     name_of_basis(1:inp_strlen(name_of_basis))
      nmo = nmo_file(1)
      if (.not. rtdb_get(rtdb, 'scf:nclosed', mt_int, 1, nocc))
     $     call errquit('nocc?',0)
      if (.not. rtdb_get(rtdb, 'scf:nopen', mt_int, 1, nopen))
     $     call errquit('nopen?',0)
      if (nopen .ne. 0) call errquit('asjdlfkadjsl',0)
      nvirt= nmo - nocc
      write(6,*) ' No. of closed shells       ', nclosed
      write(6,*) ' No. of molecular orbitals: ', nmo
      write(6,*) ' No. of basis functions:    ', nbf
c
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 0, 0, g_tmp))
     &     call errquit('scf_v_g: tmp', 0)
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf)
      if (.not. ma_push_get(mt_dbl, nbf,'eval',l_eval, k_eval))
     $     call errquit('ma eval', nbf)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'mos', l_mos, k_mos))
     $     call errquit('ma mos', nbf*nbf)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'mos', l_most, k_most))
     $     call errquit('ma mos', nbf*nbf)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_eval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0)
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
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0)
c
c     Make all AO integrals
c
      if (.not. ma_push_get(mt_dbl,nbf**4,'aoint',l_aoint,k_aoint))
     $     call errquit('allocation of AO integrals failed',nbf**4)
      call jan_all_ao_integrals(rtdb,basis,nbf,'dirac',dbl_mb(k_aoint))
      call jan_debug_print('AOINTS',dbl_mb(k_aoint), nbf,  nbf,  nbf,
     $     nbf)
c
c     Make all MO integrals in Dirac order
c
      if (.not. ma_push_get(mt_dbl,nmo**4,'aoint',l_moint,k_moint))
     $     call errquit('allocation of MO integrals failed',nbf**4)
      call jan_full_transform(
     $     rtdb, basis, 
     $     nmo, nmo, nmo, nmo,
     $     nmo, nmo, nmo, nmo,
     $     dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most),dbl_mb(k_most), 
     $     dbl_mb(k_moint), 'Dirac')
      call jan_debug_print('MOINTS',dbl_mb(k_moint), nmo,  nmo,  nmo,
     $     nmo)
c
c     do some incore cc
c     set nso=2*nmo, noso=2*nocc
c
      nso=2*nmo
      noso=2*nocc
      if (.not. ma_push_get(mt_dbl,nso**4,'t amps',l_t,k_t))
     $     call errquit('allocation of t amplitudes failed',nso**4)
      if (.not. ma_push_get(mt_dbl,nso**4,'t2 amps',l_t2,k_t2))
     $     call errquit('allocation of t2 amplitudes failed',nso**4)
      call ccsd_incore(dbl_mb(k_aoint), dbl_mb(k_moint), 
     &                 dbl_mb(k_eval), dbl_mb(k_t), dbl_mb(k_t2),
     &                 nbf, nmo, nocc, nso, noso)
c
c     Tidy up
c
      if (.not. ma_chop_stack(l_occ)) call errquit(' ma chop?', 0)
      if (.not. bas_destroy(basis)) call errquit(' bas ?',0)
      if (.not. geom_destroy(geom)) call errquit(' geom ?',0)
      jantest = .true.
c
      end
      subroutine jan_all_ao_integrals(rtdb, basis, nbf, order, ao)
      implicit none
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
     $     'ao_fock_2e: problem with call to bas_numcont', basis)
      call int_mem_2e4c(maxg2,maxs2)
      if (.not. ma_push_get(mt_dbl,maxs2,'scr',l_scr, k_scr))
     $   call errquit('ma scr',maxg2)
      if (.not. ma_push_get(mt_dbl,maxg2,'buf',l_buf, k_buf))
     $   call errquit('ma buf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'i',l_i, k_i))
     $   call errquit('ma ibuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'j',l_j, k_j))
     $   call errquit('ma jbuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'k',l_k, k_k))
     $   call errquit('ma kbuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'l',l_l, k_l))
     $   call errquit('ma lbuf',maxg2)
c
      call jan_do_all_ao_integrals(basis, dbl_mb(k_buf), dbl_mb(k_scr),
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     maxg2, maxs2, nbf, nsh, order, ao)
c
      if (.not. ma_chop_stack(l_scr)) call errquit('janallao: ma?',0)
c
      call int_terminate
c
      end
      subroutine jan_do_all_ao_integrals(
     $     basis, buf, scr, ilab, jlab, klab,
     $  llab, maxg2, maxs2, nbf, nsh, order, ao)
      implicit none
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
      if (order .eq. 'mulliken') then
         omulliken = .true.
      else if (order .eq. 'dirac') then
         omulliken = .false.
      else
         call errquit(' unknown order',0)
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
         call errquit('jan_full_trans: unkown integral option',0)
      endif
c
c     Initialize integrals and Schwarz screening
c
      if (.not. bas_geom(basis, geom))
     $     call errquit('jan_transform: basis ', basis)
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
c      
      if (.not. bas_numcont(basis, nsh)) call errquit(
     $     'jan_transform: bas_numcont', basis)
      if (.not. bas_nbf_cn_max(basis,maxbfsh)) call errquit(
     $     'jan_transform: bas_nbf_cn_max', basis)
c
      lenhalf = n3*n4*maxbfsh**2
      lenthird= n2*n3*n4*maxbfsh
c
      if (.not. ma_push_get(mt_dbl,lenhalf,'half',l_half, k_half))
     $     call errquit('ma half', lenhalf)
      if (.not. ma_push_get(mt_dbl,lenthird,'third',l_third, k_third))
     $     call errquit('ma third', lenthird)
c
      call dfill(n1*n2*n3*n4, 0.0d0, full, 1)
      do ksh = 1, nsh
         if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $        call errquit('jan_transform: bas_cn2bfr',basis)
         call dfill(n2*n3*n4*(khi-klo+1), 0.0d0, dbl_mb(k_third), 1)
         do lsh = 1, nsh
            if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $           call errquit('jan_transform: bas_cn2bfr',basis)
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
      if (.not. ma_pop_stack(l_third)) call errquit('ma third',0)
      if (.not. ma_pop_stack(l_half)) call errquit('ma half',0)
c
      call schwarz_tidy()
      call int_terminate
c
      end
      subroutine jan_asym_trans(full,n1,n2,n3,n4,side)
      implicit none
c
      integer n1, n2, n3, n4
      double precision full(n1,n2,n3,n4)
      character*(*) side
c
      integer p, q, r, s
      double precision tmp
c
      if (side .eq. 'left') then
         if (n1 .ne. n2) call errquit('jan_asym_trans: left', n1)
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
         if (n3 .ne. n4) call errquit('jan_asym_trans: right', n3)
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
     $     'jan_transform: bas_numcont', basis)
      if (.not. bas_nbf_cn_max(basis,maxbfsh)) call errquit(
     $     'jan_transform: bas_nbf_cn_max', basis)
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('jan_transform: nbf',basis)
      if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $     call errquit('jan_transform: bas_cn2bfr',basis)
      if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $     call errquit('jan_transform: bas_cn2bfr',basis)
      call int_mem_2e4c(maxg2,maxs2)
      lenaobuf = (khi-klo+1)*(lhi-llo+1)*maxbfsh*n2 ! (iq|kl)
c
c     Allocate scratch space for integrals and buffers for
c     transformation
c
      if (.not. ma_push_get(mt_dbl,maxs2,'scr',l_scr, k_scr))
     $   call errquit('ma scr',maxs2)
      if (.not. ma_push_get(mt_dbl,maxg2,'buf',l_buf, k_buf))
     $   call errquit('ma buf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'i',l_i, k_i))
     $   call errquit('ma ibuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'j',l_j, k_j))
     $   call errquit('ma jbuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'k',l_k, k_k))
     $   call errquit('ma kbuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'l',l_l, k_l))
     $   call errquit('ma lbuf',maxg2)
      if (.not. ma_push_get(mt_dbl,lenaobuf,'aobuf',l_aobuf, k_aobuf))
     $   call errquit('ma aobuf',lenaobuf)
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
     $     ('jan_transform: chopping stack', 0)
c
      end
      subroutine jan_do_half_transform(
     $     basis, 
     $     buf, scr, ilab, jlab, klab, llab, maxg2, maxs2, 
     $     nbf, nsh, maxbfsh, n1, n2, tol2e, 
     $     ksh, lsh, klo, khi, llo, lhi,
     $     half, c1t, c2t, ld1, ld2, aobuf, ochargecloud)
      implicit none
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
     $        call errquit('jan_do_half_transform',ish)
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
 1                   format(4i5,2x,f12.6)
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine ccsd_incore(g_ao, g_mo, e, t, t2, 
     &                       nbf, nmo, nocc, nso, noso)
      implicit none
      integer nbf, nmo, nocc, nso, noso
c
c     nbf = number of basis functions
c     nmo = number of molecular orbitals
c     nocc = number of occupied orbitals
c     nso = number of spin orbitals (for now set to 2*nmo)
c     noso = number of occupied spin orbitals (for now set to 2*nocc)
c
c     toy coupled cluster program based on the green book
c
      double precision g_ao(nbf,nbf,nbf,nbf)
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision e(nmo)
      double precision t(nso,nso,nso,nso)
      double precision t2(nso,nso,nso,nso)
c
c     cc stuff (greek indices denote occupied spin-orbitals)
c     cc stuff (r, s, t, u, ... denote unoccupied spin-orbitals)
c
      write(6,*)' nbf, nmo, nocc, nso, noso: ',
     &            nbf, nmo, nocc, nso, noso
c
c     First generate an initial guess for t2
c
      call dfill(nso**4, 0.0d0, t, 1) 
      call dfill(nso**4, 0.0d0, t2, 1) 

      call t2_init(g_mo, e, t, nmo, nso, noso)
c
c     next put this guess into t2 expression keeping all terms linear in t2
c
      call t2_l(g_mo, e, t, t2, nmo, nso, noso)
c
c     next put this into t2 expression keeping all terms
c
      call dcopy(nso**4, t2, 1, t, 1)
      call dfill(nso**4, 0.0d0, t2, 1) 
      call t2_l_q(g_mo, e, t, t2, nmo, nso, noso)
c
c how did we lose symmetry in t2?????
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
         write(6,*)' m, m_orb, m_spin: ',
     &               m, m_orb, m_spin
         do n = noso+1, nso
            n_orb = (n-1)/2 + 1
            n_spin = 0
            if ((n-2*n_orb).ne.0)n_spin = 1
            write(6,*)' n, n_orb, n_spin: ',
     &                  n, n_orb, n_spin
            do alpha = 1, noso
               alpha_orb = (alpha-1)/2 + 1
               alpha_spin = 0
               if ((alpha-2*alpha_orb).ne.0)alpha_spin = 1
               write(6,*)' alpha, alpha_orb, alpha_spin: ',
     &                     alpha, alpha_orb, alpha_spin
               do beta = 1, noso
                  beta_orb = (beta-1)/2 + 1
                  beta_spin = 0
                  if ((beta-2*beta_orb).ne.0)beta_spin = 1
                  write(6,*)' beta, beta_orb, beta_spin: ',
     &                        beta, beta_orb, beta_spin
c
                  denom = e(m_orb)+e(n_orb)-e(alpha_orb)-e(beta_orb)
                  write(6,*)' denom: ', denom
                  g=0.0d0
                  if (m_spin.eq.alpha_spin.and.n_spin.eq.beta_spin)then
                     g = g_mo(m_orb,n_orb,alpha_orb,beta_orb)
                  endif
                  if (m_spin.eq.beta_spin.and.n_spin.eq.alpha_spin)then
                     g = g - g_mo(m_orb,n_orb,beta_orb,alpha_orb)
                  endif
                  write(6,*)' g:',g
                  t(m,n,alpha,beta) = g/denom
                  write(6,*)' m,n,alpha,beta,t(m,n,alpha,beta):',
     &                        m,n,alpha,beta,t(m,n,alpha,beta)
               enddo
            enddo
         enddo
      enddo
      write(6,*)' Initial Guess for T2 '
      call writet(t,nso)
      call correlation(g_mo,t,nmo,nso,noso)
      return
      end
      subroutine t2_l(g_mo, e, t, t2, nmo, nso, noso)
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
                  write(6,*)' denom: ', denom
                  g=0.0d0
                  if (m_spin.eq.alpha_spin.and.n_spin.eq.beta_spin)then
                     g = g_mo(m_orb,n_orb,alpha_orb,beta_orb)
                  endif
                  if (m_spin.eq.beta_spin.and.n_spin.eq.alpha_spin)then
                     g = g - g_mo(m_orb,n_orb,beta_orb,alpha_orb)
                  endif
                  write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
                        t2(m,n,alpha,beta) = t2(m,n,alpha,beta) + 
     &                                       g*t(n,p,beta,gamma)
                     enddo
                  enddo
                  t2(m,n,alpha,beta) = t2(m,n,alpha,beta)/denom
                  write(6,*)' t2(m,n,alpha,beta):',
     &                        t2(m,n,alpha,beta)
               enddo
            enddo
         enddo
      enddo
      write(6,*)' T2 with linear terms. '
      call writet(t2,nso)
      call correlation(g_mo,t2,nmo,nso,noso)
      return
      end
      subroutine t2_l_q(g_mo, e, t, t2, nmo, nso, noso)
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
                  write(6,*)' denom: ', denom
                  g=0.0d0
                  if (m_spin.eq.alpha_spin.and.n_spin.eq.beta_spin)then
                     g = g_mo(m_orb,n_orb,alpha_orb,beta_orb)
                  endif
                  if (m_spin.eq.beta_spin.and.n_spin.eq.alpha_spin)then
                     g = g - g_mo(m_orb,n_orb,beta_orb,alpha_orb)
                  endif
                  write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                        write(6,*)' g:',g
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
                              write(6,*)' g:',g
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
                  write(6,*)' t2(m,n,alpha,beta):',
     &                        t2(m,n,alpha,beta)
               enddo
            enddo
         enddo
      enddo
      write(6,*)' T2 with all terms. '
      call writet(t2,nso)
      call correlation(g_mo,t2,nmo,nso,noso)
      return
      end
      subroutine writet(t,len)
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
      subroutine correlation(g_mo,t,nmo,nso,noso)
      implicit none
      integer nmo, nso, noso
      double precision g_mo(nmo,nmo,nmo,nmo)
      double precision t(nso,nso,nso,nso)
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
         do n = noso+1, m-1
            n_orb = (n-1)/2 + 1
            n_spin = 0
            if ((n-2*n_orb).ne.0)n_spin = 1
            do alpha = 1, noso
               alpha_orb = (alpha-1)/2 + 1
               alpha_spin = 0
               if ((alpha-2*alpha_orb).ne.0)alpha_spin = 1
               do beta = 1, alpha-1
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
                  write(6,*)' g:',g
                  e2 = e2 + g*t(m,n,alpha,beta)
c
               enddo
            enddo
         enddo
      enddo
      write(6,*)' e2 = ',e2
c
      return
      end
