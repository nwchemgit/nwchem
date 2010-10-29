      logical function rjhtest(rtdb)
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
      integer rtdb
c
      integer basis, geom, nbf
      integer l_aoint, k_aoint        ! nbf**4 array of AO integrals
      character*255 movecs      ! Name of movector file
      character*80 title, name_of_basis, scftype
      integer nbf_file, nsets, nmo_file(2)
      logical movecs_read, movecs_read_header
      external movecs_read, movecs_read_header
c
      integer nmo, nalpha, nbeta, g_tmp, l_occ, k_occ,
     $     l_aeval, k_aeval, l_beval, k_beval, l_amos, k_amos,
     $     l_bmos, k_bmos
      integer nocc, nvirt, nso, l_sos, k_sos, l_seval, k_seval,
     $     l_sospin, k_sospin, l_tmp, k_tmp, l_soint, k_soint
      integer l_t2, k_t2, l_t1, k_t1
      integer l_r2, k_r2, l_r1, k_r1
      integer i, iter
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
      if (.not.int_normalize(rtdb, basis))
     $           call errquit('scf:int_normalize failed', 0, INT_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0, BASIS_ERR)
c
c     Read the MO vectors and evals from a UHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('rjhtest: failed to read movecs header',0, DISK_ERR)
      if (nsets.ne.2) call errquit('rjhtest:UHF ONLY!',0, INPUT_ERR)
      write(6,*) ' Read movecs header from ', movecs
      write(6,*) ' Job title :                ', 
     $     title(1:inp_strlen(title))
      write(6,*) ' Basis name:                ', 
     $     name_of_basis(1:inp_strlen(name_of_basis))
      nmo = nmo_file(1)
      if (nmo_file(1).ne.nmo_file(2)) call errquit('rjhtest:Uh?',0,
     &       UNKNOWN_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nalpha', mt_int, 1, nalpha))
     $     call errquit('rjhtest:nalpha?',0, RTDB_ERR)
      if (.not. rtdb_get(rtdb, 'scf:nbeta', mt_int, 1, nbeta))
     $     call errquit('nbeta?',0, RTDB_ERR)
      nocc = nalpha+nbeta
      nso  = nmo*2
      nvirt= nso - nocc
      write(6,*) ' No. of alpha electrons:    ', nalpha
      write(6,*) ' No. of beta  electrons:    ', nbeta
      write(6,*) ' No. of molecular orbitals: ', nmo
      write(6,*) ' No. of basis functions:    ', nbf
      write(6,*) ' No. of spin-orbitals (SO)  ', nso
      write(6,*) ' No. of occupied SOs:       ', nocc
      write(6,*) ' No. of virtual SOs:        ', nvirt
c
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 1, 1, g_tmp))
     &     call errquit('scf_v_g: tmp', 0, GA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'aeval',l_aeval, k_aeval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf,'beval',l_beval, k_beval))
     $     call errquit('ma eval', nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'amos', l_amos, k_amos))
     $     call errquit('ma amos', nbf*nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'bmos', l_bmos, k_bmos))
     $     call errquit('ma bmos', nbf*nbf, MA_ERR)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_aeval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_amos), nbf)
      if (.not. movecs_read(movecs, 2, dbl_mb(k_occ), dbl_mb(k_beval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0,
     &       DISK_ERR)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_bmos), nbf)
c
      write(6,*) ' Alpha eigenvalues '
      call output(dbl_mb(k_aeval),1,nmo,1,1,nmo,1,1)
      write(6,*) ' Beta  eigenvalues '
      call output(dbl_mb(k_beval),1,nmo,1,1,nmo,1,1)
c$$$      write(6,*) ' Alpha MOs'
c$$$      call output(dbl_mb(k_amos),1,nbf,1,nmo,nbf,nmo,1)
c$$$      write(6,*) ' Beta  MOs'
c$$$      call output(dbl_mb(k_bmos),1,nbf,1,nmo,nbf,nmo,1)
c
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0, GA_ERR)
c
c     So now have nalpha, nbeta, nmo, a/beval, a/bmos
c
c     For brute force approach form the spin-orbitals in the order
c
c     beta occupied, alpha occupied, beta virtual, alpha virtual
c
      if (.not. ma_push_get(mt_dbl, nbf*nso,'sos', l_sos, k_sos))
     $     call errquit('ma sos', nbf*nso, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nso,'sos', l_seval, k_seval))
     $     call errquit('ma seval', nso, MA_ERR)
      if (.not. ma_push_get(mt_int, nso,'sospin', l_sospin, k_sospin))
     $     call errquit('ma sospin', nso, MA_ERR)
      call rjh_munge_mo_to_so(nbf, nmo, nso, nalpha, nbeta,
     $     dbl_mb(k_aeval), dbl_mb(k_beval), dbl_mb(k_amos),
     $     dbl_mb(k_bmos), dbl_mb(k_seval), dbl_mb(k_sos),
     $     int_mb(k_sospin))
c
      write(6,*) ' Spin-orbital eigenvalues '
      call output(dbl_mb(k_seval), 1, nso, 1, 1, nso, 1, 1)
c$$$      write(6,*) ' Spin-orbitals '
c$$$      call output(dbl_mb(k_sos), 1, nbf, 1, nso, nbf, nso, 1)
      write(6,*) ' Spin-orbital spins'
      do i = 1, nso
         write(6,11) i, int_mb(k_sospin+i-1)
 11      format(i5,2x,i2)
      end do
c
c     Now do a brute force transformation to the spin-orbit basis
c     zeroing out afterwards the components that are zero by spin
c
      if (.not. ma_push_get(mt_dbl, nso**4,'soint', l_soint, k_soint))
     $     call errquit('ma soint', nso**4, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nso**4,'tmp', l_tmp, k_tmp))
     $     call errquit('ma tmp', nso**4, MA_ERR)
      if (.not. ma_push_get(mt_dbl,nbf**4,'aoint',l_aoint,k_aoint))
     $     call errquit('allocation of AO integrals failed',nbf**4,
     &       MA_ERR)
c
      call rjh_all_ao_integrals(rtdb,basis,nbf,'dirac',dbl_mb(k_aoint))
      call rjh_uhf_trans(nbf, nso, int_mb(k_sospin),
     $     dbl_mb(k_sos), dbl_mb(k_aoint), 
     $     dbl_mb(k_tmp), dbl_mb(k_soint))
      if (.not. ma_chop_stack(l_tmp))call errquit('chopping from tmp',0,
     &       MA_ERR)
c
c     Compute the MP2 energy from the plain SO integrals
c
      call rjh_uhf_mp2(nocc, nso, dbl_mb(k_seval), dbl_mb(k_soint))
c
c     Antisymmetrize the SO integrals and write to file ASOINTS
c
      call rjh_so_asym(nso, dbl_mb(k_soint))
c
c     Make the first-order T2 amplitudes and recompute the MP2 energy
c     from them.  Write T2 to file T2.  Write a zero T1 to T1.
c
      if (.not. ma_push_get(mt_dbl, (nocc**2)*(nvirt**2),'t2', 
     $     l_t2, k_t2))
     $     call errquit('ma tmp', (nocc**2)*(nvirt**2),
     &       MA_ERR)
c
      call rjh_uhf_mp2_from_asoint(nocc, nvirt, nso, 
     $     dbl_mb(k_seval), dbl_mb(k_soint), dbl_mb(k_t2))
c
c     Now do some CC iterations
c
      if (.not. ma_push_get(mt_dbl, (nocc**2)*(nvirt**2),'t1', 
     $     l_t1, k_t1))
     $     call errquit('ma tmp', (nocc**2)*(nvirt**2),
     &       MA_ERR)
      call dfill(nvirt*nocc, 0.0d0, dbl_mb(k_t1), 1)
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,int_mb(k_sospin),
     $     dbl_mb(k_t1), dbl_mb(k_t2), 0)
c
      if (.not. ma_push_get(mt_dbl, (nocc**2)*(nvirt**2),'r2', 
     $     l_r2, k_r2))
     $     call errquit('ma tmp', (nocc**2)*(nvirt**2), MA_ERR)
      if (.not. ma_push_get(mt_dbl, (nocc**2)*(nvirt**2),'r1', 
     $     l_r1, k_r1))
     $     call errquit('ma tmp', (nocc**2)*(nvirt**2), MA_ERR)
      do iter = 1, 50
         call rjh_uhf_cc_energy(nocc, nvirt, nso,
     $        dbl_mb(k_t1), dbl_mb(k_t2),
     $        dbl_mb(k_seval), dbl_mb(k_soint))
         call rjh_uhf_cc_rhs(nocc, nvirt, nso, int_mb(k_sospin),
     $        dbl_mb(k_t1), dbl_mb(k_t2),
     $        dbl_mb(k_seval), dbl_mb(k_soint),
     $        dbl_mb(k_r1), dbl_mb(k_r2))
         call rjh_uhf_cc_update(nocc, nvirt, nso,
     $        dbl_mb(k_t1), dbl_mb(k_t2),
     $        dbl_mb(k_seval),
     $        dbl_mb(k_r1), dbl_mb(k_r2))
      end do
      call rjh_t1_t2_print(nocc,nvirt,nso,dbl_mb(k_t1),dbl_mb(k_t2))

c
c     Tidy up
c
      if (.not. ma_chop_stack(l_occ)) call errquit(' ma chop?', 0,
     &       MA_ERR)
      if (.not. bas_destroy(basis)) call errquit(' bas ?',0, BASIS_ERR)
      if (.not. geom_destroy(geom)) call errquit(' geom ?',0, GEOM_ERR)
c
      rjhtest = .true.
c
      end
      subroutine rjh_all_ao_integrals(rtdb, basis, nbf, order, ao)
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
      call rjh_do_all_ao_integrals(basis, dbl_mb(k_buf), dbl_mb(k_scr),
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     maxg2, maxs2, nbf, nsh, order, ao)
c
      if (.not. ma_chop_stack(l_scr)) call errquit('rjhallao: ma?',0,
     &       MA_ERR)
c
      call int_terminate
c
      end
      subroutine rjh_do_all_ao_integrals(
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
      subroutine rjh_munge_mo_to_so(nbf, nmo, nso, nalpha, nbeta,
     $     aeval, beval, amos, bmos, seval, sos, sospin)
      implicit none
#include "errquit.fh"
      integer nbf, nmo, nso, nalpha, nbeta, sospin(nso)
      double precision aeval(nmo), beval(nmo), seval(nso)
      double precision amos(nbf,nmo), bmos(nbf,nmo), sos(nbf,nso)
c
      integer i, j, cnt
c
      cnt = 0
c
      do i = 1, nbeta
         cnt = cnt + 1
         sospin(cnt)= -1
         seval(cnt) = beval(i)
         do j = 1, nbf
            sos(j,cnt) = bmos(j,i)
         end do
      end do
      do i = 1, nalpha
         cnt = cnt + 1
         sospin(cnt)= 1
         seval(cnt) = aeval(i)
         do j = 1, nbf
            sos(j,cnt) = amos(j,i)
         end do
      end do
      do i = nbeta+1,nmo
         cnt = cnt + 1
         sospin(cnt)= -1
         seval(cnt) = beval(i)
         do j = 1, nbf
            sos(j,cnt) = bmos(j,i)
         end do
      end do
      do i = nalpha+1,nmo
         cnt = cnt + 1
         sospin(cnt)= 1
         seval(cnt) = aeval(i)
         do j = 1, nbf
            sos(j,cnt) = amos(j,i)
         end do
      end do
c
      if (cnt .ne. nso) call errquit('munge: mismatch', cnt,
     &       UNKNOWN_ERR)
c
      end

      subroutine rjh_uhf_trans(nbf, nso, sospin, sos, ao, tmp, soint)
      implicit none
      integer nbf, nso
      integer sospin(nso)
      double precision sos(nbf,nso)
      double precision ao(nbf,nbf,nbf,nbf)
      double precision tmp(nso,nso,nso,nso)
      double precision soint(nso,nso,nso,nso)
c
      integer i, j, k, l        ! AO indices
      integer p, q, r, s        ! SO indices
c
c$$$      call dfill(nso**4, 0.0d0, soint, 1)
c$$$      do i = 1, nbf
c$$$         do j = 1, nbf
c$$$            do k = 1, nbf
c$$$               do l = 1, nbf
c$$$                  do s = 1, nso
c$$$                     do r = 1, nso
c$$$                        do q = 1, nso
c$$$                           do p = 1, nso
c$$$                              soint(p,q,r,s) = soint(p,q,r,s) + 
c$$$     $                             ao(i,j,k,l)*
c$$$     $                             sos(i,p)*sos(j,q)*sos(k,r)*sos(l,s)
c$$$                           end do
c$$$                        end do
c$$$                     end do
c$$$                  end do
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$c
c$$$      goto 999
c
      call dfill(nso**4, 0.0d0, tmp, 1)
      do i = 1, nbf
         do j = 1, nbf
            do k = 1, nbf
               do l = 1, nbf
                  do p = 1, nso
                     tmp(p,k,j,i) = tmp(p,k,j,i) + ao(l,k,j,i)*sos(l,p)
                  end do
               end do
            end do
         end do
      end do
c
      call dfill(nso**4, 0.0d0, soint, 1)
      do i = 1, nbf
         do j = 1, nbf
            do k = 1, nbf
               do q = 1, nso
                  do p = 1, nso
                     soint(p,q,j,i)=soint(p,q,j,i)+tmp(p,k,j,i)*sos(k,q)
                  end do
               end do
            end do
         end do
      end do
c
      call dfill(nso**4, 0.0d0, tmp, 1)
      do i = 1, nbf
         do j = 1, nbf
            do r = 1, nso
               do q = 1, nso
                  do p = 1, nso
                     tmp(p,q,r,i)=tmp(p,q,r,i)+soint(p,q,j,i)*sos(j,r)
                  end do
               end do
            end do
         end do
      end do
c
      call dfill(nso**4, 0.0d0, soint, 1)
      do i = 1, nbf
         do s = 1, nso
            do r = 1, nso
               do q = 1, nso
                  do p = 1, nso
                     soint(p,q,r,s)=soint(p,q,r,s)+tmp(p,q,r,i)*sos(i,s)
                  end do
               end do
            end do
         end do
      end do
c
c     Now zero out spin-forbidden terms ... DIRAC ORDERING <pq|rs>
c
      do s = 1, nso
         do r = 1, nso
            do q = 1, nso
               do p = 1, nso
                  if ( sospin(p).ne.sospin(r) .or. 
     $                 sospin(q).ne.sospin(s) ) soint(p,q,r,s)=0.0d0
               end do
            end do
         end do
      end do
c
c$$$ 999  write(6,*)
c$$$      write(6,*) ' Spin-orbit integrals '
c$$$      write(6,*)
c$$$      do s = 1, nso
c$$$         do r = 1, nso
c$$$            do q = 1, nso
c$$$               do p = 1, nso
c$$$                  if ( abs(soint(p,q,r,s)) .gt. 1e-6 ) 
c$$$     $                 write(6,7) p,q,r,s,
c$$$     $                 sospin(p),sospin(q),sospin(r),sospin(s), 
c$$$     $                 soint(p,q,r,s)
c$$$ 7                format(1x,4i5,2x,4i5,2x,f12.6)
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c
      end
      subroutine rjh_uhf_mp2(nocc, nso, seval, soint)
      implicit none
c
      integer nocc, nso
      double precision seval(nso), soint(nso,nso,nso,nso)
c
      integer i, j, a, b
      double precision energy, numerator, denominator
c
      energy = 0.0d0
      do i = 1, nocc
         do j = 1, i-1
            do a = nocc+1, nso
               do b = nocc+1, a
                  numerator = (soint(a,b,i,j)-soint(b,a,i,j))**2
                  denominator = seval(i)+seval(j)-seval(a)-seval(b)
                  energy = energy + numerator/denominator
               end do
            end do
         end do
      end do
c
      write(6,1) energy
 1    format(//' MP2 energy from SO integrals ', f20.10)
c
      end
      subroutine rjh_so_asym(nso, soint)
      implicit none
      integer nso
      double precision soint(nso,nso,nso,nso)
      double precision tmp
      integer i, j, k, l
c
      do l = 1, nso
         do k= 1, l
            do j = 1, nso
               do i = 1, nso
                  tmp = soint(i,j,k,l)-soint(i,j,l,k)
                  soint(i,j,k,l) =  tmp
                  soint(i,j,l,k) = -tmp
               end do
            end do
         end do
      end do
c
c$$$      open(1,file='ASOINTS',form='formatted',status='unknown')
c$$$      write(1,*) nso
c$$$      do l = 1, nso
c$$$         do k= 1, nso
c$$$            do j = 1, nso
c$$$               do i = 1, nso
c$$$                  if (abs(soint(i,j,k,l)).gt. 1d-9) 
c$$$     $                 write(1,10) i, j, k, l, soint(i,j,k,l)
c$$$ 10               format(4i5,f20.10)
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      close(1,status='keep')
c$$$      write(6,*) ' Wrote antisymmetrized spin-orbit integrals'
c
      end
      subroutine rjh_uhf_mp2_from_asoint(nocc, nvirt, nso, 
     $     seval, asoint, t2)
      implicit none
c
      integer nocc, nso, nvirt
      double precision seval(nso), asoint(nso,nso,nso,nso)
      double precision t2(nvirt,nvirt,nocc,nocc)
c
      integer i, j, a, b
      double precision energy
c
      energy = 0.0d0
      do i = 1, nocc
         do j = 1, nocc
            do a = 1, nvirt
               do b = 1, nvirt
                  t2(a,b,i,j) = asoint(a+nocc,b+nocc,i,j)/
     $                 (seval(i)+seval(j)-seval(a+nocc)-seval(b+nocc))
                  energy = energy + 
     $                 t2(a,b,i,j)*asoint(a+nocc,b+nocc,i,j)
               end do
            end do
         end do
      end do
c
      energy = energy / 4.0d0
c
      write(6,1) energy
 1    format(//' MP2 energy from ASO integrals ', f20.10)
c
c$$$      open(1,file='T2',form='formatted',status='unknown')
c$$$      write(1,*) nocc, nvirt
c$$$      do i = 1, nocc
c$$$         do j = 1, nocc
c$$$            do a = 1, nvirt
c$$$               do b = 1, nvirt
c$$$                  if (abs(t2(a,b,i,j)).gt.1d-9) 
c$$$     $                 write(1,10) i, j, a, b, t2(a,b,i,j)
c$$$ 10               format(4i5,f20.10)
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      close(1,status='keep')
c$$$      write(6,*) ' Wrote T2 amplitudes '
c$$$c
c$$$      open(1,file='EVALS',form='formatted',status='unknown')
c$$$      write(1,*) nso
c$$$      do i = 1, nso
c$$$         write(1,11) i, seval(i)
c$$$ 11      format(i5,f20.10)
c$$$      end do
c$$$      close(1,status='keep')
c$$$      write(6,*) ' Wrote eigenvalues '
c
      end
      subroutine rjh_uhf_cc_rhs(nocc,nvirt,nso,sospin,
     $     t1,t2,eval,g,r1,r2)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer nocc, nvirt, nso
      double precision t1(nvirt,nocc) ! Input amplitudes
      double precision t2(nvirt,nvirt,nocc,nocc)
      double precision eval(nso) ! Eigenvalues of spin-orbitals
      double precision g(nso,nso,nso,nso) ! Antisym, spin-orb integrals
      double precision r1(nvirt,nocc) ! Output vectors
      double precision r2(nvirt,nvirt,nocc,nocc)
      integer sospin(nso)
c
c     Given the input amplitudes t1 and t2 form the RHS of the
c     spin-orbital CCSD equations as given by Bartlett, MEST II p1067.
c
      integer l_(14), k_(14), size(14)
      character*8 name
      integer i
c
c     Allocate memory for the intermediates and zero them out
c
      size(1) = nocc*nvirt
      size(2) = nvirt*nvirt
      size(3) = nocc*nocc
      size(4) = nocc*nocc
      size(5) = (nvirt**3)*nocc
      size(6) = (nvirt**3)*nocc
      size(7) = nvirt**4
      size(8) = (nocc**3)*nvirt
      size(9) = nocc**4
      size(10)= (nocc**2)*(nvirt**2)
      size(11)= (nocc**2)*(nvirt**2)
      size(12)= (nocc**2)*(nvirt**2)
      size(13)= (nvirt**3)*nocc
      size(14)= (nocc**3)*nvirt
c
      do i = 1, 14
         write(name,'(''temp'',i4)')
         if (.not. ma_push_get(mt_dbl, size(i), name, l_(i), k_(i)))
     $        call errquit(' rjh_uhf_cc_rhs: error allocating ', i,
     &       MA_ERR)
         call dfill(size(i), 0.0d0, dbl_mb(k_(i)), 1)
      end do
c
c     Now make the intermediates
c
      call rjh_uhf_cc_inter(nocc,nvirt,nso,t1,t2,eval,g,
     $     dbl_mb(k_( 1)),dbl_mb(k_( 2)),dbl_mb(k_( 3)),dbl_mb(k_( 4)),
     $     dbl_mb(k_( 5)),dbl_mb(k_( 6)),dbl_mb(k_( 7)),dbl_mb(k_( 8)),
     $     dbl_mb(k_( 9)),dbl_mb(k_(10)),dbl_mb(k_(11)),dbl_mb(k_(12)),
     $     dbl_mb(k_(13)),dbl_mb(k_(14)))
c
c     Form the product
c
      call dfill(nocc*nvirt, 0.0d0, r1, 1)
      call dfill((nocc*nvirt)**2, 0.0d0, r2, 1)
      call rjh_uhf_cc_product(nocc,nvirt,nso,t1,t2,eval,g,
     $     dbl_mb(k_( 1)),dbl_mb(k_( 2)),dbl_mb(k_( 3)),dbl_mb(k_( 4)),
     $     dbl_mb(k_( 5)),dbl_mb(k_( 6)),dbl_mb(k_( 7)),dbl_mb(k_( 8)),
     $     dbl_mb(k_( 9)),dbl_mb(k_(10)),dbl_mb(k_(11)),dbl_mb(k_(12)),
     $     dbl_mb(k_(13)),dbl_mb(k_(14)), r1, r2, sospin)
*      call rjh_uhf_cc_product_dumb(nocc,nvirt,nso,t1,t2,eval,g,
*     $     r1, r2, sospin)
c
      if (.not. ma_chop_stack(l_(1))) call errquit
     $     ('rjh_cc_prod: ma corrupt', 0,
     &       MA_ERR)
c
      end
      subroutine rjh_uhf_cc_inter(nocc,nvirt,nso,t1,t2,eval,g,
     $     x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14)
      implicit none
#include "mafdecls.fh"
c
      integer nocc, nvirt, nso
      double precision t1(nvirt,nocc) ! Input amplitudes
      double precision t2(nvirt,nvirt,nocc,nocc)
      double precision eval(nso) ! Eigenvalues of spin-orbitals
      double precision g(nso,nso,nso,nso) ! Antisym, spin-orb integrals
c
c     In indexing the intermediates, the convention is upper indices
c     first, followed by lower indices, in order from left to right.
c
      double precision  x1(nocc,nvirt)
      double precision  x2(nvirt,nvirt)
      double precision  x3(nocc,nocc)
      double precision  x4(nocc,nocc)
      double precision  x5(nvirt,nvirt,nvirt,nocc)
      double precision  x6(nvirt,nvirt,nvirt,nocc)
      double precision  x7(nvirt,nvirt,nvirt,nvirt)
      double precision  x8(nocc,nocc,nocc,nvirt)
      double precision  x9(nocc,nocc,nocc,nocc)
      double precision x10(nocc,nvirt,nocc,nvirt)
      double precision x11(nocc,nvirt,nocc,nvirt)
      double precision x12(nocc,nvirt,nocc,nvirt)
      double precision x13(nvirt,nvirt,nocc,nvirt)
      double precision x14(nocc,nocc,nvirt,nocc)
c
      integer voff
      integer i, j, k, l, m, n, a, b, c, d, e, f
c
      voff = nocc
c
      do e = 1, nvirt
         do a = 1, nvirt
            do m = 1, nocc
               do i = 1, nocc
                  x1(i,a) = x1(i,a) + g(i,m,a+voff,e+voff)*t1(e,m)
               end do
            end do
         end do
      end do
c
      do b = 1, nvirt
         do a = 1, nvirt
            do e = 1, nvirt
               do m = 1, nocc
                  x2(a,b) = x2(a,b) + 
     $                 g(a+voff,m,b+voff,e+voff)*t1(e,m)
               end do
            end do
            do m = 1, nocc
               x2(a,b) = x2(a,b) - x1(m,b)*t1(a,m)
            end do
            do n = 1, nocc
               do m = 1, nocc
                  do e = 1, nvirt
                     x2(a,b) = x2(a,b) - 0.5d0*g(e+voff,b+voff,m,n)*
     $                    t2(e,a,m,n)
                  end do
               end do
            end do
         end do
      end do
c
      do j = 1, nocc
         do i = 1, nocc
            do m = 1, nocc
               do e = 1, nvirt
                  x3(j,i) = x3(j,i) + g(e+voff,j,m,i)*t1(e,m)
               end do
            end do
            do m = 1, nocc
               do f = 1, nvirt
                  do e = 1, nvirt
                     x3(j,i) = x3(j,i) + 0.5d0*g(e+voff,f+voff,m,i)*
     $                    t2(e,f,m,j)
                  end do
               end do
            end do
         end do
      end do
c
      do j = 1, nocc
         do i = 1, nocc
            x4(j,i) = x3(j,i)
            do e = 1, nvirt
               x4(j,i) = x4(j,i) + x1(i,e)*t1(e,j)
            end do
         end do
      end do
c
      do i = 1, nocc
         do b = 1, nvirt
            do c = 1, nvirt
               do a = 1, nvirt
                  x5(a,b,c,i) = g(a+voff,i,b+voff,c+voff)
                  do m = 1, nocc
                     x5(a,b,c,i) = x5(a,b,c,i) - 
     $                    0.5d0*g(m,i,b+voff,c+voff)*t1(a,m)
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, nocc
         do b = 1, nvirt
            do c = 1, nvirt
               do a = 1, nvirt
                  x6(a,b,c,i) = x5(a,b,c,i)
                  do m = 1, nocc
                     x6(a,b,c,i) = x6(a,b,c,i) - 
     $                    0.5d0*g(m,i,b+voff,c+voff)*t1(a,m)
                  end do
               end do
            end do
         end do
      end do
c
      do d= 1, nvirt
         do c = 1, nvirt
            do b = 1, nvirt
               do a = 1, nvirt
                  x7(a,b,c,d) = g(a+voff,b+voff,c+voff,d+voff)
               end do
            end do
            do m = 1, nocc
               do b = 1, nvirt
                  do a = 1, nvirt
                     x7(a,b,c,d) = x7(a,b,c,d) - 
     $                    (x5(a,c,d,m)*t1(b,m) - x5(b,c,d,m)*t1(a,m))
                  end do
               end do
            end do
         end do
      end do
c
      do a = 1, nvirt
         do j = 1, nocc
            do i = 1, nocc
               do k = 1, nocc
                  x8(k,i,j,a) = g(i,j,k,a+voff)
                  do e = 1, nvirt
                     x8(k,i,j,a) = x8(k,i,j,a) + 
     $                    0.5d0*g(e+voff,a+voff,i,j)*t1(e,k)
                  end do
               end do
            end do
         end do
      end do
c
      do j = 1, nocc
         do i = 1, nocc
            do l = 1, nocc
               do k = 1, nocc
                  x9(k,l,i,j) = g(k,l,i,j)
               end do
            end do
            do e = 1, nvirt
               do l = 1, nocc
                  do k = 1, nocc
                     x9(k,l,i,j) = x9(k,l,i,j) 
     $                    + x8(k,i,j,e)*t1(e,l) - x8(l,i,j,e)*t1(e,k)
                  end do
               end do
            end do
            do l = 1, nocc
               do k = 1, nocc
                  do f = 1, nvirt
                     do e = 1, nvirt
                        x9(k,l,i,j) = x9(k,l,i,j) + 
     $                       g(e+voff,f+voff,i,j)*t2(e,f,k,l)*0.5d0
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do b = 1, nvirt
         do i = 1, nocc
            do a = 1, nvirt
               do j = 1, nocc
                  x10(j,a,i,b) = g(i,a+voff,j,b+voff)
                  do m = 1, nocc
                     x10(j,a,i,b) = x10(j,a,i,b) - 
     $                    g(i,m,j,b+voff)*t1(a,m)
                  end do
                  do e = 1, nvirt
                     x10(j,a,i,b) = x10(j,a,i,b) + x5(a,b,e,i)*t1(e,j)
                  end do
               end do
            end do
         end do
      end do
c
      do b = 1, nvirt
         do i = 1, nocc
            do a = 1, nvirt
               do j = 1, nocc
                  x11(j,a,i,b) = x10(j,a,i,b)
                  do m = 1, nocc
                     x11(j,a,i,b) = x11(j,a,i,b) - 
     $                    g(m,i,b+voff,j)*t1(a,m)
                  end do
               end do
            end do
            do e = 1, nvirt
               do a = 1, nvirt
                  do j = 1, nocc
                     x11(j,a,i,b) = x11(j,a,i,b) + x6(a,b,e,i)*t1(e,j)
                  end do
               end do
            end do
         end do
      end do
c
      do b = 1, nvirt
         do i = 1, nocc
            do a = 1, nvirt
               do j = 1, nocc
                  x12(j,a,i,b) = x11(j,a,i,b)
               end do
            end do
            do m = 1, nocc
               do j = 1, nocc
                  do a = 1, nvirt
                     do e = 1, nvirt
                        x12(j,a,i,b) = x12(j,a,i,b) -
     $                       g(e+voff,b+voff,i,m)*t2(e,a,j,m)*0.5d0
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      do c = 1, nvirt
         do i = 1, nocc
            do b = 1, nvirt
               do a = 1, nvirt
                  x13(a,b,i,c) = g(a+voff,b+voff,c+voff,i)
                  do e = 1, nvirt
                     x13(a,b,i,c) = x13(a,b,i,c) + 
     $                    0.5d0*g(e+voff,c+voff,b+voff,a+voff)*t1(e,i)
                  end do
               end do
            end do
            do m = 1, nocc
               do b = 1, nvirt
                  do a = 1, nvirt
                     x13(a,b,i,c) = x13(a,b,i,c) - 
     $                    (x10(i,a,m,c)*t1(b,m)-x10(i,b,m,c)*t1(a,m))
                  end do
               end do
            end do
         end do
      end do
c
      do i = 1, nocc
         do a = 1, nvirt
            do m = 1, nocc
               do k = 1, nocc
                  do j= 1, nocc
                     x14(j,k,a,i) = x14(j,k,a,i) + 
     $                    g(j,k,i,m)*t1(a,m)
                  end do
               end do
            end do
         end do
      end do
c
      end
      subroutine rjh_uhf_cc_product(nocc,nvirt,nso,t1,t2,eval,g,
     $     x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,
     $     r1, r2, sospin)
      implicit none
#include "mafdecls.fh"
c
      integer nocc, nvirt, nso
      double precision t1(nvirt,nocc) ! Input amplitudes
      double precision t2(nvirt,nvirt,nocc,nocc)
      double precision eval(nso) ! Eigenvalues of spin-orbitals
      double precision g(nso,nso,nso,nso) ! Antisym, spin-orb integrals
c
c     Final step in computing the CC RHS
c
c     In indexing the intermediates, the convention is upper indices
c     first, followed by lower indices, in order from left to right.
c
      double precision  x1(nocc,nvirt)
      double precision  x2(nvirt,nvirt)
      double precision  x3(nocc,nocc)
      double precision  x4(nocc,nocc)
      double precision  x5(nvirt,nvirt,nvirt,nocc)
      double precision  x6(nvirt,nvirt,nvirt,nocc)
      double precision  x7(nvirt,nvirt,nvirt,nvirt)
      double precision  x8(nocc,nocc,nocc,nvirt)
      double precision  x9(nocc,nocc,nocc,nocc)
      double precision x10(nocc,nvirt,nocc,nvirt)
      double precision x11(nocc,nvirt,nocc,nvirt)
      double precision x12(nocc,nvirt,nocc,nvirt)
      double precision x13(nvirt,nvirt,nocc,nvirt)
      double precision x14(nocc,nocc,nvirt,nocc)
      double precision  r1(nvirt,nocc)
      double precision  r2(nvirt,nvirt,nocc,nocc)
c
      integer sospin(nso)
c
      integer voff
      integer i, j, m, n, a, b, e, f
c
      voff = nocc
c
      do i = 1, nocc
         do a = 1, nvirt
            do e = 1, nvirt
               r1(a,i) = r1(a,i) + x2(a,e)*t1(e,i)
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,2)
            do m = 1, nocc
               r1(a,i) = r1(a,i) - x3(i,m)*t1(a,m)
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,21)
            do m = 1, nocc
               do e = 1, nvirt
                  r1(a,i) = r1(a,i) + g(e+voff,i,m,a+voff)*t1(e,m)
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,3)
            do m = 1, nocc
               do e = 1, nvirt
                  r1(a,i) = r1(a,i) + x1(m,e)*t2(e,a,m,i)
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,33)
            do m = 1, nocc
               do f = 1, nvirt
                  do e = 1, nvirt
                     r1(a,i) = r1(a,i) +
     $                    0.5d0*g(e+voff,f+voff,m,a+voff)*t2(e,f,m,i)
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,4)
            do n = 1, nocc
               do m = 1, nocc
                  do e = 1, nvirt
                     r1(a,i) = r1(a,i) - 
     $                    0.5d0*g(e+voff,i,m,n)*t2(e,a,m,n)
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,5)
         end do
      end do
c
      do j = 1, nocc
         do i = 1, nocc
            do b = 1, nvirt
               do a = 1, nvirt
                  r2(a,b,i,j) = g(a+voff,b+voff,i,j)
                  do e = 1, nvirt
                     r2(a,b,i,j) = r2(a,b,i,j) + 
     $                    x2(b,e)*t2(a,e,i,j) - x2(a,e)*t2(b,e,i,j)
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,6)
            do m = 1, nocc
               do b = 1, nvirt
                  do a = 1, nvirt
                     r2(a,b,i,j) = r2(a,b,i,j) -
     $                    (x4(j,m)*t2(a,b,i,m) - x4(i,m)*t2(a,b,j,m))
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,7)
            do f = 1, nvirt
               do e = 1, nvirt
                  do b = 1, nvirt
                     do a = 1, nvirt
                        r2(a,b,i,j) = r2(a,b,i,j) + 
     $                       0.5d0*x7(a,b,e,f)*t2(e,f,i,j)
                     end do
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,8)
            do n = 1, nocc
               do m = 1, nocc
                  do b = 1, nvirt
                     do a = 1, nvirt
                        r2(a,b,i,j) = r2(a,b,i,j) +
     $                       0.5d0*x9(i,j,m,n)*t2(a,b,m,n)
                     end do
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,9)
            do e = 1, nvirt
               do m = 1, nocc
                  do b = 1, nvirt
                     do a = 1, nvirt
                        r2(a,b,i,j) = r2(a,b,i,j)
     $                       - x12(i,b,m,e)*t2(a,e,m,j)
     $                       + x12(j,b,m,e)*t2(a,e,m,i)
     $                       + x12(i,a,m,e)*t2(b,e,m,j)
     $                       - x12(j,a,m,e)*t2(b,e,m,i)
                     end do
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,10)
            do e = 1, nvirt
               do b = 1, nvirt
                  do a = 1, nvirt
                     r2(a,b,i,j) = r2(a,b,i,j)
     $                    + x13(a,b,j,e)*t1(e,i) - x13(a,b,i,e)*t1(e,j)
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,11)
            do m = 1, nocc
               do b = 1, nvirt
                  do a = 1, nvirt
                     r2(a,b,i,j) = r2(a,b,i,j)
     $                    -(x14(i,j,b,m)*t1(a,m) + x14(i,j,a,m)*t1(b,m))
                  end do
               end do
            end do
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,12)
         end do
      end do
c
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,1)

c
      end
      subroutine rjh_uhf_cc_energy(nocc,nvirt,nso,t1,t2,eval,g)
      implicit none
#include "mafdecls.fh"
c
      integer nocc, nvirt, nso
      double precision t1(nvirt,nocc) ! Input amplitudes
      double precision t2(nvirt,nvirt,nocc,nocc)
      double precision eval(nso) ! Eigenvalues of spin-orbitals
      double precision g(nso,nso,nso,nso) ! Antisym, spin-orb integrals
c
      double precision energy
      integer i, j, a, b
c
      energy = 0.0d0
      do j = 1, nocc
         do i = 1, nocc
            do b = 1, nvirt
               do a = 1, nvirt
                  energy = energy + g(a+nocc,b+nocc,i,j)*
     $                 (t2(a,b,i,j) + t1(a,i)*t1(b,j) - t1(b,i)*t1(a,j))
                  if (abs(t2(a,b,i,j)).gt.1e-6 .and. 
     $                 abs(g(a+nocc,b+nocc,i,j)).lt.1e-6) then
                     write(6,*) ' ?? ', i, j, a, b, 
     $                    g(a+nocc,b+nocc,i,j), t2(a,b,i,j)
                  end if
               end do
            end do
         end do
      end do
      energy = energy * 0.25d0
c
      write(6,1) energy
 1    format(' CCSD energy = ', f20.10)
c
      end
      subroutine rjh_uhf_cc_update(nocc,nvirt,nso,t1,t2,eval,r1,r2)
      implicit none
#include "mafdecls.fh"
c
      integer nocc, nvirt, nso
      double precision t1(nvirt,nocc) ! Input amplitudes
      double precision t2(nvirt,nvirt,nocc,nocc)
      double precision eval(nso) ! Eigenvalues of spin-orbitals
      double precision r1(nvirt,nocc) ! Input amplitudes
      double precision r2(nvirt,nvirt,nocc,nocc)
c
      double precision denom, t1norm, t2norm, damp
      integer i, j, a, b
c
      damp = 0.0d0
c
      t2norm = 0.0d0
      do j = 1, nocc
         do i = 1, nocc
            do b = 1, nvirt
               do a = 1, nvirt
                  denom = eval(i) + eval(j) - 
     $                 eval(a+nocc) - eval(b+nocc)
                  r2(a,b,i,j) = r2(a,b,i,j) / denom
                  t2norm = t2norm + (r2(a,b,i,j)-t2(a,b,i,j))**2
                  t2(a,b,i,j) = damp*t2(a,b,i,j) + 
     $                 (1.0d0-damp)*r2(a,b,i,j)
               end do
            end do
         end do
      end do
c
      t1norm = 0.0d0
      do i = 1, nocc
         do a = 1, nvirt
            denom = eval(i) - eval(a+nocc)
            r1(a,i) = r1(a,i) / denom
            t1norm = t1norm + (r1(a,i)-t1(a,i))**2
            t1(a,i) = damp*t1(a,i) + (1.0d0-damp)*r1(a,i)
*            t1(a,i) = 0.0d0
         end do
      end do
c
      write(6,1) sqrt(t1norm), sqrt(t2norm)
 1    format('  Residuals:  t1 ',1p,d9.1,'   t2 ',d9.1)
c
      end
      subroutine rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,code)
      implicit none
#include "errquit.fh"
      integer nocc, nvirt, nso, code, sospin(nso)
      double precision r1(nvirt,nocc)
      double precision r2(nvirt,nvirt,nocc,nocc)
      integer i, j, a, b
c
      do i = 1, nocc
         do a = 1, nvirt
            if (abs(r1(a,i)).gt.1d-10 .and. 
     $           sospin(i).ne.sospin(a+nocc)) then
               write(6,*) ' Spin forbidden ', i, a, r1(a,i), 
     $              sospin(i), sospin(nocc+a)
               call errquit(' bad t1 ', code, UNKNOWN_ERR)
            end if
         end do
      end do
      do j = 1, nocc
         do i = 1, nocc
            do b = 1, nvirt
               do a = 1, nvirt
                  if ( abs(r2(a,b,i,j)).gt.1d-10  .and.
     $                 (sospin(i)*sospin(j)).ne.
     $                 (sospin(a+nocc)*sospin(b+nocc))) 
     $                 then
                     write(6,*) ' Spin forbidden ', i, j, a, b, 
     $                    r2(a,b,i,j)
                     call errquit(' bad t2 ', code, UNKNOWN_ERR)
                  end if
                  if ((i.eq.j .or. a.eq.b) .and. 
     $                 abs(r2(a,b,i,j)).gt.1d-10) then
                     write(6,*) ' Non-zero diagonal t2 ',
     $                    i, j, a, b, r2(a,b,i,j)
                     call errquit(' bad t2 ', code, UNKNOWN_ERR)
                  end if
               end do
            end do
         end do
      end do
c
      end
      subroutine rjh_t1_t2_print(nocc,nvirt,nso,r1,r2)
      implicit none
      integer nocc, nvirt, nso
      double precision r1(nvirt,nocc)
      double precision r2(nvirt,nvirt,nocc,nocc)
      integer i, j, a, b
      do i = 1, nocc
         do j = 1, nocc
            do b = 1, nvirt
               do a = 1, nvirt
                  if(abs(r2(a,b,i,j)).gt.1e-6) then
                     write(6,1) a, b, i, j, r2(a, b, i, j)
 1                   format(4i5,f12.6)
                  end if
               end do
            end do
         end do
      end do
      do i = 1, nocc
         do a = 1, nvirt
            if(abs(r1(a,i)).gt.1e-6) then
               write(6,2) a, i, r1(a,i)
 2             format(2i5,f12.6)
            end if
         end do
      end do
c
      end
      subroutine rjh_uhf_cc_product_dumb(nocc,nvirt,nso,t1,t2,eval,g,
     $     r1, r2, sospin)
      implicit none
#include "mafdecls.fh"
c
      integer nocc, nvirt, nso
      double precision t1(nvirt,nocc) ! Input amplitudes
      double precision t2(nvirt,nvirt,nocc,nocc)
      double precision r1(nvirt,nocc) ! Input amplitudes
      double precision r2(nvirt,nvirt,nocc,nocc)
      double precision eval(nso) ! Eigenvalues of spin-orbitals
      double precision g(nso,nso,nso,nso) ! Antisym, spin-orb integrals
      integer sospin(nso)
c
c     Final step in computing the CC RHS ... DUMB version CCD only
c
      integer i, j, k, l, a, b, c, d, voff
c
      voff = nocc
c
      do j = 1, nocc
         do i = 1, nocc
            do b = 1, nvirt
               do a = 1, nvirt
                  r2(a,b,i,j) = g(a+voff,b+voff,i,j)
                  do d = 1, nvirt
                     do c = 1, nvirt
                        r2(a,b,i,j) = r2(a,b,i,j) + 0.5d0*
     $                       g(c+voff,d+voff,a+voff,b+voff)*
     $                       t2(c,d,i,j)
                     enddo
                  enddo
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,81)
                  do l = 1, nocc
                     do k = 1, nocc
                        r2(a,b,i,j) = r2(a,b,i,j) + 0.5d0*
     $                       g(k,l,i,j)*t2(a,b,k,l)
                     enddo
                  enddo
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,82)
                  do k = 1, nocc
                     do c = 1, nvirt
                        r2(a,b,i,j) = r2(a,b,i,j)
     $                       + g(c+voff,j,a+voff,k)*t2(b,c,i,k) 
     $                       - g(c+voff,j,b+voff,k)*t2(a,c,i,k) 
     $                       - g(c+voff,i,a+voff,k)*t2(b,c,j,k) 
     $                       + g(c+voff,i,b+voff,k)*t2(a,c,j,k) 
                     enddo
                  enddo
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,83)
                  do k = 1, nocc
                     do l = 1, nocc
                        do d = 1, nvirt
                           do c = 1, nvirt
                              r2(a,b,i,j) = r2(a,b,i,j) + 
     $                             g(c+voff,d+voff,k,l) * (
     $                             + t2(a,c,i,k)*t2(d,b,l,j)
     $                             - t2(b,c,i,k)*t2(d,a,l,j)

     $                             - t2(c,d,k,i)*t2(a,b,l,j)*0.5d0
     $                             + t2(c,d,k,j)*t2(a,b,l,i)*0.5d0

     $                             - t2(c,a,k,l)*t2(d,b,i,j)*0.5d0
     $                             + t2(c,b,k,l)*t2(d,a,i,j)*0.5d0

     $                             + t2(c,d,i,j)*t2(a,b,k,l)*0.25d0)
                           enddo
                        enddo
                     enddo
                  enddo
      call rjh_t1_t2_spin_check(nocc,nvirt,nso,sospin,r1,r2,84)
               enddo
            enddo
         enddo
      enddo
c
      end

