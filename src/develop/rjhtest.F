      logical function rjhtest(rtdb)
*     
*     $Id: rjhtest.F,v 1.23 2000-02-04 19:53:16 mg201 Exp $
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
      integer l_t2, k_t2
      integer i
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
      if (.not.int_normalize(rtdb,basis))
     $           call errquit('scf:int_normalize failed', 0)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0)
c
c     Read the MO vectors and evals from a UHF calculation
c
      call util_file_name('movecs',.false.,.false.,movecs)
      if (.not. movecs_read_header(movecs, title, name_of_basis,
     $     scftype, nbf_file, nsets, nmo_file, 2)) call errquit
     $     ('rjhtest: failed to read movecs header',911)
      if (nsets.ne.2) call errquit('UHF ONLY!',0)
      write(6,*) ' Read movecs header from ', movecs
      write(6,*) ' Job title :                ', 
     $     title(1:inp_strlen(title))
      write(6,*) ' Basis name:                ', 
     $     name_of_basis(1:inp_strlen(name_of_basis))
      nmo = nmo_file(1)
      if (nmo_file(1).ne.nmo_file(2)) call errquit('Uh?',0)
      if (.not. rtdb_get(rtdb, 'scf:nalpha', mt_int, 1, nalpha))
     $     call errquit('nalpha?',0)
      if (.not. rtdb_get(rtdb, 'scf:nbeta', mt_int, 1, nbeta))
     $     call errquit('nbeta?',0)
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
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'tmp', 0, 0, g_tmp))
     &     call errquit('scf_v_g: tmp', 0)
      if (.not. ma_push_get(mt_dbl, nbf,'occ',l_occ, k_occ))
     $     call errquit('ma occ', nbf)
      if (.not. ma_push_get(mt_dbl, nbf,'aeval',l_aeval, k_aeval))
     $     call errquit('ma eval', nbf)
      if (.not. ma_push_get(mt_dbl, nbf,'beval',l_beval, k_beval))
     $     call errquit('ma eval', nbf)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'amos', l_amos, k_amos))
     $     call errquit('ma amos', nbf*nbf)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'bmos', l_bmos, k_bmos))
     $     call errquit('ma bmos', nbf*nbf)
c
      if (.not. movecs_read(movecs, 1, dbl_mb(k_occ), dbl_mb(k_aeval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0)
      call ga_get(g_tmp, 1, nbf, 1, nmo, dbl_mb(k_amos), nbf)
      if (.not. movecs_read(movecs, 2, dbl_mb(k_occ), dbl_mb(k_beval), 
     $     g_tmp)) call errquit('movecs_read of amos failed ',0)
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
      if (.not. ga_destroy(g_tmp)) call errquit(' ga bad?',0)
c
c     So now have nalpha, nbeta, nmo, a/beval, a/bmos
c
c     For brute force approach form the spin-orbitals in the order
c
c     beta occupied, alpha occupied, beta virtual, alpha virtual
c
      if (.not. ma_push_get(mt_dbl, nbf*nso,'sos', l_sos, k_sos))
     $     call errquit('ma sos', nbf*nso)
      if (.not. ma_push_get(mt_dbl, nso,'sos', l_seval, k_seval))
     $     call errquit('ma seval', nso)
      if (.not. ma_push_get(mt_int, nso,'sospin', l_sospin, k_sospin))
     $     call errquit('ma sospin', nso)
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
     $     call errquit('ma soint', nso**4)
      if (.not. ma_push_get(mt_dbl, nso**4,'tmp', l_tmp, k_tmp))
     $     call errquit('ma tmp', nso**4)
      if (.not. ma_push_get(mt_dbl,nbf**4,'aoint',l_aoint,k_aoint))
     $     call errquit('allocation of AO integrals failed',nbf**4)
c
      call rjh_all_ao_integrals(rtdb,basis,nbf,'dirac',dbl_mb(k_aoint))
      call rjh_uhf_trans(nbf, nso, int_mb(k_sospin),
     $     dbl_mb(k_sos), dbl_mb(k_aoint), 
     $     dbl_mb(k_tmp), dbl_mb(k_soint))
      if (.not. ma_chop_stack(l_tmp))call errquit('chopping from tmp',0)
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
     $     call errquit('ma tmp', (nocc**2)*(nvirt**2))
c
      call rjh_uhf_mp2_from_asoint(nocc, nvirt, nso, 
     $     dbl_mb(k_seval), dbl_mb(k_soint), dbl_mb(k_t2))
c
c     Tidy up
c
      if (.not. ma_chop_stack(l_occ)) call errquit(' ma chop?', 0)
      if (.not. bas_destroy(basis)) call errquit(' bas ?',0)
      if (.not. geom_destroy(geom)) call errquit(' geom ?',0)
      rjhtest = .true.
c
      end
      subroutine rjh_all_ao_integrals(rtdb, basis, nbf, order, ao)
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
      call rjh_do_all_ao_integrals(basis, dbl_mb(k_buf), dbl_mb(k_scr),
     $     int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $     maxg2, maxs2, nbf, nsh, order, ao)
c
      if (.not. ma_chop_stack(l_scr)) call errquit('rjhallao: ma?',0)
c
      call int_terminate
c
      end
      subroutine rjh_do_all_ao_integrals(
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
      subroutine rjh_munge_mo_to_so(nbf, nmo, nso, nalpha, nbeta,
     $     aeval, beval, amos, bmos, seval, sos, sospin)
      implicit none
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
      if (cnt .ne. nso) call errquit('munge: mismatch', cnt)
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
      open(1,file='ASOINTS',form='formatted',status='unknown')
      write(1,*) nso
      do l = 1, nso
         do k= 1, nso
            do j = 1, nso
               do i = 1, nso
                  if (abs(soint(i,j,k,l)).gt. 1d-9) 
     $                 write(1,10) i, j, k, l, soint(i,j,k,l)
 10               format(4i5,f20.10)
               end do
            end do
         end do
      end do
      close(1,status='keep')
      write(6,*) ' Wrote antisymmetrized spin-orbit integrals'
c
      end
      subroutine rjh_uhf_mp2_from_asoint(nocc, nvirt, nso, 
     $     seval, asoint, t2)
      implicit none
c
      integer nocc, nso, nvirt
      double precision seval(nso), asoint(nso,nso,nso,nso)
      double precision t2(nocc,nocc,nvirt,nvirt)
c
      integer i, j, a, b
      double precision energy
c
      energy = 0.0d0
      do i = 1, nocc
         do j = 1, nocc
            do a = 1, nvirt
               do b = 1, nvirt
                  t2(i,j,a,b) = asoint(i,j,a+nocc,b+nocc)/
     $                 (seval(i)+seval(j)-seval(a+nocc)-seval(b+nocc))
                  energy = energy + 
     $                 t2(i,j,a,b)*asoint(i,j,a+nocc,b+nocc)
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
      open(1,file='T2',form='formatted',status='unknown')
      write(1,*) nocc, nvirt
      do i = 1, nocc
         do j = 1, nocc
            do a = 1, nvirt
               do b = 1, nvirt
                  if (abs(t2(i,j,a,b)).gt.1d-9) 
     $                 write(1,10) i, j, a, b, t2(i,j,a,b)
 10               format(4i5,f20.10)
               end do
            end do
         end do
      end do
      close(1,status='keep')
      write(6,*) ' Wrote T2 amplitudes '
c
      open(1,file='EVALS',form='formatted',status='unknown')
      write(1,*) nso
      do i = 1, nso
         write(1,11) i, seval(i)
 11      format(i5,f20.10)
      enddo
      close(1,status='keep')
      write(6,*) ' Wrote eigenvalues '
c
      end
