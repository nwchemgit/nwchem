      subroutine uccsdt_eterm(qO_handles,qV_handles)
c
c$Id: eterm.F,v 1.5 2002-10-08 00:17:48 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer qO_handles(0:7,8), qV_handles(0:7,8) 
      integer fO_handles(2,0:7), fV_handles(2,0:7)
C
C     GENERATION OF E-TERM
C

c (1) R(E1) aa = +(1-(ef)) t(me,na) [ f(fa)V - 0.5*q_1_O - q_2_O] m=e=n=f=a=alpha
c (1) R(E1) ab =           t(ma,nf) [ f(ea)V - 0.5*q_1_O - q_2_O] m=e=a=alpha n=f=beta
c (2) R(E1) bb = +(1-(ef)) t(me,na) [ f(fa)V - 0.5*q_4_O - q_5_O] m=e=n=f=a=beta 
c (3) R(E2) aa = -(1-(mn)) t(me,if) [ f(in)O + 0.5*q_1_V + q_2_V] m=e=n=f=i=alpha
c (3) R(E4) ab =           t(ie,nf) [ f(im)O + 0.5*q_1_V + q_2_V] i=e=m=alpha n=f=beta
c (4) R(E2) bb = -(1-(mn)) t(me,if) [ f(in)O + 0.5*q_4_V + q_5_V] m=e=n=f=i=beta      
c (4) R(E3) ab =           t(me,if) [ f(in)O + 0.5*q_4_V + q_5_V] m=e=alpha i=n=f=beta
c (5) R(E2) ab =           t(me,na) [ f(af)V - 0.5*q_4_O - q_8_O] m=e=alpha a=n=f=beta
c
c     first combine terms in [] and them "matmul" with t's
c     (1) is in qO(2,sym)
c     (2) is in qO(5,sym)
c     (3) is in qV(2,sym)
c     (4) is in qV(5,sym)
c     (5) is in qO(8,sym)
c
c     Who has the f blocks, they should be in format, let's just make them
c
      call make_fock_mo(fO_handles,fV_handles)
c
      do symf = 0, nir-1
         syma = symf
         call ga_add(-1.0d0,qO(2,symf),-0.5d0,qO(1,symf),qO(2,symf))
         call ga_add( 1.0d0,qO(2,symf), 1.0d0,fO(a,a,symf),qO(2,symf))   
         call ga_add(-1.0d0,qO(5,symf),-0.5d0,qO(4,symf),qO(5,symf))
         call ga_add( 1.0d0,qO(5,symf), 1.0d0,fO(b,b,symf),qO(5,symf))   
         call ga_add( 1.0d0,qV(2,symf), 0.5d0,qV(1,symf),qV(2,symf))
         call ga_add( 1.0d0,qV(2,symf), 1.0d0,fV(a,a,symf),qV(2,symf))  
         call ga_add( 1.0d0,qV(5,symf), 0.5d0,qV(4,symf),qV(5,symf))
         call ga_add( 1.0d0,qV(5,symf), 1.0d0,fV(b,b,symf),qV(5,symf)) 
         call ga_add(-1.0d0,qO(8,symf),-0.5d0,qO(4,symf),qO(8,symf))
         call ga_add( 1.0d0,qO(8,symf), 1.0d0,fO(a,a,symf),qO(8,symf))  
      end do
      do symnf = 0, nir-1
         symme = symnf
         symna = symnf
         alo_a = v_sym(1,0,1)
         ahi_a = v_sym(2,nir-1,1)
         alo_b = v_sym(1,0,2)
         ahi_b = v_sym(2,nir-1,2)
         call uccsdt_ampfile_read_t2(D_AMP,1,1,1,1,symna,alo_a,ahi_a,
     &                               g_t2_mena_pure,.true.,'block')
C23456789012345678901234567890123456789012345678901234567890123456789012
         do symf = 0, nir-1
            syma = symf
            symn = ieor(symf,symnf)
            symi = symn
            if (.not. ma_push_get(mt_dbl,no_sym(symn,1),'n array',
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0)
            do a = v_sym(1,syma,1), v_sym(2,syma,1)
               t_ilo = 1
               t_ihi = ov_len(symnf,1,1)
               t_jlo = ov_off(a,symnf,1,1)+1
               t_jhi = t_jlo + no_sym(symi,1)
c              matmul(g_t2_mena_pure,qO(2,symn)) into R(E1) aa   #(1)
               call ga_get(qO(symn,2),o_sym(1,symn,1),o_sym(2,symn,1),
     &                     a,a,dbl_mb(k_n),no_sym(symn,1)
               do n = 1, no_sym(symn,1)
                  t_jlo = ov_off(a,symnf,1,1)+1
                  t_jhi = t_jlo + no_sym(symi,1)
                  call ga_add_patch(dbl_mb(k_n+n-1),g_t2_mena_pure,
     &                              t_ilo,t_hi,t_jlo,t_jhi,1.0d0,
     &                              Raa,t_ilo,t_ihi,t_flo,t_fhi,
     &                              Raa,t_ilo,t_ihi,t_flo,t_fhi)
               enddo
c              matmul(g_t2_mena_pure,qV(2,symn)) into R(E2) aa   #(3)
               q_ijlo = o_sym(1,symi,1)
               q_ijhi = o_sym(2,symi,1)
               call ga_matmul_patch('n','n',1.0d0,1.0d0,
     &                g_t2_mena_pure,t_ilo,t_iho,t_jlo,t_jhi,
     &                qV(symn,2),q_ijlo,q_ijhi,q_ijlo,q_ijhi,
     &                g_raa,t_ilo,t_iho,t_jlo,t_jhi)
            end do
            if (.not. ma_pop_stack(l_n)) call
     &          errquit('eterm: n dealloc failed',0)
         end do 
         if (.not. ga_destroy(g_t2_mena_pure)) call 
     &       errquit('eterm: g_t2_mena_pure aa dealloc failed',0)
c
         call uccsdt_ampfile_read_t2(D_AMP,2,2,2,2,symna,alo_b,ahi_b,
     &                               g_t2_mena_pure,.true.,'block')
         do symf = 0, nir-1
            syma = symf
            symn = ieor(symf,symnf)
            symi = symn
            if (.not. ma_push_get(mt_dbl,no_sym(symn,1),'n array',
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0)
            do a = v_sym(1,syma,1), v_sym(2,syma,1)
               t_ilo = 1
               t_ihi = ov_len(symnf,2,2)
               t_jlo = ov_off(a,symnf,2,2)+1
               t_jhi = t_jlo + no_sym(symi,2)
c              matmul(g_t2_mena_pure,qO(5,symf)) into R(E1) bb   #(2)
               call ga_get(qO(symn,5),o_sym(1,symn,1),o_sym(2,symn,1),
     &                     a,a,dbl_mb(k_n),no_sym(symn,1)
               do n = 1, no_sym(symn,2)
                  t_jlo = ov_off(a,symnf,2,2)+1
                  t_jhi = t_jlo + no_sym(symi,2)
                  call ga_add_patch(dbl_mb(k_n+n-1),g_t2_mena_pure,
     &                              t_ilo,t_hi,t_jlo,t_jhi,1.0d0,
     &                              Raa,t_ilo,t_ihi,t_flo,t_fhi,
     &                              Raa,t_ilo,t_ihi,t_flo,t_fhi)
               enddo
c              matmul(g_t2_mena_pure,qV(5,symn)) into R(E2) bb   #(4)
               q_ijlo = o_sym(1,symi,2)
               q_ijhi = o_sym(2,symi,2)
               call ga_matmul_patch('n','n',1.0d0,1.0d0,
     &                g_t2_mena_pure,t_ilo,t_iho,t_jlo,t_jhi,
     &                qV(symn,5),q_ijlo,q_ijhi,q_ijlo,q_ijhi,
     &                g_raa,t_ilo,t_iho,t_jlo,t_jhi)
            end do
            if (.not. ma_pop_stack(l_n)) call
     &          errquit('eterm: n dealloc failed',0)
         end do
         if (.not. ga_destroy(g_t2_mena_pure)) call 
     &       errquit('eterm: g_t2_mena_pure bb dealloc failed',0)
c
         call uccsdt_ampfile_read_t2(D_AMP,1,1,2,2,symna,alo_b,ahi_b,
     &                               g_t2_mena_pure,.true.,'block')
         do symf = 0, nir-1
            syma = symf
            symn = ieor(symf,symnf)
            symi = symn
            do loop over a
               matmul(g_t2_mena_mix,qO(2,syma)) into R(E1) ab   #(1)
            end do  
            do loop over a
               matmul(g_t2_mena_mix,qV(2,symn)) into R(E4) ab   #(3)
            end do
            do a = v_sym(1,syma,2), v_sym(2,syma,2)
c              matmul(g_t2_mena_mix,qV(5,symn)) into R(E3) ab   #(4)
               t_ilo = 1
               t_ihi = ov_len(symna,1,1)
               t_jlo = ov_off(a,symna,2,2)+1
               t_jhi = t_jlo + no_sym(symi,2) 
               q_ijlo = o_sym(1,symi,2)
               q_ijhi = o_sym(2,symi,2)
               call ga_matmul_patch('n','n',1.0d0,1.0d0,
     &                g_t2_mena_pure,t_ilo,t_iho,t_jlo,t_jhi,
     &                qV(2,symn),q_ijlo,q_ijhi,q_ijlo,q_ijhi,
     &                g_raa,t_ilo,t_iho,t_jlo,t_jhi)
            end do
            do loop over i
               matmul(g_t2_mena_mix,qV(8,symn)) into R(E2) ab   #(5)
            end do
         end do
         if (.not. ga_destroy(g_t2_mena_mix)) call 
     &       errquit('eterm: g_t2_mena_mix ab dealloc failed',0)
      end do

      do symf = 0, nir-1
         do inum = 1, 8
            if (qO_handle(inum,symf) .gt. -1) then
                if (.not. ga_destroy(qO_handle(inum,symf))
     &             call errquit('eterm: qQ dealloc failed',inum)
            if (qV_handle(inum,symf) .gt. -1) then
                if (.not. ga_destroy(qV_handle(inum,symf))
     &             call errquit('eterm: qV dealloc failed',inum)
            end if
         end do
      end do
c
      end
c
      subroutine make_fock_mo(nmo, noa, nob, ha, iaa, iab, fa)
      implicit none
      integer nmo, noa, nob
      double precision ha(nmo,nmo)
      double precision iaa(nmo,nmo,nmo,nmo), iab(nmo,nmo,nmo,nmo)
      double precision fa(nmo,nmo)
C
C     MAKE FOCK MATRIX ALPHA AND BETA IN AO BASIS AND TRANSFORM WITH
C     G_PART(LEFT SIDE) AND G_HOLE(RIGHT SIDE) TO FOCK MATRIX IN MO BASIS
c
c     fij = hij + <ik||jk> = hij + <ik||jk>(same spin) + <ik|jk>(diff spin)
c
c     routine taken from test code, assumes full matrices, but we do have
c     the integrals in GA's and we should be able to quickly construct this
c
      integer i, j, k
c
      call dfill(nmo*nmo, 0.0d0, fa, 1)
c
      do j = 1, nmo
         do i = 1, nmo
            fa(i,j) = ha(i,j)
            fb(i,j) = hb(i,j)
            do k = 1, noa
               fa(i,j) = fa(i,j) + iaa(i,k,j,k)
               fb(i,j) = fb(i,j) + iba(i,k,j,k)
            end do
            do k = 1, nob
               fa(i,j) = fa(i,j) + iab(i,k,j,k)
               fb(i,j) = fb(i,j) + ibb(i,k,j,k)
            end do
         end do
      end do
c
      end
      subroutine jan_h(
     $     rtdb, basis,
     $     n1, n2,
     $     lda1, lda2,
     $     c1t, c2t,
     $     h)
      implicit none
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
     $     call errquit('jan_transform: basis ', basis)
      call schwarz_init(geom, basis)
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('jan_h: nbf',basis)
c
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'tmp1', l_tmp1, k_tmp1))
     $     call errquit('tmp1', nbf*nbf)
      if (.not. ma_push_get(mt_dbl, nbf*nbf,'tmp2', l_tmp2, k_tmp2))
     $     call errquit('tmp2', nbf*nbf)
c
      if (.not. ga_create(mt_dbl, nbf, nbf, 'tmp', 1, 1, g_tmp))
     &     call errquit('scf_v_g: tmp', 0)
      call ga_zero(g_tmp)
      call int_1e_ga(basis, basis, g_tmp, 'kinetic', .false.)
      call int_1e_ga(basis, basis, g_tmp, 'potential', .false.)
      call ga_get(g_tmp, 1, nbf, 1, nbf, dbl_mb(k_tmp1), nbf)
      if (.not. ga_destroy(g_tmp)) call errquit('ga?',0)
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
      if (.not. ma_chop_stack(l_tmp1)) call errquit('ma',0)
c
      return
      end

