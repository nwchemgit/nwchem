      subroutine uccsdt_eterm(qO,qV)
c
c$Id: eterm.F,v 1.6 2002-10-08 23:29:23 bert Exp $
c
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer qO(0:7,8), qV(0:7,8) 
      integer fO(2,0:7), fV(2,0:7)
      integer symnf,symme,symna,symf,syma,symn,symi
      integer alo_a,ahi_a,alo_b,ahi_b,q_ijlo,q_ijhi
      integer t_fhi,t_flo,t_ilo,t_ihi,t_jlo,t_jhi
      integer g_t2_mena_pure,g_t2_mena_mix,inum,n,a,l_n,k_n
      integer g_raa,g_rab,g_rbb,d_amp
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
      call make_fock_mo(fO,fV)
c
      do symf = 0, nir-1
         syma = symf
         call ga_add(-1.0d0,qO(symf,2),-0.5d0,qO(symf,1),qO(symf,2))
         call ga_add( 1.0d0,qO(symf,2), 1.0d0,fO(1,symf),qO(symf,2))
         call ga_add(-1.0d0,qO(symf,5),-0.5d0,qO(symf,4),qO(symf,5))
         call ga_add( 1.0d0,qO(symf,5), 1.0d0,fO(2,symf),qO(symf,5))
         call ga_add( 1.0d0,qV(symf,2), 0.5d0,qV(symf,1),qV(symf,2))
         call ga_add( 1.0d0,qV(symf,2), 1.0d0,fV(1,symf),qV(symf,2))
         call ga_add( 1.0d0,qV(symf,5), 0.5d0,qV(symf,4),qV(symf,5))
         call ga_add( 1.0d0,qV(symf,5), 1.0d0,fV(2,symf),qV(symf,5))
         call ga_add(-1.0d0,qO(symf,8),-0.5d0,qO(symf,4),qO(symf,8))
         call ga_add( 1.0d0,qO(symf,8), 1.0d0,fO(1,symf),qO(symf,8))
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
     &                     a,a,dbl_mb(k_n),no_sym(symn,1))
               do n = 1, no_sym(symn,1)
                  t_flo = ov_off(n,symnf,1,1)+1
                  t_fhi = t_flo + no_sym(symi,1)
                  call ga_add_patch(dbl_mb(k_n+n-1),g_t2_mena_pure,
     &                              t_ilo,t_ihi,t_jlo,t_jhi,1.0d0,
     &                              g_raa,t_ilo,t_ihi,t_flo,t_fhi,
     &                              g_raa,t_ilo,t_ihi,t_flo,t_fhi)
               enddo
c              matmul(g_t2_mena_pure,qV(2,symn)) into R(E2) aa   #(3)
               q_ijlo = o_sym(1,symi,1)
               q_ijhi = o_sym(2,symi,1)
               call ga_matmul_patch('n','n',-1.0d0,1.0d0,
     &                g_t2_mena_pure,t_ilo,t_ihi,t_jlo,t_jhi,
     &                qV(symn,2),q_ijlo,q_ijhi,q_ijlo,q_ijhi,
     &                g_raa,t_ilo,t_ihi,t_jlo,t_jhi)
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
     &                     a,a,dbl_mb(k_n),no_sym(symn,1))
               do n = 1, no_sym(symn,2)
                  t_flo = ov_off(n,symnf,2,2)+1
                  t_fhi = t_flo + no_sym(symi,2)
                  call ga_add_patch(dbl_mb(k_n+n-1),g_t2_mena_pure,
     &                              t_ilo,t_ihi,t_jlo,t_jhi,1.0d0,
     &                              g_rbb,t_ilo,t_ihi,t_flo,t_fhi,
     &                              g_rbb,t_ilo,t_ihi,t_flo,t_fhi)
               enddo
c              matmul(g_t2_mena_pure,qV(5,symn)) into R(E2) bb   #(4)
               q_ijlo = o_sym(1,symi,2)
               q_ijhi = o_sym(2,symi,2)
               call ga_matmul_patch('n','n',-1.0d0,1.0d0,
     &                g_t2_mena_pure,t_ilo,t_ihi,t_jlo,t_jhi,
     &                qV(symn,5),q_ijlo,q_ijhi,q_ijlo,q_ijhi,
     &                g_rbb,t_ilo,t_ihi,t_jlo,t_jhi)
            end do
            if (.not. ma_pop_stack(l_n)) call
     &          errquit('eterm: n dealloc failed',0)
         end do
         if (.not. ga_destroy(g_t2_mena_pure)) call 
     &       errquit('eterm: g_t2_mena_pure bb dealloc failed',0)
c
         call uccsdt_ampfile_read_t2(D_AMP,1,1,2,2,symna,alo_b,ahi_b,
     &                               g_t2_mena_mix,.true.,'block')
         do symf = 0, nir-1
            syma = symf
            symn = ieor(symf,symnf)
            symi = symn
            if (.not. ma_push_get(mt_dbl,no_sym(symn,1),'n array',
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0)
            do a = v_sym(1,syma,1), v_sym(2,syma,1)
               t_ilo = ov_off(a,symnf,1,1)+1
               t_ihi = t_jlo + no_sym(symi,1)
               t_jlo = 1
               t_jhi = ov_len(symnf,2,2)
c              matmul(g_t2_mena_mix,qO(syma,2)) into R(E1) ab   #(1)
               call ga_get(qO(symn,2),o_sym(1,symn,1),o_sym(2,symn,1),
     &                     a,a,dbl_mb(k_n),no_sym(symn,1))
               do n = 1, no_sym(symn,1) 
                  t_flo = ov_off(n,symnf,1,1)+1
                  t_fhi = t_flo + no_sym(symi,1)
                  call ga_add_patch(dbl_mb(k_n+n-1),g_t2_mena_mix,
     &                              t_ilo,t_ihi,t_jlo,t_jhi,1.0d0,
     &                              g_rab,t_flo,t_fhi,t_jlo,t_jhi,
     &                              g_rab,t_flo,t_fhi,t_jlo,t_jhi)
               enddo
c              -matmul(g_t2_mena_mix,qV(symn,2)) into R(E4) ab   #(3)
               q_ijlo = o_sym(1,symi,2)
               q_ijhi = o_sym(2,symi,2)
               call ga_matmul_patch('t','n',-1.0d0,1.0d0,
     &                qV(symn,2),q_ijlo,q_ijhi,q_ijlo,q_ijhi,
     &                g_t2_mena_mix,t_ilo,t_ihi,t_jlo,t_jhi,
     &                g_rab,t_ilo,t_ihi,t_jlo,t_jhi)
            end do
            if (.not. ma_pop_stack(l_n)) call
     &          errquit('eterm: n dealloc failed',0)
c
            if (.not. ma_push_get(mt_dbl,no_sym(symn,2),'n array',
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0)
            do a = v_sym(1,syma,2), v_sym(2,syma,2)
               t_ilo = 1
               t_ihi = ov_len(symnf,1,1)
               t_jlo = ov_off(a,symnf,2,2)+1
               t_jhi = t_jlo + no_sym(symi,2)
c              matmul(g_t2_mena_mix,qO(symn,8)) into R(E2) ab   #(5)
               call ga_get(qO(symn,8),o_sym(1,symn,2),o_sym(2,symn,2),
     &                     a,a,dbl_mb(k_n),no_sym(symn,2))
               do n = 1, no_sym(symn,2) 
                  t_flo = ov_off(n,symnf,2,2)+1
                  t_fhi = t_flo + no_sym(symi,2)
                  call ga_add_patch(dbl_mb(k_n+n-1),g_t2_mena_mix,
     &                              t_ilo,t_ihi,t_jlo,t_jhi,1.0d0,
     &                              g_rab,t_ilo,t_ihi,t_flo,t_fhi,
     &                              g_rab,t_ilo,t_ihi,t_flo,t_fhi)
               enddo
c              matmul(g_t2_mena_mix,qV(symn,5)) into R(E3) ab   #(4)
               q_ijlo = o_sym(1,symi,2)
               q_ijhi = o_sym(2,symi,2)
               call ga_matmul_patch('n','n',-1.0d0,1.0d0,
     &                g_t2_mena_mix,t_ilo,t_ihi,t_jlo,t_jhi,
     &                qV(symn,2),q_ijlo,q_ijhi,q_ijlo,q_ijhi,
     &                g_rab,t_ilo,t_ihi,t_jlo,t_jhi)
            end do
            if (.not. ma_pop_stack(l_n)) call
     &          errquit('eterm: n dealloc failed',0)
         end do
         if (.not. ga_destroy(g_t2_mena_mix)) call 
     &       errquit('eterm: g_t2_mena_mix ab dealloc failed',0)
      end do

      do symf = 0, nir-1
         do inum = 1, 2
            if (.not. ga_destroy(fO(inum,symf)))
     &          call errquit('eterm: fO dealloc failed',inum)
            if (.not. ga_destroy(fV(inum,symf)))
     &          call errquit('eterm: fV dealloc failed',inum)
         end do
         do inum = 1, 8
            if (qO(symf,inum) .gt. -1) then
                if (.not. ga_destroy(qO(inum,symf)))
     &             call errquit('eterm: qQ dealloc failed',inum)
            end if
            if (qV(symf,inum) .gt. -1) then
                if (.not. ga_destroy(qV(inum,symf)))
     &             call errquit('eterm: qV dealloc failed',inum)
            end if
         end do
      end do
c
      end
c
      subroutine make_fock_mo(fO_handles,fV_handles)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer fO_handles(2,0:7), fV_handles(2,0:7)
C
C     MAKE FOCK MATRIX ALPHA AND BETA IN AO BASIS AND TRANSFORM WITH
C     G_PART(LEFT SIDE) AND G_HOLE(RIGHT SIDE) TO FOCK MATRIX IN MO BASIS
C     ADD FINAL FOCK MATRIX PIECES TO CORRECT FO AND FV HANDLES
C     FO CONTAINS VIRTUALS AND FV CONTAINS OCCUPIED ORBITALS
C     ALL PER SYMMETRY AND SPIN-BLOCK
C
c
c     fij = hij + <ik||jk> = hij + <ik||jk>(same spin) + <ik|jk>(diff spin)
c
      return
c
      end

