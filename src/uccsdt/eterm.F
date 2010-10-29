      subroutine uccsdt_eterm(qO,qV)
c
c$Id$
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
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0,0)
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
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0,0)
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
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0,0)
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
     &          l_n,k_n)) call errquit('eterm: n alloc failed',0,0)
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
     &          call errquit('eterm: fO dealloc failed',inum,0)
            if (.not. ga_destroy(fV(inum,symf)))
     &          call errquit('eterm: fV dealloc failed',inum,0)
         end do
         do inum = 1, 8
            if (qO(symf,inum) .gt. -1) then
                if (.not. ga_destroy(qO(inum,symf)))
     &             call errquit('eterm: qQ dealloc failed',inum,0)
            end if
            if (qV(symf,inum) .gt. -1) then
                if (.not. ga_destroy(qV(inum,symf)))
     &             call errquit('eterm: qV dealloc failed',inum,0)
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
#include "geom.fh"
#include "global.fh"
#include "cuccsdtP.fh"
c
      integer fO_handles(2,0:7), fV_handles(2,0:7)
      integer d(4),f(4),nalpha,nbeta,spini,symf
      integer g_a_coul,g_a_exch,g_a_dens
      integer g_b_coul,g_b_exch,g_b_dens
      integer g_Fa,g_Fb,g_F_MO(2),dim_Fa,dim_Fb
      integer o_lo,o_hi,v_lo,v_hi,o_len,v_len
      double precision jfac(4),kfac(4),tol2e
      logical oskel
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      data jfac / 1.0d0, 0.0d0, 1.0d0, 0.0d0/
      data kfac / 0.0d0, 1.0d0, 0.0d0, 1.0d0/
c
c     Create Particle/Hole transformed MO Fock matrix:
c       1. Generate AO Fock matrix for alpha and beta
c       2. Transform to MO basis: g_part * F * g_hole
c       3. Place Fock matrix pieces into fO and fV handles
c          fO contains virtuals and fV contains occupied orbitals
c          all per symmetry and spin-block
c
c     Create Coulomb, Exchange and Density matrix for alpha and beta 
c
      g_a_coul = ga_create_atom_blocked(geom,basis,'eterm:a coul')
      g_b_coul = ga_create_atom_blocked(geom,basis,'eterm:b coul')
      g_a_exch = ga_create_atom_blocked(geom,basis,'eterm:a exch')
      g_b_exch = ga_create_atom_blocked(geom,basis,'eterm:b exch')
      g_a_dens = ga_create_atom_blocked(geom,basis,'eterm:a dens')
      g_b_dens = ga_create_atom_blocked(geom,basis,'eterm:b dens')
c
c     Create Density matrices from vectors in g_movecs(2)
c
      oskel = .false.
      tol2e = 1.0d-15
c
      nalpha = nc(1) + no(1)
      nbeta = nc(2) + no(2)
      call ga_zero(g_a_dens)
      call ga_zero(g_b_dens)
      call ga_dgemm('n','t',nbf,nbf,nalpha,1.0d0,g_movecs(1),
     &     g_movecs(1),0.0d0,g_a_dens)
      call ga_dgemm('n','t',nbf,nbf,nbeta,1.0d0,g_movecs(2),
     &     g_movecs(2),0.0d0,g_b_dens)
      if (oskel) then
         call sym_symmetrize(geom, basis, .true., g_a_dens)
         call sym_symmetrize(geom, basis, .true., g_b_dens)
      endif
c
c     Get Coulomb and Exchange matrices
c
      d(1) = g_a_dens
      d(2) = g_a_dens
      d(3) = g_b_dens
      d(4) = g_b_dens
      f(1) = g_a_coul
      f(2) = g_a_exch
      f(3) = g_b_coul
      f(4) = g_b_exch
      call fock_2e(geom, basis, 4, jfac, kfac, tol2e,
     &     oskel, d, f, .false.)
c
c     Assemble alpha and beta Fock matrices:
c     Fa = h + Ja + Jb - Ka
c     Fb = h + Ja + Jb - Kb
c
      call ga_zero(g_Fa)
      call int_1e_ga(basis, basis, g_Fa, 'kinetic', oskel)
      call int_1e_ga(basis, basis, g_Fa, 'potential', oskel)
      if (oskel) then
         call sym_symmetrize(geom, basis, .false., g_Fa)
      endif
      call ga_dadd(1.0d0,g_a_coul,1.0d0,g_Fa,g_Fa)
      call ga_dadd(1.0d0,g_b_coul,1.0d0,g_Fa,g_Fa)
      call ga_copy(g_Fa, g_Fb)
      call ga_dadd(1.0d0,g_Fa,-1.0d0,g_a_exch,g_Fa)
      call ga_dadd(1.0d0,g_Fb,-1.0d0,g_b_exch,g_Fb)
c
c     Clean up dead GA's
c
      if (.not. ga_destroy(g_a_dens)) call errquit('eterm: destroy',1,0)
      if (.not. ga_destroy(g_b_dens)) call errquit('eterm: destroy',2,0)
      if (.not. ga_destroy(g_a_exch)) call errquit('eterm: destroy',3,0)
      if (.not. ga_destroy(g_b_exch)) call errquit('eterm: destroy',4,0)
      if (.not. ga_destroy(g_a_coul)) call errquit('eterm: destroy',5,0)
      if (.not. ga_destroy(g_b_coul)) call errquit('eterm: destroy',6,0)
c
c     Transform AO Fock matrices to MO Fock matrices 
c
      dim_Fa = no(1) + nv(1)
      dim_Fb = no(2) + nv(2)
      if (.not. ga_create(MT_DBL,dim_Fa,dim_Fa,'g_F_MO 1',0,0,
     &    g_F_MO(1))) call errquit('eterm: g_F_MO 1 alloc failed',0,0)
      if (.not. ga_create(MT_DBL,dim_Fb,dim_Fb,'g_F_MO 2',0,0,
     &    g_F_MO(2))) call errquit('eterm: g_F_MO 2 alloc failed',0,0)
      call two_index_transf(g_Fa,g_part(1),g_hole(1),g_a_dens,g_F_MO(1))
      call two_index_transf(g_Fb,g_part(2),g_hole(2),g_a_dens,g_F_MO(2))
      if (.not. ga_destroy(g_Fa)) call errquit('eterm: destroy Fa',0,0)
      if (.not. ga_destroy(g_Fb)) call errquit('eterm: destroy Fb',0,0)
c
c     Store the pieces in the right places in fO and fV
c
      do spini = 1, 2
         do symf = 0, nir-1
            o_lo = o_sym(1,symf,spini)
            o_hi = o_sym(2,symf,spini)
            v_lo = v_sym(1,symf,spini)
            v_hi = v_sym(2,symf,spini)
            o_len = o_hi - o_lo + 1
            v_len = v_hi - v_lo + 1
            if (.not. ga_create(MT_DBL,o_len,o_len,'fV',0,0,
     &          fV_handles(spini,symf))) call 
     &          errquit('eterm: fV alloc failed',symf)
            call ga_copy_patch('n',g_F_MO(spini),o_lo,o_hi,o_lo,o_hi,
     &                         fV_handles(spini,symf),1,o_len,1,o_len)
            if (.not. ga_create(MT_DBL,v_len,v_len,'fO',0,0,
     &          fO_handles(spini,symf))) call 
     &          errquit('eterm: fO alloc failed',symf)
            call ga_copy_patch('n',g_F_MO(spini),v_lo,v_hi,v_lo,v_hi,
     &                         fO_handles(spini,symf),1,v_len,1,v_len)
         end do
      end do
c
c     Do some final cleanup
c
      if (.not. ga_destroy(g_F_MO(1))) 
     &              call errquit('eterm: Mdestroy',1,0)
      if (.not. ga_destroy(g_F_MO(2)))
     &              call errquit('eterm: Mdestroy',2,0)
c
      return
c
      end

