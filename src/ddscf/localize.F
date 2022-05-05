      subroutine fudgeit(rtdb)
      implicit none
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
c
      integer rtdb
      integer g_over, g_sc
      integer ga_create_atom_blocked
c
      g_over = ga_create_atom_blocked(geom, basis, 'tmp')
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
      call ga_matpow(g_over, -0.5d0, 1d-12)
      call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $     1.0d0, g_over, g_movecs, 0.0d0, g_sc)
      

      end
c
c     =================================================================
c      
      subroutine testlocalize(rtdb)      
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
      integer rtdb
      integer maxnloc
      parameter (maxnloc = 10000)
      integer iloc(maxnloc)
      integer g_over, i, g_sc, k_sc, l_sc, k_c, l_c, nloc, ncore
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      character*8 loctype
      double precision y, w, ef, scale
      logical debug
c
      integer g_uc(4), x
*      integer  g_mosc

      debug = .false. .and. ga_nodeid().eq.0
c
      if (.not. rtdb_cget(rtdb, 'scf:loctype', 1, loctype)) 
     $     loctype = 'PM'
      if (.not. rtdb_get(rtdb, 'scf:y', mt_dbl, 1, y))
     $     y = 0.0d0
      if (.not. rtdb_get(rtdb, 'scf:w', mt_dbl, 1, w))
     $     w = 1.0d0

      if (debug) then
        write(luout,*) 'hello from testlocalize'
        write(luout,*) 'y, w, =',y,w
      end if
c
      if (y.ne.0.0d0) then
c
c     Scale the MOs by the laplace quadrature factors
c
         ef = 0.5d0*(dbl_mb(k_eval+nclosed-1) + dbl_mb(k_eval+nclosed))
         if (.not. ma_push_get(mt_dbl, nbf, 'c', l_c, k_c))
     $        call errquit('ma for c', 0, MA_ERR)
         call ga_sync
         if (ga_nodeid() .eq. 0)
     $        write(6,*) ' EF ', ef, ' Y ', y, ' W ', w
         do i = ga_nodeid()+1, nclosed, ga_nnodes()
            scale = (w**0.125d0)*exp((dbl_mb(k_eval-1+i)-ef)*0.5d0*y)
*     write(6,*) i, dbl_mb(k_eval-1+i),  scale
            call ga_get(g_movecs,1,nbf,i,i,dbl_mb(k_c), 1)
            call dscal(nbf,scale,dbl_mb(k_c),1)
            call ga_put(g_movecs,1,nbf,i,i,dbl_mb(k_c), 1)
         enddo
         do i = ga_nodeid()+nclosed+1, nmo, ga_nnodes()
            scale = (w**0.125d0)*exp(-(dbl_mb(k_eval-1+i)-ef)*0.5d0*y)
*     write(6,*) i, dbl_mb(k_eval-1+i),  scale
            call ga_get(g_movecs,1,nbf,i,i,dbl_mb(k_c), 1)
            call dscal(nbf,scale,dbl_mb(k_c),1)
            call ga_put(g_movecs,1,nbf,i,i,dbl_mb(k_c), 1)
         enddo
         if (.not. ma_pop_stack(l_c)) call errquit('pop',0, MA_ERR)
      endif
      call ga_sync
c
      if (loctype .eq. 'NONE') then
*         call ga_print(g_movecs)
         return
      endif
c
      if (loctype .eq. 'PM') then
         if (ga_nodeid().eq.0)
     $        call util_print_centered(6,
     $        'Pipek-Mezey orbital localization', 40, .true.)
         g_over  = ga_create_atom_blocked(geom, basis, 'rohf_mull:over')
         call ga_zero(g_over)
         call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
c     
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'sc',
     $        nbf, 0, g_sc)) call errquit('testlocalize: sc',0, GA_ERR)
         call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $        1.0d0, g_over, g_movecs, 0.0d0, g_sc)
c     
         if (.not. ma_push_get(mt_dbl, 2*nbf, 'sc', l_sc, k_sc))
     $        call errquit('ma for sc', 0, MA_ERR)
         if (.not. ma_push_get(mt_dbl, 2*nbf, 'c', l_c, k_c))
     $        call errquit('ma for c', 0, MA_ERR)
c     
c     Localize the core orbitals
c     
         if (.not. geom_num_core(rtdb, geom, 'ddscf', ncore)) ncore = 0
         if (ncore .gt. 0) then
            do i = 1, ncore
               iloc(i) = i
            end do
            nloc = ncore
            call localizePM(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $           nloc, iloc, nbf, nmo, g_movecs, g_sc)
         end if
c     
c     Localized the occupied
c     
         do i = ncore+1, nclosed
            iloc(i-ncore) = i
         end do
         nloc = nclosed - ncore
         call localizePM(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $        nloc, iloc, nbf, nmo, g_movecs, g_sc)
c     
c     Do the virtuals
c
         do i = nclosed+1, nmo
            iloc(i-nclosed) = i
         end do
         nloc = nmo - nclosed
         call localizePM(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $        nloc, iloc, nbf, nmo, g_movecs, g_sc)
c     
*     *      call ga_print(g_movecs)
         if (.not. ga_destroy(g_over)) call errquit('ga?',0, GA_ERR)
         if (.not. ga_destroy(g_sc)) call errquit('ga?',0, GA_ERR)
         if (.not. ma_pop_stack(l_c)) call errquit('c',0, MA_ERR)
         if (.not. ma_pop_stack(l_sc)) call errquit('sc',0, MA_ERR)
c
      else                      ! Foster-Boys localization
         if (ga_nodeid().eq.0) 
     $        call util_print_centered(6,
     $        'Foster-Boys orbital localization', 40, .true.)
c
         do x= 1, 4
*ga:1:0
            if (.not. ga_create(MT_DBL, nbf, nbf, 'uc',
     $           nbf, 0, g_uc(x))) call errquit('testlocalize: uc',x,
     &       GA_ERR)
         end do
         call int_dip_ga(basis, basis, g_uc(1), g_uc(2), g_uc(3))
*ga:1:0
         if (.not. ga_create(MT_DBL, nbf, nmo, 'sc',
     $        nbf, 0, g_sc)) call errquit('testlocalize: sc',0, GA_ERR)
         do x = 1, 3
            call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $           1.0d0, g_uc(x), g_movecs, 0.0d0, g_sc)
            call ga_copy_patch('n',g_sc,1,nbf,1,nmo,g_uc(x),1,nbf,1,nmo)
         end do
         g_over  = ga_create_atom_blocked(geom, basis, 'rohf_mull:over')
         call ga_zero(g_over)
         call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
c     
         call ga_dgemm('n', 'n', nbf, nmo, nbf, 
     $        1.0d0, g_over, g_movecs, 0.0d0, g_uc(4))
c     
         if (.not. ma_push_get(mt_dbl, 8*nbf, 'sc', l_sc, k_sc))
     $        call errquit('ma for sc', 0, MA_ERR)
         if (.not. ma_push_get(mt_dbl, 8*nbf, 'c', l_c, k_c))
     $        call errquit('ma for c', 0, MA_ERR)
c     
c     Localize the core orbitals
c     
         if (.not. geom_num_core(rtdb,geom,'ddscf',ncore))
     $       ncore = 0
         if (ncore .gt. 0) then
            do i = 1, ncore
               iloc(i) = i
            end do
            nloc = ncore
            call localizeFB(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $           nloc, iloc, nbf, nmo, g_movecs, g_uc)
         end if
c     
c     Localized the occupied
c     
         do i = ncore+1, nclosed
            iloc(i-ncore) = i
         end do
         nloc = nclosed - ncore
         call localizeFB(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $        nloc, iloc, nbf, nmo, g_movecs, g_uc)
c     
c     Do the virtuals
c
         do i = nclosed+1, nmo
            iloc(i-nclosed) = i
         end do
         nloc = nmo - nclosed
         call localizeFB(basis, dbl_mb(k_c), dbl_mb(k_sc), 
     $        nloc, iloc, nbf, nmo, g_movecs, g_uc)
c     
*         call ga_print(g_movecs)
         do x = 1, 4
            if (.not. ga_destroy(g_uc(x))) call errquit('ga?',x, GA_ERR)
         end do
         if (.not. ga_destroy(g_over)) call errquit('ga?',0, GA_ERR)
         if (.not. ga_destroy(g_sc)) call errquit('ga?',0, GA_ERR)
         if (.not. ma_pop_stack(l_c)) call errquit('c',0, MA_ERR)
         if (.not. ma_pop_stack(l_sc)) call errquit('sc',0, MA_ERR)
c
      end if
c
      call ga_sync()
c
*      call moints_screen(basis, ncore+1, nclosed, g_movecs(1), g_mosc)
*      if (.not. ga_destroy(g_mosc)) call errquit('fred did it!',0)
c
      end
c
c     =================================================================
c      
      subroutine localizePM(basis, c, sc, nloc, iloc, nbf, nmo,
     $  g_c, g_sc)
c
c     Pipek-Mezey localization
c     
      implicit none
      
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "util.fh"
c
c     Localize the nloc orbitals in iloc(*) by mixing with each other
c
      integer basis, nloc, iloc(*), nbf, nmo
      double precision c(nbf, 2), sc(nbf, 2)
      integer g_c, g_sc
      integer maxat, nlist
*............................   these should be dynamically allocated ?
      parameter (maxat = nw_max_atom)
      integer list(maxat)
      double precision pop(maxat)
c
      integer iter, ss, s, tt, t, a, u, bflo, bfhi, natoms, geom
      double precision ast, bst, qast, qat, qas, gamma, cosg, sing, d,
     $     qs, dprev, tol, dmax, gamma_tol, gamma_max, tmp
      integer nrot, set, pair, neven
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('localize: basis ', 0, BASIS_ERR)
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('localize: geom',0, GEOM_ERR)
c
      if (natoms.gt.maxat) call errquit
     &      ('localize: maxat too small ', 911, UNKNOWN_ERR)
c
      tol = 1d-8
      gamma_tol = 1d-10
c
      if (ga_nodeid() .eq. 0) then
         write(6,2)
 2       format(/10x,' iter   Max. delocal   Mean delocal    Converge'/
     $        10x,' ----   ------------   ------------   ---------')
         call util_flush(6)
      end if
c
      dprev = 0.0d0
      gamma_max = 0.0d0
      do iter = 1, 100
         call ga_sync
         nrot = 0
c
c     Analyze convergence by forming P
c
         d = 0.0d0
         dmax = 0.0d0
         do ss = 1+ga_nodeid(), nloc, ga_nnodes()
            s = iloc(ss)
            call ga_get(g_c,  1, nbf, s, s, c(1,1), 1)
            call ga_get(g_sc, 1, nbf, s, s,sc(1,1), 1)
            qs = 0.0d0
            do a = 1, natoms
               if (.not. bas_ce2bfr(basis, a, bflo, bfhi))
     $              call errquit('localized: basis ', 0, BASIS_ERR)
               qas  = 0.0d0
               do u = bflo, bfhi
                  qas  = qas  + c(u,1)*sc(u,1)
               end do
               qs = qs + qas**2
            end do
*            write(6,*) ' ds ', s, 1.0d0/qs
            dmax = max(dmax, 1.0d0/qs)
            d = d + 1.0d0/qs
         end do
c
#ifdef NWCHEM_USE_GOP_ABSMAX
         call ga_dgop(1, gamma_max, 1, 'absmax')
         call ga_dgop(1, dmax, 1, 'absmax')
#else
         gamma_max = abs(gamma_max)
         dmax      = abs(dmax)
         call ga_dgop(1, gamma_max, 1, 'max')
         call ga_dgop(1, dmax, 1, 'max')
#endif
         call ga_dgop(2, d , 1, '+')
c
         if (ga_nodeid() .eq. 0) then
            write(6,1) iter, dmax, d/dble(nloc), gamma_max
 1          format(10x, i5, 2f15.10, 1p,d12.2,d12.2)
            call util_flush(6)
         end if
         call ga_sync
c
**         if (abs(d-dprev)/dble(nloc) .lt. tol) goto 1000
**         dprev = d
         if (iter.gt.1 .and. gamma_max.lt.tol) goto 1000
         gamma_max = 0.0d0
c
c     Loop over pairs with as much parallelism as possible
c
         neven = nloc + mod(nloc,2)
         do set = 1, neven-1
            do pair = 1+ga_nodeid(), neven/2, ga_nnodes()
               call localize_pairs(neven, set, pair, ss, tt)
               if (tt .le. nloc) then
                  s = iloc(ss)
                  t = iloc(tt)
*                  write(6,*) nloc, neven, set, pair, ss, tt, s, t
                  call ga_get(g_c,  1, nbf, s, s, c(1,1), 1)
                  call ga_get(g_sc, 1, nbf, s, s,sc(1,1), 1)
                  call ga_get(g_c,  1, nbf, t, t, c(1,2), 1)
                  call ga_get(g_sc, 1, nbf, t, t,sc(1,2), 1)
c     
c     Form rotation information
c     
                  ast = 0.0d0
                  bst = 0.0d0
                  do a = 1, natoms
                     if (.not. bas_ce2bfr(basis, a, bflo, bfhi))
     $                    call errquit('localized: basis ', 0,
     &       BASIS_ERR)
c     
                     qast = 0.0d0
                     qas  = 0.0d0
                     qat  = 0.0d0
                     do u = bflo, bfhi
                        qast = qast + c(u,2)*sc(u,1) + c(u,1)*sc(u,2)
                        qas  = qas  + c(u,1)*sc(u,1)
                        qat  = qat  + c(u,2)*sc(u,2)
                     end do
                     qast = qast * 0.5d0
c     
                     ast = ast + qast**2 - 0.25d0*(qas - qat)**2
                     bst = bst + qast*(qas - qat)
                  end do
c     
                  gamma = 0.25d0*acos(-ast/sqrt(ast**2+bst**2))
                  gamma = sign(gamma,bst)
                  gamma_max = max(gamma_max, abs(gamma))
*                  if (iter .eq. 1 .and. abs(gamma).lt.0.01d0) then
*                     gamma = (util_random(0)-0.5d0)*3.14d0
*                  endif
c     
                  if (abs(gamma) .gt. gamma_tol) then
                     nrot = nrot + 1
                     cosg = cos(gamma)
                     sing = sin(gamma)
c     
c     Do the rotation of C and SC
c     
                     call drot(nbf, c(1,1), 1, c(1,2), 1, cosg, sing)
                     call drot(nbf,sc(1,1), 1,sc(1,2), 1, cosg, sing)
                     call ga_put(g_c,  1, nbf, s, s, c(1,1), 1)
                     call ga_put(g_sc, 1, nbf, s, s,sc(1,1), 1)
                     call ga_put(g_c,  1, nbf, t, t, c(1,2), 1)
                     call ga_put(g_sc, 1, nbf, t, t,sc(1,2), 1)
                  end if
               end if
            end do
            call ga_sync
         end do
      end do
c
 1000 continue

      
c     Analyze localization of each mo:
c     per lmo, a list of atomic populations is printed
c     in decreasing magnitude, with the polulations in parentheses. 
c
      if (ga_nodeid() .eq. 0) then
         write(6,*)
         do ss = 1, nloc
            s = iloc(ss)
            call ga_get(g_c,  1, nbf, s, s, c(1,1), 1)
            call ga_get(g_sc, 1, nbf, s, s,sc(1,1), 1)
            nlist = 0
            do a = 1, natoms
               if (.not. bas_ce2bfr(basis, a, bflo, bfhi))
     $              call errquit('localized: basis ', 0,
     &       BASIS_ERR)
               qas  = 0.0d0
               do u = bflo, bfhi
                  qas  = qas  + c(u,1)*sc(u,1)
               end do
               if (abs(qas) .gt. 0.01d0) then
                  nlist = nlist + 1
                  list(nlist) = a
                  pop(nlist) = qas
               end if
            end do
            do u = 1, nlist
               do t = 1, u-1
                  if (abs(pop(t)).lt.abs(pop(u))) then
                     tmp = pop(u)
                     pop(u) = pop(t)
                     pop(t) = tmp
                     tt = list(u)
                     list(u) = list(t)
                     list(t) = tt
                  end if
               end do
            end do
            write(6,77) s, (list(a), pop(a), a=1,nlist)
 77         format(i5, 100(2x,i4,'(',f5.2,')'))
         end do
         call util_flush(6)
      end if
c
      call ga_sync
c
      end      
c
c     =================================================================
c      
      subroutine localizeFB(basis, c, uc, nloc, iloc, nbf, nmo,
     $  g_c, g_uc)
c
c     Foster-Boys localization
c     
      
      implicit none
      
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "util.fh"
c
c     Localize the nloc orbitals in iloc(*) by mixing with each other
c
      integer basis, nloc, iloc(*), nbf, nmo
      double precision c(nbf, 2), uc(nbf, 2, 4)
      integer g_c, g_uc(4) ! x, y, z, overlap
c
      integer iter, ss, s, tt, t, u, geom, a, bflo, bfhi
      double precision ast, bst, gamma, cosg, sing, d, tmp,
     $     qs, dprev, tol, dmax, gamma_tol, gamma_max, u1, u2, u12
      double precision u21
      integer nrot, set, pair, neven, x, natoms
      integer maxat, nlist
*............................   these should be dynamically allocated ?
      parameter (maxat = nw_max_atom)
      integer list(maxat)
      double precision pop(maxat), qas
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('localize: basis ', 0, BASIS_ERR)
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('localize: geom',0, GEOM_ERR)
c
      if (natoms.gt.maxat) call errquit
     &      ('localize: maxat too small ', 911,
     &       UNKNOWN_ERR)
c
      tol = 1d-8
      gamma_tol = 1d-10
      u21 = util_random(12345)
c
      if (ga_nodeid() .eq. 0) then
         write(6,2)
 2       format(/10x,' iter   Max. dipole2   Mean dipole2    Converge'/
     $        10x,' ----   ------------   ------------   ---------')
         call util_flush(6)
      end if
c
      dprev = 0.0d0
      gamma_max = 0.0d0
      do iter = 1, 100
         call ga_sync
         nrot = 0
c
c     Analyze convergence by forming functional
c
         d = 0.0d0
         dmax = 0.0d0
         do ss = 1+ga_nodeid(), nloc, ga_nnodes()
            s = iloc(ss)
            call ga_get(g_c,  1, nbf, s, s, c(1,1), 1)
            do x = 1, 3
               call ga_get(g_uc(x), 1, nbf, s, s,uc(1,1,x), 1)
            end do
            qs = 0.0d0
            do x = 1, 3
               u1 = 0.0d0
               do u = 1, nbf
                  u1 = u1 + c(u,1)*uc(u,1,x)
               end do
               qs  = qs  + u1**2
            end do
            dmax = max(dmax, qs)
            d = d + qs
         end do
c
#ifdef NWCHEM_USE_GOP_ABSMAX
         call ga_dgop(1, gamma_max, 1, 'absmax')
         call ga_dgop(1, dmax, 1, 'absmax')
#else
         gamma_max = abs(gamma_max)
         dmax      = abs(dmax)
         call ga_dgop(1, gamma_max, 1, 'max')
         call ga_dgop(1, dmax, 1, 'max')
#endif
         call ga_dgop(2, d , 1, '+')
c
         if (ga_nodeid() .eq. 0) then
            write(6,1) iter, dmax, d/dble(nloc), gamma_max
 1          format(10x, i5, 2f17.8, 1p,2d12.2)
            call util_flush(6)
         end if
         call ga_sync
c
*         if (abs(d-dprev)/dble(nloc) .lt. tol
*     $        .and. iter.gt.1) goto 1000
*         dprev = d
         if (iter.gt.1 .and. gamma_max.lt.tol) goto 1000
         gamma_max = 0.0d0
c
c     Loop over pairs with as much parallelism as possible
c
         neven = nloc + mod(nloc,2)
         do set = 1, neven-1
            do pair = 1+ga_nodeid(), neven/2, ga_nnodes()
               call localize_pairs(neven, set, pair, ss, tt)
               if (tt .le. nloc) then
                  s = iloc(ss)
                  t = iloc(tt)
*                  write(6,*) nloc, neven, set, pair, ss, tt, s, t
                  call ga_get(g_c,  1, nbf, s, s, c(1,1), 1)
                  call ga_get(g_c,  1, nbf, t, t, c(1,2), 1)
                  do x = 1, 4
                     call ga_get(g_uc(x), 1, nbf, s, s,uc(1,1,x), 1)
                     call ga_get(g_uc(x), 1, nbf, t, t,uc(1,2,x), 1)
                  end do
c     
c     Form rotation information
c     
                  ast = 0.0d0
                  bst = 0.0d0
                  do x = 1, 3
                     u1  = 0.0d0
                     u2  = 0.0d0
                     u12 = 0.0d0
                     u21 = 0.0d0
                     do u = 1, nbf
                        u12 = u12 + c(u,1)*uc(u,2,x)
                        u21 = u21 + c(u,2)*uc(u,1,x)
                        u1  = u1 + c(u,1)*uc(u,1,x)
                        u2  = u2 + c(u,2)*uc(u,2,x)
                     end do
                     ast = ast + u12*u12 - 0.25d0*(u1-u2)**2
                     bst = bst + u12*(u1 - u2)
                     if (abs(u12-u21)/max(1.0d0,abs(u12)).gt.1d-8) then
                        write(6,*) ' U12, U21 ', u12, u21
                        call errquit('bad u12', 0, UNKNOWN_ERR)
                     endif
                  end do
c     
                  gamma = 0.25d0*acos(-ast/sqrt(ast**2+bst**2))
                  gamma = sign(gamma,bst)
*                  gamma = 0.25*abs(atan2(bst,-ast))
                  gamma_max = max(gamma_max, abs(gamma))
                  if (iter .eq. 1 .and. abs(gamma).lt.0.01d0) then
                     gamma = (util_random(0)-0.5d0)*3.14d0
                  endif
c     
                  if (abs(gamma) .gt. gamma_tol) then
                     nrot = nrot + 1
                     cosg = cos(gamma)
                     sing = sin(gamma)
c     
c     Do the rotation of C and UC
c     
                     call drot(nbf, c(1,1), 1, c(1,2), 1, cosg, sing)
                     call ga_put(g_c,  1, nbf, s, s, c(1,1), 1)
                     call ga_put(g_c,  1, nbf, t, t, c(1,2), 1)
                     do x = 1, 4
                        call drot(nbf,uc(1,1,x), 1,uc(1,2,x), 1, 
     $                       cosg, sing)
                        call ga_put(g_uc(x), 1, nbf, s, s,uc(1,1,x), 1)
                        call ga_put(g_uc(x), 1, nbf, t, t,uc(1,2,x), 1)
                     end do
                  end if
               end if
            end do
            call ga_sync
         end do
      end do
c
 1000 continue
c
c     Analyze localization of each mo
c
      if (ga_nodeid() .eq. 0) then
         write(6,*)
         do ss = 1, nloc
            s = iloc(ss)
            call ga_get(g_c,  1, nbf, s, s, c(1,1), 1)
            call ga_get(g_uc(4), 1, nbf, s, s,uc(1,1,1), 1)
            nlist = 0
            do a = 1, natoms
               if (.not. bas_ce2bfr(basis, a, bflo, bfhi))
     $              call errquit('localized: basis ', 0, BASIS_ERR)
               qas  = 0.0d0
               do u = bflo, bfhi
                  qas  = qas  + c(u,1)*uc(u,1,1)
               end do
               if (abs(qas) .gt. 0.01d0) then
                  nlist = nlist + 1
                  list(nlist) = a
                  pop(nlist) = qas
               end if
            end do
            do u = 1, nlist
               do t = 1, u-1
                  if (abs(pop(t)).lt.abs(pop(u))) then
                     tmp = pop(u)
                     pop(u) = pop(t)
                     pop(t) = tmp
                     tt = list(u)
                     list(u) = list(t)
                     list(t) = tt
                  end if
               end do
            end do
            write(6,77) s, (list(a), pop(a), a=1,nlist)
 77         format(i5, 100(2x,i4,'(',f5.2,')'))
         end do
         call util_flush(6)
      end if
c
      call ga_sync
c
      end
c
c     =================================================================
c      
      subroutine localizeIBO(basis, minbas,
     &  c, nloc, iloc, nbf, nmo,
     &  g_c, g_smat, g_iao, mnbf, eval, occ)
c
c     IAO-IBO localization
c     
      
      implicit none
      
#include "errquit.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
c
c     Localize the nloc orbitals in iloc(*) by mixing with each other
c
      integer basis, minbas, nloc, iloc(*), nbf, nmo, mnbf
      double precision c(mnbf, 2), eval(nbf), occ(nbf)
      integer g_c, g_smat, g_iao
      integer maxat, nlist
*............................   these should be dynamically allocated ?
      parameter (maxat = nw_max_atom)
      integer list(maxat)
      double precision pop(maxat)

      integer iter, ss, s, tt, t, a, u, bflo, bfhi, natoms, geom
      double precision ast, bst, qast, qat, qas, gamma, cosg, sing, d,
     &  qs, dprev, tol, dmax, gamma_tol, gamma_max, tmp
      double precision minval, swap
      integer nrot, set, pair, neven
      integer ntmp1, ntmp2
      logical master, debug
      integer imo, jmo, kmo

      integer g_tmp1, g_tmp2, g_cibocc

      character*(11) pname

      character*(256) lmotrans
      logical file_write_ga
      external file_write_ga      

c     -----------------------------------------------------------------

      pname = 'iao-ibo loc'
      master = (ga_nodeid() .eq. 0)
      debug = .false.  ! true, during development

      if (debug.and.master) write(luout,*) 'entering '//pname
c
      if (.not. bas_geom(basis, geom)) call errquit
     $     (pname//': basis ', 0, BASIS_ERR)
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     (pname//': geom',0, GEOM_ERR)
c
      if (natoms.gt.maxat) call errquit
     &  (pname//': maxat too small ', 911, UNKNOWN_ERR)

c     note sure if we need this for sure, but later we assume this
      if (nloc.gt.mnbf) call errquit
     &  (pname//': nloc > mnbf', 66, UNKNOWN_ERR)
      if (nloc.gt.nbf) call errquit
     &  (pname//': nloc > nbf', 66, UNKNOWN_ERR)      


c     DEBUG some things related to the basis function centers
      if (debug) then
        do a = 1, natoms
          if (.not. bas_ce2bfr(basis, a, bflo, bfhi))
     $      call errquit('localized: basis ', 0, BASIS_ERR)
          if (master) then
            write(luout,'(1x,a,3i4)') 'ao basis, atom = ', a, bflo, bfhi
          end if
          if (.not. bas_ce2bfr(minbas, a, bflo, bfhi))
     &      call errquit('localized: basis ', 0, BASIS_ERR)
          if (master) then
            write(luout,'(1x,a,3i4)') 'minbas  , atom = ', a, bflo, bfhi
          end if
        end do
      end if ! debug

c     transform the set of occupied MOs to be localized to the IAO basis
c     via C(iao, trans) S C(MO)

      if (.not. ga_create(MT_DBL, mnbf, nloc , 'loc:cibocc',
     &  mnbf, 0, g_cibocc)) call errquit(pname//': cibocc',0, GA_ERR)

      if (.not. ga_create(MT_DBL, nbf, nloc , 'loc:tmp1',
     &  nbf, 0, g_tmp1)) call errquit(pname//': tmp1',0, GA_ERR)

      call ga_dgemm('n', 'n', nbf, nloc, nbf, 
     &  1.0d0, g_smat, g_c, 0.0d0, g_tmp1)

      call ga_dgemm('t', 'n', mnbf, nloc, nbf, 
     &  1.0d0, g_iao, g_tmp1, 0.0d0, g_cibocc)

      if (.not. ga_destroy(g_tmp1))
     &  call errquit(pname//': ga_destroy failed g_tmp1',0, GA_ERR)

c      if (debug) call ga_print(g_cibocc)

      if (debug .and. master) write(luout,*) 'cibocc done'

c     save a copy of the starting MOs for later

      if (.not.ga_duplicate(g_cibocc, g_tmp2, 'g_tmp2')) call
     &  errquit(pname//': ga_dup cibocc', 0, GA_ERR)

      call ga_copy(g_cibocc, g_tmp2)

c     the MOs in the orthonormal IAO basis must form an orthogonal
c     matrix. test that before localization, and again later

      if (.not. ga_create(MT_DBL, nloc, nloc , 'loc:tmp1',
     &  mnbf, 0, g_tmp1)) call errquit(pname//': tmp1',0, GA_ERR)

      call ga_dgemm('t', 'n', nloc, nloc, mnbf, 
     &  1.0d0, g_cibocc, g_cibocc, 0.0d0, g_tmp1)

      call ga_add_todiag(g_tmp1, -1.0d0)

      call ga_maxelt(g_tmp1, tmp)
      if (master) then
        write(luout,'(/1x,a,1x,f18.12)')
     &    'IBO loc: largest element of C(MO, trans) C(MO) -1: ',
     &    tmp
        write(luout,'(1x,a)')
     &    'should be zero, for CMOs in the IAO basis'
      end if

      if (.not. ga_destroy(g_tmp1))
     &  call errquit(pname//': ga_destroy failed g_tmp1',0, GA_ERR)

      
c
      tol = 1d-8
      gamma_tol = 1d-10
c
      if (master) then
         write(luout,2)
 2       format(/10x,' iter   Max. delocal   Mean delocal    Converge'/
     $        10x,' ----   ------------   ------------   ---------')
         call util_flush(luout)
      end if
c
      dprev = 0.0d0
      gamma_max = 0.0d0
      do iter = 1, 100
         call ga_sync
         nrot = 0
c
c     Analyze convergence by forming P
c
         d = 0.0d0
         dmax = 0.0d0
         do ss = 1+ga_nodeid(), nloc, ga_nnodes()
            s = iloc(ss)
            call ga_get(g_cibocc, 1, mnbf, s, s, c(1,1), 1)
            qs = 0.0d0
            do a = 1, natoms
               if (.not. bas_ce2bfr(minbas, a, bflo, bfhi))
     &              call errquit(pname//': basis ', 1, BASIS_ERR)
               qas  = 0.0d0
               do u = bflo, bfhi
                  qas  = qas  + c(u,1)*c(u,1)
               end do
               qs = qs + qas**2
            end do
*            write(luout,*) ' ds ', s, 1.0d0/qs
            dmax = max(dmax, 1.0d0/qs)
            d = d + 1.0d0/qs
         end do
c
#ifdef NWCHEM_USE_GOP_ABSMAX
         call ga_dgop(1, gamma_max, 1, 'absmax')
         call ga_dgop(1, dmax, 1, 'absmax')
#else
         gamma_max = abs(gamma_max)
         dmax      = abs(dmax)
         call ga_dgop(1, gamma_max, 1, 'max')
         call ga_dgop(1, dmax, 1, 'max')
#endif
         call ga_dgop(2, d , 1, '+')
c
         if (ga_nodeid() .eq. 0) then
            write(luout,1) iter, dmax, d/dble(nloc), gamma_max
 1          format(10x, i5, 2f15.10, 1p,d12.2,d12.2)
            call util_flush(luout)
         end if
         call ga_sync
c
**         if (abs(d-dprev)/dble(nloc) .lt. tol) goto 1000
**         dprev = d
         if (iter.gt.1 .and. gamma_max.lt.tol) goto 1000
         gamma_max = 0.0d0
c
c     Loop over pairs with as much parallelism as possible
c
         neven = nloc + mod(nloc,2)
         do set = 1, neven-1
            do pair = 1+ga_nodeid(), neven/2, ga_nnodes()
               call localize_pairs(neven, set, pair, ss, tt)
               if (tt .le. nloc) then
                  s = iloc(ss)
                  t = iloc(tt)
*                  write(luout,*) nloc, neven, set, pair, ss, tt, s, t
                  call ga_get(g_cibocc, 1, mnbf, s, s, c(1,1), 1)
                  call ga_get(g_cibocc, 1, mnbf, t, t, c(1,2), 1)
c     
c     Form rotation information
c     
                  ast = 0.0d0
                  bst = 0.0d0
                  do a = 1, natoms
                     if (.not. bas_ce2bfr(minbas, a, bflo, bfhi))
     &                    call errquit(pname//': basis ', 2,
     &                BASIS_ERR)
c     
                     qast = 0.0d0
                     qas  = 0.0d0
                     qat  = 0.0d0
                     do u = bflo, bfhi
                        qast = qast + c(u,2)*c(u,1) + c(u,1)*c(u,2)
                        qas  = qas  + c(u,1)*c(u,1)
                        qat  = qat  + c(u,2)*c(u,2)
                     end do
                     qast = qast * 0.5d0
c     
c                    ast = ast + qast**2 - 0.25d0*(qas - qat)**2
c                    bst = bst + qast*(qas - qat)
                     tmp = -1.0d0*(qas**4 + qat**4)
                     tmp = tmp + 6.0d0*((qas**2 + qat**2)* qast**2)
                     tmp = tmp + qas**3 * qat + qas * qat**3
                     ast = ast + 0.25d0*tmp
                     tmp = 4.0d0 * qast*(qas**3 - qat**3)
                     bst = bst + 0.25d0*tmp
                  end do
c     
                  gamma = 0.25d0*acos(-ast/sqrt(ast**2+bst**2))
                  gamma = sign(gamma,bst)
                  gamma_max = max(gamma_max, abs(gamma))
*                  if (iter .eq. 1 .and. abs(gamma).lt.0.01d0) then
*                     gamma = (util_random(0)-0.5d0)*3.14d0
*                  endif
c     
                  if (abs(gamma) .gt. gamma_tol) then
                     nrot = nrot + 1
                     cosg = cos(gamma)
                     sing = sin(gamma)
c     
c     Do the rotation of C 
c     
                     call drot(mnbf, c(1,1), 1, c(1,2), 1, cosg, sing)
                     call ga_put(g_cibocc, 1, mnbf, s, s, c(1,1), 1)
                     call ga_put(g_cibocc, 1, mnbf, t, t, c(1,2), 1)
                  end if
               end if
            end do
            call ga_sync
         end do
      end do
c
 1000 continue


c     the MOs in the orthonormal IAO basis must form an orthogonal
c     matrix. test that after localization

      if (.not. ga_create(MT_DBL, nloc, nloc , 'loc:tmp1',
     &  mnbf, 0, g_tmp1)) call errquit(pname//': tmp1',0, GA_ERR)

      call ga_dgemm('t', 'n', nloc, nloc, mnbf, 
     &  1.0d0, g_cibocc, g_cibocc, 0.0d0, g_tmp1)

      call ga_add_todiag(g_tmp1, -1.0d0)

      call ga_maxelt(g_tmp1, tmp)
      if (master) then
        write(luout,'(/1x,a,1x,f18.12)')
     &    'IBO loc: largest element of C(MO, trans) C(MO) -1: ',
     &    tmp
        write(luout,'(1x,a)')
     &    'should be zero, for LMOs in the IAO basis'
      end if

c     Following steps:
      
c     1. while we have array g_tmp1 allocated, also calculate the
c     localization transform. The starting MOs were saved in g_tmp2 in
c     the IAO basis

c     2. calculate LMO energies and order the LMO set by energy

c     3.  save the localization transform to a file, instead of
c     generating it again in the calling routine and saving it from
c     there

c     4. We then transform the CMOs to re-ordered LMOs in the
c     AO basis and store the transformed MOs back in g_c      


c     CMOs(iao, trans) * LMOs(iao) = localization transform -> g_tmp1
      call ga_dgemm('t', 'n', nloc, nloc, mnbf, 
     &  1.0d0, g_tmp2, g_cibocc, 0.0d0, g_tmp1)

      if (.not. ga_destroy(g_tmp2))
     &  call errquit(pname//': ga_destroy failed g_tmp2',0, GA_ERR)

c     transform MO energies. We assume that the first nloc MOs were
c     localized, and that the Fock matrix was diagonal in the basis of
c     input MOs. We also assume (this is tested upon entry) that
c     nloc is smaller than the column sizes of arays eval, c(:,1:2)
c     so we don't segfault when using c(:,:) as temp memory

      do jmo = 1,nloc
        call ga_get(g_tmp1, 1, nloc, jmo, jmo, c(1,1), 1)
        tmp = 0.0d0
        do imo = 1,nloc
          tmp = tmp + eval(imo) * c(imo,1)**2
        end do ! jmo
        c(jmo,2) = tmp ! temp storage for transformed MO energies
      end do ! imo
      do imo = 1,nloc
        eval(imo) = c(imo,2)
      end do
      
c     sort LMOs by energy. We'll have to do this on the master node, or
c     otherwise avoid having multiple nodes working on the same data

      if (master) then
      do imo = 1,nloc -1
        minval = eval(imo)
c       find lowest eval(kmo) below eval(imo)
        kmo = 0
        do jmo = imo+1, nloc
          if (eval(jmo).lt.minval) then
            kmo = jmo
            minval = eval(jmo)
          end if
        end do
c       if kmo > 0 we swap imo and jmo
        if (kmo.ne.0) then
          if (debug.and.master) write(luout,*) 'swap ',imo,kmo
c         swap eigenvalues
          swap = eval(imo)
          eval(imo) = eval(kmo)
          eval(kmo) = swap
c         swap columns of LMO transformation
          call ga_get(g_tmp1, 1, nloc, imo, imo, c(1,1), 1)
          call ga_get(g_tmp1, 1, nloc, kmo, kmo, c(1,2), 1)
          call ga_put(g_tmp1, 1, nloc, imo, imo, c(1,2), 1)
          call ga_put(g_tmp1, 1, nloc, kmo, kmo, c(1,1), 1)
c         swap LMOs in IAO basis, for the analysis below
          call ga_get(g_cibocc, 1, mnbf, imo, imo, c(1,1), 1)
          call ga_get(g_cibocc, 1, mnbf, kmo, kmo, c(1,2), 1)
          call ga_put(g_cibocc, 1, mnbf, imo, imo, c(1,2), 1)
          call ga_put(g_cibocc, 1, mnbf, kmo, kmo, c(1,1), 1)          
        end if
      end do
      end if ! master

      call ga_sync

c     write transformation to scratch file
      
      call util_file_name('lmotrans',.true.,.true.,lmotrans)
      if(.not.file_write_ga(lmotrans,g_tmp1)) call errquit
     &  (pname//': could not write lmotrans',0, DISK_ERR)

      if (debug.and.master) write (luout,*)
     &  pname//': lmos reordered and lmotrans written'

c     transform CMOs to re-ordered LMOs
      
      if (.not. ga_create(MT_DBL, nbf, nloc , 'loc:tmp2',
     &  nbf, 0, g_tmp2)) call errquit(pname//': tmp2',0, GA_ERR)      

      call ga_dgemm('n', 'n', nbf, nloc, nloc, 
     &  1.0d0, g_c, g_tmp1, 0.0d0, g_tmp2)

      call ga_copy(g_tmp2, g_c)

      if (debug.and.master) write (luout,*) pname//': g_c now has LMOs'

      if (.not. ga_destroy(g_tmp1))
     &  call errquit(pname//': ga_destroy failed g_tmp1',0, GA_ERR)

      if (.not. ga_destroy(g_tmp2))
     &  call errquit(pname//': ga_destroy failed g_tmp2',0, GA_ERR)

c     Analyze localization of each MO:
c     per LMO, a list of atomic populations is printed
c     in decreasing magnitude, with the polulations in parentheses

      if (master) then
         write(luout,'(/1x,a/)') 'IAO-IBO localized orbitals'
         do ss = 1, nloc
            s = iloc(ss)
            call ga_get(g_cibocc,  1, mnbf, s, s, c(1,1), 1)
            nlist = 0
            do a = 1, natoms
               if (.not. bas_ce2bfr(minbas, a, bflo, bfhi))
     &              call errquit('localized: basis ', 3,
     &          BASIS_ERR)
               qas  = 0.0d0
               do u = bflo, bfhi
                  qas  = qas  + c(u,1)*c(u,1)
               end do
               if (abs(qas) .gt. 0.01d0) then
                  nlist = nlist + 1
                  list(nlist) = a
                  pop(nlist) = qas
               end if
            end do
            do u = 1, nlist
               do t = 1, u-1
                  if (abs(pop(t)).lt.abs(pop(u))) then
                     tmp = pop(u)
                     pop(u) = pop(t)
                     pop(t) = tmp
                     tt = list(u)
                     list(u) = list(t)
                     list(t) = tt
                  end if
               end do
            end do
            write(luout,77) s, eval(s),
     &        occ(s),(list(a), pop(a), a=1,nlist)
 77         format(i5, 1x, f14.6,1x, f5.3, 1x,100(2x,i4,'(',f5.2,')'))
          end do
          write(luout,*)
         call util_flush(luout)
       end if

       if (.not. ga_destroy(g_cibocc))
     &   call errquit(pname//': ga_destroy failed g_cibocc',0, GA_ERR)

      if (debug.and.master) write(luout,*) 'leaving '//pname

      end      
c
c     =================================================================
c      
      subroutine localize_pairs(n, set, pair, left, right)
      implicit none
      integer n, set, pair, left, right
c
c     This routine returns maximally overlapped independent pairs.
c     Use it with the following code fragment to generate all 
c     unique pairs (n*(n-1)/2 of them).  The loop over pair can
c     be executed in parallel.
c
c     .    neven = n + mod(n,2)
c     .    do set = 1, neven-1
c     .       do pair = 1, neven/2
c     .          call pairs(neven, set, pair, i, j)
c*    .           write(6,*) ' GOT ', set, pair, i, j
c     .          if (j. le. n) then
c     .             p(i,j) = 1
c     .             p(j,i) = 1
c     .          end if
c     .       end do
c     .    end do
c
      left = pair - set + 1
      if (left .le. 0) left = n - 1 + left
      if (pair .eq. 1) then
         right = n
      else
         right = n - pair + 1 - set + 1
         if (right .le. 0) right = n - 1 + right
      end if
c
      end
