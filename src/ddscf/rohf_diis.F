      logical function rohf_diis_solve(rtdb, energy, eone, etwo, enrep)
C     $Id: rohf_diis.F,v 1.12 2000-05-22 23:27:26 mdupuis Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "crohf.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "util.fh"
#include "cscf.fh"
c     
c     Solve the ROHF equations using DIIS
c     
c     Arguments
c     
      integer rtdb
      double precision energy                  ! Return ROHF energy
      double precision eone, etwo, enrep       ! Return ROHF energy contribs
c
      integer DFLT_MAX_DIIS_BAS
      parameter(DFLT_MAX_DIIS_BAS=5)
c     
c     Local GA handles 
c     
      integer g_fock
      integer g_xorth
      integer g_tmp
      integer g_tmp2
      integer g_gdiis
      integer g_grad
      integer g_over
      integer g_err
c     
c     Local variables
c     
      double precision gnorm                   ! measures for convergence
      double precision diis_emax
      logical converged, odo_diis
      integer vlen, ndiis
      integer voff, aoff, aend
      logical oprint_parm
      integer mxdiisbas
c     
c     Functions
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked, ga_iter_lsolve
c
c     cosmo code
c
      integer g_cos1e
      dimension g_cos1e(3)
      logical odbug
      logical ocosmo
      logical osome
      double precision ecosmo
c     
      vlen = crohf_vlen
      aoff = nclosed + 1
      aend = nclosed + nopen
      voff = nclosed + nopen + 1
c
c
c
      oprint_parm = util_print('parameters', print_default)
      if (.not.rtdb_get(rtdb,'scf:diisbas',MT_INT, 1, mxdiisbas ))
     $   mxdiisbas = DFLT_MAX_DIIS_BAS
c
c   Print info
c
      if (ga_nodeid().eq.0.and. oprint_parm) then
         write(6,1) gnorm_tol, maxiter, tol2e, mxdiisbas
         write(6,111)
 1       format(//,
     $        1x,'----------------------------------------------',/
     $        1x,'               DIIS ROHF',//,
     $        1x,'Convergence threshold     :',9x,1p,e10.3,0p,/,
     $        1x,'Maximum no. of iterations :',9x,i4,/,
     $        1x,'Integral*density screening:',9x,1p,e10.3,0p,/,
     $        1x,'DIIS subspace             :',9x,i4,/)
 111     format(
     $        1x,'----------------------------------------------',/)
         call util_flush(6)
      endif
c     
c     Allocate global arrays.  
c     
*ga:1:0
      if (.not. ga_create(MT_DBL, vlen, 1, 'diis:grad', 0, 0, g_grad))
     $     call errquit('diis:grad', 0)
*ga:1:0
      if (.not. ga_create(MT_DBL, (nbf*nbf), (2*mxdiisbas),
     $     'diis:basis', 0, 0, g_gdiis)) call errquit('diis: work', 0)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:fock',
     $     0, 0, g_fock)) call
     $     errquit('rohf_diis_solve: ga failed for fock', 0)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:over',
     $     0, 0, g_over)) call
     $     errquit('rohf_diis_solve: ga failed for over', 0)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:xorth',
     $     0, 0, g_xorth)) call
     $     errquit('rohf_diis_solve: ga failed for xorth', 0)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:tmp',
     $     0, 0, g_tmp)) call
     $     errquit('rohf_diis_solve: ga failed for tmp', 0)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:tmp2',
     $     0, 0, g_tmp2)) call
     $     errquit('rohf_diis_solve: ga failed for tmp2', 0)
*ga:1:0
      if (.not. ga_create(MT_DBL, (nbf*nbf), 1, 'rohf_diis_solve:err',
     $     0, 0, g_err)) call
     $     errquit('rohf_diis_solve: ga failed for error', 0)
c
c  Generate overlap and orthonormalization matrix
c
      call ga_zero(g_over)
      call int_1e_ga( basis, basis, g_over, 'overlap', .false. )
      call rohf_ao_orthonorm( nbf, g_over, g_xorth )
c     
c  SCF-DIIS loop
c
      converged = .false.
      ndiis = 0
      call ga_zero(g_err)
      do iter=1,maxiter
c
c     ----- cosmo charges and energy contribution -----
c
        if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo)) then
           if(ocosmo) then
              odbug=.false.
              if(odbug) then
                 osome=.true.
              else
                 osome=.false.
              endif
              odbug=odbug.and.ga_nodeid().eq.0
              osome=osome.and.ga_nodeid().eq.0
c
              call cosmo_charges(rtdb,basis,geom,ecosmo,iter,osome)
c
              odbug=.false.
              osome=.false.
           endif
        endif
        call ga_sync()
c
c     ----- cosmo 1-electron potential operator -----
c
        if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo)) then
           if(ocosmo) then
c
              odbug=.true.
              if(odbug) then
                 osome=.true.
              else
                 osome=.false.
              endif
              odbug=odbug.and.ga_nodeid().eq.0
              osome=osome.and.ga_nodeid().eq.0
c
              if(odbug) then
                 write(6,*) 'cos1e-pot ... try g_over'
              endif
              call ga_zero(g_over)
              call int_1e_ga(basis,basis,g_over,'overlap',.false.)
              if(odbug) then
                 call ga_print(g_over)
              endif
c
              g_cos1e(1)  = ga_create_atom_blocked(geom, basis, 
     $                                         'rohf_diis_solve:cos1e')
              if(odbug) then
                 write(6,*) 'cos1e-pot ... try g-cos1e'
              endif
              call ga_zero(g_cos1e(1))
              call int_1e_ga(basis,basis,g_cos1e(1),'cos_chg_pot',
     $                                                    .false.)
              if(odbug) then
                 call ga_print(g_cos1e(1))
              endif
c
              odbug=.false.
              osome=.false.
c
              if (.not. ga_destroy(g_cos1e)) call errquit
     $           ('rohf_diis_solve: ga_destroy cos1e', 0)
           endif
        endif
        call ga_sync()
c
        call rohf_energy( g_movecs, eone, etwo, enrep, energy, g_grad)
c
c     ----- add cosmo energy contribution -----
c
        if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo)) then
           if(ocosmo) then
              energy=energy+ecosmo
c
              osome=.true.
              osome=osome.and.ga_nodeid().eq.0 
              osome=.false.
              if(osome) then
                 write(6,9999) (energy-ecosmo),ecosmo,energy
 9999 format(' energ0 = ',f20.10,/,' ecosmo = ',f20.10,/,
     $       ' energy = ',f20.10)
              endif
           endif
        endif
        call ga_sync()
c
c   Check convergence
c
        gnorm = sqrt(ga_ddot(g_grad, g_grad))
        converged = (gnorm.lt.gnorm_tol)
        if (converged) goto 333
        call rohf_get_fock( g_fock )
c
c   Recover AO Fock matrix
c
c                  t
c   F   =  S.C.F  C.S
c    ao         mo
c
        call ga_dgemm( 'n', 't', nbf, nbf, nbf, 1.d0, g_fock,
     $                 g_movecs, 0.d0, g_tmp )
        call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_tmp,
     $                 g_over, 0.d0, g_fock )
        call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_movecs,
     $                 g_fock, 0.d0, g_tmp )
        call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_over,
     $                 g_tmp, 0.d0, g_fock )
c
c   DIIS Error vector
c
        call rohf_errvec( g_fock, g_movecs, g_over, g_xorth,
     $                    g_tmp, g_tmp2, g_err )
        call ga_maxelt( g_err, diis_emax )
        write(6,911) iter, ndiis, energy, gnorm, diis_emax
 911    format(5x,i5,i3,5x,f20.12,5x,2e10.2)
c
c   DIIS extrapolation
c
        odo_diis = (diis_emax.lt.0.5d0)
        call rohf_diis( nbf, iter, ndiis, g_err, g_gdiis,
     $                  g_fock, odo_diis )
c
c   Diagonalize
c
#if defined(PARALLEL_DIAG)
        call ga_diag( g_fock, g_over, g_tmp, dbl_mb(k_eval))
#else 
        call ga_diag_seq( g_fock, g_over, g_tmp, dbl_mb(k_eval))
#endif
        call ga_copy(g_tmp, g_movecs)
      enddo
c     
c     End SCF minimisation
c     
 333  continue
      if (ga_nodeid().eq.0) call util_flush(6)
c     
c     *ALWAYS* return canonical MOs ... don't do aufbau
c     
      call rohf_canon(.false., .false.)
      call movecs_fix_phase(g_movecs)
      call scf_movecs_write(rtdb)
      if (olagr) call rohf_lagr()
c     
c     Free globals
c     
      if (.not. ga_destroy(g_grad)) call errquit
     $     ('rohf_diis_solve: ga_destroy grad', 0)
      if (.not. ga_destroy(g_fock)) call errquit
     $     ('rohf_diis_solve: ga_destroy fock', 0)
      if (.not. ga_destroy(g_err)) call errquit
     $     ('rohf_diis_solve: ga_destroy err', 0)
      if (.not. ga_destroy(g_gdiis)) call errquit
     $     ('rohf_diis_solve: ga_destroy gdiis', 0)
      if (.not. ga_destroy(g_over)) call errquit
     $     ('rohf_diis_solve: ga_destroy over', 0)
      if (.not. ga_destroy(g_xorth)) call errquit
     $     ('rohf_diis_solve: ga_destroy xorth', 0)
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('rohf_diis_solve: ga_destroy tmp', 0)
      if (.not. ga_destroy(g_tmp2)) call errquit
     $     ('rohf_diis_solve: ga_destroy tmp2', 0)
c     
      if (ga_nodeid().eq.0) call util_flush(6)
      call ga_sync()
c     
      rohf_diis_solve = converged
c     
      end




c
c  Generic DIIS extrapolation routine
c  Cyclic storage of error-vectors and Fock matrices
c
      subroutine rohf_diis( nbf, iter, ndiis, g_err, g_bas, g_fock,
     $                      odo_diis )
      implicit none
#include "global.fh"      
      integer nbf, iter, ndiis
      integer g_err
      integer g_fock
      integer g_bas
      logical odo_diis
c
      integer mxbas, nbas, vlen, itype
      integer i, j, ij, ji
      integer info
      double precision xx(500), yy(100), tmp(100), zz
c
c  Copy new error vector to basis
c
      call ga_inquire( g_bas, itype, vlen, mxbas )
      mxbas = mxbas/2
      if (odo_diis) then
        ndiis = ndiis + 1
        i = mod( (ndiis-1), mxbas ) + 1
        nbas = min(ndiis, mxbas)
        call ga_copy_patch( 'n', g_err, 1, vlen, 1, 1,
     $                           g_bas, 1, vlen, i, i )
        call ga_copy_patch( 'n', g_fock, 1, nbf, 1, nbf,
     $                           g_bas, 1, vlen, i+mxbas, i+mxbas )
      endif
      if ((iter.le.3).or.(.not.odo_diis)) return
c
c  Construct error metric
c
      call dfill((nbas+1)*(nbas+1),-1.d0,xx,1)
      xx(1) = 0.d0
      do i=1,nbas
        do j=1,i
          zz = ga_ddot_patch( g_bas, 'n', 1, vlen, i, i,
     $                        g_bas, 'n', 1, vlen, j, j ) 
          ij = (i-1+1)*(nbas+1) + j + 1
          ji = (j-1+1)*(nbas+1) + i + 1
          xx(ij) = zz
          xx(ji) = zz
        enddo
      enddo
      call dfill((nbas+1),0.d0,yy,1)
      yy(1) = -1.d0
      call dgesv((nbas+1), 1, xx, (nbas+1), tmp, yy, (nbas+1), info)
c
c  Construct extrapolated vector
c
      call ga_zero(g_err)
      do i=1,nbas
        call ga_dadd_patch( yy(i+1), g_bas, 1, vlen, mxbas+i, mxbas+i,
     $                      1.d0,    g_err, 1, vlen, 1, 1,
     $                               g_err, 1, vlen, 1, 1 )
      enddo
      call ga_copy_patch( 'n', g_err, 1, vlen, 1, 1,
     $                         g_fock, 1, nbf, 1, nbf )
      odo_diis = .true.
      return
      end




c
c  Generate Pulay's error vector from AO fock matrix
c
c             -1/2 t             1/2
c       e = (S    ).(FDS - SDF).S
c
c
      subroutine rohf_errvec( g_fock, g_vecs, g_over, g_xorth,
     $                        g_tmp1, g_tmp2, g_err )
      implicit none
#include "global.fh"
#include "cscf.fh"
      integer g_fock
      integer g_vecs
      integer g_over
      integer g_xorth
      integer g_tmp1
      integer g_tmp2
      integer g_err
c
c  ROHF density
c
      call ga_matmul_patch( 'n', 't', 2.d0, 0.d0,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_movecs, 1, nclosed, 1, nbf,
     $                       g_tmp2, 1, nbf, 1, nbf )
      if (nopen.gt.0) then
        call ga_matmul_patch( 'n', 't', 1.d0, 1.d0,
     $                       g_movecs, 1, nbf, nclosed+1, nclosed+nopen,
     $                       g_movecs, nclosed+1, nclosed+nopen, 1, nbf,
     $                       g_tmp2, 1, nbf, 1, nbf )
      endif
c
c                         t
c  FDS - SDF  =  FDS - FDS
c
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_fock,
     $               g_tmp2, 0.d0, g_tmp1 )
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_tmp1,
     $               g_over, 0.d0, g_tmp2 )
      call ga_transpose( g_tmp2, g_tmp1 )
      call ga_dadd( 1.d0, g_tmp2, -1.d0, g_tmp1, g_tmp2 )
c
c       -1/2          t   -1/2
c  e = S   .(FDS - FDS ).S
c
      call ga_dgemm( 't', 'n', nbf, nbf, nbf, 1.d0, g_xorth,
     $                g_tmp2, 0.d0, g_tmp1 )
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_tmp1,
     $                g_xorth, 0.d0, g_tmp2 )
      call ga_copy_patch( 'n', g_tmp2, 1, nbf, 1, nbf,
     $                         g_err, 1, (nbf*nbf), 1, 1)
c
c
c
      return
      end






c
c
c   Construct AO orthonormalization matrix
c   using symmetric orthogonalization.
c   Does *not* handle linear dependence, yet
c
c
      subroutine rohf_ao_orthonorm( nbf, g_over, g_xorth )
      implicit none
#include "mafdecls.fh"
#include "global.fh"      
      integer nbf
      integer g_over
      integer g_xorth
c
      integer l_e, k_e, i
      double precision xx
c
c
c
      if (.not.ma_push_get(MT_DBL, nbf, 'overlap', l_e, k_e ))
     $     call errquit('rohf_ao_orthonorm: cannot allocate over',0)
#if defined(PARALLEL_DIAG)
      call ga_diag_std( g_over, g_xorth, dbl_mb(k_e))
#else 
      call ga_diag_std_seq( g_over, g_xorth, dbl_mb(k_e))
#endif
c               1/2
c   X   = U  / s
c    ij    ij   j
c
      do i=ga_nodeid()+1,nbf,ga_nnodes()
        xx = 1.d0/sqrt(dbl_mb(k_e+i-1))
        call ga_dscal_patch( g_xorth, 1, nbf, i, i, xx )
      enddo
c
c
c
      if (.not.ma_pop_stack(l_e))
     $     call errquit('rohf_ao_orthonorm: cannot pop stack',0)
      return
      end
