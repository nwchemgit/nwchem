      logical function rohf_diis_solve(rtdb, energy, eone, etwo, ecosmo,
     $                                 enrep)
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "crohf.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "util.fh"
#include "cscf.fh"
#include "cosmo.fh"
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
      logical oprint_parm, nodezero
      integer mxdiisbas
c     
c     Functions
c     
      integer ga_create_atom_blocked,g_oexch,g_cfock
      external ga_create_atom_blocked, ga_iter_lsolve
c
c     cosmo code
c
      logical odbug
      logical osome
      double precision ecosmo !< [Output] The COSMO solvation energy
c
      double precision egas
      double precision esol
c
      integer itype_wfn
c     
      vlen = crohf_vlen
      aoff = nclosed + 1
      aend = nclosed + nopen
      voff = nclosed + nopen + 1
c
      odbug=.false.
      if(odbug) then
         write(6,*) 'in rohf_diis_solve ... with node id =',
     $              ga_nodeid()
      endif
c
      oprint_parm = util_print('parameters', print_default)
      if (.not.rtdb_get(rtdb,'scf:diisbas',MT_INT, 1, mxdiisbas ))
     $   mxdiisbas = DFLT_MAX_DIIS_BAS
c
c   Print info
c
      nodezero = ga_nodeid().eq.0
      if (nodezero.and. oprint_parm) then
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
     $     call errquit('diis:grad', 0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, (nbf*nbf), (2*mxdiisbas),
     $     'diis:basis', 0, 0, g_gdiis)) call errquit('diis: work', 0,
     &       GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:fock',
     $     0, 0, g_fock)) call
     $     errquit('rohf_diis_solve: ga failed for fock', 0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:over',
     $     0, 0, g_over)) call
     $     errquit('rohf_diis_solve: ga failed for over', 0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:xorth',
     $     0, 0, g_xorth)) call
     $     errquit('rohf_diis_solve: ga failed for xorth', 0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:tmp',
     $     0, 0, g_tmp)) call
     $     errquit('rohf_diis_solve: ga failed for tmp', 0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, nbf, nbf, 'rohf_diis_solve:tmp2',
     $     0, 0, g_tmp2)) call
     $     errquit('rohf_diis_solve: ga failed for tmp2', 0, GA_ERR)
*ga:1:0
      if (.not. ga_create(MT_DBL, (nbf*nbf), 1, 'rohf_diis_solve:err',
     $     0, 0, g_err)) call
     $     errquit('rohf_diis_solve: ga failed for error', 0, GA_ERR)
c
c  Generate overlap and orthonormalization matrix
c
      call ga_zero(g_over)
      call int_1e_ga( basis, basis, g_over, 'overlap', .false. )
      call rohf_ao_orthonorm( nbf, g_over, g_xorth )
      call ga_sync()
c
 1000 continue
c     
c  SCF-DIIS loop
c
      converged = .false.
      ndiis = 0
      call ga_zero(g_err)
      do iter=1,maxiter
        call rohf_energy( g_movecs, eone, etwo, enrep, ecosmo, energy,
     $                    g_grad, rtdb)
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
c
c   Print energy
c
        if(iter.eq.1) then
           if(cosmo_on.and.cosmo_phase.eq.1) then
              if(nodezero) write(6,909)
           else
              if(nodezero) write(6,910)
           endif
        endif
 909    format(2x,'   COSMO gas phase')
 910    format(2x,'   COSMO solvation phase')
c
        if(nodezero) then
           write(6,911) iter, ndiis, energy, gnorm, diis_emax
        endif
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
        if ( util_print('evp_input', print_never) ) then
            call scf_dump_evp_input(rtdb,g_fock,g_over)
        endif
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
        call ga_pdsygv(g_fock,g_over, g_tmp, dbl_mb(k_eval))
#else
        call ga_diag( g_fock, g_over, g_tmp, dbl_mb(k_eval))
#endif
#else 
        call ga_diag_seq( g_fock, g_over, g_tmp, dbl_mb(k_eval))
#endif
        call ga_copy(g_tmp, g_movecs)
      enddo
c     
c     End SCF minimisation
c     
 333  continue
      if (nodezero) call util_flush(6)
c     
c     *ALWAYS* return canonical MOs ... do not do aufbau
c     
      if (.not. ga_destroy(g_tmp)) call errquit
     $     ('rohf_diis_solve: ga_destroy tmp', 0, GA_ERR)
      if (.not. ga_destroy(g_tmp2)) call errquit
     $     ('rohf_diis_solve: ga_destroy tmp2', 0, GA_ERR)
      call rohf_canon(.false., .false.)
      call movecs_fix_phase(g_movecs)
      call scf_movecs_write(rtdb)
      if (olagr) then
         g_cfock = ga_create_atom_blocked(geom, basis, 'rohf: cfock')
         g_oexch = ga_create_atom_blocked(geom, basis, 'rohf: ocoul')
         call rohf_readfock(nopen,g_oexch,g_cfock)
         if (.not. ga_duplicate(g_movecs, g_tmp, 'scfmsa'))
     $        call errquit('rohf_lagr: gaduplicate?',0, GA_ERR)
         call two_index_transf(g_cfock, g_movecs, g_movecs, g_tmp, 
     C        crohf_g_fcv)
         
         if (nopen .gt. 0) then
            call two_index_transf(g_oexch, g_movecs, g_movecs, g_tmp, 
     C           crohf_g_fpv)
         endif
         if (.not. ga_destroy(g_oexch)) call errquit
     $        ('rohf_diis_solve: ga_destroy oexch', 0, GA_ERR)
         if (.not. ga_destroy(g_cfock)) call errquit
     $        ('rohf_diis_solve: ga_destroy cfock', 0, GA_ERR)
         if (.not. ga_destroy(g_tmp)) call errquit
     $        ('rohf_diis_solve: ga_destroy tmp', 0, GA_ERR)
         call rohf_lagr()
      endif
c
c     ----- if sol_phase=.T. go back and do solution phase -----
c
      odbug=.false.
      if(cosmo_on.and.cosmo_phase.eq.1) then
         converged =.false.
         egas      =energy
         cosmo_phase = 2
         if(odbug) then
            write(6,*) 'gas_phase done, do sol-phase now ...'
         endif
         call ga_sync()
         go to 1000
      elseif(cosmo_on.and.cosmo_phase.eq.2) then
         esol      =energy
         if(nodezero) then
            write(6,*) '                 COSMO solvation results'
            write(6,*) '                 -----------------------'
            write(6,*) 
            write(6,912) egas
            write(6,913) esol
            write(6,914) (egas-esol),(egas-esol)*627.509451d+00
 912     format('                 gas phase energy = ',f20.10)
 913     format('                 sol phase energy = ',f20.10)
 914     format(' (electrostatic) solvation energy = ',f20.10,
     $          ' (',f8.2,' kcal/mol)'                       )
         endif
      endif
      call ga_sync()
c     
c     Free globals
c     
      if (.not. ga_destroy(g_grad)) call errquit
     $     ('rohf_diis_solve: ga_destroy grad', 0, GA_ERR)
      if (.not. ga_destroy(g_fock)) call errquit
     $     ('rohf_diis_solve: ga_destroy fock', 0, GA_ERR)
      if (.not. ga_destroy(g_err)) call errquit
     $     ('rohf_diis_solve: ga_destroy err', 0, GA_ERR)
      if (.not. ga_destroy(g_gdiis)) call errquit
     $     ('rohf_diis_solve: ga_destroy gdiis', 0, GA_ERR)
      if (.not. ga_destroy(g_over)) call errquit
     $     ('rohf_diis_solve: ga_destroy over', 0, GA_ERR)
      if (.not. ga_destroy(g_xorth)) call errquit
     $     ('rohf_diis_solve: ga_destroy xorth', 0, GA_ERR)
c     
      if (nodezero) call util_flush(6)
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
      subroutine rohf_ao_orthonorm( nbf, g_over, g_xorth )
      implicit none
#include "errquit.fh"
c
c   Construct AO orthonormalization matrix
c   using symmetric orthogonalization.
c   Does *not* handle linear dependence, yet
c
#include "mafdecls.fh"
#include "global.fh"      
      integer nbf
      integer g_over
      integer g_xorth
c
      integer l_e, k_e, i
      double precision xx
c
      if (.not.ma_push_get(MT_DBL, nbf, 'overlap', l_e, k_e ))
     $     call errquit('rohf_ao_orthonorm: cannot allocate over',0,
     &       MA_ERR)
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
      call ga_pdsyev(g_over, g_xorth, dbl_mb(k_e), 0)
#else
      call ga_diag_std( g_over, g_xorth, dbl_mb(k_e))
#endif
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
     $     call errquit('rohf_ao_orthonorm: cannot pop stack',0,
     &       MA_ERR)
      return
      end
