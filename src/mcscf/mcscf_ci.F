c
c   MCSCF Secular solution front-end
c
      subroutine mcscf_cisolve( rtdb, geom, basis, nbf, nclosed,
     $                          nact, nsym, orbsym, e0, eps, g_h,
     $                          g_coul, citol, ociread, oprint,
     $                          ociprint, dm1, dm2, eci )
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "mcscfprof.fh"
c
      integer rtdb, geom, basis
      integer nbf, nclosed, nact
      integer nsym
      integer orbsym(nbf)
      double precision eps(nbf)
      double precision e0
      integer g_h
      integer g_coul
      double precision citol
      logical ociread
      logical oprint
      logical ociprint
      double precision dm1(nact,nact)
      double precision dm2(nact,nact,nact,nact)
      double precision eci
c
#ifdef OLD_DAVIDSON
      integer mxworkvec
      parameter(mxworkvec=16)
#else
      integer mxworkvec
      parameter(mxworkvec=5)
#endif      
      integer nela, nelb, nactel, multip, symstate
      integer actsym(detci_maxorb), actmap(detci_maxorb)
      integer guess_cfg(detci_maxelec_tot*detci_maxguess_cfg)
      integer itmp(detci_maxelec_tot)
      integer nguesscfg
      double precision guess_coeff(6)
      double precision ciener
      double precision acteps(detci_maxorb)
      double precision mxcnt
      integer maxiter
      integer ijmap(detci_maxorb*detci_maxorb)
      integer i, iprint, hasprint
      integer ntij, n2
      integer l_h1, k_h1
      integer l_eri, k_eri
      integer g_civec, g_wvec(mxworkvec)
      character*20 vstr
      character*255 civec_in
      logical civec_saved
c
      external dummy_mxv
      external dummy_precon
      external detci_sigma
      external detci_ciprecon
      logical file_read_ga, file_write_ga
      external file_read_ga, file_write_ga
c
      save hasprint
      data hasprint/0/
C      data civec_in/'civec'/
c
c
      iprint = 1
      if (hasprint.ne.0 .or. (.not. oprint)) iprint = 0
      hasprint = hasprint + 1
      if (.not.rtdb_get(rtdb,'mcscf:cimaxiter',MT_INT,1,maxiter))
     $  maxiter = 100
c
c   Retrieve spin-related info from database
c
      if (.not.rtdb_get(rtdb,'mcscf:multiplicity',MT_INT,1,multip))
     $  multip = 0
      if (.not.rtdb_get(rtdb,'mcscf:nactelec',MT_INT,1,nactel))
     $  nactel = 0
      nela = (nactel + multip - 1)/2
      nelb = nactel - nela
      if ((mod((nactel + multip - 1),2).ne.0).or.
     $    (nela.lt.0).or.(nelb.lt.0))
     $    call errquit('detci: incompatible elec and spin',0, INPUT_ERR)
      if ((nela.gt.detci_maxelec).or.(nelb.gt.detci_maxelec).or.
     $    ((nela+nelb).gt.detci_maxelec_tot))
     $    call errquit('mcscf_ci: too many electrons - redimension',0,
     &       INPUT_ERR)
c
c   Retrieve symmetry stuff
c
      if (.not.rtdb_get(rtdb,'mcscf:statesymmetry',MT_INT,1,symstate))
     $   symstate = 1
c
c   Remap active indices into symmetry blocks
c
      do i=1,nact
        actmap(i) = i
        actsym(i) = orbsym(nclosed+i)
        acteps(i) = eps(nclosed+i)
      enddo
      call detci_ijmap( nact, nsym, actsym, ntij, ijmap )
c
c   Allocate local space for 1e and 2e MO integrals
c   Push global MO integrals into local space
c
      n2 = nact*nact
      if (.not.ma_push_get(MT_DBL,n2,'1e hamil',l_h1,k_h1))
     $  call errquit('mcscf_cisolve: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,(n2*n2),'eri',l_eri,k_eri))
     $  call errquit('mcscf_cisolve: cannot allocate',0, MA_ERR)
      call dfill(n2,0.d0,dbl_mb(k_h1),1)
      call dfill((n2*n2),0.d0,dbl_mb(k_eri),1)
      call mcscf_getmolocal( nclosed, nact, nbf, g_h, g_coul,
     $                       ntij, ijmap, dbl_mb(k_h1),
     $                       dbl_mb(k_eri) )
c
c   Initialize DetCI internals
c
      call detci_init( nact, nela, nelb, nsym, symstate,
     $                 actsym, iprint, acteps,
     $                 dbl_mb(k_h1), dbl_mb(k_eri))
c
c   Create CI-vectors
c
      if (.not.ga_civec_create('CI vector',g_civec))
     $   call errquit('detci: cannot create CI vector',0, GA_ERR)
      do i=1,mxworkvec
        vstr = 'CI work ' // char(ichar('0')+i)
        if (.not.ga_civec_create(vstr,g_wvec(i)))
     $    call errquit('detci: cannot create CI vector',0, GA_ERR)
        call ga_zero(g_wvec(i))
      enddo
c
c   Guess CI vector
c   Can manually set the guess config. in the database
c   otherwise generate guess via Aufbau rules with correct
c   symmetry. Automatic
c
      call util_file_name('civec', .false.,.false.,civec_in)
      if (ga_nodeid() .eq. 0) inquire(file=civec_in,exist=civec_saved)
      call ga_brdcst(666, civec_saved, ma_sizeof(mt_log, 1, mt_byte), 0)
c
      if ((ociread).and.(civec_saved)) then
        if (.not.(file_read_ga(civec_in,g_civec)))
     $    call errquit('mcscf_ci: cannot read CI vector',0, DISK_ERR)
      else
        if (rtdb_get(rtdb,'mcscf:guess',MT_INT,(nela+nelb),itmp)) then
          do i=1,(nela+nelb)
            guess_cfg(i) = itmp(i)
          enddo
          call util_isort( nela, guess_cfg )
          call util_isort( nelb, guess_cfg(nela+1) )
        else
          call detci_guesscfg( nact, nsym, nela, nelb, actsym,
     $                         acteps, symstate, guess_cfg, nguesscfg )
        endif
        do i=1,nguesscfg
           guess_coeff(i) = 1.d0/nguesscfg
        enddo
        call detci_guess( nguesscfg, guess_cfg, guess_coeff,
     $                    g_civec, g_wvec(1) )
        if (util_print('ci guess',print_debug)) then
          if (ga_nodeid().eq.0) write(6,222)
 222      format(/,'Spin-adapted guess vector')
          call detci_print(g_civec,1.d-2)
        endif
      endif
c
c   Davidson diag
c
      if ((omcscfprof).and.(ociprof)) call pstat_on(ps_cidiag)
#ifdef OLD_DAVIDSON
      call detci_dvddiag0( detci_sigma, detci_ciprecon, mxworkvec,
     $                     maxiter, citol, e0, ciener, g_wvec, g_civec,
     $                     oprint )
#else
      call detci_dvddiag( detci_sigma, detci_ciprecon, maxiter,
     $                    citol, oprint, g_civec, g_wvec(1), g_wvec(2),
     $                    g_wvec(3), g_wvec(4), g_wvec(5), ciener )
#endif
      eci = e0 + ciener
      if ((omcscfprof).and.(ociprof)) call pstat_off(ps_cidiag)
c
c   Free work vectors ASAP
c
      do i=1,mxworkvec
        if (.not.ga_destroy(g_wvec(i)))
     $    call errquit('mcscf_cisolve: cannot destroy CI vector',0,
     &       GA_ERR)
      enddo
c
c   Print out converged ground state CI vector
c
      if (ociprint) then
        if (ga_nodeid().eq.0) write(6,885) 
 885    format(/,10x,'Converged CI vector')
        call detci_print( g_civec, 2.d-2 )
      endif
c
c   Check symmetry contamination of soln vector
c
      if (nsym.gt.1) then
        call detci_symmproj( g_civec, mxcnt, .false. )
        if (mxcnt.gt.1.d-14) then
          if (ga_nodeid().eq.0) write(6,881) mxcnt
 881      format(/,'Warning - CI symmetry contamination:',e12.3)
c$$$        else if (ociprint) then
c$$$          if (ga_nodeid().eq.0) write(6,882) mxcnt
c$$$ 882      format(//,10x,'CI Symmetry contamination:',e12.3)
        endif
      endif
c
c   Compute MO densities
c
      if ((omcscfprof).and.(ociprof)) call pstat_on(ps_cidens)
      call detci_density( g_civec, dm1, dm2 )
      if ((omcscfprof).and.(ociprof)) call pstat_off(ps_cidens)
c
c   Save CI vector for restart
c
      if (.not. file_write_ga(civec_in,g_civec))
     $  call errquit('mcscf_ci: error in CI vector save',0, DISK_ERR)
c
c   Release CI vectors
c
      if (.not.ga_destroy(g_civec))
     $   call errquit('mcscf_cisolve: cannot destroy CI vector',0,
     &       GA_ERR)
c
c   Release DetCI internals
c
      call detci_free()
c
c   Clean up
c
      if (.not.ma_pop_stack(l_eri))
     $  call errquit('mcscf_cisolve: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_h1))
     $  call errquit('mcscf_cisolve: cannot pop stack',0, MA_ERR)
c
c
      return
      end






      subroutine mcscf_getmolocal( nclosed, nact, nbf, g_h, g_coul,
     $                             ntij, ijmap, h1, eri )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nclosed
      integer nact
      integer nbf
      integer g_h
      integer g_coul
      integer ntij
      integer ijmap(nact,nact)
      double precision h1(ntij)
      double precision eri(ntij,ntij)
c
      integer i, j, k, l, gij, ij, kl, ii, jj
      integer coff, kk, klo, khi
      integer l_t,k_t
      double precision xx
      integer type, ncol, nrow
      logical ojact
c
c
c
      call ga_inquire(g_coul, type, ncol, nrow )
      ojact = (ncol.eq.(nact*nact))
      if ((.not.(ojact)).and.(ncol.ne.(nbf*nbf)))
     $  call errquit('mcscf_getmolocal: wrong J dimension',0, INPUT_ERR)
      if (.not.ma_push_get(MT_DBL,(nact*nact),'t',l_t,k_t))
     $  call errquit('mcscf_getmolocal: cannot allocate',0, MA_ERR)
      call dfill(ntij,0.d0,h1,1)
      call dfill((ntij*ntij),0.d0,eri,1)

      coff = nclosed*nbf + nclosed
      do i=1,nact
        do j=1,i
          gij = ((i-1)*i)/2 + j
          ij = ijmap(i,j)
          do k=1,nact
            kk = (k-1)*nact
            if (ojact) then
              klo = (k-1)*nact + 1
              khi = (k-1)*nact + nact
            else
              klo = coff + (k-1)*nbf + 1
              khi = coff + (k-1)*nbf + nact
            endif
            call ga_get(g_coul,klo,khi,gij,gij,dbl_mb(k_t+kk),1)
          enddo
          do k=1,nact
            kk = (k-1)*nact
            do l=1,k
              kl = ijmap(k,l)
              eri(ij,kl) = dbl_mb(k_t+kk+l-1)
            enddo
          enddo
          ii = nclosed + i
          jj = nclosed + j
          call ga_get(g_h,ii,ii,jj,jj,xx,1)
          h1(ij) = xx
        enddo
      enddo

      if (.not.ma_pop_stack(l_t))
     $  call errquit('mcscf_getmolocal: cannot pop stack',0, MA_ERR)
      return
      end

          




      subroutine dummy_mxv( g_a, g_b )
      implicit none
#include "global.fh"
      integer g_a, g_b

      call ga_copy( g_a, g_b )
      return
      end



      subroutine dummy_precon( g_a, g_b )
      implicit none
#include "global.fh"
      integer g_a, g_b

      return
      end


