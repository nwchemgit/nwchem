c
c   MCSCF Secular solution front-end
c
      subroutine mcscf_cisolve( rtdb, geom, basis, nbf, nclosed,
     $                          nact, nsym, orbsym, eps, g_h, g_coul,
     $                          dm1, dm2, oprint, ociprint )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "detciP.fh"
#include "detci.fh"
#include "mcscfprof.fh"
      integer rtdb
      integer geom
      integer basis
      integer nbf
      integer nclosed
      integer nact
      integer nsym
      integer orbsym(nbf)
      integer g_h
      integer g_coul
      double precision eps(nbf)
      double precision dm1(nact,nact)
      double precision dm2(nact,nact,nact,nact)
      logical oprint
      logical ociprint
c
c
c
#ifdef OLD_DAVIDSON
      integer mxworkvec
      parameter(mxworkvec=24)
#else
      integer mxworkvec
      parameter(mxworkvec=5)
#endif      
      integer nela, nelb, nactel, multip, symstate
      integer actsym(detci_maxorb), actmap(detci_maxorb)
      integer actrmap(detci_maxorb)
      integer guess_cfg(12)
      double precision guess_coeff(6)
      double precision ciener, citol
      double precision acteps(detci_maxorb)
      double precision mxcnt
      integer maxiter
      integer mul(detci_maxsy*detci_maxsy)
      integer ijmap(detci_maxorb*detci_maxorb)
      integer i, iprint, hasprint
      integer ntij, n2
      integer l_h1, k_h1
      integer l_eri, k_eri
      integer g_civec, g_wvec(mxworkvec)
      INTEGER J,II
      character*10 vstr
      character*10 civec_in
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
      data citol/1.e-7/
      data maxiter/30/
      data civec_in/'civec'/
c
c
      iprint = 1
      if (hasprint.ne.0) iprint = 0
      hasprint = hasprint + 1
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
     $    call errquit('detci: incompatible elec and spin',0)
c
c   Retrieve symmetry stuff
c
      if (.not.rtdb_get(rtdb,'mcscf:statesymmetry',MT_INT,1,symstate))
     $   symstate = 1
      call getsymm_multab( nsym, mul )
c
c   Remap active indices into symmetry blocks
c
      ii = 0
      do j=1,nsym
        do i=1,nact
          if (orbsym(nclosed+i).eq.j) then
            ii = ii + 1
            actmap(ii) = i
            actrmap(i) = ii
            actsym(ii) = orbsym(nclosed+i)
            acteps(ii) = eps(nclosed+i)
          endif
        enddo
      enddo
c
c$$$      MAXITER = 1                      ! The following removes 
c$$$      SYMSTATE = 1                     ! symmetry mapping
c$$$      NSYM = 1
c$$$      DO I=1,NACT
c$$$        ACTSYM(I) = 1
c$$$        ACTMAP(I) = I
c$$$        ACTRMAP(I) = I
c$$$        ACTEPS(I) = EPS(NCLOSED+I)
c$$$      ENDDO
c
c
c$$$      if (ga_nodeid().eq.0) then
c$$$        write(6,771) 'Orb map',(actmap(i),i=1,nact)
c$$$        write(6,771) 'Orb rev map',(actrmap(i),i=1,nact)
c$$$        write(6,771) 'Orb symm',(actsym(i),i=1,nact)
c$$$ 771    format(a20,5x,10i3)
c$$$      endif
c
c
      call detci_ijmap( nact, nsym, actsym, mul, ntij, ijmap )
c
c   Allocate local space for 1e and 2e MO integrals
c   Push global MO integrals into local space
c
      n2 = nact*nact
      if (.not.ma_push_get(MT_DBL,n2,'1e hamil',l_h1,k_h1))
     $  call errquit('mcscf_cisolve: cannot allocate',0)
      if (.not.ma_push_get(MT_DBL,(n2*n2),'eri',l_eri,k_eri))
     $  call errquit('mcscf_cisolve: cannot allocate',0)
      call dfill(n2,0.d0,dbl_mb(k_h1),1)
      call dfill((n2*n2),0.d0,dbl_mb(k_eri),1)
      call mcscf_getmolocal( nclosed, nact, nbf, actrmap,
     $                       g_h, g_coul,
     $                       ntij, ijmap, dbl_mb(k_h1),
     $                       dbl_mb(k_eri) )
c
c   Initialize DetCI internals
c
      call detci_init( nact, nela, nelb, nsym, symstate,
     $                 mul, actsym, iprint, acteps,
     $                 dbl_mb(k_h1), dbl_mb(k_eri) )
c
c   Create CI-vectors
c
      if (.not.ga_civec_create('CI vector',g_civec))
     $   call errquit('detci: cannot create CI vector',0)
      do i=1,mxworkvec
        vstr = 'CI work ' // char(ichar('0')+i)
        if (.not.ga_civec_create(vstr,g_wvec(i)))
     $    call errquit('detci: cannot create CI vector',0)
        call ga_zero(g_wvec(i))
      enddo
c
c   Guess CI vector for first time thru
c   Can manually set the guess config. in the database
c   otherwise generate guess via Aufbau rules with correct
c   symmetry
c
      inquire(file=civec_in,exist=civec_saved)
      if ((hasprint.le.1).or.(.not.(civec_saved))) then
        guess_coeff(1) = 1.d0
        if (.not.rtdb_get(rtdb,'mcscf:guess',MT_INT,
     $    (nela+nelb),guess_cfg))
     $    call detci_guesscfg( nact, nsym, nela, nelb, actsym, mul,
     $                         acteps, symstate, guess_cfg )
        call detci_guess( 1, guess_cfg, guess_coeff, g_civec )
      else
        if (.not.(file_read_ga(civec_in,g_civec)))
     $    call errquit('mcscf_ci: cannot read CI vector',0)
      endif
c
c   Davidson diag
c
      if ((omcscfprof).and.(ociprof)) call pstat_on(ps_cidiag)
#ifdef OLD_DAVIDSON
      call detci_dvddiag0( detci_sigma, detci_ciprecon, mxworkvec,
     $                     1.d-12, ciener, g_wvec, g_civec, oprint )
#else
      call detci_dvddiag( detci_sigma, detci_ciprecon, maxiter,
     $                    citol, oprint, g_civec, g_wvec(1), g_wvec(2),
     $                    g_wvec(3), g_wvec(4), g_wvec(5), ciener )
#endif
      if ((omcscfprof).and.(ociprof)) call pstat_off(ps_cidiag)
c
c   Print out converged ground state CI vector
c
      if (ociprint) call detci_print( g_civec, 2.d-2 )
c
c   Check symmetry contamination of soln vector
c
      if (nsym.gt.1) then
        call detci_symmproj( g_civec, mxcnt, .false. )
        if (mxcnt.gt.1.d-14) then
          if (ga_nodeid().eq.0) write(6,881) mxcnt
 881      format('Warning - CI symmetry contamination:',e12.3)
        else if (ociprint) then
          if (ga_nodeid().eq.0) write(6,882) mxcnt
 882      format(//,10x,'CI Symmetry contamination:',e12.3)
        endif
      endif
c
c   Compute MO densities
c
      if ((omcscfprof).and.(ociprof)) call pstat_on(ps_cidens)
      call detci_density( g_civec, actmap, dm1, dm2 )
      if ((omcscfprof).and.(ociprof)) call pstat_off(ps_cidens)
c
c   Save CI vector for restart
c
      if (.not. file_write_ga(civec_in,g_civec))
     $  call errquit('mcscf_ci: error in CI vector save',0)
c
c   Release CI vectors
c
      if (.not.ga_destroy(g_civec))
     $   call errquit('mcscf_cisolve: cannot destroy CI vector',0)
      do i=1,mxworkvec
        if (.not.ga_destroy(g_wvec(i)))
     $    call errquit('mcscf_cisolve: cannot destroy CI vector',0)
      enddo
c
c   Release DetCI internals
c
      call detci_free()
c
c   Clean up
c
      if (.not.ma_pop_stack(l_eri))
     $  call errquit('mcscf_cisolve: cannot pop stack',0)
      if (.not.ma_pop_stack(l_h1))
     $  call errquit('mcscf_cisolve: cannot pop stack',0)
c
c
c
      return
      end






      subroutine mcscf_getmolocal( nclosed, nact, nbf, rmap,
     $                             g_h, g_coul,
     $                             ntij, ijmap, h1, eri )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nclosed
      integer nact
      integer nbf
      integer rmap(nact)
      integer g_h
      integer g_coul
      integer ntij
      integer ijmap(nact,nact)
      double precision h1(ntij)
      double precision eri(ntij,ntij)
c
      integer i, j, k, l, gij, ij, kl, ii, jj
      integer coff, kk, klo, khi, n2
      integer l_t,k_t
      double precision xx
c
c
c
      if (.not.ma_push_get(MT_DBL,(nact*nact),'t',l_t,k_t))
     $  call errquit('mcscf_getmolocal: cannot allocate',0)
      call dfill(ntij,0.d0,h1,1)
      call dfill((ntij*ntij),0.d0,eri,1)

      coff = nclosed*nbf + nclosed
      do i=1,nact
        do j=1,i
          gij = ((i-1)*i)/2 + j
          ij = ijmap(rmap(i),rmap(j))
          do k=1,nact
            kk = (k-1)*nact
            klo = coff + (k-1)*nbf + 1
            khi = coff + (k-1)*nbf + nact
            call ga_get(g_coul,klo,khi,gij,gij,dbl_mb(k_t+kk),1)
          enddo
          do k=1,nact
            kk = (k-1)*nact
            do l=1,k
              kl = ijmap(rmap(k),rmap(l))
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
     $  call errquit('mcscf_getmolocal: cannot pop stack',0)
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


