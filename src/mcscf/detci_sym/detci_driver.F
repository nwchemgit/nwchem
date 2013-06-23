      program detci_driver
*
* $Id$
*
      implicit none
#include "detciP.fh"
#include "detci.fh"
#include "util.fh"
      integer norb, nela, nelb, nsym
      integer ntij
      integer civlen
      integer maxiter
      integer osym(detci_maxorb)
      integer ijmap(detci_maxorb*detci_maxorb)
      double precision escf, edetci, potnuc
      double precision eps(detci_maxorb)
      double precision h(detci_maxtij)
      double precision g(detci_maxtij*detci_maxtij)
      double precision onepdm(detci_maxorb*detci_maxorb)
      double precision twopdm(detci_maxorb*detci_maxorb*
     $                        detci_maxorb*detci_maxorb)
      double precision tol
      integer iprint
      integer i, j
      integer iioff, ngs, cfggs(2*detci_maxelec*100)
      double precision cgs(100)
      double precision xx
      data iprint/5/
      data maxiter/12/
c     
c     
c     
      integer nwork
      integer stack_req
      parameter(nwork=5)
      double precision civec
      double precision work
      common civec(detci_maxstr*detci_maxstr),
     $       work(detci_maxstr*detci_maxstr,nwork)
c
c
c     
      DOUBLE PRECISION DRAND
      EXTERNAL DRAND
      external detci_sigma, detci_ciprecon
c
c
c
      REAL TX, T2(2)
      XX = DRAND(234)
c
c read CI parameters
c     
      read(5,*) norb,nela,nelb,nsym
      read(5,*) tol
      if ((norb.gt.detci_maxorb).or.(nela.gt.detci_maxelec).or.
     $    (nelb.gt.detci_maxelec)) then
        write(6,922)
 922    format(/,'Number of orbitals/electron too large - redimension')
        stop 777
      endif
      read(5,*) (osym(i),i=1,norb)
c
c Read in orbital energies
c
 145  read(5,*) i,xx
      if (i.gt.0) then
        eps(i) = xx
        goto 145
      endif
c
c Read in starting Ground-state config 
c
      ngs = 1
      iioff = 0
 133  read(5,*) cgs(ngs),(cfggs(iioff+i),i=1,(nela+nelb))
      iioff = iioff + nela + nelb
      if (cgs(ngs).ne.0.d0) then
        ngs = ngs + 1
        goto 133
      endif
      ngs = ngs - 1
c     
c  Read integral stuff
c     
      call detci_ijmap( norb, nsym, osym, ntij, ijmap )
      call detci_readints( norb, ntij, ijmap, potnuc, escf, eps, h, g )
c
c  Initialize DETCI internals
c
      stack_req = detci_memsiz( norb, nela, nelb, nsym, 
     $                        osym, civlen )
      call detci_init( norb, nela, nelb, nsym, osym,
     $               iprint, eps, h, g )
c
c  Call guess routine
c
      call detci_guess( ngs, cfggs, cgs, civec )
c
c  Spin adapt
c
      call detci_spinadapt( civec, work(1,1) )
c
c  Test density, preconditioner, sigma
c
C$$$      CALL DETCI_DENSITY( CIVEC, ONEPDM, TWOPDM )
C$$$      CALL DETCI_ETRACE(NORB,NTIJ,IJMAP,H,G,ONEPDM,TWOPDM)
C$$$      CALL DETCI_CIPRECON( CIVEC )
c
c Time one sigma product
c     
      TX = util_cpusec()
      CALL DETCI_SIGMA( CIVEC, WORK(1,1) )
C      CALL DETCI_PRINT( WORK(1,1), 1.0D-3 )
      TX = util_cpusec() - TX
      WRITE(6,933) TX
 933  FORMAT(/,'Time for single CI vector product:',f12.3,/)
c
c  Call Davidson routine
c
C$$$      CALL DETCI_DVDDIAG0( CIVLEN, DETCI_SIGMA, DETCI_CIPRECON,
C$$$     $                  NWORK, WORK, TOL, CIVEC, EDETCI )
C$$$      CALL DETCI_GUESS( GS_ELVA, GS_ELVB, CIVEC )
c
      call detci_dvddiag( civlen, detci_sigma, detci_ciprecon,
     $                  maxiter, tol, civec, work(1,1), work(1,2),
     $                  work(1,3), work(1,4), work(1,5), edetci)
      write(6,'(/,10x,a,13x,f20.10)') '2e energy',edetci
      write(6,'(10x,a,10x,f20.10)') 'Total energy',(potnuc+edetci)
c
c  Write out CI vector
c
      write(6,'(/,10x,a20,/)') 'Solution CI vector'
      call detci_print( civec, 5.0d-2 )
c
c  Test spin
c
C      call detci_spinadapt( civec, work(1,1) )
C      call detci_print( work(1,1), 1.0d-3 )
c
c  Density matrix construction
c
      call detci_density( civec, onepdm, twopdm )
c
c  Free DETCI internals
c
      call detci_free()
c
c  Energy trace 
c  !! this can be called after detci_free !!
c
      call detci_etrace( norb, ntij, ijmap, h, g, onepdm, twopdm )
c
c
c
 1000 CONTINUE
      end

