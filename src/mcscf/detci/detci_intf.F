
#define DETCI_INTFACE
c
c  Need this macro to exclude declarations
c  in detci.fh
c

      subroutine detci_init( norb, nela, nelb, nsym, symstate,
     $                       mul, osym, iprint, eps, h, g )
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
      integer norb
      integer nela
      integer nelb
      integer nsym
      integer symstate
      integer mul(nsym,nsym)
      integer osym(norb)
      integer iprint
      double precision eps(norb)
      double precision h(*)
      double precision g(*)
c
c
c
      real tx,t2(2)
      integer nexa, nexb
      integer nstra, nstrb
      integer nekla, neklb
      integer ntij, civlen
      integer ijmap(detci_maxorb*detci_maxorb)
      integer vtab((detci_maxorb+1)*(detci_maxelec+1)*(detci_maxsy))
      integer i, j, ij
      integer memuse, memuse01, memuse02, memuse03
c
c  Reset statistics
c
      tx = util_cpusec()
      detci_init_etime = 0.d0
      detci_saa_etime = 0.d0
      detci_sbb_etime = 0.d0
      detci_sab_etime = 0.d0
      detci_sigma_etime = 0.d0
      detci_density_etime = 0.d0
      detci_sigma_calls = 0
      detci_density_calls = 0
c
c  Copy CI parameters
c
      cdetci_norb = norb
      cdetci_nela = nela
      cdetci_nelb = nelb
      cdetci_nsym = nsym
      cdetci_symstate = symstate
      do i=1,norb
        cdetci_osym(i) = osym(i)
        cdetci_eps(i) = eps(i)
      enddo
      do i=1,nsym
        do j=1,nsym
          ij = (i-1)*nsym + j
          cdetci_mul(ij) = mul(j,i)
        enddo
      enddo
c     
c  Compute CI parameters
c     
      call detci_ijmap( norb, nsym, osym, mul, ntij, ijmap )
      nexa = (norb-nela+1)*nela
      nexb = (norb-nelb+1)*nelb
      nstra = detci_binomial(norb,nela)
      nstrb = detci_binomial(norb,nelb)
      nekla = detci_binomial((norb-1),(nela-1))
      neklb = detci_binomial((norb-1),(nelb-1))
      civlen = nstra*nstrb
c
c  Default reference energy is Aufbau ordering
c  ...although maybe reset by CI-guess routine!
c
      cdetci_eref = detci_refenergy( norb, nela, nelb, eps )
c
c  Toggle spin-adaption
c  Need some way of automatically turning this on and off 
c  depending requested spin-state
c
      cdetci_spinadapt = .false.
      cdetci_squant = (nela - nelb)/2.d0
c     
c  Info print
c     
      if ((iprint.gt.0).and.(ga_nodeid().eq.0)) then
        write(6,902) norb, nsym, civlen,
     $               nela, nelb, nstra, nstrb,
     $               nexa, nexb, nekla, neklb,
     $               cdetci_symstate,
     $               cdetci_spinadapt,
     $               cdetci_squant,
     $               cdetci_eref
 902   format(10x,'Active shells:',23x,i5,/,
     $        10x,'Irreps:',30x,i5,/,
     $        10x,'CI vector length:',5x,i20,//,
     $        37x,'Alpha',6x,'Beta',/,
     $        35x,2(7('-'),3x),/,
     $        10x,'Electrons:',12x,2i10,/,
     $        10x,'Strings:',14x,2i10,/,
     $        10x,'E_ij per string:',6x,2i10,/,
     $        10x,'Strings per E_ij:',5x,2i10,//,
     $        10x,'State symmetry:',8x,i5,/,
     $        10x,'Spin adaption:',13x,l1,/,
     $        10x,'S quantum number:',3x,f10.3,/,
     $        10x,'Reference energy:',5x,f15.8)
        call util_flush(6)
      endif
c     
c  Construct arc weight tables
c     
      call detci_vatable( norb, nela, nsym, osym, mul,
     $                  vtab, cdetci_ataba )
      call detci_vatable( norb, nelb, nsym, osym, mul,
     $                  vtab, cdetci_atabb )
c     
c  Allocate and construct excitation operator table
c     
      l_detci_exa = CDETCI_INVALID
      l_detci_exb = CDETCI_INVALID
      if (nexa.gt.0) then
        if (.not.ma_push_get(MT_INT, (6*nexa*(nstra+1)), 'detci:exa',
     $                       l_detci_exa, k_detci_exa))
     $    call errquit('detci: cannot allocate',0)
        call detci_excit( norb, nela, nsym, nstra, nexa, osym,
     $                    mul, ijmap, cdetci_ataba,
     $                    int_mb(k_detci_exa) )
      endif
      if (nexb.gt.0) then
        if (.not.ma_push_get(MT_INT, (6*nexb*(nstrb+1)), 'detci:exb',
     $                       l_detci_exb, k_detci_exb))
     $    call errquit('detci: cannot allocate',0)
        call detci_excit( norb, nelb, nsym, nstrb, nexb, osym,
     $                    mul, ijmap, cdetci_atabb,
     $                    int_mb(k_detci_exb) )
      endif
c
c  Allocate ERI block
c  Copy integrals into internal block
c
      if (.not.ma_push_get(MT_DBL,ntij,'detci: h1',
     $   l_detci_h, k_detci_h))
     $   call errquit('detci: cannot allocate',0)
      if (.not.ma_push_get( MT_DBL,(ntij*ntij),'detci: eri',
     $   l_detci_g, k_detci_g))
     $   call errquit('detci: cannot allocate',0)
      do i=1,ntij
        dbl_mb(k_detci_h+i-1) = h(i)
      enddo
      do i=1,(ntij*ntij)
        dbl_mb(k_detci_g+i-1) = g(i)
      enddo
c
c  Memory for scratch space should be allocated HERE!
c
      memuse03 = 2*max(nstra,nstrb) + 3*nekla + 4*nstrb*nekla
      
      if (.not.ma_push_get(MT_DBL, max(nstra,nstrb),
     $     'detci: fscr', l_detci_fscr, k_detci_fscr))
     $     call errquit('detci: cannot allocate fscr',0)

      if (.not.ma_push_get(MT_INT, nekla,
     $     'detci: rhsscr', l_detci_rhsscr, k_detci_rhsscr))
     $     call errquit('detci: cannot allocate rhs',0)

      if (.not.ma_push_get(MT_INT, nekla,
     $     'detci: lhsscr', l_detci_lhsscr, k_detci_lhsscr))
     $     call errquit('detci: cannot allocate lhs',0)

      if (.not.ma_push_get(MT_INT, nekla,
     $     'detci: iscr', l_detci_iscr, k_detci_iscr))
     $     call errquit('detci: cannot allocate iscr',0)

      if (.not.ma_push_get(MT_DBL, (nstrb*nekla),
     $     'detci: cprime', l_detci_cprime, k_detci_cprime))
     $     call errquit('detci: cannot allocate cprime',0)

      if (.not.ma_push_get(MT_DBL, (nstrb*nekla),
     $     'detci: sprime', l_detci_sprime, k_detci_sprime))
     $     call errquit('detci: cannot allocate sprime',0)
c
c  Validate internals
c
      cdetci_valid = CDETCI_MAGIC
      detci_init_etime = util_cpusec() - tx
      return
      end






c
c  Compute internal memory required & CI-vector length
c

      integer function detci_memsiz( norb, nela, nelb, nsym,
     $                             mul, osym, vlen )
      implicit none
      integer norb
      integer nela
      integer nelb
      integer nsym
      integer mul(nsym,nsym)
      integer osym(norb)
      integer vlen
c
      integer nexa 
      integer nexb 
      integer nstra
      integer nstrb
      integer nekla
      integer neklb
      integer detci_binomial
      external detci_binomial
c
c
c
      nexa = (norb-nela+1)*nela
      nexb = (norb-nelb+1)*nelb
      nstra = detci_binomial(norb,nela)
      nstrb = detci_binomial(norb,nelb)
      nekla = detci_binomial((norb-1),(nela-1))
      neklb = detci_binomial((norb-1),(nelb-1))
      
      vlen = nstra*nstrb
      detci_memsiz = 0
      return
      end
      




      subroutine detci_free()
      implicit none
c
#include "mafdecls.fh"      
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
c
c  Check if initialized
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_free with uninitialized internals',0)
c
c Print stats
c
c$$$      write(6,901) detci_sigma_calls,
c$$$     $             detci_saa_etime,
c$$$     $             detci_sbb_etime,
c$$$     $             detci_sab_etime,
c$$$     $             detci_sigma_etime
 901  format(/,10x,'Number of sigma calls:',4x,i5,
     $       /,10x,'Elapsed time (aa):',4x,f10.2,
     $       /,23x,'(bb):',4x,f10.2,
     $       /,23x,'(ab):',4x,f10.2,
     $       /,23x,'(total)',2x,f10.2)
c
c  Free scratch memory
c
      if (.not.ma_pop_stack(l_detci_sprime))
     $   call errquit('detci: cannot pop stack',0)
      if (.not.ma_pop_stack(l_detci_cprime))
     $   call errquit('detci: cannot pop stack',0)
      if (.not.ma_pop_stack(l_detci_iscr))
     $   call errquit('detci: cannot pop stack',0)
      if (.not.ma_pop_stack(l_detci_lhsscr))
     $   call errquit('detci: cannot pop stack',0)
      if (.not.ma_pop_stack(l_detci_rhsscr))
     $   call errquit('detci: cannot pop stack',0)
      if (.not.ma_pop_stack(l_detci_fscr))
     $   call errquit('detci: cannot pop stack',0)
c
c  Free other internal memory blocks
c
      if (.not.ma_pop_stack(l_detci_g))
     $   call errquit('detci: cannot pop stack',0)
      if (.not.ma_pop_stack(l_detci_h))
     $   call errquit('detci: cannot pop stack',0)
      if (l_detci_exb.ne.CDETCI_INVALID) then
         if (.not.ma_pop_stack(l_detci_exb))
     $       call errquit('detci: cannot pop stack',0)
      endif
      if (l_detci_exa.ne.CDETCI_INVALID) then
        if (.not.ma_pop_stack(l_detci_exa))
     $       call errquit('detci: cannot pop stack',0)
      endif
c
c
c  Invalidate common block
c
      cdetci_valid = 0
c
c
      return
      end





c
c  Create a GA CI vector consistent 
c  with previously initialized parameters
c
      logical function ga_civec_create( label, g_a )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
      character*(*) label
      integer g_a
c
      integer nstra
      integer nstrb
c
c  Check if initialized
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('ga_civec_create called uninitialized',0)
c
c
c
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
      ga_civec_create = ga_create( MT_DBL, nstrb, nstra, label,
     $                             nstrb, 1, g_a )
      return
      end









      subroutine detci_guess( ngs, cfggs, cgs, g_civec )
      implicit none
c
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
c
      integer ngs
      integer cfggs(*)
      double precision cgs(ngs)
      integer g_civec
c
      integer nstra
      integer nstrb
c
c
c     
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_guess with uninitialized internals',0)
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c
c
      call detci_ciguess( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                    cdetci_nelb, nstra, nstrb, cdetci_osym,
     $                    cdetci_mul, cdetci_ataba, cdetci_atabb,
     $                    ngs, cfggs, cgs, g_civec )
c
c
c
      return
      end






c
c  Easy interface to sigma vector product
c
c
      subroutine detci_sigma( g_civec, g_sigma )
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
      integer g_civec
      integer g_sigma
c
c$$$      double precision civec(*)
c$$$      double precision sigma(*)
c
      real tx, ty, tz, t2(2)
      integer nexa, nexb
      integer nstra, nstrb
      integer nekla, neklb
      integer ntij
      integer ijmap(detci_maxorb*detci_maxorb)
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_sigma with uninitialized internals',0)
c
c
c
      detci_sigma_calls = detci_sigma_calls + 1
      tx = util_cpusec()
c     
c  CI parameters
c     
      call detci_ijmap( cdetci_norb, cdetci_nsym, cdetci_osym,
     $                cdetci_mul, ntij, ijmap )
      nexa = (cdetci_norb-cdetci_nela+1)*cdetci_nela
      nexb = (cdetci_norb-cdetci_nelb+1)*cdetci_nelb
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
      nekla = detci_binomial( (cdetci_norb-1), (cdetci_nela-1) )
      neklb = detci_binomial( (cdetci_norb-1), (cdetci_nelb-1) )
c
c  Sigma components
c
      call detci_sigmaaa( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                    cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                    nekla, neklb, cdetci_osym, cdetci_mul,
     $                    ijmap, int_mb(k_detci_exa),
     $                    int_mb(k_detci_exb),
     $                    cdetci_ataba, cdetci_atabb, ntij,
     $                    dbl_mb(k_detci_h), dbl_mb(k_detci_g),
     $                    dbl_mb(k_detci_fscr),
     $                    g_civec, g_sigma )

      ty = util_cpusec()
      detci_saa_etime = detci_saa_etime + ty - tx

      call detci_sigmabb( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                    cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                    nekla, neklb, cdetci_osym, cdetci_mul,
     $                    ijmap, int_mb(k_detci_exa),
     $                    int_mb(k_detci_exb),
     $                    cdetci_ataba, cdetci_atabb, ntij,
     $                    dbl_mb(k_detci_h), dbl_mb(k_detci_g),
     $                    dbl_mb(k_detci_fscr),
     $                    g_civec, g_sigma )

C$$$      CALL DETCI_CIVEC_PRINT( CDETCI_NORB, CDETCI_NSYM, CDETCI_NELA, CDETCI_NELB,
C$$$     $                      NSTRA, NSTRB, CDETCI_MUL, CDETCI_OSYM,
C$$$     $                      CDETCI_ATABA, CDETCI_ATABB, SIGMA, 1.D-5 )

      tz = util_cpusec()
      detci_sbb_etime = detci_sbb_etime + tz - ty

      call detci_sigmaab( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                    cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                    nekla, neklb, cdetci_osym, cdetci_mul,
     $                    ijmap, int_mb(k_detci_exa),
     $                    int_mb(k_detci_exb),
     $                    cdetci_ataba, cdetci_atabb, ntij,
     $                    dbl_mb(k_detci_g),
     $                    int_mb(k_detci_rhsscr),
     $                    int_mb(k_detci_lhsscr),
     $                    int_mb(k_detci_iscr),
     $                    dbl_mb(k_detci_fscr),
     $                    dbl_mb(k_detci_cprime),
     $                    dbl_mb(k_detci_sprime),
     $                    g_civec, g_sigma )
c
c
c$$$      CALL DETCI_CIVEC_PRINT( CDETCI_NORB, CDETCI_NSYM, CDETCI_NELA, CDETCI_NELB,
c$$$     $                      NSTRA, NSTRB, CDETCI_MUL, CDETCI_OSYM,
c$$$     $                      CDETCI_ATABA, CDETCI_ATABB, SIGMA, 1.D-5 )
c
c
c
      ty = util_cpusec()
      detci_sab_etime = detci_sab_etime + ty - tz
      detci_sigma_etime = detci_sigma_etime + ty - tx
      return
      end








      subroutine detci_ciprecon( g_civec, g_workvec )
      implicit none
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
c
      integer g_civec
      integer g_workvec
      integer nstra, nstrb
c
      double precision xx, yy
      double precision mxcnt
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_ciprecon with uninitialized common',0)
c     
c  CI parameters
c     
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c  First spin-adapt the vector
c
      if (cdetci_spinadapt) then
C        CALL DETCI_RANDOM_ERROR( G_CIVEC )       ! INTRODUCE CONTAMINATION TO DEBUG
        xx = sqrt( ga_ddot( g_civec, g_civec ) )
        call detci_spinadapt( g_civec, g_workvec )
        yy = sqrt( ga_ddot( g_civec, g_civec ) )
c$$$        WRITE(6,988) (1.d0 - YY/XX)
c$$$ 988    FORMAT('Spin contamination:',E12.3)
      endif
c
c  Apply preconditioning using orbital 
c  energies (Moller-Plesset denominators)
c
      call detci_diagscale( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                      cdetci_nelb, nstra, nstrb, cdetci_osym,
     $                      cdetci_mul, cdetci_ataba, cdetci_atabb,
     $                      cdetci_eref, cdetci_eps, g_civec )
c
c  Ensure vector is symmetry-adapted
c
      if (cdetci_nsym.gt.1) then
        call detci_symmproject( cdetci_norb, cdetci_nsym,
     $                          cdetci_nela, cdetci_nelb, nstra,
     $                          nstrb, cdetci_osym, cdetci_mul,
     $                          cdetci_ataba, cdetci_atabb,
     $                          cdetci_symstate, .true.,
     $                          mxcnt, g_civec )
      endif

      return
      end







c
c  Print CI vector coefficients larger
c  than threshold
c
      subroutine detci_print( g_civec, thresh )
      implicit none
c
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
c
      integer g_civec
C      double precision civec(*)
      double precision thresh
      integer nstra, nstrb
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_print with uninitialized common',0)
c     
c  CI parameters
c     
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c
c
      call detci_civec_print( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                        cdetci_nelb, nstra, nstrb, cdetci_mul,
     $                        cdetci_osym, cdetci_ataba, cdetci_atabb,
     $                        g_civec, thresh )
c
c
c
      return
      end








c
c  Generate 1- and 2-particle density matrices
c  for given CI vector (both RHS and LHS)
c
      subroutine detci_density( g_civec, map, onepdm, twopdm )
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
c
      integer g_civec
      integer map(*)
      double precision onepdm(*)
      double precision twopdm(*)
c
      integer nexa,nexb
      integer nstra,nstrb
      integer ntij
      integer ijmap(detci_maxorb*detci_maxorb)
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_density with uninitialized common',0)
c     
c  CI parameters
c     
      call detci_ijmap( cdetci_norb, cdetci_nsym, cdetci_osym,
     $                  cdetci_mul, ntij, ijmap )
      nexa = (cdetci_norb-cdetci_nela+1)*cdetci_nela
      nexb = (cdetci_norb-cdetci_nelb+1)*cdetci_nelb
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c
c
      call detci_onepdm( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                   cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                   cdetci_osym, cdetci_mul, int_mb(k_detci_exa),
     $                   int_mb(k_detci_exb), map, g_civec, onepdm )
c
c
c
      call detci_twopdm( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                   cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                   cdetci_osym, cdetci_mul, int_mb(k_detci_exa),
     $                   int_mb(k_detci_exb), map, g_civec, twopdm )
c
c
c
c$$$      call detci_etrace( cdetci_norb, ntij, ijmap,
c$$$     $                   dbl_mb(k_detci_h), dbl_mb(k_detci_g),
c$$$     $                   onepdm, twopdm )
c
c
c
      return
      end





c
c  Spin-adaption routine
c  Project off contaminants by Lowdin projection operator
c  S quantum number is determined by the highest component, M,
c  ((nela - nelb)/2)
c

      subroutine detci_spinadapt( g_civec, g_wvec )
      implicit none
c
#include "mafdecls.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
c$$$      double precision civec(*)
c$$$      double precision wvec(*)
      integer g_civec
      integer g_wvec
c
      integer nexa, nexb
      integer nstra, nstrb
      integer ntij
      integer ijmap(detci_maxorb*detci_maxorb)
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_spin with uninitialized internals',0)
c
c  CI parameters
c     
      call detci_ijmap( cdetci_norb, cdetci_nsym, cdetci_osym,
     $                cdetci_mul, ntij, ijmap )
      nexa = (cdetci_norb-cdetci_nela+1)*cdetci_nela
      nexb = (cdetci_norb-cdetci_nelb+1)*cdetci_nelb
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c  
c
      call detci_spadpt( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                   cdetci_nelb, nstra, nstrb, cdetci_mul,
     $                   cdetci_osym, cdetci_ataba, cdetci_atabb,
     $                   nexa, nexb, int_mb(k_detci_exa),
     $                   int_mb(k_detci_exb), g_civec, g_wvec )
c
c
      return
      end

      








      subroutine detci_symmproj( g_civec, mxcnt, oscreen )
      implicit none
#include "global.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
c
      integer g_civec
      logical oscreen
      double precision mxcnt
      integer nstra, nstrb
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_symmproj with uninitialized common',0)
c     
c  CI parameters
c     
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c  Symmetry adapt
c
      mxcnt = 0.d0
      if (cdetci_nsym.gt.1) then
        call detci_symmproject( cdetci_norb, cdetci_nsym,
     $                          cdetci_nela, cdetci_nelb, nstra,
     $                          nstrb, cdetci_osym, cdetci_mul,
     $                          cdetci_ataba, cdetci_atabb,
     $                          cdetci_symstate, oscreen,
     $                          mxcnt, g_civec )
      endif

      return
      end

