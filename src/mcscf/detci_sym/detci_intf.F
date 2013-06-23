#define DETCI_INTFACE
c
c  Need this macro to exclude declarations
c  in detci.fh
c

      subroutine detci_init( norb, nela, nelb, nsym, symstate,
     $                       osym, iprint, eps, h, g )
*
* $Id$
*
      implicit none
#include "errquit.fh"
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
      integer osym(norb)
      integer iprint
      double precision eps(norb)
      double precision h(*)
      double precision g(*)
c
c
c
      double precision tx
      integer nexa, nexb
      integer nstra, nstrb
      integer nekla, neklb
      integer ntij, civlen
      integer ijmap(detci_maxorb*detci_maxorb)
      integer vtab((detci_maxorb+1)*(detci_maxelec+1)*(detci_maxsy))
      integer i, j, ii
c$$$      integer ij
      integer memuse03
      integer mseg
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
      detci_spinadp_etime = 0.d0
      detci_symmadp_etime = 0.d0
      detci_aaff_etime = 0.d0
      detci_aadot_etime = 0.d0
      detci_aagop_etime = 0.d0
      detci_abstr_etime = 0.d0
      detci_abgath_etime = 0.d0
      detci_abdotab_etime = 0.d0
      detci_abscat_etime = 0.d0
      detci_absync_etime = 0.d0
      detci_density_onepdm = 0.0d0
      detci_density_twopdm = 0.0d0
      detci_density_twopdmab = 0.0d0
      cdetci_profprint = .true. 
c
c  Copy CI parameters
c
      cdetci_norb = norb
      cdetci_nela = nela
      cdetci_nelb = nelb
      cdetci_nsym = nsym
      cdetci_symstate = symstate
c
c   Reorder orbital indices by symmetry
c
      ii = 0
      do j=1,nsym
        do i=1,norb
          if (osym(i).eq.j) then
            ii = ii + 1
            cdetci_ixmap(ii) = i
            cdetci_irmap(i) = ii
            cdetci_osym(ii) = osym(i)
            cdetci_eps(ii) = eps(i)
          endif
        enddo
      enddo
c
c$$$      DO I=1,NORB
c$$$        WRITE(6,'(I5,A5,I5,5X,F10.4)') I,' ->',CDETCI_IRMAP(I),
c$$$     $                                 CDETCI_EPS(CDETCI_IRMAP(I))
c$$$      ENDDO
c
c  Compute CI parameters
c     
      call detci_ijmap( norb, nsym, cdetci_osym, ntij, ijmap )
      nexa = (norb-nela+1)*nela
      nexb = (norb-nelb+1)*nelb
      nstra = detci_binomial(norb,nela)
      nstrb = detci_binomial(norb,nelb)
      nekla = detci_binomial((norb-1),(nela-1))
      neklb = detci_binomial((norb-1),(nelb-1))
      civlen = nstra*nstrb
      mseg = (10*ga_nnodes())/max(min(nexa,nexb),1)
      mseg = min(mseg, 20)
      MSEG = 1                                       ! Delete for better parallel effeciency
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
      cdetci_spinadapt = .true.
      cdetci_squant = (nela - nelb)/2.d0
c     
c  Info print
c     
      if ((iprint.gt.0).and.(ga_nodeid().eq.0)) then
        write(6,902) norb, nsym, civlen,
     $               cdetci_symstate,
     $               cdetci_spinadapt,
     $               cdetci_squant,
     $               nela, nelb, nstra, nstrb,
     $               nexa, nexb, nekla, neklb
 902   format(10x,'Active shells:',23x,i5,/,
     $        10x,'Irreps:',30x,i5,/,
     $        10x,'CI vector length:',5x,i20,/,
     $        10x,'State symmetry:',25x,i2,/,
     $        10x,'Spin adaption:',27x,l1,/,
     $        10x,'S quantum number:',15x,f10.3,//,
     $        37x,'Alpha',6x,'Beta',/,
     $        35x,2(7('-'),3x),/,
     $        10x,'Electrons:',12x,2i10,/,
     $        10x,'Strings:',14x,2i10,/,
     $        10x,'E_ij per string:',6x,2i10,/,
     $        10x,'Strings per E_ij:',5x,2i10,/)
        call util_flush(6)
      endif
c     
c  Construct arc weight tables
c     
      call detci_vatable( norb, nela, nsym, cdetci_osym, 
     $                    vtab, cdetci_ataba )
      call detci_vatable( norb, nelb, nsym, cdetci_osym,
     $                    vtab, cdetci_atabb )
c     
c  Allocate and construct excitation operator table
c     
      l_detci_exa = CDETCI_INVALID
      l_detci_exb = CDETCI_INVALID
      if (nexa.gt.0) then
        if (.not.ma_push_get(MT_INT, (6*nexa*(nstra+1)), 'detci:exa',
     $                       l_detci_exa, k_detci_exa))
     $    call errquit('detci: cannot allocate',0, MA_ERR)
        call detci_excit( norb, nela, nsym, nstra, nexa, cdetci_osym,
     $                    ijmap, cdetci_ataba, int_mb(k_detci_exa) )
      endif
      if (nexb.gt.0) then
        if (.not.ma_push_get(MT_INT, (6*nexb*(nstrb+1)), 'detci:exb',
     $                       l_detci_exb, k_detci_exb))
     $    call errquit('detci: cannot allocate',0, MA_ERR)
        call detci_excit( norb, nelb, nsym, nstrb, nexb, cdetci_osym,
     $                    ijmap, cdetci_atabb, int_mb(k_detci_exb) )
      endif
c$$$      IF (GA_NODEID().EQ.0) THEN
c$$$        IJ = 6*(NEXA*(NSTRA+1) + NEXB*(NSTRB+1))
c$$$        I = MA_INQUIRE_AVAIL( MT_DBL )
c$$$        WRITE(6,913) IJ, I
c$$$ 913    FORMAT('ALLOCATED AND GENERATED EXCITATION TABLE',/,
c$$$     $          'Requested ',I10,' Free',I10,' words')
c$$$        CALL UTIL_FLUSH(6)
c$$$      ENDIF
c$$$      CALL GA_SYNC()
c
c  Allocate ERI block
c  Copy integrals into internal block
c
      if (.not.ma_push_get(MT_DBL,ntij,'detci: h1',
     $   l_detci_h, k_detci_h))
     $   call errquit('detci: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get( MT_DBL,(ntij*ntij),'detci: eri',
     $   l_detci_g, k_detci_g))
     $   call errquit('detci: cannot allocate',0, MA_ERR)
      call detci_moint_copy( norb, ntij, cdetci_ixmap, ijmap,
     $                       cdetci_osym, h, g,
     $                       dbl_mb(k_detci_h), dbl_mb(k_detci_g) )
c
c  Memory for scratch space should be allocated HERE!
c
      memuse03 = 2*max(nstra,nstrb)*mseg + 3*nekla + 4*nstrb*nekla
      if (.not.ma_push_get(MT_DBL, (max(max(nstra,nstrb)*mseg,nexb)),
     $     'detci: fscr', l_detci_fscr, k_detci_fscr))
     $     call errquit('detci: cannot allocate fscr',0, MA_ERR)

      if (.not.ma_push_get(MT_INT, nekla,
     $     'detci: rhsscr', l_detci_rhsscr, k_detci_rhsscr))
     $     call errquit('detci: cannot allocate rhs',0, MA_ERR)

      if (.not.ma_push_get(MT_INT, max(nekla,neklb,nexb,nexa),
     $     'detci: lhsscr', l_detci_lhsscr, k_detci_lhsscr))
     $     call errquit('detci: cannot allocate lhs',0, MA_ERR)

      if (.not.ma_push_get(MT_INT, nekla,
     $     'detci: iscr', l_detci_iscr, k_detci_iscr))
     $     call errquit('detci: cannot allocate iscr',0, MA_ERR)

      if (.not.ma_push_get(MT_DBL, (nstrb*nekla),
     $     'detci: cprime', l_detci_cprime, k_detci_cprime))
     $     call errquit('detci: cannot allocate cprime',0, MA_ERR)

      if (.not.ma_push_get(MT_DBL, (nstrb*nekla),
     $     'detci: sprime', l_detci_sprime, k_detci_sprime))
     $     call errquit('detci: cannot allocate sprime',0, MA_ERR)
c
c$$$      IF (GA_NODEID().EQ.0) THEN
c$$$        IJ = MAX(NSTRA,NSTRB) + 3*NEKLA + 2*(NSTRB*NEKLA)
c$$$        I = MA_INQUIRE_AVAIL( MT_DBL )
c$$$        WRITE(6,914) IJ, I
c$$$ 914    FORMAT('ALLOCATED SCRATCH ARRAYS',/,
c$$$     $         'Requested ',I10,' Free',I10,' words')
c$$$        CALL UTIL_FLUSH(6)
c$$$      ENDIF
c$$$      CALL GA_SYNC()
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
     $                               osym, vlen )
      implicit none
      integer norb
      integer nela
      integer nelb
      integer nsym
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
#include "errquit.fh"
c
#include "mafdecls.fh"      
#include "global.fh"
#include "msgids.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
c
c  Check if initialized
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_free with uninitialized internals',0,
     &       INPUT_ERR)
c
c Print stats
c
      if (cdetci_profprint) then
        if (ga_nodeid().eq.0) write(6,901) detci_sigma_calls,
     $                                     detci_saa_etime,
     $                                     detci_sbb_etime,
     $                                     detci_sab_etime,
     $                                     detci_sigma_etime,
     $                                     detci_aaff_etime,
     $                                     detci_aagop_etime,
     $                                     detci_aadot_etime,
     $                                     detci_abstr_etime,
     $                                     detci_abgath_etime,
     $                                     detci_abdotab_etime,
     $                                     detci_abscat_etime,
     $                                     detci_absync_etime,
     $                                     detci_density_etime,
     $                                     detci_density_onepdm,
     $                                     detci_density_twopdm,
     $                                     detci_density_twopdmab,
     $                                     detci_spinadp_etime,
     $                                     detci_symmadp_etime
 901    format(/,10x,'Number of sigma calls:',4x,i5,
     $         /,23x,'o',5('<'),' (aa):',7x,f10.2,
     $         /,23x,'o',5('<'),' (bb):',7x,f10.2,
     $         /,23x,'o',5('<'),' (ab):',7x,f10.2,
     $         /,23x,'o',5('<'),' (total)',5x,f10.2,
     $         /,23x,'o',5('<'),' (aa) ff',5x,f10.2,
     $         /,23x,'o',5('<'),' (aa) gop',4x,f10.2,
     $         /,23x,'o',5('<'),' (aa) dot',4x,f10.2,
     $         /,23x,'o',5('<'),' (ab) str',4x,f10.2,
     $         /,23x,'o',5('<'),' (ab) gath',3x,f10.2,
     $         /,23x,'o',5('<'),' (ab) dotab',2x,f10.2,
     $         /,23x,'o',5('<'),' (ab) scat',3x,f10.2,
     $         /,23x,'o',5('<'),' (ab) sync',3x,f10.2,
     $         /,23x,'o',5('<'),' Density',5x,f10.2,
     $         /,23x,'o',5('<'),' Density one',1x,f10.2,
     $         /,23x,'o',5('<'),' Density two',1x,f10.2,
     $         /,23x,'o',5('<'),' Density ab',2x,f10.2,
     $         /,23x,'o',5('<'),' Spin adapt',2x,f10.2,
     $         /,23x,'o',5('<'),' Symm adapt',2x,f10.2)

        call ga_dgop(msg_detci_maxsync, detci_absync_etime,
     $               1, 'max' )
        if (ga_nodeid().eq.0) write(6,933) detci_absync_etime
 933    format(/,23x,'o',5('<'),' (ab) max sync:',f10.2)
      endif
c
c  Free scratch memory
c
      if (.not.ma_pop_stack(l_detci_sprime))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_detci_cprime))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_detci_iscr))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_detci_lhsscr))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_detci_rhsscr))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_detci_fscr))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
c
c  Free other internal memory blocks
c
      if (.not.ma_pop_stack(l_detci_g))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_detci_h))
     $   call errquit('detci: cannot pop stack',0, MA_ERR)
      if (l_detci_exb.ne.CDETCI_INVALID) then
         if (.not.ma_pop_stack(l_detci_exb))
     $       call errquit('detci: cannot pop stack',0, MA_ERR)
      endif
      if (l_detci_exa.ne.CDETCI_INVALID) then
        if (.not.ma_pop_stack(l_detci_exa))
     $       call errquit('detci: cannot pop stack',0, MA_ERR)
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
#include "errquit.fh"
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
c  Hardwired for now for 1024 processors
c
      integer map1(1024)
c
c  Check if initialized
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('ga_civec_create called uninitialized',0,
     &       INPUT_ERR)
c
c  Distribute both string_b and string_a such that different
c  symmetry blocks in string_b are not split between processors.
c  String_a will be split over nproc/nsym tasks
c
      nblock   = 0
      detstore = 0
      ntasks = nproc/nsym
      if (mod(nproc,nsym).ne.0) ntasks=ntasks+1
      do isym=1,nsym
         jsym=MULX(isym,ksym)
         nblock=nblock+1
         map1(nblock)=detstore+1
         maxa=cdetci_ataba(norb,nela,jsym)-cdetci_ataba(0,0,jsym)+1
         maxb=cdetci_atabb(norb,nelb,isym)-cdetci_atabb(0,0,isym)+1
         do i=2,ntasks
            nblock=nblock+1
            map1(nblock)=detstore+(i-1)*(maxa/ntasks)+1
         enddo
         detstore=detstore+maxa*maxb
      enddo
c
c  Generate GA
c
      ga_civec_create = ga_create_irreg( MT_DBL, detstore, 1, label,
     $                             map1,nblock,1,1, g_a )
c
      return
      end






c
c  Takes a list of configurations with coefficients 
c  and generates an initial guess CI vector that
c  is both spin- and symmetry-adapted
c  Check to ensure initial guess has non-zero
c  component in requisite spin and symmetry state
c
      subroutine detci_guess( ngs, cfggs, cgs, g_civec, g_workvec )
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
      integer ngs
      integer cfggs(*)
      double precision cgs(ngs)
      integer g_civec
      integer g_workvec
c
      integer nstra
      integer nstrb
      integer rcfg(detci_maxelec_tot*detci_maxguess_cfg)
      integer i, j, jj
      double precision mxcnt, xx, tx
c
c
c     
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_guess with uninitialized internals',0,
     &       INPUT_ERR)
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c  Reorder guess indices to internal convention
c
      do j=1,ngs
        jj = (j-1)*(cdetci_nela+cdetci_nelb)
        do i=1,cdetci_nela+cdetci_nelb
          rcfg(jj+i) = cdetci_irmap(cfggs(jj+i))
        enddo
        call util_isort(cdetci_nela,rcfg(jj+1))
        call util_isort(cdetci_nelb,rcfg(jj+cdetci_nela+1))
      enddo
c
c  Generate normalized guess CI vector
c
      call detci_ciguess( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                    cdetci_nelb, nstra, nstrb, cdetci_osym,
     $                    cdetci_ataba, cdetci_atabb,
     $                    ngs, rcfg, cgs, g_civec )
c
c  Spin adapt
c
      if (cdetci_spinadapt) then
        tx = util_cpusec()
        call detci_spinadapt( g_civec, g_workvec )
        detci_spinadp_etime = detci_spinadp_etime + util_cpusec() - tx
      endif
c
c  Symmetry adapt
c
      if (cdetci_nsym.gt.1) then
        tx = util_cpusec()
        call detci_symmproject( cdetci_norb, cdetci_nsym,
     $                          cdetci_nela, cdetci_nelb, nstra,
     $                          nstrb, cdetci_osym, 
     $                          cdetci_ataba, cdetci_atabb,
     $                          cdetci_symstate, .true.,
     $                          mxcnt, g_civec )
        detci_symmadp_etime = detci_symmadp_etime + util_cpusec() - tx
      endif
c
c  Check for zero projected vector
c
      xx = ga_ddot(g_civec,g_civec)
      if (xx.lt.0.0001d0)
     $  call errquit('detci: initial guess has wrong spin/symmetry',0,
     &       INPUT_ERR)
      return
      end






c
c  Easy interface to sigma vector product
c
c
      subroutine detci_sigma( g_civec, g_sigma )
      implicit none
#include "errquit.fh"
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
      double precision t0, tx
      integer nexa, nexb
      integer nstra, nstrb
      integer nekla, neklb
      integer ntij
      integer ijmap(detci_maxorb*detci_maxorb)
      integer g_ct, g_st
      integer mseg
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_sigma with uninitialized internals',0,
     &       INPUT_ERR)
c
c
c
      detci_sigma_calls = detci_sigma_calls + 1
      t0 = util_cpusec()
c     
c  CI parameters
c     
      call detci_ijmap( cdetci_norb, cdetci_nsym, cdetci_osym,
     $                  ntij, ijmap )
      nexa = (cdetci_norb-cdetci_nela+1)*cdetci_nela
      nexb = (cdetci_norb-cdetci_nelb+1)*cdetci_nelb
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
      nekla = detci_binomial( (cdetci_norb-1), (cdetci_nela-1) )
      neklb = detci_binomial( (cdetci_norb-1), (cdetci_nelb-1) )
      mseg = (10*ga_nnodes())/max(min(nexa,nexb),1)
      mseg = min(mseg, 20)
      MSEG = 1                                 ! Delete for better parallel effeciency
c                                              ! mseg = number of simultaneous F string contructed
c
c  Sigma components
c
c   beta-beta using transpose vectors
c
      tx = util_cpusec()
*ga:1:0
      if (.not.ga_create( MT_DBL, nstra, nstrb, 'transp CI',
     $                    nstra, 0, g_ct ))
     $   call errquit('detci_sigma: cannot allocate transp ci',0,
     &       GA_ERR)
*ga:1:0
      if (.not.ga_create( MT_DBL, nstra, nstrb, 'transp sig',
     $                    nstra, 0, g_st ))
     $   call errquit('detci_sigma: cannot allocate transp sig',0,
     &       GA_ERR)
      call ga_transpose( g_civec, g_ct )
      call ga_zero(g_st)

      call detci_sigmaaa( cdetci_norb, cdetci_nsym, cdetci_nelb,
     $                    cdetci_nela, nstrb, nstra, nexb, nexa,
     $                    neklb, nekla, cdetci_osym, 
     $                    ijmap, int_mb(k_detci_exb),
     $                    int_mb(k_detci_exa),
     $                    cdetci_atabb, cdetci_ataba, ntij, mseg,
     $                    dbl_mb(k_detci_h), dbl_mb(k_detci_g),
     $                    dbl_mb(k_detci_fscr),
     $                    g_ct, g_st )

      call ga_transpose( g_st, g_sigma )
      if (.not.ga_destroy(g_ct))
     $  call errquit('detci_sigma: cannot destroy transp ci',0,
     &       GA_ERR)
      if (.not.ga_destroy(g_st))
     $  call errquit('detci_sigma: cannot destroy transp sig',0,
     &       GA_ERR)
      detci_sbb_etime = detci_sbb_etime + util_cpusec() - tx
c
c   alpha-alpha
c
      tx = util_cpusec()
      call detci_sigmaaa( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                    cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                    nekla, neklb, cdetci_osym, 
     $                    ijmap, int_mb(k_detci_exa),
     $                    int_mb(k_detci_exb),
     $                    cdetci_ataba, cdetci_atabb, ntij, mseg,
     $                    dbl_mb(k_detci_h), dbl_mb(k_detci_g),
     $                    dbl_mb(k_detci_fscr),
     $                    g_civec, g_sigma )
      detci_saa_etime = detci_saa_etime + util_cpusec() - tx
c
c   alpha-beta
c
      tx = util_cpusec()
      call detci_sigmaab( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                    cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                    nekla, neklb, cdetci_osym,
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
      detci_sab_etime = detci_sab_etime + util_cpusec() - tx
      detci_sigma_etime = detci_sigma_etime + util_cpusec() - t0
c
c
c$$$      CALL DETCI_CIVEC_PRINT( CDETCI_NORB, CDETCI_NSYM, CDETCI_NELA, CDETCI_NELB,
c$$$     $                      NSTRA, NSTRB, CDETCI_OSYM,
c$$$     $                      CDETCI_ATABA, CDETCI_ATABB, SIGMA, 1.D-5 )
c
c
c
      return
      end








      subroutine detci_ciprecon( g_civec, g_workvec )
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
      integer g_civec
      integer g_workvec
      integer nstra, nstrb
c
c$$$  double precision xx, yy
      double precision tx
      double precision mxcnt
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_ciprecon with uninitialized common',0,
     &       INPUT_ERR)
c     
c  CI parameters
c     
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c  Apply preconditioning using orbital 
c  energies (Moller-Plesset denominators)
c
      call detci_diagscale( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                      cdetci_nelb, nstra, nstrb, cdetci_osym,
     $                      cdetci_ataba, cdetci_atabb,
     $                      cdetci_eref, cdetci_eps, g_civec )
c
c  Spin-adapt the vector
c
      if (cdetci_spinadapt) then
c$$$        PRINT*,'=== ENTRY VECTOR ==='
c$$$        CALL DETCI_PRINT(G_CIVEC,1.D-3)
CCCCC        CALL DETCI_RANDOM_ERROR( G_CIVEC )       ! INTRODUCE CONTAMINATION TO DEBUG
c$$$        XX = SQRT( GA_DDOT( G_CIVEC, G_CIVEC ) )

        tx = util_cpusec()
        call detci_spinadapt( g_civec, g_workvec )
        detci_spinadp_etime = detci_spinadp_etime + util_cpusec() - tx

c$$$        YY = SQRT( GA_DDOT( G_CIVEC, G_CIVEC ) )
c$$$        WRITE(6,988) (1.d0 - YY/XX)
c$$$ 988    FORMAT('Spin contamination:',E12.3)
c$$$        PRINT*,'=== AFTER SPIN ADAPT ==='
c$$$        CALL DETCI_PRINT(G_CIVEC,1.D-3)
      endif
c
c  Ensure vector is symmetry-adapted
c
c$$$      PRINT*,'=== AFTER PRECON ==='
c$$$      CALL DETCI_PRINT(G_CIVEC,1.D-3)
      if (cdetci_nsym.gt.1) then
        tx = util_cpusec()
        call detci_symmproject( cdetci_norb, cdetci_nsym,
     $                          cdetci_nela, cdetci_nelb, nstra,
     $                          nstrb, cdetci_osym, 
     $                          cdetci_ataba, cdetci_atabb,
     $                          cdetci_symstate, .true.,
     $                          mxcnt, g_civec )
        detci_symmadp_etime = detci_symmadp_etime + util_cpusec() - tx
      endif

c$$$      PRINT*,'=== AFTER SYMMPROJ ==='
c$$$      CALL DETCI_PRINT(G_CIVEC,1.D-3)

      return
      end







c
c  Print CI vector coefficients larger
c  than threshold
c
      subroutine detci_print( g_civec, thresh )
      implicit none
#include "errquit.fh"
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
     $  call errquit('detci_print with uninitialized common',0,
     &       INPUT_ERR)
c     
c  CI parameters
c     
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c
c
      call detci_civec_print( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                        cdetci_nelb, nstra, nstrb,
     $                        cdetci_osym, cdetci_ixmap,
     $                        cdetci_ataba, cdetci_atabb,
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
      subroutine detci_density( g_civec, onepdm, twopdm )
      implicit none
#include "errquit.fh"
c
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "detciP.fh"
#include "detci.fh"
#include "cdetci.fh"
#include "cdetcistats.fh"
c
      integer g_civec
      double precision onepdm(*)
      double precision twopdm(*)
c
      integer nexa, nexb
      integer nstra, nstrb
      integer nekla, neklb
      integer nn
      integer g_civect
      double precision tx, tx2
c
c  Check internals
c
      if (cdetci_valid .ne. CDETCI_MAGIC)
     $  call errquit('detci_density with uninitialized common',0,
     &       INPUT_ERR)
c     
c  CI parameters
c     
      nexa = (cdetci_norb-cdetci_nela+1)*cdetci_nela
      nexb = (cdetci_norb-cdetci_nelb+1)*cdetci_nelb
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
      nekla = detci_binomial( (cdetci_norb-1),(cdetci_nela-1))
      neklb = detci_binomial( (cdetci_norb-1),(cdetci_nelb-1))
      nn = cdetci_norb*cdetci_norb
c
c
c
      tx = util_cpusec()
*ga:1:0
      if (.not.ga_create( MT_DBL, nstra, nstrb, 'civec transp',
     $                    nstra, 0, g_civect ))
     $    call errquit('detci_density: cannot create CI transp.',0,
     &       GA_ERR)
      call ga_transpose( g_civec, g_civect )
c
c
c
      tx2 = util_cpusec() 
      call detci_onepdm( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                   cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                   cdetci_osym, int_mb(k_detci_exa),
     $                   int_mb(k_detci_exb), cdetci_ixmap,
     $                   g_civec, g_civect, onepdm )
      detci_density_onepdm = detci_density_onepdm + util_cpusec() - tx2
c
c
c
      tx2 = util_cpusec() 
      call dfill(nn*nn, 0.0d0, twopdm, 1)
      call detci_twopdm( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                   cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                   cdetci_osym, int_mb(k_detci_exa),
     $                   int_mb(k_detci_exb), cdetci_ixmap,
     $                   g_civec, g_civect, twopdm )
      detci_density_twopdm = detci_density_twopdm + util_cpusec() - tx2
c
c
c
      tx2 = util_cpusec() 
      call detci_twopdm_ab( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                      cdetci_nelb, nstra, nstrb, nexa, nexb,
     $                      nekla, neklb, cdetci_osym, cdetci_ixmap,
     $                      int_mb(k_detci_exa), int_mb(k_detci_exb),
     $                      cdetci_ataba, cdetci_atabb,
     $                      int_mb(k_detci_rhsscr),
     $                      int_mb(k_detci_lhsscr),
     $                      int_mb(k_detci_iscr),
     $                      dbl_mb(k_detci_cprime),
     $                      dbl_mb(k_detci_sprime),
     $                      g_civec, twopdm )
      detci_density_twopdmab=detci_density_twopdmab+util_cpusec()-tx2
      call ga_dgop(1, twopdm, nn*nn, '+')
c
c
c
      if (.not.ga_destroy( g_civect ))
     $   call errquit('detci_density: cannot destroy CI transp.',0,
     &       GA_ERR)
      detci_density_etime = detci_density_etime + util_cpusec() - tx
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
#include "errquit.fh"
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
     $  call errquit('detci_spin with uninitialized internals',0,
     &       INPUT_ERR)
c
c  CI parameters
c     
      call detci_ijmap( cdetci_norb, cdetci_nsym, cdetci_osym,
     $                  ntij, ijmap )
      nexa = (cdetci_norb-cdetci_nela+1)*cdetci_nela
      nexb = (cdetci_norb-cdetci_nelb+1)*cdetci_nelb
      nstra = detci_binomial( cdetci_norb, cdetci_nela )
      nstrb = detci_binomial( cdetci_norb, cdetci_nelb )
c
c
c
      call detci_spadpt( cdetci_norb, cdetci_nsym, cdetci_nela,
     $                   cdetci_nelb, nstra, nstrb,
     $                   cdetci_osym, cdetci_ataba, cdetci_atabb,
     $                   nexa, nexb, int_mb(k_detci_exa),
     $                   int_mb(k_detci_exb), g_civec, g_wvec )
c
c
      return
      end

      








      subroutine detci_symmproj( g_civec, mxcnt, oscreen )
      implicit none
#include "errquit.fh"
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
     $  call errquit('detci_symmproj with uninitialized common',0,
     &       INPUT_ERR)
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
     $                          nstrb, cdetci_osym, 
     $                          cdetci_ataba, cdetci_atabb,
     $                          cdetci_symstate, oscreen,
     $                          mxcnt, g_civec )
      endif

      return
      end









      subroutine detci_moint_copy( norb, ntij, irmap, ijmap, orbsym,
     $                             h, g, h1, g1 )
      implicit none
#include "errquit.fh"
#include "global.fh"
      integer norb
      integer ntij
      integer irmap(norb)
      integer ijmap(norb,norb)
      integer orbsym(norb)
      double precision h(ntij)
      double precision g(ntij,ntij)
      double precision h1(ntij)
      double precision g1(ntij,ntij)
c      
      integer i, j, k, l
      integer ij, rij, kl, rkl
      integer ijsym, ijklsym, nmixed
      double precision TOL
#include "symmdef.fh"
#include "bitops.fh"
#include "symmmul.fh"
      data TOL/1.d-6/           ! rjh - was 1d-9 but this is not met - why?
c      
      nmixed = 0
      do i=1,norb
        do j=1,i
          ij = ijmap(i,j)
          rij = ijmap(irmap(i),irmap(j))
          h1(ij) = h(rij)
          ijsym = MULX(orbsym(i),orbsym(j))
          if (ijsym.ne.1) then
            if (abs(h1(ij)).gt.TOL) then
               write(6,*) ' H1 contaminated ', i, j, h1(ij)
               nmixed = nmixed + 1
            endif
            h1(ij) = 0.0d0
          endif
          do k=1,norb
            do l=1,k
              kl = ijmap(k,l)
              rkl = ijmap(irmap(k),irmap(l))
              g1(ij,kl) = g(rij,rkl)
              ijklsym = MULX(ijsym,MULX(orbsym(k),orbsym(l)))
              if (ijklsym.ne.1) then
                 if (abs(g1(ij,kl)).gt.TOL) then
                    nmixed = nmixed + 1
                    write(6,*) ' G1 contaminated ', i,j,k,l,g1(ij,kl)
                 endif
                 g1(ij,kl) = 0.0d0
              endif
            enddo
          enddo
        enddo
      enddo
      if ((nmixed.gt.0).and.(ga_nodeid().eq.0)) then
        call errquit('DETCI: MO ints symmetry contaminated',nmixed,
     &       GEOM_ERR)
      endif
c
      end

