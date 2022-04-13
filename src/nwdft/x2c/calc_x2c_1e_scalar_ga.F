      subroutine calc_x2c_1e_scalar_ga (rtdb, 
     &  nexc, g_dens, g_vxc, g_x2c, g_u, g_u2c, g_c2u, toll_s)

c     -----------------------------------------------------------------
c     Purpose: calculate scalar X2C one-electron Hamiltonian matrix

c     Integrals are supposed to be initialized by the calling
c     routine

c     Driver (c) 2012, 2015, 2022 J. Autschbach, jochena@buffalo.edu
c     'relham' code used for the original (2012) implementation
c     (c) 2012 Daoling Peng. This code is now obsolete in NWChem.
c     Contracted basis extension: J. Autschbach, 2015
c     New code using GA exclusively (c) 2022 J. Autschbach
      
c     Code is not yet (2022) well tested. Expect bugs & email me at the
c     address shown above if you find one. Please provide a minimal
c     working example (in-,output) for your bug report, if possible.

c     The relevant equations implemented here can be found in the
c     following works, among others:
      
c     Autschbach, Peng, Reiher, J. Chem. Theory Comput. 2012, 8,
c     4239–4248. (initial NWChem implementation, field gradient calcs.)

c     Liu & Peng, J. Chem. Phys. 2009, 131, 031104. (X2C formalism in
c     non-orthogonal basis)

c     Peng & Reiher, Theor. Chem. Acc. 2012, 131, 1081. (comparison of
c     X2C approaches)

c     Autschbach, J., Quantum Theory for Chemical Applications, Oxford
c     University Press (2021), Chapter 24. (further details re. the
c     formalism)

c     The workflow for mapping onto an uncontracted basis and back is
c     similar to int_1edk_ga.F. We have the problem here that anything
c     that is calculated by numerical integration needs to be done in
c     the contracted basis, because of the way the grid routines are
c     written. I.e. unlike in int_1edk_ga.F we also need a
c     back-transformation.  For the time being, we use a naive
c     left-inverse for this.  Eventually, the numerical integrations
c     will be replaced completely.
c     -----------------------------------------------------------------

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "dftpara.fh"
#include "cdft.fh"
#include "rel_consts.fh"
#include "cint1cache.fh"
#include "numerical_constants.fh"
#include "x2c_scf.fh"

c     subroutine arguments:

      integer rtdb              ! [inp] runtime database GA handle
      integer nexc              ! [inp] no. of XC terms
      integer g_dens(2)         ! [inp] atom densities for model pot.
      integer g_vxc(4)          ! scratch 
      integer g_x2c             ! [out] exact decoupling 1e operator 
      integer g_u(2)            ! [out] decoupling transformation
      integer g_u2c             ! [out] basis trafo
      integer g_c2u             ! [out] basis trafo
      double precision toll_s   ! [inp] Smat diag. tolerance

c     GA handles used for temp storage in this routine:
      integer g_t, g_s, g_w, g_v, g_temp, g_tmpv, g_X
      integer g_smix, g_densu
      integer g_tm12, g_sm12, g_h, g_m

c     local MA handles

      integer l_eval, k_eval
      integer l_scr, k_scr

c     other local variables:

      integer ubas, cbas
      integer nu, nu2, nunu, nc
      integer me, type, i, j, ntemp
      double precision sum, energy, rtemp
      double precision c_2c, c_2c2, c_4c2

      double precision ga_trace_diag
      external ga_trace_diag

      integer info
      integer ga_cholesky,ga_llt_i
      external ga_cholesky,ga_llt_i
      
      double precision tol2e
      logical master, analytic_pvp, modelpot_2e, use_ga, oskel_sav
      integer dbg
      double precision jfac(1),kfac(1)

      character*13 pname

c     =================================================================

c     dbg = 1: debug numbers and flow. 2: debug also matrices
c     dbg = 0: no debug output
      dbg=1
      me=ga_nodeid()
      master = me.eq.0
      pname = "x2c_1e_scalar"
      
      if (dbg>0.and.master) write(luout,*) 'entering '//pname
      
      call ga_sync()
      call ga_zero(g_x2c)

      tol2e=10.d0**(-itol2e)

c     src/NWints/dk/int_1edk_ga.F uses similar transformations as used
c     here. There, oskel is deliberately set to .false., so I'm guessing
c     we should do the same here but I'm not sure yet

      oskel_sav = oskel
      oskel = .false.


c     temporarily add alpha and beta model densities:
c     TO BE REPLACED WITH 4-COMPONENT MODEL DENSITIES

      if (ipol.gt.1) 
     &  call ga_dadd(1.d0, g_dens(1), 1.d0, g_dens(2), g_dens(1))

c     ---------------
c     initializations
c     ---------------

      ubas = x2c_bas            ! uncontracted basis
      cbas = ao_bas_han         ! contracted basis

c     make sure the input GAs have the correct dimension

      call ga_inquire(g_x2c,type,nc,nc)
      if (nc.ne.nbf_ao) then
        call errquit(pname//': g_x2c has wrong dimension',nc,GA_ERR)
      end if

      if (.not. bas_numbf(cbas,ntemp))
     &  call errquit(pname//': bas_ao error',ntemp, BASIS_ERR)
      if (ntemp.ne.nc)
     &  call errquit(pname//': nc.ne.ntemp',ntemp, BASIS_ERR)
      if (.not. bas_numbf(ubas,nu))
     &  call errquit(pname//': nu error',nu, BASIS_ERR)

      if (master.and.dbg>0) write (luout,*)
     &  pname//': nu, nc =',nu,nc

      if (nu.lt.nc) call errquit(
     &  'x2c uncontracted basis smaller than contracted', 1,
     &  RTDB_ERR)

      lx2c_u2c = (nu.gt.nc)
      if (master) then
        if (lx2c_u2c) then
          write (luout,*)
     &      'x2c INFO: the operator will be calculated in an'
          write (luout,*)
     &      ' uncontracted auxiliary basis'
        else
          write (luout,*)
     &      'x2c INFO: Assuming that the AO basis is fully uncontracted'
        end if
      end if

c     no of elements in the AO matrices and dimension of mDE:

      nunu = nu * nu
      nu2 = 2* nu

c     ---------------------------
c     read some options from rtdb
c     ---------------------------

c     analytic pVp integrals are default, but it can be disabled

      analytic_pvp = .true.
      if (.not. rtdb_get(rtdb, 'x2c:analyticpvp',
     &  mt_log, 1, analytic_pvp)) continue

c     model potential for 2e terms?

      modelpot_2e = .false.
      if (.not. rtdb_get(rtdb, 'x2c:modelpot2e',
     &  mt_log, 1, modelpot_2e)) continue

c     if analyticpvp is .F. or if we need the model potential,
c     we also need a contracted-to-uncontracted basis transformation
c     unless the two bases are identical

      lx2c_c2u = ((.not.analytic_pvp .or. modelpot_2e).and.lx2c_u2c)
      if (master.and.lx2c_c2u) write (luout,*)
     &  'x2c INFO: x2c options require constructing the c2u trafo'

c     check rtdb switch to see whether we should perform decoupling
c     in GA, or call relham routines with MA arrays

      use_ga = .false.
      if (.not.rtdb_get(rtdb, 'x2c:usega', mt_log, 1,
     &  use_ga)) continue
      if (master.and.use_ga) write(luout,*)
     &  pname//': will perform decoupling in GA mode'
      
c     -----------------------------------------
c     allocate GAs for T, V, S, and for W = pVp
c     -----------------------------------------

c     create S matrix GA. Do we need to use that atom-block stuff
c     here, like they do in the SCF calling routine dft_main0d.F ???

      if (.not. ga_create(MT_DBL,nu,nc,'g_s',1,1,g_s)) call
     &  errquit(pname//': ga_create failed g_s',11,GA_ERR)
      call ga_zero(g_s)

c     the matrices for T, V, W are created as duplicates of g_s:
      
      if (.not.ga_duplicate(g_s, g_t,'g_t')) call 
     &  errquit(pname//': ga_create failed g_t',12,GA_ERR)
      call ga_zero(g_t)

      if (.not.ga_duplicate(g_s,g_v, 'g_v')) call 
     &  errquit(pname//': ga_create failed g_v',13,GA_ERR)
      call ga_zero(g_v)

      if (.not.ga_duplicate(g_s, g_w, 'g_w')) call
     &  errquit(pname//': ga_create failed g_w',14,GA_ERR)  
      call ga_zero(g_w)
      
c     ---------------------------------------
c     Analytic one-electron matrices: S, T, V
c     ---------------------------------------

      if (ocache) call int_1e_uncache_ga()
      call int_1e_ga(ubas, ubas ,g_s, 'overlap'   ,oskel)
      call int_1e_ga(ubas, ubas, g_t, 'kinetic'  , oskel)
      call int_1e_ga(ubas, ubas, g_v, 'potential', oskel)

c     symmetrize, just in case ... 
      call ga_symmetrize(g_s)
      call ga_symmetrize(g_t)
      call ga_symmetrize(g_v)

      call ga_sync()

      if (master.and.dbg>0) write(luout,*) pname//': T,V,S done'

c     ------------------------------------------------------
c     later we need a transformation from the contracted (c)
c     to the uncontracted (u) basis, and vice versa for the
c     density matrix sooner. the sequence comes
c     from routine int_1edk_ga
c     ------------------------------------------------------

c     u2c transforms the basis from uncontracted to contracted
c     (really, these are just the contraction coefficients,
c     but this is how it is done in the DKH routines)

      if (lx2c_u2c) then

        if (.not. ga_create(MT_DBL,nu,nc,'u2c',1,1,g_u2c))
     &    call errquit(pname//': create u2c failed',51,GA_ERR)
        call ga_zero(g_u2c)
        
        if (.not. ga_create(MT_DBL,nu,nc,'smix',1,1,g_smix))
     &    call errquit(pname//': create mixed overlap failed',52,GA_ERR)
        call ga_zero(g_smix)
        
        if (.not. ga_create(MT_DBL,nu,nu,'temp',1,1,g_temp))
     &    call errquit(pname//': error creating temp',53,GA_ERR)
        
        if (master.and.dbg>0) write(luout,*) pname//': u2c mem done'
        
        call int_1e_ga(ubas, cbas, g_smix, 'overlap' ,oskel)
        call ga_copy(g_s,g_temp)
        call ga_matpow(g_temp, -1.0d0, 1d-12)
        call ga_dgemm('n','n',nu,nc,nu,1.0d0,
     &    g_temp, g_smix, 0.0d0, g_u2c)
        
        if (.not. ga_destroy(g_temp))
     &    call errquit(pname//': error destroy temp',53, GA_ERR)
        
        if (.not. ga_destroy(g_smix))
     &    call errquit(pname//': error destroy smix',52, GA_ERR)
        
        if (master.and.dbg>0) write(luout,*) pname//': u2c done'
        
      end if                    ! lx2c_u2c

c     c2u transforms the basis from contracted to uncontracted.
c     we have to do this as a generalized left-inverse of u2c.
c     Use [(A^T A)**(-1) A^T] A = 1 to define the left-inverse.

      if (lx2c_c2u) then

        if (.not. ga_create(MT_DBL,nc,nu,'c2u',1,1,g_c2u))
     &    call errquit(pname//': create c2u failed',54,GA_ERR)
        call ga_zero(g_c2u)
        
        if (.not. ga_create(MT_DBL,nc,nc,'temp',1,1,g_temp))
     &    call errquit(pname//': failed creating temp',55,GA_ERR)
        
        if (master.and.dbg>0) write(luout,*) pname//': c2u mem done'
        
        call ga_dgemm('t','n',nc,nc,nu,1.0d0,g_u2c,g_u2c,
     &    0.0d0,g_temp)   
        
        call ga_matpow(g_temp, -1.0d0, 1d-12)
        if (master.and.dbg>0) write(luout,*) pname//': c2u matpow done'
        
        call ga_dgemm('n','t',nc,nu,nc,1.0d0,g_temp,g_u2c,
     &    0.0d0,g_c2u)
        
        if (.not. ga_destroy(g_temp))
     &    call errquit(pname//': error destroy temp',55, GA_ERR)
        
        if (master.and.dbg>0) write(luout,*) pname//': c2u done'
        
c       check c2u * u2c
        
        if (.not. ga_create(MT_DBL,nc,nc,'temp',1,1,g_temp))
     &    call errquit(pname//': error creating temp',56,GA_ERR)
        
        call ga_dgemm('n','n',nc,nc,nu,1.0d0,g_c2u,g_u2c,
     &    0.0d0,g_temp)

        if (dbg>1) then
          call ga_print(g_u2c)
          call ga_print(g_c2u)
          call ga_print(g_temp)
        end if

c       check the quality of the left-inverse. the product
c       c2u*u2c must be a unit matrix of dimension nc

        sum = ga_trace_diag(g_temp)  
        if (master) write (luout,'(1x,a,1x,e15.7)')
     &    'x2c c2u*u2c: ',sum

        if (abs(sum-real(nc)).gt.1.0d-10) then
          if (master) write (luout,*)
     &      'x2c WARNING: left-inverse in basis trafo is inaccurate'
        end if

        if (.not. ga_destroy(g_temp))
     &    call errquit(pname//': error destroy temp',56, GA_ERR)

        if (master.and.dbg>0) write(luout,*) pname//': c2u done'

      end if                    ! lx2c_c2u

      if (lx2c_u2c) then

c       transform density matrix to the uncontracted basis:
c       c2u is how the basis would transform in this case.
c       the density matrix is based on the MO coefficients,
c       P = C n C^T, with C being a collections of MO vectors
c       in the c-basis. We do a passive transformation,
c       since the MOs are to remain the same. In this case,
c       the transformation of the MO vectors is contra-
c       gredient to the basis and we need u2c:
c       P(u) = u2c P(c) u2c^T. 

        if (.not. ga_create(MT_DBL,nc,nu,'temp',1,1,g_temp))
     &    call errquit(pname//': error creating temp',57,GA_ERR)
        
        if (.not. ga_create(MT_DBL,nu,nu,'densu',1,1,g_densu))
     &    call errquit(pname//': error creating densu',58,GA_ERR)
        
        if (master.and.dbg>0) write(luout,*) pname//': dens c2u start'
        
        call ga_zero(g_densu)
        call ga_dgemm('n','t',nc,nu,nc, 1.0d0,
     &    g_dens(1),g_u2c,0.0d0,g_temp)
        if (master.and.dbg>0) write(luout,*) pname//': dens c2u step1'
        call ga_dgemm('n','n',nu,nu,nc ,1.0d0,
     &    g_u2c, g_temp, 0.0d0,g_densu)
        call ga_symmetrize(g_densu)
        
        if (.not. ga_destroy(g_temp))
     &    call errquit(pname//': error destroy temp',57, GA_ERR) 
        
        if (master.and.dbg>0) write(luout,*) pname//': dens c2u done'

      else

c       if there is no need for the u2c trafo, we
c       just assign the GA handle of g_dens(1) to g_densu 

        g_densu = g_dens(1)

      end if                    ! lx2c_u2c

      if (master.and.dbg>0) write (luout,*)
     &  'g_densu, g_dens(1) = ',g_densu,g_dens(1)

c     ------------------------------------
c     use sum of atomic densities to check
c     that the AO matrices are OK
c     ------------------------------------

      energy = 0d0

      if (.not. ga_create(MT_DBL,nu,nu,'g_temp',1,1,g_temp)) call
     &  errquit(pname//': ga_create g_temp failed',59,GA_ERR)
      call ga_zero(g_temp)
     
      if (master.and.dbg>0) write(luout,*) pname//': debug S:'

      call ga_dgemm('n','n',nu,nu,nu, 1.0d0,
     &  g_densu,g_s,0.0d0,g_temp)
      
      sum = ga_trace_diag(g_temp)
      
      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': tr(P * S) : ',sum

      if (master.and.dbg>0) write(luout,*) pname//': debug T:'

      call ga_dgemm('n','n',nu,nu,nu, 1.0d0,
     &  g_densu,g_t,0.0d0,g_temp)
      
      sum = ga_trace_diag(g_temp)
      
      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': tr(P * T) : ',sum
      energy = energy + sum

      if (master.and.dbg>0) write(luout,*) pname//': debug V:'

      call ga_dgemm('n','n',nu,nu,nu, 1.0d0,
     &  g_densu,g_v,0.0d0,g_temp)

      sum = ga_trace_diag(g_temp)
      
      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': tr(P * V) : ',sum
      energy = energy + sum

      if (.not. ga_destroy(g_temp)) call errquit
     &  (pname//': temp ga corrupt?',59, GA_ERR)

      if (master) write (luout,'(1x,a,1x,e15.7)')
     &  pname//': 1e energy : ',energy

c     ----------------------------------------------------------
c     as this routine gets embedded in the scalar DFT framework, 
c     we make it work similar to ZORA. I.e., for compatibility 
c     we need to SUBTRACT the nonrel H from the operator. 
c     this is done here, upfront. 
c     ----------------------------------------------------------

c     if lx2c_u2c is .T., we now need to allocate an array for the
c     x2c operator in the uncontracted basis => g_X
c     otherwise, we simply use the GA handle of g_x2c.

      if (lx2c_u2c) then
        if (.not.ga_duplicate(g_s, g_X, 'g_X')) call
     &  errquit(pname//': ga_create failed g_X',61,GA_ERR)  
        call ga_zero(g_X)
      else
        g_X = g_x2c
      end if                    ! lx2c_u2c

      call ga_dadd(1.d0, g_X, -1.0d0, g_t, g_X)
      call ga_dadd(1.d0, g_X, -1.0d0, g_v, g_X)

c     ------------------------------------------------------
c     In the case of using a model Coulomb potential in the
c     transformation (V and W integrals), we need to
c     subtract the corresponding Coulomb potential from the 
c     one-electron operator and add it to V (and later in W)
c     ------------------------------------------------------

      if (modelpot_2e) then
        
        if(.not.ga_create(mt_dbl,nu,nu,'tmpv',0,0,g_tmpv))
     &    call errquit(pname//': error creating ga',141, GA_ERR)
        
        jfac(1) = 1.0d0
        kfac(1) = 0.0d0
        
        call ga_zero(g_tmpv)
        
        if (master.and.dbg>0) write (luout,*)
     &    'x2c_1e: calling fock2e for model potential'
        
        call fock_2e(geom, ubas, 1, jfac, kfac,
     &    tol2e, oskel, g_densu, g_tmpv, .false.)
        
        call ga_sync()
        
c       subtract from the 1e Hamiltonian:
        call ga_add(1.0d0, g_X, -1.0d0, g_tmpv, g_X)
        
c       add Coul. pot. to g_v:
        call ga_add(1.0d0, g_v, 1.0d0, g_tmpv, g_v)
        
        if (.not. ga_destroy(g_tmpv)) call errquit
     &    (pname//': ga corrupt?',141, GA_ERR)
        
      end if                    ! modelpot_2e
     
c     g_densu not needed past this point
      if (lx2c_u2c) then
        if (.not. ga_destroy(g_densu))
     &    call errquit(pname//': error destroy densu',58, GA_ERR) 
      else
        if (g_densu.ne.g_dens(1)) call errquit(
     &    pname//': g_densu inconsistency 1',g_densu, UNKNOWN_ERR)
      end if
      
      
      call ga_zero(g_w)
      
c     ---------------------------------------------------------------
c     nuclear potential only: calculate analytic pV dot p integrals
c     and add them to the array g_w(:).
c     the calls to the integral routines below have been adapted
c     from the douglas-kroll code in nwchem.
c     ---------------------------------------------------------------
      
      if (analytic_pvp) then
        if (master) write (luout,*)
     &    pname//': calculating analytic p.Vp integrals'
        
c       add scalar part p.Vp to g_w:
        call int_1e_ga(ubas, ubas ,g_w, 'pvp' ,oskel)
        call ga_symmetrize(g_w)
      end if ! analytic_pvp

c     --------------------------------------------------
c     if needed: calculate all or part of
c     (sig.p)V(sig.p) integrals by numerical integration,
c     similar to how the ZORA spin-orbit arrays are done.
c     We have the option to skip the nuclear part, if those
c     integrals are done analytically
c     --------------------------------------------------

      if (analytic_pvp .and. .not.modelpot_2e) then
        if (master) write(luout,*)
     &    pname//': skipping numerical W integrals'
      else
        
        if (master) write(luout,*) pname//': calculating W numerically'
        if (modelpot_2e .and. master)
     &    write(luout,*) '     (2e model pot. only)'
        
        
        if (lx2c_u2c) then
          
c         here we need to use the contracted basis ...
c         use g_tmpv to hold the array.
c         g_vxc has the right dimension already
          if (dbg>0.and.master) write (luout,*) 'g_w numint c2u start'
          if(.not.ga_create(mt_dbl,nc,nc,'tmpw',0,0,g_tmpv))
     &      call errquit(pname//': error creating tmpw',147, GA_ERR)
          call ga_zero(g_tmpv)
          call x2c_getw_scalar(rtdb, g_dens(1), g_tmpv, g_vxc, nexc,
     &      analytic_pvp, .not.modelpot_2e)
          if (dbg>0.and.master) write (luout,*) 'g_w numint c2u end'
c         transform to ubas: tmpv = c2u^T g_tmpv c2u, add to g_w

          if(.not.ga_create(mt_dbl,nc,nu,'temp',0,0,g_temp))
     &      call errquit(pname//': error creating temp',148, GA_ERR)
          call ga_zero(g_temp)
          call ga_dgemm('n','n',nc,nu,nc, 1.0d0,
     &      g_tmpv,g_c2u,0.0d0,g_temp)
          if (master.and.dbg>0) write(luout,*) pname//': g_w c2u step1'
          call ga_dgemm('t','n',nu,nu,nc ,1.0d0,
     &      g_c2u, g_temp, 1.0d0,g_w)
          if (master.and.dbg>0) write(luout,*) pname//': g_w c2u step2'
          if (.not. ga_destroy(g_tmpv)) call errquit
     &      (pname//': ga corrupt?',147, GA_ERR)
          if (.not. ga_destroy(g_temp)) call errquit
     &      (pname//': ga corrupt?',148, GA_ERR)
        else
c         cbas = ubas. Just continue as usual
          if (dbg>0.and.master) write (luout,*) 'g_w numint start'
          call x2c_getw_scalar(rtdb, g_dens(1), g_w, g_vxc, nexc,
     &      analytic_pvp, .not.modelpot_2e)
          if (dbg>0.and.master) write (luout,*) 'g_w numint end'
        end if                  ! lx2c_u2c

      end if                    ! need analytic integration
      
c     W is the pVp matrix, incl model pot if requested:
      if (dbg>1) call ga_print(g_w)

c     sanity checks:

      if (.not.lx2c_u2c) then
        if (g_X.ne.g_x2c) call errquit(
     &    pname//': g_X/g_x2c inconsistency 1',g_X-g_x2c, UNKNOWN_ERR)
      end if

c     -----------------------------------------------
c     Input data are all in place for decoupling step
c     -----------------------------------------------

c     define some constants related to speed of light:
c     (define them here directly as inverses, which we need later)
      
      c_2c  = 1.0d0/(cau*2.d0) ! 1/2c
      c_2c2 = 1.0d0/(cau*cau*2.d0) ! 1/2c**2
      c_4c2 = 1.0d0/(cau*cau*4.d0) ! 1/4c**2
      
c     scale W by 1/(4c^2)
      
      call ga_scale(g_w, c_4c2)
      
c     assemble the 'scalar 4-component' mDE Hamiltonian
      
      if (.not. ga_create(MT_DBL,nu2,nu2,'x2c_h',1,1,g_h))
     &  call errquit(pname//': ga_create failed g_h',66, GA_ERR)
      call ga_zero(g_h)
      
      call ga_copy_patch('n', ! L-L block = V
     &  g_v, 1, nu, 1, nu,
     &  g_h, 1, nu, 1, nu)
      
      call ga_copy_patch('n', ! L-S block = T
     &  g_t, 1, nu, 1, nu,
     &  g_h, 1, nu, nu+1, nu2)
      
      call ga_copy_patch('n', ! S-L block = T
     &  g_t, 1, nu, 1, nu,
     &  g_h, nu+1, nu2, 1, nu)
      
      call ga_add_patch( ! S-S block = W/(4c**2)-T
     &  1.0d0, g_w, 1, nu, 1, nu,
     &  -1.0d0, g_t, 1, nu, 1, nu,
     &  g_h, nu+1, nu2, nu+1, nu2)
      
c     assemble the metric for the mDE:

      if (.not. ga_create(MT_DBL,nu2,nu2,'x2c_m',1,1,g_m))
     &  call errquit(pname//': ga_create failed g_m',67, GA_ERR)
      
      call ga_zero(g_m)
      
      call ga_copy_patch('n', ! L-L block = S
     &  g_s, 1, nu, 1, nu,
     &  g_m, 1, nu, 1, nu)
      
c     call ga_scale(g_t, c_2c2)        
      call ga_copy_patch('n', ! S-S block = T/(2c**2)
     &  g_t, 1, nu, 1, nu,
     &  g_m, nu+1, nu2, nu+1, nu2)
      call ga_scale_patch(
     &  g_m, nu+1, nu2, nu+1, nu2, c_2c2)

c     diagonalize the mDE Hamiltonian:

      if (.not. ma_push_get(mt_dbl, nu2, 'eval', l_eval, k_eval))
     &  call errquit(pname//': failed allocation eval', nu2,
     &  MA_ERR)
      
      if (.not.ga_duplicate(g_h, g_temp,'g_temp')) call 
     &  errquit(pname//': ga_dup failed g_temp',68,GA_ERR)
      call ga_zero(g_temp)
      
      if (dbg>1) call ga_print(g_h)
      if (dbg>1) call ga_print(g_m)
      
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
c     
c     SCALAPACK piece needs code similar to what's
c     in src/NWints/dk/int_1edk_ga.F
      if (master .and. dbg>0) write (luout,*)
     &  'calling ga_pdsygv replacement sequence of routines'
c     call ga_pdsygv(g_h, g_m, g_temp, dbl_mb(k_eval))
      info= ga_cholesky('L', g_m)
      info = ga_llt_i('L', g_m, -1)
      info = ga_cholesky('L', g_m, -1)
      
      call ga_dgemm('t','n',nu2, nu2, nu2, 1.0d0,
     &  g_m, g_h, 0.0d0, g_temp)
      
      call ga_dgemm('n','n', nu2, nu2, nu2, 1.0d0,
     &  g_temp, g_m, 0.0d0, g_h)
      
      call ga_pdsyev(g_h, g_temp, dbl_mb(k_eval), 0)
      
      call ga_dgemm('n','n', nu2, nu2, nu2, 1.0d0,
     &  g_m, g_temp,   0.0d0, g_h)
      call ga_copy (g_h, g_temp)
      
#else
      if (master .and. dbg>0) write (luout,*) 'calling ga_diag'
      call ga_diag(g_h, g_m, g_temp, dbl_mb(k_eval))
#endif
#else
      if (master .and. dbg>0) write (luout,*) 'calling ga_diag_seq'
      call ga_diag_seq(g_h, g_m, g_temp, dbl_mb(k_eval))
#endif

c     we don't need the mDE eigenvalues further
      
      if (.not.ma_chop_stack(l_eval)) call
     &  errquit(pname//': ma_chop_stack failed k_eval',l_eval,MA_ERR)
      if (dbg>1) call ga_print(g_temp)
      
c     destroy arrays for mDE Hamiltonian and metric,
c     allocate the same handles again for scratch space
      
      if (.not. ga_destroy(g_h))
     &  call errquit(pname//': ga_destroy failed g_h',66, GA_ERR)
      
      if (.not. ga_destroy(g_m))
     &  call errquit(pname//': ga_destroy failed g_m',67, GA_ERR)
      
      if (.not.ga_duplicate(g_s, g_h, 'g_h (temp)')) call 
     &  errquit(pname//': ga_dup failed g_h(temp)',71,GA_ERR)
      
      if (.not.ga_duplicate(g_s, g_m, 'g_m (temp)')) call 
     &  errquit(pname//': ga_dup failed g_m(temp)',72,GA_ERR)
      
c     set up and perform the decoupling.    
      
      call ga_copy_patch('n', 
     &  g_temp, 1, nu, nu+1, nu2,
     &  g_m, 1, nu, 1, nu) ! g_m now holds array A
      
      call ga_copy_patch('n', 
     &  g_temp, nu+1, nu2, nu+1, nu2,
     &  g_h, 1, nu, 1, nu) ! g_h now holds array B

      if (dbg>1) call ga_print(g_m)
      if (dbg>1) call ga_print(g_h)
      
c     destroy temp GA and reallocate with nubas**2 elements
      
      if (.not. ga_destroy(g_temp))
     &  call errquit(pname//': ga_destroy failed temp',68, GA_ERR)
      
      if (.not.ga_duplicate(g_s, g_temp,'g_temp')) call 
     &  errquit(pname//': ga_dup failed g_temp',73,GA_ERR)    

c     solve matrix equation X A = B for X, store in g_h
c     that is, solve A' X' = B' (' = Transpose), then transpose.
c     matrix A is square, but not symmetric.
      
      call ga_transpose(g_m, g_temp)
      call ga_copy(g_temp, g_m)
      call ga_transpose(g_h, g_temp)
      call ga_copy(g_temp, g_h) ! keep copy of B' in g_temp, too
      
#ifdef SCALAPACK
      if (dbg>0.and.master) write(luout,*) 'calling ga_solve'
      info = ga_solve(g_m, g_h) ! g_h now has solution X
      if (dbg>0.and.master) write(luout,*) 'info from ga_solve:', info
#else
      if (dbg>0.and.master) write(luout,*) 'calling ma_solve'
      call ma_solve(g_m, g_h)
#endif      

c     test quality of the solution:
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_m, g_h, -1.d0, g_temp) ! g_temp = A'X'-B', should be zero
      call ga_maxelt(g_temp, rtemp)
      if (master) write(luout,*)
     &  'largest element in A''X''-B'': ',rtemp
      call ga_transpose(g_h, g_temp)
      call ga_copy(g_temp, g_h)


      if (dbg>1) call ga_print(g_h)
  
c     assemble S + X(+)TX/(2c**2) in g_s
c     assemble X(+)T + TX - X(+)TX + V in g_t
      
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_t, g_h, 0.d0, g_temp)
      
      call ga_dgemm('T', 'N', nu, nu, nu, 1.d0,
     &  g_h, g_t, 0.d0, g_m) 
      
      call ga_copy(g_m, g_t) 
      
      call ga_add(1.0d0, g_t, 1.0d0, g_temp, g_t)
      
      call ga_dgemm('T', 'N', nu, nu, nu, 1.d0,
     &  g_h, g_temp, 0.d0, g_m) 
      
      call ga_add (1.0d0, g_t, -1.0d0, g_m, g_t) 
      if (dbg>1) call ga_print(g_t)
      
      call ga_add (1.0d0, g_t,  1.0d0, g_v, g_t) 
      
      call ga_copy(g_s, g_temp) ! save a copy of S in g_temp
      
      call ga_scale(g_m, c_2c2) 
      
      call ga_add (1.0d0, g_s, 1.0d0, g_m, g_s)
      
      if (dbg>1) then
        if (master) write (luout,*) 'tilde S = S + X''TX/(2c^2):'
        call ga_print(g_s)
        if (master) write (luout,*) 'X''T + TX - X''TX + V:'
        call ga_print(g_t)
      end if
      
c     now assemble the metric correction for the FW Hamiltonian
c     g_m, g_v can be used for scratch
c     g_temp still holds S
c     g_h still holds X
      
      call ga_copy (g_temp, g_v)
      call ga_matpow(g_v, -FP_Half, 1d-12) 
      
      call ga_dgemm('T', 'N', nu, nu, nu, 1.d0, 
     &  g_v, g_s, 0.d0, g_m)      
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_m, g_v, 0.d0, g_s)
    
      call ga_matpow(g_s, -FP_Half, 1d-12)
            
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_v, g_s, 0.d0, g_m)
      
      call ga_matpow(g_temp, FP_Half, 1d-12) 
      
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_m, g_temp, 0.d0, g_s) ! g_s = = UL
      
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_h, g_s, 0.d0, g_m) ! g_m = US   g_h, g_v, g_temp now free
       
c     construct the FW Hamiltonian in g_h. 
c     we then need to ADD it to g_X since it might 
c     contain -V(coul) already. 
c     also, the g_X array already has T+V subtracted 
c     such that it ties in with the zora code
      
      call ga_dgemm('T', 'N', nu, nu, nu, 1.d0,
     &  g_s, g_t, 0.d0, g_temp)      
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_temp, g_s, 0.d0, g_h)

      call ga_dgemm('T', 'N', nu, nu, nu, 1.d0,
     &  g_m, g_w, 0.d0, g_temp)
      call ga_dgemm('N', 'N', nu, nu, nu, 1.d0,
     &  g_temp, g_m, 1.d0, g_h)
      
      if (dbg>1) call ga_print(g_h)

      call ga_add(1.0d0, g_h, 1.0d0, g_X, g_X)

c     free up some memory:
      
      if (.not. ga_destroy(g_temp))
     &  call errquit(pname//': ga_destroy failed temp',73, GA_ERR)
      
      if (.not. ga_destroy(g_h))
     &  call errquit(pname//': ga_destroy failed g_h',71, GA_ERR)
      
      if (.not. ga_destroy(g_t))
     &  call errquit(pname//': ga_destroy failed t',12, GA_ERR)
      
      if (.not. ga_destroy(g_v))
     &  call errquit(pname//': ga_destroy failed v',13, GA_ERR)
      
      if (.not. ga_destroy(g_w))
     &    call errquit(pname//': ga_destroy failed w',14, GA_ERR)
        

c     -----------------------------------------------------
c     save decoupling transformation, if requested by input
c     this is kept in the uncontracted basis.
c     if we don't keep the transformation we also don't
c     need to keep the u2c [or c2u] transformation
c     -----------------------------------------------------
      
      if (lx2c_keepU) then
        if (.not. ga_create(MT_DBL,nu,nu,'x2c_UL',1,1,g_u(1)))
     &    call errquit(pname//': ga_create failed UL',75, GA_ERR)
        if (.not. ga_create(MT_DBL,nu,nu,'x2c_US',1,1,g_u(2)))
     &    call errquit(pname//': ga_create failed US',76, GA_ERR)
        
        call ga_copy(g_s, g_u(1))
        call ga_copy(g_m, g_u(2))
        call ga_scale(g_u(2), c_2c)
      end if ! lx2c_keepU

      if (.not. ga_destroy(g_s))
     &  call errquit(pname//': ga_destroy failed s',11, GA_ERR)
      
      if (.not. ga_destroy(g_m))
     &  call errquit(pname//': ga_destroy failed m',72, GA_ERR)      
      
      if (dbg>1) call ga_print(g_X)

c     sanity check:

      if (.not.lx2c_u2c) then
        if (g_X.ne.g_x2c) call errquit(
     &    pname//': g_X/g_x2c inconsistency 2',g_X-g_x2c, UNKNOWN_ERR)
      end if

c     ----------------------------------------------
c     transform x2c operator to contracted basis and
c     store in g_x2c. Otherwise, we are already done
c     ----------------------------------------------

      if (lx2c_u2c) then
        if (dbg>0.and.master) write(luout,*)
     &    pname//': u2c trafo of result'
        if (.not. ga_create(MT_DBL,nc,nu,'temp',1,1,g_temp))
     &    call errquit(pname//': ga_create failed temp',77,
     &    GA_ERR)
        if (g_X.eq.g_x2c) call errquit(
     &    pname//': g_X/g_x2c inconsistency 3',g_X-g_x2c, UNKNOWN_ERR)
        
        call ga_dgemm('t','n',nc,nu,nu,1.0d0,
     &    g_u2c,g_X,0.0d0,g_temp)
c       note: if lx2c_u2c.eq..T. then g_x2c and g_X are not the 
c       same. Hence, g_x2c is over-written here
        call ga_dgemm('n','n',nc,nc,nu,1.0d0,
     &    g_temp,g_u2c, 0.0d0,g_x2c)

        if (.not. ga_destroy(g_temp))
     &    call errquit(pname//': ga_destroy failed temp',77,GA_ERR)
      end if                    ! lx2c_u2c

      call ga_symmetrize(g_x2c)

c     de-allocate g_X, if needed, and related arrays

      if (lx2c_u2c) then
        if (.not. ga_destroy(g_X))
     &    call errquit(pname//': ga_destroy failed g_X',61,GA_ERR)

        if (.not.lx2c_keepU) then
          if (.not. ga_destroy(g_u2c)) call errquit
     &      (pname//': ga_destroy failed u2c',51,GA_ERR)
        end if
      end if ! lx2c_u2c
      if (lx2c_c2u .and. .not.lx2c_keepU) then
        if (.not. ga_destroy(g_c2u)) call errquit
     &    (pname//': ga_destroy failed c2u',54,GA_ERR)
      end if ! lx2c_c2u

c     restore model alpha and beta densities:

      if (ipol.gt.1)
     &  call ga_dadd(1.d0, g_dens(1), -1.d0, g_dens(2), g_dens(1))


c     --------
c     all done
c     --------

      call ga_sync()

      if (dbg>1) then
        call ga_print(g_x2c)
      end if

      oskel = oskel_sav ! restore oskel value
      
      if (dbg>0.and.master) write(luout,*) 'leaving x2c_1e_scalar'

c     =================================================================

      return
      end

