      logical function dft_scf(rtdb,
     $     Etold,Enuc,iVcoul_opt,iVxc_opt,
     $     iter,
     $     g_dens,g_movecs,g_vxc,
     $     g_fock,g_svecs,isvals,
     $     g_cdinv,g_xcinv,g_2ceri)
c     
C     $Id: dft_scf.F,v 1.61 1997-02-27 02:52:48 gutowski Exp $
c     
      implicit none
      integer rtdb              ! [input]
      double precision Etold,Enuc
      integer iVcoul_opt
      integer iVxc_opt
      integer iter
      integer g_dens(*),g_movecs(*),g_vxc(*),
     $     g_fock,g_svecs,isvals,
     $     g_cdinv,g_xcinv,g_2ceri, g_scr
c     
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
c     
#if defined(DFT_TIMER)
#include "timer.fh"
#endif
      Logical movecs_write, movecs_converged
      External movecs_write, movecs_converged
      logical oprint,oprint_sym,oadapt
      double precision Exc(2), rms(2), derr(2), n3c_dbl_r, availm_r
      double precision energy
      integer noc(2),nmo(2),icall(2)
      integer n3c_dbl,n3c_int,n_batch
      integer avail, availm
      integer iwhat_max
      integer l_3cwhat,k_3cwhat,l_3cERI,k_3cERI
      integer dft_n3cint, n_semi_bufs
      external dft_n3cint
      double precision dft_n3cdbl
      external dft_n3cdbl
      Integer l_eval
      integer k_eval(2)
      Double Precision rmoments(3) ! dipole
      integer lmax              ! CMM
      integer finest            ! CMM
      integer natoms,nTotEl
      integer l_occ,k_occ
      integer i, j, iw
      integer me,nproc
      integer g_tmp,g_fockt,g_s,g_wght,g_xyz
      integer nheap,nsets,nstack
      integer ispin,nxang,idone
      integer nexc
      integer iswitc
      integer itol_max,iaoacc_max
      integer mxshbf_ao, nao2_max, batch_max, nshells_ao, icount
      integer idmtot,nshells,nprint
      integer leneval,lcd_coef,icd_coef
      double precision ecoul,ecore,aaz
      double precision pp,delta
      double precision anucl_charg,anel
      double precision anoca,anocb,onempp
      double precision etnew,tol2e
c     convergence declarations
      double precision rlshift_input,rlshift_def
      integer ndamp_input,ndamp_def
      character*10 start_10c, NA_10c, asap_10c
      character*10 strng1,strng2,strng3,strng4,strng5,strng6
c     
c     Note, damping, levelshifting, and diising logicals
c     are used to turn on/off these procedures per
c     iteration.  The alternative logicals nodamping, 
c     nolevelshifting, and nodiis are specified and held
c     for the entire convergence sequence.
c     
      logical diising, damping, levelshifting
      logical keep_damp_on,keep_levl_on,keep_diis_on
      Logical LResult,IOLGC,mulliken,do_gradients,lSPM
      logical converged,wght_GA
      logical oprint_parm,oprint_conv,oprint_vecs,
     $     oprint_eval ,oprint_syma,oprint_time,
     $     oprint_info,oprint_tol,oprint_final_vecs,
     $     oprint_energy_step, oprint_intermediate_fock
      integer itol_min,iAOacc_min
      double precision zero,onem,one
      parameter(zero=0.d0,one=1.d0,onem=-one)
      parameter(itol_min=6,iAOacc_min=15)
      double precision dft_dencvg,dabssum,dft_time
      external dft_dencvg,dabssum
      double precision homo, lumo, homo_lumo_gap
      integer stat1, stat2, l_ir, k_ir
      logical status, last_time, last_time_vecs, last_time_energy

      call ecce_print_module_entry('dft')

      dft_scf = .false.
c-----------------------------------------------------------------------
      oprint = util_print('information',print_low)
      oprint_info = util_print('common',print_debug)
      oprint_parm = util_print('parameters', print_default)
      oprint_conv = util_print('convergence', print_default)
      oprint_vecs = util_print('intermediate vectors', print_high)
      oprint_eval = util_print('intermediate evals', print_high)
      oprint_syma = util_print('interm vector symm',
     $     print_high)
      oprint_time = util_print('dft timings',print_high)
      oprint_tol = util_print('screening parameters',print_high)
      oprint_energy_step = util_print('intermediate energy info',
     $     print_high)
      oprint_intermediate_fock = util_print('intermediate fock matrix',
     $     print_high)
      oprint_final_vecs = util_print('final vectors',print_high)
C     
C     oadapt true to get symm analysis
C     
      oadapt=.true.
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c**** (from NWdft data structures, not the nwchem objects)
c     
      if ( .not. bas_numbf(ao_bas_han,nbf_ao) ) then
         call errquit('Exiting from dft',1)
      endif 
      if (.not. geom_ncent(geom, natoms))
     $     call errquit('dft_scf: geom_ncent failed',73)      
      if (.not. geom_nuc_charge(geom, anucl_charg))
     $     call errquit('scf: geom_nuc_charge failed', 0)
c     
      anel = int(anucl_charg) - rcharge
C     
C     determine pattern of orbitals' occupancy
C     
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     $     call errquit('dftscf: rtdb_put of noc failed', 0) 
      if(.not. MA_Push_Get(
     $     MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ))
     $     call errquit('dft: failed to alloc',999)
      
      if( ipol.eq.1 )then
c     
c     noc(1) = No. of electrons alpha
c     
         nTotEl = 2*noc(1)
         nmo(1)=nbf_ao
      else
c     
         nTotEl = noc(1) + noc(2)
         nmo(1)=nbf_ao
         nmo(2)=nbf_ao
      end if
      anoca=noc(1)
      anocb=noc(2)
      if ( ipol .eq. 1 ) then
c**** RHF occupations
         call dfill(nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 2.0d0
         enddo
      else  
c**** UHF occupations
         call dfill(2*nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
         do i = nbf_ao+1, nbf_ao+noc(2)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
      endif 
c     
      wght_GA=.false.
c     
c     Determine whether to fit the electronic charge density.
c     
      CDFIT = .FALSE.
      if( iVcoul_opt.eq.1 ) CDFIT = .TRUE.
      if(CDFIT) then
         if ( .not. bas_numbf(CD_bas_han,nbf_cd) ) then
            call errquit('Exiting from dft',1)
         endif 
      endif
      XCFIT = .FALSE.
      if( iVxc_opt.eq.1 ) XCFIT = .TRUE.
c     
c     Define various constants.
c     
      iDMtot = (ipol*(ipol+1))/2
c     
      itol_max = itol2e
      iaoacc_max = iaoacc
      if(CDFIT) then
c     
c     Determine how big a buffer can be allocated to 3-center
c     2e- integrals.
c     
c     - amount needed for all incore:
c     
         n3c_dbl_r = dft_n3cdbl(AO_bas_han,CD_bas_han,itol_max)
         n3c_int = dft_n3cint(AO_bas_han,CD_bas_han,itol_max)
c     
c     - (minimum) amount available on stack (all nodes):
c     
         call ga_sync
         avail = MA_inquire_stack(mt_dbl)
         call ga_igop(msg_min_stack_avail, avail, 1, 'min')
c     
c     - subtract off a few bits for later use, e.g., quadrature
c     
         availm = avail - ipol*(nbf_ao + nbf_cd + nbf_xc) - 125000
c     
         availm_r = dble(availm)
         availm_r = min(availm_r,n3c_dbl_r)
         n3c_dbl=nint(availm_r)
         if(me.eq.0 .and. oprint_parm)
     $        write(LuOut,3228)avail, availm, n3c_dbl
c     
c     Loops are parallelized over the products of AO shells
c     (check for zero ... must be at least 1).
c     
         if( .not. bas_nbf_cn_max(AO_bas_han,mxshbf_ao) )then
            call errquit('Exiting in dft_scf.',1)
         end if
         if( .not. bas_numcont(AO_bas_han,nshells_ao) )then
            call errquit('Exiting in dft_scf.',2)
         end if
c     
         icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
         icount = max(icount,1)
         if(nproc.eq.1)icount = 1
c     
         nao2_max = mxshbf_ao**2
c
c        compute size of largest block of integrals computed
c
         batch_max = nao2_max*nbf_cd*icount
c     
         if((n3c_dbl.lt.batch_max).or.direct)then
            incore = .false.
            n3c_dbl = 1
            n3c_int = 1
         endif
         if(availm_r.lt.n3c_dbl_r)then
c
c          cannot all be done incore, compute number of buffers 
c          required if semi-direct (I/O)
c
           n_semi_bufs = nint(n3c_dbl_r/availm_r) + 1
         endif
         if( .not.MA_Push_Get(
     $        MT_Dbl,n3c_dbl,'3c ERI space',l_3cERI,k_3cERI))
     $        call errquit('dft_scf: push_get failed', 12)
c     
         if( .not.MA_Push_Get(
     $        MT_int,n3c_int,'3c what space',l_3cwhat,k_3cwhat))
     $        call errquit('dftscf:push_get failed', 13)
c     
         if(incore)then
            if(me.eq.0 .and. oprint_parm)
     $           write(LuOut,3229)n3c_dbl*1.d-6
            call dft_3cincor(rtdb, n_batch, n3c_int,
     $           int_mb(k_3cwhat), dbl_mb(k_3cERI), n3c_dbl,
     $           iwhat_max, n_semi_bufs)
         else
            if(me.eq.0 .and. oprint_parm)write(LuOut,3230)
         endif
      endif
 3228 format(10x,'Minimum dble words available (all nodes) is: ',i15,
     $     /,10x,'         This is reduced (for later use) to: ',i15,
     $     /,10x,'                   Suggested buffer size is: ',i15)
 3229 format(/,10x,f6.3,' MW buffer allocated for incore 3-center '
     $     /,10x,'2e- integral storage on stack. ')
 3230 format(/,10x,'Incore memory use for 3-center 2e- integrals is ',
     $     'turned off. ')
c     
c     set initial coulomb acc to 1e-5
c     
      iswitc=0
      if(ltight) then
         itol2e=itol_max
         iAOacc=iAOacc_max
         iswitc=1
      else
         itol2e=min(itol_min,itol_max)
         iAOacc=min(iAOacc_min,iAOacc_max)
      endif
c     
      tol2e=10.d0**(-itol_max)
c     
c     get info for int2e_
c     
      call scf_get_fock_param(rtdb, tol2e)
      mulliken=.false.
      if(imull.eq.1) mulliken=.true.
      IOLGC=.TRUE.
      if(noio.eq.1) IOLGC=.FALSE.
c     
c     Energy decomposition switch
c     
      nExc    = idecomp + 1
      Etnew=0.d0
C     
C     start timer
C     
      dft_time=-util_cpusec()
c     
c     SCF energy convergence criterion. 
c     
      call gacrea(geom,ao_bas_han,
     $     nbf_ao,nbf_ao,'ga temp',g_tmp,'atom')
      if(ipol.gt.1) then
         call gacrea(geom,ao_bas_han,
     $        nbf_ao,nbf_ao,'fock tr',g_fockt,'atom')
      endif
c     
c     Set up local convergence parameters
c     
      diising = diis
      damping = damp
      levelshifting = levelshift
      keep_damp_on = .false.
      keep_levl_on = .false.
      keep_diis_on = .false.
      ndamp_input = ndamp
      rlshift_input = rlshift
      ndamp_def = 0
      rlshift_def = 0.0
      ndamp = ndamp_input
      rlshift = rlshift_def
c     
c     Three types of convergence speedup schemes:
c     1) based on differences of total energy less than some threshold,
c     2) based on number of cycles performed,
c     3) DIIS plus levelshifting if homo-lumo gap is small.
c     
c     Energetic defaults should be damping until energy converged to < 10-1, 
c     then damping turned off, followed by level shifting and DIIS 
c     until energy converged to < 10-3, then level shifting turned 
c     off, followed by DIIS only.
c     
      if(nodamping)then
         ncydp = 0
         ndamp = 0
         damping = .false.
      endif
      if(nolevelshifting)then
         ncysh = 0
         rlshift = 0.0
         levelshifting = .false.
      endif
      if(nodiis)then
         ncyds = 0
         diising = .false.
      endif
c     
c     If convergence input based upon #cycles then turn off energy constraints.
c     
      if(ncydp.ne.0)then
         dampon  = -999.9
         dampoff = -999.9
         damping = .true. 
         ndamp = ndamp_input
      endif
      if(ncysh.ne.0)then
         levlon  = -999.9
         levloff = -999.9
         levelshifting = .true.
         rlshift = rlshift_input
      endif
      if(ncyds.ne.0)then
         diison  = -999.9
         diisoff = -999.9
         diising = .true.
      endif
      if( me.eq.0 .and. oprint_parm)then
         write(LuOut,3231)
         write(LuOut,3232)ndamp_input, rlshift_input
         asap_10c  = '  ASAP    '
         start_10c = '  start   '
         NA_10c = '   N/A    '
         if(ncydp.ne.0)then
            strng1 = start_10c
            write(strng4,'(i3,7h iters )')ncydp
         elseif(nodamping)then
            strng1 = NA_10c
            strng4 = NA_10c
         else
            write(strng1,'(d10.2)')dampon
            write(strng4,'(d10.2)')dampoff
         endif
c     
         if(ncysh.ne.0)then
            strng2 = asap_10c
            write(strng5,'(i3,7h iters )')ncysh
         elseif(nolevelshifting)then
            strng2 = NA_10c
            strng5 = NA_10c
         else
            write(strng2,'(d10.2)')levlon
            write(strng5,'(d10.2)')levloff
         endif
c     
         if(ncyds.ne.0)then
            strng3 = start_10c
            write(strng6,'(i3,7h iters )')ncyds
         elseif(nodiis)then
            strng3 = NA_10c
            strng6 = NA_10c
         else
            write(strng3,'(d10.2)')diison
            write(strng6,'(d10.2)')diisoff
         endif
         write(LuOut,3233)strng1,strng2,strng3,strng4,strng5,strng6
         call flush(LuOut)
      end if
 3231 format(/,10x,'Convergence aids based upon iterative change in ',
     $     /,10x,'total energy or number of iterations. ')
 3232 format(/,10x,
     $     10x,'Damping(',i2,'%)  Levelshifting(',f3.1,')       DIIS',/,
     $     10x,8x,15('-'),1x,19('-'),1x,15('-'))
 3233 format(10x,'dE  on:',2x,a10,7x,a10,10x,a10,/,
     $     10x,'dE off:',2x,a10,7x,a10,10x,a10,/)
c     
c     Do not allow levelshifting initially until sure that the
c     transformed Fock matrix will be diagonally dominant.
c     
      levelshifting = .false.
c     
c     Initialize DIIS call counter.
c     
      icall(1) = 0
      icall(2) = 0
c     
c     Begin the SCF cycle.
c**** 
c**** allocate eigenvalue array, including second pointer to beta block
c**** 
      leneval = nbf_ao * ipol
      if(.not.MA_Push_Get(MT_Dbl,leneval,'eval',l_eval,k_eval(1)))
     $     call errquit('dft_scf: cannot allocate eval',0)
      if ( ipol .eq. 2) then  
         k_eval(2) = k_eval(1) + nbf_ao
      endif 
c     
c     Dump DFT parameters to see if they make sense
c     
      if(me.eq.0.and.oprint_info)call dft_dump_info
c     
c     Get initial density.
c     
      call dft_init_guess(rtdb, g_dens, g_movecs, g_fock)
      iter = 1
c     
c     Figure input/output MO vectors (done prev in dft_init_guess)
c     
c     See if we are already converged
c     
      if (iguess .eq. 999) then ! If MO starting guess 
         if (movecs_converged(rtdb, ao_bas_han, 'dft', movecs_in)) then
            if (.not. rtdb_get(rtdb, 'dft:energy', mt_dbl, 1, energy))
     $           call errquit('dft: failed getting converged energy',0)
            if (ga_nodeid().eq.0 .and. oprint) then
               write(6,1011) energy
 1011          format(/'  The DFT is already converged '//,
     $              '         Total DFT energy =', f20.12/)
               call util_flush(6)
            endif
            converged = .true.
            goto 1101
         endif
      endif

c**** 
c**** Top of infinite SCF iteration loop
c**** 
      last_time = .false.
      last_time_vecs = .false.
      last_time_energy = .false.

 1000 continue

      if(me.eq.0.and.oprint_tol)write(LuOut,3234)itol2e,iAOacc,iXCacc
 3234 format(10x,'itol2e=',i2,' iAOacc=',i2,' iXCacc=',i2)

      Ecoul  = ZERO
      Exc(1) = ZERO
      Exc(2) = ZERO
      call ga_zero(g_fock)
C     
C     Accumulate core hamiltonian into Fock matrix
C     compute core energy
C     
      call int_1e_ga(ao_bas_han,ao_bas_han,
     $     g_fock,'kinetic',.false.)
      call int_1e_ga(ao_bas_han,ao_bas_han,
     $     g_fock,'potential',.false.)
      Ecore = ga_ddot(g_dens(1),g_fock)
c     
c     Perform SPM analysis
c     
      LSPM=.true.
      if(.not. rtdb_get(rtdb, 'dft:cmm', mt_int, 1, lmax))
     $     LSPM=.false.
      if(LSPM) then
         lresult=geom_charge_center(geom)
         if (.not. rtdb_get(rtdb, 'dft:finest', mt_int, 1, finest))
     $        finest=4
         if(oprint_time) call dft_tstamp(' Before call to CMM.   ')
         call dft_getspm(geom,lmax,finest,g_dens,ao_bas_han)
         call errquit(' end of CMM ',1234)
      endif

#if defined(DFT_TIMER)
      call dft_second(tp3)
#endif

      if( CDFIT ) then
c     
c     Fit the electron charge density.
c     
         if(oprint_time) 
     $        call dft_tstamp(' Before call to FITCD.   ')
         if(.not.MA_Push_Get(MT_Dbl,nbf_cd,'cd_coef',lcd_coef,
     $        icd_coef))
     $        call errquit('dft_scf: cannot allocate cd_coef',0)
         
         call dft_fitcd( rtdb,ao_bas_han, CD_bas_han,geom,
     $        Dbl_MB(icd_coef),  nbf_cd,
     $        g_cdinv, g_dens, nTotEl, IOLGC,itol2e,
     $        Ecoul,incore,
     $        n_batch,n3c_int,
     $        int_mb(k_3cwhat),dbl_mb(k_3cERI),n3c_dbl,
     $        iwhat_max)
      endif
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(6)='   FitCD  '
      call tcllct(t_pgm(1,6),tp4,tp3)
      call dft_second(tp3)
#endif
      if(oprint_time)
     $     call dft_tstamp(' Before call to GETVCOUL.')
      call dft_getvc(rtdb,geom,ao_bas_han,CD_bas_han,
     $     Dbl_MB(icd_coef),g_tmp,
     $     Ecoul,nbf_cd,iVcoul_opt,
     $     g_2ceri,g_dens,IOLGC,itol2e,incore,
     $     n_batch,n3c_int,
     $     int_mb(k_3cwhat),dbl_mb(k_3cERI),n3c_dbl,
     $     iwhat_max)
c     
c     Add V coul to Fock Matrix
c     
      call ga_dadd(one,g_tmp,one,g_fock,g_fock)
      if(CDFIT)then
         if(.not.ma_pop_stack(lcd_coef))
     $        call errquit('dft_scf: cannot pop stack',0)
      endif
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(7)='   GtVcoul'
      call tcllct(t_pgm(1,7),tp4,tp3)
      call dft_second(tp3)
#endif
      if(ipol.gt.1) then
c     
C     Restore alpha and beta densities.
c     
         call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))
c     
c     Note that g_dens(1) now contains the alpha density
c     matrix and g_dens(2) contains the beta
c     
      endif
c     
c     Compute the XC potential and energy.
c     
      g_vxc(1)=g_tmp
      call ga_zero(g_vxc(1))
      if(ipol.eq.2)call ga_zero(g_vxc(2))
      if(oprint_time)call dft_tstamp(' Before call to GETVXC.  ')
      call xc_getv(rtdb,
     $     Exc,nExc,
     $     iVxc_opt,g_xcinv,
     $     g_dens,g_vxc,IOLGC,
     $     g_wght,g_xyz,wght_GA)
      call ga_sync
      call ga_dadd(1.d0,g_vxc(1),1.d0,g_fock,g_fock)
      if(ipol.gt.1) then
c     
c     Put (-alpha+beta) in g_vxc(2)
c     
         call ga_dadd(-1.d0,g_vxc(1),1.d0,g_vxc(2),g_vxc(2))
      endif
c     
      if(oprint_time)
     $     call dft_tstamp(' End of parallel region. ')
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(8)='   GetVXC '
      call tcllct(t_pgm(1,8),tp4,tp3)
#endif
c     
c     Calculate the total electronic energy.
c     
      if( nExc.eq.1 )then
         Etnew = Ecore + Ecoul + Exc(1)
      else
         Etnew = Ecore + Ecoul + Exc(1) + Exc(2)
      end if
c     
      if(last_time_energy)then
c     
c     If open shell put the total density matrix back in 
c     g_dens(1) and quit.
c     
         if( ipol.eq.2 )then
            call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
         end if
         goto 2000
      endif
      call ga_sync
c     
      delta=abs(etold-etnew)
      call ga_sync
      rms(1) = 0.d0
      rms(2) = 0.d0
      homo_lumo_gap = 100.0d0
      do ispin = 1,ipol
C     
C     Add V xc to Fock Matrix
C     
         if(ispin.eq.2) 
     *        call ga_copy(g_fockt,g_fock)
         if(ispin.eq.2) then
            call ga_dadd(one,g_vxc(ispin),one,g_fock,g_fock)
         endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp3)
#endif
c     
         if(ispin.eq.1.and.ipol.gt.1) call ga_copy(g_fock,g_fockt)
         if(diising)then
c     
c     DIIS step taken here.
c     
            if(diising) then
               call diis_driver(derr(ispin),ispin,ipol,icall,
     $              nfock,nbf_ao,geom,ao_bas_han,
     $              g_fock,g_dens(ispin),g_tmp,
     $              g_svecs,isvals,diising,nodiis)
            endif
         endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp4)
         pgmsect(9)='   DIIS   '
         call tcllct(t_pgm(1,9),tp4,tp3)
#endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp3)
#endif
         call gacrea(geom,ao_bas_han,
     $        nbf_ao,nbf_ao,'ga scr',g_scr,'atom')
c     
c     Put s-1/2 in g_scr.
c     
         iw=2
         call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     $        g_tmp,nbf_ao,iw)
c     
c     Transform Fock matrix.
c     
         call ga_dgemm('T','N',nbf_ao,nbf_ao,nbf_ao,one,
     $        g_scr,g_fock,zero,g_tmp)
         call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     $        g_tmp,g_scr,zero,g_fock)
c     
c     Level shifting is implemented here (similarity 
c     transformation before standard eigensolver).  Note,
c     levelshifting is appropriate once a transformation
c     is available which makes the resulting Fock matrix 
c     diagonally dominant, e.g., in an approximate MO basis.  
c     Also note, there are many matrix multiplies with S^+-1/2 
c     which are redundant if one is sure that the former basis
c     is orthonormal.
c     
         if(levelshifting) then
c     
c     Transform Fock matrix to diagonal dominant form with
c     old MOs.
c     
c     Build S^+1/2 * old movecs (as a transform).
c     
            iw=3
            call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     $           g_tmp,nbf_ao,iw)
            call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     $           g_scr,g_movecs(ispin),zero,g_tmp)
            call ga_copy(g_tmp, g_movecs(ispin))
c     
c     Transform Fock matrix         
c     
            call ga_dgemm('T','N',nbf_ao,nbf_ao,nbf_ao,one,
     $           g_movecs(ispin),g_fock,zero,g_scr)
            call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     $           g_scr,g_movecs(ispin),zero,g_fock)
c     
c     Levelshift here.
c     
            do j=noc(ispin)+1+me,nbf_ao,nproc
               call ga_get(g_fock,j,j,j,j,AAZ,1)
               AAZ=AAZ+rlshift
               call ga_put(g_fock,j,j,j,j,AAZ,1)
            enddo
c     write(*,*)' Fock matrix after shift '
c     call ga_print(g_fock)
         else
            rlshift = 0.0
         endif
         call ga_sync
c     
c     Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c     
         if(oprint_intermediate_fock)then     
            write(*,*)' Fock matrix before diagonalization: '
            call ga_print(g_fock)
         endif
#if defined(PARALLEL_DIAG)
         call ga_diag_std(g_fock,g_tmp,Dbl_MB(k_eval(ispin)))
#else 
         call ga_diag_std_seq(g_fock,g_tmp,Dbl_MB(k_eval(ispin)))
#endif
c     
c     Check HOMO/LUMO gap.
c     
         homo = Dbl_MB(k_eval(ispin)+noc(ispin)-1)
         lumo = Dbl_MB(k_eval(ispin)+noc(ispin))
         homo_lumo_gap = min(homo_lumo_gap,(lumo - homo))
c     write(*,*)' 1:homo,lumo,homo_lumo_gap: ',homo,lumo,homo_lumo_gap
c     
c     If levelshifting then tidy up.
c     
         if(levelshifting)then
c     
c     Back transform eigenvectors with S^+1/2 * old movecs.
c     
            call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     $           g_movecs(ispin),g_tmp,zero,g_scr)
            call ga_copy(g_scr, g_tmp)
c     
c     Put S^-1/2 back in g_scr (use g_fock as temp scr).
c     
            iw=2
            call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     $           g_fock,nbf_ao,iw)
c     
c     Put old-movecs back.
c     
            call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     $           g_scr,g_movecs(ispin),zero,g_fock)
            call ga_copy(g_fock, g_movecs(ispin))
            lumo = lumo - rlshift
            homo_lumo_gap = min(homo_lumo_gap,(lumo - homo))
c     write(*,*)' 2:lumo,homo_lumo_gap: ',lumo,homo_lumo_gap
         endif
c     
c     Back-transform eigenvectors with S^-1/2.
c     
         call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     $        g_scr,g_tmp,zero,g_fock)
         call gadest(g_scr)
         call ga_sync
c     
c     Keep orbital ordering according to principle
c     of maximum overlap with previous iteration.
c     Note, do this after MOs have settled down a
c     few iterations.
c     
         if(lmaxov.and.(delta.lt.0.01)) then
            call dft_mxovl(ao_bas_han, dbl_mb(k_eval(ispin)),
     $           g_fock, g_movecs(ispin),g_tmp)
            homo = -100.0d0
            lumo =  100.0d0
            do i = 1, noc(ispin)
               homo = max(homo,(Dbl_MB(k_eval(ispin)+i-1)))
            enddo
            do i = noc(ispin)+1, nbf
               lumo = min(lumo,(Dbl_MB(k_eval(ispin)+i-1)))
            enddo
         else
            call ga_copy(g_fock,g_movecs(ispin))
         endif
c     
         call ga_sync
#if defined(DFT_TIMER)
         call dft_second(tp4)
         pgmsect(10)='   Diag   '
         call tcllct(t_pgm(1,10),tp4,tp3)
#endif
c     
c     Save previous density for convergence check.
c     
         call ga_copy(g_dens(ispin),g_tmp)
c     
c     Form a new density matrix.
c     
         call dft_densm(g_dens(ispin),g_movecs(ispin),
     $        nbf_ao,ipol,noc(ispin))
c     
         call ga_sync
c     
c     write(*,*)' New ',ispin,' density '
c     call ga_print(g_dens(ispin))
c     
c     Check convergence on Density.
c     
         rms(ispin) = dft_dencvg(g_dens(ispin),g_tmp,nbf_ao)
         call ga_sync
c     
c     Damping implemented here.
c     
         if(damping)then
            pp=dble(ndamp)*1.d-2
            onempp=1.0d0 - pp
            call ga_dadd(pp,g_tmp,onempp,g_dens(ispin),g_dens(ispin))
         else
            ndamp = 0
         endif
         call ga_sync
      enddo                     ! end loop over ispin
c     
      if(oprint_conv.and.iter.eq.1.and.me.eq.0)then
         nheap = MA_Inquire_Heap(MT_Dbl)
         nstack = MA_Inquire_Stack(MT_Dbl)
         write(LuOut,21)
         write(LuOut,'(10x,a,f10.2,i20)')
     $        ' Heap Space remaining (MW):  ',dble(nheap)*1.D-06,nheap
         write(LuOut,'(10x,a,f10.2,i20)')
     $        'Stack Space remaining (MW):  ',dble(nstack)*1.D-06,nstack
         call flush(LuOut)
         write(LuOut,1)
      endif
 21   format(/,10x,' Memory utilization after 1st SCF pass: ')
    1 format(/,
     $     1x,'  convergence    iter        energy       DeltaE   ',
     $     'RMS-Dens  Diis-err    time'/
     $     1x,'---------------- ----- ----------------- --------- ',
     $     '--------- ---------  ------')
      if(.not.last_time.and.oprint_conv.and.me.eq.0) then
         if(diising) then
            write(LuOut,2) ndamp,rlshift,
     $           iter, Etnew+Enuc,
     $           -etold+etnew,sqrt(rms(1)),derr(1),util_cpusec()
            if(ipol.eq.2)write(LuOut,3) sqrt(rms(2)),derr(2)
         else
            write(LuOut,22) ndamp,rlshift,
     $           iter, Etnew+Enuc,
     $           -etold+etnew,sqrt(rms(1)), util_cpusec()
            if(ipol.eq.2)write(LuOut,23) sqrt(rms(2))
         endif
         call util_flush(LuOut)
      endif
    2 format(1x,'d=',i2,',ls=',f3.1,',diis',1x,i5,f18.10,
     $     1p,3d10.2,0p,f8.1)
    3 format(51x,1p,2d10.2)
 22   format(1x,'d=',i2,',ls=',f3.1,6x,i5,f18.10,
     $     1p,2d10.2,10x,0p,f8.1)
 23   format(51x,1p,1d10.2)
      call ga_sync
c     
c     save eigenvectors to movecs file
c     
      if (.not.
     $     movecs_write(rtdb, ao_bas_han, movecs_out, 'dft', 
     $     title,
     $     nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     $     dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     $     call errquit('dft_scf: movec_write failed', 0)
c     
      call ga_sync
c     
      if (me .eq. 0.and.oprint_eval) then
         if (util_print('intermediate evals', print_default)) then
            call util_print_centered(LuOut,'alpha eigenvalues',
     $           20,.true.)
            call output(dbl_mb(k_eval(1)),
     $           1, min(noc(1)+10,nbf_ao),
     $           1, 1, nbf_ao, 1, 1)
            if(ipol.eq.2) then
               write(LuOut,*)
               call util_print_centered(LuOut,'beta eigenvalues',
     $              20,.true.)
               call output(dbl_mb(k_eval(2)),
     $              1, min(noc(1)+10,nbf_ao),
     $              1, 1, nbf_ao, 1, 1)
            endif
            call util_flush(6)
         endif
      endif
      if (oprint_vecs) then
         if (me .eq. 0) then
            write(LuOut,*)
            call util_print_centered(LuOut,
     $           'Intermediate MO vectors',40,.true.)
            write(LuOut,*)
            call util_flush(LuOut)
         end if
         call ga_print(g_movecs)
         if(ipol.eq.2) call ga_print(g_movecs(2))
      end if
      if(oprint_syma) then
         call scf_movecs_sym_adapt(ao_bas_han, g_movecs,
     $        oprint_syma,nbf_ao, '- alpha', .false., 1)
         if(ipol.eq.2) 
     $        call scf_movecs_sym_adapt(ao_bas_han, g_movecs(2), 
     $        oprint_syma, nbf_ao, '- beta', .false., 1)
      endif
      
c     
c     Form the total density matrix.
c     
      if( ipol.eq.2 )then
         call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
      end if
      call ga_sync
c     
c     Check for SCF convergence.
c     
      call ga_sync
      if(.not.last_time)call dft_scfcvg(rms,
     $     Etold,Etnew,Enuc,e_conv,d_conv,g_conv,
     $     ipol,iter,iterations,idone,converged,rtdb)
      if(delta.lt.1.d-4) then
c     
c     Set coulomb acc to max (e.g., input parameter).
c     (note, may also require re-initializing DIIS)
c     
         itol2e=itol_max
         iAOacc=iAOacc_max
         iswitc=iswitc+1
      endif
c     
      if(.not.last_time)iter = iter + 1
c     
c     Check convergence parameters.
c     
      if((delta.lt.dampon.and.delta.gt.dampoff).or.iter.le.ncydp)then
         damping = .true.
         ndamp = ndamp_input
      else
         damping = .false.
         ndamp = ndamp_def
      endif
c     
      if((delta.lt.levlon.and.delta.gt.levloff).or.
     $     (iter.le.ncysh))then
c     write(*,*)' homo_lumo_gap, hl_tol: ',homo_lumo_gap, hl_tol
         if(homo_lumo_gap.lt.hl_tol)then
            levelshifting = .true.
            rlshift = rlshift_input
         else
            levelshifting = .false.
            rlshift = rlshift_def
         endif 
c     write(*,*)' 3: levelshifting = ',levelshifting
      else
         levelshifting = .false.
         rlshift = rlshift_def
      endif
c     write(*,*)' 4: levelshifting = ',levelshifting
c     
      if((delta.lt.diison.and.delta.gt.diisoff).or.
     $     iter.le.ncyds.or.keep_diis_on)then
         diising = .true.
c     
c     Once started, keep DIIS on until diisoff threshold.
c     
         keep_diis_on = .true.
      else
         diising = .false.
      endif
      if(delta.lt.diisoff.or.(ncyds.gt.0.and.iter.gt.ncyds))then
         diising = .false.
         keep_diis_on = .false.
      endif
c     
      if(nodamping)damping = .false.
      if(nolevelshifting)levelshifting = .false.
      if(nodiis)diising = .false.
c     
      Etold = Etnew
c     
c     write(*,*)'homo,lumo,lumo-homo,-hl_tol: ',
c     $           homo,lumo,lumo-homo,-hl_tol
      if ((lumo - homo).lt.-hl_tol.and.lmaxov)then
         lmaxov = .false.
         if(me.eq.0 .and. oprint_conv)write(LuOut,224)homo, lumo
 224     format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     $        /,10x,'Significant orbital reordering with',
     $        ' maximum overlap',
     $        /,10x,'turned ON.  Turning max_ovl OFF.')
      endif
c     
      if(oprint_energy_step.and.me.eq.0)then         
         if(nexc.le.1)then
            write(LuOut,222) etnew+enuc, ecore, Ecoul,Exc(1), enuc, 
     $           util_cpusec()
         else
            write(LuOut,223) etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     $           enuc, util_cpusec()
         endif
      endif
c     
      if( idone.eq.0.or.
     $     (iswitc.lt.2.and.iter.lt.iterations)) go to 1000 ! begin new iteration
      if(idone.eq.1.and.(.not.last_time_vecs))then
         last_time_vecs = .true.
         last_time = .true.
         if(me.eq.0 .and. oprint_conv)write(LuOut,225)
 225     format(/,10x,'Done. Constructing final Fock matrix, orbitals',
     $        ' and orbital energies.',/)
         diising = .false.
         levelshifting = .false.
         damping = .false.
         go to 1000             ! build final orbitals/energies
      endif
      if(idone.eq.1.and.(.not.last_time_energy))then
         last_time_energy = .true.
         last_time = .true.
         if(me.eq.0 .and. oprint_conv)write(LuOut,226)
 226     format(/,10x,'      Constructing final total energies.',/)
         diising = .false.
         levelshifting = .false.
         damping = .false.
         go to 1000             ! build final orbitals/energies
      endif
c     
 2000 continue
c     
      if(me.eq.0.and.oprint)then
         if(.not.converged)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     $           'Calculation failed to converge', 20, .true.)
            write(LuOut,*)
         endif
         dft_time=dft_time+util_cpusec()
         if(nexc.le.1)then
            write(LuOut,222) etnew+enuc, ecore, Ecoul,Exc(1), enuc, 
     $           dft_time
         else
            write(LuOut,223) etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     $           enuc, dft_time
         endif
 222     format(//
     $        '         Total DFT energy =', f20.12/
     $        '      One electron energy =', f20.12/
     $        '           Coulomb energy =', f20.12/
     $        '    Exchange-Corr. energy =', f20.12/
     $        ' Nuclear repulsion energy =', f20.12//
     $        '        Time for solution =', f9.1,'s'//)
 223     format(//
     $        '         Total DFT energy =', f20.12/
     $        '      One electron energy =', f20.12/
     $        '           Coulomb energy =', f20.12/
     $        '          Exchange energy =', f20.12/
     $        '       Correlation energy =', f20.12/
     $        ' Nuclear repulsion energy =', f20.12//
     $        '        Time for solution =', f9.1,'s'//)
         call util_flush(LuOut)
      endif
c     
c     call to Mulliken Pop Ananlysis
c     
      if(mulliken) then
         call gacrea(geom,ao_bas_han,
     $        nbf_ao,nbf_ao,'ovl ',g_s,'atom')
         call ga_zero(g_s)
         call int_1e_ga(ao_bas_han,ao_bas_han,
     $        g_s,'overlap',.false.)
         if(me.eq.0)
     $        call dft_header(
     $        ' Total Density - Mulliken Population Analysis')
         call mull_pop(geom,ao_bas_han,g_dens(1),g_s)
         if(ipol.eq.2) then
c     
c     analysis of spin density
c     
            if(me.eq.0) call dft_header(
     $           ' Spin Density - Mulliken Population Analysis')
            call ga_dadd(one,g_dens(1),-2.d0,g_dens(2),g_dens(2))
            call mull_pop(geom,ao_bas_han,g_dens(2),g_s) 
            call ga_dadd(one,g_dens(1),-1.d0,g_dens(2),g_dens(2))
            call ga_dscal(g_dens(2),0.5d0)
         endif
         call gadest(g_s)
      endif
c**** 
c**** end infinite loop for SCF iterations
c**** 
c**** 

c     
c     Store energy and convergence status ... must store before
c     write movecs since date of insertion is used.
c     
      if (.not. rtdb_put(rtdb, 'dft:energy', MT_DBL, 1, (Etold+Enuc))) 
     *     call errquit('dft: failed to store energy in rtdb', 0)
      if (.not. rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged)) 
     *     call errquit('dft: failed to store converged in rtdb', 0)


c**** (if gradients to be calculated) output energies and eigenvectors
c**** to disk
c**** 
      do_gradients = .true.
      if ( do_gradients ) then
c-----------------------------------------------------------------------
c**** 
c**** Save the mo vectors (this code from ddscf, generalized for
c**** spin-unrestricted case)
c**** 
c**** collect the eigenvectors from global storage; output
c**** to disk, allowing for RHF and UHF cases
c**** 
c**** Should distinguish occupied and virtual here...
c**** 

         if (.not.
     $        movecs_write(rtdb, ao_bas_han, movecs_out, 'dft',
     $        title,
     $        nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     $        dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     $        call errquit('rhf: movec_write failed', 0)
         call ga_sync()
      endif 
c     
c     Shut down DIIS.
c     
      if(icall(1).gt.0) then
         icall(1) = -1
         call diis_driver(derr(1),1,ipol,icall,
     $        nfock,nbf_ao,geom,ao_bas_han,
     $        g_fock,g_dens(1),g_tmp,
     $        g_svecs,isvals,diising,nodiis)
      end if
c     
c     If open shell compute overlap of alpha orbitals with beta orbitals.
c     
      if(ipol.gt.1)then
         call dft_mxspin_ovlp(ao_bas_han, g_movecs(1), g_movecs(2), 
     $        g_tmp)
      endif
c     
      if(ipol.gt.1) then
         call gadest(g_fockt)
      endif
      call gadest(g_tmp)
      if(wght_GA) then
         call gadest(g_wght)
         call gadest(g_xyz)
      endif
c     
*     if (filesize.gt.0 .or. memsize.gt.0) then
*     if (.not. int2e_file_close())
*     $        call errquit('dft_scf: closing aoints?', 0)
*     endif
c     
c-----------------------------------------------------------------------
c
C     Restore alpha and beta densities.
c
      if (ipol .gt. 1)
     $     call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))
c-----------------------------------------------------------------------
C     
C     computation of <S2> for open shell
C     
      if(ipol.gt.1)then

         call dft_s2_value(geom,AO_bas_han,.false.,noc(1),noc(2),
     $        nbf_ao,g_dens(1),g_dens(2))

      endif
c-----------------------------------------------------------------------
c     
c     computation of moments
      
      if (natoms .gt. 1) 
     $     call dft_mpole(rtdb, ao_bas_han, ipol, g_dens(1), g_dens(2))

C     
C     print stolen for uhf.F
C     
      if (util_print('schwarz',print_high).and.(.not.CDFIT)) then
         if ( .NOT. bas_numcont(AO_bas_han, nshells) ) 
     $        Call ErrQuit(
     (        'dft: problem with call to bas_numcont', AO_bas_han)
         call schwarz_print(natoms, nshells)
      endif
C     
      if (me .eq. 0) then
         if (util_print('final evals', print_default)) then
            call util_print_centered(LuOut,'Final alpha eigenvalues',
     $           20,.false.)
            call util_print_centered(LuOut,
     $           '(all occupied plus 10 virtual)',20,.true.)
            call output(dbl_mb(k_eval(1)),
     $           1, min(noc(1)+10,nbf_ao),
     $           1, 1, nbf_ao, 1, 1)
            if(ipol.eq.2) then
               write(LuOut,*)
               call util_print_centered(LuOut,'Final beta eigenvalues',
     $              20,.true.)
               call util_print_centered(LuOut,
     $              '(all occupied plus 10 virtual)',20,.true.)
               call output(dbl_mb(k_eval(2)),
     $              1, min(noc(1)+10,nbf_ao),
     $              1, 1, nbf_ao, 1, 1)
            endif
            call util_flush(6)
         endif
         if (oprint_final_vecs) then
            write(LuOut,*)
            call util_print_centered(
     (           LuOut,'Final MO vectors',40,.true.)
            write(LuOut,*)
            call util_flush(LuOut)
         end if
      end if
      if (oadapt) then
         nprint = min(noc(1)+10,nbf_ao)
         if (util_print('all vector symmetries', print_high)) then
            nprint = nbf_ao
            oprint_sym = util_print('final vector symmetries',
     $           print_default)
            call scf_movecs_sym_adapt(ao_bas_han, g_movecs,
     $           oprint_sym,nprint, '- alpha', .false., 1)
            if(ipol.eq.2) 
     $           call scf_movecs_sym_adapt(ao_bas_han, g_movecs(2), 
     $           oprint_sym, nprint, '- beta', .false., 1)
         endif
      endif
      if (oprint_final_vecs) then
         call ga_print(g_movecs)
         if(ipol.eq.2) call ga_print(g_movecs(2))
      end if
C     
C     ECCE printout (for SC '96)
C     
C     requested?
      stat1 = 0
      call ecce_print_control(stat1, stat2) ! stat2 set to print status
      stat1 = stat2
      call ecce_print_control(stat1, stat2) ! restore previous setting
C     we need the value from process 0
      call ga_brdcst(12321, stat1, ma_sizeof(mt_int,1,mt_byte), 0)
      if ( stat1 .eq. 1 ) then
C     allocate array for irreps
         status = ma_push_get ( mt_dbl, nbf_ao,'dft:irreps', 
     $        l_ir, k_ir )
         if ( .not. status ) 
     $        call errquit ('dft:error allocating ir array', 110 )
         if ( oadapt ) then
            call scf_movecs_sym_adapt( ao_bas_han, g_movecs(1),
     $           .false., nbf_ao, '- alpha', .true., int_mb(k_ir) )
         else
            call ifill ( nbf_ao, 1, int_mb(k_ir), 1 )
         endif
         call movecs_ecce(nbf_ao, nmo, 1, nmo, dbl_mb(k_eval(1)),
     $        dbl_mb(k_occ), int_mb(k_ir), 
     $        g_movecs(1), 'dft', 'alpha')
         if ( nsets.eq.2 ) then ! spin-unrestricted
            if ( oadapt ) then
               call scf_movecs_sym_adapt( ao_bas_han, g_movecs(2),
     $              .false., nbf_ao, '- alpha', .true., int_mb(k_ir) )
            endif
            call movecs_ecce(nbf_ao, nmo, 1, nmo, dbl_mb(k_eval(2)),
     $           dbl_mb(k_occ+nbf_ao), int_mb(k_ir), 
     $           g_movecs(2), 'dft', 'alpha')
         endif
         status = ma_pop_stack ( l_ir )
         call ecce_print1 ('total energy', mt_dbl, (Etold+Enuc), 1 )
         call ecce_print1 ('nuclear repulsion energy', mt_dbl, Enuc, 1 )
         call ecce_print1 ('coulomb energy', mt_dbl, Ecoul, 1 )
         call ecce_print1 ('exchange energy', mt_dbl, Exc(1), 1 )
         if ( nexc.gt. 1) then
            call ecce_print1 ('correlation energy', mt_dbl, Exc(2), 1 )
         endif
      endif
c
 1101 continue                  ! Jump here if already converged
      call fock_2e_tidy(rtdb)
c     
      if ( converged ) then
         call ecce_print_module_exit('dft', 'ok')
      else
         call ecce_print_module_exit('dft', 'failed')
      endif
c**** 
c**** eval deallocation moved here from inside iteration loop
c**** 
      if(.not.ma_pop_stack(l_eval))
     $     call errquit('dft_scf: cannot pop stack',0)
      if(CDFIT) then
         if(.not.ma_pop_stack(l_3cwhat))
     $        call errquit('dft_scf: cannot pop stack',0)
         if(.not.ma_pop_stack(l_3cERI))
     $        call errquit('dft_scf: cannot pop stack',0)
      endif
      if(.not.ma_pop_stack(l_occ))
     $     call errquit('dft_scf: cannot pop stack',0)
c     
c-----------------------------------------------------------------------

      dft_scf = converged

      return
c     
 1111 format(15x,'Core Energy:              ',f20.10)
c     
      end

      subroutine dft_dump_info
      implicit none
      integer i
c      
#include "cdft.fh"
c
c     Logicals
c
      write(*,*)' LOGICALS '
      write(*,*)' nodamping = ', nodamping
      write(*,*)' nolevelshifting = ', nolevelshifting
      write(*,*)' nodiis = ', nodiis
      write(*,*)' ltight = ', ltight
      write(*,*)' lmaxov = ', lmaxov
      write(*,*)' incore = ', incore
      write(*,*)' ldelley = ', ldelley
      write(*,*)' store_wght = ', store_wght
      write(*,*)' XCFIT = ', XCFIT
      write(*,*)' CDFIT = ', CDFIT
c
c     Integers
c
      write(*,*)' INTEGERS '
      write(*,*)' iAOacc = ', iAOacc
      write(*,*)' iXCacc = ', iXCacc
      write(*,*)' irqmax = ', irqmax
      write(*,*)' itol2e = ', itol2e
      write(*,*)' imull = ', imull
      write(*,*)' iangquad = ', iangquad
      write(*,*)' iterations = ', iterations
      write(*,*)' mult = ', mult
      write(*,*)' ipol = ', ipol
      write(*,*)' nfock = ', nfock
      write(*,*)' ncydp = ', ncydp
      write(*,*)' ncysh = ', ncysh
      write(*,*)' ncyds = ', ncyds
      write(*,*)' idirect = ', idirect
      write(*,*)' noio = ', noio
      write(*,*)' ndamp = ', ndamp
      write(*,*)' idecomp = ', idecomp
      write(*,*)' nagrid = ', nagrid
      write(*,*)' nrinc = ', nrinc
      write(*,*)' nrmax = ', nrmax
      write(*,*)' geom = ', geom
      write(*,*)' ncenters = ', ncenters
      write(*,*)' nbf = ', nbf
      write(*,*)' nbf_ao = ', nbf_ao
      write(*,*)' nbf_cd = ', nbf_cd
      write(*,*)' nbf_xc = ', nbf_xc
      write(*,*)' MXqshells = ', MXqshells
      write(*,*)' MXrad = ', MXrad
      write(*,*)' MXagrid = ', MXagrid
      write(*,*)' MXtheta = ', MXtheta
      write(*,*)' MXphi = ', MXphi
      write(*,*)' MXnang = ', MXnang
      write(*,*)' IAC_MAX = ', IAC_MAX
      write(*,*)' MAX_CHUNK = ', MAX_CHUNK
      write(*,*)' nqshells = ', nqshells
      write(*,*)' iqrot = ', iqrot
      write(*,*)' ntypes = ', ntypes
      write(*,*)' nquad_task = ', nquad_task
      write(*,*)' dft_ntags_bsmx = ', dft_ntags_bsmx
      write(*,*)' dft_max_cent = ', dft_max_cent
      write(*,*)' AO_bas_han = ', AO_bas_han
      write(*,*)' CD_bas_han = ', CD_bas_han
      write(*,*)' XC_bas_han = ', XC_bas_han
c
c     Doubles
c
      write(*,*)' DOUBLES '
      write(*,*)' rlshift = ', rlshift
      write(*,*)' dampon = ', dampon
      write(*,*)' dampoff = ', dampoff
      write(*,*)' levlon = ', levlon
      write(*,*)' levloff = ', levloff
      write(*,*)' diison = ', diison
      write(*,*)' diisoff = ', diisoff
      write(*,*)' rcharge = ', rcharge
      write(*,*)' rm = ', rm
      write(*,*)' tol_rho = ', tol_rho
      write(*,*)' radius = ', radius
      write(*,*)' rscale = ', rscale
      write(*,*)' cfac = ', (cfac(i),i=1,15)
      write(*,*)' xfac =', (xfac(i),i=1,15)
c
c     Characters
c
      write(*,*)' CHARACTERS '
      write(*,*)' title = ', title
      write(*,*)' gridtype = ', gridtype
      write(*,*)' xname = ', xname
      write(*,*)' cname = ', cname
c
      return
      end
