      Subroutine dft_scf(rtdb,
     &     Etold,Enuc,iVcoul_opt,iVxc_opt,
     &     iter,
     &     g_dens,g_movecs,g_vxc,
     &     g_fock,g_svecs,isvals,
     &     g_cdinv,g_xcinv,g_2ceri)
c
C$Id: dft_scf.F,v 1.38 1996-10-18 00:53:16 d3h449 Exp $
c
      implicit none
      integer rtdb ! [input]
      double precision Etold,Enuc
      integer iVcoul_opt
      integer iVxc_opt
      integer iter
      integer g_dens(*),g_movecs(*),g_vxc(*),
     &     g_fock,g_svecs,isvals,
     &     g_cdinv,g_xcinv,g_2ceri, g_scr
c      
c**** patches for SC95,  ECCE' parsing
***** #define SC95
c
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "global.fh"
#include "msgids.fh"
c
#if defined(DFT_TIMER)
#include "timer.fh"
#endif
      Logical movecs_write, movecs_read_header, movecs_read
      External movecs_write,movecs_read_header, movecs_read
      logical oprint,oprint_sym,oadapt
      double precision Exc(2), rms(2), derr(2)
      integer noc(2),nmo(2),icall(2)
      integer n3c_dbl,n3c_int,n_batch
      integer avail
      integer l_3cwhat,k_3cwhat,l_3cERI,k_3cERI
      integer dft_n3cdbl,dft_n3cint
      external dft_n3cdbl,dft_n3cint
      Integer l_eval,Type
      integer k_eval(2)
      Character*(nw_max_path_len) movecs_out,movecs_in
      Character*255 name, trans 
      Double Precision rmoments(3)  ! dipole
      integer lmax   ! CMM
      integer finest ! CMM
      integer natoms,nTotEl
      integer l_occ,k_occ
      integer i,ii,j, iw
      integer n_ashift,n_bshift
      integer me,nproc
      integer g_tmp,g_fockt,g_s,g_wght,g_xyz
      integer nheap,nsets,nstack
      integer ispin,nxang,idone
      integer ielem
      integer nexc
      integer l_afshift,n_afshift,i_afshift
      integer l_bfshift,n_bfshift,i_bfshift
      integer iswitc
      integer itol_max,iaoacc_max
      integer idmtot,nshells,nprint
      integer leneval,lcd_coef,icd_coef
      double precision ecoul,rb,ecore,anps,aaz
      double precision pp,delta
      double precision anucl_charg,anel
      double precision anoca,anocb,onempp
      double precision etnew,tol2e,fock_el
c convergence declarations
      double precision rlshift_input,rlshift_def
      integer ndamp_input,ndamp_def
      character*10 start_10c, NA_10c, asap_10c
      character*10 strng1,strng2,strng3,strng4,strng5,strng6
c
c     Note, damping, levelshifting, and diis logicals
c     are used to turn on/off these procedures per
c     iteration.  The alternative logicals nodamping, 
c     nolevelshifting, and nodiis are specified and held
c     for the entire convergence sequence.
c
      logical diis, damping, levelshifting
      logical keep_damp_on,keep_levl_on,keep_diis_on
      Logical LResult,IOLGC,mulliken,do_gradients,lSPM
      logical converged,wght_GA
      logical oprint_parm,oprint_conv,oprint_vecs,
     ,     oprint_eval ,oprint_syma,oprint_time
      integer itol_min,iAOacc_min
      double precision scfcon
      double precision zero,onem,one
      parameter(zero=0.d0,one=1.d0,onem=-one)
      parameter(itol_min=6,iAOacc_min=15)
      double precision dft_dencvg,dabssum,dft_time
      external dft_dencvg,dabssum
#if defined(SC95)
      double precision homo_lumo_gap
      logical movecs_to_logfile
      integer n_homo, n_lumo            ! # of HOMOs, LUMOs to output to stdout
#endif      
c-----------------------------------------------------------------------
      oprint = util_print('information',print_low)
      oprint_parm = util_print('parameters', print_default)
      oprint_conv = util_print('convergence', print_default)
      oprint_vecs = util_print('intermediate vectors', print_high)
      oprint_eval = util_print('intermediate evals', print_high)
      oprint_syma = util_print('intermediate vector symmetries',
     ,     print_high)
      oprint_time=util_print('dft timings',print_high)
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
     &     call errquit('dft_scf: geom_ncent failed',73)      
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
     &     MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ))
     &     call errquit('dft: failed to alloc',999)
      
      if( ipol.eq.1 )then
c       
c       noc(1) = No. of electrons alpha
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
c****   RHF occupations
        call dfill(nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
        do i = 1, noc(1)
          dbl_mb(i-1+k_occ) = 2.0d0
        enddo
      else  
c****   UHF occupations
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
      if(CDFIT.and.incore) then
C
C       try in-core Coulomb calculation
C
        n3c_dbl=dft_n3cdbl(AO_bas_han,CD_bas_han,itol_max)
        n3c_int=dft_n3cint(AO_bas_han,CD_bas_han,itol_max)
C        if(me.eq.0) then
C          write(0,*) ' n3c_dbl ',n3c_dbl,ga_nodeid()
C          write(0,*) ' n3c_int ',n3c_int,ga_nodeid()
C        endif
        avail = MA_inquire_avail(mt_dbl)
C        if (ga_uses_ma()) avail = avail - ga_memory_avail()
        call ga_igop(msg_int2e_file, avail, 1, 'min')
        if(avail.lt.n3c_dbl) then
          incore=.false.
        else
          if(.not.MA_Push_Get(
     &         MT_Dbl,n3c_dbl,'3c ERI space',l_3cERI,k_3cERI))
     &         call errquit('dft_scf: push_get failed', 12)
          avail = MA_inquire_avail(mt_int)
          if(avail.lt.n3c_int) then
            incore=.false.
            if(.not.ma_pop_stack(l_3cERI))
     &        call errquit('dft_scf: cannot pop stack',0)
          else
            if(.not.MA_Push_Get(
     &           MT_int,n3c_int,'3c what space',l_3cwhat,k_3cwhat))
     &           call errquit('dftscf:push_get failed', 13)
            call dft_3cincor(rtdb, geom, AO_bas_han, CD_bas_han, 
     &           nbf_cd, itol_max,n_batch,n3c_int,
     &           int_mb(k_3cwhat),dbl_mb(k_3cERI),n3c_dbl)
          endif
        endif
        if(.not.incore.and.me.eq.0) then
          write (LuOut,*) ' CANNOT HOLD  INTS in core'
          write(LuOut,*) ' add ',n3c_dbl-avail,
     & ' doubles to keep ints in core '
        elseif(incore.and.me.eq.0) then
          write (LuOut,*) '         Coulomb component computed incore.'
        endif
      endif
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
c     check if you want to cache ints to disk
c
      tol2e=10.d0**(-itol_max)
c
c     get info for int2e_
c
      call scf_get_fock_param(rtdb, tol2e*100.)
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
      scfcon = 10.D0**(-iscfcon)
      call gacrea(geom,ao_bas_han,
     &     nbf_ao,nbf_ao,'ga temp',g_tmp,'atom')
      if(ipol.gt.1) then
        call gacrea(geom,ao_bas_han,
     &       nbf_ao,nbf_ao,'fock tr',g_fockt,'atom')
      endif
c
c     Set up local convergence parameters
c
      diis = .false.
      damping = .true.
      levelshifting = .false.
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
c     Energetic defaults should be damping until energy converged to < 10-1, 
c     then damping turned off, followed by level shifting and DIIS 
c     until energy converged to < 10-3, then level shifting turned 
c     off, followed by DIIS only.
c 
      if(nodamping)then
        ncydp = 0
        ndamp = 0
      endif
      if(nolevelshifting)then
        ncysh = 0
        rlshift = 0.0
      endif
      if(nodiis)then
        ncyds = 0
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
        diis = .true.
      endif
      if( me.eq.0 )then
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
     &       /,10x,'total energy or number of iterations. ')
 3232 format(/,10x,
     &'       Damping(',i2,'%)  Levelshifting(',f3.1,')       DIIS',/,

     &10x,'     --------------- ------------------- ---------------')
 3233 format(10x,' on:',2x,a10,7x,a10,10x,a10,/,
     &       10x,'off:',2x,a10,7x,a10,10x,a10,/)
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
c     
      if( me.eq.0 )then
        nheap = MA_Inquire_Heap(MT_Dbl)
        nstack = MA_Inquire_Stack(MT_Dbl)
        write(LuOut,'(10x,a,f10.2,i20)')
     &       ' Heap Space remaining (MW):  ',dble(nheap)*1.D-06,nheap
        write(LuOut,'(10x,a,f10.2,i20)')
     &       'Stack Space remaining (MW):  ',dble(nstack)*1.D-06,nstack
	call flush(LuOut)
      end if
c**** 
c**** Top of infinite SCF iteration loop
c**** 
c**** 
c**** allocate eigenvalue array, including second pointer to beta block
c**** 
      leneval = nbf_ao * ipol
      if(.not.MA_Push_Get(MT_Dbl,leneval,'eval',l_eval,k_eval(1)))
     &  call errquit('dft_scf: cannot allocate eval',0)
      if ( ipol .eq. 2) then  
        k_eval(2) = k_eval(1) + nbf_ao
      endif 
c
c     Dump DFT parameters to see if they make sense
c
c      call dft_dump_info
c
c     Get initial density.
c
      call dft_init_guess(rtdb, g_dens, g_movecs, g_fock)
      iter = 1
c
c     Figure input/output MO vectors 
c     
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $     call util_file_name('movecs', .false.,.false.,movecs_in)
      if (.not. rtdb_cget(rtdb, 'dft:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
            movecs_out = movecs_in
      endif

 1000 continue

      Ecoul  = ZERO
      Exc(1) = ZERO
      Exc(2) = ZERO
      call ga_zero(g_fock)
C       
C     Accumulate core hamiltonian into Fock matrix
C     compute core energy
C       
      call int_1e_ga(ao_bas_han,ao_bas_han,
     &  g_fock,'kinetic',.false.)
      call int_1e_ga(ao_bas_han,ao_bas_han,
     &  g_fock,'potential',.false.)
      Ecore = ga_ddot(g_dens(1),g_fock)
c
c     Perform SPM analysis
c
      LSPM=.true.
      if(.not. rtdb_get(rtdb, 'dft:cmm', mt_int, 1, lmax))
     &  LSPM=.false.
      if(LSPM) then
        lresult=geom_charge_center(geom)
        if (.not. rtdb_get(rtdb, 'dft:finest', mt_int, 1, finest))
     &  finest=4
        if(oprint_time) call dft_tstamp(' Before call to CMM.   ')
        call dft_getspm(geom,lmax,finest,g_dens,ao_bas_han)
        call errquit(' end of CMM ',1234)
      endif

#if defined(DFT_TIMER)
      call dft_second(tp3)
#endif  

      if( CDFIT ) then
c         
c       Fit the electron charge density.
c         
        if(oprint_time) 
     &    call dft_tstamp(' Before call to FITCD.   ')
        if(.not.MA_Push_Get(MT_Dbl,nbf_cd,'cd_coef',lcd_coef,
     &    icd_coef))
     &    call errquit('dft_scf: cannot allocate cd_coef',0)
          
        call dft_fitcd( rtdb,ao_bas_han, CD_bas_han,geom,
     &    Dbl_MB(icd_coef),  nbf_cd,
     &    g_cdinv, g_dens, nTotEl, IOLGC,itol2e,
     &    Ecoul,incore,
     &    n_batch,n3c_int,
     &    int_mb(k_3cwhat),dbl_mb(k_3cERI),n3c_dbl)
      endif
c       
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(6)='   FitCD  '
      call tcllct(t_pgm(1,6),tp4,tp3)
      call dft_second(tp3)
#endif  
      if(oprint_time)
     &  call dft_tstamp(' Before call to GETVCOUL.')
      call dft_getvc(rtdb,geom,ao_bas_han,CD_bas_han,
     &  Dbl_MB(icd_coef),g_tmp,
     &  Ecoul,nbf_cd,iVcoul_opt,
     &  g_2ceri,g_dens,IOLGC,itol2e,incore,
     &  n_batch,n3c_int,
     &  int_mb(k_3cwhat),dbl_mb(k_3cERI),n3c_dbl)
c       
c     Add V coul to Fock Matrix
c       
      call ga_dadd(one,g_tmp,one,g_fock,g_fock)
      if(CDFIT)then
        if(.not.ma_pop_stack(lcd_coef))
     &    call errquit('dft_scf: cannot pop stack',0)
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
C       Restore alpha and beta densities.
c         
        call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))
c
c       Note that g_dens(1) now contains the alpha density
c       matrix and g_dens(2) contains the beta
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
     &  Exc,nExc,
     &  iVxc_opt,g_xcinv,
     &  g_dens,g_vxc,IOLGC,
     &  g_wght,g_xyz,wght_GA)
      call ga_sync
      call ga_dadd(1.d0,g_vxc(1),1.d0,g_fock,g_fock)
      if(ipol.gt.1) then
c         
c       Put (-alpha+beta) in g_vxc(2)
c         
        call ga_dadd(-1.d0,g_vxc(1),1.d0,g_vxc(2),g_vxc(2))
      endif
c       
      if(oprint_time)
     &  call dft_tstamp(' End of parallel region. ')
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
      delta=abs(etold-etnew)
      call ga_sync
      rms(1) = 0.d0
      rms(2) = 0.d0
      do ispin = 1,ipol
C         
C       Add V xc to Fock Matrix
C
        if(ispin.eq.2) 
     *       call ga_copy(g_fockt,g_fock)
        if(ispin.eq.2) then
          call ga_dadd(one,g_vxc(ispin),one,g_fock,g_fock)
        endif
c         
#if defined(DFT_TIMER)
        call dft_second(tp3)
#endif    
c
        if(ispin.eq.1.and.ipol.gt.1) call ga_copy(g_fock,g_fockt)
        if(diis)then
c
c         DIIS step taken here.
c
          if(diis) then
            call diis_driver(derr(ispin),ispin,ipol,icall,
     &        nfock,nbf_ao,geom,ao_bas_han,
     &        g_fock,g_dens(ispin),g_tmp,
     &        g_svecs,isvals)
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
     &    nbf_ao,nbf_ao,'ga scr',g_scr,'atom')
c
c       Put s-1/2 in g_scr.
c
        iw=2
        call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     &    g_tmp,nbf_ao,iw)
c
c       Transform Fock matrix.
c
        call ga_dgemm('T','N',nbf_ao,nbf_ao,nbf_ao,one,
     &     g_scr,g_fock,zero,g_tmp)
        call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &     g_tmp,g_scr,zero,g_fock)
c       
c       Level shifting is implemented here (similarity 
c       transformation before standard eigensolver).  Note,
c       levelshifting is appropriate once a transformation
c       is available which makes the resulting Fock matrix 
c       diagonally dominant, e.g., in an approximate MO basis.  
c       Also note, there are many matrix multiplies with S^+-1/2 
c       which are redundant if one is sure that the former basis
c       is orthonormal.
c       
        if(levelshifting) then
c
c         Transform Fock matrix to diagonal dominant form with
c         old MOs.
c
c         Build S^+1/2 * old movecs (as a transform).
c
          iw=3
          call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     &      g_tmp,nbf_ao,iw)
          call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &      g_scr,g_movecs(ispin),zero,g_tmp)
          call ga_copy(g_tmp, g_movecs(ispin))
c
c         Transform Fock matrix         
c
          call ga_dgemm('T','N',nbf_ao,nbf_ao,nbf_ao,one,
     &      g_movecs(ispin),g_fock,zero,g_scr)
          call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &      g_scr,g_movecs(ispin),zero,g_fock)
c
c         Levelshift here.
c
          do j=noc(ispin)+1+me,nbf_ao,nproc
            call ga_get(g_fock,j,j,j,j,AAZ,1)
            AAZ=AAZ+rlshift
            call ga_put(g_fock,j,j,j,j,AAZ,1)
          enddo
c          write(*,*)' Fock matrix after shift '
c          call ga_print(g_fock)
        else
          rlshift = 0.0
        endif
        call ga_sync
c     
c       Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c     
#if defined(PARALLEL_DIAG)
        call ga_diag_std(g_fock,g_tmp,Dbl_MB(k_eval(ispin)))
#else   
        call ga_diag_std_seq(g_fock,g_tmp,Dbl_MB(k_eval(ispin)))
#endif  
c
c       If levelshifting then tidy up.
c
        if(levelshifting)then
c
c         Back transform eigenvectors with S^+1/2 * old movecs.
c
          call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &      g_movecs(ispin),g_tmp,zero,g_scr)
          call ga_copy(g_scr, g_tmp)
c
c         Put S^-1/2 back in g_scr (use g_fock as temp scr).
c
          iw=2
          call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     &      g_fock,nbf_ao,iw)
c
c         Put old-movecs back.
c
          call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &      g_scr,g_movecs(ispin),zero,g_fock)
          call ga_copy(g_fock, g_movecs(ispin))
        endif
c
c       Back-transform eigenvectors with S^-1/2.
c
        call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &    g_scr,g_tmp,zero,g_fock)
        call gadest(g_scr)
        call ga_sync
c
c       Keep orbital ordering according to principle
c       of maximum overlap with previous iteration.
c       Note, do this after MOs have settled down a
c       few iterations.
c       
        if(lmaxov.and.(delta.lt.0.01)) then
          call dft_mxovl(ao_bas_han, dbl_mb(k_eval(ispin)),
     &      g_fock, g_movecs(ispin),g_tmp)
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
c       Save previous density for convergence check.
c       
        call ga_copy(g_dens(ispin),g_tmp)
c       
c       Form a new density matrix.
c
        call dft_densm(g_dens(ispin),g_movecs(ispin),
     &       nbf_ao,ipol,noc(ispin))
c       
        call ga_sync
c
c        write(*,*)' New ',ispin,' density '
c        call ga_print(g_dens(ispin))
c
c       Damping implemented here.
c
        call ga_sync
        if(damping)then
          pp=ndamp*1.d-2
          onempp=1.0d0 - pp
          call ga_dadd(pp,g_tmp,onempp,g_dens(ispin),g_dens(ispin))
        else
          ndamp = 0
        endif
c
c       Check convergence on Density.
c
        call ga_sync
        rms(ispin) = rms(ispin) +
     &               dft_dencvg(g_dens(ispin),g_tmp,nbf_ao)
        call ga_sync
c
      enddo     ! end loop over ispin
c
      if(oprint.and.iter.eq.1.and.me.eq.0)then
        write(LuOut,1)
      endif
 1    format(/,
     $      1x,'  convergence    iter        energy       DeltaE   ',
     $         'RMS-Dens  Diis-err    time'/
     $      1x,'---------------- ----- ----------------- --------- ',
     $         '--------- ---------  ------')
      if(oprint.and.me.eq.0) then
        if(diis) then
          write(LuOut,2) ndamp,rlshift,
     &         iter, Etnew+Enuc,
     &         -etold+etnew,sqrt(rms(1)),derr(1),util_cpusec()
          if(ipol.eq.2)write(LuOut,3) sqrt(rms(2)),derr(2)
        else
          write(LuOut,22) ndamp,rlshift,
     &         iter, Etnew+Enuc,
     &         -etold+etnew,sqrt(rms(1)), util_cpusec()
          if(ipol.eq.2)write(LuOut,23) sqrt(rms(2))
        endif
        call util_flush(LuOut)
      endif
    2 format(1x,'d=',i2,',ls=',f3.1,',diis',1x,i5,f18.10,
     &       1p,3d10.2,0p,f8.1)
    3 format(51x,1p,2d10.2)
   22 format(1x,'d=',i2,',ls=',f3.1,6x,i5,f18.10,
     &       1p,2d10.2,10x,0p,f8.1)
   23 format(51x,1p,1d10.2)
      call ga_sync
c
c     save eigenvectors to movecs file
c
      if (.not.
     &     movecs_write(movecs_out, 'Kohn-Sham orbitals',
     &     'ao basis',
     &     nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &     dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     &     call errquit('dft_scf: movec_write failed', 0)
c
      call ga_sync
c
      if (me .eq. 0.and.oprint_eval) then
        if (util_print('intermediate evals', print_default)) then
          call util_print_centered(LuOut,'alpha eigenvalues',
     $         20,.true.)
          call output(dbl_mb(k_eval(1)),
     $         1, min(noc(1)+10,nbf_ao),
     $         1, 1, nbf_ao, 1, 1)
          if(ipol.eq.2) then
            write(LuOut,*)
            call util_print_centered(LuOut,'beta eigenvalues',
     $           20,.true.)
            call output(dbl_mb(k_eval(2)),
     $           1, min(noc(1)+10,nbf_ao),
     $           1, 1, nbf_ao, 1, 1)
          endif
          call util_flush(6)
        endif
      endif
      if (oprint_vecs) then
        if (me .eq. 0) then
          write(LuOut,*)
          call util_print_centered(LuOut,
     $         'Intermediate MO vectors',40,.true.)
          write(LuOut,*)
          call util_flush(LuOut)
        end if
        call ga_print(g_movecs)
        if(ipol.eq.2) call ga_print(g_movecs(2))
      end if
      if(oprint_syma) then
        call scf_movecs_sym_adapt(ao_bas_han, g_movecs,
     $       oprint_syma,nbf_ao, '- alpha', .false., 1)
        if(ipol.eq.2) 
     $       call scf_movecs_sym_adapt(ao_bas_han, g_movecs(2), 
     $       oprint_syma, nbf_ao, '- beta', .false., 1)
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
      if(itrscf.eq.0) then
c       
        call ga_sync
c       
c       
c       Compute the contribution to the total energy from the core Hamiltonian.
c       
        Ecore = ga_ddot(g_dens(1),g_fock)
        if(oprint_conv.and.me.eq.0) write(LuOut,1111) Ecore
      else
        call ga_sync
        call dft_scfcvg(rms,
     &       Etold,Etnew,Enuc,
     &       scfcon,ipol,igcon,iter,itrscf,idone,converged,rtdb)
        if(delta.lt.1.d-4) then
c
c         Set coulomb acc to max (e.g., input parameter).
c         (note, may also require re-initializing DIIS)
c
          itol2e=itol_max
          iAOacc=iAOacc_max
          iswitc=iswitc+1
        endif
c
        iter = iter + 1
c
c       Check convergence parameters.
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
     &     (iter.le.ncysh))then
          levelshifting = .true.
          rlshift = rlshift_input
        else
          levelshifting = .false.
          rlshift = rlshift_def
        endif
c
        if((delta.lt.diison.and.delta.gt.diisoff).or.
     &    iter.le.ncyds.or.keep_diis_on)then
          diis = .true.
c
c         Once started, keep DIIS on until diisoff threshold.
c
          keep_diis_on = .true.
        else
          diis = .false.
        endif
        if(delta.lt.diisoff.or.(ncyds.gt.0.and.iter.gt.ncyds))then
          diis = .false.
          keep_diis_on = .false.
        endif
c
        if(nodamping)damping = .false.
        if(nolevelshifting)levelshifting = .false.
        if(nodiis)diis = .false.
c
        Etold = Etnew
c
        if( idone.eq.0.or.
     &   (iswitc.lt.2.and.iter.lt.itrscf)) go to 1000 ! begin new iteration
c
        if(me.eq.0.and.oprint)then
          if(.not.converged)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     $         'Calculation failed to converge', 20, .true.)
            write(LuOut,*)
          endif
          dft_time=dft_time+util_cpusec()
          if(nexc.le.1)then
            write(LuOut,222) etnew+enuc, ecore, Ecoul,Exc(1), enuc, 
     $                       dft_time
          else
            write(LuOut,223) etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     $                        enuc, dft_time
          endif
 222      format(//
     $        '          Final results '/
     $        '       ------------------ '//
     $        '         Total DFT energy =', f20.12/
     $        '      One electron energy =', f20.12/
     $        '           Coulomb energy =', f20.12/
     $        '    Exchange-Corr. energy =', f20.12/
     $        ' Nuclear repulsion energy =', f20.12//
     $        '        Time for solution =', f9.1,'s'//)
 223      format(//
     $        '          Final results '/
     $        '       ------------------ '//
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
c       call to Mulliken Pop Ananlysis
c     
        if(mulliken) then
          call gacrea(geom,ao_bas_han,
     &         nbf_ao,nbf_ao,'ovl ',g_s,'atom')
          call ga_zero(g_s)
          call int_1e_ga(ao_bas_han,ao_bas_han,
     &         g_s,'overlap',.false.)
          if(me.eq.0)
     &      call dft_header(
     &       ' Total Density - Mulliken Population Analysis')
          call mull_pop(geom,ao_bas_han,g_dens(1),g_s)
          if(ipol.eq.2) then
c         
c           analysis of spin density
c         
            if(me.eq.0) call dft_header(
     &        ' Spin Density - Mulliken Population Analysis')
            call ga_dadd(one,g_dens(1),-2.d0,g_dens(2),g_dens(2))
            call mull_pop(geom,ao_bas_han,g_dens(2),g_s) 
            call ga_dadd(one,g_dens(1),-1.d0,g_dens(2),g_dens(2))
            call ga_dscal(g_dens(2),0.5d0)
          endif
          call gadest(g_s)
        endif
c****   
c****   end infinite loop for SCF iterations
c****   
c****   
c****   (if gradients to be calculated) output energies and eigenvectors
c****   to disk
c****   
        do_gradients = .true.
        if ( do_gradients ) then
c-----------------------------------------------------------------------
c****     
c****     Save the mo vectors (this code from ddscf, generalized for
c****     spin-unrestricted case)
c****     
c****     collect the eigenvectors from global storage; output
c****     to disk, allowing for RHF and UHF cases
c****     
c****     logical function movecs_write(filename, title, basis_name, 
c****     $     nbf_ao, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
c****     occ is an occupation vector (not used here)
c****     
c****     Should distinguish occupied and virtual here...
c****     
          if (.not.
     &         movecs_write(movecs_out, 'Kohn-Sham orbitals',
     &         'ao basis',
     &         nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &         dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     &         call errquit('rhf: movec_write failed', 0)
          call ga_sync()
          if(me.eq.0) then
#if defined(SC95)
            homo_lumo_gap = min( dbl_mb( k_eval(1)+noc(1)),
     &                           dbl_mb( k_eval(2)+noc(2) ) )
     &           - max( dbl_mb( k_eval(1)+noc(1)-1),
     &           dbl_mb( k_eval(2)+noc(2)-1 ) )
            write (LuOut,*) 
            write (LuOut,'(a,e14.4)')
     &           'HOMO/LUMO Gap (Hartrees): ', homo_lumo_gap
            write (LuOut,*) 
            write (LuOut,*) 'Alpha-Occupied Eigenvalues: '
            write (LuOut,'(5x,4(e14.4))')
     &           (dbl_mb(k_eval(1)+ii),ii=0,noc(1)-1)
            write (LuOut,*) 
            write (LuOut,*) 'Alpha-Virtual Eigenvalues: '
            write (LuOut,'(5x,4(e14.4))')
     &           (dbl_mb(k_eval(1)+ii),ii=noc(1),nbf_ao-1)
            write (LuOut,*) 
            if ( ipol.eq.2 ) then
              write (LuOut,*) 'Beta-Occupied Eigenvalues: '
              write (LuOut,'(5x,4(e14.4))')
     &             (dbl_mb(k_eval(2)+ii),ii=0,noc(2)-1)
              write (LuOut,*) 
              write (LuOut,*) 'Beta-Virtual Eigenvalues: '
              write (LuOut,'(5x,4(e14.4))')
     &             (dbl_mb(k_eval(2)+ii),ii=noc(2),nbf_ao-1)
              write (LuOut,*) 
            endif 
#endif
          endif

        endif 
c-----------------------------------------------------------------------
#if defined(SC95)
c****   number of HOMOs and LUMOs to output to logfile, if any
        if ( converged ) then

          if (.not. rtdb_get(rtdb, 'dft:n_HOMO_to_logfile',
     &         mt_int, 1, n_homo)) n_homo = 0
          if (.not. rtdb_get(rtdb, 'dft:n_LUMO_to_logfile',
     &         mt_int, 1, n_lumo)) n_lumo = 0
c****     check values
          if ( nsets.eq.2 ) then ! spin-unrestricted
            n_homo = min( n_homo, noc(1), noc(2) )
            n_lumo = min( n_lumo, nbf_ao-noc(1), nbf_ao-noc(2) )
          else  ! ! spin-restricted
            n_homo = min( n_homo, noc(1) )
            n_lumo = min( n_lumo, nbf_ao-noc(1) )
          endif 

          if ( n_homo.gt.0 .or. n_lumo.gt.0 ) then ! yes, output MO coeff.

            if (.not.
     &           movecs_to_logfile( rtdb, n_homo, n_lumo,
     &           nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &           dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     &           call errquit('dft_scf: movecs_to_logfile failed', 2)
          endif 
        endif 
#endif
c         
c       Shut down DIIS.
c         
        if(icall(1).gt.0) then
          icall(1) = -1
          call diis_driver(derr(1),1,ipol,icall,
     &         nfock,nbf_ao,geom,ao_bas_han,
     &         g_fock,g_dens(1),g_tmp,
     &         g_svecs,isvals)
        end if
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
        call fock_2e_tidy(rtdb)
c
*      if (filesize.gt.0 .or. memsize.gt.0) then
*         if (.not. int2e_file_close())
*     $        call errquit('dft_scf: closing aoints?', 0)
*      endif
c
c-----------------------------------------------------------------------
c
c     computation of dipole moments
c
      
      if( .not. bas_high_angular(AO_bas_han,NXang) )
     &     call errquit('exiting in dft_scf.',2)
      if(natoms.gt.1.and.NXang.gt.1) 
     &     call dft_mpole( geom, AO_bas_han, AO_bas_han,1, 
     &     g_dens, ipol, rmoments, 3)

C
C     computation of <S2> for open shell
C
      if(ipol.gt.1)then

C         Restore alpha and beta densities.
          call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))

          call dft_s2_value(geom,AO_bas_han,.false.,noc(1),noc(2),
     &     nbf_ao,g_dens(1),g_dens(2))

        endif
c**** 
c**** this is the final section of the ddscf routine rhf()
c**** have substituted SCF with DFT
c**** 
c     Store energy and convergence status ... not sure if better to
c     do this here or have each wavefunction type do it
c     
      if (.not. bas_name(AO_bas_han, name, trans))
     $     call errquit('dft_scf: bas_name?', 0)
      if (name .eq. ' ') name = 'ao basis'

      if (.not. rtdb_put(rtdb, 'dft:energy', MT_DBL, 1, (Etold+Enuc))) 
     *     call errquit('dft: failed to store energy in rtdb', 0)
      if (.not. rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged)) 
     *     call errquit('dft: failed to store converged in rtdb', 0)
C
C     print stolen for uhf.F
C
      if (util_print('schwarz',print_high).and.(.not.CDFIT)) then
        if ( .NOT. bas_numcont(AO_bas_han, nshells) ) 
     $       Call ErrQuit(
     (       'dft: problem with call to bas_numcont', AO_bas_han)
         call schwarz_print(natoms, nshells)
      endif
C
      if (me .eq. 0) then
        if (util_print('final evals', print_default)) then
          call util_print_centered(LuOut,'Final alpha eigenvalues',
     $         20,.false.)
          call util_print_centered(LuOut,
     $         '(all occupied plus 10 virtual)',20,.true.)
          call output(dbl_mb(k_eval(1)),
     $         1, min(noc(1)+10,nbf_ao),
     $         1, 1, nbf_ao, 1, 1)
          if(ipol.eq.2) then
            write(LuOut,*)
            call util_print_centered(LuOut,'Final beta eigenvalues',
     $           20,.true.)
            call util_print_centered(LuOut,
     $           '(all occupied plus 10 virtual)',20,.true.)
            call output(dbl_mb(k_eval(2)),
     $           1, min(noc(1)+10,nbf_ao),
     $           1, 1, nbf_ao, 1, 1)
          endif
          call util_flush(6)
        endif
        if (util_print('final vectors', print_debug)) then
          write(LuOut,*)
          call util_print_centered(
     (         LuOut,'Final MO vectors',40,.true.)
          write(LuOut,*)
          call util_flush(LuOut)
        end if
      end if
      if (oadapt) then
        nprint = min(noc(1)+10,nbf_ao)
        if (util_print('all vector symmetries', print_high)) then
          nprint = nbf_ao
          oprint_sym = util_print('final vector symmetries',
     $         print_default)
          call scf_movecs_sym_adapt(ao_bas_han, g_movecs,
     $         oprint_sym,nprint, '- alpha', .false., 1)
          if(ipol.eq.2) 
     $         call scf_movecs_sym_adapt(ao_bas_han, g_movecs(2), 
     $         oprint_sym, nprint, '- beta', .false., 1)
        endif
      endif
      if (util_print('final vectors', print_debug)) then
        call ga_print(g_movecs)
        if(ipol.eq.2) call ga_print(g_movecs(2))
      end if
c****   
c****   eval deallocation moved here from inside iteration loop
c****   
        if(.not.ma_pop_stack(l_eval))
     &    call errquit('dft_scf: cannot pop stack',0)
        if(incore.and.CDFIT) then
          if(.not.ma_pop_stack(l_3cwhat))
     &      call errquit('dft_scf: cannot pop stack',0)
          if(.not.ma_pop_stack(l_3cERI))
     &      call errquit('dft_scf: cannot pop stack',0)
        endif
        if(.not.ma_pop_stack(l_occ))
     &    call errquit('dft_scf: cannot pop stack',0)
c     
c-----------------------------------------------------------------------

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
      write(*,*)' itrscf = ', itrscf
      write(*,*)' iscfcon = ', iscfcon
      write(*,*)' iguess = ', iguess
      write(*,*)' mult = ', mult
      write(*,*)' ipol = ', ipol
      write(*,*)' igcon = ', igcon
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
      write(*,*)' igeotype = ', igeotype
      write(*,*)' numgrp = ', numgrp
      write(*,*)' numset = ', numset
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
      write(*,*)' TOLL = ', TOLL
      write(*,*)' radius = ', radius
      write(*,*)' rscale = ', rscale
      write(*,*)' cdist = ', (cdist(i),i=1,3)
      write(*,*)' cang =', (cang(i),i=1,3)
      write(*,*)' cfac = ', (cfac(i),i=1,10)
      write(*,*)' xfac =', (xfac(i),i=1,10)
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
