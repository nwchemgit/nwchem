      Subroutine dft_scf(rtdb,
     &     ipol,Etold,Enuc,iVcoul_opt,iVxc_opt,
     &     iconst,nconst,iter,
     &     g_dens,g_evecs,g_olddens,g_vxc,
     &     g_fock,g_diis,g_svecs,isvals,
     &     g_cdinv,g_xcinv,g_2ceri)

C$Id: dft_scf.F,v 1.21 1995-12-27 23:30:04 og845 Exp $
      implicit none
      integer rtdb ! [input]
      integer ipol
      double precision Etold,Enuc
      integer iVcoul_opt
      integer iVxc_opt
      integer iconst(*),nconst,iter
      integer g_dens(*),g_evecs(*),g_olddens,g_vxc(*),
     &     g_fock,g_diis,g_svecs,isvals,
     &     g_cdinv,g_xcinv,g_2ceri
      
c**** patches for SC95,  ECCE' parsing
***** #define SC95

#include "bas.fh"
#include "geom.fh"

#include "context.fh"

#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cfock.fh"
#include "cdft.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"

#if defined(DFT_TIMER)
#include "timer.fh"
#endif
      Logical movecs_write, movecs_read_header, movecs_read
      External movecs_write,movecs_read_header, movecs_read
      
      double precision Exc(2)
      integer noc(2),nmo(2),nmo_read(2),icall(2)

      Integer lSbn2,Type
      integer iSbn2(2)
      Character*255 movecs_out,movecs_in
      character*255 title_vecs, basis_vecs
      Character*255 name, trans 
      Double Precision rmoments(3)  ! dipole
      logical ltight ! if TRUE always tight computing tolerances
      integer lmax   ! CMM
      integer finest ! CMM
      character*255 filename
      logical int2e_file_open,int2e_test_mem,int2e_file_close

      integer natoms,nTotEl
      integer l_occ,k_occ
      integer n_full,n_empty
      integer i,ii,j
      integer n_ashift,n_bshift
      integer me,nproc
      integer g_jeff,g_tmp,g_fockt,g_s,g_wght,g_xyz
      integer nheap,nsets,nstack
      integer nbf_vecs,ispin,nxang,idone
      integer ncysh,lshift,ndamp,itrscf
      integer igcon,idisoff,isdiag,ielem
      integer nexc,nfock,itol2e
      integer l_afshift,n_afshift,i_afshift
      integer l_bfshift,n_bfshift,i_bfshift
      integer icvgtyp,iaoacc,iswitc
      integer itol_max,iaoacc_max
      integer imull,iguess,idmtot
      integer leneval,lcd_coef,icd_coef
      double precision ecoul,rb,facto,ecore,afact,anps,aaz
      double precision rms,pp,bb,delta
      double precision anucl_charg,rcharge,anel,anew_nel
      double precision anoca,anocb,onempp
      double precision etnew,tol2e,fock_el
      Logical LResult,CDFIT,IOLGC,mulliken,do_gradients,lSPM
      logical converged,XCFIT,wght_GA
      logical lfocka,lfockb,lmaxov
      integer itol_min,iAOacc_min,iscfcon
      double precision scfcon
      double precision zero,onem,one
      parameter(zero=0.d0,one=1.d0,onem=-one)
      parameter(itol_min=6,iAOacc_min=15)
      character*3 scftype
      double precision dft_dencvg,dabssum
      external dft_dencvg,dabssum
#if defined(SC95)
      double precision homo_lumo_gap
      logical movecs_to_logfile
      integer n_homo, n_lumo            ! # of HOMOs, LUMOs to output to stdout
#endif      
c-----------------------------------------------------------------------
c     
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c**** (from NWdft data structures, not the nwchem objects)
c     
      if ( .not. bas_numbf(AO_bas_han,nbf_ao) ) then
        call errquit('Exiting from dft',1)
      endif 
      if (.not. geom_ncent(geom, natoms))
     &     call errquit('dft_scf: geom_ncent failed',73)      
      if (.not. geom_nuc_charge(geom, anucl_charg))
     $     call errquit('scf: geom_nuc_charge failed', 0)
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, rcharge)) 
     %     rcharge = 0.0d0
      anel = int(anucl_charg) - rcharge
C
C     determine pattern of orbitals' occupancy
C
      if(.not. MA_Push_Get(
     &     MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ))
     &     call errquit('dft: failed to alloc',999)
      
      if( ipol.eq.1 )then
c       
c       noc(1) = No. of electrons
c       
        noc(1) = iconst(15)
        noc(2) = 0
        nTotEl = 2*noc(1)
        nmo(1)=nbf_ao
      else
c       
c       noc(1) = No. of spin alpha electrons
c       noc(2) = No. of spin beta  electrons
c       
        noc(1) = iconst(16)
        noc(2) = iconst(17)
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
C     
C     shift on Fock matrix to change occupancy
C
      n_ashift=0
      n_bshift=0
      lfocka = RTDB_MA_Get(RTDB, 'dft:alpha empty', Type,
     $   n_afshift, l_afshift)
      lfockb = RTDB_MA_Get(RTDB, 'dft:beta empty', Type,
     $   n_bfshift, l_bfshift)
      if(lfocka) then
        if(me.eq.0) write(LuOut,*) ' Fock shifta'
        If ( Type .ne. MT_Int)  Call ErrQuit(
     $       'dft_scf: Bad MA type from RTDB_MA_Get', Type)
        If ( .NOT. MA_Get_Index(l_afshift, i_afshift) ) Call ErrQuit(
     $       'dft_scf: Bad MA handle from RTDB_MA_Get', l_afshift)
      endif
      if(lfockb) then
        if(me.eq.0) write(LuOut,*) ' Fock shiftb'
        If ( Type .ne. MT_Int)  Call ErrQuit(
     $       'dft_scf: Bad MA type from RTDB_MA_Get', Type)
        If ( .NOT. MA_Get_Index(l_bfshift, i_bfshift) ) Call ErrQuit(
     $       'dft_scf: Bad MA handle from RTDB_MA_Get', l_bfshift)
      endif


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
      if (.not. rtdb_get(rtdb, 'dft:tight', mt_log, 1, ltight))
     $     ltight=.false.
      if (.not. rtdb_get(rtdb, 'dft:iguess', mt_int, 1, iguess))
     &     call errquit('dft: rtdb_get failed',1) 
      if (.not. rtdb_get(rtdb, 'dft:imull', mt_int, 1, imull))
     $     call errquit('dft: rtdb_get failed', 2)
      if (.not. rtdb_get(rtdb, 'dft:itol2e', mt_int, 1, itol_max))
     $     call errquit('fitcd: rtdb_get failed', 0)
      if (.not. rtdb_get(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc_max))
     $     call errquit('quadvxc0: rtdb_get failed', 0)
c
c     set initial coulomb acc to 1e-5
c
      iswitc=0
      if(ltight) then
        itol2e=itol_max
        iAOacc=iAOacc_max
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
      itrscf  = iconst(8)
      iscfcon = iconst(9)
      igcon   = iconst(18)
      icvgtyp = iconst(36)
      nfock   = iconst(37)
      idisoff = iconst(40)
      lshift=iconst(50)
      ncysh=iconst(52)
      ndamp=iconst(51)
      IOLGC=.TRUE.
      if(iconst(49).eq.1) IOLGC=.FALSE.
      isdiag   = iconst(46)
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
c     Energy decomposition switch
c     
      nExc    = iconst(39) + 1
      Etnew=0.d0

c     
c     SCF energy convergence criterion. 
c     
      scfcon = 10.D0**(-iscfcon)
      call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'ga temp',g_tmp,'atom')
      if(ipol.gt.1) then
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'fock tr',g_fockt,'atom')
        if(icvgtyp.eq.2)
     &       call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'fock save',g_jeff,'atom')
      else
        g_fockt=g_fock
      endif
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
        write(LuOut,'(a,f10.2,i20)')
     &       ' Heap Space remaining (MW):  ',dble(nheap)*1.D-06,nheap
        write(LuOut,'(a,f10.2,i20)')
     &       'Stack Space remaining (MW):  ',dble(nstack)*1.D-06,nstack
	call flush(LuOut)
      end if
c     
c     Construct a temporary global array.
c     
c     
c**** 
c**** top of infinite SCF iteration loop
c**** 
c**** 
c**** allocate eigenvalue array, including second pointer to beta block
c**** 
      leneval = nbf_ao * ipol
      LResult = MA_Push_Get(MT_Dbl,leneval,'eval',lSbn2,iSbn2(1))
      if ( ipol .eq. 2) then  
        iSbn2(2) = iSbn2(1) + nbf_ao
      endif 

c
c     Figure input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c     
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $     movecs_in = 'movecs'
      if (.not. rtdb_cget(rtdb, 'dft:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
            movecs_out = movecs_in
      endif
c
c  get MO vectors from file
c
      LResult = movecs_read_header(movecs_in, title_vecs, basis_vecs,
     &     nbf_vecs, nsets, nmo_read, 2)
      if(Lresult.and.
     &     (nbf_vecs.eq.nbf_ao).and.
     &     (nmo(1).eq.nmo_read(1))) then
        LResult=.true.
        if(ipol.eq.1.and.nsets.eq.2) nsets=1
        do ispin = 1,nsets
c****     
c****     get orbitals (i.e. fock eigenvectors) and fock eigenvalues
c****     
          LResult = LResult .and.
     &         movecs_read(movecs_in, ispin, dbl_mb(k_occ),
     &         dbl_mb(iSbn2(ispin)), g_evecs(ispin))
          call ga_orthog_mos(AO_bas_han, g_evecs(ispin))
        enddo 
          if(LResult) then
            if(ipol.eq.2.and.nsets.eq.1) then
              call ga_copy(g_evecs(1),g_evecs(2))
            endif
c
c     Apply any swaps of orbitals
c
            scftype='RHF'
            if(ipol.eq.2) scftype='UHF'
            call movecs_swap(rtdb, scftype, g_evecs,
     $           dbl_mb(k_occ), dbl_mb(Isbn2(1)))
            do ispin=1,ipol
              call dft_densm(g_dens(ispin),g_evecs(ispin),
     &             nbf_ao,ipol,noc(ispin))
            enddo
          endif
      else
        LResult=.false.
      endif
      if(LResult) then
        iter=1
        iguess=999
        itol2e=itol_max
        iAOacc=iAOacc_max
        iswitc=1
        if(me.eq.0)  
     &     write(LuOut,*) ' read movecs file as starting guess '
        
c       
c       Form the total density matrix.
c       
        if( ipol.eq.2 )then
          call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
        end if
c       
c       call to Mulliken Pop Ananlysis
c       
        if (mulliken) then
          call gacrea(geom,AO_bas_han,
     &         nbf_ao,nbf_ao,'ovl ',g_s,'atom')
          call ga_zero(g_s)
          call int_1e_ga(AO_bas_han,AO_bas_han,
     &         g_s,'overlap',.false.)
          if(me.eq.0) call dft_header(
     &         ' Total Density - Mulliken Population Analysis')
          call mull_pop(geom,AO_bas_han,g_dens(1),g_s)
          if(ipol.eq.2) then
            if(me.eq.0) call dft_header(
     &           ' Beta Density - Mulliken Population Analysis')
            call mull_pop(geom,AO_bas_han,g_dens(2),g_s)
          endif
          call gadest(g_s)
        endif
      endif

 1000 continue

c     
c     Compute the contribution to the total energy from the core Hamiltonian.
c     

      if( me.eq.0 ) write(LuOut,9000) iter
      Ecoul  = ZERO
      Exc(1) = ZERO
      Exc(2) = ZERO
      if(iguess.eq.0.and.iter.eq.0) then
C
C       call to atomscf for initial density matrix
C
        call ga_zero(g_dens(1))
        call guess_dens(geom,AO_bas_han,g_dens)
        if(ipol.eq.2) then
          call ga_copy(g_dens(1),g_dens(2))
          if((lfocka).or.(lfockb)) then
            call ga_dscal(g_dens(2),.5d0)
          else
            rb=anocb/anel
            if(me.eq.0) write(luerr,*) ' scaling beta DM by',rb
            call ga_dscal(g_dens(2),rb)
          endif
        endif

        if(abs(rcharge).gt.0.001) then
          afact=anel/dble(anel+rcharge)
          call ga_dscal(g_dens(1),afact)
          if(ipol.eq.2) call ga_dscal(g_dens(2),afact)
        endif
c       
c       check normalization of DM
c       
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'ovl ',g_s,'atom')
        call ga_zero(g_s)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       g_s,'overlap',.false.)
        if(abs(rcharge).gt.0.001) then
          anPS = ga_ddot(g_dens(1),g_s)
          if(dabs(anPS-anel).gt.1.d-7) then
            if(me.eq.0) write(LuOut,*) ' renormalizing tot DM '
            facto=anel/anPS
            if(me.eq.0) write(LuOut,*) ' noc anel ',anPS,anel
            call ga_dscal(g_dens(1),facto)
          endif
        endif
        if(ipol.eq.2) then
          anPS = ga_ddot(g_dens(2),g_s)
          if(dabs(anPS-anocb).gt.1.d-7) then
            if(me.eq.0) write(luOut,*) ' renormalizing beta DM'
            facto=anocb/anPS
            if(me.eq.0) write(LuOut,*) ' noc anel ',anPS,anocb
            call ga_dscal(g_dens(2),facto)
          endif
        endif
        call gadest(g_s)
        iter=1
c       
c       call to Mulliken Pop Ananlysis
c       
        if (mulliken) then
          call gacrea(geom,AO_bas_han,
     &         nbf_ao,nbf_ao,'ovl ',g_s,'atom')
          call ga_zero(g_s)
          call int_1e_ga(AO_bas_han,AO_bas_han,
     &         g_s,'overlap',.false.)
          if(me.eq.0) call dft_header(
     &         ' Total Density - Mulliken Population Analysis')
          call mull_pop(geom,AO_bas_han,g_dens(1),g_s)
          if(ipol.eq.2) then
            if(me.eq.0) call dft_header(
     &           ' Beta Density - Mulliken Population Analysis')
            call mull_pop(geom,AO_bas_han,g_dens(2),g_s)
          endif
          call gadest(g_s)
        endif
      endif
      if( iter.gt.0 )then
        call ga_zero(g_fock)
C       
C       accumulate core hamiltonian into Fock matrix
C       compute core energy
C       
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       g_fock,'kinetic',.false.)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       g_fock,'potential',.false.)
        Ecore = ga_ddot(g_dens(1),g_fock)


c
c       perform SPM analysis
c
        LSPM=.true.
      if (.not. rtdb_get(rtdb, 'dft:cmm', mt_int, 1, lmax))
     &     LSPM=.false.
        


      if(LSPM) then
        lresult=geom_charge_center(geom)


      if (.not. rtdb_get(rtdb, 'dft:finest', mt_int, 1, finest))
     &     finest=4
          call dft_tstamp(' Before call to CMM.   ')
        call dft_getspm(geom,lmax,finest,g_dens,AO_bas_han)
        call errquit(' end of CMM ',1234)

      endif


#if defined(DFT_TIMER)
        call dft_second(tp3)
#endif  

        if( CDFIT ) then
c         
c         Fit the electron charge density.
c         
          call dft_tstamp(' Before call to FITCD.   ')
          LResult = MA_Push_Get(MT_Dbl,nbf_cd,'cd_coef',
     &         lcd_coef,icd_coef)
          
          call dft_fitcd( rtdb,AO_bas_han, CD_bas_han,geom,
     &         Dbl_MB(icd_coef),  nbf_cd,
     &         g_cdinv, g_dens, nTotEl, IOLGC,itol2e,
     &         Ecoul)
        endif
c       
#if defined(DFT_TIMER)
        call dft_second(tp4)
        pgmsect(6)='   FitCD  '
        call tcllct(t_pgm(1,6),tp4,tp3)
        call dft_second(tp3)
#endif  
        call dft_tstamp(' Before call to GETVCOUL.')
        call dft_getvc(rtdb,geom,AO_bas_han,CD_bas_han,
     &       Dbl_MB(icd_coef),g_tmp,
     &       Ecoul,nbf_cd,iVcoul_opt,
     &       g_2ceri,g_dens,IOLGC,itol2e)
c       
C       Add V coul to Fock Matrix
C       
        call ga_dadd(one,g_tmp,one,g_fock,g_fock)
        if (CDFIT)
     &       LResult = MA_Pop_Stack(lcd_coef)
c       
#if defined(DFT_TIMER)
        call dft_second(tp4)
        pgmsect(7)='   GtVcoul'
        call tcllct(t_pgm(1,7),tp4,tp3)
        call dft_second(tp3)
#endif  
        if(ipol.gt.1) then
c         
C         Restore alpha and beta densities.
c         
          call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))
        endif
c       
c       Compute the XC potential and energy.
c       
        g_vxc(1)=g_tmp
        call ga_zero(g_vxc(1))
        if(ipol.eq.2)
     &       call ga_zero(g_vxc(2))
        call dft_tstamp(' Before call to GETVXC.  ')
        call xc_getv(rtdb,geom,AO_bas_han,XC_bas_han,
     &       Exc,ipol,nExc,
     &       iVxc_opt,g_xcinv,
     &       g_dens,g_vxc,IOLGC,iAOacc,
     &       g_wght,g_xyz,wght_GA)

        call ga_dadd(1.d0,g_vxc(1),1.d0,g_fock,g_fock)
        if(ipol.gt.1) then
c         
C         
c         Put (-alpha+beta) in g_vxc(2)
c         
          call ga_dadd(-1.d0,g_vxc(1),1.d0,g_vxc(2),g_vxc(2))
        endif
c       
        call dft_tstamp(' End of parallel region. ')
c       
#if defined(DFT_TIMER)
        call dft_second(tp4)
        pgmsect(8)='   GetVXC '
        call tcllct(t_pgm(1,8),tp4,tp3)
#endif  
c       
c       Calculate the total electronic energy.
c       
        if( nExc.eq.1 )then
          Etnew = Ecore + Ecoul + Exc(1)
        else
          Etnew = Ecore + Ecoul + Exc(1) + Exc(2)
        end if
        if( me.eq.0 )then
          if( nExc.eq.1 )then
            write(LuOut,1111) Ecore
            write(LuOut,1112) Ecoul
            write(LuOut,1113) Exc(1)
          else
            write(LuOut,1111) Ecore
            write(LuOut,1112) Ecoul
            write(LuOut,1114) Exc(1)
            write(LuOut,1115) Exc(2)
          end if
          call flush(LuOut)
        end if
        call ga_sync

      end if
c     
c     Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c     
      rms =0.d0
      do ispin = 1,ipol
        if(iter.eq.0) then
          if(iguess.eq.2) then
            call errquit(' guess read option no more available',999)
          elseif(iguess.eq.1) then
C           
C           initialize fock matrix
C           
            call ga_zero(g_fock)
C           
C           accumulate core hamiltonian into Fock matrix
C           compute core energy
C           
            call int_1e_ga(AO_bas_han,AO_bas_han,
     &           g_fock,'kinetic',.false.)
            call int_1e_ga(AO_bas_han,AO_bas_han,
     &           g_fock,'potential',.false.)
          endif
          call ga_sync
        elseif(iter.gt.0) then
C         
C         Add V xc to Fock Matrix
C         
          if(icvgtyp.eq.2.and.ispin.eq.2) 
     *         call ga_copy(g_jeff,g_fock)
          if(ispin.eq.2) then
            call ga_dadd(one,g_vxc(ispin),one,g_fock,g_fock)
          endif
c         
#if defined(DFT_TIMER)
          call dft_second(tp3)
#endif    
c         
          delta=abs(etold-etnew)
C
C         set tighter conditions when de < 1.e-4
C
        if(delta.lt.1.d-4) then
c
c         set  coulomb acc to input parameter
c
          itol2e=itol_max
          iAOacc=iAOacc_max
          if(iswitc.eq.0.and.icvgtyp.eq.2.and.icall(1).gt.0)then
            icall(1) = -1
            icall(2) = -1
            call diis_driver(
     &           1,ipol,icall,
     &           nfock,nbf_ao,geom,AO_bas_han,
     &           g_fock,g_dens(1),g_tmp,
     &           g_svecs,g_diis,isvals)
            icall(1)=0
            icall(2)=0
          endif
          iswitc=iswitc+1
        endif
C
C       rape the Fock matrix
C
        if(lfocka.and.iter.lt.2.and.ispin.eq.1) then
          if(me.eq.0) then
C           
C           alpha
C           
            do ii=0,n_afshift-1
              ielem= int_mb(i_afshift+ii)
              call ga_get(g_fock, ielem, ielem,
     *             ielem, ielem, fock_el, 1)
              fock_el = fock_el +5.d0
              write(luerr,*) ' alpha ',ii,fock_el
              call ga_put(g_fock, ielem, ielem,
     *             ielem, ielem, fock_el, 1)
            enddo
          endif
        endif
C       
C       beta
C       
        if(lfockb.and.ispin.eq.2.and.iter.lt.2) then
          if(me.eq.0) then
            do ii=0,n_bfshift-1
              ielem= int_mb(i_bfshift+ii)
              call ga_get(g_fock, ielem, ielem,
     *             ielem, ielem, fock_el, 1)
              fock_el = fock_el +5.d0
              call ga_put(g_fock, ielem, ielem,
     *             ielem, ielem, fock_el, 1)
              write(luerr,*) ' beta ',ii,fock_el
            enddo
          endif
          call ga_sync
        endif

          if( icvgtyp.eq.2.and.delta.gt.10.d0**(-idisoff)) then
            if(ispin.eq.1.and.ipol.gt.1) 
     &           call ga_copy(g_fock,g_jeff)
            call diis_driver(ispin,ipol,icall,
     &           nfock,nbf_ao,geom,AO_bas_han,
     &           g_fock,g_dens(ispin),g_tmp,
     &           g_svecs,g_diis,isvals)
          endif
c         
#if defined(DFT_TIMER)
          call dft_second(tp4)
          pgmsect(9)='   DIIS   '
          call tcllct(t_pgm(1,9),tp4,tp3)
#endif    
c         
        end if 
C       
C       level shifter: similarity transformation+ standard eigensolver
C       

#if defined(DFT_TIMER)
        call dft_second(tp3)
#endif  
        call ga_dgemm('T','N',nbf_ao,nbf_ao,nbf_ao,one,
     &       g_evecs(ispin),g_fock,zero,g_tmp)
        call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &       g_tmp,g_evecs(ispin),zero,g_fockt)

        if(lshift.ne.0.and.iter.lt.ncysh) then
          bb=0.1d0*lshift
          do j=noc(ispin)+1+me,nbf_ao,nproc
            call ga_get(g_fockt,
     *           j,j,j,j,AAZ,1)
            AAZ=AAZ+bb
            call ga_put(g_fockt,
     *           j,j,j,j,AAZ,1)
          enddo
        endif
        call ga_sync
#if defined(PARALLEL_DIAG)
        call ga_diag_std(g_fockt,g_tmp,Dbl_MB(ISbn2(ispin)))
#else   
        call ga_diag_std_seq(g_fockt,g_tmp,Dbl_MB(ISbn2(ispin)))
#endif  
        call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &       g_evecs(ispin),g_tmp,zero,g_fockt)
C
C       keep orbital ordering according to principle
C       of maximum overlap with previous iteration
C
        if (rtdb_get(rtdb, 'dft:max_ovl', mt_log, 1,lmaxov )) then
          
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'ovl ',g_s,'atom')
        call ga_zero(g_s)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       g_s,'overlap',.false.)
        call dft_mxovl(AO_bas_han, dbl_mb(Isbn2(ispin)),
     ,       g_fockt, g_evecs(ispin),g_s,g_tmp)
        call gadest(g_s)
        endif
        call ga_copy(g_fockt,g_evecs(ispin))
        call ga_sync
#if defined(DFT_TIMER)
        call dft_second(tp4)
        pgmsect(10)='   Diag   '
        call tcllct(t_pgm(1,10),tp4,tp3)
#endif  
c       
c       Save previous density for convergence check.
c       
c       Save the old alpha (ispin = 1) density matrix.
c       OR
c       Sum in the old beta (ispin = 2) density matrix.
c       
        if(ispin.eq.1)then
          call ga_copy(g_dens(1),g_tmp)
        else
          call ga_copy(g_dens(2),g_olddens)
          call ga_dadd(one,g_tmp,one,g_dens(2),g_tmp)
        endif
c       
c       Form a new density matrix.
c       
        call dft_densm(g_dens(ispin),g_evecs(ispin),
     &       nbf_ao,ipol,noc(ispin))
c       
        call ga_sync
        if(ndamp.ne.0.and.iter.gt.0.and.iter.lt.ncysh) then
          pp=ndamp*1.d-2
          onempp=1.0d0 - pp
          if(ispin.eq.1)then
            call ga_dadd(
     &           pp,g_tmp,onempp,g_dens(ispin),g_dens(ispin))
          else
            call ga_dadd(
     &           pp,g_olddens,onempp,g_dens(ispin),g_dens(ispin))
          endif
        endif
c
c           check convergence on Density
c
        if(ispin.eq.1) then
            rms=rms+dft_dencvg(g_dens(ispin),g_tmp,nbf_ao)
          else
            rms=rms+dft_dencvg(g_dens(ispin),g_olddens,nbf_ao)
          endif
        call ga_sync
      enddo     ! end loop over ispin
      call ga_sync
c
c     save eigenvectors to movecs file
c
      if(iter.gt.0) then
        
        if (.not.
     &       movecs_write(movecs_out, 'Kohn-Sham orbitals',
     &       'ao basis',
     &       nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &       dbl_mb(ISbn2(1)), nbf_ao, g_evecs))
     &       call errquit('dft_scf: movec_write failed', 0)

        call ga_sync

      endif
c     
c     Form the total density matrix.
c     
      if( ipol.eq.2 )then
        call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
          
        
      end if
      call ga_sync
c     
c     call to Mulliken Pop Ananlysis
c     
      if(mulliken) then
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'ovl ',g_s,'atom')
        call ga_zero(g_s)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       g_s,'overlap',.false.)
        if(me.eq.0)
     &       call dft_header(
     &       ' Total Density - Mulliken Population Analysis')
        call mull_pop(geom,AO_bas_han,g_dens(1),g_s)
        if(ipol.eq.2) then
c         
c         analysis of spin density
c         
          if(me.eq.0) call dft_header(
     &         ' Spin Density - Mulliken Population Analysis')
          call ga_dadd(one,g_dens(1),-2.d0,g_dens(2),g_dens(2))
          call mull_pop(geom,AO_bas_han,g_dens(2),g_s) 
          call ga_dadd(one,g_dens(1),-1.d0,g_dens(2),g_dens(2))
          call ga_dscal(g_dens(2),0.5d0)
        endif
        call gadest(g_s)
      endif
      
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
        write(LuOut,1111) Ecore
      else
        call ga_sync
        rms = rms/(ipol)
        call dft_scfcvg(rms,
     &       Etold,Etnew,Enuc,
     &       scfcon,igcon,iter,itrscf,idone,converged,rtdb)

        iter = iter + 1
        Etold = Etnew
        if( idone.eq.0.or.
     .   (iswitc.lt.2.and.iter.lt.itrscf)) go to 1000 ! begin new iteration
C
C       call for symmetry analysis
C
         call scf_movecs_sym_adapt(AO_bas_han, g_evecs(1),.true.,
     $        nbf_ao, '- alpha')
         if(ipol.eq.2)
     $        call scf_movecs_sym_adapt(AO_bas_han, g_evecs(2),
     $        .true., nbf_ao, '- beta')
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
c****     
c****     collect the eigenvectors from global storage; output
c****     to disk, allowing for RHF and UHF cases
c****     


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
     &         dbl_mb(ISbn2(1)), nbf_ao, g_evecs))
     &         call errquit('rhf: movec_write failed', 0)
          call ga_sync()
          if(me.eq.0) then
#if defined(SC95)
            homo_lumo_gap = min( dbl_mb( ISbn2(1)+noc(1)),
     &                           dbl_mb( ISbn2(2)+noc(2) ) )
     &           - max( dbl_mb( ISbn2(1)+noc(1)-1),
     &           dbl_mb( ISbn2(2)+noc(2)-1 ) )
            write (LuOut,*) 
            write (LuOut,'(a,e14.4)')
     &           'HOMO/LUMO Gap (Hartrees): ', homo_lumo_gap
            write (LuOut,*) 
            write (LuOut,*) 'Alpha-Occupied Eigenvalues: '
            write (LuOut,'(5x,4(e14.4))')
     &           (dbl_mb(ISbn2(1)+ii),ii=0,noc(1)-1)
            write (LuOut,*) 
            write (LuOut,*) 'Alpha-Virtual Eigenvalues: '
            write (LuOut,'(5x,4(e14.4))')
     &           (dbl_mb(ISbn2(1)+ii),ii=noc(1),nbf_ao-1)
            write (LuOut,*) 
            if ( ipol.eq.2 ) then
              write (LuOut,*) 'Beta-Occupied Eigenvalues: '
              write (LuOut,'(5x,4(e14.4))')
     &             (dbl_mb(ISbn2(2)+ii),ii=0,noc(2)-1)
              write (LuOut,*) 
              write (LuOut,*) 'Beta-Virtual Eigenvalues: '
              write (LuOut,'(5x,4(e14.4))')
     &             (dbl_mb(ISbn2(2)+ii),ii=noc(2),nbf_ao-1)
              write (LuOut,*) 
            endif 
#else
            write (LuOut,*) 'Alpha Eigenvalues: '
            write (LuOut,'(5x,4(e14.4))')
     &           (dbl_mb(ISbn2(1)+ii),ii=0,nbf_ao-1)
            write (LuOut,*) 
            if ( ipol.eq.2 ) then
              write (LuOut,*) 'Beta Eigenvalues: '
              write (LuOut,'(5x,4(e14.4))')
     &             (dbl_mb(ISbn2(2)+ii),ii=0,nbf_ao-1)
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
     &           dbl_mb(ISbn2(1)), nbf_ao, g_evecs))
     &           call errquit('dft_scf: movecs_to_logfile failed', 2)
          endif 
        endif 
#endif

c****   
c****   eval deallocation moved here from inside iteration loop
c****   
        LResult = MA_Pop_Stack(lSbn2)
        LResult = MA_Pop_Stack(l_occ)

        if(icvgtyp.eq.2.and.icall(1).gt.0) then
c         
c         One call kills all
c         
          icall(1) = -1
          call diis_driver(
     &         1,ipol,icall,
     &         nfock,nbf_ao,geom,AO_bas_han,
     &         g_fock,g_dens(1),g_tmp,
     &         g_svecs,g_diis,isvals)
        end if
      endif
      if(ipol.gt.1) then
        call gadest(g_fockt)
        if(icvgtyp.eq.2)
     &       call gadest(g_jeff)
      endif
      call gadest(g_tmp)
      if(wght_GA) then
        call gadest(g_wght)
        call gadest(g_xyz)
      endif
      if(lfocka) then
        If ( .NOT. MA_Free_Heap(l_afshift) ) Call ErrQuit(
     $       'dft_scf: Unable to free MA handle', l_bfshift)
        endif
        if(lfockb) then
        If ( .NOT. MA_Free_Heap(l_bfshift) ) Call ErrQuit(
     $       'dft_scf: Unable to free MA handle', l_afshift)
        endif
c
c
c
      if (filesize.gt.0 .or. memsize.gt.0) then
         if (.not. int2e_file_close())
     $        call errquit('scf_tidy: closing aoints?', 0)
      endif
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

c
c
c
c**** 
c**** this is the final section of the ddscf routine rhf()
c**** have substituted SCF with DFT
c**** 
c     
c     Store energy and convergence status ... not sure if better to
c     do this here or have each wavefunction type do it
c     
      if (.not. bas_name(AO_bas_han, name, trans))
     $     call errquit('dft_scf: bas_name?', 0)
      if (name .eq. ' ') name = 'ao basis'

      if (.not. context_prefix('energy', name)) call errquit
     &     ('dft: context_prefix failed', 0)
      if (.not. rtdb_put(rtdb, name, MT_DBL, 1, (Etold+Enuc))) 
     *     call errquit('dft: failed to store energy in rtdb', 0)
      if (.not. context_prefix('converged', name)) call errquit
     $     ('dft: context_prefix failed', 0)
      if (.not. rtdb_put(rtdb, name, MT_LOG, 1, converged)) 
     *     call errquit('dft: failed to store converged in rtdb', 0)
c     
c-----------------------------------------------------------------------

      return
c     
 1111 format(15x,'Core Energy:              ',f20.10)
 1112 format(15x,'Coulomb Energy:           ',f20.10)
 1113 format(15x,'Exch-Corr Energy:         ',f20.10)
 1114 format(15x,'Exchange Energy:          ',f20.10)
 1115 format(15x,'Correlation Energy:       ',f20.10)
 9000 format(/,5x,'Iteration number',I4)
c     
      end
