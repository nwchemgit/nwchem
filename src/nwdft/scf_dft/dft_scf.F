      Subroutine dft_scf(rtdb,geom,AO_bas_han,CD_bas_han,XC_bas_han,
     &     ipol,Etold,Enuc,iVcoul_opt,iVxc_opt,
     &     iconst,nconst,iter,iExc,
     &     iga_dens,iga_evecs,iga_olddens,iga_vxc,
     &     iga_fock,iga_diis,iga_svecs,isvals,
     &     iga_cdinv,iga_xcinv,iga_2ceri)

C$Id: dft_scf.F,v 1.1 1995-04-27 03:16:20 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      Logical LResult,CDFIT,IOLGC,mulliken,do_gradients
      logical converged,XCFIT,wght_GA
      parameter(itol_min=5,iAOacc_min=15)
      
c**** 
c**** nwchem handles
c**** 
      integer rtdb, geom, AO_bas_han, CD_bas_han, XC_bas_han

#include <bas.fh>
#include <geom.fh>

#include <context.fh>

#include <mafdecls.h>
#include <stdio.fh>
#include <rtdb.fh>
#include <global.fh>
#include <tcgmsg.fh>
#include <msgids.fh>

#if defined(DFT_TIMER)
#include <timer.fh>
#endif
      Logical movecs_write, movecs_read_header, movecs_read
      External movecs_write,movecs_read_header, movecs_read

      Dimension iga_dens(2), iga_evecs(2), iga_vxc(4)
      Dimension iconst(nconst)
      Dimension Exc(2),noc(2),nmo(2),nmo_read(2),icall(2)

      Integer lSbn2
      Dimension iSbn2(2)
      Character*255 movecs_out
      character*255 title_vecs, basis_vecs
      Character*255 name, trans 
c-----------------------------------------------------------------------
c     
c     Determine the characteristics of the AO and CD Gaussian basis sets.
c**** (from NWdft data structures, not the nwchem objects)
c     
      if ( .not. bas_numbf(AO_bas_han,nbf_ao) ) then
        call errquit('Exiting from dft',1)
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
      itol2e=min(itol_min,itol_max)
      iAOacc=min(iAOacc_min,iAOacc_max)
      mulliken=.false.
      if(imull.eq.1) mulliken=.true.
      itrscf  = iconst(8)
      iscfcon = iconst(9)
      ichg    = iconst(12) 
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
c     
c     SCF energy convergence criterion. 
c     
      scfcon = 10.D0**(-iscfcon)
      zero = 0.0d0
      one = 1.0d0
      onem = -1.0d0
      call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'ga temp',iga_tmp,'atom')
      if(ipol.gt.1) then
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'fock tr',iga_fockt,'atom')
        if(icvgtyp.eq.2)
     &       call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'fock save',iga_jeff,'atom')
      else
        iga_fockt=iga_fock
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

      if (.not. rtdb_cget(rtdb, 'dft:output vectors',
     &     1, movecs_out)) then
        movecs_out = "movecs"
      endif 
c
C  get MO vectors from file
c
      LResult = movecs_read_header(movecs_out, title_vecs, basis_vecs,
     &     nbf_vecs, nsets, nmo_read, 2)
      if(Lresult.and.
     &     (ipol.eq.nsets).and.
     &     (nmo(1).eq.nmo_read(1)).and.
     &     (nbf_vecs.eq.nbf_ao))then
        if(.not. MA_Push_Get(
     &       MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ))
     &  call errquit('dft: failed to alloc',999)
        LResult=.true.
        do ispin = 1,ipol
c****     
c****     get orbitals (i.e. fock eigenvectors) and fock eigenvalues
c****     
          LResult = LResult .and.
     &         movecs_read(movecs_out, ispin, dbl_mb(k_occ),
     &         dbl_mb(iSbn2(ispin)), iga_evecs(ispin))
          call ga_orthog_mos(AO_bas_han, iga_evecs(ispin))
          if(LResult)
     &         call dft_densm(iga_dens(ispin),iga_evecs(ispin),
     &         nbf_ao,ipol,noc(ispin))
        enddo 
        if(.not.MA_Pop_Stack(l_occ))
     &       call errquit(' failed pop_stack in dft ',333)
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
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'ovl ',iga_s,'atom')
        call ga_zero(iga_s)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       iga_s,'overlap',.false.)
        
c       
c       Form the total density matrix.
c       
        if( ipol.eq.2 )then
          call ga_dadd(one,iga_dens(1),one,iga_dens(2),iga_dens(1))
        end if
c       
c       call to Mulliken Pop Ananlysis
c       
        if (mulliken) then
          if(me.eq.0) call dft_header(
     &         ' Total Density - Mulliken Population Analysis')
          call mull_pop(geom,AO_bas_han,iga_dens(1),iga_s)
          if(ipol.eq.2) then
            if(me.eq.0) call dft_header(
     &           ' Beta Density - Mulliken Population Analysis')
            call mull_pop(geom,AO_bas_han,iga_dens(2),iga_s)
          endif
        endif
        call gadest(iga_s)
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
        call ga_zero(iga_dens(1))
        call guess_dens(geom,AO_bas_han,iga_dens,ipol)
        if(ipol.eq.2) then
          call ga_copy(iga_dens(1),iga_dens(2))
          rb=dble(noc(2))/dble(nTotEl)
C          if(me.eq.0) write(6,*) ' rb ',rb
          call ga_dscal(iga_dens(2),rb)
        endif
        if(ichg.ne.0) then
          afact=dble(ntotEl)/dble(ntotEl+ichg)
C          if(me.eq.0) write(6,*) ' afact ',afact
          call ga_dscal(iga_dens(1),afact)
          if(ipol.eq.2) call ga_dscal(iga_dens(2),afact)
        endif
c       
c       check normalization of DM
c       
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'ovl ',iga_s,'atom')
        call ga_zero(iga_s)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       iga_s,'overlap',.false.)
        if(ichg.ne.0) then
          anPS = ga_ddot(iga_dens(1),iga_s)
          anel=dble(ntotEl)
          if(dabs(anPS-anel).gt.1.d-7) then
            if(me.eq.0) write(LuOut,*) ' renormalizing DM'
            facto=anel/anPS
            if(me.eq.0) write(LuOut,*) ' noc anel ',anPS,anel
            call ga_dscal(iga_dens(1),facto)
          endif
        endif
        if(ipol.eq.2) then
          anPS = ga_ddot(iga_dens(2),iga_s)
          anel=dble(noc(2))
          if(dabs(anPS-anel).gt.1.d-7) then
            if(me.eq.0) write(luOut,*) ' renormalizing DM'
            facto=anel/anPS
            if(me.eq.0) write(LuOut,*) ' noc anel ',anPS,anel
            call ga_dscal(iga_dens(2),facto)
          endif
        endif
        call gadest(iga_s)
        iter=1
c       
c       call to Mulliken Pop Ananlysis
c       
        if (mulliken) then
          call gacrea(geom,AO_bas_han,
     &         nbf_ao,nbf_ao,'ovl ',iga_s,'atom')
          call ga_zero(iga_s)
          call int_1e_ga(AO_bas_han,AO_bas_han,
     &         iga_s,'overlap',.false.)
          if(me.eq.0) call dft_header(
     &         ' Total Density - Mulliken Population Analysis')
          call mull_pop(geom,AO_bas_han,iga_dens(1),iga_s)
          if(ipol.eq.2) then
            if(me.eq.0) call dft_header(
     &           ' Beta Density - Mulliken Population Analysis')
            call mull_pop(geom,AO_bas_han,iga_dens(2),iga_s)
          endif
          call gadest(iga_s)
        endif
      endif
      if( iter.gt.0 )then
c       
c       Construct goofy scaled density for fitcd, getvcoul, etc.
c       
        call ga_copy(iga_dens(1),iga_tmp)
        do i=me+1,nbf_ao,nproc
          call ga_get(iga_tmp,
     *         i,i,i,i,aaz,1)
          aaz=aaz*0.5d0
          call ga_put(iga_tmp,
     *         i,i,i,i,aaz,1)
        enddo
        call ga_sync
        call ga_dscal(iga_tmp,2.D0)
C       
C       initialize fock matrix
C       
#ifdef DATA_SERVER
        call ga_release(iga_fock,1,nbf_ao,1,nbf_ao)
#endif  
        call ga_zero(iga_fock)
C       
C       accumulate core hamiltonian into Fock matrix
C       compute core energy
C       
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       iga_fock,'kinetic',.false.)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       iga_fock,'potential',.false.)
        Ecore = ga_ddot(iga_dens(1),iga_fock)

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
     &         iga_cdinv, iga_tmp, nTotEl, IOLGC,itol2e,
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
     &       Dbl_MB(icd_coef),iga_tmp,
     &       Ecoul,nbf_cd,iVcoul_opt,
     &       iga_2ceri,iga_dens,IOLGC,itol2e)
c       
C       Add V coul to Fock Matrix
C       
        call ga_dadd(one,iga_tmp,one,iga_fock,iga_fock)
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
          call ga_dadd(one,iga_dens(1),onem,iga_dens(2),iga_dens(1))
        endif
c       
c       Compute the XC potential and energy.
c       
        iga_vxc(1)=iga_tmp
        call ga_zero(iga_vxc(1))
        if(ipol.eq.2)
     &       call ga_zero(iga_vxc(2))
        call dft_tstamp(' Before call to GETVXC.  ')
        call xc_getv(rtdb,geom,AO_bas_han,XC_bas_han,
     &       Exc,ipol,nExc,
     &       iExc,iVxc_opt,iga_xcinv,
     &       iga_dens,iga_vxc,IOLGC,iAOacc,
     &       iga_wght,iga_xyz,wght_GA)

        call ga_dadd(1.d0,iga_vxc(1),1.d0,iga_fock,iga_fock)
        if(ipol.gt.1) then
c         
C         
c         Put (-alpha+beta) in iga_vxc(2)
c         
          call ga_dadd(-1.d0,iga_vxc(1),1.d0,iga_vxc(2),iga_vxc(2))
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
            call ga_zero(iga_fock)
C           
C           accumulate core hamiltonian into Fock matrix
C           compute core energy
C           
            call int_1e_ga(AO_bas_han,AO_bas_han,
     &           iga_fock,'kinetic',.false.)
            call int_1e_ga(AO_bas_han,AO_bas_han,
     &           iga_fock,'potential',.false.)
          endif
          call ga_sync
        elseif(iter.gt.0) then
C         
C         Add V xc to Fock Matrix
C         
          if(icvgtyp.eq.2.and.ispin.eq.2) 
     *         call ga_copy(iga_jeff,iga_fock)
          if(ispin.eq.2) then
#ifdef DATA_SERVER
            call ga_release(iga_fock,1,nbf_ao,1,nbf_ao)
#endif      
            call ga_dadd(one,iga_vxc(ispin),one,iga_fock,iga_fock)
          endif
c         
#if defined(DFT_TIMER)
          call dft_second(tp3)
#endif    
c         
          delta=abs(etold-etnew)
        if(delta.lt.1.d-5) then
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
     &           iga_fock,iga_dens(1),iga_tmp,
     &           iga_svecs,iga_diis,isvals)
            icall(1)=0
            icall(2)=0
          endif
          iswitc=iswitc+1
        endif

          if( icvgtyp.eq.2.and.delta.gt.10.d0**(-idisoff)) then
            if(ispin.eq.1.and.ipol.gt.1) 
     &           call ga_copy(iga_fock,iga_jeff)
            call diis_driver(ispin,ipol,icall,
     &           nfock,nbf_ao,geom,AO_bas_han,
     &           iga_fock,iga_dens(ispin),iga_tmp,
     &           iga_svecs,iga_diis,isvals)
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
     &       iga_evecs(ispin),iga_fock,zero,iga_tmp)
        call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &       iga_tmp,iga_evecs(ispin),zero,iga_fockt)
c        if(lshift.ne.0.and.iter.gt.1.and.iter.lt.ncysh) then
        if(lshift.ne.0.and.iter.lt.ncysh) then
          bb=0.1d0*lshift
          do j=noc(ispin)+1+me,nbf_ao,nproc
            call ga_get(iga_fockt,
     *           j,j,j,j,AAZ,1)
            AAZ=AAZ+bb
            call ga_put(iga_fockt,
     *           j,j,j,j,AAZ,1)
          enddo
        endif
        call ga_sync
#if defined(PARALLEL_DIAG)
        call ga_diag_std(iga_fockt,iga_tmp,Dbl_MB(ISbn2(ispin)))
#else   
        call ga_diag_std_seq(iga_fockt,iga_tmp,Dbl_MB(ISbn2(ispin)))
#endif  
        call ga_dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,one,
     &       iga_evecs(ispin),iga_tmp,zero,iga_fockt)
        call ga_copy(iga_fockt,iga_evecs(ispin))
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
          call ga_copy(iga_dens(1),iga_tmp)
        else
          call ga_copy(iga_dens(2),iga_olddens)
          call ga_dadd(one,iga_tmp,one,iga_dens(2),iga_tmp)
        endif
c       
c       Form a new density matrix.
c       
        call dft_densm(iga_dens(ispin),iga_evecs(ispin),
     &       nbf_ao,ipol,noc(ispin))
c       
        call ga_sync
        if(ndamp.ne.0.and.iter.gt.0.and.iter.lt.ncysh) then
          pp=ndamp*1.d-2
          onempp=1.0d0 - pp
          if(ispin.eq.1)then
            call ga_dadd(
     &           pp,iga_tmp,onempp,iga_dens(ispin),iga_dens(ispin))
          else
            call ga_dadd(
     &           pp,iga_olddens,onempp,iga_dens(ispin),iga_dens(ispin))
          endif
        endif
c
c           check convergence on Density
c
        if(ispin.eq.1) then
            rms=rms+dft_dencvg(iga_dens(ispin),iga_tmp,nbf_ao)
          else
            rms=rms+dft_dencvg(iga_dens(ispin),iga_olddens,nbf_ao)
          endif
        call ga_sync
      enddo     ! end loop over ispin
      call ga_sync
c
c     save eigenvectors to movecs file
c
      if(iter.gt.0) then
        
        LResult = MA_Push_Get(
     &       MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ)
        if ( ipol .eq. 1 ) then
c****     RHF occupations
          call dfill(nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
          do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 2.0d0
          enddo
        else  
c****     UHF occupations
          call dfill(2*nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
          do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 1.0d0
          enddo
          do i = nbf_ao+1, nbf_ao+noc(2)
            dbl_mb(i-1+k_occ) = 1.0d0
          enddo
        endif 
        if (.not.
     &       movecs_write(movecs_out, 'Kohn-Sham orbitals',
     &       'AO_basis',
     &       nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &       dbl_mb(ISbn2(1)), nbf_ao, iga_evecs))
     &       call errquit('rhf: movec_write failed', 0)

        LResult = MA_Pop_Stack(l_occ)

        call ga_sync

      endif
c     
c     Form the total density matrix.
c     
      if( ipol.eq.2 )then
        call ga_dadd(one,iga_dens(1),one,iga_dens(2),iga_dens(1))
      end if
      call ga_sync
c     
c     call to Mulliken Pop Ananlysis
c     
      if(mulliken) then
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'ovl ',iga_s,'atom')
        call ga_zero(iga_s)
        call int_1e_ga(AO_bas_han,AO_bas_han,
     &       iga_s,'overlap',.false.)
        if(me.eq.0)
     &       call dft_header(
     &       ' Total Density - Mulliken Population Analysis')
        call mull_pop(geom,AO_bas_han,iga_dens(1),iga_s)
        if(ipol.eq.2) then
c         
c         analysis of spin density
c         
          if(me.eq.0) call dft_header(
     &         ' Spin Density - Mulliken Population Analysis')
          call ga_dadd(one,iga_dens(1),-2.d0,iga_dens(2),iga_dens(2))
          call mull_pop(geom,AO_bas_han,iga_dens(2),iga_s) 
          call ga_dadd(one,iga_dens(1),-1.d0,iga_dens(2),iga_dens(2))
          call ga_dscal(iga_dens(2),0.5d0)
        endif
        call gadest(iga_s)
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
        Ecore = ga_ddot(iga_dens(1),iga_fock)
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
          LResult = MA_Push_Get(
     &         MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ)
c****     
c****     not sure what this occupation number array is good for,
c****     but filling it anyway...  NOT TESTED
c****     
          if ( ipol .eq. 1 ) then
c****       RHF occupations
            call dfill(nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
            do i = 1, noc(1)
              dbl_mb(i-1+k_occ) = 2.0d0
            enddo
          else  ! UHF
c****       UHF occupations
            call dfill(2*nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
            do i = 1, noc(1)
              dbl_mb(i-1+k_occ) = 1.0d0
            enddo
            do i = nbf_ao+1, nbf_ao+noc(2)
              dbl_mb(i-1+k_occ) = 1.0d0
            enddo
          endif 


c****     
c****     logical function movecs_write(filename, title, basis_name, 
c****     $     nbf_ao, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
c****     occ is an occupation vector (not used here)
c****     
          if (.not.
     &         movecs_write(movecs_out, 'Kohn-Sham orbitals',
     &         'Unknown',
     &         nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &         dbl_mb(ISbn2(1)), nbf_ao, iga_evecs))
     &         call errquit('rhf: movec_write failed', 0)
          call ga_sync()
          if(me.eq.0) then
            write (6,*) 'Alpha Eigenvalues: '
            write (6,'(5x,4(e14.4))')
     &           (dbl_mb(ISbn2(1)+ii),ii=0,nbf_ao-1)
            write (6,*) 
            if ( ipol.eq.2 ) then
              write (6,*) 'Beta Eigenvalues: '
              write (6,'(5x,4(e14.4))')
     &             (dbl_mb(ISbn2(2)+ii),ii=0,nbf_ao-1)
              write (6,*) 
          endif 
          endif
          LResult = MA_Pop_Stack(l_occ)

        endif 
c-----------------------------------------------------------------------
c****   
c****   eval deallocation moved here from inside iteration loop
c****   
        LResult = MA_Pop_Stack(lSbn2)

        if(icvgtyp.eq.2.and.icall(1).gt.0) then
c         
c         One call kills all
c         
          icall(1) = -1
          call diis_driver(
     &         1,ipol,icall,
     &         nfock,nbf_ao,geom,AO_bas_han,
     &         iga_fock,iga_dens(1),iga_tmp,
     &         iga_svecs,iga_diis,isvals)
        end if
      endif
      if(ipol.gt.1) then
        call gadest(iga_fockt)
        if(icvgtyp.eq.2)
     &       call gadest(iga_jeff)
      endif
      call gadest(iga_tmp)
      if(wght_GA) then
        call gadest(iga_wght)
        call gadest(iga_xyz)
      endif
c-----------------------------------------------------------------------
c**** 
c**** this is the final section of the ddscf routine rhf()
c**** have substituted SCF with DFT
c**** 
c     
c     Store energy and convergence status ... not sure if better to
c     do this here or have each wavefunction type do it
c     
      if (.not. bas_name(AO_bas_han, name, trans))
     $     call errquit('scf: bas_name?', 0)
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
c     Overwrite the input movecs option with where the final movecs
c     ended up so that a restart automatically picks them up
c     
      if (.not. rtdb_cput(rtdb, 'dft:input vectors', 1, movecs_out))
     $     call errquit('dft:rtdb_cput failed dft:input vectors',0)

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
