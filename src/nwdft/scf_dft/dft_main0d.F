      Subroutine dft_main0d( rtdb, geom, AO_bas_han, CD_bas_han,
     &                      XC_bas_han, iconst, nconst )

C$Id: dft_main0d.F,v 1.3 1995-10-17 08:54:53 d3g681 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical LResult,CDFIT,XCFIT,IOLGC

c**** 
c**** nwchem handles
c**** 
      integer rtdb, geom, AO_bas_han, CD_bas_han, XC_bas_han
c**** 

#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"

#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "schwarz.fh"
#include "msgids.fh"

#if defined(DFT_TIMER)
#include "timer.fh"
#endif

      Dimension iconst(nconst),iga_dens(2),iga_evecs(2),iga_vxc(4)

c*****************************************************************************

#if defined(DFT_TIMER)
      call t_stat(10,-2)
      call dft_second(tp1)
#endif

      me = ga_nodeid()
      nproc = ga_nnodes()
      call schwarz_init(geom,AO_bas_han)
      IOLGC=.TRUE.
      if(iconst(49).eq.1) IOLGC=.FALSE.

c     
c     Evaluate the nuclear replusion energy.
c     
      LResult = geom_nuc_rep_energy( geom, Enuc )
c     
c     Determine the characteristics of the AO Gaussian basis set.
c     This uses the old NWdft data structures, rather than tthe rtdb
c     
      if ( .not. bas_numbf(AO_bas_han,nbf_ao) ) then
        call errquit('Exiting from scfdriver',1)
      endif 

c     
c     Generate initial guess for the density matrix.
c     
      iter = 0
      Etold = 0.D0
c     
      itrscf  = iconst(8)
      iscfcon = iconst(9)
      ipol    = iconst(11)
      igcon   = iconst(18)
      icvgtyp = iconst(36)
      nfock   = iconst(37)
      nExc    = iconst(39) + 1
      isdiag  = iconst(46)
      lshift  = iconst(50)
      ncysh   = iconst(52)
      ndamp   = iconst(51)

c     
c     Select method to compute the XC potential.
c     
      if( XC_bas_han.eq.-99999 )then
        iVxc_opt = 0
        XCFIT    = .FALSE.
      else
        iVxc_opt = 1
        XCFIT    = .TRUE.
      end if

      if( XCFIT )then
c     
c       Determine the characteristics of the XC Gaussian basis set.
c     
        if ( .not. bas_numbf(XC_bas_han,nbf_xc) ) then
          call errquit('Exiting from scfdriver',4)
        endif 
        lmiss=1
        if(IOLGC)then
          call gacrea(
     &         geom,XC_bas_han,nbf_xc,nbf_xc,'XCinv',iga_xcinv,'row')
          if(me.eq.0) then
            call dft_invio('XC',iga_xcinv,nbf_xc,'READ',lmiss)
          endif
        call ga_sync
        call ga_brdcst(Msg_XC_Ex,lmiss,mitob(1),0)
        endif
        if(lmiss.eq.1) then
c
c Compute the matrix inverse of the 2-ctr OIs.
c         
          if(.not.IOLGC)             
     &         call gacrea(
     &         geom,XC_bas_han,nbf_xc,nbf_xc,'XCinv',iga_xcinv,'row')

          call gacrea(
     &         geom,XC_bas_han,nbf_xc,nbf_xc,'iga_oi',iga_oi,'row')
          call ga_zero(iga_oi)
          call xc_get2oi(XC_bas_han,iga_oi)
#if defined(PARALLEL_DIAG)         
          call ga_chol(iga_oi,iga_xcinv,info)
#else
          call ga_chol_seq(iga_oi,iga_xcinv,info)
#endif
	  if(info.ne.0.and.me.eq.0) 
     &	  write(LuErr,*) ' XC info = ',info
          if(info.eq.0) then
#if defined(PARALLEL_DIAG)         
            call ga_inverse(iga_xcinv,iga_oi)
            call ga_dgemm('T','N',nbf_xc,nbf_xc,nbf_xc,1.d0,
     &      iga_oi,iga_oi,0.d0,iga_xcinv)
#else
            call ga_copy(iga_xcinv,iga_oi)
            call ga_inv_seq(iga_oi,iga_xcinv)
#endif
          else
            call dft_invdiag(geom,XC_bas_han,iga_oi,iga_xcinv,nbf_xc)
          endif
          call ga_sync
          call gadest(iga_oi)
          if (me.eq.0.and.IOLGC) then
            call dft_invio('XC',iga_xcinv,nbf_xc,'WRITE',lmiss)
          end if
        endif
        call ga_sync
        if(IOLGC)  call gadest(iga_xcinv)
      end if

c     
#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(1)='  XCinv   '
      call tcllct(t_pgm(1,1),tp2,tp1)
      call dft_second(tp1)
#endif
c     
c     Select method to compute Coulomb potential.
c     
      if( CD_bas_han.eq.-99999 )then
        iVcoul_opt = 0
        CDFIT      = .FALSE.
      else
        iVcoul_opt = 1
        CDFIT      = .TRUE.
      end if

      if( CDFIT )then
c     
c       Determine the characteristics of the CD Gaussian basis set.
c       
        if ( .not. bas_numbf(CD_bas_han,nbf_cd) ) then
          call errquit('Exiting from scfdriver',5)
        endif 
c       
c       Compute the matrix inverse of the CD 2-ctr ERIs.
c       
        call gacrea(geom,CD_bas_han,
     &       nbf_cd,nbf_cd,'CD 2cERI',iga_2ceri,'row')
        call ga_zero(iga_2ceri)
        call dft_get2eri(CD_bas_han,iga_2ceri,geom)
        call gacrea(geom,CD_bas_han,
     &       nbf_cd,nbf_cd,'CD 2cERInv',iga_cdinv,'row')
        lmiss=1
        if( me.eq.0.and.IOLGC)then
          call dft_invio('CDI',iga_cdinv,nbf_cd,'READ',lmiss)
        endif
        call ga_brdcst(Msg_CD_Ex,lmiss,mitob(1),0)
        if(lmiss.eq.1) then
#if defined(PARALLEL_DIAG)         
        call ga_chol(iga_2ceri,iga_cdinv,info)
#else
        call ga_chol_seq(iga_2ceri,iga_cdinv,info)
#endif

	  if(info.ne.0.and.me.eq.0) 
     &         write(LuErr,*) ' CD info = ',info
          if(info.eq.0) then
            iga_tmpb=iga_2ceri
#if defined(PARALLEL_DIAG)         
            call ga_inverse(iga_cdinv,iga_tmpb)
            call ga_dgemm('T','N',nbf_cd,nbf_cd,nbf_cd,1.d0,
     &           iga_tmpb,iga_tmpb,0.d0,iga_cdinv)
#else
            call ga_copy(iga_cdinv,iga_tmpb)
            call ga_inv_seq(iga_tmpb,iga_cdinv)
#endif
c         second gaget
c
            call ga_zero(iga_2ceri)
            call dft_get2eri(CD_bas_han,iga_2ceri,geom)
          else
            call dft_invdiag(geom,CD_bas_han,iga_2ceri,iga_cdinv,nbf_cd)
          endif
          if(IOLGC.and.(me.eq.0)) then
            call dft_invio('CDI',iga_cdinv,nbf_cd,'WRITE',lmiss)
            call dft_invio('CD',iga_2ceri,nbf_cd,'WRITE',lmiss)
          endif
          call ga_sync
        endif
        if(IOLGC) then
          call gadest(iga_cdinv)
          call gadest(iga_2ceri)
        endif
        call ga_sync
      endif
c     

#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(2)='  CDinv   '
      call tcllct(t_pgm(1,2),tp2,tp1)
      call dft_second(tp1)
#endif
c     
c     Compute the AO overlap matrix.
c     
      call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'ovl ',iga_s,'atom')
      LResult = 
     &MA_Push_Get(MT_Dbl,nbf_ao,'ovl eig vals',lsvals,isvals)
      call ga_zero(iga_s)
      call int_1e_ga(AO_bas_han,AO_bas_han,
     &     iga_s,'overlap',.false.)

c     
#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(3)='  Overl   '
      call tcllct(t_pgm(1,3),tp2,tp1)
      call dft_second(tp1)
#endif
c     
c     orthog BS
c     
c     
c--   > Diag Overlap Matrix
c     
      call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'ovl eig vecs',iga_svecs,'atom')
#if defined(PARALLEL_DIAG)         
      call ga_diag_std(iga_s,iga_svecs,Dbl_MB(isvals))
#else
      call ga_diag_std_seq(iga_s,iga_svecs,Dbl_MB(isvals))
#endif
c     
#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(4)='  Xortho   '
      call tcllct(t_pgm(1,4),tp2,tp1)
      call dft_second(tp1)
#endif
c     
      call ga_sync
      call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'alpha evecs',iga_evecs(1),'atom')
      call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'alpha density',iga_dens(1),'atom')
      call ga_zero(iga_dens(1))
      if( ipol.ne.1 )then
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'beta density',iga_dens(2),'atom')
        call ga_zero(iga_dens(2))
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'beta evecs',iga_evecs(2),'atom')
        call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao,'olddens',iga_olddens,'atom')
        call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'Vxcs',iga_Vxc(2),'atom')
      end if 
      nVxc=ipol
      call gacrea(geom,AO_bas_han,
     &     nbf_ao,nbf_ao,'fock Matrix',iga_fock,'atom')
      iga_vxc(1)=iga_fock
      if(icvgtyp.eq.2) 
     &       call gacrea(geom,AO_bas_han,
     &       nbf_ao,nbf_ao*(2*nfock*ipol),
     &     'diis bufs',iga_diis,'row')
c     
#if defined(DFT_TIMER)
      call dft_second(tp2)
      pgmsect(5)='  CREA     '
      call tcllct(t_pgm(1,5),tp2,tp1)
#endif
c     
c     
c     Begin iteration sequence with S**(-1/2) in iga_evecs arrays
c     
      iw=2
      iga_tmp=iga_s
      call diis_builds12(
     &     isvals,iga_svecs,iga_evecs(1),iga_tmp,nbf_ao,iw)

      if (ipol.gt.1) call ga_copy(iga_evecs(1),iga_evecs(2))
      if(icvgtyp.ne.2) call gadest(iga_svecs)

      if(CDFIT)call gadest(iga_s)

      call dft_scf(rtdb,geom,AO_bas_han,CD_bas_han,XC_bas_han,
     &     ipol,Etold,Enuc,iVcoul_opt,iVxc_opt,
     &     iconst,nconst,iter,
     &     iga_dens,iga_evecs,iga_olddens,iga_vxc,
     &     iga_fock,iga_diis,iga_svecs,isvals,
     &     iga_cdinv,iga_xcinv,iga_2ceri)
c     
c     Clean up memory usage.
c     
      LResult = MA_Pop_Stack(lsvals)
      if(.not.CDFIT)call gadest(iga_s)
      if(me.eq.0) then
        if (CDFIT.and.IOLGC) then
          call dft_invio('CD',idum,1,'KILL',lmiss)
          call dft_invio('CDI',idum,1,'KILL',lmiss)
        endif
        if (XCFIT.and.IOLGC) then
          call dft_invio('XC',idum,1,'KILL',lmiss)
        endif
      endif
      if (CDFIT.and.(.not.IOLGC)) then
        call gadest(iga_cdinv)
        call gadest(iga_2ceri)
      endif
      if (XCFIT.and.(.not.IOLGC)) call gadest(iga_xcinv)
      if(icvgtyp.eq.2)call gadest(iga_diis)
      call gadest(iga_fock)
      if(ipol.eq.2)
     &     call gadest(iga_vxc(2))
      if(ipol.ne.1)then
        call gadest(iga_olddens)
        call gadest(iga_evecs(2))
        call gadest(iga_dens(2))
      endif
      call gadest(iga_evecs(1))
      call gadest(iga_dens(1))
      if(icvgtyp.eq.2) call gadest(iga_svecs)
      call schwarz_tidy
      call ga_sync
c     
#if defined(DFT_TIMER)
      if(me.eq.0)then
        ltime = LuOut
        write(ltime,*)
      endif  
      pgmname=' SCFDRIVER'
      call t_stat(10,2)
#endif
c     
#if defined(DFT_GASTAT)
      if(me.eq.0)call gadump
#endif
c     
      return
      end
