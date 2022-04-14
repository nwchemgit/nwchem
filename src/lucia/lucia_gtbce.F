************************************************************************
      subroutine lucia_gtbce(irefspace,itrefspc,maxit_gtbce)
************************************************************************
*
* Master routine for Generalized Two-Body operator Cluster Expansion,
* i.e. CC  expansions which allow excitations of rank +2,+1,0,-1,-2
*
************************************************************************
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'clunit.inc'
      include 'glbbas.inc'
      include 'cands.inc'
      include 'cecore.inc'
      include 'cc_exc.inc'
      include 'cicisp.inc'
      include 'cintfo.inc'
      include 'gtbce.inc'
      include 'lucinp.inc'
      include 'csmprd.inc'
      include 'multd2h.inc'
      include 'frorbs.inc'

************************************************************************
      integer*8 n_ci_det
      character*8 cctype
      dimension icascr(mxpngas)
      dimension ioccun(100)
      dimension ioff(8)
************************************************************************

      ntest = 5
      i_spin_adapt = 0
      i_do_newccv=0
      icc_exc = 0
      irefspc = irefspace
      
c now defined via lucia.f or input:
c      expg_thrsh = 1d-20
c      mxterm_expg = 200

      if (ntest.ge.5) then
        write(6,'(/,2(x,a,/),2(x,a,i3,/),x,a,3i2,/x,a,5i2)')
     &    'Generalized Two-Body operator Cluster Expansion at work',
     &    '=======================================================',
     &    ' reference space          = ', irefspc,
     &    ' space defining operators = ', itrefspc,
     &    ' SING = ', inc_sing(1:3),
     &    ' DOUB = ', inc_doub(1:5)
        if (i_mode_gtbce.eq.0) then
          write(6,'(2x,a,/,2x,a,i6)')
     &         'Trying to solve the Nakasuji-equations',
     &         'max. iterations = ',maxit_gtbce
        else if (i_mode_gtbce.eq.1) then
          write(6,'(2x,a,/,2x,a,i6)')
     &         'Trying to minimize the E-expectation value directly',
     &         'max. iterations = ',maxit_gtbce
        else
          write(6,'(2x,a,i3)')
     &         'Beeing a bit surprised by i_mode_gtbce = ', i_mode_gtbce
        end if
        if (igtbmod.eq.0) write(6,*) '***  exp(G) expansion  ***'
        if (igtbmod.eq.1) write(6,*) '***  exp(G^2) expansion  ***'
        if (igtbmod.eq.2) write(6,*)
     &       '***  exp(G) expansion with G=LL  ***'
        if (igtbmod.eq.3) write(6,*)
     &       '***  exp(G) expansion with G=U Omega U  ***'
      end if


* currently hard wired specifications:
      ionly_excop = 0
      i_ign_ovl = 1
c      icexc_rank_min = -4   ! allow rank -2 to +2 (multipl. by 2)
c      inc_sing = (/0, 0, 0/)
c      inc_doub = (/1, 1, 1, 1, 1/)
      ! Hermitian or unitary operators requested? -
      if (igtbcs.eq.-1.or.igtbcs.eq.+1) then
        ! --> then generate only operators with positive rank
        inc_sing(3) = 0
        inc_doub(4:5) = 0
      end if
      mn_crea = 1
      mn_anni = 1

* set mark in memory manager
      call memman(idum,idum,'MARK  ',idum,'GTBCE ')
* expand reference wave-function to FCI space
      isym = IREFSM
      icopy = 1
      lblk = -1

      call expciv(isym,irefspc,luc,
     &            itrefspc,lusc1,lblk,
     &            lusc2,nroot,icopy,idc,ntest)

      ! regenerate Fock matrix
      call copvec(work(kint1o),work(kint1),nint1)
      icc_exc = 0
      i_use_simtrh = 0
      call fi(work(kint1),eccx,1)
      ecore = ecore_ini

* define C and Sigma space for mv7 and sigden_cc
      icsm  = irefsm    ! the symmetry of the reference ...
      issm  = irefsm    ! ... is also the symmetry of the wavefunction
      icspc = itrefspc
      isspc = itrefspc

      call cc_ac_spaces(irefspc,ireftyp)
      nael = nelec(1)
      nbel = nelec(1)
      
      iadd_uni = 0

      call gen_ic_orbop2(1,nobex_tp,idummy,
     &              inc_sing,inc_doub,
     &              ionly_excop,i_ign_ovl,
     &              irefspc,itrefspc,iadd_uni)
*. and the orbital excitations
      call memman(kobex_tp,2*ngas*nobex_tp,'ADDL ',2,'GTBOBX')
      klobex = kobex_tp
      call gen_ic_orbop2(2,nobex_tp,work(kobex_tp),
     &              inc_sing,inc_doub,
     &              ionly_excop,i_ign_ovl,
     &              irefspc,itrefspc,iadd_uni)
      nobex_tpe = nobex_tp+1
*
      if(i_spin_adapt.eq.1) then
*
*. excitation operators will be spin adapted
*
        do jobex_tp = 1, nobex_tp
          write(6,*) ' constructing ca confs for jobex_tp = ', jobex_tp
*. integer arrays for creation and annihilation part 
          call icopve2(work(kobex_tp),1+(jobex_tp-1)*2*ngas,2*ngas,
     &                  icascr)
          nop_c = ielsum(icascr,ngas)
          nop_a = ielsum(icascr(1+ngas),ngas)
          nop_ca = nop_c + nop_a
          call get_ca_conf_for_orbex(icascr,icascr(1+ngas),
     &         ncoc_fsm(1,jobex_tp),naoc_fsm(1,jobex_tp),
     &         ibcoc_fsm(1,jobex_tp),ibaoc_fsm(1,jobex_tp),
     &         kcoc(jobex_tp),kaoc(jobex_tp),
     &         kzc(jobex_tp),kza(jobex_tp),
     &         kcreo(jobex_tp),kareo(jobex_tp))
          write(6,*) ' ncoc_fsm and naoc_fsm after get_ca ... '
          call iwrtma(ncoc_fsm,1,nsmst,1,nsmst)
          call iwrtma(naoc_fsm,1,nsmst,1,nsmst)
         
*. offsets in ca block for given symmetry of creation occ
c ioff_symblk_mat(nsmst,na,nb,itotsm,ioff,irestrict
          call ioff_symblk_mat(nsmst,ncoc_fsm(1,jobex_tp),
     &         naoc_fsm(1,jobex_tp),1,ibcaoc_fsm(1,jobex_tp),0)
c                           ndim_1el_mat(ihsm,nrpsm,ncpsm,nsm,ipack)
          ncaoc(jobex_tp) = ndim_1el_mat(1,ncoc_fsm(1,jobex_tp),
     &                      naoc_fsm(1,jobex_tp),nsmst,0)
*. and the actual configurations 
          call memman(kcaoc(jobex_tp),nop_ca*ncaoc(jobex_tp),'ADDL  ',
     &                2,'CA_OC ')
c     get_conf_for_orbex(ncoc_fsm,naoc_fsm,icoc,iaoc,
c    &           nop_c,nop_a, ibcoc_fsm,ibaoc_fsm,nsmst,iopsm,
c    &           icaoc)
          call get_conf_for_orbex(
     &         ncoc_fsm(1,jobex_tp),naoc_fsm(1,jobex_tp),
     &         work(kcoc(jobex_tp)),work(kaoc(jobex_tp)),
     &         nop_c, nop_a,
     &         ibcoc_fsm(1,jobex_tp),ibaoc_fsm(1,jobex_tp),
     &         nsmst,1,work(kcaoc(jobex_tp)) )
        end do
      end if ! i_spin_adapt
*. number of creation and annihilation operators per op
      call memman(klcobex_tp,nobex_tpe,'ADDL ',1,'COBEX ')
      call memman(klaobex_tp,nobex_tpe,'ADDL ',1,'AOBEX ')
      call get_nca_for_orbop(nobex_tp,work(kobex_tp),
     &     work(klcobex_tp),work(klaobex_tp),ngas)
*. number of spinorbital excitations
      izero = 0
      mxspox = 0
      iact_spc = 0
      iaaexc_typ = 3
      irefspcx = 0
      call obex_to_spobex2(1,work(kobex_tp),work(klcobex_tp),
     &     work(klaobex_tp),nobex_tp,idummy,nspobex_tp,ngas,
     &     nobpt,0,izero ,iaaexc_typ,iact_spc,iprcc,idummy,
     &     mxspox,work(knsox_for_ox),
     &     work(kibsox_for_ox),work(kisox_for_ox),nael,nbel,irefspcx,
     &     mn_crea,mn_anni)
      nspobex_tpe = nspobex_tp + 1
*. and the actual spinorbital excitations
      call memman(klsobex,4*ngas*nspobex_tpe,'ADDL  ',1,'SPOBEX')
*. map spin-orbital exc type => orbital exc type
      call memman(klsox_to_ox,nspobex_tpe,'ADDL  ',1,'SPOBEX')
*. first sox of given ox ( including zero operator )
      call memman(kibsox_for_ox,nobex_tpe,'ADDL  ',1,'IBSOXF')
*. number of sox's for given ox
      call memman(knsox_for_ox,nobex_tpe,'ADDL  ',1,'IBSOXF')
*. sox for given ox
      call memman(kisox_for_ox,nspobex_tpe,'ADDL  ',1,'IBSOXF')
*
      call obex_to_spobex2(2,work(kobex_tp),work(klcobex_tp),
     &     work(klaobex_tp),nobex_tp,work(klsobex),nspobex_tp,ngas,
     &     nobpt,0,mscomb_cc,iaaexc_typ,iact_spc,iprcc,
     &     work(klsox_to_ox),mxspox,work(knsox_for_ox),
     &     work(kibsox_for_ox),work(kisox_for_ox),nael,nbel,irefspcx,
     &     mn_crea,mn_anni)
*
*
      write(6,*) 'Generated excitations:'
      write(6,*) '======================'
      call wrt_spox_tp(work(klsobex),nspobex_tp) 

*
* dimension and offsets of ic operators
*
      call memman(kllsobex,nspobex_tpe,'ADDL  ',1,'LSPOBX')
      call memman(klibsobex,nspobex_tpe,'ADDL  ',1,'LSPOBX')
      call memman(klspobex_ac,nspobex_tpe,'ADDL  ',1,'SPOBAC')
*. all spinorbital excitations are initially active
      ione = 1
      call isetvc(work(klspobex_ac),ione,nspobex_tpe)
*
      itop_sm = 1
      write(6,*) ' irefspc before idim.. ', irefspc
      call idim_tcc(work(klsobex),nspobex_tp,itop_sm,
     &     mx_st_tsoso,mx_st_tsoso_blk,mx_tblk,
     &     work(kllsobex),work(klibsobex),len_t_vec,
     &     mscomb_cc,mx_tblk_as,
     &     work(kisox_for_occls),noccls,work(kibsox_for_occls),
     &     ntconf,iprcc)

      ! set up nfrobs
      call set_frobs(nfrob,nfrobs)
      ntaobs(1:nsmob) = ntoobs(1:nsmob)-nfrobs(1:nsmob)
      ntaob = ntoob-nfrob

      n_cc_amp = len_t_vec
      n_ci_det = xispsm(1,itrefspc)
      n_ci_csf = ncsf_for_cispace(itrefspc,irefsm)

      i12loc = 1
      i34loc = 1
      i1234loc = 1
      
      imode = 0
      call pnt4dm2(nh2elm_p11,imode,
     &     nsmob,nsmsx,mxpobs,ntaobs,ntaobs,ntaobs,ntaobs,
     &     itsdx,adsxa,sxdxsx,i12loc,i34loc,i1234loc,
     &     idum,idum,adasx)
      
      i12loc = 1
      i34loc = 1
      i1234loc = -1
      
      imode = 0
      call pnt4dm2(nh2elm_m11,imode,
     &     nsmob,nsmsx,mxpobs,ntaobs,ntaobs,ntaobs,ntaobs,
     &     itsdx,adsxa,sxdxsx,i12loc,i34loc,i1234loc,
     &     idum,idum,adasx)
      
      i12loc = -1
      i34loc = -1
      i1234loc = 1
      
      imode = 0
      call pnt4dm2(nh2elm_p33,imode,
     &     nsmob,nsmsx,mxpobs,ntaobs,ntaobs,ntaobs,ntaobs,
     &     itsdx,adsxa,sxdxsx,i12loc,i34loc,i1234loc,
     &     idum,idum,adasx)
      
      i12loc = -1
      i34loc = -1
      i1234loc = -1
      
      imode = 0
      call pnt4dm2(nh2elm_m33,imode,
     &     nsmob,nsmsx,mxpobs,ntaobs,ntaobs,ntaobs,ntaobs,
     &     itsdx,adsxa,sxdxsx,i12loc,i34loc,i1234loc,
     &     idum,idum,adasx)
      
c      call num_ssaa2op(nndiag,ndiag)

      write(6,*)
     &     '======================================================'
      write(6,'(x,a,i20)')
     &     ' number of amplitudes:          ',
     &     n_cc_amp
      write(6,'(x,a,2(/,x,a,2i20))')
     &     ' number of indep. two-body parameters',
     &     ' eff H (non-diagonal/diagonal): ',
     &     nh2elm_m11,nh2elm_p11-nh2elm_m11,
     &     ' G (non-diagonal/diagonal):     ',
     &     nh2elm_m11+nh2elm_m33,nh2elm_p11+nh2elm_p33
     &                          -nh2elm_m11-nh2elm_m33
      write(6,'(x,a,/,x,a,i20)')
     &     ' number of determinants/combinations ',
     &     ' in the underlying CI-Space:    ',
     &     n_ci_det
      write(6,'(x,a,/,x,a,i20)')
     &     ' number of CSFs ',
     &     ' in the underlying CI-Space:    ',
     &     n_ci_csf
      write(6,*)
     &     '======================================================'
      
      if (nh2elm_m11+nh2elm_m33.gt.n_ci_csf) then
        write(6,*)
     &  ' Well, the number of non-linear parameters is larger than the'
        write(6,*)
     &  ' the number of CI-paramters! This appears rather silly to me!'
        do ii = 1, 30
          write(6,*) '    ???????????? silly calculation ????????????'
        end do
      end if

      write(6,*) ' dimension of the various types '
      call iwrtma(work(kllsobex),1,nspobex_tp,1,nspobex_tp)
      write(6,*) ' offsets of the various types '
      call iwrtma(work(klibsobex),1,nspobex_tp,1,nspobex_tp)
*
      mx_st_tsoso_mx = mx_st_tsoso
      mx_st_tsoso_blk_mx = mx_st_tsoso_blk
      mx_tblk_mx = mx_tblk
      mx_tblk_as_mx = mx_tblk_as
      len_t_vec_mx =  len_t_vec
*. some more scratch etc
*. alpha- and beta-excitations constituting the spinorbital excitations
*. number 

      call spobex_to_abobex(work(klsobex),nspobex_tp,ngas,
     &     1,naobex_tp,nbobex_tp,idummy,idummy)
*. and the alpha-and beta-excitations
      lena = 2*ngas*naobex_tp
      lenb = 2*ngas*nbobex_tp
      call memman(klaobex,lena,'ADDL  ',2,'IAOBEX')
      call memman(klbobex,lenb,'ADDL  ',2,'IAOBEX')
      call spobex_to_abobex(work(klsobex),nspobex_tp,ngas,
     &     0,naobex_tp,nbobex_tp,work(klaobex),work(klbobex))
*. max dimensions of ccop !kstr> = !istr> maps
*. for alpha excitations
      iatp = 1
      ioctpa = ibspgpftp(iatp)
      noctpa = nspgpftp(iatp)
      call len_genop_str_map(
     &     naobex_tp,work(klaobex),noctpa,nelfspgp(1,ioctpa),
     &     nobpt,ngas,maxlena)
      ibtp = 2
      ioctpb = ibspgpftp(ibtp)
      noctpb = nspgpftp(ibtp)
      call len_genop_str_map(
     &     nbobex_tp,work(klbobex),noctpb,nelfspgp(1,ioctpb),
     &     nobpt,ngas,maxlenb)
      maxlen_i1 = max(maxlena,maxlenb)
      if(ntest.ge.5) write(6,*) ' maxlen_i1 = ', maxlen_i1

c get work space: 
c get dimensions for FCI (wow) behind the curtains
      call get_3blks_gcc(kvec1,kvec2,kvec3,mxcj)
      kc2=kvec3
      write(6,*) 'max block length from get_3blks: ', mxcj
*. and two CC vectors
c      n_sd_int = 1
      lenny = n_cc_amp ! + n_sd_int
      call memman(kcc1,lenny,'ADDL  ',2,'CC1_VE')
      call memman(kcc2,lenny,'ADDL  ',2,'CC2_VE')
*
      if (igtbcs.eq.1.or.igtbcs.eq.-1.or.isymmet_G.ne.0) 
     &    call memman(kcc3,lenny,'ADDL  ',2,'CC3_VE')
      if (isymmet_G.ne.0) 
     &    call memman(kiccvec,lenny,'ADDL  ',1,'ICCVEC')
*. and the cc diagonal 
      if (igtbmod.eq.2.or.igtbmod.eq.3) lenny = max(lenny,(2*ntoob)**2)
      call memman(kdia,lenny,'ADDL  ',2,'CC_DIA')

      if (igtbmod.lt.2) then
        imod = 1  ! Fock-matrix based on rho1
        call gencc_f_diag_m(imod,work(klsobex),nspobex_tp,work(kdia),1,
     &      xdum,idum,idum,0,
     &      work(kvec1),work(kvec2),mx_st_tsoso_mx,
     &      mx_st_tsoso_blk_mx)
c the approximate Hessian is two times the diagonal!
        call scalve(work(kdia),2d0,n_cc_amp)
        ! well, at the moment I do not know better than removing
        ! all negative and small stuff:
        if (isymmet_G.ne.0) then
          do ii = 1, n_cc_amp
            work((kdia-1)+ii) = abs(work((kdia-1)+ii))
          end do      
        end if
        xmin = 100d0
        do ii = 1, n_cc_amp
          xmin = min(xmin,work((kdia-1)+ii))
        end do
        write(6,*) 'diagonal: lowest element = ',xmin
        xsh = max(0d0,0.01d0-xmin)
        write(6,*) 'shift diagonal by ',xsh
        do ii = 1, n_cc_amp
          work((kdia-1)+ii) = work((kdia-1)+ii) + xsh
        end do      
        
        if (igtb_closed.eq.0) then
          call vec_to_disc(work(kdia),n_cc_amp,1,lblk,ludia)
        else

          call memman(kpamp, 2*nsmob**3,'ADDL  ',1,'PSMTR ')
          call memman(kpamp2,2*nsmob**3,'ADDL  ',1,'PSMTR2')

          call setup4idx(isymmet_G,n11amp,n33amp,
     &                   work(kpamp),work(kpamp2),ntaobs)

          namp_packed = n11amp + n33amp
         
c TESTING
c          work(kcc3:kcc3-1+namp_packed) = 0d0
c          idx = 0
cc          do isymq = 1, nsmob
cc            do isymp = 1, isymq
cc              isymrs = multd2h(isymp,isumq)
cc              do isymr = 1, nsmob
cc                isyms = multd2h(isymr,isymrs)
c              
c          do idxs = nfrob+1, ntoob 
c            do idxr = nfrob+1, ntoob
c              do idxq = nfrob+1, ntoob
c                do idxp = nfrob+1, ntoob
cc                  idxsr = (idxs-1)*ntoob + idxr
cc                  if (idxpq.ge.idxrs) cycle
cc                  if (idxpq.gt.idxsr) cycle
c                  iadr = i2addr2(   idxp,idxq,idxr,idxs,
c     &                              work(kpamp),1,1,-1)
c                  if (iadr.eq.-2) cycle
c                  idx = idx+1
c                  
c                  print *,'-----------------------------------------'
c                  print '(a,i5,a,4i5)',
c     &                 ' INDEX: ',idx,'  ',idxp,idxq,idxr,idxs
c                  print *,'SINGLET-SINGLET'
c                  iadr = i2addr2(   idxp,idxq,idxr,idxs,
c     &                              work(kpamp),1,1,-1)
c                  print '(4i5,a,i5,x,"S")',
c     &                 idxp,idxq,idxr,idxs,' --> ',iadr
c                  if (iadr.lt.1.or.iadr.gt.n11amp) then
c                    print *,'EVIL RANGE ERROR: ',1,iadr,namp_packed
c                  else
c                    if (work(kcc3-1+iadr).eq.0d0) then
c                      work(kcc3-1+iadr) = dble(idx)
c                    else
c                      print *,'EIEIEI, wer hat auf meinem Plaetzchen'//
c     &                     ' gesessen?',
c     &                     work(kcc3-1+iadr)
c                    end if
c                  end if
c
c                  print *,'TRIPLET-TRIPLET'
c                  iadr = i2addr2(   idxp,idxq,idxr,idxs,
c     &                              work(kpamp+nsmob**3),-1,-1,-1)
c                  print '(4i5,a,i5,x,"T")',
c     &                 idxp,idxq,idxr,idxs,' --> ',iadr
c                  if (iadr.lt.1.or.iadr.gt.n33amp) then
c                    print *,'EVIL RANGE ERROR: ',1,iadr,namp_packed
c                  else
c                    iadr = iadr+n11amp
c                    if (work(kcc3-1+iadr).eq.0d0) then
c                      work(kcc3-1+iadr) = dble(idx)
c                    else
c                      print *,'EIEIEI, wer hat auf meinem Plaetzchen'//
c     &                     ' gesessen?',
c     &                     work(kcc3-1+iadr)
c                    end if
c                  end if
c
cc                  if (idxp.ne.idxq) then
cc                    print *,'+ INDEX: ',idx
cc                    idx = idx+1
cc                    iadr3 = i2addr2(   idxp,idxq,idxs,idxr,
cc     &                   work(kpamp),1,0,-1)
cc                    print '(4i5,a,i5)',idxp,idxq,idxs,idxr,' --> ',iadr3
cc                    if (iadr3.lt.1.or.iadr3.gt.namp_packed)
cc     &                 print *,'RANGE ERROR: ',1,iadr3,namp_packed
cc                  end if
c
cc                  iadr1 = i2addr2(   idxr,idxs,idxp,idxq,
cc     &                              work(kpamp),1,0,-1)
cc                  print '(4i5,a,i5)',idxr,idxs,idxp,idxq,' --> ',iadr1
cc                  if (iadr1.ne.iadr)
cc     &                 print *,'SYM. ERROR'
cc
cc                  iadr2 = i2addr2(   idxq,idxp,idxr,idxs,
cc     &                              work(kpamp),1,0,-1)
cc                  print '(4i5,a,i5)',idxq,idxp,idxr,idxs,' --> ',iadr2
cc                  if (iadr2.lt.1.or.iadr2.gt.namp_packed)
cc     &                 print *,'RANGE ERROR: ',1,iadr2,namp_packed
cc
cc                  iadr4 = i2addr2(   idxq,idxp,idxs,idxr,
cc     &                              work(kpamp),1,0,-1)
cc                  print '(4i5,a,i5)',idxq,idxp,idxs,idxr,' --> ',iadr4
cc                  if (iadr4.ne.iadr)
cc     &                 print *,'SYM. ERROR'
cc
c                end do
c              end do
c            end do
c          end do
c
cc              end do
cc            end do
cc          end do
c          print *,'-----------------------------------------'
c
cc          call wrtmat(work(kcc3),namp_packed,1,namp_packed,1)
c          do ii = 1, namp_packed
c            if (work(kcc3-1+ii).eq.0d0) then
c              print *,ii,work(kcc3-1+ii),' <--'
c            else
c              print *,ii,work(kcc3-1+ii)
c            end if
c          end do
c
c          stop 'testing'
c TESTING


          iway = 1 ! pack (no symmetrizing, would result in 0d0's)
          idual = 3
          call pack_g(iway,idual,isymmet_G,work(kcc1),work(kdia),
     &         nspobex_tp,work(klsobex),work(klibsobex),n11amp,n33amp,
     &         work(kpamp),n_cc_amp)
          call vec_to_disc(work(kcc1),namp_packed,1,lblk,ludia)
        end if

        if (ntest.gt.100) then
          write(6,*) 'the preconditioner: '
          cctype='GEN_CC'
          call wrt_cc_vec2(work(kdia),6,cctype)
        end if

      else if (igtmode.eq.2) then
c some init for G=LL        

        idx = 0
        do ism = 1, nsmob
          ioff(ism) = idx
          idx = idx + (ntoobs(ism)+1)*ntoobs(ism)/2
        end do
        
        do ii = 1, ntoob
          do jj = 1, ntoob
            ism = ismfto(ii)
            jsm = ismfto(jj)
            idx = ireots(ii) - ibso(ism) + 1
            jdx = ireots(jj) - ibso(jsm) + 1

            iidx = ioff(ism) + (idx+1)*idx/2
            jjdx = ioff(jsm) + (jdx+1)*jdx/2
            
            work(kdia-1+(ii-1)*ntoob+jj) =
     &           work(kfiz-1+iidx)-work(kfiz-1+jjdx)

            print *,ii,jj,'->',work(kfiz-1+iidx),work(kfiz-1+jjdx)

          end do
        end do

        do ii = 1, ntoob**2
          work(kdia-1+ii) = max(.1d0,work(kdia-1+ii))
        end do
        call vec_to_disc(work(kdia),ntoob**2,1,lblk,ludia)
      else if (igtbmod.eq.3) then

        ! get memory for G= U Om U variant
        nlen = ntoob**2*4
        call memman(komvec,nlen,'ADDL  ',2,'OMVEC ')
        call memman(kurvec,nlen,'ADDL  ',2,'URVEC ')
        call memman(kuivec,nlen,'ADDL  ',2,'UIVEC ')
        call memman(komgrd,nlen,'ADDL  ',2,'OMGRD ')
        call memman(kurgrd,nlen,'ADDL  ',2,'URGRD ')
        call memman(kuigrd,nlen,'ADDL  ',2,'UIGRD ')

      end if

      i_test_fock = 0

      if (i_test_fock.ne.1) then

        call gtbce_opt(maxit_gtbce,irefspc,itrefspc,
     &               work(kcc1),work(kcc2),work(kdia),work(kcc3),
     &               work(kvec1),work(kvec2),work(kc2),
     &               nspobex_tp,work(klsobex),
     &               work(kllsobex),work(klibsobex),
     &               igtbcs,mxcj,
     &               n11amp,n33amp,work(kpamp),
     &               work(komvec),work(kurvec),work(kuivec),
     &               work(komgrd),work(kurgrd),work(kuigrd),
     &               work(kiccvec),
     &               luc,lu_ccamp,lu_ccvecf,ludia,
     &               lusc3,luhc)

      else

        call gucc_fock(irefspc,itrefspc,
     &       work(kcc1),work(kcc2),work(kdia),work(kcc3),
     &       work(kvec1),work(kvec2),work(kc2),
     &       nspobex_tp,work(klsobex),
     &       work(kllsobex),work(klibsobex),
     &       igtbcs,mxcj,
     &       luc,lu_ccamp,lu_ccvecf,ludia,
     &       lusc3,luhc)

      end if

c TESTING: copy exp(G)|ref> to |ref>
      call copvcd(lusc3,luc,work(kvec1),1,lblk)

      call memman(idum,idum,'FLUSM ',idum,'GTBCE ')

      return
      
      end
************************************************************************
      subroutine setup4idx(isymmet_G,n11amp,n33amp,
     &                     ioff_amp,isy_amp,ntaobs)
*     little slave routine to address parts of work(kpamp),
*     the curse of using self-made allocation
*     routines ....

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'lucinp.inc'
      include 'csm.inc'
      include 'csmprd.inc'

      dimension ioff_amp(nsmob*nsmob*nsmob,2)
      dimension isy_amp(nsmob*nsmob*nsmob,2)
      dimension ntaobs(*)

      ! singlet-singlet amplitudes
      i12loc = 1
      i34loc = 1      
      i1234loc = isymmet_G      ! antisymmetry between 12 and 34

      imode = 1
      call pnt4dm2(n11amp,imode,
     &     nsmob,nsmsx,mxpobs,ntaobs,ntaobs,ntaobs,ntaobs,
     &     itsdx,adsxa,sxdxsx,i12loc,i34loc,i1234loc,
     &     ioff_amp(1,1),isy_amp(1,1),adasx)

      ! triplet-triplet amplitudes
      i12loc = -1
      i34loc = -1      
      i1234loc = isymmet_G      ! antisymmetry between 12 and 34
      
      imode = 1
      call pnt4dm2(n33amp,imode,
     &     nsmob,nsmsx,mxpobs,ntaobs,ntaobs,ntaobs,ntaobs,
     &     itsdx,adsxa,sxdxsx,i12loc,i34loc,i1234loc,
     &     ioff_amp(1,2),isy_amp(1,2),adasx)
      
      return
      end
************************************************************************
* DECK: gtbce_opt
************************************************************************
      subroutine gtbce_opt(maxiter,irefspc,itrefspc,
     &                     ccvec1,ccvec2,ccvec3,ccvec4,
     &                     civec1,civec2,c2vec,
     &                     n_cc_typ,i_cc_typ,
     &                     namp_cc_typ,ioff_cc_typ,
     &                     iopsym,mxb_ci,
     &                     n11amp,n33amp,iamp_packed,
     &                     omvec,urvec,uivec,
     &                     omgrd,urgrd,uigrd,
     &                     iccvec,
     &                     luc,luamp,luomg,ludia,
     &                     luec,luhc)
************************************************************************
*
* purpose : driver for the optimization of the Generalize TwoBody
*           operator Cluster Expansion wavefunction (if it works at all)
*
*  ak, early 2004
*
************************************************************************
*
* units:
*   luc   = definition of reference function
*   luamp = amplitude vectors (also output for most recent vector)
*   luampold = scratch containing old vectors from previous iterations
*           (on input it may also be a first trial vector)
*   luomg = error vectors
*   ludia = diagonal preconditioner
*   luec  = scratch for exp(G)|ref>
*   luhc  = scratch for H exp(G)|ref>

* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cecore.inc'
      include 'gtbce.inc'
      include 'opti.inc'
      include 'glbbas.inc'
      include 'cintfo.inc'
* constants
      integer, parameter ::
     &     ntest = 5

* arrays
      integer ::
     &     ioff_cc_typ(n_cc_typ), namp_cc_typ(n_cc_typ),
     &     i_cc_typ(4*ngas,n_cc_typ), iccvec(n_cc_amp)
      real*8 ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp), ccvec3(n_cc_amp),
     &     omvec(ntoob,ntoob,2,2),urvec(ntoob,ntoob,2,2),
     &     uivec(ntoob,ntoob,2,2),
     &     omgrd(ntoob,ntoob,2,2),urgrd(ntoob,ntoob,2,2),
     &     uigrd(ntoob,ntoob,2,2)
* local
      logical ::
     &     calc_Omg, calc_gradE, tstgrad, tst_hss, comm_ops,
     &     do_eag, do_foo, do_hss, do_rdvec, did_rdvec, do_h0
      character*8 cctype
      integer ::
     &     ictp(n_cc_typ)
* external functions
      real*8, external :: inprod, inprdd

* ============================================================
* initialize : restart, set coefs to zero
* ============================================================

      call atim(cpu0,wall0)

      nprint = 1
      lblk = -1

      if (ntest.ge.5) then
        write(6,*) '======================='
        write(6,*) 'entered gtbce_opt with:'
        write(6,*) '======================='
        write(6,*) ' iopsym = ',iopsym
      end if

      calc_gradE = .true.
      calc_Omg   = .true.

* unit init
      lusc1 = iopen_nus('GTBSCR1')
      lusc2 = iopen_nus('GTBSCR2')
      lusc3 = iopen_nus('GTBSCR3')
      lusc4 = iopen_nus('GTBSCR4')
      lusc5 = iopen_nus('GTBSCR5')
      lusc6 = iopen_nus('GTBSCR6')
      lusc7 = iopen_nus('GTBSCR7')
      lusc8 = iopen_nus('GTBSCR8')
      lusc9 = iopen_nus('GTBSCR9')

      luhss = iopen_nus('GTBHESS')
      luh0  = iopen_nus('GTBH0')
      lufoo = iopen_nus('GTBFOO')

      lutrvec = iopen_nus('GTBTRVC')
      lusig   = iopen_nus('GTBSIG')

      lurdvec = iopen_nus('GTBRDVEC')

      ! our functional is variational:
      ivar = 1
      if (igtbfusk.gt.50) then
c preliminary: set common opti
        iorder = 1
        iprecnd = 1
        isubsp = 1
        ilsrch = 2
        icnjgrd = 1
        mxsp_sbspja = 10
        isbspjatyp = 1
        isbspja_start = 2       ! lowest possible iteration!
        thr_sbspja = 1d-1
        mxsp_diis = 10
        idiistyp = 2
        idiis_start = 0
        thr_diis = 1d-1
        trini = 0.140d0
        trmin = 0.025d0
        trmax = 0.5d0
        trthr1l = 0.8d0
        trthr1u = 1.2d0
        trthrfac1 = 1.2d0
        trthr2l = 0.4d0
        trthr2u = 1.6d0
        trfac1  = 1.2d0
        trfac2  = 0.8d0
        trfac3  = 0.3d0
        thrstp  = 1d-5
        thrgrd  = 1d-5
        thr_de  = 1d-8
      end if
      maxmacit = maxiter
      micifac  = 20
      maxmicit = maxmacit*micifac
*
      if (igtbmod.eq.2) then
        len = ntoob*ntoob ! very simple, to be adapted for frozen orbitals
        n_l_amp = len
c        if ((len*len+1)/2.gt.n_cc_amp) then
c          write(6,*) ' ',(len+1)*len/2,' <---> ',n_cc_amp 
c          write(6,*) 'input not appropriate for this test!'
c          stop 'ihtest'
c        end if
        call memman(khvec1,len,'ADDL  ',2,'HTEST1 ')
        call memman(khvec2,len,'ADDL  ',2,'HTEST2 ')
        do ii = 1, 30
          write(6,*) ' !!!!!!!!!!! G = LL test active !!!!!!!!!!!'
        end do
      end if

*
* set initial G
*
      if (igtbmod.eq.0) then
        
c        imode = -1
c for the moment better:
        imode = 1
        luinp = luamp

        nwfpar = n_cc_amp
        if (igtb_closed.eq.1) then
          imode = 1
          namp_packed = n11amp + n33amp
          nwfpar = namp_packed
        end if
        
        call gtbce_initG(ccvec1,
     &                 imode,luinp,
     &                 ccvec2,
     &                 ngas,igsoccx(1,1,itrefspc),
     &                 ihpvgas,nwfpar,i_cc_typ,n_cc_typ,
     &                 namp_cc_typ,ioff_cc_typ)

        if (igtb_disptt.eq.1) then
          write(6,*) ' ACCORDING TO YOUR WISHES I DISPOSE THE '//
     &         'ANTISYMMETRIC PART OF G !!!'
          ccvec1(n11amp+1:n11amp+n33amp) = 0d0
        end if


        ! not necessary for igtb_closed.eq.1
        if (isymmet_G.ne.0
!     &       .and.igtb_closed.eq.0
     &       ) then
          call conj_t_pairs(ictp,ierr,
     &         i_cc_typ,n_cc_typ,ngas)
          if (ierr.ne.0) then
            write(6,*)
     &           'The definition of the G operator is not compatible '//
     &           'with the symmetrizing option!'
            stop 'symmetrizing problem'
          end if
          if (igtb_closed.eq.0) then
            call symmet_t(isymmet_G,1,
     &                  ccvec1,ccvec2,
     &                  ictp,i_cc_typ,n_cc_typ,
     &                  namp_cc_typ,ioff_cc_typ,ngas)
          end if
        end if

c        ! project out redundant components:
c        call prjout_red(ccvec1,ccvec2,nspobex_tp,work(klsobex),
c     &       work(klibsobex))

        call vec_to_disc(ccvec1,nwfpar,1,lblk,luamp)

      else if (igtbmod.eq.1) then
        imode = -1
        luinp = luamp
        ! well, at the moment there are problems, so ...
          imode = 0
        
        call gtbce_initG(ccvec1,
     &                 imode,luinp,
     &                 ccvec2,
     &                 ngas,igsoccx(1,1,itrefspc),
     &                 ihpvgas,n_cc_amp,i_cc_typ,n_cc_typ,
     &                 namp_cc_typ,ioff_cc_typ)

        if (igtbfusk.ge.5) then
          call memman(kcan,2*ntoob**2,'ADDS  ',2,'OPCAN ')
* fusk init of operator in canonical symmetry blocked form:

c          work(kcan:kcan+2*ntoob**2-1) = 0.0d0
 
c for testing the gradient         
c          do ii = 1, ntoob**2
c            work(kcan:kcan+2*ntoob**2-1) = 1.d0/(dble(ii)+4d0)
c          end do

c some info on occ/virt orbital per symmetry would be nice here:

cc        ! init for CH2
          ioff = 0
          do ism = 1, nsmst
            if (ism.eq.1) then
              do ii = 2,3
                do jj = 4,7
                  idx = ioff+(ii-1)*ntoobs(ism)+jj
                  work(kcan+idx-1) =  0.2d0
c                  idx = ioff+(jj-1)*ntoobs(ism)+ii
c                  work(kcan+idx-1) = -0.05d0
                end do
              end do
            end if
            if (ism.eq.2) then
              do ii = 1, 1
                do jj = 2,4
                  idx = ioff+(ii-1)*ntoobs(ism)+jj
                  work(kcan+idx-1) =  0.2d0
c                  idx = ioff+(jj-1)*ntoobs(ism)+ii
c                  work(kcan+idx-1) = -0.025d0
                end do
              end do
            end if
            ioff = ioff + ntoobs(ism)*ntoobs(ism)
          end do

c a routine to get from the usual (I called it "canonical") symmetry blocked
c form to LUCIA's string ordering; just for convenience ...
          call can2str(1,work(kcan),ccvec1,
     &         nspobex_tp,i_cc_typ,ioff_cc_typ)

          call vec_to_disc(ccvec1,n_cc_amp,1,lblk,luamp)

        end if

        call vec_to_disc(ccvec1,n_cc_amp,1,lblk,luamp)


      else if (igtbmod.eq.2) then
        if (igtbfusk.ge.10) then
          ! just something but different for each element (for testing purps)
          do ii = 1, n_l_amp
            work(khvec1+ii-1) = 1d0/(dble(ii)+4d0) ! 0d0 
          end do
        else
          ! hm, everything set to a small value:
          do ii = 1, n_l_amp
            work(khvec1+ii-1) = 0.01d0
          end do
        end if

c        ! init for CH2
c        do ii = 2, 4
c          do jj = 5, 8
c            work(khvec1+(ii-1)*ntoob+jj) = 0.05d0
c            work(khvec1+(jj-1)*ntoob+ii) = 0.05d0
c          end do
c        end do

        call vec_to_disc(work(khvec1),n_l_amp,1,lblk,luamp)
c testing
c        call l2g(work(khvec1),ccvec1,nspobex_tp,work(klsobex),0,ntoob)
c        call wrt_cc_vec2(ccvec1,6,'GEN_CC')
c        stop 'brute force'
c testing
      else if (igtbmod.eq.3) then

        nlen = ntoob**2*4

        ! we need three files:
        luom = iopen_nus('OMEGA_VEC')
        luur = iopen_nus('UREAL_VEC')
        luui = iopen_nus('UIMAG_VEC')

        luomgr = iopen_nus('OMEGA_GRD')
        luurgr = iopen_nus('UREAL_GRD')
        luuigr = iopen_nus('UIMAG_GRD')

        ! find out how to set up the preconditioner:::
        call memman(idum,idum,'MARK  ',idum,'LOCAL ')
c        call memman(kfdia,nacob,'ADDL  ',2,'KFDIA ')
c        CALL GT1DIS(WORK(KFDIA),IREOTS,WORK(KPINT1),
c     &            WORK(KFI),ISMFTO,IBSO,NACOB)
c
c        ! Ur with diagonal 1d0
c        do imp = 1,2
c          do imq = 1,2
c            do ip = 1, ntoob
c              do iq = 1, ntoob
c                urvec(iq,ip,imq,imp) =
c     &             abs(2d0*(  work(kfdia + ip) - work(kfdia + iq)))
c                if (urvec(iq,ip,imq,imp).lt.1d-3)
c     &               urvec(iq,ip,imq,imp) = 10d0
c              end do
c            end do
c          end do
c        end do
c
        ! well no, take only 1d0
        urvec(1:ntoob,1:ntoob,1:2,1:2) = 1d0

        call memman(idum,idum,'FLUSM ',idum,'LOCAL ')

        call vec_to_disc(urvec,nlen,1,-1,luom)
        call vec_to_disc(urvec,nlen,1,-1,luur)
        call vec_to_disc(urvec,nlen,1,-1,luui)
        
        imode = 11
        call cmbamp(imode,luom,luur,luui,ludia,
     &       omvec,nlen,nlen,nlen)

        ! try to restart, if file luamp is present
        write(6,*) ' testing unit ',luamp
        rewind(luamp,err=100)
        read(luamp,err=100,end=100) namp_read
        if (namp_read.eq.nlen) then
          imode = 01
          call cmbamp(imode,luom,luur,luui,luamp,
     &       omvec,nlen,nlen,nlen)
          write(6,*) '================='
          write(6,*) ' RESTART SUCCESS'
          write(6,*) '================='
          goto 200
        end if

 100    continue
        ! else: we init

        ! Omega with zero
        omvec(1:ntoob,1:ntoob,1:2,1:2) = 0.d0
c        do im = 1,2
c          do ii = 1, ntoob
c            omvec(ii,ii,im,im) = 1d0
c          end do
c        end do

        ! Ur and Ui with diagonal 1d0
c        urvec(1:ntoob,1:ntoob,1,1) = 1d-3
c        urvec(1:ntoob,1:ntoob,1,2) = 0d0
c        urvec(1:ntoob,1:ntoob,2,1) = 0d0
c        urvec(1:ntoob,1:ntoob,2,2) = 1d-3
c test
        fac = 1d0/sqrt(dble(ntoob))

        urvec(1:ntoob,1:ntoob,1:2,1:2) = 0d0
        do im = 1,2
          do ii = 1, ntoob
c fusk for 1 frozen orbital:
            if (ii.eq.1) then
              urvec(ii,ii,im,im) = fac
            else
              do jj = ii,ntoob
                urvec(jj,ii,im,im) = fac
              end do
            end if
          end do
        end do

c        do im = 1,2
c          do ii = 1, ntoob
c            urvec(ii,ii,im,im) = 1d0/sqrt(2d0)
c          end do
c        end do

c        uivec(1:ntoob,1:ntoob,1,1) = -1d-3
c        uivec(1:ntoob,1:ntoob,1,2) = 0d0
c        uivec(1:ntoob,1:ntoob,2,1) = 0d0
c        uivec(1:ntoob,1:ntoob,2,2) = -1d-3

        uivec(1:ntoob,1:ntoob,1:2,1:2) = 0d0
        do im = 1,2
          do ii = 1, ntoob
c fusk for 1 frozen orbital:            
            if (ii.gt.1) then
              do jj = 1, ii-1
                uivec(jj,ii,im,im) = -fac
              end do
            end if
            uivec(ii,ii,im,im) = fac
          end do
        end do
c        do im = 1,2
c          do ii = 1, ntoob
c            uivec(ii,ii,im,im) = 1d0/sqrt(2d0)
c          end do
c        end do

c fusk        
c        ihom = 2
c        ilum = 3
c        urvec(ihom,ihom,1,1)  = 1d0/2d0
c        urvec(ihom,ilum,1,1) = 1d0/sqrt(2d0)
c        urvec(ilum,ihom,1,1) = 0d0
c        urvec(ilum,ilum,1,1)= 1d0/2d0
c        urvec(ihom,ihom,2,2)  = 1d0/2d0
c        urvec(ihom,ilum,2,2) = 1d0/sqrt(2d0)
c        urvec(ilum,ihom,2,2) = 0d0
c        urvec(ilum,ilum,2,2)= 1d0/2d0
c
c        uivec(ihom,ihom,1,1)  = 1d0/2d0
c        uivec(ihom,ilum,1,1) = 0d0
c        uivec(ilum,ihom,1,1) = -1d0/sqrt(2d0)
c        uivec(ilum,ilum,1,1)= 1d0/2d0
c        uivec(ihom,ihom,2,2)  = 1d0/2d0
c        uivec(ihom,ilum,2,2) = 0d0
c        uivec(ilum,ihom,2,2) = -1d0/sqrt(2d0)
c        uivec(ilum,ilum,2,2)= 1d0/2d0

c        urvec(3,3,1,1)  = 1d0/2d0
c        urvec(3,6,1,1) = 1d0/sqrt(2d0)
c        urvec(6,3,1,1) = 0d0
c        urvec(6,6,1,1)= 1d0/2d0
c        urvec(3,3,2,2)  = 1d0/2d0
c        urvec(3,6,2,2) = 1d0/sqrt(2d0)
c        urvec(6,3,2,2) = 0d0
c        urvec(6,6,2,2)= 1d0/2d0
c
c        uivec(3,3,1,1)  = 1d0/2d0
c        uivec(3,6,1,1) = 0d0
c        uivec(6,3,1,1) = -1d0/sqrt(2d0)
c        uivec(6,6,1,1)= 1d0/2d0
c        uivec(3,3,2,2)  = 1d0/2d0
c        uivec(3,6,2,2) = 0d0
c        uivec(6,3,2,2) = -1d0/sqrt(2d0)
c        uivec(6,6,2,2)= 1d0/2d0

c        do im = 1,2
c          do ii = 1, ntoob
c            do jj = 1, ntoob
c              urvec(ii,jj,im,im) = 0.1d0 * 
c     &             sqrt(abs((ii-1.5d0*jj)/(ii+1.5d0*jj)))
c            end do
c          end do
c        end do
c        do im = 1,2
c          do ii = 1, ntoob
c            do jj = 1, ntoob
c              uivec(ii,jj,im,im) = 0.1d0*(ii-2d0*jj)/(ii+2d0*jj)
c            end do
c          end do
c        end do

        call vec_to_disc(omvec,nlen,1,-1,luom)
        call vec_to_disc(urvec,nlen,1,-1,luur)
        call vec_to_disc(uivec,nlen,1,-1,luui)

        !
        imode = 11
        call cmbamp(imode,luom,luur,luui,luamp,
     &       omvec,nlen,nlen,nlen)

 200    continue

      end if

      ! Header for iteration info
      if (calc_Omg.and.calc_gradE) then
        write (6,'(">>>",2a/,">>>",2a)')
     &       '  iter              energy   variance     norm(G) ',
     &       '  norm(dE/dG) norm(Omega)',
     &       '--------------------------------------------------',
     &       '--------------------------'
      else if (calc_Omg) then
        write (6,'(">>>",2a/,">>>",2a)')
     &       '  iter              energy   variance     norm(G) ',
     &       '  norm(Omega)',
     &       '--------------------------------------------------',
     &       '--------------'
      else if (calc_gradE) then
        write (6,'(">>>",2a/,">>>",2a)')
     &       '  iter              energy   variance     norm(G) ',
     &       '  norm(dE/dG)',
     &       '--------------------------------------------------',
     &       '--------------'
      end if

      xngrad = 1000
      xnomg  = 1000
      itask = 0
      imacit = 0
      imicit = 0
      imicit_tot = 0
      energy = 0d0
      itask = 0
      nrdvec = 0
      did_rdvec = .false.
      do while (itask.lt.8)

        call atim(cpu0i,wall0i)

        call memchk2('b optc')

        if (igtbmod.ne.2) then
          ! usual route:
          nwfpar = n_cc_amp
          if (igtb_closed.eq.1) nwfpar = namp_packed
          call optcont(imacit,imicit,imicit_tot,iprint,
     &                   itask,iconv,
     &                   luamp,lutrvec,
     &                   energy,
     &                   ccvec1,ccvec2,nwfpar,
     &                   luomg,lusig,ludia,
     &                   nrdvec,lurdvec)
        else
          call optcont(imacit,imicit,imicit_tot,iprint,
     &                   itask,iconv,
     &                   luamp,lutrvec,
     &                   energy,
     &                   work(khvec1),work(khvec2),n_l_amp,
     &                   luomg,lusig,ludia,
     &                   0,lurdvec)
        end if
        call memchk2('a optc')

        if (igtbmod.lt.2) then
          ! the usual route:
          if (igtb_closed.eq.0) then
            call vec_from_disc(ccvec1,n_cc_amp,1,lblk,luamp)
            xnamp = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))
          else if (igtb_closed.eq.1) then
            ! expand to full spin-orbital basis, if necessary
            call vec_from_disc(ccvec2,namp_packed,1,lblk,luamp)
            xnamp = sqrt(inprod(ccvec2,ccvec2,namp_packed))
            iway = -1 ! unpack
            idual = 3
            call pack_g(iway,idual,isymmet_G,ccvec2,ccvec1,
     &                  n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                  n11amp,n33amp,iamp_packed,n_cc_amp)
          else
            write(6,*) 'igtb_closed has strange value'
            stop 'gtbce'
          end if


          if (isymmet_G.ne.0) then
            write(6,*) 'checking new T:'
c            call vec_from_disc(ccvec1,n_cc_amp,1,lblk,luamp)
            call chksym_t(isymmet_G,1,
     &           ccvec1,ccvec2,
     &           ictp,i_cc_typ,n_cc_typ,
     &           namp_cc_typ,ioff_cc_typ,ngas)
          end if

          
          if (ntest.ge.1000) then
            write(6,*) 'The new operator:'
            call wrt_cc_vec2(ccvec1,6,'GEN_CC')
          end if

        else if (igtbmod.eq.2) then
          call vec_from_disc(work(khvec1),n_cc_amp,1,lblk,luamp)
          xnamp = sqrt(inprod(work(khvec1),work(khvec1),n_l_amp))
          ! the '0' actually means, that we have so far identical
          ! alpha and beta parts for L 
          call l2g(work(khvec1),ccvec1,nspobex_tp,
     &         work(klsobex),work(klibsobex),
     &         0  ,ntoob)
        else if (igtbmod.eq.3) then
          imode=01
          call cmbamp(imode,luom,luur,luui,luamp,
     &       omvec,nlen,nlen,nlen)
          call vec_from_disc(omvec,nlen,1,-1,luom)
          call vec_from_disc(urvec,nlen,1,-1,luur)
          call vec_from_disc(uivec,nlen,1,-1,luui)

          call uou2g(omvec,urvec,uivec,ccvec1,
     &         nspobex_tp,
     &         work(klsobex),work(klibsobex),ntoob)
          
          write(6,*)
     &         '=============================================='
          write(6,*) 'calling chksym_t for the new variant:'

c          call chksym_t(isymmet_G,1,
c     &       ccvec1,ccvec2,
c     &       ictp,i_cc_typ,n_cc_typ,
c     &       namp_cc_typ,ioff_cc_typ,ngas)

          write(6,*)
     &         '=============================================='


        end if
        if (isymmet_G.ne.0) then
          call chksym_t(isymmet_G,1,
     &       ccvec1,ccvec2,
     &       ictp,i_cc_typ,n_cc_typ,
     &       namp_cc_typ,ioff_cc_typ,ngas)
        end if


        if (iand(itask,1).eq.1) then
* calculate energy ...
          call gtbce_E(igtbmod,elen,variance,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)
        end if

        if (iand(itask,2).eq.2) then
          if (calc_Omg) then
* ...  and vector function (Nakasuji CSE residual)  ...
            call gtbce_Omg(ccvec2,xnomg,
     &               elen,ovl,iopsym,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luec,luhc,lusc1,lusc2)
          end if

        ! we currently overwrite Omega if gradient is calculated
        ! ... I know, the usage of files to pass vectors would be
        ! more appropriate, but for the moment it is as it is
          if (calc_gradE) then
            inumint=1
            igrdmod=1
            npts=5
            call gtbce_gradE(
     &               isymmet_G,ccvec2,xngrad,igrdmod,
     &               inumint,npts,
     &               elen,ovl,
     &               ccvec1,iopsym,ccvec3,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_typ,i_cc_typ,ictp,
     &               namp_cc_typ,ioff_cc_typ,
     &               n_cc_amp,mxb_ci,nprint,
     &               luamp,luc,luec,luhc,
     &               lusc1,lusc2,lusc3,lusc4,lusc5,lusc6)
            if (igtbmod.eq.2) then
              ! transform into L-gradient
              call ggrad2lgrad(ccvec2,work(khvec2),work(khvec1),
     &             nspobex_tp,work(klsobex),0,ntoob)
              xngrad = sqrt(inprod(work(khvec2),work(khvec2),n_l_amp))
            else if (igtbmod.eq.3) then
              ! transform into Om-gradient
              call ggrad2omgrad(ccvec2,omgrd,omvec,urvec,uivec,
     &             nspobex_tp,work(klsobex),ntoob)
              ! transform into Ur-gradient
              irmod = 1
              call ggrad2ugrad(ccvec2,urgrd,omvec,urvec,uivec,
     &             nspobex_tp,work(klsobex),ntoob,irmod)
              ! transform into Ui-gradient
              irmod = 2
              call ggrad2ugrad(ccvec2,uigrd,omvec,uivec,urvec,
     &             nspobex_tp,work(klsobex),ntoob,irmod)

              xnom = sqrt(inprod(omgrd,omgrd,4*ntoob**2))
              xnur = sqrt(inprod(urgrd,urgrd,4*ntoob**2))
              xnui = sqrt(inprod(uigrd,uigrd,4*ntoob**2))
              
              write (6,'(">>>",i6," |grd|: ",3(2x,e10.4))')
     &             imacit,xnom,xnur,xnui
              xnom = sqrt(inprod(omvec,omvec,4*ntoob**2))
              xnur = sqrt(inprod(urvec,urvec,4*ntoob**2))
              xnui = sqrt(inprod(uivec,uivec,4*ntoob**2))
              
              write (6,'(">>>",i6," |vec|: ",3(2x,e10.4))')
     &             imacit,xnom,xnur,xnui

              if (mod(imacit,10).eq.0.or.
     &            imacit.eq.1.or.
     &            imacit.eq.maxmacit) then
                write(6,*) 'Information on vectors in iteration ',imacit
                write(6,*) 'Omega:'
                do ii = 1,2
                  do jj = 1,2
                    xnrm = sqrt(inprod(omvec(1,1,jj,ii),
     &                   omvec(1,1,jj,ii),ntoob**2))
                    write(6,*)'spin case ',ii,jj,xnrm
                    call wrtmat2(omvec(1,1,jj,ii),ntoob,ntoob,
     &                   ntoob,ntoob)
                  end do
                end do
                write(6,*) 'U(Re):'
                do ii = 1,2
                  xnrm = sqrt(inprod(urvec(1,1,ii,ii),
     &                 urvec(1,1,ii,ii),ntoob**2))
                  write(6,*)'spin case ',ii,ii,xnrm
                  call wrtmat2(urvec(1,1,ii,ii),ntoob,ntoob,
     &                 ntoob,ntoob)
                end do
                write(6,*) 'U(Im):'
                do ii = 1,2
                  xnrm = sqrt(inprod(uivec(1,1,ii,ii),
     &                 uivec(1,1,ii,ii),ntoob**2))
                  write(6,*)'spin case ',ii,ii,xnrm
                  call wrtmat2(uivec(1,1,ii,ii),ntoob,ntoob,
     &                 ntoob,ntoob)
                end do
                
                write(6,*) 'dE/dOmega:'
                do ii = 1,2
                  do jj = 1,2
                    xnrm = sqrt(inprod(omgrd(1,1,jj,ii),
     &                   omgrd(1,1,jj,ii),ntoob**2))
                    write(6,*)'spin case ',ii,jj,xnrm
                    call wrtmat2(omgrd(1,1,jj,ii),ntoob,ntoob,
     &                   ntoob,ntoob)
                  end do
                end do
                write(6,*) 'dE/dU(Re):'
                do ii = 1,2
                  xnrm = sqrt(inprod(urgrd(1,1,ii,ii),
     &                 urgrd(1,1,ii,ii),ntoob**2))
                  write(6,*)'spin case ',ii,ii,xnrm
                  call wrtmat2(urgrd(1,1,ii,ii),ntoob,ntoob,
     &                 ntoob,ntoob)
                end do
                write(6,*) 'dE/dU(Im):'
                do ii = 1,2
                  xnrm = sqrt(inprod(uigrd(1,1,ii,ii),
     &                 uigrd(1,1,ii,ii),ntoob**2))
                  write(6,*)'spin case ',ii,ii,xnrm
                  call wrtmat2(uigrd(1,1,ii,ii),ntoob,ntoob,
     &                 ntoob,ntoob)
                end do
                
              end if

            end if

          end if ! calc_gradE

          
          ! save gradient/omega
c        call vec_to_disc(ccvec1,n_cc_amp,1,lblk,luamp)
          if (igtbmod.lt.2) then
            ! the usual route:
            if (igtb_closed.eq.0) then
              call vec_to_disc(ccvec2,n_cc_amp,1,lblk,luomg)
            else
              iway = 2          ! pack and symmetrize
              idual = 3
              call pack_g(iway,idual,isymmet_G,ccvec1,ccvec2,
     &             n_cc_typ,i_cc_typ,ioff_cc_typ,
     &             n11amp,n33amp,iamp_packed,n_cc_amp)

              if (igtb_disptt.eq.1) then
                write(6,*) ' ACCORDING TO YOUR WISHES I DISPOSE THE '//
     &               'ANTISYMMETRIC PART OF dE/dG !!!'
                ccvec1(n11amp+1:n11amp+n33amp) = 0d0
              end if
              xngrad = sqrt(inprod(ccvec1,ccvec1,namp_packed))
              call vec_to_disc(ccvec1,namp_packed,1,lblk,luomg)
            end if
          else if (igtbmod.eq.2) then
            call vec_to_disc(work(khvec2),n_l_amp,1,lblk,luomg)
          else if (igtbmod.eq.3) then
            call vec_to_disc(omgrd,nlen,1,-1,luomgr)
            call vec_to_disc(urgrd,nlen,1,-1,luurgr)
            call vec_to_disc(uigrd,nlen,1,-1,luuigr)
            imode = 11
            call cmbamp(imode,luomgr,luurgr,luuigr,luomg,
     &           omvec,nlen,nlen,nlen)
            
          end if

c test and analysis routines follow:
          if (calc_gradE) then
            tstgrad = .false. !imacit.eq.3
            if (tstgrad.and.igtbmod.lt.2) then
              if (igtb_close.eq.0) then
                call copvec(ccvec2,ccvec3,n_cc_amp)
              else
                call copvec(ccvec1,ccvec3,n_cc_amp)
              end if

              ! vector is reloaded from luamp inside
              call gtbce_testgradE(igtbmod,
     &                       isymmet_G,igtb_closed,
     &                       ccvec3,ccvec2,xngrad,
     &                       ecore,
     &                       ccvec1,iopsym,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_typ,i_cc_typ,namp_cc_typ,ioff_cc_typ,
     &                       n_cc_amp,mxb_ci,
     &                       n11amp,n33amp,iamp_packed,ictp,
     &                       luamp,luomg,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)
              stop 'stop after testgradE'
            else if (tstgrad.and.igtbmod.eq.2) then
              call gtbce_testgradE_L(
     &                       work(khvec2),work(khvec1),
     &                       ecore,
     &                       ccvec1,iopsym,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_amp,n_l_amp,mxb_ci,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)              
              stop 'stop after testgradE_L'
            else if (tstgrad.and.igtbmod.eq.3.and.imacit.eq.5) then
              imode = 1
              write(6,*) 'calling test for Omega gradient'
              namp = 4*ntoob**2
              call gtbce_testgradE_UOU(imode,
     &                       omgrd,omvec,urvec,uivec,
     &                       elen,ecore,
     &                       ccvec1,iopsym,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_amp,namp,mxb_ci,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)              
              imode = 2
              write(6,*) 'calling test for U(R) gradient'
              namp = 50 !4*ntoob**2
              call gtbce_testgradE_UOU(imode,
     &                       urgrd,omvec,urvec,uivec,
     &                       elen,ecore,
     &                       ccvec1,iopsym,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_amp,namp,mxb_ci,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)              
              imode = 3
              write(6,*) 'calling test for U(I) gradient'
              namp = 50 !4*ntoob**2
              call gtbce_testgradE_UOU(imode,
     &                       uigrd,omvec,urvec,uivec,
     &                       elen,ecore,
     &                       ccvec1,iopsym,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_amp,namp,mxb_ci,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)              
              stop 'stop after testgradE_L'
            end if
          end if                ! calc_gradE (analysis mode)

        end if ! iand(itask,2)

        if (iand(itask,4).eq.4) then
          imode=1
          iomg =1
          inumint=1
          npnts = 5
          call gtbce_num2drv(igtbmod,imode,iomg,
     &                       igtb_closed,isymmet_G,
     &                       inumint,npnts,
     &                       ecore,
     &                       iccvec,nSdim,
     &                       ccvec1,iopsym,ccvec2,ccvec3,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_typ,i_cc_typ,ictp,
     &                       namp_cc_typ,ioff_cc_typ,
     &                       n_cc_amp,mxb_ci,
     &                       n11amp,n33amp,iamp_packed,
     &                       lusig,
     &                       luamp,lutrvec,luc,luec,luhc,
     &                       lusc1,lusc2,lusc3,lusc4,lusc5,lusc6,lusc7)

        end if

        do_rdvec = .false.
        if (igtb_prjout.eq.1.and.
     &         xnamp.gt.1d-6.and..not.did_rdvec
     &         .and.imicit.eq.0) do_rdvec=.true.

        if (do_rdvec) then
          did_rdvec = .true.
          if (igtbmod.ne.0) stop 'does not work'
          inumint=1
          npnts = 5
          comm_ops = .false.
c test
          irestart = 1
          if (irestart.ne.0) then
            iramp = irestart
            call mk_iccvec(isymmet_G,lufoo,iramp,
     &                    iccvec,nSdim,ccvec1,ccvec2,
     &                    n_cc_typ,i_cc_typ,ictp,
     &                    namp_cc_typ,ioff_cc_typ,ngas,
     &                    n_cc_amp)
          end if
          imode = 0
          call gtbce_h0(imode,igtb_closed,isymmet_G,
     &                  iccvec,nSdim,
     &                  ccvec1,ccvec2,ccvec3,
     &                  civec1,civec2,c2vec,
     &                  n_cc_amp,mxb_ci,
     &                  n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                  n11amp,n33amp,iamp_packed,
     &                  lufoo,ludum,
     &                  luamp,luec,luhc,
     &                  lusc1,lusc2)
c          if (iramp.lt.nsdim) then
c
c            call gtbce_foo( isymmet_G,iramp,
c     &                    inumint,npnts,
c     &                    ovl,
c     &                    iccvec,nSdim,
c     &                    ccvec1,iopsym,comm_ops,
c     &                    ccvec2,ccvec3,
c     &                    civec1,civec2,c2vec,
c     &                    n_cc_typ,i_cc_typ,ictp,
c     &                    namp_cc_typ,ioff_cc_typ,
c     &                    n_cc_amp,mxb_ci,
c     &                    lufoo,
c     &                    luamp,luc,luec,luhc,
c     &                    lusc1,lusc2,lusc3,lusc4,
c     &                    lusc5,lusc6,lusc7,lusc8,
c     &                    lusc9,lusc10)
c          end if
          call memman(idum,idum,'MARK  ',2,'FOO MA')
          lenhss=nSdim*nSdim
          call memman(khss,lenhss,'ADDL  ',2,'HSSIAN')
          istmode = 2
          call gtbce_getrdvec(isymmet_G,work(khss),lufoo,lurdvec,nrdvec,
     &                nSdim,n_cc_amp,iccvec,
     &                ccvec1,ccvec2)
          idum = 0
          call memman(idum,idum,'FLUSM ',2,'FOO MA')
        end if

        if (nrdvec.gt.0.and.iand(itask,2).eq.2) then
          call gtbce_prjout_rdvec(nrdvec,lurdvec,luomg,
     &         n_cc_amp,ccvec1,ccvec2)
          xngrad = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))
          if (isymmet_G.ne.0) then
            write(6,*) 'checking projected gradient:'
            call chksym_t(isymmet_G,1,
     &           ccvec1,ccvec2,
     &           ictp,i_cc_typ,n_cc_typ,
     &           namp_cc_typ,ioff_cc_typ,ngas)
          end if

        end if

* minimal output
        energy = elen + ecore
        if (imicit.eq.0.and..not.iand(itask,8).eq.8) then
         if (calc_Omg.and.calc_gradE) then
          write (6,'(">>>",i6,f21.12,4(2x,e10.4))')
     &             imacit,energy,variance,xnamp,xngrad,xnomg
         else if (calc_Omg) then
          write (6,'(">>>",i6,f21.12,3(2x,e10.4))')
     &             imacit,energy,variance,xnamp,xnomg
         else if (calc_gradE) then
          write (6,'(">>>",i6,f21.12,3(2x,e10.4))')
     &             imacit,energy,variance,xnamp,xngrad
         end if
         call util_flush(6)
        end if

* analysis section:
        do_eag = .false.
        do_foo = .false.
        do_hss = .false.

        if (imicit.eq.0) then
          do ii = 1, n_eag
            if (it_eag(ii).eq.imacit) do_eag = .true.
          end do

          do ii = 1, n_foo
            if (it_foo(ii).eq.imacit) do_foo = .true.
          end do

          do ii = 1, n_hss
            if (it_hss(ii).eq.imacit) do_hss = .true.
          end do
        end if

c        tst_hss = .false.
        if (do_eag) then
          if (igtbmod.ne.0) stop 'does not work'
          do ii = 1, nn_eag

c     reload amplitudes:
            if (igtb_closed.eq.0) then
              call vec_from_disc(ccvec1,n_cc_amp,1,-1,luamp)
            else
              call vec_from_disc(ccvec3,namp_packed,1,-1,luamp)
              iway = -1
              idual = 0
              call pack_g(iway,idual,isymmet_G,ccvec3,ccvec1,
     &             n_cc_typ,i_cc_typ,ioff_cc_typ,
     &             n11amp,n33amp,iamp_packed,n_cc_amp)
            end if

            write(6,'("@p",a,i4)') 'printout for amplitude ', ng_eag(ii)
            if (igtb_closed.eq.0) then
              ccvec2(1:n_cc_amp) = 0d0
              ccvec2(ng_eag(ii)) = 1d0
              if (isymmet_G.ne.0) then
                stop 'adapt this section'
              end if
            else
              ccvec3(1:n11amp+n33amp) = 0d0
              if (ng_eag(ii).ge.-1) then
                if (ng_eag(ii).ge.1) then
                  ccvec3(ng_eag(ii)) = 1d0
                else if(ng_eag(ii).eq.-1) then
                  ccvec3(1:n11amp+n33amp) = 1d0
                end if
                iway = -1
                idual = 0
                call pack_g(iway,idual,isymmet_G,ccvec3,ccvec2,
     &             n_cc_typ,i_cc_typ,ioff_cc_typ,
     &             n11amp,n33amp,iamp_packed,n_cc_amp)
              else if (ng_eag(ii).eq.-2) then
                ccvec2(1:n_cc_amp) = ccvec1(1:n_cc_amp)
              else if (ng_eag(ii).eq.-3) then
                stop 'not impl.'
c no no no
c                iramp = 0
c                call mk_iccvec(isymmet_G,lufoo,iramp,
c     &                    iccvec,nSdim,ccvec1,ccvec2,
c     &                    n_cc_typ,i_cc_typ,ictp,
c     &                    namp_cc_typ,ioff_cc_typ,ngas,
c     &                    n_cc_amp)
c                do iamp = 1, n_cc_amp
c                  if (iccvec(iamp).lt.1) then
c                    ccvec2(iamp) = -1d0
c                  else
c                    ccvec2(iamp) = 1d0
c                  end if
c                end do

              end if
            end if
            from_g = st_eag(ii)
            to_g   = en_eag(ii)
            npnts = np_eag(ii)

            call gtbce_EalongG(ccvec2,npnts,from_g,to_g,
     &               ecore,
     &               ccvec1,iopsym,ccvec3,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)
          end do
        end if

        if (do_foo) then
          if (igtbmod.ne.0) stop 'does not work'
          inumint=1
          npnts = 5
          comm_ops = .false.
c          call gtbce_foo_old(inumint,npnts,
c     &                   ovl,
c     &                   ccvec1,iopsym,comm_ops,
c     &                   ccvec2,ccvec3,
c     &                   civec1,civec2,c2vec,
c     &                   n_cc_amp,mxb_ci,
c     &                   lufoo,
c     &                   luamp,luc,luec,luhc,
c     &                   lusc1,lusc2,lusc3,lusc4,
c     &                   lusc5,lusc6,lusc7,lusc8)
c
c          stop 'test foo'
          call gtbce_foo( isymmet_G,0,
     &                    inumint,npnts,
     &                    ovl,
     &                    iccvec,nSdim,
     &                    ccvec1,iopsym,comm_ops,
     &                    ccvec2,ccvec3,
     &                    civec1,civec2,c2vec,
     &                    n_cc_typ,i_cc_typ,ictp,
     &                    namp_cc_typ,ioff_cc_typ,
     &                    n_cc_amp,mxb_ci,
     &                    lufoo,
     &                    luamp,luc,luec,luhc,
     &                    lusc1,lusc2,lusc3,lusc4,
     &                    lusc5,lusc6,lusc7,lusc8,
     &                    lusc9,lusc10)
          call memman(idum,idum,'MARK  ',2,'FOO MA')
          lenhss=nSdim*nSdim
          call memman(khss,lenhss,'ADDL  ',2,'HSSIAN')
          istmode = 2
          call gtbce_anahss(work(khss),lufoo,ludum,istmode,
     &                nSdim,n_cc_typ,i_cc_typ,
     &                namp_cc_typ,ioff_cc_typ,iopsym)
          idum = 0
          call memman(idum,idum,'FLUSM ',2,'FOO MA')
        end if

        do_h0 = 
     &       i_do_h0.ne.0.and.(xnamp.gt.1d-6)
     &       .and.(imacit.eq.2.or.mod(imacit,30).eq.0)
     &       .and.imicit.eq.0

        if (do_h0) then
          if (isymmet_G.ne.0.and.igtb_closed.eq.0) then
            iramp = 0
            call mk_iccvec(isymmet_G,lufoo,iramp,
     &                    iccvec,nSdim,ccvec1,ccvec2,
     &                    n_cc_typ,i_cc_typ,ictp,
     &                    namp_cc_typ,ioff_cc_typ,ngas,
     &                    n_cc_amp)
          else if (igtb_closed.eq.1) then
            nSdim = namp_packed
          else
            nSdim = n_cc_amp
          end if

          imode = 2
          call gtbce_h0(imode,igtb_closed,isymmet_G,
     &                  iccvec,nSdim,
     &                  ccvec1,ccvec2,ccvec3,
     &                  civec1,civec2,c2vec,
     &                  n_cc_amp,mxb_ci,
     &                  n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                  n11amp,n33amp,iamp_packed,
     &                  luh0,ludia,
     &                  luamp,luec,luhc,
     &                  lusc1,lusc2)

c          idum = 0
c          call memman(idum,idum,'MARK  ',2,'HESSMA')
c          lenhss=nSdim*nSdim
c          call memman(khss,lenhss,'ADDL  ',2,'HSSIAN')
c          istmode = 3
c          call gtbce_anahss(work(khss),luh0,ludia,istmode,
c     &                nSdim,n_cc_typ,i_cc_typ,
c     &                namp_cc_typ,ioff_cc_typ,iopsym)
c
c          idum = 0
c          call memman(idum,idum,'FLUSM ',2,'HESSMA')
c
        end if

c          call rewino(lufoo)
c          call rewino(luhss)
c          do ii = 1, n_cc_amp
c            print *,'column ',ii
c            call cmp2vcd(ccvec2,ccvec3,lufoo,luhss,1d-10,0,lblk)
c          end do
        if (do_hss) then
          if (igtbmod.ne.0.and.igtbmod.ne.2) stop 'does not work'
          if (isymmet_G.ne.0.and.igtb_closed.eq.0) then
            iramp = 0
            call mk_iccvec(isymmet_G,lufoo,iramp,
     &                    iccvec,nSdim,ccvec1,ccvec2,
     &                    n_cc_typ,i_cc_typ,ictp,
     &                    namp_cc_typ,ioff_cc_typ,ngas,
     &                    n_cc_amp)
          else if (igtb_closed.eq.1) then
            nSdim = n11amp+n33amp
          else
            nSdim = n_cc_amp
          end if

c test h0
c          call gtbce_h0(isymmet_G,
c     &                  iccvec,nSdim,
c     &                  ccvec1,ccvec2,
c     &                  civec1,civec2,c2vec,
c     &                  n_cc_amp,mxb_ci,
c     &                  luh0,
c     &                  luamp,luec,luhc,
c     &                  lusc1,lusc2)
c          idum = 0
c          call memman(idum,idum,'MARK  ',2,'HESSMA')
c          lenhss=nSdim*nSdim
c          call memman(khss,lenhss,'ADDL  ',2,'HSSIAN')
c          istmode = 3
c          call gtbce_anahss(work(khss),luh0,ludum,istmode,
c     &                nSdim,n_cc_typ,i_cc_typ,
c     &                namp_cc_typ,ioff_cc_typ,iopsym)
c
c          idum = 0
c          call memman(idum,idum,'FLUSM ',2,'HESSMA')

          imode=2
          iomg =1
          inumint=1
          npnts = 5
          call gtbce_num2drv(igtbmod,imode,iomg,
     &                       igtb_closed,isymmet_G,
     &                       inumint,npnts,
     &                       ecore,
     &                       iccvec,nSdim,
     &                       ccvec1,iopsym,ccvec2,ccvec3,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_typ,i_cc_typ,ictp,
     &                       namp_cc_typ,ioff_cc_typ,
     &                       n_cc_amp,mxb_ci,
     &                       n11amp,n33amp,iamp_packed,
     &                       luhss,
     &                       luamp,luleq,luc,luec,luhc,
     &                       lusc1,lusc2,lusc3,lusc4,lusc5,lusc6,lusc7)
          
          idum = 0
          call memman(idum,idum,'MARK  ',2,'HESSMA')
          lenhss=nSdim*nSdim
          call memman(khss,lenhss,'ADDL  ',2,'HSSIAN')
          istmode = 1
          call gtbce_anahss(work(khss),luhss,ludum,istmode,
     &                nSdim,n_cc_typ,i_cc_typ,
     &                namp_cc_typ,ioff_cc_typ,iopsym)

          idum = 0
          call memman(idum,idum,'FLUSM ',2,'HESSMA')


        end if

        call memchk2('afcalc')

        call atim(cpui,walli)
        call prtim(6,'time for current iteration',
     &       cpui-cpu0i,walli-wall0i)

      end do ! optimization loop
 
      call atim(cpu,wall)
      call prtim(6,'time in GTBCE optimization',
     &       cpu-cpu0,wall-wall0)


      ! somewhat unmotivated here, actually just for looking at
      ! the amplitudes in another way:
      if (igtbmod.eq.1) then
        call can2str(2,work(kcan),ccvec1,
     &       nspobex_tp,i_cc_typ,ioff_cc_typ)
      end if

      write (6,*) ' ANALYSIS: '
      if (igtb_closed.eq.0) then
        call vec_from_disc(ccvec1,n_cc_amp,1,-1,luamp)
        call ana_gencc(ccvec1,1)
      else
        write(6,*) ' ANALYSIS in spin-adapted basis: '
        call vec_from_disc(ccvec2,namp_packed,1,-1,luamp)
        call ana_gucc(ccvec2,n11amp,n33amp,iamp_packed,
     &                ireost,nsmob,ntoob)
        iway = -1
        idual = 3
        call pack_g(iway,idual,isymmet_G,ccvec2,ccvec1,
     &             n_cc_typ,i_cc_typ,ioff_cc_typ,
     &             n11amp,n33amp,iamp_packed,n_cc_amp)
        write(6,*) ' ANALYSIS in spin-orbital basis: '
        call ana_gencc(ccvec1,1)
      end if

      idum = 0
      call memman(idum,idum,'FLUSH  ',idum,'GTBCOP')

      return
      end
**********************************************************************
**********************************************************************
* DECK: gtbce_initG
**********************************************************************
      subroutine gtbce_initG(ccamp,
     &                       imode,luamp,
     &                       ccscr,
     &                       ngas_,iocc,ihpv,n_cc_amp,i_cc_typ,n_cc_typ,
     &                       namp_cc_typ,ioff_cc_typ)
**********************************************************************
* 
* purpose: initialize G (depending on imode) with
*     
*          -1 :    automatic
*           0 :    zero
*           1 :    a full previous G vector on luamp
*           2 :    a singles and doubles vector on luamp
*           3 :    a doubles vectors on luamp
*
*  ak, early 2004
*
**********************************************************************
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cc_exc.inc'
      include 'orbinp.inc'
      include 'cgas.inc'
      include 'csm.inc'
* input
      integer, intent(in) ::
     &     ihpv(ngas), iocc(mxpngas,2), i_cc_typ(ngas,4,n_cc_typ),
     &     ioff_cc_typ(n_cc_typ), namp_cc_typ(n_cc_typ)     
* output
      real*8, intent(out) ::
     &     ccamp(n_cc_amp), ccscr(n_cc_amp)
* constants
      integer, parameter ::
     &     ntest = 00
* local scratch
      logical ::
     &     dont, not_possible
      integer ::
     &     ioff_sing(2), ilen_sing(2), ioff_doub(3), ilen_doub(3),
     &     nph(nsmst,2)
      character*8 cctype

      if (ntest.ge.5) then
        write(6,*) '==========='
        write(6,*) 'gtbce_initG'
        write(6,*) '==========='
        write(6,*) ' imode = ',imode
        write(6,*) ' luamp = ',luamp
        write(6,*) ' mscomb_cc = ',mscomb_cc
      end if

      nsing = 2
      ndoub = 3
      if (mscomb_cc.ne.0) then
        nsing = 1  ! only alpha part
        ndoub = 2  ! only alpha and alpha/beta part
      end if

      imode_ = imode
      ! test for existence of file
      if (imode.gt.0.or.imode.eq.-1) then
        rewind(luamp,err=100)
        read(luamp,err=100,end=100) namp_read
        if (namp_read.gt.0.and.namp_read.le.n_cc_amp) goto 200

 100    write(6,*) 'no proper amplitudes found to restart from'
        imode_ = 0

 200    continue
        
      end if


      if (imode_.eq.2.or.imode_.eq.3.or.imode_.eq.-1) then
        ! get the D or SD vector
        lblk = -1
        call vec_from_disc(ccscr,n_cc_amp,1,lblk,luamp)
        ! find the matching blocks in G
        ! and hope that LUCIA keeps the ordering of the blocks
        ioff_sing(1:2) = 0  !(alpha / beta)
        ioff_doub(1:3) = 0  !(alpha-beta / alpha-alpha / beta-beta) 
        not_possible = .true.
        do itp = 1, n_cc_typ
          nca = 0
          ncb = 0
          naa = 0
          nab = 0
          dont = .false.
          do igs = 1, ngas
            if (ihpv(igs).eq.1) then ! hole space
              naa = naa + i_cc_typ(igs,3,itp) 
              nab = nab + i_cc_typ(igs,4,itp) 
              if (i_cc_typ(igs,1,itp).gt.0.or.
     &            i_cc_typ(igs,2,itp).gt.0    ) then
                dont = .true.
              end if
            else if(ihpv(igs).eq.2) then ! particle space
              nca = nca + i_cc_typ(igs,1,itp) 
              ncb = ncb + i_cc_typ(igs,2,itp) 
              if (i_cc_typ(igs,3,itp).gt.0.or.
     &            i_cc_typ(igs,4,itp).gt.0    ) then
                dont = .true.
              end if
            else if(ihpv(igs).eq.3) then ! valence space
              not_possible = .true.   ! we cannot handle this currently
              stop 'valence spaces are too difficult for me!'
            else
              stop'ihpv is inconsistent in init_gtbce'
            end if
          end do

          if (ntest.ge.100) then
            write(6,*) 'ityp = ',itp
            write(6,*) ' nca, ncb ', nca, ncb
            write(6,*) ' naa, nab ', naa, nab
            write(6,*) ' dont     ',dont
          end if

          if (.not.dont) then
            if (nca.eq.1.and.ncb.eq.0.and.
     &          naa.eq.1.and.nab.eq.0     ) then
              ioff_sing(1) = ioff_cc_typ(itp)
              ilen_sing(1) = namp_cc_typ(itp)
            else if (nca.eq.0.and.ncb.eq.1.and.
     &          naa.eq.0.and.nab.eq.1     ) then
              ioff_sing(2) = ioff_cc_typ(itp)
              ilen_sing(2) = namp_cc_typ(itp)
            else if (nca.eq.1.and.ncb.eq.1.and.
     &          naa.eq.1.and.nab.eq.1     ) then
              ioff_doub(1) = ioff_cc_typ(itp)
              ilen_doub(1) = namp_cc_typ(itp)
            else if (nca.eq.2.and.ncb.eq.0.and.
     &          naa.eq.2.and.nab.eq.0     ) then
              ioff_doub(2) = ioff_cc_typ(itp)
              ilen_doub(2) = namp_cc_typ(itp)
            else if (nca.eq.0.and.ncb.eq.2.and.
     &          naa.eq.0.and.nab.eq.2     ) then
              ioff_doub(3) = ioff_cc_typ(itp)
              ilen_doub(3) = namp_cc_typ(itp)
            end if
          end if
        end do
        
        if (mscomb_cc.ne.0) then
          ! don't worry about missing info
          ioff_sing(2) = 1
          ilen_sing(2) = 0
          ioff_doub(3) = 1
          ilen_doub(3) = 0
        end if

        if (ntest.ge.5) then
          write(6,*) 'offsets and lengthes extracted:'
          write(6,*) '(mscomb_cc = ',mscomb_cc,')'
          write(6,*) ioff_sing(1:nsing), ioff_doub(1:ndoub)
          write(6,*) ilen_sing(1:nsing), ilen_doub(1:ndoub)
        end if

        if (ilen_sing(1)*ilen_sing(nsing).eq.0) then
          ! try to guess singles size from the number of 
          ! possible holes and particles
          nph(1:nsmst,1:2) = 0
          do igs = 1, ngas
            if (igs.eq.1) then
              nelmin = iocc(1,1)
              nelmax = iocc(1,2)
            else
              nelmin = iocc(igs,1)-iocc(igs-1,2)
              nelmax = iocc(igs,2)-iocc(igs-1,1)
            end if
            ! may at least one electron be removed in this space?
            ihp = 0
            if (nelmin.lt.2*nobpt(igs).and.ihpv(igs).eq.1) ihp=1
            ! may at least one electron be added in this space
            if (nelmax.gt.0.and.ihpv(igs).eq.2) ihp=2
            if (ihp.gt.0) then
              do ism = 1, nsmst
                ! get the number of holes/particles per symmetry
                nph(ism,ihp) = nph(ism,ihp) + ngssh(ism,igs)
              end do
            end if
          end do
          lsing = 0
          do ism = 1, nsmst
            lsing = lsing + nph(ism,1)*nph(ism,2)
          end do
          ! there has to be done some more work for open-shell cases!
          ! for now:
          ilen_sing(1:nsing) = lsing 

          write(6,*) 'There seem to be no singles in your general '//
     &               'TWOBODY operator!'
          write(6,*) 'From the number of active holes and particles'//
     &               ' I guess ',ilen_sing(1:nsing)
          
        end if
        
        if (ioff_doub(1)*ioff_doub(2)*ioff_doub(3).eq.0) then
          write(6,*) 'No offsets for doubles found!!!'
          stop 'difficulties in gtbce_init'
        end if

        ! decide what to do
        if (imode.eq.-1) then
          namp_d  = ilen_doub(1)+ilen_doub(2)+ilen_doub(3)
          namp_sd = namp_d + ilen_sing(1) + ilen_sing(2)
          imode_ = 0
          if (namp_read.eq.namp_d ) imode_ = 3
          if (namp_read.eq.namp_sd) imode_ = 2
          if (namp_read.eq.n_cc_amp)imode_ = 1
          
          if (ntest.ge.5) then
            write (6,*) 'namp_read ',namp_read
            write (6,*) 'namp_d    ',namp_d
            write (6,*) 'namp_sd   ',namp_sd
            write (6,*) 'n_cc_amp  ',n_cc_amp
            write (6,*) ' imode_  =',imode_
          end if
        end if ! imode.eq.-1
  
      end if ! imode_.eq.2/3/-1

      if (imode_.eq.0) then
        ccamp(1:n_cc_amp) = 0d0
      else if (imode_.eq.1) then
        lblk = -1
        call vec_from_disc(ccamp,n_cc_amp,1,lblk,luamp)
      else if (imode_.eq.2.or.imode_.eq.3) then
        ioff1 = 0
        if (imode_.eq.2) then
          do ii = 1, nsing
            if (ilen_sing(ii).gt.0.and.ioff_sing(ii).gt.0)
     &         ccamp(ioff_sing(ii)  :ioff_sing(ii)+ilen_sing(ii)-1) =
     &         ccscr(ioff1        +1:ioff1        +ilen_sing(ii))
            ioff1 = ioff1 + ilen_sing(ii)
          end do
        end if
        do ii = 1, ndoub
          if (ilen_doub(ii).gt.0.and.ioff_doub(ii).gt.0)
     &         ccamp(ioff_doub(ii)  :ioff_doub(ii)+ilen_doub(ii)-1) =
     &         ccscr(ioff1        +1:ioff1        +ilen_doub(ii))
          ioff1 = ioff1 + ilen_doub(ii)
        end do

      else
        write(6,*) 'unknown imode in init_gtbce(', imode,') !'
        stop 'init_gtbce'
      end if

      if (ntest.ge.100) then
        write(6,*) 'Initialized G: '
        call wrt_cc_vec2(ccamp,6,'GEN_CC')
      end if

      return
      end
**********************************************************************
**********************************************************************
* DECK: gtbce_E
**********************************************************************
      subroutine gtbce_E(igtbmod_l,
     &                   elen,variance,ovl,
     &                   e_core,
     &                   ccvec1,iopsym,ccvecscr,
     &                   civec1,civec2,c2vec,
     &                   n_cc_amp,mxb_ci,
     &                   luc,luec,luhc,lusc1,lusc2)
**********************************************************************
*
* purpose: calculate the Energy of the GTBCE.
*
*    E = <0|exp(G^+) H exp(G)|0> / <0|exp(G^+)exp(G^+)|0>
*
*  input:          |0>    on luc
*                   G     on ccvec1 
*
*  output:   exp(G)|0>    on luec
*          H exp(G)|0>    on luhc
*
*          E                            on elen
*          S = <0|exp(G^+)exp(G^+)|0>   on ovl
*          v = <H^2>/S - E^2            on variance
*
*  igtbmod_l.eq.(0/2) proceed as usual
*  igtbmod_l.eq.1     use exp(G^2)
*
*  ak, early 2004
*
**********************************************************************
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
c      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cands.inc'
      include 'oper.inc'
      include 'gtbce.inc'
* debugging:
      integer, parameter :: ntest = 5

* input arrays
      real*8 ccvec1(n_cc_amp)

* local
      logical test_h1

* scratch arrays
      character*8 cctype
      real*8 civec1(mxb_ci),civec2(mxb_ci),c2vec(*),
     &       ccvecscr(n_cc_amp)
* external functions
      real*8 inprod, inprdd

      call atim(cpu0,wall0)

      ! settings for expt_ref2
      thresh=expg_thrsh
      mx_term=-mxterm_expg
      cctype='GEN_CC'

      if (ntest.ge.5) then
        write (6,*) '================='
        write (6,*) ' This is gtbce_E '
        write (6,*) '================='
        write (6,*) 
        write (6,*) 'on entry: '
        write (6,*) 'e_core   : ', e_core
        write (6,*) 'n_cc_amp,mxb_ci : ', n_cc_amp,mxb_ci
        write (6,*) 'luc,luec,luhc,lusc1,lusc2: ',
     &               luc,luec,luhc,lusc1,lusc2
        write (6,*) 'igtbmod_l: ',igtbmod_l
      end if
      if (ntest.ge.5) then
        write(6,*) ' gtbce_E > '
        xnorm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))
        write(6,*) '     n_cc_amp,norm of T: ',n_cc_amp,xnorm
      end if
      if (ntest.ge.100) then
        call wrt_cc_vec2(ccvec1,6,cctype)
      end if
      
      lblk = -1
*--------------------------------------------------------------------*
* |0tilde> = exp(G)|0>
*
*  |0> on luc, |0tilde> on luec, 
*  G is on ccvec1
*--------------------------------------------------------------------*
      if (igtbmod_l.ne.1) then
        call expt_ref2(luc,luec,luhc,lusc1,lusc2,
     &              thresh,mx_term, ccvec1, ccvecscr, civec1, civec2,
     &              n_cc_amp,cctype,iopsym)
      else
        call expt2_ref(luc,luec,luhc,lusc1,lusc2,
     &              thresh,mx_term,
     &              1d0,ccvec1, ccvecscr, civec1, civec2, n_cc_amp,
     &              iopsym)
      end if
*--------------------------------------------------------------------*
* |H0tilde> = H exp(G)|0>
*
*  |H0tilde> on luhc
*--------------------------------------------------------------------*
      if (igtb_test_h1.eq.1) i12 = 1
      call mv7(civec1,civec2,luec,luhc)
*--------------------------------------------------------------------*
* S = <0tilde|0tilde>
*--------------------------------------------------------------------*
      xs = inprdd(civec1,civec2,luec,luec,1,lblk)
      if (xs.eq.0) then
        write(6,*) 'gtbce_E > Wavefunction with zero norm!!'
        write(6,*) '          Are we trying to be funny today?'
        stop 'fatal inconsistency'
      end if
*--------------------------------------------------------------------*
* E S = <0tilde|H|0tilde>, E = <0tilde|H|0tilde>/S
*--------------------------------------------------------------------*
      xes= inprdd(civec1,civec2,luec,luhc, 1,lblk)
      elen = xes/xs
      ovl = xs
*--------------------------------------------------------------------*
* variance of <H>: <0tilde|H^2|0tilde>/S - E^2
*--------------------------------------------------------------------*
      xh2 = inprdd(civec1,civec2,luhc,luhc,1,lblk)
      variance = xh2/xs - xes*xes/(xs*xs)
      if (ntest.ge.5) then
        write(6,*) ' gtbce_E > '
        write(6,*) '       <0tilde|0tilde> = ',xs
        write(6,*) '     <0tilde|H|0tilde> = ',xes
        write(6,*) '   <0tilde|H^2|0tilde> = ',xh2
        write(6,*) '           el. energy  = ',elen
        write(6,*) '               e_core  = ',e_core
        write(6,*) '               energy  = ',elen+e_core
        write(6,*) '             variance  = ',variance
      end if
      if (ntest.ge.1000) then
        write(6,*) ' gtbce_E > '
        write(6,*) ' |0tilde>:'
        call wrtvcd(civec1,luec,1,lblk)
        write(6,*) ' H|0tilde>:'
        call wrtvcd(civec1,luhc,1,lblk)
      end if

      call atim(cpu,wall)
      call prtim(6,'time in gtbce_E',cpu-cpu0,wall-wall0)

      return
      end
*--------------------------------------------------------------------*
**********************************************************************
* DECK: gtbce_Omg
**********************************************************************
      subroutine gtbce_Omg(omg,xnomg,
     &                     elen,ovl,iopsym,
     &                     civec1,civec2,c2vec,
     &                     n_cc_amp,mxb_ci,
     &                     luec,luhc,lusc1,lusc2)
**********************************************************************
*
* purpose: calculate the Nakasuji-type
*          Vectorfunction Omega of the GTBCE (or Contracted Schroedinger
*          Equations (CSE) residual, if you will, if the operator space
*          was chosen accordingly (SING,0,0,0/DOUB,1,1,1,1,1))
*
*    Omg = 1/S <0|exp(G^+) gamma (H-E) exp(G)|0>
*
*   input:   exp(G)|0> on luec
*          H exp(G)|0> on luhc
*
*   output:   Omg      on omg
*            |Omg|     on xnomg
*
*  ak, early 2004
*
**********************************************************************
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
c      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cands.inc'
* debugging:
      integer, parameter :: ntest = 0

* input/output arrays
      real*8 omg(n_cc_amp)
* scratch arrays
      character*8 cctype
      real*8 civec1(mxb_ci),civec2(mxb_ci),c2vec(*)
* external functions
      real*8 inprod, inprdd

      call atim(cpu0,wall0)

      if (ntest.ge.5) then
        write (6,*) '========================='
        write (6,*) ' This is gtbce_Omg'
        write (6,*) '========================='
        write (6,*) 
        write (6,*) 'on entry: '
        write (6,*) 'el. energy,mxb_ci : ', elen, mxb_ci
        write (6,*) 'luec,luhc,lusc1,lusc2: ',
     &               luec,luhc,lusc1,lusc2
      end if
      
      lblk = -1
*--------------------------------------------------------------------*
* (H-E)|0tilde>
*  result on lusc1
*--------------------------------------------------------------------*
      call vecsmdp(civec1,civec2,1d0,-elen,luhc,luec,lusc1,1,lblk)
*--------------------------------------------------------------------*
* Omg_u = <0(tilde)|gamma_u(H-E)|0(tilde)>
*  result on omg
*--------------------------------------------------------------------*
      isigden=2
      omg(1:n_cc_amp) = 0d0
      call sigden_cc(civec1,civec2,luec,lusc1,omg,isigden)
      if (iopsym.eq.1.or.iopsym.eq.-1) then
        if (iopsym.eq.-1) call scalve(omg,-1d0,n_cc_amp)
        call conj_t
        call sigden_cc(civec1,civec2,luec,lusc1,omg,isigden)
        call conj_t
        if (iopsym.eq.-1) call scalve(omg,-1d0,n_cc_amp)
      end if
c      call memchk
      call scalve(omg,1d0/ovl,n_cc_amp)

      xnomg = sqrt(inprod(omg,omg,n_cc_amp))

      if (ntest.ge.5) then
c        call memchk
        write(6,*) ' gtbce_Omg > '
        write(6,*) '     n_cc_amp,norm of omega: ',n_cc_amp,xnomg
      end if
      if (ntest.ge.100) then
        cctype='GEN_CC'
        call wrt_cc_vec2(omg,6,cctype)
      end if
      
      call atim(cpu,wall)
      call prtim(6,'time in gtbce_Omg',cpu-cpu0,wall-wall0)

      return
      end
*--------------------------------------------------------------------*
**********************************************************************
* DECK: gtbce_gradE
**********************************************************************
      subroutine gtbce_gradE(!igtbmod,
     &                       isymmet_G,grad,xngrad,igradmode,
     &                       imode,npnts,
     &                       elen,ovl,
     &                       ccvec1,iopsym,ccvec2,ccvec3,
     &                       civec1,civec2,c2vec,
     &                       n_cc_typ,i_cc_typ,ictp,
     &                       namp_cc_typ,ioff_cc_typ,
     &                       n_cc_amp,mxb_ci,nprint,
     &                       luamp,luc,luec,luhc,
     &                       lusc1,lusc2,lusc3,lusc4,lusc5,lusc6)
**********************************************************************
*
* purpose: calculate the gradient of the GTBCE energy by numerical
*          integration of the Wilcox identity
*
*          Ref. van Voorhis, Head-Gordon, JCP 115(11) 5033 (2001)
*
*    gradE = 
*      2/S int_0^1 da <0|exp(G^+) (H-E) exp((1-a)G) gamma exp(aG)|0>
*
*   input:         |0> on luc
*            exp(G)|0> on luec
*          H exp(G)|0> on luhc
*              E       on elen
*              S       on ovl
*              G       on ccvec1 
*
*          imode: num. integration scheme
*          npnts: number of integration points
*
*   note on scratch vectors: ccvec3 is only needed if iopsym.eq.+/-1
*
*   output:   gradE      on grad
*            |gradE|     on xngrad
*
*     igtbmod.eq.1: use exp(G^2)
*
*  ak, early 2004
*
**********************************************************************
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
c      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cands.inc'
      include 'gtbce.inc'
* debugging:
      integer, parameter :: ntest = 005
      logical, parameter :: tstgrad = .false.

* input/output arrays
      integer, intent(in) ::
     &     igradmode, ioff_cc_typ(n_cc_typ), namp_cc_typ(n_cc_typ)
      real*8, intent(inout) ::
     &     grad(n_cc_amp)
* scratch arrays
      real*8 civec1(mxb_ci),civec2(mxb_ci),c2vec(*),
     &       ccvec1(n_cc_amp), ccvec2(n_cc_amp), ccvec3(n_cc_amp)
* local arrays
      character*8 cctype
      real*8 alp(npnts+2), wght(npnts+2)
* external functions
      real*8 inprod, inprdd

      call atim(cpu0,wall0)

      nprintl = max(ntest,nprint)

      lblk = -1
      if (ntest.ge.5) then
        write (6,*) '====================='
        write (6,*) ' This is gtbce_gradE'
        write (6,*) '====================='
        write (6,*) 
        write (6,*) 'on entry: '
        write (6,*) 'imode, npnts   : ', imode, npnts
        write (6,*) 'igradmode      : ', igradmode
        write (6,*) 'isymmet_G      : ', isymmet_G
        write (6,*) 'ovl, elen: ',ovl,elen
        write (6,*) 'n_cc_amp,mxb_ci : ', n_cc_amp,mxb_ci
        write (6,*) 'luc,luec,luhc,lusc1,lusc2: ',
     &               luc,luec,luhc,lusc1,lusc2
      end if

      if (ntest.ge.1000) then
        write(6,*) 'on entry:'
        write(6,*) 'Reference on LUC'
        call wrtvcd(civec1,luc,1,lblk)
        write(6,*) 'e^G|0> on LUEC'
        call wrtvcd(civec1,luec,1,lblk)
        write(6,*) 'H e^G|0> on LUHC'
        call wrtvcd(civec1,luhc,1,lblk)
      end if

c      if (ntest.ge.5) then
c        xnorm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))
c        write (6,*) 'norm of T: ',xnorm
c      end if
c      if (ntest.ge.100) then
c        call wrt_cc_vec2(ccvec1,6,cctype)
c      end if
      
      ! for I/O
      lblk = -1
      ! for expt_ref
      thresh=expg_thrsh
      mx_term=-mxterm_expg
      cctype='GEN_CC'
*--------------------------------------------------------------------*
* set up points and weights
*--------------------------------------------------------------------*
      select case (imode)
      case (0)  ! just testing 
        do ipnt = 1, npnts
          alp(ipnt) = dble(ipnt-1)/dble(npnts-1)
          wght(ipnt) = 1d0
        end do
      case (1)  ! Gauss-Legendre
        call gl_weights(0d0,1d0,npnts,alp,wght)
      case (2)  ! Simpson
c        if (mod(npnts,2).eq.0) npnts = npnts-1
        call s_weights(0d0,1d0,npnts,alp,wght)
      case default
        stop 'unknown imode in gtbce_gradE'
      end select
c      call test_quad(0d0,1d0,npnts,alp,wght)
c      stop 'enf stop after quad'

      mxpnts=npnts
      ! if G == 0 ...
      xnrm2 = inprod(ccvec1,ccvec1,n_cc_amp)
      ! ... things are trivial and we evaluate the formula only once
      if (xnrm2.lt.10d-20) then
        mxpnts=1
        wght(1)=1d0
        alp(1)=0d0
        if (ntest.ge.5) then
          write(6,*) 'Detected zero amplitudes: ',
     &               'only case alpha = 0 will be processed'
        end if
      else if (tstgrad) then
        ! does not work in route 3!
        mxpnts = npnts+2
        wght(npnts+1)=0d0
        wght(npnts+2)=0d0
        alp(npnts+1)=0d0
        alp(npnts+2)=1d0        
      end if
      call setvec(grad,0d0,n_cc_amp)

*--------------------------------------------------------------------*
* (H-E)|0tilde>
*  result on lusc1
*--------------------------------------------------------------------*
      if (igradmode.eq.1) then      ! (H-E)|0tilde>
        call vecsmdp(civec1,civec2,1d0,-elen,luhc,luec,lusc1,1,lblk)
      else if (igradmode.eq.2) then !     H|0tilde> only
        call copvcd(luhc,lusc1,civec1,1,lblk)
      else if (igradmode.eq.3) then !      |0tilde> only
        call copvcd(luec,lusc1,civec1,1,lblk)
      end if

**-------------------------------------------------------------------*
* loop over quadrature points
**-------------------------------------------------------------------*
      do ipnt = 1, mxpnts
        if (ntest.ge.5) then
          write(6,*) 'info for quadrature point: ', ipnt,'/',npnts
          write(6,*) 'point, weight: ', alp(ipnt), wght(ipnt)
        end if

        if (ipnt.gt.1.and.(alp(ipnt).le.alp(ipnt-1))) then
          write(6,*) 'quadrature point should be in ascending order!'
          stop 'gtbce_gradE > quadrature '
        end if

        if (ipnt.eq.1) then
          dltalp = alp(1)
        else
          dltalp = alp(ipnt)-alp(ipnt-1)
          call copvcd(lusc2,lusc1,civec1,1,lblk)
        end if
*--------------------------------------------------------------------*
* |a_i> = exp(a_i G^+) [(H-E)exp(G)|0>]
*       = exp((a_i-a_{i-1})G^+) [exp(a_{i-1}G^+) (H-E)exp(G)|0>]
*  result on lusc2
*--------------------------------------------------------------------*
        if (ntest.ge.5) then
          write(6,*)
     &         'constructing |a_i> = exp(a_i G^+) [(H-E)exp(G)|0>]'
        end if
        
        if (abs(dltalp).lt.1d-20) then
          call copvcd(lusc1,lusc2,civec1,1,lblk)
        else
          ! get the conjugate operator G^+ on ccvec2
          call conj_ccamp(ccvec1,1,ccvec2)
          if (igtbmod.ne.1) then
            ! and scale it
            call scalve(ccvec2,dltalp,n_cc_amp)
            call conj_t
            call expt_ref2(lusc1,lusc2,lusc4,lusc5,lusc6,
     &         thresh,mx_term, ccvec2, ccvec3, civec1, civec2,          
     &         n_cc_amp,cctype, iopsym)
            call conj_t
          else
            call conj_t
            call expt2_ref(lusc1,lusc2,lusc4,lusc5,lusc6,
     &         thresh,mx_term,
     &         dltalp,ccvec2, ccvec3, civec1, civec2,n_cc_amp,
     &         iopsym)
            call conj_t            
          end if
          if (ntest.ge.5) then
            xnrm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
            etest = inprdd(civec1,civec2,luc,lusc2,1,lblk)
            write(6,*) '|dlta G^+|, dlta = ',xnrm, dltalp
            write(6,*) '<ref|a_i> = ', etest,
     &                 'for alp(i) = ', alp(ipnt) 
          end if

        end if

*--------------------------------------------------------------------*
* |b_i> = exp(-a_i G)exp(G)|0> =
*       = exp(-(a_i-a_{i-1})G) [exp(-a_{i-1}G)exp(G)|0>]
*  result on lusc3
*--------------------------------------------------------------------*
        if (ipnt.eq.1) then
          call copvcd(luec,lusc1,civec1,1,lblk)
        else
          call copvcd(lusc3,lusc1,civec1,1,lblk)
        end if

        if (ntest.ge.5) then
          write(6,*) 'constructing |b_i> = exp(-a_i G) exp(G)|0>]'
        end if

        if (abs(dltalp).lt.1d-20) then
          call copvcd(lusc1,lusc3,civec1,1,lblk)          
        else
          if (igtbmod.ne.1) then
            ! get a copy of G
            call copvec(ccvec1,ccvec2,n_cc_amp)
            ! and scale it
            call scalve(ccvec2,-dltalp,n_cc_amp)
            call expt_ref2(lusc1,lusc3,lusc4,lusc5,lusc6,
     &           thresh,mx_term, ccvec2, ccvec3, civec1, civec2,
     &           n_cc_amp,cctype, iopsym)
          else
            call expt2_ref(lusc1,lusc3,lusc4,lusc5,lusc6,
     &              thresh,mx_term,
     &              -dltalp,ccvec1, ccvec3, civec1, civec2,n_cc_amp,
     &              iopsym)
          end if
          if (ntest.ge.5) then
            xnrm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
            etest = inprdd(civec1,civec2,lusc3,lusc3,1,lblk)
            etest2= inprdd(civec1,civec2,lusc2,lusc3,1,lblk)
            write(6,*) '|dltaG|, dlta = ',xnrm, dltalp
            write(6,*) '<b_i|b_i> , S = ', etest, ovl,
     &           'for alp(i) = ', alp(ipnt) 
            write(6,*) '<a_i|b_i>     = ', etest2,
     &           'for alp(i) = ', alp(ipnt) 
          end if
        end if
        
*--------------------------------------------------------------------*
* dE_u +=  w_i <a_i|gamma_u|b_i>
*  note: sigden implements ccvec2 = <lusc2|gamma_u|lusc3>
*
* for exp(G^2) we have
*
*  dE_u += w_i ( <a_i|G gamma_u|b_i> + <a_i|gamma_u G|b_i> ) 
*
*--------------------------------------------------------------------*

        if (ntest.ge.1000) then
          write(6,*) 'Before calling sigden_cc:'
          write(6,*) '|a_i> on lusc2:'
          call wrtvcd(civec1,lusc2,1,lblk)
          write(6,*) '|b_i> on lusc3:'
          call wrtvcd(civec1,lusc3,1,lblk)
        end if

        if (igtbmod.ne.1) then
          isigden=2
          ccvec2(1:n_cc_amp)=0d0
          call sigden_cc(civec1,civec2,lusc3,lusc2,ccvec2,isigden)

          call vecsum(grad,grad,ccvec2,1d0,wght(ipnt),n_cc_amp)

          if (ntest.ge.150) then
            xnorm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
            write(6,*)
     &           'non-weighted contrib to gradient: norm = ', xnorm
            if (iopsym.ne.0) write(6,*)
     &           ' (from non-conjugated exc. op.)'
            call wrt_cc_vec2(ccvec2,6,cctype)
            if (imode.eq.0) then
              ist = 1
              do
                ind = min(ist+19,n_cc_amp)
                if (ind-ist.gt.0) write(6,*) '@@ ',ist,ind,alp(ipnt),
     &               grad(ist:ind)
                if (ind.ge.n_cc_amp) exit
                ist = ist + 20
              end do
            end if
          end if

          if (iopsym.eq.1.or.iopsym.eq.-1) then
            ccvec2(1:n_cc_amp)=0d0
            call conj_t
            call sigden_cc(civec1,civec2,lusc3,lusc2,ccvec2,isigden)
            call conj_ccamp(ccvec2,1,ccvec3)
            call conj_t
            fac = wght(ipnt)
            if (iopsym.eq.-1) fac = -wght(ipnt)
            call vecsum(grad,grad,ccvec3,1d0,fac,n_cc_amp)
            
            if (ntest.ge.150) then
              xnorm = sqrt(inprod(ccvec3,ccvec3,n_cc_amp))
              write(6,*)
     &             'non-weighted contrib to gradient: norm = ', xnorm
              write(6,*)' (from conjugated exc. op.)'
              call wrt_cc_vec2(ccvec3,6,cctype)
              if (imode.eq.0) then
                ist = 1
                do
                  ind = min(ist+19,n_cc_amp)
                  if (ind-ist.gt.0) write(6,*) '@@ ',ist,ind,alp(ipnt),
     &                 grad(ist:ind)
                  if (ind.ge.n_cc_amp) exit
                  ist = ist + 20
                end do
              end if
            end if
          end if

        else ! exp(G^2) part:
*    G^+ |a> on lusc4
          isigden=1
          call conj_ccamp(ccvec1,1,ccvec2)
          call conj_t
          call sigden_cc(civec1,civec2,lusc2,lusc4,ccvec2,isigden)
          call conj_t
        
*    <a| G gamma |b> contribution:
*  note: sigden implements ccvec2 = <lusc4|gamma_u|lusc3>
          isigden=2
          ccvec2(1:n_cc_amp)=0d0
          call sigden_cc(civec1,civec2,lusc3,lusc4,ccvec2,isigden)

*  increment gradient:
          call vecsum(grad,grad,ccvec2,1d0,wght(ipnt),n_cc_amp)
        
*    G |b> on lusc4
          isigden=1
          call sigden_cc(civec1,civec2,lusc3,lusc4,ccvec1,isigden)
        
*    <a| gamma G |b> contribution:
          isigden=2
          ccvec2(1:n_cc_amp)=0d0
          call sigden_cc(civec1,civec2,lusc4,lusc2,ccvec2,isigden)

*  increment gradient:
          call vecsum(grad,grad,ccvec2,1d0,wght(ipnt),n_cc_amp)

          if (ntest.ge.150) then
            xnorm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
            write(6,*)
     &           'non-weighted contrib to gradient: norm = ', xnorm
            if (iopsym.ne.0) write(6,*)
     &           ' (from non-conjugated exc. op.)'
            call wrt_cc_vec2(ccvec2,6,cctype)
          end if

          if (iopsym.eq.1.or.iopsym.eq.-1)
     &      stop 'not prepared for iopsym.ne.0'

        end if

      end do


      if (isymmet_G.ne.0) then
        if (ntest.ge.1000) then
          write(6,*) 'The new gradient (bef. symmetrizing):'
          call wrt_cc_vec2(grad,6,'GEN_CC')
        end if
        call symmet_t(isymmet_G,1,
     &       grad,ccvec2,
     &       ictp,i_cc_typ,n_cc_typ,
     &       namp_cc_typ,ioff_cc_typ,ngas)
      end if

      if (igradmode.eq.1) then
        ! normalize gradient
        call scalve(grad,2d0/ovl,n_cc_amp)
        xngrad = sqrt(inprod(grad,grad,n_cc_amp))
      end if

      if (ntest.ge.5) then
        write(6,*) ' gtbce_gradE > '
        write(6,*) '     n_cc_amp,norm of grad: ',n_cc_amp,xngrad
      end if
      if (ntest.ge.100) then
        call wrt_cc_vec2(grad,6,'GEN_CC')
      end if

      if (nprintl.ge.1) then
        write(6,'(4(/x,a))')
     &   ' Contributions to gradient norm per operator type:',
     &   '-----------------------------------------------------------',
     &   '   type     n      norm     norm/n        max        min',
     &   '-----------------------------------------------------------'
        do itp = 1, n_cc_typ
          ist = ioff_cc_typ(itp)
          len = namp_cc_typ(itp)
          xnorm = sqrt(inprod(grad(ist),grad(ist),len))
          xmax = fndmnx(grad(ist),len,2)
          xmin = fndmnx(grad(ist),len,1)
          write(6,'(4x,i3,x,i7,4(x,e10.4))')
     &      itp,len,xnorm,xnorm/dble(len),xmax,xmin
        end do
        write(6,'(x,a,/)')
     &   '-----------------------------------------------------------'

      end if
      
      call atim(cpu,wall)
      call prtim(6,'time in gtbce_gradE',cpu-cpu0,wall-wall0)

      return
      end
*--------------------------------------------------------------------*
* DECK: gtbce_tstgradE
*--------------------------------------------------------------------*
      subroutine gtbce_testgradE(igtbmod,isymmet_G,igtb_closed,
     &                       ccvec1,ccvec2,xngrad_num,
     &                       ecore,
     &                       ccvec3,iopsym,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_typ,i_cc_typ,namp_cc_typ,ioff_cc_typ,
     &                       n_cc_amp,mxb_ci,
     &                       n11amp,n33amp,iamp_packed,ictp,
     &                       luamp,lugrd,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)
*--------------------------------------------------------------------*
*
* test gradient by numerical differentiation
* the exact gradient should be passed
*
*--------------------------------------------------------------------*
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
c      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cands.inc'
* debugging:
      integer, parameter :: ntest = 1000000

* input/output arrays
      real*8 ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp)
* scratch arrays
      real*8 ::
     &     civec1(mxb_ci),civec2(mxb_ci),c2vec(*)
      real*8 ::
     &     ccvec3(n_cc_amp), ccvec4(n_cc_amp)
* external functions
      real*8 ::
     &     inprod

      write (6,'(/,3(x,a,/))')
     &     '============================',
     &     ' Welcome to gtbce_tstgradE!',
     &     '============================'

* increment is 0.001
      xinc = 0.00001d0 

      if (igtb_closed.eq.0) then
        namp = n_cc_amp
      else
        namp_packed = n11amp+n33amp
        namp = namp_packed
      end if

      do iamp = 1, namp

        if (igtb_closed.eq.0) then
          call vec_from_disc(ccvec3,namp,1,-1,luamp)
* increment +
          ccvec3(iamp) = ccvec3(iamp) + xinc
        else
          call vec_from_disc(ccvec1,namp,1,-1,luamp)          
* increment +
          ccvec1(iamp) = ccvec1(iamp) + xinc
          iway = -1
          idual = 3
          call pack_g(iway,idual,isymmet_G,ccvec1,ccvec3,
     &                n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                n11amp,n33amp,iamp_packed,n_cc_amp)          

            call chksym_t(isymmet_G,1,
     &           ccvec3,ccvec1,
     &           ictp,i_cc_typ,n_cc_typ,
     &           namp_cc_typ,ioff_cc_typ,ngas)

        end if

        call gtbce_E(igtbmod,elenp,varp,ovl,
     &               ecore,
     &               ccvec3,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

* increment -
        if (igtb_closed.eq.0) then
          call vec_from_disc(ccvec3,namp,1,-1,luamp)          
          ccvec3(iamp) = ccvec3(iamp) - xinc
        else
          call vec_from_disc(ccvec1,namp,1,-1,luamp)          
          ccvec1(iamp) = ccvec1(iamp) - xinc
          iway = -1
          idual = 3
          call pack_g(iway,idual,isymmet_G,ccvec1,ccvec3,
     &                n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                n11amp,n33amp,iamp_packed,n_cc_amp)          
            call chksym_t(isymmet_G,1,
     &           ccvec3,ccvec1,
     &           ictp,i_cc_typ,n_cc_typ,
     &           namp_cc_typ,ioff_cc_typ,ngas)
        end if
        call gtbce_E(igtbmod,elenm,varm,ovl,
     &               ecore,
     &               ccvec3,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

* compare
        gradnum = (elenp-elenm)/(2d0*xinc)
        ccvec2(iamp) = gradnum
        call vec_from_disc(ccvec1,namp,1,-1,lugrd)          
        if (ntest.gt.150) then
          write(6,'(/,x,a,/x,a,i6,/x,a,3(/x,a,e20.13)/)')
     &       '==================================',
     &       ' RESULT FOR IAMP = ',iamp,
     &       '==================================',
     &       ' analytic ',ccvec1(iamp),
     &       ' numeric  ',gradnum,
     &       ' diff     ',ccvec1(iamp)-gradnum
          if (gradnum.ne.0d0)
     &         write(6,*)
     &       ' a/n      ',ccvec1(iamp)/gradnum
          if (ccvec1(iamp).ne.0d0)
     &         write(6,*)
     &       ' n/a      ',gradnum/ccvec1(iamp)
        end if

      end do

      write (6,*) 'comparison of analytical and numerical gradient:'
      call cmp2vc(ccvec1,ccvec2,namp,.1d-2*xinc*xinc)

      xngrad_num = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))

      return

      end
*--------------------------------------------------------------------*
* stop card: ccvec1, ccvec2, grad, grad_num
*--------------------------------------------------------------------*
* DECK: gtbce_tstgradE_L
*--------------------------------------------------------------------*
      subroutine gtbce_testgradE_L(
     &                       gradL,ampL,
     &                       ecore,
     &                       ccvec1,iopsym,ccvec2,
     &                       civec1,civec2,c2vec,
     &                       n_cc_amp,n_l_amp,mxb_ci,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)
*--------------------------------------------------------------------*
*
* test gradient by numerical differentiation
* the exact gradient should be passed
*
*--------------------------------------------------------------------*
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
c      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cands.inc'
      include 'glbbas.inc'
* debugging:
      integer, parameter :: ntest = 1000000

* input/output arrays
      real*8 ::
     &     gradL(*), ampL(*)
* scratch arrays
      real*8 ::
     &     civec1(mxb_ci),civec2(mxb_ci),c2vec(*)
      real*8 ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp)
* external functions
      real*8 ::
     &     inprod

      write (6,'(/,3(x,a,/))')
     &     '=============================',
     &     ' Welcome to gtbce_tstgradE_L',
     &     '============================='

* increment is 0.001
      xinc = 0.0001d0 

      do iamp = 1, n_l_amp
        
* increment +
        ampL(iamp) = ampL(iamp) + xinc

        call l2g(ampL,ccvec1,nspobex_tp,
     &       work(klsobex),work(klibsobex),0  ,ntoob)

        igtbmod = 2 ! obviously
        call gtbce_E(igtbmod,elenp,varp,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec2,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

* increment -
        ampL(iamp) = ampL(iamp) - 2d0*xinc

        call l2g(ampL,ccvec1,nspobex_tp,
     &       work(klsobex),work(klibsobex),0  ,ntoob)

        call gtbce_E(igtbmod,elenm,varm,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec2,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

* reset
        ampL(iamp) = ampL(iamp) + xinc

* compare
        gradnum = (elenp-elenm)/(2d0*xinc)
        if (ntest.gt.150) then
          ii = iamp/ntoob + 1
          jj = mod(iamp-1,ntoob) + 1
          write(6,'(/,x,a,/x,a,i6,x,i3,x,i3,/x,a,3(/x,a,e20.13)/)')
     &       '==================================',
     &       ' RESULT FOR IAMP = ',iamp,ii,jj,
     &       '==================================',
     &       ' analytic ',gradL(iamp),
     &       ' numeric  ',gradnum,
     &       ' diff     ',gradL(iamp)-gradnum
        end if

      end do

c      write (6,*) 'comparistion of analytical and numerical gradient:'
c      call cmp2vc(grad,grad_num,n_cc_amp,.1d-2*xinc*xinc)
c
c      xngrad_num = sqrt(inprod(grad_num,grad_num,n_cc_amp))

      return

      end
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
      subroutine gtbce_testgradE_UOU(imode,
     &                       grad,omvec,urvec,uivec,
     &                       elen,ecore,
     &                       ccvec1,iopsym,ccvec2,
     &                       civec1,civec2,c2vec,
     &                       n_cc_amp,n_l_amp,mxb_ci,
     &                       luc,luec,luhc,
     &                       lusc1,lusc2)
*--------------------------------------------------------------------*
*
* test gradient by numerical differentiation
* the exact gradient should be passed
*
*--------------------------------------------------------------------*
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
c      include 'crun.inc'
      include 'cstate.inc'
      include 'cgas.inc'
      include 'ctcc.inc'
      include 'gasstr.inc'
      include 'strinp.inc'
      include 'orbinp.inc'
      include 'cprnt.inc'
      include 'corbex.inc'
      include 'csm.inc'
      include 'cands.inc'
      include 'glbbas.inc'
* debugging:
      integer, parameter :: ntest = 1000000

* input/output arrays
      real*8 ::
     &     grad(*), omvec(*), urvec(*), uivec(*)
* scratch arrays
      real*8 ::
     &     civec1(mxb_ci),civec2(mxb_ci),c2vec(*)
      real*8 ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp)
* external functions
      real*8 ::
     &     inprod

      write (6,'(/,3(x,a,/))')
     &     '===============================',
     &     ' Welcome to gtbce_tstgradE_UOU',
     &     '==============================='
      write(6,*) ' imode = ', imode
      write(6,*) ' number of amplitudes = ', n_l_amp

      call uou2g(omvec,urvec,uivec,ccvec1,
     &         nspobex_tp,
     &         work(klsobex),work(klibsobex),ntoob)
        igtbmod = 3 ! obviously
        call gtbce_E(igtbmod,elen0,varp,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec2,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

* increment is 0.001
      xinc = 0.0001d0 

      write(6,*) ' xinc = ',xinc
      write(6,*) ' elen0 = ',elen0

      do iamp = 1, n_l_amp
        
* increment +
        if (imode.eq.1) omvec(iamp) = omvec(iamp) + xinc
        if (imode.eq.2) urvec(iamp) = urvec(iamp) + xinc
        if (imode.eq.3) uivec(iamp) = uivec(iamp) + xinc

        call uou2g(omvec,urvec,uivec,ccvec1,
     &         nspobex_tp,
     &         work(klsobex),work(klibsobex),ntoob)

        igtbmod = 3 ! obviously
        call gtbce_E(igtbmod,elenp,varp,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec2,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

* increment -
        if (imode.eq.1) omvec(iamp) = omvec(iamp) - 2d0*xinc
        if (imode.eq.2) urvec(iamp) = urvec(iamp) - 2d0*xinc
        if (imode.eq.3) uivec(iamp) = uivec(iamp) - 2d0*xinc

        call uou2g(omvec,urvec,uivec,ccvec1,
     &         nspobex_tp,
     &         work(klsobex),work(klibsobex),ntoob)

        call gtbce_E(igtbmod,elenm,varm,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec2,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

* reset
        if (imode.eq.1) omvec(iamp) = omvec(iamp) + xinc
        if (imode.eq.2) urvec(iamp) = urvec(iamp) + xinc
        if (imode.eq.3) uivec(iamp) = uivec(iamp) + xinc

* compare
        gradnum = (elenp-elenm)/(2d0*xinc)

        hessnum = (elenp+elenm - 2d0*elen0)/(xinc*xinc)

        if (ntest.gt.150) then
          ii = iamp/ntoob + 1
          jj = mod(iamp-1,ntoob) + 1
          write(6,'(/,x,a,/x,a,i6,x,i3,x,i3,/x,a,4(/x,a,e20.13)/)')
     &       '==================================',
     &       ' RESULT FOR IAMP = ',iamp,ii,jj,
     &       '==================================',
     &       ' analytic ',grad(iamp),
     &       ' numeric  ',gradnum,
     &       ' diff     ',grad(iamp)-gradnum,
     &       ' num.hess ',hessnum
        end if

      end do

c      write (6,*) 'comparistion of analytical and numerical gradient:'
c      call cmp2vc(grad,grad_num,n_cc_amp,.1d-2*xinc*xinc)
c
c      xngrad_num = sqrt(inprod(grad_num,grad_num,n_cc_amp))

      return

      end
*--------------------------------------------------------------------*
* DECK: gtbce_num2drv
*--------------------------------------------------------------------*
      subroutine gtbce_num2drv(igtbmod,imode,iomggrd,
     &                       igtb_closed,isymmet_G,
     &                       inumint,npnts,
     &                       ecore,
     &                       iccvec,nSdim,
     &                       ccvec1,iopsym,ccvec2,ccvec3,ccvec4,
     &                       civec1,civec2,c2vec,
     &                       n_cc_typ,i_cc_typ,ictp,
     &                       namp_cc_typ,ioff_cc_typ,
     &                       n_cc_amp,mxb_ci,
     &                       n11amp,n33amp,iamp_packed,
     &                       luhss,
     &                       luamp,luleqv,luc,luec,luhc,
     &                       lusc1,lusc2,lusc3,lusc4,lusc5,lusc6,lusc7)
*--------------------------------------------------------------------*
*
* purpose: calculate the numerical second derivatives of E/S or
*          the Jacobian dOmg/dG resp.ly
*
*  imode = 1   get matrix-vector product  
*          2   calculate complete H-ES matrix
*          3   calculate complete H matrix
*          4   calculate complete S matrix
*
*  iomggrd = 0  calc. Omega
*            1  calc. Gradient
*
*  ak, early 2004
*
*--------------------------------------------------------------------*
      include 'implicit.inc'

* constants
      integer, parameter ::
     &     ntest = 010

* input

* scratch
      real*8, intent(inout) ::
     &     ccvec1(n_cc_amp),ccvec2(n_cc_amp),
     &     ccvec3(n_cc_amp),ccvec4(n_cc_amp),
     &     civec1(*), civec2(*), c2vec(*)
      integer, intent(inout) ::
     &     iccvec(n_cc_amp)

      real(8), external ::
     &     inprod

      lblk = -1
      xinc = 1d-5

      if (ntest.gt.0) then
        write(6,*) '======================='
        write(6,*) ' This is gtbce_num2drv'
        write(6,*) '======================='
        write(6,*) ' imode = ',imode
        write(6,*) ' xinc = ',xinc
        write(6,*) ' luhss,luamp,luleqv: ',luhss,luamp,luleqv
        write(6,*) ' igtbmod,isymmet_G,igtb_closed,iopsym: ',
     &       igtbmod,isymmet_G,igtb_closed,iopsym
        if (igtb_closed.ne.0) then
          write(6,*) 'n11amp, n33amp: ',n11amp,n33amp          
        end if
      end if

c      if (imode.gt.1) igradmode=imode-1
      igradmode = 1

      namp = n_cc_amp
      if (igtb_closed.eq.1) namp = n11amp+n33amp
      nloops = namp
      if (imode.eq.1) nloops = 1

      ! rewind output file
      call rewino(luhss)

* loop over elements in vector
      do iloop = 1, nloops

        if (imode.ne.1.and.igtb_closed.ne.1.and.isymmet_G.ne.0) then
          if (iccvec(iloop).lt.0) cycle
        end if

        if (ntest.ge.5) then
          write(6,*) 'iloop = ',iloop,'/',nloops
        end if

        ! reload amplitudes
        if (igtb_closed.eq.0) then
          call vec_from_disc(ccvec1,namp,1,-1,luamp)
        else
          call vec_from_disc(ccvec2,namp,1,-1,luamp)
        end if

* inc + xinc
        if (ntest.ge.10) then
          write(6,*) '------------------'
          write(6,*) 'positive increment'
          write(6,*) '------------------'
        end if
        call memchk2('zzz---')
        if (imode.eq.1.and.igtb_closed.eq.0) then
          call vec_from_disc(ccvec2,namp,1,-1,luleqv)
          if (ntest.ge.100.and.imode.eq.1) then
            xnorm = sqrt(inprod(ccvec2,ccvec2,namp))
            write(6,*) ' norm of input vector = ',xnorm
          end if
          ccvec1(1:namp) =
     &         ccvec1(1:namp)+xinc*ccvec2(1:namp)
        else if (imode.eq.1.and.igtb_closed.ne.0) then
          call vec_from_disc(ccvec1,namp,1,-1,luleqv)
          if (ntest.ge.100.and.imode.eq.1) then
            xnorm = sqrt(inprod(ccvec1,ccvec1,namp))
            write(6,*) ' norm of input vector = ',xnorm
          end if
          ccvec2(1:namp) =
     &         ccvec2(1:namp)+xinc*ccvec1(1:namp)
        else if (isymmet_G.eq.0) then
          ccvec1(iloop) = ccvec1(iloop) + xinc
        else if (isymmet_G.ne.0.and.igtb_closed.ne.0) then
          call memchk2('yyy---')
          ccvec2(iloop) = ccvec2(iloop) + xinc
          call memchk2('xxx---')
        else
          iadj = abs(iccvec(iloop))
          fac = dble(isymmet_G)
          ccvec1(iloop) = ccvec1(iloop) + sqrt(2d0)*xinc
          ccvec1(iadj) = ccvec1(iadj) + fac*sqrt(2d0)*xinc
        end if

        call memchk2('aaa---')

        if (igtb_closed.ne.0) then
          iway = -1
          call pack_g(iway,idum,isymmet_G,ccvec2,ccvec1,
     &                n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                n11amp,n33amp,iamp_packed,n_cc_amp)
        end if
        call memchk2('bbb---')

        call gtbce_E(igtbmod,elen,variance,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)
        if (iomggrd.eq.0) then
          call gtbce_Omg(ccvec3,xnomg,
     &                   elen,ovl,iopsym,
     &                   civec1,civec2,c2vec,
     &                   n_cc_amp,mxb_ci,
     &                   luec,luhc,lusc1,lusc2)
        else
          ipr=0
          call gtbce_gradE(isymmet_G,ccvec3,xngrad,igradmode,
     &                 inumint,npnts,
     &                 elen,ovl,
     &                 ccvec1,iopsym,ccvec2,ccvec4,
     &                 civec1,civec2,c2vec,
     &                 n_cc_typ,i_cc_typ,ictp,
     &                 namp_cc_typ,ioff_cc_typ,
     &                 n_cc_amp,mxb_ci,ipr,
     &                 luamp,luc,luec,luhc,
     &                 lusc1,lusc2,lusc3,lusc4,lusc5,lusc6)
        end if

* save
        if (ntest.ge.1000) then
          write (6,*) 'gradient for positive increment:'
          call wrt_cc_vec2(ccvec3,6,'GEN_CC')
        end if
        call vec_to_disc(ccvec3,n_cc_amp,1,lblk,lusc7)

* inc - xinc
        if (ntest.ge.10) then
          write(6,*) '------------------'
          write(6,*) 'negative increment'
          write(6,*) '------------------'
        end if
        if (imode.eq.1.and.igtb_closed.eq.0) then
          call vec_from_disc(ccvec2,namp,1,-1,luleqv)
          ccvec1(1:namp) =
     &         ccvec1(1:namp)-2d0*xinc*ccvec2(1:namp)
        else if (imode.eq.1.and.igtb_closed.eq.1) then
          call vec_from_disc(ccvec2,namp,1,-1,luamp)
          call vec_from_disc(ccvec1,namp,1,-1,luleqv)
          ccvec2(1:namp) =
     &         ccvec2(1:namp)-2d0*xinc*ccvec1(1:namp)
        else if (isymmet_G.eq.0) then
          ccvec1(iloop) = ccvec1(iloop) - 2d0*xinc
        else if (isymmet_G.ne.0.and.igtb_closed.ne.0) then
          call vec_from_disc(ccvec2,namp,1,-1,luamp)
          ccvec2(iloop) = ccvec2(iloop) - xinc
        else
          iadj = abs(iccvec(iloop))
          fac = dble(isymmet_G)
          ccvec1(iloop) = ccvec1(iloop) - 2d0*sqrt(2d0)*xinc
          ccvec1(iadj) = ccvec1(iadj) - fac*2d0*sqrt(2d0)*xinc
        end if

        if (igtb_closed.ne.0) then
          iway = -1
          call pack_g(iway,idum,isymmet_G,ccvec2,ccvec1,
     &                n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                n11amp,n33amp,iamp_packed,n_cc_amp)
        end if

        call gtbce_E(igtbmod,elen,variance,ovl,
     &               ecore,
     &               ccvec1,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)
        if (iomggrd.eq.0) then
          call gtbce_Omg(ccvec3,xnomg,
     &                   elen,ovl,iopsym,
     &                   civec1,civec2,c2vec,
     &                   n_cc_amp,mxb_ci,
     &                   luec,luhc,lusc1,lusc2)
        else
          ipr=0
          call gtbce_gradE(isymmet_G,ccvec3,xngrad,igradmode,
     &                 inumint,npnts,
     &                 elen,ovl,
     &                 ccvec1,iopsym,ccvec2,ccvec4,
     &                 civec1,civec2,c2vec,
     &                 n_cc_typ,i_cc_typ,ictp,
     &                 namp_cc_typ,ioff_cc_typ,
     &                 n_cc_amp,mxb_ci,ipr,
     &                 luamp,luc,luec,luhc,
     &                 lusc1,lusc2,lusc3,lusc4,lusc5,lusc6)
        end if

        if (ntest.ge.1000) then
          write (6,*) 'gradient for negative increment:'
          call wrt_cc_vec2(ccvec3,6,'GEN_CC')
        end if

* get difference
        call vec_from_disc(ccvec2,n_cc_amp,1,lblk,lusc7)
        fac = 1d0/(2d0*xinc)
        call vecsum(ccvec3,ccvec3,ccvec2,-fac,fac,n_cc_amp)

        if (ntest.ge.500) then
          write(6,*) 'result for iloop = ', iloop
          call wrt_cc_vec2(ccvec3,6,'GEN_CC')
        end if

        if (isymmet_G.ne.0.and.igtb_closed.eq.0.and.imode.ne.1) then
          ! compress result vector
          idx = 0
          do ii = 1, n_cc_amp
            if (iccvec(ii).le.0) cycle
            idx = idx + 1
            ccvec2(idx) = 2d0*ccvec3(ii)
          end do
          if (idx.ne.nSdim) stop 'verdacht'
          call vec_to_disc(ccvec2,nSdim,0,lblk,luhss)
        else if (igtb_closed.ne.0) then
          iway = 2
          call pack_g(iway,idum,isymmet_G,ccvec1,ccvec3,
     &                n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                n11amp,n33amp,iamp_packed,n_cc_amp)          
          if (imode.eq.1.and.ntest.ge.100) then
            xnorm = sqrt(inprod(ccvec1,ccvec1,namp))
            write(6,*) ' norm of MV-product: ',xnorm
          end if
          call vec_to_disc(ccvec1,namp,0,lblk,luhss)
        else
* save result
          if (imode.eq.1.and.ntest.ge.100) then
            xnorm = sqrt(inprod(ccvec3,ccvec3,n_cc_amp))
            write(6,*) ' norm of MV-product: ',xnorm
          end if
c          if (imode.eq.1) call scalve(ccvec3,-1d0,n_cc_amp)
          call vec_to_disc(ccvec3,n_cc_amp,0,lblk,luhss)
        end if

      end do

      if (ntest.gt.0) then
        write(6,*) '======================'
        write(6,*) ' END OF gtbce_num2drv'
        write(6,*) '======================'
      end if

      return

      end
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
* DECK: gtbce_foo
*--------------------------------------------------------------------*
      subroutine gtbce_foo_old(inumint,npnts,
     &                     ovl,
     &                     ccvec1,iopsym,comm_ops,
     &                     ccvec2,ccvec3,
     &                     civec1,civec2,c2vec,
     &                     n_cc_amp,mxb_ci,
     &                     lufoo,
     &                     luamp,luc,luec,luhc,
     &                     lusc1,lusc2,lusc3,lusc4,
     &                     lusc5,lusc6,lusc7,lusc8)
*--------------------------------------------------------------------*
*
* purpose: Calculate the overlap of the first order wavefunction
*          change
*
*          S_ij = N <0|(d/dg_i exp(G^+))(d/dg_j exp(G))|0>
*
*  ak, early 2004
*
*--------------------------------------------------------------------*
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
      include 'gtbce.inc'
* debugging:
      integer, parameter :: ntest = 000
      logical, parameter :: tstgrad = .false.
 
* input/output arrays
      logical comm_ops
      integer, intent(in) ::
     &     inumint, npnts
* scratch arrays
      real*8 civec1(mxb_ci),civec2(mxb_ci),c2vec(*),
     &       ccvec1(n_cc_amp), ccvec2(n_cc_amp), ccvec3(n_cc_amp)
* local arrays
      character*8 cctype
      real*8 alp(npnts+2), wght(npnts+2)
* external functions
      real*8 inprod, inprdd

      call atim(cpu0,wall0)

      lblk = -1
      if (ntest.ge.5) then
        write (6,*) '====================='
        write (6,*) ' This is gtbce_foo (old)'
        write (6,*) '====================='
        write (6,*) 
        write (6,*) 'on entry: '
        write (6,*) 'inumint, npnts   : ', inumint, npnts
        write (6,*) 'ovl, elen: ',ovl
        write (6,*) 'n_cc_amp,mxb_ci : ', n_cc_amp,mxb_ci
        write (6,*) 'luc,luec,luhc,lusc1,lusc2: ',
     &               luc,luec,luhc,lusc1,lusc2
      end if
      if (ntest.ge.1000) then
        write(6,*) 'on entry:'
        write(6,*) '|0> on LUC'
        call wrtvcd(civec1,luc,1,lblk)
        write(6,*) 'e^G|0> on LUEC'
        call wrtvcd(civec1,luec,1,lblk)
      end if
      
      ! for I/O
      lblk = -1
      ! for expt_ref
      thresh=expg_thrsh
      mx_term=-mxterm_expg
      cctype='GEN_CC'
*--------------------------------------------------------------------*
* set up points and weights
*--------------------------------------------------------------------*
      select case (inumint)
      case (0)  ! just testing 
        do ipnt = 1, npnts
          alp(ipnt) = dble(ipnt-1)/dble(npnts-1)
          wght(ipnt) = 1d0
        end do
      case (1)  ! Gauss-Legendre
        call gl_weights(0d0,1d0,npnts,alp,wght)
      case (2)  ! Simpson
c        if (mod(npnts,2).eq.0) npnts = npnts-1
        call s_weights(0d0,1d0,npnts,alp,wght)
      case default
        stop 'unknown inumint in gtbce_foo'
      end select

      mxpnts=npnts
      ! if G == 0 ...
      xnrm2 = inprod(ccvec1,ccvec1,n_cc_amp)
      ! ... things are trivial and we evaluate the formula only once
c      comm_ops = .false.
      if (xnrm2.lt.10d-20.or.comm_ops) then
        mxpnts=1
        wght(1)=1d0
        alp(1)=0.0d0
        if (ntest.ge.5) then
          write(6,*) 'Detected zero amplitudes: ',
     &               'only case alpha = 0 will be processed'
        end if
      end if

      ! rewind output file
      call rewino(lufoo)

**-------------------------------------------------------------------*
* loop i over parameters
**-------------------------------------------------------------------*
      do iamp = 1, n_cc_amp
        if (ntest.ge.10) write(6,*) 'iamp = ',iamp,'/',n_cc_amp

        ! reset |0tilde> = exp(G)|0>
        call copvcd(luec,lusc1,civec1,1,lblk)

**-------------------------------------------------------------------*
* loop over quadrature points
**-------------------------------------------------------------------*
        do ipnt = 1, mxpnts
          if (ntest.ge.5) then
            write(6,*) 'info for quadrature point: ', ipnt,'/',npnts
            write(6,*) 'point, weight: ', alp(ipnt), wght(ipnt)
          end if

          if (ipnt.gt.1.and.(alp(ipnt).le.alp(ipnt-1))) then
            write(6,*) 'quadrature points should be in ascending order!'
            stop 'gtbce_foo > quadrature '
          end if

          if (ipnt.eq.1) then
            dltalp = -alp(1)
          else
            dltalp = -alp(ipnt)+alp(ipnt-1)
          end if
*--------------------------------------------------------------------*
* |a_i>(1) = exp(-a_i G) [exp(G)|0>]
*          = exp(-(a_i-a_{i-1})G) exp(-a_{i-1} G) [exp(G)|0>]
*  result on lusc2
*--------------------------------------------------------------------*
          if (ntest.ge.5) then
            write(6,*)
     &           'constructing |a_i> = exp(-a_i G^+) exp(G)|0>]'
          end if
        
          if (abs(dltalp).lt.1d-20) then
            call copvcd(lusc1,lusc2,civec1,1,lblk)
          else
            ! get G on ccvec2
            call copvec(ccvec1,ccvec2,n_cc_amp)
            ! and scale it
            call scalve(ccvec2,dltalp,n_cc_amp)
            call expt_ref2(lusc1,lusc2,lusc4,lusc5,lusc6,
     &           thresh,mx_term, ccvec2, ccvec3, civec1, civec2,          
     &           n_cc_amp,cctype, iopsym)
            if (ntest.ge.5) then
              xnrm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
              etest = inprdd(civec1,civec2,luc,lusc2,1,lblk)
              write(6,*) '|dlta G^+|, dlta = ',xnrm, dltalp
              write(6,*) '<ref|a_i> = ', etest,
     &                   'for alp(i) = ', alp(ipnt) 
            end if
            ! save for next round
            call copvcd(lusc2,lusc1,civec1,1,lblk)
          end if

*--------------------------------------------------------------------*
* |a_i(ii)>(2) = tau_ii exp(-a_i G)[exp(G)|0>]
*  result on lusc3
*--------------------------------------------------------------------*
          ccvec2(1:n_cc_amp) = 0d0
          ccvec2(iamp) = 1d0
          isigden=1
          call sigden_cc(civec1,civec2,lusc2,lusc3,ccvec2,isigden)
          if (iopsym.ne.0) then
            fac = dble(iopsym)
            call conj_ccamp(ccvec2,1,ccvec3)
            call conj_t
            call sigden_cc(civec1,civec2,lusc2,lusc4,ccvec3,isigden)
            call conj_t
            call vecsmdp(civec1,civec2,1d0,fac,lusc3,lusc4,lusc5,1,lblk)
            call copvcd(lusc5,lusc3,civec1,1,lblk)
          end if

*--------------------------------------------------------------------*
* |a_i(ii)>(3) = exp(a_i G) tau_ii exp(-a_i G)[exp(G)|0>]
*  result on lusc2 again
*--------------------------------------------------------------------*
          if (abs(alp(ipnt)).lt.1d-20) then
            call copvcd(lusc3,lusc2,civec1,1,lblk)
          else
            ! get G on ccvec2
            call copvec(ccvec1,ccvec2,n_cc_amp)
            ! and scale it
            call scalve(ccvec2,alp(ipnt),n_cc_amp)
            call expt_ref2(lusc3,lusc2,lusc4,lusc5,lusc6,
     &           thresh,mx_term, ccvec2, ccvec3, civec1, civec2,          
     &           n_cc_amp,cctype, iopsym)
            if (ntest.ge.5) then
              xnrm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
              etest = inprdd(civec1,civec2,luc,lusc2,1,lblk)
              write(6,*) '|alp(i) G^+|, alp(i) = ',xnrm, alp(ipnt)
              write(6,*) '<ref|a_i> = ', etest,
     &                 'for alp(i) = ', alp(ipnt) 
            end if
          end if

          if (ntest.ge.2000) then
            write (6,*) 'contribution to 1st derivative of ',
     &           'wavefunction, element ',iamp,alp(ipnt)
            call wrtvcd(civec1,lusc2,1,lblk)
          end if

          ! update result on lusc8
          if (ipnt.gt.1) then
            call vecsmdp(civec1,civec2,1d0,wght(ipnt),lusc8,lusc2,
     &           lusc4,1,lblk)
            call copvcd(lusc4,lusc8,civec1,1,lblk)
          else
            call sclvcd(lusc2,lusc8,wght(ipnt),civec1,1,lblk)
          end if

        end do ! loop over quadrature

c TEST compare with numerical wavefunction derivative:
        itest = 1
        if (itest.eq.1) then
          xinc = 1d-4
          call copvec(ccvec1,ccvec2,n_cc_amp)
          ccvec2(iamp) = ccvec2(iamp)+xinc
          call expt_ref2(luc,lusc2,lusc4,lusc5,lusc6,
     &         thresh,mx_term, ccvec2, ccvec3, civec1, civec2,          
     &         n_cc_amp,cctype, iopsym)
          call copvec(ccvec1,ccvec2,n_cc_amp)
          ccvec2(iamp) = ccvec2(iamp)-xinc
          call expt_ref2(luc,lusc3,lusc4,lusc5,lusc6,
     &         thresh,mx_term, ccvec2, ccvec3, civec1, civec2,          
     &         n_cc_amp,cctype, iopsym)
          fac=1d0/(2d0*xinc)
          call vecsmdp(civec1,civec2,fac,-fac,lusc2,lusc3,
     &           lusc4,1,lblk)
          print *,'==============================================='
          print *,' RESULT for iamp = ',iamp
          print *,' analytic 1st der. of WF: norm = ',
     &         sqrt(inprdd(civec1,civec2,lusc8,lusc8,1,lblk))
          print *,'  numeric 1st der. of WF: norm = ',
     &         sqrt(inprdd(civec1,civec2,lusc4,lusc4,1,lblk))
          print *,' calling compare routine:'
          call cmp2vcd(civec1,civec2,lusc4,lusc8,1d-10,1,lblk)
          print *,'==============================================='
        end if
c TEST
        if (ntest.ge.1000) then
          write (6,*) '1st derivative of wavefunction, element ',iamp
          call wrtvcd(civec1,lusc8,1,lblk)
        end if


        ! rewind file with old |a(ii)>
        call rewino(lusc7)

        ccvec2(1:iamp) = 0d0
        do jamp = 1, iamp-1
          call rewino(lusc8)
          sij = inprdd(civec1,civec2,lusc7,lusc8,0,lblk)
          if (ntest.ge.1000) write(6,*) iamp,jamp,': sij =',sij
          ccvec2(jamp) = sij
        end do
        sii = inprdd(civec1,civec2,lusc8,lusc8,1,lblk)
        if (ntest.ge.1000) write(6,*) iamp,iamp,': sii =',sii
        ccvec2(iamp) = sii
        call rewino(lusc8)
        ! append as last record
        call copvcd(lusc8,lusc7,civec1,0,lblk)
        call vec_to_disc(ccvec2,iamp,0,lblk,lufoo)

      end do ! loop over iamp

      call atim(cpu,wall)
      call prtim(6,'time in gtbce_foo',cpu-cpu0,wall-wall0)

      return
      end 
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
* DECK: gtbce_foo
*--------------------------------------------------------------------*
      subroutine gtbce_foo(igtb_closed_al,isymmet_G,irest,
c !!!!!!!!!!!!!!!!!!!!!!!!!^^^^^^^^^^^^^^!!!!!!!!!!!!!!!!!!!!!
     &                     inumint,npnts,
     &                     ovl,
     &                     iccvec,nsdim,
     &                     ccvec1,iopsym,comm_ops,
     &                     ccvec2,ccvec3,
     &                     civec1,civec2,c2vec,
     &                     n_cc_typ,i_cc_typ,ictp,
     &                     namp_cc_typ,ioff_cc_typ,
     &                     n_cc_amp,mxb_ci,
     &                     n11amp,n33amp,iamp_packed,
     &                     lufoo,
     &                     luamp,luc,luec,luhc,
     &                     lusc1,lusc2,lusc3,lusc4,
     &                     lusc5,lusc6,lusc7,lusc8,
     &                     lusc9,lusc10)
*--------------------------------------------------------------------*
*
* purpose: Calculate the overlap of the first order wavefunction
*          change
*
*          S_ij = N <0|(d/dg_i exp(G^+))(d/dg_j exp(G))|0>
*
*  ak, early 2004
*
*--------------------------------------------------------------------*
* diverse inludes with commons and paramters
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
      include 'gtbce.inc'
* debugging:
      integer, parameter :: ntest = 50
      logical, parameter :: tstrgad = .false.

* input/output arrays
      logical comm_ops
      integer iccvec(n_cc_amp)
      integer, intent(in) ::
     &     inumint, npnts
* scratch arrays
      character*8 cctype
      real*8 civec1(mxb_ci),civec2(mxb_ci),c2vec(*),
     &       ccvec1(n_cc_amp), ccvec2(n_cc_amp), ccvec3(n_cc_amp)
* local arrays
      real*8 alp(npnts+2), wght(npnts+2)
* external functions
      real*8 inprod, inprdd

      call atim(cpu0,wall0)

      iamp_rst = 0
      if (irest.gt.0) then
        iamp_rst = irest
      else
        if (isymmet_G.ne.0)
     &     iccvec(1:n_cc_amp) = 0
      end if

      lblk = -1
      if (ntest.ge.5) then
        write (6,*) '====================='
        write (6,*) ' This is gtbce_foo   '
        write (6,*) '====================='
        write (6,*) 
        write (6,*) 'on entry: '
        write (6,*) 'igtb_closed, isymmet_G: ',
     &       igtb_closed, isymmet_G
        write (6,*) 'inumint, npnts   : ', inumint, npnts
        write (6,*) 'iopsym: ',iopsym
        write (6,*) 'n_cc_amp,mxb_ci : ', n_cc_amp,mxb_ci
        write (6,*) 'luc,luec,luhc,lusc1,lusc2: ',
     &               luc,luec,luhc,lusc1,lusc2
        write (6,*) 'lusc3,lusc4,lusc5,lusc6,lusc7,lusc8,lusc9: ',
     &               lusc3,lusc4,lusc5,lusc6,lusc7,lusc8,lusc9
      end if
      if (ntest.ge.1000) then
        write(6,*) 'on entry:'
        write(6,*) 'e^G|0> on LUEC'
        call wrtvcd(civec1,luec,1,lblk)
      end if
      
      ! for I/O
      lblk = -1
      ! for expt_ref
      thresh=expg_thrsh
      mx_term=-mxterm_expg
      cctype='GEN_CC'
*--------------------------------------------------------------------*
* set up points and weights
*--------------------------------------------------------------------*
      select case (inumint)
      case (0)  ! just testing 
        do ipnt = 1, npnts
          alp(ipnt) = dble(ipnt-1)/dble(npnts-1)
          wght(ipnt) = 1d0
        end do
      case (1)  ! Gauss-Legendre
        call gl_weights(0d0,1d0,npnts,alp,wght)
      case (2)  ! Simpson
c        if (mod(npnts,2).eq.0) npnts = npnts-1
        call s_weights(0d0,1d0,npnts,alp,wght)
      case default
        stop 'unknown inumint in gtbce_foo'
      end select

      call vec_from_disc(ccvec1,n_cc_amp,1,-1,luamp)
      mxpnts=npnts
      ! if G == 0 ...
      xnrm2 = inprod(ccvec1,ccvec1,n_cc_amp)
      ! ... we have a set of commuting operators
      ! things are trivial and we evaluate the formula only once
      if (xnrm2.lt.10d-20.or.comm_ops) then
        mxpnts=1
        wght(1)=1d0
        alp(1)=0d0
        if (ntest.ge.5.and..not.comm_ops) then
          write(6,*) 'Detected zero amplitudes: ',
     &               'only case alpha = 0 will be processed'
        end if
      end if

      ! rewind output file
      call rewino(lufoo)

      if (iamp_rst.gt.0) then

        write(6,*) 'position unit ',lufoo,' after record ',iamp_rst
        call util_flush(6)
        call skpvcd(lufoo,iamp_rst,ccvec2,1,lblk)

      end if

      nsdim = 0
      mxperbatch = 50
      nbatch = n_cc_amp/mxperbatch
      if (mod(n_cc_amp,mxperbatch).gt.0) nbatch = nbatch+1
**-------------------------------------------------------------------*
* loop i over parameters in batches
**-------------------------------------------------------------------*
      do ibatch = 1, nbatch
        namp = mxperbatch
        if (ibatch.eq.nbatch) namp = n_cc_amp - mxperbatch*(nbatch-1)
        ista = (ibatch-1)*mxperbatch+1
        iend = (ibatch-1)*mxperbatch+namp

        if (iamp_rst.ne.0) then
          if (iend.lt.iamp_rst+1) then
            cycle
          else if (ista.le.iamp_rst+1) then
            ista = iamp_rst+1
            iamp_rst = 0
            write(6,*) 'restarting calculation from amplitude ',ista
          else
            write(6,*) 'error: ',ista,iend,iamp_rst
            stop 'impossible things happen sometimes ....'
          end if
        end if

        if (isymmet_G.ne.0) then
          iskip = 1
          do iamp = ista,iend
            if (iccvec(iamp).eq.0) then
              iskip = 0
              exit
            end if
          end do
          if (iskip.eq.1) cycle
        end if
c      do iamp = 1, n_cc_amp
        if (ntest.ge.10) write(6,*) 'batch, start, end ',
     &       ibatch,ista,iend

        ! reset |0tilde> = exp(G)|0>
        call copvcd(luec,lusc1,civec1,1,lblk)

**-------------------------------------------------------------------*
* loop over quadrature points
**-------------------------------------------------------------------*
        do ipnt = 1, mxpnts
          if (ntest.ge.50) then
            write(6,*) 'info for quadrature point: ', ipnt,'/',npnts
            write(6,*) 'point, weight: ', alp(ipnt), wght(ipnt)
          end if

          if (ipnt.gt.1.and.(alp(ipnt).le.alp(ipnt-1))) then
            write(6,*) 'quadrature point should be in ascending order!'
            stop 'gtbce_foo > quadrature '
          end if

          if (ipnt.eq.1) then
            dltalp = -alp(1)
          else
            dltalp = -alp(ipnt)+alp(ipnt-1)
          end if
*--------------------------------------------------------------------*
* |a_i>(1) = exp(-a_i G) [exp(G)|0>]
*          = exp(-(a_i-a_{i-1})G) exp(-a_{i-1} G) [exp(G)|0>]
*  result on lusc2
*--------------------------------------------------------------------*
          if (ntest.ge.50) then
            write(6,*)
     &           'constructing |a_i> = exp(-a_i G^+) exp(G)|0>]'
          end if
        
          if (abs(dltalp).lt.1d-20) then
            call copvcd(lusc1,lusc2,civec1,1,lblk)
          else
            ! get G on ccvec2
            call copvec(ccvec1,ccvec2,n_cc_amp)
            ! and scale it
            call scalve(ccvec2,dltalp,n_cc_amp)
            call expt_ref2(lusc1,lusc2,lusc4,lusc5,lusc6,
     &           thresh,mx_term, ccvec2, ccvec3, civec1, civec2,          
     &           n_cc_amp,cctype, iopsym)
            if (ntest.ge.100) then
              xnrm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
              etest = inprdd(civec1,civec2,luc,lusc2,1,lblk)
              write(6,*) '|dlta G^+|, dlta = ',xnrm, dltalp
              write(6,*) '<ref|a_i> = ', etest,
     &                   'for alp(i) = ', alp(ipnt) 
            end if
            ! save for next round
            call copvcd(lusc2,lusc1,civec1,1,lblk)
          end if

*--------------------------------------------------------------------*
* |a_i(ii)>(2) = tau_ii exp(-a_i G)[exp(G)|0>] for each paramter in batch
*  result on lusc3
*--------------------------------------------------------------------*
          call rewino(lusc8)
          call rewino(lusc9)
          
          ! alternate units to collect contributions
          ! the final result is on lunew in the end
          if (ipnt.eq.1) then
            lunew = lusc8
            luold = lusc9
          else
            if (lunew.eq.lusc8) then
              lunew = lusc9
              luold = lusc8
            else
              lunew = lusc8
              luold = lusc9
            end if
          end if

          if (ntest.ge.50) then
            write(6,*)
     &           'constructing |a_i> = exp((1-a_i) G^+)tau_i '//
     &           'exp(-a_i G^+) exp(G)|0>]'
          end if
          do iamp = ista, iend
            if (ntest.ge.50) then
              write (6,*) 'batch: ',ibatch,' iamp = ',iamp

              if (isymmet_G.ne.0) then
                if (ntest.ge.50)
     &               write(6,*) ' iccvec(iamp): ',
     &               iamp,iccvec(iamp)
                if (ntest.ge.50.and.iccvec(iamp).lt.0)
     &               write(6,*) ' this amplitude is skipped'
                if (iccvec(iamp).lt.0) cycle
              end if
            end if
            ccvec2(1:n_cc_amp) = 0d0
            ccvec2(iamp) = 1d0
            if (isymmet_G.ne.0) then
              ! (anti-)symmetrize
              call symmet_t(isymmet_G,1,
     &             ccvec2,ccvec3,
     &             ictp,i_cc_typ,n_cc_typ,
     &             namp_cc_typ,ioff_cc_typ,ngas)
              ! if not already marked, do that now:
              if (iccvec(iamp).eq.0) then
                ! remains non-vanishing amplitude afterwards?
                if (abs(inprod(ccvec2,ccvec2,n_cc_amp)).lt.1d-12) then
                  if (ntest.ge.100)
     &                write(6,*) ' aha, amplitude was diagonal! skipped'
                  iccvec(iamp) = -iamp
                  cycle
                end if
                if (ntest.ge.50)
     &               write(6,*) ' this amplitude is taken'

                if (abs(abs(ccvec2(iamp)-1d0)).lt.1d-12) then
                  iccvec(iamp) = iamp
                  if (ntest.ge.50) then
                    write(6,*) ' iamp, counterpart : ',iamp,iamp
                  end if
                else
                  ! mark counterpart as inactive              
                  do ii = iamp+1, n_cc_amp
                    if (abs(abs(ccvec2(ii))-0.5d0).lt.1d-12) then
                      if (ntest.ge.50) then
                        write(6,*) ' iamp, counterpart : ',iamp,ii
                      end if
                      iccvec(ii) = -iamp
                      iccvec(iamp) = ii
                      exit
                    end if
                  end do
                end if
                nsdim = nsdim + 1
              end if
            end if

            isigden=1
            if (iopsym.eq.0) then
              call sigden_cc(civec1,civec2,lusc2,lusc3,ccvec2,isigden)
            else
              call sigden_cc(civec1,civec2,lusc2,lusc4,ccvec2,isigden)
              fac = dble(iopsym)
              call conj_ccamp(ccvec2,1,ccvec3)
              call conj_t
              call sigden_cc(civec1,civec2,lusc2,lusc5,ccvec3,isigden)
              call conj_t
              call vecsmdp(civec1,civec2,1d0,fac,
     &                     lusc4,lusc5,lusc3,1,lblk)
            end if

*--------------------------------------------------------------------*
* |a_i(ii)>(3) = exp(a_i G) tau_ii exp(-a_i G)[exp(G)|0>]
*  result on lusc4
*--------------------------------------------------------------------*
            if (abs(alp(ipnt)).lt.1d-20) then
              call copvcd(lusc3,lusc4,civec1,1,lblk)
            else
              ! get G on ccvec2
              call copvec(ccvec1,ccvec2,n_cc_amp)
              ! and scale it
              call scalve(ccvec2,alp(ipnt),n_cc_amp)
              call expt_ref2(lusc3,lusc4,lusc5,lusc6,lusc10,
     &             thresh,mx_term, ccvec2, ccvec3, civec1, civec2,          
     &             n_cc_amp,cctype, iopsym)
              if (ntest.ge.100) then
                xnrm = sqrt(inprod(ccvec2,ccvec2,n_cc_amp))
                etest = inprdd(civec1,civec2,luc,lusc4,1,lblk)
                write(6,*) '|alp(i) G^+|, alp(i) = ',xnrm, alp(ipnt)
                write(6,*) '<ref|a_i> = ', etest,
     &                 'for alp(i) = ', alp(ipnt) 
              end if
            end if

            if (ntest.ge.2000) then
              write (6,*) 'contribution to 1st derivative of ',
     &                    'wavefunction, element ',iamp,alp(ipnt)
              call wrtvcd(civec1,lusc4,1,lblk)
            end if

            ! add lusc4 to luold giving lunew
            if (ipnt.eq.1) then
              call rewino(lusc4)
              call sclvcd(lusc4,lunew,wght(ipnt),civec1,0,lblk)
            else
              call rewino(lusc4)
              call vecsmdp(civec1,civec2,1d0,wght(ipnt),
     &             luold,lusc4,lunew,0,lblk) 
            end if
          end do ! loop over iamp

        end do ! loop over quadrature

        call rewino(lunew)
        do iamp = ista, iend
          ! rewind file with old |a(ii)>          
          if (isymmet_G.ne.0) then
            if (ntest.ge.100)
     &           write(6,*) ' iamp, iccvec(iamp): ',iamp,iccvec(iamp)
            if (iccvec(iamp).le.0) cycle
            if (ntest.ge.100)
     &           write(6,*) ' taken! '
          end if
          call rewino(lusc7)

c          call skpvcd(lunew,iamp-ista,civec1,1,lblk)
          call rewino(lusc3)

          call copvcd(lunew,lusc3,civec1,0,lblk)

          if (ntest.ge.1000) then
            write (6,*) '1st derivative of wavefunction, element ',iamp
            call wrtvcd(civec1,lusc3,1,lblk)
          end if

          ccvec2(1:iamp) = 0d0
          icnt = 0
          ! get the lunew/lusc7 contrib to Sij
          do jamp = 1, iamp-1
            if (isymmet_G.ne.0) then
              if (ntest.ge.100)
     &           write(6,*) ' jamp, iccvec(jamp): ',jamp,iccvec(jamp)
              if (iccvec(jamp).le.0) cycle
              if (ntest.ge.100)
     &             write(6,*) ' taken! '
            end if
            icnt = icnt + 1
            call rewino(lusc3)
            sij = inprdd(civec1,civec2,lusc3,lusc7,0,lblk)
            if (ntest.ge.100) write(6,*) iamp,jamp,': sij =',sij
            ccvec2(icnt) = sij
c            ccvec2(jamp) = sij
          end do
          ! get the lunew/lunew contrib to Sij
          sii = inprdd(civec1,civec2,lusc3,lusc3,1,lblk)
          if (ntest.ge.100) write(6,*) iamp,iamp,': sii =',sii
          icnt = icnt+1
          ccvec2(icnt) = sii
c          ccvec2(iamp) = sii
          ! append vector iamp a last record on lusc7
          call rewino(lusc3)
          call copvcd(lusc3,lusc7,civec1,0,lblk)

          call vec_to_disc(ccvec2,icnt,0,lblk,lufoo)
c          call vec_to_disc(ccvec2,iamp,0,lblk,lufoo)

        end do ! loop over iamp within batch

      end do ! loop over batches of iamp

      if (isymmet_G.eq.0) nsdim = n_cc_amp
      if (ntest.ge.50) then
        write(6,*) 'dimension: ',nsdim
      end if

      call atim(cpu,wall)
      call prtim(6,'time in gtbce_foo',cpu-cpu0,wall-wall0)

      return
      end 
*--------------------------------------------------------------------*
      subroutine mk_iccvec(isymmet_G,lufoo,irest,
     &                    iccvec,nSdim,ccvec1,ccvec2,
     &                    n_cc_typ,i_cc_typ,ictp,
     &                    namp_cc_typ,ioff_cc_typ,ngas,
     &                    n_cc_amp)
*--------------------------------------------------------------------*
*     set up iccvec array and nsdim for restarts
*--------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     lufoo,
     &     isymmet_G,n_cc_amp,ngas,n_cc_typ(*),i_cc_typ(*),ictp(*),
     &     namp_cc_typ(*),ioff_cc_typ(*)

      integer, intent(out) ::
     &     iccvec(n_cc_amp), nSdim

      integer, intent(inout) ::
     &     irest

      real(8), intent(inout) ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp)

      logical ::
     &     testrec

      integer ::
     &     iamp, ii, ierr

      real(8), external ::
     &     inprod

      nsdim = 0
      testrec = irest.ne.0
      if (irest.ne.0) call rewino(lufoo)
      iccvec(1:n_cc_amp) = 0
      do iamp = 1, n_cc_amp
        if (ntest.ge.100)
     &       write(6,*) ' iccvec(iamp): ',
     &       iamp,iccvec(iamp)
        if (iccvec(iamp).eq.0) then
          ccvec1(1:n_cc_amp) = 0d0
          ccvec1(iamp) = 1d0
          ! (anti-)symmetrize
          call symmet_t(isymmet_G,1,
     &         ccvec1,ccvec2,
     &         ictp,i_cc_typ,n_cc_typ,
     &         namp_cc_typ,ioff_cc_typ,ngas)
          ! remains non-vanishing amplitude afterwards?
          if (abs(inprod(ccvec1,ccvec1,n_cc_amp)).lt.1d-12) then
            if (ntest.ge.100)
     &           write(6,*) ' aha, amplitude was diagonal! skipped'
            iccvec(iamp) = -iamp
            cycle
          end if
          ! if requested, test whether this record is present on lufoo
          if (testrec) then
            call vec_from_disc_e(ccvec2,nsdim+1,0,-1,lufoo,ierr)
            if (ierr.eq.2) write(6,*) 'I/O-error detected :-('
            if (ierr.eq.1) write(6,*) 'EOF detected :-|'
            if (ierr.eq.0) write(6,*) 'record is fine :-)'
            if (ierr.ne.0) then
              irest = nsdim
              testrec = .false.
            else
              irest = nsdim+1
            end if
          end if

          if (abs(abs(ccvec1(iamp)-1d0)).lt.1d-12) then
            nsdim = nsdim + 1
            iccvec(iamp) = iamp
            if (ntest.ge.100) then
              write(6,*) ' iamp, counterpart : ',iamp,iamp
            end if
          else
            ! mark counterpart as inactive              
            nsdim = nsdim + 1
            do ii = iamp+1, n_cc_amp
              if (abs(abs(ccvec1(ii))-0.5d0).lt.1d-12) then
                if (ntest.ge.100) then
                  write(6,*) ' iamp, counterpart : ',iamp,ii
                end if
                iccvec(ii) = -iamp
                iccvec(iamp) = ii
                exit
              end if
            end do
          end if
        end if
                    
      end do

      if (ntest.ge.100) write(6,*) 'dimension of S: ',nSdim

      return
      
      end
*--------------------------------------------------------------------*
* DECK: gtbce_anahss
*--------------------------------------------------------------------*
      subroutine gtbce_anahss(hessi,luhss,ludia,istmode,
     &                        n_cc_amp,n_cc_typ,i_cc_typ,
     &                        namp_cc_typ,ioff_cc_typ,iopsym)
*--------------------------------------------------------------------*
*
* analyze a 2nd derivative matrix:
*  print blocks and get eigenvalues
*
*  istmode: 1 -- full matrix on file (one column per block)
*           2 -- upper triangle on file (one column up to diagonal
*                per block)
*--------------------------------------------------------------------*
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
* constants
      integer, parameter ::
     &     ntest = 100

* external functions
      real*8 inprod, inprdd

* input
      integer, intent(in) ::
     &     i_cc_typ(n_cc_typ), namp_cc_typ(n_cc_typ),
     &     ioff_cc_typ(n_cc_typ)
      real*8, intent(inout) ::
     &     hessi(n_cc_amp,n_cc_amp)

      lblk = -1

* read file luhss
      call rewino(luhss)
      hessi(1:n_cc_amp,1:n_cc_amp) = 0d0
      do iirec = 1, n_cc_amp
        if (ntest.ge.10) write (6,*) 'read rec. ',iirec
        nread = n_cc_amp
        if (istmode.eq.2) nread = iirec
        call vec_from_disc(hessi(1,iirec),nread,0,lblk,luhss)
      end do
      if (ntest.ge.100) then
        write(6,*) 'The Hessian as read in:'
        call wrtmat2(hessi,n_cc_amp,n_cc_amp,n_cc_amp,n_cc_amp)
      end if

      ! some waste of time, but for the moment much easier:
      ! get full matrix
      if (istmode.eq.2) then
        do ii = 1, n_cc_amp
          do jj = ii+1, n_cc_amp
            hessi(jj,ii)=hessi(ii,jj)
          end do
        end do
      else if (istmode.eq.3) then
        do ii = 1, n_cc_amp
          do jj = ii+1, n_cc_amp
            xel = 0.5d0*(hessi(jj,ii)+hessi(ii,jj))
            hessi(jj,ii)= xel
            hessi(ii,jj)= xel
          end do
        end do
      end if
      if (ntest.ge.100) then
        write(6,*) 'The Hessian as full matrix:'
        call wrtmat2(hessi,n_cc_amp,n_cc_amp,n_cc_amp,n_cc_amp)
      end if

* print-out of raw blocks
c      if (ntest.ge.5) then
c       do ii_tp = 1, n_cc_typ
c        iioff = ioff_cc_typ(ii_tp)
c        iilen = namp_cc_typ(ii_tp)
c        do jj_tp = 1, n_cc_typ
c          jjoff = ioff_cc_typ(jj_tp)
c          jjlen = namp_cc_typ(jj_tp)
c          write (6,*) 'block: ',ii_tp, jj_tp
c          call wrtmat(hessi(iioff,jjoff),iilen,jjlen,n_cc_amp,n_cc_amp)
c        end do
c       end do
c      end if

* diagonalize the matrix
      ltria = n_cc_amp*(n_cc_amp+1)/2
      leig  = n_cc_amp
      lscr  = 80*n_cc_amp
      idum = 0
      call memman(idum,idum,'MARK',idum,'TSTHSS')
      call memman(ktria,ltria,'ADDL',2,'HSSTRIA')
      call memman(keig,leig,'ADDL',2,'HSS EIG')
      call memman(kscr,lscr,'ADDL',2,'HSS SCR')
      
      call copdia(hessi,work(keig),n_cc_amp,0)
      write(6,*) 'the diagonal:'
      call wrtmat_ep(work(keig),n_cc_amp,1,n_cc_amp,1)

      irt = 1
      if (irt.eq.0) then
        iway = -1 ! symmetrize on the way
        call tripak(hessi,work(ktria),iway,n_cc_amp,n_cc_amp)
        call jacobi(work(ktria),hessi,n_cc_amp,n_cc_amp)
        call copdia(work(ktria),work(keig),n_cc_amp,1)
        stop 'test purpose route only'
      else if(irt.eq.1) then
        call diag_symmat_eispack(hessi,work(keig),work(ktria),
     &         n_cc_amp,iret)
        if (ntest.ge.100) then
          write(6,*) 'Eigenvector array:'
          call wrtmat2(hessi,n_cc_amp,n_cc_amp,n_cc_amp,n_cc_amp)
        end if
      else
        stop 'irt = ???'
      end if
c      hessi(1:n_cc_amp,1:n_cc_amp) = 0d0
c      do ii = 1, n_cc_amp
c        hessi(ii,ii) = 1d0
c      end do
c      work(keig:keig-1+leig) = 0d0
c      eps = 1d-14
c      call rdiag(work(ktria),hessi,work(keig),n_cc_amp,eps,work(kscr))
      
      write(6,*) 'the eigenvalues:'
      call wrtmat_ep(work(keig),n_cc_amp,1,n_cc_amp,1)

c      thrs = 1d-8
c      do ii = 1, n_cc_amp
c        if (work(keig-1+ii).gt.thrs) then
c          write(6,*) 'the eigenvector ',ii,work(keig-1+ii)
c          do ii_tp = 1, n_cc_typ
c            iioff = ioff_cc_typ(ii_tp)
c            iilen = namp_cc_typ(ii_tp)
c            xnrm = sqrt(inprod(hessi(iioff,ii),hessi(iioff,ii),iilen))
c            write (6,*) ' contributions from typ', ii_tp, xnrm
c            if (xnrm.gt.0.1*dble(iilen)) 
c     &       call wrtmat(hessi(iioff,ii),1,iilen,1,n_cc_amp)
c          end do
c        end if
c      end do
      
      imk_hinv = 0
      if (imk_hinv.eq.1) then
        ! find lowest eigenvalue and shift according to xdiag_min
        ! get column of hinv as
        !   hinv(i,j) = U(i,k) eig(k) U(j,k) 
      end if

      idum = 0
      call memman(idum,idum,'FLUSM',idum,'TSTHSS')

      return
      end
*--------------------------------------------------------------------*
*--------------------------------------------------------------------*
      subroutine gtbce_getrdvec(isymmet_G,
     &                        xsmat,lusmat,lurdvec,nrdvec,
     &                        nsmat,n_cc_amp,iccvec,
     &                        ccvec1,ccvec2)
*--------------------------------------------------------------------*
*
*  get redundant directions from smat
*  upper triangle on file (one column up to diagonal per block)
*--------------------------------------------------------------------*
c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
* constants
      integer, parameter ::
     &     ntest = 100

* external functions
      real*8 inprod, inprdd

* input
      integer, intent(in) ::
     &     iccvec(n_cc_amp)
      real*8, intent(inout) ::
     &     xsmat(nsmat,nsmat), ccvec1(n_cc_amp), ccvec2(n_cc_amp)

      lblk = -1

* read file luhss
      call rewino(lusmat)
      xsmat(1:nsmat,1:nsmat) = 0d0
      do iirec = 1, nsmat
        if (ntest.ge.10) write (6,*) 'read rec. ',iirec
        nread = iirec
c        if (istmode.eq.2) nread = iirec
        call vec_from_disc(xsmat(1,iirec),nread,0,lblk,lusmat)
      end do
      if (ntest.ge.100) then
        write(6,*) 'The S-matrix as read in:'
        call wrtmat2(xsmat,nsmat,nsmat,nsmat,nsmat)
      end if

      ! some waste of time, but for the moment much easier:
      ! get full matrix
      do ii = 1, nsmat
        do jj = ii+1, nsmat
          xsmat(jj,ii)=xsmat(ii,jj)
        end do
      end do
      if (ntest.ge.100) then
        write(6,*) 'The S-matrix as full matrix:'
        call wrtmat2(xsmat,nsmat,nsmat,nsmat,nsmat)
      end if

* diagonalize the matrix
      ltria = nsmat*(nsmat+1)/2
      leig  = nsmat
      lscr  = 80*nsmat
      idum = 0
      call memman(idum,idum,'MARK',idum,'TSTHSS')
      call memman(ktria,ltria,'ADDL',2,'HSSTRIA')
      call memman(keig,leig,'ADDL',2,'HSS EIG')
      call memman(kscr,lscr,'ADDL',2,'HSS SCR')
      
      call copdia(xsmat,work(keig),nsmat,0)
      write(6,*) 'the diagonal:'
      call wrtmat_ep(work(keig),nsmat,1,nsmat,1)

      irt = 1
      if (irt.eq.0) then
        iway = -1 ! symmetrize on the way
        call tripak(xsmat,work(ktria),iway,nsmat,nsmat)
        call jacobi(work(ktria),xsmat,nsmat,nsmat)
        call copdia(work(ktria),work(keig),nsmat,1)
        stop 'test purpose route only'
      else if(irt.eq.1) then
        call diag_symmat_eispack(xsmat,work(keig),work(ktria),
     &         nsmat,iret)
        if (ntest.ge.100) then
          write(6,*) 'Eigenvector array:'
          call wrtmat2(xsmat,nsmat,nsmat,nsmat,nsmat)
        end if
      else
        stop 'irt = ???'
      end if

      write(6,*) 'the eigenvalues:'
      call wrtmat_ep(work(keig),nsmat,1,nsmat,1)

      thrsh = 1d-12
c      thrsh = 1d-7
      nrdvec=0
      fac = dble(isymmet_G)
      call rewino(lurdvec)
      do ii = 1, nsmat
        if (work(keig-1+ii).lt.thrsh) then
          nrdvec = nrdvec+1
          ! expand this eigenvector to full aray
          ccvec1(1:n_cc_amp) = 0d0
          ismat = 0
          do iamp = 1, n_cc_amp
            if (iccvec(iamp).gt.0) then
              ismat = ismat+1
              if (ismat.gt.nsmat)
     &             stop 'inconsistency!'
              ccvec1(iamp)=xsmat(ismat,ii)
              if (isymmet_G.ne.0) then
                idx = iccvec(iamp)
                ccvec1(idx)=fac*xsmat(ismat,ii)
              end if
            end if
          end do
          ! renormalize and
          ! save as next record on lurdvec
          xnrm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))
          ccvec1(1:n_cc_amp) = 1d0/xnrm*ccvec1(1:n_cc_amp)
          call vec_to_disc(ccvec1,n_cc_amp,0,-1,lurdvec)
        end if
      end do

      write(6,*) '>> # redundant vectors:     ',nrdvec

      idum = 0
      call memman(idum,idum,'FLUSM',idum,'TSTHSS')

      return
      end
*--------------------------------------------------------------------*
      subroutine gtbce_prjout_rdvec(nrdvec,lurdvec,luvec,
     &     n_cc_amp,ccvec1,ccvec2)

      implicit none

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     nrdvec, lurdvec, luvec, n_cc_amp
      real(8), intent(inout) ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp)
      
      integer ::
     &     irdvec
      real(8) ::
     &     ovl, xnrm
      real(8), external ::
     &     inprod

      call vec_from_disc(ccvec1,n_cc_amp,1,-1,luvec)

      xnrm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))

      write(6,*) ' norm of unprojected gradient: ',xnrm

      call rewino(lurdvec)
      do irdvec = 1, nrdvec
        call vec_from_disc(ccvec2,n_cc_amp,0,-1,lurdvec)
        ovl = inprod(ccvec1,ccvec2,n_cc_amp)
        write(6,*) ' overlap with vec ',irdvec,' :',ovl
        ccvec1(1:n_cc_amp) = ccvec1(1:n_cc_amp) - ovl*ccvec2(1:n_cc_amp)
      end do

      xnrm = sqrt(inprod(ccvec1,ccvec1,n_cc_amp))

      write(6,*) ' norm of projected gradient:   ',xnrm

      call vec_to_disc(ccvec1,n_cc_amp,1,-1,luvec)

      return
      end
*--------------------------------------------------------------------*
      subroutine gtbce_EalongG(tvec,npnts,from_g,to_g,
     &               ecore,
     &               ccvec1,iopsym,ccvec3,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)
*--------------------------------------------------------------------*
*
* purpose: calculate energy along a selected direction tvec and 
*          generate plot data
*
*  ak, early 2004
*
*--------------------------------------------------------------------*
      include "implicit.inc"

* input
      real*8, intent(in) ::
     &     ccvec1(n_cc_amp), tvec(n_cc_amp)
      real*8, intent(inout) ::
     &     ccvec3(n_cc_amp), ccvec4(n_cc_amp)
* external
      real*8 ::
     &     inprod


      xdelt = to_g - from_g
      
      xinc = xdelt/dble(npnts-1)
      xnorm = sqrt(inprod(tvec,  tvec,n_cc_amp))
      xovl  =      inprod(ccvec1,tvec,n_cc_amp)

      write (6,'("@p",a,e10.4)') ' comp. of G along t: ',xovl/xnorm
      
      write (6,'("@p",a)') ' n    c    energy   variance   dnorm'
      
      do ipnt = 0, npnts-1
        
        fac = (from_g+xdelt*dble(ipnt)/dble(npnts-1))/xnorm
        ccvec3(1:n_cc_amp) = ccvec1(1:n_cc_amp)+fac*tvec(1:n_cc_amp)
        
c        igtbmod = 1
        call gtbce_E(igtbmod,elen,var,ovl,
     &               ecore,
     &               ccvec3,iopsym,ccvec4,
     &               civec1,civec2,c2vec,
     &               n_cc_amp,mxb_ci,
     &               luc,luec,luhc,lusc1,lusc2)

        write(6,'("@p",i4,e14.6,f21.12,2(2x,e10.4))')
     &           ipnt,fac,elen,var,1d0-sqrt(ovl)

      end do

      return

      end

*--------------------------------------------------------------------*
      subroutine prjout_red(gop,xrs,ntss_tp,itss_tp,ibtss_tp)
*--------------------------------------------------------------------*
*
*     project out redundant directions from input vector
*
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 1000

* input
      real(8), intent(inout) ::
     &     gop(*), xrs(*)

c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp),
     &     ibtss_tp(ntss_tp)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4),
     &     irs(ntoob*ntoob)

      if (ntest.ge.1000) then
        write(6,*) ' input amplitudes: '
        call wrt_cc_vec2(gop,6,'GEN_CC')
        write(6,*) 'ibtss_tp:'
        call iwrtma(ibtss_tp,1,ntss_tp,1,ntss_tp)
      end if
        
      ! init
      xrs(1:ntoob*ntoob) = 0d0
      irs(1:ntoob*ntoob) = 0
      do ipass = 1, 2
        !
        ! run over all operator elements and ...
        !
        ! pass 1:
        !  X^{(rs)} = sum_p G_pprs (p,p of equal spin)
        !
        ! pass 2:
        ! subtract X^{(rs)} from each entry G_pprs

        ! loop over types
        idx = 0
        do itss = 1, ntss_tp
          ! identify two-particle excitations:
          nel_ca = ielsum(itss_tp(1,1,itss),ngas)
          nel_cb = ielsum(itss_tp(1,2,itss),ngas)
          nel_aa = ielsum(itss_tp(1,3,itss),ngas)
          nel_ab = ielsum(itss_tp(1,4,itss),ngas)
          nc = nel_ca + nel_cb
          na = nel_aa + nel_ab
          if (na.ne.2) cycle
          ! transform occupations to groups
          call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
          call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
          call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
          call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)
          
          if (mscomb_cc.ne.0) then
            call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &           itss_tp(1,3,itss),itss_tp(1,4,itss),
     &           ngas,idiag)
          else
            idiag = 0
          end if
        

          ! loop over symmetry blocks
          ism = 1 ! totally symmetric operators
          do ism_c = 1, nsmst
            ism_a = multd2h(ism,ism_c)
            do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      ! have one particle and one hole operator the same index?
c TEST tabula rasa test
c                      if ((idx_s(1).eq.idx_s(2).and.
                      if ((!idx_s(1).eq.idx_s(2).and.
     &                      (idx_c(1).eq.idx_c(3).or.
     &                       idx_c(1).eq.idx_c(4).or.
     &                       idx_c(2).eq.idx_c(3).or.
     &                       idx_c(2).eq.idx_c(4))) .or.
     &                    (idx_s(1).ne.idx_s(2).and.
     &                      (idx_c(1).eq.idx_c(3).or.
     &                       idx_c(2).eq.idx_c(4)) ) ) then
                        if (idx_c(1).eq.idx_c(3))
     &                       idx_rs = (idx_c(2)-1)*ntoob
     &                               + idx_c(4)
                        if (idx_c(1).eq.idx_c(4))
     &                       idx_rs = (idx_c(2)-1)*ntoob
     &                               + idx_c(3)
                        if (idx_c(2).eq.idx_c(3))
     &                       idx_rs = (idx_c(1)-1)*ntoob
     &                               + idx_c(4)
                        if (idx_c(2).eq.idx_c(4))
     &                       idx_rs = (idx_c(1)-1)*ntoob
     &                               + idx_c(3)
                        if (ipass.eq.1) then
                          xrs(idx_rs) = xrs(idx_rs) + gop(idx)
                          irs(idx_rs) = irs(idx_rs) + 1
                        end if
c                        if (ipass.eq.2)
c     &                     gop(idx) = gop(idx)
c     &                       - 1d0/dble(irs(idx_rs))*xrs(idx_rs)
c TEST --- tabula rasa for all amplitudes with repeated indices
                        if (ipass.eq.2)
     &                       gop(idx) = 0d0


                      end if

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
            end do ! ism_ca
          end do ! ism_c 

        end do ! itss

        if (ipass.eq.1.and.ntest.ge.150) then
          write(6,*) 'The xrs array:'
          call wrtmat(xrs,ntoob,ntoob,ntoob,ntoob)
          write(6,*) 'The irs array:'
          call iwrtma(irs,ntoob,ntoob,ntoob,ntoob)
        end if

      end do ! ipass

      if (ntest.ge.1000) then
        write(6,*) ' output amplitudes: '
        call wrt_cc_vec2(gop,6,'GEN_CC')
      end if

      return
      end

*--------------------------------------------------------------------*
      subroutine ggrad2lgrad(ggrad,lgrad,lop,
     &     ntss_tp,itss_tp,nloff,ldiml)
*--------------------------------------------------------------------*
*
*
*
*--------------------------------------------------------------------*

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 100

* input
      real(8), intent(in) ::
     &     ggrad(*), lop(*)
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp)

      real(8), intent(out) ::
     &     lgrad(*)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)

      ! init
      lgrad(1:ldiml**2) = 0d0
      if (nloff.gt.0) lgrad(nloff:nloff+ldiml**2-1) = 0d0

      ! loop over types
      idx = 0
      do itss = 1, ntss_tp
        ! identify two-particle excitations:
        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab
        if (na.ne.2) cycle

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      idxpq = idx_s(1)*nloff +
     &                     (idx_c(3)-1)*ldiml + idx_c(1)
                      if (idx_s(1).ne.idx_s(3)) stop 'ups (1)'
                      idxrs = idx_s(2)*nloff +
     &                     (idx_c(4)-1)*ldiml + idx_c(2)
                      if (idx_s(2).ne.idx_s(4)) stop 'ups (2)'

c                      print *,'pq = ',idx_c(3), idx_c(1), idx_s(1)
c                      print *,'rs = ',idx_c(4), idx_c(2), idx_s(2)

c                      print '(x,a,3i4,a,3i4)',
c     &                     ' contr. ',idx_c(3),idx_c(1),idx_s(1),
c     &                         ' to ',idx_c(4),idx_c(2),idx_s(2)
                      lgrad(idxpq) = lgrad(idxpq)+ggrad(idx)*lop(idxrs)
c                      print '(x,a,3i4,a,3i4)',
c     &                     ' contr. ',idx_c(4),idx_c(2),idx_s(2),
c     &                         ' to ',idx_c(3),idx_c(1),idx_s(1)
c                      print *,' grad(',idxpq,idxrs,')=',ggrad(idx)
                      lgrad(idxrs) = lgrad(idxrs)+ggrad(idx)*lop(idxpq)

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss

      if (ntest.ge.100) then
        write(6,*) 'L gradient:'
        do ii = 1, ntoob
          do jj = 1, ntoob
            idx = (ii-1)*ntoob+jj
            write(6,*) ii,jj,lgrad(idx)
          end do
        end do
      end if

      return

      end

*--------------------------------------------------------------------*
      subroutine ggrad2omgrad(ggrad,omgrad,omop,urop,uiop,
     &     ntss_tp,itss_tp,ndim)
*--------------------------------------------------------------------*
*
*     get Omega gradient acc. to chain rule
*
*--------------------------------------------------------------------*

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 00

* input
      real(8), intent(in) ::
     &     ggrad(*),  omop(ndim,ndim,2,2),
     &     urop(ndim,ndim,2,2),  uiop(ndim,ndim,2,2)
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp)

      real(8), intent(out) ::
     &     omgrad(ndim,ndim,2,2)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)

      call atim(cpu0,wall0)

      ! init
      omgrad(1:ndim,1:ndim,1:2,1:2) = 0d0

      ! loop over types
      idx = 0
      do itss = 1, ntss_tp
        ! identify two-particle excitations:
        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab
        if (na.ne.2) cycle

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      ip = idx_c(1)
                      ir = idx_c(2)
                      iq = idx_c(3)
                      is = idx_c(4)

                      imp = idx_s(1)
                      imr = idx_s(2)
                      imq = idx_s(3)
                      ims = idx_s(4)

                      do imt = 1, 2
                        do imu = 1, 2
                          do it = 1, ndim
                            urur_t = urop(ip,it,imp,imt)*
     &                               urop(iq,it,imq,imt)
                            uiui_t = uiop(ip,it,imp,imt)*
     &                               uiop(iq,it,imq,imt)
                            uiur_t = uiop(ip,it,imp,imt)*
     &                               urop(iq,it,imq,imt)
                            urui_t = urop(ip,it,imp,imt)*
     &                               uiop(iq,it,imq,imt)


                            do iu = 1, ndim

                              urur_u = urop(ir,iu,imr,imu)*
     &                                 urop(is,iu,ims,imu)
                              uiui_u = uiop(ir,iu,imr,imu)*
     &                                 uiop(is,iu,ims,imu)
                              uiur_u = uiop(ir,iu,imr,imu)*
     &                                 urop(is,iu,ims,imu)
                              urui_u = urop(ir,iu,imr,imu)*
     &                                 uiop(is,iu,ims,imu)

                              omgrad(it,iu,imt,imu) =
     &                             omgrad(it,iu,imt,imu) +
     &                             ggrad(idx)
     &                             *((urur_t+uiui_t)*(uiur_u-urui_u)
     &                              +(uiur_t-urui_t)*(urur_u+uiui_u))
* new: update also the inversed pair
                              if (imu.ne.imt) then
                                omgrad(iu,it,imu,imt) =
     &                             omgrad(iu,it,imu,imt) +
     &                             ggrad(idx)
     &                             *((urur_t+uiui_t)*(uiur_u-urui_u)
     &                              +(uiur_t-urui_t)*(urur_u+uiui_u))
                              end if
                            
                            end do
                          end do
                        end do
                      end do

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss

      if (ntest.ge.100) then
        write(6,*) 'Omega gradient:'
        do imp = 1, 2
          do imq = 1, 2
            write(6,*) 'spin block: ',imp,imq
            call wrtmat2(omgrad(1,1,imp,imq),ndim,ndim,ndim,ndim)
          end do
        end do
      end if

      call atim(cpu,wall)
      call prtim(6,'time in ggrad2omgrad',cpu-cpu0,wall-wall0)

      return

      end

*--------------------------------------------------------------------*
      subroutine ggrad2ugrad(ggrad,urgrad,omop,urop,uiop,
     &     ntss_tp,itss_tp,ndim,irmod)
*--------------------------------------------------------------------*
*
*     get U gradient acc. to chain rule
*
*--------------------------------------------------------------------*

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 00

* input
      real(8), intent(in) ::
     &     ggrad(*),  omop(ndim,ndim,2,2),
     &     urop(ndim,ndim,2,2),  uiop(ndim,ndim,2,2)
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp)
      real(8), intent(out) ::
     &     urgrad(ndim,ndim,2,2)
* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)

      call atim(cpu0,wall0)

      ! init
      urgrad(1:ndim,1:ndim,1:2,1:2) = 0d0

      ! loop over types
      idx = 0
      do itss = 1, ntss_tp
        ! identify two-particle excitations:
        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab
        if (na.ne.2) cycle

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      idxpq = idx_s(1)*nloff +
     &                     (idx_c(3)-1)*ldiml + idx_c(1)
                      if (idx_s(1).ne.idx_s(3)) stop 'ups (1)'
                      idxrs = idx_s(2)*nloff +
     &                     (idx_c(4)-1)*ldiml + idx_c(2)
                      if (idx_s(2).ne.idx_s(4)) stop 'ups (2)'

                      ip = idx_c(1)
                      ir = idx_c(2)
                      iq = idx_c(3)
                      is = idx_c(4)

                      imp = idx_s(1)
                      imr = idx_s(2)
                      imq = idx_s(3)
                      ims = idx_s(4)

                      if (irmod.eq.1) fac =  1d0
c                      if (irmod.eq.2) fac = -1d0
                      if (irmod.eq.2) fac = 1d0


                      do imu = 1,2
                        do imw = 1, 2
                          do iu = 1, ndim
                            urur_rs = urop(ir,iu,imr,imu)*
     &                                  urop(is,iu,ims,imu)
                            uiui_rs = uiop(ir,iu,imr,imu)*
     &                                  uiop(is,iu,ims,imu)
                            uiur_rs = uiop(ir,iu,imr,imu)*
     &                                  urop(is,iu,ims,imu)
                            urui_rs = urop(ir,iu,imr,imu)*
     &                                  uiop(is,iu,ims,imu)

                            urur_pq = urop(ip,iu,imp,imu)*
     &                                  urop(iq,iu,imq,imu)
                            uiui_pq = uiop(ip,iu,imp,imu)*
     &                                  uiop(iq,iu,imq,imu)
                            uiur_pq = uiop(ip,iu,imp,imu)*
     &                                  urop(iq,iu,imq,imu)
                            urui_pq = urop(ip,iu,imp,imu)*
     &                                  uiop(iq,iu,imq,imu)


                            do iw = 1,ndim
* term 1
                              ur = urop(iq,iw,imq,imw)
                              urgrad(ip,iw,imp,imw) =
     &                                 urgrad(ip,iw,imp,imw) +
     &                                 ur * (uiur_rs-urui_rs) *
     &                                 omop(iw,iu,imw,imu) * ggrad(idx)

                              ur = urop(ip,iw,imp,imw)
                              urgrad(iq,iw,imq,imw) =
     &                                 urgrad(iq,iw,imq,imw) +
     &                                 ur * (uiur_rs-urui_rs) *
     &                                 omop(iw,iu,imw,imu) * ggrad(idx)

* term 2
                              
                              ur = urop(is,iw,ims,imw)
                              urgrad(ir,iw,imr,imw) =
     &                             urgrad(ir,iw,imr,imw) +
     &                             ur * (uiur_pq-urui_pq) *
     &                             omop(iu,iw,imu,imw)* ggrad(idx)

                              ur = urop(ir,iw,imr,imw)
                              urgrad(is,iw,ims,imw) =
     &                             urgrad(is,iw,ims,imw) +
     &                             ur * (uiur_pq-urui_pq) *
     &                             omop(iu,iw,imu,imw)* ggrad(idx)
                              
* term 3
                              ui = uiop(ir,iw,imr,imw)
                              urgrad(is,iw,ims,imw) =
     &                             urgrad(is,iw,ims,imw) + 
     &                             (urur_pq+uiui_pq) * ui * 
     &                             omop(iu,iw,imu,imw)* ggrad(idx)

                              ui = - uiop(is,iw,ims,imw)
                              urgrad(ir,iw,imr,imw) =
     &                             urgrad(ir,iw,imr,imw) +
     &                             (urur_pq+uiui_pq) * ui * 
     &                             omop(iu,iw,imu,imw)* ggrad(idx)

* term 4
                              ui = uiop(ip,iw,imp,imw)
                              urgrad(iq,iw,imq,imw) =
     &                             urgrad(iq,iw,imq,imw) + fac*
     &                             (urur_rs+uiui_rs) * ui * 
     &                             omop(iw,iu,imw,imu)* ggrad(idx)

                              ui = - uiop(iq,iw,imq,imw)
                              urgrad(ip,iw,imp,imw) =
     &                             urgrad(ip,iw,imp,imw) + fac*
     &                             (urur_rs+uiui_rs) * ui * 
     &                             omop(iw,iu,imw,imu)* ggrad(idx)

                            end do
                          end do
                        end do
                      end do

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss

c scale with 2d0
      fac = 2d0
      if (irmod.eq.2) fac = -2d0
      call scalve(urgrad,fac,4*ndim**2)
      

      if (ntest.ge.100) then
        write(6,*) 'U gradient:'
        do imp = 1, 2
          do imq = 1, 2
            write(6,*) 'spin block: ',imp,imq
            call wrtmat2(urgrad(1,1,imp,imq),ndim,ndim,ndim,ndim)
          end do
        end do
      end if

      call atim(cpu,wall)
      call prtim(6,'time in ggrad2ugrad',cpu-cpu0,wall-wall0)

      return

      end

*--------------------------------------------------------------------*
      subroutine ggrad2ugrad_old(ggrad,urgrad,omop,urop,uiop,
     &     ntss_tp,itss_tp,ndim,irmod)
*--------------------------------------------------------------------*
*
*     get U gradient acc. to chain rule
*
*--------------------------------------------------------------------*

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 100

* input
      real(8), intent(in) ::
     &     ggrad(*),  omop(ndim,ndim,2,2),
     &     urop(ndim,ndim,2,2),  uiop(ndim,ndim,2,2)
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp)
      real(8), intent(out) ::
     &     urgrad(ndim,ndim,2,2)
* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)

      call atim(cpu0,wall0)

      ! init
      urgrad(1:ndim,1:ndim,1:2,1:2) = 0d0

      ! loop over types
      idx = 0
      do itss = 1, ntss_tp
        ! identify two-particle excitations:
        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab
        if (na.ne.2) cycle

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      idxpq = idx_s(1)*nloff +
     &                     (idx_c(3)-1)*ldiml + idx_c(1)
                      if (idx_s(1).ne.idx_s(3)) stop 'ups (1)'
                      idxrs = idx_s(2)*nloff +
     &                     (idx_c(4)-1)*ldiml + idx_c(2)
                      if (idx_s(2).ne.idx_s(4)) stop 'ups (2)'

                      ip = idx_c(1)
                      ir = idx_c(2)
                      iq = idx_c(3)
                      is = idx_c(4)

                      imp = idx_s(1)
                      imr = idx_s(2)
                      imq = idx_s(3)
                      ims = idx_s(4)

                      if (irmod.eq.1) fac =  1d0
c                      if (irmod.eq.2) fac = -1d0
                      if (irmod.eq.2) fac = 1d0


                      do imu = 1,2
                        do imv = 1, 2
                          do imw = 1, 2
                            do iu = 1, ndim
                              urur_rs = urop(ir,iu,imr,imu)*
     &                                  urop(is,iu,ims,imu)
                              uiui_rs = uiop(ir,iu,imr,imu)*
     &                                  uiop(is,iu,ims,imu)
                              uiur_rs = uiop(ir,iu,imr,imu)*
     &                                  urop(is,iu,ims,imu)
                              urui_rs = urop(ir,iu,imr,imu)*
     &                                  uiop(is,iu,ims,imu)

                              urur_pq = urop(ip,iu,imp,imu)*
     &                                  urop(iq,iu,imq,imu)
                              uiui_pq = uiop(ip,iu,imp,imu)*
     &                                  uiop(iq,iu,imq,imu)
                              uiur_pq = uiop(ip,iu,imp,imu)*
     &                                  urop(iq,iu,imq,imu)
                              urui_pq = urop(ip,iu,imp,imu)*
     &                                  uiop(iq,iu,imq,imu)


                              do iv = 1, ndim
                                do iw = 1,ndim
* term 1
                                  ur = 0d0
                                  if (iv.eq.ip.and.imv.eq.imp) then
                                    ur = urop(iq,iw,imq,imw)
                                  end if
                                  if (iv.eq.iq.and.imv.eq.imq) then
                                    ur = ur + urop(ip,iw,imp,imw)
                                  end if
                                  urgrad(iv,iw,imv,imw) =
     &                                 urgrad(iv,iw,imv,imw) +
     &                                 ur * (uiur_rs-urui_rs) *
     &                                 omop(iw,iu,imw,imu) * ggrad(idx)
* term 2
                                  ur = 0d0
                                  if (iv.eq.ir.and.imv.eq.imr) then
                                    ur = urop(is,iw,ims,imw)
                                  end if
                                  if (iv.eq.is.and.imv.eq.ims) then
                                    ur = ur + urop(ir,iw,imr,imw)
                                  end if
                                  urgrad(iv,iw,imv,imw) =
     &                                 urgrad(iv,iw,imv,imw) +
     &                                 ur * (uiur_pq-urui_pq) *
     &                                 omop(iu,iw,imu,imw)* ggrad(idx)
* term 3
                                  ui = 0d0
                                  if (iv.eq.is.and.imv.eq.ims) then
                                    ui = uiop(ir,iw,imr,imw)
                                  end if
                                  if (iv.eq.ir.and.imv.eq.imr) then
                                    ui = ui - uiop(is,iw,ims,imw)
                                  end if
                                  urgrad(iv,iw,imv,imw) =
     &                                 urgrad(iv,iw,imv,imw) + fac*
     &                                 (urur_pq+uiui_pq) * ui * 
     &                                 omop(iu,iw,imu,imw)* ggrad(idx)
* term 4
                                  ui = 0d0
                                  if (iv.eq.iq.and.imv.eq.imq) then
                                    ui = uiop(ip,iw,imp,imw)
                                  end if
                                  if (iv.eq.ip.and.imv.eq.imp) then
                                    ui = ui - uiop(iq,iw,imq,imw)
                                  end if
                                  urgrad(iv,iw,imv,imw) =
     &                                 urgrad(iv,iw,imv,imw) + fac*
     &                                 (urur_rs+uiui_rs) * ui * 
     &                                 omop(iw,iu,imw,imu)* ggrad(idx)

                                end do
                              end do
                            end do
                          end do
                        end do
                      end do

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss

c scale with 2d0
      fac = 2d0
      if (irmod.eq.2) fac = -2d0
      call scalve(urgrad,fac,4*ndim**2)
      

      if (ntest.ge.100) then
        write(6,*) 'U gradient:'
        do imp = 1, 2
          do imq = 1, 2
            write(6,*) 'spin block: ',imp,imq
            call wrtmat2(urgrad(1,1,imp,imq),ndim,ndim,ndim,ndim)
          end do
        end do
      end if

      call atim(cpu,wall)
      call prtim(6,'time in ggrad2ugrad',cpu-cpu0,wall-wall0)

      return

      end

*--------------------------------------------------------------------*
      subroutine uou2g(omop,urop,uiop,gop,
     &                 ntss_tp,itss_tp,ibtss_tp,ndim)
*--------------------------------------------------------------------*
*
*     Set up elements of two-particle operator G according to
*
*     G(pq,rs)a_pqrs =  (....) a_pqrs
*
*--------------------------------------------------------------------*

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 000

* input
      real(8), intent(inout) ::
     &     omop(ndim,ndim,2,2),
     &     urop(ndim,ndim,2,2),
     &     uiop(ndim,ndim,2,2)
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp),
     &     ibtss_tp(ntss_tp)

      real(8), intent(out) ::
     &     gop(*)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)

      if (ntest.eq.1000) then
        write(6,*) '======'
        write(6,*) 'Omega:'
        write(6,*) '======'
        do imp = 1, 2
          do imq = 1, 2
            write(6,*) 'spin block: ',imp,imq
            call wrtmat2(omop(1,1,imp,imq),ndim,ndim,ndim,ndim)
          end do
        end do
        write(6,*) '======'
        write(6,*) 'U(Re):'
        write(6,*) '======'
        do imp = 1, 2
          do imq = 1, 2
            write(6,*) 'spin block: ',imp,imq
            call wrtmat2(urop(1,1,imp,imq),ndim,ndim,ndim,ndim)
          end do
        end do
        write(6,*) '======'
        write(6,*) 'U(Im):'
        write(6,*) '======'
        do imp = 1, 2
          do imq = 1, 2
            write(6,*) 'spin block: ',imp,imq
            call wrtmat2(uiop(1,1,imp,imq),ndim,ndim,ndim,ndim)
          end do
        end do

      end if

      ! loop over types
      idx = 0
      do itss = 1, ntss_tp
        if (ibtss_tp(itss).ne.idx+1) then
          write(6,*) 'problem with offset for op. ',itss
          write(6,*) '  ',ibtss_tp(itss),' != ',idx+1
        end if
        ! identify two-particle excitations:
        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab
        if (na.ne.2) cycle

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      ip = idx_c(1)
                      ir = idx_c(2)
                      iq = idx_c(3)
                      is = idx_c(4)

                      imp = idx_s(1)
                      imr = idx_s(2)
                      imq = idx_s(3)
                      ims = idx_s(4)

                      gop(idx) = 0d0

                      do imt = 1, 2
                        do imu = 1, 2
                          do it = 1, ndim
                            urur_t = urop(ip,it,imp,imt)*
     &                               urop(iq,it,imq,imt)
                            uiui_t = uiop(ip,it,imp,imt)*
     &                               uiop(iq,it,imq,imt)
                            uiur_t = uiop(ip,it,imp,imt)*
     &                               urop(iq,it,imq,imt)
                            urui_t = urop(ip,it,imp,imt)*
     &                               uiop(iq,it,imq,imt)

                            do iu = 1, ndim

                              urur_u = urop(ir,iu,imr,imu)*
     &                                 urop(is,iu,ims,imu)
                              uiui_u = uiop(ir,iu,imr,imu)*
     &                                 uiop(is,iu,ims,imu)
                              uiur_u = uiop(ir,iu,imr,imu)*
     &                                 urop(is,iu,ims,imu)
                              urui_u = urop(ir,iu,imr,imu)*
     &                                 uiop(is,iu,ims,imu)

                              gop(idx) = gop(idx) +
     &                             ((urur_t+uiui_t)*(uiur_u-urui_u)
     &                             +(uiur_t-urui_t)*(urur_u+uiui_u))
     &                             *omop(it,iu,imt,imu)

                              if (imt.ne.imu) then
                                gop(idx) = gop(idx) +
     &                               ((urur_t+uiui_t)*(uiur_u-urui_u)
     &                               +(uiur_t-urui_t)*(urur_u+uiui_u))
     &                               *omop(iu,it,imu,imt)
                              end if


                            end do
                          end do
                        end do
                      end do

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss
      if (ntest.ge.1000) then
        write(6,*) 'the two-particle operator:'
        call wrt_cc_vec2(gop,6,'GEN_CC')
      end if


      return

      end


*------------------------------------------------------------------------*
*     another clone of EXPT_REF:
*------------------------------------------------------------------------*
      SUBROUTINE EXPT2_REF(LUC,LUHC,LUSC1,LUSC2,LUSC3,
     &                    THRES_C,MX_TERM,
     &                    ALPHA,TAMP,TSCR,VEC1,VEC2,N_CC_AMP,
     &                    IOPTYP)
*
* Obtain Exp (alpha T^2) !ref> by Taylor expansion of exponential
*
* Orig. Version: Jeppe Olsen, March 1998 
*
* Extended to include general CC, summer of 99
*
* IOPTYP defines symmetry of operator: 
*
*    +1 Hermitian
*    -1 unitary
*     0 general
*
* TSCR is only needed in the first two cases.
*
c      IMPLICIT REAL*8(A-H,O-Z)
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'

      REAL*8 INPRDD, INPROD
*
      INCLUDE 'glbbas.inc'
      INCLUDE 'cprnt.inc'
*
      DIMENSION VEC1(*),VEC2(*),TAMP(*),TSCR(*)
      COMMON/CINT_CC/INT_CC
* 
      LBLK = -1
*
      NTEST = 5
      NTEST = MAX(NTEST,IPRCC)
*
      IF (IOPTYP.EQ.1) THEN
        SFAC = 1d0
      ELSE IF(IOPTYP.EQ.-1) THEN
        SFAC = -1d0
      ELSE IF (IOPTYP.NE.0) THEN
        WRITE(6,*) 'Indigestible input in EXPT_REF2!!!'
        STOP 'IOPTYP in EXPT_REF2'
      END IF
*
      IF(NTEST.GE.5) THEN
       WRITE(6,*)
       WRITE(6,*) '===================='
       WRITE(6,*) 'EXPT2_REF in action '
       WRITE(6,*) '===================='
       WRITE(6,*) ' ioptyp  = ',ioptyp
       WRITE(6,*) ' alpha   = ',alpha
       WRITE(6,*) ' mx_term = ',mx_term
       WRITE(6,*) ' thresh  = ',THRES_C
       WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' LUC,LUHC,LUSC1,LUSC2',LUC,LUHC,LUSC1,LUSC2
       WRITE(6,*) ' Initial vector on LUC '
       IF (NTEST.GE.1000) THEN
         CALL WRTVCD(VEC1,LUC,1,LBLK)
       ELSE
         CALL WRTVSD(VEC1,LUC,1,LBLK)
       END IF
      END IF
* Tell integral fetcher to fetch cc amplitudes, not integrals
      INT_CC = 1
*. Loop over orders of expansion
      N = 0
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) 'TAMP:'
        CALL WRT_CC_VEC2(TAMP,6,'GEN_CC')
      END IF

      IF (IOPTYP.NE.0) THEN
        CALL CONJ_CCAMP(TAMP,1,TSCR)
        CALL SCALVE(TSCR,SFAC,N_CC_AMP)
        IF(NTEST.GE.500) THEN
          WRITE(6,*) 'TAMP+:'
          CALL WRT_CC_VEC2(TSCR,6,'GEN_CC')
        END IF        
      END IF
*
      CALL COPVCD(LUC,LUSC1,VEC1,1,LBLK)
      CALL COPVCD(LUC,LUHC,VEC1,1,LBLK)
*
      DO
        N = N+1
        IF(NTEST.GE.5) THEN
          WRITE(6,*) ' Info for N = ', N
        END IF
*. (T^2)^N  times vector on LUSC1
C?     WRITE(6,*) ' Input vector to MV7 '
C?     CALL WRTVCD(VEC1,LUSC1,1,LBLK)
*.  T   * 1/(N-1)! (T^2)^(N-1)
        CALL SIG_GCC(VEC1,VEC2,LUSC1,LUSC2,TAMP)
*.  T^2 * 1/(N-1)! (T^2)^(N-1)
        CALL SIG_GCC(VEC1,VEC2,LUSC2,LUSC3,TAMP)
        IF(NTEST.GE.500.AND.IOPTYP.NE.0) THEN
          WRITE(6,*) ' 1/(N-1)! (T^2)**(N-1) |0> '
          WRITE(6,*) ' =================================='
          CALL WRTVCD(VEC1,LUSC3,1,LBLK)
        END IF

        FAC = ALPHA/DBLE(N)

        IF(IOPTYP.NE.0) THEN
* Part for unitary/hermitean operators:
          STOP 'NOT PREPARED FOR IOPTYPE.NE.0'
          CALL SCLVCD(LUSC2,LUSC3,FAC,VEC1,1,LBLK)
          CALL CONJ_T
          CALL SIG_GCC(VEC1,VEC2,LUSC1,LUSC2,TSCR)
          CALL CONJ_T
          IF(NTEST.GE.500) THEN
            WRITE(6,*) ' 1/(N-1)! T^+ (T +/- T^+)**(N-1) |0> '
            WRITE(6,*) ' =================================='
            IF (NTEST.GE.5000) THEN
              CALL WRTVCD(VEC1,LUSC2,1,LBLK)
            ELSE
              CALL WRTVSD(VEC1,LUSC2,1,LBLK)
            END IF
          END IF
c                                      in1   in2   res         
          CALL VECSMD(VEC1,VEC2,FAC,1d0,LUSC2,LUSC3,LUSC1,1,LBLK)
        ELSE
* Part for unsymmetric operators:
          CALL SCLVCD(LUSC3,LUSC1,FAC,VEC1,1,LBLK)
        END IF
        IF(NTEST.GE.500) THEN
          WRITE(6,*) ' 1/N! (T**2)**(N) |0> '
          WRITE(6,*) ' ================'
          IF (NTEST.GE.5000) THEN
            CALL WRTVCD(VEC1,LUSC1,1,LBLK)
          ELSE
            CALL WRTVSD(VEC1,LUSC1,1,LBLK)
          END IF
        END IF
*. Norm of this correction term
c       XNORM2 = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
c       XNORM = SQRT(XNORM2)
c I prefer the maximum-norm:
        XMXNRM = FDMNXD(LUSC1,2,VEC1,1,LBLK)
        IF(NTEST.GE.5) THEN
          WRITE(6,*) ' Max.-norm of correction ', XMXNRM
        END IF
*. Update output file with 1/N! T^N !ref>
        ONE = 1.0D0
        CALL VECSMD(VEC1,VEC2,ONE,ONE,LUSC1,LUHC,LUSC2,1,LBLK)
        CALL COPVCD(LUSC2,LUHC,VEC1,1,LBLK)
*. give up?
        IF (XMXNRM.GT.1d+100) THEN
          WRITE(6,*) 'Wavefunction blows up! Take a step back :-)'
          WRITE(6,*) ' Norm of last 1/N! T^N !ref>: ',XMXNRM,' for N=',N
          XNORM=SQRT(INPROD(TAMP,TAMP,N_CC_AMP))
          WRITE(6,*) ' Norm of T was: ', XNORM
          STOP 'WOOMM!'
        END IF
*. Finito ?
        IF (XMXNRM.LE.THRES_C .OR. N.GE.MX_TERM) EXIT

      END DO
*. NOTE: Result on LUHC
*
* Not converged ?
      IF (XMXNRM.GT.THRES_C) THEN
        WRITE(6,'(x,a,i5,a)')
     $        'Fatal: No convergence in EXPT_REF (max. iter.:',
     $        MX_TERM, ' )'
        STOP 'No convergence in EXPT_REF!'
      END IF
C      CALL COPVCD(LUSC3,LUHC,VEC1,1,LBLK)
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Convergence obtained in ', N, ' iterations'
        WRITE(6,*) ' Max.-norm of last correction ', XMXNRM
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' Exp (T^2) |ref> '
        WRITE(6,*) ' ==============='
        WRITE(6,*)
         IF (NTEST.GE.1000) THEN
           CALL WRTVCD(VEC1,LUHC,1,LBLK)
         ELSE
           CALL WRTVSD(VEC1,LUHC,1,LBLK)
         END IF
      END IF
*
      RETURN
      END 
*------------------------------------------------------------------------*
*--------------------------------------------------------------------*
      subroutine can2str(iway,gcan,gstr,ntss_tp,itss_tp,ibtss_tp)
*--------------------------------------------------------------------*
*
*     Set up elements of operator G in spinstring ordering using
*     operator G' in canonical, symmetry-blocked ordering
*
*      iway == 1  :   canonical -> string
*      iway == 2  :   canonical <- string
*
*--------------------------------------------------------------------*

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'lucinp.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 1000

* input
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp),
     &     ibtss_tp(ntss_tp)

      real(8), intent(inout) ::
     &     gcan(*), gstr(*)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4), isym_c(4), isymoff(nsmst)

      if (iway.ne.1.and.iway.ne.2) then
        write(6,*) 'can2str: illegal value for iway: ',iway
        stop 'can2str'
      end if

      if (ntest.ge.500) then
        write(6,*) ' iway = ',iway
        write(6,*) 'Input operator'
        if (iway.eq.1) then
          call aprblm2(gcan,ntoobs,ntoobs,nsmst,0)
        else if (iway.eq.2) then
          call wrt_cc_vec2(gstr,6,'GEN_CC')
        end if
      end if
      
      ! get symmetry offsets (for 1-particle operators)
      idx = 0
      do ism = 1, nsmst
        isymoff(ism) = idx
        idx = idx + ntoobs(ism)*ntoobs(ism)
      end do
      nlen = idx

      ! now we loop over the elements in string-ordered form

      ! loop over types
      idx = 0
      do itss = 1, ntss_tp
c        if (ibtss_tp(itss).ne.idx+1) then
c          write(6,*) 'problem with offset for op. ',itss
c          write(6,*) '  ',ibtss_tp(itss),' != ',idx+1
c        end if
        ! identify two-particle excitations:
        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index n-tupel
                      ! ireots: translates type-ordering to symmetry-ordering
                      ! ibso:   orbital-offset for symmetry
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = ireots(iocc_ca((ica-1)*nel_ca+iel))
     &                       -ibso(ism_ca) + 1
                        idx_s(ii) = 1
                        isym_c(ii) = ism_ca
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = ireots(iocc_cb((icb-1)*nel_cb+iel))
     &                       -ibso(ism_cb) + 1
                        idx_s(ii) = 2
                        isym_c(ii) = ism_cb
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = ireots(iocc_aa((iaa-1)*nel_aa+iel))
     &                       -ibso(ism_aa) + 1
                        idx_s(ii) = 1
                        isym_c(ii) = ism_aa
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = ireots(iocc_ab((iab-1)*nel_ab+iel))
     &                       -ibso(ism_ab) + 1
                        idx_s(ii) = 2
                        isym_c(ii) = ism_ab
                      end do
                      
                      ! lots of if's in the inner loop ...
                      if (na.eq.1) then
                        ! 1-particle operators
                        idxpq = (idx_s(1)-1)*nlen +
     &                       isymoff(isym_c(1)) +
     &                       (idx_c(2)-1)*ntoobs(isym_c(1)) + idx_c(1)
                        if (idx_s(1).ne.idx_s(2)) stop 'flip (1)'

                        if (iway.eq.1) then
                          gstr(idx) = gcan(idxpq)
                        else if (iway.eq.2) then
                          gcan(idxpq) = gstr(idx)
                        end if
                        
                      else if (na.eq.2) then
                        ! 2-particle operators
                        stop 'too lazy'

                      end if

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss
      if (ntest.ge.500) then
        write(6,*) ' iway = ',iway
        write(6,*) 'Output operator'
        if (iway.eq.1) then
          call wrt_cc_vec2(gstr,6,'GEN_CC')
        else if (iway.eq.2) then
          call aprblm2(gcan,ntoobs,ntoobs,nsmst,0)
        end if
      end if

      return

      end
*--------------------------------------------------------------------*
      subroutine l2g(lop,gop,ntss_tp,itss_tp,ibtss_tp,nloff,ldiml)
*--------------------------------------------------------------------*
*
*     Set up elements of two-particle operator G according to
*
*     G(pq,rs)a_pqrs = L(pq)L(rs)a_pqrs
*
*     G is blocked over operator types, each of these symmetry-blocked 
*     and in string ordering.
*
*     L is quadratic array p,q running over indices in type ordering
*     and includes also frozen or deleted orbitals, which are 
*     ignored when setting up G.
*
*     Probably not the most elegant routine on earth, but at least
*     it works....
*
*--------------------------------------------------------------------*

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 1000

* input
      real(8), intent(in) ::
     &     lop(*)
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp),
     &     ibtss_tp(ntss_tp)

      real(8), intent(out) ::
     &     gop(*)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)

      ! loop over types
      idx = 0
      do itss = 1, ntss_tp
        if (ibtss_tp(itss).ne.idx+1) then
          write(6,*) 'problem with offset for op. ',itss
          write(6,*) '  ',ibtss_tp(itss),' != ',idx+1
        end if
        ! identify two-particle excitations:
        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab
        if (na.ne.2) cycle

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

c                      print *,idx,'-> (',idx_c(1:4),')'
c                      print *,'    -> (',idx_s(1:4),')'
                      
                      idxpq = idx_s(1)*nloff +
     &                     (idx_c(3)-1)*ldiml + idx_c(1)
                      if (idx_s(1).ne.idx_s(3)) stop 'ups (1)'
                      idxrs = idx_s(2)*nloff +
     &                     (idx_c(4)-1)*ldiml + idx_c(2)
                      if (idx_s(2).ne.idx_s(4)) stop 'ups (2)'

                        gop(idx) = lop(idxpq)*lop(idxrs)
c                      print *,' gop(',idxpq,idxrs,')=',gop(idx)
c                      print *,' ',lop(idxpq),lop(idxrs)

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss
      if (ntest.ge.1000) then
        write(6,*) 'the two-particle operator:'
        call wrt_cc_vec2(gop,6,'GEN_CC')
      end if


      return

      end

*--------------------------------------------------------------------*
      subroutine pack_g(iway,idum,isymG,gop_pack,gop,
     &     ntss_tp,itss_tp,ibtss_tp,
     &     n11amp,n33amp,ioff_amp_pack,n_cc_amp)
*--------------------------------------------------------------------*
*
*     pack G from form defined by ntss_tp to usual lower triangle
*     used for 2-el. integrals (for closed shell cases)
*
*     iway: 2 pack and symmetrize
*           1 pack (no symmetrizations)
*          -1 unpack
*
*     be careful with changes:
*       1 and -1 should pack and unpack giving the same vector again
*       AND: a packed gradient should be exactly the gradient wrt. the
*       packed amplitudes (!!), else the optimization routines will go
*       gaga ....
*
*--------------------------------------------------------------------*

c      include 'implicit.inc'
c      include 'mxpdim.inc'
      include 'wrkspc.inc'
      include 'cgas.inc'
      include 'multd2h.inc'
      include 'orbinp.inc'
      include 'lucinp.inc'
      include 'csm.inc'
      include 'ctcc.inc'
      include 'glbbas.inc'
      include 'cc_exc.inc'
      
      integer, parameter ::
     &     ntest = 000
      real(8), parameter ::
c     &     f1 = 1d0,
c     &     f2 = 1.73205080756887729352d0  ! sqrt(3)
     &     f1 =  .70710678118654752440d0, ! sqrt(0.5)
     &     f2 = 1.22474487139158904909d0  ! sqrt(1.5)

* input
      real(8), intent(inout) ::
     &     gop(*)
c      input needed: itss_tp <-- work(klsobex), ntss_tp <-- nspobex_tp
      integer, intent(in) ::
     &     ntss_tp,
     &     itss_tp(ngas,4,ntss_tp),
     &     ibtss_tp(ntss_tp),
     &     ioff_amp_pack(*)

      real(8), intent(inout) ::
     &     gop_pack(*)

* local
      integer ::
     &     igrp_ca(mxpngas), igrp_cb(mxpngas),
     &     igrp_aa(mxpngas), igrp_ab(mxpngas),
     &     iocc_ca(mx_st_tsoso_blk_mx),
     &     iocc_cb(mx_st_tsoso_blk_mx),
     &     iocc_aa(mx_st_tsoso_blk_mx),
     &     iocc_ab(mx_st_tsoso_blk_mx),
     &     idx_c(4), idx_s(4)

      if (ntest.ge.10) then
        write(6,*) '================'
        write(6,*) ' this is pack_g'
        write(6,*) '================'
        print *,'iway = ', iway
        print *,'isymG = ', isymG
        print *,'ntss_tp:', ntss_tp
        print *,'ibtss_tp: ',ibtss_tp
        print *,'n11amp,n33amp,n_cc_amp: ',n11amp,n33amp,n_cc_amp

        if (ntest.ge.1000) then
          if (iway.gt.0) then
            print *,'input vector:'
            call wrt_cc_vec2(gop,6,'GEN_CC')
          else
            print *,'input packed vector (11 part):'
            call wrtmat(gop_pack,n11amp,1,n11amp,1)
            print *,'input packed vector (33 part):'
            call wrtmat(gop_pack(n11amp+1),n33amp,1,n11amp,1)
          end if
        end if
      end if

      iap_off = nsmob**3+1
      ittoff = n11amp+1

      if (iway.ne.1.and.iway.ne.2.and.iway.ne.3.and.iway.ne.-1) then
        write(6,*) 'strange iway = ', iway
        stop 'pack_G'
      end if
      if (isymG.ne.1.and.isymG.ne.-1) then
        write(6,*) 'pack_G called for non-symmetric G ',isymG
        stop 'pack_G'
      end if

      if (iway.ge.1)  gop_pack(1:n11amp+n33amp) = 0d0
      if (iway.le.-1) gop(1:n_cc_amp) = 0d0

      ! loop over types
      do itss = 1, ntss_tp
        idx = ibtss_tp(itss) - 1
c        if (ibtss_tp(itss).ne.idx+1) then
c          write(6,*) 'problem with offset for op. ',itss
c          write(6,*) '  ',ibtss_tp(itss),' != ',idx+1
c        end if
        ! identify two-particle excitations:

        nel_ca = ielsum(itss_tp(1,1,itss),ngas)
        nel_cb = ielsum(itss_tp(1,2,itss),ngas)
        nel_aa = ielsum(itss_tp(1,3,itss),ngas)
        nel_ab = ielsum(itss_tp(1,4,itss),ngas)
        nc = nel_ca + nel_cb
        na = nel_aa + nel_ab
        if (na.ne.2) stop 'accept only G2, not G1+G2 !'

        ! skip all aa or bb operators on packing
        ! (only bb case for gradient packing)
        if ((iway.eq.1.or.iway.eq.2).and.
     &       (nel_ca.eq.2.or.nel_cb.eq.2)) cycle

        ! transform occupations to groups
        call occ_to_grp(itss_tp(1,1,itss),igrp_ca,1)
        call occ_to_grp(itss_tp(1,2,itss),igrp_cb,1)
        call occ_to_grp(itss_tp(1,3,itss),igrp_aa,1)
        call occ_to_grp(itss_tp(1,4,itss),igrp_ab,1)

        if (mscomb_cc.ne.0) then
          call diag_exc_cc(itss_tp(1,1,itss),itss_tp(1,2,itss),
     &                     itss_tp(1,3,itss),itss_tp(1,4,itss),
     &                     ngas,idiag)
        else
          idiag = 0
        end if
        
        ! loop over symmetry blocks
        ism = 1 ! totally symmetric operators, n'est-ce pas?
        do ism_c = 1, nsmst
          ism_a = multd2h(ism,ism_c)
          do ism_ca = 1, nsmst
            ism_cb = multd2h(ism_c,ism_ca)
            do ism_aa = 1, nsmst
              ism_ab = multd2h(ism_a,ism_aa)
              ! get alpha and beta symmetry index
              ism_alp = (ism_aa-1)*nsmst+ism_ca  ! = (sym Ca,sym Aa)
              ism_bet = (ism_ab-1)*nsmst+ism_cb  ! = (sym Cb,sym Ab)
              
              ! restrict to (sym Ca,sym Aa) >= (sym Cb,sym Ab)
              if (idiag.eq.1.and.ism_bet.gt.ism_alp) cycle
              if (idiag.eq.0.or.ism_alp.gt.ism_bet) then
                irestr = 0
              else
                irestr = 1
              end if
              
              ! get the strings
              call getstr2_totsm_spgp(igrp_ca,ngas,ism_ca,nel_ca,
     &             lca,iocc_ca,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_cb,ngas,ism_cb,nel_cb,
     &             lcb,iocc_cb,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_aa,ngas,ism_aa,nel_aa,
     &             laa,iocc_aa,norb,0,idum,idum)
              call getstr2_totsm_spgp(igrp_ab,ngas,ism_ab,nel_ab,
     &             lab,iocc_ab,norb,0,idum,idum)

              ! length of strings in this symmetry block
              if (lca*lcb*laa*lab.eq.0) cycle

              do iab = 1, lab
                if (irestr.eq.1) then
                  iaa_min = iab
                else
                  iaa_min = 1
                end if
                do iaa = iaa_min, laa
                  do icb = 1, lcb
                    if (irestr.eq.1.and.iaa.eq.iab) then
                      ica_min = icb
                    else
                      ica_min = 1
                    end if
                    do ica = ica_min, lca
                      idx = idx + 1
                      ! translate into canonical index quadrupel
                      ii = 0
                      do iel = 1, nel_ca
                        ii = ii + 1
                        idx_c(ii) = iocc_ca((ica-1)*nel_ca+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_cb
                        ii = ii + 1
                        idx_c(ii) = iocc_cb((icb-1)*nel_cb+iel)
                        idx_s(ii) = 2
                      end do
                      do iel = 1, nel_aa
                        ii = ii + 1
                        idx_c(ii) = iocc_aa((iaa-1)*nel_aa+iel)
                        idx_s(ii) = 1
                      end do
                      do iel = 1, nel_ab
                        ii = ii + 1
                        idx_c(ii) = iocc_ab((iab-1)*nel_ab+iel)
                        idx_s(ii) = 2
                      end do

                      idxp = idx_c(1) 
                      idxr = idx_c(2)
                      idxq = idx_c(3) 
                      idxs = idx_c(4)

                      idxpr = (min(idxp,idxr)-1)*ntoob+max(idxp,idxr)
                      idxqs = (min(idxq,idxs)-1)*ntoob+max(idxq,idxs)

                      if (iway.eq.1) then
                        ! packing
                        ! take only triangle
                        if (idxp.gt.idxr .or. idxq.gt.idxs) cycle
 
                        iadr = i2addr2(ireots(idxp),ireots(idxr),
     &                              ireots(idxq),ireots(idxs),
     &                              ioff_amp_pack,1,1,isymG)
                        if (iadr.lt.0) cycle
                        if (iadr.eq.0) stop 'iadr error'

                        gop_pack(iadr) = gop(idx)

                        if (ntest.ge.1000)
     &                    print '(a,2i4,a,4i4,i5,2(x,e12.6))',
     &                       '1S',itss,idx-ibtss_tp(itss)+1,'->',
     &                       idxp,idxr,idxq,idxs,iadr,
     &                       gop(idx)
     &                       

                        if (idxp.eq.idxr .or. idxq.eq.idxs) cycle

                        iadr = i2addr2(ireots(idxp),ireots(idxr),
     &                              ireots(idxq),ireots(idxs),
     &                              ioff_amp_pack(iap_off),-1,-1,isymG)
                        if (iadr.lt.0) cycle
                        if (iadr.eq.0) stop 'iadr error'

                        gop_pack(ittoff+iadr) = gop(idx)

                        if (ntest.ge.1000)
     &                    print '(a,2i4,a,4i4,i5,2(x,e12.6))',
     &                       '1T',itss,idx-ibtss_tp(itss)+1,'->',
     &                       idxp,idxr,idxq,idxs,iadr,
     &                       gop(idx)

*----------------------------------------------------------------------*
*     2: packing and (anti-)symmetrizing
*----------------------------------------------------------------------*
                      else if (iway.eq.2) then
*----------------------------------------------------------------------*
*     2A: contributions to G(+):
*----------------------------------------------------------------------*

                        fac = f1 
                        if (idxpr.ge.idxqs) fac = dble(isymG)*fac

c     &                       sqrt(dble(isymfac(idxp,idxr,idxq,idxs)))
c                        if (idxp.eq.idxr) fac = fac/2d0
c                        if (idxq.eq.idxs) fac = fac/2d0
c                        if (idxp.eq.idxr.and.idxq.eq.idxs) fac = fac/2d0
C                        if (idxp.eq.idxr.or.idxq.eq.idxs) fac = fac/2d0

       write(6,*) ' Jeppe commented this out to get code running '
C                       if ( idxp.eq.idxr.xor.idxq.eq.idxs)
C    &                       fac = fac*sqrt(2d0)

                        if (idxp.le.idxr.and.idxq.le.idxs) then
                          iadr = i2addr2(ireots(idxp),ireots(idxr),
     &                         ireots(idxq),ireots(idxs),
     &                         ioff_amp_pack,1,1,isymG)
                          if (iadr.lt.0) cycle
                          if (iadr.eq.0) stop 'iadr error'

c                          sfac = 1d0
c                          if (idxp.eq.idxr.or.idxq.eq.idxs) sfac = 2d0
                          gop_pack(iadr) = gop_pack(iadr)
     &                         + fac*gop(idx)

                          if (ntest.ge.1000)
     &                      print '(a,2i4,i5,a,4i4,i5,3(x,e12.6))',
     &                       '2S1',itss,idx-ibtss_tp(itss)+1,idx,'->',
     &                       idxp,idxr,idxq,idxs,iadr,
     &                       gop(idx),gop_pack(iadr),fac

                        else if (idxp.lt.idxr.and.idxq.ne.idxs) then
                          iadr = i2addr2(ireots(idxp),ireots(idxr),
     &                         ireots(idxs),ireots(idxq),
     &                         ioff_amp_pack,1,1,isymG)
                          if (iadr.lt.0) cycle
                          if (iadr.eq.0) stop 'iadr error'

                          gop_pack(iadr) = gop_pack(iadr)
     &                         + fac*gop(idx)
                          if (ntest.ge.1000)
     &                      print '(a,2i4,i5,a,4i4,i5,3(x,e12.6))',
     &                       '2S2',itss,idx-ibtss_tp(itss)+1,idx,'->',
     &                       idxp,idxr,idxs,idxq,iadr,
     &                       gop(idx),gop_pack(iadr),fac

                        end if

*----------------------------------------------------------------------*
*     2B: contributions to G(-):
*----------------------------------------------------------------------*
                        fac = f2
                        if (idxpr.ge.idxqs) fac = dble(isymG)*fac

                        if (idxp.lt.idxr.and.idxq.lt.idxs) then
                          iadr = i2addr2(ireots(idxp),ireots(idxr),
     &                         ireots(idxq),ireots(idxs),
     &                         ioff_amp_pack(iap_off),-1,-1,isymG)
                          if (iadr.lt.0) cycle
                          if (iadr.eq.0) stop 'iadr error'
                          
                          gop_pack(ittoff+iadr) = gop_pack(ittoff+iadr)
     &                         + fac*gop(idx)

                          if (ntest.ge.1000)
     &                      print '(a,2i4,i5,a,4i4,i5,3(x,e12.6))',
     &                       '2T1',itss,idx-ibtss_tp(itss)+1,idx,'->',
     &                       idxp,idxr,idxs,idxq,iadr,
     &                       gop(idx),gop_pack(ittoff+iadr),fac

                        else if (idxp.lt.idxr.and.idxq.ne.idxs) then
                          iadr = i2addr2(ireots(idxp),ireots(idxr),
     &                         ireots(idxs),ireots(idxq),
     &                         ioff_amp_pack(iap_off),-1,-1,isymG)

                          if (iadr.lt.0) cycle
                          if (iadr.eq.0) stop 'iadr error'
                          
                          gop_pack(ittoff+iadr) = gop_pack(ittoff+iadr)
     &                         - fac*gop(idx)
                          if (ntest.ge.1000)
     &                      print '(a,2i4,i5,a,4i4,i5,3(x,e12.6))',
     &                       '2T2',itss,idx-ibtss_tp(itss)+1,idx,'->',
     &                       idxp,idxr,idxs,idxq,iadr,
     &                       gop(idx),gop_pack(ittoff+iadr),-fac

                        end if

*----------------------------------------------------------------------*
*     -1: unpacking
*----------------------------------------------------------------------*
                      else
*----------------------------------------------------------------------*
*     -1A: unpack contrib.s from G(-) to either G(aa),G(bb) or G(ab)
*----------------------------------------------------------------------*
c                        fac =
c     &                      1d0/sqrt(dble(isymfac(idxp,idxr,idxq,idxs)))
                        if (idx_s(1).eq.idx_s(2)) then
                          fac = 0.5d0/f2
                        else
                          fac = 0.5d0/f2
                          if (idxp.ne.idxr.and.idxq.ne.idxs)
     &                         fac = fac/2d0
                        end if

                        if ( idxpr.gt.idxqs)
     &                       fac = dble(isymG)*fac

c                        if (idx_s(1).ne.idx_s(2))
c     &                       fac = 0.5d0*fac

                        sfac = 1d0
                        if (idxp.gt.idxr) sfac = sfac*(-1d0)
                        if (idxq.gt.idxs) sfac = sfac*(-1d0)
                        iadr = i2addr2(ireots(idxp),ireots(idxr),
     &                              ireots(idxq),ireots(idxs),
     &                              ioff_amp_pack(iap_off),-1,-1,isymG)
                      
                        if (iadr.ge.0) then
                          if (iadr.eq.0) stop 'iadr error'
                          if (iadr.gt.n33amp) then
                            print *,'1: ',idxp,idxq,idxr,idxs
                            print *,'2: ',ioff_amp_pack(1:3)
                            stop 'error error'
                          end if

                          gop(idx) = sfac*fac*gop_pack(ittoff+iadr)

                          if (ntest.ge.1000)
     &                      print '(a,2i4,i5,a,4i4,i5,3(x,e12.6))',
     &                       '3:-',itss,idx-ibtss_tp(itss)+1,idx,'<-',
     &                       idxp,idxr,idxq,idxs,iadr,
     &                       gop(idx),gop_pack(ittoff+iadr),sfac*fac
     &                       
                        end if

                        if (idx_s(1).eq.idx_s(2)) cycle
*----------------------------------------------------------------------*
*     -1B: unpack contrib.s from G(+) to G(ab)
*----------------------------------------------------------------------*

c                        fac = 1.0d0
                        fac = 0.5d0/f1
 
       write(6,*) ' Jeppe commented this out to get code running '
C                       if ( idxp.eq.idxr.xor.idxq.eq.idxs)
C    &                       fac = fac/sqrt(2d0)

                        if ( idxpr.gt.idxqs)
     &                       fac = dble(isymG)*fac
                        if (idxp.ne.idxr.and.idxq.ne.idxs)
     &                       fac = fac/2d0

c                        if (idxp.eq.idxr) fac = fac*2d0
c                        if (idxq.eq.idxs) fac = fac*2d0

                        iadr = i2addr2(
     &                              ireots(idxp),ireots(idxr),
     &                              ireots(idxs),ireots(idxq),
     &                              ioff_amp_pack,1,1,isymG)

                        if (iadr.lt.0) cycle
                        if (iadr.eq.0) stop 'iadr error'
                        gop(idx) = gop(idx)+fac*gop_pack(iadr)

                        if (ntest.ge.1000)
     &                    print '(a,2i4,i5,a,4i4,i5,3(x,e12.6))',
     &                       '3:+',itss,idx-ibtss_tp(itss)+1,idx,'<-',
     &                       idxp,idxr,idxs,idxq,iadr,
     &                       gop(idx),gop_pack(iadr),fac
     &                       
                      end if

                    end do ! ica
                  end do ! icb
                end do ! iaa
              end do ! iab

            end do ! ism_aa
          end do ! ism_ca
        end do ! ism_c 

      end do ! itss

      if (ntest.eq.1000) then
        if (iway.gt.0) then
          print *,'packed vector (11)'
          call wrtmat(gop_pack,n11amp,1,n11amp,1)
          print *,'packed vector (33)'
          call wrtmat(gop_pack(n11amp+1),n33amp,1,n33amp,1)
        else
          print *,'unpacked vector'
          call wrt_cc_vec2(gop,6,'GEN_CC')
        end if
      end if

      return

      end

      integer function isymfac(ip,ir,iq,is)
* return the number of non-identical permutations of the index-quadruple
* unter (anti-)hermitian and particle symmetry
*                 
*  identity         (ip,ir,iq,is)
*  hermitian conj.  (iq,is,ip,ir)
*  particle perm.   (ir,ip,is,iq)
*  h. c. + p. p.    (is,iq,ir,ip)
*          
      implicit none

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ip,ir,iq,is

      integer ::
     &     ifac

      ifac = 1
      ! symmetric under herm. conj.?
      if (.not.(ip.eq.iq.and.ir.eq.is)) ifac = ifac*2
      ! symmetric under particle perm.?
      if (.not.(ip.eq.ir.and.iq.eq.is)) ifac = ifac*2
      ! symmetric under combination of both?
      if (ifac.eq.4 .and.
     &    (ip.eq.is.and.iq.eq.ir)) ifac = ifac/2

      if (ntest.ge.100)
     &     write(6,'(x,a,4i10,a,i2)')
     &     'isymfac: ',ip,ir,iq,is,' --> ',ifac

      isymfac = ifac
      return

      end

      subroutine set_frobs(nfrob,nfrobs)

      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'lucinp.inc'
      include 'csm.inc'
      include 'csmprd.inc'
      include 'cgas.inc'

      dimension nfrobs(nsmob)

      isym = 1

      nfrobs(1:nsmob) = 0
      do igas = 1, ngas-1
        ! excitations allowed in this GAS space?
        nrem = igsocc(igas,2)-igsocc(igas,1)
        if (nrem.gt.0) exit
        nfrobs(1:nsmob) = nfrobs(1:nsmob)+ngssh(1:nsmob,igas)
        print *,'1> ',igspc,nrem
        print *,'   ',nfrobs(1:nsmob)
      end do

      nfrob = sum(nfrobs,nsmob)

      print *,'final suggestion:'
      print *,'  ',nfrobs(1:nsmob)
      print *,' >',nfrob

      return
      end

      subroutine num_ssaa2op(nndiag,ndiag)

*     find the number of symmetry and spin (i.e. singlet) adapted 
*     antisymmetric two-body operators
      include 'implicit.inc'
      include 'mxpdim.inc'
      include 'lucinp.inc'
      include 'csm.inc'
      include 'csmprd.inc'
      include 'cgas.inc'

      logical lpdiag,lhdiag,lhpdiag
      dimension iact(ngas)

      isym = 1

      do igspc = 1, ngas
        ! 2-body excitations allowed in this GAS space?
        nrem = igsocc(igspc,2)-igsocc(igspc,1)
        nadd = 0
        if (igspc.gt.1)
     &       nadd = igsocc(igspc,2)-igsocc(igspc-1,1)
        if (nrem.ne.0.or.nadd.ne.0) then
          iact(igspc) = 1
        else
          iact(igspc) = 0
        end if
        print *,'1> ',igspc,nrem,nadd
      end do

      print *,'-> ',iact(1:ngas)

      isum = 0
      isumd = 0
      do ip1spc = 1, ngas
        if (iact(ip1spc).eq.0) cycle
        do ip2spc = 1, ip1spc
          if (iact(ip2spc).eq.0) cycle
          lpdiag = ip1spc.eq.ip2spc
          ipidx = (ip1spc-1)*ngas + ip2spc
          do ih1spc = 1, ngas
            if (iact(ih1spc).eq.0) cycle
            do ih2spc = 1, ih1spc
              if (iact(ih2spc).eq.0) cycle
              lhdiag = ih1spc.eq.ih2spc
              ihidx = (ih1spc-1)*ngas + ih2spc
              lhpdiag = ihidx.eq.ipidx
              if (ihidx.gt.ipidx) cycle

              ii12 = 0
              ii34 = 0
              ii1234 = 0
              if (lpdiag) ii12 = 1
              if (lhdiag) ii34 = 1
              if (lhpdiag) ii1234 = 1

              print *,'>> ',lpdiag,lhdiag,lhpdiag

              print *,' p1: ',ngssh(1:nirrep,ip1spc)
              print *,' p2: ',ngssh(1:nirrep,ip2spc)
              print *,' h1: ',ngssh(1:nirrep,ih1spc)
              print *,' h2: ',ngssh(1:nirrep,ih2spc)

              inum = ndxfsm(nsmob,nsmsx,mxpobs,
     &             ngssh(1,ip1spc),ngssh(1,ip2spc),
     &             ngssh(1,ih1spc),ngssh(1,ih2spc),
     &             isym,adsxa,sxdxsx,ii12,ii34,ii1234,0)

              idiag = 0
              if (lhpdiag) then
                do ii = 1, nirrep
                  do jj = 1, nirrep
                    idiag = idiag + ngssh(ii,ip1spc)*ngssh(jj,ih1spc)
                  end do
                end do
              end if

              print '(a,4i3,2i8)','> ',ip1spc,ip2spc,ih1spc,ih2spc,inum,
     &             idiag

              isum = isum + inum
              isumd = isumd + idiag

            end do
          end do
        end do
      end do

      ndiag = isumd
      nndiag = isum-isumd

      return
      
      end
*----------------------------------------------------------------------*
      subroutine gtbce_h0(imode,igtb_closed,isymmet_G,
     &                    iccvec,nSdim,
     &                    ccvec1,ccvec2,ccvec3,
     &                    civec1,civec2,c2vec,
     &                    n_cc_amp,mxb_ci,
     &                    n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                    n11amp,n33amp,iamp_packed,
     &                    luh0,ludia,
     &                    luamp,luec,luhc,
     &                    lusc1,lusc2)
*----------------------------------------------------------------------*
*
*     imode == 0: <ref|exp(-G)tau(mu)tau(nu)exp(G)|ref>
*     imode == 1: 2<ref|exp(-G)tau(mu)H tau(nu)exp(G)|ref>
*                -2<ref|exp(-G)H tau(mu)tau(nu)exp(G)|ref>
*     imode == 2: dto. and save diagonal on ludia
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     isymmet_G, igtb_closed, n11amp, n33amp,
     &     n_cc_amp, nsdim, mxb_ci,
     &     luamp, luec, luhc, luh0, ludia,
     &     lusc1, lusc2, iamp_packed(*), iccvec(n_cc_amp),
     &     n_cc_typ(*), i_cc_typ(*), ioff_cc_typ(*)

      real(8), intent(inout) ::
     &     ccvec1(n_cc_amp), ccvec2(n_cc_amp), ccvec3(n_cc_amp),
     &     civec1(mxb_ci), civec2(mxb_ci), c2vec(*)

      integer ::
     &     iamp, iadj, lblk, isigden, idx, ii, imode, icnt,
     &     namp, nsave, iway, idum
      real(8) ::
     &     fac, xmin, xsh,
     &     wall0, wall, cpu0, cpu
      
      real(8), external ::
     &     inprod
 
      call atim(cpu0,wall0)

      if (ntest.gt.0) then
        write(6,*) '====================='
        write(6,*) ' here comes gtbce_h0'
        write(6,*) '====================='
        write(6,*) ' isymmet_G, igtb_closed : ',isymmet_G, igtb_closed
        write(6,*) ' nSdim, n_cc_amp: ',nSdim,n_cc_amp
        write(6,*) ' luh0, luamp, luec, luhc: ',luh0, luamp, luec, luhc
      end if

      call rewino(luh0)

      icnt = 0
      lblk = -1
      fac = dble(isymmet_G)
      namp = n_cc_amp
      if (igtb_closed.eq.1) then
        namp = n11amp+n33amp
      end if
      do iamp = 1, namp

        if (ntest.ge.10) write(6,*) ' iamp = ',iamp,'/',namp

        if (igtb_closed.eq.1.and.isymmet_G.ne.0) then
          if (iccvec(iamp).lt.0) cycle
        end if
        icnt = icnt+1

        if (isymmet_G.ne.0.and.igtb_closed.eq.1) then
          ccvec1(1:namp) = 0d0
          ccvec1(iamp) = 1d0
          iway = -1
          call pack_g(iway,idum,isymmet_G,ccvec1,ccvec2,
     &                n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                n11amp,n33amp,iamp_packed,n_cc_amp)

        else if (isymmet_G.eq.0) then
          ccvec2(1:namp) = 0d0
          ccvec2(iamp) = 1d0
        else
          ccvec2(1:namp) = 0d0
          iadj = abs(iccvec(iamp))
          ccvec2(iamp) = sqrt(2d0)
          ccvec2(iadj) = fac*sqrt(2d0)
        end if

*----------------------------------------------------------------------*
*     calculate tau_(iamp)exp(G)|ref>
*----------------------------------------------------------------------*
        isigden=1
        call sigden_cc(civec1,civec2,luec,lusc1,ccvec2,isigden)


        if (imode.ge.1) then
*----------------------------------------------------------------------*
*     calculate H tau_(iamp)exp(G)|ref>
*----------------------------------------------------------------------*
          call mv7(civec1,civec2,lusc1,lusc2)

*----------------------------------------------------------------------*
*     1: <ref|exp(G)tau_(iamp) tau(nu) H exp(G)|ref>
*----------------------------------------------------------------------*
          ccvec1(1:n_cc_amp) = 0d0
          isigden = 2
          call sigden_cc(civec1,civec2,luhc,lusc1,ccvec1,isigden)

*----------------------------------------------------------------------*
*     2: <ref|exp(G)tau_(iamp) H tau(nu) exp(G)|ref>
*----------------------------------------------------------------------*
          ccvec2(1:n_cc_amp) = 0d0
          isigden = 2
          call sigden_cc(civec1,civec2,luec,lusc2,ccvec2,isigden)
          call vecsum(ccvec1,ccvec1,ccvec2,-2d0,2d0,n_cc_amp)
        else
*----------------------------------------------------------------------*
*     2: <ref|exp(G)tau_(iamp) tau(nu) exp(G)|ref>
*----------------------------------------------------------------------*
          ccvec1(1:n_cc_amp) = 0d0
          isigden = 2
          call sigden_cc(civec1,civec2,luec,lusc1,ccvec1,isigden)
        end if

        if (isymmet_G.ne.0.and.igtb_closed.eq.0) then
          ! collect diagonal
          iadj = abs(iccvec(iamp))
          ccvec3(iamp) = ccvec1(iamp)+fac*ccvec1(iadj)
          ccvec3(iadj) = ccvec3(iamp) ! we want them positive
          ! compress result vector
          idx = 0
          do ii = 1, n_cc_amp
            if (iccvec(ii).le.0) cycle
            idx = idx + 1
            iadj = abs(iccvec(ii))
            ccvec2(idx) = ccvec1(ii)+fac*ccvec1(iadj)
          end do
          nsave = nSdim
          if (imode.eq.0) nsave = icnt
          call vec_to_disc(ccvec2,nsave,0,lblk,luh0)
        else if (igtb_closed.eq.1) then
          ! pack again
          iway = 2
          call pack_g(iway,idum,isymmet_G,ccvec2,ccvec1,
     &                n_cc_typ,i_cc_typ,ioff_cc_typ,
     &                n11amp,n33amp,iamp_packed,n_cc_amp)          
          ccvec3(iamp) = ccvec2(iamp)
          nsave = nSdim
          if (imode.eq.0) nsave = icnt
          call vec_to_disc(ccvec2,nsave,0,lblk,luh0)
        else
          ccvec3(iamp) = ccvec1(iamp)
          nsave = nSdim
          if (imode.eq.0) nsave = icnt
          call vec_to_disc(ccvec1,nsave,0,lblk,luh0)
        end if

      end do ! iamp

      if (imode.eq.2) then
      ! look at diagonal
        xmin = 1000d0
        do ii = 1, namp
          xmin = min(ccvec3(ii),xmin)
        end do
        write(6,*) 'diagonal: lowest element = ',xmin
        xsh = max(0d0,0.01d0-xmin)
        write(6,*) 'shift diagonal by ',xsh
        do ii = 1, namp
          ccvec3(ii) = ccvec3(ii)+xsh
        end do
        if (isymmet_G.ne.0) then
          do ii = 1, namp
            if (iccvec(ii).eq.-ii) ccvec3(ii)=1d12
          end do
        end if
        call vec_to_disc(ccvec3,namp,1,-1,ludia)
      end if
        
      call atim(cpu,wall)
      call prtim(6,'time in gtbce_h0',cpu-cpu0,wall-wall0)

      return

      end
**********************************************************************
      subroutine ana_gucc(vec,n11amp,n33amp,iamp_packed,
     &                    ireost,nsmob,ntoob)

      implicit none

      integer, parameter ::
     &     ntest = 100, nlist = 20
      
      integer, intent(in) ::
     &     n11amp, n33amp, iamp_packed(*), ireost(*), nsmob, ntoob
      real(8), intent(in) ::
     &     vec(*)

      real(8) ::
     &     xlist(nlist), x11n, x33n

      integer ::
     &     ii, ilist(nlist), ijkllist(4,nlist)

      real(8), external ::
     &     inprod

      x11n = sqrt(inprod(vec,vec,n11amp))

      call list_asl(2,vec,n11amp,xlist,ilist,nlist)

      call ijkl2iadr(ijkllist,ilist,nlist,
     &               ntoob,ireost,iamp_packed,1,1,-1)

      write(6,*) 'singlet-singlet coupled part: '
      write(6,'(x,a,i10,a,g20.8)')' amplitudes: ',n11amp,'  norm: ',x11n
      write(6,*) 'largest amplitudes:'
      do ii = 1, nlist
        write(6,'(x,i8,x,4i5,g20.8)')
     &       ilist(ii),ijkllist(1:4,ii),xlist(ii)
      end do

      x33n = sqrt(inprod(vec(n11amp+1),vec(n11amp+1),n33amp))
      
      call list_asl(2,vec(n11amp+1),n33amp,xlist,ilist,nlist)

      call ijkl2iadr(ijkllist,ilist,nlist,
     &              ntoob,ireost,iamp_packed(nsmob**3+1),-1,-1,-1)

      write(6,*) 'triplet-triplet coupled part: '
      write(6,'(x,a,i10,a,g20.8)')' amplitudes: ',n33amp,'  norm: ',x33n
      write(6,*) 'largest amplitudes:'
      do ii = 1, nlist
        write(6,'(x,i8,x,4i5,g20.8)')
     &       ilist(ii),ijkllist(1:4,ii),xlist(ii)
      end do

      return

      end
**********************************************************************
      FUNCTION NCSF_FOR_CISPACE(ISPC,ISYM)
*
* Find number of CSF's, CONF's (and SD's) for given CISPACE
* and symmetry
*
* The CI space is defined by the integer ISPC
*
* The spin-multiplicity, 2*Ms and combination flags
* are obtained from MULTS, MS2 and PSSIGN in CSTATE.
*
* The symmetry is defined by ISYM
*
*
* A bit of modifications from CSFDIM_GAS for Andreas
*
* Jeppe Olsen, Aug 2004
*
*
* ( Spin signaled by PSSIGN in CIINFO)
*
c      INCLUDE 'implicit.inc'
c      INCLUDE 'mxpdim.inc'
      INCLUDE 'wrkspc.inc'
      INCLUDE 'orbinp.inc'
      INCLUDE 'cstate.inc'
      INCLUDE 'glbbas.inc'
      INCLUDE 'cgas.inc'
      INCLUDE 'spinfo.inc'
      INCLUDE 'cprnt.inc'
      INCLUDE 'gasstr.inc'
* Scratch for one occupation class
      INTEGER IOCCLS(MXPNGAS)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'NCSF_F')
*
      NTEST = 1
      NTEST = MAX(IPRCIX,NTEST)
      IF(NTEST.GE.10) WRITE(6,*) '  PSSIGN : ', PSSIGN
      IF(NTEST.GE.10) WRITE(6,*) ' MULTS, MS2 = ', MULTS,MS2
*. Obtain the occupation classes for this CISPACE
*. Number of occupation classes
*. Number
      IATP = 1
      IBTP = 2
      NEL = NELFTP(IATP)+NELFTP(IBTP)
      CALL OCCLSE(1,NOCCLS,IOCCLS,NEL,ISPC,0,0,NOBPT)
*. And the occupation classes
      CALL MEMMAN(KLOCCLS,NOCCLS*NGAS,'ADDL  ',1,'OCCLS ')
      CALL OCCLSE(2,NOCCLS,WORK(KLOCCLS),NEL,ISPC,0,0,NOBPT)
*. Number of occupation classes for T-operators
*
*.. Define parameters in SPINFO
*
*. Allowed number of open orbitals
      MINOP = ABS(MS2)
      CALL MAX_OPEN_ORB(MAXOP,WORK(KLOCCLS),NGAS,NOCCLS,NOBPT)
      IF( NTEST .GE. 2 )
     &WRITE(6,*) ' MINOP MAXOP ',MINOP,MAXOP
C
C.. Number of prototype sd's and csf's per configuration prototype
C
      ITP = 0
      DO IOPEN = 0, MAXOP
        ITP = IOPEN + 1
*. Unpaired electrons :
        IAEL = (IOPEN + MS2 ) / 2
        IBEL = (IOPEN - MS2 ) / 2
        IF(IAEL+IBEL .EQ. IOPEN .AND. IAEL-IBEL .EQ. MS2 .AND.
     &            IAEL .GE. 0 .AND. IBEL .GE. 0) THEN
          NPDTCNF(ITP) = IBION(IOPEN,IAEL)
          IF(PSSIGN.EQ. 0.0D0 .OR. IOPEN .EQ. 0 ) THEN
            NPCMCNF(ITP) = NPDTCNF(ITP)
          ELSE
            NPCMCNF(ITP) = NPDTCNF(ITP)/2
          END IF
          IF(IOPEN .GE. MULTS-1) THEN
            NPCSCNF(ITP) = IWEYLF(IOPEN,MULTS)
          ELSE
            NPCSCNF(ITP) = 0
          END IF
        ELSE
          NPDTCNF(ITP) = 0
          NPCMCNF(ITP) = 0
          NPCSCNF(ITP) = 0
        END IF
      END DO
*
      IF(NTEST.GE.1) THEN
      IF(PSSIGN .EQ. 0 ) THEN
        WRITE(6,*) '  (Combinations = Determinants ) '
      ELSE
        WRITE(6,*) '  (Spin combinations in use ) '
      END IF
      WRITE(6,'(/A)') ' Information about prototype configurations '
      WRITE(6,'( A)') ' ========================================== '
      WRITE(6,'(/A)')
     &'  Open orbitals   Combinations    CSFs '
      DO IOPEN = MINOP,MAXOP,2
        WRITE(6,'(5X,I3,10X,I6,7X,I6)')
     &  IOPEN,NPCMCNF(IOPEN+1),NPCSCNF(IOPEN+1)
      END DO
*
      END IF
C
C.. Number of Configurations per occupation type
C
      DO JOCCLS = 1, NOCCLS
        IF(JOCCLS.EQ.1) THEN
          INITIALIZE_CONF_COUNTERS = 1
        ELSE
          INITIALIZE_CONF_COUNTERS = 0
        END IF
*
        IDOREO = 0
        CALL ICOPVE2(WORK(KLOCCLS),(JOCCLS-1)*NGAS+1,NGAS,IOCCLS)
        IB_ORB = NINOB + 1
        CALL GEN_CONF_FOR_OCCLS(IOCCLS,
     &     IDUM,INITIALIZE_CONF_COUNTERS,
     &     NGAS,ISYM,MINOP,MAXOP,NSMST,1,NOCOB,
     &     NOBPT,NCONF_PER_OPEN(1,ISYM),NCONF_OCCLS,
     &     IB_CONF_REO,IB_CNOCC_OPEN,
     &     IDUM,IDOREO,IDUMMY,IDUMMY,NCONF_ALL_SYM,IB_ORB)
*
      END DO
*. Number of CSF's in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN(1,ISYM),NPCSCNF,
     &                   NCSF)
*. Number of SD's in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN(1,ISYM),NPDTCNF,
     &                    NSD)
*. Number of combinations in expansion
      CALL NCNF_TO_NCOMP(MAXOP,NCONF_PER_OPEN(1,ISYM),NPCMCNF,
     &                    NCM)
*
      NCSF_PER_SYM(ISYM) = NCSF
      NSD_PER_SYM(ISYM) = NSD
      NCM_PER_SYM(ISYM) = NCM
      NCONF_PER_SYM(ISYM) = IELSUM(NCONF_PER_OPEN(1,ISYM),MAXOP+1)
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Number of CSFs  ', NCSF
        WRITE(6,*) ' Number of SDs   ', NSD
        WRITE(6,*) ' Number of Confs ', NCONF_PER_SYM(ISYM)
        WRITE(6,*) ' Number of CMs   ', NCM
      END IF
*
      NCSF_FOR_CISPACE = NCSF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'NCSF_F')
*
      RETURN
      END
c $Id$
