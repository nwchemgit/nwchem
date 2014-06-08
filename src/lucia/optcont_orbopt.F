*----------------------------------------------------------------------*
      subroutine optcont_orbopt
     &                  (imacit,imicit,imicit_tot,iprint,
     &                   itask,iconv,iorbopt,ikapmod,itransf,
     &                   luamp,lulamp,lukappa,lutrmat,
     &                   lutrvec,lutrv_l,lutrv_o,
     &                   energy,
     &                   vec1,vec2,nwfpar,norbr,icmpr,nspin,
     &                   lugrvf,lulgrvf,luogrd,lusig,lusig_l,lusig_o,
     &                   ludia,luodia,nrdvec,lurdvec)
*----------------------------------------------------------------------*
*
* General optimization control routine for non-linear optimization.
* "Special edition" for orbital-optimized and Brueckner CC
*
*   iorbopt = 001  --  (   t,kappa)
*             002      (L, t,kappa)
*             003       L,(t,kappa)
*
*   ikapmod = 0   do not keep track of total kappa, i.e. only
*                 the U matrix is updated as U = exp(dkappa)U0
*                 does not allow for extrapolation techniques like DIIS
*
*   ikapmod = 1   the U matrix is updated as U = exp(dkappa)U0 and
*                 a new kappa is evaluated as ln(U)
*                 allows for DIIS, but ln(U) may not converge for large
*                 orbital rotations
*
*   ikapmod = 2   it is assumed that luogrd contains the real gradient at
*                 kappa0 (i.e. without prior change of the expansion point)
*                 then the new kappa is kappa = kappa0+dkappa and the new
*                 U matrix can be calculated as exp(kappa)
*
* The optimization strategy is set on common block /opti/, see there 
* for further comments.
*
* It holds track of the macro/micro-iterations and expects to be
* called with imacit/imicit set to zero at the first time.
*
* optcont() requests information via the itask parameter (bit flags):
*     
*     1: calculate energy
*     2: calculate gradient/vectorfunction
*     4: calculate Hessian/Jacobian(H/A) times trialvector product
*     8: exit
*
* Files: 
* a) passed to slave routines
*   luamp   -- current set of wave-function parameters  luamp 
*   lutrvec -- additional trial-vector for H/A matrix-vector products
*
* b) passed from slave routines
*   lugrvf  -- gradient or vectorfunction 
*   lusig   -- H/A matrix-vector product
*
* c) own book keeping files
*   lugrvfold -- previous gradient/vectorfunction 
*   lusigold  -- previous H/A matrix-vector product
*   lu_newstp -- current new step (scratch)
*   lust_sbsp -- subspace of previous steps
*   lugv_sbsp -- subspace of previous gradient/vector function diff.s
*
*----------------------------------------------------------------------*
* common blocks
      include "wrkspc.inc"  ! includes also implicit statement
      include "opti.inc"

* constants
      integer, parameter ::
     &     ntest = 00, mxptsk = 10

* interface:
      integer, intent(out) ::
     &     itask, itransf
      integer, intent(inout) ::
     &     imacit, imicit
      integer, intent(in) ::
     &     nwfpar, nspin, norbr(nspin), nrdvec,
     &     luamp, lutrvec, lutrv_l, lutrv_o, lugrvf,
     &     lusig, lusig_l, lusig_o, ludia, luodia, lurdvec,
     &     icmpr(*)
      real(8), intent(in) ::
     &     energy
* two scratch vectors
      real(8), intent(inout) ::
     &     vec1(nwfpar), vec2(nwfpar)

* some private static variables:
      integer, save ::
     &     lugrvfold, lusigold, lu_newstp, lu_corstp,
     &     lusig_c, lutrv_c, lugrvf_c, ludia_c,
     &     lust_sbsp, lupst_sbsp, lutpst_sbsp, lugv_sbsp, luhg_sbsp,
     &     luvec_sbsp, lumv_sbsp, lures,
     &     lu_intstp,lu_corgrvf,luscr,luhg,
     &     lulgrvfold, lulsigold, lu_newlstp, lu_corlstp,
     &     lulst_sbsp, lulpst_sbsp, lultpst_sbsp, lulgv_sbsp,
     &     luscr_l,
     &     luosigold, lu_newostp, luamp_c,
     &     luost_sbsp, luopst_sbsp, luotpst_sbsp, luogv_sbsp,
     &     luscr_o,lukappa_sum,
     &     kumat,kbmat,ksbscr,
     &     kumat_l,kbmat_l,
     &     khred,kgred,kcred,kscr1,kscr2,
     &     ngvtask, nsttask, nhgtask,
     &     igvtask(mxptsk), isttask(mxptsk), ihgtask(mxptsk),
     &     nstdim,npstdim,ntpstdim,ngvdim,nhgdim,maxsbsp,
     &     ndiisdim, ndiisdim_last,nsbspjadim, nsbspjadim_last,
     &     nlstdim,nlpstdim,nltpstdim,nlgvdim,
     &     nldiisdim, nldiisdim_last,nlsbspjadim, nlsbspjadim_last,
     &     nostdim,nopstdim,notpstdim,nogvdim,
     &     nodiisdim, nodiisdim_last,nosbspjadim, nosbspjadim_last,
     &     ibstep, imicdim, i2nd_mode, itrnsym
      real(8), save ::
     &     ener_old, trrad, xngrd, xngrd_tot, xngrd_tot_old,
     &     xnstp, de, de_pred, alpha_last, gamma,
     &     xdamp, xdamp_old, xmicthr, crate, crate_prev,
     &     xdamp_l, xdamp_l_old
      real(8) ::
     &     facs(mxptsk)

      logical ::
     &     laccept, lin_dep, ldamp,
     &     lgrconv, lstconv, ldeconv, lconv, lexit,
     &     llgrconv, llstconv,
     &     logrconv, lostconv

* some external functions
      real*8 ::
     &     inprdd

*
      lblk = -1
      iprintl = max(ntest,iprint)  


* be verbose?
      if (iprintl.ge.1) then
        write(6,'(/,x,a,/,x,a)')
     &         'Optimization control (orbopt)',
     &         '============================='
      end if
      if (ntest.ge.10) then
        write(6,*) 'entered optcont with:'
        write(6,*) ' imacit, imicit, imicit_tot: ',
     &         imacit, imicit, imicit_tot
        write(6,*) ' itask: ',itask
        write(6,*) ' energy:',energy
        write(6,*) ' nwfpar:',nwfpar
        write(6,*) ' lugrvf,lusig,luamp,lutrvec: ',
     &         lugrvf,lusig,luamp,lutrvec
      end if

* some init
      norbr_tot = norbr(1)
      if (nspin.eq.2) norbr_tot = norbr_tot + norbr(2)

*======================================================================*
* check iteration number:
* zero --> init everything
*======================================================================*
      if (imacit.eq.0) then

        if (ntest.eq.10) then
          write(6,*) 'Initialization step entered'
        end if

* print some initial information:
        if (iprintl.ge.1) then
          if (ivar.eq.0) then
            write(6,*) 'Non-variational functional optimization'
          else if (ivar.eq.1) then
            write(6,*) 'Variational functional optimization'
          else
            write(6,*) 'illegal value for ivar'
            stop 'optcont: ivar'
          end if
          if (ilin.eq.0) then
          else if (ilin.eq.1) then
            write(6,*) 'Solve linear equations'
          else
            write(6,*) 'illegal value for ilin'
            stop 'optcont: ilin'
          end if
          if (iorder.eq.1) then
            write(6,*) 'First-order method'
          else if (iorder.eq.2) then
            write(6,*) 'Second-order method'
          else
            write(6,*) 'illegal value for iorder'
            stop 'optcont: iorder'
          end if
          if (iprecnd.eq.0) then
            write(6,*) 'No preconditioner'
          else if(iprecnd.eq.1) then
            write(6,*) 'Diagonal preconditioner'
          else if (iprecnd.eq.2) then
            write(6,*) 'Subspace Jacobian'
          else
            write(6,*) 'illegal value for iprecnd'
            stop 'optcont: iprecnd'
          end if
          if (isubsp.eq.0) then
            write(6,*) 'No subspace'
          else if (isubsp.eq.1) then
            write(6,*) 'Conjugate gradient correction'
          else if (isubsp.eq.2) then
            write(6,*) 'DIIS extrapolation'
            if (idiistyp.eq.1) then
              write(6,*) ' error vector: amplitude diff.s'
            else if (idiistyp.eq.2) then
              write(6,*) ' error vector: preconditioned residual'
            else if (idiistyp.eq.3) then
              write(6,*) ' error vector: residual'
            end if
          else
            write(6,*) 'illegal value for isubsp'
            stop 'optcont: isubsp'
          end if
          if (ilsrch.eq.0) then
            write(6,*) 'Linesearch: alpha = 1'
          else if (ilsrch.eq.1) then
            write(6,*) 'Linesearch: alpha est. from diagonal'
          else if (ilsrch.eq.2) then
            write(6,*) 'Linesearch: additional E calculation'
          else
            write(6,*) 'illegal value for ilsrch'
            stop 'optcont: ilsrch'
          end if
        end if ! iprint

* internal unit numbers
        lugrvfold = iopen_nus('OLDGRD')
        lusigold  = iopen_nus('OLDSIG')
        lu_newstp  = iopen_nus('NEWSTP')
        lu_corstp  = iopen_nus('CORSTP') 
        luscr     = iopen_nus('OPTSCR')
        if (iabs(iorbopt).ge.2) then
          lu_newlstp  = iopen_nus('NEWSTP_L')
        end if
        if (iabs(iorbopt).ge.3) then
          lulgrvfold = iopen_nus('OLDGRD_L')
          lulsigold  = iopen_nus('OLDSIG_L')
          lu_corlstp  = iopen_nus('CORSTP_L')
          luscr_l     = iopen_nus('OPTSCR_L')
        end if
        if (iabs(iorbopt).ge.1) then
          lugrvf_c = iopen_nus('GRDCMB')
          ludia_c   = iopen_nus('DIACMD')
          luosigold  = iopen_nus('OLDSIG_O')
          lu_newostp  = iopen_nus('NEWSTP_O')
          luamp_c  = iopen_nus('AMPCMB')
          luscr_o     = iopen_nus('OPTSCR_O')
        end if
        if (iprecnd.eq.2.or.
     &      isubsp.eq.1 .or.
     &      (isubsp.eq.2.and.idiistyp.eq.1).or.
     &      (isubsp.eq.2.and.idiistyp.eq.4)) then
          lust_sbsp = iopen_nus('STPSBSP')
          if (iabs(iorbopt).ge.3) lulst_sbsp = iopen_nus('STPSBSP_L')
        end if
        if (isubsp.eq.2.and.idiistyp.ge.2) then
          lupst_sbsp = iopen_nus('PSTPSBSP')
          if (iabs(iorbopt).ge.3) lulpst_sbsp = iopen_nus('PSTPSBSP_L')
        end if
        if (isubsp.eq.2.and.idiistyp.eq.2.or.idiistyp.eq.3) then
          lutpst_sbsp = iopen_nus('TPSTPSBSP')
          if (iabs(iorbopt).ge.3)
     &         lultpst_sbsp = iopen_nus('TPSTPSBSP_L')
        end if
        if (iprecnd.eq.2) then
          lugv_sbsp = iopen_nus('GRDSBSP')
          if (iabs(iorbopt).ge.3) lulgv_sbsp = iopen_nus('GRDSBSP_L')
        end if
        if (iprecnd.eq.2.and.isbspjatyp.gt.10) then
          stop 'this route is not active in Lamdba edition'
        end if

        if (isubsp.eq.2.and.idiistyp.eq.4) then
          stop 'this route is not active in Lambda edition'
        end if
        
        if (iorder.eq.2.and.(iabs(iorbopt).ne.2!.or.ikapmod.ne.0
     &         )) then
          write(6,*) 'one or both conditions not met for iorder.eq.2:'
          write(6,*) ' iorbopt = ',iorbopt,' required: +2/-2'
          write(6,*) ' ikapmod = ',ikapmod,' required: 0'
          stop 'check options for optcont_orbopt'
        end if

        if (iorder.eq.2) then
          luvec_sbsp = iopen_nus('SO_V_SBSP')
          lumv_sbsp = iopen_nus('SO_MV_SBSP')
          lures = iopen_nus('SO_RESID')
          lusig_c = iopen_nus('SIGCMB')
          lutrv_c = iopen_nus('TRVCMB')
        end if

        idum = 0
        call memman(idum,idum,'MARK  ',idum,'OPTCON')

* we are still working on it:
        iconv = 0

* reset task list for subspace manager
        nsttask = 0
        ngvtask = 0
        npstdim = 0
        ntpstdim = 0
        nstdim = 0
        ngvdim = 0
        nlpstdim = 0
        nltpstdim = 0
        nlstdim = 0
        nlgvdim = 0
        nopstdim = 0
        notpstdim = 0
        nostdim = 0
        nogvdim = 0

        ndiisdim = 0
        ndiisdim_last = 0
        nldiisdim = 0
        nldiisdim_last = 0
        nodiisdim = 0
        nodiisdim_last = 0

        nsbspjadim = 0
        nsbspjadim_last = 0
        nlsbspjadim = 0
        nlsbspjadim_last = 0
        nosbspjadim = 0
        nosbspjadim_last = 0

        xdamp = 0d0
        xdamp_old = 0d0
        xdamp_l = 0d0
        xdamp_l_old = 0d0
        xdamp_o = 0d0
        xdamp_o_old = 0d0

* work space for subspace jacobian?
        maxsbsp = 0
        if (iprecnd.eq.2) then
          maxsbsp = mxsp_sbspja
          lscr_sbspja = maxsbsp**2 
          lenu= maxsbsp**2
          call memman(kumat,lenu   ,'ADDL  ',2,'U SBJA')
          if (iabs(iorbopt).ge.3)
     &         call memman(kumat_l,lenu   ,'ADDL  ',2,'ULSBJA')
        else
          lscr_sbspja = 0        
        end if

        if (isubsp.eq.2) then
* work space for DIIS?
          maxsbsp = max(mxsp_diis-1,mxsp_sbspja)
          mxdim = mxsp_diis+1
          lenb = mxdim*(mxdim+1)/2
          lscr_diis = mxdim*(mxdim+1)/2
          call memman(kbmat,lenb,'ADDL  ',2,'BMDIIS')
          if (iabs(iorbopt).ge.3)
     &         call memman(kbmat_l,lenb,'ADDL  ',2,'BMLDIS')
        else
          lscr_diis = 0
        end if

        if (iorder.eq.2) then
          lenh = (2*micifac)**2
          lenhex = (2*micifac+1)**2
          call memman(khred,lenh,'ADDL  ',2,'HREDSO')
          call memman(kscr1,lenhex,'ADDL  ',2,'SCR1SO')
          call memman(kscr2,lenhex,'ADDL  ',2,'SCR2SO')
          leng = 2*micifac
          call memman(kcred,leng,'ADDL  ',2,'GREDSO')
          call memman(kgred,leng,'ADDL  ',2,'GREDSO')
        end if

* for conj. grad.s we need a subspace of one (last step)
        if (isubsp.eq.1) then
          stop 'this route is deactivated for Lambda opt'
          maxsbsp = max(maxsbsp,1)
        end if

* init trust radius:
        trrad = trini

        lscr = max(lscr_sbspja,lscr_diis)
        if (lscr.gt.0) then
          call memman(ksbscr,lscr,'ADDL  ',2,'SB SCR')
        end if

* init 2nd-order solver
        if (iorder.eq.2) then
          ! begin with Newton-eigenvector method ...
          i2nd_mode = 2
          ! ... and a gamma of 1
          gamma = 1d0 
        end if

        if (ilin.eq.1) then
          stop 'ilin is not active for lambda equations!'
        end if ! ilin.eq.1

* set itask -- we want the energy in any way
        itask = 1
* ... and the gradient/vector function
        itask = itask + 2

        imacit = 1

* return and let the slaves work
        return

      end if ! init-part

      if (imicit.eq.0) then
*======================================================================*
* beginning of a macro iteration
*======================================================================*
        if (ntest.ge.10) then
          write(6,*) 'macro iteration part entered'
        end if

        if (iabs(iorbopt).eq.1.or.iabs(iorbopt).eq.3) then
          ! combine omega and grad_kappa
          call cmbamp(11,lugrvf,luogrd,0,lugrvf_c,
     &                vec1,nwfpar,norbr_tot,0)
          ! combine diagonals
          call cmbamp(11,ludia,luodia,0,ludia_c,
     &                vec1,nwfpar,norbr_tot,0)
          call cmbamp(11,luamp,lukappa,0,luamp_c,
     &         vec1,nwfpar,norbr_tot,0)
        else if (iabs(iorbopt).eq.2) then
          ! combine omega, grad_T and grad_kappa
          call cmbamp(11,lugrvf,lulgrvf,luogrd,lugrvf_c,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          ! combine diagonals
          call cmbamp(11,ludia,ludia,luodia,ludia_c,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          call cmbamp(11,luamp,luamp,lukappa,luamp_c,
     &         vec1,nwfpar,nwfpar,norbr_tot)
        else
          lugrvf_c = lugrvf
          ludia_c  = ludia
          luamp_c = luamp
        end if


* subspace method needing delta gradient/delta vecfun?
* --->     call subspace manager
        if (iprecnd.eq.2.and.imacit.gt.1) then
          igvtask(ngvtask+1) = 1
          igvtask(ngvtask+2) = 1
          facs(1) = 1d0
          ngvtask = ngvtask + 2
          
          ndel = max(0,ngvdim - 1 - nsbspjadim)
          if (iabs(iorbopt).ge.3) then
            ndel2 = max(0,nlgvdim - 1 - nlsbspjadim)
            ndel = min(ndel,ndel2)
          end if
          if (ndel.gt.0) then
            igvtask(ngvtask+1) = 2
            igvtask(ngvtask+2) = ndel
            ngvtask = ngvtask + 2
          end if

          if (ntest.ge.100)
     &         write(6,*) 'calling sbspman for diff. t-gradient'
            
          idiff = -1
          call optc_sbspman(lugv_sbsp,lugrvf_c,facs,
     &         lugrvfold,ngvdim,maxsbsp,
     &         igvtask,ngvtask,idiff,ndel_recent_gv,
     &         vec1,vec2)

          if (iabs(iorbopt).ge.3) then
            if (ntest.ge.100)
     &           write(6,*) 'calling sbspman for diff. tbar-gradient'

            call optc_sbspman(lulgv_sbsp,lulgrvf,facs,
     &           lulgrvfold,nlgvdim,maxsbsp,
     &           igvtask,ngvtask,idiff,ndel_recent_gv,
     &           vec1,vec2)
          end if

          ngvtask = 0
        end if
      
*----------------------------------------------------------------------*
* check trust radius criterion
*----------------------------------------------------------------------*
        if (iprintl.ge.1) then
          write(6,*) ' current trust radius: ',trrad
        end if

        if (nrdvec.gt.0)
     &         call optc_prjout(nrdvec,lurdvec,lugrvf,
     &                          vec1,vec2,nwfpar,.false.)

        xngrd   = sqrt(inprdd(vec1,vec2,lugrvf,lugrvf,1,lblk))
        xngrd_l = 0d0
        xngrd_o = 0d0
        if (iabs(iorbopt).ge.2)
     &       xngrd_l = sqrt(inprdd(vec1,vec2,lulgrvf,lulgrvf,1,lblk))
        if (iabs(iorbopt).ge.1) then
          xngrd_o = sqrt(inprdd(vec1,vec2,luogrd,luogrd,1,lblk))
        end if
        xngrd_tot_old = xngrd_tot
        xngrd_tot = sqrt(xngrd*xngrd+xngrd_l*xngrd_l+xngrd_o*xngrd_o)
        de = energy - ener_old
        ener_old = energy
        if (ivar.eq.1.and.imacit.gt.1) then 
          if (abs(de_pred).gt.1d-100) then
            rat = de / de_pred
          else
            rat = 1d100*sign(de_pred,1d0)*sign(de,1d0)
          end if
          if (rat.gt.trthr1l.and.rat.lt.trthr1u) then
            trrad = trrad * trfac1
          else if ((rat.gt.trthr2l.and.rat.le.trthr1l).or.
     &           (rat.ge.trthr1u.and.rat.lt.trthr2u)) then
            trrad = trrad * trfac2
          else
            trrad = trrad * trfac3
            ! trigger maybe line-search backsteps
          end if
          trrad = min(trrad,trmax)
          trrad = max(trrad,trmin)
          if (iprintl.ge.2) then
            write(6,*) ' delta E predicted: ',de_pred
            write(6,*) ' delta E observed:  ',de
            write(6,*) ' ratio = ',rat
            write(6,*) ' updated trust radius = ',trrad
          end if
        else if (ivar.eq.0.and.imacit.gt.1) then
          crate_prev = crate
          crate = xngrd_tot_old/xngrd_tot
          if (imacit.gt.2) then
            ratio = crate/crate_prev
          else
            ratio = 1d0
          end if
          if (crate.ge.1.5d0.and.(imacit.eq.1.or.crate_prev.ge.1d0)
     &         .and.ratio.gt.trthr1l) then
            trrad = trrad * trfac1
          else if (crate.ge.1d0.and.
     &           (ratio.gt.trthr2l)) then
            trrad = trrad
          else if (crate.ge.0.90d0) then
            trrad = trrad * trfac2
          else
            trrad = trrad * trfac2**2
c commented out for the moment
c            xngrd = xngrd_old
c            laccept = .false.
          end if
          trrad = min(trrad,trmax)
          trrad = max(trrad,trmin)
          if (iprintl.ge.2) then
            write(6,'(x,a,E10.2)')
     &           '#  current convergence rate: ',crate
            write(6,'(x,a,E10.2)')
     &           '# previous convergence rate: ',crate_prev
            write(6,'(x,a,E10.2)')
     &           '#                     ratio: ',ratio
            write(6,*) '# updated trust radius = ',trrad
c            if (.not.laccept)
c     &           write(6,*) '# backstepping triggered'
          end if

        end if

*----------------------------------------------------------------------*
* get a new step direction
*  a) taking directly the gradient
*  b) taking a diagonal preconditioner (i.e. approx. diagonal Hessian)
*  c) use a subspace Hessian/Jacobian
*----------------------------------------------------------------------*
        if (iprecnd.eq.0) then
          stop 'iprecnd.eq.0 inactive for Lambda edition'
        else if (iprecnd.eq.1) then
*   use diagonal preconditioner/quasi-Newton step with diag. Hess.
          if (iabs(iorbopt).ge.3) 
     &      call optc_diagp(lulgrvf,ludia,trrad,
     &                      lu_newlstp,xdamp_l,.true.,vec1,vec2,nwfpar)
          call optc_diagp(lugrvf_c,ludia_c,trrad,
     &                    lu_newstp,xdamp,.true.,vec1,vec2,nwfpar)
        else if (iprecnd.eq.2) then
          xdamp_old = xdamp
          if (iabs(iorbopt).ge.3) xdamp_l_old = xdamp_l
          if (iabs(iorbopt).ge.3) lulstp = lu_newlstp
          lustp = lu_newstp
          if (imacit.ge.isbspja_start) lulstp = luscr_l
          if (imacit.ge.isbspja_start) lustp = luscr
c new:
          if (iabs(iorbopt).ge.3) 
     &       call optc_diagp(lulgrvf,ludia,trrad,
     &           lulstp,xdamp_l,.true.,vec1,vec2,nwfpar)
          call optc_diagp(lugrvf_c,ludia_c,trrad,
     &         lustp,xdamp,.true.,vec1,vec2,nwfpar)
c old:
c          if (iabs(iorbopt).ge.3) 
c     &       call optc_diagp(lulgrvf,ludia,trrad,
c     &           lulstp,xdamp_l,.true.,vec1,vec2,nwfpar)
c          call optc_diagp(lugrvf_c,ludia_c,trrad,
c     &         lustp,xdamp,.true.,vec1,vec2,nwfpar)
          if (imacit.ge.isbspja_start) then
*   use (approximate) subspace Hessian/Jacobian
            if (iabs(iorbopt).ge.3) then
              nadd = 1
              nlsbspjadim = min(nlsbspjadim + nadd,mxsp_sbspja)
             call optc_sbspja_new(isbspjatyp,thr_sbspja,nlstdim,nlgvdim,
     &           nlsbspjadim,mxsp_sbspja,
     &           nadd,nlsbspjadim_last,
     &           lulgrvf,ludia,trrad,lulstp,lu_newlstp,
     &           xdamp_l,xdamp_l_old,
     &           lulst_sbsp,lulgv_sbsp,
     &           work(kumat_l),work(ksbscr),vec1,vec2,iprint)
            end if

            nadd = 1
            nsbspjadim = min(nsbspjadim + nadd,mxsp_sbspja)
            call optc_sbspja_new(isbspjatyp,thr_sbspja,nstdim,ngvdim,
     &           nsbspjadim,mxsp_sbspja,
     &           nadd,nsbspjadim_last,
     &           lugrvf_c,ludia_c,trrad,lustp,lu_newstp,
     &           xdamp,xdamp_old,
     &           lust_sbsp,lugv_sbsp,
     &           work(kumat),work(ksbscr),vec1,vec2,iprint)

          end if
        end if ! iprecnd

c        xnrm = sqrt(inprdd(vec1,vec1,lukappa,lukappa,1,lblk))

c V A: use kap2u here
c V B: skip it
        if (iorder.lt.2) then
        ! unless we had the real gradient wrt kappa at kappa!=0
        ! we have to do this trick (i.pt. for Brueckner CC):
        ! (described in Hampel et al., CPL something etc etc)
         if ((iorbopt.eq.1.or.iorbopt.eq.3).and.ikapmod.eq.1) then

          ! calc.  dkappa = ln(exp(kappa)exp(dkappa_cep))-kappa
          call cmbamp(01,luscr,luscr_o,0,lu_newstp,
     &                vec1,nwfpar,norbr_tot,0)
          call kap2u(4,luscr_o,lukappa,lutrmat,icmpr,norbr,nspin)

          call cmbamp(11,luscr,luscr_o,0,lu_newstp,
     &                vec1,nwfpar,norbr_tot,0)

         else if (iorbopt.eq.2.and.ikapmod.eq.1) then

          luscr2 = iopen_nus('OPTSCR_TMP')
          ! calc.  dkappa = ln(exp(kappa)exp(dkappa_cep))-kappa
          call cmbamp(01,luscr,luscr2,luscr_o,lu_newstp,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          call kap2u(4,luscr_o,lukappa,lutrmat,icmpr,norbr,nspin)

          call cmbamp(11,luscr,luscr2,luscr_o,lu_newstp,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          call relunit(luscr2,'delete')

         end if
        end if

        lustp = lu_newstp
        if (iabs(iorbopt).ge.3) lulstp = lu_newlstp
        if (iabs(iorbopt).ge.1) then
          luostp = lu_newostp
          luostp2 = luogrd
        end if
        if (isubsp.eq.1) then
          stop 'isubsp.eq.1 not active for Lambda edition'
        else if (isubsp.eq.2) then
* use DIIS?
          if (imacit.ge.idiis_start) then
c              nadd = 1
c              ndiisdim = min(ndiisdim + nadd,mxsp_diis)
              if (idiistyp.eq.1) then
                stop 'idiistyp.eq.1 not active for Lambda edition'
              else if (idiistyp.eq.2.or.idiistyp.eq.3) then
c                ludiis_err = lupst_sbsp
c                ludiis_bas = lutpst_sbsp
c                nsbspdim = ntpstdim
              else
                write(6,*) 'unknown diistyp = ',idiistyp
                stop 'optcont'
              end if
              if (iabs(iorbopt).ge.3) then
                nadd = 1
                nldiisdim = min(nldiisdim + nadd,mxsp_diis)
                ludiis_err = lulpst_sbsp
                ludiis_bas = lultpst_sbsp
                nsbspdim = nltpstdim
                call optc_diis(idiistyp,thr_diis,nsbspdim,
     &               nldiisdim,mxsp_diis,
     &               nadd,nldiisdim_last,1d0,
     &               lu_newlstp,lu_corlstp,ludiis_err,ludiis_bas,
     &               lulamp,lulgrvf,
     &               work(kbmat_l),work(ksbscr),vec1,vec2,iprint)
              end if

              nadd = 1
              ndiisdim = min(ndiisdim + nadd,mxsp_diis)
              ludiis_err = lupst_sbsp
              ludiis_bas = lutpst_sbsp
              nsbspdim = ntpstdim

              call optc_diis(idiistyp,thr_diis,nsbspdim,
     &             ndiisdim,mxsp_diis,
     &             nadd,ndiisdim_last,1d0,
     &             lu_newstp,lu_corstp,ludiis_err,ludiis_bas,
     &             luamp_c,lugrvf_c,
     &             work(kbmat),work(ksbscr),vec1,vec2,iprint)

              lustp = lu_corstp
              if (iabs(iorbopt).ge.3) lulstp = lu_corlstp
          end if ! imacit.ge.idiis_start

        else ! isubsp.eq.0
          lustp = lu_newstp
          if (iabs(iorbopt).ge.3) lulstp = lu_newlstp
          if (iabs(iorbopt).ge.1) luostp = lu_newostp
        end if

        if (iabs(iorbopt).eq.1.or.iabs(iorbopt).eq.3) then
          ! decompose
          call cmbamp(01,luscr_o,luostp,0,lustp,
     &                vec1,nwfpar,norbr_tot,0)
          lustp_c = lustp
          lustp = luscr_o

        else if (iabs(iorbopt).eq.2) then

          ! decompose
          call cmbamp(01,luscr_o,lu_newlstp,luostp,lustp,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          lulstp = lu_newlstp
          lustp_c = lustp
          lustp = luscr_o

        else
          lustp_c = lustp
        end if

* norm of unscaled step:
        xnstp = sqrt(inprdd(vec1,vec2,lustp,lustp,1,lblk))
        if (iabs(iorbopt).ge.2)
     &       xnstp_l = sqrt(inprdd(vec1,vec2,lulstp,lulstp,1,lblk))
        if (iabs(iorbopt).ge.1) 
     &      xnstp_o = sqrt(inprdd(vec1,vec2,luostp,luostp,1,lblk))

        if (ntest.ge.20) then
          write(6,*) 'New unscaled step length: ',xnstp,xnstp_l
          write(6,*) ' <s|s> was ',xnstp*xnstp,xnstp_l*xnstp_l
        end if

        cont1   = inprdd(vec1,vec2,lustp,lugrvf,1,lblk)
c        de_pred = 0.5d0*cont1
        de_pred = cont1
        if (iabs(iorbopt).ge.2) then
          cont2   = inprdd(vec1,vec2,lulstp,lulgrvf,1,lblk)
c          de_pred = de_pred+0.5d0*cont2
          de_pred = de_pred+cont2
        end if
        if (iabs(iorbopt).ge.1) then
          cont3 = inprdd(vec1,vec2,luostp,luogrd,1,lblk)
          de_pred = de_pred + cont3
        end if

        if (ilsrch.eq.2) then
          stop 'ilsrch.eq.2 inactive in Lambda edition'
        else
          if (iorder.eq.1) then
* 1st order method: increase macroiteration counter and request
* energy and gradient/vectorfunction
            itask = 1+2
            imicit = 0
          else if (iorder.eq.2) then
* 2nd order method: increase microiteration counter and request
* matrix-vector product
            itask = 4
            imicit = 1
            imicdim = 0
            imicit_tot = imicit_tot + 1
            itrnsym = 0
          end if
        end if

* obtain new paramter set |X> = |Xold> + alpha |d>
*  macro iteration --> see end-of-macro-iteration section
*  micro iteration:
        if (imicit.eq.1.and.ilsrch.eq.2) then
          stop 'unexpected event'
        else if (imicit.eq.1) then
          ! normalize trial vector
          xnrm = sqrt(inprdd(vec1,vec2,lustp,lustp,1,lblk)
     &               +inprdd(vec1,vec2,lulstp,lulstp,1,lblk)
     &               +inprdd(vec1,vec2,luostp,luostp,1,lblk))
          ! and let it point into the direction of the gradient again
          call sclvcd(lustp,lutrvec,-1d0/xnrm,vec1,1,lblk)
          call sclvcd(lulstp,lutrv_l,-1d0/xnrm,vec1,1,lblk)
          call sclvcd(luostp,lutrv_o,-1d0/xnrm,vec1,1,lblk)
          ! set threshold for microiterations
          ! hard: 1d-4
          ! weak: 1d1
          xmicthr = max(
     &        inprdd(vec1,vec2,lugrvf_c,lugrvf_c,1,lblk),!**(3d0/2d0),
     &                  thrgrd)
c          xmicthr = min(xmicthr,0.1d0)
          xmicthr = min(xmicthr,0.5d0)
          write(6,'(" >>",x,a,e12.6)')
     &         'micro-iterations started -- threshold = ',xmicthr

        end if

* save old gradient for next iteration
        if (isubsp.gt.0.or.iprecnd.ge.2) then 
          call copvcd(lugrvf_c,lugrvfold,vec1,1,lblk)
          if (iabs(iorbopt).ge.3)
     &         call copvcd(lulgrvf,lulgrvfold,vec1,1,lblk)
        end if

      else
*======================================================================*
* micro-iteration
*======================================================================*
        if (ntest.ge.10) then
          write(6,*) '*** micro iteration part entered ***'
        end if

        if (iorder.eq.2) then

          ! combine sigma-vectors
          call cmbamp(11,lusig,lusig_l,lusig_o,lusig_c,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          ! recombine trial-vectors
          call cmbamp(11,lutrvec,lutrv_l,lutrv_o,lutrv_c,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          ! recombine grad-vectors
          call cmbamp(11,lugrvf,lulgrvf,luogrd,lugrvf_c,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          ! combine diagonals
          call cmbamp(11,ludia,ludia,luodia,ludia_c,
     &                vec1,nwfpar,nwfpar,norbr_tot)
          call cmbamp(11,luamp,luamp,lukappa,luamp_c,
     &         vec1,nwfpar,nwfpar,norbr_tot)

          imicdim = imicdim + 1
          maxit = 2*micifac
          fac = 1d0
          isttask(1) = 1
          isttask(2) = 1
          nsttask = 2
          ! push trial vector onto subspace
          call optc_sbspman(luvec_sbsp,lutrv_c,fac,ludum,imicdim-1,
     &                                                         maxit,
     &           isttask,nsttask,0,ndum,vec1,vec2)
          ! push mv-product onto subspace
          call optc_sbspman(lumv_sbsp,lusig_c,fac,ludum,imicdim-1,maxit,
     &           isttask,nsttask,0,ndum,vec1,vec2)
          nsttask = 0

          ! update current reduced Hessian
          isymm = ivar ! enforce symmetric Hessian
     &                 ! (unless it is a Jacobian as for ivar==0)
          call optc_redh(isymm,work(khred),imicdim,imicdim-1,
     &                   lutrv_c,lusig_c,luvec_sbsp,lumv_sbsp,
     &                   vec1,vec2)

          ! set up current reduced gradient
          call optc_sbspja_prjpstp(work(kgred),imicdim,0,
     &                             lugrvf_c,luvec_sbsp,
     &                             vec1,vec2,.false.)
          
          ! loop over the next routine call to enable mode-switches
          do int = 1, 2
            ! call TR-Newton to 
            !  - get new reduced expansion coeff.s and damping
            call optc_trnewton(i2nd_mode,iret,itrnsym,
     &           work(khred),work(kgred),
     &           work(kcred),work(kscr1),work(kscr2),
     &           imicdim,xlamb,gamma,trrad,de_pred)
          
            if (iret.ne.0.and.int.eq.2) then
              ! still not satisfied?
              write(6,*) 'I hoped this never happens ....'
              stop 'wrong hopes'
            end if

            if (iret.ne.0.and.i2nd_mode.eq.2) then
              ! switch to newton step
              if (ntest.ge.10) write(6,*) '>> switched to Newton step'
              xlamb = 0d0
              i2nd_mode = 1
            else if (iret.ne.0.and.i2nd_mode.eq.1) then
              ! switch to newton-eigenvector procedure
              if (ntest.ge.10) write(6,*) '>> switched to Newton EV'
              gamma = 1d0
              i2nd_mode = 2
            else
              exit ! exit the int=1,2 loop
            end if
          end do ! int = 1,2

          !  - assemble new residual in full space
          !  - get residual norm
          call optc_trn_resid(work(kcred),imicdim,xlamb,
     &         lures,xresnrm,
     &         luvec_sbsp,lumv_sbsp,lugrvf_c,
     &         vec1,vec2)

          if (nrdvec.gt.0) then
            call optc_prjout(nrdvec,lurdvec,lures,
     &                       vec1,vec2,nwfpar,.false.)
            xresnrm = sqrt(inprdd(vec1,vec1,lures,lures,1,-1))
          end if


          if (xresnrm.lt.xmicthr.or.imicit.ge.maxit) then
            if (imicit.ge.maxit) then
              write(6,*) '>> WARNING: No convergence in microiterations'
              ! ... but we go on anyway
            end if

            ! assemble new step
            lustp_c = lu_corstp
            call mvcsmd(luvec_sbsp,work(kcred),lu_corstp,luscr,
     &           vec1,vec2,imicdim,1,lblk)
            xnstp = sqrt(inprdd(vec1,vec1,lu_corstp,lu_corstp,1,lblk))
            alpha_eff = 1d0
            alpha = 1d0

            ! unpack step
            lustp = lu_newstp
            lulstp = lu_newlstp
            luostp = lu_newostp
            call cmbamp(01,lustp,lulstp,luostp,lu_corstp,
     &           vec1,nwfpar,nwfpar,norbr_tot)

c V A: 
c            nix
c V B: 
            if (ikapmod.ne.0)
     &           call kap2u(4,luostp,lukappa,lutrmat,icmpr,norbr,nspin)

            imicit = 0
            imicdim = 0
            itask  = 1 + 2
          else
            ! new raw step by preconditioning
            call dmtvcd2(vec1,vec2,ludia_c,lures,lutrv_c,
     &           -1d0,-xlamb,1,1,lblk)

            ! orthogonalize against previous space and normalize;
            call optc_orthvec(work(kscr1),work(kscr2),
     &           imicdim,1,lin_dep,
     &           luvec_sbsp,lutrv_c,
     &           vec1,vec2)

            if (lin_dep) then
              if (iprintl.ge.2) then
                write(6,*)
     &               '>> linear dependency in subspace; restarting ...'
              end if
              ! combine all previous vectors to give the best solution
              call mvcsmd(luvec_sbsp,work(kcred),lu_newstp,luscr,
     &             vec1,vec2,imicdim,1,lblk)
              ! ... and replace
              isttask(1) = 2       ! delete ...
              isttask(2) = imicdim ! ... whole space
              isttask(3) = 1
              isttask(4) = 1
              nsttask = 4
              call optc_sbspman(luvec_sbsp,lu_newstp,fac,ludum,imicdim,
     &                                                         maxit,
     &           isttask,nsttask,0,ndum,vec1,vec2)

              ! combine all previous matrix-vector products 
              call mvcsmd(lumv_sbsp,work(kcred),lu_newstp,luscr,
     &             vec1,vec2,imicdim,1,lblk)
              ! ... and put combined matrix-vector onto subspace
              call optc_sbspman(lumv_sbsp,lu_newstp,fac,ludum,imicdim,
     &                                                         maxit,
     &           isttask,nsttask,0,ndum,vec1,vec2)
              nsttask = 0
              imicdim = 1
    
              ! orthogonalize against the combined vector
              call optc_orthvec(work(kscr1),work(kscr2),
     &           imicdim,1,lin_dep,
     &           luvec_sbsp,lutrv_c,
     &           vec1,vec2)

              if (lin_dep) then
                write(6,*) ' unresolvable linear dependency problem!'
                stop 'optcont (2nd order)'
              end if

            end if

            ! unpack trial-vectors for next micro-iteration
            call cmbamp(01,lutrvec,lutrv_l,lutrv_o,lutrv_c,
     &                vec1,nwfpar,nwfpar,norbr_tot)

            imicit = imicit + 1
            imicit_tot = imicit_tot + 1
            itask  = 4

          end if
        else
          write(6,*) 'unexpected event in micro-iteration section'
          stop 'optcont'
        end if

      end if ! macro/micro iteration switch


      lexit = .false.
      lconv = .false.
      ! convergence check:
      !   end of iteration for 1st-order methods
      !   before first micro-iteration for 2nd-order methods
      if ((iorder.eq.1.and.imicit.eq.0) .or.
     &    (iorder.eq.2.and.imicit.eq.1)) then

*======================================================================*
* end of macro-iteration (indicated by imicit.eq.0):
*======================================================================*

*----------------------------------------------------------------------*
*  check convergence and max. iterations:
*----------------------------------------------------------------------*
        lexit = .false.

        lstconv = (iorder.eq.2).or.xnstp.lt.thrstp
        lgrconv = xngrd.lt.thrgrd
        ldeconv = (iorder.eq.2).or.abs(de).lt.thr_de

        lconv = lstconv .and. lgrconv .and. ldeconv 

        if (iabs(iorbopt).ge.2) then
          llstconv = (iorder.eq.2).or.xnstp_l.lt.thrstp
          llgrconv = xngrd_l.lt.thrgrd
          lconv = lconv .and. llstconv .and. llgrconv
        end if

        if (iabs(iorbopt).ge.1) then
          lostconv = (iorder.eq.2).or.xnstp_o.lt.thrstp
          logrconv = xngrd_o.lt.thrgrd
          lconv = lconv .and. lostconv .and. logrconv
          if (lostconv.and.logrconv) then
c            itransf = 0
            itransf = 1
          else
            itransf = 1
          end if
        end if
        if (iorder.eq.2.and.imicit.eq.1) itransf = 0

        if (iprintl.ge.1) then
          if (iorder.eq.1)
     &         write(6,*) 'after iteration ',imacit
          if (iorder.eq.2)
     &         write(6,*) 'after macro-iteration ',imacit
          if (iorder.ne.2)
     &          write (6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of new t step:  ', xnstp,
     &                           '   threshold:  ', thrstp,
     &                           '   converged:  ', lstconv
          write(6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of t gradient:  ', xngrd,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', lgrconv
          if (iabs(iorbopt).ge.2) then
            if (iorder.ne.2)
     &            write (6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of new l step:  ', xnstp_l,
     &                           '   threshold:  ', thrstp,
     &                           '   converged:  ', llstconv
            write(6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of l gradient:  ', xngrd_l,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', llgrconv
          end if
          if (iabs(iorbopt).ge.1) then
            if (iorder.ne.2)
     &            write (6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of new k step:  ', xnstp_o,
     &                           '   threshold:  ', thrstp,
     &                           '   converged:  ', lostconv
            write(6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of k gradient:  ', xngrd_o,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', logrconv
          end if
          if (iorder.ne.2)
     &          write (6,'(x,2(a,e10.3),a,l)')
     &                    '   change in energy:  ', de,
     &                           '   threshold:  ', thr_de,
     &                           '   converged:  ', ldeconv

          if (lconv.and.iorder.eq.1)
     &         write(6,'(x,a,i5,a)')
     &         'CONVERGED IN ',imacit,' ITERATIONS'
          if (lconv.and.iorder.eq.2) then
            imicit_tot = imicit_tot-1
            write(6,'(x,a,i5,a,i6,a)')
     &         'CONVERGED IN ',imacit,' MACRO-ITERATIONS (',imicit_tot,
     &         ' MICRO-ITERATIONS)'
          end if
          if (lconv) iconv = 1
          if (lconv) imicit = 0
        end if

        if (.not.lconv) imacit = imacit + 1

        if (.not.lconv.and.
     &       (imacit.gt.maxmacit.or.
     &       imicit_tot.gt.maxmicit)) then
          write(6,*) 'NO CONVERGENCE OBTAINED'
          imacit = imacit - 1
          imicit_tot = imicit_tot - 1
          imicit = 0
          lexit = .true.
        end if
      end if

      ! some stuff to be done at the end of the macro-iteration:
      if (imicit.eq.0) then
*----------------------------------------------------------------------*
* clean up
*----------------------------------------------------------------------*
        if (lconv.or.lexit) then

          idum = 0
          call memman(idum,idum,'FLUSM  ',idum,'OPTCON')
          call relunit(lugrvfold,'delete')
          call relunit(lusigold,'delete')
          call relunit(lu_newstp,'delete')
          call relunit(lu_corstp,'delete')
          call relunit(luscr,'delete')
          if (iabs(iorbopt).ge.2) then
            call relunit(lu_newlstp,'delete')
          end if
          if (iabs(iorbopt).ge.3) then
            call relunit(lulgrvfold,'delete')
            call relunit(lulsigold,'delete')
            call relunit(lu_corlstp,'delete')
            call relunit(luscr_l,'delete')
          end if
          if (iprecnd.eq.2.or.isubsp.eq.1.or.
     &       (isubsp.eq.2.and.idiistyp.eq.1).or.
     &        (isubsp.eq.2.and.idiistyp.eq.4)) then
            call relunit(lust_sbsp,'delete')
            if (iabs(iorbopt).ge.3) call relunit(lulst_sbsp,'delete')
          end if
          if (isubsp.eq.2.and.idiistyp.ge.2) then
            call relunit(lupst_sbsp,'delete')
            if (iabs(iorbopt).ge.3) call relunit(lulpst_sbsp,'delete')
          end if
          if (isubsp.eq.2.and.idiistyp.eq.2.or.idiistyp.eq.3) then
            call relunit(lutpst_sbsp,'delete')
            if (iabs(iorbopt).ge.3) call relunit(lultpst_sbsp,'delete')
          end if
          if (iprecnd.eq.2) then
            call relunit(lugv_sbsp,'delete')
            if (iabs(iorbopt).ge.3) call relunit(lulgv_sbsp,'delete')
          end if
          
          if (iabs(iorbopt).ge.1) then
            call relunit(lugrvf_c,'delete')
            call relunit(ludia_c,'delete')
            call relunit(luosigold,'delete')
            call relunit(lu_newostp,'delete')
            call relunit(luamp_c,'delete')
            call relunit(luscr_o,'delete')
          end if

          if (iorder.eq.2) then
            call relunit(luvec_sbsp,'delete')
            call relunit(lumv_sbsp,'delete')
            call relunit(lures,'delete')
            call relunit(lusig_c,'delete')
            call relunit(lutrv_c,'delete')
          end if

          itask = 8 ! stop it

        else ! do some stuff for the next macro-iteration

          if (iorder.eq.2.and.imicit.eq.0) itransf = 1

* subspace method needing the steps? call subspace manager
          if (isubsp.eq.1.or.isubsp.eq.2.or.iprecnd.eq.2) then
            isttask(nsttask+1) = 1
            isttask(nsttask+2) = 1
            nsttask = nsttask + 2
            facs(1) = 1d0
            ndel = 0
! manage here requests of diis and others to delete vectors
ccc currently only diis:
            if (isubsp.eq.2.and.idiistyp.eq.1) then
              ndel = nstdim + 1 - ndiisdim
              if (iabs(iorbopt).ge.3) then
                ndel2 = nlstdim + 1 - nldiisdim
                ndel = min(ndel,ndel2)
              end if
            else if (isubsp.eq.2.and.idiistyp.gt.1) then
              ndel = ntpstdim + 1 - ndiisdim
              if (iabs(iorbopt).ge.3) then
                ndel2 = nltpstdim + 1 - nldiisdim
                ndel = min(ndel,ndel2)
              end if
            end if
ccc and subspace jac.
            if (iprecnd.eq.2) then
              ndel2 = max(0,nstdim - 1 - nsbspjadim)
              ndel = min(ndel,ndel2)
              if (iabs(iorbopt).ge.3) then
                ndel3 = max(0,nlstdim - 1 - nlsbspjadim)
                ndel = min(ndel,ndel3)
              end if
            end if

            if (ndel.gt.0) then
              isttask(nsttask+1) = 2
              isttask(nsttask+2) = ndel
              nsttask = nsttask + 2
            end if

            if (iprecnd.eq.2.or.isubsp.eq.1.or.
     &          (isubsp.eq.2.and.idiistyp.eq.1).or.
     &          (isubsp.eq.2.and.idiistyp.eq.4) ) then
              idiff = 0
              ludum = 0
              if (iabs(iorbopt).ge.3) then
                if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for L-L(last)'
                call optc_sbspman(lulst_sbsp,lulstp,facs,
     &             ludum,nlstdim,maxsbsp,
     &             isttask,nsttask,idiff,ndel_recent_st,
     &             vec1,vec2)
              end if

              if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for T-T(last)'
              call optc_sbspman(lust_sbsp,lustp_c,facs,
     &             ludum,nstdim,maxsbsp,
     &             isttask,nsttask,idiff,ndel_recent_st,
     &             vec1,vec2)
            end if
            if (isubsp.eq.2.and.idiistyp.ge.2) then
c              npstdim = ntpstdim
              if (idiistyp.ne.4) then
                idiff = 1     
                if (iabs(iorbopt).ge.3) then
                  if (ntest.ge.100)
     &               write(6,*) 'calling sbspman for L+dL(pert)'
                  call optc_sbspman(lultpst_sbsp,lulamp,facs,
     &               lu_newlstp,nltpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_tpst,
     &               vec1,vec2)
                end if
                if (ntest.ge.100)
     &               write(6,*) 'calling sbspman for T+dT(pert)'
                call optc_sbspman(lutpst_sbsp,luamp_c,facs,
     &               lu_newstp,ntpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_tpst,
     &               vec1,vec2)
              end if
              if (idiistyp.eq.2) then
                ludum = 0
                idiff = 0 
                if (iabs(iorbopt).ge.3) then
                  if (ntest.ge.100) write(6,*)
     &               'calling sbspman for dL(pert)'
                  call optc_sbspman(lulpst_sbsp,lu_newlstp,facs,
     &               ludum,nlpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
                end if
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for dT(pert)'
                call optc_sbspman(lupst_sbsp,lu_newstp,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
              else if (idiistyp.eq.3.or.idiistyp.eq.4) then
                ludum = 0
                idiff = 0 
                if (iabs(iorbopt).ge.3) then
                  if (ntest.ge.100) write(6,*)
     &               'calling sbspman for L residual'
                  call optc_sbspman(lulpst_sbsp,lulgrvf,facs,
     &               ludum,nlpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
                end if
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for T residual'
                call optc_sbspman(lupst_sbsp,lugrvf_c,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
              end if
            end if
            nsttask = 0


          end if ! isubsp.eq.1.or.isubsp.eq.2.or.iprecnd.eq2

          ! build new amplitudes for next energy and gradient evalutation
          ! the step is expected on lustp the old amplitudes on luamp

* obtain new paramter set |X> = |Xold> + alpha |d>
          if (iabs(iorbopt).ge.2) then
            call vecsmd(vec1,vec2,1d0,1d0,
     &         lulamp,lulstp,luscr,1,lblk)
            call copvcd(luscr,lulamp,vec1,1,lblk)
          end if

          call vecsmd(vec1,vec2,1d0,1d0,
     &         luamp,lustp,luscr,1,lblk)
          call copvcd(luscr,luamp,vec1,1,lblk)
          if (lustp.ne.lutrvec)
     &         call copvcd(lustp,lutrvec,vec1,1,lblk)

          if (iabs(iorbopt).ge.1.and.ikapmod.ne.0) then
            call vecsmd(vec1,vec2,1d0,1d0,
     &           lukappa,luostp,luscr,1,lblk)
            call copvcd(luscr,lukappa,vec1,1,lblk)
            
            call kap2u(3,ludum,lukappa,lutrmat,icmpr,norbr,nspin)
          else if (iabs(iorbopt).ge.1) then
            call kap2u(-3,luostp,ludum,lutrmat,icmpr,norbr,nspin)
          end if

        end if ! "prepare for next macro iteration" part

        if (ntest.ge.10) then
          write(6,*) 'at the end of optcont:'
          write(6,*) ' itask = ',itask
          write(6,*) ' imacit,imicit,imicit_tot: ',
     &         imacit,imicit,imicit_tot
        end if

      end if ! end-of-macro-iteration part

      end
*----------------------------------------------------------------------*
c $Id$
