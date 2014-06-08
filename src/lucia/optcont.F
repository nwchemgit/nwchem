*----------------------------------------------------------------------*
      subroutine optcont(imacit,imicit,imicit_tot,iprint,
     &                   itask,iconv,
     &                   luamp,lutrvec,
     &                   energy,
     &                   vec1,vec2,nwfpar,
     &                   lugrvf,lusig,ludia,
     &                   nrdvec,lurdvec)
*----------------------------------------------------------------------*
*
* General optimization control routine for non-linear optimization.
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
*              should remain connected such that in the next call
*              the trial vector still can be found here
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
* treatment of redundant vectors: if nrdvec.ne.0 a set of vectors
*     spanning the redundant space is expected on lurdvec
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
     &     itask
      integer, intent(inout) ::
     &     imacit, imicit
      integer, intent(in) ::
     &     nwfpar,
     &     luamp, lutrvec, lugrvf, lusig, ludia
      real(8), intent(in) ::
     &     energy
* two scratch vectors
      real(8), intent(inout) ::
     &     vec1(nwfpar), vec2(nwfpar)

* some private static variables:
      integer, save ::
     &     lugrvfold, lusigold, lu_newstp, lu_corstp,
     &     lust_sbsp, lupst_sbsp, lutpst_sbsp, lugv_sbsp, luhg_sbsp,
     &     luvec_sbsp, lumv_sbsp, lures,
     &     lu_intstp,lu_corgrvf,luscr,luhgam,luhg_last,
     &     ku_mat,kb_mat,ksbscr,
     &     khred,kgred,kcred,kscr1,kscr2,
     &     ngvtask, nsttask, nhgtask,
     &     igvtask(mxptsk), isttask(mxptsk), ihgtask(mxptsk),
     &     nstdim,npstdim,ntpstdim,ngvdim,nhgdim,maxsbsp,
     &     ndiisdim, ndiisdim_last,nsbspjadim, nsbspjadim_last,
     &     ibstep, imicdim, i2nd_mode
      real(8), save ::
     &     ener_old, trrad, xngrd, xngrd_old,
     &     xnstp, de, de_pred, alpha_last, gamma,
     &     xdamp, xdamp_old, xmicthr, crate, crate_prev
      real(8) ::
     &     facs(mxptsk)

      logical ::
     &     laccept, lgrconv, lstconv, ldeconv, lconv, lexit, lin_dep,
     &     ldamp

* some external functions
      real*8 ::
     &     inprdd, inprod

*
      lblk = -1
      iprintl = max(ntest,iprint)  


* be verbose?
      if (iprintl.ge.1) then
        write(6,'(/,x,a,/,x,a)')
     &         'Optimization control',
     &         '===================='
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
        if (iprecnd.eq.2.or.
     &      isubsp.eq.1 .or.
     &      (isubsp.eq.2.and.idiistyp.eq.1).or.
     &      (isubsp.eq.2.and.idiistyp.eq.4))
     &       lust_sbsp = iopen_nus('STPSBSP')
        if (isubsp.eq.2.and.idiistyp.ge.2)
     &       lupst_sbsp = iopen_nus('PSTPSBSP')
        if (isubsp.eq.2.and.idiistyp.eq.2.or.idiistyp.eq.3)
     &       lutpst_sbsp = iopen_nus('TPSTPSBSP')
        if (iprecnd.eq.2)
     &       lugv_sbsp = iopen_nus('GRDSBSP')
        if (iprecnd.eq.2.and.isbspjatyp.gt.10)
     &       luhg_sbsp = iopen_nus('HDGSBSP')
        if (iprecnd.eq.2.and.isbspjatyp.gt.10)
     &       luhgam = iopen_nus('OPTHDGAM')
        if (iprecnd.eq.2.and.isbspjatyp.gt.10)
     &       luhg_last = iopen_nus('OPTHGOLD')

        if (isubsp.eq.2.and.idiistyp.eq.4) then
          lu_intstp  = iopen_nus('DIISISTP')
          lu_corgrvf = iopen_nus('DIISCGRD')
        end if

        if (iorder.eq.2) then
          luvec_sbsp = iopen_nus('SO_V_SBSP')
          lumv_sbsp = iopen_nus('SO_MV_SBSP')
          lures = iopen_nus('SO_RESID')
        end if
        
        idum = 0
        call memman(idum,idum,'MARK  ',idum,'OPTCON')

* we are still working on it:
        iconv = 0

* reset task list for subspace manager
        nsttask = 0
        ngvtask = 0
        nhgtask = 0
        npstdim = 0
        ntpstdim = 0
        nstdim = 0
        ngvdim = 0
        nhgdim = 0

        ndiisdim = 0
        ndiisdim_last = 0

        nsbspjadim = 0
        nsbspjadim_last = 0

        xdamp = 0d0
        xdamp_old = 0d0

        ibstep = 0

* work space for subspace jacobian?
        maxsbsp = 0
        if (iprecnd.eq.2) then
          maxsbsp = mxsp_sbspja
          lscr_sbspja = maxsbsp**2
c TEST
c          lscr_sbspja = lscr_sbspja + 4*27*(maxsbsp+2)**2
          lenu= maxsbsp**2
          call memman(ku_mat,lenu   ,'ADDL  ',2,'U SBJA')
        else
          lscr_sbspja = 0        
        end if

        if (isubsp.eq.2) then
* work space for DIIS?
          maxsbsp = max(mxsp_diis-1,mxsp_sbspja)
          mxdim = mxsp_diis+1
          lenb = mxdim*(mxdim+1)/2
          lscr_diis = mxdim*(mxdim+1)/2
          call memman(kb_mat,lenb,'ADDL  ',2,'BMDIIS')
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
* linear equations:
*  the RHS resides on lugrvf; save as old gradient ...
          call copvcd(lugrvf,lugrvfold,vec1,1,lblk)

*  ... and get first trial step by dividing with diagonal on ludia.
          call optc_diagp(lugrvf,ludia,trrad,
     &         luamp,xdamp,.false.,vec1,vec2,nwfpar)

          lustp = luamp

          if (isubsp.eq.1.or.isubsp.eq.2.or.iprecnd.eq.2) then
*  save the trial step as first subspace member

            isttask(nsttask+1) = 1
            isttask(nsttask+2) = 1
            nsttask = nsttask + 2
            facs(1) = 1d0
            ndel = 0

            if (iprecnd.eq.2.or.isubsp.eq.1.or.
     &          (isubsp.eq.2.and.idiistyp.eq.1).or.
     &          (isubsp.eq.2.and.idiistyp.eq.4) ) then
              idiff = 0
              ludum = 0
              if (ntest.ge.100)
     &             write(6,*) 'initial sbspman for T-T(last)'
              call optc_sbspman(lust_sbsp,lustp,facs,
     &             ludum,nstdim,maxsbsp,
     &             isttask,nsttask,idiff,ndel_recent_st,
     &             vec1,vec2)
            end if
            if (isubsp.eq.2.and.idiistyp.ge.2) then
              idiff = 0
              if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for T+dT(pert)'
              if (idiistyp.ne.4) then
                call optc_sbspman(lutpst_sbsp,lustp,facs,
     &               ludum,ntpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_tpst,
     &               vec1,vec2)
              end if
              if (idiistyp.eq.2) then
                ludum = 0
                idiff = 0 
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for dT(pert)'
                call optc_sbspman(lupst_sbsp,lustp,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
              else if (idiistyp.eq.3.or.idiistyp.eq.4) then
                ludum = 0
                idiff = 0 
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for Omega/gradient'
                call optc_sbspman(lupst_sbsp,lugrvf,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
              end if
            end if
            nsttask = 0
          end if ! isubsp.eq.1.or.isubsp.eq.2.or.iprecnd.eq.2

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

      
*----------------------------------------------------------------------*
* check trust radius criterion
*----------------------------------------------------------------------*
        if (iprintl.ge.1) then
          write(6,*) ' current trust radius: ',trrad
        end if

        if (nrdvec.gt.0)
     &         call optc_prjout(nrdvec,lurdvec,lugrvf,
     &                          vec1,vec2,nwfpar,.false.)

        laccept = .true.
        xngrd_old = xngrd
        xngrd = sqrt(inprdd(vec1,vec2,lugrvf,lugrvf,1,lblk))
        de = energy - ener_old
        ener_old = energy
        if (ibstep.gt.0) then
          laccept = .false.
          if (ibstep.eq.3) then
            ibstep = 0
            laccept = .true.
          else if (ivar.eq.1.and.de.lt.0d0) then
            ibstep = 0
            laccept = .true.
          else if (ivar.eq.0.and.xngrd.lt.xngrd_old) then
            ibstep = 0
            laccept = .true.
          end if
        else if (ivar.eq.1.and.imacit.gt.1) then 
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
          else if (rat.gt.0d0) then
            trrad = trrad * trfac3
          else
            ! trigger maybe line-search backsteps
            trrad = trrad * trfac2
            xngrd = xngrd_old
            laccept = .false.
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
          crate = xngrd_old/xngrd
          if (imacit.gt.2) then
            ratio = crate/crate_prev
          else
            ratio = 1d0
          end if
          if (crate.gt.1d0.and.ratio.gt.trthr1l) then
            trrad = trrad * trfac1
          else if (crate.gt.1d0.and.
     &           (rat.gt.trthr2l.and.rat.le.trthr1l)) then
            trrad = trrad
          else if (crate.gt.1d0) then
            trrad = trrad * trfac2
          else
            trrad = trrad * trfac2
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
            if (.not.laccept)
     &           write(6,*) '# backstepping triggered'
          end if

        end if


        if (.not.laccept) then
          ibstep = ibstep + 1
          write(6,*) '# back-stepping by 1/2, step ',ibstep
          ! |X(new)> = |X(current)> - 0.5|dX(last)>
          call vecsmd(vec1,vec2,1d0,-0.5d0,
     &         luamp,lutrvec,luscr,1,lblk)
          ! save 0.5|dX(last)>:
          call vecsmd(vec1,vec2,1d0,-1.0d0,
     &         luamp,luscr,lutrvec,1,lblk)
          call copvcd(luscr,luamp,vec1,1,lblk)

          ! new energy
          itask = 1
          ! ... and gradient, if necessary
          if (ivar.eq.2) itask = itask + 2
          imicit = 0

          return

        end if

*-----------------------------------------------------------------------*
* subspace method needing delta gradient/delta vecfun?
* --->     call subspace manager
*-----------------------------------------------------------------------*
        if (iprecnd.eq.2.and.
     &      isbspjatyp.lt.10.and. 
     &      (imacit.gt.1.or.ilin.eq.1)) then
c        if (iprecnd.eq.2.and.imacit.gt.1) then
          igvtask(ngvtask+1) = 1
          igvtask(ngvtask+2) = 1
          facs(1) = 1d0
          ngvtask = ngvtask + 2
          
          ndel = max(0,ngvdim - 1 - nsbspjadim)
          if (ndel.gt.0) then
            igvtask(ngvtask+1) = 2
            igvtask(ngvtask+2) = ndel
            ngvtask = ngvtask + 2
          end if

c          if (ilin.eq.0) then
            if (ntest.ge.100)
     &           write(6,*) 'calling sbspman for Omg-Omg(last)'
            
            idiff = -1
            call optc_sbspman(lugv_sbsp,lugrvf,facs,
     &           lugrvfold,ngvdim,maxsbsp,
     &           igvtask,ngvtask,idiff,ndel_recent_gv,
     &           vec1,vec2)
c          else
c            if (ntest.ge.100)
c     &           write(6,*) 'calling sbspman for Jac*T'
c            
c            idiff = 0
c            call optc_sbspman(lugv_sbsp,lusig,facs,
c     &           ludum,ngvdim,maxsbsp,
c     &           igvtask,ngvtask,idiff,ndel_recent_gv,
c     &           vec1,vec2)
c          end if
          ngvtask = 0
        end if

*----------------------------------------------------------------------*
* 
*  perform Omega-DIIS, i.e. find the optimal solution in the space of
*  the previous parameters and the corresponding vectorfunction/gradient
*  the new step is then the reoptimized step within the iterative 
*  subspace plus the new step from the reoptimized residual determined
*  from the reoptimized Omega in the preconditioning/quasi-newton step
*
*  should be used with iprecnd==2
*  if used with iprecnd==1, the usual DIIS results (use idiis=2 instead)
*
*----------------------------------------------------------------------*

        if (isubsp.eq.2.and.idiistyp.eq.4) then
          if (imacit.ge.idiis_start) then
            nadd = 1
            ndiisdim = min(ndiisdim + nadd,mxsp_diis)
            ludiis_err = lupst_sbsp   ! Omega subspace
            ludiis_bas = lust_sbsp    ! step subspace  
            nsbspdim = npstdim

            call optc_diis(idiistyp,thr_diis,nsbspdim,
     &           ndiisdim,mxsp_diis,
     &           nadd,ndiisdim_last,alpha_last,
     &           lu_intstp,lu_corgrvf,ludiis_err,ludiis_bas,
     &           luamp,lugrvf,
     &           work(kb_mat),work(ksbscr),vec1,vec2,iprint)

            lustp = lu_corstp
            lugrd = lu_corgrvf
            if (ntest.ge.10) write(6,*) 'after diis: ndel_request = ',
     &           ndel_request_diis
          else
            lugrd = lugrvf
          end if
        else
          lugrd = lugrvf
        end if


*----------------------------------------------------------------------*
* get a new step direction
*  a) taking directly the gradient
*  b) taking a diagonal preconditioner (i.e. approx. diagonal Hessian)
*  c) use a subspace Hessian/Jacobian
*----------------------------------------------------------------------*
        if (iprecnd.eq.0) then
          lblk = -1
          call sclvcd(lugrd,lu_newstp,-1d0,vec1,1,lblk)
        else if (iprecnd.eq.1) then
*   use diagonal preconditioner/quasi-Newton step with diag. Hess.
          ldamp = iorder.ne.2 ! do not damp initial vector for 2nd order opt
          call optc_diagp(lugrd,ludia,trrad,
     &                    lu_newstp,xdamp,ldamp,vec1,vec2,nwfpar)
        else if (iprecnd.eq.2) then
          xdamp_old = xdamp
          lustp = lu_newstp
          if (imacit.ge.isbspja_start) lustp = luscr
          call optc_diagp(lugrd,ludia,trrad,
     &         lustp,xdamp,.true.,vec1,vec2,nwfpar)
          if (imacit.ge.isbspja_start) then
            nadd = 1
            nsbspjadim = min(nsbspjadim + nadd,mxsp_sbspja)
*   use (approximate) subspace Hessian/Jacobian
            if (isbspjatyp.lt.10) then
              ! rank n updates of diagonal Hessian
c TEST
c              xdamp = 0d0
c              if (imacit.eq.5) xdamp = 4d0
              call optc_sbspja_new(isbspjatyp,thr_sbspja,nstdim,ngvdim,
     &             nsbspjadim,mxsp_sbspja,
     &             nadd,nsbspjadim_last,
     &             lugrd,ludia,trrad,lustp,lu_newstp,
     &             xdamp,xdamp_old,
     &             lust_sbsp,lugv_sbsp,
     &             work(ku_mat),work(ksbscr),vec1,vec2,iprint)
            else
              ! rank 2 updates of previous Hessians
c              ndim_l = min(nsbspjadim - 1,mxsp_sbspja)
              ndim_l = min(nstdim,nhgdim+1)
              nrank = 2
              if (isbspjatyp.eq.15) nrank = 1
              call optc_updtja(isbspjatyp,nrank,thr_sbspja,
     &             nstdim,nhgdim,
     &             ndim_l,mxsp_sbspja,
     &             nadd,nsbspjadim_last,
     &             lugrvf,lugrvfold,
     &             ludia,trrad,lustp,lu_newstp,
     &             luhg_last,luhgam,
     &             xdamp,xdamp_old,
     &             lust_sbsp,luhg_sbsp,
     &             work(ksbscr),vec1,vec2,iprint)
            end if
          end if
        end if ! iprecnd

        if (ntest.ge.20) then
          xnrm = inprdd(vec1,vec2,lu_newstp,lu_newstp,1,lblk)
          write(6,*) 'New primitive step length: ',sqrt(xnrm)
          write(6,*) ' <s|s> was ',xnrm
        end if

        if (ntest.ge.1000) then
          write(6,*) 'New primitive step:'
          call vec_from_disc(vec1,nwfpar,1,lblk,lu_newstp)
          call wrt_cc_vec2(vec1,6,'GEN_CC')
c          call wrtvcd(vec1,lu_newstp,1,lblk)
        end if

        lustp = lu_newstp
        lustp2 = lugrvf
        if (isubsp.eq.1) then

* correct to approximate conjugate gradient direction
          if (imacit.gt.1) then
            call optc_conjgrad(icnjgrd,ilin,
     &           lu_corstp,
     &           lu_newstp,lust_sbsp,nstdim,
     &           lugrvf,lugrvfold,lusig,
     &           vec1,vec2,nwfpar,iprintl)
            lustp = lu_corstp
          end if
        else if (isubsp.eq.2) then
* use DIIS?
          if (imacit.ge.idiis_start) then
            if (idiistyp.lt.4) then
              nadd = 1
c     dbg
              print *,'last DIIS dim: ',ndiisdim
              print *,'last DIIS dim for B matrix: ', ndiisdim_last
c     dbg
              ndiisdim = min(ndiisdim + nadd,mxsp_diis)
c     dbg
              print *,'nadd : ',nadd
              print *,' ---> new DIIS dim = ',ndiisdim
c     dbg
              if (idiistyp.eq.1) then
                ludiis_err = lust_sbsp
                ludiis_bas = lust_sbsp              
                nsbspdim = nstdim
              else if (idiistyp.eq.2.or.idiistyp.eq.3) then
                ludiis_err = lupst_sbsp
                ludiis_bas = lutpst_sbsp
                nsbspdim = ntpstdim
              else
                write(6,*) 'unknown diistyp = ',idiistyp
                stop 'optcont'
              end if
              call optc_diis(idiistyp,thr_diis,nsbspdim,
     &             ndiisdim,mxsp_diis,
     &             nadd,ndiisdim_last,alpha_last,
     &             lu_newstp,lu_corstp,ludiis_err,ludiis_bas,
     &             luamp,lugrvf,
     &             work(kb_mat),work(ksbscr),vec1,vec2,iprint)

              lustp = lu_corstp
              if (ntest.ge.10) write(6,*) 'after diis: ndel_request = ',
     &             ndel_request_diis
            else if (idiistyp.eq.4) then            
c              stop 'weitermachen'
* add subspace-internal correction from Omega-DIIS
          ! add lu_newstp and lu_intstp givin lu_corstp
              if (ndiisdim.ge.2) then
                call vecsmd(vec1,vec2,1d0,1d0,lu_newstp,lu_intstp,
     &                      lu_corstp,1,lblk)
              else
                call copvcd(lu_newstp,lu_corstp,vec1,1,lblk)
              end if
              lustp = lu_corstp
            end if ! idiistyp
          end if ! imacit.ge.idiis_start

        else ! isubsp.eq.0
          lustp = lu_newstp
        end if

* norm of unscaled step:
        xnstp = sqrt(inprdd(vec1,vec2,lustp,lustp,1,lblk))

        if (ntest.ge.20) then
          xnrm = inprdd(vec1,vec2,lustp,lustp,1,lblk)
          write(6,*) 'New unscaled step length: ',sqrt(xnrm)
          write(6,*) ' <s|s> was ',xnrm
        end if

* line-search along new direction?
c        if (ilsrch.eq.0) then
c          alpha = 1d0
c        else
          ipass = 1
          call optc_linesearch(ilsrch,ilin,ivar,iprecnd,ipass,
     &         alpha,xdum,energy,de_pred,trrad,xnstp,
     &         vec1, vec2,
c     &         lustp,lu_newstp,ludia,iprint)
     &         lustp,lustp2,lusig,ludia,iprint)

          if (abs(de_pred).lt.1d-13.and.ilsrch.eq.2) then
            write(6,*) ' expected energy difference to small ',de_pred
            write(6,*) ' reverting to one-point model'
            ilsrch = 1
          end if

c        end if ! ilsrch

        alpha_last = alpha
c        if (alpha.ne.1d0.and.isubsp.eq.2) then
c          if (iprintl.ge.5) write(6,*) 'reset DIIS!'
c          ndel_request_diis=ndiisdim-nadd
c        end if

        xnstp = xnstp * alpha
        if (ntest.ge.20) then
          write(6,*) 'New unscaled step length: ',xnstp
          write(6,*) ' <s|s> was ',xnstp*xnstp
        end if


        if (ilsrch.eq.2) then
* if we need a second energy point:
* increase microiteration counter and request one more energy calculation
          itask  = 1
          imicit = 1
          imicit_tot = imicit_tot + 1
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
          end if
        end if

* obtain new paramter set |X> = |Xold> + alpha |d>
*  macro iteration --> see end-of-macro-iteration section
        alpha_eff = alpha
*  micro iteration:
        if (imicit.eq.1.and.ilsrch.eq.2) then
          call vecsmd(vec1,vec2,1d0,alpha,luamp,lustp,luscr,1,lblk)
          call copvcd(luscr,luamp,vec1,1,lblk)
          call copvcd(lustp,lutrvec,vec1,1,lblk)
c     &      call sclvcd(lustp,lutrvec,alpha,vec1,1,lblk)
        else if (imicit.eq.1) then
          ! normalize trial vector
          xnrm = sqrt(inprdd(vec1,vec2,lustp,lustp,1,lblk))
          ! and let it point into the direction of the gradient again
          call sclvcd(lustp,lutrvec,-1d0/xnrm,vec1,1,lblk)
          ! set threshold for microiterations
          ! hard: 1d-4
          ! weak: 1d1
          xmicthr = max(
     &        inprdd(vec1,vec2,lugrvf,lugrvf,1,lblk),!**(3d0/2d0),
     &                  thrgrd)
          xmicthr = min(xmicthr,0.1d0)
          write(6,'(" >>",x,a,e12.6)')
     &         'micro-iterations started -- threshold = ',xmicthr
        end if

* save old gradient for next iteration
        if (isubsp.gt.0.or.iprecnd.ge.2) 
     &       call copvcd(lugrvf,lugrvfold,vec1,1,lblk)

      else
*======================================================================*
* micro-iteration
*======================================================================*
        if (ntest.ge.10) then
          write(6,*) '*** micro iteration part entered ***'
        end if

        if (iorder.eq.1.and.ilsrch.eq.2) then
          ipass = 2
          xnstp = sqrt(inprdd(vec1,vec2,lutrvec,lutrvec,1,lblk))
          call optc_linesearch(ilsrch,ilin,ivar,iprecnd,ipass,
     &         alpha,alpha_eff,energy,de_pred,trrad,xnstp,
     &         vec1, vec2,
     &         lutrvec,lu_ucrstp,lusig,ludia,iprint)

          ! rescale step an add it to old parameter set          
          lustp = lu_corstp
          call sclvcd(lutrvec,lustp,alpha,vec1,1,lblk)
          ! well, this was done just for getting the step length
          xnstp = sqrt(inprdd(vec1,vec2,lustp,lustp,1,lblk))
          ! the really needed step length is alpha_eff (applied
          ! in the final section)
c          call sclvcd(lutrvec,lustp,dalp,vec1,1,lblk)
          lustp = lutrvec


          imicit = 0
          itask  = 1 + 2
        else if (iorder.eq.2) then

          imicdim = imicdim + 1
          maxit = 2*micifac
          fac = 1d0
          isttask(1) = 1
          isttask(2) = 1
          nsttask = 2
          ! push trial vector onto subspace
          call optc_sbspman(luvec_sbsp,lutrvec,fac,ludum,imicdim-1,
     &                                                         maxit,
     &           isttask,nsttask,0,ndum,vec1,vec2)
          ! push mv-product onto subspace
          call optc_sbspman(lumv_sbsp,lusig,fac,ludum,imicdim-1,maxit,
     &           isttask,nsttask,0,ndum,vec1,vec2)
          nsttask = 0

          ! update current reduced Hessian
          isymm = ivar ! enforce symmetric Hessian
     &                 ! (unless it is a Jacobian as for ivar==0)
          call optc_redh(isymm,work(khred),imicdim,imicdim-1,
     &                   lutrvec,lusig,luvec_sbsp,lumv_sbsp,
     &                   vec1,vec2)

          ! set up current reduced gradient
          call optc_sbspja_prjpstp(work(kgred),imicdim,0,
     &                             lugrvf,luvec_sbsp,
     &                             vec1,vec2,.false.)
          
          ! loop over the next routine call to enable mode-switches
          do int = 1, 2
            ! call TR-Newton to 
            !  - get new reduced expansion coeff.s and damping
            call optc_trnewton(i2nd_mode,iret,isymm,
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
              if (ntest.ge.100) write(6,*) '>> switched to Newton step'
              xlamb = 0d0
              i2nd_mode = 1
            else if (iret.ne.0.and.i2nd_mode.eq.1) then
              ! switch to newton-eigenvector procedure
              if (ntest.ge.100) write(6,*) '>> switched to Newton EV'
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
     &         luvec_sbsp,lumv_sbsp,lugrvf,
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
            call mvcsmd(luvec_sbsp,work(kcred),lu_newstp,luscr,
     &           vec1,vec2,imicdim,1,lblk)
            lustp = lu_newstp
            xnstp = sqrt(inprdd(vec1,vec1,lu_newstp,lu_newstp,1,lblk))
            alpha_eff = 1d0
            alpha = 1d0

            imicit = 0
            imicdim = 0
            itask  = 1 + 2
          else
            ! new raw step by preconditioning
            call dmtvcd2(vec1,vec2,ludia,lures,lutrvec,
     &           -1d0,-xlamb,1,1,lblk)

            ! orthogonalize against previous space and normalize;
            call optc_orthvec(work(kscr1),work(kscr2),
     &           imicdim,1,lin_dep,
     &           luvec_sbsp,lutrvec,
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
     &           luvec_sbsp,lutrvec,
     &           vec1,vec2)

              if (lin_dep) then
                write(6,*) ' unresolvable linear dependency problem!'
                stop 'optcont (2nd order)'
              end if

            end if

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
c        if (ntest.ge.10) then
c          write(6,*) 'end-of-macro-iteration part entered'
c        end if
*======================================================================*
* end of macro-iteration (indicated by imicit.eq.0):
*======================================================================*

*----------------------------------------------------------------------*
*  check convergence and max. iterations:
*----------------------------------------------------------------------*
        ! step criterion too expensive for 2nd order methods
        lstconv = (iorder.eq.2).or.xnstp.lt.thrstp
        lgrconv = xngrd.lt.thrgrd
        ! energy criterion makes no sense for 2nd order methods
        ldeconv = (iorder.eq.2).or.abs(de).lt.thr_de
        lconv = lstconv .and. lgrconv .and. ldeconv
        lexit = .false.

        if (iprintl.ge.1) then
          if (iorder.eq.1)
     &         write(6,*) 'after iteration ',imacit
          if (iorder.eq.2)
     &         write(6,*) 'in macro-iteration ',imacit
          if (iorder.ne.2)
     &         write(6,'(x,2(a,e10.3),a,l)')
     &                    '   norm of new step:  ', xnstp,
     &                           '   threshold:  ', thrstp,
     &                           '   converged:  ', lstconv
          if (ivar.eq.0)
     &         write(6,'(x,2(a,e10.3),a,l)')
     &                    '   norm of residual:  ', xngrd,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', lgrconv
          if (ivar.eq.1)
     &         write(6,'(x,2(a,e10.3),a,l)')
     &                    '   norm of gradient:  ', xngrd,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', lgrconv
          if (ilin.eq.0.and.iorder.ne.2)
     &        write (6,'(x,2(a,e10.3),a,l)')
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
          if (iprecnd.eq.2.or.isubsp.eq.1.or.
     &        (isubsp.eq.2.and.idiistyp.eq.1).or.
     &        (isubsp.eq.2.and.idiistyp.eq.4))
     &         call relunit(lust_sbsp,'delete')
          if (isubsp.eq.2.and.idiistyp.ge.2)
     &         call relunit(lupst_sbsp,'delete')
          if (isubsp.eq.2.and.idiistyp.eq.2.or.idiistyp.eq.3)
     &         call relunit(lutpst_sbsp,'delete')
          if (iprecnd.eq.2)
     &         call relunit(lugv_sbsp,'delete')
          if (iprecnd.eq.2.and.isbspjatyp.gt.10)
     &         call relunit(luhg_sbsp,'delete')
          if (iprecnd.eq.2.and.isbspjatyp.gt.10)
     &         call relunit(luhgam,'delete')
          if (iprecnd.eq.2.and.isbspjatyp.gt.10)
     &         call relunit(luhg_last,'delete')

          if (isubsp.eq.2.and.idiistyp.eq.4) then
            call relunit(lu_intstp,'delete')
            call relunit(lu_corgrvf,'delete')
          end if
          
          if (iorder.eq.2) then
            call relunit(luvec_sbsp,'delete')
            call relunit(lumv_sbsp,'delete')
            call relunit(lures,'delete')
          end if

          itask = 8 ! stop it

        else ! do some stuff for the next macro-iteration
* subspace method needing the steps? call subspace manager
          if (isubsp.eq.1.or.isubsp.eq.2.or.iprecnd.eq.2) then
            isttask(nsttask+1) = 1
            isttask(nsttask+2) = 1
            nsttask = nsttask + 2
            facs(1) = alpha
            ndel = 0
! manage here requests of diis and others to delete vectors
ccc currently only diis:
            if (isubsp.eq.2.and.idiistyp.eq.1) then
              ndel = nstdim + 1 - ndiisdim
            else if (isubsp.eq.2.and.idiistyp.gt.1) then
              ndel = ntpstdim + 1 - ndiisdim
            end if
ccc and subspace jac.
            if (iprecnd.eq.2) then
              ndel2 = max(0,nstdim - 1 - nsbspjadim)
              ndel = min(ndel,ndel2)
            end if

            if (ndel.gt.0) then
              isttask(nsttask+1) = 2
              isttask(nsttask+2) = ndel
              nsttask = nsttask + 2
            end if

c            if (lustp.ne.lu_newstp.and.lustp.ne.lu_corstp) then
c              write(6,*) 'unexpected unit before call to sbspman:',
c     &             lustp, ' expected: ',lu_newstp, lu_corstp
c              stop 'unexpected unit number'
c            end if

            if (iprecnd.eq.2.or.isubsp.eq.1.or.
     &          (isubsp.eq.2.and.idiistyp.eq.1).or.
     &          (isubsp.eq.2.and.idiistyp.eq.4) ) then
              idiff = 0
              ludum = 0
              if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for T-T(last)'
              call optc_sbspman(lust_sbsp,lustp,facs,
     &             ludum,nstdim,maxsbsp,
     &             isttask,nsttask,idiff,ndel_recent_st,
     &             vec1,vec2)
            end if
            if (isubsp.eq.2.and.idiistyp.ge.2) then
c              npstdim = ntpstdim
              idiff = 1     
              if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for T+dT(pert)'
              if (idiistyp.ne.4) then
                call optc_sbspman(lutpst_sbsp,luamp,facs,
     &               lu_newstp,ntpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_tpst,
     &               vec1,vec2)
              end if
              if (idiistyp.eq.2) then
                ludum = 0
                idiff = 0 
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for dT(pert)'
                call optc_sbspman(lupst_sbsp,lu_newstp,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
              else if (idiistyp.eq.3.or.idiistyp.eq.4) then
                ludum = 0
                idiff = 0 
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for Omega/gradient'
                call optc_sbspman(lupst_sbsp,lugrvf,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
              end if
              if (idiistyp.ne.4.and.(npstdim.ne.ntpstdim.or.
     &            ndel_recent_pst.ne.ndel_recent_tpst)) then
                write(6,*) 'dimension problem with p-step subspace:'
                write(6,*) 'dim: ',npstdim,ntpstdim
                write(6,*) 'del: ',ndel_recent_pst,ndel_recent_tpst
                stop 'optcont: dimension problem'
              end if
            end if
            nsttask = 0


          end if ! isubsp.eq.1.or.isubsp.eq.2.or.iprecnd.eq2

          if (iprecnd.eq.2.and.isbspjatyp.gt.10.and.
     &         imacit.gt.isbspja_start+1) then
* recursive update methods will need the previous H|delta g> product:
            ihgtask(nhgtask+1) = 1
            ihgtask(nhgtask+2) = 1
            nhgtask = nhgtask + 2
            ndel = 0
            facs(1) = 1d0
! manage here requests to delete vectors
            ndel = max(0,nhgdim - 1 - (nsbspjadim-1))

            if (ndel.gt.0) then
              ihgtask(nhgtask+1) = 2
              ihgtask(nhgtask+2) = ndel
              nhgtask = nhgtask + 2
            end if

            idiff = 0
            ludum = 0
            if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for H|g_i-g_i-1>'
            call optc_sbspman(luhg_sbsp,luhgam,facs,
     &           ludum,nhgdim,maxsbsp-1,
     &           ihgtask,nhgtask,idiff,ndel_recent_gh,
     &           vec1,vec2)
            nhgtask = 0

          end if ! iprecond.eq.2.and.isbspjatyp.gt.10

          ! build new amplitudes for next energy and gradient evalutation
          ! the step is expected on lustp the old amplitudes on luamp

          if (ntest.ge.1000) then
            write(6,*) 'New optimal step:'
            call vec_from_disc(vec1,nwfpar,1,lblk,lustp)
            call wrt_cc_vec2(vec1,6,'GEN_CC')
cc            call wrtvcd(vec1,lustp,1,lblk)
          end if

* obtain new paramter set |X> = |Xold> + alpha |d>
          call vecsmd(vec1,vec2,1d0,alpha_eff,
     &         luamp,lustp,luscr,1,lblk)
          call copvcd(luscr,luamp,vec1,1,lblk)
          if (lustp.ne.lutrvec)
     &         call copvcd(lustp,lutrvec,vec1,1,lblk)

          if (ntest.ge.1000) then
            write(6,*) 'New amplitudes:'
            call vec_from_disc(vec1,nwfpar,1,lblk,luamp)
            call wrt_cc_vec2(vec1,6,'GEN_CC')
c            call wrtvcd(vec1,luamp,1,lblk)
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
