*----------------------------------------------------------------------*
      subroutine optcont_lambda
     &                  (imacit,imicit,imicit_tot,iprint,
     &                   itask,iconv,i_oocc,itransf,
     &                   luamp,lulamp,lukappa,lutrvec,
     &                   energy,
     &                   vec1,vec2,nwfpar,norbr,
     &                   lugrvf,lulgrvf,luogrd,lusig,ludia,luodia)
*----------------------------------------------------------------------*
*
* General optimization control routine for non-linear optimization.
* "Special edition" for Lambda functional based CC
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
     &     ntest = 100, mxptsk = 10

* interface:
      integer, intent(out) ::
     &     itask, itransf
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
     &     lu_intstp,lu_corgrvf,luscr,luhg,
     &     lulgrvfold, lulsigold, lu_newlstp, lu_corlstp,
     &     lulst_sbsp, lulpst_sbsp, lultpst_sbsp, lulgv_sbsp,
     &     luscr_l,
     &     luogrvfold, luosigold, lu_newostp, lu_corostp,
     &     luost_sbsp, luopst_sbsp, luotpst_sbsp, luogv_sbsp,
     &     luscr_o,lukappa_sum,
     &     kumat,kbmat,ksbscr,
     &     kumat_l,kbmat_l,
     &     kumat_o,kbmat_o,
     &     ngvtask, nsttask, nhgtask,
     &     igvtask(mxptsk), isttask(mxptsk), ihgtask(mxptsk),
     &     nstdim,npstdim,ntpstdim,ngvdim,nhgdim,maxsbsp,
     &     ndiisdim, ndiisdim_last,nsbspjadim, nsbspjadim_last,
     &     nlstdim,nlpstdim,nltpstdim,nlgvdim,
     &     nldiisdim, nldiisdim_last,nlsbspjadim, nlsbspjadim_last,
     &     nostdim,nopstdim,notpstdim,nogvdim,
     &     nodiisdim, nodiisdim_last,nosbspjadim, nosbspjadim_last
      real(8), save ::
     &     ener_old, trrad, xngrd, xnstp, de, de_pred, alpha_last,
     &     xdamp, xdamp_old,
     &     xdamp_l, xdamp_l_old,
     &     xdamp_o, xdamp_o_old
      real(8) ::
     &     facs(mxptsk)

      logical ::
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
     &         'Optimization control (Lambda)',
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
        lulgrvfold = iopen_nus('OLDGRD_L')
        lulsigold  = iopen_nus('OLDSIG_L')
        lu_newlstp  = iopen_nus('NEWSTP_L')
        lu_corlstp  = iopen_nus('CORSTP_L')
        luscr     = iopen_nus('OPTSCR')
        luscr_l     = iopen_nus('OPTSCR_L')
        if (i_oocc.eq.1) then
          luogrvfold = iopen_nus('OLDGRD_O')
          luosigold  = iopen_nus('OLDSIG_O')
          lu_newostp  = iopen_nus('NEWSTP_O')
          lu_corostp  = iopen_nus('CORSTP_O')
          luscr_o     = iopen_nus('OPTSCR_O')
          lukappa_sum = iopen_nus('OPTKAPSM')
          vec1(1:norbr) = 0d0
          call vec_to_disc(vec1,norbr,1,lblk,lukappa_sum)
        end if
        if (iprecnd.eq.2.or.
     &      isubsp.eq.1 .or.
     &      (isubsp.eq.2.and.idiistyp.eq.1).or.
     &      (isubsp.eq.2.and.idiistyp.eq.4)) then
          lust_sbsp = iopen_nus('STPSBSP')
          lulst_sbsp = iopen_nus('STPSBSP_L')
          if (i_oocc.eq.1) then
            luost_sbsp = iopen_nus('STPSBSP_O')
          end if
        end if
        if (isubsp.eq.2.and.idiistyp.ge.2) then
          lupst_sbsp = iopen_nus('PSTPSBSP')
          lulpst_sbsp = iopen_nus('PSTPSBSP_L')
          if (i_oocc.eq.1) then
            luopst_sbsp = iopen_nus('PSTPSBSP_O')
          end if
        end if
        if (isubsp.eq.2.and.idiistyp.eq.2.or.idiistyp.eq.3) then
          lutpst_sbsp = iopen_nus('TPSTPSBSP')
          lultpst_sbsp = iopen_nus('TPSTPSBSP_L')
          if (i_oocc.eq.1) then
            luotpst_sbsp = iopen_nus('TPSTPSBSP_O')
          end if
        end if
        if (iprecnd.eq.2) then
          lugv_sbsp = iopen_nus('GRDSBSP')
          lulgv_sbsp = iopen_nus('GRDSBSP_L')
          if (i_oocc.eq.1) then
            luogv_sbsp = iopen_nus('GRDSBSP_O')
          end if
        end if
        if (iprecnd.eq.2.and.isbspjatyp.gt.10) then
          stop 'this route is not active in Lamdba edition'
        end if

        if (isubsp.eq.2.and.idiistyp.eq.4) then
          stop 'this route is not active in Lambda edition'
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

        ndiisdim = 0
        ndiisdim_last = 0
        nldiisdim = 0
        nldiisdim_last = 0

        nsbspjadim = 0
        nsbspjadim_last = 0
        nlsbspjadim = 0
        nlsbspjadim_last = 0

        xdamp = 0d0
        xdamp_old = 0d0
        xdamp_l = 0d0
        xdamp_l_old = 0d0

* work space for subspace jacobian?
        maxsbsp = 0
        if (iprecnd.eq.2) then
          maxsbsp = mxsp_sbspja
          lscr_sbspja = maxsbsp**2
          lenu= maxsbsp**2
          call memman(kumat,lenu   ,'ADDL  ',2,'U SBJA')
          call memman(kumat_l,lenu   ,'ADDL  ',2,'ULSBJA')
          call memman(kumat_o,lenu   ,'ADDL  ',2,'UOSBJA')
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
          call memman(kbmat_l,lenb,'ADDL  ',2,'BMLDIS')
          call memman(kbmat_o,lenb,'ADDL  ',2,'BMODIS')
        else
          lscr_diis = 0
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

* subspace method needing delta gradient/delta vecfun?
* --->     call subspace manager
        if (iprecnd.eq.2.and.imacit.gt.1) then
          igvtask(ngvtask+1) = 1
          igvtask(ngvtask+2) = 1
          facs(1) = 1d0
          ngvtask = ngvtask + 2
          
          ndel1 = max(0,ngvdim - 1 - nsbspjadim)
          ndel2 = max(0,nlgvdim - 1 - nlsbspjadim)
          ndel3 = max(0,nogvdim - 1 - nosbspjadim)
          ndel = min(ndel1,ndel2,ndel3)
          if (ndel.gt.0) then
            igvtask(ngvtask+1) = 2
            igvtask(ngvtask+2) = ndel
            ngvtask = ngvtask + 2
          end if

          if (ntest.ge.100)
     &         write(6,*) 'calling sbspman for diff. t-gradient'
            
          idiff = -1
          call optc_sbspman(lugv_sbsp,lugrvf,facs,
     &         lugrvfold,ngvdim,maxsbsp,
     &         igvtask,ngvtask,idiff,ndel_recent_gv,
     &         vec1,vec2)

          if (ntest.ge.100)
     &         write(6,*) 'calling sbspman for diff. tbar-gradient'

          call optc_sbspman(lulgv_sbsp,lulgrvf,facs,
     &         lulgrvfold,nlgvdim,maxsbsp,
     &         igvtask,ngvtask,idiff,ndel_recent_gv,
     &         vec1,vec2)

          ngvtask = 0
        end if
      
*----------------------------------------------------------------------*
* check trust radius criterion
*----------------------------------------------------------------------*
        if (iprintl.ge.1) then
          write(6,*) ' current trust radius: ',trrad
        end if

        xngrd   = sqrt(inprdd(vec1,vec2,lugrvf,lugrvf,1,lblk))
        xngrd_l = sqrt(inprdd(vec1,vec2,lulgrvf,lulgrvf,1,lblk))
        if (i_oocc.eq.1) then
          xngrd_o = sqrt(inprdd(vec1,vec2,luogrd,luogrd,1,lblk))
        end if
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
*   I you didn't know by now: the t gradient determines the new tbar
*                             the tbar gradient the new t
          call optc_diagp(lugrvf,ludia,trrad,
     &                    lu_newlstp,xdamp_l,.true.,vec1,vec2,nwfpar)
          call optc_diagp(lulgrvf,ludia,trrad,
     &                    lu_newstp,xdamp,.true.,vec1,vec2,nwfpar)
          if (i_oocc.eq.1) then
            call optc_diagp(luogrd,luodia,trrad,
     &           lu_newostp,xdamp_o,.true.,vec1,vec2,norbr)
          end if
        else if (iprecnd.eq.2) then
          xdamp_old = xdamp
          xdamp_l_old = xdamp_l
          lulstp = lu_newlstp
          lustp = lu_newstp
          if (imacit.ge.isbspja_start) lulstp = luscr_l
          if (imacit.ge.isbspja_start) lustp = luscr
          call optc_diagp(lugrvf,ludia,trrad,
     &         lulstp,xdamp_l,.false.,vec1,vec2,nwfpar)
          call optc_diagp(lulgrvf,ludia,trrad,
     &         lustp,xdamp,.false.,vec1,vec2,nwfpar)
          if (i_oocc.eq.1) then
            xdamp_o_old = xdamp_o
            luostp = lu_newostp
            if (imacit.ge.isbspja_start) luostp = luscr_o
            call optc_diagp(luogrd,luodia,trrad,
     &           luostp,xdamp_o,.false.,vec1,vec2,norbr)
          end if
          if (imacit.ge.isbspja_start) then
c            nadd = 1
c            nlsbspjadim = min(nlsbspjadim + nadd,mxsp_sbspja)
c*   use (approximate) subspace Hessian/Jacobian
c            call optc_sbspja(isbspjatyp,thr_sbspja,nlstdim,ngvdim,
c     &           nlsbspjadim,mxsp_sbspja,
c     &           nadd,nlsbspjadim_last,
c     &           lugrvf,ludia,trrad,lulstp,lu_newlstp,
c     &           xdamp_l,xdamp_l_old,
c     &           lulst_sbsp,lugv_sbsp,
c     &           work(kumat_l),work(ksbscr),vec1,vec2,iprint)
c            nadd = 1
c            nsbspjadim = min(nsbspjadim + nadd,mxsp_sbspja)
c            call optc_sbspja(isbspjatyp,thr_sbspja,nstdim,nlgvdim,
c     &           nsbspjadim,mxsp_sbspja,
c     &           nadd,nsbspjadim_last,
c     &           lulgrvf,ludia,trrad,lustp,lu_newstp,
c     &           xdamp,xdamp_old,
c     &           lust_sbsp,lulgv_sbsp,
c     &           work(kumat),work(ksbscr),vec1,vec2,iprint)
            if (i_oocc.eq.1) then
              nadd = 1
              nosbspjadim = min(nosbspjadim + nadd,mxsp_sbspja)
              call optc_sbspja(isbspjatyp,thr_sbspja,nostdim,nogvdim,
     &           nosbspjadim,mxsp_sbspja,
     &           nadd,nosbspjadim_last,
     &           luogrd,luodia,trrad,luostp,lu_newostp,
     &           xdamp_o,xdamp_o_old,
     &           luost_sbsp,luogv_sbsp,
     &           work(kumat_o),work(ksbscr),vec1,vec2,iprint)
            end if

          end if
        end if ! iprecnd

c        if (ntest.ge.20) then
c          xnrm = inprdd(vec1,vec2,lu_newstp,lu_newstp,1,lblk)
c          write(6,*) 'New primitive step length: ',sqrt(xnrm)
c          write(6,*) ' <s|s> was ',xnrm
c        end if

        lustp = lu_newstp
        lustp2 = lulgrvf
        lulstp = lu_newlstp
        lulstp2 = lugrvf
        if (i_oocc.eq.1) then
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
              nadd = 1
              nldiisdim = min(nldiisdim + nadd,mxsp_diis)
              ludiis_err = lulpst_sbsp
              ludiis_bas = lultpst_sbsp
              nsbspdim = nltpstdim
              call optc_diis(idiistyp,thr_diis,nsbspdim,
     &             nldiisdim,mxsp_diis,
     &             nadd,nldiisdim_last,alpha_last,
     &             lu_newlstp,lu_corlstp,ludiis_err,ludiis_bas,
     &             lulamp,lugrvf,
     &             work(kbmat_l),work(ksbscr),vec1,vec2,iprint)
              nadd = 1
              ndiisdim = min(ndiisdim + nadd,mxsp_diis)
              ludiis_err = lupst_sbsp
              ludiis_bas = lutpst_sbsp
              nsbspdim = ntpstdim
              call optc_diis(idiistyp,thr_diis,nsbspdim,
     &             ndiisdim,mxsp_diis,
     &             nadd,ndiisdim_last,alpha_last,
     &             lu_newstp,lu_corstp,ludiis_err,ludiis_bas,
     &             luamp,lulgrvf,
     &             work(kbmat),work(ksbscr),vec1,vec2,iprint)
              if (i_oocc.eq.1) then
                nadd = 1
                nodiisdim = min(nodiisdim + nadd,mxsp_diis)
                ludiis_err = luopst_sbsp
                ludiis_bas = luotpst_sbsp
                nsbspdim = notpstdim
                call optc_diis(idiistyp,thr_diis,nsbspdim,
     &             nodiisdim,mxsp_diis,
     &             nadd,nodiisdim_last,alpha_last,
     &             lu_newostp,lu_corostp,ludiis_err,ludiis_bas,
     &             lukappa_sum,luogrd,
     &             work(kbmat_o),work(ksbscr),vec1,vec2,iprint)

              end if

              lustp = lu_corstp
              lulstp = lu_corlstp
              if (i_oocc.eq.1) luostp = lu_corostp
          end if ! imacit.ge.idiis_start

        else ! isubsp.eq.0
          lustp = lu_newstp
          lulstp = lu_newlstp
          if (i_oocc.eq.1) luostp = lu_newostp
        end if

* norm of unscaled step:
        xnstp = sqrt(inprdd(vec1,vec2,lustp,lustp,1,lblk))
        xnstp_l = sqrt(inprdd(vec1,vec2,lulstp,lulstp,1,lblk))
        if (i_oocc.eq.1) 
     &      xnstp_o = sqrt(inprdd(vec1,vec2,luostp,luostp,1,lblk))

        if (ntest.ge.20) then
          write(6,*) 'New unscaled step length: ',xnstp,xnstp_l
          write(6,*) ' <s|s> was ',xnstp*xnstp,xnstp_l*xnstp_l
        end if

* line-search along new direction?
c        if (ilsrch.eq.0) then
c          alpha = 1d0
c        else
c          ipass = 1
c          call optc_linesearch(ilsrch,ilin,ivar,iprecnd,ipass,
c     &         alpha,xdum,energy,de_pred,trrad,xnstp,
c     &         vec1, vec2,
cc     &         lustp,lu_newstp,ludia,iprint)
c     &         lustp,lustp2,lusig,ludia,iprint)
c        end if ! ilsrch

        cont1   = inprdd(vec1,vec2,lustp,lugrvf,1,lblk)
        cont2   = inprdd(vec1,vec2,lulstp,lulgrvf,1,lblk)
        de_pred = 0.5d0*(cont1+cont2)
        if (i_oocc.eq.1) then
          cont3 = inprdd(vec1,vec2,luostp,luogrd,1,lblk)
          de_pred = de_pred + cont3
        end if


        alpha = 1d0
        alpha_last = alpha
c        if (alpha.ne.1d0.and.isubsp.eq.2) then
c          if (iprintl.ge.5) write(6,*) 'reset DIIS!'
c          ndel_request_diis=ndiisdim-nadd
c        end if


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
            imicit_tot = imicit_tot + 1
          end if
        end if


* obtain new paramter set |X> = |Xold> + alpha |d>
*  macro iteration --> see end-of-macro-iteration section
        alpha_eff = alpha
*  micro iteration:
        if (imicit.eq.1.and.ilsrch.eq.2) then
          stop 'unexpected event'
        else if (imicit.eq.1) then
          stop 'unexpected event 2'
        end if

* save old gradient for next iteration
        if (isubsp.gt.0.or.iprecnd.ge.2) then 
          call copvcd(lugrvf,lugrvfold,vec1,1,lblk)
          call copvcd(lulgrvf,lulgrvfold,vec1,1,lblk)
          if (i_oocc.eq.1) call copvcd(luogrd,luogrvfold,vec1,1,lblk)
        end if

      else
*======================================================================*
* micro-iteration
*======================================================================*
        if (ntest.ge.10) then
          write(6,*) 'micro iteration part entered'
        end if

        stop 'under construction'

      end if ! macro/micro iteration switch

      if (imicit.eq.0) then
        if (ntest.ge.10) then
          write(6,*) 'end-of-macro-iteration part entered'
        end if
*======================================================================*
* end of macro-iteration (indicated by imicit.eq.0):
*======================================================================*

*----------------------------------------------------------------------*
*  check convergence and max. iterations:
*----------------------------------------------------------------------*
        lstconv = xnstp.lt.thrstp
        lgrconv = xngrd.lt.thrgrd
        llstconv = xnstp_l.lt.thrstp
        llgrconv = xngrd_l.lt.thrgrd
        ldeconv = abs(de).lt.thr_de
        lconv = lstconv .and. lgrconv .and. ldeconv .and.
     &          llstconv .and. llgrconv
        lexit = .false.

        if (i_oocc.eq.1) then
          lostconv = xnstp_o.lt.thrstp
          logrconv = xngrd_o.lt.thrgrd
          lconv = lconv .and. lostconv .and. logrconv
          if (lostconv.and.logrconv) then
            itransf = 0
          else
            itransf = 1
          end if
        end if

        if (iprintl.ge.1) then
          if (iorder.eq.1)
     &         write(6,*) 'after iteration ',imacit
          if (iorder.eq.2)
     &         write(6,*) 'after macro-iteration ',imacit
          write (6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of new t step:  ', xnstp,
     &                           '   threshold:  ', thrstp,
     &                           '   converged:  ', lstconv
          write(6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of t gradient:  ', xngrd,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', lgrconv
          write (6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of new l step:  ', xnstp_l,
     &                           '   threshold:  ', thrstp,
     &                           '   converged:  ', llstconv
          write(6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of l gradient:  ', xngrd_l,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', llgrconv
          if (i_oocc.eq.1) then
            write (6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of new k step:  ', xnstp_o,
     &                           '   threshold:  ', thrstp,
     &                           '   converged:  ', lostconv
            write(6,'(x,2(a,e10.3),a,l)')
     &                    ' norm of k gradient:  ', xngrd_o,
     &                           '   threshold:  ', thrgrd,
     &                           '   converged:  ', logrconv
          end if
          write (6,'(x,2(a,e10.3),a,l)')
     &                    '   change in energy:  ', de,
     &                           '   threshold:  ', thr_de,
     &                           '   converged:  ', ldeconv
          if (lconv) write(6,*) 'CONVERGED'
          if (lconv) iconv = 1
        end if

        if (.not.lconv) imacit = imacit + 1

        if (.not.lconv.and.
     &       (imacit.gt.maxmacit.or.
     &       imicit_tot.gt.maxmicit)) then
          write(6,*) 'NO CONVERGENCE OBTAINED'
          imacit = imacit - 1
          imicit = imicit - 1
          lexit = .true.
        end if

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
          call relunit(lulgrvfold,'delete')
          call relunit(lulsigold,'delete')
          call relunit(lu_newlstp,'delete')
          call relunit(lu_corlstp,'delete')
          call relunit(luscr_l,'delete')
          if (iprecnd.eq.2.or.isubsp.eq.1.or.
     &       (isubsp.eq.2.and.idiistyp.eq.1).or.
     &        (isubsp.eq.2.and.idiistyp.eq.4)) then
            call relunit(lust_sbsp,'delete')
            call relunit(lulst_sbsp,'delete')
          end if
          if (isubsp.eq.2.and.idiistyp.ge.2) then
            call relunit(lupst_sbsp,'delete')
            call relunit(lulpst_sbsp,'delete')
          end if
          if (isubsp.eq.2.and.idiistyp.eq.2.or.idiistyp.eq.3) then
            call relunit(lutpst_sbsp,'delete')
            call relunit(lultpst_sbsp,'delete')
          end if
          if (iprecnd.eq.2) then
            call relunit(lugv_sbsp,'delete')
            call relunit(lulgv_sbsp,'delete')
          end if
          
          if (i_oocc.eq.1) then
            call relunit(luogrvfold,'delete')
            call relunit(luosigold,'delete')
            call relunit(lu_newostp,'delete')
            call relunit(lu_corostp,'delete')
            call relunit(luscr_o,'delete')
            call relunit(lukappa_sum,'delete')
            if (iprecnd.eq.2.or.isubsp.eq.1.or.
     &           (isubsp.eq.2.and.idiistyp.eq.1).or.
     &           (isubsp.eq.2.and.idiistyp.eq.4)) then
              call relunit(luost_sbsp,'delete')
            end if
            if (isubsp.eq.2.and.idiistyp.ge.2) then
              call relunit(luopst_sbsp,'delete')
            end if
            if (isubsp.eq.2.and.idiistyp.eq.2.or.idiistyp.eq.3) then
              call relunit(luotpst_sbsp,'delete')
            end if
            if (iprecnd.eq.2) then
              call relunit(luogv_sbsp,'delete')
            end if
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
              ndel2 = nlstdim + 1 - nldiisdim
              ndel = min(ndel,ndel2)
              if (i_oocc.eq.1) then
                ndel3 = nostdim + 1 -nodiisdim
                ndel = min(ndel,ndel3)
              end if
            else if (isubsp.eq.2.and.idiistyp.gt.1) then
              ndel = ntpstdim + 1 - ndiisdim
              ndel2 = nltpstdim + 1 - nldiisdim
              ndel = min(ndel,ndel2)
              if (i_oocc.eq.1) then
                ndel3 = notpstdim + 1 - nodiisdim
                ndel = min(ndel,ndel3)
              end if
            end if
ccc and subspace jac.
            if (iprecnd.eq.2) then
              ndel2 = max(0,nstdim - 1 - nsbspjadim)
              ndel3 = max(0,nlstdim - 1 - nlsbspjadim)
              ndel = min(ndel,ndel2,ndel3)
              if (i_oocc.eq.1) then
                ndel4 = nostdim - 1 - nosbspjadim
                ndel = min(ndel,ndel4)
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
              if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for L-L(last)'
              call optc_sbspman(lulst_sbsp,lulstp,facs,
     &             ludum,nlstdim,maxsbsp,
     &             isttask,nsttask,idiff,ndel_recent_st,
     &             vec1,vec2)
              if (ntest.ge.100)
     &             write(6,*) 'calling sbspman for T-T(last)'
              call optc_sbspman(lust_sbsp,lustp,facs,
     &             ludum,nstdim,maxsbsp,
     &             isttask,nsttask,idiff,ndel_recent_st,
     &             vec1,vec2)
              if (i_oocc.eq.1) then
                if (ntest.ge.100)
     &               write(6,*) 'calling sbspman for kappa-kappa(last)'
                call optc_sbspman(luost_sbsp,luostp,facs,
     &               ludum,nostdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_st,
     &               vec1,vec2)
              end if
            end if
            if (isubsp.eq.2.and.idiistyp.ge.2) then
c              npstdim = ntpstdim
              if (idiistyp.ne.4) then
                idiff = 1     
                if (ntest.ge.100)
     &               write(6,*) 'calling sbspman for L+dL(pert)'
                call optc_sbspman(lultpst_sbsp,lulamp,facs,
     &               lu_newlstp,nltpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_tpst,
     &               vec1,vec2)
                if (ntest.ge.100)
     &               write(6,*) 'calling sbspman for T+dT(pert)'
                call optc_sbspman(lutpst_sbsp,luamp,facs,
     &               lu_newstp,ntpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_tpst,
     &               vec1,vec2)
                if (i_oocc.eq.1) then
                  if (ntest.ge.100)
     &              write(6,*) 'calling sbspman for kappa+dkappa(pert)'
                  call optc_sbspman(luotpst_sbsp,lukappa_sum,facs,
     &                 lu_newostp,notpstdim,maxsbsp,
     &                 isttask,nsttask,idiff,ndel_recent_tpst,
     &                 vec1,vec2)
                end if
              end if
              if (idiistyp.eq.2) then
                ludum = 0
                idiff = 0 
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for dL(pert)'
                call optc_sbspman(lulpst_sbsp,lu_newlstp,facs,
     &               ludum,nlpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for dT(pert)'
                call optc_sbspman(lupst_sbsp,lu_newstp,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
                if (i_oocc.eq.1) then
                  if (ntest.ge.100) write(6,*)
     &                 'calling sbspman for dkappa(pert)'
                  call optc_sbspman(luopst_sbsp,lu_newostp,facs,
     &                 ludum,nopstdim,maxsbsp,
     &                 isttask,nsttask,idiff,ndel_recent_pst,
     &                 vec1,vec2)
                end if
              else if (idiistyp.eq.3.or.idiistyp.eq.4) then
                ludum = 0
                idiff = 0 
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for T gradient = L residual'
                call optc_sbspman(lulpst_sbsp,lugrvf,facs,
     &               ludum,nlpstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
                if (ntest.ge.100) write(6,*)
     &               'calling sbspman for L gradient = T residual'
                call optc_sbspman(lupst_sbsp,lulgrvf,facs,
     &               ludum,npstdim,maxsbsp,
     &               isttask,nsttask,idiff,ndel_recent_pst,
     &               vec1,vec2)
                if (i_oocc.eq.1) then
                  if (ntest.ge.100) write(6,*)
     &                 'calling sbspman for kappa gradient'
                  call optc_sbspman(luopst_sbsp,luogrvf,facs,
     &                 ludum,nopstdim,maxsbsp,
     &                 isttask,nsttask,idiff,ndel_recent_pst,
     &                 vec1,vec2)
                end if
              end if
            end if
            nsttask = 0


          end if ! isubsp.eq.1.or.isubsp.eq.2.or.iprecnd.eq2

          ! build new amplitudes for next energy and gradient evalutation
          ! the step is expected on lustp the old amplitudes on luamp

* obtain new paramter set |X> = |Xold> + alpha |d>
c TEST TEST TEST -- only MP2 update of amplitudes
c          if (imacit.le.2) then
          call vecsmd(vec1,vec2,1d0,alpha_eff,
     &         lulamp,lulstp,luscr,1,lblk)
          call copvcd(luscr,lulamp,vec1,1,lblk)
          
          call vecsmd(vec1,vec2,1d0,alpha_eff,
     &         luamp,lustp,luscr,1,lblk)
          call copvcd(luscr,luamp,vec1,1,lblk)
c          else 
c            do i = 1,10
c              print *,' ==== TEST active ==== '
c            end do
c          end if

          ! for kappa we pass only the step
          if (i_oocc.eq.1) then
            call vecsmd(vec1,vec2,1d0,alpha_eff,
     &           lukappa_sum,luostp,luscr,1,lblk)
            call copvcd(luscr,lukappa_sum,vec1,1,lblk)
            call copvcd(luostp,lukappa,vec1,1,lblk)
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
