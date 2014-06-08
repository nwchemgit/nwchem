*----------------------------------------------------------------------*
      subroutine optc_sbspman(lusbsp,lunew,facs,lunew2,nvec,nmaxvec,
     &                        itask,ntask,iaddif,ndel_out,
     &                        vec1,vec2)
*----------------------------------------------------------------------*
*
* purpose: manage subspace of vectors on file lusbsp. at most nmaxvec
*          vectors will be kept on file. itask contains the tasks:
*
*          1: add n vectors (where n is contained in the succeeding
*                            field of itask)
*          2: delete n vectors (n: see above). the number of vectors
*             that need to be deleted anyway to keep the subspace
*             dimensions is subtracted from this
*          3: restart (all vectors will be deleted)
*
*     nvec is current size of subspace
*
*   iaddif = 1: add sums   lunew + lunew2
*   iaddif = 1: add diffs  lunew - lunew2
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00
      
* input/output
      integer, intent(in) ::
     &     lusbsp, lunew, lunew2,
     &     nmaxvec, ntask,
     &     itask(ntask), iaddif
      integer, intent(inout) ::
     &     nvec
      integer, intent(out) ::
     &     ndel_out
      real(8), intent(in) ::
     &     facs(*)
* scratch
      real(8), intent(inout) ::
     &     vec1(*), vec2(*)
* external function
      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd
* local
      integer ::
     &     i, ierr, ndel, nadd, nnew, lblk,
     &     luscr, irec
      real(8) ::
     &     xnrm

      lblk = -1

      if (ntest.ge.10) then
        write(6,*) 'optc_sbspman:'
        write(6,*) '=============='
        write(6,*) ' nvec, nmaxvec : ', nvec, nmaxvec
        write(6,*) ' ntask,itask :', ntask,'(',itask(1:ntask),')'
        write(6,*) ' lusbsp, lunew: ',lusbsp, lunew
        if (iabs(iaddif).eq.1) 
     &    write(6,*) ' lnew2 :',lunew2
        if (ntest.ge.100) then
          write(6,*) 'Initial contents on subspace-file:'
          write(6,*) '   record  norm(vector)'
          call rewino(lusbsp)
          do irec = 1, nvec
            xnrm = sqrt(inprdd(vec1,vec1,lusbsp,lusbsp,0,lblk))
            write(6,'(4x,i5,2x,e15.7)') irec, xnrm
          end do
        end if
      end if

      i = 0
      ierr = 0
      ndel = 0
      nadd = 0
      do while(i.lt.ntask)
        i = i+1
        if (itask(i).eq.1) then
          i = i+1
          if (i.gt.ntask) then
            ierr = 1
            exit
          end if
          nadd = nadd + itask(i)
        else if (itask(i).eq.2) then
          i = i+1
          if (i.gt.ntask) then
            ierr = 2
            exit
          end if
          ndel = min(ndel+itask(i),nvec)
        else if (itask(i).eq.3) then
          ndel = nvec
        else
          ierr = 4
          exit
        end if
      end do

      if (ierr.ne.0) then
        write(6,*) 'Error in optc_sbspman: ',ierr
        stop 'error in optc_sbspman'
      end if

      nnew = nvec + nadd - ndel
      if (nnew.gt.nmaxvec) then
        ndel = ndel + nnew-nmaxvec
      end if

      if (ntest.ge.10) then
        write(6,*) 'nadd, ndel: ',nadd, ndel
      end if

      if (ndel.gt.0.and.nvec-ndel.gt.0) then
        luscr = iopen_nus('sbspmanscr')
        call rewino(luscr)
        call skpvcd(lusbsp,ndel,vec1,1,lblk)
        do i = ndel+1,nvec
          call copvcd(lusbsp,luscr,vec1,0,lblk)
        end do
        call rewino(luscr)
        call rewino(lusbsp)
        do i = ndel+1,nvec
          call copvcd(luscr,lusbsp,vec1,0,lblk)
        end do
        call relunit(luscr,'delete')
      else if (ndel.eq.nvec) then
        call rewino(lusbsp)
      else if (ndel.eq.0) then
        call skpvcd(lusbsp,nvec,vec1,1,lblk)
      else
        write(6,*) 'Inconsistency in optc_sbspman'
        stop 'Inconsistency in optc_sbspman'
      end if

      if (nadd.gt.0) then
        call rewino(lunew)
        if (iabs(iaddif).eq.1) call rewino(lunew2)
        do i = 1, nadd
          if (iaddif.eq.1) then
            call vecsmd(vec1,vec2,facs(i),facs(i),
     &           lunew,lunew2,lusbsp,0,lblk)
          else if (iaddif.eq.-1) then
            call vecsmd(vec1,vec2,facs(i),-facs(i),
     &           lunew,lunew2,lusbsp,0,lblk)
          else
            call sclvcd(lunew,lusbsp,facs(i),vec1,0,lblk)
          end if
        end do
      end if

      ! return new subspace size ...
      nvec = nvec - ndel + nadd

      ! ... and how many vectors were deleted
      ndel_out = ndel

      if (ntest.ge.100) then
        write(6,*) 'Final contents on subspace-file:'
        write(6,*) '   record  norm(vector)'
        call rewino(lusbsp)
        do irec = 1, nvec
          xnrm = sqrt(inprdd(vec1,vec1,lusbsp,lusbsp,0,lblk))
          write(6,'(4x,i5,2x,e15.7)') irec, xnrm
        end do
      end if

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_diagp(lugrvf,ludia,xmaxstp,
     &                      luamp,xdamp,ldamp,
     &                      vec1,vec2,nwfpar)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     ldamp
      integer, intent(in) ::
     &     lugrvf, ludia, luamp,
     &     nwfpar
      real(8), intent(in) ::
     &     vec1(nwfpar), vec2(nwfpar),
     &     xmaxstp
      real(8), intent(out) ::
     &     xdamp
      integer ::
     &     inv, lblk

      logical ::
     &     lconv
      integer ::
     &     maxiter, iter
      real(8) ::
     &     xnrm, xgg, xval, xder, xthr, xmx, xdamp_min,
     &     xdamp_old, xval_old, diamin, xstep, xhigh, xlow,
     &     xvhigh, xvlow
      
      real(8), external ::
     &     inprdd, inprdd3, fdmnxd

      inv=1
      lblk = -1
      diamin = fdmnxd(ludia,-1,vec1,1,lblk)
      xdamp_min = max(-diamin,0d0)
      if (ntest.ge.100) write(6,*) 'smallest diagonal element: ',diamin

      call dmtvcd2(vec1,vec2,ludia,lugrvf,luamp,-1d0,0d0,1,inv,lblk)

      xnrm = sqrt(inprdd(vec1,vec1,luamp,luamp,1,lblk))

      if (ntest.ge.100) write(6,*) 'norm of predicted step: ',xnrm
      xdamp = 0d0

      if ((xnrm.gt.xmaxstp.or.diamin.lt.0d0).and..not.ldamp) then
        xdamp = 100d0 ! just used as a flag, another routine will
                      ! evaluate the actual damping
      end if

      xlow = xdamp_min
      xvlow = +100
      xhigh= 1d13
      xvhigh= -100

* step length OK?
      if ((xnrm.gt.xmaxstp.or.diamin.lt.0d0).and.ldamp) then
        ! else we have to search for the lowest positive value of
        ! xdamp such that that condition is fulfilled
        xgg = inprdd(vec1,vec1,lugrvf,lugrvf,1,lblk)
        ! a guess, assuming sqrt(<1/H^2> ~ 0.4d0)
        xdamp = max(sqrt(xgg)/xmaxstp - 0.4d0,0.1d0,-diamin+0.1d0)
c        xdamp = 0
        maxiter = 100
        xmx  = 1.d0
        xthr = 1d-12
        do iter = 1, maxiter
          ! get the value
          xval=inprdd3(vec1,vec2,lugrvf,lugrvf,ludia,xdamp,-2d0,1,lblk)
     &         - xmaxstp*xmaxstp

c?          if (iter.gt.1.and.sign(1d0,xval).ne.sign(1d0,xval_old)
c?     &         .and.abs(xval_old)/10d0.lt.abs(xval)) then
c?            xdamp = xdamp_old -
c?c     &           xval_old*(xdamp_old-xdamp)/(xval_old-xval)
c?     &            0.5d0*(xdamp_old-xdamp)
c?
c?           xval=inprdd3(vec1,vec2,lugrvf,lugrvf,ludia,xdamp,-2d0,1,lblk)
c?     &         - xmaxstp*xmaxstp
c?          end if

c          if (iter.gt.1.and.abs(xval).gt.abs(xval_old)) then
c            xdamp = xdamp_old + (xdamp_old - xdamp)
c           xval=inprdd3(vec1,vec2,lugrvf,lugrvf,ludia,xdamp,-2d0,1,lblk)
c     &         - xmaxstp*xmaxstp
c          end if


          ! test convergence
          lconv = abs(xval).lt.xthr .or.
     &            xval.lt.0d0 .and. abs(xdamp-xdamp_min).lt.xthr
          if (lconv) exit
          if (xval.gt.0d0) then
            xlow = xdamp
            xvlow = xval
          else
            xhigh = xdamp
            xvhigh = xval
          end if

          if (xlow.gt.xhigh) stop 'AHA!'

          ! get the derivative
          xder =
     &     -2d0*inprdd3(vec1,vec2,lugrvf,lugrvf,ludia,xdamp,-3d0,1,lblk)
          if (ntest.ge.100) then
            write(6,'(x,a,i4,3(x,e15.6))')
     &           'iter, damp, val, der ',iter, xdamp, xval, xder
          end if

          xdamp_old = xdamp
          ! newton step:
          xstep = - xval/xder

          if (xdamp+xstep.ge.xhigh) then
            xdamp = xdamp-xval*(xdamp-xhigh)/(xval-xvhigh)
c            xdamp = (xdamp+xhigh)/2d0    
          else if (xdamp+xstep.lt.xlow) then
            if (xlow.eq.xdamp_min) xvlow=-xval
            xdamp = xdamp-xval*(xdamp-xlow)/(xval-xvlow)
c            xdamp = (xdamp+xlow)/2d0
          else
            xdamp = xdamp + xstep 
          end if

          xval_old = xval

        end do

        if (.not.lconv) then
          write(6,*) 'problem with damping!!'
          stop 'optc_diagp'
        end if
       
        if (ntest.ge.100) write(6,*) 'final damping: ',xdamp

        call dmtvcd2(vec1,vec2,ludia,lugrvf,luamp,-1d0,xdamp,1,inv,lblk)
 
      end if
      
      return
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_conjgrad(itype,ilin,
     &                         lu_corstp,     
     &                         lu_ucrstp,lust_sbsp,nrec,
     &                         lugrvf,lugrvfold,lusig,
     &                         vec1,vec2,nwfpar,iprint)
*----------------------------------------------------------------------*
*
* apply conjugate gradient correction to new search direction:
*
*     |d> = |d_unc> - beta |d_last>
*
*   |d_uncorr> is initially on    lu_ucrstp
*   |d>        will finally be on lu_corstp
*
*   |d_last>   is the last step (record nrec on file lust_sbsp)
*
* beta can be calculated according to the following approximations:
*
*   1) orthodox:
*                   < g - g_last | d_unc >  
*           beta = ------------------------
*                   < g - g_last | d_last>
*
*   2) Polak-Ribiere:
*                    < g - g_last| d_unc >
*           beta = -------------------------
*                     < g_last | d_last >
*
*   3) Fletcher-Reeves:
*                       <  g  | d_unc >
*           beta = ------------------------
*                    < g_last | d_last >
*
*  if ilin==1:
*                   < d_last A | d_unc >
*           beta = ------------------------
*                   < d_last A | d_last >
*
*    where < d_last A | is on lusig
*
*----------------------------------------------------------------------*

      implicit none

* constants
      integer, parameter ::
     &     ntest = 00

* input/output
      integer, intent(in) ::
     &     itype,ilin,
     &     lu_ucrstp,lu_corstp,lust_sbsp,
     &     lugrvf,lugrvfold,lusig,
     &     nrec,nwfpar,iprint
      real(8), intent(inout) ::
     &     vec1(nwfpar), vec2(nwfpar)

* local variables
      integer ::
     &     lblk, luscr, iprintl
      real(8) ::
     &     xnum, xdnm, beta
      
*external functions
      real(8), external ::
     &     inprdd
      integer, external ::
     &     iopen_nus

      lblk = -1
      iprintl = max(iprint,ntest)

      if (ntest.ge.10) then
        write(6,*) 'optc_conjgrad:'
        write(6,*) '=============='
        write(6,*) ' itype, ilin : ', itype, ilin
        write(6,*) ' lu_corstp, lu_ucrstp, lust_sbsp, nrec: ',
     &       lu_corstp, lu_ucrstp, lust_sbsp, nrec
        write(6,*) ' lugrvf, lugrvfold, lusig: ',
     &               lugrvf, lugrvfold, lusig
      end if

c      if (ilin.eq.0.and.(itype.eq.1.or.itype.eq.2)) then
      if (itype.eq.1.or.itype.eq.2) then
        luscr = iopen_nus('DLTGRD')
        call vecsmd(vec1,vec2,1d0,-1d0,lugrvf,lugrvfold,luscr,1,lblk)
        ! calc <dlt g|d_unc>
        xnum = inprdd(vec1,vec2,luscr,lu_ucrstp,1,lblk)
        if (itype.eq.1) then
          call skpvcd(lust_sbsp,nrec-1,vec1,1,lblk)
          call rewino(luscr)
          ! calc <dlt g|d_last>
          xdnm = inprdd(vec1,vec2,luscr,lust_sbsp,0,lblk)
        end if
        call relunit(luscr,'delete')
      end if

c      if (ilin.eq.0.and.(itype.eq.2.or.itype.eq.3)) then
      if (itype.eq.2.or.itype.eq.3) then
        ! calc <g_last|d_last>
        call skpvcd(lust_sbsp,nrec-1,vec1,1,lblk)
        call rewino(lugrvfold)
        xdnm = inprdd(vec1,vec2,lugrvfold,lust_sbsp,0,lblk)
        if (itype.eq.3) then
          ! calc < g | d_unc>
          xnum = inprdd(vec1,vec2,lugrvf,lu_ucrstp,1,lblk)
        end if
      end if

c      if (ilin.eq.1) then
c        ! calc < d_last A | d_unc >
c        xnum = inprdd(vec1,vec2,lusig,lu_ucrstp,1,lblk)
cc???        ! calc < d_unc | d_last>
cc???        call skpvcd(lust_sbsp,nrec-1,vec1,1,lblk)
cc???        call rewino(lu_ucrstp)
cc???        xnum = xnum - inprdd(vec1,vec2,lu_ucrstp,lust_sbsp,0,lblk)
c        ! calc < d_last A | d_last>
c        call skpvcd(lust_sbsp,nrec-1,vec1,1,lblk)
c        call rewino(lusig)
c        xdnm = inprdd(vec1,vec2,lusig,lust_sbsp,0,lblk)
c      end if

      if (ntest.ge.10) then
        write(6,*) 'ilin  = ', ilin
        write(6,*) 'itype = ', itype
        write(6,*) ' xnum = ', xnum
        write(6,*) ' xdnm = ', xdnm
      end if

      beta =  xnum / xdnm

c??
      if (itype.eq.2.or.itype.eq.3) beta = -beta

      if (iprintl.ge.1) then
        write (6,'(x,a)') 'conjugate gradient correction:'
        if (ilin.eq.1)
     &       write(6,'(4x,a,e10.4)') '(lin. equations):  beta = ', beta
        if (ilin.eq.0.and.itype.eq.1)
     &       write(6,'(4x,a,e10.4)') '(orthodox):        beta = ', beta
        if (ilin.eq.0.and.itype.eq.2)
     &       write(6,'(4x,a,e10.4)') '(Polak-Ribiere):   beta = ', beta
        if (ilin.eq.0.and.itype.eq.3)
     &       write(6,'(4x,a,e10.4)') '(Fletcher-Reeves): beta = ', beta
      end if

      ! obtain corrected step:
      call rewino(lu_ucrstp)
      call rewino(lu_corstp)
      call skpvcd(lust_sbsp,nrec-1,vec1,1,lblk)
      call vecsmd(vec1,vec2,1d0,-beta,lu_ucrstp,lust_sbsp,
     &            lu_corstp,0,lblk)


      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_linesearch(ilsrch,ilin,ivar,iprecnd,ipass,
     &         alpha,dalp,energy,de_pred,trrad,xnstp,
     &         vec1,vec2,
     &         lu_corstp,lu_ucrstp,lusig,ludia,iprint)
*----------------------------------------------------------------------*
*
* purpose: obtain step length from an approximate second-order model
*          along the given search direction.
*
* One point model:
*
*  The local model is
*
*   E(alpha) = E_0 - alpha <g|d> + 1/2 alpha^2 <d|A|d>
*
*  where the diagonal preconditioner C is used as approximation to A     
*
*
*                     < d | g >            < d | g >
*        alpha = - ---------------  ~ - ---------------
*                   < d | A | d >        < d | C | d >
*
*  if ilin==1: exact < d A |  on lusig
*
*     for prec. steepest descend alpha then always is 1.0
*
*  | d > is always on lu_corstp
*  | g > is found on lu_ucrstp
*
*
*  Two point model: (needs two passes)
*
* We predict an initial alpha as above, and set up the local function 
* (along our new step direction):
*
*   E(alpha) = E_0 + alpha <g|d> + 1/2 alpha^2 h
*
* to find the optimum alpha as
*
*   alpha = |g|/h
*
* where h can be estimated from our extra point E(alpha) as
*
*   h = 2 (E(alpha) - E_0 - alpha <g|d>)/alpha^2
*     = 2 (E(alpha) - E_0 - dE(pred) ) / alpha^2
*
* dE(pred) being the energy estimate from the linear model 
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00

* input/output
      integer, intent(in) ::
     &     ilsrch,ilin,ivar,iprecnd,ipass,
     &     lu_corstp,lu_ucrstp,ludia,lusig,
     &     iprint
      real(8), intent(in) ::
     &     energy, xnstp, trrad
      real(8), intent(out) ::
     &     alpha, de_pred, dalp
      real(8), intent(inout) ::
     &     vec1(*), vec2(*)

* local save variables
      integer, save ::
     &     last_ipass
      real(8), save ::
     &     de_lin, de_qdr, e_old, alpha_old, alpha_max, xdg
* local
      integer ::
     &     lblk, inv, luscr, iprintl
      real(8) ::
     &     xdad, xdad2, xh, de

* external functions
      real(8), external ::
     &     inprdd
      integer, external ::
     &     iopen_nus

      lblk = -1

      iprintl = max(ntest,iprint)

      if (iprintl.ge.1) write (6,*) 'line-search:'

      if (ntest.ge.10) then
        write(6,*) ' entered optc_linesearch with:'
        write(6,*) '  ivar,iprecnd, ipass :  ',ivar,iprecnd, ipass
        write(6,*) '  energy,trrad,xnstp  :  ',energy,trrad,xnstp
        write(6,*) '  lu_corstp,lu_ucrstp,ludia: ',
     &       lu_corstp,lu_ucrstp,ludia
      end if

      if (ipass.eq.1) then

        last_ipass = 1
        
        luscr = iopen_nus('DADSCR')

        ! calc < d | g > or < d | C^-1 g >, respectively
        xdg = inprdd(vec1,vec2,lu_corstp,lu_ucrstp,1,lblk)

        inv = 0
        ! calc < d | A | d > 
        call dmtvcd(vec1,vec2,ludia,lu_corstp,luscr,0d0,1,inv,lblk)
        xdad = inprdd(vec1,vec2,lu_corstp,luscr,1,lblk)
        alpha = 1d0
        if (ilsrch.gt.0)
     &     alpha = - xdg / xdad  ! interestingly, a factor of 2 is helpful
                                 ! with conjugate gradient methods ...
        ! get max alpha

        if (ntest.ge.10) then
          write(6,*) ' xdad, xdg: ', xdad, xdg
        end if

c        if (alpha.lt.0d0) stop 'PANIC: alpha.lt.0d0!!!'
        if (alpha.lt.0d0)write(6,*) '>>> WARNING: alpha.lt.0d0!!!'

        if (ilsrch.eq.2) then
          if (iprintl.ge.2)
     &         write(6,*) ' reduced trust radius for trial step (.75)'
          alpha_max = 0.75*trrad/xnstp
        else
          alpha_max = trrad/xnstp
        end if
        if (abs(alpha).gt.alpha_max) then
          if (iprintl.ge.2) then
            write (6,*) ' alpha = ',alpha 
            write (6,*) ' norm step = ',alpha*xnstp
            write (6,*) ' trust radius = ',trrad
          end if
          alpha = alpha_max*sign(1d0,alpha)
        end if
c TEST
        alpha = 1d0
        write(6,*) 'alpha set to 1d0'
c
        if (ivar.eq.1) then
          alpha_old = alpha
          e_old  = energy
          de_lin = alpha*xdg
          de_qdr = alpha*xdg + .5d0 * alpha*alpha*xdad
          de_pred = de_qdr
        end if

        if (iprintl.ge.1) then
          write (6,*)   '  alpha from one-point model: ',alpha
          if (ivar.eq.1) then
            write (6,*) '  predicted energy change:    ',de_pred
          end if
        end if

        call relunit(luscr,'delete')

      else if (ipass.eq.2) then
        if (ivar.eq.0) then
          write (6,*)
     &    'Error: 2 point line-search entered for non-variational model'
          stop 'line-search'
        end if
        if (last_ipass.ne.1) then
          write (6,*)
     &    'Error: pass 2 without previous pass 1'
          stop 'line-search'
        end if

        if (ntest.ge.10) then
          write(6,*) 'energy, e_old, de_lin, alpha_old: ',
     &         energy, e_old, de_lin, alpha_old
        end if

        ! values have hopefully been set in previous pass
        de = energy-e_old
        if (abs(de).lt.1d-12) then
          write(6,*) 'Energy difference too small : ', de
          write(6,*) 'reverting to old alpha'
          alpha = alpha_old
        else
          if (de.lt.0d0) then
            alpha_max = min(alpha_max,5d0*alpha_old)
          else
            alpha_max = 0.5*alpha_old
          end if
          if (ntest.ge.10) then
            write(6,*) 'alpha_max set to ',alpha_max
          end if
          xh = 2d0*(de - de_lin) / (alpha_old*alpha_old)
          if (ntest.ge.10) then
            write(6,*) ' xh, xdg: ', xh, xdg
          end if
          
          if (xh.le.0d0) then
            alpha = alpha_max
            if (ntest.ge.5) then
              write(6,*)'h <= 0, so we step to the trust radius border'
              write(6,*) 'alpha opt   = ',alpha
            end if
          else
            alpha = -xdg/xh
c            if (alpha.lt.0d0) stop 'PANIC: alpha.lt.0d0!!! (2)'
            if (alpha.lt.0d0)
     &           write(6,*) '>>>WARNING: alpha.lt.0d0!!! (2)'
            if (ntest.ge.5) then
              write(6,*) 'alpha opt    = ',alpha
              write(6,*) 'norm of step = ',alpha*xnstp
              write(6,*) 'trust rad    = ',trrad
            end if
          end if

        end if

        alpha = sign(min(alpha_max,abs(alpha)),alpha)
        if (ntest.ge.5) then
          write(6,*) 'alpha taken  = ',alpha
        end if

        de_pred = alpha * xdg + 0.5d0*xh*alpha*alpha

        if (iprintl.ge.1) then
          write (6,*) '  alpha from two-point model: ',alpha
          write (6,*) '  predicted energy change:    ',de_pred
        end if

        ! return also difference to old alpha
        dalp = alpha - alpha_old
        
        if (ntest.ge.10) then
          write(6,*) 'returned alpha, alpha-alpha_old ',alpha,dalp
        end if

      end if

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_diis(itype,thrsh,nvec_sbsp,
     &     navec,maxvec,
     &     nadd_in,navec_last,alpha_last,
     &     lu_step,lu_corstep,lu_sbsp,lu_sbsp2,
     &     luamp,lugrvf,
     &     bmat,scr,vec1,vec2,iprint)
*----------------------------------------------------------------------*
*
*  perform DIIS extrapolation of wavefunction vectors.
*
*  itype = 1    take actual t_k+1-t_k as error vector for t_k+1
*               and interpolate in (t_k) basis
*               (*not* recommended!!)
*  
*  itype = 2    take perturbation steps delta t_k (= preconditioned 
*               vector function) as error vectors (on lu_sbsp)
*               and interpolate in (t_k + delta t_k) basis
*
*  itype = 3    take vector function as error vectors (on lu_sbsp)
*               and interpolate in (t_k + delta t_k) basis
*
*  itype = 4    take vector function as error vectors (on lu_sbsp)
*               and interpolate (t_k) and (Omega_k) separately
*
*----------------------------------------------------------------------*

      implicit none

* constants
      integer, parameter ::
     &     ntest = 00

* input/output
      integer, intent(in) ::
     &     itype,
     &     nvec_sbsp,      ! number of vectors in subspace
     &     maxvec,         ! max number of vectors allowed
     &     nadd_in,        ! number of new vectors on lu_step
     &     lu_step, lu_corstep, lu_sbsp,
     &     lu_sbsp2, luamp, lugrvf,       ! unit numbers
     &     iprint                         ! print level
c      integer, intent(out) ::
c     &     ndel_out        ! request to sbspman: delete these vectors
      integer, intent(inout) ::
     &     navec,          ! number active vectors in subspace
     &                     ! (skip the first nvec-navec vectors)
     &     navec_last      ! number active vectors in subspace in 
                           ! previous call
      real(8), intent(in) ::
     &     alpha_last,     ! scaling factor for previous vector
     &     thrsh           ! threshold for accepting DIIS step
      real(8), intent(inout) ::
     &     bmat(*), scr(*), vec1(*), vec2(*) ! DIIS-matrix and scratch

* local save
c      integer, save ::

* local O(N) scratch
      integer ::
     &     kpiv(navec+1)
      real(8) ::
     &     scrvec(navec+1)

* local
      logical ::
     &     lincore, again
      integer ::
     &     lblk, iprintl, navec_l,
     &     nskip, nold, ndel, ii, jj, iioff, iioff2,
     &     lu_sbsp_un, lu_curerr, lust, itype_l
      real(8) ::
     &     cond, xcorsum

      iprintl = max(iprint,ntest)
      lblk = -1
      lincore = .false.

      if (iprintl.ge.5) then
        write(6,*) 'DIIS-extrapolation'
        write(6,*) '=================='
        if (itype.eq.1) then
          write(6,*) ' error vector: c^(k+1)-c^(k)'
        else if (itype.eq.2) then
          write(6,*) ' error vector: c^(k+1)_pert'
        else if (itype.eq.3.or.itype.eq.4) then
          write(6,*) ' error vector: grad_c^(k+1)'
        end if
      end if

      if (ntest.ge.10) then
        write(6,*) 'on entry in optc_diis:'
        write(6,*) 'nvec_sbsp, navec, maxvec, nadd_in, navec_last: ',
     &              nvec_sbsp, navec, maxvec, nadd_in, navec_last
        write(6,*) 'lu_step, lu_corstep, lu_sbsp: ',
     &              lu_step, lu_corstep, lu_sbsp
      end if

      if (itype.lt.1.or.itype.gt.4) then
        write(6,*) 'DIIS: unknown itype = ',itype
        stop 'DIIS arguments'
      end if


c      if (nvec_sbsp+nadd_in.ne.navec) then
c        write(6,*) 'Panic 1: DIIS subspace dimensioning inconsistent!'
c        stop 'DIIS dimensions'
c      end if
      if (navec.gt.maxvec) then
        write(6,*) 'Panic 2: DIIS subspace dimensioning inconsistent!'
        stop 'DIIS dimensions'
      end if

      nskip = nvec_sbsp+nadd_in-navec
                                 ! vectors that are to be skipped on Bmat
                                 ! (after deleting), normally 0
      ndel  = navec_last+nadd_in-navec
                                 ! vectors that are on Bmat but are no more
                                 ! needed
      nold  = navec_last-ndel    ! vectors that are already on Bmat and are
                                 ! still needed

      navec_last = navec         ! save for next iteration

      if (nold.lt.0.or.nskip.lt.0) then
        write(6,*) 'inconsistent dimensions in optc_diis!'
        write(6,*) 'nold , nskip : ', nold, nskip
        stop 'DIIS dimensions'
      end if

      if (ntest.ge.10) then
        write(6,*) 'vectors to be deleted in B mat (ndel): ',ndel
        write(6,*) 'old vectors on B matrix: (nold) ',nold
        write(6,*) 'skipped vectors on B mat:(nskip)',nskip
      end if

* vectors deleted? move B matrix by the corresp. number of col.s and rows
      if (ndel.gt.0) then
c        if (ndel.gt.nskip) then
c          write(6,*) 'Panic 2: DIIS subspace dimensioning inconsistent!'
c          stop 'DIIS dimensions'
c        end if
        ! nold is the future dimension of the already initialized
        ! DIIS matrix
        ! we now move all entries according to
        !   B(i,j) = B(i+ndel,j+ndel), where ij actually refer to a
        !                              upper triangular matrix
        if (ntest.ge.100) then
          write(6,*) 'DIIS matrix before deleting:'
          call prtrlt(bmat,nold+ndel)
        end if
        
        do ii = 1, nold
          iioff  = ii*(ii-1)/2
          iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
          do jj = 1, ii
            bmat(iioff+jj) = bmat(iioff2+jj)
          end do
        end do
        if (ntest.ge.100) then
          write(6,*) 'DIIS matrix after deleting:'
          call prtrlt(bmat,nold)
        end if

      end if

* modify most recent vector (if alpha was != 1.0)
      if (itype.eq.1.and.nold.gt.0.and.alpha_last.ne.1d0) then
        iioff = nold*(nold-1)/2
        do ii = 1, nold-1
          bmat(iioff+ii) = alpha_last*bmat(iioff+ii)
        end do
        bmat(iioff+nold) = alpha_last*alpha_last*bmat(iioff+nold)
        if (ntest.ge.100) then
          write(6,*) 'DIIS matrix after modifying last vector:'
          call prtrlt(bmat,nold)
        end if
      end if

* add new vectors
      ndel = 0
      if (itype.eq.1.or.itype.eq.2) then
        lu_curerr = lu_step
      else
        lu_curerr = lugrvf
      end if
      call optc_diis_bmat(bmat,nadd_in,nold,ndel,nskip,
     &     lu_sbsp,lu_curerr,vec1,vec2,lincore)

      again = .true.
      navec_l = navec
      do while(again)
* transfer to scratch array and augment
        do ii = 1, navec_l
          iioff  = ii*(ii-1)/2
          iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
          do jj = 1, ii
            scr(iioff+jj) = bmat(iioff2+jj)
          end do
        end do
        iioff = navec_l*(navec_l+1)/2
        do jj = 1, navec_l
          scr(iioff+jj) = -1d0
        end do
        scr(iioff+navec_l+1) = 0d0
        
        if (ntest.ge.100) then
          write(6,*) 'Augmented matrix passed to dspco:'
          call prtrlt(scr,navec_l+1)
        end if
* solve DIIS problem to obtain new weights
        ! first factorize the DIIS matrix and get its condition number
        call dspco(scr,navec_l+1,kpiv,cond,scrvec)
c        call dppco(scr,navec+1,cond,scrvec,info)
        if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &       navec_l,cond
        if (ntest.ge.100) then
          write(6,*) 'Factorized matrix from dspco:'
          call prtrlt(scr,navec_l+1)
          write(6,*) 'Pivot array:'
          call iwrtma(kpiv,1,navec_l+1,1,navec_l+1)
        end if

        again = .false.
c        if (itype.eq.1) again = cond.lt.1d-14

        ! if everything is fine ...
        if (.not.again) then
          ! ... we set up the RHS vector and solve the DIIS equations
          scrvec(1:navec_l) =  0d0
          scrvec(navec_l+1) = -1d0
          call dspsl(scr,navec_l+1,kpiv,scrvec)
c          call dppsl(scr,navec_l+1,kpiv,scrvec)
          if (ntest.ge.10) then
            write(6,*) 'result of dspsl:'
            write(6,*) 'w:', scrvec(1:navec_l)
            write(6,*) 'l:', scrvec(navec_l+1)
          end if
        else
          ! ... else we have to reduce the number of active dimensions
          nskip = nskip+1
          ndel  = ndel +1
          navec_l = navec_l-1
        end if
* analyze solution and event. request deletion of vectors
        if (.not.again) then
          xcorsum = 0d0
          do ii = 1, navec_l-1
            xcorsum = xcorsum + abs(scrvec(ii))
          end do
          if (xcorsum/abs(scrvec(navec_l)).gt.1.2d0) then
            again = .true.
            nskip = nskip + 1
            ndel  = ndel  + 1
            navec_l = navec_l -1
          end if
        end if

      end do

* get new step according to weigths
      if (itype.eq.1.or.itype.eq.4) lu_sbsp_un = lu_sbsp
      if (itype.eq.2.or.
     &    itype.eq.3) lu_sbsp_un = lu_sbsp2
      if (itype.eq.1.or.itype.eq.2.or.itype.eq.3) lust = lu_step
      if (itype.eq.4) lust = lugrvf
      call optc_diis_nstp(itype,thrsh,scrvec,navec_l,nskip,lu_sbsp_un,
     &                    lust,lu_corstep,luamp,vec1,vec2,lincore)

      if (itype.eq.4) then
        ! extrapolate here new internal step
        itype_l = 5
        lust = 0
        call optc_diis_nstp(itype_l,thrsh,scrvec,navec_l,nskip,lu_sbsp2,
     &       lust,lu_step,luamp,vec1,vec2,lincore)

      end if

      if (itype.eq.1) then
* reupdate DIIS B-matrix with actual new step
        call optc_diis_bmat(bmat,nadd_in,nold,ndel,nskip,
     &       lu_sbsp,lu_corstep,vec1,vec2,lincore)
      end if

* return the actual last DIIS dimension
      navec = navec_l

      return
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_diis_bmat(bmat,nadd,nold,nskipb,nskip,
     &                          lu_sbsp,lu_step,
     &                          vec1,vec2,lincore)
*----------------------------------------------------------------------*
*
*  Update the DIIS B matrix ( B = <e_i|e_j> ) with nadd new vectors
*  from lu_step.
*  If lincore.eq..true., the vectors vec1 and vec2 can hold a complete
*  vector from disc.
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     lincore
      integer, intent(in) ::
     &     nadd,nold,lu_sbsp,lu_step,nskip,nskipb
      real(8), intent(inout) ::
     &     bmat(*),vec1(*),vec2(*)

      real(8) ::
     &     bij, bii
      integer ::
     &     ii, jj, iioff, iirec, lblk

      real(8), external ::
     &     inprdd

      lblk = -1
      if (lincore) stop 'no lincore'
      if (nadd.gt.1.and..not.lincore) then
        write(6,*) 'optc_diis not completely read to add more than '//
     &             'vector at a time'
        stop 'DIIS: nadd.gt.1'
      end if

      if (ntest.ge.10) then
        write(6,*) 'lu_sbsp, lu_step: ',lu_sbsp, lu_step
        write(6,*) 'updating B matrix with ', nadd,' vector(s)'
        write(6,*) 'old dimension was ',nold
      end if

      do ii = nold + 1, nold+nadd
        iirec = ii-nold
        iioff = ii*(ii-1)/2
        call rewino(lu_sbsp)
        if (nskip.gt.0) call skpvcd(lu_sbsp,nskip,vec1,1,lblk)
        do jj = 1, nskipb
          bmat(iioff+jj) = 0d0
        end do
c        if (lincore) call vec_from_disc(vec1,)
        do jj = nskipb + 1, nold
          call rewino(lu_step)
          if (iirec.gt.1) call skpvcd(lu_step,iirec-1,vec1,1,lblk)
          bij = inprdd(vec1,vec2,lu_sbsp,lu_step,0,lblk)
          bmat(iioff+jj) = bij
        end do
        ! well, here adding more than one vector becomes a bit tricky
        ! unless we can do that incore. we leave it as is for a while
        call rewino(lu_step)
        if (iirec.gt.1) call skpvcd(lu_step,iirec-1,vec1,1,lblk)
        bii = inprdd(vec1,vec2,lu_step,lu_step,0,lblk)
        bmat(iioff+nold+1) = bii
      end do
*
      if (ntest.ge.100) then
        write(6,*) 'DIIS matrix after adding new vector:'
        call prtrlt(bmat,nold+nadd)
      end if

      return
      end

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_diis_nstp(itype,thrsh,wvec,navec,nskip,lusbsp,
     &                          lu_step,lu_corstep,luamp,
     &                          vec1,vec2,lincore)
*----------------------------------------------------------------------*
*
* obtain new step from the DIIS-weigths on wvec, the old steps on
* lu_sbsp and the current perturbation estimate on lu_step
*
*  if |diis_correction|/|perturbation_step| > thrsh, the DIIS is rejected
*
*----------------------------------------------------------------------*

      implicit none
* constants
      integer, parameter ::
     &     ntest = 00

* input/output
      logical, intent(in) ::
     &     lincore
      integer, intent(in) ::
     &     itype,navec, nskip,
     &     lusbsp,lu_step,lu_corstep,luamp
      real(8), intent(in) ::
     &     thrsh,
     &     wvec(navec), vec1(*), vec2(*)

      integer ::
     &     ii, jj, lblk, ndim,
     &     luscr
      real(8) ::
     &     swvec(1:navec), xfac, xnrm

      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd

      if (lincore) stop 'no lincore in optc_diis_newstp' 

      lblk = -1
      luscr = iopen_nus('DIISCR')

      if (itype.eq.1.or.itype.eq.5) then
        swvec(1:navec) = 0d0
        do jj = 1, navec-1
          do ii = jj,navec      ! jj+1 -> jj ?
            swvec(jj) = swvec(jj) + wvec(ii)
          end do
          swvec(jj) = swvec(jj) - 1d0
        end do
        swvec(navec) = wvec(navec)
        if (itype.eq.5) then
          swvec(navec) = swvec(navec) - 1d0
          do ii = 1, navec
            swvec(ii) = swvec(ii+1)
          end do
        end if
      else if (itype.eq.2.or.itype.eq.3.or.itype.eq.4) then
        swvec(1:navec) = wvec(1:navec)
      else
        write(6,*) 'unexpected itype in optc_diis_nstep: ',itype
        stop 'optc_diis_nstep'
      end if

      if (ntest.ge.10) then
        write(6,*)
     &       'optc_diis_nstep: weight vector for new step: (itype = ',
     &       itype,')'
        ndim = navec
        if (itype.eq.5) ndim = navec-1
        call wrtmat(swvec,ndim,1,ndim,1)
      end if

      call skpvcd(lusbsp,nskip,vec1,1,lblk)
      call rewino(lu_corstep)
      call rewino(luscr)

      if (navec-1.gt.0) then
        call mvcsmd(lusbsp,swvec,luscr,lu_corstep, ! lu_corstep is only scratch
     &            vec1,vec2,navec-1,0,lblk)
        if (lu_step.ne.0.and.itype.ne.5) then
          call vecsmd(vec1,vec2,1d0,swvec(navec),luscr,lu_step,
     &         lu_corstep,1,lblk)
        else
          call copvcd(luscr,lu_corstep,vec1,1,lblk)
        end if
        if (itype.eq.2.or.itype.eq.3) then
* to be consistent, we add only (w_n - 1)times the current amplitudes
          call vecsmd(vec1,vec2,1d0,swvec(navec)-1d0,
     &         lu_corstep,luamp,
     &         luscr,1,lblk)
* too much copying, but for the moment it is OK
          call copvcd(luscr,lu_corstep,vec1,1,lblk)
        end if
      else
        if (lu_step.ne.0.and.itype.ne.5)
     &     call sclvcd(lu_step,lu_corstep,swvec(navec),vec1,1,lblk)
      end if

      ! get norm of DIIS correction:
      if (lu_step.ne.0.and.itype.ne.5.and.itype.ne.4) then
        call vecsmd(vec1,vec2,1d0,-1d0,lu_corstep,lu_step,luscr,1,lblk)
        xnrm = sqrt(inprdd(vec1,vec2,luscr,luscr,1,lblk))
        if (ntest.ge.10) write(6,*) '|DIIS correction:| = ',xnrm
        if (xnrm.gt.thrsh) then
          if (ntest.ge.10) write(6,*) 'DIIS step not accepted!'
          call copvcd(lu_step,lu_corstep,vec1,1,lblk)
        end if
      end if

      call relunit(luscr,'delete')

      return
      end 
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja(itype,thrsh,nstdim,ngvdim,
     &           navec,maxvec,
     &           nadd,navec_last,
     &           lugrvf,ludia,trrad,lu_pertstp,lu_newstp,
     &           xdamp,xdamp_last,
     &           lust_sbsp,lugv_sbsp,
     &           umat,scr,vec1,vec2,iprint)
*----------------------------------------------------------------------*
*
* subspace jacobian: A = A_0(1-P) + A_ex P
*
*  if xdamp is 0d0, we try the direct inversion procedure; if that
*  step is too large, we go over to the iterative solver to find the
*  optimal damping (where we end up directly, if xdamp was != 0d0) 
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00
      
* input/output
      integer, intent(in) ::
     &     itype,                    ! type (s.above)
     &     nstdim,ngvdim,            ! # records on sbsp files
     &     maxvec,                   ! max dimension of sbsp
     &     nadd,                     ! number of new vectors (usually 1)
     &     lugrvf,ludia,             ! current gradient, diagonal Hess/Jac.
     &     lu_pertstp,               ! step from diagonal Hess/Jac.
     &     lu_newstp,                ! new step (output)
     &     lust_sbsp,lugv_sbsp,      ! sbsp files
     &     iprint                    ! print level
      integer, intent(inout) ::
     &     navec,                    ! current dimension of sbsp
     &     navec_last                ! previous sbsp dimension
      real(8), intent(in) ::
     &     thrsh,                    ! thresh for accepting low-rank correction
     &     trrad,                    ! trust radius for step
     &     xdamp_last                ! previous damping
      real(8), intent(out) ::
     &     xdamp                     ! current damping
      real(8), intent(inout) ::
     &     umat(*),                  ! low rank Hessian/Jacobian
     &     scr(*), vec1(*), vec2(*)  ! scratch

* local O(N) scratch
      integer ::
     &     kpiv(navec)
      real(8) ::
     &     scrvec(navec), scrvec2(navec)

* local
      logical ::
     &     lincore, again, accept, converged
      integer ::
     &     lblk, iprintl, navec_l, isym, job,
     &     nskipst, nskipgv, nnew, nold, nold2, ndel,
     &     ii, jj, iioff, iioff2, irhsoff,
     &     lu_sbsp_un, lu_curerr,
     &     luvec, luAvec, nadd_l, nold_l, nskip_l, ludum,
     &     iter
      real(8) ::
     &     cond, xcorsum
      integer, external ::
     &     iopen_nus

      iprintl = max(iprint,ntest)
      lblk = -1
      lincore = .false.
      if (itype.eq.1) isym = 0  ! non-symmetric jacobian
      if (itype.eq.2) isym = 1  ! symmetric jacobian

      if (iprintl.ge.5) then
        write(6,*) 'Subspace-Hessian/Jacobian'
        write(6,*) '========================='
        if (itype.eq.1) then
          write(6,*) ' asymmetric update'
        else if (itype.eq.2) then
          write(6,*) ' symmetric update'
        end if
      end if

      if (ntest.ge.10) then
        write(6,*) 'on entry in optc_sbspja:'
        write(6,*)
     &       'nstdim, ngvdim, navec, maxvec, nadd, navec_last: ',
     &        nstdim, ngvdim, navec, maxvec, nadd, navec_last
        write(6,*)
     &       'lugrvf, lu_newstp, lu_pertstp, lust_sbsp, lugv_sbsp: ',
     &        lugrvf, lu_newstp, lu_pertstp, lust_sbsp, lugv_sbsp
      end if

      if (itype.lt.1.or.itype.gt.2) then
        write(6,*) 'SBSP hessian/jacobian: unknown itype = ',itype
        stop 'SBSPJA arguments'
      end if

      if (navec.gt.maxvec) then
        write(6,*)
     &       'Panic 2: Hess./Jac. subspace dimensioning inconsistent!'
        stop 'SBSP dimensions'
      end if

      nskipst = nstdim-navec
      nskipgv = ngvdim-navec
      ndel  = navec_last+nadd-navec
      nold  = navec_last-ndel

      navec_last = navec         ! save for next iteration

      if (nold.lt.0.or.nskipst.lt.0.or.nskipgv.lt.0) then
        write(6,*) 'inconsistent dimensions in optc_sbspja!'
        write(6,*) 'nold , nskipst, nskipgv : ', nold, nskipst, nskipgv
        stop 'SBSPJA dimensions'
      end if

      if (ntest.ge.10) then
        write(6,*) 'vectors to be deleted in U mat (ndel): ',ndel
        write(6,*) 'old vectors on U matrix: (nold) ',nold
        write(6,*) 'skipped vectors on step file:(nskipst)',nskipst
        write(6,*) 'skipped vectors on grad file:(nskipgv)',nskipgv
      end if

* vectors deleted? move B matrix by the corresp. number of col.s and rows
      if (ndel.gt.0) then
        if (isym.eq.1) then ! symmetric U matrix
        ! nold is the future dimension of the already initialized
        ! low-rank matrix
        ! we now move all entries according to
        !   U(i,j) = U(i+ndel,j+ndel), where ij actually refer to a
        !                              upper triangular matrix
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix before deleting:'
            call prtrlt(umat,nold+ndel)
          end if
        
          do ii = 1, nold
            iioff  = ii*(ii-1)/2
            iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
            do jj = 1, ii
              umat(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix after deleting:'
            call prtrlt(umat,nold)
          end if
        else ! non-symmetric (and therefore full) matrix
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix before deleting:'
            call wrtmat2(umat,nold+ndel,maxvec,nold+ndel,maxvec)
          end if
        
          do ii = 1, nold
            iioff  = (ii-1)*maxvec
            iioff2 = (ii+ndel-1)*maxvec + ndel
            do jj = 1, nold
              umat(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix after deleting:'
            call wrtmat2(umat,nold,maxvec,nold,maxvec)
          end if

        end if ! isym

      end if ! ndel.gt.0
      
* add new vectors
      ndel = 0
      nnew = nadd
      nold2 = nold
      call optc_sbspja_umat(itype,
     &     umat,nnew,nold2,ndel,maxvec,
     &     nskipst,nskipgv,
     &     lust_sbsp,lugv_sbsp,ludia,xdamp,
     &     vec1,vec2,lincore)

* get projection of perturbation step in current subspace
      call optc_sbspja_prjpstp(scrvec2,navec,nskipst,
     &                         lu_pertstp,lust_sbsp,
     &                         vec1,vec2,lincore)

      navec_l = navec

* direct-inversion section:
* =========================
      again = .true.
      do while(again)
* transfer to scratch array
        if (isym.eq.1) then
          do ii = 1, navec_l
            iioff  = ii*(ii-1)/2
            iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
            do jj = 1, ii
              scr(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix passed to dspco:'
            call prtrlt(scr,navec_l)
          end if
        else
          do ii = 1, navec_l
            iioff = (ii-1)*navec_l
            iioff2 = (ii+ndel-1)*maxvec + ndel
            do jj = 1, navec_l
              scr(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix passed to dgeco:'
            call wrtmat2(scr,navec_l,navec_l,navec_l,navec_l)
          end if
        end if

* invert low-rank matrix:
* factorize and get condition
        if (isym.eq.1) then
          call dspco(scr,navec_l,kpiv,cond,scrvec)
          if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &         navec_l,cond
          if (ntest.ge.100) then
            write(6,*) 'Factorized matrix from dspco:'
            call prtrlt(scr,navec_l)
            write(6,*) 'Pivot array:'
            call iwrtma(kpiv,1,navec_l,1,navec_l)
          end if
        else
          call dgeco(scr,navec_l,navec_l,kpiv,cond,scrvec)
          if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &         navec_l,cond
          if (ntest.ge.100) then
            write(6,*) 'Factorized matrix from dgeco:'
            call wrtmat2(scr,navec_l,navec_l,navec_l,navec_l)
            write(6,*) 'Pivot array:'
            call iwrtma(kpiv,1,navec_l,1,navec_l)
          end if
        end if ! isym
        
        again = cond.lt.1d-14

        if (.not.again) then
          irhsoff = navec - navec_l + 1
          scrvec(1:navec_l) = scrvec2(irhsoff:irhsoff-1+navec_l)
* if OK, solve U x = rhs
          if (ntest.ge.100) then
            write(6,*) 'RHS: ',scrvec(1:navec_l)
          end if
          if (isym.eq.1) then
            call dspsl(scr,navec_l,kpiv,scrvec)
          else
            job = 0
            call dgesl(scr,navec_l,navec_l,kpiv,scrvec,job)
          end if
          if (ntest.ge.100) then
            write(6,*) 'Result for x:'
            call wrtmat(scrvec,navec_l,1,navec_l,1)
          end if
        end if

        if (.not.again) then
* get low-rank contribution to new vector
          call optc_sbspja_lrstep(itype,thrsh,accept,scrvec,navec_l,
     &         nskipst,nskipgv,
     &         lust_sbsp,lugv_sbsp,ludia,xdamp,
     &         lu_pertstp,lu_newstp,
     &         vec1,vec2,lincore)
* and analyze it
          if (.not.accept.and.navec_l.gt.1) then
            again = .true.
          else if (.not.accept) then
            navec_l = 0
          end if
        end if

        if (again) then
          nskipst = nskipst+1
          nskipgv = nskipgv+1
          ndel = ndel+1
          navec_l = navec_l-1
        end if

      end do

* return the actual last subspace dimension
      navec = navec_l

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja_new(itype,thrsh,nstdim,ngvdim,
     &           navec,maxvec,
     &           nadd,navec_last,
     &           lugrvf,ludia,trrad,lu_pertstp,lu_newstp,
     &           xdamp,xdamp_last,
     &           lust_sbsp,lugv_sbsp,
     &           umat,scr,vec1,vec2,iprint)
*----------------------------------------------------------------------*
*
* subspace jacobian: A(xdamp) = A_0(1-P) + A_ex P + xdamp*1
*
*  xdamp is adjusted according to the trust radius on trrad
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00, mxd_iter = 100
      real(8), parameter ::
     &     thrdamp = 1d-4, xinc = 1d-5, xfailfac = 1d5,
     &     xmxdxdamp = 0.1d0
      
* input/output
      integer, intent(in) ::
     &     itype,                    ! type (s.above)
     &     nstdim,ngvdim,            ! # records on sbsp files
     &     maxvec,                   ! max dimension of sbsp
     &     nadd,                     ! number of new vectors (usually 1)
     &     lugrvf,ludia,             ! current gradient, diagonal Hess/Jac.
     &     lu_pertstp,               ! step from diagonal Hess/Jac.
     &     lu_newstp,                ! new step (output)
     &     lust_sbsp,lugv_sbsp,      ! sbsp files
     &     iprint                    ! print level
      integer, intent(inout) ::
     &     navec,                    ! current dimension of sbsp
     &     navec_last                ! previous sbsp dimension
      real(8), intent(in) ::
     &     thrsh,                    ! thresh for accepting low-rank correction
     &     trrad,                    ! trust radius for step
     &     xdamp_last                ! previous damping
      real(8), intent(out) ::
     &     xdamp                     ! current damping
      real(8), intent(inout) ::
     &     umat(*),                  ! low rank Hessian/Jacobian
     &     scr(*), vec1(*), vec2(*)  ! scratch

* local O(N) scratch
      integer ::
     &     kpiv(navec)
      real(8) ::
     &     scrvec(navec), scrvec2(navec), xd(mxd_iter), xv(mxd_iter)

* local
      logical ::
     &     lincore, again, accept, converged
      integer ::
     &     lblk, iprintl, navec_l, isym, job,
     &     nskipst, nskipgv, nnew, nold, nold2, ndel,
     &     ii, jj, iioff, iioff2, irhsoff,
     &     lu_sbsp_un, lu_curerr,
     &     luvec, luAvec, luvec_sbsp, luAvec_sbsp, lures,
     &     nadd_l, nold_l, nskip_l, ludum, iter, maxiter,
     &     kend, ksmat, ksred, kared, kevec, kscr,
     &     itype_evp, itask(2), ntask, ndum,
     &     ixd_iter,isubcnt
      real(8) ::
     &     cond, xcorsum, xlam, xresnrm, fac, thrsh_l, xnrm,
     &     xdiff, xgr, xhs, f0, f1, f2, a1, a2, xdis, dxdamp, xmxdxdmp2,
     &     xval, xlow, xhigh, xvlow, xvhigh, xdamp_min
      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd, fdmnxd

      iprintl = max(iprint,ntest)
      lblk = -1
      lincore = .false.
      if (itype.eq.1) isym = 0  ! non-symmetric jacobian
      if (itype.eq.2) isym = 1  ! symmetric jacobian

      if (iprintl.ge.5) then
        write(6,*) 'Subspace-Hessian/Jacobian'
        write(6,*) '========================='
        if (itype.eq.1) then
          write(6,*) ' asymmetric update'
        else if (itype.eq.2) then
          write(6,*) ' symmetric update'
        end if
      end if

      if (ntest.ge.10) then
        write(6,*) 'on entry in optc_sbspja:'
        write(6,*)
     &       'nstdim, ngvdim, navec, maxvec, nadd, navec_last: ',
     &        nstdim, ngvdim, navec, maxvec, nadd, navec_last
        write(6,*)
     &       'lugrvf, lu_newstp, lu_pertstp, lust_sbsp, lugv_sbsp: ',
     &        lugrvf, lu_newstp, lu_pertstp, lust_sbsp, lugv_sbsp
      end if

      if (itype.lt.1.or.itype.gt.2) then
        write(6,*) 'SBSP hessian/jacobian: unknown itype = ',itype
        stop 'SBSPJA arguments'
      end if

      if (navec.gt.maxvec) then
        write(6,*)
     &       'Panic 2: Hess./Jac. subspace dimensioning inconsistent!'
        stop 'SBSP dimensions'
      end if

      nskipst = nstdim-navec
      nskipgv = ngvdim-navec
      ndel  = navec_last+nadd-navec
      nold  = navec_last-ndel

      navec_last = navec         ! save for next iteration

      if (nold.lt.0.or.nskipst.lt.0.or.nskipgv.lt.0) then
        write(6,*) 'inconsistent dimensions in optc_sbspja!'
        write(6,*) 'nold , nskipst, nskipgv : ', nold, nskipst, nskipgv
        stop 'SBSPJA dimensions'
      end if

      if (ntest.ge.10) then
        write(6,*) 'vectors to be deleted in U mat (ndel): ',ndel
        write(6,*) 'old vectors on U matrix: (nold) ',nold
        write(6,*) 'skipped vectors on step file:(nskipst)',nskipst
        write(6,*) 'skipped vectors on grad file:(nskipgv)',nskipgv
      end if

* vectors deleted? move B matrix by the corresp. number of col.s and rows
      if (ndel.gt.0) then
        if (isym.eq.1) then ! symmetric U matrix
        ! nold is the future dimension of the already initialized
        ! low-rank matrix
        ! we now move all entries according to
        !   U(i,j) = U(i+ndel,j+ndel), where ij actually refer to a
        !                              upper triangular matrix
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix before deleting:'
            call prtrlt(umat,nold+ndel)
          end if
        
          do ii = 1, nold
            iioff  = ii*(ii-1)/2
            iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
            do jj = 1, ii
              umat(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix after deleting:'
            call prtrlt(umat,nold)
          end if
        else ! non-symmetric (and therefore full) matrix
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix before deleting:'
            call wrtmat2(umat,nold+ndel,maxvec,nold+ndel,maxvec)
          end if
        
          do ii = 1, nold
            iioff  = (ii-1)*maxvec
            iioff2 = (ii+ndel-1)*maxvec + ndel
            do jj = 1, nold
              umat(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix after deleting:'
            call wrtmat2(umat,nold,maxvec,nold,maxvec)
          end if

        end if ! isym

      end if ! ndel.gt.0

      ixd_iter = 0  ! iteration counter for xdamp optimization
      isubcnt = 0   ! counter for optimizer steps
      converged = .false.

      xdamp_min = max(-fdmnxd(ludia,-1,vec1,1,lblk),0d0)
      xlow = xdamp_min
      xhigh = 1d13
      xvlow = 100d0
      xvhigh = -100d0
      
c test: if we damp, start from a somewhat larger value
CCC      if (xdamp.gt.0d0) xdamp = xdamp+0.5d0

      do while(.not.converged)
      
* add new vectors
        ndel = 0
        if (xdamp.eq.0d0) then
c          nnew = nadd
c          nold2 = nold
c for the moment: rebuild everything, everytime
          nnew = nold+nadd
          nold2 = 0
        else
          nnew = nold+nadd
          nold2=0
          call optc_smat(umat,nnew,0,nold2,
     &         lust_sbsp,nskipst,nold+nadd,
     &         ludum,0,0,
     &         vec1,vec2,lincore)
        end if

        call optc_sbspja_umat(itype,
     &       umat,nnew,nold2,ndel,maxvec,
     &       nskipst,nskipgv,
     &       lust_sbsp,lugv_sbsp,ludia,xdamp,
     &       vec1,vec2,lincore)

* for xdamp.ne.0d0, we produce the perturbation step here:
c        if (xdamp.ne.0d0) then
        ! to be sure, rebuild the perturbation step always here:
          call dmtvcd2(vec1,vec2,ludia,lugrvf,lu_pertstp,
     &         -1d0,xdamp,1,1,lblk)
c        end if

* get projection of perturbation step in current subspace
        call optc_sbspja_prjpstp(scrvec2,navec,nskipst,
     &                         lu_pertstp,lust_sbsp,
     &                         vec1,vec2,lincore)

        navec_l = navec

      again = .true.
      do while(again)
* transfer to scratch array
        if (isym.eq.1) then
          do ii = 1, navec_l
            iioff  = ii*(ii-1)/2
            iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
            do jj = 1, ii
              scr(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix passed to dspco:'
            call prtrlt(scr,navec_l)
          end if
        else
          do ii = 1, navec_l
            iioff = (ii-1)*navec_l
            iioff2 = (ii+ndel-1)*maxvec + ndel
            do jj = 1, navec_l
              scr(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix passed to dgeco:'
            call wrtmat2(scr,navec_l,navec_l,navec_l,navec_l)
          end if
        end if

* invert low-rank matrix:
* factorize and get condition
        if (isym.eq.1) then
          call dspco(scr,navec_l,kpiv,cond,scrvec)
          if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &         navec_l,cond
          if (ntest.ge.100) then
            write(6,*) 'Factorized matrix from dspco:'
            call prtrlt(scr,navec_l)
            write(6,*) 'Pivot array:'
            call iwrtma(kpiv,1,navec_l,1,navec_l)
          end if
        else
          call dgeco(scr,navec_l,navec_l,kpiv,cond,scrvec)
          if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &         navec_l,cond
          if (ntest.ge.100) then
            write(6,*) 'Factorized matrix from dgeco:'
            call wrtmat2(scr,navec_l,navec_l,navec_l,navec_l)
            write(6,*) 'Pivot array:'
            call iwrtma(kpiv,1,navec_l,1,navec_l)
          end if
        end if ! isym

c we do not care for the moment        
c        again = cond.lt.1d-14
        again = .false.

        if (.not.again) then
          irhsoff = navec - navec_l + 1
          scrvec(1:navec_l) = scrvec2(irhsoff:irhsoff-1+navec_l)
* if OK, solve U x = rhs
          if (ntest.ge.100) then
            write(6,*) 'RHS: ',scrvec(1:navec_l)
          end if
          if (isym.eq.1) then
            call dspsl(scr,navec_l,kpiv,scrvec)
          else
            job = 0
            call dgesl(scr,navec_l,navec_l,kpiv,scrvec,job)
          end if
          if (ntest.ge.100) then
            write(6,*) 'Result for x:'
            call wrtmat(scrvec,navec_l,1,navec_l,1)
          end if
        end if

        if (.not.again) then
* get low-rank contribution to new vector
c TEST unset thrsh
          thrsh_l = 1d10
          call optc_sbspja_lrstep(itype,thrsh_l,accept,scrvec,navec_l,
     &         nskipst,nskipgv,
     &         lust_sbsp,lugv_sbsp,ludia,xdamp,
     &         lu_pertstp,lu_newstp,
     &         vec1,vec2,lincore)

          xnrm = sqrt(inprdd(vec1,vec1,lu_newstp,lu_newstp,1,lblk))

          xval = xnrm-trrad
          converged = xnrm.le.trrad.and.(xdamp-xdamp_min).le.thrdamp
     &            .or.( xdamp.ne.0d0 .and. (abs(xnrm-trrad).lt.thrdamp))

          if (xval.gt.0d0) then
            xlow = xdamp
            xvlow = xval
          else
            xhigh = xdamp
            xvhigh = xval
          end if
          write(6,'(">>",4x,2i4,3(2x,e20.5))')
     &         ixd_iter+1,isubcnt,xdamp,xnrm,xnrm-trrad

* convergence control for xdamp:
*  the step-length function should for any xdamp larger than the
*  negative of the lowest eigenvalue of the approximate Jacobian
*  be a monotonically decreasing konvex function. If our Quasi-Newton
*  algorithm encounters problems, that can only mean that we are still
*  not beyond the lowest eigenvalue.
          if (.not.converged) then
            if (xdamp.lt.0d0) then ! that should never ever happen!!
              write(6,*) 'What did you do? xdamp = ',xdamp
              stop 'optc_sbspja: xdamn'
            end if
            ixd_iter = ixd_iter+1
            if (ixd_iter.ge.mxd_iter) then
              write(6,*) 'unsolvable problems finding xdamp'
              stop 'optc_sbspja: xdamp'
            end if
            xd(ixd_iter) = xdamp
            xv(ixd_iter) = xnrm-trrad
            if (isubcnt.eq.0) then
              ! restart at a very new point
              ! make a small increment to get the gradient              
              if (xv(ixd_iter).gt.0d0) then
                dxdamp = + xinc
              else
                dxdamp = - xinc
              end if
              if (ntest.ge.150)
     &             write(6,*)
     &             'xdamp> initialized num. gradient with xinc = ',
     &                        xinc
              ! next time we see us in step two
              isubcnt = 2
            else if (isubcnt.eq.1) then
              ! this is step one of a numerical Newton optimization
              ! make a small increment to get the gradient              
              if (xv(ixd_iter).gt.0d0) then
                dxdamp = + xinc
              else
                dxdamp = - xinc
              end if
              if (ntest.ge.150)
     &             write(6,*)
     &             'xdamp> initialized num. gradient with xinc = ',
     &                        xinc
              ! next time we see us in step two
              isubcnt = 2
            else if (isubcnt.eq.2) then
              ! get gradient from finite difference
              if (ntest.ge.150) then
                write(6,*) 'xdamp> linear model '
                write(6,*) ' points used:'
                do ii = 0, 1
                  write(6,'(3x,i2,2(2x,e25.8))')
     &                 ii, xd(ixd_iter-ii), xv(ixd_iter-ii)
                end do
              end if
              xdiff = xd(ixd_iter)-xd(ixd_iter-1)
              xgr = (xv(ixd_iter)-xv(ixd_iter-1))/xdiff
              if (ntest.ge.150)
     &             write(6,*) 'xdamp> current num. gradient: ',xgr    
              ! gradient has to be negative
              if (xgr.lt.0d0) then
                if (ntest.ge.150)
     &             write(6,*) 'xdamp> accepting step'  
                ! make a Newton step
c                dxdamp = - xv(ixd_iter-1)/xgr
                dxdamp = - xv(ixd_iter)/xgr
                isubcnt = 1
                if (ntest.eq.150)
     &               write(6,*) 'xdamp> step = ',dxdamp
c switch to three point formula disabled
c                if (dxdamp.lt.0.5d0) isubcnt = 3
                if (dxdamp.lt.0.5d0) isubcnt = 2
              else
                if (ntest.ge.150)
     &             write(6,*) 'xdamp> search at a new place'  
                ! retry with a larger xdamp
                dxdamp = + xfailfac*xinc
                ixd_iter = ixd_iter-1
                isubcnt = 0
              end if
            else if (isubcnt.eq.3) then
              if (ntest.ge.150) then
                write(6,*) 'xdamp> quadratic model'  
                write(6,*) ' points used:'
                do ii = 0, 2
                  write(6,'(3x,i2,2(2x,e25.8))')
     &                 ii, xd(ixd_iter-ii), xv(ixd_iter-ii)
                end do
              end if
              ! for convenience:
              f0 = xv(ixd_iter)
              f1 = xv(ixd_iter-1)
              f2 = xv(ixd_iter-2)
              a1 = xd(ixd_iter-1)-xd(ixd_iter)
              a2 = xd(ixd_iter-2)-xd(ixd_iter)
              ! the points should correspond to a monotonic
              ! decreasing function
              if (.not.(a1*(f0-f1).gt.0d0.and.a2*(f0-f2).gt.0d0)) then
                if (ntest.ge.150)
     &             write(6,*) 'xdamp> search at a new place'  
                ! search somewhere else
                dxdamp = + xfailfac*xinc
                xhigh = 1d13
                xvhigh = -100d0
                isubcnt = 0
              else
                ! get gradient and hessian from last three values
                ! gradient
                xgr = ((f1-f0)-(a1*a1)/(a2*a2)*(f2-f0))/(a1-(a1*a1)/a2)
                if (ntest.ge.150)
     &             write(6,*) 'xdamp> num. gradient ',xgr  
                ! hessian
                xhs = 2d0*((f1-f0)-a1/a2*(f2-f0))/(a1*a1-a1*a2)
                if (ntest.ge.150)
     &               write(6,*) 'xdamp> num. hessian ',xhs  
                ! discriminant
                xdis = (xgr*xgr)/(xhs*xhs) - 2d0*f0/xhs
                if (ntest.ge.150)
     &               write(6,*) 'xdamp> discriminant ',xdis  
                ! hessian positive?
                if (xhs.gt.0d0.and.xdis.gt.0d0) then
                  if (ntest.ge.150)
     &                 write(6,*) 'xdamp> accepting solution ' 
                  ! take the lower solution
                  dxdamp = -xgr/xhs - sqrt(xdis)
c                  dxdamp = -xgr/xhs + sqrt(xdis)
                  ! we go on with the 3-point search
                  isubcnt = 3
                else
                  ! gradient negative (and to be trusted)
                  if (xgr.lt.0d0.and.abs(a2).lt.xinc) then
                    if (ntest.ge.150)
     &                   write(6,*) 'xdamp> follow gradient step' 
                    ! take only the gradient step
                    dxdamp = - f0/xgr
                    ! we go on with the 3-point search
                    isubcnt = 3
                  else
                    ! xgr.gt.0 does not really mean something:
                    ! we try to get a new gradient here
                    if (ntest.ge.150)
     &                   write(6,*)
     &                     'xdamp> prepare new gradient evaluation' 
                    if (f0.gt.0d0) then
                      dxdamp = + xinc
                    else
                      dxdamp = - xinc
                    end if
                    isubcnt = 2
                  end if
                end if ! hessian > 0?
              end if ! monotonic function?
            else
              write(6,*) 'unknown isubcnt = ',isubcnt
              stop 'optc_sbspja'
            end if

c            if (isubcnt.ne.0) then
c              ! we use a simple step restriction
c              xmxdxdmp2 = max(xmxdxdamp,xdamp-1d0/xdamp)
c              dxdamp = sign(min(abs(dxdamp),xmxdxdmp2),dxdamp)
c            end if
            ! unless isubcnt was reset to 0
            if (isubcnt.ne.0) then
             if (xdamp+dxdamp.ge.xhigh) then
              if (xval.eq.xvhigh) stop 'strange: xdamp (1)'
              xdamp = xdamp-xval*(xdamp-xhigh)/(xval-xvhigh)
             else if (xdamp+dxdamp.le.xlow) then
              if (xlow.eq.xdamp_min) xvlow=-xval
              if (xval.eq.xvlow) stop 'strange: xdamp (2)'
              xdamp = xdamp-xval*(xdamp-xlow)/(xval-xvlow)
             else
              xdamp = xdamp+dxdamp
             end if
            else
              xdamp = xdamp+dxdamp
            end if

c            if (isubcnt.eq.0) then
c              xdamp = xdamp + dxdamp
c            else
c              ! we use a simple step restriction
c              xmxdxdmp2 = max(xmxdxdamp,xdamp-1d0/xdamp)
c              dxdamp = sign(min(abs(dxdamp),xmxdxdmp2),dxdamp)
c              xdamp = max(0d0,xdamp + dxdamp)
c            end if

            if (ntest.ge.150)
     &           write(6,*)'xdamp> next xdamp = ',xdamp

c            if (xdamp.eq.0d0) xdamp = 10d0
c            xdamp = xdamp - 0.009d0
c            if (xdamp.lt.-4d0) stop 'test'
          end if ! .not.converged
        

c* and analyze it
c          if (.not.accept.and.navec_l.gt.1) then
c            again = .true.
c          else if (.not.accept) then
c            navec_l = 0
c          end if
        end if ! again

        if (again) then
          nskipst = nskipst+1
          nskipgv = nskipgv+1
          ndel = ndel+1
          navec_l = navec_l-1
        end if

      end do ! again

      end do ! xdamp - optimization

* return the actual last subspace dimension
      navec = navec_l

* make an energy prediction (for variational cases)
c
c      RHS : scrvec2(irhsoff:irhsoff-1+navec_l)
c      c   : scrvec
c
c      call matvcb(work(khred),work(kscr1),work(kscr2),
c     &            imicdim,imicdim,0)
c      de_pred = inprod(work(kgred),work(kscr1),imicdim) +
c     &      0.5d0*inprod(work(kscr2),work(kscr1),imicdim)
c

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja_new_iter(itype,thrsh,nstdim,ngvdim,
     &           navec,maxvec,
     &           nadd,navec_last,
     &           lugrvf,ludia,trrad,lu_pertstp,lu_newstp,
     &           xdamp,xdamp_last,
     &           lust_sbsp,lugv_sbsp,
     &           umat,scr,vec1,vec2,iprint)
*----------------------------------------------------------------------*
*
* subspace jacobian: A = A_0(1-P) + A_ex P
*
*  if xdamp is 0d0, we try the direct inversion procedure; if that
*  step is too large, we go over to the iterative solver to find the
*  optimal damping (where we end up directly, if xdamp was != 0d0) 
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00
      
* input/output
      integer, intent(in) ::
     &     itype,                    ! type (s.above)
     &     nstdim,ngvdim,            ! # records on sbsp files
     &     maxvec,                   ! max dimension of sbsp
     &     nadd,                     ! number of new vectors (usually 1)
     &     lugrvf,ludia,             ! current gradient, diagonal Hess/Jac.
     &     lu_pertstp,               ! step from diagonal Hess/Jac.
     &     lu_newstp,                ! new step (output)
     &     lust_sbsp,lugv_sbsp,      ! sbsp files
     &     iprint                    ! print level
      integer, intent(inout) ::
     &     navec,                    ! current dimension of sbsp
     &     navec_last                ! previous sbsp dimension
      real(8), intent(in) ::
     &     thrsh,                    ! thresh for accepting low-rank correction
     &     trrad,                    ! trust radius for step
     &     xdamp_last                ! previous damping
      real(8), intent(out) ::
     &     xdamp                     ! current damping
      real(8), intent(inout) ::
     &     umat(*),                  ! low rank Hessian/Jacobian
     &     scr(*), vec1(*), vec2(*)  ! scratch

* local O(N) scratch
      integer ::
     &     kpiv(navec)
      real(8) ::
     &     scrvec(navec), scrvec2(navec)

* local
      logical ::
     &     lincore, again, accept, converged, lin_dep
      integer ::
     &     lblk, iprintl, navec_l, isym, job,
     &     nskipst, nskipgv, nnew, nold, nold2, ndel,
     &     ii, jj, iioff, iioff2, irhsoff,
     &     lu_sbsp_un, lu_curerr,
     &     luvec, luAvec, luvec_sbsp, luAvec_sbsp, lures,
     &     nadd_l, nold_l, nskip_l, ludum, iter, maxiter,
     &     kend, ksmat, ksred, kared, kevec, kscr,
     &     itype_evp, itask(2), ntask, ndum
      real(8) ::
     &     cond, xcorsum, xlam, xresnrm, fac, thrres, xnrm
      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd

      iprintl = max(iprint,ntest)
      lblk = -1
      lincore = .false.
      if (itype.eq.1) isym = 0  ! non-symmetric jacobian
      if (itype.eq.2) isym = 1  ! symmetric jacobian

      if (iprintl.ge.5) then
        write(6,*) 'Subspace-Hessian/Jacobian'
        write(6,*) '========================='
        if (itype.eq.1) then
          write(6,*) ' asymmetric update'
        else if (itype.eq.2) then
          write(6,*) ' symmetric update'
        end if
      end if

      if (ntest.ge.10) then
        write(6,*) 'on entry in optc_sbspja:'
        write(6,*)
     &       'nstdim, ngvdim, navec, maxvec, nadd, navec_last: ',
     &        nstdim, ngvdim, navec, maxvec, nadd, navec_last
        write(6,*)
     &       'lugrvf, lu_newstp, lu_pertstp, lust_sbsp, lugv_sbsp: ',
     &        lugrvf, lu_newstp, lu_pertstp, lust_sbsp, lugv_sbsp
      end if

      if (itype.lt.1.or.itype.gt.2) then
        write(6,*) 'SBSP hessian/jacobian: unknown itype = ',itype
        stop 'SBSPJA arguments'
      end if

      if (navec.gt.maxvec) then
        write(6,*)
     &       'Panic 2: Hess./Jac. subspace dimensioning inconsistent!'
        stop 'SBSP dimensions'
      end if

      nskipst = nstdim-navec
      nskipgv = ngvdim-navec
      ndel  = navec_last+nadd-navec
      nold  = navec_last-ndel

      navec_last = navec         ! save for next iteration

      if (nold.lt.0.or.nskipst.lt.0.or.nskipgv.lt.0) then
        write(6,*) 'inconsistent dimensions in optc_sbspja!'
        write(6,*) 'nold , nskipst, nskipgv : ', nold, nskipst, nskipgv
        stop 'SBSPJA dimensions'
      end if

      if (ntest.ge.10) then
        write(6,*) 'vectors to be deleted in U mat (ndel): ',ndel
        write(6,*) 'old vectors on U matrix: (nold) ',nold
        write(6,*) 'skipped vectors on step file:(nskipst)',nskipst
        write(6,*) 'skipped vectors on grad file:(nskipgv)',nskipgv
      end if

* vectors deleted? move B matrix by the corresp. number of col.s and rows
      if (ndel.gt.0) then
        if (isym.eq.1) then ! symmetric U matrix
        ! nold is the future dimension of the already initialized
        ! low-rank matrix
        ! we now move all entries according to
        !   U(i,j) = U(i+ndel,j+ndel), where ij actually refer to a
        !                              upper triangular matrix
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix before deleting:'
            call prtrlt(umat,nold+ndel)
          end if
        
          do ii = 1, nold
            iioff  = ii*(ii-1)/2
            iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
            do jj = 1, ii
              umat(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix after deleting:'
            call prtrlt(umat,nold)
          end if
        else ! non-symmetric (and therefore full) matrix
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix before deleting:'
            call wrtmat2(umat,nold+ndel,maxvec,nold+ndel,maxvec)
          end if
        
          do ii = 1, nold
            iioff  = (ii-1)*maxvec
            iioff2 = (ii+ndel-1)*maxvec + ndel
            do jj = 1, nold
              umat(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix after deleting:'
            call wrtmat2(umat,nold,maxvec,nold,maxvec)
          end if

        end if ! isym

      end if ! ndel.gt.0
      
* add new vectors
      ndel = 0
      nnew = nadd
      nold2 = nold
      call optc_sbspja_umat(itype,
     &     umat,nnew,nold2,ndel,maxvec,
     &     nskipst,nskipgv,
     &     lust_sbsp,lugv_sbsp,ludia,xdamp,
     &     vec1,vec2,lincore)

* get projection of perturbation step in current subspace
      call optc_sbspja_prjpstp(scrvec2,navec,nskipst,
     &                         lu_pertstp,lust_sbsp,
     &                         vec1,vec2,lincore)

      navec_l = navec

      if (xdamp.eq.0d0) then 

* direct-inversion section:
* =========================
      again = .true.
      do while(again)
* transfer to scratch array
        if (isym.eq.1) then
          do ii = 1, navec_l
            iioff  = ii*(ii-1)/2
            iioff2 = (ii+ndel)*(ii+ndel-1)/2 + ndel
            do jj = 1, ii
              scr(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix passed to dspco:'
            call prtrlt(scr,navec_l)
          end if
        else
          do ii = 1, navec_l
            iioff = (ii-1)*navec_l
            iioff2 = (ii+ndel-1)*maxvec + ndel
            do jj = 1, navec_l
              scr(iioff+jj) = umat(iioff2+jj)
            end do
          end do
          if (ntest.ge.100) then
            write(6,*) 'Low-rank matrix passed to dgeco:'
            call wrtmat2(scr,navec_l,navec_l,navec_l,navec_l)
          end if
        end if

* invert low-rank matrix:
* factorize and get condition
        if (isym.eq.1) then
          call dspco(scr,navec_l,kpiv,cond,scrvec)
          if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &         navec_l,cond
          if (ntest.ge.100) then
            write(6,*) 'Factorized matrix from dspco:'
            call prtrlt(scr,navec_l)
            write(6,*) 'Pivot array:'
            call iwrtma(kpiv,1,navec_l,1,navec_l)
          end if
        else
          call dgeco(scr,navec_l,navec_l,kpiv,cond,scrvec)
          if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &         navec_l,cond
          if (ntest.ge.100) then
            write(6,*) 'Factorized matrix from dgeco:'
            call wrtmat2(scr,navec_l,navec_l,navec_l,navec_l)
            write(6,*) 'Pivot array:'
            call iwrtma(kpiv,1,navec_l,1,navec_l)
          end if
        end if ! isym
        
        again = cond.lt.1d-14

        if (.not.again) then
          irhsoff = navec - navec_l + 1
          scrvec(1:navec_l) = scrvec2(irhsoff:irhsoff-1+navec_l)
* if OK, solve U x = rhs
          if (ntest.ge.100) then
            write(6,*) 'RHS: ',scrvec(1:navec_l)
          end if
          if (isym.eq.1) then
            call dspsl(scr,navec_l,kpiv,scrvec)
          else
            job = 0
            call dgesl(scr,navec_l,navec_l,kpiv,scrvec,job)
          end if
          if (ntest.ge.100) then
            write(6,*) 'Result for x:'
            call wrtmat(scrvec,navec_l,1,navec_l,1)
          end if
        end if

        if (.not.again) then
* get low-rank contribution to new vector
          call optc_sbspja_lrstep(itype,thrsh,accept,scrvec,navec_l,
     &         nskipst,nskipgv,
     &         lust_sbsp,lugv_sbsp,ludia,xdamp,
     &         lu_pertstp,lu_newstp,
     &         vec1,vec2,lincore)

          xnrm = sqrt(inprdd(vec1,vec1,lu_newstp,lu_newstp,1,lblk))
          accept = xnrm.le.trrad

c* and analyze it
c          if (.not.accept.and.navec_l.gt.1) then
c            again = .true.
c          else if (.not.accept) then
c            navec_l = 0
c          end if
        end if

        if (again) then
          nskipst = nskipst+1
          nskipgv = nskipgv+1
          ndel = ndel+1
          navec_l = navec_l-1
        end if

      end do

      ! abuse xdamp as flag
      if (.not.accept) xdamp = 100d0

* end of direct-inversion section
      end if

* here starts the iterative inversion section with damping:
* =========================================================
      if (xdamp.ne.0d0) then

        kend = 1
        
        ksmat = kend
        kend = kend + navec_l*(navec_l+2)/2

        kscr = kend
        kend = kend + 27*(navec_l+2)**2

        ksred = kend
        kend = kend + 27*(navec_l+2)**2

        kared = kend
        kend = kend + 27*(navec_l+2)**2

        kevec = kend
        kend = kend + 27*(navec_l+2)**2

        luvec  = iopen_nus('SSJA_VEC')
        luAvec = iopen_nus('SSJAAVEC')
        luvec_sbsp  = iopen_nus('SSJA_V_SP')
        luAvec_sbsp = iopen_nus('SSJAAV_SP')
        lures  = iopen_nus('SSJARVEC')

        ! set up the S matrix
        nadd_l = navec_l
        nold_l = 0
        nskip_l = 0
        call optc_smat(scr(ksmat),navec_l,nold_l,nskip_l,
     &       lust_sbsp,nskipst,nstdim,
     &       ludum,0,0,
     &       vec1,vec2,lincore)

        ! decompose it
        call dspco(scr(ksmat),navec_l,kpiv,cond,scrvec)
        if (ntest.ge.10) write(6,*)'navec_l,condition number: ',
     &       navec_l,cond
        if (ntest.ge.100) then
          write(6,*) 'Factorized matrix from dspco:'
          call prtrlt(scr(ksmat),navec_l)
          write(6,*) 'Pivot array:'
          call iwrtma(kpiv,1,navec_l,1,navec_l)
        end if
        ! delete vectors in subspace if S becomes singular
        ! to be done later ...

        ! microiterations:
        converged = .false.
        iter = 0
        maxiter = 3*(navec_l+1)
        thrres = 1d-6

        call sclvcd(lu_pertstp,luvec,-1d0,vec1,1,lblk)

        do while(.not.converged)
          iter = iter+1
          ! orthonormalize trial vector
          call optc_orthvec(scr(kscr),scr(ksred),
     &         iter-1,1,lin_dep,
     &         luvec_sbsp,luvec,
     &         vec1,vec2)

          ! get Mv product
          call optc_sbspja_mvp(itype,
     &         luvec,luAvec,
     &         scr(ksmat),kpiv,navec_l,
     &         lust_sbsp,nskipst,lugv_sbsp,nskipgv,
     &         ludia,
     &         vec1,vec2,lincore)

          ! solve reduced EVP
          itype_evp = 0 ! trust-radius method
          !itype_evp = 1 ! augmented hessian method
          call optc_sbspja_redevp(itype_evp,
     &         scr(kared),scr(ksred),scr(kscr),scr(kevec),trrad,
     &         xlam,lures,xresnrm,
     &         iter,
     &         lugrvf,
     &         luvec,luAvec,luvec_sbsp,luAvec_sbsp,
     &         vec1,vec2,lincore)
          
          ! test convergence
          if (abs(xresnrm).lt.thrres) converged = .true.

          write(6,'(">>",i4,2(2x,e12.2),2x,l)')
     &         iter,xlam,xresnrm,converged

          if (iter.eq.maxiter) exit

          if (.not.converged) then
            ! put luvec on subspace-file
            fac = 1d0
            itask(1) = 1
            itask(2) = 1
            ntask = 2
            call optc_sbspman(luvec_sbsp,luvec,fac,ludum,iter-1,maxiter,
     &           itask,ntask,0,ndum,vec1,vec2)

            ! put luAvec on subspace file
            call optc_sbspman(luAvec_sbsp,luAvec,fac,
     &                                             ludum,iter-1,maxiter,
     &           itask,ntask,0,ndum,vec1,vec2)

            ! precondition new direction:
            call dmtvcd2(vec1,vec2,ludia,lures,luvec,1d0,0d0,1,1,lblk)

          end if

        end do

        if (.not.converged) then
          write(6,*) 'WARNING: Microiterations not converged'
          write(6,*) '         final residual : ',xresnrm
        end if

        ! assemble new step from first vector in scr(kevec) and
        ! the subspace vectors in luvec_sbsp and the last
        ! trial vector (which we didn't put on the subspace file)
        if (iter.gt.1) then
          ! luAvec and luAvec_sbsp used as scratch
          call mvcsmd(luvec_sbsp,scr(kevec),luAvec,luAvec_sbsp,
     &              vec1,vec2,iter-1,1,lblk)
c TEST
c          call vecsmd(vec1,vec2,-1d0,-scr(kevec-1+iter),
          call vecsmd(vec1,vec2,1d0,scr(kevec-1+iter),
     &         luAvec,luvec,lu_newstp,1,lblk)
        else
c TEST
c          call sclvcd(luvec,lu_newstp,-scr(kevec),vec1,1,lblk)
          call sclvcd(luvec,lu_newstp,scr(kevec),vec1,1,lblk)
        end if

        call relunit(lures,'delete')
        call relunit(luvec,'delete')
        call relunit(luAvec,'delete')
        call relunit(luvec_sbsp,'delete')
        call relunit(luAvec_sbsp,'delete')

      end if

* return the actual last subspace dimension
      navec = navec_l

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_orthvec(scr1,scr2,ndim1,ndim2,lin_dep,
     &                        luvec1,luvec2,
     &                        vec1,vec2)
*----------------------------------------------------------------------*
*
*     luvec1 contains a set of ndim1 orthonormal vectors,
*     luvec2 a set of ndim2 non-orthonormal vectors. On output luvec2
*     is orthonormalized against luvec1 and among itself.
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ndim1, ndim2,
     &     luvec1, luvec2
      logical, intent(out) ::
     &     lin_dep(ndim2)
      real(8), intent(inout) ::
     &     scr1(ndim1+ndim2,ndim1+ndim2),
     &     scr2(ndim1+ndim2,ndim1+ndim2),
     &     vec1(*), vec2(*)

      integer ::
     &     lblk, ii, jj, ndim,
     &     lusc1, lusc2, lusc3, lusc4
      real(8) ::
     &     xs, xnorm

      real(8) ::
     &     scrvec(ndim1+ndim2)

      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd, inprod

      lblk = -1

      lusc1 = iopen_nus('OPTC_ORTHSC1')
      lusc2 = iopen_nus('OPTC_ORTHSC2')
      lusc3 = iopen_nus('OPTC_ORTHSC3')
      lusc4 = iopen_nus('OPTC_ORTHSC4')

      ! set up metric on s
      ndim = ndim1+ndim2
      scr1(1:ndim,1:ndim) = 0d0
      do ii = 1, ndim1
        scr1(ii,ii) = 1d0
      end do

      call rewino(luvec2)
      do jj = ndim1+1, ndim
        call rewino(lusc1)
        call copvcd(luvec2,lusc1,vec1,0,lblk)
        call rewino(luvec1)
        do ii = 1, ndim1
          call rewino(lusc1)
          xs = inprdd(vec1,vec2,lusc1,luvec1,0,lblk)
          scr1(ii,jj) = xs
          scr1(jj,ii) = xs
        end do
      end do

      ! get first element
      scr1(ndim1+1,ndim1+1) = inprdd(vec1,vec1,luvec2,luvec2,1,lblk)
      do ii = ndim1+2, ndim
        ! we are now at element ii, copy to scratch
        call copvcd(luvec2,lusc1,vec1,0,lblk)
        ! and rewind
        call rewino(luvec2)
        do jj = ndim1+2, ii
          call rewino(lusc1)
          xs = inprdd(vec1,vec2,luvec2,lusc1,0,lblk)
          scr1(ii,jj) = xs
          scr1(jj,ii) = xs
        end do
      end do

      if (ntest.ge.100) then
        write(6,*) 'overlap matrix:'
        call wrtmat2(scr1,ndim,ndim,ndim,ndim)
      end if

      ! call modified Gram-Schmidt orthonormalization
      call mgs3(scr2,scr1,ndim,scrvec)

      if (ntest.ge.100) then
        write(6,*) 'orthonormalization matrix:'
        call wrtmat2(scr2,ndim,ndim,ndim,ndim)
      end if

      do ii = ndim1+1, ndim
        xnorm = inprod(scr2(1,ii),scr2(1,ii),ndim)
        if (xnorm.lt.epsilon(1d0)) then
          if (ntest.ge.10) write(6,*)
     &         ' linear dependency detected for vector ',ii
          lin_dep(ii-ndim1) = .true.
          ! leave this vector unmodified
          scr2(ii,ii) = 1d0
        else
          lin_dep(ii-ndim1) = .false.
        end if
      end do

      ! update luvec2
      call rewino(luvec2)
      call rewino(lusc2)
      do ii = 1, ndim2
        call copvcd(luvec2,lusc2,vec1,0,lblk)
      end do

      call rewino(luvec2)
      do ii = ndim1+1, ndim
        ! assemble new vector ii for luvec2
        ! contributions from luvec1
        if (ndim1.gt.0)
     &    call mvcsmd(luvec1,scr2(1,ii),lusc1,lusc3,
     &         vec1,vec2,ndim1,1,lblk)

        ! contributions from luvec2 (now on lusc2)
        if (ndim2.gt.0)
     &    call mvcsmd(lusc2,scr2(ndim1+1,ii),lusc4,lusc3,
     &         vec1,vec2,ndim2,1,lblk)
        
        if (ndim1.gt.0.and.ndim2.gt.0) then
          call rewino(lusc1)
          call rewino(lusc4)
          call vecsmd(vec1,vec2,1d0,1d0,lusc1,lusc4,luvec2,0,lblk)
        else if (ndim1.gt.0) then
          call rewino(lusc1)
          call copvcd(lusc1,luvec2,vec1,0,lblk)
        else if (ndim2.gt.0) then
          call rewino(lusc4)
          call copvcd(lusc4,luvec2,vec1,0,lblk)
        end if

      end do

      call relunit(lusc1,'delete')
      call relunit(lusc2,'delete')
      call relunit(lusc3,'delete')
      call relunit(lusc4,'delete')

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja_redevp(itype,
     &     amat,smat,scr,xvec,xs,
     &     xlam,lures,xresnrm,
     &     ndim,
     &     lugvec,
     &     luvec,luAvec,luvec_sbsp,luAvec_sbsp,
     &     vec1,vec2,lincore)
*----------------------------------------------------------------------*
*
*     solve
*
*            A_aug   |x_aug>  =
*
*      /                            \ /    \            /    \
*      |  <1|A|1> .. <1|A|n> <1|g>  | | x1 |            | x1 |
*      |     .          .      .    | | .  |            | .  |
*      |     .          .      .    | | .  | = lambda S | .  |
*      |  <n|A|1> .. <n|A|n> <n|g>  | | xn |            | xn |
*      |   <g|1>  ..  <g|n>    0    | | 1  |            | 1  |
*      \                            / \    /            \    /
*
*     (i.e. |x_aug> will be normalized such that its last coeff. is +1)
*
*     where for itype==0 (trust radius method)
*
*          /                 \
*          | 1 0 . 0    0    |
*          | 0 .        .    |
*      S = | .   . 0    .    |
*          | 0 . 0 1    0    |
*          | 0 ... 0   -s**2 |
*          \                 /
*
*     or for itype==1 (rational function method)
*
*          /               \
*          | s 0 . 0    0  |
*          | 0 .        .  |
*      S = | .   . 0    .  |
*          | 0 . 0 s    0  |
*          | 0 ... 0    1  |
*          \               /
*
*                    on input:                  on output:
*     amat       A_aug in (n-1)-dim subspc. A_aug in n-dim subspace
*     smat       scratch holding --->       S in n-dim subspace
*     xvec                                  eigenvectors of reduced EVP
*     xs         scalar defining the metric
*     xlam                                  lowest eigenvalue
*     xresnrm                               residual norm in full space
*     lures                                 file: residual vector
*     ndim       current dimension of 
*                subspace
*     lugvec     file: g in full space
*     luvec      file: last trial vector |n>
*     luAvec     file: last A|n>
*     luvec_sbsp  file: (n-1)-dim subspace {|n>} 
*     luAvec_sbsp file: (n-1)-dim subspace {A|n>}
*     vec1,vec2  scratch (size: at least largest block on vector file)
*     lincore    .true. if vec1/vec2 can hold one complete vector
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00
      
* interface:
      logical, intent(in) ::
     &     lincore
      integer, intent(in) ::
     &     itype, ndim,
     &     lures,
     &     lugvec,luvec,luAvec,luvec_sbsp,luAvec_sbsp
      real(8), intent(in) ::
     &     xs
      real(8), intent(out) ::
     &     xlam, xresnrm
      real(8), intent(inout) ::
     &     amat(*), smat(ndim+1,*), scr(*),
     &     xvec(ndim+1,*), vec1(*), vec2(*)

* local
      integer ::
     &     lblk, ndimold, ndimp1, ndimp1old, matz, ierr,
     &     iad, iadold, iminev,
     &     ii, jj, iiold, jjold,
     &     lusc1, lusc2

      real(8) ::
     &     xx, xminev

* local order(N)-scratch, allocated on the stack:
      real(8) ::
     &     evr(ndim+1),evi(ndim+1),evd(ndim+1)

* external
      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd, inprod

c TEST
      real(8) xtest
          
      lblk = -1
      lusc1 = iopen_nus('OPTCEVPSC1')
      lusc2 = iopen_nus('OPTCEVPSC2')

      ndimp1 = ndim+1
      ! set up projected metric matrix
      smat(1:ndim+1,1:ndim+1) = 0d0
      if (itype.eq.0) then
        xx = 1d0
      else
        xx = 1d0/xs
      end if
      do ii = 1, ndim
        smat(ii,ii) = xx
      end do
      if (itype.eq.0) then
        smat(ndim+1,ndim+1) = 1d0!-xs*xs
      else
        smat(ndim+1,ndim+1) = 1d0
      end if

      if (ntest.ge.100) then
        write(6,*) 'The S matrix: '
        call wrtmat(smat,ndimp1,ndimp1,ndimp1,ndimp1)
      end if

      ! update A matrix
      ! the ndim+1,ndim+1 element is zero
      amat((ndimp1-1)*ndimp1+ndimp1) = 0d0
      ! first rearrange from (ndim+1)-1 to (ndim+1)
      !  a) the g column 
      ndimold = ndim-1
      ndimp1old = ndimp1-1
      if (ntest.ge.100.and.ndimp1old.gt.1) then
        write(6,*) 'The previous augmented A matrix: '
        call wrtmat(amat,ndimp1old,ndimp1old,ndimp1old,ndimp1old)
      end if
      jjold = ndimp1old
      jj    = ndimp1
      do iiold = ndimold, 1, -1
        iadold = ((jjold-1)*ndimp1old)+iiold
        iad    = ((jj-1)*ndimp1)+iiold
        amat(iad) = amat(iadold)
      end do

      do jjold = ndimold, 1, -1
        ! b) one element of the g row
        iadold = ((jjold-1)*ndimp1old)+ndimp1old
        iad    = ((jjold-1)*ndimp1old)+ndimp1
        amat(iad) = amat(iadold)

        ! c) the elements of A
        do iiold = ndimold, 1, -1
          iadold = ((jjold-1)*ndimp1old)+iiold
          iad    = ((jjold-1)*ndimp1)+iiold
          amat(iad) = amat(iadold)
        end do
          
      end do

      if (ntest.ge.100) then
        write(6,*) 'The previous augm. A matrix after rearranging: '
        call wrtmat(amat,ndimp1,ndimp1,ndimp1,ndimp1)
      end if

      ! update the new row in A_aug
      call rewino(luvec_sbsp)
      do ii = 1, ndim-1
        call rewino(luAvec)
        amat(((ndim-1)*ndimp1)+ii)
     &       = inprdd(vec1,vec2,luvec_sbsp,luAvec,0,lblk)
      end do

      ! update the new row in A_aug
      call rewino(luAvec_sbsp)
      do jj = 1, ndim-1
        call rewino(luvec)
        amat(((jj-1)*ndimp1)+ndim)
     &       = inprdd(vec1,vec2,luvec,luAvec_sbsp,0,lblk)
      end do

      ! the (ndim,ndim)-element:
      amat((ndim-1)*ndimp1+ndim) = inprdd(vec1,vec2,luvec,luAvec,1,lblk)

      ! update with <g|x_n>
      xx = inprdd(vec1,vec2,lugvec,luvec,1,lblk)
c TEST
      amat((ndimp1-1)*ndimp1+ndim) = 0d0!xx
      amat((ndim-1)*ndimp1+ndimp1) = 0d0!xx

      if (ntest.ge.100) then
        write(6,*) 'The new augmented A matrix: '
        call wrtmat(amat,ndimp1,ndimp1,ndimp1,ndimp1)
      end if

      ! solve reduced EVP
      matz = 1
      ! transfer amat to scratch array (as rgg destroys it)
      scr(1:ndimp1**2) = amat(1:ndimp1**2)
      ! call the general EVP-solver from eispack:
c TEST
      call rgg(ndimp1,ndim,scr,smat,evr,evi,evd,
c      call rgg(ndimp1,ndimp1,scr,smat,evr,evi,evd,
     &         matz,xvec,ierr)
      if (ierr.ne.0) then
        write(6,*) 'error code from rgg: ',ierr
        stop 'optc_sbspja_redevp'
      end if

      xminev = 1d100
      iminev = 0
c TEST
      do ii = 1, ndim
c      do ii = 1, ndimp1
        evr(ii) = evr(ii)/evd(ii)
        if (evr(ii).lt.xminev) then
          xminev = evr(ii)
          iminev = ii
        end if
        evi(ii) = evi(ii)/evd(ii)
      end do 
      if (iminev.eq.0.or.iminev.gt.ndimp1) then
        write(6,*) 'internal error'
        stop 'optc_sbspja_redevp'
      end if

      if (ntest.ge.100) then
        write(6,*) 'The eigenvalues from RGG:'
        do ii = 1, ndimp1
          write(6,'(x,i4,2f12.6)') ii,evr(ii),evi(ii)
        end do
        write(6,*) 'The eigenvector with lowest real component: ',
     &       evr(iminev)
        call wrtmat(xvec(1,iminev),ndimp1,1,ndimp1,1)
        if (evi(1).ne.0d0) then
          write(6,*) 'imaginary component:'
          call wrtmat(xvec(2,iminev),ndimp1,1,ndimp1,1)
        end if
      end if

      xlam = evr(iminev)
      if (evi(iminev).ne.0d0) then
        write(6,*) 'Huah! Imaginary lowest eigenvalue encountered!'
        write(6,*) evr(iminev),' +/- i * ',evi(iminev)
        write(6,*) 'Don''t know what to do and give up ...'
        stop 'optc_sbspja_redevp'
      end if
      ! renormalize
c TEST
      xx = inprod(xvec(1,iminev),xvec(1,iminev),ndim)
c      xx = xvec(ndimp1,iminev)
      do ii = 1, ndimp1
        xvec(ii,iminev) = xvec(ii,iminev)/xx
      end do

      if (ntest.ge.100) then
        write(6,*) 'The lowest eigenvector (renormalized):'
        call wrtmat(xvec(1,iminev),ndimp1,1,ndimp1,1)
      end if

      ! get residual in full space
      if (itype.eq.0) then
        xx = 1d0
      else
        xx = 1d0/xs
      end if

      if (ndim.gt.1) then
        ! (A - lambda S) x ...
        call mvcsmd(luAvec_sbsp,xvec(1,iminev),lusc1,lusc2,
     &              vec1,vec2,ndim-1,1,lblk)
        call vecsmd(vec1,vec2,1d0,xvec(ndim,iminev),
     &              lusc1,luAvec,lusc2,1,lblk)
        ! lures used as scratch
        call mvcsmd(luvec_sbsp,xvec(1,iminev),lusc1,lures,
     &              vec1,vec2,ndim-1,1,lblk)

        call vecsmd(vec1,vec2,-xlam*xx,1d0,
     &              lusc1,lusc2,lures,1,lblk)
        call vecsmd(vec1,vec2,1d0,-xlam*xx*xvec(ndim,iminev),
     &              lures,luvec,lusc2,1,lblk)
      else
        ! (A - lambda S) x ...
        call vecsmd(vec1,vec2,
     &              xvec(ndim,iminev),-xlam*xx*xvec(ndim,iminev),
     &              luAvec,luvec,lusc2,1,lblk)
      end if

      call copvcd(lusc2,lures,vec1,1,lblk)

      xresnrm = sqrt(inprdd(vec1,vec1,lures,lures,1,lblk))

      call relunit(lusc1,'delete')
      call relunit(lusc2,'delete')

      return
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja_mvp(itype,
     &         luvec,luAvec,
     &         smat,kpiv,ndim,
     &         lust_sbsp,nskipst,lujt_sbsp,nskipjt,
     &         ludia,
     &         vec1,vec2,lincore)
*----------------------------------------------------------------------*
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     lincore

      integer, intent(in) ::
     &     itype, luvec, luAvec,
     &     ndim, lust_sbsp, lujt_sbsp,
     &     nskipst, nskipjt, ludia

      integer, intent(in) ::
     &     kpiv(ndim)

      real(8), intent(in) ::
     &     smat(*)

      real(8), intent(inout) ::
     &     vec1(*), vec2(*) 

      integer ::
     &     lusc1, lusc2, lblk,
     &     ii, jj

      real(8) ::
     &     v(ndim)

      real(8), external ::
     &     inprdd
      integer, external ::
     &     iopen_nus

      lblk = -1

      lusc1 = iopen_nus('MATVSCR1')
      lusc2 = iopen_nus('MATVSCR2')

      call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
      ! v(j) = <x_j|v>
      do jj = 1, ndim
        call rewino(luvec)
        v(jj) = inprdd(vec1,vec2,luvec,lust_sbsp,0,lblk)
      end do

      if (ntest.ge.100) then
        write(6,*) 'v:'
        call wrtmat(v,ndim,1,ndim,1)
        write(6,*) 'smat:'
        call prtrlt(smat,ndim)
      end if

      ! S_{ij}^{-1}v_j
      call dspsl(smat,ndim,kpiv,v)

      if (ntest.ge.100) then
        write(6,*) 'vbar:'
        call wrtmat(v,ndim,1,ndim,1)
      end if

      call skpvcd(lujt_sbsp,nskipjt,vec1,1,lblk)
      call rewino(luAvec)
      call rewino(lusc2)
      ! v(i)*|Jx_i> --> luAvec
      call mvcsmd(lujt_sbsp,v,luAvec,lusc2,vec1,vec2,
     &     ndim,0,lblk)

      if (ntest.ge.1000) then
        write(6,*) 'v(i)|Jx_i>:'
        call wrtvcd(vec1,luAvec,1,lblk)
      end if

      call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
      do ii = 1, ndim
        ! D|x_i> --> lusc2
        call rewino(ludia)
        call rewino(lusc2)
        call dmtvcd2(vec1,vec2,ludia,lust_sbsp,lusc2,1d0,0d0,0,0,lblk)
        ! add -v(i)D|x_i> to luAvec
        call vecsmd(vec1,vec2,1d0,-v(ii),luAvec,lusc2,lusc1,1,lblk)
        call copvcd(lusc1,luAvec,vec1,1,lblk)
      end do

      if (ntest.ge.1000) then
        write(6,*) 'v(i)(|Jx_i>-D|x_i>):'
        call wrtvcd(vec1,luAvec,1,lblk)
      end if

      ! D|v>
      call dmtvcd2(vec1,vec2,ludia,luvec,lusc2,1d0,0d0,1,0,lblk)
      ! add D|v> to luAvec
      call vecsmd(vec1,vec2,1d0,1d0,luAvec,lusc2,lusc1,1,lblk)
      call copvcd(lusc1,luAvec,vec1,1,lblk)

      if (ntest.ge.1000) then
        write(6,*) 'D|v> + v(i)(|Jx_i>-D|x_i>):'
        call wrtvcd(vec1,luAvec,1,lblk)
      end if

      call relunit(lusc1,'delete')
      call relunit(lusc2,'delete')

      return
      
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja_umat(itype,
     &                      umat,nnew,nold,nskip_u,maxvec,
     &                      nskipst,nskipjt,
     &                      lust_sbsp,lujt_sbsp,ludia,xdamp,
     &                      vec1,vec2,lincore)
*----------------------------------------------------------------------*
*
*  update the low-rank hessian/jacobian
*
*   itype = 1:       U_ij = <dt^(i)|(D+damp)^-1|Adt^(j)>
*
*   itype = 2:       dto. but symmetrized
*
*  |dt^(i)>   is on lust_sbsp
*  |Adt^(i)>  is on lujt_sbsp (usually gradient differences)
*   D         is on ludia
*
*  if (xdamp.ne.0d0) the complete matrix needs to be rebuild. we
*  expect the overlap matrix in packed triangular form on input.
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00
      logical, parameter ::
     &     ldo_extra = .true.

* input/output
      logical, intent(in) ::
     &     lincore
      integer, intent(in) ::
     &     itype,
     &     nnew, nold, nskip_u, maxvec,
     &     nskipst, nskipjt,
     &     lust_sbsp, lujt_sbsp, ludia
      real(8), intent(in) ::
     &     xdamp
      real(8), intent(inout) ::
     &     umat(*),
     &     vec1(*), vec2(*)

* local
      integer ::
     &     ii, jj, ijdx, isym, lblk, luscr
      real(8) ::
     &     uij, fac

* external
      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd, inprdd3

      lblk = -1

      if (itype.lt.1.or.itype.gt.2) then
        write(6,*) 'optc_sbspja_umat: not prepared for itype = ',itype
        stop 'optc_sbspja_umat'
      end if
      
      luscr = iopen_nus('UMATSCR')
      
      if (itype.eq.1) isym = 0
      if (itype.eq.2) isym = 1

      if (itype.eq.1) fac = 1d0
      if (itype.eq.2) fac = .5d0

      if (xdamp.ne.0d0.and.ldo_extra) then
        if (nold.ne.0) then
          write(6,*) 'error: nold must be 0 if xdamp.ne.0d0!'
          stop 'optc_sbspja_umat'
        end if
        if (ntest.ge.100) then
          write(6,*) 'on entry: S packed'
          call prtrlt(umat,nnew)
        end if
        if (isym.eq.0) then
          call uptripak(umat,umat,2,nnew,maxvec)
          if (ntest.ge.100) then
            write(6,*) 'S unpacked'
            call wrtmat(umat,nnew,nnew,maxvec,maxvec)
          end if
        end if
      else

        ! zero the new elements:
        if (isym.eq.0) then
          do jj = nskip_u+nold+1, nskip_u+nold+nnew
            do ii = nskip_u+1, nskip_u+nold+nnew
              ijdx = (jj-1)*maxvec + ii
              umat(ijdx) = 0d0
            end do
          end do
          do ii = nskip_u+nold+1, nskip_u+nold+nnew
            do jj = nskip_u+1, nskip_u+nold
              ijdx = (jj-1)*maxvec + ii
              umat(ijdx) = 0d0
            end do
          end do
        else
          do jj = nskip_u+nold+1, nskip_u+nold+nnew
            do ii = nskip_u+1, jj
              ijdx = jj*(jj-1)/2 + ii
              umat(ijdx) = 0d0
            end do
          end do
        end if

      end if ! if (xdamp.ne.0d0)
            
c a) update columns jj 
c first let´s do the damping terms (to see if they are significant)
c    u(i,j) -= <dt^i|(A_0+xdamp)^-1 A_0 |dt^j>       
      if (xdamp.ne.0d0.and.ldo_extra) then
        do jj = nskip_u+nold+1, nskip_u+nold+nnew
          call rewino(ludia)
          call rewino(luscr)
          ! well, not ideal ... :
          call skpvcd(lust_sbsp,nskipst-1+jj,vec1,1,lblk)
          call dmtvcd(vec1,vec2,ludia,lust_sbsp,luscr,0d0,0,0,lblk)
          call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
          do ii = nskip_u+1, jj
            call rewino(luscr)
            call rewino(ludia)
            uij = inprdd3(vec1,vec2,lust_sbsp,luscr,ludia,
     &                    xdamp,-1d0,0,lblk)
            if (isym.eq.0) then
              ijdx = (jj-1)*maxvec + ii
              umat(ijdx) = umat(ijdx) - uij
              if (ii.ne.jj) then
                ijdx = (ii-1)*maxvec + jj
                umat(ijdx) = umat(ijdx) - uij
              end if
            else
              ijdx = jj*(jj-1)/2 + ii
              umat(ijdx) = umat(ijdx) - uij
            end if
          end do
        end do
      end if

      if (ntest.ge.100.and.xdamp.ne.0d0.and.ldo_extra) then
        write(6,*) 'S_ij-<i|A_0(l)^-1 A_0|j>  contrib. to U matrix: '
        if (isym.eq.1) then
          call prtrlt(umat,nold+nnew)
        else
          call wrtmat(umat,nold+nnew,nold+nnew,maxvec,maxvec)
        end if
      end if

      call skpvcd(lujt_sbsp,nskipjt+nold,vec1,1,lblk)
c and now the important contribution:
c    u(i,j) += <dt^i|A_0^-1|Adt^j>       
      do jj = nskip_u+nold+1, nskip_u+nold+nnew
        call rewino(ludia)
        call rewino(luscr)
        call dmtvcd(vec1,vec2,ludia,lujt_sbsp,luscr,xdamp,0,1,lblk)
        call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
        do ii = nskip_u+1, nskip_u+nold+nnew
          call rewino(luscr)
          uij = inprdd(vec1,vec2,lust_sbsp,luscr,0,lblk)
          if (isym.eq.0) then
            ijdx = (jj-1)*maxvec + ii
          else
            if (ii.le.jj) then
              ijdx = jj*(jj-1)/2 + ii
            else
              ijdx = ii*(ii-1)/2 + jj
            end if
          end if
          umat(ijdx) = umat(ijdx) + fac*uij
          if (isym.eq.1.and.ii.eq.jj) umat(ijdx) = umat(ijdx) + fac*uij
        end do
      end do
* b) update rows ii
*    (as we assume nold.eq.0 for xdamp.ne.0d0, here nothing is 
*    to be done for this case)
      if (nold.gt.0) then
        call skpvcd(lust_sbsp,nskipst+nold,vec1,1,lblk)
        do ii = nskip_u+nold+1, nskip_u+nold+nnew
          call rewino(ludia)
          call rewino(luscr)
          call dmtvcd(vec1,vec2,ludia,lust_sbsp,luscr,xdamp,0,1,lblk)
          call skpvcd(lujt_sbsp,nskipjt,vec1,1,lblk)
          do jj = nskip_u+1, nskip_u+nold ! do not doubly visit the elements
                                          ! we already updated in a)
            call rewino(luscr)
            uij = inprdd(vec1,vec2,lujt_sbsp,luscr,0,lblk)
            if (isym.eq.0) then
              ijdx = (jj-1)*maxvec + ii
            else
              if (ii.le.jj) then
                ijdx = jj*(jj-1)/2 + ii
              else
                ijdx = ii*(ii-1)/2 + jj
              end if
            end if
            umat(ijdx) = umat(ijdx) + fac*uij
          end do
        end do

      end if ! nold.gt.0

      if (ntest.ge.100) then
        write(6,*) 'Updated U matrix: ',isym
        if (isym.eq.1) then
          call prtrlt(umat,nold+nnew)
        else
          call wrtmat(umat,nold+nnew,nold+nnew,maxvec,maxvec)
        end if
      end if
      call relunit(luscr,'delete')

      return
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja_prjpstp(vecout,navec,nskipst,
     &                           lu_pertstp,lust_sbsp,
     &                           vec1,vec2,lincore)
*----------------------------------------------------------------------*
*
*     get projection of vector on lu_pertstp in subspace lust_sbsp
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     lincore
      integer, intent(in) ::
     &     navec, nskipst,
     &     lu_pertstp,lust_sbsp
      real(8), intent(out) ::
     &     vecout(navec)
      real(8), intent(inout) ::
     &     vec1(*), vec2(*)

      integer ::
     &     ii, lblk

      real(8), external ::
     &     inprdd

      lblk = -1

      call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
      do ii = 1, navec
        call rewino(lu_pertstp)
        vecout(ii) = inprdd(vec1,vec2,lust_sbsp,lu_pertstp,0,lblk)
      end do

      if (ntest.ge.100) then
        write(6,*) 'projected RHS:'
        call wrtmat(vecout,navec,1,navec,1)
      end if

      return
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_sbspja_lrstep(itype,thrsh,accept,vecin,navec,
     &         nskipst,nskipjt,
     &         lust_sbsp,lujt_sbsp,ludia,xdamp,
     &         lu_pertstp,lu_newstp,
     &         vec1,vec2,lincore)
*----------------------------------------------------------------------*
      implicit none
      
      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     lincore
      logical, intent(out) ::
     &     accept
      integer, intent(in) ::
     &     itype, navec, nskipst, nskipjt,
     &     lust_sbsp, lujt_sbsp, ludia,
     &     lu_pertstp, lu_newstp
      real(8), intent(in) ::
     &     vecin(navec), xdamp, thrsh
      real(8), intent(inout) ::
     &     vec1(*), vec2(*)

      integer ::
     &     lblk, luscr, luscr2
      real(8) ::
     &     xnrm
      
      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd

      lblk = -1

      if (ntest.ge.100) then
        write(6,*) 'assembling low-rank step'
        write(6,*) '------------------------'
        write(6,*) ' navec, nskipst, nskipjt: ',navec, nskipst, nskipjt
        write(6,*) ' lust_sbsp, lujt_sbsp, ludia: ',
     &               lust_sbsp, lujt_sbsp, ludia
        write(6,*) ' lu_pertstp, lu_newstp: ',lu_pertstp, lu_newstp
      end if

      luscr  = iopen_nus('LRCSCR1')
      luscr2 = iopen_nus('LRCSCR2')

      call rewino(luscr)
      call rewino(luscr2)
      call skpvcd(lujt_sbsp,nskipjt,vec1,1,lblk)
      ! Omg = vec(i) * Omg(i) --> luscr
      call mvcsmd(lujt_sbsp,vecin,luscr,luscr2,vec1,vec2,
     &     navec,0,lblk)
      ! D^-1 Omg --> luscr2
      call dmtvcd(vec1,vec2,ludia,luscr,luscr2,xdamp,1,1,lblk)

      call rewino(luscr)
      call rewino(luscr2)
      call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
      ! dt = vec(i) * dt(i) --> lu_newstp (only intermediate scratch)
      call mvcsmd(lust_sbsp,vecin,lu_newstp,luscr,vec1,vec2,
     &     navec,0,lblk)
* in the case of damping, we have to modify the result with A_0/(A_0+xdamp)
      if (xdamp.ne.0d0) then
        call dmtvcd(vec1,vec2,ludia,lu_newstp,luscr,xdamp,1,1,lblk)
        call dmtvcd(vec1,vec2,ludia,luscr,lu_newstp,0d0,1,0,lblk)
      end if
      ! assemble total low-rank correction on luscr
      call vecsmd(vec1,vec2,1d0,-1d0,lu_newstp,luscr2,luscr,1,lblk)

      xnrm = sqrt(inprdd(vec1,vec1,luscr,luscr,1,lblk))

      if (ntest.ge.10) write(6,*) '|low-rank correction:| =',xnrm
      if (xnrm.gt.thrsh) then
        if (ntest.ge.10) write(6,*) 'low-rank correction not accepted!'
        call copvcd(lu_pertstp,lu_newstp,vec1,1,lblk)
        accept = .false.
      else
        call vecsmd(vec1,vec2,1d0,1d0,lu_pertstp,luscr,lu_newstp,1,lblk)
        accept = .true.
      end if

      call relunit(luscr,'delete')
      call relunit(luscr2,'delete')
      
      return
      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_smat(smat,ndim,nold,nskips,
     &                     lu_sbsp1,nskip1,nvec1,
     &                     lu_sbsp2,nskip2,nvec2,
     &                     vec1,vec2,lincore)
*----------------------------------------------------------------------*
*     update or build S matrix = <i|j> as upper triangular matrix
*     ndim     the new S dimension
*     nold     vectors already reside in the smat array
*     nskips   vectors for that the new <i|new> contributions are not
*              needed
*     lu_sbsp1 the first file containing (usually) the old vectors
*              of the previously initialized subspace
*     nskip1   number of vectors to be skipped on that file
*     nvec1    total (!) number of vectors on that file
*     lu_sbsp2 the second file containing (usually) the new vector(s)
*              to be added to the subspace
*     nskip2   same as nskip1 but for second file
*     nvec2    dto.
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     lincore
      integer, intent(in) ::
     &     ndim,nold,nskips,
     &     lu_sbsp1,lu_sbsp2,nskip1,nskip2,nvec1,nvec2
      real(8), intent(inout) ::
     &     smat(*),vec1(*),vec2(*)

      integer ::
     &     ii, jj, idxii, lblk, lusc1,
     &     nadd, nadd1, nadd2, nvecs

      real(8), external ::
     &     inprdd
      integer, external ::
     &     iopen_nus

      if (ntest.ge.100) then
        write(6,*) 'Welcome to optc_smat:'
        write(6,*) '====================='
        write(6,*) ' ndim,nold,nskips: ',ndim,nold,nskips
        write(6,*) ' lu_sbsp1,nskip1,nvec1: ',lu_sbsp1,nskip1,nvec1
        write(6,*) ' lu_sbsp2,nskip2,nvec2: ',lu_sbsp2,nskip2,nvec2
        write(6,*) ' infos on units:'
        if (nvec1.gt.0) then
          print *,'lu_sbsp1:'
          call unit_info(lu_sbsp1)
        end if
        if (nvec2.gt.0) then
          print *,'lu_sbsp2:'
          call unit_info(lu_sbsp2)
        end if
      end if

      lblk = -1
      ! consistency check:
      nvecs = nvec1+nvec2 - (nskip1+nskip2)
      if (nvecs.ne.ndim) then
        write(6,*) 'consistency error in optc_smat: '
        write(6,*) ' nvec1,  nvec2,  sum1 : ',nvec1, nvec2, nvec1+nvec2
        write(6,*) ' nskip1, nskip2, sum2 : ',
     &       nskip1, nskip2, nskip1+nskip2
        write(6,*) ' sum1 - sum2, ndim    : ',nvecs,ndim 
        stop 'optc_smat'
      end if

      lusc1 = iopen_nus('SMATSCR')

      if (ntest.ge.100) then
        write(6,*) 'S on input:'
        call prtrlt(smat,ndim)
      end if

      ! number of rows to be added
      nadd = ndim - nold
      nadd2 = nvec2-nskip2
      nadd1 = nadd-nadd2

      ! contributions within lu_sbsp1
      do ii = nold+1, nold+nadd1
        ! goto the appropriate vector on lu_sbsp1
        call skpvcd(lu_sbsp1,nskip1+ii-1,vec1,1,lblk)
        ! get a copy
        call rewino(lusc1)
        call copvcd(lu_sbsp1,lusc1,vec1,0,lblk)
        ! rewind
        call skpvcd(lu_sbsp1,nskip1,vec1,1,lblk)
        idxii = (ii-1)*ii/2
        do jj = 1+nskips, ii
          call rewino(lusc1)
          smat(idxii+jj) = inprdd(vec1,vec2,lusc1,lu_sbsp1,0,lblk)
        end do
      end do

      if (ntest.ge.100) then
        write(6,*) 'S after 1st block:'
        call prtrlt(smat,ndim)
      end if

      ! contributions from lu_sbsp1/lu_sbsp2
      if (nadd2.gt.0)
     &     call skpvcd(lu_sbsp2,nskip2,vec1,1,lblk)
      do ii = nold+nadd1+1, nold+nadd1+nadd2
        call rewino(lusc1)
        call copvcd(lu_sbsp2,lusc1,vec1,0,lblk)
        call skpvcd(lu_sbsp1,nskip1,vec1,1,lblk)
        idxii = (ii-1)*ii/2
        do jj = 1+nskips, nold+nadd1
          call rewino(lusc1)
          smat(idxii+jj) = inprdd(vec1,vec2,lusc1,lu_sbsp1,0,lblk)
        end do
      end do
      
      if (ntest.ge.100) then
        write(6,*) 'S after 2nd block:'
        call prtrlt(smat,ndim)
      end if

      ! contributions within lu_sbsp2
      do ii = nold+nadd1+1, nold+nadd1+nadd2
        ! go to the appropriate vector on lu_sbsp2
        call skpvcd(lu_sbsp2,nskip1+ii-1,vec1,1,lblk)
        ! get a copy
        call rewino(lusc1)
        call copvcd(lu_sbsp2,lusc1,vec1,0,lblk)
        ! rewind
        call skpvcd(lu_sbsp2,nskip1,vec1,1,lblk)
        idxii = (ii-1)*ii/2
        do jj = nold+nadd1+1, ii
          call rewino(lusc1)
          smat(idxii+jj) = inprdd(vec1,vec2,lusc1,lu_sbsp2,0,lblk)
        end do
      end do

      if (ntest.ge.100) then
        write(6,*) 'S after 3rd block:'
        call prtrlt(smat,ndim)
      end if

      call relunit(lusc1,'delete')

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine optc_updtja(itype,nrank,thrsh,
     &           nstdim,nhgdim,
     &           navec,maxvec,
     &           nadd,navec_last,
     &           lugrvf,lugrvf_last,
     &           ludia,trrad,lu_pertstp,lu_newstp,
     &           luhg_last,luhgam_new,
     &           xdamp,xdamp_last,
     &           lust_sbsp,luhg_sbsp,
     &           hkern,vec1,vec2,iprint)
*----------------------------------------------------------------------*
*
*     calculate the matrix-vector product of the current gradient passed
*     in with an updated inverse Hessian matrix H^-1(k):
*
*      |w> = H^-1(k)|v> = H_0^-1|v> + sum(i=2,k) E_i |v>
* 
*     where the matrices E_i are obtained according to the update
*     formulae of
* 
*      11-14: Broyden-Family:
*      11  Broyden-Fletcher-Goldfarb-Shanno
*      12  Davidon-Fletcher-Powell
*      13  Broyden Rank 1
*      14  Hoshino
*      
*      15  Broyden assymetric rank 1
*
*     |v> comes in on lugrvf
*     |w> goes into the wide world through unit lu_newstp
*
*     |delta_i>    previous steps on lust_sbsp
*     |H gamma_i>  previous grad. diffs. times prev. H on luhg_sbsp
*                  for Powell: |gamma - H^{-1}delta>, instead
*     H_0 = H_1    on ludia
*     H_0^-1|g> is assumed to be on lu_pertstp
*     H_(k-1)^-1|gamma_(k-1)> will leave on luhgnew 
*                  for Powell: |gamma - H_(k-1)^{-1}delta>
*
*     H_(k-1)^-1|g_(k-1)> comes from previous iteration thru luhg_last
*     H_(k)^-1|g_(k)> goes to next iteration thru luhg_last
*                  ( not used for Powell )
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00
      
* input/output
      integer, intent(in) ::
     &     itype,                    ! type (s.above)
     &     nrank,                    ! rank of update
     &     nstdim,nhgdim,            ! # records on sbsp files
     &     maxvec,                   ! max dimension of sbsp
     &     nadd,                     ! number of new vectors (usually 1)
     &     lugrvf,ludia,             ! current gradient, diagonal Hess/Jac.
     &     lugrvf_last,              ! previous gradient
     &     lu_pertstp,               ! step from diagonal Hess/Jac.
     &     lu_newstp,luhgam_new,     ! new step (output)
     &     luhg_last,                ! H_(k-1)^-1|g_(k-1)> from last iteration
     &     lust_sbsp,                ! sbsp files
     &     luhg_sbsp,
     &     iprint                    ! print level
      integer, intent(inout) ::
     &     navec,                    ! current dimension of sbsp
     &     navec_last                ! previous sbsp dimension
      real(8), intent(in) ::
     &     thrsh,                    ! thresh for accepting low-rank correction
     &     trrad,                    ! trust radius for step
     &     xdamp_last                ! previous damping
      real(8), intent(inout) ::
     &     xdamp                     ! current damping
      real(8), intent(inout) ::
     &     hkern(nrank*(nrank+1)/2,navec),
     &                          ! low-rank kernels of previous updates
     &     vec1(*), vec2(*)     ! scratch

* local O(N) scratch
c      integer ::
c     &     kpiv(navec)
      real(8) ::
     &     v1(navec), v2(navec),
     &     hv1(navec), hv2(navec)

* local
c      logical ::
c     &     lincore, again, accept
      integer ::
     &     ii, iprintl, ndel1, ndel2, nn,
     &     nskipst, nskiphg, lblk,
     &     luscr1, luscr2, luscr3, ludlcnt, lupnt1, lupnt2
      real(8) ::
     &     x4, x2, x3, x1, phi,
     &     f1, f2, fac

      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd

      luscr1 = iopen_nus('RUHSCR1') 
      luscr2 = iopen_nus('RUHSCR2') 
      luscr3 = iopen_nus('RUHSCR3') 

      nskipst = nstdim - navec
      nskiphg = nhgdim - (navec-1)

      iprintl = max(iprint,ntest)
      lblk = -1

      if (itype.lt.11.or.itype.gt.15) then
        write(6,*) 'unknown update method!'
        stop 'optc_updtja'
      end if

      if ((itype.ge.11.and.itype.le.14.and.nrank.ne.2).or.
     &    (itype.eq.15                .and.nrank.ne.1)) then
        write(6,*) 'illegal combination of type and rank: ',itype,nrank
        stop 'optc_updtja'
      end if

      if (iprintl.ge.5) then
        write(6,*) 'Updated-Hessian/Jacobian'
        write(6,*) '========================'
        if (itype.eq.11) then
          write(6,*) ' BFGS update'
        else if (itype.eq.12) then
          write(6,*) ' DFP update'
        else if (itype.eq.13) then
          write(6,*) ' Broyden rank 1 update'
        else if (itype.eq.14) then
          write(6,*) ' Hoshino update'
        else if (itype.eq.15) then
          write(6,*) ' Broyden asymmetric update'
        end if
      end if

      if (ntest.ge.10) then
        write(6,*) 'on entry in optc_sbspja:'
        write(6,*)
     &      'nstdim, nhgdim, navec, maxvec, nadd, navec_last: ',
     &       nstdim, nhgdim, navec, maxvec, nadd, navec_last
        write(6,*)
     &       'lugrvf, lugrvf_last, lu_newstp, lu_pertstp: ',
     &        lugrvf, lugrvf_last, lu_newstp, lu_pertstp
        write(6,*)
     &       'luhgam_new, luhg_last: ',luhgam_new, luhg_last
        write(6,*)
     &       'lust_sbsp, luhg_sbsp: ', lust_sbsp, luhg_sbsp
      end if

      if (ntest.ge.10) then
        write(6,*) 'nskipst, nskiphg: ',
     &       nskipst, nskiphg
      end if

      if (ntest.ge.100) then
        write(6,*) 'previous low-rank matrices:'
        do ii = 1, navec-1
          write(6,*) '(',ii,')'
          call prtrlt(hkern(1,ii),nrank)
        end do
      end if

      ndel1 = 0
      ndel2 = 0
      if (navec.eq.maxvec) then
        ! be prepared to remove the first vector from the subspace
        ndel1 = 1  ! for Hinv_{k}g_{k} contributions
        ndel2 = 1  ! for moving h^{k} around
c        if (itype.eq.15) ndel1 = 0 ! not necessary
        if (ndel1.ne.0) ludlcnt = iopen_nus('RUHDLCNT')
      end if

* 1) calculate
*
*    Hinv_{k-1} gamma_{k-1} = Hinv_{k-1} g_{k} - Hinv_{k-1} g_{k-1}
*
*    where
*    
*    Hinv_{k-1} g_(k) = Hinv_{1}g_k 
*                                                    
*      + sum_{i=2,k-1} (delta_{i-1}, Hgamma_{i-1}) x 
*                                                    
*                         / h11^(i) h12^(i) \ / <delta_{i-1}|g_(k)>\
*                         |                 | |                    |
*                         \ h12^(i) h22^(i) / \<Hgamma_{i-1}|g_(k)>/
*
*    the low-rank kernel h^(i) is on hkern(1..3,i)
*    -Hinv_{1}g_k is on lu_pertstp
*    -Hinv_{k-1}g_{k-1} is on luhg_last
*
      ! get this file into position
      if (navec.gt.0) then
        call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
      end if
        
      if (navec.gt.1.and.itype.ge.11.and.itype.le.14) then

        ! <delta_{i}|g_(k)>
        ! this has been done outside the 'if':
        ! call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
        do ii = 1, navec-1
          call rewino(lugrvf)
          v1(ii) = inprdd(vec1,vec2,lust_sbsp,lugrvf,0,lblk)
        end do

        ! <Hinv gamma_{i}|g_(k)>
        call skpvcd(luhg_sbsp,nskiphg,vec1,1,lblk)
        do ii = 1, navec-1
          call rewino(lugrvf)
          v2(ii) = inprdd(vec1,vec2,luhg_sbsp,lugrvf,0,lblk)
        end do

        if (ntest.ge.100) then
          write(6,*) '<delta_{i}|g_(k)>:'
          call wrtmat(v1,navec-1,1,navec-1,1)
          write(6,*) '<Hgamma_{i}|g_(k)>:'
          call wrtmat(v2,navec-1,1,navec-1,1)
        end if

        do ii = 1, navec-1
          hv1(ii) = hkern(1,ii)*v1(ii) + hkern(2,ii)*v2(ii)
          hv2(ii) = hkern(2,ii)*v1(ii) + hkern(3,ii)*v2(ii)
        end do

        if (ntest.ge.100) then
          write(6,*) 'hv1:'
          call wrtmat(hv1,navec-1,1,navec-1,1)
          write(6,*) 'hv2:'
          call wrtmat(hv2,navec-1,1,navec-1,1)
        end if

        ! bring subspace-files into position:
        call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
        call skpvcd(luhg_sbsp,nskiphg,vec1,1,lblk)
        if (ndel1.gt.0) then
          ! store contributions from 'to-be-deleted' vectors on
          ! separate file ludlcnt
          call mvcsmd(lust_sbsp,hv1(1),luscr1,luscr2,
     &                vec1,vec2,ndel1,0,lblk)
          call mvcsmd(luhg_sbsp,hv2(1),luscr2,luscr3,
     &                vec1,vec2,ndel1,0,lblk)
          call vecsmd(vec1,vec2,1d0,1d0,luscr1,luscr2,ludlcnt,1,lblk)
        end if

        ! assemble delta contributions
        call rewino(luscr1)
        call rewino(luscr2)
        call mvcsmd(lust_sbsp,hv1(1+ndel1),luscr1,luscr2,
     &              vec1,vec2,navec-1-ndel1,0,lblk)
c        call mvcsmd(lust_sbsp,hv1,luscr1,luscr2,
c     &              vec1,vec2,navec-1,0,lblk)
        ! add to Hinv_1g_k on lu_pertstp --> luscr3
        call vecsmd(vec1,vec2,1d0,-1d0,lu_pertstp,luscr1,luscr3,1,lblk)

        ! assemble Hinv gamma contributions
        call rewino(luscr1)
        call rewino(luscr2)
        call mvcsmd(luhg_sbsp,hv2(1+ndel1),luscr1,luscr2,
     &              vec1,vec2,navec-1-ndel1,0,lblk)
c        call mvcsmd(luhg_sbsp,hv2,luscr1,luscr2,
c     &              vec1,vec2,navec-1,0,lblk)
        ! add to luscr3 --> lu_newstp (used as intermediate scratch)
        if (ndel1.eq.0) then
          call vecsmd(vec1,vec2,1d0,-1d0,luscr3,luscr1,lu_newstp,1,lblk)
        else
          ! add contribution from 'to-be-deleted' vectors here
          call vecsmd(vec1,vec2,1d0,1d0,luscr1,ludlcnt,luscr2,1,lblk)
          call vecsmd(vec1,vec2,1d0,-1d0,luscr3,luscr2,lu_newstp,1,lblk)
        end if

c        if (itype.ge.11.and.itype.le.14) then
c          ! BFGS and DFP:
c          ! so far, we have built -Hinv_{k-1}g_k on lu_newstp
c          ! combine with -Hinv_{k-1}g_{k-1} on luhg_last to get
c          ! Hinv_{k-1}gamma_{k-1} --> luhgam_new
c          call vecsmd(vec1,vec2,-1d0,1d0,
c     &              lu_newstp,luhg_last,luhgam_new,1,lblk)
c        end if

      else if (navec.gt.1.and.itype.eq.15) then
*
*     assymetric Broyden: we have to recursively update Hinv g
*        
        call skpvcd(luhg_sbsp,nskiphg,vec1,1,lblk)
        ! point to lu_newstp as starting point
        lupnt1 = luscr2
        lupnt2 = lu_pertstp
        ! initial prefactor of -1 as lu_newstp carries -|Hinv0 g>
        fac = -1d0
        do ii = 1, navec-1
          ! <delta^{ii}|Hinv{ii}g>
          call rewino(lupnt2)
          x1 = fac*inprdd(vec1,vec2,lust_sbsp,lupnt2,0,lblk)
          fac = -1d0
          ! multiply with rank-1 kernel
          x2 = hkern(1,ii)*x1

          if (ntest.ge.100) then
            write(6,*) ii,'/',navec-1, '  x1 = ',x1,'  x2 = ',x2
          end if

          if (ii.eq.navec-1) lupnt1 = lu_newstp

          call rewino(lupnt1)
          call rewino(lupnt2)
          call vecsmd(vec1,vec2,1d0,-x2,lupnt2,luhg_sbsp,lupnt1,0,lblk)
          if (mod(ii,2).eq.1) then
            lupnt1 = luscr1
            lupnt2 = luscr2
          else
            lupnt1 = luscr2
            lupnt2 = luscr1
          end if

        end do
        ! now we should have -Hinv^{k-1}g^{k} on lu_newstp (as before)
c        ! combine with old -Hinv^{k-1}g^{k-1}
c        call vecsmd(vec1,vec2,-1d0,1d0,
c     &              lu_newstp,luhg_last,luhgam_new,1,lblk)

      else ! navec.gt.1

        ! first two iterations: just copy
        call copvcd(lu_pertstp,lu_newstp,vec1,1,lblk)
      end if

      if (navec.gt.0) then
        ! so far, we have built -Hinv_{k-1}g_k on lu_newstp
        ! combine with -Hinv_{k-1}g_{k-1} on luhg_last to get
        ! Hinv_{k-1}gamma_{k-1} --> luhgam_new
        call vecsmd(vec1,vec2,-1d0,1d0,
     &              lu_newstp,luhg_last,luhgam_new,1,lblk)

      else
        ! else Hinv_{k-1}gamma{k-1} is identical with Hinv_1 g_{k}
        call sclvcd(lu_pertstp,luhgam_new,-1d0,vec1,1,lblk)
      end if
*
*     what is expected at this position: 
*                 -Hinv^{k-1}g^{k} on lu_newstp
*                 lust_sbsp positioned on delta^{k-1}
*
      if (navec.gt.0) then
        ! we take advantage of the file lust_sbsp being 
        ! the correct postion for this
        ! delta^{k-1} --> luscr2
        call rewino(luscr2)
        call copvcd(lust_sbsp,luscr2,vec1,0,lblk)

      end if

*
* 2) calculate new h^{k}
*
*     we expect that delta^(k-1) is still save(d) on luscr2
*
      if (navec.gt.0.and.(itype.ge.11.and.itype.le.14))
     &     then

        ! calculate the current gamma from current and previous gradient
        ! --> luscr3
        call vecsmd(vec1,vec2,1d0,-1d0,lugrvf,lugrvf_last,luscr3,1,lblk)

        ! <gamma|Hgamma>
        x1 = inprdd(vec1,vec2,luscr3,luhgam_new,1,lblk)
        ! <delta^(k-1)|gamma^{k-1}>
        x2 = inprdd(vec1,vec2,luscr2,luscr3,1,lblk)
        ! <delta^(k-1)|g^{k}>
        x3 = inprdd(vec1,vec2,luscr2,lugrvf,1,lblk)
        ! <Hgamma|g^{k}>
        x4 = inprdd(vec1,vec2,luhgam_new,lugrvf,1,lblk)

        if (itype.eq.11) then
          phi = 1d0 ! BFGS
        else if (itype.eq.12) then
          phi = 0d0 ! DFP
        else if (itype.eq.13) then
          phi = 1d0/(1d0 - x1/x2) ! Broyden rank 1
        else if (itype.eq.14) then
          phi = 1d0/(1d0 + x1/x2) ! Hoshino
        end if

        hkern(1,navec) = (1d0+phi*x1/x2)/x2
        hkern(2,navec) = -phi/x2
        hkern(3,navec) = (phi-1d0)/x1

        if (ntest.ge.100) then
          write(6,*) '<gamma^{k-1}|Hgamma^{k-1}> = ',x1
          write(6,*) ' <delta^(k-1)|gamma^{k-1}> = ',x2
          write(6,*) '       <delta^(k-1)|g^{k}> = ',x3
          write(6,*) '      <Hgamma^{k-1}|g^{k}> = ',x4
        end if
      else if (navec.gt.0.and.itype.eq.15) then
        ! <delta^{k-1}|Hgamma^{k-1}>
        x1 = inprdd(vec1,vec2,luscr2,luhgam_new,1,lblk)
        ! |delta^{k-1}> - |Hgamma^{k-1}> --> luhgam_new
        call vecsmd(vec1,vec2,1d0,-1d0,
     &              luscr2,luhgam_new,luscr3,1,lblk)
        call copvcd(luscr3,luhgam_new,vec1,1,lblk)

        ! <delta^{k-1}|Hg^{k}>
        x3 = -inprdd(vec1,vec2,luscr2,lu_newstp,1,lblk)

        if (ntest.ge.100) then
          print *,'<delta^{k-1}|Hgamma^{k-1}> = ',x1
          print *,'<delta^{k-1}|Hg^{k}>       = ',x3 
        end if

        hkern(1,navec) = 1d0/x1
      end if

      if (navec.gt.0.and.ntest.ge.100) then
        write(6,*) 'new rank n kernel (rank = ',nrank,'):'
        call prtrlt(hkern(1,navec),nrank)
      end if
*
* 3) complete Hinv^{k} g^{k}
*
      if (navec.gt.0.and.nrank.eq.2) then
        f1 = hkern(1,navec) * x3 + hkern(2,navec) * x4
        f2 = hkern(2,navec) * x3 + hkern(3,navec) * x4

        if (ntest.ge.100) then
          write(6,*) 'f1 = ', f1
          write(6,*) 'f2 = ', f2
        end if

        ! add to what we already have on lu_newstp
        call vecsmd(vec1,vec2,1d0,-f1,lu_newstp,luscr2,luscr3,1,lblk)
        call vecsmd(vec1,vec2,1d0,-f2,
     &              luscr3,luhgam_new,lu_newstp,1,lblk)
      else if (navec.gt.0.and.nrank.eq.1) then
        f1 = hkern(1,navec)*x3

        if (ntest.ge.100) then
          write(6,*) 'f1 = ', f1
        end if

        call vecsmd(vec1,vec2,1d0,-f1,
     &       lu_newstp,luhgam_new,luscr3,1,lblk)
        call copvcd(luscr3,lu_newstp,vec1,1,lblk)

      end if ! navec.gt.0
     
      ! save Hinv_k|g_k> for next iteration (only BFGS, DFP)
      if (itype.ge.11.and.itype.le.15) then
        if (ndel1.eq.0) then
          call copvcd(lu_newstp,luhg_last,vec1,1,lblk)
        else
          ! remove contributions from 'to-be-deleted' vectors
          ! such that in the next iteration we really calculate
          ! H^{k-1}gamma{k-1} = H^{k-1}g^k - H^{k-1}g^{k-1}
          ! as this ------------^^^ won't have these vectors anymore
          ! yes, the signs are ok:
          call vecsmd(vec1,vec2,1d0,1d0,
     &              lu_newstp,ludlcnt,luhg_last,1,lblk)
        end if
      end if

      ! if necessary, move h matrices down by ndel2
      if (ndel2.gt.0) then
        nn = nrank*(nrank+1)/2
        do ii = ndel2+1, navec
          hkern(1:nn,ii-ndel2) = hkern(1:nn,ii)
        end do
      end if

      ! remove scratch files
      call relunit(luscr1,'delete')
      call relunit(luscr2,'delete')
      call relunit(luscr3,'delete')
      if (ndel1.gt.0) call relunit(ludlcnt,'delete')

      return

      end
*----------------------------------------------------------------------*
* the version with wrong Powell:
*----------------------------------------------------------------------*
      subroutine optc_updtja_old(itype,nrank,thrsh,
     &           nstdim,nhgdim,
     &           navec,maxvec,
     &           nadd,navec_last,
     &           lugrvf,lugrvf_last,
     &           ludia,trrad,lu_pertstp,lu_newstp,
     &           luhg_last,luhgam_new,
     &           xdamp,xdamp_last,
     &           lust_sbsp,luhg_sbsp,
     &           hkern,vec1,vec2,iprint)
*----------------------------------------------------------------------*
*
*     calculate the matrix-vector product of the current gradient passed
*     in with an updated inverse Hessian matrix H^-1(k):
*
*      |w> = H^-1(k)|v> = H_0^-1|v> + sum(i=2,k) E_i |v>
* 
*     where the matrices E_i are obtained according to the update
*     formulae of
* 
*      11  Broyden-Fletcher-Goldfarb-Shanno
*      12  Davidon-Fletcher-Powell
*      13  Powell symmetric Broyden
*
*     |v> comes in on lugrvf
*     |w> goes into the wide world through unit lu_newstp
*
*     |delta_i>    previous steps on lust_sbsp
*     |H gamma_i>  previous grad. diffs. times prev. H on luhg_sbsp
*                  for Powell: |gamma - H^{-1}delta>, instead
*     H_0 = H_1    on ludia
*     H_0^-1|g> is assumed to be on lu_pertstp
*     H_(k-1)^-1|gamma_(k-1)> will leave on luhgnew 
*                  for Powell: |gamma - H_(k-1)^{-1}delta>
*
*     H_(k-1)^-1|g_(k-1)> comes from previous iteration thru luhg_last
*     H_(k)^-1|g_(k)> goes to next iteration thru luhg_last
*                  ( not used for Powell )
*
*----------------------------------------------------------------------*
      implicit none

* constants
      integer, parameter ::
     &     ntest = 00
      
* input/output
      integer, intent(in) ::
     &     itype,                    ! type (s.above)
     &     nrank,                    ! rank of update
     &     nstdim,nhgdim,            ! # records on sbsp files
     &     maxvec,                   ! max dimension of sbsp
     &     nadd,                     ! number of new vectors (usually 1)
     &     lugrvf,ludia,             ! current gradient, diagonal Hess/Jac.
     &     lugrvf_last,              ! previous gradient
     &     lu_pertstp,               ! step from diagonal Hess/Jac.
     &     lu_newstp,luhgam_new,     ! new step (output)
     &     luhg_last,                ! H_(k-1)^-1|g_(k-1)> from last iteration
     &     lust_sbsp,                ! sbsp files
     &     luhg_sbsp,
     &     iprint                    ! print level
      integer, intent(inout) ::
     &     navec,                    ! current dimension of sbsp
     &     navec_last                ! previous sbsp dimension
      real(8), intent(in) ::
     &     thrsh,                    ! thresh for accepting low-rank correction
     &     trrad,                    ! trust radius for step
     &     xdamp_last                ! previous damping
      real(8), intent(inout) ::
     &     xdamp                     ! current damping
      real(8), intent(inout) ::
     &     hkern(nrank*(nrank+1)/2,navec),
     &                          ! low-rank kernels of previous updates
     &     vec1(*), vec2(*)     ! scratch

* local O(N) scratch
c      integer ::
c     &     kpiv(navec)
      real(8) ::
     &     v1(navec), v2(navec),
     &     hv1(navec), hv2(navec)

* local
c      logical ::
c     &     lincore, again, accept
      integer ::
     &     ii, iprintl, ndel1, ndel2, nn,
     &     nskipst, nskiphg, lblk,
     &     luscr1, luscr2, luscr3, ludlcnt
      real(8) ::
     &     x4, x2, x3, x1, phi,
     &     f1, f2

      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd

      luscr1 = iopen_nus('RUHSCR1') 
      luscr2 = iopen_nus('RUHSCR2') 
      luscr3 = iopen_nus('RUHSCR3') 

      nskipst = nstdim - navec
      nskiphg = nhgdim - (navec-1)

      iprintl = max(iprint,ntest)
      lblk = -1

      if (iprintl.ge.5) then
        write(6,*) 'Updated-Hessian/Jacobian'
        write(6,*) '========================'
        if (itype.eq.11) then
          write(6,*) ' BFGS update'
        else if (itype.eq.12) then
          write(6,*) ' DFP update'
        else if (itype.eq.13) then
          write(6,*) ' Broyden rank 1 update'
        else if (itype.eq.14) then
          write(6,*) ' Hoshino update'
        else if (itype.eq.15) then
          write(6,*) ' PSB update'
        end if
      end if

      if (ntest.ge.10) then
        write(6,*) 'on entry in optc_sbspja:'
        write(6,*)
     &      'nstdim, nhgdim, navec, maxvec, nadd, navec_last: ',
     &       nstdim, nhgdim, navec, maxvec, nadd, navec_last
        write(6,*)
     &       'lugrvf, lugrvf_last, lu_newstp, lu_pertstp: ',
     &        lugrvf, lugrvf_last, lu_newstp, lu_pertstp
        write(6,*)
     &       'luhgam_new, luhg_last: ',luhgam_new, luhg_last
        write(6,*)
     &       'lust_sbsp, luhg_sbsp: ', lust_sbsp, luhg_sbsp
      end if

      if (ntest.ge.10) then
        write(6,*) 'nskipst, nskiphg: ',
     &       nskipst, nskiphg
      end if

      if (ntest.ge.100) then
        write(6,*) 'previous low-rank matrices:'
        do ii = 1, navec-1
          write(6,*) '(',ii,')'
          call prtrlt(hkern(1,ii),nrank)
        end do
      end if

      ndel1 = 0
      ndel2 = 0
      if (navec.eq.maxvec) then
        ! be prepared to remove the first vector from the subspace
        ndel1 = 1  ! for Hinv_{k}g_{k} contributions
        ndel2 = 1  ! for moving h^{k} around
        if (itype.eq.15) ndel1 = 0 ! not necessary
        if (ndel1.ne.0) ludlcnt = iopen_nus('RUHDLCNT')
      end if

* 1) calculate
*
*    Hinv_{k-1} gamma_{k-1} = Hinv_{k-1} g_{k} - Hinv_{k-1} g_{k-1}
*
*    where
*    
*    Hinv_{k-1} g_(k) = Hinv_{1}g_k 
*                                                    
*      + sum_{i=2,k-1} (delta_{i-1}, Hgamma_{i-1}) x 
*                                                    
*                         / h11^(i) h12^(i) \ / <delta_{i-1}|g_(k)>\
*                         |                 | |                    |
*                         \ h12^(i) h22^(i) / \<Hgamma_{i-1}|g_(k)>/
*
*    the low-rank kernel h^(i) is on hkern(1..3,i)
*    -Hinv_{1}g_k is on lu_pertstp
*    -Hinv_{k-1}g_{k-1} is on luhg_last
*
      ! get this file into position
      if (navec.gt.0) then
        call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
      end if
        
      if (navec.gt.1) then

        ! <delta_{i}|g_(k)>
        ! this has been done outside the 'if':
        ! call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
        do ii = 1, navec-1
          call rewino(lugrvf)
          v1(ii) = inprdd(vec1,vec2,lust_sbsp,lugrvf,0,lblk)
        end do

        ! <Hinv gamma_{i}|g_(k)>
        call skpvcd(luhg_sbsp,nskiphg,vec1,1,lblk)
        do ii = 1, navec-1
          call rewino(lugrvf)
          v2(ii) = inprdd(vec1,vec2,luhg_sbsp,lugrvf,0,lblk)
        end do

        if (ntest.ge.100) then
          write(6,*) '<delta_{i}|g_(k)>:'
          call wrtmat(v1,navec-1,1,navec-1,1)
          write(6,*) '<Hgamma_{i}|g_(k)>:'
          call wrtmat(v2,navec-1,1,navec-1,1)
        end if

        do ii = 1, navec-1
          hv1(ii) = hkern(1,ii)*v1(ii) + hkern(2,ii)*v2(ii)
          hv2(ii) = hkern(2,ii)*v1(ii) + hkern(3,ii)*v2(ii)
        end do

        if (ntest.ge.100) then
          write(6,*) 'hv1:'
          call wrtmat(hv1,navec-1,1,navec-1,1)
          write(6,*) 'hv2:'
          call wrtmat(hv2,navec-1,1,navec-1,1)
        end if

        ! bring subspace-files into position:
        call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
        call skpvcd(luhg_sbsp,nskiphg,vec1,1,lblk)
        if (ndel1.gt.0) then
          ! store contributions from 'to-be-deleted' vectors on
          ! separate file ludlcnt
          call mvcsmd(lust_sbsp,hv1(1),luscr1,luscr2,
     &                vec1,vec2,ndel1,0,lblk)
          call mvcsmd(luhg_sbsp,hv2(1),luscr2,luscr3,
     &                vec1,vec2,ndel1,0,lblk)
          call vecsmd(vec1,vec2,1d0,1d0,luscr1,luscr2,ludlcnt,1,lblk)
        end if

        ! assemble delta contributions
        call rewino(luscr1)
        call rewino(luscr2)
        call mvcsmd(lust_sbsp,hv1(1+ndel1),luscr1,luscr2,
     &              vec1,vec2,navec-1-ndel1,0,lblk)
c        call mvcsmd(lust_sbsp,hv1,luscr1,luscr2,
c     &              vec1,vec2,navec-1,0,lblk)
        ! add to Hinv_1g_k on lu_pertstp --> luscr3
        call vecsmd(vec1,vec2,1d0,-1d0,lu_pertstp,luscr1,luscr3,1,lblk)

        ! assemble Hinv gamma contributions
        call rewino(luscr1)
        call rewino(luscr2)
        call mvcsmd(luhg_sbsp,hv2(1+ndel1),luscr1,luscr2,
     &              vec1,vec2,navec-1-ndel1,0,lblk)
c        call mvcsmd(luhg_sbsp,hv2,luscr1,luscr2,
c     &              vec1,vec2,navec-1,0,lblk)
        ! add to luscr3 --> lu_newstp (used as intermediate scratch)
        if (ndel1.eq.0) then
          call vecsmd(vec1,vec2,1d0,-1d0,luscr3,luscr1,lu_newstp,1,lblk)
        else
          ! add contribution from 'to-be-deleted' vectors here
          call vecsmd(vec1,vec2,1d0,1d0,luscr1,ludlcnt,luscr2,1,lblk)
          call vecsmd(vec1,vec2,1d0,-1d0,luscr3,luscr2,lu_newstp,1,lblk)
        end if

        if (itype.ge.11.and.itype.le.14) then
          ! BFGS and DFP:
          ! so far, we have built -Hinv_{k-1}g_k on lu_newstp
          ! combine with -Hinv_{k-1}g_{k-1} on luhg_last to get
          ! Hinv_{k-1}gamma_{k-1} --> luhgam_new
          call vecsmd(vec1,vec2,-1d0,1d0,
     &              lu_newstp,luhg_last,luhgam_new,1,lblk)
        end if

      else ! navec.gt.1
        
        ! else Hinv_{k-1}gamma{k-1} is identical with Hinv_1 g_{k}
        if (itype.ge.11.and.itype.le.14)
     &       call sclvcd(lu_pertstp,luhgam_new,-1d0,vec1,1,lblk)
        call copvcd(lu_pertstp,lu_newstp,vec1,1,lblk)

      end if
*
*     what is expected at this position: 
*                 -Hinv^{k-1}g^{k} on lu_newstp
*                 lust_sbsp positioned on delta^{k-1}
*
      if (navec.gt.0) then
        ! we take advantage of the file lust_sbsp being 
        ! the correct postion for this
        ! delta^{k-1} --> luscr2
        call rewino(luscr2)
        call copvcd(lust_sbsp,luscr2,vec1,0,lblk)

      end if
*
* 1) a.  for Powell-Update:
*
*      calculate H^{k-1} |delta^{k-1}>
*
      if (navec.gt.1.and.itype.eq.15) then

        ! <delta_{i}|delta_(k-1)>
        call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
        do ii = 1, navec-1
          call rewino(luscr2)
          v1(ii) = inprdd(vec1,vec2,lust_sbsp,luscr2,0,lblk)
        end do

        ! <tau_{i}|delta_(k-1)>
        call skpvcd(luhg_sbsp,nskiphg,vec1,1,lblk)
        do ii = 1, navec-1
          call rewino(luscr2)
          v2(ii) = inprdd(vec1,vec2,luhg_sbsp,luscr2,0,lblk)
        end do

        if (ntest.ge.100) then
          write(6,*) '<delta_{i}|delta_(k-1)>:'
          call wrtmat(v1,navec-1,1,navec-1,1)
          write(6,*) '<tau_{i}|delta_(k-1)>:'
          call wrtmat(v2,navec-1,1,navec-1,1)
        end if

        do ii = 1, navec-1
          hv1(ii) = hkern(1,ii)*v1(ii) + hkern(2,ii)*v2(ii)
          hv2(ii) = hkern(2,ii)*v1(ii) + hkern(3,ii)*v2(ii)
        end do

        if (ntest.ge.100) then
          write(6,*) 'hv1:'
          call wrtmat(hv1,navec-1,1,navec-1,1)
          write(6,*) 'hv2:'
          call wrtmat(hv2,navec-1,1,navec-1,1)
        end if

        ! bring subspace-files into position again:
        call skpvcd(lust_sbsp,nskipst,vec1,1,lblk)
        call skpvcd(luhg_sbsp,nskiphg,vec1,1,lblk)

        ! assemble delta contributions
        call rewino(luscr1)
        call rewino(luscr2)
        call mvcsmd(lust_sbsp,hv1,luscr1,luscr3,
     &              vec1,vec2,navec-1,0,lblk)

        ! assemble gamma - Hinv delta contributions
        call rewino(luscr1)
        call rewino(luscr2)
        ! luhgam_new used as scratch:
        call mvcsmd(luhg_sbsp,hv2,luscr3,luhgam_new,
     &              vec1,vec2,navec-1,0,lblk)
        ! add luscr1 and luscr3 --> luhgam_new
        call vecsmd(vec1,vec2,1d0,1d0,luscr1,luscr3,luhgam_new,1,lblk)

        ! luhgam_new is still missing the contributions
        !   H_0^-1 delta^{k-1}  and  gamma^{k-1}
        ! coming soon ...

      end if

*
* 2) calculate new h^{k}
*
*     we expect that delta^(k-1) is still save(d) on luscr2
*
      if (navec.gt.0.and.(itype.ge.11.and.itype.le.14))
     &     then

        ! calculate the current gamma from current and previous gradient
        ! --> luscr3
        call vecsmd(vec1,vec2,1d0,-1d0,lugrvf,lugrvf_last,luscr3,1,lblk)

        ! <gamma|Hgamma>
        x1 = inprdd(vec1,vec2,luscr3,luhgam_new,1,lblk)
        ! <delta^(k-1)|gamma^{k-1}>
        x2 = inprdd(vec1,vec2,luscr2,luscr3,1,lblk)
        ! <delta^(k-1)|g^{k}>
        x3 = inprdd(vec1,vec2,luscr2,lugrvf,1,lblk)
        ! <Hgamma|g^{k}>
        x4 = inprdd(vec1,vec2,luhgam_new,lugrvf,1,lblk)

        if (itype.eq.11) then
          phi = 1d0 ! BFGS
        else if (itype.eq.12) then
          phi = 0d0 ! DFP
        else if (itype.eq.13) then
          phi = 1d0/(1d0 - x1/x2) ! Broyden rank 1
        else if (itype.eq.14) then
          phi = 1d0/(1d0 + x1/x2) ! Hoshino
        end if

        hkern(1,navec) = (1d0+phi*x1/x2)/x2
        hkern(2,navec) = -phi/x2
        hkern(3,navec) = (phi-1d0)/x1

        if (ntest.ge.100) then
          write(6,*) '<gamma^{k-1}|Hgamma^{k-1}> = ',x1
          write(6,*) ' <delta^(k-1)|gamma^{k-1}> = ',x2
          write(6,*) '       <delta^(k-1)|g^{k}> = ',x3
          write(6,*) '      <Hgamma^{k-1}|g^{k}> = ',x4
        end if
      else if (navec.gt.0.and.itype.eq.15) then

        if (navec.gt.1) then
          ! Hinv^{1} delta^{k-1} on luscr3
          call dmtvcd(vec1,vec2,ludia,luscr2,luscr3,0d0,1,1,lblk)
          ! add to luhgam_new --> luscr1
          call vecsmd(vec1,vec2,1d0,1d0,luhgam_new,luscr3,luscr1,1,lblk)
        else
          ! no previous contributions?
          ! Hinv^{1} delta^{k-1} directly on luscr1
          call dmtvcd(vec1,vec2,ludia,luscr2,luscr1,0d0,1,1,lblk)
        end if

        ! gamma^{k-1} on luscr3
        call vecsmd(vec1,vec2,1d0,-1d0,lugrvf,lugrvf_last,luscr3,1,lblk)

        ! subtract luscr1 --> final result on luhgam_new
        call vecsmd(vec1,vec2,-1d0,1d0,luscr1,luscr3,luhgam_new,1,lblk)

        !   <delta^{k-1}|gamma^{k-1}-Hinv^{k-1}delta^{k-1}>
        ! = <delta^{k-1}|tau^{k-1}>
        x1 = inprdd(vec1,vec2,luscr2,luhgam_new,1,lblk)
        ! <delta^{k-1}|delta^{k-1}>
        x2 = inprdd(vec1,vec1,luscr2,luscr2,1,lblk)
        ! <delta^{k-1}|g^k>
        x3 = inprdd(vec1,vec2,luscr2,lugrvf,1,lblk)
        ! <tau^{k-1}|g^k>
        x4 = inprdd(vec1,vec2,luhgam_new,lugrvf,1,lblk)

c        hkern(1,navec) = -x1/(x2*x2)
c        hkern(2,navec) = 1d0/x2
c        hkern(3,navec) = 0d0
c TEST
        hkern(1,navec) = -x1
        hkern(2,navec) = x2
        hkern(3,navec) = 0d0

        if (ntest.ge.100) then
          write(6,*) '  <delta^{k-1}|tau^{k-1}> = ',x1
          write(6,*) '<delta^(k-1)|delta^{k-1}> = ',x2
          write(6,*) '      <delta^(k-1)|g^{k}> = ',x3
          write(6,*) '        <tau^{k-1}|g^{k}> = ',x4
        end if

      end if

      if (navec.gt.0.and.ntest.ge.100) then
        write(6,*) 'new rank n kernel (rank = ',nrank,'):'
        call prtrlt(hkern(1,navec),nrank)
      end if
*
* 3) complete Hinv^{k} g^{k}
*
      if (navec.gt.0) then
        f1 = hkern(1,navec) * x3 + hkern(2,navec) * x4
        f2 = hkern(2,navec) * x3 + hkern(3,navec) * x4

        if (ntest.ge.100) then
          write(6,*) 'f1 = ', f1
          write(6,*) 'f2 = ', f2
        end if

        ! add to what we already have on lu_newstp
        call vecsmd(vec1,vec2,1d0,-f1,lu_newstp,luscr2,luscr3,1,lblk)
        call vecsmd(vec1,vec2,1d0,-f2,
     &              luscr3,luhgam_new,lu_newstp,1,lblk)
      end if ! navec.gt.0
     
      ! save Hinv_k|g_k> for next iteration (only BFGS, DFP)
      if (itype.ge.11.and.itype.le.14) then
        if (ndel1.eq.0) then
          call copvcd(lu_newstp,luhg_last,vec1,1,lblk)
        else
          ! remove contributions from 'to-be-deleted' vectors
          ! such that in the next iteration we really calculate
          ! H^{k-1}gamma{k-1} = H^{k-1}g^k - H^{k-1}g^{k-1}
          ! as this ------------^^^ won't have these vectors anymore
          ! yes, the signs are ok:
          call vecsmd(vec1,vec2,1d0,1d0,
     &              lu_newstp,ludlcnt,luhg_last,1,lblk)
        end if
      end if

      ! if necessary, move h matrices down by ndel2
      if (ndel2.gt.0) then
        nn = nrank*(nrank+1)/2
        do ii = ndel2+1, navec
          hkern(1:nn,ii-ndel2) = hkern(1:nn,ii)
        end do
      end if

      ! remove scratch files
      call relunit(luscr1,'delete')
      call relunit(luscr2,'delete')
      call relunit(luscr3,'delete')
      if (ndel1.gt.0) call relunit(ludlcnt,'delete')

      return

      end
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine cmbamp(imode,lu1,lu2,lu3,lu_comb,
     &     vec,namp1,namp2,namp3)
*----------------------------------------------------------------------*
*
*  imode = 11 : combine vectors on lu1, lu2, lu3 into one single vector
*  imode = 01 : reverse procedure
*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     imode,
     &     lu1, lu2, lu3, lu_comb,
     &     namp1, namp2, namp3
      real(8), intent(inout) ::
     &     vec(*)

      logical ::
     &     l_end
      integer ::
     &     ipass,
     &     lblk, lu, iamzero, iampacked, namp, namp_r, namp_sum!, imone
      real(8), external ::
     &     inprod

      lblk = -1

      if (ntest.ge.10) then
        write(6,*) 'cmbamp at work!'
        write(6,*) '==============='
        write(6,*) 'imode = ',imode
        if (imode.lt.12)write(6,*) 'lu1, lu2, lu3 = ', lu1, lu2, lu3
        write(6,*) 'lu_comb = ',lu_comb
      end if

      if (imode.ge.11) then
* lu1, lu2, lu3         --> lu_comb

        call rewino(lu_comb)
        do ipass = 1, 3
          if (ipass.eq.1) lu = lu1
          if (ipass.eq.2) lu = lu2
          if (ipass.eq.3) lu = lu3
          if (ipass.eq.1) namp = namp1
          if (ipass.eq.2) namp = namp2
          if (ipass.eq.3) namp = namp3
          if (ntest.ge.50) then
            write(6,*) 'ipass, lu, namp: ', ipass, lu, namp
          end if
          if (lu.gt.0) then
            call rewino(lu)
            namp_sum = 0
            l_end = .false.
            do while(.not.l_end)
              call ifrmds(namp_r,1,lblk,lu)
              if (namp_r.eq.-1) then
                if (ntest.ge.50) write(6,*) 'end of vector on lu ',lu
                l_end = .true.
              else
                namp_sum = namp_sum + namp_r
                if (ntest.ge.50)
     &               write(6,*) 'transfer ',namp_r,
     &               ' words from lu ',lu,' to lu_comb ',lu_comb
                call frmdsc(vec,namp_r,lblk,lu,iamzero,iampacked)
                if (ntest.ge.100)
     &               write(6,*) ' norm of that block: ',
     &               sqrt(inprod(vec,vec,namp_r))
                call itods(namp_r,1,lblk,lu_comb)
                call todsc(vec,namp_r,lblk,lu_comb)
              end if
            end do
            if (namp_sum.ne.namp) then
              write(6,*) 'WARNING from cmbamp:'
              write(6,'(x,a,i2,2(a,i10))') ' File ',ipass,
     &             ': read ',namp_sum,' expected: ',namp
            end if
          end if
        end do ! ipass
        ! mark end of record
        call itods(-1,1,lblk,lu_comb)

      else
* lu_comb --> lu1, lu2

        call rewino(lu_comb)
        pass_loop: do ipass = 1, 3
          if (ipass.eq.1) lu = lu1
          if (ipass.eq.2) lu = lu2
          if (ipass.eq.3) lu = lu3
          if (ipass.eq.1) namp = namp1
          if (ipass.eq.2) namp = namp2
          if (ipass.eq.3) namp = namp3

          if (lu.gt.0) then
            call rewino(lu)
            namp_sum = 0
            l_end = .false.
            do while(.not.l_end)
              call ifrmds(namp_r,1,lblk,lu_comb)
              if (namp_r.eq.-1) then
                write(6,*)
     &               'WARNING from cmbamp: Unexpected end of vector'
                write(6,'(x,a,i2,2(a,i10))') ' File ',ipass,
     &             ': read ',namp_sum,' expected: ',namp
                exit pass_loop
              else if (namp_r+namp_sum.gt.namp) then
                write(6,*)
     &               'WARNING from cmbamp: Too long new record ',namp_r
                write(6,'(x,a,i2,2(a,i10))') ' File ',ipass,
     &             ': read ',namp_r,' expected: ',namp
                ! well, this is fatal
                stop 'cmbamp'
                l_end = .true.
              else
                namp_sum = namp_sum + namp_r
                if (ntest.ge.50)
     &               write(6,*) 'transfer ',namp_r,
     &               ' words from lu_comb ',lu_comb,' to lu ',lu
                call frmdsc(vec,namp_r,lblk,lu_comb,iamzero,iampacked)
                if (ntest.ge.100)
     &               write(6,*) ' norm of that block: ',
     &               sqrt(inprod(vec,vec,namp_r))
                call itods(namp_r,1,lblk,lu)
                call todsc(vec,namp_r,lblk,lu)
                if (namp_sum.eq.namp) then
                  if (ntest.ge.50)
     &               write(6,*) 'vector on lu ',lu,' complete'
                  l_end = .true.
                  call itods(-1,1,lblk,lu)
                end if
              end if
            end do ! while (.not.l_end)
          end if ! lu.gt.0
        end do pass_loop
        ! skip the -1
c        call ifrmds(imone,1,lblk,lu_comb)

      end if

      return
      end
*----------------------------------------------------------------------*
* follow: routines for 2nd-order optimization
*----------------------------------------------------------------------*
      subroutine optc_redh(isymm,hred,ndim,ndim_o,
     &     luvec,lumv,luvec_sbsp,lumv_sbsp,
     &     vec1,vec2)
*----------------------------------------------------------------------*
*     update/set up reduced H matrix
*
*      isymm != 0 :  symmetrize (for more numerical stability)
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     isymm,
     &     ndim, ndim_o,
     &     luvec, lumv, luvec_sbsp, lumv_sbsp

      real(8), intent(inout) ::
     &     hred(ndim*ndim), vec1(*), vec2(*)

      real(8) ::
     &     xel, fac
      integer ::
     &     ii, jj,
     &     lblk, luscr

      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd

      lblk = -1

      if (ntest.ge.10) then
        write(6,*) '-----------'
        write(6,*) ' optc_redh'
        write(6,*) '-----------'
        write(6,*) ' ndim, ndim_o : ',ndim, ndim_o
        write(6,*) ' luvec,lumv,luvec_sbsp,lumv_sbsp: ',
     &                luvec,lumv,luvec_sbsp,lumv_sbsp
      end if
      
      if (ndim_o.gt.0) then
        if (ntest.ge.100) then
          write(6,*) 'H(red) on input:'
          call wrtmat2(hred,ndim_o,ndim_o,ndim_o,ndim_o)
        end if
      
        do ii = ndim_o, 1, -1
          do jj = ndim_o, 1, -1
            hred((ii-1)*ndim+jj) = hred((ii-1)*ndim_o+jj)
          end do
        end do

        if (ntest.ge.100) then
          write(6,*) 'H(red) on input (resorted):'
          call wrtmat2(hred,ndim,ndim,ndim,ndim)
        end if

      end if

      if (ndim-ndim_o.gt.1) then
        stop 'i shall not enter this section!'
        luscr = iopen_nus('HRED_SCR')

        call skpvcd(luvec_sbsp,ndim_o,vec1,1,lblk)
        do ii = ndim_o+1, ndim-1
          call copvcd(luvec_sbsp,luscr,vec1,0,lblk)
          call rewino(lumv_sbsp)
          do jj = 1, ndim-1
            call rewino(luscr)
            xel = inprdd(vec1,vec2,luscr,lumv_sbsp,0,lblk)
c            hred((ii-1)*ndim+jj) = xel
            hred((jj-1)*ndim+ii) = xel
          end do
        end do
        
        call skpvcd(lumv_sbsp,ndim_o,vec1,1,lblk)
        do jj = ndim_o+1, ndim-1
          call copvcd(lumv_sbsp,luscr,vec1,0,lblk)
          call rewino(luvec_sbsp)
          do ii = 1, ndim_o
            call rewino(luscr)
            xel = inprdd(vec1,vec2,luscr,luvec_sbsp,0,lblk)
c            hred((ii-1)*ndim+jj) = xel
            hred((jj-1)*ndim+ii) = xel
          end do
        end do

        call relunit(luscr,'delete')

      end if

      if (ndim.gt.1) then
        fac = 1d0
        if (isymm.ne.0) fac = 0.5d0
        call rewino(lumv_sbsp)
        do jj = 1, ndim-1
          call rewino(luvec)
          xel = inprdd(vec1,vec2,luvec,lumv_sbsp,0,lblk)
c          hred((ndim-1)*ndim+jj) = xel
          hred((jj-1)*ndim+ndim) = fac*xel
          if (isymm.ne.0)
     &         hred((ndim-1)*ndim+jj) = fac*xel
        end do
        if (isymm.eq.0) hred((ndim-1)*ndim+1:(ndim-1)*ndim+ndim) = 0d0
        call rewino(luvec_sbsp)
        do ii = 1, ndim-1
          call rewino(lumv)
          xel = inprdd(vec1,vec2,lumv,luvec_sbsp,0,lblk)
          hred((ndim-1)*ndim+ii) = hred((ndim-1)*ndim+ii) + fac*xel
          if (isymm.ne.0)
     &         hred((ii-1)*ndim+ndim) = hred((ii-1)*ndim+ndim) +fac*xel
        end do
      end if

      xel = inprdd(vec1,vec2,lumv,luvec,1,lblk)
      hred(ndim*ndim) = xel

      if (ntest.ge.100) then
        write(6,*) 'Final H(red):'
        call wrtmat2(hred,ndim,ndim,ndim,ndim)
      end if

      return
      end
*----------------------------------------------------------------------*
      subroutine optc_trnewton(imode,iret,isymmet,
     &     hred,gred,cred,scr,scr2,ndim,xlamb,gamma,trrad,de_pred)
*----------------------------------------------------------------------*
*     solve reduced trust-radius newton equations
*
*      imode = 1:  try bare newton step first
*                  the subspace Hessian will be tested for negative
*                  eigenvalues; if so, iret!=0 will signal that this
*                  does not work
*      imode = 2:  solve augmented Hessian equations with adapted damping
*                  (reminiscent of the NEO, see "the book") If a too
*                  low gamma is found (which is, together with the lowest
*                  eigenvalue of the augmented Hessian, nu, related to
*                  the usual damping parameter lambda = nu*gamma**2)
*                  the routine will use the lowest possible gamma and
*                  suggest to take the newton step by setting iret!=0
*
*----------------------------------------------------------------------*
      
      implicit none

      integer, parameter ::
     &     ntest = 00, maxit = 400
      real(8), parameter ::
     &     thrsh = 1d-6, xmxstp = 1.0d0, gamma_min = 1d-1

      integer, intent(in) ::
     &     ndim, imode
      real(8), intent(in) ::
     &     hred(ndim*ndim),gred(ndim),
     &     trrad

      integer, intent(inout) ::
     &     isymmet
      integer, intent(out) ::
     &     iret
      real(8), intent(out) ::
     &     de_pred,xlamb, cred(ndim)
      real(8), intent(inout) ::
     &     gamma,
     &     scr((ndim+1)*(ndim+1)),scr2((ndim+1)*(ndim+1))

      logical ::
     &     opt, conv
      integer ::
     &     iter, ii, jj, idx_min, ierr, irg
      real(8) ::
     &     eig_min, dlt, fac, xel,
     &     xnrm, xnrm_prev, xdum, gamma_prev, gamma_new,
     &     gamma_high, gamma_low
      ! local O(N) scratch
      integer ::
     &     iscr(ndim+1)
      real(8) ::
     &     xscr(ndim+1), eigr(ndim+1), eigi(ndim+1)

      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd, inprod

      if (ntest.ge.10) then
        write(6,*) '---------------'
        write(6,*) ' optc_trnewton'
        write(6,*) '---------------'
        write(6,*) ' imode = ',imode
        write(6,*) ' ndim  = ',ndim
        write(6,*) ' trrad = ',trrad
      end if

      iret = 0 ! start optimistic

      if (imode.eq.1) then
        if (ntest.ge.10) then
          write(6,*) '--------------------'
          write(6,*) ' trying Newton step'
          write(6,*) '--------------------'
        end if

        ! get a copy of hred
        scr(1:ndim*ndim) = hred(1:ndim*ndim)
          
        ! and solve for eigenvalues
        irg = 0 
        call rg(ndim,ndim,scr,eigr,eigi,irg,xdum,iscr,xscr,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Internal ERROR: rg gives ierr = ',ierr
          stop 'optc_trnewton: rg'
        end if
        
        if (ntest.ge.100) then
          write(6,*) 'eigenvalues of Hred:'
          do ii = 1, ndim
            write(6,*) ii,eigr(ii),eigi(ii)
          end do
        end if
        
        eig_min = 0d0
        do ii = 1, ndim
          eig_min = min(eig_min,eigr(ii))
        end do

        if (eig_min.lt.0d0) then
          write(6,*) 'detected negative eigenvalue in subspace'
          write(6,*) 'you are clearly not in the local region!'
          write(6,*) 'I will switch to NEO algorithm'
          iret = 1
        end if

        scr(1:ndim*ndim) = hred(1:ndim*ndim)
        xscr(1:ndim) = -gred(1:ndim)
        ! solve reduced LEQ
        call dgefa(scr,ndim,ndim,iscr,ierr)
        if (ierr.ne.0) then
          write(6,*) 'Internal ERROR: dgefa gives ierr = ',ierr
          stop 'optc_trnewton: dgefa'
        end if
        call dgesl(scr,ndim,ndim,iscr,xscr,0)
        
        xnrm = sqrt(inprod(xscr,xscr,ndim))

        if (ntest.ge.100) then
          write(6,*) 'the Newton step: '
          call wrtmat(xscr,ndim,1,ndim,1)
          write(6,*) ' norm = ',xnrm
          write(6,*) ' trrad = ',trrad
        end if

        if (xnrm.gt.trrad) then
          write(6,*) 'detected too long step!'
          write(6,*) 'you are clearly not in the local region!'
          write(6,*) 'I will switch to NEO algorithm'
          iret = 1
        end if

        if (iret.eq.0) then
          cred(1:ndim) = xscr(1:ndim)
          ! energy prediction
          call matvcb(hred,cred,xscr,
     &         ndim,ndim,0)
          de_pred = inprod(gred,cred,ndim) +
     &        0.5d0*inprod(xscr,cred,ndim)
        end if
        xlamb = 0d0

      end if

      if (imode.eq.2) then

        if (ntest.ge.10) then
          write(6,*) '--------------------------------------'
          write(6,*) ' solving Newton eigenvector equations'
          write(6,*) '--------------------------------------'
        end if
        
        ! gamma is set outside
c TEST
c          gamma = 20d0
c TEST
        ! signal that no previous gamma is available
        gamma_prev = -1d0
        gamma_low = gamma_min
        gamma_high = 1d20
        iter = 0
        do
          iter = iter + 1

          do ! trial loop: 1 w/o symmetrization, 2 w/ symmetrization
            ! set up augmented Hessian (g scaled with gamma)
            do jj = 1, ndim
              do ii = 1, ndim
                scr((jj-1)*(ndim+1)+ii) = hred((jj-1)*ndim+ii)
              end do
              scr((jj-1)*(ndim+1)+ndim+1) = gamma*gred(jj)
            end do
            do ii = 1, ndim
              scr(ndim*(ndim+1)+ii) = gamma*gred(ii)  
            end do
            scr((ndim+1)*(ndim+1))=0d0

            if (isymmet.eq.1) then
              do ii = 1, ndim
                do jj = 1, ii-1
                  xel = 0.5d0*(scr((ii-1)*(ndim+1)+jj)+
     &                 scr((jj-1)*(ndim+1)+ii))
                  scr((ii-1)*(ndim+1)+jj)=xel
                  scr((jj-1)*(ndim+1)+ii)=xel
                end do
              end do
            end if

            if (ntest.ge.100) then
              write(6,*) 'gradient scaled augmented Hessian (gamma = ',
     &                                                     gamma,')'
              call wrtmat2(scr,ndim+1,ndim+1,ndim+1,ndim+1)
            end if

            irg = 1
c          call test_symmat(scr,ndim+1,ndim+1)
            call rg(ndim+1,ndim+1,scr,eigr,eigi,irg,scr2,iscr,xscr,ierr)
            if (ierr.ne.0) then
              write(6,*) 'Internal ERROR: rg gives ierr = ',ierr
              stop 'optc_trnewton: rg'
            end if
          
            if (ntest.ge.100) then
              write(6,*) 'eigenvalues of Hred:'
              do ii = 1, ndim+1
                write(6,*) ii,eigr(ii),eigi(ii)
              end do
            end if

            eig_min = 100d0
            idx_min = 1
            do ii = 1, ndim+1
              if (eig_min.gt.eigr(ii)) then
                eig_min = eigr(ii)
                idx_min = ii
              end if
            end do
c            if (abs(eigi(idx_min)).gt.epsilon(1d0)) then
c            write(6,*) 'Help, imaginary lowest eigenvalue!'
c            stop 'trnewton'
c            end if
            if (abs(eigi(idx_min)).lt.epsilon(1d0)) exit
                      
            if (abs(eigi(idx_min)).gt.epsilon(1d0)
     &               .and.isymmet.eq.1) then
              ! then something is REALLY wrong
              write(6,*)
     &             'Strange, imaginary lowest eigenvalue persists!'
              stop 'trnewton'
            end if

            write(6,*) '>> symmetrizing trick (auweia)'
            isymmet = 1 ! OK, so we try to symmetrize (desperate trick)

          end do
          isymmet = 0

          ! rescale eigenvector
          xscr(1:ndim) = (1d0/gamma)*scr2((idx_min-1)*(ndim+1)+1:
     &                                    (idx_min-1)*(ndim+1)+ndim)
          xscr(ndim+1) = scr2((idx_min-1)*(ndim+1)+ndim+1)
          if (ntest.ge.100) then
            write(6,*) 'the raw eigenvector to ',eig_min
            call wrtmat(xscr,1,ndim+1,1,ndim+1)
          end if

          ! normalize such that the last element is 1
          fac = xscr(ndim+1)
          call scalve(xscr,1d0/fac,ndim+1)

          if (ntest.ge.100) then
            write(6,*) 'the renormalized eigenvector to ',eig_min
            call wrtmat(xscr,1,ndim+1,1,ndim+1)
          end if

          ! get step length (only the first ndim elements)
          xnrm = sqrt(inprod(xscr,xscr,ndim))

          if (ntest.ge.100) then
            write(6,'(x,a,i4,2(x,e12.6))')
     &           ' iter, gamma, step length: ', iter, gamma, xnrm
            write(6,'(x,a,18x,e12.6)')
     &           ' trrad =                   ', trrad
          end if

          if (abs(xnrm-trrad).gt.thrsh) then
            if (xnrm.gt.trrad) then
              ! save highest gamma for that so far xnrm.gt.trrad
              ! i.e. the low limit on gamma
              gamma_low = max(gamma,gamma_low)
              if (gamma_prev.eq.-1d0) then
                gamma_new = gamma + 0.1d0
              else
                dlt = -(gamma-gamma_prev)/(xnrm-xnrm_prev)*(xnrm-trrad)
                if (dlt.lt.0d0) dlt = xmxstp
                gamma_new = gamma + min(xmxstp,dlt)
              end if
            else
              ! save lowest gamma for that so far xnrm.lt.trrad
              ! i.e. the high limit on gamma
              gamma_high = min(gamma,gamma_high)
              if (gamma.eq.gamma_min) then
                write(6,*) 'arrived at lowest possible gamma;'
                write(6,*) 'suggesting to take newton step instead!'
                iret = 1
                exit
              end if
c              if (abs(xnrm-trrad).lt.0.1d0.and.gamma_prev.ne.-1d0) then
              if (gamma_prev.ne.-1d0) then
                dlt = -(gamma-gamma_prev)/(xnrm-xnrm_prev)*(xnrm-trrad)
                if (dlt.gt.0d0) dlt = -xmxstp
                dlt = max(-xmxstp,dlt)
              else
                dlt = -0.1d0
              end if
              gamma_new = max(gamma_min,gamma + dlt)
            end if
            ! stabilize by bracketing
            if (gamma_new.gt.gamma_high.or.
     &           gamma_new.lt.gamma_low) then
              gamma_new = 0.5*(gamma_high+gamma_low)
            end if
            if (ntest.ge.100) then
              write(6,*) 'optimization of gamma:'
              write(6,*) ' gamma_high, gamma_low: ',gamma_high,gamma_low
              write(6,*) ' previous gamma, step: ',gamma_prev,xnrm_prev
              write(6,*) ' current gamma, step : ',gamma,xnrm
              write(6,*) ' new gamma: ',gamma_new
            end if
            gamma_prev = gamma
            gamma = gamma_new
            xnrm_prev = xnrm
          else
            if (ntest.ge.100) write(6,*) 'exited optim. loop'
            exit
          end if

          if (iter.gt.maxit) then
            write(6,*) ' problem in optimizing gamma'
            stop 'trnewton'
          end if

        end do
        
        cred(1:ndim) = xscr(1:ndim)
c        de_pred = 0.5d0*eig_min/(gamma*gamma)
c        print *,' energy prediction 1: ',0.5d0*eig_min/(gamma*gamma)
        call matvcb(hred,cred,xscr,
     &       ndim,ndim,0)
        de_pred = inprod(gred,cred,ndim) +
     &       0.5d0*inprod(xscr,cred,ndim)
c        print *,' energy prediction 2: ',de_pred

        xlamb = eig_min!*gamma*gamma

      end if

      if (ntest.ge.100) then
        write(6,*) 'final solution for lambda = ',xlamb
        call wrtmat(cred,ndim,1,ndim,1)
        write(6,*) ' norm = ',xnrm
        write(6,*) ' trrad = ',trrad
      end if

      return
      end
*----------------------------------------------------------------------*
      subroutine optc_trn_resid(cred,ndim,xlamb,
     &     lures,xresnrm,
     &     luvec_sbsp,lumv_sbsp,lugrvf,
     &     vec1,vec2)
*----------------------------------------------------------------------*
*     set up residual vector in full space
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00
      
      integer, intent(in) ::
     &     ndim, lures, luvec_sbsp, lumv_sbsp, lugrvf
      real(8), intent(in) ::
     &     xlamb, cred(ndim)
      real(8), intent(out) ::
     &     xresnrm
      real(8), intent(inout) ::
     &     vec1(*), vec2(*)

      integer ::
     &     lblk, luscr1, luscr2

      integer, external ::
     &     iopen_nus
      real(8), external ::
     &     inprdd

      lblk = -1

      luscr1 = iopen_nus('RES_SCR1')
      
      call mvcsmd(lumv_sbsp,cred,luscr1,lures,vec1,vec2,ndim,1,lblk)

      if (abs(xlamb).gt.epsilon(1d0)) then
        luscr2 = iopen_nus('RES_SCR2')
        call mvcsmd(luvec_sbsp,cred,lures,luscr2,vec1,vec2,ndim,1,lblk)
        call vecsmd(vec1,vec2,-xlamb,1d0,lures,luscr1,luscr2,1,lblk)
        call vecsmd(vec1,vec2,1d0,1d0,luscr2,lugrvf,lures,1,lblk)
        call relunit(luscr2,'delete')
      else
        call vecsmd(vec1,vec2,1d0,1d0,luscr1,lugrvf,lures,1,lblk)
      end if

      xresnrm = sqrt(inprdd(vec1,vec1,lures,lures,1,lblk))

      call relunit(luscr1,'delete')

      return
      end

      subroutine optc_prjout(nrdvec,lurdvec,luvec,
     &     vec1,vec2,nwfpar,lincore)

      implicit none

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     lincore
      integer, intent(in) ::
     &     nrdvec, lurdvec, luvec, nwfpar
      real(8), intent(inout) ::
     &     vec1(*), vec2(*)
      
      integer ::
     &     irdvec, luscr1, luscr2
      real(8) ::
     &     ovl(nrdvec), xnrm
      real(8), external ::
     &     inprod, inprdd
      integer, external ::
     &     iopen_nus

      if (lincore) then
        call vec_from_disc(vec1,nwfpar,1,-1,luvec)
      end if
      if (ntest.ge.100) then
        if (lincore) then
          xnrm = sqrt(inprod(vec1,vec1,nwfpar))
        else
          xnrm = sqrt(inprdd(vec1,vec1,luvec,luvec,1,-1))
        end if
        write(6,*) ' norm of unprojected vector: ',xnrm
      end if

      call rewino(lurdvec)
      do irdvec = 1, nrdvec
        if (lincore) then
          call vec_from_disc(vec2,nwfpar,0,-1,lurdvec)
          ovl(irdvec) = inprod(vec1,vec2,nwfpar)
        else
          call rewino(luvec)
          ovl(irdvec) = inprdd(vec1,vec2,luvec,lurdvec,0,-1)
        end if
        if (ntest.ge.100)
     &       write(6,*) ' overlap with vec ',irdvec,' :',ovl(irdvec)
        if (lincore)
     &       vec1(1:nwfpar) =
     &       vec1(1:nwfpar)-ovl(irdvec)*vec2(1:nwfpar)
      end do

      if (lincore) then
        call vec_to_disc(vec1,nwfpar,1,-1,luvec)
      else
        luscr1 = iopen_nus('PRJOUT_SCR1')
        luscr2 = iopen_nus('PRJOUT_SCR2')

        call mvcsmd(lurdvec,ovl,luscr1,luscr2,vec1,vec2,nrdvec,1,-1)
        call vecsmd(vec1,vec2,1d0,-1d0,luvec,luscr1,luscr2,1,-1)
        call copvcd(luscr2,luvec,vec1,1,-1)

        call relunit(luscr1,'delete')
        call relunit(luscr2,'delete')
      end if

      if (ntest.ge.100) then
        if (lincore) then
          xnrm = sqrt(inprod(vec1,vec1,nwfpar))
        else
          xnrm = sqrt(inprdd(vec1,vec1,luvec,luvec,1,-1))
        end if
        write(6,*) ' norm of projected vector:   ',xnrm

      end if

      return
      end
c $Id$
