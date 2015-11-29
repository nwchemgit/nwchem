C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Perform a linesearch for the minimum gradient with the state
C> and search direction specified as described in `wfn1_next_step.dox`
C>
C> This subroutine is in many ways similar to `wfn1_linesearch`.
C> The main differences are that:
C>
C> 1. The wavefunction and search direction are stored differently.
C>
C> 2. Different routines are invoked to compute the energies and
C>    gradients
C>
C> 3. This routine seeks to minimize the gradient instead of the
C>    energy
C>
C> This routine performs a line search the find a minimum of the
C> gradient. In this implementation a second order approximation is
C> used in that a parabola is fitted to 3 data elements. The minimum
C> of the parabola is taken as a guess for the true minimum.
C> This general idea is repeated in a cycle until the gradient of the
C> parabola at the last point is small enough.
C>
C> The algorithm behind this uses a start up phase and an iterative
C> process. In both phases a parabola is fitted to available data and
C> its minimum established. In the start up phase a parabola is
C> fitted to f(0), df(0) and f(x1). In the iterative phases a different
C> approach based on known function values is used. I.e. the input
C> data is f(x1), f(x2), and f(x3).
C>
C> In short the start up phase works as follows:
C> - Given is f(0) and df(0), calc |df(0)|
C> - x1 = df(0), calc f(x1), df(x1) and |df(x1)|
C> - x2 = -df(0), calc f(x2), df(x2) and |df(x2)|
C> - fit parabola to |df(x1)|, |df(0)|, and |df(x2)|, and find its
C>   minimum xm
C> - Calculate f(xm), df(xm) and |df(xm)|
C> - Add (0,|df(0)|), (x1,|df(x1)|), (x2,|df(x2)|) and (xm,|df(xm)|) to
C>   an ordered table
C>
C> The iterative process then proceeds:
C> - Find the point xg with the smallest gradient df(xg)
C> - Find two points closest to xg in position
C> - Call the sequences of these three points x1, x2, and x3, along
C>   with associated values |df(x1)|, |df(x2)|, and |df(x3)|
C> - Fit a parabola to the three data points and establish the
C>   corresponding minum xm
C> - Calculate f(xm), df(xm) and |df(xm)|
C> - Add (xm,|df(xm)|) to the table
C>
      subroutine wfn1_nxt2_linesearch_g(rtdb,geom,nbf,nea,neb,step,
     &           h1,eri,erix,ov,ov12,ovm12,state,fock,grad,etot,tol)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfn1_nxt.fh"
c
      integer rtdb !< [Input] The RTDB handle
      integer geom !< [Input] The geometry handle
      integer nbf  !< [Input] The number of basis functions
      integer nea  !< [Input] The number of alpha electrons
      integer neb  !< [Input] The number of beta electrons
c
      double precision h1(nbf,nbf) !< [Input] The 1-electron Hamiltonian
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                            !< Coulomb integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                             !< eXchange integrals
      double precision ov(nbf,nbf) !< [Input] The overlap integrals
      double precision ov12(nbf,nbf) !< [Input] The sqrt(overlap)
      double precision ovm12(nbf,nbf) !< [Input] The 1/sqrt(overlap)
c
      double precision state(8*nbf*nbf) !< [Input] The state vector
      double precision fock(8*nbf*nbf)  !< [Input] The Fock matrices
      double precision grad(8*nbf*nbf)  !< [Input] The gradient
c
      double precision step !< [Output] The step to the minimum
      double precision etot !< [In/Output] The total energy at the minimum
      double precision tol  !< [Input] The tolerance on the gradient
c
      integer mxpts
      parameter(mxpts = 100)
      integer          npts      !< The number of points
      double precision xx(mxpts) !< The coordinate
      double precision fx(mxpts) !< The function value fx(i) = f(xx(i))
      double precision trst      !< The trust region
      double precision xm        !< The coordinate of the predicted 
                                 !< minimum
      double precision fxm       !< The predicted minimum value
      double precision a         !< Polynomial coefficient
      double precision b         !< Polynomial coefficient
      double precision c         !< Polynomial coefficient
      integer          ix1       !< The point before the minimum
      integer          ix2       !< The point at the minimum
      integer          ix3       !< The point after the minimum
      integer          indxmin   !< The minimal index
      double precision df        !< Estimate of gradient at a point
      double precision dnrm      !< The norm
      double precision dnrm2     !< Another copy of the norm
c
      double precision en0, e0    !< Energies
      double precision levelshift !< A levelshifter (not used)
c
      integer l_lstate, k_lstate !< Local state vector
      integer l_lgrad,  k_lgrad  !< Local gradient vector
      integer l_lfock,  k_lfock  !< Local Fock vector
c
      integer          wfn1_idmin
      external         wfn1_idmin
      double precision wfn1_nxt_dot, ddot
      external         wfn1_nxt_dot, ddot
c
      levelshift = -1.0d0
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"lstate",
     +                     l_lstate,k_lstate))
     +   call errquit("wfn1_nxt_linesearch: could not allocate lstate",
     +   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"lgrad",
     +                     l_lgrad,k_lgrad))
     +   call errquit("wfn1_nxt_linesearch: could not allocate lgrad",
     +   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,wfn1_nxt_size1(),"lfock",
     +                     l_lfock,k_lfock))
     +   call errquit("wfn1_nxt_linesearch: could not allocate lfock",
     +   ma_sizeof(MT_DBL,wfn1_nxt_size1(),MT_BYTE),MA_ERR)
c
      dnrm = 
     +   ddot(nbf*nbf,grad(wfn1_nxt_anat()),1,grad(wfn1_nxt_anat()),1)
     +  +ddot(nbf*nbf,grad(wfn1_nxt_bnat()),1,grad(wfn1_nxt_bnat()),1)
     +  +ddot(nbf*nbf,grad(wfn1_nxt_aocc()),1,grad(wfn1_nxt_aocc()),1)
     +  +ddot(nbf*nbf,grad(wfn1_nxt_bocc()),1,grad(wfn1_nxt_bocc()),1)
      dnrm = sqrt(dnrm)
cDEBUG
      write(*,*)"wfn1_nxt2_linesearch_g: step: ",0.0d0," grad: ",dnrm
cDEBUG
c
      trst     = 1.0d0
      npts     = 1
      xx(npts) = 0.0d0
      fx(npts) = dnrm
c
      call dcopy(wfn1_nxt_size1(),grad,1,dbl_mb(k_lgrad),1)
      call dscal(wfn1_nxt_size1(),1.0d0,dbl_mb(k_lgrad),1)
      call dcopy(wfn1_nxt_size1(),fock,1,dbl_mb(k_lfock),1)
      call wfn1_nxt2_step_g(nea,neb,nbf,1,state,dbl_mb(k_lfock),
     +     dbl_mb(k_lgrad),dbl_mb(k_lstate))
cDEBUG
c     write(*,*)'step: ',1.0d0
c     call wfn1_nxt_print_state(dbl_mb(k_lstate),nbf)
cDEBUG
      call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,h1,eri,erix,
     +     ov,ov12,ovm12,dbl_mb(k_lstate),dbl_mb(k_lfock),
     +     dbl_mb(k_lgrad),en0,levelshift)
      k_lgrad = k_lgrad-1
      dnrm2 = 
     +   ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_anat()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_anat()),1)
     +  +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_bnat()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_bnat()),1)
     +  +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_aocc()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_aocc()),1)
     +  +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_bocc()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_bocc()),1)
      dnrm2 = sqrt(dnrm2)
      k_lgrad = k_lgrad+1
cDEBUG
      write(*,*)"wfn1_nxt2_linesearch_g: step: ",1.0d0," grad: ",dnrm2
cDEBUG
      call wfn1_insert(mxpts,npts,xx,fx,dnrm,dnrm2)
c
      call dcopy(wfn1_nxt_size1(),grad,1,dbl_mb(k_lgrad),1)
      call dscal(wfn1_nxt_size1(),-1.0d0,dbl_mb(k_lgrad),1)
      call dcopy(wfn1_nxt_size1(),fock,1,dbl_mb(k_lfock),1)
      call wfn1_nxt2_step_g(nea,neb,nbf,1,state,dbl_mb(k_lfock),
     +     dbl_mb(k_lgrad),dbl_mb(k_lstate))
cDEBUG
c     write(*,*)'step: ',-1.0d0
c     call wfn1_nxt_print_state(dbl_mb(k_lstate),nbf)
cDEBUG
      call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,h1,eri,erix,
     +     ov,ov12,ovm12,dbl_mb(k_lstate),dbl_mb(k_lfock),
     +     dbl_mb(k_lgrad),en0,levelshift)
      k_lgrad = k_lgrad-1
      dnrm2 = 
     +   ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_anat()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_anat()),1)
     +  +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_bnat()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_bnat()),1)
     +  +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_aocc()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_aocc()),1)
     +  +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_bocc()),1,
     +                dbl_mb(k_lgrad+wfn1_nxt_bocc()),1)
      dnrm2 = sqrt(dnrm2)
      k_lgrad = k_lgrad+1
cDEBUG
      write(*,*)"wfn1_nxt2_linesearch_g: step: ",-1.0d0," grad: ",dnrm2
cDEBUG
      call wfn1_insert(mxpts,npts,xx,fx,-dnrm,dnrm2)
c
      indxmin = wfn1_idmin(npts,fx,1)
      call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
      call wfn1_f1f2f3(xx(ix1),xx(ix2),xx(ix3),
     +                 fx(ix1),fx(ix2),fx(ix3),trst,a,b,c,xm,fxm)
      df = 2.0d0*a*xx(indxmin)+b
      step = xm/dnrm
cDEBUG
      write(*,*)"*** f1f2f3 a: a,b,c,xx,df= ",a,b,c,xx(indxmin),df
c     write(*,*)"*** f0df0f1 a: etot,en0,a,b,c= ",etot,en0,a,b,c
c     write(*,*)"*** f0df0f1 b: dnrm,xm,fxm,df= ",dnrm,xm,fxm,df
cDEBUG
c
      do while (abs(df).gt.a*tol.and.abs(xm).ge.1.0d-16.and.
     +          npts.lt.mxpts)
        call dcopy(wfn1_nxt_size1(),grad,1,dbl_mb(k_lgrad),1)
        call dscal(wfn1_nxt_size1(),xm/dnrm,dbl_mb(k_lgrad),1)
        call dcopy(wfn1_nxt_size1(),fock,1,dbl_mb(k_lfock),1)
        call wfn1_nxt2_step_g(nea,neb,nbf,1,state,dbl_mb(k_lfock),
     +       dbl_mb(k_lgrad),dbl_mb(k_lstate))
        call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,h1,
     +       eri,erix,ov,ov12,ovm12,dbl_mb(k_lstate),dbl_mb(k_lfock),
     +       dbl_mb(k_lgrad),en0,levelshift)
        dnrm2 = 
     +     ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_anat()),1,
     +                  dbl_mb(k_lgrad+wfn1_nxt_anat()),1)
     +    +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_bnat()),1,
     +                  dbl_mb(k_lgrad+wfn1_nxt_bnat()),1)
     +    +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_aocc()),1,
     +                  dbl_mb(k_lgrad+wfn1_nxt_aocc()),1)
     +    +ddot(nbf*nbf,dbl_mb(k_lgrad+wfn1_nxt_bocc()),1,
     +                  dbl_mb(k_lgrad+wfn1_nxt_bocc()),1)
        dnrm2 = sqrt(dnrm2)
cDEBUG
        write(*,*)"wfn1_nxt2_linesearch_g: step: ",xm/dnrm,
     +                                   " grad: ",dnrm2,df
cDEBUG
        call wfn1_insert(mxpts,npts,xx,fx,xm,dnrm2)
        indxmin = wfn1_idmin(npts,fx,1)
        call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
        call wfn1_f1f2f3(xx(ix1),xx(ix2),xx(ix3),
     +                   fx(ix1),fx(ix2),fx(ix3),trst,a,b,c,xm,fxm)
        df = 2.0d0*a*xx(indxmin)+b
        step = xm/dnrm
cDEBUG
c       write(*,'("*** ix1,ix2,ix3",3i22)')ix1,ix2,ix3
c       write(*,'("*** xx1,xx2,xx3",3f22.14)')xx(ix1),xx(ix2),xx(ix3)
c       write(*,'("*** fx1,fx2,fx3",3f22.14)')fx(ix1),fx(ix2),fx(ix3)
c       write(*,*)"*** f1f2f3 b: a,b,c,xx,df= ",a,b,c,xx(indxmin),df
c       write(*,*)"*** npts,en0,xm,fxm,df=",npts,en0,xm,fxm,df
cDEBUG
      enddo
cDEBUG
      write(*,*)"*** df,xm,npts, step: ",df,xm,npts,step
      write(*,*)
      write(*,*)
c     step = max(step,1.0d-12)
cDEBUG
c
      if (.not.ma_pop_stack(l_lfock))
     +  call errquit("wfn1_nxt_linesearch: could not deallocate lfock",
     +       0,MA_ERR)
      if (.not.ma_pop_stack(l_lgrad))
     +  call errquit("wfn1_nxt_linesearch: could not deallocate lgrad",
     +       0,MA_ERR)
      if (.not.ma_pop_stack(l_lstate))
     +  call errquit("wfn1_nxt_linesearch: could not deallocate lstate",
     +       0,MA_ERR)
c
      end
C>
C> @}
