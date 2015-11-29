C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Perform a linesearch for the minimum energy with the state
C> and search direction specified as described in `wfn1_next_step.dox`
C>
C> This subroutine is in many ways similar to `wfn1_linesearch`.
C> The main differences are that:
C>
C> 1. The wavefunction and search direction are stored differently.
C>
C> 2. Different routines are invoked to compute the energies.
C>
C> This routine performs a line search the find a minimum of the
C> energy. In this implementation a second order approximation is
C> used in that a parabola is fitted to 3 data elements. The minimum
C> of the parabola is taken as a guess for the true minimum.
C> This general idea is repeated in a cycle until the gradient at the
C> last point is small enough.
C>
C> The algorithm behind this uses a start up phase and an iterative
C> process. In both phases a parabola is fitted to available data and
C> its minimum established. In the start up phase a parabola is
C> fitted to f(0), df(0) and f(x1). In the iterative phases a different
C> approach based on known function values is used. I.e. the input
C> data is f(x1), f(x2), and f(x3).
C>
C> In short the start up phase works as follows:
C> - Given is f(0) and df(0)
C> - x1 = -df(0), calc f(x1)
C> - fit parabola to f(0), df(0), and f(x1), and find its minimum xm
C> - Calculate f(xm)
C> - Add (0,f(0)), (x1,f(x1)), and (xm,f(xm)) to an ordered table
C>
C> The iterative process then proceeds:
C> - Find the point xg with the lowest energy f(xg)
C> - Find two points closest to xg in position
C> - Call the sequences of these three points x1, x2, and x3, along
C>   with associated values f(x1), f(x2), and f(x3)
C> - Fit a parabola to the three data points and establish the
C>   corresponding minum xm
C> - Calculate f(xm)
C> - Add (xm,f(xm)) to the table
C>
      subroutine wfn1_nxt2_linesearch_e(rtdb,geom,nbf,nea,neb,step,
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
      double precision state(8*nbf*nbf+2*nbf) !< [Input] The state vector
      double precision fock(8*nbf*nbf+2*nbf)  !< [Input] The Fock matrices
      double precision grad(8*nbf*nbf+2*nbf)  !< [Input] The gradient
c
      double precision step !< [Output] The step to the minimum
      double precision etot !< [In/Output] The total energy at the minimum
      double precision tol  !< [Input] The tolerance on the gradient
c
      integer mxpts
      parameter(mxpts = 200)
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
c
      double precision en0, e0    !< Energies
      double precision levelshift !< Levelshift (not actually used)
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
      trst     = 1.0d0
c
c     The next three lines are essentially equivalent to a wfn1_insert
c     call.
c
      npts     = 1
      xx(npts) = 0.0d0
      fx(npts) = etot
c
      call dcopy(wfn1_nxt_size1(),grad,1,dbl_mb(k_lgrad),1)
      call dscal(wfn1_nxt_size1(),-1.0d0,dbl_mb(k_lgrad),1)
      call dcopy(wfn1_nxt_size1(),fock,1,dbl_mb(k_lfock),1)
      call wfn1_nxt2_step_ls(nea,neb,nbf,100,state,dbl_mb(k_lfock),
     +     dbl_mb(k_lgrad),dbl_mb(k_lstate),levelshift)
      dnrm = sqrt(ddot(wfn1_nxt_size1(),grad,1,grad,1))
c
      call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,h1,eri,erix,
     +     ov,ov12,ovm12,dbl_mb(k_lstate),dbl_mb(k_lfock),
     +     dbl_mb(k_lgrad),en0,levelshift)
      call wfn1_insert(mxpts,npts,xx,fx,dnrm,en0)
      call wfn1_f0df0f1(dnrm,etot,-dnrm,en0,trst,a,b,c,xm,fxm)
      df = 2.0d0*a*dnrm+b
      step = 1.0d0
cDEBUG
c     write(*,*)"*** f0df0f1 a: etot,en0,a,b,c= ",etot,en0,a,b,c
c     write(*,*)"*** f0df0f1 b: dnrm,xm,fxm,df= ",dnrm,xm,fxm,df
c     write(*,*)"*** linesearch tol = ",tol
cDEBUG
c
      do while (abs(df).gt.1.0d-12.and.xm.ge.1.0d-15.and.npts.lt.mxpts)
        call dcopy(wfn1_nxt_size1(),grad,1,dbl_mb(k_lgrad),1)
        call dscal(wfn1_nxt_size1(),-xm/dnrm,dbl_mb(k_lgrad),1)
        call wfn1_nxt2_step_ls(nea,neb,nbf,100,state,dbl_mb(k_lfock),
     +       dbl_mb(k_lgrad),dbl_mb(k_lstate),levelshift)
        call wfn1_nxt2_energy_gradient(rtdb,geom,nbf,nea,neb,h1,
     +       eri,erix,ov,ov12,ovm12,dbl_mb(k_lstate),dbl_mb(k_lfock),
     +       dbl_mb(k_lgrad),en0,levelshift)
        call wfn1_insert(mxpts,npts,xx,fx,xm,en0)
        indxmin = wfn1_idmin(npts,fx,1)
        call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
        if (indxmin.eq.1) then
          call wfn1_f0df0f1(xx(ix2),etot,-dnrm,fx(ix2),trst,a,b,c,
     +                      xm,fxm)
        else
          call wfn1_f1f2f3(xx(ix1),xx(ix2),xx(ix3),
     +                     fx(ix1),fx(ix2),fx(ix3),trst,a,b,c,xm,fxm)
        endif
        df = 2.0d0*a*xx(indxmin)+b
        step = xx(indxmin)/dnrm
cDEBUG
c       write(*,*)"*** npts,en0,xm,fxm,df=",npts,en0,xm,fxm,df
cDEBUG
      enddo
cDEBUG
c     write(*,*)"*** npts, step: ",npts,step
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
