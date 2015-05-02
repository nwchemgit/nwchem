      Subroutine xc_hcth(tol_rho, xfac, lxfac, nlxfac,
     ,                    cfac, lcfac, nlcfac,rho, delrho, 
     &                    Amat, Cmat, nq, ipol, Ex, Ec,  qwght,
     &                    ldew,func,funcname)     
c
c$Id$
c
      Implicit none
c
#include "dft2drv.fh"
c
      logical ldew ! [input]
      logical lcfac, nlcfac,  lxfac, nlxfac ! [input]
      double precision func(*) ![input/output]
      double precision cfac, xfac ![input]
      character*4 funcname ! functional name [input]
c
      integer ipol  ! no. of spin states [input]
      integer nq    ! no. of quadrature pts [input]
      double precision tol_rho! [input]!threshold on density
      double precision Ec ! Correlation energy [input/output] 
      double precision Ex ! Exchange    energy [input/output] 
      double precision rho(nq,ipol*(ipol+1)/2)! Charge Density [input] 
      double precision delrho(nq,3,ipol) ! Charge Density Gradient[input] 
      double precision qwght(nq) ! Quadrature Weights [input]
      double precision Amat(nq,ipol)  !Sampling Matrices for the XC [output]
      double precision Cmat(nq,*)!Potential & Energy [output]
c     
c References:
c    F.A.Hamprecht, A.J.Cohen, D.J.Tozer and N.C.Handy, 
c    J. Chem. Phys. 109, 6264-6271 (1998)
c    subroutine supplied by Fred Hamprecht fah@igc.phys.chem.ethz.ch
c
      integer n
      double precision gammaval
c to hcth
      double precision rhoa 
      double precision rhob 
      double precision za   
      double precision zb   
      double precision hE_x, hE_c 
      double precision dfdrax, dfdrac
      double precision dfdrbx, dfdrbc
      double precision dfdzac,dfdzax,dfdza
      double precision dfdzbc,dfdzbx,dfdzb
c
      if(ipol.eq.1) then
        do n=1,nq
          if(rho(n,1).gt.tol_rho) then 
            rhoa=0.d0
            rhoa=0.5d0*rho(n,1)
            rhob=0.d0
            rhob=rhoa
            za=0.d0
            gammaval=0.25d0*(delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1))
            if(gammaval.gt.tol_rho**2)  za=sqrt(gammaval)
            zb=za
            call hcth(ipol,funcname,
     *           dfdrax,dfdrac, dfdzax,dfdzac, 
     *           dfdrbx,dfdrbc, dfdzbx,dfdzbc, 
     .           rhoa, rhob, 
     &           za, zb, hE_x, hE_c, tol_rho)
            if(ldew) func(n)=func(n)+hE_c*cfac+hE_x*xfac
            Ec=Ec+hE_c*qwght(n)*cfac
            Ex=Ex+hE_x*qwght(n)*xfac
            Amat(n,1) = Amat(n,1)+dfdrac*cfac+dfdrax*xfac
            dfdza=(dfdzac*cfac+dfdzax*xfac)*0.5d0 
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dfdza
          endif
        enddo
      else
        do n=1,nq
          if(rho(n,1).gt.tol_rho) then
          rhoa=rho(n,2)
          rhob=rho(n,3)
            za=0.d0
            gammaval=delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
            if(gammaval.gt.tol_rho**2) za=sqrt(gammaval)
            zb=0.d0
            gammaval=delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &           delrho(n,3,2)*delrho(n,3,2) 
            if(gammaval.gt.tol_rho**2) zb=sqrt(gammaval)
            call hcth(ipol,funcname,
     *           dfdrax,dfdrac, dfdzax,dfdzac, 
     *           dfdrbx,dfdrbc, dfdzbx,dfdzbc, 
     .           rhoa, rhob, 
     &           za, zb, hE_x, hE_c, tol_rho)
            if(ldew) func(n)=func(n)+hE_c*cfac+hE_x*xfac 
            Ec=Ec+hE_c*qwght(n)*cfac
            Ex=Ex+hE_x*qwght(n)*xfac
            Amat(n,1) = Amat(n,1)+dfdrac*cfac+dfdrax*xfac
            Amat(n,2) = Amat(n,2)+dfdrbc*cfac+dfdrbx*xfac
            dfdza=dfdzac*cfac+dfdzax*xfac
            dfdzb=dfdzbc*cfac+dfdzbx*xfac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dfdza*0.5d0
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dfdzb*0.5d0
          endif
        enddo
      endif
      return
      end

      Subroutine xc_hcth_d2(tol_rho, xfac, lxfac, nlxfac,
     &                      cfac, lcfac, nlcfac,rho, delrho,
     &                      Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &                      Ex, Ec,  qwght, ldew,func,funcname)
c
c$Id$
c
      Implicit none
#include "errquit.fh"
c
#include "dft2drv.fh"
#include "2ndDerivB97.h"
#include "mafdecls.fh"
c
      logical ldew ! [input]
      logical lcfac, nlcfac,  lxfac, nlxfac ! [input]
      double precision func(*) ![input/output]
      double precision cfac, xfac ![input]
      character*4 funcname ! functional name [input]
c
      integer ipol  ! no. of spin states [input]
      integer nq    ! no. of quadrature pts [input]
      double precision tol_rho! [input]!threshold on density
      double precision Ec ! Correlation energy [input/output] 
      double precision Ex ! Exchange    energy [input/output] 
      double precision rho(nq,ipol*(ipol+1)/2)! Charge Density [input] 
      double precision delrho(nq,3,ipol) ! Charge Density Gradient[input] 
      double precision qwght(nq) ! Quadrature Weights [input]
      double precision Amat(nq,ipol)  !Sampling Matrices for the XC [output]
      double precision Cmat(nq,*)!Potential & Energy [output]
      double precision Amat2(nq,NCOL_AMAT2) ! XC functional seconds [output]
      double precision Cmat2(nq,NCOL_CMAT2) ! XC functional seconds [output]
c
      integer i,max_pow_u
      double precision rho_a(0:3),rho_b(0:3)
      double precision FX(0:_FXC_NUMDERI),FC(0:_FXC_NUMDERI)
#include "xc_hcth.fh"
      double precision sol((limpow+1)*3)
      call xc_htch_loadcf(funcname,sol,max_pow_u)
c      if(max_pow_u.gt.2) call errquit(' not ready for maxpow=4',0,0)
      do i=1,nq
         if(ipol.eq.2) then
            rho_a(0)=rho(i,2)
            rho_b(0)=rho(i,3)
            rho_a(1)=delrho(i,1,1)**2
            rho_a(2)=delrho(i,2,1)**2
            rho_a(3)=delrho(i,3,1)**2
            rho_b(1)=delrho(i,1,2)**2
            rho_b(2)=delrho(i,2,2)**2
            rho_b(3)=delrho(i,3,2)**2
         else
            rho_a(0)=rho(i,1)*.5d0
            rho_b(0)=rho_a(0)
            rho_a(1)=(delrho(i,1,1)*.5d0)**2
            rho_a(2)=(delrho(i,2,1)*.5d0)**2
            rho_a(3)=(delrho(i,3,1)*.5d0)**2
            rho_b(1)=rho_a(1)
            rho_b(2)=rho_a(2)
            rho_b(3)=rho_a(3)
         endif
         if(rho_a(0).gt.tol_rho.or.rho_b(0).gt.tol_rho) then
            call dft_xckernel_xb97(rho_a, rho_b, 1d0, tol_rho, 
     F           FX, max_pow_u, sol)
            call dscal(_FXC_NUMDERI+1,xfac,FX,1)
            Ex = Ex + FX(_FXC_E)*qwght(i)
            if(ldew) func(i)=func(i) + FX(_FXC_E)
            call xc_fxc2acmat(nq,i,ipol,FX,
     A           Amat,Cmat,Amat2,Cmat2)
            call dft_xckernel_cb97(rho_a, rho_b, 1d0, tol_rho, 
     F           FC, max_pow_u, sol)
            call dscal(_FXC_NUMDERI+1,cfac,FC,1)
            Ec = Ec + FC(_FXC_E)*qwght(i)
            if(ldew) func(i)=func(i) + FC(_FXC_E)
            call xc_fxc2acmat(nq,i,ipol,FC,
     A           Amat,Cmat,Amat2,Cmat2)
         endif
      enddo

      return
      end
      Subroutine xc_hcth_d2num(tol_rho, xfac, lxfac, nlxfac,
     &                      cfac, lcfac, nlcfac,rho, delrho,
     &                      Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &                      Ex, Ec,  qwght, ldew,func,funcname)
c
c$Id$
c
      Implicit none
#include "errquit.fh"
c
#include "dft2drv.fh"
#include "mafdecls.fh"
c
      logical ldew ! [input]
      logical lcfac, nlcfac,  lxfac, nlxfac ! [input]
      double precision func(*) ![input/output]
      double precision cfac, xfac ![input]
      character*4 funcname ! functional name [input]
c
      integer ipol  ! no. of spin states [input]
      integer nq    ! no. of quadrature pts [input]
      double precision tol_rho! [input]!threshold on density
      double precision Ec ! Correlation energy [input/output] 
      double precision Ex ! Exchange    energy [input/output] 
      double precision rho(nq,ipol*(ipol+1)/2)! Charge Density [input] 
      double precision delrho(nq,3,ipol) ! Charge Density Gradient[input] 
      double precision qwght(nq) ! Quadrature Weights [input]
      double precision Amat(nq,ipol)  !Sampling Matrices for the XC [output]
      double precision Cmat(nq,*)!Potential & Energy [output]
      double precision Amat2(nq,NCOL_AMAT2) ! XC functional seconds [output]
      double precision Cmat2(nq,NCOL_CMAT2) ! XC functional seconds [output]
c
c     Local variables
c
      integer l_storage, i_prho, i_pdelrho, i_pAmat, i_pCmat, i_pfunc,
     &     i_qwght_copy, npert
      double precision ExDum, EcDum
c
c     First get the functional and first derivative values
c
      call xc_hcth(tol_rho, xfac, lxfac, nlxfac, cfac, lcfac, nlcfac,
     &     rho, delrho, Amat, Cmat, nq, ipol, Ex, Ec, qwght, ldew, func,
     &     funcname)
c
c     Compute the second derivative values by finite difference
c
      call xc_setup_fd(tol_rho, rho, delrho, qwght, nq, ipol, .true.,
     &     l_storage, i_prho, i_pdelrho, i_pAmat, i_pCmat, i_pfunc,
     &     i_qwght_copy)
c
c     Compute functional first derivatives at perturbed density parameter
c     values - note that the number of points is nq*2*npert and that the
c     routine is called as unrestricted
c
      npert = 5
      call xc_hcth(tol_rho, xfac, lxfac, nlxfac, cfac, lcfac, nlcfac,
     &     dbl_mb(i_prho), dbl_mb(i_pdelrho), dbl_mb(i_pAmat),
     &     dbl_mb(i_pCmat), nq*2*npert, ipol, ExDum, EcDum,
     &     dbl_mb(i_qwght_copy), ldew, dbl_mb(i_pfunc), funcname)
      call xc_make_fd(Amat2, Cmat2, nq, .true., dbl_mb(i_pAmat),
     &     dbl_mb(i_pCmat))
c
c     Free temporary storage allocated by xc_setup_fd
c
      if (.not.ma_free_heap(l_storage))
     &     call errquit('xc_hcth_d2: cannot pop stack',0, MA_ERR)
c
      return
      end

      SUBROUTINE hcth(ipol,functional,
     *     dfdrax,dfdrac, dfdzax,dfdzac, 
     *     dfdrbx,dfdrbc, dfdzbx,dfdzbc, 
     1     rhoa, rhob, za, zb, hE_x, hE_c, 
     2     tol_rho)

c    subroutine supplied by Fred Hamprecht fah@igc.phys.chem.ethz.ch
C    SUPPLIED TO THE ROUTINE:
C    
C    rhoa   -- value of rhoalpha at a given grid point 
C    rhob   -- value of rhobeta at a given grid point
C    za     -- zeta_alpha, as defined in the TH1 paper (JCP 108 2545), 
C              that is mod(grad(rhoalpha)), a scalar quantity.
C    zb     -- mod(grad(rhobeta)) 
C    zab    -- zeta_{alpha beta} as defined in the TH1 paper, that is
C              grad(rhoalpha).grad(rhobeta) 
C    energy -- a boolean variable deciding whether to compute the energy 
C              contribution at the point in space (true) or the
C              appropriate derivatives (false) needed for the KS matrix
C              _and_ the energy contribution.

C    RETURNED FROM THE ROUTINE:

C    hE_x --   the contribution to the energy at this point in space.
C    hE_c --   the contribution to the energy at this point in space.
C    dfdra  -- partial functional derivative of F_xc with respect to 
C              rhoalpha
C    dfdrb  -- partial functional derivative of F_xc with respect to   
C              rhobeta
C    dfdza  -- partial functional derivative of F_xc with respect to   
C              mod(grad(rhoalpha)), divided by za !!!!!!!!!
C              i.e.  1    d f 
C                   --- * ---- 
C                   za    d za
C              This is a consequence of the Cadpac implementation
C    dfdzb  -- partial functional derivative of F_xc with respect to   
C              mod(grad(rhobeta)), divided by zb !!!!!!!!!


      implicit none
#include "errquit.fh"
      
      double precision rhoa ![input]
      double precision rhob ![input]
      double precision za   ![input]
      double precision zb   ![input]
      integer ipol ![input]
      double precision hE_x ![output]
      double precision hE_c ![output]
      double precision tol_rho ![input]
      double precision dfdrax    ![output]
      double precision dfdrac    ![output]
      double precision dfdrbx    ![output]
      double precision dfdrbc    ![output]
      double precision dfdzax    ![output]
      double precision dfdzac    ![output]
      double precision dfdzbx    ![output]
      double precision dfdzbc    ![output]
      character*4 functional
      double precision pi
      PARAMETER (PI=3.1415926535898D0)
#include "xc_hcth.fh"
c
c     variables passed to hcderiv
c
      integer nofunc,max_pow_u
      double precision sol((limpow+1)*3), F((limpow+1)*3,4), 
C     &          FF((limpow+1)*3,5,4),
     &          F_xc((limpow+1)*3)

Cfah sol -- contains the coefficients of the terms in F_xc
Cfah        convention: sol(1) = c_{x alpha, 0}, c_{x beta, 0}
Cfah                    sol(2) = c_{c alpha alpha, 0}, c_{c beta beta, 0} 
Cfah                    sol(3) = c_{c alpha beta, 0} 
Cfah                    sol(4) = c_{x alpha, 1}, c_{x beta, 1}
Cfah                    sol(5) = c_{c alpha alpha, 1}, c_{c beta beta, 1} 
Cfah                    sol(6) = c_{c alpha beta, 1} 
Cfah                           
Cfah                           etc.
Cfah 
Cfah f(5) -- contains the partial first functional derivatives of F_xc with 
Cfah respect to 
Cfah the four quantities (IN THIS ORDER): ra, rb, za, zb
Cfah 
Cfah ff(5,5) contains the second derivatives with
Cfah respect to the same five quantities

Cfah F_xa -- contains the alpha exchange bit containing the various powers 
Cfah         of u_{x alpha} (eq. (18) of Becke V paper) 
Cfah F_xb --              beta       
Cfah            u_{x beta} 
Cfah F_caa -- contains the alpha parallel spin correlation bit with the powers
Cfah          of u_{c alpha alpha} 
Cfah F_cbb --              beta 
Cfah             u_{c beta beta} 
Cfah F_cab -- contains the anti-parallel spin correlation bit with the powers 
Cfah          of u_{c alpha beta} 

C     Initialise

      dfdrac = 0.D0
      dfdrax = 0.D0
      dfdrbc = 0.D0
      dfdrbx = 0.D0
      dfdzac = 0.D0
      dfdzax = 0.D0
      dfdzbc = 0.D0
      dfdzbx = 0.D0
      hE_c = 0.D0
      hE_x = 0.D0

      IF (rhoa .LT. tol_rho.and.rhob.lt.tol_rho) RETURN
Cfah numerical cutoff: if the density is too low, its contribution is 
Cfah neglectable. 
      call xc_htch_loadcf(functional,sol,max_pow_u)
      
      CALL hcderiv(max_pow_u,ipol,
     &     F,
CFF,
     &     F_xc,
     &     rhoa,rhob,za,zb,
     &     sol,tol_rho)

c     if(ipol.eq.2) then
c       DO n = 1, (max_pow_u+1)*3 
c         dfdra = dfdra + F(n,1) 
c         dfdrb = dfdrb + F(n,2) 
c         if(za.gt.tol_rho) dfdza = dfdza + F(n,3) / za
c         if(zb.gt.tol_rho) dfdzb = dfdzb + F(n,4) / zb  
c       ENDDO
c     else
c        DO n = 1, (max_pow_u+1)*3 
c          dfdra = dfdra + F(n,1) 
c        enddo
c        if(za.gt.tol_rho) then
c          DO n = 1, (max_pow_u+1)*3 
c            dfdza = dfdza + F(n,3) / za
c          enddo
c        endif
c     endif
Cfah big thanks to NCH: cadpac requires df/(za * dza), NOT 
Cfah                                    df/dza 
cDEC$ NOVECTOR
      DO n = 0, max_pow_u 
        hE_x = hE_x + F_xc (n*3 + 1) 
        hE_c = hE_c + F_xc (n*3 + 2) + F_xc (n*3 + 3)
        dfdrax = dfdrax + F(n*3+1,1) 
        dfdrac = dfdrac + F(n*3+2,1) + F(n*3+3,1) 
        if(za.gt.tol_rho) then
            dfdzax = dfdzax + F(n*3+1,3) / za
            dfdzac = dfdzac + (F(n*3+2,3)+F(n*3+3,3)) / za
        endif
      if(ipol.eq.2) then
        dfdrbx = dfdrbx + F(n*3+1,2) 
        dfdrbc = dfdrbc + F(n*3+2,2) + F(n*3+3,2) 
        if(zb.gt.tol_rho) then
            dfdzbx = dfdzbx + F(n*3+1,4) / zb
            dfdzbc = dfdzbc + (F(n*3+2,4)+F(n*3+3,4)) / zb
        endif
      endif
      ENDDO
      RETURN
      END
      SUBROUTINE hcderiv(max_pow_u,ipol,
     &     F,
CFF, 
     &     F_xc,
     &             rhoa,rhob,za,zb,
     &             sol,tol_rho)

c    subroutine supplied by Fred Hamprecht fah@igc.phys.chem.ethz.ch
      implicit none
      INTEGER max_pow_u,ipol
      integer limpow
      parameter (limpow=4)
      double precision f_xc((limpow + 1)*3)
      double precision f((limpow+1)*3,4)
C, ff((limpow+1)*3,5,4)
      double precision rhoa, rhob, za, zb,tol_rho

Cfah  COMMON/special/h_atom
Cfah  LOGICAL h_atom
      DOUBLE PRECISION sol((limpow+1)*3)
    
      
      DOUBLE PRECISION dF_xa(4)
      DOUBLE PRECISION dF_xb(4)
      DOUBLE PRECISION dF_caa(4)
      DOUBLE PRECISION dF_cbb(4)
      DOUBLE PRECISION dF_cab(4)
Cfah these are the first derivatives of the terms of F_xc with respect to 
Cfah the 4 quantities. the index
Cfah runs over the particular partial derivatives of each term.  
Cfah More explicitly: these are the partial functional derivatives of 
Cfah F_XXX with respect to rhoa, rhob, za and zb. 

c      DOUBLE PRECISION d2F_xa(4,4)
c      DOUBLE PRECISION d2F_xb(4,4)
c      DOUBLE PRECISION d2F_caa(4,4)
c      DOUBLE PRECISION d2F_cbb(4,4)
c      DOUBLE PRECISION d2F_cab(4,4)

Cfah these are the first derivatives of the different transformed variables 
Cfah u with respect to rhoa, rhob, za and zb. These different derivatives 
Cfah with respect to these 4 quantities named above are stored in these 
Cfah arrays.

      DOUBLE PRECISION Pi 
      PARAMETER (Pi = 3.1415926535898D0)
      double precision rho
      DOUBLE PRECISION s_a2, s_b2, s_avg2, u_caa, u_cbb, u_cab
      DOUBLE PRECISION du_caa_by_drhoa, du_caa_by_dza, du_cbb_by_drhob 
      DOUBLE PRECISION du_cbb_by_dzb, du_cab_by_drhoa, du_cab_by_drhob 
      DOUBLE PRECISION du_cab_by_dza, du_cab_by_dzb
C, du_caa_by_drhoa_dza 
C      DOUBLE PRECISION du_caa_by_dza_dza, du_cbb_by_dzb_dzb
C      DOUBLE PRECISION du_cbb_by_drhob_dzb, du_cab_by_drhoa_dza 
C      DOUBLE PRECISION du_cab_by_drhoa_dzb, du_cab_by_drhob_dza 
C      DOUBLE PRECISION du_cab_by_drhob_dzb, du_cab_by_dza_dza, 
C     ,du_cab_by_dza_dzb 
C      DOUBLE PRECISION du_cab_by_dzb_dzb 
      DOUBLE PRECISION rsa, rsa12, rsa32, rsa21, rsb, 
     ,rsb12, rsb32, rsb21 
      DOUBLE PRECISION rsab, rsab12, rsab32, rsab21 
      DOUBLE PRECISION drsa_by_drhoa, drsb_by_drhob, drsab_by_drhoa
      DOUBLE PRECISION drsab_by_drhob 
      DOUBLE PRECISION zeta, dzeta_by_drhoa, dzeta_by_drhob 
      DOUBLE PRECISION fzeta, dfzeta_by_dzeta, 
     ,     e_crsa1, e_crsb1
      DOUBLE PRECISION e_crsab1, e_crsab0, a_crsab 
      DOUBLE PRECISION e_crsabzeta, de_crsa1_by_drsa, de_crsb1_by_drsb 
      DOUBLE PRECISION da_crsab_by_drsab, de_crsab0_by_drsab 
      DOUBLE PRECISION de_crsab1_by_drsab, de_crsabzeta_by_drsab 
      DOUBLE PRECISION de_crsabzeta_by_dzeta, e_caa, e_cbb, e_cab, 
     & de_caa_by_drhoa, de_cbb_by_drhob, de_cab_by_drhoa, 
     & de_cab_by_drhob,
     & c_naa, c_nbb, c_nab
      DOUBLE PRECISION F_xs,F_xs0 ! this is a function which is called. 
      DOUBLE PRECISION dF_xs_by_drhos, dF_xs_by_dzs,
     ,dF_xs_by_drhos0,dF_xs_by_drhos1,  dF_xs_by_dzs1
      INTEGER i, j,  n
      integer n1
      double precision x1,x2,x3,x4
      double precision e_crs1
      double precision drsbydrh
      double precision d2ez
      double precision decrsdrs
C
C     F_xs computes HCTH contribution to exchange Energy
C     using Dirac functional as LDA part
C     usage of F_xs
C     F_xs(n, sol(), rhoa, za)
C
      F_xs(n1, x1, x2, x3) =
     = (-3.D0*x1*(3.D0/Pi)**(1.D0/3.D0)*x2**(4.D0/3.D0)*
     -     ((0.004D0*x3**2.D0)/(0.004D0*x3**2.D0 +
     -     x2**(8.D0/3.D0)))**n1)/(2.D0*2.D0**(2.D0/3.D0))
      F_xs0( x1, x2, x3) =
     = (-3.D0*x1*(3.D0/Pi)**(1.D0/3.D0)*x2**(4.D0/3.D0)
     -     )/(2.D0*2.D0**(2.D0/3.D0))
C
C     dF_xs_by_drhos computes dE_x/drho derivative
Cfah  computes the derivative of the term with u^n of the exchange part of
Cfah  F_xc with respect to rho of the same spin.
Cfah  n     -- the power of u involved in this term
Cfah  c_xs  -- the coefficient c_xs(n) of the term of spin s with the
Cfah           power n of u; is NOT passed over as an array.
Cfah  rhos -- rhosigma, that is, either rhoalpha or rhobeta
Cfah  zs    -- mod(grad(rhosigma)), again for alpha or beta
C     usage dF_xs_by_drhos(n, c_xs, rhos, zs)
C
       dF_xs_by_drhos(n1, x1, x2, x3) =  -(x1*(6.D0/Pi)**
     *     (1.D0/3.D0)*x2**(1.D0/3.D0)*((0.004D0*x3*x3)/
     /     (x2**(8.D0/3.D0)+0.004D0*x3*x3))**n1)+
     +     (x1*0.008D0*n1*(6.D0/Pi)**(1.D0/3.D0)*
     *     x2**3.D0*x3*x3*((0.004D0*x3*x3)/
     /     (x2**(8.D0/3.D0)+0.004D0*x3*x3))**(-1 + n1))/
     /     (x2**(8.D0/3.D0)+0.004D0*x3*x3)**2
       dF_xs_by_drhos0(x1, x2, x3) =  -(x1*(6.D0/Pi)**
     *     (1.D0/3.D0)*x2**(1.D0/3.D0))
       dF_xs_by_drhos1(x1, x2, x3) =  -(x1*(6.D0/Pi)**
     *     (1.D0/3.D0)*x2**(1.D0/3.D0)*((0.004D0*x3*x3)/
     /     (x2**(8.D0/3.D0)+0.004D0*x3*x3)))+
     +     (x1*0.008D0*(6.D0/Pi)**(1.D0/3.D0)*
     *     x2**3.D0*x3*x3)/
     /     (x2**(8.D0/3.D0)+0.004D0*x3*x3)**2
C
C     F_xc with respect to zs
Cfah  see above (function dF_xs_by_drhos) for definition of the
Cfah  other variables
C     usage  dF_xs_by_dzs (n, c_xs, rhos, zs)
c
      dF_xs_by_dzs(n1, x1, x2, x3) = 
     =      (-3.d0*x1*n1*(3.D0/Pi)**(1.D0/3.D0)*
     *     x2**(4.D0/3.D0)*((0.004D0*x3*x3)/
     /     (x2**(8.D0/3.D0)+ 0.004D0*x3*x3))**(-1+n1)*
     *     ((-0.000032D0*x3*x3*x3)/(x2**(8.D0/3.D0)+
     +     0.004D0*x3*x3)**2+(0.008D0*x3)/
     /     (x2**(8.D0/3.D0)+0.004D0*x3*x3)))/
     /     (2.D0**(5.D0/3.D0))
      dF_xs_by_dzs1(x1, x2, x3) = 
     =      (-3.d0*x1*(3.D0/Pi)**(1.D0/3.D0)*
     *     x2**(4.D0/3.D0)*
     *     ((-0.000032D0*x3*x3*x3)/(x2**(8.D0/3.D0)+
     +     0.004D0*x3*x3)**2+(0.008D0*x3)/
     /     (x2**(8.D0/3.D0)+0.004D0*x3*x3)))/
     /     (2.D0**(5.D0/3.D0))
c      dF_xs_by_dzs0(x1, x2, x3) = 0d0
c
c     usage
c     e_crsa1 = e_crs1(rsa12,rsa,rsa32,rsa21)
c
      e_crs1(x1,x2,x3,x4) = -0.03108999999999999d0*
     *  dlog(1.d0 + 32.16468317787069d0/
     /     (14.1189d0*x1+6.1977d0*x2 + 3.3662d0*x3 +
     +     0.6251699999999999d0*x4))*(1.d0 + 0.20548d0*x2)
      drsbydrh(x1) = -((1.d0/x1)**(4.D0/3.D0)/
     -    (6.d0**(2.D0/3.D0)*Pi**(1.D0/3.D0)))
c     usage decrsdrs(rsa,rsa12,rsa21,rsa32)
      decrsdrs(x1,x2,x3,x4) = ((1.d0 + 0.20548d0*x1)*
     -     (6.1977d0 + 7.05945d0/x2 + 1.25034d0*x1+5.0493d0*x2))/
     -     ((6.1977d0*x1+14.1189d0*x2+0.6251699999999999d0*x3 + 
     +     3.3662d0*x4)**2d0*(1.d0 + 32.16468317787069d0/
     -     (6.1977d0*x1+14.1189d0*x2+0.6251699999999999d0*x3+
     +     3.3662d0*x4))) - 0.006388373199999999d0*
     -   dlog(1.d0 + 32.16468317787069d0/(6.1977d0*x1 + 14.1189d0*x2 + 
     -        0.6251699999999999d0*x3 + 3.3662d0*x4)) 
c
      DO j = 1, 4
        DO n = 1, (max_pow_u+1)*3
          F(n,j) = 0.D0
Cfah  later on, n has a different meaning: n as power of u, not 
Cfah  as number of the coefficient. 
        ENDDO
        dF_xa(j) = 0.D0
        dF_xb(j) = 0.D0
        dF_caa(j) = 0.D0
        dF_cbb(j) = 0.D0
        dF_cab(j) = 0.D0
C        DO k = 1, 4
C          DO n = 1, (max_pow_u+1)*3
C            FF(n,j,k) = 0.D0
C          ENDDO
C          d2F_xa(j,k) = 0.D0
C          d2F_xb(j,k) = 0.D0
C          d2F_caa(j,k) = 0.D0
C          d2F_cbb(j,k) = 0.D0
C          d2F_cab(j,k) = 0.D0
C        ENDDO
      ENDDO
      DO j = 1, (max_pow_u+1)*3
        F_xc(j) = 0.D0
      ENDDO 

Cfah --------------------------------------------------------------

Cfah call the expensive correlation parts here just once, and store their
Cfah values in a temporary variable. Then compute the actual F_c derivatives
Cfah with the various powers of u.  

      rho = rhoa + rhob

      s_a2=0.d0
      if(za.gt.tol_rho.and.rhoa.gt.tol_rho)
     +   s_a2 = za**2.D0 / rhoa**(8.D0/3.D0)
      s_b2=0.d0
      if(zb.gt.tol_rho.and.rhob.gt.tol_rho)
     +   s_b2 = zb**2.D0 / rhob**(8.D0/3.D0)
      s_avg2 = 0.5D0*(s_a2 + s_b2)

      u_caa = 0.2D0*s_a2/(1.D0+0.2D0*s_a2) 
      u_cbb = 0.2D0*s_b2/(1.D0+0.2D0*s_b2) 

      u_cab = 0.006D0*s_avg2/(1.d0+0.006D0*s_avg2)
      if(rhoa.gt.tol_rho) then
         rsa = ((3.d0/Pi)**(1.D0/3.D0)*
     -        (1.d0/rhoa)**(1.D0/3.D0))/
     -        2**(2.D0/3.D0)
         rsa12 = rsa**(1.D0/2.D0)
         rsa32 = rsa**(3.D0/2.D0)
         rsa21 = rsa**2.D0
      else
         rsa=0d0
         rsa12=0d0
         rsa32=0d0
         rsa21=0d0
      endif

      if(rhob.gt.tol_rho) then
         rsb = ((3.d0/Pi)**(1.D0/3.D0)*
     -        (1.d0/rhob)**(1.D0/3.D0))/
     -        2**(2.D0/3.D0)
         rsb12 = rsb**(1.D0/2.D0)
         rsb32 = rsb**(3.D0/2.D0)
         rsb21 = rsb**2.D0
C     
C     pw91 LDA Ecorr
C
         if(rhob.gt.tol_rho) then
            call xc_pw91ldag(rsb, 2, e_crsb1 ,de_crsb1_by_drsb, d2ez)
         else
            e_crsb1 = 0d0
            de_crsb1_by_drsb = 0d0
         endif
            
         du_cbb_by_drhob = (-1.6D0*zb**2*rhob**(5.D0/3.D0))/
     -        (3.d0*(0.2D0*zb**2 + 
     -        rhob**(8.D0/3.D0))**2)
         du_cbb_by_dzb = (2*0.2D0*zb*rhob**(8.D0/3.D0))/
     -        (0.2D0*zb**2 + rhob**(8.D0/3.D0))**2
         if(rhoa.gt.0d0) then
         du_cab_by_drhoa = (-16*0.006D0*za*za*rhoa**(5.D0/3.D0)*
     -        rhob**(16.D0/3.D0))/
     -        (3.d0*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -        0.006D0*za*za*rhob**(8.D0/3.D0) + 
     -        2.d0*rhoa**(8.D0/3.D0)*
     -        rhob**(8.D0/3.D0))**2) 
         du_cab_by_dza = (4*0.006D0*za*rhoa**(8.D0/3.D0)*
     -        rhob**(16.D0/3.D0))/
     -        (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -        0.006D0*za*za*rhob**(8.D0/3.D0) + 
     -        2.d0*rhoa**(8.D0/3.D0)*
     -        rhob**(8.D0/3.D0))**2
         du_cab_by_drhob = (-16*0.006D0*zb**2*rhob**(5.D0/3.D0)*
     -        rhoa**(16.D0/3.D0))/
     -        (3.d0*(0.006D0*za*za*rhob**(8.D0/3.D0) + 
     -        0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -        2.d0*rhob**(8.D0/3.D0)*
     -        rhoa**(8.D0/3.D0))**2) 
         du_cab_by_dzb = (4*0.006D0*zb*rhoa**(16.D0/3.D0)*
     -        rhob**(8.D0/3.D0))/
     -        (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -        0.006D0*za*za*rhob**(8.D0/3.D0) + 
     -        2.d0*rhoa**(8.D0/3.D0)*
     -        rhob**(8.D0/3.D0))**2
         else
            du_cab_by_drhoa = 0d0
            du_cab_by_dza = 0d0
            du_cab_by_drhob = 0d0
            du_cab_by_dzb = 0d0
         endif
         drsb_by_drhob = drsbydrh(rhob)
      else
         e_crsb1 = 0d0
         de_crsb1_by_drsb = 0d0
         du_cbb_by_drhob = 0d0
         du_cbb_by_dzb = 0d0
         du_cab_by_drhoa = 0d0
         du_cab_by_drhob = 0d0
         du_cab_by_dza = 0d0
         du_cab_by_dzb = 0d0
         drsb_by_drhob = 0d0
      endif

      rsab = ((3.d0/Pi)**(1.D0/3.D0)*
     -    (1.d0/rho)**(1.D0/3.D0))/
     -  2**(2.D0/3.D0)
      rsab12 = rsab**(1.D0/2.D0)
      rsab32 = rsab**(3.D0/2.D0)
      rsab21 = rsab**2.D0

      zeta = (rhoa-rhob)/rho
      if(zeta.lt.-1d0) zeta=-1d0
      if(zeta.gt.1d0) zeta=1d0

      if(abs(1d0-zeta).gt.tol_rho) then
         fzeta = (-2.d0 + sign(1d0,1.d0 - zeta)*(abs(1.d0 - zeta))**
     -        (4.D0/3.D0) +
     -        (1.d0 + zeta)**(4.D0/3.D0))/
     -        (-2.d0 + 2.d0*2.d0**(1.D0/3.D0))
         else
            fzeta = 1d0
         endif

C     
C     pw91 LDA Ecorr
C
      if(rhoa.gt.tol_rho) then
         call  xc_pw91ldag(rsa, 2, e_crsa1,
     D        de_crsa1_by_drsa, d2ez)
      else
         e_crsa1 = 0d0
      endif
      if(rho.gt.tol_rho) then
         call xc_pw91ldag(rsab, 2, e_crsab1,
     D        de_crsab1_by_drsab, d2ez)
c eps_c(rs,0)
         call xc_pw91ldag(rsab, 1, e_crsab0,
     D        de_crsab0_by_drsab, d2ez)
c     alpha(rs)
         call xc_pw91ldag(rsab, 3, a_crsab,
     D        da_crsab_by_drsab, d2ez)
      else
         e_crsab1 = 0d0
         e_crsab0 = 0d0
         a_crsab=0d0
      endif

      e_crsabzeta = e_crsab0+a_crsab*fzeta*(1.d0-zeta**4)/1.709921D0+
     -  (e_crsab1-e_crsab0)*fzeta*zeta**4

      e_caa = rhoa*e_crsa1
      e_cbb = rhob*e_crsb1
      e_cab = rho*e_crsabzeta - rhoa*e_crsa1 - rhob*e_crsb1
      if(rhoa.gt.tol_rho) then
         du_caa_by_drhoa = (-1.6D0*za*za*rhoa**(5.D0/3.D0))/
     -        (3.*(0.2D0*za*za + 
     -        rhoa**(8.D0/3.D0))**2)
         du_caa_by_dza = (2*0.2D0*za*rhoa**(8.D0/3.D0))/
     -        (0.2D0*za*za + rhoa**(8.D0/3.D0))**2
      else
         du_caa_by_drhoa = 0d0
         du_caa_by_dza = 0d0
      endif



Cfah Second derivatives are not required by cadpac. 
Cfah   du_caa_by_drhoa_dza = (16*0.2D0*za*rhoa**(5.D0/3.D0)*
Cfah -    (0.2D0*za**2 - rhoa**(8.D0/3.D0)))/
Cfah -  (3.*(0.2D0*za**2 + 
Cfah -       rhoa**(8.D0/3.D0))**3)
Cfah   du_cbb_by_drhob_dzb = (16*0.2D0*zb*rhob**(5.D0/3.D0)*
Cfah -    (0.2D0*zb**2 - rhob**(8.D0/3.D0)))/
Cfah -  (3.*(0.2D0*zb**2 + 
Cfah -       rhob**(8.D0/3.D0))**3)
Cfah
Cfah   du_caa_by_dza_dza = (2*0.2D0*rhoa**(8.D0/3.D0)*
Cfah -    (-3*0.2D0*za**2 + rhoa**(8.D0/3.D0))
Cfah -    )/
Cfah -  (0.2D0*za**2 + rhoa**(8.D0/3.D0))**3
Cfah   du_cbb_by_dzb_dzb = (2*0.2D0*rhob**(8.D0/3.D0)*
Cfah -    (-3*0.2D0*zb**2 + rhob**(8.D0/3.D0))
Cfah -    )/
Cfah -  (0.2D0*zb**2 + rhob**(8.D0/3.D0))**3
Cfah
Cfah   du_cab_by_drhoa_dza = (-32*0.006D0*rhoa**(5.D0/3.D0)*
Cfah -    (0.006D0*za*zb**2*
Cfah -       rhoa**(8.D0/3.D0)*
Cfah -       rhob**(16.D0/3.D0) - 
Cfah -      0.006D0*za**3*rhob**8 + 
Cfah -      2*za*rhoa**(8.D0/3.D0)*rhob**8))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_drhoa_dzb = (64*0.006D0**2*za**2*zb*
Cfah -    rhoa**(13.D0/3.D0)*
Cfah -    rhob**(16.D0/3.D0))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_drhob_dza = (64*0.006D0**2*za*zb**2*
Cfah -    rhoa**(16.D0/3.D0)*
Cfah -    rhob**(13.D0/3.D0))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_drhob_dzb = (-32*0.006D0*rhob**(5.D0/3.D0)*
Cfah -    (-(0.006D0*zb**3*rhoa**8) + 
Cfah -      0.006D0*za**2*zb*
Cfah -       rhoa**(16.D0/3.D0)*
Cfah -       rhob**(8.D0/3.D0) + 
Cfah -      2*zb*rhoa**8*rhob**(8.D0/3.D0)))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_dza_dza = (4*0.006D0*(0.006D0*zb**2*
Cfah -       rhoa**(16.D0/3.D0)*
Cfah -       rhob**(16.D0/3.D0) - 
Cfah -      3*0.006D0*za**2*rhoa**(8.D0/3.D0)*
Cfah -       rhob**8 + 2*rhoa**(16.D0/3.D0)*rhob**8
Cfah -      ))/
Cfah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -     2*rhoa**(8.D0/3.D0)*
Cfah -      rhob**(8.D0/3.D0))**3 
Cfah   du_cab_by_dza_dzb = (-16*0.006D0**2*za*zb*
Cfah -    rhoa**(16.D0/3.D0)*
Cfah -    rhob**(16.D0/3.D0))/
Cfah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -     2*rhoa**(8.D0/3.D0)*
Cfah -      rhob**(8.D0/3.D0))**3 
Cfah   du_cab_by_dzb_dzb = (4*0.006D0*rhoa**(16.D0/3.D0)*
Cfah -    rhob**(8.D0/3.D0)*
Cfah -    (-3*0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -      0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -      2*rhoa**(8.D0/3.D0)*
Cfah -       rhob**(8.D0/3.D0)))/
Cfah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -     2*rhoa**(8.D0/3.D0)*
Cfah -      rhob**(8.D0/3.D0))**3 

      if(rhoa.gt.tol_rho) then
         drsa_by_drhoa = drsbydrh(rhoa)
      else
         drsa_by_drhoa =0d0
      endif
      if(rho.gt.tol_rho) then
         drsab_by_drhoa = drsbydrh(rho) 
         drsab_by_drhob = drsab_by_drhoa 
      else
         drsab_by_drhoa = 0d0
         drsab_by_drhob = 0d0
      endif

      dzeta_by_drhoa = 2.d0*rhob/rho**2
      dzeta_by_drhob = -2.d0*rhoa/rho**2

      dfzeta_by_dzeta = ((-4.d0*sign(1d0,1.d0 - zeta)*
     *     (abs(1.d0 - zeta))**(1.D0/3.D0))/
     -     3.d0 + (4.d0*(1.d0 + zeta)**
     -        (1.D0/3.D0))/3.)/
     -  (-2.d0 + 2d0*2**(1.D0/3.D0))

      if(rhoa.gt.tol_rho) then
         call  xc_pw91ldag(rsa, 2, e_crsa1,
     D        de_crsa1_by_drsa, d2ez)
      else
         de_crsa1_by_drsa = 0d0
         de_crsab1_by_drsab = 0d0
      endif


      de_crsabzeta_by_drsab = 1.124999956683108D0*(1.d0 - zeta**4)*
     -   (-2.d0+sign(1d0,1.d0-zeta)*(abs(1.d0-zeta))**(4.D0/3.D0)+ 
     -     (1.d0 + zeta)**(4.D0/3.D0))*
     -   da_crsab_by_drsab + 
     -  de_crsab0_by_drsab + 
     -  fzeta*zeta**4*
     -   (- de_crsab0_by_drsab + 
     -      de_crsab1_by_drsab  )

      de_crsabzeta_by_dzeta = 1.499999942244144D0*(-1.d0 + zeta**4)*
     -   (sign(1d0,1.d0 - zeta)*(abs(1.d0 - zeta))**(1.D0/3.D0) - 
     -     (1.d0 + zeta)**(1.D0/3.D0))*a_crsab - 
     -  4.499999826732434D0*zeta**3*
     -   (-2.d0 + sign(1d0,1.d0-zeta)*abs(1.d0-zeta)**(4.D0/3.D0)+ 
     -     (1.d0 + zeta)**(4.D0/3.D0))*a_crsab + 
     -  (2*zeta**4*sign(1d0,1.d0-zeta)*(abs(1.d0-zeta)**(1.D0/3.D0)- 
     -       (1.d0 + zeta)**(1.D0/3.D0))*
     -     (e_crsab0 - e_crsab1))/
     -   (3.*(-1.d0 + 2**(1.D0/3.D0))) + 
     -  4*fzeta*(-e_crsab0 + 
     -     e_crsab1)*zeta**3

Cfah this is with application of the chain rule; I keep it that general
Cfah because this way, I only have to define one "G". 
      de_caa_by_drhoa = e_crsa1 + rhoa*de_crsa1_by_drsa*drsa_by_drhoa 
      de_cbb_by_drhob = e_crsb1 + rhob*de_crsb1_by_drsb*drsb_by_drhob 

      de_cab_by_drhoa = -e_crsa1 + 
     -  e_crsabzeta - 
     -  rhoa*de_crsa1_by_drsa* 
     -  drsa_by_drhoa + 
     -  rho*(de_crsabzeta_by_drsab*
     -  drsab_by_drhoa + 
     -  de_crsabzeta_by_dzeta*
     -  dzeta_by_drhoa)

      de_cab_by_drhob = -e_crsb1 + 
     -  e_crsabzeta - 
     -  rhob*de_crsb1_by_drsb* 
     -  drsb_by_drhob + 
     -  rho*(de_crsabzeta_by_drsab*
     -  drsab_by_drhob + 
     -  de_crsabzeta_by_dzeta*
     -  dzeta_by_drhob)



Cfah Here starts the big outer loop over the powers u 
      DO n = 0, max_pow_u 
        c_naa = sol((n*3) + 2)
        c_nbb = c_naa
        c_nab = sol((n*3) + 3) 

Cfah construction of the F_xc itself
Cfah -------------------------------
        IF (rhoa.GT.tol_rho) THEN
          if(n.eq.0) then
            F_xc(1) = F_xs0 (sol(1), rhoa, za)
          else
            F_xc(n*3+1) = F_xs (n, sol((n*3) + 1), rhoa, za)
          endif
        ENDIF
          if(u_caa.gt.tol_rho)then
             if(n.eq.0) then
             F_xc(2) = e_caa*c_naa
             else
             F_xc(n*3+2) = e_caa*u_caa**n*c_naa
             endif
          endif

        IF (rhob.GT.tol_rho) THEN
           if(n.eq.0) then
          F_xc(1) = F_xc(1)+F_xs0(sol(1), rhob, zb)
           else
          F_xc(n*3+1) = F_xc(n*3+1)+F_xs(n, sol((n*3) + 1), rhob, zb)
          endif
        ENDIF
          if(u_cbb.gt.tol_rho) then
             if(n.eq.0) then
                F_xc(2) = F_xc(2)+e_cbb*c_nbb
             else
                F_xc(n*3+2) = F_xc(n*3+2)+e_cbb*u_cbb**n*c_nbb
             endif
          endif

          if(u_cab.gt.tol_rho) then
             if(n.eq.0) then
                F_xc(3) = e_cab*c_nab
             else
                F_xc(n*3+3) = e_cab*u_cab**n*c_nab
             endif
          endif

Cfah       print*, 'in deriv:', e_cab, u_cab, c_nab

Cfah    First Derivatives
Cfah ---------------------

        if(za.gt.tol_rho)then 
           if(n.eq.0) then
            dF_xa(1) = dF_xs_by_drhos0 ( sol(1), rhoa, za) 
            dF_xa(3) = 0d0
         elseif(n.eq.1) then
            dF_xa(1) = dF_xs_by_drhos1 (sol(4), rhoa, za) 
            dF_xa(3) = dF_xs_by_dzs1 ( sol(4), rhoa, za)
         else
            dF_xa(1) = dF_xs_by_drhos (n, sol((n*3) + 1), rhoa, za) 
            dF_xa(3) = dF_xs_by_dzs (n, sol((n*3) + 1), rhoa, za)
         endif
          endif

        if(zb.gt.tol_rho) then
           if(n.eq.0) then
             dF_xb(2) = dF_xs_by_drhos0(sol(1), rhob, zb)
             dF_xb(4) = 0d0
          elseif(n.eq.1) then
             dF_xb(2) = dF_xs_by_drhos1(sol(4), rhob, zb)
             dF_xb(4) = dF_xs_by_dzs1 ( sol(4), rhob, zb)
          else
             dF_xb(2) = dF_xs_by_drhos (n, sol((n*3) + 1), rhob, zb)
             dF_xb(4) = dF_xs_by_dzs (n, sol((n*3) + 1), rhob, zb)
          endif
           endif

        if(u_caa.gt.tol_rho) then
           if(n.eq.0) then
          dF_caa(1) = c_naa*de_caa_by_drhoa
          dF_caa(3) = 0d0
           elseif(n.eq.1) then
          dF_caa(1) = c_naa*u_caa*
     *       de_caa_by_drhoa+c_naa*e_caa*du_caa_by_drhoa
          dF_caa(3) = c_naa*e_caa*du_caa_by_dza
          else
          dF_caa(1) = c_naa*u_caa**n*
     *       de_caa_by_drhoa+c_naa*n*e_caa*u_caa**(-1+n)*
     *       du_caa_by_drhoa
          dF_caa(3) = c_naa*n*e_caa*u_caa**(-1+n)*du_caa_by_dza
          endif
        endif

        if(u_cbb.gt.tol_rho) then
           if(n.eq.0) then
          dF_cbb(2) = c_nbb*de_cbb_by_drhob 
          dF_cbb(4) = 0d0
           elseif(n.eq.1) then
          dF_cbb(2) = c_nbb*u_cbb*de_cbb_by_drhob + 
     -   c_nbb*e_cbb*du_cbb_by_drhob
          dF_cbb(4) = c_nbb*e_cbb*du_cbb_by_dzb
           else
          dF_cbb(2) = c_nbb*u_cbb**n*de_cbb_by_drhob + 
     -   c_nbb*n*e_cbb*u_cbb**(-1 + n)*du_cbb_by_drhob
          dF_cbb(4) = c_nbb*n*e_cbb*u_cbb**(-1+n)*du_cbb_by_dzb
          endif
        endif


        if(u_cab.gt.tol_rho) then
           if(n.eq.0) then
          dF_cab(1) = c_nab*de_cab_by_drhoa
          dF_cab(2) = c_nab*de_cab_by_drhob
          dF_cab(3) = 0d0
          dF_cab(4) = 0d0
           elseif(n.eq.1) then
          dF_cab(1) = c_nab*u_cab*
     *         de_cab_by_drhoa+c_nab*n*e_cab*du_cab_by_drhoa 
          dF_cab(2) = c_nab*u_cab*
     -         de_cab_by_drhob+c_nab*n*e_cab*du_cab_by_drhob
          dF_cab(3) = c_nab*e_cab*du_cab_by_dza
          dF_cab(4) = c_nab*n*e_cab*du_cab_by_dzb
          else
          dF_cab(1) = c_nab*u_cab**n*
     *         de_cab_by_drhoa+c_nab*n*e_cab*u_cab**(-1+n)*
     *         du_cab_by_drhoa 
          dF_cab(2) = c_nab*u_cab**n*
     -         de_cab_by_drhob+c_nab*n*e_cab*u_cab**(-1+n)*
     -         du_cab_by_drhob
          dF_cab(3) = c_nab*n*e_cab*
     -         u_cab**(-1+n)*du_cab_by_dza
          dF_cab(4) = c_nab*n*e_cab*
     -         u_cab**(-1+n)*du_cab_by_dzb
          endif
        endif

Cfah Second Derivatives
Cfah ------------------

Cfah         d2F_xa(1,1) = d2F_xs_by_drhos_drhos (n, sol((n*3) + 1), 
Cfah      &                                         rhoa, za)
Cfah see comment below, for the (2,2) term. 
Cfah    d2F_xa(1,2) = 0 
Cfah    d2F_xa(1,3) = d2F_xs_by_drhos_dzs (n, sol((n*3) + 1), rhoa, za)
Cfah    d2F_xa(1,4) = 0 
Cfah    d2F_xa(2,2) = 0  
Cfah    d2F_xa(2,3) = 0  
Cfah    d2F_xa(2,4) = 0 
Cfah    d2F_xa(3,3) = d2F_xs_by_dzs_dzs (n, sol((n*3) + 1), rhoa, za)
Cfah    d2F_xa(3,4) = 0 
Cfah    d2F_xa(4,4) = 0 

Cfah for alpha spin, elements are non-zero when both indices are odd; 
Cfah for beta spin, elements are non-zero when both indices are even. 
Cfah the matrix is symmetric, and the upper triangle contains the 
Cfah 10 elements given above and below. 

Cfah    d2F_xb(1,1) = 0
Cfah    d2F_xb(1,2) = 0  
Cfah    d2F_xb(1,3) = 0           
Cfah    d2F_xb(1,4) = 0 
Cfah        d2F_xb(2,2) = d2F_xs_by_drhos_drhos (n, sol((n*3) + 1), 
Cfah     &                                         rhob, zb)
Cfah this term is NOT zero, but needs not be evaluated since we don't 
Cfah need it for the construction of v (cf. routine "va" in the fit 
Cfah program) 
Cfah    d2F_xb(2,3) = 0  
Cfah    d2F_xb(2,4) = d2F_xs_by_drhos_dzs (n, sol((n*3) + 1), rhob, zb)
Cfah    d2F_xb(3,3) = 0
Cfah    d2F_xb(3,4) = 0 
Cfah    d2F_xb(4,4) = d2F_xs_by_dzs_dzs (n, sol((n*3) + 1), rhob, zb)


Cfah    d2F_caa(1,1) = !=0, but not needed 
Cfah    d2F_caa(1,2) = 0.D0 (not needed)  
Cfah    d2F_caa(1,3) = c_naa*n*u_caa**(-1 + n)*
Cfah -   de_caa_by_drhoa*
Cfah -   du_caa_by_dza + 
Cfah -  c_naa*(-1 + n)*n*e_caa*
Cfah -   u_caa**(-2 + n)*
Cfah -   du_caa_by_dza*
Cfah -   du_caa_by_drhoa + 
Cfah -  c_naa*n*e_caa*u_caa**(-1 + n)*
Cfah -   du_caa_by_drhoa_dza
Cfah    d2F_caa(1,4) = 0.D0
Cfah    d2F_caa(2,2) = 0.D0 (not needed)
Cfah    d2F_caa(2,3) = 0.D0 
Cfah    d2F_caa(2,4) = 0.D0 
Cfah    d2F_caa(3,3) = c_naa*n*e_caa*u_caa**(-2 + n)*
Cfah -  ((-1 + n)*du_caa_by_dza**
Cfah -      2 + u_caa*
Cfah -     du_caa_by_dza_dza)
Cfah    d2F_caa(3,4) = 0.D0 
Cfah    d2F_caa(4,4) = 0.D0 


Cfah    d2F_cbb(1,1) = 0.D0 (not needed)
Cfah    d2F_cbb(1,2) = 0.D0 (not needed)
Cfah    d2F_cbb(1,3) = 0.D0
Cfah    d2F_cbb(1,4) = 0.D0
Cfah    d2F_cbb(2,2) = !=0, but not needed 
Cfah    d2F_cbb(2,3) = 0.D0
Cfah    d2F_cbb(2,4) = c_nbb*n*u_cbb**(-1 + n)*
Cfah -   de_cbb_by_drhob*
Cfah -   du_cbb_by_dzb +
Cfah -  c_nbb*(-1 + n)*n*e_cbb*
Cfah -   u_cbb**(-2 + n)*
Cfah -   du_cbb_by_dzb*
Cfah -   du_cbb_by_drhob +
Cfah -  c_nbb*n*e_cbb*u_cbb**(-1 + n)*
Cfah -   du_cbb_by_drhob_dzb
Cfah    d2F_cbb(3,3) = 0.D0
Cfah    d2F_cbb(3,4) = 0.D0
Cfah    d2F_cbb(4,4) =  c_nbb*n*e_cbb*u_cbb**(-2 + n)*
Cfah -  ((-1 + n)*du_cbb_by_dzb**
Cfah -      2 + u_cbb*
Cfah -     du_cbb_by_dzb_dzb)

Cfah    d2F_cab(1,1) = not needed
Cfah    d2F_cab(1,2) = not needed
Cfah    d2F_cab(1,3) = c_nab*n*u_cab**(-2 + n)*
Cfah -  ((-1 + n)*e_cab*
Cfah -     du_cab_by_dza
Cfah -      *du_cab_by_drhoa + 
Cfah -    u_cab*(de_cab_by_drhoa*
Cfah -        du_cab_by_dza + e_cab*
Cfah -        du_cab_by_drhoa_dza
Cfah -  ))
Cfah    d2F_cab(1,4) = c_nab*n*u_cab**(-2 + n)*
Cfah -  (  (-1 + n)*e_cab*
Cfah -     du_cab_by_dzb
Cfah -      *du_cab_by_drhoa + 
Cfah -    u_cab*
Cfah -     (de_cab_by_drhoa*
Cfah -        du_cab_by_dzb + 
Cfah -       e_cab*
Cfah -        du_cab_by_drhoa_dzb))
Cfah    d2F_cab(2,2) = not needed
Cfah    d2F_cab(2,3) = c_nab*n*u_cab**(-2 + n)*
Cfah -  ((-1 + n)*e_cab*
Cfah -     du_cab_by_dza
Cfah -      *du_cab_by_drhob +
Cfah -    u_cab*
Cfah -     (de_cab_by_drhob*
Cfah -        du_cab_by_dza + 
Cfah -       e_cab*
Cfah -        du_cab_by_drhob_dza))
Cfah    d2F_cab(2,4) = c_nab*n*u_cab**(-2 + n)*
Cfah -  ((-1 + n)*e_cab*
Cfah -     du_cab_by_dzb
Cfah -      *du_cab_by_drhob + 
Cfah -    u_cab*(de_cab_by_drhob*
Cfah -        du_cab_by_dzb + e_cab*
Cfah -        du_cab_by_drhob_dzb ))
Cfah    d2F_cab(3,3) = c_nab*n*e_cab*
Cfah -  u_cab**(-2 + n)*
Cfah -  ((-1 + n)*du_cab_by_dza**2 + 
Cfah -    u_cab*
Cfah -     du_cab_by_dza_dza)
Cfah    d2F_cab(3,4) = c_nab*n*e_cab*
Cfah -  u_cab**(-2 + n)*
Cfah -  ((-1 + n)*du_cab_by_dzb*
Cfah -     du_cab_by_dza
Cfah -       + u_cab*
Cfah -     du_cab_by_dza_dzb)
Cfah    d2F_cab(4,4) = c_nab*n*e_cab*
Cfah -  u_cab**(-2 + n)*
Cfah -  ((-1 + n)*du_cab_by_dzb**2 + 
Cfah -    u_cab*
Cfah -     du_cab_by_dzb_dzb)
Cfah
Cfah here, the second derivatives are completed (Schwartz's rule: 
Cfah df/(dadb) = df/(dbda) 
Cfah    DO i = 1, 4
Cfah      DO j = i, 4
Cfah        d2F_xa(j,i) = d2F_xa(i,j)  
Cfah        d2F_xb(j,i) = d2F_xb(i,j)
Cfah        d2F_caa(j,i) = d2F_caa(i,j) 
Cfah        d2F_cbb(j,i) = d2F_cbb(i,j) 
Cfah        d2F_cab(j,i) = d2F_cab(i,j) 
Cfah      ENDDO
Cfah    ENDDO

Cfah test for zero densities (as in beta part of H atom):
        IF (rhob.LT.tol_rho) THEN 
          DO i = 1, 4
            dF_xb(i) = 0.D0
            dF_cbb(i) = 0.D0
            dF_cab(i) = 0.D0
Cfah        DO j = 1, 4
Cfah          d2F_xb(i,j) = 0.D0
Cfah          d2F_cbb(i,j) = 0.D0
Cfah          d2F_cab(i,j) = 0.D0
Cfah        ENDDO 
          ENDDO 
        ENDIF

        IF (rhoa.LT.tol_rho) THEN
          DO i = 1, 4
            dF_xa(i) = 0.D0
            dF_caa(i) = 0.D0
            dF_cab(i) = 0.D0
Cfah        DO j = 1, 4
Cfah          d2F_xa(i,j) = 0.D0
Cfah          d2F_caa(i,j) = 0.D0
Cfah          d2F_cab(i,j) = 0.D0
Cfah        ENDDO
          ENDDO
        ENDIF


Cfah Sum up all the partial derivatives with respect to the same function
Cfah of terms containing different powers of u with the help of the big outer 
Cfah loop: 

Cfah have the partial derivative 

        DO i = 1, 4
          F(n*3+1,i) = dF_xa(i) + dF_xb(i) 
          F(n*3+2,i) = dF_caa(i) +dF_cbb(i) 
          F(n*3+3,i) = dF_cab(i) 
Cfah      DO j = 1, 4
Cfah        FF(n*3+1,i,j) = d2F_xa(i,j) + d2F_xb(i,j) 
Cfah        FF(n*3+2,i,j) = d2F_caa(i,j) + d2F_cbb(i,j) 
Cfah        FF(n*3+3,i,j) = d2F_cab(i,j)  
Cfah      ENDDO 
        ENDDO 

Cfah these partial derivatives have not been computed because they are
Cfah zero since we don't have a gradrhoagradrhob term in the Becke V functional
C        F(n*3+1,5) = 0
C        F(n*3+2,5) = 0
C        F(n*3+3,5) = 0
Cfah    DO i = 1, 5
Cfah      FF(n*3+1,i,5) = 0
Cfah      FF(n*3+2,i,5) = 0
Cfah      FF(n*3+3,i,5) = 0
Cfah      FF(n*3+1,5,i) = 0
Cfah      FF(n*3+2,5,i) = 0
Cfah      FF(n*3+3,5,i) = 0
Cfah    ENDDO

      ENDDO

      RETURN

      END

Cfah-----------------------------------------------------------

Cfah-----------------------------------------------------------

ccc   DFT_XCKernel_PWLDA(ra, rb, FCLDA);
      subroutine DFT_XCKernel_PWLDA(ra, rb, FCLDA)
      implicit none
#include "2ndDerivB97.h"
      double precision ra, rb
      double precision FCLDA(0:_FXC_RARB)
c
      double precision Amat(2)
      double precision Amat2(3)
      double precision rho(3)
      double precision ec,qwght
c
      double precision func
c
      rho(1)=ra+rb
      rho(2)=ra
      rho(3)=rb
cinitialize      
      ec = 0d0
      qwght=1d0
      Amat(1)=0d0
      Amat(2)=0d0
      Amat2(1)=0d0
      Amat2(2)=0d0
      Amat2(3)=0d0

      call xc_pw91lda_d2(1d-20, 1d0, .true., .false., rho,
     &     Amat, Amat2, 1, 2, Ec, qwght, .false., func)

      FCLDA(_FXC_E)=ec
      FCLDA(_FXC_RA)=Amat(1)
      FCLDA(_FXC_RB)=Amat(2)
c
      FCLDA(_FXC_RARA)=Amat2(D2_RA_RA)
      FCLDA(_FXC_RBRB)=Amat2(D2_RB_RB)
      FCLDA(_FXC_RARB)=Amat2(D2_RA_RB)
c      rubbish
      return
      end
      subroutine xc_fxc2acmat(nq,i,ipol,FX,
     A     Amat,Cmat,Amat2,Cmat2)
      implicit none
#include "2ndDerivB97.h"
#include "mafdecls.fh"
      integer nq,ipol
      integer i
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
      double precision FX(0:_FXC_NUMDERI)
      Amat(i,1) = Amat(i,1) + FX(_FXC_RA) 
      Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FX(_FXC_GAA)
      Amat2(i,D2_RA_RA) = Amat2(i,D2_RA_RA) + FX(_FXC_RARA)
      Amat2(i,D2_RA_RB) = Amat2(i,D2_RA_RB) + FX(_FXC_RARB)
      Cmat2(i,D2_GAA_GAA) = Cmat2(i,D2_GAA_GAA) +
     +     FX(_FXC_GAAGAA)
      Cmat2(i,D2_GAA_GBB) = Cmat2(i,D2_GAA_GBB) +
     +     FX(_FXC_GAAGBB)
      Cmat2(i,D2_RA_GAA) = Cmat2(i,D2_RA_GAA) +
     +     FX(_FXC_RAGAA)
      Cmat2(i,D2_RA_GBB) = Cmat2(i,D2_RA_GBB) +
     +     FX(_FXC_RAGBB)
      Cmat2(i,D2_RB_GAA) = Cmat2(i,D2_RB_GAA) +
     +     FX(_FXC_RBGAA)
      if(ipol.eq.2) then
         Amat(i,2) = Amat(i,2) + FX(_FXC_RB) 
         Cmat(i,D1_GBB) = Cmat(i,D1_GBB) + FX(_FXC_GBB)
         Amat2(i,D2_RB_RB) = Amat2(i,D2_RB_RB) + FX(_FXC_RBRB)
         Cmat2(i,D2_GBB_GBB) = Cmat2(i,D2_GBB_GBB) +
     +        FX(_FXC_GBBGBB)
         Cmat2(i,D2_RB_GBB) = Cmat2(i,D2_RB_GBB) +
     +        FX(_FXC_RBGBB)
      endif
      return
      end
      subroutine xc_htch_loadcf(functional,sol,max_pow_u)
      implicit none
#include "errquit.fh"
      character*4 functional ! functional name [input]
      integer max_pow_u ! functional name [input]
#include "xc_hcth.fh"
      double precision sol((limpow+1)*3) ! [out]
c
      integer nofunc
c
      nofunc = -1               ! take care of compiler warnings
      do n=1,numfunc
        if(functional.eq.funcnam(n)) nofunc=n
      enddo
      if(nofunc.eq.-1) call errquit('xchcth: cant pair funcname ',0,
     &       UNKNOWN_ERR)
      max_pow_u=maxpow(nofunc)
      do n=1,3*(limpow+1)
        sol(n)=coeffs(n,nofunc)
      enddo

C please refer to these coeffs as THCH1/iterate-e750-g500-v1-m4-n4

c      sol( 1) =     0.109320D+01
c      sol( 2) =     0.222601D+00
c      sol( 3) =     0.729974D+00
c      sol( 4) =    -0.744056D+00
c      sol( 5) =    -0.338622D-01
c      sol( 6) =     0.335287D+01
c      sol( 7) =     0.559920D+01
c      sol( 8) =    -0.125170D-01
c      sol( 9) =    -0.115430D+02
c      sol(10) =    -0.678549D+01
c      sol(11) =    -0.802496D+00
c      sol(12) =     0.808564D+01
c      sol(13) =     0.449357D+01
c      sol(14) =     0.155396D+01
c      sol(15) =    -0.447857D+01
      return
      end
