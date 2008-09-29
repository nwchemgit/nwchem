      Subroutine xc_hcth(tol_rho, xfac, lxfac, nlxfac,
     ,                    cfac, lcfac, nlcfac,rho, delrho, 
     &                    Amat, Cmat, nq, ipol, Ex, Ec,  qwght,
     &                    ldew,func,funcname)     
c
c$Id: xc_hcth.F,v 1.35 2008-09-29 20:25:06 edo Exp $
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
            if(gammaval.gt.tol_rho)  za=sqrt(gammaval)
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
            if(gammaval.gt.tol_rho) za=sqrt(gammaval)
            zb=0.d0
            gammaval=delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &           delrho(n,3,2)*delrho(n,3,2) 
            if(gammaval.gt.tol_rho) zb=sqrt(gammaval)
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
c$Id: xc_hcth.F,v 1.35 2008-09-29 20:25:06 edo Exp $
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

      integer limpow,n
      PARAMETER(limpow = 4)
Cfah limpow is equivalent to "m" in the Becke V paper, that is, the greatest 
Cfah power of u appearing in the power expansion. 
Cfah The "Becke V" paper is: 
Cfah Becke A. D.  Density-functional thermochemistry. V.
Cfah Systematic optimization of exchange-correlation functionals, 
Cfah J. Chem. Phys., 1997, 107, 8554-8560
c
c     variables passed to hcderiv
c
      integer numfunc,nofunc,max_pow_u
      parameter(numfunc=14)
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

Cfah these transformed variables u will be defined and given short-cut names 
Cfah below. 
      double precision coeffs(3*(limpow+1),numfunc)
      character*4 funcnam(numfunc)
      integer maxpow(numfunc)
      data maxpow  /     2,   2 ,      2,    4,      4,    4,     4,
     ,              2 ,    4 ,    4    ,   2  ,  4  ,  4 ,        2/

      data funcnam/'b970','b980','b971','hcth','hcta','h120','h147',
     ,             'b97g','h407','hp14','b972','407p','b973','b97d'/
C               B97        B98   B97-1          HCTH        HCTH-A
C m max          2,            2,            4,           4,
Cc X s,0
      data (coeffs(1,n),n=1,numfunc)/
     /     +0.80940d+00,0.790194d0, +0.789518d+00,+0.109320d+01,
     ,     +0.109878d+01,1.09163d0,  1.09025d0, 1.1068d0,   1.08184d0,
     ,     +0.103161d+01,+0.827642D+00,+1.08018D0,+7.334648D-01,
     ,     +1.086620d+00/

c C ss,0 
      data (coeffs(2,n),n=1,numfunc)/
     ,     +0.17370d+00,-0.120163d0,+0.820011d-01,+0.222601d+00,
     ,     +0.136823d-01, 0.48951d0,  0.56258d0, 0.4883d0,  1.18777d0,
     ,     +0.282414d+01,+0.585808D+00,+0.80302D0,+5.623649D-01,
     ,     +0.22340d+00/
c C ab,0 
      data (coeffs(3,n),n=1,numfunc)/
     ,     +0.94540d+00,0.934715d0,+0.955689d+00,+0.729974d+00,
     ,     +0.836897d+00,0.51473d0,  0.54235d0,  0.7961d0,  0.58908d0,
     ,     +0.821827d-01,+0.999849D+00,+0.73604D0,+1.133830D+00,
     ,     +0.690410d+00/
c X s,1 
      data (coeffs(4,n),n=1,numfunc)/
     ,     +0.50730d+00,0.400271d0,+0.573805d+00,-0.744056d+00,
     ,     -0.251173d+01,-0.74720d0, -0.79920d0, -0.8765d0, -0.5183d0,
     ,     -0.360781d+00,+0.478400D-01,-0.4117D0,+2.925270D-01,
     ,     -0.521270d+00/
c C ss,1 
      data (coeffs(5,n),n=1,numfunc)/
     ,     +0.23487d+01,2.82332d0,+0.271681d+01,-0.338622d-01,
     ,     +0.268920d+00,-0.26070d0, -0.01710d0, -2.117d0, -2.4029d0,
     ,     +0.318843d-01,-0.691682D+00,-1.0479D0,-1.322980D+00,
     ,     -1.562080d+00/
c C ab,1 
      data (coeffs(6,n),n=1,numfunc)/
     ,     +0.74710d+00,1.14105d0,+0.788552d+00,+0.335287d+01,
     ,     +0.172051d+01,6.92980d0,  7.01460d0, 5.7060d0,  4.4237d0,
     ,     +0.456466d+01,+0.140626D+01,+3.0270D0,-2.811967D+00,
     ,     +6.302700d+00/
c X s,2 
      data (coeffs(7,n),n=1,numfunc)/
     ,     +0.74810d+00,0.832857d0,+0.660975d+00,+0.559920d+01,
     ,     +0.156233d-01,5.07830d0,  5.57210d0, 4.2639d0,   3.4256d0,
     ,     +0.351994d+01,+0.176125D+01,+2.4368D0,+3.338789D+00,
     ,     +3.254290d+00/
c C ss,2
      data (coeffs(8,n),n=1,numfunc)/
     ,     -0.24868d+01,-2.59412d0,-0.287103d+01,-0.125170d-01,
     ,     -0.550769d+00,0.43290d0, -1.30640d0, 2.3235d0,   5.6174d0,
     ,     -0.178512d+01,+0.394796D+00,+4.9807D0,+6.359191D+00,
     ,     +1.942930d+00/
c C ab,2
      data (coeffs(9,n),n=1,numfunc)/
     ,     -0.45961d+01,-5.33398d0,-0.547869d+01,-0.115430d+02,
     ,     -0.278498d+01,-24.7070d0, -28.3820d0,-14.9820d0,-19.222d0,
     ,     -0.135529d+02,-0.744060D+01,-10.075D0,+7.431302D+00,
     ,     -14.97120d+00/
c X s,3 
      data (coeffs(10,n),n=1,numfunc)/
     ,      0.0000000d0,0.0d000000, 0.00000000d0,-0.678549d+01,
     ,     0.00000000d0,-4.10750d0, -5.86760d0 ,0d0       , -2.6290d0,
     ,     -0.495944d+01,0.d0000000000,+1.3890D0,-1.051158D+01,
     ,      0.0000000d0/
c C ss,3
      data (coeffs(11,n),n=1,numfunc)/
     ,      0.0000000d0,0.0d000000, 0.00000000d0,-0.802496d+00,
     ,     +0.103947d+01,-1.99250d0,  1.05750d0,0d0       , -9.1792d0,
     ,     +0.239795d+01,0.d0000000000,-12.890D0,-7.464002D+00,
     ,      0.0000000d0/
c C ab,3
      data (coeffs(12,n),n=1,numfunc)/
     ,      0.0000000d0,0.0d000000, 0.00000000d0,+0.808564d+01,
     ,     -0.457504d+01, 23.1100d0,  35.0330d0,0d0       , 42.572d0,
     ,     +0.133820d+02,0.000000d0,+20.611D0,-1.969342D+00,
     ,      0.0000000d0/
c X s,4 
      data (coeffs(13,n),n=1,numfunc)/
     ,      0.0000000d0,0.0d000000, 0.00000000d0,+0.449357d+01,
     ,     0.00000000d0, 1.17170d0,  3.04540d0, 0d0       , 2.2886d0,
     ,     +0.241165d+01,0.000000d0,-1.3529D0,+1.060907D+01,
     ,      0.0000000d0/
c C ss,4
      data (coeffs(14,n),n=1,numfunc)/
     ,      0.0000000d0,0.0d000000, 0.00000000d0,+0.155396d+01,
     ,     0.00000000d0, 2.48530d0,  0.88540d0, 0d0       , 6.2480d0,
     ,     -0.876909d+00,0.000000d0,9.6446D0,+1.827082D+00,
     ,      0.0000000d0/
c C ab,4
      data (coeffs(15,n),n=1,numfunc)/
     ,      0.0000000d0,0.0d000000, 0.00000000d0,-0.447857d+01,
     ,     0.00000000d0,-11.3230d0, -20.4280d0, 0d0, -42.005d0,
     ,     -0.317493d+01,0.000000d0,-29.418D0,-1.174423D+01,
     ,      0.0000000d0/
c
c X +0.1! coeffs for HF exchange
c
C     ,     +0.19430d+00,+0.210000d+00, 0.00000000d0, 0.00000000d0/

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
      double precision decrsdrs
      double precision eps
      parameter (eps=1d-19)
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
      if(za.gt.tol_rho) s_a2 = za**2.D0 / rhoa**(8.D0/3.D0)
      s_b2=0.d0
      if(zb.gt.tol_rho) s_b2 = zb**2.D0 / rhob**(8.D0/3.D0)
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
            e_crsb1 = e_crs1(rsb12,rsb,rsb32,rsb21)
            de_crsb1_by_drsb =  decrsdrs(rsb,rsb12,rsb21,rsb32)
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

      if(abs(1d0-zeta).gt.eps) then
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
         e_crsa1 = e_crs1(rsa12,rsa,rsa32,rsa21)
      else
         e_crsa1 = 0d0
      endif
      if(rho.gt.tol_rho) then
         e_crsab1 = e_crs1(rsab12,rsab,rsab32,rsab21)
         e_crsab0 = -0.062182d0*dlog(1.d0 +
     -        16.0818243221511d0/
     -        (7.595699999999999d0*rsab12 +
     -        3.5876d0*rsab +
     -        1.6382d0*rsab32 +
     -        0.49294d0*rsab21))*
     -        (1.d0 + 0.2137d0*rsab)

         a_crsab = 0.03377399999999999d0*
     -        dlog(1.d0 + 29.60857464321667/
     -        (10.35699999999999d0*rsab12 +
     -        3.623099999999999d0*rsab +
     -        0.88026d0*rsab32 +
     -        0.49671d0*rsab21))*
     -        (1.d0 + 0.11125d0*rsab)
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


      da_crsab_by_drsab = (-1.d0*(1.d0 + 0.11125d0*rsab)*
     -     (3.623099999999999d0 + 
     -       5.178499999999999d0/rsab12 + 
     -       0.99342d0*rsab + 1.32039d0*rsab12))/
     -   ((1.d0 + 29.60857464321667d0/
     -        (3.623099999999999d0*rsab + 
     -          10.35699999999999d0*rsab12 + 
     -          0.49671d0*rsab21 + 0.88026d0*rsab32))*
     -     (3.623099999999999d0*rsab + 
     -        10.35699999999999d0*rsab12 + 
     -        0.49671d0*rsab21 + 0.880260*rsab32)**2) + 
     -  0.003757357499999999d0*
     -   dlog(1.d0 + 29.60857464321667/
     -      (3.623099999999999d0*rsab + 
     -        10.35699999999999d0*rsab12 + 
     -        0.49671d0*rsab21 + 0.88026d0*rsab32)) 

      de_crsab0_by_drsab = (1.d0*(1.d0 + 0.2137d0*rsab)*
     -     (3.5876d0 + 3.797849999999999d0/rsab12 + 
     -       0.98588d0*rsab + 2.4573*rsab12))/
     -   ((3.5876d0*rsab + 7.595699999999999d0*rsab12 + 
     -        0.49294d0*rsab21 + 1.6382d0*rsab32)**2*
     -     (1.d0 + 16.0818243221511d0/
     -        (3.5876d0*rsab + 7.595699999999999d0*rsab12 + 
     -          0.49294d0*rsab21 + 1.6382d0*rsab32))) - 
     -  0.01328829339999999d0*
     -   dlog(1.d0 + 16.0818243221511d0/
     -      (3.5876*rsab + 7.595699999999999d0*rsab12 + 
     -        0.49294d0*rsab21 + 1.6382d0*rsab32))
      if(rhoa.gt.tol_rho) then
         de_crsa1_by_drsa = decrsdrs(rsa,rsa12,rsa21,rsa32)
         de_crsab1_by_drsab = decrsdrs(rsab,rsab12,rsab21,rsab32)
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
          if(abs(u_caa).gt.tol_rho)then
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
          if(abs(u_cbb).gt.tol_rho) then
             if(n.eq.0) then
                F_xc(2) = F_xc(2)+e_cbb*c_nbb
             else
                F_xc(n*3+2) = F_xc(n*3+2)+e_cbb*u_cbb**n*c_nbb
             endif
          endif

          if(abs(u_cab).gt.tol_rho) then
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

        if(abs(u_caa).gt.tol_rho) then
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

        if(abs(u_cbb).gt.tol_rho) then
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


        if(abs(u_cab).gt.tol_rho) then
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

