
      subroutine cs_minres_qlp(ov1,ov2,factor,x,b,npoles,maxpoles,n,
     $                         ipol)
      implicit none
#include "mafdecls.fh"
      integer,intent(in) :: maxpoles,npoles(*),n,ipol
      double precision,intent(in) :: ov1(n,*),ov2(n,*),b(n)
      double complex,intent(in) :: factor(maxpoles,*)
      double complex,intent(inout) :: x(n)

      double complex, parameter :: zzero = (0d0,0d0)
      integer, parameter :: maxiter = 30
      double precision, parameter :: rtol = 1d-4, maxxnorm=1d7,
     $                               Acondlim = 1d15, TranCond=1d7 
      logical,parameter :: likeLS = .false.

      double precision, dimension(n) :: resvec(maxiter),aresvec(maxiter)
      double complex, dimension(n) :: r1,r2,r3,v,w,wl,wl2,xl,xl2,xbak

      double precision :: beta1,beta,betal,betan,gmin,cs,cr1,cr2
      double precision :: Acond,Acondl,Anorm,relres,relresl
      double precision :: xnorm,t2,t1,relAresl,rnorm,xnorm_tmp
      double precision :: xl2norm,Axnorm,abs_gama,epsx,pnorm
      double precision :: gminl,gminl2,Anorml,rootl,realmin
      double precision :: eps,Arnorml,Arnorm,rnorml,xnorml
      double complex :: tau,taul,phi,sn,sr1,sr2,dltan
      double complex :: gama,gamal,gamal2,eta,etal,etal2,vepln
      double complex :: veplnl,veplnl2,ul3,ul2,ul,u,ul4
      double complex :: gamal_QLP,gamal_tmp,vepln_QLP,gama_QLP
      double complex :: ul_QLP,u_QLP,alpha,gbar,dlta,dlta_tmp
      double complex :: gamal3,taul2,dbar,eplnn,dlta_QLP,epln
      double complex :: gama_tmp
      double precision, external :: dznrm2
      integer,external :: ga_nodeid

      integer flag, flag0, iter, lines, headlines, QLPiter


      call ga_sync()

      realmin = tiny(1d0)
      eps = epsilon(0d0)

      if (ipol.eq.1) then
       call gw_aprod_cmplx(OV1,factor,x,r2,4d0,0d0,npoles(1),n)
      else
       call gw_aprod_cmplx(OV1,factor,x,r2,2d0,0d0,npoles(1),n)
       call gw_aprod_cmplx(OV2,factor(1,2),x,r2,2d0,1d0,
     $                      npoles(2),n)
      endif
      call ga_dgop(1038,r2,n,'+')
      r2 = dcmplx(b,0d0) - (r2 + x)
      xbak = x

c      r2(:) = b(:)
      r3(:) = r2(:)
      beta1 = dsqrt(dble(dot_product(r2,r2)))

      if (beta1.lt.rtol) then
        return
      endif

      flag0 = -2
      flag = flag0
      iter = 0
      QLPiter = 0
      lines = 1
      headlines = 20
      beta = 0d0
      tau  = zzero
      taul = zzero
      phi = dcmplx(beta1,0d0)
      betan = beta1
      gmin = 0d0
      cs = -1d0
      sn = zzero
      cr1 = 1d0
      sr1 = zzero
      cr2 = -1d0
      sr2 = zzero
      dltan = zzero
      eplnn = zzero
      gama = zzero
      gamal = zzero
      gamal2 = zzero
      eta = zzero
      etal = zzero
      etal2 = zzero
      vepln = zzero
      veplnl = zzero
      veplnl2 = zzero
      ul3 = zzero
      ul2 = zzero
      ul = zzero
      u = zzero
      rnorm = betan
      xnorm = 0d0
      xl2norm = 0d0
      Axnorm = 0d0
      Anorm = 0d0
      Acond = 1d0
      relres = rnorm/(beta1 + 1d-50)
      x = zzero
      w = zzero
      wl = zzero
      r1 = zzero

      resvec(1) = beta1

      do while ((iter.lt.maxiter).and.(flag.eq.flag0))

        iter = iter + 1
        betal = beta
        beta = betan

        v(:) = (1d0/beta)*r3(:)

        if (ipol.eq.1) then
         call gw_aprod_cmplx(OV1,factor,conjg(v),r3,4d0,0d0,npoles(1),n)
        else
         call gw_aprod_cmplx(OV1,factor,conjg(v),r3,2d0,0d0,npoles(1),n)
         call gw_aprod_cmplx(OV2,factor(1,2),conjg(v),r3,2d0,1d0,
     $                        npoles(2),n)
        endif
        call ga_dgop(1038,r3,n,'+')

        r3(:) = r3(:) + conjg(v(:))

        if (iter.gt.1) r3 = r3 - (beta/betal)*r1
        alpha = dot_product(v,r3)
        r3 = r3 - (alpha/beta)*r2
        r1 = r2
        r2 = r3

        betan = dsqrt(dble(dot_product(r3,r3)))
        if (iter.eq.1) then
          if (betan.eq.0d0) then
            if (alpha.eq.zzero) then
              flag = 0
              exit
            else
              flag = -1
              x = dcmplx(b(:),0d0)/alpha
              exit
            endif
          endif
        endif
        pnorm = dsqrt(betal**2 + dble(alpha*conjg(alpha)) + betan**2)

        dbar = dltan
        dlta = cs*dbar + sn*alpha
        epln = eplnn
        gbar = conjg(sn)*dbar - cs*alpha
        eplnn = sn*betan
        dltan = -cs*betan
        dlta_QLP = dlta

        gamal3 = gamal2
        gamal2 = gamal
        gamal = gama
        call symortho(gbar, dcmplx(betan,0d0), cs, sn, gama)
        gama_tmp = gama
        taul2 = tau
        taul = tau
        tau = cs*phi

        Axnorm = dsqrt(Axnorm**2 + dble(tau*conjg(tau)))

        phi = conjg(sn)*phi

        if (iter.gt.2) then
          veplnl2 = veplnl
          etal2 = etal
          etal = eta
          dlta_tmp = sr2*vepln - cr2*dlta
          veplnl = cr2*vepln + conjg(sr2)*dlta
          dlta = dlta_tmp
          eta = conjg(sr2)*gama
          gama = -cr2*gama
        endif

        if (iter.gt.1) then
          call symortho(conjg(gamal), conjg(dlta), cr1, sr1, gamal)
          gamal = conjg(gamal)
          vepln = conjg(sr1)*gama
          gama = -cr1*gama
        endif

        xnorml = xnorm
        ul4 = ul3
        ul3 = ul2
        if (iter.gt.2) ul2 = (taul2-etal2*ul4-veplnl2*ul3)/gamal2
        if (iter.gt.1) ul = (taul-etal*ul3-veplnl*ul2)/gamal
        xnorm_tmp = dznrm2(3,
     $      [dcmplx(xl2norm,0d0),ul2,ul],1)
        if (abs(gama) > pnorm*realmin .and. xnorm_tmp .lt.maxxnorm)then
          u = (tau - eta*ul2 - vepln*ul)/gama
          if (dsqrt(xnorm_tmp**2 + dble(u*conjg(u))).gt.maxxnorm )then
            u = zzero
            flag = 6
          endif
        else
          u = zzero
          flag = 9
        endif
        xl2norm = dsqrt( xl2norm**2 + dble(ul2*conjg(ul2)))
        xnorm = dznrm2(3,[dcmplx(xl2norm,0d0),ul,u],1)

        if ((Acond.lt.TranCond).and.(flag.eq.flag0).and.
     $      (QLPiter.eq.0)) then
          wl2 = wl
          wl = w
          w = (conjg(v) - epln*wl2 - dlta_QLP*wl) * (1d0/gama_tmp)
          if (xnorm .lt. maxxnorm) then
            x = x + tau*w
          else
            flag = 6
          endif
        else
          QLPiter = QLPiter + 1
          if (QLPiter.eq.1) then
            xl2 = zzero
            if (iter.gt.1) then
              if (iter.gt.3) wl2 = gamal3*wl2 + veplnl2*wl + etal*w
              if (iter.gt.2) wl = gama_QLP*wl + vepln_QLP*w
              w = gama_QLP*w
              xl2 = x - wl*ul_QLP - w*u_QLP
            endif
          endif
          if (iter.eq.1) then
            wl2 = wl
            wl  = conjg(v)*conjg(sr1)
            w   = -conjg(v)*cr1
          elseif (iter.eq.2) then
            wl2 = wl
            wl  = w*cr1 + conjg(v)*conjg(sr1)
            w   = w*sr1 - conjg(v)*cr1
          else
            wl2 = wl
            wl = w
            w = wl2*sr2 - conjg(v)*cr2
            wl2 = wl2*cr2 + conjg(v)*conjg(sr2)
            v = wl*cr1 + w *conjg(sr1)
            w = wl*sr1 - w*cr1
            wl = v
          endif
          xl2 = xl2 + wl2*ul2
          x = xl2 + wl*ul + w*u
        endif

        gamal_tmp = gamal
        call symortho(conjg(gamal_tmp),conjg(eplnn),cr2,sr2,gamal)
        gamal = conjg(gamal)

        gamal_QLP = gamal_tmp
        vepln_QLP = vepln
        gama_QLP = gama
        ul_QLP = ul
        u_QLP = u

        abs_gama = abs(gama)
        Anorml = Anorm
        Anorm = maxval( (/ Anorm, abs(gamal), abs_gama /) )
        if (iter.eq.1) then
          gmin = abs(gama)
          gminl = gmin
        else
          gminl2 = gminl
          gminl = gmin
          gmin = minval( (/ abs(gminl2), abs(gamal), abs_gama /) )
        endif
        Acondl = Acond
        Acond = Anorm/gmin
        rnorml = rnorm
        relresl = relres
        if (flag.ne.9) rnorm = abs(phi)
        relres = rnorm/(Anorm*xnorm + beta1)
        rootl = dznrm2( 2, [ gbar, dltan ], 1 )
        Arnorml = rnorml*rootl
        relAresl = rootl/Anorm
        epsx = Anorm*xnorm*eps


        if ((flag.eq.flag0) .or. (flag.eq.9)) then
          t1 = 1d0 + relres
          t2 = 1d0 + relAresl
          if (iter .ge. maxiter) flag = 8
          if (Acond .ge. Acondlim) flag = 7
          if (xnorm .ge. maxxnorm) flag = 6
          if (epsx .ge. beta1) flag = 5
          if (t2 .le. 1d0) flag = 4
          if (t1 .le. 1d0) flag = 3
          if (relAresl .le. rtol) flag = 2
          if (relres .le. rtol) flag = 1
          if (rnorm .le. rtol) flag = 1
        endif

        if ((flag.eq.2).or.(flag.eq.4).or.((flag.eq.6).and.likeLS).or.
     $      (flag.eq.7)) then
          iter = iter - 1
          Acond = Acondl
          rnorm = rnorml
          relres = relresl
        else
          resvec(iter+1) = rnorm
          Aresvec(iter) = Arnorml
        endif
      enddo

      x = x + xbak

      end subroutine

      subroutine symortho(a, b, c, s, r)
      double complex,intent(in)    :: a, b
      double complex,intent(out)   :: s, r
      double precision,intent(out) :: c

      double precision absa, absb, t
      double complex signa, signb
      double complex, parameter :: zzero = (0d0, 0d0), zone = (1d0, 0d0)

      signa = 1d0
      signb = 1d0
      
      absa = abs(a)
      absb = abs(b)
      if (absa.gt.0d0) signa = a/absa
      if (absb.gt.0d0) signb = b/absb

      if (b.eq.zzero) then
        c = 1d0
        s = (0d0, 0d0)
        r = a
        return
      endif

      if (a.eq.zzero) then
        c = 0d0
        s = (1d0, 0d0)
        r = b
        return
      endif

      if (absb .gt. absa) then
        t = absa/absb
        c = 1d0/dsqrt(1d0 + t**2)
        s = c*conjg(signb/signa)
        c = c*t
        r = b/conjg(s)
      else
        t = absb/absa
        c = 1d0/dsqrt(1d0 + t**2)
        s = c*t*conjg(signb/signa)
        r = a/c
      endif

      end subroutine

      subroutine gw_aprod_cmplx(ov,factor,x,y,scal,mult,npoles,n)
      implicit none
      integer,intent(in) :: npoles, n
      double precision, intent(in) :: scal, mult
      double precision, intent(in) :: ov(n,npoles)
      double complex, intent(in) :: factor(npoles), x(n)
      double complex, intent(out) :: y(n)

      integer :: i
      double complex :: tmp(npoles)

      call dgemm('n','n',2,npoles,n,scal,x,2,ov,n,0d0,tmp,2)

      tmp(:) = tmp(:)*factor(:)

      call dgemm('n','t',2,n,npoles,1d0,tmp,2,ov,n,mult,y,2)

      end subroutine

