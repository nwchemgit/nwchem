c
C$Id$
c
      subroutine xc_kop(tol_rho,whichf,x,
     &     kalpha, dkadxa)
      implicit none
c
      character*4 whichf
      double precision tol_rho,x
      double precision kalpha,dkadxa
      if(whichf.eq.'be88') then
         call xc_kbecke88(tol_rho,x,
     &     kalpha, dkadxa)
      endif
      if(whichf.eq.'pb96') then
         call xc_kpbe96(tol_rho,x,
     &     kalpha, dkadxa)
      endif
      return
      end
      subroutine xc_kbecke88(tol_rho,x,
     &     kalpha, dkadxa)
      implicit none
c
      double precision tol_rho,x
      double precision kalpha,dkadxa
c
      double precision BETA, C
      Parameter (BETA = 0.0042D0)
      double precision g,gdenom,dgdenom,dg
      double precision arcsinh, darcsinh
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
c
c
c     Uniform electron gas constant
c
      C =  3d0*(0.75d0/acos(-1d0))**(1d0/3d0)

      if (x.gt.0d0)then
         gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
         g = 2d0*BETA*x*x / gdenom
         dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
         dg = g*(2d0/x-dgdenom/gdenom)
         
         kalpha= C + g
         dkadxa = dg
         
      else
         kalpha= C
         dkadxa = 0d0
      endif
      return
      end
      subroutine xc_kpbe96(tol_rho,x,
     &     kalpha, dkadxa)
      implicit none
c
      double precision tol_rho,x
      double precision kalpha,dkadxa
c
      double precision pi,um, uk, umk
      parameter(um=0.21951d0, uk=0.804d0, umk=um/uk)
      double precision C
      double precision forty8,deno
c
c
c     Uniform electron gas constant
c
      pi = acos(-1.d0)
      C =  3d0*(0.75d0/pi)**(1d0/3d0)

      if (x.gt.0d0)then
         forty8=1d0/((48d0*pi*pi)**(2d0/3d0))
         deno=1d0/(1d0+um*x*x*forty8/uk)
         kalpha= C * (1d0 + uk - uk *deno)
         dkadxa = C * (2d0*um*x*deno*deno*
     *        forty8)
         
      else
         kalpha= C
         dkadxa = 0d0
      endif
      return
      end
