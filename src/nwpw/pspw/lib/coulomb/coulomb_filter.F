      

*     **********************************************
*     *                                            *
*     *              nwpw_pochhammer               *
*     *                                            *
*     ********************************************** 
      real*8 function nwpw_pochhammer(a,n)
      implicit none
      real*8 a
      integer n

*     **** local variables ****
      integer i
      real*8 fac

      fac = 1.0d0
      do i=0,n-1
         fac = fac * (a+dble(i))
      end do

      nwpw_pochhammer = fac
      return 
      end

*     **********************************************
*     *                                            *
*     *                nwpw_Qfac                   *
*     *                                            *
*     ********************************************** 
      real*8 function nwpw_Qfac(k,m,n)
      implicit none
      integer k,m,n

*     **** local variables ****
      real*8 nfac,num,den

*     **** external functions ****
      real*8   nwpw_pochhammer
      external nwpw_pochhammer

      if (abs(m).gt.0) then
         nfac = 2.0d0*nwpw_pochhammer(1.0d0,n-2*k-abs(m)+1)
     >               /nwpw_pochhammer(1.0d0,n-2*k+abs(m)+1)
      else
         nfac = 1.0d0
      end if
      num = nwpw_pochhammer(0.5d0,k)
     >     *nwpw_pochhammer(1.0d0,n-k) 
     >     *(1.0d0+2.0d0*n - 4.0d0*k)
      den = nwpw_pochhammer(1.0d0,k)*nwpw_pochhammer(1.5d0,n-k)

      nwpw_Qfac = nfac*num/den
      return
      end 


*     ***************************************************
*     *                                                 *
*     *                 nwpw_plgndrs                    *
*     *                                                 *
*     ***************************************************

*  This routine calculates all the associated Legendre
*  polynomials thru lmax at the value x.
*
*  P_lm(x) = (1-x^2)^(m/2) * d^m/dx^m Pl(x)
*
*  and 
*
*  P_l(x) = 1/(2^l*l!) * d^l/dx^l (x^2-1)^l
*
*  In this routine the following relationships are used to calculate P_lm(x)
*
*   P_mm(x)    = (2*m-1)!! (1-x^2)^(m/2)
*   P_m+1,m(x) = x*(2*m+1)*P_mm(x)
*
* and
*
*    (l-m)*P_lm(x) = x*(2*l-1)*P_(l-1),m(x) - (l+m-1)*P_(l-2),m(x)

      subroutine nwpw_plgndrs(lmax,x,Plm)
      implicit none
      integer lmax
      real*8 x,Plm(*)

*     **** local variables ****
      integer i,l,m, lm,lm0,lm1,lm2,lmi,lmi0,lmi1,lmi2
      real*8  fact,pmm,somx2

      !*** compute P_mm(x) thru lmax ***
      somx2 = ((1.0d0-x)*(1.0d0+x))
      if (somx2.lt.0.0d0) somx2=0.0d0
      somx2 = dsqrt(somx2)
      pmm   = 1.0d0
      fact  = 1.0d0
      lm = 2
      Plm(0+1) = 1.0d0
      do l=1,lmax
         pmm = pmm*fact*somx2
         fact = fact + 2.0d0
         Plm(lm+1)=pmm
         lm = lm + (l+2)
      end do

      !*** compute P_(m+1),m(x) thru lmax ***
      lm  = 0
      lm2 = 1
      do m=0,(lmax-1)
         Plm(lm2+1) = x*(2*m+1)*Plm(lm+1)
         lm  = lm  + (m+2)
         lm2 = lm2 + (m+3)
      end do

      !*** compute the rest of the P_lm's ***
      lmi0 = 0
      lmi1 = 1
      lmi  = 3
      do m=0,lmax-2
         lm0 = lmi0
         lm1 = lmi1
         lm  = lmi
         do l=m+2,lmax
            Plm(lm+1) = ( x*dble(2*l-1)*Plm(lm1+1)
     >                  - dble(l+m-1)*Plm(lm0+1))/dble(l-m)
            lm0 = lm1 
            lm1 = lm
            lm  = lm + (l+1)
         end do
         lmi0 = lmi0 + (m+2)
         lmi1 = lmi1 + (m+3)
         lmi  = lmi  + (m+4)
      end do

      return
      end


*     ***************************************************
*     *                                                 *
*     *                 nwpw_Tesserals                  *
*     *                                                 *
*     ***************************************************/

*  this routine calculates Tesseral harmonics thru lmax, where
*
*  Tlm(theta,phi) =                            Pl|m|(cos(theta))            for m =0
*  Tlm(theta,phi) = sqrt(2*(l-|m|)!/(l+|m|)) * Pl|m|(cos(theta))*cos(m*phi) for m> 0
*  Tlm(theta,phi) = sqrt(2*(l-|m|)!/(l+|m|)) * Pl|m|(cos(theta))*sin(m*phi) for m<0
*
*  Tlm(theta,phi) = Pl|m|(cos(theta))            for m =0
*  Tlm(theta,phi) = Pl|m|(cos(theta))*cos(m*phi) for m> 0
*  Tlm(theta,phi) = Pl|m|(cos(theta))*sin(m*phi) for m<0
*   
*  Plm(x) = (1-x^2)^(m/2) d^m/dx^m Pl(x)
*
*  where -l<= m <= l  and -1 <= x <= 1
*
*  Note that the (-1)^m factor is not included
*
      subroutine nwpw_Tesserals(lmax,costheta,cosphi,sinphi,Plm,Tlm)
      implicit none
      integer lmax
      real*8 costheta,cosphi,sinphi
      real*8 Plm(*),Tlm(*)

*     **** local variables ****
      integer l,m
      integer lm,lmc,lms,lmi,lmci,lmsi
      real*8  cs0,cs1,cs2,sn0,sn1,sn2,fact

      call nwpw_plgndrs(lmax,costheta,Plm)

      !*** m==0 copy ***
      lm  = 0
      lmc = 0
      do l=0,lmax
         Tlm(lmc+1) = Plm(lm+1)
         lm  = lm  + (l+1)
         lmc = lmc + (2*l+2)
      end do

      !*** m>0 copy ***
      cs0 = 0.0d0
      sn0 = 0.0d0
      cs1 = 1.0d0
      sn1 = 0.0d0
      cs2 = cosphi
      sn2 = sinphi
      lmi  = 2
      lmci = 3
      lmsi = 1
      do m=1,l
         lm  = lmi
         lmc = lmci
         lms = lmsi
         do l=m,lmax
            Tlm(lmc+1) = cs2*Plm(lm+1)
            Tlm(lms+1) = sn2*Plm(lm+1)
            lm  = lm  + (l+1)
            lmc = lmc + (2*l+2)
            lms = lms + (2*l+2)
         end do
         cs0 = cs1
         cs1 = cs2
         cs2 = (2.0d0*cosphi)*cs1 - cs0
         sn0 = sn1
         sn1 = sn2
         sn2 = (2.0d0*cosphi)*sn1 - sn0
         lmi  = lmi  + (m+2)
         lmci = lmci + (2*m+3)
         lmsi = lmsi + (2*m+1)
      end do
      return
      end


*     **********************************************
*     *                                            *
*     *              f3_to_theta_min               *
*     *                                            *
*     ********************************************** 
      real*8 function f3_to_theta_min(f3)
      implicit none
      real*8 f3
 
*     **** local variables ****
      real*8 pi,theta_min

      pi = 4.0*datan(1.0d0)
      theta_min = 0.0d0

      if (dabs(f3).lt.1.0d-6) then
         theta_min = 0.0d0
      else if (f3.gt.0.0d0)  then
         theta_min = 0.0d0
      else 
         theta_min = 0.5d0*pi
      end if
    
      f3_to_theta_min = theta_min 
      return
      end

*     **********************************************
*     *                                            *
*     *              f3_to_theta_max               *
*     *                                            *
*     ********************************************** 
      real*8 function f3_to_theta_max(f3)
      implicit none
      real*8 f3
 
*     **** local variables ****
      real*8 pi,theta_max

      pi = 4.0*datan(1.0d0)
      theta_max = 0.0d0

      if (dabs(f3).lt.1.0d-6) then
         theta_max = pi
      else if (f3.gt.0.0d0)  then  
         theta_max = 0.5d0*pi
      else                 
         theta_max = pi
      end if
    
      f3_to_theta_max = theta_max 
      return
      end


*     **********************************************
*     *                                            *
*     *              f1f2_to_phi_min               *
*     *                                            *
*     ********************************************** 
      real*8 function f1f2_to_phi_min(f1,f2)
      implicit none
      real*8 f1,f2

*     **** local variables ****
      real*8 phi_min,pi

      pi = 4.0d0*datan(1.0d0)

      phi_min = 0.0d0
      if ((dabs(f1).lt.1.0d-6).and.(dabs(f2).lt.1.0d-6)) then
         phi_min = 0.0d0
      else if ((f1.gt.0.0d0).and.(dabs(f2).lt.1.0d-6)) then
         phi_min = -0.5d0*pi
      else if ((f1.lt.0.0d0).and.(dabs(f2).lt.1.0d-6)) then
         phi_min =  0.5d0*pi
      else if ((dabs(f1).lt.1.0d-6).and.(f2.gt.0.0)) then
         phi_min =  0.0d0
      else if ((f1.gt.0.0d0).and.(f2.gt.0.0d0)) then
         phi_min =  0.0d0
      else if ((f1.lt.0.0d0).and.(f2.gt.0.0d0)) then
         phi_min =  0.5d0*pi
      else if ((dabs(f1).lt.1.0d-6).and.(f2.lt.0.0d0)) then
         phi_min =  pi
      else if ((f1.gt.0.0d0).and.(f2.lt.0.0d0)) then
         phi_min =  1.5d0*pi
      else if ((f1.lt.0.0d0).and.(f2.lt.0.0d0)) then
         phi_min =  pi
      end if

      f1f2_to_phi_min = phi_min
      return
      end

*     **********************************************
*     *                                            *
*     *              f1f2_to_phi_max               *
*     *                                            *
*     ********************************************** 
      real*8 function f1f2_to_phi_max(f1,f2)
      implicit none
      real*8 f1,f2

*     **** local variables ****
      real*8 phi_max,phi,pi

      pi=4.0d0*datan(1.0d0)

      phi_max = 0.0d0
      if ((dabs(f1).lt.1.0d-6).and.(dabs(f2).lt.1.0d-6)) then
         phi_max = 2.0d0*pi
      else if ((f1.gt.0.0d0).and.(dabs(f2).lt.1.0d-6)) then
         phi_max =  0.5d0*pi
      else if ((f1.lt.0.0d0).and.(dabs(f2).lt.1.0d-6)) then
         phi_max =  1.5d0*pi
      else if ((dabs(f1).lt.1.0d-6).and.(f2.gt.0.0d0)) then
         phi_max =  pi
      else if ((f1.gt.0.0d0).and.(f2.gt.0.0d0)) then
         phi_max =  0.5d0*pi
      else if ((f1.lt.0.0d0).and.(f2.gt.0.0d0)) then
         phi_max =  pi
      else if ((dabs(f1).lt.1.0d-6).and.(f2.lt.0.0d0)) then
         phi_max =  2.0d0*pi
      else if ((f1.gt.0.0d0).and.(f2.lt.0.0d0)) then
         phi_max =  2.0d0*pi
      else if ((f1.lt.0.0d0).and.(f2.lt.0.0d0)) then
         phi_max =  1.5d0*pi
      end if

      f1f2_to_phi_max = phi_max
      return
      end 


*     ***************************************************
*     *                                                 *
*     *    setup_integrate_oneoverG2_BZ_Gengenbauer2    *
*     *                                                 *
*     ***************************************************
      subroutine setup_integrate_oneoverG2_BZ_Gengenbauer2(Nt3,
     >                gammasdxdydz,lk,nmax,Mpack,
     >                QQ,Plm,Tlm)
      implicit none
      integer Nt3
      real*8 gammasdxdydz(*),lk(3,*)
      integer nmax
      real*8  Mpack(*)
      real*8  QQ(*),Plm(*),Tlm(*)

*     **** local variables ****
      integer i,n,k,m,l,lm,indx
      real*8  gx,gy,gz,g,gxy,gn,costheta,cosphi,sinphi

*     **** external functions ****
      real*8   nwpw_Qfac
      external nwpw_Qfac

      indx = 0
      do n=0,nmax
         do k=0,n/2
            do m=(-n+2*k),(n-2*k)
               Mpack(indx+1) = 0.0d0
               QQ(indx+1)    = nwpw_Qfac(k,m,n)
               indx = indx + 1
            end do
         end do
      end do

      do i=1,Nt3
         gx = lk(1,i)
         gy = lk(2,i)
         gz = lk(3,i)
         g   = dsqrt(gx*gx + gy*gy + gz*gz)
         gxy = dsqrt(gx*gx + gy*gy)
         if (g.gt.1.0d-9) then
            costheta = gz/g
            if (gxy.lt.1.0d-9) then
               cosphi = 1.0
               sinphi = 0.0
            else
               cosphi = gx/gxy
               sinphi = gy/gxy
            end if

            call nwpw_Tesserals(nmax,costheta,cosphi,sinphi,Plm,Tlm)

            indx = 0
            gn = 1.0d0
            do n=0,nmax
               do k=0,n/2
                  l = n-2*k
                  do m=(-n+2*k),(n-2*k)
                     lm = (l+1)*(l+1) - l + m - 1
                     Mpack(indx+1) = Mpack(indx+1) 
     >                             + QQ(indx+1)*gn*Tlm(lm+1)
     >                              *gammasdxdydz(i)
                     indx = indx + 1
                  end do
               end do
               gn = gn*g
            end do
         end if
      end do
      return
      end 

*     ***************************************************
*     *                                                 *
*     *        integrate_oneoverG2_BZ_Gengenbauer2      *
*     *                                                 *
*     ***************************************************

      real*8 function integrate_oneoverG2_BZ_Gengenbauer2(nmax,Mpack,
     >                    gx,gy,gz,Plm,Tlm)
      implicit none
      integer nmax
      real*8 Mpack(*),gx,gy,gz
      real*8 Plm(*),Tlm(*)

*     **** local variables ****
      integer indx,n,k,m,l,lm
      real*8  sum0,sum,costheta,cosphi,sinphi
      real*8  g,gxy,gn

      g   = dsqrt(gx*gx + gy*gy + gz*gz)
      gxy = dsqrt(gx*gx + gy*gy)
      gn = g*g

      costheta = gz/g
      if (gxy.lt.1.0d-9) then
         cosphi = 1.0d0
         sinphi = 0.0d0
      else
         cosphi = gx/gxy
         sinphi = gy/gxy
      end if
   
      call nwpw_Tesserals(nmax,costheta,cosphi,sinphi,Plm,Tlm)
      !if ((dabs(gx).lt.1.0d-9).and.(dabs(gz).lt.1.0d-9)) then
      !   call errquit("gx==0 and gz==0",0,0)
      !end if

      indx = 0
      sum  = 0.0d0
      do n=0,nmax
         sum0 = 0.0d0
         do k=0,(n/2)
            l = n-2*k
            do m=(-n+2*k),(n-2*k)
               lm = (l+1)*(l+1) -  l + m - 1
               sum0 = sum0 + Mpack(indx+1)*Tlm(lm+1)
               indx = indx + 1
            end do
         end do
         sum = sum + sum0/gn
         gn = gn*g
      end do

      integrate_oneoverG2_BZ_Gengenbauer2 = sum
      return 
      end


*     **********************************************
*     *                                            *
*     *    integrate_spherical_oneoveralphaG2_BZ   *
*     *                                            *
*     **********************************************

      real*8 function integrate_spherical_oneoveralphaG2_BZ(af1,af2,af3,
     >                                     Nr,r0,r1,
     >                                     Ntheta,theta0,theta1,
     >                                     Nphi,phi0,phi1,
     >                                     unitg,alpha)
      implicit none
      real*8 af1,af2,af3
      integer Nr
      real*8  r0,r1
      integer Ntheta
      real*8  theta0,theta1
      integer Nphi
      real*8  phi0,phi1
      real*8 unitg(3,3),alpha

*     **** local variables ****
      integer ir,itheta,iphi
      real*8  dr,dtheta,dphi,r,theta,phi,dx,dy,dz
      real*8  gamma1,gamma2,gamma3
      real*8  tgamma1,tgamma2,tgamma3,tR2,tlkx,tlky,tlkz,w
      real*8  fourpi,asum

      fourpi = 16.0d0*datan(1.0d0)

      dr     = (r1-r0)        /dble(Nr)
      dtheta = (theta1-theta0)/dble(Ntheta)
      dphi   = (phi1-phi0)    /dble(Nphi)

      asum = 0.0d0
      do iphi=0,Nphi
      do itheta=0,Ntheta
      do ir=0,Nr
         r     = r0     + ir*dr
         theta = theta0 + itheta*dtheta
         phi   = phi0   + iphi*dphi
         gamma1 = r*dcos(phi)*dsin(theta)
         gamma2 = r*dsin(phi)*dsin(theta)
         gamma3 = r*dcos(theta)

         tgamma1 = dcos(phi)*dsin(theta)
         tgamma2 = dsin(phi)*dsin(theta)
         tgamma3 = dcos(theta)
         tlkx=(tgamma1*unitg(1,1)+tgamma2*unitg(1,2)+tgamma3*unitg(1,3))
         tlky=(tgamma1*unitg(2,1)+tgamma2*unitg(2,2)+tgamma3*unitg(2,3))
         tlkz=(tgamma1*unitg(3,1)+tgamma2*unitg(3,2)+tgamma3*unitg(3,3))
         tR2 = tlkx*tlkx + tlky*tlky + tlkz*tlkz

         if (mod(ir,2).eq.0) then
            dx = 2.0d0*dr/3.0d0
         else
            dx = 4.0d0*dr/3.0d0
         end if

         if (mod(itheta,2).eq.0) then
            dy = 2.0d0*dtheta/3.0d0
         else 
            dy = 4.0d0*dtheta/3.0d0
         end if

         if (mod(iphi,2).eq.0) then 
            dz = 2.0d0*dphi/3.0d0
         else  
            dz = 4.0d0*dphi/3.0d0
         end if

         if ((ir.eq.0)    .or.(ir.eq.Nr))         dx = dr/3.0d0
         if ((itheta.eq.0).or.(itheta.eq.Ntheta)) dy = dtheta/3.0d0
         if ((iphi.eq.0)  .or.(iphi.eq.Nphi))     dz = dphi/3.0d0

         w =   ((1.0d0-af1)*(1.0d0-dabs(gamma1)) + af1*dabs(gamma1))
         w = w*((1.0d0-af2)*(1.0d0-dabs(gamma2)) + af2*dabs(gamma2))
         w = w*((1.0d0-af3)*(1.0d0-dabs(gamma3)) + af3*dabs(gamma3))

         asum = asum + fourpi*(dexp(-alpha*r*r*tR2)/(tR2))
     >                     *w*dsin(theta)*dx*dy*dz
      end do
      end do
      end do

      integrate_spherical_oneoveralphaG2_BZ = asum
      return
      end 

*     **********************************************
*     *                                            *
*     *         integrate_oneoveralphaG2_BZ        *
*     *                                            *
*     **********************************************

      real*8 function integrate_oneoveralphaG2_BZ(Nt3,
     >                             gammasdxdydz,lk,
     >                             alpha,gx,gy,gz)
      implicit none
      integer Nt3
      real*8  gammasdxdydz(*),lk(3,*)
      real*8  alpha,gx,gy,gz

*     **** local variables ****
      integer i
      real*8 glkx,glky,glkz,gg
      real*8 fourpi,sum

      sum = 0.0d0
      do i=1,Nt3
         glkx = (gx+lk(1,i))
         glky = (gy+lk(2,i))
         glkz = (gz+lk(3,i))
         gg = glkx*glkx + glky*glky + glkz*glkz
         if (gg.lt.1.0d-9) then
            sum = sum + gammasdxdydz(i)*alpha
         else
            sum = sum + gammasdxdydz(i)*(1.0d0-dexp(-alpha*gg))/gg
         end if
      end do

      integrate_oneoveralphaG2_BZ = sum
      return
      end 



*     **********************************************
*     *                                            *
*     *          coulomb_filter_setup              *
*     *                                            *
*     ********************************************** 
      subroutine coulomb_filter_setup(N1,a1,b1,
     >                                N2,a2,b2,
     >                                N3,a3,b3,
     >                                unitg,gammasdxdydz,lk)
      implicit none
      integer N1
      real*8  a1,b1
      integer N2
      real*8  a2,b2
      integer N3
      real*8  a3,b3
      real*8  unitg(3,3)
      real*8  gammasdxdydz(*)
      real*8  lk(3,*)

*     **** local variables ****
      real*8 dx,dy,dz
      real*8 dgamma1,dgamma2,dgamma3,gamma1,gamma2,gamma3
      real*8 fourpi
      integer i,i1,i2,i3

      fourpi = 16.0d0*datan(1.0d0)

      dgamma1 = (b1-a1)/dble(N1)
      dgamma2 = (b2-a2)/dble(N2)
      dgamma3 = (b3-a3)/dble(N3)

      i = 1
      do i3=0,N3
      do i2=0,N2
      do i1=0,N1
         gamma1 = a1 + i1*dgamma1
         gamma2 = a2 + i2*dgamma2
         gamma3 = a3 + i3*dgamma3
         if (mod(i1,2).eq.0) then
            dx = 2.0d0*dgamma1/3.0d0
         else  
            dx = 4.0d0*dgamma1/3.0d0
         end if

         if (mod(i2,2).eq.0) then
            dy = 2.0d0*dgamma2/3.0d0
         else
            dy = 4.0d0*dgamma2/3.0d0
         end if

         if (mod(i3,2).eq.0) then
            dz = 2.0d0*dgamma3/3.0d0
         else
            dz = 4.0d0*dgamma3/3.0d0
         end if

         if ((i1.eq.0).or.(i1.eq.N1)) dx = dgamma1/3.0d0
         if ((i2.eq.0).or.(i2.eq.N2)) dy = dgamma2/3.0d0
         if ((i3.eq.0).or.(i3.eq.N3)) dz = dgamma3/3.0d0

         gammasdxdydz(i) = fourpi*(1.0d0-dabs(gamma1))
     >                           *(1.0d0-dabs(gamma2))
     >                           *(1.0d0-dabs(gamma3))*dx*dy*dz
         lk(1,i) = gamma1*unitg(1,1) 
     >           + gamma2*unitg(1,2) 
     >           + gamma3*unitg(1,3)

         lk(2,i) = gamma1*unitg(2,1) 
     >           + gamma2*unitg(2,2) 
     >           + gamma3*unitg(2,3)

         lk(3,i) = gamma1*unitg(3,1)
     >           + gamma2*unitg(3,2)
     >           + gamma3*unitg(3,3)
         i = i + 1
      end do
      end do
      end do

      return 
      end




*     **********************************************
*     *                                            *
*     *              coulomb_filter                *
*     *                                            *
*     ********************************************** 
*
*  This code computes the filtered coulomb potential using
*  filon integration strategy given in "An Accurate Integration 
*  Strategy for Calculating Exact Exchange in Periodic Boundary 
*  Conditions: A Plane-Wave DFT Implementation."  
*
*  The I1,I2,and I3 integrals are computed using Simpson integration.
*
*  Warning: While this routine should be good enough for testing, it 
*  needs to be restructured to improve parallelization and 
*  the use of memory. 
*
      subroutine coulomb_filter(unita,unitg,nfft3d,G,Vfilter)
      implicit none
      real*8 unita(3,3),unitg(3,3)
      integer nfft3d
      real*8 G(nfft3d,3)
      real*8 Vfilter(nfft3d)

#include "stdio.fh"
#include "util.fh"
#include "bafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,found_datafile,oprint
      integer Nt1,Nt2,Nt3,Ntt3,Ns1,Ns2,Ns3
      parameter (Nt1=228, Nt2=228, Nt3=228,Ntt3=(Nt1+1)*(Nt2+1)*(Nt3+1))
      parameter (Ns1=432, Ns2=432, Ns3=432)

      character*255 datafile

      integer nmax,i,j,l,ifound,ngrid(3),tngrid(3)
      integer gammasdxdydz(2),lk(2),tmp(2),tmp0(2)
      integer Mpack(2),QQ(2),Plm(2),Tlm(2)
      real*8 gx,gy,gz,f1,f2,f3
      real*8 rcut,scut,scut1,salpha,alpha,alpha0
      real*8 twopi,oneovertwopi
      real*8 tunita(3,3),tunitg(3,3)

*     **** external functions ****
      real*8   integrate_oneoveralphaG2_BZ
      external integrate_oneoveralphaG2_BZ
      real*8   f3_to_theta_min,f3_to_theta_max
      external f3_to_theta_min,f3_to_theta_max
      real*8   f1f2_to_phi_min,f1f2_to_phi_max
      external f1f2_to_phi_min,f1f2_to_phi_max
      real*8   integrate_spherical_oneoveralphaG2_BZ
      external integrate_spherical_oneoveralphaG2_BZ
      real*8   integrate_oneoverG2_BZ_Gengenbauer2
      external integrate_oneoverG2_BZ_Gengenbauer2


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      call D3dB_nx(1,ngrid(1))
      call D3dB_ny(1,ngrid(2))
      call D3dB_nz(1,ngrid(3))

*     *******************************
*     **** check for kernel file ****
*     *******************************
      call util_file_name_noprefix('filtered_kernel.dat',
     >                                .false.,
     >                                .false.,
     >                                datafile)
      if (taskid.eq.MASTER) then
         ifound = 0
         inquire(file=datafile,exist=found_datafile)
         if (found_datafile) ifound = 1
      end if
      call Parallel_Brdcst_ivalue(MASTER,ifound)
      found_datafile = .false.
      if (ifound.eq.1) found_datafile=.true.

      if (found_datafile) then
         l = index(datafile,' ') -1
         ifound = 0
         if (taskid.eq.MASTER) then
            call openfile(5,datafile,l,'r',l)
            call iread(5,tngrid,3)
            call dread(5,tunita,9)
            call dread(5,tunitg,9)
            call closefile(5)
            found_datafile = ((tngrid(1).eq.ngrid(1)).and.
     >                        (tngrid(2).eq.ngrid(2)).and.
     >                        (tngrid(3).eq.ngrid(3)))
            do j=1,3
               do i=1,3
                  found_datafile = found_datafile.and.
     >                 (dabs(tunita(i,j)-unita(i,j)).lt.1.0d-6).and.
     >                 (dabs(tunitg(i,j)-unitg(i,j)).lt.1.0d-6)
               end do
            end do
            if (found_datafile) ifound = 1
         end if
         call Parallel_Brdcst_ivalue(MASTER,ifound)
         found_datafile = .false.
         if (ifound.eq.1) found_datafile=.true.

      end if


*     *******************************
*     **** read kernel from file ****
*     *******************************
      if (found_datafile) then
         if (oprint) then
            l = index(datafile,' ') - 1
            write(luout,*)
            write(luout,*) "Reading filtered kernel filename:",
     >                     datafile(1:l)

         end if
         value = BA_push_get(mt_dbl,nfft3d,'tmp',tmp(2),tmp(1))
         if (.not.value) 
     >      call errquit('coulomb_filter: out of stack memory',9,MA_ERR)
         if (taskid.eq.MASTER) then
            call openfile(5,datafile,l,'r',l)
            call iread(5,tngrid,3)
            call dread(5,tunita,9)
            call dread(5,tunitg,9)
         end if
         call D3dB_t_read(1,5,Vfilter,dbl_mb(tmp(1)),-1)
         if (taskid.eq.MASTER) then
            call closefile(5)
         end if
         value = BA_pop_stack(tmp(2))
         if (.not.value) 
     >      call errquit('coulomb_filter:popping stack memory',9,MA_ERR)
         !write(*,*) "Vfilter=",(Vfilter(i),i=1,10)


*     *******************************
*     **** generating kernel     ****
*     *******************************
      else

         if (oprint) then
            l = index(datafile,' ') - 1
            write(luout,*)
            write(luout,*) "Generating filtered kernel filename:",
     >                     datafile(1:l)

         end if

         twopi = 8.0d0*datan(1.0d0)
         oneovertwopi = 1.0d0/twopi

         nmax = 10

*        **** allocate memory ****
         value = BA_push_get(mt_dbl,3*Ntt3,'gammasdxdydz',
     >                    gammasdxdydz(2),gammasdxdydz(1))
         value = value.and.
     >        BA_push_get(mt_dbl,3*Ntt3,'lk',
     >                    lk(2),lk(1))
         value = value.and.
     >        BA_push_get(mt_dbl,5*(nmax+1)*(nmax+1),'Mpack',
     >                    Mpack(2),Mpack(1))
         value = value.and.
     >        BA_push_get(mt_dbl,5*(nmax+1)*(nmax+1),'QQ',
     >                    QQ(2),QQ(1))
         value = value.and.
     >        BA_push_get(mt_dbl,5*(nmax+1)*(nmax+1),'Plm',
     >                    Plm(2),Plm(1))
         value = value.and.
     >        BA_push_get(mt_dbl,5*(nmax+1)*(nmax+1),'Tlm',
     >                    Tlm(2),Tlm(1))
         if (.not.value) 
     >   call errquit('coulomb_filter: out of stack memory',0,MA_ERR)


         call coulomb_filter_setup(Nt1,-1.0d0,1.0d0,
     >                          Nt2,-1.0d0,1.0d0,
     >                          Nt3,-1.0d0,1.0d0,
     >                          unitg,dbl_mb(gammasdxdydz(1)),
     >                          dbl_mb(lk(1)))
         call setup_integrate_oneoverG2_BZ_Gengenbauer2(Ntt3,
     >                    dbl_mb(gammasdxdydz(1)),
     >                    dbl_mb(lk(1)),
     >                    nmax,dbl_mb(Mpack(1)),
     >                    dbl_mb(QQ(1)),dbl_mb(Plm(1)),dbl_mb(Tlm(1)))

         scut  = dsqrt( unitg(1,1)*unitg(1,1) 
     >                + unitg(2,1)*unitg(2,1) 
     >                + unitg(3,1)*unitg(3,1))
         scut1 = dsqrt( unitg(1,2)*unitg(1,2)
     >                + unitg(2,2)*unitg(2,2)
     >                + unitg(3,2)*unitg(3,2))
         if (scut1.lt.scut) scut = scut1
         scut1 = dsqrt( unitg(1,3)*unitg(1,3)
     >                + unitg(2,3)*unitg(2,3)
     >                + unitg(3,3)*unitg(3,3))
         if (scut1.lt.scut) scut = scut1
         rcut = 0.3d0*scut
         !rcut = 0.2d0*scut
         salpha = 1.0d0/rcut
         alpha  = salpha*salpha

         do i=1,nfft3d
            gx = G(i,1)
            gy = G(i,2)
            gz = G(i,3)
         
            f1=oneovertwopi*(unita(1,1)*gx+unita(1,2)*gy+unita(1,3)*gz)
            f2=oneovertwopi*(unita(2,1)*gx+unita(2,2)*gy+unita(2,3)*gz)
            f3=oneovertwopi*(unita(3,1)*gx+unita(3,2)*gy+unita(3,3)*gz)
            if ((dabs(f1).lt.2.0d0).and.
     >          (dabs(f2).lt.2.0d0).and.
     >          (dabs(f3).lt.2.0)) then
               Vfilter(i) = integrate_oneoveralphaG2_BZ(Ntt3,
     >                            dbl_mb(gammasdxdydz(1)),dbl_mb(lk(1)),
     >                            alpha,gx,gy,gz) 
               alpha0 = integrate_spherical_oneoveralphaG2_BZ(
     >                  dabs(f1),dabs(f2),dabs(f3),
     >                  Ns1,0.0d0,1.0d0,
     >                  Ns2,f3_to_theta_min(f3),f3_to_theta_max(f3),
     >                  Ns3,f1f2_to_phi_min(f1,f2),
     >                  f1f2_to_phi_max(f1,f2),
     >                  unitg,alpha) 
               Vfilter(i) = Vfilter(i) + alpha0
            else
               Vfilter(i) = integrate_oneoverG2_BZ_Gengenbauer2(nmax,
     >                            dbl_mb(Mpack(1)),gx,gy,gz,
     >                            dbl_mb(Plm(1)),dbl_mb(Tlm(1)))
            end if
         end do
         value = BA_pop_stack(Tlm(2))
         value = value.and.BA_pop_stack(Plm(2))
         value = value.and.BA_pop_stack(QQ(2))
         value = value.and.BA_pop_stack(Mpack(2))
         value = value.and.BA_pop_stack(lk(2))
         value = value.and.BA_pop_stack(gammasdxdydz(2))
         if (.not.value) 
     >      call errquit('coulomb_filter:popping stack memory',0,MA_ERR)


*        *****************************
*        **** save kernel to file ****
*        *****************************
         l = index(datafile,' ') -1
         value = BA_push_get(mt_dbl,nfft3d,'tmp',tmp(2),tmp(1))
         value = value.and.
     >           BA_push_get(mt_dbl,nfft3d,'tmp',tmp0(2),tmp0(1))
         if (.not.value)
     >      call errquit('coulomb_filter: out of stack',99,MA_ERR)
         if (taskid.eq.MASTER) then
            call openfile(6,datafile,l,'w',l)
            call iwrite(6,ngrid,3)
            call dwrite(6,unita,9)
            call dwrite(6,unitg,9)
         end if
         call D3dB_t_Copy(1,Vfilter,dbl_mb(tmp0(1)))
         call D3dB_t_write(1,6,dbl_mb(tmp0(1)),dbl_mb(tmp(1)),0)
         if (taskid.eq.MASTER) then
            call closefile(6)
         end if
         value =           BA_pop_stack(tmp0(2))
         value = value.and.BA_pop_stack(tmp(2))
         if (.not.value)
     >      call errquit('coulomb_filter:popping stack',99,MA_ERR)

         !write(*,*) "Vfilter=",(Vfilter(i),i=1,10)
      end if

      return
      end
