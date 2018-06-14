
      subroutine vdw_DF_init()
      implicit none

c     **** Langreth kernel common block ****
      integer Na
      parameter (Na=256)
      !real*8 amin,amax,a(Na),aweights(Na),a2(Na),cos_a(Na),sin_a(Na)
      !real*8 nu(Na),nu1(Na)
      !real*8 Wab(Na,Na)
      real*8 amin,amax
      integer a(2),aweights(2),a2(2),cos_a(2),sin_a(2)
      integer nu(2),nu1(2),Wab(2)
      common /lang_phi_common/ Wab,a,aweights,a2,cos_a,sin_a,nu,nu1,
     >                         amin,amax

      return
      end

      subroutine vdw_DF_init0(Na,amin,amax,a,a2,aweights,
     >                        cos_a,sin_a,nu,nu1,Wab)
      implicit none
      integer Na
      real*8 amin,amax,a(Na),aweights(Na),a2(Na),cos_a(Na),sin_a(Na)
      real*8 nu(Na),nu1(Na),Wab(Na,Na)

      integer, parameter :: Nqs = 20

       real(dp), dimension(Nqs):: q_mesh = (/  1.00D-5, 0.0449420825586261D0, 0.0975593700991365D0, &
       0.159162633466142D0, 0.231286496836006D0, 0.315727667369529D0, 0.414589693721418D0, &
       0.530335368404141D0, 0.665848079422965D0, 0.824503639537924D0, 1.010254382520950D0, &
       1.227727621364570D0, 1.482340921174910D0, 1.780437058359530D0, 2.129442028133640D0, &
       2.538050036534580D0, 3.016440085356680D0, 3.576529545442460D0, 4.232271035198720D0, &
       5.0D0 /)


      integer i,j
      real*8  D,delta,d1,d2,fourpi
      real*8  phi0,phi5,phi9

*     **** external functions ****
      real*8   phi_value
      external phi_value

      fourpi = 16.0d0*datan(1.0d0)
      amin = 0.0d0
      amax = 64.0d0

      call GaussLegendre(datan(amin),datan(amax),Na,a,aweights)
      do i=1,Na
         a(i) = dtan(a(i))
         a2(i) = a(i)**2
         aweights(i) = aweights(i)*(1.0d0+a2(i))
         cos_a(i) = dcos(a(i))
         sin_a(i) = dsin(a(i))
      end do
      do i = 1, Na
      do j = 1, Na
         Wab(i,j) = 2.0d0 * aweights(i)*aweights(j) *  (
     >       (3.0d0-a2(i))*a(j)*cos_a(j)*sin_a(i)  +
     >       (3.0d0-a2(j))*a(i)*cos_a(i)*sin_a(j)  +
     >       (a2(i)+a2(j)-3.0d0)*sin_a(i)*sin_a(j) -
     >       3.0d0*a(i)*a(j)*cos_a(i)*cos_a(j) )   / (a(i)*a(j))
      end do
      end do


      ! The total number of phi_alpha_beta functions that have to be calculated
      Ntotal = (Nqs**2 + Nqs)/2

      dr = (r_max)/(Nr_points)
      do idx = 1, my_Nqs
      do i = 1, Nr_points

        d1 = q_mesh(indices(idx+my_start_q-1, 1)) * (dr * i)
        d2 = q_mesh(indices(idx+my_start_q-1, 2)) * (dr * i)

        phi(r_i, idx) = phi_value(d1, d2)

      end do
      end do
     



      return
      end

      subroutine vdw_DF_energy(dn,gdn)
      implicit none
      real*8 dn(*),gdn(*)
