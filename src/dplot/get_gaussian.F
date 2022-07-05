!>    @brief Evaluate a Gaussian centered on an atom at a grid point
!>
!>    We have a normalized spherical Gausian centered at atomic position
!>    RA with a width determined by the blur parameter and scale by the
!>    nuclear charge Z, and we evaluate its value at grid point ri.
!>    \f{eqnarray*}{
!>      g(r_i) &=& \frac{Z}{\sqrt{\sigma\pi}^3}
!>                 \exp{-\frac{|r_i-R_a|^2}{\sigma}}
!>    \f}
!>    where \f$\sigma\f$ is equal to blur.
!>
!>    @param[in] RA   The atom position
!>    @param[in] ri   The grid point position
!>    @param[in] blur The blurring factor of the Gaussian
!>    @param[in] Z    The nuclear charge
!>
!>    @return The function value
!>
      double precision function get_gaussian(RA,ri,blur,Z)
      implicit none
      double precision :: RA(3)
      double precision :: ri(3)
      double precision :: blur
      double precision :: Z
!
      double precision :: xAi
      double precision :: yAi
      double precision :: zAi
      double precision :: dAi
      double precision :: pi
      double precision :: A
      double precision :: E
!
      xAi = RA(1)-ri(1)
      yAi = RA(2)-ri(2)
      zAi = RA(3)-ri(3)
      dAi = dsqrt(xAi*xAi+yAi*yAi+zAi*zAi)
      pi  = dacos(-1.0d0)
      A   = Z/(dsqrt(blur*pi)**3)
      E   = -dAi*dAi/blur
      get_gaussian = A*exp(E)
      return
      end
