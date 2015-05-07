C> \brief Program to demonstrate the numerical implications of chi
C>
C> In 1988 Becke (doi: 10.1103/PhysRevA.38.3098) introduced the
C> dimensionless parameter chi:
C> \f{eqnarray}{
C>    \chi &=& \frac{|\nabla\rho(r)|}{\rho(r)^{4/3}} \\\\
C>         &=& \frac{\sqrt{\gamma}{\rho^{4/3}}
C> \f}
C> Since then this parameter has been used in a number of functionals.
C> In later functionals the quantity that is often used is chi squared
C> (e.g. see Becke97, doi: 10.1063/1.475007). In actual codes, however,
C> chi squared is still often calculated by first calculating chi and
C> then squaring it. I.e. first taking the square root of gamma,
C> multiplying it with some factors, and then squaring it up again. In
C> formal mathematics there is nothing wrong with this. 
C>
C> On finite precision machines and in particular for derivatives
C> (especially second order and higher) this approach has consequences
C> for the accuracy of the derivatives. 
C>
C> This program demonstrates the impact of calculating chi and its
C> derivatives, and then calculating gamma and its derivatives from
C> chi, versus calculating the derivatives of gamma directly. Obviously
C> the second derivative of gamma wrt. itself should be 0. Calculating
C> this derivative via chi leads to errors of the order of epsilon/gamma
C> (where epsilon is the machine precision). I.e. for small gammas
C> the second order derivatives are seriously in error for no good
C> reason.
C>
C>
C> THE RECOMMENDATION: therefore is that if the functional uses only
C> even powers of chi then eliminate all explicit references to chi
C> and express the equations directly in terms of gamma instead.
C> For functionals with odd orders of chi there is no hope and we have
C> to accept that derivatives will become singular for small gammas.
C>
C> 
C> Inputs:
C>
C> This program takes two floating point numbers on one line as input.
C> The first one is the top value for gamma and the second one is
C> the lowest value for gamma. It automatically generates a number of
C> intermediate values as well. Example:
C>
C>    10.0 1.0e-21
C>
C>
C> Outputs:
C>
C> The program prints a table of pairs of lines. The top line of each
C> pair lists the row number, gamma, the first derivative of gamma wrt.
C> to itself (i.e. 1.0), and the second derivative of gamma wrt. itself
C> (i.e. 0.0). The bottom line of each pair lists the same quantities
C> but calculated via chi. On an infinite precision machine both results
C> would be identical but on a finite precision machine they are not.
C> Example:
C>
C>    100   0.100000000000000E-20   0.100000000000000E+01   0.000000000000000E+00
C>    100   0.100000000000000E-20   0.100000000000000E+01  -0.655360000000000E+05
C> 
      program numerics_chi
      implicit none
      integer mxp
      parameter(mxp = 7)
      double precision gamma(mxp)
      double precision dgamma(mxp)
      double precision d2gamma(mxp)
      double precision d3gamma(mxp)
      double precision chi(mxp)
      double precision dchi(mxp)
      double precision d2chi(mxp)
      double precision d3chi(mxp)
      double precision gamma2(mxp)
      double precision dgamma2(mxp)
      double precision d2gamma2(mxp)
      double precision d3gamma2(mxp)

      integer ii

      double precision top, bot, fact
      read(*,*)top,bot
      fact = exp(-log(top/bot)/(mxp-1))

      gamma(1) = top
      do ii = 2, mxp
        gamma(ii) = gamma(ii-1)*fact
      enddo
c
      do ii = 1, mxp
        chi(ii) = sqrt(gamma(ii))
      enddo
c
      do ii = 1, mxp
        dgamma(ii) = 1.0d0
      enddo
c
      do ii = 1, mxp
        dchi(ii) = 0.5d0/sqrt(gamma(ii))
      enddo
c
      do ii = 1, mxp
        d2gamma(ii) = 0.0d0
      enddo
c
      do ii = 1, mxp
        d2chi(ii) = -0.25d0/(gamma(ii)**1.5d0)
      enddo
c
      do ii = 1, mxp
        gamma2(ii) = chi(ii)*chi(ii)
      enddo
c
      do ii = 1, mxp
        dgamma2(ii) = dchi(ii)*chi(ii) + chi(ii)*dchi(ii)
      enddo
c
      do ii = 1, mxp
        d2gamma2(ii) = d2chi(ii)*chi(ii) + 2*dchi(ii)*dchi(ii)
     +               + chi(ii)*d2chi(ii)
      enddo
c
      do ii = 1, mxp
        d3chi(ii) = 1.5d0*0.25d0/(gamma(ii)**2.5d0)
      enddo
c
      do ii = 1, mxp
        d3gamma(ii) = 0.0d0
      enddo
c
      do ii = 1, mxp
        d3gamma2(ii) = d3chi(ii)*chi(ii) + 3*d2chi(ii)*dchi(ii)
     +               + 3*dchi(ii)*d2chi(ii) + chi(ii)*d3chi(ii)
      enddo
c
      do ii = 1, mxp
        write(*,'(i6,4e24.15)')ii,gamma(ii),dgamma(ii),
     +                            d2gamma(ii),d3gamma(ii)
        write(*,'(i6,4e24.15)')ii,gamma2(ii),dgamma2(ii),
     +                            d2gamma2(ii),d3gamma2(ii)
        write(*,*)
      enddo
c
      end
