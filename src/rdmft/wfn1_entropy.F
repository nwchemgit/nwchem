c-----------------------------------------------------------------------
c
      subroutine wfn1_entropy1(p,s,ds)
      implicit none
C>
C> \brief Calculate the 1-electron entropy associated with an 
C> orbital and its derivative
C>
C> The entropy of a natural orbital with occupation number \f$p\f$ is
C> given by \f$S^1 = p\ln(p)+(1-p)\ln(1-p)\f$. As a result the
C> derivative is given by
C> \f$\frac{\mathrm{d}S_1}{\mathrm{d}p}=\ln(p)-\ln(1-p)\f$.
C> 
      double precision, intent(in)  :: p  ! The occupation number
      double precision, intent(out) :: s  ! The entropy
      double precision, intent(out) :: ds ! The entropy derivative
c
      s  = p*log(p) + (1.0d0-p)*log(1.0d0-p)
      ds = log(p) - log(1.0d0-p)
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_entropy2(p1,p2,s,ds1,ds2)
      implicit none
C>
C> \brief Calculate the 2-electron entropy associated with an 
C> orbital pair and its derivatives
C>
C> The entropy of a natural orbital pair with occupation numbers
C> \f$p_1\f$ and \f$p_2\f$ is
C> given by
C> \f{eqnarray}{
C>    \f$S^2 &=& p_1 p_2\ln(p_1 p_2) 
C>               +(1-p_1)p_2\ln((1-p_1)p_2) \\\\
C>           &&  +p_1(1-p_2)\ln(p_1(1-p_2))
C>               +(1-p_1)(1-p_2)\ln((1-p_1)(1-p_2))
C> \f}.
C> 
      double precision, intent(in)  :: p1  ! The occupation number of
                                           ! orbital 1.
      double precision, intent(in)  :: p2  ! The occupation number of
                                           ! orbital 2.
      double precision, intent(out) :: s   ! The entropy
      double precision, intent(out) :: ds1 ! The entropy derivative for
                                           ! orbital 1.
      double precision, intent(out) :: ds2 ! The entropy derivative for
                                           ! orbital 2.
c
      s =   p1*p2*log(p1*p2)
     &  +   (1.0d0-p1)*p2*log((1.0d0-p1)*p2)
     &  +   p1*(1.0d0-p2)*log(p1*(1.0d0-p2))
     &  +   (1.0d0-p1)*(1.0d0-p2)*log((1.0d0-p1)*(1.0d0-p2))
      ds1 = p2*log(p1*p2) + (2*p2-1.0d0)
     &    - (1.0d0-p2)*log((1.0d0-p1)*(1.0d0-p2))
      ds2 = p1*log(p1*p2) + (2*p1-1.0d0)
     &    - (1.0d0-p1)*log((1.0d0-p1)*(1.0d0-p2))
c
      end
c
c-----------------------------------------------------------------------
c
c Note: x*ln(x) is the Sterling approximation to ln(Gamma(x+1)) which 
c       the generalization of ln(x!) to real numbers. The latter is
c       referred to as the log-gamma function and its derivative is
c       called the digamma function. The Sterling approximation is very
c       accurate for large x, but we are interested in 0 <= x <= 1.
c       In that domain the Sterling approximation is not very accurate.
c       Worse the digamma function has a finite value at x=0 whereas the
c       derivative of the Stirling approximation has an asymptote at
c       x=0. The latter makes optimizations very hard. Nevertheless
c       we cannot carelessly switch to the log-gamma and digamma 
c       functions as there are few implementations of these functions.
c       The most commonly available one resides in the GNU scientific
c       library (GSL) which is released under the GPL. Hence we cannot
c       use that implementation in NWChem. So we are between a rock and
c       a hard place here.
c
c-----------------------------------------------------------------------
