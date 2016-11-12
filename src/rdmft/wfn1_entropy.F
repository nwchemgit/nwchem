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
      subroutine wfn1_entropy2_ab(pa,pb,s,dsa,dsb)
      implicit none
C>
C> \brief Calculate the 2-electron entropy associated with an 
C> orbital pair and its derivatives
C>
C> The entropy of a natural orbital pair with occupation numbers
C> \f$p_\alpha\f$ and \f$p_\beta\f$ is
C> given by \f$S^2_{ab} = p_\alpha p_\beta\ln(p_\alpha p_\beta)+(1-p_\alpha)(1-p_\beta)\ln((1-p_\alpha)(1-p_\beta))\f$.
C> 
      double precision, intent(in)  :: pa  ! The \f$\alpha\f$ occupation
                                           ! number
      double precision, intent(in)  :: pb  ! The \f$\beta\f$ occupation
                                           ! number
      double precision, intent(out) :: s   ! The entropy
      double precision, intent(out) :: dsa ! The entropy derivative for
                                           ! the \f$\alpha\f$-spin
      double precision, intent(out) :: dsb ! The entropy derivative for
                                           ! the \f$\beta\f$-spin
c
      s =   pa*pb*log(pa*pb)
     &  +   (1.0d0-pa)*(1.0d0-pb)*log((1.0d0-pa)*(1.0d0-pb))
      dsa = pb*log(pa*pb) + (2*pb-1.0d0)
     &    - (1.0d0-pb)*log((1.0d0-pa)*(1.0d0-pb))
      dsb = pa*log(pa*pb) + (2*pa-1.0d0)
     &    - (1.0d0-pa)*log((1.0d0-pa)*(1.0d0-pb))
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_entropy2_aa(pa,s,dsa)
      implicit none
C>
C> \brief Calculate the 2-electron entropy associated with an 
C> equal spin orbital pair and its derivatives
C>
C> The entropy of an equal spin natural orbital pair with occupation
C> numbers \f$p_{\sigma\sigma}\f$ is given by:
C> \f{eqnarray}{
C>   p &=& \sqrt{p_{\sigma\sigma}} \\\\
C>   S^2_{\sigma\sigma} &=& p^2\ln(p^2)+(1-p)^2\ln((1-p)^2)
C> \f}.
C> 
      double precision, intent(in)  :: pa  ! The \f$\alpha\f$ occupation
                                           ! number
      double precision, intent(out) :: s   ! The entropy
      double precision, intent(out) :: dsa ! The entropy derivative for
                                           ! the \f$\alpha\f$-spin
c
      double precision :: pa12 ! The square root of the occupation 
                               ! number
      double precision :: ds1  ! A derivative component
      double precision :: ds2  ! A derivative component
c
      pa12 = sqrt(pa)
      call wfn1_entropy2_ab(pa12,pa12,s,ds1,ds2)
      dsa = 0.5d0*(ds1+ds2)/pa12
c
      end
c
c-----------------------------------------------------------------------
