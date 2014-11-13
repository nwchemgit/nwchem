C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the Heaviside function
C>
C> The Heaviside function [1] is defined as
C> \f{eqnarray*}{
C>   H(x) &=& \int_{-\infty}^{x} \delta(s)\mathrm{d}s
C> \f}
C> Essentially this boils down to a step function that can be defined
C> as
C> \f{eqnarray*}{
C>   H(x) &=& \left\{\begin{array}{cc}
C>              0, & \iff x<0 \\\\ 
C>              \frac{1}{2}, & \iff x=0 \\\\
C>              1, & \iff x>0
C>            \end{array}\right.
C> \f}
C>
C> ### References ###
C>
C> [1] Oliver Heaviside (1850-1925), invented this function to 
C>     model the current in an electric circuit, exact reference
C>     unknown.
C>
      double precision function nwxcm_heaviside(x)
      implicit none
      double precision x !< [Input] The coordinate
c
      if (x<0) then
        nwxcm_heaviside = 0.0d0
      else if (x>0) then
        nwxcm_heaviside = 1.0d0
      else
        nwxcm_heaviside = 0.5d0
      endif
c
      return
      end
C> @}
