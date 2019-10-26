*    ************************************
*    *                                  *
*    *      gen_PBE96_x_unrestricted    *
*    *                                  *
*    ************************************
*
*    This function returns the PBE96 exchange
*  energy density, xe, and its derivatives with respect
*  to nup or ndn and |grad nup| or |grad ndn|.
*
*   Entry - dn_in  : spin densites nup or ndn
*           agr_in : |grad nup| or |grad ndn|
*
*   Exit - xe  : PBE96 energy density
*        - fn  : d(n*xe)/dnup or d(n*xe)/dndn
*        - fdn : d(n*xe)/d|grad nup| or d(n*xe)/d|grad ndn|

      subroutine gen_PBE96_x_unrestricted(dn_in,agr_in,
     >                           xe,fn,fdn)
      implicit none
      
      real*8 dn_in
      real*8 agr_in
      real*8 xe
      real*8 fn
      real*8 fdn
      
*     **** Density cutoff parameter ****
      real*8 DNS_CUT
      parameter (DNS_CUT = 1.0d-20)

c     ***** PBE96 GGA exchange constants ******
      real*8 MU,KAPPA
      parameter (MU    = 0.2195149727645171d0)
      parameter (KAPPA = 0.8040000000000000d0)
 
c     **** other constants ****
      real*8 onethird,fourthird
      parameter (onethird=1.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)

c     **** local variables ****
      integer i
      real*8 n,agr
      real*8 kf,ks,s,P0,n_onethird,pi
      real*8 F
      real*8 Fs                ! dF/ds
      real*8 ex,ex_lda
      real*8 fnx,fdnx          ! d(n*ex)/dnup, d(n*ex)/dndn
      real*8 fdnx_const

      pi = 4.0d0*datan(1.0d0)
      fdnx_const = -3.0d0/(8.0d0*pi)

      
      n     = dn_in
      agr   = agr_in
 
c     ****************************************************************
c     ***** calculate polarized Exchange energies and potentials *****
c     ****************************************************************

      n     = 2.0d0*n
      agr   = 2.0d0*agr

      n_onethird = (3.0d0*n/pi)**onethird
      ex_lda     = -0.75d0*n_onethird

      kf = (3.0d0*pi*pi*n)**onethird
      s  = agr/(2.0d0*kf*n)
      P0 = 1.0d0 + (MU/KAPPA)*s*s

      F   = (1.0d0 + KAPPA - KAPPA/P0)
      Fs  = 2.0d0*MU/(P0*P0)*s

      ex = ex_lda*F
      fnx = fourthird*(ex - ex_lda*Fs*s)
      fdnx = fdnx_const*Fs

      xe  = ex
      fn  = fnx  
      fdn = fdnx 

      
      return
      end
      
*    ************************************
*    *					*
*    *	    gen_PBE96_x_restricted	*
*    *					*
*    ************************************
*
*   This routine calculates the PBE96 exchange 
*   potential(xp) and energy density(xe).
*
*
*   Entry - rho_in : density (nup+ndn)
*           agr_in : |grad rho_in|
*
*     Exit  - xe  : PBE96 exchange energy density
*             fn  : d(n*xe)/dn
*             fdn : d(n*xe/d|grad n|
*
      subroutine gen_PBE96_x_restricted(rho_in,agr_in,
     >                                xe,fn,fdn)
      implicit none

      real*8     rho_in
      real*8     agr_in
      real*8     xe
      real*8     fn
      real*8     fdn

      
*     **** Density cutoff parameter ****
      real*8 DNS_CUT
      parameter (DNS_CUT = 1.0d-20)

c     ***** PBE96 GGA exchange constants ******
      real*8 MU,KAPPA
      parameter (MU    = 0.2195149727645171d0)
      parameter (KAPPA = 0.8040000000000000d0)
 
c     **** other constants ****
      real*8 onethird,fourthird
      parameter (onethird=1.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)

c     **** local variables ****
      real*8 n,agr
      real*8 kf,ks,s,P0,n_onethird,pi
      real*8 fdnx_const
      real*8 F,Fs
      real*8 ex_lda
      real*8 ex
      real*8 fnx,fdnx


      pi         = 4.0d0*datan(1.0d0)
      fdnx_const = -3.0d0/(8.0d0*pi)
      
      n     = rho_in
      agr   = agr_in
      
c     ***** calculate unpolarized Exchange energies and potentials *****
      n_onethird = (3.0d0*n/pi)**onethird
      ex_lda     = -0.75d0*n_onethird

      kf = (3.0d0*pi*pi*n)**onethird
      s  = agr/(2.0d0*kf*n)
      P0 = 1.0d0 + (MU/KAPPA)*s*s

c     if (n.gt.DNS_CUT) then
c        F   = (1.0d0 + KAPPA - KAPPA/P0)
c        Fs  = 2.0d0*MU/(P0*P0)*s
c     else
c        F   = 1.0d0
c        Fs  = 0.0d0
c     end if
      F   = (1.0d0 + KAPPA - KAPPA/P0)
      Fs  = 2.0d0*MU/(P0*P0)*s

      ex   = ex_lda*F
      fnx  = fourthird*(ex - ex_lda*Fs*s) 
      fdnx = fdnx_const*Fs


      xe  = ex   
      fn  = fnx  
      fdn = fdnx 
         

c       write(*,*) "pbe96:",i,ec,fnc,fdnc


      return
      end
 
