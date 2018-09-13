C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_query.F
C> Routines to query the properties of the current functional
C> 
C> @}
C
C> \ingroup nwpwxc_api
C> @{
C>
C> \brief Is this library used to evaluate the functional
C>
C> Certainly in interim stages one might want to choose to use
C> this library or another mechanism to evaluate the functional.
C> This function returns .true. if this library has been activated.
C>
      logical function nwpwxc_is_on()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_is_on = nwpwxc_on
      return
      end
C>
C> \brief Does this functional have Hartree-Fock exchange?
C>
C> A function to query the functional for the presence of a
C> Hartree-Fock exchange term. This function just queries the 
C> functional definition. In practice the Hartree-Fock exchange term
C> must be implemented in the application code as this library cannot
C> provide that.
C>
C> \return Returns .true. if the functional contains a Hartree-Fock 
C> exchange term (i.e. is a hybrid functional), and .false. otherwise.
C>
      logical function nwpwxc_has_hfx()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_has_hfx = nwpwxc_wght_hf .gt. 0.0d0
      return
      end
C>
C> \brief Does this functional have MP2?
C>
C> A function to query the functional for the presence of an MP2
C> correlation term. This function just queries the 
C> functional definition. In practice the MP2 correlation term
C> must be implemented in the application code as this library cannot
C> provide that.
C>
C> \return Returns .true. if the functional contains an MP2 correlation
C> term (i.e. is a double hybrid functional), and .false. otherwise.
C>
      logical function nwpwxc_has_mp2c()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_has_mp2c = nwpwxc_wght_mp2 .gt. 0.0d0
      return
      end
C>
C> \brief Does this functional have Coulomb Attenuation?
C>
C> A function to query the functional for the use Coulomb Attentuated
C> operators. With Coulomb Attentuation part of the expression is 
C> incorporate in the Density Functional component and part needs to
C> be done by the application code.
C>
C> \return Returns .true. if the functional use Coulomb Attenuation,
C> and .false. otherwise.
C>
      logical function nwpwxc_has_cam()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_has_cam = nwpwxc_cam_gamma .ne. 0.0d0 .or.
     +               nwpwxc_cam_beta  .gt. 0.0d0 .or.
     +               nwpwxc_cam_alpha .gt. 0.0d0
      return
      end
C>
C> \brief Does this functional have dispersion corrections?
C>
C> A function to query the functional for the presence of a
C> dispersion correction term. This function just queries the 
C> functional definition. In practice the dispersion correction term
C> must be implemented in the application code as this library cannot
C> provide that.
C>
C> \return Returns .true. if the functional contains a dispersion
C> correction term, and .false. otherwise.
C>
      logical function nwpwxc_has_disp()
      implicit none
#include "nwpwxcP.fh"
#include "nwpwxc_param.fh"
      nwpwxc_has_disp = nwpwxc_disp_id.ne.NWXC_DISP_NONE
      return
      end
C>
C> \brief Retrieve the density tolerance
C>
C> Retrieve the current setting of the density tolerance.
C> In the evaluation of the functional densities less 
C> than this tolerance will be ignored.
C>
      double precision function nwpwxc_get_rho_tolerance()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_get_rho_tolerance = nwpwxc_rho_tol
      return
      end
C>
C> \brief Retrieve the Coulomb Attenuation parameters
C>
C> Using the Coulomb Attenuation or range separation the Coulomb
C> operator is written as:
C> \f{eqnarray*}{
C>   \frac{1}{r_{12}} &=& \frac{1-[\alpha+\beta\mathrm{erf}(\gamma r_{12})]}{r_{12}} 
C>                     +  \frac{\alpha+\beta\mathrm{erf}(\gamma r_{12})}{r_{12}}
C> \f}
C> This routine retrieves the parameters \f$ \alpha \f$, \f$ \beta \f$,
C> and \f$ \gamma \f$. In addition the argument \f$ srhfx \f$ 
C> indicates whether the short-range term should be evaluated using
C> the Hartree-Fock formulation of exchange (`srhfx = .true.`) or the
C> long range part (`srhfx = .false.`).
C>
C> To enable returning multiple values at once this routine is
C> implemented as a subroutine.
C>
      subroutine nwpwxc_get_cam(alpha,beta,gamma,srhfx)
      implicit none
#include "nwpwxcP.fh"
      double precision alpha !< [Output] The alpha parameter
      double precision beta  !< [Output] The beta parameter
      double precision gamma !< [Output] The gamma parameter
      logical srhfx !< [Output] T if short-range exchange is 
                    !< Hartree-Fock exchange
      alpha = nwpwxc_cam_alpha
      beta  = nwpwxc_cam_beta
      gamma = nwpwxc_cam_gamma
      srhfx = nwpwxc_cam_srhfex
      end
C>
C> \brief Retrieve the dispersion correction parameters
C>
C> A variety of dispersion corrections have been suggested in the 
C> literature. These corrections consist of \f$\frac{S_6}{R^6}\f$
C> and \f$\frac{S_8}{R^8}\f$ factors combined with a damping factor
C> \f$f(\alpha;R)\f$. The exact form of the dispersion correction
C> comes in a few different flavors. So to specify which expression
C> should be used this subroutine returns the \f$S\f$-factors, the
C> \f$\alpha\f$-damping parameter as well as the DispID. It
C> is assumed that the host code implements the actual expression
C> and has tables of the atomic \f$C_6\f$, and if applicable, \f$C_8\f$
C> coeffcients.
C>
      subroutine nwpwxc_get_disp(dispid,s6,s8,sr6,sr8,alpha)
      implicit none
#include "nwpwxcP.fh"
      integer dispid      !< [Output] The dispersion identifier
      double precision s6 !< [Output] The \f$S_6\f$ scale factor
      double precision s8 !< [Output] The \f$S_8\f$ scale factor
      double precision sr6 !< [Output] The \f$S_{r6}\f$ scale factor
      double precision sr8 !< [Output] The \f$S_{r8}\f$ scale factor
      double precision alpha !< [Output] The \f$\alpha\f$ damping factor
      dispid = nwpwxc_disp_id
      s6     = nwpwxc_disp_s6
      s8     = nwpwxc_disp_s8
      sr6    = nwpwxc_disp_sr6
      sr8    = nwpwxc_disp_sr8
      alpha  = nwpwxc_disp_alpha
      end
C>
C> \brief Retrieve the Hartree-Fock exchange weight
C>
C> Simply look up and return the weight for the Hartree-Fock exchange.
C>
C> \return Returns the weight of the Hartree-Fock exchange term
C>
      double precision function nwpwxc_wght_hfx()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_wght_hfx = nwpwxc_wght_hf
      return
      end
C>
C> \brief Retrieve the MP2 correlation weight
C>
C> Simply look up and return the weight for the MP2 correlation term.
C>
C> \return Returns the weight of the MP2 correlation term
C>
      double precision function nwpwxc_wght_mp2c()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_wght_mp2c = nwpwxc_wght_mp2
      return
      end
C>
C> \brief Does this functional have any exchange-correlation terms?
C>
C> A function to query the functional for the presence of any 
C> exchange-correlation terms. With the current hybrid and double hybrid
C> functionals it has become possible to specify functionals that 
C> produce just pure Hartree-Fock or MP2 without any DFT 
C> exchange-correlation at all. Although unusual the code has to be able
C> to handle these cases correctly. Hence there is a need to be able to
C> query whether any exchange-correlation terms need to be evaluated.
C>
C> \return Return .true. if the functional has an exchange-correlation
C> term of any kind, and .false. otherwise.
C>
      logical function nwpwxc_has_xc()
      implicit none
      logical  nwpwxc_is_lda, nwpwxc_is_gga, nwpwxc_is_mgga
      external nwpwxc_is_lda, nwpwxc_is_gga, nwpwxc_is_mgga
      nwpwxc_has_xc = nwpwxc_is_lda() .or. nwpwxc_is_gga() .or. 
     >                nwpwxc_is_mgga()
      return
      end
C>
C> \brief Is this functional an LDA functional?
C>
C> An LDA functional is a functional that depends on the local density
C> only. 
C>
C> The recommended way to query the nature of a functional (e.g. for
C> memory allocation purposes) is:
C> \code
C>    if (nwpwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwpwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwpwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density only,
C> and .false. otherwise.
C>
      logical function nwpwxc_is_lda()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_is_lda = nwpwxc_depends_on_rho .and. 
     +              (.not. nwpwxc_depends_on_gamma) .and.
     +              (.not. nwpwxc_depends_on_tau)
      return
      end
C>
C> \brief Is this functional a GGA functional?
C>
C> A GGA functional is a functional that depends on the local density
C> and the density gradient, but not on the kinetic energy density.
C>
C> The recommended way to query the nature of a functional (e.g. for
C> memory allocation purposes) is:
C> \code
C>    if (nwpwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwpwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwpwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density and
C> the density gradient, and .false. otherwise.
C>
      logical function nwpwxc_is_gga()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_is_gga = nwpwxc_depends_on_rho .and. 
     +              nwpwxc_depends_on_gamma .and.
     +              (.not. nwpwxc_depends_on_tau)
      return
      end
C>
C> \brief Is this functional a meta-GGA functional?
C>
C> A meta-GGA functional is a functional that depends on the local
C> density, the density gradient, and on the kinetic energy density.
C>
C> The recommended way to query the nature of a functional (e.g. for
C> memory allocation purposes) is:
C> \code
C>    if (nwpwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwpwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwpwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density,
C> the density gradient, and the kinetic energy density, and .false.
C> otherwise.
C>
      logical function nwpwxc_is_mgga()
      implicit none
#include "nwpwxcP.fh"
      nwpwxc_is_mgga = nwpwxc_depends_on_rho .and. 
     +               nwpwxc_depends_on_gamma .and.
     +               nwpwxc_depends_on_tau
      return
      end
C> @}
C
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Extract information of specific density functional terms from the
C> user list
C>
C> This function queries the current density functional. It extracts the
C> functional identifier and the corresponding weight for a given entry in
C> the user list. The list entry is incremented upon return. This allows
C> one to extract the whole functional with the following loop:
C> \code
C>    ientry = 1
C>    do while (nwpwxc_get_df(ientry,func_id,weight))
C>      write(*,*)"functional id: ",func_id," weight=",weight
C>    enddo
C> \endcode
C>
C> \return Return .true. if the entry number is greater than or equal to 1 and
C> and less than or equal to the number of density functional terms, and .false.
C> otherwise.
C>
       logical function nwpwxc_get_df(ientry,func_id,weight)
       implicit none
#include "nwpwxcP.fh"
       integer ientry  !< [In/Output] On input the entry number to extract,
                       !< on output the input value incremented by one
       integer func_id !< [Output] The functional identifier
       double precision weight !< [Output] The functional weight
c
       nwpwxc_get_df = .true.
       if (.not.((ientry.ge.1).and.(ientry.le.nwpwxc_usr_num))) then
         nwpwxc_get_df = .false.
         return
       endif
       func_id = nwpwxc_usr_ids(ientry)
       weight  = nwpwxc_usr_wghts(ientry)
       ientry  = ientry + 1
c
       return
       end
C>
C> \brief Extract information about a particular functional id
C>
C> Retrieve information about a functional id. This information consists
C> of the name and literature reference of that functional. This information
C> is returned in a string so as to be convenient for printing, for example.
C>
C> \return Return .true. if the data was found and could be copied successfully
C> to the buffer provided, and return .false. otherwise.
C>
       logical function nwpwxc_get_info(func_id,cdata)
       implicit none
#include "nwpwxcP.fh"
       integer func_id     !< [Input] The functional id
       character*(*) cdata !< [Output] Information on the functional
c
       integer lendat      !< the length of the cdata buffer
       logical stat        !< the function status
       integer i           !< counter
c
       stat   = .true.
       lendat = len(cdata)
       select case (func_id)
         case (NWXCP_HFXC)
           cdata = "Hartree-Fock Exchange"
         case (NWXCP_MP2C)
           cdata = "Moller-Plesset 2nd Order Correlation"
         case (NWXCP_X_SLATER)
           cdata = "Slater Exchange (JC Slater, Phys.Rev. 81, 385 "
     +             //"(1951) doi:10.1103/PhysRev.81.385)"
         case (NWXCP_X_B0)
           cdata = "Becke0 Exchange (ACM0 variant of Becke88)"
         case (NWXCP_X_B86B)
           cdata = "Becke86B Exchange (AD Becke, "
     +             //"J.Chem.Phys. 85, 7184 (1986)"
     +             //" doi:10.1063/1.451353)"
         case (NWXCP_X_B88)
           cdata = "Becke88 Exchange (AD Becke, "
     +             //"Phys.Rev.A 38, 3098 (1988)"
     +             //" doi:10.1103/PhysRevA.38.3098)"
         case (NWXCP_X_B97)
           cdata = "B97 Exchange (AD Becke, "
     +             //"J.Chem.Phys. 107, 8554 (1997)"
     +             //" doi:10.1063/1.475007)"
         case (NWXCP_X_B97_1)
           cdata = "B97-1 Exchange (FA Hamprecht, A Cohen, DJ Tozer, "
     +             //"NC Handy, J.Chem.Phys. 109, 6264 (1998)"
     +             //" doi:10.1063/1.477267)"
         case (NWXCP_X_B97_2)
           cdata = "B97-2 Exchange (PJ Wilson, TJ Bradley, DJ Tozer, "
     +             //"J.Chem.Phys. 115, 9233 (2001)"
     +             //" doi:10.1063/1.1412605)"
         case (NWXCP_X_B97_3)
           cdata = "B97-3 Exchange (TW Keal, DJ Tozer, "
     +             //"J.Chem.Phys. 123, 121103 (2005)"
     +             //" doi:10.1063/1.2061227)"
         case (NWXCP_X_B97_D)
           cdata = "B97-D Exchange (S Grimme, "
     +             //"J.Comp.Chem. 27, 1787 (2006)"
     +             //" doi:10.1002/jcc.20495)"
         case (NWXCP_X_B97_G)
           cdata = "B97-GGA1 Exchange (AJ Cohen, NC Handy, "
     +             //"Chem.Phys.Lett. 316, 160 (2000)"
     +             //" doi:10.1016/S0009-2614(99)01273-7)"
         case (NWXCP_X_B98)
           cdata = "B98 Exchange (HL Schmider, AD Becke, "
     +             //"J.Chem.Phys. 108, 9624 (1998)"
     +             //" doi:10.1063/1.476438)"
         case (NWXCP_X_BNL)
           cdata = "BNL Exchange (R Baer, D Neuhauser, E Livshits, "
     +             //"Phys.Chem.Chem.Phys. 9, 2932 (2007)"
     +             //" doi:10.1039/B617919C)"
         case (NWXCP_X_CAMB88)
           cdata = "CAM-B88 Exchange (T Yanai, DP Tew, "
     +             //"NC Handy, Chem.Phys.Lett. 393, 51 (2004)"
     +             //" doi:10.1016/j.cplett.2004.06.011)"
         case (NWXCP_X_CAMLSD)
           cdata = "CAM-LSD Exchange (T Yanai, DP Tew, "
     +             //"NC Handy, Chem.Phys.Lett. 393, 51 (2004)"
     +             //" doi:10.1016/j.cplett.2004.06.011)"
         case (NWXCP_X_CAMPBE)
           cdata = "CAM-PBE Exchange (OA Vydrov, GE Scuseria, "
     +             //"JP Perdew, J.Chem.Phys. 126, 154109 (2007)"
     +             //" doi:10.1063/1.2723119)"
         case (NWXCP_X_CAMS12G)
           cdata = "CAM-S12g Exchange (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_X_CAMS12H)
           cdata = "CAM-S12h Exchange (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_X_DLDF)
           cdata = "DLDF Exchange (K Pernal, R Podeszwa, "
     +             //"K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, "
     +             //"263201 (2009) doi:10.1103/PhysRevLett.103.263201)"
         case (NWXCP_X_FT97)
           cdata = "FT97 B Exchange (M Filatov, W Thiel, "
     +             //"Mol.Phys. 91, "
     +             //"847 (1997) doi:10.1080/002689797170950)"
         case (NWXCP_X_GILL)
           cdata = "GILL Exchange (PMW Gill, "
     +             //"Mol.Phys. 89, 433 (1996) "
     +             //"doi:10.1080/002689796173813)"
         case (NWXCP_X_HCTH)
           cdata = "HCTH Exchange (FA Hamprecht, A Cohen, "
     +             //"DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)"
     +             //" doi:10.1063/1.477267)"
         case (NWXCP_X_HCTH120)
           cdata = "HCTH120 Exchange (AD Boese, NL Doltsinis, "
     +             //"NC Handy, M Sprik, J.Chem.Phys. 112, 1670 (2000)"
     +             //" doi:10.1063/1.480732)"
         case (NWXCP_X_HCTH147)
           cdata = "HCTH147 Exchange (AD Boese, NL Doltsinis, "
     +             //"NC Handy, M Sprik, J.Chem.Phys. 112, 1670 (2000)"
     +             //" doi:10.1063/1.480732)"
         case (NWXCP_X_HCTH407)
           cdata = "HCTH407 Exchange (AD Boese, "
     +             //"NC Handy, J.Chem.Phys. 114, 5497 (2001)"
     +             //" doi:10.1063/1.1347371)"
         case (NWXCP_X_HCTH407P)
           cdata = "HCTH407P Exchange (AD Boese, A Chandra, JML Martin,"
     +             //" D Marx, J.Chem.Phys. 119, 5965 (2003)"
     +             //" doi:10.1063/1.1599338)"
         case (NWXCP_X_HCTH_A)
           cdata = "HCTH(AC) Exchange (DJ Tozer, NC Handy, "
     +             //"J.Chem.Phys. 109, 10180 (1998)"
     +             //" doi:10.1063/1.477711)"
         case (NWXCP_X_HCTHP14)
           cdata = "HCTHP14 Exchange (G Menconi, PJ Wilson, "
     +             //"DJ Tozer, J.Chem.Phys. 114, 3958 (2001)"
     +             //" doi:10.1063/1.1342776)"
         case (NWXCP_X_M05)
           cdata = "M05 Exchange (Y Zhao, NE Schultz, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
         case (NWXCP_X_M05_2X)
           cdata = "M05-2X Exchange (Y Zhao, NE Schultz, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
         case (NWXCP_X_M06)
           cdata = "M06 Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_X_M06_HF)
           cdata = "M06-HF Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_X_M06_L)
           cdata = "M06-L Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_X_M06_2X)
           cdata = "M06-2X Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_X_M08_HX)
           cdata = "M08-HX Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Chem.TheoryComput. 4, 1849 (2008)"
     +             //" doi:10.1021/ct800246v)"
         case (NWXCP_X_M08_SO)
           cdata = "M08-SO Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Chem.TheoryComput. 4, 1849 (2008)"
     +             //" doi:10.1021/ct800246v)"
         case (NWXCP_X_M11)
           cdata = "M11 Exchange (R Peverati, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 2810 (2011)"
     +             //" doi:10.1021/jz201170d)"
         case (NWXCP_X_M11_L)
           cdata = "M11-L Exchange (R Peverati, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 2810 (2011)"
     +             //" doi:10.1021/jz201170d)"
         case (NWXCP_X_MPW91)
           cdata = "MPW91 Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem.A 109, 5656 (2005)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_X_OPT)
           cdata = "OPT Exchange (NC Handy, AJ Cohen, "
     +             //"Mol.Phys. 99, 403 (2001)"
     +             //" doi:10.1080/00268970010018431)"
         case (NWXCP_X_PW6B95)
           cdata = "PW6B95 Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem. A 109, 5656 (2005)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_X_PWB6K)
           cdata = "PWB6K Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem. A 109, 5656 (2005)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_X_PW91)
           cdata = "PW91 Exchange (JP Perdew, JA Chevary, "
     +             //"SH Vosko, KA Jackson, MR Pederson, DJ Singh, "
     +             //"C. Fiolhais, Phys.Rev. B 46, 6671 (1992)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_X_PBE)
           cdata = "PBE Exchange (JP Perdew, K Burke, M Ernzerhof, "
     +             //"Phys.Rev.Lett. 77, 3865 (1996)"
     +             //" doi:10.1103/PhysRevLett.77.3865)"
         case (NWXCP_X_PBE0)
           cdata = "PBE0 Exchange (C Adamo, V Barone, "
     +             //"J.Chem.Phys. 110, 6158 (1999)"
     +             //" doi:10.1063/1.478522)"
         case (NWXCP_X_PBESOL)
           cdata = "PBESOL Exchange (JP Perdew, A Ruzsinszky, "
     +             //"GI Csonka, et al., "
     +             //"Phys.Rev.Lett. 100, 136406 (2008)"
     +             //" doi:10.1103/PhysRevLett.100.136406)"
         case (NWXCP_X_REVPBE)
           cdata = "revPBE Exchange (Y Zhang, W Yang, "
     +             //"Phys.Rev.Lett. 80, 890 (1998)"
     +             //" doi:10.1103/PhysRevLett.80.890)"
         case (NWXCP_X_REVPBE0)
           cdata = "revPBE0 Exchange (ACM0 variant of REVPBE Exchange)"
         case (NWXCP_X_RPBE)
           cdata = "RPBE Exchange (B Hammer, LB Hanssen, JK Norskov, "
     +             //"Phys.Rev.B 59, 7413 (1999)"
     +             //" doi:10.1103/PhysRevB.59.7413)"
         case (NWXCP_X_PKZB)
           cdata = "PKZB Exchange (JP Perdew, S Kurth, A Zupan, "
     +             //"P Blaha, Phys.Rev.Lett. 82, "
     +             //"2544 (1999) doi:10.1103/PhysRevLett.82.2544)"
         case (NWXCP_X_S12G)
           cdata = "S12g Exchange (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_X_S12H)
           cdata = "S12h Exchange (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_X_SOGGA)
           cdata = "SOGGA Exchange (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 128, 184109 (2008)"
     +             //" doi:10.1063/1.2912068)"
         case (NWXCP_X_SOGGA11)
           cdata = "SOGGA11 Exchange (R Peverati, Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 1991 (2011)"
     +             //" doi:10.1021/jz200616w)"
         case (NWXCP_X_SOGGA11_X)
           cdata = "SOGGA11-X Exchange (R Peverati, DG Truhlar, "
     +             //"J.Chem.Phys. 135, 191102 (2011)"
     +             //" doi:10.1063/1.3663871)"
         case (NWXCP_X_SSB_D)
           cdata = "SSB-D Exchange (M Swart, M Sola, FM Bickelhaupt, "
     +             //"J.Chem.Phys. 131, 094103 (2009)"
     +             //" doi:10.1063/1.3213193)"
         case (NWXCP_X_TPSS)
           cdata = "TPSS Exchange (J Tao, JP Perdew, VN Staveroverov, "
     +             //"GE Scuseria, Phys.Rev.Lett. 91, 146401 (2003) "
     +             //"doi:10.1103/PhysRevLett.91.146401)"
         case (NWXCP_X_TPSSH)
           cdata = "TPSSH Exchange (VN Staveroverov, GE Scuseria, "
     +             //"J Tao, JP Perdew, J.Chem.Phys. 119, "
     +             //"12129 (2003) doi:10.1063/1.1626543)"
         case (NWXCP_X_WPBE)
           cdata = "WPBE Exchange (E Weintraub, TM Henderson, "
     +             //"GE Scuseria, J.Chem.Theory.Comput. 5, 754 (2009)"
     +             //" doi:10.1021/ct800530u)"
         case (NWXCP_X_VS98)
           cdata = "VS98 Exchange (T van Voorhis, GE Scuseria, "
     +             //"J.Chem.Phys. 109, 400 (1998)"
     +             //" doi:10.1063/1.3005348)"
         case (NWXCP_X_SCAN)
           cdata = "SCAN Exchange (T van Voorhis, GE Scuseria, "
     +             //"J.Chem.Phys. 109, 400 (1998)"
     +             //" doi:10.1063/1.3005348)"
         case (NWXCP_C_B95)
           cdata = "B95 Correlation (AD Becke, "
     +             //"J.Chem.Phys. 104, 1040 (1996)"
     +             //" doi:10.1063/1.4708298)"
         case (NWXCP_C_B97)
           cdata = "B97 Correlation (AD Becke, "
     +             //"J.Chem.Phys. 107, 8554 (1997)"
     +             //" doi:10.1063/1.475007)"
         case (NWXCP_C_B97_1)
           cdata = "B97-1 Correlation (FA Hamprecht, A Cohen, "
     +             //"DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)"
     +             //" doi:10.1063/1.477267)"
         case (NWXCP_C_B97_2)
           cdata = "B97-2 Correlation (PJ Wilson, TJ Bradley, "
     +             //"DJ Tozer, J.Chem.Phys. 115, 9233 (2001)"
     +             //" doi:10.1063/1.1412605)"
         case (NWXCP_C_B97_3)
           cdata = "B97-3 Correlation (TW Keal, DJ Tozer, "
     +             //"J.Chem.Phys. 123, 121103 (2005)"
     +             //" doi:10.1063/1.2061227)"
         case (NWXCP_C_B97_D)
           cdata = "B97-D Correlation (S Grimme, "
     +             //"J.Comp.Chem. 27, 1787 (2006)"
     +             //" doi:10.1002/jcc.20495)"
         case (NWXCP_C_B97_G)
           cdata = "B97-GGA1 Correlation (AJ Cohen, NC Handy, "
     +             //"Chem.Phys.Lett. 316, 160 (2000)"
     +             //" doi:10.1016/S0009-2614(99)01273-7)"
         case (NWXCP_C_B98)
           cdata = "B98 Correlation (HL Schmider, AD Becke, "
     +             //"J.Chem.Phys. 108, 9624 (1998)"
     +             //" doi:10.1063/1.476438)"
         case (NWXCP_C_CAMS12G)
           cdata = "CAM-S12g Correlation (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_C_CAMS12H)
           cdata = "CAM-S12h Correlation (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_C_DLDF)
           cdata = "DLDF Correlation (K Pernal, R Podeszwa, "
     +             //"K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, "
     +             //"263201 (2009) doi:10.1103/PhysRevLett.103.263201)"
         case (NWXCP_C_FT97)
           cdata = "FT97 Correlation (M Filatov, W Thiel, "
     +             //"Int.J.Quant.Chem. 62, 603 (1997) doi:10.1002/(SIC"
     +             //"I)1097-461X(1997)62:6<603::AID-QUA4>3.0.CO;2-#)"
         case (NWXCP_C_HCTH)
           cdata = "HCTH Correlation (FA Hamprecht, A Cohen, "
     +             //"DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)"
     +             //" doi:10.1063/1.477267)"
         case (NWXCP_C_HCTH120)
           cdata = "HCTH120 Correlation (AD Boese, NL Doltsinis, "
     +             //"NC Handy, M Sprik, J.Chem.Phys. 112, 1670 (2000)"
     +             //" doi:10.1063/1.480732)"
         case (NWXCP_C_HCTH147)
           cdata = "HCTH147 Correlation (AD Boese, NL Doltsinis, "
     +             //"NC Handy, M Sprik, J.Chem.Phys. 112, 1670 (2000)"
     +             //" doi:10.1063/1.480732)"
         case (NWXCP_C_HCTH407)
           cdata = "HCTH407 Correlation (AD Boese, "
     +             //"NC Handy, J.Chem.Phys. 114, 5497 (2001)"
     +             //" doi:10.1063/1.1347371)"
         case (NWXCP_C_HCTH407P)
           cdata = "HCTH407P Correlation (AD Boese, A Chandra, "
     +             //"JML Martin, D Marx, J.Chem.Phys. 119, 5965 (2003)"
     +             //" doi:10.1063/1.1599338)"
         case (NWXCP_C_HCTH_A)
           cdata = "HCTH(AC) Correlation (DJ Tozer, NC Handy, "
     +             //"J.Chem.Phys. 109, 10180 (1998)"
     +             //" doi:10.1063/1.477711)"
         case (NWXCP_C_HCTHP14)
           cdata = "HCTHP14 Correlation (G Menconi, PJ Wilson, "
     +             //"DJ Tozer, J.Chem.Phys. 114, 3958 (2001)"
     +             //" doi:10.1063/1.1342776)"
         case (NWXCP_C_OP)
           cdata = "OP Correlation (T Tsuneda, T Suzumura, K Hirao, "
     +             //"J.Chem.Phys. 110, 10664 (1999)"
     +             //" doi:10.1063/1.479012)"
         case (NWXCP_C_OPT)
           cdata = "OPT Correlation (NC Handy, AJ Cohen, "
     +             //"Mol.Phys. 99, 607 (2001)"
     +             //" doi:10.1080/00268970010023435)"
         case (NWXCP_C_PZ81)
           cdata = "Perdew81 Correlation (JP Perdew, A Zunger, "
     +             //"Phys.Rev. B 23, 5048 (1981)"
     +             //" doi:10.1103/PhysRevB.23.5048)"
         case (NWXCP_C_P86)
           cdata = "Perdew86 Correlation (JP Perdew, "
     +             //"Phys.Rev. B 33, 8822 (1986)"
     +             //" doi:10.1103/PhysRevB.33.8822)"
         case (NWXCP_C_P91)
           cdata = "Perdew91 Correlation (JP Perdew, JA Chevary, "
     +             //"SH Vosko, KA Jackson, MR Pederson, DJ Singh, "
     +             //"C. Fiolhais, Phys.Rev. B 46, 6671 (1992)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_C_P91_VWN5)
           cdata = "Perdew91 with VWN5 Correlation"
         case (NWXCP_C_PW6B95)
           cdata = "PW6B95 Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem. A 109, 5656 (2005)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_C_PWB6K)
           cdata = "PWB6K Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem. A 109, 5656 (2005)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_C_PW91LDA)
           cdata = "PW91 LDA Correlation (JP Perdew, Y Wang, "
     +             //"Phys.Rev.B 45, 13244 (1992)"
     +             //" doi:10.1103/PhysRevB.45.13244)"
         case (NWXCP_C_LYP)
           cdata = "LYP Correlation (C Lee, W Yang, RG Parr, "
     +             //"Phys.Rev.B 37, 785 (1988)"
     +             //" doi:10.1103/PhysRevB.37.785)"
         case (NWXCP_C_M05)
           cdata = "M05 Correlation (Y Zhao, NE Schultz, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
         case (NWXCP_C_M05_2X)
           cdata = "M05-2X Correlation (Y Zhao, NE Schultz, DG Truhlar,"
     +             //" J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
         case (NWXCP_C_M06)
           cdata = "M06 Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_C_M06_HF)
           cdata = "M06-HF Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_C_M06_L)
           cdata = "M06-L Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_C_M06_2X)
           cdata = "M06-2X Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_C_M08_HX)
           cdata = "M08-HX Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Chem.TheoryComput. 4, 1849 (2008)"
     +             //" doi:10.1021/ct800246v)"
         case (NWXCP_C_M08_SO)
           cdata = "M08-SO Correlation (Y Zhao, DG Truhlar, "
     +             //"J.Chem.TheoryComput. 4, 1849 (2008)"
     +             //" doi:10.1021/ct800246v)"
         case (NWXCP_C_M11)
           cdata = "M11 Correlation (R Peverati, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 2810 (2011)"
     +             //" doi:10.1021/jz201170d)"
         case (NWXCP_C_M11_L)
           cdata = "M11-L Correlation (R Peverati, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 2810 (2011)"
     +             //" doi:10.1021/jz201170d)"
         case (NWXCP_C_PBE)
           cdata = "PBE Correlation (JP Perdew, K Burke, M Ernzerhof, "
     +             //"Phys.Rev.Lett. 77, 3865 (1996)"
     +             //" doi:10.1103/PhysRevLett.77.3865)"
         case (NWXCP_C_PBESOL)
           cdata = "PBESOL Correlation (JP Perdew, A Ruzsinszky, "
     +             //"GI Csonka, et al., "
     +             //"Phys.Rev.Lett. 100, 136406 (2008)"
     +             //" doi:10.1103/PhysRevLett.100.136406)"
         case (NWXCP_C_PKZB)
           cdata = "PKZB Correlation (JP Perdew, S Kurth, A Zupan, "
     +             //"P Blaha, Phys.Rev.Lett. 82, "
     +             //"2544 (1999) doi:10.1103/PhysRevLett.82.2544)"
         case (NWXCP_C_S12G)
           cdata = "S12g Correlation (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_C_S12H)
           cdata = "S12h Correlation (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_C_SPBE)
           cdata = "SSB-D PBE Correlation (M Swart, M Sola, "
     +             //"FM Bickelhaupt, J.Chem.Phys. 131, 094103 (2009)"
     +             //" doi:10.1063/1.3213193)"
         case (NWXCP_C_SOGGA11)
           cdata = "SOGGA11 Correlation (R Peverati, Y Zhao, "
     +             //"DG Truhlar, J.Phys.Chem.Lett. 2, 1991 (2011)"
     +             //" doi:10.1021/jz200616w)"
         case (NWXCP_C_SOGGA11_X)
           cdata = "SOGGA11-X Correlation (R Peverati, DG Truhlar, "
     +             //"J.Chem.Phys. 135, 191102 (2011)"
     +             //" doi:10.1063/1.3663871)"
         case (NWXCP_C_TPSS)
           cdata = "TPSS Correlation (J Tao, JP Perdew, VN "
     +             //"Staveroverov, GE Scuseria, Phys.Rev.Lett. 91, "
     +             //"146401 (2003) doi:10.1103/PhysRevLett.91.146401)"
         case (NWXCP_C_VS98)
           cdata = "VS98 Correlation (T van Voorhis, GE Scuseria, "
     +             //"J.Chem.Phys. 109, 400 (1998)"
     +             //" doi:10.1063/1.3005348)"
         case (NWXCP_C_SCAN)
           cdata = "SCAN Correlation (T van Voorhis, GE Scuseria, "
     +             //"J.Chem.Phys. 109, 400 (1998)"
     +             //" doi:10.1063/1.3005348)"
         case (NWXCP_C_VWN1)
           cdata = "VWN1 Correlation (SH Vosko, L Wilk, M Nusair, "
     +             //"Can.J.Phys. 58, 1200 (1980) doi:10.1139/p80-159)"
         case (NWXCP_C_VWN1_RPA)
           cdata = "VWN1RPA Correlation (SH Vosko, L Wilk, M Nusair, "
     +             //"Can.J.Phys. 58, 1200 (1980) doi:10.1139/p80-159)"
         case (NWXCP_C_VWN2)
           cdata = "VWN2 Correlation (SH Vosko, L Wilk, M Nusair, "
     +             //"Can.J.Phys. 58, 1200 (1980) doi:10.1139/p80-159)"
         case (NWXCP_C_VWN3)
           cdata = "VWN3 Correlation (SH Vosko, L Wilk, M Nusair, "
     +             //"Can.J.Phys. 58, 1200 (1980) doi:10.1139/p80-159)"
         case (NWXCP_C_VWN4)
           cdata = "VWN4 Correlation (SH Vosko, L Wilk, M Nusair, "
     +             //"Can.J.Phys. 58, 1200 (1980) doi:10.1139/p80-159)"
         case (NWXCP_C_VWN5)
           cdata = "VWN5 Correlation (SH Vosko, L Wilk, M Nusair, "
     +             //"Can.J.Phys. 58, 1200 (1980) doi:10.1139/p80-159)"
         case (NWXCP_ACM)
           cdata = "ACM (AD Becke, "
     +             //"J.Chem.Phys. 98, 5648 (1993) "
     +             //"doi:10.1063/1.464913)"
         case (NWXCP_B1B95)
           cdata = "B1B95 (AD Becke, "
     +             //"J.Chem.Phys. 104, 1040 (1996) "
     +             //"doi:10.1063/1.470829)"
         case (NWXCP_B2PLYP)
           cdata = "B2PLYP (S Grimme, "
     +             //"J.Chem.Phys. 124, 034108 (2006) "
     +             //"doi:10.1063/1.2148954)"
         case (NWXCP_B3LYP)
           cdata = "B3LYP (PJ Stephens, FJ Devlin, CF Chabalowski, "
     +             //"MJ Frisch, J.Phys.Chem. 98, 11623 (1994) "
     +             //"doi:10.1021/j100096a001)"
         case (NWXCP_B3P86)
           cdata = "B3P86 (unpublished analog of B3LYP)"
         case (NWXCP_B3PW91)
           cdata = "B3PW91 (unpublished analog of B3LYP)"
         case (NWXCP_BLYP)
           cdata = "BLYP (unpublished)"
         case (NWXCP_B0LYP)
           cdata = "B0LYP (ACM0 variant of BLYP)"
         case (NWXCP_B97)
           cdata = "B97 (AD Becke, "
     +             //"J.Chem.Phys. 107, 8554 (1997)"
     +             //" doi:10.1063/1.475007)"
         case (NWXCP_B97_1)
           cdata = "B97-1 (FA Hamprecht, A Cohen, "
     +             //"DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)"
     +             //" doi:10.1063/1.477267)"
         case (NWXCP_B97_2)
           cdata = "B97-2 (PJ Wilson, TJ Bradley, "
     +             //"DJ Tozer, J.Chem.Phys. 115, 9233 (2001)"
     +             //" doi:10.1063/1.1412605)"
         case (NWXCP_B97_3)
           cdata = "B97-3 (TW Keal, DJ Tozer, "
     +             //"J.Chem.Phys. 123, 121103 (2005)"
     +             //" doi:10.1063/1.2061227)"
         case (NWXCP_B97_D)
           cdata = "B97-D (S Grimme, "
     +             //"J.Comp.Chem. 27, 1787 (2006)"
     +             //" doi:10.1002/jcc.20495)"
         case (NWXCP_B97_G)
           cdata = "B97-GGA1 (AJ Cohen, NC Handy, "
     +             //"Chem.Phys.Lett. 316, 160 (2000)"
     +             //" doi:10.1016/S0009-2614(99)01273-7)"
         case (NWXCP_B98)
           cdata = "B98 (HL Schmider, AD Becke, "
     +             //"J.Chem.Phys. 108, 9624 (1998)"
     +             //" doi:10.1063/1.476438)"
         case (NWXCP_BB1K)
           cdata = "BB1K (Y Zhao, BJ Lynch, DG Truhlar, "
     +             //"J.Phys.Chem. A 108, 2715 (2004) "
     +             //"doi:10.1021/jp049908s)"
         case (NWXCP_BHANDH)
           cdata = "BECKEHANDH (AD Becke, "
     +             //"J.Chem.Phys. 98, 1372 (1993)"
     +             //" doi:10.1063/1.464304)"
         case (NWXCP_BHLYP)
           cdata = "BHLYP (BHANDHLYP, Gaussian inc.)"
         case (NWXCP_BOP)
           cdata = "BOP (T Tsuneda, T Suzumura, K Hirao, "
     +             //"J.Chem.Phys. 110, 10664 (1999)"
     +             //" doi:10.1063/1.479012)"
         case (NWXCP_CAMS12G)
           cdata = "CAM-S12g (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_CAMS12H)
           cdata = "CAM-S12h (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_DLDF)
           cdata = "DLDF (K Pernal, R Podeszwa, "
     +             //"K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, "
     +             //"263201 (2009) doi:10.1103/PhysRevLett.103.263201)"
         case (NWXCP_FT97)
           cdata = "FT97 (M Filatov, W Thiel, "
     +             //"Mol.Phys. 91, "
     +             //"847 (1997) doi:10.1080/002689797170950)"
         case (NWXCP_HCTH)
           cdata = "HCTH (FA Hamprecht, A Cohen, "
     +             //"DJ Tozer, NC Handy, J.Chem.Phys. 109, 6264 (1998)"
     +             //" doi:10.1063/1.477267)"
         case (NWXCP_HCTH120)
           cdata = "HCTH120 (AD Boese, NL Doltsinis, "
     +             //"NC Handy, M Sprik, J.Chem.Phys. 112, 1670 (2000)"
     +             //" doi:10.1063/1.480732)"
         case (NWXCP_HCTH147)
           cdata = "HCTH147 (AD Boese, NL Doltsinis, "
     +             //"NC Handy, M Sprik, J.Chem.Phys. 112, 1670 (2000)"
     +             //" doi:10.1063/1.480732)"
         case (NWXCP_HCTH407)
           cdata = "HCTH407 (AD Boese, "
     +             //"NC Handy, J.Chem.Phys. 114, 5497 (2001)"
     +             //" doi:10.1063/1.1347371)"
         case (NWXCP_HCTH407P)
           cdata = "HCTH407P (AD Boese, A Chandra, "
     +             //"JML Martin, D Marx, J.Chem.Phys. 119, 5965 (2003)"
     +             //" doi:10.1063/1.1599338)"
         case (NWXCP_HCTH_A)
           cdata = "HCTH(AC) (DJ Tozer, NC Handy, "
     +             //"J.Chem.Phys. 109, 10180 (1998)"
     +             //" doi:10.1063/1.477711)"
         case (NWXCP_HCTHP14)
           cdata = "HCTHP14 (G Menconi, PJ Wilson, "
     +             //"DJ Tozer, J.Chem.Phys. 114, 3958 (2001)"
     +             //" doi:10.1063/1.1342776)"
         case (NWXCP_KT1)
           cdata = "KT1 (TW Keal, "
     +             //"DJ Tozer, J.Chem.Phys. 119, 3015 (2003)"
     +             //" doi:10.1063/1.1590634)"
         case (NWXCP_KT2)
           cdata = "KT2 (TW Keal, "
     +             //"DJ Tozer, J.Chem.Phys. 119, 3015 (2003)"
     +             //" doi:10.1063/1.1590634)"
         case (NWXCP_LDA)
           cdata = "LDA SVWN5 (SH Vosko, L Wilk, M Nusair, "
     +             //"Can.J.Phys. 58, 1200 (1980) doi:10.1139/p80-159)"
         case (NWXCP_M05)
           cdata = "M05 (Y Zhao, NE Schultz, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
         case (NWXCP_M05_2X)
           cdata = "M05-2X (Y Zhao, NE Schultz, DG Truhlar,"
     +             //" J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
         case (NWXCP_M06)
           cdata = "M06 (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_M06_HF)
           cdata = "M06-HF (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_M06_L)
           cdata = "M06-L (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_M06_2X)
           cdata = "M06-2X (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 194101 (2006)"
     +             //" doi:10.1063/1.2370993)"
         case (NWXCP_M08_HX)
           cdata = "M08-HX (Y Zhao, DG Truhlar, "
     +             //"J.Chem.TheoryComput. 4, 1849 (2008)"
     +             //" doi:10.1021/ct800246v)"
         case (NWXCP_M08_SO)
           cdata = "M08-SO (Y Zhao, DG Truhlar, "
     +             //"J.Chem.TheoryComput. 4, 1849 (2008)"
     +             //" doi:10.1021/ct800246v)"
         case (NWXCP_M11)
           cdata = "M11 (R Peverati, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 2810 (2011)"
     +             //" doi:10.1021/jz201170d)"
         case (NWXCP_M11_L)
           cdata = "M11-L (R Peverati, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 2810 (2011)"
     +             //" doi:10.1021/jz201170d)"
         case (NWXCP_MPW1K)
           cdata = "MPW1K (BJ Lynch, PL Fast, M Harris, DG Truhlar, "
     +             //"J.Phys.Chem. A 104, 4811 (2000)"
     +             //" doi:10.1021/jp000497z)"
         case (NWXCP_MPW1B95)
           cdata = "MPW1B95 (AD Becke, "
     +             //"J.Chem.Phys. 104, 1040 (1996)"
     +             //" doi:10.1063/1.470829)"
         case (NWXCP_MPWB1K)
           cdata = "MPWB1K (AD Becke, "
     +             //"J.Chem.Phys. 104, 1040 (1996)"
     +             //" doi:10.1063/1.470829)"
         case (NWXCP_HSE)
           cdata = "HSE (J Heyd, GE Scuseria, M Ernzerhof, "
     +             //"J.Chem.Phys. 118, 8207 (2003)"
     +             //" doi:10.1063/1.2204597)"
         case (NWXCP_OPT)
           cdata = "OPT (NC Handy, AJ Cohen, "
     +             //"Mol.Phys. 99, 607 (2001)"
     +             //" doi:10.1080/00268970010023435)"
         case (NWXCP_PBE)
           cdata = "PBE (JP Perdew, K Burke, M Ernzerhof, "
     +             //"Phys.Rev.Lett. 77, 3865 (1996)"
     +             //" doi:10.1103/PhysRevLett.77.3865)"
         case (NWXCP_PBE0)
           cdata = "PBE0 (C Adamo, V Barone, "
     +             //"J.Chem.Phys. 110, 6158 (1999)"
     +             //" doi:10.1063/1.478522)"
         case (NWXCP_PBEOP)
           cdata = "PBEOP (T Tsuneda, T Suzumura, K Hirao, "
     +             //"J.Chem.Phys. 111, 5656 (1999)"
     +             //" doi:10.1063/1.479954)"
         case (NWXCP_PBESOL)
           cdata = "PBESOL (JP Perdew, A Ruzsinszky, "
     +             //"GI Csonka, et al., "
     +             //"Phys.Rev.Lett. 100, 136406 (2008)"
     +             //" doi:10.1103/PhysRevLett.100.136406)"
         case (NWXCP_REVPBE)
           cdata = "revPBE (Y Zhang, W Yang, "
     +             //"Phys.Rev.Lett. 80, 890 (1998)"
     +             //" doi:10.1103/PhysRevLett.80.890)"
         case (NWXCP_REVPBE0)
           cdata = "revPBE0 (ACM0 variant of revPBE)"
         case (NWXCP_RPBE)
           cdata = "RPBE (B Hammer, LB Hanssen, JK Norskov, "
     +             //"Phys.Rev.B 59, 7413 (1999)"
     +             //" doi:10.1103/PhysRevB.59.7413)"
         case (NWXCP_PKZB)
           cdata = "PKZB (JP Perdew, S Kurth, A Zupan, "
     +             //"P Blaha, Phys.Rev.Lett. 82, "
     +             //"2544 (1999) doi:10.1103/PhysRevLett.82.2544)"
         case (NWXCP_PW6B95)
           cdata = "PW6B95 (Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem. A 109, 5656 (2005)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_PWB6K)
           cdata = "PWB6K (Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem. A 109, 5656 (2005)"
     +             //" doi:10.1021/jp050536c)"
         case (NWXCP_S12G)
           cdata = "S12g (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_S12H)
           cdata = "S12h (M Swart, "
     +             //"Chem.Phys.Lett. 580, 166 (2013)"
     +             //" doi:10.1016/j.cplett.2013.06.045)"
         case (NWXCP_SOGGA)
           cdata = "SOGGA (Y Zhao, DG Truhlar, "
     +             //"J.Chem.Phys. 128, 184109 (2008)"
     +             //" doi:10.1063/1.2912068)"
         case (NWXCP_SOGGA11)
           cdata = "SOGGA11 (R Peverati, Y Zhao, DG Truhlar, "
     +             //"J.Phys.Chem.Lett. 2, 1991 (2011)"
     +             //" doi:10.1021/jz200616w)"
         case (NWXCP_SOGGA11_X)
           cdata = "SOGGA11-X (R Peverati, DG Truhlar, "
     +             //"J.Chem.Phys. 135, 191102 (2011)"
     +             //" doi:10.1063/1.3663871)"
         case (NWXCP_SOP)
           cdata = "SOP (T Tsuneda, T Suzumura, K Hirao, "
     +             //"J.Chem.Phys. 111, 5656 (1999)"
     +             //" doi:10.1063/1.479954)"
         case (NWXCP_SSB_D)
           cdata = "SSB-D (M Swart, M Sola, FM Bickelhaupt, "
     +             //"J.Chem.Phys. 131, 094103 (2009)"
     +             //" doi:10.1063/1.3213193)"
         case (NWXCP_TPSS)
           cdata = "TPSS (J Tao, JP Perdew, VN "
     +             //"Staveroverov, GE Scuseria, Phys.Rev.Lett. 91, "
     +             //"146401 (2003) doi:10.1103/PhysRevLett.91.146401)"
         case (NWXCP_TPSSH)
           cdata = "TPSSH (VN Staveroverov, GE Scuseria, J Tao, JP "
     +             //"Perdew, J.Chem.Phys. 119, "
     +             //"12129 (2003) doi:10.1063/1.1626543)"
         case (NWXCP_VS98)
           cdata = "VS98 (T van Voorhis, GE Scuseria, "
     +             //"J.Chem.Phys. 109, 400 (1998)"
     +             //" doi:10.1063/1.3005348)"
         case (NWXCP_SCAN)
           cdata = "SCAN (T van Voorhis, GE Scuseria, "
     +             //"J.Chem.Phys. 109, 400 (1998)"
     +             //" doi:10.1063/1.3005348)"
         case (NWXCP_PBE_D1)
           cdata = "PBE-D1 (S Grimme, "
     +             //"J.Comput.Chem. 25, 1463 (2004)"
     +             //" doi:10.1002/jcc.20078)"
         case default
           stat = .false.
       end select
       i = lendat
       do while (cdata(i:i).eq." ")
         i = i - 1
       enddo
       stat = stat .and. (i.lt.lendat)
       nwpwxc_get_info = stat
c
       return
       end
C>
C> @}
c $Id$
