C> \ingroup nwxc
C> @{
C>
C> \file nwxc_query.F
C> Routines to query the properties of the current functional
C> 
C> @}
C
C> \ingroup nwxc_api
C> @{
C>
C> \brief Is this library used to evaluate the functional
C>
C> Certainly in interim stages one might want to choose to use
C> this library or another mechanism to evaluate the functional.
C> This function returns .true. if this library has been activated.
C>
      logical function nwxc_is_on()
      implicit none
#include "nwxcP.fh"
      nwxc_is_on = nwxc_on
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
      logical function nwxc_has_hfx()
      implicit none
#include "nwxcP.fh"
      nwxc_has_hfx = nwxc_wght_hf .gt. 0.0d0
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
      logical function nwxc_has_mp2c()
      implicit none
#include "nwxcP.fh"
      nwxc_has_mp2c = nwxc_wght_mp2 .gt. 0.0d0
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
      logical function nwxc_has_cam()
      implicit none
#include "nwxcP.fh"
      nwxc_has_cam = nwxc_cam_gamma .ne. 0.0d0 .or.
     +               nwxc_cam_beta  .gt. 0.0d0 .or.
     +               nwxc_cam_alpha .gt. 0.0d0
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
      logical function nwxc_has_disp()
      implicit none
#include "nwxcP.fh"
#include "nwxc_param.fh"
      nwxc_has_disp = nwxc_disp_id.ne.NWXC_DISP_NONE
      return
      end
C>
C> \brief Retrieve the density tolerance
C>
C> Retrieve the current setting of the density tolerance.
C> In the evaluation of the functional densities less 
C> than this tolerance will be ignored.
C>
      double precision function nwxc_get_rho_tolerance()
      implicit none
#include "nwxcP.fh"
      nwxc_get_rho_tolerance = nwxc_rho_tol
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
      subroutine nwxc_get_cam(alpha,beta,gamma,srhfx)
      implicit none
#include "nwxcP.fh"
      double precision alpha !< [Output] The alpha parameter
      double precision beta  !< [Output] The beta parameter
      double precision gamma !< [Output] The gamma parameter
      logical srhfx !< [Output] T if short-range exchange is 
                    !< Hartree-Fock exchange
      alpha = nwxc_cam_alpha
      beta  = nwxc_cam_beta
      gamma = nwxc_cam_gamma
      srhfx = nwxc_cam_srhfex
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
      subroutine nwxc_get_disp(dispid,s6,s8,alpha)
      implicit none
#include "nwxcP.fh"
      integer dispid      !< [Output] The dispersion identifier
      double precision s6 !< [Output] The \f$S_6\f$ scale factor
      double precision s8 !< [Output] The \f$S_8\f$ scale factor
      double precision alpha !< [Output] The \f$\alpha\f$ damping factor
      dispid = nwxc_disp_id
      s6     = nwxc_disp_s6
      s8     = nwxc_disp_s8
      alpha  = nwxc_disp_alpha
      end
C>
C> \brief Retrieve the Hartree-Fock exchange weight
C>
C> Simply look up and return the weight for the Hartree-Fock exchange.
C>
C> \return Returns the weight of the Hartree-Fock exchange term
C>
      double precision function nwxc_wght_hfx()
      implicit none
#include "nwxcP.fh"
      nwxc_wght_hfx = nwxc_wght_hf
      return
      end
C>
C> \brief Retrieve the MP2 correlation weight
C>
C> Simply look up and return the weight for the MP2 correlation term.
C>
C> \return Returns the weight of the MP2 correlation term
C>
      double precision function nwxc_wght_mp2c()
      implicit none
#include "nwxcP.fh"
      nwxc_wght_mp2c = nwxc_wght_mp2
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
      logical function nwxc_has_xc()
      implicit none
      logical  nwxc_is_lda, nwxc_is_gga, nwxc_is_mgga
      external nwxc_is_lda, nwxc_is_gga, nwxc_is_mgga
      nwxc_has_xc = nwxc_is_lda() .or. nwxc_is_gga() .or. nwxc_is_mgga()
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
C>    if (nwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density only,
C> and .false. otherwise.
C>
      logical function nwxc_is_lda()
      implicit none
#include "nwxcP.fh"
      nwxc_is_lda = nwxc_depends_on_rho .and. 
     +              (.not. nwxc_depends_on_gamma) .and.
     +              (.not. nwxc_depends_on_tau)
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
C>    if (nwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwxc_is_lda()) then
C>      ... do something ...
C>    else
C>      ... bomb ...
C>    endif
C> \endcode
C>
C> \return Return .true. if the functional depends on the density and
C> the density gradient, and .false. otherwise.
C>
      logical function nwxc_is_gga()
      implicit none
#include "nwxcP.fh"
      nwxc_is_gga = nwxc_depends_on_rho .and. 
     +              nwxc_depends_on_gamma .and.
     +              (.not. nwxc_depends_on_tau)
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
C>    if (nwxc_is_mgga()) then
C>      ... do something ...
C>    else if (nwxc_is_gga()) then
C>      ... do something ...
C>    else if (nwxc_is_lda()) then
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
      logical function nwxc_is_mgga()
      implicit none
#include "nwxcP.fh"
      nwxc_is_mgga = nwxc_depends_on_rho .and. 
     +               nwxc_depends_on_gamma .and.
     +               nwxc_depends_on_tau
      return
      end
C> @}
C
C> \ingroup nwxc_priv
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
C>    do while (nwxc_get_df(ientry,func_id,weight))
C>      write(*,*)"functional id: ",func_id," weight=",weight
C>    enddo
C> \endcode
C>
C> \return Return .true. if the entry number is greater than or equal to 1 and
C> and less than or equal to the number of density functional terms, and .false.
C> otherwise.
C>
       logical function nwxc_get_df(ientry,func_id,weight)
       implicit none
#include "nwxcP.fh"
       integer ientry  !< [In/Output] On input the entry number to extract,
                       !< on output the input value incremented by one
       integer func_id !< [Output] The functional identifier
       double precision weight !< [Output] The functional weight
c
       nwxc_get_df = .true.
       if (.not.((ientry.ge.1).and.(ientry.le.nwxc_usr_num))) then
         nwxc_get_df = .false.
         return
       endif
       func_id = nwxc_usr_ids(ientry)
       weight  = nwxc_usr_wghts(ientry)
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
       logical function nwxc_get_info(func_id,cdata)
       implicit none
#include "nwxcP.fh"
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
         case (NWXCP_X_B88)
           cdata = "Becke88 Exchange (AD Becke, "
     +             //"Phys.Rev.A 38, 3098 (1988)"
     +             //" doi:10.1103/PhysRevA.38.3098)"
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
         case (NWXCP_X_DLDF)
           cdata = "DLDF Exchange (K Pernal, R Podeszwa, "
     +             //"K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, "
     +             //"263201 (2009) doi:10.1103/PhysRevLett.103.263201)"
         case (NWXCP_X_M05)
           cdata = "M05 Exchange (Y Zhao, NE Schultz, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
         case (NWXCP_X_M05_2X)
           cdata = "M05-2X Exchange (Y Zhao, NE Schultz, DG Truhlar, "
     +             //"J.Chem.Phys. 123, 161103 (2005)"
     +             //" doi:10.1063/1.2126975)"
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
         case (NWXCP_X_TPSS)
           cdata = "TPSS Exchange (J Tao, JP Perdew, VN Staveroverov, "
     +             //"GE Scuseria, Phys.Rev.Lett. 91, 146401 (2003) "
     +             //"doi:10.1103/PhysRevLett.91.146401)"
         case (NWXCP_X_WPBE)
           cdata = "WPBE Exchange (E Weintraub, TM Henderson, "
     +             //"GE Scuseria, J.Chem.Theory.Comput. 5, 754 (2009)"
     +             //" doi:10.1021/ct800530u)"
         case (NWXCP_X_VS98)
           cdata = "VS98 Exchange (T van Voorhis, GE Scuseria, "
     +             //"J.Chem.Phys. 109, 400 (1998)"
     +             //" doi:10.1063/1.3005348)"
         case (NWXCP_C_B95)
           cdata = "B95 Correlation (AD Becke, "
     +             //"J.Chem.Phys. 104, 1040 (1996)"
     +             //" doi:10.1063/1.4708298)"
         case (NWXCP_C_DLDF)
           cdata = "DLDF Correlation (K Pernal, R Podeszwa, "
     +             //"K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, "
     +             //"263201 (2009) doi:10.1103/PhysRevLett.103.263201)"
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
         case (NWXCP_C_PBE)
           cdata = "PBE Correlation (JP Perdew, K Burke, M Ernzerhof, "
     +             //"Phys.Rev.Lett. 77, 3865 (1996)"
     +             //" doi:10.1103/PhysRevLett.77.3865)"
         case (NWXCP_C_PBESOL)
           cdata = "PBESOL Correlation (JP Perdew, A Ruzsinszky, "
     +             //"GI Csonka, et al., "
     +             //"Phys.Rev.Lett. 100, 136406 (2008)"
     +             //" doi:10.1103/PhysRevLett.100.136406)"
         case (NWXCP_C_TPSS)
           cdata = "TPSS Correlation (J Tao, JP Perdew, VN "
     +             //"Staveroverov, GE Scuseria, Phys.Rev.Lett. 91, "
     +             //"146401 (2003) doi:10.1103/PhysRevLett.91.146401)"
         case (NWXCP_C_VS98)
           cdata = "VS98 Correlation (T van Voorhis, GE Scuseria, "
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
         case (NWXCP_B3LYP)
           cdata = "B3LYP (PJ Stephens, FJ Devlin, CF Chabalowski, "
     +             //"MJ Frisch, J.Phys.Chem. 98, 11623 (1994) "
     +             //"doi:10.1021/j100096a001)"
         case (NWXCP_B3PW91)
           cdata = "B3PW91 (AD Becke, J.Chem.Phys. 98, 5648 (1993) "
     +             //"doi:10.1063/1.464913)"
         case (NWXCP_BLYP)
           cdata = "BLYP"
         case (NWXCP_B0LYP)
           cdata = "B0LYP (ACM0 variant of BLYP)"
         case (NWXCP_DLDF)
           cdata = "DLDF (K Pernal, R Podeszwa, "
     +             //"K Patkowski, K Szalewicz, Phys.Rev.Lett. 103, "
     +             //"263201 (2009) doi:10.1103/PhysRevLett.103.263201)"
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
         case (NWXCP_HSE)
           cdata = "HSE (J Heyd, GE Scuseria, M Ernzerhof, "
     +             //"J.Chem.Phys. 118, 8207 (2003)"
     +             //" doi:10.1063/1.2204597)"
         case (NWXCP_PBE)
           cdata = "PBE (JP Perdew, K Burke, M Ernzerhof, "
     +             //"Phys.Rev.Lett. 77, 3865 (1996)"
     +             //" doi:10.1103/PhysRevLett.77.3865)"
         case (NWXCP_PBE0)
           cdata = "PBE0 (C Adamo, V Barone, "
     +             //"J.Chem.Phys. 110, 6158 (1999)"
     +             //" doi:10.1063/1.478522)"
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
         case (NWXCP_TPSS)
           cdata = "TPSS (J Tao, JP Perdew, VN "
     +             //"Staveroverov, GE Scuseria, Phys.Rev.Lett. 91, "
     +             //"146401 (2003) doi:10.1103/PhysRevLett.91.146401)"
         case (NWXCP_VS98)
           cdata = "VS98 (T van Voorhis, GE Scuseria, "
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
       nwxc_get_info = stat
c
       return
       end
C>
C> @}
