C> \ingroup nwxc
C> @{
C>
C> \file nwxc_add.F
C> Routines to add new terms to the functional
C>
C> @}
C
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Add a functional term to the evaluation list
C>
C> This function is for internal use of the library.
C>
C> The evaluation list is the list of functional terms as they will
C> be used to compute the functional values. This routine adds a term
C> identified by an integer to this list. The routine also takes a 
C> weight for this term.
C>
C> It checks whether the term is already on the list. If so the weight
C> for the entry on list is increased otherwise a new entry is added.
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function nwxc_add_func(iterm,weight)
      implicit none
#include "nwxcP.fh"
      integer iterm           !< [Input] The functional term identifier
      double precision weight !< [Input] The functional term weight
c
      integer i
c
      nwxc_add_func = .false.
c
      i = 1 
      do while (i.le.nwxc_num_f.and.nwxc_ids(i).ne.iterm)
        i = i + 1
      enddo
      if (i.gt.nwxc_num_f) then
c       the functional was not found on the list
        if (nwxc_num_f.ge.NWXC_MAX_F) then
c         we are out of space
          return
        endif
        nwxc_num_f = nwxc_num_f + 1
        nwxc_ids(nwxc_num_f)   = iterm
        nwxc_wghts(nwxc_num_f) = weight
      else
c       the functional was found on the list
        nwxc_wghts(i) = nwxc_wghts(i) + weight
      endif
      nwxc_add_func = .true.
c
      return
      end
C> 
C>
C> \brief Add a functional term to the user list
C>
C> This function is for internal use of the library.
C>
C> The user list is the list of functional terms as they were provided
C> by the user in the input. This routine adds a term
C> identified by an integer to this list. The routine also takes a 
C> weight for this term.
C>
C> It checks whether the term is already on the list. If so the weight
C> for the entry on list is increased otherwise a new entry is added.
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function nwxc_usr_add_func(iterm,weight)
      implicit none
#include "nwxcP.fh"
      integer iterm           !< [Input] The functional term identifier
      double precision weight !< [Input] The functional term weight
c
      integer i
c
      nwxc_usr_add_func = .true.
c
      if (iterm.eq.NWXCP_HFXC) then
        return
      else if (iterm.eq.NWXCP_MP2C) then
        return
      endif
c
      i = 1 
      do while (i.le.nwxc_usr_num.and.nwxc_usr_ids(i).ne.iterm)
        i = i + 1
      enddo
      if (i.gt.nwxc_usr_num) then
c       the functional was not found on the list
        if (nwxc_usr_num.ge.NWXC_MAX_F) then
c         we are out of space
          nwxc_usr_add_func = .false.
          return
        endif
        nwxc_usr_num = nwxc_usr_num + 1
        nwxc_usr_ids(nwxc_usr_num)   = iterm
        nwxc_usr_wghts(nwxc_usr_num) = weight
      else
c       the functional was found on the list
        nwxc_usr_wghts(i) = nwxc_usr_wghts(i) + weight
      endif
c
      return
      end
C> 
C> \brief Add a functional term to both the user and evaluation list
C>
C> This function is for internal use of the library.
C>
C> If the user specifies terms that are not aliases the same term
C> has to be added to both the user and the evaluation list. As this
C> is likely a common occurance this routine provided the convenience
C> of doing this with one call.
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function nwxc_all_add_func(iterm,weight)
      implicit none
      integer iterm           !< [Input] The functional term identifier
      double precision weight !< [Input] The functional term weight
c
      logical  nwxc_add_func, nwxc_usr_add_func
      external nwxc_add_func, nwxc_usr_add_func
c
      nwxc_all_add_func = nwxc_add_func(iterm,weight).and.
     +                    nwxc_usr_add_func(iterm,weight)
c
      return
      end
C>
C> @}
C
C> \ingroup nwxc_api
C> @{
C>
C> \brief Add a functional by name to the current functional definition
C>
C> This function adds a named functional to the current functional
C> definition. It takes both the name and the corresponding weight 
C> as arguments. If necessary this routine will resolve aliases.
C>
C> \return Returns .true. if successfull, and .false. otherwise.
C>
      logical function nwxc_add_df(dfname,dfweight)
      implicit none
#include "nwxcP.fh"
#include "nwxc_param.fh"
      character*(*)    dfname   !< [Input] The name of the functional
      double precision dfweight !< [Input] The weight of the functional
c
      logical  nwxc_add_func, nwxc_usr_add_func, nwxc_all_add_func
      external nwxc_add_func, nwxc_usr_add_func, nwxc_all_add_func
c
      nwxc_add_df = .false.
      if (dfname.eq."hfexch") then
        nwxc_wght_hf = nwxc_wght_hf + dfweight
        nwxc_add_df = .true.
c       Hartree-Fock exchange and MP2 correlation are handled different
c       nwxc_add_df = nwxc_usr_add_func(NWXCP_HFXC,dfweight)
      else if (dfname.eq."mp2") then
        nwxc_wght_mp2 = nwxc_wght_mp2 + dfweight
        nwxc_add_df = .true.
c       Hartree-Fock exchange and MP2 correlation are handled different
c       nwxc_add_df = nwxc_usr_add_func(NWXCP_MP2C,dfweight)
      else if (dfname.eq."slater") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_SLATER,dfweight)
      else if (dfname.eq."becke88") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B88,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,dfweight)
      else if (dfname.eq."becke0") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B0,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.75d0*dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.25d0*dfweight
      else if (dfname.eq."optx") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_OPT,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,1.05151d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_OPT,1.43169d0*dfweight)
      else if (dfname.eq."xbnl07".or.dfname.eq."bnl") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_BNL,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 1.0d0
        nwxc_cam_gamma  = 0.5d0
      else if (dfname.eq."xcamb88") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_CAMB88,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_CAMB88,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.19d0
        nwxc_cam_beta   = 0.46d0
        nwxc_cam_gamma  = 0.33d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xcamlsd") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_CAMLSD,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.5d0
        nwxc_cam_beta   = 0.5d0
        nwxc_cam_gamma  = 0.3d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xcampbe".or.dfname.eq."xcampbe96") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_CAMPBE,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + 0.25d0*dfweight
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 1.0d0
        nwxc_cam_gamma  = 0.33d0
        nwxc_cam_srhfex = .true.
      else if (dfname.eq."xdldf") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_DLDF,dfweight)
      else if (dfname.eq."xm05") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M05,dfweight)
      else if (dfname.eq."xm05-2x") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M05_2X,dfweight)
      else if (dfname.eq."xm06") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M06,dfweight)
      else if (dfname.eq."xm06-hf") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M06_HF,dfweight)
      else if (dfname.eq."xm06-l") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M06_L,dfweight)
      else if (dfname.eq."xm06-2x") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M06_2X,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.54*dfweight
      else if (dfname.eq."xm08-hx") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M08_HX,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.5223d0*dfweight
      else if (dfname.eq."xm08-so") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M08_SO,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.5679d0*dfweight
      else if (dfname.eq."xm11") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M11,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 1.0000d0*dfweight
        nwxc_cam_alpha  = 0.428d0
        nwxc_cam_beta   = 0.572d0
        nwxc_cam_gamma  = 0.25d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xm11-l") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_M11_L,dfweight)
      else if (dfname.eq."xpbe96") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,dfweight)
      else if (dfname.eq."xpbe0") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_PBE0,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,0.75d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.25d0*dfweight
      else if (dfname.eq."xpbesol") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_PBESOL,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBESOL,dfweight)
      else if (dfname.eq."xrevpbe") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_REVPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_REVPBE,dfweight)
      else if (dfname.eq."xrevpbe0") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_REVPBE0,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_REVPBE,0.75d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.25d0*dfweight
      else if (dfname.eq."xrpbe") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_RPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_RPBE,dfweight)
      else if (dfname.eq."xpkzb99") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_PKZB,dfweight)
      else if (dfname.eq."xtpss03") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_TPSS,dfweight)
      else if (dfname.eq."xwpbe") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_WPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_WPBE,dfweight)
c       nwxc_wght_hf    = 1.00d0
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 1.0d0
        nwxc_cam_gamma  = 0.30d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xvs98") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_VS98,dfweight)
      else if (dfname.eq."cdldf") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_DLDF,dfweight)
      else if (dfname.eq."cpw91lda".or.dfname.eq."pw91lda") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PW91LDA,dfweight)
      else if (dfname.eq."lyp") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_LYP,dfweight)
      else if (dfname.eq."cm05") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M05,dfweight)
      else if (dfname.eq."cm05-2x") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M05_2X,dfweight)
      else if (dfname.eq."cm06") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M06,dfweight)
      else if (dfname.eq."cm06-hf") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M06_HF,dfweight)
      else if (dfname.eq."cm06-l") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M06_L,dfweight)
      else if (dfname.eq."cm06-2x") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M06_2X,dfweight)
      else if (dfname.eq."cm08-hx") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M08_HX,dfweight)
      else if (dfname.eq."cm08-so") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M08_SO,dfweight)
      else if (dfname.eq."cm11") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M11,dfweight)
      else if (dfname.eq."cm11-l") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_M11_L,dfweight)
      else if (dfname.eq."cpbe96") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."cpbesol") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PBESOL,dfweight)
      else if (dfname.eq."cpkzb99") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PKZB,dfweight)
      else if (dfname.eq."ctpss03") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_TPSS,dfweight)
      else if (dfname.eq."cvs98") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VS98,dfweight)
      else if (dfname.eq."perdew91") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_C_P91,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW91LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_P91,dfweight)
      else if (dfname.eq."vwn_1") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VWN1,dfweight)
      else if (dfname.eq."vwn_1_rpa") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VWN1_RPA,dfweight)
      else if (dfname.eq."vwn_2") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VWN2,dfweight)
      else if (dfname.eq."vwn_3") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VWN3,dfweight)
      else if (dfname.eq."vwn_4") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VWN4,dfweight)
      else if (dfname.eq."vwn_5") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VWN5,dfweight)
      else if (dfname.eq."bc95") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B95,dfweight)
      else if (dfname.eq."cpw6b95") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PW6B95,dfweight)
      else if (dfname.eq."cpwb6k") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PWB6K,dfweight)
      else if (dfname.eq."blyp") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_BLYP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_LYP,dfweight)
      else if (dfname.eq."b0lyp".or.dfname.eq."blyp0") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B0LYP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_LYP,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.25d0*dfweight
      else if (dfname.eq."b3lyp") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B3LYP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.80d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VWN1_RPA,0.19d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_LYP,0.81d0*dfweight)
        nwxc_wght_hf    = 0.20d0
c     else if (dfname.eq."b3pw91") then ! alias
c       nwxc_add_df = nwxc_usr_add_func(NWXCP_B3PW91,dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_X_SLATER,0.80d0*dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_X_B88,0.72d0*dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_C_VWN5,0.19d0*dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_C_PW91,0.81d0*dfweight)
c       nwxc_wght_hf    = 0.20d0
      else if (dfname.eq."dldf") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_DLDF,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_DLDF,0.3855871d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_DLDF,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.6144129d0*dfweight
      else if (dfname.eq."lda") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VWN5,dfweight)
      else if (dfname.eq."m05") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M05,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M05,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M05,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.28d0*dfweight
      else if (dfname.eq."m05-2x") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M05_2X,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M05_2X,0.44d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M05_2X,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.56d0*dfweight
      else if (dfname.eq."m06") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M06,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M06,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M06,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.27d0*dfweight
      else if (dfname.eq."m06-hf") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M06_HF,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M06_HF,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M06_HF,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + dfweight
      else if (dfname.eq."m06-l") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M06_L,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M06_L,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M06_L,dfweight)
      else if (dfname.eq."m06-2x") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M06_2X,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M06_2X,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M06_2X,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.54d0*dfweight
      else if (dfname.eq."m08-hx") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M08_HX,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M08_HX,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M08_HX,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.5223d0*dfweight
      else if (dfname.eq."m08-so") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M08_SO,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_X_M08_SO,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     +                nwxc_add_func(NWXCP_C_M08_SO,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.5679d0*dfweight
      else if (dfname.eq."m11") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M11,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_M11,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_M11,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 1.0000d0*dfweight
        nwxc_cam_alpha  = 0.428d0
        nwxc_cam_beta   = 0.572d0
        nwxc_cam_gamma  = 0.25d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."m11-l") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_M11_L,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_M11_L,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_M11_L,dfweight)
      else if (dfname.eq."hse03".or.dfname.eq."hse") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_HSE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_CAMPBE,-0.25d0*dfweight)
        nwxc_wght_hf    = 0.25d0
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 1.0d0
        nwxc_cam_gamma  = 0.33d0
        nwxc_cam_srhfex = .true.
      else if (dfname.eq."pbe96") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."pbe0") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PBE0,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.25d0*dfweight
      else if (dfname.eq."pbesol") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PBESOL,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBESOL,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBESOL,dfweight)
      else if (dfname.eq."revpbe") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_REVPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_REVPBE,dfweight)
c       In the NWChem DFT code REVPBE does not include correlation
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."revpbe0") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_REVPBE0,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.75d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_REVPBE,0.75d0*dfweight)
c       In the NWChem DFT code REVPBE does not include correlation
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.25d0*dfweight
      else if (dfname.eq."rpbe") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_RPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_RPBE,dfweight)
c       In the NWChem DFT code RPBE does not include correlation
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."pkzb99".or.dfname.eq."pkzb") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PKZB,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PKZB,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PKZB,dfweight)
      else if (dfname.eq."tpss03".or.dfname.eq."tpss") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_TPSS,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_TPSS,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_TPSS,dfweight)
      else if (dfname.eq."vs98") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_VS98,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_VS98,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VS98,dfweight)
c
c     Dispersion corrected functionals
c
      else if (dfname.eq."pbe96-d1") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PBE_D1,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_disp_id    = NWXC_DISP_D1
        nwxc_disp_s6    = 0.7d0
        nwxc_disp_s8    = 0.0d0
        nwxc_disp_sr6   = 0.0d0
        nwxc_disp_sr8   = 0.0d0
        nwxc_disp_alpha = 23.0d0
      else
c       Oh dear, the functional was not recognized
      endif
c
      return
      end
C>
C> \brief Set the density tolerance
C>
C> Sets the tolerance that is used to screen small density values.
C> Valid tolerances are nonnegative double precision values.
C>
C> \return Return .true. if a valid tolerance value is passed in,
C> and .false. otherwise.
C>
      logical function nwxc_set_rho_tolerance(tol)
      implicit none
#include "nwxcP.fh"
      double precision tol
      nwxc_set_rho_tolerance = .false.
      if (tol.lt.0.0d0) return
      nwxc_rho_tol = tol
      return
      end
C>
C> \brief Detect the dependencies for the functional
C>
C> An important aspect in DFT is to keep track of what quantities
C> a functional depends on as this relates to many expressions a code
C> needs to evaluate. This function updates the appropriate data
C> in accordance to the current definition of the functional.
C>
      subroutine nwxc_detect_dependencies()
      implicit none
#include "nwxcP.fh"
      integer i      !< Counter
      integer maxpt  !< The number of points
      parameter (maxpt = 3)
      integer maxpt2  !< The number of points
      parameter (maxpt2 = 2*maxpt)
      integer maxpt3  !< The number of points
      parameter (maxpt3 = 3*maxpt)
c     
      double precision rho(maxpt2)    !< The density
      double precision rgamma(maxpt3) !< The norm of the gradient
      double precision tau(maxpt2)    !< The kinetic energy density
c
      double precision func(maxpt)    !< The functional value
      double precision dfdr(maxpt2)   !< Derivative wrt rho
      double precision dfdg(maxpt3)   !< Derivative wrt rgamma
      double precision dfdt(maxpt2)   !< Derivative wrt tau
c
      double precision abssum
c
c     Data taken from the quartet Nitrogen test case nwxc_nwdft_4n;
c
      data rho   /2.2175d-3, 8.2309d-3, 1.4192d-2,
     +            1.8849d-4, 1.1904d-3, 2.4799d-3/
      data rgamma/1.7078d-5, 2.6544d-4, 8.3431d-4,
     +            2.1459d-6, 5.1828d-5, 1.9672d-4,
     +            2.6964d-7, 1.0119d-5, 4.6388d-5/
      data tau   /1.1802d-3, 5.2247d-3, 9.8571d-3,
     +            1.7882d-4, 1.0629d-3, 2.3387d-3/
      data func/maxpt*0.0d0/
      data dfdr/maxpt2*0.0d0/
      data dfdg/maxpt3*0.0d0/
      data dfdt/maxpt2*0.0d0/
c
      call nwxc_eval_df(2,maxpt,rho,rgamma,tau,func,dfdr,dfdg,dfdt)
c
      abssum = 0.0d0
      do i = 1, maxpt2
        abssum = abssum + abs(dfdr(i))
      enddo
      nwxc_depends_on_rho = (abssum.gt.0.0d0)
c
      abssum = 0.0d0
      do i = 1, maxpt3
        abssum = abssum + abs(dfdg(i))
      enddo
      nwxc_depends_on_gamma = (abssum.gt.0.0d0)
c
      abssum = 0.0d0
      do i = 1, maxpt2
        abssum = abssum + abs(dfdt(i))
      enddo
      nwxc_depends_on_tau = (abssum.gt.0.0d0)
c
      end

C> @}
c $Id$
