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
      else if (dfname.eq."srhfexch") then
c       Represent the short-range exchange with exact exchange
        nwxc_cam_srhfex = .true.
      else if (dfname.eq."lrhfexch") then
c       Represent the long-range exchange with exact exchange
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."slater") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_SLATER,dfweight)
      else if (dfname.eq."becke86b") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B86B,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B86B,dfweight)
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
      else if (dfname.eq."xbecke97") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B97,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.1943d0*dfweight
      else if (dfname.eq."xbecke97-1") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B97_1,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_1,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.2100d0*dfweight
      else if (dfname.eq."xbecke97-2") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B97_2,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_2,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.2100d0*dfweight
      else if (dfname.eq."xbecke97-3") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B97_3,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_3,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.269288d0*dfweight
      else if (dfname.eq."xbecke97-d") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_B97_D,dfweight)
      else if (dfname.eq."xbecke97gga1") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_B97_G,dfweight)
      else if (dfname.eq."xbecke98") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_B98,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B98,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.219847d0*dfweight
      else if (dfname.eq."mpw91".or.dfname.eq."xmpw91") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_MPW91,dfweight)
      else if (dfname.eq."xperdew91".or.dfname.eq."xpw91") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_PW91,dfweight)
      else if (dfname.eq."xpw6b95") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_PW6B95,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,1.00d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PW6B95,1.00d0*dfweight)
c       nwxc_wght_hf = nwxc_wght_hf + 0.28d0*dfweight
      else if (dfname.eq."xpwb6k") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_PWB6K,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,1.00d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PWB6K,1.00d0*dfweight)
c       nwxc_wght_hf = nwxc_wght_hf + 0.46d0*dfweight
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
        nwxc_cam_srhfex = .false.
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
        nwxc_wght_hf    = nwxc_wght_hf + 1.00d0*dfweight
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 1.0d0
        nwxc_cam_gamma  = 0.33d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xcamrevpbe") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_CAMREVPBE,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + 1.00d0*dfweight
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 1.0d0
        nwxc_cam_gamma  = 0.33d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xcamrpbe") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_CAMRPBE,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + 1.00d0*dfweight
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 1.0d0
        nwxc_cam_gamma  = 0.33d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xcams12g".or.dfname.eq."xcam-s12g") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_CAMS12G,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 0.34485046d0
        nwxc_cam_gamma  = 1.52420731d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xcams12h".or.dfname.eq."xcam-s12h") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_CAMS12H,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.25d0
        nwxc_cam_beta   = 0.10897845d0
        nwxc_cam_gamma  = 0.48516891d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."xdldf") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_DLDF,dfweight)
      else if (dfname.eq."xft97") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_FT97,dfweight)
        nwxc_add_df = nwxc_add_df .and. 
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and. 
     &                nwxc_add_func(NWXCP_X_FT97,dfweight)
      else if (dfname.eq."xgill".or.dfname.eq."gill96") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_GILL,dfweight)
      else if (dfname.eq."xhcth".or.dfname.eq."xhcth93") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_HCTH,dfweight)
      else if (dfname.eq."xhcth120") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_HCTH120,dfweight)
      else if (dfname.eq."xhcth147") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_HCTH147,dfweight)
      else if (dfname.eq."xhcth407") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_HCTH407,dfweight)
      else if (dfname.eq."xhcth407p") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_HCTH407P,dfweight)
c     else if (dfname.eq."xhcth-ac") then
c       nwxc_add_df = nwxc_all_add_func(NWXCP_X_HCTH_A,dfweight)
      else if (dfname.eq."xhcthp14") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_HCTHP14,dfweight)
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
      else if (dfname.eq."xs12g") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_S12G,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,
     &                              dfweight*1.03842032d0)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_S12G,dfweight)
      else if (dfname.eq."xs12h") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_S12H,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight*
     &                              0.75d0*1.02543951d0)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_S12H,dfweight*
     &                              0.75d0)
        nwxc_wght_hf= nwxc_wght_hf + 0.25d0*dfweight
      else if (dfname.eq."xsogga") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_SOGGA,dfweight)
      else if (dfname.eq."xsogga11") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_SOGGA11,dfweight)
      else if (dfname.eq."xsogga11-x") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_SOGGA11_X,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.4015d0*dfweight
      else if (dfname.eq."xssb-d") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_SSB_D,dfweight)
      else if (dfname.eq."xtpss03") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_TPSS,dfweight)
      else if (dfname.eq."xtpssh") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_TPSSH,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_TPSS,0.90d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.10d0*dfweight
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
      else if (dfname.eq."cbecke97") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B97,dfweight)
      else if (dfname.eq."cbecke97-1") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B97_1,dfweight)
      else if (dfname.eq."cbecke97-2") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B97_2,dfweight)
      else if (dfname.eq."cbecke97-3") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B97_3,dfweight)
      else if (dfname.eq."cbecke97-d") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B97_D,dfweight)
        nwxc_disp_id    = NWXC_DISP_D2
        nwxc_disp_s6    = 1.25d0
        nwxc_disp_s8    = 0.0d0
        nwxc_disp_sr6   = 0.0d0
        nwxc_disp_sr8   = 0.0d0
        nwxc_disp_alpha = 20.0d0
      else if (dfname.eq."cbecke97gga1") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B97_G,dfweight)
      else if (dfname.eq."cbecke98") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_B98,dfweight)
      else if (dfname.eq."ccams12g".or.dfname.eq."ccam-s12g") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_C_CAMS12G,dfweight)
        nwxc_add_df = nwxc_add_df .and. 
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.86124355d0
        nwxc_disp_sr6   = 1.20250451d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."ccams12h".or.dfname.eq."ccam-s12h") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_C_CAMS12H,dfweight)
        nwxc_add_df = nwxc_add_df .and. 
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.37999939d0
        nwxc_disp_sr6   = 1.08034183d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."cdldf") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_DLDF,dfweight)
      else if (dfname.eq."cft97") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_FT97,dfweight)
      else if (dfname.eq."chcth".or.dfname.eq."chcth93") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_HCTH,dfweight)
      else if (dfname.eq."chcth120") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_HCTH120,dfweight)
      else if (dfname.eq."chcth147") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_HCTH147,dfweight)
      else if (dfname.eq."chcth407") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_HCTH407,dfweight)
      else if (dfname.eq."chcth407p") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_HCTH407P,dfweight)
c     else if (dfname.eq."chcth-ac") then
c       nwxc_add_df = nwxc_all_add_func(NWXCP_C_HCTH_A,dfweight)
      else if (dfname.eq."chcthp14") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_HCTHP14,dfweight)
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
      else if (dfname.eq."cmpbe") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_MPBE,dfweight)
      else if (dfname.eq."cop".or.dfname.eq."op") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_OP,dfweight)
      else if (dfname.eq."optc") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_OPT,dfweight)
      else if (dfname.eq."cpbe96") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."cpbesol") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PBESOL,dfweight)
      else if (dfname.eq."cpkzb99") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PKZB,dfweight)
      else if (dfname.eq."cs12g") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_C_S12G,dfweight)
        nwxc_add_df = nwxc_add_df .and. 
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.84432515d0
        nwxc_disp_sr6   = 1.17755954d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."cs12h") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_C_S12H,dfweight)
        nwxc_add_df = nwxc_add_df .and. 
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.37705816d0
        nwxc_disp_sr6   = 1.07735222d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."cspbe96") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_SPBE,dfweight)
      else if (dfname.eq."csogga11") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_SOGGA11,dfweight)
      else if (dfname.eq."csogga11-x") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_SOGGA11_X,dfweight)
      else if (dfname.eq."ctpss03") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_TPSS,dfweight)
      else if (dfname.eq."cvs98") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VS98,dfweight)
      else if (dfname.eq."perdew81".or.dfname.eq."cperdew81") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PZ81,dfweight)
      else if (dfname.eq."perdew86".or.dfname.eq."cperdew86") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_P86,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PZ81,dfweight)
      else if (dfname.eq."perdew91".or.dfname.eq."cperdew91") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_C_P91,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW91LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_P91,dfweight)
      else if (dfname.eq."cp91vwn5") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_C_P91_VWN5,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VWN5,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_P91_VWN5,dfweight)
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
      else if (dfname.eq."acm") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_ACM,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.80d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VWN5,1.00d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_P91_VWN5,0.81d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.20d0*dfweight
      else if (dfname.eq."beckehandh") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_BHANDH,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.5d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW91LDA,0.5d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.5d0*dfweight
      else if (dfname.eq."bhlyp") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_BHLYP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.5d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.5d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_LYP,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.5d0*dfweight
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
      else if (dfname.eq."b1b95") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B1B95,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B95,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.28d0*dfweight
      else if (dfname.eq."b2plyp") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B2PLYP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.47d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.47d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_LYP,0.73d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.53d0*dfweight
        nwxc_wght_mp2 = nwxc_wght_mp2 + 0.27d0*dfweight
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
        nwxc_wght_hf= nwxc_wght_hf + 0.20d0*dfweight
      else if (dfname.eq."b3p86") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B3P86,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.80d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VWN1_RPA,1.00d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_P86,0.81d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.20d0*dfweight
      else if (dfname.eq."b3pw91") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B3PW91,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.80d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW91LDA,1.00d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_P91,0.81d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.20d0*dfweight
      else if (dfname.eq."becke97") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B97,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B97,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.1943d0*dfweight
      else if (dfname.eq."becke97-1") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B97_1,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_1,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B97_1,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.2100d0*dfweight
      else if (dfname.eq."becke97-2") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B97_2,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_2,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B97_2,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.2100d0*dfweight
      else if (dfname.eq."becke97-3") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B97_3,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_3,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B97_3,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.269288d0*dfweight
      else if (dfname.eq."becke97-d") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B97_D,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_D,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B97_D,dfweight)
        nwxc_disp_id    = NWXC_DISP_D2
        nwxc_disp_s6    = 1.25d0
        nwxc_disp_s8    = 0.0d0
        nwxc_disp_sr6   = 0.0d0
        nwxc_disp_sr8   = 0.0d0
        nwxc_disp_alpha = 20.0d0
      else if (dfname.eq."becke97gga1") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B97_G,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B97_G,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B97_G,dfweight)
      else if (dfname.eq."becke98") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_B98,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B98,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B98,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.219847d0*dfweight
      else if (dfname.eq."bb1k") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_BB1K,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.58d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,0.58d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B95,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.42d0*dfweight
      else if (dfname.eq."bop") then ! alias
        nwxc_add_df = nwxc_all_add_func(NWXCP_BOP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_B88,dfweight)
      else if (dfname.eq."camb3lyp".or.dfname.eq."cam-b3lyp") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_CAMB3LYP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_CAMB88,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VWN5,0.19d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_LYP,0.81d0*dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.19d0
        nwxc_cam_beta   = 0.46d0
        nwxc_cam_gamma  = 0.33d0
        nwxc_cam_srhfex = .false.
      else if (dfname.eq."cams12g".or.dfname.eq."cam-s12g") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_CAMS12G,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_CAMS12G,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.0d0
        nwxc_cam_beta   = 0.34485046d0
        nwxc_cam_gamma  = 1.52420731d0
        nwxc_cam_srhfex = .false.
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.86124355d0
        nwxc_disp_sr6   = 1.20250451d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."cams12h".or.dfname.eq."cam-s12h") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_CAMS12H,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_CAMS12H,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_wght_hf    = nwxc_wght_hf + dfweight
        nwxc_cam_alpha  = 0.25d0
        nwxc_cam_beta   = 0.10897845d0
        nwxc_cam_gamma  = 0.48516891d0
        nwxc_cam_srhfex = .false.
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.37999939d0
        nwxc_disp_sr6   = 1.08034183d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."dldf") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_DLDF,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_DLDF,0.3855871d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_DLDF,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.6144129d0*dfweight
      else if (dfname.eq."ft97") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_FT97,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_FT97,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_FT97,dfweight)
      else if (dfname.eq."hcth".or.dfname.eq."hcth93") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_HCTH,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_HCTH,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_HCTH,dfweight)
      else if (dfname.eq."hcth120") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_HCTH120,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_HCTH120,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_HCTH120,dfweight)
      else if (dfname.eq."hcth147") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_HCTH147,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_HCTH147,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_HCTH147,dfweight)
      else if (dfname.eq."hcth407") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_HCTH407,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_HCTH407,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_HCTH407,dfweight)
      else if (dfname.eq."hcth407p") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_HCTH407P,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_HCTH407P,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_HCTH407P,dfweight)
c     else if (dfname.eq."hcth-ac") then
c       nwxc_add_df = nwxc_usr_add_func(NWXCP_HCTH_A,dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_X_HCTH_A,dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_C_HCTH_A,dfweight)
      else if (dfname.eq."hcthp14") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_HCTHP14,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_HCTHP14,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_HCTHP14,dfweight)
      else if (dfname.eq."kt1") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_KT1,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW91LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_KT1,dfweight)
      else if (dfname.eq."kt2") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_KT2,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,1.07173d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW91LDA,0.576727d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_KT1,dfweight)
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
      else if (dfname.eq."mpw1k") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_MPW1K,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW91LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_P91,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_MPW91,0.572d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.4280d0*dfweight
      else if (dfname.eq."mpw1b95") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_MPW1B95,dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_C_PW91LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B95,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_MPW91,0.69d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.31d0*dfweight
      else if (dfname.eq."mpwb1k") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_MPWB1K,dfweight)
c       nwxc_add_df = nwxc_add_df .and.
c    &                nwxc_add_func(NWXCP_C_PW91LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_B95,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_MPW91,0.56d0*dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.44d0*dfweight
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
      else if (dfname.eq."opt") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_OPT,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_OPT,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_OPT,dfweight)
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
      else if (dfname.eq."pbeop") then ! alias
        nwxc_add_df = nwxc_all_add_func(NWXCP_PBEOP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,dfweight)
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
      else if (dfname.eq."pw6b95") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PW6B95,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PW6B95,0.72d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PW6B95,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.28d0*dfweight
      else if (dfname.eq."pwb6k") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PWB6K,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,0.54d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PWB6K,0.54d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PWB6K,dfweight)
        nwxc_wght_hf = nwxc_wght_hf + 0.46d0*dfweight
      else if (dfname.eq."s12g") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_S12G,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,
     &                              dfweight*1.03842032d0)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_S12G,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.84432515d0
        nwxc_disp_sr6   = 1.17755954d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."s12h") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_S12H,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight*
     &                              0.75d0*1.02543951d0)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_S12H,dfweight*
     &                              0.75d0)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.25d0*dfweight
        nwxc_disp_id    = NWXC_DISP_D3
        nwxc_disp_s6    = 1.0d0
        nwxc_disp_s8    = 0.37705816d0
        nwxc_disp_sr6   = 1.07735222d0
        nwxc_disp_sr8   = 1.0d0
        nwxc_disp_alpha = 14.0d0
      else if (dfname.eq."sogga") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_SOGGA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SOGGA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."sogga11") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_SOGGA11,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SOGGA11,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_SOGGA11,dfweight)
      else if (dfname.eq."sogga11-x") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_SOGGA11_X,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SOGGA11_X,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_SOGGA11_X,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.4015d0*dfweight
      else if (dfname.eq."sop") then ! alias
        nwxc_add_df = nwxc_all_add_func(NWXCP_SOP,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
      else if (dfname.eq."ssb-d") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_SSB_D,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,1.079966d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SSB_D,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_SPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_KT1,0.358830d0*dfweight)
      else if (dfname.eq."tpss03".or.dfname.eq."tpss") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_TPSS,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_TPSS,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_TPSS,dfweight)
      else if (dfname.eq."xctpssh".or.dfname.eq."tpssh") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_TPSSH,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_TPSS,0.90d0*dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_TPSS,dfweight)
        nwxc_wght_hf= nwxc_wght_hf + 0.10d0*dfweight
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
c     Due to the Automatic Differentiation stuff nwxc_eval_df
c     now does a lot more behind the scenes. The complication is 
c     that the routine now wants to know whether we are using 
c     an LDA, GGA, or meta-GGA functional, which is what we want
c     to detect. 
c
c     The resolution of this problem is to assume, for now, that the
c     functional depends on everything, run nwxc_eval_df, and then
c     check what the functional really depends on.
c
      nwxc_depends_on_rho   = .true.
      nwxc_depends_on_gamma = .true.
      nwxc_depends_on_tau   = .true.
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
      nwxc_depends_on_gamma = nwxc_depends_on_gamma .or.
     +                        nwxc_depends_on_tau
      nwxc_depends_on_rho   = nwxc_depends_on_rho .or.
     +                        nwxc_depends_on_gamma
c
      end
C>
C> \brief Detect whether the functional supports 2nd derivatives
C>
C> For capabilities such as TDDFT and Hessians 2nd order partial
C> derivatives of the functional are required. This function detects
C> whether the current functional supports these derivatives. 
C>
C> \return
C> - `.TRUE.` if 2nd order derivatives are supported
C> - otherwise the function does not return (it will kill the code
C>   instead)
C>
      logical function nwxc_has_2nd()
      implicit none
#include "nwxcP.fh"
      integer i      !< Counter
      integer maxpt  !< The number of points
      parameter (maxpt = 3)
      integer maxpt2  !< The number of points
      parameter (maxpt2 = 2*maxpt)
      integer maxpt3  !< The number of points
      parameter (maxpt3 = 3*maxpt)
      integer maxpt12 !< The number of points
      parameter (maxpt12 = 12*maxpt)
      integer maxpt13 !< The number of points
      parameter (maxpt13 = 13*maxpt)
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
      double precision dfdr2(maxpt3)   !< 2nd derivative wrt rho
      double precision dfdg2(maxpt12)  !< 2nd derivative wrt rgamma
      double precision dfdt2(maxpt13)  !< 2nd derivative wrt tau
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
      data dfdr2/maxpt3*0.0d0/
      data dfdg2/maxpt12*0.0d0/
      data dfdt2/maxpt13*0.0d0/
c
      nwxc_has_2nd = .false.
c
      call nwxc_eval_df2(2,maxpt,rho,rgamma,tau,func,dfdr,dfdr2,
     &                   dfdg,dfdg2,dfdt,dfdt2)
c
      nwxc_has_2nd = .true.
c
      return
      end
C>
C> \brief Detect whether the functional supports 3rd derivatives
C>
C> For capabilities such as TDDFT gradients 3rd order partial
C> derivatives of the functional are required. This function detects
C> whether the current functional supports these derivatives. 
C>
C> \return
C> - `.TRUE.` if 3rd order derivatives are supported
C> - otherwise the function does not return (it will kill the code
C>   instead)
C>
      logical function nwxc_has_3nd()
      implicit none
#include "nwxcP.fh"
      integer i      !< Counter
      integer maxpt  !< The number of points
      parameter (maxpt = 3)
      integer maxpt2  !< The number of points
      parameter (maxpt2 = 2*maxpt)
      integer maxpt3  !< The number of points
      parameter (maxpt3 = 3*maxpt)
      integer maxpt4  !< The number of points
      parameter (maxpt4 = 4*maxpt)
      integer maxpt12 !< The number of points
      parameter (maxpt12 = 12*maxpt)
      integer maxpt13 !< The number of points
      parameter (maxpt13 = 13*maxpt)
      integer maxpt16 !< The number of points
      parameter (maxpt16 = 16*maxpt)
      integer maxpt31 !< The number of points
      parameter (maxpt31 = 31*maxpt)
      integer maxpt33 !< The number of points
      parameter (maxpt33 = 33*maxpt)
c     
      double precision rho(maxpt2)    !< The density
      double precision rgamma(maxpt3) !< The norm of the gradient
      double precision tau(maxpt2)    !< The kinetic energy density
c
      double precision func(maxpt)    !< The functional value
c
      double precision dfdr(maxpt2)   !< Derivative wrt rho
      double precision dfdg(maxpt3)   !< Derivative wrt rgamma
      double precision dfdt(maxpt2)   !< Derivative wrt tau
c
      double precision dfdr2(maxpt3)   !< 2nd derivative wrt rho
      double precision dfdg2(maxpt12)  !< 2nd derivative wrt rgamma
      double precision dfdt2(maxpt13)  !< 2nd derivative wrt tau
c
      double precision dfdr3(maxpt4)   !< 3rd derivative wrt rho
      double precision dfdg3(maxpt31)  !< 3rd derivative wrt rgamma
      double precision dfdt3(maxpt16)  !< 3rd derivative wrt rho tau
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
      data dfdr2/maxpt3*0.0d0/
      data dfdg2/maxpt12*0.0d0/
      data dfdt2/maxpt13*0.0d0/
      data dfdr3/maxpt4*0.0d0/
      data dfdg3/maxpt31*0.0d0/
      data dfdt3/maxpt16*0.0d0/
c
      nwxc_has_3nd = .false.
c
      call nwxc_eval_df3(2,maxpt,rho,rgamma,tau,func,dfdr,dfdr2,dfdr3,
     &                   dfdg,dfdg2,dfdg3,dfdt,dfdt2,dfdt3)
c
      nwxc_has_3nd = .true.
c
      return
      end

C> @}
c $Id$
