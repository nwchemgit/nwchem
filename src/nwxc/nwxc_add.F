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
      logical  nwxc_set_dependencies
      external nwxc_set_dependencies
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
        nwxc_usr_add_func = nwxc_set_dependencies(iterm)
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
      else if (dfname.eq."xpbe96") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,dfweight)
      else if (dfname.eq."xrevpbe") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_REVPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_REVPBE,dfweight)
      else if (dfname.eq."xrpbe") then
        nwxc_add_df = nwxc_usr_add_func(NWXCP_X_RPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_RPBE,dfweight)
      else if (dfname.eq."xvs98") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_X_VS98,dfweight)
      else if (dfname.eq."cpw91lda".or.dfname.eq."pw91lda") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PW91LDA,dfweight)
      else if (dfname.eq."cpbe96") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."cvs98") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VS98,dfweight)
      else if (dfname.eq."cvwn5") then
        nwxc_add_df = nwxc_all_add_func(NWXCP_C_VWN5,dfweight)
      else if (dfname.eq."lda") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_LDA,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_VWN5,dfweight)
      else if (dfname.eq."pbe96") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_PBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."revpbe") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_REVPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_REVPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
      else if (dfname.eq."rpbe") then ! alias
        nwxc_add_df = nwxc_usr_add_func(NWXCP_RPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_SLATER,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_X_RPBE,dfweight)
        nwxc_add_df = nwxc_add_df .and.
     &                nwxc_add_func(NWXCP_C_PBE,dfweight)
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
C> \brief Set the dependencies for the new functional term
C>
C> An important aspect in DFT is to keep track of what quantities
C> a functional depends on as this relates to many expressions a code
C> needs to evaluate. This function updates the appropriate data
C> in accordance to the term being added.
C>
      logical function nwxc_set_dependencies(iterm)
      implicit none
#include "nwxcP.fh"
      integer iterm
      logical stat
      stat = .true.
      select case (iterm)
        case (NWXCP_HFXC,NWXCP_MP2C)
c         do not do anything
        case (NWXCP_X_SLATER)
          nwxc_depends_on_rho   = .true.
        case (NWXCP_X_PBE)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
        case (NWXCP_X_REVPBE)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
        case (NWXCP_X_RPBE)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
        case (NWXCP_X_VS98)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
          nwxc_depends_on_tau   = .true.
        case (NWXCP_C_PW91LDA)
          nwxc_depends_on_rho   = .true.
        case (NWXCP_C_PBE)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
        case (NWXCP_C_VS98)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
          nwxc_depends_on_tau   = .true.
        case (NWXCP_C_VWN5)
          nwxc_depends_on_rho   = .true.
        case (NWXCP_LDA)
          nwxc_depends_on_rho   = .true.
        case (NWXCP_PBE)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
        case (NWXCP_REVPBE)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
        case (NWXCP_RPBE)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
        case (NWXCP_VS98)
          nwxc_depends_on_rho   = .true.
          nwxc_depends_on_gamma = .true.
          nwxc_depends_on_tau   = .true.
        case default
          stat = .false.
      end select
      nwxc_set_dependencies = stat
      return
      end
C> @}
