C> \ingroup nwxc
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
      nwxc_usr_add_func = .false.
c
      i = 1 
      do while (i.le.nwxc_usr_num.and.nwxc_usr_ids(i).ne.iterm)
        i = i + 1
      enddo
      if (i.gt.nwxc_usr_num) then
c       the functional was not found on the list
        if (nwxc_usr_num.ge.NWXC_MAX_F) then
c         we are out of space
          return
        endif
        nwxc_usr_num = nwxc_usr_num + 1
        nwxc_usr_ids(nwxc_usr_num)   = iterm
        nwxc_usr_wghts(nwxc_usr_num) = weight
      else
c       the functional was found on the list
        nwxc_usr_wghts(i) = nwxc_usr_wghts(i) + weight
      endif
      nwxc_usr_add_func = .true.
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
      character*(*)    dfname   !< [Input] The name of the functional
      double precision dfweight !< [Input] The weight of the functional
c
      logical  nwxc_add_func, nwxc_usr_add_func, nwxc_all_add_func
      external nwxc_add_func, nwxc_usr_add_func, nwxc_all_add_func
c
      nwxc_add_df = .false.
c
      if (name.eq."hfex") then
        nwxc_wght_hf = nwxc_wght_hf + weight
        nwxc_add_df = nwxc_usr_add_func(NWXC_HFXC,dfweight)
      else if (name.eq."mp2") then
        nwxc_wght_mp2 = nwxc_wght_mp2 + weight
        nwxc_add_df = nwxc_usr_add_func(NWXC_MP2C,dfweight)
      else if (name.eq."slater") then
        nwxc_add_df = nwxc_all_add_func(NWXC_SLATER,dfweight)
      else
c       Oh dear, the functional was not recognized
      endif

      return
      end
C> @}
