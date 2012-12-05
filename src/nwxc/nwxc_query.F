C> \ingroup nwxc
C> @{
C>
C> \file nwxc_query.F
C> Routines to query the properties of the current functional
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
C> correction term (i.e. is a hybrid functional), and .false. otherwise.
C>
      logical function nwxc_has_disp()
      implicit none
#include "nwxcP.fh"
      call errquit("not sure how to do this yet",0,0)
      nwxc_has_disp = nwxc_wght_hf .gt. 0.0d0
      return
      end
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
         case (NWXC_HFXC)
           cdata = "Hartree-Fock Exchange"
         case (NWXC_MP2C)
           cdata = "Moller-Plesset 2nd Order Correlation"
         case (NWXC_SLATER)
           cdata = "Slater Exchange ()"
         case default
           stat = .false.
       end select
       i = lendat
       do while (cdata(i:i).ne." ")
         i = i - 1
       enddo
       stat = stat .and. (i.lt.lendat)
       nwxc_get_info = stat
c
       return
       end
C> @}
