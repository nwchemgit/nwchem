C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_nwchem.F
C> Routines for the integration of the NWXC library into NWChem
C>
C> The library is designed to rely on as few NWChem specific things 
C> as possible. Nevertheless to make it work with NWChem we need
C> routines to handle 
C>
C> - the input, 
C>
C> - printing output,
C>
C> - storing the functional on the RTDB
C>
C> - retrieving the functional from the RTDB
C>
C> - deleting the functional from the RTDB
C>
C> This file contains these integration routines.
C>
C> @}
C
C> \ingroup nwpwxc_api
C> @{
C>
C> \brief Parse the input of the functional
C>
C> Read the keywords specifying the density functional and set the
C> appropriate items in the code.
C>
C> \return Return .true. if the functional was successfully parsed,
C> and .false. otherwise.
C>
      logical function nwpwxc_input()
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"
#include "nwpwxcP.fh"
      logical       stat    !< The overall status
      character*32  func    !< Functional tag
      character*32  local   !< "Locality" of the functional 
                            !< (local/nonlocal)
      double precision wght !< Functional weight
      integer ilen          !< Length of functional name
      logical olocal        !< If .true. use local part
      logical ononlocal     !< If .true. use non-local part
c
c     Functions
c
      logical  nwpwxc_add_df, nwpwxc_usr_add_func, nwpwxc_add_func
      external nwpwxc_add_df, nwpwxc_usr_add_func, nwpwxc_add_func
c
      call nwpwxc_init()
      nwpwxc_on = .true.
      stat = .true.
      do while (inp_a(func))
        olocal = .true.
        ononlocal = .true.
c
        if (inp_a(local)) then
          if (inp_compare(.false.,'local',local)) then
            ononlocal = .false.
          else if (inp_compare(.false.,'nonlocal',local)) then
            olocal = .false.
          else
            call inp_prev_field()
          endif
        endif
        if (.not.inp_f(wght)) then
          wght = 1.0d0
c         call inp_prev_field() ! not necessary inp_f does not move
c                               ! cursor if failed
        endif
        ilen = inp_strlen(func)
        select case (func(1:ilen))
          case ("xpbe96")
            if (olocal.and.ononlocal) then
              stat = stat .and. nwpwxc_add_df(func(1:ilen),wght)
            else if ((.not.olocal).and.ononlocal) then
              stat = stat .and. nwpwxc_usr_add_func(NWXCP_X_PBE,wght)
              stat = stat.and.nwpwxc_usr_add_func(NWXCP_X_SLATER,-wght)
              stat = stat .and. nwpwxc_add_func(NWXCP_X_PBE,wght)
            else if (olocal.and.(.not.ononlocal)) then
              write(LuOut,*)' xpbe96 "local" is equivalent to the',
     &                ' slater exchange local functional. '
              write(LuOut,*)' Resubmit using slater exchange. '
              call errquit(' xc_inp - functional description ',0,
     &                     INPUT_ERR)
            else
              write(LuOut,*)' Want to add nothing',olocal,ononlocal
              call errquit(' xc_inp - functional description ',0,
     &                     INPUT_ERR)
            endif
          case default
            stat = stat .and. nwpwxc_add_df(func(1:ilen),wght)
        end select
      enddo
      call nwpwxc_detect_dependencies()
      nwpwxc_input = stat
      return
      end
C>
C> \brief Print the current density functional
C>
C> Extract various components of the current functional and write them
C> to the NWChem output file. The functional is printed as specified
C> on the user list, so as to match the user input as closely as 
C> possible.
C>
      subroutine nwpwxc_print()
      implicit none
#include "stdio.fh"
#include "inp.fh"
#include "global.fh"
#include "errquit.fh"
#include "nwpwxc_param.fh"
#include "nwpwxcP.fh"
c
      double precision wght  !< The functional weight
      character*255 term     !< The functional term details
      integer iterm          !< The functional identifier
      integer ientry         !< Counter
      integer dispid         !< Dipersion type ID
      double precision s6    !< The C6 coefficient
      double precision s8    !< The C8 coefficient
      double precision sr6   !< The Cr,6 coefficient
      double precision sr8   !< The Cr,8 coefficient
      double precision alpha !< The damping parameter
c
c     Functions
c
      logical  nwpwxc_has_hfx, nwpwxc_has_mp2c, 
     >         nwpwxc_has_cam, nwpwxc_has_disp
      external nwpwxc_has_hfx, nwpwxc_has_mp2c,
     >         nwpwxc_has_cam, nwpwxc_has_disp
c
      logical  nwpwxc_get_df, nwpwxc_get_info
      external nwpwxc_get_df, nwpwxc_get_info
c
      double precision nwpwxc_wght_mp2c, nwpwxc_wght_hfx
      external         nwpwxc_wght_mp2c, nwpwxc_wght_hfx

c
      nwpwxc_oroot = ga_nodeid().eq.0
      if (nwpwxc_oroot) then
        write(luout,*)" =================================="
        write(luout,*)" === Current Density Functional ==="
        write(luout,*)" =================================="
        write(luout,*)
        if (nwpwxc_has_hfx()) then
          wght = nwpwxc_wght_hfx()
          if (.not.nwpwxc_get_info(NWXCP_HFXC,term)) then
            call errquit("illegal functional",NWXCP_HFXC,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        endif
        if (nwpwxc_has_mp2c()) then
          wght = nwpwxc_wght_mp2c()
          if (.not.nwpwxc_get_info(NWXCP_MP2C,term)) then
            call errquit("illegal functional",NWXCP_MP2C,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        endif
        ientry = 1
        do while (nwpwxc_get_df(ientry,iterm,wght))
          if (.not.nwpwxc_get_info(iterm,term)) then
            call errquit("illegal functional",iterm,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        enddo
c       if (nwpwxc_has_cam()) then
c         if (ga_nodeid().eq.0) then
c           write(luout,*)"must think about CAM"
c         endif
c       endif
        if (nwpwxc_has_disp()) then
          if (ga_nodeid().eq.0) then
            call nwpwxc_get_disp(dispid,s6,s8,sr6,sr8,alpha)
            select case (dispid)
              case (NWXC_DISP_D1)
                write(luout,1010)"GRIMME D1","Correction type"
                write(luout,1000)s6,"C6 coefficient"
                write(luout,1000)alpha,"Alpha damping parameter"
              case (NWXC_DISP_D2)
                write(luout,1010)"GRIMME D2","Correction type"
                write(luout,1000)s6,"C6 coefficient"
                write(luout,1000)alpha,"Alpha damping parameter"
              case (NWXC_DISP_D3)
                write(luout,1010)"GRIMME D3","Correction type"
                write(luout,1000)s6,"C6 coefficient"
                write(luout,1000)s8,"C8 coefficient"
                write(luout,1000)sr6,"Cr,6 coefficient"
                write(luout,1000)sr8,"Cr,8 coefficient"
                write(luout,1000)alpha,"Alpha damping parameter"
              case default
                call errquit("illegal dispersion correction",
     +                       dispid,UERR)
            end select
          endif
        endif
      endif
cDEBUG
c     call nwpwxc_printP()
cDEBUG
c
 1000 format(1x,f14.8,"  ",a)
 1010 format(1x,a14,"  ",a)
      end
C>
C> \brief Store the density functional on the RTDB
C>
C> Store the current settings of the density functional on the RTDB
C> under a specific name. The companion routine `nwpwxc_rtdb_load`
C> may be used to retrieve the data from the RTDB.
C> We integrate the variable names in the RTDB keys to eliminate any
C> source of ambiguity for what the data elements are.
C>
C> This routine is for convenience reasons implemented as what would
C> be a "friend" in C++. Otherwise a zoo of additional routines would 
C> need to query and set all sorts of variables.
C>
C> \return Return .true. if the data was successfully stored, and
C> .false. otherwise.
C> 
      logical function nwpwxc_rtdb_store(rtdb,name)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "nwpwxcP.fh"
#include "inp.fh"
      integer rtdb       !< [Input] The RTDB handle
      character*(*) name !< [Input] The name under which the current
                         !< functional will be stored
c
      character*256 rtdb_key !< Key to access data on the RTDB
      integer lenname        !< Length of the functional name
      logical status         !< Status of the store operation
      logical  nwpwxc_has_cam
      external nwpwxc_has_cam
cHACK
      logical junk
cHACK
c
      lenname = inp_strlen(name)
      status = .true.
      rtdb_key = "dft:xc:"//name(1:lenname)//":nwpwxc_num_f"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_INT,1,nwpwxc_num_f)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_usr_num"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_INT,1,nwpwxc_usr_num)
c
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_alpha"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,
     >                           nwpwxc_cam_alpha)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_beta"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,nwpwxc_cam_beta)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_gamma"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,
     >                           nwpwxc_cam_gamma)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_srhfex"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_LOG,1,
     >                           nwpwxc_cam_srhfex)
c
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_id"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_INT,1,nwpwxc_disp_id)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_s6"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_s6)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_s8"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_s8)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_sr6"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_sr6)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_sr8"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_sr8)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_alpha"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,
     >                           nwpwxc_disp_alpha)
c
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_wght_hf"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,nwpwxc_wght_hf)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_wght_mp2"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,nwpwxc_wght_mp2)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_usr_wghts"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,
     +                           nwpwxc_usr_num,nwpwxc_usr_wghts)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_wghts"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,
     +                           nwpwxc_num_f,nwpwxc_wghts)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_ids"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_INT,
     +                           nwpwxc_num_f,nwpwxc_ids)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_usr_ids"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_INT,
     +                           nwpwxc_usr_num,nwpwxc_usr_ids)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_depends_on_rho"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_LOG,1,
     +                           nwpwxc_depends_on_rho)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_depends_on_gamma"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_LOG,1,
     +                           nwpwxc_depends_on_gamma)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_depends_on_tau"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_LOG,1,
     +                           nwpwxc_depends_on_tau)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_rho_tol"
      status=status.and.rtdb_put(rtdb,rtdb_key,MT_DBL,1,
     +                           nwpwxc_rho_tol)
cHACK
c     This stuff is to be removed when the integration is finalized
      if (nwpwxc_has_cam()) then
        status=status.and.
     &     rtdb_put(rtdb,'int:cando_txs',mt_log,1,.false.)
        status=status.and.
     &     rtdb_put(rtdb,'dft:cam_exch',mt_log,1,.true.)
        status=status.and.
     &     rtdb_put(rtdb,'dft:direct',mt_log,1,.true.)
      else
        status=status.and.
     &     rtdb_put(rtdb,'dft:cam_exch',mt_log,1,.false.)
        if (rtdb_get(rtdb,'int:cando_txs',mt_log,1,junk)) then
          status=status.and.
     &       rtdb_delete(rtdb,'int:cando_txs')
        endif
      endif
cHACK
      nwpwxc_rtdb_store = status
c
      return
      end
C>
C> \brief Load the density functional from the RTDB
C>
C> Load the current settings of the density functional stored with
C> a particular name from the RTDB. The companion routine
C> `nwpwxc_rtdb_store` may be used to store the data on the RTDB.
C> We integrate the variable names in the RTDB keys to eliminate any
C> source of ambiguity for what the data elements are.
C>
C> This routine is for convenience reasons implemented as what would
C> be a "friend" in C++. Otherwise a zoo of additional routines would 
C> need to query and set all sorts of variables.
C>
C> \return Return .true. if the data was successfully loaded, and
C> .false. if no functional with the specified name was found, abort 
C> if only partial information on a particular functional is found.
C> 
      logical function nwpwxc_rtdb_load(rtdb,name)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "nwpwxcP.fh"
      integer rtdb       !< [Input] The RTDB handle
      character*(*) name !< [Input] The name under which the functional
                         !< is stored
c
      character*256 rtdb_key !< Key to access data on the RTDB
      integer lenname        !< Length of the functional name
      logical status         !< Status of the store operation
c
      call nwpwxc_init()
      lenname = inp_strlen(name)
      status = .true.
      rtdb_key = "dft:xc:"//name(1:lenname)//":nwpwxc_num_f"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_INT,1,nwpwxc_num_f)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_usr_num"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_INT,1,nwpwxc_usr_num)
c
      nwpwxc_rtdb_load = status
      nwpwxc_on = status
      if (.not. status) then
c
c       The specified functional was not found so return .false.
c
        return
      endif
c
      if (nwpwxc_num_f.gt.NWXC_MAX_F.or.
     >   nwpwxc_usr_num.gt.NWXC_MAX_F) then
c
c       The RTDB was created by a code with different dimensions. Hence
c       the RTDB is incompatible with this code. So, produce an error
c       message and bomb out before something really bad happens...
c
        call errquit("RTDB is incompatible with this version",0,UERR)
      endif
c
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_alpha"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,
     >                           nwpwxc_cam_alpha)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_beta"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,nwpwxc_cam_beta)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_gamma"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,
     >                           nwpwxc_cam_gamma)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_cam_srhfex"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_LOG,1,
     >                          nwpwxc_cam_srhfex)
c
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_id"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_INT,1,nwpwxc_disp_id)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_s6"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_s6)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_s8"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_s8)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_sr6"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_sr6)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_sr8"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,nwpwxc_disp_sr8)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_disp_alpha"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,
     >                           nwpwxc_disp_alpha)
c
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_wght_hf"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,nwpwxc_wght_hf)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_wght_mp2"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,nwpwxc_wght_mp2)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_usr_wghts"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,
     +                           nwpwxc_usr_num,nwpwxc_usr_wghts)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_wghts"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,
     +                           nwpwxc_num_f,nwpwxc_wghts)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_ids"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_INT,
     +                           nwpwxc_num_f,nwpwxc_ids)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_usr_ids"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_INT,
     +                           nwpwxc_usr_num,nwpwxc_usr_ids)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_depends_on_rho"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_LOG,1,
     +                           nwpwxc_depends_on_rho)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_depends_on_gamma"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_LOG,1,
     +                           nwpwxc_depends_on_gamma)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_depends_on_tau"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_LOG,1,
     +                           nwpwxc_depends_on_tau)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_rho_tol"
      status=status.and.rtdb_get(rtdb,rtdb_key,MT_DBL,1,
     +                           nwpwxc_rho_tol)
c
      nwpwxc_rtdb_load = status
      if (.not.status) then
c
c       We found nwpwxc_num_f and nwpwxc_usr_num so the functional must have
c       existed but some of the other components were not found. So
c       the functional data on the RTDB is corrupt and we would better
c       terminate the program here.
c
        call errquit("Functional data on the RTDB corrupt",0,UERR)
      endif
c
      return
      end
C>
C> \brief Delete the density functional from the RTDB
C>
C> Delete the current settings of the density functional stored with
C> a particular name from the RTDB. We integrate the variable names
C> in the RTDB keys to eliminate any source of ambiguity for what
C> the data elements are. If the functional could not be deleted 
C> the routine returns .false.. This allows one to attempt to
C> remove a functional even if no functional was previously defined.
C> Such behavior is OK depending on the coding conventions in the
C> host code.
C>
C> \return Return .true. if the functional was successfully
C> deleted, and .false. otherwise.
C> 
      logical function nwpwxc_rtdb_delete(rtdb,name)
      implicit none
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb       !< [Input] The RTDB handle
      character*(*) name !< [Input] The name under which the functional
                         !< is stored

c
      character*256 rtdb_key !< Key to access data on the RTDB
      integer lenname        !< Length of the functional name
      logical status         !< Status of the store operation
c
      lenname = inp_strlen(name)
      status = .true.
      rtdb_key = "dft:xc:"//name(1:lenname)//":nwpwxc_num_f"
      status=status.and.rtdb_delete(rtdb,rtdb_key)
      rtdb_key="dft:xc:"//name(1:lenname)//":nwpwxc_usr_num"
      status=status.and.rtdb_delete(rtdb,rtdb_key)

      nwpwxc_rtdb_delete = status
      return
      end
C> @}
C
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Private printing routine to investigate the state of the
C> evaluation list
C>
      subroutine nwpwxc_printP()
      implicit none
#include "stdio.fh"
#include "inp.fh"
#include "global.fh"
#include "errquit.fh"
#include "nwpwxcP.fh"
c
      double precision wght !< The functional weight
      character*255 term    !< The functional term details
      integer iterm         !< The functional identifier
      integer ientry        !< Counter
c
c     Functions
c
      logical  nwpwxc_has_hfx, nwpwxc_has_mp2c, 
     >         nwpwxc_has_cam, nwpwxc_has_disp
      external nwpwxc_has_hfx, nwpwxc_has_mp2c, 
     >         nwpwxc_has_cam, nwpwxc_has_disp
c
      logical  nwpwxc_get_df, nwpwxc_get_info
      external nwpwxc_get_df, nwpwxc_get_info
c
      double precision nwpwxc_wght_mp2c, nwpwxc_wght_hfx
      external         nwpwxc_wght_mp2c, nwpwxc_wght_hfx

c
 20   format(1x,a20,i4)
      nwpwxc_oroot = ga_nodeid().eq.0
      if (nwpwxc_oroot) then
        call util_flush(luout)
#ifdef dlflfglf
        write(luout,*)" =================================="
        write(luout,*)" === Current Density Functional ==="
        write(luout,*)" =================================="
        write(luout,*)
        write(luout,20)"NWXCP_HFXC",NWXCP_HFXC
        write(luout,20)"NWXCP_MP2C",NWXCP_MP2C
        write(luout,20)"NWXCP_X_SLATER",NWXCP_X_SLATER
        write(luout,20)"NWXCP_X_B88",NWXCP_X_B88
        write(luout,20)"NWXCP_X_B0",NWXCP_X_B0
        write(luout,20)"NWXCP_X_BNL",NWXCP_X_BNL
        write(luout,20)"NWXCP_X_CAMB88",NWXCP_X_CAMB88
        write(luout,20)"NWXCP_X_CAMLSD",NWXCP_X_CAMLSD
        write(luout,20)"NWXCP_X_CAMPBE",NWXCP_X_CAMPBE
        write(luout,20)"NWXCP_X_DLDF",NWXCP_X_DLDF
        write(luout,20)"NWXCP_X_M05",NWXCP_X_M05
        write(luout,20)"NWXCP_X_M05_2X",NWXCP_X_M05_2X
        write(luout,20)"NWXCP_X_M06",NWXCP_X_M06
        write(luout,20)"NWXCP_X_M06_HF",NWXCP_X_M06_HF
        write(luout,20)"NWXCP_X_M06_L",NWXCP_X_M06_L
        write(luout,20)"NWXCP_X_M06_2X",NWXCP_X_M06_2X
        write(luout,20)"NWXCP_X_M08_HX",NWXCP_X_M08_HX
        write(luout,20)"NWXCP_X_M08_SO",NWXCP_X_M08_SO
        write(luout,20)"NWXCP_X_M11",NWXCP_X_M11
        write(luout,20)"NWXCP_X_M11_L",NWXCP_X_M11_L
        write(luout,20)"NWXCP_X_PBE",NWXCP_X_PBE
        write(luout,20)"NWXCP_X_PBE0",NWXCP_X_PBE0
        write(luout,20)"NWXCP_X_PBESOL",NWXCP_X_PBESOL
        write(luout,20)"NWXCP_X_REVPBE",NWXCP_X_REVPBE
        write(luout,20)"NWXCP_X_RPBE",NWXCP_X_RPBE
        write(luout,20)"NWXCP_X_TPSS",NWXCP_X_TPSS
        write(luout,20)"NWXCP_X_WPBE",NWXCP_X_WPBE
        write(luout,20)"NWXCP_X_VS98",NWXCP_X_VS98
        write(luout,20)"NWXCP_X_SCAN",NWXCP_X_SCAN
        write(luout,20)"NWXCP_C_B95",NWXCP_C_B95
        write(luout,20)"NWXCP_C_DLDF",NWXCP_C_DLDF
        write(luout,20)"NWXCP_C_LYP",NWXCP_C_LYP
        write(luout,20)"NWXCP_C_M05",NWXCP_C_M05
        write(luout,20)"NWXCP_C_M05_2X",NWXCP_C_M05_2X
        write(luout,20)"NWXCP_C_M06",NWXCP_C_M06
        write(luout,20)"NWXCP_C_M06_HF",NWXCP_C_M06_HF
        write(luout,20)"NWXCP_C_M06_L",NWXCP_C_M06_L
        write(luout,20)"NWXCP_C_M06_2X",NWXCP_C_M06_2X
        write(luout,20)"NWXCP_C_M08_HX",NWXCP_C_M08_HX
        write(luout,20)"NWXCP_C_M08_SO",NWXCP_C_M08_SO
        write(luout,20)"NWXCP_C_M11",NWXCP_C_M11
        write(luout,20)"NWXCP_C_M11_L",NWXCP_C_M11_L
        write(luout,20)"NWXCP_C_PW6B95",NWXCP_C_PW6B95
        write(luout,20)"NWXCP_C_PWB6K",NWXCP_C_PWB6K
        write(luout,20)"NWXCP_C_PW91LDA",NWXCP_C_PW91LDA
        write(luout,20)"NWXCP_C_PBE",NWXCP_C_PBE
        write(luout,20)"NWXCP_C_PBESOL",NWXCP_C_PBESOL
        write(luout,20)"NWXCP_C_TPSS",NWXCP_C_TPSS
        write(luout,20)"NWXCP_C_VS98",NWXCP_C_VS98
        write(luout,20)"NWXCP_C_SCAN",NWXCP_C_SCAN
        write(luout,20)"NWXCP_C_VWN1",NWXCP_C_VWN1
        write(luout,20)"NWXCP_C_VWN1_RPA",NWXCP_C_VWN1_RPA
        write(luout,20)"NWXCP_C_VWN2",NWXCP_C_VWN2
        write(luout,20)"NWXCP_C_VWN3",NWXCP_C_VWN3
        write(luout,20)"NWXCP_C_VWN4",NWXCP_C_VWN4
        write(luout,20)"NWXCP_C_VWN5",NWXCP_C_VWN5
        write(luout,20)"NWXCP_B3LYP",NWXCP_B3LYP
        write(luout,20)"NWXCP_B3PW91",NWXCP_B3PW91
        write(luout,20)"NWXCP_B0LYP",NWXCP_B0LYP
        write(luout,20)"NWXCP_BLYP",NWXCP_BLYP
        write(luout,20)"NWXCP_DLDF",NWXCP_DLDF
        write(luout,20)"NWXCP_LDA",NWXCP_LDA
        write(luout,20)"NWXCP_M05",NWXCP_M05
        write(luout,20)"NWXCP_M05_2X",NWXCP_M05_2X
        write(luout,20)"NWXCP_M06",NWXCP_M06
        write(luout,20)"NWXCP_M06_HF",NWXCP_M06_HF
        write(luout,20)"NWXCP_M06_L",NWXCP_M06_L
        write(luout,20)"NWXCP_M06_2X",NWXCP_M06_2X
        write(luout,20)"NWXCP_M08_HX",NWXCP_M08_HX
        write(luout,20)"NWXCP_M08_SO",NWXCP_M08_SO
        write(luout,20)"NWXCP_M11",NWXCP_M11
        write(luout,20)"NWXCP_M11_L",NWXCP_M11_L
        write(luout,20)"NWXCP_PBE",NWXCP_PBE
        write(luout,20)"NWXCP_PBE0",NWXCP_PBE0
        write(luout,20)"NWXCP_PBESOL",NWXCP_PBESOL
        write(luout,20)"NWXCP_REVPBE",NWXCP_REVPBE
        write(luout,20)"NWXCP_REVPBE0",NWXCP_REVPBE0
        write(luout,20)"NWXCP_RPBE",NWXCP_RPBE
        write(luout,20)"NWXCP_TPSS",NWXCP_TPSS
        write(luout,20)"NWXCP_HSE",NWXCP_HSE
        write(luout,20)"NWXCP_VS98",NWXCP_VS98
        write(luout,20)"NWXCP_SCAN",NWXCP_SCAN
        write(luout,20)"NWXCP_PBE_D1",NWXCP_PBE_D1
#endif
        write(luout,*)
        if (nwpwxc_has_hfx()) then
          wght = nwpwxc_wght_hfx()
          if (.not.nwpwxc_get_info(NWXCP_HFXC,term)) then
            call errquit("illegal functional",NWXCP_HFXC,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        endif
        if (nwpwxc_has_mp2c()) then
          wght = nwpwxc_wght_mp2c()
          if (.not.nwpwxc_get_info(NWXCP_MP2C,term)) then
            call errquit("illegal functional",NWXCP_MP2C,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        endif
        ientry = 1
        do ientry = 1, nwpwxc_num_f
          if (.not.nwpwxc_get_info(nwpwxc_ids(ientry),term)) then
            call errquit("illegal functional",NWXCP_HFXC,UERR)
          endif
          write(luout,1010)nwpwxc_wghts(ientry),nwpwxc_ids(ientry),
     +                     term(1:inp_strlen(term))
        enddo
        if (nwpwxc_has_cam()) then
c         call errquit("must think about CAM",0,UERR)
          write(luout,1000)nwpwxc_cam_alpha,"NWXC_CAM_ALPHA"
          write(luout,1000)nwpwxc_cam_beta ,"NWXC_CAM_BETA"
          write(luout,1000)nwpwxc_cam_gamma,"NWXC_CAM_GAMMA"
        endif
        if (nwpwxc_has_disp()) then
c         call errquit("must think about DISP",0,UERR)
        endif
        call util_flush(luout)
      endif
c
 1000 format(1x,f12.6,"  ",a)
 1010 format(1x,f12.6,"  ",i4,"  ",a)
      end
C>
C> \brief Print the current density functional in nwpw format
C>
C> Extract various components of the current functional and write them
C> to the NWChem output file. The functional is printed as specified
C> on the user list, so as to match the user input as closely as 
C> possible.
C>
      subroutine nwpwxc_print_nwpw()
      implicit none
#include "stdio.fh"
#include "inp.fh"
#include "global.fh"
#include "errquit.fh"
#include "nwpwxcP.fh"
c
      double precision wght !< The functional weight
      character*255 term    !< The functional term details
      integer iterm         !< The functional identifier
      integer ientry        !< Counter
c
c     Functions
c
      logical  nwpwxc_has_hfx, nwpwxc_has_mp2c, 
     >         nwpwxc_has_cam, nwpwxc_has_disp
      external nwpwxc_has_hfx, nwpwxc_has_mp2c, 
     >         nwpwxc_has_cam, nwpwxc_has_disp
c
      logical  nwpwxc_get_df, nwpwxc_get_info
      external nwpwxc_get_df, nwpwxc_get_info
c
      double precision nwpwxc_wght_mp2c, nwpwxc_wght_hfx
      external         nwpwxc_wght_mp2c, nwpwxc_wght_hfx

c
      nwpwxc_oroot = ga_nodeid().eq.0
      if (nwpwxc_oroot) then
        write(luout,*)
        write(luout,*)"     exchange correlation = "
        if (nwpwxc_has_hfx()) then
          wght = nwpwxc_wght_hfx()
          if (.not.nwpwxc_get_info(NWXCP_HFXC,term)) then
            call errquit("illegal functional",NWXCP_HFXC,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        endif
        if (nwpwxc_has_mp2c()) then
          wght = nwpwxc_wght_mp2c()
          if (.not.nwpwxc_get_info(NWXCP_MP2C,term)) then
            call errquit("illegal functional",NWXCP_MP2C,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        endif
        ientry = 1
        do while (nwpwxc_get_df(ientry,iterm,wght))
          if (.not.nwpwxc_get_info(iterm,term)) then
            call errquit("illegal functional",iterm,UERR)
          endif
          write(luout,1000)wght,term(1:inp_strlen(term))
        enddo
c       if (nwpwxc_has_cam()) then
c         if (ga_nodeid().eq.0) then
c           write(luout,*)"must think about CAM"
c         endif
c       endif
        if (nwpwxc_has_disp()) then
          if (ga_nodeid().eq.0) then
            write(luout,*)"must think about DISP"
          endif
          call errquit("must think about DISP",0,UERR)
        endif
        write(luout,*)
      endif
c
 1000 format(2x,f14.8,"  ",a)
      end



C> @}
c $Id$
