      logical function nwchem_libxc_init()
      use,intrinsic :: iso_c_binding
#include "cdft.fh"
#include "libxc.fh"
      libxcon = .false.
      libxc_nfuncs = 0
      nwchem_libxc_init = .true.
      end function

      logical function is_libxcon()
      use,intrinsic :: iso_c_binding
#include "cdft.fh"
      is_libxcon = libxcon
      end function


      logical function nwchem_libxc_family(family)
      use,intrinsic :: iso_c_binding
#ifdef USE_LIBXC
      use xc_f03_lib_m
#endif

      implicit none

      character*(*) family
      integer ifunc

#include "libxc.fh"

      nwchem_libxc_family = .false.

#ifdef USE_LIBXC
      do ifunc=1,libxc_nfuncs
        if (family.eq."gga") then
          select case(libxc_family(ifunc))
          case(xc_family_gga,xc_family_hyb_gga,xc_family_mgga,
     $         xc_family_hyb_mgga)
            nwchem_libxc_family = .true.
          case default
            nwchem_libxc_family = .false.
          end select            
        else if (familY.eq."mgga") then
          select case(libxc_family(ifunc))
          case(xc_family_mgga,xc_family_hyb_mgga)
            nwchem_libxc_family = .true.
          case default
            nwchem_libxc_family = .false.
          end select
        elseif (family.eq."lapl") then
          nwchem_libxc_family = iand(libxc_flags(ifunc),
     $                           xc_flags_needs_laplacian).eq.
     $                           xc_flags_needs_laplacian
        elseif (family.eq."exc") then
          nwchem_libxc_family = iand(libxc_flags(ifunc),
     $                           xc_flags_have_exc).eq.
     $                           xc_flags_have_exc
          if (.not.nwchem_libxc_family) return
        elseif (family.eq."fxc") then
          nwchem_libxc_family = iand(libxc_flags(ifunc),
     $                           xc_flags_have_fxc).eq.
     $                           xc_flags_have_fxc
          if (.not.nwchem_libxc_family) return
        elseif (family.eq."kxc") then
          nwchem_libxc_family = iand(libxc_flags(ifunc),
     $                           xc_flags_have_kxc).eq.
     $                           xc_flags_have_kxc
          if (.not.nwchem_libxc_family) return
        endif
      enddo
#endif
      end function

      subroutine nwchem_libxc_print
      use,intrinsic :: iso_c_binding
#ifdef USE_LIBXC
      use xc_f03_lib_m, only: xc_f03_functional_get_name
#endif
      implicit none

#include "libxc.fh"
#include "stdio.fh"
#include "global.fh"

      character(len=1024) :: namexc,version
      character(len=5) :: blanks
      integer :: ifunc
      double precision :: fact


#ifdef USE_LIBXC
      if (ga_nodeid().eq.0) then

        do ifunc=1,libxc_nfuncs
          namexc = xc_f03_functional_get_name(libxc_funcs(ifunc))
          fact = libxc_facts(ifunc)
          write(blanks,'(I5)') 50-len_trim(namexc)
          write(luOut,'(a'//adjustl(blanks)//')',advance='no') 
          write(luOut,100) trim(namexc),fact
        enddo

      endif

  100 format(A50,1x,f6.3)
#endif
      end subroutine

      subroutine nwchem_libxc_print_header
      use,intrinsic :: iso_c_binding
#ifdef USE_LIBXC
      use xc_f03_lib_m, only: xc_f03_version_string
#endif
      implicit none

#include "libxc.fh"
#include "stdio.fh"
#include "global.fh"

      character(len=1024) :: namexc,version
      character(len=5) :: blanks
      integer :: ifunc
      double precision :: fact


#ifdef USE_LIBXC
      if (ga_nodeid().eq.0) then

        call xc_f03_version_string(version)
        write(luOut,100) trim(adjustl(version))

      endif

  100 format(10X,'Using LibXC version ',A10)
#endif
      end subroutine

      subroutine nwchem_libxc_rtdbput(rtdb,modname)
      use, intrinsic :: iso_c_binding
      implicit none
#include "libxc.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      character(len=*) modname
      integer :: itemp(maxfunc)
      integer,parameter :: one = 1

      if (.not.rtdb_put(rtdb,modname//':libxc_nfuncs',mt_int,one,
     $    libxc_nfuncs))
     $  call errquit(' libxc_rtdbput: rtdb_put failed',1,RTDB_ERR)

      itemp(:) = libxc_funcs(:)
      if (.not.rtdb_put(rtdb,modname//':libxc_funcs',mt_int,
     $    libxc_nfuncs,itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_put failed',1,RTDB_ERR)

      itemp(:) = libxc_family(:)
      if (.not.rtdb_put(rtdb,modname//':libxc_family',mt_int,
     $    libxc_nfuncs,itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_put failed',2,RTDB_ERR)

      itemp(:) = libxc_kind(:)
      if (.not.rtdb_put(rtdb,modname//':libxc_kind',mt_int,
     $    libxc_nfuncs,itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_put failed',3,RTDB_ERR)

      itemp(:) = libxc_flags(:)
      if (.not.rtdb_put(rtdb,modname//':libxc_flags',mt_int,
     $    libxc_nfuncs,itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_put failed',4,RTDB_ERR)

      if (.not.rtdb_put(rtdb,modname//':libxc_facs',mt_dbl,libxc_nfuncs,
     $    libxc_facts))
     $  call errquit(' libxc_rtdbput: rtdb_put failed',5,RTDB_ERR)

      end subroutine


      subroutine nwchem_libxc_rdinput(rtdb,modname)
      use, intrinsic :: iso_c_binding
      implicit none
#include "libxc.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      character(len=*) modname
      integer :: itemp(maxfunc)

      if (.not.rtdb_get(rtdb,modname//':libxc_nfuncs',mt_int,1,
     $    libxc_nfuncs))
     $  call errquit(' libxc_rtdbput: rtdb_get failed',1,RTDB_ERR)

      if (.not.rtdb_get(rtdb,modname//':libxc_funcs',mt_int,
     $    libxc_nfuncs,itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_get failed',1,RTDB_ERR)
      libxc_funcs(:libxc_nfuncs) = int(itemp(:libxc_nfuncs),kind=c_int)

      if (.not.rtdb_get(rtdb,modname//':libxc_family',mt_int,
     $    libxc_nfuncs,itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_get failed',2,RTDB_ERR)
      libxc_family(:libxc_nfuncs) = int(itemp(:libxc_nfuncs),kind=c_int)

      if (.not.rtdb_get(rtdb,modname//':libxc_kind',mt_int,libxc_nfuncs,
     $    itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_get failed',3,RTDB_ERR)
      libxc_kind(:libxc_nfuncs) = int(itemp(:libxc_nfuncs),kind=c_int)

      if (.not.rtdb_get(rtdb,modname//':libxc_flags',mt_int,
     $    libxc_nfuncs,itemp(1)))
     $  call errquit(' libxc_rtdbput: rtdb_get failed',4,RTDB_ERR)
      libxc_flags(:libxc_nfuncs) = int(itemp(:libxc_nfuncs),kind=c_int)

      if (.not.rtdb_get(rtdb,modname//':libxc_facs',mt_dbl,libxc_nfuncs,
     $    libxc_facts))
     $  call errquit(' libxc_rtdbput: rtdb_get failed',5,RTDB_ERR)

      end subroutine
