*
* $Id$
*
      subroutine hnd_propty(rtdb,basis,geom)
c
      implicit none
c
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "int_nbf.fh"
#include "zora.fh"
c
      integer     rtdb   ! runtime database handle
      integer     basis  ! basis set "object" handle
      integer     geom   ! geometry "object" handle
      logical     status
      logical     some
      logical     out 
      integer     iopt
      integer     i,ipol
      character*20 scftype
c
      integer
     &     nodip,noqdp,nootp,nopop,nospin,nodpl,nodhp,nod2hp,
     1     nofpl,nofhp,nof2hp,nosos,noelp,noelf,noelfg,noelfgZ4,noden,
     2     nogiao,nogshift,noston,nospinspin,noaoresponse,noboysloc,
     3     nopmloc,noiboloc,noelpiso
      data nodip,noqdp,nootp,nopop,nospin,nodpl,nodhp,nod2hp,
     1     nofpl,nofhp,nof2hp,nosos,noelp,noelf,noelfg,noelfgZ4,noden,
     2     nogiao,nogshift,noston,nospinspin,noaoresponse,noboysloc,
     3     nopmloc,noiboloc,noelpiso
     4    /  1,    1,    1,    1,    1,    1,    1,    1,    1,
     5       1,    1,    1,    1,    1,    1,    1,    1,    1,
     6       1,    1,    1,    1,    1,    1,    1,    1 /
c
      out   =.false.
      some  =.true.        
      some  =some.or.out
      some  =some.and.ga_nodeid().eq.0

c     Determine scftype
      if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) 
     &     scftype = "RHF"
      if (scftype.eq."RHF") then
           ipol = 1
      else 
           ipol = 2
      endif
c
      status=rtdb_get(rtdb,'prop:dipole',MT_INT,1,nodip)            !done
      status=rtdb_get(rtdb,'prop:quadrupole',MT_INT,1,noqdp)        !done
      status=rtdb_get(rtdb,'prop:octupole',MT_INT,1,nootp)          !done
      status=rtdb_get(rtdb,'prop:mulliken',MT_INT,1,nopop)          !done
      status=rtdb_get(rtdb,'prop:hyperfine',MT_INT,1,nospin)        !done
      status=rtdb_get(rtdb,'prop:loc_boys',MT_INT,1,noboysloc)

      status=rtdb_get(rtdb,'prop:loc_pm',MT_INT,1,nopmloc)

      status=rtdb_get(rtdb,'prop:loc_ibo',MT_INT,1,noiboloc) 

      status=rtdb_get(rtdb,'prop:esp',MT_INT,1,noelp) !done

      status=rtdb_get(rtdb,'prop:espiso',MT_INT,1,noelpiso)         !done
      status=rtdb_get(rtdb,'prop:efield',MT_INT,1,noelf)            !done
      status=rtdb_get(rtdb,'prop:efieldgrad',MT_INT,1,noelfg)       !done
      status=rtdb_get(rtdb,'prop:efieldgradZ4',MT_INT,1,noelfgZ4)   !done
      status=rtdb_get(rtdb,'prop:electrondensity',MT_INT,1,noden)   !done
      status=rtdb_get(rtdb,'prop:stoneanalysis',MT_INT,1,noston)
      status=rtdb_get(rtdb,'prop:giao',MT_INT,1,nogiao)             !done
      status=rtdb_get(rtdb,'prop:gshift',MT_INT,1,nogshift)         !done
      status=rtdb_get(rtdb,'prop:spinspin',MT_INT,1,nospinspin)     !done
      status=rtdb_get(rtdb,'prop:aoresponse',MT_INT,1,noaoresponse) !done
c
      if(nodhp.eq.0.or.nod2hp.eq.0) nodpl=0
      if(nofhp.eq.0.or.nof2hp.eq.0) nofpl=0
      if(nospinspin.eq.0) call hnd_spinspin(rtdb,basis,geom)

c     ----- orbital localization ----

c     Pipek-Mezey localization
      if (nopmloc.eq.0) then
        call localization_driver(rtdb, basis, geom, 1)
      end if
c
c     Boys localization
      if (noboysloc.eq.0) then
        call localization_driver(rtdb, basis, geom, 2)
      end if

c     IBO localization
      if (noiboloc.eq.0) then
        call localization_driver(rtdb, basis, geom, 3)
      end if      
c
c     ----- aoresponse: frequency dependent response -----
c
      iopt = 1 ! default
      if (.not. rtdb_get(rtdb,'prop:newaoresp',MT_INT,1,iopt))  
     &   iopt = 1 ! toggle for old driver for closed shells
      if (ipol.gt.1.or.scftype.eq."UHF") iopt = 1  ! always set new driver for open shells
c
      if (noaoresponse.eq.0) then
        if (iopt.eq.0) then
          if (ga_nodeid().eq.0) then
               write (luout,*)
     &   '*** CALLING OLD AORESP DRIVER FOR CLOSED SHELLS ***'
          endif
          call aoresponse_driver(rtdb, basis, geom)     ! old driver for closed shells only
        else ! iopt .ne. 0
          if (ga_nodeid().eq.0) then
               write (luout,*)
     &   '*** CALLING NEW AORESP DRIVER FOR CLOSED AND OPEN SHELLS ***'
          endif
          call aoresponse_driver_new(rtdb, basis, geom) ! new driver for closed and open shells
        end if
      end if
c
c     ----- dipole moment -----
c
      if(nodip.eq.0.or.nopmloc.eq.0.or.noboysloc.eq.0.or.nodpl.eq.0) 
     &  call hnd_mtpole(rtdb,basis,geom,1)
c
c     ----- quadrupole moment -----
c
      if(noqdp.eq.0) 
     &  call hnd_mtpole(rtdb,basis,geom,2)
c
c     ----- octupole moment -----
c
      if(nootp.eq.0) 
     &  call hnd_mtpole(rtdb,basis,geom,3)
c
c     ----- (hyper)polarizabilities by -sos- -----
c
*     if(nosos.eq.0) call hnd_sosx
c
c     ----- dipole polarizability -----
c     ----- symmetry turned off   -----
c
*     if(nosym.ne.0) call hnd_symoff
c
*     if(nodpl.eq.0.and.nofpl.ne.0) call hnd_dpl
c
c     ----- dipole hyperpolarizability -----
c
*     if(nodhp.eq.0.and.nofhp.ne.0) call hnd_dhp
c
c     ----- dipole second hyperpolarizability -----
c
*     if(nod2hp.eq.0.and.nof2hp.ne.0) call hnd_d2hp
c
c     ----- frequency dependent dipole polarizability -----
c
*     if(nofpl.eq.0) call hnd_fpl
c
c     ----- freqency dependent dipole hyperpolarizability -----
c
*     if(nofhp.eq.0) call hnd_fhp
c
c     ----- frequency dependent dipole second hyperpolarizability -----
c
*     if(nof2hp.eq.0) call hnd_f2hp
c
c     ----- symmetry turned on again -----
c
*      if(nosym.ne.0) call hnd_symon
c
c     ----- mulliken population analysis -----
c
      if(nopop.eq.0) call hnd_mulken(rtdb,basis,geom)
c
c     ----- atomic spin density -----
c
      if(nospin.eq.0.and. (scftype.eq.'UHF')) then  ! only perform for open-shell
       if (.not.do_zora) then
          call hnd_spind(rtdb,basis,geom) ! non-relativistic
       else
        if (ga_nodeid().eq.0) then
        write(LuOut,*)
        call util_print_centered(LuOut,'ZORA NMR Hyperfine', 23, .true.)
        write(LuOut,*)
        endif
        call hnd_hyperfine_ZORA(rtdb,basis,geom)
       end if
      end if
c
c     ----- electrostatic potential -----
c
      if(noelp.eq.0) call hnd_elpmap(rtdb,basis,geom)
c
c     ----- electrostatic potential on a specific electron density surface ---
c
      if(noelpiso.eq.0) call hnd_elpiso(rtdb,basis,geom)
 
c     ----- electric field -----
c
      if(noelf.eq.0) call hnd_elfmap(rtdb,basis,geom)
c
c     ----- electric field gradient -----
c
      if(noelfg.eq.0) call hnd_efgmap(rtdb,basis,geom)
c
c     ----- electric field gradient with ZORA picture change effects -----
c
      if (do_zora. and. (noelfgZ4.eq.0)) then
        if (ga_nodeid().eq.0) then
         write(LuOut,*)
         call util_print_centered(LuOut,'ZORA EFG-Z4', 23, .true.)
         write(LuOut,*)
        endif
        call hnd_efgmap(rtdb,basis,geom)  ! normal EFG
        call hnd_efgmap_Z4(rtdb,basis,geom)  ! with picture change (Z4) effects
      end if
c
c     ----- electron, spin, alpha and beta densities -----
c
      if(noden.eq.0) then
        call hnd_eldmap(rtdb,basis,geom)
        call hnd_eldmap_grid(rtdb,basis,geom)
      end if
c
c     ----- giao -nmr- chemical shifts -----
c
      if((nogiao.eq.0)) then
        if (ga_nodeid().eq.0 .and. (do_zora)) then
         write(LuOut,*)
           call util_print_centered(LuOut,'Scalar ZORA NMR Shielding', 
     &                   23, .true.)
         write(LuOut,*)
        endif
        call hnd_giaox_ZORA(rtdb,basis,geom) ! handles both w/wo zora
      end if
c
c     ------- gshifts are a relativistic property -----
c
      if(nogshift.eq.0 .and. do_zora) then
        if (ga_nodeid().eq.0) then
         write(LuOut,*)
         call util_print_centered(LuOut,'ZORA g-Shift', 23, .true.)
         write(LuOut,*)
        endif
        call hnd_gshift_ZORA(rtdb,basis,geom) ! zora g-shift
      end if
c
c     === deallocate used ga's ===
      if (.not. rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &   ipol=1

      if (nogshift.eq.0 .or. nogiao  .eq.0 .or.
     &    nospin  .eq.0 .or. noelfgZ4.eq.0) then
       if (do_zora) then 
        do i=1,ipol
           if (.not. ga_destroy(g_zora_Kinetic(i))) call errquit(
     &          'prp: ga_destroy failed g_zora_K',0, GA_ERR)
           if (.not. ga_destroy(g_Cifull(i))) call errquit(
     &         'prp: ga_destroy failed ',0, GA_ERR)
        enddo ! ipol
c Do not destroy g_AtNr if doing g-shift because it was not created!!
        status=rtdb_get(rtdb,'prop:gshift'      ,MT_INT,1,nogshift) 
        if (nogshift.ne.0) then
         if (.not. ga_destroy(g_AtNr)) call errquit( ! defined in dft_zora_NMR.F
     &    'prp: ga_destroy failed g_AtNr',0, GA_ERR)   
        endif
       if (.not. ga_destroy(g_Ci)) call errquit( ! defined in dft_zora_utils.F
     &    'prp: ga_destroy failed g_Ci',0, GA_ERR)                 
       end if
      end if ! check-nogshift-nogiao-nospin-noelfgZ4
c
c     ----- stone's distributed multipole analysis -----
c
      noston=1
!     if(noston.eq.0) call hnd_stonex
c
c     ----- boys orbital localization -----
c
      noboysloc=1
!     if(noboysloc.eq.0) call hnd_boyloc
!     if(noboysloc.eq.0) call hnd_mulken(rtdb,basis,geom)
c
      return
      end
