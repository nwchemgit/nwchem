*
* $Id: prp.F,v 1.65 2006-11-16 18:23:01 jochen Exp $
*
      subroutine hnd_propty(rtdb,basis,geom)

c ... jochen: added some aoresponse stuff here

      implicit double precision (a-h,o-z)
#include "nwc_const.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "int_nbf.fh"
      integer     rtdb
      logical     status
      logical     some
      logical     out 
      data nodip,noqdp,nootp,nopop,nospin,noloc,nodpl,nodhp,nod2hp,
     1     nofpl,nofhp,nof2hp,nosos,noelp,noelf,noelfg,noden,
     2     nogiao,noston,nospinspin, noaoresponse
     3    /  1,    1,    1,    1,    1,    1,    1,    1,    1,
     4       1,    1,    1,    1,    1,    1,    1,    1, 
     5       1,    1,    1,    1/
c
      out   =.false.
      some  =.true.        
      some  =some.or.out
      some  =some.and.ga_nodeid().eq.0
c
      status=rtdb_get(rtdb,'prop:dipole',MT_INT,1,nodip)            !done
      status=rtdb_get(rtdb,'prop:quadrupole',MT_INT,1,noqdp)        !done
      status=rtdb_get(rtdb,'prop:octupole',MT_INT,1,nootp)          !done
      status=rtdb_get(rtdb,'prop:mulliken',MT_INT,1,nopop)          !done
      status=rtdb_get(rtdb,'prop:hyperfine',MT_INT,1,nospin)        !done
      status=rtdb_get(rtdb,'prop:boyslocalization',MT_INT,1,noloc)
      status=rtdb_get(rtdb,'prop:esp',MT_INT,1,noelp)               !done
      status=rtdb_get(rtdb,'prop:efield',MT_INT,1,noelf)            !done
      status=rtdb_get(rtdb,'prop:efieldgrad',MT_INT,1,noelfg)       !done
      status=rtdb_get(rtdb,'prop:electrondensity',MT_INT,1,noden)   !done
      status=rtdb_get(rtdb,'prop:stoneanalysis',MT_INT,1,noston)
      status=rtdb_get(rtdb,'prop:giao',MT_INT,1,nogiao)             !done
      status=rtdb_get(rtdb,'prop:spinspin',MT_INT,1,nospinspin)     !done
c ... jochen:
      status=rtdb_get(rtdb,'prop:aoresponse',MT_INT,1,noaoresponse) !done
c
      if(nodhp.eq.0.or.nod2hp.eq.0) nodpl=0
      if(nofhp.eq.0.or.nof2hp.eq.0) nofpl=0
      if(nospinspin.eq.0) call hnd_spinspin(rtdb,basis,geom)

c     ----- aoresponse: frequency dependent response -----
c
      if (noaoresponse.eq.0)
     +   call aoresponse_driver(rtdb, basis, geom)
c
c     ----- dipole moment -----
c
      if(nodip.eq.0.or.noloc.eq.0.or.nodpl.eq.0) 
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
      if(nospin.eq.0) call hnd_spind(rtdb,basis,geom)
c
c     ----- electrostatic potential -----
c
      if(noelp.eq.0) call hnd_elpmap(rtdb,basis,geom)
 
c     ----- electric field -----
c
      if(noelf.eq.0) call hnd_elfmap(rtdb,basis,geom)
c
c     ----- electric field gradient -----
c
      if(noelfg.eq.0) call hnd_efgmap(rtdb,basis,geom)
c
c     ----- electron and spin densities -----
c
      if(noden.eq.0) call hnd_eldmap(rtdb,basis,geom)
c
c     ----- giao -nmr- chemical shifts -----
c
      if(nogiao.eq.0) call hnd_giaox(rtdb,basis,geom)
c
c     ----- stone's distributed multipole analysis -----
c
      noston=1
!     if(noston.eq.0) call hnd_stonex
c
c     ----- boys orbital localization -----
c
      noloc=1
!     if(noloc.eq.0) call hnd_boyloc
!     if(noloc.eq.0) call hnd_mulken(rtdb,basis,geom)
c
      return
      end
