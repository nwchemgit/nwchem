      subroutine cosmo_input(rtdb)
*
* $Id: cosmo.F,v 1.1 2004-07-21 17:28:26 bert Exp $
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer    mxatom
      parameter (mxatom=nw_max_atom)
c
      integer rtdb
      integer iat
      integer nfield
      integer nrad   
      integer irad   
      character*255 field
      logical ocosmo
      logical do_cosmo
      integer nap
      integer iap
c
c
      double precision radius
      double precision rad    
      double precision dielec
      double precision rsolv
      integer          iscren
      integer          ifscrn
      integer          minbem
      integer          maxbem
      integer          ificos
      integer          lineq
c
      logical odbug
      logical osome
c
      common/hnd_lcapid/nap,iap
      common/hnd_cosrad/radius(mxatom)
c
c     ----- /lcapid/ -----
c
      nap = ga_nnodes()
      iap = ga_nodeid() + 1
c
      odbug=.false.
      odbug=odbug.and.(iap.eq.1)
c
c     ----- read input data -----
c
c     ----- defaults -----
c
      do_cosmo = .true.
      dielec   =78.4d+00
      iscren   =0
      ifscrn   =2
      minbem   =2
      maxbem   =3
      ificos   =0
      lineq    =1
      rsolv    =0.0d+00
      do iat=1,mxatom
         radius(iat)=0.0d+00
      enddo
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('cosmo_input: unexpected eof',911, INPUT_ERR)
      nfield = inp_n_field()
 20   if (.not. inp_a(field))
     $     call errquit('cosmo_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         go to 40
c
c   --- dielec
c
      else if(inp_compare(.false.,'dielec',field)) then
         if(.not. inp_f(dielec))
     $        call errquit('cosmo_input: dielec',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:dielec',mt_dbl,1,dielec))
     $        call errquit('cosmo_input: rtdb put failed',911,
     &       RTDB_ERR)
c
c   --- off  
c
      else if(inp_compare(.false.,'off',field)) then
         do_cosmo=.false.
c
c   --- rsolv
c
      else if(inp_compare(.false.,'rsolv',field)) then
         if(.not. inp_f(rsolv))
     $        call errquit('cosmo_input: rsolv',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:rsolv',mt_dbl,1,rsolv))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- radius
c
      else if(inp_compare(.false.,'radius',field)) then
         irad=0
         nfield=nfield-1
   30    if(nfield.gt.0) then
            if(odbug) then
               write(Luout,*) 'nfield = ',nfield
            endif
            if(inp_f(rad)) then
               irad=irad+1
               radius(irad)=rad
               if(odbug) then
                  write(Luout,*) 'radius value found = ',rad
               endif
               nfield=nfield-1
               go to 30
            else
               if(odbug) then
                  write(Luout,*) 'not a radius value'
               endif
               nrad=irad  
               if (.not. rtdb_put(rtdb,'cosmo:radius',
     $                            mt_dbl,mxatom,radius))
     $            call errquit('cosmo_input: rtdb put failed',911,
     &       RTDB_ERR)
               if(odbug) then
                  write(Luout,*) 'nrad = ',nrad
               endif
               go to 20
            endif
         else
            if(.not.inp_read())
     $         call errquit('cosmo_input: unexpected eof',911,
     &       INPUT_ERR)
            nfield = inp_n_field()
            if(odbug) then
               write(Luout,*) 'new line with # fields = ',nfield
            endif
            go to 30
         endif
c
c   --- iscren
c
      else if(inp_compare(.false.,'iscren',field)) then
         if(.not. inp_i(iscren))
     $        call errquit('cosmo_input: iscren',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:iscren',mt_int,1,iscren))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- ifscrn
c
      else if(inp_compare(.false.,'ifscrn',field)) then
         if(.not. inp_i(ifscrn))
     $        call errquit('cosmo_input: ifscrn',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:ifscrn',mt_int,1,ifscrn))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- minbem
c
      else if(inp_compare(.false.,'minbem',field)) then
         if(.not. inp_i(minbem))
     $        call errquit('cosmo_input: minbem',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:minbem',mt_int,1,minbem))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- maxbem
c
      else if(inp_compare(.false.,'maxbem',field)) then
         if(.not. inp_i(maxbem))
     $        call errquit('cosmo_input: maxbem',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:maxbem',mt_int,1,maxbem))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- ificos
c
      else if(inp_compare(.false.,'ificos',field)) then
         if(.not. inp_i(ificos))
     $        call errquit('cosmo_input: ificos',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:ificos',mt_int,1,ificos))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- lineq
c
      else if(inp_compare(.false.,'lineq',field)) then
         if(.not. inp_i(lineq))
     $        call errquit('cosmo_input: lineq',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:lineq',mt_int,1,lineq))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
      endif
c
      go to 10
c
   40 continue
c
c     Let SCF/DFT know cosmo has to be calculated
c
      if (do_cosmo) then 
         if (.not. rtdb_put(rtdb,'slv:cosmo',mt_log,1,.true.))
     $        call errquit('cosmo_input: rtdb_put failed',0, RTDB_ERR)
      else
         if (.not. rtdb_put(rtdb,'slv:cosmo',mt_log,1,.false.))
     $        call errquit('cosmo_input: rtdb_put failed',0, RTDB_ERR)
      endif
c
c     ----- print ... -----
c
      if(odbug) then
         write(Luout,*) '-cosmo input-'
         write(Luout,*) 'dielec     = ',dielec
         write(Luout,*) 'iscren     = ',iscren
         write(Luout,*) 'ifscrn     = ',ifscrn
         write(Luout,*) 'minbem     = ',minbem
         write(Luout,*) 'maxbem     = ',maxbem
         write(Luout,*) 'ificos     = ',ificos
         write(Luout,*) 'lineq      = ',lineq
         write(Luout,*) 'rsolv      = ',rsolv
         write(Luout,*) 'radius     = ',(radius(irad),irad=1,nrad)
         write(Luout,*) 'do cosmo   = ',do_cosmo
c
         call cosmo_input_print(rtdb)
      endif
c
      call util_flush(Luout)
c
      return
      end
      subroutine cosmo_tidy()
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "prop.fh"
#include "geom.fh"
c
      if (.not. geom_destroy(cosmo_geom_efc)) call errquit
     $     ('cosmo_tidy: geom destroy failed', 0, GEOM_ERR)
c
      return
      end
c
      subroutine cosmo_input_print(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
c
      integer    mxatom
      parameter (mxatom=nw_max_atom)
c
      integer rtdb
      integer iat
      integer irad
      integer nrad
      logical cosmo
      integer nap
      integer iap
c
      double precision dielec
      double precision rsolv
      double precision radius
      integer          iscren
      integer          ifscrn
      integer          minbem
      integer          maxbem
      integer          ificos
      integer          lineq
c
      logical odbug
c
      common/hnd_lcapid/nap,iap
      common/hnd_cosrad/radius(mxatom)
c
      odbug=.false.
      odbug=odbug.and.(iap.eq.1)
c
      if(odbug) then
         write(Luout,*) 'in cosmo_input_print ...'
      endif
c
c     ----- defaults -----
c
      dielec   =78.4d+00
      iscren   =0
      ifscrn   =2
      minbem   =2
      maxbem   =3
      ificos   =0
      lineq    =1
      rsolv    =0.0d+00
      do iat=1,mxatom
         radius(iat)=0.0d+00
      enddo
c
      if(odbug) then
         write(Luout,*) '-cosmo input defaults-'
         write(Luout,*) 'dielec     = ',dielec
         write(Luout,*) 'iscren     = ',iscren
         write(Luout,*) 'ifscrn     = ',ifscrn
         write(Luout,*) 'minbem     = ',minbem
         write(Luout,*) 'maxbem     = ',maxbem
         write(Luout,*) 'ificos     = ',ificos
         write(Luout,*) 'lineq      = ',lineq
         write(Luout,*) 'rsolv      = ',rsolv
      endif
c
c   --- cosmo
c
      if (.not. rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo))
     $     call errquit('cosmo_input_print: rtdb_get failed',911,
     &       RTDB_ERR)
         if(cosmo) then
            if(odbug) then
               write(Luout,*) 'cosmo  ... found on rtdb'
            endif
         else
            if(odbug) then
               write(Luout,*) 'cosmo  not found on rtdb'
            endif
         endif
c
c   --- dielec
c
      if (.not. rtdb_get(rtdb,'cosmo:dielec',mt_dbl,1,dielec)) then
         if(odbug) then
            write(Luout,*) 'dielec ... from defaults = ',dielec
         endif
      else
         if(odbug) then
            write(Luout,*) 'dielec ... from  -rtdb-  = ',dielec
         endif
      endif
c
c   --- rsolv
c
      if (.not. rtdb_get(rtdb,'cosmo:rsolv',mt_dbl,1,rsolv)) then
         if(odbug) then
            write(Luout,*) 'rsolv  ... from defaults = ',rsolv 
         endif
      else
         if(odbug) then
            write(Luout,*) 'rsolv  ... from  -rtdb-  = ',rsolv 
         endif
      endif
c
c   --- iscren
c
      if (.not. rtdb_get(rtdb,'cosmo:iscren',mt_int,1,iscren)) then
         if(odbug) then
            write(Luout,*) 'iscren ... from defaults = ',iscren 
         endif
      else
         if(odbug) then
            write(Luout,*) 'iscren ... from  -rtdb-  = ',iscren
         endif
      endif
c
c   --- ifscrn
c
      if (.not. rtdb_get(rtdb,'cosmo:ifscrn',mt_int,1,ifscrn)) then
         if(odbug) then
            write(Luout,*) 'ifscrn ... from defaults = ',ifscrn 
         endif
      else
         if(odbug) then
            write(Luout,*) 'ifscrn ... from  -rtdb-  = ',ifscrn
         endif
      endif
c
c   --- minbem
c
      if (.not. rtdb_get(rtdb,'cosmo:minbem',mt_int,1,minbem)) then
         if(odbug) then
            write(Luout,*) 'minbem ... from defaults = ',minbem
         endif
      else
         if(odbug) then
            write(Luout,*) 'minbem ... from  -rtdb-  = ',minbem
         endif
      endif
c
c   --- maxbem
c
      if (.not. rtdb_get(rtdb,'cosmo:maxbem',mt_int,1,maxbem)) then
         if(odbug) then
            write(Luout,*) 'maxbem ... from defaults = ',maxbem
         endif
      else
         if(odbug) then
            write(Luout,*) 'maxbem ... from  -rtdb-  = ',maxbem
         endif
      endif
c
c   --- ificos
c
      if (.not. rtdb_get(rtdb,'cosmo:ificos',mt_int,1,ificos)) then
         if(odbug) then
            write(Luout,*) 'ificos ... from defaults = ',ificos
         endif
      else
         if(odbug) then
            write(Luout,*) 'ificos ... from  -rtdb-  = ',ificos
         endif
      endif
c
c   --- lineq
c
      if (.not. rtdb_get(rtdb,'cosmo:lineq',mt_int,1,lineq)) then
         if(odbug) then
            write(Luout,*) 'lineq  ... from defaults = ',lineq 
         endif
      else
         if(odbug) then
            write(Luout,*) 'lineq  ... from  -rtdb-  = ',lineq 
         endif
      endif
c
c   --- radius
c
      nrad=mxatom
      if (.not. rtdb_get(rtdb,'cosmo:radius',mt_dbl,nrad,radius)) then
         if(odbug) then
            write(Luout,*) 'radius ... from defaults = ',
     $                     (radius(irad),irad=1,nrad)
         endif
      else
         if(odbug) then
            write(Luout,*) 'radius ... from  -rtdb-  = ',
     $                     (radius(irad),irad=1,nrad)
         endif
      endif
c
c   --- print
c
      if(odbug) then
         write(Luout,*) '-cosmo input from -rtdb- -'
         write(Luout,*) 'dielec     = ',dielec
         write(Luout,*) 'iscren     = ',iscren
         write(Luout,*) 'ifscrn     = ',ifscrn
         write(Luout,*) 'minbem     = ',minbem
         write(Luout,*) 'maxbem     = ',maxbem
         write(Luout,*) 'ificos     = ',ificos
         write(Luout,*) 'lineq      = ',lineq
         write(Luout,*) 'rsolv      = ',rsolv
      endif
c
      call util_flush(Luout)
c
      return
      end
      subroutine cosmo_initialize(rtdb,geom,basis,oprint)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "prop.fh"
#include "geom.fh"
#include "bas.fh"
c
      integer    mxatom
      parameter (mxatom=nw_max_atom)
      parameter (mxefc =3*7*7*7)
c
      logical status
      integer rtdb
      integer geom
      integer basis
      integer iat
      logical cosmo
c
      double precision dielec
      double precision rsolv
      double precision radius
      integer          iscren
      integer          ifscrn
      integer          minbem
      integer          maxbem
      integer          ificos
      integer          lineq
c
      logical odbug
      logical osome
      logical oprint
      integer me
c
      character*4 keyblk,keycos,icos,iefc
      character*8 blank
      character*8 efclab
      common/hnd_iofile/ir,iw,ip
      common/hnd_lcapid/nap,iap
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_cosrad/radius(mxatom)
      common/hnd_cosprt/iprt
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_efcpar/efcc(3,mxefc),efcz(mxefc),efclab(mxefc),nefc
      common/hnd_efcopt/iefc	
      character*16 atom_tags(mxatom)
      dimension vdwr(103)
      dimension vander(36)
c
      data blank  /'        '/
      data keyblk /'    '/
      data keycos /' COS'/
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pt5    /0.5d+00/
      data radtol /0.1d+00/
      data mxpass /6/
C
C     ----- van der waals radii in angstroms	-----
C
      data def  /1.80D+00/
      data vdwr /
     1   0.80,0.49,0.00,0.00,0.00,1.65,1.55,1.50,1.50,0.00,
     2   2.30,1.70,2.05,2.10,1.85,1.80,1.80,0.00,2.80,2.75,
     3   0.00,0.00,1.20,0.00,0.00,0.00,2.70,0.00,0.00,0.00,
     4   0.00,0.00,0.00,1.90,1.90,0.00,0.00,0.00,0.00,1.55,
     5   0.00,1.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     6   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     7   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     8   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     9   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.65,
     2   0.00,0.00,0.00/
C
C     vdw radii from Merz/Kollman/Singh
C
      data (vander(i),i=1,36)
     1                  /1.20D+00,1.20D+00,1.37D+00,1.45D+00,
     2 1.45D+00,1.50D+00,1.50D+00,1.40D+00,1.35D+00,1.30D+00,
     3 1.57D+00,1.36D+00,1.24D+00,1.17D+00,1.80D+00,1.75D+00,
     4 1.70D+00,19*0.0D+00/
c
c     VDW RADII FROM GAVEZZOTTI (J.AM.CHEM.SOC. 105, 5220 (1983))
c     SUPPLEMENTED BY VALUES FOR BORON AND THIRD AND FOURTH ROWS
c     FROM M.A.SPACKMAN (J.CHEM.PHYS. 85, 6579 (1986))
c
c     DATA (VANDER(I),I=1,36)
c    1                  /1.20D+00,0.00D+00,0.00D+00,0.00D+00,
c    2 1.85D+00,1.50D+00,1.50D+00,1.40D+00,1.35D+00,0.00D+00,
c    3 0.00D+00,0.00D+00,2.07D+00,2.05D+00,1.96D+00,1.89D+00,
c    4 1.80D+00,19*0.0D+00/
c
c     VDW RADII FROM BRENEMAN & WIBERG
c
c     DATA (VANDER(I),I=1,36)
c    1                  /1.45D+00,1.45D+00,1.50D+00,1.50D+00,
c    2 1.50D+00,1.50D+00,1.70D+00,1.70D+00,1.70D+00,1.70D+00,
c    3 2.00D+00,2.00D+00,2.00D+00,2.00D+00,2.00D+00,2.00D+00,
c    4 2.00D+00,19*0.0D+00/
c
c
c     ----- /iofile/ -----
c
      ir=Luin
      iw=Luout
      ip=7
c
c     ----- /lcapid/ -----
c
      me=ga_nodeid()
      nap = ga_nnodes()
      iap = ga_nodeid() + 1
c
      odbug=.false.
      odbug=odbug.and.(iap.eq.1)
      osome=oprint
      osome=osome.or.odbug
      osome=osome.and.(iap.eq.1)
c
      if(odbug) then
         write(Luout,*) 'in cosmo_initialize ...'
      endif
c
c     Check if spherical basis sets are requested
c
      if(bas_is_spherical(basis)) call errquit(
     .   ' cosmo does not handle spherical basis yet',0, BASIS_ERR)
c
c     Create a geometry for the COSMO efc charges, to be used later on
c
      if (.not. geom_create(cosmo_geom_efc, 'cosmo efc geom'))
     $     call errquit('cosmo_initialize: geom_create failed', 0,
     &       GEOM_ERR)
c
c     ----- /molxyz/ /molnuc/ -----
c
      status=geom_ncent(geom,nat)
      if(.not.status) call errquit
     $    (' cosmo_initialize: geom_ncent    failed.',911, GEOM_ERR)
      status=geom_cart_get(geom,nat,atom_tags,c,zan)
      if(.not.status) call errquit
     $    (' cosmo_initialize: geom_cart_get failed.',911, GEOM_ERR)
      do iat=1,nat
         nuc(iat)=nint(zan(iat))
      enddo
      if(odbug) then
         do iat=1,nat
            write(iw,*) ' ... i,Z,c = ',iat,zan(iat),c(1,iat),c(2,iat),
     $                                               c(3,iat)
         enddo
      endif
c
c     ----- defaults -----
c
      dielec   =78.4d+00
      iscren   =0
      ifscrn   =2
      minbem   =2
      maxbem   =3
      ificos   =0
      lineq    =1
      rsolv    =0.0d+00
      do iat=1,mxatom
         radius(iat)=0.0d+00
      enddo
      do i=1,36
         vdwr(i)=vander(i) 
      enddo
      do iat=1,nat
         rad=vdwr(nuc(iat))
         if (geom_any_ecp(geom).and.geom_ecp_get(geom,iat)) then
            rad = zero
            status = geom_tag_to_default_mass(atom_tags(iat),rad)
         endif
         if(rad.eq.zero) then
            rad=def
         endif
         radius(iat)=rad
      enddo
c
      if(odbug) then
         write(Luout,*) '-cosmo input defaults-'
         write(Luout,*) 'dielec     = ',dielec
         write(Luout,*) 'iscren     = ',iscren
         write(Luout,*) 'ifscrn     = ',ifscrn
         write(Luout,*) 'minbem     = ',minbem
         write(Luout,*) 'maxbem     = ',maxbem
         write(Luout,*) 'ificos     = ',ificos
         write(Luout,*) 'lineq      = ',lineq
         write(Luout,*) 'rsolv      = ',rsolv
      endif
c
c   --- cosmo
c
      if (.not. rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo))
     $     call errquit('cosmo_initialize: rtdb_get failed',911,
     &       RTDB_ERR)
         if(cosmo) then
            if(odbug) then
               write(Luout,*) 'cosmo  ... found on rtdb'
            endif
         else
            if(odbug) then
               write(Luout,*) 'cosmo  not found on rtdb'
            endif
         endif
c
c   --- dielec
c
      if (.not. rtdb_get(rtdb,'cosmo:dielec',mt_dbl,1,dielec)) then
         if(odbug) then
            write(Luout,*) 'dielec ... from defaults = ',dielec
         endif
      else
         if(odbug) then
            write(Luout,*) 'dielec ... from  -rtdb-  = ',dielec
         endif
      endif
c
c   --- rsolv
c
      if (.not. rtdb_get(rtdb,'cosmo:rsolv',mt_dbl,1,rsolv)) then
         if(odbug) then
            write(Luout,*) 'rsolv  ... from defaults = ',rsolv 
         endif
      else
         if(odbug) then
            write(Luout,*) 'rsolv  ... from  -rtdb-  = ',rsolv 
         endif
      endif
c
c   --- iscren
c
      if (.not. rtdb_get(rtdb,'cosmo:iscren',mt_int,1,iscren)) then
         if(odbug) then
            write(Luout,*) 'iscren ... from defaults = ',iscren 
         endif
      else
         if(odbug) then
            write(Luout,*) 'iscren ... from  -rtdb-  = ',iscren
         endif
      endif
c
c   --- ifscrn
c
      if (.not. rtdb_get(rtdb,'cosmo:ifscrn',mt_int,1,ifscrn)) then
         if(odbug) then
            write(Luout,*) 'ifscrn ... from defaults = ',ifscrn 
         endif
      else
         if(odbug) then
            write(Luout,*) 'ifscrn ... from  -rtdb-  = ',ifscrn
         endif
      endif
c
c   --- minbem
c
      if (.not. rtdb_get(rtdb,'cosmo:minbem',mt_int,1,minbem)) then
         if(odbug) then
            write(Luout,*) 'minbem ... from defaults = ',minbem
         endif
      else
         if(odbug) then
            write(Luout,*) 'minbem ... from  -rtdb-  = ',minbem
         endif
      endif
c
c   --- maxbem
c
      if (.not. rtdb_get(rtdb,'cosmo:maxbem',mt_int,1,maxbem)) then
         if(odbug) then
            write(Luout,*) 'maxbem ... from defaults = ',maxbem
         endif
      else
         if(odbug) then
            write(Luout,*) 'maxbem ... from  -rtdb-  = ',maxbem
         endif
      endif
c
c   --- ificos
c
      if (.not. rtdb_get(rtdb,'cosmo:ificos',mt_int,1,ificos)) then
         if(odbug) then
            write(Luout,*) 'ificos ... from defaults = ',ificos
         endif
      else
         if(odbug) then
            write(Luout,*) 'ificos ... from  -rtdb-  = ',ificos
         endif
      endif
c
c   --- lineq
c
      if (.not. rtdb_get(rtdb,'cosmo:lineq',mt_int,1,lineq)) then
         if(odbug) then
            write(Luout,*) 'lineq  ... from defaults = ',lineq 
         endif
      else
         if(odbug) then
            write(Luout,*) 'lineq  ... from  -rtdb-  = ',lineq 
         endif
      endif
c
c   --- radius
c
      nrad=mxatom
      if (.not. rtdb_get(rtdb,'cosmo:radius',mt_dbl,nrad,radius)) then
         if(odbug) then
            write(Luout,*) 'radius ... from defaults = ',
     $                     (radius(irad),irad=1,nrad)
         endif
      else
         if(odbug) then
            write(Luout,*) 'radius ... from  -rtdb-  = ',
     $                     (radius(irad),irad=1,nrad)
         endif
      endif
c
c     ----- checks ... -----
c
      if(dielec.le.zero) then
         dielec=78.4d+00
      endif
c
      IF(OSOME) THEN
         WRITE(IW,9999)
         WRITE(IW,9997) DIELEC
      ENDIF
C
C     ----- CHARGE SCREENING APPROACH ... -----
C
      IF(OSOME) THEN
         WRITE(IW,9983) IFSCRN
      ENDIF

      IF(IFSCRN.NE.1.AND.IFSCRN.NE.2) THEN
         ISCREN=-1
      ENDIF
C
C     ----- SCREENING FACTOR -----
C
      IF(ISCREN.LT.0) THEN
         SCREEN=ONE
         IF (OSOME) WRITE(IW,9984) SCREEN
      ELSEIF(ISCREN.GT.0) THEN
         SCREEN=(DIELEC-ONE)/(DIELEC+PT5)
         IF (OSOME) WRITE(IW,9996) SCREEN
      ELSE
         SCREEN=(DIELEC-ONE)/(DIELEC    )
         IF (OSOME) WRITE(IW,9995) SCREEN
      ENDIF
C
C     ----- TESSELATION -----
C
      IF(MINBEM.LT.1.OR.MINBEM.GT.MXPASS) THEN
         MINBEM=2
      ENDIF
      IF(MAXBEM.LT.1.OR.MAXBEM.GT.MXPASS) THEN
         MAXBEM=4
      ENDIF
      IF(MAXBEM.LT.MINBEM) THEN
         MAXBEM=MINBEM+2
      ENDIF
C
      IFICOS=ABS(IFICOS)
C
      ICOS=KEYCOS
C
      IF(OSOME) THEN
         WRITE(IW,9993) LINEQ
         WRITE(IW,9991) MINBEM
         WRITE(IW,9992) MAXBEM
         IF(IFICOS.EQ.0) THEN
            WRITE(IW,9990)
         ELSE
            WRITE(IW,9989)
         ENDIF
         WRITE(IW,9994) RSOLV
         WRITE(IW,9988)
      ENDIF
      DO IAT=1,NAT
         RADIUS(IAT)=ABS(RADIUS(IAT))
         IF(RADIUS(IAT).LT.RADTOL) THEN
            RADIUS(IAT)=ZERO
         ENDIF
         IF(OSOME) THEN
            WRITE(IW,9987) IAT,NUC(IAT),RADIUS(IAT)
         ENDIF
      ENDDO
C
C     ----- CREATE BOUNDARY ELEMENTS OF CAVITY -----
C
      IF(OSOME) THEN
         IPRT=1
      ELSE
         IPRT=0
      ENDIF

*     **********************************************
*     **** do call in parallel ? -changed by EJB****
*     **********************************************
c-    IF(IAP.EQ.1) THEN
c-        CALL HND_COSSET(RTDB)
c-    ENDIF
      CALL HND_COSSET(RTDB)
      IPRT=0
*     **********************************************
C
C     ----- FOR NOW ... NO INDUCED CHARGES ... -----
C
      NEFC=0
      DO IEF=1,MXEFC
         EFCLAB(IEF)=BLANK
         EFCZ(  IEF)=ZERO
         DO I=1,3
            EFCC(I,IEF)=ZERO
         ENDDO
      ENDDO
      IEFC=KEYBLK
c
c     ----- done -----
c
      if(osome) then
         write(iw,9986)
         write(iw,9985)
      endif
c
      if(me.eq.0) call util_flush(Luout)
c
      return
 9999 FORMAT(/,10X,15(1H-),
     1       /,10X,'-cosmo- solvent',
     2       /,10X,15(1H-))
 9998 FORMAT(' no -cosmo- solvent.')
 9997 FORMAT(' dielectric constant -eps-  = ',F6.2)
 9996 FORMAT(' screen = (eps-1)/(eps+1/2) = ',F9.5)
 9995 FORMAT(' screen = (eps-1)/(eps    ) = ',F9.5)
 9994 FORMAT(' solvent radius (ang.)      = ',F7.3)
 9993 FORMAT(' -lineq- algorithm          = ',I3)
 9992 FORMAT(' -bem- high level           = ',I3)
 9991 FORMAT(' -bem- low  level           = ',I3)
 9990 FORMAT(' -bem- from -octahedral-')
 9989 FORMAT(' -bem- from -icosahedral-')
 9988 FORMAT(' atomic radii = ',/,' --------------')
 9987 FORMAT(i5,i5,f7.3)
 9986 FORMAT(' ...... end of -cosmo- initialization ......')
 9985 FORMAT(/)
 9984 FORMAT(' screen =  one              = ',F9.5)
 9983 FORMAT(' charge screening approach  = ',I3)
      end
      SUBROUTINE HND_COSSET(RTDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
C
C              ----- STARTING FROM -ICOSAHEDRON- -----
C
C     PASS, NAPEX, NFACE, ERROR =   0      12      20      20
C     PASS, NAPEX, NFACE, ERROR =   1      42      80     100    0.4982
C     PASS  NAPEX, NFACE, ERROR =   2     162     320     420    0.1848
C     PASS  NAPEX, NFACE, ERROR =   3     642    1280    1700    0.0523
C     PASS  NAPEX, NFACE, ERROR =   4    2562    5120    6820    0.0135
C     PASS  NAPEX, NFACE, ERROR =   5   10242   20480   27300    0.0034
C
C              ----- STARTING FROM -OCTAHEDRON-  -----
C
C     PASS, NAPEX, NFACE, ERROR =   0       6       8       8
C     PASS, NAPEX, NFACE, ERROR =   1      18      32      40    0.8075
C     PASS  NAPEX, NFACE, ERROR =   2      66     128     168    0.4557
C     PASS  NAPEX, NFACE, ERROR =   3     258     512     680    0.1619
C     PASS  NAPEX, NFACE, ERROR =   4    1026    2048    2728    0.0451
C     PASS  NAPEX, NFACE, ERROR =   5    4098    8192   10920    0.0116
C     PASS  NAPEX, NFACE, ERROR =   6   16386   32768   43688    0.0029
C
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      integer rtdb
      integer geom
      logical status
c
      PARAMETER (MXFACE=10920)
      PARAMETER (MXAPEX= 4098)
      PARAMETER (MXEFC =3*7*7*7)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      CHARACTER*8 EFCLAB
      CHARACTER*8 ERRMSG
      LOGICAL     DBUG
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_COSMOD/LINEQ,MINBEM,MAXBEM,IFICOS,IFSCRN
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-COSSET-'/
c
      character*16 atom_tags(mxatom)
C
C     XYZATM(3,MXATOM)
C       RATM(  MXATOM)
C       NSPA(  MXATOM)
C       NPPA(  MXATOM)
C     IJKFAC(3,MXFACE)
C     XYZSEG(3,MXFACE)
C     IJKSEG(  MXFACE)
C     INSSEG(  MXFACE,MXATOM)
C     XYZSPA(3,MXFACE,MXATOM)
C     IJKSPA(  MXFACE,MXATOM)
C     NUMPPS(  MXFACE,MXATOM)
C       APEX(3,MXAPEX)
C
      DBUG=.FALSE.
      DBUG=DBUG.AND.(IAP.EQ.1)
C
      IF(DBUG) THEN
         WRITE(IW,9999)
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
      IF(IFICOS.EQ.0.AND.MAXBEM.GT.6) THEN
         WRITE(IW,*) '-MAXBEM- TOO LARGE FOR PARAMETERS IN -COSSET-'
CBERT    CALL HND_FIFLSH(IW)
CBERT    CALL HND_HNDERR(3,ERRMSG)
      ELSEIF(IFICOS.NE.0.AND.MAXBEM.GT.7) THEN
         WRITE(IW,*) '-MAXBEM- TOO LARGE FOR PARAMETERS IN -COSSET-'
CBERT    CALL HND_FIFLSH(IW)
CBERT    CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      NFAC=MXFACE
      NAPX=MXAPEX
C
C     ----- PARTITION MEMORY -----
C
      I10 =1      
      I20 =I10 +3*NAT
      I30 =I20 +  NAT
      I40 =I30 +  NAT
      I50 =I40 +  NAT
      I60 =I50 +3*NFAC
      I70 =I60 +3*NFAC
      I80 =I70 +  NFAC
      I90 =I80 +  NFAC*NAT
      I100=I90 +3*NFAC*NAT
      I110=I100+  NFAC*NAT
      I120=I110+  NFAC*NAT
      I130=i120+3*NAPX
      LAST=I130-1
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      NEED = LAST
      IF(.NOT.MA_PUSH_GET(MT_DBL,NEED,'MEM INIT:cosmo:hnd_cosset:1',
     &    I_INIT,INIT))
     & CALL ERRQUIT('HND_COSSET, MALLOC OF INIT  FAILED',911)
C
      I10 =INIT   
      I20 =I10 +3*NAT
      I30 =I20 +  NAT
      I40 =I30 +  NAT
      I50 =I40 +  NAT
      I60 =I50 +3*NFAC
      I70 =I60 +3*NFAC
      I80 =I70 +  NFAC
      I90 =I80 +  NFAC*NAT
      I100=I90 +3*NFAC*NAT
      I110=I100+  NFAC*NAT
      I120=I110+  NFAC*NAT
      I130=i120+3*NAPX
      LAST=I130-1
C
C     ----- GET -COSMO- SURFACE -----
C
      CALL HND_COSSRF(NAT,C,NUC,NAT,NFAC,NAPX,
     1                X(I10),X(I20),X(I30),X(I40),
     2                X(I50),X(I60),X(I70),
     3                X(I80),X(I90),X(I100),X(I110),X(I120),rtdb)
C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_COSSET, MA_POP_STACK OF INIT FAILED',911)
C
      RETURN
 9999 FORMAT(/,10X,15(1H-),
     1       /,10X,'-COSMO- SURFACE',
     2       /,10X,15(1H-))
      END
      SUBROUTINE HND_COSSRF(NAT,C,NUC,MXATM,MXFAC,MXAPX,
     1                  XYZATM,RATM,NSPA,NPPA,
     2                  IJKFAC,XYZSEG,IJKSEG,INSSEG,
     3                  XYZSPA,IJKSPA,NUMPPS,APEX,rtdb)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
#include "nwc_const.fh"
#include "rtdb.fh"
      integer rtdb
c
      PARAMETER (BOHR=0.529177249D+00)
      PARAMETER (MXATOM=NW_MAX_ATOM)
      LOGICAL    SOME
      LOGICAL    OUT
      LOGICAL    MORE
      LOGICAL    DBUG
      LOGICAL    INSSEG
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_COSDAT/ADIAG,DSURF,DVOL,SRFMOL,VOLMOL 
      COMMON/HND_COSPAR/DIELEC,SCREEN,RSOLV
      COMMON/HND_COSRAD/RADIUS(MXATOM)
      DIMENSION      C(3,*    )
      DIMENSION    NUC(  *    )
      DIMENSION XYZATM(3,MXATM)
      DIMENSION   RATM(  MXATM)
      DIMENSION   NSPA(  MXATM)
      DIMENSION   NPPA(  MXATM)
      DIMENSION IJKFAC(3,MXFAC)
      DIMENSION XYZSEG(3,MXFAC)
      DIMENSION IJKSEG(  MXFAC)
      DIMENSION INSSEG(  MXFAC,MXATM)
      DIMENSION XYZSPA(3,MXFAC,MXATM)
      DIMENSION IJKSPA(  MXFAC,MXATM)
      DIMENSION NUMPPS(  MXFAC,MXATM)
      DIMENSION   APEX(3,MXAPX)
      DATA ZERO /0.0D+00/
C
      DBUG=.FALSE.
      MORE=.FALSE.
      MORE=MORE.OR.DBUG
      OUT =.FALSE. 
      OUT =OUT.OR.MORE
      SOME=.FALSE.
      SOME=SOME.OR.OUT
C
      DBUG=DBUG.AND.(IAP.EQ.1)
      MORE=MORE.AND.(IAP.EQ.1)
      OUT =OUT .AND.(IAP.EQ.1)
      SOME=SOME.AND.(IAP.EQ.1)
C
C     ----- APPROXIMATE SPHERE WITH SEGMENTS AND POINTS -----
C
      CALL HND_COSSPH(NSEG,NFAC,NDIV,
     1                IJKFAC,XYZSEG,IJKSEG,MXFAC,APEX,MXAPX,
     2                DSURF,DVOL,ADIAG)
C
C     ----- DEBUG PRINTING -----
C
      IF(OUT) THEN
         WRITE(IW,9999) NSEG,NFAC,NDIV,DSURF,DVOL
         WRITE(IW,9995) ADIAG
         IF(MORE) THEN
            WRITE(IW,9998)
            DO LSEG=1,NSEG
               WRITE(IW,9997) LSEG,XYZSEG(1,LSEG),XYZSEG(2,LSEG),
     1                             XYZSEG(3,LSEG),IJKSEG(  LSEG)
            ENDDO
         ENDIF
         IF(DBUG) THEN
            WRITE(IW,9996)
            DO LFAC=1,NFAC
               MFAC=LFAC+NSEG
               WRITE(IW,9997) MFAC,XYZSEG(1,MFAC),XYZSEG(2,MFAC),
     1                             XYZSEG(3,MFAC),IJKSEG(  MFAC)
            ENDDO
         ENDIF
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- SET MOLECULE -----
C
      DO IAT=1,NAT
         DO I=1,3
            XYZATM(I,IAT)=C(I,IAT)
         ENDDO
      ENDDO
      DO IAT=1,NAT
         IF(RADIUS(IAT).EQ.ZERO) THEN
            RATM(IAT)=ZERO
         ELSE
            RATM(IAT)=(RADIUS(IAT)+RSOLV)/BOHR
         ENDIF
      ENDDO
C
C     ----- CREATE -SOLVENT ACCESSIBLE SURFACE- OF THE MOLECULE -----
C
      CALL HND_COSSAS(NAT,XYZATM,RATM,MXATM,
     1                NSPA,NPPA,XYZSPA,IJKSPA,
     2                NSEG,NFAC,XYZSEG,IJKSEG,INSSEG,
     3                NUMPPS,MXFAC,rtdb)
C
      RETURN
 9999 FORMAT(' NSEG,NFAC,NDIV=NFAC/NSEG,DSURF,DVOL = ',3I7,2F10.6)
 9998 FORMAT('  SEG  ','      X     ','      Y     ','      Z     ',
     1       ' SEG ',/,1X,47(1H-))
 9997 FORMAT(I7,3F12.8,I5,F12.8)
 9996 FORMAT('  FAC  ','      X     ','      Y     ','      Z     ',
     1       ' SEG ',/,1X,47(1H-))
 9995 FORMAT(' ADIAG           = ',F12.6)
      END
      SUBROUTINE HND_COSSAS(NAT,XYZATM,RATM,MXATOM,
     1                  NSPA,NPPA,XYZSPA,IJKSPA,
     2                  NSEG,NFAC,XYZSEG,IJKSEG,INSSEG,
     3                  NUMPPS,MXFACE,rtdb)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
c
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
c
      PARAMETER  (MXEFC =3*7*7*7)
      PARAMETER  (BOHR=0.529177249D+00)
      CHARACTER*8 EFCLAB
      CHARACTER*8 ERRMSG
      CHARACTER*8 BLANK
      LOGICAL     SOME
      LOGICAL     OUT
      LOGICAL     DBUG
      LOGICAL     INSSEG
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_COSPRT/IPRT      
      COMMON/HND_COSOPT/ICOS,NCOS
      COMMON/HND_COSPAR/DIELEC,SCREEN,RSOLV
      COMMON/HND_COSDAT/ADIAG,DSURF,DVOL,SRFMOL,VOLMOL 
      COMMON/HND_COSNRG/ATMEFC,ELCEFC,EFCEFC,ALLEFC,SOLNRG,CAVDSP
      COMMON/HND_COSBEM/EFCS(MXEFC),IATEFC(MXEFC)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      DIMENSION XYZATM(3,MXATOM)
      DIMENSION   RATM(  MXATOM)
      DIMENSION   NSPA(  MXATOM)
      DIMENSION   NPPA(  MXATOM)
      DIMENSION XYZSEG(3,MXFACE)
      DIMENSION IJKSEG(  MXFACE)
      DIMENSION INSSEG(  MXFACE,MXATOM)
      DIMENSION XYZSPA(3,MXFACE,MXATOM)
      DIMENSION IJKSPA(  MXFACE,MXATOM)
      DIMENSION NUMPPS(  MXFACE,MXATOM)
      DIMENSION ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-COSSAS-'/
      DATA BLANK   /'        '/
      DATA ZERO    /0.0D+00/
C
      DIST(XI,YI,ZI,XJ,YJ,ZJ)=SQRT((XJ-XI)**2+(YJ-YI)**2+(ZJ-ZI)**2)
C
      DBUG=.FALSE.
      OUT =.FALSE.  
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
      SOME=SOME.OR.IPRT.NE.0
C
      DBUG=DBUG.AND.(IAP.EQ.1)
      OUT =OUT .AND.(IAP.EQ.1)
      SOME=SOME.AND.(IAP.EQ.1)
C
      IF(SOME) THEN
         WRITE(IW,9999)
      ENDIF
C
C     ----- PRINT ATOMIC CENTERS -----
C
      IF(OUT) THEN
         WRITE(IW,9998)
         DO IAT=1,NAT
            WRITE(IW,9997) IAT,XYZATM(1,IAT),XYZATM(2,IAT),
     1                                       XYZATM(3,IAT),
     2                    (RATM(IAT)*BOHR-RSOLV)
         ENDDO
      ENDIF
C
C     ----- CLEAR ARRAYS ..... -----
C
      DO IAT=1,NAT
         DO I=1,MXFACE
            IJKSPA(I,IAT)=0
            NUMPPS(I,IAT)=0
         ENDDO
      ENDDO
C
C     ----- SIFT THROUGH ATOMIC CENTERS AND DECIDE IF A FACE -----
C           BELONGS TO THE -SAS- OR IS INSIDE THE MOLECULE.
C
      DO IAT=1,NAT
C
         IF(RATM(IAT).NE.ZERO) THEN
C
            DO ISEG=1,NSEG
               IJKSPA(       ISEG,IAT)=IJKSEG(       ISEG)
               DO M=1,3
                  XYZSPA(M,     ISEG,IAT)=XYZSEG(M,     ISEG)*RATM(IAT)
     1                                   +XYZATM(M,IAT)
               ENDDO
            ENDDO
            DO IFAC=1,NFAC
               IJKSPA(  IFAC+NSEG,IAT)=IJKSEG(  IFAC+NSEG)
               DO M=1,3
                  XYZSPA(M,IFAC+NSEG,IAT)=XYZSEG(M,IFAC+NSEG)*RATM(IAT)
     1                                   +XYZATM(M,IAT)
               ENDDO
            ENDDO
            IF(DBUG) THEN
               WRITE(IW,9996) IAT
               WRITE(IW,9995) (IJKSPA(IFAC+NSEG,IAT),IFAC=1,NFAC)
            ENDIF
            DO JAT=1,NAT
               DIJ=DIST(XYZATM(1,IAT),
     1                  XYZATM(2,IAT),
     2                  XYZATM(3,IAT),
     3                  XYZATM(1,JAT),
     4                  XYZATM(2,JAT),
     5                  XYZATM(3,JAT))
               IF(JAT.NE.IAT.AND.(DIJ.LT.(RATM(IAT)+RATM(JAT)))) THEN
                  DO IFAC=1,NFAC
                     DUM=DIST(XYZSPA(1,IFAC+NSEG,IAT),
     1                        XYZSPA(2,IFAC+NSEG,IAT),
     2                        XYZSPA(3,IFAC+NSEG,IAT),
     3                        XYZATM(1,JAT),
     4                        XYZATM(2,JAT),
     5                        XYZATM(3,JAT))
                     IF(DUM.LT.RATM(JAT)) THEN
                        IJKSPA(IFAC+NSEG,IAT)=0
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            IF(DBUG) THEN
               WRITE(IW,9996) IAT
               WRITE(IW,9995) (IJKSPA(IFAC+NSEG,IAT),IFAC=1,NFAC)
            ENDIF
C
C     ----- CHECK SEGMENTS BELONGING TO -SAS- -----
C
            DO IFAC=1,NSEG+NFAC
               INSSEG(IFAC,IAT)=.TRUE.
            ENDDO
            DO IFAC=1,NFAC
               ISEG=IJKSPA(IFAC+NSEG,IAT)
               IF(ISEG.NE.0) THEN
                  INSSEG(IFAC+NSEG,IAT)=.FALSE.
                  INSSEG(     ISEG,IAT)=.FALSE.
               ENDIF
            ENDDO
            IF(DBUG) THEN
               WRITE(IW,9994) IAT
               WRITE(IW,9993) (INSSEG(IFAC,IAT),IFAC=1,NSEG+NFAC)
            ENDIF
            MSEG=0
            DO ISEG=1,NSEG
               IF(.NOT.INSSEG(     ISEG,IAT)) THEN
                  MSEG=MSEG+1
               ENDIF
            ENDDO
            MFAC=0
            DO IFAC=1,NFAC
               IF(.NOT.INSSEG(IFAC+NSEG,IAT)) THEN
                  MFAC=MFAC+1
               ENDIF
            ENDDO
            NSPA(IAT)=MSEG
            NPPA(IAT)=MFAC
C
C     ----- SURFACE AREA OF SEGMENTS -----
C
            DO ISEG=1,NSEG
               NUMPPS(ISEG,IAT)=0
            ENDDO
            DO IFAC=1,NFAC
               ISEG=IJKSPA(IFAC+NSEG,IAT)
               IF(ISEG.NE.0) THEN
                  NUMPPS(ISEG,IAT)=NUMPPS(ISEG,IAT)+1
               ENDIF
            ENDDO
C
         ENDIF
C
      ENDDO
C
      IF(SOME) THEN
         WRITE(IW,9985) NSEG,NFAC
         WRITE(IW,9992)
         DO IAT=1,NAT
            NPP=0
            DO ISEG=1,NSEG
               NPP=NPP+NUMPPS(ISEG,IAT)
            ENDDO
            WRITE(IW,9991) IAT,NSPA(IAT),NPPA(IAT),NPP
         ENDDO
      ENDIF
      IF(DBUG) THEN
         WRITE(IW,9987)
         DO IAT=1,NAT
            DO ISEG=1,NSEG
               WRITE(IW,9986) IAT,ISEG,NUMPPS(ISEG,IAT)
            ENDDO
         ENDDO
      ENDIF
C
C     ----- SAVE COORDINATES OF SURFACE POINTS -----
C           SAVE SEGMENT SURFACES
C           SAVE SEGMENT TO ATOM MAPPING
C
      SRFMOL=ZERO
      VOLMOL=ZERO
      IEF   =0
      DO IAT=1,NAT
         IF(RATM(IAT).NE.ZERO) THEN
            DO ISEG=1,NSEG
               IF(.NOT.INSSEG(ISEG,IAT)) THEN
                  IEF=IEF+1
                  IF(IEF.GT.MXEFC) THEN
                     WRITE(IW,*) 'TOO MANY -COSMO- SURFACE POINTS'
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  DO I=1,3
                     EFCC(I,IEF)=XYZSEG(I,ISEG)*RATM(IAT)
     1                          +XYZATM(I,IAT )
                  ENDDO
                  IPP=NUMPPS(ISEG,IAT)
                  EFCS(IEF)=         DBLE(IPP)*DSURF*RATM(IAT)**2
                  SRFMOL   =SRFMOL + DBLE(IPP)*DSURF*RATM(IAT)**2
                  VOLMOL   =VOLMOL + DBLE(IPP)*DVOL *RATM(IAT)**3
                  IATEFC(IEF)=IAT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      SRFMOL=SRFMOL*(BOHR**2)
      VOLMOL=VOLMOL*(BOHR**3)
C
      NEFC=IEF
      NCOS=IEF
      IF(SOME) THEN
         WRITE(IW,9990) NEFC
         WRITE(IW,9984) SRFMOL 
         WRITE(IW,9983) VOLMOL 
      ENDIF 
C
C     ----- Cavity/Dispersion FREE ENERGY ---
C           Sitkoff, Sharp, and Honig,
C           J.Phys.Chem. 98, 1978 (1994)
C
      CAVDSP=0.860+0.005*SRFMOL
      IF(SOME) THEN
         WRITE(IW,9981) CAVDSP
      ENDIF
C
C     ----- PRINT SEGMENT SURFACES -----
C
      IF(DBUG) THEN
         WRITE(IW,9989)
         DO IEF=1,NEFC
            WRITE(IW,9988) IEF,EFCS(IEF),IATEFC(IEF)
         ENDDO
      ENDIF
C
      DO IEF=1,NEFC
         EFCZ(IEF)=ZERO
      ENDDO
      DO IEF=1,NEFC
         EFCLAB(IEF)=BLANK
      ENDDO
c
c     ----- write out to -rtdb- -----
c
      if(.not.rtdb_put(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     $   call errquit('hnd_cossas: rtdb put failed for nefc  ',911,
     &       RTDB_ERR)
      if(.not.rtdb_put(rtdb,'cosmo:efcc',mt_dbl,3*nefc,efcc))
     $   call errquit('hnd_cossas: rtdb put failed for efcc  ',912,
     &       RTDB_ERR)
      if(.not.rtdb_put(rtdb,'cosmo:efcz',mt_dbl,  nefc,efcz))
     $   call errquit('hnd_cossas: rtdb put failed for efcz  ',913,
     &       RTDB_ERR)
      if(.not.rtdb_put(rtdb,'cosmo:efcs',mt_dbl,  nefc,efcs))
     $   call errquit('hnd_cossas: rtdb put failed for efcs  ',914,
     &       RTDB_ERR)
      if(.not.rtdb_put(rtdb,'cosmo:iatefc',mt_int,nefc,iatefc))
     $   call errquit('hnd_cossas: rtdb put failed for iatefc',915,
     &       RTDB_ERR)
      if(.not.rtdb_cput(rtdb,'character variable',nefc,efclab))
     $   call errquit('hnd_cossas: rtdb put failed for efclab',916,
     &       RTDB_ERR)
C
      RETURN
 9999 FORMAT(/,1X,'solvent accessible surface',/,1X,26(1H-))
 9998 FORMAT(/,1X,'---------- ATOMIC COORDINATES (A.U.) ----------',
     1            '-- VDWR(ANG.) --')
 9997 FORMAT(  1X,I5,3F14.8,F10.3)
 9996 FORMAT(/,1X,'---------- SEGMENTS FOR -IAT- = ',I5)
 9995 FORMAT(16I4)
 9994 FORMAT(/,1X,'-INSSEG- FACES FOR IAT = ',I5)
 9993 FORMAT(16L4)
 9992 FORMAT(  1X,'atom',' ( ','  nspa',',','  nppa',' )',/,1X,22(1H-))
 9991 FORMAT(  1X,   I4 ,' ( ',     I6 ,',',     I6 ,' )',I8)
 9990 FORMAT(  1X,'number of -cosmo- surface points = ',I8)
 9989 FORMAT(  1X,'SEGMENT SURFACES =',/,1X,18(1H-))
 9988 FORMAT(I8,F10.5,I5)
 9987 FORMAT(  1X,'NUMBER OF FACES / SEGMENT =',/,1X,27(1H-))
 9986 FORMAT(3I5)
 9985 FORMAT(' number of segments per atom = ',I10,/,
     1       ' number of   points per atom = ',I10)
 9984 FORMAT(' molecular surface = ',F10.3,' angstrom**2')
 9983 FORMAT(' molecular volume  = ',F10.3,' angstrom**3')
 9981 FORMAT(' G(cav/disp)       = ',F10.3,' kcal/mol')
      END
      SUBROUTINE HND_COSSPH(NSEG,NFAC,NDIV,
     1                  IJKFAC,XYZSEG,IJKSEG,MXFACE,APEX,MXAPEX,
     2                  DSURF,DVOL,ADIAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C              ----- STARTING FROM -ICOSAHEDRON- -----
C
C     PASS, NAPEX, NFACE, ERROR =   0      12      20      20
C     PASS, NAPEX, NFACE, ERROR =   1      42      80     100    0.4982
C     PASS  NAPEX, NFACE, ERROR =   2     162     320     420    0.1848
C     PASS  NAPEX, NFACE, ERROR =   3     642    1280    1700    0.0523
C     PASS  NAPEX, NFACE, ERROR =   4    2562    5120    6820    0.0135
C     PASS  NAPEX, NFACE, ERROR =   5   10242   20480   27300    0.0034
C
C              ----- STARTING FROM -OCTAHEDRON-  -----
C
C     PASS, NAPEX, NFACE, ERROR =   0       6       8       8
C     PASS, NAPEX, NFACE, ERROR =   1      18      32      40    0.8075
C     PASS  NAPEX, NFACE, ERROR =   2      66     128     168    0.4557
C     PASS  NAPEX, NFACE, ERROR =   3     258     512     680    0.1619
C     PASS  NAPEX, NFACE, ERROR =   4    1026    2048    2728    0.0451
C     PASS  NAPEX, NFACE, ERROR =   5    4098    8192   10920    0.0116
C     PASS  NAPEX, NFACE, ERROR =   6   16386   32768   43688    0.0029
C
      PARAMETER (MXPASS=    7)
      LOGICAL ICOS
      LOGICAL OCTA
      LOGICAL SOME,OUT,DBUG
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_COSMOD/LINEQ,MINBEM,MAXBEM,IFICOS,IFSCRN
      DIMENSION   APEX(3,*)
      DIMENSION IJKFAC(3,*)
      DIMENSION IJKSEG(  *)
      DIMENSION XYZSEG(3,*)
      DIMENSION MINFAC(MXPASS)
      DIMENSION MAXFAC(MXPASS)
      DIMENSION MINICO(MXPASS)
      DIMENSION MAXICO(MXPASS)
      DIMENSION MINOCT(MXPASS)
      DIMENSION MAXOCT(MXPASS)
      DIMENSION IJKNEW(3)
      DIMENSION IJKOLD(3)
      EQUIVALENCE (IJKOLD(1),IOLD),(IJKOLD(2),JOLD),(IJKOLD(3),KOLD)
      EQUIVALENCE (IJKNEW(1),INEW),(IJKNEW(2),JNEW),(IJKNEW(3),KNEW)
      DATA MINICO /    1,   21,  101,  421, 1701, 6821,    0/
      DATA MAXICO /   20,  100,  420, 1700, 6820,27300,    0/
      DATA MINOCT /    1,    9,   41,  169,  681, 2729,10921/
      DATA MAXOCT /    8,   40,  168,  680, 2728,10920,43688/
      DATA ZERO  /0.0D+00/
      DATA ONE   /1.0D+00/
      DATA TWO   /2.0D+00/
      DATA THREE /3.0D+00/
      DATA FOUR  /4.0D+00/
C
      DIST(XI,YI,ZI,XJ,YJ,ZJ)=SQRT((XJ-XI)**2+(YJ-YI)**2+(ZJ-ZI)**2)
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
C
      DBUG=DBUG.AND.(IAP.EQ.1)
      OUT =OUT .AND.(IAP.EQ.1)
      SOME=SOME.AND.(IAP.EQ.1)
C
      PI=FOUR*ATAN(ONE)
      RAD=ONE
      CIR= TWO*PI*RAD
      SRF=FOUR*PI*RAD**2
      VOL=FOUR*PI*RAD**3/THREE
C
      NPASS=MAXBEM
C
C     ----- DEFINE  HEDRON  -----
C           DEFINE -MINFAC- 
C           DEFINE -MAXFAC- 
C
      ICOS=IFICOS.NE.0
      OCTA=.NOT.ICOS
      IF(ICOS) THEN
         CALL HND_SPHICO(APEX,NAPEX,IJKFAC,IJKSEG,NFACE)
         DO IPASS=1,MXPASS
            MINFAC(IPASS)=MINICO(IPASS)
            MAXFAC(IPASS)=MAXICO(IPASS)
         ENDDO
      ENDIF
      IF(OCTA) THEN
         CALL HND_SPHOCT(APEX,NAPEX,IJKFAC,IJKSEG,NFACE)
         DO IPASS=1,MXPASS
            MINFAC(IPASS)=MINOCT(IPASS)
            MAXFAC(IPASS)=MAXOCT(IPASS)
         ENDDO
      ENDIF
      IF(SOME.OR.OUT.OR.DBUG) THEN
         IF(ICOS) THEN
            WRITE(IW,9994)
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
         IF(OCTA) THEN
            WRITE(IW,9982)
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
         IF(OUT) THEN
            WRITE(IW,*) '-MINBEM- = ',MINBEM
            WRITE(IW,*) '-MAXBEM- = ',MAXBEM
            WRITE(IW,*) '-MINFAC- = ',MINFAC
            WRITE(IW,*) '-MAXFAC- = ',MAXFAC
            WRITE(IW,*) '-NPASS - = ',NPASS
            WRITE(IW,9999)
            DO IAPEX=1,NAPEX
               WRITE(IW,9998) IAPEX,APEX(1,IAPEX),
     1                              APEX(2,IAPEX),
     2                              APEX(3,IAPEX)
            ENDDO
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
      ENDIF
C
C     ----- LOOP OVER DIVISIONS TO CREATE SPHERE -----
C
      MXFAC=0
      IPASS=1
  100 IPASS=IPASS+1
         MNFAC=MXFAC+1
         MXFAC=NFACE
         IF(OUT) THEN
            WRITE(IW,9996) IPASS,NAPEX,NFACE,MNFAC,MXFAC
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
C
         DMIN =ONE
         MAPEX=NAPEX
         MFACE=NFACE
         DO LFACE=MNFAC,MXFAC
            IOLD=IJKFAC(1,LFACE)
            JOLD=IJKFAC(2,LFACE)
            KOLD=IJKFAC(3,LFACE)
            CALL HND_SPHAPX(APEX,MAPEX,IJKFAC,IJKSEG,MFACE,LFACE,
     1                      IJKOLD,IJKNEW,DIJK)
            DMIN=MIN(DMIN,DIJK)
         ENDDO
         NAPEX=MAPEX
         NFACE=MFACE
         IF(OUT) THEN
            WRITE(IW,9995) NAPEX,NFACE
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
C
C     ----- PRINT OUT NEW APECES -----
C
         IF(DBUG) THEN
            DO IAPEX=1,NAPEX
               WRITE(IW,9998) IAPEX,APEX(1,IAPEX),APEX(2,IAPEX),
     1                              APEX(3,IAPEX)
            ENDDO
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
C
C     ----- PRINT APPROXIMATE VOLUME -----
C
         RADAPP=    DMIN
         RADRAT=    DMIN
         RADERR=ONE-RADRAT
         SRFAPP=SRF*DMIN**2
         SRFRAT=    DMIN**2
         SRFERR=ONE-SRFRAT
         VOLAPP=VOL*DMIN**3
         VOLRAT=    DMIN**3
         VOLERR=ONE-VOLRAT
         IF(OUT) THEN
            WRITE(IW,9997) VOL,VOLAPP,VOLRAT,VOLERR
            WRITE(IW,9992) SRF,SRFAPP,SRFRAT,SRFERR
            WRITE(IW,9991) RAD,RADAPP,RADRAT,RADERR
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
C
C     ----- ASSIGN EARLY SEGMENT TO EACH FACE -----
C
         IF(IPASS.GT.(MINBEM+1)) THEN
            IF(DBUG) THEN
               WRITE(IW,9981) IPASS
               WRITE(IW,9980) (MINFAC(I),I=1,IPASS)
               WRITE(IW,9979) (MAXFAC(I),I=1,IPASS)
CBERT          CALL HND_FIFLSH(IW)
            ENDIF
            MAXSEG=MAXFAC(MINBEM)
            LFACMN=MINFAC(IPASS)
            LFACMX=MAXFAC(IPASS)
            IF(DBUG) THEN
               WRITE(IW,9990) IPASS
               WRITE(IW,9988) (IJKSEG(LFACE),LFACE=LFACMN,LFACMX)
CBERT          CALL HND_FIFLSH(IW)
            ENDIF
            DO LFACE=LFACMN,LFACMX
               IJKSEG(LFACE)=IJKSEG(IJKSEG(LFACE))
               IF(IJKSEG(LFACE).GT.MAXSEG) THEN
                  WRITE(IW,9987) LFACE,IJKSEG(LFACE)
               ENDIF
            ENDDO
            IF(DBUG) THEN
               WRITE(IW,9989) IPASS
               WRITE(IW,9988) (IJKSEG(LFACE),LFACE=LFACMN,LFACMX)
CBERT          CALL HND_FIFLSH(IW)
            ENDIF
         ENDIF
C
      IF(IPASS.LT.NPASS) GO TO 100
C
C     ----- END OF LOOP OVER TESSALATING PASSES -----
C
      IF(DBUG) THEN
         DO IPASS=1,NPASS
            LFACMN=MINFAC(IPASS)
            LFACMX=MAXFAC(IPASS)
            WRITE(IW,9989) IPASS
            WRITE(IW,*) '-LFACMN- = ',LFACMN
            WRITE(IW,*) '-LFACMX- = ',LFACMX
            WRITE(IW,9988) (IJKSEG(LFACE),LFACE=LFACMN,LFACMX)
CBERT       CALL HND_FIFLSH(IW)
         ENDDO
      ENDIF
      IF(SOME.OR.OUT.OR.DBUG) THEN
         WRITE(IW,9993) NPASS,NAPEX,MINFAC(NPASS),MAXFAC(NPASS),
     1                  RADAPP,RADERR,SRFAPP,SRFERR,VOLAPP,VOLERR
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- AT THIS POINT EACH OF THE FACES IS ASSIGNED TO ONE -----
C           SEGMENT. NOW DEFINE CENTERS OF SEGMENTS ...
C
      THIRD =ONE/THREE
      LFACMN= MINFAC(MINBEM)
      LFACMX= MAXFAC(MINBEM)
      DO LFACE=LFACMN,LFACMX
         MFACE=LFACE-LFACMN+1
         IJKSEG(MFACE)=MFACE
         I=IJKFAC(1,LFACE)
         J=IJKFAC(2,LFACE)
         K=IJKFAC(3,LFACE)
         DO M=1,3
            XYZSEG(M,MFACE)=(APEX(M,I)+APEX(M,J)+APEX(M,K))*THIRD
         ENDDO
         DSEG=ONE/DIST(XYZSEG(1,MFACE),XYZSEG(2,MFACE),XYZSEG(3,MFACE),
     1                 ZERO,ZERO,ZERO)
         DO M=1,3
            XYZSEG(M,MFACE)=XYZSEG(M,MFACE)*DSEG
         ENDDO
      ENDDO
      NSEG=(LFACMX-LFACMN+1)
C
      IF(DBUG) THEN
         LFACMN=1
         LFACMX=NSEG
         WRITE(IW,*)    'SEGMENT TO SEGMENT MAPPING = '
         WRITE(IW,9988) (IJKSEG(LFACE),LFACE=LFACMN,LFACMX)
      ENDIF
C
C     ----- NOW THE FACES ... -----
C
      IF(NPASS.GT.MINBEM) THEN
         LFACMN=MINFAC(MINBEM+1)
         LFACMX=MAXFAC(NPASS   )
         DO LFACE=LFACMN,LFACMX
            MFACE=LFACE-LFACMN+1    
     1                        +(MAXFAC(MINBEM)-MINFAC(MINBEM)+1)
            IJKSEG(MFACE)=IJKSEG(LFACE)
     1                        -(               MINFAC(MINBEM)-1)
            I=IJKFAC(1,LFACE)
            J=IJKFAC(2,LFACE)
            K=IJKFAC(3,LFACE)
            DO M=1,3
               XYZSEG(M,MFACE)=(APEX(M,I)+APEX(M,J)+APEX(M,K))*THIRD
            ENDDO
            DSEG=ONE/DIST(XYZSEG(1,MFACE),
     1                    XYZSEG(2,MFACE),
     2                    XYZSEG(3,MFACE),ZERO,ZERO,ZERO)
            DO M=1,3
               XYZSEG(M,MFACE)=XYZSEG(M,MFACE)*DSEG
            ENDDO
         ENDDO
         NFAC=(LFACMX-LFACMN+1)
      ELSE
         DO ISEG=1,NSEG
            IFAC=ISEG+NSEG
            IJKSEG(IFAC)=IJKSEG(ISEG)
            DO M=1,3
               XYZSEG(M,IFAC)=XYZSEG(M,ISEG)
            ENDDO
         ENDDO
         NFAC=NSEG
      ENDIF
C
      IF(DBUG) THEN
         LFACMN=NSEG+1
         LFACMX=NSEG+NFAC
         WRITE(IW,*)    ' FACE   TO SEGMENT MAPPING = '
         WRITE(IW,9988) (IJKSEG(LFACE),LFACE=LFACMN,LFACMX)
      ENDIF
C
C     ----- CALCULATE -DSURF DVOL- FOR THE -COSMO- THEORY -----
C
      NDIV =NFAC/NSEG
      DSURF=SRF/DBLE(NFAC)
      DVOL =VOL/DBLE(NFAC)
      IF(SOME.OR.OUT.OR.DBUG) THEN
         WRITE(IW,9986) NSEG,NFAC,NDIV,DSURF,DVOL
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
      IF(OUT) THEN
         WRITE(IW,9985)
         DO I=1,NSEG
            DONE=DIST(XYZSEG(1,I),XYZSEG(2,I),XYZSEG(3,I),
     1                ZERO,ZERO,ZERO)
            WRITE(IW,9984) I,
     1                     XYZSEG(1,I),XYZSEG(2,I),XYZSEG(3,I),
     2                     IJKSEG(I),DONE
         ENDDO
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
      IF(DBUG) THEN
         WRITE(IW,9985)
         DO I=NSEG+1,NSEG+NFAC
            DONE=DIST(XYZSEG(1,I),XYZSEG(2,I),XYZSEG(3,I),
     1                ZERO,ZERO,ZERO)
            WRITE(IW,9984) (I-NSEG),
     1                     XYZSEG(1,I),XYZSEG(2,I),XYZSEG(3,I),
     2                     IJKSEG(I),DONE
         ENDDO
      ENDIF
C
C     ----- CALCULATE -ADIAG- OF THE -COSMO- THEORY -----
C
      AVGDIA=ZERO
      AVGFAC=ZERO
      DO MSEG=1,NSEG
         SUM=ZERO
         DO LSEG=1,NSEG
            IF(LSEG.NE.MSEG) THEN
               L1=MSEG
               L2=LSEG
         SUM=SUM+RAD/DIST(XYZSEG(1,L2),XYZSEG(2,L2),XYZSEG(3,L2),
     1                    XYZSEG(1,L1),XYZSEG(2,L1),XYZSEG(3,L1))
            ENDIF
         ENDDO
         FAC=(DBLE(NSEG)-SUM)/SQRT(DBLE(NSEG))
         ADIAG=SQRT(FOUR*PI)*FAC
         IF(SOME.OR.OUT.OR.DBUG) THEN
            WRITE(IW,9983) MSEG,ADIAG,FAC,DBLE(NSEG),SUM
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
         AVGDIA=AVGDIA+ADIAG
         AVGFAC=AVGFAC+FAC
      ENDDO
      ADIAG=AVGDIA/DBLE(NSEG)
      FAC  =AVGFAC/DBLE(NSEG)
      IF(SOME.OR.OUT.OR.DBUG) THEN
         WRITE(IW,9978)      ADIAG,FAC               
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
      RETURN
 9999 FORMAT(/,1X,'  APEX',5X,'X',6X,5X,'Y',6X,5X,'Z',6X,/,1X,42(1H-))
 9998 FORMAT(1X,I6,F12.8,F12.8,F12.8)
 9997 FORMAT(' VOL, APPROX., RATIO, ERROR = ',2F12.8,2 F8.4)
 9996 FORMAT(' PASS, NAPEX, NFACE, MNFAC, MXFAC = ',I3,4I8)
 9995 FORMAT('       NAPEX, NFACE               = ',3X,2I8)
 9994 FORMAT(1X,'SPHERE FROM -ICOSAHEDRON-',/,1X,25(1H-))
 9993 FORMAT(' NPASS = ',I2,' NAPEX = ',I8,
     1       ' MINFAC = ',I8,' MAXFAC = ',I8,/,
     2       ' RAD = ',F10.6,' ERROR = ',F8.4,/,
     3       ' SRF = ',F10.6,' ERROR = ',F8.4,/,
     4       ' VOL = ',F10.6,' ERROR = ',F8.4)
 9992 FORMAT(' SRF, APPROX., RATIO, ERROR = ',2F12.8,2 F8.4)
 9991 FORMAT(' RAD, APPROX., RATIO, ERROR = ',2F12.8,2 F8.4)
 9990 FORMAT(' ABSOLUTE -IJKSEG- , FOR -IPASS- = ',I5)
 9989 FORMAT(' RELATIVE -IJKSEG- , FOR -IPASS- = ',I5)
 9988 FORMAT(12I6)
 9987 FORMAT(' ASSIGNED SEGMENT FOR -LFACE- = ',I7,
     1       ' IS = ',I7,' ( GREATER THAN -MAXSEG- = ',I4,' )')
 9986 FORMAT(' NSEG,NFAC,NDIV=NFAC/NSEG,DSURF,DVOL = ',3I7,2F10.6)
 9985 FORMAT('   PT  ','      X     ','      Y     ','      Z     ',
     1       ' SEG ','    NORM    ',/,1X,59(1H-))
 9984 FORMAT(I7,3F12.8,I5,F12.8)
 9983 FORMAT(' MSEG,ADIAG,FAC,M,SUM = ',I7,4F12.6)
 9982 FORMAT(1X,'SPHERE FROM -OCTAHEDRON-',/,1X,24(1H-))
 9981 FORMAT(' PASS # = ',I5)
 9980 FORMAT(' MINFAC = ',10I5)
 9979 FORMAT(' MAXFAC = ',10I5)
 9978 FORMAT('      ADIAG,FAC       = ',   2F12.6)
      END
      SUBROUTINE HND_SPHOCT(APEX,NAPEX,IJKFAC,IJKSEG,NFACE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION    XYZ(3,6)
      DIMENSION    IJK(3,8)
      DIMENSION   APEX(3,*)
      DIMENSION IJKFAC(3,*)
      DIMENSION IJKSEG(  *)
      DATA XYZ / 1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00,
     1          -1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00,
     2           0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00/
      DATA IJK / 5, 1, 2, 5, 2, 3, 5, 3, 4, 5, 4, 1,
     1           6, 1, 2, 6, 2, 3, 6, 3, 4, 6, 4, 1/
C
      OUT=.FALSE.
      OUT=OUT.AND.(IAP.EQ.1)
C
      IF(OUT) THEN
         WRITE(IW,9997)
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- SET THE 6 APECES OF AN OCTAHEDRON -----
C
C     1     1.     0.     0.
C     2     0.     1.     0.
C     3    -1.     0.     0.
C     4     0.    -1.     0.
C     5     0.     0.     1.
C     6     0.     0.    -1.
C
      NAPEX=6
      DO IAPEX=1,NAPEX
         APEX(1,IAPEX)=XYZ(1,IAPEX)
         APEX(2,IAPEX)=XYZ(2,IAPEX)
         APEX(3,IAPEX)=XYZ(3,IAPEX)
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9999)
         DO IAPEX=1,NAPEX
            WRITE(IW,9998) IAPEX,APEX(1,IAPEX),APEX(2,IAPEX),
     1                           APEX(3,IAPEX)
         ENDDO
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
      NFACE=8
      DO LFACE=1,NFACE
         IJKFAC(1,LFACE)=IJK(1,LFACE)
         IJKFAC(2,LFACE)=IJK(2,LFACE)
         IJKFAC(3,LFACE)=IJK(3,LFACE)
         IJKSEG(  LFACE)=      LFACE
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,*) '...... END OF -SPHOCT- ......'
      ENDIF
      RETURN
 9999 FORMAT(/,1X,'  APEX',5X,'X',6X,5X,'Y',6X,5X,'Z',6X,/,1X,42(1H-))
 9998 FORMAT(1X,I6,F12.8,F12.8,F12.8)
 9997 FORMAT(/,1X,'OCTAHEDRON',/,1X,10(1H-))
      END
      SUBROUTINE HND_SPHICO(APEX,NAPEX,IJKFAC,IJKSEG,NFACE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION      C(3,12)
      DIMENSION      S(3,12)
      DIMENSION    IJK(3,20)
      DIMENSION   APEX(3,*)
      DIMENSION IJKFAC(3,*)
      DIMENSION IJKSEG(  *)
      DATA C   / 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00,
     1           0.0D+00,-1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00,
     2           0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, 1.0D+00,
     3           0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00,
     4           1.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00,
     5          -1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00, 0.0D+00/
      DATA S   / 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00,
     1           0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00,
     2           1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00, 0.0D+00,
     3           1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00, 0.0D+00,
     4           0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00,
     5           0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00,-1.0D+00, 0.0D+00/
      DATA IJK / 1, 2, 9, 1, 9, 5, 1, 5, 6, 1, 6,11, 1,11, 2,
     1                    2, 9, 7, 2, 7, 8, 2, 8,11,
     2           3, 4,10, 3,10, 5, 3, 5, 6, 3, 6,12, 3,12, 4,
     3                    4,10, 7, 4, 7, 8, 4, 8,12,
     4           9,10, 7, 9, 5,10,11,12, 8,11, 6,12/
      DATA ONE  /1.0D+00/
      DATA TWO  /2.0D+00/
      DATA FIVE /5.0D+00/
C
      OUT=.FALSE.
      OUT=OUT.AND.(IAP.EQ.1)
C
      IF(OUT) THEN
         WRITE(IW,9997)
      ENDIF
C
C     ----- SET THE 12 APECES OF AN ICOSAHEDRON -----
C
C     1     0.     COSA   SINA
C     2     0.     COSA  -SINA
C     3     0.    -COSA   SINA
C     4     0.    -COSA  -SINA
C     5     SINA   0.     COSA
C     6    -SINA   0.     COSA
C     7     SINA   0.    -COSA
C     8    -SINA   0.    -COSA
C     9     COSA   SINA   0.
C    10     COSA  -SINA   0.
C    11    -COSA   SINA   0.
C    12    -COSA  -SINA   0.
C
      ANG=ACOS(ONE/SQRT(FIVE))/TWO
      COSA=COS(ANG)
      SINA=SIN(ANG)
      NAPEX=12
      DO IAPEX=1,NAPEX
         APEX(1,IAPEX)=COSA*C(1,IAPEX)+SINA*S(1,IAPEX)
         APEX(2,IAPEX)=COSA*C(2,IAPEX)+SINA*S(2,IAPEX)
         APEX(3,IAPEX)=COSA*C(3,IAPEX)+SINA*S(3,IAPEX)
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9999)
         DO IAPEX=1,NAPEX
            WRITE(IW,9998) IAPEX,APEX(1,IAPEX),APEX(2,IAPEX),
     1                           APEX(3,IAPEX)
         ENDDO
      ENDIF
C
      NFACE=20
      DO LFACE=1,NFACE
         IJKFAC(1,LFACE)=IJK(1,LFACE)
         IJKFAC(2,LFACE)=IJK(2,LFACE)
         IJKFAC(3,LFACE)=IJK(3,LFACE)
         IJKSEG(  LFACE)=      LFACE
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,*) '...... END OF -SPHICO- ......'
      ENDIF
      RETURN
 9999 FORMAT(/,1X,'  APEX',5X,'X',6X,5X,'Y',6X,5X,'Z',6X,/,1X,42(1H-))
 9998 FORMAT(1X,I6,F12.8,F12.8,F12.8)
 9997 FORMAT(/,1X,'ICOSAHEDRON',/,1X,11(1H-))
      END
      SUBROUTINE HND_SPHAPX(APEX,MAPEX,IJKFAC,IJKSEG,MFACE,LFACE,
     1                             IJKOLD,IJKNEW,DMIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT
      LOGICAL DUPLIC
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION   APEX(3,*)
      DIMENSION IJKFAC(3,*)
      DIMENSION IJKSEG(  *)
      DIMENSION IJKOLD(3)
      DIMENSION IJKNEW(3)
      DIMENSION    XYZ(3,3)
      DIMENSION      D(3)
      DIMENSION XYZIJK(3)
      EQUIVALENCE (XYZ(1,1),XIJ),(XYZ(2,1),YIJ),(XYZ(3,1),ZIJ),
     1            (XYZ(1,2),XJK),(XYZ(2,2),YJK),(XYZ(3,2),ZJK),
     2            (XYZ(1,3),XKI),(XYZ(2,3),YKI),(XYZ(3,3),ZKI)
      DATA ZERO  /0.0D+00/
      DATA PT5   /0.5D+00/
      DATA ONE   /1.0D+00/
      DATA THREE /3.0D+00/
      DATA TOL   /1.0D-04/
C
      DIST(X1,Y1,Z1,X2,Y2,Z2)=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
C
      OUT=.FALSE.
      OUT=OUT.AND.(IAP.EQ.1)
C
C     ----- CREATE MID-POINT OF THE 3 EDGES -----
C
      IOLD=IJKOLD(1)
      JOLD=IJKOLD(2)
      KOLD=IJKOLD(3)
      DO M=1,3
         XYZ(M,1)=(APEX(M,IOLD)+APEX(M,JOLD))*PT5
         XYZ(M,2)=(APEX(M,JOLD)+APEX(M,KOLD))*PT5
         XYZ(M,3)=(APEX(M,KOLD)+APEX(M,IOLD))*PT5
      ENDDO
C
C     ----- PROJECT ONTO SPHERE -----
C
      D(1)=DIST(XIJ,YIJ,ZIJ,ZERO,ZERO,ZERO)
      D(2)=DIST(XJK,YJK,ZJK,ZERO,ZERO,ZERO)
      D(3)=DIST(XKI,YKI,ZKI,ZERO,ZERO,ZERO)
      D(1)=ONE/D(1)
      D(2)=ONE/D(2)
      D(3)=ONE/D(3)
      DO L=1,3
         DO M=1,3
            XYZ(M,L)=XYZ(M,L)*D(L)
         ENDDO
      ENDDO
C
C     ----- CHECK FOR DUPLICATE APECES -----
C
      NEWAPX=0
      DO IAPX=1,3
         DUPLIC=.FALSE.
         DO LAPEX=1,MAPEX
            DD=DIST(XYZ(1,  IAPX),XYZ(2,  IAPX),XYZ(3,  IAPX),
     1              APEX(1,LAPEX),APEX(2,LAPEX),APEX(3,LAPEX))
            IF(ABS(DD).LT.TOL) THEN
               DUPLIC=.TRUE.
               LDUPLC=LAPEX
            ENDIF
         ENDDO
         IF(DUPLIC) THEN
            IJKNEW(IAPX)=LDUPLC
            IF(OUT) THEN
               WRITE(IW,9999) IAPX,IJKOLD,LDUPLC
            ENDIF
         ELSE
            NEWAPX=NEWAPX+1
            JAPX=MAPEX+NEWAPX
            IJKNEW(IAPX)=JAPX
            DO M=1,3
               APEX(M,JAPX)=XYZ(M,IAPX)
            ENDDO
            IF(OUT) THEN
               WRITE(IW,9998) IAPX,IJKOLD,JAPX,
     1                        APEX(1,JAPX),APEX(2,JAPX),APEX(3,JAPX)
            ENDIF
         ENDIF
      ENDDO
      MAPEX=MAPEX+NEWAPX
C
C     ----- MAKE UP NEW FACES AND THEIR CENTERS -----
C
      THIRD=ONE/THREE
      DMIN =ONE
C
      MFACE=MFACE+1
      IJKSEG(  MFACE)=LFACE
      IJKFAC(1,MFACE)=IJKOLD(1)
      IJKFAC(2,MFACE)=IJKNEW(1)
      IJKFAC(3,MFACE)=IJKNEW(3)
      DO M=1,3
         XYZIJK(M)=(APEX(M,IOLD)+APEX(M,JOLD)+APEX(M,KOLD))*THIRD
      ENDDO
      DIJK=DIST(XYZIJK(1),XYZIJK(2),XYZIJK(3),ZERO,ZERO,ZERO)
      DMIN=MIN(DMIN,DIJK)
C
      MFACE=MFACE+1
      IJKSEG(  MFACE)=LFACE
      IJKFAC(1,MFACE)=IJKOLD(2)
      IJKFAC(2,MFACE)=IJKNEW(1)
      IJKFAC(3,MFACE)=IJKNEW(2)
      DO M=1,3
         XYZIJK(M)=(APEX(M,IOLD)+APEX(M,JOLD)+APEX(M,KOLD))*THIRD
      ENDDO
      DIJK=DIST(XYZIJK(1),XYZIJK(2),XYZIJK(3),ZERO,ZERO,ZERO)
      DMIN=MIN(DMIN,DIJK)
C
      MFACE=MFACE+1
      IJKSEG(  MFACE)=LFACE
      IJKFAC(1,MFACE)=IJKOLD(3)
      IJKFAC(2,MFACE)=IJKNEW(2)
      IJKFAC(3,MFACE)=IJKNEW(3)
      DO M=1,3
         XYZIJK(M)=(APEX(M,IOLD)+APEX(M,JOLD)+APEX(M,KOLD))*THIRD
      ENDDO
      DIJK=DIST(XYZIJK(1),XYZIJK(2),XYZIJK(3),ZERO,ZERO,ZERO)
      DMIN=MIN(DMIN,DIJK)
C
      MFACE=MFACE+1
      IJKSEG(  MFACE)=LFACE
      IJKFAC(1,MFACE)=IJKNEW(1)
      IJKFAC(2,MFACE)=IJKNEW(2)
      IJKFAC(3,MFACE)=IJKNEW(3)
      DO M=1,3
         XYZIJK(M)=(APEX(M,IOLD)+APEX(M,JOLD)+APEX(M,KOLD))*THIRD
      ENDDO
      DIJK=DIST(XYZIJK(1),XYZIJK(2),XYZIJK(3),ZERO,ZERO,ZERO)
      DMIN=MIN(DMIN,DIJK)
C
      IF(OUT) THEN
         WRITE(IW,9997) DMIN,MFACE
      ENDIF
C
      RETURN
 9999 FORMAT(' DUPLICATED APEX =',I2,' FOR FACE ',3I5,'. SAME AS = ',I5)
 9998 FORMAT('    NEW     APEX =',I2,' FOR FACE ',3I5,'.  NEWAPX = ',I5,
     1       /,7X,3F12.8)
 9997 FORMAT(' --- DMIN = ',F12.8,' --- MFACE = ',I10)
      END
      subroutine cosmo_charges(rtdb,basis,geom,ecos,iter,some,itype_wfn,
     1                         g_movecs_dft)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      integer rtdb
      integer basis
      integer geom
      logical status
c
      parameter (mxatom=nw_max_atom)
      parameter (mxefc =3*7*7*7)

      parameter (zero=0.0d+00)
c
      integer itype_wfn

      integer   g_movecs_dft
      dimension g_movecs_dft(2)
c
      logical dbug
      logical out
      logical some
      character*8 wfntyp
      character*8 scftyp
      character*8 efclab
c
      common/hnd_iofile/ir,iw,ip
      common/hnd_lcapid/nap,iap
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_basnum/num
      common/hnd_cosrad/radius(mxatom)
      common/hnd_cosprt/iprt
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosbem/efcs(mxefc),iatefc(mxefc)
      common/hnd_efcpar/efcc(3,mxefc),efcz(mxefc),efclab(mxefc),nefc
      common/hnd_efcopt/iefc
      common/hnd_wfnopt/wfntyp
      common/hnd_scfopt/scftyp
      character*16 atom_tags(mxatom)
c
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=some.or.out
c
c     dbug=dbug.and.(iap.eq.1)
c     out =out .and.(iap.eq.1)
c     some=some.and.(iap.eq.1)
c
      if(some) then
         write(iw,9999)
CBERT    call hnd_fiflsh(iw)
      endif
c
c     ----- nwchem/hondo interface -----
c
CBERT call hnd_prp_init(geom,basis)
CBERT call hnd_prp_nw_hnd
CBERT call hnd_atlabl(rtdb,geom,basis)
      if(out) then
         write(iw,*) 'in cosmo_charges ... geom = ',geom
      endif
c
c     ----- read -efc- coordinates from -rtdb- -----
c
      if(.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     &   call errquit('cosmo_charges: rtdb get failed for nefc  ',911,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,efcc))
     &   call errquit('cosmo_charges: rtdb get failed for efcc  ',912,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,  nefc,efcz))
     &   call errquit('cosmo_charges: rtdb get failed for efcz  ',913,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcs',mt_dbl,  nefc,efcs))
     &   call errquit('cosmo_charges: rtdb get failed for efcs  ',914,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:iatefc',mt_int,nefc,iatefc))
     &   call errquit('cosmo_charges: rtdb get failed for iatefc',915,
     &       RTDB_ERR)
      if(out) then
         write(iw,*) 'in cosmo_charges, nefc = ',nefc
         do ief=1,nefc
            write(iw,*) efcc(1,ief),efcc(2,ief),efcc(3,ief)
         enddo
      endif
c
c     ----- partition memory -----
c
      l1= num
      l2=(num*(num+1))/2
      l3= num* num
c
      i10 =1      
      i20 =i10+l2
      i30 =i20+l2
      i40 =i30+l3
      i50 =i40+l3
      last=i50-1 
c
c     ----- allocate memory block -----
c
      need = last
      if(.not.ma_push_get(mt_dbl,need,
     &   'mem init:cosmo:cosmo_charges:1',i_init,init))
     &    call errquit('cosmo_charges, malloc of init failed',911,
     &       MA_ERR)
c
      i10 =init   
      i20 =i10+l2
      i30 =i20+l2
      i40 =i30+l3
      i50 =i40+l3
      last=i50-1 
c
      wfntyp='SCF     '
      scftyp='RHF     '
c
c     ----- get density matrix -----
c
      do i=1,l2
         x(i+i10-1)=zero
         x(i+i20-1)=zero
      enddo
c
      if(out) then
         write(iw,*) 'in -cosmo_charges, wfntyp = ',wfntyp
         write(iw,*) 'in -cosmo_charges, scftyp = ',scftyp
         write(iw,*) 'itype_wfn = ',itype_wfn
         write(iw,*) 'calling hnd_cos_get_dens_... '
         call util_flush(iw)
      endif

      if(itype_wfn.eq.1) then
CBERT    call hnd_cos_get_dens_scf(rtdb,geom,basis,
CBERT1                            x(i10),x(i20),x(i30),x(i40),
CBERT2                            scftyp)
      elseif(itype_wfn.eq.2) then
CBERT    call hnd_cos_get_dens_dft(rtdb,geom,basis,
CBERT1                            x(i10),x(i20),x(i30),x(i40),
CBERT2                            scftyp,g_movecs_dft)
      else
         call errquit('cosmo_charges, wrong wfn flag',911,
     &       INPUT_ERR)
      endif
c
      if(out) then
         write(iw,*) 'in -cosmo_charges, wfntyp = ',wfntyp
         write(iw,*) 'in -cosmo_charges, scftyp = ',scftyp
         if(dbug) then
            write(iw,*) 'in -cosmo_charges, -da- ...'
CBERT       call hnd_prtr(x(i10),l1)
            if(scftyp.eq.'UHF     ') then
               write(iw,*) 'in -cosmo_charges, -db- ...'
CBERT          call hnd_prtr(x(i20),l1)
            endif
         endif
      endif
c
c     ----- get -cosmo- charges -----
c
      call hnd_coschg(X(i10),X(i20),l2,l1,ecos,iter,rtdb,geom)
c
c     ----- release memory block -----
c
      if(.not.ma_pop_stack(i_init))
     & call errquit('cosmo_charges, ma_pop_stack of init failed',911,
     &       MA_ERR)
c
      return
 9999 format(/,10X,15(1H-),
     1       /,10X,'-cosmo- charges',
     2       /,10X,15(1H-))
      end
      SUBROUTINE HND_COSCHG(DA,DB,L2,L1,ECOS,ITER,RTDB,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "prop.fh"
#include "geom.fh"
      INTEGER  RTDB
      INTEGER  GEOM
      LOGICAL  STATUS
C
      PARAMETER  (MXATOM=NW_MAX_ATOM)
      PARAMETER  (MXEFC =3*7*7*7)
      LOGICAL     DBUG
      LOGICAL     MORE
      LOGICAL     OUT
      LOGICAL     DIRECT
      LOGICAL     NOALL
      LOGICAL     ALL
      LOGICAL     ELEC
      LOGICAL     NUCL
      CHARACTER*8 ERRMSG
      CHARACTER*8 EFCLAB
      CHARACTER*8 SCFTYP
      CHARACTER*8 WFNTYP
      CHARACTER*8 SCF
      CHARACTER*8 UHF
      CHARACTER*8 MC
      character*4 KEYEFC,IEFC
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_COSMOD/LINEQ,MINBEM,MAXBEM,IFICOS,IFSCRN
      COMMON/HND_COSDAT/ADIAG,DSURF,DVOL,SRFMOL,VOLMOL 
      COMMON/HND_COSPAR/DIELEC,SCREEN,RSOLV
      COMMON/HND_COSNRG/ATMEFC,ELCEFC,EFCEFC,ALLEFC,SOLNRG,CAVDSP
      COMMON/HND_COSBEM/EFCS(MXEFC),IATEFC(MXEFC)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
      DIMENSION   DA(*)
      DIMENSION   DB(*)
      DIMENSION   ERRMSG(3)
      EXTERNAL    HND_COSAX,HND_COSXA
      EXTERNAL    HND_COSAXD,HND_COSXAD
      DATA ERRMSG /'PROGRAM ','STOP IN ','-COSCHG-'/
      DATA ZERO   /0.0D+00/
      DATA PT5    /0.5D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA MXDIM  /1024/
      DATA SCF    /'SCF     '/
      DATA UHF    /'UHF     '/
      DATA MC     /'MC      '/
      DATA KEYEFC /' EFC'/
C
      MORE=.FALSE.
      DBUG=.FALSE.
      DBUG=DBUG.OR.MORE
      OUT =.FALSE. 
      OUT =OUT.OR.DBUG
C
      DBUG=DBUG.AND.(IAP.EQ.1)
      MORE=MORE.AND.(IAP.EQ.1)
      OUT =OUT .AND.(IAP.EQ.1)
C
      IF(OUT) THEN
         WRITE(IW,9999)
         WRITE(IW,9986) WFNTYP
         WRITE(IW,9985) SCFTYP
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
      IF(NEFC.LE.0) THEN
         WRITE(IW,*) 'NO -COSMO- SURFACE POINTS. STOP.'
         CALL HND_HNDERR(3,ERRMSG)
      ELSE
         IF(OUT) THEN
            WRITE(IW,*) '-NEFC- = ',NEFC
            WRITE(IW,*) 'ITER   = ',ITER
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
      ENDIF
C
C     ----- DENSITY MATRIX -----
C
      IF(WFNTYP.EQ.SCF) THEN
         IF(SCFTYP.EQ.UHF) THEN
            IF(OUT) THEN
               WRITE(IW,*) 'ALPHA DENSITY'
               IF(DBUG) THEN
                  CALL HND_PRTR(DA,L1)
               ENDIF
               WRITE(IW,*) 'BETA  DENSITY'
               IF(DBUG) THEN
                  CALL HND_PRTR(DB,L1)
               ENDIF
CBERT          CALL HND_FIFLSH(IW)
            ENDIF
            DO I=1,L2
               DA(I)=DA(I)+DB(I)
            ENDDO
         ENDIF
         IF(OUT) THEN
            WRITE(IW,*) 'TOTAL DENSITY'
            IF(DBUG) THEN
               CALL HND_PRTR(DA,L1)
            ENDIF
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
      ENDIF
C
      ECOS=ZERO
C
C     ----- ALLOCATE MEMORY BLOCK -----
C
      IF(NEFC.LE.MXDIM.AND.LINEQ.NE.1) THEN
         NEED=3*NEFC+6*NEFC+NEFC*NEFC+NEFC+NEFC
         IF(.NOT.MA_PUSH_GET(MT_DBL,NEED,
     &      'MEM INIT:cosmo:hnd_coschg:1',I_INIT,INIT)) THEN
            NEED=NEFC+NEFC+NEFC+NEFC+NEFC+NEFC+NEFC
            LINEQ=1
            IF(.NOT.MA_PUSH_GET(MT_DBL,NEED,
     &         'MEM INIT:cosmo:hnd_coschg:1',I_INIT,INIT))
     &         CALL ERRQUIT(
     &              'HND_COSCHG, MALLOC OF INIT  FAILED',911)
            IF(OUT) THEN
               WRITE(IW,*) 'USE ITERATIVE METHOD FOR -LINEQ-'
            ENDIF
         ELSE
            IF(OUT) THEN
               WRITE(IW,*) 'USE IN-MEMORY METHOD FOR -LINEQ-'
            ENDIF
         ENDIF
      ELSE
         NEED=3*NEFC+3*NEFC+NEFC+NEFC+NEFC+NEFC+NEFC
         IF(.NOT.MA_PUSH_GET(MT_DBL,NEED,
     &      'MEM INIT:cosmo:hnd_coschg:2',I_INIT,INIT))
     &    CALL ERRQUIT('HND_COSCHG, MALLOC OF INIT  FAILED',912)
         LINEQ=1
         IF(OUT) THEN
            WRITE(IW,*) 'USE ITERATIVE METHOD FOR -LINEQ-'
         ENDIF
      ENDIF
C
C     ----- GET ELECTROSTATIC POTENTIAL AT SURFACE POINTS -----
C
      I10=INIT
      I11=I10+NEFC
      I12=I11+NEFC
      I20=I12+NEFC

      CALL HND_COSELP(DA,X(I10),X(I11),X(I12))
C
C     ----- ( NUCL, ELEC, ALL ) PASSES ... -----
C
      IF(IFSCRN.NE.2) THEN
         NOALL=.TRUE.
      ELSE
         NOALL=.FALSE.
      ENDIF

      IF(NOALL) THEN
         ALL =.TRUE.
         ELEC=.FALSE.
         NUCL=.FALSE.
      ELSE
         ALL  =.FALSE.
         ELEC =.FALSE.
         NUCL =.TRUE.
      ENDIF
C
   10 CONTINUE

      IF(OUT) THEN
         WRITE(IW,*) 'IN -COSCHG- ... NOALL, ALL, ELEC, NUCL = ',
     1                                NOALL, ALL, ELEC, NUCL
      ENDIF
C
C     ----- GET SURFACE CHARGES -----
C
      I10=INIT
      I11=I10+NEFC
      I12=I11+NEFC
      I20=I12+NEFC
      I21=I20+NEFC
      I22=I21+NEFC
      I30=I22+NEFC
      I40=I30+NEFC*NEFC
      I50=I40+NEFC
      I60=I50+NEFC
C
      IF(LINEQ.EQ.0) THEN
         I10=INIT
         I11=I10+NEFC
         I12=I11+NEFC
         I20=I12+NEFC
         I21=I20+NEFC
         I22=I21+NEFC
         I30=I22+NEFC
         I40=I30+NEFC*NEFC
         I50=I40+NEFC
         I60=I50+NEFC
C
C     ----- IN MEMORY -----
C
         IF(OUT) THEN
            WRITE(IW,*) '-LINEQ- IN MEMORY'
CBERT       CALL HND_FIFLSH(IW)
         ENDIF

         IF(ALL) THEN
            IEP=I10
         ELSEIF(NUCL) THEN
            IEP=I11
         ELSEIF(ELEC) THEN
            IEP=I12
         ELSE
            WRITE(IW,*) '... SHOULD NOT BE HERE ... ALL,ELEC,NUCL = ',
     1                                              ALL,ELEC,NUCL
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF

         DO IEF=1,NEFC
            X(IEF+I20-1)=X(IEF+IEP-1)
         ENDDO

         DO JEF=1,NEFC
            XJ=EFCC(1,JEF)
            YJ=EFCC(2,JEF)
            ZJ=EFCC(3,JEF)
            DO IEF=1,NEFC
               IF(IEF.EQ.JEF) THEN
                  AII=ADIAG/SQRT(EFCS(IEF))
                  X(IEF+NEFC*(JEF-1)+I30-1)=AII
               ELSE
                  XI=EFCC(1,IEF)
                  YI=EFCC(2,IEF)
                  ZI=EFCC(3,IEF)
                  DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
                  X(IEF+NEFC*(JEF-1)+I30-1)=ONE/DIJ
               ENDIF
            ENDDO
         ENDDO
         IF(OUT) THEN
            WRITE(IW,*) '-B-'
            IF(MORE) THEN
               CALL HND_PRSQ(X(I20),1,NEFC,NEFC)
            ENDIF
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
         IF(OUT) THEN
            WRITE(IW,*) '-A-'
            IF(MORE) THEN
               CALL HND_PRSQ(X(I30),NEFC,NEFC,NEFC)
            ENDIF
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
C
         NODCMP=1
         CALL HND_LINEQU(X(I30),NEFC,X(I20),NEFC,
     1                   X(I40),X(I50),DETA,IERR,NODCMP)
C
         IF(OUT) THEN
            WRITE(IW,*) '-Q-'
            IF(DBUG) THEN
               CALL HND_PRSQ(X(I20),1,NEFC,NEFC)
            ENDIF
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
      ELSE
         I10=INIT
         I11=I10+NEFC
         I12=I11+NEFC
         I20=I12+NEFC
         I21=I20+NEFC
         I22=I21+NEFC
         I30=I22+NEFC
         I40=I30+NEFC
         I50=I40+NEFC
         I60=I50+NEFC
         I70=I60+NEFC
         I80=I70+NEFC
C
C     ----- ITERATIVE -----
C
         IF(OUT) THEN
            WRITE(IW,*) 'ITERATIVE -LINEQ-'
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
C
         DO IEF=1,NEFC
            X(IEF+I20-1)=ZERO
         ENDDO
C
         DIRECT=.TRUE.
C
C     ----- FIRST TIME AROUND CALCULATE MATRIX -----
C
         IF(.NOT.DIRECT) THEN
            IF(WFNTYP.EQ.MC) THEN
               IF(ITER.EQ.1) THEN
                  CALL HND_COSMAT(X(I30),NEFC)
               ENDIF
            ELSE
               IF(ITER.EQ.0) THEN
                  CALL HND_COSMAT(X(I30),NEFC)
               ENDIF
            ENDIF
         ENDIF
C
C     ----- SOLVE ... -----
C
         IF(ALL) THEN
            IEP=I10
         ELSEIF(NUCL) THEN
            IEP=I11
         ELSEIF(ELEC) THEN
            IEP=I12
         ENDIF

         CALL HND_COSEQU(X(IEP),X(I20),X(I30),NEFC,
     1                   X(I40),X(I50),X(I60),X(I70),
     2                   HND_COSAX,HND_COSXA,
     3                   HND_COSAXD,HND_COSXAD,DIRECT)
C
         IF(OUT) THEN
            WRITE(IW,*) '-Q-'
            IF(MORE) THEN
               CALL HND_PRSQ(X(I20),1,NEFC,NEFC)
            ENDIF
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
      ENDIF
C
C     ----- DEBUG CHECK -----
C
      IF(MORE) THEN
         DO JEF=1,NEFC
            XJ=EFCC(1,JEF)
            YJ=EFCC(2,JEF)
            ZJ=EFCC(3,JEF)
            DUM=ZERO
            DO IEF=1,NEFC
               IF(IEF.EQ.JEF) THEN
                  AII=ADIAG/SQRT(EFCS(IEF))
                  DUM=DUM+AII*X(IEF+I20-1)      
               ELSE
                  XI=EFCC(1,IEF)
                  YI=EFCC(2,IEF)
                  ZI=EFCC(3,IEF)
                  DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
                  AIJ=ONE/DIJ
                  DUM=DUM+AIJ*X(IEF+I20-1)      
               ENDIF
            ENDDO
            WRITE(IW,9992) JEF,DUM
         ENDDO
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- CHARGE SCREENING ISSUE ... -----
C
      CHGCOS=ZERO
      DO IEF=1,NEFC
         CHGIEF=-ONE*X(IEF+I20-1)
         CHGCOS=CHGCOS+CHGIEF
      ENDDO
      CHGSCR=CHGCOS*SCREEN
      IF(OUT) THEN
         WRITE(IW,*) '-COSMO- CHARGE SCREENING ISSUE ...'
         WRITE(IW,*) ' COSMO  CHARGE UN-SCREENED    = ',CHGCOS
         WRITE(IW,*) '                  SCREENING   = ',SCREEN
         WRITE(IW,*) ' COSMO  CHARGE    SCREENED    = ',CHGSCR
      ENDIF

c
c     ----- number of electrons and charge from rtdb ... -----
c     Figure out the number of electrons from the required total
c     charge and the sum of nuclear charges
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, charge_nuclear))
     $     call errquit('scf: geom_nuc_charge failed', 0, GEOM_ERR)
      nelec = nint(charge_nuclear - charge)
      if (nelec .le. 0) call errquit
     $     ('scf: negative no. of electrons ?', nelec, INPUT_ERR)
      if (abs(charge_nuclear - charge - dble(nelec)) .gt. 1d-8)
     $     call errquit('scf: non-integral no. of electrons ?', 0,
     &       INPUT_ERR)

      if(out) then
         write(iw,*) 'from -rtdb- nelec  = ',nelec
         write(iw,*) 'from -rtdb- charge = ',charge
         write(iw,*) 'from -rtdb- charge_nuclear = ',charge_nuclear
      endif

      NE=nelec

      IF(ALL) THEN
         CHGGSS=(         -charge)*SCREEN
      ELSEIF(ELEC) THEN
         CHGGSS=( dble(ne)       )*SCREEN
      ELSEIF(NUCL) THEN
         CHGGSS=(-charge_nuclear )*SCREEN
      ENDIF
      ERRCOS=CHGCOS-CHGGSS
      ERRSCR=CHGSCR-CHGGSS
      IF(OUT) THEN
         IF(ALL) THEN
            WRITE(IW,*) ' SOLUTE TOTAL CHARGE          = ',
     $                  charge
         ELSEIF(ELEC) THEN
            WRITE(IW,*) ' SOLUTE TOTAL ELECTRON CHARGE = ',
     $                 -dble(ne)
         ELSEIF(NUCL) THEN
            WRITE(IW,*) ' SOLUTE TOTAL NUCLEAR  CHARGE = ',
     $                  charge_nuclear
         ENDIF
         WRITE(IW,*) ' GAUSS  CHARGE                = ',CHGGSS
         WRITE(IW,*) ' ERROR  (UN-SCREENED - GAUSS) = ',ERRCOS
         WRITE(IW,*) ' ERROR  (   SCREENED - GAUSS) = ',ERRSCR
         WRITE(IW,*) ' RATIO  (GAUSS / UN-SCREENED) = ',CHGGSS/CHGCOS
         WRITE(IW,*) ' RATIO  (GAUSS /    SCREENED) = ',CHGGSS/CHGSCR
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- SET SCREENING FACTOR -----
C
      IF(ALL) THEN
         IF(IFSCRN.EQ.0) THEN
            CHGFAC=ONE
         ELSEIF(IFSCRN.EQ.1) THEN
            CHGFAC=SCREEN
            CHGFAC=ONE
         ELSEIF(IFSCRN.EQ.2) THEN
            CHGFAC=CHGGSS/CHGCOS
            CHGFAC=ONE
         ELSE
            WRITE(IW,*) 'ILLEGAL -IFSCRN- ... STOP. -IFSCRN- = ',IFSCRN
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         IF(IFSCRN.EQ.0) THEN
            CHGFAC=ONE
         ELSEIF(IFSCRN.EQ.1) THEN
            CHGFAC=SCREEN
         ELSEIF(IFSCRN.EQ.2) THEN
            CHGFAC=CHGGSS/CHGCOS
         ELSE
            WRITE(IW,*) 'ILLEGAL -IFSCRN- ... STOP. -IFSCRN- = ',IFSCRN
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF

      IF(OUT) THEN
         WRITE(IW,*) ' IFSCRN                       = ',IFSCRN
         WRITE(IW,9984) CHGFAC
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- APPLY SCREENING FACTOR -----
C
      DO IEF=1,NEFC
         X(IEF+I20-1)=(-CHGFAC)*X(IEF+I20-1)
      ENDDO
      IF(NUCL) THEN
         DO IEF=1,NEFC
            X(IEF+I21-1)=X(IEF+I20-1)
         ENDDO
      ENDIF
      IF(ELEC) THEN
         DO IEF=1,NEFC
            X(IEF+I22-1)=X(IEF+I20-1)
         ENDDO
      ENDIF
C
C     ----- STORE EFFECTIVE CHARGES IN -EFCZ- ... -----
C           CHECK CONVERGENCE ...
C
      IF(IEFC.EQ.KEYEFC) THEN
         IF(ALL) THEN
            CHGCVG=ZERO
            DO IEF=1,NEFC
               OLDIEF=EFCZ(IEF)
               CHGIEF=   X(IEF+I20-1)
               DELCHG=ABS(CHGIEF-OLDIEF)
               IF(DELCHG.GT.CHGCVG) THEN
                  CHGCVG=DELCHG
               ENDIF
               EFCZ(IEF)=CHGIEF
            ENDDO
         ELSE
            DO IEF=1,NEFC
               EFCZ(IEF)=X(IEF+I20-1)
            ENDDO
         ENDIF
      ELSE
         IF(ALL) THEN
            CHGCVG=ZERO
            DO IEF=1,NEFC
               OLDIEF=   ZERO
               CHGIEF=   X(IEF+I20-1)
               DELCHG=ABS(CHGIEF-OLDIEF)
               IF(DELCHG.GT.CHGCVG) THEN
                  CHGCVG=DELCHG
               ENDIF
               EFCZ(IEF)=CHGIEF
            ENDDO
            IEFC=KEYEFC
         ELSE
            DO IEF=1,NEFC
               EFCZ(IEF)=X(IEF+I20-1)
            ENDDO
         ENDIF
      ENDIF

      IF(ALL) THEN
         IF(OUT) THEN
            WRITE(IW,9987) CHGCVG
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
      ENDIF

      IF(DBUG) THEN
         WRITE(IW,9998)
         DO IEF=1,NEFC
            WRITE(IW,9997) IEF,(EFCC(I,IEF),I=1,3),EFCZ(IEF)
         ENDDO
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- CALCULATE ENERGY TERMS FROM -COSMO- CHARGES -----
C
      ALLEFC=ZERO
      ATMEFC=ZERO
      EFCEFC=ZERO
      DO JEF=1,NEFC
         XJ=EFCC(1,JEF)
         YJ=EFCC(2,JEF)
         ZJ=EFCC(3,JEF)
         QJ=EFCZ(  JEF)
C
         ALLEFC=ALLEFC+QJ*X(JEF+I10-1)
C
         DO IAT=1,NAT
            XI=C(1,IAT)
            YI=C(2,IAT)
            ZI=C(3,IAT)
            QI=ZAN(IAT)
            DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
            BIJ=ONE/DIJ
            ATMEFC=ATMEFC+QI*BIJ*QJ
         ENDDO
C
         DO IEF=1,NEFC
            XI=EFCC(1,IEF)
            YI=EFCC(2,IEF)
            ZI=EFCC(3,IEF)
            QI=EFCZ(  IEF)
            IF(IEF.EQ.JEF) THEN
               AIJ=ADIAG/SQRT(EFCS(IEF))
            ELSE
               DIJ=SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
               AIJ=ONE/DIJ
            ENDIF
            EFCEFC=EFCEFC+QI*AIJ*QJ
         ENDDO
      ENDDO
      EFCEFC= EFCEFC/(TWO*CHGFAC)
      ELCEFC= ALLEFC-ATMEFC
      SOLNRG= ALLEFC+EFCEFC
      ECOS  = ATMEFC+EFCEFC
      IF(OUT) THEN
         WRITE(IW,9991) ATMEFC
         WRITE(IW,9990) ELCEFC
         WRITE(IW,9995) EFCEFC
         WRITE(IW,9988) SOLNRG
         WRITE(IW,9989) ALLEFC,(-TWO*EFCEFC)
         WRITE(IW,9994) ECOS  
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- OTHER FORM OF THE SOLVATION ENERGY ... -----
C
      ALLEFC=ZERO
      ATMEFC=ZERO
      ELCEFC=ZERO
      DO JEF=1,NEFC
         XJ=EFCC(1,JEF)
         YJ=EFCC(2,JEF)
         ZJ=EFCC(3,JEF)
         QJ=EFCZ(  JEF)
C
         ALLEFC=ALLEFC+QJ*X(JEF+I10-1)
         ATMEFC=ATMEFC+QJ*X(JEF+I11-1)
         ELCEFC=ELCEFC+QJ*X(JEF+I12-1)
      ENDDO
      SOLNRG= PT5* ALLEFC
      ECOS  = PT5*(ATMEFC-ELCEFC)
      IF(OUT) THEN
         WRITE(IW,9991) ATMEFC
         WRITE(IW,9990) ELCEFC
         WRITE(IW,9989) ALLEFC
         WRITE(IW,9988) SOLNRG
         WRITE(IW,9994) ECOS
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
C
C     ----- AFTER 'ELEC' PASS, GET THE TOTAL ... -----
C
      IF(NUCL) THEN
         ATMEFCN=ATMEFC
         ELCEFCN=ELCEFC
         EFCEFCN=EFCEFC
         ALLEFCN=ALLEFC
         SOLNRGN=SOLNRG
           ECOSN=  ECOS
      ENDIF
      IF(ELEC) THEN
         ATMEFCE=ATMEFC
         ELCEFCE=ELCEFC
         EFCEFCE=EFCEFC
         ALLEFCE=ALLEFC
         SOLNRGE=SOLNRG
           ECOSE=  ECOS
         ATMEFCT=ATMEFCN+ATMEFCE
         ELCEFCT=ELCEFCN+ELCEFCE
         ALLEFCT=ALLEFCN+ALLEFCE
         SOLNRGT=SOLNRGN+SOLNRGE
           ECOST=  ECOSN+  ECOSE
         IF(OUT) THEN
            WRITE(IW,*) '-ATMEFC (NUCL     )  - = ',ATMEFCN
            WRITE(IW,*) '-ATMEFC (     ELEC)  - = ',ATMEFCE
            WRITE(IW,*) '-ATMEFC (NUCL+ELEC)  - = ',ATMEFCT
            WRITE(IW,*) '-ELCEFC (NUCL     )  - = ',ELCEFCN
            WRITE(IW,*) '-ELCEFC (     ELEC)  - = ',ELCEFCE
            WRITE(IW,*) '-ELCEFC (NUCL+ELEC)  - = ',ELCEFCT
            WRITE(IW,*) '-ALLEFC (NUCL     )  - = ',ALLEFCN
            WRITE(IW,*) '-ALLEFC (     ELEC)  - = ',ALLEFCE
            WRITE(IW,*) '-ALLEFC (NUCL+ELEC)  - = ',ALLEFCT
            WRITE(IW,*) '-SOLNRG (NUCL     )  - = ',SOLNRGN
            WRITE(IW,*) '-SOLNRG (     ELEC)  - = ',SOLNRGE
            WRITE(IW,*) '-SOLNRG (NUCL+ELEC)  - = ',SOLNRGT
            WRITE(IW,*) '-  ECOS (NUCL     )  - = ',  ECOSN
            WRITE(IW,*) '-  ECOS (     ELEC)  - = ',  ECOSE
            WRITE(IW,*) '-  ECOS (NUCL+ELEC)  - = ',  ECOST
            WRITE(IW,*) '-EFCEFC (NUCL     )  - = ',EFCEFCN
            WRITE(IW,*) '-EFCEFC (     ELEC)  - = ',EFCEFCE
         ENDIF
      ENDIF
C
C     ----- AFTER 'ALL' PASS, RESET CHARGES AND ENERGY ... -----
C
      IF(.NOT.NOALL.AND.ALL) THEN
         ATMEFC=ATMEFCT
         ELCEFC=ELCEFCT
         SOLNRG=SOLNRGT
           ECOS=  ECOST
         EFCEFC=ZERO
         DO IEF=1,NEFC
            EFCZ(IEF)=X(IEF+I21-1)+X(IEF+I22-1)
         ENDDO

         IF(DBUG) THEN
            WRITE(IW,9998)
            DO IEF=1,NEFC
               WRITE(IW,9997) IEF,(EFCC(I,IEF),I=1,3),EFCZ(IEF)
            ENDDO
CBERT       CALL HND_FIFLSH(IW)
         ENDIF
      ENDIF
C
C     ----- -NUCL- FIRST, -ELEC- SECOND, THEN -ALL- ... -----
C
      IF(.NOT.ALL) THEN
         IF(.NOT.ELEC) THEN
            ELEC=.TRUE.
            NUCL=.FALSE.
            GO TO 10
         ELSE
            IF(.NOT.NOALL) THEN
               ALL =.TRUE.
               ELEC=.FALSE.
               GO TO 10
            ENDIF
         ENDIF
      ENDIF
C
C     ----- RESET DENSITY MATRIX IF -UHF- BEFORE RETURNING -----
C
      IF(SCFTYP.EQ.UHF) THEN
         DO I=1,L2
            DA(I)=DA(I)-DB(I)
         ENDDO
      ENDIF
C
      IF(OUT) THEN
         WRITE(IW,9993)
CBERT    CALL HND_FIFLSH(IW)
      ENDIF
c
c     ----- save -cosmo- charges and energy to -rtdb- -----
c
      if (.not. rtdb_put(rtdb,'cosmo:energy',mt_dbl,1,ecos))
     &   call errquit('hnd_coschg: rtdb put failed for ecos',911,
     &       RTDB_ERR)
      if(.not.rtdb_put(rtdb,'cosmo:efcz',mt_dbl,  nefc,efcz))
     $   call errquit('hnd_coschg: rtdb put failed for efcz',912,
     &       RTDB_ERR)
c
c     ----- for the time being, save in 'geometry' object -----
c
      if(out) then
         write(iw,*) 'in -hnd_coschg ... for -efc- geom = ',geom
      endif
      status=geom_efc_cart_set(cosmo_geom_efc,nefc,efcc,efcz)
      if (.not.status) then
         call errquit('hnd_coschg: geom_efc_cart_set failed !', 0,
     &       GEOM_ERR)
      endif

C
C     ----- RELEASE MEMORY BLOCK -----
C
      IF(.NOT.MA_POP_STACK(I_INIT))
     & CALL ERRQUIT('HND_COSCHG, MA_POP_STACK OF INIT FAILED',911)
C
      RETURN
 9999 FORMAT(/,10X,15(1H-),
     1       /,10X,'-COSMO- CHARGES',
     2       /,10X,15(1H-))
 9998 FORMAT(4X,'IEFC',6X,'X',5X,6X,'Y',5X,6X,'Z',5X,5X,'Q',4X,
     1     /,1X,53(1H-))
 9997 FORMAT(1X,I7,3F12.6,F10.6)
 9995 FORMAT(' -EFCEFC- ENERGY = ',F20.12)
 9994 FORMAT(' -ECOS  - ENERGY = ',F20.12)
 9993 FORMAT(' ...... END OF -COSCHG- ......')
 9992 FORMAT(' COSMO POTENTIAL AT -IEF = ',I6,F16.10)
 9991 FORMAT(' -ATMEFC- ENERGY = ',F20.12)
 9990 FORMAT(' -ELCEFC- ENERGY = ',F20.12)
 9989 FORMAT(' -ALLEFC- ENERGY = ',F20.12,F20.12)
 9988 FORMAT(' -SOLNRG- ENERGY = ',F20.12)
 9987 FORMAT(' -COSMO- CHARGES CONVERGENCE = ',F10.6)
 9986 FORMAT(' -WFNTYP- = ',A8)
 9985 FORMAT(' -SCFTYP- = ',A8)
 9984 FORMAT(' APPLIED SCREENING FACTOR -CHGFAC- = ',F10.6)
      END
      SUBROUTINE HND_COSAX(A,X,AX,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      COMMON/HND_IOFILE/IR,IW,IP
      DIMENSION A(N),X(N),AX(N)
      DATA NFT17 /17/
C
CBERT CALL HND_REWFIL(NFT17)
      DO I=1,N
         READ(NFT17) A
         DUM=ZERO
         DO J=1,N
            AIJ=A(J)
             XJ=X(J)
            DUM=DUM+AIJ*XJ
         ENDDO
         AX(I)=DUM
      ENDDO
CBERT CALL HND_REWFIL(NFT17)
C
      RETURN
      END
      SUBROUTINE HND_COSXA(A,X,XA,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION A(N),X(N),XA(N)
      DATA NFT17 /17/
C
CBERT CALL HND_REWFIL(NFT17)
      DO I=1,N
         READ(NFT17) A
         DUM=ZERO
         DO J=1,N
            AIJ=A(J)
             XJ=X(J)
            DUM=DUM+AIJ*XJ
         ENDDO
         XA(I)=DUM
      ENDDO
CBERT CALL HND_REWFIL(NFT17)
C
      RETURN
      END
      SUBROUTINE HND_COSMAT(A,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (MXEFC =3*7*7*7)
      CHARACTER*8 EFCLAB
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_COSDAT/ADIAG,DSURF,DVOL,SRFMOL,VOLMOL 
      COMMON/HND_COSBEM/EFCS(MXEFC),IATEFC(MXEFC)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      DIMENSION A(N)
      DATA NFT17 /17/
C
      RIJ(XI,YI,ZI,XJ,YJ,ZJ)=
     1                    SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
C
CBERT CALL HND_REWFIL(NFT17)
      DO IEF=1,NEFC
         XI=EFCC(1,IEF)
         YI=EFCC(2,IEF)
         ZI=EFCC(3,IEF)
         DO JEF=1,NEFC
            XJ=EFCC(1,JEF)
            YJ=EFCC(2,JEF)
            ZJ=EFCC(3,JEF)
            IF(JEF.NE.IEF) THEN
               AIJ=ONE/RIJ(XI,YI,ZI,XJ,YJ,ZJ)
            ELSE
               AIJ=ADIAG/SQRT(EFCS(IEF))
            ENDIF
            A(JEF)=AIJ 
         ENDDO
         WRITE(NFT17) A
      ENDDO
CBERT CALL HND_REWFIL(NFT17)
C
      RETURN
      END
      SUBROUTINE HND_COSAXD(A,X,AX,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (MXEFC =3*7*7*7)
      CHARACTER*8 EFCLAB
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_COSDAT/ADIAG,DSURF,DVOL,SRFMOL,VOLMOL
      COMMON/HND_COSBEM/EFCS(MXEFC),IATEFC(MXEFC)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      DIMENSION A(N),X(N),AX(N)
C
      B(I,J)=ONE/SQRT((EFCC(1,I)-EFCC(1,J))**2+
     1                (EFCC(2,I)-EFCC(2,J))**2+
     2                (EFCC(3,I)-EFCC(3,J))**2)
      D(I  )=ADIAG/SQRT(EFCS(I))
C
      DO I=1,N
         DUM=ZERO
         DO J=1,N
            IF(J.EQ.I) THEN
               AIJ=D(I)
            ELSE
               AIJ=B(I,J)
            ENDIF
             XJ=X(J)
            DUM=DUM+AIJ*XJ
         ENDDO
         AX(I)=DUM
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_COSXAD(A,X,XA,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
      PARAMETER (MXEFC =3*7*7*7)
      CHARACTER*8 EFCLAB
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_COSDAT/ADIAG,DSURF,DVOL,SRFMOL,VOLMOL
      COMMON/HND_COSBEM/EFCS(MXEFC),IATEFC(MXEFC)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      DIMENSION A(N),X(N),XA(N)
C
      B(I,J)=ONE/SQRT((EFCC(1,I)-EFCC(1,J))**2+
     1                (EFCC(2,I)-EFCC(2,J))**2+
     2                (EFCC(3,I)-EFCC(3,J))**2)
      D(I  )=ADIAG/SQRT(EFCS(I))
C
      DO I=1,N
         DUM=ZERO
         DO J=1,N
            IF(J.EQ.I) THEN
               AIJ=D(I)
            ELSE
               AIJ=B(I,J)
            ENDIF
             XJ=X(J)
            DUM=DUM+AIJ*XJ
         ENDDO
         XA(I)=DUM
      ENDDO
C
      RETURN
      END
      SUBROUTINE HND_COSEQU(B,X,A,N,G,H,XI,XJ,
     1                      HND_COSAX,HND_COSXA,
     2                      HND_COSAXD,HND_COSXAD,DIRECT)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- SOLVE A * X = B , USING AN ITERATIVE PROCEDURE       -----
C
C     ----- NUMERICAL RECIPES (P.70), CAMBRIDGE UNIVERSITY PRESS -----
C          W.H.PRESS, B.P.FLANNERY, S.A.TEUKOLSKY, W.T.VETTERLING
C
      CHARACTER*8 ERRMSG
      LOGICAL     DBUG
      LOGICAL     DIRECT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      DIMENSION   ERRMSG(3)
      DIMENSION   A(N),X(N),B(N)
      DIMENSION   G(N),H(N),XI(N),XJ(N)
      EXTERNAL    HND_COSAX,HND_COSXA
      EXTERNAL    HND_COSAXD,HND_COSXAD
      DATA ERRMSG /'PROGRAM ','STOP IN ','-COSEQU-'/
      DATA ZERO   /0.0D+00/
      DATA EPS    /1.0D-07/
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,*) 'IN -COSEQU-'
         DO I=1,N
            WRITE(IW,9999) I,B(I)
         ENDDO
      ENDIF
C
      EPS2=N*EPS*EPS
      IRST=0
   10 IRST=IRST+1
      IF(DIRECT) THEN
         CALL HND_COSAXD(A,X,XI,N)
      ELSE
         CALL HND_COSAX(A,X,XI,N)
      ENDIF
      RP=ZERO
      BSQ=ZERO
      DO 11 J=1,N
         BSQ=BSQ+B(J)*B(J)
         XI(J)=XI(J)-B(J)
         RP=RP+XI(J)*XI(J)
   11 CONTINUE
      IF(DIRECT) THEN
         CALL HND_COSXAD(A,XI,G,N)
      ELSE
         CALL HND_COSXA(A,XI,G,N)
      ENDIF
      DO 12 J=1,N
      G(J)=-G(J)
      H(J)= G(J)
   12 CONTINUE
      DO 19 ITER=1,10*N
         IF(DIRECT) THEN
            CALL HND_COSAXD(A,H,XI,N)
         ELSE
            CALL HND_COSAX(A,H,XI,N)
         ENDIF
         ANUM=ZERO
         ADEN=ZERO
         DO 13 J=1,N
            ANUM=ANUM+G(J)*H(J)
            ADEN=ADEN+XI(J)*XI(J)
   13    CONTINUE
         IF(ADEN.EQ.ZERO) THEN
            WRITE(IW,*) 'VERY SINGULAR MATRIX'
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         ANUM=ANUM/ADEN
         DO 14 J=1,N
            XI(J)=X(J)
            X(J)=X(J)+ANUM*H(J)
   14    CONTINUE
         IF(DIRECT) THEN
            CALL HND_COSAXD(A,X,XJ,N)
         ELSE
            CALL HND_COSAX(A,X,XJ,N)
         ENDIF
         RSQ=ZERO
         DO 15 J=1,N
            XJ(J)=XJ(J)-B(J)
            RSQ=RSQ+XJ(J)*XJ(J)
   15    CONTINUE
         IF(RSQ.EQ.RP.OR.RSQ.LE.BSQ*EPS2) RETURN
         IF(RSQ.GT.RP) THEN
            DO 16 J=1,N
               X(J)=XI(J)
   16       CONTINUE
            IF(IRST.GE.3) RETURN
            GO TO 10
         ENDIF
         RP=RSQ
         IF(DIRECT) THEN
            CALL HND_COSXAD(A,XJ,XI,N)
         ELSE
            CALL HND_COSXA(A,XJ,XI,N)
         ENDIF
         GG=ZERO
         DGG=ZERO
         DO 17 J=1,N
            GG=GG+G(J)*G(J)
            DGG=DGG+(XI(J)+G(J))*XI(J)
   17    CONTINUE
         IF(GG.EQ.ZERO) RETURN
         GAM=DGG/GG
         DO 18 J=1,N
            G(J)=-XI(J)
            H(J)=G(J)+GAM*H(J)
   18    CONTINUE
   19 CONTINUE
      WRITE(IW,*) 'TOO MANY ITERATIONS'
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(I8,F16.8)
      END
      SUBROUTINE HND_COSELP(D,ELP,ELPN,ELPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "nwc_const.fh"
      PARAMETER  (MXATOM=NW_MAX_ATOM)
      PARAMETER  (MXEFC =3*7*7*7)
      CHARACTER*8 EFCLAB
      LOGICAL     SOME
      LOGICAL     OUT
      LOGICAL     DBUG
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFCOPT/IEFC
      DIMENSION      D(*)
      DIMENSION    ELP(*)
      DIMENSION   ELPN(*)
      DIMENSION   ELPE(*)
      DATA ZERO   /0.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
C
C     ----- CALCULATE ELECTROSTATIC POTENTIAL -----
C
      IF(SOME) THEN
         WRITE(IW,9999)
         WRITE(IW,9998)
      ENDIF
C
C     ----- ELECTRON DENSITY CONTIBUTION -----
C
CBERT CALL HND_ELPOTC(NEFC,EFCC,ELPE,D)
      DO IEF=1,NEFC
         ELPE(IEF)=-ELPE(IEF)
      ENDDO
C
C     ----- ADD NUCLEAR CONTRIBUTION -----
C
      IF(DBUG) THEN
         WRITE(IW,9997)
      ENDIF
C
      DO IEF=1,NEFC
         XP = EFCC(1,IEF)
         YP = EFCC(2,IEF)
         ZP = EFCC(3,IEF)
C
C     ----- NUCLEAR CONTRIBUTION ... -----
C
         ELPOTN = ZERO
         DO IAT = 1,NAT
            XN = XP - C(1,IAT) 
            YN = YP - C(2,IAT) 
            ZN = ZP - C(3,IAT) 
            RR = SQRT(XN*XN + YN*YN + ZN*ZN)
            DUM    = ZAN(IAT)/RR
            ELPOTN = ELPOTN +           DUM        
         ENDDO    
         ELPN(IEF) = ELPOTN
C
C     ----- NUCLEAR + ELECTRON CONTRIBUTION ... -----
C
         DUME    =  ELPE(IEF)
         DUMN    =  ELPN(IEF)
         DUM     =  DUMN   + DUME    
         ELP(IEF)=  DUM

         IF(DBUG) THEN
            WRITE(IW,9996) IEF,XP,YP,ZP,DUM,DUMN,DUME 
         ENDIF

      ENDDO
C
      IF(SOME) THEN
         WRITE(IW,9994)
      ENDIF
CBERT CALL HND_FIFLSH(IW)
C
      RETURN
 9999 FORMAT(/,10X,24(1H-),/,10X,'-ELP- AT -COSMO- SURFACE',
     1       /,10X,24(1H-),/)
 9998 FORMAT(' 1 A.U. = 9.07618 ESU/CM ( OR STATVOLTS ) ')
 9997 FORMAT(3X,'POINT',6X,'X',9X,'Y',9X,'Z',5X,'POTENTIAL(A.U.)',/,
     1       83(1H-))
 9996 FORMAT(1X,I5,2X,3F10.5,3F15.6)
 9994 FORMAT(' ...... END OF -COSELP- ......')
      END
      SUBROUTINE HND_LINEQU(A,LDA,B,N,IB,T,DETA,IERR,NODCMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(LDA,1),B(N),IB(N),T(N)
C
C     ----- SOLVE A * X = B , WITH X RETURNED IN B -----
C
      IF(NODCMP.NE.1) GO TO 20
C
CBERT CALL HND_LUDCMP(A,LDA,N,IB,T,DETA,IERR)
      DO 10 J=1,N
   10 DETA=DETA*A(J,J)
C
   20 CONTINUE
CBERT CALL HND_LUBKSB(A,LDA,N,IB,B)
C
      RETURN
      END
      subroutine hnd_cos_get_dens_scf(rtdb,geom_handle,basis_handle,
     1                               da,db,d1nw,d2nw,
     2                               scftyp)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"
#include "cscf.fh"
      integer       rtdb
      integer       geom_handle
      integer       basis_handle
c
      character*8   scftyp
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer g_dens      ! density matrices (up to 3)
      dimension g_dens(3) ! max # of density matrices
 
      integer ndens

c
      integer ir, iw, ip
      common/hnd_iofile/ir,iw,ip

      integer nap, iap
      common/hnd_lcapid/nap,iap
c
      double precision da(*),db(*),d1nw(*),d2nw(*)

      integer i
c
      logical dbug
      dbug=.false.
      dbug=dbug.and.(iap.eq.1)
c
c     ----- ... scftyp ... -----
c
      if (scftype.eq.'UHF') then
         scftyp  = 'UHF     '
      else
         scftyp  = 'RHF     '
      endif
c
      ndens = 3
c
c     handles for density matrices
c
      do i=1, ndens
        g_dens(i) = ga_create_atom_blocked ( geom_handle, 
     $              basis_handle,'density matrix' )
      enddo
c
      if(dbug) then
         write(iw,*) 'in hnd_cos_dens_read'
         write(iw,*) 'scftype = ',scftype
         write(iw,*) 'scftyp  = ',scftyp 
         write(iw,*) 'nbf     = ',nbf
         write(iw,*) 'nmo     = ',nmo
         write(iw,*) 'ndens   = ',ndens
         write(iw,*) 'calling hnd_cos_dens_make_scf ... '
         call util_flush(iw)
      endif
c
c     ----- make density matrices -----
c
      call hnd_cos_dens_make_scf(g_dens, ndens, 
     $                           da, db, d1nw, d2nw)
c
      return
      end
      subroutine hnd_cos_dens_make_scf(g_dens, ndens,
     $                                 da, db, d1nw, d2nw)
      implicit none
#include "global.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "cscf.fh"

      logical status

      integer g_dens, ndens
      dimension g_dens(ndens)

      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      integer ir, iw, ip
      common/hnd_iofile/ir,iw,ip

      integer nap, iap
      common/hnd_lcapid/nap,iap

      double precision da, db, d1nw, d2nw
      dimension da(*), db(*), d1nw(*), d2nw(*)

c     local variables

      integer i, j, ij

c     calculate the density matrices from the g_movecs

      logical dbug
      dbug=.false.
      dbug=dbug.and.(iap.eq.1)

      if(dbug) then
         write(iw,*) 'in hnd_cos_dens_make'
         write(iw,*) 'nclosed, nopen, nbf = ',nclosed,nopen,nbf
         write(iw,*) 'forming density ...'
         call util_flush(iw)
      endif

c       density matrices ...

      if (scftype .eq. 'RHF' .or. scftype .eq. 'ROHF') then
c
c     ----- rhf -----
c
        if(nclosed.eq.0) then
           call ga_zero(g_dens(1))
        else
           if(dbug) then
              write(iw,*) 'density from closed shells ...'
              call util_flush(iw)
           endif

           call ga_matmul_patch('n', 't', 2.0D0, 0.0d0  ,
     $          g_movecs(1), 1, nbf,            1, nclosed,
     $          g_movecs(1), 1, nclosed,        1, nbf    ,
     $          g_dens(1)  , 1, nbf,            1, nbf    )
           if(nopen.gt.0) then
              if(dbug) then
                 write(iw,*) 'density from open   shells ...'
                 call util_flush(iw)
              endif

              call ga_matmul_patch('n', 't', 1.0D0, 0.0d0       ,
     $             g_movecs(1), 1, nbf, nclosed+1, nclosed+nopen,
     $             g_movecs(1), nclosed+1, nclosed+nopen, 1, nbf,
     $             g_dens(2)  , 1, nbf,            1, nbf       )
              call ga_dadd ( 1.D0, g_dens(1), 1.D0, g_dens(2),
     $                                              g_dens(3) )
              call ga_dadd ( 1.D0, g_dens(3), 0.D0, g_dens(2),
     $                                              g_dens(1) )
           endif
           if(dbug) then
              write(iw,*) 'density formed ...'
              call util_flush(iw)
           endif
        endif

        call ga_symmetrize(g_dens(1))

        call ga_get(g_dens(1),1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

CBERT   call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              da(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
CBERT      call hnd_prtr(da,nbf)
           call util_flush(iw)
        endif

      else  ! UHF
c
c     ----- uhf -----
c
        call ga_zero(g_dens(1))
        call ga_zero(g_dens(2))

        call ga_dgemm('n', 't', nbf, nbf, nalpha, 1.D0, 
     $                g_movecs(1), g_movecs(1), 0.D0, g_dens(1))
        if (nbeta .gt. 0) then
           call ga_dgemm('n', 't', nbf, nbf, nbeta, 1.D0, 
     $                   g_movecs(2), g_movecs(2), 0.D0, g_dens(2))
        endif

C       density total

        call ga_dadd ( 1.D0, g_dens(1), 1.D0, g_dens(2), g_dens(3) )
        call ga_symmetrize(g_dens(1))
        call ga_symmetrize(g_dens(2))
        call ga_symmetrize(g_dens(3))

        call ga_get(g_dens(1),1,nbf,1,nbf,d1nw,nbf)

c       hondo to nwchem mapping

CBERT   call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              da(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
CBERT      call hnd_prtr(da,nbf)
        endif

        call ga_get(g_dens(2),1,nbf,1,nbf,d1nw,nbf)
c       hondo to nwchem mapping

CBERT   call hnd_nwhnd_tran(d1nw,d2nw,nbf)

c       square to triangular format

        ij=0
        do i=1,nbf
           do j=1,i
              ij=ij+1
              db(ij)=d2nw(i+nbf*(j-1))
           enddo
        enddo
        if(dbug) then
           call hnd_prtr(db,nbf)
        endif

c
      endif

c
      do i=1,3
         status = ga_destroy ( g_dens(i) )
      enddo
      return
      end
      SUBROUTINE HND_COSDER(RTDB,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
      INTEGER  RTDB
      INTEGER  GEOM
      LOGICAL  GEOM_CREATE
      EXTERNAL GEOM_CREATE
      LOGICAL  GEOM_DESTROY
      EXTERNAL GEOM_DESTROY
      LOGICAL  GEOM_EFC_CART_SET
      EXTERNAL GEOM_EFC_CART_SET
C
      PARAMETER   (MXATOM=NW_MAX_ATOM)
      PARAMETER   (MXEFC =3*7*7*7)
c     CHARACTER*8  ERRMSG
      CHARACTER*8  EFCLAB
      CHARACTER*8 SCFTYP
      CHARACTER*8 WFNTYP
      CHARACTER*8 SCF
      CHARACTER*8 UHF
      character*4 keyblk,iefd
      integer init,i_init,need
      LOGICAL      DBUG
      LOGICAL      OUT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_LCAPID/NAP,IAP
      COMMON/HND_MEMORY/MAXCOR,MAXLCM
      COMMON/HND_COSMOD/LINEQ,MINBEM,MAXBEM,IFICOS,IFSCRN
      COMMON/HND_COSDAT/ADIAG,DSURF,DVOL,SRFMOL,VOLMOL 
      COMMON/HND_COSPAR/DIELEC,SCREEN,RSOLV
      COMMON/HND_COSNRG/ATMEFC,ELCEFC,EFCEFC,ALLEFC,SOLNRG,CAVDSP
      COMMON/HND_COSBEM/EFCS(MXEFC),IATEFC(MXEFC)
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_WFNOPT/WFNTYP
      COMMON/HND_SCFOPT/SCFTYP
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_EFDOPT/IEFD
      COMMON/HND_GRAD12/DE(3,MXATOM)
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
c     DIMENSION ERRMSG(3)
      DIMENSION DE1(3,MXATOM)
c     DATA ERRMSG  /'PROGRAM ','STOP IN ','-COSDER-'/
      DATA SCF     /'SCF     '/
      DATA UHF     /'UHF     '/
      DATA KEYBLK  /'    '/
      DATA ZERO    /0.0D+00/
      DATA TWO     /2.0D+00/
C
      DBUG=.FALSE.
      OUT =.TRUE. 
      OUT =OUT.OR.DBUG
C
      DBUG=DBUG.AND.(IAP.EQ.1)
      OUT =OUT .AND.(IAP.EQ.1)
C
      IF(OUT) THEN
         WRITE(IW,9999)           
         WRITE(IW,9998) IEFC
      ENDIF
C
C     ----- GET CORE MEMORY -----
C
c     CALL HND_CMEM(LOADCM)
c     LOCX1=LOCCM(X(1))
c     NGOTMX=MAXCOR-LOCX1
C
      I10=1
      I20=I10+(NUM*(NUM+1))/2
      I30=I20+(NUM*(NUM+1))/2
C
      LAST=I20+1
      IF(WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF) LAST=I30
      LAST=LAST-1
c     IF(LAST.GT.NGOTMX) THEN
c        NEED=LOCX1+LAST
c        CALL HND_NEEDCM(LAST,NEED)
c        CALL HND_HNDERR(3,ERRMSG)
c     ENDIF
C
C     GET MEMORY BLOCK
C
      need = last
      if(.not.ma_push_get(mt_dbl,need,
     &   'mem init:cosmo:hnd_cosder:1',i_init,init))
     &    call errquit('hnd_cosder, malloc of init failed',911, MA_ERR)
      I10=init
      I20=I10+(NUM*(NUM+1))/2
      I30=I20+(NUM*(NUM+1))/2
      IF(WFNTYP.EQ.SCF.AND.SCFTYP.EQ.UHF) LAST=I30
C
C     ----- GET DENSITY MATRIX -----
C
c--   CALL HND_DENDDM(X(I10),X(I20),NUM)
C
C     ----- GET GRADIENT OF FIELD INTERACTION -----
C
      IEFD=KEYBLK
      CALL HND_COSDIN(X(I10))
C
C     ----- ATM / EFC -COSMO- GRADIENT TERM -----
C
      DO IAT=1,NAT
         DO I=1,3
            DE1(I,IAT)=ZERO
         ENDDO
      ENDDO
      ATMEFC=ZERO
      DO JEF=1,NEFC
         XJ=EFCC(1,JEF)
         YJ=EFCC(2,JEF)
         ZJ=EFCC(3,JEF)
         QJ=EFCZ(  JEF)
         JAT=IATEFC(JEF)
         DO IAT=1,NAT
            XI=C(1,IAT)
            YI=C(2,IAT)
            ZI=C(3,IAT)
            QI=ZAN(IAT)
            RIJ =SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
            RIJ2=RIJ*RIJ
            RIJ3=RIJ*RIJ2
            QIJ =QI*QJ
            DUM =QIJ/RIJ3
            DDXJ= DUM*(XI-XJ)
            DDYJ= DUM*(YI-YJ)
            DDZJ= DUM*(ZI-ZJ)
C
C     -----  THIS IS THE : z*B*q TERM    -----
C           MATCHING TERM FROM -NUCDER- ! 
C
            DE1(1,JAT)=DE1(1,JAT)+DDXJ
            DE1(2,JAT)=DE1(2,JAT)+DDYJ
            DE1(3,JAT)=DE1(3,JAT)+DDZJ
C
            ATMEFC=ATMEFC+QIJ/RIJ
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9997) ATMEFC
         WRITE(IW,9990)
         DO IAT=1,NAT
            WRITE(IW,9993) IAT,(DE1(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- ADD TO TOTAL GRADIENT -----
C
      DO IAT=1,NAT
         DO I=1,3
            DE(I,IAT)=DE(I,IAT)+DE1(I,IAT)
         ENDDO
      ENDDO
C
C     ----- PRINT CURRENT TOTAL GRADIENT -----
C
      IF(OUT) THEN
         WRITE(IW,9992)
         DO IAT=1,NAT
            WRITE(IW,9993) IAT,(DE(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- EFC / EFC -COSMO- GRADIENT TERM -----
C
      DO IAT=1,NAT
         DO I=1,3
            DE1(I,IAT)=ZERO
         ENDDO
      ENDDO
      EFCEFC=ZERO
      DO JEF=1,NEFC
         XJ=EFCC(1,JEF)
         YJ=EFCC(2,JEF)
         ZJ=EFCC(3,JEF)
         QJ=EFCZ(  JEF)
         JAT=IATEFC(JEF)
         DO IEF=1,NEFC
            XI=EFCC(1,IEF)
            YI=EFCC(2,IEF)
            ZI=EFCC(3,IEF)
            QI=EFCZ(  IEF)
            IAT=IATEFC(IEF)
            IF(IEF.EQ.JEF) THEN
               AIJ =ADIAG/SQRT(EFCS(IEF))
               QIJ =QI*QJ/(TWO*SCREEN)
               EFCEFC=EFCEFC+QIJ*AIJ
            ELSE
               RIJ =SQRT((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
               RIJ2=RIJ*RIJ
               RIJ3=RIJ*RIJ2
               QIJ =QI*QJ/(TWO*SCREEN)
               DUM =QIJ/RIJ3 
               DDXJ= DUM*(XI-XJ)
               DDYJ= DUM*(YI-YJ)
               DDZJ= DUM*(ZI-ZJ)
               DDXI=-DDXJ
               DDYI=-DDYJ
               DDZI=-DDZJ
C
C     -----    1/2 * f *  q * A * q term -----
C
               DE1(1,IAT)=DE1(1,IAT)+DDXI
               DE1(2,IAT)=DE1(2,IAT)+DDYI
               DE1(3,IAT)=DE1(3,IAT)+DDZI
               DE1(1,JAT)=DE1(1,JAT)+DDXJ
               DE1(2,JAT)=DE1(2,JAT)+DDYJ
               DE1(3,JAT)=DE1(3,JAT)+DDZJ
C            
               EFCEFC=EFCEFC+QIJ/RIJ   
            ENDIF
         ENDDO
      ENDDO
      ECOS=ATMEFC+EFCEFC
      IF(OUT) THEN
         WRITE(IW,9997) ATMEFC
         WRITE(IW,9996) EFCEFC
         WRITE(IW,9995) ECOS
         WRITE(IW,9994)
         DO IAT=1,NAT
            WRITE(IW,9993) IAT,(DE1(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- ADD TO TOTAL GRADIENT -----
C
      DO IAT=1,NAT
         DO I=1,3
            DE(I,IAT)=DE(I,IAT)+DE1(I,IAT)
         ENDDO
      ENDDO
C
C     ----- PRINT CURRENT TOTAL GRADIENT -----
C
      IF(OUT) THEN
         WRITE(IW,9992)
         DO IAT=1,NAT
            WRITE(IW,9993) IAT,(DE(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- RESET CORE MEMORY -----
C
      if(.not.ma_pop_stack(i_init))
     & call errquit('hnd_cosder, ma_pop_stack of init failed',911,
     &       MA_ERR)
c     CALL HND_CMEM(NGOTCM)
c     IF(NGOTCM.NE.LOADCM) CALL HND_SETC(LOADCM)
C
      IF(OUT) THEN
         WRITE(IW,9991)
      ENDIF
C
      RETURN
 9999 FORMAT(' -COSDER- ',/,' -------- ')
 9998 FORMAT(' -IEFC- = ',A4)
 9997 FORMAT(' -ATMEFC- ENERGY = ',F20.12)
 9996 FORMAT(' -EFCEFC- ENERGY = ',F20.12)
 9995 FORMAT(' -ECOS  - ENERGY = ',F20.12)
 9994 FORMAT(' EFC / EFC -COSMO- GRADIENT CONTRIBUTIONS')
 9993 FORMAT(I5,3F15.8)
 9992 FORMAT(' CURRENT TOTAL GRADIENT')
 9991 FORMAT(' ...... END OF -COSDER- ......')
 9990 FORMAT(' ATM / EFC -COSMO- GRADIENT CONTRIBUTIONS')
      END
      SUBROUTINE HND_COSDIN(DAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
C
#include "hnd_rys.fh"
#include "hnd_tol.fh"
C
      PARAMETER   (MXATOM=NW_MAX_ATOM)
      PARAMETER   (MXEFC =3*7*7*7)
      PARAMETER   (MXPRIM=2048)
      PARAMETER   (MXSHEL=1024)
      CHARACTER*8 EFCLAB
      CHARACTER*8 ERRMSG
      character*4 keyefd,iefd
      LOGICAL     IANDJ
      LOGICAL     NORM
      LOGICAL     DOUBLE
      LOGICAL     DBUG
      LOGICAL     OUT
      LOGICAL     SOME
      LOGICAL     EFD
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW,IP
      COMMON/HND_EFCOPT/IEFC
      COMMON/HND_EFDOPT/IEFD
      COMMON/HND_EFCPAR/EFCC(3,MXEFC),EFCZ(MXEFC),EFCLAB(MXEFC),NEFC
      COMMON/HND_EFDPAR/EFDD(3,MXEFC),EFDALF(MXEFC),NEFD
      COMMON/HND_COSBEM/EFCS(MXEFC),IATEFC(MXEFC)
      COMMON/HND_NSHEL/EX(MXPRIM),CS(MXPRIM),CP(MXPRIM),
     1                 CD(MXPRIM),CF(MXPRIM),CG(MXPRIM),
     2                 KSTART(MXSHEL),KATOM(MXSHEL),KTYPE(MXSHEL),
     3                 KNG(MXSHEL),KLOC(MXSHEL),KMIN(MXSHEL),
     4                 KMAX(MXSHEL),NSHELL
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_BASNUM/NUM
      COMMON/HND_XYZDER/XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     1                  NI,NJ
     1                              ,CX,CY,CZ
      COMMON/HND_GRAD12/DE(3,MXATOM)
      COMMON/HND_IJPAIR/IA(1)
      DIMENSION DAB(*)
      DIMENSION DE1(3,MXATOM)
      DIMENSION DE2(3)
      DIMENSION DE3(3,MXATOM)
      DIMENSION DE4(3,MXATOM)
      DIMENSION DE5(3)
      DIMENSION DIJ(225)
      DIMENSION IJX(35),IJY(35),IJZ(35)
      DIMENSION   XV(5,5,5),  YV(5,5,5),  ZV(5,5,5)
      DIMENSION  DXV(5,5,5), DYV(5,5,5), DZV(5,5,5)
      DIMENSION DDXV(5,5,5),DDYV(5,5,5),DDZV(5,5,5)
      DIMENSION W2(5),W4(5)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-COSDIN-'/
      DATA RLN10  /2.30258D+00/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA FOUR   /4.0D+00/
      DATA PI212  /1.1283791670955D+00/
      DATA SQRT3  /1.73205080756888D+00/
      DATA SQRT5  /2.23606797749979D+00/
      DATA SQRT7  /2.64575131106459D+00/
      DATA IJX    / 1, 2, 1, 1, 3, 1, 1, 2, 2, 1,
     1              4, 1, 1, 3, 3, 2, 1, 2, 1, 2,
     2              5, 1, 1, 4, 4, 2, 1, 2, 1, 3,
     3              3, 1, 3, 2, 2/
      DATA IJY    / 1, 1, 2, 1, 1, 3, 1, 2, 1, 2,
     1              1, 4, 1, 2, 1, 3, 3, 1, 2, 2,
     2              1, 5, 1, 2, 1, 4, 4, 1, 2, 3,
     3              1, 3, 2, 3, 2/
      DATA IJZ    / 1, 1, 1, 2, 1, 1, 3, 1, 2, 2,
     1              1, 1, 4, 1, 2, 1, 2, 3, 3, 2,
     2              1, 1, 5, 1, 2, 1, 2, 4, 4, 1,
     3              3, 3, 2, 2, 3/
      DATA KEYEFD /' EFD'/
C
C     ----- ELECTRON/EFC CONTRIBUTION TO EFC GRADIENTS -----
C     ----- ELECTRON/EFD CONTRIBUTION TO EFC GRADIENTS -----
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.NPRINT.EQ. 6.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
      IF(SOME) THEN
         WRITE(IW,9992)
      ENDIF
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
      DO IAT=1,NAT 
         DO I=1,3
            DE1(I,IAT)=ZERO
         ENDDO
      ENDDO
C
      DO IAT=1,NAT 
         DO I=1,3
            DE3(I,IAT)=ZERO
            DE4(I,IAT)=ZERO
         ENDDO
      ENDDO
C
      EFD=IEFD.EQ.KEYEFD
      IF(SOME) THEN
         WRITE(IW,*) 'EFD = ',EFD
      ENDIF
C
      NDER=1
      IF(EFD) THEN
         NDER=2
      ENDIF
C
C     ----- ISHELL -----
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL -----
C
      DO 8000 JJ=1,II
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
C
      IANDJ=II.EQ.JJ
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      NROOTS=(LIT+LJT+NDER-2)/2+1
      IF(NROOTS.GT.MAXRYS) THEN
         WRITE(IW,9997) MAXRYS,LIT,LJT,NROOTS
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- I PRIMITIVE -----
C
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
      CGI=CG(IG)
C
C     ----- J PRIMITIVE -----
C
      JGMAX=J2
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GO TO 6000
      FAC= EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      CGJ=CG(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR -----
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      IJ=0
      DO 360 I=MINI,MAXI
      GO TO (110,120,220,220,130,220,220,140,220,220,
     1       150,220,220,160,220,220,220,220,220,170,
     2       180,220,220,190,220,220,220,220,220,200,
     3       220,220,210,220,220),I
  110 DUM1=CSI*FAC
      GO TO 220
  120 DUM1=CPI*FAC
      GO TO 220
  130 DUM1=CDI*FAC
      GO TO 220
  140 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  150 DUM1=CFI*FAC
      GO TO 220
  160 IF(NORM) DUM1=DUM1*SQRT5
      GO TO 220
  170 IF(NORM) DUM1=DUM1*SQRT3
      GO TO 220
  180 DUM1=CGI*FAC
      GO TO 220
  190 IF(NORM) DUM1=DUM1*SQRT7
      GO TO 220
  200 IF(NORM) DUM1=DUM1*SQRT5/SQRT3
      GO TO 220
  210 IF(NORM) DUM1=DUM1*SQRT3
  220 CONTINUE
C
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 360 J=MINJ,JMAX
      GO TO (230,250,350,350,260,350,350,270,350,350,
     1       280,350,350,290,350,350,350,350,350,300,
     2       310,350,350,320,350,350,350,350,350,330,
     3       350,350,340,350,350),J
  230 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GO TO 350
      IF(I.GT.1) GO TO 240
      DUM2=DUM2+DUM2
      GO TO 350
  240 DUM2=DUM2+CSI*CPJ*FAC
      GO TO 350
  250 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  260 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  270 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  280 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  290 IF(NORM) DUM2=DUM2*SQRT5
      GO TO 350
  300 IF(NORM) DUM2=DUM2*SQRT3
      GO TO 350
  310 DUM2=DUM1*CGJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GO TO 350
  320 IF(NORM) DUM2=DUM2*SQRT7
      GO TO 350
  330 IF(NORM) DUM2=DUM2*SQRT5/SQRT3
      GO TO 350
  340 IF(NORM) DUM2=DUM2*SQRT3
  350 CONTINUE
C
      LOCII=LOCI+I
      LOCJJ=LOCJ+J
      NN=(LOCII*(LOCII-1))/2+LOCJJ
      DEN=DAB(NN)
      IF(.NOT.IANDJ.OR.J.NE.I) DEN=DEN+DEN
      IJ=IJ+1
  360 DIJ(IJ)=DUM2*DEN
C
C     ----- ELECTRIC FIELD          TERM -----
C           ELECTRIC FIELD GRADIENT TERM
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DO 500 IEF=1,NEFC
      ZNUC=-EFCZ(IEF)
      CX=EFCC(1,IEF)
      CY=EFCC(2,IEF)
      CZ=EFCC(3,IEF)
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL HND_RT123
      IF(NROOTS.EQ.4) CALL HND_ROOT4
      IF(NROOTS.EQ.5) CALL HND_ROOT5
      DO 420 IROOT=1,NROOTS
CBERT UU=U(IROOT)*AA
CBERT U2=UU
CBERT U4=UU*UU
CBERT WW=W(IROOT)
      W2(IROOT)=WW*U2*TWO
      W4(IROOT)=WW*U4*FOUR
      TT=ONE/(AA+UU)
      T = SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
C
      DO 410 J=1,LJT
      NJ=J
      DO 410 I=1,LIT
      NI=I
CBERT CALL HND_DSXYZ
      XV(I,J,IROOT)=XINT
      YV(I,J,IROOT)=YINT
      ZV(I,J,IROOT)=ZINT
CBERT CALL HND_DVXYZ
      DXV(I,J,IROOT)=XINT
      DYV(I,J,IROOT)=YINT
      DZV(I,J,IROOT)=ZINT
      IF(EFD) THEN
CBERT    CALL HND_DDVXYZ
         DDXV(I,J,IROOT)=XINT
         DDYV(I,J,IROOT)=YINT
         DDZV(I,J,IROOT)=ZINT
      ENDIF
  410 CONTINUE
  420 CONTINUE
C
      IJ=0
      DO 450 I=MINI,MAXI
      IX=IJX(I)
      IY=IJY(I)
      IZ=IJZ(I)
      JMAX=MAXJ
      IF(IANDJ) JMAX=I
      DO 440 J=MINJ,JMAX
      JX=IJX(J)
      JY=IJY(J)
      JZ=IJZ(J)
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      DUMXX=ZERO
      DUMYX=ZERO
      DUMZX=ZERO
      DUMXY=ZERO
      DUMYY=ZERO
      DUMZY=ZERO
      DUMXZ=ZERO
      DUMYZ=ZERO
      DUMZZ=ZERO
      DO 430 IROOT=1,NROOTS
      DUMX =DUMX +
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMY =DUMY +
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
      DUMZ =DUMZ +
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W2(IROOT)
      IF(EFD) THEN
         DUM  =
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W2(IROOT)
         DUMXX=DUMXX-DUM+
     1   DDXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMYX=DUMYX+
     1    DXV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMZX=DUMZX+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMXY=DUMXY+
     1    DXV(IX,JX,IROOT)* DYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMYY=DUMYY-DUM+
     1     XV(IX,JX,IROOT)*DDYV(IY,JY,IROOT)*  ZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMZY=DUMZY+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMXZ=DUMXZ+
     1    DXV(IX,JX,IROOT)*  YV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMYZ=DUMYZ+
     1     XV(IX,JX,IROOT)* DYV(IY,JY,IROOT)* DZV(IZ,JZ,IROOT)*W4(IROOT)
         DUMZZ=DUMZZ-DUM+
     1     XV(IX,JX,IROOT)*  YV(IY,JY,IROOT)*DDZV(IZ,JZ,IROOT)*W4(IROOT)
      ENDIF
  430 CONTINUE
      IJ=IJ+1
      DUM=PI212*AA1*DIJ(IJ)
      TEMPX=DUMX*(DUM*ZNUC)
      TEMPY=DUMY*(DUM*ZNUC)
      TEMPZ=DUMZ*(DUM*ZNUC)
      IAT=IATEFC(IEF)
       DE1(1,IAT)= DE1(1,IAT)+TEMPX
       DE1(2,IAT)= DE1(2,IAT)+TEMPY
       DE1(3,IAT)= DE1(3,IAT)+TEMPZ
       DE3(1,IAT)= DE3(1,IAT)+TEMPX
       DE3(2,IAT)= DE3(2,IAT)+TEMPY
       DE3(3,IAT)= DE3(3,IAT)+TEMPZ
      IF(EFD) THEN
         TEMPX=     -DUMXX*(DUM*EFDD(1,IEF))
         TEMPY=     -DUMYX*(DUM*EFDD(1,IEF))
         TEMPZ=     -DUMZX*(DUM*EFDD(1,IEF))
         TEMPX=TEMPX-DUMXY*(DUM*EFDD(2,IEF))
         TEMPY=TEMPY-DUMYY*(DUM*EFDD(2,IEF))
         TEMPZ=TEMPZ-DUMZY*(DUM*EFDD(2,IEF))
         TEMPX=TEMPX-DUMXZ*(DUM*EFDD(3,IEF))
         TEMPY=TEMPY-DUMYZ*(DUM*EFDD(3,IEF))
         TEMPZ=TEMPZ-DUMZZ*(DUM*EFDD(3,IEF))
          DE1(1,IAT)= DE1(1,IAT)+TEMPX
          DE1(2,IAT)= DE1(2,IAT)+TEMPY
          DE1(3,IAT)= DE1(3,IAT)+TEMPZ
          DE4(1,IAT)= DE4(1,IAT)+TEMPX
          DE4(2,IAT)= DE4(2,IAT)+TEMPY
          DE4(3,IAT)= DE4(3,IAT)+TEMPZ
      ENDIF
  440 CONTINUE
  450 CONTINUE
C
  500 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
      IF(DBUG) THEN
         DO I=1,3
            DE2(I)=ZERO
            DE5(I)=ZERO
            DO IAT=1,NEFC
               DE2(I)=DE2(I)+DE3(I,IAT)
               DE5(I)=DE5(I)+DE4(I,IAT)
            ENDDO
         ENDDO
         WRITE(IW,9995) II,JJ
         WRITE(IW,*)    '-DE3-'
         DO IAT=1,NAT 
            WRITE(IW,9994) IAT,DE3(1,IAT),DE3(2,IAT),DE3(3,IAT)
         ENDDO
         IF(EFD) THEN
            WRITE(IW,9995) II,JJ
            WRITE(IW,*)    '-DE4-'
            DO IAT=1,NAT 
               WRITE(IW,9994) IAT,DE4(1,IAT),DE4(2,IAT),DE4(3,IAT)
            ENDDO
         ENDIF
         WRITE(IW,9995) II,JJ
         WRITE(IW,*) 'TRANSLATIONAL INVARIANCE'
         WRITE(IW,*) '-X- = ',DE2(1)
         WRITE(IW,*) '-Y- = ',DE2(2)
         WRITE(IW,*) '-Z- = ',DE2(3)
         WRITE(IW,9995) II,JJ
         DO IAT=1,NAT
            WRITE(IW,9994) IAT,DE1(1,IAT),DE1(2,IAT),DE1(3,IAT)
         ENDDO
      ENDIF
C
 8000 CONTINUE
 9000 CONTINUE
C
      IF(OUT) THEN
         WRITE(IW,*) '-DE3- = ELC / EFC GRADIENT TERM FROM -COSMO-'
         DO IAT=1,NAT
            WRITE(IW,9994) IAT,DE3(1,IAT),DE3(2,IAT),DE3(3,IAT)
         ENDDO
         IF(EFD) THEN
            WRITE(IW,*) '-DE4- = ELC / EFD GRADIENT TERM FROM -COSMO-'
            DO IAT=1,NAT 
               WRITE(IW,9994) IAT,DE4(1,IAT),DE4(2,IAT),DE4(3,IAT)
            ENDDO
         ENDIF
         IF(EFD) THEN
            WRITE(IW,*) '-DE3- + -DE4- = '
            DO IAT=1,NAT
               WRITE(IW,9994) IAT,DE1(1,IAT),DE1(2,IAT),DE1(3,IAT)
            ENDDO
         ENDIF
      ENDIF
C
C     ----- ADD TO FINAL GRADIENT -----
C
      DO IAT=1,NAT
         DO I=1,3
            DE(I,IAT)=DE(I,IAT)+DE1(I,IAT)
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,*) 'TOTAL CURRENT GRADIENT'
         DO IAT=1,NAT
            WRITE(IW,9994) IAT,DE(1,IAT),DE(2,IAT),DE(3,IAT)
         ENDDO
      ENDIF
C
      RETURN
 9997 FORMAT(' IN -COSDIN- , THE RYS QUADRATURE IS NOT IMPLEMENTED',
     1       ' BEYOND -NROOTS- = ',I3,/,
     2       ' LIT,LJT,NROOTS= ',3I3)
 9996 FORMAT(' II,JJ,IG,JG,IEF,I,J= ',7I3,/,4F15.8)
 9995 FORMAT(' SHELLS II,JJ = ',2I5)
 9994 FORMAT(I5,3F15.8)
 9993 FORMAT(3F15.8,/,3F15.8)
 9992 FORMAT(' COSDIN',/,' ------')
 9987 FORMAT(' I,J,DUMX ,DUMY ,DUMZ = ',2I3,3F15.8)
 9986 FORMAT(' I,J,DUMXX,DUMXY,DUMXZ= ',2I3,3F15.8)
 9985 FORMAT(' I,J,DUMYX,DUMYY,DUMYZ= ',2I3,3F15.8)
 9984 FORMAT(' I,J,DUMZX,DUMZY,DUMZZ= ',2I3,3F15.8)
      END
