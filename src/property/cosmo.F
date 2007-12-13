      subroutine cosmo_input(rtdb)
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
      integer rtdb
c
      integer iat, nfield, nrad, irad 
      integer iscren, ifscrn, minbem, maxbem, ificos, lineq
      integer l_rad, k_rad
      double precision dielec, rsolv, rad
      character*255 field
      logical do_cosmo
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
      if(.not.ma_push_get(mt_dbl,nw_max_atom,'cosmo rads',l_rad,k_rad))
     & call errquit('cosmo_input malloc k_rad failed',911,MA_ERR)
      do iat=1,nw_max_atom
         dbl_mb(k_rad+iat-1)=0.0d+00
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
            if(inp_f(rad)) then
               irad=irad+1
               dbl_mb(k_rad+irad-1)=rad
               nfield=nfield-1
               go to 30
            else
               nrad=irad  
               if (.not.rtdb_put(rtdb,'cosmo:nrrad',mt_int,1,nrad)) call
     $            errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
               if (.not. rtdb_put(rtdb,'cosmo:radius',mt_dbl,nrad,
     $            dbl_mb(k_rad))) call 
     $            errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
               go to 20
            endif
         else
            if(.not.inp_read()) call
     $         errquit('cosmo_input: unexpected eof',911,INPUT_ERR)
            nfield = inp_n_field()
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
     $       call errquit('cosmo_input: rtdb_put failed',0, RTDB_ERR)
      else
         if (.not. rtdb_put(rtdb,'slv:cosmo',mt_log,1,.false.))
     $       call errquit('cosmo_input: rtdb_put failed',0, RTDB_ERR)
      endif
c
      if(.not.ma_pop_stack(l_rad)) call
     &  errquit('cosmo_input, ma_pop_stack of l_rad failed',911,MA_ERR)
c
      return
      end
c
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
#include "msgids.fh"
c
      integer rtdb
      integer geom
      integer basis
      logical oprint
c
      logical status
      logical cosmo
c
      integer iat, nrad
      integer iscren, ifscrn, minbem, maxbem, ificos, lineq
      integer l_rad, k_rad, l_costag, k_costag, l_coszan, k_coszan
      integer l_coscoor, k_coscoor
c
      logical odbug, osome
      integer me
c
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
      double precision vdwr(102)
      double precision vander(102)
c
      integer i
      integer atn(203)
      double precision arad(203)
      character*16 aname,symb,elem
      integer fn,gn
      character*255 map
      logical util_io_unit
      external util_io_unit
c
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pt5    /0.5d+00/
      data radtol /0.1d+00/
      data mxpass /6/
C
C     ----- van der waals radii in angstroms	-----
C
      data def  /1.80D+00/

c      data vdwr /
c     1   0.80,0.49,0.00,0.00,0.00,1.65,1.55,1.50,1.50,0.00,
c     2   2.30,1.70,2.05,2.10,1.85,1.80,1.80,0.00,2.80,2.75,
c     3   0.00,0.00,1.20,0.00,0.00,0.00,2.70,0.00,0.00,0.00,
c     4   0.00,0.00,0.00,1.90,1.90,0.00,0.00,0.00,0.00,1.55,
c     5   0.00,1.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     6   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     7   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     8   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     9   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.65,
c     2   0.00,0.00,0.00/
c
c    default radii provided by Andreas Klamt (Cosmologic)
c    vdw radii: 1.17 (+/- 0.02) * Bondi radius 
c    (Bondi, J. Phys. Chem., 68, 441, 1964)
c
c    optimal vdw radii for H, C, N, O, F, S, Cl, Br, I
c    (Klamt et al, J. Phys. Chem. A, 102, 5074 (1998)
c
c    for heavy elements: 1.17*1.9
c
      data (vander(i),i=1,102)
     1 / 1.300,1.638,1.404,1.053,2.0475,2.00,
     2   1.830,1.720,1.720,1.8018,1.755,1.638,
     3   1.404,2.457,2.106,2.160,2.05,2.223,
     4   2.223,2.223,2.223,2.223,2.223,2.223,
     5   2.223,2.223,2.223,2.223,2.223,2.223,
     6   2.223,2.223,2.223,2.223,2.160,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223,
     8   2.223,2.223,2.223,2.223,2.223,2.223,
     9   2.223,2.223,2.223,2.223,2.320,2.223,
     1   2.223,2.223,2.223,2.223,2.223,2.223,
     2   2.223,2.223,2.223,2.223,2.223,2.223,
     3   2.223,2.223,2.223,2.223,2.223,2.223,
     4   2.223,2.223,2.223,2.223,2.223,2.223,
     5   2.223,2.223,2.223,2.223,2.223,2.223,
     6   2.223,2.223,2.223,2.223,2.223,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223/
c
c     vdw radii from Merz/Kollman/Singh
c
c      data (vander(i),i=1,36)
c     1                  /1.20D+00,1.20D+00,1.37D+00,1.45D+00,
c     2 1.45D+00,1.50D+00,1.50D+00,1.40D+00,1.35D+00,1.30D+00,
c     3 1.57D+00,1.36D+00,1.24D+00,1.17D+00,1.80D+00,1.75D+00,
c     4 1.70D+00,19*0.0D+00/
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
      me=ga_nodeid()
c
      odbug=.false.
      osome=oprint
c
      if(odbug.and.me.eq.0) then
         write(Luout,*) 'in cosmo_initialize ...'
      endif
c
      do i=1,102
         vdwr(i)=vander(i)*1.d0 
      enddo
c
c     start processing cosmo map
c     --------------------------
      nmap = 0
      if(.not.rtdb_cget(rtdb,'cosmo:map',1,map))
     +  map = "none"
      
      if(map.ne."none") then
        call util_file_name_resolve(map,.false.)
        if(me.eq.0) then
        if(.not.util_io_unit(80,90,fn))
     +     call errquit('cosmo_initialize:cannot get free unit', 0,
     +       0)

        open(unit=fn,status="old",form="formatted",file=map)
        write(luout,9982) map(1:inp_strlen(map))
10      continue
        read(fn,*,end=11) aname                 
        if(aname.ne." ") then
         nmap = nmap + 1
         goto 10
        end if
11      continue
        rewind(fn)
        do i=1,nmap
         read(fn,*) aname,arad(i)
         if(.not.geom_tag_to_element(aname, symb, elem, atn(i))) 
     +     call errquit('cosmo_initialize:cannot resolve map'
     +                  //aname, 0,0)
          write(luout,9981) symb,arad(i)
        end do
        close(fn)
        end if
        call ga_brdcst(msg_cosmo0,nmap,
     >                 ma_sizeof(mt_int,1,mt_byte),0)
        call ga_brdcst(msg_cosmo1,atn,
     >                 nmap*ma_sizeof(mt_int,1,mt_byte),0)
        call ga_brdcst(msg_cosmo2,arad,
     >                 nmap*ma_sizeof(mt_dbl,1,mt_byte),0)
        call ga_sync()

        do i=1,nmap
          vdwr(atn(i)) = arad(i)
        end do

      end if
c
c     end processing cosmo map
c     --------------------------

c
c     Create a geometry for the COSMO efc charges, to be used later on
c
      if (.not. geom_create(cosmo_geom_efc, 'cosmo efc geom'))
     $     call errquit('cosmo_initialize: geom_create failed', 0,
     &       GEOM_ERR)
c
c     Get data from the molecule itself
c
      if(.not.geom_ncent(geom,nat)) call errquit
     $    (' cosmo_initialize: geom_ncent    failed.',911, GEOM_ERR)
c
      if(.not.ma_push_get(mt_dbl,nat,'cosmo rads',l_rad,k_rad))
     & call errquit('cosmo_init malloc k_rad failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nat*3,'coord',l_coscoor,k_coscoor))
     & call errquit('cosmo_init malloc k_coscoor failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nat,'cosmo z',l_coszan,k_coszan)) call
     &  errquit('cosmo_init malloc k_coszan failed',911,MA_ERR)
      if(.not.ma_push_get(mt_byte,nat*16,'tags',l_costag,k_costag)) call
     &  errquit('cosmo_init malloc k_costag failed',911,MA_ERR)
c
      if(.not.geom_cart_get(geom,nat,byte_mb(k_costag),
     &   dbl_mb(k_coscoor),dbl_mb(k_coszan))) call errquit
     $    (' cosmo_initialize: geom_cart_get failed.',911, GEOM_ERR)
c
      do iat=1,nat
         rad=vdwr(nint(dbl_mb(k_coszan+iat-1)))
         if (geom_any_ecp(geom).and.geom_ecp_get(geom,iat)) then
            if (.not.geom_tag_to_element(byte_mb(k_costag+16*(iat-1)), 
     &         symb, elem, gn))
     &          call errquit('cosmo_init:geom_tag_to_element',0,0)
            rad=vdwr(gn)
         endif
         if(rad.eq.zero) then
            rad=def
         endif
         dbl_mb(k_rad+iat-1)=rad
      enddo
c
      if(.not.ma_pop_stack(l_costag)) call
     &  errquit('cosmo_init,ma_pop_stack of l_tag failed',911, MA_ERR)
c
      dielec   =78.4d+00
      iscren   =0
      ifscrn   =2
      minbem   =2
      maxbem   =3
      ificos   =0
      lineq    =1
      rsolv    =0.0d+00
c
      status = rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo)
      status = rtdb_get(rtdb,'cosmo:dielec',mt_dbl,1,dielec) 
      status = rtdb_get(rtdb,'cosmo:rsolv',mt_dbl,1,rsolv)
      status = rtdb_get(rtdb,'cosmo:iscren',mt_int,1,iscren) 
      status = rtdb_get(rtdb,'cosmo:ifscrn',mt_int,1,ifscrn)
      status = rtdb_get(rtdb,'cosmo:minbem',mt_int,1,minbem) 
      status = rtdb_get(rtdb,'cosmo:maxbem',mt_int,1,maxbem)
      status = rtdb_get(rtdb,'cosmo:ificos',mt_int,1,ificos) 
      status = rtdb_get(rtdb,'cosmo:lineq',mt_int,1,lineq) 
      nrad = 0
      if (.not. rtdb_get(rtdb,'cosmo:nrrad',mt_int,1,nrad)) nrad = 0
      if (nrad.gt.0) then
         status=rtdb_get(rtdb,'cosmo:radius',mt_dbl,nrad,dbl_mb(k_rad))
      endif
c
c     ----- checks ... -----
c
      if(dielec.le.zero) then
         dielec=78.4d+00
      endif
c
      if(osome.and.me.eq.0) then
         write(Luout,9999)
         write(Luout,9997) dielec
      endif
c
c     ----- charge screening approach ... -----
c
      if(osome.and.me.eq.0) then
         write(Luout,9983) ifscrn
      endif

      if(ifscrn.ne.1.and.ifscrn.ne.2) then
         iscren=-1
      endif
c
c     ----- screening factor -----
c
      if(iscren.lt.0) then
         screen=one
         if (osome.and.me.eq.0) write(luout,9984) screen
      elseif(iscren.gt.0) then
         screen=(dielec-one)/(dielec+pt5)
         if (osome.and.me.eq.0) write(luout,9996) screen
      else
         screen=(dielec-one)/(dielec    )
         if (osome.and.me.eq.0) write(luout,9995) screen
      endif
c
c     ----- tesselation -----
c
      if(minbem.lt.1.or.minbem.gt.mxpass) then
         minbem=2
      endif
      if(maxbem.lt.1.or.maxbem.gt.mxpass) then
         maxbem=4
      endif
      if(maxbem.lt.minbem) then
         maxbem=minbem+2
      endif
c
      ificos=abs(ificos)
c
      if(osome.and.me.eq.0) then
         write(luout,9993) lineq
         write(luout,9991) minbem
         write(luout,9992) maxbem
         if(ificos.eq.0) then
            write(luout,9990)
         else
            write(luout,9989)
         endif
         write(luout,9994) rsolv
         write(luout,9988)
      endif
      do iat=1,nat
         dbl_mb(k_rad+iat-1)=abs(dbl_mb(k_rad+iat-1))
         if(dbl_mb(k_rad+iat-1).lt.radtol) then
            dbl_mb(k_rad+iat-1)=zero
         endif
         if(osome.and.me.eq.0) then
           write(luout,9987) iat,dbl_mb(k_coszan+iat-1),
     &                           dbl_mb(k_rad+iat-1)
         endif
      enddo
c
      if(.not.ma_pop_stack(l_coszan))
     & call errquit('cosmo_initialize,ma_pop_stack of l_zan failed',
     & 911,MA_ERR)
c
c     ----- create boundary elements of cavity -----
c
      call hnd_cosset(rtdb,nat,dbl_mb(k_coscoor),dbl_mb(k_rad))
c
c     ----- done -----
c
      if(osome.and.me.eq.0) then
         write(luout,9986)
         write(luout,9985)
      endif
c
      if(me.eq.0) call util_flush(Luout)
c
      if(.not.ma_pop_stack(l_coscoor)) call
     & errquit('cosmo_initialize,ma_pop_stack l_coor failed',
     & 911,MA_ERR)
      if(.not.ma_pop_stack(l_rad)) call
     & errquit('cosmo_initialize,ma_pop_stack l_rad failed',911,MA_ERR)
c
      return
 9999 format(/,10x,15(1h-),
     1       /,10x,'-cosmo- solvent',
     2       /,10x,15(1h-))
 9998 format(' no -cosmo- solvent.')
 9997 format(' dielectric constant -eps-  = ',f6.2)
 9996 format(' screen = (eps-1)/(eps+1/2) = ',f9.5)
 9995 format(' screen = (eps-1)/(eps    ) = ',f9.5)
 9994 format(' solvent radius (ang.)      = ',f7.3)
 9993 format(' -lineq- algorithm          = ',i3)
 9992 format(' -bem- high level           = ',i3)
 9991 format(' -bem- low  level           = ',i3)
 9990 format(' -bem- from -octahedral-')
 9989 format(' -bem- from -icosahedral-')
 9988 format(' atomic radii = ',/,' --------------')
 9987 format(i5,2f7.3)
 9986 format(' ...... end of -cosmo- initialization ......')
 9985 format(/)
 9984 format(' screen =  one              = ',f9.5)
 9983 format(' charge screening approach  = ',i3)
 9982 format(' processing cosmo map from ',a)
 9981 format(a16,4x,f9.5)
      end
c
      subroutine cosmo_tidy()
      implicit none
#include "errquit.fh"
#include "prop.fh"
#include "geom.fh"
c
      if (.not. geom_destroy(cosmo_geom_efc)) call errquit
     $     ('cosmo_tidy: geom destroy failed', 0, GEOM_ERR)
      return
      end
c
      subroutine hnd_cosset(rtdb,nat,c,radius)
      implicit double precision (a-h,o-z)
c
c              ----- starting from -icosahedron- -----
c
c     pass, napex, nface, error =   0      12      20      20
c     pass, napex, nface, error =   1      42      80     100    0.4982
c     pass  napex, nface, error =   2     162     320     420    0.1848
c     pass  napex, nface, error =   3     642    1280    1700    0.0523
c     pass  napex, nface, error =   4    2562    5120    6820    0.0135
c     pass  napex, nface, error =   5   10242   20480   27300    0.0034
c
c              ----- starting from -octahedron-  -----
c
c     pass, napex, nface, error =   0       6       8       8
c     pass, napex, nface, error =   1      18      32      40    0.8075
c     pass  napex, nface, error =   2      66     128     168    0.4557
c     pass  napex, nface, error =   3     258     512     680    0.1619
c     pass  napex, nface, error =   4    1026    2048    2728    0.0451
c     pass  napex, nface, error =   5    4098    8192   10920    0.0116
c     pass  napex, nface, error =   6   16386   32768   43688    0.0029
c
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
c
      integer rtdb
      integer nat
      double precision c(3,nat)
      double precision radius(nat)
c
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
c
      parameter (mxface=10920)
      parameter (mxapex= 4098)
      logical     dbug
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      dbug=.false.
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9999)
      endif
c
      if(ificos.eq.0.and.maxbem.gt.6) then
         write(luout,*) '-maxbem- too large for parameters in -cosset-'
         call errquit('hnd_cosset, -maxbem- too large',911,0)
      elseif(ificos.ne.0.and.maxbem.gt.7) then
         write(luout,*) '-maxbem- too large for parameters in -cosset-'
         call errquit('hnd_cosset, -maxbem- too large',911,0)
      endif
c
c     ----- partition memory -----
c
      need = 6*nat + 7*mxface + 6*mxface*nat + 3*mxapex
c
c     ----- allocate memory block -----
c
      if(.not.ma_push_get(mt_dbl,need,'mem init:cosmo:hnd_cosset:1',
     &    i_init,init))
     & call errquit('hnd_cosset, malloc of init  failed',911,MA_ERR)
c
      i10 =init                    ! xyzatm(3,nat)
      i20 =i10 +3*nat              !   ratm(  nat)
      i30 =i20 +  nat              !   nspa(  nat)
      i40 =i30 +  nat              !   nppa(  nat)
      i50 =i40 +  nat              ! ijkfac(3,mxface)
      i60 =i50 +3*mxface             ! xyzseg(3,mxface)
      i70 =i60 +3*mxface             ! ijkseg(  mxface)
      i80 =i70 +  mxface             ! insseg(  mxface,nat)
      i90 =i80 +  mxface*nat         ! xyzspa(3,mxface,nat)
      i100=i90 +3*mxface*nat         ! ijkspa(  mxface,nat)
      i110=i100+  mxface*nat         ! numpps(  mxface,nat)
      i120=i110+  mxface*nat         ! apex(3,mxapex)
c
c     ----- get -cosmo- surface -----
c
      call hnd_cossrf(nat,c,radius,nat,mxface,mxapex,
     1                x(i10),x(i20),x(i30),x(i40),
     2                x(i50),x(i60),x(i70),
     3                x(i80),x(i90),x(i100),x(i110),x(i120),rtdb)
c
c     ----- release memory block -----
c
      if(.not.ma_pop_stack(i_init)) call
     &  errquit('hnd_cosset, ma_pop_stack of init failed',911,MA_ERR)
c
      return
 9999 format(/,10x,15(1h-),
     1       /,10x,'-cosmo- surface',
     2       /,10x,15(1h-))
      end
c
      subroutine hnd_cossrf(nat,c,radius,mxatm,mxfac,mxapx,
     1                  xyzatm,ratm,nspa,nppa,
     2                  ijkfac,xyzseg,ijkseg,insseg,
     3                  xyzspa,ijkspa,numpps,apex,rtdb)
      implicit double precision (a-h,o-z)
c
#include "nwc_const.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb, nat
      dimension      c(3,nat  )
      dimension radius(    nat)
      dimension xyzatm(3,mxatm)
      dimension   ratm(  mxatm)
      dimension   nspa(  mxatm)
      dimension   nppa(  mxatm)
      dimension ijkfac(3,mxfac)
      dimension xyzseg(3,mxfac)
      dimension ijkseg(  mxfac)
      dimension insseg(  mxfac,mxatm)
      dimension xyzspa(3,mxfac,mxatm)
      dimension ijkspa(  mxfac,mxatm)
      dimension numpps(  mxfac,mxatm)
      dimension   apex(3,mxapx)
c
      parameter (bohr=0.529177249d+00)
      logical    some
      logical    out
      logical    more
      logical    dbug
      logical    insseg
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
c
      dbug=.false.
      more=.false.
      more=more.or.dbug
      out =.false. 
      out =out.or.more
      some=.false.
      some=some.or.out
c
c     ----- approximate sphere with segments and points -----
c
      nseg = 0
      nfac = 0
      ndiv = 0
      call hnd_cossph(nseg,nfac,ndiv,
     1                ijkfac,xyzseg,ijkseg,mxfac,apex,mxapx,
     2                dsurf,dvol,adiag)
c
c     ----- debug printing -----
c
      if(out.and.ga_nodeid().eq.0) then
         write(luout,9999) nseg,nfac,ndiv,dsurf,dvol
         write(luout,9995) adiag
         if(more) then
            write(luout,9998)
            do lseg=1,nseg
               write(luout,9997) lseg,xyzseg(1,lseg),xyzseg(2,lseg),
     1                             xyzseg(3,lseg),ijkseg(  lseg)
            enddo
         endif
         if(dbug) then
            write(luout,9996)
            do lfac=1,nfac
               mfac=lfac+nseg
               write(luout,9997) mfac,xyzseg(1,mfac),xyzseg(2,mfac),
     1                             xyzseg(3,mfac),ijkseg(  mfac)
            enddo
         endif
      endif
c
c     ----- set molecule -----
c
      do iat=1,nat
         do i=1,3
            xyzatm(i,iat)=c(i,iat)
         enddo
      enddo
      do iat=1,nat
         if(radius(iat).eq.0.0d0) then
            ratm(iat)=0.0d0
         else
            ratm(iat)=(radius(iat)+rsolv)/bohr
         endif
      enddo
c
c     ----- create -solvent accessible surface- of the molecule -----
c
      call hnd_cossas(nat,xyzatm,ratm,mxatm,
     1                nspa,nppa,xyzspa,ijkspa,
     2                nseg,nfac,xyzseg,ijkseg,insseg,
     3                numpps,mxfac,rtdb)
c
      return
 9999 format(' nseg,nfac,ndiv=nfac/nseg,dsurf,dvol = ',3i7,2f10.6)
 9998 format('  seg  ','      x     ','      y     ','      z     ',
     1       ' seg ',/,1x,47(1h-))
 9997 format(i7,3f12.8,i5,f12.8)
 9996 format('  fac  ','      x     ','      y     ','      z     ',
     1       ' seg ',/,1x,47(1h-))
 9995 format(' adiag           = ',f12.6)
      end
c
      subroutine hnd_cossas(nat,xyzatm,ratm,mxatom,
     1                  nspa,nppa,xyzspa,ijkspa,
     2                  nseg,nfac,xyzseg,ijkseg,insseg,
     3                  numpps,mxface,rtdb)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb
c
      parameter  (bohr=0.529177249d+00)
      logical     dbug
      logical     insseg
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
      dimension xyzatm(3,mxatom)
      dimension   ratm(  mxatom)
      dimension   nspa(  mxatom)
      dimension   nppa(  mxatom)
      dimension xyzseg(3,mxface)
      dimension ijkseg(  mxface)
      dimension insseg(  mxface,mxatom)
      dimension xyzspa(3,mxface,mxatom)
      dimension ijkspa(  mxface,mxatom)
      dimension numpps(  mxface,mxatom)
      data zero    /0.0d+00/
      integer l_efcc, k_efcc, l_efcs, k_efcs, l_efcz, k_efcz
      integer l_efclb, k_efclb, k_efciat, l_efciat
c
      dist(xi,yi,zi,xj,yj,zj)=sqrt((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)
c
      dbug=.false.
c
      if(ga_nodeid().eq.0) then
         write(luout,9999)
      endif
c
c     ----- print atomic centers -----
c
      if(ga_nodeid().eq.0) then
         write(luout,9998)
         do iat=1,nat
            write(luout,9997) iat,xyzatm(1,iat),xyzatm(2,iat),
     1                                       xyzatm(3,iat),
     2                    (ratm(iat)*bohr-rsolv)
         enddo
      endif
c
c     ----- clear arrays ..... -----
c
      do iat=1,nat
         do i=1,mxface
            ijkspa(i,iat)=0
            numpps(i,iat)=0
         enddo
      enddo
c
c     ----- sift through atomic centers and decide if a face -----
c           belongs to the -sas- or is inside the molecule.
c
      do iat=1,nat
c
         if(ratm(iat).ne.zero) then
            do iseg=1,nseg
               ijkspa(iseg,iat)=ijkseg(iseg)
               do m=1,3
               xyzspa(m,iseg,iat)=xyzseg(m,iseg)*ratm(iat)+xyzatm(m,iat)
               enddo
            enddo
            do ifac=1,nfac
               ijkspa(ifac+nseg,iat)=ijkseg(ifac+nseg)
               do m=1,3
                  xyzspa(m,ifac+nseg,iat)=xyzseg(m,ifac+nseg)*ratm(iat)
     1                                   +xyzatm(m,iat)
               enddo
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9996) iat
               write(luout,9995) (ijkspa(ifac+nseg,iat),ifac=1,nfac)
            endif
            do jat=1,nat
               dij=dist(xyzatm(1,iat),
     1                  xyzatm(2,iat),
     2                  xyzatm(3,iat),
     3                  xyzatm(1,jat),
     4                  xyzatm(2,jat),
     5                  xyzatm(3,jat))
               if(jat.ne.iat.and.(dij.lt.(ratm(iat)+ratm(jat)))) then
                  do ifac=1,nfac
                     dum=dist(xyzspa(1,ifac+nseg,iat),
     1                        xyzspa(2,ifac+nseg,iat),
     2                        xyzspa(3,ifac+nseg,iat),
     3                        xyzatm(1,jat),
     4                        xyzatm(2,jat),
     5                        xyzatm(3,jat))
                     if(dum.lt.ratm(jat)) then
                        ijkspa(ifac+nseg,iat)=0
                     endif
                  enddo
               endif
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9996) iat
               write(luout,9995) (ijkspa(ifac+nseg,iat),ifac=1,nfac)
            endif
c
c     ----- check segments belonging to -sas- -----
c
            do ifac=1,nseg+nfac
               insseg(ifac,iat)=.true.
            enddo
            do ifac=1,nfac
               iseg=ijkspa(ifac+nseg,iat)
               if(iseg.ne.0) then
                  insseg(ifac+nseg,iat)=.false.
                  insseg(     iseg,iat)=.false.
               endif
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9994) iat
               write(luout,9993) (insseg(ifac,iat),ifac=1,nseg+nfac)
            endif
            mseg=0
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) mseg=mseg+1
            enddo
            mfac=0
            do ifac=1,nfac
               if(.not.insseg(ifac+nseg,iat)) mfac=mfac+1
            enddo
            nspa(iat)=mseg
            nppa(iat)=mfac
c
c     ----- surface area of segments -----
c
            do iseg=1,nseg
               numpps(iseg,iat)=0
            enddo
            do ifac=1,nfac
               iseg=ijkspa(ifac+nseg,iat)
               if(iseg.ne.0) numpps(iseg,iat)=numpps(iseg,iat)+1
            enddo
c
         endif
c
      enddo
c
      if(ga_nodeid().eq.0) then
         write(luout,9985) nseg,nfac
         write(luout,9992)
         do iat=1,nat
            npp=0
            do iseg=1,nseg
               npp=npp+numpps(iseg,iat)
            enddo
            write(luout,9991) iat,nspa(iat),nppa(iat),npp
         enddo
      endif
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9987)
         do iat=1,nat
            do iseg=1,nseg
               write(luout,9986) iat,iseg,numpps(iseg,iat)
            enddo
         enddo
      endif
c
c    Count the number of surface points, i.e. number of point charges
c    and generate memory to store them
c
      nefc = 0
      do iat=1,nat
         if(ratm(iat).ne.zero) then
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) nefc = nefc+1
            enddo
         endif
      enddo
c
c     Allocate memory for point charges
c
      if(.not.ma_push_get(mt_dbl,nefc*3,'cosmo efcc',l_efcc,k_efcc))
     & call errquit('cosmo_cossas malloc k_efcc failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcs',l_efcs,k_efcs))
     & call errquit('cosmo_cossas malloc k_efcs failed',911,MA_ERR)
      if(.not.ma_push_get(mt_int,nefc,'cosmo efciat',l_efciat,k_efciat))
     & call errquit('cosmo_cossas malloc k_efciat failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcz',l_efcz,k_efcz))
     & call errquit('cosmo_cossas malloc k_efcz failed',911,MA_ERR)
      if(.not.ma_push_get(mt_byte,nefc*8,'cosmo tags',l_efclb,k_efclb))
     & call errquit('cosmo_cossas malloc k_tag failed',911,MA_ERR)
c
c     ----- save coordinates of surface points -----
c           save segment surfaces
c           save segment to atom mapping
c
      srfmol=zero
      volmol=zero
      ief   =0
      do iat=1,nat
         if(ratm(iat).ne.zero) then
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) then
                  ief=ief+1
                  do i=1,3
                     dbl_mb(k_efcc+3*(ief-1)+i-1)=xyzatm(i,iat)
     1                          +xyzseg(i,iseg)*ratm(iat)
                  enddo
                  ipp=numpps(iseg,iat)
                  dbl_mb(k_efcs+ief-1) = dble(ipp)*dsurf*ratm(iat)**2
                  srfmol   = srfmol + dble(ipp)*dsurf*ratm(iat)**2
                  volmol   = volmol + dble(ipp)*dvol *ratm(iat)**3
                  int_mb(k_efciat+ief-1)=iat
               endif
            enddo
         endif
      enddo
      srfmol=srfmol*(bohr**2)
      volmol=volmol*(bohr**3)
c
      if(ga_nodeid().eq.0) then
         write(luout,9990) nefc
         write(luout,9984) srfmol 
         write(luout,9983) volmol 
      endif 
c
c     ----- Cavity/Dispersion free energy ---
c           Sitkoff, Sharp, and Honig,
c           J.Phys.Chem. 98, 1978 (1994)
c
      cavdsp=0.860+0.005*srfmol
      if(ga_nodeid().eq.0) then
         write(luout,9981) cavdsp
      endif
c
c     ----- print segment surfaces -----
c
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9989)
         do ief=1,nefc
            write(luout,9988) ief,dbl_mb(k_efcs+ief-1),
     &                        int_mb(k_efciat+ief-1)
         enddo
      endif
c
      do ief=1,nefc
         dbl_mb(k_efcz+ief-1)=zero
      enddo
      do ief=1,nefc
         byte_mb(k_efclb+(ief-1)*8)='        '
      enddo
c
c     ----- write out to -rtdb- -----
c
      if(.not.rtdb_put(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     $   call errquit('hnd_cossas: rtdb put failed for nefc  ',911,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcc',mt_dbl,3*nefc,dbl_mb(k_efcc)))
     $   call errquit('hnd_cossas: rtdb put failed for efcc  ',912,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcz',mt_dbl,  nefc,dbl_mb(k_efcz)))
     $   call errquit('hnd_cossas: rtdb put failed for efcz  ',913,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcs',mt_dbl,  nefc,dbl_mb(k_efcs)))
     $   call errquit('hnd_cossas: rtdb put failed for efcs  ',914,
     &       rtdb_err)
c     if(.not.rtdb_put(rtdb,'cosmoiatefc',mt_int,nefc,int_mb(k_efciat)))
c    $   call errquit('hnd_cossas: rtdb put failed for iatefc',915,
c    &       rtdb_err)
c     if(.not.rtdb_cput(rtdb,'char variable',nefc,byte_mb(k_efclb)))
c    $   call errquit('hnd_cossas: rtdb put failed for efclab',916,
c    &       rtdb_err)
c
      if(.not.ma_pop_stack(l_efclb)) call
     &   errquit('cosmo_cossas chop stack k_efclb failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efcz)) call
     &   errquit('cosmo_cossas chop stack k_efcz failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efciat)) call
     &   errquit('cosmo_cossas chop stack k_efciat failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efcs)) call
     &   errquit('cosmo_cossas chop stack k_efcs failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efcc)) call
     &   errquit('cosmo_cossas chop stack k_efcc failed',911,MA_ERR)
c
      return
 9999 format(/,1x,'solvent accessible surface',/,1x,26(1h-))
 9998 format(/,1x,'---------- ATOMIC COORDINATES (A.U.) ----------',
     1            '-- VDWR(ANG.) --')
 9997 format(  1x,i5,3f14.8,f10.3)
 9996 format(/,1x,'---------- SEGMENTS FOR -IAT- = ',i5)
 9995 format(16i4)
 9994 format(/,1x,'-INSSEG- FACES FOR IAT = ',i5)
 9993 format(16l4)
 9992 format(  1x,'atom',' ( ','  nspa',',','  nppa',' )',/,1x,22(1h-))
 9991 format(  1x,   i4 ,' ( ',     i6 ,',',     i6 ,' )',i8)
 9990 format(  1x,'number of -cosmo- surface points = ',i8)
 9989 format(  1x,'SEGMENT SURFACES =',/,1x,18(1h-))
 9988 format(i8,f10.5,i5)
 9987 format(  1x,'NUMBER OF FACES / SEGMENT =',/,1x,27(1h-))
 9986 format(3i5)
 9985 format(' number of segments per atom = ',i10,/,
     1       ' number of   points per atom = ',i10)
 9984 format(' molecular surface = ',f10.3,' angstrom**2')
 9983 format(' molecular volume  = ',f10.3,' angstrom**3')
 9981 format(' G(cav/disp)       = ',f10.3,' kcal/mol')
      enD
c
      subroutine hnd_cossph(nseg,nfac,ndiv,
     1                  ijkfac,xyzseg,ijkseg,mxface,apex,mxapex,
     2                  dsurf,dvol,adiag)
      implicit double precision (a-h,o-z)
#include "global.fh"
#include "stdio.fh"
c
c              ----- starting from -icosahedron- -----
c
c     pass, napex, nface, error =   0      12      20      20
c     pass, napex, nface, error =   1      42      80     100    0.4982
c     pass  napex, nface, error =   2     162     320     420    0.1848
c     pass  napex, nface, error =   3     642    1280    1700    0.0523
c     pass  napex, nface, error =   4    2562    5120    6820    0.0135
c     pass  napex, nface, error =   5   10242   20480   27300    0.0034
c
c              ----- starting from -octahedron-  -----
c
c     pass, napex, nface, error =   0       6       8       8
c     pass, napex, nface, error =   1      18      32      40    0.8075
c     pass  napex, nface, error =   2      66     128     168    0.4557
c     pass  napex, nface, error =   3     258     512     680    0.1619
c     pass  napex, nface, error =   4    1026    2048    2728    0.0451
c     pass  napex, nface, error =   5    4098    8192   10920    0.0116
c     pass  napex, nface, error =   6   16386   32768   43688    0.0029
c
      dimension   apex(3,*)
      dimension ijkfac(3,*)
      dimension ijkseg(  *)
      dimension xyzseg(3,*)
c
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
c
      parameter (mxpass=    7)
      dimension minfac(mxpass)
      dimension maxfac(mxpass)
      dimension minico(mxpass)
      dimension maxico(mxpass)
      dimension minoct(mxpass)
      dimension maxoct(mxpass)
      dimension ijknew(3)
      dimension ijkold(3)
      equivalence (ijkold(1),iold),(ijkold(2),jold),(ijkold(3),kold)
      equivalence (ijknew(1),inew),(ijknew(2),jnew),(ijknew(3),knew)
      logical icos
      logical octa
      logical some,out,dbug
      data minico /    1,   21,  101,  421, 1701, 6821,    0/
      data maxico /   20,  100,  420, 1700, 6820,27300,    0/
      data minoct /    1,    9,   41,  169,  681, 2729,10921/
      data maxoct /    8,   40,  168,  680, 2728,10920,43688/
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data two   /2.0d+00/
      data three /3.0d+00/
      data four  /4.0d+00/
c
      dist(xi,yi,zi,xj,yj,zj)=sqrt((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.
      some=some.or.out
c
      pi=four*atan(one)
      rad=one
      cir= two*pi*rad
      srf=four*pi*rad**2
      vol=four*pi*rad**3/three
c
      npass=maxbem
c
c     ----- define  hedron  -----
c           define -minfac- 
c           define -maxfac- 
c
      icos=ificos.ne.0
      octa=.not.icos
      if(icos) then
         call hnd_sphico(apex,napex,ijkfac,ijkseg,nface)
         do ipass=1,mxpass
            minfac(ipass)=minico(ipass)
            maxfac(ipass)=maxico(ipass)
         enddo
      endif
      if(octa) then
         call hnd_sphoct(apex,napex,ijkfac,ijkseg,nface)
         do ipass=1,mxpass
            minfac(ipass)=minoct(ipass)
            maxfac(ipass)=maxoct(ipass)
         enddo
      endif
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         if(icos) then
            write(luout,9994)
         endif
         if(octa) then
            write(luout,9982)
         endif
         if(out) then
            write(luout,*) '-minbem- = ',minbem
            write(luout,*) '-maxbem- = ',maxbem
            write(luout,*) '-minfac- = ',minfac
            write(luout,*) '-maxfac- = ',maxfac
            write(luout,*) '-npass - = ',npass
            write(luout,9999)
            do iapex=1,napex
               write(luout,9998) iapex,apex(1,iapex),
     1                              apex(2,iapex),
     2                              apex(3,iapex)
            enddo
         endif
      endif
c
c     ----- loop over divisions to create sphere -----
c
      mxfac=0
      ipass=1
  100 ipass=ipass+1
         mnfac=mxfac+1
         mxfac=nface
         if(out.and.ga_nodeid().eq.0) then
            write(luout,9996) ipass,napex,nface,mnfac,mxfac
         endif
c
         dmin =one
         mapex=napex
         mface=nface
         do lface=mnfac,mxfac
            iold=ijkfac(1,lface)
            jold=ijkfac(2,lface)
            kold=ijkfac(3,lface)
            call hnd_sphapx(apex,mapex,ijkfac,ijkseg,mface,lface,
     1                      ijkold,ijknew,dijk)
            dmin=min(dmin,dijk)
         enddo
         napex=mapex
         nface=mface
         if(out.and.ga_nodeid().eq.0) then
            write(luout,9995) napex,nface
         endif
c
c     ----- print out new apeces -----
c
         if(dbug.and.ga_nodeid().eq.0) then
            do iapex=1,napex
               write(luout,9998) iapex,apex(1,iapex),apex(2,iapex),
     1                              apex(3,iapex)
            enddo
         endif
c
c     ----- print approximate volume -----
c
         radapp=    dmin
         radrat=    dmin
         raderr=one-radrat
         srfapp=srf*dmin**2
         srfrat=    dmin**2
         srferr=one-srfrat
         volapp=vol*dmin**3
         volrat=    dmin**3
         volerr=one-volrat
         if(out.and.ga_nodeid().eq.0) then
            write(luout,9997) vol,volapp,volrat,volerr
            write(luout,9992) srf,srfapp,srfrat,srferr
            write(luout,9991) rad,radapp,radrat,raderr
         endif
c
c     ----- assign early segment to each face -----
c
         if(ipass.gt.(minbem+1)) then
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9981) ipass
               write(luout,9980) (minfac(i),i=1,ipass)
               write(luout,9979) (maxfac(i),i=1,ipass)
            endif
            maxseg=maxfac(minbem)
            lfacmn=minfac(ipass)
            lfacmx=maxfac(ipass)
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9990) ipass
               write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
            endif
            do lface=lfacmn,lfacmx
               ijkseg(lface)=ijkseg(ijkseg(lface))
               if(ijkseg(lface).gt.maxseg.and.ga_nodeid().eq.0) then
                  write(luout,9987) lface,ijkseg(lface)
               endif
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9989) ipass
               write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
            endif
         endif
c
      if(ipass.lt.npass) go to 100
c
c     ----- end of loop over tessalating passes -----
c
      if(dbug.and.ga_nodeid().eq.0) then
         do ipass=1,npass
            lfacmn=minfac(ipass)
            lfacmx=maxfac(ipass)
            write(luout,9989) ipass
            write(luout,*) '-lfacmn- = ',lfacmn
            write(luout,*) '-lfacmx- = ',lfacmx
            write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
         enddo
      endif
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9993) npass,napex,minfac(npass),maxfac(npass),
     1                  radapp,raderr,srfapp,srferr,volapp,volerr
      endif
c
c     ----- at this point each of the faces is assigned to one -----
c           segment. now define centers of segments ...
c
      third =one/three
      lfacmn= minfac(minbem)
      lfacmx= maxfac(minbem)
      do lface=lfacmn,lfacmx
         mface=lface-lfacmn+1
         ijkseg(mface)=mface
         i=ijkfac(1,lface)
         j=ijkfac(2,lface)
         k=ijkfac(3,lface)
         do m=1,3
            xyzseg(m,mface)=(apex(m,i)+apex(m,j)+apex(m,k))*third
         enddo
         dseg=one/dist(xyzseg(1,mface),xyzseg(2,mface),xyzseg(3,mface),
     1                 zero,zero,zero)
         do m=1,3
            xyzseg(m,mface)=xyzseg(m,mface)*dseg
         enddo
      enddo
      nseg=(lfacmx-lfacmn+1)
c
      if(dbug.and.ga_nodeid().eq.0) then
         lfacmn=1
         lfacmx=nseg
         write(luout,*)    'segment to segment mapping = '
         write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
      endif
c
c     ----- now the faces ... -----
c
      if(npass.gt.minbem) then
         lfacmn=minfac(minbem+1)
         lfacmx=maxfac(npass   )
         do lface=lfacmn,lfacmx
            mface=lface-lfacmn+1    
     1                        +(maxfac(minbem)-minfac(minbem)+1)
            ijkseg(mface)=ijkseg(lface)
     1                        -(               minfac(minbem)-1)
            i=ijkfac(1,lface)
            j=ijkfac(2,lface)
            k=ijkfac(3,lface)
            do m=1,3
               xyzseg(m,mface)=(apex(m,i)+apex(m,j)+apex(m,k))*third
            enddo
            dseg=one/dist(xyzseg(1,mface),
     1                    xyzseg(2,mface),
     2                    xyzseg(3,mface),zero,zero,zero)
            do m=1,3
               xyzseg(m,mface)=xyzseg(m,mface)*dseg
            enddo
         enddo
         nfac=(lfacmx-lfacmn+1)
      else
         do iseg=1,nseg
            ifac=iseg+nseg
            ijkseg(ifac)=ijkseg(iseg)
            do m=1,3
               xyzseg(m,ifac)=xyzseg(m,iseg)
            enddo
         enddo
         nfac=nseg
      endif
c
      if(dbug.and.ga_nodeid().eq.0) then
         lfacmn=nseg+1
         lfacmx=nseg+nfac
         write(luout,*)    ' face   to segment mapping = '
         write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
      endif
c
c     ----- calculate -dsurf dvol- for the -cosmo- theory -----
c
      ndiv =nfac/nseg
      dsurf=srf/dble(nfac)
      dvol =vol/dble(nfac)
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9986) nseg,nfac,ndiv,dsurf,dvol
      endif
      if(out.and.ga_nodeid().eq.0) then
         write(luout,9985)
         do i=1,nseg
            done=dist(xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     1                zero,zero,zero)
            write(luout,9984) i,
     1                     xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     2                     ijkseg(i),done
         enddo
      endif
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9985)
         do i=nseg+1,nseg+nfac
            done=dist(xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     1                zero,zero,zero)
            write(luout,9984) (i-nseg),
     1                     xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     2                     ijkseg(i),done
         enddo
      endif
c
c     ----- calculate -adiag- of the -cosmo- theory -----
c
      avgdia=zero
      avgfac=zero
      do mseg=1,nseg
         sum=zero
         do lseg=1,nseg
            if(lseg.ne.mseg) then
               l1=mseg
               l2=lseg
         sum=sum+rad/dist(xyzseg(1,l2),xyzseg(2,l2),xyzseg(3,l2),
     1                    xyzseg(1,l1),xyzseg(2,l1),xyzseg(3,l1))
            endif
         enddo
         fac=(dble(nseg)-sum)/sqrt(dble(nseg))
         adiag=sqrt(four*pi)*fac
         if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
            write(luout,9983) mseg,adiag,fac,dble(nseg),sum
         endif
         avgdia=avgdia+adiag
         avgfac=avgfac+fac
      enddo
      adiag=avgdia/dble(nseg)
      fac  =avgfac/dble(nseg)
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9978)      adiag,fac               
      endif
c
      return
 9999 format(/,1x,'  apex',5x,'x',6x,5x,'y',6x,5x,'z',6x,/,1x,42(1h-))
 9998 format(1x,i6,f12.8,f12.8,f12.8)
 9997 format(' vol, approx., ratio, error = ',2f12.8,2 f8.4)
 9996 format(' pass, napex, nface, mnfac, mxfac = ',i3,4i8)
 9995 format('       napex, nface               = ',3x,2i8)
 9994 format(1x,'sphere from -icosahedron-',/,1x,25(1h-))
 9993 format(' npass = ',i2,' napex = ',i8,
     1       ' minfac = ',i8,' maxfac = ',i8,/,
     2       ' rad = ',f10.6,' error = ',f8.4,/,
     3       ' srf = ',f10.6,' error = ',f8.4,/,
     4       ' vol = ',f10.6,' error = ',f8.4)
 9992 format(' srf, approx., ratio, error = ',2f12.8,2 f8.4)
 9991 format(' rad, approx., ratio, error = ',2f12.8,2 f8.4)
 9990 format(' absolute -ijkseg- , for -ipass- = ',i5)
 9989 format(' relative -ijkseg- , for -ipass- = ',i5)
 9988 format(12i6)
 9987 format(' assigned segment for -lface- = ',i7,
     1       ' is = ',i7,' ( greater than -maxseg- = ',i4,' )')
 9986 format(' nseg,nfac,ndiv=nfac/nseg,dsurf,dvol = ',3i7,2f10.6)
 9985 format('   pt  ','      x     ','      y     ','      z     ',
     1       ' seg ','    norm    ',/,1x,59(1h-))
 9984 format(i7,3f12.8,i5,f12.8)
 9983 format(' mseg,adiag,fac,m,sum = ',i7,4f12.6)
 9982 format(1x,'sphere from -octahedron-',/,1x,24(1h-))
 9981 format(' pass # = ',i5)
 9980 format(' minfac = ',10i5)
 9979 format(' maxfac = ',10i5)
 9978 format('      adiag,fac       = ',   2f12.6)
      end

      subroutine hnd_sphoct(apex,napex,ijkfac,ijkseg,nface)
      implicit double precision (a-h,o-z)
#include "global.fh"
#include "stdio.fh"
c
      logical out
      dimension    xyz(3,6)
      dimension    ijk(3,8)
      dimension   apex(3,*)
      dimension ijkfac(3,*)
      dimension ijkseg(  *)
      data xyz / 1.0d+00, 0.0d+00, 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00,
     1          -1.0d+00, 0.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00,
     2           0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00/
      data ijk / 5, 1, 2, 5, 2, 3, 5, 3, 4, 5, 4, 1,
     1           6, 1, 2, 6, 2, 3, 6, 3, 4, 6, 4, 1/
c
      out=.false.
      out=out.and.ga_nodeid().eq.0
c
      if(out) then
         write(luout,9997)
      endif
c
c     ----- set the 6 apeces of an octahedron -----
c
c     1     1.     0.     0.
c     2     0.     1.     0.
c     3    -1.     0.     0.
c     4     0.    -1.     0.
c     5     0.     0.     1.
c     6     0.     0.    -1.
c
      napex=6
      do iapex=1,napex
         apex(1,iapex)=xyz(1,iapex)
         apex(2,iapex)=xyz(2,iapex)
         apex(3,iapex)=xyz(3,iapex)
      enddo
      if(out) then
         write(luout,9999)
         do iapex=1,napex
            write(luout,9998) iapex,apex(1,iapex),apex(2,iapex),
     1                           apex(3,iapex)
         enddo
      endif
c
      nface=8
      do lface=1,nface
         ijkfac(1,lface)=ijk(1,lface)
         ijkfac(2,lface)=ijk(2,lface)
         ijkfac(3,lface)=ijk(3,lface)
         ijkseg(  lface)=      lface
      enddo
c
      if(out) then
         write(luout,*) '...... end of -sphoct- ......'
      endif
      return
 9999 format(/,1x,'  apex',5x,'x',6x,5x,'y',6x,5x,'z',6x,/,1x,42(1h-))
 9998 format(1x,i6,f12.8,f12.8,f12.8)
 9997 format(/,1x,'octahedron',/,1x,10(1h-))
      end
c
      subroutine hnd_sphico(apex,napex,ijkfac,ijkseg,nface)
      implicit double precision (a-h,o-z)
#include "global.fh"
#include "stdio.fh"
c
      logical out
      dimension      c(3,12)
      dimension      s(3,12)
      dimension    ijk(3,20)
      dimension   apex(3,*)
      dimension ijkfac(3,*)
      dimension ijkseg(  *)
      data c   / 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00,
     1           0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00,
     2           0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00, 1.0d+00,
     3           0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00,
     4           1.0d+00, 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,
     5          -1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00/
      data s   / 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00,
     1           0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00,
     2           1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,
     3           1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,
     4           0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00,
     5           0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00/
      data ijk / 1, 2, 9, 1, 9, 5, 1, 5, 6, 1, 6,11, 1,11, 2,
     1                    2, 9, 7, 2, 7, 8, 2, 8,11,
     2           3, 4,10, 3,10, 5, 3, 5, 6, 3, 6,12, 3,12, 4,
     3                    4,10, 7, 4, 7, 8, 4, 8,12,
     4           9,10, 7, 9, 5,10,11,12, 8,11, 6,12/
      data one  /1.0d+00/
      data two  /2.0d+00/
      data five /5.0d+00/
c
      out=.false.
      out=out.and.ga_nodeid().eq.0
c
      if(out) then
         write(luout,9997)
      endif
c
c     ----- set the 12 apeces of an icosahedron -----
c
c     1     0.     cosa   sina
c     2     0.     cosa  -sina
c     3     0.    -cosa   sina
c     4     0.    -cosa  -sina
c     5     sina   0.     cosa
c     6    -sina   0.     cosa
c     7     sina   0.    -cosa
c     8    -sina   0.    -cosa
c     9     cosa   sina   0.
c    10     cosa  -sina   0.
c    11    -cosa   sina   0.
c    12    -cosa  -sina   0.
c
      ang=acos(one/sqrt(five))/two
      cosa=cos(ang)
      sina=sin(ang)
      napex=12
      do iapex=1,napex
         apex(1,iapex)=cosa*c(1,iapex)+sina*s(1,iapex)
         apex(2,iapex)=cosa*c(2,iapex)+sina*s(2,iapex)
         apex(3,iapex)=cosa*c(3,iapex)+sina*s(3,iapex)
      enddo
      if(out) then
         write(luout,9999)
         do iapex=1,napex
            write(luout,9998) iapex,apex(1,iapex),apex(2,iapex),
     1                           apex(3,iapex)
         enddo
      endif
c
      nface=20
      do lface=1,nface
         ijkfac(1,lface)=ijk(1,lface)
         ijkfac(2,lface)=ijk(2,lface)
         ijkfac(3,lface)=ijk(3,lface)
         ijkseg(  lface)=      lface
      enddo
c
      if(out) then
         write(luout,*) '...... end of -sphico- ......'
      endif
      return
 9999 format(/,1x,'  apex',5x,'x',6x,5x,'y',6x,5x,'z',6x,/,1x,42(1h-))
 9998 format(1x,i6,f12.8,f12.8,f12.8)
 9997 format(/,1x,'icosahedron',/,1x,11(1h-))
      end
c
      subroutine hnd_sphapx(apex,mapex,ijkfac,ijkseg,mface,lface,
     1                             ijkold,ijknew,dmin)
      implicit double precision (a-h,o-z)
#include "global.fh"
#include "stdio.fh"
c
      logical out
      logical duplic
      dimension   apex(3,*)
      dimension ijkfac(3,*)
      dimension ijkseg(  *)
      dimension ijkold(3)
      dimension ijknew(3)
      dimension    xyz(3,3)
      dimension      d(3)
      dimension xyzijk(3)
      equivalence (xyz(1,1),xij),(xyz(2,1),yij),(xyz(3,1),zij),
     1            (xyz(1,2),xjk),(xyz(2,2),yjk),(xyz(3,2),zjk),
     2            (xyz(1,3),xki),(xyz(2,3),yki),(xyz(3,3),zki)
      data zero  /0.0d+00/
      data pt5   /0.5d+00/
      data one   /1.0d+00/
      data three /3.0d+00/
      data tol   /1.0d-04/
c
      dist(x1,y1,z1,x2,y2,z2)=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
c
      out=.false.
      out=out.and.ga_nodeid().eq.0
c
c     ----- create mid-point of the 3 edges -----
c
      iold=ijkold(1)
      jold=ijkold(2)
      kold=ijkold(3)
      do m=1,3
         xyz(m,1)=(apex(m,iold)+apex(m,jold))*pt5
         xyz(m,2)=(apex(m,jold)+apex(m,kold))*pt5
         xyz(m,3)=(apex(m,kold)+apex(m,iold))*pt5
      enddo
c
c     ----- project onto sphere -----
c
      d(1)=dist(xij,yij,zij,zero,zero,zero)
      d(2)=dist(xjk,yjk,zjk,zero,zero,zero)
      d(3)=dist(xki,yki,zki,zero,zero,zero)
      d(1)=one/d(1)
      d(2)=one/d(2)
      d(3)=one/d(3)
      do l=1,3
         do m=1,3
            xyz(m,l)=xyz(m,l)*d(l)
         enddo
      enddo
c
c     ----- check for duplicate apeces -----
c
      newapx=0
      do iapx=1,3
         duplic=.false.
         do lapex=1,mapex
            dd=dist(xyz(1,  iapx),xyz(2,  iapx),xyz(3,  iapx),
     1              apex(1,lapex),apex(2,lapex),apex(3,lapex))
            if(abs(dd).lt.tol) then
               duplic=.true.
               lduplc=lapex
            endif
         enddo
         if(duplic) then
            ijknew(iapx)=lduplc
            if(out) then
               write(luout,9999) iapx,ijkold,lduplc
            endif
         else
            newapx=newapx+1
            japx=mapex+newapx
            ijknew(iapx)=japx
            do m=1,3
               apex(m,japx)=xyz(m,iapx)
            enddo
            if(out) then
               write(luout,9998) iapx,ijkold,japx,
     1                        apex(1,japx),apex(2,japx),apex(3,japx)
            endif
         endif
      enddo
      mapex=mapex+newapx
c
c     ----- make up new faces and their centers -----
c
      third=one/three
      dmin =one
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijkold(1)
      ijkfac(2,mface)=ijknew(1)
      ijkfac(3,mface)=ijknew(3)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijkold(2)
      ijkfac(2,mface)=ijknew(1)
      ijkfac(3,mface)=ijknew(2)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijkold(3)
      ijkfac(2,mface)=ijknew(2)
      ijkfac(3,mface)=ijknew(3)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijknew(1)
      ijkfac(2,mface)=ijknew(2)
      ijkfac(3,mface)=ijknew(3)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      if(out) then
         write(luout,9997) dmin,mface
      endif
c
      return
 9999 format(' duplicated apex =',i2,' for face ',3i5,'. same as = ',i5)
 9998 format('    new     apex =',i2,' for face ',3i5,'.  newapx = ',i5,
     1       /,7x,3f12.8)
 9997 format(' --- dmin = ',f12.8,' --- mface = ',i10)
      end
c
      subroutine cosmo_charges(rtdb,basis,geom,ecos,some,itype_wfn,
     1                         vectors,nclosed,nopen,nbf,nmo)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "prop.fh"
c
      integer rtdb
      integer basis
      integer geom, itype_wfn
      double precision ecos
      integer vectors(2), nclosed(2), nopen(2)

      logical dbug, out, some
      character*8 scfruntyp
c
      common/hnd_cospar/dielec,screen,rsolv
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
c
      integer l_efcc, k_efcc, l_efcs, k_efcs, l_efcz, k_efcz
      integer g_dens(3),ndens,nvirt(2)
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
c
      dbug=.false..and.ga_nodeid().eq.0
      out =.false..and.ga_nodeid().eq.0
      out =out.or.dbug
      some=some.or.out.and.ga_nodeid().eq.0
c
      if(some) then
         write(luout,9999)
      endif
c
      if(out) then
         write(luout,*) 'in cosmo_charges ... geom = ',geom
      endif
c
c     ----- read -efc- coordinates from -rtdb- -----
c
      if(.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     &   call errquit('cosmo_charges: rtdb get failed for nefc  ',911,
     &       RTDB_ERR)
c
c     ----- allocate memory for efc's -----
c
      if(.not.ma_push_get(mt_dbl,nefc*3,'cosmo efcc',l_efcc,k_efcc))
     & call errquit('cosmo_cossas malloc k_efcc failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcs',l_efcs,k_efcs))
     & call errquit('cosmo_cossas malloc k_efcs failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcz',l_efcz,k_efcz))
     & call errquit('cosmo_cossas malloc k_efcz failed',911,MA_ERR)
c
      if(.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,dbl_mb(k_efcc)))
     &   call errquit('cosmo_charges: rtdb get failed for efcc  ',912,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,  nefc,dbl_mb(k_efcz)))
     &   call errquit('cosmo_charges: rtdb get failed for efcz  ',913,
     &       RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cosmo:efcs',mt_dbl,  nefc,dbl_mb(k_efcs)))
     &   call errquit('cosmo_charges: rtdb get failed for efcs  ',914,
     &       RTDB_ERR)
c
      if(out) then
         write(luout,*) 'in cosmo_charges, nefc = ',nefc
         do ief=1,nefc
            write(luout,*) dbl_mb(k_efcc+(ief-1)*3), 
     &            dbl_mb(k_efcc+(ief-1)*3+1),dbl_mb(k_efcc+(ief-1)*3+2)
         enddo
      endif
cc
cc     ----- get density matrix -----
cc
c      scfruntyp='RHF'
c      nocc=nclosed(1)+nopen(1)
c      if (itype_wfn.eq.2) then
c         scfruntyp='UHF'
c         nocc=max(nocc,nclosed(2)+nopen(2))
c      endif
c      if(.not.ma_push_get(mt_dbl,nocc*2,'cosmo occ',l_occ,k_occ))
c     &   call errquit('cosmo_charges malloc k_occ failed',911,MA_ERR)
c      do i=1,nocc*2
c         dbl_mb(k_occ+i-1)=1.0d0
c      enddo

c
c     changes to make cosmo work for open shell DFT (MV)
c
      if(.not.ma_push_get(mt_dbl,nbf*2,'cosmo occ',l_occ,k_occ))
     &   call errquit('cosmo_charges malloc k_occ failed',911,MA_ERR)
         call dfill(2*nbf, 0.0d0, dbl_mb(k_occ), 1)
c
c     Only need to set occupation numbers for UHF
c     occupation numbers for RHF are done inside hnd_prop_dens_make
c
      if (itype_wfn.eq.2) then
         scfruntyp='UHF'
         do i = 1, nopen(1)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
         do i = nbf+1, nbf+nopen(2)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
      else if (itype_wfn.eq.1) then
         scfruntyp='RHF'
      else
         call errquit("unknown function type",0,0)
      endif
c
c     end of changes (MV)
c

      call hnd_prop_dens_make(rtdb,geom,basis,nbf,nmo,nclosed,nopen,
     &                       nvirt,scfruntyp,vectors,dbl_mb(k_occ),
     &                       g_dens,ndens)
c
      if(out) then
         write(luout,*) 'in -cosmo_charges, wfntyp = ',itype_wfn
         write(luout,*) 'in -cosmo_charges, scftyp = ',scfruntyp
         if(dbug) then
            write(luout,*) 'in -cosmo_charges, -da- ...'
            call ga_print(g_dens(1))
            if(scfruntyp.eq.'UHF     ') then
               write(luout,*) 'in -cosmo_charges, -db- ...'
               call ga_print(g_dens(2))
            endif
         endif
      endif
c
c     ----- get -cosmo- charges -----
c
      call hnd_coschg(g_dens,ndens,rtdb,geom,basis,nefc,dbl_mb(k_efcc),
     &                dbl_mb(k_efcs),dbl_mb(k_efcz),ecos)
c
c     ----- release memory block -----
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('hnd_coschg: ga_destroy failed g_dens',0,GA_ERR)
      enddo
      if(.not.ma_chop_stack(l_efcc))
     & call errquit('cosmo_charges, ma_chop_stack of l_efcc failed',911,
     &       ma_err)
c
      return
 9999 format(/,10X,15(1H-),
     1       /,10X,'-cosmo- charges',
     2       /,10X,15(1H-))
      end
c
      subroutine hnd_coschg(g_dens,ndens,rtdb,geom,basis,nefc,
     &                      efcc,efcs,efcz,ecos)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "stdio.fh"
#include "prop.fh"
#include "util.fh"
c
      integer  rtdb
      integer  geom
      integer  basis
      logical  status
c
      logical  dbug,more,out,direct,noall,all,elec,nucl,iefc_done
      character*16 at_tag
      character*255 cosmo_file
      integer fn
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
      common/hnd_cospar/dielec,screen,rsolv
      dimension efcc(3,nefc),efcz(nefc),efcs(nefc)
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      integer g_dens(ndens)
c
      logical util_io_unit
      external util_io_unit
c
      more=.false.
      dbug=.false.
      dbug=dbug.or.more
      out =.false.
      out =out.or.dbug
c
      dbug=dbug.and.ga_nodeid().eq.0
      more=more.and.ga_nodeid().eq.0
      out =out .and.ga_nodeid().eq.0
c
      iefc_done=.false.
      ecos=zero
c
c     ----- allocate memory block -----
c
      if(lineq.ne.1) then
         need=11*nefc+nefc*nefc
         if(.not.ma_push_get(mt_dbl,need,
     &      'mem init:cosmo:hnd_coschg:1',i_init,init)) then
            need=nefc+nefc+nefc+nefc+nefc+nefc+nefc
            lineq=1
            if(.not.ma_push_get(mt_dbl,need,
     &         'mem init:cosmo:hnd_coschg:1',i_init,init))
     &         call errquit(
     &              'hnd_coschg, malloc of init  failed',911,MA_ERR)
            if(out) then
               write(luout,*) 'use iterative method for -lineq-'
            endif
         else
            if(out) then
               write(luout,*) 'use in-memory method for -lineq-'
            endif
         endif
      else
         need=11*nefc
         if(.not.ma_push_get(mt_dbl,need,
     &      'mem init:cosmo:hnd_coschg:2',i_init,init))
     &    call errquit('hnd_coschg, malloc of init  failed',912,MA_ERR)
         lineq=1
         if(out) then
            write(luout,*) 'use iterative method for -lineq-'
         endif
      endif
c
c     ----- get electrostatic potential at surface points -----
c
      i10=init
      i11=i10+nefc
      i12=i11+nefc
c
c     ----- calculate electronic contribution at all points -----
c
      call hnd_elfcon(basis,geom,g_dens(ndens),efcc,nefc,x(i12),0)
c
c     ----- nuclear contribution -----
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_coschg: geom_ncent',911,GEOM_ERR)
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'epot pnt',l_epot,k_epot))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_coschg: geom_cent_get',911,GEOM_ERR)
   30 continue
c
      do 230  ipt=1,nefc
         xp = efcc(1,ipt)
         yp = efcc(2,ipt)
         zp = efcc(3,ipt)
         x(i11+ipt-1) = 0.0d0
         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            x(i11+ipt-1) = x(i11+ipt-1) + zan/rr
  210    continue
         x(i12+ipt-1) = -x(i12+ipt-1)
         x(i10+ipt-1) = x(i11+ipt-1) + x(i12+ipt-1)
  230 continue
c
c     ----- ( nucl, elec, all ) passes ... -----
c
      if(ifscrn.ne.2) then
         noall=.true.
      else
         noall=.false.
      endif

      if(noall) then
         all =.true.
         elec=.false.
         nucl=.false.
      else
         all  =.false.
         elec =.false.
         nucl =.true.
      endif
c
   10 continue

      if(out) then
         write(luout,*) 'in -coschg- ... noall, all, elec, nucl = ',
     1                                noall, all, elec, nucl
      endif
c
c     ----- get surface charges -----
c
      if(lineq.eq.0) then
         i10=init
         i11=i10+nefc
         i12=i11+nefc
         i20=i12+nefc
         i21=i20+nefc
         i22=i21+nefc
         i30=i22+nefc
         i40=i30+nefc*nefc
         i50=i40+nefc
         i60=i50+nefc
c
c     ----- in memory -----
c
         if(out) then
            write(luout,*) '-lineq- in memory'
         endif

         if(all) then
            iep=i10
         elseif(nucl) then
            iep=i11
         elseif(elec) then
            iep=i12
         else
            write(luout,*)'... should not be here ... all,elec,nucl = ',
     1                                              all,elec,nucl
            call errquit('hnd_coschg: should not be here...',911,UERR)
         endif

         do ief=1,nefc
            x(ief+i20-1)=x(ief+iep-1)
         enddo

         do jef=1,nefc
            xj=efcc(1,jef)
            yj=efcc(2,jef)
            zj=efcc(3,jef)
            do ief=1,nefc
               if(ief.eq.jef) then
                  aii=adiag/sqrt(efcs(ief))
                  x(ief+nefc*(jef-1)+i30-1)=aii
               else
                  xi=efcc(1,ief)
                  yi=efcc(2,ief)
                  zi=efcc(3,ief)
                  dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
                  x(ief+nefc*(jef-1)+i30-1)=one/dij
               endif
            enddo
         enddo
c
         nodcmp=1
         call hnd_linequ(x(i30),nefc,x(i20),nefc,
     1                   x(i40),x(i50),deta,ierr,nodcmp)
c
      else
         i10=init
         i11=i10+nefc
         i12=i11+nefc
         i20=i12+nefc
         i21=i20+nefc
         i22=i21+nefc
         i30=i22+nefc
         i40=i30+nefc
         i50=i40+nefc
         i60=i50+nefc
         i70=i60+nefc
         i80=i70+nefc
c
c     ----- iterative -----
c
         if(out) then
            write(luout,*) 'iterative -lineq-'
         endif
c
         do ief=1,nefc
            x(ief+i20-1)=zero
         enddo
c
         direct=.true.
c
c     ----- solve ... -----
c
         if(all) then
            iep=i10
         elseif(nucl) then
            iep=i11
         elseif(elec) then
            iep=i12
         endif

         call hnd_cosequ(x(iep),x(i20),x(i30),nefc,
     1                   x(i40),x(i50),x(i60),x(i70),
     2                   efcc,efcs)
      endif
c
c     ----- debug check -----
c
      if(more) then
         do jef=1,nefc
            xj=efcc(1,jef)
            yj=efcc(2,jef)
            zj=efcc(3,jef)
            dum=zero
            do ief=1,nefc
               if(ief.eq.jef) then
                  aii=adiag/sqrt(efcs(ief))
                  dum=dum+aii*x(ief+i20-1)      
               else
                  xi=efcc(1,ief)
                  yi=efcc(2,ief)
                  zi=efcc(3,ief)
                  dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
                  aij=one/dij
                  dum=dum+aij*x(ief+i20-1)      
               endif
            enddo
            write(luout,9992) jef,dum
         enddo
      endif
c
c     ----- charge screening issue ... -----
c
      chgcos=zero
      do ief=1,nefc
         chgief=-one*x(ief+i20-1)
         chgcos=chgcos+chgief
      enddo
      chgscr=chgcos*screen
      if(out) then
         write(luout,*) '-cosmo- charge screening issue ...'
         write(luout,*) ' cosmo  charge un-screened    = ',chgcos
         write(luout,*) '                  screening   = ',screen
         write(luout,*) ' cosmo  charge    screened    = ',chgscr
      endif
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
         write(luout,*) 'from -rtdb- nelec  = ',nelec
         write(luout,*) 'from -rtdb- charge = ',charge
         write(luout,*) 'from -rtdb- charge_nuclear = ',charge_nuclear
      endif

      ne=nelec

      if(all) then
         chggss=(         -charge)*screen
      elseif(elec) then
         chggss=( dble(ne)       )*screen
      elseif(nucl) then
         chggss=(-charge_nuclear )*screen
      endif
      errcos=chgcos-chggss
      errscr=chgscr-chggss
      if(out) then
         if(all) then
            write(luout,*) ' solute total charge          = ',
     $                  charge
         elseif(elec) then
            write(luout,*) ' solute total electron charge = ',
     $                 -dble(ne)
         elseif(nucl) then
            write(luout,*) ' solute total nuclear  charge = ',
     $                  charge_nuclear
         endif
         write(luout,*) ' gauss  charge                = ',chggss
         write(luout,*) ' error  (un-screened - gauss) = ',errcos
         write(luout,*) ' error  (   screened - gauss) = ',errscr
         write(luout,*) ' ratio  (gauss / un-screened) = ',chggss/chgcos
         write(luout,*) ' ratio  (gauss /    screened) = ',chggss/chgscr
      endif
c
c     ----- set screening factor -----
c
      if(all) then
         if(ifscrn.eq.0) then
            chgfac=one
         elseif(ifscrn.eq.1) then
            chgfac=screen
            chgfac=one
         elseif(ifscrn.eq.2) then
            chgfac=chggss/chgcos
            chgfac=one
         else
            write(luout,*) 'illegal -ifscrn- ... stop. = ',ifscrn
            call errquit('hnd_coschg: illegal -ifscrn-...',911,UERR)
         endif
      else
         if(ifscrn.eq.0) then
            chgfac=one
         elseif(ifscrn.eq.1) then
            chgfac=screen
         elseif(ifscrn.eq.2) then
            chgfac=chggss/chgcos
         else
            write(luout,*) 'illegal -ifscrn- ... stop. = ',ifscrn
            call errquit('hnd_coschg: illegal -ifscrn-...',911,UERR)
         endif
      endif

      if(out) then
         write(luout,*) ' ifscrn                       = ',ifscrn
         write(luout,9984) chgfac
      endif
c
c     ----- apply screening factor -----
c
      do ief=1,nefc
         x(ief+i20-1)=(-chgfac)*x(ief+i20-1)
      enddo
      if(nucl) then
         do ief=1,nefc
            x(ief+i21-1)=x(ief+i20-1)
         enddo
      endif
      if(elec) then
         do ief=1,nefc
            x(ief+i22-1)=x(ief+i20-1)
         enddo
      endif
c
c     ----- store effective charges in -efcz- ... -----
c           check convergence ...
c
      if(iefc_done) then
         if(all) then
            chgcvg=zero
            do ief=1,nefc
               oldief=efcz(ief)
               chgief=   x(ief+i20-1)
               delchg=abs(chgief-oldief)
               if(delchg.gt.chgcvg) then
                  chgcvg=delchg
               endif
               efcz(ief)=chgief
            enddo
         else
            do ief=1,nefc
               efcz(ief)=x(ief+i20-1)
            enddo
         endif
      else
         if(all) then
            chgcvg=zero
            do ief=1,nefc
               oldief=   zero
               chgief=   x(ief+i20-1)
               delchg=abs(chgief-oldief)
               if(delchg.gt.chgcvg) then
                  chgcvg=delchg
               endif
               efcz(ief)=chgief
            enddo
            iefc_done=.true.
         else
            do ief=1,nefc
               efcz(ief)=x(ief+i20-1)
            enddo
         endif
      endif

      if(all) then
         if(out) then
            write(luout,9987) chgcvg
         endif
      endif

      if(dbug) then
         write(luout,9998)
         do ief=1,nefc
            write(luout,9997) ief,(efcc(i,ief),i=1,3),efcz(ief)
         enddo
      endif
c
c     ----- calculate energy terms from -cosmo- charges -----
c
      allefc=zero
      atmefc=zero
      efcefc=zero
      do jef=1,nefc
         xj=efcc(1,jef)
         yj=efcc(2,jef)
         zj=efcc(3,jef)
         qj=efcz(  jef)
c
         allefc=allefc+qj*x(jef+i10-1)
c
         do iat=1,nat
            xi=dbl_mb(k_xyzpt  +3*(iat-1))
            yi=dbl_mb(k_xyzpt+1+3*(iat-1))
            zi=dbl_mb(k_xyzpt+2+3*(iat-1))
            qi=dbl_mb(k_zanpt+iat-1)
            dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
            bij=one/dij
            atmefc=atmefc+qi*bij*qj
         enddo
c
         do ief=1,nefc
            xi=efcc(1,ief)
            yi=efcc(2,ief)
            zi=efcc(3,ief)
            qi=efcz(  ief)
            if(ief.eq.jef) then
               aij=adiag/sqrt(efcs(ief))
            else
               dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
               aij=one/dij
            endif
            efcefc=efcefc+qi*aij*qj
         enddo
      enddo
      efcefc= efcefc/(two*chgfac)
      elcefc= allefc-atmefc
      solnrg= allefc+efcefc
      ecos  = atmefc+efcefc
      if(out) then
         write(luout,9991) atmefc
         write(luout,9990) elcefc
         write(luout,9995) efcefc
         write(luout,9988) solnrg
         write(luout,9989) allefc,(-two*efcefc)
         write(luout,9994) ecos  
      endif
c
c     ----- other form of the solvation energy ... -----
c
      allefc=zero
      atmefc=zero
      elcefc=zero
      do jef=1,nefc
         xj=efcc(1,jef)
         yj=efcc(2,jef)
         zj=efcc(3,jef)
         qj=efcz(  jef)
c
         allefc=allefc+qj*x(jef+i10-1)
         atmefc=atmefc+qj*x(jef+i11-1)
         elcefc=elcefc+qj*x(jef+i12-1)
      enddo
      solnrg= pt5* allefc
      ecos  = pt5*(atmefc-elcefc)
      if(out) then
         write(luout,9991) atmefc
         write(luout,9990) elcefc
         write(luout,9989) allefc
         write(luout,9988) solnrg
         write(luout,9994) ecos
      endif
c
c     ----- after 'elec' pass, get the total ... -----
c
      if(nucl) then
         atmefcn=atmefc
         elcefcn=elcefc
         efcefcn=efcefc
         allefcn=allefc
         solnrgn=solnrg
           ecosn=  ecos
      endif
      if(elec) then
         atmefce=atmefc
         elcefce=elcefc
         efcefce=efcefc
         allefce=allefc
         solnrge=solnrg
           ecose=  ecos
         atmefct=atmefcn+atmefce
         elcefct=elcefcn+elcefce
         allefct=allefcn+allefce
         solnrgt=solnrgn+solnrge
           ecost=  ecosn+  ecose
         if(out) then
            write(luout,*) '-atmefc (nucl     )  - = ',atmefcn
            write(luout,*) '-atmefc (     elec)  - = ',atmefce
            write(luout,*) '-atmefc (nucl+elec)  - = ',atmefct
            write(luout,*) '-elcefc (nucl     )  - = ',elcefcn
            write(luout,*) '-elcefc (     elec)  - = ',elcefce
            write(luout,*) '-elcefc (nucl+elec)  - = ',elcefct
            write(luout,*) '-allefc (nucl     )  - = ',allefcn
            write(luout,*) '-allefc (     elec)  - = ',allefce
            write(luout,*) '-allefc (nucl+elec)  - = ',allefct
            write(luout,*) '-solnrg (nucl     )  - = ',solnrgn
            write(luout,*) '-solnrg (     elec)  - = ',solnrge
            write(luout,*) '-solnrg (nucl+elec)  - = ',solnrgt
            write(luout,*) '-  ecos (nucl     )  - = ',  ecosn
            write(luout,*) '-  ecos (     elec)  - = ',  ecose
            write(luout,*) '-  ecos (nucl+elec)  - = ',  ecost
            write(luout,*) '-efcefc (nucl     )  - = ',efcefcn
            write(luout,*) '-efcefc (     elec)  - = ',efcefce
         endif
      endif
c
c     ----- after 'all' pass, reset charges and energy ... -----
c
      if(.not.noall.and.all) then
         atmefc=atmefct
         elcefc=elcefct
         solnrg=solnrgt
           ecos=  ecost
         efcefc=zero
         do ief=1,nefc
            efcz(ief)=x(ief+i21-1)+x(ief+i22-1)
         enddo

         if(dbug) then
            write(luout,9998)
            do ief=1,nefc
               write(luout,9997) ief,(efcc(i,ief),i=1,3),efcz(ief)
            enddo
         endif
      endif
c
c     ----- -nucl- first, -elec- second, then -all- ... -----
c
      if(.not.all) then
         if(.not.elec) then
            elec=.true.
            nucl=.false.
            go to 10
         else
            if(.not.noall) then
               all =.true.
               elec=.false.
               go to 10
            endif
         endif
      endif
c
      if(out) then
         write(luout,9993)
      endif
c
c     ----- save -cosmo- charges and energy to -rtdb- -----
c
      if (.not. rtdb_put(rtdb,'cosmo:energy',mt_dbl,1,ecos))
     &   call errquit('hnd_coschg: rtdb put failed for ecos',911,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcz',mt_dbl,  nefc,efcz))
     $   call errquit('hnd_coschg: rtdb put failed for efcz',912,
     &       rtdb_err)
c
c     ----- for the time being, save in 'geometry' object -----
c
      if(out) then
         write(luout,*) 'in -hnd_coschg ... for -efc- geom = ',geom
      endif
      status=geom_efc_cart_set(cosmo_geom_efc,nefc,efcc,efcz)
      if (.not.status) then
         call errquit('hnd_coschg: geom_efc_cart_set failed !', 0,
     &       geom_err)
      endif
     
c
c     printing cosmo charges for bq module
c     -------------------------------------
      cosmo_file = "cosmo.xyzq"
      call util_file_name_resolve(cosmo_file,.false.)
      if(ga_nodeid().eq.0) then
        if(.not.util_io_unit(80,90,fn))
     +     call errquit('cannot get free unit', 0,
     +       0)

        open(unit=fn,form="formatted",file=cosmo_file)
        write(*,*) "printing cosmo charges for bq module",
     +     cosmo_file

        do ief=1,nefc
           write(fn,*) 
     +      efcz(  ief),
     +      efcc(1,ief),
     +      efcc(2,ief),
     +      efcc(3,ief)
        end do
        close(fn)
      end if
c
c     ----- release memory block -----
c
      if(.not.ma_chop_stack(i_init)) call
     &  errquit('hnd_coschg, ma_pop_stack of init failed',911,MA_ERR)
c
      return
 9999 format(/,10x,15(1h-),
     1       /,10x,'-cosmo- charges',
     2       /,10x,15(1h-))
 9998 format(4x,'iefc',6x,'x',5x,6x,'y',5x,6x,'z',5x,5x,'q',4x,
     1     /,1x,53(1h-))
 9997 format(1x,i7,3f12.6,f10.6)
 9995 format(' -efcefc- energy = ',f20.12)
 9994 format(' -ecos  - energy = ',f20.12)
 9993 format(' ...... end of -coschg- ......')
 9992 format(' cosmo potential at -ief = ',i6,f16.10)
 9991 format(' -atmefc- energy = ',f20.12)
 9990 format(' -elcefc- energy = ',f20.12)
 9989 format(' -allefc- energy = ',f20.12,f20.12)
 9988 format(' -solnrg- energy = ',f20.12)
 9987 format(' -cosmo- charges convergence = ',f10.6)
 9986 format(' -wfntyp- = ',a8)
 9985 format(' -scftyp- = ',a8)
 9984 format(' applied screening factor -chgfac- = ',f10.6)
      end
c
      subroutine hnd_cosaxd(a,x,ax,n,efcc,efcs)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d+00, one=1.0d+00)
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
      dimension efcc(3,n), efcs(n)
      dimension a(n),x(n),ax(n)
c
      b(i,j)=one/sqrt((efcc(1,i)-efcc(1,j))**2+
     1                (efcc(2,i)-efcc(2,j))**2+
     2                (efcc(3,i)-efcc(3,j))**2)
      d(i  )=adiag/sqrt(efcs(i))
c
      do i=1,n
         dum=zero
         do j=1,n
            if(j.eq.i) then
               aij=d(i)
            else
               aij=b(i,j)
            endif
             xj=x(j)
            dum=dum+aij*xj
         enddo
         ax(i)=dum
      enddo
c
      return
      end
c
      subroutine hnd_cosxad(a,x,xa,n,efcc,efcs)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d+00, one=1.0d+00)
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol 
      dimension efcc(3,n), efcs(n)
      dimension a(n),x(n),xa(n)
c
      b(i,j)=one/sqrt((efcc(1,i)-efcc(1,j))**2+
     1                (efcc(2,i)-efcc(2,j))**2+
     2                (efcc(3,i)-efcc(3,j))**2)
      d(i  )=adiag/sqrt(efcs(i))
c
      do i=1,n
         dum=zero
         do j=1,n
            if(j.eq.i) then
               aij=d(i)
            else
               aij=b(i,j)
            endif
             xj=x(j)
            dum=dum+aij*xj
         enddo
         xa(i)=dum
      enddo
c
      return
      end
c
      subroutine hnd_cosequ(b,x,a,n,g,h,xi,xj,efcc,efcs)
      implicit real*8 (a-h,o-z)
c
c     ----- solve a * x = b , using an iterative procedure       -----
c
c     ----- numerical recipes (p.70), cambridge university press -----
c          w.h.press, b.p.flannery, s.a.teukolsky, w.t.vetterling
c
#include "errquit.fh"
c
      logical     dbug
      double precision efcc(3,n),efcs(n)
      dimension   a(n),x(n),b(n)
      dimension   g(n),h(n),xi(n),xj(n)
      data zero   /0.0d+00/
      data eps    /1.0d-07/
c
      dbug=.false.
      if(dbug) then
         write(luout,*) 'in -cosequ-'
         do i=1,n
            write(luout,9999) i,b(i)
         enddo
      endif
c
      eps2=n*eps*eps
      irst=0
   10 irst=irst+1
      call hnd_cosaxd(a,x,xi,n,efcc,efcs)
      rp=zero
      bsq=zero
      do 11 j=1,n
         bsq=bsq+b(j)*b(j)
         xi(j)=xi(j)-b(j)
         rp=rp+xi(j)*xi(j)
   11 continue
      call hnd_cosxad(a,xi,g,n,efcc,efcs)
      do 12 j=1,n
      g(j)=-g(j)
      h(j)= g(j)
   12 continue
      do 19 iter=1,10*n
         call hnd_cosaxd(a,h,xi,n,efcc,efcs)
         anum=zero
         aden=zero
         do 13 j=1,n
            anum=anum+g(j)*h(j)
            aden=aden+xi(j)*xi(j)
   13    continue
         if(aden.eq.zero) then
            write(luout,*) 'very singular matrix'
            call errquit('hnd_cosequ: singular matrix',911,UERR)
         endif
         anum=anum/aden
         do 14 j=1,n
            xi(j)=x(j)
            x(j)=x(j)+anum*h(j)
   14    continue
         call hnd_cosaxd(a,x,xj,n,efcc,efcs)
         rsq=zero
         do 15 j=1,n
            xj(j)=xj(j)-b(j)
            rsq=rsq+xj(j)*xj(j)
   15    continue
         if(rsq.eq.rp.or.rsq.le.bsq*eps2) return
         if(rsq.gt.rp) then
            do 16 j=1,n
               x(j)=xi(j)
   16       continue
            if(irst.ge.3) return
            go to 10
         endif
         rp=rsq
         call hnd_cosxad(a,xj,xi,n,efcc,efcs)
         gg=zero
         dgg=zero
         do 17 j=1,n
            gg=gg+g(j)*g(j)
            dgg=dgg+(xi(j)+g(j))*xi(j)
   17    continue
         if(gg.eq.zero) return
         gam=dgg/gg
         do 18 j=1,n
            g(j)=-xi(j)
            h(j)=g(j)+gam*h(j)
   18    continue
   19 continue
      write(luout,*) 'too many iterations'
      call errquit('hnd_cosequ: too many iters',911,UERR)
      return
 9999 format(i8,f16.8)
      end
c
      subroutine hnd_linequ(a,lda,b,n,ib,t,deta,ierr,nodcmp)
      implicit real*8 (a-h,o-z)
      dimension a(lda,1),b(n),ib(n),t(n)
c
c     ----- solve a * x = b , with x returned in b -----
c
      if(nodcmp.ne.1) go to 20
c
      call hnd_ludcmp(a,lda,n,ib,t,deta,ierr)
      do 10 j=1,n
   10 deta=deta*a(j,j)
c
   20 continue
      call hnd_lubksb(a,lda,n,ib,b)
c
      return
      end
c
      subroutine hnd_lubksb(a,lda,n,ib,b)
      implicit real*8 (a-h,o-z)
c
      dimension a(lda,1),ib(n),b(n)
      data zero /0.0d+00/
c
      ii=0
      do 12 i=1,n
         ll=ib(i)
         sum=b(ll)
         b(ll)=b(i)
         if(ii.ne.0) then
            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
   11       continue
         else if (sum.ne.zero) then
            ii=i
         endif
         b(i)=sum
   12 continue
      do 14 i=n,1,-1
         sum=b(i)
         if(i.lt.n) then
            do 13 j=i+1,n
               sum=sum-a(i,j)*b(j)
   13       continue
         endif
         b(i)=sum/a(i,i)
   14 continue
      return
      end
c
      subroutine hnd_ludcmp(a,lda,n,ib,vv,d,ierr)
      implicit real*8 (a-h,o-z)
c
      dimension a(lda,1),ib(n),vv(n)
      data tiny     /1.0d-20/
      data zero,one /0.0d+00,1.0d+00/
c
      ierr=0
      d=one
      do 12 i=1,n
         aamax=zero
         do 11 j=1,n
            if( abs(a(i,j)).gt.aamax) aamax= abs(a(i,j))
   11    continue
         if(aamax.eq.zero) then
            ierr=1
            return
         endif
         vv(i)=one/aamax
   12 continue
      do 19 j=1,n
         if(j.gt.1) then
            do 14 i=1,j-1
               sum=a(i,j)
               if(i.gt.1) then
                  do 13 k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
   13             continue
                  a(i,j)=sum
               endif
   14       continue
         endif
         aamax=zero
         do 16 i=j,n
            sum=a(i,j)
            if(j.gt.1) then
               do 15 k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
   15          continue
               a(i,j)=sum
            endif
            dum=vv(i)* abs(sum)
            if(dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
   16    continue
         if(j.ne.imax) then
            do 17 k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
   17       continue
            d=-d
            vv(imax)=vv(j)
         endif
         ib(j)=imax
         if(j.ne.n) then
            if(a(j,j).eq.zero) a(j,j)=tiny
            dum=one/a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
   18       continue
         endif
   19 continue
      if(a(n,n).eq.zero) a(n,n)=tiny
      return
      end
