C>
C> \ingroup cosmo
C> @{
C>
C> \file cosmo_initialize.F
C> The COSMO initialization
C>
C> \brief Initialize COSMO data structures
C>
C> Initialize the COSMO data structures, in particular this involves
C> computing the positions of the COSMO charges.
C>
      subroutine cosmo_initialize(rtdb,geom,basis,oprint)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "prop.fh"
#include "geom.fh"
#include "bq.fh"
#include "bas.fh"
#include "msgids.fh"
#include "cosmoP.fh"
#include "cosmo_params.fh"
#include "util_params.fh"
c
      integer rtdb   !< [Input] the RTDB handle
      integer geom   !< [Input] the geometry handle
      integer basis  !< [Input] the basis set handle
      logical oprint !< [Input] a print flag
c
      logical status
      logical cosmo
c
      integer iat, nrad, cosmo_screen
      integer l_rad, k_rad, l_costag, k_costag, l_coszan, k_coszan
      integer l_coscoor, k_coscoor
      integer invnuc, nat, nefc, nmap
c
      character*16 ctag(203)
      double precision crad(203)
c
      logical odbug, osome
      integer me
c
c
c MN solvation models -->
c
      logical do_cosmo_smd
      integer icds
      integer do_cosmo_vem
      double precision sola,solb,solc,solg,solh,soln
      double precision pol_cosmo_vem(2), gcds, areacds, disp_cosmo_vem
c
      character*8 solvname_short ! short name of solvent
      character*35 solvname_long ! long name of solvent
c
      character*32 theory
c
c <-- MN solvation models
c
      integer i,ixmem
      double precision arad(203),rad
      logical ofile
      character*255 parfile
      character*16 aname,symb,elem
      integer fn,gn
      character*255 map
      character*255 token
      character*(rtdb_max_file) geom_geom  ! The "geometry" value
      character*(rtdb_max_file) cosmo_geom ! The geometry for which the COSMO
                                           ! parameters were set
      logical util_io_unit,solv_data
      external util_io_unit,solv_data
c
      double precision zero, one, pt5, radtol, def
      integer          mxpass
      data zero   /0.0d+00/
      data one    /1.0d+00/
      data pt5    /0.5d+00/
      data radtol /0.1d+00/
      data mxpass /6/
c
      me=ga_nodeid()
c
      odbug=.false.
      osome=oprint
      call util_print_push
      call util_print_rtdb_load(rtdb, 'cosmo')
c
      if(odbug.and.me.eq.0) then
         write(Luout,*) 'in cosmo_initialize ...'
      endif
c
c     start processing cosmo map
c     --------------------------
      nmap = 0
      if(rtdb_cget(rtdb,'cosmo:parfile',1,token)) then
         map = token
      else if(rtdb_cget(rtdb,'cosmo:map',1,token)) then
         map = token
      else
        map = "none"
      end if
      
      if(map.ne."none") then
  
        if(me.eq.0) then
          call util_getenv('NWCHEM_COSMO_LIBRARY',parfile)
          if(parfile .ne. "") then
            parfile = TRIM(parfile)//"/"//map
            inquire(file=parfile,exist=ofile)
            write(*,*) "looking for parfile in ",TRIM(parfile)
          else
            ofile = .false.
          end if
          if(.not. ofile) then
            parfile = map
            call util_file_name_resolve(parfile,.false.)
            inquire(file=parfile,exist=ofile)
            write(*,*) "looking for parfile in ",TRIM(parfile)
          end if
          if(.not. ofile) then
            parfile = "./"//map
            inquire(file=parfile,exist=ofile)
            write(*,*) "looking for parfile in ",TRIM(parfile)
          end if
          
          if(ofile) then
             write(*,*) "found parfile in ",TRIM(parfile)
          else
             call errquit('cosmo_initialize:cannot find parameter file',
     +       0,0)
             write(*,*) "cannot find parfile"
          end if
  
          call util_flush(6)

          if(.not.util_io_unit(80,90,fn))
     +       call errquit('cosmo_initialize:cannot get free unit', 0,
     +         0)

          open(unit=fn,status="old",form="formatted",file=parfile)
          write(luout,9982) parfile(1:inp_strlen(parfile))
10        continue
          read(fn,*,end=11) aname                 
          if(aname.ne." ") then
           nmap = nmap + 1
           goto 10
          end if
11        continue
          rewind(fn)
          do i=1,nmap
           read(fn,*) ctag(i),crad(i)
           write(luout,9981) ctag(i),crad(i)
          end do
          close(fn)
        end if
        call ga_brdcst(msg_cosmo0,nmap,
     >                 ma_sizeof(mt_int,1,mt_byte),0)
        call ga_brdcst(msg_cosmo1,crad,
     >                 nmap*ma_sizeof(mt_dbl,1,mt_byte),0)
        call ga_brdcst(msg_cosmo2,ctag,
     >                 nmap*16*ma_sizeof(mt_byte,1,mt_byte),0)
        call ga_sync()

        if(.not.rtdb_put(rtdb,'cosmo:crad',mt_dbl,  nmap,crad))
     $     call errquit('cosmo_initialize: rtdb put crad  ',914,
     &         rtdb_err)
c
        if(.not.rtdb_cput(rtdb,'cosmo:ctag',nmap,ctag(1:nmap)))
     $     call errquit('cosmo_initialize: rtdb put ctag',917,
     &         rtdb_err)

      end if
c     end processing cosmo map
c     --------------------------

c
c     Create a Bq instance for the COSMO efc charges, to be used later
c     on
c
      if (.not. bq_create('cosmo efc bq', cosmo_bq_efc))
     $     call errquit('cosmo_initialize: bq_create failed', 0,
     &       GEOM_ERR)
c
c     The next Bq instance is abused to hold a zero array, this
c     array is needed in int_1eefc to keep the ECP contributions
c     coming out right.
c
      if (.not. bq_create('cosmo invnuc bq', cosmo_bq_invnuc))
     $     call errquit('cosmo_initialize: bq_create failed', 0,
     &       GEOM_ERR)
c
c     Try and load the COSMO charges. 
c     1. If loading fails then we hope that this data structure will
c        be properly initialized later on.
c     2. If loading succeeds then we avoid problems if the SCF 
c        wavefunction is converged already. In that case cosmo_charges
c        is not called and the 'cosmo efc geom' geometry is not 
c        initialized. This leads to a catastrophic overwrite of the
c        COSMO geometry. So having it loaded is the best next thing.
c
      if (.not. bq_rtdb_load(rtdb,cosmo_bq_efc))
     $then
c       Hope for the best
      endif
c
c     Get data from the molecule itself
c
      if(.not.geom_ncent(geom,nat)) call errquit
     $    (' cosmo_initialize: geom_ncent    failed.',911, GEOM_ERR)
c
      if(.not.ma_push_get(mt_dbl,nat,'cosmo rads',l_rad,k_rad))
     & call errquit('cosmo_init malloc k_rad failed',911,MA_ERR)
c
      call cosmo_def_radii(rtdb,geom,nat,dbl_mb(k_rad),.true.)
c
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
C      atomloop: do iat=1,nat
C         call util_set_ma_char_arrays(16,
C     $                                byte_mb(k_costag+16*(iat-1)),
C     $                                aname)
C         do i=1,nmap
C            if(inp_compare(.false.,
C     &                     aname,
C     &                     ctag(i))) then
C              dbl_mb(k_rad+iat-1) = crad(i)
C              if(me.eq.0) 
C     &            write(*,9977) ctag(i),crad(i)
C              cycle atomloop
C             
C            endif
C         enddo
C      enddo atomloop
c
      if(.not.ma_pop_stack(l_costag)) call
     &  errquit('cosmo_init,ma_pop_stack of l_tag failed',911, MA_ERR)
c
      do_cosmo_model = DO_COSMO_YK
      dielec   =78.4d+00 ! default water dielectric
      cosmo_screen = COSMO_SCREEN_ST
      cosmo_sccor  = COSMO_SCCOR_LAGRA ! surface charge correction adapted for geometry optimzations
      minbem   =2
      maxbem   =4 ! used in Klamt-Schuurmann model but not in 
                  ! York-Karplus model
      ificos   =0
      lineq    =-1 ! 0: direct solver, 1: iterative solver
      rsolv    =0.50d+00
      zeta     =4.50d+00 ! apparently for Lebedev grids zeta=4.9
      gammas   =1.00d+00 ! full switching, see York, Karplus Eq.86.
      swtol    =1.00d-04 ! switching tolerance, see York, Karplus
                         ! Eq.62-88.
      iangleb  =-1
      fibonacci=-1
      thomson  =0
c
c MN solvation models -->
c
c     reset some defaults if do_cosmo_smd is true
c
      do_cosmo_smd = .false.
      sola = -1.d0
      solb = -1.d0
      solc = -1.d0
      solg = -1.d0
      solh = -1.d0
      soln = -1.d0
      if (dielec.eq.78.4d+00) then
         solvname_short  = 'h2o'
         solvname_long  = 'water'
         sola =  0.82d0
         solb =  0.35d0
         solc =  0.d0
         solg =  103.62d0
         solh =  0.d0
         soln =  1.333d0
         dielecinf = soln*soln ! n**2 (refractive index squarred)
      end if
      status = 
     $ rtdb_get(rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd)
c
      if (do_cosmo_smd) then
          minbem   =3
          maxbem   =4 
          ificos   =0
      end if
c
      status = rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo)
      status = rtdb_get(rtdb,'cosmo:rsolv',mt_dbl,1,rsolv)
      status = rtdb_get(rtdb,'cosmo:do_cosmo_model',mt_int,1,
     +                  do_cosmo_model) 
      if (do_cosmo_model.eq.DO_COSMO_KS) then
        cosmo_sccor  = COSMO_SCCOR_SCALE
        cosmo_screen = COSMO_SCREEN_KS
      else if (do_cosmo_model.eq.DO_COSMO_YK) then
        cosmo_sccor  = COSMO_SCCOR_LAGRA
        cosmo_screen = COSMO_SCREEN_ST
      endif
      status = rtdb_get(rtdb,'cosmo:screen',mt_int,1,cosmo_screen) 
      status = rtdb_get(rtdb,'cosmo:sccor',mt_int,1,cosmo_sccor)
      status = rtdb_get(rtdb,'cosmo:minbem',mt_int,1,minbem) 
      status = rtdb_get(rtdb,'cosmo:maxbem',mt_int,1,maxbem)
      status = rtdb_get(rtdb,'cosmo:ificos',mt_int,1,ificos) 
      status = rtdb_get(rtdb,'cosmo:lineq',mt_int,1,lineq) 
      status = rtdb_get(rtdb,'cosmo:zeta',mt_dbl,1,zeta) 
      status = rtdb_get(rtdb,'cosmo:gamma_s',mt_dbl,1,gammas) 
      status = rtdb_get(rtdb,'cosmo:sw_tol',mt_dbl,1,swtol) 
      status = rtdb_cget(rtdb,'cosmo:solvent',1,solvname_short)
      status = rtdb_get(rtdb,'cosmo:iangleb',mt_int,1,iangleb)
      status = rtdb_get(rtdb,'cosmo:fibonacci',mt_int,1,fibonacci)
      status = rtdb_get(rtdb,'cosmo:thomson',mt_int,1,thomson)
c
c      nrad = 0 ! this block was moved below
c      if (.not. rtdb_get(rtdb,'cosmo:nrrad',mt_int,1,nrad)) nrad = 0
c      if (nrad.gt.nat) call 
c     &   errquit('cosmo_init:nr radii gt nr atoms?',911,MA_ERR)
c      if (nrad.gt.0) then
c         status=rtdb_get(rtdb,'cosmo:radius',mt_dbl,nrad,dbl_mb(k_rad))
c      endif
c
c     get sola,solb,solc,solg,solh,soln parameters, only used if cosmo_smd is true
c     get dielecinf, only used if cosmo_vem is true
c     get dielec, always used 
c
      if(odbug.and.me.eq.0) 
     &     write(luout,*) "calling solv_data"
      if(.not.solv_data(solvname_short,solvname_long, !short and long solvent names
     &          dielec,dielecinf, ! slow and fast components
     &          sola,solb,solc,solg,solh,soln,icds)) then
c
c     if set by the user
c
      status = rtdb_get(rtdb,'cosmo:dielec',mt_dbl,1,dielec)
      status = rtdb_get(rtdb,'cosmo:dielecinf',mt_dbl,1,dielecinf)
      status = rtdb_get(rtdb,'cosmo:sola',mt_dbl,1,sola)
      status = rtdb_get(rtdb,'cosmo:solb',mt_dbl,1,solb)
      status = rtdb_get(rtdb,'cosmo:solc',mt_dbl,1,solc)
      status = rtdb_get(rtdb,'cosmo:solg',mt_dbl,1,solg)
      status = rtdb_get(rtdb,'cosmo:solh',mt_dbl,1,solh)
      status = rtdb_get(rtdb,'cosmo:soln',mt_dbl,1,soln)
      status = rtdb_get(rtdb,'cosmo:icds',mt_int,1,icds)
      endif
c
c     set sola,solb,solc,solg,solh,soln,icds parameters (either from
c     solv_data or user-defined)
c
      if (.not. rtdb_put(rtdb,'cosmo:dielec',mt_dbl,1,dielec))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:dielecinf',mt_dbl,1,dielecinf))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:sola',mt_dbl,1,sola))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:solb',mt_dbl,1,solb))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:solc',mt_dbl,1,solc))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:solg',mt_dbl,1,solg))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:solh',mt_dbl,1,solh))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:soln',mt_dbl,1,soln))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'cosmo:icds',mt_int,1,icds))
     $    call errquit('cosmo_initialize: rtdb_put failed',0, RTDB_ERR)
c
      if (me.eq.0) then
c       write(luout,'(A)') " after solv_data"
       write(luout,*) 
        write(luout,'(A)')   " solvent parameters"
        write(luout,'(A,A)')" solvname_short: ",
     $ solvname_short
        write(luout,'(A,A)')" solvname_long:  ",solvname_long
         if (dielec.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for dielec',0,0)
         if (dielecinf.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for dielecinf',0,0)
        write(luout,'(A,F8.4)') " dielec:      ",dielec
        write(luout,'(A,F8.4)') " dielecinf:   ",dielecinf
        if (do_cosmo_smd) then
         if (icds.eq.1) then
          write(luout,'(A)')   " aqueous SMD model solvent descriptors"
          write(luout,'(A,F8.4)') " dielec:      ",dielec
          write(luout,'(A,F8.4)') " sola:        ",sola
         else if (icds.eq.2) then
          write(luout,'(A)')" nonaqueous SMD model solvent descriptors"
          write(luout,'(A,F8.4)') " dielec:      ",dielec
          write(luout,'(A,F8.4)') " sola:        ",sola
          write(luout,'(A,F8.4)') " solb:        ",solb
          write(luout,'(A,F8.4)') " solc:        ",solc
          write(luout,'(A,F8.4)') " solg:        ",solg
          write(luout,'(A,F8.4)') " solh:        ",solh
          write(luout,'(A,F8.4)') " soln:        ",soln
          if (sola.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for sola',0,0)
          if (solb.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for solb',0,0)
          if (solc.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for solc',0,0)
          if (solg.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for solg',0,0)
          if (solh.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for solh',0,0)
          if (soln.lt.zero) call errquit
     $ ('cosmo_initialize: illegal/undefined value for soln',0,0)
         else
          write(luout,'(A,I1)')   " icds:          ",icds
          call errquit
     $ ('cosmo_initialize: illegal value for icds',0,0)
         endif
        endif
      endif
c 
c     set cosmo radii to smd radii if do_cosmo_smd is true; if the radii
c     are provided by the user then we will use those
c
      if (do_cosmo_smd) then
       call cosmo_def_radii(rtdb,geom,nat,dbl_mb(k_rad),.true.)
      endif
c
      if (.not.
     $rtdb_get(rtdb,'cosmo:do_cosmo_vem',mt_int,1,do_cosmo_vem))
     $  call errquit(
     $ 'cosmo_initialize: cannot get do_cosmo_vem from rtdb',
     $ 0,rtdb_err)
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit(
     $ 'cosmo_initialize: cannot get theory from rtdb',
     $ 0,rtdb_err)
c    
      if (do_cosmo_vem.ne.0.or.do_cosmo_smd) then
        if (theory.ne.'dft'.and.theory.ne.'tddft'.and.
     $      theory.ne.'sodft' .and.theory.ne.'pspw')
     $    call errquit('cosmo-smd can be used only with dft and sodft. 
     $ cosmo-vem can only be used with tddft', 0,0)
      endif
c
      if (do_cosmo_vem.ne.0) then
         if (
     $   rtdb_get(rtdb,'cosmo:pol_cosmo_vem',mt_dbl,2,pol_cosmo_vem))
     $   then
          if (me.eq.0) write(luout,9974)
     $     pol_cosmo_vem(1),pol_cosmo_vem(2)
          disp_cosmo_vem = (dielecinf-1.d0)/(dielecinf+2.d0) * 
     $   (pol_cosmo_vem(2)-pol_cosmo_vem(1)) * 1.850d0 * 0.0015936011d0
          if (.not. 
     $ rtdb_put(rtdb,'cosmo:disp_cosmo_vem',mt_dbl,1,disp_cosmo_vem))
     $     call errquit(
     $ 'cosmo_initialize: cannot put disp_cosmo_vem in rtdb',
     $ 0,rtdb_err)
         endif
      endif
c
c <-- MN solvation models
c
      nrad = 0
      if (.not. rtdb_get(rtdb,'cosmo:nrrad',mt_int,1,nrad)) nrad = 0
      if (nrad.gt.nat) call
     &   errquit('cosmo_init:nr radii gt nr atoms?',911,MA_ERR)
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
         if (do_cosmo_model.eq.DO_COSMO_KS) then
           write(Luout,9970)
         else if (do_cosmo_model.eq.DO_COSMO_YK) then
           write(Luout,9971)
         endif
         write(Luout,9997) dielec
      endif
c
c     ----- screening factor -----
c
      if (cosmo_screen.eq.COSMO_SCREEN_ID) then
         screen=one
         if (osome.and.me.eq.0) write(luout,9984) screen
      else if (cosmo_screen.eq.COSMO_SCREEN_KS) then
         screen=(dielec-one)/(dielec+pt5)
         if (osome.and.me.eq.0) write(luout,9996) screen
      else if (cosmo_screen.eq.COSMO_SCREEN_ST) then
         screen=(dielec-one)/(dielec    )
         if (osome.and.me.eq.0) write(luout,9995) screen
      else
         call errquit("cosmo_initialize: invalid screening method",
     &                cosmo_screen,UERR)
      endif
c
c     ----- surface charge correction -----
c
      if(osome.and.me.eq.0) then
        if (cosmo_sccor.eq.COSMO_SCCOR_NO) then
          write(Luout,9972)"no"
        else if (cosmo_sccor.eq.COSMO_SCCOR_SCALE) then
          write(Luout,9972)"scale"
        else if (cosmo_sccor.eq.COSMO_SCCOR_LAGRA) then
          write(Luout,9972)"lagrangian"
        else
          write(Luout,9972)"*** invalid ***"
        endif
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
         maxbem=min(mxpass,minbem+2)
      endif
c
      ificos=abs(ificos)
c
      do iat=1,nat
         dbl_mb(k_rad+iat-1)=abs(dbl_mb(k_rad+iat-1))
         if(dbl_mb(k_rad+iat-1).lt.radtol) then
            dbl_mb(k_rad+iat-1)=zero
         endif
      enddo
c
c     ---- lebedev & fibonacci ----
c
      if ((iangleb.gt.0).and.(fibonacci.gt.0)) then
        call errquit('cosmo_initialize: Lebedev or Fibonacci?',0,
     $                INPUT_ERR)
      endif

      if ((do_cosmo_model.eq.do_cosmo_ks)) then
        if (iangleb.gt.0) then
          iangleb = -1
        elseif (fibonacci.gt.0) then
          fibonacci = -1
        elseif (thomson.gt.0) then
          thomson = 0
        endif
      endif

      if (zeta.eq.0d0) then
        if (iangleb.gt.0) then
          zeta = 4.907d0
        else
          zeta = 4.500d0
        endif
      endif
c
c     ----- create boundary elements of cavity -----
c
      call hnd_cosset(rtdb,nat,dbl_mb(k_coscoor),dbl_mb(k_rad),geom)
c
      if(.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     $   call errquit('cosmo_initialize: rtdb get failed for nefc',911,
     &       rtdb_err)
c
c     Adjust lineq according to nefc
c
      if ((nefc.ge.7500).and.(lineq.eq.-1)) then
        lineq = 1
      elseif (lineq.ne.1) then
        lineq = 0
      endif
c
      if(osome.and.me.eq.0) then
         write(luout,9993) lineq
         if ((iangleb.lt.1).and.(fibonacci.lt.1)) then
           write(luout,9991) minbem
           if (do_cosmo_model.eq.DO_COSMO_KS) then
             write(luout,9992) maxbem
           endif
           if(ificos.eq.0) then
              write(luout,9990)
           else
              write(luout,9989)
           endif
         elseif (fibonacci.lt.1) then
           write(luout,9969)
           write(luout,9968) iangleb
         else
           write(luout,9967)
           write(luout,9966) fibonacci
         endif
         if (thomson.gt.1) then
           write(luout,9965)
         endif
         if (do_cosmo_model.eq.DO_COSMO_KS) then
           write(luout,9994) rsolv
         endif
         if (do_cosmo_model.eq.DO_COSMO_YK) then
           write(luout,9980) zeta
           write(luout,9979) gammas
           write(luout,9978) swtol
         endif
         write(luout,9988)
      endif
      do iat=1,nat
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
c     ----- done -----
c
      if(osome.and.me.eq.0) then
         write(luout,9986)
         write(luout,9985)
      endif
c
c MN solvation models -->
c
c     ----- add in the smd model non-electrostatic corrections -----
c
      if (do_cosmo_smd) then
c
c        calculate the x memory
         call mnsol_xmem(nat,ixmem)
         if (odbug.and.me.eq.0) write(luout,*) "nat: ", nat
         if (odbug.and.me.eq.0) write(luout,*) "ixmem: ", ixmem
         if (odbug.and.me.eq.0) 
     &       write(luout,*) "calling mnsol_interface()"
c
c        call smd model
         call mnsol_interface(rtdb,geom,nat,ixmem,
     &        sola,solb,solc,solg,solh,soln,icds)
          gcds=0.d0
          areacds=0.d0
          if(.not.rtdb_get(rtdb,'cosmo:gcds',mt_dbl,1,gcds))
     $   call errquit('cosmo_initialize: cannot get gcds from rtdb',
     $ 911,rtdb_err)
          if(.not.rtdb_get(rtdb,'cosmo:areacds',mt_dbl,1,areacds))
     $   call errquit(
     $ 'cosmo_initialize: cannot get areacds from rtdb',
     $ 911,rtdb_err)
         if (me.eq.0) then
          write(luout,9976) gcds*cau2kcal
          write(luout,9975) areacds   
         endif 
      end if
c
c <-- MN solvation models
c
      if(me.eq.0) call util_flush(Luout)
c
      if(.not.ma_pop_stack(l_coscoor)) call
     & errquit('cosmo_initialize,ma_pop_stack l_coor failed',
     & 911,MA_ERR)
      if(.not.ma_pop_stack(l_rad)) call
     & errquit('cosmo_initialize,ma_pop_stack l_rad failed',911,MA_ERR)
c
c     Create invnuc array and zero it
c
      if (.not.bq_alloc(cosmo_bq_invnuc,nefc))
     &  call errquit("cosmo_initialize: could not allocate invnuc",
     &               nefc,MA_ERR)
      if (.not.bq_index_charge(cosmo_bq_invnuc,invnuc))
     &  call errquit("cosmo_initialize: could not get index",
     &               cosmo_bq_invnuc,UERR)
      call dfill(nefc,0.0d0,dbl_mb(invnuc),1)
c
c     Wipe the raw COSMO charges as we potentially have changed the
c     number of charges, the atoms they are associated with, etc.
c
      status = rtdb_delete(rtdb,'cosmo:qraw')
c
      return
 9999 format(/,10x,15(1h-),
     1       /,10x,'-cosmo- solvent',
     2       /,10x,15(1h-))
 9998 format(' no -cosmo- solvent.')
 9997 format(' dielectric constant -eps-     = ',f6.2)
 9996 format(' screen = (eps-1)/(eps+1/2)    = ',f9.5)
 9995 format(' screen = (eps-1)/(eps    )    = ',f9.5)
 9994 format(' solvent radius (ang.)         = ',f7.3)
 9993 format(' -lineq- algorithm             = ',i3)
 9992 format(' -bem- high level              = ',i3)
 9991 format(' -bem- low  level              = ',i3)
 9990 format(' -bem- from -octahedral-')
 9989 format(' -bem- from -icosahedral-')
 9988 format(' atomic radii = ',/,' --------------')
 9987 format(i5,2f7.3)
 9986 format(' ...... end of -cosmo- initialization ......')
 9985 format(/)
 9984 format(' screen =  one                 = ',f9.5)
 9983 format(' charge screening approach     = ',a)
 9982 format(' processing cosmo parameter file ',a)
 9981 format(a19,4x,f9.5)
 9980 format(' gaussian surface charge width = ',f9.5)
 9979 format(' degree of switching           = ',f9.5)
 9978 format(' switching function tolerance  = ',f9.5)
 9977 format(' setting custom COSMO radius: ', a6,f7.3)
c
c MN solvation models -->
c
 9976 format(' G(SMD-CDS) energy (kcal/mol)  = ',f9.3)
 9975 format(' SMD-CDS SASA (angstrom**2)    = ',f9.3,/)
 9974 format(' GS polarizability (angstrom**3) = ',f10.4,/,
     $        ' ES polarizability (angstrom**3) = ',f10.4,/)
c
c <-- MN solvation models
c
 9972 format(' surface charge correction     = ',a)
 9971 format(' Cosmo: York-Karplus, doi: 10.1021/jp992097l')
 9970 format(' Cosmo: Klamt-Schuurmann, doi: 10.1039/P29930000799')

 9969 format(' Using -Lebedev- grid')
 9968 format(' -iangleb',I4)
 9967 format(' Using -Fibonacci spiral- grid')
 9966 format(' -Number of points:',I4)
 9965 format(' Minimizing -Thomson- energy')

      end
C>
C> @}
c $Id$
