C>
C> \defgroup cosmo COnductor-like Screening MOdel (COSMO)
C>
C> \ingroup cosmo
C> @{
C> 
C> \file cosmo_input.F
C> The COSMO input
C>
C> \brief The COSMO dielectric continuum solvation model
C>
C> COSMO is a dielectric continuum solvation model proposed by
C> Klamt et al. [1] and extended by York et al. [2]. 
C>
C> ### References ###
C>
C> [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="https://doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
C> [2] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="https://doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
C> @}
C>
C> \ingroup cosmo
C> @{
C>
C> \brief Read the COSMO input section
C>
C> COSMO is a dielectric continuum solvation model proposed by
C> Klamt et al. [1] and extended by York et al. [2]. This routine reads
C> the input section and stores the data in the RTDB.
C>
C> ### References ###
C>
C> [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="https://doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
C> [2] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="https://doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
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
#include "cosmoP.fh"
c
      integer rtdb !< [Input] the RTDB handle
c
      integer iat, nfield, nrad, irad
      integer iscren, ifscrn, minbem, maxbem, ificos, lineq
      integer cosmo_screen, cosmo_sccor
      integer iangleb
      integer fibonacci
      integer thomson
      integer l_rad, k_rad
      double precision dielec,dielecinf,rsolv,rad,zeta,gammas,swtol
      character*255 field
      logical do_cosmo
      logical do_gasphase
      logical status
      character*30 tag
      character*255 token
      character*(rtdb_max_file) geom_geom  ! The "geometry" value
      character*(rtdb_max_file) cosmo_geom ! The geometry for which the COSMO
                                           ! parameters were set
c
c     COSMO models
c
      integer do_cosmo_model
c
c MN solvation models -->
c
      character*8 solvname_short ! short name of solvent
      character*35 solvname_long ! long name of solvent
c
c     smd model parameters
c
      logical do_cosmo_smd
      double precision sola
      double precision solb
      double precision solc
      double precision solg
      double precision solh
      double precision soln
      double precision pol_cosmo_vem(2)
      double precision polgs_cosmo_vem, poles_cosmo_vem
      integer icds
c
c     vem model parameters
c
      integer do_cosmo_vem
      integer istep_cosmo_vem
      logical do_cosmo_smd_false
c
c <-- MN solvation models
c
c     ----- defaults -----
c     if solvent name is not given then we will use water solvent
c     descriptors by default
c
c     if solvent name is given but it is not on the list of solvents
c     then the program will be terminated unless all necessary solvent
c     descriptors are found in the input file
c
c     if solvent name is valid then database solvent descriptors will be
c     used which can be overridden by user-provided values
c
      do_cosmo = .true.
      do_cosmo_model = DO_COSMO_YK 
      dielec   =78.4d+00 ! default water dielectric
      solvname_short = 'h2o'
      solvname_long = 'water'
      cosmo_screen = COSMO_SCREEN_ST
      cosmo_sccor  = COSMO_SCCOR_LAGRA ! surface charge correction adapted for geometry optimzations
      minbem   =2
      maxbem   =4 ! used with the Klamt-Schuurmann model
      ificos   =0
      lineq    =0  ! 0 fast direct solver, 1 slow iterative solver
      rsolv    =0.50d+00
      zeta     =4.50d+00 ! apparently for Lebedev grids zeta=4.9
      gammas   =1.00d+00 ! full switching, see York, Karplus Eq.86.
      swtol    =1.00d-04 ! switching tolerance, see York, Karplus Eq.62-88.
      iangleb  =-1
      fibonacci=-1
      thomson  =0
c
c MN solvation models -->
c
c     smd model defaults
c
      do_cosmo_smd = .false.
      if (dielec.eq.78.4d+00) then
         solvname_short  = 'h2o'
         solvname_long  = 'water'
         sola =  0.82d0
         solb =  0.35d0
         solc =  0.d0
         solg =  103.62d0
         solh =  0.d0
         soln =  1.333d0
         dielecinf = soln*soln ! n**2 (refractive index squared)
      end if
c
c     vem model defaults
c
      do_cosmo_vem = 0
      istep_cosmo_vem = 0
      do_cosmo_smd_false = .false.
      pol_cosmo_vem(1) = -99999.9d0
      pol_cosmo_vem(2) = -99999.9d0
c
c     set do_cosmo_smd and do_cosmo_vem flags
      if (.not. 
     $ rtdb_put(rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $ call errquit('cosmo_input: cannot put do_cosmo_smd in rtdb',
     $ 0,rtdb_err)
      if (.not.
     $ rtdb_put(rtdb,'cosmo:do_cosmo_vem',mt_int,1,do_cosmo_vem))
     $ call errquit('cosmo_input: cannot put do_cosmo_vem in rtdb',
     $ 0,rtdb_err)
c
c <-- MN solvation models
c
      if(.not.ma_push_get(mt_dbl,nw_max_atom,'cosmo rads',l_rad,k_rad))
     & call errquit('cosmo_input malloc k_rad failed',nw_max_atom,
     &              MA_ERR)
      call dfill(nw_max_atom,0d0,dbl_mb(k_rad),1)
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
     $     call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
c
c MN solvation models -->
c
c   --- dielecinf
c
      else if(inp_compare(.false.,'dielecinf',field)) then
         if(.not. inp_f(dielecinf))
     $        call errquit('cosmo_input: dielecinf',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:dielecinf',mt_dbl,1,dielecinf))
     $     call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
c
c   --- solvent explicitly specified via short name
c
      else if(inp_compare(.false.,'solvent',field)) then
         if(.not. inp_a(solvname_short))
     $     call errquit('cosmo_input: solvent',911,INPUT_ERR)
         if (solvname_short.eq.'water') then
             solvname_short = 'h2o'
             solvname_long = 'water'
         end if
         if (.not.rtdb_cput(rtdb,'cosmo:solvent',1,solvname_short))
     $     call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
c
c   --- do_cosmo_vem
c
      else if(inp_compare(.false.,'do_cosmo_vem',field)) then
         if(.not. inp_i(do_cosmo_vem))
     $        call errquit('cosmo_input: do_cosmo_vem',911, INPUT_ERR)
         if (do_cosmo_vem.lt.0.or.do_cosmo_vem.gt.2) do_cosmo_vem = 0
         if (.not. 
     $      rtdb_put(rtdb,'cosmo:do_cosmo_vem',mt_int,1,do_cosmo_vem))
     $      call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- pol_cosmo_vem
c
      else if(inp_compare(.false.,'polgs_cosmo_vem',field)) then
         if(.not. inp_f(polgs_cosmo_vem))
     $        call errquit(
     $ 'cosmo_input: polgs_cosmo_vem',911, INPUT_ERR)
         pol_cosmo_vem(1) = polgs_cosmo_vem
      else if(inp_compare(.false.,'poles_cosmo_vem',field)) then
         if(.not. inp_f(poles_cosmo_vem))
     $        call errquit(
     $ 'cosmo_input: poles_cosmo_vem',911, INPUT_ERR)
         pol_cosmo_vem(2) = poles_cosmo_vem
         if (pol_cosmo_vem(1).lt.-9999.d0) 
     $        call errquit(
     $ 'cosmo_input: error in polgs_cosmo_vem',911,INPUT_ERR)
         if (pol_cosmo_vem(2).lt.-9999.d0)
     $        call errquit(
     $ 'cosmo_input: error in poles_cosmo_vem',911,INPUT_ERR)
         if (.not.
     $      rtdb_put(rtdb,'cosmo:pol_cosmo_vem',mt_dbl,2,pol_cosmo_vem))
     $      call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- do_cosmo_smd
c
      else if(inp_compare(.false.,'do_cosmo_smd',field)) then
         if (.not. inp_l(do_cosmo_smd)) do_cosmo_smd = .true.
         if (.not. rtdb_put(rtdb,'cosmo:do_cosmo_smd',mt_log,1,
     $    do_cosmo_smd))
     $    call errquit('cosmo_input: rtdb_put failed',0, RTDB_ERR)
         if (.not.do_cosmo_smd) do_cosmo_smd_false = .true.
c
c  --- smd definitions: sola,solb,solc,solg,solh,soln
c
      else if(inp_compare(.false.,'sola',field)) then
        status = inp_f(sola)
        if (.not. rtdb_put(rtdb,'cosmo:sola',mt_dbl,1,sola))
     $    call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
      else if(inp_compare(.false.,'solb',field)) then
        status = inp_f(solb)
        if (.not. rtdb_put(rtdb,'cosmo:solb',mt_dbl,1,solb))
     $    call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
      else if(inp_compare(.false.,'solc',field)) then
        status = inp_f(solc)
        if (.not. rtdb_put(rtdb,'cosmo:solc',mt_dbl,1,solc))
     $    call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
      else if(inp_compare(.false.,'solg',field)) then
        status = inp_f(solg)
        if (.not. rtdb_put(rtdb,'cosmo:solg',mt_dbl,1,solg))
     $    call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
      else if(inp_compare(.false.,'solh',field)) then
        status = inp_f(solh)
        if (.not. rtdb_put(rtdb,'cosmo:solh',mt_dbl,1,solh))
     $    call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
      else if(inp_compare(.false.,'soln',field)) then
        status = inp_f(soln)
        if (.not. rtdb_put(rtdb,'cosmo:soln',mt_dbl,1,soln))
     $    call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
      else if(inp_compare(.false.,'icds',field)) then
        status = inp_i(icds)
        if (.not. rtdb_put(rtdb,'cosmo:icds',mt_int,1,icds))
     $    call errquit('cosmo_input: rtdb put failed',912,RTDB_ERR)
c
c <-- MN solvation models
c
c   --- off  
c
      else if(inp_compare(.false.,'off',field)) then
         do_cosmo=.false.
c
c   --- do_cosmo_ks: switch to Klamt-Schuurmann model
c
      else if(inp_compare(.false.,'do_cosmo_ks',field)) then
         do_cosmo_model=DO_COSMO_KS
         if (.not. rtdb_put(rtdb,'cosmo:do_cosmo_model',mt_int,1,
     $                      do_cosmo_model))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- do_cosmo_yk: switch to York-Karplus model
c
      else if(inp_compare(.false.,'do_cosmo_yk',field)) then
         do_cosmo_model=DO_COSMO_YK
         if (.not. rtdb_put(rtdb,'cosmo:do_cosmo_model',mt_int,1,
     $                      do_cosmo_model))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- rsolv: Used in Klamt-Schuurmann, but not in York-Karplus
c   ---        approach
c
      else if(inp_compare(.false.,'rsolv',field)) then
         if(.not. inp_f(rsolv))
     $        call errquit('cosmo_input: rsolv',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:rsolv',mt_dbl,1,rsolv))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- parameter file
c
      else if(inp_compare(.false.,'parameters',field)) then
         if(.not.inp_a(token)) call errquit('no parameters',0,
     &         INPUT_ERR)
         if (.not.rtdb_cput(rtdb,'cosmo:parfile',1,token))
     >        call errquit('failed to store cosmo:parfile',0,
     >         RTDB_ERR)
c
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
c   --- iscren !!! Deprecated: use SCREEN instead
c
      else if(inp_compare(.false.,'iscren',field)) then
         if(.not. inp_i(iscren))
     $        call errquit('cosmo_input: iscren',911, INPUT_ERR)
         if (iscren.eq.0) then
           cosmo_screen = COSMO_SCREEN_ST
         else if (iscren.gt.0) then 
           cosmo_screen = COSMO_SCREEN_KS
         else
           cosmo_screen = COSMO_SCREEN_ID
         endif
         if (.not. rtdb_put(rtdb,'cosmo:screen',mt_int,1,cosmo_screen))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- SCREEN: Dielectric screening to use
c       SCREEN IDEAL : No screening to be used f(e)=1
c       SCREEN KS    : Klamt-Schuurmann approach to be used
c                      f(e)=(e-1)/(e+0.5)
c       SCREEN ST    : Stefanovich-Truong approach to be used
c                      f(e)=(e-1)/e
c       the default is: 
c       - SCREEN KS : with the Klamt-Schuurmann model
c       - SCREEN ST : with the York-Karplus model
c
      else if(inp_compare(.false.,'screen',field)) then
        if (.not.inp_a(token))
     $        call errquit('cosmo_input: screen',911, INPUT_ERR)
        if (inp_compare(.false.,'ideal',token)) then
          cosmo_screen = COSMO_SCREEN_ID
        else if (inp_compare(.false.,'ks',token)) then
          cosmo_screen = COSMO_SCREEN_KS
        else if (inp_compare(.false.,'st',token)) then
          cosmo_screen = COSMO_SCREEN_ST
        else
          write(luout,*)"invalid option for SCREEN, valid options are:"
          write(luout,*)"- ST    for Stefanovich-Truong f(e)=(e-1)/e"
          write(luout,*)"- KS    for Klamt-Schuurmann   "
     $                //"f(e)=(e-1)/(e+1/2)"
          write(luout,*)"- IDEAL for ideal screening    f(e)=1"
          call util_flush(luout)
          call errquit("cosmo_input: invalid option for SCREEN",911,
     $                 INPUT_ERR)
        endif
        if (.not. rtdb_put(rtdb,'cosmo:screen',mt_int,1,cosmo_screen))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- ifscrn !!! Deprecated use CHARGE_CORRECTION instead
c
      else if(inp_compare(.false.,'ifscrn',field)) then
         if(.not. inp_i(ifscrn))
     $        call errquit('cosmo_input: ifscrn',911, INPUT_ERR)
         if (ifscrn.eq.2) then
           cosmo_sccor = COSMO_SCCOR_LAGRA
         else if (ifscrn.eq.1) then
           cosmo_sccor = COSMO_SCCOR_SCALE
         else if (ifscrn.eq.0) then
           cosmo_sccor = COSMO_SCCOR_NO
         endif
         if (.not. rtdb_put(rtdb,'cosmo:sccor',mt_int,1,cosmo_sccor))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- CHARGE_CORRECTION: Surface charge correction to use
c       CHARGE_CORRECTION NO         : No surface charge corrections to
c                                      be used
c       CHARGE_CORRECTION SCALE      : Scale the electron surface
c                                      charges and the nuclear surface
c                                      charges
c       CHARGE_CORRECTION LAGRANGIAN : Use a Lagrangian constraint to
c                                      variationally correct the charges
c       the default is: CHARGE_CORRECTION LAGRANGIAN unless the
c       Klamt-Schuurmann COSMO model is used when the default is
c       CHARGE_CORRECTION SCALE
c
      else if(inp_compare(.false.,'charge_correction',field)) then
        if (.not.inp_a(token))
     $        call errquit('cosmo_input: charge_correction',911,
     $                     INPUT_ERR)
        if (inp_compare(.false.,'no',token)) then
          cosmo_sccor = COSMO_SCCOR_NO
        else if (inp_compare(.false.,'off',token)) then
          cosmo_sccor = COSMO_SCCOR_NO
        else if (inp_compare(.false.,'scale',token)) then
          cosmo_sccor = COSMO_SCCOR_SCALE
        else if (inp_compare(.false.,'lagrangian',token)) then
          cosmo_sccor = COSMO_SCCOR_LAGRA
        else
          write(luout,*)"invalid option for CHARGE_CORRECTION, valid "
     $                //"options are:"
          write(luout,*)"- LAGRANGIAN for a variational constraint"
          write(luout,*)"- SCALE      for a plain charge scaling"
          write(luout,*)"- NO         for no correction"
          call util_flush(luout)
          call errquit("cosmo_input: invalid option for "
     $               //"CHARGE_CORRECTION,",911,INPUT_ERR)
        endif
        if (.not. rtdb_put(rtdb,'cosmo:sccor',mt_int,1,cosmo_sccor))
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
c   --- maxbem: used in the Klamt-Schuurmann approach but not in 
c   ---         the York-Karplus
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
c   --- lebedev
c
      else if(inp_compare(.false.,'lebedev',field)) then
         if(.not. inp_i(iangleb)) iangleb = 4
         if (iangleb.lt.1)
     $     call errquit('cosmo_input: lebedev <1',911,INPUT_ERR)
         if(.not. rtdb_put(rtdb,'cosmo:iangleb',mt_int,1,iangleb))
     $     call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
c
c   --- fibonacci
c
      else if(inp_compare(.false.,'fibonacci',field)) then
         if(.not. inp_i(fibonacci)) fibonacci = 86
         if (fibonacci.lt.1)
     $     call errquit('cosmo_input: fibonacci <1',911,INPUT_ERR)
         if(.not. rtdb_put(rtdb,'cosmo:fibonacci',mt_int,1,fibonacci))
     $     call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
c
c   --- Thomson
c
      else if(inp_compare(.false.,'thomson',field)) then
        thomson = 1
        if(.not. rtdb_put(rtdb,'cosmo:thomson',mt_int,1,thomson))
     $     call errquit('cosmo_input: rtdb put failed',911,RTDB_ERR)
c
c   --- lineq
c
      else if(inp_compare(.false.,'lineq',field)) then
         if(.not. inp_i(lineq))
     $        call errquit('cosmo_input: lineq',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:lineq',mt_int,1,lineq))
     $        call errquit('cosmo_input: rtdb put failed',911, RTDB_ERR)
c
c   --- do_gasphase
c
      else if(inp_compare(.false.,'do_gasphase',field)) then
         if (.not. inp_l(do_gasphase)) do_gasphase = .true.
         if (do_gasphase) then
           if (.not. rtdb_put(rtdb,'cosmo_phase',mt_int,1,1))
     $        call errquit('cosmo_input: put cosmo_phase failed',
     $                     911, RTDB_ERR)
         else
           if (.not. rtdb_put(rtdb,'cosmo_phase',mt_int,1,2))
     $        call errquit('cosmo_input: put cosmo_phase failed',
     $                     911, RTDB_ERR)
         endif
c
c   --- zeta ! the exponent of the Gaussian charge distributions
c            ! that represent the surface charges
c
      else if(inp_compare(.false.,'zeta',field)) then
         if (.not. inp_f(zeta))
     $      call errquit('cosmo_input: zeta not a floating point '
     $                   //'number',911,INPUT_ERR)
         if (zeta.le.0.0d0)
     $      call errquit('cosmo_input: zeta must be a positive '
     $                   //'number',911,INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:zeta',mt_dbl,1,zeta))
     $      call errquit('cosmo_input: put zeta failed',
     $                   911, RTDB_ERR)
c
c   --- gamma_s ! the degree of switching
c
      else if(inp_compare(.false.,'gamma_s',field)) then
         if (.not.inp_f(gammas))
     $      call errquit('cosmo_input: gamma_s not a floating point '
     $                   //'number',911,INPUT_ERR)
         if (gammas.lt.0.0d0) 
     $      call errquit('cosmo_input: gamma_s must be >= 0',911,
     $                   INPUT_ERR)
         if (gammas.gt.1.0d0) 
     $      call errquit('cosmo_input: gamma_s must be =< 1',911,
     $                   INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:gamma_s',mt_dbl,1,gammas))
     $      call errquit('cosmo_input: put cosmo:gamma_s failed',
     $                   911, RTDB_ERR)
c
c   --- sw_tol ! the switching tolerance
c
      else if(inp_compare(.false.,'sw_tol',field)) then
         if (.not.inp_f(swtol))
     $      call errquit('cosmo_input: sw_tol not a floating point '
     $                   //'number',911,INPUT_ERR)
         if (swtol.lt.0.0d0) 
     $      call errquit('cosmo_input: sw_tol must be >= 0',911,
     $                   INPUT_ERR)
         if (swtol.gt.1.0d0) 
     $      call errquit('cosmo_input: sw_tol must be =< 1',911,
     $                   INPUT_ERR)
         if (.not. rtdb_put(rtdb,'cosmo:sw_tol',mt_dbl,1,swtol))
     $      call errquit('cosmo_input: put cosmo:sw_tol failed',
     $                   911, RTDB_ERR)
c
      else if(inp_compare(.false.,'print',field)) then
         call util_print_input(rtdb, 'cosmo')
      else
         if (ga_nodeid().eq.0) then
           write(LuOut,*)'cosmo_input: read unknown keyword: ',field
           call util_flush(LuOut)
         endif
         call errquit('cosmo_input: unknown keyword',911,INPUT_ERR)
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
c MN solvation models -->
c
c if VEM is called and do_cosmo_smd above is not set explicitly 
c to .false. then set do_cosmo_smd = .true. by default
c
      if (do_cosmo_vem.ne.0.and..not.do_cosmo_smd_false) then
        do_cosmo_smd = .true.
        if (.not. 
     $ rtdb_put(rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $  call errquit('cosmo_input: cannot put do_cosmo_smd in rtdb',
     $  0,rtdb_err)
      endif
c
      if (.not. rtdb_put(rtdb,'cosmo:istep_cosmo_vem',mt_int,1,0))
     $ call errquit('cosmo_input: cannot put istep_cosmo_vem in rtdb',
     $ 0,rtdb_err)
c
c <-- MN solvation models
c

      return
      end
c $Id$
