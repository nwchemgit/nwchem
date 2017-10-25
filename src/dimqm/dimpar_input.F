      subroutine dimpar_input(rtdb)
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input module for DIM parameters
c     -------------------------------
c
c     Called from: input_parse.F
c
c     Calls: CubicSpline.F
c
c     This subroutine parses and stores all the parameteres for
c     each element within the DIM system.  All global variables
c     are defined and stored in common blocks located within
c     dimqm.fh
c
c     ----- Input Variable -----
c
      integer rtdb ! run-time database
c
c     ----- Global Variables
c
c      logical lcpim
c      integer nDIMTypes
c      integer nDIM
c      integer nOrder
c
c     ----- Local Variables -----
c      
      integer nField, i, j
      integer max_atoms
      double precision xyz(3)
      character*2 element
      character aelement(2)
      integer dims(3)
      integer chunk(3)
      logical stat
      character*255 field
      logical xyzfile
      integer iexp
      double precision expdata(3, 10000)
      double precision splined(2, 10000)
      integer nAtoms                  
      integer nLrtz

c
c     ----- Common variables from dimqm.fh -----
c
c
c     ----- Element Parameters -----
c
      integer          n_types           ! Number of Types
      character*2      typename(92)      ! Symbol
      character        ctemp(2)          ! Temp storage for element symbol
      double precision typepol(92)       ! Polarizability
      double precision typecap(92)       ! Capacitance
      double precision typechar(92)      ! Charge
      double precision typeom1(92)       ! 1st resonant frequency
      double precision typeom2(92)       ! 2nd resonant frequency
      double precision typegm1(92)       ! 1st oscillator width
      double precision typegm2(92)       ! 2nd oscillator width
      double precision typesize(92)      ! Size-dependent parameter
      double precision typerad(92)       ! Atomic Radius
      double precision typebound(92)     ! Bound term in dielectric
      double precision typefermi(92)     ! Fermi Velocity
      double precision typedrude(2,92)   ! Drude term
      character*12     units             ! For unit conversions
      double precision typelrtz(3,50,92) ! Lorentzian term(s), up to 50
      logical          hasExpDie(92)     ! Experimental dielectric defined

c     ----- System Parameters -----
c     These are the handles to the global arrays storing the actual
c     values.
c
      integer l_coords, k_coords
      integer l_typeindex, k_typeindex
c
c
      if (ga_nodeid() .ne. 0) return
c
c     ----- Initialize Parameters -----
c
      nDIMTypes    = 0
      nDIM      = 0
      typename  = ''
      typepol   = -ONE
      typecap   = -ONE
      typechar  = -ONE
      typeom1   = -ONE
      typeom2   = -ONE
      typegm1   = -ONE
      typegm2   = -ONE
      typesize  = -ONE
      typerad   = -ONE
      typedrude = -ONE
      typelrtz  = -ONE
      typebound = ONE
      typefermi = -ONE
      nLrtz    = 1
      xyzfile = .false.
      hasExpDie = .false.
c
c     ----- Read values from input -----
c
      write(LuOut,*) "Start Input Read"
 10   if(.not.inp_read())
     $  call errquit('dimqm_input: unexpected eof',911, INPUT_ERR)
      nField = inp_n_field()
 20   if (.not. inp_a(field))
     $  call errquit('dimqm_input: failed to read field',911,INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         go to 50
      end if
c
c     ----- Atom Parameters -----
c
      if (index(field, ' ') .eq. 3) then
        nDIMTypes = nDIMTypes + 1
        typename(nDIMTypes) = field
c
c       --- Collect values from keywords and check for errors
c
c        write(LuOut,*) "Found atom block ", field
 30     if(.not.inp_read())
     $    call errquit('dimpar: unexpected eof',911,INPUT_ERR)
        if (.not. inp_a(field))
     $    call errquit('dimpar: failed to read field',911,INPUT_ERR)
c
c       END          
        if(inp_compare(.false.,'end',field)) then
          go to 10
        else if(inp_compare(.false.,'subend',field)) then
          go to 10
c
c       POL   
        else if(inp_compare(.false.,'pol',field)) then
          if(.not. inp_f(typepol(nDIMTypes)))
     $      call errquit('dimpar:pol',911,INPUT_ERR)
c
c       CAP
        else if(inp_compare(.false.,'cap',field)) then
          if(.not. inp_f(typecap(nDIMTypes)))
     $      call errquit('dimpar:cap',911,INPUT_ERR)
c
c       CHAR
        else if(inp_compare(.false.,'char',field)) then
          if(.not. inp_f(typechar(nDIMTypes)))
     $      call errquit('dimpar:char',911,INPUT_ERR)
c
c       OM1
        else if(inp_compare(.false.,'om1',field)) then
          if(.not. inp_f(typeom1(nDIMTypes)))
     $      call errquit('dimpar:om1',911,INPUT_ERR)
c
c       OM2
        else if(inp_compare(.false.,'om2',field)) then
          if(.not. inp_f(typeom2(nDIMTypes)))
     $      call errquit('dimpar:om2',911,INPUT_ERR)
c
c       GM1
        else if(inp_compare(.false.,'gm1',field)) then
          if(.not. inp_f(typegm1(nDIMTypes)))
     $      call errquit('dimpar:gm1',911,INPUT_ERR)
c
c       GM2
        else if(inp_compare(.false.,'gm2',field)) then
          if(.not. inp_f(typegm2(nDIMTypes)))
     $      call errquit('dimpar:gm2',911,INPUT_ERR)
c
c       SIZE
        else if(inp_compare(.false.,'size',field)) then
          if(.not. inp_f(typesize(nDIMTypes)))
     $      call errquit('dimpar:size',911,INPUT_ERR)
c
c       RAD
        else if(inp_compare(.false.,'rad',field)) then
          if(.not. inp_f(typerad(nDIMTypes)))
     $      call errquit('dimpar:rad',911,INPUT_ERR)
c
c       DRUDE
        else if(inp_compare(.false.,'drude',field)) then
          stat = inp_f(typedrude(1,nDIMTypes))
          stat = stat .and. inp_f(typedrude(2,nDIMTypes))
          if(inp_n_field() > 3) then
            stat = stat .and. inp_a(units)
c       ADD UNIT CONVERSION HERE
          end if
          if(.not. stat) call errquit('dimpar:drude',911,INPUT_ERR)
c
c       LRTZ
        else if(inp_compare(.false.,'lrtz',field)) then
c          pass
c
c       BOUND
        else if(inp_compare(.false.,'bound',field)) then
          if(.not. inp_f(typebound(nDIMTypes)))
     $      call errquit('dimpar:bound', 911, INPUT_ERR)
c
c       FERMI
        else if(inp_compare(.false.,'fermi',field)) then
          stat = inp_f(typefermi(nDIMTypes))
          if(inp_n_field() > 2) then
            stat = stat .and. inp_a(units)
c       ADD UNIT CONVERSION HERE
          endif
          if(.not. stat) call errquit('dimpar:fermi',911,INPUT_ERR)
c       EXP DIELECTRIC
        else if(inp_compare(.false.,'exp',field)) then
          if(.not.inp_a(field))
     $      call errquit('dimpar: exp missing file path', 9, INPUT_ERR)
          write(luout,*) "Field:", field
          call inp_save_state()
          open(unit=54, file=field, status='old')
          call inp_init(54, luout)
          iexp = 0
          do
            if(.not.inp_read()) exit
            iexp = iexp + 1
            if(iexp > 10000)
     $        call errquit('exp die: too many data points',1,INPUT_ERR)
            if(.not.inp_f(expdata(1,iexp)))
     $        call errquit('exp die: missing wavelength',iexp,INPUT_ERR)
            if(.not.inp_f(expdata(2,iexp)))
     $        call errquit('exp die: missing real die',iexp,INPUT_ERR)
            if(.not.inp_f(expdata(3,iexp)))
     $        call errquit('exp die: missing imag die',iexp,INPUT_ERR)
          end do  

          call inp_restore_state

          close(unit=54)

          call spline (iexp, expdata(1,1:iexp),
     $                 expdata(2,1:iexp), expdata(3,1:iexp),
     $                 splined(1,1:iexp), splined(2,1:iexp))
c jbecca TODO: for some reason, no functionals but unspecified LDA will
c        run with exp dielectric if this following iexp print is not
c        here. Tried 5-4, 5-5, and 5-6 toolkits. Maybe this line allows
c        for badly parallelized loop to complete before storing values?
c          write(luout,*) "iexp"
c          write(luout,*) expdata(1,1:iexp)
c          write(luout,*) expdata(2,1:iexp)
c          write(luout,*) expdata(3,1:iexp)
c          write(luout,*) "Splined data"
c          write(luout,*) splined(1,1:iexp)
c          write(luout,*) splined(2,1:iexp)
          hasExpDie(nDIMTypes) = .true.
          if(.not.rtdb_put(rtdb,'dimpar:expdata'//char(nDIMTypes),
     $                     mt_dbl, 3*iexp, expdata(1:3, 1:iexp)))
     $      call errquit('expdata: rtdb put failed',nDIMTypes,RTDB_ERR)
          if(.not.rtdb_put(rtdb,'dimpar:splinedata'//char(nDIMTypes),
     $                     mt_dbl, 2*iexp, splined(1:2, 1:iexp)))
     $      call errquit('splined: rtdb put failed',nDIMTypes,RTDB_ERR)
          if(.not.rtdb_put(rtdb,'dimpar:iexp'//char(nDIMTypes),
     $                     mt_int, 1, iexp))
     $      call errquit('iexp: rtdb put failed',nDIMTypes,RTDB_ERR)
c
c       Quit on unknown keyword
        else
          call errquit('dimpar: bad keyword!',911,INPUT_ERR)
        endif
        go to 30
c    
c     ----- Coordinates Block -----
      else if(inp_compare(.false.,'xyz',field)) then
c        write(LuOut,*) "Found XYZ block"
        if(.not. inp_read())
     $    call errquit('dimpar:xyz: unexpected eof',911,INPUT_ERR)
c
c       If the next line is the # of atoms, then the coordinate block
c       is in the input file.  Otherwise, check to see if an xyz file
c       was given.
        if(.not. inp_i(nAtoms)) then
         if(.not. inp_a(field))
     $      call errquit('dimpar_input failed xyz file', 9, INPUT_ERR)
          write(luout,*) "XYZ path:", field
          xyzfile = .true.
          open(unit=55, file=field, status='old')
          call inp_save_state()
          call inp_init(55, LuOut)
          ! Get first line (which should be # of atoms)
          if(.not. inp_read())
     $      call errquit('xyzfile: unexpected eof',911,INPUT_ERR)
          if(.not. inp_i(nAtoms))
     $      call errquit('xyzfile: missing nAtoms',911,INPUT_ERR)
        end if
c
c       Create MAs beased off of # of atoms
c       --- Coordinates
          if(.not.ma_push_get(mt_dbl, nAtoms*3, 'dim coords',
     $                    l_coords, k_coords))
     $    call errquit('dimpar_input malloc k_coords failed',1,MA_ERR)
c 
c       --- Type Index
          if(.not.ma_push_get(mt_int, nAtoms, 'dim:type index',
     $                    l_typeindex, k_typeindex))
     $      call errquit('type index malloc failed', 1, MA_ERR)
c
c       Read in line by line, stopping at the next end or EOF for XYZ file
 40       if(.not. inp_read()) then
            if(xyzfile) then
              call inp_restore_state
              close(unit=55)
              ! Take care of END to xyz
              if(.not. inp_read())
     $          call errquit('xyzfile: unexpected eof',911,INPUT_ERR)
              go to 10
            else
              call errquit('dimpar:xyz: unexpected eof',911,INPUT_ERR)
            end if
          end if
          if (.not. inp_a(field))
     $      call errquit('dimpar: failed to read field',911,INPUT_ERR)
c
          if(inp_compare(.false.,'end',field)) then
            go to 10
          else if(inp_compare(.false.,'subend',field)) then
            go to 10
          end if
      
c
          nDIM = nDIM + 1
          element = field
          do i=1,3
            if(.not. inp_f(xyz(i))) then
              nDIM = nDIM - 1
              if(xyzfile) goto 40
              call errquit('dimpar:xyz bad coords',911,INPUT_ERR)
            end if
          end do
c
          dbl_mb(k_coords+3*(nDIM-1))   = xyz(1)*ANGSTROM2BOHR
          dbl_mb(k_coords+3*(nDIM-1)+1) = xyz(2)*ANGSTROM2BOHR
          dbl_mb(k_coords+3*(nDIM-1)+2) = xyz(3)*ANGSTROM2BOHR
c
c       Loop over atom types and when a match is found, put all the
c       data for that atom type into the GA for this atom
          do j=1,nDIMTypes
            if(typename(j) .eq. element) then
              int_mb(k_typeindex+nDIM-1) = j
              go to 40
            endif
          end do ! Loop over nDIMtypes
          ! If we get here, this element doesn't have parameters so soft
          ! exit
          call errquit(element // " has not been defined in input!",
     $                 1, INPUT_ERR)
      end if
      go to 10
c
 50   continue
c
c     Check number of atoms specified versus number read
      if(nDIM .ne. nAtoms)
     $   call errquit("atoms read doesn't match atoms specified",
     $                 1, INPUT_ERR)
c
      nOrder = 3 * nDIM
c
c   Convert Units
      typerad = typerad * ANGSTROM2BOHR
      typedrude = typedrude * EV2HART
      do i = 1, nDIMTypes
        if(typepol(i) < ZERO) then
          typepol(i) = (SIX / PI) * typerad(i)**3
        end if
      end do
c
c   Store in RTDB
      stat = .true.
c      write(luout,*) "Storing"
      stat=stat.and.rtdb_put(rtdb,'dimpar:ndim',mt_int,1,
     $                       ndim)
      stat=stat.and.rtdb_put(rtdb,'dimpar:ndimtypes',mt_int,1,
     $                       ndimtypes)
      stat=stat.and.rtdb_put(rtdb,'dimpar:norder',mt_int,1,
     $                       norder)
! Is there a reason the following is rtdb put twice???
      stat=stat.and.rtdb_put(rtdb,'dimpar:name',mt_byte,nDIMTypes*2,
     $                       typename(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:name',mt_byte,nDIMTypes*2,
     $                       typename(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:sPol',mt_dbl,nDIMTypes,
     $                       typepol(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:cap',mt_dbl,nDIMTypes,
     $                       typecap(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:char',mt_dbl,nDIMTypes,
     $                       typechar(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:om1',mt_dbl,nDIMTypes,
     $                       typeom1(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:om2',mt_dbl,nDIMTypes,
     $                       typeom2(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:gm1',mt_dbl,nDIMTypes,
     $                       typegm1(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:gm2',mt_dbl,nDIMTypes,
     $                       typegm2(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:size',mt_dbl,nDIMTypes,
     $                       typesize(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:rad',mt_dbl,nDIMTypes,
     $                       typerad(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:drude',mt_dbl,nDIMTypes*2,
     $                       typeDrude(1:2, 1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:bound',mt_dbl,nDIMTypes,
     $                       typebound(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:fermi',mt_dbl,nDIMTypes,
     $                       typefermi(1:nDIMTypes))
      stat=stat.and.rtdb_put(rtdb,'dimpar:lrtz',mt_dbl,
     $                       nDIMTypes*3*nLrtz,
     $                       typelrtz(1:3, 1:nLrtz, 1:nDIMTypes))
      stat = stat .and. rtdb_put(rtdb,'dimpar:coords',mt_dbl, nAtoms*3,
     $                           dbl_mb(k_coords))
      stat = stat.and.rtdb_put(rtdb, 'dimpar:typeindex',mt_int, nAtoms,
     $                         int_mb(k_typeindex))
      stat = stat.and.rtdb_put(rtdb, 'dimpar:hasexp',mt_log, nDIMTypes,
     $                         hasExpDie(1:nDIMTypes))
      if(.not.stat) call errquit('error in dimpar_input', 1, RTDB_ERR)
c
c   Deallocate
      stat = ma_chop_stack(l_coords)
c
      end subroutine dimpar_input
