C> \ingroup stpr_api
C> @{
      subroutine stepper_input(rtdb)
C $Id$
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c::passed
      integer rtdb ! run-time-data-base handle
c
c no input required
c stepper input must follow geometry input!!!!!!
c stepper 
c  [min] find a minimum energy point 
c  [ts]  find a transition state
c  [std] do a stepest descent walk to ts
c  [rstd real] maximum step for stepest descent 
c  [trot] allows scailing of step near a transition state
c  [trtr] allows scailing of step during entire walk
c         this is dangerous since you can miss the minimum
c  [mwc] use mass weighted coordinates
c  [seco] analytic second derivatives in the first step only
c  [swsec] switch to second derivatives and stop
c  [swfir] switch to first derivatives and stop
c  [trust real] initial trust radius 
c  [nmode int] mode to follow
c  [track] track specified or default mode
c  [select 
c     int (3*number of atoms)] mode to select
c  [rpath] reaction path
c  [rconvg real] convergence for reaction path
c  [setmass 
c   int,real (number of atoms [maximum]) 
c   endsetmass]  atomic number and atomic mass
c  [convggm real] maximum gradient convergence tolerance
c  [convgg  real] gradient norme convergence tolerance
c  [convge  real] energy convergence tolerance
c  [maxiter int] maximun number of geometry steps
c  [fdat]  generate finite difference hessian at current 
c          geometry using fd basis
c  [fdopt] generate finite difference hessian at current 
c          geometry using fd basis and optimize geometry in 
c          fd basis.  
c  
c::functions
      integer stpr_get_nat
      external stpr_get_nat
c::local
      logical omin, ts, std, trot, trtr, mwc, seco, swfir, swsec
      logical track, rpath, setm
      logical FF, FT
      parameter (FF = .false., FT = .true.)
      logical oselect, osetmass
      integer hsel, ksel, hm, him
*needs a fix via geom mass interface:rak:      integer km, kim
      integer nmode, nat, maxiter
      double precision rconvg, rstd, trustc, trustm, trust
      double precision convggm, convgg, convge
      character*255 stpr_token
c
      if (nodeid() .ne. 0) return
      oselect  = .false.
      osetmass = .false.
c
      call inp_set_field(0)   ! goto start of line
c
      if (.not. inp_a(stpr_token))
     &    call errquit('stepper_input: no input present',911, INPUT_ERR)
      if (.not. inp_compare(FF,'stepper',stpr_token))
     &    call errquit('stepper_input: not stepper input', 911,
     &       INPUT_ERR)
c
c If stepper module input then turn driver off
c
      if (.not.rtdb_put(rtdb, 'opt:driver', mt_log, 1, FF))
     &    call errquit
     &    ('stepper_input: unable to turn opt:driver off',911, RTDB_ERR)
c
c set defaults
c
      call stpr_set_defaults(omin,ts,std,trot,trtr,mwc,seco,swfir,
     &    swsec,rstd,trustc,trustm,trust,nmode,track,rpath,rconvg,
     &    setm,convggm,convgg,convge, maxiter)
c
 1    continue
      if( inp_cur_field() .lt. inp_n_field() ) then
        if (.not. inp_a(stpr_token))
     &      call errquit
     &      ('stepper_input: expected ascii stpr_token failue',911,
     &       INPUT_ERR)
      else
        if(inp_read()) then
          goto 1
        else
          if (inp_eof()) goto 3         ! all input parsed?
        endif
      endif
c
      if (inp_compare(FF,'end',stpr_token)) then
c     input is finished
        goto 3
      else if (inp_compare(FF,'print',stpr_token)) then
         call util_print_input(rtdb, 'stepper')
      else if (inp_compare(FF,'noprint',stpr_token)) then
         call util_print_input(rtdb, 'stepper')
      else if(inp_compare(FF,'min',stpr_token)) then
        write(luout,*)"Stepper: WARNING: Always use MIN with "
     +                //"TASK ... OPTIMIZE"
        omin = FT
      else if(inp_compare(FF,'ts',stpr_token)) then
        write(luout,*)"Stepper: WARNING: Always use TS with "
     +                //"TASK ... SADDLE"
        ts = FT
        omin = FF
      else if(inp_compare(FF,'std',stpr_token)) then
        omin = FF
        std  = FT
      else if(inp_compare(FF,'rstd',stpr_token)) then
        if(.not.inp_f(rstd))
     &      call errquit
     &      ('stepper_input: error parsing value of rstd',911,
     &       INPUT_ERR)
      else if(inp_compare(FF,'trot',stpr_token)) then
        trot = FT
      else if(inp_compare(FF,'trtr',stpr_token)) then
        trtr = FT
      else if(inp_compare(FF,'mwc',stpr_token)) then
        mwc = FT
      else if(inp_compare(FF,'seco',stpr_token)) then
        seco = FT
      else if(inp_compare(FF,'swsec',stpr_token)) then
        swsec = FT
      else if(inp_compare(FF,'swfir',stpr_token)) then
        swfir = FT
      else if(inp_compare(FF,'trust',stpr_token)) then
        if(.not.inp_f(trust))
     &      call errquit
     &      ('stepper_input: error parsing trust value',911, INPUT_ERR)
      else if(inp_compare(FF,'track',stpr_token)) then
        track = FT
      else if(inp_compare(FF,'nmode',stpr_token)) then
        if(.not.inp_i(nmode))
     &      call errquit
     &      ('stepper_input: error parsing nmode value',911, INPUT_ERR)
      else if(inp_compare(FF,'maxiter',stpr_token)) then
        if(.not.inp_i(maxiter))
     &      call errquit
     &      ('stepper_input: error parsing maxiter value',911,
     &       INPUT_ERR)
      else if(inp_compare(FF,'fdat',stpr_token)) then
c
c.. store fd at flag as true on runtime data base.
c.. absence on rtdb is false.
        if (.not.rtdb_put(rtdb,'fdat',MT_LOG,1,FT))
     &      call errquit
     &      ('stepper_input: failed to store fdat logical',911,
     &       RTDB_ERR)
c
      else if(inp_compare(FF,'fdopt',stpr_token)) then
c
c.. store fd opt flag as true on runtime data base.
c.. absence on rtdb is false.
        if (.not.rtdb_put(rtdb,'fdopt',MT_LOG,1,FT))
     &      call errquit
     &      ('stepper_input: failed to store fdopt logical',911,
     &       RTDB_ERR)
c
      else if(inp_compare(FF,'select',stpr_token)) then
        oselect = FT
        nat = stpr_get_nat(rtdb)
        if (.not.MA_Push_Get(MT_INT,3*nat,'select scratch',hsel,ksel))
     &      call errquit
     &      ('stepper_input: MA failed for select scratch',911, MA_ERR)
        call stpr_parse_select(nat,int_mb(ksel))
*
* write to rtdb and free scratch
*
        call stpr_select_srtdb(nat,rtdb,int_mb(ksel))
        if (.not.MA_Pop_Stack(hsel))
     &      call errquit
     &      ('stepper_input: error freeing select scratch',911, MA_ERR)
c
      else if(inp_compare(FF,'rpath',stpr_token)) then
        rpath = FT
      else if(inp_compare(FF,'rconvg',stpr_token)) then
        if(.not.inp_f(rconvg))
     &      call errquit
     &      ('stepper_input: error parsing rconvg value',911, INPUT_ERR)
      else if(inp_compare(FF,'setmass',stpr_token)) then
        write(luout,*)' SETMASS is no longer supported ' 
        write(luout,*)' use the mass directive in the GEOMETRY object'
        call errquit('stepper_input: obsolete directive ',911,
     &       INPUT_ERR)
*needs a fix via geom mass interface:rak:        osetmass = FT
*needs a fix via geom mass interface:rak:        if (.not.MA_Push_Get
*needs a fix via geom mass interface:rak:     &      (MT_INT,nat,'setmass integer scratch',him,kim))
*needs a fix via geom mass interface:rak:     &      call errquit
*needs a fix via geom mass interface:rak:     &      ('stepper_input: ma error on setmass int',911)
*needs a fix via geom mass interface:rak:        if (.not.MA_Push_Get
*needs a fix via geom mass interface:rak:     &      (MT_DBL,nat,'setmass double  scratch',hm,km))
*needs a fix via geom mass interface:rak:     &      call errquit
*needs a fix via geom mass interface:rak:     &      ('stepper_input: ma error on setmass dble',911)
*needs a fix via geom mass interface:rak:        call stpr_parse_setmass(nat,int_mb(kim),dbl_mb(km))
*needs a fix via geom mass interface:rak:*
*needs a fix via geom mass interface:rak:* write to rtdb and free scratch
*needs a fix via geom mass interface:rak:*
*needs a fix via geom mass interface:rak:        call stpr_setmass_srtdb(nat,rtdb,int_mb(kim),dbl_mb(km))
        if(.not.MA_Pop_Stack(hm))
     &      call errquit
     &      ('stepper_input: error freeing double setmass',911, MA_ERR)
        if(.not.MA_Pop_Stack(him))
     &      call errquit
     &      ('stepper_input: error freeing int setmass',911, MA_ERR)
      else if(inp_compare(FF,'convggm',stpr_token)) then
        if(.not.inp_f(convggm))
     &      call errquit
     &      ('stepper_input: error parsing value of convggm',911,
     &       INPUT_ERR)
      else if(inp_compare(FF,'convgg',stpr_token)) then
        if(.not.inp_f(convgg))
     &      call errquit
     &      ('stepper_input: error parsing value of convgg',911,
     &       INPUT_ERR)
      else if(inp_compare(FF,'convge',stpr_token)) then
        if(.not.inp_f(convge))
     &      call errquit
     &      ('stepper_input: error parsing value of convge',911,
     &       INPUT_ERR)
      else
        write(6,*)'unknown stpr_token: ',stpr_token
        call errquit('stepper_input: fatal error',911, INPUT_ERR)
      endif
      goto 1
c
 3    continue
c
      call stpr_store_rtdb(rtdb,omin,ts,std,trot,trtr,mwc,seco,swfir,
     &    swsec,track,rpath,setm,oselect,osetmass,nmode,rconvg, rstd,
     &    trustc, trustm, trust, convggm, convgg, convge,maxiter)
c
      return
      end
C> @}
*needs a fix via geom mass interface:rak:      subroutine stpr_parse_setmass(nat,imass,mass)
*needs a fix via geom mass interface:rak:      implicit none
*needs a fix via geom mass interface:rak:#include "mafdecls.fh"
*needs a fix via geom mass interface:rak:#include "inp.fh"
*needs a fix via geom mass interface:rak:c
*needs a fix via geom mass interface:rak:c routine to read stepper mass information from input
*needs a fix via geom mass interface:rak:c      
*needs a fix via geom mass interface:rak:c
*needs a fix via geom mass interface:rak:c::passed      
*needs a fix via geom mass interface:rak:      integer nat, imass(nat)
*needs a fix via geom mass interface:rak:      double precision mass(nat)
*needs a fix via geom mass interface:rak:c::local      
*needs a fix via geom mass interface:rak:      integer icount, where
*needs a fix via geom mass interface:rak:      integer idum
*needs a fix via geom mass interface:rak:      double precision xdum
*needs a fix via geom mass interface:rak:      logical result
*needs a fix via geom mass interface:rak:      character*255 stpr_token
*needs a fix via geom mass interface:rak:c
*needs a fix via geom mass interface:rak:      icount = 0
*needs a fix via geom mass interface:rak:      result = .true.
*needs a fix via geom mass interface:rak: 1    continue
*needs a fix via geom mass interface:rak:      if (result) then
*needs a fix via geom mass interface:rak:        if(.not.(inp_cur_field() .lt. inp_n_field())) then
*needs a fix via geom mass interface:rak:          if(.not.inp_read())
*needs a fix via geom mass interface:rak:     &        call errquit
*needs a fix via geom mass interface:rak:     &        ('stpr_parse_setmass: error reading next line',911)          
*needs a fix via geom mass interface:rak:        endif
*needs a fix via geom mass interface:rak:        where = inp_cur_field()
*needs a fix via geom mass interface:rak:        result = inp_i(idum)
*needs a fix via geom mass interface:rak:        result = inp_i(imass(icount))
*needs a fix via geom mass interface:rak:
*needs a fix via geom mass interface:rak:        result = result .and. inp_f(xdum)
*needs a fix via geom mass interface:rak:        if (.not.result) then
*needs a fix via geom mass interface:rak:          call inp_set_field(where)
*needs a fix via geom mass interface:rak:          result = inp_a(stpr_token)
*needs a fix via geom mass interface:rak:          if(inp_compare(.false.,'endsetmass',stpr_token))then
*needs a fix via geom mass interface:rak:            goto 2
*needs a fix via geom mass interface:rak:          else
*needs a fix via geom mass interface:rak:            call errquit
*needs a fix via geom mass interface:rak:     &          ('stpr_parse_setmass: mass input pair not contiguous',
*needs a fix via geom mass interface:rak:     &          911)
*needs a fix via geom mass interface:rak:          endif
*needs a fix via geom mass interface:rak:        endif
*needs a fix via geom mass interface:rak:        icount = icount + 1
*needs a fix via geom mass interface:rak:        imass(icount) = idum
*needs a fix via geom mass interface:rak:        mass(icount)  = xdum
*needs a fix via geom mass interface:rak:        goto 1
*needs a fix via geom mass interface:rak:      else
*needs a fix via geom mass interface:rak:        call errquit('stpr_parse_setmass: error parsing masses ',911)
*needs a fix via geom mass interface:rak:      endif
*needs a fix via geom mass interface:rak:c
*needs a fix via geom mass interface:rak: 2    continue
*needs a fix via geom mass interface:rak:c
*needs a fix via geom mass interface:rak:      end
C> \ingroup stpr_priv
C> @{
      integer function stpr_get_nat(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer nat, rtdb, geom
c
c get number of atoms
c
      if(.not. geom_create(geom, 'geometry'))
     &    call errquit('stpr_get_nat: geom_create?',911, GEOM_ERR)
      if(.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &    call errquit
     &    ('stpr_get_nat: geom input must be before stepper input !!',
     &    911, RTDB_ERR)
      if(.not. geom_ncent(geom,nat))
     &    call errquit('stpr_get_nat: geom_ncent failed?',911, GEOM_ERR)
      if(.not.geom_destroy(geom))
     &    call errquit('stpr_get_nat: geom_destroy failed?',911,
     &       GEOM_ERR)
      stpr_get_nat = nat
      end
      subroutine stpr_parse_select(nat,isel)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
c::passed
      integer nat, isel(3*nat)
c::local
      integer i
c
c routine to read in the select mode.  always (3*nat) in length
c
      do i=1,(3*nat)
        if (.not.(inp_cur_field().lt.inp_n_field())) then
          if(.not.inp_read())
     &        call errquit
     &        ('stpr_parse_select: error reading next line',911,
     &       INPUT_ERR)
        endif
c
        if (.not.inp_i(isel(i)))
     &      call errquit
     &      ('stpr_parse_select: error reading NAT values',911,
     &       INPUT_ERR)
      enddo
c
      end
C> @}
