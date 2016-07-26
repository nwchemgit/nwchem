      subroutine dimqm_input(rtdb)
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
c     Input module for DIM/QM control parameters
c     ------------------------------------------
c          
c     Called from: input_parse.F
c
c     Calls: None
c
c     This subroutine reads the user defined parameters from the
c     dimqm block of the input file.  Defaults are set for
c     optional commands.  All global variables are defined and
c     stored in common blocks located within dimqm.fh
c
c     ----- Input -----
      integer rtdb ! Run time database
c
c     ----- Global Variables -----
c
c      logical ldimqm
c      logical dimqm_go
c      logical lfreq
c      logical lclfld 
c      logical lrsp
c      logical lefield
c      logical ldebug
c      logical latmdip
c      integer alg
c      double precision dimtol
c      double precision scrnFactor
c      double precision finiteField     
c      logical lfirst
c      logical lcmplx
c
c     ----- Local Variables ----
c
      integer nfield
      character*255 field
      character*4 scrnTmp
c
c     ----- Defaults -----
c
      ldimqm            = .true.
      dimqm_go          = .false.
      lfreq             = .false. ! Frequency-dependent DIM
      lclfld            = .false. ! Local Fields
      lrsp              = .true.  ! DIM/QM in linear response
      lefield           = .false. ! Static electric field to be included
      ldebug            = .false. ! Print off DIM/QM debug info
      latmdip           = .true.  ! Print atomic dipoles and charges
      alg               = -1      ! Algorithm for solving DIM equations
      scrnType          = NOSCRN  ! DIM/QM screening function to use
      dimtol            = 1.0d-6  ! Iterative solver cutoff
      scrnFactor        = 1.0d0   ! DIM/QM screening factor
      lfirst            = .true.  ! This is the first SCF cycle
      lcmplx            = .false. ! The FD matrix will be complex
      lagrange          = 0.0d0   ! Lagrange multiplier (CPIM only)
      finiteFieldX      = ZERO
      finiteFieldY      = ZERO
      finiteFieldZ      = ZERO
      dimqm_noseed      = .false. ! Seeding is on by default
c
c     ----- Read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('dimqm_input: unexpected eof',911, INPUT_ERR)
      nfield = inp_n_field()
 20   if (.not. inp_a(field))
     $     call errquit('dimqm_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         go to 40
c
c   --- Methods
c
      else if(inp_compare(.false.,'off', field)) then
         ldimqm = .false.
      else if(inp_compare(.false.,'frequency',field)) then
         lfreq = .true.
         lcmplx = .true.
      else if(inp_compare(.false.,'algorithm',field)) then
         if(.not. inp_i(alg))
     $      call errquit('dimqm_input: alg',911, INPUT_ERR)
      else if(inp_compare(.false.,'tolerance',field)) then
         if(.not. inp_f(dimtol))
     $      call errquit('dimqm_input: tol',911, INPUT_ERR)
c
c   --- Interaction Options
c
      else if(inp_compare(.false.,'efield',field)) then
         lefield = .true.
         if(.not. inp_f(finiteFieldX))
     $     call errquit('dimqm_input: efieldX val', 1, INPUT_ERR)
         if(.not. inp_f(finiteFieldY))
     $     call errquit('dimqm_input: efieldY val', 1, INPUT_ERR)
         if(.not. inp_f(finiteFieldZ))
     $     call errquit('dimqm_input: efieldZ val', 1, INPUT_ERR)

      else if(inp_compare(.false.,'localfield',field)) then
         lclfld = .true.
      else if(inp_compare(.false.,'noresp',field)) then
         lrsp = .false.
      else if(inp_compare(.false.,'screen',field)) then
         if(.not. inp_a_trunc(scrnTmp))
     $      call errquit('dimqm_input: screen type',911, INPUT_ERR)
         if(scrnTmp .eq. "NONE")  scrnType = NOSCRN
         if(scrnTmp .eq. "EXP ")  then
             scrnType = EXPSCRN
             if(.not. inp_f(scrnFactor))
     $         call errquit('dimqm_input: screen factor',911, INPUT_ERR)
         end if
         if( scrnTmp .eq. "ERF " )  then
           scrnType = ERFSCRN
           if( .not. inp_f(scrnFactor) )
     $       call errquit('dimqm_input: screen factor',911, INPUT_ERR)
         end if
c      Seeding
      else if(inp_compare(.false.,'noseed',field)) then
         dimqm_noseed = .true.
c
c   --- Printing Options
c
      else if(inp_compare(.false.,'debug',field)) then
         ldebug = .true.
      else if(inp_compare(.false.,'printatomicdipoles',field)) then
         latmdip = .true.
      end if
      go to 10
c
 40   continue
c
c   --- Valid algorithm choices are -1 through 3
c
      if(alg < -1 .or. alg > 3) then
         write(luout,*) 'Invalid choise for ALGORITHM : ', alg
         call errquit('dimqm_input: alg',911, INPUT_ERR)
      end if
c     Some DIM/QM global variables are needed in the rtdb to
c     reduce module dependency contamination
      write(6,*) ga_nodeid(),' ldimqm value before rt_put ',ldimqm
      if (.not. rtdb_put(rtdb, 'dimqm:ldimqm', mt_log, 1, ldimqm))
     $  call errquit('dimqm_input: ldimqm put failed', 1, RTDB_ERR)         
      if (.not. rtdb_put(rtdb, 'dimqm:lrsp', mt_log, 1, lrsp))
     $  call errquit('dimqm_input: lrsp put failed', 1, RTDB_ERR)         
      if (.not. rtdb_put(rtdb, 'dimqm:lfirst', mt_log, 1, lfirst))
     $  call errquit('dimqm_input: lfirst put failed', 1, RTDB_ERR)         
      if (.not. rtdb_put(rtdb, 'dimqm:alg', mt_int, 1, alg))
     $  call errquit('dimqm_input: alg put failed', 1, RTDB_ERR)         
      if (.not. rtdb_put(rtdb, 'dimqm:dimtol', mt_dbl, 1, dimtol))
     $  call errquit('dimqm_input: dimtol put failed', 1, RTDB_ERR)         
      return
      end subroutine dimqm_input
