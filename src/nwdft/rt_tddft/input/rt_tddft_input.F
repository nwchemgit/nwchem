C
C     rt_tddft_input.F
C     
C     Parses input deck for rt-tddft parameters.
C     
C
      subroutine rt_tddft_input (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in) :: rtdb


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_input: "
      integer, parameter      :: maxnum = 999999999 ! dummy # when using "*" as # prints

C     == Varibles ==
      logical done
      character*255 test, curr_popt, loadopts, loadtarget
      logical got_opt

      double precision tmin, tmax, dt
      
      logical prof, noprop, static, nodisk, matrix_checks
C      logical dplot_do, subgs
C      character*20 dplot_opts
      integer nchecks, nprints, nrestarts
C      integer nsnapshots
      character*20 num_str
      integer checklvl
      character*16 field_name
      character*255 tag_in, method
      integer prop_method, exp_method
      logical use_dmat

      
C     (parameters that must be supplied--no defaults)
      logical got_tmin, got_tmax, got_dt
      integer nfields, nexcites
      
      

C     == External ==
      integer, external :: atoi


C     (values which have no defaults and must be supplied)
      got_tmin = .false.
      got_tmax = .false.
      got_dt = .false.

      nfields = 0
      nexcites = 0

      if (.not.rtdb_put(rtdb,'rt_tddft:nfields',
     $     mt_int,1,nfields))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:nexcites',
     $     mt_int,1,nexcites))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      


C
C     Dump defaults into rtdb (will be replaced later if values found in
C     input deck).
C
      call rt_tddft_input_put_defaults (rtdb)


      
C     
C     Parse the input; we will put in rtdb later after checking.
C     
      done = .false.
      do while (.not. done)

         
         if (.not. inp_read()) 
     $        call errquit(pname//'Read failed input',0, INPUT_ERR)
         if (.not. inp_a(test)) 
     $        call errquit(pname//'Read failed keyword',0, INPUT_ERR)



C
C     Tag/title for the run.
C
         if (inp_compare(.false.,test,"tag")) then
            if (.not.inp_a(tag_in)) then
               call errquit (pname//"failed to parse tag",0,0)
            endif
            
            if (len_trim(tag_in) .gt. 24)
     $           call errquit (pname//"tag too long, max length is 24",
     $           0,0)
            
            if (.not. rtdb_cput (rtdb, "rt_tddft:tag",
     $           1, trim(tag_in)))
     $           call errquit (pname//
     $           "failed to put 'tag' target into rtdb",0,0)



C
C     Look for starting state to load.
C
         elseif (inp_compare(.false.,test,"load")) then
            if (.not.inp_a(loadopts)) then
               call errquit (pname//"failed to parse load option",0,0)
            endif
            
            if (loadopts .eq. "vectors") then
               if (.not.inp_a(loadtarget)) then
                  call errquit (pname//
     $                 "failed to parse 'load vectors' target",0,0)
               endif
               if (.not. rtdb_cput (rtdb, "rt_tddft:init_movecs",
     $              1, loadtarget)) call errquit (pname//
     $              "failed to put 'init_movecs' target into rtdb",0,0)

            elseif (loadopts .eq. "density") then
               use_dmat=.true.
               if (.not.rtdb_put(rtdb,'rt_tddft:use_dmat',mt_log,
     $              1,use_dmat))
     $            call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c              Read total (if CS) or alpha (if OS) dmat
               if (.not.inp_a(loadtarget)) then
                  call errquit (pname//
     $                 "failed to parse 'load density' target",0,0)
               endif
               if (.not. rtdb_cput (rtdb, "rt_tddft:init_dmat1",
     $              1, loadtarget)) call errquit (pname//
     $              "failed to put 'init_density' target into rtdb",0,0)
c              Read beta dmat if given (i.e. doing OS)
               if (inp_a(loadtarget)) then
                 if (.not. rtdb_cput (rtdb, "rt_tddft:init_dmat2",
     $               1, loadtarget)) call errquit (pname//
     $              "failed to put 'init_density' target into rtdb",0,0)
               endif

            elseif (loadopts .eq. "scf") then
C     (no need to do anything--will use SCF vectors as starting point if nothing in rtdb)
               
            elseif (loadopts .eq. "restart") then ! will look for *.rt_restart file
               if (.not.rtdb_put (rtdb, "rt_tddft:restart",
     $              mt_log, 1, .true.))
     $              call errquit(pname//
     $              'Write failed to rtdb',0,RTDB_ERR)
               
            else
               call errquit (pname//"invalid 'load' option: "
     $              //loadopts,0,0)
            endif
                  
            
C
C     applied field sub-block
C
         elseif (inp_compare(.false.,test,'field')) then
            
            if (.not. inp_a (field_name))
     $           call errquit (pname//"failed to read field name",0,0)

            nfields = nfields + 1
            
            if (.not.rtdb_put(rtdb,'rt_tddft:nfields',mt_int,1,nfields))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            call rt_tddft_input_field (rtdb, field_name, nfields)


            

C
C     visualization sub-block
C
         elseif (inp_compare(.false.,test,'visualization')) then
            call rt_tddft_input_visualization (rtdb)


C
C     MO CAP sub-block
C
         elseif (inp_compare(.false.,test,'mocap')) then
            call rt_tddft_input_mocap (rtdb)


C 
C     Spatial CAP sub-block
C
         elseif (inp_compare(.false.,test,'spatialcap')) then
            call rt_tddft_input_spatialcap (rtdb)
            
            
C     
C     tmin
C     
C     XXX HARDCODED TO -5*dt
C
c$$$         elseif (inp_compare(.false.,test,'tmin')) then
c$$$            if (.not.inp_f(tmin)) call errquit (pname//
c$$$     $           "tmin takes a float", 0, 0)
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:tmin',mt_dbl,1,tmin))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c$$$            got_tmin = .true.


            
C     
C     tmax
C     
         elseif (inp_compare(.false.,test,'tmax')) then
            if (.not.inp_f(tmax)) call errquit (pname//
     $           "tmax takes a float", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_tddft:tmax',mt_dbl,1,tmax))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            got_tmax = .true.
            
            
C     
C     dt
C     
         elseif (inp_compare(.false.,test,'dt')) then
            if (.not.inp_f(dt)) call errquit (pname//
     $           "dt takes a float >= 0", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_tddft:dt',mt_dbl,1,dt))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            got_dt = .true.


C     
C     checklvl
C     
         elseif (inp_compare(.false.,test,'checklvl')) then
            if (.not.inp_i(checklvl)) call errquit (pname//
     $           "checklvl takes a value of 1, 2, or 3", 0, 0)
            
            if (.not.rtdb_put(rtdb,'rt_tddft:checklvl',
     $           mt_int,1,checklvl))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)





         elseif (inp_compare(.false.,test,'tolerances')) then
            call rt_tddft_input_tolerances (rtdb)
            
C     
C     nchecks
C     
         elseif (inp_compare(.false.,test,'nchecks')) then
            if (.not.inp_a(num_str))
     $           call errquit (pname//
     $           "nchecks takes an int >= 0 (or *)", 0, 0)

            if ( trim(num_str) .eq. "*") then
               nchecks = maxnum
            else
               nchecks = atoi (num_str)
               if (nchecks < 0)
     $              call errquit (pname//
     $              "nchecks takes an int >= 0 (or *)", 0, 0)
            endif

            if (.not.rtdb_put(rtdb,'rt_tddft:nchecks',mt_int,1,nchecks))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)




C     
C     nprints (note * means every step, uses large dummy number)
C     
         elseif (inp_compare(.false.,test,'nprints')) then
            
            if (.not.inp_a(num_str))
     $           call errquit (pname//
     $           "nprints takes an int >= 0 (or *)", 0, 0)

            if ( trim(num_str) .eq. "*") then
               nprints = maxnum
            else
               nprints = atoi (num_str)
               if (nprints < 0)
     $              call errquit (pname//
     $              "nprints takes an int >= 0 (or *)", 0, 0)
            endif

            if (.not.rtdb_put(rtdb,'rt_tddft:nprints',mt_int,1,nprints))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)



C     
C     nsnapshots
C     
c$$$         elseif (inp_compare(.false.,test,'nsnapshots')) then
c$$$            if (.not.inp_a(num_str))
c$$$     $           call errquit (pname//
c$$$     $           "nsnapshots takes an int >= 0 (or *)", 0, 0)
c$$$
c$$$            if ( trim(num_str) .eq. "*") then
c$$$               nsnapshots = maxnum
c$$$            else
c$$$               nsnapshots = atoi (num_str)
c$$$               if (nsnapshots < 0)
c$$$     $              call errquit (pname//
c$$$     $              "nsnapshots takes an int >= 0 (or *)", 0, 0)
c$$$            endif
c$$$            
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:nsnapshots',
c$$$     $           mt_int,1,nsnapshots))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C     
C     nrestarts
C     
         elseif (inp_compare(.false.,test,'nrestarts')) then
            if (.not.inp_a(num_str))
     $           call errquit (pname//
     $           "nrestarts takes an int >= 0 (or *)", 0, 0)

            if ( trim(num_str) .eq. "*") then
               nrestarts = maxnum
            else
               nrestarts = atoi (num_str)
               if (nrestarts < 0)
     $              call errquit (pname//
     $              "nrestarts takes an int >= 0 (or *)", 0, 0)
            endif

            if (.not.rtdb_put(rtdb,'rt_tddft:nrestarts',
     $           mt_int,1,nrestarts))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            


C
C     Propagator
C
         elseif (inp_compare(.false.,test,'propagator')) then

            if (.not.inp_a(method)) then
               call errquit (pname//
     $              "failed to parse propagator method",0,0)
            endif

            if (method.eq."euler") then
               prop_method = 1  ! euler

            elseif (method.eq."rk4") then
               prop_method = 2  ! 4th order runge-kutta
               
            elseif (method.eq."magnus") then
               prop_method = 3  ! 2nd order magnus
               
            else
               call errquit (pname//"invalid propagator: "//method)
            endif

            if (.not.rtdb_put(rtdb,'rt_tddft:prop_method',
     $           mt_int,1,prop_method))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)



C
C     Exponentiation method
C
         elseif (inp_compare(.false.,test,'exp')) then
            
            if (.not.inp_a(method)) then
               call errquit (pname//
     $              "failed to parse exponentiation method",0,0)
            endif

            if (method.eq."diag") then
               exp_method = 2  ! diagonalization

            elseif (method.eq."pseries") then
               exp_method = 1  ! power series
               
            elseif (method.eq."magnus") then
               exp_method = 3   ! baker-campbell-hausdorff (disabled) !!since assumes e^X A e^-X
               call errquit (pname//"BCH exp disabled"//method)
               
            else
               call errquit (pname//
     $              "invalid exponentiation method: "//method)
            endif
            
            if (.not.rtdb_put(rtdb,'rt_tddft:exp_method',
     $           mt_int,1,exp_method))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            


C
C     profiling
C
         elseif (inp_compare(.false.,test,'prof')) then
            prof=.true.

            if (.not.rtdb_put(rtdb,'rt_tddft:prof',mt_log,1,prof))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     profiling
C
         elseif (inp_compare(.false.,test,'matrix_checks')) then
            matrix_checks=.true.

            if (.not.rtdb_put(rtdb,'rt_tddft:matrix_checks',
     $           mt_log,1,matrix_checks))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C     
C     Override DFT settings so that there is *absolutely* no I/O.  If
C     "usedisk" is present, we will use DFT rtdb options, otherwise we
C     disable everything (default).
C     
         elseif (inp_compare(.false.,test,'usedisk')) then
            call errquit (pname//"usedisk disabled",0,0)

c$$$            nodisk=.false.
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:nodisk',mt_log,1,nodisk))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            

C
C     post process using dplot? (i.e. make density plots)
C
c$$$         elseif (inp_compare(.false.,test,'dplot')) then
c$$$            dplot_do=.true.
c$$$            
c$$$            if (inp_a(dplot_opts)) then
c$$$               if (trim(dplot_opts).eq."subgs") then
c$$$                  subgs = .true.
c$$$                  if (.not.rtdb_put(rtdb,'rt_tddft:subgs',
c$$$     $                 mt_log,1,subgs))
c$$$     $                 call errquit(pname//'Write failed to rtdb',
c$$$     $                 0,RTDB_ERR)
c$$$               else
c$$$                  call errquit (pname//"invalid dplot option: "
c$$$     $                 //trim(dplot_opts), 0, 0)
c$$$               endif
c$$$            endif
c$$$
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:dplot', mt_log,1,dplot_do))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c$$$
            
            
C
C     noprop directive (dont propagate)
C
         elseif (inp_compare(.false.,test,'noprop')) then
            noprop=.true.
            if (.not.rtdb_put(rtdb,'rt_tddft:noprop',mt_log,1,noprop))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     static directive (never rebuild Fock matrix)
C
         elseif (inp_compare(.false.,test,'static')) then
            static = .true.
            if (.not.rtdb_put(rtdb,'rt_tddft:static',mt_log,1,static))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            
            
C
C     Rules for applying fields to geometries
C
         elseif (inp_compare(.false.,test,"excite")) then
            nexcites = nexcites + 1
            
            if (.not.rtdb_put(rtdb,'rt_tddft:nexcites',
     $           mt_int,1,nexcites))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            call rt_tddft_input_excite (rtdb, nexcites, nfields)


C
C     List of time-dependent system properties to print.
C
         elseif (inp_compare(.false.,test,'print')) then
            call rt_tddft_input_print (rtdb)

C     
C     end
C     
         else if (inp_compare(.false.,test,'end')) then
            done = .true.
         else
            call errquit(pname//'Unknown directive: '//trim(test),
     $           0, INPUT_ERR)
         endif
      enddo  !end main parsing loop

      

C      if (.not.got_tmin) call errquit(pname//"must supply tmin",0,0)
      if (.not.got_tmax) call errquit(pname//"must supply tmax",0,0)
      if (.not.got_dt) call errquit(pname//"must supply dt",0,0)


C     HARDCODED TMIN
C      tmin = -5.0*dt
      tmin = 0d0


      if (.not.rtdb_put(rtdb,'rt_tddft:tmin',mt_dbl,1,tmin))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Check that parameters are valid and compatible with each other.
C     
C      if (tmin.lt.0d0) call errquit(pname//"tmin must be > 0", 0, 0)
      if (tmax.lt.0d0) call errquit(pname//"tmax must be > 0", 0, 0)
      if (tmax.lt.tmin) call errquit(pname//"tmax must be > tmin", 0, 0)
      if (dt.lt.0d0) call errquit(pname//"dt must be > 0", 0, 0)
      
      end subroutine
      
c $Id$
