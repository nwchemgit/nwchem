c
c     == task raman driver ==
c
C> \ingroup task
C> @{
      logical function task_raman(rtdb)
C
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"
#include "eaf.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
c
      logical int_normalize
      logical hyperraman
      logical vroa
      logical status
      logical freq_done
      logical debug
      logical task_freq
      logical task_property
      logical task_energy
c
      external task_property
      external task_energy
      external task_freq
      external int_normalize
c
      integer rtdb, geom ! [input] handles for rtdb
      integer nat,nc ! [input] # of atoms
      integer steps ! # of steps used in finite differencing
      integer rmmodes ! # of raman active modes
      integer l_freq, k_freq
      integer l_nmode, k_nmode
      integer k_m, l_m
      integer fd
      integer nfreq, l_fr, k_fr, response_order, i, first, last
      data nfreq/0/  ! number of frequencies, static if zero
c
      character*255 filename 
      character*32 theory
      character*16 plot ! normal or resonance
      character*16 line !  lorentzian (l) or gaussian (g) lineshape
c
      double precision offset
      double precision low, high
      double precision width ! full-width at half maximum (FWHM) in 1/cm
      double precision cpu, wall
      double precision step_size ! multiplictive factor for step along normal mode
      double precision laser(1)
c
      character*32 pname
c
      pname = 'task_raman: ' 
c
      debug =  ( .false. .and. ga_nodeid().eq.0 )
c
c --------------------------------------------------------------------
c
      call ecce_print_module_entry('task raman')
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c --------------task:theory call and response--------------------------
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit(pname//'theory not specified',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit(pname//'failed to invalidate status',0,
     &       RTDB_ERR)
      if (ga_nodeid().eq.0 .and.
     $    util_print('task_raman', print_low)) then
        write(LuOut,*)
        call util_print_centered(LuOut,'Raman Analysis',40,.true.)
        write(LuOut,*)
      endif
c
c ------First Hessian(frequencies) Either read or calculate--------------------
c
      if (.not. rtdb_get(rtdb,'raman:freq_done',mt_log,1,freq_done))
     $      freq_done=.false.
      status = .true.
      if (.not.freq_done) status = task_freq(rtdb)
      if (.not.status) then
         call errquit(pname//'task_freq failed',911, CALC_ERR)
      else
         if (.not. rtdb_put(rtdb,'raman:freq_done',mt_log,1,freq_done))
     $   call errquit(pname//'rtdb_put freq_done',911, RTDB_ERR)
      endif
c
c --------Create/load reference geometry to get the number of atoms------------    
c
      if (.not.geom_create(geom,'geometry')) call errquit
     $  (pname//'geom_create failed?',1, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $  call errquit(pname//'geom_rtdb_load failed?',2, RTDB_ERR)
      if (.not. geom_ncent(geom,nat)) 
     &  call errquit(pname//'geom_ncent failed?',3, GEOM_ERR)
      nc = nat*3
      rmmodes = nc-6
c
c -------- get masses for weighting coords ------------------------------------
c
      if (.not. ma_push_get(mt_dbl, nat, 'mass', l_m,k_m))
     & call errquit(pname//'mass not allocated',nat,MA_ERR)
      if (.not.geom_masses_get(geom,nat,dbl_mb(k_m)))
     & call errquit(pname//'geom_get_mass failed',911, GEOM_ERR)
      if (.not. geom_destroy(geom)) 
     & call errquit(pname//'geom_destroy failed?',911, GEOM_ERR)
c
c ----------Load parameters from rtdb or assign defaults-----------------------    
c
      if (.not. rtdb_cget(rtdb,'raman:plot',1,plot))
     $      plot = 'normal' ! set default "normal"
      if (.not. rtdb_cget(rtdb,'raman:line',1,line))
     $      line = 'lorentzian' !  lorentzian (l) or gaussian (g) lineshape
      if (.not. rtdb_get(rtdb,'raman:width',mt_dbl,1,width))
     $      width = 20.0D+00 ! full-width at half maximum (FWHM) in 1/cm
      if (.not. rtdb_get(rtdb,'raman:step_size',mt_dbl,1,step_size))
     $      step_size = 0.010D+00 ! srtep size for displacement along modes
      if (.not. rtdb_get(rtdb,'raman:steps',mt_int,1,steps))
     $      steps = 2 ! # steps to use for FD
      if (.not. rtdb_get(rtdb,'raman:low',mt_dbl,1,low))
     $      low  = 0.0D+00 ! lowest wavenumber  normal mode to use
      if (.not. rtdb_get(rtdb,'raman:high',mt_dbl,1,high))
     $      high = 100000.0D+00 ! Highest wavenumber normal mode to use
      if (.not. rtdb_get(rtdb,'raman:first',mt_int,1,first))
     $      first = 7 ! first normal mode to use
      if (.not. rtdb_get(rtdb,'raman:last',mt_int,1,last))
     $      last = 10000 ! last normal mode to use
      if (.not. rtdb_get(rtdb,'raman:hyperraman',mt_log,1,hyperraman))
     $      hyperraman = .false. ! # flag to calculate hyperraman
      if (.not. rtdb_get(rtdb,'raman:vroa',mt_log,1,vroa))
     $      vroa = .false. ! # flag to calculate vibrational 
      if (hyperraman) then
         first = 0
         rmmodes = nc
c
c --- in case we want overide the defaults for modes to include ---
c
         if (.not. rtdb_get(rtdb,'raman:first',mt_int,1,first))
     $      first = 7 ! srtep size for displacement along modes
      endif 
c
c ----------alocate space for freq and normal modes----------------------------
c
      if (.not. ma_push_get(mt_dbl, nc,'frequencies',l_freq,k_freq))
     &  call errquit(pname//'could not allocate l_freq',nc,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nc*nc,'normal modes',l_nmode,k_nmode))
     &  call errquit(pname//'could not allocate l_nmode',nc*nc,MA_ERR)
c
c --------attempt to read in frequencies and normal modes---------------------- 
c
      if (ga_nodeid().eq.0) then
        call util_file_name('nmode',.false.,.false.,filename)
        if (eaf_open(filename,eaf_rw,fd).ne.0)
     &    call errquit(pname//'eaf_open failed',555, DISK_ERR)
        offset = 16 
        if (eaf_read(fd,offset,dbl_mb(k_freq),nc*8).ne.0)
     &    call errquit(pname//'error reading frequencies',555,
     &       DISK_ERR)
        offset = offset + nc*8
        if (eaf_read(fd,offset,dbl_mb(k_nmode),nc*nc*8).ne.0)
     &    call errquit(pname//'error reading modes',555, DISK_ERR)
        if (eaf_close(fd).ne.0)
     &    call errquit(pname//'file close problem',555, DISK_ERR)
      endif  
c
c -------------grab laser freq(s) from rtdb---------------------------
c          (blatanly stolen from aoresponse_driver.F)
c
      if (.not.rtdb_get(rtdb, 'aoresponse:nfreq', mt_int, 1, nfreq))
     &  nfreq = 0
      if (debug) write (luout,*) 'nfreq',nfreq
c
c     at present we will tie the response order to the number
c     of frequencies that the user inputs, and calculate the 
c     response only for one set of frequencies. i.e. we calculate
c     alpha(-sum_i \omega_i; omega_1, omega_2, omega_3, ...).
c
      if (nfreq .eq. 0) then
        response_order = 1
      elseif(nfreq.gt.0) then
        response_order = nfreq
      else
        call errquit (pname//'nfreq .lt. zero',nfreq,RTDB_ERR)
      end if
c
c     exit if response order is larger than one (TEMPORARY)
      if (response_order.gt.1)
     &   call errquit (pname//'order>1',response_order,INPUT_ERR)
c
c     allocate memory for frequencies list and read values 
c     into dbl_mb(k_fr). If nfreq = 0, assign a single frequency = 0
c
      if (nfreq.ne.0) then
c       dynamic response: get frequencies
        if (.not. ma_push_get(mt_dbl,nfreq,'freqs',l_fr,k_fr)) 
     &    call errquit(pname//'ma_push_get failed k_fr',l_fr,MA_ERR)
        if (.not. rtdb_get(rtdb, 'aoresponse:freqlist', mt_dbl, nfreq,
     &     dbl_mb(k_fr)))
     &     call errquit(pname//'freqlist failed', 0, RTDB_ERR)
      else
c       static response: set all frequencies to zero
        if (.not. ma_push_get(mt_dbl,response_order,'freqs',l_fr,k_fr)) 
     &     call errquit(pname//'ma_push_get failed k_fr',0,MA_ERR) 
        do i=1, response_order
          dbl_mb(k_fr + (i-1)) = 0d0
        end do
      endif                     ! nfreq.ne.0
c
c     force nfreq = 1 if response order equals 1 (no freq. input)
      if (response_order.eq.1 .and. nfreq.ne.1) nfreq = 1
      if (debug) write (luout,*) 'freqs',dbl_mb(k_fr)
      nfreq = 1
      laser(1) = 0.088559D+00  ! 514 nm
      if (debug) write (luout,*) "freq",laser(1)
c
c ------------allocate freq and nmode and broadcast from node 0----------------
c
      call ga_brdcst(msg_frq,dbl_mb(k_freq),
     &               nc*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(msg_eigv,dbl_mb(k_nmode),
     &               nc*nc*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(msg_frq,dbl_mb(k_fr),
     &               nfreq*ma_sizeof(mt_dbl,1,mt_byte),0)
c
c ------------enough setup really do the calculation------------------------
c
      call task_raman_doit(rtdb,geom,nc,nat,first,last,rmmodes,steps,
     &     nfreq,plot,line,width,step_size,low,high,dbl_mb(k_nmode),
     &     dbl_mb(k_freq),dbl_mb(k_m),dbl_mb(k_fr) )
c
c -------------get rid of extra memory-----------------------------------------
c
      if (.not.ma_chop_stack(l_m)) 
     &  call errquit(pname//'ma_chop_stack failed k_m',l_m,MA_ERR)
c
c ---------------collect and print timings data--------------------------------
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - cpu
c
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit(pname//'failed storing cputime',0,
     &     MA_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit(pname//'failed storing walltime',0,
     &     MA_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit(pname//'failed to set status',0,
     &     MA_ERR)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      task_raman = status
      if (task_raman) then
        call ecce_print_module_exit('task raman', 'ok')
      else
        call ecce_print_module_exit('task raman', 'failed')
      endif
c
      return
      end
C> @}
c
c     == perform raman calculation ==
c
      subroutine task_raman_doit(rtdb,geom,nc,nat,first,last,rmmodes,
     &         steps,nfreq,plot,line,width,step_size,low,high,eigenvecs,
     &         eigenvals,mass,laser )
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"
#include "eaf.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "bas.fh"
c
      logical status
      logical hyperraman ! flag to calculate hyperraman (currently not supported)
      logical debug
      logical task_energy
      logical task_property
      logical int_normalize
c
      external int_normalize
      external task_energy
      external task_property
c
      integer rtdb, geom ! [input] handles for rtdb
      integer nat! [input] # of atoms
      integer nc ! [input] # of coordinates
      integer i, ii, iii ! countine indexes
      integer rmmodes ! # of raman active modes
      integer steps ! [input] # of steps used in Finited differencing
      integer l_coord, k_coord
      integer nfreq ! [input] number of frequencies from property
      integer first ! [input] first mode to calculate for aoresponse
      integer last  ! [input] last mode to calculate for aoresponse
      integer begin ! first=begin if not restarting, else modified
      integer start ! [input] first mode to use in calculation of plot
c
      character*16 plot ! normal or resonance
      character*16 line !  lorentzian (l) or gaussian (g) lineshape
c
      double precision rminfo(rmmodes,4) ! data for raman spec
      double precision width ! [input] full-width at half maximum (FWHM) in 1/cm
      double precision eigenvecs(nc,nc) ! [input](xyz&atom,mode)
      double precision eigenvals(nc)      ! [input] (mode)
      double precision mass(nat)      ! [input] (masses)
      double precision step_size ! multiplictive factor for step along normal mode
      double precision master(3,nat) ! saved initial coordinates
      double precision coord(3,nat)  ! displaced coordinate passed to raman_modestep
      double precision junk3(6,3), junk4(6,3) ! temp. arrays for FD of polarizanbility
      double precision laser(nfreq) ! [input] frequency of laser
      double precision low, high ! wavenumber bounds for normal modes allowed in raman calc
c
      character*32 pname
c
      pname = 'task_raman_doit: ' 
      debug =  ( .false. .and. ga_nodeid().eq.0 )
c
      begin = 2
      start = 1
      call dfill(rmmodes*4,0.0D+00,rminfo,1)  ! initialize rminfo
      if (last .eq. 10000) then  ! if last is default set to the number of normal modes
        last = first + rmmodes - 1 
      endif
      if (.not. rtdb_get(rtdb, 'raman:rminfo ', mt_dbl,
     &    rmmodes*4, rminfo)) then
        if (ga_nodeid().eq.0) then
          write(LuOut,*)pname//"no rminfo found"
c         call errquit('aoresponse:failed to get rminfo', 0, RTDB_ERR)
        endif
      endif
      do i = 1,rmmodes
         if (rminfo(i,1) .le. low ) then
             goto 100
         endif
      enddo
      begin = i
 100  continue
      begin = first + begin - 2
      if (debug) then
              write(6,*)"rmmodes", rmmodes 
              write(6,*)"f l b", first, last, begin
              call output(rminfo,1,rmmodes,1,4,rmmodes,4,1)
      endif
c
c ------open/load/get coords/ get cartesian coordinates for master geom--------
c
      if (.not.geom_create(geom,'geometry')) 
     &  call errquit(pname//'geom create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) 
     &  call errquit(pname//'geom_rtdb_load failed',911, GEOM_ERR)
      if (.not.geom_cart_coords_get(geom,master)) 
     &  call errquit(pname//'geom_get_cart_coords failed',911, GEOM_ERR)
      if (.not.geom_destroy(geom)) 
     &  call errquit(pname//'geom_destroy failed',911, GEOM_ERR)
c
c -------------------actually do the deed.-------------------------------------
c
       if (begin .ge. last)  go to 200             
       do ii=begin,last ! 2*rmmodes in the end. shorten for tests.
        call dfill(18,0.0D+00,junk3,1)  ! real junk3 a real (FD)
        call dfill(18,0.0D+00,junk4,1)  ! real junk4 a imaginary (FD)
        if (rminfo(ii,1) .ge. high) then ! skips normal modes over a certain wavenumber
          last = ii -1  ! modify last from the high [input] for plot
          goto 200
        endif 
c
        do iii=1,steps
c
c      ----------------open/load/get coords/ get cartesian coordinates---------
c
         if (.not.geom_create(geom,'geometry')) 
     &      call errquit(pname//'geom create failed',911, GEOM_ERR)
         if (.not.geom_rtdb_load(rtdb,geom,'geometry')) 
     &      call errquit(pname//'geom_rtdb_load failed',911, GEOM_ERR)
         if (.not.geom_cart_coords_get(geom,coord)) 
     &      call errquit(pname//'geom_get_cart_coords failed',911, 
     &           GEOM_ERR)
c
c      -------------Step along normal mode ii ---------------------------------
c
          call raman_modestep (rtdb,nc,nat,geom,rmmodes,ii,iii,first,
     &    eigenvecs,eigenvals,coord,rminfo,step_size)
c      store displaced geom in rtb
         if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &     call errquit(pname//'geom_rtdb_store failed',5)
c      destroy geom
         if (.not. geom_destroy(geom))
     &     call errquit(pname//'geom corrupt',0, GEOM_ERR)
c
c      --------------------run properties (AORESPONSE) on displaced geom-------
c
       status = task_property(rtdb)
c
c      -------------------open/load/get coords/ get cartesian coordinates------
c
         if (.not.geom_create(geom,'geometry')) 
     &     call errquit(pname//'geom create failed',911, GEOM_ERR)
         if (.not.geom_rtdb_load(rtdb,geom,'geometry')) 
     &     call errquit(pname//'geom_rtdb_load failed',911, GEOM_ERR)
c      restore inital geom "master"
         if (.not. geom_cart_coords_set(geom,master))
     $     call errquit(pname//'bad geom',0, GEOM_ERR)
c     --- coord now stored in rtdb ---
         if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &     call errquit(pname//'geom_rtdb_store failed',5)
c      --- destroy geom ---
         if (.not. geom_destroy(geom))
     &   call errquit(pname//'geom corrupt',0,GEOM_ERR)
c      --- save response data from aoresponse
         call raman_save(rtdb, iii, junk3, junk4 )
        enddo ! iii
c
c ------------------------do finite differenced for mode ii--------------------
c
         call fd_raman(rtdb, ii, rmmodes, nat, nc, junk3, junk4,
     &           step_size, rminfo, eigenvecs, mass)
       enddo ! ii 
c
c ---------------------calculate ramam spectra---------------------------------
c
 200  continue ! come here on a restart to recalculate a plot / from the high variable
c
      if (ga_nodeid().eq.0 ) then
        call raman_scattering(rtdb,start,last,rmmodes,nfreq,plot,
     &                      line,width,laser,rminfo)
      end if
c
      return
      end
c $Id$
