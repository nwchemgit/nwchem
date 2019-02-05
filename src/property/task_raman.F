c
c     == task raman driver ==

      logical function task_raman(rtdb)
c
c Authors: Jonathan Mullin, Northwestern University (ver 1: Jan. 2011)
c          Fredy W. Aquino, Northwestern University (ver 2: Oct. 2012)

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

      logical int_normalize
      logical hyperraman
      logical vroa
      logical status
      logical freq_done
      logical debug
      logical task_freq
      logical task_property
      logical task_energy

      external task_property
      external task_energy
      external task_freq
      external int_normalize

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

      character*255 filename 
      character*32 theory
      character*16 plot ! normal or resonance
      character*16 line !  lorentzian (l) or gaussian (g) lineshape

      double precision offset
      double precision low, high
      double precision width ! full-width at half maximum (FWHM) in 1/cm
      double precision cpu, wall
      double precision step_size ! multiplictive factor for step along normal mode
      double precision laser(1)
      double precision gamwidth
c      logical lifetime

      integer j,pos,first0 ! FA-06-15-12
      logical preraman     ! FA-06-18-12
      logical linear
      logical lgarbage

      character*32 pname

      pname = 'task_raman: ' 

      debug =  ( .false. .and. ga_nodeid().eq.0 )
c
c --------------------------------------------------------------------

      call ecce_print_module_entry('task raman')
      cpu  = util_cpusec()
      wall = util_wallsec()

c
c --------------task:theory call and response--------------------------

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

      if (.not. rtdb_get(rtdb,'raman:freq_done',mt_log,1,freq_done))
     $      freq_done=.false.

c     if (ga_nodeid().eq.0) then
c      write(*,*) 'FA-raman-2 freq_done=',freq_done
c     endif

      status = .true.
      if (.not.freq_done) then
        status = task_freq(rtdb)
      endif
      if (.not.status) then
         call errquit(pname//'task_freq failed',911, CALC_ERR)
      else
         if (.not. rtdb_put(rtdb,'raman:freq_done',mt_log,1,freq_done))
     $   call errquit(pname//'rtdb_put freq_done',911, RTDB_ERR)
      endif
c
c --------Figure out if molecule is linear------------    

c     if vib module doesn't list molecule as linear, assume it is not
      if (.not. rtdb_get(rtdb,'vib:linear',mt_log,1,linear))
     $      linear=.false.
c
c --------Create/load reference geometry to get the number of atoms------------    

      if (.not.geom_create(geom,'geometry')) call errquit
     $  (pname//'geom_create failed?',1, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $  call errquit(pname//'geom_rtdb_load failed?',2, RTDB_ERR)
      if (.not. geom_ncent(geom,nat)) 
     &  call errquit(pname//'geom_ncent failed?',3, GEOM_ERR)
      nc = nat*3
      if (linear) then
        rmmodes = nc-5
      else
        rmmodes = nc-6
      end if

c      if (ga_nodeid().eq.0) then
c       write(*,1) nat,nc,rmmodes
c 1     format('FA-(nat,nc,rmmodes)=(',i4,',',i4,',',i4,')')
c      endif
c
c -------- get masses for weighting coords ------------------------------------
      if (.not. ma_push_get(mt_dbl, nat, 'mass', l_m,k_m))
     & call errquit(pname//'mass not allocated',nat,MA_ERR)
      if (.not.geom_masses_get(geom,nat,dbl_mb(k_m)))
     & call errquit(pname//'geom_get_mass failed',911, GEOM_ERR)
      if (.not. geom_destroy(geom)) 
     & call errquit(pname//'geom_destroy failed?',911, GEOM_ERR)
c
c ----------Load parameters from rtdb or assign defaults-----------------------    
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
      if (.not. rtdb_get(rtdb,'raman:first',mt_int,1,first)) then
            if (linear) then
              first = 6 ! first normal mode to use
            else
              first = 7 ! first normal mode to use
            end if
      end if
      if (.not. rtdb_get(rtdb,'raman:last',mt_int,1,last))
     $      last = 10000 ! last normal mode to use
      if (.not. rtdb_get(rtdb,'raman:hyperraman',mt_log,1,hyperraman))
     $      hyperraman = .false. ! # flag to calculate hyperraman
      if (.not. rtdb_get(rtdb,'raman:vroa',mt_log,1,vroa))
     $      vroa = .false. ! # flag to calculate vibrational 
      if (.not. rtdb_get(rtdb,'raman:preraman',mt_log,1,preraman))
     $      preraman = .false. ! # flag to do task_freq() and leave
      if (linear) then
        first0=6 ! constant
      else
        first0=7 ! constant
      end if
c
c     force lifetime/damping when plot=resonance
c     default damping value 0.004 (if not assigned)
c
      if (plot.eq.'resonance') then
         if (.not.rtdb_put(rtdb, 'aoresponse:damping', mt_int, 1,0)) 
     &        call errquit('aor_driver: damping failed ',0, RTDB_ERR)
         gamwidth=0d0
         lgarbage=rtdb_get(rtdb, 'aoresponse:gamwidth',mt_dbl, 1,
     G        gamwidth)
        if(abs(gamwidth).lt.1d-6) then
           gamwidth=0.004d0
           if(ga_nodeid().eq.0) write(6,1964) gamwidth
 1964 format(' Resonance Raman: Using a nonzero Damping value of',f10.4)
           if (.not.rtdb_put(rtdb, 'aoresponse:gamwidth',
     &          mt_dbl, 1,gamwidth)) 
     &          call errquit('aor_driver: gamwidth failed ',0, RTDB_ERR)
        endif
      endif

c ======== FA-debug =============== START
c      if (ga_nodeid().eq.0) then
c       write(*,2) plot,line,width,step_size,steps
c 2     format('FA-(plot,line,width,step_size,steps)=(',
c     &         a,',',a,',',f15.8,',',f15.8,',',i4,')')
c       write(*,22) low,high,first,last,hyperraman,vroa,preraman
c 22    format('FA-(low,high,first,last,hyperraman,vroa,preraman)=(',
c     &        f15.8,',',f15.8,',',i8,',',i8,',',L,',',L,',',L,')')
c      endif 
c ======== FA-debug =============== END
      if (hyperraman) then
         first = 0
         rmmodes = nc
c
c --- in case we want overide the defaults for modes to include ---
         if (.not. rtdb_get(rtdb,'raman:first',mt_int,1,first)) then
            if (linear) then
              first = 6 ! srtep size for displacement along modes
            else
              first = 7 ! srtep size for displacement along modes
            end if
         end if
      endif 
c
c ----------alocate space for freq and normal modes----------------------------

      if (.not. ma_push_get(mt_dbl, nc,'frequencies',l_freq,k_freq))
     &  call errquit(pname//'could not allocate l_freq',nc,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nc*nc,'normal modes',l_nmode,k_nmode))
     &  call errquit(pname//'could not allocate l_nmode',nc*nc,MA_ERR)
c
c --------attempt to read in frequencies and normal modes---------------------- 

      if (ga_nodeid().eq.0) then
        call util_file_name('nmode',.false.,.false.,filename)
    
c        write(*,3) filename
c 3      format('FA-Reading filename=',a)

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
c ============== FA-debug ==================== START
c      if (ga_nodeid().eq.0) then
c       write(*,*) 'Listing frequencies------- START'
c       do i=1,nc
c        write(*,44) i,dbl_mb(k_freq+i-1)
c 44     format('freq(',i4,')=',f15.8)
c       enddo ! end-loop-i
c       write(*,*) 'Listing frequencies------- END'
c       write(*,*) 'Listing modes------- START'
c       pos=1
c       do i=1,nc
c         do j=1,nc
c           write(*,4) i,j,dbl_mb(k_nmode+pos-1)
c  4        format('mode('i4,',',i4,')=',f15.8)
c           pos=pos+1
c         enddo ! end-loop-j
c       enddo ! end-loop-i
c       write(*,*) 'Listing modes------- END'
c      endif 
c ============== FA-debug ==================== END
c
c -------------grab laser freq(s) from rtdb---------------------------
c          (blatanly stolen from aoresponse_driver.F)

      if (.not.rtdb_get(rtdb, 'aoresponse:nfreq', mt_int, 1, nfreq))then
         nfreq=0
         if (.not. rtdb_put(rtdb, 'aoresponse:nfreq', mt_int, 1,0))
     $     call errquit('task_raman: rtdb_put failed', 0, RTDB_ERR)
      endif
      if (debug) write (luout,*) 'nfreq',nfreq

c
c     at present we will tie the response order to the number
c     of frequencies that the user inputs, and calculate the 
c     response only for one set of frequencies. i.e. we calculate
c     alpha(-sum_i \omega_i; omega_1, omega_2, omega_3, ...).

      if (nfreq .eq. 0) then
        response_order = 1
      elseif(nfreq.gt.0) then
        response_order = nfreq
      else
        call errquit (pname//'nfreq .lt. zero',nfreq,RTDB_ERR)
      end if
c
c     allocate memory for frequencies list and read values 
c     into dbl_mb(k_fr). If nfreq = 0, assign a single frequency = 0

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
c     exit if response order is larger than one (TEMPORARY)
      if (response_order.gt.1)
     &   call errquit (pname//'order>1',response_order,INPUT_ERR)
c     force nfreq = 1 if response order equals 1 (no freq. input)
      if (response_order.eq.1 .and. nfreq.ne.1) nfreq = 1
      if (debug) write (luout,*) 'freqs',dbl_mb(k_fr)
      nfreq = 1
      laser(1) = 0.088559D+00  ! 514 nm
      if (debug) write (luout,*) "freq",laser(1)
c ============== FA-debug ==================== START
c     if (ga_nodeid().eq.0) then
c      do i=1,nfreq
c       write(*,6) i,dbl_mb(k_fr+i-1)
c6      format('FA-freq(',i4,')=',f15.8)
c      enddo
c     endif
c ============== FA-debug ==================== END
c
c ------------allocate freq and nmode and broadcast from node 0----------------

      call ga_brdcst(msg_frq,dbl_mb(k_freq),
     &               nc*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(msg_eigv,dbl_mb(k_nmode),
     &               nc*nc*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(msg_frq,dbl_mb(k_fr),
     &               nfreq*ma_sizeof(mt_dbl,1,mt_byte),0)
c
c ------------enough setup really do the calculation------------------------
      if (.not.preraman) then
      call task_raman_doit(rtdb,geom,nc,nat,
     &                     first0, ! = 6 or 7
     &                     first,last,rmmodes,
     &                     steps,nfreq,plot,line,width,
     &                     step_size,
     &                     low,high,
     &                     dbl_mb(k_nmode),dbl_mb(k_freq),
     &                     dbl_mb(k_m),dbl_mb(k_fr))
      endif
c -------------get rid of extra memory-----------------------------------------

      if (.not.ma_chop_stack(l_m)) 
     &  call errquit(pname//'ma_chop_stack failed k_m',l_m,MA_ERR)
c
c ---------------collect and print timings data--------------------------------

      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - cpu
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit(pname//'failed storing cputime',0,
     &     MA_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit(pname//'failed storing walltime',0,
     &     MA_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit(pname//'failed to set status',0,
     &     MA_ERR)
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
c
c     == perform raman calculation ==
      subroutine task_raman_doit(rtdb,geom,nc,nat,
     &                           first0, ! = 6 or 7
     &                           first,last,rmmodes,
     &                           steps,nfreq,
     &                           plot,line,width,
     &                           step_size,
     &                           low,high,
     &                           eigenvecs,eigenvals,
     &                           mass,laser)
c
c Authors: Jonathan Mullin, Northwestern University (ver 1: Jan. 2011)
c          Fredy W. Aquino, Northwestern University (ver 2: Oct. 2012)

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

      logical status
      logical hyperraman ! flag to calculate hyperraman (currently not supported)
      logical debug
      logical ramanspc,! to store aoresp data in files
     &        ramanlist! to read rminfo from files and create list of missing modes            
      integer rtdb,geom! [input] handles for rtdb
      integer nat      ! [input] # of atoms
      integer nc       ! [input] # of coordinates
      integer i,ii,iii ! countine indexes
      integer rmmodes  ! # of raman active modes
      integer steps    ! [input] # of steps used in Finited differencing
      integer l_coord, 
     &        k_coord
      integer nfreq    ! [input] number of frequencies from property
      integer first,   ! [input] first mode to calculate for aoresponse
     &        first0  
      integer last     ! [input] last mode to calculate for aoresponse
      integer begin    ! first=begin if not restarting, else modified
      integer start    ! [input] first mode to use in calculation of plot
      integer mode_ini,! window of normal modes to process
     &        mode_end
      character*16 plot ! normal or resonance
      character*16 line !  lorentzian (l) or gaussian (g) lineshape
      double precision rminfo(rmmodes,4)! data for raman spec
      double precision width            ! [input] full-width at half maximum (FWHM) in 1/cm
      double precision eigenvecs(nc,nc) ! [input] (xyz&atom,mode)
      double precision eigenvals(nc)    ! [input] (mode)
      double precision mass(nat)        ! [input] (masses)
      double precision step_size        ! factor for step along normal mode
      double precision master(3,nat)    ! saved initial coordinates
      double precision coord(3,nat)     ! displaced coordinate passed to raman_modestep
      double precision junk3(6,3),      ! temp. arrays for FD of polarizability   
     &                 junk4(6,3)       ! temp. arrays for FD of polarizability
      double precision laser(nfreq)     ! [input] frequency of laser
      double precision low, high        ! wavenumber bounds for normal modes allowed in raman calc
      character*32 pname
c --------- FA vars ---------- START
      character*5 lbl_raman
      parameter(lbl_raman='raman')
      character*255 lbl_raman1
      integer nmodes_missing,
     &        modes_missing(rmmodes),
     &        ignore_missingmode 
      character*255 fname_raman
      integer modecc,cmode,nmodes_slc,
     &        l_rmodes,k_rmodes,ramantype
      logical raman_modes_list,raman_done
      integer unitno
      parameter (unitno = 77)
      logical raman_write,
     &        raman_read_list  
      logical oldmode
      external raman_write,
     &         raman_read_list
c --------- FA vars ---------- END
      raman_done=.false.
      pname = 'task_raman_doit: ' 
      debug =  ( .false. .and. ga_nodeid().eq.0 )

      begin = 2
      start = 1
      call dfill(rmmodes*4,0.0D+00,rminfo,1)  ! initialize rminfo
       if (last .eq. 10000) then  ! if last is default set to the number of normal modes
         last = first + rmmodes - 1 
       endif
      mode_ini=first-first0+1 ! in : mode_ini
      mode_end=last -first0+1 ! in : mode_end
      if (ga_nodeid().eq.0) then
       write(luout,1) mode_ini,mode_end
 1     format('(mode_ini,mode_end)=(',i5,',',i5,')')
      endif
      if (.not. rtdb_get(rtdb, 'raman:rminfo ', mt_dbl,
     &    rmmodes*4, rminfo))
     & write(6,*)pname//"no rminfo found"
c     &   call errquit('aoresponse:failed to get rminfo', 0, RTDB_ERR)
       do i = 1,rmmodes
         if (rminfo(i,1) .le. low ) then
             goto 100
         endif
       enddo
       begin = i
 100   continue
        begin = first + begin - 2

      if (debug) then
              write(6,*)"rmmodes", rmmodes 
              write(6,*)"f l b", first, last, begin
              call output(rminfo,1,rmmodes,1,4,rmmodes,4,1)
      endif
c
c ------open/load/get coords/ get cartesian coordinates for master geom--------

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
c ============== FA-debug ============= START
      if (ga_nodeid().eq.0) then
       write(luout,151) begin,last,steps
 151   format('(begin,last,steps)=(',i4,',',i4,',',i4,')')
      endif
c ---- ramanspc=.true. will save aoresponse data in files
      ramanspc=.false.
      status=rtdb_get(rtdb,'raman:aores0',mt_log,1,ramanspc)
c ============== FA-debug ============= END
c 0000000 Reading modes from input script 000000000 START
      ramantype=1 ! defaul value
      status=rtdb_get(rtdb,'raman:ramantype',
     &                mt_int,1,ramantype)

      select case (ramantype)
       case (1)  ! Default
                   if (ga_nodeid().eq.0)
     &               write(luout,*) 'ENTER raman default'
                   call raman_default(
     &             modecc,     ! in: mode counter
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             lbl_raman,  ! in: raman label
     &             begin,      ! in: 
     &             last,       ! in: 
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom)  ! in: (rtdb,geom) handles
       case (2)  ! raman list
                   if (ga_nodeid().eq.0)
     &               write(luout,*) 'ENTER raman list'
                   ignore_missingmode=0 ! default do nothing
                   status=rtdb_get(rtdb,'raman:ramanskipmode',
     &                             mt_int,1,ignore_missingmode)
                   call raman_inp_files_rmodes(
     &             modecc,     ! in: mode counter
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             lbl_raman,  ! in: raman label
     &             mode_ini,   ! in: 
     &             mode_end,   ! in: 
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom,  ! in: (rtdb,geom) handles
     &             ignore_missingmode) ! in: =1 ignore calc missing mode
       case (3)  ! raman modes list
                   if (ga_nodeid().eq.0)
     &               write(luout,*) 'ENTER raman modes list'
                   call raman_inp_list_rmodes(
     &             modecc,     ! in: mode counter
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             lbl_raman,  ! in: raman label
     &             begin,      ! in: starting mode
     &             last,       ! in: ending   mode
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom)  ! in: (rtdb,geom) handles
       case (:0) ! All other values
       case (4:)
                   call errquit(
     &                  'task_raman_doit: ramantype value not 1,2 or 3',
     $                  555, RTDB_ERR)        
      end select

      oldmode =  rtdb_parallel(.false.)
      if (ga_nodeid().eq.0 ) then
        write(luout,*) 'Showing rminfo:'
        write(luout,204) start,last,rmmodes,nfreq,plot,width
 204    format('(start,last,rmmodes,nfreq,plot,width)=(',
     &         i4,',',i4,',',i4,',',i4,',',a,',',f15.8,')')
        do i=1,rmmodes
         write(luout,201) i,
     &                rminfo(i,1),rminfo(i,2),
     &                rminfo(i,3),rminfo(i,4)         
 201     format('rminfo(',i4,')=(',
     &          f15.8,',',f15.8,',',f15.8,',',f15.8,')')
        enddo ! end-loop-i 
        if (ramantype.ne.3) then
         call raman_scattering(
     &                       rtdb,
     &                       mode_ini,
     &                       mode_end,
     &                       rmmodes,
     &                       nfreq,
     &                       plot,
     &                       line,
     &                       width,
     &                       laser,
     &                       rminfo)
        endif ! enter-only-if-raman_type=1,2
      endif
      oldmode =  rtdb_parallel(oldmode)
      return
      end
c 00000000000000000000000000000000000000000000000000000000000000
c 000000000000000 FA-06-19-12 extra raman routines 0000000 START
c 00000000000000000000000000000000000000000000000000000000000000
      subroutine raman_per_mode(
     &             ii,         ! in: ii-th raman mode
     &             modecc,     ! in: mode counter
     &             fname_raman,! in: filename to store raman data
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom,  ! in: (rtdb,geom) handles
     W       oprint_dalphadq, filedadq)
c
c Authors: Jonathan Mullin, Northwestern University (ver 1: Jan. 2011)
c          Fredy W. Aquino, Northwestern University (ver 2: Oct. 2012)
     
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
      logical status
      logical debug
      logical task_property
      logical ramanspc ! FA-06-18-12 to store aoresp data in files
      external task_property,
     &         raman_modestep,
     &         fd_raman,raman_save
      integer rtdb,geom! [input] handles for rtdb
      integer nat      ! [input] # of atoms
      integer nc       ! [input] # of coordinates
      integer i,ii,iii ! countine indexes
      integer rmmodes  ! # of raman active modes
      integer steps    ! [input] # of steps used in Finited differencing
      integer nfreq    ! [input] number of frequencies from property
      integer first0   ! [input] first mode to calculate for aoresponse
      double precision rminfo(rmmodes,4) ! data for raman spec
      double precision eigenvecs(nc,nc)  ! [input] (xyz&atom,mode)
      double precision eigenvals(nc)     ! [input] (mode)
      double precision mass(nat)         ! [input] (masses)
      double precision step_size         ! [input] multiplictive factor for step along normal mode
      character*(*) filedadq       ! [in]
      logical oprint_dalphadq      ! [in]
      double precision master(3,nat)     ! saved initial coordinates
      double precision coord(3,nat)      ! displaced coordinate passed to raman_modestep
      double precision junk3(6,3),       ! temp. arrays for FD of polarizanbility
     &                 junk4(6,3)        ! temp. arrays for FD of polarizanbility
      double precision alfare(3,3),
     &                 alfaim(3,3)
      character*32 pname
c --------- FA vars ---------- START
      character*255 lbl_ramanpol1
      character*255 ramanpolfile
      character*8 lbl_ramanpol
      parameter(lbl_ramanpol='ramanpol') 
      character*255 fname_raman
      integer modecc
      integer unitno
      parameter (unitno = 77)
      logical raman_write,ramanpol_read,status1
      external raman_write,ramanpol_read,
     &         raman_save1

c --------- FA vars ---------- END
        if (ramanspc) then
         if (.not. rtdb_put(rtdb,'raman:aores1',mt_int, 1,ii)) call ! FA-06-18-12
     $     errquit('task_raman: failed to put ii',0,RTDB_ERR)      
        endif
        call dfill(18,0.0D+00,junk3,1)  ! real junk3 a real (FD)
        call dfill(18,0.0D+00,junk4,1)  ! real junk4 a imaginary (FD)
        do iii=1,steps
          if (ramanspc) then  ! FA-06-18-12
           if (.not. rtdb_put(rtdb,'raman:aores2',
     &                        mt_int,1,iii)) 
     $      call  errquit('task_raman: failed to put iii',0,RTDB_ERR)      
          endif
c      ----------------open/load/get coords/ get cartesian coordinates---------
          if (.not.geom_create(geom,'geometry')) 
     &      call errquit(pname//'geom create failed'   ,911, GEOM_ERR)
          if (.not.geom_rtdb_load(rtdb,geom,'geometry')) 
     &      call errquit(pname//'geom_rtdb_load failed',911, GEOM_ERR)
          if (.not.geom_cart_coords_get(geom,coord)) 
     &      call errquit(pname//'geom_get_cart_coords failed',911, 
     &           GEOM_ERR)
c      -------------Step along normal mode ii ---------------------------------
          call raman_modestep(rtdb,nc,nat,geom,
     &                        rmmodes,ii,iii,first0,
     &                        eigenvecs,eigenvals,
     &                        coord,rminfo,step_size)
c      store displaced geom in rtb
         if (.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &     call errquit(pname//'geom_rtdb_store failed',5)
c      destroy geom
         if (.not. geom_destroy(geom))
     &     call errquit(pname//'geom corrupt',0, GEOM_ERR)
c 000000000 run properties (AORESPONSE) on displaced geom 000 START
         write(lbl_ramanpol1,'(a8,"_",i4.4,"-",i1)') 
     &                        lbl_ramanpol,
     &                        ii-first0+1,iii 
         call util_file_name(lbl_ramanpol1,
     &                      .false.,.false.,ramanpolfile)  
         if (ga_nodeid().eq.0)
     &    write(luout,*) 'ramanpolfile=',ramanpolfile

         status1=ramanpol_read(
     &           ramanpolfile,! in: filename
     &           ii,          ! in: mode-nr
     &           rmmodes,     ! in: total modes
     &           iii,         ! in: step counter
     &           alfare,      ! ou: RE of alpha (polarization)
     &           alfaim)      ! ou: IM of alpha (polarization)
         if (.not.status1) then ! 00000000 Calc-polarization-step-iii-th --START
          status = task_property(rtdb) 
          status = rtdb_parallel(.true.)
c  -------------get raman ploarizability from rtdb ----------------------------
          if (.not. rtdb_get(rtdb, 'raman:alfare ', mt_dbl, 9, alfare))
     &    call errquit('raman_save:failed to get alfare', 1, RTDB_ERR)
          if (.not. rtdb_get(rtdb, 'raman:alfaim ', mt_dbl, 9, alfaim)) 
     &    call errquit('raman_save:failed to get alfaim', 2, RTDB_ERR)
         endif                  ! 00000000 Calc-polarization-step-iii-th --END
c 000000000 run properties (AORESPONSE) on displaced geom 000 END
c --------open/load/get coords/ get cartesian coordinates------
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

c         call raman_save(rtdb, iii, junk3, junk4 )

         call raman_save1(
     &            alfare,      ! in: RE alpha, polarization
     &            alfaim,      ! in: IM alpha, polarization
     &            ii-first0+1, ! in: mode nr
     &            rmmodes,     ! in: total raman modes
     &            iii,         ! in: step index, = 1,2 for steps=2
     &            ramanpolfile,! in: filename to store polariz tensor
     &            junk3,       ! ou: two sets of polarization stored RE
     &            junk4)       ! ou: two sets of polarization stored IM

        enddo ! iii
c ------------------------do finite defernced for mode ii----------------------
         call fd_raman(rtdb,ii,rmmodes,
     &                 nat,nc,
     &                 junk3, ! in : (+/-) polarization tensor RE (6,3)
     &                 junk4, ! in : (+/-) polarization tensor IM (6,3)
     &                 step_size,
     &                 rminfo,
     &                 eigenvecs,
     &                 mass,
     &       first0,
     W       oprint_dalphadq, filedadq)
c ======== FA: Writing to file rminfo ========= START
c         if (ga_nodeid().eq.0)
c     &     write(*,*) 'BEF raman_write() ...'
         status=raman_write(
     &           fname_raman,! in: filename
     &           ii-first0+1,! in: mode-nr
     &           modecc,     ! in: mode counter
     &           rmmodes,    ! in: total modes
     &           rminfo)     ! in: raman info (to store only [ii-6]th mode)
       modecc=modecc+1 ! counter
c ======== FA: Writing to file rminfo ========= END
      return
      end
c $Id$
c 00000000000000000000000000000000000000000000000000000000000000
c 000000000000000 FA-06-19-12 extra raman routines 0000000 END
c 00000000000000000000000000000000000000000000000000000000000000
      subroutine raman_inp_list_rmodes(
     &             modecc,     ! in: mode counter
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             lbl_raman,  ! in: raman label
     &             begin,      ! in: starting mode
     &             last,       ! in: ending   mode
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom)  ! in: (rtdb,geom) handles
c
c Author: Fredy W. Aquino, Northwestern University (Oct. 2012)

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
      logical status
      logical debug
      logical ramanspc ! FA-06-18-12 to store aoresp data in files
      integer rtdb,geom! [input] handles for rtdb
      integer nat      ! [input] # of atoms
      integer nc       ! [input] # of coordinates
      integer i,ii,iii ! countine indexes
      integer nmodes_slc,
     &        cmode,
     &        l_rmodes,k_rmodes,
     &        begin,last
      integer rmmodes  ! # of raman active modes
      integer steps    ! [input] # of steps used in Finited differencing
      integer nfreq    ! [input] number of frequencies from property
      integer first0   ! [input] first mode to calculate for aoresponse
      double precision rminfo(rmmodes,4) ! data for raman spec
      double precision eigenvecs(nc,nc)  ! [input] (xyz&atom,mode)
      double precision eigenvals(nc)     ! [input] (mode)
      double precision mass(nat)         ! [input] (masses)
      double precision step_size         ! [input] multiplictive factor for step along normal mode
      double precision master(3,nat)     ! saved initial coordinates
      double precision coord(3,nat)      ! displaced coordinate passed to raman_modestep
      character*32 pname
      character*(*) lbl_raman
      character*255 lbl_raman1
c --------- FA vars ---------- START
      character*255 fname_raman
      integer modecc
      character*255 filedadq
      logical oprint_dalphadq

       status=rtdb_get(rtdb,'raman:nmodesslc',
     &                mt_int,1,nmodes_slc)   
        if (.not. ma_push_get(mt_int,nmodes_slc,'rmodes',
     &                        l_rmodes,k_rmodes))
     &   call errquit('task_raman ma_push_get fail k_rmodes',
     &                0,MA_ERR)
       status=rtdb_get(rtdb,'raman:modeslist',
     &                mt_int,nmodes_slc,int_mb(k_rmodes)) 
 
       do i=1,nmodes_slc
        if (ga_nodeid().eq.0) then
         write(luout,21) i,int_mb(k_rmodes+i-1)
 21      format('raman-mode-slc(',i4,')=',i4)
        endif
        if (.not.(int_mb(k_rmodes+i-1).ge. begin .and.
     &            int_mb(k_rmodes+i-1).le. last)) then
         write(luout,42) i,int_mb(k_rmodes+i-1),begin,last
 42      format('Error: raman-mode-slc(',i4,')=',
     &          i4,' not in range=[',
     &          i4,',',i4,']')
         call errquit('prop_input-raman: rtdb_put failed',
     $                555, RTDB_ERR)
        endif
       enddo ! end-loop-i
c ======== construct filename in listing mode ====== START
        write(lbl_raman1,'(a5,"_listslc_",i4.4,"-",i4.4)') 
     &                    lbl_raman,1,nmodes_slc       
        call util_file_name(lbl_raman1,
     &                     .false.,.false.,fname_raman)
        if (ga_nodeid().eq.0)
     &   write(luout,131) fname_raman
 131     format('Generating fname_raman1=',a)         
c ======== construct filename in listing mode ====== END
         oprint_dalphadq=util_print('raman_dalphadq', print_high).and.
     c        ga_nodeid().eq.0
      if(oprint_dalphadq) then
         call util_file_name('dadq', .false., .false.,filedadq)
         call raman_open_dadq(filedadq)
      endif
        modecc=1 ! reset counter         
        do cmode=1,nmodes_slc
         ii=int_mb(k_rmodes+cmode-1)
         if (ga_nodeid().eq.0) then
          write(luout,173) ii,rminfo(ii,1)
 173      format('(ii,rminfo(ii,1))=(',
     &            i5,',',f15.8,')')
         endif
         if (ga_nodeid().eq.0)
     &    write(luout,*) 'Calculating slc-raman mode ii=',ii

         call raman_per_mode(
     &             ii,         ! in: ii-th raman mode
     &             modecc,     ! in: mode counter
     &             fname_raman,! in: filename to store raman data
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom,  ! in: (rtdb,geom) handles
     &        oprint_dalphadq ,
     F        filedadq)
        enddo ! end-loop-cmode
      if(oprint_dalphadq) then
         call raman_close_dadq(filedadq)
      endif
        if (.not.ma_pop_stack(l_rmodes)) call
     &        errquit('task_raman: ma_pop_stack failed k_rmodes',
     &                0,MA_ERR)   
      return
      end

      subroutine raman_inp_files_rmodes(
     &             modecc,     ! in: mode counter
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             lbl_raman,  ! in: raman label
     &             mode_ini,   ! in: 
     &             mode_end,   ! in: 
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom,  ! in: (rtdb,geom) handles
     &             ignore_missingmode) ! in: =1 ignore calc missing mode
c
c Author: Fredy W. Aquino, Northwestern University (Oct. 2012)

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
      logical status
      logical debug
      logical ramanspc ! FA-06-18-12 to store aoresp data in files
      integer rtdb,geom! [input] handles for rtdb
      integer nat      ! [input] # of atoms
      integer nc       ! [input] # of coordinates
      integer i,ii,iii ! countine indexes
      integer mode_ini,mode_end
      integer rmmodes  ! # of raman active modes
      integer nmodes_slc,
     &        cmode,
     &        nmodes_missing,
     &        modes_missing(rmmodes),
     &        ignore_missingmode
      integer steps    ! [input] # of steps used in Finited differencing
      integer nfreq    ! [input] number of frequencies from property
      integer first0   ! [input] first mode to calculate for aoresponse
      double precision rminfo(rmmodes,4) ! data for raman spec
      double precision eigenvecs(nc,nc)  ! [input] (xyz&atom,mode)
      double precision eigenvals(nc)     ! [input] (mode)
      double precision mass(nat)         ! [input] (masses)
      double precision step_size         ! [input] multiplictive factor for step along normal mode
      double precision master(3,nat)     ! saved initial coordinates
      double precision coord(3,nat)      ! displaced coordinate passed to raman_modestep
      character*32 pname
      character*(*) lbl_raman
      character*255 lbl_raman1
c --------- FA vars ---------- START
      character*255 fname_raman
      integer modecc
      logical raman_read_list
      logical oprint_dalphadq
      character*255 filedadq
      external raman_read_list

c ====== Step 0: Read rminfo file(s) ====== START
c Note.- Determine missing modes : modes_missing() and 
c        do calc. to generate them.
c ---- Reading all lbl_ramanXXX files ---- START
         status=raman_read_list(
     &           rmmodes,        ! in : total nr of modes
     &           mode_ini,       ! in : mode_ini
     &           mode_end,       ! in : mode_end
     &           rminfo,         ! in:  raman info updated
     &           nmodes_missing, ! out: nr of modes missing
     &           modes_missing)  ! out: list of missing modes

         if (ga_nodeid().eq.0) then
          write(luout,*) 'List of modes missing: nmodes_missing=',
     &                nmodes_missing
           do i=1,nmodes_missing
            write(luout,112) i,modes_missing(i)
 112        format('mode(',i4,')=',i4)
           enddo ! end-loop-i
         endif
c ======== construct filename in listing mode ====== START
        write(lbl_raman1,'(a5,"_list_",i4.4,"-",i4.4)') 
     &                    lbl_raman,1,nmodes_missing       
        call util_file_name(lbl_raman1,
     &                     .false.,.false.,fname_raman)
        if (ga_nodeid().eq.0)
     &   write(luout,121) fname_raman
 121     format('Generating fname_raman2=',a)         
c ======== construct filename in listing mode ====== END
c ************************************************************
c ************ skipping calculations of missing modes *** START
        if (ignore_missingmode.eq.1) then
          if (ga_nodeid().eq.0) 
     &    write(luout,*) 'SKIP calc modes missing!!'
          goto 123
        endif
        oprint_dalphadq=util_print('raman_dalphadq', print_high).and.
     c       ga_nodeid().eq.0
      if(oprint_dalphadq) then
         call util_file_name('dadq', .false., .false.,filedadq)
         call raman_open_dadq(filedadq)
      endif
c ************ skipping calculations of missing modes *** END
c ************************************************************
        modecc=1 ! reset counter         
        do cmode=1,nmodes_missing 
         ii=modes_missing(cmode)+first0-1

         if (ga_nodeid().eq.0) then
          write(luout,174) ii,rminfo(ii,1)
 174      format('(ii,rminfo(ii,1))=(',
     &            i5,',',f15.8,')')
         endif

         if (ga_nodeid().eq.0) then
           write(luout,113) ii
 113       format("Calculating mode=",i4)
         endif

         call raman_per_mode(
     &             ii,         ! in: ii-th raman mode
     &             modecc,     ! in: mode counter
     &             fname_raman,! in: filename to store raman data
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom,  ! in: (rtdb,geom) handles
     &        oprint_dalphadq ,
     F        filedadq)
        enddo ! end-loop-cmode
 123  continue ! skipping
      if(oprint_dalphadq) then
         call raman_close_dadq(filedadq)
      endif
      return
      end
      subroutine raman_default(
     &             modecc,     ! in: mode counter
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             lbl_raman,  ! in: raman label
     &             begin,      ! in: 
     &             last,       ! in: 
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom)  ! in: (rtdb,geom) handles
c
c Authors: Jonathan Mullin, Northwestern University (ver 1: Jan. 2011)
c          Fredy W. Aquino, Northwestern University (ver 2: Oct. 2012)

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
      logical status
      logical debug
      logical ramanspc ! FA-06-18-12 to store aoresp data in files
      integer rtdb,geom! [input] handles for rtdb
      integer nat      ! [input] # of atoms
      integer nc       ! [input] # of coordinates
      integer i,ii,iii ! countine indexes
      integer begin,last
      integer nmodes_slc,
     &        cmode
      integer rmmodes  ! # of raman active modes
      integer steps    ! [input] # of steps used in Finited differencing
      integer nfreq    ! [input] number of frequencies from property
      integer first0   ! [input] first mode to calculate for aoresponse
      double precision rminfo(rmmodes,4) ! data for raman spec
      double precision eigenvecs(nc,nc)  ! [input] (xyz&atom,mode)
      double precision eigenvals(nc)     ! [input] (mode)
      double precision mass(nat)         ! [input] (masses)
      double precision step_size         ! [input] multiplictive factor for step along normal mode
      double precision master(3,nat)     ! saved initial coordinates
      double precision coord(3,nat)      ! displaced coordinate passed to raman_modestep
      character*32 pname
      character*(*) lbl_raman
      character*255 lbl_raman1
c --------- FA vars ---------- START
      character*255 fname_raman
      integer modecc
      character*255 filedadq
      logical oprint_dalphadq
c ======== construct filename in partitioning mode ====== START
        write(lbl_raman1,'(a5,"_",i4.4,"-",i4.4)') 
     &                    lbl_raman,begin,last       
        call util_file_name(lbl_raman1,
     &                     .false.,.false.,fname_raman)
        if (ga_nodeid().eq.0)
     &   write(luout,111) fname_raman
 111     format('Generating fname_raman3=',a)
         oprint_dalphadq=util_print('raman_dalphadq', print_high).and.
     c        ga_nodeid().eq.0
      if(oprint_dalphadq) then
         call util_file_name('dadq', .false., .false.,filedadq)
         call raman_open_dadq(filedadq)
      endif
c ======== construct filename in partitioning mode ====== END
c ====== Step 0: Read rminfo file(s) ====== END
c        if (begin .gt. last)  go to 200    
        modecc=1 ! reset counter         
        do ii=begin,last ! 2*rmmodes in the end. shorten for tests.

         if (ga_nodeid().eq.0) then
          write(luout,175) ii,rminfo(ii,1)
 175      format('(ii,rminfo(ii,1))=(',
     &            i5,',',f15.8,')')
         endif
         call raman_per_mode(
     &             ii,         ! in: ii-th raman mode
     &             modecc,     ! in: mode counter
     &             fname_raman,! in: filename to store raman data
     &             rmmodes,    ! in: total nr. modes
     &             rminfo,     ! in: stores raman info
     &             nc,nat,     ! in: (nc,nat)=(nr coord,nr atoms)
     &             first0,     ! in: = 6 or 7
     &             eigenvecs,  ! in: hessian data (modes)
     &             eigenvals,  ! in: hessian data (frequencies)
     &             mass,       ! in: mass(i) i=1,nat
     &             step_size,  ! in: step size for displacing geom
     &             steps,      ! in: steps for finite difference
     &             pname,      ! in: label
     &             master,     ! in: master geometry
     &             ramanspc,   ! in: logical = T save in files aoresp data
     &             rtdb,geom,  ! in: (rtdb,geom) handles
     &        oprint_dalphadq ,
     F        filedadq)
        enddo ! ii 
 9      if(oprint_dalphadq) then
         call raman_close_dadq(filedadq)
      endif
      return
      end
