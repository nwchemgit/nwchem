c
c     == Raman input handler ==

      subroutine raman_input(rtdb)
c
c     Jonathan Mullin, Northwestern University (ver 1: Jan. 2011)
c     Fredy W. Aquino, Northwestern University (ver 2: Oct. 2012)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "rtdb.fh"
#include "mafdecls.fh"

      integer rtdb
      integer rmmodes,first,last

      character*32 rtdb_id
      character*255 test
      character*255 field
      character*16 plot ! normal or resonance
      character*16 line !  lorentzian (l) or gaussian (g) lineshapec

      double precision width ! full-width at half maximum (FWHM) in 1/cm
      double precision step_size ! multiplictive factor for step along normal mode
      double precision frequency
      double precision low, high ! lowest and highest wavenumber normal modes to use.

      logical hyperraman
      logical vroa
c ---------- FA-definitions -------- START
      integer ignore_missingmode
      logical freq_done, ! FA-06-16-12
     &        ramanspc,
     &        ramanlist,
     &        preraman
      integer icount,ii,nmodes_slc,
     &        l_rmodes,k_rmodes
      logical raman_modes_list
      integer cphf3write
      logical status
c ---------- FA-definitions -------- END

      character*32 pname    
      pname = 'raman_input: '
c
c  set some defaults
c
      field=' '
      plot = 'normal' ! normal or resonance
      line = 'lorentzian' !  lorentzian (l) or gaussian (g) lineshape
      width = 20.0D+00 ! full-width at half maximum (FWHM) in 1/cm
      step_size = 0.010D+00
      hyperraman = .false. ! flag to calculate hyperaman terms
      vroa = .false. ! flag to calculate vibrational raman spec
      rmmodes = 0
      last = 10000
      low = 0.0D+00
      high = 100000.0D+00
      freq_done = .false. ! flag to avoid recalc *.nmode if exists FA-06-16-12
      ramanspc  = .false. ! flag to save aoresp data from raman calc.
      ramanlist = .false. ! flag to read rminfo from files and determine 
                          ! missing modes to calculate
      preraman  = .false. ! flag to do task_freq() creating XXX.nmode file
      raman_modes_list=.false.  ! default for input selected raman modes
c
c ----------------------- read values from raman input ------------------------------
c
c     set prop:newaoresp 1  ! set the new solver in the rtdb
      if (.not.rtdb_put(rtdb, 'prop:newaoresp', mt_int, 1,1))
     & call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c
 10   if(.not.inp_read() )
     $  call errquit(pname//'unexpected eof',911, INPUT_ERR)
      if (.not. inp_a(test) )
     $  call errquit(pname//'failed to read field',911, INPUT_ERR)
c        ---- have we ended ? ---
      if(inp_compare(.false.,'end',test)) then
         return
c        --- determine plot style for output ---
      else if(inp_compare(.false.,'normal',test)) then
         plot="normal"
         if (.not. rtdb_cput(rtdb,'raman:plot',1,plot))
     $    call errquit(pname//'rtdb put failed',0, RTDB_ERR)
      else if(inp_compare(.false.,'resonance',test)) then
         plot="resonance"
         if (.not. rtdb_cput(rtdb,'raman:plot',1,plot))
     $    call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c          --- are we calculating hyper raman spec. --- 
      else if(inp_compare(.false.,'hyperraman',field)) then
              hyperraman = .true.
       if (.not.rtdb_put(rtdb,'raman:hyperraman',mt_log,1,hyperraman))
     $    call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c          --- are we calculating VROA spec. --- 
      else if(inp_compare(.false.,'vroa',field)) then
              vroa = .true.
       if (.not.rtdb_put(rtdb,'raman:vroa',mt_log,1,vroa))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c        --- determine lineshape to be used to broaden the peaks ---
      else if(inp_compare(.false.,'lorentzian',test)) then
         if (.not. rtdb_cput(rtdb,'raman:line',1,line))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
      else if(inp_compare(.false.,'gaussian',test)) then
         if (.not. rtdb_cput(rtdb,'raman:line',1,line))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c        --- determine width of the peaks ---
      else if(inp_compare(.false.,'width',test)) then
         if(.not. inp_f(width))
     $     call errquit(pname//'missing value for width',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'raman:width',mt_dbl,1,width))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c        --- determine step size along normal mode ---
      else if(inp_compare(.false.,'dq',test)) then
         if(.not. inp_f(step_size))
     $     call errquit(pname//'missing value for dq',911, INPUT_ERR)
         if (.not.rtdb_put(rtdb,'raman:step_size',mt_dbl,1,step_size))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c        --- determine highest wavernumber normal mode to use ---
      else if(inp_compare(.false.,'high',test)) then
         if(.not. inp_f(high))
     $     call errquit(pname//'missing value for high',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'raman:high',mt_dbl,1,high))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c        --- determine lowest wavernumber normal mode to use ---
      else if(inp_compare(.false.,'low',test)) then
         if(.not. inp_f(low))
     $     call errquit(pname//'missing value for low',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'raman:low',mt_dbl,1,low))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c        --- determine first normal mode to use ---
      else if(inp_compare(.false.,'first',test)) then
         if(.not. inp_i(first))
     $     call errquit(pname//'missing value for first',911, INPUT_ERR)
c        --- not setting default here, it will be set later after
c            frequency calculation has been done so we know if we have
c            a linear molecule or not
         if (.not. rtdb_put(rtdb,'raman:first',mt_int,1,first))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
      else if(inp_compare(.false.,'last',test)) then
         if(.not. inp_i(last)) ! FA-06-16-12 bug-fixed (BEF: first AFT: last)
     $     call errquit(pname//'missing value for last',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'raman:last',mt_int,1,last))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c========== Added by FA-06-16-12 ========== START
      else if(inp_compare(.false.,'freqdone',test)) then
              freq_done = .true.
       if (.not.rtdb_put(rtdb,'raman:freq_done',mt_log,1,freq_done))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
      else if(inp_compare(.false.,'saveaoresp',test)) then
              ramanspc = .true.
       if (.not.rtdb_put(rtdb,'raman:aores0',mt_log,1,ramanspc))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
       cphf3write=2 ! = 2 it will write aoresp files from cphf3 
                    !     but only after it converge to a solution
       if (.not. rtdb_put(rtdb,'cphf:cphf3write',mt_int,1,
     $                         cphf3write))
     $     call errquit('prop_input: rtdb_put failed', 0, RTDB_ERR)
      else if(inp_compare(.false.,'ramanlist',test)) then
       if (.not.rtdb_put(rtdb,'raman:ramantype',mt_int,1,2))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
       if (inp_i(ignore_missingmode)) then
        if (ignore_missingmode.eq.1) then
           if (.not.rtdb_put(rtdb,'raman:ramanskipmode',mt_int,1,
     &                    ignore_missingmode))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
        endif
       endif 
      else if(inp_compare(.false.,'raman_modes_list',test)) then
       if (.not.rtdb_put(rtdb,'raman:ramantype',mt_int,1,3))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c 00000000000000000000000000000000000000000000000000000000000000
         if (inp_i(nmodes_slc)) then ! xxxxxxxx if-nmodes-START
          if (.not. ma_push_get(mt_int,nmodes_slc,'rmodes',
     &                         l_rmodes,k_rmodes))
     &      call errquit('prp_input ma_push_get fail k_pr',0,MA_ERR)
          icount = 0
          do ii = 1, nmodes_slc
           status = inp_i(int_mb(k_rmodes+ii-1))
           if (status) then
            icount = icount+1
           else
            write(6,'(a,a,i4,a)') 
     &           'Raman property input found',
     &           ' fewer modes than expected, only',icount,
     &           ' will be used'
            nmodes_slc = icount
            if (icount.eq.0) then 
             if (.not.ma_pop_stack(l_rmodes)) call errquit(
     &           'prop_input: ma_pop_stack failed k_rmodes',
     &           0,MA_ERR)
            endif
            goto 13
           endif
          enddo ! end-loop-ii
         endif ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx if-nmodes-END
   13    continue      
c --------- storing selected raman modes in rtdb ------- START
          if (nmodes_slc.gt.0) then
           if (.not. rtdb_put(rtdb, 'raman:nmodesslc',mt_int,
     $                        1,nmodes_slc))
     $      call errquit('prop_input-raman: rtdb_put failed',
     $                   555, RTDB_ERR)
           if (.not. rtdb_put(rtdb, 'raman:modeslist',mt_int,
     $                        nmodes_slc,int_mb(k_rmodes)))
     $      call errquit('prop_input-raman: rtdb_put failed',
     $                   555, RTDB_ERR)
          if (.not.ma_pop_stack(l_rmodes)) call
     &        errquit('prop_input: ma_pop_stack failed k_rmodes',
     &                0,MA_ERR)
          endif
c --------- storing selected raman modes in rtdb ------- END
c 00000000000000000000000000000000000000000000000000000000000000
      else if(inp_compare(.false.,'preraman',test)) then
              preraman = .true.
       if (.not.rtdb_put(rtdb,'raman:preraman',mt_log,1,preraman))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c========== Added by FA-06-16-12 ========== END
      else
           call errquit(pname//'unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end
c $Id$
