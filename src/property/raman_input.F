c
c     == Raman input handler ==
c
      subroutine raman_input(rtdb)
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer rmmodes,first,last
c
      character*32 rtdb_id
      character*255 test
      character*255 field
      character*16 plot ! normal or resonance
      character*16 line !  lorentzian (l) or gaussian (g) lineshapec
c
      double precision width ! full-width at half maximum (FWHM) in 1/cm
      double precision step_size ! multiplictive factor for step along normal mode
      double precision frequency
      double precision low, high ! lowest and highest wavenumber normal modes to use.
c
      logical hyperraman
      logical vroa
      logical status
c
      character*32 pname    
      pname = 'raman_input: '
c
c  set some defaults
c
      plot = 'normal' ! normal or resonance
      line = 'lorentzian' !  lorentzian (l) or gaussian (g) lineshape
      width = 20.0D+00 ! full-width at half maximum (FWHM) in 1/cm
      step_size = 0.010D+00
      hyperraman = .false. ! flag to calculate hyperaman terms
      vroa = .false. ! flag to calculate vibrational raman spec
      rmmodes = 0
      first = 7
      last = 10000
      low = 0.0D+00
      high = 100000.0D+00
c
c ----------------------- read values from raman input ------------------------------
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
         if (.not. rtdb_put(rtdb,'raman:first',mt_int,1,first))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
c        --- determine first normal mode to use ---
      else if(inp_compare(.false.,'last',test)) then
         if(.not. inp_i(first))
     $     call errquit(pname//'missing value for last',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'raman:last',mt_int,1,last))
     $     call errquit(pname//'rtdb put failed',0, RTDB_ERR)
      else
           call errquit(pname//'unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end
