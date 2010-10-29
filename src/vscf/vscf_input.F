      subroutine vscf_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
c#include "coptopt.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb, ngrid, ncoup 
      integer iexc, ilen
      character*255 field
      double precision vcfct
      logical dmdr
c
c  set some defaults
c
      ngrid = 16      ! number of grid points per dimension
      ncoup = 2       ! 1 for diagonal terms, 2 for pair couplings
      vcfct = 1.0d+0  ! possible scale factor for pair-coupling potential
      iexc = 1        ! level of excitation; 1 - fundamentals
      dmdr = .false.  ! don't just use the dipole derivative to get intensities
 
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('vscf_input: unexpected eof',911, INPUT_ERR)
      if (.not. inp_a(field))
     $     call errquit('vscf_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'ngrid',field)) then
         if(.not. inp_i(ngrid))
     $        call errquit('vscf_input: missing grid size',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'vscf:ngrid',mt_int,1,ngrid))
     $        call errquit('vscf_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'iexcite',field)) then
         if(.not. inp_i(iexc))
     $        call errquit('vscf_input: missing grid size',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'vscf:iexc',mt_int,1,iexc))
     $        call errquit('vscf_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'coupling',field)) then
         if(.not. inp_a(field))
     $        call errquit('vscf_input: missing coupling type',911,
     &       INPUT_ERR)
         if(inp_compare(.false.,'diagonal',field).or.
     $      inp_compare(.false.,'diag',field)) then
              ncoup = 1
         else if(inp_compare(.false.,'pair',field)) then
              ncoup = 2
         else if(inp_compare(.false.,'triplet',field)) then
              ncoup = 3
         else
              call errquit('vscf_input: unrecognized coupling type',0,
     &       INPUT_ERR)
         endif
         if (.not. rtdb_put(rtdb,'vscf:ncoup',mt_int,1,ncoup))
     $        call errquit('vscf_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'scale_factor',field)) then
         if(.not. inp_f(vcfct))
     $        call errquit('vscf_input: missing scale factor',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'vscf:vcfct',mt_dbl,1,vcfct))
     $        call errquit('vscf_input: rtdb put failed',0, RTDB_ERR)
c
      else
         call errquit('vscf_input: unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end
