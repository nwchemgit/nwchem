      logical function drdy_nwchem(rtdb)
*
* $Id: drdy.F,v 1.1 1999-08-16 22:54:33 d3e129 Exp $
*
C
C DRDYNWChem - DiRect DYnamics with the NWChem electronic structure code
C   This code prepares the file30 input for POLYRATE from electronic 
C   structure calculations of gradients
C   and hessians at the reactant, product, and saddle point geometries, 
C   and along the minimum energy path.
C   Cartesian geometries for the reactants, products, and saddle points 
C   need to be input to this code;
C   optimization for these geometries is not performed in the code. 
C   Points along the minimum energy path
C   are calculated here.
C
C Read coordinate X, gradient DX (negative of the force), and the
C   Force constant matrix F from the runtime database and hessian file of 
C   NWChem predict the next step along the IRC, generate input file 
C   (via the RTDB) for NWChem initiate NWChem ab initio run
C
C  File useage
*?*C     5  - standard input
*?*C     6  - standard output 
*?*C     9  - new Gaussian input file
*?*C    10  - formatted checkpoint file
C    30  - RPH data
C
      Implicit none             ! Double Precision (A-H,O-Z)
C
#include "stdio.fh"
#include "drdyP.fh"
#include "nwc_drdyP.fh"
#include "global.fh"
C
*::passed:
      Integer rtdb ! [input] RTDB handle
*::local:
      Logical lhead,lsad,lskip
      integer jop
      integer ns
      integer jopmx
      integer i
      integer i_dummy
*
      me_drdy = ga_nodeid()
      my_rtdb = rtdb
      if (me_drdy.ne.0) then
        i_dummy = 0
        call drdy_synch(i_dummy)
        goto 99999
      endif
C
C  Read in input data
      Call Input_for_drdy(rtdb)
C
C  Set up for constructing nwchem info on rtdb
      call potset_nwchem(rtdb)
C
      jop = 0
      lsad = .false.
      ns = 0
c
      if (lgs(8).eq.0) then
C
         call util_file_unlink('fort.30')
C  Initiate file 30
         Call drdy_rphwrt(0)
      else
C
C  Read restart data from file 30
         Call drdy_sort(lhead,jop,lsad,ns,0)
         if (.not.lhead) then
            write (luout,*) ' trouble with restart, file unreadable'
            call errquit('drdy: fatal error',911)
         endif
      endif
      jopmx = 3
      if (lgs(6).eq.1.or.lgs(6).eq.3) jopmx = 4
C
C  Reactant and product calculations
      if (jop.lt.jopmx) then
         Call drdy_react(jop)
      endif
C
C  Saddle point calculations
      if (lsad) then
         do i = 1,n3
            xspsv(i) = xr(i,5)
         enddo
      else
         Call drdy_saddle
         Call drdy_rphwrt(6)
      endif
      if (slp.ne.0.0d0.or.slm.ne.0.0d0) then
C
C  Check whether to skip the saddle point set up
C
        if (lgs(1).eq.0) then
C
C    no saddle point so skip
C
           lskip = .true.
        else if (ns.gt.0) then
C
C restart calculation; last s values on reactant and product sides should be in 
C sgrid(1) and sgrid(ns) - 
C skip if either the reactant or product side has not be done
C
           lskip = sgrid(1).lt.0.0d00 .and. sgrid(ns).gt.0.0d00
        else
C for saddle point with no previous calculation, don't skip
C
           lskip = .false.
        endif
C
C  Set up to take first step off saddle point
C
        if(.not.lskip) then
           Call drdy_setpth
        endif
C
C  Calculate steepest descent path
C
        Call drdy_path(ns)
C
C  Sort data along reaction path and rewrite file30
C
        Call drdy_sort(lhead,jop,lsad,ns,1)
      endif
      write (luout,*) 'drdy_NWChem has finished corretly, bye !'
99999 continue
      call ga_sync()
      drdy_nwchem = .true.
      End
      
