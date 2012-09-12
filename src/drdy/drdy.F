C> \ingroup task
C> @{
      logical function task_drdy(rtdb)
      implicit none
*
      logical drdy_nwchem
      external drdy_nwchem
*
      integer rtdb
*
      task_drdy = drdy_nwchem(rtdb)
      end
C> @}
      logical function drdy_nwchem(rtdb)
*
* $Id$
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
C     6  - standard output 
C    30  - RPH data
C
      Implicit none             ! Double Precision (A-H,O-Z)
#include "errquit.fh"
C
#include "stdio.fh"
#include "drdyP.fh"
#include "nwc_drdyP.fh"
#include "global.fh"
#include "rtdb.fh"
C
*::passed:
      Integer rtdb ! [input] RTDB handle
*::local:
      Logical status_rtdb
      Logical lhead,lsad,lskip
      integer jop
      integer ns
      integer jopmx
      integer i, ii, j
      integer i_dummy
*
      me_drdy = ga_nodeid()
      my_rtdb = rtdb
      fu6 = luout
      fu30 = 30
*      call setdbg(1)
      
      status_rtdb = rtdb_parallel(.false.)
      if (me_drdy.ne.0) then
        i_dummy = 0
        call drdy_synch(i_dummy,'drdy')
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
            write (fu6,*) ' trouble with restart, file unreadable'
            call errquit('drdy: fatal error',911, UNKNOWN_ERR)
         endif
      endif
      jopmx = 3
      if (lgs(6).eq.1.or.lgs(6).eq.3) jopmx = 4
C
C  Reactant and product calculations
      if (jop.lt.jopmx) Call drdy_react(jop)
C
      if (lgs(1).ne.0) then
C  Saddle point calculations
         if (lsad) then
C  Saddle point information successful read in from fort.30
            do i = 1,n3
               xspsv(i) = xr(i,5)
            enddo
         else
C  Calculate potential and hessian at saddle point geometry
            Call drdy_saddle
            Call drdy_rphwrt(5)
            Call drdy_rphwrt(6)
         endif
      else
         if (ns.le.0) then
C  No saddle point, if ns<=0 compute v, dx, h at starting geometry
            call drdy_saddle
C  lsad is set true in sort if the two lines beginning the MEP data
C   are already written
            if(.not.lsad) call drdy_rphwrt(6)
            call drdy_rphwrt(7)
            ns = 1
         else
            ii = 0
            vspsv = vgrid(1) + vzero
            vspspsv = vgrdsp(1) + vzerosp
            do i = 1,n3
               xspsv(i) = xgrid(i,1)
               dxspsv(i) = dxgrid(i,1)
               do j = 1,i
                  ii = ii + 1
                  fspsv(j,i) = hgrid(ii,1)
                  fspsv(i,j) = hgrid(ii,1)
               enddo
            enddo
         endif
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
        if(.not.lskip) Call drdy_setpth
C
C  Calculate steepest descent path
C
        Call drdy_path(ns)
C
C  Sort data along reaction path and rewrite file30
C
        Call drdy_sort(lhead,jop,lsad,ns,1)
      endif
      write (fu6,*) 'drdy_NWChem has finished correctly, bye !'
      call drdy_synch(DRDY_CODE_DONE,'drdy 2')
99999 continue
*      call setdbg(0)
      call ga_sync()
      status_rtdb = rtdb_parallel(status_rtdb)
      drdy_nwchem = .true.
      End
      
