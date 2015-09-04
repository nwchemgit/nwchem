c $Id: xc_getv.F 19821 2010-12-14 07:46:49Z d3y133 $
c
      Subroutine dim_getv(rtdb, Edim, ecoul,nExc, iVxc_opt, g_xcinv,
     &                   g_dens, g_vxc, IOLGC, g_wght, g_xyz,g_nq, 
     &                   wght_GA, rho_n, rdens_atom,
     &                   cetobfr, natoms)
c
      implicit none
#include "errquit.fh"
c      
      integer nExc
      integer iVxc_opt
      integer g_xcinv, g_dens(2), g_vxc(4), g_wght, g_xyz,g_nq
      integer natoms
      logical IOLGC, wght_GA
      integer rtdb
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "oep.fh"
#include "dftpara.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"
#include "case.fh"
c
      integer cetobfr(2,natoms)
      double precision rho_n, rdens_atom(ipol*natoms*natoms)
      double precision jfac(4),kfac(4)
      integer g_jk(4), g_d(4)
      logical havehfxc
c     
      integer  ga_create_atom_blocked
      logical xc_gotxc
      external ga_create_atom_blocked,xc_gotxc
c
c--> XC Energy
c
      double precision Edim(2)
      double precision ecoul ! [output]
c
c This driver routine solves for the XC energy and potential (Vxc) via 
c numerical quadrature methods. The results are obtained either by direct 
c numerical integration or by means of a LSQ fit of the Vxc to a set of 
c Gaussian functions. This fitted function can be used to evaluate Vxc
c via a summation of a series of 3-center overlap integrals (3OIs). The
c algorithms are formulated in terms of matrix products. See subsequent
c subroutines for further explanation.
c 
c              XC Energy and Potential Index Key, Vxc(pq,i)
c
c              Value of     |     Definition of index "i"
c            ipol     nExc  |    1        2        3       4
c           --------------------------------------------------
c              1        1   |   Vxc      
c              2        1   |   Vxc^up   Vxc^dw   
c              1        2   |   Vxc      
c              2        2   |   Vxc^up   Vxc^dw   
c
c           nTcols = ipol 
c
      integer me,nTrows,nTcols
      integer lTmat,iTmat,g_truevxc(2)
      double precision zero,one,onem
      logical oprint_intermediate_xc, oprint_time, 
     ,     oprint_oep
      parameter(zero=0.d0,one=1.d0,onem=-1.d0)
      double precision tol2e
      integer g_tmp(2)
c
c     timings
c
      double precision time1_2e,time2_2e
      double precision time1_xc,time2_xc
c
c******************************************************************************
c
c Compute the matrix elements for the XC potential and energy.
c
c      write(luout,*) "dim_getv start"
      oprint_intermediate_xc = util_print('intermediate XC matrix',
     $     print_debug)
      oprint_time = util_print('dft timings', print_high)
      oprint_oep = util_print('oep', print_high)
      Edim(1)=0.d0
      Edim(2)=0.d0
c
      me=ga_nodeid()
c
c     Get the DFT exchange-correlation contribution
      if(util_print('dft timings', print_high)) 
     &      time1_xc=util_cpusec()   ! start xc build time
c
      call dim_grid_quadv0(rtdb, g_vxc, nexc, rho_n,  Edim)
c
      if(util_print('dft timings', print_high)) 
     &         time2_xc=util_cpusec()   ! end xc build time
c
c     print fock_xc build time
      if(util_print('dft timings', print_high)) then
        if (me.eq.0) then
          write(*,"(4x,'Fock_xc Build Time:',F13.1,'s')") 
     &                 time2_xc-time1_xc
        endif
      end if
c
c      write(luout,*) "dim_getv end"
      return
      end
