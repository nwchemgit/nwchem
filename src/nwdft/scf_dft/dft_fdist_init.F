      subroutine dft_fdist_init()
      implicit none
#include "dist.fh"
#include "dft_fdist.fh"
c
c     xc_tabcd
      monitor_tabcd = getdist( 'xc_tabcd')
      monitor_comm_ga_acc1 = getdist( 'comm tabcd ga_acc1')
      monitor_comm_ga_acc2 = getdist( 'comm tabcd ga_acc2')
      monitor_size_ga_acc1 = getdist( 'size tabcd ga_acc1')
      monitor_screen0 = getdist( 'comp tabcd screen0')
      monitor_screen1 = getdist( 'comp tabcd screen1')
      monitor_screen2 = getdist( 'comp tabcd screen2')
      monitor_screen3 = getdist( 'comp tabcd screen3')
      monitor_screen4 = getdist( 'comp tabcd screen4')
      monitor_screen5 = getdist( 'comp tabcd screen5')
      monitor_comp_transp = getdist( 'comp tabcd transp')
      monitor_mult1 = getdist( 'comp tabcd mult1')
      monitor_mult2 = getdist( 'comp tabcd mult2')
      monitor_wait1 = getdist( 'comm tabcd wait1')
      monitor_wait2 = getdist( 'comm tabcd wait2')
c     xcrhogen
      monitor_xcrho = getdist( 'xc_rhogen')
      monitor_gaget = getdist( 'comm rhogen ga_get')
      monitor_ganbget = getdist( 'comm rhogen ga_nbget')
      monitor_wait3 = getdist( 'comm rhogen wait')
      monitor_size_ga_get = getdist( 'size rhogen ga_get')
      monitor_rscreen0 = getdist( 'comp rhogen screen0')
      monitor_rhocomp = getdist( 'comp rhogen mult1')
      monitor_rhocomp2 = getdist( 'comp rhogen mult2')
c
      end
c $Id$
