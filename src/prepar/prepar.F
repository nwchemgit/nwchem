       logical function prepar(irtdb0)
c
c $Id:
c
c ********************************************************
c ********************************************************
c **                                                    **
c **  nwprep is the NWChem molecular dynamics prepare   **
c **  module                                            **
c **                                                    **
c **  Author:                                           **
c **    Dr. T. P. Straatsma                             **
c **    High Performance Computational Chemistry        **
c **    Environmental Molecular Sciences Laboratory     **
c **    Pacific Northwest National Laboratory           **
c **  Copyright 1998                                    **
c **              Pacific Northwest National Laboratory **
c **                                                    **
c **  This module is not parallelized                   **
c **  This module is based in part on the utility codes **
c **  of the molecular simulation package ARGOS,        ** 
c **  copyright 1989,1995 Dr.T.P.Straatsma              **
c **                                                    **
c ********************************************************
c ********************************************************
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
c
      logical pre_rtdbin,pre_mkseq,pre_mktop
      external pre_rtdbin,pre_mkseq,pre_mktop
c
      integer irtdb0,irtdb
c
      character*255 filpdb,filtyp,filseq,filtop,filpar
      character*255 dir_s,dir_x,dir_u,dir_t
      character*80 sysnam,ffield,slvmdl
      character*3 slvnam
c
      integer lfnpdb,lfnout,lfnfrg,lfntyp,lfnseq,lfnsgm,lfntop,lfnpar
      integer len,lend,lprint
      integer newtop,newseq
c
      irtdb=irtdb0
c
      call util_print_push
      call util_print_rtdb_load(irtdb,'prepar')
c
c     node 0 only please
c     ------------------
c
      prepar=.true.
      if(ga_nodeid().ne.0) return
c
c     header
c     ------
c
      write(6,*)
      write(6,*)
      call util_print_centered(6,'nwARGOS Preparation Phase',40,.true.)
      write(6,*)
      write(6,*)
      call ecce_print_module_entry('prepare')
c
c     get info from rtdb
c     ------------------
c
      if(.not.pre_rtdbin(irtdb,ffield,dir_s,dir_x,dir_u,dir_t,sysnam,
     + slvnam,slvmdl,newtop,newseq))
     + call errquit('pre_rtdbin failed',9999)
c
c     construct file names
c     --------------------
c
      len=index(ffield,' ')-1
      lend=index(dir_s,' ')-1
      filtyp=dir_s(1:lend)//ffield(1:len)//'.typ '
      filpar=ffield(1:len)//'.par '
c
      len=index(sysnam,' ')-1
      lend=index(dir_t,' ')-1
      filpdb=sysnam(1:len)//'.pdb '
      filtop=sysnam(1:len)//'.top '
      filseq=dir_t(1:lend)//sysnam(1:len)//'.seq '
c
c     set logical file numbers to be used
c     -----------------------------------
c
      lfnout=6
      lfnpdb=12
      lfnfrg=13
      lfntyp=14
      lfnseq=15
      lfnsgm=16
      lfntop=17
      lfnpar=18
c
      lprint=10
      if(util_print('sequence',print_none)) lprint=0
      if(util_print('sequence',print_low)) lprint=1
      if(util_print('sequence',print_medium)) lprint=2
      if(util_print('sequence',print_high)) lprint=4
      if(util_print('sequence',print_debug)) lprint=5
c
c     check if the topology file exists
c     ---------------------------------
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=1)
      if(newtop.eq.0) then
      close(unit=lfntop)
      write(lfnout,1000) filtop(1:index(filtop,' ')-1)
 1000 format('USING EXISTING TOPOLOGY FILE ',a)
      goto 4
      else
      close(unit=lfntop,status='delete')
      write(lfnout,1005) filtop(1:index(filtop,' ')-1)
 1005 format('DELETED EXISTING TOPOLOGY FILE ',a)
      endif
c
    1 continue
c
c     check if the sequence file exists
c     ---------------------------------
c
      open(unit=lfnseq,file=filseq(1:index(filseq,' ')-1),
     + form='formatted',status='old',err=2)
      if(newseq.eq.0) then
      close(unit=lfnseq)
      write(lfnout,1001) filseq(1:index(filseq,' ')-1)
 1001 format('USING EXISTING SEQUENCE FILE ',a)
      goto 3
      else
      close(unit=lfnseq,status='delete')
      write(lfnout,1006) filseq(1:index(filseq,' ')-1)
 1006 format('DELETED EXISTING SEQUENCE FILE ',a)
      endif
c
    2 continue
c
c     generate sequence & segments from coordinates
c     ---------------------------------------------
c
      if(.not.pre_mkseq(lfnout,lprint,ffield,
     + lfnpdb,filpdb,lfnseq,filseq,lfntyp,filtyp,
     + lfnfrg,lfnsgm,dir_s,dir_x,dir_u,dir_t,slvnam))
     + call errquit('pre_mkseq failed',9999)
      write(lfnout,1002) filseq(1:index(filseq,' ')-1)
 1002 format('CREATED SEQUENCE FILE ',a)
c
    3 continue
c
c     generate topology
c     -----------------
c
      if(.not.pre_mktop(lfnout,lprint,ffield,
     + lfnseq,filseq,lfntop,filtop,lfnsgm,lfnpar,filpar,
     + dir_s,dir_x,dir_u,dir_t))
     + call errquit('pre_mktop failed',9999)
      write(lfnout,1003) filtop(1:index(filtop,' ')-1)
 1003 format('CREATED TOPOLOGY FILE ',a)
c
    4 continue
c
c
      call util_print_pop
c
      prepar=.true.
      return
      end
