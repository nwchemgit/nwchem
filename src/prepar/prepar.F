       logical function prepar(irtdb0)
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
#include "inp.fh"
c
      logical pre_rtdbin,pre_mkseq,pre_mktop
      external pre_rtdbin,pre_mkseq,pre_mktop
c
      integer irtdb0,irtdb
c
      character*80 filpdb,ffield,filtyp,filseq,filtop,filpar
      character*80 dir_s,dir_x,dir_u,dir_t,sysnam
      character*3 namslv
c
      integer lfnpdb,lfnout,lfnfrg,lfntyp,lfnseq,lfnsgm,lfntop,lfnpar
      integer len,lend
c
      ffield='amber '
      filpdb='/usr/people/d3j191/crown.pdb '
      filseq='/usr/people/d3j191/test/crown.seq '
      filtop='/usr/people/d3j191/test/crown.top '
      filtyp='/usr/people/d3j191/nwdata/amber.typ '
      filpar='amber.par'
      dir_s='/usr/people/d3j191/data/amber_s/ '
      dir_x='/usr/people/d3j191/data/amber_x/ '
      dir_u='/usr/people/d3j191/data/amber_u/ '
      dir_t='/usr/people/d3j191/test/ '
      namslv='HOH'
      lfnout=6
      lfnpdb=12
      lfnfrg=13
      lfntyp=14
      lfnseq=15
      lfnsgm=16
      lfntop=17
      lfnpar=18
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
      irtdb=irtdb0
      if(.not.pre_rtdbin(irtdb,ffield,dir_s,dir_x,dir_u,dir_t,sysnam))
     + call errquit('pre_rtdbin failed',9999)
c
c     construct file names
c     --------------------
c
      len=index(sysnam,' ')-1
      lend=index(dir_t,' ')-1
      filpdb=sysnam(1:len)//'.pdb '
      filtop=sysnam(1:len)//'.top '
      filseq=dir_t(1:lend)//sysnam(1:len)//'.seq '
c
c     check if the topology file exists
c     ---------------------------------
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=1)
      close(unit=lfntop)
      write(lfnout,1000) filtop(1:index(filtop,' ')-1)
 1000 format('TOPOLOGY FILE ',a,' EXISTS')
      goto 4
c
    1 continue
c
c     check if the sequence file exists
c     ---------------------------------
c
      open(unit=lfnseq,file=filseq(1:index(filseq,' ')-1),
     + form='formatted',status='old',err=2)
      close(unit=lfnseq)
      write(lfnout,1001) filseq(1:index(filseq,' ')-1)
 1001 format('SEQUENCE FILE ',a,' EXISTS')
      goto 3
c
    2 continue
c
c     generate sequence & segments from coordinates
c     ---------------------------------------------
c
      if(.not.pre_mkseq(lfnout,ffield,
     + lfnpdb,filpdb,lfnseq,filseq,lfntyp,filtyp,
     + lfnfrg,lfnsgm,dir_s,dir_x,dir_u,dir_t,namslv))
     + call errquit('pre_mkseq failed',9999)
      write(lfnout,1002) filseq(1:index(filseq,' ')-1)
 1002 format('SEQUENCE FILE ',a,' CREATED')
c
    3 continue
c
c     generate topology
c     -----------------
c
      if(.not.pre_mktop(lfnout,.true.,ffield,
     + lfnseq,filseq,lfntop,filtop,lfnsgm,lfnpar,filpar,
     + dir_s,dir_x,dir_u,dir_t))
     + call errquit('pre_mktop failed',9999)
      write(lfnout,1003) filtop(1:index(filtop,' ')-1)
 1003 format('TOPOLOGY FILE ',a,' CREATED')
c
    4 continue
c
c     generate topology
c     -----------------
c
      prepar=.true.
      return
      end
