       logical function prepar(irtdb0)
c
c $Id: prepar.F,v 1.16 1998-04-29 18:35:31 d3j191 Exp $
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
#include "msgids.fh"
c
      logical pre_rtdbin,pre_mkseq,pre_mktop,pre_task,pre_mkrst
      external pre_rtdbin,pre_mkseq,pre_mktop,pre_task,pre_mkrst
c
      integer irtdb0,irtdb,itask
c
      character*255 filpdb,filseq,filtop,filrst,filpar,filtmp,filcmd
      character*255 dir_s,dir_x,dir_u,dir_t
      character*80 sysnam,calc,ffield
      character*10 slvmdl
      character*3 slvnam
      character*80 source
c
      integer lfnpdb,lfnout,lfnfrg,lfnseq,lfnsgm,lfntop,lfnrst,lfnpar
      integer len,lenc,lend,lprint,lfntmp,lfncmd,lfnslv
      integer newtop,newseq,newrst,mgrid,mnoe,maxscf
      real*8 gdist,qscale
      logical lstate
c
      irtdb=irtdb0
c
      call util_print_push
      call util_print_rtdb_load(irtdb,'prepar')
      lprint=10
      if(util_print('prepare',print_none)) lprint=0
      if(util_print('prepare',print_low)) lprint=1
      if(util_print('prepare',print_medium)) lprint=2
      if(util_print('prepare',print_high)) lprint=4
      if(util_print('prepare',print_debug)) lprint=5
c
      lstate=rtdb_parallel(.false.)
c
c     set logical file numbers to be used
c     -----------------------------------
c
      lfnout=6
      lfnpdb=12
      lfnfrg=13
      lfnseq=14
      lfnsgm=15
      lfntop=16
      lfnrst=17
      lfnpar=18
      lfntmp=19
      lfncmd=20
      lfnslv=21
c
c     node 0 only please
c     ------------------
c
      prepar=.true.
      if(ga_nodeid().ne.0) then
  100 continue
      if(.not.pre_task(irtdb,itask))
     + call errquit('pre_task failed',9999)
      if(itask.eq.0) then
      lstate=rtdb_parallel(.true.)
      return
      endif
      goto 100
      endif
c
c     header
c     ------
c
      if(lprint.ge.3) then
      call util_print_centered(lfnout,'Prepare Module',40,.true.)
      endif
      call ecce_print_module_entry('prepare')
c
c     get info from rtdb
c     ------------------
c
      if(.not.pre_rtdbin(irtdb,ffield,dir_s,dir_x,dir_u,dir_t,source,
     + sysnam,calc,slvnam,slvmdl,newtop,newseq,newrst,mgrid,gdist,mnoe))
     + call errquit('pre_rtdbin failed',9999)
c
c     directories
c     -----------
c
      if(lprint.ge.4) then
      write(lfnout,'(a,a)') 'ffield=',ffield(1:index(ffield,' ')-1)
      write(lfnout,'(a,a)') 'dir_s =',dir_s(1:index(dir_s,' ')-1)
      write(lfnout,'(a,a)') 'dir_x =',dir_x(1:index(dir_x,' ')-1)
      write(lfnout,'(a,a)') 'dir_u =',dir_u(1:index(dir_u,' ')-1)
      write(lfnout,'(a,a)') 'dir_t =',dir_t(1:index(dir_t,' ')-1)
      endif
c
c     construct file names
c     --------------------
c
      len=index(ffield,' ')-1
      lend=index(dir_s,' ')-1
      filpar=ffield(1:len)//'.par '
c
      len=index(sysnam,' ')-1
      lenc=index(calc,' ')-1
      lend=index(dir_t,' ')-1
      filpdb=sysnam(1:len)//'.pdb '
      filtop=sysnam(1:len)//'.top '
      filrst=sysnam(1:len)//'_'//calc(1:lenc)//'.rst '
      filseq=dir_t(1:lend)//sysnam(1:len)//'.seq '
      filtmp=sysnam(1:len)//'_tmp_pdb.nam '
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,filcmd))
     + call errquit('rtdb_cget failed on file_prefix',9999)
      filcmd=filcmd(1:index(filcmd,' ')-1)//'.cmd '
c
c     check if the source file exists
c     -------------------------------
c
      if(source(1:1).eq.' '.or.source(1:3).eq.'pdb') then
      open(unit=lfnpdb,file=filpdb(1:index(filpdb,' ')-1),
     + form='formatted',status='old',err=7)
      if(source(1:1).eq.' ') source='pdb'
      close(unit=lfnpdb)
      goto 8
    7 continue
      if(source(1:3).eq.'pdb')
     + call errquit('pdb file could not be found',9999)
    8 continue
      endif
      if(source(1:1).eq.' ') source='geometry'
      if(source(1:4).eq.'rtdb') source='geometry'
c
c
c     check if the topology file exists
c     ---------------------------------
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=1)
      if(newtop.eq.0) then
      close(unit=lfntop)
      write(lfnout,1001) filtop(1:index(filtop,' ')-1)
 1001 format('Using existing topology file ',a)
      goto 4
      else
      close(unit=lfntop,status='delete')
      write(lfnout,1002) filtop(1:index(filtop,' ')-1)
 1002 format('Deleted existing topology file ',a)
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
      write(lfnout,1003) filseq(1:index(filseq,' ')-1)
 1003 format('Using existing sequence file ',a)
      goto 3
      else
      close(unit=lfnseq,status='delete')
      write(lfnout,1004) filseq(1:index(filseq,' ')-1)
 1004 format('Deleted existing sequence file ',a)
      endif
c
    2 continue
c
c     generate sequence & segments from coordinates
c     ---------------------------------------------
c
      if(.not.pre_mkseq(irtdb,lfnout,lprint,source,ffield,
     + lfnpdb,filpdb,lfnseq,filseq,lfnpar,filpar,lfnfrg,lfnsgm,
     + lfntmp,filtmp,dir_s,dir_x,dir_u,dir_t,slvnam,slvmdl,maxscf,
     + qscale))
     + call errquit('pre_mkseq failed',9999)
      write(lfnout,1005) filseq(1:index(filseq,' ')-1)
 1005 format('Created sequence file ',a)
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
      write(lfnout,1006) filtop(1:index(filtop,' ')-1)
 1006 format('Created topology file ',a)
      newrst=1
c
    4 continue
c
c     check if the command file exists
c     --------------------------------
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=10)
      close(unit=lfncmd)
      goto 5
   10 continue
c
c     check if the restart file exists
c     --------------------------------
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + form='formatted',status='old',err=5)
      if(newrst.eq.0) then
      close(unit=lfnrst)
      write(lfnout,1007) filrst(1:index(filrst,' ')-1)
 1007 format('Using existing restart file ',a)
      goto 6
      else
      close(unit=lfnrst,status='delete')
      write(lfnout,1008) filrst(1:index(filrst,' ')-1)
 1008 format('Deleted existing restart file ',a)
      endif
c
    5 continue
c
c     generate restart file
c     ---------------------
c
      if(.not.pre_mkrst(irtdb,source,lfnout,lprint,lfntop,filtop,
     + lfnpdb,filpdb,lfnrst,filrst,lfntmp,filtmp,lfncmd,filcmd,
     + lfnslv,slvmdl,slvnam,mgrid,gdist,mnoe))
     + call errquit('pre_mkrst failed',9999)
c
    6 continue
c
c     remove temporary file if they exist
c     -----------------------------------
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=11)
      close(unit=lfncmd,status='delete')
   11 continue
      open(unit=lfntmp,file=filtmp(1:index(filtmp,' ')-1),
     + form='formatted',status='old',err=12)
      close(unit=lfntmp,status='delete')
   12 continue
c
c     reset print
c     -----------
c
      call util_print_pop
c
c     signal other nodes to return
c     ----------------------------
c
      itask=0
      if(.not.pre_task(irtdb,itask))
     + call errquit('pre_task failed',9999)
c
      lstate=rtdb_parallel(.true.)
c
      prepar=.true.
      return
      end
