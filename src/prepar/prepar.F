       logical function prepar(irtdb0)
c
c $Id: prepar.F,v 1.44 2001-11-21 21:39:30 d3j191 Exp $
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
      logical pre_mknoe
      external pre_rtdbin,pre_mkseq,pre_mktop,pre_task,pre_mkrst
      external pre_mknoe
c
      integer irtdb0,irtdb,itask
c
      integer mcount
      parameter(mcount=10)
      integer mfract
      parameter(mfract=250)
c
      character*255 filpdb,filseq,filtop,filrst,filpar,filtmp,filcmd
      character*255 filmod,prefix,filxyz,filqqq,filnoe,filpmf,filpov
      character*255 filmrg
      character*255 dir_s,dir_x,dir_q,dir_u,dir_t,dir_c
      character*255 par_s,par_x,par_q,par_u,par_t,par_c
      character*80 sysnam,calc,ffield
      character*10 slvmdl
      character*3 slvnam
      character*80 source
      character*1 altloc,chain,his
      character*4 scount(mcount)
c
      integer lfnpdb,lfnout,lfnfrg,lfnseq,lfnsgm,lfntop,lfnrst,lfnpar
      integer lfnxyz,lfnqqq,lfnnoe,lfnpmf,lfnpov,lfnmat,lfnmrg
      integer len,lenc,lend,lfntmp,lfncmd,lfnslv,lfnmod,nxlnk,mdold
      integer newtop,newseq,newrst,mgrid,mnoe,maxscf,icount(mcount)
      integer ncount,icyren,model,iconst
      integer nfract,ifract(mfract)
      integer ignore
      real*8 gdist,qscale
      logical lstate
      real*8 scale,cpk
c
      irtdb=irtdb0
c
      call util_print_push
      call util_print_rtdb_load(irtdb,'prepar')
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
      lfnmod=22
      lfnxyz=23
      lfnqqq=24
      lfnnoe=25
      lfnpmf=26
      lfnpov=27
      lfnmat=28
      lfnmrg=29
c
c     node 0 only please
c     ------------------
c
      prepar=.true.
      if(ga_nodeid().ne.0) then
  100 continue
      if(.not.pre_task(irtdb,itask,lfnout))
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
      if(util_print('information',print_low)) then
      call util_print_centered(lfnout,'Prepare Module',40,.true.)
      endif
      call ecce_print_module_entry('prepare')
c
c     get info from rtdb
c     ------------------
c
      if(.not.pre_rtdbin(irtdb,ffield,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + par_s,par_x,par_q,par_u,par_t,par_c,
     + source,
     + sysnam,calc,slvnam,slvmdl,newtop,newseq,newrst,mcount,ncount,
     + icount,mgrid,
     + gdist,mnoe,maxscf,qscale,altloc,chain,icyren,iconst,
     + model,nxlnk,mdold,
     + ignore,scount,mfract,nfract,ifract,scale,cpk,his))
     + call errquit('pre_rtdbin failed',9999)
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
      lend=index(dir_c,' ')-1
      filpdb=sysnam(1:len)//'.pdb '
      filpov=sysnam(1:len)//'.pov '
      if(source.ne.'geometry'.and.source.ne.'pdb'.and.
     + source(1:1).ne.' ') filpdb=source
      filxyz=sysnam(1:len)//'.xyz '
      filqqq=sysnam(1:len)//'.q '
      filnoe=sysnam(1:len)//'.drs '
      filpmf=sysnam(1:len)//'.pmf '
      filtop=sysnam(1:len)//'.top '
      filrst=sysnam(1:len)//'_'//calc(1:lenc)//'.rst '
      filseq=dir_c(1:lend)//sysnam(1:len)//'.seq '
      filtmp=sysnam(1:len)//'.nam '
      filmrg=' '
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,prefix))
     + call errquit('rtdb_cget failed on file_prefix',9999)
      filcmd=prefix(1:index(prefix,' ')-1)//'.cmd '
      filmod=prefix(1:index(prefix,' ')-1)//'.mod '
c
      if(par_s(1:1).eq.' '.and.dir_s(1:1).ne.' ') par_s=filpar
      if(par_x(1:1).eq.' '.and.dir_x(1:1).ne.' ') par_x=filpar
      if(par_q(1:1).eq.' '.and.dir_q(1:1).ne.' ') par_q=filpar
      if(par_u(1:1).eq.' '.and.dir_u(1:1).ne.' ') par_u=filpar
      if(par_t(1:1).eq.' '.and.dir_t(1:1).ne.' ') par_t=filpar
      if(par_c(1:1).eq.' '.and.dir_c(1:1).ne.' ') par_c=filpar
c
c     directories
c     -----------
c
      if(util_print('files',print_default)) then
      write(lfnout,2000) ffield(1:index(ffield,' ')-1)
      if(dir_s(1:1).ne.' ')
     + write(lfnout,2001) dir_s(1:index(dir_s,' ')-1),
     + par_s(1:index(par_s,' ')-1)
      if(dir_x(1:1).ne.' ')
     + write(lfnout,2002) dir_x(1:index(dir_x,' ')-1),
     + par_x(1:index(par_x,' ')-1)
      if(dir_q(1:1).ne.' ')
     + write(lfnout,2003) dir_q(1:index(dir_q,' ')-1),
     + par_q(1:index(par_q,' ')-1)
      if(dir_u(1:1).ne.' ')
     + write(lfnout,2004) dir_u(1:index(dir_u,' ')-1),
     + par_u(1:index(par_u,' ')-1)
      if(dir_t(1:1).ne.' ')
     + write(lfnout,2005) dir_t(1:index(dir_t,' ')-1),
     + par_t(1:index(par_t,' ')-1)
      if(dir_c(1:1).ne.' ')
     + write(lfnout,2006) dir_c(1:index(dir_c,' ')-1),
     + par_c(1:index(par_c,' ')-1)
 2000 format(/,' Force field',t40,a,/)
 2001 format(' Directory standard files',t40,a,/,
     + ' File name standard parameters',t40,a,/)
 2002 format(' Directory extensions files',t40,a,/,
     + ' File name extensions parameters',t40,a,/)
 2003 format(' Directory QA test files',t40,a,/,
     + ' File name QA test parameters',t40,a,/)
 2004 format(' Directory user preferences files',t40,a,/,
     + ' File name user preference parameters',t40,a,/)
 2005 format(' Directory temporary files',t40,a,/,
     + ' File name temporary parameters',t40,a,/)
 2006 format(' Directory current files',t40,a,/,
     + ' File name current parameters',t40,a,//)
      endif
c
c     check if the source file exists
c     -------------------------------
c
      if(source(1:8).ne.'geometry') then
      open(unit=lfnpdb,file=filpdb(1:index(filpdb,' ')-1),
     + form='formatted',status='old',err=7)
      source='pdb'
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
c     check if the topology file exists
c     ---------------------------------
c
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=1)
      if(newtop.eq.0) then
      close(unit=lfntop)
      if(util_print('files',print_default)) then
      write(lfnout,2015) filtop(1:index(filtop,' ')-1)
 2015 format(' Using existing topology',t40,a)
      endif
      goto 4
      else
      close(unit=lfntop,status='delete')
      if(util_print('files',print_default)) then
      write(lfnout,2016) filtop(1:index(filtop,' ')-1)
 2016 format(' Deleted existing topology',t40,a)
      endif
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
      if(util_print('files',print_default)) then
      write(lfnout,2007) filseq(1:index(filseq,' ')-1)
 2007 format(' Using existing sequence',t40,a,/)
      endif
      goto 3
      else
      close(unit=lfnseq,status='delete')
      if(util_print('files',print_default)) then
      write(lfnout,2008) filseq(1:index(filseq,' ')-1)
 2008 format(' Deleted existing sequence',t40,a,/)
      endif
      endif
c
    2 continue
c
c     generate sequence & segments from coordinates
c     ---------------------------------------------
c
      if(.not.pre_mkseq(irtdb,lfnout,source,ffield,sysnam,
     + lfnpdb,filpdb,lfnseq,filseq,lfnpar,lfnfrg,lfnsgm,
     + lfntmp,filtmp,lfnmod,filmod,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + par_s,par_x,par_q,par_u,par_t,par_c,
     + slvnam,slvmdl,maxscf,qscale,altloc,chain,icyren,model,nxlnk,
     + mcount,ncount,icount,scount,mfract,nfract,ifract,iconst,his))
     + call errquit('pre_mkseq failed',9999)
      if(util_print('files',print_default)) then
      write(lfnout,2009) filseq(1:index(filseq,' ')-1)
 2009 format(/,' Created sequence',t40,a,/)
      endif
c
    3 continue
c
c     generate topology
c     -----------------
c
      if(.not.pre_mktop(lfnout,ffield,
     + lfnseq,filseq,lfntop,filtop,lfnsgm,lfnpar,lfnmod,filmod,lfnmat,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c,
     + par_s,par_x,par_q,par_u,par_t,par_c,
     + ignore,slvnam))
     + call errquit('pre_mktop failed',9999)
      if(util_print('files',print_default)) then
      write(lfnout,2010) filtop(1:index(filtop,' ')-1)
 2010 format(/,' Created topology',t40,a,/)
      endif
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
      if(util_print('files',print_default)) then
      write(lfnout,2011) filrst(1:index(filrst,' ')-1)
 2011 format(' Using existing restart',t40,a,/)
      endif
      goto 6
      else
      close(unit=lfnrst,status='delete')
      if(util_print('files',print_default)) then
      write(lfnout,2012) filrst(1:index(filrst,' ')-1)
 2012 format(' Deleted existing restart',t40,a,/)
      endif
      endif
c
    5 continue
c
c     generate restart file
c     ---------------------
c
      if(.not.pre_mkrst(irtdb,mdold,source,lfnout,lfntop,filtop,
     + lfnpdb,filpdb,lfnrst,filrst,lfntmp,filtmp,lfncmd,filcmd,
     + lfnslv,slvmdl,slvnam,mgrid,gdist,mnoe,altloc,chain,model,
     + lfnxyz,filxyz,lfnqqq,filqqq,lfnpov,filpov,scale,cpk,lfnmat,
     + lfnmrg,filmrg))
     + call errquit('pre_mkrst failed',9999)
c
    6 continue
c
c     handle other commands
c     ---------------------
c
      if(.not.pre_mknoe(lfnout,lfntop,filtop,lfncmd,filcmd,
     + lfnnoe,filnoe,mnoe,lfnpmf,filpmf,slvnam))
     + call errquit('re_mknoe failed',9999)
c
c     remove temporary files if they exist
c     ------------------------------------
c
      open(unit=lfnmod,file=filmod(1:index(filmod,' ')-1),
     + form='formatted',status='old',err=11)
      close(unit=lfnmod,status='delete')
   11 continue
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=12)
      close(unit=lfncmd,status='delete')
   12 continue
c      open(unit=lfntmp,file=filtmp(1:index(filtmp,' ')-1),
c     + form='formatted',status='old',err=13)
c      close(unit=lfntmp,status='delete')
c   13 continue
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
      if(.not.pre_task(irtdb,itask,lfnout))
     + call errquit('pre_task failed',9999)
c
      lstate=rtdb_parallel(.true.)
c
      prepar=.true.
      return
      end
