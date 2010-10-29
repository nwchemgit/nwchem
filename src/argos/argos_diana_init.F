      subroutine argos_diana_init
c
c $Id$
c
c
c     initialize variables
c     --------------------
c
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer i
      character*255 filcmd
c
c     set sequential access to rtdb
c     and ignore current state
c
      ignore=rtdb_parallel(.false.)
c
      lfnram=60
      lfndat=61
      lfnexp=62
      lfntcl=63
      lfnpca=64
      lfnavg=65
      lfnhol=66
      lfnpin=67
      lfnpou=68
      lfncmd=69
      lfnloc=70
      lfnhba=71
      lfnecc=72
      lfnrdf=73
      lfnref=76
      lfntrj=77
      lfncop=78
      lfnsup=79
      lfnchg=80
      lfnplt=81
      lfnxyz=82
      lfnana=83
      lfnrms=84
      lfnprj=85
      lfnval=86
      lfncov=87
      lfnvec=88
      lfnmin=89
      lfnmax=90
      lfntop=91
      lfnhis=92
      lfnepz=93
      lfngrp=94
      lfnord=95
      lfnsel=96
      lfndie=97
      lfnamp=98
c
      ntval=0
      ntvec=0
c
      ltop=.false.
c
      do 44 i=1,mhis
      idhis(i,3)=40+i
   44 continue
c      do 46 i=1,mord
c      rord(i,1)=zero
c      rord(i,2)=zero
c   46 continue
c
      lsonly=.false.
      lsuper=.false.
c
      ifrfr=1
      ifrto=0
      ifrsk=1
      ifrst=0
      ilast=0
c
      ndist=0
      nangl=0
      ntors=0
      nimps=0
      nord=0
c
      ndef=0
c
      nhis=0
      lhis=0
c
      nrot=0
c
      lsel=.false.
      lsels=.false.
      lselo=.false.
      active=.false.
      lrms=.false.
      lana=.false.
      lrama=.false.
      lesd=.false.
      lloc=.false.
      lhbond=.false.
      lrdf=.false.
      lesppt=.false.
c
      nave=0
      ndata=0
      time=0.0d0
      timr=0.0d0
      ldist=0
      icent=0
      jcent=0
      scale=1.0d0
      cpk=1.0d0
      stick=1.0d0
c
      ngrid=10
      iesppb=0
      rcut=1.0d0
c
      nfrtot=0
      nfrme=0
c
      indx=1
c
      rhbmin=0.0d0
      rhbmax=0.2d0
      hbdmin=1.9d0
      hbdmax=3.1416d0
      hbamin=1.9d0
      hbamax=3.1416d0
      numwhb=0
c
      ngroup=0
      ngroups=0
      if(.not.ma_push_get(mt_int,maxgrp*mxdef*mxnum,'iwrk',l_wrk,i_wrk))
     + call md_abort('Could not allocate wrk',0)
c
c
c     open the command file
c     ---------------------
c
      if(me.eq.0) then
      if(.not.rtdb_cget(irtdb,'ana:filcmd',1,filcmd))
     + call md_abort('Error rtdb_get filcmd',0)
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=9999)
      rewind(unit=lfncmd)
      endif
c
      return
 9999 continue
      call md_abort('Could not open command file',0)
      return
      end
