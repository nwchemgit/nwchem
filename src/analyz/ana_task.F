      subroutine ana_task
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      logical ignore
c
      integer lfncmd
      character*255 filcmd,fil,string
c
      character*80 card
      integer nsai,msai,nwmi,mwmi,nwai,mwai,it,ndx
      integer nframe,nread,num,i,j,k,ivec,ia,ja
      real*8 rms0,rms1,val,range
      character*3 fmt
      logical active
c
      logical ana_rdfram,ana_wrtvec
      real*8 ana_bond,ana_angle,ana_torsion
      external ana_rdfram,ana_wrtvec
      external ana_bond,ana_angle,ana_torsion
      integer idpdb,i_temp,l_temp
      logical lsuper,lsaver
      character*8 option
c
c     set sequential access to rtdb
c     and ignore current state
c
      ignore=rtdb_parallel(.false.)
c
      lfntcl=64
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
      ltop=.false.
c
      do 44 i=1,mhis
      idhis(i,3)=40+i
   44 continue
      do 46 i=1,mord
      rord(i,1)=zero
      rord(i,2)=zero
   46 continue
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
      nbnds=0
      nangs=0
      ntors=0
      nimps=0
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
      ifr=1
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
      if(me.eq.0) then
      if(.not.rtdb_cget(irtdb,'ana:filcmd',1,filcmd))
     + call md_abort('Error rtdb_get filcmd',0)
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=9999)
      rewind(unit=lfncmd)
      endif
c
    1 continue
      lesd=.false.
c
      if(me.eq.0) then
      read(lfncmd,1000,end=9999) card
 1000 format(a)
      endif
      call util_char_ga_brdcst(mag_d00,card,0)
c
c     read reference coordinates
c     --------------------------
c
      if(card(1:6).eq.'refer ') then
      filref=card(8:80)
c
c     get the size from the reference file
c
      call ana_sizref(nsai,msai,nwmi,mwmi,nwai,mwai)
c
c     allocate memory and initialize
c
      call ana_init(nsai,msai,nwmi,mwmi,nwai,mwai,.true.)
c
c     read the reference coordinates
c
      call ana_rdref(dbl_mb(i_xref))
c
c
c     determine which atoms will be handled by this node
c
      num=nsa/np+1
      ifr=me*num+1
      ito=(me+1)*num
      if(ifr.gt.nsa) then
      num=0
      ifr=nsa
      else
      if(ito.gt.nsa) ito=nsa
      num=ito-ifr+1
      endif
c
      call ana_all(int_mb(i_ssel),1,dbl_mb(i_swt),nsel)
      lsel=.false.
      call ana_all(int_mb(i_isel),1,dbl_mb(i_wt),nsel)
      call ana_wall(int_mb(i_wsel),1,int_mb(i_ndxw))
c
      if(me.eq.0) then
      write(*,'(a,a)') ' Reference coordinates read from ',
     + filref(1:index(filref,' ')-1)
      write(*,'(/,a,i5)') ' Number of atoms is ',nsai
      endif
      goto 1
      endif
c
c     read topology file
c     ------------------
c
      if(card(1:6).eq.'topol ') then
      filtop=card(8:80)
      call ana_siztop()
      if(.not.ma_push_get(mt_int,7*msgm,'iram',l_iram,i_iram))
     + call md_abort('Could not allocate iram',0)
      if(.not.ma_push_get(mt_int,msa,'imol',l_imol,i_imol))
     + call md_abort('Could not allocate imol',0)
      if(.not.ma_push_get(mt_int,2*msb,'sbnd',l_sbnd,i_sbnd))
     + call md_abort('Could not allocate sbnd',0)
      call ana_rdtop(byte_mb(i_snam),dbl_mb(i_qdat),dbl_mb(i_qwdat),
     + int_mb(i_iram),int_mb(i_imol),int_mb(i_sbnd))
c
      if(me.eq.0) then
      if(ltop) then
      write(*,'(a,a)') ' Topology read from ',
     + filtop(1:index(filtop,' ')-1)
      else
      write(*,'(a,a)') ' Topology could not be read from ',
     + filtop(1:index(filtop,' ')-1)
      endif
      endif
      goto 1
      endif
c
c     read trajectory file header
c     ---------------------------
c
      if(card(1:6).eq.'file  ') then
      read(card(8:17),'(2i5)') ifrst,ilast
      filtrj=card(18:80)
c
      fmttrj='trj'
      if(index(filtrj,'.sco').gt.0) fmttrj='sco'
      if(index(filtrj,'.coo').gt.0) fmttrj='coo'
c
      call ana_rdhdr(byte_mb(i_snam))
c
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
      call ana_all(int_mb(i_isel),1,dbl_mb(i_wt),nsel)
      call ana_wall(int_mb(i_wsel),1,int_mb(i_ndxw))
      lsel=.false.
c
      if(me.eq.0) then
      write(*,'(/,a,a)') ' Trajectory file header from ',
     + filtrj(1:index(filtrj,' ')-1)
      endif
      goto 1
      endif
c
      if(card(1:6).eq.'scale ') then
      read(card(7:18),'(f12.6)') scale
      goto 1
      endif
c
      if(card(1:6).eq.'cpk   ') then
      read(card(7:18),'(f12.6)') cpk
      goto 1
      endif
c
      if(card(1:6).eq.'stick ') then
      read(card(7:18),'(f12.6)') stick
      goto 1
      endif
c
c     atom selection
c     --------------
c
      if(card(1:6).eq.'select') then
      if(.not.lsel) then
      call ana_all(int_mb(i_isel),0,dbl_mb(i_wt),nsel)
      lsel=.true.
      endif
      call ana_select(card,byte_mb(i_snam),int_mb(i_isel),dbl_mb(i_wt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,3003) nsel,nsa
 3003 format(' Selected',i5,' out of',i5,' atoms')
      endif
      goto 1
      endif
      if(card(1:6).eq.'selecs') then
      if(.not.lsels) then
      call ana_all(int_mb(i_ssel),0,dbl_mb(i_swt),nsel)
      lsels=.true.
      endif
      call ana_select(card,byte_mb(i_snam),int_mb(i_ssel),dbl_mb(i_swt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,3533) nsel,nsa
 3533 format(' Selected',i5,' out of',i5,' atoms for superposition')
      endif
      goto 1
      endif
      if(card(1:6).eq.'seleco') then
      if(.not.lselo) then
      call ana_all(int_mb(i_osel),0,dbl_mb(i_owt),nselo)
      lselo=.true.
      endif
      call ana_select(card,byte_mb(i_snam),int_mb(i_osel),dbl_mb(i_owt),
     + dbl_mb(i_xref),nselo)
      if(me.eq.0) then
      write(*,3633) nselo,nsa
 3633 format(' Selected',i5,' out of',i5,' atoms for output')
      endif
      goto 1
      endif
c
c     save atom selection
c     -------------------
c
      if(card(1:6).eq.'savsel') then
      filsel=card(8:80)
      call ana_savsel(lfnsel,filsel,int_mb(i_wsel),int_mb(i_isel))
      goto 1
      endif
      if(card(1:6).eq.'savsls') then
      filsel=card(8:80)
      call ana_savsel(lfnsel,filsel,int_mb(i_wsel),int_mb(i_ssel))
      goto 1
      endif
c
c     read atom selection
c     -------------------
c
      if(card(1:5).eq.'rdsel') then
      filsel=card(8:80)
      call ana_readsel(lfnsel,filsel,int_mb(i_wsel),int_mb(i_isel))
      call ana_wcount(int_mb(i_wsel),int_mb(i_ndxw))
      goto 1
      endif
      if(card(1:6).eq.'rdsels') then
      filsel=card(8:80)
      call ana_readsel(lfnsel,filsel,int_mb(i_wsel),int_mb(i_ssel))
      call ana_wcount(int_mb(i_wsel),int_mb(i_ndxw))
      goto 1
      endif
c
c     index
c     -----
c
      if(card(1:6).eq.'index ') then
      read(card(7:16),'(i10)') indx
      goto 1
      endif
c
c     set hydrogen bond criteria
c     --------------------------
c
      if(card(1:7).eq.'hbdist ') then
      read(card(8:31),'(2f12.6)') rhbmin,rhbmax
      goto 1
      endif
      if(card(1:7).eq.'hbdang ') then
      read(card(8:31),'(2f12.6)') hbdmin,hbdmax
      goto 1
      endif
      if(card(1:7).eq.'hbaang ') then
      read(card(8:31),'(2f12.6)') hbamin,hbamax
      goto 1
      endif
      if(card(1:7).eq.'numwhb ') then
      read(card(8:12),'(i5)') numwhb
      goto 1
      endif
      if(card(1:6).eq.'hbond  ') then
      lhbond=.true.
      call ana_hbond_init(byte_mb(i_snam),int_mb(i_sbnd),int_mb(i_hbnd))
      goto 1
      endif
c
c     solvent selection
c     -----------------
c
      if(card(1:6).eq.'wselct') then
      read(card(8:19),'(f12.6)') range
c
      call ana_wall(int_mb(i_wsel),0,int_mb(i_ndxw))
c
      if(nsel.eq.0) call md_abort('No solute atoms selected',0)
c
      call ana_rdhdr(byte_mb(i_snam))
      nread=0
      nwrit=0
      nframe=0
   20 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 20
      if(nframe.eq.(nframe/ifrsk)*ifrsk.and.lsx.and.lxw) then
c
      call ana_wselect(int_mb(i_isel),dbl_mb(i_xdat),
     + int_mb(i_wsel),dbl_mb(i_wdat),range,int_mb(i_ndxw))
c
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 20
      call ga_igop(mag_d13,int_mb(i_wsel),mwm,'+')
      endif
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
c
      call ana_wcount(int_mb(i_wsel),int_mb(i_ndxw))
      call ana_all(int_mb(i_isel),1,dbl_mb(i_wt),nsel)
      lsel=.false.
      write(*,3303) nwsel,nwm
 3303 format(/,' Selected',i5,' out of',i5,' solvent molecules')
      endif
c
      goto 1
      endif
c
c     solute atom group definition
c     ----------------------------
c
      if(card(1:6).eq.'define') then
      call ana_define(card,byte_mb(i_snam),int_mb(i_isel),dbl_mb(i_wt),
     + dbl_mb(i_xref))
      if(me.eq.0) then
      write(*,3004) ndef
 3004 format(' Defined',i5,' atom groups')
      endif
      goto 1
      endif
c
c     solvent atom group definition
c     -----------------------------
c
      if(card(1:6).eq.'wdefin') then
      call ana_wdefin(card)
      if(me.eq.0) then
      write(*,3004) ndef
      endif
      goto 1
      endif
c
c     atom group analysis
c     -------------------
c
      if(card(1:6).eq.'group ') then
      call ana_group(card,byte_mb(i_snam),int_mb(i_imol),
     + int_mb(i_isel),dbl_mb(i_wt),dbl_mb(i_xref),int_mb(i_wrk))
      goto 1
      endif
c
c     atom groups directives
c     ----------------------
c
      if(card(1:6).eq.'groups') then
      call ana_groups(card,byte_mb(i_snam),int_mb(i_imol),
     + int_mb(i_isel),dbl_mb(i_wt),dbl_mb(i_xref))
      goto 1
      endif
c
c     atom groups directives
c     ----------------------
c
      if(card(1:6).eq.'histo ') then
      if(card(22:26).eq.'    1') then
      nhis=nhis+1
      rhis=-0.6d0*box(3)
      if(nhis.gt.mhis)
     + call md_abort('Too many histogram directives ',nhis)
      read(card(8:21),'(2i7)') idhis(nhis,1),idhis(nhis,2)
      if(idhis(nhis,2).eq.0) idhis(nhis,2)=mlhis
      do 45 i=1,mlhis
      ihis(i,nhis)=0
   45 continue
      fil=card(54:80)
      if(fil(1:1).eq.' ')
     + call md_abort('No histogram filename given',nhis)
      open(unit=idhis(nhis,3),file=fil(1:index(fil,' ')-1),
     + form='formatted',status='unknown')
      endif
      goto 1
      endif
c
c     radial distribution function
c     ----------------------------
c
      if(card(1:6).eq.'rdf   ') then
      read(card,'(7x,i7,f12.6)') nrdf,rrdf
      lrdf=.true.
      goto 1
      endif
c
c     electrostatic potentials at selected points
c     -------------------------------------------
c
      if(card(1:6).eq.'esppnt') then
      string=card(8:80)
      ndx=index(string,' ')-1
      filpin=string(1:ndx)
      string=string(ndx+2:80)
      ndx=index(string,' ')-1
      filhol=string(1:index(string,' ')-1)
      string=string(ndx+2:80)
      ndx=index(string,' ')-1
      filpou=string(1:index(string,' ')-1)
      string=string(ndx+2:80)
      filavg=string(1:index(string,' ')-1)
      write(*,'(a,a,a,a)') ' esppnt ',
     + filpin(1:index(filpin,' ')),filhol(1:index(filhol,' ')),
     + filpou(1:index(filpou,' ')-1),filavg(1:index(filavg,' ')-1)
      call ana_esppt(int_mb(i_isel),dbl_mb(i_xdat),dbl_mb(i_wdat),
     + dbl_mb(i_qdat),dbl_mb(i_qwdat))
      endif
c
c     order parameters
c     ----------------
c
      if(card(1:7).eq.'order  ') then
      read(card(8:26),'(i5,2i7)') it,ia,ja
      if(it.eq.nord+1) nord=nord+1
      if(it.gt.nord) call md_abort('Order index too large',it)
      if(it.gt.mord) call md_abort('Order index too large',it)
      idord(it,1)=ia
      idord(it,2)=ja
      if(me.eq.0) then
      do 56 i=3,4
      read(lfncmd,1000,end=9999) card
      call ana_index(card,byte_mb(i_snam),idord(it,i),ia)
   56 continue
      endif
      call ga_brdcst(mag_d08,idord,4*mord*ma_sizeof(mt_int,1,mt_byte),0)
      goto 1
      endif
c
c     root mean square deviation
c     --------------------------
c
      if(card(1:6).eq.'rmsdev') then
      lrms=.true.
      ndata=0
      goto 1
      endif
c
c     ramachandran analysis
c     ---------------------
c
      if(card(1:6).eq.'diel  ') then
      ldiel=.true.
      goto 1
      endif
c
c     ramachandran analysis
c     ---------------------
c
      if(card(1:6).eq.'rama  ') then
      lrama=.true.
      goto 1
      endif
c
c     bond specification
c     ------------------
c
      if(card(1:6).eq.'bond  ') then
      read(card(8:12),'(i5)') it
      if(it.eq.nbnds+1) nbnds=nbnds+1
      if(it.gt.nbnds) call md_abort('Bond index too large',it)
      if(it.gt.mxbnds) call md_abort('Bond index too large',it)
      if(me.eq.0) then
      do 7 i=1,2
      read(lfncmd,1000,end=9999) card
      call ana_index(card,byte_mb(i_snam),ibnds(it,i,1),ibnds(it,i,2))
    7 continue
c      write(*,'(2i5)') (ibnds(it,i),i=1,2)
      endif
      call ga_brdcst(mag_d07,ibnds,
     + 4*mxbnds*ma_sizeof(mt_int,1,mt_byte),0)
      goto 1
      endif
c
c     angle specification
c     -------------------
c
      if(card(1:5).eq.'angle') then
      read(card(8:12),'(i5)') it
      if(it.eq.nangs+1) nangs=nangs+1
      if(it.gt.nangs) call md_abort('Angle index too large',it)
      if(it.gt.mxangs) call md_abort('Angle index too large',it)
      if(me.eq.0) then
      do 8 i=1,3
      read(lfncmd,1000,end=9999) card
      call ana_index(card,byte_mb(i_snam),iangs(it,i,1),iangs(it,i,2))
    8 continue
c      write(*,'(3i5)') (iangs(it,i),i=1,3)
      endif
      call ga_brdcst(mag_d08,iangs,
     + 6*mxangs*ma_sizeof(mt_int,1,mt_byte),0)
      goto 1
      endif
c
c     torsion specification
c     ---------------------
c
      if(card(1:6).eq.'torsio') then
      read(card(8:12),'(i5)') it
      if(it.eq.ntors+1) ntors=ntors+1
      if(it.gt.ntors) call md_abort('Torsion index skipped',it)
      if(it.gt.mxtors) call md_abort('Torsion index too large',it)
      if(me.eq.0) then
      do 9 i=1,4
      read(lfncmd,1000,end=9999) card
      call ana_index(card,byte_mb(i_snam),itors(it,i,1),itors(it,i,2))
    9 continue
c      write(*,'(4i5)') (itors(it,i),i=1,4)
      endif
      call ga_brdcst(mag_d09,itors,
     + 8*mxtors*ma_sizeof(mt_int,1,mt_byte),0)
      goto 1
      endif
c
c     improper specification
c     ----------------------
c
      if(card(1:6).eq.'improp') then
      read(card(8:12),'(i5)') it
      if(it.eq.nimps+1) nimps=nimps+1
      if(it.gt.nimps) call md_abort('Improper index too large',it)
      if(it.gt.mximps) call md_abort('Improper index too large',it)
      if(me.eq.0) then
      do 12 i=1,4
      read(lfncmd,1000,end=9999) card
      call ana_index(card,byte_mb(i_snam),iimps(it,i,1),iimps(it,i,2))
   12 continue
c      write(*,'(4i5)') (iimps(it,i),i=1,4)
      endif
      call ga_brdcst(mag_d10,iimps,
     + 8*mximps*ma_sizeof(mt_int,1,mt_byte),0)
      goto 1
      endif
c
c     charge distribution
c     -------------------
c
      if(card(1:6).eq.'chgdst') then
      read(card(8:12),'(i5)') ldist
      goto 1
      endif
c
c     center fraction
c     ---------------
      if(card(1:6).eq.'center') then
      read(card(8:17),'(2i5)') icent,jcent
      if(icent.gt.0) then
c
c     center-fold the reference coordinates
c
      call ana_center_ref(dbl_mb(i_xref))
c
      endif
      goto 1
      endif
c
c     rotations
c
      if(card(1:6).eq.'rotat0') then
      nrot=0
      goto 1
      endif
      if(card(1:6).eq.'rotatx') then
      nrot=nrot+1
      if(nrot.gt.100) call md_abort('Increase dimension irot',0)
      irot(nrot)=1
      read(card(8:19),'(f12.6)') arot(nrot)
      goto 1
      endif
      if(card(1:6).eq.'rotaty') then
      nrot=nrot+1
      if(nrot.gt.100) call md_abort('Increase dimension irot',0)
      irot(nrot)=2
      read(card(8:19),'(f12.6)') arot(nrot)
      goto 1
      endif
      if(card(1:6).eq.'rotatz') then
      nrot=nrot+1
      if(nrot.gt.100) call md_abort('Increase dimension irot',0)
      irot(nrot)=3
      read(card(8:19),'(f12.6)') arot(nrot)
      goto 1
      endif
c
c     atom tag
c     --------
c
      if(card(1:6).eq.'label ') then
      call ana_tag(card,byte_mb(i_snam),dbl_mb(i_xref))
      goto 1
      endif
c
c     select which frames to read
c     ---------------------------
c
      if(card(1:6).eq.'frame ') then
      read(card(8:37),'(3i10)') ifrfr,ifrto,ifrsk
      if(me.eq.0) then
      write(*,'(/,3(a,i10))') ' Selected frames ',ifrfr,' to ',ifrto,
     + ' by',ifrsk
      endif
      timoff=0.0d0
      time=0.0d0
      timr=0.0d0
      goto 1
      endif
c
c     write pdb
c     ---------
c
      if(card(1:6).eq.'wrtpdb') then
      lsuper=card(15:17).eq.'sup'
      lsonly=card(19:21).eq.'sol'
      read(card(7:14),'(i8)') idpdb
      read(card(23:80),'(a)') filcop
      i=index(filcop,'.')
      fmt=filcop(i+1:i+3)
c
      call ana_rdhdr(byte_mb(i_snam))
c
      nread=0
      nwrit=0
      nframe=0
      if(idpdb.gt.0) then
  222 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.idpdb) goto 222
      if(lsuper) then
      call super2(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),
     + num,msa,dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms),
     + nsa,msa,.true.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa)
      call ana_xcollect(dbl_mb(i_xdat),ifr,num,msa)
      endif
      if(fmt.eq.'pdb') then
      call ana_wrtpdb(lfncop,filcop,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xdat),int_mb(i_wsel),dbl_mb(i_wdat),.false.,.not.lsonly,
     + int_mb(i_ndxw))
      elseif(fmt.eq.'xyz') then
      call ana_wrtxyz(lfncop,filcop,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xdat),dbl_mb(i_qdat),int_mb(i_wsel),dbl_mb(i_wdat),
     + .false.,.not.lsonly,int_mb(i_ndxw))
      else
      call md_abort('Illegal file format specified',0)
      endif
      endif
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      close(unit=lfncop)
      write(*,'(/,a,i10,a,a)') ' Written frame ',idpdb,' to ',
     + filcop(1:index(filcop,' ')-1)
      lsonly=.false.
      endif
      else
      if(fmt.eq.'pdb') then
      call ana_wrtpdb(lfncop,filcop,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xref),int_mb(i_wsel),dbl_mb(i_wdat),.false.,.false.,
     + int_mb(i_ndxw))
      elseif(fmt.eq.'xyz') then
      call ana_wrtxyz(lfncop,filcop,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xdat),dbl_mb(i_qdat),int_mb(i_wsel),dbl_mb(i_wdat),
     + .false.,.false.,int_mb(i_ndxw))
      else
      call md_abort('Illegal file format specified',0)
      endif
      if(me.eq.0) then
      close(unit=lfncop)
      write(*,'(/,a,a)') ' Written reference solute coordinates to ',
     + filcop(1:index(filcop,' ')-1)
      lsonly=.false.
      endif
      endif
c
      goto 1
      endif
c
c     copy frames
c     -----------
c
      if(card(1:6).eq.'copy  '.or.card(1:6).eq.'copys ') then
      lsonly=card(1:6).eq.'copys '
      read(card(8:12),'(i5)') mcopf
      read(card(13:24),'(f12.6)') rangle
      read(card(25:80),'(a)') filcop
      i=index(filcop,'.')
      fmt=filcop(i+1:i+3)
      root=filcop(1:i)
      ncop=1
      icopf=1
c
c      if(nsel.eq.0) call md_abort('No atoms selected',0)
c
c     write header
c
      call ana_rdhdr(byte_mb(i_snam))
      call ana_wthdr(lfncop,fmt,byte_mb(i_snam),byte_mb(i_tag),
     + int_mb(i_isel),.not.lsonly)
c
      nread=0
      nwrit=0
      nframe=0
    2 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 2
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      call ana_wtfram(lfncop,fmt,dbl_mb(i_xdat),byte_mb(i_snam),
     + dbl_mb(i_wdat),int_mb(i_isel),byte_mb(i_tag),dbl_mb(i_val),
     + int_mb(i_wsel),int_mb(i_ndxw))
      nwrit=nwrit+1
      indx=indx+1
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 2
      endif
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      close(unit=lfncop)
      write(*,'(/,a,i10,a,a)') ' Copied ',nwrit,' frames to ',
     + filcop(1:index(filcop,' ')-1)
      lsonly=.false.
      endif
c
      goto 1
      endif
c
c     clean frames
c     ------------
c
      if(card(1:6).eq.'clean ') then
      read(card(8:13),'(i6)') nclean
      lclean=nclean.gt.0
      goto 1
      endif
c
c     super frames
c     -----------
c
      if(card(1:6).eq.'super '.or.card(1:6).eq.'supers') then
      lsonly=card(1:6).eq.'supers'
      read(card(8:19),'(f12.6)') rangle
      read(card(20:80),'(a)') filsup
      i=index(filsup,'.')
      fmt=filsup(i+1:i+3)
      root=filsup(1:i)
c
      if(nsel.eq.0) call md_abort('No atoms selected',0)
c
c     write header
c
      call ana_rdhdr(byte_mb(i_snam))
      call ana_wthdr(lfnsup,fmt,byte_mb(i_snam),byte_mb(i_tag),
     + int_mb(i_isel),.not.lsonly)
c
      nread=0
      nwrit=0
      nframe=0
    3 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 3
      if(nframe.eq.(nframe/ifrsk)*ifrsk.and.lsx) then
      call super2(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms),
     + nsa,msa,.true.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa)
      call ana_xcollect(dbl_mb(i_xdat),ifr,num,msa)
      call ana_wtfram(lfnsup,fmt,dbl_mb(i_xdat),byte_mb(i_snam),
     + dbl_mb(i_wdat),int_mb(i_isel),byte_mb(i_tag),dbl_mb(i_val),
     + int_mb(i_wsel),int_mb(i_ndxw))
      nwrit=nwrit+1
      indx=indx+1
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 3
      endif
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      close(unit=lfnsup)
      write(*,'(/,a,i10,a,a)') ' Superimposed ',nwrit,' frames to ',
     + filsup(1:index(filsup,' ')-1)
      lsonly=.false.
      endif
c
      goto 1
      endif
c
c     scan frames
c     -----------
c
      if(card(1:6).eq.'scan  '.or.card(1:6).eq.'scans ') then
      lsuper=card(1:6).eq.'scans '
      read(card(8:80),'(a)') string
      filram=string(1:index(string,' ')-1)//'.ram'
      filrms=string(1:index(string,' ')-1)//'.rms'
      filana=string(1:index(string,' ')-1)//'.ana'
      filhis=string(1:index(string,' ')-1)//'.his'
      filgrp=string(1:index(string,' ')-1)//'.grp'
      filord=string(1:index(string,' ')-1)//'.ord'
      filloc=string(1:index(string,' ')-1)//'.loc'
      filhba=string(1:index(string,' ')-1)//'.hba'
      fildie=string(1:index(string,' ')-1)//'.die'
      time=0.0d0
c
      if(lrdf) then
      call ana_rdf_init()
      print*,'rdf_init'
      endif
c
c     read file header
c
      call ana_rdhdr(byte_mb(i_snam))
c
      lana=nbnds+nangs+ntors+nimps.gt.0
c
      if(me.eq.0) then
      if(lana) then
      open(unit=lfnana,file=filana(1:index(filana,' ')-1),
     + form='formatted',status='unknown',err=9999)
      write(*,3333) filana(1:index(filana,' ')-1)
 3333 format(/,' Opening scn file ',a)
      call ana_anahdr(byte_mb(i_snam))
      endif
      if(lrama) then
      open(unit=lfnram,file=filram(1:index(filram,' ')-1),
     + form='formatted',status='unknown')
      call ana_ramhdr(byte_mb(i_snam),int_mb(i_iram))
      endif
      if(ldiel) then
      open(unit=lfndie,file=fildie(1:index(fildie,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(lrms) then
      open(unit=lfnrms,file=filrms(1:index(filrms,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(lhbond) then
      open(unit=lfnhba,file=filhba(1:index(filhba,' ')-1),
     + form='formatted',status='unknown')
      call ana_hbond_header(byte_mb(i_snam),int_mb(i_isel),
     + int_mb(i_hbnd))
      if(.not.ma_push_get(mt_int,mwm,'temp',l_temp,i_temp))
     + call md_abort('Could not allocate temp',0)
      endif
      if(ngroups.gt.0) then
      open(unit=lfngrp,file=filgrp(1:index(filgrp,' ')-1),
     + form='formatted',status='unknown')
      write(lfngrp,4001) ngroups
 4001 format(i5)
      do 4002 i=1,ngroups
      write(lfngrp,4003) igroups(i,5),
     + ldef(igroups(i,1)),ldef(igroups(i,2))
 4003 format(3i5)
 4002 continue
      endif
      if(nord.gt.0) then
      open(unit=lfnord,file=filord(1:index(filord,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(ngroup.gt.0) then
      open(unit=lfnloc,file=filloc(1:index(filloc,' ')-1),
     + form='formatted',status='unknown')
      call ana_lochdr(byte_mb(i_snam),int_mb(i_imol),int_mb(i_isel))
      endif
      endif
c
      if(ldist.gt.0) then
      if(.not.ma_push_get(mt_dbl,2*ldist,'hist',l_hist,i_hist))
     + call md_abort('Failed to allocate hist',0)
      open(unit=lfnhis,file=filhis(1:index(filhis,' ')-1),
     + form='formatted',status='unknown')
      endif
c
      nread=0
      nwrit=0
      nframe=0
   10 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 10
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
c
      if(lsuper) then
      call super2(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),
     + num,msa,dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms),
     + nsa,msa,.true.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.true.)
      call ana_xcollect(dbl_mb(i_xdat),ifr,num,msa)
      endif
c
c     evaluate bond lenghts
c
      do 11 i=1,nbnds
      rbnds(i)=ana_bond(dbl_mb(i_xdat),msa,dbl_mb(i_wdat),mwm,mwa,
     + ibnds(i,1,1),ibnds(i,1,2),ibnds(i,2,1),ibnds(i,2,2),box)
   11 continue
c
c     evaluate angles
c
      do 13 i=1,nangs
      rangs(i)=ana_angle(dbl_mb(i_xdat),msa,dbl_mb(i_wdat),mwm,mwa,
     + iangs(i,1,1),iangs(i,1,2),iangs(i,2,1),iangs(i,2,2),
     + iangs(i,3,1),iangs(i,3,2))
   13 continue
c
c     evaluate torsions
c
      do 14 i=1,ntors
      rtors(i)=ana_torsion(dbl_mb(i_xdat),msa,dbl_mb(i_wdat),mwm,mwa,
     + itors(i,1,1),itors(i,1,2),itors(i,2,1),itors(i,2,2),
     + itors(i,3,1),itors(i,3,2),itors(i,4,1),itors(i,4,2))
   14 continue
c
c     charge histogram
c
      if(ldist.gt.0) then
      if(.not.ltop) call md_abort('Topology needed',0)
      call ana_chgdst(dbl_mb(i_hist),dbl_mb(i_xdat),dbl_mb(i_qdat))
      endif
c
c     write bonded interaction values to file
c
      if(me.eq.0.and.lana) then
      write(lfnana,2323) time,(rbnds(i),i=1,nbnds),(rangs(i),i=1,nangs),
     + (rtors(i),i=1,ntors)
 2323 format(5f12.6)
      nwrit=nwrit+1
      endif
c
c     ramachandran analysis
c
      if(lrama) then
      call ana_rama(dbl_mb(i_xdat),int_mb(i_iram))
      endif
c
c     dielectric analysis
c
      if(ldiel) then
      call ana_diel(dbl_mb(i_xdat),dbl_mb(i_qdat),int_mb(i_isel),
     + dbl_mb(i_wdat),dbl_mb(i_qwdat),int_mb(i_wsel))
      endif
c
c     root mean square deviation
c
      if(lrms) then
      if(np.gt.1) then
      call ana_rmsdev(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,
     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms))
      else
      call ana_rmsdev(dbl_mb(i_xdat),int_mb(i_idat),nsa,
     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms))
      endif
      endif
c
c     hydrogen bond analysis
c
      if(lhbond) then
      call ana_hbond(dbl_mb(i_xdat),dbl_mb(i_wdat),int_mb(i_isel),
     + int_mb(i_hbnd),int_mb(i_temp))
      endif
c
      if(ngroup.gt.0) then
      do 3307 i=1,ngroup
      if(igroup(i,2).eq.1)
     + call ana_grploc(byte_mb(i_snam),int_mb(i_imol),int_mb(i_isel),
     + dbl_mb(i_wt),dbl_mb(i_xdat),i,int_mb(i_wrk))
 3307 continue
      endif
c
      if(lrdf) then
      call ana_rdf(int_mb(i_isel),dbl_mb(i_xdat),dbl_mb(i_wdat),
     + int_mb(i_rdf))
      endif
c
      if(ngroups.gt.0) then
      if(me.eq.0) write(lfngrp,4004) time
 4004 format('frame',/,f12.6)
      do 3301 i=1,ngroups
      if(igroups(i,5).eq.1)
     + call ana_grpcogdis(byte_mb(i_snam),int_mb(i_imol),int_mb(i_isel),
     + dbl_mb(i_wt),dbl_mb(i_xdat),i)
      if(igroups(i,5).eq.2)
     + call ana_grpdis(byte_mb(i_snam),int_mb(i_imol),int_mb(i_isel),
     + dbl_mb(i_wt),dbl_mb(i_xdat),i)
      if(igroups(i,5).eq.4)
     + call ana_grpcogang(byte_mb(i_snam),int_mb(i_imol),int_mb(i_isel),
     + dbl_mb(i_wt),dbl_mb(i_xdat),i)
      if(igroups(i,5).eq.5)
     + call ana_grpvectors(byte_mb(i_snam),int_mb(i_imol),
     + int_mb(i_isel),dbl_mb(i_wt),dbl_mb(i_xdat),i)
 3301 continue
      endif
c
      if(nhis.gt.0) then
      do 3302 i=1,nhis
      call ana_histo(dbl_mb(i_xdat),dbl_mb(i_wdat),i)
 3302 continue
      endif
c
      if(nord.gt.0) then
      call ana_order(dbl_mb(i_xdat))
      endif
c
      ndata=ndata+1
c
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 10
      endif
c
      if(ldist.gt.0) then
      if(.not.ma_pop_stack(l_hist))
     + call md_abort('Failed to deallocate hist',0)
      close(unit=lfnhis)
      endif
c
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      if(lrdf) call ana_rdfwrt(int_mb(i_isel),int_mb(i_rdf))
      if(lana) then
      close(unit=lfnana)
      write(*,'(a)') ' Closing scn file '
      write(*,'(/,a,i10,a,a)') ' Scanned ',nwrit,' frames to ',
     + filana(1:index(filana,' ')-1)
      endif
      if(lrama) close(unit=lfnram)
      if(ldiel) close(unit=lfndie)
      if(lrms) then
      write(*,'(a)') ' Closing scn file '
      write(*,'(/,a,i10,a,a)') ' Scanned ',ndata,' frames to ',
     + filrms(1:index(filrms,' ')-1)
      endif
      if(ngroups.gt.0) then
      close(unit=lfngrp)
      write(*,'(a)') ' Closing grp file '
      write(*,'(/,a,i10,a,a)') ' Scanned ',ndata,' frames to ',
     + filgrp(1:index(filgrp,' ')-1)
      endif
      if(nord.gt.0) then
      do 3308 i=1,nord
      write(lfnord,3309) i,abs(half*(three*rord(i,2)/dble(ndata)-one))
 3309 format(i5,f12.6)
 3308 continue
      close(unit=lfnord)
      write(*,'(a)') ' Closing ord file '
      write(*,'(/,a,i10,a,a)') ' Scanned ',ndata,' frames to ',
     + filord(1:index(filord,' ')-1)
      endif
      if(nhis.gt.0) then
      do 3304 i=1,nhis
      dhis=1.2d0*box(3)/real(idhis(i,2))
      do 3305 j=1,idhis(i,2)
      write(idhis(i,3),3306) rhis+(j-1)*dhis,real(ihis(j,i))/real(ndata)
 3306 format(f12.6,f12.6)
 3305 continue
      close(unit=idhis(i,3))
 3304 continue
      endif
      endif
c
      if(lrdf) then
      print*,'deallocating rdf'
      if(.not.ma_pop_stack(l_rdf))
     + call md_abort('Failed to deallocate rdf',0)
      endif
c
      if(lhbond) then
      if(.not.ma_pop_stack(l_temp))
     + call md_abort('Failed to deallocate temp',0)
      endif
c
      nbnds=0
      nangs=0
      ntors=0
      nimps=0
      lrms=.false.
      lana=.false.
      ldist=0
c
      goto 1
      endif
c
c     calculate average coordinates
c     -----------------------------
c
      if(card(1:6).eq.'saver '.or.card(1:6).eq.'aver  ') then
c
      lsaver=card(1:6).eq.'saver '
c
      if(card(7:7).eq.' ') then
      filcop=' '
      else
      read(card(7:80),'(a)') filcop
      endif
      i=index(filcop,'.')
      if(i.le.0) then
      fmt='pdb'
      else
      fmt=filcop(i+1:i+3)
      endif
c
      if(nsel.eq.0) call md_abort('No atoms selected',0)
c
      call ana_rdhdr(byte_mb(i_snam))
      nave=0
      nread=0
      nframe=0
    4 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 4
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      if(lsaver) then
      call super2(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms),
     + nsa,msa,.true.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.true.)
      call ana_xcollect(dbl_mb(i_xdat),ifr,num,msa)
      endif
c
      call ana_xaver(dbl_mb(i_xdat+ifr-1),dbl_mb(i_xadat+ifr-1),
     + num,msa,nave)
c
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 4
      endif
c
      if(me.eq.0) then
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
      call ana_xmean(dbl_mb(i_xadat+ifr-1),num,msa,nave)
      if(np.gt.1) call ana_update(dbl_mb(i_xadat),msa,ifr,ifr+num-1)
c
      if(me.eq.0) then
      write(*,'(/,a,i10,a)') ' Averaged ',nave,' frames'
      endif
c
      if(filcop(1:1).ne.' ') then
      if(fmt.eq.'pdb') then
      call ana_wrtpdb(lfncop,filcop,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xadat),int_mb(i_wsel),dbl_mb(i_wdat),.false.,.false.,
     + int_mb(i_ndxw))
      elseif(fmt.eq.'xyz') then
      call ana_wrtxyz(lfncop,filcop,byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_xadat),dbl_mb(i_qdat),int_mb(i_wsel),dbl_mb(i_wdat),
     + .false.,.false.,int_mb(i_ndxw))
      endif
      if(lsaver) then
      write(*,'(/,a,a)')
     +  ' Averaged superimposed coordinates written to ',
     + filcop(1:index(filcop,' ')-1)
      else
      write(*,'(/,a,a)')
     +  ' Averaged coordinates written to ',
     + filcop(1:index(filcop,' ')-1)
      endif
      endif
      goto 1
      endif
c
c     build covariance matrix
c     -----------------------
c
      if(card(1:6).eq.'covar ') then
      fil=card(8:80)
      if(nsel.eq.0) call md_abort('No atoms selected',0)
c
      if(active) call ana_edfinal()
      active=.true.
c
      call ana_rdhdr(byte_mb(i_snam))
      nave=0
      nread=0
      nframe=0
      call ana_edinit()
    5 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 5
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      call super2(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms),
     + nsa,msa,.true.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.true.)
      call ana_xcollect(dbl_mb(i_xdat),ifr,num,msa)
c
      if(np.gt.1) call ana_update(dbl_mb(i_xdat),msa,ifr,ifr+num-1)
c
      call ana_covar(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
     + dbl_mb(i_xdat),dbl_mb(i_xadat))
c
      nave=nave+1
c
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 5
      endif
c
      call ga_sync()
      call ga_scale(ga_cov,1.0d0/dble(nave))
c
      if(me.eq.0) then
      if(util_print('frames_read',print_debug)) then
      write(*,'(/,a)') ' Covariance matrix generation completed'
      endif
      endif
c
      call ga_sync()
c      call ga_print(ga_cov)
      call ga_diag_std(ga_cov,ga_vec,dbl_mb(i_cov))
c
      call ga_sync()
      do 19 i=1,3*nsel
      int_mb(i_ord+i-1)=i
   19 continue
      do 17 i=1,3*nsel-1
      do 18 j=i+1,3*nsel
      if(abs(dbl_mb(i_cov-1+int_mb(i_ord+i-1))).lt.
     + abs(dbl_mb(i_cov-1+int_mb(i_ord+j-1)))) then
      k=int_mb(i_ord+i-1)
      int_mb(i_ord+i-1)=int_mb(i_ord+j-1)
      int_mb(i_ord+j-1)=k 
      endif
   18 continue
   17 continue
c
      if(me.eq.0) then
      write(*,'(/,a,i10,a)') ' Covariance matrix diagonalized for ',
     + nave,' frames'
      write(*,'(/,a,//,(5x,6e12.5))') ' Eigenvalues',
     + (dbl_mb(i_cov-1+int_mb(i_ord+i-1)),i=1,3*nsel)
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
      goto 1
      endif
c
c     projection
c     ----------
c
      if(card(1:6).eq.'projec') then
      lesd=.true.
      read(card(8:14),'(i7)') ivec
      fil=card(15:80)
      i=index(fil,'.')
      if(i.eq.0) then
      fmt='trj'
      filprj=fil(1:index(fil,' ')-1)//'.trj'
      filval=fil(1:index(fil,' ')-1)//'.val'
      filvec=fil(1:index(fil,' ')-1)//'.vec'
      filtcl=fil(1:index(fil,' ')-1)//'.tcl'
      filecc=fil(1:index(fil,' ')-1)//'.ecce.out'
      filmin=fil(1:index(fil,' ')-1)//'_min.pdb'
      filmax=fil(1:index(fil,' ')-1)//'_max.pdb'
      filamp=fil(1:index(fil,' ')-1)//'.amp'
      else
      filprj=fil
      fmt=filprj(i+1:i+3)
      filval=fil(1:i-1)//'.val'
      filvec=fil(1:i-1)//'.vec'
      filtcl=fil(1:i-1)//'.tcl'
      filecc=fil(1:i-1)//'.ecce.out'
      filmin=fil(1:i-1)//'_min.pdb'
      filmax=fil(1:i-1)//'_max.pdb'
      filamp=fil(1:i-1)//'.amp'
      endif
c
      if(nsel.eq.0) call md_abort('No atoms selected',0)
      if(.not.active) call md_abort('No ED matrix active',0)
c
      if(me.eq.0) then
      open(unit=lfnprj,file=filprj(1:index(filprj,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnval,file=filval(1:index(filval,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
c
      call ana_wthdr(lfnprj,fmt,byte_mb(i_snam),byte_mb(i_tag),
     + int_mb(i_isel),.false.)
c
      k=int_mb(i_ord+ivec-1)
      call ga_get(ga_vec,1,3*nsel,k,k,dbl_mb(i_cov),1)
c
      call ana_rdhdr(byte_mb(i_snam))
      nave=0
      nread=0
      nframe=0
c
      valmin=1.0d0
      valmax=-1.0d0
c
    6 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 6
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
c      call super2(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
c     + dbl_mb(i_xref),dbl_mb(i_swt),dbl_mb(i_wt),dbl_mb(i_xrms),
c     + nsa,msa,.true.,rms0,rms1,dbl_mb(i_wdat),mwm,nwm,mwa,nwa,.true.)
      call ana_xcollect(dbl_mb(i_xdat),ifr,num,msa)
c
      if(np.gt.1) call ana_update(dbl_mb(i_xdat),msa,ifr,ifr+num-1)
c
      call ana_projec(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
     + dbl_mb(i_xdat),dbl_mb(i_xadat),msa,val,
     + dbl_mb(i_xp))
c
      lxw=.false.
      call ana_wtfram(lfnprj,fmt,dbl_mb(i_xp),byte_mb(i_snam),
     + dbl_mb(i_wdat),int_mb(i_isel),byte_mb(i_tag),dbl_mb(i_val),
     + int_mb(i_wsel),int_mb(i_ndxw))
c
      nave=nave+1
c
      if(me.eq.0) then
      write(lfnval,'(i5,f12.6)') nave,val
      endif
c
      endif
      nframe=nframe+1
      if(ifrto.lt.ifrfr.or.nread.lt.ifrto) goto 6
      endif
c
      if(.not.ana_wrtvec(ivec,fil,int_mb(i_ndx),dbl_mb(i_xadat)))
     + goto 9999
c
      lesd=.false.
c
      goto 1
      endif
c
c     make report from file
c     ---------------------
c
      if(card(1:6).eq.'report') then
      option=card(8:15)
      fil=card(16:80)
      call ana_report_file(option,fil,lfnloc)
      endif
c
c     generate esp grid from an xyzq file
c     -----------------------------------
c
      if(card(1:6).eq.'espplt') then
c
      string=card(8:80)
      read(string(1:4),'(i4)') ngrid
      read(string(5:8),'(i4)') iesppb
      read(string(9:20),'(f12.6)') rcut
      string=card(28:80)
      i=index(string,' ')
      filchg=string(1:i)
      filplt=string(i+1:80)
      filepz=filplt(1:index(filplt,'.'))//'epz '
c
      if(filchg(1:1).ne.' ') call ana_pltgrd(.true.,int_mb(i_isel))
c
      goto 1
      endif
c
      if(me.eq.0) close(unit=lfncmd,status='delete')
c      if(me.eq.0) close(unit=lfncmd,status='keep')
c
      if(active) call ana_edfinal()
c
      if(.not.ma_pop_stack(l_sbnd))
     + call md_abort('Could not deallocate sbnd',0)
      if(.not.ma_pop_stack(l_imol))
     + call md_abort('Could not deallocate imol',0)
      if(.not.ma_pop_stack(l_iram))
     + call md_abort('Could not deallocate iram',0)
c
      call ana_finish(byte_mb(i_snam),.true.)
c
      if(.not.ma_pop_stack(l_wrk))
     + call md_abort('Could not deallocate iwrk',0)
c
      ignore=rtdb_parallel(.true.)
c
      return
c
 9999 continue
      call md_abort('Error opening command file',0)
c
      return
      end
      subroutine ana_update(x,m,ifr,ito)
c
      implicit none
c
      integer m,ifr,ito
      real*8 x(m,3)
c
#include "global.fh"
#include "msgids.fh"
c
      integer i,j
c
      if(ifr.gt.1) then
      do 1 j=1,3
      do 2 i=1,ifr-1
      x(i,j)=0.0d0
    2 continue
    1 continue
      endif
      if(ito.lt.m) then
      do 3 j=1,3
      do 4 i=ito+1,m
      x(i,j)=0.0d0
    4 continue
    3 continue
      endif
c
      call ga_dgop(mag_d12,x,3*m,'+')
c
      return
      end
      subroutine ana_anahdr(sgmnam)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
c
      character*16 sgmnam(nsa)
      character*11 atom(4)
      integer i,j
c
      write(lfnana,1000) nbnds,nangs,ntors,nimps
 1000 format(4i10)
c
      do 1 i=1,nbnds
      do 11 j=1,2
      if(ibnds(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(ibnds(i,j,1))(12:16),':',
     + sgmnam(ibnds(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+ibnds(i,j,1),':',
     + wnam(ibnds(i,j,2))(6:10)
      endif
   11 continue
c
      write(lfnana,1001) (atom(j),j=1,2)
 1001 format(a11,'-',a11)
    1 continue
c
      do 2 i=1,nangs
      do 21 j=1,3
      if(iangs(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(iangs(i,j,1))(12:16),':',
     + sgmnam(iangs(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+iangs(i,j,1),':',
     + wnam(iangs(i,j,2))(6:10)
      endif
   21 continue
      write(lfnana,1002) (atom(j),j=1,3)
 1002 format(a11,'-',a11,'-',a11)
    2 continue
c
      do 3 i=1,ntors
      do 31 j=1,4
      if(itors(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(itors(i,j,1))(12:16),':',
     + sgmnam(itors(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+itors(i,j,1),':',
     + wnam(itors(i,j,2))(6:10)
      endif
   31 continue
      write(lfnana,1003) (atom(j),j=1,4)
 1003 format(a11,'-',a11,'-',a11,'-',a11)
    3 continue
c
      do 4 i=1,nimps
      do 41 j=1,4
      if(iimps(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(iimps(i,j,1))(12:16),':',
     + sgmnam(iimps(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+iimps(i,j,1),':',
     + wnam(iimps(i,j,2))(6:10)
      endif
   41 continue
      write(lfnana,1003) (atom(j),j=1,4)
    4 continue
c
      return
      end
      logical function ana_wrtvec(ivec,fil,ndx,xa)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ivec
      character*255 fil
      integer i,j,k,l
      real*8 value,xa(msa,3)
      integer ndx(nsa)
c
      if(me.eq.0) then
      open(unit=lfntcl,file=filtcl(1:index(filtcl,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnvec,file=filvec(1:index(filvec,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnecc,file=filecc(1:index(filecc,' ')-1),
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnamp,file=filamp(1:index(filamp,' ')-1),
     + form='formatted',status='unknown',err=9999)
      k=0
      l=0
      do 15 i=1,nsel
      value=0.0d0
      do 16 j=1,3
      value=value+(dbl_mb(i_cov+k))**2
      k=k+1
   16 continue
      l=l+1
      write(lfntcl,5555)
     + (10.0*xa(ndx(i),j),j=1,3),
     + (dbl_mb(i_cov+k-4+j),j=1,3)
 5555 format('draw vector2 {',3f12.6,'} {',3f12.6,'} 1.0 10')
      write(lfnamp,'(i5,f12.6)') i,sqrt(value)
   15 continue
c
      write(lfnecc,3013)
 3013 format('%begin%input file%1%char',/,
     + 'geometry noautosym units nm')
      do 314 i=1,nsa
      if(int_mb(i_osel+i-1).eq.1) then
      if(byte_mb(i_snam+(i-1)*16+5).eq.'2'.or.
     + byte_mb(i_snam+(i-1)*16+5).eq.'3'.or.
     + byte_mb(i_snam+(i-1)*16+5).eq.'4') then
      write(lfnecc,3014) ' ',byte_mb(i_snam+(i-1)*16+6),
     + (dbl_mb(i_xadat+i-1+j*msa),j=0,2)
      else
      write(lfnecc,3014) (byte_mb(i_snam+(i-1)*16+j),j=5,6),
     + (dbl_mb(i_xadat+i-1+j*msa),j=0,2)
      endif
 3014 format(2a1,1x,3f12.6)
      endif
  314 continue
      write(lfnecc,3012) nselo
 3012 format('end',/,'basis "ao basis" spherical print',/,
     + ' C library "6-31+G*"',/,
     + ' N library "6-31+G*"',/,
     + ' O library "6-31+G*"',/,
     + ' H library "6-31+G*"',/,
     + ' S library "6-31+G*"',/,
     + ' P library "6-31+G*"',/,
     + ' Ca library "6-31+G*"',/,
     + ' Na library "6-31+G*"',/,
     + 'end',/,
     + 'task scf gradient',/,
     + '%end%input file%1%char',/,
     + 'task_gradient%begin%total gradient%3 ',i7,'%double')
c
      k=0
      do 315 i=1,nsa
      if(int_mb(i_osel+i-1).eq.1) then
      if(int_mb(i_isel+i-1).eq.0) then
      write(lfnecc,3015) 0.0d0,0.0d0,0.0d0
 3015 format(3f12.6)
      else
      write(lfnecc,3015) (dbl_mb(i_cov+k+j),j=0,2)
      k=k+3
      endif
      endif
  315 continue
      write(lfnecc,3011) nselo
 3011 format('task_gradient%end%total gradient%3 ',i7,'%double')
c
      close(unit=lfntcl)
      close(unit=lfnvec)
      close(unit=lfnecc)
      close(unit=lfnval)
      close(unit=lfnprj)
      close(unit=lfnamp)
      write(*,'(/,a,i5,a,a)') ' Projection of vector ',ivec,
     +' written to ',fil(1:index(fil,' ')-1)
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
      ana_wrtvec=.true.
      return
 9999 continue
      ana_wrtvec=.false.
      return
      end
      subroutine ana_xcollect(x,ifrom,nx,mx)
c
      implicit none
c
#include "ana_common.fh"
#include "ana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ifrom,nx,mx
      real*8 x(mx,3)
c
      integer i,j
c
      if(np.gt.1) then
c
      do 1 j=1,3
      do 2 i=1,ifrom-1
      x(i,j)=zero
    2 continue
    1 continue
c
      do 3 j=1,3
      do 4 i=ifrom+nx+1,mx
      x(i,j)=zero
    4 continue
    3 continue
c
      call ga_dgop(mag_d08,x,3*mx,'+')
c
      endif
c
      return
      end
     
