      subroutine ana_task
c
c $Id: ana_task.F,v 1.16 2000-03-21 22:23:50 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
c
      logical ignore
c
      integer lfncmd,lfnprj,lfnval,lfncov
      character*255 filcmd,fil,string
c
      character*80 card
      integer nsai,msai,ifrsk,it
      integer nframe,nread,nwrit,ifr,ito,num,i,ivec
      real*8 rms0,rms1,val
      character*3 fmt
      logical active
c
      logical ana_rdfram
      integer ana_index
      real*8 ana_bond,ana_angle,ana_torsion
      external ana_rdfram,ana_index
      external ana_bond,ana_angle,ana_torsion
c
c
c     set sequential access to rtdb
c     and ignore current state
c
      ignore=rtdb_parallel(.false.)
c
      lfncmd=69
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
      lfncov=86
c
      ifrfr=1
      ifrto=1
      ifrsk=1
      ifrst=0
      ilast=0
c
      nbnds=0
      nangs=0
      ntors=0
      nimps=0
c
      lsel=.false.
      active=.false.
      lrms=.false.
      lana=.false.
c
      nave=0
      ndata=0
      time=0.0d0
      timr=0.0d0
c
      ifr=1
c
      if(me.eq.0) then
      if(.not.rtdb_cget(irtdb,'ana:filcmd',1,filcmd))
     + call errquit('Error rtdb_get filcmd',0)
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=9999)
      rewind(unit=lfncmd)
      endif
c
    1 continue
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
      call ana_sizref(nsai,msai)
c
c     allocate memory and initialize
c
      call ana_init(nsai,msai,.true.)
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
      if(me.eq.0) then
      write(*,'(a,a)') ' Reference coordinates read from ',
     + filref(1:index(filref,' ')-1)
      write(*,'(/,a,i5)') ' Number of atoms is ',nsai
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
      call ana_rdhdr(byte_mb(i_snam))
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      call ana_all(int_mb(i_isel),1,dbl_mb(i_wt))
      lsel=.false.
c
      if(me.eq.0) then
      write(*,'(/,a,a)') ' Trajectory file header from ',
     + filtrj(1:index(filtrj,' ')-1)
      endif
      goto 1
      endif
c
c     atom selection
c     --------------
c
      if(card(1:6).eq.'select') then
      if(.not.lsel) then
      call ana_all(int_mb(i_isel),0,dbl_mb(i_wt))
      lsel=.true.
      endif
      call ana_select(card,byte_mb(i_snam),int_mb(i_isel),dbl_mb(i_wt),
     + dbl_mb(i_xref))
      if(me.eq.0) then
      write(*,3003) nsel,nsa
 3003 format(' Selected',i5,' out of',i5,' atoms')
      endif
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
c     bond specification
c     ------------------
c
      if(card(1:6).eq.'bond  ') then
      read(card(8:12),'(i5)') it
      if(it.eq.nbnds+1) nbnds=nbnds+1
      if(it.gt.nbnds) call errquit('Bond index too large',it)
      if(it.gt.mxbnds) call errquit('Bond index too large',it)
      do 7 i=1,2
      read(lfncmd,1000,end=9999) card
      ibnds(it,i)=ana_index(card,byte_mb(i_snam))
    7 continue
      write(*,'(2i5)') (ibnds(it,i),i=1,2)
      goto 1
      endif
c
c     angle specification
c     -------------------
c
      if(card(1:5).eq.'angle') then
      read(card(8:12),'(i5)') it
      if(it.eq.nangs+1) nangs=nangs+1
      if(it.gt.nangs) call errquit('Angle index too large',it)
      if(it.gt.mxangs) call errquit('Angle index too large',it)
      do 8 i=1,3
      read(lfncmd,1000,end=9999) card
      iangs(it,i)=ana_index(card,byte_mb(i_snam))
    8 continue
      write(*,'(3i5)') (iangs(it,i),i=1,3)
      goto 1
      endif
c
c     torsion specification
c     ---------------------
c
      if(card(1:6).eq.'torsio') then
      read(card(8:12),'(i5)') it
      if(it.eq.ntors+1) ntors=ntors+1
      if(it.gt.ntors) call errquit('Torsion index too large',it)
      if(it.gt.mxtors) call errquit('Torsion index too large',it)
      do 9 i=1,4
      read(lfncmd,1000,end=9999) card
      itors(it,i)=ana_index(card,byte_mb(i_snam))
    9 continue
      write(*,'(4i5)') (itors(it,i),i=1,4)
      goto 1
      endif
c
c     improper specification
c     ----------------------
c
      if(card(1:6).eq.'improp') then
      read(card(8:12),'(i5)') it
      if(it.eq.nimps+1) nimps=nimps+1
      if(it.gt.nimps) call errquit('Improper index too large',it)
      if(it.gt.mximps) call errquit('Improper index too large',it)
      do 12 i=1,4
      read(lfncmd,1000,end=9999) card
      iimps(it,i)=ana_index(card,byte_mb(i_snam))
   12 continue
      write(*,'(4i5)') (iimps(it,i),i=1,4)
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
      goto 1
      endif
c
c     copy frames
c     -----------
c
      if(card(1:6).eq.'copy  ') then
      read(card(8:80),'(a)') filcop
      i=index(filcop,'.')
      fmt=filcop(i+1:i+3)
c
c     write header
c
      call ana_rdhdr(byte_mb(i_snam))
      call ana_wthdr(lfncop,fmt,byte_mb(i_snam),byte_mb(i_tag),
     + int_mb(i_isel))
c
      nread=0
      nwrit=0
      nframe=0
    2 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nframe.lt.ifrfr) goto 2
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      call ana_wtfram(lfncop,fmt,dbl_mb(i_xdat),byte_mb(i_snam),
     + dbl_mb(i_wdat),int_mb(i_isel))
      nwrit=nwrit+1
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 2
      endif
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
c
      if(me.eq.0) then
      close(unit=lfncop)
      write(*,'(/,a,i10,a,a)') ' Copied ',nwrit,' frames to ',
     + filcop(1:index(filcop,' ')-1)
      endif
c
      goto 1
      endif
c
c     super frames
c     -----------
c
      if(card(1:6).eq.'super ') then
      read(card(8:80),'(a)') filcop
      i=index(filcop,'.')
      fmt=filcop(i+1:i+3)
c
      if(nsel.eq.0) call errquit('No atoms selected',0)
c
c     write header
c
      call ana_rdhdr(byte_mb(i_snam))
      call ana_wthdr(lfncop,fmt,byte_mb(i_snam),byte_mb(i_tag),
     + int_mb(i_isel))
c
      nread=0
      nwrit=0
      nframe=0
    3 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 3
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      call super(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
     + dbl_mb(i_xref),dbl_mb(i_wt),dbl_mb(i_xrms),nsa,msa,.true.,
     + rms0,rms1)
      call ana_wtfram(lfncop,fmt,dbl_mb(i_xdat),byte_mb(i_snam),
     + dbl_mb(i_wdat),int_mb(i_isel))
      nwrit=nwrit+1
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 3
      endif
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
c
c
      if(me.eq.0) then
      close(unit=lfncop)
      write(*,'(/,a,i10,a,a)') ' Superimposed ',nwrit,' frames to ',
     + filcop(1:index(filcop,' ')-1)
      endif
c
      goto 1
      endif
c
c     scan frames
c     -----------
c
      if(card(1:6).eq.'scan  ') then
      read(card(8:80),'(a)') string
      filrms=string(1:index(string,' ')-1)//'.rms'
      filana=string(1:index(string,' ')-1)//'.ana'
      write(*,'(a,a)') 'filrms=',filrms(1:index(filrms,' ')-1)
      time=0.0d0
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
      write(lfnana,2324) nbnds,nangs,ntors,nimps
 2324 format(4i10)
      endif
      if(lrms) then
      open(unit=lfnrms,file=filrms(1:index(filrms,' ')-1),
     + form='formatted',status='unknown')
      endif
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
c     evaluate bond lenghts
c
      do 11 i=1,nbnds
      rbnds(i)=ana_bond(dbl_mb(i_xdat),msa,
     + ibnds(i,1),ibnds(i,2))
   11 continue
c
c     evaluate angles
c
      do 13 i=1,nangs
      rangs(i)=ana_angle(dbl_mb(i_xdat),msa,
     + iangs(i,1),iangs(i,2),iangs(i,3))
   13 continue
c
c     evaluate torsions
c
      do 14 i=1,ntors
      rtors(i)=ana_torsion(dbl_mb(i_xdat),msa,
     + itors(i,1),itors(i,2),itors(i,3),itors(i,4))
   14 continue
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
c     root mean square deviation
c
      if(lrms) then
      if(np.gt.1) then
      call ana_rmsdev(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,
     + dbl_mb(i_xref),dbl_mb(i_wt),dbl_mb(i_xrms))
      else
      call ana_rmsdev(dbl_mb(i_xdat),int_mb(i_idat),nsa,
     + dbl_mb(i_xref),dbl_mb(i_wt),dbl_mb(i_xrms))
      endif
      ndata=ndata+1
      endif
c
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 10
      endif
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
c
      if(me.eq.0) then
      close(unit=lfnana)
      write(*,'(a)') ' Closing scn file '
      write(*,'(/,a,i10,a,a)') ' Scanned ',nwrit,' frames to ',
     + filcop(1:index(filana,' ')-1)
      endif
c
      nbnds=0
      nangs=0
      ntors=0
      nimps=0
      lrms=.false.
      lana=.false.
c
      goto 1
      endif
c
c     calculate average coordinates
c     -----------------------------
c
      if(card(1:6).eq.'xaver ') then
c
      if(nsel.eq.0) call errquit('No atoms selected',0)
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
      call super(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
     + dbl_mb(i_xref),dbl_mb(i_wt),dbl_mb(i_xrms),nsa,msa,.true.,
     + rms0,rms1)
c
      call ana_xaver(dbl_mb(i_xdat+ifr-1),dbl_mb(i_xadat+ifr-1),
     + num,msa,nave)
c
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 4
      endif
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      call ana_xmean(dbl_mb(i_xadat+ifr-1),num,msa,nave)
c
      if(me.eq.0) then
      write(*,'(/,a,i10,a)') ' Averaged ',nave,' frames'
      endif
c
      goto 1
      endif
c
c     build covariance matrix
c     -----------------------
c
      if(card(1:6).eq.'covar ') then
      fil=card(8:80)
      if(nsel.eq.0) call errquit('No atoms selected',0)
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
      call super(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
     + dbl_mb(i_xref),dbl_mb(i_wt),dbl_mb(i_xrms),nsa,msa,.true.,
     + rms0,rms1)
c
      call ana_covar(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
     + dbl_mb(i_xdat+ifr-1),dbl_mb(i_xadat+ifr-1))
c
      nave=nave+1
c
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 5
      endif
c
      call ga_scale(ga_cov,1.0d0/dble(nave))
c
      call ga_print(ga_cov)
      call ga_diag_std(ga_cov,ga_vec,dbl_mb(i_cov))
      call ga_print(ga_vec)
c
      if(me.eq.0) then
      write(*,'(/,a,i10,a)') ' Covariance matrix diagonalized for ',
     + nave,' frames'
      write(*,'(/,a,//,(5x,6e12.5))') ' Eigenvalues',
     + (dbl_mb(i_cov-1+i),i=3*nsel,1,-1)
c      if(fil(1:1).ne.' ') then
c      open(unit=lfncov,file=fil(1:index(fil,' ')-1),
c     + form='formatted',status='unknown',err=9999)
c      write(lfncov,'(i5,6e12.5)') (i,dbl_mb(i_cov-1+i),i=3*nsel,1,-1)
c      close(unit=lfncov)
c      endif
      endif
c
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
c
      goto 1
      endif
c
c     projection
c     ----------
c
      if(card(1:6).eq.'projec') then
      read(card(8:14),'(i7)') ivec
      fil=card(15:80)
      if(nsel.eq.0) call errquit('No atoms selected',0)
      if(.not.active) call errquit('No ED matrix active',0)
c
      if(me.eq.0) then
      open(unit=lfnprj,file=fil(1:index(fil,' ')-1)//'.prj',
     + form='formatted',status='unknown',err=9999)
      open(unit=lfnval,file=fil(1:index(fil,' ')-1)//'.val',
     + form='formatted',status='unknown',err=9999)
      endif
      fmt='trj'
      call ana_wthdr(lfnprj,fmt,byte_mb(i_snam),byte_mb(i_tag),
     + int_mb(i_isel))
c
      call ga_get(ga_vec,1,3*nsel,ivec,ivec,dbl_mb(i_cov),nsel)
c
      call ana_rdhdr(byte_mb(i_snam))
      nave=0
      nread=0
      nframe=0
    6 continue
      if(ana_rdfram(dbl_mb(i_xdat),int_mb(i_idat),dbl_mb(i_wdat))) then
      nread=nread+1
      if(nread.lt.ifrfr) goto 6
      if(nframe.eq.(nframe/ifrsk)*ifrsk) then
      call super(dbl_mb(i_xdat+ifr-1),int_mb(i_idat+ifr-1),num,msa,
     + dbl_mb(i_xref),dbl_mb(i_wt),dbl_mb(i_xrms),nsa,msa,.true.,
     + rms0,rms1)
c
      call ana_projec(int_mb(i_ndx),dbl_mb(i_cov),int_mb(i_isel),
     + dbl_mb(i_xdat+ifr-1),dbl_mb(i_xadat+ifr-1),val,
     + dbl_mb(i_xp+ifr-1))
c
      lxw=.false.
      call ana_wtfram(lfnprj,fmt,dbl_mb(i_xp),byte_mb(i_snam),
     + dbl_mb(i_wdat),int_mb(i_isel))
c
      nave=nave+1
c
      if(me.eq.0) then
      write(lfnval,'(i5,f12.6)') nave,val
      endif
c
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 6
      endif
      if(me.eq.0) then
      close(unit=lfnval)
      close(unit=lfnprj)
      write(*,'(/,a,i5,a,a)') ' Projection of vector ',ivec,
     +' written to ',fil(1:index(fil,' ')-1)
      endif
c
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      endif
c
c     generate esp grid from an xyzq file
c     -----------------------------------
c
      if(card(1:6).eq.'espplt') then
      call ana_pltgrd(card)
      goto 1
      endif
c
      close(unit=lfncmd,status='delete')
c
      if(active) call ana_edfinal()
      call ana_finish(byte_mb(i_snam),.true.)
c
      ignore=rtdb_parallel(.true.)
c
      return
c
 9999 continue
      call errquit('Error opening command file',0)
c
      return
      end

