      subroutine ana_task
c
c $Id: ana_task.F,v 1.8 2000-01-25 19:53:35 d3j191 Exp $
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
      integer lfncmd,lfn
      character*255 filcmd,fil
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
      real*8 ana_torsion
      external ana_rdfram,ana_index,ana_torsion
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
      lfn=83
c
      ifrfr=1
      ifrto=1
      ifrsk=1
      ifrst=0
      ilast=0
c
      ntors=0
c
      lsel=.false.
      active=.false.
c
      nave=0
      time=0.0d0
      timr=0.0d0
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
      do 7 i=1,4
      read(lfncmd,1000,end=9999) card
      itors(it,i)=ana_index(card,byte_mb(i_snam))
    7 continue
      write(*,'(4i5)') (itors(it,i),i=1,4)
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
      read(card(8:22),'(3i5)') ifrfr,ifrto,ifrsk
      if(me.eq.0) then
      write(*,'(/,3(a,i5))') ' Selected frames ',ifrfr,' to ',ifrto,
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
      call ana_wthdr(lfncop,fmt,byte_mb(i_snam),byte_mb(i_tag))
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
     + dbl_mb(i_wdat))
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
      call ana_wthdr(lfncop,fmt,byte_mb(i_snam),byte_mb(i_tag))
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
     + dbl_mb(i_wdat))
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
      read(card(8:80),'(a)') filcop
      i=index(filcop,'.')
      fmt=filcop(i+1:i+3)
c
c     read file header
c
      call ana_rdhdr(byte_mb(i_snam))
c
      if(me.eq.0) then
      open(unit=lfncop,file=filcop(1:index(filcop,' ')-1),
     + form='formatted',status='unknown',err=9999)
      write(*,3333) filcop(1:index(filcop,' ')-1)
 3333 format(/,' Opening scn file ',a)
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
      do 11 i=1,ntors
      rtors(i)=ana_torsion(dbl_mb(i_xdat),msa,
     + itors(i,1),itors(i,2),itors(i,3),itors(i,4))
   11 continue
      write(lfncop,2323) time,(rtors(i),i=1,ntors)
 2323 format(5f12.6)
      nwrit=nwrit+1
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 10
      endif
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
c
      if(me.eq.0) then
      close(unit=lfncop)
      write(*,'(a)') ' Closing trj file '
      write(*,'(/,a,i10,a,a)') ' Scanned ',nwrit,' frames to ',
     + filcop(1:index(filcop,' ')-1)
      endif
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
      call ga_diag_std(ga_cov,ga_vec,dbl_mb(i_cov))
c
      if(me.eq.0) then
      write(*,'(/,a,i10,a)') ' Covariance matrix diagonalized for ',
     + nave,' frames'
      write(*,'(/,a,//,(5x,6e12.5))') ' Eigenvalues',
     + (dbl_mb(i_cov-1+i),i=3*nsel,1,-1)
      if(fil(1:1).ne.' ') then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='unknown',err=9999)
      write(lfn,'(i5,6e12.5)') (i,dbl_mb(i_cov-1+i),i=3*nsel,1,-1)
      close(unit=lfn)
      endif
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
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
c
      call ga_get(ga_vec,1,nsel,nsel-ivec+1,nsel-ivec+1,
     + dbl_mb(i_cov),nsel)
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
     + dbl_mb(i_xdat+ifr-1),dbl_mb(i_xadat+ifr-1),val)
c
      nave=nave+1
c
      if(me.eq.0) then
      write(lfn,'(i5,f12.6)') nave,val
      endif
c
      endif
      nframe=nframe+1
      if(ifrto.le.ifrfr.or.nread.lt.ifrto) goto 6
      endif
      if(me.eq.0) then
      close(unit=lfn)
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

