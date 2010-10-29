      subroutine dia_init
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
#include "dia_common.fh"
#include "dia_params.fh"
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
      lfncub=59
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
      subroutine dia_getcmd()
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      lesd=.false.
c
      if(me.eq.0) then
      read(lfncmd,1000,end=9999) card
 1000 format(a)
      endif
      call util_char_ga_brdcst(mag_d00,card,0)
      cmd=card(1:6)
c
      return
 9999 continue
      call md_abort('Error readin command file',0) 
      return
      end
      logical function dia_getcard()
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(me.eq.0) then
      read(lfncmd,1000,end=9999) card
 1000 format(a)
      endif
c
      dia_getcard=.true.
      return
c
 9999 continue
      dia_getcard=.false.
c
      return
      end
      subroutine dia_final()
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
c 
      logical lopen
c
      if(me.eq.0) close(unit=lfncmd,status='delete')
c      if(me.eq.0) close(unit=lfncmd,status='keep')
c
      if(active) call dia_edfinal()
c
      if(.not.ma_pop_stack(l_sbnd))
     + call md_abort('Could not deallocate sbnd',0)
      if(.not.ma_pop_stack(l_iequi))
     + call md_abort('Could not deallocate iequi',0)
      if(.not.ma_pop_stack(l_imol))
     + call md_abort('Could not deallocate imol',0)
      if(.not.ma_pop_stack(l_iram))
     + call md_abort('Could not deallocate iram',0)
c
      call dia_finish(byte_mb(i_snam),.true.)
c
      if(.not.ma_pop_stack(l_wrk))
     + call md_abort('Could not deallocate iwrk',0)
c
      inquire(unit=lfnpca,opened=lopen)
      if(lopen) close(lfnpca,status='keep')
      inquire(unit=lfnord,opened=lopen)
      if(lopen) close(lfnord,status='keep')
      inquire(unit=lfndat,opened=lopen)
      if(lopen) close(lfndat,status='keep')
c
      ignore=rtdb_parallel(.true.)
c
      return
      end
      subroutine dia_update(x,m,ifr,ito)
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
      subroutine dia_anahdr(sgmnam)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
c
      character*16 sgmnam(nsa)
      character*11 atom(4)
      integer i,j
c
      write(lfnana,1000) ndist,nangl,ntors,nimps
 1000 format(4i10)
c
      do 1 i=1,ndist
      do 11 j=1,2
      if(idist(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(idist(i,j,1))(12:16),':',
     + sgmnam(idist(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+idist(i,j,1),':',
     + wnam(idist(i,j,2))(6:10)
      endif
   11 continue
c
      write(lfnana,1001) (atom(j),j=1,2)
 1001 format(a11,'-',a11)
    1 continue
c
      do 2 i=1,nangl
      do 21 j=1,3
      if(iangl(i,j,2).eq.0) then
      write(atom(j),'(a5,a,a5)') sgmnam(iangl(i,j,1))(12:16),':',
     + sgmnam(iangl(i,j,1))(6:10)
      else
      write(atom(j),'(i5,a,a5)') nsgm+iangl(i,j,1),':',
     + wnam(iangl(i,j,2))(6:10)
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
      logical function dia_wrtvec(ivec,fil,ndex,xa)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
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
      integer ndex(nsa)
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
     + (10.0*xa(ndex(i),j),j=1,3),
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
      dia_wrtvec=.true.
      return
 9999 continue
      dia_wrtvec=.false.
      return
      end
      subroutine dia_xcollect(x,ifrom,nx,mx)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
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
      subroutine dia_test(x,nsa)
      integer nsa
      real*8 x(nsa,3)
c
      integer i,j
c
      do 1 i=1,10
      write(*,'(i5,3f12.6)') i,(x(i,j),j=1,3)
    1 continue
c
      return
      end
      subroutine dia_test2(xtest,nxdim,nfr)
#include "dia_common.fh"
      integer nxdim
      real*8 xtest(nxdim,3,nfr)
c
      integer i,j,k
c
      do 2 k=1,nfr
      write(*,'(a,i5)') 'FRAME ',me*nfr+k
      do 1 i=1,10
      write(*,'(i5,3f12.6)') i,(xtest(i,j,k),j=1,3)
    1 continue
    2 continue
c
      return
      end
