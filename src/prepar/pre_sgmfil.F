      logical function pre_sgmfil(sgmnam,lfnsgm,filsgm,
     + dir_s,dir_x,dir_q,dir_u,dir_t,dir_c)
c
c $Id: pre_sgmfil.F,v 1.5 2001-11-15 16:45:57 d3j191 Exp $
c
      implicit none
c
      integer lfnsgm
      character*10 sgmnam
      character*255 filsgm,dir_s,dir_x,dir_q,dir_u,dir_t,dir_c
c
      integer lens,lend,lenf
c
      lens=index(sgmnam,' ')-1
      if(lens.le.0) lens=10
c
c     try dir_c
c     ---------
c
      lend=index(dir_c,' ')-1
      if(lend.gt.0) then
      filsgm=dir_c(1:lend)//sgmnam(1:lens)//'.sgm '
      else
      filsgm=sgmnam(1:lens)//'.sgm '
      endif
      lenf=index(filsgm,' ')-1
c
      open(unit=lfnsgm,file=filsgm(1:lenf),form='formatted',
     + status='old',err=1)
      close(unit=lfnsgm)
      pre_sgmfil=.true.
      return
c
    1 continue
c
c     try dir_t
c     ---------
c
      lend=index(dir_t,' ')-1
      if(lend.le.0) goto 2
c
      filsgm=dir_t(1:lend)//sgmnam(1:lens)//'.sgm '
      lenf=index(filsgm,' ')-1
c
      open(unit=lfnsgm,file=filsgm(1:lenf),form='formatted',
     + status='old',err=2)
      close(unit=lfnsgm)
      pre_sgmfil=.true.
      return
c
    2 continue
c
c     try dir_u
c     ---------
c
      lend=index(dir_u,' ')-1
      if(lend.le.0) goto 3
c
      filsgm=dir_u(1:lend)//sgmnam(1:lens)//'.sgm '
      lenf=index(filsgm,' ')-1
c
      open(unit=lfnsgm,file=filsgm(1:lenf),form='formatted',
     + status='old',err=3)
      close(unit=lfnsgm)
      pre_sgmfil=.true.
      return
c
    3 continue
c
c     try dir_q
c     ---------
c
      lend=index(dir_q,' ')-1
      if(lend.le.0) goto 4
c
      filsgm=dir_q(1:lend)//sgmnam(1:lens)//'.sgm '
      lenf=index(filsgm,' ')-1
c
      open(unit=lfnsgm,file=filsgm(1:lenf),form='formatted',
     + status='old',err=4)
      close(unit=lfnsgm)
      pre_sgmfil=.true.
      return
c
    4 continue
c
c     try dir_x
c     ---------
c
      lend=index(dir_x,' ')-1
      if(lend.le.0) goto 5
c
      filsgm=dir_x(1:lend)//sgmnam(1:lens)//'.sgm '
      lenf=index(filsgm,' ')-1
c
      open(unit=lfnsgm,file=filsgm(1:lenf),form='formatted',
     + status='old',err=5)
      close(unit=lfnsgm)
      pre_sgmfil=.true.
      return
c
    5 continue
c
c     try dir_s
c     ---------
c
      lend=index(dir_s,' ')-1
      if(lend.le.0) goto 6
c
      filsgm=dir_s(1:lend)//sgmnam(1:lens)//'.sgm '
      lenf=index(filsgm,' ')-1
c
      open(unit=lfnsgm,file=filsgm(1:lenf),form='formatted',
     + status='old',err=6)
      close(unit=lfnsgm)
      pre_sgmfil=.true.
      return
c
    6 continue
c
c     return .false. with dir_t filname
c     ---------------------------------
c
      lend=index(dir_t,' ')-1
c
      if(lend.gt.0) then
      filsgm=dir_t(1:lend)//sgmnam(1:lens)//'.sgm '
      else
      filsgm=sgmnam(1:lens)//'.sgm '
      endif
c
      pre_sgmfil=.false.
      return
      end

