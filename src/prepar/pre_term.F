      logical function pre_term(lfnout,lprint,lfnsgm,
     + dir_s,dir_x,dir_u,dir_t,
     + lseq,cseq,mseq,nseq,llnk,clnk,mlnk,nlnk)
c
c     function to determine the terminal segments
c
      implicit none
c
      logical pre_sgmfil,pre_atmscn
      external pre_sgmfil,pre_atmscn
c
      integer lfnout,lfnsgm,mseq,nseq,mlnk,nlnk
      integer lseq(4,mseq),llnk(2,mlnk)
      character*10 cseq(mseq)
      character*4 clnk(2,mlnk)
      character*255 dir_s,dir_x,dir_u,dir_t
      integer lprint
c
      character*255 filsgm
      integer i,isgm,iu,j,jsgm,num,link
c
c     process the explicit link list
c     ------------------------------
c
      do 1 i=1,nlnk
c
      isgm=0
      jsgm=0
      do 2 j=1,nseq
      if(lseq(1,j).eq.llnk(1,i)) isgm=j
      if(lseq(1,j).eq.llnk(2,i)) jsgm=j
    2 continue
      if(isgm.eq.0.or.jsgm.eq.0) call errquit('error in pre_term',9999)
c
      if(.not.pre_sgmfil(cseq(isgm),lfnsgm,filsgm,
     + dir_s,dir_x,dir_u,dir_t))
     + call errquit('pre_sgmfil failed',9999)
c
      if(.not.pre_atmscn(lfnsgm,filsgm,clnk(1,i),num,link))
     + call errquit('pre_atmscn failed',9999)
c
      if(link.eq.1) then
      if(lseq(2,isgm).ge.0)  call errquit('multiple links',9999)
      if(lseq(2,isgm).eq.-1) lseq(2,isgm)=0
      if(lseq(2,isgm).eq.-2) lseq(2,isgm)=1
      endif
      if(link.eq.2) then
      if(lseq(2,isgm).eq.-1) call errquit('multiple links',9999)
      if(lseq(2,isgm).eq.0)  call errquit('multiple links',9999)
      if(lseq(2,isgm).eq.-2) lseq(2,isgm)=-1
      if(lseq(2,isgm).gt.0) lseq(2,isgm)=0
      endif
c
      if(.not.pre_sgmfil(cseq(jsgm),lfnsgm,filsgm,
     + dir_s,dir_x,dir_u,dir_t))
     + call errquit('pre_sgmfil failed',9999)
c
      if(.not.pre_atmscn(lfnsgm,filsgm,clnk(2,i),num,link))
     + call errquit('pre_atmscn failed',9999)
c
      if(link.eq.1) then
      if(lseq(2,jsgm).ge.0)  call errquit('multiple links',9999)
      if(lseq(2,jsgm).eq.-1) lseq(2,jsgm)=0
      if(lseq(2,jsgm).eq.-2) lseq(2,jsgm)=1
      endif
      if(link.eq.2) then
      if(lseq(2,jsgm).eq.-1) call errquit('multiple links',9999)
      if(lseq(2,jsgm).eq.0)  call errquit('multiple links',9999)
      if(lseq(2,jsgm).eq.-2) lseq(2,jsgm)=-1
      if(lseq(2,jsgm).gt.0) lseq(2,jsgm)=0
      endif
c
    1 continue
c
c     adjust the segment names
c     ------------------------
c
      do 3 i=1,nseq+1
      if(lseq(2,i).ne.0) then
      iu=index(cseq(i),' ')
      if(lseq(2,i).eq.-2) cseq(i)(iu:iu+1)='_M'
      if(lseq(2,i).eq.-1) cseq(i)(iu:iu+1)='_N'
      if(lseq(2,i).gt.0) cseq(i)(iu:iu+1)='_C'
      endif
    3 continue
c
      pre_term=.true.
      return
c
 9999 continue
      pre_term=.false.
      return
      end
