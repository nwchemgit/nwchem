      logical function pre_dimens(lfnout,lprint,lfnsgm,
     + dir_s,dir_x,dir_u,dir_t,
     + lseq,cseq,mseq,nseq,natm,nbnd,nang,ndih,nimp)
c
c     function to determine total system array dimensions
c
      implicit none
c
      logical pre_sgmfil,pre_sgmsiz
      external pre_sgmfil,pre_sgmsiz
c
      integer lfnout,lfnsgm,mseq,nseq,natm,nbnd,nang,ndih,nimp
      integer lseq(4,mseq)
      character*10 cseq(mseq)
      character*255 dir_s,dir_x,dir_u,dir_t
      integer lprint
c
      character*255 filsgm
      integer i,nsatm,nsbnd,nsang,nsdih,nsimp
c
      natm=0
      nbnd=0
      nang=0
      ndih=0
      nimp=0
c
      pre_dimens=.true.
c
c     loop over sequence list
c     -----------------------
c
      do 1 i=1,nseq+1
      if(i.lt.nseq+1.or.cseq(i)(1:1).ne.' ') then
c
      if(lprint.ge.5) write(lfnout,1000) cseq(i)
 1000 format('SEARCHING FOR SEGMENT FILE FOR ',a10)
c
c     find segment file
c     -----------------
c
      if(pre_sgmfil(cseq(i),lfnsgm,filsgm,
     + dir_s,dir_x,dir_u,dir_t)) then
c
      if(lprint.ge.5) write(lfnout,1001) filsgm(1:index(filsgm,' ')-1)
 1001 format('FOUND SEGMENT FILE ',a)
c
c     find segment dimensions
c     -----------------------
c
      if(.not.pre_sgmsiz(lfnsgm,filsgm,nsatm,nsbnd,nsang,nsdih,nsimp))
     + call errquit('pre_sgmsiz failed',9999)
c
      natm=natm+nsatm
      nbnd=nbnd+nsbnd
      nang=nang+nsang
      ndih=ndih+nsdih
      nimp=nimp+nsimp
c
      else
      if(lprint.ge.1) write(lfnout,1002) cseq(i)
 1002 format('A SEGMENT FILE COULD NOT BE FOUND FOR ',a10)
      pre_dimens=.false.
      endif
      endif
    1 continue
c
      return
      end
