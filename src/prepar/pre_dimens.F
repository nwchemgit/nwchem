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
      character*80 dir_s,dir_x,dir_u,dir_t
      integer lprint
c
      character*80 filsgm
      integer i,nsatm,nsbnd,nsang,nsdih,nsimp
c
      natm=0
      nbnd=0
      nang=0
      ndih=0
      nimp=0
c
c     loop over sequence list
c     -----------------------
c
      do 1 i=1,nseq
c
c     find segment file
c     -----------------
c
      if(.not.pre_sgmfil(cseq(i),lfnsgm,filsgm,
     + dir_s,dir_x,dir_u,dir_t))
     + call errquit('pre_sgmfil failed',9999)
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
    1 continue
c
      pre_dimens=.true.
      return
c
 9999 continue
      pre_dimens=.false.
      return
      end
