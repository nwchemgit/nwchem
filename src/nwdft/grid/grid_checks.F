      subroutine grid_checks(acc_in,iat,alphaa,lla,xa,
     ,     ibt,alphab,llb,xb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
      double precision alphaa,alphab,xa(3),xb(3)
      integer iat,ibt
      integer lla,llb
      double precision acc_in ! [in]
c
      double precision pi
      integer nptstot,mxangpp,mxnangpp,mxangpp1
      parameter (mxnangpp=iac_max,mxangpp=mx_l_ang,mxangpp1=5294)
      parameter (nptstot=mxangpp)
      integer itype,jtype
      double precision s_err
      double precision rij
      double precision range_a,range_b,value,value_ex,maxrad
      integer nrada,nradb,nanga,nangb,me,
     ,     dum,niter,mxiter
      parameter(maxrad=131d0,mxiter=3)
      logical grid_nangl,grid_nradl,same_ctr
      external grid_nangl,grid_nradl
      logical lsafe
      integer l_pts,k_pts,l_wgt,k_wgt,k_xyzw,l_xyzw,
     ,     k_iscr,l_iscr,k_rq,l_rq,iag_old
      double precision reqacc
c
      itype=iatype(iat)
      jtype=iatype(ibt)
      same_ctr=iat.eq.ibt
      if(nrad(itype).eq.0.or.nrad(jtype).eq.0) return
      rij=sqrt((xa(1)-xb(1))*(xa(1)-xb(1))+
     .     (xa(2)-xb(2))*(xa(2)-xb(2))+
     .     (xa(3)-xb(3))*(xa(3)-xb(3)))
      pi=dacos(-1d0)
      niter=0
c
c     allocate arrays
c
      if (.not.MA_Push_Get(MT_Dbl,mxrad,'pts',l_pts,k_pts))
     &     call errquit('grchk: cannot allocate pts',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,mxrad,'wgts',l_wgt,k_wgt))
     &     call errquit('grchk: cannot allocate wgt',1, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,4*nptstot,'xyzw',l_xyzw,k_xyzw))
     &     call errquit('grchk: cannot allocate xyzw',2, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nptstot*2,'rq',l_rq,k_rq))
     &     call errquit('grchk: cannot allocate rq',3, MA_ERR)
      if (.not.MA_Push_Get(MT_int,nptstot,'iscr',l_iscr,k_iscr))
     &     call errquit('grchk: cannot allocate iscr',4, MA_ERR)
c
c     preliminary call to lebedev
c      
      if(leb) then
         iag_old=iangquad
         iangquad=mxnangpp
      else
         iangquad=mxnangpp
      endif
      call xc_lebedev(dum,.false.)
      iangquad=iag_old
      


c
c     radial
c
      lsafe=grid_nradl(itype,mxrad,nrada)
      lsafe=grid_nradl(jtype,mxrad,nradb)
c
c     ang
c
      if(leb) then
         lsafe=grid_nangl(itype,mxangpp,nanga)
         lsafe=grid_nangl(jtype,mxangpp,nangb)
      else
         nanga=1202
         nangb=1202
      endif
c
c     range
c
      range_a=dble(rad_cutoff(1,itype))
      range_b=dble(rad_cutoff(1,jtype))
      reqacc=acc_in
      me=ga_nodeid()
      if(me.eq.0) write(0,*) ' required accuracy ',reqacc
      if(me.eq.0) write(0,*) ' alpha ',alphaa,alphab
      if(me.eq.0) write(0,*) ' lla llb rij ',lla,llb,rij
c
c     initial test
c
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,
     ,     same_ctr,
     ,     value,value_ex,s_err,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq))
      if(s_err.lt.reqacc) goto 2001
      if(.not.same_ctr) then
c
c     first <a|a> at r=0
c

      call grid_checks1(acc_in,alphaa,lla,alphaa,lla,0d0,
     ,     nrada,nrada,nanga,nanga,range_a,range_a,.true.,
     ,     value,value_ex,s_err,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq),
     ,     niter,mxiter,mxangpp,mxangpp1,mxrad,mxnangpp,maxrad)
c
c     2nd <b|b> at r=0
c
      call grid_checks1(acc_in,alphab,llb,alphab,llb,0d0,
     ,     nradb,nradb,nangb,nangb,range_b,range_b,.true.,
     ,     value,value_ex,s_err,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq),
     ,     niter,mxiter,mxangpp,mxangpp1,mxrad,mxnangpp,maxrad)
      endif
c
c     3rd <a|b> at r
c
      call grid_checks1(acc_in,alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,same_ctr,
     ,     value,value_ex,s_err,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq),
     ,     niter,mxiter,mxangpp,mxangpp1,mxrad,mxnangpp,maxrad)

 2001 if(me.eq.0) then
         write(0,*) ' **********'
         write(0,*) ' #iter:',niter,
     ,        ' reached accuracy ? ',abs(s_err).lt.acc_in
         write(0,1) alphaa,alphab,lla,llb,
     ,        '  err_rel',s_err/value_ex, ' err_abs',s_err
         write(0,*) ' **********'
      endif
      rad_cutoff(1,itype)=range_a
      rad_cutoff(1,jtype)=range_b
      nrad(itype)=nrada
      nrad(jtype)=nradb
      if(leb) then
         nang(itype)=nanga
         nang(jtype)=nangb
      endif
      if (.not. ma_chop_stack(l_pts)) call errquit
     $     ('grchks: ma pop?', 0, MA_ERR)
      return
 1    format(' exps',2d10.2,' l ',2i4,A,d10.2,A,d12.4)
      end
      subroutine grid_checks1(acc_in,alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,same_ctr,
     ,     value,value_ex,s_err,
     ,     pts,wghts,xyzw,iscr,rq,
     ,niter,mxiter,mxangpp,mxangpp1,mxrad,mxnangpp,maxrad)
      implicit none
#include "global.fh"
      double precision alphaa,alphab
      integer lla,llb
      double precision acc_in ! [in]
      double precision reqacc
      double precision s_err,s_err0
      double precision rij
      double precision range_a,range_b,delta_a,delta_b,
     .     range_a0,range_b0, value,value_ex,maxrad
      integer nrada,nradb,nanga,nangb,nanga0,nangb0,
     ,     deltara,deltarb,niter,mxiter,nrada0,nradb0
      integer mxangpp,mxangpp1,mxrad,mxnangpp
      double precision pts(*),wghts(*),xyzw(4,*),rq(*) ! [scratch]
      integer iscr(*) ! [scratch]
      integer me
      logical same_ctr,lnotacc,lsafe,limprove,lsafea,lsafeb
      logical grid_nangl,grid_nradl,grid_incang,grid_incrad
      external grid_nangl,grid_nradl,grid_incang,grid_incrad
c
      niter=0
c
c     initial test
c
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,
     ,     same_ctr,
     ,     value,value_ex,s_err,
     ,     pts,wghts,xyzw,
     ,     iscr,rq)
      lnotacc=s_err.gt.acc_in
      if(.not.lnotacc) return
      me = ga_nodeid()

      range_a0=range_a
      range_b0=range_b
      nanga0=nanga
      nangb0=nangb
      nrada0=nrada
      nradb0=nradb
      deltara=nrada/6
      deltarb=nradb/6
      delta_a=range_a/6d0
      delta_b=range_b/6d0
c
c     test range
c
 12   niter=niter+1
      if(me.eq.0) write(0,*) ' 0-ranges ',range_a,range_b
      s_err0=1d99
      s_err=0d0
      reqacc=acc_in*dble(mxiter)/dble(niter)
      range_a=range_a0
      range_b=range_b0
 1688 continue
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,max(mxangpp1,nanga),
     ,     max(mxangpp1,nangb),range_a,range_b,
     ,     same_ctr,
     ,     value,value_ex,s_err,
     ,     pts,wghts,xyzw,
     ,     iscr,rq)
      if(same_ctr) then
         lnotacc=s_err.gt.(reqacc*0.025d0)
      else
         lnotacc=s_err.gt.(reqacc*0.025d0)
      endif
      lsafe=range_a.gt.maxrad.and.range_b.gt.maxrad
      limprove=(s_err0-s_err).gt.reqacc*0.001d0
      if(lnotacc.and.(.not.lsafe).
     .     and.limprove) then
         range_b=min(range_b+delta_b,maxrad)
         range_a=min(range_a+delta_a,maxrad)
         s_err0=s_err
         goto 1688
      endif
      range_a=min(maxrad,range_a)
      range_b=min(maxrad,range_b)
      if(me.eq.0) write(0,*) ' f-ranges ',range_a,range_b
      if(.not.lnotacc) then
         if(me.eq.0) write(0,*) '@radial extension OK'
      else
cedo         if(.not.limprove) then
cedo            range_a=range_a-delta_a
cedo            range_b=range_b-delta_b
cedo         endif
         reqacc=(s_err0+reqacc)*0.5d0
      endif
c
c     test nrad
c
      if(me.eq.0) write(0,*) ' 0-nrad ',nrada,nradb
      s_err0=1d99
      s_err=0d0
      nrada=nrada0
      nradb=nradb0
 2688 continue
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,
     ,     same_ctr,
     ,     value,value_ex,s_err,
     ,     pts,wghts,xyzw,
     ,     iscr,rq)
      lsafe=nrada.gt.mxrad.and.nradb.gt.mxrad
      lnotacc=s_err.gt.reqacc
      limprove=(s_err0-s_err).gt.reqacc*0.01d0
      s_err0=s_err
      if(lnotacc.and.(.not.lsafe).and.limprove) then
         lsafea=grid_incrad(mxrad,nrada,deltara)
         lsafeb=grid_incrad(mxrad,nradb,deltarb)
         if(lsafea.and.lsafeb) goto 2689
         goto 2688
      endif
 2689 continue
      if(me.eq.0) write(0,*) ' f-nrad ',nrada,nradb
      if(.not.lnotacc.and.me.eq.0) write(0,*) '@ nrad  OK '
      if(same_ctr) goto 2000
c
c     test nang
c
      nanga=nanga0
      nangb=nangb0
      if(me.eq.0) write(0,*) ' 0-nang ',nanga,nangb
      s_err0=1d99
      s_err=0d0
 3688 continue
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,
     ,     same_ctr,
     ,     value,value_ex,s_err,
     ,     pts,wghts,xyzw,
     ,     iscr,rq)
      lsafe=nanga.gt.mxangpp.and.nanga.gt.mxangpp
      lnotacc=s_err.gt.reqacc
      s_err0=s_err
      if(lnotacc.and.(.not.lsafe)) then
         lsafea=grid_incang(mxnangpp,nangb)
         lsafeb=grid_incang(mxnangpp,nanga)
         if(lsafea.and.lsafeb) goto 3689
         goto 3688
      endif
 3689 continue
      if(me.eq.0) write(0,*) ' f-nang ',nanga,nangb
      if(.not.lnotacc.and.me.eq.0)  write(0,*) '@ nang OK '
 2000 reqacc=acc_in
      lnotacc=s_err.gt.reqacc
      if(lnotacc.and.(niter.le.mxiter)) then
         if(me.eq.0) write(0,*) '#iter',niter,' acc= ',s_err
         s_err0=1d99
         s_err=0d0
         goto 12
      endif
      return
      end
      subroutine grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,lsame,
     ,     value,value_ex,s_err,
     ,     pts,wghts,xyzw,iscr,rq)
      implicit none
      double precision alphaa,alphab,rij ! [in]
      integer nrada,nradb,nanga,nangb,lla,llb ! [in]
      double precision range_a,range_b ! [in]
      logical lsame
      double precision value,value_ex,s_err ! [out]
      double precision pts(*),wghts(*),
     ,     xyzw(4,*),rq(*)      ! [scratch]
      integer iscr(*) ! [scratch]
c
      double precision range_mk
      double precision grid_sint,grid_norm,grid_intovl,
     .     grid_s2int0
      external grid_sint,grid_norm,grid_intovl,
     .     grid_s2int0
c

      if(lsame) then
         range_mk=range_a
         call xc_murakn(pts, 
     &        wghts,nrada,range_mk,.true.)
         value=grid_sint(nrada,pts,wghts,
     ,        alphaa+alphab,lla+llb)
         value_ex= grid_norm(lla+llb,alphaa+alphab)
      else
         range_mk=range_a
         value=grid_s2int0(pts,wghts,nrada,nanga,
     .        range_mk,xyzw,
     ,        alphaa,lla,alphab,llb,rij,
     ,        iscr,rq)
c     
         range_mk=range_b
         value=value+grid_s2int0(pts,wghts,nradb,nangb,
     ,        range_mk,xyzw,
     ,        alphab,llb,alphaa,lla,(-rij),
     ,        iscr,rq)
         value_ex=grid_intovl(alphaa,lla,alphab,llb,rij)
      endif
      s_err=abs(value-value_ex)
      return
      end
      double precision function grid_sint(n0,pts,wghts,alpha,l)
      implicit none
      integer n0,l
      double precision pts(*)
      double precision wghts(*)
      double precision alpha

c
      integer i
c
      grid_sint=0d0
      do i=1,n0
         if(l.eq.0) then
         grid_sint=grid_sint+wghts(i)*exp(-alpha*pts(i)*pts(i))
         else
         grid_sint=grid_sint+wghts(i)*exp(-alpha*pts(i)*pts(i))*
     .        pts(i)**l
         endif
!      write(0,11) i,wghts(i),pts(i),grid_sint
      enddo
 11   format(' int ',i3,d12.4,F20.16,d12.4)
      return
      end
      double precision function grid_s2int0(
     ,     pts,wghts,n0,nangtst,range_mk,xyzw,
     ,     alphaa,lla,
     .     alphab,llb,rij,iscr,rq)
      implicit none
#include "cdft.fh"
#include "global.fh"
#include "msgids.fh"
      integer n0,nangtst ! [in]
      double precision xyzw(4,nangtst),pts(n0),wghts(n0)
      double precision rij,alphaa,alphab,range_mk 
      double precision rq(2,nangtst) ! [scratch]
      integer iscr(nangtst)
      double precision p(2),rdum(3)
      integer lla,llb
c
      external xc_leb_nrd2acc
      double precision grid_s2int
      integer nq,j,k,acc_leb,xc_leb_nrd2acc
      double precision xyz(3,2),rrad,ww
      integer me,nproc
      external grid_s2int
c              call eumac(2d0, pts, 
c     &             wghts, n0, 2d0)
      call xc_murakn(pts, 
     &     wghts,n0,range_mk,.true.)
      me = ga_nodeid()
      nproc = ga_nnodes()
      xyz(1,1)=0d0
      xyz(2,1)=0d0
      xyz(3,1)=0d0
      xyz(1,2)=rij
      xyz(2,2)=0d0
      xyz(3,2)=0d0
      acc_leb=xc_leb_nrd2acc(nangtst)  
      grid_s2int0=0d0
      do j=1+me,n0,nproc
         nq=0
         rrad=pts(j)
         ww=wghts(j)
         do k = 1, nangtst
            nq = nq + 1
            xyzw(1,nq) =rrad*za_leb(1,k,acc_leb)
            xyzw(2,nq) =rrad*za_leb(2,k,acc_leb)
            xyzw(3,nq) =rrad*za_leb(3,k,acc_leb)
            xyzw(4,nq) = wa_leb(k,acc_leb)*ww
         enddo
         call grid_ssw(xyzw, xyz, Rdum, rq, p, 1, 
     &        2,nq,iscr, 0, 'erf1')
      grid_s2int0=grid_s2int0 +
     +        grid_s2int(nq,xyzw,alphaa,lla,alphab,llb,rij)
      enddo
      call ga_dgop(msg_igop_tru1,grid_s2int0,1,'+') 
      return
      end
                
      double precision function grid_s2int(n0,xyzw,alphaa,lla,
     .     alphab,llb,rij)
      implicit none
      integer n0,lla,llb
      double precision alphaa,alphab
      double precision xyzw(4,n0)
      double precision rij
c
      integer i
      double precision xx,xr,x,xmr
c
      grid_s2int=0d0
cedo      write(0,*) ' n0 ',n0, ' rij ',rij,' alpha ',alphaa,alphab,
cedo     , ' ll ',lla,llb
      do i=1,n0
         x=xyzw(1,i)
         xmr=xyzw(1,i)-rij
         xx=xyzw(1,i)**2+xyzw(2,i)**2+xyzw(3,i)**2
         xr=(xyzw(1,i)-rij)**2+xyzw(2,i)**2+xyzw(3,i)**2
         grid_s2int=grid_s2int+xyzw(4,i)*(exp(-alphaa*xx)*x**lla*
     *        exp(-alphab*xr)*xmr**llb)
cedo      write(0,11) i,xyzw(4,i),x,grid_s2int
cedo 11   format(' int2 ',i3,d12.4,f20.16,d12.4)
      enddo
      return
      end
      
      double precision function grid_intovl(alphaa,lla,
     .     alphab,llb,rij)
      integer lla,llb
      double precision alphaa,alphab,rij
      integer lscr,lstv
      parameter(lscr=20000,lstv=100)
      double precision xa(3),xb(3),dum,S(lstv)
      double precision scr(lscr)
      integer ndum
      xa(1)=0d0
      xa(2)=0d0
      xa(3)=0d0
      xb(1)=rij
      xb(2)=0d0
      xb(3)=0d0
      call hf1(
     &     xa,alphaa,
     &     1d0,
     &     1, 1, lla,
     &     xb,alphab,
     &     1d0,
     &     1, 1, llb,
     &     dum,dum,dum,ndum,
     &     S,dum,dum,lstv,.true.,.false.,.false.,
     .     .false.,.false.,scr,lscr)
cedo      call hnd_stvint(
cedo     &     xa,alphaa,
cedo     &     1d0,0
cedo     &     1, 1, lla,
cedo     &     xb,alphab,
cedo     &     1d0,
cedo     &     1, 1, llb,
cedo     &     dum,dum,ndum,
cedo     &     S,dum,dum,lstv,.true.,.false.,.false.,scr,lscr)
      grid_intovl=S(1)
c      write(0,'(A,d12.4,/" X",6F10.6)') ' hnd S',grid_intovl,xa,xb
      return
      end
      double precision function grid_norm(l,alpha)
      implicit none
      integer l
      double precision alpha
c
      logical even
      integer nff,n,n_in
      double precision sqrtpi
      integer grid_fact,grid_ffact
      external grid_fact,grid_ffact
c
      sqrtpi=sqrt(dacos(-1d0))
      even=mod(l,2).eq.0
      if(even) then
         n=(l+2)/2
c     
c     compute 2n-1 factorial ... 
c     
         n_in=2*n-1
         nff=grid_ffact(n_in)
cedo         write(0,*) ' 2n-1',2*n-1
cedo         write(0,*) l,' even fact ',nff
         grid_norm=sqrtpi*alpha**(-n-0.5d0)*nff*2d0**(-n-1)
      else
c     
c     compute n factorial ... 
c     
         n=(l+1)/2
         nff=grid_fact(n)
cedo         write(0,*) l,' odd fact ',nff
         grid_norm=alpha**(-n-1)*nff*0.5d0
      endif
      return
      end
      integer function grid_fact(n)
      integer n
      integer i,pp
      pp=1
      if(n.gt.1) then
         pp=1
         do i=2,n
            pp=i*pp
         enddo
      endif
      grid_fact=pp
cedo      write(0,*) ' fl ',grid_fact
      return
      end
      integer function grid_ffact(n)
      integer n
      integer i,pp
      pp=1
cedo      write(0,*) ' ff: n ',n
      if(n.gt.1) then
         pp=1
         do i=3,n,2
            pp=i*pp
         enddo
      endif
      grid_ffact=pp
cedo      write(0,*) ' ffl ',grid_ffact
      return
      end
      logical function grid_nangl(itype,mxangpp,nang_out)
      implicit none
#include "cdft.fh"     
      integer itype,nang_out
      integer mxangpp
c     return false if we exceeded mxang
      grid_nangl=.true.
      if(leb) then
         nang_out=nang(itype)
         if(nang_out.gt.mxangpp) then
c     edo         write(0,*) ' testing angular grid up to ',mxangpp,' pts'
            grid_nangl=.false.
            nang_out=mxangpp
         endif
      else
         nang_out=1202
      endif
      return
      end
      logical function grid_incang(mxnangpp,nang_out)
      implicit none
#include "cdft.fh"     
      integer nang_out,nang_in
      integer mxnangpp,acc
      integer xc_leb_nrd2acc,xc_leb_acc2nrd
      external  xc_leb_nrd2acc,xc_leb_acc2nrd
c     return false if we exceeded mxang
      grid_incang=.false.
      nang_in=nang_out
      acc=xc_leb_nrd2acc(nang_out)
      acc=acc+1
      if(acc.gt.mxnangpp) then
cedo         write(0,*) ' testing angular grid up to level ',mxangpp
         grid_incang=.true.
         nang_out=nang_in
         return
      endif
      nang_out=xc_leb_acc2nrd(acc)
      return
      end
      logical function grid_incrad(mxrad,nrad_out,deltar)
      implicit none
      integer nrad_out,nrad_in,deltar
      integer mxrad
c     return false if we exceeded mxrad
      grid_incrad=.false.
      nrad_in=nrad_out
      nrad_out=nrad_out+deltar
      if(nrad_out.gt.mxrad) then
cedo         write(0,*) ' testing radial grid up to ',mxrad, ' pts'
         grid_incrad=.true.
         nrad_out=mxrad
      endif
      return
      end
      logical function grid_nradl(itype,mxradpp,nrad_out)
      implicit none
#include "cdft.fh"     
      integer itype,nrad_out
      integer mxradpp
c     return false if we exceeded mxang
      grid_nradl=.true.
      nrad_out=nrad(itype)
      if(nrad_out.gt.mxradpp) then
cedo         write(0,*) ' testing radial grid up to ',mxradpp
         grid_nradl=.false.
         nrad_out=mxradpp
      endif
      return
      end
         
