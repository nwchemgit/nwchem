      subroutine grid_checks(acc_in,iat,alphaa,lla,xa,
     ,     ibt,alphab,llb,xb)
      implicit none
#include "cdft.fh"
#include "mafdecls.fh"
      double precision alphaa,alphab,xa(3),xb(3)
      integer iat,ibt
      integer lla,llb
      double precision acc_in ! [in]
c
      double precision pi
      integer nptstot,mxangpp,mxnangpp
      parameter (mxnangpp=11,mxangpp=5810)
      parameter (nptstot=mxangpp)
c      parameter (nptstot=mxrad*mxangpp)
      integer itype,jtype
      double precision s_err,s_err0
      double precision rij
      double precision range_a,range_b,delta_a,delta_b,
     .     range_a0,range_b0, value,value_ex,maxrad
      integer nrada,nradb,nanga,nangb,nanga0,nangb0,
     ,     deltara,deltarb,dum,niter,mxiter,nrada0,nradb0
      parameter(maxrad=131d0,mxiter=2)
      logical grid_nangl,grid_nradl,grid_incang,grid_incrad
      external grid_nangl,grid_nradl,grid_incang,grid_incrad
      logical lsafe,limprove,lnotacc,lsafea,lsafeb
      integer l_pts,k_pts,l_wgt,k_wgt,k_xyzw,l_xyzw,
     ,     k_iscr,l_iscr,k_rq,l_rq
      double precision reqacc
c
      do_interp=.false.
      rij=sqrt((xa(1)-xb(1))**2+
     .     (xa(2)-xb(2))**2+
     .     (xa(3)-xb(3))**2)
      pi=dacos(-1d0)
      niter=0
c
c     allocate arrays
c
      if (.not.MA_Push_Get(MT_Dbl,mxrad,'pts',l_pts,k_pts))
     &     call errquit('grchk: cannot allocate pts',0)
      if (.not.MA_Push_Get(MT_Dbl,mxrad,'wgts',l_wgt,k_wgt))
     &     call errquit('grchk: cannot allocate wgt',0)
      if (.not.MA_Push_Get(MT_Dbl,4*nptstot,'xyzw',l_xyzw,k_xyzw))
     &     call errquit('grchk: cannot allocate xyzw',0)
      if (.not.MA_Push_Get(MT_Dbl,nptstot*2,'rq',l_rq,k_rq))
     &     call errquit('grchk: cannot allocate rq',0)
      if (.not.MA_Push_Get(MT_int,nptstot,'iscr',l_iscr,k_iscr))
     &     call errquit('grchk: cannot allocate iscr',0)
c
c     preliminary call to lebedev
c      
      iangquad=mxnangpp
      call xc_lebedev(dum,.false.)


      itype=iatype(iat)
      jtype=iatype(ibt)
c
c     radial
c
      lsafe=grid_nradl(itype,mxrad,nrada)
      lsafe=grid_nradl(jtype,mxrad,nradb)
      nrada0=nrada
      nradb0=nradb
c
c     ang
c
      lsafe=grid_nangl(itype,mxangpp,nanga)
      lsafe=grid_nangl(jtype,mxangpp,nangb)
      nanga0=nanga
      nangb0=nangb
c
c     range
c
      range_a=dble(rad_cutoff(1,itype))
      range_b=dble(rad_cutoff(1,jtype))
      range_a0=range_a
      range_b0=range_b
      if(iat.eq.ibt) then
         deltara=1
         deltarb=1
         delta_a=1d0
         delta_b=1d0
      else
         deltara=nrada/6
         deltarb=nradb/6
         delta_a=range_a/6d0
         delta_b=range_b/6d0
      endif
      reqacc=acc_in
      write(0,*) ' required accuracy ',reqacc
      write(0,*) ' alpha ',alphaa,alphab
      write(0,*) ' lla llb rij ',lla,llb,rij
c
c     initial test
c
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,
     ,     iat.eq.ibt,
     ,     value,value_ex,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq))
      s_err = value - value_ex
      lnotacc=abs(s_err).gt.reqacc
      if(.not.lnotacc) goto 2001
c
c     test range
c
      write(0,*) ' 0-ranges ',range_a,range_b
      s_err0=1d99
      s_err=0d0
 12   niter=niter+1
      reqacc=acc_in*mxiter/niter
      range_a=range_a0
      range_b=range_b0
 1688 continue
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,mxangpp,mxangpp,range_a,range_b,
     ,     iat.eq.ibt,
     ,     value,value_ex,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq))
      s_err = value - value_ex
      lnotacc=abs(s_err).gt.reqacc
c      lnotacc=abs(s_err).gt.reqacc*0.05d0
      lsafe=range_a.gt.maxrad.and.range_b.gt.maxrad
      limprove=abs(s_err0-s_err).gt.reqacc*0.1d0
      s_err0=s_err
c      write(0,*) ' rr err ',range_a,range_b,s_err
      if(lnotacc.and.(.not.lsafe).
     .     and.limprove) then
         range_b=range_b+delta_b
         range_a=range_a+delta_a
         goto 1688
      endif
      range_a=min(maxrad,range_a)
      range_b=min(maxrad,range_b)
      rad_cutoff(1,itype)=range_a
      if(iat.ne.ibt) rad_cutoff(1,jtype)=range_b
      write(0,*) ' f-ranges ',range_a,range_b
      if(.not.lnotacc) then
         write(0,*) '@radial extension OK'
      else
         write(0,*) ' cant go beyond acc ',s_err, ' ?? '
         reqacc=(abs(s_err)+abs(reqacc))*0.5d0
         write(0,*) ' new acc ',reqacc
      endif
c
c     test nrad
c
c .99997 for UO2      deltara=nrada/4
c      deltarb=nradb/4
      write(0,*) ' 0-nrad ',nrada,nradb
      s_err0=1d99
      s_err=0d0
      nrada=nrada0
      nradb=nradb0
 2688 continue
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,
     ,     iat.eq.ibt,
     ,     value,value_ex,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq))
      s_err = value - value_ex
      lsafe=nrada.gt.mxrad.and.nradb.gt.mxrad
      lnotacc=abs(s_err).gt.reqacc
      s_err0=s_err
      limprove=(s_err0-s_err).gt.reqacc*0.025d0
      write(0,*) ' nr acc ',nrada,nradb,abs(s_err)
      if(lnotacc.and.(.not.lsafe).and.limprove) then
         lsafea=grid_incrad(mxrad,nrada,deltara)
         lsafeb=grid_incrad(mxrad,nradb,deltarb)
         if(lsafea.and.lsafeb) goto 2689
         goto 2688
      endif
 2689 nrad(itype)=nrada
      nrad(jtype)=nradb
      write(0,*) ' f-nrada nradb',nrada,nradb
      if(.not.lnotacc) write(0,*) '@ nrad  OK '
      if(iat.eq.ibt) goto 2000
c
c     test nang
c
      write(0,*) ' 0-nang ',nanga,nangb
      nanga=nanga0
      nangb=nangb0
      s_err0=1d99
      s_err=0d0
 3688 continue
      call grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,
     ,     iat.eq.ibt,
     ,     value,value_ex,
     ,     dbl_mb(k_pts),dbl_mb(k_wgt),dbl_mb(k_xyzw),
     ,     int_mb(k_iscr),dbl_mb(k_rq))
      s_err = value - value_ex
      lsafe=nanga.gt.mxangpp.and.nanga.gt.mxangpp
      lnotacc=abs(s_err).gt.reqacc
      s_err0=s_err
      if(lnotacc.and.(.not.lsafe)) then
         lsafea=grid_incang(mxnangpp,nangb)
         lsafeb=grid_incang(mxnangpp,nanga)
         if(lsafea.and.lsafeb) goto 3689
         goto 3688
      endif
 3689 nang(itype)=nanga
      nang(jtype)=nangb
      write(0,*) ' f-nang ',nanga,nangb
      if(.not.lnotacc)  write(0,*) '@ nang OK '
 2000 reqacc=acc_in
      lnotacc=abs(s_err).gt.reqacc
      if(lnotacc.and.(niter.le.mxiter)) then
         write(0,*) '#iter',niter,' acc= ',s_err
         s_err0=1d99
         s_err=0d0
         goto 12
      endif
 2001 write(0,*) ' **********'
      write(0,*) ' #iter:',niter,
     ,     ' reached accuracy ? ',abs(s_err).lt.acc_in
      write(0,1) alphaa,alphab,lla,llb,
     ,     '  err_rel',s_err/value_ex, ' err_abs',s_err
      write(0,*) ' **********'
 1    format(' exps',2d10.2,' l ',2i4,A,d10.2,A,d12.4)
      if (.not. ma_chop_stack(l_pts)) call errquit
     $     ('grchks: ma pop?', 0)                 
      return
      end
      subroutine grid_checkscore(alphaa,lla,alphab,llb,rij,
     ,     nrada,nradb,nanga,nangb,range_a,range_b,lsame,
     ,     value,value_ex,
     ,     pts,wghts,xyzw,iscr,rq)
      implicit none
      double precision alphaa,alphab,rij ! [in]
      integer nrada,nradb,nanga,nangb,lla,llb ! [in]
      double precision range_a,range_b ! [in]
      logical lsame
      double precision value,value_ex ! [out]
      double precision pts(*),wghts(*),xyzw(*),rq(*) ! [scratch]
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
     ,        alphab,llb,alphaa,lla,-rij,
     ,        iscr,rq)
         value_ex=grid_intovl(alphaa,lla,alphab,llb,rij)
      endif
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
      integer n0,nangtst ! [in]
      double precision xyzw(4,*),pts(*),wghts(*)
      double precision rij,alphaa,alphab,range_mk 
      double precision rq(*) ! [scratch]
      integer iscr(*)
      double precision p(2),rdum(2)
      integer lla,llb
c
      external xc_leb_nrd2acc
      double precision grid_s2int
      integer nq,j,k,acc_leb,xc_leb_nrd2acc
      double precision xyz(3,2),rrad,tx,ty,tz
      external grid_s2int
c              call eumac(2d0, pts, 
c     &             wghts, n0, 2d0)
      call xc_murakn(pts, 
     &     wghts,n0,range_mk,.true.)
      xyz(1,1)=0d0
      xyz(2,1)=0d0
      xyz(3,1)=0d0
      xyz(1,2)=rij
      xyz(2,2)=0d0
      xyz(3,2)=0d0
      acc_leb=xc_leb_nrd2acc(nangtst)  
      grid_s2int0=0d0
      do j=1,n0
         nq=0
         rrad=pts(j)
         do k = 1, ntheta(acc_leb)
            tx=za_leb(1,k,acc_leb)
            ty=za_leb(2,k,acc_leb)
            tz=za_leb(3,k,acc_leb)
            nq = nq + 1
            xyzw(1,nq) =rrad*tx
            xyzw(2,nq) =rrad*ty
            xyzw(3,nq) =rrad*tz
            xyzw(4,nq) = wa_leb(k,acc_leb)*wghts(j)
         enddo
         call grid_ssw(xyzw, xyz, Rdum, rq, p, 1, 
     &        2,nq,iscr, 'erf1')
      grid_s2int0=grid_s2int0 +
     +        grid_s2int(nq,xyzw,alphaa,lla,alphab,llb,rij)
      enddo
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
      parameter(lscr=10000,lstv=50)
      double precision xa(3),xb(3),dum,S(lstv)
      double precision scr(lscr)
      integer ndum
      xa(1)=0d0
      xa(2)=0d0
      xa(3)=0d0
      xb(1)=rij
      xb(2)=0d0
      xb(3)=0d0
!     call hf1(
!    &     xa,alphaa,
!    &     1d0,
!    &     1, 1, lla,
!    &     xb,alphab,
!    &     1d0,
!    &     1, 1, llb,
!    &     dum,dum,dum,ndum,
!    &     S,dum,dum,lstv,.true.,.false.,.false.,
!    .     .false.,.false.,scr,lscr,'int_1eov')
      call hnd_stvint(
     &     xa,alphaa,
     &     1d0,0
     &     1, 1, lla,
     &     xb,alphab,
     &     1d0,
     &     1, 1, llb,
     &     dum,dum,ndum,
     &     S,dum,dum,lstv,.true.,.false.,.false.,scr,lscr)
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
      nang_out=nang(itype)
      if(nang_out.gt.mxangpp) then
cedo         write(0,*) ' testing angular grid up to ',mxangpp,' pts'
         grid_nangl=.false.
         nang_out=mxangpp
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
         
