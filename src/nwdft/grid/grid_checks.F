      subroutine grid_checks(iat,alphaa,lla,xa,
     ,     ibt,alphab,llb,xb)
      implicit none
#include "cdft.fh"
      double precision alphaa,alphab,xa(3),xb(3)
      integer iat,ibt
      integer lla,llb
c
      double precision pi
      double precision pts(MXrad)
      double precision wghts(MXrad),xyzw(4,MXrad*1202)
      integer itype,n0,L
      double precision s_err,alpha,gg,nn,sfact,s_err0
      double precision grid_sint,range_mk,grid_norm,grid_intovl,
     .     rij,grid_s2int0
      external grid_sint,grid_norm,grid_intovl,grid_s2int0
      double precision range_a,range_b,delta_a,delta_b,
     .     range_a0,range_b0
      double precision reqacc,maxrad
      parameter(reqacc=1d-5,maxrad=91d0)
      logical lsafe,limprove
c
      do_interp=.false.
      rij=sqrt((xa(1)-xb(1))**2+
     .     (xa(2)-xb(2))**2+
     .     (xa(3)-xb(3))**2)
      pi=dacos(-1d0)
c
      itype=iatype(iat)
      sfact=1d0
      s_err0=1d99
      s_err=0d0
      range_a=rad_cutoff(1,itype)
      range_b=rad_cutoff(1,iatype(ibt))
      range_a0=range_a
      range_b0=range_b
      delta_a=0.25*range_a
      delta_b=0.25*range_b
      write(0,*) ' alpha ',alphaa,alphab
      write(0,*) ' lla llb ',lla,llb
 1688 continue
      write(0,*) ' ranges ',range_a,range_b
             if(iat.eq.ibt) then
                n0=nrad(itype)
                range_mk=range_a
                call xc_murakn(pts, 
     &               wghts,n0,range_mk,.true.)
                alpha=alphaa+alphab
                l=lla+llb
                gg=grid_sint(n0,pts,wghts,alpha,l)
                sfact=1d0/acos(-1d0)*.25d0*(lla+llb+1d0)
             else
                n0=nrad(itype)
                range_mk=range_a
                
                gg=grid_s2int0(pts,wghts,n0,range_mk,xyzw,
     ,               nang(itype),
     ,               alphaa,lla,alphab,llb,rij)
                
                itype=iatype(ibt)
                n0=nrad(itype)
                range_mk=range_b
                gg=gg+grid_s2int0(pts,wghts,n0,range_mk,xyzw,
     ,               nang(itype),
     ,               alphab,llb,alphaa,lla,-rij)

             endif
             nn= grid_intovl(alphaa,lla,alphab,llb,rij)
             s_err = gg/sfact - nn
             lsafe=range_a.gt.maxrad.and.range_b.gt.maxrad
             limprove=(s_err0-s_err).gt.reqacc*0.01d0
             s_err0=s_err
             if(abs(s_err).gt.reqacc.and.(.not.lsafe).
     .            and.limprove) then
                write(0,*) ' insuff acc ',s_err
                if(range_a.gt.range_b) then
                   range_b0=range_b
                   range_b=range_b+delta_b
                else
                   range_a0=range_a
                   range_a=range_a+delta_a
                endif
                goto 1688
             endif
             range_a=range_a0
             range_b=range_b0
             rad_cutoff(1,iatype(iat))=range_a
             if(iat.ne.ibt) rad_cutoff(1,iatype(ibt))=range_b
             write(0,*) ' **********'
             write(0,1) alphaa,alphab,l,
     ,            '  err_rel',s_err/nn, ' err_abs',s_err
             write(0,*) ' **********'
 1           format(' exps',2d10.2,' l ',i2,A,d10.2,A,d12.4)
             write(0,*) ' cuts for ',iat,ibt,
     ,            rad_cutoff(1,iatype(ibt)),
     ,            rad_cutoff(1,iatype(iat))
      return
      end
      
      double precision function grid_sint(n0,pts,wghts,alpha,l)
      implicit none
      integer n0,l
      double precision pts(n0)
      double precision wghts(n0)
      double precision alpha

c
      integer i
c
      grid_sint=0d0
c      write(0,*) ' n0 ',n0
      do i=1,n0
         if(l.eq.0) then
         grid_sint=grid_sint+wghts(i)*exp(-alpha*pts(i)*pts(i))
         else
         grid_sint=grid_sint+wghts(i)*exp(-alpha*pts(i)*pts(i))*
     .        pts(i)**l
         endif
cedo      write(0,11) i,wghts(i),pts(i),grid_sint
      enddo
 11   format(' int ',i3,d12.4,F20.16,d12.4)
      return
      end
      double precision function grid_s2int0(
     ,     pts,wghts,n0,range_mk,xyzw,nangtst,
     ,     alphaa,lla,
     .     alphab,llb,rij)
      implicit none
#include "cdft.fh"
      integer n0,nangtst ! [in]
      double precision xyzw(4,*),pts(*),wghts(*)
      double precision rij,alphaa,alphab,range_mk,iscr(MXrad*1202),
     ,     rq(Mxrad*1202,2),p(2),rdum(2)
      integer lla,llb,iag_old,nangpts
c
      double precision grid_s2int
      integer nq,j,k,xc_leb_nrd2acc
      double precision xyz(3,2),rrad,tx,ty,tz
      external grid_s2int,xc_leb_nrd2acc
c              call eumac(2d0, pts, 
c     &             wghts, n0, 2d0)
      call xc_murakn(pts, 
     &     wghts,n0,range_mk,.true.)
      iag_old=iangquad
      iangquad= xc_leb_nrd2acc(nangtst)
      xyz(1,1)=0d0
      xyz(2,1)=0d0
      xyz(3,1)=0d0
      xyz(1,2)=rij
      xyz(2,2)=0d0
      xyz(3,2)=0d0
      call xc_lebedev(nangpts,.false.)
      nq=0
      do j=1,n0
         rrad=pts(j)
            do k = 1, ntheta(iangquad)
               tx=za_leb(1,k,iangquad)
               ty=za_leb(2,k,iangquad)
               tz=za_leb(3,k,iangquad)
               nq = nq + 1
               xyzw(1,nq) =rrad*tx
               xyzw(2,nq) =rrad*ty
               xyzw(3,nq) =rrad*tz
               xyzw(4,nq) = wa_leb(k,iangquad)*wghts(j)
            enddo
         enddo
cedo         write(0,*) ' nq ',nq
      iangquad=iag_old
      if(abs(rij).gt.0.3d0)
     .     call grid_ssw(xyzw, xyz, Rdum, rq, p, 1, 
     &     2,nq,iscr, 'erf1')
      grid_s2int0=grid_s2int(nq,xyzw,alphaa,lla,alphab,llb,rij)
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
      
