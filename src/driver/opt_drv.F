      logical function drv_opt(rtdb)
C$Id: opt_drv.F,v 1.8 1996-12-04 23:14:23 mdupuis Exp $
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
      integer  rtdb
      logical  hnd_optx
      external hnd_optx
      integer  ir, iw
      logical  some, dbug
      common/hnd_iofile/ir,iw
c
      dbug=.false.
      some=.false.
      some=some.or.dbug
      ir=LuIn
      iw=LuOut  
      if(some) then
         write(iw,*) ' drv_opt: calling hnd_optx . '
      endif
c
      drv_opt = hnd_optx(rtdb)
c
      end
      logical function hnd_optx(rtdb)
      implicit double precision (a-h,o-z)
#include "rtdb.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer  rtdb
      integer  geom
      logical  status
      logical  opt_input
      logical  some
      logical  dbug
      logical  cvged
      logical  zcoord
      logical  rstart
      logical  excess
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c                mxcoor=max(mxcart,mxzmat)
      common/hnd_iofile/ir,iw
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      common/hnd_optrun/npts,nserch,nupdat,cvged
      common/hnd_optopt/cvgopt,nptopt,inhess,linopt,modupd,moddir
      common/hnd_optfrz/icfrz(mxcoor),iatfrz(mxatom)
      common/hnd_optfun/e,eg(mxcart)
      common/hnd_optvar/zcoord,ncoord,mcoord
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(mxcart),zan(mxatom),nat
      dimension t0(mxcoor)
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data alphmn /0.4d+00/
      data alphmx /2.5d+00/
c
      dbug=.false.  
      some=.false.
      some=some.or.dbug
      if (ga_nodeid().eq.0.or.some)  then
         write(iw,9999)
      endif
c
c     ----- parameters for optimization driver -----
c           may have been read in already. if
c           not, use default.
c
      status = rtdb_parallel(.false.)
      if (ga_nodeid().eq.0)  then
         status=rtdb_get(rtdb,'driver:opt',mt_log,1,opt_input)
         if(.not.status) then
            cvgopt=8.0d-04
            nptopt=20
            linopt=10
            inhess=0
            modupd=1
            do i=1,mxcoor
               icfrz(i)=0
            enddo
            do i=1,mxatom
               iatfrz(i)=0
            enddo
         endif
         call ga_brdcst( 1+msgdbl,cvgopt,mdtob(1     ),0)
         call ga_brdcst( 2+msgint,nptopt,mitob(1     ),0)
         call ga_brdcst( 3+msgint,linopt,mitob(1     ),0)
         call ga_brdcst( 4+msgint,inhess,mitob(1     ),0)
         call ga_brdcst( 5+msgint,modupd,mitob(1     ),0)
         call ga_brdcst( 6+msgint,icfrz ,mitob(mxcoor),0)
         call ga_brdcst( 7+msgint,iatfrz,mitob(mxatom),0)
      else
         call ga_brdcst( 1+msgdbl,cvgopt,mdtob(1     ),0)
         call ga_brdcst( 2+msgint,nptopt,mitob(1     ),0)
         call ga_brdcst( 3+msgint,linopt,mitob(1     ),0)
         call ga_brdcst( 4+msgint,inhess,mitob(1     ),0)
         call ga_brdcst( 5+msgint,modupd,mitob(1     ),0)
         call ga_brdcst( 6+msgint,icfrz ,mitob(mxcoor),0)
         call ga_brdcst( 7+msgint,iatfrz,mitob(mxatom),0)
      endif
      call ga_sync()
      status = rtdb_parallel(.true.) 
c
      if (ga_nodeid().eq.0.or.dbug)  then
         write(iw,9996) cvgopt,nptopt,linopt,inhess,modupd
      endif
c
c     ----- get going ... -----
c
      if (ga_nodeid().eq.0)  then
         call ga_brdcst( 8+msgint,nzmat ,mitob(1     ),0)
         call ga_brdcst( 9+msgint,nzvar ,mitob(1     ),0)
         call ga_brdcst(10+msgint,nvar  ,mitob(1     ),0)
         call ga_brdcst(11+msgint,nat   ,mitob(1     ),0)
         call ga_brdcst(12+msgdbl,c     ,mdtob(mxcart),0)
         call ga_brdcst(13+msgdbl,zan   ,mdtob(mxatom),0)
      else
         call ga_brdcst( 8+msgint,nzmat ,mitob(1     ),0)
         call ga_brdcst( 9+msgint,nzvar ,mitob(1     ),0)
         call ga_brdcst(10+msgint,nvar  ,mitob(1     ),0)
         call ga_brdcst(11+msgint,nat   ,mitob(1     ),0)
         call ga_brdcst(12+msgdbl,c     ,mdtob(mxcart),0)
         call ga_brdcst(13+msgdbl,zan   ,mdtob(mxatom),0)
      endif
c
      rstart=.false.
      excess=.false.
      zcoord=.false.
      zcoord=nzmat.gt.0
      if(zcoord) then
         call geom_zmt_icoord
         ncart =3*nat
         ncoord=nvar                  
         mcoord=nzvar
      else
         nzmat =0
         nzvar =0
         nvar  =3*nat
         ncart =3*nat
         ncoord=ncart
         mcoord=ncart
      endif
      if(dbug)  then
         write(iw,9997) nzvar,nvar,ncart,ncoord,mcoord
      endif
c
c
c     ----- fresh start : initialize hessian matrix -----
c
      call hnd_opt_hss_init(rtdb)
c
c     ----- set up ----  
c
      e     =zero
      func  =zero
      func0 =zero
      gmax  =zero
      gmax0 =zero
      curv  =zero
      do i=1,mxcoor
         dx(i)=zero
         ds(i)=zero
         eg(i)=zero
         g0(i)=zero
         x0(i)=zero
         t0(i)=zero
          x(i)=zero
          g(i)=zero
      enddo
      alpha =zero
      npts  =0
      nserch=0
      nupdat=0
      cvged =.false.
c
c     ----- first point -----  
c
      npts=1
      do i=1,ncart
         x(i)=c(i)
      enddo
c
      call hnd_opt_energy(rtdb)
c
      func=e
c
      call hnd_opt_gradient(rtdb)
c
      do i=1,ncart
         g(i)=eg(i)
      enddo
      if(zcoord) then
         call hnd_opt_tfgx(rtdb)
      endif
      call hnd_opt_gmax(rtdb)
      call hnd_opt_converge(rtdb,cvged)
      call hnd_opt_print(rtdb)
      do i=1,ncart
         x0(i)=x(i)
      enddo
      do i=1,ncoord
         g0(i)=g(i)
         t0(i)=g(i)
      enddo
c
c     ----- exit if converged -----                     
c
      if(cvged) then
         if (ga_nodeid().eq.0.or.dbug)  then
            write(iw,9998)
         endif
         hnd_optx=.true.
         return
      endif
c
  100 continue
c
c     ----- search direction -----
c
      nserch=nserch+1
      call hnd_opt_search_dir(rtdb)
      do i=1,ncoord
         ds(i)=dx(i)
      enddo
      if(zcoord) then
         call hnd_opt_tfds(rtdb)               
      endif
c
c     ----- linear search -----
c
      if(alpha.eq.zero) then
         if(zcoord) then
            alpha=one
         else
            alpha=pt5
         endif
      endif
      call hnd_opt_search_lin(rtdb,geom,rstart,excess)
      if(excess) then
         if (ga_nodeid().eq.0.or.dbug)  then
            write(iw,9995)
         endif
         return
      endif
      do i=1,ncoord
         dx(i)=dx(i)*alpha
         ds(i)=ds(i)*alpha
      enddo
c
c     ----- improved point found -----
c
      if(zcoord) then
         call hnd_opt_tfgx(rtdb)
      endif
      nupdat=nupdat+1
      call hnd_opt_gmax(rtdb)
      call hnd_opt_converge(rtdb,cvged)
      call hnd_opt_print(rtdb)
      if(cvged) then
         if (ga_nodeid().eq.0.or.dbug)  then
            write(iw,9998)
         endif
         hnd_optx=.true.
         return
      endif
c
c     ----- update hessian matrix -----
c
      do i=1,ncoord
         g0(i)=t0(i)
         dx(i)=ds(i)
      enddo
      call hnd_opt_hss_update(rtdb)
c
      alpha=sqrt(abs(alpha))
      if(alpha.gt.alphmx) alpha=alphmx
      if(alpha.lt.alphmn) alpha=alphmn
c
c     ----- reset pt -0- -----
c
      do i=1,ncart
         x0(i)=x(i)
      enddo
      do i=1,ncoord
         t0(i)=g(i)
      enddo
c
c     ----- go back to next search direction -----
c
      go to 100
c
      hnd_optx=.false.
      return
 9999 format(/,10x,21(1h-),
     1       /,10x,'Geometry Optimization',
     2       /,10x,21(1h-))
 9998 format(/,10x,22(1h-),
     1       /,10x,'Optimization converged',
     2       /,10x,22(1h-))
 9997 format(' in hnd_optx, nzvar,nvar,ncart,ncoord,mcoord = ',5i5)
 9996 format(' cvgopt = ',f10.6,/,' nptopt = ',i4,/,' linopt = ',i4,
     1                          /,' inhess = ',i4,/,' modupd = ',i4,/)
 9995 format(' too many points in geometry optimization. stop')
      end
      subroutine optdrv_input(rtdb)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer  rtdb
      logical  status
      logical  opt_input
c
      call hnd_opt_input
c
      opt_input=.true.
      status=rtdb_put(rtdb,'driver:opt',mt_log,1,opt_input)
      if(.not.status) then
         call errquit('optdrv_input : rtdb_put failed',911)
      endif
c
      return
      end
      subroutine hnd_opt_input
      implicit double precision (a-h,o-z)
#include "inp.fh"
c
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c                mxcoor=max(mxcart,mxzmat)
      character*255 field
      logical       status
      logical       ready
      logical       dbug
      common/hnd_iofile/ir,iw
      common/hnd_optopt/cvgopt,nptopt,inhess,linopt,modupd,moddir
      common/hnd_optfrz/icfrz(mxcoor),iatfrz(mxatom)
c
      dbug=.false.
c
c     ----- default parameters for optimization driver -----
c
      cvgopt=8.0d-04
      nptopt=20
      linopt=10
      inhess=0
      modupd=1
      do i=1,mxcoor
         icfrz(i)=0
      enddo
      do i=1,mxatom
         iatfrz(i)=0
      enddo
c
c     ----- read values from input -----
c
      ifield=inp_cur_field()
      if(dbug) then
         write(iw,*) 'in hnd_opt_input, ifield = ',ifield
      endif
      ready=.false.
c
   10 continue
      if(.not.inp_a(field)) then
c
         if(.not.inp_read()) then
            write(iw,*)  'unexpected eof while reading ',
     1                   'hnd_opt_input data'
            call errquit('hnd_opt_input: unexpected eof',911)
         else
            if(dbug) then
               write(iw,*) 'new line read in hnd_opt_input'
            endif
            go to 10
         endif
c
      else       
         if(dbug) then
            write(iw,*) 'in hnd_opt_input, field = ',field
         endif
         if(.not.ready) then
            if(inp_compare(.false.,'opt',field)) then
               ready=.true.
               go to 10
            else
               write(iw,*) 'un-expected value of field'
               call errquit('hnd_opt_input: misposition',911)
            endif
         endif
c
         if(inp_compare(.false.,'end',field)) then
            go to 100
         else
            if(inp_compare(.false.,'cvgopt',field)) then
               if(inp_f(dum)) then
                  cvgopt=dum
                  go to 10
               else
                  call errquit('hnd_opt_input: cvgopt',911)
               endif
            elseif(inp_compare(.false.,'nptopt',field)) then
               if(inp_i(idum)) then
                  nptopt=idum
                  go to 10
               else
                  call errquit('hnd_opt_input: nptopt',911)
               endif
            elseif(inp_compare(.false.,'linopt',field)) then
               if(inp_i(idum)) then
                  linopt=idum
                  go to 10
               else
                  call errquit('hnd_opt_input: linopt',911)
               endif
            elseif(inp_compare(.false.,'inhess',field)) then
               if(inp_i(idum)) then
                  inhess=idum
                  go to 10
               else
                  call errquit('hnd_opt_input: inhess',911)
               endif
            elseif(inp_compare(.false.,'modupd',field)) then
               if(inp_i(idum)) then
                  modupd=idum
                  go to 10
               else
                  call errquit('hnd_opt_input: modupd',911)
               endif
            endif
         endif
c
      endif
c
  100 continue
      if(dbug) then
         write(iw,*) 'cvgopt',cvgopt
         write(iw,*) 'nptopt',nptopt
         write(iw,*) 'linopt',linopt
         write(iw,*) 'inhess',inhess
         write(iw,*) 'modupd',modupd
      endif
c
      return
      end
      subroutine hnd_opt_gmax(rtdb)
      implicit double precision (a-h,o-z)
      integer  rtdb
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      logical    zcoord
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      common/hnd_optopt/cvgopt,nptopt,inhess,linopt,modupd,moddir
      common/hnd_optfrz/icfrz(mxcoor),iatfrz(mxatom)
      common/hnd_optvar/zcoord,ncoord,mcoord
      common/hnd_zmtpar/nzmat,nzvar,nvar
      data zero /0.0d+00/
c
      gmax =zero
      do i=1,ncoord 
         if(icfrz(i).eq.0) then
            if(abs(g(i)).gt.gmax) gmax=abs(g(i))
         endif
      enddo
c
      gnorm=zero
      do i=1,ncoord
         gnorm=gnorm+g(i)*g(i)
      enddo
      gnorm=sqrt(gnorm)
c
      return
      end
      subroutine hnd_opt_converge(rtdb,cvged)
      implicit double precision (a-h,o-z)
#include "rtdb.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer    rtdb
      logical    cvged
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      common/hnd_iofile/ir,iw
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      common/hnd_optopt/cvgopt,nptopt,inhess,linopt,modupd,moddir
c
      if( ga_nodeid().eq.0) then
         write(iw,9999) gmax,cvgopt,gmax0
      endif
c
      cvged=gmax.lt.cvgopt
c
      return
 9999 format('          largest component of the gradient = ',F14.7,
     1       ' ( cvgopt = ',F14.7,' ) ',/,
     2       ' previous largest component of the gradient = ',F14.7)
      end
      subroutine hnd_opt_energy(rtdb)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      logical  task_energy
      external task_energy
      logical  status
      logical  cvged
      logical  dbug
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      common/hnd_iofile/ir,iw
      common/hnd_optrun/npts,nserch,nupdat,cvged
      common/hnd_optfun/e,eg(mxcart)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      dbug=.false.
c
      if( ga_nodeid().eq.0.or.dbug) then
         write(iw,9999) nserch,npts
         write(iw,9996)
         do iat=1,nat
            write(iw,9995) iat,(c(i,iat),i=1,3)
         enddo
         write(iw,9997)
      endif
c
      if(nzmat.gt.0) then
         call hnd_bandbi(nat,c)
         call hnd_zmtprt
         write(iw,9997)
      endif
c
      if (task_energy(rtdb)) then
         status=rtdb_get(rtdb,'task:energy',MT_DBL,1,e)
         if(dbug) then
            write(iw,9998) e
         endif
      else
        write(iw,*) 'scf_optimize: scf failed'
      endif
c
      return
 9999 format(1h1,' nserch',i3,//,'   point',I3)
 9998 format(' in opt_energy, e = ',f15.10)
 9997 format(/)
 9996 format(21x,21(1h-),/,
     1       21x,'cartesian coordinates',/,
     2       21x,21(1h-))
 9995 format(9x,i5,3f15.8) 
      end
      subroutine hnd_opt_gradient(rtdb)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      logical  task_gradient
      external task_gradient
      logical  status
      logical  dbug
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      common/hnd_iofile/ir,iw
      common/hnd_optfun/e,eg(mxcart)
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      data tol  /1.0d-08/
      data zero /0.0d+00/
c
      dbug=.false.
c
      ncart=3*nat
      if (task_gradient(rtdb)) then
         status=rtdb_get(rtdb,'task:gradient',MT_DBL,ncart,eg)
         do i=1,3*nat
            if(abs(eg(i)).lt.tol) eg(i)=zero
         enddo
         if(dbug) then
            write(iw,9999)
            write(iw,9998) (eg(i),i=1,ncart)
         endif
      else
         write(iw,*) 'scf_optimize: gradients failed'
      endif
c
      return
 9999 format(' in opt_gradient, eg = ')
 9998 format(3f12.6)
      end
      subroutine hnd_opt_tfgx(rtdb)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      logical  gradients
      external gradients
      logical  geom_zmt_get_binv
      logical  geom_hnd_tfgx
      logical  status
      logical  dbug
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      common/hnd_iofile/ir,iw
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      dimension gs(mxcoor)
      data zero /0.0d+00/
c
      dbug=.false.
c
c     ----- get -b- and -bi- -----
c
      call hnd_bandbi(nat,x)
c
c     ----- transform -gx- to -gs- -----
c
      ncart=3*nat
      status=ma_push_get(MT_DBL,ncart*nzvar,'mem bi',i_ibi,ibi)
      status=geom_zmt_get_binv(dbl_mb(ibi),ncart,nzvar,ncart)
      status=geom_hnd_tfgx(gs,g,dbl_mb(ibi),nvar,ncart,ncart)
      status=ma_pop_stack(i_ibi)
      do i=1,nvar
         g(i)=gs(i)
      enddo
      do i=nvar+1,ncart
         g(i)=zero
      enddo
c
      return
      end
      subroutine hnd_opt_tfds(rtdb)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      logical  gradients
      external gradients
      logical  geom_zmt_get_binv
      logical  geom_zmt_get_zsym
      logical  geom_hnd_tfds
      logical  status
      logical  dbug
      logical  iterate
      parameter (mxiter=5 )
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      common/hnd_iofile/ir,iw
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      dimension d(mxcoor)
      dimension z(mxzmat)
      dimension t(mxzmat)
      data zero /0.0D+00/
      data tol  /1.0d-03/
c
      dbug=.false.
c
c     ----- store internal coordinate displacement in -d- -----
c           -x - = cart. coord. of current point
c           -ds- = int. coord. displacement
c
      ncart=3*nat
      do i=1,nvar
         d(i)=ds(i)
      enddo
c
      iterate=.false.
      if(.not.iterate) then
c
c     ----- get -b- and -bi- and transform -ds- to -dx- -----              
c
         call hnd_bandbi(nat,x)
c
         status=ma_push_get(MT_DBL,ncart*nzvar,'mem bi',i_ibi,ibi)
         status=geom_zmt_get_binv(dbl_mb(ibi),ncart,nzvar,ncart)
         status=geom_hnd_tfds(d,dx,dbl_mb(ibi),nvar,ncart,ncart)
         status=ma_pop_stack(i_ibi)
         if(dbug) then
            write(iw,*) ' in -tfds- , ds = '
            write(iw,9999) (ds(i),i=1,nvar)
            write(iw,*) ' in -tfds- , dx = '
            write(iw,9999) (dx(i),i=1,ncart)
         endif
      else
c
c     ----- if accuracy is requested, iterate the transformation -----
c           of -ds- to -dx- to get the correct -dx- displacement
c
      endif      
c
      return
 9999 format(f15.10)
      end
      subroutine hnd_opt_tfdx(rtdb)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      logical  gradients
      external gradients
      logical  geom_zmt_get_bmat
      logical  geom_hnd_tfdx
      logical  status
      logical  dbug
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      common/hnd_iofile/ir,iw
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      dimension d(mxcoor)
c
      dbug=.false.
c
c     ----- get -b- and -bi- -----
c
      call hnd_bandbi(nat,x)
c
c     ----- transform -dx- to -ds- -----              
c
      ncart=3*nat
      do i=1,nvar
         d(i)=dx(i)
      enddo
      status=ma_push_get(MT_DBL,ncart*nzvar,'mem b ',i_ib,ib)
      status=geom_zmt_get_bmat(dbl_mb(ib),ncart,nzvar,ncart)
      status=geom_hnd_tfdx(d,dx,dbl_mb(ib),nvar,ncart,ncart)
      status=ma_pop_stack(i_ib)
      if(dbug) then
         write(iw,*) ' in -tfdx- , dx = '
         write(iw,9999) (dx(i),i=1,ncart)
         write(iw,*) ' in -tfdx- , ds = '
         write(iw,9999) (d(i),i=1,nvar)
      endif
c
      return
 9999 format(f15.10)
      end
      logical function scf_opt(rtdb)
      implicit none
      integer  rtdb
      scf_opt=.true.
      return
      end
      logical function scf_freq(rtdb)
      implicit none
      integer  rtdb
      scf_freq=.true.
      return
      end
      subroutine hnd_opt_hss_init(rtdb)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "context.fh"
c
      integer     rtdb
      logical     status
c
c parameters read from rtdb
c prefix:energy (later) rohf and uhf
c
      integer geom               ! geometry handle
      integer nat                ! number of atoms
      character*255 name         ! name buffer for rtdb context
      character*255 title        ! title buffer
      double precision energy    ! energy (scf, rohf, uhf)
      logical converged          ! energy convergence flag
c
      integer h_first, k_first   ! MA handle/index of first mem location
      integer h_hss  , k_hss     ! MA handle/index 
      integer h_hssx , k_hssx    ! MA handle/index 
      integer h_hssq , k_hssq    ! MA handle/index 
      integer h_binv , k_binv    ! MA handle/index 
      integer h_zsym , k_zsym    ! MA handle/index 
      integer h_t    , k_t       ! MA handle/index 
      integer h_zmat , k_zmat    ! MA handle/index 
      integer h_izmat, k_izmat   ! MA handle/index 
      integer h_h    , k_h       ! MA handle/index 
      integer h_g    , k_g       ! MA handle/index
      integer h_v    , k_v       ! MA handle/index
      integer h_ia   , k_ia      ! MA handle/index
      integer h_last , k_last    ! MA handle/index of last  mem location
c
      logical geom_zmt_get_binv
      logical geom_zmt_get_zmat
      logical geom_zmt_get_zsym
      logical geom_zmt_get_izmat
      logical dbug
      logical out 
c
      character*16 atom_tags
c
      parameter (zero=0.0d+00)
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)     
      parameter (mxcoor=1500)
      logical zcoord
      common/hnd_iofile/ir,iw
      common/hnd_optvar/zcoord,ncoord,mcoord
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      dimension atom_tags(mxatom)
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(dbug)  then
         write(iw,9999)
      endif
c
      status = rtdb_parallel(.false.)
      if (ga_nodeid().eq.0)  then
         if (.not. geom_create(geom, 'geometry'))
     &       call errquit('hnd_opt: geom_create?', 911)
         if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &       call errquit('hnd_opt: no geometry ', 911)
         if (.not. geom_ncent(geom,nat))
     &       call errquit('hnd_opt: geom_ncent?',911)
         call ga_brdcst(14+msgint,nat   ,mitob(1     ),0)
      else
         call ga_brdcst(14+msgint,nat   ,mitob(1     ),0)
      endif
      call ga_sync()
      status = rtdb_parallel(.true.) 
      if (dbug)  then
          write(iw,*)'   in hnd_opt_hss_init, nat =',nat
      endif
c
c     ----- define memory pointers -----
c
      ndimx=3*nat
      ndimq=nzvar
      ndim =mcoord           
c
      need=     ndim *ndim 
      need=need+ndimx*ndimx
      need=need+ndimq*ndimq
      need=need+ndimq
      need=need+abs(nzmat)
      if (.not.
     &   MA_Push_Get(mt_dbl,need       ,'mem first',h_first ,k_first ))
     &   call errquit('hnd_opt_hss_init: malloc of first  failed?',911)
      k_hss  =k_first
      k_hssx =k_hss  +ndim *ndim 
      k_hssq =k_hssx +ndimx*ndimx
      k_zmat =k_hssq +ndimq*ndimq
      k_izmat=k_zmat +ndimq
      k_last =k_izmat+abs(nzmat)
c
      if(dbug) then
         write(iw,9997) k_first,k_last
         write(iw,*) ' k_first   = ',k_first
         write(iw,*) ' k_hss     = ',k_hss  
         write(iw,*) ' k_hssx    = ',k_hssx 
         write(iw,*) ' k_hssq    = ',k_hssq 
         write(iw,*) ' k_zmat    = ',k_zmat 
         write(iw,*) ' k_izmat   = ',k_izmat 
         write(iw,*) ' k_last    = ',k_last 
         do iat=1,mxatom
            zan(iat)=zero
            do ixyz=1,3
               c(ixyz,iat)=zero
            enddo
         enddo
      endif
c
c     ----- get cartesian coordinates -----  
c
      status = rtdb_parallel(.false.)
      if (ga_nodeid().eq.0)  then
         if(.not.
     &   geom_cart_get(geom,nat,atom_tags,c,zan))
     &   call errquit('hnd_opt_hss_init: geom_cart_get failed?',911)
         call ga_brdcst(15+msgint,nat   ,mitob(1     ),0)
         call ga_brdcst(16+msgdbl,c     ,mdtob(mxcart),0)
         call ga_brdcst(17+msgdbl,zan   ,mdtob(mxatom),0)
      else
         call ga_brdcst(15+msgint,nat   ,mitob(1     ),0)
         call ga_brdcst(16+msgdbl,c     ,mdtob(mxcart),0)
         call ga_brdcst(17+msgdbl,zan   ,mdtob(mxatom),0)
      endif
      call ga_sync()
      status = rtdb_parallel(.true.) 
      if(dbug) then
         write(iw,*) 'in hnd_opt_hss_init , z and c = '
         write(iw,9996) (zan(iat),(c(ixyz,iat),ixyz=1,3),iat=1,nat)
      endif
c
c     ----- create -zmat- data if requested -----
c
      if(nzmat.gt.0) then
         call hnd_bandbi(nat,c)
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            if(.not.geom_zmt_get_zmat(dbl_mb(k_zmat),nzvar))                 
     &      call errquit('hnd_opt_hss_init:zmt_get_zmat failed?',911)
            if(.not.geom_zmt_get_izmat(dbl_mb(k_izmat),nzmat))                 
     &      call errquit('hnd_opt_hss_init:zmt_get_izmat failed?',911)
            call ga_brdcst(24+msgdbl,dbl_mb(k_zmat ),mdtob(nzvar),0)
            call ga_brdcst(25+msgdbl,dbl_mb(k_izmat),mdtob(nzmat),0)
         else
            call ga_brdcst(24+msgdbl,dbl_mb(k_zmat ),mdtob(nzvar),0)
            call ga_brdcst(25+msgdbl,dbl_mb(k_izmat),mdtob(nzmat),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.)
      endif
c
c     ----- initial guess hessian ( cartesian or internal ) -----
c
      call hnd_opt_hss_init_x(dbl_mb(k_hss),ndim,
     &             dbl_mb(k_hssx),ndimx,dbl_mb(k_hssq),ndimq,          
     &             dbl_mb(k_zmat),dbl_mb(k_izmat))
c
c     ----- release memory -----
c
      if (.not.MA_Pop_Stack(h_first))
     &   call errquit('hnd_opt_hss_init: ma_pop_stack k_first   ?',911)
      status=rtdb_parallel(.false.)
      if( ga_nodeid().eq.0) then
         if (.not.geom_destroy(geom))
     &   call errquit('hnd_opt: geom_destroy?', 911)
      endif
      call ga_sync()
      status=rtdb_parallel(.true.)
c
      if(dbug) then
         write(iw,9998)
      endif
      return 
 9999 format(' start of hnd_opt_hss_init ')
 9998 format('   end of hnd_opt_hss_init ')
 9997 format('       in hnd_opt_hss_init , k_first, k_last = ',2i10)
 9996 format(' z = ',f6.2,' c = ',3f15.8)
      end
      SUBROUTINE HND_OPT_HSS_INIT_X(
     1           HSS,NCOOR,HSSX,NCART,HSSQ,NINT,ZMAT,IZMAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
C
C      -----        INITIALIZE THE HESSIAN MATRIX           -----
C      ----- IF -INHESS- EQUALS 1, INPUT THE HESSIAN MATRIX -----
C                 OTHERWISE MAKE UP INITIAL GUESS
C              IN CARTESIAN OR INTERNAL COORDINATES
C
      INTEGER   GA_HESS
      LOGICAL   DRIVER_OPT_PUT_HESS
      LOGICAL   DRIVER_OPT_GET_HESS
      LOGICAL   STATUS
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      CHARACTER*8 ERRMSG
      LOGICAL OUT
      LOGICAL SOME
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_GA_OPT/GA_HESS
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION HSSX(NCART,*),HSSQ(NINT,*)       
      DIMENSION HSS (NCOOR,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-OPTHSS-'/
      DATA ZERO   /0.0D+00/
      DATA HDIAGX /0.5D+00/
C
      OUT =.FALSE.
      SOME=.FALSE.
      SOME=SOME.OR.OUT    
C
      IF(OUT) THEN
         WRITE(IW,9999) NAT,NCART,NZVAR,NVAR,NINT,NCOOR,NZMAT
      ENDIF
      DO J=1,NCOOR
         DO I=1,NCOOR
            HSS(I,J)=ZERO
         ENDDO
      ENDDO
C
      IF(NZMAT.NE.0) THEN
C
C     ----- INTERNALLY GENERATED INITIAL HESSIAN -----
C
         IF(OUT) THEN
            WRITE(IW,9999) NAT,NCART,NZVAR,NVAR,NINT,NCOOR,NZMAT
         ENDIF
         CALL HND_OPT_HSS_INIT_ZMT(HSSQ,NINT,ZMAT,IZMAT,OUT)
         DO J=1,NZVAR
            DO I=1,NZVAR 
               HSS(I,J)=HSSQ(I,J) 
            ENDDO
         ENDDO
         IF(SOME) THEN
            WRITE(IW,9998)
            CALL HND_PRSQ(HSS,NZVAR,NZVAR,NCOOR)
         ENDIF
      ELSE
C
C     ----- CARTESIAN HESSIAN SET TO 0.5 * UNIT MATRIX -----
C
         DO J=1,NCART 
            DO I=1,NCART 
               HSSX(I,J)=ZERO
            ENDDO
            HSSX(J,J)=HDIAGX
         ENDDO
         DO J=1,NCART 
            DO I=1,NCART
               HSS(I,J)=HSSX(I,J) 
            ENDDO
         ENDDO
         IF(SOME) THEN
            WRITE(IW,9997)
            CALL HND_PRSQ(HSS,NCART,NCART,NCOOR)
         ENDIF
      ENDIF
C
C     ----- SAVE INITIAL HESSIAN -----
C
      STATUS=GA_CREATE(MT_DBL,NCOOR,NCOOR,'DRIVER:HESS',1,1,GA_HESS)
      IF(.NOT.STATUS) THEN
         CALL ERRQUIT('HSS_INIT_X:GA_CREATE?',911)
      ENDIF
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF( GA_NODEID().EQ.0) THEN
         STATUS=DRIVER_OPT_PUT_HESS(HSS,NCOOR,NCOOR)
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
      RETURN
 9999 FORMAT(/,' IN  OPT_HSS-INIT_X , ', 
     1         ' NAT,NCART,NZVAR,NVAR,NINT,NCOOR,NZMAT = ',7I4)
 9998 FORMAT(/,10X,38(1H-),/,10X,'INITIAL -INTERN.COORD.- ',
     1                           'HESSIAN MATRIX',/,10X,38(1H-))
 9997 FORMAT(/,10X,38(1H-),/,10X,'INITIAL - CART.COORD. - ',
     1                           'HESSIAN MATRIX',/,10X,38(1H-))
      END
      SUBROUTINE HND_OPT_HSS_INIT_ZMT(HSSQ,NINT,ZMAT,IZMAT,SOME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (ZERO=0.0D+00,      ONE=1.0D+00,
     1           QUARTR=0.25D+00,   PT5=0.5D+00,
     2           THIRD=ONE/3.0D+00, TOHB2=6.423D-02,
     3           TOANGS=0.52917724924D+00)
      PARAMETER (MXATOM=500)
      LOGICAL OUT 
      LOGICAL DBUG
      LOGICAL BADZ
      LOGICAL NATURL 
      LOGICAL SYMTRC
      LOGICAL INDEP
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION HSSQ(NINT,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION IROW(86),DPARM(6,6),APARM(6,6),RCOV(86)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-HSSGSS-'/
C
C     --- DATA TABLES FOR BADGER'S RULES (SEE J.C.P. 2, 128(1934)) ---
C        A GENERALIZED BADGER'S RULE, AND THESE PARAMETERS ARE FROM
C        D.R.HERSCHBACH, V.W.LAURIE, J.CHEM.PHYS. 35, 458-463(1961).
C
      DATA IROW  /2*1,8*2,8*3,18*4,18*5,32*6/
      DATA APARM /1.26D+00,1.66D+00,1.84D+00,1.98D+00,2.03D+00,2.03D+00,
     1            1.66D+00,1.91D+00,2.28D+00,2.35D+00,2.33D+00,2.50D+00,
     2            1.84D+00,2.28D+00,2.41D+00,2.52D+00,2.61D+00,2.60D+00,
     3            1.98D+00,2.35D+00,2.52D+00,2.58D+00,2.66D+00,2.75D+00,
     4            2.03D+00,2.33D+00,2.61D+00,2.66D+00,2.85D+00,2.76D+00,
     5            2.03D+00,2.50D+00,2.60D+00,2.75D+00,2.75D+00,3.00D+00/
      DATA DPARM /2.50D-02,0.30D+00,0.38D+00,0.49D+00,0.51D+00,0.81D+00,
     1            0.30D+00,0.68D+00,0.74D+00,0.85D+00,0.68D+00,0.97D+00,
     2            0.38D+00,0.74D+00,1.18D+00,1.02D+00,1.28D+00,0.84D+00,
     3            0.49D+00,0.85D+00,1.02D+00,1.41D+00,0.86D+00,1.14D+00,
     4            0.51D+00,0.68D+00,1.28D+00,0.86D+00,1.62D+00,1.25D+00,
     5            0.81D+00,0.97D+00,0.84D+00,1.14D+00,1.25D+00,1.83D+00/
C
C      COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
C      EXCEPT VAN DER WAALS RADII FOR HE,NE,AR,KR (SAME SOURCE),
C      AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT,RN
C
      DATA (RCOV(NUCZ),NUCZ=1,2)/0.30D+00,1.22D+00/
      DATA (RCOV(NUCZ),NUCZ=3,10)
     1  /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     2   0.70D+00,0.66D+00,0.58D+00,1.60D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     1  /1.66D+00,1.36D+00,1.25D+00,1.17D+00,
     2   1.10D+00,1.04D+00,0.99D+00,1.91D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     1  /2.03D+00,1.74D+00,
     2   1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     3  1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     4   1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.98D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     1  /2.22D+00,1.92D+00,
     2   1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     3   1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     4   1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,2.09D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     1  /2.35D+00,1.98D+00,
     2   1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D+00,
     3   1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D+00,
     4   1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     5   1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     6   1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,2.20D+00/
C
      DATA GUESS,GUESS1,GUESS2 /8HGUESS   ,8HGUESS1  ,8HGUESS2  /
C
      DIST(IAT,JAT)=SQRT((C(1,IAT)-C(1,JAT))**2+
     1                   (C(2,IAT)-C(2,JAT))**2+
     1                   (C(3,IAT)-C(3,JAT))**2) * TOANGS
      BADZ(NUCZ) = NUCZ.LE.0.OR.NUCZ.GT.86
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      IF(OUT) THEN
         WRITE(IW,9168)
         WRITE(IW,9188) NAT,NVAR,NZVAR,NINT,NZMAT
      ENDIF
C
C     ----- THIS ROUTINE MAKES A DIAGONAL GUESS OF THE HESSIAN -----
C                     IN INTERNAL COORDINATE SPACE.
C
      DO J=1,NZVAR
         DO I=1,NZVAR
            HSSQ(I,J)=ZERO
         ENDDO
      ENDDO
C
C     -----          INTERNAL COORDINATES                      -----
C
C     FOR MORE INFORMATION ON -GUESS1- OR -GUESS2- SEE
C     H.B.SCHLEGEL, THEORET.CHIM.ACTA, 66, 333-340 (1984).
C     T.H.FISCHER, J.ALMLOF, J.PHYS.CHEM. 96, 9786-9774 (1992).
C     -GUESS1- AND -GUESS2- ARE UNDOCUMENTED WAYS TO GET PURE
C     BADGER'S RULES, OR PURE FISCHER/ALMLOF RULES.
C     THE DEFAULT -GUESS- USES F/A FOR STRETCHES AND BENDS,
C     BUT NOT FOR TORSIONS OR OUT-OF-PLANE BENDS.
C
      IF(DBUG) THEN
         WRITE(IW,9208) (IZMAT(IZ),IZ=1,NZMAT)
         WRITE(IW,9128) ( ZMAT(IZ),IZ=1,NZVAR)
      ENDIF
      HGUESS=GUESS2
      IF(OUT) THEN
         WRITE(IW,9008) HGUESS
      ENDIF
C
      IDEGEN=0
      IZ=1
      DO 190 I=1,NZVAR
         ITYPE=IZMAT(IZ)
         IZ1 = 0
         IZ2 = 0
         IZ3 = 0
         IZ4 = 0
         IZA = 0
         IZB = 0
         IZC = 0
         IZD = 0
C
C                  BOND STRETCH, B-A
C
         IF(ITYPE.GT.1) GO TO 120
         IZ1 = IZMAT(IZ+1)
         IZ2 = IZMAT(IZ+2)
         IZA = NUC(IZ1)
         IZB = NUC(IZ2)
         IF(BADZ(IZA)  .OR.  BADZ(IZB)) THEN
            HSSQ(I,I) = THIRD
            GO TO 180
         ENDIF
C
         IROWA = IROW(IZA)
         IROWB = IROW(IZB)
         AIJ = APARM(IROWA,IROWB)
         DIJ = DPARM(IROWA,IROWB)
         RRR = TOANGS*ZMAT(I)
         HII = (AIJ-DIJ)/(RRR-DIJ)
         HESS1 = TOHB2 * HII * HII * HII
C
         RAB = DIST(IZ1,IZ2)
         CAB = RCOV(IZA) + RCOV(IZB)
         HESS2 = 0.3601 * EXP(-1.944*(RAB - CAB))
C
         IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS2
         IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
         IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
         GO TO 180
C
C                  ANGLE BEND, B-A-C
C
  120    CONTINUE
         IF(ITYPE.GT.2) GO TO 130
         IZ1 = IZMAT(IZ+1)
         IZ2 = IZMAT(IZ+2)
         IZ3 = IZMAT(IZ+3)
         IZB = NUC(IZ1)
         IZA = NUC(IZ2)
         IZC = NUC(IZ3)
         IF(BADZ(IZA)  .OR.  BADZ(IZB)  .OR.  BADZ(IZC)) THEN
            HSSQ(I,I) = QUARTR
            GO TO 180
         ENDIF
C
         HESS1 = QUARTR
         IF(IZB.EQ.1  .OR.  IZC.EQ.1) HESS1 = 0.16D+00
C
         RAB = DIST(IZ2,IZ1)
         RAC = DIST(IZ2,IZ3)
         CAB = RCOV(IZA) + RCOV(IZB)
         CAC = RCOV(IZA) + RCOV(IZC)
         HESS2 = 0.089 + 0.11 * (CAB*CAC)**0.42
     1                   * EXP(-0.44*(RAB + RAC - CAB - CAC))
C
         IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS2
         IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
         IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
         GO TO 180
C
C                  TORSION, D-A-B-C
C
  130    CONTINUE
         IF(ITYPE.GT.3) GO TO 140
         IZ1 = IZMAT(IZ+1)
         IZ2 = IZMAT(IZ+2)
         IZ3 = IZMAT(IZ+3)
         IZ4 = IZMAT(IZ+4)
         IZD = NUC(IZ1)
         IZA = NUC(IZ2)
         IZB = NUC(IZ3)
         IZC = NUC(IZ4)
         IF(BADZ(IZA) .OR. BADZ(IZB) .OR. BADZ(IZC) .OR. BADZ(IZD)) THEN
            HSSQ(I,I) = QUARTR
            GO TO 180
         ENDIF
C
         HESS1 = QUARTR
C
C            COUNT BONDS TO ATOMS A AND B, EXCEPT FOR A-B ITSELF
C
         NBA = 0
         NBB = 0
         DO 135 J=1,NAT
            IZJ = NUC(J)
            IF(IZJ.LE.0) GO TO 135
            IF(IZJ.LE.86) THEN
               CAJ = RCOV(IZA) + RCOV(IZJ)
               CBJ = RCOV(IZB) + RCOV(IZJ)
            ELSE
               CAJ = RCOV(IZA) + 1.70D+00
               CBJ = RCOV(IZB) + 1.70D+00
            ENDIF
            IF(DIST(IZ2,J)  .LT.  1.2D+00*CAJ  .AND.
     1         J.NE.IZ2  .AND.  J.NE.IZ3) NBA = NBA+1
            IF(DIST(IZ3,J)  .LT.  1.2D+00*CBJ  .AND.
     1         J.NE.IZ3  .AND.  J.NE.IZ2) NBB = NBB+1
  135    CONTINUE
C
C     THE FACTOR 14.0 IN PREEXP SEEMS LARGE, SO BRACKET THE RESULT
C
         RAB = DIST(IZ2,IZ3)
         CAB = RCOV(IZA) + RCOV(IZB)
         PREEXP = 14.0 * (NBA+NBB)**0.57 / (RAB*CAB)**4
         HESS2 = 0.0015 + PREEXP * EXP(-2.85*(RAB-CAB))
         IF(HESS2.GT.0.30D+00) HESS2 = 0.30D+00
         IF(HESS2.LT.0.05D+00) HESS2 = 0.05D+00
C
         IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS1
         IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
         IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
         GO TO 180
C
C                  OUT OF PLANE BEND, X-A-B-C, X IS HOOKED TO A
C
  140    CONTINUE
         IF(ITYPE.GT.4) GO TO 150
         IZ1 = IZMAT(IZ+1)
         IZ2 = IZMAT(IZ+2)
         IZ3 = IZMAT(IZ+3)
         IZ4 = IZMAT(IZ+4)
         IZX = NUC(IZ1)
         IZA = NUC(IZ2)
         IZB = NUC(IZ3)
         IZC = NUC(IZ4)
         IF(BADZ(IZA) .OR. BADZ(IZB) .OR. BADZ(IZC) .OR. BADZ(IZX)) THEN
            HSSQ(I,I) = QUARTR
            GO TO 180
         ENDIF
C
         HESS1 = QUARTR
C
         CAX = RCOV(IZA) + RCOV(IZX)
         CAB = RCOV(IZA) + RCOV(IZB)
         CAC = RCOV(IZA) + RCOV(IZC)
         RAX = DIST(IZ1,IZ2)
         PREEXP = 0.0061 * (CAB*CAC)**0.80 * (COS(ZMAT(I)))**4.00
         HESS2 = 0.0025 + PREEXP * EXP(-3.00*(RAX-CAX))
C
         IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS1
         IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
         IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
         GO TO 180
C
C                  ALL OTHER COORDINATE TYPES...
C
  150    CONTINUE
         HSSQ(I,I)=QUARTR
C
  180    CONTINUE
         IF(ITYPE.EQ.1) IZ=IZ+3
         IF(ITYPE.EQ.2) IZ=IZ+4
         IF(ITYPE.EQ.3) IZ=IZ+5
         IF(ITYPE.EQ.4) IZ=IZ+5
         IF(ITYPE.EQ.5) THEN
            IDEGEN=IDEGEN+1
            IF(IDEGEN.EQ.2) IZ=IZ+4
            IF(IDEGEN.EQ.2) IDEGEN=0
         ENDIF
         IF(ITYPE.EQ.6) IZ=IZ+6
         IF(ITYPE.EQ.7) IZ=IZ+7
         IF(OUT) THEN
            WRITE(IW,9148) I,ITYPE,IZ1,IZ2,IZ3,IZ4,
     1                             IZA,IZB,IZC,IZD,HSSQ(I,I)
         ENDIF
  190 CONTINUE
C
      IF(OUT) THEN
         WRITE(IW,9048)
         CALL HND_PRSQ(HSSQ,NZVAR,NZVAR,NINT)
      ENDIF
C
      RETURN
 9208 FORMAT(' IZMAT = ',12I4)
 9188 FORMAT(' IN  OPT_HSS-INIT_ZMT ,', 
     1       ' NAT,NVAR,NZVAR,NINT,NZMAT = ',5I4)
 9168 FORMAT(' IN  OPT_HSS_INIT_ZMT , GUESS INTERNAL HESSIAN ')
 9148 FORMAT(' I,ITYPE,HSSQ(I,I) = ',2I5,8I3,F12.8)
 9128 FORMAT('  ZMAT = ',5F7.3)
 9048 FORMAT(/,1X,38(1H-),
     1       /,1X,'INITIAL -INTERN.COORD.- HESSIAN MATRIX',
     2       /,1X,38(1H-))
 9008 FORMAT(' IN  OPT_HSS_INIT_ZMT , HGUESS = ',A8) 
      END
      logical function driver_opt_put_hess(hess,nint,mint)   
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
c
      integer ga_hess
      logical status
c
      integer nint 
      integer mint
      integer i
      integer j
      integer ir
      integer iw
      integer ip
      logical dbug
      common/hnd_iofile/ir,iw
      common/hnd_ga_opt/ga_hess
      double precision hess(mint,*)  
c
      dbug=.false.
c
      if(dbug) then
         write(iw,*) ' _put_hess '
         call hnd_prsq(hess,nint,nint,mint)    
      endif
c
      call ga_put(ga_hess,1,nint,1,nint,hess,mint)
c
      driver_opt_put_hess=.true.
      return
      end
      logical function driver_opt_get_hess(hess,nint,mint)   
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
c
      integer ga_hess
      logical status
c    
      integer nint
      integer mint
      integer i
      integer j
      integer ir
      integer iw
      integer ip
      logical dbug
      common/hnd_iofile/ir,iw
      common/hnd_ga_opt/ga_hess
      double precision hess(mint,*) 
c
      dbug=.false.
c
      call ga_get(ga_hess,1,nint,1,nint,hess,mint)
c
      if(dbug) then
         write(iw,*) ' _get_hess'
         call hnd_prsq(hess,nint,nint,mint)    
      endif
c
      driver_opt_get_hess=.true.
      return
      end
      subroutine hnd_opt_search_dir(rtdb)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "context.fh"
c
      integer    rtdb
      logical    status
      logical    driver_opt_get_hess
      logical    driver_opt_put_hess
      logical    geom_zmt_get_zsym
c
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      logical dbug
      logical out
      logical baker
      logical zcoord
      common/hnd_iofile/ir,iw
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_optvar/zcoord,ncoord,mcoord
      dimension x(1)
      equivalence (x(1),dbl_mb(1))
      data zero /0.0d+00/
c
      dbug =.false.
      out  =.false.
      out  =out.or.dbug
c
c     ----- set broadcast off -----
c
      baker=.false.
      if(out) then
         write(iw,*) 'start of hnd_opt_search_dir, baker = ',
     1                                             baker
      endif
c
      ndimx=3*nat
      ndimq=nzvar
      ndimi=ncoord
      ndim =mcoord           
      if(baker) then
         ndim1=ndimi+1
      else
         ndim1=ndimi
      endif
      if(dbug) then
         write(iw,*) 'in search_dir,',
     1   ' ndimx,ndimq,ndimi,ndim,ndim1 = ', 
     2     ndimx,ndimq,ndimi,ndim,ndim1
      endif
c
      need=      ndim1* ndim1
      need=need+ ndim1* ndim1
      need=need+(ndim1*(ndim1+1))/2
      need=need+ ndim1
      need=need+ ndim1
      need=need+ ndim1
      need=need+ ndim1
      if(zcoord) then
         need=need+ ndim *ndim 
         need=need+ ndim *ndim 
         need=need+ ndim 
      endif
c
c     ----- get memory -----
c
      if(.not.ma_push_get(mt_dbl,
     1               need,'mem first',i_ifirst,ifirst))
     2   call errquit('hnd_opt_search_dir: get first failed?',911)
c
      i10=ifirst 
      i20=i10+ ndim1* ndim1
      i30=i20+ ndim1* ndim1
      i40=i30+(ndim1*(ndim1+1))/2
      i50=i40+ ndim1
      i60=i50+ ndim1
      i70=i60+ ndim1
      i80=i70+ ndim1
      if(zcoord) then
         j10=i80
         j20=j10+ ndim *ndim 
         j30=j20+ ndim *ndim 
         j40=j30+ ndim 
      endif
c
c     ----- get current hessian matrix -----
c
      if(zcoord) then
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_get_hess(x(j10),ndimq,ndim)     
            call ga_brdcst(26+msgdbl,x(j10),mdtob(ndim*ndim),0)
         else
            call ga_brdcst(26+msgdbl,x(j10),mdtob(ndim*ndim),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.) 
         if(dbug) then
            write(iw,*) 'in -search_dir- , int.coord. hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
         endif
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=geom_zmt_get_zsym(x(j20),ndimq,ndimq,ndim)
            call ga_brdcst(27+msgdbl,x(j20),mdtob(ndim*ndim),0)
         else
            call ga_brdcst(27+msgdbl,x(j20),mdtob(ndim*ndim),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.) 
         if(dbug) then
            write(iw,*) 'in -search_dir- , zsym = '                 
            call hnd_prsq(x(j20),ndimq,ndimq,ndim)
         endif
         call hnd_tfhs(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
         if(out) then
            write(iw,*) 'in -search_dir- , transformed hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
         endif
         do j=1,ndimi
            do i=1,ndimi
               x(i+ndim1*(j-1)+i10-1)=x(i+ndim*(j-1)+j10-1)
            enddo
         enddo
      else
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_get_hess(x(i10),ndimx,ndim1)       
            call ga_brdcst(28+msgdbl,x(i10),mdtob(ndim1*ndim1),0)
         else
            call ga_brdcst(28+msgdbl,x(i10),mdtob(ndim1*ndim1),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.) 
         if(dbug) then
            write(iw,*) 'in -search_dir- , cartesian hessian = '
            call hnd_prsq(x(i10),ndimx,ndimx,ndim1)
         endif
      endif
c
c     --- search direction (quasi-newton or augmented hessian) ---
c
         if(baker) then
            call hnd_opt_bkr_search(ncoord,
     1      x(i10),x(i20),x(i30),x(i40),x(i50),x(i60),x(i70),
     2      ndim1)
         else
            call hnd_opt_dir_search(ncoord,
     1      x(i10),x(i20),x(i30),x(i40),x(i50),x(i60),x(i70),
     2      ndim1)
         endif
c
c     ----- transform back if needed and store hessian -----
c
      if(zcoord) then
         do j=1,ndimq
            do i=1,ndimq
               x(i+ndim*(j-1)+j10-1)=Zero                  
            enddo
         enddo
         do j=1,ndimi
            do i=1,ndimi
               x(i+ndim*(j-1)+j10-1)=x(i+ndim1*(j-1)+i10-1)
            enddo
         enddo
         if(dbug) then
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
         endif
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=geom_zmt_get_zsym(x(j20),ndimq,ndimq,ndim)
            call ga_brdcst(29+msgdbl,x(j20),mdtob(ndim*ndim),0)
         else
            call ga_brdcst(29+msgdbl,x(j20),mdtob(ndim*ndim),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.) 
c
         call hnd_tfhsi(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
c
         if(out) then
            write(iw,*) 'in -search_dir- , ',
     1                  'back-transformed hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
            if(dbug) then
               call hnd_tfhs(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
               call hnd_prsq(x(j10),ndimq,ndimq,ndim)
            endif
         endif
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_put_hess(x(j10),ndimq,ndim)       
         endif
         call ga_sync()
         status=rtdb_parallel(.true.) 
      else
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_put_hess(x(i10),ndimx,ndim1)       
         endif
         call ga_sync()
         status=rtdb_parallel(.true.) 
         if(dbug) then
            write(iw,*) 'in -search_dir- , cartesian hessian = '
            call hnd_prsq(x(i10),ndimx,ndimx,ndim1)
         endif
      endif
c
c     ----- release memory -----
c
      if(.not.ma_pop_stack(i_ifirst))
     1   call errquit('hnd_opt_search_dir: pop first failed?',911)
c
      if(dbug) then
         write(iw,*) '   end of hnd_opt_search_dir '
      endif
c
      return
      end
      SUBROUTINE HND_OPT_BKR_SEARCH(NVAR,
     1           HESS,HESVEC,HESST,HESEIG,HESEDM,T,IA,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- SELECT SEARCH DIRECTION -----
C
      LOGICAL   OUT
      LOGICAL   DBUG
      LOGICAL   CVGED
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED
      COMMON/HND_OPTOPT/CVGOPT,NPTOPT,INHESS,LINOPT,MODUPD,MODDIR
      COMMON/HND_OPTFRZ/ICFRZ(MXCOOR),IATFRZ(MXATOM)
      DIMENSION HESS(NDIM,*),HESVEC(NDIM,*),HESST(*),HESEIG(*)
      DIMENSION HESEDM(*),T(*),IA(*)
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TINY   /1.0D-06/
      DATA SMALL  /1.0D-04/
      DATA DAMP1  /5.0D+00/
      DATA DAMP2  /2.0D+00/
      DATA HUNDRD /1.0D+02/
      DATA TENM1  /1.0D-01/
C
      DBUG=.FALSE.
      OUT =.TRUE. 
      OUT =OUT.OR.DBUG
      IF(OUT) THEN
         WRITE(IW,9999) NVAR,NDIM
      ENDIF
C
      DO I=1,NDIM
         IA(I)=(I*(I-1))/2
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,9988)
         DO J=1,NVAR
            WRITE(IW,9995) J,G(J)                 
         ENDDO
      ENDIF
C
C     ----- ASSEMBLE AUGMENTED HESSIAN -----
C
      DO I=1,NVAR
         HESS(I,NVAR+1)=G(I)
         HESS(NVAR+1,I)=G(I)
      ENDDO
      HESS(NVAR+1,NVAR+1)=ZERO
      IF(DBUG) THEN
         WRITE(IW,9997)
         CALL HND_PRSQ(HESS,NVAR+1,NVAR+1,NDIM)
      ENDIF
C
C     ----- DIAGONALIZE THE AUGMENTED HESSIAN MATRIX -----
C
      IJ=0
      DO I=1,NVAR+1
         DO J=1,I
            IJ=IJ+1
            HESST(IJ)=HESS(I,J)
         ENDDO
      ENDDO
      CALL HND_DIAGIV(HESST,HESVEC,HESEIG,IA,NVAR+1,NVAR+1,NDIM)
      IF(DBUG) THEN
         WRITE(IW,9996)
         CALL HND_PREV(HESVEC,HESEIG,NVAR+1,NVAR+1,NDIM)
      ELSEIF(OUT) THEN
         WRITE(IW,9996)
         DO J=1,NVAR+1
            WRITE(IW,9995) J,HESEIG(J)                 
         ENDDO
      ENDIF
C
C     ----- THIS IS THE AUGMENTED HESSIAN DIRECTION -----
C
      DUM=ONE/HESVEC(NVAR+1,1)
      DO I=1,NVAR+1
         HESVEC(I,1)=HESVEC(I,1)*DUM
      ENDDO
      DO I=1,NVAR
         DX(I)=HESVEC(I,1)
      ENDDO
      IF(DBUG) THEN
         CALL HND_PRSQ(HESVEC,1,NVAR+1,NDIM)
      ENDIF
      IF(OUT) THEN
         WRITE(IW,9994)
         DO I=1,NVAR
            WRITE(IW,9995) I,DX(I)
         ENDDO
      ENDIF
C
      RETURN
 9999 FORMAT(' IN BKR_SEARCH, NVAR,NDIM = ',2I5)            
 9998 FORMAT(' IN BKR_SEARCH, CURRENT HESSIAN ')
 9997 FORMAT(' IN BKR_SEARCH, AUGMENTED HESSIAN ')
 9996 FORMAT(' IN BKR_SEARCH, HESSIAN EIGENMODES = ')
 9995 FORMAT(I5,2F12.8)
 9994 FORMAT(' IN BKR_SEARCH, AUGMENTED HESSIAN STEP = ')
 9988 FORMAT(' IN BKR_SEARCH, CURRENT GRADIENT ')
      END
      SUBROUTINE HND_OPT_DIR_SEARCH(NVAR,
     1           HESS,HESVEC,HESST,HESEIG,HESEDM,T,IA,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- SELECT SEARCH DIRECTION -----
C
      LOGICAL   OUT
      LOGICAL   DBUG
      LOGICAL   CVGED
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED
      COMMON/HND_OPTOPT/CVGOPT,NPTOPT,INHESS,LINOPT,MODUPD,MODDIR
      COMMON/HND_OPTFRZ/ICFRZ(MXCOOR),IATFRZ(MXATOM)
      DIMENSION HESS(NDIM,*),HESVEC(NDIM,*),HESST(*),HESEIG(*)
      DIMENSION HESEDM(*),T(*),IA(*)
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TINY   /1.0D-06/
      DATA SMALL  /1.0D-04/
      DATA DAMP1  /5.0D+00/
      DATA DAMP2  /2.0D+00/
      DATA HUNDRD /1.0D+02/
      DATA TENM1  /1.0D-01/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(DBUG) THEN
         WRITE(IW,9999) NVAR,NDIM
      ENDIF
C
      DO I=1,NDIM
         IA(I)=(I*(I-1))/2
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9988)
         DO J=1,NVAR
            WRITE(IW,9995) J,G(J)                 
         ENDDO
      ENDIF
C
C     ----- DIAGONALIZE THE HESSIAN MATRIX -----
C
      IJ=0
      DO I=1,NVAR
         DO J=1,I
            IJ=IJ+1
            HESST(IJ)=HESS(I,J)
         ENDDO
      ENDDO
      CALL HND_DIAGIV(HESST,HESVEC,HESEIG,IA,NVAR,NVAR,NDIM)
      IF(DBUG) THEN
         WRITE(IW,9997)
         CALL HND_PREV(HESVEC,HESEIG,NVAR,NVAR,NDIM)
      ENDIF
C
C     ----- PROJECT GRADIENT VECTOR ON TO HESSIAN EIGEN MODES -----
C
      DO J=1,NVAR
         DUM=ZERO
         DO K=1,NVAR
            DUM=DUM+G(K)*HESVEC(K,J)
         ENDDO
         HESEDM(J)=DUM
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9996)
         DO J=1,NVAR
            WRITE(IW,9995) J,HESEIG(J),HESEDM(J)
         ENDDO
      ENDIF
C
C     --- CHECK FOR NEGATIVE EIGENVALUES AND CONVERT TO POSITIVE ----
C               AFTER DAMPING. ALSO, DAMP ANY ROTATIONAL 
C                      OR TRANSLATIONAL MODES.
C
      NEGEIG=0
      DO I=1,NVAR
         IF(ABS(HESEIG(I)).LT.TINY) THEN
            HESEIG(I)=ZERO
            HESEDM(I)=ZERO
         ELSE 
            IF(HESEIG(I).LT.ZERO) THEN
               NEGEIG=NEGEIG+1
               EIGNEG=HESEIG(I)
               EIGMOD=ABS(HESEIG(I))*DAMP1
               HESEIG(I)=EIGMOD
               IF(OUT) THEN
                  WRITE(IW,9994) I,EIGNEG,EIGMOD
               ENDIF
            ENDIF
            HESEDM(I)=HESEIG(I)
            IF((HESEIG(I).NE.ZERO ).AND.
     1         (HESEIG(I).LT.SMALL)     ) THEN
               EIGSML=HESEIG(I)
               EIGMOD=MAX(HESEIG(I)*DAMP2,SMALL)
               HESEDM(I)=EIGMOD
               IF(OUT) THEN
                  WRITE(IW,9993) I,EIGSML,EIGMOD
               ENDIF
            ENDIF
         ENDIF
      ENDDO     
C
C     ----- INVERT THE HESSIAN MATRIX -----
C
      DO I=1,NVAR
         IF(HESEDM(I).NE.ZERO) THEN
            HESEDM(I)=ONE/HESEDM(I)
         ELSE
            HESEDM(I)=ZERO
         ENDIF
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,9989)
         CALL HND_PREV(HESVEC,HESEDM,NVAR,NVAR,NDIM)
      ENDIF
      DO I=1,NVAR
         DO J=1,NVAR
            DUM=ZERO
            DO K=1,NVAR
               DUM=DUM+HESVEC(I,K)*HESEDM(K)*HESVEC(J,K)
            ENDDO   
            HESS(I,J)=DUM
         ENDDO   
      ENDDO   
      IF(DBUG) THEN
         WRITE(IW,9991)
         CALL HND_PRSQ(HESS,NVAR,NVAR,NDIM)
      ENDIF
C
C      ----- SELECT THE NEW DIRECTION ( NEWTON-RAPHSON'S ) -----
C
      DO I=1,NVAR
         DUM=ZERO
         DO K=1,NVAR
            DUM=DUM+HESS(I,K)*G(K)
         ENDDO
         IF(ICFRZ(I).NE.0) THEN
            DX(I)=ZERO
         ELSE
            DX(I)=-DUM
         ENDIF
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9992)
         DO I=1,NVAR
            WRITE(IW,9995) I,DX(I)
         ENDDO
      ENDIF
C
C     ----- REBUILD THE HESSIAN MATRIX -----
C
      DO I=1,NVAR
         DO J=1,NVAR
            DUM=ZERO
            DO K=1,NVAR
               DUM=DUM+HESVEC(I,K)*HESEIG(K)*HESVEC(J,K)
            ENDDO   
            HESS(I,J)=DUM
         ENDDO   
      ENDDO   
      IF(DBUG) THEN
         WRITE(IW,9990)
         CALL HND_PRSQ(HESS,NVAR,NVAR,NDIM)
      ENDIF
C
      RETURN
 9999 FORMAT(' IN DIR_SEARCH, NVAR,NDIM = ',2I5)            
 9998 FORMAT(' IN DIR_SEARCH, CURRENT HESSIAN ')
 9997 FORMAT(' IN DIR_SEARCH, HESSIAN EIGENMODES = ')
 9996 FORMAT(' IN DIR_SEARCH, GRADIENT PROJECTIONS ON EIGENMODES = ')
 9995 FORMAT(I5,2F12.8)
 9994 FORMAT(' IN DIR_SEARCH, -NEG EIG- CHANGED, I,EIGNEG,EIGMOD = ',
     1       I5,2F15.8)
 9993 FORMAT(' IN DIR_SEARCH, -SML EIG- CHANGED, I,EIGSML,EIGMOD = ',
     1       I5,2F15.8)
 9992 FORMAT(' IN DIR_SEARCH, NEWTON-RAPHSON STEP = ')
 9991 FORMAT(' IN DIR_SEARCH, HESSIAN INVERSE = ')
 9990 FORMAT(' IN DIR_SEARCH, REBUILT HESSIAN = ')
 9989 FORMAT(' IN DIR_SEARCH, -HESEDM-        = ')
 9988 FORMAT(' IN DIR_SEARCH, CURRENT GRADIENT ')
      END
      SUBROUTINE HND_OPT_PRINT(RTDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      INTEGER   RTDB
      PARAMETER (MXBFN=3072)
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL      DBUG
      LOGICAL      CVGED
      CHARACTER*16 ATMNAM
      CHARACTER*8  BFLAB 
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLLAB/ATMNAM(MXATOM),BFLAB(MXBFN)
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED
      COMMON/HND_OPTOPT/CVGOPT,NPTOPT,INHESS,LINOPT,MODUPD,MODDIR
      COMMON/HND_OPTFRZ/ICFRZ(MXCOOR),IATFRZ(MXATOM)
      COMMON/HND_OPTFUN/E,EG(MXCART)
C
      DBUG=.FALSE.
C
C     ----- PRINT OPTIMIZATION SUMMARY -----
C
      DELF=FUNC-FUNC0
      IF( GA_NODEID().EQ.0.OR.DBUG) THEN
         WRITE(IW,9999) NSERCH,NUPDAT,NPTS,FUNC,
     1                  GNORM,ALPHA,CVGED,FUNC0,DELF 
      ENDIF
      RETURN
 9999 FORMAT(1H1,/,1X,
     1       'nserch  nupdat   npts       func             ',
     2       'gnorm',8X,'  alpha   cvged',
     3       /,1X,I5,I8,I7,F17.8,F15.8,F14.5,L4,
     4       /,1X,'     previous energy',F17.8,
     5       /,1X,'  energy convergence',F17.8)
      END
      subroutine hnd_tfhs(hs,zs,t,ns,ms,ls)    
      implicit double precision (a-h,o-z)
      logical out
      logical dbug
      common/hnd_iofile/ir,iw
      dimension hs(ls,*),zs(ls,*),t(*)
      data zero   /0.0d+00/       
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(out) then
         write(iw,*) 'in -tfhs- '
         call hnd_prsq(hs,ms,ms,ls)
         call hnd_prsq(zs,ms,ms,ls)
      endif
c
      do i=1,ms
         do j=1,ns
            dum=zero
            do k=1,ms
               dum=dum+hs(i,k)*zs(k,j)
            enddo
            t(j)=dum
         enddo
         do j=1,ns
            hs(i,j)=t(j)
         enddo
      enddo
c
      do j=1,ns
         do i=1,ns
            dum=zero
            do k=1,ms
               dum=dum+zs(k,i)*hs(k,j)
            enddo
            t(i)=dum
         enddo
         do i=1,ns
            hs(i,j)=t(i)
         enddo
      enddo
c
      if(out) then
         call hnd_prsq(hs,ms,ms,ls)
      endif
c
      return
      end
      subroutine hnd_tfhsi(hs,zs,t,ns,ms,ls)    
      implicit double precision (a-h,o-z)
      logical out
      logical dbug
      common/hnd_iofile/ir,iw
      dimension hs(ls,*),zs(ls,*),t(*)
      data zero   /0.0d+00/       
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      do i=1,ns
         do j=1,ms
            dum=zero
            do k=1,ns
               dum=dum+hs(i,k)*zs(j,k)
            enddo
            t(j)=dum
         enddo
         do j=1,ms
            hs(i,j)=t(j)
         enddo
      enddo
c
      do j=1,ms
         do i=1,ms
            dum=zero
            do k=1,ns
               dum=dum+zs(i,k)*hs(k,j)
            enddo
            t(i)=dum
         enddo
         do i=1,ms
            hs(i,j)=t(i)
         enddo
      enddo
c
      return
      end
      logical function geom_hnd_tfgx(gs,gx,bi,ns,nx,mx)
      implicit double precision (a-h,o-z)
      logical out
      logical dbug
      common/hnd_iofile/ir,iw
      dimension gs(*),gx(*),bi(mx,*)
      data zero   /0.0d+00/       
      data tenm08 /1.0d-08/
c
      dbug=.false.
      out =.false.
c
      do j=1,ns
         dum=zero
         do i=1,nx
            dum=dum+gx(i)*bi(i,j)
         enddo
         gs(j)=dum
      enddo
c
      if(out) then
         write(iw,*) ' in -tfgx- gx = '
         do i=1,nx
            write(iw,9999) i,gx(i)
         enddo
         if(dbug) then
            call hnd_prsq(bi,ns,nx,nx)
         endif
         write(iw,*) ' in -tfgx- gs = '
         do j=1,ns
            write(iw,9999) j,gs(j)
         enddo
      endif
c
      do j=1,ns
         if(abs(gs(j)).lt.tenm08) gs(j)=zero
      enddo
c
      geom_hnd_tfgx=.true.
      return
 9999 format(i5,f15.10)
      end
      logical function geom_hnd_tfds(ds,dx,bi,ns,nx,mx)
      implicit double precision (a-h,o-z)
      logical out
      logical dbug
      common/hnd_iofile/ir,iw
      dimension ds(*),dx(*),bi(mx,*)
      data zero   /0.0d+00/       
      data tenm09 /1.0d-09/
c
      dbug=.false.
      out =.false.
c
      do i=1,nx
         dum=zero
         do j=1,ns
            dum=dum+ds(j)*bi(i,j)
         enddo
         dx(i)=dum
      enddo
c
      if(out) then
         write(iw,*) ' in -tfds- ds = '
         do i=1,ns
            write(iw,9999) i,ds(i)
         enddo
         if(dbug) then
            call hnd_prsq(bi,ns,nx,nx)
         endif
         write(iw,*) ' in -tfds- dx = '
         do j=1,nx
            write(iw,9999) j,dx(j)
         enddo
      endif
c
      do i=1,nx
         if(abs(dx(i)).lt.tenm09) dx(i)=zero
      enddo
c
      geom_hnd_tfds=.true.
      return
 9999 format(i5,f15.10)
      end
      logical function geom_hnd_tfdx(ds,dx,b,ns,nx,mx)
      implicit double precision (a-h,o-z)
      logical out
      logical dbug
      common/hnd_iofile/ir,iw
      dimension ds(*),dx(*),b(mx,*)
      data zero /0.0d+00/       
c
      dbug=.false.
      out =.false.
c
      do j=1,ns
         dum=zero
         do i=1,nx
            dum=dum+dx(i)*b(i,j)
         enddo
         ds(j)=dum
      enddo
c
      if(out) then
         write(iw,*) ' in -tfdx- dx = '
         do i=1,nx
            write(iw,9999) i,dx(i)
         enddo
         if(dbug) then
            call hnd_prsq(b,ns,nx,nx)
         endif
         write(iw,*) ' in -tfdx- ds = '
         do j=1,ns
            write(iw,9999) j,ds(j)
         enddo
      endif
c
      geom_hnd_tfdx=.true.
      return
 9999 format(i5,f15.10)
      end
      SUBROUTINE HND_OPT_SEARCH_LIN(RTDB,GEOM,RSTART,EXCESS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "context.fh"
C
C     ----- ONE DIMENSIONAL SEARCH.  -----
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
C
      INTEGER RTDB
      INTEGER GEOM
      LOGICAL STATUS
      CHARACTER*16 TAGS_NW
      DIMENSION COORDS_NW(MXCART)
      DIMENSION CHARGE_NW(MXATOM)
      DIMENSION   TAGS_NW(MXATOM)
C
      LOGICAL RSTART
      LOGICAL EXCESS
      LOGICAL LOWER,ILLCON,RECALC,RETRY
      LOGICAL CVGED
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/RESTAR/NREST
      COMMON/HND_MOLXYZ/C(MXCART),ZAN(MXATOM),NAT         
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED 
      COMMON/HND_OPTOPT/CVGOPT,NPTOPT,INHESS,LINOPT,MODUPD,MODDIR
      COMMON/HND_OPTFRZ/ICFRZ(MXCOOR),IATFRZ(MXATOM)
      COMMON/HND_OPTFUN/E,EG(MXCART)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      DIMENSION ALP(5),FUNL(5)
      DATA ZERO,TWO,THREE       /0.0D+00,2.0D+00,3.0D+00/
      DATA TOLALP,TOLFUN        /1.0D-01,1.5D-01/
      DATA FALMIN,FALNOR,FALMAX /1.0D-01,2.0D+00,5.0D+00/
      DATA DALNOR,DALMAX        /0.9D+00,0.5D+00/
C
      DBUG=.FALSE.
      OUT =.TRUE.
      OUT =OUT.OR.DBUG
      IF(OUT) THEN
         WRITE(IW,9994)
         WRITE(IW,9992) ALPHA
      ENDIF
C
      LOWER =.FALSE.
      ILLCON=.FALSE.
      RECALC=.FALSE.
      RETRY =.FALSE.
C
      NCART =3*NAT
C
C     ----- NORMAL START -----
C
      ALPH = ALPHA
  100 CONTINUE
      FUNC0 = FUNC
      GMAX0 = GMAX 
      DO I = 1,NCART 
         X0(I) = X(I)
         G0(I) =EG(I)
         G (I) =EG(I)
      ENDDO
      DUM=ZERO
      DO I=1,NCART 
         DUM=DUM+G0(I)*DX(I)
      ENDDO
      GS0=DUM
C
           NLPTS  = 1
       ALP(NLPTS) = ZERO
      FUNL(NLPTS) = FUNC0
C
C     ----- SEARCH ALONG -DX-DIRECTION. -----
C
  200 CONTINUE
      IF(OUT) THEN
         WRITE(IW,9991) ALPH
      ENDIF
      DO I = 1,NCART 
         X(I) = X0(I)+ALPH*DX(I)
         C(I) = X(I)
      ENDDO
      NPTS=NPTS+1
      IF(NPTS.GT.NPTOPT) THEN
         EXCESS=.TRUE.
         RETURN
      ENDIF
C
C     ----- WRITE TO -NWCHEM- -----
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF(GA_NODEID().EQ.0) THEN
         STATUS=GEOM_CREATE(GEOM,'geometry')
         STATUS=GEOM_RTDB_LOAD(RTDB,GEOM,'geometry')
         STATUS=GEOM_CART_GET(GEOM,NAT_NW,TAGS_NW,COORDS_NW,
     1                                            CHARGE_NW)
         CALL GA_BRDCST(30+MSGINT,  NAT_NW ,MITOB(1    ),0)
         CALL GA_BRDCST(31+MSGDBL,  TAGS_NW,   16*NCART ,0)
         CALL GA_BRDCST(32+MSGDBL,COORDS_NW,MDTOB(NCART),0)
         CALL GA_BRDCST(33+MSGDBL,CHARGE_NW,MDTOB(NCART),0)
      ELSE
         CALL GA_BRDCST(30+MSGINT,  NAT_NW ,MITOB(1    ),0)
         CALL GA_BRDCST(31+MSGDBL,  TAGS_NW,   16*NCART ,0)
         CALL GA_BRDCST(32+MSGDBL,COORDS_NW,MDTOB(NCART),0)
         CALL GA_BRDCST(33+MSGDBL,CHARGE_NW,MDTOB(NCART),0)
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
      IF(DBUG) THEN
         WRITE(IW,9997)
         WRITE(IW,9998) (       X0(I),I=1,NCART)
         WRITE(IW,9996)
         WRITE(IW,9998) (       DX(I),I=1,NCART)
         WRITE(IW,9995) ALPH
         WRITE(IW,9998) (        X(I),I=1,NCART)
         WRITE(IW,9999)
         WRITE(IW,9998) (COORDS_NW(I),I=1,NCART)
      ENDIF
      DO I = 1,NCART
         COORDS_NW(I)=X(I)
      ENDDO
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF(GA_NODEID().EQ.0) THEN
         STATUS=GEOM_CART_SET(GEOM,NAT_NW,TAGS_NW,COORDS_NW,
     1                                            CHARGE_NW)
         STATUS=GEOM_RTDB_STORE(RTDB,GEOM,'geometry')
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
      IF(OUT) THEN
         WRITE(IW,9999)
         WRITE(IW,9998) (COORDS_NW(I),I=1,NCART)
      ENDIF
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF(GA_NODEID().EQ.0) THEN
         STATUS=GEOM_DESTROY(GEOM)
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
C     ----- CALL FUNCTION EVALUATION -----
C
      CALL HND_OPT_ENERGY(RTDB)
      FUNC=E
C
      RSTART = .FALSE.
      IF(RECALC) THEN
         GO TO 900
      ENDIF
C
C     ----- SAVE UP TO THE LAST FOUR STEPS -----
C
          NLPTS  = NLPTS+1
           NPOS  = MIN0(4,NLPTS)
       ALP(NPOS) = ALPH
      FUNL(NPOS) = FUNC
      LOWER      = FUNC.LT.FUNC0
C
C     ----- ORDER THE FUNCTION VALUES AND THEIR ALPHAS -----
C
      NLPOS = NPOS
      DO I = 2,NPOS
         J = NPOS-I+1
         IF(FUNC.GT.FUNL(J)) THEN
            GO TO 400
         ELSE
            FUNL(J+1) = FUNL(J)
            FUNL(J  ) = FUNC
             ALP(J+1) = ALP(J)
             ALP(J  ) = ALPH
            NLPOS     = J
         ENDIF
      ENDDO    
C
  400 IF(NLPTS.GT.2) GO TO 440
C
C     ----- FIT A PARABOLA -----      
C
      AQUAD = (FUNC-GS0*ALPH-FUNC0)/(ALPH*ALPH)
      BQUAD = GS0
      CQUAD = FUNC0
      IF(AQUAD.LT.ZERO) THEN
         IF(LOWER) THEN
            GO TO 900
         ELSE
            IF(ILLCON) THEN
               GO TO 900
            ELSE
               ILLCON=.TRUE.
               ALPH=ALPH*FALMIN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF
      ALPHI = -BQUAD/(TWO*AQUAD)
      FESTI =  AQUAD*ALPHI**2+BQUAD*ALPHI+CQUAD
C
      IF(LOWER) THEN
         ALPHN = DALNOR*(ALPHI-ALPH)+ALPH
      ELSE
         ALPHN = ALPHI
      ENDIF
      IF(ALPHN.GT.FALNOR*ALPH) THEN
         ALPHN = DMIN1(FALMAX*ALPH,
     1                 DALMAX*(ALPHN-FALNOR*ALPH)+FALNOR*ALPH)
      ENDIF
      IF(ALPHN.LT.ZERO) THEN
         IF(LOWER) THEN
            GO TO 900
         ELSE
            IF(ILLCON) THEN
               GO TO 900
            ELSE
               ILLCON=.TRUE.
               ALPH=ALPH*FALMIN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF
      FEST = AQUAD*ALPHN**2+BQUAD*ALPHN+CQUAD
      IF(FEST .GT. FUNL(1)) THEN
         IF(LOWER) THEN
            GO TO 900
         ELSE
            IF(ILLCON) THEN
               GO TO 900
            ELSE
               ILLCON=.TRUE.
               ALPH=ALPH*FALMIN
               GO TO 200
            ENDIF
         ENDIF
      ENDIF
      IF(.NOT.LOWER) THEN
         ALPH = ALPHN
         GO TO 200
      ENDIF
      IF(ABS((ALPHN-ALPH)/ALPH).LT.TOLALP) THEN
         GO TO 900
      ENDIF
      IF(ABS((FEST-FUNL(1))/(FUNC0-FUNL(1))).LT.TOLFUN) THEN
         GO TO 900
      ENDIF
      ALPH = ALPHN
      GO TO 200
C
C     --- CALCULATE THE PARABOLIC MINIMUM WITH THE LOWEST THREE POINTS ---
C
  440 IF(.NOT.LOWER) THEN
         GO TO 500
      ENDIF
      IF(NLPOS.NE.1) THEN
         GO TO 600
      ENDIF
      ALDIF1 = ALP(2)-ALP(3)
      ALDIF2 = ALP(3)-ALP(1)
      ALDIF3 = ALP(1)-ALP(2)
      AQUAD  = -(FUNL(1)*ALDIF1+FUNL(2)*ALDIF2+FUNL(3)*ALDIF3)
     1         /(ALDIF1*ALDIF2*ALDIF3)
      BQUAD = (FUNL(1)-FUNL(2))/ALDIF3-AQUAD*(ALP(1)+ALP(2))
      CQUAD = FUNL(1)-ALP(1)*(ALP(1)*AQUAD+BQUAD)
      IF(AQUAD.LT.ZERO) THEN
         GO TO 900
      ENDIF
C
      ALPHI = -BQUAD/(TWO*AQUAD)
      FESTI = (AQUAD*ALPHI+BQUAD)*ALPHI+CQUAD
      ALPHN = DALNOR*(ALPHI-ALPH)+ALPH
      IF(ALPHN.GT.FALNOR*ALPH) THEN
         ALPHN = DMIN1(FALMAX*ALPH,
     1                 DALMAX*(ALPHN-FALNOR*ALPH)+FALNOR*ALPH)
      ENDIF
      IF(ALPHN.LT.ZERO) THEN
         GO TO 900
      ENDIF
      FEST = (AQUAD*ALPHN+BQUAD)*ALPHN+CQUAD
      IF(FEST.GT.FUNL(1)) THEN
         GO TO 900
      ENDIF
      IF(ABS((ALPHN-ALPH)/ALPH).LT.TOLALP) THEN
         GO TO 900
      ENDIF
      IF(ABS((FEST-FUNL(1))/(FUNC0-FUNL(1))).LT.TOLFUN*NLPTS) THEN
         GO TO 900
      ENDIF
      ALPH = ALPHN
      GO TO 200
C
C     ----- IF WE ALREADY HAVE A LOWER ENERGY, TAKE IT -----
C
  480 IF(LOWER) THEN
         GO TO 900
      ENDIF
C
C     ----- FUNCTION IS NOT PARABOLIC ALONG SEARCH DIRECTION -----
C
  500 IF(ILLCON) THEN
         GO TO 900
      ENDIF
      ILLCON = .TRUE.
      ALPH   = ALPH*FALMIN
      GO TO 200
C
C     ----- THE LAST POINT MUST BE RECALCULATED -----
C
  600 RECALC = .TRUE.
      ALPH   = ALP(1)
      GO TO 200
C
C     ----- CALCULATE THE GRADIENT FOR THE FINAL POINT -----
C
  900 CONTINUE
      CALL HND_OPT_GRADIENT(RTDB)
      DO I=1,NCART 
         G(I) =EG(I)
      ENDDO
C
      RSTART = .FALSE.
C
      ALPHA  = ALPH
C
      DUM=ZERO
      DO I=1,NCART
         DUM=DUM+G(I)*G(I)
      ENDDO
      GNORM= SQRT(DUM)
      DUM=ZERO
      DO I=1,NCART
         DUM=DUM+G(I)*DX(I)
      ENDDO
      GS=DUM
      IF((GS.GT.GS0).OR.RETRY) THEN      
         GO TO 1000
      ENDIF
      RETRY = .TRUE.
      ACUBE = (-TWO*(FUNC-FUNC0)+(GS+GS0)*ALPH)/ALPH**3
      BCUBE = -ACUBE*ALPH+(FUNC-FUNC0-GS0*ALPH)/ALPH**2
      DSCRM =  BCUBE*BCUBE-THREE*ACUBE*GS0
      IF(DSCRM.LT.ZERO) THEN
         RETURN
      ENDIF
      ALPH  = (-BCUBE+ SQRT(DSCRM))/(THREE*ACUBE)
      IF((ALPH.LT.ALPHA).AND.(ALPH.GT.ZERO)) THEN
         GO TO 100
      ENDIF
C
C     ----- RETURN AFTER SUCCESSFUL SEARCH -----
C
 1000 CONTINUE
      IF(DBUG) THEN
         WRITE(IW,9993)
      ENDIF
      RETURN
 9999 FORMAT(' IN SEARCH_LIN, COORDS = ')
 9998 FORMAT(F12.7)
 9997 FORMAT(' IN SEARCH_LIN, X0 = ')
 9996 FORMAT(' IN SEARCH_LIN, DX = ')
 9995 FORMAT(' IN SEARCH_LIN, ALPH = ',F8.4,' X = ')
 9994 FORMAT(' STARTING SEARCH_LIN ')
 9993 FORMAT('   ENDING SEARCH_LIN ')
 9992 FORMAT(' IN SEARCH_LIN, ALPHA = ',F8.4)
 9991 FORMAT(' IN SEARCH_LIN, ALPH  = ',F8.4)
      END
      SUBROUTINE HND_OPT_HSS_UPDATE(RTDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "geomP.fh"
#include "rtdb.fh"
#include "context.fh"
      INTEGER RTDB
      LOGICAL STATUS
      LOGICAL DRIVER_OPT_PUT_HESS
      LOGICAL DRIVER_OPT_GET_HESS
      LOGICAL GEOM_ZMT_GET_ZSYM
      DIMENSION X(1)
      EQUIVALENCE (X(1),DBL_MB(1))
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL   DBUG
      LOGICAL   OUT
      LOGICAL   ZCOORD
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_OPTVAR/ZCOORD,NCOORD,MCOORD
      DATA ZERO /0.0D+00/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(out) then
         write(iw,*) 'start of hnd_opt_hss_update'
         write(iw,*) 'nat, nzvar, ncoord, mcoord = ',
     1                nat, nzvar, ncoord, mcoord
      endif
c
      ndimx=3*nat
      ndimq=nzvar
      ndimi=ncoord
      ndim =mcoord
c
c     ----- get memory -----
c
      need=      ndimi*ndimi
      need=need +ndimi
      if(zcoord) then
         need=need +ndim *ndim
         need=need +ndim *ndim
         need=need +ndim
      endif
      status=ma_push_get(mt_dbl,need,'mem update',i_ifirst,ifirst)
      i10=ifirst   
      i20=i10+ndimi*ndimi
      i30=i20+ndimi
      if(zcoord) then
         j10=i30
         j20=j10+ndim *ndim
         j30=j20+ndim *ndim
         j40=j30+ndim
      endif
c
c     ----- get current hessian matrix -----
c
      if(zcoord) then
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_get_hess(x(j10),ndimq,ndim) 
            status=geom_zmt_get_zsym(x(j20),ndimq,ndimq,ndim)
            call ga_brdcst(35+msgdbl,x(j10),mdtob(ndim*ndim),0)
            call ga_brdcst(36+msgdbl,x(j20),mdtob(ndim*ndim),0)
         else
            call ga_brdcst(35+msgdbl,x(j10),mdtob(ndim*ndim),0)
            call ga_brdcst(36+msgdbl,x(j20),mdtob(ndim*ndim),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.)
c
         if(dbug) then
            write(iw,*) 'in hss_update, current hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
            write(iw,*) 'in hss_update, zsym = '
            call hnd_prsq(x(j20),ndimq,ndimq,ndim)
         endif
         call hnd_tfhs(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
         do j=1,ndimi
            do i=1,ndimi
               x(i+ndimi*(j-1)+i10-1)=x(i+ndim*(j-1)+j10-1)
            enddo
         enddo
         if(out) then
            write(iw,*) 'in hss_update, transformed hessian = '
            call hnd_prsq(x(i10),ndimi,ndimi,ndimi)
         endif
      else
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_get_hess(x(i10),ndimx,ndim)
            call ga_brdcst(37+msgdbl,x(i10),mdtob(ndim*ndim),0)
         else
            call ga_brdcst(37+msgdbl,x(i10),mdtob(ndim*ndim),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.)
         if(dbug) then
         if( ga_nodeid().eq.0) then
            write(iw,*) 'in hss_update, cartesian hessian = '
            call hnd_prsq(x(i10),ndimx,ndimx,ndim)
         endif
         call ga_sync()
         if( ga_nodeid().eq.1) then
            write(iw,*) 'in hss_update, cartesian hessian = '
            call hnd_prsq(x(i10),ndimx,ndimx,ndim)
         endif
         call ga_sync()
         endif
      endif
c
c     ----- update hessian matrix -----
c
      call hnd_opt_hss_update_x(x(i10),x(i20),ndimi,ndimi) 
c
c     ----- transform back updated hessian and store -----
c
      if(zcoord) then
         if(dbug) then
            write(iw,*) 'in hss_update, updated hessian = '
            call hnd_prsq(x(i10),ndimi,ndimi,ndimi)
         endif
         do j=1,ndimq
            do i=1,ndimq
               x(i+ndim*(j-1)+j10-1)=zero
            enddo
         enddo
         do j=1,ndimi
            do i=1,ndimi
               x(i+ndim*(j-1)+j10-1)=x(i+ndimi*(j-1)+i10-1)
            enddo
         enddo
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=geom_zmt_get_zsym(x(j20),ndimq,ndimq,ndim)
            call ga_brdcst(38+msgdbl,x(j20),mdtob(ndim),0)
         else
            call ga_brdcst(38+msgdbl,x(j20),mdtob(ndim),0)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.)
c
         call hnd_tfhsi(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
c
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_put_hess(x(j10),ndimq,ndim)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.)
         if(out) then
            write(iw,*) 'in hss_update, transformed updated hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
         endif
      else
         status=rtdb_parallel(.false.)
         if( ga_nodeid().eq.0) then
            status=driver_opt_put_hess(x(i10),ndimx,ndim)
         endif
         call ga_sync()
         status=rtdb_parallel(.true.)
         if(dbug) then
            write(iw,*) 'in hss_update, updated hessian = '
            call hnd_prsq(x(i10),ndimx,ndimx,ndim)
         endif
      endif
c
c     ----- release memory -----
c
      status=ma_pop_stack(i_ifirst)
c
      if(dbug) then
         write(iw,*) '  end of hnd_opt_hss_update'
      endif
c
      RETURN
      END
      SUBROUTINE HND_OPT_HSS_UPDATE_X(HESS,HESSDX,NVAR,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL DBUG
      LOGICAL OUT 
      LOGICAL CVGED
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED
      COMMON/HND_OPTOPT/CVGOPT,NPTOPT,INHESS,LINOPT,MODUPD,MODDIR
      COMMON/HND_OPTFRZ/ICFRZ(MXCOOR),IATFRZ(MXATOM)
      DIMENSION HESS(NDIM,*),HESSDX(*)
      DATA ZERO /0.0D+00/
      DATA TOL  /1.0D-08/
C
      DBUG=.FALSE.
      OUT =.TRUE. 
      OUT =OUT.OR.DBUG
      IF(DBUG) THEN
         WRITE(IW,9999)
         CALL HND_PRSQ(HESS,NVAR,NVAR,NDIM)
      ELSEIF(OUT) THEN
         WRITE(IW,9999)
      ENDIF
      IF(OUT) THEN
         WRITE(IW,9997) 
         WRITE(IW,9994) (G0(I),I=1,NVAR)
         WRITE(IW,9996) 
         WRITE(IW,9994) ( G(I),I=1,NVAR)
         WRITE(IW,9995) 
         WRITE(IW,9994) (DX(I),I=1,NVAR)
      ENDIF
C
C     ----- DAMP OUT FROZEN COORDINATES -----
C
      DO I=1,NVAR
         IF(ICFRZ(I).NE.0) THEN      
            G(I)=G0(I)
            DO J=1,NVAR
               HESS(I,J)=ZERO
               HESS(J,I)=ZERO
            ENDDO
         ENDIF
      ENDDO    
C
      DO I=1,NVAR
         DUM=ZERO
         DO J=1,NVAR
            DUM=DUM+HESS(I,J)*DX(J)
         ENDDO
         HESSDX(I)=DUM
      ENDDO    
      DUMG=ZERO
      DUMX=ZERO
      DO I=1,NVAR
         DUMX=DUMX+DX(I)*DX(I)
         DUMG=DUMG+DX(I)*(G(I)-G0(I))
      ENDDO
      DXDG=DUMG
      DXDX=DUMX
      DUM=ZERO
      DO I=1,NVAR
         DUM=DUM+DX(I)*HESSDX(I)
      ENDDO
      DXHDX=DUM
      IF(OUT) THEN
         WRITE(IW,*) 'DXDG  = ',DXDG
         WRITE(IW,*) 'DXDX  = ',DXDX
         WRITE(IW,*) 'DXHDX = ',DXHDX
      ENDIF
C
      IF(MODUPD.LE.1) THEN       
C
C     ----- -BFGS- UPDATE -----
C
         IF(ABS(DXDG).GT.TOL) THEN
            DO I=1,NVAR
               DO J=1,NVAR
                  HESS(I,J)=HESS(I,J)+(G(I)-G0(I))*(G(J)-G0(J))/DXDG
     1                               -  HESSDX(I) *  HESSDX(J) /DXHDX
               ENDDO   
            ENDDO   
         ENDIF
      ELSE      
C
C     ----- -PSB- UPDATE -----
C
         DO I=1,NVAR
            DO J=1,NVAR
               HESS(I,J)=HESS(I,J)+((G(I)-G0(I))-HESSDX(I))*DX(J)/DXDX
     1                            +((G(J)-G0(J))-HESSDX(J))*DX(I)/DXDX
     2                            -DX(I)*DX(J)*(DXDG-DXHDX)/(DXDX*DXDX)
            ENDDO    
         ENDDO    
      ENDIF    
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         CALL HND_PRSQ(HESS,NVAR,NVAR,NDIM)
      ENDIF
      RETURN
 9999 FORMAT(' STARTING HSS_UPDATE ')
 9998 FORMAT('   ENDING HSS_UPDATE ')
 9997 FORMAT('       IN HSS_UPDATE, G0 = ')
 9996 FORMAT('       IN HSS_UPDATE, G  = ')
 9995 FORMAT('       IN HSS_UPDATE, DX = ')
 9994 FORMAT(F10.5)
      END
