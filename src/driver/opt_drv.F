      logical function drv_opt(rtdb)
C$Id: opt_drv.F,v 1.66 1999-05-30 21:50:28 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer  rtdb
      integer  idriver
      logical  hnd_optx
      external hnd_optx
      integer  ir, iw
      logical  some, dbug
      common/hnd_iofile/ir,iw
c
      call util_print_push
      call util_print_rtdb_load(rtdb, 'driver')
      call ecce_print_module_entry('driver')
      if (util_print('brdcst', print_never)) call setdbg(1)
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
c     ----- optimization -----
c
      idriver=0
c
      drv_opt = hnd_optx(rtdb,idriver)
c
      if (drv_opt) then
         call ecce_print_module_exit('driver', 'ok')
      else
         call ecce_print_module_exit('driver', 'failed')
      endif
c
      call movecs_ecce_print_on() ! Restore MO printing
      if (util_print('brdcst', print_never)) call setdbg(0)
      call util_print_pop
c
      end
      logical function hnd_optx(rtdb,idriver)
      implicit double precision (a-h,o-z)
#include "rtdb.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "geom.fh"
#include "util.fh"
#include "coptopt.fh"
#include "inp.fh"
      logical       geom_zmt_get_nizmat
      logical       geom_zmt_get_nzfrz 
      logical       geom_zmt_get_izfrz 
      integer       rtdb
      integer       geom
      integer       idriver
      logical       hnd_sadx
      logical       some
      logical       out
      logical       dbug
      logical       cvged
      logical       zcoord
      logical       rstart
      logical       excess
      logical       zmatrix     ! True if user specified zmatrix
      character*8   source
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c     mxcoor=max(mxcart,mxzmat)
c     
      logical geom_zmt_get_varinfo, geom_zmt_get_source, 
     $     geom_zmt_get_nzvar
      character*8 zvarname(mxcoor),test
      double precision zvarsign(mxcoor)

      common/hnd_iofile/ir,iw
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1     g0(mxcoor),g(mxcoor),ds(mxcoor),
     2     func,func0,gmax,gmax0,curv,alpha,gnorm
      common/hnd_optrun/npts,nserch,nupdat,cvged
      common/hnd_optsad/sadstp
      common/hnd_optfrz/nzfrz,izfrz(mxcoor),iatfrz(mxatom)
      common/hnd_optfun/e,eg(mxcart)
      common/hnd_optvar/zcoord,ncoord,mcoord
      common/hnd_opttim/energy_time,gradient_time
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(mxcart),zan(mxatom),nat
      character*16 atmlab
      common/hnd_mollab/atmlab(mxatom)
      dimension t0(mxcoor)
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data alphmn /0.4d+00/
      data alphmx /2.5d+00/
      data four   /4.0d+00/
c     
      hnd_optx = .false.
c     
c     Initialize HONDO common blocks with junk for paranoia
c     since it was relying on (node 0 only for last geometry read!) 
c     passing data from the input routines ... explains some modes 
c     of failure that were hard to understand.
c
      call dfill(6*mxcoor+7,1d300,x0,1)
      call ifill(4,-99999999,npts,1)
      nzfrz = 9999
      ncoord = -99999
      mcoord = -99099
      nzmat = 98888888
      nzvar = 89898989
      nvar = 81818181
      sadstp = 1d300
      call dfill(mxcart+mxatom,1d300,c,1)
      nat = 99999999
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.
      some=some.or.out
      if ((ga_nodeid().eq.0.or.some) .and.
     $     util_print('information', print_low))  then
         write(iw,9999)
      endif
c     
c     ----- default parameters for optimization driver -----
c     
      do i=1,mxcoor
         izfrz(i)=0
      enddo
      do i=1,mxatom
         iatfrz(i)=0
      enddo
c     
c     ----- parameters for optimization driver -----
c     
      if (.not. rtdb_get(rtdb,'driver:cvgopt',mt_dbl,1,cvgopt))
     $     cvgopt=8.0d-04
      if (.not. rtdb_get(rtdb,'driver:nptopt',mt_int,1,nptopt))
     $     nptopt=20
      if (.not. rtdb_get(rtdb,'driver:inhess',mt_int,1,inhess))
     $     inhess=0
      if (.not. rtdb_get(rtdb,'driver:linopt',mt_int,1,linopt)) 
     $     linopt=10
      if (.not. rtdb_get(rtdb,'driver:modupd',mt_int,1,modupd))
     $     modupd=1
      if (.not. rtdb_get(rtdb,'driver:moddir',mt_int,1,moddir)) 
     $     moddir=1
      if (.not. rtdb_get(rtdb,'driver:modsad',mt_int,1,modsad)) 
     $     modsad=0
      if (.not. rtdb_get(rtdb,'driver:sadstp',mt_dbl,1,sadstp)) 
     $     sadstp=0.35d+00
      if (.not. rtdb_get(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $     oqstep = .true.
c     
c     ----- get going ... -----
c     
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('hnd_opt: geom_create?', 911)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('hnd_opt: no geometry ', 911)
      if (.not. geom_cart_get(geom, ncent, atmlab, c, zan))
     &     call errquit('hnd_opt: geom_get ', 911)
      nat = ncent
c     
c     Zmatrix info
c
      if (.not.geom_zmt_get_source(geom,source))
     $     call errquit('hnd_optx: zmat source ',geom)
      zcoord  = (source .ne. ' ')    ! Any kind of zmatrix
      zmatrix = (source .eq. 'user') ! zmatrix from user
      if (zcoord) then
         if(.not.geom_zmt_get_nizmat(geom,nzmat))
     &        call errquit('geom_input: geom_zmt_get_nizmat failed',0)
         if(.not.geom_zmt_get_nzvar(geom,nzvar))
     &        call errquit('geom_input: geom_zmt_get_nzvar failed',0)
c     This next call fills in nvar, nizmat and other bits of HONDO common
         call geom_zmt_icoord(rtdb,geom)
      endif
c
      if (.not. zmatrix) oqstep = .false.
      if (nzvar .ne. nvar) oqstep = .false.
c     
      if (util_print('information', print_low) .and. 
     $     (ga_nodeid().eq.0.or.dbug))  then
         write(iw,9996) cvgopt,nptopt,linopt,inhess,modupd,
     1        modsad,moddir,
     2        sadstp, source, oqstep
      endif
c     
      nzfrz=0
      do i=1,mxcoor           
         izfrz(i)=0
      enddo
      do i=1,mxatom
         iatfrz(i)=0
      enddo
c     
c     ----- branch to transition state search if demanded -----
c     
      if(idriver.ne.0.or.modsad.ne.0) then
         hnd_optx=hnd_sadx(rtdb,idriver)
         return
      else
         if ((ga_nodeid().eq.0.or.some) .and.
     $        util_print('information', print_low))  then
            write(iw,9994)
         endif
      endif
c     
c     ----- frozen coordinates ? if so .... -----
c     
      if(.not.geom_zmt_get_nzfrz(geom,nzfrz)) then
         nzfrz=0
      else
         if(.not.geom_zmt_get_izfrz(geom,izfrz,nzfrz)) then
            call errquit('geom_input: geom_zmt_get_izfrz  failed',0)
         endif
      endif
c     
      if ((ga_nodeid().eq.0.or.dbug) .and.
     $     util_print('information', print_low)) then
         if(nzfrz.gt.0) then
            write(iw,9993) nzfrz
            write(iw,9992) (izfrz(i),i=1,nzfrz) 
            if (.not. oqstep) write(6,*) 'freezing of variables',
     $           ' will not be exact since qstep is not used'
         else
            write(iw,9991)
         endif
         if (zmatrix) then
            if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $           call errquit('hnd_optx: varinfo?',0)
            write(6,*)
            write(6,*) '  names of Z-matrix variables '
            do i = 1, nzvar, 5
               write(6,5432) (j,zvarname(j),j=i,min(nzvar,i+4))
 5432          format(1x,5(i4,1x,a8,2x))
            enddo
            write(6,*)
            write(6,*) 'variables with the same non-blank name',
     $           ' are constrained to be equal'
            if (.not. oqstep) write(6,*) 'constraints will be not be',
     $           ' exact since qstep is not used'
            write(6,*)
            do i = 1, nzvar
               if (zvarname(i).ne.' ') then
                  num = 0
                  test = zvarname(i)
                  do j = 1, nzvar
                     if (zvarname(j).eq.test) then 
                        num = num + 1
                     endif
                  enddo
                  if (num .gt. 1) then
                     write(6,5433) test(1:inp_strlen(test))
 5433                format(1x,'    variables named ',a8,' = ',$)
                     do j = 1, nzvar
                        if (zvarname(j).eq.test) then 
                           if (zvarsign(j).gt.0d0) then
                              write(6,5434) j
                           else
                              write(6,5434) -j
                           endif
 5434                      format(i5,$)
                           zvarname(j) = ' '
                        endif
                     enddo
                     write(6,*)
                  endif
               endif
            enddo
         endif
      endif
c     
c     ----- get ready to start now -----
c     if internal coordinates ...
c     nzvar  = # of (redundant) internal coordinates
c     ncoord = # of independent internal coordinates
c     mcoord = max( 3*nat , nzvar )
c     if cartesian coordinates ...
c     ncoord = 3*nat
c     mcoord = 3*nat
c     
      rstart=.false.
      excess=.false.
      if(zcoord) then
         ncart =3*nat
         ncoord=nvar                  
         mcoord=max(ncart,nzvar)
      else
         nzmat =0
         nzvar =0
         nvar  =3*nat
         ncart =3*nat
         ncoord=ncart
         mcoord=ncart
      endif
c     
c     ----- initialize hessian matrix -----
c     
      call hnd_opt_hss_init(rtdb,geom)
      if(zcoord) then
         if(inhess.eq.1) then
            call hnd_opt_tfhx(rtdb,geom)
         endif
      endif
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
c     ----- energy -----
c     
      time_start=util_wallsec()
      call hnd_opt_energy(rtdb,geom)
      time_end  =util_wallsec()
      energy_time=time_end-time_start
      if(out) then
         write(iw,*) 'energy_time = ',energy_time
      endif
c     
      func=e
c     
c     ----- gradient -----
c     
      time_start=util_wallsec()
      call hnd_opt_gradient(rtdb,geom)
      time_end  =util_wallsec()
      gradient_time=time_end-time_start
      if(out) then
         write(iw,*) 'gradient_time = ',gradient_time
      endif
c     
      do i=1,ncart
         g(i)=eg(i)
      enddo
      call ga_sync()
c
      if(zcoord) then
         call hnd_opt_tfgx(rtdb,geom)
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
      if (cvged) goto 3301
c     
 100  continue
c     
c     ----- search direction -----
c     
      nserch=nserch+1
      call hnd_opt_search_dir(rtdb,geom)
      do i=1,ncoord
         ds(i)=dx(i)
      enddo
c     
      if(alpha.eq.zero) then
         if(zcoord) then
            if(abs(gmax).lt.pt5) then
               alpha=one
            else
               alpha=one/four
            endif
         else
            if(abs(gmax).lt.pt5) then
               alpha=pt5
            else
               alpha=pt5/four
            endif
         endif
      endif
c     
      if(zcoord) call hnd_opt_tfds(rtdb,geom)               
      if (oqstep) then
         call hnd_opt_search_lin_int(rtdb,geom,rstart,excess)
      else
         call hnd_opt_search_lin(rtdb,geom,rstart,excess)
      endif
c     
      if(excess) then
         if (ga_nodeid().eq.0.or.dbug)  then
            write(iw,9995) nptopt
         endif
         if (.not.geom_destroy(geom))
     &        call errquit('hnd_opt: geom_destroy?', 911)
         call ga_sync()
         hnd_optx=.false.
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
         call hnd_opt_tfgx(rtdb,geom)
      endif
      nupdat=nupdat+1
      call hnd_opt_gmax(rtdb)
      call hnd_opt_converge(rtdb,cvged)
      call hnd_opt_print(rtdb)
c
      if(cvged) goto 3301
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
c     Converged ... print info and clean up
c     
 3301 if ((ga_nodeid().eq.0.or.dbug)
     $     .and. util_print('geometry',print_low))  then
         write(iw,9998)
         if (.not. geom_print(geom)) call errquit
     $        ('hnd_opt_drv: geom_print?',0)
      endif
      if ((dbug) .and. util_print('geometry',print_low))  then
         if (zcoord)  then 
            call hnd_zmtprt(geom) ! Must be called in parallel
            write(6,*)
         endif
      endif
      if ((ga_nodeid().eq.0.or.dbug)
     $     .and. util_print('geometry',print_low))  then
         if (util_print('bonds',print_default)) then
            if (.not.geom_print_distances(geom)) call errquit(
     &           'hnd_opt_drv: geom_print_distances failed',911)
         endif
         if (util_print('angles',print_default)) then
            if (.not.geom_print_angles(geom)) call errquit(
     &           'hnd_opt_drv: geom_print_angles failed',911)
         endif
      endif
      if (.not.geom_destroy(geom))
     &     call errquit('hnd_opt: geom_destroy?', 911)

      call ga_sync()
      hnd_optx=.true.
      return
c     
 9999 format(/,10x,21(1h-),
     1     /,10x,'Geometry Optimization',
     2     /,10x,21(1h-))
 9998 format(/,10x,22(1h-),
     1     /,10x,'Optimization converged',
     2     /,10x,22(1h-))
 9997 format(' in hnd_optx, nzvar,nvar,ncart,ncoord,mcoord = ',5i5)
 9996 format(' optimization convergence threshold    (cvgopt) = ',f10.6,
     1     /,' maximum number of energy evaluations  (nptopt) = ',i4,/,
     2     ' max.num. of points in linear search   (linopt) = ',i4,/,
     3     ' initial guess option for opt.hessian  (inhess) = ',i4,/,
     4     ' hessian update selection option       (modupd) = ',i4,/,
     5     ' optimization to saddle point          (modsad) = ',i4,/,
     6     ' initial eigen-mode followed to saddle (moddir) = ',i4,/,
     7     ' maximum step size to saddle           (sadstp) = ',f6.2/
     7     ' source of zmatrix                              =   ',a8/
     8     ' line search in internals              (oqstep) =    ',l1)
 9995 format(' maximum number of energy points reached. stop',
     1     ' currently -nptopt- = ',i5)
 9994 format(/,10x,19(1h-),
     1       /,10x,'Energy Minimization',
     2     /,10x,19(1h-),/)
 9993 format(' number of frozen internal coordinates = ',i4)
 9992 format('           frozen internal coordinates = ',12i4)
 9991 format('        no frozen internal coordinates.')
      end
      subroutine driver_input(rtdb)
      implicit none
#include "inp.fh"
#include "coptopt.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c                mxcoor=max(mxcart,mxzmat)
      character*255 field
c
      double precision sadstp
      common/hnd_optsad/sadstp
      integer nzfrz,izfrz,iatfrz
      common/hnd_optfrz/nzfrz,izfrz(mxcoor),iatfrz(mxatom)
c
c     So that task knows that driver is in control
c
      if (.not. rtdb_put(rtdb,'opt:driver',mt_log,1,.true.))
     $     call errquit('driver_input: rtdb_put failed',0)
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('driver_input: unexpected eof',911)
      if (.not. inp_a(field))
     $     call errquit('driver_input: failed to read field',911)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'cvgopt',field) .or. 
     $        inp_compare(.false.,'thresh',field)) then
         if(.not. inp_f(cvgopt))
     $        call errquit('driver_input: cvgopt/thresh',911)
         if (.not. rtdb_put(rtdb,'driver:cvgopt',mt_dbl,1,cvgopt))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'sadstp',field)) then
         if(.not. inp_f(sadstp)) 
     $        call errquit('driver_input: sadstp',911)
         if (.not. rtdb_put(rtdb,'driver:sadstp',mt_dbl,1,sadstp))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'nptopt',field) .or. 
     $        inp_compare(.false.,'maxiter',field)) then
         if(.not. inp_i(nptopt))
     $        call errquit('driver_input: nptopt/maxiter',911)
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,nptopt))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'linopt',field)) then
         if(.not. inp_i(linopt))
     $        call errquit('driver_input: linopt',911)
         if (.not. rtdb_put(rtdb,'driver:linopt',mt_int,1,linopt))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if (inp_compare(.false.,'clear',field)) then
         call driver_opt_del_hess()
c
      else if(inp_compare(.false.,'inhess',field)) then
         if(.not. inp_i(inhess)) 
     $        call errquit('driver_input: inhess',911)
         if (.not. rtdb_put(rtdb,'driver:inhess',mt_int,1,inhess))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'modupd',field)) then
         if(.not. inp_i(modupd))
     $        call errquit('driver_input: modupd',911)
         if (.not. rtdb_put(rtdb,'driver:modupd',mt_int,1,modupd))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'modsad',field)) then
         if(.not. inp_i(modsad))
     $        call errquit('driver_input: modsad',911)
         if (.not. rtdb_put(rtdb,'driver:modsad',mt_int,1,modsad))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'moddir',field)) then
         if(.not. inp_i(moddir)) 
     $        call errquit('driver_input: moddir',911)
         if (.not. rtdb_put(rtdb,'driver:moddir',mt_int,1,moddir))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'print',field)) then
         call util_print_input(rtdb,'driver')
c
      else if(inp_compare(.false.,'noprint',field)) then
         call util_print_input(rtdb,'driver')
c
      else if (inp_compare(.false.,'qstep',field)) then
         oqstep = .true.
         if (.not. rtdb_put(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if (inp_compare(.false.,'noqstep',field)) then
         oqstep = .false.
         if (.not. rtdb_put(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $        call errquit('driver_input: rtdb put failed',0)
      endif
c
      goto 10                   ! While loop
c
      end
      subroutine hnd_opt_gmax(rtdb)
      implicit double precision (a-h,o-z)
      integer  rtdb
#include "global.fh"
#include "coptopt.fh"
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      logical    zcoord
      logical    frozen
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      common/hnd_optfrz/nzfrz,izfrz(mxcoor),iatfrz(mxatom)
      common/hnd_optvar/zcoord,ncoord,mcoord
      common/hnd_zmtpar/nzmat,nzvar,nvar
      data zero /0.0d+00/
c
      gmax =zero
      gnorm=zero
      do i=1,ncoord 
         frozen=.false.
         do iz=1,nzfrz
            frozen=frozen.or.(i.eq.izfrz(iz))
         enddo
         if(.not.frozen) then
            if(abs(g(i)).gt.gmax) gmax=abs(g(i))
            gnorm=gnorm+g(i)*g(i)
         endif
      enddo
c
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
#include "util.fh"
#include "coptopt.fh"
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
c
      if((ga_nodeid().eq.0) .and.
     $     util_print('optimization',print_low)) then
         write(iw,9999) gmax,cvgopt,gmax0
      endif
c
      cvged=gmax.lt.cvgopt
c
      return
 9999 format('          largest component of the gradient =',F11.7,
     1       ' (cvgopt =',F11.7,') ',/,
     2       ' previous largest component of the gradient =',F11.7)
      end
      subroutine hnd_opt_energy(rtdb,geom)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "util.fh"
      integer  rtdb
      integer  geom
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
      if((ga_nodeid().eq.0.or.dbug) .and.
     $     util_print('optimization',print_low)) then
         write(iw,9999) nserch,npts
         if(.not.geom_print(geom)) 
     1     call errquit('hnd_opt_energy: print error',911)
      endif
c
      if(nzmat.gt.0) then
         call hnd_bandbi(geom,nat,c)
         if (util_print('optimization', print_low)) then
            call hnd_zmtprt(geom)
            if (ga_nodeid().eq.0) write(iw,9997)
         endif
      endif
c
      if (task_energy(rtdb)) then
         status=rtdb_get(rtdb,'task:energy',MT_DBL,1,e)
         if(dbug) then
            write(iw,9998) e
         endif
      else
         call errquit('optimize: energy failed', 0)
      endif
c
c     Disable printing to ecce of movecs after first point
c
      call movecs_ecce_print_off()
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
      subroutine hnd_opt_gradient(rtdb,geom)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      integer  geom
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
         call errquit('optimize: gradient failed', 0)
      endif
      call movecs_ecce_print_off() ! Disable MO printing after first time
c
      return
 9999 format(' in opt_gradient, eg = ')
 9998 format(3f12.6)
      end
      subroutine hnd_opt_tfhx(rtdb,geom)
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      integer  geom
      logical  gradients
      external gradients
      logical  geom_hnd_tfhx
      logical  dbug
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      integer ir, iw
      common/hnd_iofile/ir,iw
c
      integer nzmat, nzvar, nvar
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      logical zcoord
      integer ncoord, mcoord
      common/hnd_optvar/zcoord,ncoord,mcoord
c
      double precision c, zan
      integer nat
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c
      integer ibi, i_ibi, i, ncart
      double precision zero
      data zero /0.0d+00/
c
c
      dbug=.false.
c
c     ----- get -b- and -bi- -----
c
      call hnd_bandbi(geom,nat,c)
c
c     ----- transform -hx- to -hs- -----
c
      ncart=3*nat
      memsiz=ncart*nzvar
      memsiz=memsiz+mcoord*mcoord
      memsiz=memsiz+mcoord*mcoord
      if (.not. ma_push_get(MT_DBL,memsiz,'mem tfhx',i_ibi,ibi))
     $     call errquit('hnd_opt_tfhx: ma', memsiz )
      ihx=ibi+ncart*nzvar
      ihs=ihx+mcoord*mcoord
c
c     ----- only for node 0 -----
c     rjh ... NOW everyone does this
c
      do i=1,ncart*nzvar
         dbl_mb(i+ibi-1)=zero
      enddo
      call geom_hnd_get_data('binv', dbl_mb(ibi), ncart*nzvar)
      if(dbug) then
         write(iw,*) '-binv-'
         call hnd_prsq(dbl_mb(ibi),nzvar,ncart,ncart)
      endif
      do i=1,mcoord*mcoord
         dbl_mb(i+ihx-1)=zero
      enddo
      call geom_hnd_get_data('drv.hess',dbl_mb(ihx),mcoord*mcoord)
c
c     ----- now transform -----
c
      do i=1,mcoord*mcoord
         dbl_mb(i+ihs-1)=zero
      enddo
      if (.not. geom_hnd_tfhx(dbl_mb(ihs),dbl_mb(ihx),dbl_mb(ibi),
     $     nvar,ncart,ncart,mcoord))
     $     call errquit('hnd_opt_tfhx: hnd_tfhx failed', 911)
      if(dbug) then
         write(iw,*) '-hs-'
         call hnd_prsq(dbl_mb(ihs),nvar,nvar,mcoord)
      endif
c
c     ----- write out transformed hessian -----
c
      call geom_hnd_put_data('drv.hess',dbl_mb(ihs),mcoord*mcoord)
c
      if (.not. ma_pop_stack(i_ibi))
     $     call errquit('hnd_opt_tfhx: ma pop', memsiz )
c
      return
      end
      subroutine hnd_opt_tfgx(rtdb,geom)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      integer  geom
      logical  gradients
      external gradients
      logical  geom_hnd_tfgx
      logical  dbug
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      integer ir, iw
      common/hnd_iofile/ir,iw
c
      integer nzmat, nzvar, nvar
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      integer nzfrz, izfrz, iatfrz
      common/hnd_optfrz/nzfrz,izfrz(mxcoor),iatfrz(mxatom)
c
      double precision c, zan
      integer nat
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c
      double precision x0, x, dx, g0, g, ds, func, func0, gmax, gmax0,
     $     curv, alpha, gnorm
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
c
      logical geom_zmt_get_varinfo
      character*8 zvarname(mxcoor),test
      double precision zvarsign(mxcoor)
      double precision gs(mxcoor)
      integer ibi, i_ibi, i, ncart, j, num
      double precision sum
      double precision zero
      data zero /0.0d+00/
c
c
      dbug=.false.
c
c     ----- get -b- and -bi- -----
c
      call hnd_bandbi(geom,nat,x)
c
c     ----- transform -gx- to -gs- -----
c
      ncart=3*nat
      if (.not. ma_push_get(MT_DBL,ncart*nzvar,'mem bi',i_ibi,ibi))
     $     call errquit('hnd_opt_tfgx: ma', ncart*nzvar)
      if(dbug) then
         write(iw,*) 'in opt_tfgx'
         write(iw,*) '-gx-'
         call hnd_prsq(g,1,ncart,ncart)
      endif
      call geom_hnd_get_data('binv',dbl_mb(ibi),ncart*nzvar)
      if(dbug) then
         write(iw,*) '-binv-', ga_nodeid()
         call hnd_prsq(dbl_mb(ibi),nzvar,ncart,ncart)
         call util_flush(6)
      endif
      if (.not. geom_hnd_tfgx(gs,g,dbl_mb(ibi),nvar,ncart,ncart))
     $     call errquit('hnd_opt_tfgx: tfgx', ncart*nzvar)
      if(dbug) then
         write(iw,*) '-gs-',ga_nodeid()
         call hnd_prsq(gs,1,nvar,ncart)
         call util_flush(6)
      endif
      if (.not. ma_pop_stack(i_ibi))
     $     call errquit('hnd_opt_tfgx: ma pop', ncart*nzvar)
      do i=1,nvar
         g(i)=gs(i)
      enddo
      if(ncart.gt.nvar) then
         do i=nvar+1,ncart
            g(i)=zero
         enddo
      endif
      if(dbug) then
         write(iw,*) '-gs-',ga_nodeid(),nvar,ncart
         call hnd_prsq(g,1,ncart,ncart)
         call util_flush(6)
      endif
c
c     Apply constants and constraints to the gradient before
c     it is used for upating the Hessian.
c
      if (nzfrz .gt. 0) then
         do i = 1, nzfrz
            if (dbug) write(6,*) ' Freezing ', izfrz(i)
            g(izfrz(i)) = 0d0
         enddo
      endif
c
      if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $     call errquit('hnd_opt_tfgx: varinfo?',0)
      do i = 1, nzvar
         if (zvarname(i).ne.' ') then
            sum = 0d0
            num = 0
            test = zvarname(i)
            do j = 1, nzvar
               if (zvarname(j).eq.test) then 
                  sum = sum + g(j)*zvarsign(j)
                  num = num + 1
               endif
            enddo
            if (num .gt. 1) then
               sum = sum / dble(num)
               do j = 1, nzvar
                  if (zvarname(j).eq.test) then 
                     if(dbug)write(6,77) j, test, zvarsign(j), g(j),sum
 77                  format(' Constraining : ',i3,2x,a,2x,f3.1,2f10.6)
                     g(j) = sum*zvarsign(j)
                     zvarname(j) = ' '
                  endif
               enddo
            endif
         endif
      enddo
c
      if(dbug) then
         write(iw,*) '-gs- at end ',ga_nodeid(),nvar,ncart
         call hnd_prsq(g,1,ncart,ncart)
         call util_flush(6)
      endif
c
      return
      end
      subroutine hnd_opt_tfds(rtdb,geom)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      integer  geom
      logical  gradients
      external gradients
      logical  geom_hnd_tfds
      logical  dbug
      logical  iterate
      integer mxiter, mxatom, mxcart, mxzmat, mxcoor
      parameter (mxiter=5 )
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c
      integer ir, iw
      common/hnd_iofile/ir,iw
c
      integer nzmat, nzvar, nvar
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      double precision c, zan
      integer nat
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c
      double precision x0, x, dx, g0, g, ds, func, func0, gmax, gmax0,
     $     curv, alpha, gnorm
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      double precision d(mxcoor)
      integer ncart, i, i_ibi, ibi
      double precision zero, tol
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
         call hnd_bandbi(geom,nat,x)
c
         if (.not. ma_push_get(MT_DBL,ncart*nzvar,'mem bi',i_ibi,ibi))
     $        call errquit('hnd_opt_tfds: ma', ncart*nzvar)
         call geom_hnd_get_data('binv',dbl_mb(ibi),ncart*nzvar)
         if (.not. geom_hnd_tfds(d,dx,dbl_mb(ibi),nvar,ncart,ncart))
     $        call errquit('hnd_opt_tfds: tfds', ncart*nzvar)
         if (.not. ma_pop_stack(i_ibi))
     $        call errquit('hnd_opt_tfds: ma pop', ncart*nzvar)
c
         if(dbug) then
            write(iw,*) ' in -tfds- , ds = '
            write(iw,9999) (ds(i),i=1,nvar)
            write(iw,*) ' in -tfds- , bi = '
            call hnd_prsq(dbl_mb(ibi),nzvar,ncart,ncart)
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
      subroutine hnd_opt_tfdx(rtdb,geom)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer  rtdb
      integer  geom 
      logical  gradients
      external gradients
      logical  geom_hnd_tfdx
      logical  dbug
c
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c
      integer ir, iw
      common/hnd_iofile/ir,iw
c
      integer nzmat, nzvar, nvar
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      double precision c, zan
      integer nat
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c
      double precision x0, x, dx, g0, g, ds, func, func0, gmax, gmax0,
     $     curv, alpha, gnorm
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      double precision d(mxcoor)
c
      integer i_ib, ib, ncart, i
c
      dbug=.false.
c
c     ----- get -b- and -bi- -----
c
      call hnd_bandbi(geom,nat,x)
c
c     ----- transform -dx- to -ds- -----              
c
      ncart=3*nat
      do i=1,nvar
         d(i)=dx(i)
      enddo
      if (.not. ma_push_get(MT_DBL,ncart*nzvar,'mem b ',i_ib,ib))
     $     call errquit('hnd_opt_tfdx: ma ', ncart*nzvar)
      call geom_hnd_get_data('bmat',dbl_mb(ib),ncart*nzvar)
      if (.not. geom_hnd_tfdx(d,dx,dbl_mb(ib),nvar,ncart,ncart))
     $     call errquit('hnd_opt_tfdx: tfdx ', ncart*nzvar)
      if (.not. ma_pop_stack(i_ib))
     $     call errquit('hnd_opt_tfdx: ma pop ', ncart*nzvar)
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
      subroutine hnd_opt_hss_init(rtdb,geom)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
c
      integer     rtdb
      integer     geom
c
c parameters read from rtdb
c prefix:energy (later) rohf and uhf
c
c      character*255 name         ! name buffer for rtdb context never used
c      character*255 title        ! title buffer never used
c      double precision energy    ! energy (scf, rohf, uhf)
c      logical converged          ! energy convergence flag
c
      integer h_first, k_first   ! MA handle/index of first mem location
c
      logical geom_zmt_get_izmat
      logical geom_zmt_get_nizmat
      external geom_zmt_get_izmat
      external geom_zmt_get_nizmat
      logical dbug
      logical out 
c
c
      double precision zero
      parameter (zero=0.0d+00)
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)     
      parameter (mxcoor=1500)
c
      integer ir, iw
      common/hnd_iofile/ir,iw
c
      logical zcoord
      integer ncoord, mcoord
      common/hnd_optvar/zcoord,ncoord,mcoord
c
      integer nzmat, nzvar, nvar
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      double precision c, zan
      integer nat
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c
      character*16 atom_tags(mxatom)
      integer ndimx, ndimq, ndim, need
      integer k_hss, k_hssx, k_hssq, k_zsym, k_zmat, k_izmat, k_last
      integer iat, ixyz
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(dbug)  then
         write(iw,9999)
      endif
c
      if (.not. geom_ncent(geom,nat))
     1       call errquit('hnd_opt: geom_ncent?',911)
      if (dbug)  then
          write(iw,*)'   in hnd_opt_hss_init, nat =',nat
      endif
c
c     ----- define memory pointers -----
c     if internal coordinates ...
c           nzvar  = # of (redundant) internal coordinates
c           ncoord = # of independent internal coordinates
c           mcoord = max( 3*nat , nzvar )
c     if cartesian coordinates ...
c           ncoord = 3*nat
c           mcoord = 3*nat
c
      ndimx=3*nat
      ndimq=nzvar
      ndim =mcoord           
c
      need=     ndim *ndim 
      need=need+ndimx*ndimx
      need=need+ndimq*ndimq
      need=need+ndim *ndim 
      need=need+ndimq
      need=need+abs(nzmat)
      if (.not.
     1   MA_Push_Get(mt_dbl,need       ,'mem first',h_first ,k_first ))
     2   call errquit('hnd_opt_hss_init: malloc of first  failed?',911)
      k_hss  =k_first
      k_hssx =k_hss  +ndim *ndim 
      k_hssq =k_hssx +ndimx*ndimx
      k_zsym =k_hssq +ndimq*ndimq
      k_zmat =k_zsym +ndim *ndim 
      k_izmat=k_zmat +ndimq
      k_last =k_izmat+abs(nzmat)
c
      if(dbug) then
         write(iw,9997) k_first,k_last
         write(iw,*) ' k_first   = ',k_first
         write(iw,*) ' k_hss     = ',k_hss  
         write(iw,*) ' k_hssx    = ',k_hssx 
         write(iw,*) ' k_hssq    = ',k_hssq 
         write(iw,*) ' k_zsym    = ',k_zsym 
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
      if(.not.
     1   geom_cart_get(geom,nat,atom_tags,c,zan))
     2   call errquit('hnd_opt_hss_init: geom_cart_get failed?',911)
      if(dbug) then
         write(iw,*) 'in hnd_opt_hss_init , z and c = '
         write(iw,9996) (zan(iat),(c(ixyz,iat),ixyz=1,3),iat=1,nat)
      endif
c
c     ----- create -zmat- data if requested -----
c
      if(nzmat.gt.0) then
         if(dbug) then
            write(iw,9997) k_first,k_last
            write(iw,*) ' k_first   = ',k_first
            write(iw,*) ' k_hss     = ',k_hss
            write(iw,*) ' k_hssx    = ',k_hssx
            write(iw,*) ' k_hssq    = ',k_hssq
            write(iw,*) ' k_zsym    = ',k_zsym 
            write(iw,*) ' k_zmat    = ',k_zmat
            write(iw,*) ' k_izmat   = ',k_izmat
            write(iw,*) ' k_last    = ',k_last
         endif
         call hnd_bandbi(geom,nat,c)
         call geom_hnd_get_data('zsym',dbl_mb(k_zsym),ndim*ndim)
         call geom_hnd_get_data('zmat',dbl_mb(k_zmat),nzvar)
         if(.not.geom_zmt_get_nizmat(geom,nzmat))                 
     1      call errquit('hnd_opt_hss_init:zmt_get_nizmat failed?',911)
         if(.not.geom_zmt_get_izmat(geom, dbl_mb(k_izmat),nzmat))                 
     1      call errquit('hnd_opt_hss_init:zmt_get_izmat failed?',911)
      endif
c
c     ----- initial guess hessian ( cartesian or internal ) -----
c
      if(dbug) then
         write(iw,9997) k_first,k_last
         write(iw,*) ' k_first   = ',k_first
         write(iw,*) ' k_hss     = ',k_hss
         write(iw,*) ' k_hssx    = ',k_hssx
         write(iw,*) ' k_hssq    = ',k_hssq
         write(iw,*) ' k_zsym    = ',k_zsym 
         write(iw,*) ' k_zmat    = ',k_zmat
         write(iw,*) ' k_izmat   = ',k_izmat
         write(iw,*) ' k_last    = ',k_last
      endif
      call hnd_opt_hss_init_x(dbl_mb(k_hss  ),ndim ,
     1                        dbl_mb(k_hssx ),ndimx,
     2                        dbl_mb(k_hssq ),ndimq,          
     3                        dbl_mb(k_zsym ),
     4                        dbl_mb(k_zmat ),
     5                        dbl_mb(k_izmat))
      if(dbug) then
         write(iw,9997) k_first,k_last
         write(iw,*) ' k_first   = ',k_first
         write(iw,*) ' k_hss     = ',k_hss
         write(iw,*) ' k_hssx    = ',k_hssx
         write(iw,*) ' k_hssq    = ',k_hssq
         write(iw,*) ' k_zsym    = ',k_zsym 
         write(iw,*) ' k_zmat    = ',k_zmat
         write(iw,*) ' k_izmat   = ',k_izmat
         write(iw,*) ' k_last    = ',k_last
      endif
c
c     ----- release memory -----
c
      if (.not.MA_Pop_Stack(h_first))
     1   call errquit('hnd_opt_hss_init: ma_pop_stack k_first   ?',911)
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
      subroutine hnd_opt_hss_read(hess,n)
      implicit none
#include "global.fh"
c:: reads file in vib_vib format using vib_vib filename default
c:: Note the default filename is set in task_freq
c:: filenames must be made identical.
c
c format of file is ascii lower triangular elements only.
c
c::passed
      integer n ! [input] the rank of the hessian (3*number of atoms)
      double precision hess (n,n) ! [input] the matrix
c::local
      integer h_unit
      parameter (h_unit=47)
      character*255 fname
      integer i,j
      if (ga_nodeid().eq.0) then
c     
c::   -- open default file
         call util_file_name('hess',.false.,.false.,fname)
         open(unit=h_unit,file=fname,form='formatted',status='unknown',
     $        err=99990,access='sequential')
c
         rewind h_unit
c::   -- read  information
         do i = 1,n
            do j = 1,i
               read(h_unit,10000,err=99992,end=99992) hess(i,j)
            enddo
         enddo
         close(unit=h_unit,status='keep')
      endif
      call ga_sync()
      return
10000 format(f30.15)
99990 write(6,*)' could not open <',fname,'> as unknown file'
      call errquit('hnd_opt_hss_read: fatal error', 911)
99991 write(6,*)' could not open <',fname,'> as new file'
      call errquit('hnd_opt_hss_read: fatal error', 911)
99992 write(6,*)' error in reading <',fname,'> as hessian file'
      call errquit('hnd_opt_hss_read: fatal error', 911)
      end
      SUBROUTINE HND_OPT_HSS_INIT_X(
     1           HSS,NCOOR,HSSX,NCART,HSSQ,NINT,ZSYM,
     2           ZMAT,IZMAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "coptopt.fh"
C
C      -----        INITIALIZE THE HESSIAN MATRIX           -----
C      ----- IF -INHESS- EQUALS 1, INPUT THE HESSIAN MATRIX -----
C                 OTHERWISE MAKE UP INITIAL GUESS
C              IN CARTESIAN OR INTERNAL COORDINATES
C
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL       STATUS
      LOGICAL       DBUG
      LOGICAL       OUT
      LOGICAL       SOME
      LOGICAL       OLD_HESSIAN
      CHARACTER*8   ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION HSS (NCOOR,*)
      DIMENSION HSSX(NCART,*)
      DIMENSION HSSQ(NINT ,*)       
      DIMENSION ZSYM(NCOOR,*)       
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-OPTHSS-'/
      DATA ZERO   /0.0D+00/
      DATA HDIAGX /0.5D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT    
      IF(OUT) THEN
         WRITE(IW,9999) NAT,NCART,NZVAR,NVAR,NINT,NCOOR,NZMAT
      ENDIF
C
C     ----- CHECK IF THERE IS ALREADY A HESSIAN ON DISK -----
C
      OLD_HESSIAN=.FALSE.
      IF(INHESS.EQ.0) THEN
         call driver_opt_check_hess(ncoor, old_hessian)
         if (old_hessian)
     $        call geom_hnd_get_data('drv.hess',hss,ncoor*ncoor)
         IF(GA_NODEID().EQ.0) THEN
            IF(OLD_HESSIAN) THEN
               if (util_print('information', print_low))
     $              WRITE(IW,9993)
               IF(SOME) THEN
                  CALL HND_PRSQ(HSS,NCOOR,NCOOR,NCOOR)
               ENDIF
            ELSE
               if (util_print('information', print_low))
     $              WRITE(IW,9992)
            ENDIF
         ENDIF
      ENDIF
      CALL GA_SYNC()
C
C     ----- GENERATE A GUESS ... BUT DO NOT NECESSARILY USE IT -----
C     ----- IF drv.hess EXISTS, AND -INHESS- = 0               -----
C     ----- THEN THAT WILL BE THE ONE USED                     -----
C
      DO J=1,NCOOR
         DO I=1,NCOOR
            HSS(I,J)=ZERO
         ENDDO
      ENDDO
C
      IF(NZMAT.NE.0) THEN
         IF(OUT) THEN
            WRITE(IW,9999) NAT,NCART,NZVAR,NVAR,NINT,NCOOR,NZMAT
         ENDIF
         IF(INHESS.EQ.1) THEN
            if (ga_nodeid().eq.0 .and. 
     $           util_print('information', print_low)) write(6,*) 
     $           ' Reading Hessian from previous frequency calculation'
C
C     ----- READ CARTESIAN HESSIAN FROM FILE AND TRANSFORM -----
C
            DO J=1,NCART
               DO I=1,NCART
                  HSSX(I,J)=ZERO
               ENDDO
            ENDDO
            CALL HND_OPT_HSS_READ(HSSX,NCART)
            DO I=1,NCART
               DO J=1,I
                  HSS(I,J)=HSSX(I,J)
                  HSS(J,I)=HSSX(I,J)
               ENDDO
            ENDDO
            IF(SOME) THEN
               WRITE(IW,9997)
               CALL HND_PRSQ(HSS,NCART,NCART,NCOOR)
            ENDIF
C
C     ----- TRANSFORMATION TO INTERNAL COORDINATES IS DONE LATER -----
C
         ELSE
C
C     ----- INTERNALLY GENERATED (DIAGONAL) INITIAL HESSIAN -----
C
            CALL HND_OPT_HSS_INIT_ZMT(HSSQ,NINT,ZMAT,IZMAT,OUT)
C
C     ----- TRANSFORM TO SYMMETRY INTERNAL COORDINATES -----
C
            IF(OUT) THEN
               WRITE(IW,9996)
               CALL HND_PRSQ(ZSYM,NZVAR,NZVAR,NCOOR)
            ENDIF
            IF(SOME) THEN
               WRITE(IW,9995)
               CALL HND_PRSQ(HSSQ,NZVAR,NZVAR,NINT)
            ENDIF
            DO J=1,NZVAR
               DO I=1,NZVAR 
                  DUM=ZERO
                  DO K=1,NZVAR
                     DUM=DUM+ZSYM(K,I)*HSSQ(K,K)*ZSYM(K,J)
                  ENDDO
                  HSS(I,J)=DUM
               ENDDO
            ENDDO
            IF(SOME) THEN
               IF(OUT) THEN
                  WRITE(IW,9994)
               ENDIF
               WRITE(IW,9998)
               CALL HND_PRSQ(HSS,NZVAR,NZVAR,NCOOR)
            ENDIF
         ENDIF
      ELSE
         IF(INHESS.EQ.0) THEN
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
         ELSE
C
C     ----- READ CARTESIAN HESSIAN FROM FILE -----
C
            DO J=1,NCART
               DO I=1,NCART
                  HSSX(I,J)=ZERO
               ENDDO
            ENDDO
            CALL HND_OPT_HSS_READ(HSSX,NCART)
            DO I=1,NCART
               DO J=1,I
                  HSS(I,J)=HSSX(I,J)
                  HSS(J,I)=HSSX(I,J)
               ENDDO
            ENDDO
            IF(SOME) THEN
               WRITE(IW,9997)
               CALL HND_PRSQ(HSS,NCART,NCART,NCOOR)
            ENDIF
         ENDIF
      ENDIF
C
C     ----- SAVE INITIAL HESSIAN IF THERE IS NOT ONE ALREADY -----
C           IF THERE IS ONE, KEEP IT, DO NOT OVERWRITE IT
C           NOTE :
C     NCOOR = mcoord = max(3*nat,ncoord) with ncoord=nzvar
C
      if (.not. old_hessian)
     $     call geom_hnd_put_data('drv.hess',hss,ncoor*ncoor)
C
      RETURN
 9999 FORMAT(/,' IN  OPT_HSS-INIT_X , ', 
     1         ' NAT,NCART,NZVAR,NVAR,NINT,NCOOR,NZMAT = ',7I4)
 9998 FORMAT(/,10X,38(1H-),/,10X,'INITIAL -INTERN.COORD.- ',
     1                           'HESSIAN MATRIX',/,10X,38(1H-))
 9997 FORMAT(/,10X,38(1H-),/,10X,'INITIAL - CART.COORD. - ',
     1                           'HESSIAN MATRIX',/,10X,38(1H-))
 9996 FORMAT(/,' -ZSYM- ')
 9995 FORMAT(/,' -HSSQ- BEFORE SYMMETRIZATION ')
 9994 FORMAT(/,' -HSSQ- AFTER  SYMMETRIZATION ')
 9993 FORMAT(/,' hessian found on  -drv.hess- will be used')
 9992 FORMAT(/,' hessian internally generated or read from',
     1         ' external file will be used.')
      END
      SUBROUTINE HND_OPT_HSS_INIT_ZMT(HSSQ,NINT,ZMAT,IZMAT,SOME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (ZERO=0.0D+00,      ONE=1.0D+00,
     1     QUARTR=0.25D+00,   PT5=0.5D+00,
     2     THIRD=ONE/3.0D+00, TOHB2=6.423D-02,
     3     TOANGS=0.52917724924D+00)
      PARAMETER (MXATOM=500)
      LOGICAL OUT 
      LOGICAL SOME
      LOGICAL DBUG
      LOGICAL BADZ
      CHARACTER*8 HGUESS
      CHARACTER*8 GUESS
      CHARACTER*8 GUESS1
      CHARACTER*8 GUESS2
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
C     A GENERALIZED BADGER'S RULE, AND THESE PARAMETERS ARE FROM
C     D.R.HERSCHBACH, V.W.LAURIE, J.CHEM.PHYS. 35, 458-463(1961).
C     
      DATA IROW  /2*1,8*2,8*3,18*4,18*5,32*6/
      DATA APARM /1.26D+00,1.66D+00,1.84D+00,1.98D+00,2.03D+00,2.03D+00,
     1     1.66D+00,1.91D+00,2.28D+00,2.35D+00,2.33D+00,2.50D+00,
     2     1.84D+00,2.28D+00,2.41D+00,2.52D+00,2.61D+00,2.60D+00,
     3     1.98D+00,2.35D+00,2.52D+00,2.58D+00,2.66D+00,2.75D+00,
     4     2.03D+00,2.33D+00,2.61D+00,2.66D+00,2.85D+00,2.76D+00,
     5     2.03D+00,2.50D+00,2.60D+00,2.75D+00,2.75D+00,3.00D+00/
      DATA DPARM /2.50D-02,0.30D+00,0.38D+00,0.49D+00,0.51D+00,0.81D+00,
     1     0.30D+00,0.68D+00,0.74D+00,0.85D+00,0.68D+00,0.97D+00,
     2     0.38D+00,0.74D+00,1.18D+00,1.02D+00,1.28D+00,0.84D+00,
     3     0.49D+00,0.85D+00,1.02D+00,1.41D+00,0.86D+00,1.14D+00,
     4     0.51D+00,0.68D+00,1.28D+00,0.86D+00,1.62D+00,1.25D+00,
     5     0.81D+00,0.97D+00,0.84D+00,1.14D+00,1.25D+00,1.83D+00/
C     
C     COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
C     EXCEPT VAN DER WAALS RADII FOR HE,NE,AR,KR (SAME SOURCE),
C     AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT,RN
C     
      DATA (RCOV(NUCZ),NUCZ=1,2)/0.30D+00,1.22D+00/
      DATA (RCOV(NUCZ),NUCZ=3,10)
     1     /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     2     0.70D+00,0.66D+00,0.58D+00,1.60D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     1     /1.66D+00,1.36D+00,1.25D+00,1.17D+00,
     2     1.10D+00,1.04D+00,0.99D+00,1.91D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     1     /2.03D+00,1.74D+00,
     2     1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     3     1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     4     1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.98D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     1     /2.22D+00,1.92D+00,
     2     1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     3     1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     4     1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,2.09D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     1     /2.35D+00,1.98D+00,
     2     1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D0,
     3     1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D0,
     4     1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     5     1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     6     1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,2.20D+00/
C     
      DATA GUESS,GUESS1,GUESS2 /'GUESS   ','GUESS1  ','GUESS2  '/
C     
      DIST(IAT,JAT)=SQRT((C(1,IAT)-C(1,JAT))**2+
     1     (C(2,IAT)-C(2,JAT))**2+
     1     (C(3,IAT)-C(3,JAT))**2) * TOANGS
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
C     IN INTERNAL COORDINATE SPACE.
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
C     BOND STRETCH, B-A
C     
         if (itype .eq. 1) then
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
C     
C     ANGLE BEND, B-A-C
C     
         else if (itype.eq.2 .or. itype.eq.7) then
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
     1           * EXP(-0.44*(RAB + RAC - CAB - CAC))
C     
            IF(HGUESS.EQ.GUESS ) HSSQ(I,I) = HESS2
            IF(HGUESS.EQ.GUESS1) HSSQ(I,I) = HESS1
            IF(HGUESS.EQ.GUESS2) HSSQ(I,I) = HESS2
C     
C     TORSION, D-A-B-C
C     
         else if (itype .eq. 3) then
            IZ1 = IZMAT(IZ+1)
            IZ2 = IZMAT(IZ+2)
            IZ3 = IZMAT(IZ+3)
            IZ4 = IZMAT(IZ+4)
            IZD = NUC(IZ1)
            IZA = NUC(IZ2)
            IZB = NUC(IZ3)
            IZC = NUC(IZ4)
            IF(BADZ(IZA).OR.BADZ(IZB).OR.BADZ(IZC).OR.BADZ(IZD))THEN
               HSSQ(I,I) = QUARTR
               GO TO 180
            ENDIF
C     
            HESS1 = QUARTR
C     
C     COUNT BONDS TO ATOMS A AND B, EXCEPT FOR A-B ITSELF
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
     1              J.NE.IZ2  .AND.  J.NE.IZ3) NBA = NBA+1
               IF(DIST(IZ3,J)  .LT.  1.2D+00*CBJ  .AND.
     1              J.NE.IZ3  .AND.  J.NE.IZ2) NBB = NBB+1
 135        CONTINUE
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
C     
C     OUT OF PLANE BEND, X-A-B-C, X IS HOOKED TO A
C     
         else if (itype .eq. 4) then
            IZ1 = IZMAT(IZ+1)
            IZ2 = IZMAT(IZ+2)
            IZ3 = IZMAT(IZ+3)
            IZ4 = IZMAT(IZ+4)
            IZX = NUC(IZ1)
            IZA = NUC(IZ2)
            IZB = NUC(IZ3)
            IZC = NUC(IZ4)
            IF(BADZ(IZA) .OR. BADZ(IZB) .OR. BADZ(IZC) .OR. 
     $           BADZ(IZX)) THEN
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
C     
C     ALL OTHER COORDINATE TYPES...
C     
         else 
            HSSQ(I,I)=QUARTR
         endif
C     
 180     CONTINUE
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
         IF(ITYPE.EQ.7) IZ=IZ+5 ! RJH WAS 7 ... only place 7 appeared????
         IF(OUT) THEN
            WRITE(IW,9148) I,ITYPE,IZ1,IZ2,IZ3,IZ4,
     1           IZA,IZB,IZC,IZD,HSSQ(I,I)
         ENDIF
 190  CONTINUE
C     
c--   
c     DO J=1,NZVAR
c     DO I=1,NZVAR
c     HSSQ(I,J)=ZERO
c     ENDDO
c     HSSQ(J,J)=0.5D+00
c     ENDDO
c--   
      IF(OUT) THEN
         WRITE(IW,9048)
         CALL HND_PRSQ(HSSQ,NZVAR,NZVAR,NINT)
      ENDIF
C     
      RETURN
 9208 FORMAT(' IZMAT = ',12I4)
 9188 FORMAT(' IN  OPT_HSS-INIT_ZMT ,', 
     1     ' NAT,NVAR,NZVAR,NINT,NZMAT = ',5I4)
 9168 FORMAT(' IN  OPT_HSS_INIT_ZMT , GUESS INTERNAL HESSIAN ')
 9148 FORMAT(' I,ITYPE,HSSQ(I,I) = ',2I5,8I3,F12.8)
 9128 FORMAT('  ZMAT = ',5F7.3)
 9048 FORMAT(/,1X,38(1H-),
     1     /,1X,'INITIAL -INTERN.COORD.- HESSIAN MATRIX',
     2     /,1X,38(1H-))
 9008 FORMAT(' IN  OPT_HSS_INIT_ZMT , HGUESS = ',A8) 
      END
      subroutine driver_opt_check_hess(ncoor, old_hessian)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
c
      integer ncoor
      logical old_hessian
      character*255 filename
c
      integer m
c
c     Look at an existing hessian file and verify it
c
      call util_file_name('drv.hess',
     1     .false.,.false.,filename)
c
      if (ga_nodeid() .eq. 0) then
         open(32,file=filename,form='unformatted',status='old',err=10)
         write(6,*) ' CHECK: OPENED OK',filename
         read(32,err=11) m
         write(6,*) ' CHECK: READ OK ',m,ncoor*ncoor
         if (m.ne.ncoor*ncoor) goto 11
         close(32)
         old_hessian = .true.
         goto 20
c
 11      close(32)
 10      old_hessian = .false.
      endif
c
 20   call ga_brdcst(323, old_hessian, mitob(1), 0)
c
      end
      subroutine driver_opt_del_hess()
      implicit none
#include "util.fh"
c
c     Delete the Hessian information restart file
c     
      character*255 opt_hess_fil
c
      call util_file_name('drv.hess',
     1     .false.,.false.,opt_hess_fil)
c
      call util_file_unlink(opt_hess_fil)
c
      if (util_print('information',print_low)) then
         write(6,*)
         write(6,*) ' Deleted DRIVER Hessian restart file '
         write(6,*)
      endif
c
      end
c$$$      logical function driver_opt_put_hess(hess,nint,mint)   
c$$$      implicit none
c$$$#include "mafdecls.fh"
c$$$#include "tcgmsg.fh"
c$$$#include "msgtypesf.h"
c$$$#include "global.fh"
c$$$#include "rtdb.fh"
c$$$#include "util.fh"
c$$$c
c$$$      integer nint 
c$$$      integer mint
c$$$      integer i
c$$$      integer j
c$$$      integer ir
c$$$      integer iw
c$$$      integer ift
c$$$      logical dbug
c$$$      logical out
c$$$      common/hnd_iofile/ir,iw
c$$$      double precision hess(mint,mint)  
c$$$      character*255 opt_hess_fil
c$$$      data ift /31/
c$$$c
c$$$      dbug=.false.
c$$$      dbug=dbug.or.util_print('hess_fil', print_debug)
c$$$      out =.false.
c$$$      out =out.or.dbug
c$$$c
c$$$      call util_file_name('drv.hess',
c$$$     1     .false.,.false.,opt_hess_fil)
c$$$      open(unit=ift, file=opt_hess_fil, form='unformatted',
c$$$     1     access='sequential', status='unknown', err=911)
c$$$c
c$$$      if(out) then
c$$$         write(iw,*) ' _put_hess', opt_hess_fil
c$$$         write(iw,*) 'nint, mint = ',nint,mint 
c$$$      endif
c$$$      if(dbug) then
c$$$         call hnd_prsq(hess,nint,nint,mint)    
c$$$      endif
c$$$c
c$$$      rewind ift
c$$$      write(ift) mint
c$$$      write(ift) ((hess(i,j),i=1,mint),j=1,mint)
c$$$c
c$$$      close(ift,status='keep')
c$$$c
c$$$      driver_opt_put_hess=.true.
c$$$      return
c$$$c
c$$$ 911  call errquit('driver_opt_put_hess: open failed',0)
c$$$c
c$$$      end
c$$$      logical function driver_opt_get_hess(hess,nint,mint)   
c$$$      implicit none
c$$$#include "mafdecls.fh"
c$$$#include "tcgmsg.fh"
c$$$#include "msgtypesf.h"
c$$$#include "global.fh"
c$$$#include "rtdb.fh"
c$$$#include "util.fh"
c$$$c
c$$$      integer nint
c$$$      integer mint
c$$$      integer i
c$$$      integer j
c$$$      integer ir
c$$$      integer iw
c$$$      integer ift
c$$$      logical dbug
c$$$      logical out
c$$$c
c$$$      integer mmint
c$$$c
c$$$      common/hnd_iofile/ir,iw
c$$$      double precision hess(mint,mint) 
c$$$      character*255 opt_hess_fil
c$$$      data ift /31/
c$$$c
c$$$      dbug=.false.
c$$$      dbug=dbug.or.util_print('hess_fil', print_debug)
c$$$      out =.false.
c$$$      out =out.or.dbug
c$$$c
c$$$      call util_file_name('drv.hess',
c$$$     1     .false.,.false.,opt_hess_fil)
c$$$      open(unit=ift, file=opt_hess_fil, form='unformatted',
c$$$     1     access='sequential', status='unknown', err=911)
c$$$      if(out) then
c$$$         write(iw,*) ' _get_hess', opt_hess_fil
c$$$         write(iw,*) 'nint, mint = ',nint,mint 
c$$$      endif
c$$$c
c$$$      rewind ift
c$$$      read(ift,end=910,err=910) mmint
c$$$      if (mint .ne. mmint) goto 910
c$$$      read(ift,end=910,err=910) ((hess(i,j),i=1,mint),j=1,mint)
c$$$c
c$$$      close(ift,status='keep')
c$$$c
c$$$      if(dbug) then
c$$$         call hnd_prsq(hess,nint,nint,mint)    
c$$$      endif
c$$$c
c$$$      driver_opt_get_hess=.true.
c$$$      return
c$$$c
c$$$c     ----- unable to complete the read ---
c$$$c
c$$$ 910  close(ift,status='keep')
c$$$      if(out) then
c$$$         write(iw,*) 'unable to complete the read in _get_hess'
c$$$      endif
c$$$      driver_opt_get_hess=.false.
c$$$      return 
c$$$c
c$$$c     ----- unable to open ; stop ... -----
c$$$c
c$$$ 911  call errquit('driver_opt_get_hess: open failed',0)
c$$$      return
c$$$c
c$$$      end
      subroutine hnd_opt_search_dir(rtdb,geom)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "util.fh"
c
      integer    rtdb, geom
      logical    status
c
      integer mxatom, mxcart, mxzmat, mxcoor
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c
      logical dbug
      logical out
      logical baker
c
      integer ir, iw
      common/hnd_iofile/ir,iw
c
      double precision c, zan
      integer nat
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
c
      integer nzmat, nzvar, nvar
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      logical zcoord
      integer ncoord, mcoord
      common/hnd_optvar/zcoord,ncoord,mcoord
c
      integer ndimx, ndimq, ndimi, ndim, ndim1
      integer need, i, j, i_ifirst, ifirst
      integer i10, i20, i30, i40, i50, i60, i70, i80
      integer j10, j20, j30, j40
c
      double precision x(1)
      equivalence (x(1),dbl_mb(1))
      double precision zero
      data zero /0.0d+00/
c
      dbug =.false.
      dbug = dbug.or.util_print('search dir debug',print_debug)
      out  =.false.
      out  = out .or.util_print('search dir out', print_debug)
      out  = out .or.dbug
c
      baker=.false.
c
      if(out) then
         write(iw,*) 'start of hnd_opt_search_dir, baker = ',
     1                                             baker
      endif
c     if internal coordinates ...
c           nzvar  = # of (redundant) internal coordinates
c           ncoord = # of independent internal coordinates
c           mcoord = max( 3*nat , nzvar )
c     if cartesian coordinates ...
c           ncoord = 3*nat
c           mcoord = 3*nat
c
      ndimx=3*nat
      ndimq=nzvar
      ndimi=ncoord
      ndim =mcoord           
      if(baker) then
         ndim1=ndim +1
      else
         ndim1=ndim 
      endif
      if(out) then
         write(iw,*) 'in search_dir,',
     1   ' ndimx,ndimq,ndimi,ndim,ndim1,nzvar,nvar = ',
     2     ndimx,ndimq,ndimi,ndim,ndim1,nzvar,nvar
      endif
c
      need=      ndim1* ndim1
      need=need+ ndim1* ndim1
      need=need+max(6,(ndim1*(ndim1+1))/2)
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
      i40=i30+max(6,(ndim1*(ndim1+1))/2)
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
         call geom_hnd_get_data('drv.hess',x(j10),ndim*ndim)
         if(dbug) then
            write(iw,*) 'in -search_dir- , int.coord. hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
         endif
         call geom_hnd_get_data('zind',x(j20),ndim*ndim)
         if(dbug) then
            write(iw,*) 'in -search_dir- , zind = '                 
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
         if(dbug) then
            write(iw,*) 'in -search_dir- , hessian = '
            call hnd_prsq(x(i10),ndimq,ndimq,ndim1)
         endif
      else
c
         call geom_hnd_get_data('drv.hess',x(i10),ndim1*ndim1)
         if(dbug) then
            write(iw,*) 'in -search_dir- , cartesian hessian = '
            call hnd_prsq(x(i10),ndimx,ndimx,ndim1)
         endif
c--      call hnd_opt_eckart(x(i10),x(i20),x(i30),x(i40),ndim1)
      endif
c
c     --- search direction (quasi-newton or augmented hessian) ---
c
         if(baker) then
            call hnd_opt_bkr_search(ncoord,
     1      x(i10),x(i20),x(i30),x(i40),x(i50),x(i60),x(i70),
     2      ndim1,zcoord)
         else
            call hnd_opt_dir_search(ncoord,
     1      x(i10),x(i20),x(i30),x(i40),x(i50),x(i60),x(i70),
     2      ndim1,zcoord,geom)
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
         call geom_hnd_get_data('zind',x(j20),ndim*ndim)
c
         call hnd_tfhsi(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
c
         if(out) then
            write(iw,*) 'in -search_dir- , ',
     1                  'back-transformed hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
c--
c--      these are for dbug only. one may not continue after those
c--         if(dbug) then
c--            call hnd_tfhs(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
c--            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
c--            stop
c--         endif
         endif
         call geom_hnd_put_data('drv.hess',x(j10),ndim*ndim)       
      else
         call geom_hnd_put_data('drv.hess',x(i10),ndim1*ndim1)       
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
      subroutine hnd_opt_eckart(hess,proj,tr,t,ndim)
      implicit double precision (a-h,o-z)
      parameter (zero=0.0d+00)
      parameter ( one=1.0d+00)
      parameter ( tol=1.0d-10)
      parameter (mxatom=500)
      logical      dbug
      logical      out
      character*24 errmsg
      common/hnd_iofile/ir,iw
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      dimension hess(ndim,*)
      dimension proj(ndim,*)
      dimension   tr(ndim,*)
      dimension    t(*)
      dimension    g(3)
      data errmsg /'program stop in -eckart-'/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      ncart=3*nat
c
      if(out) then
         write(iw,*) 'in -hnd-opt-eckart- ....'
         write(iw,*) 'un-projected hessian'
         call hnd_prsq(hess,ncart,ncart,ndim)
      endif
c
      do jxyz=1,6
         do ixyz=1,ncart
            tr(ixyz,jxyz)=zero
         enddo
      enddo
c
c     ----- center of mass -----
c
      do ixyz=1,3
         dum=zero
         do iat=1,nat
            dum=dum+c(ixyz,iat)
         enddo
         g(ixyz)=dum/dble(nat)
      enddo
c
c     ----- define translations -----
c
      do jxyz=1,3
         do iat=1,nat
            ixyz=jxyz+3*(iat-1)
            tr(ixyz,jxyz)=one
         enddo
      enddo
c
c     ----- define rotations -----
c 
      do jxyz=4,6
         if(jxyz.eq.4) then
            do iat=1,nat
               ixyz=3*(iat-1)
               tr(2+ixyz,jxyz)=-c(3,iat)+g(3)
               tr(3+ixyz,jxyz)= c(2,iat)-g(2)
            enddo
         elseif(jxyz.eq.5) then
            do iat=1,nat
               ixyz=3*(iat-1)
               tr(1+ixyz,jxyz)= c(3,iat)-g(3)
               tr(3+ixyz,jxyz)=-c(1,iat)+g(1)
            enddo
         elseif(jxyz.eq.6) then
            do iat=1,nat
               ixyz=3*(iat-1)
               tr(1+ixyz,jxyz)=-c(2,iat)+g(2)
               tr(2+ixyz,jxyz)= c(1,iat)-g(1)
            enddo
         endif
      enddo
      if(dbug) then
         write(iw,*) 'translations and rotations vectors'
         call hnd_prsq(tr,6,ncart,ndim)
      endif
c
c     ----- orthonormalize rot. + trans. -----
c
      nxyz=0
      do jxyz=1,6
         dum=dnrm2(ncart,tr(1,jxyz),1)
         if(abs(dum).gt.tol) then
            nxyz=nxyz+1
            dum=one/dum
            call dscal(ncart,dum,tr(1,jxyz),1)
            if(jxyz.lt.6) then
               do kxyz=jxyz+1,6
                  dum= -ddot(ncart,    tr(1,jxyz),1,tr(1,kxyz),1)
                  call daxpy(ncart,dum,tr(1,jxyz),1,tr(1,kxyz),1)
               enddo
            endif
         else
            do ixyz=1,ncart
               tr(ixyz,jxyz)=zero
            enddo
         endif
      enddo
      if(out) then
         write(iw,*) 'translations and rotations vectors'
         call hnd_prsq(tr,6,ncart,ndim)
      endif
      if(out.and.nxyz.ne.6) then
         write(iw,*) '# of transl.+rot. found',nxyz
      endif
c
c     ----- eckart transform ----
c
      do jxyz=1,ncart
         do ixyz=1,ncart
            proj(ixyz,ncart+1-jxyz)=zero
         enddo
         if(jxyz.le.nxyz) then
            do ixyz=1,ncart
               proj(ixyz,ncart+1-jxyz)=tr(ixyz,jxyz)
            enddo
         endif
      enddo
      mxyz=nxyz+1
      do jxyz=1,ncart
         if(mxyz.le.ncart) then
            do ixyz=1,ncart
               proj(ixyz,ncart+1-mxyz)=zero
            enddo
            proj(jxyz,ncart+1-mxyz)=one
            do kxyz=1,mxyz-1 
               dum=zero
               do ixyz=1,ncart
                  dum=dum+proj(ixyz,ncart+1-kxyz)
     1                   *proj(ixyz,ncart+1-mxyz)
               enddo
               do ixyz=1,ncart
                  proj(ixyz,ncart+1-mxyz)=proj(ixyz,ncart+1-mxyz)-
     1                                    proj(ixyz,ncart+1-kxyz)*dum
               enddo   
            enddo
            dum=zero
            do ixyz=1,ncart
               dum=dum+proj(ixyz,ncart+1-mxyz)*proj(ixyz,ncart+1-mxyz)
            enddo
            if(dum.gt.tol) then
               dum=one/sqrt(dum)
               do ixyz=1,ncart
                  proj(ixyz,ncart+1-mxyz)=proj(ixyz,ncart+1-mxyz)*dum
               enddo
               mxyz=mxyz+1
            endif
         endif
      enddo
      if(dbug) then
         write(iw,*) 'eckart transformation'
         call hnd_prsq(proj,ncart,ncart,ndim)
      endif
      do jxyz=1,nxyz 
         do ixyz=1,ncart
            proj(ixyz,ncart+1-jxyz)=zero
         enddo
      enddo
      if(out) then
         write(iw,*) 'eckart transformation'
         call hnd_prsq(proj,ncart,ncart,ndim)
      endif
      if(mxyz.le.ncart) then
         write(iw,*) 'something wrong in -eckart- transform'
         write(iw,*) 'nxyz,mxyz,ncart = ',nxyz,mxyz,ncart
         call hnd_hnderr(3,errmsg)
      endif
c
c     ----- transform hessian -----
c
      if(dbug) then
         write(iw,*) 'un-transformed hessian'
         call hnd_prsq(hess,ncart,ncart,ndim)
      endif
c
      do ixyz=1,ncart
         do jxyz=1,ncart
            dum=zero
            do kxyz=1,ncart
               dum=dum+hess(ixyz,kxyz)*proj(kxyz,jxyz)
            enddo
            t(jxyz)=dum
         enddo
         do jxyz=1,ncart
            hess(ixyz,jxyz)=t(jxyz)
         enddo
      enddo
c
      if(dbug) then
         write(iw,*) 'half-transformed hessian'
         call hnd_prsq(hess,ncart,ncart,ndim)
      endif
c
      do jxyz=1,ncart
         do ixyz=1,ncart
            dum=zero
            do kxyz=1,ncart
               dum=dum+proj(kxyz,ixyz)*hess(kxyz,jxyz)
            enddo
            t(ixyz)=dum
         enddo
         do ixyz=1,ncart
            hess(ixyz,jxyz)=t(ixyz)
         enddo
      enddo
      if(out) then
         write(iw,*) '   transformed hessian'
         call hnd_prsq(hess,ncart,ncart,ndim)
      endif
c
c     ----- project gradient -----
c

c
      return
      end
      SUBROUTINE HND_OPT_BKR_SEARCH(NVAR,
     1           HESS,HESVEC,HESST,HESEIG,HESEDM,T,IA,NDIM,ZCOORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "coptopt.fh"
C
C     ----- SELECT SEARCH DIRECTION -----
C
      LOGICAL   OUT
      LOGICAL   DBUG
      LOGICAL   CVGED
      LOGICAL   ZCOORD
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      DIMENSION HESS(NDIM,*),HESVEC(NDIM,*),HESST(*),HESEIG(*)
      DIMENSION HESEDM(*),T(*),IA(*)
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TINY   /1.0D-06/
      DATA SMALL  /1.0D-04/
      DATA DAMP1  /5.0D+00/
      DATA DAMP2  /2.0D+00/
      DATA TENM1  /1.0D-01/
C
      DBUG=.FALSE.
      OUT =.FALSE.
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
     1           HESS,HESVEC,HESST,HESEIG,HESEDM,T,IA,NDIM,ZCOORD,
     $     geom)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "coptopt.fh"
C
C     ----- SELECT SEARCH DIRECTION -----
C
      integer geom
      LOGICAL   OUT
      LOGICAL   DBUG
      LOGICAL   CVGED
      LOGICAL   ZCOORD
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      character*8 zvarname(mxcoor),test
      double precision zvarsign(mxcoor)
      logical geom_zmt_get_varinfo
      DIMENSION HESS(NDIM,*),HESVEC(NDIM,*),HESST(*),HESEIG(*)
      DIMENSION HESEDM(*),T(*),IA(*)
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TINY   /1.0D-06/
      DATA DAMP1  /5.0D+00/
      DATA DAMP2  /2.0D+00/
      DATA HUNDRD /1.0D+02/
      DATA TENM1  /1.0D-01/
      DATA TENM8  /1.0D-08/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(DBUG) THEN
         WRITE(IW,9999) NVAR,NDIM
      ENDIF
C
      IF(ZCOORD) THEN
         SMALL=1.0D-04
      ELSE
         SMALL=1.0D-01
      ENDIF
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
C     ----- DAMP OUT FROZEN COORDINATES -----
C
      IF(NZFRZ.GT.0) THEN
         DO IZ=1,NZFRZ
            I=IZFRZ(IZ)
            DO J=1,NVAR
               HESS(I,J)=ZERO
               HESS(J,I)=ZERO
            ENDDO
            HESS(I,I)=  1d9  ! RJH ... was HUNDRD
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
      IF(DBUG) THEN
         CALL HND_PRTR(HESST,NVAR)
      ENDIF
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
         IF(ABS(HESEIG(I)).LT.TENM8) THEN
            HESEIG(I)=ZERO
            HESEDM(I)=ZERO
            DO J=1,NVAR
               HESVEC(J,I)=ZERO
            ENDDO
         ELSEIF(ABS(HESEIG(I)).LT.TINY) THEN
            HESEIG(I)=SMALL
            HESEDM(I)=SMALL
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
         DX(I)=-DUM
      ENDDO
C
c
c     Apply constants and constraints to the gradient before
c     it is used for upating the Hessian.
c

      IF(NZFRZ.GT.0) THEN
         DO I=1,NZFRZ
            IZ=IZFRZ(I)
            DX(IZ)=ZERO
         ENDDO
      ENDIF
c
      if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $     call errquit('hnd_opt_tfgx: varinfo?',0)
      if(dbug) then
         do i = 1, nvar
            write(6,*) ' zvarinfo ', zvarname(i),zvarsign(i)
         enddo
      endif
      do i = 1, nvar
         if (zvarname(i).ne.' ') then
            sum = 0d0
            num = 0
            test = zvarname(i)
            do j = 1, nvar
               if (zvarname(j).eq.test) then 
                  sum = sum + dx(j)*zvarsign(j)
                  num = num + 1
               endif
            enddo
            if (num .gt. 1) then
               sum = sum / dble(num)
               do j = 1, nvar
                  if (zvarname(j).eq.test) then 
                     if(dbug)write(6,77) j, test, zvarsign(j), sum
 77                  format(' XConstraining : ',i3,2x,a,2x,f3.1,f10.6)
                     dx(j) = sum*zvarsign(j)
                     zvarname(j) = ' '
                  endif
               enddo
            endif
         endif
      enddo



C
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
#include "rtdb.fh"
#include "global.fh"
#include "util.fh"
#include "coptopt.fh"
      INTEGER   RTDB
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL      DBUG
      LOGICAL      CVGED
      CHARACTER*16 ATMNAM
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLLAB/ATMNAM(MXATOM)
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      COMMON/HND_OPTFUN/E,EG(MXCART)
C
      DBUG=.FALSE.
C
C     ----- PRINT OPTIMIZATION SUMMARY -----
C
      DELF=FUNC-FUNC0
      IF((GA_NODEID().EQ.0.OR.DBUG) .and.
     $     util_print('options', print_low)) THEN
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
      logical function geom_hnd_tfhx(hs,hx,bi,ns,nx,mx,mh)
      implicit double precision (a-h,o-z)
      logical out
      logical dbug
      common/hnd_iofile/ir,iw
      dimension hs(mh,*),hx(mh,*),bi(mx,*)
      data zero   /0.0d+00/
      data tenm08 /1.0d-08/
c
      dbug=.false.
      out =.false.
      if(dbug) then
         write(iw,*) ' in -tfhx- bi = '
         call hnd_prsq(bi,ns,nx,nx)
         write(iw,*) ' in -tfhx- hs = '
         call hnd_prsq(hx,nx,nx,mh)
      endif
c
      do k=1,nx
         do j=1,ns
            dum=zero
            do i=1,nx
               dum=dum+hx(k,i)*bi(i,j)
            enddo
            hs(k,j)=dum
         enddo
      enddo
      do j=1,mh
         do i=1,mh
            hx(i,j)=hs(i,j)
         enddo
      enddo
      do k=1,ns
         do i=1,ns
            dum=zero
            do j=1,nx
               dum=dum+bi(j,i)*hx(j,k)
            enddo
            hs(i,k)=dum
         enddo
      enddo
c
      if(out) then
         write(iw,*) ' in -tfhx- hs = '
         call hnd_prsq(hs,ns,ns,mh)
      endif
c
      do j=1,ns
         do i=1,ns
            if(abs(hs(i,j)).lt.tenm08) hs(i,j)=zero
         enddo
      enddo
c
      geom_hnd_tfhx=.true.
      return
 9999 format(i5,f15.10)
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
      data tenm15 /1.0d-15/
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
         if(abs(dx(i)).lt.tenm15) dx(i)=zero
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
#include "rtdb.fh"
#include "coptopt.fh"
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
      LOGICAL  UTIL_TEST_TIME_REMAINING
      EXTERNAL UTIL_TEST_TIME_REMAINING
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
      COMMON/HND_MOLXYZ/C(MXCART),ZAN(MXATOM),NAT         
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED 
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      COMMON/HND_OPTFUN/E,EG(MXCART)
      COMMON/HND_OPTTIM/ENERGY_TIME,GRADIENT_TIME
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      DIMENSION ALP(5),FUNL(5)
      DATA ZERO,TWO,THREE       /0.0D+00,2.0D+00,3.0D+00/
      DATA TOLALP,TOLFUN        /1.0D-01,1.5D-01/
      DATA FALMIN,FALNOR,FALMAX /1.0D-01,2.0D+00,5.0D+00/
      DATA DALNOR,DALMAX        /0.9D+00,0.5D+00/
      DATA ONEPT5               /1.5D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
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
      dxnorm = 0d0
      DO I = 1,NCART 
         X(I) = X0(I)+ALPH*DX(I)
         dxnorm = dxnorm + dx(i)*dx(i)
         C(I) = X(I)
      ENDDO
      dxnorm = alph*sqrt(dxnorm)
      NPTS=NPTS+1
      IF(NPTS.GT.NPTOPT) THEN
         EXCESS=.TRUE.
         RETURN
      ENDIF
C
C     ----- WRITE TO -NWCHEM- -----
C
      STATUS=GEOM_CART_GET(GEOM,NAT_NW,TAGS_NW,COORDS_NW,
     1                                         CHARGE_NW)
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
      STATUS=GEOM_CART_SET(GEOM,NAT_NW,TAGS_NW,COORDS_NW,
     1                                         CHARGE_NW)

c
c     Enforce symmetry on the geometry both inside the
c     NWChem geometry object and the **&@$ing hondo duplicates.
c
      call sym_geom_project(geom, dxnorm)
      STATUS=GEOM_CART_GET(GEOM,NAT_NW,TAGS_NW,COORDS_NW,
     1                                         CHARGE_NW)

      DO I = 1,NCART 
         X(I) = coords_nw(i)
         C(I) = X(I)
      ENDDO

      STATUS=GEOM_RTDB_STORE(RTDB,GEOM,'geometry')
C
      IF(OUT) THEN
         WRITE(IW,9999)
         WRITE(IW,9998) (COORDS_NW(I),I=1,NCART)
      ENDIF
C
C     ----- CALL FUNCTION EVALUATION -----
C
         STATUS=UTIL_TEST_TIME_REMAINING(RTDB,INT(ENERGY_TIME*ONEPT5))
         IF(OUT) THEN
            WRITE(IW,*) 'ENOUGH TIME REMAINING FOR ENERGY = ? ',STATUS
         ENDIF
         IF(.NOT.STATUS) THEN
            if (ga_nodeid().eq.0) 
     $           WRITE(IW,*) 'Not enough time remaining,',
     $           ' shutting down ... '
            NPTOPT=NPTS-1
            EXCESS=.TRUE.
            RETURN
         ENDIF
C
      CALL HND_OPT_ENERGY(RTDB,GEOM)
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
         STATUS=UTIL_TEST_TIME_REMAINING(RTDB,INT(GRADIENT_TIME*ONEPT5))
         IF(OUT) THEN
            WRITE(IW,*) 'ENOUGH TIME REMAINING FOR GRADIENT = ? ',STATUS
         ENDIF
         IF(.NOT.STATUS) THEN
            if (ga_nodeid().eq.0) 
     $           WRITE(IW,*) 'Not enough time remaining,',
     $           ' shutting down ... '
            NPTOPT=NPTS-1
            EXCESS=.TRUE.
            RETURN
         ENDIF
C
      CALL HND_OPT_GRADIENT(RTDB,GEOM)
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
      ALPHA=ALPH
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
      SUBROUTINE HND_OPT_SEARCH_LIN_INT(RTDB,GEOM,RSTART,EXCESS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "coptopt.fh"
C
C     ----- ONE DIMENSIONAL SEARCH.  -----
c
c     RJH ... modified so search happens in internals.
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
C
      INTEGER RTDB
      INTEGER GEOM
      LOGICAL STATUS
      LOGICAL  UTIL_TEST_TIME_REMAINING
      EXTERNAL UTIL_TEST_TIME_REMAINING
c
      double precision s(mxcart),news(mxcart)
C
      LOGICAL RSTART
      LOGICAL EXCESS
      LOGICAL LOWER,ILLCON,RECALC,RETRY
      LOGICAL CVGED
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(MXCART),ZAN(MXATOM),NAT         
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_OPTRUN/NPTS,NSERCH,NUPDAT,CVGED 
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      COMMON/HND_OPTFUN/E,EG(MXCART)
      COMMON/HND_OPTTIM/ENERGY_TIME,GRADIENT_TIME
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      DIMENSION ALP(5),FUNL(5)
      DATA ZERO,TWO,THREE       /0.0D+00,2.0D+00,3.0D+00/
      DATA TOLALP,TOLFUN        /1.0D-01,1.5D-01/
      DATA FALMIN,FALNOR,FALMAX /1.0D-01,2.0D+00,5.0D+00/
      DATA DALNOR,DALMAX        /0.9D+00,0.5D+00/
      DATA ONEPT5               /1.5D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE.
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
c
*      write(6,*) ' Step in internals'
*      call output(ds,1,nzvar,1,1,nzvar,1,1)
      gs0 = ddot(nzvar,ds,1,g,1)
**      write(6,*) ' HERE GS0 ', gs0
*     rescale internal step into appropriate units
      call geom_zmat_ico_scale(ds,nzvar,
     $     0.52917715d0, 0.52917715d0*180d0/(4d0*atan(1d0)))
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
         G (I) =EG(I)           ! Overwrites int. grad with cart. grad
      ENDDO
      DUM=ZERO
      DO I=1,NCART 
        DUM=DUM+G0(I)*DX(I)
      ENDDO
**      write(6,*) ' HERE DUM ', dum
c$$$      GS0=DUM
      call geom_hnd_get_data('zmat',s,nzvar)
c
           NLPTS  = 1
       ALP(NLPTS) = ZERO
      FUNL(NLPTS) = FUNC0
C
C     ----- SEARCH ALONG -DX-DIRECTION. -----
C
  200 CONTINUE
*      write(6,*) ' HERE alpha ', alph
      IF(OUT) THEN
         WRITE(IW,9991) ALPH
      ENDIF
*      DO I = 1,NCART 
*         X(I) = X0(I)+ALPH*DX(I)
*         C(I) = X(I)
*      ENDDO
      NPTS=NPTS+1
      IF(NPTS.GT.NPTOPT) THEN
         EXCESS=.TRUE.
         goto 8080              ! return
      ENDIF
C
**      do alph = -0.1d0,0.1d0,0.05d0
      do i = 1, nzvar
         news(i) = s(i) + alph*ds(i)
      enddo
      status = geom_cart_coords_set(geom,x0)
      call opt_internal_to_cartesian(rtdb,geom,news,nzvar,nvar,ncart)
      if (.not. GEOM_RTDB_STORE(RTDB,GEOM,'geometry'))
     $     call errquit('line search internals: geom?',0)
C
C
C     ----- CALL FUNCTION EVALUATION -----
C
         STATUS=UTIL_TEST_TIME_REMAINING(RTDB,INT(ENERGY_TIME*ONEPT5))
         IF(OUT) THEN
            WRITE(IW,*) 'ENOUGH TIME REMAINING FOR ENERGY = ? ',STATUS
         ENDIF
         IF(.NOT.STATUS) THEN
            if (ga_nodeid().eq.0) 
     $           WRITE(IW,*) 'Not enough time remaining,',
     $           ' shutting down ... '
            NPTOPT=NPTS-1
            EXCESS=.TRUE.
            goto 8080           ! return
         ENDIF
C
      CALL HND_OPT_ENERGY(RTDB,GEOM)
      FUNC=E
**      write(6,*) ' HERE ',alph,e
**      enddo
**      call errquit('done',0)
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
         STATUS=UTIL_TEST_TIME_REMAINING(RTDB,INT(GRADIENT_TIME*ONEPT5))
         IF(OUT) THEN
            WRITE(IW,*) 'ENOUGH TIME REMAINING FOR GRADIENT = ? ',STATUS
         ENDIF
         IF(.NOT.STATUS) THEN
            if (ga_nodeid().eq.0) 
     $           WRITE(IW,*) 'Not enough time remaining,',
     $           ' shutting down ... '
            NPTOPT=NPTS-1
            EXCESS=.TRUE.
            goto 8080
         ENDIF
C
      CALL HND_OPT_GRADIENT(RTDB,GEOM)
      DO I=1,NCART 
         G(I) =EG(I)
      ENDDO
      status = geom_cart_coords_get(geom,x) ! Get latest cartesians
      status = geom_cart_coords_get(geom,c) ! 
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
      IF((GS.GT.GS0).OR.RETRY) THEN ! Needs modifying for internal search????
         GO TO 1000
      ENDIF
      RETRY = .TRUE.
      ACUBE = (-TWO*(FUNC-FUNC0)+(GS+GS0)*ALPH)/ALPH**3
      BCUBE = -ACUBE*ALPH+(FUNC-FUNC0-GS0*ALPH)/ALPH**2
      DSCRM =  BCUBE*BCUBE-THREE*ACUBE*GS0
      IF(DSCRM.LT.ZERO) THEN
         goto 8080
      ENDIF
      ALPH  = (-BCUBE+ SQRT(DSCRM))/(THREE*ACUBE)
*     RJH cannot do the following with internals (doubly so with
*     constants/constraints) unless first transform gradient to
*     internal space to compute gs0 ... might as well do a full 
*     geometry optimization step.
c$$$      IF((ALPH.LT.ALPHA).AND.(ALPH.GT.ZERO)) THEN
c$$$         write(6,*) ' Going to 100 '
c$$$         GO TO 100
c$$$      ENDIF
C
C     ----- RETURN AFTER SUCCESSFUL SEARCH -----
C
 1000 CONTINUE
      ALPHA=ALPH
      IF(DBUG) THEN
         WRITE(IW,9993)
      ENDIF
c
c     return statement for every exit so can unscale internal step
c
 8080 continue
*     undo rescale internal step into appropriate units
      call geom_zmat_ico_scale(ds,nzvar,
     $     1d0/0.52917715d0, 4d0*atan(1d0)/(0.52917715d0*180d0))
c
      RETURN
c
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
#include "rtdb.fh"
      INTEGER RTDB
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
      if (.not. ma_push_get(mt_dbl,need,'mem update',i_ifirst,ifirst))
     $     call errquit('hnd_opt_hss_update: ma ', need)
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
         call geom_hnd_get_data('zind',x(j20),ndim*ndim)
         call geom_hnd_get_data('drv.hess',x(j10),ndim*ndim)
c
         if(dbug) then
            write(iw,*) 'in hss_update, current hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
            write(iw,*) 'in hss_update, zind = '
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
         call geom_hnd_get_data('drv.hess',x(i10),ndim*ndim)
         if(dbug) then
            if( ga_nodeid().eq.0) then
               write(iw,*) 'in hss_update, cartesian hessian = '
               call hnd_prsq(x(i10),ndimx,ndimx,ndim)
            endif
         endif
         call ga_sync()
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
         call geom_hnd_get_data('zind',x(j20),ndim*ndim)
c
         call hnd_tfhsi(x(j10),x(j20),x(j30),ndimq,ndimq,ndim)
c
         call geom_hnd_put_data('drv.hess',x(j10),ndim*ndim)
         if(out) then
            write(iw,*) 'in hss_update, ',
     $                  'back-transformed updated hessian = '
            call hnd_prsq(x(j10),ndimq,ndimq,ndim)
         endif
      else
         call geom_hnd_put_data('drv.hess',x(i10),ndim*ndim)
         if(dbug) then
            write(iw,*) 'in hss_update, updated hessian = '
            call hnd_prsq(x(i10),ndimx,ndimx,ndim)
         endif
      endif
c
c     ----- release memory -----
c
      if (.not. ma_pop_stack(i_ifirst))
     $     call errquit('hnd_opt_hss_update: ma pop ', 0)

c
      if(dbug) then
         write(iw,*) '  end of hnd_opt_hss_update'
      endif
c
      RETURN
      END
      SUBROUTINE HND_OPT_HSS_UPDATE_X(HESS,HESSDX,NVAR,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "coptopt.fh"
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
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      DIMENSION HESS(NDIM,*),HESSDX(*)
      DATA ZERO /0.0D+00/
      DATA TOL  /1.0D-08/
C
      DBUG=.FALSE.
      OUT =.FALSE.
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
      DO IZ=1,NZFRZ
         I=IZFRZ(IZ)     
         G(I)=G0(I)
         DO J=1,NVAR
            HESS(I,J)=ZERO
            HESS(J,I)=ZERO
         ENDDO
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
