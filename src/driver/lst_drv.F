      logical function drv_lst(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer  rtdb
      logical  hnd_lstx
      external hnd_lstx
      integer  ir, iw
      logical  some, dbug
      common/hnd_iofile/ir,iw
c
      call util_print_push
      call util_print_rtdb_load(rtdb, 'driver')
      call ecce_print_module_entry('driver')
      if (util_print('brdcst', print_never)) call setdbg(1)
c
      dbug=.true.
      some=.false.
      some=some.or.dbug
      ir=LuIn
      iw=LuOut  
      if(some) then
         write(iw,*) ' drv_lst: calling hnd_lstx . '
      endif
c
c     ----- lst path -----
c
      drv_lst = hnd_lstx(rtdb)
c
      if (drv_lst) then
         call ecce_print_module_exit('driver', 'ok')
      else
         call ecce_print_module_exit('driver', 'failed')
      endif
c
      if (util_print('brdcst', print_never)) call setdbg(0)
      call util_print_pop
c
      end
      logical function hnd_lstx(rtdb)
      implicit double precision (a-h,o-z)
#include "rtdb.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "geom.fh"
#include "util.fh"
      logical       geom_zmt_get_nizmat
      logical       geom_usr_get_nzfrz 
      logical       geom_usr_get_izfrz 
      integer       rtdb
      integer       geom
      integer       idriver
      logical       hnd_sadx
      logical       some
      logical       out
      logical       dbug
      logical       status
      logical       cvged
      logical       zcoord
      logical       rstart
      logical       excess
      character*255 opt_hess_fil
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
      common/hnd_optopt/cvgopt,nptopt,inhess,linopt,modupd,moddir,
     1                                                     modsad,
     2                  convggm,convgg,convge,maxiter
      common/hnd_optsad/sadstp
      common/hnd_optfrz/nzfrz,izfrz(mxcoor),iatfrz(mxatom)
      common/hnd_optfun/e,eg(mxcart)
      common/hnd_optvar/zcoord,ncoord,mcoord
      common/hnd_optfil/opt_hess_fil
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
c
      hnd_lstx = .false.
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.
      some=some.or.out
      if (ga_nodeid().eq.0.or.some)  then
         write(iw,9999)
      endif
c--
      if(some.or..not.some) then
         hnd_lstx = .true.  
         return
      endif
c--
c
c     ----- parameters for optimization driver -----
c
      if (.not. rtdb_get(rtdb,'driver:cvgopt',mt_dbl,1,cvgopt)) then
           cvgopt=8.0d-04
           if(dbug) then
              write(iw,*) 'using default value of -cvgopt-',cvgopt
           endif
      endif
      if (.not. rtdb_get(rtdb,'driver:nptopt',mt_int,1,nptopt)) then
           nptopt=20
           if(dbug) then
              write(iw,*) 'using default value of -nptopt-',nptopt
           endif
      endif
      if (.not. rtdb_get(rtdb,'driver:inhess',mt_int,1,inhess)) then
           inhess=0
           if(dbug) then
              write(iw,*) 'using default value of -inhess-',inhess
           endif
      endif
      if (.not. rtdb_get(rtdb,'driver:linopt',mt_int,1,linopt)) then
           linopt=10
           if(dbug) then
              write(iw,*) 'using default value of -linopt-',linopt
           endif
      endif
      if (.not. rtdb_get(rtdb,'driver:modupd',mt_int,1,modupd)) then
           modupd=1
           if(dbug) then
              write(iw,*) 'using default value of -modupd-',modupt
           endif
      endif
      if (.not. rtdb_get(rtdb,'driver:moddir',mt_int,1,moddir)) then
           moddir=1
           if(dbug) then
              write(iw,*) 'using default value of -moddir-',moddir
           endif
      endif
      if (.not. rtdb_get(rtdb,'driver:modsad',mt_int,1,modsad)) then
           modsad=0
           if(dbug) then
              write(iw,*) 'using default value of -modsad-',modsad
           endif
      else
           if(dbug) then
              write(iw,*) 'using -rtdb-  value of -modsad-',modsad
           endif
      endif
      if (.not. rtdb_get(rtdb,'driver:sadstp',mt_dbl,1,sadstp)) then
           sadstp=0.35d+00
           if(dbug) then
              write(iw,*) 'using default value of -sadstp-',sadstp
           endif
      endif
c
      if (ga_nodeid().eq.0.or.dbug)  then
         write(iw,9996) cvgopt,nptopt,linopt,inhess,modupd,
     1                                       modsad,moddir,
     2                                       sadstp
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
         hnd_lstx=hnd_sadx(rtdb,idriver)
         return
      else
         if (ga_nodeid().eq.0.or.some)  then
            write(iw,9994)
         endif
      endif
c
c     ----- get going ... -----
c
      if (.not. geom_create(geom, 'geometry'))
     &       call errquit('hnd_opt: geom_create?', 911)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &       call errquit('hnd_opt: no geometry ', 911)
      if (.not. geom_cart_get(geom, ncent, atmlab, c, zan))
     &       call errquit('hnd_opt: geom_get ', 911)
       nat = ncent
      if(.not.geom_zmt_get_nizmat(geom,nzmat))
     &       call errquit('geom_input: geom_zmt_get_nizmat failed',0)
c
c     ----- frozen coordinates ? if so .... -----
c
      status = rtdb_parallel(.false.)
      if(ga_nodeid().eq.0) then
         if(.not.geom_usr_get_nzfrz(geom,nzfrz)) then
            if(dbug) then
               write(iw,*) '-nzfrz- not found in -geomP- object'
            endif
            nzfrz=0
         else
            if(dbug) then
               write(iw,*) '-nzfrz-     found in -geomP- object',nzfrz
               write(iw,*) 'now read -izfrz-                   '
            endif
            if(.not.geom_usr_get_izfrz(geom,izfrz,nzfrz)) then
               call errquit('geom_input: geom_usr_get_izfrz  failed',0)
            endif
            if(dbug) then
               if(nzfrz.gt.0) then
                  write(iw,*) 'izfrz = ',(izfrz(i),i=1,nzfrz)
               endif
            endif
         endif
         call ga_brdcst(24+msgint,nzfrz,mitob(1     ),0)
         call ga_brdcst(25+msgint,izfrz,mitob(mxcoor),0)
      else
         call ga_brdcst(24+msgint,nzfrz,mitob(1     ),0)
         call ga_brdcst(25+msgint,izfrz,mitob(mxcoor),0)
      endif
      call ga_sync()
      status = rtdb_parallel(.true.)
      if (ga_nodeid().eq.0.or.dbug)  then
         if(nzfrz.gt.0) then
            write(iw,9993) nzfrz
            write(iw,9992) (izfrz(i),i=1,nzfrz) 
         else
            write(iw,9991)
         endif
      endif
c
c     ----- get ready to start now -----
c     if internal coordinates ...
c           nzvar  = # of (redundant) internal coordinates
c           ncoord = # of independent internal coordinates
c           mcoord = max( 3*nat , nzvar )
c     if cartesian coordinates ...
c           ncoord = 3*nat
c           mcoord = 3*nat
c
      rstart=.false.
      excess=.false.
      zcoord=.false.
      zcoord=nzmat.gt.0
      if(dbug)  then
         write(iw,*) 'nzmat,zcoord = ',nzmat,zcoord    
         write(iw,*) 'nzfrz        = ',nzfrz
         if(nzfrz.gt.0) then
            write(iw,*) 'izfrz = ',(izfrz(i),i=1,nzfrz)
         endif
      endif
      if(zcoord) then
         call geom_zmt_icoord(rtdb,geom)
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
      if(dbug)  then
         write(iw,9997) nzvar,nvar,ncart,ncoord,mcoord
      endif
c
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
      if(cvged) then
         if (ga_nodeid().eq.0.or.dbug)  then
            write(iw,9998)
            if (.not. geom_print(geom)) call errquit
     $           ('hnd_opt_drv: geom_print?',0)
*rak:            if (util_print('bonds',print_default))
*rak:     $           call distan(nat,c)
*rak:            if (util_print('angles',print_default))
*rak:     $           call angle(nat,c)
            if (util_print('bonds',print_default)) then
              if (.not.geom_print_distances(geom)) call errquit(
     &            'hnd_opt_drv: geom_print_distances failed',911)
            endif
            if (util_print('angles',print_default)) then
              if (.not.geom_print_angles(geom)) call errquit(
     &            'hnd_opt_drv: geom_print_angles failed',911)
            endif
         endif
         call ga_sync()
         if (.not.geom_destroy(geom))
     &      call errquit('hnd_opt: geom_destroy?', 911)
         hnd_lstx=.true.
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
         call hnd_opt_tfds(rtdb,geom)               
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
            write(iw,9995) nptopt
         endif
         if (.not.geom_destroy(geom))
     &      call errquit('hnd_opt: geom_destroy?', 911)
         call ga_sync()
         hnd_lstx=.false.
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
      if(cvged) then
         if (ga_nodeid().eq.0.or.dbug)  then
            write(iw,9998)
            if (.not. geom_print(geom)) call errquit
     $           ('hnd_opt_drv: geom_print?',0)
*rak:            if (util_print('bonds',print_default))
*rak:     $           call distan(nat,c)
*rak:            if (util_print('angles',print_default))
*rak:     $           call angle(nat,c)
            if (util_print('bonds',print_default)) then
              if (.not.geom_print_distances(geom)) call errquit(
     &            'hnd_opt_drv: geom_print_distances failed',911)
            endif
            if (util_print('angles',print_default)) then
              if (.not.geom_print_angles(geom)) call errquit(
     &            'hnd_opt_drv: geom_print_angles failed',911)
            endif
         endif
         if (.not.geom_destroy(geom))
     &      call errquit('hnd_opt: geom_destroy?', 911)
         call ga_sync()
         hnd_lstx=.true.
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
 9999 format(/,10x,8(1h-),
     1       /,10x,'LST Path',
     2       /,10x,8(1h-))
 9998 format(/,10x,18(1h-),
     1       /,10x,'LST Path completed',
     2       /,10x,18(1h-))
 9997 format(' in hnd_lstx, nzvar,nvar,ncart,ncoord,mcoord = ',5i5)
 9996 format(' optimization convergence threshold    (cvgopt) = ',f10.6,
     1     /,' maximum number of energy evaluations  (nptopt) = ',i4,/,
     2       ' max.num. of points in linear search   (linopt) = ',i4,/,
     3       ' initial guess option for opt.hessian  (inhess) = ',i4,/,
     4       ' hessian update selection option       (modupd) = ',i4,/,
     5       ' optimization to saddle point          (modsad) = ',i4,/,
     6       ' initial eigen-mode followed to saddle (moddir) = ',i4,/,
     7       ' maximum step size to saddle           (sadstp) = ',f6.2)
 9995 format(' maximum number of energy points reached. stop',
     1       ' currently -nptopt- = ',i5)
 9994 format(/,10x,19(1h-),
     1       /,10x,'Energy Minimization',
     2       /,10x,19(1h-),/)
 9993 format(' number of frozen internal coordinates = ',i4)
 9992 format('           frozen internal coordinates = ',12i4)
 9991 format('        no frozen internal coordinates.')
      end
