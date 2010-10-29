      logical function drv_lst(rtdb)
*
* $Id$
*
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
      dbug=.false.
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
#include "errquit.fh"
#include "rtdb.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "geom.fh"
#include "util.fh"
      logical       geom_zmt_get_nzvar
      logical       geom_zmt_get_nizmat
      external      geom_zmt_get_nizmat
      logical       geom_zmt_get_nzfrz 
      external      geom_zmt_get_nzfrz 
      logical       geom_zmt_get_izfrz 
      external      geom_zmt_get_izfrz 
      logical       geom_lst_get_coord
      external      geom_lst_get_coord
      integer       rtdb
      integer       geom
      logical       do_gradient
      logical       some
      logical       out
      logical       dbug
      logical       status
      logical       done
      logical       rstart
      logical       zcoord
      character*16  atmlab
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c                mxcoor=max(mxcart,mxzmat)
      common/hnd_iofile/ir,iw
      common/hnd_optmiz/x0(mxcoor),x(mxcoor),dx(mxcoor),
     1                  g0(mxcoor),g(mxcoor),ds(mxcoor),
     2                  func,func0,gmax,gmax0,curv,alpha,gnorm
      common/hnd_optfrz/nzfrz,izfrz(mxcoor),iatfrz(mxatom)
      common/hnd_optvar/zcoord,ncoord,mcoord
      common/hnd_optfun/e,eg(mxcart)
      common/hnd_lstrun/nptlst,iptlst
      common/hnd_lsttim/energy_time,gradient_time
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_molxyz/c(mxcart),zan(mxatom),nat
      common/hnd_mollab/atmlab(mxatom)
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
c
      hnd_lstx = .false.
c
      do_gradient = .false.
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.
      some=some.or.out
      if (ga_nodeid().eq.0.or.some)  then
         write(iw,9999)
      endif
c
c     ----- get going ... -----
c
      if (.not. geom_create(geom, 'geometry'))
     &       call errquit('hnd_lst: geom_create?', 911, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &       call errquit('hnd_lst: no geometry ', 911, RTDB_ERR)
      if (.not. geom_cart_get(geom, ncent, atmlab, c, zan))
     &       call errquit('hnd_lst: geom_get ', 911, GEOM_ERR)
       nat = ncent
      if(.not.geom_zmt_get_nizmat(geom,nzmat))
     &       call errquit('geom_input: geom_zmt_get_nizmat failed',0,
     &       GEOM_ERR)
      if(.not.geom_zmt_get_nzvar(geom,nzvar))
     &     call errquit('geom_input: geom_zmt_get_nzvar failed',0,
     &       GEOM_ERR)
c
c     ----- scan -lst- coordinates file -----
c
      nptlst=1
   10 continue
         status=geom_lst_get_coord(x(1),x(1+nat),x(1+nat),
     &                             nat,nptlst)
         if(status) then
            nptlst=nptlst+1
            go to 10
         else
            nptlst=nptlst-1
         endif
c
      if (ga_nodeid().eq.0.or.dbug)  then
         write(iw,9996) nptlst
      endif
      if(nptlst.eq.0) then
         write(iw,*) 'no -lst- points found. return and continue ..'
         return
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
c     ----- frozen coordinates ? if so .... -----
c
      nzfrz=0
      if (geom_zmt_get_nzfrz(geom,nzfrz)) then
         if(.not.geom_zmt_get_izfrz(geom,izfrz,nzfrz)) 
     $        call errquit('geom_input: geom_zmt_get_izfrz  failed',0,
     &       GEOM_ERR)
      endif
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
c
c     if internal coordinates ...
c           nzvar  = # of (redundant) internal coordinates
c           ncoord = # of independent internal coordinates
c           mcoord = max( 3*nat , nzvar )
c     if cartesian coordinates ...
c           ncoord = 3*nat
c           mcoord = 3*nat
c
      rstart     =.false.
c
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
         call geom_bandbi(geom)
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
          x(i)=zero
          g(i)=zero
      enddo
c
c     ----- first point -----  
c
      iptlst=1
      do i=1,ncart
         x(i)=c(i)
      enddo
c
c     ----- energy -----
c
         time_start=util_wallsec()
      call hnd_lst_energy(rtdb,geom)
         time_end  =util_wallsec()
         energy_time=time_end-time_start
         if(out) then
            write(iw,*) 'energy_time = ',energy_time
         endif
c
      func=e
      if(do_gradient) then
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
      endif
c
      call hnd_lst_print(rtdb)
      do i=1,ncart
         x0(i)=x(i)
      enddo
      do i=1,ncoord
         g0(i)=g(i)
      enddo
c
      if (ga_nodeid().eq.0.and.dbug)  then
         if (.not. geom_print(geom)) call errquit
     $        ('hnd_lst_drv: geom_print?',0, GEOM_ERR)
         if (util_print('bonds',print_default)) then
           if (.not.geom_print_distances(geom)) call errquit(
     &         'hnd_lst_drv: geom_print_distances failed',911,
     &       GEOM_ERR)
         endif
         if (util_print('angles',print_default)) then
           if (.not.geom_print_angles(geom)) call errquit(
     &         'hnd_lst_drv: geom_print_angles failed',911, GEOM_ERR)
         endif
      endif
      call ga_sync()
c
c     ----- lst path -----
c
  100 continue
      iptlst=iptlst+1
c
      call hnd_lst_path(rtdb,geom,rstart)
c
      if(zcoord) then
         call hnd_opt_tfgx(rtdb,geom)
      endif
      call hnd_lst_print(rtdb)
c
      if (ga_nodeid().eq.0.and.dbug)  then
         if (.not. geom_print(geom)) call errquit
     $        ('hnd_lst_drv: geom_print?',0, GEOM_ERR)
         if (util_print('bonds',print_default)) then
           if (.not.geom_print_distances(geom)) call errquit(
     &         'hnd_lst_drv: geom_print_distances failed',911, GEOM_ERR)
         endif
         if (util_print('angles',print_default)) then
           if (.not.geom_print_angles(geom)) call errquit(
     &         'hnd_lst_drv: geom_print_angles failed',911,
     &       GEOM_ERR)
         endif
      endif
      call ga_sync()
c
c     ----- go back to next -lst- point -----
c
      done=iptlst.eq.nptlst
      if(.not.done) then
         go to 100
      else
         if (ga_nodeid().eq.0.or.dbug)  then
            write(iw,9998)
         endif
         call ga_sync()
         if (.not.geom_destroy(geom))
     &      call errquit('hnd_lst: geom_destroy?',911, GEOM_ERR)
         hnd_lstx=.true. 
      endif       
c
 9999 format(/,10x,8(1h-),
     1       /,10x,'LST Path',
     2       /,10x,8(1h-))
 9998 format(/,10x,18(1h-),
     1       /,10x,'LST Path completed',
     2       /,10x,18(1h-))
 9997 format(' in hnd_lstx, nzvar,nvar,ncart,ncoord,mcoord = ',5i5)
 9996 format(' number of points on -lst- path (nptlst) = ',i4)    
 9993 format(' number of frozen internal coordinates = ',i4)
 9992 format('           frozen internal coordinates = ',12i4)
 9991 format('        no frozen internal coordinates.')
      end
      SUBROUTINE HND_LST_PRINT(RTDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      INTEGER   RTDB
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL      DBUG
      CHARACTER*16 ATMNAM
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLLAB/ATMNAM(MXATOM)
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,CURV,ALPHA,GNORM
      COMMON/HND_LSTRUN/NPTLST,IPTLST
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      COMMON/HND_OPTFUN/E,EG(MXCART)
C
      DBUG=.FALSE.
C
C     ----- PRINT OPTIMIZATION SUMMARY -----
C
      IF( GA_NODEID().EQ.0.OR.DBUG) THEN
         WRITE(IW,9999) NPTLST,IPTLST,FUNC       
      ENDIF
      RETURN
 9999 FORMAT(1H1,/,1X,
     1       'nptlst  iptlst      func             ',1X,I5,I8,F17.8)
      END
      SUBROUTINE HND_LST_PATH(RTDB,GEOM,RSTART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
C
C     ----- ONE DIMENSIONAL SEARCH.  -----
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
C
      INTEGER   RTDB
      INTEGER   GEOM
      LOGICAL   STATUS
      LOGICAL   UTIL_TEST_TIME_REMAINING
      EXTERNAL  UTIL_TEST_TIME_REMAINING
      LOGICAL   GEOM_LST_GET_COORD
      EXTERNAL  GEOM_LST_GET_COORD
      CHARACTER*16 TAGS_NW
      DIMENSION COORDS_NW(MXCART)
      DIMENSION CHARGE_NW(MXATOM)
      DIMENSION   TAGS_NW(MXATOM)
C
      LOGICAL DO_GRADIENT
      LOGICAL RSTART
      LOGICAL DBUG
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(MXCART),ZAN(MXATOM),NAT         
      COMMON/HND_OPTMIZ/X0(MXCOOR),X(MXCOOR),DX(MXCOOR),
     1                  G0(MXCOOR),G(MXCOOR),DS(MXCOOR),
     2                  FUNC,FUNC0,GMAX,GMAX0,DELSQ,ALPHA,GNORM
      COMMON/HND_LSTRUN/NPTLST,IPTLST
      COMMON/HND_OPTFRZ/NZFRZ,IZFRZ(MXCOOR),IATFRZ(MXATOM)
      COMMON/HND_OPTTIM/ENERGY_TIME,GRADIENT_TIME
      COMMON/HND_OPTFUN/E,EG(MXCART)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      DATA ZERO,TWO,THREE       /0.0D+00,2.0D+00,3.0D+00/
      DATA ONEPT5               /1.5D+00/
C
      DO_GRADIENT=.FALSE.
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(OUT) THEN
         WRITE(IW,9994)
      ENDIF
C
      NCART =3*NAT
C
C     ----- NORMAL SETUP -----
C
      FUNC0 = FUNC
      DO I = 1,NCART 
         X0(I) = X(I)
         G0(I) =EG(I)
         G (I) =EG(I)
         X (I) =ZERO
      ENDDO
C
C     ----- NEW -LST- POINT -----
C
      STATUS=GEOM_LST_GET_COORD(X(1),X(1+NAT),X(1+NAT*2),
     1                          NAT,IPTLST)
      DO IAT = 1,NAT
         C(1+3*(IAT-1)) = X(IAT)
         C(2+3*(IAT-1)) = X(IAT+NAT)
         C(3+3*(IAT-1)) = X(IAT+NAT*2)
      ENDDO
      DO I = 1,NCART
         X(I)=ZERO
         X(I)=C(I)
      ENDDO
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
         WRITE(IW,9998) (        X(I),I=1,NCART)
         WRITE(IW,9999)
         WRITE(IW,9998) (COORDS_NW(I),I=1,NCART)
      ENDIF
      DO I = 1,NCART
         COORDS_NW(I)=X(I)
      ENDDO
      STATUS=GEOM_CART_SET(GEOM,NAT_NW,TAGS_NW,COORDS_NW,
     1                                         CHARGE_NW)
      STATUS=GEOM_RTDB_STORE(RTDB,GEOM,'geometry')
C
      IF(OUT) THEN
         WRITE(IW,9999)
         WRITE(IW,9998) (COORDS_NW(I),I=1,NCART)
      ENDIF
C
C
C     ----- CALL FUNCTION EVALUATION -----
C
         STATUS=UTIL_TEST_TIME_REMAINING(RTDB,INT(ENERGY_TIME*ONEPT5))
         IF(OUT) THEN
            WRITE(IW,*) 
     1      'ENOUGH TIME REMAINING FOR ENERGY = ? ',STATUS
         ENDIF
         IF(.NOT.STATUS) THEN
            WRITE(IW,*) 
     1      'NOT ENOUGH TIME REMAINING, SHUTTING DOWN ... '
            RETURN
         ENDIF
C
      CALL HND_LST_ENERGY(RTDB,GEOM)
      FUNC=E
C
C     ----- CALCULATE THE GRADIENT FOR THE FINAL POINT -----
C
  900 CONTINUE
      IF(DO_GRADIENT) THEN
         STATUS=UTIL_TEST_TIME_REMAINING(RTDB,INT(GRADIENT_TIME*ONEPT5))
         IF(OUT) THEN
            WRITE(IW,*) 
     1      'ENOUGH TIME REMAINING FOR GRADIENT = ? ',STATUS
         ENDIF
         IF(.NOT.STATUS) THEN
            WRITE(IW,*) 
     1      'NOT ENOUGH TIME REMAINING, SHUTTING DOWN ... '
            RETURN
         ENDIF
C
         CALL HND_OPT_GRADIENT(RTDB,GEOM)
         DO I=1,NCART 
            G(I) =EG(I)
         ENDDO
      ENDIF
C
C     ----- RETURN AFTER SUCCESSFUL POINT -----
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
 9995 FORMAT(' IN SEARCH_LIN,  X = ')
 9994 FORMAT(' STARTING SEARCH_LIN ')
 9993 FORMAT('   ENDING SEARCH_LIN ')
      END
      subroutine hnd_lst_energy(rtdb,geom)
      implicit double precision (a-h,o-z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
      integer  rtdb
      integer  geom
      logical  task_energy
      external task_energy
      logical  status
      logical  dbug
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
      common/hnd_iofile/ir,iw
      common/hnd_lstrun/nptlst,iptlst
      common/hnd_optfun/e,eg(mxcart)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_zmtpar/nzmat,nzvar,nvar
c
      logical geom_print_zmatrix
c
      dbug=.false.
c
      if( ga_nodeid().eq.0.or.dbug) then
         write(iw,9999) nptlst,iptlst
         if(dbug) then
            write(iw,9996)
            do iat=1,nat
               write(iw,9995) iat,(c(i,iat),i=1,3)
            enddo
            write(iw,9997)
         endif
         if(.not.geom_print(geom)) 
     1     call errquit('hnd_lst_energy: print error',911, GEOM_ERR)
      endif
c
      if(nzmat.gt.0) then
         call geom_bandbi(geom)
         if (.not. geom_print_zmatrix(geom,0d0,.false.))
     $        call errquit('hnd_lst: print zmat failed',0, GEOM_ERR)
         write(iw,9997)
      endif
c
      if (task_energy(rtdb)) then
         status=rtdb_get(rtdb,'task:energy',MT_DBL,1,e)
         if(dbug) then
            write(iw,9998) e
         endif
      else
         call errquit('lst path: energy failed', 0, GEOM_ERR)
      endif
c
      return
 9999 format(1h1,' nptlst',i3,//,'  iptlst',I3)
 9998 format(' in lst_energy, e = ',f15.10)
 9997 format(/)
 9996 format(21x,21(1h-),/,
     1       21x,'cartesian coordinates',/,
     2       21x,21(1h-))
 9995 format(9x,i5,3f15.8) 
      end
