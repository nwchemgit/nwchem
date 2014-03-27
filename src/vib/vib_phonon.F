
*     *********************************************
*     *                                           *
*     *              vib_phonon                   *
*     *                                           *
*     *********************************************
      subroutine vib_phonon(rtdb,nion,rion,nion3,hess,mass)
      implicit none
      integer rtdb,nion
      real*8  rion(3,nion)
      integer nion3
      real*8  hess(nion3,nion3)
      real*8  mass(nion3)

#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "errquit.fh"

      double precision ams, wave
      integer numans
      COMMON /cvib_SETCON/ AMS(36),wave,numans      ! setup parameters

*     **** local variables ****
      logical value,oprint
      integer dosgrid(3),isize,i1,i2,i3,i
      real*8  omega,unita(3,3),unitg(3,3),lattice(6)
      real*8  xx,yy,zz,xxx,yyy,zzz,ks(3),kq(3)
      complex*16 cfac

      integer weight_dos(2),eigs_dos(2),dmat(2),wp(2),hessadjust(2)
      integer ierr,lwork,work(2),rwork(2),w2(1)

      write(*,*) "calculate phonon spectrum"

      oprint = util_print('phonon',print_medium)

      if (oprint) then
         call util_print_centered(luout,
     >        'NWChem Phonon Spectra Analysis ',
     >        40,.true.)
         write(luout,*)
      end if


*     **** get lattice vectors, reciprocal lattice vectors ****
      call vib_getlattice(rtdb,omega,unita,unitg,lattice)

      if (oprint) then
        write(luout,*)
        write(luout,*) '     Lattice Parameters '
        write(luout,*) '     ------------------ '
        write(luout,*)
        write(luout,5) 
 5    format('      lattice vectors in a.u.')
        write(luout,*)
        write(luout,1241) unita(1,1),
     >                    unita(2,1),
     >                    unita(3,1)
        write(luout,1242) unita(1,2),
     >                    unita(2,2),
     >                    unita(3,2)
        write(luout,1243) unita(1,3),
     >                    unita(2,3),
     >                    unita(3,3)
        write(luout,1232) lattice(1),
     >                    lattice(2),
     >                    lattice(3),
     >                    lattice(4),
     >                    lattice(5),
     >                    lattice(6)
        write(luout,1231) omega
        write(luout,*)
        write(luout,6)
 6      format('      reciprocal lattice vectors in a.u.')
        write(luout,*)
        write(luout,1244) unitg(1,1),
     >                    unitg(1,2),
     >                    unitg(1,3)
        write(luout,1245) unitg(2,1),
     >                    unitg(2,2),
     >                    unitg(2,3)
        write(luout,1246) unitg(3,1),
     >                    unitg(3,2),
     >                    unitg(3,3)


      end if


*     **** allocate  dos grid ****
      if (.not.rtdb_get(rtdb,'phonon:dos-grid',mt_int,3,dosgrid)) then
         dosgrid(1) = 5
         dosgrid(2) = 5
         dosgrid(3) = 5
      end if
      if (oprint) write(luout,1250) dosgrid

*     **** allocate memory for DOS calculation ****
      isize = dosgrid(1)*dosgrid(2)*dosgrid(3)*nion3
      value =           MA_push_get(mt_dbl,isize,
     >                       'weight_dos',weight_dos(2),weight_dos(1))
      value = value.and.MA_push_get(mt_dbl,isize,
     >                       'eigs_dos',eigs_dos(2),eigs_dos(1))

      value = value.and.MA_push_get(mt_dcpl,nion3*nion3,
     >                       'dmat',dmat(2),dmat(1))

      value = value.and.MA_push_get(mt_dbl,nion3*nion3,
     >                    'hessadjust',hessadjust(2),hessadjust(1))

      value = value.and.MA_push_get(mt_dbl,nion3,
     >                             'wpphonon',wp(2),wp(1))

      value = value.and.MA_push_get(mt_dbl,3*nion3,
     >                             'rworkphonon',rwork(2),rwork(1))

      lwork = nion3*nion3
      if (lwork.lt.(3*nion3)) lwork = 3*nion3
      value = value.and.MA_push_get(mt_dcpl,lwork,
     >                    'workphonon',work(2),work(1))
      value = value.and.MA_push_get(mt_dbl,nion3*nion3,
     >                    'w2phonon',w2(2),w2(1))

      if (.not.value)
     >   call errquit('vib_phonon:out of stack',0,MA_ERR)


      call vib_adjustfordynamic(unita,nion,rion,
     >                    nion3,mass,dbl_mb(rwork(1)),
     >                    dcpl_mb(work(1)),dbl_mb(w2(1)),
     >                    hess,dbl_mb(hessadjust(1)))

      value = MA_pop_stack(w2(2))
      if (.not.value)
     >   call errquit('vib_phonon:popping stack',0,MA_ERR)



      xxx = 1.0d0/(1.0d0*dosgrid(1))
      yyy = 1.0d0/(1.0d0*dosgrid(2))
      zzz = 1.0d0/(1.0d0*dosgrid(3))
      do i3=0,(dosgrid(3)-1)
      do i2=0,(dosgrid(2)-1)
      do i1=0,(dosgrid(1)-1)
          ks(1) = xx*xxx
          ks(2) = yy*yyy
          ks(3) = zz*zzz
          kq(1) = unitg(1,1)*ks(1)+unitg(1,2)*ks(2)+unitg(1,3)*ks(3)
          kq(2) = unitg(2,1)*ks(1)+unitg(2,2)*ks(2)+unitg(2,3)*ks(3)
          kq(3) = unitg(2,1)*ks(1)+unitg(2,2)*ks(2)+unitg(2,3)*ks(3)

          call vib_gendynamicmatrix(unita,kq,nion,rion,
     >                          nion3,dbl_mb(hessadjust(1)),
     >                          dcpl_mb(dmat(1)),dbl_mb(wp(1)),
     >                          dcpl_mb(work(1)),lwork,dbl_mb(rwork(1)))

      write(*,*) 
      write(*,*) "kvec=",kq," eigs="
      do i=1,nion3
         write(*,*) "w(",i,")=",dbl_mb(wp(1)+i-1),
     >               dsqrt(dabs(dbl_mb(wp(1)+i-1)))*wave
      end do
      write(*,*) 
      write(*,*) 

      end do
      end do
      end do
c      kq(1)=0.1d0
c      kq(2)=0.2d0
c      kq(3)=0.3d0
c      call vib_gendynamicmatrix(unita,kq,nion,rion,
c     >                          nion3,dbl_mb(hessadjust(1)),
c     >                          dcpl_mb(dmat(1)),dbl_mb(wp(1)),
c     >                          dcpl_mb(work(1)),lwork,dbl_mb(rwork(1)))



*     **** deallocate stack memory ****
      value =           MA_pop_stack(work(2))
      value = value.and.MA_pop_stack(rwork(2))
      value = value.and.MA_pop_stack(wp(2))
      value = value.and.MA_pop_stack(hessadjust(2))
      value = value.and.MA_pop_stack(dmat(2))
      value = value.and.MA_pop_stack(eigs_dos(2))
      value = value.and.MA_pop_stack(weight_dos(2))
      if (.not.value)
     >   call errquit('vib_phonon:popping stack',1,MA_ERR)

      write(*,*) "done phonon spectrum"
      return

 1231 FORMAT(5x,' omega=',f8.1)
 1232 FORMAT(5x,' a=    ',f8.3,' b=   ',f8.3,' c=    ',f8.3,
     >      /5x,' alpha=',f8.3,' beta=',f8.3,' gamma=',f8.3)
 1241 FORMAT(5x,' a1=<',3f8.3,' >')
 1242 FORMAT(5x,' a2=<',3f8.3,' >')
 1243 FORMAT(5x,' a3=<',3f8.3,' >')
 1244 FORMAT(5x,' b1=<',3f8.3,' >')
 1245 FORMAT(5x,' b2=<',3f8.3,' >')
 1246 FORMAT(5x,' b3=<',3f8.3,' >')
 1250 FORMAT(/5x,' dos grid=',3I4,' ')
      end

*     *************************************************
*     *                                               *
*     *              vib_adjustfordynamic             *
*     *                                               *
*     *************************************************

      subroutine vib_adjustfordynamic(unita,nion,rion,
     >                                nion3,mass,vc,w1,w2,
     >                                hess,hessout)
      implicit none
      real*8 unita(3,3)
      integer nion
      real*8  rion(3,nion)
      integer nion3
      real*8  mass(nion3),vc(nion3,3)
      real*8  w1(nion3,nion3),w2(nion3,nion3)
      real*8  hess(nion3,nion3),hessout(nion3,nion3)

*     **** local variables ****
      integer ii,jj,i,j,xyz,rst,ierr
      integer i1,i2,i3,dcount
      real*8  dx,dy,dz,x,y,z,d2,d2min,fac

      double precision ams, wave
      integer numans
      COMMON /cvib_SETCON/ AMS(36),wave,numans      ! setup parameters


      do jj=1,nion

*        **** diagonal term ****
         do xyz=1,3
         do rst=1,3
           i = 3*(jj-1)+xyz
           j = 3*(jj-1)+rst
           hessout(i,j) = hess(i,j)
         end do
         end do

         do ii=jj+1,nion
            dx = rion(1,jj)-rion(1,ii)
            dy = rion(1,jj)-rion(1,ii)
            dz = rion(1,jj)-rion(1,ii)

            d2min=9.99d12
            do i3=-1,1
            do i2=-1,1
            do i1=-1,1
               x = dx + unita(1,1)*i1 + unita(1,2)*i2 + unita(1,3)*i3
               y = dy + unita(2,1)*i1 + unita(2,2)*i2 + unita(2,3)*i3
               z = dz + unita(3,1)*i1 + unita(3,2)*i2 + unita(3,3)*i3
               d2 = x**2 + y**2 + z**2
               if (d2.lt.d2min) d2min = d2
            end do
            end do
            end do

            dcount = 0
            do i3=-1,1
            do i2=-1,1
            do i1=-1,1
               x = dx + unita(1,1)*i1 + unita(1,2)*i2 + unita(1,3)*i3
               y = dy + unita(2,1)*i1 + unita(2,2)*i2 + unita(2,3)*i3
               z = dz + unita(3,1)*i1 + unita(3,2)*i2 + unita(3,3)*i3
               d2 = x**2 + y**2 + z**2
               if (dabs(d2-d2min).lt.1.0d-4) dcount = dcount + 1
            end do
            end do
            end do
            fac = 1.0d0/dble(dcount)
            do xyz=1,3
            do rst=1,3
               i = 3*(ii-1)+xyz
               j = 3*(jj-1)+rst
               hessout(i,j) = fac*hess(i,j)
               hessout(j,i) = fac*hess(j,i)
            end do
            end do

            write(*,*) ii,jj,dcount

         end do !**ii**

      end do !**jj**

*     **** remove translations - hessout = (I-vc*vc')*hessout*(I-vc*vc') ****
      call dcopy(3*nion3,0.0d0,0,vc,1)
      x = 1.0d0/dsqrt(dble(nion))
      do i=1,3 ! unit translation vector for x-, y-, and z-directions
        do ii=1,nion
           vc(3*(ii-1)+i,i) = x
        end do
      end do
      call dcopy(nion3*nion3,0.0d0,0,w1,1)
      do i=1,nion3
         w1(i,i) = 1.0d0
      end do
      call dgemm('n','t',nion3,nion3,3,-1.0d0,
     >           vc,nion3,
     >           vc,nion3,
     >           1.0d0,w1,nion3)
      call dgemm('n','n',nion3,nion3,nion3,1.0d0,
     >           w1,nion3,
     >           hessout,nion3,
     >           0.0d0,w2,nion3)
      call dgemm('n','n',nion3,nion3,nion3,1.0d0,
     >           w2,nion3,
     >           w1,nion3,
     >           0.0d0,hessout,nion3)


*     **** divide by sqrt of masses ****
      do jj=1,nion
         do ii=1,nion
            fac = 1.0d0/dsqrt(mass(ii)*mass(jj))
            do xyz=1,3
            do rst=1,3
               i = 3*(ii-1)+xyz
               j = 3*(jj-1)+rst
               hessout(i,j) = hessout(i,j)*fac
            end do
            end do
         end do
      end do

      !*** call Gainsville matrix output ***
      write(*,*) 
      write(*,*) "hess="
      call output(hess,1,nion3,1,nion3,nion3,nion3,1)
      write(*,*) 
      write(*,*) "hessout="
      call output(hessout,1,nion3,1,nion3,nion3,nion3,1)
      write(*,*) 

c      call dscal(nion3*nion3,1.0d3,hessout,1)
c
c      dcount = nion3*nion3
c      call dsyev('N','U',nion3,hessout,nion3,vc,w1,dcount,ierr)
c
c      call dscal(nion3,1.0d-3,vc,1)
c      
c      write(*,*)"ierr=",ierr
c      write(*,*)
c      do i=1,nion3
c         write(*,*) "w(",i,")=",vc(i,1),dsqrt(dabs(vc(i,1)))*wave
c
c      end do


      return
      end

*     *************************************************
*     *                                               *
*     *                 vib_gendynamicmatrix          *
*     *                                               *
*     *************************************************

      subroutine vib_gendynamicmatrix(unita,q,nion,rion,
     >                                nion3,hess,dmat,wp,
     >                                work,lwork,rwork)
      implicit none
      real*8     unita(3,3),q(3)
      integer    nion
      real*8     rion(3,nion)
      integer    nion3
      real*8     hess(nion3,nion3)
      complex*16 dmat(nion3,nion3)
      real*8     wp(nion3)
      integer    lwork
      complex*16 work(lwork)
      real*8     rwork(3*nion3)

*     **** local variables ****
      integer ii,jj,i,j,xyz,rst,ierr
      integer i1,i2,i3,dcount
      real*8  dx,dy,dz,x,y,z,d2,d2min,qrmin,fac
      complex*16 cfac
      
      write(*,*) "generating dyanmical matrix"
      do jj=1,nion

*        **** diagonal term ****
         do xyz=1,3
         do rst=1,3
           i = 3*(jj-1)+xyz
           j = 3*(jj-1)+rst
           dmat(i,j) = hess(i,j)
         end do
         end do

         do ii=jj+1,nion
            dx = rion(1,jj)-rion(1,ii)
            dy = rion(1,jj)-rion(1,ii)
            dz = rion(1,jj)-rion(1,ii)

            d2min=9.99d12
            do i3=-1,1
            do i2=-1,1
            do i1=-1,1
               x = dx + unita(1,1)*i1 + unita(1,2)*i2 + unita(1,3)*i3
               y = dy + unita(2,1)*i1 + unita(2,2)*i2 + unita(2,3)*i3
               z = dz + unita(3,1)*i1 + unita(3,2)*i2 + unita(3,3)*i3
               d2 = x**2 + y**2 + z**2
               if (d2.lt.d2min) then
                  d2min = d2
                  qrmin = q(1)*x + q(2)*y + q(3)*z
               end if
            end do
            end do
            end do
            cfac = dcmplx(dcos(qrmin),dsin(qrmin))
      
            do xyz=1,3
            do rst=1,3
               i = 3*(ii-1)+xyz
               j = 3*(jj-1)+rst
               dmat(i,j) = cfac        *hess(i,j)
               dmat(j,i) = dconjg(cfac)*hess(j,i)
            end do
            end do

         end do !**ii**

      end do !**jj**

      call dscal(2*nion3*nion3,1.0d3,dmat,1)
      call zheev('N','U',nion3,dmat,nion3,wp,work,lwork,rwork,ierr)
      call dscal(nion3,1.0d-3,wp,1)

      return
      end

*     *************************************************
*     *                                               *
*     *                 vib_getlattice                *
*     *                                               *
*     *************************************************

      subroutine vib_getlattice(rtdb,omega,unita,unitg,lattice)
      implicit none
      integer rtdb
      real*8 omega
      real*8 unita(3,3)
      real*8 unitg(3,3)
      real*8 lattice(6)

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer geom
      real*8  twopi

      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('vib_getlattice: error creating geometry',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('vib_getlattice: error loading geometry',911, RTDB_ERR)

      if (.not. geom_amatrix_get(geom,unita))
     $  call errquit('driver: failed to get lattice',0,0)

      if (.not.geom_destroy(geom)) call errquit
     &    ('vib_getlattice: geom_destroy failed',911, GEOM_ERR)

      twopi = 8.0d0*datan(1.0d0)

*     -----------------------------------------
*     primitive vectors in the reciprocal space 
*     -----------------------------------------
      unitg(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      unitg(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      unitg(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      unitg(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      unitg(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      unitg(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      unitg(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      unitg(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      unitg(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)

*     ---------------------
*     volume of a unit cell
*     ---------------------
      omega = unita(1,1)*unitg(1,1)
     >      + unita(2,1)*unitg(2,1)
     >      + unita(3,1)*unitg(3,1)
      call dscal(9,twopi/omega,unitg,1)

      call xlattice_abc_abg(lattice(1),lattice(2),lattice(3),
     >                      lattice(4),lattice(5),lattice(6),unita)

      return
      end

*     *************************************************
*     *                                               *
*     *                 vib_isystype                  *
*     *                                               *
*     *************************************************

* returns the system type from geometry

      integer function vib_isystype(rtdb)
      implicit none
      integer rtdb

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"

      integer geom,isystype

      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('vib_isystype: error creating geometry',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('vib_isystype: error loading geometry',911, RTDB_ERR)

      if (.not. geom_systype_get(geom,isystype))
     $     call errquit('vib_isystype: systype?',0, GEOM_ERR)

      if (.not.geom_destroy(geom)) call errquit
     &    ('vib_isystype: geom_destroy failed',911, GEOM_ERR)

      vib_isystype = isystype
      return
      end

