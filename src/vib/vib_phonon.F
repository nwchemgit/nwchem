
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

c      double precision ams, wave
c      integer numans
c      COMMON /cvib_SETCON/ AMS(36),wave,numans      ! setup parameters
      real*8   autocm,autoTHz
      parameter (autocm=219474.6313705d0,autoTHz=6579.683920729d0)

*     **** local variables ****
      logical value,oprint,hprint
      integer dosgrid(3),isize,i1,i2,i3,i,k,ind,num_kvectors,icount,indx
      integer npoints
      real*8  omega,unita(3,3),unitg(3,3),lattice(6),emin,emax
      real*8  xx,yy,zz,xxx,yyy,zzz,ks(3),kq(3),kold(3),pathlength,dist
      complex*16 cfac
      character*50 zone_name,rtdb_name
      character*255 full_filename

      integer weight_dos(2),eigs_dos(2),dmat(2),wp(2),hessadjust(2)
      integer ierr,lwork,work(2),rwork(2),w2(1),kvs(2)


      oprint = util_print('phonon',print_medium)
      hprint = util_print('phonon',print_debug)

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
        write(luout,'(A,I5)') '    number of atoms =',nion 
        write(luout,'(A,I5)') '    number of bands =',nion3
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
         dosgrid(1) = 11
         dosgrid(2) = 11
         dosgrid(3) = 11
      end if
      if (oprint) write(luout,1250) dosgrid

*     **** allocate memory calculation ****
      isize = dosgrid(1)*dosgrid(2)*dosgrid(3)*nion3
      value = MA_alloc_get(mt_dbl,isize,
     >                     'eigs_dos',eigs_dos(2),eigs_dos(1))

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



*     ************************************************************
*     ******************* Density of States **********************
*     ************************************************************

      if (oprint) write(luout,*) 
      if (oprint) write(luout,*) "     Calculating Density of States"
      icount = 1
      xxx = 1.0d0/(1.0d0*dosgrid(1))
      yyy = 1.0d0/(1.0d0*dosgrid(2))
      zzz = 1.0d0/(1.0d0*dosgrid(3))
      do i3=0,(dosgrid(3)-1)
      do i2=0,(dosgrid(2)-1)
      do i1=0,(dosgrid(1)-1)
          ks(1) = i1*xxx
          ks(2) = i2*yyy
          ks(3) = i3*zzz
          kq(1) = unitg(1,1)*ks(1)+unitg(1,2)*ks(2)+unitg(1,3)*ks(3)
          kq(2) = unitg(2,1)*ks(1)+unitg(2,2)*ks(2)+unitg(2,3)*ks(3)
          kq(3) = unitg(2,1)*ks(1)+unitg(2,2)*ks(2)+unitg(2,3)*ks(3)

          call vib_gendynamicmatrix(unita,kq,nion,rion,
     >                          nion3,dbl_mb(hessadjust(1)),
     >                          dcpl_mb(dmat(1)),dbl_mb(wp(1)),
     >                          dcpl_mb(work(1)),lwork,dbl_mb(rwork(1)))

          do i=1,nion3
             indx = eigs_dos(1)
     >            + (icount-1)
     >            + (i-1)*dosgrid(1)*dosgrid(2)*dosgrid(3)
             dbl_mb(indx) = dsqrt(dabs(dbl_mb(wp(1)+i-1)))
          end do

          if (hprint) then
             write(luout,*) 
             write(luout,*) icount," out of ",
     >                      dosgrid(1)*dosgrid(2)*dosgrid(3)
             write(luout,*) "kvec= ",kq
             write(luout,*) "ksvec=",ks
             do i=1,nion3
                write(luout,'(5x,A3,I4,A2,E10.3,F10.3)') "wp(",i,")=",
     >             dsqrt(dabs(dbl_mb(wp(1)+i-1))),
     >             dsqrt(dabs(dbl_mb(wp(1)+i-1)))*autocm
             end do
          end if
          icount = icount + 1

      end do
      end do
      end do


*     ************************************************************
*     ********************   Dispersion **************************
*     ************************************************************

*     ***** writeout dispersion curves to nwchem_name.dispersion ****
      if (.not.rtdb_cget(rtdb,'vib:zone_name',1,zone_name))
     >   zone_name = 'vib_zone_default'

      ind = index(zone_name,' ') -1
      rtdb_name = zone_name(1:ind)//':number_kvectors'
      if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >   num_kvectors = 0
      if (num_kvectors.gt.0) then

         if (.not.MA_push_get(mt_dbl,(4*num_kvectors),
     >                       'kvs',kvs(2),kvs(1)))
     >      call errquit('vib_phonon: out of stack memory', 0,
     >       MA_ERR)

         rtdb_name = zone_name(1:ind)//':kvectors'
         if (.not.rtdb_get(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1))))
     >    call errquit('vib_phonon: failed to get kvs',0,
     &                 RTDB_ERR)

         call util_file_name('dispersion',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         ind = index(full_filename,' ') -1
         if (oprint) 
     >      write(luout,*) "     Generating Dispersion Curves: ",
     >                     full_filename(1:ind)
         open(unit=57,file=full_filename,form='formatted')
         call util_file_name('cm.dispersion',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         ind = index(full_filename,' ') -1
         if (oprint) 
     >      write(luout,*) "     Generating Dispersion Curves: ",
     >                     full_filename(1:ind)
         open(unit=58,file=full_filename,form='formatted')
         call util_file_name('THz.dispersion',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         ind = index(full_filename,' ') -1
         if (oprint) 
     >      write(luout,*) "     Generating Dispersion Curves: ",
     >                     full_filename(1:ind)
         open(unit=59,file=full_filename,form='formatted')

         do i=1,num_kvectors
           ks(1) = dbl_mb(kvs(1)+4*(i-1))
           ks(2) = dbl_mb(kvs(1)+4*(i-1)+1)
           ks(3) = dbl_mb(kvs(1)+4*(i-1)+2)
           kq(1) = unitg(1,1)*ks(1)+unitg(1,2)*ks(2)+unitg(1,3)*ks(3)
           kq(2) = unitg(2,1)*ks(1)+unitg(2,2)*ks(2)+unitg(2,3)*ks(3)
           kq(3) = unitg(2,1)*ks(1)+unitg(2,2)*ks(2)+unitg(2,3)*ks(3)

           call vib_gendynamicmatrix(unita,kq,nion,rion,
     >                          nion3,dbl_mb(hessadjust(1)),
     >                          dcpl_mb(dmat(1)),dbl_mb(wp(1)),
     >                          dcpl_mb(work(1)),lwork,dbl_mb(rwork(1)))

            if (i.eq.1) then
               pathlength = 0.0
            else
               dist=dsqrt((kq(1)-kold(1))**2
     >                   +(kq(2)-kold(2))**2
     >                   +(kq(3)-kold(3))**2)
               pathlength = pathlength + dist
            end if
            kold(1) = kq(1)
            kold(2) = kq(2)
            kold(3) = kq(3)
            write(57,'(1000E14.6)') pathlength,
     >          (dsqrt(dabs(dbl_mb(wp(1)+i1-1))),i1=1,nion3)
            write(58,'(1000E14.6)') pathlength,
     >          (dsqrt(dabs(dbl_mb(wp(1)+i1-1)))*autocm,i1=1,nion3)
            write(59,'(1000E14.6)') pathlength,
     >          (dsqrt(dabs(dbl_mb(wp(1)+i1-1)))*autoTHz,i1=1,nion3)

            if (hprint) then
               write(luout,*) 
               write(luout,*) i," out of ",num_kvectors
               write(luout,*) "kvec= ",kq
               write(luout,*) "ksvec=",ks
               do i1=1,nion3
                 write(luout,'(5x,A3,I4,A2,E10.3,F10.3)') "wp(",i1,")=",
     >               dsqrt(dabs(dbl_mb(wp(1)+i1-1))),
     >               dsqrt(dabs(dbl_mb(wp(1)+i1-1)))*autocm
               end do
            end if
         end do
         close(57)
         close(58)
         close(59)

         if (.not.MA_pop_stack(kvs(2)))
     >   call errquit('vib_phonon: failed to free stack',1,MA_ERR)

      end  if

*     **** deallocate stack memory ****
      value =           MA_pop_stack(work(2))
      value = value.and.MA_pop_stack(rwork(2))
      value = value.and.MA_pop_stack(wp(2))
      value = value.and.MA_pop_stack(hessadjust(2))
      value = value.and.MA_pop_stack(dmat(2))
      if (.not.value)
     >   call errquit('vib_phonon:popping stack',1,MA_ERR)


*     ************************************************************
*     ******************** DOS plotting **************************
*     ************************************************************

*     **** allocate memory for DOS calculation ****
      isize = dosgrid(1)*dosgrid(2)*dosgrid(3)*nion3
c      if (.not.MA_push_get(mt_dbl,isize,
c     >                    'weight_dos',weight_dos(2),weight_dos(1)))
c     >   call errquit('vib_phonon: out of stack memory',2,MA_ERR)


      if (.not.rtdb_get(rtdb,'vib:dos:npoints',mt_int,1,npoints)) then
        npoints = 1000
      end if

      if (.not.rtdb_get(rtdb,'vib:dos:emin',mt_dbl,1,emin)) emin=0.0d0
      if (.not.rtdb_get(rtdb,'vib:dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do i=1,(nion3)*dosgrid(1)*dosgrid(2)*dosgrid(3)
           if (dbl_mb(eigs_dos(1)+i-1).gt.emax)
     >       emax = dbl_mb(eigs_dos(1)+i-1)
         end do
         emax = emax + 1.0d-3
      end if

      call util_file_name('vib.dos',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      ind = index(full_filename,' ') -1
      if (oprint) write(luout,*) 
      if (oprint) write(luout,*) "     Generating Density of States: ",
     >                           full_filename(1:ind)
      open(unit=58,file=full_filename,form='formatted')
      call vib_dos_generate(58,dosgrid(1),dosgrid(2),dosgrid(3),
     >                      dbl_mb(eigs_dos(1)),nion3,
     >                      1.0d0,1.0d0,
     >                      npoints,emin,emax,unitg)
      close(58)

      call util_file_name('vib.cm.dos',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      ind = index(full_filename,' ') -1
      if (oprint) write(luout,*) "     Generating Density of States: ",
     >                           full_filename(1:ind)
      open(unit=58,file=full_filename,form='formatted')
      call vib_dos_generate(58,dosgrid(1),dosgrid(2),dosgrid(3),
     >                      dbl_mb(eigs_dos(1)),nion3,
     >                      autocm,1.0d0,
     >                      npoints,emin,emax,unitg)
      close(58)

      call util_file_name('vib.THz.dos',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      ind = index(full_filename,' ') -1
      if (oprint) write(luout,*) "     Generating Density of States: ",
     >                           full_filename(1:ind)
      open(unit=58,file=full_filename,form='formatted')
      call vib_dos_generate(58,dosgrid(1),dosgrid(2),dosgrid(3),
     >                      dbl_mb(eigs_dos(1)),nion3,
     >                      autoTHz,1.0d0,
     >                      npoints,emin,emax,unitg)
      close(58)



*     ************************************************************
*     ************************* THERMO ***************************
*     ************************************************************
      if (oprint) write(luout,*) 
      if (oprint) write(luout,*) "     Generating Thermo"


*     **** deallocate memory for DOS and THERMO calculations ****
c      value =           MA_pop_stack(weight_dos(2))
      value = MA_free_heap(eigs_dos(2))
      if (.not.value)
     >   call errquit('vib_phonon:popping stack',0,MA_ERR)

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

c      double precision ams, wave
c      integer numans
c      COMMON /cvib_SETCON/ AMS(36),wave,numans      ! setup parameters


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
            fac = fac/1822.89d0
            do xyz=1,3
            do rst=1,3
               i = 3*(ii-1)+xyz
               j = 3*(jj-1)+rst
               hessout(i,j) = hessout(i,j)*fac
            end do
            end do
         end do
      end do

c      !*** call Gainsville matrix output ***
c      write(*,*) 
c      write(*,*) "hess="
c      call output(hess,1,nion3,1,nion3,nion3,nion3,1)
c      write(*,*) 
c      write(*,*) "hessout="
c      call output(hessout,1,nion3,1,nion3,nion3,nion3,1)
c      write(*,*) 

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

      call dscal(2*nion3*nion3,1.0d6,dmat,1)
      call zheev('N','U',nion3,dmat,nion3,wp,work,lwork,rwork,ierr)
      call dscal(nion3,1.0d-6,wp,1)

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




**************** Definitions of Cubes and Tetrahedrons *****************
*                                                                      *
*                                                                      *
*     (011)------------(111)                    ( 3 )------------( 7 ) *
*       +                +                        +                +   *
*      /.               /|                       /.               /|   *
*     / .              / |                      / .              / |   *
*    /  .             /  |                     /  .             /  |   *
*   /   .            /   |                    /   .            /   |   *
* (001)------------(101) |      <====>      ( 1 )------------( 5 ) |   *
*   |   .            |   |                    |   .            |   |   *
*   | (010)..........|.(110)                  | ( 2 )..........|.( 6 ) *
*   |   .            |   /                    |   .            |   /   *
*   |  .             |  /                     |  .             |  /    *
*   | .              | /                      | .              | /     *
*   |.               |/                       |.               |/      *
*   +                +                        +                +       *
* (000)------------(100)                    ( 0 )------------( 4 )     *
*                                                                      *
*                                                                      *
* Algorithm to find diagaonals                                         *
*                                                                      *
*  Given a cube vertice d1                                             *
*  then d2 = d1^(111) = d1^7                                           *
*                                                                      *
*   Where the cOR bit operator "^" is defined as follows:              *
*      0^0 = 0                                                         *
*      1^1 = 0                                                         *
*      1^0 = 1                                                         *
*      0^1 = 1                                                         *
*                                                                      *
* The four possible cube diagonals are                                 *
*     (000) --- (111)                              (0, 7)              *
*     (001) --- (110)           <====>  2-tuple    (1, 6)              *
*     (010) --- (101)                   rep.       (2, 5)              *
*     (011) --- (100)                              (3, 4)              *
*                                                                      *
* Given a 2-tuple (d1,d2) that defines the diagonal of the cube,       *
* six tetrahedrons are defined, e.g.                                   *
*                                                                      *
*                      (111)                                           *
*                     .  / .                                           *
*                   .   /  .                                           *
*                 .    /   .                                           *
*                .    /   .                                            *
*              .     /    .                                            *
*             .    (101)  .                                            *
*           .     /  |   .        <====> 4-tuple (0, 7, 4, 5)          *
*          .    /    |   .               rep.                          *
*        .    /      |   .                                             *
*       .   /        |  .                                              *
*     .   /          |  .                                              *
*    .  /            |  .                                              *
*  .  /              | .                                               *
* (000)------------(100)                                               *
*                                                                      *
*                                                                      *
* Algorithm to find the six tetradedrons                               *
*                                                                      *
*  Given the diagonals vertices d1 and d2 such that d2=d1^7, the six   * 
*  tetradedrons (six 4-tuples) can be found using the following        *
*  algorithm:                                                          *
*                                                                      *
*   shift(0) = (001) = 1                                               *
*   shift(1) = (010) = 2                                               *
*   shift(2) = (100) = 4                                               *
*   tcount = 0                                                         *
*   For i=0,2                                                          *
*   For j=0,2                                                          *
*     c1 = d1^shift(i)                                                 *
*     c2 = c1^shift(j)                                                 *
*     If (c1 != d1) and (c1 != d2) and (c2!=d1) and (c2!=d2) Then      *
*       tetra(tcount) = (d1,d2,c1,c2)                                  *
*       tcount = tcount + 1                                            *
*     End If                                                           *
*   End For                                                            *
*   End For                                                            *
*                                                                      *
**************** Definitions of Cubes and Tetrahedrons *****************


*     *********************************************
*     *                                           *
*     *             vib_dos_generate              *
*     *                                           *
*     *********************************************

      subroutine vib_dos_generate(unit,idx,idy,idz,eigs,neigs,
     >                            ss,sgn,npoints,emin,emax,l_unitg)
      implicit none
      integer unit
      integer idx,idy,idz
      real*8 eigs(idx,idy,idz,*)
      integer neigs
      real*8  ss,sgn
      integer npoints
      real*8 emin,emax
      real*8 l_unitg(3,3)

*     **** local variables ****
      integer dosgrid(3)
      integer i,j,k,ii,jj,kk,ncubes,ntetra,count
      integer ishft,jshft,kshft
      integer k1_d(4),k2_d(4),k3_d(4),k1_dd(4),k2_dd(4),k3_dd(4)
      integer id,d1(4),d2(4)
      integer itetra(4,6)
      real*8  VT,VG
      real*8  B(3,3),unitg(3,3),e,ecube(8),f,g,de
      real*8  k1,k2,k3,kx,ky,kz,kkx,kky,kkz,r,rmax

*     **** external functions ****
      real*8   vib_Dstates_Cube,vib_Nstates_Cube
      external vib_Dstates_Cube,vib_Nstates_Cube

       dosgrid(1) = idx
       dosgrid(2) = idy
       dosgrid(3) = idz

c      write(unit,*) "dosgrid:",dosgrid
c      write(unit,*) "neigs:     ",neigs
c      write(unit,*) "sign:      ",sign
c      write(unit,*) "npoints:   ",npoints
c      write(unit,*) "emin:      ", emin
c      write(unit,*) "emax:      ", emax

      do j=1,3
      do i=1,3
        B(i,j) = l_unitg(i,j)
      end do  
      end do  

*     **** volume of reciprocal unit cell, VG ****
      unitg(1,1) = B(2,2)*B(3,3) - B(3,2)*B(2,3)
      unitg(2,1) = B(3,2)*B(1,3) - B(1,2)*B(3,3)
      unitg(3,1) = B(1,2)*B(2,3) - B(2,2)*B(1,3)

      unitg(1,2) = B(2,3)*B(3,1) - B(3,3)*B(2,1)
      unitg(2,2) = B(3,3)*B(1,1) - B(1,3)*B(3,1)
      unitg(3,2) = B(1,3)*B(2,1) - B(2,3)*B(1,1)

      unitg(1,3) = B(2,1)*B(3,2) - B(3,1)*B(2,2)
      unitg(2,3) = B(3,1)*B(1,2) - B(1,1)*B(3,2)
      unitg(3,3) = B(1,1)*B(2,2) - B(2,1)*B(1,2)
      VG = B(1,1)*unitg(1,1)
     >   + B(2,1)*unitg(2,1)
     >   + B(3,1)*unitg(3,1)
      
      ncubes = dosgrid(1)*dosgrid(2)*dosgrid(3)
      ntetra = ncubes*6
      VT = VG/dble(ntetra)

c      write(unit,*) "VG:     ",VG
c      write(unit,*) "number of cubes:",ncubes
c      write(unit,*) "number of tetra:",ntetra
c      write(unit,*) "VT:     ",VT
c
c      count = 0
c      do k=0,dosgrid(3)-1
c      do j=0,dosgrid(2)-1
c      do i=0,dosgrid(1)-1
c         count = count + 1
c         k1 = (dble(i)/dble(dosgrid(1)))
c         k2 = (dble(j)/dble(dosgrid(2)))
c         k3 = (dble(k)/dble(dosgrid(3)))
c         kx = k1*B(1,1) + k2*B(1,2) + k3*B(1,3)
c         ky = k1*B(2,1) + k2*B(2,2) + k3*B(2,3)
c         kz = k1*B(3,1) + k2*B(3,2) + k3*B(3,3)
c         write(unit,*) i,j,k
c         write(unit,3508) count,k1,k2,k3,kx,ky,kz
c         write(unit,*)
c      end do
c      end do
c      end do

*     ********************************
*     **** find shortest diagonal ****
*     ********************************

*     **** (000) ---- (111) ****
      k1_d(1) = 0
      k2_d(1) = 0
      k3_d(1) = 0
      k1_dd(1) = 1
      k2_dd(1) = 1
      k3_dd(1) = 1
      d1(1) = 0
      d2(1) = 7

*     **** (001) ---- (110) ****
      k1_d(2) = 1
      k2_d(2) = 0
      k3_d(2) = 0
      k1_dd(2) = 0
      k2_dd(2) = 1
      k3_dd(2) = 1
      d1(2) = 1
      d2(2) = 6

*     **** (010) ---- (101) ****
      k1_d(3) = 0
      k2_d(3) = 1
      k3_d(3) = 0
      k1_dd(3) = 1
      k2_dd(3) = 0
      k3_dd(3) = 1
      d1(3) = 2
      d2(3) = 5

*     **** (011) ---- (100) ****
      k1_d(4) = 1
      k2_d(4) = 1
      k3_d(4) = 0
      k1_dd(4) = 0
      k2_dd(4) = 0
      k3_dd(4) = 1
      d1(4) = 3
      d2(4) = 4

      id = 1
      rmax = 9.99d9
      do i=1,4
         kx = k1_d(i)*B(1,1) + k2_d(i)*B(1,2) + k3_d(i)*B(1,3)
         ky = k1_d(i)*B(2,1) + k2_d(i)*B(2,2) + k3_d(i)*B(2,3)
         kz = k1_d(i)*B(3,1) + k2_d(i)*B(3,2) + k3_d(i)*B(3,3)

         kkx = k1_dd(i)*B(1,1) + k2_dd(i)*B(1,2) + k3_dd(i)*B(1,3)
         kky = k1_dd(i)*B(2,1) + k2_dd(i)*B(2,2) + k3_dd(i)*B(2,3)
         kkz = k1_dd(i)*B(3,1) + k2_dd(i)*B(3,2) + k3_dd(i)*B(3,3)
         r = (kx-kkx)**2 + (ky-kky)**2 + (kz-kkz)**2
         !write(unit,*) "diagonal distance:",i,r,rmax
         if (r.lt.rmax) then
           rmax = r
           id = i
         end if
      end do

      !write(unit,*) "diagonal d1,d2 =",d1(id),d2(id)

*     **** define six tetradrons - clunky but don't know defn of cOR in fortran ****
      if (id.eq.1) then
        itetra(1,1) = 0   +1
        itetra(2,1) = 7   +1
        itetra(3,1) = 1   +1
        itetra(4,1) = 3   +1

        itetra(1,2) = 0   +1
        itetra(2,2) = 7   +1
        itetra(3,2) = 1   +1
        itetra(4,2) = 5   +1

        itetra(1,3) = 0   +1
        itetra(2,3) = 7   +1
        itetra(3,3) = 2   +1
        itetra(4,3) = 3   +1

        itetra(1,4) = 0   +1
        itetra(2,4) = 7   +1
        itetra(3,4) = 2   +1
        itetra(4,4) = 6   +1

        itetra(1,5) = 0   +1
        itetra(2,5) = 7   +1
        itetra(3,5) = 4   +1
        itetra(4,5) = 5   +1

        itetra(1,6) = 0   +1
        itetra(2,6) = 7   +1
        itetra(3,6) = 4   +1
        itetra(4,6) = 6   +1
      else if (id.eq.2) then
        itetra(1,1) = 1   +1
        itetra(2,1) = 6   +1
        itetra(3,1) = 0   +1
        itetra(4,1) = 2   +1

        itetra(1,2) = 1   +1
        itetra(2,2) = 6   +1
        itetra(3,2) = 0   +1
        itetra(4,2) = 4   +1

        itetra(1,3) = 1   +1
        itetra(2,3) = 6   +1
        itetra(3,3) = 3   +1
        itetra(4,3) = 2   +1

        itetra(1,4) = 1   +1
        itetra(2,4) = 6   +1
        itetra(3,4) = 3   +1
        itetra(4,4) = 7   +1

        itetra(1,5) = 1   +1
        itetra(2,5) = 6   +1
        itetra(3,5) = 5   +1
        itetra(4,5) = 4   +1

        itetra(1,6) = 1   +1
        itetra(2,6) = 6   +1
        itetra(3,6) = 5   +1
        itetra(4,6) = 7   +1
      else if (id.eq.3) then
        itetra(1,1) = 2   +1
        itetra(2,1) = 5   +1
        itetra(3,1) = 3   +1
        itetra(4,1) = 1   +1

        itetra(1,2) = 2   +1
        itetra(2,2) = 5   +1
        itetra(3,2) = 3   +1
        itetra(4,2) = 7   +1

        itetra(1,3) = 2   +1
        itetra(2,3) = 5   +1
        itetra(3,3) = 0   +1
        itetra(4,3) = 1   +1

        itetra(1,4) = 2   +1
        itetra(2,4) = 5   +1
        itetra(3,4) = 0   +1
        itetra(4,4) = 4   +1

        itetra(1,5) = 2   +1
        itetra(2,5) = 5   +1
        itetra(3,5) = 6   +1
        itetra(4,5) = 7   +1

        itetra(1,6) = 2   +1
        itetra(2,6) = 5   +1
        itetra(3,6) = 6   +1
        itetra(4,6) = 4   +1
      else if (id.eq.4) then
        itetra(1,1) = 3   +1
        itetra(2,1) = 4   +1
        itetra(3,1) = 2    +1
        itetra(4,1) = 0   +1

        itetra(1,2) = 3   +1
        itetra(2,2) = 4   +1
        itetra(3,2) = 2   +1
        itetra(4,2) = 6   +1

        itetra(1,3) = 3   +1
        itetra(2,3) = 4   +1
        itetra(3,3) = 1   +1
        itetra(4,3) = 0   +1

        itetra(1,4) = 3   +1
        itetra(2,4) = 4   +1
        itetra(3,4) = 1   +1
        itetra(4,4) = 5   +1

        itetra(1,5) = 3   +1
        itetra(2,5) = 4   +1
        itetra(3,5) = 7   +1
        itetra(4,5) = 6   +1

        itetra(1,6) = 3   +1
        itetra(2,6) = 4   +1
        itetra(3,6) = 7   +1
        itetra(4,6) = 5   +1
      end if


c      do i=1,6
c        write(unit,*) id,"tetra :",i,"(",(itetra(j,i),j=1,4),")"
c      end do

      de = (emax-emin)/dble(npoints-1)
      do k=1,npoints
        e = emin + (k-1)*de

        f = 0.0d0
        g = 0.0d0
        do kk=1,dosgrid(3)
        do jj=1,dosgrid(2)
        do ii=1,dosgrid(1)
          ishft = ii+1
          jshft = jj+1
          kshft = kk+1
          if (ishft.gt.dosgrid(1)) ishft=1
          if (jshft.gt.dosgrid(2)) jshft=1
          if (kshft.gt.dosgrid(3)) kshft=1
          do i=1,neigs
            ecube(1) = eigs(ii,       jj,    kk, i)  ! (000)
            ecube(2) = eigs(ishft,    jj,    kk, i)  ! (001)
            ecube(3) = eigs(ii,    jshft,    kk, i)  ! (010)
            ecube(4) = eigs(ishft, jshft,    kk, i)  ! (011)
            ecube(5) = eigs(ii,       jj, kshft, i)  ! (100)
            ecube(6) = eigs(ishft,    jj, kshft, i)  ! (101)
            ecube(7) = eigs(   ii, jshft, kshft, i)  ! (110)
            ecube(8) = eigs(ishft, jshft, kshft, i)  ! (111)

           ! write(*,*) "ecube=",e,ii,jj,kk,i,ecube
           
            f = f + vib_Dstates_Cube(e,itetra,ecube)
            g = g + vib_Nstates_Cube(e,itetra,ecube)

           ! write(*,*) "f,g  =",e,ii,jj,kk,i,f,g,VT/VG
          end do
        end do
        end do
        end do
        f = f*(VT/VG)
        g = g*(VT/VG)

        write(unit,1310) ss*e,f*sgn,g*sgn
      end do


      return
 1310 FORMAT(3E15.5)
 3508 FORMAT(/' Brillouin zone point: ',i5,
     >       /'    k     =<',3f8.3,'> . <b1,b2,b3> ',
     >       /'          =<',3f8.3,'>')
      end


      real*8 function vib_Dstates_Cube(e,itetra,ecube)
      implicit none
      real*8  e
      integer itetra(4,6)
      real*8  ecube(8)

*     **** local variables ****
      integer i,j,k
      real*8 ds,etetra(4),swap

      real*8   vib_Dstates_Tetra
      external vib_Dstates_Tetra

*     **** sum over 6 tetrahedrons ****
      ds = 0.0d0
      do k=1,6
        etetra(1) = ecube(itetra(1,k))
        etetra(2) = ecube(itetra(2,k))
        etetra(3) = ecube(itetra(3,k))
        etetra(4) = ecube(itetra(4,k))

*       **** bubble sort ****
        do j=1,3
        do i=j+1,4
          if (etetra(j).gt.etetra(i)) then
            swap      = etetra(i)
            etetra(i) = etetra(j)
            etetra(j) = swap
          end if
        end do
        end do
        ds = ds + vib_Dstates_Tetra(e,etetra)
      end do

      vib_Dstates_cube = ds
      return
      end


      real*8 function vib_Dstates_Tetra(e,ee)
      implicit none
      real*8 e
      real*8 ee(4)

*     **** local variables ****
      real*8 ds
      real*8 e1,e2,e4
      real*8 e21,e31,e41,e32,e42,e43

      if ((ee(1).le.e).and.(e.lt.ee(2))) then
        e1 = e-ee(1)
        e21 = ee(2) - ee(1)
        e31 = ee(3) - ee(1)
        e41 = ee(4) - ee(1)
        ds = 3.0d0*e1*e1/(e21*e31*e41)
      else if ((ee(2).le.e).and.(e.lt.ee(3))) then
        e2 = e-ee(2) 
        e21 = ee(2) - ee(1) 
        e31 = ee(3) - ee(1) 
        e41 = ee(4) - ee(1) 
        e32 = ee(3) - ee(2) 
        e42 = ee(4) - ee(2) 
        ds = (3.0d0*e21+6.0d0*e2-3.0d0*(e31+e42)*e2*e2/(e32*e42))
     >       /(e31*e41)
      else if ((ee(3).le.e).and.(e.lt.ee(4))) then
        e4 = ee(4)-e 
        e41 = ee(4) - ee(1) 
        e42 = ee(4) - ee(2) 
        e43 = ee(4) - ee(3) 
        ds = 3.0d0*e4*e4/(e41*e42*e43)
      else
        ds = 0.0d0
      end if


      vib_Dstates_Tetra = ds
      return
      end
      

      real*8 function vib_Nstates_Cube(e,itetra,ecube)
      implicit none
      real*8  e
      integer itetra(4,6)
      real*8  ecube(8)

*     **** local variables ****
      integer i,j,k
      real*8 ds,etetra(4),swap

      real*8   vib_Nstates_Tetra
      external vib_Nstates_Tetra

*     **** sum over 6 tetrahedrons ****
      ds = 0.0d0
      do k=1,6
        etetra(1) = ecube(itetra(1,k))
        etetra(2) = ecube(itetra(2,k))
        etetra(3) = ecube(itetra(3,k))
        etetra(4) = ecube(itetra(4,k))

*       **** bubble sort ****
        do j=1,3
        do i=j+1,4
          if (etetra(j).gt.etetra(i)) then
            swap      = etetra(i)
            etetra(i) = etetra(j)
            etetra(j) = swap
          end if
        end do
        end do
        ds = ds + vib_Nstates_Tetra(e,etetra)
      end do

      vib_Nstates_cube = ds
      return
      end


      real*8 function vib_Nstates_Tetra(e,ee)
      implicit none
      real*8 e
      real*8 ee(4)

*     **** local variables ****
      real*8 ds
      real*8 e1,e2,e4
      real*8 e21,e31,e41,e32,e42,e43

      if ((ee(1).le.e).and.(e.lt.ee(2))) then
        e1 = e-ee(1)
        e21 = ee(2) - ee(1)
        e31 = ee(3) - ee(1)
        e41 = ee(4) - ee(1)
        ds = e1*e1*e1/(e21*e31*e41)
      else if ((ee(2).le.e).and.(e.lt.ee(3))) then
        e2 = e-ee(2) 
        e21 = ee(2) - ee(1) 
        e31 = ee(3) - ee(1) 
        e41 = ee(4) - ee(1) 
        e32 = ee(3) - ee(2) 
        e42 = ee(4) - ee(2) 
        ds = (e21*e21 
     >        + 3.0d0*e21*e2
     >        + 3.0d0*e2*e2
     >        - (e31+e42)*e2*e2*e2/(e32*e42))
     >       /(e31*e41)
      else if ((ee(3).le.e).and.(e.lt.ee(4))) then
        e4 = ee(4)-e 
        e41 = ee(4) - ee(1) 
        e42 = ee(4) - ee(2) 
        e43 = ee(4) - ee(3) 
        ds = 1.0d0 - e4*e4*e4/(e41*e42*e43)
      else if (e.ge.ee(4)) then
        ds = 1.0d0
      else
        ds = 0.0d0
      end if


      vib_Nstates_Tetra = ds
      return
      end

