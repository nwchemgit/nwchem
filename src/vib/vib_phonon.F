
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
      logical value,oprint,hprint,animate
      integer dosgrid(3),isize,i1,i2,i3,i,k,ind,num_kvectors,icount,indx
      integer npoints
      real*8  omega,unita(3,3),unitg(3,3),lattice(6),emin,emax,gmin,gmax
      real*8  xx,yy,zz,xxx,yyy,zzz,ks(3),kq(3),kold(3),pathlength,dist
      real*8  animate_step
      complex*16 cfac
      character*50 zone_name,rtdb_name
      character*255 full_filename

      integer efn_dos(2),eigs_dos(2),dmat(2),wp(2),hessadjust(2)
      integer ierr,lwork,work(2),rwork(2),w2(2),kvs(2)


      oprint = util_print('phonon',print_medium)
      hprint = util_print('phonon',print_debug)

      if (oprint) then
         call util_print_centered(luout,
     >        'NWChem Phonon Spectra Analysis ',
     >        40,.true.)
         write(luout,*)
      end if

      if (.not.rtdb_get(rtdb,'vib:animate',mt_log,1,animate))
     >    animate=.false.

      if (animate) then
          if (.not.rtdb_get(rtdb,'vib:animate:step_size',
     >                      mt_dbl,1,animate_step)) 
     >       animate_step = 0.15d0
      end if


*     **** get lattice vectors, reciprocal lattice vectors ****
      call vib_getlattice(rtdb,omega,unita,unitg,lattice)

      if (oprint) then
        write(luout,*)
        write(luout,'(A,I5)') '      number of atoms =',nion 
        write(luout,'(A,I5)') '      number of bands =',nion3
        write(luout,*)
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
      if (oprint) write(luout,*) 
      if (oprint) write(luout,*) 
      if (oprint) write(luout,*) 
     > "     Calculating Vibrational Density of States"
      if (oprint) write(luout,*) 
     > "     -----------------------------------------"
      if (oprint) write(luout,1250) dosgrid
      if (oprint) write(luout,*) 

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
          kq(3) = unitg(3,1)*ks(1)+unitg(3,2)*ks(2)+unitg(3,3)*ks(3)

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

         call util_file_name('vib.dispersion',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         ind = index(full_filename,' ') -1
         if (oprint) 
     >      write(luout,*) "     Phonon Dispersion Curves: ",
     >                     full_filename(1:ind)
         open(unit=57,file=full_filename,form='formatted')
         call util_file_name('vib.cm.dispersion',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         ind = index(full_filename,' ') -1
         if (oprint) 
     >      write(luout,*) "     Phonon Dispersion Curves: ",
     >                     full_filename(1:ind)
         open(unit=58,file=full_filename,form='formatted')
         call util_file_name('vib.THz.dispersion',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         ind = index(full_filename,' ') -1
         if (oprint) 
     >      write(luout,*) "     Phonon Dispersion Curves: ",
     >                     full_filename(1:ind)
         open(unit=59,file=full_filename,form='formatted')


         do i=1,num_kvectors
           ks(1) = dbl_mb(kvs(1)+4*(i-1))
           ks(2) = dbl_mb(kvs(1)+4*(i-1)+1)
           ks(3) = dbl_mb(kvs(1)+4*(i-1)+2)
           kq(1) = unitg(1,1)*ks(1)+unitg(1,2)*ks(2)+unitg(1,3)*ks(3)
           kq(2) = unitg(2,1)*ks(1)+unitg(2,2)*ks(2)+unitg(2,3)*ks(3)
           kq(3) = unitg(3,1)*ks(1)+unitg(3,2)*ks(2)+unitg(3,3)*ks(3)

           call vib_gendynamicmatrix(unita,kq,nion,rion,
     >                          nion3,dbl_mb(hessadjust(1)),
     >                          dcpl_mb(dmat(1)),dbl_mb(wp(1)),
     >                          dcpl_mb(work(1)),lwork,dbl_mb(rwork(1)))

            if (animate) then
               if (oprint) 
     >          write(luout,'(A,F6.3,A,3F12.6)') 
     >          "      Animating Mode - animate_step_size=",
     >          animate_step," -  kvec= ",kq

               call vib_writedmat(rtdb,i,animate_step,unita,kq,
     >                          nion,rion,
     >                          nion3,dcpl_mb(dmat(1)))
            end if

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

c     **** animate at kq=<0,0,0> if no dispersion specified ****
      if ((animate).and.(num_kvectors.lt.1)) then
         kq(1) = 0.0d0
         kq(2) = 0.0d0
         kq(3) = 0.0d0
         call vib_gendynamicmatrix(unita,kq,nion,rion,
     >                      nion3,dbl_mb(hessadjust(1)),
     >                      dcpl_mb(dmat(1)),dbl_mb(wp(1)),
     >                      dcpl_mb(work(1)),lwork,dbl_mb(rwork(1)))
         call vib_writedmat(rtdb,0,animate_step,unita,kq,
     >                      nion,rion,
     >                      nion3,dcpl_mb(dmat(1)))
         if (oprint) 
     >     write(luout,'(A,F6.3,A,3F12.6)') 
     >     "      Animating Mode - animate_step_size=",
     >     animate_step," -  kvec= ",kq
        end if



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
      if (.not.rtdb_get(rtdb,'vib:dos:npoints',mt_int,1,npoints)) then
        npoints = 1001
      end if

*     **** allocate memory for DOS calculation ****
      isize = 3*npoints
      if (.not.MA_push_get(mt_dbl,isize,
     >                    'efn_dos',efn_dos(2),efn_dos(1)))
     >   call errquit('vib_phonon: out of stack memory',2,MA_ERR)



      if (.not.rtdb_get(rtdb,'vib:dos:emin',mt_dbl,1,emin)) emin=0.0d0
      if (.not.rtdb_get(rtdb,'vib:dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do i=1,(nion3)*dosgrid(1)*dosgrid(2)*dosgrid(3)
           if (dbl_mb(eigs_dos(1)+i-1).gt.emax)
     >       emax = dbl_mb(eigs_dos(1)+i-1)
         end do
         emax = emax + 1.0d-3
      end if
      call vib_dos_generate(dosgrid(1),dosgrid(2),dosgrid(3),
     >                      dbl_mb(eigs_dos(1)),nion3,
     >                      npoints,emin,emax,unitg,dbl_mb(efn_dos(1)))

      call util_file_name('vib.vdos',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      ind = index(full_filename,' ') -1
      if (oprint) write(luout,*) 
      if (oprint) write(luout,*) 
     > "     Vibrational Density of States: ",
     >                           full_filename(1:ind)
      open(unit=58,file=full_filename,form='formatted')
      do k=1,npoints
         write(58,'(3E15.5)') dbl_mb(efn_dos(1)          +k-1),
     >                        dbl_mb(efn_dos(1)+npoints  +k-1),
     >                        dbl_mb(efn_dos(1)+2*npoints+k-1)
      end do 
      close(58)

      call util_file_name('vib.cm.vdos',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      ind = index(full_filename,' ') -1
      if (oprint) write(luout,*) 
     > "     Vibrational Density of States: ",
     >                           full_filename(1:ind)
      open(unit=58,file=full_filename,form='formatted')
      do k=1,npoints
         write(58,'(3E15.5)') dbl_mb(efn_dos(1)          +k-1)*autocm,
     >                        dbl_mb(efn_dos(1)+npoints  +k-1)/autocm,
     >                        dbl_mb(efn_dos(1)+2*npoints+k-1)
      end do 
      close(58)

      call util_file_name('vib.THz.vdos',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      ind = index(full_filename,' ') -1
      if (oprint) write(luout,*) 
     > "     Vibrational Density of States: ",
     >                           full_filename(1:ind)
      open(unit=58,file=full_filename,form='formatted')
      do k=1,npoints
         write(58,'(3E15.5)') dbl_mb(efn_dos(1)          +k-1)*autoTHz,
     >                        dbl_mb(efn_dos(1)+npoints  +k-1)/autoTHz,
     >                        dbl_mb(efn_dos(1)+2*npoints+k-1)
      end do 
      close(58)


*     ************************************************************
*     ************************* THERMO ***************************
*     ************************************************************
      if (oprint) 
     >   call vib_vdos_thermo(rtdb,nion,npoints,dbl_mb(efn_dos(1)))


*     **** deallocate memory for DOS and THERMO calculations ****
      value =           MA_pop_stack(efn_dos(2))
      value = value.and.MA_free_heap(eigs_dos(2))
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
*     *              vib_vdos_thermo                  *
*     *                                               *
*     *************************************************

      subroutine vib_vdos_thermo(rtdb,nion,npoints,efn)
      implicit none
      integer rtdb,nion,npoints
      real*8 efn(npoints,3)

#include "errquit.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"


*     ****  conversions and constants ****
      double precision AUKCAL,c,h,kgas,Rgas,autocm,autoTHz
      parameter (AUKCAL=627.5093314d0)
      parameter (c=2.998e+10,h=6.626e-27,kgas=1.3807e-16)   ! cgs units
      parameter (Rgas = 1.9863d0/1000.0d0/AUKCAL)           ! atomic units
      parameter (autocm=219474.6313705d0,autoTHz=6579.683920729d0)

*     **** local variables ****
      integer k,it,ntemp,h_actlist,l_actlist,ma_type
      real*8 temp,scalefreq,threshfreq,RT,RT2,RT1
      real*8 ezero,ezero1,ezero2
      real*8 ethermal,ethermal1,ethermal2
      real*8 fhelmholtz,fhelmholtz1,fhelmholtz2
      real*8 Svib,Svib1,Svib2
      real*8 Cv_vib,Cv_vib1,Cv_vib2
      real*8 sa,sa1,sa2,sb,sb1,sb2,de,g,xdum,hbarw
      real*8 nmin,nmax,gmin,gmax,emin,emax,xx,yy,zz


*     ************************************************************
*     *************** ascii plot of the density of states ********
*     ************************************************************
      write(luout,*)
      write(luout,*)
      write(luout,*)
      write(luout,*) 
     > '     g(w) and N(w) Used For Calculating Thermodynamic Functions'
      write(luout,*) 
     > '     ----------------------------------------------------------'
      emin =  9999999.99d0
      gmin =  9999999.99d0
      nmin =  9999999.99d0
      emax = -9999999.99d0
      gmax = -9999999.99d0
      nmax = -9999999.99d0
      do k=1,npoints
         xx = efn(k,1)
         yy = efn(k,2)
         zz = efn(k,3)
         if (xx.lt.emin) emin=xx
         if (yy.lt.gmin) gmin=yy
         if (zz.lt.nmin) nmin=zz
         if (xx.gt.emax) emax=xx
         if (yy.gt.gmax) gmax=yy
         if (zz.gt.nmax) nmax=zz
      end do
      write(luout,10012) 
      write(luout,10013) npoints
      write(luout,10014) emin*autocm
      write(luout,10015) emax*autocm
      write(luout,10016) gmin/autocm
      write(luout,10017) gmax/autocm
      write(luout,10018) nmin
      write(luout,10019) nmax
      write(luout,*)
      write(luout,*)

      !**** scale efn ****
      do k=1,npoints
         efn(k,1) = efn(k,1)*autocm
         efn(k,2) = efn(k,2)/autocm
      end do

      if (dabs(gmax-gmin).lt.1.0d-3)  gmax = gmin+1.0d-3
      call util_ascii_setwindow(emin*autocm,emax*autocm,
     >                          gmin/autocm,gmax/autocm)
      !write(luout,'(A,E10.3,A)') "gmax=",gmax/autocm," states/cm-1"
      call util_ascii_plotter("",luout,npoints,efn(1,1),efn(1,2),'*',
     >                     "g(w) - Vibrational Density of States",
     >                      "w (cm-1)","g(w) (states/cm-1)")
      !write(luout,'(A,E10.3,A)') "gmin=",gmin/autocm," states/cm-1"
      !write(luout,*)
      !write(luout,900)emin*autocm,0.5d0*(emax+emin)*autocm,emax*autocm
      write(luout,*)
      write(luout,*)
      call util_ascii_setwindow(emin*autocm,emax*autocm,nmin,nmax)
      !write(luout,'(A,F7.1,A)') "Nmax=",nmax," states"
      call util_ascii_plotter("",luout,npoints,efn(1,1),efn(1,3),'*',
     >                     "N(w) - Number of Vibrational States",
     >                     "w (cm-1)", "N(w) (states)")
      !write(luout,'(A,F7.1,A)') "Nmin=",nmin," states"
      !write(luout,*)
      !write(luout,900)emin*autocm,0.5d0*(emax+emin)*autocm,emax*autocm
      write(luout,*)
c  900 format('wmin=',f6.1,'cm-1',18x,f6.1,'cm-1',
c     >                           17x,'wmax=',f6.1,'cm-1')

      !**** unscale efn ****
      do k=1,npoints
         efn(k,1) = efn(k,1)/autocm
         efn(k,2) = efn(k,2)*autocm
      end do

*     ************************************************************
*     *************** get the other input from rtdb       ********
*     ************************************************************

*     *** define scalefreq ***
      if(.not. rtdb_get(rtdb,'vib:scalefreq',MT_DBL,1,scalefreq))
     &    scalefreq = 1.0d0

      if(.not. rtdb_get(rtdb,'vib:threshfreq',MT_DBL,1,threshfreq))
     &    threshfreq = 100.0d0/autocm

*     **** define temperatures
      if (.not. rtdb_ma_get(rtdb, 'vib:temperature', ma_type,
     >        ntemp, h_actlist)) then
         ntemp = 1

         if (.not.MA_alloc_get(mt_dbl,(ntemp),'temp',
     >                            h_actlist,l_actlist))
     >   call errquit(
     >       'vib_WRTFREQ: ma_get_index failed for actlist',911, MA_ERR)

         dbl_mb(l_actlist) = 298.15d0
      else
         if (.not.MA_get_index(h_actlist,l_actlist))
     >   call errquit(
     >       'vib_WRTFREQ: ma_get_index failed for actlist',911, MA_ERR)
      end if

*     ************************************************************
*     *************** Compute the Thermodynamic Functions ********
*     ************************************************************

      write(luout,*)
      write(luout,*)
      write(luout,*) 
     > '     Thermodynamic Functions from Vibrational Density of States'
      write(luout,*) 
     > '     ----------------------------------------------------------'
      do it=1,ntemp
         temp = dbl_mb(l_actlist+it-1)
         RT = Rgas*temp
         RT2 = 1.0d0/(2.0d0*RT)
         RT1 = 1.0d0/(1.0d0*RT)

*        **********************************************************************
*        **** Write out thermal Corrections to Energies and entropies      ****
*        **** frequencies are assumed to be in units of a.u.               ****
*        **** Formulas taken from Herzberg                                 ****
*        **********************************************************************
         ezero      = 0.0d0
         ezero1     = 0.0d0
         ezero2     = 0.0d0
         ethermal   = 0.0d0
         ethermal1  = 0.0d0
         ethermal2  = 0.0d0
         fhelmholtz = 0.0d0
         fhelmholtz1= 0.0d0
         fhelmholtz2= 0.0d0
         sa         = 0.0d0
         sa1        = 0.0d0
         sa2        = 0.0d0
         sb         = 0.0d0
         sb1        = 0.0d0
         sb2        = 0.0d0
         Cv_vib     = 0.0d0
         Cv_vib1    = 0.0d0
         Cv_vib2    = 0.0d0
         de = efn(2,1)-efn(1,1)
         do k=2,npoints-1
            !g = efn(k,2)
            hbarw = scalefreq*efn(k,1)
            g     = 0.5d0*(efn(k+1,3)-efn(k-1,3))/de  !*** central difference derivative ****
            if (temp.gt.0.0d0) then
               xdum = dexp(-hbarw*RT1)
            else
               xdum = 0.0d0
            end if
            xdum = xdum/(1.0d0-xdum)
            ezero      = ezero      + g*hbarw*(0.5d0)
            ethermal   = ethermal   + g*hbarw*(0.5d0+xdum)
            fhelmholtz = fhelmholtz + g*RT*dlog(2.0d0*dsinh(hbarw*RT2))

            sa     = sa     + g*hbarw*xdum
            sb     = sb     + g*dlog(1.0d0-dexp(-hbarw*RT1))
            Cv_vib = Cv_vib + g*dexp(hbarw*RT1)*(hbarw*RT1*xdum)**2
            if (hbarw.lt.threshfreq) then
               ezero1      = ezero1   + g*hbarw*(0.5d0)
               ethermal1   = ethermal1+ g*hbarw*(0.5d0+xdum)
               fhelmholtz1=fhelmholtz1+g*RT*dlog(2.0d0*dsinh(hbarw*RT2))
               sa1    = sa1    + g*hbarw*xdum
               sb1    = sb1    + g*dlog(1.0d0-dexp(-hbarw*RT1))
               Cv_vib1= Cv_vib1+ g*dexp(hbarw*RT1)*(hbarw*RT1*xdum)**2
            else
               ezero2      = ezero2   + g*hbarw*(0.5d0)
               ethermal2   = ethermal2+ g*hbarw*(0.5d0+xdum)
               fhelmholtz2=fhelmholtz2+g*RT*dlog(2.0d0*dsinh(hbarw*RT2))
               sa2    = sa2    + g*hbarw*xdum
               sb2    = sb2    + g*dlog(1.0d0-dexp(-hbarw*RT1))
               Cv_vib2= Cv_vib2+ g*dexp(hbarw*RT1)*(hbarw*RT1*xdum)**2
            end if
 
         end do
         ezero      = ezero*de
         ezero1     = ezero1*de
         ezero2     = ezero2*de
         ethermal   = ethermal*de
         ethermal1  = ethermal1*de
         ethermal2  = ethermal2*de
         fhelmholtz = fhelmholtz*de
         fhelmholtz1= fhelmholtz1*de
         fhelmholtz2= fhelmholtz2*de
         if (temp.gt.0.0d0) then
            Svib = sa*de/temp  - Rgas*sb*de
            Svib1= sa1*de/temp - Rgas*sb1*de
            Svib2= sa2*de/temp - Rgas*sb2*de
         else
            Svib = 0.0d0
            Svib1= 0.0d0
            Svib2= 0.0d0
         end if
         Cv_vib = Rgas*Cv_vib*de
         Cv_vib1= Rgas*Cv_vib1*de
         Cv_vib2= Rgas*Cv_vib2*de

         write(luout,10023) temp,RT
         write(luout,10022) scalefreq
         write(luout,10021) threshfreq*autocm

         write(luout,10026) fhelmholtz*AUKCAL,fhelmholtz
         write(luout,10027) 
     >         threshfreq*autocm,fhelmholtz1*AUKCAL,fhelmholtz1
         write(luout,10028) 
     >         threshfreq*autocm,fhelmholtz2*AUKCAL,fhelmholtz2

         write(luout,10024) ezero*AUKCAL,ezero
         write(luout,10027) threshfreq*autocm,ezero1*AUKCAL,ezero1
         write(luout,10028) threshfreq*autocm,ezero2*AUKCAL,ezero2

         write(luout,10025) ethermal*AUKCAL,ethermal
         write(luout,10027) threshfreq*autocm,ethermal1*AUKCAL,ethermal2
         write(luout,10028) threshfreq*autocm,ethermal2*AUKCAL,ethermal2


         write(luout,10030) Svib   *AUKCAL*1000.0d0
         write(luout,10032) threshfreq*autocm,Svib1*AUKCAL*1000.0d0
         write(luout,10033) threshfreq*autocm,Svib2*AUKCAL*1000.0d0

         write(luout,10031) Cv_vib*AUKCAL*1000.0d0
         write(luout,10032) threshfreq*autocm,Cv_vib1*AUKCAL*1000.0d0
         write(luout,10033) threshfreq*autocm,Cv_vib2*AUKCAL*1000.0d0
         write(luout,*)

      end do !it loop

*     *** free heap ***
      if (.not. MA_free_heap(h_actlist))
     >call errquit('h_actlist:error freeing heap memory',0, MA_ERR)

      return
10012 FORMAT(/6x,'Density of States Details')
10013 FORMAT(6x,'   - grid size (npoints)         = ',I10)
10014 FORMAT(6x,'   - minimum frequency w         = ',f10.1,' cm-1')
10015 FORMAT(6x,'   - maximum frequency w         = ',f10.1,' cm-1')
10016 FORMAT(6x,'   - mininum g(w)                = ',e10.3,
     >       ' states/cm-1')
10017 FORMAT(6x,'   - maximum g(w)                = ',e10.3,
     >       ' states/cm-1')
10018 FORMAT(6x,'   - mininum N(w)                = ',f10.1,' states')
10019 FORMAT(6x,'   - maximum N(w)                = ',f10.1,' states')
10020 FORMAT(/6x,
     > 'Thermodynamic Functions from Vibrational Density of States')
10021 FORMAT(6x,'frequency dividing parameter     = ',f8.1,' cm-1')
10022 FORMAT(6x,'frequency scaling parameter      = ',f8.4)
10023 FORMAT(/6x,'Temperature                      = ',f8.2,'K',
     >       ' (RT=',E13.6,' au)')
10024 FORMAT(/6x,'Zero-Point correction to Energy  = ',f8.3,
     .     ' kcal/mol  (',f10.6,' au)')
10025 FORMAT(/6x,'Thermal correction to Energy     = ',f8.3,
     .     ' kcal/mol  (',f10.6,' au)')
10026 FORMAT(/6x,'Helmholtz Free Energy            = ',f8.3,
     .     ' kcal/mol  (',f10.6,' au)')
10027 FORMAT( 6x,'  - (less than   ',f8.1,' cm-1)  = ',f8.3,
     .     ' kcal/mol  (',f10.6,' au)')
10028 FORMAT( 6x,'  - (greater than',f8.1,' cm-1)  = ',f8.3,
     .     ' kcal/mol  (',f10.6,' au)')


10030 FORMAT(/6x,'Entropy                            = ',f8.3,
     .       ' cal/mol-K')
10031 FORMAT(/6x,'Cv (constant volume heat capacity) = ',f8.3,
     .       ' cal/mol-K')
10032 FORMAT( 6x,'  - (less than   ',f8.1,' cm-1)    = ',f8.3,
     .     ' cal/mol-K')
10033 FORMAT( 6x,'  - (greater than',f8.1,' cm-1)    = ',f8.3,
     .     ' cal/mol-K')
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
            dy = rion(2,jj)-rion(2,ii)
            dz = rion(3,jj)-rion(3,ii)

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

      character*7 function vib_index_name(i)
      integer i

      integer itmp,j0,j1,j2,j3,j4,j5
      character*7 name

      itmp = i

      j5 = itmp/100000
      itmp = itmp - j5*100000
      j4 = itmp/10000
      itmp = itmp - j4*10000
      j3 = itmp/1000
      itmp = itmp - j3*1000
      j2 = itmp/100
      itmp = itmp - j2*100
      j1 = itmp/10
      itmp = itmp - j1*10
      j0 = itmp/1

      name(1:1) = '_'
      name(2:2) = CHAR(j5+ICHAR('0'))
      name(3:3) = CHAR(j4+ICHAR('0'))
      name(4:4) = CHAR(j3+ICHAR('0'))
      name(5:5) = CHAR(j2+ICHAR('0'))
      name(6:6) = CHAR(j1+ICHAR('0'))
      name(7:7) = CHAR(j0+ICHAR('0'))
      vib_index_name = name
      return
      end

*     *************************************************
*     *                                               *
*     *                 vib_writedmat                 *
*     *                                               *
*     *************************************************

* generates phonon animations

      subroutine vib_writedmat(rtdb,ik,step,unita,q,
     >                         nion,rion,nion3,dmat)
      implicit none
      integer   rtdb
      integer   ik
      real*8    step
      real*8    unita(3,3),q(3)
      integer    nion
      real*8     rion(3,nion)
      integer    nion3
      complex*16 dmat(nion3,nion3)

      integer    jjj,ii,i1,i2,i3,k
      real*8     x,y,z,qr,dx,dy,dz,dd
      complex*16 cfacp,cfacm,cpx,cpy,cpz,cmx,cmy,cmz
      character*80  fname
      character*255 full_filename
      character*2 symbol

      character*2 vib_getsymbol
      external    vib_getsymbol
      character*7 vib_index_name
      external    vib_index_name

      do jjj=1,nion3
         fname="kcos"//vib_index_name(ik)//".mode"//vib_index_name(jjj)
     >               //".xyz"
         call util_file_name(fname,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         open(unit=77,file=full_filename,form='formatted')
         fname="ksin"//vib_index_name(ik)//".mode"//vib_index_name(jjj)
     >               //".xyz"
         call util_file_name(fname,
     >                    .false.,
     >                    .false.,
     >                    full_filename)
         open(unit=78,file=full_filename,form='formatted')

         do k=-10,10
            dd = k*step
            write(77,*) nion*27
            write(78,*) nion*27
            write(77,*) 
            write(78,*) 
            do ii=1,nion
              symbol = vib_getsymbol(rtdb,ii)
              do i3=-1,1
              do i2=-1,1
              do i1=-1,1
                x = rion(1,ii)+unita(1,1)*i1+unita(1,2)*i2+unita(1,3)*i3
                y = rion(2,ii)+unita(2,1)*i1+unita(2,2)*i2+unita(2,3)*i3
                z = rion(3,ii)+unita(3,1)*i1+unita(3,2)*i2+unita(3,3)*i3
                qr = q(1)*x + q(2)*y + q(3)*z
                cfacp = dcmplx(dcos(qr),dsin(qr))
                cpx = dmat(3*(ii-1)+1,jjj)*cfacp
                cpy = dmat(3*(ii-1)+2,jjj)*cfacp
                cpz = dmat(3*(ii-1)+3,jjj)*cfacp
                dx = dble(cpx)*dd
                dy = dble(cpy)*dd
                dz = dble(cpz)*dd
                write(77,'(A2,5x,3F12.6)') symbol,
     >                                    (x+dx)*0.529177, 
     >                                    (y+dy)*0.529177, 
     >                                    (z+dz)*0.529177
                dx = dimag(cpx)*dd
                dy = dimag(cpy)*dd
                dz = dimag(cpz)*dd
                write(78,'(A2,5x,3F12.6)') symbol,
     >                                    (x+dx)*0.529177, 
     >                                    (y+dy)*0.529177, 
     >                                    (z+dz)*0.529177
              end do
              end do
              end do
            end do
         end do
         close(77)
         close(78)
      end do

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
            dy = rion(2,jj)-rion(2,ii)
            dz = rion(3,jj)-rion(3,ii)

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
*     *                 vib_getsymbol                 *
*     *                                               *
*     *************************************************

      character*2 function vib_getsymbol(rtdb,ii)
      implicit none
      integer rtdb,ii

#include "stdio.fh"
#include "geom.fh"
#include "errquit.fh"

*     **** local variables ****
      integer      geom
      character*2  symbol
      character*16 t,tname
      real*8 q

      if (.not.geom_create(geom,'geometry')) call errquit
     >    ('vib_getsymbol: error creating geometry',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     >    ('vib_getsymbol: error loading geometry',911, RTDB_ERR)
      
      if (.not. geom_cent_tag(geom,ii,t))
     >        call errquit('vib_getsymbol failed ',ii, GEOM_ERR)
      if (.not. geom_tag_to_element(t,symbol,tname,q))
     >        call errquit('vib_getsymbol failed ',ii, GEOM_ERR)

      if (.not.geom_destroy(geom)) call errquit
     >    ('vib_getsymbol: geom_destroy failed',911, GEOM_ERR)


      vib_getsymbol = symbol
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

      subroutine vib_dos_generate(idx,idy,idz,eigs,neigs,
     >                            npoints,emin,emax,l_unitg,efn)
      implicit none
      integer idx,idy,idz
      real*8 eigs(idx,idy,idz,*)
      integer neigs
      integer npoints
      real*8 emin,emax
      real*8 l_unitg(3,3)
      real*8 efn(npoints,3)

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

        efn(k,1) = e
        efn(k,2) = f
        efn(k,3) = g
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

c $Id$
