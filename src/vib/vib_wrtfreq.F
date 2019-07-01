      SUBROUTINE vib_WRTFREQ(rtdb,EIGVAL,NAT3,ZEROPE,NPRI)
* $Id$
C
C This routine scales the eigenvalues properly and then converts them
C      to wavenumbers
C
      IMPLICIT none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer nat3
      LOGICAL ZEROPE,isatom
      integer npri
      double precision EIGVAL(NAT3)    ! eigenvalues
c
      integer maxfrq
      PARAMETER (MAXFRQ = 2000)
      double precision d1, d2
      PARAMETER(D1=1.0D00,D2=2.0D00)
      double precision ams, wave
      integer numans
      COMMON /cvib_SETCON/ AMS(36),WAVE,NUMANS      ! setup parameters
      double precision PRFREQ(MAXFRQ)    ! freqencies to be written to tape10
c
      double precision scale, dd1, zero, sumfreq, rate, xdum, freq, zpe
      double precision ezero,ethermal,hthermal,temp,thetav
      double precision inertia_A,inertia_B,inertia_C,thetaA,thetaB,
     ,     thetaC,M
      double precision Stotal,Strans,Srot,Svib,scalefreq
      double precision Cv_total,Cv_trans,Cv_rot,Cv_vib
      integer sigma
      integer i,it,ntemp,h_actlist,l_actlist,ma_type
      logical linear

*     ****  conversions and constants ****
      double precision AUKCAL,c,h,kgas,Rgas,inertprod
      parameter (AUKCAL=627.5093314d0)
      parameter (c=2.998e+10,h=6.626e-27,kgas=1.3807e-16)   ! cgs units
      parameter (Rgas = 1.9863d0/1000.0d0/AUKCAL)           ! atomic units

c
      isatom = (nat3.eq.3)

*rak:      ITAP96 = 96
c: error check :
      if (nat3.gt.maxfrq) then
c        dimension problems
        write(luout,*)' maxfrq = ',maxfrq
        write(luout,*)' nat3   = ',nat3
        write(luout,*)' increase maxfrq to ',nat3
        call errquit('vib_wrtfreq:',nat3, INPUT_ERR)
        stop ' error'
      endif
C
C Scale eigenvalues down by 1000
C
      SCALE = 1.0D-03   ! Scaling factor after diagonaization
      call dscal(nat3,scale,eigval,1) ! Scale eigenvalues
C
      DD1 = D1          ! used with dsign function
C      WRITE(6,10000)    ! write header
10000 FORMAT(//,1X,31('-'),' Vibrations ',31('-'),/,
     +5X,'#',6X,'Eigenvalue',14X,'Frequency (cm**-1)',/,1X,74('-'))
10001 FORMAT(1X,I5,D20.10,F20.3)
      ZERO=0.0D00
      SUMFREQ=zero
      RATE=2.860D-03    ! Note that there is a slight error in this conversion - RATE = R*h*c/k = 2.857778587 kcal/mol 
      DO 00100 I=1,NAT3 ! loop over coordinates
         XDUM = EIGVAL(I) ! eigenvalue
         FREQ = SIGN(DD1,XDUM)*SQRT(ABS(XDUM))*WAVE ! form frequency
         PRFREQ(I) = FREQ
C         WRITE(6,10001,ERR=00099)I,XDUM,FREQ ! write it out
         EIGVAL(I) = FREQ
         IF(ZEROPE) THEN
            IF (FREQ.GT.ZERO)THEN
              SUMFREQ=SUMFREQ+FREQ
            ENDIF
         END IF
         GO TO 00100
C00099    WRITE(6,*)EIGVAL(I),FREQ  ! in case of error with format
00100 CONTINUE
      IF(ZEROPE)THEN
        ZPE=0.5D00*SUMFREQ*RATE

         !**** removed bad ZPE printout ...EJB***
c        WRITE (6,10002)ZPE
c        call ecce_print1('zero point energy',mt_dbl,zpe,1)
c        if(.not. rtdb_put(rtdb,'vib:zpe',MT_DBL,1,zpe))
c     &      call errquit('vib_wrtfreq:rtdb_put of ZPE failed',
c     &      555, RTDB_ERR)
      ENDIF
10002 FORMAT(1X,//,'From the projected analysis ',/,
     &    'The Zero-Point Energy (Kcal/mol) = ',
     &    f20.8)


*     ********************************************************
*     **** Cannoncial Partition Equations for Ideal       ****
*     **** polyatomic gas                                 ****
*     ********************************************************
      IF(ZEROPE)THEN

*       *** define rotational constants ****
        call vib_geom_parameters_get(rtdb,M,
     >     inertia_A,inertia_B,inertia_C,sigma,linear)
        thetaA = inertia_A*(h*c/kgas)
        thetaB = inertia_B*(h*c/kgas)
        thetaC = inertia_C*(h*c/kgas)
        write(6,10010) 
        write(6,10011) inertia_A,thetaA
        write(6,10012) inertia_B,thetaB
        write(6,10013) inertia_C,thetaC
10010 FORMAT(/1X,'Rotational Constants',
     >       /1X,'--------------------')
10011 FORMAT( 1X,'A= ',f10.6,' cm-1  (',f10.6,' K)')
10012 FORMAT( 1X,'B= ',f10.6,' cm-1  (',f10.6,' K)')
10013 FORMAT( 1X,'C= ',f10.6,' cm-1  (',f10.6,' K)')
       
*       *** define scalefreq ***
        if(.not. rtdb_get(rtdb,'vib:scalefreq',MT_DBL,1,scalefreq))
     &      scalefreq = 1.0d0

*        **** define temperatures
         if (.not. rtdb_ma_get(rtdb, 'vib:temperature', ma_type,
     >        ntemp, h_actlist)) then
            ntemp = 1

            if (.not.MA_alloc_get(mt_dbl,(ntemp),'temp',
     >                            h_actlist,l_actlist))
     >      call errquit(
     >       'vib_WRTFREQ: ma_get_index failed for actlist',911, MA_ERR)

            dbl_mb(l_actlist) = 298.15d0
         else        

            if (.not.MA_get_index(h_actlist,l_actlist))
     >      call errquit(
     >       'vib_WRTFREQ: ma_get_index failed for actlist',911, MA_ERR)
         end if

        do it=1,ntemp
        temp = dbl_mb(l_actlist+it-1)

*       ********************************************************
*       **** Write out thermal Corrections to Energies      ****
*       **** frequencies are assumed to be in units of cm-1 ****
*       ********************************************************
        ezero    = 0.0d0
        ethermal = 0.0d0

*       **** add vibrational thermal corrections ***
        do i=1,nat3
            
          freq = scalefreq*eigval(i) ! eigenvalue

          if (freq.gt.0.1d0) then
            thetav = freq*(h*c/kgas)    !freqency temperature in Kelvin from cm-1
            if (temp.gt.zero) then
              xdum   = exp(-thetav/temp)
            else
              xdum = 0.0d0
            end if
            xdum   = xdum/(1.0d0-xdum)
            ezero    = ezero    + 0.5d0*thetav
            ethermal = ethermal + thetav*(0.5d0 + xdum)
          end if
        end do
        ezero    = ezero*Rgas
        ethermal = ethermal*Rgas
 
        if (isatom) then
              ethermal = ethermal + 1.5d0*Rgas*temp ! add only translational(3/2RT) thermal corrections
        else
           if (linear) then
              ethermal = ethermal + 2.5d0*Rgas*temp ! add translational(3/2RT) and rotation(2/2RT) thermal corrections
           else
              ethermal = ethermal + 3.0d0*Rgas*temp ! add translational(3/2RT) and rotation(3/2RT) thermal corrections
           end if
        end if

        hthermal = ethermal + Rgas*temp       ! htermal = ethermal+pV=ethermal+RT, since pV=RT
       
        WRITE (6,10023) temp
        WRITE (6,10022) scalefreq
        if (isatom) then
           write(6,10028)
        else
           if (linear) then
              write(6,10027)
              if (.not.rtdb_put(rtdb,'vib:linear',mt_log,1,.true.))
     &           call errquit('vib: failed to set linear',555, RTDB_ERR)
           end if
        end if
 
        if(.not. rtdb_put(rtdb,'vib:zpe',MT_DBL,1,ezero))
     &       call errquit('vib_wrtfreq:rtdb_put of ZPE failed',
     &       555, RTDB_ERR)
        if(.not. rtdb_put(rtdb,'vib:ethermal',MT_DBL,1,ethermal))
     &       call errquit('vib_wrtfreq:rtdb_put of ethermal failed',
     &       555, RTDB_ERR)
        if(.not. rtdb_put(rtdb,'vib:hthermal',MT_DBL,1,hthermal))
     &       call errquit('vib_wrtfreq:rtdb_put of hthermal failed',
     &       555, RTDB_ERR)
 
        WRITE (6,10024) ezero*AUKCAL,   ezero
        WRITE (6,10025) ethermal*AUKCAL,ethermal
        WRITE (6,10026) hthermal*AUKCAL,hthermal
10022 FORMAT(1X,'frequency scaling parameter      = ',f8.4)
10023 FORMAT(//1X,'Temperature                      = ',f8.2,'K')
10024 FORMAT(/1X,'Zero-Point correction to Energy  = ',f8.3, 
     .     ' kcal/mol  (',f10.6,' au)')
10025 FORMAT( 1X,'Thermal correction to Energy     = ',f8.3, 
     .     ' kcal/mol  (',f10.6,' au)')
10026 FORMAT( 1X,'Thermal correction to Enthalpy   = ',f8.3, 
     .     ' kcal/mol  (',f10.6,' au)')
10027 FORMAT(//1X,'Linear Molecule                   ')
10028 FORMAT(//1X,'Atom                              ')

*     ********************************************************
*     **** Write out Entropies                            ****
*     **** frequencies are assumed to be in units of cm-1 ****
*     **** Formulas taken from Herzberg                   ****
*     ********************************************************

        Cv_trans = (3.0d0/2.0d0)*Rgas
        Strans = (5.0d0/2.0d0)*Rgas*log(temp)
     >         + (3.0d0/2.0d0)*Rgas*log(M) 
     >         - 2.3135d0/AUKCAL/1000.0d0

        if (isatom) then
           Cv_rot = 0.0d0
           Srot = 0.0d0
        else
           inertprod=inertia_A*inertia_B*inertia_C
           if (dabs(inertprod).lt.1.0d-9) then
              Cv_rot = (2.0d0/2.0d0)*Rgas
              Srot=0.0d0
              if (linear) then
                  Srot   = log(temp) 
     >                           - log(dble(sigma))
                 if(inertia_C.gt.1d-9) 
     >              Srot= Srot - log(inertia_C) 
                  Srot   = Rgas*Srot
     >                    + 1.2639d0/AUKCAL/1000.0d0
              end if
           else
              Cv_rot = (3.0d0/2.0d0)*Rgas
              Srot   = (Rgas/2.0d0)*( 3.0d0*log(temp) -
     >             log(inertprod) -
     >             2.0d0*log(dble(sigma))) +
     >              3.0327d0/AUKCAL/1000.0d0
           endif
        end if

*       **** add up vibrational contributions ****
        Cv_vib = 0.0d0
        Svib   = 0.0d0
        if (temp.gt.0.0d0) then
        do i=1,nat3

          freq = scalefreq*eigval(i) ! eigenvalue

          if (freq.gt.(0.1d0)) then
            thetav = freq*(h*c/kgas)    !freqency temperature in Kelvin from cm-1
            xdum   = thetav/temp
c  fix to take avoid abnormal values
            if(xdum.lt.20d0*2.3025d0) then
            Cv_vib = Cv_vib + exp(xdum)* (xdum/(exp(xdum)-1.0d0))**2 
            Svib   = Svib + xdum/(exp(xdum)-1.0d0)   
     >             - log(1.0d0-exp(-xdum))
            endif
          end if
        end do
        end if
        Cv_vib = Cv_vib*Rgas
        Svib   = Svib*Rgas


        Stotal = (Strans+Srot+Svib)*AUKCAL*1000.0d0
        Cv_total = (Cv_trans+Cv_rot+Cv_vib)*AUKCAL*1000.0d0
        
        if(.not. rtdb_put(rtdb,'vib:entropy',MT_DBL,1,Stotal))
     &       call errquit('vib_wrtfreq:rtdb_put of entropy failed',
     &       555, RTDB_ERR)

        if(.not. rtdb_put(rtdb,'vib:Cv',MT_DBL,1,Cv_total))
     &       call errquit('vib_wrtfreq:rtdb_put of Cv failed',
     &       555, RTDB_ERR)
        
        write(6,10030) Stotal
*        write(6,10030) (Strans+Srot+Svib)*AUKCAL*1000.0d0
        write(6,10031) Strans*AUKCAL*1000.0d0,M
        write(6,10032) Srot*AUKCAL*1000.0d0,sigma
        write(6,10033) Svib*AUKCAL*1000.0d0
10030   FORMAT(/1X,'Total Entropy                    = ',f8.3,
     .       ' cal/mol-K')
10031   FORMAT( 1X,'  - Translational                = ',f8.3,
     .       ' cal/mol-K (mol. weight =',f9.4,')')
10032   FORMAT( 1X,'  - Rotational                   = ',f8.3,
     .       ' cal/mol-K (symmetry #  =',I9,  ')')
10033   FORMAT( 1X,'  - Vibrational                  = ',f8.3,
     .       ' cal/mol-K')

        write(6,10040) Cv_total
*        write(6,10040) (Cv_trans+Cv_rot+Cv_vib)*AUKCAL*1000.0d0
        write(6,10041) Cv_trans*AUKCAL*1000.0d0
        write(6,10042) Cv_rot*AUKCAL*1000.0d0
        write(6,10043) Cv_vib*AUKCAL*1000.0d0
10040   FORMAT(/1X,'Cv (constant volume heat capacity) = ',f8.3,
     .       ' cal/mol-K')
10041   FORMAT( 1X,'  - Translational                  = ',f8.3,
     .       ' cal/mol-K')
10042   FORMAT( 1X,'  - Rotational                     = ',f8.3,
     .       ' cal/mol-K')
10043   FORMAT( 1X,'  - Vibrational                    = ',f8.3,
     .       ' cal/mol-K')

       end do !it loop

*      *** free heap ***
       if (.not. MA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)

      END IF
*rak:      IF (NPRI.EQ.0) WRITE (ITAP96,10002)ZPE
*rak:      GO TO 00096
*rak:      ELSE
*rak:      END IF
C
00096 RETURN                    ! go home
      END

      subroutine vib_geom_parameters_get(rtdb,M,AA,BB,CC,sigma,linear)
      implicit none
#include "errquit.fh"
      integer rtdb
      double precision M,AA,BB,CC
      integer          sigma
      logical          linear


#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "sym.fh"

*     **** local variables ****
      integer          geom,natoms,ii,nops
      double precision amass,meig(3),mtensor(3,3),q,qt
      integer          r1(2)
      integer          lwork,ierr
      parameter (lwork=9)
      double precision work(lwork)
      character*16 t
      double precision fconv,eps

*     ****  conversions and constants ****
      double precision AUKCAL,c,h,kgas,Rgas
      parameter (AUKCAL=627.5093314d0)
      parameter (c=2.998e+10,h=6.626e-27,kgas=1.3807e-16)   ! cgs units
      parameter (Rgas = 1.9863d0/1000.0d0/AUKCAL)           ! atomic units
      parameter (fconv=4.6497791d-41, eps=1d-9)

*     **** external functions ****
      double precision deter3
      external         deter3


      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('vib_geom_parameters_get: error creating geometry',911,
     &       GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('vib_geom_parameters_get: error loading geometry',911,
     &       RTDB_ERR)
      if (.not.geom_ncent(geom,natoms)) call errquit
     &    ('vib_geom_parameters_get: errorgetting number of atoms',911,
     &       GEOM_ERR)

*     **** calculate the symmetry number ****
      nops  = sym_number_ops(geom)
      if (.not. MA_push_get(mt_dbl,(nops*3*4),'r1',r1(2),r1(1)))
     >   call errquit('vib_geom_parameters_set: push stack',0, MA_ERR)
      call sym_ops_get(geom,nops,dbl_mb(r1(1)))
      sigma = 1                                                 ! Identity not included
      do ii=1,nops
        mtensor(1,1) = dbl_mb(r1(1)+ (ii-1)*3 + 0 + 0*(nops*3)) !crazy format
        mtensor(2,1) = dbl_mb(r1(1)+ (ii-1)*3 + 1 + 0*(nops*3))
        mtensor(3,1) = dbl_mb(r1(1)+ (ii-1)*3 + 2 + 0*(nops*3))
        mtensor(1,2) = dbl_mb(r1(1)+ (ii-1)*3 + 0 + 1*(nops*3))
        mtensor(2,2) = dbl_mb(r1(1)+ (ii-1)*3 + 1 + 1*(nops*3))
        mtensor(3,2) = dbl_mb(r1(1)+ (ii-1)*3 + 2 + 1*(nops*3))
        mtensor(1,3) = dbl_mb(r1(1)+ (ii-1)*3 + 0 + 2*(nops*3))
        mtensor(2,3) = dbl_mb(r1(1)+ (ii-1)*3 + 1 + 2*(nops*3))
        mtensor(3,3) = dbl_mb(r1(1)+ (ii-1)*3 + 2 + 2*(nops*3))
        q = deter3(mtensor) 
        if (abs(q-1.0d0) .lt. 1.0d-6) sigma = sigma+1
      end do

      if (.not. MA_pop_stack(r1(2)))
     >  call errquit('vib_geom_parameters_set: pop stack',1, MA_ERR)


*     ***** calculate molecular weight ****
      M = 0.0d0     
      do ii=1,natoms
         if(.not.geom_mass_get(geom, ii, amass)) call
     &        errquit(' mass_get  failed ',ii, GEOM_ERR)
         M = M + amass
      enddo

*     ***** calculate inertia tensor ****
      if (.not. MA_push_get(mt_dbl,(3*natoms),'r1',r1(2),r1(1)))
     >   call errquit('vib_geom_parameters_set: push stack',2, MA_ERR)

      call dcopy(3*natoms,0.0d0,0,dbl_mb(r1(1)),1)
      do ii=1,natoms
         if (.not.geom_cent_get(geom,ii,t,dbl_mb(r1(1)+(ii-1)*3),q))
     >   call errquit('vib_geom_parameters_set: error reading ions',0,
     &       GEOM_ERR)
      end do

      call geom_momint0(geom,dbl_mb(r1(1)),natoms,meig,mtensor,.true.,
     .     .false.,.false.)

      if (.not. MA_pop_stack(r1(2)))
     >  call errquit('vib_geom_parameters_set: pop stack',3, MA_ERR)
      call DSYEV('N','U',3,mtensor,3,meig,work,lwork,ierr)

      !meig(1) = fconv*meig(1)  ! convert to cgs units
      !meig(2) = fconv*meig(2)
      !meig(3) = fconv*meig(3)

      q = 4.0d0*datan(1.0d0)
      q = 8.0d0*q*q             ! 8*pi^2
      q = h/(q*c)               ! h/(8*pi^2*c) - there is an error in old McQuarrie

      q = q/fconv               ! convert to cgs units

      AA = 0.0d0
      BB = 0.0d0
      CC = 0.0d0
      if(abs(meig(1)).gt.eps) AA = q/meig(1)    ! in units of cm-1
      if(abs(meig(2)).gt.eps) BB = q/meig(2)
      if(abs(meig(3)).gt.eps) CC = q/meig(3)

      linear = .false.
      if ((AA.le.eps).and.((BB-CC).le.eps)) linear = .true.


*     **** recalculate the symmetry number  - Dinfh = 2, Cinfh = 1 ****
      if (linear) then

        nops  = sym_number_ops(geom)
        if (.not. MA_push_get(mt_dbl,(nops*3*4),'r1',r1(2),r1(1)))
     >     call errquit('vib_geom_parameters_set: push stack',4, MA_ERR)
        call sym_ops_get(geom,nops,dbl_mb(r1(1)))
        sigma = 1                                                 ! Identity not included
        do ii=1,nops
          mtensor(1,1) = dbl_mb(r1(1)+ (ii-1)*3 + 0 + 0*(nops*3)) !crazy format
          mtensor(2,1) = dbl_mb(r1(1)+ (ii-1)*3 + 1 + 0*(nops*3))
          mtensor(3,1) = dbl_mb(r1(1)+ (ii-1)*3 + 2 + 0*(nops*3))
          mtensor(1,2) = dbl_mb(r1(1)+ (ii-1)*3 + 0 + 1*(nops*3))
          mtensor(2,2) = dbl_mb(r1(1)+ (ii-1)*3 + 1 + 1*(nops*3))
          mtensor(3,2) = dbl_mb(r1(1)+ (ii-1)*3 + 2 + 1*(nops*3))
          mtensor(1,3) = dbl_mb(r1(1)+ (ii-1)*3 + 0 + 2*(nops*3))
          mtensor(2,3) = dbl_mb(r1(1)+ (ii-1)*3 + 1 + 2*(nops*3))
          mtensor(3,3) = dbl_mb(r1(1)+ (ii-1)*3 + 2 + 2*(nops*3))
          q = deter3(mtensor)
          qt = mtensor(1,1)+mtensor(2,2)+mtensor(3,3)
          if ((abs(q+1.0d0).lt.1.0d-6).and.(abs(qt+3.0d0).lt.1.0d-6))
     >       sigma = 2
        end do
        if (.not. MA_pop_stack(r1(2)))
     >    call errquit('vib_geom_parameters_set: pop stack',5, MA_ERR)

      end if
*
      if (.not.geom_destroy(geom)) call errquit
     &    ('vib_geom_parameters_get: geom_destroy failed',911, GEOM_ERR)

      return
      end
