      SUBROUTINE vib_WRTFREQ(rtdb,EIGVAL,NAT3,ZEROPE,NPRI)
* $Id: vib_wrtfreq.F,v 1.8 2002-10-08 00:36:35 edo Exp $
C
C This routine scales the eigenvalues properly and then converts them
C      to wavenumbers
C
      IMPLICIT none
c
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer nat3
      LOGICAL ZEROPE
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
      double precision Strans,Srot,Svib
      integer sigma
      integer i

*     ****  conversions and constants ****
      double precision AUKCAL,c,h,kgas,Rgas
      parameter (AUKCAL=627.5093314d0)
      parameter (c=2.998e+10,h=6.626e-27,kgas=1.3807e-16)   ! cgs units
      parameter (Rgas = 1.9863d0/1000.0d0/AUKCAL)           ! atomic units

c
*rak:      ITAP96 = 96
c: error check :
      if (nat3.gt.maxfrq) then
c        dimension problems
        write(luout,*)' maxfrq = ',maxfrq
        write(luout,*)' nat3   = ',nat3
        write(luout,*)' increase maxfrq to ',nat3
        call errquit('vib_wrtfreq:',nat3)
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
        WRITE (6,10002)ZPE
        call ecce_print1('zero point energy',mt_dbl,zpe,1)
        if(.not. rtdb_put(rtdb,'vib:zpe',MT_DBL,1,zpe))
     &      call errquit('vib_wrtfreq:rtdb_put of ZPE failed',
     &      555)
      ENDIF
10002 FORMAT(1X,//,'From the projected analysis ',/,
     &    'The Zero-Point Energy (Kcal/mol) = ',
     &    f20.8)


*     ********************************************************
*     **** Write out thermal Corrections to Energies      ****
*     **** frequencies are assumed to be in units of cm-1 ****
*     ********************************************************
      IF(ZEROPE)THEN

*       *** define temperature ***
        if(.not. rtdb_get(rtdb,'vib:temperature',MT_DBL,1,temp))
     &      temp = 298.15d0

        ezero    = 0.0d0
        ethermal = 0.0d0

*       **** add vibrational thermal corrections ***
        do i=1,nat3
            
          freq = eigval(i) ! eigenvalue

          if (freq.gt.zero) then
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

        ethermal = ethermal + 3.0d0*Rgas*temp ! add translational(3/2RT) and rotation(3/2RT) thermal corrections
        hthermal = ethermal + Rgas*temp       ! htermal = ethermal+pV=ethermal+RT, since pV=RT
       
        WRITE (6,10003) temp
        WRITE (6,10004) ezero*AUKCAL,   ezero
        WRITE (6,10005) ethermal*AUKCAL,ethermal
        WRITE (6,10006) hthermal*AUKCAL,hthermal
10003 FORMAT(/1X,'Temperature                      = ',f8.2,'K')
10004 FORMAT(/1X,'Zero-Point correction to Energy  = ',f8.3, 
     .     ' kcal/mol  (',f10.6,' au)')
10005 FORMAT( 1X,'Thermal correction to Energy     = ',f8.3, 
     .     ' kcal/mol  (',f10.6,' au)')
10006 FORMAT( 1X,'Thermal correction to Enthalpy   = ',f8.3, 
     .     ' kcal/mol  (',f10.6,' au)')

*     ********************************************************
*     **** Write out Entropies                            ****
*     **** frequencies are assumed to be in units of cm-1 ****
*     **** Formulas taken from Herzberg                   ****
*     ********************************************************
        call vib_geom_parameters_get(rtdb,M,
     ,     inertia_A,inertia_B,inertia_C,sigma)
        thetaA = inertia_A*(h*c/kgas)
        thetaB = inertia_B*(h*c/kgas)
        thetaC = inertia_C*(h*c/kgas)

        Strans = (5.0d0/2.0d0)*Rgas*log(temp) + 
     +       (3.0d0/2.0d0)*Rgas*log(M) 
     >         - 2.3135d0/AUKCAL/1000.0d0

        Srot   = (Rgas/2.0d0)*( 3.0d0*log(temp) 
     >                        - log(inertia_A*inertia_B*inertia_C) 
     >                        - 2.0d0*log(dble(sigma))) 
     >         + 3.0327d0/AUKCAL/1000.0d0

*       **** add up vibrational contributions ****
        Svib   = 0.0d0
        if (temp.gt.0.0d0) then
        do i=1,nat3

          freq = eigval(i) ! eigenvalue

          if (freq.gt.zero) then
            thetav = freq*(h*c/kgas)    !freqency temperature in Kelvin from cm-1
            xdum   = thetav/temp
            Svib   = Svib + xdum/(exp(xdum)-1.0d0) - 
     -           log(1.0d0-exp(-xdum))
          end if
        end do
        end if
        Svib = Svib*Rgas


        write(6,10010) (Strans+Srot+Svib)*AUKCAL*1000.0d0
        write(6,10011) Strans*AUKCAL*1000.0d0,M
        write(6,10012) Srot*AUKCAL*1000.0d0,sigma
        write(6,10013) Svib*AUKCAL*1000.0d0
10010   FORMAT(/1X,'Total Entropy                    = ',f8.3,
     .       ' cal/mol-K')
10011   FORMAT( 1X,'  - Translational                = ',f8.3,
     .       ' cal/mol-K (mol. weight =',f9.4,')')
10012   FORMAT( 1X,'  - Rotational                   = ',f8.3,
     .       ' cal/mol-K (symmetry #  =',I9,  ')')
10013   FORMAT( 1X,'  - Vibrational                  = ',f8.3,
     .       ' cal/mol-K')

      END IF
*rak:      IF (NPRI.EQ.0) WRITE (ITAP96,10002)ZPE
*rak:      GO TO 00096
*rak:      ELSE
*rak:      END IF
C
00096 RETURN                    ! go home
      END

      subroutine vib_geom_parameters_get(rtdb,M,A,B,C,sigma)
      implicit none
      integer rtdb
      double precision M,A,B,C
      integer          sigma


#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

*     **** local variables ****
      integer          geom,natoms,ii
      double precision amass,mtensor 

      sigma = 1

      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('vib_geom_parameters_get: error creating geometry',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('vib_geom_parameters_get: error loading geometry',911)
      if (.not.geom_ncent(geom,natoms)) call errquit
     &    ('vib_geom_parameters_get: errorgetting number of atoms',911)

c
*     ***** calculate molecular weight ****
      M = 0.0d0     ! M = 32.0261d0
      do ii=1,natoms
         if(.not.geom_mass_get(geom, ii, amass)) call
     &        errquit(' mass_get  failed ',ii)
         M = M + amass
      enddo

*     ***** calculate inertia tensor ****
      A = 4.194766472d0  ! not implemented yet
      B = 0.826501251d0
      C = 0.797463164d0
*
      if (.not.geom_destroy(geom)) call errquit
     &    ('vib_geom_parameters_get: geom_destroy failed',911)

      return
      end
