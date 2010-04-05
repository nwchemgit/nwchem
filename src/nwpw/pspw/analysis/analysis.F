c
c $Id: analysis.F,v 1.12 2008-11-03 06:06:11 bylaska Exp $
c

      subroutine pspw_analysis(flag,rtdb,ispin,ne,psi,eig)
      implicit none
      integer flag
      integer rtdb
      integer ispin,ne(2)
      complex*16 psi(*)
      real*8    eig(*)

#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "rtdb.fh"
#include "errquit.fh"



* ================================================
*   This code is a rewrite of an earlier analysis
* code written by R. Kawai.
*
*     VERSION  MPI-1.00   11/15/96 by Eric Bylaska
*
* ================================================
 
*     **** parallel variables ****
      integer  taskid
      integer  MASTER
      parameter(MASTER=0)


*     **** electronic variables ****
      integer npack1
      integer n1(2),n2(2)
    


*     **** local variables ****
      integer i,k,l,n,ms,l1,l2,j
      integer ll,spin,ind
      real*8  ttl1,ttl2,subttl

      integer h_actlist,l_actlist,nactive_atoms,ma_type
      integer npoints,ii,lmax
      integer weight(2)
      real*8 emin,emax,alpha
      character*255 filename

      logical value,mulliken_kawai,fixatoms
      character*28 DD
      character*255 id,test
      integer npsp,nion,nemax
      integer lorb(2)    ! integer lorb(npsp)
      integer b0(2)      ! real*8 b0(0:5,npsp)
      integer a(2)       ! real*8 a(36,nemax,nion)
      integer total(2)   ! real*8 total(nion)
      integer sum(2)     ! real*8 sum(nemax)
      integer subtl(2)   ! real*8 subtl(0:5,3)
 
      character*4 spn(2)
      DATA SPN / 'UP  ', 'DOWN' /


*     **** external functions ****
      logical  control_DOS,nwpw_filefind
      logical  control_mulliken_kawai,aorbs_init,aorbs_readall
      external control_DOS,nwpw_filefind
      external control_mulliken_kawai,aorbs_init,aorbs_readall
      character   spdf_name
      external    spdf_name
      character*4 ion_atom_qm
      external    ion_atom_qm
      integer  ion_nion_qm,ion_katm_qm,ion_nkatm_qm
      external ion_nion_qm,ion_katm_qm,ion_nkatm_qm
      real*8   ion_rion
      external ion_rion
      real*8      lattice_omega,lattice_ecut,lattice_unita
      external    lattice_omega,lattice_ecut,lattice_unita


      call Parallel_taskid(taskid)
      call Pack_npack(1,npack1)
      mulliken_kawai = control_mulliken_kawai()

      npsp = ion_nkatm_qm()
      nion = ion_nion_qm()
      nemax = ne(1)+ne(2)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

      value = MA_alloc_get(mt_int,npsp,'lorb',lorb(2),lorb(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,6*npsp,'b0',b0(2),b0(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,36*nemax*nion,'a',a(2),a(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nion,'total',total(2),total(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nemax,'sum',sum(2),sum(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,6*3,'subtl',subtl(2),subtl(1))


*      ***********************************************
*      **** create psp1 files if they don't exist ****
*      ***********************************************
       do k=1,npsp
         DD = '                          '
         DD = ion_atom_qm(k)
         ind = index(DD,' ') - 1
         test = DD(1:ind)//'.psp1'
         if (.not.nwpw_filefind(test))
     >      call aorbs_formatter_auto(ion_atom_qm(k))
       end do

*      ****************************************
*      **** read in expansion coefficients ****
*      ****************************************
       do k=1,npsp
          id = 'analysis:lorb'//ion_atom_qm(k)
          if (.not. rtdb_get(rtdb,id,mt_int,1,int_mb(lorb(1)+k-1))) then
            DD = '                          '
            DD = ion_atom_qm(k)
            ind = index(DD,' ') - 1
            test = DD(1:ind)//'.psp1'

            !write(*,*) "test:",test,ind
            value = rtdb_parallel(.false.)
            if (taskid.eq.MASTER) then
            ind = index(test,' ') - 1
            !write(*,*) "test:",test,ind
            call analysis_expansion_coef(test,-1,rtdb)
            end if
            value = rtdb_parallel(.true.)
            call ga_sync()

            if (.not. rtdb_get(rtdb,id,mt_int,1,int_mb(lorb(1)+k-1))) 
     >        call errquit(
     >        'analysis: rtdb_get lorb failed', 0, RTDB_ERR)
          end if

          id = 'analysis:expansion'//ion_atom_qm(k)
          if (.not. rtdb_get(rtdb,id,mt_dbl,(int_mb(lorb(1)+k-1)+1),
     >                                  dbl_mb(b0(1)+(k-1)*6))) then
            DD = '                          '
            DD = ion_atom_qm(k)
            ind = index(DD,' ') -1
            test = DD(1:ind)//'.psp1'
            call analysis_expansion_coef(test,-1,rtdb)

            if (.not. rtdb_get(rtdb,id,mt_dbl,(int_mb(lorb(1)+k-1)+1),
     >                                       dbl_mb(b0(1)+(k-1)*6))) 
     >       call errquit(
     >       'analysis: rtdb_get failed', 0, RTDB_ERR)
          end if
       end do

      if (.not.mulliken_kawai) then
         value = aorbs_init()
         value = value.and.aorbs_readall()
         if (.not.value) go to 1901
      end if



      if (taskid.eq.MASTER) then
         call util_date(DD)
 
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*) 
     >  '*************************************************************'
         WRITE(6,*) 
     >  '**                                                         **'
         WRITE(6,*) 
     >  '**                  PSPW Mulliken analysis                 **'
         WRITE(6,*) 
     >  '**                                                         **'
         if (flag.eq.1)
     >   WRITE(6,*)
     >  '**                   (Virtual Orbitals)                    **'
         WRITE(6,*) 
     >  '** Population analysis algorithm devloped by Ryoichi Kawai **'
         WRITE(6,*) 
     >  '**                                                         **'
         WRITE(6,1000) DD
 1000    FORMAT(
     >  ' **                   ',A16,'                      **')
         WRITE(6,*)
     >   '**                                                         **'
         WRITE(6,*) 
     >  '*************************************************************'
      end if


c     **** ouput xyz format ****
      call ion_Print_XYZ(6)


      if (taskid.eq.MASTER) then
         if (mulliken_kawai) then
         WRITE(6,1300)
         WRITE(6,1305) 'ATOM ','S','P','D','F'
         DO k=1,npsp
           WRITE(6,1310) ion_atom_qm(K),
     >                  (dbl_mb(b0(1)+l+(k-1)*6), 
     >                   l=0,int_mb(lorb(1)+k-1))
         END DO
         call util_flush(6)
         end if
 1300    FORMAT(//'== Kawai Expansion Coefficients =='/)
 1305    FORMAT(A5,6X,A,12X,A,12X,A,12X,A)
 1310    FORMAT(A4,' : ',4E13.5)
      end if
 

 
      call util_file_name('ORBOUT',
     >                     .true.,
     >                     .false.,
     >                      id)
      if (taskid.eq.MASTER) 
     > OPEN(UNIT=65,FILE=id,FORM='FORMATTED')
      
      call Orb_Analysis(65,flag,ispin,ne,npack1,nemax,psi,
     >                        int_mb(lorb(1)),
     >                        dbl_mb(b0(1)),
     >                        dbl_mb(a(1)),
     >                        dbl_mb(sum(1)))
     
      if (taskid.eq.MASTER) close(unit=65)



      if (taskid.eq.MASTER) then
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '====================================================='
      if (flag.eq.0)
     >WRITE(6,*) '| POPULATION ANALYSIS OF FILLED MOLECULAR ORBITALS  |'
      if (flag.eq.1)
     >WRITE(6,*) '| POPULATION ANALYSIS OF VIRTUAL MOLECULAR ORBITALS |'
      WRITE(6,*) '====================================================='
      if (.not.mulliken_kawai) then
        WRITE(6,1311)
      else
        WRITE(6,1312)
      end if
c     WRITE(6,1313)
c     WRITE(6,1314)
c     WRITE(6,1315)
 1311 FORMAT(//'== Using pseudoatomic orbital expansion          ==')
 1312 FORMAT(//'== Using Kawai projected atomic expansion        ==')
c1313 FORMAT('== order of orbitals: s                          ==')
c1314 FORMAT('==                    px,  py,     pz            ==')
c1315 FORMAT('==                    dzz, dx2-y2, dxy, dyz, dzx ==')
 
      DO SPIN=1,ISPIN
        DO N=N1(SPIN),N2(SPIN)
          WRITE(6,1500)
          IF(ISPIN.EQ.2) THEN
            WRITE(6,1510) N,SPN(SPIN),dbl_mb(sum(1)+N-1),eig(n),
     >                    eig(n)*27.2116d0
          ELSE
            WRITE(6,1515) N,dbl_mb(SUM(1)+N-1),eig(n),
     >                    eig(n)*27.2116d0
          ENDIF
          !write(6,1519)
          WRITE(6,1520) 'NO','ATOM','L','POPULATION'

          DO L=0,5
            dbl_mb(SUBTL(1)+L)=0.0d0
          END DO
          DO I=1,nion
            dbl_mb(TOTAL(1)+I-1)=0.0d0
            DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
              L1=L**2+1
              L2=(L+1)**2
              SUBTTL=0.0d0
              DO LL=L1,L2
                dbl_mb(TOTAL(1)+I-1)=dbl_mb(TOTAL(1)+I-1)+
     >                dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
                SUBTTL=SUBTTL+
     >                dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
              END DO
              dbl_mb(SUBTL(1)+L)=dbl_mb(SUBTL(1)+L)+SUBTTL
              if (l.eq.0) write(6,1516)
              if (l.eq.1) write(6,1517)
              if (l.eq.2) write(6,1518)
              if (l.eq.3) write(6,1519)
              WRITE(6,1530) I,ion_atom_qm(ion_katm_qm(I)),L,SUBTTL,
     >               (dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax),
     >                        LL=L1,L2)

            END DO
          END DO

          WRITE(6,1540)
          WRITE(6,1550) (I,ion_atom_qm(ion_katm_qm(I)),
     >                   dbl_mb(TOTAL(1)+I-1),I=1,nion)
          WRITE(6,1555)
          WRITE(6,1560)  's','p','d','f'
          WRITE(6,1570)  (dbl_mb(SUBTL(1)+L),L=0,3)
        END DO
      END DO
      call util_flush(6)

 1500 FORMAT(//'------------------------------------------------',
     >         '------------------------------'//)
 1510 FORMAT('*** ORBITAL=',I4,'***  SPIN=',A4,4X,'SUM=',E12.5,
     >       ' E=',E12.5,' (',F8.3,'eV)'/)
 1515 FORMAT('*** ORBITAL=',I4,'***  SPIN=BOTH',4X,'SUM=',E12.5,
     >       ' E=',E12.5,' (',F8.3,'eV)'/)
 1516 FORMAT(30x,'  s')
 1517 FORMAT(30x,' px         py         pz')
 1518 FORMAT(27x,'d3z2-1        dxy        dyz        dzx     dx2-y2')
 1519 FORMAT(30x,' f1         f2         f3         f4         f5',
     >           '         f6         f7',
     >      /27x,' x(3x2  y(3(1-z2)   z(x2-y2)        xyz   x(5z2-1)',
     >           '   y(5z2-1)   z(5z2-3)'
     >      /24x,'-3(1-z2))      -4y2)' ) 

c 1519 FORMAT(30x,'  s',
c     >      /30x,' px         py        pz '
c     >      /27x,'d3z2-1        dxy        dyz        dzx     dx2-y2',
c     >      /30x,' f1         f2        f3          f4         f5', 
c     >           '         f6         f7')

 1520 FORMAT(A2,2X,A4,2X,A1,2X,A10)
 1530 FORMAT(I2,3X,A4,3X,I1,8F11.5)
 1531 FORMAT(8F11.5)
 1540 FORMAT(//'=== DISTRIBUTION ==='/)
 1550 FORMAT(4(I6,'(',A4,')',F9.4))
 1555 FORMAT(//'== ANGULAR MOMENTUM POPULATIONS ==='/)
 1560 FORMAT(6X,A1,3(9X,A1))
 1570 FORMAT(4F10.4)
 
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '========================================'
      WRITE(6,*) '|   POPULATION ANALYSIS ON EACH ATOM   |'
      WRITE(6,*) '========================================'
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,1610) 'NO','ATOM','SPIN','TOTAL','s','p','d','f'
      DO I=1,nion
        TTL1=0.0d0
        TTL2=0.0d0

        DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
          DO ms=1,ISPIN
            dbl_mb(SUBTL(1)+L+(ms-1)*6)=0.0d0
            DO N=N1(ms),N2(ms)
              L1=L**2+1
              L2=(L+1)**2
              DO LL=L1,L2
                dbl_mb(SUBTL(1)+L+(ms-1)*6)=
     >            dbl_mb(SUBTL(1)+L+(ms-1)*6) +
     >            dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2

              END DO
            END DO
          END DO
        END DO
     
        TTL1=0.0d0
        DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
          TTL1=TTL1+dbl_mb(SUBTL(1)+L)
        END DO
        WRITE(6,1620) I,ion_atom_qm(ion_katm_qm(I)),SPN(1),TTL1,
     >                 ( dbl_mb(SUBTL(1)+L),
     >                   L=0,int_mb(lorb(1)+ion_katm_qm(I)-1) )
        TTL1=0.0d0
        DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
          TTL1=TTL1+dbl_mb(SUBTL(1)+L+(ispin-1)*6)
        END DO
        WRITE(6,1620) I,ion_atom_qm(ion_katm_qm(I)),SPN(2),TTL1,
     >                 ( dbl_mb(SUBTL(1)+L+(ispin-1)*6),
     >                   L=0,int_mb(lorb(1)+ion_katm_qm(I)-1) )

      END DO
      call util_flush(6)

 1610 FORMAT(A2,2X,A4,2X,A4,4X,A5,7X,A,10X,A,10X,A,10X,A)
 1620 FORMAT(I2,3X,A4,3X,A4,5F11.5)
 
      DO L=0,3
        dbl_mb(SUBTL(1)+L)  =0.0d0
        dbl_mb(SUBTL(1)+L+6)=0.0d0
      END DO
      DO I=1,nion
        DO SPIN=1,ISPIN
          DO N=N1(SPIN),N2(SPIN)
            DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
              L1=L**2+1
              L2=(L+1)**2
              DO LL=L1,L2
                dbl_mb(SUBTL(1)+L+(SPIN-1)*6)=
     >          dbl_mb(SUBTL(1)+L+(SPIN-1)*6)+
     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
              END DO
            END DO
          END DO
        END DO
      END DO

      DO L=0,3
        dbl_mb(SUBTL(1)+L+2*6)=
     >  (dbl_mb(SUBTL(1)+L)+dbl_mb(SUBTL(1)+L+(ISPIN-1)*6))
     >   *100.d0/(NE(1)+NE(ISPIN))
        dbl_mb(SUBTL(1)+L)=dbl_mb(SUBTL(1)+L)*100.0d0/dble(NE(1))
        IF((ISPIN.EQ.2).and.(NE(2).gt.0)) 
     >    dbl_mb(SUBTL(1)+L+6)=dbl_mb(SUBTL(1)+L+6)*100.0d0/dble(NE(2))
      END DO

      WRITE(6,1700)
      WRITE(6,1710) ' SPIN ','s','p','d','f'
      WRITE(6,1720) SPN(1),(dbl_mb(SUBTL(1)+L),L=0,3)
      WRITE(6,1720) SPN(ISPIN),(dbl_mb(SUBTL(1)+L+(ISPIN-1)*6),L=0,3)
      WRITE(6,1720) ' TOTAL',(dbl_mb(SUBTL(1)+L+(3-1)*6),L=0,3)
      call util_flush(6)
 1700 FORMAT(///'=== TOTAL ANGULAR MOMENTUM POPULATION ==='/)
 1710 FORMAT(A6,6X,A1,3(11X,A1))
 1720 FORMAT(A6,4(F10.2,'% '))
 
      end if


*     ***********************************************
*     **** generate projected DENSITY OF STATES *****
*     ***********************************************
      if (control_DOS()) then

      value = MA_push_get(mt_dbl,(nemax),'weight',weight(2),weight(1))
      if (.not. value)
     >  call errquit('analysis:out of stack memory',0, MA_ERR)
      call dcopy(nemax,1.0d0,0,dbl_mb(weight(1)),1)


      if (.not.rtdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
        alpha = 0.05d0/27.2116d0
      end if

      if (.not.rtdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
        npoints = 500
      end if

      if (.not.rtdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
         emin = 99999.0d0
         do ii=1,ne(1)+ne(2)
           if (eig(ii).lt.emin) emin = eig(ii)
         end do
         emin = emin - 0.1d0
      end if

      if (.not.rtdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
         emax = -99999.0d0
         do ii=1,ne(1)+ne(2)
           if (eig(ii).gt.emax) emax = eig(ii)
         end do
         emax = emax + 0.1d0
      end if

*     **** explicit number of atoms have been requested ****
      fixatoms = .false.
      if (rtdb_ma_get(rtdb, 'nwpw:dos:actlist', ma_type,
     >        nactive_atoms, h_actlist)) then

         if (.not.MA_get_index(h_actlist,l_actlist))
     >      call errquit(
     >       'analysis: ma_get_index failed for actlist',911,
     &       MA_ERR)

           fixatoms = .true.
      end if

      write(6,1800)
      if (.not.fixatoms) write(6,1801)
      if (fixatoms) then
        write(6,1802)
        write(6,1803) (int_mb(l_actlist+j-1),j=1,nactive_atoms)
      end if
 1800 FORMAT(///'=== PROJECTED DENSITY OF STATES ==='/)
 1801 FORMAT('  All atoms were used to determine weights')
 1802 FORMAT('  The following atoms were used to determine weights:')
 1803 FORMAT(2x,8I6)

*     **** angular momentum decomposition *****
      lmax = -1
      do k=1,npsp
        if (lmax.le.int_mb(lorb(1)+k-1)) lmax = int_mb(lorb(1)+k-1)
      end do

      do L=0,lmax

         call dcopy(nemax,0.0d0,0,dbl_mb(weight(1)),1)

         if (.not.fixatoms) then
         DO I=1,nion
           DO SPIN=1,ISPIN
           DO N=N1(SPIN),N2(SPIN)
            if (L.le.int_mb(lorb(1)+ion_katm_qm(I)-1)) then
              L1=L**2+1
              L2=(L+1)**2
              DO LL=L1,L2
                dbl_mb(weight(1)+n-1)=
     >          dbl_mb(weight(1)+n-1)+
     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
              END DO
            end if
           END DO
           END DO
         END DO
         else
          DO j=1,nactive_atoms
           I=int_mb(l_actlist+j-1)
           DO SPIN=1,ISPIN
           DO N=N1(SPIN),N2(SPIN)
            if (L.le.int_mb(lorb(1)+ion_katm_qm(I)-1)) then
              L1=L**2+1
              L2=(L+1)**2
              DO LL=L1,L2
                dbl_mb(weight(1)+n-1)=
     >          dbl_mb(weight(1)+n-1)+
     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
              END DO
            end if
           END DO
           END DO
          END DO
         end if

         if (ispin.eq.1) then
           if (flag.eq.0) filename = "DOS_both_"//spdf_name(l)
           if (flag.eq.1) filename = "VDOS_both_"//spdf_name(l)
           call densityofstates(filename,
     >                     eig,dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
         end if
         if (ispin.eq.2) then
           if (flag.eq.0) filename = "DOS_alpha_"//spdf_name(l)
           if (flag.eq.1) filename = "VDOS_alpha_"//spdf_name(l)
           call densityofstates(filename,
     >                     eig,dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
           if (flag.eq.0) filename = "DOS_beta_"//spdf_name(l)
           if (flag.eq.1) filename = "VDOS_beta_"//spdf_name(l)
           call densityofstates(filename,
     >               eig(1+ne(1)),dbl_mb(weight(1)+ne(1)),ne(2),
     >               -1.0d0,alpha,npoints,emin,emax)
         end if
      end do

*     **** combined angular momentum decomposition *****
         call dcopy(nemax,0.0d0,0,dbl_mb(weight(1)),1)
         if (.not.fixatoms) then
         DO I=1,nion
           DO SPIN=1,ISPIN
           DO N=N1(SPIN),N2(SPIN)
            DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
              L1=L**2+1
              L2=(L+1)**2
              DO LL=L1,L2
                dbl_mb(weight(1)+n-1)=
     >          dbl_mb(weight(1)+n-1)+
     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
              END DO
            END DO
           END DO
           END DO
         END DO
         else
          DO j=1,nactive_atoms
           I=int_mb(l_actlist+j-1)
           DO SPIN=1,ISPIN
           DO N=N1(SPIN),N2(SPIN)
            DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
              L1=L**2+1
              L2=(L+1)**2
              DO LL=L1,L2
                dbl_mb(weight(1)+n-1)=
     >          dbl_mb(weight(1)+n-1)+
     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
              END DO
            END DO
           END DO
           END DO
          END DO
         end if

         if (ispin.eq.1) then
           if (flag.eq.0) filename = "DOS_both_all"
           if (flag.eq.1) filename = "VDOS_both_all"
           call densityofstates(filename,
     >                     eig,dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
         end if
         if (ispin.eq.2) then
           if (flag.eq.0) filename = "DOS_alpha_all"
           if (flag.eq.1) filename = "VDOS_alpha_all"
           call densityofstates(filename,
     >                     eig,dbl_mb(weight(1)),ne(1),
     >                     1.0d0,alpha,npoints,emin,emax)
           if (flag.eq.0) filename = "DOS_beta_all"
           if (flag.eq.1) filename = "VDOS_beta_all"
           call densityofstates(filename,
     >               eig(1+ne(1)),dbl_mb(weight(1)+ne(1)),ne(2),
     >               -1.0d0,alpha,npoints,emin,emax)
         end if


*     *** free heap ***
      if(fixatoms) then
        if (.not. MA_free_heap(h_actlist))
     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)
      end if

c*     **** atom decomposition *****
c      do I=1,nion
c
c         call dcopy(nemax,0.0d0,0,dbl_mb(weight(1)),1)
c           DO SPIN=1,ISPIN
c           DO N=N1(SPIN),N2(SPIN)
c            DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
c              L1=L**2+1
c              L2=(L+1)**2
c              DO LL=L1,L2
c                dbl_mb(weight(1)+n-1)=
c     >          dbl_mb(weight(1)+n-1)+
c     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
c              END DO
c            END DO
c           END DO
c           END DO
c
c         if (ispin.eq.1) then
c           write(filename,1801) I
c           call densityofstates(filename,
c     >                     eig,dbl_mb(weight(1)),ne(1),
c     >                     1.0d0,alpha,npoints,emin,emax)
c         end if
c         if (ispin.eq.2) then
c           write(filename,1802) I
c           call densityofstates(filename,
c     >                     eig,dbl_mb(weight(1)),ne(1),
c     >                     1.0d0,alpha,npoints,emin,emax)
c           write(filename,1803) I
c           call densityofstates(filename,
c     >               eig(1+ne(1)),dbl_mb(weight(1)+ne(1)),ne(2),
c     >               -1.0d0,alpha,npoints,emin,emax)
c         end if
c 1801 FORMAT('DOS_both_atom_',I4.4)
c 1802 FORMAT('DOS_alpha_atom_',I4.4)
c 1803 FORMAT('DOS_beta_atom_',I4.4)
c      end do

      value = MA_pop_stack(weight(2))
      if (.not. value)
     >  call errquit('analysis: error freeing stack',0, MA_ERR)

      end if !*** control_DOS ***





*     **** free heap space ****
 1901 continue
      if (.not.mulliken_kawai) call aorbs_end()
      value = MA_free_heap(lorb(2))
      value = MA_free_heap(b0(2))
      value = MA_free_heap(a(2))
      value = MA_free_heap(total(2))
      value = MA_free_heap(sum(2))
      value = MA_free_heap(subtl(2))

      return
      end
 

