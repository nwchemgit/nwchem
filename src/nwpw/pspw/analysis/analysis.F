      subroutine pspw_analysis(rtdb,ispin,ne,psi,eig)
      implicit none
      integer rtdb
      integer ispin,ne(2)
      complex*16 psi(*)
      real*8    eig(*)

#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "rtdb.fh"



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
      integer i,k,l,n,ms,l1,l2
      integer ll,spin,ind
      real*8  ttl1,ttl2,subttl

      logical value
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
      character*2 psp_atom
      external    psp_atom
      integer  ion_nion,psp_npsp,ion_katm
      external ion_nion,psp_npsp,ion_katm
      real*8   ion_rion
      external ion_rion
      real*8      lattice_omega,lattice_ecut,lattice_unita
      external    lattice_omega,lattice_ecut,lattice_unita


      call Parallel_taskid(taskid)
      call Pack_npack(1,npack1)

      npsp = psp_npsp()
      nion = ion_nion()
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

*     ****************************************
*     **** read in expansion coefficients ****
*     ****************************************
      do k=1,npsp
          id = 'analysis:lorb'//psp_atom(k)
          if (.not. rtdb_get(rtdb,id,mt_int,1,int_mb(lorb(1)+k-1))) then
            DD = '                          '
            DD = psp_atom(k)
            ind = index(DD,' ') - 1
            test = DD(1:ind)//'.psp'

            write(*,*) "test:",test,ind
            value = rtdb_parallel(.false.)
            if (taskid.eq.MASTER) then
            ind = index(test,' ') - 1
            write(*,*) "test:",test,ind
            call analysis_expansion_coef(test,-1,rtdb)
            end if
            value = rtdb_parallel(.true.)
            call ga_sync()

            if (.not. rtdb_get(rtdb,id,mt_int,1,int_mb(lorb(1)+k-1))) 
     >        call errquit(
     >        'analysis: rtdb_get lorb failed', 0)
          end if

          id = 'analysis:expansion'//psp_atom(k)
          if (.not. rtdb_get(rtdb,id,mt_dbl,(int_mb(lorb(1)+k-1)+1),
     >                                  dbl_mb(b0(1)+(k-1)*6))) then
            DD = '                          '
            DD = psp_atom(k)
            ind = index(DD,' ') -1
            test = DD(1:ind)//'.psp'
            call analysis_expansion_coef(test,-1,rtdb)

            if (.not. rtdb_get(rtdb,id,mt_dbl,(int_mb(lorb(1)+k-1)+1),
     >                                       dbl_mb(b0(1)+(k-1)*6))) 
     >       call errquit(
     >       'analysis: rtdb_get failed', 0)
          end if
      end do



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
         WRITE(6,1300)
         WRITE(6,1305) 'ATOM ','S','P','D','F'
         DO k=1,npsp
           WRITE(6,1310) psp_atom(K),
     >                  (dbl_mb(b0(1)+l+(k-1)*6), 
     >                   l=0,int_mb(lorb(1)+k-1))
         END DO
         call util_flush(6)
 1300    FORMAT(//'== EXPANSION COEFFICIENTS =='/)
 1305    FORMAT(A5,6X,A,12X,A,12X,A,12X,A)
 1310    FORMAT(A2,' : ',4E13.5)
      end if
 

 
      call util_file_name('ORBOUT',
     >                     .true.,
     >                     .false.,
     >                      id)
      if (taskid.eq.MASTER) 
     > OPEN(UNIT=65,FILE=id,FORM='FORMATTED')

      call Orb_Analysis(65,ispin,ne,npack1,nemax,psi,
     >                        int_mb(lorb(1)),
     >                        dbl_mb(b0(1)),
     >                        dbl_mb(a(1)),
     >                        dbl_mb(sum(1)))
     
      if (taskid.eq.MASTER) close(unit=65)



      if (taskid.eq.MASTER) then
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '==============================================='
      WRITE(6,*) '|  POPULATION ANALYSIS OF MOLECULAR ORBITALS  |'
      WRITE(6,*) '==============================================='
 
      DO SPIN=1,ISPIN
        DO N=N1(SPIN),N2(SPIN)
          WRITE(6,1500)
          IF(ISPIN.EQ.2) THEN
            WRITE(6,1510) N,SPN(SPIN),dbl_mb(sum(1)+N-1),eig(n)
          ELSE
            WRITE(6,1515) N,dbl_mb(SUM(1)+N-1),eig(n)
          ENDIF
          WRITE(6,1520) 'NO','ATOM','L','POPULATION'

          DO L=0,5
            dbl_mb(SUBTL(1)+L)=0.0d0
          END DO
          DO I=1,nion
            dbl_mb(TOTAL(1)+I-1)=0.0d0
            DO L=0,int_mb(lorb(1)+ion_katm(I)-1)
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
              WRITE(6,1530) I,psp_atom(ion_katm(I)),L,SUBTTL,
     >               (dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax),
     >                        LL=L1,L2)

            END DO
          END DO

          WRITE(6,1540)
          WRITE(6,1550) (I,psp_atom(ion_katm(I)),
     >                   dbl_mb(TOTAL(1)+I-1),I=1,nion)
          WRITE(6,1555)
          WRITE(6,1560)  's','p','d','f'
          WRITE(6,1570)  (dbl_mb(SUBTL(1)+L),L=0,3)
        END DO
      END DO
      call util_flush(6)

 1500 FORMAT(//'------------------------------------------------'//)
 1510 FORMAT('*** ORBITAL=',I2,'***  SPIN=',A4,4X,'SUM=',E12.5,
     >       ' E=',E12.5/)
 1515 FORMAT('*** ORBITAL=',I2,'***  SPIN=BOTH',4X,'SUM=',E12.5,
     >       ' E=',E12.5/)
 1520 FORMAT(A2,2X,A4,2X,A1,5X,A10)
 1530 FORMAT(I2,3X,A2,3X,I1,8F11.5)
 1531 FORMAT(8F11.5)
 1540 FORMAT(//'=== DISTRIBUTION ==='/)
 1550 FORMAT(4(I6,'(',A2,')',F9.4))
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

        DO L=0,int_mb(lorb(1)+ion_katm(I)-1)
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
        DO L=0,int_mb(lorb(1)+ion_katm(I)-1)
          TTL1=TTL1+dbl_mb(SUBTL(1)+L)
        END DO
        WRITE(6,1620) I,psp_atom(ion_katm(I)),SPN(1),TTL1,
     >                 ( dbl_mb(SUBTL(1)+L),
     >                   L=0,int_mb(lorb(1)+ion_katm(I)-1) )
        TTL1=0.0d0
        DO L=0,int_mb(lorb(1)+ion_katm(I)-1)
          TTL1=TTL1+dbl_mb(SUBTL(1)+L+(ispin-1)*6)
        END DO
        WRITE(6,1620) I,psp_atom(ion_katm(I)),SPN(2),TTL1,
     >                 ( dbl_mb(SUBTL(1)+L+(ispin-1)*6),
     >                   L=0,int_mb(lorb(1)+ion_katm(I)-1) )

      END DO
      call util_flush(6)

 1610 FORMAT(A2,2X,A4,2X,A4,4X,A5,7X,A,10X,A,10X,A,10X,A)
 1620 FORMAT(I2,3X,A2,3X,A4,5F11.5)
 
      DO L=0,3
        dbl_mb(SUBTL(1)+L)  =0.0d0
        dbl_mb(SUBTL(1)+L+6)=0.0d0
      END DO
      DO I=1,nion
        DO SPIN=1,ISPIN
          DO N=N1(SPIN),N2(SPIN)
            DO L=0,int_mb(lorb(1)+ion_katm(I)-1)
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
        IF(ISPIN.EQ.2) dbl_mb(SUBTL(1)+L+6)=dbl_mb(SUBTL(1)+L+6)
     >                                      *100.0d0/dble(NE(2))
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


*     **** free heap space ****
      value = MA_free_heap(lorb(2))
      value = MA_free_heap(b0(2))
      value = MA_free_heap(a(2))
      value = MA_free_heap(total(2))
      value = MA_free_heap(sum(2))
      value = MA_free_heap(subtl(2))

      return
      end
 

