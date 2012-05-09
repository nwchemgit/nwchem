c
c $Id: band_analysis.F 22397 2012-05-07 23:08:15Z bylaska $
c

      subroutine band_analysis(flag,rtdb,ispin,ne,
     >                         nbrill,nbrillq,psi_tag,eig_tag)
      implicit none
      integer flag
      integer rtdb
      integer ispin,ne(2),nbrill,nbrillq
      integer psi_tag
      integer eig_tag

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"

*     **** parallel variables ****
      integer  taskid,taskid_k
      integer  MASTER
      parameter(MASTER=0)


*     **** electronic variables ****
      integer npack1
      integer n1(2),n2(2)
    
*     **** local variables ****
      integer i,k,l,n,ms,l1,l2,j,nbq,nb,pk
      integer ll,spin,ind
      real*8  ttl1,ttl2,subttl

      integer h_actlist,l_actlist,nactive_atoms,ma_type
      integer npoints,ii,lmax
      integer weight(2)
      real*8 emin,emax,alpha,lmbda,rcut
      character*255 filename

      logical value,fixatoms
      character*28 DD
      character*255 id,test
      integer npsp,nion,nemax
      integer lorb(2)    ! integer lorb(npsp)
      integer b0(2)      ! real*8 b0(0:5,npsp)
      integer total(2)   ! real*8 total(nion)
      integer subtl(2)   ! real*8 subtl(0:5,3)
      integer a(2)       ! real*8 a(36,nemax,nion)
      integer sum(2)     ! real*8 sum(nemax)
      integer eig(2)     ! real*8 sum(nemax)
      integer A_tag,sum_tag
      integer A_shift,sum_shift,psi_shift,eig_shift
 
      character*4 spn(2)
      DATA SPN / 'UP  ', 'DOWN' /


*     **** external functions ****
      logical  control_DOS,nwpw_filefind
      logical  borbs_init,borbs_readall
      external control_DOS,nwpw_filefind
      external borbs_init,borbs_readall
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
      real*8   borbs_rcut,borbs_lmbda
      external borbs_rcut,borbs_lmbda
      integer  cpsi_data_alloc,cpsi_data_get_chnk
      external cpsi_data_alloc,cpsi_data_get_chnk


      call Parallel_taskid(taskid)
      call Parallel3d_taskid_k(taskid_k)
      call Cram_max_npack(npack1)

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
     >        MA_alloc_get(mt_dbl,nion,'total',total(2),total(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,6*3,'subtl',subtl(2),subtl(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,36*nemax*nion,'a',a(2),a(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nemax,'sum',sum(2),sum(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nemax,'eig',eig(2),eig(1))

      A_tag   = cpsi_data_alloc(nbrillq,36*nemax*nion,1)
      sum_tag = cpsi_data_alloc(nbrillq,nemax,1)


*      ***********************************************
*      **** create psp1 files if they don't exist ****
*      ***********************************************
       do k=1,npsp
         DD = '                          '
         DD = ion_atom_qm(k)
         ind = index(DD,' ') - 1
         test = DD(1:ind)//'.psp1'
         id   = DD(1:ind)//'.borb'
         call control_mullikenparameters(ion_atom_qm(k),rcut,lmbda)
         if ((.not.nwpw_filefind(test)).or.(.not.nwpw_filefind(id)))
     >      call borbs_formatter_auto(ion_atom_qm(k),rcut,lmbda)
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


      value = borbs_init()
      value = value.and.borbs_readall()
      if (.not.value) go to 1901

      value = .true.
      do k=1,npsp
         call control_mullikenparameters(ion_atom_qm(k),rcut,lmbda)
         if ((dabs(borbs_rcut(k)-rcut).gt.1.0d-6).or.
     >       (dabs(borbs_lmbda(k)-lmbda).gt.1.0d-6)) then
            call borbs_formatter_auto(ion_atom_qm(k),rcut,lmbda)
            value = .false.
         end if
      end do
      if (.not.value) value = borbs_readall()
      if (.not.value) go to 1901


      if (taskid.eq.MASTER) then
         call util_date(DD)
 
         WRITE(luout,*)
         WRITE(luout,*)
         WRITE(luout,*)
         WRITE(luout,*) 
     >  '*************************************************************'
         WRITE(luout,*) 
     >  '**                                                         **'
         WRITE(luout,*) 
     >  '**                  BAND Mulliken analysis                 **'
         WRITE(luout,*) 
     >  '**                                                         **'
         if (flag.eq.1)
     >   WRITE(luout,*)
     >  '**                   (Virtual Orbitals)                    **'
         WRITE(luout,*) 
     >  '** Population analysis algorithm developed by Ryoichi Kawai**'
         WRITE(luout,*) 
     >  '**                                                         **'
         WRITE(luout,1000) DD
 1000    FORMAT(
     >  ' **                   ',A16,'                      **')
         WRITE(luout,*)
     >   '**                                                         **'
         WRITE(luout,*) 
     >  '*************************************************************'
      end if


c     **** ouput xyz format ****
      call ion_Print_XYZ(luout)


      value = rtdb_parallel(.false.)
      if (taskid.eq.MASTER) then
         write(luout,1307)
         do k=1,npsp
           call control_mullikenparameters(ion_atom_qm(k),rcut,lmbda)
           if (lmbda.gt.0.0d0) then
              write(luout,1308) ion_atom_qm(k),"damping",rcut,lmbda
           else
              write(luout,*) ion_atom_qm(k)," nodamping"
           end if
         end do
         call util_flush(luout)
 1307    FORMAT(//'== Atomic Orbital Expansion =='/)
 1308    FORMAT(A5,A10,4x,"rcut=",F8.3,2x,"lmbda=",F8.3)
      end if
      value = rtdb_parallel(.true.)
 

 
c      call util_file_name('ORBOUT',
c     >                     .true.,
c     >                     .false.,
c     >                      id)
c      if (taskid.eq.MASTER) 
c     > OPEN(UNIT=65,FILE=id,FORM='FORMATTED')
c      call Orb_Analysis(65,flag,ispin,ne,npack1,nemax,psi,
c     >                        int_mb(lorb(1)),
c     >                        dbl_mb(b0(1)),
c     >                        dbl_mb(a(1)),
c     >                        dbl_mb(sum(1)))
c      if (taskid.eq.MASTER) close(unit=65)

      do nbq=1,nbrillq
        psi_shift = cpsi_data_get_chnk(psi_tag,nbq)
        A_shift   = cpsi_data_get_chnk(A_tag,nbq)
        sum_shift = cpsi_data_get_chnk(sum_tag,nbq)

        call Orb_pop_borb(flag,nbq,ispin,ne,
     >                  npack1,nemax,dbl_mb(psi_shift),
     >                  int_mb(lorb(1)),
     >                  dbl_mb(A_shift),dbl_mb(sum_shift))
      end do


      if (taskid.eq.MASTER) then
      WRITE(luout,*)
      WRITE(luout,*)
      WRITE(luout,*) 
     > '====================================================='
      if (flag.eq.0)
     >WRITE(luout,*) 
     > '| POPULATION ANALYSIS OF FILLED MOLECULAR ORBITALS  |'
      if (flag.eq.1)
     >WRITE(luout,*) 
     > '| POPULATION ANALYSIS OF VIRTUAL MOLECULAR ORBITALS |'
      WRITE(luout,*) 
     > '====================================================='
      WRITE(luout,1311)
c     WRITE(6,1313)
c     WRITE(6,1314)
c     WRITE(6,1315)
 1311 FORMAT(//'== Using pseudoatomic orbital expansion          ==')
      end if

 
      do nb=1,nbrill
         call dcopy(nemax,0.0d0,0,dbl_mb(eig(1)),1)
         call dcopy(nemax,0.0d0,0,dbl_mb(sum(1)),1)
         call dcopy(36*nion*nemax,0.0d0,0,dbl_mb(a(1)),1)

         call K1dB_ktoqp(nb,nbq,pk)
         if (pk.eq.taskid_k) then
            eig_shift = cpsi_data_get_chnk(eig_tag,nbq)
            A_shift   = cpsi_data_get_chnk(A_tag,nbq)
            sum_shift = cpsi_data_get_chnk(sum_tag,nbq)
            call dcopy(nemax,dbl_mb(eig_shift),1,dbl_mb(eig(1)),1)
            call dcopy(nemax,dbl_mb(sum_shift),1,dbl_mb(sum(1)),1)
            call dcopy(36*nion*nemax,dbl_mb(A_shift),1,dbl_mb(a(1)),1)
         end if
         call K1dB_Vector_SumAll(nemax,dbl_mb(eig(1)))
         call K1dB_Vector_SumAll(nemax,dbl_mb(sum(1)))
         call K1dB_Vector_SumAll(36*nion*nemax,dbl_mb(a(1)))

         if (taskid.eq.MASTER) then
         DO SPIN=1,ISPIN
         DO N=N1(SPIN),N2(SPIN)
          WRITE(luout,1500)
          IF(ISPIN.EQ.2) THEN
            WRITE(luout,1510) nb,N,SPN(SPIN),dbl_mb(sum_shift+N-1),
     >                    dbl_mb(eig(1)+n-1),
     >                    dbl_mb(eig(1)+n-1)*27.2116d0
          ELSE
            WRITE(luout,1515) nb,N,dbl_mb(sum(1)+N-1),
     >                   dbl_mb(eig(1)+n-1),
     >                   dbl_mb(eig(1)+n-1)*27.2116d0
          ENDIF
          !write(6,1519)
          WRITE(luout,1520) 'NO','ATOM','L','POPULATION'

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
     >                dbl_mb(a(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
                SUBTTL=SUBTTL+
     >                dbl_mb(a(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
              END DO
              dbl_mb(SUBTL(1)+L)=dbl_mb(SUBTL(1)+L)+SUBTTL
              if (l.eq.0) write(luout,1516)
              if (l.eq.1) write(luout,1517)
              if (l.eq.2) write(luout,1518)
              if (l.eq.3) write(luout,1519)
              WRITE(luout,1530) I,ion_atom_qm(ion_katm_qm(I)),L,SUBTTL,
     >               (dbl_mb(a(1)+LL-1+(N-1)*36+(I-1)*36*nemax),
     >                        LL=L1,L2)

            END DO
          END DO

          WRITE(luout,1540)
          WRITE(luout,1550) (I,ion_atom_qm(ion_katm_qm(I)),
     >                   dbl_mb(TOTAL(1)+I-1),I=1,nion)
          WRITE(luout,1555)
          WRITE(luout,1560)  's','p','d','f'
          WRITE(luout,1570)  (dbl_mb(SUBTL(1)+L),L=0,3)
        END DO
        END DO
        call util_flush(luout)
        end if
      end do

 1500 FORMAT(//'------------------------------------------------',
     >         '------------------------------'//)
 1510 FORMAT('*** NB=',I4,' ORBITAL=',I4,'***  SPIN=',A4,
     >       4X,'SUM=',E12.5,
     >       ' E=',E12.5,' (',F8.3,'eV)'/)
 1515 FORMAT('*** NB=',I4,' ORBITAL=',I4,'***  SPIN=BOTH',
     >       4X,'SUM=',E12.5,
     >       ' E=',E12.5,' (',F8.3,'eV)'/)
 1516 FORMAT(27x,' s')
 1517 FORMAT(27x,' px         pz         py')
 1518 FORMAT(27x,' dx2-y2     dzx        d3z2-1     dyz        dxy')
 1519 FORMAT(27x,
     > ' fx(x2-3y2) fz(5z2-1)  fx(5z2-1)  fz(5z2-3)  fy(5z2-1) ',
     > ' fxyz       fy(3x2-y2)')

c 1519 FORMAT(30x,'  s',
c     >      /30x,' px         py        pz '
c     >      /27x,'d3z2-1        dxy        dyz        dzx     dx2-y2',
c     >      /30x,' fy(3x2-y2) fxyz      fy(5z2-1)   fz(5z2-3)  fx(5z2-1)', 
c     >           '  fz(5z2-1)  fx(x2-3y2)')

 1520 FORMAT(A2,2X,A4,2X,A1,2X,A10)
 1530 FORMAT(I2,3X,A4,3X,I1,8F11.5)
 1531 FORMAT(8F11.5)
 1540 FORMAT(//'=== DISTRIBUTION ==='/)
 1550 FORMAT(4(I6,'(',A4,')',F9.4))
 1555 FORMAT(//'== ANGULAR MOMENTUM POPULATIONS ==='/)
 1560 FORMAT(6X,A1,3(9X,A1))
 1570 FORMAT(4F10.4)
 
      if (taskid.eq.MASTER) then
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '========================================'
      WRITE(6,*) '|   POPULATION ANALYSIS ON EACH ATOM   |'
      WRITE(6,*) '========================================'
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,1610) 'NO','ATOM','SPIN','TOTAL','s','p','d','f'
      end if

c      DO I=1,nion
c        TTL1=0.0d0
c        TTL2=0.0d0
c
c        DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
c          DO ms=1,ISPIN
c            dbl_mb(SUBTL(1)+L+(ms-1)*6)=0.0d0
c            DO N=N1(ms),N2(ms)
c              L1=L**2+1
c              L2=(L+1)**2
c              DO LL=L1,L2
c                dbl_mb(SUBTL(1)+L+(ms-1)*6)=
c     >            dbl_mb(SUBTL(1)+L+(ms-1)*6) +
c     >            dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
c
c              END DO
c            END DO
c          END DO
c        END DO
c     
c        TTL1=0.0d0
c        DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
c          TTL1=TTL1+dbl_mb(SUBTL(1)+L)
c        END DO
c        WRITE(luout,1620) I,ion_atom_qm(ion_katm_qm(I)),SPN(1),TTL1,
c     >                 ( dbl_mb(SUBTL(1)+L),
c     >                   L=0,int_mb(lorb(1)+ion_katm_qm(I)-1) )
c        TTL1=0.0d0
c        DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
c          TTL1=TTL1+dbl_mb(SUBTL(1)+L+(ispin-1)*6)
c        END DO
c        WRITE(luout,1620) I,ion_atom_qm(ion_katm_qm(I)),SPN(2),TTL1,
c     >                 ( dbl_mb(SUBTL(1)+L+(ispin-1)*6),
c     >                   L=0,int_mb(lorb(1)+ion_katm_qm(I)-1) )
c
c      END DO
c      call util_flush(luout)

 1610 FORMAT(A2,2X,A4,2X,A4,4X,A5,7X,A,10X,A,10X,A,10X,A)
 1620 FORMAT(I2,3X,A4,3X,A4,5F11.5)
 
c      DO L=0,3
c        dbl_mb(SUBTL(1)+L)  =0.0d0
c        dbl_mb(SUBTL(1)+L+6)=0.0d0
c      END DO
c      DO I=1,nion
c        DO SPIN=1,ISPIN
c          DO N=N1(SPIN),N2(SPIN)
c            DO L=0,int_mb(lorb(1)+ion_katm_qm(I)-1)
c              L1=L**2+1
c              L2=(L+1)**2
c              DO LL=L1,L2
c                dbl_mb(SUBTL(1)+L+(SPIN-1)*6)=
c     >          dbl_mb(SUBTL(1)+L+(SPIN-1)*6)+
c     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
c              END DO
c            END DO
c          END DO
c        END DO
c      END DO
c
c      DO L=0,3
c        dbl_mb(SUBTL(1)+L+2*6)=
c     >  (dbl_mb(SUBTL(1)+L)+dbl_mb(SUBTL(1)+L+(ISPIN-1)*6))
c     >   *100.d0/(NE(1)+NE(ISPIN))
c        dbl_mb(SUBTL(1)+L)=dbl_mb(SUBTL(1)+L)*100.0d0/dble(NE(1))
c        IF((ISPIN.EQ.2).and.(NE(2).gt.0)) 
c     >    dbl_mb(SUBTL(1)+L+6)=dbl_mb(SUBTL(1)+L+6)*100.0d0/dble(NE(2))
c      END DO

      if (taskid.eq.MASTER) then
      WRITE(luout,1700)
      WRITE(luout,1710) ' SPIN ','s','p','d','f'
      WRITE(luout,1720) SPN(1),(dbl_mb(SUBTL(1)+L),L=0,3)
      WRITE(luout,1720) SPN(ISPIN),
     >                  (dbl_mb(SUBTL(1)+L+(ISPIN-1)*6),L=0,3)
      WRITE(luout,1720) ' TOTAL',(dbl_mb(SUBTL(1)+L+(3-1)*6),L=0,3)
      call util_flush(luout)
 1700 FORMAT(///'=== TOTAL ANGULAR MOMENTUM POPULATION ==='/)
 1710 FORMAT(A6,6X,A1,3(11X,A1))
 1720 FORMAT(A6,4(F10.2,'% '))
      end if
 


*     ***********************************************
*     **** generate projected DENSITY OF STATES *****
*     ***********************************************
      if (control_DOS()) then

c      value = MA_push_get(mt_dbl,(nemax),'weight',weight(2),weight(1))
c      if (.not. value)
c     >  call errquit('analysis:out of stack memory',0, MA_ERR)
c      call dcopy(nemax,1.0d0,0,dbl_mb(weight(1)),1)
c
c
c      if (.not.rtdb_get(rtdb,'dos:alpha',mt_dbl,1,alpha)) then
c        alpha = 0.05d0/27.2116d0
c      end if
c
c      if (.not.rtdb_get(rtdb,'dos:npoints',mt_int,1,npoints)) then
c        npoints = 500
c      end if
c
c      if (.not.rtdb_get(rtdb,'dos:emin',mt_dbl,1,emin)) then
c         emin = 99999.0d0
c         do ii=1,ne(1)+ne(2)
c           if (eig(ii).lt.emin) emin = eig(ii)
c         end do
c         emin = emin - 0.1d0
c      end if
c
c      if (.not.rtdb_get(rtdb,'dos:emax',mt_dbl,1,emax)) then
c         emax = -99999.0d0
c         do ii=1,ne(1)+ne(2)
c           if (eig(ii).gt.emax) emax = eig(ii)
c         end do
c         emax = emax + 0.1d0
c      end if
c
c*     **** explicit number of atoms have been requested ****
c      fixatoms = .false.
c      if (rtdb_ma_get(rtdb, 'nwpw:dos:actlist', ma_type,
c     >        nactive_atoms, h_actlist)) then
c
c         if (.not.MA_get_index(h_actlist,l_actlist))
c     >      call errquit(
c     >       'analysis: ma_get_index failed for actlist',911,
c     &       MA_ERR)
c
c           fixatoms = .true.
c      end if
c
c      write(6,1800)
c      if (.not.fixatoms) write(6,1801)
c      if (fixatoms) then
c        write(6,1802)
c        write(6,1803) (int_mb(l_actlist+j-1),j=1,nactive_atoms)
c      end if
 1800 FORMAT(///'=== PROJECTED DENSITY OF STATES ==='/)
 1801 FORMAT('  All atoms were used to determine weights')
 1802 FORMAT('  The following atoms were used to determine weights:')
 1803 FORMAT(2x,8I6)

c*     **** angular momentum decomposition *****
c      lmax = -1
c      do k=1,npsp
c        if (lmax.le.int_mb(lorb(1)+k-1)) lmax = int_mb(lorb(1)+k-1)
c      end do
c
c      do L=0,lmax
c
c         call dcopy(nemax,0.0d0,0,dbl_mb(weight(1)),1)
c
c         if (.not.fixatoms) then
c         DO I=1,nion
c           DO SPIN=1,ISPIN
c           DO N=N1(SPIN),N2(SPIN)
c            if (L.le.int_mb(lorb(1)+ion_katm_qm(I)-1)) then
c              L1=L**2+1
c              L2=(L+1)**2
c              DO LL=L1,L2
c                dbl_mb(weight(1)+n-1)=
c     >          dbl_mb(weight(1)+n-1)+
c     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
c              END DO
c            end if
c           END DO
c           END DO
c         END DO
c         else
c          DO j=1,nactive_atoms
c           I=int_mb(l_actlist+j-1)
c           DO SPIN=1,ISPIN
c           DO N=N1(SPIN),N2(SPIN)
c            if (L.le.int_mb(lorb(1)+ion_katm_qm(I)-1)) then
c              L1=L**2+1
c              L2=(L+1)**2
c              DO LL=L1,L2
c                dbl_mb(weight(1)+n-1)=
c     >          dbl_mb(weight(1)+n-1)+
c     >          dbl_mb(A(1)+LL-1+(N-1)*36+(I-1)*36*nemax)**2
c              END DO
c            end if
c           END DO
c           END DO
c          END DO
c         end if
c
c         if (ispin.eq.1) then
c           if (flag.eq.0) filename = "DOS_both_"//spdf_name(l)
c           if (flag.eq.1) filename = "VDOS_both_"//spdf_name(l)
c           call densityofstates(filename,
c     >                     eig,dbl_mb(weight(1)),ne(1),
c     >                     1.0d0,alpha,npoints,emin,emax)
c         end if
c         if (ispin.eq.2) then
c           if (flag.eq.0) filename = "DOS_alpha_"//spdf_name(l)
c           if (flag.eq.1) filename = "VDOS_alpha_"//spdf_name(l)
c           call densityofstates(filename,
c     >                     eig,dbl_mb(weight(1)),ne(1),
c     >                     1.0d0,alpha,npoints,emin,emax)
c           if (flag.eq.0) filename = "DOS_beta_"//spdf_name(l)
c           if (flag.eq.1) filename = "VDOS_beta_"//spdf_name(l)
c           call densityofstates(filename,
c     >               eig(1+ne(1)),dbl_mb(weight(1)+ne(1)),ne(2),
c     >               -1.0d0,alpha,npoints,emin,emax)
c         end if
c      end do
c
c*     **** combined angular momentum decomposition *****
c         call dcopy(nemax,0.0d0,0,dbl_mb(weight(1)),1)
c         if (.not.fixatoms) then
c         DO I=1,nion
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
c         END DO
c         else
c          DO j=1,nactive_atoms
c           I=int_mb(l_actlist+j-1)
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
c          END DO
c         end if
c
c         if (ispin.eq.1) then
c           if (flag.eq.0) filename = "DOS_both_all"
c           if (flag.eq.1) filename = "VDOS_both_all"
c           call densityofstates(filename,
c     >                     eig,dbl_mb(weight(1)),ne(1),
c     >                     1.0d0,alpha,npoints,emin,emax)
c         end if
c         if (ispin.eq.2) then
c           if (flag.eq.0) filename = "DOS_alpha_all"
c           if (flag.eq.1) filename = "VDOS_alpha_all"
c           call densityofstates(filename,
c     >                     eig,dbl_mb(weight(1)),ne(1),
c     >                     1.0d0,alpha,npoints,emin,emax)
c           if (flag.eq.0) filename = "DOS_beta_all"
c           if (flag.eq.1) filename = "VDOS_beta_all"
c           call densityofstates(filename,
c     >               eig(1+ne(1)),dbl_mb(weight(1)+ne(1)),ne(2),
c     >               -1.0d0,alpha,npoints,emin,emax)
c         end if
c
c
c*     *** free heap ***
c      if(fixatoms) then
c        if (.not. MA_free_heap(h_actlist))
c     >   call errquit('h_actlist:error freeing heap memory',0, MA_ERR)
c      end if
c
c      value = MA_pop_stack(weight(2))
c      if (.not. value)
c     >  call errquit('analysis: error freeing stack',0, MA_ERR)

      end if !*** control_DOS ***




*     **** free heap space ****
 1901 continue
      call borbs_end()
      call cpsi_data_dealloc(A_tag)
      call cpsi_data_dealloc(sum_tag)
      value = MA_free_heap(a(2))
      value = MA_free_heap(sum(2))
      value = MA_free_heap(eig(2))
      value = MA_free_heap(lorb(2))
      value = MA_free_heap(b0(2))
      value = MA_free_heap(total(2))
      value = MA_free_heap(subtl(2))

      return
      end
 

