      subroutine ducc(rtdb,d_t1,k_t1_offset,d_t2,k_t2_offset,
     &    d_v2,k_v2_offset,ene_orb,nos,nvs,noas,nobs,nvas,nvbs,nactv,
     D        transform,actindex1,printindex1,
     D        actindex2,printindex2,
     D        h,horb,v,
     D        t1,t2,
     D        ht,vt)
#ifndef OPTIMIZE_DUCC
CDIR$ OPTIMIZE: 2
#endif
c
c  Only for RHF
c
c all varaible below have to be defined in the tce_energy.f
c the nact should be defined in the input set tce:nactv
c
c All occupied orbitals are treated as active
c nactv defines the number of active virtual orbitals
c please do not freeze the occupied orbitals for now
c
c nactv - number of active virtual orbitals
c nos   - number of occupied spinorbitals
c nvs   - number of virtual spinorbitals
c noas  - number of occupied alpha spinorbitals
c nobs  - number of occupied beta spinorbitals
c nvas  - number of virtual  alhoa spinorbitals
c nvbs  - number of virtual beta spinorbitals
c
c spinorbital convention
c
c | noas | nobs | nvas | nvbs |
c
c ene_orb contains orbital energies
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_t1    ! GA handle for t1
      integer k_t1_offset
      integer d_t2    ! GA handle for t2
      integer k_t2_offset
      integer d_v2    ! GA handle for v2
      integer k_v2_offset
      integer nos,nvs ! # of occupied/virtual spinorbitals
      integer noas,nobs,nvas,nvbs ! # of occupied/virtual alphas betas
      integer nactv
      double precision ene_orb(nos+nvs)  ! a copy of dbl_mb(k_sorted)
      integer size
      integer rtdb
      integer l_aux1,k_aux1 ! for local memory allocator: loc. mem. buffer 1
      integer i,j,k,l         ! auxiliary indices
      integer ia,ib,ic,id,ie,if !auxiliary indices
      integer nact              !number of active virtual orbitals
c nicks transform.
      integer g,m,n,p,q,r,s
      DOUBLE PRECISION PTHRESH  !PRINTING THRESHOLD
      INTEGER TRANSFORM(NOS+NVS)  !RETURNS THE CORRESPONDING ORBITAL LABEL FOR A GIVEN SPIN ORBITAL

      DOUBLE PRECISION ZZZ
C INDEXING ARRAY TO DETERMINE IF A SET OF INDICES CORRESPONDS TO AN
C INTERNAL OR EXTERNAL SET.
      INTEGER ACTINDEX1(NOS+NVS,NOS+NVS)
      INTEGER PRINTINDEX1(NOS+NVS,NOS+NVS)
      logical*1 ACTINDEX2(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)
      logical*1 PRINTINDEX2(NOS+NVS,NOS+NVS,NOS+NVS,NOS+NVS)
c matrices with original integrals
c h created here from "orbital" matrix horb
      double precision h(nos+nvs,nos+nvs)
c horb valid only for RHF case
      double precision horb((nos+nvs)/2,(nos+nvs)/2)
      double precision v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
c t2 in a nice representation : be careful here (you may want to reindex
c virtuals
      double precision t1(nos,nos+1:nos+nvs)
      double precision t2(nos,nos,nos+1:nos+nvs,nos+1:nos+nvs)
c transformed matrices
      double precision ht(nos+nvs,nos+nvs)
!      double precision FOCKt(nos+nvs,nos+nvs)
      double precision vt(nos+nvs,nos+nvs,nos+nvs,nos+nvs)  
c correlation energies
      double precision eccsd
      double precision xxx
c
      logical nodezero
      logical oprint_qa
c
c
      nodezero = (ga_nodeid().eq.0)
      oprint_qa=util_print('ducc_qa', print_high)
c
c one electron integrals (original core Hamiltonian)    
c      call kinetic_hcore_1(rtdb,horb)
c       call kinetic_hcore(rtdb,horb,noas+nvas,spust)
c mapping from horb-->h --------------
c
c
      do i=1,(nos+nvs)/2
      do j=1,(nos+nvs)/2
        horb(i,j)=0.0d0
      enddo
      enddo
c
      do i=1,nos+nvs
      do j=1,nos+nvs
        h(i,j)=0.0d0
        ht(i,j)=0.0d0
      enddo
      enddo
c
      call kinetic_hcore_1(rtdb,horb,h,nos,nvs,noas,nobs,nvas,nvbs)
c
c     
c forming v matrix
c
      do i=1,nos+nvs
      do j=1,nos+nvs
      do k=1,nos+nvs
      do l=1,nos+nvs
        v(i,j,k,l)=0.0d0 
        vt(i,j,k,l)=0.0d0
      enddo
      enddo
      enddo
      enddo
      call mapping_v2_m(rtdb,d_v2,k_v2_offset,v,nos,nvs)
c
      do i=1,nos
      do ia=nos+1,nos+nvs
        t1(i,ia)=0.0d0
      enddo
      enddo
c
      do i=1,nos
        do j=1,nos
          do ia=nos+1,nos+nvs
            do ib=nos+1,nos+nvs
               t2(i,j,ia,ib)=0.0d0
            enddo
          enddo
        enddo
      enddo
c
      call mapping_t1(rtdb,d_t1,k_t1_offset,t1,nos,nvs)
c
c
      call mapping_t2(rtdb,d_t2,k_t2_offset,t2,nos,nvs)
c
c
c check point: calculating CCSD correlation energy 
c 
      eccsd=0.0d0
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        eccsd=eccsd+0.250d0*v(i,j,ia,ib)*t2(i,j,ia,ib)
        eccsd=eccsd+0.50d0*v(i,j,ia,ib)*t1(i,ia)*t1(j,ib)
      enddo
      enddo
      enddo
      enddo
c 
      if(nodezero) then 
        write(6,*)'From DUCC CCSD corr. ene.',eccsd
        call util_flush(6)
      endif
c Attention: t1e and t2e calculated for specific SES g(no,nactv)
c Attention: for effective interactions different types of SESs
c            are needed!
c transformation t1-->t1e  
      do i=1,nos
c      a->a 
       do ia=nos+1,nos+nactv
        t1(i,ia)=0.0d0
       enddo
c      b->b
       do ia=nos+nvas+1,nos+nvas+nactv
        t1(i,ia)=0.0d0
       enddo
      enddo !t1e
c transformation t2-->t2e  ! in QDK only t2e is used
c we are using the same t1 and t2 tensors to store t1e and t2e
      do i=1,nos
      do j=1,nos
c virt a a
        do ia=nos+1,nos+nactv
        do ib=nos+1,nos+nactv
         t2(i,j,ia,ib)=0.0d0
        enddo
        enddo
c virt a b  & b a
        do ia=nos+1,nos+nactv
        do ib=nos+nvas+1,nos+nvas+nactv
         t2(i,j,ia,ib)=0.0d0
         t2(i,j,ib,ia)=0.0d0
        enddo
        enddo
c virt b b 
        do ia=nos+nvas+1,nos+nvas+nactv
        do ib=nos+nvas+1,nos+nvas+nactv
         t2(i,j,ia,ib)=0.0d0
        enddo
        enddo
c
       enddo  !j t2->t2e
       enddo  !i t2->t2e

c  Nick's work starts here 
c  one electron  - only: alpha-alpha & beta-beta
c  two electron  - only: vt(alpha,beta,alpha,beta)
c
c  calculate ht  (print only ht(alpha,alpha)

c  calculate vt  (print only vt(alpha,beta,alpha,beta)

! ***|######## ########     ###    ##    ##  ######  ########  #######  ########  ##     ##
! ***|   ##    ##     ##   ## ##   ###   ## ##    ## ##       ##     ## ##     ## ###   ###
! ***|   ##    ##     ##  ##   ##  ####  ## ##       ##       ##     ## ##     ## #### ####
! ***|   ##    ########  ##     ## ## ## ##  ######  ######   ##     ## ########  ## ### ##
! ***|   ##    ##   ##   ######### ##  ####       ## ##       ##     ## ##   ##   ##     ##
! ***|   ##    ##    ##  ##     ## ##   ### ##    ## ##       ##     ## ##    ##  ##     ##
! ***|   ##    ##     ## ##     ## ##    ##  ######  ##        #######  ##     ## ##     ##

      HT = 0.0d0
      VT = 0.0d0
! START 11111111111111111111111111111111111111111111111111111111111111111111111
! ***| HT_{A}^{B}
      DO ia = NOS+1, NOS+NVS
        DO ib = NOS+1, NOS+NVS

! ***| H_N
          IF(ia.EQ.ib) HT(ia,ib) = HT(ia,ib) + ENE_ORB(ia)

! ***| (H_N(T_EXT))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(ia,ib) = HT(ia,ib) + V(ie,ia,M,ib)*T1(M,ie)
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                HT(ia,ib) = HT(ia,ib) - 0.5*V(ie,ia,M,N)*T2(M,N,ie,ib)
              ENDDO
            ENDDO
          ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(ib,ia) = HT(ib,ia) + V(M,ib,ie,ia)*T1(M,ie)
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                HT(ib,ia) = HT(ib,ia) - 0.5*V(M,N,ie,ia)*T2(M,N,ie,ib)
              ENDDO
            ENDDO
          ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
          DO M = 1, NOS
            DO N = 1, NOS
              IF(M.EQ.N) THEN
                HT(ib,ia) = HT(ib,ia) + 0.5*T1(M,ib)*ENE_ORB(M)*T1(N,ia)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              IF(ib.EQ.ie) THEN
                HT(ib,ia)=HT(ib,ia)-0.5d0*T1(M,ie)*ENE_ORB(ib)*T1(M,ia)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(ia.EQ.ie) THEN
                    HT(ib,ia) = HT(ib,ia) +
     &              0.25*T2(M,N,if,ib)*ENE_ORB(ia)*T2(M,N,ie,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(ie.EQ.if) THEN
                    HT(ib,ia) = HT(ib,ia) +
     &              0.25*T2(M,N,ie,ib)*ENE_ORB(ie)*T2(M,N,ia,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO K = 1, NOS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(K.EQ.N) THEN
                    HT(ib,ia) = HT(ib,ia) -
     &              0.5*T2(M,K,ie,ib)*ENE_ORB(K)*T2(N,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
          DO M = 1, NOS
            DO N = 1, NOS
              IF(M.EQ.N) THEN
                HT(ia,ib) = HT(ia,ib) + 0.5*T1(M,ib)*ENE_ORB(M)*T1(N,ia)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              IF(ib.EQ.ie) THEN
                HT(ia,ib)=HT(ia,ib) -0.5d0*T1(M,ie)*ENE_ORB(ib)*T1(M,ia)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(ia.EQ.ie) THEN
                    HT(ia,ib) = HT(ia,ib) +
     &              0.25*T2(M,N,if,ib)*ENE_ORB(ia)*T2(M,N,ie,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(ie.EQ.if) THEN
                    HT(ia,ib) = HT(ia,ib) +
     &              0.25*T2(M,N,ie,ib)*ENE_ORB(ie)*T2(M,N,ia,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO K = 1, NOS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(K.EQ.N) THEN
                    HT(ia,ib) = HT(ia,ib) -
     &              0.5*T2(M,K,ie,ib)*ENE_ORB(K)*T2(N,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   DO IA = NOS+1, NOS+NVS
      !     DO IB = NOS+1, NOS+NVS
      !       WRITE(6,*)HT(ia,ib), HT(ib,ia), HT(ia,ib)-HT(ib,ia)
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'PROOF THAT H (COLUMN 1) CHANGED TO HT (COLUMN 2)'

      !   DO IA = NOS+1, NOS+NVS
      !     DO IB = NOS+1, NOS+NVS
      !       WRITE(6,*)H(ia,ib), HT(ia,ib), H(ia,ib)-HT(ia,ib)
      !     ENDDO
      !   ENDDO
      ! ENDIF
! ***| END DEBUG
!  END  11111111111111111111111111111111111111111111111111111111111111111111111

! START 22222222222222222222222222222222222222222222222222222222222222222222222
! ***| HT_{I}^{J}
      DO I = 1, NOS
        DO J = 1, NOS
! ***| H_N
          IF(I.EQ.J) HT(I,J) = HT(I,J) + ENE_ORB(I)

! ***| (H_N(T_EXT))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(I,J) = HT(I,J) + V(ie,I,M,J)*T1(M,ie)
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                HT(I,J) = HT(I,J) + 0.5*V(ie,if,M,J)*T2(M,I,ie,if)
              ENDDO
            ENDDO
          ENDDO
! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(J,I) = HT(J,I) + V(M,J,ie,I)*T1(M,ie)
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                HT(J,I) = HT(J,I) + 0.5*V(M,J,ie,if)*T2(M,I,ie,if)
              ENDDO
            ENDDO
          ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              IF(ie.EQ.if) THEN
                HT(I,J) = HT(I,J) +
     &          0.5*T1(J,ie)*ENE_ORB(ie)*T1(I,if)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              IF(J.EQ.M) THEN
                HT(I,J) = HT(I,J) -
     &          0.5*T1(M,ie)*ENE_ORB(J)*T1(I,ie)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N =  1, NOS
                  IF(I.EQ.N) THEN
                    HT(I,J) = HT(I,J) -
     &              0.25*T2(M,J,ie,if)*ENE_ORB(I)*T2(M,N,ie,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO G = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ie.EQ.if) THEN
                    HT(I,J) = HT(I,J) +
     &              0.5*T2(M,J,ie,G)*ENE_ORB(ie)*T2(M,I,if,G)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(M.EQ.N) THEN
                    HT(I,J) = HT(I,J) +
     &              0.25*T2(I,M,ie,if)*ENE_ORB(M)*T2(N,J,ie,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              IF(ie.EQ.if) THEN
                HT(J,I) = HT(J,I) +
     &          0.5*T1(J,ie)*ENE_ORB(ie)*T1(I,if)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              IF(J.EQ.M) THEN
                HT(J,I) = HT(J,I) -
     &          0.5*T1(M,ie)*ENE_ORB(J)*T1(I,ie)
              ENDIF
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N =  1, NOS
                  IF(I.EQ.N) THEN
                    HT(J,I) = HT(J,I) -
     &              0.25*T2(M,J,ie,if)*ENE_ORB(I)*T2(M,N,ie,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO G = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ie.EQ.if) THEN
                    HT(J,I) = HT(J,I) +
     &              0.5*T2(M,J,ie,G)*ENE_ORB(ie)*T2(M,I,if,G)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(M.EQ.N) THEN
                    HT(J,I) = HT(J,I) +
     &              0.25*T2(I,M,ie,if)*ENE_ORB(M)*T2(N,J,ie,if)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO

        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       WRITE(6,*)HT(I,J), HT(J,I), HT(I,J)-HT(J,I)
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'PROOF THAT H (COLUMN 1) CHANGED TO HT (COLUMN 2)'

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       WRITE(6,*)H(I,J), HT(I,J), H(I,J)-HT(I,J)
      !     ENDDO
      !   ENDDO
      ! ENDIF
! ***| END DEBUG
!  END  22222222222222222222222222222222222222222222222222222222222222222222222

! START 33333333333333333333333333333333333333333333333333333333333333333333333
! ***| HT_{A}^{I}
      DO IA = NOS+1, NOS+NVS
        DO I = 1, NOS
! ***| THERE IS NO H_N TERM SINCE WE ARE ASSUMING AN RHF REFERENCE.
! ***| (H_N(T_EXT))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(ia,I) = HT(ia,I) + V(ie,ia,M,I)*T1(M,ie)
            ENDDO
          ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
          DO IE = NOS+1, NOS+NVS
            IF(ia.EQ.ie) THEN
              HT(ia,I) = HT(ia,I) + ENE_ORB(ia)*T1(I,ie)
            ENDIF
          ENDDO

          DO M = 1, NOS
            IF(I.EQ.M) THEN
              HT(ia,I) = HT(ia,I) - ENE_ORB(I)*T1(M,ia)
            ENDIF
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(ia,I) = HT(ia,I) + V(M,ia,ie,I)*T1(M,ie)
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                HT(ia,I) = HT(ia,I) - 0.5*V(M,N,ie,I)*T2(M,N,ie,ia)
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                HT(ia,I) = HT(ia,I) + 0.5*V(M,ia,ie,if)*T2(M,I,ie,if)
              ENDDO
            ENDDO
          ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                IF(M.EQ.N) THEN
                  HT(ia,I) = HT(ia,I) -
     &            0.5*T2(m,i,ie,ia)*ENE_ORB(M)*T1(n,ie)
                ENDIF
              ENDDO
            ENDDO
          ENDDO


          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                IF(ie.EQ.if) THEN
                  HT(ia,I) = HT(ia,I) +
     &            0.5*T2(m,i,if,ia)*ENE_ORB(ie)*T1(M,ie)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                IF(ia.EQ.if) THEN
                  HT(ia,I) = HT(ia,I) +
     &            0.5*T1(M,ie)*ENE_ORB(ia)*T2(M,I,ie,if)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                IF(I.EQ.N) THEN
                  HT(ia,I) = HT(ia,I) -
     &            0.5*T1(M,ie)*ENE_ORB(I)*T2(M,N,ie,ia)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                IF(M.EQ.N) THEN
                  HT(ia,I) = HT(ia,I) -
     &            0.5*T1(M,ie)*ENE_ORB(M)*T2(I,N,ia,ie)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                IF(ie.EQ.if) THEN
                  HT(ia,I) = HT(ia,I) +
     &            0.5*T1(M,ie)*ENE_ORB(ie)*T2(I,M,ia,if)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   WRITE(6,*) 'PROOF THAT H (COLUMN 1) CHANGED TO HT (COLUMN 2)'

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       WRITE(6,*)H(ia,I), HT(ia,I), H(ia,I)-HT(ia,I)
      !     ENDDO
      !   ENDDO
      ! ENDIF
! ***| END DEBUG
!  END  33333333333333333333333333333333333333333333333333333333333333333333333

! START 44444444444444444444444444444444444444444444444444444444444444444444444
! ***| HT_{I}^{A}
      DO I = 1, NOS
        DO IA = NOS+1, NOS+NVS
! ***| THERE IS NO H_N TERM SINCE WE ARE ASSUMING AN RHF REFERENCE.
! ***| (H_N(T_EXT))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            IF(ia.EQ.ie) THEN
              HT(I,ia) = HT(I,ia) + ENE_ORB(ia)*T1(I,ie)
            ENDIF
          ENDDO

          DO M = 1, NOS
            IF(I.EQ.M) THEN
              HT(I,ia) = HT(I,ia) - ENE_ORB(I)*T1(M,ia)
            ENDIF
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(I,ia) = HT(I,ia) + V(ie,I,M,ia)*T1(M,ie)
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                HT(I,ia) = HT(I,ia) - 0.5*V(ie,I,M,N)*T2(M,N,ie,ia)
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                HT(I,ia) = HT(I,ia) + 0.5*V(ie,if,M,ia)*T2(M,I,ie,if)
              ENDDO
            ENDDO
          ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              HT(I,ia) = HT(I,ia) + V(M,I,ie,ia)*T1(M,ie)
            ENDDO
          ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                IF(ia.EQ.if) THEN
                  HT(I,ia) = HT(I,ia) +
     &            0.5*T1(M,ie)*ENE_ORB(ia)*T2(M,I,ie,if)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                IF(I.EQ.N) THEN
                  HT(I,ia) = HT(I,ia) -
     &            0.5*T1(M,ie)*ENE_ORB(I)*T2(M,N,ie,ia)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                IF(M.EQ.N) THEN
                  HT(I,ia) = HT(I,ia) -
     &            0.5*T1(M,ie)*ENE_ORB(M)*T2(I,N,ia,ie)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                IF(ie.EQ.if) THEN
                  HT(I,ia) = HT(I,ia) +
     &            0.5*T1(M,ie)*ENE_ORB(ie)*T2(I,M,ia,if)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
          DO IE = NOS+1, NOS+NVS
            DO M = 1, NOS
              DO N = 1, NOS
                IF(M.EQ.N) THEN
                  HT(I,ia) = HT(I,ia) -
     &            0.5*T2(m,i,ie,ia)*ENE_ORB(M)*T1(n,ie)
                ENDIF
              ENDDO
            ENDDO
          ENDDO


          DO IE = NOS+1, NOS+NVS
            DO IF = NOS+1, NOS+NVS
              DO M = 1, NOS
                IF(ie.EQ.if) THEN
                  HT(I,ia) = HT(I,ia) +
     &            0.5*T2(m,i,if,ia)*ENE_ORB(ie)*T1(M,ie)
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   DO IA = NOS+1, NOS+NVS
      !     DO I = 1, NOS
      !       WRITE(6,*)HT(ia,I), HT(I,ia), HT(ia,I)-HT(I,ia)
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'PROOF THAT H (COLUMN 1) CHANGED TO HT (COLUMN 2)'

      !   DO IA = NOS+1, NOS+NVS
      !     DO I = 1, NOS
      !       WRITE(6,*)H(I,ia), HT(I,ia), H(I,ia)-HT(I,ia)
      !     ENDDO
      !   ENDDO
      ! ENDIF
! ***| END DEBUG
!  END  44444444444444444444444444444444444444444444444444444444444444444444444

! START 55555555555555555555555555555555555555555555555555555555555555555555555
! ***| VT_{IA}^{BC} = VT_{AI}^{CB}
      DO I = 1, NOS
        DO IA = NOS+1, NOS+NVS
          DO IB = NOS+1, NOS+NVS
            DO IC = NOS+1, NOS+NVS
! ***| H_N
              VT(I,ia,ib,ic) = VT(I,ia,ib,ic) + V(I,ia,ib,ic)

! ***| (H_N(T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(I,ia,ib,ic) = VT(I,ia,ib,ic) +
     &          V(ie,ia,ib,ic)*T1(I,ie)
              ENDDO

              DO M = 1, NOS
                VT(I,ia,ib,ic) = VT(I,ia,ib,ic) -
     &          V(ia,I,M,ib)*T1(M,ic)
              ENDDO

              DO M = 1, NOS
                VT(I,ia,ib,ic) = VT(I,ia,ib,ic) +
     &          V(ia,I,M,ic)*T1(M,ib)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,ia,ib,ic) = VT(I,ia,ib,ic) -
     &            V(ie,ia,M,ib)*T2(M,I,ie,ic)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,ia,ib,ic) = VT(I,ia,ib,ic) +
     &            V(ie,ia,M,ic)*T2(M,I,ie,ib)
                ENDDO
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  VT(I,ia,ib,ic) = VT(I,ia,ib,ic) +
     &            0.5*V(I,ia,M,N)*T2(M,N,ib,ic)
                ENDDO
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
! ***| THIS IS WRITTEN AS GAMMA_{CI}^{AB} = GAMMA_{IC}^{BA} IN THE PAPER
              DO M = 1, NOS
                VT(I,ic,ib,ia) = VT(I,ic,ib,ia) -
     &          V(M,I,ia,ib)*T1(M,ic)
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
! ***| THIS IS WRITTEN AS GAMMA_{CI}^{AB} = GAMMA_{IC}^{BA} IN THE PAPER
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(M.EQ.N) THEN
                    VT(I,ic,ib,ia) = VT(I,ic,ib,ia) +
     &              0.5*T1(M,ic)*ENE_ORB(M)*T2(N,I,ia,ib)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ib.EQ.ie) THEN
                    VT(I,ic,ib,ia) = VT(I,ic,ib,ia) -
     &              0.5*T1(M,ic)*ENE_ORB(ib)*T2(I,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ia.EQ.ie) THEN
                    VT(I,ic,ib,ia) = VT(I,ic,ib,ia) +
     &              0.5*T1(M,ic)*ENE_ORB(ia)*T2(I,M,ie,ib)
                  ENDIF
                ENDDO
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  IF(I.EQ.N) THEN
                    VT(I,ic,ib,ia) = VT(I,ic,ib,ia) +
     &              0.5*T1(M,ic)*ENE_ORB(I)*T2(M,N,ia,ib)
                  ENDIF
                ENDDO
              ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
! ***| THIS IS WRITTEN AS GAMMA_{CI}^{AB} = GAMMA_{IC}^{BA} IN THE PAPER
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ic.EQ.ie) THEN
                    VT(I,ic,ib,ia) = VT(I,ic,ib,ia) -
     &              0.5*T2(M,I,ia,ib)*ENE_ORB(ic)*T1(M,ie)
                  ENDIF
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| FILLING IN THE ARRAY
! ***| VT_{IA}^{BC} = VT_{AI}^{CB}
      DO I = 1, NOS
        DO IA = NOS+1, NOS+NVS
          DO IB = NOS+1, NOS+NVS
            DO IC = NOS+1, NOS+NVS
              VT(ia,I,ic,ib) = VT(I,ia,ib,ic)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

      !   DO I = 1, NOS
      !     DO IA = NOS+1, NOS+NVS
      !       DO IB = NOS+1, NOS+NVS
      !         DO IC = NOS+1, NOS+NVS
      !           WRITE(6,*)V(I,ia,ib,ic), VT(I,ia,ib,ic), V(I,ia,ib,ic)-VT(I,ia,ib,ic)
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG
!  END  55555555555555555555555555555555555555555555555555555555555555555555555

! START 66666666666666666666666666666666666666666666666666666666666666666666666
! ***| VT_{IJ}^{KA} = VT_{JI}^{AK}
      DO I = 1, NOS
        DO J = 1, NOS
          DO K = 1, NOS
            DO IA = NOS+1, NOS+NVS
! ***| H_N
              VT(I,J,K,ia) = VT(I,J,K,ia) + V(I,J,K,ia)

! ***| (H_N(T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(I,J,K,ia) = VT(I,J,K,ia) +
     &          V(ie,J,K,ia)*T1(I,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                VT(I,J,K,ia) = VT(I,J,K,ia) -
     &          V(ie,I,K,ia)*T1(J,ie)
              ENDDO

              DO M = 1, NOS
                VT(I,J,K,ia) = VT(I,J,K,ia) +
     &          V(I,J,M,K)*T1(M,ia)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  VT(I,J,K,ia) = VT(I,J,K,ia) +
     &            0.5*V(ie,if,K,ia)*T2(I,J,ie,if)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,J,K,ia) = VT(I,J,K,ia) -
     &            V(ie,J,M,K)*T2(M,I,ie,ia)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,J,K,ia) = VT(I,J,K,ia) +
     &            V(ie,I,M,K)*T2(M,J,ie,ia)
                ENDDO
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(I,J,K,ia) = VT(I,J,K,ia) +
     &          V(I,J,ie,ia)*T1(K,ie)
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  IF(ie.EQ.if) THEN
                    VT(I,J,K,ia) = VT(I,J,K,ia) +
     &              0.5*T1(K,ie)*ENE_ORB(ie)*T2(I,J,if,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(J.EQ.M) THEN
                    VT(I,J,K,ia) = VT(I,J,K,ia) -
     &              0.5*T1(K,ie)*ENE_ORB(J)*T2(I,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(I.EQ.M) THEN
                    VT(I,J,K,ia) = VT(I,J,K,ia) +
     &              0.5*T1(K,ie)*ENE_ORB(I)*T2(J,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  IF(ia.EQ.if) THEN
                    VT(I,J,K,ia) = VT(I,J,K,ia) +
     &              0.5*T1(K,ie)*ENE_ORB(ia)*T2(I,J,ie,if)
                  ENDIF
                ENDDO
              ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(K.EQ.M) THEN
                    VT(I,J,K,ia) = VT(I,J,K,ia) -
     &              0.5*T2(I,J,ie,ia)*ENE_ORB(K)*T1(M,ie)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  IF(ie.EQ.if) THEN
                    VT(I,J,K,ia) = VT(I,J,K,ia) +
     &              0.5*T2(j,i,ia,if)*ENE_ORB(ie)*T1(k,ie)
                  ENDIF
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| FILLING IN THE ARRAY
! ***| VT_{IJ}^{KA} = VT_{JI}^{AK}
      DO I = 1, NOS
        DO J = 1, NOS
          DO K = 1, NOS
            DO IA = NOS+1, NOS+NVS
              VT(J,I,ia,K) = VT(I,J,K,ia)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       DO K = 1, NOS
      !         DO IA = NOS+1, NOS+NVS
      !           WRITE(6,*)V(I,J,K,ia), VT(I,J,K,ia), V(I,J,K,ia)-VT(I,J,K,ia)
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG
!  END  66666666666666666666666666666666666666666666666666666666666666666666666

! START 77777777777777777777777777777777777777777777777777777777777777777777777
! ***| VT_{AB}^{CI} =  VT_{BA}^{IC}
      DO IA = NOS+1, NOS+NVS
        DO IB = NOS+1, NOS+NVS
          DO IC = NOS+1, NOS+NVS
            DO I = 1, NOS
! ***| H_N
              VT(ia,ib,ic,I) = VT(ia,ib,ic,I) + V(ia,ib,ic,I)


! ***| (H_N(T_EXT))_C,OPEN
              DO M = 1, NOS
                VT(ia,ib,ic,I) = VT(ia,ib,ic,I) -
     &          V(ia,ib,M,I)*T1(M,ic)
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
! ***| THIS IS WRITTEN AS GAMMA_{BC}^{IA} = GAMMA_{CB}^{AI} IN THE PAPER
              DO IE = NOS+1, NOS+NVS
                VT(ic,ib,ia,I) = VT(ic,ib,ia,I) +
     &          V(ib,ic,ie,ia)*T1(I,ie)
              ENDDO

              DO M = 1, NOS
                VT(ic,ib,ia,I) = VT(ic,ib,ia,I) -
     &          V(M,ib,ia,I)*T1(M,ic)
              ENDDO

              DO M = 1, NOS
                VT(ic,ib,ia,I) = VT(ic,ib,ia,I) +
     &          V(M,ic,ia,I)*T1(M,ib)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(ic,ib,ia,I) = VT(ic,ib,ia,I) -
     &            V(M,ib,ie,ia)*T2(M,I,ie,ic)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(ic,ib,ia,I) = VT(ic,ib,ia,I) +
     &            V(M,ic,ie,ia)*T2(M,I,ie,ib)
                ENDDO
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  VT(ic,ib,ia,I) = VT(ic,ib,ia,I) +
     &            0.5*V(M,N,I,ia)*T2(M,N,ib,ic)
                ENDDO
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ic.EQ.ie) THEN
                    VT(ia,ib,ic,I) = VT(ia,ib,ic,I) -
     &              0.5*T2(M,I,ia,ib)*ENE_ORB(ic)*T1(M,ie)
                  ENDIF
                ENDDO
              ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
              DO M = 1, NOS
                DO N = 1, NOS
                  IF(M.EQ.N) THEN
                    VT(ia,ib,ic,I) = VT(ia,ib,ic,I) +
     &              0.5*T1(M,ic)*ENE_ORB(M)*T2(N,I,ia,ib)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ib.EQ.ie) THEN
                    VT(ia,ib,ic,I) = VT(ia,ib,ic,I) -
     &              0.5*T1(M,ic)*ENE_ORB(ib)*T2(I,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(ia.EQ.ie) THEN
                    VT(ia,ib,ic,I) = VT(ia,ib,ic,I) +
     &              0.5*T1(M,ic)*ENE_ORB(ia)*T2(I,M,ie,ib)
                  ENDIF
                ENDDO
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  IF(I.EQ.N) THEN
                    VT(ia,ib,ic,I) = VT(ia,ib,ic,I) +
     &              0.5*T1(M,ic)*ENE_ORB(I)*T2(M,N,ia,ib)
                  ENDIF
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| FILLING IN THE ARRAY
! ***| VT_{AB}^{CI} =  VT_{BA}^{IC}
      DO IA = NOS+1, NOS+NVS
        DO IB = NOS+1, NOS+NVS
          DO IC = NOS+1, NOS+NVS
            DO I = 1, NOS
              VT(ib,ia,I,ic) = VT(ia,ib,ic,I)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
    !   IF(NODEZERO) THEN

    !     DO IA = NOS+1, NOS+NVS
    !       DO IB = NOS+1, NOS+NVS
    !         DO IC = NOS+1, NOS+NVS
    !           DO I = 1, NOS
    !             WRITE(6,*)VT(I,ia,ib,ic),VT(ic,ib,ia,I),
    !  &          VT(I,ia,ib,ic)-VT(ic,ib,ia,I)
    !           ENDDO
    !         ENDDO
    !       ENDDO
    !     ENDDO

    !     WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

    !     DO IA = NOS+1, NOS+NVS
    !       DO IB = NOS+1, NOS+NVS
    !         DO IC = NOS+1, NOS+NVS
    !           DO I = 1, NOS
    !             WRITE(6,*)V(ia,ib,ic,I), VT(ia,ib,ic,I), V(ia,ib,ic,I)-VT(ia,ib,ic,I)
    !           ENDDO
    !         ENDDO
    !       ENDDO
    !     ENDDO

    !   ENDIF
! ***| END DEBUG
!  END  77777777777777777777777777777777777777777777777777777777777777777777777

! START 88888888888888888888888888888888888888888888888888888888888888888888888
! ***| VT_{KA}^{IJ} = VT_{AK}^{JI}
      DO K = 1, NOS
        DO IA = NOS+1, NOS+NVS
          DO I = 1, NOS
            DO J = 1, NOS
! ***| H_N
              VT(K,ia,I,J) = VT(K,ia,I,J) + V(K,ia,I,J)

! ***| (H_N(T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(K,ia,I,J) = VT(K,ia,I,J) +
     &          V(ie,ia,I,J)*T1(K,ie)
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(K,ia,I,J) = VT(K,ia,I,J) +
     &          V(K,ia,ie,J)*T1(I,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                VT(K,ia,I,J) = VT(K,ia,I,J) -
     &          V(K,ia,ie,I)*T1(J,ie)
              ENDDO

              DO M = 1, NOS
                VT(K,ia,I,J) = VT(K,ia,I,J) +
     &          V(M,K,I,J)*T1(M,ia)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  VT(K,ia,I,J) = VT(K,ia,I,J) +
     &            0.5*V(K,ia,ie,if)*T2(I,J,ie,if)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(K,ia,I,J) = VT(K,ia,I,J) -
     &            V(M,K,ie,J)*T2(M,I,ie,ia)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(K,ia,I,J) = VT(K,ia,I,J) +
     &            V(M,K,ie,I)*T2(M,J,ie,ia)
                ENDDO
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(K.EQ.M) THEN
                    VT(K,ia,I,J) = VT(K,ia,I,J) -
     &              0.5*T2(I,J,ie,ia)*ENE_ORB(K)*T1(M,ie)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  IF(ie.EQ.if) THEN
                    VT(K,ia,I,J) = VT(K,ia,I,J) +
     &              0.5*T2(j,i,ia,if)*ENE_ORB(ie)*T1(k,ie)
                  ENDIF
                ENDDO
              ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  IF(ie.EQ.if) THEN
                    VT(K,ia,I,J) = VT(K,ia,I,J) +
     &              0.5*T1(K,ie)*ENE_ORB(ie)*T2(I,J,if,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(J.EQ.M) THEN
                    VT(K,ia,I,J) = VT(K,ia,I,J) -
     &              0.5*T1(K,ie)*ENE_ORB(J)*T2(I,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  IF(I.EQ.M) THEN
                    VT(K,ia,I,J) = VT(K,ia,I,J) +
     &              0.5*T1(K,ie)*ENE_ORB(I)*T2(J,M,ie,ia)
                  ENDIF
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  IF(ia.EQ.if) THEN
                    VT(K,ia,I,J) = VT(K,ia,I,J) +
     &              0.5*T1(K,ie)*ENE_ORB(ia)*T2(I,J,ie,if)
                  ENDIF
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| FILLING IN THE ARRAY
! ***| VT_{KA}^{IJ} = VT_{AK}^{JI}
      DO K = 1, NOS
        DO IA = NOS+1, NOS+NVS
          DO I = 1, NOS
            DO J = 1, NOS
              VT(ia,K,J,I) = VT(K,ia,I,J)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
    !   IF(NODEZERO) THEN

    !     DO K = 1, NOS
    !       DO IA = NOS+1, NOS+NVS
    !         DO I = 1, NOS
    !           DO J = 1, NOS
    !             WRITE(6,*)VT(K,ia,I,J),VT(I,J,K,ia),
    !  &          VT(K,ia,I,J)-VT(I,J,K,ia)
    !           ENDDO
    !         ENDDO
    !       ENDDO
    !     ENDDO

    !     WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

    !     DO K = 1, NOS
    !       DO IA = NOS+1, NOS+NVS
    !         DO I = 1, NOS
    !           DO J = 1, NOS
    !             WRITE(6,*)V(K,ia,I,J), VT(K,ia,I,J), V(K,ia,I,J)-VT(K,ia,I,J)
    !           ENDDO
    !         ENDDO
    !       ENDDO
    !     ENDDO

    !   ENDIF
! ***| END DEBUG
!  END  88888888888888888888888888888888888888888888888888888888888888888888888

! START 99999999999999999999999999999999999999999999999999999999999999999999999
! ***| VT_{IA}^{JB} = VT_{AI}^{BJ} = -VT_{IA}^{BJ} = -VT_{AI}^{JB}
      DO I = 1, NOS
        DO IA = NOS+1, NOS+NVS
          DO J = 1, NOS
            DO IB = NOS+1, NOS+NVS
! ***| H_N
              VT(I,ia,J,ib) = VT(I,ia,J,ib) + V(I,ia,J,ib)

! ***| (H_N(T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(I,ia,J,ib) = VT(I,ia,J,ib) +
     &          V(ie,ia,J,ib)*T1(I,ie)
              ENDDO

              DO M = 1, NOS
                VT(I,ia,J,ib) = VT(I,ia,J,ib) +
     &          V(I,ia,M,J)*T1(M,ib)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,ia,J,ib) = VT(I,ia,J,ib) -
     &            V(ie,ia,M,J)*T2(M,I,ie,ib)
                ENDDO
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(J,ib,I,ia) = VT(J,ib,I,ia) +
     &          V(J,ib,ie,ia)*T1(I,ie)
              ENDDO

              DO M = 1, NOS
                VT(J,ib,I,ia) = VT(J,ib,I,ia) +
     &          V(M,J,I,ia)*T1(M,ib)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(J,ib,I,ia) = VT(J,ib,I,ia) -
     &            V(M,J,ie,ia)*T2(M,I,ie,ib)
                ENDDO
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(J.EQ.N) THEN
                      VT(J,ib,I,ia) = VT(J,ib,I,ia) +
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(J)*T2(M,N,ie,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(ia.EQ.if) THEN
                      VT(J,ib,I,ia) = VT(J,ib,I,ia) -
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(ia)*T2(M,J,ie,if)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(ie.EQ.if) THEN
                      VT(J,ib,I,ia) = VT(J,ib,I,ia) -
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(ie)*T2(M,J,if,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M =  1, NOS
                  DO N = 1, NOS
                    IF(M.EQ.N) THEN
                      VT(J,ib,I,ia) = VT(J,ib,I,ia) +
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(M)*T2(N,J,ie,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(J.EQ.N) THEN
                      VT(I,ia,J,ib) = VT(I,ia,J,ib) +
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(J)*T2(M,N,ie,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(ia.EQ.if) THEN
                      VT(I,ia,J,ib) = VT(I,ia,J,ib) -
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(ia)*T2(M,J,ie,if)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(ie.EQ.if) THEN
                      VT(I,ia,J,ib) = VT(I,ia,J,ib) -
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(ie)*T2(M,J,if,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M =  1, NOS
                  DO N = 1, NOS
                    IF(M.EQ.N) THEN
                      VT(I,ia,J,ib) = VT(I,ia,J,ib) +
     &                0.5*T2(M,I,ie,ib)*ENE_ORB(M)*T2(N,J,ie,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| FILLING IN THE ARRAY
! ***| VT_{IA}^{JB} = VT_{AI}^{BJ} = -VT_{IA}^{BJ} = -VT_{AI}^{JB}
      DO I = 1, NOS
        DO IA = NOS+1, NOS+NVS
          DO J = 1, NOS
            DO IB = NOS+1, NOS+NVS
              VT(ia,I,J,ib) = -VT(I,ia,J,ib)
              VT(I,ia,ib,J) = -VT(I,ia,J,ib)
              VT(ia,I,ib,J) =  VT(I,ia,J,ib)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
    !   IF(NODEZERO) THEN

    !     DO I = 1, NOS
    !       DO IA = NOS+1, NOS+NVS
    !         DO J = 1, NOS
    !           DO IB = NOS+1, NOS+NVS
    !             WRITE(6,*)VT(I,ia,J,ib),VT(J,ib,I,ia),
    !  &          VT(I,ia,J,ib)-VT(J,ib,I,ia)
    !           ENDDO
    !         ENDDO
    !       ENDDO
    !     ENDDO

    !     WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

    !     DO I = 1, NOS
    !       DO IA = NOS+1, NOS+NVS
    !         DO J = 1, NOS
    !           DO IB = NOS+1, NOS+NVS
    !             WRITE(6,*)V(I,ia,J,ib), VT(I,ia,J,ib), V(I,ia,J,ib)-VT(I,ia,J,ib)
    !           ENDDO
    !         ENDDO
    !       ENDDO
    !     ENDDO

    !   ENDIF
! ***| END DEBUG
!  END  99999999999999999999999999999999999999999999999999999999999999999999999

! START 10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10
! ***| VT_{AB}^{CD}
      DO IA = NOS+1, NOS+NVS
        DO IB = NOS+1, NOS+NVS
          DO IC = NOS+1, NOS+NVS
            DO ID= NOS+1, NOS+NVS
! ***| H_N
              VT(ia,ib,ic,id) = VT(ia,ib,ic,id) + V(ia,ib,ic,id)

! ***| (H_N(T_EXT))_C,OPEN
              DO M = 1, NOS
                VT(ia,ib,ic,id) = VT(ia,ib,ic,id) +
     &          V(ia,ib,M,ic)*T1(M,id)
              ENDDO

              DO M = 1, NOS
                VT(ia,ib,ic,id) = VT(ia,ib,ic,id) -
     &          V(ia,ib,M,id)*T1(M,ic)
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  VT(ia,ib,ic,id) = VT(ia,ib,ic,id) +
     &            0.5*V(ia,ib,M,N)*T2(M,N,ic,id)
                ENDDO
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
              DO M = 1, NOS
                VT(ic,id,ia,ib) = VT(ic,id,ia,ib) +
     &          V(M,ic,ia,ib)*T1(M,id)
              ENDDO

              DO M = 1, NOS
                VT(ic,id,ia,ib) = VT(ic,id,ia,ib) -
     &          V(M,id,ia,ib)*T1(M,ic)
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  VT(ic,id,ia,ib) = VT(ic,id,ia,ib) +
     &            0.5*V(M,N,ia,ib)*T2(M,N,ic,id)
                ENDDO
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(ia.EQ.ie) THEN
                      VT(ic,id,ia,ib) = VT(ic,id,ia,ib) +
     &                0.25*T2(M,N,ic,id)*ENE_ORB(ia)*T2(M,N,ie,ib)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(ib.EQ.ie) THEN
                      VT(ic,id,ia,ib) = VT(ic,id,ia,ib) -
     &                0.25*T2(M,N,ic,id)*ENE_ORB(ib)*T2(M,N,ie,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO K = 1, NOS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(M.EQ.N) THEN
                      VT(ic,id,ia,ib) = VT(ic,id,ia,ib) -
     &                0.5*T2(M,K,ic,id)*ENE_ORB(M)*T2(N,K,ia,ib)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(ia.EQ.ie) THEN
                      VT(ia,ib,ic,id) = VT(ia,ib,ic,id) +
     &                0.25*T2(M,N,ic,id)*ENE_ORB(ia)*T2(M,N,ie,ib)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(ib.EQ.ie) THEN
                      VT(ia,ib,ic,id) = VT(ia,ib,ic,id) -
     &                0.25*T2(M,N,ic,id)*ENE_ORB(ib)*T2(M,N,ie,ia)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO K = 1, NOS
                DO M = 1, NOS
                  DO N = 1, NOS
                    IF(M.EQ.N) THEN
                      VT(ia,ib,ic,id) = VT(ia,ib,ic,id) -
     &                0.5*T2(M,K,ic,id)*ENE_ORB(M)*T2(N,K,ia,ib)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   DO IA = NOS+1, NOS+NVS
      !     DO IB = NOS+1, NOS+NVS
      !       DO IC = NOS+1, NOS+NVS
      !         DO ID= NOS+1, NOS+NVS
      !           ZZZ=VT(ia,ib,ic,id)-VT(ic,id,ia,ib)
      !           IF(abs(ZZZ).GT.(1.0d-17)) THEN
      !             WRITE(6,*)VT(ia,ib,ic,id),VT(ic,id,ia,ib),ZZZ
      !           ENDIF
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

      !   DO IA = NOS+1, NOS+NVS
      !     DO IB = NOS+1, NOS+NVS
      !       DO IC = NOS+1, NOS+NVS
      !         DO ID= NOS+1, NOS+NVS
      !           WRITE(6,*)V(ia,ib,ic,id), VT(ia,ib,ic,id), V(ia,ib,ic,id)-VT(ia,ib,ic,id)
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG
!  END  10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10-10

! START 11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11
! ***| VT_{IJ}^{KL}
      DO I = 1, NOS
        DO J = 1, NOS
          DO K = 1, NOS
            DO L = 1, NOS
! ***| H_N
              VT(I,J,K,L) = VT(I,J,K,L) + V(I,J,K,L)

! ***| (H_N(T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(I,J,K,L) = VT(I,J,K,L) +
     &          V(ie,J,K,L)*T1(I,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                VT(I,J,K,L) = VT(I,J,K,L) -
     &          V(ie,I,K,L)*T1(J,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  VT(I,J,K,L) = VT(I,J,K,L) +
     &            0.5*V(ie,if,K,L)*T2(I,J,ie,if)
                ENDDO
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(K,L,I,J) = VT(K,L,I,J) +
     &          V(K,L,ie,J)*T1(I,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                VT(K,L,I,J) = VT(K,L,I,J) -
     &          V(K,L,ie,I)*T1(J,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  VT(K,L,I,J) = VT(K,L,I,J) +
     &            0.5*V(K,L,ie,if)*T2(I,J,ie,if)
                ENDDO
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(L.EQ.M) THEN
                      VT(K,L,I,J) = VT(K,L,I,J) +
     &                0.25*T2(I,J,ie,if)*ENE_ORB(L)*T2(M,K,ie,if)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(K.EQ.M) THEN
                      VT(K,L,I,J) = VT(K,L,I,J) -
     &                0.25*T2(I,J,ie,if)*ENE_ORB(K)*T2(M,L,ie,if)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO G = NOS+1, NOS+NVS
                    IF(ie.EQ.if) THEN
                      VT(K,L,I,J) = VT(K,L,I,J) +
     &                0.5*T2(I,J,ie,G)*ENE_ORB(ie)*T2(K,L,if,G)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(L.EQ.M) THEN
                      VT(I,J,K,L) = VT(I,J,K,L) +
     &                0.25*T2(I,J,ie,if)*ENE_ORB(L)*T2(M,K,ie,if)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO M = 1, NOS
                    IF(K.EQ.M) THEN
                      VT(I,J,K,L) = VT(I,J,K,L) -
     &                0.25*T2(I,J,ie,if)*ENE_ORB(K)*T2(M,L,ie,if)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  DO G = NOS+1, NOS+NVS
                    IF(ie.EQ.if) THEN
                      VT(I,J,K,L) = VT(I,J,K,L) +
     &                0.5*T2(I,J,ie,G)*ENE_ORB(ie)*T2(K,L,if,G)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       DO K = 1, NOS
      !         DO L = 1, NOS
      !           ZZZ=VT(I,J,K,L)-VT(K,L,I,J)
      !           IF(abs(ZZZ).GT.(1.0d-17)) THEN
      !             WRITE(6,*)VT(I,J,K,L),VT(K,L,I,J),ZZZ
      !           ENDIF
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       DO K = 1, NOS
      !         DO L = 1, NOS
      !           WRITE(6,*)V(I,J,K,L), VT(I,J,K,L), V(I,J,K,L)-VT(I,J,K,L)
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG
!  END  11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11-11

! START 12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12
! ***| VT_{AB}^{IJ}
      DO IA = NOS+1, NOS+NVS
        DO IB = NOS+1, NOS+NVS
          DO I = 1, NOS
            DO J = 1, NOS
! ***| H_N
              VT(ia,ib,I,J) = VT(ia,ib,I,J) + V(ia,ib,I,J)

! ***| (H_N(T_EXT))_C,OPEN
! ***| THERE IS NO CORRESPONDING TERM

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &          V(ia,ib,ie,J)*T1(I,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                VT(I,J,ia,ib) = VT(I,J,ia,ib) -
     &          V(ia,ib,ie,I)*T1(J,ie)
              ENDDO

              DO M = 1, NOS
                VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &          V(M,ia,I,J)*T1(M,ib)
              ENDDO

              DO M = 1, NOS
                VT(I,J,ia,ib) = VT(I,J,ia,ib) -
     &          V(M,ib,I,J)*T1(M,ia)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                IF(ia.EQ.ie) THEN
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &            ENE_ORB(ia)*T2(I,J,ie,ib)
                ENDIF
              ENDDO

              DO IE = NOS+1, NOS+NVS
                IF(ib.EQ.ie) THEN
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) -
     &            ENE_ORB(ib)*T2(I,J,ie,ia)
                ENDIF
              ENDDO

              DO M = 1, NOS
                IF(J.EQ.M) THEN
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &            ENE_ORB(J)*T2(M,I,ia,ib)
                ENDIF
              ENDDO

              DO M = 1, NOS
                IF(I.EQ.M) THEN
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) -
     &            ENE_ORB(I)*T2(M,J,ia,ib)
                ENDIF
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &            0.5*V(ia,ib,ie,if)*T2(I,J,ie,if)
                ENDDO
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &            0.5*V(M,N,I,J)*T2(M,N,ia,ib)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) -
     &            V(M,ia,ie,J)*T2(M,I,ie,ib)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &            V(M,ia,ie,I)*T2(M,J,ie,ib)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) +
     &            V(M,ib,ie,J)*T2(M,I,ie,ia)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(I,J,ia,ib) = VT(I,J,ia,ib) -
     &            V(M,ib,ie,I)*T2(M,J,ie,ia)
                ENDDO
              ENDDO

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
! ***| THERE IS NO CORRESPONDING TERM

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
! ***| THERE IS NO CORRESPONDING TERM

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

      !   DO IA = NOS+1, NOS+NVS
      !     DO IB = NOS+1, NOS+NVS
      !       DO I = 1, NOS
      !         DO J = 1, NOS
      !           WRITE(6,*)V(ia,ib,I,J), VT(ia,ib,I,J), V(ia,ib,I,J)-VT(ia,ib,I,J)
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG
!  END  12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12-12

! START 13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13
! ***| VT_{IJ}^{AB}
      DO I = 1, NOS
        DO J = 1, NOS
          DO IA = NOS+1, NOS+NVS
            DO IB = NOS+1, NOS+NVS
! ***| H_N
              VT(I,J,ia,ib) = VT(I,J,ia,ib) + V(I,J,ia,ib)

! ***| (H_N(T_EXT))_C,OPEN
              DO IE = NOS+1, NOS+NVS
                VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &          V(ie,J,ia,ib)*T1(I,ie)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                VT(ia,ib,I,J) = VT(ia,ib,I,J) -
     &          V(ie,I,ia,ib)*T1(J,ie)
              ENDDO

              DO M = 1, NOS
                VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &          V(I,J,M,ia)*T1(M,ib)
              ENDDO

              DO M = 1, NOS
                VT(ia,ib,I,J) = VT(ia,ib,I,J) -
     &          V(I,J,M,ib)*T1(M,ia)
              ENDDO

              DO IE = NOS+1, NOS+NVS
                IF(ia.EQ.ie) THEN
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &            ENE_ORB(ia)*T2(I,J,ie,ib)
                ENDIF
              ENDDO

              DO IE = NOS+1, NOS+NVS
                IF(ib.EQ.ie) THEN
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) -
     &            ENE_ORB(ib)*T2(I,J,ie,ia)
                ENDIF
              ENDDO

              DO M = 1, NOS
                IF(J.EQ.M) THEN
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &            ENE_ORB(J)*T2(M,I,ia,ib)
                ENDIF
              ENDDO

              DO M = 1, NOS
                IF(I.EQ.M) THEN
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) -
     &            ENE_ORB(I)*T2(M,J,ia,ib)
                ENDIF
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO IF = NOS+1, NOS+NVS
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &            0.5*V(ie,if,ia,ib)*T2(I,J,ie,if)
                ENDDO
              ENDDO

              DO M = 1, NOS
                DO N = 1, NOS
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &            0.5*V(I,J,M,N)*T2(M,N,ia,ib)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) -
     &            V(ie,J,M,ia)*T2(M,I,ie,ib)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &            V(ie,I,M,ia)*T2(M,J,ie,ib)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) +
     &            V(ie,J,M,ib)*T2(M,I,ie,ia)
                ENDDO
              ENDDO

              DO IE = NOS+1, NOS+NVS
                DO M = 1, NOS
                  VT(ia,ib,I,J) = VT(ia,ib,I,J) -
     &            V(ie,I,M,ib)*T2(M,J,ie,ia)
                ENDDO
              ENDDO

! ***| ((T^{DAGGER}_EXT)H_N)_C,OPEN
! ***| THERE IS NO CORRESPONDING TERM

! ***| 1/2((T^{DAGGER}_EXT) * (F_N (T_EXT))_C,OPEN))_C,OPEN
! ***| THERE IS NO CORRESPONDING TERM

! ***| 1/2(((T^{DAGGER}_EXT)F_N)_C,OPEN * (T_EXT))_C,OPEN
! ***| THERE IS NO CORRESPONDING TERM

            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       DO IA = NOS+1, NOS+NVS
      !         DO IB = NOS+1, NOS+NVS
      !           ZZZ=VT(I,J,ia,ib)-VT(ia,ib,I,J)
      !           IF(abs(ZZZ).GT.(1.0d-17)) THEN
      !             WRITE(6,*)VT(I,J,ia,ib),VT(ia,ib,I,J),ZZZ
      !           ENDIF
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'PROOF THAT V (COLUMN 1) CHANGED TO VT (COLUMN 2)'

      !   DO I = 1, NOS
      !     DO J = 1, NOS
      !       DO IA = NOS+1, NOS+NVS
      !         DO IB = NOS+1, NOS+NVS
      !           WRITE(6,*)V(I,J,ia,ib), VT(I,J,ia,ib), V(I,J,ia,ib)-VT(I,J,ia,ib)
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG
!  END  13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13-13

!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************

! ***| FINAL CHECK FOR UNITARY and ANTISYMMETRY ****
! ***| UNCOMMENTING THE FOLLOWING LOOPS WILL TEST THE UNITARY AND ANTISYMETRY
! ***| CHARACTER OF THE TRANSFORMED MATRIX ELEMENTS. IDEALLY, THE FOLLOWING OUGHT
! ***| TO BE ZERO OR NEAR ZERO. IT WILL PRINT OUT OTHERWISE. THE THRESHOLD IS SO
! ***| THAT IT PRINT OUT 'SOMETHING', BUT THIS WILL BE SMALL IF THERE IS NO PROBLEM.
! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   WRITE(6,*) 'ONE-BODY UNITARY CHECK'
      !   DO P = 1, NOS+NVS
      !     DO Q = 1, NOS+NVS
      !       ZZZ=HT(P,Q)-HT(Q,P)
      !       IF(abs(ZZZ).GT.(1.0d-16)) THEN
      !         WRITE(6,*)P,Q,ZZZ
      !       ENDIF
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'TWO-BODY UNITARY CHECK'

      !   DO P = 1, NOS+NVS
      !     DO Q = 1, NOS+NVS
      !       DO R = 1, NOS+NVS
      !         DO S = 1, NOS+NVS
      !           ZZZ=VT(P,Q,R,S)-VT(R,S,P,Q)
      !           IF(abs(ZZZ).GT.(1.0d-16)) THEN
      !             WRITE(6,*)P,Q,R,S,ZZZ
      !           ENDIF
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      !   WRITE(6,*) 'TWO-BODY ANTISYMMETRY CHECK'
      !   DO P = 1, NOS+NVS
      !     DO Q = 1, NOS+NVS
      !       DO R = 1, NOS+NVS
      !         DO S = 1, NOS+NVS
      !           ZZZ=VT(P,Q,R,S)+VT(Q,P,R,S)
      !           IF(abs(ZZZ).GT.(1.0d-16)) THEN
      !             WRITE(6,*)P,Q,R,S,ZZZ
      !           ENDIF
      !         ENDDO
      !       ENDDO
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG

!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
! ***| DEBUG
      ! IF(NODEZERO) THEN

      !   WRITE(6,*) 'DIAGONALS OF THE TRANSFORMED FOCK MATRIX'
      !   DO P = 1, NOS+NVS
      !       WRITE(6,*)P,HT(P,P)
      !   ENDDO

      ! ENDIF
! ***| END DEBUG
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************

      !  FOCKT=HT

       DO IA = NOS+1, NOS+NVS
         DO IB = NOS+1, NOS+NVS
           DO M = 1, NOS
             HT(ia,ib) = HT(ia,ib) - VT(M,ia,M,ib)
           ENDDO
         ENDDO
       ENDDO

       DO I = 1, NOS
         DO IA = NOS+1, NOS+NVS
           DO M = 1, NOS
             HT(I,ia) = HT(I,ia) - VT(M,I,M,ia)
           ENDDO
         ENDDO
       ENDDO

       DO IA = NOS+1, NOS+NVS
         DO I = 1, NOS
           DO M = 1, NOS
             HT(ia,I) = HT(ia,I) - VT(M,ia,M,I)
           ENDDO
         ENDDO
       ENDDO

       DO I = 1, NOS
         DO J = 1, NOS
           DO M = 1, NOS
             HT(I,J) = HT(I,J) - VT(M,I,M,J)
           ENDDO
         ENDDO
       ENDDO

! ***| DEBUG
! ***| THERE IS AN EASY TEST TO CHECK THAT HT(P,Q) = HT(P,Q) - VT(M,P,M,Q)
! ***| TERMS ABOVE ARE CORRECT (FOR THE MOST PART). MAKE ALL VIRTUAL ORBITALS ACTIVE.
! ***| THE ABOVE LINES THEN CORRESPOND TO FOCK MATRIX ELEMENTS MINUS THE EFFECTIVE
! ***| ONE-ELECTRON MATRIX ELEMENTS WHICH SHOULD GIVE THE TRUE ONE-BODY MATRIX
! ***| ELEMENTS (HT(P,Q)-H(P,Q) ~= 0)
      ! IF(NODEZERO) THEN

      !   WRITE(6,*) 'DIAGONALS OF THE BARE FOCK MATRIX'
      !   DO P = 1, NOS+NVS
      !       WRITE(6,*)P,ENE_ORB(P)
      !   ENDDO

      !   WRITE(6,*) ' HT (COLUMN 1) VS H (COLUMN 2) AND THE DIFFERENCE'
      !   DO P = 1, NOS+NVS
      !     DO Q = 1, NOS+NVS
      !       WRITE(6,*)P,Q,HT(P,Q),H(P,Q),HT(P,Q)-H(P,Q)
      !     ENDDO
      !   ENDDO

      ! ENDIF
! ***| END DEBUG



! ***|  _______   _______   ______  __    __  ________
! ***| |       \ |       \ |      \|  \  |  \|        \
! ***| | $$$$$$$\| $$$$$$$\ \$$$$$$| $$\ | $$ \$$$$$$$$
! ***| | $$__/ $$| $$__| $$  | $$  | $$$\| $$   | $$
! ***| | $$    $$| $$    $$  | $$  | $$$$\ $$   | $$
! ***| | $$$$$$$ | $$$$$$$\  | $$  | $$\$$ $$   | $$
! ***| | $$      | $$  | $$ _| $$_ | $$ \$$$$   | $$
! ***| | $$      | $$  | $$|   $$ \| $$  \$$$   | $$
! ***|  \$$       \$$   \$$ \$$$$$$ \$$   \$$    \$$

! ***|  THIS IS SET UP WITH ALL OCCUPIED ORBITALS BEING ACTIVE ALONG WITH
! ***|  A SUBSET OF VIRTUALS, WHICH IS INDICATED BY THE INPUT VARIABLE NACTV
! ***|  AND COUNTED STARTING AT THE LOWEST ENERGY UNOCCUPIED ORBITAL.

! ***| PRINTING ONLY HT(ALPHA,ALPHA)
! ***| PRINTING ONLY VT(ALPHA,BETA,ALPHA,BETA)

! ***| PRINTING THRESHOLD
      PTHRESH = 1.0d-10

! ***| SPIN ORBITAL TO ORBITAL TRANSFORM. USED WHEN PRINTING.

      Q = 0

      DO P = 1, NOAS
        Q=Q+1
        TRANSFORM(Q) = P
      ENDDO

      DO P = NOAS+1, NOS
        Q=Q+1
        TRANSFORM(Q) = P-NOAS
      ENDDO

      DO P = NOS+1, NOS+NVAS
        Q=Q+1
        TRANSFORM(Q) = P-NOS+NOAS
      ENDDO

      DO P = NOS+NVAS+1, NOS+NVS
        Q=Q+1
        TRANSFORM(Q) = P-NOS-NVAS+NOAS
      ENDDO

! ***|  ACTINDEX1 IS AN ARRAY TO TEST IF A SET OF INDICES {P,Q} CONTAIN
! ***|  AT LEAST ONE EXTERNAL(INACTIVE) INDEX. IF SO, THEN THAT SET OF
! ***|  INDICES IS 'EXTERNAL' AND ACTINDEX1=1, ELSE ACTINDEX1=0.

      ACTINDEX1 = 0

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          IF(P.GT.(NOS+NACTV).AND.P.LE.(NOS+NVAS))ACTINDEX1(P,Q)=1
          IF(Q.GT.(NOS+NACTV).AND.Q.LE.(NOS+NVAS))ACTINDEX1(P,Q)=1
          IF(P.GT.(NOS+NVAS+NACTV).AND.P.LE.(NOS+NVS))ACTINDEX1(P,Q)=1
          IF(Q.GT.(NOS+NVAS+NACTV).AND.Q.LE.(NOS+NVS))ACTINDEX1(P,Q)=1
        ENDDO
      ENDDO

      PRINTINDEX1 = 0

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
! ***| SKIP PRINTING IF P OR Q IS EXTERNAL
          IF(ACTINDEX1(P,Q).EQ.1)PRINTINDEX1(P,Q)=1
! ***| SKIP PRINTING IF P OR Q IS BETA
          IF(P.GT.NOAS.AND.P.LE.NOS)PRINTINDEX1(P,Q)=1
          IF(Q.GT.NOAS.AND.Q.LE.NOS)PRINTINDEX1(P,Q)=1
          IF(P.GT.NOS+NVAS)PRINTINDEX1(P,Q)=1
          IF(Q.GT.NOS+NVAS)PRINTINDEX1(P,Q)=1
        ENDDO
      ENDDO

!      IF(NODEZERO) THEN
!       WRITE(6,*)'TRANSFORMED ~FOCK~ MATRIX'
!       CALL UTIL_FLUSH(6)
!      ENDIF
!
!      DO P = 1, NOS+NVS
!        DO Q = 1, NOS+NVS
!          IF(PRINTINDEX1(P,Q).EQ.1)CYCLE
!          IF(abs(FOCKT(P,Q)).GT.PTHRESH)THEN
!             IF(NODEZERO) THEN
!               WRITE(6,'(2I5,3X,F18.10)')
!     &         TRANSFORM(P),TRANSFORM(Q),FOCKT(P,Q)
!               CALL UTIL_FLUSH(6)
!             ENDIF
!          ENDIF
!        ENDDO
!      ENDDO

      IF(NODEZERO) THEN
       WRITE(6,*)'begin_one_electron_integrals'
       CALL UTIL_FLUSH(6)
      ENDIF

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          IF(PRINTINDEX1(P,Q).EQ.1)CYCLE
          IF(abs(HT(P,Q)).GT.PTHRESH)THEN
            if(oprint_qa) then
             IF(NODEZERO) THEN
               WRITE(6,'(A,2I5,3X,F18.10)')
     &         '1-e int ',TRANSFORM(P),TRANSFORM(Q),HT(P,Q)
           else
               WRITE(6,'(2I5,3X,F18.10)')
     &         TRANSFORM(P),TRANSFORM(Q),HT(P,Q)
           endif
               CALL UTIL_FLUSH(6)
             ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF(NODEZERO) THEN
       WRITE(6,*)'end_one_electron_integrals'
       CALL UTIL_FLUSH(6)
      ENDIF

! ***|  ACTINDEX2 IS AN ARRAY TO TEST IF A SET OF INDICES {P,Q,R,S} CONTAIN
! ***|  AT LEAST ONE EXTERNAL(INACTIVE) INDEX. IF SO, THEN THAT SET OF
! ***|  INDICES IS 'EXTERNAL' AND ACTINDEX2=1, ELSE ACTINDEX2=0.

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          DO R = 1, NOS+NVS
            DO S = 1, NOS+NVS
      actindex2(p,q,r,s)=(P.GT.(NOS+NACTV).AND.P.LE.(NOS+NVAS)).or.
     O              (Q.GT.(NOS+NACTV).AND.Q.LE.(NOS+NVAS)).or.
     O              (R.GT.(NOS+NACTV).AND.R.LE.(NOS+NVAS)).or.
     O              (S.GT.(NOS+NACTV).AND.S.LE.(NOS+NVAS)).or.
     O              (P.GT.(NOS+NVAS+NACTV).AND.P.LE.(NOS+NVS)).or.
     O              (Q.GT.(NOS+NVAS+NACTV).AND.Q.LE.(NOS+NVS)).or.
     O              (R.GT.(NOS+NVAS+NACTV).AND.R.LE.(NOS+NVS)).or.
     O              (S.GT.(NOS+NVAS+NACTV).AND.S.LE.(NOS+NVS))
            ENDDO
          ENDDO
        ENDDO
      ENDDO


      DO P = 1, NOS+NVS
        DO R = 1, NOS+NVS
          DO Q = 1, NOS+NVS
            DO S = 1, NOS+NVS
!     ***| SKIP PRINTING IF P, Q, R, OR S IS EXTERNAL
               printindex2(p,q,r,s)=ACTINDEX2(P,Q,R,S).or.
! ***| SKIP PRINTING IF P OR R IS BETA
     O              (P.GT.NOAS.AND.P.LE.NOS).or.
     O              (R.GT.NOAS.AND.R.LE.NOS).or.
     O              (P.GT.NOS+NVAS).or.
     O              (R.GT.NOS+NVAS).or.
! ***| SKIP PRINTING IF Q OR S IS ALPHA
     O              (Q.LE.NOAS).or.
     O              (S.LE.NOAS).or.
     O              (Q.GT.NOS.AND.Q.LE.NOS+NVAS).or.
     O              (S.GT.NOS.AND.S.LE.NOS+NVAS)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

! ***| THE ORDER OF INDICES IS SPECIFIC QDK YAML FILE

      IF(NODEZERO) THEN
        WRITE(6,*)'begin_two_electron_integrals'
        CALL UTIL_FLUSH(6)
      ENDIF

      DO P = 1, NOS+NVS
        DO Q = 1, NOS+NVS
          DO R = 1, NOS+NVS
            DO S = 1, NOS+NVS
              IF(PRINTINDEX2(P,Q,R,S))CYCLE
              IF(abs(VT(P,Q,R,S)).GT.PTHRESH)THEN
                 IF(NODEZERO) THEN
                   WRITE(6,'(4I5,2X,F18.10)')TRANSFORM(P),
     &             TRANSFORM(R),TRANSFORM(Q),TRANSFORM(S),
     &             VT(P,Q,R,S)
                   CALL UTIL_FLUSH(6)
                 ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(NODEZERO) THEN
        WRITE(6,*)'end_two_electron_integrals'
        CALL UTIL_FLUSH(6)
      ENDIF






! ! ***| NICK'S PRINTING FOR CCSDTQ
!       IF(NODEZERO) THEN
!         WRITE(6,*)'onebody.inp'
!         CALL UTIL_FLUSH(6)
!       ENDIF
!       I=0
!       DO P = 1, NOS+NVS
!         DO Q = 1, P
!           IF(PRINTINDEX1(P,Q).EQ.1)CYCLE
!           I=I+1
!              IF(NODEZERO) THEN
!                WRITE(6,'(F18.10,3X,I3)')
!      &         HT(P,Q),I
!                CALL UTIL_FLUSH(6)
!              ENDIF
!         ENDDO
!       ENDDO

!       IF(NODEZERO) THEN
!         WRITE(6,*)'twobody.inp'
!         CALL UTIL_FLUSH(6)
!       ENDIF
!       DO P = 1, NOS+NVS
!           DO R = 1, NOS+NVS
!         DO Q = 1, NOS+NVS
!             DO S = 1, NOS+NVS
!               IF(PRINTINDEX2(P,Q,R,S))CYCLE
!               IF(abs(VT(P,Q,R,S)).GT.PTHRESH)THEN
!                  IF(NODEZERO) THEN
!                    WRITE(6,'(4I5,2X,F18.10)')TRANSFORM(P),
!      &             TRANSFORM(R),TRANSFORM(Q),TRANSFORM(S),VT(P,Q,R,S)
!                    CALL UTIL_FLUSH(6)
!                  ENDIF
!               ENDIF
!             ENDDO
!           ENDDO
!         ENDDO
!       ENDDO

! ! ***| PRINTING FOR DMRG
!       IF(NODEZERO) THEN
!         WRITE(6,*)'fcidump'
!         CALL UTIL_FLUSH(6)
!       ENDIF
!       DO P = 1, NOS+NVS
!         DO R = 1, NOS+NVS
!           DO Q = 1, NOS+NVS
!             DO S = 1, NOS+NVS
!               IF(PRINTINDEX2(P,Q,R,S))CYCLE
!               IF(abs(VT(P,Q,R,S)).GT.PTHRESH)THEN
!                  IF(NODEZERO) THEN
!                    WRITE(6,'(F18.10,2X,4I4)')VT(P,Q,R,S),TRANSFORM(P),
!      &             TRANSFORM(R),TRANSFORM(Q),TRANSFORM(S)
!                    CALL UTIL_FLUSH(6)
!                  ENDIF
!               ENDIF
!             ENDDO
!           ENDDO
!         ENDDO
!       ENDDO
!       DO P = 1, NOS+NVS
!         DO Q = 1, NOS+NVS
!           IF(PRINTINDEX1(P,Q).EQ.1)CYCLE
!           IF(abs(HT(P,Q)).GT.PTHRESH)THEN
!              IF(NODEZERO) THEN
!                WRITE(6,'(F18.10,2X,4I4)')
!      &         HT(P,Q),TRANSFORM(P),TRANSFORM(Q),0,0
!                CALL UTIL_FLUSH(6)
!              ENDIF
!           ENDIF
!         ENDDO
!       ENDDO
c 
      return
      end
c
c
c
c
c
c
c

