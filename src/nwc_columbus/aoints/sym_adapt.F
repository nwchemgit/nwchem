      SUBROUTINE sym_adapt(ibas)
        ! Routine to adapt the AOs to symmetrized orbitals
      USE nwc_sym
      IMPLICIT NONE
#include "errquit.fh"
!#include "nwc_const.fh"
!#include "sym_adapt.fh"
#include "sym.fh"
#include "geom.fh"
#include "bas.fh"
#include "basP.fh"

      ! IN
      INTEGER, INTENT(IN) :: ibas ! basis handle
      ! LOCAL
      INTEGER :: geom ! geometry handle
      INTEGER :: i,j,k,ii,jj,kk,mm ! counters
      INTEGER :: l(mxaqn),m(mxaqn),n(mxaqn)
      INTEGER :: Lvar, Mvar, Nvar
      INTEGER :: ang_max !maximum L for basis ibas
      INTEGER :: iclass(mxsym),otype
      INTEGER :: naorb, naorbd, nsorb, ivarb, norbpsh, ncent
      INTEGER :: jkb, cent
      INTEGER :: ireplm ! integer function
      INTEGER :: nprim, ngeno, spcrt ! dummies for bas_continfo
      INTEGER :: mult(0:7)
      INTEGER :: lorb, lex_least
      DOUBLE PRECISION :: chars(mxsym*mxsym)
      CHARACTER*8 :: zclass(mxsym)

      LOGICAL :: status


      ! Get geometry (geom)
      IF (.NOT. bas_geom(ibas,geom)) CALL errquit
     $     ('sym_adapt: bas_geom failed for ibas', ibas,
     &       BASIS_ERR)

      ! Get number of centers (ncent)
      IF (.NOT. geom_ncent(geom,ncent)) CALL errquit
     $     ('sym_adapt: geom_ncent failed for geom', geom,
     &       GEOM_ERR)
      WRITE(*,*)"C LB, ncent=",ncent
      
      ! Get highest angular momentum for basis (ang_max)
      status = bas_high_angular(ibas,ang_max)

      ! Get group name (grp_name) 
      CALL sym_group_name(geom,grp_name)
      IF (.NOT. sym_abelian_group(geom)) CALL errquit
     $     ('sym_adapt: Columbus can only process abelian point
     & groups',geom, GEOM_ERR)

      ! Get group info A
      status = sym_char_table(grp_name,nop,nir,iclass,zir,zclass,chars)
      WRITE(*,*)"C LB, nop =",nop

      ! Get group info B
      CALL sym_bas_irreps(ibas,.true.,nbpsy(1:nir))

      ! Get group info C
      oprint=.TRUE.
      CALL sym_abelian_bas_info(ibas,
     $     sh_op_map, bf_op_map, bf_phase,
     $     sh_n_uq_op, bf_n_uq_op,
     $     sh_uq_op, bf_uq_op,
     $     sh_n_uq, bf_n_uq,
     $     sh_uq, bf_uq,
     $     sh_nat, bf_nat, sh_uq_bf,
     $     char_tab,
     $     bf_per_ir,
     $     bf_per_ir_cum,
     $     bf_so_ir,
     $     zir,
     $     so_uq_to_lab, so_lab1, so_lab2, oprint)

!      WRITE(*,'(a)')"C LB apply operator to shell"
!      WRITE(*,'(a)')"C LB sh_op_map ="
!      WRITE(*,'(8i4)')sh_op_map(:,1:6)
!      WRITE(*,'(a)')"C LB apply operator to basis function"
!      WRITE(*,'(a)')"C LB bf_op_map ="
!      WRITE(*,'(8i4)')bf_op_map(:,1:10)
!      WRITE(*,'(a)')"C LB apply operator to basis function (phase)"
!      WRITE(*,'(a)')"C LB bf_phase ="
!      WRITE(*,'(8f6.1)')bf_phase(:,1:10)
!      WRITE(*,'(a)')"C LB number of operations that produce unique
!     & shells"
!      WRITE(*,'(a)')"C LB sh_n_uq_op ="
!      WRITE(*,'(6i4)')sh_n_uq_op(1:6)
!      WRITE(*,'(a)')"C LB number of operations that produce unique
!     & basis functions"
!      WRITE(*,'(a)')"C LB bf_n_uq_op ="
!      WRITE(*,'(10i4)')bf_n_uq_op(1:10)
!      WRITE(*,'(a)')"C LB list of operations that produce unique
!     & shells"
!      WRITE(*,'(a)')"C LB sh_uq_op ="
!      WRITE(*,'(8i8)')sh_uq_op(:,1:6)
!      WRITE(*,'(a)')"C LB list of operations that produce unique
!     & bfns"
!      WRITE(*,'(a)')"C LB bf_uq_op ="
!      WRITE(*,'(8i4)')bf_uq_op(:,1:10)
!      WRITE(*,'(a)')"C LB number of unique shells"
!      WRITE(*,'(a)')"C LB sh_n_uq ="
!      WRITE(*,'(i4)')sh_n_uq
!      WRITE(*,'(a)')"C LB number of unique bfns"
!      WRITE(*,'(a)')"C LB bf_n_uq ="
!      WRITE(*,'(i4)')bf_n_uq
!      WRITE(*,'(a)')"C LB lexically highest sym-related shell"
!      WRITE(*,'(a)')"C LB sh_uq ="
!      WRITE(*,'(6i4)')sh_uq(1:6)
!      WRITE(*,'(a)')"C LB lexically highest sym-related bfn"
!      WRITE(*,'(a)')"C LB bf_uq ="
!      WRITE(*,'(10i4)')bf_uq(1:10)
!      WRITE(*,'(a)')"C LB map shells from nwchem to natural AO order"
!      WRITE(*,'(a)')"C LB sh_nat1 ="
!      WRITE(*,'(6i4)')sh_nat(1,1:6)
!      WRITE(*,'(a)')"C LB map shells from natural to nwchem AO order"
!      WRITE(*,'(a)')"C LB sh_nat2 ="
!      WRITE(*,'(6i4)')sh_nat(2,1:6)
!      WRITE(*,'(a)')"C LB map bfn from nwchem to natural AO order"
!      WRITE(*,'(a)')"C LB bf_nat ="
!      WRITE(*,'(10i4)')bf_nat(1,1:10)
!      WRITE(*,'(a)')"C LB map bfn from natural to nwchem AO order"
!      WRITE(*,'(a)')"C LB bf_nat ="
!      WRITE(*,'(10i4)')bf_nat(2,1:10)
!      WRITE(*,'(a)')"1st natural bfn from unique shell"
!      WRITE(*,'(a)')"C LB sh_uq_bf ="
!      WRITE(*,'(6i4)')sh_uq_bf(1,1:6)
!      WRITE(*,'(a)')"last natural bfn from unique shell"
!      WRITE(*,'(a)')"C LB sh_uq_bf ="
!      WRITE(*,'(6i4)')sh_uq_bf(2,1:6)
!      WRITE(*,'(a)')"1st unique bfn from unique shell"
!      WRITE(*,'(a)')"C LB sh_uq_bf ="
!      WRITE(*,'(6i4)')sh_uq_bf(3,1:6)
!      WRITE(*,'(a)')"last unique bfn from unique shell"
!      WRITE(*,'(a)')"C LB sh_uq_bf ="
!      WRITE(*,'(6i4)')sh_uq_bf(4,1:6)
!      WRITE(*,'(a)')"C LB bf_per_ir ="
!      WRITE(*,'(8i4)')bf_per_ir(:)
!      WRITE(*,'(a)')"list of irreps generated from unique bfn"
!      WRITE(*,'(a)')"C LB bf_so_ir ="
!      WRITE(*,'(8i4)')bf_so_ir(:,1:10)
!      WRITE(*,'(a)')"C LB so_uq_to_lab 1 ="
!      WRITE(*,'(8i12)')so_uq_to_lab(1,:,1:10)
!      WRITE(*,'(a)')"C LB so_uq_to_lab 2 ="
!      WRITE(*,'(8i12)')so_uq_to_lab(2,:,1:10)
!      WRITE(*,'(a)')"irrep of given SO"
!      WRITE(*,'(a)')"C LB so_lab1 ="
!      WRITE(*,'(10i4)')so_lab1(1,1:10)
!      WRITE(*,'(a)')"map twixt labelling schemes"
!      WRITE(*,'(a)')"C LB so_lab1 ="
!      WRITE(*,'(10i4)')so_lab1(2,1:10)
!      WRITE(*,'(a)')"unique bfn that is generator"
!      WRITE(*,'(a)')"C LB so_lab1 ="
!      WRITE(*,'(10i4)')so_lab1(3,1:10)
!      WRITE(*,'(a)')"map from sym-blocked SO order to natural SO order"
!      WRITE(*,'(a)')"C LB so_lab2 ="
!      WRITE(*,'(10i4)')so_lab2(1:10)
!      WRITE(*,*)"C LB, before isymax, grp_name=",grp_name

      ! number of generating operators
      ngen = INT(LOG(DBLE(nop))/LOG(2.0))

      ! Initialize isymax, reinitialize irrep and operator maps 
      ! column 1 of isymax is irreps of  x, y, z
      ! column 2 of isymax is irreps of Rx,Ry,Rz
      ! Note that for the math here, the irreps are in DALTON order!  
      ! This will be fixed with the irrep and operator maps.
      !
      ! Group Irrep order                     Operator Order
      ! ----- ------------------------------  ---------------------------
      ! C1    A                               E
      ! Cs    A'  A''                         E Oxy
      ! Ci    Ag  Au                          E i
      ! C2    A   B                           E C2z
      ! D2    A   B2  B1  B3                  E C2y C2x C2z
      ! C2v   A1  B1  B2  A2                  E Oyz Oxz C2z
      ! C2h   Ag  Bu  Au  Bg                  E Oxy C2z i
      ! D2h   Ag  B3u B2u B1g B1u B2g B3g Au  E Oyz Oxz C2z Oxy C2y C2x i
      SELECT CASE (grp_name)
        CASE ('C1')
          ISYMAX(1,1) = 0
          ISYMAX(2,1) = 0
          ISYMAX(3,1) = 0
          ISYMAX(1,2) = 0
          ISYMAX(2,2) = 0
          ISYMAX(3,2) = 0

        CASE ('Cs')
          ISYMAX(1,1) = 0
          ISYMAX(2,1) = 0
          ISYMAX(3,1) = 1
          ISYMAX(1,2) = 1
          ISYMAX(2,2) = 1
          ISYMAX(3,2) = 0

        CASE ('Ci')
          ISYMAX(1,1) = 1
          ISYMAX(2,1) = 1
          ISYMAX(3,1) = 1
          ISYMAX(1,2) = 0
          ISYMAX(2,2) = 0
          ISYMAX(3,2) = 0

        CASE ('C2')
          ISYMAX(1,1) = 1
          ISYMAX(2,1) = 1
          ISYMAX(3,1) = 0
          ISYMAX(1,2) = 1
          ISYMAX(2,2) = 1
          ISYMAX(3,2) = 0

        CASE ('D2')
          ISYMAX(1,1) = 3
          ISYMAX(2,1) = 1
          ISYMAX(3,1) = 2
          ISYMAX(1,2) = 3
          ISYMAX(2,2) = 1
          ISYMAX(3,2) = 2

          irrepmap(1)=2
          irrepmap(2)=1
          irrepmap(3)=3

        CASE ('C2v')
          ISYMAX(1,1) = 1
          ISYMAX(2,1) = 2
          ISYMAX(3,1) = 0
          ISYMAX(1,2) = 2
          ISYMAX(2,2) = 1
          ISYMAX(3,2) = 3

          irrepmap(1)=3
          irrepmap(2)=1
          irrepmap(3)=2

          opmap(1)=3
          opmap(2)=2
          opmap(3)=1

        CASE ('C2h')
          ISYMAX(1,1) = 1
          ISYMAX(2,1) = 1
          ISYMAX(3,1) = 2
          ISYMAX(1,2) = 3
          ISYMAX(2,2) = 3
          ISYMAX(3,2) = 0

          irrepmap(1)=2
          irrepmap(2)=3
          irrepmap(3)=1

          opmap(1)=1
          opmap(2)=3
          opmap(3)=2

        CASE ('D2h')
          ISYMAX(1,1) = 1
          ISYMAX(2,1) = 2
          ISYMAX(3,1) = 4
          ISYMAX(1,2) = 6
          ISYMAX(2,2) = 5
          ISYMAX(3,2) = 3

          irrepmap(1)=7
          irrepmap(2)=3
          irrepmap(3)=4
          irrepmap(4)=5
          irrepmap(5)=2
          irrepmap(6)=6
          irrepmap(7)=1
           
          opmap(1)=3
          opmap(2)=5
          opmap(3)=6
          opmap(4)=7
          opmap(5)=4
          opmap(6)=2
          opmap(7)=1

        CASE DEFAULT
          CALL errquit ('sym_adapt: symmetry not recognized',
     &          grp_name, GEOM_ERR)
      END SELECT
!      WRITE(*,*)"C LB, isymax="
!      WRITE(*,'(3i4)')isymax(1:3,1:2)

      ! Initialize the stabilizer list
      CALL stab_cent(geom,ncent)

!      ! initialize mult
!      DO i = 0,7
!        mult(i)=2**MAX(0,ngen-is(i))
!      ENDDO
!      WRITE(*,*)"C LB, ibas =",ibas
!      WRITE(*,*)"C LB, bas_is_spherical =",bas_is_spherical(ibas)
!      WRITE(*,*)"C LB, ang_max=",ang_max
!
      ! Initialize isymao
      IF (.NOT. bas_is_spherical(ibas)) THEN
        ! ----- using cartesian basis -----
        DO i = 1, ang_max + 1
          CALL lmnval(i,i*(i + 1)/2,l,m,n)
            DO j = 1, i*(i + 1)/2
                Lvar = MOD(l(j),2)*isymax(1,1)
                Mvar = MOD(m(j),2)*isymax(2,1)
                Nvar = MOD(n(j),2)*isymax(3,1)
                isymao(i,j) = IEOR(Lvar,IEOR(Mvar,Nvar))
            ENDDO
        ENDDO

      ELSE
        ! ----- using spherical basis -----
        DO i = 0, ang_max
          ii = i + 1
          IF (i .EQ. 0) THEN
            isymao(ii,1) = ireplm(0,0)
          ELSEIF (i .EQ. 1) THEN
            isymao(ii,1) = ireplm(1, 1)
            isymao(ii,2) = ireplm(1,-1)
            isymao(ii,3) = ireplm(1, 0)
          ELSE
            DO j = -i, i
              jj = j + i + 1
              isymao(ii,jj) = ireplm(i,j)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
       
      ! ----- Calculate transformation matrix from AO to SO -----

      nsorb = 0 ! number of symmetry orbitals

      ! --- Loop over irreps
      DO ii = 0, nir-1
      i=irrepmap(ii)
!      WRITE(*,'(a,1i4)')"  C LB, irrep= ",i
        naorbd = 0

        ! --- Loop over unique shells
        DO j = 1, sh_n_uq
          jj=sh_uq(j)

          IF (.NOT. bas_cn2ce(ibas,jj,cent)) CALL errquit
     &         ("sym_adapt: can't find center for bfn",jj,GEOM_ERR)
          status = bas_continfo(ibas,jj,otype,nprim,ngeno,spcrt)
          norbpsh = sh_uq_bf(4,j)-sh_uq_bf(3,j)+1

          ! --- Loop over orbitals
          DO k = 1, norbpsh
          naorbd = naorbd + 1
!      WRITE(*,'(a,2i4)')"      C LB, i=",i
!      WRITE(*,'(a,2i4)')"      C LB, shell,orbital=",jj,k
!      WRITE(*,'(a,2i4)')"      C LB, otype+1=",otype+1
!      WRITE(*,'(a,2i4)')"      C LB, isymao=",isymao(otype+1,k)
         !   naorb = naorb + 1
            ivarb = IEOR(i,isymao(otype+1,k))
!      WRITE(*,'(a,i4)')"ivarb =",ivarb
            IF (IAND(stabilizer(cent),ivarb).EQ.0) THEN
!      WRITE(*,*)"        C LB, orbital contributes!"
              nsorb = nsorb + 1
!      WRITE(*,'(a,1i4)')
!     &"         C LB, this will be symmetrized orbital",nsorb
              jkb = 0
!      WRITE(*,'(a,1i4)')"C LB, naorbd=",naorbd

              ! --- Loop over symmetry operations
              DO mm = 1, nop
                kk = opmap(mm - 1)
!      WRITE(*,'(a,2i4)')"      C LB, mm=",mm
!      WRITE(*,'(a,2i4)')"C LB, mm,bf_uq_op=",
!     $  mm,bf_uq_op(mm,bf_uq(naorbd))
                IF(bf_uq_op(mm,bf_uq(naorbd)) .EQ. 0) EXIT
                jkb = jkb + 1
!      WRITE(*,'(a,2i4)')"C LB, KK, IVARB =",kk,ivarb
!      WRITE(*,'(a,2i4)')"C LB, IAND =",iand(kk,ivarb)
!      WRITE(*,'(a,2i4)')"C LB, parbit =",parbit(iand(kk,ivarb))

                lorb = lex_least(bf_uq(naorbd))
!      WRITE(*,'(a,1i4)')"C LB, aorb =",bf_uq(naorbd)
!      WRITE(*,'(a,1i4)')"C LB, lorb =",lorb
                itran(nsorb,jkb)=
     &      bf_op_map(kk+1,lorb)*parbit(IAND(mm-1,ivarb))
!     &      bf_op_map(mm,bf_uq(naorbd))*parbit(IAND(kk,ivarb))
!     &      bf_op_map(mm,lorb)*bf_phase(mm,lorb)*parbit(IAND(kk,ivarb))
!      WRITE(*,*)"C LB, itran =",itran(nsorb,jkb)
              ENDDO !mm (operators)
              wt(nsorb,1)=nsorb
              wt(nsorb,2)=jkb

      !       ELSE
      !         naorbd = naorbd + mult(stabilizer(cent))
            ENDIF
          ENDDO !k (orbitals)
        ENDDO !j (unique shells)
      ENDDO !ii (irreps)

      WRITE(*,'(a)')" "
      WRITE(*,'(a)')"Symmetrized Orbitals"
      WRITE(*,'(a)')"--------------------"
      WRITE(*,'(a)')"          Sym Orb       Atomic Orbs"
      k = 0 
      DO j = 1, nir
        WRITE(*,'(a)')" "
        WRITE(*,'(a,a)')"Group ",zir(j)
        DO i = 1, nbpsy(j)
          k = k + 1
          WRITE(*,'(a,1i4,a,8i4)')'          ',k,'     ',
     &          itran(k,1:wt(k,2))
        ENDDO ! i (sym orbitals)
      ENDDO ! j (irreps)

      RETURN
      END !sym_adapt




!-----------------------------------------------------------------      
      INTEGER FUNCTION lex_least(orb)
      ! Find the lexically lowest equivalent orbital

      USE nwc_sym, ONLY : mxbf, mxshel, bf_uq, bf_op_map, nop
      IMPLICIT NONE

      ! IN
      INTEGER, INTENT(IN) :: orb ! orbital in question

      ! LOCAL
      INTEGER :: k ! counter
      INTEGER :: lorb !marks where lowest eq orb is in the map

      lorb = 1
      DO k = 1, nop
        IF (bf_op_map(k,orb) .EQ. 0) EXIT
        IF (bf_op_map(k,orb) .LT. bf_op_map(lorb,orb)) THEN
          lorb = k
        ENDIF
      ENDDO ! k
      
      lex_least = bf_op_map(lorb,orb)

      RETURN
      END !FUNCTION lex_least




!------------------------------------------------------------------
      INTEGER FUNCTION ireplm(L,M)
      ! Taken from DALTON (hersol.F)
      ! Symmetry of RLM integrals for spherical basis

      USE nwc_sym, ONLY: isymax
      IMPLICIT NONE
!#include "nwc_const.fh"
!#include "sym_adapt.fh"

      ! IN
      INTEGER, INTENT(IN) :: L, M

      ireplm = 0
      IF (MOD(L + M,2) .EQ. 1)    ireplm = isymax(3,1)
      IF (MOD(ABS(M),2) .EQ. 1) ireplm = IEOR(ireplm,isymax(1,1))
      IF (M.LT.0) ireplm=IEOR(ireplm,IEOR(isymax(1,1),isymax(2,1)))

      RETURN
      END !FUNCTION ireplm





!------------------------------------------------------------------
      SUBROUTINE lmnval(nhkta,khkta,l,m,n)
      ! Taken from DALTON (hergam.F and carpow.F)
      ! determine l, m, n for cartesian basis

      USE nwc_sym, ONLY : mxaqnm, mxqnm
      IMPLICIT NONE
!#include "nwc_const.fh"
!#include "sym_adapt.fh"

      ! IN
      INTEGER, INTENT(IN) :: nhkta, khkta
      ! OUT
      INTEGER, INTENT(OUT) :: l(khkta),m(khkta),n(khkta)      
      ! LOCAL
      INTEGER :: icomp,istep(mxaqnm),mval(mxaqnm),nval(mxaqnm)
      INTEGER :: i,j,ij ! counters

      ! Calculate cartesian powers
      ij=0
      DO i=1, mxqnm
        DO j=1, i
          ij=ij + 1
          istep(ij) = i
          mval(ij)  = i - j
          nval(ij)  = j - 1
        ENDDO
      ENDDO

      icomp = 0
      DO i = 1, khkta
        icomp = icomp + 1
        l(icomp) = nhkta - istep(i)
        m(icomp) = mval(i)
        n(icomp) = nval(i)
      ENDDO

      RETURN
      END ! lmnval



!------------------------------------------------------------------
      SUBROUTINE stab_cent(geom,ncent)
      ! Taken from DALTON (herrdn.F)
      ! Map stabilizers of each atom center

      USE nwc_sym, ONLY: ngen, isymax, stabilizer, coords
      IMPLICIT NONE
!#include "nwc_const.fh"
!#include "sym_adapt.fh"

      ! IN
      INTEGER, INTENT(IN) :: ncent  ! number of unique atom centers
      INTEGER, INTENT(IN) :: geom   ! geom handle

      ! LOCAL
      INTEGER :: i,l,ll,m ! counters
      INTEGER :: mul



      ! --- Loop over atoms
      DO i = 1, ncent
        mul = 0
        ll   = 1

        ! --- Loop over generating operators
        DO 10 l = 1, ngen
!      WRITE(*,*)"C LB, L=",L
          DO m = 1, 3
!      WRITE(*,*)"C LB, M=",M
!      WRITE(*,*)"C LB, LL=",LL
!      WRITE(*,*)"C LB, isymax=",isymax(m,1)
            IF(IAND(ll,isymax(m,1)).NE.0) THEN
!      WRITE(*,*)"C LB, i,coords=",i,coords(m,i,geom)
              IF(ABS(coords(m,i,geom)).GE.1D-6)
     &  GOTO 10
            ENDIF
          ENDDO
          mul = mul + ll
!      WRITE(*,*)"C LB, mul=",mul
10      ll = 2 * ll
        stabilizer(i) = mul
!      WRITE(*,*)"C LB, i, stabilizer=",i,stabilizer(i)
      ENDDO

      RETURN
      END ! stab_cent




!-----------------------------------------------------------------      
      SUBROUTINE sym_1int(num,numints,intbuf,labbuf,val,lab)
      ! Calculate the SOs from the AOs

      USE nwc_sym, ONLY: nir, bf_per_ir_cum, nbpsy, wt, itran
      IMPLICIT NONE
!#include "nwc_const.fh"
!#include "sym_adapt.fh"

      ! IN
      INTEGER :: num  !# of values to place in packed buffer
      DOUBLE PRECISION :: intbuf(1:num) ! AO integer buffer
      INTEGER :: labbuf(1:2,1:num) ! AO label buffer

      ! LOCAL
      INTEGER :: i,j,k,l,m !counters
      INTEGER :: offset
      INTEGER :: i_lab,j_lab ! SO integral labels
      INTEGER :: ao_i_lab,ao_j_lab ! AO integral labels
      INTEGER :: sgni, sgnj ! sign of AO

      ! OUT
      INTEGER :: numints ! number of SO ints 
      DOUBLE PRECISION :: val(1:num) ! SO integrals
      INTEGER :: lab(1:2,1:num) ! SO labels

      WRITE(*,*)"C LB, now in sym_1int"
      numints = 0

      ! irreps
      DO i = 1, nir
      WRITE(*,'(a,2i4)')"C LB, irrep =",i
        offset = bf_per_ir_cum(i-1)
      WRITE(*,'(a,2i4)')"C LB, offset =",offset

        ! SO integral labels
        DO i_lab = offset+1, offset+nbpsy(i)
          DO j_lab = i_lab,offset+nbpsy(i)
!      WRITE(*,'(a,2i4)')"  C LB, SO integral label = ",i_lab,j_lab
            numints = numints + 1
            val(numints) = 0.0D0
            lab(1,numints) = i_lab
            lab(2,numints) = j_lab

            ! AO linear combination for i_lab
            DO j = 1, wt(i_lab,2)
              ao_i_lab = itran(i_lab,j)
              sgni = SIGN(1,ao_i_lab)
              ao_i_lab = ABS(ao_i_lab)

              ! AO linear combination for j_lab
              DO k = 1, wt(j_lab,2)
                ao_j_lab = itran(j_lab,k)
                sgnj = SIGN(1,ao_j_lab)
                ao_j_lab = ABS(ao_j_lab)
!      WRITE(*,'(a,2i4)')"     C LB, AO integral label = ",
!     &                          ao_i_lab,ao_j_lab

                ! Find the AO integral that corresponds to ao_i_lab
                ! ao_j_lab, multiply by sgni*sgnj, and add to val
                DO l = 1, num
!      WRITE(*,'(a,2i4)')"       C LB, labbuf ="
!     &                   ,labbuf(1:2,l)
                  IF (((labbuf(1,l) .EQ. ao_i_lab) .AND.
     &                 (labbuf(2,l) .EQ. ao_j_lab)) .OR.
     &                ((labbuf(1,l) .EQ. ao_j_lab) .AND.
     &                 (labbuf(2,l) .EQ. ao_i_lab))) THEN
                    ! match found!
      WRITE(*,'(a,1e15.3)')"     C LB, integral found:",
     &                               intbuf(l)
                    val(numints) = val(numints) + intbuf(l)*sgni*sgnj
                  ENDIF


                ENDDO ! l
              ENDDO !k
            ENDDO !j
!      WRITE(*,'(a,3i4,1e15.3)')"  C LB, SO integral:",
!     &                 numints,lab(1:2,numints),val(numints)
          ENDDO !j_lab (SO lab j)
        ENDDO !i_lab (SO lab i)
      ENDDO !i (irreps)
      RETURN
      END ! sym_1int
