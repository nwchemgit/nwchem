*
* $Id$
*
      subroutine geom_bandbi(geom)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"
      integer geom
c     
c     New version of bandbi
c
c     For the geometry referred to by the geometry handle form the
c     B, Binv=(G^-1)B, P=(G^-1)G, C
c
c     B and Binv are stored as (ncart,nzvar)
c     P is stored as (nzvar,nzvar)
c     C is stored as (nzvar,nzvar)
c
c     C is the matrix that applies constraints/constants to a
c     step or gradient ... C*dq = constrained gradient.
c
c     B and Binv refer to the UNconstrained redundant internal variables
c     P projects on the the constrained redundant internal variables
c     
      integer nzmat, nzvar, nvar, nat, ndim, ncart, nconstraint
      logical status
      integer k_x, k_y, k_z, k_c, k_e, k_iz
      integer l_x, l_y, l_z, l_c, l_e, l_iz
c     
      logical geom_hnd_check_data, geom_zmt_get_nzvar,
     $     geom_zmt_get_nizmat
      external geom_hnd_check_data, geom_zmt_get_nzvar,
     $     geom_zmt_get_nizmat
c     
      status = geom_ncent(geom, nat)
      status = geom_zmt_get_nzvar(geom, nzvar) .and. status
      status = geom_zmt_get_nizmat(geom, nzmat) .and. status
      if (.not. status) call errquit
     $     ('geom_bandbi: failed getting info ', geom, GEOM_ERR)
c     
c     ----- return if no internal coordinates -----
c     
      if(nzmat.eq.0) return
c
      if (nat .eq. 1) then
         nvar = 0
      else if (nat .eq. 2) then
         nvar = 1
      else
         nvar = 3*nat - 6       ! Linear molecules not done here
      endif
c
      ncart = 3*nat
      ndim = max(ncart,nzvar)
c
c     Allocate 3 matrices all (ndim,ndim)
c     
      if(.not.ma_push_get(mt_dbl,ndim**2    ,'mem x ',l_x,k_x))
     &     call errquit('hnd_bandbi, malloc of x  failed',911,
     &       MA_ERR)
      if(.not.ma_push_get(mt_dbl,ndim**2        ,'mem y ',l_y,k_y))
     &     call errquit('hnd_bandbi, malloc of y  failed',911, MA_ERR)
      if(.not.ma_push_get(mt_dbl,ncart ,'c',l_c,k_c))
     &     call errquit('hnd_bandbi, malloc of c  failed',ncart, MA_ERR)
      if(.not.ma_push_get(mt_dbl,ndim ,'e',l_e,k_e))
     &     call errquit('hnd_bandbi, malloc of e  failed',nzvar, MA_ERR)
      if(.not.ma_push_get(mt_int,ndim**2   ,'mem z ',l_iz,k_iz))
     &     call errquit('hnd_bandbi, malloc of iz  failed',911, MA_ERR)
c
c     ----- calculate -c- matrix ----
c
      call geom_make_constraint_matrix(geom, nzvar, dbl_mb(k_x),
     $     nconstraint)
c     
c     ----- calculate -b- matrix -----
c     
      call hnd_bmat(geom,nzvar,ncart,
     1     dbl_mb(k_x),dbl_mb(k_y),int_mb(k_iz),
     2     nzmat,dbl_mb(k_c))

      if(.not.ma_chop_stack(l_iz))
     &     call errquit('geom_bandbi, ma_chop_stack of iz  failed',911,
     &       MA_ERR)

      if(.not.ma_push_get(mt_dbl,ndim**2   ,'mem z ',l_z,k_z))
     &     call errquit('hnd_bandbi, malloc of z  failed',911, MA_ERR)
c
c     compute b^-1, p
c     
      call geom_binvr(nvar,nzvar,ncart,ndim,
     $     dbl_mb(k_x),dbl_mb(k_y),dbl_mb(k_z),dbl_mb(k_e))
c     
c     ----- release memory pointers -----
c     
      if(.not.ma_chop_stack(l_x))
     &     call errquit('geom_bandbi, ma_chop_stack of i01  failed',911,
     &       MA_ERR)
c     
      end
      SUBROUTINE HND_BMAT(GEOM,NZVAR,NCART,B,ZMAT,IZMAT,NZMAT,C)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "geom.fh"
#include "util_params.fh"
      integer geom
C     
C     ----- CONSTRUCT THE B MATRIX -----
C     
      integer na,nzvar,ncart,nzmat
      PARAMETER (    NA=10 )
      LOGICAL     DBUG
      LOGICAL     OUT
      LOGICAL     GEOM_ZMT_GET_IZMAT
      LOGICAL     GEOM_ZMT_GET_NIZMAT
      CHARACTER*8 ERRMSG
      double precision B(NCART,*)
      double precision ZMAT(*)
      integer          IZMAT(*)
      double precision A(3,NA)
      double precision C(3,*)
      DIMENSION ERRMSG(3)
      double precision eqval,pi2,degree,bohr,zero,four,pi,pio2,rtod
      integer iw,iadd,i,itype,izvar,j,ndum
      DATA ERRMSG /'PROGRAM ','STOP IN ','- BMAT -'/
c     DATA PI2    /6.28318530717958D+00/
      DATA DEGREE /360.0D+00/
      DATA BOHR   /cau2ang/
      DATA ZERO   /0.0D+00/
      DATA FOUR   /4.0D+00/
C     
      IW = LUOUT
      DBUG  =.FALSE.
      OUT   =.FALSE. 
      OUT   =OUT.OR.DBUG
C     
      pi    =acos(-1.0d0)
      pi2   =2.0d0*pi
      PIO2  =PI2/FOUR
      RTOD  =PI2/DEGREE
C     
C     ----- READ -IZMAT- -----
C     -A-     AT -90-
C     -IZMAT- AT -96-
C     
      NDUM=3*NA
C     
      if (.not. geom_cart_coords_get(geom, c))
     $     call errquit('hnd_bmat: failed getting coords',0, GEOM_ERR)
      if(.not.geom_zmt_get_nizmat(geom, nzmat))
     1     call errquit('geom_bmat: geom_zmt_get_nizmat failed',0,
     &       GEOM_ERR)
      if(.not.geom_zmt_get_izmat(geom, izmat, nzmat))
     1     call errquit('geom_bmat: geom_zmt_get_izmat failed',0,
     &       GEOM_ERR)
C     
      IF(OUT) THEN
         WRITE(IW,9997) 
         WRITE(IW,9996) (IZMAT(I),I=1,NZMAT)
      ENDIF
C     
C     ----- ZERO OUT THE B MATRIX -----
C     
      DO J = 1,NZVAR
         DO I = 1,NCART
            B(I,J) = ZERO
         ENDDO
      ENDDO
C     
      IADD=1
      I   =0
      DO IZVAR=1,NZVAR
         I    =I+1            
         ITYPE=IZMAT(IADD)
         IF(OUT) THEN
            WRITE(IW,9999) IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART
         ENDIF
         IF(ITYPE.EQ.1) THEN
C     
C     ----- BOND STRETCH -----
C     
            CALL HND_BSTR(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1           C,B,NCART,BOHR)
            ZMAT(I) = EQVAL
            IADD = IADD + 3
         ELSEIF(ITYPE.EQ.2 .or. itype.eq.7) THEN
C     
C     ----- ANGLE BEND -----
C     
            CALL HND_BEND(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1           IZMAT(IADD+3),
     2           C,B,NCART,RTOD,BOHR)
            ZMAT(I) = EQVAL
            IADD = IADD + 4
            if (itype .eq. 7) iadd = iadd + 1 ! Additional orientation
         ELSEIF(ITYPE.EQ.3) THEN
C     
C     ----- TORSION -----
C     
            CALL HND_TORS(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1           IZMAT(IADD+3),IZMAT(IADD+4),
     2           C,B,NCART,RTOD,BOHR)
            ZMAT(I) = EQVAL
            IADD = IADD + 5
         ELSEIF(ITYPE.EQ.4) THEN
C     
C     ----- OUT OF PLANE ANGLE BEND -----
C     
            CALL HND_OPLA(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1           IZMAT(IADD+3),IZMAT(IADD+4),
     2           C,B,NCART,RTOD,PIO2,BOHR)
            ZMAT(I) = EQVAL
            IADD = IADD + 5
         ELSEIF(ITYPE.EQ.5) THEN
C     
C     ----- LINEAR ANGLE BEND -----
C     
            CALL HND_LIBE(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1           IZMAT(IADD+3),IZMAT(IADD+4),
     2           C,B,NCART,RTOD,A,NA,BOHR)
            ZMAT(I  ) = EQVAL
            ZMAT(I+1) = EQVAL
            I=I+1
            IADD=IADD+5
         ELSEIF(ITYPE.EQ.6) THEN
C     
C     ----- DIHEDRAL ANGLE BETWEEN TWO PLANES SHARING ONE ATOM -----
C     
            CALL HND_DIHPLA(EQVAL,I,IZMAT(IADD+1),C,B,NCART,RTOD)
            ZMAT(I)=EQVAL
            IADD=IADD+6
C     
         ELSE
            WRITE(IW,9995) IZVAR,IADD,ITYPE
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF    
C     
      ENDDO
C     
      IF(OUT) THEN
**         CALL HND_PRTZMT(NZVAR,IZMAT,ZMAT,DUM,0)
         WRITE(IW,9998)
         IF(DBUG) THEN
            CALL HND_PRSQ(B,NZVAR,NCART,NCART)
         ENDIF
      ENDIF
C     
C     ----- SAVE -B- AND -Z- ON -IDAF- -----
C     -A-                                 AT -90-
C     -B-  (pristine)                     AT -91-
C     -S-  (pristine)                     AT -92-
C     -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C     AT -94-
C     -ZMAT-                              AT -95-
C     -IZMAT-                             AT -96-
C     
      NDUM=NCART*NZVAR
C     
C     ----- STORE FOR -NWCHEM- -----
C     
      call geom_hnd_put_data('zmat',zmat,nzvar)
      call geom_hnd_put_data('b',b,ncart*nzvar)
C     
      RETURN
 9999 FORMAT(' IN -BMAT- IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART = ',
     1     6I5)
 9998 FORMAT(/,' -B- MATRIX (ROWS = CARTESIAN COORDINATE.',
     1     ' COLUMNS = INTERNAL COORDINATES.)')
 9997 FORMAT(' IN -BMAT- IZMAT = ')
 9996 FORMAT(12I5)
 9995 FORMAT(' ERROR IN -BMAT- , IZVAR,IADD,ITYPE = ',3I5)
      END
      subroutine geom_binvr(nvar,nzvar,ncart,ndim,
     $     x, y, z, e)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer nvar              ! [input] Expected #indep variables
      integer nzvar             ! [input] No. of (redundant) internals
      integer ncart             ! [input] Number of cartesians
      integer ndim              ! [input] max(ncart,nzvar)
      double precision x(ndim*ndim) ! [scratch]
      double precision y(ndim*ndim) ! [scratch]
      double precision z(ndim*ndim) ! [scratch]
      double precision e(ndim)  ! [scratch]
c
      integer info, i, j, ioff, nindep, ij, ji
      integer ii,jj
      logical odebug
      external dabssum
      double precision dabssum
c
c     Form generalized inverse of the B matrix
c
      odebug = (ga_nodeid().eq.0) .and. 
     $     util_print('geom_binvr',print_never)
      call geom_hnd_get_data('b', x, ncart*nzvar) ! B -> x
      if (odebug) then
         write(LuOut,*) ' B matrix '
         call output(x, 1, ncart, 1, nzvar, ncart, nzvar, 1)
      endif
      call dgemm('t', 'n', nzvar, nzvar, ncart, 1d0,  
     $     x, ncart, x, ncart, 0d0, y, nzvar) ! BT.B -> y ... should symmetrize
      do i = 1, nzvar
         do j = 1, nzvar
            ij = i + (j-1)*nzvar
            ji = j + (i-1)*nzvar
            y(ij) = 0.5d0*(y(ij) + y(ji))
            y(ji) = y(ij)
         end do
      end do
      call dsyev('v', 'u', nzvar, y, nzvar, e, z, ndim**2, info)
      if (info .ne. 0) call errquit('geom_binvr: dsyev failed', info,
     &       MEM_ERR)
      if (odebug) then
         write(LuOut,*) ' Eigenvalues of G '
         call doutput(e, 1, nzvar, 1, 1, nzvar, 1, 1)
         write(LuOut,*) ' Eigenvectors of G'
         call output(y, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
      endif
c
c     Now have B in x, evecs of G=BT.B in y, evals in e
c
c     Count the significant eigenvalues to make sure we have enough,
c     take the reciprocal of the eigenvalue, zero out small ones.
c
      nindep = 0
      do i = 1, nzvar
         if (abs(e(i)) .gt. 1d-4) then
            nindep = nindep + 1
            e(i) = 1d0 / e(i)
         else
            e(i) = 0.0d0
         endif
      enddo
      if (nindep .ne. nvar) then
         if (ga_nodeid() .eq. 0) then
            write(6,9990) nindep, nvar
 9990       format(' !! There are insufficient internal variables:',
     $           ' expected ', i5, ' got ', i5,/,
     $             ' !! Either AUTOZ failed or your geometry has',
     $           ' changed so much that the',/,
     $           ' !! coordinates should be regenerated.')
            call util_flush(6)
         end if
         call ga_sync()
         call errquit
     $     ('geom_binvr: #indep variables incorrect',100000*nindep+nvar,
     &       GEOM_ERR)
      endif
c
c     Form G- = VT (e^-1) V
c
      do i = 1, nzvar
         ioff = (i-1)*nzvar
         do j = 1, nzvar
            z(j+ioff) = y(j+ioff) * e(i)
         enddo
      enddo
      call dgemm('n', 't', nzvar, nzvar, nzvar, 1d0, 
     $     z, nzvar, y, nzvar, 0d0, x, nzvar)
c
c     Now have G^-1 in x.  Form G again and then P = (G^-1)G.
c
      call geom_hnd_get_data('b', y, ncart*nzvar) ! B -> y
      call dgemm('t', 'n', nzvar, nzvar, ncart, 1d0,  
     $     y, ncart, y, ncart, 0d0, z, nzvar) ! G = BT.B -> z
      call dgemm('n', 'n', nzvar, nzvar, nzvar, 1d0, 
     $     z, nzvar, x, nzvar, 0d0, y, nzvar) ! P = (G^-1)G -> y
      call geom_hnd_put_data('p', y, nzvar*nzvar)
      if (odebug) then
         write(LuOut,*) ' Projector P without constraints'
         call output(y, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
      endif
c
c     Finally B^-1 = G^-1 B
c
      call geom_hnd_get_data('b', y, ncart*nzvar) ! B -> y
      call dgemm('n', 'n', ncart, nzvar, nzvar, 1d0, 
     $     y, ncart, x, nzvar, 0d0, z, ncart) ! B^-1 = B.(G^-1) -> z
      call geom_hnd_put_data('b^-1', z, ncart*nzvar)
c
c     Verify inverse property
c
      if (odebug) then
         write(LuOut,*) ' B inverse'
         call output(z, 1, ncart, 1, nzvar, ncart, nzvar, 1)
         call dgemm('n', 't', ncart, ncart, nzvar, 1d0, 
     $        y, ncart, z, ncart, 0d0, x, ncart) ! B*B^-1
         write(LuOut,*) ' B*B^-1 '
         call output(x, 1, ncart, 1, ncart, ncart, ncart, 1)
      endif
c
c     Above made P without constraints ... now add constraints
c
c     P <- P - PC((CPC)^-1)CP
c
c     where C is the projector on to the constraints (note that
c     the matrix in 'c' projects out the constraints so is the 
c     complement of what we need here).
c
      call geom_hnd_get_data('p', x, nzvar*nzvar) ! P -> x
      call geom_hnd_get_data('c', y, nzvar*nzvar) ! complement(C) -> y
      if (odebug) then
         write(LuOut,*) ' Complement(C)', nzvar
         call output(y, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
      endif
      call dscal(nzvar*nzvar, -1d0, y, 1)
      do i = 1, nzvar
         y(i + (i-1)*nzvar) = y(i + (i-1)*nzvar) + 1d0
      enddo
      if (odebug) then
         write(LuOut,*) ' C'
         call output(y, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
      endif
      if(dabssum(nzvar*nzvar, y).lt.1d-12) then
         if(ga_nodeid().eq.0) write(6,*) 
     W        ' no constraints, skipping ',dabssum(nzvar*nzvar,y)
         call dfill(nzvar*nzvar,0d0,z,1)
      else
c
c     Now have P in x, C in y.  Form CPC and then do generalized
c     inverse.  Both P and C are symmetric.
c
      call dgemm('n', 'n', nzvar, nzvar, nzvar, 1d0, 
     $     y, nzvar, x, nzvar, 0d0, z, nzvar) ! C*P -> z
      call dgemm('n', 'n', nzvar, nzvar, nzvar, 1d0, 
     $     z, nzvar, y, nzvar, 0d0, x, nzvar) ! C*P*C -> x
c
      call dsyev('v', 'u', nzvar, x, nzvar, e, z, ndim**2, info)

      if (info .ne. 0) call errquit('geom_binvr: dsyev2 failed', info,
     &       GEOM_ERR)
      if (odebug) then
         write(LuOut,*) ' Eigenvalues of CPC '
         call doutput(e, 1, nzvar, 1, 1, nzvar, 1, 1)
         write(LuOut,*) ' Eigenvectors of CPC'
         call output(x, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
      endif
c
      nindep = 0
      do i = 1, nzvar
         if (e(i) .gt. 1d-12) then
            nindep = nindep + 1
            e(i) = 1d0 / sqrt(e(i))
         else
            if (e(i) .lt. -1d-12) write(LuOut,*)
     $           ' warning: negative eval in CPC ', i, e(i)
            e(i) = 0.0d0
         endif
         call dscal(nzvar, e(i), x(1+(i-1)*nzvar), 1)
      enddo
      call dgemm('n', 't', nzvar, nzvar, nzvar, 1d0, 
     $     x, nzvar, x, nzvar, 0d0, z, nzvar) ! (CPC)^-1 -> z
c
c     Have (CPC)^-1 in z, C in y.  Transform with C then 
c     get P back and do it with that.
c
      call dgemm('n', 'n', nzvar, nzvar, nzvar, 1d0, 
     $     y, nzvar, z, nzvar, 0d0, x, nzvar) ! C*(CPC)^-1 -> x
      call dgemm('n', 'n', nzvar, nzvar, nzvar, 1d0, 
     $     x, nzvar, y, nzvar, 0d0, z, nzvar) ! C*(CPC)^-1C -> z
      endif
c
      call geom_hnd_get_data('p', y, nzvar*nzvar) ! P -> y
      call dgemm('n', 'n', nzvar, nzvar, nzvar, 1d0, 
     $     y, nzvar, z, nzvar, 0d0, x, nzvar) ! P*C*(CPC)^-1*C -> x
      call dgemm('n', 'n', nzvar, nzvar, nzvar, 1d0, 
     $     x, nzvar, y, nzvar, 0d0, z, nzvar) ! P*C*(CPC)^-1*C*P -> z
c
      if (odebug) then
         write(LuOut,*) ' Correction to P from constraints '
         call output(z, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
      endif
      call daxpy(nzvar*nzvar, -1d0, z, 1, y, 1)
      if (odebug) then
         write(LuOut,*) ' P corrected for constraints '
         call output(y, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
      endif
      call geom_hnd_put_data('p', y, nzvar*nzvar)
c
      END
      SUBROUTINE HND_BSTR(EQVAL,NOINT,I,J,C,B,NCART,BOHR)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A
C          BOND STRETCH AS DEFINED BY WILSON (SEE WDC P.55) -----
C
      double precision eqval,bohr
      integer          noint,i,j,ncart
      double precision C(3,1),B(NCART,1)
      integer          nocol1,nocol2,m
      double precision dijsq
      double precision RIJ(3)
      double precision ZERO
      DATA ZERO /0.0D+00/
      DIJSQ = ZERO
      DO 100 M = 1,3
      RIJ(M) = C(M,J)-C(M,I)
  100 DIJSQ = DIJSQ+RIJ(M)*RIJ(M)
      DO 120 M = 1,3
      NOCOL1 = 3*(I-1)+M
      NOCOL2 = 3*(J-1)+M
      B(NOCOL1,NOINT) = -RIJ(M)/ SQRT(DIJSQ)
  120 B(NOCOL2,NOINT) = -B(NOCOL1,NOINT)
      EQVAL = SQRT(DIJSQ)*BOHR
      RETURN
      END
      SUBROUTINE HND_BEND(EQVAL,NOINT,I,J,K,C,B,NCART,RTOD,BOHR)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "errquit.fh"
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS OF A
C          VALENCE ANGLE BENDING COORDINATE AS DEFINED BY WILSON.
C          SEE WDC P. 56 -----
C
C     -----I AND K ARE THE NUMBERS OF THE END ATOMS.  J IS THE
C          NUMBER OF THE CENTRAL ATOM -----
C
      double precision eqval,rtod,bohr
      integer          noint,i,j,k,ncart
      integer          m,nocol1,nocol2,nocol3
      double precision C(3,1),B(NCART,1)
      double precision RJI(3),RJK(3),EJI(3),EJK(3)
      double precision zero,one,tol,pideg
      double precision djisq,djksq
      double precision dotj,dji,djk,dot
      double precision sinj
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DO 100 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      RJK(M) = C(M,K)-C(M,J)
      DJISQ = DJISQ+RJI(M)*RJI(M)
  100 DJKSQ = DJKSQ+RJK(M)*RJK(M)
      DJI =  SQRT(DJISQ)
      DJK =  SQRT(DJKSQ)
      DOTJ = ZERO
      DO 120 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
  120 DOTJ = DOTJ+EJI(M)*EJK(M)
      IF (ONE- ABS(DOTJ)) 200,140,140
  140 continue
      SINJ = SQRT(ONE-DOTJ*DOTJ)
      DJI=DJI*BOHR
      DJK=DJK*BOHR
      DO 160 M = 1,3
      NOCOL1 = 3*(I-1)+M
      if(abs(dotj).eq.1) then
      B(NOCOL1,NOINT) = DOTJ/DJI
      else
      B(NOCOL1,NOINT) = ((DOTJ*EJI(M)-EJK(M)))/(DJI*SINJ)
      endif
      NOCOL2 = 3*(K-1)+M
      if(abs(dotj).eq.1) then
      B(NOCOL2,NOINT) = DOTJ/DJK
      else
      B(NOCOL2,NOINT) = ((DOTJ*EJK(M)-EJI(M)))/(DJK*SINJ)
      endif
      NOCOL3 = 3*(J-1)+M
      B(NOCOL3,NOINT) = -B(NOCOL1,NOINT)-B(NOCOL2,NOINT)
 160  continue
      DOT = ZERO
      DO 180 M = 1,3
  180 DOT = DOT+EJI(M)*EJK(M)
      EQVAL = ACOS(DOT)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
      GO TO 220
  200 CONTINUE
      CALL ERRQUIT('HND_BEND: ROUNDOFF ERROR',0,0)
  220 CONTINUE
      RETURN
      END
      SUBROUTINE HND_TORS(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,BOHR)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
C
C     ----- THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR THE
C           TORSION AS DEFINED BY WILSON.  SEE WDC P60.
C           NOTE THAT THE SIGN OF THE SECOND TERM ON THE RIGHT HAND
C           SIDE IN EQUATION (22) IS INCORRECT AS PRINTED ON P60 OF
C           WDC.
C
C     ----- THE ORDER OF ATOMS IS I--J--K--L. IF THE OBSERVER LOOKS
C           ALONG THE VECTOR J-->K WITH J NEARER THE OBSERVER, THEN
C           THE CLOCKWISE ROTATION OF J-->I WHICH SUPERPOSES J--> I
C           WITH K-->L IS A NEGATIVE TORSION ANGLE.
C
      double precision eqval,rtod,bohr
      integer          noint,i,j,k,l,ncart
      CHARACTER*8 ERRMSG
      double precision C(3,1),B(NCART,1)
      double precision RIJ(3),RJK(3),RKL(3)
      double precision EIJ(3),EJK(3),EKL(3)
      double precision CR1(3),CR2(3),D(3)
      double precision DKLSQ,DKL,DIJ,DIJSQ,DJKSQ,DJK,DOTPK,DOT
      double precision DOTPJ,F2,F1,SML
      double precision SINPJ,SMI,SINPK,SMJ
      integer          iw,nocol1,nocol2,nocol3,nocol4,m
      DIMENSION ERRMSG(3)
      double precision ZERO,ONE,TOLRD,TOL,PIDEG
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-TORS  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOLRD     /1.0001D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
C     ----- CALCULATING UNIT VECTORS ALONG I-->J, J-->K AND K-->L.
C
      IW = 6
      DIJSQ = ZERO
      DJKSQ = ZERO
      DKLSQ = ZERO
      DO 120 M = 1,3
         RIJ(M) = C(M,J)-C(M,I)
         DIJSQ  = DIJSQ+RIJ(M)*RIJ(M)
         RJK(M) = C(M,K)-C(M,J)
         DJKSQ  = DJKSQ+RJK(M)*RJK(M)
         RKL(M) = C(M,L)-C(M,K)
         DKLSQ  = DKLSQ+RKL(M)*RKL(M)
  120 CONTINUE
      DIJ =  SQRT(DIJSQ)
      DJK =  SQRT(DJKSQ)
      DKL =  SQRT(DKLSQ)
      DO 180 M = 1,3
         EIJ(M) = RIJ(M)/DIJ
         EJK(M) = RJK(M)/DJK
         EKL(M) = RKL(M)/DKL
  180 CONTINUE
C
C                   -->   --->   --->     -->   --->   --->
C     ----- FORMING CR1 = I--J X J--K AND CR2 = J--K X K--L
C
      CR1(1) = EIJ(2)*EJK(3)-EIJ(3)*EJK(2)
      CR1(2) = EIJ(3)*EJK(1)-EIJ(1)*EJK(3)
      CR1(3) = EIJ(1)*EJK(2)-EIJ(2)*EJK(1)
      CR2(1) = EJK(2)*EKL(3)-EJK(3)*EKL(2)
      CR2(2) = EJK(3)*EKL(1)-EJK(1)*EKL(3)
      CR2(3) = EJK(1)*EKL(2)-EJK(2)*EKL(1)
C
C     ----- CALCULATING COS AND SIN OF ANGLES I-J-K AND J-K-L.
C
      DOTPJ = ZERO
      DOTPK = ZERO
      DO 220 M = 1,3
         DOTPJ = DOTPJ-EIJ(M)*EJK(M)
         DOTPK = DOTPK-EJK(M)*EKL(M)
  220 CONTINUE
      IF (ONE- ABS(DOTPJ)) 320,340,240
  240 continue
      if((one-abs(dotpk)).lt.1d-10) goto 380
      IF (ONE- ABS(DOTPK)) 360,380,260
  260 continue
      SINPJ =  SQRT(ONE-DOTPJ*DOTPJ)
      SINPK =  SQRT(ONE-DOTPK*DOTPK)
C
C     ----- CALCULATING THE B MATRIX ELEMENTS.
C
      DIJ=DIJ*BOHR
      DJK=DJK*BOHR
      DKL=DKL*BOHR
      DO 280 M = 1,3
         SMI = -CR1(M)/(DIJ*SINPJ*SINPJ)
         NOCOL1 = 3*(I-1)+M
         B(NOCOL1,NOINT) = SMI
         F1 = (CR1(M)*(DJK-DIJ*DOTPJ))/(DJK*DIJ*SINPJ*SINPJ)
         F2 = (DOTPK*CR2(M))/(DJK*SINPK*SINPK)
         SMJ = F1-F2
         NOCOL2 = 3*(J-1)+M
         B(NOCOL2,NOINT) = SMJ
         SML = CR2(M)/(DKL*SINPK*SINPK)
         NOCOL3 = 3*(L-1)+M
         B(NOCOL3,NOINT) = SML
         NOCOL4 = 3*(K-1)+M
         B(NOCOL4,NOINT) = -(SMI+SMJ+SML)
  280 CONTINUE
C
C     ----- CALCULATING THE MAGNITUDE OF THE TORSION ANGLE, WHICH IS
C           BETWEEN 0.0 AND 180.0 DEGREES.
C
      DOT = ZERO
      DO 300 M = 1,3
         DOT = DOT+(CR1(M)*CR2(M))/(SINPJ*SINPK)
  300 CONTINUE
      IF( ABS(DOT) .GT. TOLRD) GO TO 400
      IF( ABS(DOT) .GT. ONE  ) DOT = SIGN(ONE,DOT)
      EQVAL = ACOS(DOT)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
C
C     ----- CALCULATING THE TORSION ANGLE WITH THE CONVENTIONAL SIGN.
C
C                       -PI .LE. TORSION ANGLE .LE. +PI
C
C               -->  -->   -->
C           HERE D = CR1 X CR2
C
      D(1) = CR1(2)*CR2(3)-CR1(3)*CR2(2)
      D(2) = CR1(3)*CR2(1)-CR1(1)*CR2(3)
      D(3) = CR1(1)*CR2(2)-CR1(2)*CR2(1)
      EQVAL =   SIGN(EQVAL,(D(1)*EJK(1)+D(2)*EJK(2)+D(3)*EJK(3)))
      RETURN
C
C     ----- ERROR SITUATIONS.
C
  320 WRITE (IW,9999) I,J,K,DOTPJ
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  340 continue
      WRITE (IW,9998) I,J,K
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  360 WRITE (IW,9999) J,K,L,DOTPK
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  380 eqval=0d0               
c     torsion angle equal to zero when atoms are collinear
cedo  WRITE (IW,9998) J,K,L
cedo  CALL HND_HNDERR(3,ERRMSG)
      RETURN
  400 WRITE (IW,9997) I,J,K,L
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(' COS ',I3,' - ',I3,' - ',I3,' = ',1PE12.5,'?',
     1       ' ROUNDOFF ERROR. RUN STOPPED.')
 9998 FORMAT(' ATOMS ',I3,1X,I3,1X,I3,' ARE COLLINEAR. RUN STOPPED.')
 9997 FORMAT(' ROUNDOFF ERROR IN CALCULATING TORSION ANGLE FOR',
     1       I3,' - ',I3,' - ',I3,' - ',I3,'. RUN STOPPED.')
      END
      SUBROUTINE HND_OPLA(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,PIO2,BOHR)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "errquit.fh"
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR AN
C          OUT OF PLANE WAGGING COORDINATE AS DEFINED BY WILSON.
C          SEE WDC P58. -----
C     -----I IS THE END ATOM, J IS THE APEX ATOM, AND
C          K AND L ARE THE ANCHOR ATOMS -----
C
      double precision eqval,rtod,pio2,bohr
      integer          noint,i,j,k,l,ncart
      double precision B(NCART,1),C(3,1)
      double precision RJI(3),RJK(3),RJL(3)
      double precision EJI(3),EJK(3),EJL(3)
      double precision C1(3),C2(3),C3(3)
      double precision zero,one,tol,pideg
      double precision DJISQ,DJKSQ,DJLSQ,DJI,DJK,DJL,DET,COST,DOTI,DOT
      double precision SMK,SINI,SMI,SINT,TANT,SML,THETA
      integer          m,nocol1,nocol2,nocol3,nocol4
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DJLSQ = ZERO
      DO 100 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      DJISQ  = DJISQ+RJI(M)*RJI(M)
      RJK(M) = C(M,K)-C(M,J)
      DJKSQ  = DJKSQ+RJK(M)*RJK(M)
      RJL(M) = C(M,L)-C(M,J)
  100 DJLSQ  = DJLSQ+RJL(M)*RJL(M)
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DJL = SQRT(DJLSQ)
      DO 120 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
  120 EJL(M) = RJL(M)/DJL
      C1(1) = EJK(2)*EJL(3)-EJK(3)*EJL(2)
      C1(2) = EJK(3)*EJL(1)-EJK(1)*EJL(3)
      C1(3) = EJK(1)*EJL(2)-EJK(2)*EJL(1)
      C2(1) = EJL(2)*EJI(3)-EJL(3)*EJI(2)
      C2(2) = EJL(3)*EJI(1)-EJL(1)*EJI(3)
      C2(3) = EJL(1)*EJI(2)-EJL(2)*EJI(1)
      C3(1) = EJI(2)*EJK(3)-EJI(3)*EJK(2)
      C3(2) = EJI(3)*EJK(1)-EJI(1)*EJK(3)
      C3(3) = EJI(1)*EJK(2)-EJI(2)*EJK(1)
      DOTI = ZERO
      DO 140 M = 1,3
  140 DOTI = DOTI+EJK(M)*EJL(M)
      IF (ONE- ABS(DOTI)) 240,240,160
  160 SINI = SQRT(ONE-DOTI*DOTI)
      DET = EJI(1)*C1(1)+EJI(2)*C1(2)+EJI(3)*C1(3)
      SINT = DET/SINI
      IF (ONE- ABS(SINT)) 260,260,180
  180 COST = SQRT(ONE-SINT*SINT)
      TANT = SINT/COST
      DJI=DJI*BOHR
      DJK=DJK*BOHR
      DJL=DJL*BOHR
      DO 200 M = 1,3
      NOCOL1 = 3*(I-1)+M
      SMI = ((C1(M)/(COST*SINI))-(TANT*EJI(M)))/DJI
      B(NOCOL1,NOINT) = SMI
      NOCOL2 = 3*(K-1)+M
      SMK = ((C2(M)/(COST*SINI))-((TANT*(EJK(M)-DOTI*EJL(M)))/(SINI*
     1     SINI)))/DJK
      B(NOCOL2,NOINT) = SMK
      NOCOL3 = 3*(L-1)+M
      SML = ((C3(M)/(COST*SINI))-((TANT*(EJL(M)-DOTI*EJK(M)))/(SINI*
     1     SINI)))/DJL
      B(NOCOL3,NOINT) = SML
      NOCOL4 = 3*(J-1)+M
      B(NOCOL4,NOINT) = -(SMI+SMK+SML)
  200 CONTINUE
      DOT = ZERO
      DO 220 M = 1,3
  220 DOT = DOT+C1(M)*EJI(M)/SINI
      THETA = ACOS(DOT)
      THETA = PIO2-THETA
      EQVAL = THETA/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
      GO TO 260
  240 CONTINUE
      call errquit('hnd_opla: roundoff error',0,
     &       GEOM_ERR)
  260 CONTINUE
      RETURN
      END
      SUBROUTINE HND_LIBE(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,AA,NA,BOHR)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
C
C     ----- THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A
C           PAIR OF PERPENDICULAR LINEAR BENDING COORDINATES.  SEE
C           FERIGLE AND MEISTER, J CHEM PHYS 19,982(1951) -----
C
C     ----- I AND K ARE THE END ATOMS AND J IS THE CENTRAL ATOM.
C           NO2 IS THE NUMBER OF THE SECOND INTERNAL COORDINATE,
C           I.E. THE ONE IN A PLANE PERPENDICULAR TO THE FIRST. A(I)
C           GIVES THE CARTESIAN COORDINATES OF A POINT IN SPACE SUCH
C           THAT THE VECTOR FROM ATOM J TO POINT A DEFINES
C           THE ORIENTATION OF THE COORDINATES IN SPACE.
C           THE FIRST OF THE TWO INTERNAL COORDINATES IS IN THE I-A-K
C           PLANE AND THE SECOND IN A PLANE PERPENDICULAR TO THE FIRST
C           THROUGH POINTS I,J, AND K -----
C
#include "stdio.fh"
      double precision eqval,rtod,bohr
      integer noint,i,j,k,l,ncart,na
      CHARACTER*8 ERRMSG
      double precision B(NCART,1),C(3,1)
      double precision AA(3,1),A(3)
      double precision RJI(3),RJK(3)
      double precision UNIT(3),UP(3),UN(3)
      double precision EJI(3),EJK(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-LIBE  -'/
      double precision zero,one,tol,pideg,pt0001
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
      DATA PT0001    /0.0001D+00/
      integer iw,m,nocol1,nocol2,nocol3,no2
      double precision DJISQ,DJKSQ,DAJSQ,DJI,DJK,DAJ,DOTJ,DOTP,TEST
      double precision DUM
c
      IW = LuOut
C
      IF(L.LE.NA) GO TO 120
      WRITE(IW,9048) L
      CALL HND_HNDERR(3,ERRMSG)
  120 CONTINUE
      A(1)=AA(1,L)
      A(2)=AA(2,L)
      A(3)=AA(3,L)
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DAJSQ = ZERO
      DO 140 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      DJISQ  = DJISQ+RJI(M)*RJI(M)
      RJK(M) = C(M,K)-C(M,J)
      DJKSQ  = DJKSQ+RJK(M)*RJK(M)
      UN(M)  = A(M)-C(M,J)
  140 DAJSQ  = DAJSQ+UN(M)*UN(M)
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DAJ = SQRT(DAJSQ)
      DOTJ = ZERO
      DOTP = ZERO
      DO 160 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
      UNIT(M) = UN(M)/DAJ
      DOTJ = DOTJ+EJI(M)*EJK(M)
  160 DOTP = DOTP+EJI(M)*UNIT(M)
      TEST = ( ABS(DOTJ)-ONE)
      IF (PT0001- ABS(TEST)) 280,180,180
  180 CONTINUE
      TEST = ( ABS(DOTP)-ONE)
      IF (PT0001- ABS(TEST)) 200,300,300
  200 CONTINUE
      UP(1) = EJK(2)*UNIT(3)-EJK(3)*UNIT(2)
      UP(2) = EJK(3)*UNIT(1)-EJK(1)*UNIT(3)
      UP(3) = EJK(1)*UNIT(2)-EJK(2)*UNIT(1)
      DUM=ONE/ SQRT(UP(1)**2+UP(2)**2+UP(3)**2)
      UP(1)=UP(1)*DUM
      UP(2)=UP(2)*DUM
      UP(3)=UP(3)*DUM
      UNIT(1) = UP(2)*EJK(3)-UP(3)*EJK(2)
      UNIT(2) = UP(3)*EJK(1)-UP(1)*EJK(3)
      UNIT(3) = UP(1)*EJK(2)-UP(2)*EJK(1)
      DJI=DJI*BOHR
      DJK=DJK*BOHR
C
C     ----- FIRST COMPONENT -----
C
      DO 220 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NOINT) = -UNIT(M)/DJI
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NOINT) = -UNIT(M)/DJK
      NOCOL3 = 3*(J-1)+M
  220 B(NOCOL3,NOINT) = (ONE/DJI+ONE/DJK)*UNIT(M)
      EQVAL = ACOS(DOTJ)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
C
C     ----- SECOND COMPONENT -----
C
      NO2=NOINT+1
      DO 260 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NO2) = -UP(M)/DJI
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NO2) = -UP(M)/DJK
      NOCOL3 = 3*(J-1)+M
  260 B(NOCOL3,NO2) = (ONE/DJI+ONE/DJK)*UP(M)
      GO TO 320
  280 CONTINUE
      WRITE (IW,9008) I,J,K,DOTJ
      CALL HND_HNDERR(3,ERRMSG)
  300 CONTINUE
      WRITE (IW,9028) I,J,K,L
      CALL HND_HNDERR(3,ERRMSG)
  320 CONTINUE
      RETURN
 9008 FORMAT(' IN -LIBE- I,J,K ARE NOT COLLINEAR. I,J,K,DOTJ = ',
     1 3I5,F15.9)
 9028 FORMAT(' IN -LIBE- DUMMY ATOM IS COLLINEAR WITH I,J,K.',
     1 ' I,J,K,L = ',4I5)
 9048 FORMAT(' IN -LIBE- DUMMY ATOM IS OUT OF RANGE. L = ',I5)
      END
      SUBROUTINE HND_DIHPLA(DIHANG,NOINT,IZMAT,CARTC,BMAT,NCART,
     1                      DTORAD)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "stdio.fh"
      double precision dihang,dtorad
      integer noint,ncart
      CHARACTER*8 ERRMSG
      integer          IZMAT(5)
      double precision CARTC(3,1), BMAT(NCART,1),
     1          A(3), B(3), C(3), D(3), E1(3), E2(3), E3(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG   /'PROGRAM ','STOP IN ','-DIHPLA-'/
      double precision zero,one,tol
      DATA ZERO,ONE /0.0D+00,1.0D+00/
      DATA TOL      /1.0D-06/
      integer          i,j,k,l,m,n,iw,iatom,jatom,katom,latom,matom,ixyz
      double precision b1,b2,b3,b4,b5,ADOTE2,F1,F2,F4,F5,E1DE2,CDOTE1
      double precision E2MAG,E1MAGI,E1MAG,E2MAGI,SINDI
C
C
C     COMPUTE THE B MATRIX AND DIHEDRAL ANGLE BETWEEN 5 ATOMS
C     FORMING TWO PLANES HAVING ONLY ONE ATOM IN COMMON.
C
C     WRITTEN BY STEPHEN T. ELBERT, AMES LABORATORY, ISU, 2/23/83
C
C     ON ENTRY,
C     IW     - PRINT MESSAGE FILE
C     NOINT  - NUMBER OF THE INTERNAL COORDINATE BEING EVALUATED
C     NCART  - NUMBER OF COORDINATES (3 * NATOMS)
C     IZMAT  - LIST OF THE FIVE ATOMS INVOLVED
C     CARTC  - CARTESIAN COORDINATES OF THE ATOMS
C     DTORAD - CONVERSION FACTOR (RADIANS TO DEGREES)
C
C     ON EXIT,
C     DIHANG - DIHEDRAL ANGLE, IN DEGREES
C     BMAT   - THE NOINTth COLUMN IS COMPUTED BY THIS ROUTINE
C
C     THE FIRST THREE ATOMS DEFINE THE FIRST PLANE,
C     THE LAST THREE ATOMS DEFINE THE SECOND PLANE,
C     THE MIDDLE ATOM IS IN BOTH PLANES.
C
      IW = LuOut
C
      IATOM = IZMAT(1)
      JATOM = IZMAT(2)
      KATOM = IZMAT(3)
      LATOM = IZMAT(4)
      MATOM = IZMAT(5)
C
C     DEFINE VECTORS FROM THE PIVOT ATOM TO THE OTHER ATOMS
C     A= I-K, B=J-K, C=L-K, D=M-K
C
      DO 110 IXYZ = 1,3
         A(IXYZ) = CARTC(IXYZ,IATOM) - CARTC(IXYZ,KATOM)
         B(IXYZ) = CARTC(IXYZ,JATOM) - CARTC(IXYZ,KATOM)
         C(IXYZ) = CARTC(IXYZ,LATOM) - CARTC(IXYZ,KATOM)
         D(IXYZ) = CARTC(IXYZ,MATOM) - CARTC(IXYZ,KATOM)
  110 CONTINUE
C
C     E1 = A X B   IS PERPENDICULAR TO FIRST PLANE
C
      E1(1) = A(2)*B(3) - A(3)*B(2)
      E1(2) = A(3)*B(1) - A(1)*B(3)
      E1(3) = A(1)*B(2) - A(2)*B(1)
      E1MAG =  SQRT( E1(1)**2 + E1(2)**2 + E1(3)**2)
      IF(E1MAG.LT.TOL) GO TO 810
C
C     E2 = C X D   IS PERPENDICULAR TO SECOND PLANE
C
      E2(1) = C(2)*D(3) - C(3)*D(2)
      E2(2) = C(3)*D(1) - C(1)*D(3)
      E2(3) = C(1)*D(2) - C(2)*D(1)
      E2MAG =  SQRT( E2(1)**2 + E2(2)**2 + E2(3)**2)
      IF(E2MAG.LT.TOL) GO TO 820
C
C     NORMALIZE E1 AND E2
C
      E1MAGI = ONE/E1MAG
      E2MAGI = ONE/E2MAG
      DO 120 IXYZ=1,3
         E1(IXYZ) = E1(IXYZ) * E1MAGI
         E2(IXYZ) = E2(IXYZ) * E2MAGI
  120 CONTINUE
C
C     THE DIRECTION OF E2 IS DEFINED SUCH THAT WHEN THE PLANES ARE
C     ROTATED ALONG THE INTERSECTION E3 = E1 X E2 SO THEY COINCIDE
C     (DIHEDRAL ANGLE=0), THE ATOMS 1 AND 4 LIE ON THE SAME SIDE OF E3.
C     (IF ATOMS 1 AND/OR 4 LIE ON E3, THEN ATOMS 2 AND/OR 5 ARE USED)
C
      ADOTE2 = A(1)*E2(1) + A(2)*E2(2) + A(3)*E2(3)
      IF(ADOTE2.LT.TOL) ADOTE2 = B(1)*E2(1) + B(2)*E2(2) + B(3)*E2(3)
      CDOTE1 = C(1)*E1(1) + C(2)*E1(2) + C(3)*E1(3)
      IF(CDOTE1.LT.TOL) CDOTE1 = D(1)*E1(1) + D(2)*E1(2) + D(3)*E1(3)
C
C     POSSIBLY CHANGE SENSE OF E2 SO A (OR B) WILL BE ON THE SAME
C     SIDE OF E1 X E3 AS C (OR D)
C
      IF(ADOTE2*CDOTE1 .GT.ZERO) GO TO 140
         ADOTE2 = -ADOTE2
         DO 130 IXYZ=1,3
            E2(IXYZ) = -E2(IXYZ)
  130    CONTINUE
  140 CONTINUE
C
C      EVALUATE DIHEDRAL ANGLE    0 < ABS(DIHANG) < 180
C      ANGLE IS POSITIVE FOR SECOND PLANE 'ABOVE' ATOM I
C      ANGLE IS NEGATIVE FOR SECOND PLANE 'BELOW' ATOM I
C
      E1DE2 = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
      IF((      ABS(E1DE2)-ONE)  .GT.ZERO) GO TO 840
      IF(( ABS( ABS(E1DE2)-ONE)) .LT.TOL ) GO TO 830
      DIHANG = ACOS(E1DE2)
      IF(ADOTE2.GT.ZERO) DIHANG = -DIHANG
      SINDI = ONE / SIN(DIHANG)
      DIHANG = DIHANG/DTORAD
C
C     E3 = E1 X E2 IS ALONG INTERSECTION LINE OF THE PLANES
C
      E3(1) = SINDI * (E1(2)*E2(3) - E1(3)*E2(2))
      E3(2) = SINDI * (E1(3)*E2(1) - E1(1)*E2(3))
      E3(3) = SINDI * (E1(1)*E2(2) - E1(2)*E2(1))
C
      F1 =  E1MAGI * (B(1)*E3(1) + B(2)*E3(2) + B(3)*E3(3))
      F2 = -E1MAGI * (A(1)*E3(1) + A(2)*E3(2) + A(3)*E3(3))
      F4 = -E2MAGI * (D(1)*E3(1) + D(2)*E3(2) + D(3)*E3(3))
      F5 =  E2MAGI * (C(1)*E3(1) + C(2)*E3(2) + C(3)*E3(3))
C
      I = 3*IATOM - 3
      J = 3*JATOM - 3
      K = 3*KATOM - 3
      L = 3*LATOM - 3
      M = 3*MATOM - 3
C
      DO 160 N=1,3
         B1 = E1(N) * F1
         B2 = E1(N) * F2
         B4 = E2(N) * F4
         B5 = E2(N) * F5
         B3 = -(B1+B2+B4+B5)
         BMAT(NOINT,I+N) = B1
         BMAT(NOINT,J+N) = B2
         BMAT(NOINT,K+N) = B3
         BMAT(NOINT,L+N) = B4
         BMAT(NOINT,M+N) = B5
  160 CONTINUE
      RETURN
C                                        ERROR MESSAGES
  810 WRITE(IW,910)
      GO TO 880
  820 WRITE(IW,920)
      GO TO 880
  830 WRITE(IW,930)
      GO TO 880
  840 WRITE(IW,940) E1DE2
  880 CONTINUE
      WRITE(IW,950) NOINT,IATOM,JATOM,KATOM,LATOM,MATOM
      CALL HND_HNDERR(3,ERRMSG)
C
  910 FORMAT(1X,'---- ERROR, FIRST THREE ATOMS ARE COLLINEAR')
  920 FORMAT(1X,'---- ERROR, LAST THREE ATOMS ARE COLLINEAR ')
  930 FORMAT(1X,'---- ERROR, ALL FIVE ATOMS ARE COPLANAR ',
     1          '     THE -B- MATRIX COLUMN CAN NOT BE DEFINED. STOP')
  940 FORMAT(1X,'---- COSIN OF DIHEDRAL ANGLE=',F16.10,
     1          ' IS GREATER THAN ONE')
  950 FORMAT(1X,'THIS MESSAGE COMES TO YOU FROM -DIHPLA-, WHERE ',
     1          'INTERNAL COORDINATE',I4,' USES ATOMS',5I4)
      END
      logical function geom_zmtmak(rtdb, geom, oprint)
c     implicit double precision (a-h,o-z)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "util.fh"
#include "nwc_const.fh"
#include "rtdb.fh"
#include "util_params.fh"
      integer  rtdb
      integer  geom
      logical  oprint
      integer  max_zcoord
      parameter (max_zcoord=160)
c     
c     Make a Z-matrix for the current geometry.  If successful store
c     it in the geometry and return true.  Otherwise return false.
c     
c     This version does not use the old HONDO common blocks so
c     it may be called without preamble or fear of side-effects.
c     
      double precision toangs
      integer          mxatom,mxcart,mxzmat,mxcoor,mxizmt,mxbond,mxbnds
      integer          mxangs,mxtors,mxoopa,mxlinb,mxlnba,mxseg
      parameter (toangs=cau2ang)
      parameter (mxatom=nw_max_atom)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=nw_max_zmat)    
      parameter (mxcoor=nw_max_coor)    
      parameter (mxizmt=nw_max_izmat)
      parameter (mxbond= 384)   ! for specifying additional bonds
      parameter (mxbnds= 8*mxatom)
      parameter (mxangs= 8*mxatom)
      parameter (mxtors= 8*mxatom)
      parameter (mxoopa= 8*mxatom)
      parameter (mxlinb= 8*mxatom)
      parameter (mxlnba=10       )
      parameter (mxseg=64)
      logical okay
      logical remove
      logical nubond
      logical modbnd
      logical modang
      logical modtor
      logical endatm
      logical endmod
      logical  done
      logical idone
      logical zdone
      logical geom_autoz_info_get
      external geom_autoz_info_get
      integer iw

*     common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
*     common/hnd_zmtpar/nzmat,nzvar,nvar
*     common/hnd_zmtdat/zmat(mxcoor),izmat(mxizmt),nizmat

      character*16 tags(mxatom)
      double precision cvr(mxatom) ! Gathered radii
      double precision zan(mxatom)
      double precision c(3,mxatom)
      integer nat
      integer nzmat, nzvar, izmat(mxizmt), nizmat
      double precision xyzlnb(3,mxlnba)
*
      integer ijbond(2,mxbond)
      integer ijkang(3,mxangs)
      integer ijklto(4,mxtors)
      integer ijklop(4,mxoopa)
      integer ijklnb(4,mxlinb)
*
      double precision ijbond_val(mxbond)
      double precision ijkang_val(mxangs)
      double precision ijklto_val(mxtors)
      double precision ijklop_val(mxoopa)
      double precision ijklnb_val(mxlinb)
*
      logical ijbond_frz(mxbond)
      logical ijkang_frz(mxangs)
      logical ijklto_frz(mxtors)
      logical ijklop_frz(mxoopa)
      logical ijklnb_frz(mxlinb)
*
      character*8 ijbond_nam(mxbond)
      character*8 ijkang_nam(mxangs)
      character*8 ijklto_nam(mxtors)
      character*8 ijklop_nam(mxoopa)
      character*8 ijklnb_nam(mxlinb)
*
      integer ijbnds(2,mxbnds)
      integer ijmods(2,mxbnds), ijxtra(2,mxbond)
      integer nibond(mxatom),iibond(mxatom)
      integer nimods(mxatom),iimods(mxatom)
      integer izfrz(mxcoor), nzfrz
      double precision izfrz_val(mxcoor),zmat(mxcoor)
      external geom_compute_zmatrix
      logical  geom_compute_zmatrix, frz_bnd, frz_ang, frz_tor
*
      dimension endatm(mxatom)
      dimension endmod(mxatom)
      dimension  idone(mxatom)
      integer   iatseg(mxatom,mxseg)
      integer   namseg(       mxseg)
      integer   lenseg(       mxseg)
c
      character*8 zvarname(mxcoor)
      double precision zvarsign(mxcoor)
c
      double precision zero
      parameter (zero = 0.0d0)
c
      integer iat,jat,kat,nseg,iadd,i,j,k,l,lcon,ianchr,ii,iatom,jatom
      integer ipass,itype,izvar,icon,jcon,kcon
      integer m,mlnba,mbond,mangs,max_tor_per_bond,mbnds,mlinb,ncvrpass
      integer moopa,mseg,mtors,mxconn
      integer nbnds,npass,nlnba,numang,nuconn,numods,numbnd,numlnb
      integer numijxtra,numtor,numoop,nuser
      double precision dist,cvr_scaling,cvfac,rij,rcv,radius
c     
      dist(iat,jat)=dsqrt((c(1,iat)-c(1,jat))**2+
     1     (c(2,iat)-c(2,jat))**2+
     2     (c(3,iat)-c(3,jat))**2) * toangs
c     
      iw=luout
      zdone = .false.
      if(oprint) write(iw,8888)
c     
c     Get the data originally passed via the common blocks
c     
      if (.not. geom_cart_get(geom, nat, tags, c, zan))
     $     call errquit('geom_zmtmak: failed getting coords',0,
     &       GEOM_ERR)
      do iat = 1, nat
         if (.not. geom_tag_to_covalent_radius(tags(iat),radius))
     $        radius = 2d0/toangs ! Radius = 2 angstrom if no other info
         cvr(iat) = radius * toangs
         cvr(iat) = max(0.4d0,cvr(iat)) ! Default too small for H
*         write(LuOut,*) ' Radius ', iat, cvr(iat)
      enddo
c     
c     these bonds are used to join disconnected
c     clusters in addition to the user defined bonds
c     
      numijxtra = 0
      call ifill(2*mxbond, 0, ijxtra, 1)
c     
      ncvrpass = 0
10101 ncvrpass = ncvrpass + 1
c     
      nzvar =0
      nizmat=0
      do i=1,mxizmt
         izmat(i)=0
      enddo
      do mbnds=1,mxbnds
         ijbnds(1,mbnds)=0
         ijbnds(2,mbnds)=0
      enddo
      do mbond=1,mxbond
         ijbond(1,mbond)=0
         ijbond(2,mbond)=0
      enddo
      do mangs=1,mxangs
         do i=1,3
            ijkang(i,mangs)=0
         enddo
      enddo
      do mtors=1,mxtors
         do i=1,4
            ijklto(i,mtors)=0
         enddo
      enddo
      do moopa=1,mxoopa
         do i=1,4
            ijklop(i,moopa)=0
         enddo
      enddo
      do mlinb=1,mxlinb
         do i=1,4
            ijklnb(i,mlinb)=0
         enddo
      enddo
      do mlnba=1,mxlnba
         do i=1,3
            xyzlnb(i,mlnba)=zero
         enddo
      enddo
      cvr_scaling=1.2d+00
      max_tor_per_bond = 100
c     
c     ----- read imposed bonds, angles, torsions -----    
c     ----- o-o-plane, linear bends, ... if any  -----
c     
c     read in zcoord data from -rtdb-
c     cvr_scaling,ijbond,ijkang,ijklto,ijklop,ijklnb
c     
      if (.not. geom_autoz_info_get(geom,
     $     cvr_scaling, max_tor_per_bond,
     $     ijbond,ijkang,ijklto,ijklop,ijklnb,
     $     ijbond_nam,ijkang_nam,ijklto_nam,ijklop_nam,ijklnb_nam,
     $     ijbond_val,ijkang_val,ijklto_val,ijklop_val,ijklnb_val,
     $     ijbond_frz,ijkang_frz,ijklto_frz,ijklop_frz,ijklnb_frz))then
         continue               ! no zcoord data input
      endif
c     
      if (ncvrpass .eq. 1) cvfac=cvr_scaling ! use default on first pass
c     
      numbnd=0
      do mbond=1,mxbond
         icon=ijbond(1,mbond)
         jcon=ijbond(2,mbond)
         if(icon.gt.0.and.jcon.gt.0) then
            numbnd=numbnd+1
         elseif(icon.eq.0.and.jcon.eq.0) then
         else
            call errquit('autoz: atom label in ijbond incorrect',mbond,
     &       GEOM_ERR)
         endif
      enddo
c     
      modbnd=numbnd.ne.0
c     
      numang=0
      do mangs=1,mxangs
         icon=ijkang(1,mangs)
         jcon=ijkang(2,mangs)
         kcon=ijkang(3,mangs)
         if(icon.gt.0.and.jcon.gt.0.and.kcon.gt.0) then
            numang=numang+1
         elseif(icon.eq.0.and.jcon.eq.0.and.kcon.eq.0) then
         else
            call errquit('autoz: atom label in ijkang incorrect',mangs,
     &       GEOM_ERR)
         endif
      enddo
      modang=numang.ne.0
c     
      numtor=0
      do mtors=1,mxtors
         icon=ijklto(1,mtors)
         jcon=ijklto(2,mtors)
         kcon=ijklto(3,mtors)
         lcon=ijklto(4,mtors)
         if(icon.gt.0.and.jcon.gt.0.and.
     1        kcon.gt.0.and.lcon.gt.0     ) then
            numtor=numtor+1
         elseif(icon.eq.0.and.jcon.eq.0.and. 
     1           kcon.eq.0.and.lcon.eq.0     ) then
         else
            call errquit('autoz: atom label in ijklto incorrect',mtors,
     &       GEOM_ERR)
         endif
      enddo
      modtor=numtor.ne.0
c     
c$$$      numoop=0
c$$$      do moopa=1,mxoopa
c$$$         icon=ijklop(1,moopa)
c$$$         jcon=ijklop(2,moopa)
c$$$         kcon=ijklop(3,moopa)
c$$$         lcon=ijklop(4,moopa)
c$$$         if(icon.gt.0.and.jcon.gt.0.and.
c$$$     1        kcon.gt.0.and.lcon.gt.0     ) then
c$$$            numoop=numoop+1
c$$$         elseif(icon.eq.0.and.jcon.eq.0.and.
c$$$     1           kcon.eq.0.and.lcon.eq.0     ) then
c$$$         else
c$$$            call errquit('autoz: atom label in ijklop incorrect',0)
c$$$         endif
c$$$      enddo
c$$$      modoop=numoop.ne.0
c$$$c     
c$$$      numlnb=0
c$$$      do mlinb=1,mxlinb
c$$$         icon=ijklnb(1,mlinb)
c$$$         jcon=ijklnb(2,mlinb)
c$$$         kcon=ijklnb(3,mlinb)
c$$$         lcon=ijklnb(4,mlinb)
c$$$         if(icon.gt.0.and.jcon.gt.0.and.
c$$$     1        kcon.gt.0.and.lcon.gt.0     ) then
c$$$            numlnb=numlnb+1
c$$$         elseif(icon.eq.0.and.jcon.eq.0.and.
c$$$     1           kcon.eq.0.and.lcon.eq.0     ) then
c$$$         else
c$$$            call errquit('autoz: atom label in ijklnb incorrect',0)
c$$$         endif
c$$$      enddo
c$$$      modlnb=numlnb.ne.0
c     
c     ----- set up connectivity table, including forced bonds -----
c     
      mbnds=0
      do iat=1,nat
         iibond(iat)=mbnds
c     
c     ----- first check forced bonds -----
c     
         do mbond=1,numijxtra   ! these forced by connection of clusters
            nubond=.false.
            if(ijxtra(1,mbond).eq.iat) then
               mbnds=mbnds+1
               ijbnds(1,mbnds)=iat
               ijbnds(2,mbnds)=ijxtra(2,mbond)
               nubond=.true.
            elseif(ijxtra(2,mbond).eq.iat) then
               mbnds=mbnds+1
               ijbnds(1,mbnds)=iat
               ijbnds(2,mbnds)=ijxtra(1,mbond)
               nubond=.true.
            endif
c     
c     ----- check for duplicates and eliminate -----
c     
            if(nubond.and.mbnds.gt.1) then
               m=mbnds-1
               do i=1,m
                  if(ijbnds(1,i).eq.ijbnds(1,mbnds).and.
     1                 ijbnds(2,i).eq.ijbnds(2,mbnds)     ) then
                     mbnds=mbnds-1
                  endif
               enddo
            endif
         enddo
c     
         if(modbnd) then        ! these forced by the user
            do mbond=1,mxbond
               nubond=.false.
               if(ijbond(1,mbond).eq.iat) then
                  mbnds=mbnds+1
                  ijbnds(1,mbnds)=iat
                  ijbnds(2,mbnds)=ijbond(2,mbond)
                  nubond=.true.
               elseif(ijbond(2,mbond).eq.iat) then
                  mbnds=mbnds+1
                  ijbnds(1,mbnds)=iat
                  ijbnds(2,mbnds)=ijbond(1,mbond)
                  nubond=.true.
               endif
c     
c     ----- check for duplicates and eliminate -----
c     
               if(nubond.and.mbnds.gt.1) then
                  m=mbnds-1
                  do i=1,m
                     if(ijbnds(1,i).eq.ijbnds(1,mbnds).and.
     1                    ijbnds(2,i).eq.ijbnds(2,mbnds)     ) then
                        mbnds=mbnds-1
                     endif
                  enddo
               endif
            enddo
         endif
c     
c     ----- discover connectivity to other atoms -----
c     
         do jat=1,nat
            if(jat.ne.iat) then
               nubond=.false.
               rij=dist(iat,jat)
               rcv=(cvr(iat)+cvr(jat))
               if(rij.le.rcv*cvfac) then
                  if(mbnds.eq.mxbnds) call errquit
     $                 ('autoz: too many bonds found',mbnds, GEOM_ERR)
                  mbnds=mbnds+1
                  ijbnds(1,mbnds)=iat
                  ijbnds(2,mbnds)=jat
                  nubond=.true.
c     
c     ----- check for duplicates and eliminate -----
c     
                  if(nubond.and.mbnds.gt.1) then
                     remove=.false.
                     m=mbnds-1
                     do i=1,m
                        if(ijbnds(1,i).eq.ijbnds(1,mbnds).and.
     1                       ijbnds(2,i).eq.ijbnds(2,mbnds)     ) then
                           remove=.true.
                        endif
                     enddo
                     if(remove) then
                        ijbnds(1,mbnds)=0
                        ijbnds(2,mbnds)=0
                        mbnds=mbnds-1
                     endif
                  endif
               endif
            endif
         enddo
c     
         nibond(iat)=mbnds-iibond(iat)
         if(nibond(iat).gt.1) then
            endatm(iat)=.false.
            if (nibond(iat).gt.6 .and. oprint) 
     $           write(LuOut,30304) nibond(iat),iat
30304       format(1x,'     warning. autoz generated',i5,
     $           ' bonds for atom',i5)
         else 
            endatm(iat)=.true.
c$$$  else if (nibond(iat).le.0) then ! check for unconnected atoms
c$$$  write(LuOut,30303) iat
c$$$  30303       format(1x,'      autoz did not connect atom ',i5)
c$$$  goto 55555
         endif
      enddo
      nbnds=mbnds
c     
      if(nbnds.le.0) then
**         write(iw,9999)
         goto 55555             ! increase the radius
      endif
c     
c     -----    terminal atoms have been detected   -----
c     try for segments ( rings and bridges )
c     
      do i=1,nbnds
         ijmods(1,i)=ijbnds(1,i)
         ijmods(2,i)=ijbnds(2,i)
      enddo
      mxconn=0
      do iat=1,nat
         endmod(iat)=endatm(iat)
         iimods(iat)=iibond(iat)
         nimods(iat)=nibond(iat)
         if(nimods(iat).gt.mxconn) then
            mxconn=nimods(iat)
         endif
         do mseg=1,mxseg
            iatseg(iat,mseg)=0
         enddo
      enddo
      do mseg=1,mxseg
         namseg(mseg)=0
         lenseg(mseg)=0
      enddo
c     
c     -----     take away -end atom- connectivity       -----
c     ring and bridge connectivity should be left
c     
      npass=mxconn
      do ipass=1,npass
         do iat=1,nat
            if(nimods(iat).eq.1) then
               do i=iibond(iat)+1,iibond(iat)+nibond(iat)
                  kat=ijmods(1,i)
                  jat=ijmods(2,i)
                  ijmods(1,i)=0
                  ijmods(2,i)=0
                  if(kat.ne.0) then
                     if(nimods(jat).ge.1) then
                        do j=iibond(jat)+1,iibond(jat)+nibond(jat)
                           if(ijmods(2,j).eq.iat) then
                              ijmods(1,j)=0
                              ijmods(2,j)=0
                              nimods(jat)=nimods(jat)-1
                           endif
                        enddo
                     else
                        write(iw,9978) jat,iat
                        call errquit('autoz: code logic error',0,
     &       GEOM_ERR)
                     endif
                     if(nimods(jat).eq.1) then
                        endmod(jat)=.true.
                     endif
                  endif
               enddo
               nimods(iat)=nimods(iat)-1
            endif
         enddo
      enddo
c     
c     ----- check if there will be more to do ... -----
c     
      done=.true.
      do iat=1,nat
         if(nimods(iat).gt.0) then
            idone(iat)=.false.
         else
            idone(iat)=.true. 
         endif
         done=done.and.idone(iat)
      enddo
      done=.true.
      if(done) then
c     
c     ----- create -izmat- -----
c     
         nzvar =0
         nizmat=0
         do i=1,mxizmt
            izmat(i)=0
         enddo
         do i=1,nbnds
            ijmods(1,i)=ijbnds(1,i)
            ijmods(2,i)=ijbnds(2,i)
         enddo
         do iat=1,nat
            endmod(iat)=endatm(iat)
            iimods(iat)=iibond(iat)
            nimods(iat)=nibond(iat)
         enddo
c     
c     check for disconnected clusters and attempt
c     to join them by putting extra bonds into ijxtra
c     
         call geom_check_all_connected(nat,nbnds,ijmods,
     $        c, numijxtra, ijxtra, mxbond, 3d0*cvfac, 
     $        cvr, okay, oprint)
c     
         if (.not. okay) then
*            if (cvfac .lt. 1.19d0) then
*               goto 55555       ! prefer to first try increasing cvfac
*            else if (numijxtra .eq. 0) then
            if (numijxtra .eq. 0) then
               goto 55555       ! reconnection failed
            else
               if (oprint) write(LuOut,*) 
     $              ' autoz: regenerating connections with new bonds'
               goto 10101
            endif
         endif
c     
         call geom_zmtbld(zdone,nzvar,nizmat,izmat,mxizmt,c,
     1        nat,nbnds,ijmods,iimods,nimods,numbnd,
     2        ijkang,numang,ijklto,numtor,
     3        ijklop,numoop,ijklnb,numlnb,xyzlnb,
     $        max_tor_per_bond, oprint)
c     
         nzmat=nizmat
         nlnba=3*mxlnba
         if (.not. zdone) goto 55555 ! attempt to recover
c
cold         if (nzvar .gt. 10*(max(6,3*nat-6))) then
         if (nzvar .gt. 42*nat) then
c
c     Made a z-matrix but it is asburdly big.  Nothing yet to fix this
c
            write(LuOut,*) ' autoz: excessive number of variables ',
     $           nzvar, 42*nat
            zdone = .false.
            goto 66666          ! abandon autoz
         endif
c     
c     This is the succesful return
c
c     Set variable names and signs, mark constants, track values of constants.
c
c     ----- save -izmat- and -zmat- -----
c
         do i = 1, nzvar
            zvarname(i) = ' '
            zvarsign(i) = 1d0
         enddo
c
         if (.not.geom_autoz_info_get(geom,
     $        cvr_scaling, max_tor_per_bond,
     $        ijbond,ijkang,ijklto,ijklop,ijklnb,
     $        ijbond_nam,ijkang_nam,ijklto_nam,ijklop_nam,ijklnb_nam,
     $        ijbond_val,ijkang_val,ijklto_val,ijklop_val,ijklnb_val,
     $        ijbond_frz,ijkang_frz,ijklto_frz,ijklop_frz,ijklnb_frz))
     $        then
            modbnd = .false.
            modang = .false.
            modtor = .false.
            continue            ! no zcoord data input
         endif
c
c     Match generated variables against the user input to propagate
c     constraint information into the final zmatrix
c
         nzfrz = 0
         if (.not.rtdb_get(rtdb, 'zmat:frz_bnd', MT_LOG, 1, frz_bnd))
     $       frz_bnd=.false.
         if (.not.rtdb_get(rtdb, 'zmat:frz_ang', MT_LOG, 1, frz_ang))
     $       frz_ang=.false.
         if (.not.rtdb_get(rtdb, 'zmat:frz_tor', MT_LOG, 1, frz_tor))
     $       frz_tor=.false.
         if (frz_bnd.or.frz_ang.or.frz_tor) then
            if (.not. geom_zmt_put_data(geom, nizmat, izmat,
     $           nzfrz, izfrz, izfrz_val, nzvar, zvarname, zvarsign,
     $           'autoz')) call
     $           errquit('geom_zmtmak: put data failed',0, GEOM_ERR)
            done=geom_compute_zmatrix(geom, zmat)
            iadd = 1
            do izvar = 1, nzvar
               itype = izmat(iadd)
               if (itype .eq. 1) then
                  iadd = iadd + 3
                  if (frz_bnd) then
                     nzfrz = nzfrz + 1
                     izfrz(nzfrz) = izvar
                     izfrz_val(nzfrz) = zmat(izvar)
                  endif
               else if (itype .eq. 2) then
                  iadd = iadd + 4
                  if (frz_ang) then
                     nzfrz = nzfrz + 1
                     izfrz(nzfrz) = izvar
                     izfrz_val(nzfrz) = zmat(izvar)
                  endif
               else if (itype .eq. 3) then
                  iadd = iadd + 5
                  if (frz_tor) then
                     nzfrz = nzfrz + 1
                     izfrz(nzfrz) = izvar
                     izfrz_val(nzfrz) = zmat(izvar)
                  endif
               endif
            enddo
         endif
c
         if (modbnd .or. modang .or. modtor) then
            nuser = 0
            iadd = 1
            do izvar = 1, nzvar
               itype = izmat(iadd)
               i = izmat(iadd+1)
               j = izmat(iadd+2)
               k = izmat(iadd+3)
               l = izmat(iadd+4)
               if (itype .eq. 1) then
                  iadd = iadd + 3
                  do ii = 1, max_zcoord
                     if ( (i.eq.ijbond(1,ii).and.j.eq.ijbond(2,ii)).or.
     $                    (j.eq.ijbond(1,ii).and.i.eq.ijbond(2,ii)) )
     $                    then
                        if (ijbond_nam(ii)(1:1) .eq. '-') 
     $                       zvarsign(izvar) = -1d0
                        if ( (ijbond_nam(ii)(1:1) .eq. '-') .or.
     $                       (ijbond_nam(ii)(1:1) .eq. '+'))
     $                       ijbond_nam(ii) = ijbond_nam(ii)(2:)
                        zvarname(izvar) = ijbond_nam(ii)
                        nuser = nuser + 1
                        if (ijbond_frz(ii)) then
                           nzfrz = nzfrz + 1
                           izfrz(nzfrz) = izvar
                           izfrz_val(nzfrz) = ijbond_val(ii)
                        endif
                     endif
                  enddo
               else if (itype .eq. 2) then
                  iadd = iadd + 4
                  do ii = 1, max_zcoord
                     if (((i.eq.ijkang(1,ii).and.k.eq.ijkang(3,ii)).or.
     $                    (k.eq.ijkang(1,ii).and.i.eq.ijkang(3,ii)))
     $                    .and. (j.eq.ijkang(2,ii))) then
                        if (ijkang_nam(ii)(1:1) .eq. '-') 
     $                       zvarsign(izvar) = -1d0
                        if ( (ijkang_nam(ii)(1:1) .eq. '-') .or.
     $                       (ijkang_nam(ii)(1:1) .eq. '+'))
     $                       ijkang_nam(ii) = ijkang_nam(ii)(2:)
                        zvarname(izvar) = ijkang_nam(ii)
                        nuser = nuser + 1
                        if (ijkang_frz(ii)) then
                           nzfrz = nzfrz + 1
                           izfrz(nzfrz) = izvar
                           izfrz_val(nzfrz) = ijkang_val(ii)
                        endif
                     endif
                  enddo
               else if (itype .eq. 3) then
                  iadd = iadd + 5
                  do ii = 1, max_zcoord
                     if( (i.eq.ijklto(1,ii).and.j.eq.ijklto(2,ii).and.
     $                    k.eq.ijklto(3,ii).and.l.eq.ijklto(4,ii)) .or.
     $                   (i.eq.ijklto(1,ii).and.j.eq.ijklto(2,ii).and.
     $                    k.eq.ijklto(3,ii).and.l.eq.ijklto(4,ii)) )then
                        if (ijklto_nam(ii)(1:1) .eq. '-') 
     $                       zvarsign(izvar) = -1d0
                        if ( (ijklto_nam(ii)(1:1) .eq. '-') .or.
     $                       (ijklto_nam(ii)(1:1) .eq. '+'))
     $                       ijklto_nam(ii) = ijklto_nam(ii)(2:)
                        zvarname(izvar) = ijklto_nam(ii)
                        nuser = nuser + 1
                        if (ijklto_frz(ii)) then
                           nzfrz = nzfrz + 1
                           izfrz(nzfrz) = izvar
                           izfrz_val(nzfrz) = ijklto_val(ii)
                        endif
                     endif
                  enddo
               else
                  call errquit('geom_zmtmak: izmat?',i, GEOM_ERR)
               endif
            enddo
         endif
         if (.not. geom_zmt_put_data(geom, nizmat, izmat,
     $        nzfrz, izfrz, izfrz_val, nzvar, zvarname, zvarsign,
     $        'autoz'))
     $        call errquit('geom_zmtmak: put data failed',0, GEOM_ERR)
c
         geom_zmtmak = .true.
         return
      endif
c     
c     ----- now compress modified connectivity table -----
c     
      do iat=1,nat
         if(nimods(iat).gt.0) then
            numods=0
            do i=iibond(iat)+1,iibond(iat)+nibond(iat)
               if(ijmods(1,i).ne.0.and.ijmods(2,i).ne.0) then
                  numods=numods+1
                  ijmods(1,numods)=ijbnds(1,i)
                  ijmods(2,numods)=ijbnds(2,i)
               endif
            enddo
            if(nimods(iat).ne.numods) then
               write(iw,9985) iat,nimods(iat),numods
               call errquit('autoz: code logic error',0, GEOM_ERR)
            endif
            iibond(iat)=iimods(iat)
            nibond(iat)=iimods(iat)
            endatm(iat)=endmod(iat)
         endif
      enddo
c     
c     ----- now we have only atoms with remaining connectivity -----
c     of 2 or more. they belong to rings and bridges.
c     
      nuconn=0
      do iat=1,nat
         if(nimods(iat).gt.nuconn) then
            nuconn=nimods(iat)
         endif
      enddo
c     
      if(nuconn.lt.2) then
c     
c     ----- we should have only ( > 2 ) connectivities -----
c     
         call errquit('autoz: nuconn < 2', nuconn, GEOM_ERR)
      elseif(nuconn.gt.3) then
c     
c     ----- -nuconn- > 3 .... not treated yet -----
c     
         call errquit('autoz: nuconn > 2', nuconn, GEOM_ERR)
c     
      elseif(nuconn.eq.3) then
c     
c     ----- -nuconn- = 3 .... several rings , fused or not -----
c     bridges as well.
c     first , list out all the connectivity segments
c     
         nseg  =0
         ianchr=0
 100     ianchr=ianchr+1
         if(ianchr.gt.nat) then
            go to 120
         endif
         if(nimods(ianchr).ne.3) then
            go to 100
         endif
         do i=1,3
            nseg=nseg+1
            if(nseg.gt.mxseg) call errquit
     $           ('autoz:too many segments from 3-connected atoms',nseg,
     &       GEOM_ERR)
c     
            jatom=1
            iatseg(jatom,nseg)=ianchr
c     
            iat=ianchr
            jat=ijmods(2,iimods(iat)+i)
 110        jatom=jatom+1
            iatseg(jatom,nseg)=jat
            if(nimods(jat).eq.2) then
               if(ijmods(2,iimods(jat)+1).ne.iat) then
                  kat=ijmods(2,iimods(jat)+1)
               else
                  kat=ijmods(2,iimods(jat)+2)
               endif
               iat=jat
               jat=kat
               go to 110
            elseif(nimods(jat).eq.3) then
               lenseg(nseg)=jatom
               if(jat.eq.ianchr) then
                  namseg(nseg)=1
               else
                  namseg(nseg)=2
               endif
            else
               call errquit('autoz : should have >= 2 connections',jat,
     &       GEOM_ERR)
            endif
         enddo
c     
 120     continue
c     
      elseif(nuconn.eq.2) then
c     
c     ----- -nuconn- = 2 .... disconnected rings -----
c     
         nseg  =0
         ianchr=0
 200     ianchr=ianchr+1
         if(ianchr.gt.nat) then
            go to 220
         endif
         if(nimods(ianchr).ne.2) then
            go to 200
         endif
         nseg =nseg+1
         iatom=1
         iatseg(iatom,nseg)=ianchr
c     
         iat=ianchr
         jat=ijmods(2,iimods(iat)+1)
 210     iatom=iatom+1
         iatseg(iatom,nseg)=jat
         if(nimods(jat).eq.2) then
            if(ijmods(2,iimods(jat)+1).ne.iat) then
               kat=ijmods(2,iimods(jat)+1)
            else
               kat=ijmods(2,iimods(jat)+2)
            endif
            if(kat.ne.ianchr) then
               iat=jat
               jat=kat
               go to 210
            endif
         else
            call errquit('autoz: should have only two connections',jat,
     &       GEOM_ERR)
         endif
         lenseg(nseg)=iatom
c     
 220     continue
c     
      endif
c     
      geom_zmtmak = .false.
      return                    ! failure
c     
c     get here if something went wrong and we are trying to recover
c     by adjusting cvr_scaling.
c     
55555 continue
      if (cvfac .lt. (cvr_scaling + 0.5d0)) then ! used to be 2.5 upper bound
         if (oprint) write(LuOut,9876) ncvrpass, cvfac, cvfac+0.1d0
 9876    format(1x,i5,' autoz failed with cvr_scaling =',f4.1,
     $        ' changing to',f4.1)
         cvfac = cvfac + 0.1d0
         numijxtra = 0          ! redefine extra bonds
         call ifill(2*mxbond, 0, ijxtra, 1)
         goto 10101
      endif
c
c     Give up.  Punt to cartesians.
c
66666 zdone = .false.
      geom_zmtmak = .false.
      if (oprint) write(LuOut,5505) 
 5505 format(/' AUTOZ failed to generate good internal coordinates.',
     $     /' Cartesian coordinates will be used in optimizations.'/)
      return                    ! failure
c     
 9999 format(' autoz: something is wrong. no bonds could be detected.')
 9998 format(' total number of bonds = ',i5)
 9997 format(1x,'atom   num. of bonds   first bond  endatm',/,
     1     1x,41(1h-))
 9996 format(i5,6x,i5,5x,6x,i5,2x,4x,l1,3x)
 9995 format(' connectivity for atom -iat- = ',i4,' num.bonds = ',i5)
 9994 format(6(1h(,i4,1h,,i4,1h)))
 9993 format(' no covalent radius specified for one of these atoms = ',
     1     2i5)
 9992 format(' iat,jat,iz,jz,rij,rcv,mbnds= ',2i4,2i4,2f10.3,3i5,l4)
 9991 format(/,' pass num. = ',i3,' out of ',i3,' pass(es) ')
 9990 format(/,' searching for a single ring ...')
 9988 format(/,' this routine does not handle -mxconn- .gt. 3. stop.')
 9987 format(/,' a ',i2,'-member ring was found. atoms are = ',/,20i4)
 9986 format(/,' atom ',i3,' should have two connections,',
     1     ' and only two. stop.')
 9985 format(/,' error during compression of connectivity table',
     1     ' for atom = ',i4,/,' nimods(iat), numods = ',2i5)
 9984 format(/,' searching for multiple rings and bridges ...')
 9981 format(/,' searching segment no. = ',i4)
 9980 format(/,' all 1-connectivity completed.')
 9979 format(16i4)
 9978 format(/,' atom ',i4,' has no connectivity left,',/,
     1     ' but the code is trying to remove connection',
     2     ' to atom ',i4,' . stop.')
 9976 format(2i5,3f10.6)
 9975 format(' forced bond inserted = ',3i5)
 9974 format(' bond removed = ',i5,' equivalent to = ',i5)
 8888 format(/,10x,6(1h-),/,10x,'auto-z',/,10x,6(1h-))
 8887 format(' some imposed bonds           in -zmtmak- ')
 8886 format('   no imposed bonds           in -zmtmak- ')
 8885 format(' some imposed angles          in -zmtmak- ')
 8884 format('   no imposed angles          in -zmtmak- ')
 8883 format(' some imposed torsions        in -zmtmak- ')
 8882 format('   no imposed torsions        in -zmtmak- ')
 8881 format(' some imposed o-o-plane-bends in -zmtmak- ')
 8880 format('   no imposed o-o-plane-bends in -zmtmak- ')
 8879 format(' some imposed linear-bends    in -zmtmak- ')
 8878 format('   no imposed linear-bends    in -zmtmak- ')
 8877 format(' zdone = ',l4)
      end
      subroutine geom_zmtbld(zdone,nzvar,nizmat,izmat,mxizmt,c,
     1                      nat,nbnds,ijbnds,iibnds,nibnds,numbnd,
     2                      ijkang,numang,ijklto,numtor,
     3                      ijklop,numoop,
     4                      ijklnb,numlnb,xyzlnb,
     $     max_tor_per_bond, oprint)
c     implicit double precision (a-h,o-z)
      implicit none
#include "stdio.fh"
#include "nwc_const.fh"
      integer nzvar,nizmat,mxizmt
      integer nat,nbnds,numbnd,numang,numtor,numoop,numlnb
      integer max_tor_per_bond
      integer i
      integer mxatom,mxbnds,mxangs,mxtors,mxoopa,mxlinb
      parameter (mxatom=nw_max_atom)        
      parameter (mxbnds= 8*mxatom)
      parameter (mxangs= 8*mxatom)
      parameter (mxtors= 8*mxatom)
      parameter (mxoopa= 8*mxatom)
      parameter (mxlinb= 8*mxatom)
      logical     dbug
      logical     out
      logical     zdone, oprint
      logical     dolinb
      integer iw
      double precision c(3,*)
      integer   izmat(*)   
      integer   ijbnds(2,*)
      integer   iibnds(  *)
      integer   nibnds(  *)
      integer   ijkang(3,*),ijklto(4,*),ijklop(4,*)
      integer   ijklnb(4,*)
      double precision xyzlnb(3,*)
c
      iw=luout
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(dbug) then
         write(iw,9999)
      endif
c
      if(nat.le.1) then
         zdone=.true.
         return
      endif
c
c     ----- bonds -----
c
      call hnd_zmtyp1(nzvar,nizmat,izmat,mxizmt,
     1                nat,nbnds,ijbnds,iibnds,nibnds,numbnd)
      if(dbug) then
         write(iw,9998)
         write(iw,9997) (izmat(i),i=1,nizmat)
         write(iw,9996)
         write(iw,9995) nzvar
      endif
c
c     ----- some checks for diatomics and triatomics -----
c
      if(nat.eq.2) then
         if(nzvar.eq.1) then
            zdone=.true. 
            return
         else
            nzvar=0
            zdone=.false.
            return
         endif
      endif
      if(nat.eq.3) then
         if(nzvar.ne.2.and.nzvar.ne.3) then
            nzvar=0
            zdone=.false.
            return
         endif
      endif
c
c     ----- angles -----
c
      dolinb=.false.
      zdone = .true.
      call hnd_zmtyp2(nzvar,nizmat,izmat,mxizmt,c,
     1                nat,nbnds,ijbnds,iibnds,nibnds,ijkang,numang,
     2                dolinb, zdone, oprint)
c
      if (.not. zdone) return
c
      if(dbug) then
         write(iw,9998)
         write(iw,9997) (izmat(i),i=1,nizmat)
         write(iw,9996)
         write(iw,9994) nzvar
      endif
c
c     ----- more checks for triatomics -----
c
      if(nat.eq.3) then
         if(dolinb) then
            nzvar=0
            zdone=.false.
            return
         endif
      endif
c
c     ----- torsions -----
c
      call hnd_zmtyp3(nzvar,nizmat,izmat,mxizmt,c,
     1                nat,nbnds,ijbnds,iibnds,nibnds,ijklto,numtor,
     $     max_tor_per_bond, zdone, oprint)
      if (.not. zdone) return
c
      if(dbug) then
         write(iw,9998)
         write(iw,9997) (izmat(i),i=1,nizmat)
         write(iw,9996)
         write(iw,9993) nzvar
      endif
c
c     rjh.  michel says oop were working, but generated far
c     .     too many degrees of freedom to be useful.
c
c     .     the linear bends are not working but a fix 
c     .     could be to include a point(s) at some distance
c     .     and just use the regular bends.
c
c     .     in the meantime just punt after the torsions
c     .     and then analyze the generated coordinates to
c     .     make sure that we have at least (3n-6) linearly
c     .     independent coordinates.
c
c$$$      dooopa=.false.
c$$$      if(dooopa) then
c$$$c
c$$$c     ----- out-of-plane bends -----
c$$$c
c$$$         call hnd_zmtyp4(nzvar,nizmat,izmat,mxizmt,
c$$$     1                   nat,nbnds,ijbnds,iibnds,nibnds,ijklop,numoop)
c$$$         if(dbug) then
c$$$            write(iw,9998)
c$$$            write(iw,9997) (izmat(i),i=1,nizmat)
c$$$            write(iw,9996)
c$$$            write(iw,9992) nzvar
c$$$         endif
c$$$      endif
c$$$c
c$$$c     ----- linear bends -----
c$$$c
c$$$      if(dolinb) then
c$$$         call hnd_zmtyp5(nzvar,nizmat,izmat,mxizmt,c,
c$$$     1                   nat,nbnds,ijbnds,iibnds,nibnds,
c$$$     2                   ijklnb,numlnb,xyzlnb,ijkang,numang)
c$$$         if(dbug) then
c$$$            write(iw,9998)
c$$$            write(iw,9997) (izmat(i),i=1,nizmat)
c$$$            write(iw,9996)
c$$$            write(iw,9991) nzvar
c$$$         endif
c$$$         zdone=.false.
c$$$         return
c$$$      endif
c
c     ----- done -----
c
      if (nzvar .lt. max(1,3*nat-6)) then
         write(LuOut,*) ' autoz: insufficient variables ',nzvar,3*nat-6
         zdone=.false.
      else
         zdone=.true.
      endif
c
      return
 9999 format(' in -zmtbld- ',/,' ----------- ')
 9998 format(' $zmat  ',/,' izmat= ')
 9997 format(12(i4,1h,))
 9996 format(' $end   ')
 9995 format(' after -zmtyp1- , nzvar = ',i5)
 9994 format(' after -zmtyp2- , nzvar = ',i5)
 9993 format(' after -zmtyp3- , nzvar = ',i5)
 9992 format(' after -zmtyp4- , nzvar = ',i5)
 9991 format(' after -zmtyp5- , nzvar = ',i5)
      end
      subroutine hnd_zmtyp1(nzvar,nizmat,izmat,mxizmt,
     1                      nat,nbnds,ijbnds,iibnds,nibnds,numbnd)
c     implicit double precision (a-h,o-z)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "nwc_const.fh"
      integer nzvar,nizmat,mxizmt,nat,nbnds,numbnd
      integer ione,mxatom,mxbnds
      parameter (ione=1)
      parameter (mxatom=nw_max_atom)        
      parameter (mxbnds= 8*mxatom)
      logical     dbug
      integer iw,ij,iat,i,ii,ji,jj,ibnd,jbnd,ijbnd,nbnd,nbond
      integer izmat(*)   
      integer ijbond(2,mxbnds)
      integer iibond(  mxatom)
      integer nibond(  mxatom)
      integer ijbnds(2,*)
      integer iibnds(  *)
      integer nibnds(  *)
c
      iw=luout
      dbug=.false.
      if(dbug) then
         write(iw,9999)
      endif
c
      nbond=nbnds
      do ijbnd=1,nbond
         ijbond(1,ijbnd)=ijbnds(1,ijbnd)
         ijbond(2,ijbnd)=ijbnds(2,ijbnd)
      enddo
      do iat=1,nat
         iibond(iat)=iibnds(iat)
         nibond(iat)=nibnds(iat)
      enddo
c
      if(nat.le.1) return
c
c     ----- bonds -----
c
      nbnd=nbond
      do ibnd=1,nbnd
         ii=ijbond(1,ibnd)
         ij=ijbond(2,ibnd)
         if(ii.ne.0) then
            if(ibnd.lt.nbnd) then
               do jbnd=ibnd+1,nbnd
                  ji=ijbond(1,jbnd)
                  jj=ijbond(2,jbnd)
                  if(ji.ne.0) then
                     if((ji.eq.ii.and.jj.eq.ij).or.        
     1                  (ji.eq.ij.and.jj.eq.ii)    ) then
                        ijbond(1,jbnd)=0
                        ijbond(2,jbnd)=0
                     endif
                  endif
               enddo
            endif
            if (nizmat+3 .gt. mxizmt) call errquit
     $           ('autoz bonds: too many internals',0,
     &       GEOM_ERR)
             izmat(nizmat+1)=ione
             izmat(nizmat+2)=ii
             izmat(nizmat+3)=ij   
            nizmat=nizmat+3
            nzvar =nzvar +1
            numbnd=numbnd+1
         endif
      enddo
c
      if(dbug) then
         write(iw,9998)
         write(iw,9997) (izmat(i),i=1,nizmat)
         write(iw,9996)
         write(iw,*)    'nzvar = ',nzvar
      endif
c
      return
 9999 format(' in -zmtyp1- ',/,' ----------- ')
 9998 format(' $zmat  ',/,' izmat= ')
 9997 format(12(i4,1h,))
 9996 format(' $end   ')
      end
      subroutine hnd_zmtyp2(nzvar,nizmat,izmat,mxizmt,c,
     1     nat,nbnds,ijbnds,iibnds,nibnds,ijkang,numang,
     2     dolinb, zdone, oprint)
c     implicit double precision (a-h,o-z)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "nwc_const.fh"
      integer nzvar,nizmat,nat,nbnds,numang,mxizmt
      integer itwo,mxatom,mxbnds,mxangs
      parameter (itwo=2)
      parameter (mxatom=nw_max_atom)        
      parameter (mxbnds= 8*mxatom)
      parameter (mxangs= 8*mxatom)
      logical     dbug
      logical     dolinb
      logical     zdone
      logical     linbnd
      logical     oprint
      integer iw
      double precision c(3,*)
      integer          izmat(*)   
      integer          ijkang(3,*)     
      integer          ijbond(2,mxbnds)
      integer          iibond(  mxatom)
      integer          nibond(  mxatom)
      integer          ijbnds(2,*)
      integer          iibnds(  *)
      integer          nibnds(  *)
      double precision tol,two
      data tol    /1.0d-01/     ! 5.7 degrees
*     data tol    /1.0d-08/
      data two    /2.0d+00/
      integer iat,jat,kat,nbond,ii,jk,i,icon,iang,ij
      integer ijbnd,ik,jcon,jang,jj,ji,klbnd,kcon,jkbnd,lmbnd,lat,mat
      integer mang,nang,nfound,nuser_specified
      double precision cosb,rijsq,rij,rik,riksq,rjksq,rilsq,rjk,rklsq
      double precision rkl
      double precision distsq
c     
      distsq(iat,jat)=(c(1,iat)-c(1,jat))**2+
     1     (c(2,iat)-c(2,jat))**2+
     2     (c(3,iat)-c(3,jat))**2
c     
      iw=luout
      dbug=.false.
      if(dbug) then
         write(iw,9999)
      endif
c     
c     ----- angles -----
c     
      nbond=nbnds
      do ijbnd=1,nbond
         ijbond(1,ijbnd)=ijbnds(1,ijbnd)
         ijbond(2,ijbnd)=ijbnds(2,ijbnd)
      enddo
      do iat=1,nat
         iibond(iat)=iibnds(iat)
         nibond(iat)=nibnds(iat)
      enddo
c     
      zdone = .true.
      if(nat.le.2) return
c     
      nang=numang
      do iat=1,nat
         do ijbnd=iibond(iat)+1,iibond(iat)+nibond(iat)
            icon=ijbond(1,ijbnd)
            jcon=ijbond(2,ijbnd)
            if(icon.ne.iat)
     $           call errquit('autoz: invalid connection',1, GEOM_ERR)
            jat=jcon
            do jkbnd=iibond(jat)+1,iibond(jat)+nibond(jat)
               jcon=ijbond(1,jkbnd)
               kcon=ijbond(2,jkbnd)
               if(jcon.ne.jat)
     $              call errquit('autoz: invalid connection',2,
     &       GEOM_ERR)
               kat=kcon
               if(kat.ne.iat) then
                  rijsq=distsq(iat,jat)
                  rij  =sqrt(rijsq)                          
                  rjksq=distsq(jat,kat)
                  rjk  =sqrt(rjksq)                           
                  riksq=distsq(iat,kat)
                  cosb =(rijsq+rjksq-riksq)/(two*rij*rjk)
                  linbnd = sin(acos(min(1d0,abs(cosb)))) .le. tol
                  if (linbnd .and. nibond(jat).eq.2) then
c     
c     ijk is (nearly) linear.  if j is only bonded to i & k
c     then try to generate ikl where l is bonded to k.
c     note we'll generate kji at another time and also kim, where m
c     is bonded to i.
c     
                     do klbnd=iibond(kat)+1,iibond(kat)+nibond(kat)
                        lat = ijbond(2,klbnd)
                        if (lat.ne.iat .and. lat.ne.jat) then
                           riksq=distsq(iat,kat)
                           rik  =sqrt(riksq)                          
                           rklsq=distsq(kat,lat)
                           rkl  =sqrt(rklsq)                           
                           rilsq=distsq(iat,lat)
                           cosb =(riksq+rklsq-rilsq)/(two*rik*rkl)
                           linbnd = sin(acos(min(1d0,abs(cosb)))) 
     $                          .le. tol
                           if (linbnd) then
                              if (nibond(kat) .eq. 2) then
c     
c     darn it.  ikl is also (nearly) linear and there is only one
c     choice for l.  try ilm.  could be even more open-ended but not
c     necessary (?).
c     
                                 nfound = 0
                                 do lmbnd=iibond(lat)+1,
     $                                iibond(lat)+nibond(lat)
                                    mat = ijbond(2,lmbnd)
                                    if (mat.ne.iat .and. mat.ne.jat
     $                                   .and. mat.ne.kat) then
c     
c     logic needed here to check for linear just in case !!!!!!!!!!!!
c     
                                       if(nang.eq.mxangs) goto 808
                                       nfound = nfound + 1
                                       nang=nang+1
                                       ijkang(1,nang)=iat
                                       ijkang(2,nang)=lat
                                       ijkang(3,nang)=mat
                                    endif
                                 enddo
                                 if ( nfound.eq.0 .and. 
     $                                nibond(iat).eq.1 ) then
c     
c     No connection was found and no alternatives for I
c     
                                    if (oprint) write(LuOut,*)
     $                                   ' autoz: angle generation',
     $                                   ' failed due to 4 colinear',
     $                                   ' atoms ', iat, jat, kat, lat
                                    zdone = .false.
                                    return
                                 endif
                              endif
                           else 
                              if(nang.eq.mxangs) goto 808
                              nang=nang+1
                              ijkang(1,nang)=iat
                              ijkang(2,nang)=kat
                              ijkang(3,nang)=lat
                           endif
                        endif
                     enddo
                  else
                     if(nang.eq.mxangs) goto 808
                     nang=nang+1
                     ijkang(1,nang)=iat
                     ijkang(2,nang)=jat
                     ijkang(3,nang)=kat
                  endif
               endif
            enddo
         enddo
      enddo
      if(dbug) then
         if(nang.gt.1) then
            write(iw,9994)
            do mang=1,nang
               write(iw,9993) mang,ijkang(1,mang),
     1              ijkang(2,mang),ijkang(3,mang)
            enddo
         else
            write(iw,9992)
         endif
      endif
c     
c     ----- eliminate duplicates and clean up -----
c     
      if(nang.ge.2) then
         do iang=1,nang-1
            ii=ijkang(1,iang)
            ij=ijkang(2,iang)
            ik=ijkang(3,iang)
            if(ijkang(1,iang).ne.0) then
               do jang=iang+1,nang
                  ji=ijkang(1,jang)
                  jj=ijkang(2,jang)
                  jk=ijkang(3,jang)
                  if(jj.eq.ij) then
                     if((ii.eq.ji.and.ik.eq.jk).or.
     1                    (ii.eq.jk.and.ik.eq.ji)    ) then
                        ijkang(1,jang)=0
                        ijkang(2,jang)=0
                        ijkang(3,jang)=0
                     endif
                  endif
               enddo
            endif
         enddo
      endif
c     
      nuser_specified = numang  ! Don't remove user defined angles
      numang = 0
      iang=0                  
 100  iang=iang+1
      if(iang.gt.nang) go to 110
      if(ijkang(1,iang).eq.0) then
         if(iang.lt.nang) then
            do jang=iang+1,nang
               ijkang(1,jang-1)=ijkang(1,jang)
               ijkang(2,jang-1)=ijkang(2,jang)
               ijkang(3,jang-1)=ijkang(3,jang)
            enddo
            iang=iang-1
         endif
         nang=nang-1
      else
         rijsq=distsq(ijkang(1,iang),ijkang(2,iang))
         rij  =sqrt(rijsq)                          
         rjksq=distsq(ijkang(2,iang),ijkang(3,iang))
         rjk  =sqrt(rjksq)                           
         riksq=distsq(ijkang(1,iang),ijkang(3,iang))
         cosb =(rijsq+rjksq-riksq)/(two*rij*rjk)
         linbnd = (sin(acos(min(1d0,abs(cosb)))) .le. tol)
         if (linbnd) dolinb = .true.
         if((.not.linbnd) .or. (iang.le.nuser_specified)) then
            if (nizmat+4 .gt. mxizmt) goto 808
            izmat(nizmat+1)=itwo
            izmat(nizmat+2)=ijkang(1,iang)
            izmat(nizmat+3)=ijkang(2,iang)
            izmat(nizmat+4)=ijkang(3,iang)
            nizmat=nizmat+4
            nzvar =nzvar +1
            numang=numang+1
         else
*            write(LuOut,303) ijkang(1,iang), ijkang(2,iang), 
*     &           ijkang(3,iang), 
*     $           linbnd
* 303        format(' rejecting bond angle due to linearity ',
*     $           (2x,3i4,2x,l1))
         endif
      endif
      go to 100
 110  continue
c     
      if(dbug) then
         if(nang.gt.1) then
            write(iw,9994)
            do mang=1,nang
               write(iw,9993) mang,ijkang(1,mang),
     1              ijkang(2,mang),ijkang(3,mang)
            enddo
            if(dolinb) then
               write(iw,9991)
            endif
         else
            write(iw,9992)
         endif
      endif
c     
      if(dbug) then
         write(iw,9998)
         write(iw,9997) (izmat(i),i=1,nizmat)
         write(iw,9996)
         write(iw,*)    'nzvar = ',nzvar
      endif
c     
      zdone = .true.
      return
c     
c     Jump here if have too many of something
 808  if (oprint) write(LuOut,*)' autoz: too many angles'
      zdone = .true.
      return
c     
 9999 format(' in -zmtyp2- ',/,' ----------- ')
 9998 format(' $zmat  ',/,' izmat= ')
 9997 format(12(i4,1h,))
 9996 format(' $end   ')
 9995 format(' something wrong with bond ',i2,' of atom ',i4,
     1     ' with -icon-, -jcon- = ',2i5)
 9994 format(' table of angles = ')
 9993 format(4i5)
 9992 format(' no bond angles found ')
 9991 format(' collinear atoms found. do -linear bends- ... ')
      end
      subroutine hnd_zmtyp3(nzvar,nizmat,izmat,mxizmt,c,
     1     numatom,nbnds,ijbnds,iibnds,nibnds, ! note rename nat -> numatom
     2     ijklto,numtor, max_tor_per_bond, zdone, oprint)
c     implicit double precision (a-h,o-z)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "nwc_const.fh"
      integer nzvar,nizmat,numatom,nbnds,numtor, max_tor_per_bond
      integer ithree,mxatom,mxbnds,mxtors,mxizmt
      parameter (ithree=3)
      parameter (mxatom=nw_max_atom)        
      parameter (mxbnds= 8*mxatom)
      parameter (mxtors= 8*mxatom)
      logical     dbug
      logical     linijk,linjkl
      logical odoit
      logical zdone, oprint, oopb, oi3d, oj3d
      logical geom_add_tor, bonds_span3d
      integer iw
      double precision c(3,*)
      integer          izmat(*)   
      integer          ijklto(4,*)      
      integer          ijbond(2,mxbnds)
      integer          iibond(  mxatom)
      integer          nibond(  mxatom)
      integer          ijbnds(2,*)
      integer          iibnds(  *)
      integer          nibnds(  *)
      double precision tol,two
      data tol    /1.0d-01/     ! sin(5.7 degrees) = 0.1
*     data tol    /1.0d-08/
c$$$  data zero   /0.0d+00/
c$$$  data one    /1.0d+00/
      data two    /2.0d+00/
      integer iat,jat,kat
      double precision distsq
      double precision cosbk,cosbj,rjksq,rik,rijsq,rij,rjk,riksq
      double precision rjl,rjlsq,rklsq,rkl
      integer i,ii,ij,ijbnd,ik,il,jtor,jj,itor,ji,jk,jkbnd,jlbnd,jl
      integer lat,klbnd,lmbnd,mat,mnbnd,mtor,nat,nbond,ntornonlin
      integer ntor_for_ij,ntor

c     
      distsq(iat,jat)=(c(1,iat)-c(1,jat))**2+
     1     (c(2,iat)-c(2,jat))**2+
     2     (c(3,iat)-c(3,jat))**2
c     
      iw=luout
      dbug=.false.
      if(dbug) then
         write(iw,9999)
      endif
c     
c     ----- torsions -----
c     
      nbond=nbnds
      do ijbnd=1,nbond
         ijbond(1,ijbnd)=ijbnds(1,ijbnd)
         ijbond(2,ijbnd)=ijbnds(2,ijbnd)
      enddo
      do iat=1,numatom
         iibond(iat)=iibnds(iat)
         nibond(iat)=nibnds(iat)
      enddo
c     
      oopb = .false.            ! Initially don't do oopb
c
55555 ntor=numtor               ! Jump to here to redo everything
      ntornonlin = 0            ! No. torsions not fudged for linearity
      if(numatom.le.3) return
c     
      if (max_tor_per_bond .eq. 0) return
c     
c     currently construct the torsions from the bond info.
c     it might be easier if we instead drove off the angle info
c     since we already know that the angles are good and much
c     of the same logic is reproduced here.
c     
      do iat=1,numatom
         do ijbnd=iibond(iat)+1,iibond(iat)+nibond(iat)
            jat=ijbond(2,ijbnd)
            if(ijbond(1,ijbnd).ne.iat) 
     $           call errquit('zmtyp3: iat not my bond?',ijbnd,
     &       GEOM_ERR)
            ntor_for_ij = 0
            do jkbnd=iibond(jat)+1,iibond(jat)+nibond(jat)
               kat=ijbond(2,jkbnd)
               if(ijbond(1,jkbnd).ne.jat) 
     $              call errquit('zmtyp3: jat not my bond?',jkbnd,
     &       GEOM_ERR)
               if(kat.ne.iat) then
                  rij=sqrt(distsq(iat,jat))
                  rjk=sqrt(distsq(jat,kat))
                  rik=sqrt(distsq(iat,kat))
                  cosbj=(rij**2+rjk**2-rik**2)/(two*rij*rjk)
                  linijk = sin(acos(min(1d0,abs(cosbj)))) .le. tol
c     
                  do klbnd=iibond(kat)+1,iibond(kat)+nibond(kat)
                     lat=ijbond(2,klbnd)
                     if(ijbond(1,klbnd).ne.kat) 
     $                    call errquit('zmtyp3: kat not my bond?',klbnd,
     &       GEOM_ERR)
c     
c     (h2o)5 test case showed some logic was missing.  also, can avoid
c     generating some duplicates by restricting i>l unless are restricting 
c     to a subset of the torsions or have linear angles.
c     
                     odoit = ( 
     $                    iat.ne.jat .and. iat.ne.kat .and. 
     $                    iat.ne.lat .and. jat.ne.kat .and.
     $                    jat.ne.lat .and. kat.ne.lat )
c     
                     if (odoit) then
                        rkl=sqrt(distsq(kat,lat))
                        rjl=sqrt(distsq(jat,lat))
                        cosbk=(rjk**2+rkl**2-rjl**2)/(two*rjk*rkl)
                        linjkl = sin(acos(min(1d0,abs(cosbk)))) 
     $                       .le. tol
c     
                        if (linijk .and. linjkl) then
c     
c     both bonds are linear.  go fishing for ilmn.  also add ijml
c     which is needed for h-c-c-h-X (x not colinear)
c     
*     write(LuOut,*) ' four colinear centers ',
*     $                          iat, jat, kat, lat
                           do lmbnd=iibond(lat)+1,
     $                          iibond(lat)+nibond(lat)
                              mat=ijbond(2,lmbnd)
                              if (mat.ne.iat .and. mat.ne.jat
     $                             .and. mat.ne.kat) then
*     write(LuOut,*) ' IJML ',iat,jat,mat,lat
                                 if (.not. geom_add_tor(
     $                                iat,jat,mat,lat,
     $                                ntor, ijklto, mxtors, 
     $                                ntor_for_ij)) goto 8000
                                 do mnbnd=iibond(mat)+1,
     $                                iibond(mat)+nibond(mat)
                                    nat=ijbond(2,mnbnd)
                                    if (nat.ne.iat .and. nat.ne.jat
     $                                   .and. nat.ne.kat .and.
     $                                   nat.ne.lat) then
*     write(LuOut,*) ' colinear ijkl -',
*     $                                      ' use ilmn', 
*     $                                   iat, jat, kat, lat, mat, nat
                                       if (.not. geom_add_tor(
     $                                      iat,lat,mat,nat,
     $                                      ntor, ijklto, mxtors, 
     $                                      ntor_for_ij)) goto 8000
                                    endif
                                 enddo
                              endif
                           enddo
                        else if (linijk) then
c     
c     ijk is linear.  instead of ijkl use ijlk.  don't do this
c     now since it will be done as lkji.
c     
                           continue
c$$$  write(LuOut,*) ' linear ijk ', iat, jat, kat, lat
c$$$  iii = iat
c$$$  jjj = jat
c$$$  kkk = lat
c$$$  lll = kat
                        else if (linjkl) then
c     
c     jkl is linear.  first try ijlm where m is bonded to l.
c     also add jikl.  fails jeff.nw without this.
c     
                           do lmbnd=iibond(lat)+1,
     $                          iibond(lat)+nibond(lat)
                              mat=ijbond(2,lmbnd)
                              if (mat.ne.iat .and. mat.ne.jat
     $                             .and. mat.ne.kat) then
*     write(LuOut,*) ' linear jkl - use ijlm', 
*     $                                iat, jat, kat, lat, mat
                                 if (.not. geom_add_tor(iat,jat,lat,mat,
     $                                ntor, ijklto, mxtors,ntor_for_ij))
     $                                goto 8000
                              endif
                           enddo
                           if (.not. geom_add_tor(jat,iat,kat,lat,
     $                          ntor, ijklto, mxtors, ntor_for_ij))
     $                          goto 8000
                        else
c     
c     nothing is linear.  wow.
c     
                           ntornonlin = ntornonlin + 1
                           if (.not. geom_add_tor(iat,jat,kat,lat,
     $                          ntor, ijklto, mxtors, ntor_for_ij))
     $                          goto 8000
                        endif
                        if (ntor_for_ij.eq.max_tor_per_bond)
     $                       goto 1212 ! next ij
c     
                     endif      ! if (odoit)
                  enddo         ! do klbnd
               endif            ! if (kat.ne.iat)
            enddo               ! do jkbnd
c     
            oi3d = bonds_span3d(iat,c,nibond,iibond,ijbond)
            oj3d = bonds_span3d(jat,c,nibond,iibond,ijbond)
            if (ntor_for_ij.eq.0 .and. numatom.gt.3 .and. oopb .and.
     $           (.not. (oi3d .or. oj3d))) then 
c     
c     sigh.  we have not found a valid torsion for this bond.
c     this can happen because of linear angles or even for
c     simple things like formaldehyde.  we've looked at all
c     sequences i-j-k-l ... try now 
c     .           k
c     .          /
c     .     i - j
c     .          \
c     .           l
c     specifying the torsion as i-j-k-l.  could do this as an 
c     out-of-plane bend, but this will work also.
c     
c     another catch is that for four such atoms there is only
c     one independent out-of-plane bend ... so comparison
c     for duplicates must be done without regard to order.
c     
c     but only putting in the unique one will break symmetry
c     so bite the big enchilada and do them all ... the user can always
c     restrict the number per bond.  Must check on linear bonds since
c     if restricting the number must be sure each one counts.
c     
c     using i as the central atom will be considered later.
c     
               do jkbnd=iibond(jat)+1,iibond(jat)+nibond(jat)
                  kat=ijbond(2,jkbnd)
                  if(kat.ne.iat) then
                     rij=sqrt(distsq(iat,jat))
                     rjk=sqrt(distsq(jat,kat))
                     rik=sqrt(distsq(iat,kat))
                     cosbj=(rij**2+rjk**2-rik**2)/(two*rij*rjk)
                     linijk = sin(acos(min(1d0,abs(cosbj)))) .le. tol
                     if (.not. linijk) then
                        do jlbnd=iibond(jat)+1,iibond(jat)+nibond(jat)
*                        do jlbnd=jkbnd+1,iibond(jat)+nibond(jat)
                           lat=ijbond(2,jlbnd)
                           if(lat.ne.iat .and. lat.ne.kat) then
                              rkl=sqrt(distsq(kat,lat))
                              rjl=sqrt(distsq(jat,lat))
                              cosbk=(rjk**2+rkl**2-rjl**2)/(two*rjk*rkl)
                              linjkl = sin(acos(min(1d0,abs(cosbk)))) 
     $                             .le. tol
                              if (.not. linjkl) then
*                                 write(LuOut,*) ' OOPB ',
*     $                                iat, jat, kat, lat
                                 if (.not. geom_add_tor(iat, jat, kat,
     $                                lat, ntor, ijklto, mxtors,
     $                                ntor_for_ij)) goto 8000
                                 if (ntor_for_ij .gt. max_tor_per_bond) 
     $                                goto 1212
                              endif
                           endif ! if (4 distinct atoms)
                        enddo   ! do jlbnd
                     endif      ! if (not linear ijk)
                  endif         ! if (kat .ne. iat)
               enddo            ! do jkbnd
            endif               ! if (missing a torsion)
 1212       continue            ! jump here for next ij
         enddo                  ! do ijbnd
      enddo                     ! do iat
c
      if (numatom.gt.3 .and. (.not. oopb) .and. (ntornonlin.eq.0)) then
         if (oprint) write(LuOut,*) ' Looking for out-of-plane bends'
         oopb = .true.
         goto 55555
      endif
c     
      if(dbug) then
         if(ntor.gt.1) then
            write(iw,9994)
            do mtor=1,ntor
               write(iw,9993) mtor,ijklto(1,mtor),ijklto(2,mtor),
     1              ijklto(3,mtor),ijklto(4,mtor)
            enddo
         else
            write(iw,9992)
         endif
      endif
c     
c     ----- eliminate duplicates and clean up -----
c     
      if(ntor.ge.2) then
         do itor=1,ntor-1
            ii=ijklto(1,itor)
            ij=ijklto(2,itor)
            ik=ijklto(3,itor)
            il=ijklto(4,itor) 
            if(ijklto(1,itor).ne.0) then
               do jtor=itor+1,ntor
                  ji=ijklto(1,jtor)
                  jj=ijklto(2,jtor)
                  jk=ijklto(3,jtor)
                  jl=ijklto(4,jtor)
                  if(jj.eq.ij.and.jk.eq.ik) then
                     if(ii.eq.ji.and.il.eq.jl) then
                        ijklto(1,jtor)=0
                        ijklto(2,jtor)=0
                        ijklto(3,jtor)=0
                        ijklto(4,jtor)=0
                     endif
                  elseif(jj.eq.ik.and.jk.eq.ij) then
                     if(ii.eq.jl.and.il.eq.ji) then
                        ijklto(1,jtor)=0
                        ijklto(2,jtor)=0
                        ijklto(3,jtor)=0
                        ijklto(4,jtor)=0
                     endif
                  endif
               enddo
            endif
         enddo
      endif
c     
      itor=0
 100  itor=itor+1
      if(itor.gt.ntor) go to 110
      if(ijklto(1,itor).eq.0) then
         if(itor.lt.ntor) then
            do jtor=itor+1,ntor
               ijklto(1,jtor-1)=ijklto(1,jtor)
               ijklto(2,jtor-1)=ijklto(2,jtor)
               ijklto(3,jtor-1)=ijklto(3,jtor)
               ijklto(4,jtor-1)=ijklto(4,jtor)
            enddo
            itor=itor-1
         endif
         ntor=ntor-1
      else
         rijsq=distsq(ijklto(1,itor),ijklto(2,itor))
         rij  =sqrt(rijsq)
         rjksq=distsq(ijklto(2,itor),ijklto(3,itor))
         rjk  =sqrt(rjksq)
         rklsq=distsq(ijklto(3,itor),ijklto(4,itor))
         rkl  =sqrt(rklsq)
         riksq=distsq(ijklto(1,itor),ijklto(3,itor))
         rjlsq=distsq(ijklto(2,itor),ijklto(4,itor))
         cosbj=(rijsq+rjksq-riksq)/(two*rij*rjk)
         cosbk=(rjksq+rklsq-rjlsq)/(two*rjk*rkl)
         linijk = sin(acos(min(1d0,abs(cosbj)))) .le. tol
         linjkl = sin(acos(min(1d0,abs(cosbk)))) .le. tol
c$$$  dum  = abs(cosbj)-one
c$$$  if(dum.gt.zero) then
c$$$  if(dum.le.tol) then
c$$$  linijk=.true.
c$$$  else
c$$$  write(iw,*) 'something is wrong in -zmtyp3-'
c$$$  call hnd_hnderr(3,errmsg)
c$$$  endif
c$$$  else
c$$$  dum=-dum
c$$$  if(dum.le.tol) then
c$$$  linijk=.true.
c$$$  else
c$$$  linijk=.false.
c$$$  endif
c$$$  endif
c$$$  dum  =abs(cosbk)-one
c$$$  if(dum.gt.zero) then
c$$$  if(dum.le.tol) then
c$$$  linjkl=.true.
c$$$  else
c$$$  write(iw,*) 'something is wrong in -zmtyp3-'
c$$$  call hnd_hnderr(3,errmsg)
c$$$  endif
c$$$  else
c$$$  dum=-dum
c$$$  if(dum.le.tol) then
c$$$  linjkl=.true.
c$$$  else
c$$$  linjkl=.false.
c$$$  endif
c$$$  endif
         if(.not.linijk.and..not.linjkl) then
            if (nizmat+5 .gt. mxizmt) goto 8000
            izmat(nizmat+1)=ithree
            izmat(nizmat+2)=ijklto(1,itor)
            izmat(nizmat+3)=ijklto(2,itor)
            izmat(nizmat+4)=ijklto(3,itor)
            izmat(nizmat+5)=ijklto(4,itor)
            nizmat=nizmat+5
            nzvar =nzvar +1
         else
*     write(LuOut,303) ijklto(1,itor), ijklto(2,itor), ijklto(3,itor), 
*     $           linijk, ijklto(2,itor), ijklto(3,itor), ijklto(4,itor),
*     $           linjkl
*     303        format(' rejecting torsion due to linear angle ',
*     $           2(2x,3i4,2x,l1))
         endif
      endif
      go to 100
 110  continue
c     
      if(dbug) then
         if(ntor.gt.1) then
            write(iw,9994)
            do mtor=1,ntor
               write(iw,9993) mtor,ijklto(1,mtor),ijklto(2,mtor),
     1              ijklto(3,mtor),ijklto(4,mtor)
            enddo
         else
            write(iw,9992)
         endif
      endif
c     
      if(dbug) then
         write(iw,9998)
         write(iw,9997) (izmat(i),i=1,nizmat)
         write(iw,9996)
         write(iw,*)    'nzvar = ',nzvar
      endif
c     
      zdone = .true.
      return
c     
 8000 if (oprint) write(LuOut,*) ' autoz: too many torsions '
      zdone = .false.
c     
 9999 format(' in -zmtyp3- ',/,' ----------- ')
 9998 format(' $zmat  ',/,' izmat= ')
 9997 format(12(i4,1h,))
 9996 format(' $end   ')
 9995 format(' something wrong with bond ',i2,' of atom ',i4,
     1     ' with -icon-, -jcon- = ',2i5)
 9994 format(' table of torsions = ')
 9993 format(5i5)
 9992 format(' no torsions found ')
      end
c$$$      SUBROUTINE HND_ZMTYP4(NZVAR,NIZMAT,IZMAT,
c$$$     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
c$$$     2                                       IJKLOP,NUMOOP)
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      PARAMETER (IFOUR=4)
c$$$      PARAMETER (MXATOM=500)        
c$$$      PARAMETER (MXBNDS= 8*MXATOM)
c$$$      PARAMETER (MXOOPA= 8*MXATOM)
c$$$      LOGICAL     DBUG
c$$$      CHARACTER*8 ERRMSG
c$$$      COMMON/HND_IOFILE/IR,IW
c$$$      DIMENSION IZMAT(*)   
c$$$      DIMENSION IJKLOP(4,*)      
c$$$      DIMENSION IJBOND(2,MXBNDS)
c$$$      DIMENSION IIBOND(  MXATOM)
c$$$      DIMENSION NIBOND(  MXATOM)
c$$$      DIMENSION IJBNDS(2,*)
c$$$      DIMENSION IIBNDS(  *)
c$$$      DIMENSION NIBNDS(  *)
c$$$      DIMENSION ERRMSG(3)        
c$$$      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP4-'/
c$$$C
c$$$      DBUG=.FALSE.
c$$$      IF(DBUG) THEN
c$$$         WRITE(IW,9999)
c$$$      ENDIF
c$$$C
c$$$C     ----- OUT-OF-PLANE-BENDS -----
c$$$C
c$$$      NBOND=NBNDS
c$$$      DO IJBND=1,NBOND
c$$$         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
c$$$         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
c$$$      ENDDO
c$$$      DO IAT=1,NAT
c$$$         IIBOND(IAT)=IIBNDS(IAT)
c$$$         NIBOND(IAT)=NIBNDS(IAT)
c$$$      ENDDO
c$$$C
c$$$      IF(NAT.LE.3) RETURN
c$$$C
c$$$      NOOP=0
c$$$      DO IAT=1,NAT
c$$$         DO IJBND=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
c$$$            ICON=IJBOND(1,IJBND)
c$$$            JCON=IJBOND(2,IJBND)
c$$$            IF(ICON.NE.IAT) THEN
c$$$               WRITE(IW,9995) IJBND,IAT,ICON,JCON
c$$$               CALL HND_HNDERR(3,ERRMSG)
c$$$            ENDIF
c$$$            JAT=JCON
c$$$            IF(NIBOND(JAT).GT.2) THEN
c$$$               DO JKBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
c$$$                  JCON=IJBOND(1,JKBND)
c$$$                  KCON=IJBOND(2,JKBND)
c$$$                  IF(JCON.NE.JAT) THEN
c$$$                     WRITE(IW,9995) JKBND,JAT,JCON,KCON
c$$$                     CALL HND_HNDERR(3,ERRMSG)
c$$$                  ENDIF
c$$$                  KAT=KCON
c$$$                  IF(KAT.NE.IAT) THEN
c$$$                     DO JLBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
c$$$                        JCON=IJBOND(1,JLBND)
c$$$                        LCON=IJBOND(2,JLBND)
c$$$                        IF(JCON.NE.JAT) THEN
c$$$                           WRITE(IW,9995) JLBND,JAT,JCON,LCON
c$$$                           CALL HND_HNDERR(3,ERRMSG)
c$$$                        ENDIF
c$$$                        LAT=LCON
c$$$                        IF(LAT.NE.IAT.AND.LAT.NE.KAT) THEN
c$$$                           NKCON=NIBOND(KAT)
c$$$                           NLCON=NIBOND(LAT)
c$$$                           IF((NKCON.EQ.1.AND.NLCON.EQ.1).OR.
c$$$     1                        (NKCON.GT.1.AND.NLCON.GT.1)    ) THEN
c$$$                              IF(NOOP.EQ.MXOOPA) THEN
c$$$                              WRITE(IW,*) 'TOO MANY OUT-OF-PLANE BENDS',
c$$$     1                                    ' NOOP, MXOOPA = ',NOOP,MXOOPA
c$$$                                 CALL HND_HNDERR(3,ERRMSG)
c$$$                              ENDIF
c$$$                              NOOP=NOOP+1
c$$$                              IJKLOP(1,NOOP)=IAT
c$$$                              IJKLOP(2,NOOP)=JAT
c$$$                              IJKLOP(3,NOOP)=KAT
c$$$                              IJKLOP(4,NOOP)=LAT
c$$$                           ENDIF
c$$$                        ENDIF
c$$$                     ENDDO
c$$$                  ENDIF
c$$$               ENDDO
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDDO
c$$$      IF(DBUG) THEN
c$$$         IF(NOOP.GT.1) THEN
c$$$            WRITE(IW,9994)
c$$$            DO MOOP=1,NOOP
c$$$               WRITE(IW,9993) MOOP,IJKLOP(1,MOOP),IJKLOP(2,MOOP),
c$$$     1                             IJKLOP(3,MOOP),IJKLOP(4,MOOP)
c$$$            ENDDO
c$$$         ELSE
c$$$            WRITE(IW,9992)
c$$$         ENDIF
c$$$      ENDIF
c$$$C
c$$$C     ----- ELIMINATE DUPLICATES AND CLEAN UP -----
c$$$C
c$$$      IF(NOOP.GE.2) THEN
c$$$         DO IOOP=1,NOOP-1
c$$$            II=IJKLOP(1,IOOP)
c$$$            IJ=IJKLOP(2,IOOP)
c$$$            IK=IJKLOP(3,IOOP)
c$$$            IL=IJKLOP(4,IOOP) 
c$$$            IF(IJKLOP(1,IOOP).NE.0) THEN
c$$$               DO JOOP=IOOP+1,NOOP
c$$$                  JI=IJKLOP(1,JOOP)
c$$$                  JJ=IJKLOP(2,JOOP)
c$$$                  JK=IJKLOP(3,JOOP)
c$$$                  JL=IJKLOP(4,JOOP)
c$$$                  IF(JI.EQ.II.AND.JJ.EQ.IJ) THEN
c$$$                     IF((JK.EQ.IK.AND.JL.EQ.IL).OR.
c$$$     1                  (JK.EQ.IL.AND.JL.EQ.IK)    ) THEN
c$$$                        IJKLOP(1,JOOP)=0
c$$$                        IJKLOP(2,JOOP)=0
c$$$                        IJKLOP(3,JOOP)=0
c$$$                        IJKLOP(4,JOOP)=0
c$$$                     ENDIF
c$$$                  ENDIF
c$$$               ENDDO
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDIF
c$$$C
c$$$      IOOP=0
c$$$  100 IOOP=IOOP+1
c$$$      IF(IOOP.GT.NOOP) GO TO 110
c$$$         IF(IJKLOP(1,IOOP).EQ.0) THEN
c$$$            IF(IOOP.LT.NOOP) THEN
c$$$               DO JOOP=IOOP+1,NOOP
c$$$                  IJKLOP(1,JOOP-1)=IJKLOP(1,JOOP)
c$$$                  IJKLOP(2,JOOP-1)=IJKLOP(2,JOOP)
c$$$                  IJKLOP(3,JOOP-1)=IJKLOP(3,JOOP)
c$$$                  IJKLOP(4,JOOP-1)=IJKLOP(4,JOOP)
c$$$               ENDDO
c$$$               IOOP=IOOP-1
c$$$            ENDIF
c$$$            NOOP=NOOP-1
c$$$         ELSE
c$$$             IZMAT(NIZMAT+1)=IFOUR 
c$$$             IZMAT(NIZMAT+2)=IJKLOP(1,IOOP)
c$$$             IZMAT(NIZMAT+3)=IJKLOP(2,IOOP)
c$$$             IZMAT(NIZMAT+4)=IJKLOP(3,IOOP)
c$$$             IZMAT(NIZMAT+5)=IJKLOP(4,IOOP)
c$$$            NIZMAT=NIZMAT+5
c$$$            NZVAR =NZVAR +1
c$$$         ENDIF
c$$$         GO TO 100
c$$$  110 CONTINUE
c$$$C
c$$$      IF(DBUG) THEN
c$$$         IF(NOOP.GT.1) THEN
c$$$            WRITE(IW,9994)
c$$$            DO MOOP=1,NOOP
c$$$               WRITE(IW,9993) MOOP,IJKLOP(1,MOOP),IJKLOP(2,MOOP),
c$$$     1                             IJKLOP(3,MOOP),IJKLOP(4,MOOP)
c$$$            ENDDO
c$$$         ELSE
c$$$            WRITE(IW,9992)
c$$$         ENDIF
c$$$      ENDIF
c$$$C
c$$$      IF(DBUG) THEN
c$$$         WRITE(IW,9998)
c$$$         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
c$$$         WRITE(IW,9996)
c$$$         WRITE(IW,*)    'NZVAR = ',NZVAR
c$$$      ENDIF
c$$$C
c$$$      RETURN
c$$$ 9999 FORMAT(' IN -ZMTYP4- ',/,' ----------- ')
c$$$ 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
c$$$ 9997 FORMAT(12(I4,1H,))
c$$$ 9996 FORMAT(' $END   ')
c$$$ 9995 FORMAT(' SOMETHING WRONG WITH BOND ',I2,' OF ATOM ',I4,
c$$$     1       ' WITH -ICON-, -JCON- = ',2I5)
c$$$ 9994 FORMAT(' TABLE OF OUT-OF-PLANE-BENDS = ')
c$$$ 9993 FORMAT(5I5)
c$$$ 9992 FORMAT(' NO OUT-OF-PLANE-ANGLE FOUND ')
c$$$      END
c$$$      SUBROUTINE HND_ZMTYP5(NZVAR,NIZMAT,IZMAT,C,
c$$$     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
c$$$     2                      IJKLNB,NUMLNB,XYZLNB,IJKANG,NUMANG)
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      PARAMETER (IFOUR=4)
c$$$      PARAMETER (MXATOM=500)        
c$$$      PARAMETER (MXBNDS= 8*MXATOM)
c$$$      PARAMETER (MXOOPA= 8*MXATOM)
c$$$      LOGICAL     DBUG
c$$$      CHARACTER*8 ERRMSG
c$$$      COMMON/HND_IOFILE/IR,IW
c$$$      DIMENSION C(3,*)
c$$$      DIMENSION IZMAT(*)   
c$$$      DIMENSION IJKLNB(4,*)      
c$$$      DIMENSION XYZLNB(3,*)
c$$$      DIMENSION IJKANG(3,*)           
c$$$      DIMENSION IJBNDS(2,*)
c$$$      DIMENSION IIBNDS(  *)
c$$$      DIMENSION NIBNDS(  *)
c$$$      DIMENSION ERRMSG(3)        
c$$$      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP5-'/
c$$$C
c$$$      DBUG=.FALSE.
c$$$      IF(DBUG) THEN
c$$$         WRITE(IW,9999)
c$$$         WRITE(IW,9991)
c$$$      ENDIF
c$$$C
c$$$C     ----- LINEAR-BENDS -----
c$$$C
c$$$      NLBND=0
c$$$C
c$$$      IF(NAT.LE.3) RETURN
c$$$C
c$$$      IF(DBUG) THEN
c$$$         IF(NLBND.GT.1) THEN
c$$$            WRITE(IW,9994)
c$$$            DO MLBND=1,NLBND
c$$$               WRITE(IW,9993) MLBND,IJKLNB(1,MLBND),IJKLNB(2,MLBND),
c$$$     1                              IJKLNB(3,MLBND),IJKLNB(4,MLBND)
c$$$            ENDDO
c$$$         ELSE
c$$$            WRITE(IW,9992)
c$$$         ENDIF
c$$$      ENDIF
c$$$C
c$$$      IF(DBUG) THEN
c$$$         WRITE(IW,9998)
c$$$         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
c$$$         WRITE(IW,9996)
c$$$      ENDIF
c$$$C
c$$$      RETURN
c$$$ 9999 FORMAT(' IN -ZMTYP5- ',/,' ----------- ')
c$$$ 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
c$$$ 9997 FORMAT(12(I4,1H,))
c$$$ 9996 FORMAT(' $END   ')
c$$$ 9994 FORMAT(' TABLE OF LINEAR-BENDS = ')
c$$$ 9993 FORMAT(5I5)
c$$$ 9992 FORMAT(' NO LINEAR-BENDS FOUND ')
c$$$ 9991 FORMAT(' LINEAR-BENDS NOT SET UP CURRENTLY . STOP . ')
c$$$      END
      subroutine hnd_dparsc(a,la,c,lc)
c     implicit double precision (a-h,o-z)
      implicit none
      integer la,lc,ic
      character*(*) a
      character*(*) c
      character*1   blk
      data  blk /' '/
      do ic=1,lc
         c(ic:ic)=blk
      enddo
      if(lc.le.(la-2)) then
         c(1:lc)=a(2:lc+1)
      else
         c(1:lc)=a(2:la-1)
      endif
      return
      end
      subroutine hnd_dparsi(a,la,n)
c     implicit double precision (a-h,o-z)
      implicit none
      integer la,i1,i2,isign,n,ia,ib,i
      character*(*) a
      character*1   char(12)
      data char     /'0','1','2','3','4','5','6','7','8','9',
     1               '+','-'/
c
      i1=1
      i2=la
c
c     ----- sign -----
c
      isign=1
      if(a(i1:i1).eq.char(12)) then
         isign=-1
         i1=i1+1
      elseif(a(i1:i1).eq.char(11)) then
         isign= 1
         i1=i1+1
      else
         isign= 1
      endif
c     na=i2-i1+1
c
      n=0
      do ia=i1,i2
         ib=i2-ia
         do i=1,10
            if(a(ia:ia).eq.char(i)) then
               n=n+(i-1)*10**(ib)
            endif
         enddo
      enddo
      n=n*isign
c
      return
      end
      subroutine hnd_dparsr(a,la,x)
c     implicit double precision (a-h,o-z)
      implicit none
      integer la,i1,i2,isign,ie1,ie2,ie,iexp,itmp,i,j
      double precision dum,zero,ten,x
      logical rep
      character*(*) a
      character*1   char(17)
      data char     /'0','1','2','3','4','5','6','7','8','9',
     1               '+','-','.','e','d','E','D'/
      data zero     /0.0d+00/
      data ten      /1.0d+01/
c     
      i1 =1            
      i2 =la
c
c     ----- sign -----
c
      isign=1
      if(a(i1:i1).eq.char(12)) then
         isign=-1
         i1=i1+1
      elseif(a(i1:i1).eq.char(11)) then
         isign= 1
         i1=i1+1
      else
         isign= 1
      endif
c
c     ----- exponent -----
c
      ie2=i2
      do ie=i1+1,i2
         if(a(ie:ie).eq.char(14).or.a(ie:ie).eq.char(15).or. 
     1      a(ie:ie).eq.char(16).or.a(ie:ie).eq.char(17)) go to 10
      enddo
      iexp=0
      go to 50
c
 10   ie2=i2
      ie1=ie+1
      i2 =ie-1
      iexp=1
      if(a(ie1:ie1).eq.char(12)                           ) then
         iexp=-1
      endif
      if(a(ie1:ie1).eq.char(12).or.a(ie1:ie1).eq.char(11)) then
         ie1=ie1+1
      endif
      itmp=0
      do i=ie1,ie2
         do j=1,10
            if(a(i:i).eq.char(j)) go to 30
         enddo
         go to 100
 30      itmp=itmp*10+j-1
      enddo
      iexp=iexp*itmp
c
c     ----- the number itself -----
c
 50   continue
      rep=.false.
      dum=zero
      do i=i1,i2
         if(a(i:i).ne.char(13)) then
            do j=1,10
               if(a(i:i).eq.char(j)) go to 70
            enddo
            go to 100
 70         dum=dum*ten+dble(j-1)
         else
            if(rep) go to 100
            iexp=iexp+i-i2
            rep=.true.
         endif
      enddo
      dum=dum*dble(isign)*ten**iexp
      x  =dum
      return
c
 100  continue       
      return
      end
      SUBROUTINE HND_HNDERR(LERR,ERRMSG)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
      integer lerr,iw
      CHARACTER*8  ERRMSG
      DIMENSION    ERRMSG(LERR)
c
      iw = LuOut
C
      WRITE(IW,9999)
      WRITE(IW,9998) ERRMSG
      WRITE(IW,9997)
C
      CALL ERRQUIT(' CALLS IT QUIT FROM HND_HNDERR ',0,0)
 9999 FORMAT(/,1X,72(1H-),/,31X,11HJOB STOPPED)
 9998 FORMAT(/,1X,9A8)
 9997 FORMAT(/,1X,72(1H-))
      END
      SUBROUTINE HND_GEOCLS(NFT)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "stdio.fh"
      integer nft,iw
      CHARACTER*8  ERRMSG
      DIMENSION    ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-GEOCLS-'/
      IW=LuOut
      CLOSE(UNIT=NFT,STATUS='KEEP',ERR=10)
      RETURN
   10 WRITE(IW,9999) NFT
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(/,' ----- ERROR CLOSING UNIT ',I3,' IN -GEOCLS- . STOP .')
      END
      SUBROUTINE HND_GEOOPN(NFT,GEOFIL)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "stdio.fh"
      integer nft,iw
      CHARACTER*80 GEOFIL
      CHARACTER*8  ERRMSG
      DIMENSION    ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-GEOOPN-'/
      IW = LuOut
      OPEN(UNIT=NFT,FILE=GEOFIL,STATUS='OLD',
     1     ACCESS='SEQUENTIAL',FORM=  'FORMATTED',ERR=10)
      RETURN
   10 WRITE(IW,9999) NFT
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(/,' ----- ERROR OPENING FILE -',A80,'-',/,
     1         '       AS UNIT ',I3,' IN -GEOOPN- . STOP.')
      END
      SUBROUTINE HND_ZDAT(IAT,NAT,NVAR,ZVAL,ZLST,
     1                ZMT,NUMZMT,PRSZMT,FLGZMT,ZMTCHR,
     2                NUMVAR,PRSVAR,VARCHR,FRZVAR,FRZVAL,LST,
     3                IZMAT,IZ,IZFRZ,MXIZMT,NZMOD,DBUG,
     $     zvarname, zvarsign)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
C
C     THIS ROUTINE TAKES THE CHARACTER STRING VALUES FOR THE
C     -Z- MATRIX INPUT FOR THE INDICES, AND THE VALUES AND TRANSFORMS
C     THEM INTO INTEGER AND REAL VALUES. THE CHARACTER STRING
C     VALUES CAN BE EITHER NUMERIC (INTEGER AND REAL), OR A VARIABLE
C     THAT GETS REPLACED BY A VALUE FOR THE VARIABLE.
C
C     INPUT PARAMETERS:
C        IAT     - INTEGER    -  ATOM NUMBER BEING CONSIDERED
C        NAT     - INTEGER    -  TOTAL NUMBER OF ATOMS
C        NVAR    - INTEGER    -  TOTAL NUMBER OF VARIABLES
C        NUMZMT  - INTEGER    -  ARRAY OF NUMBERS OF WORDS PER LINE
C        PRSZMT  - CHARACTER*80 -  ARRAY OF PARSED -Z- DATA
C        FLGZMT  - INTEGER    -  ARRAY OF TYPES FOR -Z- DATA
C                                  0 ... BLANK
C                                  1 ... FLOATING
C                                  2 ... INTEGER
C                                  3 ... ALPHANUMERIC
C        ZMTCHR  - INTEGER    -  ARRAY OF LENGTHS FOR -Z- DATA
C        NUMVAR  - INTEGER    -  ARRAY OF NUM. OF PIECES OF INFO
C                                  FOR EACH VARIABLE FOR SUBSTITUTION
C        PRSVAR  - CHARACTER*80 -  ARRAY OF VARIABLES FOR SUBSTITUTION
C        VARCHR  - INTEGER    -  ARRAY OF LENGTHS FOR VARIABLES
C
C     OUTPUT PARAMETERS:
C        ZVAL    - REAL*8       -  ARRAY OF -Z- MATRIX VALUES
C        ZMT     - INTEGER    -  ARRAY OF -Z- MATRIX INDICES
c        zvarname- char*8     -  array of -z- matrix variable names
c     .                          or blank if numeric value given.
c     .  zvarsign- real*8     -  array of -z-matrix variables signs
c     .                          (+1 or -1 as input, 1 if numeric given)
C
C     S. CHIN: 11/08/90 - IBM KINGSTON, NY
C
#include "stdio.fh"
#include "nwc_const.fh"
      integer       IAT,NAT,NVAR,NZMOD,IZ
      integer       MXATOM,MAXGEO,MAXZMT,MAXVAR,MXIZMT
      PARAMETER     (MXATOM=nw_max_atom)
      PARAMETER     (MAXGEO=MXATOM+1,MAXZMT=40,MAXVAR=nw_max_zmat)
      LOGICAL       DBUG
      LOGICAL       LST
      LOGICAL       READY
      CHARACTER*80  PRSZMT
      CHARACTER*80  PRSVAR
      LOGICAL       CART
      CHARACTER*8   ATNAME
      CHARACTER*1   PLUS
      CHARACTER*1   MINUS
      CHARACTER*8   ERRMSG
      double precision xx,yy,zz,atnum
      COMMON/HND_XYZGEO/XX(MAXGEO),YY(MAXGEO),ZZ(MAXGEO),
     1                  ATNAME(MAXGEO),ATNUM(MAXGEO),CART(MAXGEO)
      DIMENSION         ERRMSG(3)
      integer           NUMZMT(       MAXGEO)
      DIMENSION         PRSZMT(MAXZMT,MAXGEO)
      integer           FLGZMT(MAXZMT,MAXGEO)
      integer           ZMTCHR(MAXZMT,MAXGEO)
      DIMENSION         PRSVAR(MAXZMT,MAXVAR)
      integer           VARCHR(MAXZMT,MAXVAR)
      logical           FRZVAR(       MAXVAR)
      integer           NUMVAR(       MAXVAR)
      logical           FRZVAL(     3,MAXGEO)
      double precision    ZVAL(     3,MAXGEO)
      integer              ZMT(     5,MAXGEO)
      double precision    ZLST(     3,MAXGEO)
      integer            IZMAT(MXIZMT)
      logical            IZFRZ(MXIZMT)

      integer          i,j,k,JVAR,JAT
      character*(*) zvarname(*)
      double precision zvarsign(*)

      DATA ERRMSG /'PROGRAM ','STOP IN ','- ZDAT -'/
      double precision ZERO,TWO,PIDEG
      integer          iw
      DATA ZERO   /0.0D+00/
      DATA TWO    /2.0D+00/
      DATA PIDEG  /180.0D+00/
      DATA PLUS   /'+'/
      DATA MINUS  /'-'/
c
      IW = luout
c
      do i = nzmod+1,nzmod+3
         zvarname(i) = ' '       ! Default for numeric parameter
         zvarsign(i) = 1d0       ! Default for numeric parameter
      end do
C
C     ---- THIS ATOM IS BEING INPUTED WITH CARTESIAN COORDINATES ----
C
      IF(FLGZMT(3,IAT).EQ.1.AND.FLGZMT(4,IAT).EQ.1) THEN
             IF(NUMZMT(  IAT).EQ.4.AND.FLGZMT(1,IAT).EQ.3.AND.
     1          FLGZMT(2,IAT).EQ.1                            ) THEN
c           READ(PRSZMT(1,IAT),*) ATNAME(IAT)
c           READ(PRSZMT(2,IAT),*) XX(IAT)
c           READ(PRSZMT(3,IAT),*) YY(IAT)
c           READ(PRSZMT(4,IAT),*) ZZ(IAT)
            call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
            call hnd_dparsr(prszmt(2,iat),zmtchr(2,iat),xx(iat)      )
            call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),yy(iat)      )
            call hnd_dparsr(prszmt(4,iat),zmtchr(4,iat),zz(iat)      )
                 ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
         ELSEIF(NUMZMT(  IAT).EQ.5.AND.FLGZMT(1,IAT).EQ.3.AND.
     1          FLGZMT(2,IAT).EQ.2.AND.FLGZMT(5,IAT).EQ.1     ) THEN
c           READ(PRSZMT(1,IAT),*) ATNAME(IAT)
c           READ(PRSZMT(3,IAT),*) XX(IAT)
c           READ(PRSZMT(4,IAT),*) YY(IAT)
c           READ(PRSZMT(5,IAT),*) ZZ(IAT)
            call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
            call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),xx(iat)      )
            call hnd_dparsr(prszmt(4,iat),zmtchr(4,iat),yy(iat)      )
            call hnd_dparsr(prszmt(5,iat),zmtchr(5,iat),zz(iat)      )
                 ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
         ELSE
            WRITE(IW,9994) IAT
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         CART(IAT) = .TRUE.
         RETURN
      ENDIF
C
C     ----- CONSTRUCT THE -Z- MATRIX -----
C
      IF(IAT.EQ.1.AND.(NUMZMT(IAT).NE.1.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9989)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IAT.EQ.2.AND.(NUMZMT(IAT).NE.3.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9988)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IAT.EQ.3.AND.(NUMZMT(IAT).NE.5.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9987)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- ANCHOR ATOM              -----
C
      ZMT(1,IAT) = IAT
      IF(FLGZMT(1,IAT) .EQ. 3) THEN
c        READ (PRSZMT(1,IAT),*) ATNAME(IAT)
         call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
               ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
      ELSE
         WRITE(IW,9999)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 1) GO TO 100
C
C     ----- BOND LENGTH -----
C
      IF(FLGZMT(2,IAT) .EQ. 2) THEN
c        READ (PRSZMT(2,IAT),*) ZMT(2,IAT)
         call hnd_dparsi(prszmt(2,iat),zmtchr(2,iat),zmt(2,iat)   )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(2,IAT) (1:8) .EQ. PRSZMT(1,JAT) (1:8)) THEN
               ZMT(2,IAT)=JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(3,IAT) .EQ. 1) THEN
c        READ (PRSZMT(3,IAT),*) ZVAL(1,IAT)
         call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),zval(1,iat)  )
                                ZLST(1,IAT)=ZVAL(1,IAT)
      ELSEIF(FLGZMT(3,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(3, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(3, IAT) (1:ZMTCHR(3, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))   ) THEN
               zvarname(nzmod+1) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
c              READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                       ZLST(1,IAT)=ZVAL(1,IAT)
               FRZVAL(1,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(3, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(3, IAT) (3:ZMTCHR(3, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))   ) THEN
               IF(PRSZMT(3, IAT) (2:2).EQ.PLUS) THEN
                  zvarname(nzmod+1) = PRSVAR(1,JVAR)(2:VARCHR(1,JVAR)-1)
c                 READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                          ZLST(1,IAT)=ZVAL(1,IAT)
                  FRZVAL(1,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(3, IAT) (2:2).EQ.MINUS) THEN
                  zvarname(nzmod+1) = PRSVAR(1,JVAR)(2:VARCHR(1,JVAR)-1)
                  zvarsign(nzmod+1) = -1d0
c                 READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                          ZLST(1,IAT)=ZVAL(1,IAT)
                  FRZVAL(1,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
                  ENDIF
                  ZVAL(1,IAT)=-ZVAL(1,IAT)
                  ZLST(1,IAT)=-ZLST(1,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9998) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9998) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 2) GO TO 100
C
C     ----- BOND ANGLE -----
C
      IF(FLGZMT(4,IAT) .EQ. 2) THEN
c        READ (PRSZMT(4,IAT),*) ZMT(3,IAT)
         call hnd_dparsi(prszmt(4,iat),zmtchr(4,iat),zmt(3,iat)  )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(4,IAT) (1:8).EQ.PRSZMT(1,JAT) (1:8)) THEN
               ZMT(3,IAT) = JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(5,IAT) .EQ. 1) THEN
c        READ (PRSZMT(5,IAT),*) ZVAL(2,IAT)
         call hnd_dparsr(prszmt(5,iat),zmtchr(5,iat),zval(2,iat)  )
                                ZLST(2,IAT)=ZVAL(2,IAT)
      ELSEIF(FLGZMT(5,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(5, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(5, IAT) (1:ZMTCHR(5, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))    ) THEN
               zvarname(nzmod+2) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
c              READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                       ZLST(2,IAT)=ZVAL(2,IAT)
               FRZVAL(2,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(5, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(5, IAT) (3:ZMTCHR(5, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))) then
               IF(PRSZMT(5, IAT) (2:2).EQ.PLUS) THEN
                  zvarname(nzmod+2) = PRSVAR(1,JVAR)(2:VARCHR(1,JVAR)-1)
c                 READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                          ZLST(2,IAT)=ZVAL(2,IAT)
                  FRZVAL(2,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(5, IAT) (2:2).EQ.MINUS) THEN
                  zvarname(nzmod+2) = PRSVAR(1,JVAR)(2:VARCHR(1,JVAR)-1)
                  zvarsign(nzmod+2) = -1d0
c                 READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                          ZLST(2,IAT)=ZVAL(2,IAT)
                  FRZVAL(2,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
                  ENDIF
                  ZVAL(2,IAT)=-ZVAL(2,IAT)
                  ZLST(2,IAT)=-ZLST(2,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9997) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9997) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
c
c     this avoid weird rounding problems
c
      if(abs(zval(2,iat)-pideg).lt.1d-8) zval(2,iat)=pideg
      IF(ZVAL(2,IAT).LT.ZERO.OR.ZVAL(2,IAT).GT.PIDEG) THEN
         WRITE(IW,9986) ZVAL(2,IAT),IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 3) GO TO 100
C
C     ----- TORSION -----
C
      IF(FLGZMT(6,IAT) .EQ. 2) THEN
c        READ (PRSZMT(6,IAT),*) ZMT(4,IAT)
         call hnd_dparsi(prszmt(6,iat),zmtchr(6,iat),zmt(4,iat)  )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(6,IAT) (1:8).EQ.PRSZMT(1,JAT) (1:8)) THEN
               ZMT(4,IAT)=JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(7,IAT) .EQ. 1) THEN
c        READ (PRSZMT(7,IAT),*) ZVAL(3,IAT)
         call hnd_dparsr(prszmt(7,iat),zmtchr(7,iat),zval(3,iat)  )
                                ZLST(3,IAT)=ZVAL(3,IAT)
      ELSEIF(FLGZMT(7,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(7, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(7, IAT) (1:ZMTCHR(7, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))    ) THEN
c              READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
               zvarname(nzmod+3) = PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)-1)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                       ZLST(3,IAT)=ZVAL(3,IAT)
               FRZVAL(3,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(7, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(7, IAT) (3:ZMTCHR(7, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))    ) THEN
               IF(PRSZMT(7, IAT) (2:2).EQ.PLUS) THEN
                  zvarname(nzmod+3) = PRSVAR(1,JVAR)(2:VARCHR(1,JVAR)-1)
c                 READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                          ZLST(3,IAT)=ZVAL(3,IAT)
                  FRZVAL(3,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(7, IAT) (2:2).EQ.MINUS) THEN
                  zvarname(nzmod+3) = PRSVAR(1,JVAR)(2:VARCHR(1,JVAR)-1)
                  zvarsign(nzmod+3) = -1d0
c                 READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                          ZLST(3,IAT)=ZVAL(3,IAT)
                  FRZVAL(3,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
                  ENDIF
                  ZVAL(3,IAT)=-ZVAL(3,IAT)
                  ZLST(3,IAT)=-ZLST(3,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9996) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9996) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- LAST INTEGER PARAMETER -----
C
      IF(NUMZMT(IAT).EQ.8) THEN
         IF(FLGZMT(8,IAT).EQ.2) THEN
c           READ(PRSZMT(8,IAT),*) ZMT(5,IAT)
         call hnd_dparsi(prszmt(8,iat),zmtchr(8,iat),zmt(5,iat)  )
         ELSE
            WRITE(IW,9993) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
C     ----- POSSIBLE ERRORS -----
C
      IF(ZMT(5,IAT).EQ.0) THEN
         IF(ZVAL(3,IAT).LT.-TWO*PIDEG.OR.
     1      ZVAL(3,IAT).GT. TWO*PIDEG) THEN
            WRITE(IW,9985) ZVAL(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ZVAL(3,IAT).LT.-PIDEG) ZVAL(3,IAT)=ZVAL(3,IAT)+TWO*PIDEG
         IF(ZVAL(3,IAT).GT. PIDEG) ZVAL(3,IAT)=ZVAL(3,IAT)-TWO*PIDEG
C
         IF(ZLST(3,IAT).LT.-TWO*PIDEG.OR.
     1      ZLST(3,IAT).GT. TWO*PIDEG) THEN
            WRITE(IW,9985) ZLST(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZLST(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ZLST(3,IAT).LT.-PIDEG) ZLST(3,IAT)=ZLST(3,IAT)+TWO*PIDEG
         IF(ZLST(3,IAT).GT. PIDEG) ZLST(3,IAT)=ZLST(3,IAT)-TWO*PIDEG
      ELSE
         IF(ZVAL(3,IAT).LT. ZERO .OR.
     1      ZVAL(3,IAT).GT.PIDEG) THEN
            WRITE(IW,9985) ZVAL(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
C
         IF(ZLST(3,IAT).LT. ZERO .OR.
     1      ZLST(3,IAT).GT.PIDEG) THEN
            WRITE(IW,9985) ZLST(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZLST(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
  100 CONTINUE
      IF(DBUG) THEN
         WRITE(IW,9984) IAT,ATNAME(IAT)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(  ZVAL(I,IAT),I=1,3),
     1                                         (FRZVAL(I,IAT),I=1,3)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(  ZLST(I,IAT),I=1,3),
     1                                         (FRZVAL(I,IAT),I=1,3)
      ENDIF
      DO K=1,4
        IF (ZMT(K,IAT).NE.0) THEN
          DO J=1,K-1
            IF (ZMT(K,IAT).EQ.ZMT(J,IAT)) THEN
              WRITE(IW,9983)
              WRITE(IW,9984) IAT,ATNAME(IAT)
              WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3),
     1                           (FRZVAL(I,IAT),I=1,3)
              WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZLST(I,IAT),I=1,3),
     1                           (FRZVAL(I,IAT),I=1,3)
              CALL HND_HNDERR(3,ERRMSG)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C     ----- CREATE -IZMAT- FOR -HONDO- INPUT -----
C
      IF(IAT.GT.1) THEN
         NZMOD=NZMOD+1
         IZ=IZ+1
         IZMAT(IZ)=1
         IZ=IZ+1
         IZMAT(IZ)=ZMT(1,IAT)
         IZ=IZ+1
         IZMAT(IZ)=ZMT(2,IAT)
         IZFRZ(NZMOD)=FRZVAL(1,IAT)
C
         IF(IAT.GT.2) THEN
            NZMOD=NZMOD+1
            IZ=IZ+1
            IZMAT(IZ)=2
            IZ=IZ+1
            IZMAT(IZ)=ZMT(1,IAT)
            IZ=IZ+1
            IZMAT(IZ)=ZMT(2,IAT)
            IZ=IZ+1
            IZMAT(IZ)=ZMT(3,IAT)
            IZFRZ(NZMOD)=FRZVAL(2,IAT)
C
            IF(IAT.GT.3) THEN
               IF(ZMT(5,IAT).EQ.0) THEN
                  NZMOD=NZMOD+1
                  IZ=IZ+1
                  IZMAT(IZ)=3
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(1,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(2,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(3,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(4,IAT)
                  IZFRZ(NZMOD)=FRZVAL(3,IAT)
               ELSE
                  NZMOD=NZMOD+1
                  IZ=IZ+1
                  IZMAT(IZ)=7   ! RJH was 2 but discarded zmt(5,iat)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(1,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(2,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(4,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(5,IAT)
                  IZFRZ(NZMOD)=FRZVAL(3,IAT)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(IZ.GT.MXIZMT) THEN
         WRITE(IW,9991) IZ,MXIZMT,IAT
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      RETURN
 9999 FORMAT(' FIRST PARSED -Z- DATA IS NOT A CHARACTER STRING. STOP')
 9998 FORMAT(' THE 3-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9997 FORMAT(' THE 5-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9996 FORMAT(' THE 7-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9995 FORMAT(' IAT=',I5,' ZMAT=',5I3,3F12.5,3(2X,L4))
 9994 FORMAT(' SOMETHING IS WRONG WITH THE CARTESIAN COORDINATES INPUT',
     1       ' DATA FOR THIS ATOM, IAT = ',I3,' . STOP')
 9993 FORMAT(' THE 8-TH (INTEGER) PARAMETER FOR ATOM = ',I4,
     1       ' IS OF THE WRONG TYPE. STOP')
 9992 FORMAT(' ATOM REFERED TO FOR IAT = ',I3,' IS INVALID. STOP')
 9991 FORMAT(' TOO MANY -IZMAT- ENTRIES, IZ,MXIZMT = ',2I5,
     1       ' FOR -IAT- = ',I5,' STOP.')
 9989 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -1- IS INCORRECT. STOP')
 9988 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -2- IS INCORRECT. STOP')
 9987 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -3- IS INCORRECT. STOP')
 9986 FORMAT(' BOND ANGLE = ',F10.2,' FOR ATOM = ',I4,' IS OUT OF',
     1       ' ALLOWED RANGE ... 0.0 TO 180.0 . STOP')
 9985 FORMAT(' TORSION ANGLE OR BOND ANGLE = ',F10.2,' FOR ATOM = ',I4,
     1       ' IS OUT OF ALLOWED RANGE ... -360.0 TO 360.0 . STOP')
 9984 FORMAT(' IAT=',I5,' ATNAME=',1X,A8)
 9983 FORMAT(' SAME ATOM REFERRED TO MORE THAN ONCE ON THIS LINE')
      END
      SUBROUTINE HND_ZXYZ(NAT,NZMT,ZVAL)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
C
C     ROUTINE COMPUTES :
C
C     THE CARTESIAN COORDINATES OF A MOLECULE FROM BONDING COORDINATES
C     (BOND LENGTHS, BOND AND DIHEDRAL ANGLES).
C
C     ATOMS NA, NB, NC, HAVE KNOWN COORDINATES AND ARE NOT COLLINEAR.
C     RCD   =ZVAL(1,I) : C-D BOND LENGTH (IN ANGSTROM)
C     THBCD =ZVAL(2,I) : B-C-D BOND ANGLE (IN DEGREES)
C     PHABCD=ZVAL(3,I) : DIHEDRAL ANGLE (IN DEGREES) OF C-D
C                        RELATIVE TO A-B MEASURED CLOCKWISE
C                        ( + ) ALONG THE DIRECTION B TO C
C
C     MEANING OF PARAMETER NZMT(5,I):
C          =0     ANGLES FROM DATA
C          =+/- 1 ZVAL(3,I) IS BOND ANGLE A-C-D
C          =2     BOND ANGLE B-C-D IS TETRAHEDRAL
C          =3     ATOMS B-C-D ARE COLLINEAR
C
#include "stdio.fh"
#include "nwc_const.fh"
      integer       nat
      integer       MXATOM,MAXGEO,MAXWRD,MAXVAR
      PARAMETER     (MXATOM=nw_max_atom)
      PARAMETER     (MAXGEO=MXATOM+1,MAXWRD=40,MAXVAR=nw_max_zmat)
      LOGICAL       DBUG
      LOGICAL       CART
      CHARACTER*8   ATNAME
      CHARACTER*8   ERRMSG
      double precision xx,yy,zz,atnum
      COMMON/HND_XYZGEO/XX(MAXGEO),YY(MAXGEO),ZZ(MAXGEO),
     1                  ATNAME(MAXGEO),ATNUM(MAXGEO),CART(MAXGEO)
      integer            NZMT(5,MAXGEO)
      double precision   ZVAL(3,MAXGEO)
      DIMENSION   ERRMSG(3)
      double precision numd,zero,one,two,three,rcd,ccosp,ssinp
      double precision costd,sintd,ccos,ssin,ccost,ssint
      integer iw,na,nb,nc,nd,iat
      DATA ERRMSG    /'PROGRAM ','STOP IN ','- ZXYZ -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TWO,THREE /2.0D+00,3.0D+00/
      double precision t11,t12,t13,t21,t22,t23,t31,t32,t33
      double precision phi,bet,alp,gam,dum,dot,rca,pifac,rab,rcb
      double precision theta,thbcd
      double precision xa,xb,xxd,yyd,zzd,ya,yb,za,zb
C
      DBUG=.FALSE.
      IW = luout
C
      PIFAC=acos(-1.0d0)/180.0D+00
      COSTD=-ONE/THREE
      SINTD= TWO/THREE* SQRT(TWO)
      NA   =0
      NB   =0
      NC   =0
      ND   =0
      RCD  =ZERO
      CCOSP= 0d0
      SSIN= 0d0
      CCOS= 0d0
      SSINP= 0d0
      CCOST= 0d0
      SSINT= 0d0
C     
C     ----- ATOM -1- AND -2- -----
C
      IAT=1
      IF(CART(IAT)) THEN
         IF(DBUG) THEN
            WRITE(IW,9997) IAT
         ENDIF
      ELSE
         ND=NZMT(1,IAT)
         XX(ND)=ZERO
         YY(ND)=ZERO
         ZZ(ND)=ZERO
         CART(IAT)=.TRUE.
         IF(DBUG) THEN
            WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
         ENDIF
      ENDIF
C
      IF(NAT.GT.1) THEN
         IAT=2
         IF(CART(IAT)) THEN
            IF(DBUG) THEN
               WRITE(IW,9997) IAT
            ENDIF
         ELSE
            ND=NZMT(1,IAT)
            NC=NZMT(2,IAT)
            XX(ND)=ZERO
            YY(ND)=ZERO
            ZZ(ND)=ZVAL(1,IAT)
            RCD   =ZVAL(1,IAT)
            CART(IAT)=.TRUE.
            IF(DBUG) THEN
               WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
            ENDIF
         ENDIF
C
         IF(NAT.GT.2) THEN
C
C     ----- ATOM -3- -----
C
            IAT=3
            IF(CART(IAT)) THEN
               IF(DBUG) THEN
                  WRITE(IW,9997) IAT
               ENDIF
            ELSE
                   IF(NZMT(5,IAT).EQ.3) THEN
                  CCOS=-ONE
                  SSIN= ZERO
               ELSEIF(NZMT(5,IAT).EQ.2) THEN
                  CCOS=COSTD
                  SSIN=SINTD
               ELSEIF(NZMT(5,IAT).EQ.0) THEN
                  THBCD=ZVAL(2,IAT)*PIFAC
                  CCOS= COS(THBCD)
                  SSIN= SIN(THBCD)
               ELSE
                  WRITE(IW,9999) IAT
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               ND =NZMT(1,IAT)
               NC =NZMT(2,IAT)
               NB =NZMT(3,IAT)
               RCD=ZVAL(1,IAT)
               IF(DBUG) THEN
                  WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
               ENDIF
               XX(ND)=         RCD*SSIN
               YY(ND)=ZERO
               ZZ(ND)=ZZ(NC) - RCD*CCOS
               IF(ZZ(NC).LT.ZZ(NB)) THEN
                  ZZ(ND)=ZZ(NC) + RCD*CCOS
               ENDIF
               CART(IAT) = .TRUE.
            ENDIF
C
C     ----- ATOM -4- AND HIGHER -----
C
            IF(NAT.GT.3) THEN
               DO 20 IAT=4,NAT
                  IF(CART(IAT)) THEN
                     IF(DBUG) THEN
                        WRITE(IW,9997) IAT
                     ENDIF
                     GO TO 20
                  ENDIF
                  NA =NZMT(4,IAT)
                  NB =NZMT(3,IAT)
                  NC =NZMT(2,IAT)
                  ND =NZMT(1,IAT)
                  RCD=ZVAL(1,IAT)
                  IF(DBUG) THEN
                     WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
                  ENDIF
C
C     ----- LINEAR CONNECTION FOR THIS ATOM : D...C...B -----
C
                  IF(NZMT(5,IAT).EQ.3) THEN
                     RCB= SQRT((XX(NB)-XX(NC))**2+
     1                         (YY(NB)-YY(NC))**2+(ZZ(NB)-ZZ(NC))**2)
                     XX(ND)=XX(NC)+RCD*(XX(NC)-XX(NB))/RCB
                     YY(ND)=YY(NC)+RCD*(YY(NC)-YY(NB))/RCB
                     ZZ(ND)=ZZ(NC)+RCD*(ZZ(NC)-ZZ(NB))/RCB
                  ELSE
C
C     ----- OTHER CONNECTIONS FOR THIS ATOM : DEFINE LOCAL FRAME -----
C
                     XB=XX(NB)-XX(NC)
                     YB=YY(NB)-YY(NC)
                     ZB=ZZ(NB)-ZZ(NC)
                     RCB= SQRT(XB*XB+YB*YB+ZB*ZB)
                     XB=XB/RCB
                     YB=YB/RCB
                     ZB=ZB/RCB
                     T11=XB
                     T21=YB
                     T31=ZB
                     XA=XX(NA)-XX(NC)
                     YA=YY(NA)-YY(NC)
                     ZA=ZZ(NA)-ZZ(NC)
                     RCA= SQRT(XA*XA+YA*YA+ZA*ZA)
                     XA=XA/RCA
                     YA=YA/RCA
                     ZA=ZA/RCA
                     DOT=XA*XB+YA*YB+ZA*ZB
                     XA=XA-DOT*XB
                     YA=YA-DOT*YB
                     ZA=ZA-DOT*ZB
                     DUM= SQRT(XA*XA+YA*YA+ZA*ZA)
                     XA=XA/DUM
                     YA=YA/DUM
                     ZA=ZA/DUM
                     T12=XA
                     T22=YA
                     T32=ZA
                     T13=T21*T32-T31*T22
                     T23=T31*T12-T11*T32
                     T33=T11*T22-T21*T12
C
C     ----- GET CARTESIAN COORDINATES IN LOCAL  FRAME -----
C
                         IF(NZMT(5,IAT).EQ.0) THEN
                        PHI  = ZVAL(2,IAT)*PIFAC
                        THETA=-ZVAL(3,IAT)*PIFAC
                        CCOSP= COS(PHI)
                        SSINP= SIN(PHI)
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSEIF(NZMT(5,IAT).EQ.-1.OR.
     1                      NZMT(5,IAT).EQ. 1    ) THEN
                        RAB  =SQRT((XX(NB)-XX(NA))**2+(YY(NB)-YY(NA))**2
     1                            +(ZZ(NB)-ZZ(NA))**2)
                        PHI=ZVAL(2,IAT)*PIFAC
                        BET=ZVAL(3,IAT)*PIFAC
                        ALP=PHI
                        GAM  =ACOS((RCB**2+RCA**2-RAB**2)/(TWO*RCB*RCA))
                     numd=(COS(BET)- COS(ALP)* COS(GAM))/
     /                    (SIN(ALP)* SIN(GAM))
                     if(numd.gt.1d0) numd=1d0
                     if(numd.lt.-1d0) numd=-1d0
                       THETA=ACOS(numd)
                        IF(NZMT(5,IAT).EQ.-1) THEN
                           THETA=-THETA
                        ENDIF
                        CCOSP= COS(PHI)
                        SSINP= SIN(PHI)
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSEIF(NZMT(5,IAT).EQ.2) THEN
                        CCOSP= COSTD
                        SSINP= SINTD
                        THETA=ZVAL(3,IAT)*PIFAC
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSE
                        WRITE(IW,9999) IAT
                        CALL HND_HNDERR(3,ERRMSG)
                     ENDIF
                     XXD=RCD*CCOSP
                     YYD=RCD*SSINP*CCOST
                     ZZD=RCD*SSINP*SSINT
C
C     ----- GET CARTESIAN COORDINATES IN MASTER FRAME -----
C
                     XX(ND)=T11*XXD+T12*YYD+T13*ZZD + XX(NC)
                     YY(ND)=T21*XXD+T22*YYD+T23*ZZD + YY(NC)
                     ZZ(ND)=T31*XXD+T32*YYD+T33*ZZD + ZZ(NC)
                  ENDIF
                  CART(IAT)=.TRUE.
   20             CONTINUE
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
 9999 FORMAT(' ILLEGAL 5-TH -ZMT- PARAMETER FOR -IAT- = ',I3,' STOP.')
 9998 FORMAT(' IAT,NA,NB,NC,ND,RCD = ',5I4,F12.5)
 9997 FORMAT(' IAT = ',I4,' IS ALREADY SPECIFIED IN CARTESIAN SPACE.')
      END
      SUBROUTINE HND_PRSQ(V,M,N,NDIM)
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
#include "stdio.fh"
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C
***      COMMON/HND_LISTNG/LIST
      integer m,n,ndim
      double precision V(NDIM,1)
      integer IC(5)
      double precision C(5)
      double precision vtol
      integer icmax,iw,list,imax,imin,i,j,ii,idum,max
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
c
      IW = LuOut
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      IF(LIST.EQ.0) WRITE (IW,9048) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      WRITE(IW,9328) J
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
 9328 FORMAT(/,' VEC =',I5,/,1X,10(1H-))
 9348 FORMAT(5(I5,F11.5))
      END
      SUBROUTINE HND_PREV(V,E,M,N,NDIM)
c     IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
#include "stdio.fh"
C
C     ----- PRINT OUT E AND V-MATRICES
C
***      COMMON/HND_LISTNG/LIST
      integer m,n,ndim
      double precision V(NDIM,1),E(1)
      integer IC(5)
      double precision C(5)
      double precision vtol
      integer icmax
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
      integer iw,list,imax,imin,max,i,j,ii,idum
C
      IW = LuOut
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9068) (E(I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9168) (E(I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9268) (E(I),I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      IF(LIST.EQ.0) WRITE (IW,9048) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      WRITE(IW,9328) J,E(J)
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(15X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9068 FORMAT(15X,10F11.5)
 9128 FORMAT(15X,7(6X,I3,6X))
 9148 FORMAT(10X,I5,7F15.10)
 9168 FORMAT(15X,7F15.10)
 9228 FORMAT(15X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
 9268 FORMAT(15X,7E15.8)
 9328 FORMAT(/,' VEC =',I5,' EIG/OCC =',F11.5,/,1X,31(1H-))
 9348 FORMAT(5(I5,F11.5))
      END
      SUBROUTINE HND_PRTR(D,N)
c     IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
#include "stdio.fh"
C
C     ----- PRINT OUT A TRIANGULAR MATRIX -----
C
***      COMMON/HND_LISTNG/LIST
      integer n,iw,list,max,imax,imin,i,j,k,ii,jj,ij
      double precision D(1),DD(10)
C
      IW = LuOut
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 160 J = 1,N
      K = 0
      DO 140 I = IMIN,IMAX
      K = K+1
      II = MAX0( I, J)
      JJ = MIN0( I, J)
      IJ = (II*(II-1))/2 + JJ
  140 DD(K) = D(IJ)
      IF(LIST.EQ.0) WRITE (IW,9048) J,(DD(I),I = 1,K)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(DD(I),I = 1,K)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(DD(I),I = 1,K)
  160 CONTINUE
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
      END
      SUBROUTINE HND_TFTR(H,F,Q,T,IA,M,N,NDIM)
c     IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
#include "util.fh"
C
C     ----- H(M,M) = Q(DAGGER)(N,M) * F(N,N) * Q(N,M) -----
C
      integer m,n,ndim
      double precision H(1),F(1),Q(NDIM,1),T(1)
      integer IA(1)
      integer ij,j,ik,max,i,k
      double precision small,zero,dum,qij,hij
      double precision ddot
      external ddot
      DATA SMALL /1.0D-11/
      DATA ZERO  /0.0D+00/
      IJ = 0
      DO 180 J = 1,M
      IK = 0
      DO 140 I = 1,N
      DUM = ZERO
      QIJ = Q(I,J)
      MAX = I-1
      IF (MAX .EQ. 0) GO TO 120
      DO 100 K = 1,MAX
      IK = IK+1
      T(K) = T(K)+ F(IK)*QIJ
      DUM  = DUM + F(IK)*Q(K,J)
  100 CONTINUE
  120 IK = IK+1
      T(I) = DUM + F(IK)*QIJ
  140 CONTINUE
      DO 160 I = 1,J
      IJ = IJ+1
      HIJ = ddot(N,Q(1,I),1,T,1)
      IF( ABS(HIJ).LT.SMALL) HIJ=ZERO
      H(IJ)=HIJ
  160 CONTINUE
  180 CONTINUE
      RETURN
      END
      SUBROUTINE HND_DIAGIV(A,VVEC,EIG,IA,NVEC,N,NDIM)
c     IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
C
C     ----- GENERAL CALLING  ROUTINE FOR -----
C               -GIVEIS- OR -GIVDIA-
C
#include "mafdecls.fh"
#include "global.fh"
      logical status
      integer nvec,n,ndim
      double precision A(1),VVEC(NDIM,1),EIG(1)
      integer IA(1)
      integer i_idbl,idbl,i_iint,k_iint,iw
c
c     ----- memory -----
c
      IW = LuOut
      status=ma_push_get(MT_DBL,4*ndim,'mem_dbl dia',i_idbl,idbl)
      status=ma_push_get(MT_INT,3*ndim,'mem_int dia',i_iint,k_iint)
c
c     ----- call -----
c
      call hnd_givdia(a,vvec,eig,ia,n,ndim,
     1 dbl_mb(idbl),dbl_mb(idbl+ndim),
     . dbl_mb(idbl+2*ndim),dbl_mb(idbl+3*ndim),
     2                 dbl_mb(idbl+2*ndim),dbl_mb(idbl+2*ndim),
     3    int_mb(k_iint),int_mb(k_iint+ndim),int_mb(k_iint+2*ndim))
c
c     ----- memory -----
c
      status=ma_pop_stack(i_iint)
      if(.not.status) then
         call errquit('hnd_diagiv,ma_pop_stack of -iint- failed',911,
     &       MA_ERR)
      endif
      status=ma_pop_stack(i_idbl)
      if(.not.status) then
         call errquit('hnd_diagiv,ma_pop_stack of -idbl- failed',911,
     &       MA_ERR)
      endif
c
      RETURN
      END
      SUBROUTINE HND_GIVDIA(A,VEC,EIG,IA,N,NDIM,
     1           w,gamma,beta,betasq,p,q,iposv,ivpos,iord)
c     IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
C        -----  A GIVENS HOUSHOLDER MATRIX DIAGONALIZATION  -----
C        -----  ROUTINE SAME AS EIGEN BUT WORKS WITH A      -----
C        -----  LINEAR ARRAY.                               -----
      integer n,ndim,nr,ik
      integer*4 n1, n2
      integer IA(1)
      integer IPOSV(*),IVPOS(*),IORD(*)
      double precision A(1),VEC(NDIM,1),EIG(1)
      double precision W(*),GAMMA(*),BETA(*),BETASQ(*)
      double precision P(*),Q(*)
      double precision zero,pt5,one,two,rhosq
      double precision b,s,sum,sgn,temp,d,sqrts,g,a2,cosa,COSAP,DIF
      double precision DIA,QJ,PPBR,PP,PPBS,SINA2,SINA,R1,R2,R12,SHIFT
      double precision U,WTAW,WJ
      integer i1,i,il,ij,ii,j,itemp,jj,k,lv,l,m,nrr,np,npas,nv,nt
      DATA ZERO,PT5,ONE,TWO /0.0D+00,0.5D+00,1.0D+00,2.0D+00/
      DATA RHOSQ /1.0D-22/
C
      IF(N.EQ.0) GO TO 560
      N1=N-1
      N2=N-2
      GAMMA(1)=A(1)
      IF(N2) 280,270,120
  120 DO 260 NR=1,N2
      IK=IA(NR+1)+NR
      B=A(IK)
      S=ZERO
      DO 130 I=NR,N2
      IJ=IA(I+2)+NR
  130 S=S+A(IJ)**2
C
C        -----  PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION ----
C
      A(IK)=ZERO
      IF(S.LE.ZERO) GO TO 250
      S=S+B*B
      SGN=+ONE
      IF(B.GE.ZERO) GO TO 160
      SGN=-ONE
  160 SQRTS= SQRT(S)
      D=SGN/(SQRTS+SQRTS)
      TEMP= SQRT(PT5+B*D)
      W(NR)=TEMP
      A(IK)=TEMP
      D=D/TEMP
      B=-SGN*SQRTS
C
C        -----  -D- IS FACTOR OF PROPORTIONALITY. NOW       -----
C        -----  COMPUTE AND SAVE -W- VECTOR. EXTRA SINGLY   -----
C        -----  SUBSCRIPTED -W- VECTOR FOR SPEED.           -----
C
      DO 170 I=NR,N2
      IJ=IA(I+2)+NR
      TEMP=D*A(IJ)
      W(I+1)=TEMP
  170 A(IJ)=TEMP
C
C        -----  PREMULTIPLY VECTOR -W- BY MATRIX -A- TO     -----
C        -----  OBTAIN -P- VECTOR. SIMULTANEOUSLY ACCUMULATE ----
C        -----  DOT PRODUCT -WP- -- SCALR -K-.              -----
C
      WTAW=ZERO
      DO 220 I=NR,N1
      SUM=ZERO
      II=IA(I+1)
      DO 180 J=NR,I
      IJ=II+J+1
  180 SUM=SUM+A(IJ)*W(J)
      I1=I+1
      IF(N1.LT.I1) GO TO 210
      DO 200 J=I1,N1
      IJ=IA(J+1)+I+1
  200 SUM=SUM+A(IJ)*W(J)
  210 P(I)=SUM
  220 WTAW=WTAW+SUM*W(I)
      DO 230 I=NR,N1
  230 Q(I)=P(I)-WTAW*W(I)
C
C        -----  NOW FORM -PAP- MATRIX, REQUIRED PART        -----
C
      DO 240 J=NR,N1
      QJ=Q(J)
      WJ=W(J)
      JJ=J+1
      DO 240 I=J,N1
      IJ=IA(I+1)+JJ
  240 A(IJ)=A(IJ)-TWO*(W(I)*QJ+WJ*Q(I))
  250 BETA(NR)=B
      BETASQ(NR)=B*B
      IL=IK+1
  260 GAMMA(NR+1)=A(IL)
  270 IJ=IA(N)+N-1
      B=A(IJ)
      BETA(N-1)=B
      BETASQ(N-1)=B*B
      IJ=IJ+1
      GAMMA(N)=A(IJ)
  280 BETASQ(N)=ZERO
C
C        -----  ADJOIN AN IDENTYTY MATRIX TO BE POST-       -----
C        -----  MULTIPLIED BY ROTATIONS                     -----
C
      DO 300 I=1,N
      DO 299 J=1,N
  299 VEC(I,J)=ZERO
  300 VEC(I,I)=ONE
      M=N
      SUM=ZERO
      NPAS=1
      GO TO 400
  310 SUM=SUM+SHIFT
      COSA=ONE
      G=GAMMA(1)-SHIFT
      PP=G
      PPBS=PP*PP+BETASQ(1)
      PPBR= SQRT(PPBS)
      DO 370 J=1,M
      COSAP=COSA
      IF(PPBS.NE.ZERO) GO TO 320
      SINA=ZERO
      SINA2=ZERO
      COSA=ONE
      GO TO 350
  320 SINA=BETA(J)/PPBR
      SINA2=BETASQ(J)/PPBS
      COSA=PP/PPBR
C
C        -----  POSTMULTIPLY IDENTITY BY -P- TRANSPOSE .    -----
C
      NT=J+NPAS
      IF(NT.LT.N) GO TO 330
      NT=N
  330 CONTINUE
      DO 340 I=1,NT
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)
  340 VEC(I,J)=TEMP
  350 DIA=GAMMA(J+1)-SHIFT
      U=SINA2*(G+DIA)
      GAMMA(J)=G+U
      G=DIA-U
      PP=DIA*COSA-SINA*COSAP*BETA(J)
      IF(J.NE.M) GO TO 360
      BETA(J)=SINA*PP
      BETASQ(J)=SINA2*PP*PP
      GO TO 380
  360 PPBS=PP*PP+BETASQ(J+1)
      PPBR= SQRT(PPBS)
      BETA(J)=SINA*PPBR
  370 BETASQ(J)=SINA2*PPBS
  380 GAMMA(M+1)=G
C
C        -----  TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT ----
C
      NPAS=NPAS+1
      IF(BETASQ(M).GT.RHOSQ) GO TO 410
  390 EIG(M+1)=GAMMA(M+1)+SUM
  400 BETA(M)=ZERO
      BETASQ(M)=ZERO
      M=M-1
      IF(M.EQ.0) GO TO 430
      IF(BETASQ(M).LE.RHOSQ) GO TO 390
C
C        -----  TAKE ROOT OF CORMER 2 BY 2 NEAREST TO       -----
C        -----  LOWER DIAGONAL IN VALUE AS ESTIMATE OF      -----
C        -----  EIGENVALUE TO USE FOR SHIFT                 -----
C
  410 A2=GAMMA(M+1)
      R2=PT5*A2
      R1=PT5*GAMMA(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP= SQRT(DIF*DIF+BETASQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF=  ABS(A2-R1)-  ABS(A2-R2)
      IF(DIF.ge.ZERO) then
         SHIFT=R2
      else
         SHIFT=R1
      endif
      GO TO 310
  430 EIG(1)=GAMMA(1)+SUM
      DO 440 J=1,N
      IPOSV(J)=J
      IVPOS(J)=J
  440  IORD(J)=J
      M=N
      GO TO 470
  450 DO 460 J=1,M
      IF(EIG(J).LE.EIG(J+1)) GO TO 460
      TEMP=EIG(J)
      EIG(J)=EIG(J+1)
      EIG(J+1)=TEMP
      ITEMP=IORD(J)
      IORD(J)=IORD(J+1)
      IORD(J+1)=ITEMP
  460 CONTINUE
  470 M=M-1
      IF(M.NE.0) GO TO 450
      IF(N1.EQ.0) GO TO 500
      DO 490 L=1,N1
      NV=IORD(L)
      NP=IPOSV(NV)
      IF(NP.EQ.L) GO TO 490
      LV=IVPOS(L)
      IVPOS(NP)=LV
      IPOSV(LV)=NP
      DO 480 I=1,N
      TEMP=VEC(I,L)
      VEC(I,L)=VEC(I,NP)
  480 VEC(I,NP)=TEMP
  490 CONTINUE
C
C        -----  BACK TRANSFORM THE VECTORS OF THE TRIPLE    -----
C        -----  DIAGONAL MATRIX.                            -----
C
  500 DO 550 NRR=1,N
      K=N1
  510 K=K-1
      IF(K.LE.0) GO TO 550
      SUM=ZERO
      DO 520 I=K,N1
      IJ=IA(I+1)+K
  520 SUM=SUM+VEC(I+1,NRR)*A(IJ)
      SUM=SUM+SUM
      DO 530 I=K,N1
      IJ=IA(I+1)+K
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(IJ)
      GO TO 510
  550 CONTINUE
  560 CONTINUE
      RETURN
      END
      subroutine geom_check_all_connected(nat,nbnds,ijmods,
     $     c, numijxtra, ijxtra, mxbond, cvrfac, cvr, 
     $     okay, oprint)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
      integer nat, nbnds, numijxtra, mxbond
      integer ijmods(2,nbnds), ijxtra(2,mxbond)
      double precision c(3,nat), cvrfac, cvr(*)
      logical okay              ! return true if all connected
      logical oprint
c     
c     Examine the connectivity generated by autoz and ensure
c     that it connects everything.  If it does not then
c     whine to the user, attempt to add extra bonds to 
c     connect the clusters and return not okay.
c
c     Returns
c     
c     okay=.true. and numijxtra as on input ... ijmods is good
c
c     okay=.false. and numijxtra non-zero ... ijmods is not good
c     .           rerun connection algorithm with extra bonds
c
c     okay=.false. and numijxtra=0 ... ijmods is not good and
c     .           too many extra bonds ... modify cvr_scaling
c     .           or some other parameter.
c     
      integer i, j, k, k_m, l_m, k_clus, l_clus, ibig, test, nclus
      integer ninclus, iclus, jclus
      integer ind, imin, jmin, nnew
      double precision toangs, r, dist, rmin, close
      parameter (toangs=0.52917724924d+00)
      ind(i,j) = k_m + (i-1) + (j-1)*nat
      dist(i,j)=sqrt((c(1,i)-c(1,j))**2+
     1     (c(2,i)-c(2,j))**2+
     2     (c(3,i)-c(3,j))**2) * toangs
c     
      if (nbnds .eq. 0) then
         okay = nat .eq. 1
         return
      endif
c     
      if (.not. ma_push_get(mt_int, nat*nat,'check all0',l_m, k_m))
     $     call errquit('geom_check_all_connected: ma?',nat*nat,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, nat,'check all1',l_clus, k_clus))
     $     call errquit('geom_check_all_connected: ma?',nat,
     &       MA_ERR)
c     
      ibig = 99999999
      call ifill(nat*nat, ibig, int_mb(k_m), 1)
c     
c     Run Floyd's algorithm to find minimal path from each
c     node to any other
c     
      do k = 1, nbnds
         i = ijmods(1,k)
         j = ijmods(2,k)
         int_mb(ind(i,j)) = 1
         int_mb(ind(j,i)) = 1
      enddo
      do i = 1, nat
         int_mb(ind(i,i)) = 0
      enddo
c     
      do k = 1, nat
         do i = 1, nat
            do j = 1, nat
               test = int_mb(ind(i,k)) + int_mb(ind(k,j))
               if (test .lt. int_mb(ind(i,j))) 
     $              int_mb(ind(i,j)) = test
            enddo
         enddo
      enddo
c     
      okay = .true.
      do j = 1, nat
         do i = 1, nat
            if (int_mb(ind(i,j)) .eq. ibig) then
               int_mb(ind(i,j)) = -1
               okay = .false.
            endif
         enddo
      enddo
c     
      if (util_print('connection matrix',print_debug)) then
         write(LuOut,*) ' No. of bonds between each pair of atoms'
         do j = 1, nat
            write(LuOut,1) (int_mb(ind(i,j)),i=1,nat)
 1          format(1x,150i3)
         enddo
      endif
c     
      if (okay) goto 10000      ! Hurray!
c     
c     First identify the connected fragments
c     
      if (oprint) write(LuOut,*) 
     $     ' autoz: The atoms group into disjoint clusters'
      call ifill(nat, 0, int_mb(k_clus), 1)
      nclus = 0
      do j = 1, nat
         if (int_mb(k_clus+j-1) .eq. 0) then
            nclus = nclus + 1
            ninclus = 0
            if (oprint) write(LuOut,2) nclus
 2          format(' cluster ',i3,':',$)
            do i = j, nat
               if (int_mb(k_clus+i-1) .eq. 0) then
                  if (int_mb(ind(i,j)).ne.-1) then
                     ninclus = ninclus + 1
                     int_mb(k_clus+i-1) = nclus
                     if (oprint) write(LuOut,3) i
 3                   format(i5,$)
                     if (ninclus .gt. 10) then
                        if (oprint) write(LuOut,*)
                        if (oprint) write(LuOut,22)
 22                     format('             ',$)
                        ninclus = 0
                     endif
                  endif
               endif
            enddo
            if (oprint) write(LuOut,*)
         endif
      enddo
c
c     Now have identified the clusters.  Loop thru pairs
c     of clusters and find the closest atoms.  If this
c     distance is less than cvrfac*(r1+r2) then attempt to join
c     the clusters via this pair of atoms.
c     
c     If we have too many extra bonds then return with
c     no extra bonds at all and the outer code will
c     increase cvr_scaling to attempt to compensate.
c     
      nnew = 0
      do iclus = 1, nclus
         do jclus = iclus+1, nclus
            rmin = 1d300
            do i = 1, nat
               if (int_mb(k_clus+i-1).eq.iclus) then
                  do j = 1, nat
                     if (int_mb(k_clus+j-1).eq.jclus) then
                        r = dist(i,j)
                        if (r .lt. rmin) then
                           imin = i
                           jmin = j
                           rmin = r
                        endif
                     endif
                  enddo
               endif
            enddo
            close = cvrfac*(cvr(imin)+cvr(jmin))
            if (rmin .le. close) then
               if (numijxtra .eq. mxbond) then
                  numijxtra = 0
                  okay = .false.
                  goto 10000
               endif
               if (oprint) write(LuOut,44) iclus, jclus, 
     &                     imin, jmin, rmin
 44            format(' Connecting clusters',2i4,' via atoms',
     $              2i5,' r =',f5.2)
               numijxtra = numijxtra + 1
               ijxtra(1,numijxtra) = imin
               ijxtra(2,numijxtra) = jmin
               nnew = nnew + 1
            else
               if (oprint .and. rmin.lt.(close+2d0))
     $              write(LuOut,45) iclus, jclus, imin, jmin, rmin
 45            format(' Clusters',2i4,' closest atoms are',
     $              2i5,' r =',f5.2, ' ... close but too far?')
            endif
         enddo
      enddo
c
      if (nnew .eq. 0) then
c     Something is wrong but we did not make new bonds
         okay = .false.
         numijxtra = 0
      else
         okay = .false.
      endif
c     
10000 if (.not. ma_chop_stack(l_m)) 
     $     call errquit('geom_check_all_connected: ma pop?',nat,
     &       MA_ERR)
c     
      end
      logical function geom_hnd_check_data(name,n)
      implicit none
#include "global.fh"
#include "geomhnd.fh"
#include "mafdecls.fh"
      character*(*) name
      integer n
c
      character*255 filename
c     
      integer m
      logical status
c     
c     Verify that data of the given name and dimension exists.
c     It is assumed that the file is permanent.
c     
      call util_file_name(name,.false.,.false.,filename)
c     
      if (ga_nodeid() .eq. 0) then
         open(32,file=filename,form='unformatted',status='old',err=10)
         read(32,err=11) m
         if (m.ne.n) goto 11
         close(32)
         status = .true.
         goto 20
c     
 11      close(32)
 10      status = .false.
      endif
c     
 20   if (parallel) call ga_brdcst(323, status, 
     ,     ma_sizeof(mt_log,1,mt_byte), 0)
      geom_hnd_check_data = status
c     
      end
      subroutine geom_hnd_put_data(name, a, n)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
c
      character*(*) name
      integer n
      double precision a(n)
c
c     Put the data of given name and dimension from a()
c     into the appropriate file.
c
      character*255 filename
      logical odebug, oscratch
c
      oscratch = (name .ne. 'drv.hess' .and. name.ne.'zsym')
c
      odebug = util_print('geom_hnd',print_never)
      if (ga_nodeid() .eq. 0) then
         call util_file_name(name,oscratch,.false.,filename)
         open(32, file=filename, form='unformatted', 
     $        status='unknown', err=100)
         write(32) n
         call swrite(32, a, n)
         close(32)
         if (odebug) then
            write(LuOut,*) ' WROTE ', name, n, 
     $           filename(1:inp_strlen(filename))
            call doutput(a,1,n,1,1,n,1,1)
         endif
      endif
c
      return
c
 100  write(LuOut,*) ' name = ', name
      call errquit('geom_hnd_put_data: open failed',0, GEOM_ERR)
c
      end
      subroutine geom_hnd_get_data(name, a, n)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
#include "geomhnd.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
      character*(*) name
      integer n
      double precision a(n)
c
c     Get the data of given name and dimension into a()
c     from the appropriate file.
c
      character*255 filename
      logical odebug, oscratch
      integer m
c
      oscratch = (name .ne. 'drv.hess' .and. name.ne.'zsym')
      odebug = util_print('geom_hnd',print_never)
      if (ga_nodeid() .eq. 0) then
         call util_file_name(name,oscratch,.false.,filename)
         open(32, file=filename, form='unformatted', 
     $        status='old', err=100)
         if (odebug) write(LuOut,*) ' READ  ', name, n
         read(32,err=101) m
         if (n .ne. m) goto 102
         call sread(32, a, n)
         close(32)
         if (odebug) call doutput(a,1,n,1,1,n,1,1)
      endif
      if (parallel) call ga_brdcst(32, a, 
     $  ma_sizeof(MT_DBL,n,MT_BYTE), 0)
c
      return
c
 100  write(LuOut,*) ' name = ', name
      call errquit('geom_hnd_get_data: open failed',0, DISK_ERR)
 101  write(LuOut,*) ' name = ', name
      call errquit('geom_hnd_get_data: reading m failed',0, DISK_ERR)
 102  write(LuOut,*) ' name = ', name, ' n = ', n, ' m = ', m
      call errquit('geom_hnd_get_data: dimension mismatch',0,
     &       GEOM_ERR)
c
      end
      subroutine geom_hnd_parallel(mode)
      implicit none
#include "geomhnd.fh"
      logical mode
      parallel = mode
      end
      subroutine geom_sort4(iii,jjj,kkk,lll,i1,i2,i3,i4)
      implicit none
      integer iii,jjj,kkk,lll,i1,i2,i3,i4
c
c     Return in i[1-4] iii...lll sorted into ascending order
c
      integer i, j, k(4), itmp
c
      k(1) = iii
      k(2) = jjj
      k(3) = kkk
      k(4) = lll
      do i = 1, 4
         do j = i, 4
            if (k(i).gt.k(j)) then
               itmp = k(i)
               k(i) = k(j)
               k(j) = itmp
            endif
         enddo
      enddo
      i1 = k(1)
      i2 = k(2)
      i3 = k(3)
      i4 = k(4)
c
      end
      logical function geom_add_tor(iat, jat, kat, lat,
     $     ntor, ijklto, mxtors, ntor_for_ij)
      implicit none
c
      integer iat, jat, kat, lat
      integer ntor
      integer ijklto(4,*)
      integer mxtors, ntor_for_ij
c
c     Add the torsion to the list checking for duplicates.
c     Return false if there is no room left.
c
      integer iii, jjj, kkk, lll
      integer itmp, itor
c
      geom_add_tor = .true.
c     
c     Put indices into a standard order
c     
      iii = iat
      jjj = jat
      kkk = kat
      lll = lat
c
      if (iii .gt. lll) then
         itmp = iii
         iii  = lll
         lll  = itmp
         itmp = jjj
         jjj  = kkk
         kkk  = itmp
      endif
c     
c     Figure out if this is a duplicate.  Must remove duplicates
c     now otherwise we will not know if something is missing.
c     
c     But this should count in the overall number for this bond.  Really?
c     
      ntor_for_ij = ntor_for_ij + 1
      do itor = 1, ntor
         if ( iii.eq.ijklto(1,itor) .and.
     $        jjj.eq.ijklto(2,itor) .and.
     $        kkk.eq.ijklto(3,itor) .and.
     $        lll.eq.ijklto(4,itor) ) goto 666
      enddo
      if(ntor.eq.mxtors) then
         geom_add_tor = .false.
      else
         ntor=ntor+1
         ijklto(1,ntor)=iii
         ijklto(2,ntor)=jjj
         ijklto(3,ntor)=kkk
         ijklto(4,ntor)=lll
      endif
c
 666  return
c
      end
      subroutine geom_make_constraint_matrix(geom, nzvar, cmat,
     $     nconstraint)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "stdio.fh"
      integer geom, nzvar       ! input
      double precision cmat(nzvar, nzvar) ! scratch
      integer nconstraint       ! output
c
c     Make a matrix that projects a step or gradient onto the 
c     constrained variables.  Store in file 'cmat'.
c
c     Return in nconstraint the number of constraints so that
c     the number of indep. variables can be reduced.
c
c     E.g.,
c
c     c(i,i) = 1d0 ... i is an indep. variable
c     c(i,i) = 0d0 ... i is frozen.
c     c(i,i) = c(i,j) = c(j,i) = c(j,j) = 1/2 ... i=j constraint
c
c     C*g = sum(j) c(i,j)*g(j) = constrained g.
c
      integer i, ii, j, k, num
      double precision sum
      logical ofrozen
      logical geom_zmt_get_varinfo, geom_zmt_get_nzfrz, 
     $     geom_zmt_get_izfrz
      character*8 test
c
      integer mxcoor
      parameter (mxcoor = nw_max_coor)
      character*8 zvarname(mxcoor)
      integer nzfrz, izfrz(mxcoor)
      double precision zvarsign(mxcoor)
c
c     Start with identity
c
      nconstraint = 0
      call dfill(nzvar*nzvar, 0d0, cmat, 1)
      call dfill(nzvar, 1d0, cmat, nzvar+1)
c
c     Zero diagonals for constants
c
      if(.not.geom_zmt_get_nzfrz(geom,nzfrz))
     $     call errquit('geom_const: geom_zmt_get_nzfrz failed',0,
     &       GEOM_ERR)
      if(.not.geom_zmt_get_izfrz(geom,izfrz))
     $     call errquit('geom_const: geom_zmt_get_izfrz  failed',0,
     &       GEOM_ERR)
      do ii = 1, nzfrz
         i = izfrz(ii)
         nconstraint = nconstraint + 1
         cmat(i,i) = 0d0
      enddo
c
c     Apply constraints
c
      if (.not. geom_zmt_get_varinfo(geom,zvarname,zvarsign))
     $     call errquit('geom_const: varinfo?',0, GEOM_ERR)
      do i = 1, nzvar
         ofrozen = .false.
         do ii = 1, nzfrz
            if (i .eq. izfrz(ii)) ofrozen = .true.
         enddo
         if (zvarname(i).ne.' ' .and. (.not. ofrozen)) then
            num = 0
            test = zvarname(i)
            do j = 1, nzvar
               if (zvarname(j).eq.test) then 
                  num = num + 1
               endif
            enddo
            nconstraint = nconstraint + num - 1
            if (num .gt. 1) then
               sum = 1d0 / dble(num)
               do j = 1, nzvar
                  if (zvarname(j).eq.test) then 
                     do k = 1, nzvar
                        if (zvarname(k).eq.test) 
     $                       cmat(k,j) = zvarsign(k)*zvarsign(j)*sum
                     enddo
                  endif
               enddo
               do j = 1, nzvar
                  if (zvarname(j).eq.test) zvarname(j) = ' '
               enddo
            endif
         endif
      enddo
c
*      write(LuOut,*) ' Constraint matrix ', nzvar, nconstraint
*      call output(cmat, 1, nzvar, 1, nzvar, nzvar, nzvar, 1)
c
      call geom_hnd_put_data('c',cmat, nzvar**2)
c
      end
      logical function bonds_span3d(iat,x,nibond,iibond,ijbond)
      implicit none
      integer iat, nibond(*), iibond(*), ijbond(2,*)
      double precision x(3,*)
c
c     Return true if the bonds connecting to iat have strong
c     components in all 3D.  (use to determine if need to 
c     add out-of-plane bends for a bond.
c
c     A tetrahedral bond is 50 degrees out of plane
c     Ammonia is 34 degrees.  
c
      integer jat, kat, lat, ijbnd, ikbnd, ilbnd, i
      double precision out_of_plane_angle
      double precision a(3), b(3), c(3), angle
c
      bonds_span3d = .false.
      if (nibond(iat).lt.3) return
c      
      do ijbnd=iibond(iat)+1,iibond(iat)+nibond(iat)
         jat=ijbond(2,ijbnd)
         do ikbnd=ijbnd+1,iibond(iat)+nibond(iat)
            kat=ijbond(2,ikbnd)
            do ilbnd=ikbnd+1,iibond(iat)+nibond(iat)
               lat=ijbond(2,ilbnd)
               do i = 1, 3
                  a(i) = x(i,iat) - x(i,jat)
                  b(i) = x(i,iat) - x(i,kat)
                  c(i) = x(i,iat) - x(i,lat)
               enddo
               angle = out_of_plane_angle(a,b,c)
               if (angle .gt. 15d0) then
                  bonds_span3d = .true.
                  return
               endif
            enddo
         enddo
      enddo
c
      end
      double precision function out_of_plane_angle(a,b,c)
      implicit none
      double precision a(3),b(3),c(3)
c
c     a, b, & c are three position vectors from the origin.
c     Return the angle in degrees that c makes with the plane 
c     containing a and b.  The sign is determined in a left-handed
c     sense ... so that if a=x, b=y then for c=+z angle=+90.
c
      double precision d(3), abc, theta
      double precision ddot
      external ddot
c
      call cross_product(a,b,d)
      abc = ddot(3,c,1,d,1)/sqrt(
     $     ddot(3,a,1,a,1)*ddot(3,b,1,b,1)*ddot(3,c,1,c,1))
      theta = 90d0 - acos(min(1d0,abs(abc)))*45d0/atan(1d0)
      if (abc .lt. 0d0) theta = -theta
      out_of_plane_angle = theta
c
      end
