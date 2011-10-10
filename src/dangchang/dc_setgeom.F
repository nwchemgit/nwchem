c        1         2         3         4         5         6
C23456789012345678901234567890123456789012345678901234567890123456789
      SUBROUTINE dc_setgeom(rtdb)

      implicit none

#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "dc_common.fh"
#include "util.fh"

C     General variables

C      DOUBLE PRECISION R ![output] coodinates of atoms
      INTEGER rtdb  ![input] rtdb handle
C      INTEGER nmols  ![input] number of molecules
      INTEGER natom ![output] number of atoms
      INTEGER geom  !geometry handle
      DOUBLE PRECISION scale ![output] bohr->angstrom factor

      integer l_c, k_c, l_q, k_q, l_m, k_m, l_t, k_t

C     Variables specific to Mpoint calculation

      DOUBLE PRECISION directionM !vector from Oxygen to M
      DOUBLE PRECISION midptH !vector gets to the bisector
      DOUBLE PRECISION directionMMagnitude !magn of directionM

      INTEGER I    !indicates a particular water molecule
      INTEGER iii  !indicates component
C      INTEGER MXMOLS
C      INTEGER MXATOMS

C      DIMENSION R(MXMOLS,MXATOM,3)
      DIMENSION midptH(MXMOLS,3)
      DIMENSION directionM(MXMOLS,3)
      DIMENSION directionMMagnitude(MXMOLS)

C     MAKE THE GEOM. OBJECT

      if (.not. geom_create(geom, 'geometry'))
     & call errquit('dangchang_getgeom: failed creating geometry',
     & 0,GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     & call errquit('dangchang_getgeom: no geometry',
     & 0,RTDB_ERR)
      if (.not. geom_ncent(geom,natom))
     & call errquit('dangchang_getgeom: ncent error',
     & 0,GEOM_ERR)

C     CHECKS THAT NATOM IS MULTIPLE OF 3

      if ( MOD (natom,3) .NE. 0 )
     & call errquit('dangchang_getgeom: natom not div by 3',
     & 0,GEOM_ERR)

C     DEFINES NMOLS

      NMOLS = natom / 3

C     ALLOCATE LOCAL MEMORY

      if (.not. ma_push_get(mt_dbl, 3*natom, 'coordinate' ,
     & l_c,k_c))
     & call errquit('dangchang_getgeom: coordinate not allocated',
     & 3*natom,MA_ERR)
      if (.not. ma_push_get(mt_dbl, natom, 'charge', l_q,k_q))
     & call errquit('dangchang_getgeom: charge not allocated',
     & natom,MA_ERR)
      if (.not.ma_push_get(mt_dbl, natom, 'mass', l_m,k_m))
     & call errquit('dangchang_getgeom: mass not allocated',
     & natom,MA_ERR)
      if (.not.ma_push_get(mt_byte, 16*natom, 'tags', l_t,k_t))
     & call errquit('dangchang_getgeom: tags not allocated',
     & natom,MA_ERR)

      if (.not. geom_cart_get(geom, natom, byte_mb(k_t),
     & dbl_mb(k_c), dbl_mb(k_q)))
     & call errquit('dangchang_getgeom: coord error',0,GEOM_ERR)
      if (.not. geom_masses_get(geom, natom, dbl_mb(k_m)))
     & call errquit('dangchang_getgeom:  mass error',0,GEOM_ERR)

 
C     CHECKS THAT ONLY ATOMS O and H PRESENT BY TAG
C     w/ THE PATTERN OHHOHHOHH.....
      if (util_print('geometry',print_medium)) then
      if (.not.geom_print(geom))
     & call errquit('dangcheng_getgeom: print error',0,GEOM_ERR)
      endif

      CALL DANGCHANG_TAGCHECK(natom, byte_mb(k_t))

C      WRITE(6,*)'geom after geom_print before geom_get', geom

      if (.not. geom_get_user_scale(geom,scale))
     & call errquit('dangchang_geom_get_user_scale:
     & scale factor unknown',0,GEOM_ERR)
 

C     PUT INFO INTO DANGCHANG VARIABLES

      CALL DANGCHANG_FILLCOORD(scale,natom,dbl_mb(k_c))

C     In DANGCHANG VARIABLES, Computes the Mpoints for the H2O's

      DO I = 1,NMOLS

        DO iii = 1,3
          midptH(I,iii) = R(I,2,iii) + .5D+00*
     $                  ( R(I,3,iii) - R(I,2,iii) )
          directionM(I,iii) =  midptH(I,iii) - R(I,1,iii)
        END DO

        directionMMagnitude(I) = DSQRT(
     $                           directionM(I,1)**2
     $                         + directionM(I,2)**2
     $                         + directionM(I,3)**2 )

        DO iii = 1,3
          R(I,4,iii) = R(I,1,iii) + 0.215D+00*
     $               directionM(I,iii)/directionMMagnitude(I)

C          WRITE(6,*) 'i,4,iii',i,iii,R(I,4,iii)

        END DO

      END DO ! 1st I loop


C     GET RID OF MEMORY WITH MA_CHOP
C     DEALLOCATE THE FIRST/LOWEST PART OF THE STACK
C     THIS ACTION DEALLOCATES ALL ADDITIONAL
C     MEMORY ADDED TO THE STACK AFTER LOWEST PART

      if (.not. ma_chop_stack(l_c))
     & call errquit('dangchang_getgeom: chop of l_c did not work',
     & 0,MA_ERR)

C     DESTROY GEOM. OBJECT

C      WRITE(6,*) 'GEOM before destory', geom

      if (.not. geom_destroy(geom))
     & call errquit('dangchang_getgeom:  geom_destroy failed',
     & 911,GEOM_ERR)

      END !SUBROUTINE DANGCHANG_SETGEOM

C********************************************************************
      SUBROUTINE DANGCHANG_FILLCOORD(scale,natom,Rinit)
C********************************************************************
      implicit none

#include "dc_common.fh"

      INTEGER i, j, k

      INTEGER natom
      DOUBLE PRECISION Rinit
      DOUBLE PRECISION scale

      DIMENSION Rinit(3,natom) 


      DO i = 1, natom, 3
        DO j = 1, 3
          R(INT((i-1)/3)+1,1,j) = Rinit(j,i)/scale
        END DO
      END DO


      DO i = 2, natom, 3
        DO j = 1, 3
          R(INT((i-1)/3)+1,2,j) = Rinit(j,i)/scale
        END DO
      END DO

      DO i = 3, natom, 3
        DO j = 1, 3
          R(i/3,3,j) = Rinit(j,i)/scale
        END DO
      END DO


c      DO i = 1, natom
c        DO j = 1, 3

c          WRITE(6,*) 'i,j=',i,j
c          WRITE(6,*) 'INT((i-1)/3)+1 =', INT((i-1)/3)+1

c          WRITE(6,*) 'Rinit=', Rinit(j,i)/scale

c        END DO
c      END DO

      END !SUBROUTINE FILLCOORD
C**********************************************************************
      SUBROUTINE DANGCHANG_TAGCHECK(natom, tag)
C**********************************************************************
      implicit none

#include "geom.fh"
#include "errquit.fh"

      INTEGER i
      INTEGER atn
      INTEGER natom
      CHARACTER*16 tag
      CHARACTER*16 element
      CHARACTER*2 symbol

      DIMENSION tag(natom) 

      DO i = 1, natom, 3
C        write(6,*) 'i in the first loop is ',i
        if (.not. geom_tag_to_element(tag(i), symbol, element, atn))
     &   call errquit('DANGCHANG_TAGCHECK:  failed to interpret tags',
     &   0,GEOM_ERR)

C        write(6,*) 'atn is ',atn
        if (atn .NE. 8)
     &  call errquit('DANGCHANG_TAGCHECK:  Os not correct',
     &  0,GEOM_ERR)
      END DO

      DO i = 2, natom, 3
C        write(6,*) 'i in the second loop is ',i
        if (.not. geom_tag_to_element(tag(i), symbol, element, atn))
     &   call errquit('DANGCHANG_TAGCHECK:  failed to interpret tags',
     &   0,GEOM_ERR)

C        write(6,*) 'atn is ',atn
        if (atn .NE. 1)
     &  call errquit('DANGCHANG_TAGCHECK:  Hs not correct',
     &  0,GEOM_ERR)
      END DO

      DO i = 3, natom, 3
C        write(6,*) 'i in the third loop is ',i
        if (.not. geom_tag_to_element(tag(i), symbol, element, atn))
     &   call errquit('DANGCHANG_TAGCHECK:  failed to interpret tags',
     &   0,GEOM_ERR)

C        write(6,*) 'atn is ',atn
        if (atn .NE. 1)
     &  call errquit('DANGCHANG_TAGCHECK:  Hs not correct',
     &  0,GEOM_ERR)
      END DO

      END !SUBROUTINE TAGCHECK
c $Id$
