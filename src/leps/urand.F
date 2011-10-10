C      program main
C      integer I
C      double precision DYNRND,seed,R
C      COMMON/DYSEED/SEED,NSEED
C      seed=524.3
C      R=dynrnd(1)
C      do i=1,5000
C         write(*,*) 'dynrnd ', dynrnd(0)
C      enddo 

C      end
C      double precision FUNCTION DYNRND(N,SEED)
      FUNCTION DYNRND(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL  DBUG
      COMMON/DYSEED/SEED,NSEED
      CHARACTER*8 ERRMSG
      DIMENSION V(97)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-DYNRND-'/

C     REAL*8 (A-H,O-Z)
C      integer i,N,nseed
C      double precision dum,dummy,rnd,V,SEED,urand
C      LOGICAL  DBUG
C      CHARACTER*8 ERRMSG
C      DIMENSION V(97)
C
      SAVE V
      SAVE RND
C
      DBUG=.FALSE.
C
      IF(N.EQ.1) THEN
         NSEED=IDINT(SEED)
         DUMMY=URAND(NSEED)
         DO I=1,97
            DUM=URAND(NSEED)
         ENDDO
         DO I=1,97
            V(I)=URAND(NSEED)
         ENDDO
         RND=URAND(NSEED)
      ENDIF
      I=1+INT(DBLE(97)*RND)
      IF(I.LT.1.OR.I.GT.97) THEN
         WRITE(*,9999)
C         CALL HNDERR(3,ERRMSG)
      ENDIF
         RND=V(I)
      DYNRND=RND
        V(I)=URAND(NSEED)
C
      IF(DBUG) THEN
         WRITE(*,9998) DYNRND,NSEED
      ENDIF
C
      RETURN
 9999 FORMAT(' SOMETHING WRONG WITH -DYNRND- . STOP. ')
 9998 FORMAT(' In -Function DYNRND, -DYNRND- = ',F12.8,I12)

      END
      FUNCTION URAND(IY)
      DOUBLE PRECISION URAND
c --------------------------------------------------------------
c
c     double precision function urand(iy)
c
c       Routine from the book Computer Methods for Mathematical
c       Computations, by Forsythe, Malcolm, and Moler.
c       Prentice-Hall, Englewood CLiffs, N.J., 1997
c       Developed by George Forsythe, Mike Malcolm, and Cleve Moler.
c       As copied from http://www.netlib.org/fmm/
c
c
c  urand is a uniform random number generator based  on  theory  and
c  suggestions  given  in  D.E. Knuth (1969),  Vol  2.   The integer  iy
c  should be initialized to an arbitrary integer prior to the first call
c  to urand.  The calling program should  not  alter  the  value  of  iy
c  between  subsequent calls to urand.  Values of urand will be returned
c  in the interval (0,1).
c

      integer  iy
      integer  ia,ic,itwo,m2,m,mic
      double precision  halfm
      real  s
      double precision  datan,dsqrt
c
      data m2/0/,itwo/2/
c
      save
c
      if (m2 .ne. 0) go to 20
c
c  if first entry, compute machine integer word length
c
      m = 1
   10 m2 = m
      m = itwo*m2
      if (m .gt. m2) go to 10
      halfm = m2
c
c  compute multiplier and increment for linear congruential method
c
      ia = 8*idint(halfm*datan(1.d0)/8.d0) + 5
      ic = 2*idint(halfm*(0.5d0-dsqrt(3.d0)/6.d0)) + 1
      mic = (m2 - ic) + m2
c
c  s is the scale factor for converting to floating point
c
      s = 0.5/halfm
c
c  compute next random number
c
   20 iy = iy*ia
c
c  the following statement is for computers which do not allow
c  integer overflow on addition
c
      if (iy .gt. mic) iy = (iy - m2) - m2
c
      iy = iy + ic
c
c  the following statement is for computers where the
c  word length for addition is greater than for multiplication
c
      if (iy/2 .gt. m2) iy = (iy - m2) - m2
c
c  the following statement is for computers where integer
c  overflow affects the sign bit
c
      if (iy .lt. 0) iy = (iy + m2) + m2
      urand = dble(iy)*dble(s)
      return
      END

c $Id$
