#ifdef SNGLPR
      Subroutine SGeCop (TRANSA, M, N, A, LDA, B, LDB)
#else
      Subroutine DGeCop (TRANSA, M, N, A, LDA, B, LDB)
#endif
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: dgecop.F,v 1.1 1994-06-14 21:54:10 gg502 Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     dgecop, sgecop -- copies one matrix into another
C
C SYNOPSIS
      Character*1 TRANSA
      Integer M, N, LDA, LDB
#ifdef SNGLPR
      Real A(LDA,*), B(LDB,*)
#else
      Double precision A(LDA,*), B(LDB,*)
#endif
C
C ARGUMENTS
C     TRANSA  Consider the source to be A or A' (input)
C     M       Rows of B (input)
C     N       Columns of B (input)
C     A       Source matrix (input)
C     LDA     Leading dimension of A (input)
C     B       Desitnation matrix (output)
C     LDB     Leading dimension of B (input)
C
C DESCRIPTION
C     Copies one general matrix into another
C
C ROUTINES REQUIRED
C     From BLAS -- Support routines
      Logical LSame
      External LSame, XErBLA
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Integer i, j, Info
      Logical NotA
C
C     Test input parameters
C
      NotA = LSame( TransA, 'N')
C
      Info = 0
      If ( .not. (NotA .OR. LSame(TransA, 'T')
     $   .OR. LSame(TransA, 'C')) ) then
         Info = 1
      ElseIf (M .lt. 0) then
         Info = 2
      ElseIf (N .lt. 0) then
         Info = 3
      ElseIf ((LDA .lt. Max(1, M) .AND. NotA)
     $      .OR. ( LDA .lt. Max(1, N) .AND. .NOT. NotA)) then
         Info = 5
      ElseIf (LDB .lt. Max(1, M)) then
         Info = 7
      EndIf
      If (Info .ne. 0) then
#ifdef SNGLPR
         Call XErBLA('SGECOP', Info)
#else
         Call XErBLA('DGECOP', Info)
#endif
         Return
      EndIf
C
C     See if we can get away easily
C
      If (M .eq. 0 .OR. N .eq. 0) Return
C
C     Gotta do the work...
C
      If (NotA) then
         Do 100 j = 1, N
            Do 110 i = 1, M
               B(i,j) = A(i,j)
 110        Continue
 100     Continue
      Else
         Do 200 j = 1, N
            Do 210 i = 1, M
               B(i,j) = A(j,i)
 210        Continue
 200     Continue
      EndIf
C
      Return
      End
