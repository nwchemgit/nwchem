*fordeck triprt $Revision: 7.7 $
      Subroutine TriPrt(Title,FmtIn,A,N)
************************************************************
*
*   <DOC>
*     <Name>TriPrt</Name>
*     <Syntax>Call TriPrt(Title,FmtIn,A,N)</Syntax>
*     <Arguments>
*       \Argument{Title}{String containing a title}{Character*(*)}{in}
*       \Argument{FmtIn}{String containing a format}{Character*(*)}{in}
*       \Argument{A}{Triangular matrix to be printed}{Array Real*8}{in}
*       \Argument{N}{Dimension of matrix A}{Integer}{in}
*     </Arguments>
*     <Purpose>Print a triangular matrix</Purpose>
*     <Dependencies></Dependencies>
*     <Author></Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     </Description>
*    </DOC>
*
************************************************************
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Print a square matrix stored in packed, lower triangular         *
*     storage mode                                                     *
*                                                                      *
*     calling arguments                                                *
*     Title  : character string containing a title                     *
*              If the string is empty no title will be printed         *
*     A      : triangular matrix of double precision reals             *
*     N      : dimension of matrix A                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Character*(*) Title
      Character*(*) FmtIn
      Dimension A(N*(N+1)/2)
      Integer StrnLn
      Parameter (lPaper=120)
      Character*(lPaper) Line
      Character*20 FMT
*----------------------------------------------------------------------*
      If (N.le.0) Return
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     print the title                                                  *
*----------------------------------------------------------------------*
      lTitle=StrnLn(Title)
      If ( lTitle.gt.0 ) then
         Do 10 i=1,lPaper
             Line(i:i)=' '
10       Continue
         lLeft=1
         Do 20 i=lTitle,1,-1
            If ( Title(i:i).ne.' ' ) lLeft=i
20       Continue
         lLeft=lLeft-1
         Do 25 i=1,lPaper
            If ( i+lLeft.le.lTitle ) Line(i:i)=Title(i+lLeft:i+lLeft)
25       Continue
         Write(6,*)
         Write(6,'(2X,A)') Line
         Do 30 i=1,StrnLn(Line)
            Line(i:i)='-'
30       Continue
         Write(6,'(2X,A)') Line
         Write(6,'(2X,A,I4,A,I4)') 'mat. size = ',N,'x',N
      End If
*----------------------------------------------------------------------*
*     determine the printing format                                    *
*----------------------------------------------------------------------*
      lFmt=StrnLn(FmtIn)
      If ( lFmt.ne.0 ) then
         FMT=FmtIn
      Else
         Amax=A(1)
         Amin=A(1)
         Do 40 i=1,N*(N+1)/2
            Amax=Max(Amax,A(i))
            Amin=Min(Amin,A(i))
40       Continue
         If (Amax.ne.0.0D0) Then
           Pmax=Log10(Abs(Amax))
           iPmax=Int(1d0+Pmax)
           iPmax=Max(1,iPmax)
         Else
           iPmax=1
         End If
         If (Amin.ne.0.0D0) Then
           Pmin=Log10(Abs(Amin))
           iPmin=Int(1d0+Pmin)
           iPmin=Max(1,iPmin)
         Else
           iPmin=1
         End If
         nDigit=14
         nDecim=Min(8,ABS(nDigit-Max(iPmin,iPmax)))
         If ( Amax.lt.0d0 ) iPmax=iPmax+1
         If ( Amin.lt.0d0 ) iPmin=iPmin+1
         lNumbr=Max(iPmin,iPmax)+nDecim+2
         nCols=10
         lLine=nCols*lNumbr
         If ( lLine.gt.lPaper ) then
            If ( lLine.le.lPaper+nCols .and. nDecim.gt.1 ) then
               nDecim=nDecim-1
               lNumbr=Max(iPmin,iPmax)+nDecim
               lItem=Max(lNumbr,lPaper/nCols)
            Else
               nCols=5
               lItem=Max(lNumbr,lPaper/nCols)
            End If
         Else
            lItem=lNumbr
         End If
         Write(FMT,'(A,   I4.4,  A, I4.4,  A, I4.4,   A)')
     &             '(2X,',nCols,'F',lItem,'.',nDecim,')'
      End if
*----------------------------------------------------------------------*
*     print the data                                                   *
*----------------------------------------------------------------------*
      Write(6,*)
      jEnd=0
      Do 60 i=1,N
         jStart=jEnd+1
         jEnd=jEnd+i
         Write(6,FMT)(A(j),j=jStart,jEnd)
60    Continue
*----------------------------------------------------------------------*
*     End procedure                                                    *
*----------------------------------------------------------------------*
      Return
      End
c $Id$
