*fordeck strnln.f $Revision: 7.7 $
      Integer Function StrnLn(String)
************************************************************************
*                                                                      *
*   (c) Copyright 1993 by the authors of MOLCAS. All rights reserved   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Determine the position of the last nonblank character in         *
*     the input string.                                                *
*                                                                      *
*     calling arguments:                                               *
*     String : Character string                                        *
*                                                                      *
*     return value:                                                    *
*     StrnLn  : Integer, position of last nonblank character           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P.O. Widmark                                  *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Character*(*) String
      Character*(1) Blank,Null
*----------------------------------------------------------------------*
*     set characters                                                   *
*----------------------------------------------------------------------*
      Blank = Char(32)
      Null  = Char(0)
*----------------------------------------------------------------------*
*     get the length of the line                                       *
*----------------------------------------------------------------------*
      lString=LEN(String)
*----------------------------------------------------------------------*
*     loop over the string                                             *
*----------------------------------------------------------------------*
      StrnLn=0
      Do 100 i=1,lString
        If ( String(i:i).ne.Blank .and. String(i:i).ne.Null ) StrnLn=i
100   Continue
*----------------------------------------------------------------------*
*     Normal exit                                                      *
*----------------------------------------------------------------------*
      Return
      End
c $Id$
