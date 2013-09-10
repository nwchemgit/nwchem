*fordeck switch.inc $Revision: 7.7 $
*----------------------------------------------------------------------*
* Internal switches for AIX/IO routines.                               *
*----------------------------------------------------------------------*
      Parameter (eEof   = 1024,      eNtOpn = eEof   +1)
      Parameter (eInErr = eNtOpn +1, eTmF   = eInErr +1)
      Parameter (eTlFn  = eTmF   +1, eBlNme = eTlFn  +1)
      Parameter (eNoMsg = eBlNme +1)
