*fordeck ctl.inc $Revision: 7.7 $
*----------------------------------------------------------------------*
* Control blocks for AIX I/O routines.                                 *
*----------------------------------------------------------------------*
      Include 'MxFile.inc'
      Parameter (pHndle =  1)
      Parameter (pWhere =  pHndle + 1)
      Parameter (pDesc  =  pWhere + 1)
      Parameter (pStat  =  pDesc  + 1)
      Parameter (pEnd   =  pStat)
      Dimension CtlBlk(pEnd,MxFile)
      Character*80 FCtlBlk(MxFile)
      Common /AixCmn/ CtlBlk,FCtlBlk
