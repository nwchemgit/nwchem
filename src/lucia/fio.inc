*comdeck fio.inc $Revision: 7.7 Patch(7.7): 326_io $
************************************************************************
*                                                                      *
*     (c) Copyright. All rights reserved                               *
*                                                                      *
*     This code or parts thereof may not be copdied or redistributed,  *
*     without the written permission of the author. The use is re -    *
*     stricted to research purposes only, and the matrial may not be   *
*     included in any commercial product.                              *
*                                                                      *
************************************************************************
*                                                                      *
*     This file defines the common blocks /FIO1/ - /FIO4/              *
*     including all entries required for fast I/O system.              *
*                                                                      *
*     Following the list of entries and their usage:                   *
*     LuName    : Name of the file                                       *
*     isOpen  : open/close flag                                        *
*     FSCB    : file descriptors (C language)                          *
*     Addr    : pointer to the current position                        *
*     MxAddr  : largest address written                                *
*     Count   : counters for the I/O statistics                        *
*               (# of I/O startups, amount of data transferred)        *
*     Trace   : enable/disable debugging output                        *
*     Query   : enable/disable the traceback facilities                *
*     FirstCall : enable/disable multifile partitiong                  *
*     RefBuf  : anchor of the gather/scatter process                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L Serrano-Andres,University of Lund, Sweden, 1996                *
*                                                                      *
************************************************************************
      Integer MaxSplitFile
      Parameter (MaxSplitFile=20)
      Include 'MxFile.inc'
      Integer isOpen,FSCB,Addr,MxAddr,RefBuf,MaxFileSize,MPUnit,MBL
      Common /FIO1/ isOpen(MxFile),
     &              FSCB(MxFile),
     &              Addr(MxFile),
     &              MxAddr(MxFile),
     &              RefBuf(8),
     &              MPUnit(0:MaxSplitFile-1,MxFile),
     &              MBL(MxFile),
     &              MaxFileSize

      Logical     Trace,Query,FirstCall,Multi_File
      Common /FIO2/ Trace,
     &              Query,
     &              FirstCall,
     &              Multi_File(MxFile)

      Real*8      Count
      Common /FIO3/ Count(4,MxFile)

      Character*8 LuName
      Common /FIO4/ LuName(MxFile)
      Integer LuMark(MxFile)
      Common /FIO5/ LuMark
