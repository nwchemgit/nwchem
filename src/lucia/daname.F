*fordeck daname $Revision: 7.7 $
      Subroutine DaName2(Lu,String)
************************************************************************
*                                                                      *
*     (c) Copyright. All rights reserved                               *
*                                                                      *
*     This code or parts thereof may not be copdied or redistributed,  *
*     without the written permission of the author. The use is re -    *
*     stricted to research purposes only, and the matrial may not be   *
*     included in any commercial product.                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     purpose:                                                         *
*     Open unit Lu for direct access I/O and link the data strean to   *
*     the logical file name Name.                                      *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number                                    *
*     Name    : character string, input                                *
*               logocal file name                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L Serrano-Andres,University of Lund, Sweden, 1996                *
*                                                                      *
************************************************************************

      Character*(*) String
      Integer Lu
      Logical mf,wa
      mf=.false.
      wa=.false.
      Call Daname_Main(Lu, String, mf, wa)

      Return
      End
c $Id$
