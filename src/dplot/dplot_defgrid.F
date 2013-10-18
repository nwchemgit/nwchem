      subroutine dplot_defgrid(volume,angle,
     ,     spread,step,const_m,const_p,
     ,     cooto,coofrom,ngridxyz)
      implicit none
#include "errquit.fh"
#include "util_params.fh"
      double precision volume
      double precision angle(3),spread(3),step(3)
      integer const_m(3),const_p(3)
      double precision cooto(3),coofrom(3) ! [in]
      integer ngridxyz(3) ! [in]
c
      integer i
      double precision scale,frac_p
c
      Volume = 1.0d0
      Scale  = 1d0/cau2ang
c
      Do i = 1, 3
         Angle(i)  = 90.0d0
         If (nGridXYZ(i).gt.0) Then
            Spread(i) = CooTo(i)-CooFrom(i)
            Step(i)   = Spread(i)/dble( nGridXYZ(i) )
            Volume    = Volume * Step(i)
            If (CooFrom(i).ge.0.0d0) Then
               Const_M(i) = 0
               Const_P(i) = nGridXYZ(i)
            Else If (CooTo(i).le.0.0d0) Then
               Const_M(i) =-nGridXYZ(i)
               Const_P(i) = 0
            Else
               Frac_P     = (CooTo(i)-0.5d0*Step(i))/Spread(i)
               Const_P(i) = Int(Frac_P*nGridXYZ(i)) + 1
               Const_M(i) = Const_P(i) - nGridXYZ(i)
            End If
            Spread(i) = Spread(i)/Scale
         Else If (nGridXYZ(i).eq.0) Then
            Volume     =0.d0
            CooFrom(i) = (CooFrom(i)+CooTo(i))/2.d0
            CooTo(i)   = CooFrom(i)
            Spread(i)  = 0.d0
            Const_M(i) = 0
            Const_P(i) = 0
         Else
            Call ErrQuit('Get_Grid: Negative no_of_spacings', 0,
     &       INPUT_ERR)
         End If
      End Do
      return
      end
c $Id$
