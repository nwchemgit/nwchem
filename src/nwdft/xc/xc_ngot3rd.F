c $Id$
c
c     return true for functionals without 3rds
c
      logical function xc_ngot3rd()
      implicit none
c
#include "cdft.fh"
c
      xc_ngot3rd = (
c
     .     xfac(4).ne.0d0.or.   ! HCTH
     .     xfac(5).ne.0d0.or.   ! HCTH (b970)
     .     xfac(6).ne.0d0.or.   ! HCTH (b971)
     .     xfac(8).ne.0d0.or.   ! Gill 96
     .     xfac(9).ne.0d0.or.   ! PW91
     .     xfac(10).ne.0d0.or.  ! HCTH (h120)
     .     xfac(11).ne.0d0.or.  ! HCTH (h147)
     .     xfac(12).ne.0d0.or.  ! HCTH (b980)
     .     xfac(13).ne.0d0.or.  ! HCTH (b97g)
     .     xfac(14).ne.0d0.or.  ! HCTH (h407)
     .     xfac(15).ne.0d0.or.  ! HCTH (hp14)
     .     xfac(16).ne.0d0.or.  ! OPTX
     .     xfac(17).ne.0d0.or.  ! HCTH (b972)
     .     xfac(18).ne.0d0.or.  ! PKZB99-EXC
     .     xfac(19).ne.0d0.or.  ! mPW91
     .     xfac(20).ne.0d0.or.  ! HCTH (407p) 
     .     xfac(21).ne.0d0.or.  ! TPSS03-EXC
     .     xfac(22).ne.0d0.or.  ! HCTH (b973)
     .     xfac(25).ne.0d0.or.  ! ft97 
     .     xfac(26).ne.0d0.or.  ! PW6B95
     .     xfac(27).ne.0d0.or.  ! PW6BK
     .     xfac(28).ne.0d0.or.  ! M05 
     .     xfac(29).ne.0d0.or.  ! M05-2X
     .     xfac(32).ne.0d0.or.  ! dlDF
     .     xfac(33).ne.0d0.or.  ! VSXC
     .     xfac(34).ne.0d0.or.  ! M06-L
     .     xfac(35).ne.0d0.or.  ! M06-HF
     .     xfac(36).ne.0d0.or.  ! M06
     .     xfac(37).ne.0d0.or.  ! M06-2X
     .     xfac(39).ne.0d0.or.  ! HCTH (b97d)
     .     xfac(46).ne.0d0.or.  ! SSB-D
     .     xfac(47).ne.0d0.or.  ! kt1
     .     xfac(48).ne.0d0.or.  ! M08-HX 
     .     xfac(49).ne.0d0.or.  ! M08-SO
     .     xfac(50).ne.0d0.or.  ! M11
     .     xfac(51).ne.0d0.or.  ! M11-L
     .     xfac(52).ne.0d0.or.  ! SOGGA
     .     xfac(53).ne.0d0.or.  ! SOGGA11
     .     xfac(54).ne.0d0.or.  ! SOGGA11-X
c
     .     cfac(5).ne.0d0.or.   ! PW91
     .     cfac(13).ne.0d0.or.  ! HCTH
     .     cfac(14).ne.0d0.or.  ! HCTH (b970)
     .     cfac(15).ne.0d0.or.  ! HCTH (b971)
     .     cfac(16).ne.0d0.or.  ! HCTH (h120)
     .     cfac(17).ne.0d0.or.  ! HCTH (h147)
     .     cfac(18).ne.0d0.or.  ! HCTH (b980)
     .     cfac(19).ne.0d0.or.  ! HCTH (b97g)
     .     cfac(20).ne.0d0.or.  ! HCTH (h407)
     .     cfac(21).ne.0d0.or.  ! HCTH (hp14)
     .     cfac(22).ne.0d0.or.  ! OPTC
     .     cfac(23).ne.0d0.or.  ! HCTH (b972)
     .     cfac(24).ne.0d0.or.  ! ft97 
     .     cfac(25).ne.0d0.or.  ! PKZB99-COR
     .     cfac(26).ne.0d0.or.  ! HCTH (407p)
     .     cfac(27).ne.0d0.or.  ! TPSS03-COR
     .     cfac(28).ne.0d0.or.  ! HCTH (b973)
     .     cfac(29).ne.0d0.or.  ! cvs98
     .     cfac(30).ne.0d0.or.  ! M06-L
     .     cfac(31).ne.0d0.or.  ! Bc95-COR
     .     cfac(32).ne.0d0.or.  ! PW6B95
     .     cfac(33).ne.0d0.or.  ! PWB6K
     .     cfac(34).ne.0d0.or.  ! M05
     .     cfac(35).ne.0d0.or.  ! M05-2X
     .     cfac(36).ne.0d0.or.  ! Be88 
     .     cfac(37).ne.0d0.or.  ! OP 
     .     cfac(38).ne.0d0.or.  ! M06-HF
     .     cfac(39).ne.0d0.or.  ! M06
     .     cfac(40).ne.0d0.or.  ! M06-2X
     .     cfac(41).ne.0d0.or.  ! HCTH (b97d)
     .     cfac(42).ne.0d0.or.  ! dlDF
     .     cfac(46).ne.0d0.or.  ! SSB-D
     .     cfac(48).ne.0d0.or.  ! M08-HX
     .     cfac(49).ne.0d0.or.  ! M08-SO
     .     cfac(50).ne.0d0.or.  ! M11
     .     cfac(51).ne.0d0.or.  ! M11-L
     .     cfac(52).ne.0d0.or.  ! SOGGA
     .     cfac(53).ne.0d0.or.  ! SOGGA11
     .     cfac(54).ne.0d0      ! SOGGA11-X
     .                )
c
      return
      end
