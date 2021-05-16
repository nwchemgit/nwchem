      subroutine util_poltensor(alpha,alpha_ev,alpha_iso,alpha_ani,
     S     string_in)
      implicit none
#include "global.fh"
#include "stdio.fh"      
#include "util_params.fh"      
      double precision alpha(3,3)           ! polarizability [in]
      double precision alpha_ev(3)          ! polarizability evec  [out]
      double precision alpha_iso            ! polarizability isotr [out]
      double precision alpha_ani ! polarizability aniso [out]
      character*(*) string_in
c
      double precision alphad(3,3)
      double precision au3_2_angs3
c     diagonalize
      call dcopy(9,alpha,1,alphad,1)
      call hnd_diag(alphad,alpha_ev,3,.false.,.false.) 
c
      alpha_iso = (alpha_ev(1)+alpha_ev(2)+alpha_ev(3))/3.0d0 ! isotropic part
      alpha_ani = (alpha_ev(1)-alpha_ev(2))*(alpha_ev(1)-alpha_ev(2))
     1  + (alpha_ev(1)-alpha_ev(3))*(alpha_ev(1)-alpha_ev(3))
     2  + (alpha_ev(2)-alpha_ev(3))*(alpha_ev(2)-alpha_ev(3))
      alpha_ani = dsqrt(alpha_ani/2.d0) ! anisotropic part
      if (ga_nodeid().eq.0) write(luout,1) 
     1     string_in//" polarizability / au ",
     2  alpha(1,1),alpha(1,2),alpha(1,3),
     3  alpha(2,1),alpha(2,2),alpha(2,3),
     4  alpha(3,1),alpha(3,2),alpha(3,3),
     5  alpha_ev(1),alpha_ev(2),alpha_ev(3),
     6     alpha_iso,alpha_ani
      au3_2_angs3=cau2ang*cau2ang*cau2ang
      call dscal(9,au3_2_angs3,alpha,1)
      call dscal(3,au3_2_angs3,alpha_ev,1)
      alpha_iso = (alpha_ev(1)+alpha_ev(2)+alpha_ev(3))/3.0d0 ! isotropic part
      alpha_ani = (alpha_ev(1)-alpha_ev(2))*(alpha_ev(1)-alpha_ev(2))
     1  + (alpha_ev(1)-alpha_ev(3))*(alpha_ev(1)-alpha_ev(3))
     2  + (alpha_ev(2)-alpha_ev(3))*(alpha_ev(2)-alpha_ev(3))
      alpha_ani = dsqrt(alpha_ani/2.d0) ! anisotropic part
      if (ga_nodeid().eq.0) write(luout,1) 
     1     string_in//" polarizability / angs^3 ",
     2  alpha(1,1),alpha(1,2),alpha(1,3),
     3  alpha(2,1),alpha(2,2),alpha(2,3),
     4  alpha(3,1),alpha(3,2),alpha(3,3),
     5  alpha_ev(1),alpha_ev(2),alpha_ev(3),
     6     alpha_iso,alpha_ani
c     back to au units
      call dscal(9,1d0/au3_2_angs3,alpha,1)
      call dscal(3,1d0/au3_2_angs3,alpha_ev,1)
 1    format(/,1x, a,/
     2  1x,'             X              Y              Z',/
     3  1x,'-----------------------------------------------',/
     4  1x,'X ',3f15.7,/
     5  1x,'Y ',3f15.7,/
     6  1x,'Z ',3f15.7,/
     3  1x,'-----------------------------------------------',/
     6  1x,'Eigenvalues = ',3f15.7,/
     6  1x,'Isotropic   = ',1f15.7,/
     6  1x,'Anisotropic = ',1f15.7,/
     7  1x,'-----------------------------------------------')
      
      return
      end
