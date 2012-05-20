      subroutine smd_ewald_get_alpha(aalpha)
      implicit none
#include "smd_ewald_data.fh"

      double precision aalpha

      aalpha = alpha

      end 

      subroutine smd_ewald_get_ralphsq(aralphsq)
      implicit none
#include "smd_ewald_data.fh"

      double precision aralphsq

      aralphsq = ralphsq

      end 

      subroutine smd_ewald_get_rksqmax(arksqmax)
      implicit none
#include "smd_ewald_data.fh"

      double precision arksqmax

      arksqmax = rksqmax

      end 

      subroutine smd_ewald_get_kmax(akmax)
      implicit none
#include "smd_ewald_data.fh"

      integer akmax(3)
      integer i

      do i=1,3
        akmax(i) = kmax(i)
      end do

      end 

      subroutine smd_ewald_set_alpha(aalpha)
      implicit none
#include "smd_ewald_data.fh"

      double precision aalpha

      alpha = aalpha

      end 

      subroutine smd_ewald_set_ralphsq(aralphsq)
      implicit none
#include "smd_ewald_data.fh"

      double precision aralphsq

      ralphsq = aralphsq

      end 

      subroutine smd_ewald_set_rksqmax(arksqmax)
      implicit none
#include "smd_ewald_data.fh"

      double precision arksqmax

      rksqmax = arksqmax

      end 

      subroutine smd_ewald_set_kmax(akmax)
      implicit none
#include "smd_ewald_data.fh"

      integer akmax(3)
      integer i

      do i=1,3
        kmax(i) = akmax(i)
      end do

      end 

      subroutine smd_ewald_print(un)
      implicit none
#include "smd_ewald_data.fh"

      integer un

      write(un,*) "ewald data"
      write(un,*) "kmax",kmax(1),kmax(2),kmax(3)
      write(un,*) "alpha",alpha
      write(un,*) "ralphsq",ralphsq
      write(un,*) "rksqmax",rksqmax

      end 
c $Id$
