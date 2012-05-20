c      subroutine smd_ewald_get_alpha(aalpha)
c      implicit none
c#include "smd_ewald_data.fh"
c
c      double precision aalpha
c
c      aalpha = alpha
c
c      end 
c
c      subroutine smd_ewald_get_ralphsq(aralphsq)
c      implicit none
c#include "smd_ewald_data.fh"
c
c      double precision aralphsq
c
c      aralphsq = ralphsq
c
c      end 
c
c      subroutine smd_ewald_get_rksqmax(arksqmax)
c      implicit none
c#include "smd_ewald_data.fh"
c
c      double precision arksqmax
c
c      arksqmax = rksqmax
c
c      end 
c
c      subroutine smd_ewald_get_kmax(akmax)
c      implicit none
c#include "smd_ewald_data.fh"
c
c      integer akmax(3)
c      integer i
c
c      do i=1,3
c        akmax(i) = kmax(i)
c      end do
c
c      end 
c
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
