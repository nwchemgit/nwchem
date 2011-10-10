      subroutine chelp_grid_data_alloc_put(nr,r)
      implicit none
#include "chelp_grid_data.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer nr
      double precision r(nr)

      if(.not.ma_alloc_get(mt_dbl,nr,'chelp grid',
     + h_rchelp,i_rchelp))
     + call errquit('Failed to allocate memory chelp grid',nr,
     &       MA_ERR)
   
      call dcopy(nr,r,1,dbl_mb(i_rchelp),1)
      nrchelp=nr
      end 
c      
      subroutine chelp_grid_data_get_prt(i_r)
      implicit none
#include "chelp_grid_data.fh"
      integer i_r
      i_r=i_rchelp
      end 

      subroutine chelp_grid_data_print()
      implicit none
#include "chelp_grid_data.fh"
#include "mafdecls.fh"
#include "global.fh"
c      
      integer ir,k
      integer un
c      
      if(ga_nodeid().eq.0) then
      call util_file_open(un,.true.,"chelp.grid.xyz")
      write(un,*) nrchelp/3
      write(un,*) "grid"
      do ir=1,nrchelp/3
        write(un,*) "Bq", 
     +    (dbl_mb(i_rchelp+(ir-1)*3+k-1),k=1,3)
      end do
      end if
      call ga_sync()
c                                                                     
      return                                                          
      end                                                             

c $Id$
