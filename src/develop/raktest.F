      subroutine raktest(rtdb)
      implicit none
      integer rtdb
c
      write(6,*)' run rjh test '
c      call rjhtest(rtdb)
      write(6,*)' rjh test done '
      write(6,*) ' rak test'
      call dointsum(rtdb,'geometry', 'mo basis')
c
      end
      subroutine dointsum(rtdb,geom_name, basis_name)
      implicit none
#include "geom.fh"
#include "bas.fh"
      logical int_normalize
      external int_normalize
c
      integer rtdb
      integer basis, geom
      character*(*) geom_name, basis_name
c
      write(6,*)' rtdb      :',rtdb
      write(6,*)' geom_name : <',geom_name,'>'
      write(6,*)' basis_name: <',basis_name,'>'
c
      if (.not. geom_create(geom, geom_name))
     $     call errquit('dointsum: geom_create failed?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
     $     call errquit('dointsum: geom_load failed', 0)
      if (.not. geom_print(geom))
     $     call errquit('dointsum: geom_print failed', 0)
      if (.not. bas_create(basis, basis_name))
     $     call errquit('dointsum: basis create failed', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
     $      call errquit('dointsum: basis load failed', 0)
      if (.not. bas_print(basis))
     $      call errquit('dointsum: basis print failed', 0)
      if (.not. gbs_map_print(basis))
     $     call errquit('dointsum: gbs_map_print failed', 0)
c
c normalize basis set
c
      if (.not.int_normalize(basis))
     &       call errquit('dointsum: basis norm. failed', 0)
      write(6,*)' after normalization'
      if (.not. bas_print(basis))
     $      call errquit('dointsum: basis print failed', 0)
c
      call int_chk_sum(basis)
      if (.not.(
     &       (bas_destroy(basis)).and.(geom_destroy(geom))
     &       ))
     &       call errquit('rak:error destroying gbs pair',0)
      end
