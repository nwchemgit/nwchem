      logical function raktask_fullsc(rtdb)
      implicit none
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
c
      double precision ddot
      external ddot
      logical rak_is_unit_matrix
      external rak_is_unit_matrix
c
      integer basis, geom
      integer ma_hp , ma_kp
      integer ma_hpp, ma_kpp
      integer ma_hpp2, ma_kpp2
      integer ma_hdiff, ma_kdiff
      integer ma_hptrans, ma_kptrans
      integer ma_hpinv, ma_kpinv
      double precision norm
      integer ncart, nsph
      integer ga_hp
c
      raktask_fullsc = .false.
c
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('geom create failed',911)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('geom_rtdb_load failed',911)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &    ('bas_create failed',911)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &    ('bas_rtdb_load failed',911)
c
      write(luout,*)' geom/basis loaded'
      if (.not. bas_print(basis))
     $    call errquit(' basis print failed', 0)
      call util_flush(luout)
c
* note byproduct of resetting ncart and nsph with each call
      call spcart_get_full_cs('n',ma_hp,ncart,nsph,
     &    .true.,.false.,basis)
      if (.not.ma_get_index(ma_hp,ma_kp))
     &    call errquit('rak29: ma get index failed',911)
      if (.not.ma_alloc_get(mt_dbl,(ncart*ncart),
     &    'product difference  matrix',ma_hdiff,ma_kdiff))
     &    call errquit('rak29: ma alloc failed ',911)
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kdiff),1)
*-ga test
      call spcart_get_full_cs('n',ga_hp,ncart,nsph,
     &    .false.,.true.,basis)
      call ga_get(ga_hp,1,ncart,1,nsph,dbl_mb(ma_kdiff),ncart)
      if (.not.ga_destroy(ga_hp)) stop ' ga destroy error'
      call daxpy((ncart*nsph),-1.0d00,dbl_mb(ma_kp),1,
     &    dbl_mb(ma_kdiff),1)
      norm = ddot((ncart*ncart),dbl_mb(ma_kdiff),1,dbl_mb(ma_kdiff),1)
      write(luout,*)'rak29: new/old from ga difference norm ',norm
* note byproduct of resetting ncart and nsph with each call
      call spcart_get_full_cs('t',ma_hptrans,ncart,nsph,
     &    .true.,.false.,basis)
      if (.not.ma_get_index(ma_hptrans,ma_kptrans))
     &    call errquit('rak29: ma get index failed',911)
*-ga test
      call spcart_get_full_cs('t',ga_hp,ncart,nsph,
     &    .false.,.true.,basis)
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kdiff),1)
      call ga_get(ga_hp,1,nsph,1,ncart,dbl_mb(ma_kdiff),nsph)
      if (.not.ga_destroy(ga_hp)) stop ' ga destroy error'
      call daxpy((ncart*nsph),-1.0d00,dbl_mb(ma_kptrans),1,
     &    dbl_mb(ma_kdiff),1)
      norm = ddot((ncart*ncart),dbl_mb(ma_kdiff),1,dbl_mb(ma_kdiff),1)
      write(luout,*)'rak29: new/old trans from ga difference norm ',norm
* note byproduct of resetting ncart and nsph with each call
      call spcart_get_full_cs('i',ma_hpinv,ncart,nsph,
     &    .true.,.false.,basis)
      if (.not.ma_get_index(ma_hpinv,ma_kpinv))
     &    call errquit('rak29: ma get index failed',911)
*-ga test
      call spcart_get_full_cs('i',ga_hp,ncart,nsph,
     &    .false.,.true.,basis)
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kdiff),1)
      call ga_get(ga_hp,1,nsph,1,ncart,dbl_mb(ma_kdiff),nsph)
      call daxpy((ncart*nsph),-1.0d00,dbl_mb(ma_kpinv),1,
     &    dbl_mb(ma_kdiff),1)
      norm = ddot((ncart*ncart),dbl_mb(ma_kdiff),1,dbl_mb(ma_kdiff),1)
      write(luout,*)'rak29: new/old inverse from ga difference norm ',
     &    norm
*
      if (.not.ma_alloc_get(mt_dbl,(ncart*ncart),
     &    'product cart,sph*sph,cart matrix',ma_hpp,ma_kpp))
     &    call errquit('rak29: ma alloc failed ',911)
      if (.not.ma_alloc_get(mt_dbl,(ncart*ncart),
     &    'second product cart,sph*sph,cart matrix',ma_hpp2,ma_kpp2))
     &    call errquit('rak29: ma alloc failed ',911)
c
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kpp),1)
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kpp2),1)
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kdiff),1)
      call dgemm('N','T',ncart,ncart,nsph,1.0d00,
     &    dbl_mb(ma_kp),ncart,
     &    dbl_mb(ma_kp),ncart,
     &    0.0d00,dbl_mb(ma_kpp),ncart)
      call dgemm('N','N',ncart,ncart,nsph,1.0d00,
     &    dbl_mb(ma_kp),ncart,
     &    dbl_mb(ma_kptrans),nsph,
     &    0.0d00,dbl_mb(ma_kpp2),ncart)
c
      call dcopy((ncart*ncart),dbl_mb(ma_kpp),1,dbl_mb(ma_kdiff),1)
      call daxpy((ncart*ncart),-1.0d00,dbl_mb(ma_kpp2),1,
     &    dbl_mb(ma_kdiff),1)
      norm = ddot((ncart*ncart),dbl_mb(ma_kdiff),1,dbl_mb(ma_kdiff),1)
c
      write(luout,*)'rak29: product matrix'
      call output(dbl_mb(ma_kpp),1,ncart,1,ncart,ncart,ncart,1)
      write(luout,*)'rak29: second product matrix'
      call output(dbl_mb(ma_kpp2),1,ncart,1,ncart,ncart,ncart,1)
      write(luout,*)'rak29: product difference matrix'
      call output(dbl_mb(ma_kdiff),1,ncart,1,ncart,ncart,ncart,1)
      write(luout,*)'rak29: product difference norm ',norm
c
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kpp),1)
      call dfill((ncart*ncart),0.0d00,dbl_mb(ma_kpp2),1)
      call dgemm('N','N',ncart,ncart,nsph,1.0d00,
     &    dbl_mb(ma_kp),ncart,
     &    dbl_mb(ma_kpinv),nsph,
     &    0.0d00,dbl_mb(ma_kpp),ncart)
      write(luout,*)' matrix*(matrix**-1)'
      call output(dbl_mb(ma_kpp),1,ncart,1,ncart,ncart,ncart,1)
      if (rak_is_unit_matrix(dbl_mb(ma_kpp),nsph,nsph)) then
        write(luout,*)' is unit matrix'
      else
        write(luout,*)' is NNOOTT unit matrix'
      endif
      call dgemm('n','n',nsph,nsph,ncart,1.0d00,
     &    dbl_mb(ma_kpinv),nsph,
     &    dbl_mb(ma_kp),ncart,
     &    0.0d00,dbl_mb(ma_kpp2),nsph)
      write(luout,*)' (matrix**-1)*matrix'
      call output(dbl_mb(ma_kpp2),1,nsph,1,nsph,nsph,nsph,1)
      if (rak_is_unit_matrix(dbl_mb(ma_kpp2),nsph,nsph)) then
        write(luout,*)' is unit matrix'
      else
        write(luout,*)' is NNOOTT unit matrix'
      endif
c
      raktask_fullsc = ma_free_heap(ma_hp)
      write(luout,*)' raktask_fullsc 1    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.ga_destroy(ga_hp)
      write(luout,*)' raktask_fullsc 2    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.ma_free_heap(ma_hptrans)
      write(luout,*)' raktask_fullsc 3    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.ma_free_heap(ma_hpinv)
      write(luout,*)' raktask_fullsc 4    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.ma_free_heap(ma_hpp)
      write(luout,*)' raktask_fullsc 5    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.ma_free_heap(ma_hpp2)
      write(luout,*)' raktask_fullsc 6    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.ma_free_heap(ma_hdiff)
      write(luout,*)' raktask_fullsc 7    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.bas_destroy(basis)
      write(luout,*)' raktask_fullsc 8    ',raktask_fullsc
      raktask_fullsc = raktask_fullsc.and.geom_destroy(geom)
      write(luout,*)' raktask_fullsc last ',raktask_fullsc
c
      end
      logical function rak_is_unit_matrix(mat,nrow,ncol)
      implicit none
      integer nrow, ncol
      double precision mat(nrow,ncol)
*
      integer ir, ic
      double precision mytest
*
      double precision thresh
      parameter (thresh=1.0d-6)
c::-statement function
      logical is_it_close_to  
      double precision value,test
*---          is value close to test?
      is_it_close_to(value,test) = (abs(value-test).lt.thresh)
*
      rak_is_unit_matrix = .true.
      do ir = 1,nrow
        do ic = 1,ncol
          if (ir.eq.ic) then
            mytest = 1.0d00
          else
            mytest = 0.0d00
          endif
          if (.not.(is_it_close_to(mat(ir,ic),mytest))) then
            rak_is_unit_matrix = .false.
            return
          endif
        enddo
      enddo
      end
