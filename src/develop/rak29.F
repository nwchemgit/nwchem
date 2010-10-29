      logical function raktask_fullsc(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
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
      call rak_hondomapnwchem()
c
      if (.not.geom_create(geom,'geometry')) call errquit
     &    ('geom create failed',911, GEOM_ERR)
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) call errquit
     &    ('geom_rtdb_load failed',911, RTDB_ERR)
c
      if (.not.bas_create(basis,'ao basis')) call errquit
     &    ('bas_create failed',911, BASIS_ERR)
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis')) call errquit
     &    ('bas_rtdb_load failed',911, RTDB_ERR)
c
      write(luout,*)' geom/basis loaded'
      if (.not. bas_print(basis))
     $    call errquit(' basis print failed', 0, BASIS_ERR)
      call util_flush(luout)
c
* note byproduct of resetting ncart and nsph with each call
      call spcart_get_full_cs('n',ma_hp,ncart,nsph,
     &    .true.,.false.,basis)
      if (.not.ma_get_index(ma_hp,ma_kp))
     &    call errquit('rak29: ma get index failed',911, MA_ERR)
      if (.not.ma_alloc_get(mt_dbl,(ncart*ncart),
     &    'product difference  matrix',ma_hdiff,ma_kdiff))
     &    call errquit('rak29: ma alloc failed ',911, MA_ERR)
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
     &    call errquit('rak29: ma get index failed',911, MA_ERR)
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
     &    call errquit('rak29: ma get index failed',911, MA_ERR)
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
     &    call errquit('rak29: ma alloc failed ',911, MA_ERR)
      if (.not.ma_alloc_get(mt_dbl,(ncart*ncart),
     &    'second product cart,sph*sph,cart matrix',ma_hpp2,ma_kpp2))
     &    call errquit('rak29: ma alloc failed ',911, MA_ERR)
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
      subroutine rak_hondomapnwchem()
      implicit none
#include "mafdecls.fh"
      logical rak_is_unit_matrix
      external rak_is_unit_matrix

      integer map_hton(35)
      integer map_ntoh(35)
      double precision n2h_s(1,1)
      double precision n2h_p(3,3)
      double precision n2h_sp(4,4)
      double precision n2h_d(6,6)
      double precision n2h_f(10,10)
      double precision n2h_g(15,15)
      double precision h2n_s(1,1)
      double precision h2n_p(3,3)
      double precision h2n_sp(4,4)
      double precision h2n_d(6,6)
      double precision h2n_f(10,10)
      double precision h2n_g(15,15)
      integer i,j,ang, shell_off
      integer l2ang
      integer h_tmp, k_tmp
      integer sizei, sized
      integer h_n2hind, k_n2hind
      integer h_n2hcmp, k_n2hcmp
      integer coffset, ioffset, count
      integer ii,jj,rank,square
c
c     ----- mapping of shell components from nwchem to hondo -----
c
*0 1 2  3  4
*s p d  f  g
*1,3,6,10,15
      data map_hton /1,     2,     3,     4,
     1               5,     8,    10,     6,     7,    9,
     2              11,    17,    20,    12,    13,
     3              14,    18,    16,    19,    15,
     4              21,    31,    35,    22,    23,
     5              27,    32,    30,    34,    24,
     6              26,    33,    25,    28,    29 /
*
      square(ii,jj,rank) = (ii-1)*rank + jj - 1 ! 0 to (rank*rank)-1)
*
      do i = 1,35
        write(6,10000)i,map_hton(i)
      enddo
      do i = 1,35
        map_ntoh(map_hton(i)) = i
      enddo
      do i = 1,35
        write(6,10001)i,map_ntoh(i)
      enddo
* zero fill matrices 
      call dfill(3*3,0.0d00,n2h_p,1)
      call dfill(4*4,0.0d00,n2h_sp,1)
      call dfill(6*6,0.0d00,n2h_d,1)
      call dfill(10*10,0.0d00,n2h_f,1)
      call dfill(15*15,0.0d00,n2h_g,1)
      call dfill(3*3,0.0d00,h2n_p,1)
      call dfill(4*4,0.0d00,h2n_sp,1)
      call dfill(6*6,0.0d00,h2n_d,1)
      call dfill(10*10,0.0d00,h2n_f,1)
      call dfill(15*15,0.0d00,h2n_g,1)
* fill matrices as nwchem 2 hondo order e.g., unit matrix with map
      n2h_s(1,map_ntoh(1)) = 1.0d00
      h2n_s(1,map_hton(1)) = 1.0d00
      do i=1,3
        n2h_p(i,map_ntoh(i+1)-1)  = 1.0d00
        h2n_p(i,map_hton(i+1)-1)  = 1.0d00
      enddo
      do i=1,4
        n2h_sp(i,map_ntoh(i)) = 1.0d00
        h2n_sp(i,map_hton(i)) = 1.0d00
      enddo
      do i=1,6
        n2h_d(i,map_ntoh(i+4)-4)  = 1.0d00
        h2n_d(i,map_hton(i+4)-4)  = 1.0d00
      enddo
      do i=1,10
        n2h_f(i,map_ntoh(i+10)-10)  = 1.0d00
        h2n_f(i,map_hton(i+10)-10)  = 1.0d00
      enddo
      do i=1,15
        n2h_g(i,map_ntoh(i+20)-20)  = 1.0d00
        h2n_g(i,map_hton(i+20)-20)  = 1.0d00
      enddo
      write(6,*)' n2h_s '
      call output(n2h_s,1,1,1,1,1,1,1)
      write(6,*)' n2h_p '
      call output(n2h_p,1,3,1,3,3,3,1)
      write(6,*)' n2h_sp '
      call output(n2h_sp,1,4,1,4,4,4,1)
      write(6,*)' n2h_d '
      call output(n2h_d,1,6,1,6,6,6,1)
      write(6,*)' n2h_f '
      call output(n2h_f,1,10,1,10,10,10,1)
      write(6,*)' n2h_g '
      call output(n2h_g,1,15,1,15,15,15,1)

      write(6,*)' h2n_s '
      call output(h2n_s,1,1,1,1,1,1,1)
      write(6,*)' h2n_p '
      call output(h2n_p,1,3,1,3,3,3,1)
      write(6,*)' h2n_sp '
      call output(h2n_sp,1,4,1,4,4,4,1)
      write(6,*)' h2n_d '
      call output(h2n_d,1,6,1,6,6,6,1)
      write(6,*)' h2n_f '
      call output(h2n_f,1,10,1,10,10,10,1)
      write(6,*)' h2n_g '
      call output(h2n_g,1,15,1,15,15,15,1)

*allocate memory
      if (.not.ma_alloc_get(mt_dbl,(15*15),'rak29:tmp',
     &      h_tmp, k_tmp)) stop 'h_tmp alloc failed'
* products of p
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',3,3,3,1.0d00,n2h_p,3,h2n_p,3,
     &      0.0d00,dbl_mb(k_tmp),3)
      write(6,*)' product matrix p,nhhn',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),3,3)
      call output(dbl_mb(k_tmp),1,3,1,3,3,3,1)
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',3,3,3,1.0d00,h2n_p,3,n2h_p,3,
     &      0.0d00,dbl_mb(k_tmp),3)
      write(6,*)' product matrix p,hnnh',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),3,3)
      call output(dbl_mb(k_tmp),1,3,1,3,3,3,1)
* products of sp
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',4,4,4,1.0d00,n2h_sp,4,h2n_sp,4,
     &      0.0d00,dbl_mb(k_tmp),4)
      write(6,*)' product matrix sp,nhhn',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),4,4)
      call output(dbl_mb(k_tmp),1,4,1,4,4,4,1)
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',4,4,4,1.0d00,h2n_sp,4,n2h_sp,4,
     &      0.0d00,dbl_mb(k_tmp),4)
      write(6,*)' product matrix sp,hnnh',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),4,4)
      call output(dbl_mb(k_tmp),1,4,1,4,4,4,1)
* products of d
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',6,6,6,1.0d00,n2h_d,6,h2n_d,6,
     &      0.0d00,dbl_mb(k_tmp),6)
      write(6,*)' product matrix d,nhhn',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),6,6)
      call output(dbl_mb(k_tmp),1,6,1,6,6,6,1)
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',6,6,6,1.0d00,h2n_d,6,n2h_d,6,
     &      0.0d00,dbl_mb(k_tmp),6)
      write(6,*)' product matrix d,hnnh',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),6,6)
      call output(dbl_mb(k_tmp),1,6,1,6,6,6,1)
* products of f
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',10,10,10,1.0d00,n2h_f,10,h2n_f,10,
     &      0.0d00,dbl_mb(k_tmp),10)
      write(6,*)' product matrix f,nhhn',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),10,10)
      call output(dbl_mb(k_tmp),1,10,1,10,10,10,1)
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',10,10,10,1.0d00,h2n_f,10,n2h_f,10,
     &      0.0d00,dbl_mb(k_tmp),10)
      write(6,*)' product matrix f,hnnh',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),10,10)
      call output(dbl_mb(k_tmp),1,10,1,10,10,10,1)
* products of g
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',15,15,15,1.0d00,n2h_g,15,h2n_g,15,
     &      0.0d00,dbl_mb(k_tmp),15)
      write(6,*)' product matrix g,nhhn',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),15,15)
      call output(dbl_mb(k_tmp),1,15,1,15,15,15,1)
      call dfill((15*15),0.0d00,dbl_mb(k_tmp),1)
      call dgemm('n','n',15,15,15,1.0d00,h2n_g,15,n2h_g,15,
     &      0.0d00,dbl_mb(k_tmp),15)
      write(6,*)' product matrix g,hnnh',
     &      rak_is_unit_matrix(dbl_mb(k_tmp),15,15)
      call output(dbl_mb(k_tmp),1,15,1,15,15,15,1)
* build ma based index and compressed matrix arrays
*                   0 1  2  3  4  5
*                   1 2  3  4  5  6
*                  -1 0  1  2  3  4
      sizei = 6 ! sp, s, p, d, f, g
      sized = 4*4 + 1*1 + 3*3 + 6*6 + 10*10 + 15*15
      if (.not.ma_alloc_get(mt_int,sizei,'n2h:index',
     &      h_n2hind,k_n2hind)) stop 'error'
      if (.not.ma_alloc_get(mt_dbl,sized,'n2h:compressed',
     &      h_n2hcmp,k_n2hcmp)) stop ' error'
      call ifill(sizei,0,int_mb(k_n2hind),1)
      call dfill(sized,0.0d00,dbl_mb(k_n2hcmp),1)
      coffset = 0
      ioffset = 0
      int_mb(k_n2hind+ioffset) = coffset
*:sp
      do i=1,4
        write(6,*)' i = ',i,' square = ',square(i,map_ntoh(i),4),
     &      ' map =',map_ntoh(i)
        dbl_mb(k_n2hcmp+coffset+square(i,map_ntoh(i),4)) = 1.0d00
*        n2h_sp(i,map_ntoh(i)) = 1.0d00
*        h2n_sp(i,map_hton(i)) = 1.0d00
      enddo
      coffset = coffset + 4*4
      ioffset = ioffset + 1
*:s
      int_mb(k_n2hind+ioffset) = coffset
      write(6,*)' coffset ',coffset
      dbl_mb(k_n2hcmp+coffset+square(1,1,1)) = 2.0d00
*      n2h_s(1,map_ntoh(1)) = 1.0d00
*      h2n_s(1,map_hton(1)) = 1.0d00
      coffset = coffset + 1
      ioffset = ioffset + 1
      int_mb(k_n2hind+ioffset) = coffset
      write(6,*)' coffset ',coffset
*:p
      do i=1,3
        write(6,*)' i = ',i,' square = ',square(i,map_ntoh(i),3),
     &      ' map =',map_ntoh(i)
        dbl_mb(k_n2hcmp+coffset+square(i,map_ntoh(i),3)) = 3.0d00        
*        n2h_p(i,map_ntoh(i+1)-1)  = 1.0d00
*        h2n_p(i,map_hton(i+1)-1)  = 1.0d00
      enddo
      coffset = coffset + 3*3
      ioffset = ioffset + 1
      int_mb(k_n2hind+ioffset) = coffset
      write(6,*)' coffset ',coffset
*:d
      do i=1,6
        write(6,*)' i = ',i,' square = ',square(i,(map_ntoh(i+4)-4),6),
     &      ' map =',(map_ntoh(i+4)-4)
        dbl_mb(k_n2hcmp+coffset+square(i,(map_ntoh(i+4)-4),6)) = 4.0d00        
*        n2h_d(i,map_ntoh(i+4)-4)  = 1.0d00
*        h2n_d(i,map_hton(i+4)-4)  = 1.0d00
      enddo
      coffset = coffset + 6*6
      ioffset = ioffset + 1
      int_mb(k_n2hind+ioffset) = coffset
      write(6,*)' coffset ',coffset
*:f
      do i=1,10
        write(6,*)' i = ',i,' square = ',
     &      square(i,(map_ntoh(i+10)-10),10),
     &      ' map =',(map_ntoh(i+10)-10)
        dbl_mb(k_n2hcmp+coffset+square(i,(map_ntoh(i+10)-10),10)) =
     &      5.0d00        
*        n2h_f(i,map_ntoh(i+10)-10)  = 1.0d00
*        h2n_f(i,map_hton(i+10)-10)  = 1.0d00
      enddo
      coffset = coffset + 10*10
      ioffset = ioffset + 1
      int_mb(k_n2hind+ioffset) = coffset
      write(6,*)' coffset ',coffset
*:g
      do i=1,15
        write(6,*)' i = ',i,' square = ',
     &      square(i,(map_ntoh(i+20)-20),15),
     &      ' map =',(map_ntoh(i+20)-20)
        dbl_mb(k_n2hcmp+coffset+square(i,(map_ntoh(i+20)-20),15)) =
     &      6.0d00        
*        n2h_g(i,map_ntoh(i+20)-20)  = 1.0d00
*        h2n_g(i,map_hton(i+20)-20)  = 1.0d00
      enddo
      

      write(6,*)' fixed array sp ',int_mb(k_n2hind+0)
      call output(dbl_mb(k_n2hcmp+int_mb(k_n2hind+0)),1,4,1,4,4,4,1)
      write(6,*)' fixed array s  ',int_mb(k_n2hind+1)
      call output(dbl_mb(k_n2hcmp+int_mb(k_n2hind+1)),1,1,1,1,1,1,1)
      write(6,*)' fixed array p  ',int_mb(k_n2hind+2)
      call output(dbl_mb(k_n2hcmp+int_mb(k_n2hind+2)),1,3,1,3,3,3,1)
      write(6,*)' fixed array d  ',int_mb(k_n2hind+3)
      call output(dbl_mb(k_n2hcmp+int_mb(k_n2hind+3)),1,6,1,6,6,6,1)
      write(6,*)' fixed array f  ',int_mb(k_n2hind+4)
      call output(dbl_mb(k_n2hcmp+int_mb(k_n2hind+4)),1,10,1,10,10,10,1)
      write(6,*)' fixed array g  ',int_mb(k_n2hind+5)
      call output(dbl_mb(k_n2hcmp+int_mb(k_n2hind+5)),1,15,1,15,15,15,1)
      if (.not.ma_free_heap(h_tmp)) stop 'h_tmp free failed'
      if (.not.ma_free_heap(h_n2hind)) stop 'h_tmp free failed'
      if (.not.ma_free_heap(h_n2hcmp)) stop 'h_tmp free failed'
10000 format(1x,'hondo  component ',i5,' maps to nwchem component ',i5)
10001 format(1x,'nwchem component ',i5,' maps to hondo  component ',i5)
*
      end
