*     ***********************
*     *                     *
*     *    set_katm_list    *
*     *                     *
*     ***********************
      subroutine set_katm_list(tag,nion,aname_function,
     >                         nkatm,katm,natm,atom)
      implicit none
      integer tag
      integer nion
      character*4 aname_function
      external    aname_function
      integer nkatm,katm(2),natm(2),atom(2)

#include "bafdecls.fh"
#include "errquit.fh"

      logical value
      integer i,ka,nkatm_tmp
      integer atom_list(2)
      character*4 aname
      character*20 rtdb_name

      integer     pspw_compare_atom
      external    pspw_compare_atom
      integer tid,Parallel_threadid
      external    Parallel_threadid

      tid = Parallel_threadid()

*     **** define nkatm, and atom ****
      value = BA_push_get(mt_byte,4*nion,
     >                    'atom_list',atom_list(2),atom_list(1))
      if (.not. value) call errquit('out of stack memory',0,
     &       MA_ERR)
      nkatm_tmp = 0
      do i=1,nion
         aname = aname_function(i)
         if (
     >    pspw_compare_atom(aname,nkatm_tmp,byte_mb(atom_list(1))).eq.0)
     >   then
           nkatm_tmp = nkatm_tmp + 1
           call pspw_insert_atom(aname,nkatm_tmp,byte_mb(atom_list(1)))
         end if
      end do
      nkatm = nkatm_tmp

      rtdb_name = 'atom'//CHAR(tag+ICHAR('a'))
      value = BA_alloc_get(mt_byte,(4*nkatm),rtdb_name,atom(2),atom(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)
      call pspw_copy_atom(nkatm,byte_mb(atom_list(1)),
     >                          byte_mb(atom(1)))

      value = BA_pop_stack(atom_list(2))
      if (.not. value) call errquit('popping stack memory',0, MA_ERR)

*     **** define katm, and natm ****
      rtdb_name = 'katm'//CHAR(tag+ICHAR('a'))
      value = BA_alloc_get(mt_int,(nion),rtdb_name,katm(2),katm(1))
      rtdb_name = 'natm'//CHAR(tag+ICHAR('a'))
      value = value.and.
     >        BA_alloc_get(mt_int,(nkatm),rtdb_name,natm(2),natm(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

!$OMP MASTER
      do i=1,nkatm
         int_mb(natm(1)+i-1) = 0
      end do
!!$OMP END MASTER
!!$cOMP BARRIER
      do i=1,nion
         aname = aname_function(i)
         ka=pspw_compare_atom(aname,nkatm,byte_mb(atom(1)))
         int_mb(katm(1)+i-1) = ka
!!$OMP MASTER
         int_mb(natm(1)+ka-1) = int_mb(natm(1)+ka-1)+1
!!$OMP END MASTER
      end do
!$OMP END MASTER
!$OMP BARRIER

      return
      end


      subroutine pspw_copy_atom(nkatm,alist1,alist2)
      implicit none
      integer nkatm
      character*4 alist1(nkatm)
      character*4 alist2(nkatm)
      integer i
      do i=1,nkatm
         alist2(i) = alist1(i)
      end do
      return
      end
      subroutine pspw_insert_atom(aname,nkatm,alist)
      implicit none
      character*4 aname
      integer nkatm
      character*4 alist(nkatm)
      alist(nkatm) = aname
      return
      end
      integer function pspw_compare_atom(aname,nkatm,alist)
      implicit none
      character*4 aname
      integer nkatm
      character*4 alist(nkatm)
      integer i,j
      j=0
      do i=1,nkatm
         if (aname.eq.alist(i))  j=i
      end do
      pspw_compare_atom = j
      return
      end



c $Id$
