*     ***********************
*     *                     *
*     *    set_katm_list    *
*     *                     *
*     ***********************
      subroutine set_katm_list(nion,aname_function,
     >                         nkatm,katm,natm,atom)
      implicit none
      integer nion
      character*2 aname_function
      external    aname_function
      integer nkatm,katm(2),natm(2),atom(2)

#include "mafdecls.fh"

      logical value
      integer i,ka
      integer atom_list(2)
      character*2 aname

      integer     pspw_compare_atom
      external    pspw_compare_atom


*     **** define nkatm, and atom ****
      value = MA_push_get(mt_byte,2*nion,
     >                    'atom_list',atom_list(2),atom_list(1))
      if (.not. value) call errquit('out of stack memory',0)
      nkatm = 0
      do i=1,nion
         aname = aname_function(i)
         if (
     >    pspw_compare_atom(aname,nkatm,byte_mb(atom_list(1))).eq.0)
     >   then
           nkatm = nkatm + 1
           call pspw_insert_atom(aname,nkatm,byte_mb(atom_list(1)))
         end if
      end do
      value = MA_alloc_get(mt_byte,(2*nkatm),'atom',atom(2),atom(1))
      if (.not. value) call errquit('out of heap memory',0)
      call pspw_copy_atom(nkatm,byte_mb(atom_list(1)),
     >                          byte_mb(atom(1)))

      value = MA_pop_stack(atom_list(2))
      if (.not. value) call errquit('popping stack memory',0)

*     **** define katm, and natm ****
      value = MA_alloc_get(mt_int,(nion),'katm',katm(2),katm(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(nkatm),'natm',natm(2),natm(1))
      if (.not. value) call errquit('out of heap memory',0)

      do i=1,nkatm
         int_mb(natm(1)+i-1) = 0
      end do
      do i=1,nion
         aname = aname_function(i)
         ka=pspw_compare_atom(aname,nkatm,byte_mb(atom(1)))
         int_mb(katm(1)+i-1) = ka
         int_mb(natm(1)+ka-1) = int_mb(natm(1)+ka-1)+1
      end do

      return
      end


      subroutine pspw_copy_atom(nkatm,alist1,alist2)
      implicit none
      integer nkatm
      character*2 alist1(nkatm)
      character*2 alist2(nkatm)
      integer i
      do i=1,nkatm
         alist2(i) = alist1(i)
      end do
      return
      end
      subroutine pspw_insert_atom(aname,nkatm,alist)
      implicit none
      character*2 aname
      integer nkatm
      character*2 alist(nkatm)
      alist(nkatm) = aname
      return
      end
      integer function pspw_compare_atom(aname,nkatm,alist)
      implicit none
      character*2 aname
      integer nkatm
      character*2 alist(nkatm)
      integer i,j
      j=0
      do i=1,nkatm
         if (aname.eq.alist(i))  j=i
      end do
      pspw_compare_atom = j
      return
      end



