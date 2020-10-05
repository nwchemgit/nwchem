
      subroutine cpsi_printeigall(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint
      integer nb,i
      integer ne(2),ne_excited(2),neall(2),isize
      integer eigs(2)
      real*8  f1,f2,f3,f4,f5,f6,ee
      character*255 full_filename

*     **** external functions ****
      integer  cpsi_ispin,cpsi_ne,control_excited_ne
      external cpsi_ispin,cpsi_ne,control_excited_ne
      integer  brillioun_nbrillioun
      external brillioun_nbrillioun
      real*8   brillioun_k_brdcst,brillioun_ks_brdcst
      external brillioun_k_brdcst,brillioun_ks_brdcst
      real*8   cpsi_eigenvalue_brdcst,cpsi_virtual_brdcst
      external cpsi_eigenvalue_brdcst,cpsi_virtual_brdcst

      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      ne_excited(1) = 0
      ne_excited(2) = 0
      ne_excited(1) = control_excited_ne(1)
      if (cpsi_ispin().eq.2) ne_excited(2) = control_excited_ne(2)
      neall(1) = ne(1) + ne_excited(1)
      neall(2) = ne(2) + ne_excited(2)

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

      call util_file_name('eigenvalue_grid',
     >                    .false.,
     >                    .false.,
     >                    full_filename)
      if (taskid.eq.MASTER) 
     >   open(unit=58,file=full_filename,form='formatted')


*     **** allocate eigs ****
      isize = neall(1)+neall(2)
      if (.not.BA_push_get(mt_dbl,isize,'eigs',eigs(2),eigs(1)))
     >  call errquit('cpsi_printeigall:out of stack memory',0,MA_ERR)


      do nb=1,brillioun_nbrillioun()
         !f0 = brillioun_weight_brdcst(nb)
         f1 = brillioun_ks_brdcst(1,nb)
         f2 = brillioun_ks_brdcst(2,nb)
         f3 = brillioun_ks_brdcst(3,nb)
         f4 = brillioun_k_brdcst(1,nb)
         f5 = brillioun_k_brdcst(2,nb)
         f6 = brillioun_k_brdcst(3,nb)
         do i=1,neall(1)
            if (i.le.ne(1)) then
               ee = cpsi_eigenvalue_brdcst(nb,1,i)
            else
               ee = cpsi_virtual_brdcst(nb,1,i-ne(1))
            end if
            dbl_mb(eigs(1)+i-1) = ee
         end do
         do i=1,neall(2)
            if (i.le.ne(2)) then
               ee = cpsi_eigenvalue_brdcst(nb,2,i)
            else
               ee = cpsi_virtual_brdcst(nb,2,i-ne(2))
            end if
            dbl_mb(eigs(1)+neall(1)+i-1) = ee
         end do
         if (taskid.eq.MASTER) then
            write(58,*) nb,f1,f2,f3,f4,f5,f6,
     >                (dbl_mb(eigs(1)+i-1),i=1,neall(1)+neall(2))
            write(58,*) nb,-f1,-f2,-f3,-f4,-f5,-f6,
     >                (dbl_mb(eigs(1)+i-1),i=1,neall(1)+neall(2))
         end if
      end do


      if (taskid.eq.MASTER) close(58)


      if (.not.BA_pop_stack(eigs(2)))
     >  call errquit('cpsi_printeigall:popping stack memory',0,MA_ERR)

      return
      end 
