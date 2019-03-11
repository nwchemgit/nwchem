*
* $Id$
*

*     ***************************
*     *          	   	*
*     *       cpsi_read		*
*     *             		*
*     ***************************
      subroutine cpsi_read(spin_orbit,ispin,ne,nbrillioun,psi2_tag)
      implicit none 
      logical spin_orbit
      integer ispin,ne(2),nbrillioun
      integer psi2_tag

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"

*    *** local variables ***
      integer rtdb
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

*     **** local variables ****
      logical hprint
      integer MASTER,taskid,taskid_j,taskid_k
      parameter(MASTER=0)
      logical brillioun_check,psi_nolattice
      integer n,nb,nbrill,occupation,ms
      integer msglen,psiold,psiold1
      integer q,pj,nbq,pk

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2),tmpe(2)
      logical value

*     ***** local functions ****
      logical control_spin_orbit,control_print
      external control_spin_orbit,control_print
      character*50 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb,control_fractional_orbitals
      external control_ngrid,control_rtdb,control_fractional_orbitals
      integer  cpsi_data_get_allptr,cpsi_data_get_next
      integer  cpsi_data_get_chnk
      external cpsi_data_get_allptr,cpsi_data_get_next
      external cpsi_data_get_chnk

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
       if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,(ne(1)+ne(2)),'tmpe',tmpe(2),tmpe(1))
       if (.not. value) call errquit('out of stack memory',0, MA_ERR)
     
*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_psi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                        full_filename)
         l = index(full_filename,' ') -1
         if (control_print(print_medium)) then
           write(luout,1210) full_filename(1:l)
         end if
 1210    FORMAT(/' input psi filename:',A)

         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call iread(5,ispin,1)
         call iread(5,ne,2)
         call iread(5,nbrill,1)
         call iread(5,occupation,1)
      end if
      

c     **** send header to all nodes ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,ispin)
      msglen = 2
      call Parallel_Brdcst_ivalues(MASTER,msglen,ne)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,nbrill)
      call Parallel_Brdcst_ivalues(MASTER,msglen,occupation)


c     **** spin_orbit - ispin=3 --> ispin=2,spin_orbit=.true ****
      if (ispin.eq.3) then
         ispin = 2
         spin_orbit = .true.
      else
         spin_orbit = .false.
      end if


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        call errquit('Error reading wavefunctions', 0, CALC_ERR)
      end if

      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nolattice',
     >                  mt_log,1,psi_nolattice))
     >   psi_nolattice = .true.

      if (.not.psi_nolattice) then
      if ( (dabs(unita(1,1)-control_unita(1,1)).gt.1.0d-6) .or.
     >     (dabs(unita(2,1)-control_unita(2,1)).gt.1.0d-6) .or.
     >     (dabs(unita(3,1)-control_unita(3,1)).gt.1.0d-6) .or.
     >     (dabs(unita(1,2)-control_unita(1,2)).gt.1.0d-6) .or.
     >     (dabs(unita(2,2)-control_unita(2,2)).gt.1.0d-6) .or.
     >     (dabs(unita(3,2)-control_unita(3,2)).gt.1.0d-6) .or.
     >     (dabs(unita(1,3)-control_unita(1,3)).gt.1.0d-6) .or.
     >     (dabs(unita(2,3)-control_unita(2,3)).gt.1.0d-6) .or.
     >     (dabs(unita(3,3)-control_unita(3,3)).gt.1.0d-6) ) then
            if (hprint) then
            write(luout,*) "unita(psi)      unita(lattice)"
            write(luout,*) unita(1,1), control_unita(1,1)
            write(luout,*) unita(2,1), control_unita(2,1)
            write(luout,*) unita(3,1), control_unita(3,1)
            write(luout,*) unita(1,2), control_unita(1,2)
            write(luout,*) unita(2,2), control_unita(2,2)
            write(luout,*) unita(3,2), control_unita(3,2)
            write(luout,*) unita(1,3), control_unita(1,3)
            write(luout,*) unita(2,3), control_unita(2,3)
            write(luout,*) unita(3,3), control_unita(3,3)
            endif
        call errquit('Error reading wavefunctions', 1, CALC_ERR)
      end if
      end if

      if (.not.btdb_get(rtdb,'nwpw:psi_brillioun_check',
     >                  mt_log,1,brillioun_check))
     >   brillioun_check = .false.

      if (brillioun_check) then
      if ( nbrill.ne.nbrillioun) then
        call errquit(
     >   'Error reading wavefunctions - wrong number of k-points', 1,
     &       CALC_ERR)
      end if
      end if


*     *************************
*     *** read in 3d blocks ***
*     *************************

*     * file may contain more data then necessary *
      if (nbrillioun.le.nbrill) then

        do nb=1,nbrillioun
           call K1dB_ktoqp(nb,nbq,pk)
           do n=1,(ne(1)+ne(2))
              pj=0
              q=n
              call C3dB_c_Read(1,5,dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp(1)),pj,pk)

            if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
               call Cram_c_pack(nbq,dcpl_mb(tmp2(1)))
               call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp2(1)))
            end if
           end do
        end do

*     * file contains less data then necessary *
      else 

        do nb=1,nbrill-1
           call K1dB_ktoqp(nb,nbq,pk)
           do n=1,(ne(1)+ne(2))
              pj=0
              q=n

             call C3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)),pj,pk)

             if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
                call Cram_c_pack(nbq,dcpl_mb(tmp2(1)))
                call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp2(1)))
             end if
        end do
        end do
*       **** use the last zone for filling in unknown zones ****
        do n=1,(ne(1)+ne(2))
          pj=0
          q=n
          call C3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)),pj,-1)
          do nb=nbrill,nbrillioun
             call K1dB_ktoqp(nb,nbq,pk)
             if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
                call dcopy(2*nfft3d,dcpl_mb(tmp2(1)),1,
     >                              dcpl_mb(tmp(1)),1)
                call Cram_c_pack(nbq,dcpl_mb(tmp(1)))
                call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp(1)))
             end if
          end do
        end do

      end if

*     **** read the occupations - still a bit hacky ****
      if (occupation.gt.0) then
         l = cpsi_data_get_next(psi2_tag)
*        * file may contain more data then necessary *
c        if (nbrillioun.eq.nbrill) then
         if (nbrillioun.le.nbrill) then

            do nb=1,nbrillioun
               call K1dB_ktoqp(nb,nbq,pk)
               if (taskid.eq.MASTER) then
                  call dread(5,dbl_mb(tmpe(1)),(ne(1)+ne(2)))
               end if
               !** big broadcast - should be reduced ***
               call Parallel_Brdcst_values(MASTER,
     >                       (ne(1)+ne(2)),dbl_mb(tmpe(1)))
               if (pk.eq.taskid_k) then
                  do n=1,(ne(1)+ne(2))
                     call cpsi_data_add(l,nbq,n,dbl_mb(tmpe(1)+n-1))
                  end do
               end if
            end do
*        **** reset to zeros and ones ****
*       * file contains less data then necessary *
         else
           do nb=1,nbrill
             call K1dB_ktoqp(nb,nbq,pk) 
             if (taskid.eq.MASTER) then
               call dread(5,dbl_mb(tmpe(1)),(ne(1)+ne(2)))
             end if
             call Parallel_Brdcst_values(MASTER,
     >                                   (ne(1)+ne(2)),dbl_mb(tmpe(1)))
             if (pk.eq.taskid_k) then 
               do n =1,ne(1)+ne(2)
                   call cpsi_data_add(l,nbq,n,dbl_mb(tmpe(1)+n-1))
                 !if (n.le.(ne(ms)-control_fractional_orbitals(ms))) then
                 !  unita(1,1) = 1.0d0  
                 !else
                 !  unita(1,1) = 0.0d0  
                 !end if
                 !call cpsi_data_add(l,nbq,(n+(ms-1)*ne(1)),unita)
               end do
             end if
           end do
           do nb=nbrill+1,nbrillioun
             call K1dB_ktoqp(nb,nbq,pk) 
             if (pk.eq.taskid_k) then 
               do n =1,ne(1)+ne(2)
                   call cpsi_data_add(l,nbq,n,dbl_mb(tmpe(1)+n-1))
               end do
             end if
           end do
         end if

      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value =           BA_pop_stack(tmpe(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
*     end if

      return
      end

*

*     ***********************************
*     *          	   		*
*     *       cpsi_read_filename	*
*     *             			*
*     ***********************************
      subroutine cpsi_read_filename(filename,
     >                spin_orbit,ispin,ne,nbrillioun,psi2_tag)
      implicit none 
      character*(*) filename
      logical spin_orbit
      integer ispin,ne(2),nbrillioun
      integer psi2_tag

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*    *** local variables ***
      integer rtdb
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*255 full_filename

*     **** local variables ****
      integer MASTER,taskid,taskid_j,taskid_k
      parameter(MASTER=0)
      logical brillioun_check,psi_nolattice
      integer n,nb,nbrill,occupation
      integer msglen,psiold,psiold1
      integer q,pj,nbq,pk

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value

*     ***** local functions ****
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb,cpsi_data_get_chnk
      external control_ngrid,control_rtdb,cpsi_data_get_chnk

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
       if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                        full_filename)
         l = index(full_filename,' ') -1
         write(*,1210) full_filename(1:l)
 1210    FORMAT(/' input psi filename:',A)

         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call iread(5,ispin,1)
         call iread(5,ne,2)
         call iread(5,nbrill,1)
         call iread(5,occupation,1)
      end if

c     **** send header to all nodes ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,ispin)
      msglen = 2
      call Parallel_Brdcst_ivalues(MASTER,msglen,ne)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,nbrill)
      call Parallel_Brdcst_ivalues(MASTER,msglen,occupation)

c     **** spin_orbit - ispin=3 --> ispin=2,spin_orbit=.true ****
      if (ispin.eq.3) then
         ispin = 2
         spin_orbit = .true.
      else
         spin_orbit = .false.
      end if

*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        call errquit('Error reading wavefunctions', 0, CALC_ERR)
      end if

      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nolattice',
     >                  mt_log,1,psi_nolattice))
     >   psi_nolattice = .true.

      if (.not.psi_nolattice) then
      if ( (dabs(unita(1,1)-control_unita(1,1)).gt.1.0d-6) .or.
     >     (dabs(unita(2,1)-control_unita(2,1)).gt.1.0d-6) .or.
     >     (dabs(unita(3,1)-control_unita(3,1)).gt.1.0d-6) .or.
     >     (dabs(unita(1,2)-control_unita(1,2)).gt.1.0d-6) .or.
     >     (dabs(unita(2,2)-control_unita(2,2)).gt.1.0d-6) .or.
     >     (dabs(unita(3,2)-control_unita(3,2)).gt.1.0d-6) .or.
     >     (dabs(unita(1,3)-control_unita(1,3)).gt.1.0d-6) .or.
     >     (dabs(unita(2,3)-control_unita(2,3)).gt.1.0d-6) .or.
     >     (dabs(unita(3,3)-control_unita(3,3)).gt.1.0d-6) ) then
            write(*,*) "unita(psi)      unita(lattice)"
            write(*,*) unita(1,1), control_unita(1,1)
            write(*,*) unita(2,1), control_unita(2,1)
            write(*,*) unita(3,1), control_unita(3,1)
            write(*,*) unita(1,2), control_unita(1,2)
            write(*,*) unita(2,2), control_unita(2,2)
            write(*,*) unita(3,2), control_unita(3,2)
            write(*,*) unita(1,3), control_unita(1,3)
            write(*,*) unita(2,3), control_unita(2,3)
            write(*,*) unita(3,3), control_unita(3,3)
        call errquit('Error reading wavefunctions', 1, CALC_ERR)
      end if
      end if

      if (.not.btdb_get(rtdb,'nwpw:psi_brillioun_check',
     >                  mt_log,1,brillioun_check))
     >   brillioun_check = .false.

      if (brillioun_check) then
      if ( nbrill.ne.nbrillioun) then
        call errquit(
     >   'Error reading wavefunctions - wrong number of k-points', 1,
     &       CALC_ERR)
      end if
      end if


*     *************************
*     *** read in 3d blocks ***
*     *************************

*     * file may contain more data then necessary *
      if (nbrillioun.le.nbrill) then

        do nb=1,nbrillioun
           call K1dB_ktoqp(nb,nbq,pk)
           do n=1,(ne(1)+ne(2))
              pj=0
              q=n
              call C3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp(1)),pj,pk)
              if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
                 call Cram_c_pack(nbq,dcpl_mb(tmp2(1)))
                 call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp2(1)))
              endif
           end do
        end do

*     * file contains less data then necessary *
      else 
        do nb=1,nbrill-1
           call K1dB_ktoqp(nb,nbq,pk)
           do n=1,(ne(1)+ne(2))
              pj=0
              q=n
              call C3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)),pj,pk)
              if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
                 call Cram_c_pack(nbq,dcpl_mb(tmp2(1)))
                 call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp2(1)))
              end if
           end do
        end do
*       **** use the last zone for filling in unknown zones ****
        do n=1,(ne(1)+ne(2))
           pj=0
           q=n
           call C3dB_c_Read(1,5,dcpl_mb(tmp2(1)),
     >                     dcpl_mb(tmp(1)),pj,-1)
           do nb=nbrill,nbrillioun
              call K1dB_ktoqp(nb,nbq,pk)
              if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
                 call dcopy(2*nfft3d,dcpl_mb(tmp2(1)),1,
     >                               dcpl_mb(tmp(1)),1)
                 call Cram_c_pack(nbq,dcpl_mb(tmp(1)))
                 call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp(1)))
              end if
           end do
        end do

      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))

*     end if

      return
      end


*     ***************************
*     *          	   	*
*     *       ecpsi_read	*
*     *             		*
*     ***************************
      subroutine ecpsi_read(spin_orbit,ispin,ne,nbrillioun,psi2_tag)
      implicit none 
      logical spin_orbit
      integer ispin,ne(2),nbrillioun
      integer psi2_tag

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"

*    *** local variables ***
      integer rtdb
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

*     **** local variables ****
      logical hprint
      integer MASTER,taskid,taskid_j,taskid_k
      parameter(MASTER=0)
      logical brillioun_check,psi_nolattice
      integer n,nb,nbrill,occupation,ms
      integer msglen,psiold,psiold1
      integer q,pj,nbq,pk

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2),tmpe(2)
      logical value

*     ***** local functions ****
      logical control_spin_orbit,control_print
      external control_spin_orbit,control_print
      character*50 control_input_epsi
      external     control_input_epsi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb,control_fractional_orbitals
      external control_ngrid,control_rtdb,control_fractional_orbitals
      integer  cpsi_data_get_allptr,cpsi_data_get_next
      integer  cpsi_data_get_chnk
      external cpsi_data_get_allptr,cpsi_data_get_next
      external cpsi_data_get_chnk

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      hprint = (taskid.eq.MASTER).and.control_print(print_high)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
       if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = BA_push_get(mt_dbl,(ne(1)+ne(2)),'tmpe',tmpe(2),tmpe(1))
       if (.not. value) call errquit('out of stack memory',0, MA_ERR)
     
*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_epsi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                        full_filename)
         l = index(full_filename,' ') -1
         if (control_print(print_medium))
     >     write(luout,1210) full_filename(1:l)
 1210    FORMAT(/' input psi filename:',A)

         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call iread(5,ispin,1)
         call iread(5,ne,2)
         call iread(5,nbrill,1)
         call iread(5,occupation,1)
      end if
      

c     **** send header to all nodes ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,ispin)
      msglen = 2
      call Parallel_Brdcst_ivalues(MASTER,msglen,ne)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,nbrill)
      call Parallel_Brdcst_ivalues(MASTER,msglen,occupation)


c     **** spin_orbit - ispin=3 --> ispin=2,spin_orbit=.true ****
      if (ispin.eq.3) then
         ispin = 2
         spin_orbit = .true.
      else
         spin_orbit = .false.
      end if


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        call errquit('Error reading wavefunctions', 0, CALC_ERR)
      end if

      rtdb = control_rtdb()
      if (.not.btdb_get(rtdb,'nwpw:psi_nolattice',
     >                  mt_log,1,psi_nolattice))
     >   psi_nolattice = .true.

      if (.not.psi_nolattice) then
      if ( (dabs(unita(1,1)-control_unita(1,1)).gt.1.0d-6) .or.
     >     (dabs(unita(2,1)-control_unita(2,1)).gt.1.0d-6) .or.
     >     (dabs(unita(3,1)-control_unita(3,1)).gt.1.0d-6) .or.
     >     (dabs(unita(1,2)-control_unita(1,2)).gt.1.0d-6) .or.
     >     (dabs(unita(2,2)-control_unita(2,2)).gt.1.0d-6) .or.
     >     (dabs(unita(3,2)-control_unita(3,2)).gt.1.0d-6) .or.
     >     (dabs(unita(1,3)-control_unita(1,3)).gt.1.0d-6) .or.
     >     (dabs(unita(2,3)-control_unita(2,3)).gt.1.0d-6) .or.
     >     (dabs(unita(3,3)-control_unita(3,3)).gt.1.0d-6) ) then
            if (hprint) then
            write(luout,*) "unita(psi)      unita(lattice)"
            write(luout,*) unita(1,1), control_unita(1,1)
            write(luout,*) unita(2,1), control_unita(2,1)
            write(luout,*) unita(3,1), control_unita(3,1)
            write(luout,*) unita(1,2), control_unita(1,2)
            write(luout,*) unita(2,2), control_unita(2,2)
            write(luout,*) unita(3,2), control_unita(3,2)
            write(luout,*) unita(1,3), control_unita(1,3)
            write(luout,*) unita(2,3), control_unita(2,3)
            write(luout,*) unita(3,3), control_unita(3,3)
            endif
        call errquit('Error reading wavefunctions', 1, CALC_ERR)
      end if
      end if

      if (.not.btdb_get(rtdb,'nwpw:psi_brillioun_check',
     >                  mt_log,1,brillioun_check))
     >   brillioun_check = .false.

      if (brillioun_check) then
      if ( nbrill.ne.nbrillioun) then
        call errquit(
     >   'Error reading wavefunctions - wrong number of k-points', 1,
     &       CALC_ERR)
      end if
      end if


*     *************************
*     *** read in 3d blocks ***
*     *************************

*     * file may contain more data then necessary *
      if (nbrillioun.le.nbrill) then

        do nb=1,nbrillioun
           call K1dB_ktoqp(nb,nbq,pk)
           do n=1,(ne(1)+ne(2))
              pj=0
              q=n
              call C3dB_c_Read(1,5,dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp(1)),pj,pk)

            if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
               call Cram_c_pack(nbq,dcpl_mb(tmp2(1)))
               call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp2(1)))
            end if
           end do
        end do

*     * file contains less data then necessary *
      else 

        do nb=1,nbrill-1
           call K1dB_ktoqp(nb,nbq,pk)
           do n=1,(ne(1)+ne(2))
              pj=0
              q=n

             call C3dB_c_Read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)),pj,pk)

             if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
                call Cram_c_pack(nbq,dcpl_mb(tmp2(1)))
                call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp2(1)))
             end if
        end do
        end do
*       **** use the last zone for filling in unknown zones ****
        do n=1,(ne(1)+ne(2))
          pj=0
          q=n
          call C3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)),pj,-1)
          do nb=nbrill,nbrillioun
             call K1dB_ktoqp(nb,nbq,pk)
             if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
                call dcopy(2*nfft3d,dcpl_mb(tmp2(1)),1,
     >                              dcpl_mb(tmp(1)),1)
                call Cram_c_pack(nbq,dcpl_mb(tmp(1)))
                call cpsi_data_add(psi2_tag,nbq,q,dcpl_mb(tmp(1)))
             end if
          end do
        end do

      end if


*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value =           BA_pop_stack(tmpe(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
*     end if

      return
      end

