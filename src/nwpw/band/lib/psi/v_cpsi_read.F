*
* $Id$
*

*     ***************************
*     *          	   	*
*     *       v_cpsi_read	*
*     *             		*
*     ***************************
      subroutine v_cpsi_read(spin_orbit0,ispin0,ne0,nbrillioun0,
     >                       psi0_tag)
      implicit none 
      logical spin_orbit0
      integer ispin0,ne0(2),nbrillioun0
      integer psi0_tag

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
      logical hprint,readblock
      integer MASTER,taskid,taskid_j,taskid_k
      parameter(MASTER=0)
      logical brillioun_check,psi_nolattice
      integer n,nb,nbrill,occupation,ms
      integer msglen,psiold,psiold1
      integer q,pj,nbq,pk,i
      integer ispin,ne(2)
      logical spin_orbit

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2),tmpe(2)
      logical value

*     ***** local functions ****
      logical control_spin_orbit,control_print
      external control_spin_orbit,control_print
      character*50 control_input_v_psi
      external     control_input_v_psi
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

      value = BA_push_get(mt_dbl,(ne0(1)+ne0(2)),'tmpe',tmpe(2),tmpe(1))
       if (.not. value) call errquit('out of stack memory',0, MA_ERR)
     
*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_v_psi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                        full_filename)
         l = index(full_filename,' ') -1
         if (control_print(print_medium))
     >     write(luout,1210) full_filename(1:l)
 1210    FORMAT(/' input vpsi filename:',A)

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
      readblock = .true.

      if (spin_orbit) then
        if (.not.spin_orbit0) readblock = .false.
      end if
      if (.not.spin_orbit) then
        if (spin_orbit0) readblock = .false.
      end if

      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
         if (taskid.eq.MASTER) then
            write(*,*) "Error reading v_psi, v_psi set to zero"
            write(*,*) " - bad ngrid "
            write(*,*) " - vpsi_ngrid(1)=",nfft(1),
     >                 " ngrid(1)=",control_ngrid(1)
            write(*,*) " - vpsi_ngrid(2)=",nfft(2),
     >                 " ngrid(2)=",control_ngrid(2)
            write(*,*) " - vpsi_ngrid(3)=",nfft(3),
     >                 " ngrid(3)=",control_ngrid(3)
         end if
      end if

      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) .or.
     >     (ispin     .ne.ispin0)             .or.
     >     (nbrill    .ne.nbrillioun0)        .or.
     >     (ne(1)     .ne.ne0(1))             .or.
     >     (ne(2)     .ne.ne0(2))  ) then
         readblock = .false.
         if (taskid.eq.MASTER) then
            write(6,*) "Error reading v_psi, v_psi set to zero"
            write(6,*) " - incorrect unitcell,ispin, or ne"
         end if
      end if

*     *************************
*     *** read in 3d blocks ***
*     *************************

*     * file may contain more data then necessary *
      if (readblock) then
        do nb=1,nbrill
           call K1dB_ktoqp(nb,nbq,pk)
           do n=1,(ne(1)+ne(2))
              pj=0
              q=n
              call C3dB_c_Read(1,5,dcpl_mb(tmp2(1)),
     >                             dcpl_mb(tmp(1)),pj,pk)

            if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
               call Cram_c_pack(nbq,dcpl_mb(tmp2(1)))
               call cpsi_data_add(psi0_tag,nbq,q,dcpl_mb(tmp2(1)))
            end if
           end do
        end do
        if (occupation.gt.0) then
           l = cpsi_data_get_next(psi0_tag)
           if (nbrillioun0.eq.nbrill) then
  
              do nb=1,nbrill
                 call K1dB_ktoqp(nb,nbq,pk)
                 if (taskid.eq.MASTER) then
                    call dread(5,dbl_mb(tmpe(1)),(ne(1)+ne(2)))
                 end if
                 !** big broadcast - should be reduced ***
                 call Parallel_Brdcst_values(MASTER,
     >                         (ne(1)+ne(2)),dbl_mb(tmpe(1)))
                 if (pk.eq.taskid_k) then
                    do n=1,(ne(1)+ne(2))
                       call cpsi_data_add(l,nbq,n,dbl_mb(tmpe(1)+n-1))
                    end do
                 end if
              end do
*          **** reset to zeros and ones ****
           else
             do nb=1,nbrillioun0
               call K1dB_ktoqp(nb,nbq,pk) 
               if (pk.eq.taskid_k) then 
                 do ms=1,ispin
                 do n =1,ne(ms)
                 if (n.le.(ne(ms)-control_fractional_orbitals(ms))) then
                   unita(1,1) = 1.0d0  
                 else
                   unita(1,1) = 0.0d0  
                 end if
                 call cpsi_data_add(l,nbq,(n+(ms-1)*ne(1)),unita)
                 end do
                 end do
               end if
             end do
           end if

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

