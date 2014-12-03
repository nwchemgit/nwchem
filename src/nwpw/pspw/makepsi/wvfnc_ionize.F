*
* $Id$
*

*     *******************************************
*     *						*
*     *	 	   wvfnc_ionize			*
*     *						*
*     *******************************************

      logical function wvfnc_ionize(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

      logical value
      integer version
 
      integer ierr

      integer ne(2),ispin,nne(2)

      character*50 new_wavefunction_filename,filename
      character*50 old_wavefunction_filename
      character*255 full_filename,full_filename2


      integer ngrid(3)
      integer dngrid(3)
      integer cfull(2),dcfull(2)
      integer nfft3d,n2ft3d
      integer dnfft3d,dn2ft3d
      integer ms,n,l,ma_type,occupation


      double precision unita(3,3)
      integer norbs_up,up_orbital_list(2)
      integer norbs_dn,dn_orbital_list(2)

      logical  ionize_orbital_notfound
      external ionize_orbital_notfound

      value = .false.
      version = 3

*     **** get wavefunction information ****
      value = btdb_cget(rtdb,'ionize:old_wavefunction_filename',
     >                  1,old_wavefunction_filename)
      value = btdb_cget(rtdb,'ionize:new_wavefunction_filename',
     >                  1,new_wavefunction_filename)

      norbs_up = 0
      norbs_dn = 0
      if (rtdb_ma_get(rtdb, 'ionize:up_orbitals', ma_type,
     >                norbs_up, up_orbital_list(2))) then

         if (.not.BA_get_index(up_orbital_list(2),up_orbital_list(1)))
     >      call errquit(
     >       'wvfnc_ionize: ma_get_index failed for actlist',911,
     &       MA_ERR)
      end if
      if (rtdb_ma_get(rtdb, 'ionize:dn_orbitals', ma_type,
     >                norbs_dn, dn_orbital_list(2))) then

         if (.not.BA_get_index(dn_orbital_list(2),dn_orbital_list(1)))
     >      call errquit(
     >       'wvfnc_ionize: ma_get_index failed for actlist',911,
     &       MA_ERR)
      end if


 
      call util_file_name_noprefix(old_wavefunction_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename)

      l = index(full_filename,' ') - 1
      call openfile(5,full_filename,l,'r',l)
      call iread(5,version,1)
      call iread(5,ngrid,3)
      call dread(5,unita,9)
      call iread(5,ispin,1)
      call iread(5,ne,2)
      call iread(5,occupation,1)
 
      nne(1) = ne(1) - norbs_up
      nne(2) = ne(2) - norbs_dn
      
      l = index(new_wavefunction_filename,' ') - 1
      filename = new_wavefunction_filename(1:l)//".xpander"
      call util_file_name(filename,
     >                    .true.,
     >                    .false.,
     >                    full_filename)
      l = index(full_filename,' ') - 1
      call openfile(6,full_filename,l,'w',l)
      call iwrite(6,version,1)
      call iwrite(6,ngrid,3)
      call dwrite(6,unita,9)
      call iwrite(6,ispin,1)
      call iwrite(6,nne,2)
      call iwrite(6,occupation,1)

       
       nfft3d = ( ngrid(1)/2+1)* ngrid(2)* ngrid(3)
       n2ft3d = 2* nfft3d

      write(*,109) old_wavefunction_filename
      write(*,110) new_wavefunction_filename
      write(*,111) ne(1), ne(2),
     >             nne(1),nne(2)
  109 format(' old_filename: ',A)
  110 format(' new_filename: ',A)
  111 format(' converting  : ',I4,'x',I4,' --> ', 
     >                         I4,'x',I4)
 
*     ***** allocate wavefunction memory ****
      value = BA_alloc_get(mt_dcpl,nfft3d,
     >                     'cfull',cfull(2),cfull(1))


      do ms=1,ispin
        do n=1,ne(ms)
          call dread(5,dcpl_mb(cfull(1)),n2ft3d)

          if (ionize_orbital_notfound(ms,n,
     >                       int_mb(up_orbital_list(1)),norbs_up,
     >                       int_mb(dn_orbital_list(1)),norbs_dn)) then
            write(*,*) "copying .... psi:", n," spin:",ms
            call dwrite(6,dcpl_mb(cfull(1)),n2ft3d)
          else
            write(*,*) "removing .... psi:", n," spin:",ms
          end if

        end do
      end do
      call closefile(5)
      call closefile(6)

c     *** copy wvfnc_expander to new_wavefunction_filename ****
      call util_file_name_noprefix(new_wavefunction_filename,
     >                    .false.,
     >                    .false.,
     >                    full_filename2)
      call util_file_copy(full_filename,full_filename2)
      call util_file_unlink(full_filename)
      IERR=0
      GO TO 9999
 
 9110 IERR=10
      GO TO 9999
 9111 IERR=11
      GO TO 9999
 
 9999 value = BA_free_heap(cfull(2))
      if (norbs_up.gt.0) value = BA_free_heap(up_orbital_list(2))
      if (norbs_dn.gt.0) value = BA_free_heap(dn_orbital_list(2))
      !IF(IERR.EQ.0) THEN
      !  WRITE(6,*) ' JOB HAS BEEN COMPLETED.  CODE=',IERR
      !ELSE
      !  WRITE(6,*) ' JOB HAS BEEN TERMINATED DUE TO CODE=',IERR
      !  value = .false.
      !ENDIF
      !call nwpw_message(4)
      
      wvfnc_ionize = value
      return 
      end
 
*     *******************************************
*     *                                         *
*     *            ionize_orbital_notfound      *
*     *                                         *
*     *******************************************

      logical function ionize_orbital_notfound(ms,n,
     >                                         up_orbs,nup,
     >                                         dn_orbs,ndn)
      implicit none
      integer ms,n
      integer up_orbs(*),nup
      integer dn_orbs(*),ndn

      logical notfound
      integer i

      notfound = .true.

      if (ms.eq.1) then
         do i=1,nup
           if (up_orbs(i).eq.n) notfound = .false.
         end do
      end if
      if (ms.eq.2) then
         do i=1,ndn
           if (dn_orbs(i).eq.n) notfound = .false.
         end do
      end if

      ionize_orbital_notfound = notfound
      return
      end
