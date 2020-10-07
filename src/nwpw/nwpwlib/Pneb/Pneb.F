*
* $Id$
*
#define NBLOCKS 4


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C        Pneb - (P)arallel (n) (e)lectron with (b)rillioun zone
C
C            (P)atrick J. (N)ichols and (E)ric J. (B)ylaska


*     The Pne (parallel (n) electrons and brillioun zone) library is to
* be used for handling three kinds of data structures.  The first
* data structure, denoted by "m", represents a set of double precision matrices of
* [nexne,1:nbrillion].  The second data structure, denoted by "w" represents a set of 
* double complex matrix of size [nexne, 1:nbrillioun].  The third data structure, 
* denoted by "f" represents a double complex array of size [1:nwave,ne,1:nbrillioun].
*
* Where 
*    ne = { 1:ne(1)                     if ispin==1                        || 
*           1:ne(1),ne(1)+1:ne(1)+ne(2) if ispin==2 && spin_orbit==.false. ||
*           1:ne(1),ne(1)+1:ne(1)+ne(1) if ispin==2 && spin_orbit==.true.
*         }
*
* nexne={[1:ne(1),1:ne(1)]                                   if ispin==1                       ||
*        [1:ne(1),1:ne(1)],[1:ne(2),1:ne(2)]                 if ispin==2 && spin_orbit==.false.||
*        [1:ne(1),1:ne(1)],[1:ne(2),1:ne(2)] && ne(2)==ne(1) if ispin==2 && spin_orbit==.false.||
*       }        

* When spin_orbit==.true. the up and down compnents of the f array are assumed to be stored 
* as follows.
*
*  The up components of the ne(1) electrons are stored in
*
*      f[1:nwave,1:ne(1),1:nbrillion]
*
*  and the down components of the ne(1) electrons are stored in
*
*      f[1:nwave,ne(1)+1:ne(1)+ne(1),1:nbrillioun]
*
* When spin_orbit==.false. && ispin==1 the restricted components of the f array
* are stored as follows
*
* The restricted (up and down) ne(1) electrons are stored in
*
*      f[1:nwave,1:ne(1),1:nbrillioun]
*
* When spin_orbit==.false. && ispin==2 the up and down compnents of the f array are stored 
* as follows.
*
*  The ne(1) up electrons are stored in
*
*      f[1:nwave,1:ne(1),1:nbrillion]
*
*  and the ne(2) down electrons stored in
*
*      f[1:nwave,ne(1)+1:ne(1)+ne(2),1:nbrillioun]



* Implementation Details:
*
* All the operations on Pneb objects are controlled by the mb and nb parameters.
* 
* mb can have values of 0, 1 or 2.
*
*   mb==0 means that the 'm' array contains the data for both the up and down matrices.
*
*   mb==1 means that the 'm' array only contains the data for the up matrices.
*   
*   mb==2 means that the 'm' array only contains the data for the down matrices. 

c
c if mb==0 means m array contains all data for both up and down matrices
c if mb==1 means the

c if (mb==0 && nb==0)


* Note when spin_orbit==.true. Pneb has the following properties  
*
*     - the value for ispin is always equal to 2.
*
*     - The matrices 'm' and 'w' can be represented in block diagonal form 
*
*       m = { mup  0;          w = { wup  0;
*             0    mdown}            0    wdown}
*
*         where
*              m[i,j]== mup[i,j] == mdown[i,j]
*              w[i,j]== wup[i,j] == wdown[i,j]
               
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



C     ***********************************
c     *                                 *
c     *          Pneb_is_init           *
c     *                                 *
c     ***********************************

      logical function Pneb_is_init()
      implicit none
#include "Pneb.fh"
      Pneb_is_init=im_init
      return
      end

c     ***********************************
c      Reset the number of eigenvalues
c      to a new value. 
c     ***********************************
      subroutine Pneb_set_ne(ne_in)
      implicit none
#include "Pneb.fh"
      integer ne_in(2)
      if (spin_orbit) then
         if (ne_in(1).ne.ne_in(2)) then
           write(*,*)"WARNING EXCITED NE(1) NOT EQUAL EXCITED NE(2)"
           write(*,*)"AND SPIN ORBIT EQUAL TRUE!"
           write(*,*)"USING MAX EXCITED NE!"
           call flush(6)
           if (ne_in(1).gt.ne_in(2)) then
             ne_in(2)=ne_in(1)
           else
             ne_in(1)=ne_in(2)
           end if
         end if          
      end if 
      ne(1)=ne_in(1)
      ne(2)=ne_in(2)
cccccccccccccccccccccccccc
cccccc temporary fix !
cccccc fix me!!!!!!!
cccccccccccccccccccccccccc
      neq(1)=ne_in(1)
      neq(2)=ne_in(2)
      return
      end

*     ***********************************
*     *					*
*     *	        Pneb_init		*	
*     *					*
*     ***********************************

      subroutine Pneb_init(ispin_in,ne_in,nbrill_in,
     >                     map_in,spin_orbit_in)
      implicit none
      integer ispin_in
      integer ne_in(2)
      integer nbrill_in
      integer map_in
      logical spin_orbit_in

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

*     **** local variables ****
      logical value
      integer ms,k,i,j,npack1,n2ft3d,nthr
      integer nework,ncqmax0
      integer ii,jj,icur,jcur,indx0,indx1,indx2,nn


*     **** external functions ****
      logical  control_mparallelized
      external control_mparallelized
      integer  brillioun_nbrillq
      external brillioun_nbrillq
      integer  Parallel3d_comm_i,Parallel3d_comm_j,Parallel3d_comm_k
      external Parallel3d_comm_i,Parallel3d_comm_j,Parallel3d_comm_k
      integer  control_mreplicate_size
      external control_mreplicate_size


      ispin   = ispin_in
      ne(1)   = ne_in(1)
      ne(2)   = ne_in(2)
      nbrill  = nbrill_in
      nbrillq =  brillioun_nbrillq()
      spin_orbit = spin_orbit_in

      neq(1)  = 0
      neq(2)  = 0
      call C1dB_init(1,ne(1),map_in)
      call C1dB_nq(1,neq(1))
      if ((ispin.eq.2).and.(ne(2).gt.0)) then
         call C1dB_init(2,ne(2),map_in)
         call C1dB_nq(2,neq(2))
      end if

      call Parallel3d_np_i(np_i)
      call Parallel3d_np_j(np_j)
      call Parallel3d_np_k(np_k)
      comm_i = Parallel3d_comm_i()
      comm_j = Parallel3d_comm_j()
      comm_k = Parallel3d_comm_k()
      parallelized_j = (np_j.gt.1)
      parallelized_k = (np_k.gt.1)

      mparallelized = (parallelized_j.and.control_mparallelized())
      mreplicate = ((np_i*np_j).gt.control_mreplicate_size())


      if (spin_orbit) then
         ispinq = 1
      else
         ispinq = ispin
      end if

      if (parallelized_j) then
         value = .true.
         mcq(1) = 0
         mcq(2) = 0
         ncq(1) = 0
         ncq(2) = 0
         mcqmax(1) = 0
         mcqmax(2) = 0
         ncqmax(1) = 0
         ncqmax(2) = 0
         ncqmax0 = 0
         do ms=1,ispin
             value = value.and.
     >               BA_alloc_get(mt_int,np_i,'ma',ma(2,ms),ma(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_i,'ma2',ma2(2,ms),ma2(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_i,'mc',mc(2,ms),mc(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_j,'na',na(2,ms),na(1,ms))
             value = value.and.
     >               BA_alloc_get(mt_int,np_j,'nc',nc(2,ms),nc(1,ms))
             if (.not.value) then
               call errquit('Pneb_init: out of heap memory',0,MA_ERR)
             end if

             call icopy(np_i,0,0,int_mb(ma(1,ms)),1)
             call icopy(np_i,0,0,int_mb(ma2(1,ms)),1)
             call icopy(np_i,0,0,int_mb(mc(1,ms)),1)
             call icopy(np_j,0,0,int_mb(na(1,ms)),1)
             call icopy(np_j,0,0,int_mb(nc(1,ms)),1)

             i = 0
             j = 0
             do k=1,ne(ms)
                int_mb(mc(1,ms)+i) = int_mb(mc(1,ms)+i) + 1

                int_mb(nc(1,ms)+j) = int_mb(nc(1,ms)+j) + 1
                int_mb(na(1,ms)+j) = int_mb(na(1,ms)+j) + 1
                i = mod(i+1,np_i)
                j = mod(j+1,np_j)
             end do

             call Cram_max_npack(npack1)
             call C3dB_n2ft3d(1,n2ft3d)
             int_mb(ma(1,ms) +taskid_i) = 2*npack1
             int_mb(ma2(1,ms)+taskid_i) = n2ft3d
             call C3dB_Vector_iSumAll(np_i,int_mb(ma(1,ms)))
             call C3dB_Vector_iSumAll(np_i,int_mb(ma2(1,ms)))

             mcq(ms) = int_mb(mc(1,ms)+taskid_i)
             ncq(ms) = int_mb(nc(1,ms)+taskid_j)
             mcqmax(ms) = 0
             do i=0,np_i-1
                if (int_mb(mc(1,ms)+i).gt.mcqmax(ms))
     >             mcqmax(ms) = int_mb(mc(1,ms)+i)
             end do
             ncqmax(ms) = 0
             do j=0,np_j-1
                if (int_mb(nc(1,ms)+j).gt.ncqmax(ms))
     >             ncqmax(ms) = int_mb(nc(1,ms)+j)
             end do

             ncqmax0 = 0
             do j=0,np_j-1
                if (int_mb(nc(1,1)+j).gt.ncqmax0)
     >             ncqmax0 = int_mb(nc(1,1)+j)
             end do

             npack1_all = 0
             n2ft3d_all = 0
             do i=0,np_i-1
                npack1_all = npack1_all + int_mb(ma(1,ms) +i)
                n2ft3d_all = n2ft3d_all + int_mb(ma2(1,ms)+i)
             end do

         end do

         mdq = 0
         ndq = 0
         if ((ispin.eq.2).and.(ne(2).gt.0)) then
            value = BA_alloc_get(mt_int,np_i,'md',md(2),md(1))
            value = value.and.
     >              BA_alloc_get(mt_int,np_j,'nd',nd(2),nd(1))
            if (.not.value) then
               call errquit('Dne_init: out of heap memory',1,MA_ERR)
            end if
            call icopy(np_i,0,0,int_mb(md(1)),1)
            call icopy(np_j,0,0,int_mb(nd(1)),1)

            i = 0
            do k=1,ne(1)
               int_mb(md(1)+i) = int_mb(md(1)+i) + 1
               i = mod(i+1,np_i)
            end do

            j = 0
            do k=1,ne(2)
               int_mb(nd(1)+j) = int_mb(nd(1)+j) + 1
               j = mod(j+1,np_j)
            end do
            mdq = int_mb(md(1)+taskid_i)
            ndq = int_mb(nd(1)+taskid_j)
         end if


         value = value.and.BA_alloc_get(mt_dcpl,
     >                                 2*2*64*int_mb(ma2(1,1)+taskid_i),
     >                                  'work1',work1(2),work1(1))

         nework = 2*64*int_mb(nc(1,1)+taskid_j)
         if (nework.lt.mcq(1)*ncq(1)) nework = mcq(1)*ncq(1)
         value = value.and.BA_alloc_get(mt_dcpl,
     >                                  3*nework,
     >                                  'work2',work2(2),work2(1))

         value = value.and.BA_alloc_get(mt_dcpl,
     >                              ne(1)*ncqmax0,
     >                              'bcolwork',bcolwork(2),bcolwork(1))
         value = value.and.BA_alloc_get(mt_dcpl,
     >                              ne(1)*ncqmax0,
     >                              'bwork2',bwork2(2),bwork2(1))
         nework = int_mb(na(1,1))
         value = value.and.BA_alloc_get(mt_dcpl,
     >                              2*nework*int_mb(ma2(1,1)+taskid_i),
     >                              'rwork1',rwork1(2),rwork1(1))
         value = value.and.BA_alloc_get(mt_dcpl,
     >                              2*nework*int_mb(ma2(1,1)+taskid_i),
     >                              'rwork2',rwork2(2),rwork2(1))
         nework = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         value = value.and.BA_alloc_get(mt_dcpl,3*nework,
     >                              'mat_tmp',mat_tmp(2),mat_tmp(1))

         nework = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'mindx0',mindx(2,0),mindx(1,0))
         nework = mcq(1)*ncq(1)
         value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'mindx1',mindx(2,1),mindx(1,1))
         if (ispin.eq.2) then
         nework = mcq(2)*ncq(2)
         value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'mindx2',mindx(2,2),mindx(1,2))
         endif
         if (.not.value)
     >      call errquit('Pneb_init: out of heap memory',1,MA_ERR)

         mall(0)  = ne(1)*ne(1) + ne(2)*ne(2)
         mall(1)  = ne(1)*ne(1)
         mall(2)  = ne(2)*ne(2)
         mpack(0) = mcq(1)*ncq(1) + mcq(2)*ncq(2)
         mpack(1) = mcq(1)*ncq(1)
         mpack(2) = mcq(2)*ncq(2)
         indx0=0
         indx1=0
         indx2=0
         jj   = 1
         jcur = 0
         do j=1,ne(1)
            ii   = 1
            icur = 0
            do i=1,ne(1)
               if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  int_mb(mindx(1,0)+indx0) = i+(j-1)*ne(1)
                  int_mb(mindx(1,1)+indx1) = i+(j-1)*ne(1)
                  indx0 = indx0 + 1
                  indx1 = indx1 + 1
               end if
               ii = ii+1
               if (ii.gt.int_mb(mc(1,1)+icur)) then
                  icur = icur + 1
                  ii   = 1
               end if
            end do
            jj = jj+1
            if (jj.gt.int_mb(nc(1,1)+jcur)) then
               jcur = jcur + 1
               jj   = 1
            end if
         end do
         if (ispin.eq.2) then
         jj   = 1
         jcur = 0
         do j=1,ne(2)
            ii   = 1
            icur = 0
            do i=1,ne(2)
               if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  int_mb(mindx(1,0)+indx0) = i+(j-1)*ne(2) + ne(1)*ne(1)
                  int_mb(mindx(1,2)+indx2) = i+(j-1)*ne(2)
                  indx0 = indx0 + 1
                  indx2 = indx2 + 1
               end if
               ii = ii+1
               if (ii.gt.int_mb(mc(1,2)+icur)) then
                  icur = icur + 1
                  ii   = 1
               end if
            end do
            jj = jj+1
            if (jj.gt.int_mb(nc(1,2)+jcur)) then
               jcur = jcur + 1
               jj   = 1
            end if
         end do
         end if

         if ((ispin.eq.2).and.(ne(2).gt.0)) then
            nework = mdq*ndq
            value = value.and.BA_alloc_get(mt_int,nework,
     >                                 'smindx',smindx(2),smindx(1))
            small  = ne(1)*ne(2)
            smpack = mdq*ndq
            indx2=0
            jj   = 1
            jcur = 0
            do j=1,ne(2)
               ii   = 1
               icur = 0
               do i=1,ne(1)
               if ((icur.eq.taskid_i).and.(jcur.eq.taskid_j)) then
                  int_mb(smindx(1)+indx2) = i+(j-1)*ne(1)
                  indx2 = indx2 + 1
               end if
               ii = ii+1
               if (ii.gt.int_mb(md(1)+icur)) then
                  icur = icur + 1
                  ii   = 1
               end if
               end do
               jj = jj+1
               if (jj.gt.int_mb(nd(1)+jcur)) then
                  jcur = jcur + 1
                  jj   = 1
               end if
            end do
         end if

      end if


      im_init= .true.
      return
      end


*     ***********************************
*     *					*
*     *	          Pneb_end   		*	
*     *					*
*     ***********************************

      subroutine Pneb_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

*     ***** local variables ****
      logical value
      integer ms

      call C1dB_end(1)
      if ((ispin.eq.2).and.(ne(2).gt.0)) call C1dB_end(2)

      if (parallelized_j) then
         value = .true.
         do ms=1,ispin
            value = value.and.BA_free_heap(ma(2,ms))
            value = value.and.BA_free_heap(ma2(2,ms))
            value = value.and.BA_free_heap(na(2,ms))
            value = value.and.BA_free_heap(mc(2,ms))
            value = value.and.BA_free_heap(nc(2,ms))
         end do
         if ((ispin.eq.2).and.(ne(2).gt.0)) then
            value = value.and.BA_free_heap(md(2))
            value = value.and.BA_free_heap(nd(2))
            value = value.and.BA_free_heap(smindx(2))
         end if
         value = value.and.BA_free_heap(work1(2))
         value = value.and.BA_free_heap(work2(2))
         value = value.and.BA_free_heap(bcolwork(2))
         value = value.and.BA_free_heap(bwork2(2))
         value = value.and.BA_free_heap(rwork1(2))
         value = value.and.BA_free_heap(rwork2(2))
         value = value.and.BA_free_heap(mat_tmp(2))
         value = value.and.BA_free_heap(mindx(2,0))
         value = value.and.BA_free_heap(mindx(2,1))
         if (ispin.eq.2) value = value.and.BA_free_heap(mindx(2,2))
         if (.not.value) then
           call errquit('Pneb_end:error deallocating heap',1,MA_ERR)
         end if
      end if

      im_init=.false.
      return
      end
         


*     ***********************************
*     *					*
*     *	        Pneb_neq		*	
*     *					*
*     ***********************************

      subroutine Pneb_neq(nqtmp)
      implicit none
      integer nqtmp(2)

#include "Pneb.fh"

      nqtmp(1) = neq(1)
      nqtmp(2) = neq(2)

      return
      end


*     ***********************************
*     *					*
*     *	        Pneb_nbrill 		*	
*     *					*
*     ***********************************

      integer function Pneb_nbrill()
      implicit none

#include "Pneb.fh"

      Pneb_nbrill = nbrill
      return 
      end


*     ***********************************
*     *					*
*     *	        Pneb_ne 		*	
*     *					*
*     ***********************************

      subroutine Pneb_ne(netmp)
      implicit none
      integer netmp(2)

#include "Pneb.fh"

      netmp(1) = ne(1)
      netmp(2) = ne(2)
      return 
      end



*     ***********************************
*     *					*
*     *	        Pneb_nbrillq 		*	
*     *					*
*     ***********************************

      integer function Pneb_nbrillq()
      implicit none

#include "Pneb.fh"

      Pneb_nbrillq = nbrillq
      return 
      end


*     ***********************************
*     *					*
*     *	        Pneb_spin_orbit 	*	
*     *					*
*     ***********************************

      logical function Pneb_spin_orbit()
      implicit none

#include "Pneb.fh"

      Pneb_spin_orbit = spin_orbit
      return 
      end

*     ****************************************
*     *                                      *
*     *            Pneb_ispin                *
*     *                                      *
*     ****************************************

      integer function Pneb_ispin()
      implicit none

#include "Pneb.fh"

      Pneb_ispin = ispin
      return 
      end

*     ****************************************
*     *                                      *
*     *            Pneb_convert_nb           *
*     *                                      *
*     ****************************************
      integer function Pneb_convert_nb(nbq)
      implicit none
      integer nbq,nb

      call K1dB_qtok(nbq,nb)

      Pneb_convert_nb = nb
      return
      end 


*     ****************************************
*     *                                      *
*     *           Pneb_ispinq                *
*     *                                      *
*     ****************************************

      integer function Pneb_ispinq()
      implicit none

#include "Pneb.fh"

      Pneb_ispinq = ispinq
      return
      end


c     ****************************************
c     *                                      *
c     *        Pneb_m_size                   *
c     *                                      *
c     ****************************************

c   This function returns the size of the 'm' array
c with options mb and nb.
c
c if mb==0 means m array contains all data for both up and down matrices
c if mb==1 means the
c if (mb==0 && nb==0)

      integer function Pneb_m_size(mb,nb)
      implicit none
      integer mb,nb

#include "Pneb.fh"
     
      integer mbsize,nbsize

      if (nb.eq.0) then
        nbsize = nbrillq
      else
        nbsize = 1
      end if
 
      if (mb.eq.0) then
         if (.not.spin_orbit) then 
            mbsize = ne(1)*ne(1) + ne(2)*ne(2)       !*** Changing back to a more packed space! ***
         else
            mbsize = ne(1)*ne(1)
         end if
      else
         mbsize = ne(mb)*ne(mb)
      end if
     
      Pneb_m_size = mbsize*nbsize
      return
      end



c
c     ****************************************
c     *                                      *
c     *        Pneb_w_size                   *
c     *                                      *
c     ****************************************
      integer function Pneb_w_size(mb,nb)
      implicit none
      integer mb,nb

#include "Pneb.fh"
     
      integer mbsize,nbsize

      if (nb.eq.0) then
        nbsize = nbrillq
      else
        nbsize = 1
      end if
 
      if (mb.eq.0) then
         if (.not.spin_orbit) then 
            mbsize = ne(1)*ne(1) + ne(2)*ne(2)       !*** Changing back to a more packed space! ***
         else
            mbsize = ne(1)*ne(1)
         end if
      else
         mbsize = ne(mb)*ne(mb)
      end if
     
      Pneb_w_size = mbsize*nbsize
      return
      end



c
c     ****************************************
c     *                                      *
c     *        Pneb_wne_size                 *
c     *                                      *
c     ****************************************
      integer function Pneb_wne_size(mb,nb,nein)
      implicit none
      integer mb,nb,nein(2)

#include "Pneb.fh"

      integer mbsize,nbsize

      if (nb.eq.0) then
        nbsize = nbrillq
      else
        nbsize = 1
      end if

      if (mb.eq.0) then
         if (.not.spin_orbit) then
            mbsize = nein(1)*nein(1) + nein(2)*nein(2)       !*** Changing back to a more packed space! ***
         else
            mbsize = nein(1)*nein(1)
         end if
      else
         mbsize = nein(mb)*nein(mb)
      end if

      Pneb_wne_size = mbsize*nbsize
      return
      end



c     ****************************************
c     *                                      *
c     *        Pneb_m_allocate               *
c     *                                      *
c     ****************************************

      logical function Pneb_m_allocate(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "bafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size
     
      size = Pneb_m_size(mb,nb)
      Pneb_m_allocate = BA_alloc_get(mt_dbl,size,'hml',hml(2),hml(1))

      return
      end

c     ****************************************
c     *                                      *
c     *        Pneb_w_allocate               *
c     *                                      *
c     ****************************************

      logical function Pneb_w_allocate(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "bafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size
     
      size = Pneb_w_size(mb,nb)
      Pneb_w_allocate = BA_alloc_get(mt_dcpl,size,'hml',hml(2),hml(1))
      return
      end



c     ****************************************
c     *                                      *
c     *        Pneb_wne_allocate_block       *
c     *                                      *
c     ****************************************

      logical function Pneb_wne_allocate_block(mb,nb,nblocks,nein,hml)
      implicit none
      integer mb,nb,nblocks,nein(2)
      integer hml(2)

#include "bafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_wne_size
      external Pneb_wne_size

      size = nblocks*Pneb_wne_size(mb,nb,nein)
      Pneb_wne_allocate_block 
     > = BA_alloc_get(mt_dcpl,size,'hml',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *        Pneb_m_allocate_block         *
c     *                                      *
c     ****************************************

      logical function Pneb_m_allocate_block(mb,nb,nblocks,hml)
      implicit none
      integer mb,nb,nblocks
      integer hml(2)

#include "bafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size
     
      size = nblocks*Pneb_m_size(mb,nb)
     
      Pneb_m_allocate_block =
     >   BA_alloc_get(mt_dbl,size,'hml',hml(2),hml(1))
      return
      end



c     ****************************************
c     *                                      *
c     *        Pneb_w_allocate_block         *
c     *                                      *
c     ****************************************

      logical function Pneb_w_allocate_block(mb,nb,nblocks,hml)
      implicit none
      integer mb,nb,nblocks
      integer hml(2)

#include "bafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size
     
      size = nblocks*Pneb_w_size(mb,nb)
     
      Pneb_w_allocate_block =
     >   BA_alloc_get(mt_dcpl,size,'hml',hml(2),hml(1))
      return
      end


c
c     ****************************************
c     *                                      *
c     *          Pneb_m_free                 *
c     *                                      *
c     ****************************************
      logical function Pneb_m_free(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Pneb_m_free = BA_free_heap(hml(2))
      return
      end

c     ****************************************
c     *                                      *
c     *          Pneb_w_free                 *
c     *                                      *
c     ****************************************
      logical function Pneb_w_free(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Pneb_w_free = BA_free_heap(hml(2))
      return
      end



c     ****************************************
c     *                                      *
c     *          Pneb_m_push_get             *
c     *                                      *
c     ****************************************
      logical function Pneb_m_push_get(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "Pneb.fh"
#include "bafdecls.fh"

*     **** local variables
      integer size

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size

      size = Pneb_m_size(mb,nb)
      Pneb_m_push_get = BA_push_get(mt_dbl,size,'hml',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *          Pneb_w_push_get             *
c     *                                      *
c     ****************************************
      logical function Pneb_w_push_get(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "Pneb.fh"
#include "bafdecls.fh"

*     **** local variables
      integer size

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      size = Pneb_w_size(mb,nb)
      Pneb_w_push_get = BA_push_get(mt_dcpl,size,'hml',hml(2),hml(1))
      return
      end




c     ****************************************
c     *                                      *
c     *          Pneb_m_push_get_block       *
c     *                                      *
c     ****************************************
      logical function Pneb_m_push_get_block(mb,nb,nblock,hml)
      implicit none
      integer mb,nb,nblock
      integer hml(2)

#include "bafdecls.fh"
#include "Pneb.fh"

*     **** local variables
      integer nsize

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size

      nsize = nblock*Pneb_m_size(mb,nb)
      Pneb_m_push_get_block = 
     >    BA_push_get(mt_dbl,nsize,'hml',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *          Pneb_w_push_get_block       *
c     *                                      *
c     ****************************************
      logical function Pneb_w_push_get_block(mb,nb,nblock,hml)
      implicit none
      integer mb,nb,nblock
      integer hml(2)

#include "bafdecls.fh"
#include "Pneb.fh"

*     **** local variables
      integer nsize

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      nsize = nblock*Pneb_w_size(mb,nb)
      Pneb_w_push_get_block = 
     >    BA_push_get(mt_dcpl,nsize,'hml',hml(2),hml(1))
      return
      end



c     ****************************************
c     *                                      *
c     *          Pneb_m_pop_stack            *
c     *                                      *
c     ****************************************
      logical function Pneb_m_pop_stack(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Pneb_m_pop_stack = BA_pop_stack(hml(2))
      return
      end


c     ****************************************
c     *                                      *
c     *          Pneb_w_pop_stack            *
c     *                                      *
c     ****************************************
      logical function Pneb_w_pop_stack(hml)
      implicit none
      integer hml(2)

#include "bafdecls.fh"

      Pneb_w_pop_stack = BA_pop_stack(hml(2))
      return
      end


c     ****************************************
c     *                                      *
c     *             Pneb_w_zero              *
c     *                                      *
c     ****************************************
c      Zeros out a complex matrix, w.

      subroutine Pneb_w_zero(mb,nbb,w)
      implicit none
      integer mb,nbb
      complex*16 w(*)

#include "Pneb.fh"

c     *** local variables ***
      integer nsize

c     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      nsize = Pneb_w_size(mb,nbb)
      call dcopy(2*nsize,0.0d0,0,w,1)
      return
      end



c     ****************************************
c     *                                      *
c     *             Pneb_m_zero              *
c     *                                      *
c     ****************************************
c      Zeros out a real matrix, m.

      subroutine Pneb_m_zero(mb,nbb,m)
      implicit none
      integer mb,nbb
      real*8 m(*)

#include "Pneb.fh"

c     *** local variables ***
      integer nsize

c     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size

      nsize = Pneb_m_size(mb,nbb)
      call dcopy(nsize,0.0d0,0,m,1)

      return
      end




c     ****************************************
c     *                                      *
c     *          Pneb_ffm_Multiply           *
c     *                                      *
c     ****************************************

      subroutine Pneb_ffm_Multiply(mb,nbb,Y1,Y2,npack1,hml)
      implicit none
      integer    mb,nbb
      complex*16 Y1(*),Y2(*)
      integer    npack1
      real*8     hml(*)

#include "bafdecls.fh"   
#include "errquit.fh"
#include "Pneb.fh"

*     **** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,nsize,n,np,nn
      integer shift,shiftso,shift2,ishift2,ishift3
      integer tid,nthr

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call Parallel3d_np_i(np)

      if (parallelized_j) then
         nsize = Pneb_m_size(mb,nbb)
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
            !mm = mcq(1)*ncq(1) + mcq(2)*ncq(2)
            nn = mcqmax(1)*ncqmax(1)+mcqmax(2)*ncqmax(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            !mm = mcq(mb)*ncq(mb) 
            nn = mcqmax(mb)*ncqmax(mb)
         end if

         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if

         do nb=nb1,nb2
            if (mparallelized) then
               do ms=ms1,ms2
                  if (ne(ms).eq.0) go to 20
                  shift  = 1 + (ms-1)*neq(1)*npack1
     >                   + ((nb-nb1))*(neq(1)+neq(2))*npack1
                  shiftso = shift + neq(1)*npack1
                  shift2  = 1 + (ms-1)*ishift2 + (nb-nb1)*ishift3

                  call DMatrix_dgemm2(ne(ms),ne(ms),npack1_all,128,
     >                1.0d0,
     >                Y1(shift),int_mb(ma(1,ms)+taskid_i),
     >                          int_mb(ma(1,ms)),
     >                          int_mb(na(1,ms)),
     >                Y2(shift),int_mb(ma(1,ms)+taskid_i),
     >                          int_mb(ma(1,ms)),
     >                          int_mb(na(1,ms)),
     >                0.0d0,
     >                hml(shift2),int_mb(mc(1,ms)+taskid_i),
     >                            int_mb(mc(1,ms)),
     >                            int_mb(nc(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i, comm_j,
     >                dcpl_mb(work1(1)),dcpl_mb(work2(1)))

                  if (spin_orbit) then
                     call DMatrix_dgemm2(ne(ms),ne(ms),npack1_all,128,
     >                   1.0d0,
     >                   Y1(shiftso),int_mb(ma(1,ms)+taskid_i),
     >                             int_mb(ma(1,ms)),
     >                             int_mb(na(1,ms)),
     >                   Y2(shiftso),int_mb(ma(1,ms)+taskid_i),
     >                             int_mb(ma(1,ms)),
     >                             int_mb(na(1,ms)),
     >                   1.0d0,
     >                   hml(shift2),int_mb(mc(1,ms)+taskid_i),
     >                               int_mb(mc(1,ms)),
     >                               int_mb(nc(1,ms)),
     >                   taskid_i,taskid_j,
     >                   np_i,np_j,
     >                   comm_i, comm_j,
     >                   dcpl_mb(work1(1)),dcpl_mb(work2(1)))
                  end if
  20              continue
               end do
            else
               do ms=ms1,ms2
                  if (ne(ms).eq.0) go to 21
                  shift  = 1 + (ms-1)*neq(1)*npack1
     >                   + ((nb-nb1))*(neq(1)+neq(2))*npack1
                  shiftso = shift + neq(1)*npack1
                  shift2  = 1 + (ms-1)*ishift2 + (nb-nb1)*ishift3
                  call CMatrix_dgemm2c_omp(ne(ms),ne(ms),npack1_all,128,
     >                 Y1(shift),Y2(shift),int_mb(ma(1,ms)+taskid_i),
     >                                     int_mb(ma(1,ms)),
     >                                     int_mb(na(1,ms)),
     >                 dcpl_mb(mat_tmp(1)+shift2),
     >                                     int_mb(mc(1,ms)+taskid_i),
     >                                     int_mb(mc(1,ms)),
     >                                     int_mb(nc(1,ms)),
     >                 taskid_i,taskid_j,
     >                 np_i,np_j,
     >                 comm_i, comm_j,
     >                 tid,nthr,
     >                 dcpl_mb(work1(1)),dcpl_mb(work2(1)),
     >                 dcpl_mb(thrwork1(1)+tid*nn))
                  if (spin_orbit) then
                     call errquit(
     >                     "Pneb_ffm_Multiply: spin orbit not finished",
     >                      0,0)
                  end if
  21              continue
               end do
               call Pneb_m_gather(mall(mb),mpack(mb),
     >                            int_mb(mindx(1,mb)),
     >                            dcpl_mb(mat_tmp(1)),hml)

            end if

         end do


      else
         nsize = Pneb_m_size(mb,nbb)
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispinq
            ishift2 = ne(1)*ne(1)  !** note that ne rather than neq is used **
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
      
         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if

         do nb=nb1,nb2
         do ms=ms1,ms2
           shift  = 1 + (ms-1)*ne(1)*npack1
     >            + ((nb-nb1))*(ne(1)+ne(2))*npack1
           shiftso = shift + ne(1)*npack1
           shift2 = 1 + (ms-1)*ishift2      +(nb-nb1)*ishift3
           n = ne(ms)
           if (n.eq.0) go to 30

              call Cram_ccm_idgemm(nb,n,n,Y1(shift),Y2(shift),
     >                            1.0d0,0.0d0,hml(shift2))

              if (spin_orbit) then
                 call Cram_ccm_idgemm(nb,n,n,Y1(shiftso),Y2(shiftso),
     >                               1.0d0,1.0d0,hml(shift2))
              end if

 30        continue
         end do
         end do
         if (np.gt.1) call C3dB_Vector_SumAll(nsize,hml)
      end if

      return
      end



c     ****************************************
c     *                                      *
c     *          Pneb_ffw_Multiply           *
c     *                                      *
c     ****************************************

      subroutine Pneb_ffw_Multiply(mb,nbb,Y1,Y2,npack1,hml)
      implicit none
      integer    mb,nbb
      complex*16 Y1(*),Y2(*)
      integer    npack1
      complex*16 hml(*)

#include "bafdecls.fh"   
#include "errquit.fh"
#include "Pneb.fh"

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0))

      integer nb,nb1,nb2,ms,ms1,ms2,nsize,n,np,nn
      integer shift,shiftso,shift2,ishift2,ishift3
      integer tid,nthr

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()


      call Parallel3d_np_i(np)
      nsize = Pneb_w_size(mb,nbb)

      if (parallelized_j) then
         nsize = Pneb_w_size(mb,nbb)
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
            !mm = mcq(1)*ncq(1) + mcq(2)*ncq(2)
            nn = mcqmax(1)*ncqmax(1)+mcqmax(2)*ncqmax(2)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            !mm = mcq(mb)*ncq(mb) 
            nn = mcqmax(mb)*ncqmax(mb)
         end if

         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if
         do nb=nb1,nb2
            if (mparallelized) then
               do ms=ms1,ms2
                  if (ne(ms).eq.0) go to 20
                  shift  = 1 + (ms-1)*neq(1)*npack1
     >                   + ((nb-nb1))*(neq(1)+neq(2))*npack1
                  shiftso = shift + neq(1)*npack1
                  shift2  = 1 + (ms-1)*ishift2 + (nb-nb1)*ishift3

                  call CMatrix_zgemm2(ne(ms),ne(ms),npack1_all,128,
     >                1.0d0,
     >                Y1(shift),int_mb(ma(1,ms)+taskid_i),
     >                          int_mb(ma(1,ms)),
     >                          int_mb(na(1,ms)),
     >                Y2(shift),int_mb(ma(1,ms)+taskid_i),
     >                          int_mb(ma(1,ms)),
     >                          int_mb(na(1,ms)),
     >                0.0d0,
     >                hml(shift2),int_mb(mc(1,ms)+taskid_i),
     >                            int_mb(mc(1,ms)),
     >                            int_mb(nc(1,ms)),
     >                taskid_i,taskid_j,
     >                np_i,np_j,
     >                comm_i, comm_j,
     >                dcpl_mb(work1(1)),dcpl_mb(work2(1)))

                  if (spin_orbit) then
                     call CMatrix_zgemm2(ne(ms),ne(ms),npack1_all,128,
     >                   1.0d0,
     >                   Y1(shiftso),int_mb(ma(1,ms)+taskid_i),
     >                             int_mb(ma(1,ms)),
     >                             int_mb(na(1,ms)),
     >                   Y2(shiftso),int_mb(ma(1,ms)+taskid_i),
     >                             int_mb(ma(1,ms)),
     >                             int_mb(na(1,ms)),
     >                   1.0d0,
     >                   hml(shift2),int_mb(mc(1,ms)+taskid_i),
     >                               int_mb(mc(1,ms)),
     >                               int_mb(nc(1,ms)),
     >                   taskid_i,taskid_j,
     >                   np_i,np_j,
     >                   comm_i, comm_j,
     >                   dcpl_mb(work1(1)),dcpl_mb(work2(1)))
                  end if
  20              continue
               end do
            else
               do ms=ms1,ms2
                  if (ne(ms).eq.0) go to 21
                  shift  = 1 + (ms-1)*neq(1)*npack1
     >                   + ((nb-nb1))*(neq(1)+neq(2))*npack1
                  shiftso = shift + neq(1)*npack1
                  shift2  = 1 + (ms-1)*ishift2 + (nb-nb1)*ishift3
                  call CMatrix_zgemm2c_omp(ne(ms),ne(ms),npack1_all,128,
     >                 Y1(shift),Y2(shift),int_mb(ma(1,ms)+taskid_i),
     >                                     int_mb(ma(1,ms)),
     >                                     int_mb(na(1,ms)),
     >                 dcpl_mb(mat_tmp(1)+shift2),
     >                                     int_mb(mc(1,ms)+taskid_i),
     >                                     int_mb(mc(1,ms)),
     >                                     int_mb(nc(1,ms)),
     >                 taskid_i,taskid_j,
     >                 np_i,np_j,
     >                 comm_i, comm_j,
     >                 tid,nthr,
     >                 dcpl_mb(work1(1)),dcpl_mb(work2(1)),
     >                 dcpl_mb(thrwork1(1)+tid*nn))
                  if (spin_orbit) then
                     call errquit(
     >                     "Pneb_ffw_Multiply: spin orbit not finished",
     >                      0,0)
                  end if
  21              continue
               end do
              call Pneb_w_gather(mall(mb),mpack(mb),int_mb(mindx(1,mb)),
     >                      dcpl_mb(mat_tmp(1)),hml)

            end if

         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispinq
            ishift2 = ne(1)*ne(1)  !** note that ne rather than neq is used **
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0 
         end if
         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if

         do nb=nb1,nb2
         do ms=ms1,ms2
           shift  = 1 + (ms-1)*ne(1)*npack1 +
     >                     (nb-nb1)*(ne(1)+ne(2))*npack1
           shiftso = shift + ne(1)*npack1
           shift2 = 1 + (ms-1)*ishift2      +(nb-nb1)*ishift3
           n = ne(ms)
           if (n.eq.0) go to 30

              call Cram_ccm_izgemm(nb,n,n,Y1(shift),Y2(shift),
     >                             one,zero,hml(shift2))

              if (spin_orbit) then
                call Cram_ccm_izgemm(nb,n,n,Y1(shiftso),Y2(shiftso),
     >                               one,one,hml(shift2))
              end if

 30        continue
         end do
         end do
         if (np.gt.1) call C3dB_Vector_SumAll(2*nsize,hml)
      end if

      return
      end



c     ****************************************
c     *                                      *
c     *       Pneb_ffw_hermit_Multiply       *
c     *                                      *
c     ****************************************

      subroutine Pneb_ffw_hermit_Multiply(mb,nbb,Y1,Y2,npack1,hml)
      implicit none
      integer    mb,nbb
      complex*16 Y1(*),Y2(*)
      integer    npack1
      complex*16 hml(*)

#include "Pneb.fh"

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0))

      integer nb,nb1,nb2,ms,ms1,ms2,nsize,n,np
      integer shift,shiftso,shift2,ishift2,ishift3

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      call Parallel3d_np_i(np)
      if (parallelized_j) then
      else
         nsize = Pneb_w_size(mb,nbb)
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispinq
            ishift2 = ne(1)*ne(1)  !** note that ne rather than neq is used **
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
         end if
         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if

         do nb=nb1,nb2
            do ms=ms1,ms2
              shift  = 1 + (ms-1)*ne(1)*npack1 
     >                   + (nb-nb1)*(ne(1)+ne(2))*npack1
              shiftso = shift + ne(1)*npack1
              shift2 = 1 + (ms-1)*ishift2 
     >                   +(nb-nb1)*ishift3
              n = ne(ms)
              if (n.eq.0) go to 30

              call Cram_ccm_sym_izgemm(nb,n,Y1(shift),Y2(shift),
     >                                 one,zero,hml(shift2))

              if (spin_orbit) then
                 call Cram_ccm_sym_izgemm(nb,n,Y1(shiftso),Y2(shiftso),
     >                                    one,one,hml(shift2))
              end if

 30           continue
            end do
         end do

         if (np.gt.1) call C3dB_Vector_SumAll(2*nsize,hml)
      end if

      return
      end




      
      
c     ****************************************
c     *                                      *
c     *          Pneb_fwf_Multiply           *
c     *                                      *
c     ****************************************

      subroutine Pneb_fwf_Multiply(mb,nbb,alpha,Y1,npack1,hml,beta,Y2)
      implicit none
      integer mb,nbb
      complex*16 Y1(*),Y2(*),hml(*)
      integer npack1
      complex*16 alpha,beta

      logical value
#include "bafdecls.fh"   
#include "Pneb.fh"

*     **** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,n
      integer shift,shiftso,shift2,ishift2,ishift3

      if (parallelized_j) then
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
            ishift2 = mcq(1)*ncq(1)
            ishift3 = ne(1)*ne(1)
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0
            ishift3 = 0
         end if
         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if
         do nb=nb1,nb2
            if (mparallelized) then
            else
            end if

         end do

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispinq
            ishift2 = ne(1)*ne(1)  !** note that ne rather than neq is used **
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0 
         end if
         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if


         do nb=nb1,nb2
            do ms=ms1,ms2           
               shift = 1 + (ms-1)*ne(1)*npack1
     >                   +  (nb-nb1)*(ne(1)+ne(2))*npack1
               shiftso = shift + ne(1)*npack1
               shift2  = 1 + (ms-1)*ishift2 
     >                     + (nb-nb1)*ishift3
               n = ne(ms)
               if (n.eq.0) go to 30

                 call Cram_cmc_zgemm(nb,n,n,Y1(shift),Y2(shift),
     >                               alpha,beta,hml(shift2))

                 if (spin_orbit) then
                    call Cram_cmc_zgemm(nb,n,n,Y1(shiftso),Y2(shiftso),
     >                                  alpha,beta,hml(shift2))
                 end if

 30            continue
            end do
         end do
      end if

      return
      end
      

c     ****************************************
c     *                                      *
c     *          Pneb_fwf_multiplyAdd        *
c     *                                      *
c     ****************************************

      subroutine Pneb_fwf_multiplyAdd(mb,nbb,Y1,hml,npack1,Y2)
      implicit none
      integer mb,nbb
      complex*16 Y1(*),Y2(*),hml(*)
      integer npack1

#include "Pneb.fh"

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0))

      integer nb,nb1,nb2,ms,ms1,ms2,n
      integer shift,shiftso,shift2,ishift2,ishift3

      if (parallelized_j) then

      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispinq
            ishift2 = ne(1)*ne(1)  !** note that ne rather than neq is used **
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0 
         end if
         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if


         do nb=nb1,nb2
         do ms=ms1,ms2           
            shift = 1 + (ms-1)*ne(1)*npack1 +
     >               (nb-nb1)*(ne(1)+ne(2))*npack1
            shiftso = shift + ne(1)*npack1
            shift2  = 1 + (ms-1)*ishift2      +(nb-nb1)*ishift3
            n = ne(ms)
            if (n.eq.0) go to 30

              call Cram_cmc_zgemm(nb,n,n,Y1(shift),Y2(shift),
     >                            one,one,hml(shift2))

              if (spin_orbit) then
                 call Cram_cmc_zgemm(nb,n,n,Y1(shiftso),Y2(shiftso),
     >                               one,one,hml(shift2))
              end if

 30         continue
         end do
         end do

      end if

      return
      end
      
c     ****************************************
c     *                                      *
c     *          Pneb_fwf_multiply           *
c     *          alternative                 *
c     ****************************************

      subroutine Pneb_fwf_multiplyAlt(mb,nbb,Y1,hml,npack1,Y2)
      implicit none
      integer mb,nbb
      complex*16 Y1(*),Y2(*),hml(*)
      integer npack1

#include "Pneb.fh"

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0))

      integer nb,nb1,nb2,ms,ms1,ms2,n
      integer shift,shiftso,shift2,ishift2,ishift3

      if (parallelized_j) then
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispinq
            ishift2 = ne(1)*ne(1)  !** note that ne rather than neq is used **
         else
            ms1 = mb
            ms2 = mb
            ishift2 = 0 
         end if
         if (nbb.eq.0) then
           nb1 = 1
           nb2 = nbrillq
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         else
           nb1 = nbb
           nb2 = nbb
           ishift3 = 0
         end if

         do nb=nb1,nb2
         do ms=ms1,ms2           
            shift = 1 + (ms-1)*ne(1)*npack1 +
     >               (nb-nb1)*(ne(1)+ne(2))*npack1
            shiftso = shift + ne(1)*npack1
            shift2  = 1 + (ms-1)*ishift2      +((nb-nb1))*ishift3
            n = ne(ms)
            if (n.eq.0) go to 30

              call Cram_cmc_zgemm(nb,n,n,Y1(shift),Y2(shift),
     >                            one,zero,hml(shift2))

              if (spin_orbit) then
                 call Cram_cmc_zgemm(nb,n,n,Y1(shiftso),Y2(shiftso),
     >                               one,zero,hml(shift2))
              end if

 30         continue
         end do
         end do
      end if

      return
      end
      


c     ****************************************
c     *                                      *
c     *          Pneb_w_axpy                 *
c     *                                      *
c     ****************************************

      subroutine Pneb_w_Axpy(mb,nbb,alpha,fml,hml)
      implicit none
      integer mb,nbb
      complex*16 fml(*),hml(*),alpha

#include "Pneb.fh"

*     **** local variables ****
      integer nsize

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      nsize = Pneb_w_size(mb,nbb)
      call zaxpy(nsize,alpha,fml,1,hml,1) 
      return
      end



c     ****************************************
c     *                                      *
c     *          Pneb_w_scal                 *
c     *                                      *
c     ****************************************

      subroutine Pneb_w_scal(mb,nbb,alpha,hml)
      implicit none
      integer mb,nbb
      real*8 alpha
      complex*16 hml(*)

#include "Pneb.fh"

*     **** local variables ****
      integer nsize

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      nsize = Pneb_w_size(mb,nbb)
      call dscal(2*nsize,alpha,hml,1)
      return
      end



c     ccccccccccccccccccccccccccccccccccccccccc
c     *                                       *
c     *            Pneb_w_copy                *
c     *                                       *
c     ccccccccccccccccccccccccccccccccccccccccc

      subroutine  Pneb_w_copy(mbb,nbb,hml,gml)
      implicit none
      integer    nbb,mbb
      complex*16 hml(*),gml(*)
     
#include "Pneb.fh"

c     **** local variables ****
      integer nsize

c     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size
    
      nsize = Pneb_w_size(mbb,nbb)
      call dcopy(2*nsize,hml,1,gml,1)

      return 
      end 



c     ccccccccccccccccccccccccccccccccccccccccc
c     *                                       *
c     *            Pneb_w_copyto0             *
c     *                                       *
c     ccccccccccccccccccccccccccccccccccccccccc
c
c A special copy
c
      subroutine  Pneb_w_copyto0(mbb,nbb,hml,gml)
      implicit none
      integer    nbb,mbb
      complex*16 hml(*),gml(*)

#include "Pneb.fh"

c     **** local variables ****
      integer nsize,shift

c     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      nsize = Pneb_w_size(mbb,nbb)
      shift = 1 + (mbb-1)*ne(1)*ne(1)
      call dcopy(2*nsize,hml,1,gml(shift),1)

      return
      end

c     cccccccccccccccccccccccccccccccccccccccc
c     *                                      *
c     *           Pneb_m_copy                *
c     *                                      *
c     cccccccccccccccccccccccccccccccccccccccc

      subroutine  Pneb_m_copy(mbb,nbb,hml,gml)
      implicit none
      integer nbb,mbb
      real*8 hml(*),gml(*)

#include "Pneb.fh"
     
c     **** local variables ****
      integer nsize

c     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size

      nsize = Pneb_m_size(mbb,nbb)
      call dcopy(2*nsize,hml,1,gml,1)

      return 
      end 

c     ccccccccccccccccccccccccccccccccccc
c     *                                 *
c     *         Pneb_w_diag             *
c     *                                 *
c     ccccccccccccccccccccccccccccccccccc

      subroutine Pneb_w_diag(mbb,nbb,S,V)
      implicit none
      integer    mbb,nbb
      real*8     S(*)
      complex*16 V(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

c     *** local variables ***
      integer MASTER
      parameter (MASTER=0)
      logical value
      integer nb,nb1,nb2,ms,ms1,ms2,i,j,k
      integer rwork(2),work(2),lwork,ierr
      integer vindex,sindex
      integer W(2),VL(2),VR(2)

      if (nbb.eq.0) then
        nb1=1
        nb2=nbrillq
      else
        nb1=nbb
        nb2=nbb
      end if   

      if (mbb.eq.0) then
         ms1 = 1
         ms2 = ispinq
      else
         ms1 = mbb
         ms2 = mbb
      end if


c     **** get temp memory from stack ****
      lwork=3*ne(1)*ne(1)
      value =BA_push_get(mt_dbl,lwork,'rwork',rwork(2),rwork(1))
      value =value.and.BA_push_get(mt_dcpl,lwork,'work',work(2),work(1))
      if (.not.value) call errquit('Pneb_w_diag: out of stack',0,MA_ERR)

      sindex = 1
      vindex = 1
      do nb=nb1,nb2
      do ms=ms1,ms2
        if (ne(ms).gt.0) then
           ierr = 0
           call ZHEEV('V','L',ne(ms),
     >                 V(vindex),ne(ms),
     >                 S(sindex),
     >                 dcpl_mb(work(1)),lwork,
     >                 dbl_mb(rwork(1)),
     >                 ierr)
           if (ierr.ne.0) 
     >       call errquit("Pneb_w_diag: ZHEEV failed",ierr,0)

        end if

        vindex=vindex + ne(ms)*ne(ms)
        sindex=sindex + ne(ms)       
      end do
      end do
      call Parallela_Brdcst_values(1,MASTER,2*(vindex-1),V)
      call Parallela_Brdcst_values(1,MASTER,(sindex-1),S)

c     ***** pop stack memory *****
      value=          BA_pop_stack(work(2))
      value=value.and.BA_pop_stack(rwork(2))
      if (.not.value) then
        call errquit('Pneb_w_diag:error popping stack',0,MA_ERR)
      end if

      return
      end      



c     ccccccccccccccccccccccccccccccccccc
c     *                                 *
c     *         Pneb_wne_diag           *
c     *                                 *
c     ccccccccccccccccccccccccccccccccccc

      subroutine Pneb_wne_diag(mbb,nbb,nein,S,V)
      implicit none
      integer    mbb,nbb,nein(2)
      real*8     S(*)
      complex*16 V(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

c     *** local variables ***
      integer MASTER
      parameter (MASTER=0)
      logical value
      integer nb,nb1,nb2,ms,ms1,ms2
      integer rwork(2),work(2),lwork,ierr
      integer vindex,sindex

      if (nbb.eq.0) then
        nb1=1
        nb2=nbrillq
      else
        nb1=nbb
        nb2=nbb
      end if   

      if (mbb.eq.0) then
         ms1 = 1
         ms2 = ispinq
      else
         ms1 = mbb
         ms2 = mbb
      end if


c     **** get temp memory from stack ****
      lwork=3*nein(1)*nein(1)
      value =BA_push_get(mt_dbl,lwork,'rwork',rwork(2),rwork(1))
      value =value.and.BA_push_get(mt_dcpl,lwork,'work',work(2),work(1))
      if (.not.value) 
     >  call errquit('Pneb_wne_diag: out of stack',0,MA_ERR)

      sindex = 1
      vindex = 1
      do nb=nb1,nb2
      do ms=ms1,ms2
        if (nein(ms).gt.0) then
           ierr = 0
           call ZHEEV('V','L',nein(ms),
     >                 V(vindex),nein(ms),
     >                 S(sindex),
     >                 dcpl_mb(work(1)),lwork,
     >                 dbl_mb(rwork(1)),
     >                 ierr)
           if (ierr.ne.0) 
     >       call errquit("Pneb_wne_diag: ZHEEV failed",ierr,0)
        endif

        vindex=vindex + nein(ms)*nein(ms)
        sindex=sindex + nein(ms)       
      end do
      end do
      call Parallela_Brdcst_values(1,MASTER,2*(vindex-1),V)
      call Parallela_Brdcst_values(1,MASTER,(sindex-1),S)
             

c     ***** pop stack memory *****
      value=          BA_pop_stack(work(2))
      value=value.and.BA_pop_stack(rwork(2))
      if (.not.value) then
        call errquit('Pneb_wne_diag:error popping stack',0,MA_ERR)
      end if

      return
      end      



c     ccccccccccccccccccccccccccccccccccccccccccc
c     *                                         *
c     *               Pneb_SVD                  *
c     *                                         *
c     ccccccccccccccccccccccccccccccccccccccccccc

      subroutine Pneb_SVD(mbb,nbb,npack1,Y,U,S,V)
      implicit none
      integer mbb,nbb,npack1    
      real*8     S(*)  
      complex*16 Y(*),U(*),V(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

c     ***** local variables ******
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0))

      integer soshift,shift1,shift2,n,neall
      integer ms,ms1,ms2,nb,nb1,nb2
      real*8  tsum
      integer tmp(2)

      neall=ne(1)+ne(2)*(ispinq-1)
      if (mbb.eq.0) then
         ms1 = 1
         ms2 = ispinq
      else
         ms1 = mbb
         ms2 = mbb
         neall = ne(mbb)
      end if
      
      if (nbb.eq.0) then
         nb1 = 1
         nb2 = nbrillq
         neall = neall * nbrillq
      else
         nb1 = nbb
         nb2 = nbb
       end if

c     ***  compute C = Y'*Y ( ne x ne )= V*S*S*V' ****
      call Pneb_ffw_hermit_Multiply(mbb,nbb,Y,Y,npack1,V)


c     ***  compute S2=S*S and V from C  **** 
      call Pneb_w_diag(mbb,nbb,S,V)


c     **** sqrt(S*S) ****
      do n=1,neall
         S(n) = dsqrt(dabs(S(n)))
      end do

c     *** compute U*S=Y*V ****
      call Pneb_fwf_Multiply(mbb,nbb,one,Y,npack1,V,zero,U)

*     **** normalize U*S ****
      if (.not.BA_push_get(mt_dbl,neq(1),'tmp',tmp(2),tmp(1)))
     > call errquit('Dneb_f_SVD: out of stack memory',0,MA_ERR)

      soshift = neq(1)*npack1
      do nb=nb1,nb2
         do ms=ms1,ms2
            if (neq(ms).le.0) go to 30
               shift1 = 1 + (ms-1)* neq(1)        *npack1 
     >                    + (nb-nb1)*(neq(1)+neq(2))*npack1
               shift2 = shift1
               do n=1,neq(ms)
                  call Cram_cc_idot(nb,
     >                              U(shift1),
     >                              U(shift1),
     >                              dbl_mb(tmp(1)+n-1))

                  if (spin_orbit) then
                     call Cram_cc_idot(nb,
     >                                U(shift1+soshift),
     >                                U(shift1+soshift),
     >                                tsum)
                     dbl_mb(tmp(1)+n-1) = dbl_mb(tmp(1)+n-1) + tsum
                  end if
                  shift1 = shift1 + npack1
               end do
               call C3dB_Vector_SumAll(neq(ms),dbl_mb(tmp(1)))

               do n=1,neq(ms)
                  dbl_mb(tmp(1)+n-1) = 1.0d0/dsqrt(dbl_mb(tmp(1)+n-1))
               end do

               do n=1,neq(ms)
                  call Cram_c_SMul1(nb,dbl_mb(tmp(1)+n-1),U(shift2))
                  if (spin_orbit) then
                     call Cram_c_SMul1(nb,dbl_mb(tmp(1)+n-1),
     >                                 U(shift2+soshift))
                  end if
                  shift2 = shift2 + npack1
               end do
                 
 30         continue
         end do
      end do

      if (.not.BA_pop_stack(tmp(2)))
     >   call errquit('Pneb_f_SVD: popping stack memory',0,MA_ERR)
      return 
      end  





c     ccccccccccccccccccccccccccccccccccccccccc
c     *                                       *
c     *            Pneb_makeOrtho             *
c     *                                       *
c     ccccccccccccccccccccccccccccccccccccccccc

      subroutine Pneb_makeOrtho(mb,nbb,npack1,A)
      implicit none
      integer mb,nbb,npack1
      complex*16 A(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

c     ** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,index1,index2,index3,index4
      integer j,k,nemax
      integer bindex,nbshift,soshift
      complex*16 wc,wc1,wc2       
      real*8 w,w1,w2
      logical value


      if (nbb.eq.0) then
        nb1=1
        nb2=nbrillq
      else
        nb1=nbb
        nb2=nbb
      end if
      if (spin_orbit) then
        nemax=ne(1)
        soshift=npack1*nemax
        nbshift=2*soshift
        bindex=1
        do nb=nb1,nb2
         bindex=((nb-nb1))*nbshift
         do k=1,nemax
           index1=bindex+(k-1)*npack1+1
           call Cram_cc_dot(nb,A(index1),A(index1),w1)
           index2=index1+soshift
           call Cram_cc_dot(nb,A(index2),A(index2),w2)
           w=w1+w2
           if (w.gt.1.d-9) then
             w=1.d0/dsqrt(w)
           else
             write(*,*)"Pneb SO make ortho"
             write(*,*)"Pneb_makeOrtho::this should not happen"
             write(*,*)"Pneb_makeOrtho: < psi | psi > = 0 K= ",k,w     
             w=1.d0/dsqrt(2.0d0*w1)
             call Cram_c_copy(nb,A(index1),A(index2))
           end if
c           call Cram_c_SMul(nb,w,A(index1),A(index1))
c           call Cram_c_SMul(nb,w,A(index2),A(index2))
           call Cram_c_SMul1(nb,w,A(index1))
           call Cram_c_SMul1(nb,w,A(index2))
           do j=k+1,nemax
             index3=bindex+(j-1)*npack1+1
             index4=index3+soshift
             call Cram_cc_zdot(nb,A(index1),A(index3),wc1)
             call Cram_cc_zdot(nb,A(index2),A(index4),wc2)
             wc=wc1+wc2
             wc=-wc
             call Cram_cc_zaxpy(nb,wc,A(index1),A(index3))
             call Cram_cc_zaxpy(nb,wc,A(index2),A(index4))
           end do
         end do
        end do
      else
        if (mb.eq.0) then
          ms1=1
          ms2=ispin
        else
          ms1=mb
          ms2=mb
        end if

        nemax=ne(1)+ne(2)
        nbshift=npack1*nemax

        do nb=nb1,nb2
        do ms=ms1,ms2
          bindex=(nb-nb1)*nbshift+(ms-1)*ne(1)*npack1
          do k=1,ne(ms)
           index1=bindex+(k-1)*npack1+1
           call Cram_cc_dot(nb,A(index1),A(index1),w)

           if (w.gt.1.0d-9) then
             w=1.0d0/dsqrt(w)
           else
            write(*,*)"Pneb_makeOrtho::this should not happen"
            write(*,*)"Pneb_makeOrtho: < psi | psi > = 0 K= ",k,w     
            w=1.0d0
           end if
c           call Cram_c_SMul(nb,w,A(index1),A(index1))
           call Cram_c_SMul1(nb,w,A(index1))

           do j=k+1,ne(ms)
             index3=bindex+(j-1)*npack1+1
             call Cram_cc_zdot(nb,A(index1),A(index3),wc)
             wc=-wc
             call Cram_cc_zaxpy(nb,wc,A(index1),A(index3))
           end do
          end do
        end do
        end do
      end if
      return 
      end



c     ccccccccccccccccccccccccccccccccccccccccccc
c     *                                         *
c     *           Pneb_orthoCheckMake           *
c     *                                         *
c     ccccccccccccccccccccccccccccccccccccccccccc

      subroutine Pneb_orthoCheckMake(coutput,mb,nbb,npack1,A)
      implicit none
      logical coutput
      integer mb,nbb,npack1
      complex*16 A(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "Pneb.fh"

c     ccccc local variables ccccc
      integer taskid,MASTER
      PARAMETER(MASTER=0)

      logical oprint
      integer asize,indx
      integer nb,nb1,nb2,ms,ms1,ms2,sm2(2),sm3(2)
      integer indx2,indx3
      logical gram,value
      real*8  sum1,sum2,sum3

*     **** external functions ****
      logical  control_print
      external control_print
      
      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_high))
      if (nbb.eq.0) then
         nb1=1
         nb2=nbrillq
         asize=2*nbrillq
      else
         nb1=nbb
         nb2=nbb
         asize=2
      end if
      if (mb.eq.0) then
         ms1=1
         ms2=ispin
      else
         ms1=1
         ms2=1
      end if

*     **** allocate tmp space from the stack ****
      value =           BA_push_get(mt_dbl,asize,'sm2',sm2(2),sm2(1))
      value = value.and.BA_push_get(mt_dbl,asize,'sm3',sm3(2),sm3(1))
      if (.not.value) then
       call errquit('Pneb_orthCheckMake: out of stack memory',0,MA_ERR)
      end if


      call Pneb_ff_trace(mb,nbb,A,A,npack1,dbl_mb(sm2(1)))

      gram = .false.  
      if (spin_orbit) then
         indx=sm2(1)
         do nb=nb1,nb2
            sum1=dble(ne(1))
            sum2=dbl_mb(indx)
            indx=indx+1
            if (dabs(sum2-sum1).gt.1.0d-10) gram= .true.
         end do 
      else
         indx=sm2(1)
         do nb=nb1,nb2
         do ms=ms1,ms2
            sum1=dble(ne(ms))
            sum2=dbl_mb(indx)
            indx=indx+1         
            if (dabs(sum2-sum1).gt.1.0d-12) gram = .true.
         end do
         end do
      end if

      if (gram) then
         call Pneb_makeOrtho(mb,nbb,npack1,A)
         call Pneb_ff_trace(mb,nbb,A,A,npack1,dbl_mb(sm3(1)))

         if (oprint) then
         if (spin_orbit) then
            indx2=sm2(1)
            indx3=sm3(1)
            do nb=nb1,nb2
               sum1=dble(ne(1))
               sum2=dbl_mb(indx2)
               sum3=dbl_mb(indx3)
               indx2=indx2+1
               indx3=indx3+1
               write(*,1400) 3,nb
               write(*,1401) dabs(sum2-sum1),sum2
               write(*,1402) dabs(sum3-sum1),sum3
            end do          
         else
            indx2=sm2(1)
            indx3=sm3(1)
            do nb=nb1,nb2
            do ms=ms1,ms2
               sum1=dble(ne(ms))
               sum2=dbl_mb(indx2)
               sum3=dbl_mb(indx3)
               indx2=indx2+1
               indx3=indx3+1
               write(*,1400) ms,nb
               write(*,1401) dabs(sum2-sum1),sum2
               write(*,1402) dabs(sum3-sum1),sum3
            end do
            end do
         end if
         end if
      end if

*     **** pop stack ****
      value = BA_pop_stack(sm3(2))
      value = value.and.BA_pop_stack(sm2(2))
      if (.not. value)
     > call errquit('pneb_orthocheckmake:popping stack memory',0,0)

 1400 FORMAT('Warning - Orthogonalization performed, spin,zone:',I2,I4)
 1401 FORMAT(8x,'- error(before)=', E14.6,' (',E14.6,')')
 1402 FORMAT(8x,'- error(after)= ', E14.6,' (',E14.6,')')
      return
      end



c     ccccccccccccccccccccccccccccccccccccccccccc
c     *                                         *
c     *           Pneb_orthoCheckMake_tag       *
c     *                                         *
c     ccccccccccccccccccccccccccccccccccccccccccc

      subroutine Pneb_orthoCheckMake_tag(coutput,mb,nbb,npack1,A_tag)
      implicit none
      logical coutput
      integer mb,nbb,npack1
      integer A_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "Pneb.fh"

c     ccccc local variables ccccc
      integer taskid,MASTER
      PARAMETER(MASTER=0)

      logical oprint,value
      integer nb,nb1,nb2,ms,ms1,ms2,nbq
      integer a_shift,nshift,l,i
      integer sm2(2),sm3(2)
      real*8  sum1,sum2,sum3

*     **** external functions ****
      logical  control_print
      integer  cpsi_data_nsize,cpsi_data_get_chnk
      external control_print
      external cpsi_data_nsize,cpsi_data_get_chnk

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_high))
      oprint = (taskid.eq.MASTER)

      if (nbb.eq.0) then
         nb1=1
         nb2=nbrillq
      else
         nb1=nbb
         nb2=nbb
      end if
      if (mb.eq.0) then
         ms1=1
         ms2=ispinq
      else
         ms1=1
         ms2=1
      end if

      value =           BA_push_get(mt_dbl,2*nbrill,'sm2',sm2(2),sm2(1))
      value = value.and.BA_push_get(mt_dbl,2*nbrill,'sm3',sm3(2),sm3(1))
      if (.not.value) 
     >   call errquit('Pneb_orthCheckMake_tag:out of stack',0,MA_ERR)
      call dcopy(2*nbrill,0.0d0,0,dbl_mb(sm2(1)),1)
      call dcopy(2*nbrill,0.0d0,0,dbl_mb(sm3(1)),1)

      nshift = cpsi_data_nsize(A_tag)

      do nbq=nb1,nb2
         a_shift = cpsi_data_get_chnk(A_tag,nbq)
         do ms=ms1,ms2
            sum1=dble(ne(ms))
            call Pneb_ff_trace(ms,nbq,
     >                         dbl_mb(a_shift),
     >                         dbl_mb(a_shift),npack1,
     >                         sum2)
            if (dabs(sum2-sum1).gt.1.d-10) then
              call Pneb_makeOrtho(ms,nbq,npack1,dbl_mb(a_shift))
              call Pneb_ff_trace(ms,nbq,
     >                         dbl_mb(a_shift),
     >                         dbl_mb(a_shift),npack1,
     >                         sum3)
              call K1dB_qtok(nbq,nb)
              i = (ms-ms1)+2*(nb-1)
              dbl_mb(sm2(1)+i) = sum2
              dbl_mb(sm3(1)+i) = sum3
            end if
         end do
      end do
      call K1dB_Vector_SumAll(2*nbrill,dbl_mb(sm2(1)))
      call K1dB_Vector_SumAll(2*nbrill,dbl_mb(sm3(1)))
      if (oprint) then
      do nb=1,nbrill
         do ms=ms1,ms2
            sum1=dble(ne(ms))
            i = (ms-ms1)+2*(nb-1)
            if (dabs(dbl_mb(sm2(1)+i)).gt.1.d-10) then
              if (spin_orbit) then
                 write(*,1400) nb
              else
                 write(*,1401) ms,nb
              end if
              write(*,1402) dabs(dbl_mb(sm2(1)+i)-sum1),dbl_mb(sm2(1)+i)
              write(*,1403) dabs(dbl_mb(sm3(1)+i)-sum1),dbl_mb(sm3(1)+i)
            end if
         end do
      end do
      end if
      value =           BA_pop_stack(sm3(2))
      value = value.and.BA_pop_stack(sm2(2))
      if (.not. value)
     >   call errquit('Pneb_orthCheckMake_tag:popping stack',0,MA_ERR)

 1400 FORMAT(
     >  'Warning - Orthogonalization performed, spin_orbit,zone:',I4)
 1401 FORMAT('Warning - Orthogonalization performed, spin,zone:',I2,I4)
 1402 FORMAT(8x,'- error(before)=', E14.6,' (',E14.6,')')
 1403 FORMAT(8x,'- error(after)= ', E14.6,' (',E14.6,')')

      return
      end






c     ccccccccccccccccccccccccccccccccccccccccccccccc
c     *                                             *
c     *             Pneb_orthoCheck                 *
c     *                                             *
c     ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine Pneb_orthoCheck(ix1,mb,nbb,npack1,A)
      implicit none
      integer mb,nbb,npack1,ix1
      complex*16 A(*)

#include "errquit.fh"
#include "Pneb.fh"

c     ***** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      integer nb,nb1,nb2,ms,ms1,ms2,j,k,bindex
      integer index11,index21,index22,index12,bshift,soshift,nemax
      real*8     error
      complex*16 w,wc1,wc2
      character*255 full_filename
      character*7   p_index_name
      external      p_index_name

c     ccccccccc first line ccccccccc
      if (nbb.eq.0) then
         nb1=1
         nb2=nbrillq
      else
         nb1=nbb
         nb2=nbb
      end if
      if (mb.eq.0) then
         ms1=1
         ms2=ispin
      else
         ms1=mb
         ms2=mb
      end if

*     **** produce CHECK FILE ****
      call util_file_name('ORTHOCHECK'//p_index_name(ix1),
     >                             .true.,
     >                             .false.,
     >                             full_filename)
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         open(unit=17,file=full_filename,form='formatted')
      end if

      if (.not.im_init) then
         call errquit("PNEB : I need to be initialized",0,0)
      end if

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) THEN
         write(17,1350)
         write(17,*)ne(1),ne(2)
         call flush(17)
      end if
      if (spin_orbit) then
        nemax=ne(1)
        bshift=npack1*nemax*2
        soshift=npack1*nemax
        do nb=nb1,nb2
          bindex=bshift*((nb-nb1))
          do k=1,nemax
            index11=1+((k-1)+((nb-nb1))*nemax*2)*npack1
            index12=index11+soshift
            do j=1,nemax
              index21=1+((j-1)+((nb-nb1))*nemax*2)*npack1
              index22=index21+soshift
              call Cram_cc_zdot(nb,A(index11),A(index21),wc1)
              call Cram_cc_zdot(nb,A(index12),A(index22),wc2)
              w=wc1+wc2
              if (k.eq.j) then
                error=dabs(1.-dble(w))+dabs(aimag(w))
              else
                error=dabs(dble(w))+dabs(dimag(w))
              end if
              if (taskid.eq.MASTER) then
                write(17,1360) nb,3,k,j,w
              end if
            end do
          end do
        end do
        if (taskid.eq.MASTER) write(17,1370) error
        if (taskid.eq.MASTER) write(17,1380)
      else
        bshift=(ne(1)+ne(2))*npack1
        do nb=nb1,nb2
        do ms=ms1,ms2
          bindex=bshift*((nb-nb1))+(ms-1)*npack1*ne(1)
          do k=1,ne(ms)
            index11=1+bindex+(k-1)*npack1
            do j=1,ne(ms)
              index21=1+bindex+(j-1)*npack1
              call Cram_cc_zdot(nb,A(index11),A(index21),w)
              if (k.eq.j) then
                 error=dabs(1.-dble(w))+dabs(aimag(w))
              else
                 error=dabs(dble(w))+dabs(dimag(w))
              end if
              if (taskid.eq.MASTER) then
                write(17,1360) nb,ms,k,j,w
              end if
            end do
          end do
        end do
        end do
        if (taskid.eq.MASTER) write(17,1370) error
        if (taskid.eq.MASTER) write(17,1380)
      end if
      if (taskid.eq.MASTER) close(17)
      return
 1350 FORMAT(/'******** band structure orthonormality **********')
 1360 FORMAT(2I3,2I3,'(',2E18.7,')')
 1370 FORMAT('ERROR = ',E18.7)
 1380 FORMAT(/'*************************************************')
      end





c     ccccccccccccccccccccccccccccccccccccccccccccccc
c     *                                             *
c     *             Pneb_orthoCheck_tag             *
c     *                                             *
c     ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine Pneb_orthoCheck_tag(ix1,mb,nbb,npack1,A_tag)
      implicit none
      integer mb,nbb,npack1,ix1
      integer A_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

c     ***** local variables ****
      integer taskid,MASTER
      parameter(MASTER=0)

      integer nb,nb1,nb2,ms,ms1,ms2,j,k,bindex
      integer index11,index21,index22,index12,a_shift,soshift,nemax
      real*8     error
      complex*16 w,wc1,wc2
      character*255 full_filename

*     **** external functions ***
      character*7   p_index_name
      integer       cpsi_data_get_chnk
      external      p_index_name
      external      cpsi_data_get_chnk

c     ccccccccc first line ccccccccc
      if (nbb.eq.0) then
         nb1=1
         nb2=nbrillq
      else
         nb1=nbb
         nb2=nbb
      end if
      if (mb.eq.0) then
         ms1=1
         ms2=ispin
      else
         ms1=mb
         ms2=mb
      end if

*     **** produce CHECK FILE ****
      call util_file_name('ORTHOCHECK'//p_index_name(ix1),
     >                             .true.,
     >                             .false.,
     >                             full_filename)
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         open(unit=17,file=full_filename,form='formatted')
      end if

      if (.not.im_init) then
         call errquit("PNEB : I need to be initialized",0,0)
      end if

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) THEN
         write(17,1350)
         write(17,*)ne(1),ne(2)
         call flush(17)
      end if
      if (spin_orbit) then
        nemax=ne(1)
        soshift=2*npack1*nemax
        do nb=nb1,nb2
          a_shift = cpsi_data_get_chnk(A_tag,nb)
          do k=1,nemax
            index11=a_shift+(k-1)*2*npack1
            index12=index11+soshift
            do j=1,nemax
              index21=a_shift+(j-1)*2*npack1
              index22=index21+soshift
              call Cram_cc_zdot(nb,dbl_mb(index11),dbl_mb(index21),wc1)
              call Cram_cc_zdot(nb,dbl_mb(index12),dbl_mb(index22),wc2)
              w=wc1+wc2
              if (k.eq.j) then
                error=dabs(1.-dble(w))+dabs(aimag(w))
              else
                error=dabs(dble(w))+dabs(dimag(w))
              end if
              if (taskid.eq.MASTER) then
                write(17,1360) nb,3,k,j,w
              end if
            end do
          end do
        end do
        if (taskid.eq.MASTER) write(17,1370) error
        if (taskid.eq.MASTER) write(17,1380)
      else
        do nb=nb1,nb2
          a_shift = cpsi_data_get_chnk(A_tag,nb)
          do ms=ms1,ms2
          bindex=(ms-1)*2*npack1*ne(1)
          do k=1,ne(ms)
            index11=a_shift+bindex+(k-1)*2*npack1
            do j=1,ne(ms)
              index21=a_shift+bindex+(j-1)*2*npack1
              call Cram_cc_zdot(nb,dbl_mb(index11),dbl_mb(index21),w)
              if (k.eq.j) then
                 error=dabs(1.-dble(w))+dabs(dimag(w))
              else
                 error=dabs(dble(w))+dabs(dimag(w))
              end if
              if (taskid.eq.MASTER) then
                write(17,1360) nb,ms,k,j,w
              end if
            end do
          end do
        end do
        end do
        if (taskid.eq.MASTER) write(17,1370) error
        if (taskid.eq.MASTER) write(17,1380)
      end if
      if (taskid.eq.MASTER) close(17)
      return
 1350 FORMAT(/'******** band structure orthonormality **********')
 1360 FORMAT(2I3,2I3,'(',2E18.7,')')
 1370 FORMAT('ERROR = ',E18.7)
 1380 FORMAT(/'*************************************************')
      end



c     ****************************************
c     *                                      *
c     *          Pneb_ff_traceall            *
c     *                                      *
c     ****************************************
      subroutine Pneb_ff_traceall(mb,nbb,Y1,Y2,npack1,sum)
      implicit none
      integer mb,nbb
      complex*16 Y1(*),Y2(*)
      real*8 sum
      integer npack1

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

*     **** local variables ****
      integer nb1,nb2,nb,ms1,ms2,ms
      integer index,tmp2(2)
      logical value

c     ** external functions
      real*8 brillioun_weight
      external brillioun_weight

      sum=0.0d0
      if (nbb.eq.0) then
        nb1 = 1
        nb2 = nbrillq
      else
        nb1 = nbb
        nb2 = nbb
      end if

      value=BA_push_get(mt_dbl,(2*nbrillq),'tmp2',tmp2(2),tmp2(1))
      if (.not.value) then
       call errquit('out of stack memory',0,MA_ERR)
      end if

      call Pneb_ff_trace(mb,nbb,Y1,Y2,npack1,dbl_mb(tmp2(1)))

      if (spin_orbit) then
        sum=0.0d0
        index=tmp2(1)
        do nb=nb1,nb2
          sum=sum+brillioun_weight(nb)*dbl_mb(index)
          index=index+1
        end do
      else
        if (mb.eq.0) then
          ms1=1
          ms2=ispin
        else
           ms1=mb
           ms2=mb
        end if
        sum=0.0d0
        index=tmp2(1)
        do nb=nb1,nb2
          do ms=ms1,ms2
            sum=sum+brillioun_weight(nb)*dbl_mb(index)
            index=index+1
          end do
        end do
      end if
      value=BA_pop_stack(tmp2(2))
      if (.not.value) then
        call errquit('error popping the stack',0,MA_ERR)
      end if

      call K1dB_SumAll(sum)

      return
      end


c     ****************************************
c     *                                      *
c     *          Pneb_ff_traceall_tag        *
c     *                                      *
c     ****************************************
      subroutine Pneb_ff_traceall_tag(mb,nbb,Y1_tag,Y2_tag,sum)
      implicit none
      integer mb,nbb
      integer Y1_tag,Y2_tag
      real*8 sum

#include "bafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

*     **** local variables ****
      integer nb1,nb2,nb,ms1,ms2,ms
      integer npack1,nshift1,y1_shift,y2_shift
      real*8  tmp

c     *** external functions ***
      integer  cpsi_data_nsize,cpsi_data_get_chnk
      real*8   brillioun_weight
      external cpsi_data_nsize,cpsi_data_get_chnk
      external brillioun_weight

      if (nbb.eq.0) then
        nb1 = 1
        nb2 = nbrillq
      else
        nb1 = nbb
        nb2 = nbb
      end if

      nshift1 = cpsi_data_nsize(Y1_tag)
      npack1  = nshift1/2

      if (mb.eq.0) then
        ms1=1
        ms2=ispinq
      else
         ms1=mb
         ms2=mb
      end if
      sum=0.0d0
      do nb=nb1,nb2
        y1_shift = cpsi_data_get_chnk(Y1_tag,nb)
        y2_shift = cpsi_data_get_chnk(Y2_tag,nb)
        do ms=ms1,ms2
          call Pneb_ff_trace(ms,nb,
     >                       dbl_mb(y1_shift),
     >                       dbl_mb(y2_shift),npack1,
     >                       tmp)
          sum=sum+brillioun_weight(nb)*tmp
          !y1_shift = y1_shift + nshift1
          !y2_shift = y2_shift + nshift1
        end do
      end do
      call K1dB_SumAll(sum)

      return
      end




c     ****************************************
c     *                                      *
c     *          Pneb_ff_trace               *
c     *                                      *
c     ****************************************

      subroutine Pneb_ff_trace(mb,nbb,Y1,Y2,npack1,hml)
      implicit none
      integer mb,nbb
      complex*16 Y1(*),Y2(*)
      integer npack1
      real*8 hml(*)

#include "Pneb.fh"

*     **** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,nemax,nsize
      integer k,indx1,indexm
      integer nbshift,msshift,soshift
      integer np
      real*8 m1,m2,sum

      call Parallel3d_np_i(np)

      nemax = ne(1)+ne(2)

      if (nbb.eq.0) then
        nb1 = 1
        nb2 = nbrillq
      else
        nb1 = nbb
        nb2 = nbb
      end if

c     ***** spin-orbit dft *****
      if (spin_orbit) then
         soshift = ne(1)*npack1
         sum=0.0d0
         indexm=1
         do nb = nb1,nb2
            nbshift = ((nb-nb1))*2*soshift
            indx1 = 1 + nbshift
            sum=0.0d0
            do k=1,ne(1)
              call Cram_cc_idot(nb,Y1(indx1),
     >                             Y2(indx1),
     >                             m1)
              call Cram_cc_idot(nb,Y1(indx1+soshift),
     >                             Y2(indx1+soshift),
     >                             m2)
               sum = sum + m1 + m2
               indx1 = indx1 + npack1
            end do
            hml(indexm)=sum
            indexm=indexm+1
         end do
c     **** restricted and unrestricted dft ****
      else
         if (mb.eq.0) then
            ms1 = 1
            ms2 = ispin
         else
            ms1 = mb
            ms2 = mb
         end if

         indexm=1
         do nb = nb1,nb2
            nbshift = ((nb-nb1))*nemax*npack1
            do ms = ms1,ms2
              msshift = (ms-1)*ne(1)*npack1
              indx1 = 1 + msshift + nbshift
              sum=0.0d0
              do k = 1,ne(ms)
                call Cram_cc_idot(nb,Y1(indx1),Y2(indx1),m1)
                sum = sum + m1
                indx1 = indx1 + npack1
              end do
              hml(indexm)=sum
              indexm=indexm+1
            end do
         end do
      end if
      nsize=indexm-1
      if (np.gt.1) call C3dB_Vector_SumAll(nsize,hml)
      return
      end


c     cccccccccccccccccccccccccccccccccccc
c     *                                  *
c     *            p_index_name          *
c     *                                  *
c     cccccccccccccccccccccccccccccccccccc

      character*7 function p_index_name(i)
      implicit none
      integer i
      integer itmp,j0,j1,j2,j3,j4,j5
      character*7 name
      itmp = i
      j5 = itmp/100000
      itmp = itmp - j5*100000
      j4 = itmp/10000
      itmp = itmp - j4*10000
      j3 = itmp/1000
      itmp = itmp - j3*1000
      j2 = itmp/100
      itmp = itmp - j2*100
      j1 = itmp/10
      itmp = itmp - j1*10
      j0 = itmp/1
      name(1:1) = '_'
      name(2:2) = CHAR(j5+ICHAR('0'))
      name(3:3) = CHAR(j4+ICHAR('0'))
      name(4:4) = CHAR(j3+ICHAR('0'))
      name(5:5) = CHAR(j2+ICHAR('0'))
      name(6:6) = CHAR(j1+ICHAR('0'))
      name(7:7) = CHAR(j0+ICHAR('0'))
      p_index_name = name
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ccccccccccccccccccccccccccccccccccc
c     *                                 *
c     *             prand               *
c     *                                 *
c     ccccccccccccccccccccccccccccccccccc 
c     returns a psuedorandom number   between 0 and 1   
c     Park Miller algo. not that great but ...            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function prand(idum)
      implicit none
      integer idum,ia,im,iq,ir,mask
      real*8 am
      parameter(ia=16807,im=2147483647,AM=1.0D0/IM,IQ=127773,
     >      IR=2836,MASK=123456789)
      integer k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      PRAND=DBLE(IDUM)/DBLE(IM)
      IDUM=IEOR(IDUM,MASK)
      RETURN 
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   assigns the variables for loop over brillioun 
c     points
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine assign_nb(nbb,nb1,nb2)
      implicit none
      integer nbb,nb1,nb2
#include "Pneb.fh"
      if (nbb.eq.0) then
         nb1=1
         nb2=nbrillq
      else
         nb1=nbb
         nb2=nbb
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  assigns the variables for loop over spin states 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine assign_ms(mbb,ms1,ms2)
      implicit none
#include "Pneb.fh"
      integer mbb,ms1,ms2
      if (mbb.eq.0) then
        ms1=1 
        ms2=ispinq
      else
        ms1=mbb
        ms2=mbb
      end if
      return
      end            

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c   Copies the Hermitian Conjugate of A into B
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Pneb_w_copy_dagger(mbb,nbb,A,B)
      implicit none
      integer mbb,nbb
      complex*16 A(*),B(*)

#include "Pneb.fh"

      integer nb,nb1,nb2,ms,ms1,ms2,indx,indext,i,j,xindx

      call assign_nb(nbb,nb1,nb2) 
      call assign_ms(mbb,ms1,ms2)
      indx=1
      do nb=nb1,nb2
        do ms=ms1,ms2
          xindx=indx-1
          do i=1,ne(ms)
            do j=1,ne(ms)
              indext=xindx+i+(j-1)*ne(ms)
              B(indext)=Dconjg(A(indx))
              indx=indx+1
            end do  
          end do
        end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Pneb_www_multiply(mbb,nbb,A,B,C)
      implicit none
      integer mbb,nbb
      complex*16 A(*),B(*),C(*)
#include "Pneb.fh"
      integer nb,nb1,nb2,ms,ms1,ms2
      integer indx  
      complex*16 zero,one
      zero=dcmplx(0.0d0,0.0d0)
      one=dcmplx(1.0d0,0.0d0)
      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      indx=1
      do nb=nb1,nb2
      do ms=ms1,ms2
        call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),one,
     >    A(indx),ne(ms),B(indx),ne(ms),zero,C(indx),ne(ms))
        indx=indx+ne(ms)*ne(ms)               
      end do
      end do
      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     returns a specified value for the w array
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex*16 function Pneb_w_value(nb,ms,i,j,W)
      implicit none
      integer nb,ms,i,j
      complex*16 W(*)
#include "Pneb.fh"
      integer indx
     
      indx=i+(j-1)*ne(ms)
      if (spin_orbit) then
        indx=indx+(nb-1)*ne(1)*ne(1)
      else
        indx=indx+ (ms-1)*ne(1)*ne(1)
     >        + (nb-1)*(ne(1)*ne(1)+ne(2)*ne(2))
      end if
      Pneb_w_value=W(indx)
      return
      end  

c     *******************************************
c     *                                         *
c     *           Pneb_www_Multiply1            *
c     *                                         *
c     *******************************************

c    returns C = beta*C + alpha*A*B

      subroutine Pneb_www_Multiply1(mbb,nbb,alpha,A,B,beta,C)
      implicit none
      integer mbb,nbb
      complex*16 A(*),B(*),C(*)
      complex*16 alpha,beta

#include "bafdecls.fh"
#include "Pneb.fh"

*     **** local variables ****
      integer ms1,ms2,nb1,nb2,indx,ms,nb

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      indx=1
      do nb=nb1,nb2
      do ms=ms1,ms2
         if (ne(ms).gt.0) then
            call ZGEMM('N','N',ne(ms),ne(ms),ne(ms),
     >                 alpha,
     >                 A(indx),ne(ms),
     >                 B(indx),ne(ms),
     >                 beta,
     >                 C(indx),ne(ms))
            indx = indx + ne(ms)*ne(ms)
         end if
      end do
      end do
      return
      end


c     ccccccccccccccccccccccccccccccccccccccccccc
c     *                                         *
c     *           Pneb_www_Multiply2            *
c     *                                         *
c     ccccccccccccccccccccccccccccccccccccccccccc
  
c    returns C = beta*C + alpha*conjg(A')*B
c
      subroutine Pneb_www_Multiply2(mbb,nbb,alpha,A,B,beta,C)
      implicit none
      integer mbb,nbb
      complex*16 A(*),B(*),C(*)
      complex*16 alpha,beta

#include "bafdecls.fh"
#include "Pneb.fh"

*     **** local variables ****
      integer ms1,ms2,nb1,nb2,indx,ms,nb

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      indx=1
      do nb=nb1,nb2
      do ms=ms1,ms2
         if (ne(ms).gt.0) then
            call ZGEMM('C','N',ne(ms),ne(ms),ne(ms),
     >                 alpha,
     >                 A(indx),ne(ms),
     >                 B(indx),ne(ms),
     >                 beta,
     >                 C(indx),ne(ms))
            indx = indx + ne(ms)*ne(ms)
         end if
      end do
      end do
      return
      end


c     ccccccccccccccccccccccccccccccccccccccccccc
c     *                                         *
c     *           Pneb_www_Multiply3            *
c     *                                         *
c     ccccccccccccccccccccccccccccccccccccccccccc

c    returns C = beta*C + alpha*A*conj(B')

      subroutine Pneb_www_Multiply3(mbb,nbb,alpha,A,B,beta,C)
      implicit none
      integer mbb,nbb
      complex*16 A(*),B(*),C(*)
      complex*16 alpha,beta

#include "bafdecls.fh"
#include "Pneb.fh"

*     **** local variables ****
      integer ms1,ms2,nb1,nb2,indx,ms,nb,ishift2,ishift3
      integer shift2

      if (mparallelized) then
         
         if (mbb.eq.0) then
            ms1=1
            ms2=ispinq
            ishift2 = mcq(1)*ncq(1)
         else
            ms1=mbb
            ms2=mbb
            ishift2 = 0 
         end if
         if (nbb.eq.0) then
            nb1=1
            nb2=nbrillq
            ishift3 = neq(1)*neq(1) + neq(2)*neq(2)*(ispinq-1)
         else
            nb1=nbb
            nb2=nbb
            ishift3 = 0
         endif
      

         do nb=nb1,nb2
         do ms=ms1,ms2
            if (ne(ms).le.0) go to 20
            shift2 = 1 + (ms-1)*ishift2 + (nb-nb1)*ishift3
            call CMatrix_zgemm3(ne(ms),ne(ms),ne(ms),64,
     >             1.0d0,
     >             A(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             B(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             0.0d0,
     >             C(shift2),mcq(ms),
     >                       int_mb(mc(1,ms)),
     >                       int_mb(nc(1,ms)),
     >             taskid_i,taskid_j,
     >             np_i,np_j,
     >             comm_i, comm_j,
     >             dcpl_mb(work1(1)),dcpl_mb(work2(1)))

   20       continue
         end do
         end do

      else
         call assign_nb(nbb,nb1,nb2)
         call assign_ms(mbb,ms1,ms2)
         indx=1
         do nb=nb1,nb2
         do ms=ms1,ms2
            if (ne(ms).gt.0) then
               call ZGEMM('N','C',ne(ms),ne(ms),ne(ms),
     >                    alpha,
     >                    A(indx),ne(ms),
     >                    B(indx),ne(ms),
     >                    beta,
     >                    C(indx),ne(ms))
               indx = indx + ne(ms)*ne(ms)
            end if
         end do
         end do
      end if
      return
      end


c     ccccccccccccccccccccccccccccccccccccccccccc
c     *                                         *
c     *               Pneb_w_dmax               *
c     *                                         *
c     ccccccccccccccccccccccccccccccccccccccccccc

      real*8 function Pneb_w_dmax(mbb,nbb,A)
      implicit none
      integer mbb,nbb
      complex*16 A(*)

#include "Pneb.fh"

*     **** local variables ****
      integer ms1,ms2,nb1,nb2,ishift2,ishift3,shift2,nb,ms
      real*8 adiff1,adiff2
      complex*16 wt

*     **** external functions ****
      integer  izamax
      external izamax

      call assign_ms(mbb,ms1,ms2)
      call assign_nb(nbb,nb1,nb2)
      if (mbb.eq.0) then
         ishift2 = ne(1)*ne(1)
      else
         ishift2 = 0
      end if
      if (nbb.eq.0) then
         if (spin_orbit) then
            ishift3 = ne(1)*ne(1)
         else
            ishift3 = ne(1)*ne(1)+ne(2)*ne(2)
         end if
      else
         ishift3 = 0
      end if

      adiff1 = 0.0d0
      adiff2 = 0.0d0
      do nb=nb1,nb2
      do ms=ms1,ms2
         if (ne(ms).le.0) go to 30

         shift2 = 1 + (ms-ms1)*ishift2 + (nb-nb1)*ishift3

         wt     = A(shift2-1+izamax(ne(ms)*ne(ms),A(shift2),1))
         adiff1 = adiff2
         adiff2 = dble(wt)**2 + dimag(wt)**2
         if (adiff2.gt.adiff1) adiff1 = adiff2
  30     continue
      end do
      end do

      Pneb_w_dmax = adiff1
      return
      end

c     *********************************************
c     *                                           *
c     *          Pneb_SCVTrans1                   *
c     *                                           *
c     *********************************************

c     return cos(S*t)V' = A and sin(S*t)*V'=B

      subroutine Pneb_SCVtrans1(mbb,nbb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer    mbb,nbb
      real*8     t
      real*8     S(*)
      complex*16 A(*),B(*),Vt(*)
      real*8     SA(*),SB(*)

#include "Pneb.fh"

c     **** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,j,k,nj
      integer ishift11,ishift12,shift1
      integer ishift21,ishift22,shift2
      integer indx1,indx2


      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      if (mbb.eq.0) then
         if (spin_orbit) then
            nj = ne(1)
            ishift11 = ne(1)
            ishift21 = ne(1)*ne(1)
         else
            nj = ne(1) + ne(2)
            ishift11 = ne(1)
            ishift21 = ne(1)*ne(1)
         end if
      else
        nj = ne(mbb)
        ishift11 = 0
        ishift21 = 0
      end if

      if (nbb.eq.0) then
         nj = nj*nbrillq
         if (spin_orbit) then
            ishift12 = ne(1)
            ishift22 = ne(1)*ne(1)
         else
            ishift12 = ne(1)+ne(2)
            ishift22 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         end if
      else
         ishift12 = 0
         ishift22 = 0
      end if


      do j=1,nj
         SA(j) = dcos(S(j)*t)
         SB(j) = dsin(S(j)*t)
      end do

      do nb=nb1,nb2
      do ms=ms1,ms2

         if (ne(ms).le.0) go to 30
         shift1 = 1 + (ms-ms1)*ishift11 + (nb-nb1)*ishift12
         shift2 = 1 + (ms-ms1)*ishift21 + (nb-nb1)*ishift22

         indx2 = shift2
         do k=1,ne(ms)
            indx1 = shift1
            do j=1,ne(ms)
               A(indx2)=SA(indx1)*Vt(indx2)
               B(indx2)=SB(indx1)*Vt(indx2)
               indx2 = indx2 + 1
               indx1 = indx1 + 1
            end do
         end do

  30     continue
      end do
      end do

      return
      end

c     *********************************************
c     *                                           *
c     *          Pneb_SCVTrans2                   *
c     *                                           *
c     *********************************************

c     returns S*sin(S*t)*V' = A 
c         and S*cos(S*t)*V' = B

      subroutine Pneb_SCVtrans2(mbb,nbb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer    nbb,mbb
      real*8     t
      real*8     S(*)
      complex*16 Vt(*),A(*),B(*)
      real*8     SA(*),SB(*)

#include "Pneb.fh"

c     **** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,j,k,nj
      integer ishift11,ishift12,shift1
      integer ishift21,ishift22,shift2
      integer indx1,indx2
      

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      if (mbb.eq.0) then
         if (spin_orbit) then
            nj = ne(1)
            ishift11 = ne(1)
            ishift21 = ne(1)*ne(1)
         else
            nj = ne(1) + ne(2)
            ishift11 = ne(1)
            ishift21 = ne(1)*ne(1)
         end if
      else
        nj = ne(mbb)
        ishift11 = 0
        ishift21 = 0
      end if

      if (nbb.eq.0) then
         nj = nj*nbrillq
         if (spin_orbit) then
            ishift12 = ne(1)
            ishift22 = ne(1)*ne(1) 
         else
            ishift12 = ne(1)+ne(2)
            ishift22 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         end if
      else
         ishift12 = 0
         ishift22 = 0
      end if

      do j=1,nj
         SA(j) = S(j)*dsin(S(j)*t)
         SB(j) = S(j)*dcos(S(j)*t)
      end do
      
      do nb=nb1,nb2
      do ms=ms1,ms2

         if (ne(ms).le.0) go to 30
         shift1 = 1 + (ms-ms1)*ishift11 + (nb-nb1)*ishift12
         shift2 = 1 + (ms-ms1)*ishift21 + (nb-nb1)*ishift22

         indx2 = shift2
         do k=1,ne(ms)
            indx1 = shift1
            do j=1,ne(ms)
               A(indx2)=SA(indx1)*Vt(indx2) 
               B(indx2)=SB(indx1)*Vt(indx2)
               indx2 = indx2 + 1
               indx1 = indx1 + 1
            end do
         end do

  30     continue
      end do
      end do

      return
      end


c     *********************************************
c     *                                           *
c     *          Pneb_SCVTrans3                   *
c     *                                           *
c     *********************************************

c     Returns (1-cos(S*t))*V' in A and
c              sin(St)V' in Bc

      subroutine Pneb_SCVtrans3(mbb,nbb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer    nbb,mbb
      real*8     t
      real*8     S(*)
      complex*16 Vt(*),A(*),B(*)
      real*8     SA(*),SB(*)

#include "Pneb.fh"

c     **** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,j,k,nj
      integer ishift11,ishift12,shift1
      integer ishift21,ishift22,shift2
      integer indx1,indx2
      real*8 tmp

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      if (mbb.eq.0) then
         nj=ne(1)+ne(2)
         if (spin_orbit) then
            nj = ne(1)
         end if
         ishift11 = ne(1)
         ishift21 = ne(1)*ne(1)
      else
        nj = ne(mbb)
        ishift11 = 0
        ishift21 = 0
      end if

      if (nbb.eq.0) then
         nj = nj*nbrillq
         if (spin_orbit) then
            ishift12 = ne(1)
            ishift22 = ne(1)*ne(1) 
         else
            ishift12 = ne(1)+ne(2)
            ishift22 = ne(1)*ne(1) + ne(2)*ne(2)*(ispinq-1)
         end if
      else
         ishift12 = 0
         ishift22 = 0
      end if


      do j=1,nj
         tmp=S(j)*t
         SA(j) = dsin(tmp)
         SB(j) = 1.0d0-dcos(tmp)
      end do


      do nb=nb1,nb2
      do ms=ms1,ms2

         if (ne(ms).le.0) go to 30
         shift1 = 1 + (ms-ms1)*ishift11 + (nb-nb1)*ishift12
         shift2 = 1 + (ms-ms1)*ishift21 + (nb-nb1)*ishift22

         indx2 = shift2
         do k=1,ne(ms)
            indx1 = shift1
            do j=1,ne(ms)
               A(indx2)=SA(indx1)*Vt(indx2) 
               B(indx2)=SB(indx1)*Vt(indx2)
               indx2 = indx2 + 1
               indx1 = indx1 + 1
            end do
         end do

  30     continue
      end do
      end do

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Return the size of real diagonal s matrix 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function Pneb_s_size(mbb,nbb)
      implicit none
      integer mbb,nbb
#include "Pneb.fh"
      integer bsize
      if (mbb.eq.0) then
        bsize=ne(1)+ne(2)*(ispinq-1)
      else
        bsize=ne(mbb)
      end if
      if (nbb.eq.0) bsize=bsize*nbrillq
      Pneb_s_size=bsize
      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     returns the indices for a block of s sized memory
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function Pneb_s_push_get(mbb,nbb,sxi)
      implicit none
      integer mbb,nbb,sxi(2)
#include "Pneb.fh"
#include "bafdecls.fh"
      integer ss,Pneb_s_size
      external Pneb_s_size  
      ss=Pneb_s_size(mbb,nbb)
      Pneb_s_push_get=BA_push_get(mt_dbl,ss,'Sblk',sxi(2),sxi(1))
      return
      end
cccccccccccccccccc
c     Free some memory
cccccccccccccccccccccc
      logical function Pneb_s_pop_stack(sx)
      implicit none
      integer sx(2)
#include "Pneb.fh"
#include "bafdecls.fh"
      Pneb_s_pop_stack=BA_pop_stack(sx(2))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccc
c     Pneb w expand
c     place the ne by ne matrix, A, in the proper place
c       of "large" matrix B
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Pneb_w_expand(mbb,nbb,A,B)
      implicit none
      integer mbb,nbb
      complex*16 A(*),B(*)
#include "Pneb.fh"
      integer ms,ms1,ms2,bshift,indx2,indx1       
      if (mbb.eq.0) then
         ms1=1
         ms2=ispinq
      else
         ms1=mbb
         ms2=mbb
      end if 
      if (nbb.eq.0) then
        write(*,*)"PROGRAMMER ERROR IN PNEB_W_EXPAND!"
        return
      end if
      bshift=ne(1)*ne(1)+ne(2)*ne(2)*(ispinq-1)
      indx1=1
      do ms=ms1,ms2
        indx2=(nbb-1)*bshift+(ms-ms1)*ne(1)*ne(1)
        call Cram_c_copy(nbb,A(1),B(indx2))
      end do
      return 
      end 	
cccccccccccccccc

*     **********************************************
*     *                                            *
*     *             Pneb_w_scale_s22               *
*     *                                            *
*     **********************************************
      subroutine Pneb_w_scale_s22(mbb,nbb,dte,s22)
      implicit none
      integer mbb,nbb
      real*8 dte
      complex*16 s22(*)

#include "Pneb.fh"

*     **** local variables ****
      integer nb,ms,nb1,nb2,ms1,ms2,ishift2,ishift3,shift2,shift3
      integer j,k,indx,indxt

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      if (mbb.eq.0) then
         ishift2 = ne(1)*ne(1)
      else
         ishift2 = 0
      end if
      if (nbb.eq.0) then
         if (spin_orbit) then
           ishift3 = ne(1)*ne(1)
         else
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)
         end if
      else
         ishift3 = 0
      end if
      do nb=nb1,nb2
        shift3 = (nb-nb1)*ishift3
        do ms=ms1,ms2
          if (ne(ms).le.0) go to 30
          shift2 = (ms1-1)*ishift2 + shift3
          do k=1,ne(ms)
             indx = k + (k-1)*ne(ms) + shift2
             s22(indx) = (1.0d0-s22(indx))*0.5d0/dte

             do j=k+1,ne(ms)
                indx  = j + (k-1)*ne(ms) + shift2
                indxt = k + (j-1)*ne(ms) + shift2
                s22(indx) = -s22(indx)*0.5d0/dte
                s22(indxt) = dconjg(s22(indx))
             end do
          end do

 30     continue
        end do
      end do
      return
      end 




*     **********************************************
*     *                                            *
*     *             Pneb_w_scale_s21               *
*     *                                            *
*     **********************************************
      subroutine Pneb_w_scale_s21(mbb,nbb,dte,s21)
      implicit none
      integer mbb,nbb
      real*8 dte
      complex*16 s21(*)

#include "Pneb.fh"

*     **** local variables ****
      integer nb,ms,nb1,nb2,ms1,ms2,ishift2,ishift3,shift2,shift3
      integer j,k,indx,indxt

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      if (mbb.eq.0) then
         ishift2 = ne(1)*ne(1)
      else
         ishift2 = 0
      end if
      if (nbb.eq.0) then
         if (spin_orbit) then
           ishift3 = ne(1)*ne(1)
         else
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)
         end if
      else
         ishift3 = 0
      end if
      do nb=nb1,nb2
        shift3 = (nb-nb1)*ishift3
        do ms=ms1,ms2
          if (ne(ms).le.0) go to 30
          shift2 = (ms1-1)*ishift2 + shift3
          do k=1,ne(ms)
             indx = k + (k-1)*ne(ms) + shift2
             s21(indx) = (1.0d0-s21(indx))*0.5d0

             do j=k+1,ne(ms)
                indx  = j + (k-1)*ne(ms) + shift2
                indxt = k + (j-1)*ne(ms) + shift2
                s21(indx) = -s21(indx)*0.5d0
                s21(indxt) = dconjg(s21(indx))
             end do
          end do

 30     continue
        end do
      end do
      return
      end

*     **********************************************
*     *                                            *
*     *             Pneb_w_scale_s11               *
*     *                                            *
*     **********************************************
      subroutine Pneb_w_scale_s11(mbb,nbb,dte,s11)
      implicit none
      integer mbb,nbb
      real*8 dte
      complex*16 s11(*)

#include "Pneb.fh"

*     **** local variables ****
      integer nb,ms,nb1,nb2,ms1,ms2,ishift2,ishift3,shift2,shift3
      integer j,k,indx,indxt

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      if (mbb.eq.0) then
         ishift2 = ne(1)*ne(1)
      else
         ishift2 = 0
      end if
      if (nbb.eq.0) then
         if (spin_orbit) then
           ishift3 = ne(1)*ne(1)
         else
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)
         end if
      else
         ishift3 = 0
      end if
      do nb=nb1,nb2
        shift3 = (nb-nb1)*ishift3
        do ms=ms1,ms2
          if (ne(ms).le.0) go to 30
          shift2 = (ms1-1)*ishift2 + shift3
          do k=1,ne(ms)
             indx = k + (k-1)*ne(ms) + shift2
             s11(indx) = -s11(indx)*0.5d0*dte

             do j=k+1,ne(ms)
                indx  = j + (k-1)*ne(ms) + shift2
                indxt = k + (j-1)*ne(ms) + shift2
                s11(indx) = -s11(indx)*0.5d0*dte
                s11(indxt) = dconjg(s11(indx))
             end do
          end do

 30     continue
        end do
      end do
      return
      end


*     **********************************************
*     *                                            *
*     *             Pneb_w_trace               *
*     *                                            *
*     **********************************************
      subroutine Pneb_w_trace(mbb,nbb,hml,val)
      implicit none
      integer mbb,nbb
      complex*16 hml(*)
      real*8 val

#include "Pneb.fh"

*     **** local variables ****
      integer nb,ms,nb1,nb2,ms1,ms2,ishift2,ishift3,shift2,shift3
      integer j,k,indx,indxt

      val = 0.0d0
      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      if (mbb.eq.0) then
         ishift2 = ne(1)*ne(1)
      else
         ishift2 = 0
      end if
      if (nbb.eq.0) then
         if (spin_orbit) then
           ishift3 = ne(1)*ne(1)
         else
           ishift3 = ne(1)*ne(1) + ne(2)*ne(2)
         end if
      else
         ishift3 = 0
      end if
      do nb=nb1,nb2
        shift3 = (nb-nb1)*ishift3
        do ms=ms1,ms2
          if (ne(ms).le.0) go to 30
          shift2 = (ms1-1)*ishift2 + shift3
          do k=1,ne(ms)
             indx = k + (k-1)*ne(ms) + shift2
             val = val + dble(hml(indx))
          end do
 30     continue
        end do
      end do
      return
      end



*     ****************************
*     *                          *
*     *     Pneb_f_SOSpins_tag   *
*     *                          *
*     ****************************

*   Returns the  Spin Amplitudes for eigenvalue printout.
*  
*   This routine is done separately for c each BZ point?
*
      subroutine Pneb_f_SOSpins_tag(Y_tag,svec_tag)
      implicit none
      integer Y_tag,svec_tag

#include "bafdecls.fh"
#include "stdio.fh"
#include "Pneb.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      integer nbq,npack1,Y_shift,svec_shift

*     **** external functions ****
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      external cpsi_data_get_chnk,cpsi_data_nsize


*     **** simple error checking - make sure we are doing a spin-orbit calculation ***
      if (.not.spin_orbit) then
         call Parallel_taskid(taskid)
         if (taskid.eq.MASTER)
     >      write(luout,*)
     >      'Warning: Pneb_f_SOSpins_tag called when spin-orbit not set'
         return
      end if

      npack1 = cpsi_data_nsize(Y_tag)/2 !*** divide by 2 for complex ***

      call cpsi_data_update(svec_tag)
      do nbq=1,nbrillq
         Y_shift    = cpsi_data_get_chnk(Y_tag,   nbq)
         svec_shift = cpsi_data_get_chnk(svec_tag,nbq)

         call Pneb_f_SOSpins_sub(nbq,ne,npack1,
     >                           dbl_mb(Y_shift),
     >                           dbl_mb(svec_shift))
      end do
      call cpsi_data_noupdate(svec_tag)
      return
      end

      subroutine Pneb_f_SOSpins_sub(nbq,neq,npack1,Y,S)
      implicit none
      integer nbq,neq(2),npack1
      complex*16 Y(*)
      real*8     S(*)

*     **** local variables ****
      integer k,indx,indx2,sindx
      real*8     s1,s2
      complex*16 z1

      indx  = 1
      indx2 = neq(1)*npack1
      sindx = 1
      do k=1,neq(1)
         call Cram_cc_dot(nbq, Y(indx), Y(indx), s1)
         call Cram_cc_dot(nbq, Y(indx2),Y(indx2),s2)
         call Cram_cc_zdot(nbq,Y(indx), Y(indx2),z1)
         S(sindx)   = dble(z1)
         S(sindx+1) = dimag(z1)
         S(sindx+2) = 0.50d0*(s1-s2)
         indx  = indx  + npack1
         indx2 = indx2 + npack1
         sindx = sindx + 3
      end do
      return
      end

      subroutine Pneb_w_gather(nall,npack,indx,A,B)
      implicit none
      integer nall,npack,indx(*)
      complex*16 A(*),B(*)
      integer i
      call Parallel_shared_vector_zero(.true.,2*nall,B)
!$OMP DO
      do i=1,npack
         B(indx(i)) = A(i)
      end do
!$OMP END DO
      call C3dB_Vector_SumAll(2*nall,B)
      call C1dB_Vector_SumAll(2*nall,B)
      return
      end

      subroutine Pneb_w_scatter(npack,indx,A,B)
      implicit none
      integer npack,indx(*)
      complex*16 A(*),B(*)
      integer i
!$OMP DO
      do i=1,npack
         B(i) = A(indx(i))
      end do
!$OMP END DO
      return
      end


      subroutine Pneb_m_gather(nall,npack,indx,A,B)
      implicit none
      integer nall,npack,indx(*)
      real*8 A(*),B(*)
      integer i
      !call dcopy(nall,0.0d0,0,B,1)
      call Parallel_shared_vector_zero(.true.,nall,B)
!$OMP DO
      do i=1,npack
         B(indx(i)) = A(i)
      end do
!$OMP END DO
      !call Parallel_Vector_SumAll(nall,B)
      call C3dB_Vector_SumAll(nall,B)
      call C1dB_Vector_SumAll(nall,B)
      return
      end

      subroutine Pneb_m_scatter(npack,indx,A,B)
      implicit none
      integer npack,indx(*)
      real*8 A(*),B(*)
      integer i
!$OMP DO
      do i=1,npack
         B(i) = A(indx(i))
      end do
!$OMP END DO
      return
      end




