*
* $Id: Pneb.F,v 1.5 2006-12-12 20:12:42 d3p708 Exp $
*

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

*     ***********************************
*     *					*
*     *	        Pneb_init		*	
*     *					*
*     ***********************************

      subroutine Pneb_init(ispin_in,ne_in,nbrill_in,spin_orbit_in)
      implicit none
      integer ispin_in
      integer ne_in(2)
      integer nbrill_in
      logical spin_orbit_in

#include "Pneb.fh"

      ispin   = ispin_in
      ne(1)   = ne_in(1)
      ne(2)   = ne_in(2)
      nbrill  = nbrill_in
      spin_orbit = spin_orbit_in
  
      neq(1)  = ne(1)
      neq(2)  = ne(2)
      nbrillq = nbrill
      if (spin_orbit) then
         ispinq = 1
      else
         ispinq = ispin
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

#include "Pneb.fh"

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
        nbsize = nbrill
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
        nbsize = nbrill
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

c     ****************************************
c     *                                      *
c     *        Pneb_m_allocate               *
c     *                                      *
c     ****************************************

      logical function Pneb_m_allocate(mb,nb,hml)
      implicit none
      integer mb,nb
      integer hml(2)

#include "mafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size
     
      size = Pneb_m_size(mb,nb)
      Pneb_m_allocate = MA_alloc_get(mt_dbl,size,'hml',hml(2),hml(1))

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

#include "mafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size
     
      size = Pneb_w_size(mb,nb)
      Pneb_w_allocate = MA_alloc_get(mt_dcpl,size,'hml',hml(2),hml(1))
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

#include "mafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size
     
      size = nblocks*Pneb_m_size(mb,nb)
     
      Pneb_m_allocate_block =
     >   MA_alloc_get(mt_dbl,size,'hml',hml(2),hml(1))
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

#include "mafdecls.fh"

*     **** local variables ****
      integer size

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size
     
      size = nblocks*Pneb_w_size(mb,nb)
     
      Pneb_w_allocate_block =
     >   MA_alloc_get(mt_dcpl,size,'hml',hml(2),hml(1))
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

#include "mafdecls.fh"

      Pneb_m_free = MA_free_heap(hml(2))
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

#include "mafdecls.fh"

      Pneb_w_free = MA_free_heap(hml(2))
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
#include "mafdecls.fh"

*     **** local variables
      integer size

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size

      size = Pneb_m_size(mb,nb)
      Pneb_m_push_get = MA_push_get(mt_dbl,size,'hml',hml(2),hml(1))
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
#include "mafdecls.fh"

*     **** local variables
      integer size

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      size = Pneb_w_size(mb,nb)
      Pneb_w_push_get = MA_push_get(mt_dcpl,size,'hml',hml(2),hml(1))
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

#include "mafdecls.fh"
#include "Pneb.fh"

*     **** local variables
      integer nsize

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size

      nsize = nblock*Pneb_m_size(mb,nb)
      Pneb_m_push_get_block = 
     >    MA_push_get(mt_dbl,nsize,'hml',hml(2),hml(1))
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

#include "mafdecls.fh"
#include "Pneb.fh"

*     **** local variables
      integer nsize

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size

      nsize = nblock*Pneb_w_size(mb,nb)
      Pneb_w_push_get_block = 
     >    MA_push_get(mt_dcpl,nsize,'hml',hml(2),hml(1))
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

#include "mafdecls.fh"

      Pneb_m_pop_stack = MA_pop_stack(hml(2))
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

#include "mafdecls.fh"

      Pneb_w_pop_stack = MA_pop_stack(hml(2))
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

#include "Pneb.fh"

*     **** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,nsize,n,np
      integer shift,shiftso,shift2,ishift2,ishift3

*     **** external functions ****
      integer  Pneb_m_size
      external Pneb_m_size


      call Parallel_np(np)
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
        ishift3 = ne(1)*ne(1) + ne(2)*ne(2)
      else
        nb1 = nbb
        nb2 = nbb
        ishift3 = 0
      end if

      do nb=nb1,nb2
         do ms=ms1,ms2
           shift  = 1 + (ms-1)*ne(1)*npack1 +(nb-1)*(ne(1)+ne(2))*npack1
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

#include "Pneb.fh"

*     **** local variables ****
      complex*16 zero,one
      parameter (zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0))

      integer nb,nb1,nb2,ms,ms1,ms2,nsize,n,np
      integer shift,shiftso,shift2,ishift2,ishift3

*     **** external functions ****
      integer  Pneb_w_size
      external Pneb_w_size


      call Parallel_np(np)
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
        ishift3 = ne(1)*ne(1) + ne(2)*ne(2)
      else
        nb1 = nbb
        nb2 = nbb
        ishift3 = 0
      end if

      do nb=nb1,nb2
         do ms=ms1,ms2
           shift  = 1 + (ms-1)*ne(1)*npack1 +(nb-1)*(ne(1)+ne(2))*npack1
           shiftso = shift + ne(1)*npack1
           shift2 = 1 + (ms-1)*ishift2      +(nb-nb1)*ishift3
           n = ne(ms)
           if (n.eq.0) go to 30

              call Cram_ccm_izgemm(nb,n,n,Y1(shift),Y2(shift),
     >                            one,zero,hml(shift2))

              if (spin_orbit) then
                 call Cram_ccm_izgemm(nb,n,n,Y1(shiftso),Y2(shiftso),
     >                               one,one,hml(shift2))
              end if

 30        continue
         end do
      end do

      if (np.gt.1) call C3dB_Vector_SumAll(2*nsize,hml)
      return
      end


      
      
c     ****************************************
c     *                                      *
c     *          Pneb_fwf_multiply           *
c     *                                      *
c     ****************************************

      subroutine Pneb_fwf_multiply(mb,nbb,Y1,hml,npack1,Y2)
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
        ishift3 = ne(1)*ne(1) + ne(2)*ne(2)
      else
        nb1 = nbb
        nb2 = nbb
        ishift3 = 0
      end if


      do nb=nb1,nb2
         do ms=ms1,ms2           
            shift = 1 + (ms-1)*ne(1)*npack1 +(nb-1)*(ne(1)+ne(2))*npack1
            shiftso = shift + ne(1)*npack1
            shift2  = 1 + (ms-1)*ishift2      +(nb-nb1)*ishift3
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

c     *** local variables ***
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
      lwork=3*ne(1)*ne(1)
      value =MA_push_get(mt_dbl,lwork,'rwork',rwork(2),rwork(1))
      value =value.and.MA_push_get(mt_dcpl,lwork,'work',work(2),work(1))
      if (.not.value) call errquit('Pneb_w_diag: out of stack',0,MA_ERR)

      sindex = 1
      vindex = 1
      do nb=nb1,nb2
      do ms=ms1,ms2
        call ZHEEV('V','L',ne(ms),
     >              V(vindex),ne(ms),
     >              S(sindex),
     >              dcpl_mb(work(1)),lwork,
     >              dbl_mb(rwork(1)),
     >              ierr)
        if (ierr.ne.0) then
          call errquit("Pneb_w_diag: ZHEEV failed",ierr,0)
        end if

        vindex=vindex + ne(ms)*ne(ms)
        sindex=sindex + ne(ms)       
      end do
      end do
             

c     ***** pop stack memory *****
      value=          MA_pop_stack(work(2))
      value=value.and.MA_pop_stack(rwork(2))
      if (.not.value) then
        call errquit('Pneb_w_diag:error popping stack',0,MA_ERR)
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "Pneb.fh"

c     ***** local variables ******
      integer soshift,shift1,shift2,n,neall
      integer ms,ms1,ms2,nb,nb1,nb2
      real*8  tsum
      integer tmp(2)
      integer i,j

      if (mbb.eq.0) then
         ms1 = 1
         ms2 = ispinq
         neall = ne(1) + ne(2)
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

c     ***  compute C = Y'*Y ( ne x ne )= V*S2*V' ****
      call Pneb_ffw_Multiply(mbb,nbb,Y,Y,npack1,V)


c     ***  compute S2 and V from C  **** 
      call Pneb_w_diag(mbb,nbb,S,V)


c     **** sqrt(S*S) ****
      do n=1,neall
         if (S(n).lt.0.0d0) S(n) = dabs(S(n))
         S(n) = dsqrt(S(n))
      end do

c     *** compute U*S=Y*V ****
      call Pneb_fwf_multiply(mbb,nbb,Y,V,npack1,U)

*     **** normalize U*S ****
      if (.not.MA_push_get(mt_dbl,neq(1),'tmp',tmp(2),tmp(1)))
     > call errquit('Dneb_f_SVD: out of stack memory',0,MA_ERR)

      soshift = neq(1)*npack1
      do nb=nb1,nb2
         do ms=ms1,ms2
            if (neq(ms).le.0) go to 30
               shift1 = 1 + (ms-1)* neq(1)        *npack1 
     >                    + (nb-1)*(neq(1)+neq(2))*npack1
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
                  call dscal(2*npack1,
     >                       dbl_mb(tmp(1)+n-1),
     >                       U(shift2),1)
                  if (spin_orbit) then
                     call dscal(2*npack1,
     >                          dbl_mb(tmp(1)+n-1),
     >                          U(shift2+soshift),1)
                  end if
                  shift2 = shift2 + npack1
               end do
                 
 30         continue
         end do
      end do

      if (.not.MA_pop_stack(tmp(2)))
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

#include "errquit.fh"
#include "Pneb.fh"

c     ** local variables ****
      integer nb,nb1,nb2,ms,ms1,ms2,index1,index2,index3,index4
      integer j,k,nemax
      integer bindex,nbshift,soshift
      complex*16 wc,wc1,wc2       
      real*8 w,w1,w2

      if (nbb.eq.0) then
        nb1=1
        nb2=nbrill
      else
        nb1=nb
        nb2=nb
      end if
      if (spin_orbit) then
        nemax=ne(1)
        soshift=npack1*nemax
        nbshift=2*soshift
        bindex=1
        do nb=nb1,nb2
         bindex=(nb-1)*nbshift
         do k=1,nemax
           index1=bindex+(k-1)*npack1+1
           call Cram_cc_dot(nb,A(index1),A(index1),w1)
           index2=index1+soshift
           call Cram_cc_dot(nb,A(index2),A(index2),w2)
           w=w1+w2
           if (w.gt.1.d-9) then
             w=1.d0/dsqrt(w)
           else
             write(0,*)"Pneb SO make ortho"
             write(0,*)"Pneb_makeOrtho::this should not happen"
             write(0,*)"Pneb_makeOrtho: < psi | psi > = 0 K= ",k,w     
             w=1.d0/dsqrt(2.0d0*w1)
             call Cram_c_copy(nb,A(index1),A(index2))
           end if
           call Cram_c_SMul(nb,w,A(index1),A(index1))
           call Cram_c_SMul(nb,w,A(index2),A(index2))
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
          bindex=(nb-1)*nbshift+(ms-1)*ne(1)*npack1
          do k=1,ne(ms)
           index1=bindex+(k-1)*npack1+1
           call Cram_cc_dot(nb,A(index1),A(index1),w)
           if (w.gt.1.d-9) then
             w=1.d0/dsqrt(w)
           else
            write(0,*)"Pneb_makeOrtho::this should not happen"
            write(0,*)"Pneb_makeOrtho: < psi | psi > = 0 K= ",k,w     
            w=1.0d0
           end if
           call Cram_c_SMul(nb,w,A(index1),A(index1))
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "Pneb.fh"

c     ccccc local variables ccccc
      integer taskid,MASTER
      PARAMETER(MASTER=0)

      logical oprint
      integer asize,indx,i
      integer nb,nb1,nb2,ms,ms1,ms2,nemax,sm2(2),sm3(2)
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
         nb2=nbrill
         asize=2*nbrill
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
      value =           MA_push_get(mt_dbl,asize,'sm2',sm2(2),sm2(1))
      value = value.and.MA_push_get(mt_dbl,asize,'sm3',sm3(2),sm3(1))
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
            if (dabs(sum2-sum1).gt.1.d-10) gram= .true.
         end do 
      else
         indx=sm2(1)
         do nb=nb1,nb2
         do ms=ms1,ms2
            sum1=dble(ne(ms))
            sum2=dbl_mb(indx)
            indx=indx+1         
            if (dabs(sum2-sum1).gt.1.d-10) gram = .true.
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
      value = MA_pop_stack(sm3(2))
      value = value.and.MA_pop_stack(sm2(2))
      if (.not. value)
     > call errquit('pneb_orthocheckmake:popping stack memory',0,0)

 1400 FORMAT('Warning - Orthogonalization performed, spin,zone:',I2,I4)
 1401 FORMAT(8x,'- error(before)=', E14.6,' (',E14.6,')')
 1402 FORMAT(8x,'- error(after)= ', E14.6,' (',E14.6,')')
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

      integer taskid,MASTER,last_ix
      parameter(MASTER=0)
      integer nb,nb1,nb2,ms,ms1,ms2,j,k,bindex
      real*8 error
      integer index11,index21,index22,index12,bshift,soshift,nemax
      complex*16 w,wc1,wc2
      character*255 full_filename
      character*7 p_index_name
      external p_index_name

c     ccccccccc first line ccccccccc
      if (nbb.eq.0) then
         nb1=1
         nb2=nbrill
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
          bindex=bshift*(nb-1)
          do k=1,nemax
            index11=1+((k-1)+(nb-1)*nemax*2)*npack1
            index12=index11+soshift
            do j=1,nemax
              index21=1+((j-1)+(nb-1)*nemax*2)*npack1
              index22=index21+soshift
              call Cram_cc_zdot(nb,A(index11),A(index21),wc1)
              call Cram_cc_zdot(nb,A(index12),A(index22),wc2)
              w=wc1+wc2
              if (k.eq.j) then
                error=dabs(1.-dble(w))+dabs(imag(w))
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
          bindex=bshift*(nb-1)+(ms-1)*npack1*ne(1)
          do k=1,ne(ms)
            index11=1+bindex+(k-1)*npack1
            do j=1,ne(ms)
              index21=1+bindex+(j-1)*npack1
              call Cram_cc_zdot(nb,A(index11),A(index21),w)
              if (k.eq.j) then
                 error=dabs(1.-dble(w))+dabs(imag(w))
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
c     *          Pneb_ff_traceall              *
c     *                                      *
c     ****************************************

      subroutine Pneb_ff_traceall(mb,nbb,Y1,Y2,npack1,sum)
      implicit none
#include "errquit.fh"
      integer mb,nbb
      complex*16 Y1(*),Y2(*)
      real*8 sum
      integer npack1

#include "Pneb.fh"
#include "mafdecls.fh"

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
        nb2 = nbrill
      else
        nb1 = nbb
        nb2 = nbb
      end if

      value=MA_push_get(mt_dbl,(2*nbrill),'tmp2',tmp2(2),tmp2(1))
      if (.not.value) then
       call errquit('out of stack memory',0,MA_ERR)
      end if

      call Pneb_ff_Trace(mb,nbb,Y1,Y2,npack1,dbl_mb(tmp2(1)))

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
      value=MA_pop_stack(tmp2(2))
      if (.not.value) then
        call errquit('error popping the stack',0,MA_ERR)
      end if

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
      integer j,k,indx1,indx2,indexm
      integer nbshift,msshift,soshift
      integer np
      real*8 m1,m2,sum

      call Parallel_np(np)

      nemax = ne(1)+ne(2)

      if (nbb.eq.0) then
        nb1 = 1
        nb2 = nbrill
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
            nbshift = (nb-1)*2*soshift
            indx2 = 1 + nbshift
            sum=0.0d0
            do k=1,ne(1)
               indx1 = 1 + nbshift
               do j=1,ne(1)
                  call Cram_cc_idot(nb,Y1(indx1),
     >                                 Y2(indx2),
     >                                 m1)
                  call Cram_cc_idot(nb,Y1(indx1+soshift),
     >                                 Y2(indx2+soshift),
     >                                 m2)
                  sum = sum + m1 + m2
                  call flush(6)

                  indx1 = indx1 + npack1
               end do
               indx2 = indx2 + npack1
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
            nbshift = (nb-1)*nemax*npack1
            do ms = ms1,ms2
               msshift = (ms-1)*ne(1)*npack1
               indx2 = 1 + msshift + nbshift
               sum=0.0d0
               do k = 1,ne(ms)
                  indx1 = 1 + msshift + nbshift
                  do j = 1,ne(ms)
                     call Cram_cc_idot(nb,Y1(indx1),
     >                                    Y2(indx2),
     >                                    m1)
                     sum=sum+m1
                     indx1 = indx1 + npack1
                  end do
                  indx2 = indx2 + npack1
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     return cos(S*t)V' = A and sin(S*t)*V'=B
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Pneb_SCVtrans1(t,mbb,nbb,S,Vt,A,B)
      implicit none
      real*8 t,S(*)
      integer mbb,nbb
      complex*16 A(*),B(*),Vt(*)
#include "Pneb.fh"
      real*8 s1,s2
      integer nb,nb1,nb2,ms,ms1,ms2,j,k,sindx,vindx

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      sindx=1
      vindx=1
      do nb=nb1,nb2
        do ms=ms1,ms2
          do j=1,ne(ms)
            s1=dcos(S(sindx)*t)
            s2=dsin(S(sindx)*t)
            sindx=sindx+1
            do k=1,ne(ms)
              A(vindx)=s1*Vt(vindx) 
              B(vindx)=s2*Vt(vindx)
              vindx=vindx+1
            end do
          end do
        end do
      end do
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     returns S*sin(S*t)*V' = A 
c         and S*cos(S*t)*V' = B
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Pneb_SCVtrans2(t,mbb,nbb,S,Vt,A,B)
      implicit none
      real*8 t,S(*)
      integer nbb,mbb
      complex*16 A(*),B(*),Vt(*)
#include "Pneb.fh"
      real*8 s1,s2
      integer nb,nb1,nb2,ms,ms1,ms2,j,k,sindx,vindx

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      vindx=1
      sindx=1
      do nb=nb1,nb2
        do ms=ms1,ms2
          do j=1,ne(ms)
            s2=S(sindx)
            s1=s2*dsin(s2*t)
            s2=s2*dcos(s2*t)
            sindx=sindx+1
            do k=1,ne(ms)
              A(vindx)=s1*Vt(vindx) 
              B(vindx)=s2*Vt(vindx)
              vindx=vindx+1
            end do
          end do
        end do
      end do
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
c     Returns (1-cos(S*t))*V' in A and
c              sin(St)V' in Bc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Pneb_SCVtrans3(t,mbb,nbb,S,Vt,A,B)
      implicit none
      real*8 t,S(*)
      integer mbb,nbb
      complex*16 A(*),B(*),Vt(*)
#include "Pneb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      real*8 s1,s2
      integer nb,nb1,nb2,ms,ms1,ms2,j,k,sindx,vindx

      call assign_nb(nbb,nb1,nb2)
      call assign_ms(mbb,ms1,ms2)
      sindx=1
      vindx=1
      do nb=nb1,nb2
        do ms=ms1,ms2
          do j=1,ne(ms)
            s1=dsin(S(sindx)*t)
            s2=dcos(S(sindx)*t)
            sindx=sindx+1
            do k=1,ne(ms)
              A(vindx)=s1*Vt(vindx) 
              B(vindx)=Vt(vindx)-s2*Vt(vindx)
              vindx=vindx+1
            end do
          end do
       end do
      end do
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
          xindx=indx
          do i=1,ne(ms)
            do j=1,ne(ms)
              indext=xindx+i+(j-1)*ne(ms)
              B(indext)=Dconjg(A(indx))
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
      call assign_mb(mbb,ms1,ms2)
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
cccccccccccccccccccc


                  

 
 
        








                        
