      subroutine qshell_sort
c
c$Id$
c
c
c     Sort quadrature shells, most work to least, based on radius 
c     (assuming largest radius has most angular pieces).
c
      implicit none
#include "cdft.fh"
      integer temp1, temp2, temp3, temp4
      integer i, hsize
c
c     build heap
c
      do i = nqshells/2, 1, -1
         call heapify(nqshells, i)
      enddo
c
c     main part of sort algorithm
c
      hsize = nqshells
      do i    = nqshells, 2, -1
c
c        swap element i and 1
c
         temp1 = iqshell(1,1)
         temp2 = iqshell(2,1)
         temp3 = iqshell(3,1)
         temp4 = iqshell(4,1)
c
         iqshell(1,1) = iqshell(1,i)
         iqshell(2,1) = iqshell(2,i)
         iqshell(3,1) = iqshell(3,i)
         iqshell(4,1) = iqshell(4,i)
c
         iqshell(1,i) = temp1
         iqshell(2,i) = temp2
         iqshell(3,i) = temp3
         iqshell(4,i) = temp4
c
c        maintain heap property from element 1 down
c
         hsize= hsize - 1
         call heapify(hsize, 1)
c
      enddo
      return
      end 
      subroutine heapify(n, elem)
c
c     establish heap property for a tree branch rooted at elem 
c
      implicit none
c
#include "cdft.fh"
c
      integer n, elem 
      integer left, right, smallest, i
      integer temp1, temp2, temp3, temp4
      integer ictr_left, irsh_left, ictr_i, irsh_i
      double precision rpts_left, rpts_i
      integer ictr_right, irsh_right, ictr_smallest, irsh_smallest
      double precision rpts_right, rpts_smallest
c
      i = elem
c
c     Main Loop
c
100   continue 
        left  = 2*i
        right = 2*i + 1
        if (left. gt. n .and. right .gt. n) return   !we traversed entire branch
c
c       check heap property among element i and its children
c
        ictr_left = iatype(iqshell(3,left))
        irsh_left = iqshell(1,left)
        rpts_left = rpts(irsh_left,ictr_left)
        ictr_i = iatype(iqshell(3,i))
        irsh_i = iqshell(1,i)
        rpts_i = rpts(irsh_i,ictr_i)
        if (left .le. n .and. rpts_left .lt. rpts_i) then
           smallest = left
        else
           smallest = i
        endif
        ictr_right = iatype(iqshell(3,right))
        irsh_right = iqshell(1,right)
        rpts_right = rpts(irsh_right,ictr_right)
        ictr_smallest = iatype(iqshell(3,smallest))
        irsh_smallest = iqshell(1,smallest)
        rpts_smallest = rpts(irsh_smallest,ictr_smallest)
        if (right .le. n .and. rpts_right .lt. rpts_smallest) 
     &     smallest = right
c
        if (smallest .ne. i) then
c
c          swap array elements if smallest is not i
c
           temp1 = iqshell(1,i)
           temp2 = iqshell(2,i)
           temp3 = iqshell(3,i)
           temp4 = iqshell(4,i)
c
           iqshell(1,i) = iqshell(1,smallest)
           iqshell(2,i) = iqshell(2,smallest)
           iqshell(3,i) = iqshell(3,smallest)
           iqshell(4,i) = iqshell(4,smallest)
c
           iqshell(1,smallest) = temp1
           iqshell(2,smallest) = temp2
           iqshell(3,smallest) = temp3
           iqshell(4,smallest) = temp4
c
c          traverse down the tree
c
           i = smallest
      goto 100
      endif
      return
      end
      Subroutine mbf_ao_max(rtdb,  rq0, zprim, coord)
c
C$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "bas.fh"
c
      integer rtdb
c
c     Cartesian Coordinates of Integration Center
c
      double precision coord(3,ncenters)
c
c     Compute the quadrature points for a given 
c     set of radial shells.
c
      integer irsh, iang, iqsh,  l, ia_ictr
      double precision r
c
c     Distance Squared between Sampling Points and Centers
c
      double precision rq0(ncenters)
c
      double precision zprim(nbf_ao_mxprim)
      integer ncontrset, n1, icset, ictr, itype, nprimo, ncontr, 
     &        isphere, nbf_ang, iprimo
      double precision zmin
      double precision acc_AO_gauss
      integer mbf, nq, mbfnq
      integer m, me,nproc
      integer avail
      double precision x, y, z, r2
      integer nxtask,nt1,nt2
      external nxtask
c
c     determine node id
c
      me = ga_nodeid()
      nproc=ga_nnodes()
c
c
c     Loop over all the radial shells.
c
      acc_AO_gauss = dble(iAOacc)
      max_mbf = 0
      max_pr_mbf = 0
      max_pr_nq = 0
      max_pr_mbfnq = 0
      nt1 = 0
      nt2 = nxtask(nproc,1)
      do 70 iqsh = 1, nqshells
        if(nt1.eq.nt2) then
c
         irsh = iqshell(1,iqsh)
         ictr = iqshell(3,iqsh)
         iang = iqshell(4,iqsh)
c
         ia_ictr = iatype(ictr)
         r = rpts(irsh,ia_ictr)
c         
         nq = 0
c
         if (leb) then
c         
           nq=nq+ntheta(iang)
         else
           nq=nq+ntheta(iang)*nphi(iang)
c
         endif
         do 40 m = 1, ncenters
               x = coord(1,ictr) - coord(1,m)
               y = coord(2,ictr) - coord(2,m)
               z = coord(3,ictr) - coord(3,m)
               r2 = sqrt(x*x + y*y + z*z)
               rq0(m)=(r2-r)*(r2-r)
   40    continue
c
         if (.not.bas_numcont(AO_bas_han, ncontrset))
     &      call errquit('Exiting in mbf_ao_max.',1, BASIS_ERR)
c
         n1 = 0
c
         do 60 icset = 1,ncontrset
            if (.not.bas_cn2ce(AO_bas_han, icset, ictr))
     &         call errquit('Exiting in mbf_ao_max.',2, BASIS_ERR)
c         
c           get info about current contraction set
c         
            if (.not.bas_continfo(AO_bas_han, icset,
     &         itype, nprimo, ncontr, isphere))
     &         call errquit('Exiting in mbf_ao_max.',3, BASIS_ERR)
c         
c           angular momentum
c         
            l = 0
            if (itype .lt. 0)then
#if 0
               call errquit('mbf_ao_max: sp-type orbital not coded', 5,
     &       BASIS_ERR)
#else
               l=1
#endif
            else 
               l = itype
            endif 
c         
c           cartesian/spherical harmonic
c         
            nbf_ang = 0
            if (isphere .eq. 0)then !  cartesian set
               nbf_ang = (l+1)*(l+2)/2
            elseif (isphere .eq. 1)then !  spherical harmonic
               nbf_ang = 2*l+1
            else 
               call errquit('mbf_ao_max: illegal isphere value', 6,
     &       BASIS_ERR)
            endif 
c         
c           get exponents and contraction coefficients for this contraction set
c         
            if (.not.bas_get_exponent(AO_bas_han, icset, zprim))
     &         call errquit('Exiting in mbf_ao_max.',7, BASIS_ERR)
c   
c           Determine the minimum Gaussian exponent.
c   
            zmin = 1.D+06
            do 50 iprimo = 1,nprimo
               zmin = min(zprim(iprimo),zmin)
   50       continue
c
c           Only include those basis functions that are "non-zero" for at least one 
c           point in the sampling set.
c
            if (zmin*rq0(ictr).gt.acc_AO_gauss)goto 60
            if (l.eq.0)then
c
c              =============>  S Contractions  <============= 
c
               n1 = n1 + ncontr
            elseif (l.eq.1)then
c
c              =============>  P Contractions  <============= 
c
               n1 = n1+ncontr*3
            elseif (l.eq.2)then
c
c              =============>  D Contractions  <============= 
c
               n1 = n1 + ncontr*6
            else
c
c              =============>  General Case  <============= 
c
               n1 = n1 + ncontr*nbf_ang
            endif
c
  60     continue
c
         mbf = n1
c
c        need to determine max_mbf and max(mbf,nq) pair
c
         if (mbf.gt.max_mbf)then
            max_mbf = mbf
         endif
         mbfnq = mbf*nq
         if (mbfnq.gt.max_pr_mbfnq)then
            max_pr_mbfnq = mbfnq
            max_pr_mbf = mbf
            max_pr_nq = nq
         endif

c
            nt1 = nt1 + 1
            nt2 = nxtask(nproc,1)
         else
            nt1 = nt1 + 1
         endif
c
c
   70 continue
      nt2 = nxtask(-nproc,1)
      call ga_igop(dft_mxmbfnq, max_pr_mbfnq, 1, 'max')
      call ga_igop(dft_mxmbf, max_mbf, 1, 'max')
c
c     This can be further reduced by blocking the computed 
c     angular grid.
c
c     Assume we want to use no more than 1/3 of physical memory (stack + heap) for
c     the quadrature.  So, lets put chi(nq,mbf) and delchi(nq,3,mbf)
c     in roughly 1/3 - 8Mb by chunking up nq.
c
c     find - (minimum)amount local available memory on all nodes
c     
      call ga_sync
      avail = MA_inquire_avail(mt_dbl)
      call ga_igop(msg_min_stack_avail, avail, 1, 'min')

c      write(6,*)' avail = ',avail
      avail = avail/3
c      write(6,*)' avail/3 = ',avail
      avail = avail - 1024*1024
      if(avail.lt.0)
     &         call errquit('xc_setquad: out of memory',avail, MEM_ERR)
c      write(6,*)' amt to be used for xc = ',avail
c     
      if ( (4*max_pr_mbfnq) .gt. avail )then
         nq_chunk = avail/4/max_mbf
c
c        redefine max_pr_mbfnq
c
         max_pr_mbfnq = max_mbf*nq_chunk
c
c        reset store_wght to false (not working yet)
c
         store_wght=.false.
      else
c
c        everything fits no chunking necessary
c
         nq_chunk = 0
c
c        redefine nq_task
c
Cedo         if(nquad_task.eq.1) then
          if (.not. rtdb_get(rtdb, 'dft:nquad_task', mt_int, 1,
     &        nquad_task))then
          nquad_task = min(avail/2/(4*max_pr_mbfnq),6)
          if(nquad_task.lt.1) nquad_task=1
           nqmax=nqmax*nquad_task
           if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1, 
     &          nquad_task))
     &          call errquit('mbf_ao_max: rtdb_put failed', 911,
     &       RTDB_ERR)
         endif
         
      endif
c      write(6,*)' nq_chunk = ',nq_chunk
      if (rtdb_get(rtdb, 'dft:nq_chunk', mt_int, 1, nq_chunk))then
         if (me.eq.0)write(LuOut,*)' nq_chunk input override= ',nq_chunk
c
c        redefine max_pr_mbfnq
c
         max_pr_mbfnq = max_mbf*nq_chunk
      endif
c
c     check nquad_task (if .ne. 1 makes no sense with chunking turned on ... reset)
c
      if (nq_chunk.ne.0)then
         if (nquad_task.gt.1)then
            nquad_task = 1
c
         endif
      endif
      if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1, 
     &     nquad_task))
     &     call errquit('mbf_ao_max: rtdb_put failed', 911,
     &       RTDB_ERR)
c
c      write(6,*)' iAOacc, acc_AO_gauss: ',iAOacc, acc_AO_gauss
c      write(6,*)' max_mbf, max_pr_mbfnq, max_pr_mbf, max_pr_nq: ',
c     &            max_mbf, max_pr_mbfnq, max_pr_mbf, max_pr_nq
c
      return
      end
      double precision function dft_gaussian_range(n, alpha, eps)
      implicit none
c
      integer n
      double precision alpha, eps
c
c     Return an approximation to the outer solution of 
c     .     r^n*exp(-ar^2) = eps
c     .     r = (n*ln(-ln(eps)) - n*ln(a) - 4*ln(eps)) /
c     .         4*sqrt(-alpha*ln(eps))
c
c     Accuracy improves with smaller eps.
c
      double precision logeps
c
      logeps = log(eps)
c      
      dft_gaussian_range = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      end

      double precision function r_neglected(k, alpha, eps)
      implicit none
c
      integer k
      double precision alpha, eps
c
c     For a function f(r) = r^k*exp(-alpha*r^2) determine
c     the radial distance r such that the fraction of the 
c     function norm that is neglected if the 3D volume 
c     integration is terminated at a distance r is less
c     than or equal to eps.
c
      double precision r, test, neglected, step
c
      step = 0.5d0
      r = 1.0d0
 10   test = neglected(k,alpha,r)
      if (test .gt. eps) then
         r = r + step
      else 
         r = r - step
         if (r .lt. 0.0d0) r = 0.0d0
         step = step*0.5d0
         r = r + step
      endif
      if (step .gt. 0.01d0) goto 10
c
      r_neglected = r
c
      end
      double precision function neglected(k,alpha,r)
      implicit none
c
      integer k
      double precision alpha, r
c
c     For a function f(r) = r^k*exp(-alpha*r^2) determine
c     the fraction of the function norm that is neglected
c     if the 3D volume integration is terminated at a 
c     distance r.
c
c     neglected = int(t^2*f(t),t=r..infinity)/int(t^2*f(t),t=0..infinity)
c
      double precision ik
c
      neglected = ik(k+2,alpha,r)/ik(k+2,alpha,0.0d0)
c
      end
      double precision function ik(k,alpha,r)
      implicit none
c
      integer k
      double precision alpha, r
c
c     I(k) = int(t^k exp(-alpha*t^2), t=0..infinity)
c
c     I(k) = [(k-1)*I(k-2) + r^(k-1)*exp(-alpha*r^2)]/(2*alpha)
c
      integer i, ilo
      double precision value
#if defined(WIN32) ||defined(LINUX)
! DERFC is a GNU extension
      double precision derfc
#else
! ERFC is Fortran 2008
      double precision erfc
#endif
c
      ilo = mod(k,2)
c
      if (ilo .eq. 0) then
#if defined(WIN32) ||defined(LINUX)
         value = 0.5d0*sqrt(4.0d0*atan(1.0d0)/alpha)*
     $        derfc(sqrt(alpha)*r)
#else
         value = 0.5d0*sqrt(4.0d0*atan(1.0d0)/alpha)*
     $        erfc(sqrt(alpha)*r)
#endif
      else
         value = exp(-alpha*r*r)/(2.0d0*alpha)
      endif
c
      do i = ilo+2,k,2
         value = ((i-1)*value + r**(i-1)*exp(-alpha*r*r))/(2.0d0*alpha)
      enddo
c
      ik = value
c
      end

      Subroutine grid_repack(xyzw, qxyz, qwght, nq,rad,istep)
c
C$Id$
c
      implicit none
c
      integer nq
      double precision xyzw(4,nq), qxyz(3,nq), qwght(nq),
     .     rad
c
      integer n,istep
c
      istep=istep+1
      nq=dble(xyzw(1,istep))
      rad=xyzw(2,istep)
      do 30 n = 1, nq
c
        qxyz(1,n) = xyzw(1,n+istep)
        qxyz(2,n) = xyzw(2,n+istep)
        qxyz(3,n) = xyzw(3,n+istep)
c
        qwght(n) = xyzw(4,n+istep)
c
!            write(0,'(A,i4,4F16.9)') 
!     .           ' RE ',n+istep,
!     ,       qxyz(1,n),qxyz(2,n),qxyz(3,n),qwght(n)
   30 continue
      istep=istep+nq
c      write(6,*)' repacked buffer '
      return
      end
