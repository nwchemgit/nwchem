C $Id$
      Subroutine drdy_sort(lhead,jop,lsad,ns,iflg)
C
C Sort RPH interpolation information to give points along MEP in 
C  increasing value of s (iflg#0, final sort and print out summary)
C
      Implicit None
C
#include "drdyP.fh"
#include "stdio.fh"
C
*::passed::
      Integer jop,iflg
      Logical lhead,lsad
*::local:
      Integer i,ii,ierr,iop,is,j,jopmx,js,ks,llopt(40),nbar,nend,
     * nfreq,ns,indx
      Character*80 line
      Character*255 filename
      Double Precision t
      Logical isopen
C
      nend  = 0
      jopmx = 0
      nfreq = 0
      nbar  = 0
C
c Get the filename for fu30, make sure it is the same as the one in rphwrt.F
c
      call util_file_name('file30',.false.,.false.,filename)
      ierr = 0
      write (fu6,600)
      isopen = .false.
      inquire(unit=fu30,opened=isopen)
      if (isopen) then
        close (unit=fu30, status='keep')
      endif
      open (unit=fu30,status='old',err=1000,form='formatted',
     &      file=filename)
C Position fu30 file at beginning of MEP data
      call drdy_restor(lhead,jop,lsad,llopt)
      jopmx = 3
      if (lgs(6).eq.1.or.lgs(6).eq.3) jopmx = 4
      ns = 0
      if (.not.lhead .or. .not.lsad .or. jop.lt.jopmx) then
         write (fu6,*) ' trouble reading fu30 file, can not sort'
         return
      endif
      is = 0
      if (lgs(1).eq.0) lsad = .false.
      read (fu30,3000,end=100) line
C      write (fu6,3000) line
      read (fu30,3000) line
C      write (fu6,3000) line
      if (lgs(1).eq.0) lsad = .true.
      nfreq = nf(5)
      nbar = (n3*(n3+1))/2
C Read and sort MEP data
      do while (is.le.nsmax+1)
         is = is + 1
         read (fu30,3000,end=100) line
C         write(fu6,*) line
         if (lgs(39).eq.0) then
            read (fu30,*,end=100,err=100) s,v
         else
            read (fu30,*,end=100,err=100) s,v,vspc
         endif
C         write(fu6,*) s,v
         read (fu30,3000,end=100) line
C         write(fu6,*) line
         read (fu30,*,end=100,err=100) (x(i),i=1,n3)
C         write(fu6,*) (x(i),i=1,n3)
C Find insertion point, js, for newly read s value
         js = 1
         do while (js.lt.is .and. s.gt.sgrid(js))
            js = js + 1
         enddo
         ks = is
C Insertion point, js, is less than the total number of points, is
C   Shift grid from js and above up one
         do while (ks.gt.js)
            sgrid(ks) = sgrid(ks-1)
            vgrid(ks) = vgrid(ks-1)
            vgrdsp(ks) = vgrdsp(ks-1)
            do i = 1,n3
               xgrid(i,ks) = xgrid(i,ks-1)
               dxgrid(i,ks) = dxgrid(i,ks-1)
               ifrgrd(i,ks) = ifrgrd(i,ks-1)
               freqg(i,ks) = freqg(i,ks-1)
            enddo
            do i = 1,nbar
               hgrid(i,ks) = hgrid(i,ks-1)
            enddo
            ks = ks - 1
         enddo
C Read new point in at insertion point
         sgrid(js) = s
         vgrid(js) = v
         vgrdsp(js) = vspc
         do i = 1,n3
            xgrid(i,js) = x(i)
         enddo
         read (fu30,3000,end=100) line
C         write(fu6,*) line
         read (fu30,*,end=100,err=100) (dxgrid(i,js),i=1,n3)
C         write(fu6,*) (dxgrid(i,js),i=1,n3)
         read (fu30,*) (ifrgrd(i,js),i=1,nfreq)
C         write(fu6,*) (ifrgrd(i,js),i=1,nfreq)
         read (fu30,3000) line
C         write(fu6,*) line
         read (fu30,*) (hgrid(i,js),i=1,nbar)
C         write(fu6,*) (hgrid(i,js),i=1,nbar)
C Project and diagonalize hessian
         indx = 0
         do i = 1,n3
            dx(i) = dxgrid(i,js)
            do j = 1,i
               f(i,j) = hgrid(indx+j,js)
            enddo
            indx = indx + i
         enddo
         do i = 2,n3
            do j = 1,i-1
               f(j,i) = f(i,j)
            enddo
         enddo
         call drdy_projct(ierr)
         if(ierr.eq.0) call drdy_fdiag(n3,nf(5),ierr,0)
         if (ierr.eq.0) then
C Store frequencies in freqg array
            do i = 1,nf(5)
               freqg(i,js) = scr2(i)
            enddo
         endif
C
      enddo
      if (is.gt.nsmax) then
         write (fu6,*) ' maximum number of grid points is', nsmax
         write (fu6,*) ' unable to complete sort of data'
         close (fu30)
         call errquit('drdy_sort: fatal error:1',911,0)
      endif
100   continue
      close (fu30)
      ns = is - 1
C
      Write(fu6,601)
      if (ns.le.0) then
         if (jop.ge.1) write (fu6,602) sr,vr
         if (jop.ge.3) write (fu6,603) sp,vp
         if (lgs(1).ne.0.and.lsad) write (fu6,604) 0.0,vspsv,
     *    (xr(i,5),i=1,n3)
      else
C Calculate s at reactants and products
         sr = 0.0d00
         sp = 0.0d00
         do i = 1,n3
            t = xgrid(i,1)-xr(i,1)
            sr = sr + t*t
            t = xgrid(i,ns)-xr(i,3)
            sp = sp + t*t
         enddo
         sr = sgrid(1) - sqrt(sr)
         sp = sgrid(ns) + sqrt(sp)
         write (fu6,602) sr,vr,(xr(i,1),i=1,n3)
         write (fu6,603) sp,vp,(xr(i,3),i=1,n3)
C
C Rewrite MEP data
         open (unit=fu30,form='formatted',file=filename)
         write (fu30,3000) ctitle
         write (fu30,3005) llopt
         iop = 1
         do while (iop.le.4)
C Reactants and products
            write (fu30,3010) iop
            if (iop.eq.1) then
               if (lgs(39).eq.0) then
                  write (fu30,3020) sr,vr,vzero
               else
                  write (fu30,3021) sr,vr,vzero,vrsp,vzerosp
               endif
            else if (iop.eq.2) then
               if (lgs(39).eq.0) then
                  write (fu30,3020) sr,0.0
               else
                  write (fu30,3020) sr,0.0,0.0
               endif
            else if (iop.eq.3) then
               if (lgs(39).eq.0) then
                  write (fu30,3020) sp,vp
               else
                  write (fu30,3020) sp,vp,vpsp
               endif
            else
               if (lgs(39).eq.0) then
                  write (fu30,3020) sp,0.0
               else
                  write (fu30,3020) sp,0.0,0.0
               endif
            endif
            if (icode(iop).eq.1) then
               write (fu30,*) ' '
            elseif (icode(iop).eq.2) then
               write (fu6,6000) iop
               call errquit('drdy_sort: fatal error:2',911,0)
            else
               nfreq = nf(iop)
               nend = 3*nratom(iop)
               write (fu30,3005) (ifrgrp(i,iop),i=1,nfreq)
               write (fu30,3040)
               nbar = (nend*(nend+1))/2
               ii = 0
               do i = 1,nend
                  write (fu30,3020) (hgrp(ii+j,iop),j=1,i)
                  ii = ii + i
               enddo
            endif
            if ((iop.eq.1 .and. lgs(6).le.2) .or. 
     *       iop.eq.2. or.
     *       (iop.eq.3 .and. (lgs(6).eq.1.or.lgs(6).eq.3)) .or.
     *       iop.eq.4) then
               iop = iop + 1
            else
               iop = iop + 2
            endif
         enddo
         nfreq = nf(5)
         if(lgs(1).ne.0) then
            write (fu30,3050)
            s = 0.0
            write (fu6,604) s,vspsv,(xr(i,5),i=1,n3)
            if(lgs(39).eq.0) then
               write (fu30,3020) s,vspsv
            else
               write (fu30,3020) s,vspsv,vspspsv
            endif
            write (fu30,3005) (ifrsp(i),i=1,nfreq)
C  write out packed hessian in mass-weighted coordinates
            write (fu30,3040)
            do i = 1, n3
               write (fu30,3020) (fspsv(i,j),j=1,i)
            enddo
         endif
         write (fu30, 3060) ns
         nbar = (n3*(n3+1))/2
         write (fu6,*) ' '
         do is = 1,ns
            write (fu6,605) sgrid(is),vgrid(is),(xgrid(i,is),i=1,n3)
            if (lgs(39).eq.0) then
               write (fu30,3070) sgrid(is),vgrid(is)
            else
               write (fu30,3070) sgrid(is),vgrid(is),vgrdsp(is)
            endif
c  write out mass-weighted x
            write (fu30,3030) (xgrid(i,is),i=1,n3)
c  write out mass-weighted dx
            write (fu30,3080) (dxgrid(i,is),i=1,n3)
            write (fu30,3005) (ifrgrd(i,is),i=1,nfreq)
c  write out packed hessian in mass-weighted coordinates
            write (fu30,3040)
            ii = 0
            do i = 1,n3
               write (fu30,3020) (hgrid(ii+j,is),j=1,i)
               ii = ii + i
            enddo
         enddo
         if (iflg.ne.0) then
            write (fu30,3085)
            write (fu30,3090) 0,sr,vr
            write (fu30,3095)
            write (fu30,3090) 0,sp,vp
            call drdy_sumout(ns)
         endif
      endif
      close (fu30)
      return
1000  continue
      write (fu6,*) ' error opening fu30 file'
      call errquit('drdy_sort: fatal error:3',911,0)
600   Format (//,' Sort data in file 30 so points along MEP',
     *  ' have the reaction coordinate in ascending order.',/,
     *  ' Also calculate the values of the reaction coordinate',
     *  ' at the reactant and product geometries',/,
     *  '  and shift energies so reactant v is zero.')
601   Format (/,' Summary of sorted data (all units are atomic)',/,5x,
     * 'R = reactants, P = products, * = saddle point or starting',
     * ' geometry',//,t10,'s',11x,'V',T38,'x (mass weighted)')
602   Format (' R',t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
603   Format (' P',t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
604   Format (' *',t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
605   Format (t4,0pf12.5,1pe13.5,(T35,0p9f9.4))
3000  Format (a80)
3005  Format (1x,i2,39i3)
3010  Format (' iop =', i5, ', s,v=')
3020  Format (1x, 1pe19.10, 3e20.10)
3021  Format (1x, 1pe19.10, 4e20.10)
3030  Format (' x='/ (1x, 1pe19.10, 3e20.10))
3040  Format (' f matrix (packed)=')
3050  Format (' saddle point,  s,v=')
3060  Format (' begin points along reaction coordinate, ns=',/,1x,i5)
3070  Format (' s,v='/ (1x, 1pe19.10, 3e20.10))
3080  Format (' dx='/ (1x, 1pe19.10, 3e20.10))
3085  Format (' option for fitting potential in reactant valley,',
     *  ' ifitvr,s,v=')
3090  Format (1x,i5,1p2e13.5)
3095  Format (' option for fitting potential in product valley,',
     *  ' ifitvp,s,v=')
6000  Format (' icode=2 was chosen for iop=',i5,', but it is not',
     *   ' supported')
      END
