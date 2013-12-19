
      logical function sym_mol(zname, nop, nir, iclass_dim,
     +     zir, zclass, chars, ipntmol)
c     
C$Id$
      implicit double precision (a-h,p-w),integer (i-n),logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "errquit.fh"
      parameter(maxireps=20)
#include "util.fh"
c     
      dimension iclass_dim(maxireps), ipointer(maxireps), zir(maxireps),
     +   zclass(maxireps), chars(maxireps*maxireps), ipntmol(maxireps)
c
c     Ipointer() used to return mapping of classes from gamess order
c     to standard order ... seems of no use now but I have left the code
c     lying around ... note that some groups lack this info
c
      integer ic1
      double precision chc1
      character*10 yc1
      character*4 yrc1
      parameter(nc1=1)
      dimension ic1(5), chc1(1,1), yc1(1), yrc1(1)
c     
c     
      integer ic2
      double precision chc2
      character*10 yc2
      character*4 yrc2
      parameter(nc2=2)
      dimension ic2(8), chc2(2,2), yc2(2), yrc2(2)
c     
c     
      integer ic2h
      double precision chc2h
      character*10 yc2h
      character*4 yrc2h
      parameter(nc2h=4)
      dimension ic2h(14), chc2h(4,4), yc2h(4), yrc2h(4)
c     
c     
      integer ic2v
      double precision chc2v
      character*10 yc2v
      character*4 yrc2v
      parameter(nc2v=4)
      dimension ic2v(14), chc2v(4,4), yc2v(4), yrc2v(4)
c     
c     
      integer ici
      double precision chci
      character*10 yci
      character*4 yrci
      parameter(nci=2)
      dimension ici(8), chci(2,2), yci(2), yrci(2)
c     
c     
      integer ics
      double precision chcs
      character*10 ycs
      character*4 yrcs
      parameter(ncs=2)
      dimension ics(8), chcs(2,2), ycs(2), yrcs(2)
c     
c     
      integer id2
      double precision chd2
      character*10 yd2
      character*4 yrd2
      parameter(nd2=4)
      dimension id2(14), chd2(4,4), yd2(4), yrd2(4)
c     
c     
      integer id2h
      double precision chd2h
      character*10 yd2h
      character*4 yrd2h
      parameter(nd2h=8)
      dimension id2h(26), chd2h(8,8), yd2h(8), yrd2h(8)
c     
      logical please_print_ecce
      data please_print_ecce / .true. /
c
      data ic1/1,1,1
     +     ,1
     +     ,1/
      data yrc1/'a'/
      data yc1/'e'/
      data chc1/
     +     1.000/
c     
      data ic2/2,2,1,1
     +     ,1,2
     +     ,1,2/
      data yrc2/'a','b'/
      data yc2/'e','c2'/
      data chc2/
     +     1.000, 1.000,
     +     1.000,-1.000/
c     
      data ic2h/4,4,1,1,1,1
     +     ,1,2,3,4
     +     ,1,3,2,4/
      data yrc2h/'ag','au','bg','bu'/
      data yc2h/'e','c2','si.h','i'/
      data chc2h/
     +     1.000, 1.000, 1.000, 1.000,
     +     1.000, 1.000,-1.000,-1.000,
     +     1.000,-1.000,-1.000, 1.000,
     +     1.000,-1.000, 1.000,-1.000/
c     
      data ic2v/4,4,1,1,1,1
     +     ,1,2,4,3
     +     ,1,3,2,4/
      data yrc2v/'a1','a2','b1','b2'/
      data yc2v/'e','c2','si.xz','si.yz'/
      data chc2v/
     +     1.000, 1.000, 1.000, 1.000,
     +     1.000, 1.000,-1.000,-1.000,
     +     1.000,-1.000, 1.000,-1.000,
     +     1.000,-1.000,-1.000, 1.000/
c     
      data ici/2,2,1,1
     +     ,1,2
     +     ,1,2/
      data yrci/'ag','au'/
      data yci/'e','i=s2'/
      data chci/
     +     1.000, 1.000,
     +     1.000,-1.000/
c     
      data ics/2,2,1,1
     +     ,1,2
     +     ,1,2/
      data yrcs/'a''','a"'/
      data ycs/'e','si.h'/
      data chcs/
     +     1.000, 1.000,
     +     1.000,-1.000/
c     
      data id2/4,4,1,1,1,1
     +     ,1,2,3,4
     +     ,1,2,3,4/
      data yrd2/'a1','b1','b2','b3'/
      data yd2/'e','c2z','c2y','c2x'/
      data chd2/
     +     1.000, 1.000, 1.000, 1.000,
     +     1.000, 1.000,-1.000,-1.000,
     +     1.000,-1.000, 1.000,-1.000,
     +     1.000,-1.000,-1.000, 1.000/
c     
      data id2h/8,8,1,1,1,1,1,1,1,1
     +     ,1,2,5,6,4,3,8,7
     +     ,1,8,4,5,6,3,7,2/
      data yrd2h/'ag','au','b1g','b1u','b2g','b2u','b3g','b3u'/
      data yd2h/'e','c2z','c2y','c2x','i','si.xy','si.xz','si.yz'/
      data chd2h/
     +     1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
     +     1.000, 1.000, 1.000, 1.000,-1.000,-1.000,-1.000,-1.000,
     +     1.000, 1.000,-1.000,-1.000, 1.000, 1.000,-1.000,-1.000,
     +     1.000, 1.000,-1.000,-1.000,-1.000,-1.000, 1.000, 1.000,
     +     1.000,-1.000, 1.000,-1.000, 1.000,-1.000, 1.000,-1.000,
     +     1.000,-1.000, 1.000,-1.000,-1.000, 1.000,-1.000, 1.000,
     +     1.000,-1.000,-1.000, 1.000, 1.000,-1.000,-1.000, 1.000,
     +     1.000,-1.000,-1.000, 1.000,-1.000, 1.000, 1.000,-1.000/
c     
      zgroup = zname
      call inp_lcase(zgroup)
c     
c     here is beginning of instructions - also machine generated
c     
      if(zgroup.eq.'c1') then
         k=0
         nir=nc1
         nop=ic1(2)
         do i=1,nir
            iclass_dim(i)=ic1(i+2)
            ipointer(i)=ic1(i+2+nir)
            ipntmol(i)=ic1(i+3+nir)
            zir(i)=yrc1(i)
            zclass(i)=yc1(i)
            do j=1,nir
               k=k+1
               chars(k)=chc1(j,i)
            enddo
         enddo
         goto 10
      endif
c     
c     
      if(zgroup.eq.'c2') then
         k=0
         nir=nc2
         nop=ic2(2)
         do i=1,nir
            iclass_dim(i)=ic2(i+2)
            ipointer(i)=ic2(i+2+nir)
            ipntmol(i)=ic2(i+4+nir)
            zir(i)=yrc2(i)
            zclass(i)=yc2(i)
            do j=1,nir
               k=k+1
               chars(k)=chc2(j,i)
            enddo
         enddo
         goto 10
      endif
c     
c     
      if(zgroup.eq.'c2h') then
         k=0
         nir=nc2h
         nop=ic2h(2)
         do i=1,nir
            iclass_dim(i)=ic2h(i+2)
            ipointer(i)=ic2h(i+2+nir)
            ipntmol(i)=ic2h(i+6+nir)
            zir(i)=yrc2h(i)
            zclass(i)=yc2h(i)
            do j=1,nir
               k=k+1
               chars(k)=chc2h(j,i)
            enddo
         enddo
         goto 10
      endif
c     
c     
      if(zgroup.eq.'c2v') then
         k=0
         nir=nc2v
         nop=ic2v(2)
         do i=1,nir
            iclass_dim(i)=ic2v(i+2)
            ipointer(i)=ic2v(i+2+nir)
            ipntmol(i)=ic2v(i+6+nir)
            zir(i)=yrc2v(i)
            zclass(i)=yc2v(i)
            do j=1,nir
               k=k+1
               chars(k)=chc2v(j,i)
            enddo
         enddo
         goto 10
      endif
c     
c     
      if(zgroup.eq.'ci') then
         k=0
         nir=nci
         nop=ici(2)
         do i=1,nir
            iclass_dim(i)=ici(i+2)
            ipointer(i)=ici(i+2+nir)
            ipntmol(i)=ici(i+4+nir)
            zir(i)=yrci(i)
            zclass(i)=yci(i)
            do j=1,nir
               k=k+1
               chars(k)=chci(j,i)
            enddo
         enddo
         goto 10
      endif
c     
c     
      if(zgroup.eq.'cs') then
         k=0
         nir=ncs
         nop=ics(2)
         do i=1,nir
            iclass_dim(i)=ics(i+2)
            ipointer(i)=ics(i+2+nir)
            ipntmol(i)=ics(i+4+nir)
            zir(i)=yrcs(i)
            zclass(i)=ycs(i)
            do j=1,nir
               k=k+1
               chars(k)=chcs(j,i)
            enddo
         enddo
         goto 10
      endif
c     
c     
      if(zgroup.eq.'d2h') then
         k=0
         nir=nd2h
         nop=id2h(2)
         do i=1,nir
            iclass_dim(i)=id2h(i+2)
            ipointer(i)=id2h(i+2+nir)
            ipntmol(i)=id2h(i+10+nir)
            zir(i)=yrd2h(i)
            zclass(i)=yd2h(i)
            do j=1,nir
               k=k+1
               chars(k)=chd2h(j,i)
            enddo
         enddo
c        if(ga_nodeid().eq.0)write(6,'(8I3)')(ipntmol(i),i=1,8)
         goto 10
      endif
c     
c     here it will come if group was not recognized
c     
      sym_mol = .false.
      return
c
c     here it comes on success
c
 10   continue
      sym_mol = .true.
c
c     Munge fractional characters to their precise values
c
      twopi = 2.0d0*4.0d0*atan(1.0d0)
      do i = 1, nir*nir
c     
c     Integers, rational fractions and square-roots thereof
c
         do k = 0, 10
            do l = 1, 10
               dval = dble(k)/dble(l)
               dval = sign(dval,chars(i))
               if (abs(chars(i)-dval).lt.0.0006) goto 20
               dval = sqrt(abs(dval))
               dval = sign(dval,chars(i))
               if (abs(chars(i)-dval).lt.0.0006) goto 20
            enddo
         enddo
c
c     Sines and cosines and multiples thereof
c
         do n = 1, 24
            do k = 1, (n-1)/4+1
               angle = dble(k)*twopi/dble(n)
               do multiple = 1, 2
                  dval = multiple*sin(angle)
                  dval = sign(abs(dval),chars(i))
                  if (abs(chars(i)-dval).lt.0.0006) goto 20
                  dval = multiple*cos(angle)
                  dval = sign(abs(dval),chars(i))
                  if (abs(chars(i)-dval).lt.0.0006) goto 20
               enddo
            enddo
         enddo
c
c     Magic numbers ... shouldn't be any?
c
         write(6,*) chars(i)
         call errquit('sym_char_table: unmunged character',i, INPUT_ERR)
c
 20      chars(i) = dval
      enddo
c
      if (util_print('sym_char_table', print_never)) then
         write(6,*) ' sym_char_table ', zname
         do i = 0, nir-1
            write(6,111) (chars(i*nir+j),j=1,nir)
         enddo
 111     format(1x,20f10.6)
         call util_flush(6)
      endif
c
      if (please_print_ecce) then
         please_print_ecce = .false.
         call ecce_print1_char('group class names', zclass, nir)
         call ecce_print1_char('group irrep names', zir, nir)
      endif
c
      end
