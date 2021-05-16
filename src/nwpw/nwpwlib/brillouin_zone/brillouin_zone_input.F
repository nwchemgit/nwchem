c
c $Id$
c

*     *************************************************
*     *                                               *
*     *             brillouin_zone_input              *
*     *                                               *
*     *************************************************
      subroutine brillouin_zone_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 7)
      character*22 dirs(num_dirs)
      character*255 test

      character*50 rtdb_name
      character*50 zone_name
      integer          num_kvectors,k,l,num_kvectors_print
      double precision kvector(4)
      integer kvs(2),kvs_new(2)

      data dirs / 'zone_name:','zone_name',
     >            'kvector:','kvector',
     >            'path',
     >            'max_kpoints_print',
     >            'end'/

c     
c     **** intialize stuff ****
      num_kvectors       = 0
      num_kvectors_print = 0
      zone_name          = 'zone_default'
c
c
 10   if (.not. inp_read()) 
     >     call errquit(
     >           'brillouine_zone_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'brillouine_zone_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'brillouine_zone_input: unknown directive', 0,
     &       INPUT_ERR)


      goto ( 100,100,200,200,300,400,
     >      9999) ind
      call errquit(
     >      'brillouine_zone_input: unimplemented directive', ind,
     &       INPUT_ERR)

c     
c  zone_name
c     
 100  if (.not. inp_a(zone_name))
     >     call errquit(
     >          'brillouine_zone_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c kvector
c
 200  if (.not. inp_f(kvector(1)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0,
     &       INPUT_ERR)
      if (.not. inp_f(kvector(2)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0,
     &       INPUT_ERR)
      if (.not. inp_f(kvector(3)))
     >     call errquit(
     >          'brillouine_zone_input: failed to kvector', 0,
     &       INPUT_ERR)

      if (.not. inp_f(kvector(4))) kvector(4) = -1.0d0
     

      num_kvectors = num_kvectors + 1
      if (.not. BA_alloc_get(mt_dbl,(4*num_kvectors), 
     >        'kvs_new',kvs_new(2),kvs_new(1)))
     >     call errquit(
     >          'brillouine_zone_input: heap failed 1', 0,
     &       INPUT_ERR)

      if (num_kvectors.gt.1) then
        do k=1,(num_kvectors-1)
         dbl_mb(kvs_new(1) + 4*(k-1))   = dbl_mb(kvs(1) + 4*(k-1))
         dbl_mb(kvs_new(1) + 4*(k-1)+1) = dbl_mb(kvs(1) + 4*(k-1)+1)
         dbl_mb(kvs_new(1) + 4*(k-1)+2) = dbl_mb(kvs(1) + 4*(k-1)+2) 
         dbl_mb(kvs_new(1) + 4*(k-1)+3) = dbl_mb(kvs(1) + 4*(k-1)+3) 
        end do
        if (.not.BA_free_heap(kvs(2)))
     >     call errquit(
     >          'brillouine_zone_input: heap failed 2', 0, MA_ERR)
      end if
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1))   = kvector(1) 
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1)+1) = kvector(2) 
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1)+2) = kvector(3) 
      dbl_mb(kvs_new(1) + 4*(num_kvectors-1)+3) = kvector(4) 
      kvs(1) = kvs_new(1)
      kvs(2) = kvs_new(2)

      goto 10

c
c path
c
 300  call band_path_set(rtdb,num_kvectors,kvs)

      goto 10


c
c max_kpoints_print
c
 400  if (.not.inp_i(num_kvectors_print)) num_kvectors_print = 0

      goto 10


*     ***** add brillouin_zone to rtdb ****
 9999 continue

      l = index(zone_name,' ') -1

      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >  call errquit(
     >  'brillouin_zone_input: rtdb_put(number_kvectors) failed', 0,
     &       RTDB_ERR)

      if (num_kvectors_print.gt.0) then
         rtdb_name = zone_name(1:l)//':number_kvectors_print'
         if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_kvectors_print))
     >     call errquit(
     >     'brillouin_zone_input: rtdb_put(number_kvectors) failed', 0,
     >       RTDB_ERR)
      end if

      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1))))
     >  call errquit(
     >  'brillouin_zone_input: rtdb_put(number_kvectors) failed', 0,
     &       RTDB_ERR)


      if (.not.BA_free_heap(kvs(2)))
     >  call errquit(
     >       'brillouine_zone_input: heap failed 3', 0, MA_ERR)
      return
      end


*     *************************************************
*     *                                               *
*     *             monkhost_pack_set                  *
*     *                                               *
*     *************************************************

      subroutine monkhorst_pack_set(rtdb,zone_name,nx,ny,nz)
      implicit none
      integer rtdb
      character*(*) zone_name
      integer nx,ny,nz

#include "errquit.fh"
#include "bafdecls.fh"
#include "rtdb.fh"

      !**** local variables ****
      logical timereverse
      character*50 rtdb_name
      integer i1,i2,i3,k,l,num_kvectors
      integer kvs(2),nn(3)
      real*8  kx,ky,kz,xx,yy,zz,xxx,yyy,zzz,weight

      timereverse = .true.
      if (nx.lt.0) then
         nx = -nx
         timereverse = .false.
      end if
      if (ny.lt.0) then
         ny = -ny
         timereverse = .false.
      end if
      if (nz.lt.0) then
         nz = -nz
         timereverse = .false.
      end if

      num_kvectors = nx*ny*nz
      weight       = 1.0d0/dble(num_kvectors)

      if(.not.BA_alloc_get(mt_dbl,(4*num_kvectors),'kvs',kvs(2),kvs(1)))
     >     call errquit(
     >          'monkhorst_pack_set: heap failed 1', 0,
     >       MA_ERR)

      xxx = 1.0d0/(2.0d0*nx)
      yyy = 1.0d0/(2.0d0*ny)
      zzz = 1.0d0/(2.0d0*nz)
      k = 0
      do i3=0,(nz-1)
      do i2=0,(ny-1)
      do i1=0,(nx-1)
         xx = 1.0d0 + 2*i1 - nx
         yy = 1.0d0 + 2*i2 - ny
         zz = 1.0d0 + 2*i3 - nz
 
         kx = xx*xxx
         ky = yy*yyy
         kz = zz*zzz
 
         dbl_mb(kvs(1) + 4*(k))   = kx
         dbl_mb(kvs(1) + 4*(k)+1) = ky
         dbl_mb(kvs(1) + 4*(k)+2) = kz
         dbl_mb(kvs(1) + 4*(k)+3) = weight

         k = k + 1
      end do
      end do
      end do

      if (timereverse) 
     > call monkhorst_pack_timereversal_prune(num_kvectors,
     >                                        dbl_mb(kvs(1)))


*     ***** add brillouin_zone to rtdb ****
      l = index(zone_name,' ') -1

      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >  call errquit(
     >  'monkhorst_pack_set: rtdb_put(number_kvectors) failed', 0,
     &       RTDB_ERR)

      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1))))
     >  call errquit(
     >  'monkhorst_pack_set: rtdb_put(kvectors) failed', 0,
     &       RTDB_ERR)

      if (timereverse) then
         nn(1) = nx
         nn(2) = ny
         nn(3) = nz
      else
         nn(1) = -nx
         nn(2) = -ny
         nn(3) = -nz
      end if
      rtdb_name = zone_name(1:l)//':monkhorst-pack'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,3,nn))
     >  call errquit(
     >  'monkhorst_pack_set: rtdb_put(monkhorst-pack) failed', 0,
     &       RTDB_ERR)


      if (.not.BA_free_heap(kvs(2)))
     >  call errquit(
     >       'monkhorst_pack_set: heap failed 3', 0, MA_ERR)

      return 
      end

*     *************************************************
*     *                                               *
*     *        monkhorst_pack_timereversal_prune      *
*     *                                               *
*     *************************************************
      subroutine monkhorst_pack_timereversal_prune(nks,ks)
      implicit none
      integer nks
      real*8 ks(4,*)
      integer i,j,nks2
      real*8 norm
      
      do i=1,nks-1
         if (ks(4,i).gt.0.0d0) then
         do j=i+1,nks
           if ((dabs(ks(1,i)+ks(1,j)).lt.1.0d-9).and.
     >         (dabs(ks(2,i)+ks(2,j)).lt.1.0d-9).and.
     >         (dabs(ks(3,i)+ks(3,j)).lt.1.0d-9)) then
              ks(4,i) = ks(4,i) + ks(4,j)
              ks(4,j) = 0.0d0
           end if
         end do
         end if
      end do
      nks2 = 0
      do i=1,nks
         if (ks(4,i).gt.0.0d0) nks2 = nks2+1
      end do
      i = 1
      do while (i.lt.nks2)
         if (ks(4,i).eq.0.0d0) then
            do j=i+1,nks
               ks(1,j-1) = ks(1,j)
               ks(2,j-1) = ks(2,j)
               ks(3,j-1) = ks(3,j)
               ks(4,j-1) = ks(4,j)
            end do
         else
           i = i+1
         end if
      end do
      do i=1,nks2
        norm = ks(1,i) + ks(2,i) + ks(3,i)
        if (norm.lt.0.0d0) then
           ks(1,i) = -ks(1,i)
           ks(2,i) = -ks(2,i)
           ks(3,i) = -ks(3,i)
        end if
      end do
      nks = nks2
      return
      end


*     *************************************************
*     *                                               *
*     *             dos_states_set                    *
*     *                                               *
*     *************************************************

      subroutine dos_states_set(rtdb,zone_name,nx,ny,nz)
      implicit none
      integer rtdb
      character*(*) zone_name
      integer nx,ny,nz

#include "errquit.fh"
#include "bafdecls.fh"
#include "rtdb.fh"

      !**** local variables ****
      character*50 rtdb_name
      integer i1,i2,i3,k,l,num_kvectors,nn(3)
      integer kvs(2)
      real*8  kx,ky,kz,xx,yy,zz,xxx,yyy,zzz


      num_kvectors = nx*ny*nz

      if(.not.BA_alloc_get(mt_dbl,(4*num_kvectors),'kvs',kvs(2),kvs(1)))
     >     call errquit(
     >          'dos_states_set: heap failed 1', 0,
     >       MA_ERR)

      xxx = 1.0d0/(1.0d0*nx)
      yyy = 1.0d0/(1.0d0*ny)
      zzz = 1.0d0/(1.0d0*nz)
      k = 0
      do i3=0,(nz-1)
      do i2=0,(ny-1)
      do i1=0,(nx-1)
         xx = i1
         yy = i2
         zz = i3
 
         kx = xx*xxx
         ky = yy*yyy
         kz = zz*zzz
 
         dbl_mb(kvs(1) + 4*(k))   = kx
         dbl_mb(kvs(1) + 4*(k)+1) = ky
         dbl_mb(kvs(1) + 4*(k)+2) = kz
         dbl_mb(kvs(1) + 4*(k)+3) = -1.0d0

         k = k + 1
      end do
      end do
      end do


*     ***** add brillouin_zone to rtdb ****
      l = index(zone_name,' ') -1

      rtdb_name = zone_name(1:l)//':number_kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,num_kvectors))
     >  call errquit(
     >  'dos_states_set: rtdb_put(number_kvectors) failed', 0,
     >       RTDB_ERR)

      rtdb_name = zone_name(1:l)//':kvectors'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,
     >                   (4*num_kvectors),
     >                    dbl_mb(kvs(1))))
     >  call errquit(
     >  'dos_states_set: rtdb_put(number_kvectors) failed', 0,
     >       RTDB_ERR)

      nn(1) = nx
      nn(2) = ny
      nn(3) = nz
      rtdb_name = zone_name(1:l)//':dos-grid'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,3,nn))
     >  call errquit(
     >  'dos_states_set: rtdb_put(monkhorst-pack) failed', 0,
     &       RTDB_ERR)


      if (.not.BA_free_heap(kvs(2)))
     >  call errquit(
     >       'dos_states_set: heap failed 3', 0, MA_ERR)

      return 
      end

*     *************************************************
*     *                                               *
*     *             band_path_set                     *
*     *                                               *
*     *************************************************

      subroutine band_path_set(rtdb,nk,kvs)
      implicit none
      integer rtdb
      integer nk
      integer kvs(2)

#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"

*     *** local variables ***
      integer i,k,ii,ind,count,count2
      character*50 test

      integer num_cells         ! No. of known directives
      parameter (num_cells = 8)
      character*22 cells(num_cells)
      data cells / 'sc','fcc','bcc',
     >             'rhombohedral','hexagonal',
     >             'simple-tetragonal',
     >             'simple-orthorhombic',
     >             'bct' /

c CUB,FCC,BCC,TET,BCT1,BCT2,ORC,ORCF1,ORCF2,ORCF3,ORCI,ORCC,HEX,RHL1,RHL2,MCL,
c MCLC1,MCLC2,MCLC3,MCLC4,MCLC5,TRI1a,TRI2a,TRI1b,TRI2b
c FCCconventional,BCCconventional,BCT1conventional,BCT2conventional
c ORCF1conventional,ORCF2conventional,ORCF3conventional,ORCIconventional
c ORCCconventional
c MCLC1conventional,MCLC2conventional,MCLC3conventional,MCLC4conventional,MCLC5conventional
c
c Standard unit cell types
c Simple cubic
c Face-centered cubic - primitive,conventional
c Body-centered cubic - primitive,conventional
c Tetragonal
c Body-centered tegragonal - primitive,conventional
c Orthorhombic
c Face-centered orthorhombic - primitive,conventional
c Body-centered orthorhombic - primitive,conventional
c C-centered orthorhombic - primitive,conventional
c Hexagonal
c Rhombohedral
c Monoclinic
c C-centered monoclinic - primitive,conventional
c Triclinic

      integer pcount,  plist(100)
      integer pcount2, plist2(100)
      real*8  kvector(4),kvector1(4),kvector2(4),dist,dk
      logical goodpoint
      integer kvs_new(2)


      if (.not. inp_f(dk)) dk = 0.05d0  !** read in dk increment

      if (.not. inp_a(test))
     >     call errquit(
     >           'band_path_set: inp_read failed', 0, INPUT_ERR)

      if (.not. inp_match(num_cells, .false., test, cells, ind))
     >     call errquit(
     >           'band_path_set: unknown cell', 0,
     >       INPUT_ERR)
     
      call band_path_plist(pcount,plist)


*     *** remove non-special points ***
      pcount2 = 0
      do ii=1,pcount
        call band_path_point(ind,plist(ii),kvector,goodpoint)
        if (goodpoint) then
          pcount2=pcount2+1
          plist2(pcount2) = plist(ii)
        end if
      end do

*     *** add segments ***
      do ii=1,pcount2-1
        call band_path_point(ind,plist2(ii),  kvector1,goodpoint)
        call band_path_point(ind,plist2(ii+1),kvector2,goodpoint)
        dist = dsqrt( (kvector2(1)-kvector1(1))**2
     >              + (kvector2(2)-kvector1(2))**2
     >              + (kvector2(3)-kvector1(3))**2)
        count  = dist/dk
        count2 = count
        if (ii.eq.(pcount2-1)) count2=count2+1
        do i=1,count2
          kvector(1) = kvector1(1) 
     >               + (i-1)*((kvector2(1)-kvector1(1))/dble(count))
          kvector(2) = kvector1(2) 
     >               + (i-1)*((kvector2(2)-kvector1(2))/dble(count))
          kvector(3) = kvector1(3) 
     >               + (i-1)*((kvector2(3)-kvector1(3))/dble(count))

          nk = nk + 1
          if (.not. BA_alloc_get(mt_dbl,(4*nk),
     >        'kvs_new',kvs_new(2),kvs_new(1)))
     >        call errquit(
     >       'band_path_set: heap failed 1', 0, MA_ERR)

           if (nk.gt.1) then
            do k=1,(nk-1)
            dbl_mb(kvs_new(1) + 4*(k-1))   = dbl_mb(kvs(1) + 4*(k-1))
            dbl_mb(kvs_new(1) + 4*(k-1)+1) = dbl_mb(kvs(1) + 4*(k-1)+1)
            dbl_mb(kvs_new(1) + 4*(k-1)+2) = dbl_mb(kvs(1) + 4*(k-1)+2)
            dbl_mb(kvs_new(1) + 4*(k-1)+3) = dbl_mb(kvs(1) + 4*(k-1)+3)
            end do
            if (.not.BA_free_heap(kvs(2)))
     >      call errquit(
     >          'band_path_set: heap failed 2', 0, MA_ERR)
           end if
           dbl_mb(kvs_new(1) + 4*(nk-1))   = kvector(1)
           dbl_mb(kvs_new(1) + 4*(nk-1)+1) = kvector(2)
           dbl_mb(kvs_new(1) + 4*(nk-1)+2) = kvector(3)
           dbl_mb(kvs_new(1) + 4*(nk-1)+3) = kvector(4)
           kvs(1) = kvs_new(1)
           kvs(2) = kvs_new(2)
       end do 
      end do

      return
      end


      subroutine band_path_plist(pcount,plist)
      implicit none
      integer pcount,plist(*)

#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"

*     *** local variables ***
      integer ind
      character*50 test

      integer num_points
      parameter (num_points = 24)
      character*6 points(num_points)
      data points / 'gamma',
     >              'a','b','c','d','e',
     >              'f','g','h','k','l','m',
     >              'n','p','q','r','s',
     >              't','u','v','w','x',
     >              'y','z'/

      pcount = 0
      do while (inp_a(test)) 

        if (inp_match(num_points, .false., test, points, ind)) then
          pcount = pcount + 1
          plist(pcount) = ind
        end if

      end do

      return
      end



      subroutine band_path_point(cell,p,kvector,goodpoint)
      implicit none
      integer cell,p
      real*8 kvector(4)
      logical goodpoint

#include "errquit.fh"

      goodpoint = .false.
      kvector(4) = -1.0d0

      goto ( 100, 200, 300, 400, 500, 600, 700, 800) cell
      call errquit(
     >      'band_path_point: unimplemented cell',cell,INPUT_ERR)

c
c sc, special points gamma, m, r,x
c
  100 if (p.eq.1) then     !gamma
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.12) then    !m
        kvector(1) = 0.5d0
        kvector(2) = 0.5d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.16) then    !r
        kvector(1) = 0.5d0
        kvector(2) = 0.5d0
        kvector(3) = 0.5d0
        goodpoint = .true.
      end if
      if (p.eq.22) then    !x
        kvector(1) = 0.0d0
        kvector(2) = 0.5d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if

      goto 9999


c
c fcc, special points gamma,k,l,u,w,x
c
  200 if (p.eq.1) then     !gamma
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.10) then    !k
        kvector(1) = 0.375d0
        kvector(2) = 0.375d0
        kvector(3) = 0.750d0
        goodpoint = .true.
      end if
      if (p.eq.11) then   !l
        kvector(1) = 0.5d0
        kvector(2) = 0.5d0
        kvector(3) = 0.5d0
        goodpoint = .true.
      end if
      if (p.eq.19) then   !u
        kvector(1) =  0.625d0
        kvector(2) =  0.250d0
        kvector(3) =  0.625d0
        goodpoint = .true.
      end if
      if (p.eq.21) then   !w
        kvector(1) = 0.50d0
        kvector(2) = 0.25d0
        kvector(3) = 0.75d0
        goodpoint = .true.
      end if
      if (p.eq.22) then   !x
        kvector(1) = 0.5d0
        kvector(2) = 0.0d0
        kvector(3) = 0.5d0
        goodpoint = .true.
      end if

      goto 9999

c
c bcc, special points gamma,h,n,p
c
  300 if (p.eq.1) then     !gamma
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.9) then   !h
        kvector(1) =  0.5d0
        kvector(2) = -0.5d0
        kvector(3) =  0.5d0
        goodpoint = .true.
      end if
      if (p.eq.13) then   !n
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.5d0
        goodpoint = .true.
      end if
      if (p.eq.14) then   !p
        kvector(1) = 0.25d0
        kvector(2) = 0.25d0
        kvector(3) = 0.25d0
        goodpoint = .true.
      end if

      goto 9999

c
c rhombohedral
c
  400 write(*,*) "Rhombohedral path not implemented"
      goto 9999

c
c hexagonal, special points gamma,a,h,k,l,m
c
  500 if (p.eq.1) then     !gamma
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.2) then   !a
        kvector(1) = 0.00d0
        kvector(2) = 0.00d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.9) then   !h
        kvector(1) = 1.0d0/3.0d0
        kvector(2) = 1.0d0/3.0d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.10) then   !k
        kvector(1) = 1.0d0/3.0d0
        kvector(2) = 1.0d0/3.0d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.11) then   !l
        kvector(1) = 0.50d0
        kvector(2) = 0.00d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.12) then   !m
        kvector(1) = 0.50d0
        kvector(2) = 0.00d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      goto 9999

c
c simple tetragonal, special points gamma,a,m,r,x,z
c
  600 if (p.eq.1) then     !gamma
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.2) then   !a
        kvector(1) = 0.50d0
        kvector(2) = 0.50d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.12) then   !m
        kvector(1) = 0.50d0
        kvector(2) = 0.50d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.16) then   !r
        kvector(1) = 0.00d0
        kvector(2) = 0.50d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.22) then   !x
        kvector(1) = 0.00d0
        kvector(2) = 0.50d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.24) then   !z
        kvector(1) = 0.00d0
        kvector(2) = 0.00d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if

      goto 9999
c
c simple orthorhombic, special points gamma,r,s,t,u,x,y,z
c
  700 if (p.eq.1) then     !gamma
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.16) then   !r
        kvector(1) = 0.50d0
        kvector(2) = 0.50d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.17) then   !s
        kvector(1) = 0.50d0
        kvector(2) = 0.50d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.18) then   !t
        kvector(1) = 0.00d0
        kvector(2) = 0.50d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.19) then   !u
        kvector(1) = 0.50d0
        kvector(2) = 0.00d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if
      if (p.eq.22) then   !x
        kvector(1) = 0.50d0
        kvector(2) = 0.00d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.23) then   !y
        kvector(1) = 0.00d0
        kvector(2) = 0.50d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.24) then   !z
        kvector(1) = 0.00d0
        kvector(2) = 0.00d0
        kvector(3) = 0.50d0
        goodpoint = .true.
      end if

      goto 9999
c
c bct body-centerd tetragonal, special points gamma,m,n,p,x
c
  800 if (p.eq.1) then     !gamma
        kvector(1) = 0.0d0
        kvector(2) = 0.0d0
        kvector(3) = 0.0d0
        goodpoint = .true.
      end if
      if (p.eq.12) then   !m
        kvector(1) = 0.00d0
        kvector(2) = 0.00d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.13) then   !n
        kvector(1) = 0.00d0
        kvector(2) = 0.00d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.14) then   !p
        kvector(1) = 0.00d0
        kvector(2) = 0.00d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if
      if (p.eq.22) then   !x
        kvector(1) = 0.00d0
        kvector(2) = 0.00d0
        kvector(3) = 0.00d0
        goodpoint = .true.
      end if

      goto 9999


 9999 continue



      return
      end
