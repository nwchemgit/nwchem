      subroutine simulation_cell_input(rtdb)
c
c $Id$
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      double precision rad
      parameter(rad=57.295779513082343d0)


      logical bset,lset,nset,nset_small
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 19)
      character*22 dirs(num_dirs)
      character*255 test

      character*50 rtdb_name
      character*50 name,boundry
      double precision unita(3,3)
      double precision lat_a,lat_b,lat_c
      double precision alpha,beta,gamma
      double precision scale,ll

      logical box_orient,b1set,b2set,b3set
      integer box_type
      double precision box_delta

      integer ntokens
      parameter (ntokens=7)
      character*10 tokens(ntokens)

      integer ngrid(3),ngrid_small(3),l
      data dirs / 'cell_name:',
     >            'boundry_conditions:',
     >            'lattice_vectors:',
     >            'ngrid:',
     >            'lattice',
     >            'cell_name',
     >            'boundary_conditions',
     >            'lattice_vectors',
     >            'ngrid',
     >            'boundary_conditions:',
     >            'boundry_conditions',
     >            'sc',
     >            'fcc',
     >            'bcc',
     >            'box_delta',
     >            'box_orient',
     >            'box_different_lengths',
     >            'ngrid_small',
     >            'end'/

      data tokens /'lat_a', 'lat_b', 'lat_c',
     &     'alpha', 'beta', 'gamma', 'end'/
c

c
c     
c     **** intialize stuff ****
      call get_scalefrominput(scale)
      ngrid(1) = 32
      ngrid(2) = 32
      ngrid(3) = 32
      unita(1,1) = 0.0d0
      unita(2,1) = 0.0d0
      unita(3,1) = 0.0d0
      unita(1,2) = 0.0d0
      unita(2,2) = 0.0d0
      unita(3,2) = 0.0d0
      unita(1,3) = 0.0d0
      unita(2,3) = 0.0d0
      unita(3,3) = 0.0d0
      boundry = 'periodic'
      name    = 'cell_default'
      box_type = 0
      box_orient = .false.
      box_type = 0

      nset = .false.
      lset = .false.
      bset = .false.
      b1set = .false.
      b2set = .false.
      b3set = .false.
      nset_small = .false.
c
c
c
 10   if (.not. inp_read()) 
     >     call errquit(
     >           'simulation_cell_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'simulation_cell_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'simulation_cell_input: unknown directive', 0,
     &       INPUT_ERR)


      goto ( 100,  200, 300, 400, 500, 
     >       100,  200, 300, 400,
     >       200,  200,
     >       600, 700, 800, 
     >       1000, 1100, 1200, 1300,
     >      9999) ind
      call errquit(
     >      'simulation_cell_input: unimplemented directive', ind,
     &       INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(name))
     >     call errquit(
     >          'simulation_cell_input: failed to read keyword', 0,
     &       INPUT_ERR)
      goto 10


c
c boundry_condition
c
 200  if (.not. inp_a(boundry))
     >     call errquit(
     >          'simulation_cell_input: failed to read keyword', 0,
     &       INPUT_ERR)
      bset = .true.
      goto 10


c
c  lattice_vectors
c 
 300  if (.not. inp_read()) 
     >     call errquit(
     >     'simulation_cell_input: inp_read(unita(*,1)) failed', 0,
     &       INPUT_ERR)
      if (.not. inp_f(unita(1,1)))
     >     call errquit(
     >           'simulation_cell_input: failed to read unita', 1,
     &       INPUT_ERR)
      if (.not. inp_f(unita(2,1)))
     >     call errquit(
     >          'simulation_cell_input: failed to read unita', 2,
     &       INPUT_ERR)
      if (.not. inp_f(unita(3,1)))
     >     call errquit(
     >          'simulation_cell_input: failed to read unita', 3,
     &       INPUT_ERR)

      if (.not. inp_read()) 
     >     call errquit(
     >     'simulation_cell_input: inp_read(unita(*,2)) failed', 4,
     &       INPUT_ERR)
      if (.not. inp_f(unita(1,2)))
     >     call errquit(
     >          'simulation_cell_input: failed to read unita', 5,
     &       INPUT_ERR)
      if (.not. inp_f(unita(2,2)))
     >     call errquit(
     >           'simulation_cell_input: failed to read unita', 6,
     &       INPUT_ERR)
      if (.not. inp_f(unita(3,2)))
     >     call errquit(
     >           'simulation_cell_input: failed to read unita', 7,
     &       INPUT_ERR)


      if (.not. inp_read()) 
     >     call errquit(
     >     'simulation_cell_input: inp_read(unita(*,3)) failed', 8,
     &       INPUT_ERR)
      if (.not. inp_f(unita(1,3)))
     >     call errquit(
     >           'simulation_cell_input: failed to read unita', 9,
     &       INPUT_ERR)
      if (.not. inp_f(unita(2,3)))
     >     call errquit(
     >           'simulation_cell_input: failed to read unita', 10,
     &       INPUT_ERR)
      if (.not. inp_f(unita(3,3)))
     >     call errquit(
     >           'simulation_cell_input: failed to read unita', 11,
     &       INPUT_ERR)
      lset = .true.
      goto 10


c
c  ngrid
c
 400  if (.not. inp_i(ngrid(1)))
     >     call errquit(
     >          'simulation_cell_input: failed to read ngrid', 0,
     &       INPUT_ERR)
      if (.not. inp_i(ngrid(2)))
     >     call errquit(
     >          'simulation_cell_input: failed to read ngrid', 0,
     &       INPUT_ERR)
      if (.not. inp_i(ngrid(3)))
     >     call errquit(
     >          'simulation_cell_input: failed to read ngrid', 0,
     &       INPUT_ERR)
      nset = .true.
      goto 10

c
c  lattice
c 
 500  lat_a = 20.0d0
      lat_b = 20.0d0
      lat_c = 20.0d0
      alpha = 90.0d0
      beta  = 90.0d0
      gamma = 90.0d0

c     ***** Read new line of input *****
   34 if (.not. inp_read()) call errquit
     $     ('simulation_cell_input: premature end of file', 0,
     &       INPUT_ERR)

   33 if (.not. inp_a(test)) goto 34
     
      if (.not. inp_match(ntokens, .false., test, tokens, ind))
     >     call errquit(
     >           'simulation_cell_input: unknown directive', 0,
     &       INPUT_ERR)

          goto (11,12,13,14,15,16,
     >          35) ind
c
  11      if (.not. inp_f(lat_a)) call errquit
     >    ('simulation_cell_input: error reading lat_a', 0,
     &       INPUT_ERR)
          goto 33
  12      if (.not. inp_f(lat_b)) call errquit
     >    ('simulation_cell_input: error reading lat_b', 0,
     &       INPUT_ERR)
          goto 33
  13      if (.not. inp_f(lat_c)) call errquit
     >    ('simulation_cell_input: error reading lat_c', 0, INPUT_ERR)
          goto 33
  14      if (.not. inp_f(alpha)) call errquit
     >    ('simulation_cell_input: error reading alpha', 0, INPUT_ERR)
          goto 33
  15      if (.not. inp_f(beta)) call errquit
     >    ('simulation_cell_input: error reading beta', 0, INPUT_ERR)
          goto 33
  16      if (.not. inp_f(gamma)) call errquit
     >    ('simulation_cell_input: error reading gamma', 0, INPUT_ERR)
          goto 33

  35  alpha = alpha/rad
      beta  = beta/rad
      gamma = gamma/rad
      call generate_amatrix(lat_a,lat_b,lat_c,
     >                      alpha,beta,gamma,
     >                      unita)
      lset = .true.
      goto 10

c
c  simple cubic lattice vectors
c 
 600  if (.not. inp_f(ll)) 
     >     call errquit(
     >     'simulation_cell_input: SC inp_read(L) failed', 0, INPUT_ERR)
      unita(1,1) = ll
      unita(2,1) = 0.0d0
      unita(3,1) = 0.0d0
      unita(1,2) = 0.0d0
      unita(2,2) = ll
      unita(3,2) = 0.0d0
      unita(1,3) = 0.0d0
      unita(2,3) = 0.0d0
      unita(3,3) = ll
      lset = .true.
      goto 10

c
c  face-centered cubic lattice vectors
c 
 700  if (.not. inp_f(ll)) 
     >     call errquit(
     >     'simulation_cell_input: FCC inp_read(L) failed', 0,
     &       INPUT_ERR)
      unita(1,1) = ll/2.0d0
      unita(2,1) = ll/2.0d0
      unita(3,1) = 0.0d0
      unita(1,2) = ll/2.0d0
      unita(2,2) = 0.0d0
      unita(3,2) = ll/2.0d0
      unita(1,3) = 0.0d0
      unita(2,3) = ll/2.0d0
      unita(3,3) = ll/2.0d0
      lset = .true.
      goto 10

c
c  body-centered cubic lattice vectors
c 
 800  if (.not. inp_f(ll)) 
     >     call errquit(
     >     'simulation_cell_input: BCC inp_read(L) failed', 0,
     &       INPUT_ERR)
      unita(1,1) = -ll/2.0d0
      unita(2,1) =  ll/2.0d0
      unita(3,1) =  ll/2.0d0
      unita(1,2) =  ll/2.0d0
      unita(2,2) = -ll/2.0d0
      unita(3,2) =  ll/2.0d0
      unita(1,3) =  ll/2.0d0
      unita(2,3) =  ll/2.0d0
      unita(3,3) = -ll/2.0d0
      lset = .true.
      goto 10

c
c  box_delta
c
 1000 if (.not. inp_f(box_delta))
     >     call errquit(
     >     'simulation_cell_input: box_delta failed', 0,
     &       INPUT_ERR)
      b1set = .true.
      goto 10

c
c  box_orient
c
 1100 box_orient=.true.
      b2set = .true.
      goto 10

c
c  box_different_lengths
c
 1200 box_type = 1
      b3set = .true.
      goto 10

c
c  ngrid_small
c
 1300 if (.not. inp_i(ngrid_small(1)))
     >     call errquit(
     >          'simulation_cell_input: failed to read ngrid_small', 0,
     &       INPUT_ERR)
      if (.not. inp_i(ngrid_small(2)))
     >     call errquit(
     >          'simulation_cell_input: failed to read ngrid_small', 0,
     &       INPUT_ERR)
      if (.not. inp_i(ngrid_small(3)))
     >     call errquit(
     >          'simulation_cell_input: failed to read ngrid_small', 0,
     &       INPUT_ERR)
      nset_small = .true.
      goto 10



*     ***** add simulation_cell to rtdb ****
 9999 continue
      if (name.eq.' ') 
     >     call errquit('simulation_cell_input: nead a name', 0,
     &       INPUT_ERR)
      l = index(name,' ') -1

      if (bset) then
      rtdb_name = name(1:l)//':boundry'
      ind = index(boundry,' ') - 1
      if (.not. rtdb_cput(rtdb,rtdb_name,1,boundry(1:ind)))
     >     call errquit(
     >     'simulation_cell_input: rtdb_cput(boundry) failed', 0,
     &       INPUT_ERR)
      end if

*     **** convert to bohrs using scale variable ****
      unita(1,1) = unita(1,1) * scale
      unita(2,1) = unita(2,1) * scale
      unita(3,1) = unita(3,1) * scale

      unita(1,2) = unita(1,2) * scale
      unita(2,2) = unita(2,2) * scale
      unita(3,2) = unita(3,2) * scale

      unita(1,3) = unita(1,3) * scale
      unita(2,3) = unita(2,3) * scale
      unita(3,3) = unita(3,3) * scale

      if (lset) then
      rtdb_name = name(1:l)//':unita'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,9,unita))
     >     call errquit(
     >     'simulation_cell_input: rtdb_put(unita) failed', 0,
     &       INPUT_ERR)
      end if

      if (nset) then
      rtdb_name = name(1:l)//':ngrid'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,3,ngrid))
     >     call errquit(
     >     'simulation_cell_input: rtdb_put(ngrid) failed', 0,
     &       INPUT_ERR)
      end if

      if (nset_small) then
      rtdb_name = name(1:l)//':ngrid_small'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,3,ngrid_small))
     >     call errquit(
     >     'simulation_cell_input: rtdb_put(ngrid_small) failed', 0,
     &       INPUT_ERR)
      end if



      if (b1set) then
      box_delta = box_delta*scale
      rtdb_name = name(1:l)//':box_delta'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_dbl,1,box_delta))
     >     call errquit(
     >     'simulation_cell_input: rtdb_put(box_delta) failed', 0,
     &       INPUT_ERR)
      end if
      if (b2set) then
      rtdb_name = name(1:l)//':box_orient'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_log,1,box_orient))
     >     call errquit(
     >     'simulation_cell_input: rtdb_put(box_orient) failed', 0,
     &       INPUT_ERR)
      end if
      if (b3set) then
      rtdb_name = name(1:l)//':box_type'
      if (.not. rtdb_put(rtdb,rtdb_name,mt_int,1,box_type))
     >     call errquit(
     >     'simulation_cell_input: rtdb_put(box_type) failed', 0,
     &       INPUT_ERR)
      end if


      

      end

*     *************************************************
*     *                                               *
*     *               get_scalefrominput              *
*     *                                               *
*     *************************************************
      subroutine get_scalefrominput(scale)
      implicit none
      double precision scale

#include "inp.fh"
#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ***
      double precision bohr
      parameter (bohr=0.52917715d0) 

      integer ind
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 1)
      character*22 dirs(num_dirs)
      character*255 test

      integer nunits
      parameter (nunits=9)
      character*10 units(nunits)

      data dirs / 'units' /

      data units / 'au','atomic','bohr',
     >             'angstroms', 'an', 
     >             'nanometers','nm',
     >             'picometers','pm' /


      scale = 1.0d0
 51   if (inp_a(test)) then
      if (inp_match(num_dirs, .false., test, dirs, ind)) then

         if (.not. inp_a(test)) call errquit
     >      ('simulation_cell_input: error reading units', 0,
     &       INPUT_ERR)

        if (.not.inp_match(nunits, .false., test, units, ind)) 
     >    call errquit('simulation_cell_input: error reading units', 0,
     &       INPUT_ERR)


        goto ( 52,52,52,53,53,54,54,55,55 ) ind

*       **** atomic units ****
  52    scale = 1.0d0
        goto 60

*       **** angstrom units ****
  53    scale = 1.0d0/bohr
        goto 60

*       **** nanometer units ****
  54    scale = 10.0d0/bohr
        goto 60

*       **** picometer units ****
  55    scale = 1.0d0/bohr/100.0d0
        goto 60

      end if
      end if

  60  continue
      return
      end
     


*     *************************************************
*     *                                               *
*     *               generate_amatrix                *
*     *                                               *
*     *************************************************

*  Adapted from geom_3d code.
c
c-----------> build a-matrix
c
c Chosen convention of two axis systems.
c    Orient oblique crystal system to external rectilinear laboratory
c    coordinate system such that
c    b-axis is along j (unit rectilinear vector), ab-plane contains i (rect),
c    i.e., j along b, i in  ab-plane
c
c  Note this matrix is implicitly built to convert the input coordinates
c  from angstroms to bohr; i.e., amat is in atomic units.
c
c  the Jth column of this matrix contains the jth lattice vector
c aij is the ith Cartesian component of the jth lattice vector.
c
c
c     Entry - lat_a, lat_b, and lat_c, lattice constants in au
c             alpha,beta,gamma, lattice angles in radians
c
c     Exit  - lattice_a, lattice vectors
c

      subroutine generate_amatrix(lat_a,lat_b,lat_c,
     >                       alpha,beta,gamma,
     >                       lattice_a)
      double precision lat_a,lat_b,lat_c
      double precision alpha,beta,gamma
      double precision lattice_a(3,3)

*     **** local variables ****
      integer i,j,iang
      double precision c1,c2,c3,s3
      double precision cdist(3),cang(3)
      double precision gmat(3,3)


*     **** external functions ****
      double precision deter3
      external deter3


      cdist(1) = lat_a
      cdist(2) = lat_b
      cdist(3) = lat_c
      cang(1) = alpha
      cang(2) = beta
      cang(3) = gamma

c    **** build the metrical matrix (atomic units) ****
      do i=1,3
        gmat(i,i)=cdist(i)**2
      end do
      iang=3
      do i=1,3
        do j=i+1,3
          gmat(i,j)=cdist(i)*cdist(j)*dcos(cang(iang))
          gmat(j,i)=gmat(i,j)
          iang=iang-1
        end do  
      end do


*     **** get volume ***
      vol=dsqrt(deter3(gmat))


*     **** generate lattice_a ***
      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
      s3=dsin(cang(3))
      lattice_a(1,1) = cdist(1)*s3
      lattice_a(1,2) = 0.0d0
      lattice_a(1,3) = (cdist(3)*(c2-c1*c3)/s3)
      lattice_a(2,1) = cdist(1)*c3
      lattice_a(2,2) = cdist(2)
      lattice_a(2,3) = cdist(3)*c1
      lattice_a(3,1) = 0.0d+00
      lattice_a(3,2) = 0.0d+00
      lattice_a(3,3) = (vol/(cdist(1)*cdist(2)*s3))


      return
      end
