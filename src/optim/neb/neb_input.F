c $Id$

*     *************************************
*     *                                   *
*     *           neb_input               *
*     *                                   *
*     *************************************

      subroutine neb_input(rtdb)
      implicit none
      integer rtdb

#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"

*     ***** local variables ****
      logical impose,frze
      integer ind,i
      character*255 test,test2
      integer nbeads
      double precision timestep,kbeads
      character*255 geom_name

*     **** external functions ****
      logical  neb_parse_boolean,neb_setgeom_xyzfilename
      external neb_parse_boolean,neb_setgeom_xyzfilename
      integer  neb_nbeads_xyzfilename
      external neb_nbeads_xyzfilename
      character*7 bead_index_name
      external    bead_index_name

c
 100  if (.not. inp_read()) call errquit('neb_input: inp_read failed',0,
     &       0)
c
      if (.not. inp_a(test))
     $     call errquit('neb_input: failed to read keyword', 0,
     &       0)
c     
c     print level is read
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, "neb")
c
      goto 100

c
c     algorihm - 
c     ----------------------------
      else if (inp_compare(.false.,'algorithm', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: algorithm',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:algorithm',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',1,
     &     RTDB_ERR)
      goto 100

c
c     nhist = m - number histories
c     ----------------------------
      else if (inp_compare(.false.,'nhist', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: algorithm',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:m',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',2,
     &     RTDB_ERR)
      goto 100


c
c     stepsize
c     ----------------------------
      else if (inp_compare(.false.,'stepsize', test)) then
       if(.not. inp_f(timestep))
     &    call errquit('neb_input: timestep',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:stepsize',mt_dbl,1,timestep))
     $     call errquit('neb_input: rtdb put failed',3,
     &     RTDB_ERR)

      goto 100

c
c     trust
c     ----------------------------
      else if (inp_compare(.false.,'trust', test)) then
       if(.not. inp_f(timestep))
     &    call errquit('neb_input: trust',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:trust',mt_dbl,1,timestep))
     $     call errquit('neb_input: rtdb put failed',3,
     &     RTDB_ERR)

      goto 100

c
c     maxiter
c     ----------------------------
      else if (inp_compare(.false.,'maxiter', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: maxiter',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:steps',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',4,
     &     RTDB_ERR)

      goto 100

c
c     nbeads - number of beads
c     ----------------------------
      else if (inp_compare(.false.,'nbeads', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: nbeads',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:nbeads',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',5,
     &     RTDB_ERR)

      goto 100

c
c     kbeads - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'kbeads', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: kbeads',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:kbeads',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',6,
     &     RTDB_ERR)

      goto 100

c
c     reset  
c     ----------------------------
      else if (inp_compare(.false.,'reset', test)) then
       if (.not. rtdb_put(rtdb,'neb:nebnew',mt_log,1,.true.))
     $     call errquit('neb_input: rtdb put failed',7,
     &     RTDB_ERR)

      goto 100

c
c     gmax - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'gmax', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: gmax',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',8,
     &     RTDB_ERR)
      goto 100
c
c     grms - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'grms', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: gmax',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',9,
     &     RTDB_ERR)
      goto 100
c
c     xmax - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'xmax', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: xmax',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',10,
     &     RTDB_ERR)
      goto 100
c
c     xrms - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'xrms', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: xrms',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',11,
     &     RTDB_ERR)
      goto 100

c
c    loose - 
c     ----------------------------
      else if (inp_compare(.false.,'loose', test)) then
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,0.0045d0))
     >     call errquit('neb_input: rtdb put failed',12,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,0.0030d0))
     >     call errquit('neb_input: rtdb put failed',13,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,0.0054d0))
     >     call errquit('neb_input: rtdb put failed',14,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,0.0036d0))
     >     call errquit('neb_input: rtdb put failed',15,
     >     RTDB_ERR)
      goto 100
c
c    default -
c     ----------------------------
      else if (inp_compare(.false.,'default', test)) then
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,0.00045d0))
     >     call errquit('neb_input: rtdb put failed',16,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,0.00030d0))
     >     call errquit('neb_input: rtdb put failed',17,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,0.0018d0))
     >     call errquit('neb_input: rtdb put failed',18,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,0.0012d0))
     >     call errquit('neb_input: rtdb put failed',19,
     >     RTDB_ERR)
      goto 100
c
c    tight -
c     ----------------------------
      else if (inp_compare(.false.,'tight', test)) then
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,0.00015d0))
     >     call errquit('neb_input: rtdb put failed',20,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,0.00001d0))
     >     call errquit('neb_input: rtdb put failed',21,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,0.00006d0))
     >     call errquit('neb_input: rtdb put failed',22,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,0.00004d0))
     >     call errquit('neb_input: rtdb put failed',23,
     >     RTDB_ERR)
      goto 100

c
c    impose 
c     ----------------------------
      else if (inp_compare(.false.,'impose', test)) then
         impose = .true.
         if(inp_a(test2)) impose = neb_parse_boolean(test2,.true.)
         if (.not. rtdb_put(rtdb,'neb:impose',mt_log,1,impose))
     >     call errquit('neb_input: rtdb put failed',24,
     >     RTDB_ERR)
      
      goto 100

c
c    hasmiddle
c     ----------------------------
      else if (inp_compare(.false.,'hasmiddle', test)) then
         impose = .true.
         if(inp_a(test2)) impose = neb_parse_boolean(test2,.true.)
         if (.not. rtdb_put(rtdb,'neb:hasmiddle',mt_log,1,impose))
     >     call errquit('neb_input: rtdb put failed',25,
     >     RTDB_ERR)
      
      goto 100

c
c    freeze1
c     ----------------------------
      else if (inp_compare(.false.,'freeze1', test)) then
         if (.not. inp_l(frze)) frze = .false.
         if (.not. rtdb_put(rtdb,'neb:freeze1',mt_log,1,frze))
     >      call errquit('neb_input: rtdb put failed',0, RTDB_ERR)

      goto 100

c
c    freezen
c     ----------------------------
      else if (inp_compare(.false.,'freezen', test)) then
         if (.not. inp_l(frze)) frze = .false.
         if (.not. rtdb_put(rtdb,'neb:freezen',mt_log,1,frze))
     >      call errquit('neb_input: rtdb put failed',0, RTDB_ERR)

      goto 100

c
c    xyz_path
c     ----------------------------
      else if (inp_compare(.false.,'xyz_path', test)) then
         if(inp_a(test2)) then

            nbeads = neb_nbeads_xyzfilename(test2)
            if (.not. rtdb_put(rtdb,'neb:nbeads',mt_int,1,nbeads))
     >        call errquit('neb_input: rtdb put failed',5,
     >        RTDB_ERR)

            do i=1,nbeads
               geom_name = 'bead'//bead_index_name(i)//':geom'
               if (.not.neb_setgeom_xyzfilename(rtdb,test2,i,geom_name))
     >           call errquit('neb_input: rtdb put failed',6,
     >           RTDB_ERR)
            end do

            if (.not. rtdb_put(rtdb,'neb:custom_path',mt_log,1,.true.))
     >        call errquit('neb_input: rtdb put failed',5,
     >        RTDB_ERR)

         end if
      
      goto 100

c
c    print_shift
c     ----------------------------
      else if (inp_compare(.false.,'print_shift', test)) then
       if(.not. inp_i(ind))
     &    call errquit('neb_input: print_shift',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:print_shift',mt_int,1,ind))
     $     call errquit('neb_input: rtdb put failed',27,
     &     RTDB_ERR)

      goto 100

c
c     end of input block
c     -------------------
      else if (inp_compare(.false.,'end', test)) then
       goto 101
      endif
c     
c
101   continue
      
      return
      end
c



*     *********************************************
*     *                                           *
*     *          neb_parse_boolean                *
*     *                                           *  
*     *********************************************  
      logical function neb_parse_boolean(word,default)
      implicit none
      character*(*) word
      logical default

#include "inp.fh"
#include "errquit.fh"
      
      logical value
      value = default
      if (inp_compare(.false.,word,'on'))      value = .true.
      if (inp_compare(.false.,word,'yes'))     value = .true.
      if (inp_compare(.false.,word,'.true.'))  value = .true.
      if (inp_compare(.false.,word,'T'))       value = .true.
      if (inp_compare(.false.,word,'off'))     value = .false.
      if (inp_compare(.false.,word,'no'))      value = .false.
      if (inp_compare(.false.,word,'.false.')) value = .false.
      if (inp_compare(.false.,word,'F'))       value = .false.
      if (inp_compare(.false.,word,'fix'))     value = .false.

      neb_parse_boolean = value
      return
      end 



*     ************************************************
*     *                                              *
*     *            neb_setgeom_xyzfilename           *
*     *                                              *
*     ************************************************
      logical function neb_setgeom_xyzfilename(irtdb,in_xyz,nf,geomname)
      implicit none
      integer irtdb
      character*(*) in_xyz
      integer nf
      character*(*) geomname

#include "mafdecls.fh"
#include "errquit.fh"
#include "msgids.fh"
#include "global.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"

c     local variables
      integer ns
      integer i,j
      integer k
      logical otitle
      integer i_t,h_t
      integer i_m,h_m
      integer i_q,h_q
      integer i_ctmp,h_ctmp
      integer atn
      character*32 pname
      character*72 title
      character*16 tag
      character*16 buf
      character*255 filename
      character*255 xyzfile
      character*255 trjfile
      character*255 message

      integer fn_xyz,fn_trj
      logical end_of_file
      integer geom              ! handle for geometry
      real*8 cau2ang

      logical geom_create,geom_print
      external geom_create,geom_print
      logical geom_tag_to_element
      external geom_tag_to_element
      logical geom_cart_set,geom_masses_set
      external geom_cart_set,geom_masses_set
      logical geom_tag_to_default_mass
      external geom_tag_to_default_mass
      logical geom_rtdb_store,geom_destroy
      external geom_rtdb_store,geom_destroy
 

      cau2ang = 0.529177d0
      pname = "neb_setgeom_xyz"
      neb_setgeom_xyzfilename = .false.
c
c     we assume that xyz file has a title
c     -----------------------------------
      otitle = .true.
c
c      if(.not.util_xyz_nframes(in_xyz,fn_xyz))
c     >       call errquit("cannot get number of frames",0,0)
c      write(*,*) "number of frames ", fn_xyz

      xyzfile = in_xyz(1:inp_strlen(in_xyz))
      call util_file_name_resolve(xyzfile, .false.)
c
      filename = xyzfile(1:inp_strlen(xyzfile))

c
c     prepare files for reading/writing
c     ---------------------------------
      if(.not.util_get_io_unit(fn_xyz)) 
     >       call errquit("cannot get file number",0,0)
      filename = xyzfile
      open(fn_xyz,file=filename,form='formatted',status='old',
     $          err=133)
c
c
c     seek frame
c     ----------
      if (nf.gt.1) then
      if(.not.util_xyz_seek(fn_xyz,nf))
     >       call errquit("cannot get frame",0,0)
      end if

c
c     get number of atoms
c     ------------------
      message = " number of atoms "
      read(fn_xyz,*,err=134) ns
c
c     temporary stack memory
c     ----------------------
      if(.not.ma_push_get(mt_byte,16*ns,'t',h_t,i_t))
     + call errquit(pname//'Failed to allocate memory for t',ns,
     &       MA_ERR)

      if(.not.ma_push_get(mt_dbl,3*ns,'ctmp',h_ctmp,i_ctmp))
     + call errquit( pname//'Failed to allocate memory for ctmp',
     + 3*ns, MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'q',h_q,i_q))
     + call errquit(pname//'Failed to allocate memory for q',ns,
     &       MA_ERR)

      if(.not.ma_push_get(mt_dbl,ns,'m',h_m,i_m))
     + call errquit('qmmm: Failed to allocate memory for m',ns,
     &       MA_ERR)


c     read the coords 
c     --------------------------------
      message = " title field"
      if(otitle) 
     +  read(fn_xyz,*,err=134,end=135) title

      do i=1,ns
        tag = " "
        read(fn_xyz,*,err=134,end=135) tag,
     +        (dbl_mb(i_ctmp+3*(i-1)+k-1),k=1,3)
        do j=1,16
            byte_mb(i_t+16*(i-1)+j-1)=tag(j:j)
        end do

        if (.not.
     &      geom_tag_to_default_mass(tag,dbl_mb(i_m+i-1)))
     &      call errquit(pname//'default mass failed',
     &      911, INPUT_ERR)

        if (.not.
     &      geom_tag_to_element(tag,buf,buf,atn))
     &      call errquit(pname//'default atn failed',
     &      911, INPUT_ERR)

        dbl_mb(i_q+i-1)=atn
      end do
      call dscal(3*ns,1.0d0/cau2ang,dbl_mb(i_ctmp),1)



      if (.not. geom_create(geom, geomname)) call errquit
     $     (pname//'geom_create failed !', 0, GEOM_ERR)
c
      if(.not.geom_cart_set(geom,ns,byte_mb(i_t),
     + dbl_mb(i_ctmp),dbl_mb(i_q)))
     + call errquit('qmmm: Failed to initialize geometry',0, GEOM_ERR)
c
      if(.not.geom_masses_set(geom,ns,dbl_mb(i_m)))
     + call errquit('qmmm: Failed to initialize masses',0, GEOM_ERR)
      call geom_compute_values(geom)
c
      if(.not.geom_print(geom))
     +   call errquit('qmmm: Failed to print geom',0, RTDB_ERR)
c
      if(.not.geom_rtdb_store(irtdb,geom,geomname))
     + call errquit('qmmm: Failed to store geom to rtdb',0, RTDB_ERR)

       if(.not.geom_destroy(geom))
     + call errquit('qmmm: Failed to destroy geometry',0, GEOM_ERR)

c
      if(.not.ma_pop_stack(h_m))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_q))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_ctmp))
     & call errquit(pname//' 
     >              Failed to deallocate stack c_tmp',ns,
     &       MA_ERR)

      if(.not.ma_pop_stack(h_t))
     & call errquit(pname//'
     >              Failed to deallocate stack i_itmp',ns,
     &       MA_ERR)

      close(fn_xyz)
      neb_setgeom_xyzfilename = .true.
      return

 133  call errquit(pname//'error opening/closing '//filename,0, 0)
 134  call errquit(pname//'error reading xyz file'//message,0, 0)
 135  call errquit(pname//'error end of file at'//message,0, 0)

      end

*     ************************************************
*     *                                              *
*     *            neb_nbeads_xyzfilename            *
*     *                                              *
*     ************************************************
      integer function neb_nbeads_xyzfilename(xyzfilename)
      implicit none
      character*(*) xyzfilename

#include "inp.fh"
#include "msgtypesf.h"
#include "errquit.fh"

      integer ii,nion,nbeads
      character ch_tmp
      character*255 filename

      filename = xyzfilename(1:inp_strlen(xyzfilename))
      call util_file_name_resolve(filename, .false.)
c

      open(unit=18,file=filename(1:inp_strlen(filename)),
     >     form='formatted',status='old')
      nbeads = 0
      do while(.true.)
         read(18,*,ERR=30,END=30) nion
         read(18,*,ERR=30,END=30) ch_tmp
         do ii=1,nion
            read(18,*,ERR=30,END=30) ch_tmp
         end do
         nbeads = nbeads + 1
      end do
  30  continue
      close(18)

      neb_nbeads_xyzfilename = nbeads
      return
      end

