      subroutine bsse_input(rtdb)
      implicit none
c
c     This subroutine read the input lines in bsse     
c     to get info about monomers that build the
c     super molecule.  nmon identify monomers, mon(nmon), 
c     each monomer have mon_atm(nmon). Name of monomers is
c     allocated in mon_name(nmon)
c
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "geom.fh"
c here are defined common variables
#include "bsse_common.fh"      
c
      integer rtdb              ! [input]
c
      integer geom
      integer atm_tot           ! to check total atoms
      integer qtot              ! to print total charge
      integer nopen             ! to print spin multiplicity
      integer i, j,l,k
      integer nfield
      integer nopt
      parameter (nopt = 9)
      integer ind
c
      character*80 buf
      character*255 field
      character*18 opt (nopt)
c
      logical status
      logical do_bsse
c
c:
      logical bsse_rtdb_store
      external bsse_rtdb_store
      logical bsse_rtdb_load
      external bsse_rtdb_load
c
      data opt /'end', 'on', 'off','tidy','charge', 'input',
     &  'input_wghost', 'mon', 'mult'/

c
c     ------------------welcome------------------------
c
      buf = ' '
      write(buf,*) ' Input BSSE Module - Counter Poise Approach'
      write(LuOut,*)
      write(LuOut,*)
      call util_print_centered(LuOut,buf,40,.true.)
      write(LuOut,*)
c
c     -------------------------------------------------
c     -----   get info supermolecule geometry   -------
c     -------------------------------------------------
c
      if (.not. rtdb_cget(rtdb,'geometry', 1, spr_name))
     $     spr_name = 'geometry'
c

      if (.not. geom_create(geom, spr_name)) 
     $ call errquit('bsse_input: geom_create failed !', 0,GEOM_ERR)
c

      if (.not. geom_rtdb_load(rtdb, geom, spr_name))
     $ call errquit('bsse_input: no geometry load form rtdb', 0,
     $        GEOM_ERR)

      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('bsse_input: geom_ncent ?', 0, GEOM_ERR)
c     -------------------------------------------------
c     -----        reading  input  fields      --------
c     -------------------------------------------------
c

      call inp_set_field(0)     ! goto the begin of line

c
c:    preliminaries
c
      qtot    = 0
      atm_tot = 0
      nmon    = 0 
      do_bsse = .true.

c

      call ifill(mx_atm,1,mmon,1)     ! multiplicity default
      call dfill(mx_atm,0.0d0,qmon,1) ! charge default
      call cfill(mx_atm,' ',input,1) ! charge default
c
 100  if(.not.inp_read())
     $  call errquit('bsse_input: unexpected eof ',911,INPUT_ERR )
c
      nfield = inp_n_field()
c
 150  if (.not.inp_a(field))
     $  call errquit('bsse_input: failed to read field',911,INPUT_ERR )
c
      if (.not. inp_match(nopt, .false., field, opt, ind)) 
     $    goto 10 

      goto (900,  850, 800, 700, 600, 500, 400, 300, 200) ind

c
c:    none 
c

  10  write(LuOut,20)

  20  format

     $(/' valid bsse structure input : '/
     $  ' bsse                     '/
     $  '   mon <character name fragment> <integer list atoms>'/    
     $  '     input <input line>  '/
     $  '     input_wghost <input line>  '/
     $  '     charge <double charge> '/
     $  '     mult <integer multiplicity> '/
     $  '   off                    '/
     $  '   on                     '/
     $  ' end                      '/
     $/)
c
      call errquit('bsse_input: unknown directive', 911,INPUT_ERR)
c
c
c:    mon
c

 300  continue 
c
      nmon  = nmon + 1
c   
      mon_name(nmon) = field  ! took the first field however use to name
c
      if (.not.inp_a(mon_name(nmon)))
     $  call errquit
     $          ('bsse_input: failed to read name field',911,INPUT_ERR)
c
c     Read the atom numbers and count the number of atoms as we go along.
c     If we read something else than an integer it might be the next
c     keyword on the line. So, leave the loop, and check if we have reached
c     the end of the line. If we are not at the end of the line goto 150
c     to read the next keyword, otherwise goto 100 to read the next line.
c
      i = 0
      do while (inp_i(mon(nmon,i+1)))
        i = i + 1
        mon_atm(nmon) = i
        atm_tot = atm_tot + 1
      enddo
      if (inp_cur_field().lt.nfield) goto 150
c
      go to 100
c
c:    mult
c
 200  continue
        if(nmon.eq.0) goto 10
        if (.not. inp_i(mmon(nmon) )) call errquit
     $    ('bsse_input: failed reading monomer multiplicity',
     $     nmon,INPUT_ERR)
        if (mmon(nmon).eq.0) call errquit
     $    ('bsse_input: invalid multiplicity ',mmon(nmon),
     $     INPUT_ERR)
        if (inp_cur_field().lt.nfield) goto 150
      goto 100
c
c:    input_wghost
c
 400  continue
c
      if(nmon.eq.0) goto 10

      i=(nmon)*2

      if (.not. inp_a(input(i) )) call errquit
     $  ('bsse_input: failed reading input [input]',911,INPUT_ERR)
      if (inp_cur_field().lt.nfield) goto 150

      go to 100
c
c:    input
c
 500  continue
      if(nmon.eq.0) goto 10

      i=(nmon-1)*2+1

      if (.not. inp_a(input(i))) 
     $  call errquit
     $     ('bsse_input: failed reading input [input]',911,INPUT_ERR)
c
      if (inp_cur_field().lt.nfield) goto 150
c
      go to 100
c

 600  continue
c
c:    charge
c
      if(nmon.eq.0) goto 10
c
      if(.not. inp_f( qmon(nmon)))
     $  call errquit('bsse_input: reading monomer charge',911,INPUT_ERR)
c
      if (inp_cur_field().lt.nfield) goto 150

      go to 100
c
c:    tidy
c       clean database of any bsse info
  700 continue
c
c:    off
c
      if(.not. rtdb_delete(rtdb, 'bsse'))
     $  call errquit('bsse_input: cannot clean database',911,RTDB_ERR)
      goto 100
c 
  800 continue

      do_bsse=.false.
      buf = ' '
      write(buf,*) ' Any BSSE operations are off'
      write(LuOut,*)
      write(LuOut,*)
      call util_print_centered(LuOut,buf,40,.true.)
      write(LuOut,*)

      goto 100
c     
c
  850 continue
  
      do_bsse=.true.  
      if(.not.bsse_rtdb_load(rtdb))
     $  call errquit('bsse_input: load data input in db',911,RTDB_ERR)
c 
      atm_tot= natoms
c
      goto 100
c
c:    end
c
  900 continue

      if(do_bsse) then

        goto 1000

      else

        goto 1100

      endif
c
c:    done
c
 1000 continue

c     ------------------------------------------------------
c     check : total atoms
c     ------------------------------------------------------

      if (atm_tot.ne.natoms) 
c     $ goto 10
     $ call errquit
     $ ('bsse_input: number of atoms is wrong',911,INPUT_ERR) 

c     ------------------------------------------------------
c     check : dont repeat atoms
c     ------------------------------------------------------

      do j = 1, nmon

        do i = 1, mon_atm(j)

          do l = j, nmon

            do k = 1, mon_atm(l)

              if((i.eq.k).and.(j.eq.l)) then

                status=.true.

              elseif (mon(j,i).eq.mon(l,k)) then
  
              call errquit
     $       ('bsse_input: there are some atoms repeated',911,INPUT_ERR)
              endif
            enddo
          enddo
        enddo
      enddo

c----------------------------------------------------
c        check the atoms are correct
c----------------------------------------------------
      do j = 1, nmon
        do i = 1, mon_atm(j)

              if( mon(j,i).gt.natoms) then
              call errquit
     $       ('bsse_input: incorrect atom number',mon(j,i),INPUT_ERR)
              endif
        enddo
      enddo

c     ------------------------------------------------------
c     check : total charge
c     ------------------------------------------------------

      do j = 1, nmon

        qtot = qtot + qmon(j)

      enddo
c
c     ------------------------------------------------------
c     check : unpaired electrons
c     ------------------------------------------------------
c
      nopen = 0
      do j = 1, nmon
        nopen = nopen + abs(mmon(j)) - 1
      enddo
c
      write(LuOut, 60)  spr_name, nmon,
     $                  mod(nopen,2)+1, nopen+1, qtot
 60   format(/
     $  '    supermolecule geometry name = ', a50/
     $  '             number of monomers = ', i4/
     $  '             total multiplicity = ', i4,' to ',i4/
     $  '                   total charge = ', i4/)
c

      write(LuOut, *) 
     $  '       atoms for each monomer  ' 

      do j=1, nmon

        write(LuOut, 70) mon_name(j)
 70     format(/ 
     $      '              monomer     ',  a10,' : ',$)
c

        do i=1, mon_atm(j) 

          write(LuOut, 90)  mon(j,i)

 90       format (i3,$) 

        enddo
      enddo

      write(LuOut,*)
      call util_flush(luout)

c
      go to 1100
c
c: db is done and geometry is destroyed
c
c 1100 continue
c

 1100 continue
c
      if (do_bsse) then
c
        if (.not. rtdb_put(rtdb,'bsse',mt_log,1,.true.))
     $    call errquit('bsse_input: rtdb_put failed',911,RTDB_ERR)
c
      else
c
         if (.not. rtdb_put(rtdb,'bsse',mt_log,1,.false.))
     $     call errquit('bsse_input: rtdb_put failed',911,RTDB_ERR)
c
      endif
c
      if(.not. geom_destroy(geom))
     $  call errquit('bsse_input: geom_destroy failed', 911,RTDB_ERR)
c
      if(.not.bsse_rtdb_store(rtdb))
     $  call errquit('bsse_input: store data input in db',911,RTDB_ERR)
c
      return
c
      end
C>
C> \brief Initialize monomer calculation 
C>
C> Initialize the monomer calculation by modifying the contents of
C> the RunTime Data Base for the current calculation.
C>
      subroutine bsse_param(rtdb, mult, charge, j_mon_name, 
     &                      i_input,theory)
      implicit none
      integer rtdb !< [Input] The RTDB handle
      character*(*) j_mon_name !< [Input] Monomer name
      character*(*) i_input !< [Input] Line of input for monomer calculation
      character*(*) theory !< [Input] The theory to apply
      integer mult !< [Input] Monomer spin multiplicity
      double precision charge !< [Input] Monomer charge
      logical first_j
      character*255 vec_dbi, vec_dbo,tmp
      integer lentheo, lenname
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
c
      logical nw_inp_from_character
      external nw_inp_from_character
c
      if (.not. rtdb_get(rtdb,'bsse:first_j',mt_log,1,first_j))
     $  first_j=.true.
c
c     if (.not. rtdb_cget(rtdb,'theory', 1, theory))
c    &  call errquit('bsse_param: get theory',0)
c      
      lenname = inp_strlen(j_mon_name)
c
      lentheo = inp_strlen(theory)
c
c:    multiplicity
c:      density methods
      if ( theory(1:lentheo).eq.'dft' .or.
     $     theory(1:lentheo).eq.'tddft') then
        if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, mult))
     $    call errquit('bsse_param: rtdb_put of mult failed',
     $                 0,RTDB_ERR )
c:      wavefuntion methods

      elseif( theory(1:lentheo).ne.'dft' .and.
     $        theory(1:lentheo).ne.'tddft') then
        if (.not. rtdb_put(rtdb,'scf:nopen', MT_INT, 1, mult-1))
     $    call errquit('bsse_param: rtdb_put of nopen failed',
     $                 0,RTDB_ERR)
      endif
c
      if (charge .ne. -999.0d0) then
c
         if (.not. rtdb_put(rtdb,'charge',mt_dbl,1,charge))
     $        call errquit('bsse_param:setting charge?',911, RTDB_ERR)
c
      end if
c
      if (j_mon_name .ne. ' ') then
c
         if (.not. rtdb_cput(rtdb,'geometry',1,j_mon_name))
     $        call errquit('bsse_param: setting geometry?',911,RTDB_ERR)
c
      end if
c
        if (first_j) then
c
          tmp = 'atomic'
c
        else
c
          tmp = ' '
          tmp = j_mon_name(1:lenname)//'.bsse.movecs'
c
        endif
c MP2 SCF DFT vectors
c
        vec_dbi = ' '
        vec_dbo = ' '
c
c       if(theory(1:lentheo).ne.'scf' .or. theory(1:lentheo).ne.'dft'
c    &    .or. theory(1:lentheo).ne.'mcscf') then
        if(theory(1:lentheo).eq.'dft') then
          vec_dbi= theory(1:lentheo)//':input vectors'
          vec_dbo= theory(1:lentheo)//':output vectors'
c
       elseif (theory(1:lentheo).eq.'mcscf') then
          vec_dbi= theory(1:lentheo)//':input vectors'
          vec_dbo= theory(1:lentheo)//':output vectors'
c
        else
          vec_dbi = 'scf:input vectors'
          vec_dbo = 'scf:output vectors'

        endif
c
c       
        if (.not. rtdb_cput(rtdb, vec_dbi, 1, tmp))
     &    call errquit('bsse_param: input_vectors',0,RTDB_ERR)
c
        tmp = ' '
        tmp = j_mon_name(1:lenname)//'.bsse.movecs'
c
        if (.not. rtdb_cput(rtdb, vec_dbo, 1,tmp))
     &   call errquit('bsse_param: output_vectors',0,RTDB_ERR)
c
c
      if (i_input .ne. ' ') then

         if (.not. nw_inp_from_character(rtdb,i_input))
     $     call errquit('bsse_param: error processing input string',
     &     060,INPUT_ERR)

      endif
c
      end
c
c
      logical function bsse_rtdb_store(rtdb)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
c#include "tcgmsg.fh"
#include "bsse_common.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      
c      the propose is store bsse's sets into db 

      integer rtdb              ![input]
c     integer itmp ,k
      character*255 ctmp
c
      ctmp = 'bsse:natoms'
       if (.not. rtdb_put( rtdb, ctmp, mt_int, 1,natoms ))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
c
      ctmp = 'bsse:nmon'
       if (.not. rtdb_put( rtdb, ctmp, mt_int, 1,nmon ))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
       ctmp = 'bsse:mon_name'
       if(.not.rtdb_cput(rtdb, ctmp, nmon,  mon_name ))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c     
       ctmp = 'bsse:spr_name'
       if(.not.rtdb_cput(rtdb, ctmp, 1,  spr_name ))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
      ctmp = 'bsse:mon_atm'
       if (.not. rtdb_put( rtdb, ctmp, mt_int, nmon, mon_atm ))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
      ctmp = 'bsse:mon'
       if (.not. rtdb_put( rtdb, ctmp, mt_int,mx_atm*mx_atm,  mon))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
      ctmp = 'bsse:qmon'
       if (.not. rtdb_put( rtdb, ctmp, mt_dbl ,nmon, qmon))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
      ctmp = 'bsse:mmon'
       if (.not. rtdb_put( rtdb, ctmp, mt_int ,nmon, mmon))
     $   call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
      ctmp= 'bsse:input'
       if(.not.rtdb_cput(rtdb, ctmp, nmon*2, input))
     $    call errquit('bsse_rtdb_store: rtdb_put failed',0,RTDB_ERR)
c
      bsse_rtdb_store = .true.
      return
      end
c
      logical function bsse_rtdb_load(rtdb)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
c#include "tcgmsg.fh"
#include "bsse_common.fh"
#include "mafdecls.fh"
      integer rtdb              ![input]
      character*255 ctmp
c
c      the propose is load bsse's sets from db 
c
      ctmp = 'bsse:natoms'
       if (.not.rtdb_get( rtdb, ctmp, mt_int, 1, natoms))
     $   call errquit('bsse_rtdb_load: rtdb_put failed',0,RTDB_ERR)
c
      ctmp = 'bsse:nmon'
       if(.not.rtdb_get( rtdb, ctmp, mt_int,   1 , nmon))
     $   call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
        ctmp = 'bsse:mon_name'
        if(.not.rtdb_cget( rtdb, ctmp, nmon,   mon_name))
     $    call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
        ctmp = 'bsse:spr_name'
        if(.not.rtdb_cget( rtdb, ctmp, 1,   spr_name))
     $    call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
      ctmp = 'bsse:mon_atm'
       if(.not.rtdb_get( rtdb, ctmp, mt_int, nmon,mon_atm))
     $   call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
      ctmp = 'bsse:mon'
       if (.not.rtdb_get( rtdb, ctmp, mt_int, mx_atm*mx_atm,mon))
     $   call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
      ctmp = 'bsse:qmon'
       if (.not.rtdb_get( rtdb, ctmp, mt_dbl, nmon,    qmon ))
     $   call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
      ctmp = 'bsse:mmon'
       if (.not. rtdb_get( rtdb, ctmp, mt_int ,nmon, mmon))
     $   call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
      ctmp = 'bsse:input'
        if(.not.rtdb_cget( rtdb, ctmp, nmon*2,  input))
     $    call errquit('bsse_rtdb_load: rtdb_get failed',0,RTDB_ERR)
c
      bsse_rtdb_load = .true.
      return
      end
      logical function bsse_create_geom(rtdb)
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "geomP.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bsse_common.fh"       
c     
c     This logical function generate the geometries of monomers
c     with ghost atoms and store them into DB.
c
      integer rtdb              ![input]
      integer geom              ![input]
c
      character*255 name           
      character*16 tag_old  (mx_atm)
      character*16 tag_new  (mx_atm)
c
      integer mon_hnd
      integer k,j,l,m
c
      logical is_atm
      integer atn ![output]
      double precision mass(mx_atm)
      character*2 symbol
      character*16 element

c
      double precision q_old (mx_atm), q_new (mx_atm)
      double precision c(3,mx_atm)
      double precision ctmp(3,mx_atm)
      logical allbq
      integer lenb
c
      logical bsse_rtdb_store
      external bsse_rtdb_store
c
      lenname= inp_strlen(spr_name) 
c
      if (.not. geom_create(geom,spr_name(1:lenname))) 
     $   call errquit('bsse_create_geom: geom_create failed !',
     $                0,GEOM_ERR)
c
      if (.not. geom_rtdb_load(rtdb, geom, spr_name))
     $   call errquit('bsse_create_geom: no geometry load form rtdb',
     $                0,GEOM_ERR)
c      
      if (.not. geom_cart_get(geom, natoms, tag_old, c, q_old))
     $   call errquit('bsse_create_geom: get geom info fail !',
     $                0,GEOM_ERR)
c
c     -------------------------------------------
c     rename atoms tag to ghost and create geoms
c     -------------------------------------------
c
      do j = 1, nmon
        do l = 1, natoms
c
c         initialize all centers as bqX
c
           lenb=inp_strlen(tag_old(l))
           tag_new(l)= '  '
           tag_new(l)(1:lenb+2) = 'bq' // tag_old(l)
     T          (1:lenb)
           if(lenb.lt.14) tag_new(l)(lenb+3:lenb+3)=' '
           q_new(l) = 0.0d0
c
          do k = 1, mon_atm(j)
c
c           only do bsse with atoms
c
             is_atm = geom_tag_to_element(tag_old(l),symbol,
     E            element,atn)
            if ((.not. is_atm) .and. symbol.ne.'bq')
     $        call errquit('bsse_create_geom: not atom or bq',
     $                     0,GEOM_ERR)
c
c           compare and if this is center is an atom then
c           set its tag and charge to the proper values
c
            if (l.eq.(mon(j,k))) then
              tag_new(l) = tag_old(l)
              q_new(l)   = q_old(l)
              go to 100
            endif
          enddo
 100      continue 
        enddo
c
c       Sort out the masses of the various centers
c
        allbq=.true.
        do l = 1, natoms
           is_atm = geom_tag_to_element(tag_new(l),
     $          symbol,element,atn)
          if(symbol.ne.'bq') allbq=.false.
          if ((.not. is_atm) .and. symbol.ne.'bq')
     $       call errquit('bsse_create_geom: not atom or bq',
     $                    0,GEOM_ERR)
          if (.not.geom_atn_to_default_mass(atn,mass(l)))
     $       call errquit('bsse_create_geom: default mass failed',
     $                    911,GEOM_ERR)
       enddo

       if(allbq) then
          write(6,*) ' !!!! BSSE generated geometry made of all bqs'
          call errquit(' task_bsse: wrong geometry ',0,0)
       endif

c
        lenname= inp_strlen(mon_name(j)) 
c
        if (.not. geom_create(mon_hnd, mon_name(j)(1:lenname)//'g')) 
     $    call errquit('bsse_create_geom: geom_create failed !',
     $                 0,GEOM_ERR)
c
        if (.not. geom_cart_set(mon_hnd, natoms, tag_new, c, q_new))
     $    call errquit('bsse_create_geom: geom_cart_set failed',
     $                 0,GEOM_ERR)
c
        if (.not. geom_masses_set(mon_hnd, natoms, mass))
     $    call errquit('bsse_create_geom: geom_masses_set failed',
     $                 0,GEOM_ERR)
c:debug proposes
c     if (nodeid().eq. 0) then
c       if (.not. geom_print(mon_hnd))
c    $        call errquit('geom_input: print failed ', 0)
c     endif
c 
        if (.not.geom_rtdb_store(rtdb, mon_hnd,
     $                           mon_name(j)(1:lenname)//'g'))
     $    call errquit('bsse_create_geom: geom_rtdb_store failed',
     $                 60,GEOM_ERR)
c
        if (.not. geom_destroy(mon_hnd))
     $    call errquit('bsse_create_geom: geom_destroy failed',
     $                 60,GEOM_ERR)
c
c       creating and storing monomers without  ghost  
c
        if (.not. geom_create(mon_hnd, mon_name(j)(1:lenname)))
     $    call errquit('bsse_create_geom: geom_create failed!',0,
     $                 GEOM_ERR)
c
        m = 0
        do l = 1, natoms
          do k = 1, mon_atm(j)
            if (l.eq.(mon(j,k))) then
              m = m + 1
              tag_new(m) = tag_old(l)
              q_new  (m) = q_old  (l)
              ctmp(1,m) = c(1,l)
              ctmp(2,m) = c(2,l)
              ctmp(3,m) = c(3,l)
            endif
          enddo
        enddo
c      
        if (.not. geom_cart_set(mon_hnd,mon_atm(j),tag_new,ctmp,q_new))
     $    call errquit('bsse_create_geom: geom_cart_set failed',
     $                 0,GEOM_ERR)
c
c       Work out the masses of the atoms
c
        do l=1, mon_atm(j) 
          is_atm = geom_tag_to_element(tag_new(l),symbol,element,atn)
          if ((.not. is_atm) .and. symbol.ne.'bq')
     $      call errquit('bsse_create_geom: center is neither atom '//
     $                   'nor bq',0,GEOM_ERR)
c..       set default mass
          if(.not.geom_atn_to_default_mass(atn,mass(l)))
     &      call errquit('bsse_create_geom: default mass failed',
     &                   911,GEOM_ERR)
        enddo
c
        if (.not. geom_masses_set(mon_hnd, mon_atm(j), mass))
     $    call errquit('bsse_create_geom: geom_masses_set failed',
     $                 0,GEOM_ERR)
c
        if (nodeid() .eq. 0) then
          if(.not. geom_print(mon_hnd))
     $      call errquit('bsse_create_geom: print failed ', 0, GEOM_ERR)
        endif
c
        if (.not.geom_rtdb_store(rtdb, mon_hnd, mon_name(j)(1:lenname)))
     $    call errquit('bsse_create_geom: geom_rtdb_store failed',
     $                 0,GEOM_ERR)
c
        if (.not. geom_destroy(mon_hnd))
     $    call errquit('bsse_create_geom: geom_destroy failed',
     $                 0,GEOM_ERR)
c
c:debug proposes
c     if (.not. rtdb_print(rtdb, .true.)) call errquit('print failed',0)
c
c ======================================================
c  create basis
c     
c     status=bsse_create_basis(rtdb,geom, tag_new,natoms)
c========================================================
      enddo
c
      if(.not. geom_destroy(geom))
     $  call errquit('geom_input: geom_destroy failed', 0, GEOM_ERR)
c 
      bsse_create_geom = .true.
c
      return
      end
      logical function bsse_energy(rtdb,theory,final_spr_energy)
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "bsse_common.fh"
c
      integer rtdb              ![input]
c
      integer j ! run over monomers
      integer i 
      integer m_spr
c tmp
      character*(*) theory
      character*255 vec_spr
      character*255 vec_dbi, vec_dbo
c
      character*255 tmp
c
      double precision q_spr
c
      logical task_energy_doit
      external task_energy_doit
c
      logical bsse_rtdb_load
      external bsse_rtdb_load
c
      logical bsse_create_geom
      external bsse_create_geom
c
c
      bsse_energy=.false.
c
      if(ga_nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*)
        call util_print_centered(LuOut,
     $      'BSSE Energy Correction',
     $      40,.true.)
        write(LuOut,*)
        write(LuOut,*)
      endif
c
        if (.not. rtdb_get(rtdb,'bsse:first_j',mt_log,1,first_j))
     $  first_j=.true.
c
c
      if(.not.bsse_rtdb_load(rtdb))
     $  call errquit('bsse_energy: load data input in db', 911,RTDB_ERR)
c
c
      lentheo = inp_strlen(theory)
      lenname = inp_strlen(spr_name)
c
      if (.not. task_energy_doit(rtdb,theory,spr_energy))
     $  call errquit('bsse_energy: no geometry ',0,UNKNOWN_ERR)
c
c: take supermolecule total charge
      if(.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, q_spr))
     $  q_spr = 0.0d0


c: multiplicity
c:    density methods
        if ( theory(1:lentheo).eq.'dft' .or.
     $       theory(1:lentheo).eq.'tddft') then
          if (.not. rtdb_get(rtdb, 'dft:mult', mt_int, 1, m_spr))
     $      call errquit('bsse_energy: rtdb_get of mult failed',
     $                   0,RTDB_ERR )
c:   wavefuntion methods

        elseif( theory(1:lentheo).ne.'dft' .and.
     $          theory(1:lentheo).ne.'tddft') then
          if (.not. rtdb_get(rtdb,'scf:nopen', MT_INT, 1, m_spr))
     $      call errquit('bsse_energy: rtdb_get of nopen failed',
     $                   0,RTDB_ERR)

        endif

c
c: name of the original movecs
        if(theory(1:lentheo).eq.'dft' .or.
     $     theory(1:lentheo).eq.'tddft') then
          vec_dbo= 'dft:output vectors'
          vec_dbi= 'dft:input vectors'
        elseif (theory(1:lentheo).eq.'mcscf') then
          vec_dbo= 'mcscf:output vectors'
          vec_dbi= 'mcscf:input vectors'
        else
          vec_dbo = 'scf:output vectors'
          vec_dbi = 'scf:input vectors'
        endif
c
      if (.not. rtdb_cget(rtdb, vec_dbo, 1, vec_spr))
     $ call errquit('bsse_energy: get vectors file failed', 0,RTDB_ERR)
c
c: create geom for monomers within supermolecular geom


      if(.not.bsse_create_geom(rtdb))
     $  call errquit('bsse_energy: bsse_create_geom',911,UNKNOWN_ERR)


c
c:Obtain monomers energies from frozen geometries;
c:it makes a couple jobs for each monomer (no ghost, ghost)
c

      j = 1

      do i = 1, nmon*2 

        j_mon_name =  mon_name(j)
        lenname = inp_strlen(mon_name(j))
c

        if (mod(i,2).eq.0) then

          j_mon_name = j_mon_name(1:lenname)//'g'
          lenname = lenname + 1

        endif
c      
        call bsse_param(rtdb, mmon(j), qmon(j), j_mon_name, input(i),
     $                  theory)
c
c:evaluate energy
        if (.not. task_energy_doit(rtdb,theory, mon_energy(i)))
     $    call 
     $ errquit('bsse_energy: failed calling task_energy',0,UNKNOWN_ERR)

c
        if (mod(i,2).eq.0) then
          j = j +1
        endif

      enddo


c-----------------------------------------------------------------------
c      Evaluate BSSE error for supermolecular geometry
c-----------------------------------------------------------------------


        bsse_error = 0.0d0

        do j = 1, nmon,2

          m_error(j) = mon_energy(j) - mon_energy(j+1)
          bsse_error = bsse_error + m_error(j)

       enddo

      final_spr_energy = spr_energy + bsse_error    
c
c:return to original active geom


      lenname = inp_strlen(spr_name)


      if (.not. rtdb_cput(rtdb,'geometry',1,spr_name(1:lenname)))
     $  call errquit('bsse_energy: no geometry ',0, RTDB_ERR)

c:return to original output vectors
c
      if (.not. rtdb_cput(rtdb, vec_dbi, 1, vec_spr))
     &  call errquit('bsse_energy: input_vectors',0, RTDB_ERR)
c

      if (.not. rtdb_cput(rtdb, vec_dbo, 1,vec_spr))
     &  call errquit('bsse_energy: output_vectors',0, RTDB_ERR)


c
c:return to original  charge


      if (.not. rtdb_put(rtdb, 'charge', MT_DBL, 1, q_spr))
     $  call errquit
     $  ('bsse_energy: failed to write charge to rtdb', 0, RTDB_ERR)


c: multiplicity

c:    density methods

        if ( theory(1:lentheo).eq.'dft' .or.
     $       theory(1:lentheo).eq.'tddft') then

          if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, m_spr))
     $      call errquit('bsse_energy: rtdb_put of mult failed',
     $                   0, RTDB_ERR)

c:   wavefuntion methods

        elseif( theory(1:lentheo).ne.'dft' .and.
     $          theory(1:lentheo).ne.'tddft') then

          if (.not. rtdb_put(rtdb,'scf:nopen', MT_INT, 1, m_spr))
     $      call errquit('bsse_energy: rtdb_put of nopen failed',
     $                   0, RTDB_ERR)


        endif
c
c:put into db final energy associated with theory


      tmp = theory(1:lentheo)//':energy'
c
      if (.not. rtdb_put(rtdb,tmp,
     $                                      MT_DBL,1,final_spr_energy))
     $  call 
     $ errquit('bsse_energy: failed to write charge to rtdb',0,RTDB_ERR)


c:debug proposes vama
c     if (.not. rtdb_print(rtdb, .true.)) call errquit('print failed',0)
c
c     if (nodeid().eq.0) then
c       do j = 1, nmon
c         write(LuOut,10) j, m_error(j)
c 10      format
c    $      (/' error  for monomer ', i4, ' =', f20.12/)
c       enddo
c

      if(ga_nodeid().eq.0) then

        write(LuOut,20) bsse_error, spr_energy, final_spr_energy
  20  format (/'             BSSE error = ', f20.12/
     $         '  Supermolecular energy = ', f20.12/
     $         '       Corrected energy = ', f20.12/)

        call util_flush(luout)

      endif
c

        if (.not. rtdb_put(rtdb,'bsse:first_j',mt_log,1,.false.))
     $    call errquit('bsse_energy: rtdb_put failed',0,RTDB_ERR)

c
      call ga_sync
c
      bsse_energy = .true.
      end
      logical function bsse_gradient(rtdb,theory,final_spr_energy, 
     $                                gx_spr)
c:debug
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "inp.fh"
#include "stdio.fh"
#include "bsse_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "tcgmsg.fh"
c
c
      integer rtdb              ![input]
c
      integer j,i,l,k,n
      integer geom
c
      character*255 vec_dbi, vec_dbo
      character*255 vec_spr
      character*(*) theory
      character*16 tag
      character*255 tmp
c
      double precision q_spr
      double precision gx_spr(3,*)
      double precision gx_mon(3,mx_atm)
      double precision crd(3), q
c
      logical task_gradient_doit
      external task_gradient_doit
      logical bsse_rtdb_load
      external bsse_rtdb_load
      logical bsse_create_geom
      external bsse_create_geom
c
c
      bsse_gradient=.false.
c
      if (nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*)
        call util_print_centered(LuOut,
     $      'BSSE Correction to Energy Gradient',
     $      40,.true.)
        write(LuOut,*)
        write(LuOut,*)
      endif
c
        if (.not. rtdb_get(rtdb,'bsse:firt_j',MT_LOG,1,first_j))
     $  first_j=.true.
c
      if(.not.bsse_rtdb_load(rtdb))
     $  call 
     $ errquit('bsse_gradient: load data input in db',911,UNKNOWN_ERR)
c

      lentheo = inp_strlen(theory)
      lenname = inp_strlen(spr_name)
c

      n_cart= 3*natoms
c
c: get supermolecular gradient and energy


      if (.not.task_gradient_doit(rtdb,theory,spr_energy, gx_spr))
     $   call 
     $ errquit('bsse_gradient: no task gradient do it ',0,UNKNOWN_ERR)


c: get supermolecule charge

      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, q_spr))
     $   q_spr = 0.0d0

      
c
c: multiplicity
c
c:    wavefuntion methods
        if( theory(1:lentheo).ne.'dft') then
          if (.not. rtdb_get(rtdb,'scf:nopen', MT_INT, 1, m_spr))
     $      call 
     $  errquit('bsse_gradient: rtdb_get of nopen failed',0,RTDB_ERR)
c
c:    density methods
        elseif ( theory(1:lentheo).eq.'dft') then
          if (.not. rtdb_get(rtdb, 'dft:mult', mt_int, 1, m_spr))
     $      call 
     $  errquit('bsse_gradient: rtdb_get failed', 0,RTDB_ERR)
        endif
c
c: name of the original movecs

        if(theory(1:lentheo).eq.'dft') then
          vec_dbo= 'dft:output vectors'
          vec_dbi= 'dft:input vectors'
        elseif (theory(1:lentheo).eq.'mcscf') then
          vec_dbo= 'mcscf:output vectors'
          vec_dbi= 'mcscf:input vectors'
        else
          vec_dbo = 'scf:output vectors'
          vec_dbi = 'scf:input vectors'
        endif
c
      if (.not. rtdb_cget(rtdb, vec_dbo, 1, vec_spr))
     $ call 
     $ errquit('bsse_gradient: get vectors file failed', 0,RTDB_ERR)
c
c
c: create geom for monomers within supermolecular geom


      if(.not.bsse_create_geom(rtdb))
     $  call errquit('bsse_gradient: bsse_create_geom',60,UNKNOWN_ERR)

c
c: Obtain monomers energies from forzen geometries
c: it makes a couple job for each monomer (no ghost, ghost)
      j = 1
      do i = 1, nmon*2 
        j_mon_name =  mon_name(j)
        lenname = inp_strlen(mon_name(j))

c
        if (mod(i,2).eq.0) then
          j_mon_name = j_mon_name(1:lenname)//'g'
          lenname = lenname + 1 
        endif
c      
        call bsse_param(rtdb, mmon(j), qmon(j), j_mon_name, input(i),
     $                  theory)
c
c: evaluate gradient

      if (.not.task_gradient_doit(rtdb,theory,mon_energy(i),gx_mon))
     $  call 
     $  errquit('bsse_gradient: error calling task gradient spr',
     $   0,UNKNOWN_ERR)
c
c: mix monomers grads with supermolecule grad
        if (mod(i,2).eq.0) then

          do k = 1, natoms

            do n=1 ,3
              gx_spr(n,k) = gx_spr(n,k) - gx_mon(n,k)
            enddo

          enddo
c: go for next monomer
          j = j +1
c
        else
c
          do k = 1 , natoms
            do l = 1, mon_atm(j)
              if (k.eq.mon(j,l)) then
                do n = 1,  3
                 gx_spr(n,k)=gx_spr(n,k) + gx_mon(n,l)
                enddo
              endif
            enddo
          enddo
        endif
c
      enddo

c-----------------------------------------------------------------------      
c      Evaluate supermolecular energy free of BSSE
c-----------------------------------------------------------------------      

      bsse_error = 0.0d0

      do j = 1, nmon,2
        bsse_error = bsse_error+(mon_energy(j) - mon_energy(j+1))
      enddo
c

      final_spr_energy = spr_energy + bsse_error
c

      lenname = inp_strlen(spr_name)
c
c
      if (.not. rtdb_cput(rtdb,'geometry', 1, spr_name(1:lenname)))
     $  call errquit('bsse_gradient: no geometry ',0,RTDB_ERR)


c:return to original  charge 

      if (.not. rtdb_put(rtdb, 'charge', MT_DBL, 1, q_spr))
     $  call errquit
     $  ('bsse_gradientt: failed to write charge to rtdb', 0,RTDB_ERR)


c:return to original output vectors
c
      if (.not. rtdb_cput(rtdb, vec_dbi, 1, vec_spr))
     &  call errquit('bsse_gradient: store input_vectors',60,RTDB_ERR)
c

      if (.not. rtdb_cput(rtdb, vec_dbo, 1,vec_spr))
     &  call errquit('bsse_gradient: store output_vectors',60,RTDB_ERR)



c:   wavefuntion methods

        if( theory(1:lentheo).ne.'dft') then
          if (.not. rtdb_put(rtdb,'scf:nopen', MT_INT, 1, m_spr))
     $      call 
     $  errquit('bsse_gradient: rtdb_put of nopen failed',0,RTDB_ERR)

c:    density methods

        elseif ( theory(1:lentheo).eq.'dft') then
          if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, m_spr))
     $      call errquit('bsse_gradient: rtdb_put failed', 0,RTDB_ERR)
        endif
c
c:put into db final gradien associated with theory


      tmp=theory(1:lentheo)//':gradient'     
      if (.not. rtdb_put(rtdb,tmp  ,mt_dbl,
     $    3*natoms, gx_spr)) call errquit
     $     ('bsse_gradient: could not store gradients',1, RTDB_ERR)
c
c:put into db final energy associated with theory
c
      tmp=theory(1:lentheo)//':energy'     
c
      if (.not. rtdb_put(rtdb,tmp,
     $                             MT_DBL,1,final_spr_energy))
     $  call 
     $  errquit('bsse_gradient: failed to write charge to rtdb'
     $   ,0,RTDB_ERR)

c
        if (.not. geom_create(geom, spr_name(1:lenname))) call errquit
     $   ('bsse_gradient: geom_create final failed !', 0,GEOM_ERR)
c
c
        if (.not.geom_rtdb_load(rtdb, geom, spr_name(1:lenname)))
     $   call errquit('bsse_gradient: load geom failed !',0, GEOM_ERR)

c:print final energy gradient

      if(nodeid().eq.0) then
c 
        write(luout,1000) theory(1:inp_strlen(theory)),
     $        'x','y','z','x','y','z'

        do  i=1, natoms
          if (.not. geom_cent_get(geom, i, tag, crd, q)) call errquit
     $           ('bsse_gradient: geometry corrupt?',0, GEOM_ERR)
          write(luout,2000) i, tag,(crd(j),j=1,3),
     $           (gx_spr(j,i),j=1,3)
        enddo

        write(luout,*)
        write(luout,*)

 1000   format(/,/,25X,A,' BSSE ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $        'coordinates',
     $        24X,'gradient',/,6X,2(1X,(3(10X,A1))))

 2000   format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
c
c
        write(LuOut,20) bsse_error, spr_energy, final_spr_energy

  20    format (/'             BSSE error = ', f20.12/
     $           '  Supermolecular energy = ', f20.12/
     $           '       Corrected energy = ', f20.12/)
        call util_flush(luout)

      endif


        if(.not. geom_destroy(geom))
     $   call errquit('bsse_gradient: geom_destroy failed', 0,GEOM_ERR)
c
c      

        if (.not. rtdb_put(rtdb,'bsse:first_j',mt_log,1,.false.))
     $    call errquit('bsse_gradient: rtdb_put failed',0,RTDB_ERR)
c
      call ga_sync
c

      bsse_gradient = .true.
      end
      logical function bsse_hessian(rtdb)
      implicit none
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "bsse_common.fh"
#include "stdio.fh"
#include "tcgmsg.fh"
#include "errquit.fh"
c
      integer rtdb              ![input]
c
c     Call task_hessian_doit for each monomer
c     and supermolecule and store each one hessian matrix
c     in different files. As out put write .hess file with
c     hessian taking account the BSSE
c
      logical status
      integer j, i
c
      character*255 vec_dbi, vec_dbo
      character*255 vec_spr
      character*32  theory
c
      character*16 mult(8)
      character*(nw_max_path_len) filehess
      character*(nw_max_path_len) filehesstmp
      integer m_spr
      double precision q_spr
      integer lenhess
      integer nopen
c
      logical  task_hessian_doit
      external task_hessian_doit
c
      logical bsse_rtdb_load
      external bsse_rtdb_load
c
      logical bsse_create_geom
      external bsse_create_geom
c
c     double precision final_spr_energy
c
      data mult/'singlet','doublet','triplet','quartet',
     $ 'quintet','sextet','septet','octet'/
c
c     call ecce_print_module_entry('task hessian')
c
      bsse_hessian= .false.
c
      if (nodeid().eq.0) then
        write(LuOut,*)
        write(LuOut,*)
        call util_print_centered(LuOut,
     $        'BSSE Hessian Correction',
     $         40,.true.)
        write(LuOut,*)
        write(LuOut,*)
      endif
c

        if (.not. rtdb_get(rtdb,'bsse:firt_j',mt_log,1,first_j))
     $  first_j=.true.
c

      if(.not.bsse_rtdb_load(rtdb))
     $  call errquit('bsse_hessian: load data input in db',
     $    911,UNKNOWN_ERR )
c
c

        if (.not. rtdb_cget(rtdb,'task:theory', 1, theory))
     &    call errquit('bsse_hessian: get input_vectors',0,RTDB_ERR)


      lentheo = inp_strlen(theory)
      lenname = inp_strlen(spr_name)
c
c

      if(.not.task_hessian_doit(rtdb))
     $  call errquit('bsse_hessian: error calling hessian doit failed',
     $      0,UNKNOWN_ERR)
c
      if (.not.rtdb_get
     $  (rtdb,'task:energy',mt_dbl,1,spr_energy))
     $  call errquit('bsse_hessian: cannot get energy ',0,RTDB_ERR)


c  store original sets
c: take supermolecule total charge

      if(.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, q_spr))
     $  q_spr = 0.0d0

c: multiplicity
c:   wavefuntion methods

      if( theory(1:lentheo).ne.'dft') then

        if (.not. rtdb_get(rtdb,'scf:nopen', MT_INT, 1, m_spr))
     $    call errquit('bsse_hessian: rtdb_put of nopen failed',
     $      0,RTDB_ERR)

c:    density methods

      elseif ( theory(1:lentheo).eq.'dft') then

        if (.not. rtdb_get(rtdb, 'dft:mult', mt_int, 1, m_spr))
     $   call errquit('bsse_hessian: rtdb_put mult failed', 0, RTDB_ERR)

      endif
c
c: name of the original movecs
        if(theory(1:lentheo).eq.'dft') then
          vec_dbo= 'dft:output vectors'
          vec_dbi= 'dft:input vectors'
        elseif (theory(1:lentheo).eq.'mcscf') then
          vec_dbo= 'mcscf:output vectors'
          vec_dbi= 'mcscf:input vectors'
        else
          vec_dbo = 'scf:output vectors'
          vec_dbi = 'scf:input vectors'
        endif
c
      if (.not. rtdb_cget(rtdb, vec_dbo, 1, vec_spr))
     $ call util_file_name('movecs',.false.,.false.,vec_spr)
c    $ call errquit('bsse_hessian: get vectors file failed', 0)
c
c: name of hessian file

      call util_file_name('hess',  .false., .false.,filehess)
c:copy to operate the sum at the end|
c
      if (nodeid().eq.0) then

        lenhess= inp_strlen(filehess)

        filehesstmp = ' '
        filehesstmp= filehess(1:lenhess)//'.spr'

        call util_file_copy(filehess, filehesstmp)

      endif

c: create geom for monomers within supermolecular geom

      if(.not.bsse_create_geom(rtdb))
     $  call errquit('bsse_hessian: bsse_create_geom',911,UNKNOWN_ERR)

c
c:  Obtain monomers energies from forzen geometries;
c:  it makes a couple jobs for each monomer (no ghost, ghost)
c:monomers

      j=1

      do i=1, nmon*2   

        j_mon_name =  mon_name(j)
        lenname = inp_strlen(mon_name(j))
c
        if (mod(i,2).eq.0) then

          j_mon_name = j_mon_name(1:lenname)//'g'
          lenname = lenname + 1

        endif
c

        if (.not. rtdb_cput(rtdb,'geometry', 1, j_mon_name(1:lenname)))
     $    call errquit('bsse_hessian: no geometry ',0,RTDB_ERR)
c
c

        call bsse_param(rtdb, mmon(j), qmon(j), j_mon_name, input(i),
     $                  theory)
c
c
c
        if(.not.task_hessian_doit(rtdb))
     $    call errquit('bsse_hessian: error calling hessiandoit',
     $ 0,UNKNOWN_ERR)

c:get evaluated energy

        if (.not.rtdb_get
     $    (rtdb,'task:energy', mt_dbl, 1, mon_energy(i)))
     $     call errquit('bsse_hessian: get monomer energy',0,RTDB_ERR)

c:copy to operate the sum at the end|
c
        if(nodeid().eq.0) then
c

          filehesstmp= filehess(1:lenhess)//'.'//
     &      j_mon_name(1:lenname)

          call util_file_copy(filehess, filehesstmp)

        endif
c
        if (mod(i,2).eq.0) then
          j = j +1
        endif

      enddo
c
      bsse_error = 0.0d0

      do j = 1, nmon,2

        m_error(j) = mon_energy(j) - mon_energy(j+1)
        bsse_error = bsse_error + m_error(j)

      enddo
c
      final_spr_energy = spr_energy + bsse_error
c
c sum matrix

      if(ga_nodeid().eq.0) then

        call  hessian_matrix_sum (rtdb)

      endif
c
c: return original geom and parameters
c: store supermolecule parameters
c
      if (.not. rtdb_cput(rtdb, vec_dbi, 1, vec_spr))
     &  call errquit('bsse_hessian: put input_vectors',0,RTDB_ERR)
c

      if (.not. rtdb_cput(rtdb, vec_dbo, 1,vec_spr))
     &  call errquit('bsse_hessian: output_vectors',0,RTDB_ERR)



c:return to original active geom

c
      lenname = inp_strlen(spr_name)

      if (.not. rtdb_cput(rtdb,'geometry', 1, spr_name(1:lenname)))
     $   call errquit ('bsse_hessian: no geometry ',0,RTDB_ERR)

c: store supermolecule total charge

      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, q_spr))
     $   q_spr = 0.0d0

c: store multiplicity
c:   wavefuntion methods

      if( theory(1:lentheo).ne.'dft') then

        if (.not. rtdb_put(rtdb,'scf:nopen', MT_INT, 1, m_spr))
     $    call errquit('bsse_hessian: rtdb_put of nopen failed'
     $  , 0, RTDB_ERR)

c:    density methods

      elseif ( theory(1:lentheo).eq.'dft') then

        if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, m_spr))
     $    call errquit('bsse_gradient:cant put m_spr rtdb',0, RTDB_ERR)

      endif

c

      if(ga_nodeid().eq.0) then

        write(LuOut,20) bsse_error, spr_energy, final_spr_energy
  20  format (/'             BSSE error = ', f20.12/
     $         '  Supermolecular energy = ', f20.12/
     $         '       Corrected energy = ', f20.12/)

         call util_flush(luout)

      endif


      if (.not. rtdb_put(rtdb,'bsse:first_j',mt_log,1,.false.))
     $  call errquit('bsse_hessian: rtdb_put failed',0, RTDB_ERR)

c


      bsse_hessian = .true.

c
      call ga_sync
c
      end
      subroutine hessian_matrix_sum(rtdb)
c
c  This function makes the apropiate sum to mix the hessian of monomers matrixes
C  (with ghost and without ghost) to supermolecular hessian matrix
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
#include "bsse_common.fh"
#include "errquit.fh"
c
      integer rtdb
      integer j,i
      integer iii
      integer k
      integer nat2
      integer lenhess
      double precision dval

      integer nhesst
      character*(nw_max_path_len) filehess
c
      integer k_exy_mont
      integer k_exy_mon
      integer k_exy_spr
      integer k_exycp
      integer k_exybt
      integer k_exybs
c
      integer l_exycp
      integer l_exybt
      integer l_exy_mon
      integer l_exy_mont
      integer l_exybs
      integer n3xyz
c
      n3xyz=3*natoms
      nhesst=n3xyz*(n3xyz+1)/2
      nat2=n3xyz*n3xyz
c
      if (.not.ma_push_get(mt_dbl,nat2,'bsse exybs ',
     *  l_exybs,k_exybs))
     *  call errquit('bsse_hessian_exy: cannot allocate',650, MA_ERR)
      call dfill(nat2,0.0d00,dbl_mb(k_exybs),1)
c
      if (.not.ma_push_get(mt_dbl,nhesst,'bsse exybt ',
     *  l_exybt,k_exybt))
     *  call errquit('bsse_hessian_exy: cannot allocate',650, MA_ERR)
      call dfill(nhesst,0.0d00,dbl_mb(k_exybt),1)
c
      if (.not.ma_push_get(mt_dbl,nat2,'bsse exycp ',
     *  l_exycp,k_exycp))
     *  call errquit('bsse_hessian_exy: cannot allocate',650, MA_ERR)
      call dfill(nat2,0.0d00,dbl_mb(k_exycp),1)
c
      if (.not.ma_push_get(mt_dbl,nat2,'bsse:hessian:exy mon',
     *  l_exy_mon,k_exy_mon))
     *  call errquit('bsse_hessian_exy: cannot allocate',650, MA_ERR)
c
      if (.not.ma_push_get(MT_DBL,nhesst,'bsse:hessian:exyt mon',
     *  l_exy_mont,k_exy_mont))
     *  call errquit('bsse_hessian_exy: cannot allocate',650, MA_ERR)
c
c:read triangle matrix super and do square


      call util_file_name('hess', .false., .false.,filehess)
c

      lenhess=inp_strlen(filehess)
      filehess=filehess(1:lenhess)//'.spr'
c

      open(unit=69,file=filehess,form='formatted',status='old',
     &    err=99900,access='sequential')
      do iii = 0,(nhesst-1)
        read(69,*,err=99901,end=99902) dval
        dbl_mb(k_exybt+iii) = dval
      enddo
      close(unit=69,status='keep')
c

      call trin2squa(dbl_mb(k_exybs),dbl_mb(k_exybt),n3xyz)
c

      j=1

      do i=1, nmon*2

        j_mon_name =  mon_name(j)
        lenname = inp_strlen(mon_name(j))

        if (mod(i,2).eq.0) then

          j_mon_name = j_mon_name(1:lenname)//'g'

          lenname = lenname + 1
          filehess=filehess(1:lenhess)//'.'//j_mon_name

          n3xyz=3*natoms

          nhesst = n3xyz*(n3xyz+1)/2
c:read triangle matrix monomer and do square
          open(unit=69,file=filehess,form='formatted',status='old',
     &    err=99900,access='sequential')

          do iii = 0,(nhesst - 1)
            read(69,*,err=99901,end=99902) dval

            dbl_mb(k_exy_mont+iii) = dval

          enddo

          close(unit=69,status='keep')
    
          call trin2squa(dbl_mb(k_exy_mon),dbl_mb(k_exy_mont),n3xyz)
c:substract entire ghost matrix case
          do k = 0,(nat2- 1)

            dval=dbl_mb(k_exybs+k)-dbl_mb(k_exy_mon+k)
            dbl_mb(k_exybs+k)=dval

          enddo
c:visit other monomer

          j = j +1

        else

          n3xyz=3*mon_atm(j)

          nhesst = n3xyz*(n3xyz+1)/2
c:read triangle matrix monomer 

          filehess=filehess(1:lenhess)//'.'//j_mon_name
c

          open(unit=69,file=filehess,form='formatted',status='old',

     &    err=99900,access='sequential')

          do iii = 0,(nhesst - 1)

            read(69,*,err=99901,end=99902) dval
            dbl_mb(k_exy_mont+iii) = dval

          enddo

          close(unit=69,status='keep')

c: do square
          call trin2squa(dbl_mb(k_exy_mon),dbl_mb(k_exy_mont),n3xyz)

c:sum of no-ghost monomer case
          call dcopy(nat2,dbl_mb(k_exybs),1,dbl_mb(k_exycp),1)

c
          call sum_small_matrix(dbl_mb(k_exybs),dbl_mb(k_exycp),
     $                             dbl_mb(k_exy_mon),3*natoms, n3xyz, j)
        endif
      enddo
c
c     n3xyz=3*natoms
c
c:write the final hessian
c
      call util_file_name('hess', .false., .false.,filehess)
c
      call  stpr_wrt_fd_from_sq(dbl_mb(k_exybs),n3xyz,filehess)

c:clean memory 

      if (.not.ma_chop_stack(l_exy_mont))
     *   call errquit('hess_read: cannot deallocate Exyt',555, MA_ERR)

      if (.not.ma_chop_stack(l_exy_mon))
     *   call errquit('hess_read: cannot deallocate Exyt',555, MA_ERR)

      if (.not.ma_chop_stack(l_exycp))
     *  call errquit('hess_read: cannot deallocate l_Exycp',555, MA_ERR)

      if (.not.ma_chop_stack(l_exybt))
     *  call errquit('hess_read: cannot deallocate l_Exybt',555, MA_ERR)

      if (.not.ma_chop_stack(l_exybs))
     *  call errquit('hess_read: cannot deallocate l_Exybs',555, MA_ERR)
c
      return
99900 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('error opening file: "filehess"',
     1             DISK_ERR,911)
99901 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('error reading file: "filehess"',
     1             DISK_ERR,911)
99902 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('unexpected EOF reading file: "filehess"',
     1             DISK_ERR,911)
      end
c
      subroutine sum_small_matrix(d,a,b,dim_a,dim_b,j)
c
c     this routine sums a monomer matrix in fractions 3x3  of each atom 
c     to the supermolecular matrix
c vama
#include "bsse_common.fh"
      integer i,k,m,n,l,g,r,c
      integer dim_a  ! dimension matrix super
      integer dim_b  ! dimension matrix monomer
      integer j      ! monomer j, used to take the list of atoms
      double precision d(dim_a,dim_a) ! matrix super resul
      double precision a(dim_a,dim_a) ! matrix super
      double precision b(dim_b,dim_b) ! matrix monomer
      double precision dval 
c
      do i = 1 , natoms
        do k = 1 , natoms
c
          do m= 1, mon_atm(j)
            do n= 1, mon_atm(j)
c
              l= mon(j,n)
              g= mon(j,m)
c
              if(l.eq.k.and.g.eq.i) then
c
                do r=0, 2
                  do c=0,2
                    dval=a((i-1)*3+1+r,(k-1)*3+1 +c)+ 
     $                                        b((m-1)*3+1+r,(n-1)*3+1+c)

                    d((i-1)*3+1+r,(k-1)*3+1+c)=dval 
                    d((k-1)*3+1+c,(i-1)*3+1+r)=dval
                  enddo
                enddo
              endif   
            enddo
          enddo
        enddo
      enddo
      end
c
      subroutine trin2squa(a,b,dime)
c
c     convert a triangular matrix to a simetric square matrix
c
      integer dime
      integer j, l, i
      double precision a(dime,dime) !squ
      double precision b((dime*(dime+1))/2) !tri
      double precision dval
c
      i = 1
      do j=1, dime
        do l= 1, j
          dval=b(i)
          a(j,l)=dval
          a(l,j)=dval
          i=i+1
        enddo
      enddo
      end
      subroutine cfill(n,val,a,ia)
      implicit none
      integer n, ia
      character*(*) val, a(*)
      integer i
c
c     initialise characters array
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end

c $Id$
