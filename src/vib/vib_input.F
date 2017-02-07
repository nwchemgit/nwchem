      subroutine freq_vib_input(rtdb)
      implicit none
#include "errquit.fh"
*
* $Id$
*
c
c  {freq||vib}
c    reuse [<string> hessian_filename]
c    mass <integer> lexical_index <real> new_mass
c    mass <string> tag_identifier <real> new_mass
c    animate [<real> step_size_for_animation]     
c    fd_delta [<real> step_size_for_fd_hessian]     
c    filename <string> file_set_name [overwrite]
c  end
c  
#include "stdio.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
*::passed
      integer rtdb 
*::local      
      logical o_reuse, o_hessfile, o_masses
      logical reported
      character*255 test
      character*32 rtdb_id
      character*255 filename
      character*50 zone_name
      character*16 tag
      integer itest,ind
      integer lexical_index
      integer old_num_mod_mass
      integer num_mod_mass
      double precision new_mass, animate_step,fd_step
      integer len, num_on_line
      integer numtemps, l_temp, k_temp, itemp
      double precision valuetemp,kv(3)


      num_mod_mass = 0
      old_num_mod_mass = 0
      reported = .false.
      o_reuse = .false.
      o_hessfile = .false.
      o_masses = .false.
      call inp_set_field(0) ! reset to beginning of current line
      if (.not.inp_a(test)) call errquit
     &    ('freq_vib_input: could not parse keyword ',911, INPUT_ERR)
      if (.not.((inp_compare(.false.,test(1:3),'vib')).or.
     &          (inp_compare(.false.,test(1:4),'freq'))))
     &   call errquit
     &   ('freq_vib_input: input block identifier not "vib" or "freq" ',
     &   911, INPUT_ERR)
      
00010 if (.not.inp_read()) goto 01001
      if (inp_a(test)) then
        if (inp_compare(.false.,test,'end')) then
          goto 01000
        else if (inp_compare(.false.,test,'animate')) then
          if (.not.rtdb_put(rtdb,'vib:animate',mt_log,1,.true.))
     &        call errquit('freq_vib_input: rtdb_put failed',911,
     &       RTDB_ERR)
          num_on_line = inp_n_field()
          if (num_on_line .eq. 1) then
            continue ! okay just set animation on
          else if (num_on_line.eq.2) then
            if (inp_f(animate_step)) then
              if (.not.rtdb_put(rtdb,
     &            'vib:animate:step_size',mt_dbl,1,animate_step))
     &            call errquit('freq_vib_input: rtdb_put failed',911,
     &       RTDB_ERR)
            else
              write(luout,*)' could not read step size '
              call inp_outrec()
              goto 01001
            endif
          else
            write(luout,*)' malformed animate line '
            write(luout,*)
     &          ' should be 1 or 2 elements on line there are ',
     &          num_on_line
            call inp_outrec()
            goto 01001
          endif
          goto 00010
        else if (inp_compare(.false.,test,'fd_delta')) then
            if (inp_f(fd_step)) then
              if (.not.rtdb_put(rtdb,
     &            'stpr_gen:delta',mt_dbl,1,fd_step))
     &            call errquit('freq_vib_input: rtdb_put failed',911,
     &       RTDB_ERR)
            else
              write(luout,*)' could not read fd_delta step size '
              call inp_outrec()
              goto 01001
            endif
          goto 00010
        else if (inp_compare(.false.,test,'reuse')) then
          o_reuse = .true.
          if (inp_a(test)) then
            len = inp_strlen(test)
            filename = test(1:len)
            o_hessfile = .true.
          else
            o_hessfile = .false.
          endif
          goto 00010


        else if((inp_compare(.false.,test,'zone_name')).or.
     &          (inp_compare(.false.,test,'zone_structure_name'))) then
           if (.not. inp_a(zone_name))
     >     call errquit(
     >          'vib_input: failed to read zone_name',0,INPUT_ERR)
            ind = index(zone_name,' ') - 1
            if (.not.rtdb_cput(rtdb,'vib:zone_name',1,zone_name(1:ind)))
     >      call errquit(
     >       'vib_input: error writing to rtdb',1806,RTDB_ERR)
          goto 00010

        else if ((inp_compare(.false.,test,'dispersion')).or.
     &           (inp_compare(.false.,test,'brillouin_zone'))) then
           call vib_zone_input(rtdb)
          goto 00010

        else if ((inp_compare(.false.,test,'temperature')).or.
     &           (inp_compare(.false.,test,'temp'))) then
          if (.not. inp_i(numtemps)) 
     &      call errquit(
     &         'freq_vib_input: Forgot the number of temperatures',
     &         911, INPUT_ERR)
          if (numtemps.gt.0) then
            if (.not.MA_push_get(MT_dbl, numtemps, 'Temps',l_temp,
     &           k_temp))
     &        call errquit('freq_vib_input: failed to allocate temps',
     &          911, INPUT_ERR)
            do itemp = 1, numtemps
              if (.not. inp_f(valuetemp))
     &          call errquit(
     &            'freq_vib_input: failed to read a temperature',
     &            911, INPUT_ERR)
              dbl_mb(k_temp - 1 + itemp) = valuetemp
            enddo
            if (.not.rtdb_put(rtdb,'vib:temperature',mt_dbl,
     &        numtemps, dbl_mb(k_temp)))
     &        call errquit(
     &          'freq_vib_input: failed to write temps to rtdb',
     &          911, RTDB_ERR)
            if (.not.ma_pop_stack(l_temp))
     &        call errquit(
     &          'freq_vib_input: failed to release l_temp', 911,
     &          MA_ERR)
          else
            call errquit('freq_vib_input: Sorry, the nTemp <= 0',
     &        911, INPUT_ERR)
          endif
          goto 00010

        else if (inp_compare(.false.,test,'mass')) then
          if (.not.rtdb_get
     &        (rtdb,'vib:remass:count',mt_int,1,old_num_mod_mass))
     &        old_num_mod_mass = 0
          num_on_line = inp_n_field()
          if (num_on_line .ne. 3) then
            write(luout,*)'freq_vib_input: Invalid mass line '
            call inp_outrec()
            goto 01001
          endif
          if (inp_i(itest)) then
            lexical_index = itest
            if (inp_f(new_mass)) then
              num_mod_mass = num_mod_mass + 1
              write(rtdb_id,10001)num_mod_mass
              if (.not.rtdb_put(rtdb,rtdb_id,mt_int,1,lexical_index))
     &            call errquit('freq_vib_input: rtdb_put failed',911,
     &       RTDB_ERR)
              write(rtdb_id,10003)num_mod_mass
              if (.not.rtdb_put(rtdb,rtdb_id,mt_dbl,1,new_mass))
     &            call errquit('freq_vib_input: rtdb_put failed',911,
     &       RTDB_ERR)
            else
              call errquit('freq_vib_input: could not read mass',911,
     &       INPUT_ERR)
            endif
          else
            if (inp_a(test)) then
              tag = ' '
              len = inp_strlen(test)
              tag = test(1:len)
              if (inp_f(new_mass)) then
                num_mod_mass = num_mod_mass + 1
                write(rtdb_id,10002)num_mod_mass
                if (.not.rtdb_cput(rtdb,rtdb_id,1,tag))
     &              call errquit('freq_vib_input: rtdb_put failed',911,
     &       RTDB_ERR)
                write(rtdb_id,10003)num_mod_mass
                if (.not.rtdb_put(rtdb,rtdb_id,mt_dbl,1,new_mass))
     &              call errquit('freq_vib_input: rtdb_put failed',911,
     &       RTDB_ERR)
              else
                call errquit('freq_vib_input: could not read mass',911,
     &       INPUT_ERR)
              endif
            else
              write(luout,*)'freq_vib_input: could not parse line'
              call inp_outrec()
              goto 01001
            endif
          endif
          if ((old_num_mod_mass.gt.0).and.
     &        (.not.(reported)).and.
     &        (.not.(o_masses))) then
            write(luout,*)
     &          ' WARNING: old mass redefinitions being replaced'
            reported = .true.
          endif
          o_masses = .true.
          goto 00010
        else if (inp_compare(.false.,test,'filename')) then
          if (.not.rtdb_put(rtdb,'vib:animate:overwrite_files',mt_log,1,
     &                      .false.))
     &        call errquit('freq_vib_input: rtdb_put failed',911,
     &                     RTDB_ERR)
          if (inp_a(test)) then
            if (inp_compare(.false.,test,'overwrite')) then
              if (.not.rtdb_put(rtdb,'vib:animate:overwrite_files',
     &                          mt_log,1,.true.))
     &              call errquit('freq_vib_input: rtdb_put failed',911,
     &                           RTDB_ERR)
            else
              if (.not.rtdb_cput(rtdb,'vib:animate:filename',1,
     &                           test(1:inp_strlen(test))))
     &              call errquit('freq_vib_input: rtdb_cput failed',911,
     &                           RTDB_ERR)
              if (inp_a(test)) then
                if (inp_compare(.false.,test,'overwrite')) then
                  if (.not.rtdb_put(rtdb,'vib:animate:overwrite_files',
     &                              mt_log,1,.true.))
     &              call errquit('freq_vib_input: rtdb_put failed',911,
     &                           RTDB_ERR)
                else
                    write(luout,*)'freq_vib_input: could not parse line'
                    call inp_outrec()
                    goto 01001
                endif
              endif
            endif
          else
            write(luout,*)'freq_vib_input: could not parse line'
            call inp_outrec()
            goto 01001
          endif
          goto 00010
        else
          write(luout,10000)
          write(luout,*)'unrecognized keyword',test
          call errquit('freq_vib_input: fatal parsing error',911,
     &       INPUT_ERR)
        endif
      else
        write(luout,10000)
        call errquit('freq_vib_input: malformed input block ',911,
     &       INPUT_ERR)
      endif
01000 continue
      if (o_masses) then
        if (.not.rtdb_put
     &      (rtdb,'vib:remass:count',mt_int,1,num_mod_mass))
     &      call errquit ('freq_vib_input: could not put to rtdb ',911,
     &       RTDB_ERR)
      endif
      if (o_reuse) then
        if (.not.(rtdb_put(rtdb,'vib:reuse',mt_log,1,.true.)))
     &      call errquit ('freq_vib_input: could not put to rtdb ',911,
     &       RTDB_ERR)
      endif
      if (o_hessfile) then
        if (.not.rtdb_cput(rtdb,'vib:reuse_hessian_file',1,filename))
     &      call errquit ('freq_vib_input: could not cput to rtdb ',911,
     &       RTDB_ERR)
      endif
*debug-s
*     if (.not.rtdb_print(rtdb,.true.)) call errquit('911',911)
*debug-e
      return
01001 continue
      write(luout,10000)
      call errquit ('freq_vib_input: fatal input parsing error',911,
     &       INPUT_ERR)
10000 format(' frequency or vib input block format is: ',/,
     &       ' {freq || vib} ',/,
     &       '   reuse [<string> hessian_filename] ',/,
     &       '   mass <integer> lexical_index  <real> new_mass',/,
     &       '   mass <string>  tag_identifier <real> new_mass',/,
     &       '   animate [<real> step_size_for_animation ',
     &       '<default 0.15au>]',/,
     &       ' end ',/,/,
     &       ' NOTE: you must use a separate line for EACH', 
     &       ' redefined mass')
*             123456789012345  678901
10001 format('vib:remass:',i4,':lexi:')
10002 format('vib:remass:',i4,':tags:')
10003 format('vib:remass:',i4,':mass:')
      end



