C> \ingroup nwxc
C> @{
C>
C> \brief A utility program to print Maxima expressions of functionals
C>
C> In order to use Maxima to implement derivatives of the density
C> functionals we need the corresponding equations as Maxima
C> expressions. Given that the equations are already implemented
C> using an automatic differentiation library we can modify the
C> overloaded operators and functions to print what they do in the
C> form of Maxima expressions. After that we just need to run the code
C> to get the Maxima formulation. Doing this as part of the mainline
C> NWChem package is clunky as it brings in a tonne of code that we
C> do not need, and in normal NWChem runs you definitely do not want
C> to generate Maxima stuff. 
C>
C> Hence this separate little program is provided specifically to
C> generate the Maxima code in a form that is as good as possible for
C> having Maxima differentiate it and produce Fortran. The code
C> simply initializes the NWXC library and specifies the desired 
C> functional. It opens an output file for the Maxima expression and
C> hooks the file unit up to print instrumented automatic
C> differentiation library. Then the functional is evaluated for a few
C> points: one regular point, and one point with no beta electron
C> density (for G). After that the output file is closed and the
C> program terminates.
C>
C> The input for this program is almost trivial and tries to mimick
C> the regular NWChem input as much as possible. So the structures is
C>
C>   print_df
C>     file <filename>
C>     xc <functional definition>
C>   end
C>   task print_df
C>
C> The task "print_df" is the only task the program recognizes 
C> (everything else will cause the program to fail). The "print_df"
C> block specifies what to print and where to put the result. 
C> 
      program nwxc_fnc_print
C     use nwadp1
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
      logical  nwxc_fnc_print_input_parse
      external nwxc_fnc_print_input_parse
      character*(nw_max_path_len) input_filename
      integer stack, heap
      stack = 1000000
      heap  = 1000000
c
      call get_input_filename(input_filename)
      open(LuIn, file=input_filename, status='old', form='formatted',
     $     err=10000)
      goto 10001
10000 call errquit('nwxc_fnc_print: failed to open the input file', 0,
     &             INPUT_ERR)
10001 continue
      if (.not.ma_init(MT_DBL,stack,heap))
     &  call errquit("nwxc_fnc_print: ma_init failed",911,MA_ERR)
c
 10   if (nwxc_fnc_print_input_parse()) then
        call nwxc_fnc_print_task()
        goto 10
      endif
c
      end
C>
C> \brief Write the functional expression
C>
      subroutine nwxc_fnc_print_task
      use nwadp1
      implicit none
#include "inp.fh"
      double precision rho(2),    tr(2), dr(2)
      double precision rgamma(3), tg(3), dg(3)
      double precision tau(2),    tt(2), dt(2)
      double precision f, weight
c
      integer ii, func_id
c
      integer irefbegin, irefend
      integer idoibegin, idoiend
c
      character*(1024) description
      character*(1)    quote
c
      logical  nwxc_get_df, nwxc_get_info
      external nwxc_get_df, nwxc_get_info
c
      data rho    / 8.2309d-3, 1.1904d-3 /
      data rgamma / 2.6544d-4, 5.1828d-5, 1.0119d-5 /
      data tau    / 5.2247d-3, 1.0629d-3 /
      quote = '"'
      ii = 1
      do while (nwxc_get_df(ii,func_id,weight))
        if (weight.gt.0.0d0) then
          if (nwxc_get_info(func_id,description)) then
c
c           description now contains a string of the structure:
c           "<name> (<reference> doi:<doi>)"
c           now we want to find and print the reference and the doi.
c
            if (inp_contains(.false.,"(",description,irefbegin)) then
              irefbegin = irefbegin+1
              if (inp_contains(.false.,"doi:",description(irefbegin:),
     +                         idoibegin)) then
c
c               adjust idoibegin for the bit we cut of "description" in
c               inp_contains
c
                idoibegin = idoibegin + irefbegin - 1
c
c               irefend is of course the character before "doi:" starts
c
                irefend = idoibegin-1
c
c               we want to skip the "doi:" bit for the real DOI
c
                idoibegin = idoibegin + 4
                do while (description(idoibegin:idoibegin).eq." ")
                  idoibegin = idoibegin + 1
                enddo
c
                idoiend = inp_strlen(description)
                if (description(idoiend:idoiend).eq.")")
     +              idoiend=idoiend-1
                if (idoibegin.gt.0.and.idoiend.gt.idoibegin) then
                  write(nwad_iunit,'("doi: ",a,a,a,";")')
     +                  quote,description(idoibegin:idoiend),quote
                endif
              else
                irefend = inp_strlen(description)
                if (description(irefend:irefend).eq.")")
     +              irefend=irefend-1
              endif
              if (irefbegin.gt.0.and.irefend.gt.irefbegin) then
                write(nwad_iunit,'("reference: ",a,a,a,";")')
     +                quote,description(irefbegin:irefend),quote
              endif
            endif
          endif
        endif
      enddo
      call set_counter(-1)
      do ii = 1, 2
        tr(ii) = rho(ii)
        tt(ii) = tau(ii)
        dr(ii) = 0.0d0
        dt(ii) = 0.0d0
      enddo
      do ii = 1, 3
        tg(ii) = rgamma(ii)
        dg(ii) = 0.0d0
      enddo
      call nwxc_eval_df_print(2,1,tr,tg,tt,f,dr,dg,dt)
      tr(2) = 0.0d0
      tg(2) = 0.0d0
      tg(3) = 0.0d0
      call nwxc_eval_df_print(2,1,tr,tg,tt,f,dr,dg,dt)
      end
C>
C> \brief Parse the local mini input
C>
      logical function nwxc_fnc_print_input_parse()
      use nwadp1
      implicit none
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "nwxc.fh"
      integer nline
      integer ntask
      character*(nw_max_path_len) test, test1, func, filename
      double precision wght
      integer ilen
      nline = 0
      ntask = 0
      close(60)
 10   if (inp_read()) then
        nline = nline + 1
        if (.not. inp_a(test)) then
          call errquit("nwxc_fnc_print_input_parse: failed to read "
     &                //"keyword",911,INPUT_ERR)
        endif
        if (inp_compare(.false.,test,'task')) then
c
c         Read the "TASK" directive
c
          ntask = ntask + 1
          if (.not. inp_a(test1)) then
            call errquit("nwxc_fnc_print_input_parse: failed to read "
     &                  //"keyword",911,INPUT_ERR)
          endif
          if (inp_compare(.false.,test1,'print_df')) then
            go to 9999
          else
            call errquit("nwxc_fnc_print_input_parse: unknown task:"
     &                  //test1,911,INPUT_ERR)
          endif
c
        else if (inp_compare(.false.,test,'print_df')) then
c
c         Deal with the "PRINT_DF" input block which may contain
c         - file <filename> (to specify the output file)
c         - xc <functional> 
c         - end (to terminate the input block)
c
          call nwxc_init()
 20       if (inp_read()) then
            nline = nline + 1
            if (.not.inp_a(test)) then
              call errquit("nwxc_fnc_print_input_parse: failed to read "
     &                    //"keyword",911,INPUT_ERR)
            endif
            if (inp_compare(.false.,test,'end')) then
              go to 10
            else if (inp_compare(.false.,test,'file')) then
              if (.not.inp_a(filename)) then
                call errquit("nwxc_fnc_print_input_parse: failed to "
     &                     //"read filename",911,INPUT_ERR)
              endif
              open(UNIT=60,FILE=filename,ERR=30)
              call set_unitno(60)
            else if (inp_compare(.false.,test,'xc')) then
 40           if (.not.inp_a(func)) then
                go to 20
              endif
              ilen = inp_strlen(func)
              if (.not.inp_f(wght)) then
                wght = 1.0d0
              endif
              if (.not.nwxc_add_df(func(1:ilen),wght)) then
                call errquit("nwxc_fnc_print_input_parse: unknown "
     &                     //"functional: "//func(1:ilen),911,INPUT_ERR)
              endif
              go to 40
            else
              call errquit("nwxc_fnc_print_input_parse: unknown "
     &                   //"directive:"//test,911,INPUT_ERR)
            endif
            go to 20
          endif
c
        else 
          call errquit("nwxc_fnc_print_input_parse: unknown directive:"
     &                //test,911,INPUT_ERR)
        endif
        go to 10
      endif
 9999 continue
      call nwxc_printP()
      call nwxc_detect_dependencies()
      nwxc_fnc_print_input_parse = ntask.gt.0
      return
 30   call errquit("nwxc_fnc_print_input_parse: could not open file: "
     &            //filename,911,INPUT_ERR)
      end
C> @}
