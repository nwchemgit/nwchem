
*     *********************************************
*     *                                           *
*     *         nwpw_parse_boolean                *
*     *                                           *
*     *********************************************
      logical function nwpw_parse_boolean(word,default)
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
      if (inp_compare(.false.,word,'True'))  value = .true.
      if (inp_compare(.false.,word,'T'))       value = .true.
      if (inp_compare(.false.,word,'off'))     value = .false.
      if (inp_compare(.false.,word,'no'))      value = .false.
      if (inp_compare(.false.,word,'.false.')) value = .false.
      if (inp_compare(.false.,word,'False')) value = .false.
      if (inp_compare(.false.,word,'F'))       value = .false.
      if (inp_compare(.false.,word,'fix'))     value = .false.

      nwpw_parse_boolean = value
      return
      end
c $Id$
