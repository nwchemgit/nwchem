      subroutine util_getenv(env_name, env_value)
* $Id: util_getenv.F,v 1.4 1998-05-18 18:15:55 d3e129 Exp $
      implicit none
#include "stdio.fh"      
      character*(*) env_name  ! [input] the name of the environment variable to check for
      character*(*) env_value ! [output] the value of the named environement variable
#if defined(CRAY)
      integer lname, lvalue, lerror
#endif
*
#if defined(CRAY)
      lname  = 0
      lvalue = 0
      lerror = 0
      call pxfgetenv(env_name,lname,env_value,lvalue,lerror)
      if (lerror.eq.0) then  ! 0 ok found env_name
	return
      else if (lerror.eq.22) then  ! 22 env_name is not in environment
        env_value = ' '
        return
      else
        call errquit('util_getenv:pxfgetenv fatal error',lerror)
      endif
#else
      call getenv(env_name, env_value)
#endif
*
      end
