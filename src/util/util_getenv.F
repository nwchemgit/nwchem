      subroutine util_getenv(env_name, env_value)
* $Id: util_getenv.F,v 1.7 2004-01-19 18:19:11 edo Exp $
      implicit none
#include "errquit.fh"
#include "stdio.fh"      
      character*(*) env_name  ! [input] the name of the environment variable to check for
      character*(*) env_value ! [output] the value of the named environement variable
#if defined(CRAY)
      integer*4 lname, lvalue, lerror
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
        call errquit('util_getenv:pxfgetenv fatal error',
     *        lerror, UNKNOWN_ERR)
      endif
#else
      call getenv(env_name, env_value)
#endif
*
      end
