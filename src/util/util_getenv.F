      subroutine util_getenv(env_name, env_value)
* $Id: util_getenv.F,v 1.2 1998-02-12 20:56:25 d3e129 Exp $
      implicit none
#include "stdio.fh"      
      character*(*) env_name  ! [input] the name of the environment variable to check for
      character*(*) env_value ! [output] the value of the named environement variable
*
      call getenv(env_name, env_value)
*
      end
