      subroutine util_getenv(env_name, env_value)
      implicit none
#include "stdio.fh"      
      character*(*) env_name  ! [input] the name of the environment variable to check for
      character*(*) env_value ! [output] the value of the named environement variable
*
      call getenv(env_name, env_value)
*
      end
