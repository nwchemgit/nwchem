      function get_doit_time()
      implicit none
      double precision get_doit_time
#include "process_task.fh"

      get_doit_time = t_doit
      end !get_doit_time


      function get_gm_crt_time()
      implicit none
      double precision get_gm_crt_time
#include "process_task.fh"

      get_gm_crt_time = t_gm_crt
      end !get_gm_crt_time

      function get_rtdb_time()
      implicit none
      double precision get_rtdb_time
#include "process_task.fh"

      get_rtdb_time = t_rtdb
      end !get_rtdb_time


      function get_esp_time()
      implicit none
      double precision get_esp_time
#include "process_task.fh"

      get_esp_time = t_esp
      end !get_esp_time

      function get_gm_push_time()
      implicit none
      double precision get_gm_push_time
#include "process_task.fh"

      get_gm_push_time = t_gm_push
      end !get_push_time

      function get_chrg_set_time()
      implicit none
      double precision get_chrg_set_time
#include "process_task.fh"

      get_chrg_set_time = t_chrg_set
      end !get_chrg_set_time

      function get_im_end_time()
      implicit none
      double precision get_im_end_time
#include "process_task.fh"

      get_im_end_time = t_im_end
      end !get_im_end_time


c $Id$
