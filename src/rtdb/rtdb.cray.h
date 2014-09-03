
/*$Id$*/
#if (defined(CRAY) || defined(WIN32)) &&!defined(__crayx1) &&!defined(__MINGW32__)
#define  context_pop_              CONTEXT_POP 
#define  context_prefix_            CONTEXT_PREFIX 
#define  context_push_              CONTEXT_PUSH 
#define  context_rtdb_load_         CONTEXT_RTDB_LOAD 
#define  context_rtdb_match_        CONTEXT_RTDB_MATCH 
#define  context_rtdb_store_        CONTEXT_RTDB_STORE 
#define  context_set_               CONTEXT_SET
#define  context_get_               CONTEXT_GET
#define  rtdb_cget_                 RTDB_CGET
#define  rtdb_close_                RTDB_CLOSE
#define  rtdb_cput_                 RTDB_CPUT
#define  rtdb_delete_               RTDB_DELETE
#define  rtdb_first_                RTDB_FIRST
#define  rtdb_get_                  RTDB_GET
#define  rtdb_get_info_             RTDB_GET_INFO
#define  rtdb_ma_get_               RTDB_MA_GET
#define  rtdb_next_                 RTDB_NEXT
#define  rtdb_open_                 RTDB_OPEN
#define  rtdb_parallel_             RTDB_PARALLEL
#define  rtdb_put_                  RTDB_PUT
#define  rtdb_print_                RTDB_PRINT
#define  rtdb_print_usage_          RTDB_PRINT_USAGE
#endif



