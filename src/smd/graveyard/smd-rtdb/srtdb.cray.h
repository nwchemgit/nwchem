
/*$Id$*/
#if (defined(CRAY) || defined(WIN32)) &&!defined(__crayx1)
#define  context_pop_              CONTEXT_POP 
#define  context_prefix_            CONTEXT_PREFIX 
#define  context_push_              CONTEXT_PUSH 
#define  context_srtdb_load_         CONTEXT_SRTDB_LOAD 
#define  context_srtdb_match_        CONTEXT_SRTDB_MATCH 
#define  context_srtdb_store_        CONTEXT_SRTDB_STORE 
#define  context_set_               CONTEXT_SET
#define  context_get_               CONTEXT_GET
#define  srtdb_cget_                 SRTDB_CGET
#define  srtdb_close_                SRTDB_CLOSE
#define  srtdb_cput_                 SRTDB_CPUT
#define  srtdb_delete_               SRTDB_DELETE
#define  srtdb_first_                SRTDB_FIRST
#define  srtdb_get_                  SRTDB_GET
#define  srtdb_get_info_             SRTDB_GET_INFO
#define  srtdb_ma_get_               SRTDB_MA_GET
#define  srtdb_next_                 SRTDB_NEXT
#define  srtdb_open_                 SRTDB_OPEN
#define  srtdb_parallel_             SRTDB_PARALLEL
#define  srtdb_put_                  SRTDB_PUT
#define  srtdb_print_                SRTDB_PRINT
#define  srtdb_print_usage_          SRTDB_PRINT_USAGE
#endif



