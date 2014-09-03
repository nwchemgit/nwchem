#if (defined(CRAY) || defined(WIN32) || defined(HITACHI))&&!defined(__crayx1)&&!defined(__MINGW32__)

#define print_range_  PRINT_RANGE
#define copy_range_   COPY_RANGE
#define add_range_    ADD_RANGE
#define dot_range_    DOT_RANGE
#define init_array_   INIT_ARRAY
#define scale_patch_  SCALE_PATCH
#define compare_patches_  COMPARE_PATCHES
#define util_mitob_   UTIL_MITOB
#define util_mdtob_   UTIL_MDTOB
#define util_drand_   UTIL_DRAND
#define util_timer_   UTIL_TIMER
#define register_ext_memory_ REGISTER_EXT_MEMORY

#elif defined(F2C2_)

#define print_range_  print_range__  
#define copy_range_   copy_range__   
#define add_range_    add_range__    
#define dot_range_    dot_range__    
#define init_array_   init_array__   
#define scale_patch_  scale_patch__  
#define compare_patches_ compare_patches__
#define util_mitob_   util_mitob__  
#define util_mdtob_   util_mdtob__   
#define util_drand_   util_drand__   
#define util_timer_   util_timer__   
#define register_ext_memory_ register_ext_memory__

#endif
extern void get_range( int ndim, int dims[], int lo[], int hi[]);
extern void new_range(int ndim, int dims[], int lo[], int hi[],
                             int new_lo[], int new_hi[]);
extern void print_range(char *pre,int ndim, int lo[], int hi[], char* post);
extern void print_subscript(char *pre,int ndim, int subscript[], char* post);
extern void print_distribution(int g_a);


/* $Id$ */
