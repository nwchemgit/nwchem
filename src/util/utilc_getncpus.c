#ifdef USE_HWLOC
#include <hwloc.h>
#endif
#ifdef XLFLINUX
void utilc_getncpus_(
#else
void utilc_getncpus(
#endif
		    int *ncpus) {
#ifdef USE_HWLOC  
  hwloc_topology_t t;
    
  hwloc_topology_init(&t);
  
  hwloc_topology_load(t);
  
  unsigned nbcores_type = hwloc_get_nbobjs_by_type(t, HWLOC_OBJ_CORE);
  
  hwloc_topology_destroy(t);
  *ncpus =  (int) nbcores_type;
#else
  *ncpus = 0;
#endif
}
