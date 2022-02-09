#!/usr/bin/env bash
rm -f nodelist.patch
cat > nodelist.patch <<EOF
--- $1/comex/src-armci/armci.c
+++ $1/comex/src-armci/armci.c
@@ -91,6 +91,7 @@ void armci_init_domains(MPI_Comm comm)
       nodelist[i] = _my_node_id*_number_of_procs_per_node+i;
     comex_group_create(_number_of_procs_per_node, nodelist,
         COMEX_GROUP_WORLD, &ARMCI_Node_group);
+    free(nodelist);
   }
 }
 

EOF
patch -p0 -s -N < nodelist.patch
echo nodelist.patch applied
