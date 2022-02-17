#!/usr/bin/env bash
rm -f ptstride.patch
cat > ptstride.patch <<EOF
--- $1/comex/src-mpi-pt/comex.c
+++ $1/comex/src-mpi-pt/comex.c
@@ -2432,7 +2432,7 @@ STATIC void _put_packed_handler(header_t *header, int proc)
 
     unpack(packed_buffer, mapped_offset,
             stride->stride, stride->count, stride->stride_levels);
-
+    free(stride);
     if ((unsigned)header->length > COMEX_STATIC_BUFFER_SIZE) {
         free(packed_buffer);
     }

EOF
patch -p0 -s -N < ptstride.patch
echo ptstride.patch applied
