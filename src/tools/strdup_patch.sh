#!/usr/bin/env bash
rm -f strdup.patch
cat > strdup.patch <<EOF
--- $1/tcgmsg/fapi.c
+++ $1/tcgmsg/fapi.c
@@ -200,6 +200,7 @@ void FATR _PBEGINF_()
 
     argv[argc] = 0;
     tcgi_pbegin(argc, argv);
+    for (i=0; i<argc; i++) free(argv[i]);
     free(argv);
 }
 

EOF
patch -p0 -s -N < strdup.patch
echo strdup.patch applied

