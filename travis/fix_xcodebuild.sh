#!/bin/bash
rm -f xcodebuild
INTEL_OSXSDK_VER=`xcodebuild -sdk macosx -version | grep SDKVersion`
INTEL_OSXSDK_PATH=`xcodebuild -sdk macosx -version Path`
echo $INTEL_OSXSDK_PATH
cat > xcodebuild <<EOF
#!/bin/bash
#https://community.intel.com/t5/Intel-oneAPI-HPC-Toolkit/slow-execution-of-ifort-icpc-on-MacOSX-catalina/m-p/1203190
INTEL_OSXSDK_VER="$INTEL_OSXSDK_VER"
INTEL_OSXSDK_PATH="$INTEL_OSXSDK_PATH"
case "\$4" in
    "")
      echo \$INTEL_OSXSDK_VER;;
     *)
      echo \$INTEL_OSXSDK_PATH;;
esac
EOF
chmod +x xcodebuild
