
get_elpa(){
    if [[ $# -eq 0 ]] ; then
	elpa_shortv=2025.06.002
    else
	elpa_shortv=$1
    fi
    echo "Parameter #1 is $1"
#    VERSION=new_release_${elpa_shortv}
    VERSION=${elpa_shortv}
    echo ELPA VERSION is $VERSION
    if [ -f  elpa-${VERSION}.tar.gz ]; then
	echo "using existing"  elpa-${VERSION}.tar.gz
    else
	rm -rf elpa*
	ELPA_URL=("https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/${VERSION}/elpa-${VERSION}.tar.gz" \
		      "https://web.archive.org/web/20260214233634/https://gitlab.mpcdf.mpg.de/elpa/elpa/-/archive/new_release_2025.06.002/elpa-new_release_2025.06.002.tar.gz")
	for url in  "${ELPA_URL[@]}"
	do
	    echo using $url
	    tries=1
	    until [ "$tries" -ge 3 ]
	    do
		if [ "$tries" -gt 1 ]; then echo sleeping for 9s ;sleep 9; echo attempt no.  $tries ; fi
		curl -L --progress-bar  $url -o elpa-${VERSION}.tar.gz
		gzip -t elpa-${VERSION}.tar.gz >&  /dev/null
		if [ $? -eq 0 ]; then break ;  fi
		tries=$((tries+1))
	    done
	    gzip -t elpa-${VERSION}.tar.gz >&  /dev/null
	    if [ $? -eq 0 ]; then return 0 ;  fi
	done
    fi
    return 1
}
get_scalapack(){
    if [[ $# -eq 0 ]] ; then
	#	version=2.1.0
 	COMMIT=a23c2cdc6586c427686f6097ae66bb54ef693571
	#COMMIT=bd1768b91262b4cdc7dd5f87b373b9b18eda4636
	#COMMIT=b935167ca4d244735abc04a3cd4f6d56699702a0
    else
	COMMIT=$1
    fi
    echo "Parameter #1 is $1"
    echo Scalapack commit is $COMMIT
    rm -rf scalapack 
    if [[ -f "scalapack-$COMMIT.tar.gz" ]]; then
	echo "using existing"  "scalapack-$COMMIT.tar.gz"
    else
	echo "downloading"  "scalapack-$COMMIT.tar.gz"
	rm -f scalapack-$COMMIT.tar.gz
	tries=1
	until [ "$tries" -ge 6 ]
	do
	    if [ "$tries" -gt 1 ]; then sleep 9; echo attempt no.  $tries ; fi
	    curl -L https://github.com/Reference-ScaLAPACK/scalapack/archive/$COMMIT.tar.gz -o scalapack-$COMMIT.tar.gz
	    # check tar.gz integrity
	    gzip -t scalapack-$COMMIT.tar.gz >&  /dev/null
	    if [ $? -eq 0 ]; then return 0 ;  fi
	    tries=$((tries+1)) ;  done
    fi
    return 1
}

get_openblas()
{
    if [[ $# -eq 0 ]] ; then
	VERSION=0.3.29
    else
	VERSION=$1
    fi
    echo "Parameter #1 is $1"
    echo OpenBLAS VERSION is $VERSION
    if [ -f  OpenBLAS-${VERSION}.tar.gz ]; then
	echo "using existing"  OpenBLAS-${VERSION}.tar.gz
    else
	rm -rf OpenBLAS* openblas*
	tries=1
	until [ "$tries" -ge 6 ]
	do
	    if [ "$tries" -gt 1 ]; then sleep 9; echo attempt no.  $tries ; fi
	    curl -L https://github.com/OpenMathLib/OpenBLAS/archive/v${VERSION}.tar.gz -o OpenBLAS-${VERSION}.tar.gz ;
	    # check tar.gz integrity
	    gzip -t OpenBLAS-${VERSION}.tar.gz >&  /dev/null
	    if [ $? -eq 0 ]; then return 0 ;  fi
	    tries=$((tries+1)) ;  done
    fi
    if [ $? -ne 0 ]; then echo  "openBLAS tarball not ready"; rm -f OpenBLAS-${VERSION}.tar.gz; return 1 ; fi
    return 0
}
