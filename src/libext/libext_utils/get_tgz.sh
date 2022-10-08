check_tgz() {
    gunzip -t $1 2> /dev/null
    myexit=$?
    echo $myexit
}


get_tgz(){
    echo "Parameter #1 is $1"
    echo "Parameter #2 is $2"
    URL=$1
    TGZ=$2
    if [[ `check_tgz $TGZ` == 0 ]]; then
	echo "using existing $TGZ"
    else
        tries=0 ; until [ "$tries" -ge 5 ] ; do \
		      echo "downloading $TGZ " \
			  && rm -f  $TGZ \
			  && curl -L $URL -o $TGZ  \
			  && if [[ `check_tgz $TGZ` == 0 ]]; then break; fi ;\
		      tries=$((tries+1)) ; echo attempt no. $tries failed  ; sleep 10 ;  done
	if [[ `check_tgz $TGZ` != 0 ]]; then
	    echo
	    echo  $TGZ download failed
	    echo
	fi
    fi
}
