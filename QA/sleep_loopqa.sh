#!/bin/bash
#echo "starting sleep_loop.sh for command: " ${@}
#"${@}"  &
narg="${#}"
if [ "$narg" -eq 10 ]; then
    outfile=../testoutputs/"${10}"
 $1 $2 $3 $4 $5 $6 $7 $8 $9 >& $outfile &
elif [ "$narg" -eq 9 ]; then
 outfile=../testoutputs/$9
 $1 $2 $3 $4 $5 $6 $7 $8    >& $outfile &
elif [ "$narg" -eq 8 ]; then
 outfile=../testoutputs/$8
 $1 $2 $3 $4 $5 $6 $7       >& $outfile &
elif [ "$narg" -eq 7 ]; then
 outfile=../testoutputs/$7
 $1 $2 $3 $4 $5 $6          >& $outfile &
else
 outfile=../testoutputs/$6
 $1 $2 $3 $4 $5     >& $outfile &
fi
pid=$!
echo "sleep_loopqa got pid" $pid
ps -p "$pid" > /dev/null
if  [[ "${?}" -ne 0 ]];
then
    echo ' job stopped ' "${?}"
    exit 1
else
    while :
    do
	ps -p "$pid" > /dev/null
	ps_exit="${?}"
        if [[ "$ps_exit" -eq 0 ]]; then
            echo ' ==== ' `date` ' ==== '
            tail -1 $outfile
	else
            break          # exit loop.
	fi
        sleep 30s
    done
    exit 0
fi
#tail -1 make.log
