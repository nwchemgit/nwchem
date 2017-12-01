#!/bin/bash
echo "starting sleep_loop.sh for command: " ${@}
#"${@}"  &
outfile=../testoutputs/"$6"
echo "output files is "$outfile
$1 $2 $3 $4 $5 >& $outfile &
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
	elif wait "$pid"; then
	    echo "ps_exit code" $ps_exit
            break          # exit loop.
	fi
        sleep 30s
    done
    exit 0
fi
#tail -1 make.log
