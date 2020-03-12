#!/bin/bash
#echo "starting sleep_loop.sh for command: " ${@} >& /tmp/out
for last_arg in $@; do :; done
# echo "last arg of ${#} is $last_arg" >& /tmp/out
#"${@}"  &
set -- "${@:1:$(($#-1))}"
outfile=../testoutputs/$last_arg
"${@}" >& $outfile &
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
            tail -3 $outfile
	else
            break          # exit loop.
	fi
        sleep 30s
    done
    exit 0
fi
#tail -1 make.log
