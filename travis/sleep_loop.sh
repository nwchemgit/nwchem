#!/bin/bash
rm -f make.log
echo "starting sleep_loop.sh for command: " ${@}
"${@}" >& make.log &
pid=$!
echo "sleep_loop got pid" $pid
sleep 1
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
            tail -1 make.log
        elif [[ "$ps_exit" -ne 0 ]]; then
            break          # exit loop.
	fi
        sleep 10s
    done
    exit 0
fi
echo 'exited loop '
tail -1 make.log
