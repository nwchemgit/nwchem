#!/bin/bash
rm -f make.log
echo "starting sleep_loop.sh for command: " ${@}
"${@}" >& make.log &
pid=$!
echo "sleep_loop got pid" $pid
ps -p "$pid" > /dev/null

if  [[ "${?}" -ne 0 ]];
then
    echo ' job stopped ' "${?}"
    exit 1
else
    while :
    do
        if kill -0 "$pid" 2>/dev/null; then
            echo ' ==== ' `date` ' ==== '
            tail -1 make.log
        elif wait "$pid"; then
            break          # exit loop.
        fi
        sleep 60s
    done
    exit 0
fi
echo 'exited loop '
tail -1 make.log
