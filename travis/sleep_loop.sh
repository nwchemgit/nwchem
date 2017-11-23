#!/bin/bash -f
$* &
pid=$!
echo "TT got pid" $pid
while :
do
    sleep 60
    if kill -0 "$pid" 2>/dev/null; then
	echo ' ==== ' `date` ' ==== '
    elif wait "$pid"; then
	break          #Abandon the loop.
    fi
done
echo 'exited loop '
exit 0
