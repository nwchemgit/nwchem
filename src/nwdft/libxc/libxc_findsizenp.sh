#!/usr/bin/env bash
npstring=$(egrep np "$1"/xc.h|head -1)
np_array=($npstring)
count=-1
sizenp_string=' '
for substring in "${np_array[@]}"; do
    [[ $substring == "np," ]] && sizenp_string=${np_array[$count]} && break
    ((++count))
done
if [[ "$sizenp_string" == 'int' ]]; then
    size_np=4	       	       
elif [[ "$sizenp_string" == 'size_t' ]]; then
    size_np=8
else
    echo "unexpected sizenp_string $sizenp_string"
    exit 1
fi
echo $size_np
