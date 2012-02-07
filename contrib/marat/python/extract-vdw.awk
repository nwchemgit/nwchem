#!/usr/bin/awk -f 
BEGIN{
i=0
j=0
}
NF==1 {
j++
fn[j]=$1
}
NF==8 {
i++
a[i]=$1
an[i]=$5
w[i]=$6
r[i]=$7
vdw[i]=$NF
}
END{
n=i
for (i=1; i<=n; i++)
{
 printf "fullname['%s']='%s'\n",a[i],fn[i] > "fullname.txt"
 printf "vdw['%s']=%f\n",a[i],vdw[i] > "vdw.txt"
 printf "rc['%s']=%f\n",a[i],r[i] > "r.txt"
 printf "weight['%s']=%f\n",a[i],w[i] > "w.txt"
 printf "number['%s']=%d\n",a[i],an[i] > "an.txt"
}
}
