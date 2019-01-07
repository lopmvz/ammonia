##!/bin/bash
#
for i in *out
do
python2 ~/Workdir/ammonia/read.py  ${i} > ${i}.txt
if [ -s "${i}.txt" ]
then
python2 ~/Workdir/ammonia/conv.py ${i}.txt > ${i}.dat
fi
done
