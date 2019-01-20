##!/bin/bash
#
#for i in {1..22}
#do
#python ~/Workdir/NH3NH4/read.py  PE_NH4_4H2O_step${i}00_631G_crit2.out > PE_NH4_4H2O_step${i}00_631G_crit2.txt
#if [ -s "PE_NH4_4H2O_step${i}00_631G_crit2.txt" ]
#then
#python ~/Workdir/NH3NH4/conv.py PE_NH4_4H2O_step${i}00_631G_crit2.txt > PE_NH4_4H2O_step${i}00_631G_crit2.dat
#fi
#done
#
for i in {222..608..2}
do
python2 ~/Workdir/ammonia/read.py  PE_NH4_4H2O_step${i}0_nd.out > PE_NH4_4H2O_step${i}0.txt
if [ -s "PE_NH4_4H2O_step${i}0.txt" ]
then
python2 ~/Workdir/ammonia/conv.py PE_NH4_4H2O_step${i}0.txt > PE_NH4_4H2O_step${i}0.dat
fi
if [ -s "PE_IP_NH4_4H2O_step${i}0_IP.out" ]
then
python ~/Workdir/ammonia/get_ip.py PE_IP_NH4_4H2O_step${i}0_IP.out > PE_IP_NH4_4H2O_step${i}0_IP.dat
fi
done
#
#for i in {62..65}
#do
#python ~/Workdir/NH3NH4/read.py  PE_NH4_4H2O_step${i}00_nd.out > PE_NH4_4H2O_step${i}00_631G_crit2.txt
#if [ -s "PE_NH4_4H2O_step${i}00_631G_crit2.txt" ]
#then
#python ~/Workdir/NH3NH4/conv.py PE_NH4_4H2O_step${i}00_631G_crit2.txt > PE_NH4_4H2O_step${i}00_631G_crit2.dat
#fi
#done
#
