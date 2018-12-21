##!/bin/bash
#
#for i in {1..22}
#do
#python ~/Workdir/NH3NH4/read.py  PE_NH3_4H2O_step${i}00_631G_crit2.out > PE_NH3_4H2O_step${i}00_631G_crit2.txt
#if [ -s "PE_NH3_4H2O_step${i}00_631G_crit2.txt" ]
#then
#python ~/Workdir/NH3NH4/conv.py PE_NH3_4H2O_step${i}00_631G_crit2.txt > PE_NH3_4H2O_step${i}00_631G_crit2.dat
#fi
#done
#
for i in {222..410..10}
do
python2 ~/Workdir/NH3NH4/read.py  PE_NH3_4H2O_step${i}0_nd.out > PE_NH3_4H2O_step${i}0.txt
if [ -s "PE_NH3_4H2O_step${i}0.txt" ]
then
python2 ~/Workdir/NH3NH4/conv.py PE_NH3_4H2O_step${i}0.txt > NH3_4H2O_step${i}0.dat
fi
if [ -s "PE_IP_NH3_4H2O_step${i}0_IP.out" ]
then
python ~/Workdir/ammonia/get_ip.py PE_IP_NH3_4H2O_step${i}0_IP.out > PE_IP_NH3_4H2O_step${i}0_IP.dat
fi
done
#
#for i in {62..65}
#do
#python ~/Workdir/NH3NH4/read.py  PE_NH3_4H2O_step${i}00_nd.out > PE_NH3_4H2O_step${i}00_631G_crit2.txt
#if [ -s "PE_NH3_4H2O_step${i}00_631G_crit2.txt" ]
#then
#python ~/Workdir/NH3NH4/conv.py PE_NH3_4H2O_step${i}00_631G_crit2.txt > PE_NH3_4H2O_step${i}00_631G_crit2.dat
#fi
#done
#
