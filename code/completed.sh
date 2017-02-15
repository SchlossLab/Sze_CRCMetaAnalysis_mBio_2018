#!/bin/sh

PROC=data/process
STUDIES="brim geng"


for s in $STUDIES
do
	echo $s >> $PROC/complete.txt
	ls $PROC/$s/*.seqs | wc -l >> $PROC/complete.txt
	echo "" >> $PROC/complete.txt
done


 
