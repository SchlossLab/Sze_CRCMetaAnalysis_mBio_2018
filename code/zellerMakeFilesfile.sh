#!/bin/bash 


# Make a files file for mothur using bash commands
#Get read 1 fastq
ls *_1.fastq > read1.txt
# Get read 2 fastq
ls *_2.fastq > read2.txt
# Get Names to be used
ls *_1.fastq | cut -c1-9 > names.txt
#Combine everything together
paste -d '\t' names.txt read1.txt read2.txt > zeller.files