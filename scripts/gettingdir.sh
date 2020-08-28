#!/bin/bash
###############################################################################
#This script obtains all the func stats from the main directory 
#Author Arturo Vera
###############################################################################

wd=$(pwd)

for i in *.dir
        do
        a=$(basename $i .fasta.rps-blast.out.dir)
        cd $i/results
        cp func_stats.txt $wd/$a.func_stats.tab
        cd $wd
done
