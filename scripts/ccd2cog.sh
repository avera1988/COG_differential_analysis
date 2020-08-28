#!/bin/bash
############################################################################
#	This script uses cdd2go.pl from https://github.com/aleimba/bac-genomics-scripts/tree/master/cdd2cog
#	to annotate COG categories after rps-blast
#	Author Arturo Vera
#	Apr 2020
#	avera@ccg.unam.mx
############################################################################

print_usage(){
	echo "Usage: $0 COG_database"
}

if [ #$ -le 0 ]
then
	print_usage
	exit 1
fi

database=$1
wd=$(pwd)
for i in *.out;
        do
        mkdir $i.dir
        cd $i.dir
        ln -s ../$i .
        perl ~/scripts/COG/cdd2go.pl -r $i -c $database/cddid.tbl -f $database/fun.txt -w $database/whog
        cd $wd
done
